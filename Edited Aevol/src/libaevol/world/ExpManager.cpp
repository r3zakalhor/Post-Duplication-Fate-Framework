// ****************************************************************************
//
//          Aevol - An in silico experimental evolution platform
//
// ****************************************************************************
//
// Copyright: See the AUTHORS file provided with the package or <www.aevol.fr>
// Web: http://www.aevol.fr/
// E-mail: See <http://www.aevol.fr/contact/>
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ****************************************************************************




// =================================================================
//                              Includes
// =================================================================
#include "ExpManager.h"

#include "7/Dna_7.h"
#include "7/Individual_7.h"
#include "7/ExpManager_7.h"
#include "7/Stats_7.h"
#include "Individual.h"

#include <cerrno>
#include <err.h>
#include <iostream>
#include <list>
#include <sys/stat.h>
#include <zlib.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef __CUDACC__
#include "CUDA_Individual.h"
#include<cuda_profiler_api.h>
#endif

#ifdef __REGUL
#include "raevol/cuda_struct.h"
#include "raevol/Individual_R.h"
#endif
#include<chrono>

#include <iostream>
#include <unordered_map>
using namespace std;
using namespace std::chrono;

#ifdef __TRACING__
#include "ae_logger.h"


unordered_map<int,unordered_multiset<string>> ae_logger::logMap;
string ae_logger::logFile = "logger_csv.log";
mutex ae_logger::loggerMtx;
#endif

using std::cout;
using std::endl;


namespace aevol {
//##############################################################################
//                                                                             #
//                             Class ExpManager                            #
//                                                                             #
//##############################################################################

// ===========================================================================
//                        Definition of static attributes
// ===========================================================================

// ===========================================================================
//                                 Constructors
// ===========================================================================
ExpManager::ExpManager()
{
  // ------------------------------------------------------ Experimental setup
  exp_s_ = new ExpSetup(this);

  // ------------------------------------------------------------------- World
  world_ = nullptr;

  // ---------------------------------------------------------- Output manager
  output_m_ = new OutputManager(this);


  // -------------------------------- Timestep up to which we want to simulate
  t_end_ = 0;


  // ------------------------------------------------------------- Quit signal
  quit_signal_received_ = false;


  // ------------------------ Which Individuals will reproduce in each GridCell


#ifdef __PROXY_POW_APPROX
#ifdef __REGUL
  Rna_R::load_lookup_table();
#endif
#endif
}

// ===========================================================================
//                                  Destructor
// ===========================================================================
ExpManager::~ExpManager() noexcept
{
  if (ExpManager_7::standalone_simd) {
    delete exp_m_7_;
  }
  delete output_m_;
  delete world_;


  delete exp_s_;

  delete [] dna_mutator_array_;

  delete [] next_generation_reproducer_;

  // Several ExpManager instances can delete the fuzzy factory
  if (FuzzyFactory::fuzzyFactory != nullptr) {
    delete FuzzyFactory::fuzzyFactory;
    FuzzyFactory::fuzzyFactory = nullptr;
  }
}

// ===========================================================================
//                                   Algorithms
// ===========================================================================

// ===========================================================================
//                                 Public Methods
// ===========================================================================
void ExpManager::InitializeWorld(int16_t grid_width,
                                     int16_t grid_height,
                                     std::shared_ptr<JumpingMT> prng,
                                     std::shared_ptr<JumpingMT> mut_prng,
                                     std::shared_ptr<JumpingMT> stoch_prng,
                                     Habitat& habitat,
                                     bool share_phenotypic_target)
{
  world_ = new World();
  world_->set_prng(prng);
  world_->set_mut_prng(mut_prng);
  world_->set_stoch_prng(stoch_prng);
  world_->InitGrid(grid_width, grid_height, habitat, share_phenotypic_target);
  next_generation_reproducer_ = new int32_t[nb_indivs()];
}

/*!
  \brief Save the experiment
*/
void ExpManager::Save(bool create) const
{
  WriteSetupFiles(create);
  output_m_->write_current_generation_outputs(create);
}

/*!
  \brief Save all the static data (the data that is constant throughout
  a simulation)

  Save
    * The experimental setup
    * The output profile

  These are written in the current directory in gzipped binary files
  as well as plain text files
  This data should be saved once and for all when the experiment is created.
  There is no need to save this data in the middle of a simulation since it
  is constant throughout the experiment.

  \see load(int64_t t0,
            char* exp_setup_file_name,
            char* out_prof_file_name,
            char* pop_file_name,
            char* sel_file_name,
            char* world_file_name,
            bool verbose)
  \see WriteDynamicFiles()
  \see save_copy(char* dir, int64_t t)
*/
void ExpManager::WriteSetupFiles(bool first) const
{
  // 1) Create missing directories
  create_missing_directories();

  // 2) Open setup files (experimental setup and output profile)
  gzFile exp_s_file, out_p_file;

  if (first) {
    #ifdef HAVE_MPI
    for (int rank = 0; rank < nb_rank(); rank++) {
      open_setup_files(exp_s_file, out_p_file, AeTime::time(), "w", ".", rank);
      // 4) Write setup data
      exp_s_->write_setup_file(exp_s_file);
      output_m_->WriteSetupFile(out_p_file);

      // 5) Close setup files
      close_setup_files(exp_s_file, out_p_file);
    }
    #else
    open_setup_files(exp_s_file, out_p_file, AeTime::time(), "w");
    exp_s_->write_setup_file(exp_s_file);
    output_m_->WriteSetupFile(out_p_file);
    close_setup_files(exp_s_file, out_p_file);
    #endif
  } else {
    #ifdef HAVE_MPI
    open_setup_files(exp_s_file, out_p_file, AeTime::time(), "w", ".", rank_);
    #else
    open_setup_files(exp_s_file, out_p_file, AeTime::time(), "w");
    #endif

    // 4) Write setup data
    exp_s_->write_setup_file(exp_s_file);
    output_m_->WriteSetupFile(out_p_file);

    // 5) Close setup files
    close_setup_files(exp_s_file, out_p_file);
  }
  


}

/*!
  \brief Save all the dynamic data (the data that changes over the simulation)

  Save the state of
    * The population
    * The world
    * The PRNG (random generator) used for the selection

  These are written in the current directory in gzipped binary files
  whose names are appended with the current generation number.

  WARNING: The experimental setup and output profile are *not* saved. These
  should be saved once and for all when the experiment is created.

  \see load(int64_t t0,
            char* exp_setup_file_name,
            char* out_prof_file_name,
            char* pop_file_name,
            char* sel_file_name,
            char* world_file_name,
            bool verbose)
  \see WriteSetupFiles()
  \see save_copy(char* dir, int64_t t)
*/
void ExpManager::WriteDynamicFiles(bool first) const
{
  // Create missing directories
  create_missing_directories();

  
  // Save experiment
  if (first) {
  #ifdef HAVE_MPI
  for (int rank = 0; rank < nb_rank(); rank++) {
    int32_t rank_x = exp_s()->rank_width();
    int32_t rank_y = exp_s()->rank_height();

    int32_t local_rank_x = rank / rank_y;
    int32_t local_rank_y = rank % rank_y;

    int32_t x_offset = (local_rank_x) * grid_width(); 
    int32_t y_offset = (local_rank_y) * grid_height();

    for (int16_t x = 0; x < grid_width(); x++) {
      for (int16_t y = 0; y < grid_height(); y++) {
        int global_x = (x + x_offset) % exp_s()->global_grid_width();
        int global_y = (y + y_offset) % exp_s()->global_grid_height();
        // printf("Getting seed for %05d (%d %d -- %d %d): %d R %d (%d %d) Offset %d %d\n",global_x*exp_s()->global_grid_height()+global_y,
        //         x,y,global_x,global_y,
        //         prng_seed_[global_x*exp_s()->global_grid_height()+y],rank,
        //         local_rank_x,local_rank_y,x_offset,y_offset);
        int32_t seed = prng_seed_[global_x*exp_s()->global_grid_height()+global_y];

  #if __cplusplus == 201103L
        world()->grid(x,y)->set_reprod_prng(make_unique<JumpingMT>(seed));
        world()->grid(x,y)->set_reprod_prng_simd(make_unique<JumpingMT>(seed));
  #else
        world()->grid(x,y)->set_reprod_prng(std::make_unique<JumpingMT>(seed));
        world()->grid(x,y)->set_reprod_prng_simd(std::make_unique<JumpingMT>(seed));
  #endif
        int32_t mut_seed = mut_prng_seed_[global_x*exp_s()->global_grid_height()+global_y];
        int32_t stoch_seed = stoch_prng_seed_[global_x*exp_s()->global_grid_height()+global_y];

        world()->grid(x,y)->set_mut_prng(std::make_shared<JumpingMT>(mut_seed));
        world()->grid(x,y)->set_stoch_prng(std::make_shared<JumpingMT>(stoch_seed));

        // printf("Seed %d : %d\n",global_x*exp_s()->global_grid_height()+global_y,seed);
        // printf("Seed Grid %d : %d %d\n",global_x*exp_s()->global_grid_height()+global_y,mut_seed,stoch_seed);
      }
    }
#endif

    // Open backup files
    gzFile sel_file, world_file;
    #ifdef HAVE_MPI
    open_backup_files(sel_file, world_file, AeTime::time(), "w", ".", rank);
    #else
    open_backup_files(sel_file, world_file, AeTime::time(), "w");
    #endif


    sel()->save(sel_file);
    
    world_->save(world_file);
    
    // Close backup files
    close_backup_files(sel_file, world_file);
    #ifdef HAVE_MPI
  }
    #endif
  } else {
    // Open backup files
    gzFile sel_file, world_file;
    #ifdef HAVE_MPI
    open_backup_files(sel_file, world_file, AeTime::time(), "w", ".", rank_);
    #else
    open_backup_files(sel_file, world_file, AeTime::time(), "w");
    #endif


    sel()->save(sel_file);
    
    world_->save(world_file);
    
    // Close backup files
    close_backup_files(sel_file, world_file);
  }


}

    void ExpManager::WriteDynamicFiles(int64_t gen, SaveWorld* world, bool create) const
    {
      // Create missing directories
      create_missing_directories();

      // Open backup files
      gzFile sel_file, world_file;
      open_backup_files(sel_file, world_file, gen, "w");

      // Save experiment
      world->save(world_file, create);

      // Close backup files
      close_backup_files(sel_file, world_file);
    }

/*!
  \brief Saves a complete copy of the experiment at the provided location.

  Save both the setup (constant) and the backup (dynamic) files to the <dir>
  directory. Dynamic file names will be appended with the <num_gener>

  Static files (saved both as gzipped binary and plain text)
    * The experimental setup
    * The output profile

  Dynamic files (saved as gzipped binary)
    * The world
    * The PRNG (random generator) used for the selection

  \see load(int32_t first_gener,
             char* exp_setup_file_name,
             char* out_prof_file_name,
             char* pop_file_name,
             char* sel_file_name,
             char* world_file_name,
             bool verbose)
  \see WriteSetupFiles()
  \see WriteDynamicFiles()
*/
void ExpManager::save_copy(char* dir, int64_t time) const
{
  // Set time to time
  AeTime::set_time(time);

  // Create missing directories
  create_missing_directories(dir);

  // Open setup files and backup files
  gzFile exp_s_file, out_p_file, sel_file, world_file;
  open_setup_files(exp_s_file, out_p_file, time, "w", dir);
  open_backup_files(sel_file, world_file, time, "w", dir);

  // Write setup data
  exp_s_->write_setup_file(exp_s_file);
  output_m_->WriteSetupFile(out_p_file);

  // Write the state of selection and world into the backups
  sel()->save(sel_file);
  world_->save(world_file);

  // Close setup and backup files
  close_setup_files(exp_s_file, out_p_file);
  close_backup_files(sel_file, world_file);

  // Copy stats
  output_m_->CopyStats(dir, time);

  // Write last gener file
  output_m_->WriteLastGenerFile(dir);
}



void ExpManager::step_to_next_generation() {
  // TODO <david.parsons@inria.fr> Apply phenotypic target variation and noise
#ifdef __PERF_LOG__
    auto t1 = high_resolution_clock::now();
    auto t2 = high_resolution_clock::now();
    auto duration_2 = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    auto ta = high_resolution_clock::now();
#endif

#pragma omp single
    {
      #ifdef __REGUL
      if (!exp_m_7_->standalone() || (exp_m_7_->standalone()&&check_simd())) {
        world_->ApplyHabitatVariation();
      }

      if (exp_m_7_->standalone()) 
        exp_m_7_->phenotypic_target_handler_->ApplyVariation();

      #else
        world_->ApplyHabitatVariation();
      #endif

        // Take a step in time
        AeTime::plusplus();

        // Create the corresponding new generation
#ifdef __PERF_LOG__
        auto t1 = high_resolution_clock::now();

        //auto     s_t1 = high_resolution_clock::now();
#endif

#ifdef __CUDACC__
        if (first_gen) {
          cudaProfilerStop();
        }
        transfer_in(this, first_gen);


        auto t2 = high_resolution_clock::now();
        auto duration_A = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

        t1 = high_resolution_clock::now();

        run_a_step(nb_indivs(), (float) best_indiv()->w_max(),
                     selection_pressure(),first_gen);

        t2 = high_resolution_clock::now();
        auto duration_B = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

          /** Debug START **/
         // print_debug_promoters_start(this);
      //    print_debug_rna(this);
      //    print_debug_protein(this);
      //    print_debug_phenotype(this);
          //print_debug_fitness(this);
          /** Debug END **/

        t1 = high_resolution_clock::now();
          // Transfer data out of the GPU
        transfer_out(this,false);
        t2 = high_resolution_clock::now();
        auto duration_C = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - s_t1 ).count();
#endif

/*  bool simd_first = false;
  if (exp_m_7_ == nullptr) {
    simd_first = true;
  }*/

        //}
#ifdef __PERF_LOG__
        t1 = high_resolution_clock::now();
#endif
  }
        if (!ExpManager_7::standalone() || (ExpManager_7::standalone() && check_simd())) {
            exp_s_->step_to_next_generation();
        }

#pragma omp single
        {
#ifdef __CUDACC__
          t2 = high_resolution_clock::now();
#else
  #ifdef __PERF_LOG__
          auto t2 = high_resolution_clock::now();
  #endif
#endif

#ifdef __PERF_LOG__
          auto duration_2 =
              std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)
                  .count();

          auto ta = high_resolution_clock::now();
#endif
        }

    if (ExpManager_7::standalone())
      exp_m_7_->run_a_step(w_max_, selection_pressure());

    if (check_simd_) {
      exp_m_7_->check_result();
    }
  //}
#pragma omp single
    {
#ifdef __PERF_LOG__
        auto tb = high_resolution_clock::now();
        auto duration_simd = std::chrono::duration_cast<std::chrono::microseconds>(tb - ta).count();
#endif

        if (!ExpManager_7::standalone() || (ExpManager_7::standalone() && check_simd())) {
#ifdef __PERF_LOG__
            t1 = high_resolution_clock::now();
#endif
            // Write statistical data and store phylogenetic data (tree)
            {
                //printf("----------------> Must not be here\n");
                output_m_->write_current_generation_outputs();
            }

#ifdef __PERF_LOG__
            t2 = high_resolution_clock::now();
            duration_2 += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
#endif
        }

#ifdef __CUDACC__
        std::cout<<"PERFLOG,"<<AeTime::time()<<","<<duration_2<<","<<duration<<","<<duration_A<<","<<duration_B<<","<<duration_C<<std::endl;

        if (first_gen) {
          cudaProfilerStart();
          first_gen = false;
        }
#else

#ifdef __PERF_LOG__
        std::cout << "PERFLOG," << AeTime::time() << "," << duration_simd << ","<< exp_m_7_->cumulate_size << "," << exp_m_7_->cumulate_diff <<std::endl;
#endif

#endif
    }
}

/*!
  \brief Load an experiment with the provided files
 */
// TODO <david.parsons@inria.fr> check verbose (what does it do ?, is it consistent ?)
void ExpManager::load(gzFile& exp_s_file,
                          gzFile& exp_backup_file,
                          gzFile& world_file,
                          gzFile& out_p_file,
                          bool verbose,
                          bool to_be_run /*  = true */
                          #ifdef HAVE_MPI
                          , int32_t rank
                          #endif
                          )
{
  #ifdef HAVE_MPI
  rank_ = rank;
  #endif

  // ---------------------------------------- Retrieve experimental setup data
  printf("  Loading experimental setup...");
  fflush(stdout);
  exp_s_->load(exp_s_file, exp_backup_file, verbose);
  printf(" OK\n");

  if (FuzzyFactory::fuzzyFactory == NULL)
    FuzzyFactory::fuzzyFactory = new FuzzyFactory();


  // ---------------------------------------------------------- Retrieve world
  printf("  Loading world...");
  fflush(stdout);
  world_ = new World();
  world_->load(world_file, this);
  printf(" OK\n");

  //printf("-------------- Nb Gaussians %d\n",world()->phenotypic_target_handler()->gaussians().size());

  sel()->set_unique_id(  grid_height()*grid_width()+1 );

  dna_mutator_array_ = new DnaMutator*[grid_height()*grid_width()];
    for (int i = 0; i < grid_height()*grid_width(); i++) {
        dna_mutator_array_[i] = nullptr;
    }


  // --------------------------------------------- Retrieve output profile data
  printf("  Loading output profile...");
  fflush(stdout);
  output_m_->load(out_p_file, verbose, to_be_run);
  printf(" OK\n");

    if (ExpManager_7::standalone_simd) {
      // if (to_be_run) {
        exp_m_7_                     = new ExpManager_7(this);
        exp_m_7_->protein_grain_size = grain_size;
        exp_m_7_->rna_grain_size     = grain_size;
      // }
  }


  // -------------------------------------------- Link world and output profile
  if (record_tree()) {
      if (ExpManager_7::standalone_simd) {
          if (to_be_run) {
            exp_m_7_->addObserver(tree(), NEW_INDIV);
              for (int16_t x = 0; x < grid_width(); x++) {
                  for (int16_t y = 0; y < grid_height(); y++) {
                    exp_m_7_->current_individuals[x * grid_height() + y]->addObserver(
                              tree(),
                              END_REPLICATION);
                  }
              }
              exp_m_7_->addObserver(tree(), END_GENERATION);
          }
      } else {
          sel()->addObserver(tree(), NEW_INDIV);
          for (int16_t x = 0; x < grid_width(); x++) {
              for (int16_t y = 0; y < grid_height(); y++) {
                  world_->indiv_at(x, y)->addObserver(
                          tree(),
                          END_REPLICATION);
              }
          }
          sel()->addObserver(tree(), END_GENERATION);
      }
  }

  if (record_light_tree()){
    if (ExpManager_7::standalone_simd) {
        if (to_be_run) {
          exp_m_7_->addObserver(light_tree(), NEW_INDIV);
            for (int16_t x = 0; x < grid_width(); x++) {
                for (int16_t y = 0; y < grid_height(); y++) {
                  exp_m_7_->current_individuals[x * grid_height() + y]->addObserver(
                            light_tree(),
                            END_REPLICATION);
                }
            }
            exp_m_7_->addObserver(light_tree(), END_GENERATION);
        }
    } else {
      sel()->addObserver(light_tree(), NEW_INDIV);
      for (int16_t x = 0; x < grid_width(); x++) {
        for (int16_t y = 0; y < grid_height(); y++) {
          world_->indiv_at(x, y)->addObserver(
                  light_tree(),
                  END_REPLICATION);
        }
      }
      sel()->addObserver(light_tree(), END_GENERATION);
    }
  }

  // --------------------------------------------------- Recompute unsaved data
  world_->evaluate_individuals();
  //update_best();
}


/*!
  \brief Load an experiment with default files from a given directory
 */
// TODO <david.parsons@inria.fr> check verbose (what does it do ?, is it consistent ?)
void ExpManager::load(const char* dir,
    int64_t t0, bool verbose, bool to_be_run /*  = true */
    #ifdef HAVE_MPI
    , int32_t rank
    #endif
    )
{
  AeTime::set_time(t0);

  // -------------------------------------------------------------------------
  // Open setup files and backup files
  // -------------------------------------------------------------------------
  gzFile exp_s_file, out_p_file;
  gzFile exp_backup_file, world_file;
  open_setup_files(exp_s_file, out_p_file, t0, "r", dir
  #ifdef HAVE_MPI
  ,rank
  #endif
  );
  open_backup_files(exp_backup_file, world_file, t0, "r", dir
  #ifdef HAVE_MPI
  ,rank
  #endif
  );


  // -------------------------------------------------------------------------
  // Load data from backup and parameter files
  // -------------------------------------------------------------------------
  load(exp_s_file, exp_backup_file,
       world_file, out_p_file, verbose, to_be_run
       #ifdef HAVE_MPI
       , rank
       #endif
       );


  // -------------------------------------------------------------------------
  // Close setup and backup files
  // -------------------------------------------------------------------------
  close_setup_files(exp_s_file, out_p_file);
  close_backup_files(exp_backup_file, world_file);


  if (FuzzyFactory::fuzzyFactory == NULL)
    FuzzyFactory::fuzzyFactory = new FuzzyFactory();
}


/**
 * \brief Load an experiment with the provided constitutive files
 */
// TODO <david.parsons@inria.fr> check verbose (what does it do ?, is it consistent ?)
void ExpManager::load(int64_t t0,
                          char* exp_setup_file_name,
                          char* exp_backup_file_name,
                          char* world_file_name,
                          char* out_prof_file_name,
                          bool verbose /*= false*/,
                          bool to_be_run /*= true*/
                          #ifdef HAVE_MPI
                          , int32_t rank
                          #endif
                          )
{
  AeTime::set_time(t0);

  // ---------------------------------------------------------------------------
  // Open files and check them
  // ---------------------------------------------------------------------------
  gzFile exp_setup_file = gzopen(exp_setup_file_name, "r");
  gzFile out_prof_file = gzopen(out_prof_file_name, "r");
  gzFile exp_backup_file = gzopen(exp_backup_file_name, "r");
  gzFile world_file  = gzopen(world_file_name, "r");


  if (exp_setup_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, exp_setup_file_name);
    exit(EXIT_FAILURE);
  }
  if (out_prof_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, out_prof_file_name);
    exit(EXIT_FAILURE);
  }
  if (exp_backup_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, exp_backup_file_name);
    exit(EXIT_FAILURE);
  }
  if (world_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, world_file_name);
    exit(EXIT_FAILURE);
  }


  // ---------------------------------------------------------------------------
  // Actually load data
  // ---------------------------------------------------------------------------
  load(exp_setup_file, exp_backup_file, world_file, out_prof_file,
       verbose, to_be_run
       #ifdef HAVE_MPI
       , rank
       #endif
       );


  // ---------------------------------------------------------------------------
  // Close setup and backup files
  // ---------------------------------------------------------------------------
  gzclose(exp_setup_file);
  gzclose(exp_backup_file);
  gzclose(world_file);
  gzclose(out_prof_file);
}



//int mutator = 0;

/**
 * Run the simulation
 */
void ExpManager::run_evolution() {
  next_generation_reproducer_ = new int32_t[nb_indivs()];

  if (ExpManager_7::standalone()) {
    for (int indiv_id = 0; indiv_id < nb_indivs(); indiv_id++) {
      next_generation_reproducer_[indiv_id] = indiv_id;
    }

  }
  // We are running a simulation.
  // Save the setup files to keep track of the setup history
  WriteSetupFiles();
#ifdef __TRACING__
  ae_logger::init("logger_csv.log");

    printf("Launching TRACING...\n");
#else
  printf("Launching NOT TRACING...\n");
#endif

#ifdef __TRACING__
  high_resolution_clock::time_point t_t1 = high_resolution_clock::now();
    high_resolution_clock::time_point t_t2,t1,t2;
#endif

#ifdef __CUDACC__
  cuda_init();

  cudaProfilerStart();
#endif

  printf("Tree activated %d\n",record_tree());

        w_max_ = world()->indiv_at(0,0)->w_max();

  bool first_run = true;

      //"Post Treatment"
      if(anc_stat_ && record_light_tree()) {
        output_m_->light_tree()->setup_anc_stat();
      }

      //printf("Add indivs at gen %d\n",AeTime::time());

      //output_m_->stats()->add_indivs(AeTime::time(), indivs());

      if (ExpManager_7::standalone_simd) {
        exp_m_7_->setup_individuals(w_max_, selection_pressure());

        if (check_simd_)
          exp_m_7_->check_result();
      }

      bool finished=false;
//      int64_t time = AeTime::time();
        // For each generation
#pragma omp parallel default(shared)
                while (!finished) {
#pragma omp single
                        {
                            #ifdef HAVE_MPI
                            // MPI_Barrier(MPI_COMM_WORLD);

                          if (rank_==0) {
                             #endif

                            if (AeTime::time() % 1 == 0) {
                                printf(
                                        "============================== %" PRId64 " ==============================\n",
                                        AeTime::time());
                                if (!first_run) {
                                    if (ExpManager_7::standalone_simd) {
                                      exp_m_7_->dna_factory_->stats();
                                      exp_m_7_->fuzzy_factory_->stats();
                                      #ifndef HAVE_MPI
                                        printf(
                                                "  Best individual's distance to target (metabolic) : %f (clones %d %ld)\n",
                                                //exp_m_7_->best_indiv->indiv_id,
                                          exp_m_7_->best_indiv->metaerror,
                                          exp_m_7_->nb_clones_,exp_m_7_->mutant_list_.size());
                                      #endif
                                    } else {
                                        printf("  Aevol_6 : Best individual's distance to target (metabolic) : %f\n",
                                                //best_indiv()->id(),
                                               best_indiv()->dist_to_target_by_feature(METABOLISM));
                                    }
                                } else {
                                    printf("  Aevol_6 : Best individual's distance to target (metabolic) : %f\n",
                                            //best_indiv()->id(),
                                           best_indiv()->dist_to_target_by_feature(METABOLISM));
                                }
                            }

                            first_run = false;


#ifdef __X11
                            display();
#endif
#ifdef HAVE_MPI
                          }
#endif

#ifdef __TRACING__
                            t1 = high_resolution_clock::now();
#endif


                        }

#pragma omp single
                  {
                    if (with_mrca_ && record_light_tree()) {
                      /*if (AeTime::time() == t_end_) {
                          output_m_->light_tree()->keep_indivs(indivs());
                      }*/
                      if (output_m_->mrca_time() >= t_end_ or quit_signal_received())
                        finished = true;
                    } else if (AeTime::time() >= t_end_ or quit_signal_received())
                      finished = true;
                        }

        if (!finished)
                        // Take one step in the evolutionary loop
                        step_to_next_generation();
                        #ifdef HAVE_MPI
                        // TODO/ REMOVE §!!
                            // MPI_Barrier(MPI_COMM_WORLD);
                        #endif



                }
                                            #ifdef HAVE_MPI
                            MPI_Barrier(MPI_COMM_WORLD);
                            #endif
#ifdef __TRACING__
  t_t2 = high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t_t2 - t_t1 ).count();
        ae_logger::addLog(TOTAL,duration);
        ae_logger::flush(AeTime::time());
#endif

#ifdef __CUDACC__
  cudaProfilerStop();
#endif
   //
   // I do not think this is required, the world destructor must clean up this
   //
  // if (ExpManager_7::standalone() && check_simd_) {
  //   for (int indiv_id = 0; indiv_id < (int)nb_indivs(); indiv_id++) {
  //     int x = indiv_id / world_->height();
  //     int y = indiv_id % world_->height();

  //     world_->grid(x, y)
  //         ->individual()
  //         ->clear_everything_except_dna_and_promoters();
  //     world_->grid(x, y)->individual()->genetic_unit_list_nonconst().clear();
  //     delete world_->grid(x, y)->individual();
  //   }
  // }

  //  output_m_->flush();
//  if(with_mrca_ && record_light_tree())
//    output_m_->light_tree()->save_mrca_indiv();
//  if(anc_stat_ && record_light_tree())
//    output_m_->light_tree()->close_anc_stat();

  printf("================================================================\n");
  printf("  The run is finished. \n");
  printf("  Printing the final best individual into " BEST_LAST_ORG_FNAME "\n");

  if (!exp_m_7_->standalone()) {
    FILE* org_file = fopen(BEST_LAST_ORG_FNAME, "w");
    fputs(best_indiv()->genetic_unit_sequence(0), org_file);
    fclose(org_file);
  } else {
    // FILE* org_file = fopen(BEST_LAST_ORG_FNAME, "w");
    // printf("Last best indiv size %d\n",exp_m_7_->best_indiv->dna_->length_);
    // fputs(exp_m_7_->best_indiv->dna_->data_, org_file);
    // fclose(org_file);
  }
}

void ExpManager::update_best()
{
  world_->update_best();
}

// ===========================================================================
//                                Protected Methods
// ===========================================================================
void ExpManager::create_missing_directories(const char* dir /*= "."*/) const
{
  char cur_dir_name[255];
  int status;

  // Base directory
  status = mkdir(dir, 0755);
  if ((status == -1) && (errno != EEXIST))
  {
    err(EXIT_FAILURE, cur_dir_name, errno);
  }

  // Experimental setup
  sprintf(cur_dir_name, "%s/" EXP_S_DIR, dir);
  status = mkdir(cur_dir_name, 0755);
  if ((status == -1) && (errno != EEXIST))
  {
    err(EXIT_FAILURE, cur_dir_name, errno);
  }
  // Output profile
  sprintf(cur_dir_name, "%s/" OUT_P_DIR, dir);
  status = mkdir(cur_dir_name, 0755);
  if ((status == -1) && (errno != EEXIST))
  {
    err(EXIT_FAILURE, cur_dir_name, errno);
  }
  // World
  sprintf(cur_dir_name, "%s/" WORLD_DIR, dir);
  status = mkdir(cur_dir_name, 0755);
  if (status == -1 && errno != EEXIST)
  {
    err(EXIT_FAILURE, cur_dir_name, errno);
  }
  // Stats
  sprintf(cur_dir_name, "%s/" STATS_DIR, dir);
  status = mkdir(cur_dir_name, 0755);
  if (status == -1 && errno != EEXIST)
  {
    err(EXIT_FAILURE, cur_dir_name, errno);
  }
}

void ExpManager::open_backup_files(gzFile& exp_backup_file,
                                       gzFile& world_file,
                                       int64_t t,
                                       const char mode[3],
                                       const char* dir /*= "."*/
                                      #ifdef HAVE_MPI
                                       , int rank /*= -1*/
                                      #endif
                                       ) const
{
  assert(strcmp(mode, "w") == 0 or strcmp(mode, "r") == 0);

  // -------------------------------------------------------------------------
  // Generate backup file names for mandatory files.
  // -------------------------------------------------------------------------
  char exp_backup_file_name[255];
  char world_file_name[255];


#ifdef HAVE_MPI
  sprintf(world_file_name, "%s/" WORLD_FNAME_FORMAT, dir, t, rank);
  sprintf(exp_backup_file_name, "%s/" EXP_S_FNAME_FORMAT, dir, t, rank);
#else
  sprintf(world_file_name, "%s/" WORLD_FNAME_FORMAT, dir, t);
  sprintf(exp_backup_file_name, "%s/" EXP_S_FNAME_FORMAT, dir, t);
#endif

  // -------------------------------------------------------------------------
  // Open backup files
  // -------------------------------------------------------------------------
  exp_backup_file = gzopen(exp_backup_file_name, mode);
  world_file = gzopen(world_file_name, mode);


  // -------------------------------------------------------------------------
  // Check that files were correctly opened
  // -------------------------------------------------------------------------
  if (exp_backup_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n",
            __FILE__, __LINE__, exp_backup_file_name);
    exit(EXIT_FAILURE);
  }
  if (world_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n",
            __FILE__, __LINE__, world_file_name);
    exit(EXIT_FAILURE);
  }
}

void ExpManager::close_backup_files(gzFile& exp_backup_file,
                                        gzFile& world_file) const
{
  gzclose(exp_backup_file);
  gzclose(world_file);
}

void ExpManager::open_setup_files(
    gzFile& exp_s_file,
    gzFile& out_p_file,
    int64_t t,
    const char mode[3],
    const char* dir /*= "."*/
    #ifdef HAVE_MPI
    , int rank /*= -1*/
    #endif
    ) const
{
  // 1) Generate setup file names
  char exp_s_file_name[255];
  char out_p_file_name[255];

#ifdef HAVE_MPI
  sprintf(exp_s_file_name, "%s/" EXP_S_CONST_FNAME_FORMAT, dir, rank);
  sprintf(out_p_file_name, "%s/" OUT_P_FNAME_FORMAT, dir, rank);
#else
  sprintf(exp_s_file_name, "%s/" EXP_S_CONST_FNAME_FORMAT, dir);
  sprintf(out_p_file_name, "%s/" OUT_P_FNAME_FORMAT, dir);
#endif

  // 2) Open backup files
  exp_s_file = gzopen(exp_s_file_name, mode);
  out_p_file = gzopen(out_p_file_name, mode);

  // 3) Check that files were correctly opened
  if (exp_s_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n",
            __FILE__, __LINE__, exp_s_file_name);
    exit(EXIT_FAILURE);
  }
  if (out_p_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n",
            __FILE__, __LINE__, out_p_file_name);
    exit(EXIT_FAILURE);
  }
}

void ExpManager::close_setup_files(gzFile& exp_s_file,
                                       gzFile& out_p_file) const
{
  gzclose(exp_s_file);
  gzclose(out_p_file);
}

// ===========================================================================
//                               Non inline accessors
// ===========================================================================
Individual* ExpManager::indiv_by_id(int32_t id) const {
  return world_->indiv_by_id(id);
}

Individual* ExpManager::indiv_by_rank(int32_t rank) const {
  return world_->indiv_by_rank(rank);
}
Individual* ExpManager::indiv_by_position(int16_t x, int16_t y) const {
  return world_->indiv_at(x,y);
}

/**
 * Returns a list of all the individuals with their replication report
 */
std::list<std::pair<Individual*, ReplicationReport*>>
    ExpManager::indivs_annotated() const {
  std::list<std::pair<Individual*, ReplicationReport*>> annotated_list;
  for (const auto& indiv : indivs()) {
    annotated_list.emplace_back(indiv, tree() ?
        tree()->report_by_index(AeTime::time(), indiv->id()) : nullptr);
  }
  return annotated_list;
}


} // namespace aevol
