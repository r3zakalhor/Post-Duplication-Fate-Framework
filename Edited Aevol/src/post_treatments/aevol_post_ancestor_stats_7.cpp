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

// The input file is produced by the lineage post-treatment, please refer to it
// for e.g. the file format/content

// ============================================================================
//                                   Includes
// ============================================================================
#include <cinttypes>
#include <getopt.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <zlib.h>
#include <err.h>
#include <cerrno>
#include <sys/stat.h>
#include <unistd.h>
#include <list>
#include <iostream>

#include <cstdint>
#include <fstream>
#include <limits>
#include <string>

#include "aevol.h"
#ifdef _OPENMP
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0;}
inline omp_int_t omp_get_max_threads() { return 1;}
#endif
#include "Stats_7.h"
#include "List_Metadata.h"
using namespace aevol;

// Helper functions
void interpret_cmd_line_options(int argc, char* argv[]);
void print_help(char* prog_path);

// Command-line option variables
static char* lineage_file_name = nullptr;
static bool verbose = false;
static bool full_check = false;
static bool trace_mutations = false;

#ifdef _OPENMP
static bool run_in_parallel = false;
#endif

int main(int argc, char* argv[]) {
  interpret_cmd_line_options(argc, argv);

  printf("\n"
         "WARNING : Parameter change during simulation is not managed in general.\n"
         "          Only changes in environmental target done with aevol_modify are handled.\n"
         "\n");

  // =======================
  //  Open the lineage file
  // =======================
  gzFile lineage_file = gzopen(lineage_file_name, "r");
  if (lineage_file == Z_NULL) {
    fprintf(stderr, "ERROR : Could not read the lineage file %s\n", lineage_file_name);
    delete [] lineage_file_name;
    exit(EXIT_FAILURE);
  }

  delete [] lineage_file_name;

  int64_t t0 = 0;
  int64_t t_end = 0;
  int32_t final_indiv_index = 0;
  int32_t final_indiv_rank  = 0;

  gzread(lineage_file, &t0, sizeof(t0));
  gzread(lineage_file, &t_end, sizeof(t_end));
  gzread(lineage_file, &final_indiv_index, sizeof(final_indiv_index));
  gzread(lineage_file, &final_indiv_rank,  sizeof(final_indiv_rank));

  if (verbose) {
    printf("\n\n""===============================================================================\n");
    printf(" Statistics of the ancestors of indiv. %" PRId32
           " (rank %" PRId32 ") from time %" PRId64 " to %" PRId64 "\n",
           final_indiv_index, final_indiv_rank, t0, t_end);
    printf("================================================================================\n");
  }


#ifdef __REGUL
    int nb_iteration = 4;
#endif
  // =============================
  //  Open the experiment manager
  // =============================
  ExpManager_7::standalone_simd = true;
  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(t0, true, false);

  // The current version doesn't allow for phenotypic variation nor for
  // different phenotypic targets among the grid
  if (not exp_manager->world()->phenotypic_target_shared())
    Utils::ExitWithUsrMsg("sorry, ancestor stats has not yet been implemented "
                              "for per grid-cell phenotypic target");
  auto phenotypicTargetHandler =
      exp_manager->world()->phenotypic_target_handler();
    /*if (not (phenotypicTargetHandler->var_method() == NO_VAR))
      Utils::ExitWithUsrMsg("sorry, ancestor stats has not yet been implemented "
                                "for variable phenotypic targets");*/

  int64_t backup_step = exp_manager->backup_step();


  // =========================
  //  Open the output file(s)
  // =========================
  // Create missing directories
  int status;
  status = mkdir("stats/ancestor_stats/", 0755);
  if ((status == -1) && (errno != EEXIST)) {
    err(EXIT_FAILURE, "stats/ancestor_stats/");
  }

  // Open main output files (uses the Stats utility class)
  char prefix[255] = "stats/ancestor_stats";
  char postfix[255];
  snprintf(postfix, 255,
      "-b" TIMESTEP_FORMAT "-e" TIMESTEP_FORMAT "-i%" PRId32 "-r%" PRId32,
      t0, t_end, final_indiv_index, final_indiv_rank);
  bool best_indiv_only = true;
  bool addition_old_stats = false;
  bool delete_old_stats = true;

  Stats_7* stats_anc = new Stats_7(exp_manager->exp_m_7_, 
                                AeTime::time(), true, prefix,true);
  // Optional additional outputs

  // Open optional output files
  FILE* fixed_mutations_file = nullptr;
  if (trace_mutations) {
    char fixed_mutations_file_name[255];
    snprintf(fixed_mutations_file_name, 255,
             "stats/ancestor_stats/fixedmut-b" TIMESTEP_FORMAT
             "-e" TIMESTEP_FORMAT "-i%" PRId32 "-r%" PRId32 ".out",
             t0, t_end, final_indiv_index, final_indiv_rank);
    fixed_mutations_file = fopen(fixed_mutations_file_name, "w");
    if (fixed_mutations_file == nullptr) {
      Utils::ExitWithUsrMsg(std::string("Could not create the output file ") +
                            fixed_mutations_file_name);
    }

    // Write the header
    fprintf(fixed_mutations_file, "# #################################################################\n");
    fprintf(fixed_mutations_file, "#  Mutations in the lineage of the best indiv at generation %" PRId64 "\n", t_end);
    fprintf(fixed_mutations_file, "# #################################################################\n");
    fprintf(fixed_mutations_file, "#  1.  Generation       (mut. occurred when producing the indiv. of this generation)\n");
    fprintf(fixed_mutations_file, "#  2.  Genetic unit     (which underwent the mutation, 0 = chromosome) \n");
    fprintf(fixed_mutations_file, "#  3.  Mutation type    (0: switch, 1: smallins, 2: smalldel, 3:dupl, 4: del, 5:trans, 6:inv, 7:insert, 8:ins_HT, 9:repl_HT) \n");
    fprintf(fixed_mutations_file, "#  4.  pos_0            (position for the small events, begin_segment for the rearrangements, begin_segment of the inserted segment for ins_HT, begin_segment of replaced segment for repl_HT) \n");
    fprintf(fixed_mutations_file, "#  5.  pos_1            (-1 for the small events, end_segment for the rearrangements, end_segment of the inserted segment for ins_HT, begin_segment of donor segment for repl_HT) \n");
    fprintf(fixed_mutations_file, "#  6.  pos_2            (reinsertion point for duplic., cutting point in segment for transloc., insertion point in the receiver for ins_HT, end_segment of the replaced segment for repl_HT, -1 for other events)\n");
    fprintf(fixed_mutations_file, "#  7.  pos_3            (reinsertion point for transloc., breakpoint in the donor for ins_HT, end_segment of the donor segment for repl_HT, -1 for other events)\n");
    fprintf(fixed_mutations_file, "#  8.  invert           (transloc, was the segment inverted (0/1)?, sense of insertion for ins_HT (0=DIRECT, 1=INDIRECT), sense of the donor segment for repl_HT (0=DIRECT, 1=INDIRECT),-1 for other events)\n");
    fprintf(fixed_mutations_file, "#  9.  align_score      (score that was needed for the rearrangement to occur, score of the first alignment for ins_HT and repl_HT)\n");
    fprintf(fixed_mutations_file, "#  10. align_score2     (score for the reinsertion for transloc, score of the second alignment for ins_HT and repl_HT)\n");
    fprintf(fixed_mutations_file, "#  11. seg_len          (segment length for rearrangement, donor segment length for ins_HT and repl_HT)\n");
    fprintf(fixed_mutations_file, "#  12. repl_seg_len     (replaced segment length for repl_HT, -1 for the others)\n");
    fprintf(fixed_mutations_file, "#  13. GU_length        (before the event)\n");
    fprintf(fixed_mutations_file, "#  14. Impact of the mutation on the metabolic error (negative value = smaller gap after = beneficial mutation) \n");
    fprintf(fixed_mutations_file, "#  15. Impact of the mutation on fitness (positive value = higher fitness after = beneficial mutation) \n");
    fprintf(fixed_mutations_file, "#  16. Selection coefficient of the mutation (s > 1 = higher fitness after = beneficial mutation) \n");
    fprintf(fixed_mutations_file, "####################################################################################################################\n");
    fprintf(fixed_mutations_file, "#\n");
    fprintf(fixed_mutations_file, "# Header for R\n");
    fprintf(fixed_mutations_file, "gener gen_unit mut_type pos_0 pos_1 pos_2 pos_3 invert align_score align_score_2 seg_len repl_seg_len GU_len impact_on_gap impact_on_fitness sel_coeff\n");

  }


  // Init Factory (Fuzzy/Dna)
  DnaFactory* dna_factory_ = new DnaFactory(DnaFactory_Policy::FIRSTFIT,32,5000);
  FuzzyFactory_7* fuzzy_factory_ = new FuzzyFactory_7(exp_manager->exp_s()->get_fuzzy_flavor(),exp_manager->nb_indivs()*4,
                        exp_manager->world()->phenotypic_target_handler()->sampling());

  // ==================================================
  //  Prepare the initial ancestor and write its stats
  // ==================================================
  GridCell* grid_cell = new GridCell(lineage_file, exp_manager, nullptr);

  double w_max = exp_manager->best_indiv()->w_max();
  double selection_pressure = exp_manager->selection_pressure();

  // TODO Transform to Individual_7
  Individual_7* indiv = new Individual_7(exp_manager, grid_cell->individual()->w_max(),dna_factory_,fuzzy_factory_);
  indiv->dna_ = dna_factory_->get_dna(grid_cell->individual()->genetic_unit_seq_length(0));
  indiv->dna_->set_indiv(grid_cell->individual()->genetic_unit(0).dna(),dna_factory_);
  indiv->dna_->set_indiv(indiv);
  indiv->indiv_id = 0;
  indiv->parent_id = 0;

#ifdef __REGUL
    std::ofstream fitmeta;
    fitmeta.open("ancestor_composed_fitness.csv",std::ofstream::trunc);
    fitmeta<<"Generation,EnvID,Composed,Fitness"<<std::endl;

    SIMD_PhenotypicTargetHandler_R* phenotypic_target_handler_ = new SIMD_PhenotypicTargetHandler_R(
      std::dynamic_pointer_cast<PhenotypicTargetHandler_R>(exp_manager->world()->phenotypic_target_handler()),
      exp_manager->exp_s(),fuzzy_factory_,exp_manager->check_simd());
    phenotypic_target_handler_->ApplyVariation();
    
    exp_manager->exp_m_7_->evaluate(indiv,w_max,selection_pressure,phenotypic_target_handler_);

    printf("Initial fitness     = %e\n", indiv->fitness);
    printf("Initial genome size = %" PRId32 "\n", indiv->dna_->length());

    double base_metaerror = 0;
    double base_fitness = 0;

    int32_t max_it_env = nb_iteration > phenotypic_target_handler_->nb_env_ ? nb_iteration : phenotypic_target_handler_->nb_env_;


    SIMD_PhenotypicTargetHandler_R** pth_array_0;
    pth_array_0 = new SIMD_PhenotypicTargetHandler_R*[max_it_env];
    for (int i = 0; i < max_it_env; i++) {
      pth_array_0[i] = new SIMD_PhenotypicTargetHandler_R(phenotypic_target_handler_,fuzzy_factory_,exp_manager);

      pth_array_0[i]->var_prng_ = std::make_shared<JumpingMT>(phenotypic_target_handler_->var_prng_->random(100000000));
    }

    #pragma omp parallel for
    for (int i = 0; i < nb_iteration; i++) {
        pth_array_0[i]->ApplyVariation();
        Individual_7* cloned = new Individual_7(exp_manager,indiv,dna_factory_,fuzzy_factory_,true);
        
        exp_manager->exp_m_7_->evaluate(cloned,w_max,selection_pressure,pth_array_0[i]);


#pragma omp atomic
        base_metaerror += cloned->metaerror;

#pragma omp atomic
        base_fitness += cloned->fitness;

#pragma omp critical
        fitmeta << t0 << ",RANDOM," << i << "," << cloned->metaerror << ","
                << cloned->fitness << std::endl;

        //printf("Iteration RANDOM at %d : %d/%d : %lf %e\n", t, i, nb_iteration,
        //       cloned->dist_to_target_by_feature(METABOLISM), cloned->fitness());
        delete cloned;
    }

    indiv->fitness = base_fitness / (double) nb_iteration;
    indiv->metaerror = base_metaerror / (double) nb_iteration;

    exp_manager->exp_m_7_->write_stat(stats_anc,indiv,AeTime::time(),true);

    // Single env
//     for (int i = 0; i < phenotypic_target_handler_->nb_env_; i++) {
//         pth_array_0[i]->set_single_env(i);
//         Individual_7* cloned = new Individual_7(exp_manager,indiv,dna_factory_,fuzzy_factory_,true);
        
//         exp_manager->exp_m_7_->evaluate(cloned,w_max,selection_pressure,pth_array_0[i]);


// #pragma omp atomic
//         base_metaerror += cloned->metaerror;

// #pragma omp atomic
//         base_fitness += cloned->fitness;

// #pragma omp critical
//         fitmeta << t0 << ",TARGET," << i << "," << cloned->metaerror << ","
//                 << cloned->fitness << std::endl;

//         //printf("Iteration RANDOM at %d : %d/%d : %lf %e\n", t, i, nb_iteration,
//         //       cloned->dist_to_target_by_feature(METABOLISM), cloned->fitness());
//         delete cloned;
//     }
#else
  exp_manager->exp_m_7_->evaluate(indiv,w_max,selection_pressure);
  exp_manager->exp_m_7_->write_stat(stats_anc,indiv,AeTime::time(),true);
#endif

#ifdef __REGUL
  for (int i = 0; i < max_it_env; i++) {
      delete pth_array_0[i];
  }

  delete [] pth_array_0;
#endif
  if (verbose) {
    printf("Initial fitness     = %f\n", indiv->fitness);
    printf("Initial genome size = %" PRId32 "\n", indiv->dna_->length());
  }

  // ==========================================================================
  //  Replay the mutations to get the successive ancestors and analyze them
  // ==========================================================================
  ReplicationReport* rep = nullptr;
  int32_t index;
  ExpManager* exp_manager_backup = nullptr;
  char mut_descr_string[255];

  bool check_now = false;

  aevol::AeTime::plusplus();
  #pragma omp parallel
  {
  
  while (time() <= t_end)
  {
    #pragma omp single
    {
    rep = new ReplicationReport(lineage_file, nullptr);
    index = rep->id(); // who we are building...

    // Check now?
    check_now = time() == t_end ||
        (full_check && Utils::mod(time(), backup_step) == 0);

    if (verbose)
        printf("Rebuilding ancestor at generation %" PRId64
            " (index %" PRId32 ")...", time(), index);

    indiv->reset_metadata();
    exp_manager->exp_m_7_->evaluate(indiv,w_max,selection_pressure);

    // 2) Replay replication (create current individual's child)
    GeneticUnit* stored_gen_unit = nullptr;
    Individual* stored_indiv = nullptr;

    if (check_now)
    {
      exp_manager_backup = new ExpManager();
      exp_manager_backup->load(time(), true, true);
        // Copy the ancestor from the backup

        int l_x = index/exp_manager_backup->world()->height();
        int l_y = index%exp_manager_backup->world()->height();
        stored_indiv = exp_manager_backup->world()->indiv_at(l_x,l_y);

        stored_gen_unit = &(stored_indiv->genetic_unit_nonconst(0));
    }

    // For each genetic unit, replay the replication (undergo all mutations)
    // TODO <david.parsons@inria.fr> disabled for multiple GUs
    const auto& dnarep = rep->dna_replic_report();

    // TODO(dpa) The following 3 for loops should be factorized.
    // However, this is not as easy as it sounds :-D
    // see std::list::splice
    // TG: done :) !

    dnarep.iter_muts([&](const auto& mut) {
      int32_t unitlen_before = 0;
      double metabolic_error_before = 0.0;
      double fitness_before = 0.0;

      if (trace_mutations) {
        // Store initial values before the mutation
        metabolic_error_before = indiv->metaerror;
        fitness_before         = indiv->fitness;
        unitlen_before         = indiv->dna_->length();
      }

      // Apply mutation
      indiv->dna_->undergo_mutation(*mut);

      if (trace_mutations) {
        indiv->reset_metadata();
        exp_manager->exp_m_7_->evaluate(indiv,w_max,selection_pressure);

        // Compute the metabolic impact of the mutation
        double impact_on_metabolic_error =
            indiv->metaerror -
            metabolic_error_before;

        double impact_on_fitness = indiv->fitness - fitness_before;

        double selection_coefficient = indiv->fitness / fitness_before - 1.0;

        mut->generic_description_string(mut_descr_string);
        fprintf(fixed_mutations_file,
                "%" PRId64 " %" PRId32 " %s %" PRId32 " %.15e %.15e %.15e\n",
                time(), 0, mut_descr_string, unitlen_before,
                impact_on_metabolic_error, impact_on_fitness,
                selection_coefficient);
      }
    });

    if (check_now) {
      if (verbose) {
        printf("Checking the sequence of the unit...");
        fflush(NULL);
      }

      char * str1 = new char[indiv->dna_->length() + 1];
      memcpy(str1, indiv->dna_->data_, \
             indiv->dna_->length()*sizeof(char));
      str1[indiv->dna_->length()] = '\0';

      char * str2 = new char[(stored_gen_unit->dna())->length() + 1];
      memcpy(str2, (stored_gen_unit->dna())->data(),
             (stored_gen_unit->dna())->length()*sizeof(char));
      str2[(stored_gen_unit->dna())->length()] = '\0';

      if (strncmp(str1, str2, stored_gen_unit->dna()->length()) == 0) {
        if (verbose)
          printf(" OK\n");
      }
      else {
        if (verbose) printf(" ERROR !\n");
        fprintf(stderr, "Error: the rebuilt genetic unit is not the same as \n");
        fprintf(stderr, "the one saved at generation %" PRId64 "... ", time());
        fprintf(stderr, "Rebuilt unit : %" PRId32 " bp\n %s\n", (int32_t)strlen(str1), str1);
        fprintf(stderr, "Stored unit  : %" PRId32 " bp\n %s\n", (int32_t)strlen(str2), str2);

        delete [] str1;
        delete [] str2;
        gzclose(lineage_file);
        if (trace_mutations)
          fclose(fixed_mutations_file);
        delete exp_manager_backup;
        delete exp_manager;
        exit(EXIT_FAILURE);
      }

      delete [] str1;
      delete [] str2;
    }

#ifdef __REGUL
    if (AeTime::time() % 100 == 0) {
#endif
    // 3) All the mutations have been replayed, we can now evaluate the new individual
    Individual_7* cloned_generation = new Individual_7(exp_manager, indiv, dna_factory_,fuzzy_factory_);
    int32_t current_generation = AeTime::time();
    
    #ifdef __REGUL
    #pragma omp task
    #else
    #pragma omp task
    #endif
    {
#ifdef __REGUL
      cloned_generation->reset_metadata();
      SIMD_PhenotypicTargetHandler_R** pth_array;

      
      #pragma omp critical
      {
	int32_t max_it_env = nb_iteration > phenotypic_target_handler_->nb_env_ ? nb_iteration : phenotypic_target_handler_->nb_env_;

        phenotypic_target_handler_->ApplyVariation();
        pth_array = new SIMD_PhenotypicTargetHandler_R*[max_it_env];
        for (int i = 0; i < max_it_env; i++) {
          pth_array[i] = new SIMD_PhenotypicTargetHandler_R(phenotypic_target_handler_,fuzzy_factory_,exp_manager);

          pth_array[i]->var_prng_ = std::make_shared<JumpingMT>(phenotypic_target_handler_->var_prng_->random(100000000));
          pth_array[i]->ApplyVariation();
	}

      }

      exp_manager->exp_m_7_->evaluate(cloned_generation,w_max,selection_pressure,pth_array[0]);

        double base_metaerror = 0;
        double base_fitness = 0;
        //#pragma omp parallel for
        for (int i = 0; i < nb_iteration; i++) {
            //#pragma omp critical
            pth_array[i]->ApplyVariation();
            
            Individual_7* cloned = new Individual_7(exp_manager,cloned_generation,dna_factory_,fuzzy_factory_,true);
            
            exp_manager->exp_m_7_->evaluate(cloned,w_max,selection_pressure,pth_array[i]);


    //#pragma omp atomic
            base_metaerror += cloned->metaerror;

    //#pragma omp atomic
            base_fitness += cloned->fitness;

    //#pragma omp critical
            // fitmeta << current_generation << ",RANDOM," << i << "," << cloned->metaerror << ","
            //         << cloned->fitness << std::endl;

            //printf("Iteration RANDOM at %d : %d/%d : %lf %e\n", current_generation, i, nb_iteration,
            //       cloned->metaerror, cloned->fitness);
            delete cloned;
        }

        cloned_generation->fitness = base_fitness / (double) nb_iteration;
        cloned_generation->metaerror = base_metaerror / (double) nb_iteration;

      #pragma omp critical
      {
        printf("Write stats for generation %d\n",current_generation);
        exp_manager->exp_m_7_->write_stat(stats_anc,cloned_generation,current_generation,true);
      }
    

    // Single env
//     for (int i = 0; i < phenotypic_target_handler_->nb_env_; i++) {
//         pth_array[i]->set_single_env(i);
//         Individual_7* cloned = new Individual_7(exp_manager,indiv,dna_factory_,fuzzy_factory_,true);
        
//         exp_manager->exp_m_7_->evaluate(cloned,w_max,selection_pressure,pth_array[i]);


// #pragma omp atomic
//         base_metaerror += cloned->metaerror;

// #pragma omp atomic
//         base_fitness += cloned->fitness;

// #pragma omp critical
//         fitmeta << current_generation << ",TARGET," << i << "," << cloned->metaerror << ","
//                 << cloned->fitness << std::endl;

//         //printf("Iteration RANDOM at %d : %d/%d : %lf %e\n", t, i, nb_iteration,
//         //       cloned->dist_to_target_by_feature(METABOLISM), cloned->fitness());
//         delete cloned;
//     }

#else
    cloned_generation->reset_metadata();
    exp_manager->exp_m_7_->evaluate(cloned_generation,w_max,selection_pressure);
    
    #pragma omp critical
    exp_manager->exp_m_7_->write_stat(stats_anc,cloned_generation,current_generation,true);
#endif


#ifdef __REGUL
  for (int i = 0; i < max_it_env; i++) {
      delete pth_array[i];
  }

  delete [] pth_array;
#endif

    delete cloned_generation;
    }
#ifdef __REGUL
    }
#endif

    if (verbose) printf(" OK\n");

    delete rep;

    if (check_now) {
      delete exp_manager_backup;
    }

    aevol::AeTime::plusplus();
  }
  }

  printf("%d -- LOOP finished\n",omp_get_thread_num());
  
  printf("%d -- Finished\n",omp_get_thread_num());
  }

  gzclose(lineage_file);

  delete grid_cell;
  delete exp_manager;
  delete stats_anc;

  return EXIT_SUCCESS;
}





void interpret_cmd_line_options(int argc, char* argv[]) {
  // =====================
  //  Parse command line
  // =====================
  const char * short_options = "p:hVFvM";
  static struct option long_options[] = {
    {"help",                no_argument, NULL, 'h'},
    {"version",             no_argument, NULL, 'V'},
    {"full-check",          no_argument, NULL, 'F'},
    {"trace-mutations",     no_argument, NULL, 'M'},
    {"verbose",             no_argument, NULL, 'v'},
    {"parallel",      required_argument, nullptr, 'p'},
    {0, 0, 0, 0}
  };

  int option;
  while((option = getopt_long(argc, argv, short_options,
                              long_options, nullptr)) != -1) {
    switch(option) {
      case 'h':
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      case 'V':
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      case 'v':
        verbose = true;
        break;
      case 'F':
        full_check = true;
        break;
      case 'M':
        trace_mutations = true;
        break;
      case 'p' : {
        #ifdef _OPENMP
        run_in_parallel = true;
          if (atoi(optarg) > 0) {
            omp_set_num_threads(atoi(optarg));
          } else
            omp_set_num_threads(1);
        #endif
        break;
      }
      default:
        // An error message is printed in getopt_long, we just need to exit
        exit(EXIT_FAILURE);
    }
  }
  
  #ifdef _OPENMP
  if (not run_in_parallel)
    omp_set_num_threads(1);
  #endif
  // There should be only one remaining arg: the lineage file
  if (optind != argc - 1) {
    Utils::ExitWithUsrMsg("please specify a lineage file");
  }

  lineage_file_name = new char[strlen(argv[optind]) + 1];
  sprintf(lineage_file_name, "%s", argv[optind]);
}

void print_help(char* prog_path) {
  // Get the program file-name in prog_name (strip prog_path of the path)
  char* prog_name; // No new, it will point to somewhere inside prog_path
  if ((prog_name = strrchr(prog_path, '/'))) {
    prog_name++;
  }
  else {
    prog_name = prog_path;
  }

  printf("******************************************************************************\n");
  printf("*                                                                            *\n");
  printf("*                        aevol - Artificial Evolution                        *\n");
  printf("*                                                                            *\n");
  printf("* Aevol is a simulation platform that allows one to let populations of       *\n");
  printf("* digital organisms evolve in different conditions and study experimentally  *\n");
  printf("* the mechanisms responsible for the structuration of the genome and the     *\n");
  printf("* transcriptome.                                                             *\n");
  printf("*                                                                            *\n");
  printf("******************************************************************************\n");
  printf("\n");
  printf("%s: create an experiment with setup as specified in PARAM_FILE.\n",
  prog_name);
  printf("\n");
  printf("Usage : %s -h or --help\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("   or : %s LINEAGE_FILE [-FMv]\n",
  prog_name);
  printf("\nOptions\n");
  printf("  -h, --help\n\tprint this help, then exit\n");
  printf("  -V, --version\n\tprint version number, then exit\n");
  printf("  -F, --full-check\n");
  printf("\tperform genome checks whenever possible\n");
  printf("  -M, --trace-mutations\n");
  printf("\toutputs the fixed mutations (in a separate file)\n");
  printf("  -v, --verbose\n\tbe verbose\n");
  printf("  -p, --parallel NB_THREADS\n");
  printf("\trun on NB_THREADS threads (use -1 for system default)\n");
}
