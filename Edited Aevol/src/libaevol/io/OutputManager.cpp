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
#include "OutputManager.h"

#include <err.h>
#include <errno.h>
#include <sys/stat.h>
#include <fstream>
#include <string>

#include "ExpManager.h"
#include "7/ExpManager_7.h"
#include "AeTime.h"

using std::string;
using std::endl;

namespace aevol {



//##############################################################################
//                                                                             #
//                           Class OutputManager                           #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
OutputManager::OutputManager(ExpManager * exp_m)
{
  exp_m_  = exp_m;
  stats_  = nullptr;
  tree_   = nullptr;
  light_tree_ = nullptr;
  dump_   = nullptr;
  compute_phen_contrib_by_GU_ = false;
  record_tree_ = false;
  record_light_tree_ = false;
  make_dumps_ = false;
  dump_step_ = 0;
  logs_  = new Logging();

  backup_step_ = -1;
}

// =================================================================
//                             Destructors
// =================================================================
OutputManager::~OutputManager()
{
  delete stats_;
  delete tree_;
  delete light_tree_;
  delete dump_;
  delete logs_;
}

// =================================================================
//                            Public Methods
// =================================================================
void OutputManager::InitStats() {
  stats_ = new Stats(exp_m_);
}

void OutputManager::WriteSetupFile(gzFile setup_file) const
{
  // Write the backup steps
  gzwrite(setup_file, &backup_step_,      sizeof(backup_step_));
  gzwrite(setup_file, &big_backup_step_,  sizeof(big_backup_step_));

  // Stats
  gzwrite(setup_file, &compute_phen_contrib_by_GU_,  sizeof(compute_phen_contrib_by_GU_));

  // Tree
  int8_t record_tree = record_tree_;
  gzwrite(setup_file, &record_tree, sizeof(record_tree));


    // LightTree
    int8_t record_light_tree = (int8_t) record_light_tree_;

    gzwrite(setup_file, &record_light_tree, sizeof(record_light_tree));

    // Tree step
  if (record_tree_)
  {
      int64_t tmp_tree_step = tree_->tree_step();
    gzwrite(setup_file, &tmp_tree_step, sizeof(tmp_tree_step));
  }

    if (record_light_tree)
    {
        int64_t tmp_tree_step = light_tree_->tree_step();
        gzwrite(setup_file, &tmp_tree_step, sizeof(tmp_tree_step));
    }

    //printf("Save light tree %d %d %d %d\n",record_light_tree_,record_tree_,burp,burp2);

    // Dumps
  int8_t make_dumps = make_dumps_;
  gzwrite(setup_file, &make_dumps,  sizeof(make_dumps));
  gzwrite(setup_file, &dump_step_,  sizeof(dump_step_));

  // Logs
  int8_t logs = logs_->logs();
  gzwrite(setup_file, &logs,  sizeof(logs));
}

void OutputManager::CopyStats(const std::string& outdir, int64_t time) const {
  stats_->CreateTmpFiles(time);
  stats_->MoveTmpFiles(outdir);
}

void OutputManager::load(gzFile setup_file, bool verbose, bool to_be_run)
{
  // Write the backup steps
  gzread(setup_file, &backup_step_,      sizeof(backup_step_));
  gzread(setup_file, &big_backup_step_,  sizeof(big_backup_step_));
  // printf("Backup step %d\n",backup_step_);
  // Stats
  #ifndef HAVE_MPI
  if (to_be_run)
  {
    delete stats_;
    stats_ = new Stats(exp_m_, AeTime::time());
  }
  #endif
  gzread(setup_file, &compute_phen_contrib_by_GU_,  sizeof(compute_phen_contrib_by_GU_));

  // Tree
  int8_t record_tree;
  gzread(setup_file, &record_tree, sizeof(record_tree));
  record_tree_ = record_tree;


  // LightTree
  int8_t record_light_tree;
  gzread(setup_file, &record_light_tree, sizeof(record_light_tree));
  record_light_tree_ = record_light_tree;

  // Tree step
  if (record_tree_)
  {
    int64_t tmp_tree_step;
    gzread(setup_file, &tmp_tree_step, sizeof(tmp_tree_step));

    tree_ = new Tree(exp_m_, tmp_tree_step);
  }


  if (record_light_tree_ && to_be_run) {
    int64_t tmp_tree_step;
    gzread(setup_file, &tmp_tree_step, sizeof(tmp_tree_step));

    light_tree_ = new LightTree(exp_m_,tmp_tree_step);
    light_tree_->init_tree(AeTime::time(), exp_m_->indivs());
  }

  // Dumps
  int8_t make_dumps;
  gzread(setup_file, &make_dumps, sizeof(make_dumps));
  make_dumps_ = make_dumps;
  gzread(setup_file, &dump_step_,  sizeof(dump_step_));
  if(make_dumps_ == true)
  {
    dump_ = new Dump(exp_m_);
  }

  // Logs
  int8_t logs;
  gzread(setup_file, &logs, sizeof(logs));
  if (to_be_run)
  {
    logs_->load(logs, AeTime::time());
  }
}

void OutputManager::write_current_generation_outputs(bool create) const
{

  // we use the variable t because of the parallelisation.
  int64_t t = AeTime::time();
  std::list<Individual*> indivs = exp_m_->indivs();

  //printf("Add indivs %d\n",indivs.size());

  //stats_->add_indivs(AeTime::time(), indivs);

//  SaveWorld* backup_world;
//  //JumpingMT* backup_prng;
//  if (t % backup_step_ == 0) {
//    backup_world = exp_m_->world()->make_save(exp_m_, indivs,exp_m_->world()->phenotypic_target_shared());
//  }

//#ifdef __OPENMP_TASK
//#pragma omp task depend(out: dep[d_LT_LT])
//#endif
  // LightTree
//  if (record_light_tree_ && t > 0) {
//    light_tree_->update_tree(t, nullptr);
//    if(t % backup_step_ == 0) {
//      //std::cout << "writing light tree for gen : " << t << '\n';
//      write_light_tree(t);
//    }
//  }

//#ifdef __OPENMP_TASK
//  #pragma omp task  depend(out: dep[d_S_S])
//{
//#endif
  stats_->write_current_generation_statistics(t);

//#ifdef __OPENMP_TASK
//#pragma omp task
//#endif
  if (record_tree_ &&
      t > 0 &&
      (t % tree_->tree_step() == 0)) {
    std::cout << "writing tree for gen : " << t << '\n';
    write_tree(t);
  }

//  if(t > 0 && ((t-1) % backup_step_ != 0)) {
//    stats_->delete_indivs(t-1);
//  }
//#ifdef __OPENMP_TASK
//  }
//#endif

//#ifdef __OPENMP_TASK
//  #pragma omp task  depend(out: dep[d_B_B])
//{
//#endif
//  if ((t-1) % backup_step_ == 0)
//    stats_->delete_indivs(t-1);

  // Write backup
  if (t % backup_step_ == 0) {
    std::cout << "writing backup for gen : " << t << '\n';
    stats_->flush();
    exp_m_->WriteDynamicFiles(create);
//    if (!exp_m_->check_simd() || (exp_m_->check_simd() && AeTime::time()==0))
//      exp_m_->WriteDynamicFiles(t, backup_world, create);

    WriteLastGenerFile(".", t);
    //delete backup_prng;
//    delete backup_world;
    }
//#ifdef __OPENMP_TASK
//  }
//#endif

//#ifdef __OPENMP_TASK
//  if(omp_get_num_threads() < 2) {
//  #pragma omp taskwait
//}
//#endif

  // Write dumps
  if (make_dumps_) {
    if(AeTime::time() % dump_step_ == 0) {
      dump_->write_current_generation_dump();
    }
  }
}

// TODO <david.parsons@inria.fr> we need an output_dir attribute in this class !
    void OutputManager::WriteLastGenerFile(const string& output_dir /* = "." */, int64_t gen /* = -1 */) const {
      if(gen < 0) {
        gen = AeTime::time();
      }
      std::ofstream last_gener_file(output_dir + "/" + LAST_GENER_FNAME,
                                    std::ofstream::out);
      if (last_gener_file.fail()) {
        Utils::ExitWithUsrMsg(string("could not open file ") + LAST_GENER_FNAME);
      }
      else {
        last_gener_file << gen << endl;
        last_gener_file.close();
      }
    }

// TODO <david.parsons@inria.fr> we need an input_dir attribute in this class !
int64_t OutputManager::last_gener() {
  int64_t time = 0;
  FILE* lg_file = fopen(LAST_GENER_FNAME, "r");
  if (lg_file != NULL) {
    if (fscanf(lg_file, "%" PRId64 "\n", &time) == EOF) {
      Utils::ExitWithDevMsg("failed to read last generation", __FILE__,
                            __LINE__);
    }
    fclose(lg_file);
  }
  return time;
}
    void OutputManager::flush() {
      stats_->flush();

      // Tree
      if (record_tree_ &&
          AeTime::time() > 0 &&
          (AeTime::time() % tree_->tree_step() != 0)) {
          if (!ExpManager_7::standalone_simd) {
              write_tree(AeTime::time());
          }
      }

      // LightTree
      if (record_light_tree_ &&
          AeTime::time() > 0 &&
          (AeTime::time() % backup_step_ != 0)) {
        write_light_tree(AeTime::time());
      }

      // Write backup
      if (AeTime::time() % backup_step_ != 0) {
        stats_->flush();
        if (!ExpManager_7::standalone_simd) {
            exp_m_->WriteDynamicFiles();
            WriteLastGenerFile();
        }
      }

      // Write dumps
      if (make_dumps_) {
        if(AeTime::time() % dump_step_ != 0) {
          dump_->write_current_generation_dump();
        }
      }
    }

// =================================================================
//                           Protected Methods
// =================================================================
    void OutputManager::write_tree(int64_t gen) const
    {
      // Create the tree directory if it doesn't exist

      // debug std::cout << "Writing regular tree file ...";

      int status;
      status = mkdir(TREE_DIR, 0755);
      if ((status == -1) && (errno != EEXIST))
      {
        err(EXIT_FAILURE, "Impossible to create the directory %s", TREE_DIR);
      }

      char tree_file_name[50];

        #ifdef HAVE_MPI
        sprintf(tree_file_name, "tree/tree_" TIMESTEP_FORMAT "_%d.ae", AeTime::time(),exp_m_->rank());
        #else
        sprintf(tree_file_name, "tree/tree_" TIMESTEP_FORMAT ".ae", AeTime::time());
        #endif
        
//      gzFile tree_file = gzopen( tree_file_name, "w" );
      // Write phylogenetic data (tree)
      tree_->write_to_tree_file(tree_file_name);

//      gzclose(tree_file);

      // debug std::cout << "OK" << '\n';
    }


    void OutputManager::write_light_tree(int64_t gen) const
    {
      // debug std::cout << "Writing light tree file ...";

      int status = mkdir(LIGHTTREE_DIR, 0755);
      if ((status == -1) && (errno != EEXIST))
      {
        err(EXIT_FAILURE, "Impossible to create the directory %s", LIGHTTREE_DIR);
      }

      static char branches_file_name[50];
      sprintf(branches_file_name, "lightTree/tree_branches.ae");

      char trunc_file_name[50];
      sprintf(trunc_file_name, "lightTree/tree_trunc" TIMESTEP_FORMAT ".ae", gen);

      gzFile trunc_file = gzopen( trunc_file_name, "w" );
      gzFile branches_file = gzopen( branches_file_name, "w" );
      // Write phylogenetic data (tree)
      light_tree_->write_to_tree_file(gen, trunc_file, branches_file);

      gzclose(trunc_file);
      gzclose(branches_file);

      light_tree_->write_tree(gen);

      // debug std::cout << "OK" << '\n';
    }

// =================================================================
//                          Non inline accessors
// =================================================================
} // namespace aevol
