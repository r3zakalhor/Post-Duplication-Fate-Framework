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

// ============================================================================
//                                   Includes
// ============================================================================
#include "aevol.h"

#include <cerrno>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <getopt.h>
#include <list>
#include <queue>
#include <sys/stat.h>
#include <zlib.h>
//#ifdef HAVE_MPI
#include "ExpManager_7.h"
//#endif

enum check_type { FULL_CHECK = 0, LIGHT_CHECK = 1, NO_CHECK = 2 };

using namespace aevol;

class MullerStruct {
 public:
  MullerStruct(uint64_t l_strain_id,
               uint64_t l_parent_strain_id,
               uint64_t l_first_generation) {
    strain_id           = l_strain_id;
    parent_strain_id    = l_parent_strain_id;
    since_nb_generation = 1;
    first_generation    = l_first_generation;
  };

  uint64_t strain_id        = 0;
  uint64_t parent_strain_id = 0;
  std::list<uint64_t> children_strain_id;
  MullerStruct* parent_mstruct = nullptr;

  uint64_t since_nb_generation = 0;
  uint64_t first_generation    = 0;
  std::map<uint64_t, uint64_t> population_size_at;

  bool previously_output = false;
  bool to_delete         = false;

  std::list<int64_t> previous_generation_indexes;
  std::list<int64_t> current_generation_indexes;

  static uint64_t next_strain_id;
};

// Helper functions
void print_help(char* prog_path);
void interpret_cmd_line_options(int argc, char* argv[]);

// Command-line option variables
static bool verbose  = false;
static int64_t t0    = 0;
static int64_t t_end = -1;
static char tree_file_name[255];  // TODO(dpa) remove magic number
static uint64_t min_nb_generation = 0;
static double min_pop_percent     = 0;
static int64_t recording_end      = 0;
static int64_t recording_start    = 0;

uint64_t MullerStruct::next_strain_id = 2;

int main(int argc, char** argv) {

  // The output file (lineage.ae) contains the following information:
  // You may check that this information is up-to-date by searching
  // "lineage_file" in this source file
  //
  // - t0
  // - t_end
  // - final individual index
  // - final individual rank
  // - initial_ancestor information (including its genome)
  // - replication report of ancestor at generation t0+1
  // - replication report of ancestor at generation t0+2
  // - replication report of ancestor at generation t0+3
  // - ...
  // - replication report of ancestor at generation t_end

  interpret_cmd_line_options(argc, argv);

  printf("\n  WARNING : Parameter change in the middle of a simulation is not "
         "managed.\n");
  printf("min pop percent:  %f min nb gens:  %ld \n", min_pop_percent,
         min_nb_generation);
  

  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(t_end, true, false
#ifdef HAVE_MPI
                    ,
                    current_rank
#endif
  );

#ifdef HAVE_MPI
  exp_manager->exp_m_7_ = new ExpManager_7(exp_manager);
#endif

  // Check that the tree was recorded
  if (not exp_manager->record_tree()) {
    Utils::ExitWithUsrMsg("The phylogenetic tree wasn't recorded during "
                          "evolution, could not reconstruct the lineage");
  }

  int64_t tree_step = exp_manager->tree_step();
// Load the simulation
  //
  recording_end   = t_end;
  recording_start = t0;

  if (t0 % tree_step != 0) {
    t0 = t0 - t0 % tree_step;
  }

  if (t_end % tree_step != 0) {
    t_end = t_end - t_end % tree_step + tree_step;
  }

  if (verbose) {
    printf("Tree Step: %ld ", tree_step);
    printf("Recording start: %ld ", recording_start);
    printf("Recording end: %ld \n", recording_end);
  }

// The tree
#ifdef HAVE_MPI
  Tree** tree = new Tree*[exp_manager->nb_rank()];
#else
  Tree* tree = nullptr;
#endif

  // Muller plot data structure
  std::list<MullerStruct*> muller_list;

  // Muller plot data output file
  std::ofstream mullerfile;
  mullerfile.open("stats/mullerplot.csv", std::ofstream::trunc);
  mullerfile << "Generation,Identity,Population,Individuals" << std::endl;

  std::ofstream mullerphenofile;
  mullerphenofile.open("stats/mullerphenoplot.csv", std::ofstream::trunc);
  mullerphenofile << "Parent,Identity" << std::endl;
  // =========================
  //  Load the last tree file
  // =========================

  if (verbose) {
    printf("\n\n");
    printf("====================================\n");
    printf(" Loading the last tree file ... ");
    fflush(stdout);
  }

  // Example for ae_common::rec_params->tree_step() == 100 :
  //
  // tree_000100.ae ==>  timesteps 1 to 100.
  // tree_000200.ae ==>  timesteps 101 to 200.
  // tree_000300.ae ==>  timesteps 201 to 300.
  // etc.
  //

#ifdef HAVE_MPI
  for (int32_t current_rank = 0; current_rank < exp_manager->nb_rank();
       current_rank++) {
    exp_manager->exp_m_7_->setRank(current_rank);

    sprintf(tree_file_name, "tree/tree_" TIMESTEP_FORMAT "_%d.ae",
            t0 + tree_step, current_rank);
    tree[current_rank] = new Tree(exp_manager, tree_file_name);

#else
  sprintf(tree_file_name, "tree/tree_" TIMESTEP_FORMAT ".ae", t0 + tree_step);
  tree = new Tree(exp_manager, tree_file_name);
#endif

    for (int64_t index = 0; index < exp_manager->nb_indivs(); index++) {
      MullerStruct* mstruct =
          new MullerStruct(MullerStruct::next_strain_id++, 1, t0);

#ifdef HAVE_MPI
      int32_t x = indices / exp_manager->grid_height();
      int32_t y = indices % exp_manager->grid_height();

      int32_t global_x = exp_manager->exp_m_7_->localXtoGlobalX(x);
      int32_t global_y = exp_manager->exp_m_7_->localYtoGlobalY(y);
      int64_t global_id =
          global_x * exp_manager->exp_s()->global_grid_height() + global_y;

      mstruct->previous_generation_indexes.push_back(global_id);
#else
    mstruct->previous_generation_indexes.push_back(index);
#endif

      mstruct->population_size_at[t0] = 1;
      muller_list.push_back(mstruct);
    }

#ifdef HAVE_MPI
  }
#endif

  if (verbose) {
    printf("OK\n");
    printf("====================================\n");
  }

  // ===================================================
  //  Retrieve the replication reports of the ancestors
  // ===================================================

  if (verbose) {
    printf("\n\n\n");
    printf("==================================================================="
           "===\n");
    printf(" Parsing tree files to retrieve the muller plot data... \n");
    printf("==================================================================="
           "===\n");
  }

  // For each generation (going backwards), retrieve the index of the parent and
  // the corresponding replication report

  for (int64_t i = recording_start; i < recording_end - recording_start; i++) {
    // Where are we in time...
    int64_t t = t0 + i + 1;

    // We want to fill reports[i], that is to say, how the ancestor
    // at generation begin_gener + i + 1  was created
    if (verbose)
      printf("Getting the replication report for the ancestor at generation "
             "%" PRId64 " \n",
             t);

    // If we've exhausted the current tree file, load the next one
    if ((Utils::mod(t, tree_step) == 0) && (t + tree_step < t_end)) {
// Change the tree file
#ifdef HAVE_MPI
      for (int32_t current_rank = 0; current_rank < exp_manager->nb_rank();
           current_rank++) {
        delete tree[current_rank];

        exp_manager->exp_m_7_->setRank(current_rank);

        sprintf(tree_file_name, "tree/tree_" TIMESTEP_FORMAT "_%d.ae",
                t + tree_step, current_rank);
        tree[current_rank] = new Tree(exp_manager, tree_file_name);
      }
#else
      delete tree;

      sprintf(tree_file_name, "tree/tree_" TIMESTEP_FORMAT ".ae",
              t + tree_step);
      tree                    = new Tree(exp_manager, tree_file_name);
#endif
    }

#ifdef HAVE_MPI
    for (int32_t current_rank = 0; current_rank < exp_manager->nb_rank();
         current_rank++) {
      exp_manager->exp_m_7_->setRank(current_rank);
#endif

      for (int64_t index = 0; index < exp_manager->nb_indivs(); index++) {
#ifdef HAVE_MPI
        ReplicationReport* repl = tree[current_rank]->report_by_index(t, index);
#else
      ReplicationReport* repl = tree->report_by_index(t, index);
#endif

        //need to get the fitness of these mutations
        int32_t nb_mut = repl->nb(H_T) + repl->nb(REARR) + repl->nb(S_MUT);

        if (nb_mut > 0) {
          uint64_t parent_strain_id = 0;

#ifdef HAVE_MPI
          int32_t x = repl->parent_id() / exp_manager->grid_height();
          int32_t y = repl->parent_id() % exp_manager->grid_height();

          int32_t global_x = exp_manager->exp_m_7_->localXtoGlobalX(x);
          int32_t global_y = exp_manager->exp_m_7_->localYtoGlobalY(y);
          int64_t parent_id =
              global_x * exp_manager->exp_s()->global_grid_height() + global_y;
#else
        int64_t parent_id = repl->parent_id();
#endif
          int32_t count                = 0;
          MullerStruct* parent_mstruct = nullptr;
          for (auto& mstruct: muller_list) {

            if (std::find(mstruct->previous_generation_indexes.begin(),
                          mstruct->previous_generation_indexes.end(),
                          parent_id) !=
                mstruct->previous_generation_indexes.end()) {
              parent_strain_id = mstruct->strain_id;
              parent_mstruct   = mstruct;
            }
          }

          MullerStruct* mstruct = new MullerStruct(
              MullerStruct::next_strain_id++, parent_strain_id, t);
          mstruct->parent_mstruct         = parent_mstruct;
          mstruct->population_size_at[t0] = 1;
          parent_mstruct->children_strain_id.push_back(mstruct->strain_id);

#ifdef HAVE_MPI
          int32_t x = indices / exp_manager->grid_height();
          int32_t y = indices % exp_manager->grid_height();

          int32_t global_x = exp_manager->exp_m_7_->localXtoGlobalX(x);
          int32_t global_y = exp_manager->exp_m_7_->localYtoGlobalY(y);
          int64_t global_id =
              global_x * exp_manager->exp_s()->global_grid_height() + global_y;

          mstruct->current_generation_indexes.push_back(global_id);
#else
        mstruct->current_generation_indexes.push_back(index);
#endif

          muller_list.push_back(mstruct);
        } else {
#ifdef HAVE_MPI
          int32_t x = repl->parent_id() / exp_manager->grid_height();
          int32_t y = repl->parent_id() % exp_manager->grid_height();

          int32_t global_x = exp_manager->exp_m_7_->localXtoGlobalX(x);
          int32_t global_y = exp_manager->exp_m_7_->localYtoGlobalY(y);
          int64_t parent_id =
              global_x * exp_manager->exp_s()->global_grid_height() + global_y;
#else
        int64_t parent_id = repl->parent_id();
#endif

          for (auto& mstruct: muller_list) {
            if (std::find(mstruct->previous_generation_indexes.begin(),
                          mstruct->previous_generation_indexes.end(),
                          parent_id) !=
                mstruct->previous_generation_indexes.end()) {
#ifdef HAVE_MPI
              int32_t x = indices / exp_manager->grid_height();
              int32_t y = indices % exp_manager->grid_height();

              int32_t global_x = exp_manager->exp_m_7_->localXtoGlobalX(x);
              int32_t global_y = exp_manager->exp_m_7_->localYtoGlobalY(y);
              int64_t global_id =
                  global_x * exp_manager->exp_s()->global_grid_height() +
                  global_y;

              mstruct->current_generation_indexes.push_back(global_id);
#else
            mstruct->current_generation_indexes.push_back(index);
#endif

              break;
            }
          }
        }
      }

#ifdef HAVE_MPI
    }
#endif

    // CLEAN UP
    for (std::list<MullerStruct*>::iterator it = muller_list.begin();
         it != muller_list.end();) {
      if (((*it)->current_generation_indexes.size() == 0) &&
          ((*it)->children_strain_id.size() == 0)) {
        // delete within father
        if ((*it)->parent_mstruct != nullptr) {
          for (std::list<uint64_t>::iterator tit =
                   (*it)->parent_mstruct->children_strain_id.begin();
               tit != (*it)->parent_mstruct->children_strain_id.end();) {
            if ((*tit) == (*it)->strain_id) {
              tit = (*it)->parent_mstruct->children_strain_id.erase(tit);
            } else
              ++tit;
          }
        }

        (*it)->to_delete = true;
      }
      ++it;
    }

    for (std::list<MullerStruct*>::iterator it = muller_list.begin();
         it != muller_list.end();) {
      if ((*it)->to_delete) {
        MullerStruct* ptr = (*it);
        it                = muller_list.erase(it);
        delete ptr;
      } else {
        ++it;
      }
    }

    // SWAP LIST
    for (auto& mstruct: muller_list) {
      mstruct->since_nb_generation++;
      mstruct->population_size_at[t] =
          mstruct->current_generation_indexes.size();
      mstruct->previous_generation_indexes.swap(
          mstruct->current_generation_indexes);
      mstruct->current_generation_indexes.clear();
    }

    // Write to muller plot data file
    for (auto& mstruct: muller_list) {
#ifdef HAVE_MPI
      int64_t nb_indiv = exp_manager->exp_s()->global_grid_width() *
                         exp_manager->exp_s()->global_grid_height();
#else
      int64_t nb_indiv = exp_manager->nb_indivs();
#endif

      double pop_percent =
          (((double)mstruct->population_size_at[t]) / ((double)nb_indiv)) *
          100.0;

      if (((mstruct->since_nb_generation >= min_nb_generation) &&
           (pop_percent >= min_pop_percent)) ||
          mstruct->previously_output) {
        if (mstruct->population_size_at[t] > 0) {
          mullerfile << t << "," << mstruct->strain_id << "," << pop_percent
                     << "," << mstruct->population_size_at[t] << std::endl;
          if (verbose)
            printf(
                "%ld -- Strain %lu (exists since %lu generations, first gen "
                "%lu) : Pop Percent %lf -- Number of individual in Current / "
                "Previous generation %lu / %lu -- Cache Pop at %ld : %lu\n",
                t, mstruct->strain_id, mstruct->since_nb_generation,
                mstruct->first_generation, pop_percent,
                mstruct->current_generation_indexes.size(),
                mstruct->previous_generation_indexes.size(), t,
                mstruct->population_size_at[t]);

          if (!mstruct->previously_output) {
            // Go back up the tree
            MullerStruct* parent_mstruct = mstruct->parent_mstruct;

            for (uint64_t tmp_t = mstruct->first_generation; tmp_t < t;
                 tmp_t++) {
              double t_pop_percent =
                  (((double)mstruct->population_size_at[tmp_t]) /
                   ((double)nb_indiv)) *
                  100.0;
              mullerfile << tmp_t << "," << mstruct->strain_id << ","
                         << t_pop_percent << ","
                         << mstruct->population_size_at[tmp_t] << std::endl;
              if (verbose)
                printf("%ld -- REVERSING_TIME Strain %lu from %lu (exists "
                       "since %lu generations, first gen %lu) : Pop Percent "
                       "%lf -- Number of individual in Current / Previous "
                       "generation %lu / %lu -- Cache Pop at %ld : %lu\n",
                       t, mstruct->strain_id, mstruct->parent_strain_id,
                       mstruct->since_nb_generation, mstruct->first_generation,
                       pop_percent, mstruct->current_generation_indexes.size(),
                       mstruct->previous_generation_indexes.size(), tmp_t,
                       mstruct->population_size_at[tmp_t]);
            }

            MullerStruct* prev_mstruct = mstruct;

            while (parent_mstruct != nullptr) {
              if ((!parent_mstruct
                        ->previously_output)) {  // && (parent_mstruct->children_strain_id.size() == 1)) {
                // MERGE
                for (uint64_t tmp_t = parent_mstruct->first_generation;
                     tmp_t < prev_mstruct->first_generation; tmp_t++) {
                  double t_pop_percent =
                      (((double)parent_mstruct->population_size_at[tmp_t]) /
                       ((double)nb_indiv)) *
                      100.0;
                  mullerfile << tmp_t << "," << mstruct->strain_id << ","
                             << t_pop_percent << ","
                             << parent_mstruct->population_size_at[tmp_t]
                             << std::endl;
                  if (verbose)
                    printf("%ld -- PARENT Strain %lu from %lu (CHILD %lu from "
                           "%lu) (exists since %lu generations, first gen %lu) "
                           ": Pop Percent %lf -- Number of individual in "
                           "Current / Previous generation %lu / %lu -- Cache "
                           "Pop at %ld : %lu\n",
                           t, parent_mstruct->strain_id,
                           parent_mstruct->parent_strain_id, mstruct->strain_id,
                           mstruct->parent_strain_id,
                           parent_mstruct->since_nb_generation,
                           parent_mstruct->first_generation, pop_percent,
                           parent_mstruct->current_generation_indexes.size(),
                           parent_mstruct->previous_generation_indexes.size(),
                           tmp_t, parent_mstruct->population_size_at[tmp_t]);
                }

                prev_mstruct              = parent_mstruct;
                mstruct->parent_mstruct   = parent_mstruct->parent_mstruct;
                mstruct->parent_strain_id = parent_mstruct->parent_strain_id;
              } /*else if (!parent_mstruct->previously_output) {
              printf("INCOMPLETE IMPLEMENTATION\n");
              exit(-42);
            }*/
              else
                break;

              parent_mstruct = parent_mstruct->parent_mstruct;
            }

            mullerphenofile << mstruct->parent_strain_id << ","
                            << mstruct->strain_id << std::endl;
          }

          mstruct->previously_output = true;
        }
      }
    }
  }

#ifdef HAVE_MPI
  for (int32_t current_rank = 0; current_rank < exp_manager->nb_rank();
       current_rank++) {
    delete tree[current_rank];
  }
#else
  delete tree;  // delete the last Tree
#endif

  delete exp_manager;

  if (verbose)
    printf("OK\n");

  exit(EXIT_SUCCESS);
}

/**
 * \brief print help and exist
 */
void print_help(char* prog_path) {
  // Get the program file-name in prog_name (strip prog_path of the path)
  char* prog_name;  // No new, it will point to somewhere inside prog_path
  if ((prog_name = strrchr(prog_path, '/'))) {
    prog_name++;
  } else {
    prog_name = prog_path;
  }

  printf("*********************************************************************"
         "*********\n");
  printf("*                                                                    "
         "        *\n");
  printf("*                        aevol - Artificial Evolution                "
         "        *\n");
  printf("*                                                                    "
         "        *\n");
  printf("* Aevol is a simulation platform that allows one to let populations "
         "of       *\n");
  printf("* digital organisms evolve in different conditions and study "
         "experimentally  *\n");
  printf("* the mechanisms responsible for the structuration of the genome and "
         "the     *\n");
  printf("* transcriptome.                                                     "
         "        *\n");
  printf("*                                                                    "
         "        *\n");
  printf("*********************************************************************"
         "*********\n");
  printf("\n");
  printf("%s:\n", prog_name);
  printf(
      "\tReconstruct the lineage of a given individual from the tree files\n");
  printf("\n");
  printf("Usage : %s -h or --help\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("   or : %s [-b TIMESTEP] [-e TIMESTEP] [-m NB_GENERATIONS] [-p "
         "POPULATION_PERCENT] [-F] [-v]\n",
         prog_name);
  printf("\nOptions\n");
  printf("  -h, --help\n\tprint this help, then exit\n");
  printf("  -V, --version\n\tprint version number, then exit\n");
  printf("  -b, --begin TIMESTEP\n");
  printf("\tspecify time t0 up to which to reconstruct the lineage\n");
  printf("  -e, --end TIMESTEP\n");
  printf("\tspecify time t_end of the indiv whose lineage is to be "
         "reconstructed\n");
  printf("  -m, --min-nb-generation NB_GENERATIONS\n");
  printf("\tspecify the minimum number of generations a strain must last to "
         "appear in the output file\n");
  printf("  -p, --population-percent POPULATION_PERCENT\n");
  printf("\tspecify the minimum percent of the overall population a strain "
         "must cover\n");

  printf("  -v, --verbose\n\tbe verbose\n");
}

void interpret_cmd_line_options(int argc, char* argv[]) {
  // Define allowed options
  const char* short_options           = "hVb:e:m:p:v";
  static struct option long_options[] = {
      {"help", no_argument, nullptr, 'h'},
      {"version", no_argument, nullptr, 'V'},
      {"begin", required_argument, nullptr, 'b'},
      {"end", required_argument, nullptr, 'e'},
      {"min-nb-generation", required_argument, nullptr, 'm'},
      {"population-percent", required_argument, nullptr, 'p'},
      {"verbose", no_argument, nullptr, 'v'},
      {0, 0, 0, 0}};

  // Get actual values of the command-line options
  int option;
  while ((option = getopt_long(argc, argv, short_options, long_options,
                               nullptr)) != -1) {
    switch (option) {
    case 'h': {
      print_help(argv[0]);
      exit(EXIT_SUCCESS);
    }
    case 'V': {
      Utils::PrintAevolVersion();
      exit(EXIT_SUCCESS);
    }
    case 'b': {
      if (strcmp(optarg, "") == 0) {
        printf("%s: error: Option -b or --begin : missing argument.\n",
               argv[0]);
        exit(EXIT_FAILURE);
      }
      t0 = atol(optarg);
      break;
    }
    case 'e': {
      if (strcmp(optarg, "") == 0) {
        printf("%s: error: Option -e or --end : missing argument.\n", argv[0]);
        exit(EXIT_FAILURE);
      }
      t_end = atol(optarg);
      break;
    }
    case 'm': {
      if (strcmp(optarg, "") == 0) {
        printf(
            "%s: error: Option -m or --min-nb-generation : missing argument.\n",
            argv[0]);
        exit(EXIT_FAILURE);
      }
      min_nb_generation = atof(optarg);
      break;
    }
    case 'p': {
      if (strcmp(optarg, "") == 0) {
        printf("%s: error: Option -p or --population-percent : missing "
               "argument.\n",
               argv[0]);
        exit(EXIT_FAILURE);
      }
      min_pop_percent = atol(optarg);
      break;
    }
    case 'v': {
      verbose = true;
      break;
    }
    default: {
      // An error message is printed in getopt_long, we just need to exit
      exit(EXIT_FAILURE);
    }
    }
  }

  // If t_end wasn't provided, use default
  if (t_end < 0) {
    t_end = OutputManager::last_gener();
  }
}
