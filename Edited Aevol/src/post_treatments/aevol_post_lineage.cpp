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
#include <cerrno>
#include <cinttypes>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>

#include <zlib.h>
#include <sys/stat.h>
#include <getopt.h>

#include <list>
#include <queue>

#include "aevol.h"
//#ifdef HAVE_MPI
#include "ExpManager_7.h"
//#endif

enum check_type
{
    FULL_CHECK  = 0,
    LIGHT_CHECK = 1,
    NO_CHECK    = 2
};

using namespace aevol;

// Helper functions
void print_help(char* prog_path);
void interpret_cmd_line_options(int argc, char* argv[]);

// Command-line option variables
static bool full_check = false;
static bool best = false;
static bool verbose = false;
static int64_t t0 = 0;
static int64_t t_end = -1;
static int32_t final_indiv_index = -1;
static int32_t final_indiv_rank  = -1;
static char tree_file_name[255]; // TODO(dpa) remove magic number

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

  printf("\n  WARNING : Parameter change in the middle of a simulation is not managed.\n");
  
  #ifdef HAVE_MPI
  int32_t current_rank = 0;
  #endif

  // Load the simulation
  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(t_end, true, false
  #ifdef HAVE_MPI
  , current_rank
  #endif
  );

  #ifdef HAVE_MPI
  exp_manager->exp_m_7_                     = new ExpManager_7(exp_manager);
  #endif

  // Check that the tree was recorded
  if (not exp_manager->record_tree()) {
    Utils::ExitWithUsrMsg("The phylogenetic tree wasn't recorded during "
                          "evolution, could not reconstruct the lineage");
  }

  int64_t tree_step = exp_manager->tree_step();

  //delete exp_manager;


  // The tree
  Tree* tree = nullptr;

  // Indices, ranks and replication reports of the individuals in the lineage
  int32_t* indices = new int32_t[t_end - t0 + 1];
  ReplicationReport** reports = new ReplicationReport*[t_end - t0];
  // NB: we do not need the report of the ancestor at time t0 since we have
  // retrieved the individual itself from the initial backup
  // (plus it might be the generation 0, for which we have no reports)
  // reports[0] = how ancestor at t0 + 1 was created
  // reports[i] = how ancestor at t0 + i + 1 was created
  // reports[t_end - t0 - 1] = how the final individual was created
  //
  //           ---------------------------------------------------------------
  //  reports |  t0 => t1   |  t1 => t2   |...| t_n-1 => t_n   | XXXXXXXXXXXX |
  //           ---------------------------------------------------------------
  //  indices | index at t0 | index at t1 |...| index at t_n-1 | index at t_n |
  //           ---------------------------------------------------------------



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
    sprintf(tree_file_name,"tree/tree_" TIMESTEP_FORMAT "_%d.ae", t_end, current_rank);
#else
    sprintf(tree_file_name,"tree/tree_" TIMESTEP_FORMAT ".ae", t_end);
#endif

  printf("Loading %lld\n",t_end);

  tree = new Tree(exp_manager, tree_file_name);

  if (verbose) {
    printf("OK\n");
    printf("====================================\n");
  }


  // ============================================================================
  //  Find the index of the final individual and retrieve its replication report
  // ============================================================================
  if (best) {
    reports[t_end - t0 - 1] =
        new ReplicationReport(*(tree->report_by_index(t_end,
                                                      exp_manager->best_indiv()->id())));
    // final_indiv_rank = reports[t_end - t0 - 1]->rank();

    indices[t_end - t0]  = exp_manager->best_indiv()->id();
  } else if (final_indiv_index != -1) {
    // The index was directly provided, get the replication report and update the indices and ranks tables
    reports[t_end - t0 - 1] =
        new ReplicationReport(*(tree->report_by_index(t_end,
                                                      final_indiv_index)));
    // final_indiv_rank = reports[t_end - t0 - 1]->rank();

    indices[t_end - t0]  = final_indiv_index;
  }
  else {
    if (final_indiv_rank == -1) {
      // No index nor rank was given in the command line.
      // By default, we construct the lineage of the best individual, the rank of which
      // is simply the number of individuals in the population.
      printf("Loading individual GID %llu PID %llu\n",
              tree->report_by_index(t_end, 3)->id_,
              tree->report_by_index(t_end, 3)->parent_id_);

        reports[t_end - t0 - 1] =
                new ReplicationReport(*(tree->report_by_index(t_end, 3)));
    } else {
        reports[t_end - t0 - 1] =
                new ReplicationReport(*(tree->report_by_rank(t_end, final_indiv_rank)));
    }

    // Retrieve the replication report of the individual of interest (at t_end)
    // final_indiv_rank = reports[t_end - t0 - 1]->rank();
    final_indiv_index = reports[t_end - t0 - 1]->id();

    #ifdef HAVE_MPI
    reports[t_end - t0 - 1]->rank_ = current_rank;
    #endif

    indices[t_end - t0]  = final_indiv_index;
    //~ ranks[end_gener - begin_gener]    = final_indiv_rank;
  }

  if (verbose) {
    printf("The final individual has index %" PRId32
           " (rank %" PRId32 ")\n", final_indiv_index, final_indiv_rank);
  }


  // =======================
  //  Open the output file
  // =======================
  char output_file_name[101];

    snprintf(output_file_name, 100,
        "lineage-b" TIMESTEP_FORMAT "-e" TIMESTEP_FORMAT "-i%" PRId32 "-r%" PRId32 ".ae",
        t0, t_end, final_indiv_index, final_indiv_rank);

  gzFile lineage_file = gzopen(output_file_name, "w");
  if (lineage_file == nullptr) {
    fprintf(stderr, "File %s could not be created.\n", output_file_name);
    fprintf(stderr, "Please check your permissions in this directory.\n");
    exit(EXIT_FAILURE);
  }




  // ===================================================
  //  Retrieve the replication reports of the ancestors
  // ===================================================

  if (verbose) {
    printf("\n\n\n");
    printf("======================================================================\n");
    printf(" Parsing tree files to retrieve the ancestors' replication reports... \n");
    printf("======================================================================\n");
  }


  // Retrieve the index of the first ancestor from the last replication report
  indices[t_end - t0 -1] = reports[t_end - t0 - 1]->parent_id();

#ifdef HAVE_MPI
  int32_t global_grid_width = exp_manager->exp_s()->global_grid_width();
  int32_t global_grid_height = exp_manager->exp_s()->global_grid_height();

  int32_t grid_width = exp_manager->grid_width();
  int32_t grid_height = exp_manager->grid_height();

  int32_t rank_x = exp_manager->exp_s()->rank_width();
  int32_t rank_y = exp_manager->exp_s()->rank_height();

  int32_t TREE_CACHE_SIZE = 3;
  bool loadAll = false;

  if (TREE_CACHE_SIZE > exp_manager->exp_s()->nb_rank()) {
    TREE_CACHE_SIZE = exp_manager->exp_s()->nb_rank();
    loadAll = true;
  }

  Tree** list_tree = new Tree*[TREE_CACHE_SIZE];
  int32_t* list_rank = new int32_t[TREE_CACHE_SIZE];

  std::queue<int32_t> queue_rank;

  for (int32_t i_rank = 0; i_rank < TREE_CACHE_SIZE; i_rank++) {
    list_tree[i_rank] = nullptr;
    list_rank[i_rank] = -1;
  }
#endif

  // For each generation (going backwards), retrieve the index of the parent and
  // the corresponding replication report
  for (int64_t i = t_end - t0 - 2 ; i >= 0 ; i--) {
    int64_t t = t0 + i + 1;

    // We want to fill reports[i], that is to say, how the ancestor
    // at generation begin_gener + i + 1  was created
    if (verbose)
      printf("Getting the replication report for the ancestor at generation %" PRId64 " : %d\n", t,indices[i+1]);
    
    #ifdef HAVE_MPI
    int32_t owner_rank = exp_manager->exp_m_7_->rankOf(indices[i + 1] / global_grid_height,indices[i + 1] % global_grid_height);
    bool is_at_border = owner_rank != exp_manager->rank();

    printf("Looking for %d (%d %d) Current Rank %d Next One %d\n",indices[i + 1],
    indices[i + 1] / global_grid_height,indices[i + 1] % global_grid_height,
    exp_manager->rank(),owner_rank);

    // if (is_at_border) {

    // }
    #endif

    // If we've exhausted the current tree file, load the next one
    #ifdef HAVE_MPI
    if ((Utils::mod(t, tree_step) == 0) || (is_at_border))
    #else
    if (Utils::mod(t, tree_step) == 0)
    #endif
    {
      // Change the tree file
      
      #ifndef HAVE_MPI
      delete tree;
      #endif

      #ifdef HAVE_MPI

      current_rank = owner_rank;

      int to_load = -1;
      if (Utils::mod(t, tree_step) != 0) {
        to_load = (t / tree_step) * tree_step + tree_step;
        sprintf(tree_file_name,"tree/tree_" TIMESTEP_FORMAT "_%d.ae", to_load,current_rank);
      } else {
        AeTime::set_time(AeTime::time()-tree_step);
        sprintf(tree_file_name,"tree/tree_" TIMESTEP_FORMAT "_%d.ae", t,current_rank);
      }
        
        printf("Loading tree file %s\n",tree_file_name);
      #else
        sprintf(tree_file_name,"tree/tree_" TIMESTEP_FORMAT ".ae", t);
      #endif


      if (Utils::mod(t, tree_step) == 0)
        AeTime::set_time(AeTime::time()-tree_step);


      // printf("Loading another ExpManager (Switch from rank %d to %d\n",current_rank,owner_rank);
      //current_rank = owner_rank;

  #ifdef HAVE_MPI
  //     delete exp_manager->exp_m_7_;
  //     delete exp_manager;
  //     exp_manager = new ExpManager();
  //     exp_manager->load(t_end, true, false
  //     #ifdef HAVE_MPI
  //     , current_rank
  //     #endif
  //     );

  // exp_manager->exp_m_7_                     = new ExpManager_7(exp_manager);
      exp_manager->exp_m_7_->setRank(current_rank);
  

      if (is_at_border) {
        if (loadAll) {
          if (list_tree[current_rank]!=nullptr)
            tree = list_tree[current_rank];
          else {
            list_tree[current_rank] = new Tree(exp_manager, tree_file_name);
            tree = list_tree[current_rank];
            list_rank[current_rank] = current_rank;
          }
        } else {
          int32_t cell_id = -1;
          bool alreadyLoaded = false;

          for (int32_t i_cell = 0; i_cell < TREE_CACHE_SIZE; i_cell++) {
            if (list_rank[i_cell] == current_rank) {
              cell_id = i_cell;
              alreadyLoaded = true;
              break;
            }
          }

          if (!alreadyLoaded) {
            if (queue_rank.size() < TREE_CACHE_SIZE) {
              for (int32_t i_cell = 0; i_cell < TREE_CACHE_SIZE; i_cell++) {
                if (list_rank[i_cell] == -1) {
                  cell_id = i_cell;
                  break;
                }
              }
            } else {
              cell_id = queue_rank.front();
              queue_rank.pop();
            } 

            delete list_tree[cell_id];

            list_tree[cell_id] = new Tree(exp_manager, tree_file_name);
            list_rank[cell_id] = current_rank;
            queue_rank.push(cell_id);
          }

          tree = list_tree[cell_id];
        }
      } 
      
      if (Utils::mod(t, tree_step) == 0) {
        for (int32_t i_cell = 0; i_cell < TREE_CACHE_SIZE; i_cell++) {
            delete list_tree[i_cell];
            list_tree[i_cell] = nullptr;
            list_rank[i_cell] = -1;
            if (!queue_rank.empty()) queue_rank.pop();
        }

        list_tree[current_rank] = new Tree(exp_manager, tree_file_name);
        tree = list_tree[current_rank];
        list_rank[current_rank] = current_rank;
      }
    #else
      tree = new Tree(exp_manager, tree_file_name);
    #endif
    }



#ifdef HAVE_MPI
    int32_t local_rank_x = current_rank / rank_y;
    int32_t local_rank_y = current_rank % rank_y;

    int32_t local_x = exp_manager->exp_m_7_->globalXtoLocalX(indices[i + 1] / global_grid_height);
    int32_t local_y = exp_manager->exp_m_7_->globalYtoLocalY(indices[i + 1] % global_grid_height);
    int32_t local_id_ = local_x * grid_height + local_y;
#endif


    
    // Copy the replication report of the ancestor
    reports[i] =
        new ReplicationReport(*(tree->report_by_index(t, 
        #ifdef HAVE_MPI
        local_id_
        #else
        indices[i + 1]
        #endif
        )));
    // printf("Report for %d is %d => %d -- I %d\n",reports[i]->parent_id(),reports[i]->id(),indices[i+1]);

#ifdef HAVE_MPI
    printf("%d -- Loading RR for %d -- %d %d (Local ID %d -- %d %d) -- i : PID %d GID %d -- i+1 : PID %d GID %d\n",t,
           indices[i+1],indices[i + 1] / global_grid_height,indices[i + 1] % global_grid_height,
           local_id_,local_x,local_y,
           reports[i]->parent_id(),reports[i]->id(),reports[i+1]->parent_id(),reports[i+1]->id());
    #endif

    #ifdef HAVE_MPI
    reports[i]->rank_ = current_rank;
    #endif

    // Retreive the index and rank of the next ancestor from the report
    indices[i] = reports[i]->parent_id();
  }

  delete tree; // delete the last Tree
  //delete exp_manager->exp_m_7_;
  delete exp_manager;


  if (verbose) printf("OK\n");


  // =============================================================================
  //  Get the initial genome from the backup file and write it in the output file
  // =============================================================================

  if (verbose) {
    printf("\n\n\n");
    printf("=============================================== \n");
    printf(" Getting the initial genome sequence... ");
    fflush(nullptr);
  }

  // Load the simulation
  exp_manager = new ExpManager();
  exp_manager->load(t0, true, false
  #ifdef HAVE_MPI
  , reports[0]->rank_
  #endif
  );

      #ifdef HAVE_MPI
      exp_manager->exp_m_7_                     = new ExpManager_7(exp_manager);
      #endif

  // Copy the initial ancestor
  // NB : The list of individuals is sorted according to the index
  #ifdef HAVE_MPI
    int32_t orig_x = exp_manager->exp_m_7_->globalXtoLocalX(indices[0] / global_grid_height);
    int32_t orig_y = exp_manager->exp_m_7_->globalYtoLocalY(indices[0] % global_grid_height);
    printf("Loading GID %d : Cell %d %d on Rank %d\n",indices[0],orig_x,orig_y,reports[0]->rank_);
  #else
    int32_t orig_x = indices[0] / exp_manager->world()->height();
    int32_t orig_y = indices[0] % exp_manager->world()->height();
  #endif
      // printf("Loading GID %d : Cell %d %d\n",indices[0],orig_x,orig_y);

  const Individual& initial_ancestor = *(exp_manager->world()->indiv_at(orig_x, orig_y));

  // Write file "header"
  // printf("Write from %d to %d\n",t0,t_end);
  gzwrite(lineage_file, &t0, sizeof(t0));
  gzwrite(lineage_file, &t_end, sizeof(t_end));
  gzwrite(lineage_file, &final_indiv_index, sizeof(final_indiv_index));
  gzwrite(lineage_file, &final_indiv_rank, sizeof(final_indiv_rank));

  initial_ancestor.grid_cell()->save(lineage_file);


  if (verbose) {
    printf("OK\n");
    printf("=============================================== \n");
  }


  // ===============================================================================
  //  Write the replication reports of the successive ancestors in the output file
  //  (and, optionally, check that the rebuilt genome is correct each time a backup
  //  is available)
  // ===============================================================================

  if (verbose) {
    printf("\n\n\n");
    printf("============================================================ \n");
    printf(" Write the replication reports in the output file... \n");
    printf("============================================================ \n");
  }

  std::list<GeneticUnit>::const_iterator unit;

  Individual* stored_indiv = nullptr;
  std::list<GeneticUnit>::const_iterator stored_gen_unit;

  ExpManager* exp_manager_backup = nullptr;

  // NB: I must keep the genome encapsulated inside an Individual, because
  // replaying the mutations has side effects on the list of promoters,
  // which is stored in the individual
  bool check_genome_now = false;

  std::ofstream of;
  of.open("lignée_old.txt");

  for (int64_t i = 0 ; i < t_end - t0 ; i++) {
    // Where are we in time...
    int64_t t = t0 + i + 1;

    // Do we need to check the genome now?
    check_genome_now = t == t_end ||
        (full_check && Utils::mod(t, exp_manager->backup_step()) == 0);

    // Write the replication report of the ancestor for current generation
    if (verbose) {
      printf("Writing the replication report for t= %" PRId64
             " (built from indiv %" PRId32 " %" PRId32" at t= %" PRId64 ") %ld => %ld\n",
             t, indices[i], indices[i+1], t-1, reports[i]->parent_id(), reports[i]->id());
    }
    of << t << " : " << reports[i]->parent_id() << "  " << reports[i]->id() << std::endl;
    reports[i]->write_to_tree_file(lineage_file);
    if (verbose) printf(" OK\n");

      char *str_next;

    if (check_genome_now) {
      #ifdef HAVE_MPI
      printf("Loading ExpManager %d Rank %d\n",t,reports[i]->rank_);
      #endif
      
      // Load the simulation
      exp_manager_backup = new ExpManager();
      exp_manager_backup->load(t, true, false  
      #ifdef HAVE_MPI
      , reports[i]->rank_
      #endif
      );
  
      #ifdef HAVE_MPI
      exp_manager_backup->exp_m_7_                     = new ExpManager_7(exp_manager_backup);
      #endif

      // Copy the ancestor from the backup
      #ifdef HAVE_MPI
      int32_t local_x = exp_manager_backup->exp_m_7_->globalXtoLocalX(indices[i + 1] / global_grid_height);
      int32_t local_y = exp_manager_backup->exp_m_7_->globalYtoLocalY(indices[i + 1] % global_grid_height);
      stored_indiv = exp_manager_backup->world()->indiv_at(local_x,local_y);
      printf("Loading Individual %d (%d %d) : Rank %d Local %d %d\n",indices[i + 1],
              indices[i + 1] / global_grid_height,
              indices[i + 1] % global_grid_height, reports[i]->rank_, local_x,local_y);
      for (int i = 0; i < 3; i++)
        for (int j=0; j < 3; j++) {
          printf("%d -- DNA %d\n",exp_manager_backup->exp_m_7_->localXtoGlobalX(i) * global_grid_height + exp_manager_backup->exp_m_7_->localYtoGlobalY(j),
          
                      exp_manager_backup->world()->indiv_at(i,j)->genetic_unit_seq_length(0));
        }
      #else
      int x = indices[i+1] / exp_manager_backup->world()->height();
      int y = indices[i+1] % exp_manager_backup->world()->height();
      stored_indiv = exp_manager_backup->world()->indiv_at(x,y);
      #endif



      stored_gen_unit = stored_indiv->genetic_unit_list().cbegin();
    }


    // Warning: this portion of code won't work if the number of units changes
    // during the evolution

    // Replay the mutations stored in the current replication report on the
    // current genome

    // printf("Mutation numbers %d\n",reports[i]->dna_replic_report().mutations().size()+reports[i]->dna_replic_report().rearrangements().size());
    unit = initial_ancestor.genetic_unit_list().cbegin();
    reports[i]->dna_replic_report().iter_muts([&](const auto& mut) {
      // printf("Mutation at %d : %d (DNASize %d)\n",i,mut->mut_type(),unit->seq_length());
      (unit->dna())->undergo_this_mutation(*mut);
    });

    if (check_genome_now) {
      if (verbose) {
        printf("%lld -- Checking the sequence of the unit...",t);
        fflush(stdout);
      }
      assert(stored_gen_unit != stored_indiv->genetic_unit_list().cend());

        char * str1 = new char[unit->dna()->length() + 1];
      memcpy(str1, unit->dna()->data(),
             unit->dna()->length() * sizeof(char));
      str1[unit->dna()->length()] = '\0';

      char * str2 = new char[stored_gen_unit->dna()->length() + 1];
      memcpy(str2, stored_gen_unit->dna()->data(),
             stored_gen_unit->dna()->length() * sizeof(char));
      str2[stored_gen_unit->dna()->length()] = '\0';

      if (strncmp(str1, str2, stored_gen_unit->dna()->length()) == 0) {
        if (verbose) printf(" OK\n");
      }
      else {
        if (verbose) printf(" ERROR !\n");
        fprintf(stderr, "Error: the rebuilt unit is not the same as \n");
        fprintf(stderr, "the one stored in backup file at %" PRId64 "\n", t);
        fprintf(stderr, "Rebuilt unit : %" PRId32 " bp\n %s\n",
                (int32_t)strlen(str1), str1);
        fprintf(stderr, "Stored unit  : %" PRId32 " bp\n %s\n",
                (int32_t)strlen(str2), str2);

        delete [] str1;
        delete [] str2;
        gzclose(lineage_file);
        of.close();
        //delete exp_manager_backup->exp_m_7_;
        //delete exp_manager->exp_m_7_;
        delete exp_manager_backup;
        delete exp_manager;
        for (int64_t i_rep = 0 ; i_rep < t_end - t0 ; i_rep++) {
          delete reports[i_rep];
        }
        delete [] reports;
        delete [] indices;
        fflush(stdout);
        exit(EXIT_FAILURE);
      }

      delete [] str1;
      delete [] str2;

      ++stored_gen_unit;
    }
    ++unit;

    assert(unit == initial_ancestor.genetic_unit_list().cend());
    if (check_genome_now) {
      assert(stored_gen_unit == stored_indiv->genetic_unit_list().cend());
      //delete exp_manager_backup->exp_m_7_;
      delete exp_manager_backup;
    }
  }


  gzclose(lineage_file);
  of.close();
  for (int64_t i_rep = 0 ; i_rep < t_end - t0 ; i_rep++) {
    delete reports[i_rep];
  }
  delete [] reports;
  delete [] indices;
  //delete exp_manager->exp_m_7_;
  delete exp_manager;

  exit(EXIT_SUCCESS);
}

/**
 * \brief print help and exist
 */
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
  printf("%s:\n", prog_name);
  printf("\tReconstruct the lineage of a given individual from the tree files\n");
  printf("\n");
  printf("Usage : %s -h or --help\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("   or : %s [-b TIMESTEP] [-e TIMESTEP] [-I INDEX | -R RANK] [-F] [-v]\n",
         prog_name);
  printf("\nOptions\n");
  printf("  -h, --help\n\tprint this help, then exit\n");
  printf("  -V, --version\n\tprint version number, then exit\n");
  printf("  -b, --begin TIMESTEP\n");
  printf("\tspecify time t0 up to which to reconstruct the lineage\n");
  printf("  -e, --end TIMESTEP\n");
  printf("\tspecify time t_end of the indiv whose lineage is to be reconstructed\n");
  printf("  -I, --index INDEX\n");
  printf("\tspecify the index of the indiv whose lineage is to be reconstructed\n");
  printf("  -R, --rank RANK\n");
  printf("\tspecify the rank of the indiv whose lineage is to be reconstructed\n");
  printf("  -F, --full-check\n");
  printf("\tperform genome checks whenever possible\n");
  printf("  -v, --verbose\n\tbe verbose\n");
}

void interpret_cmd_line_options(int argc, char* argv[]) {
  // Define allowed options
  const char * short_options = "hVb:e:FI:R:vB";
  static struct option long_options[] = {
      {"help",      no_argument,       nullptr, 'h'},
      {"version",   no_argument,       nullptr, 'V'},
      {"begin",     required_argument, nullptr, 'b'},
      {"end",       required_argument, nullptr, 'e'},
      {"fullcheck", no_argument,       nullptr, 'F'},
      {"index",     required_argument, nullptr, 'I'},
      {"rank",      required_argument, nullptr, 'R'},
      {"verbose",   no_argument,       nullptr, 'v'},
      {"best",      no_argument,       nullptr, 'B'},
      {0, 0, 0, 0}
  };

  // Get actual values of the command-line options
  int option;
  while((option = getopt_long(argc, argv, short_options,
                              long_options, nullptr)) != -1) {
    switch(option) {
      case 'h' : {
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      }
      case 'V' : {
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      }
      case 'b' : {
        if (strcmp(optarg, "") == 0) {
          printf("%s: error: Option -b or --begin : missing argument.\n",
                 argv[0]);
          exit(EXIT_FAILURE);
        }
        t0  = atol(optarg);
        break;
      }
      case 'e' : {
        if (strcmp(optarg, "") == 0) {
          printf("%s: error: Option -e or --end : missing argument.\n",
                 argv[0]);
          exit(EXIT_FAILURE);
        }
        t_end = atol(optarg);
        break;
      }
      case 'F' : {
        full_check = true;
        break;
      }
      case 'B' : {
        best = true;
        break;
      }
      case 'I' : {
        final_indiv_index  = atoi(optarg);
        break;
      }
      case 'R' : {
        final_indiv_rank  = atoi(optarg);
        break;
      }
      case 'v' : {
        verbose = true;
        break;
      }
      default : {
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
