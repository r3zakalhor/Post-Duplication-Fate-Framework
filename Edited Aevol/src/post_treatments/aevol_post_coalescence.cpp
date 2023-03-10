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
//                              Libraries
// =================================================================
#include <algorithm>
#include <errno.h>
#include <fstream>
#include <getopt.h>
#include <inttypes.h>
#include <iostream>
#include <list>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unordered_map>
#include <vector>
#include <zlib.h>
// =================================================================
//                            Project Files
// =================================================================
#include "aevol.h"
#include "ExpManager_7.h"
using namespace aevol;

// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);

int main(int argc, char** argv) {
  printf("\n  WARNING : Parameters' change in the middle of a simulation is "
         "not managed.\n");

  // =====================
  //  Parse command line
  // =====================

  // Default values
  bool verbose = false;

  char tree_file_name[50];

  const char* short_options           = "hVv::";
  static struct option long_options[] = {{"help", no_argument, NULL, 'h'},
                                         {"version", no_argument, NULL, 'V'},
                                         {"verbose", no_argument, NULL, 'v'},
                                         {0, 0, 0, 0}};

  int option;
  while ((option = getopt_long(argc, argv, short_options, long_options,
                               NULL)) != -1) {
    switch (option) {
    case 'h': {
      print_help(argv[0]);
      exit(EXIT_SUCCESS);
    }
    case 'V': {
      Utils::PrintAevolVersion();
      exit(EXIT_SUCCESS);
    }
    case 'v':
      verbose = true;
      break;
      //case 'n' : check_genome = NO_CHECK;           break;
      //case 'c' : check_genome = FULL_CHECK;         break;
      //case 'b' : t0  = atol(optarg);                break;
      //case 'i' : final_indiv_index  = atol(optarg); break;
      //case 'r' : final_indiv_rank  = atol(optarg);  break;
    }
  }

  // There should be only one remaining arg: the lineage file
  if (optind != argc - 1) {
    Utils::ExitWithUsrMsg("please specify a lineage file");
  }

  char* lineage_file_name = new char[strlen(argv[optind]) + 1];
  sprintf(lineage_file_name, "%s", argv[optind]);

  // =======================
  //  Open the lineage file
  // =======================
  gzFile lineage_file = gzopen(lineage_file_name, "r");
  if (lineage_file == Z_NULL) {
    fprintf(stderr, "ERROR : Could not read the lineage file %s\n",
            lineage_file_name);
    exit(EXIT_FAILURE);
  }

  int64_t t0                = 0;
  int64_t t_end             = 0;
  int32_t final_indiv_index = 0;
  int32_t final_indiv_rank  = 0;

  gzread(lineage_file, &t0, sizeof(t0));
  gzread(lineage_file, &t_end, sizeof(t_end));
  gzread(lineage_file, &final_indiv_index, sizeof(final_indiv_index));
  gzread(lineage_file, &final_indiv_rank, sizeof(final_indiv_rank));

  // Load the simulation
  #ifdef HAVE_MPI
  int32_t current_rank = 0;
  #endif

  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(t0, true, false
  #ifdef HAVE_MPI
  , current_rank
  #endif
  );

  #ifdef HAVE_MPI
  exp_manager->exp_m_7_                     = new ExpManager_7(exp_manager);

  int32_t global_grid_width = exp_manager->exp_s()->global_grid_width();
  int32_t global_grid_height = exp_manager->exp_s()->global_grid_height();

  int32_t grid_width = exp_manager->grid_width();
  int32_t grid_height = exp_manager->grid_height();

  current_rank = exp_manager->exp_m_7_->rankOf(final_indiv_index / global_grid_height,final_indiv_index % global_grid_height);
  exp_manager->exp_m_7_->setRank(current_rank);
  #endif



  // Check that the tree was recorded
  if (not exp_manager->record_tree()) {
    Utils::ExitWithUsrMsg("The phylogenetic tree wasn't recorded during "
                          "evolution, could not reconstruct the lineage");
  }

  int64_t tree_step = exp_manager->tree_step();

  //delete exp_manager;

  // The tree
  Tree* tree = NULL;

  // =========================
  //  Load the last tree file
  // =========================

  GridCell* grid_cell = new GridCell(lineage_file, exp_manager, nullptr);
  auto* indiv         = grid_cell->individual();  //the best individual
  #ifdef HAVE_MPI
  int32_t local_index = indiv->id();
  int32_t local_x = local_index / grid_height;
  int32_t local_y = local_index % grid_height;
  int32_t x = exp_manager->exp_m_7_->localXtoGlobalX(local_x);
  int32_t y = exp_manager->exp_m_7_->localYtoGlobalY(local_y);
  int32_t index       = x*global_grid_height+y;
  #else
  int32_t index       = indiv->id();
  #endif

  int nb_muts = 0;

  ReplicationReport* rep_f = nullptr;

  World* world          = exp_manager->world();
  #ifdef HAVE_MPI
  unsigned int global_pop_size = global_grid_height * global_grid_width;
  unsigned int pop_size = grid_height * grid_width;
  #else
  int16_t grid_width    = world->width();
  int16_t grid_height   = world->height();
  unsigned int pop_size = grid_height * grid_width;
  #endif

  std::vector<int> coalescence_time;
  coalescence_time.resize(t_end);

  aevol::AeTime::set_time(t0);
  std::ofstream coalescence_file;
  coalescence_file.open("coalescence.csv", std::ofstream::trunc);
  coalescence_file << "Generation,"
                   << "Coalescence" << std::endl;

  #ifdef HAVE_MPI
  map<int, std::map<int,Tree*>> map_tree;
  #else
  map<int, Tree*> map_tree;
  printf("t0: %ld tend: %ld tree step: %ld pop size:  %u \n", t0, t_end,
         tree_step, pop_size);
  #endif

  while (time() <
         t_end) {  //time() is the number of generations. While less than max
    if (verbose)
      printf("Computing Coalescence at generation %" PRId64
             " for the lineage (index %" PRId32
             ")...",  //we calculate for the best individual
             time(), index);

    if (time() !=
        t0) {  //if it isnt the first go, get the index of the best individual
#ifdef __REGUL
      rep_f = new ReplicationReport(lineage_file,
                                    dynamic_cast<Individual_R*>(indiv));
#else
      rep_f = new ReplicationReport(lineage_file, indiv);
#endif

      index = rep_f->id();  // who we are building...
                            //printf("Update index %d\n",index);

      // For each genetic unit, replay the replication (undergo all mutations)
      const auto& dnarep = rep_f->dna_replic_report();

      nb_muts = dnarep.rearrangements().size() + dnarep.mutations().size();
    }

    //  if (nb_muts >= 1) {  //if the best individual has mutations/rearrangements
    #ifdef HAVE_MPI
    printf("nb_muts of %d (local %d): %d \n", index, local_index, nb_muts);
    #else
    printf("nb_muts of %d : %d \n", index, nb_muts);
    #endif
    // Search the coalescence time for this individual

    std::vector<int> previous;
    previous.push_back(
              index
        );  //add the index of best indiv to the end of previous

    std::vector<int> current;

    bool coal_found    = false;
    int coal_time      = 1;
    int64_t local_time = time() + 1;  //move along a generation

    //for each tree....why is this happening?
    for (auto t: map_tree) {  //seg fault is happening at this for loop

      if (t.first % 100 == 0) {
        break;
      }

      if (t.first)
        if (t.first <=
            time()) {  //if the first point in the tree is <= current gen

          #ifdef HAVE_MPI
          for (auto ptr : t.second) {
              delete ptr.second;
          }
          t.second.clear(); //deletes the value at the first point
          #else
          delete t.second;
          #endif

          map_tree.erase(t.first);  //gets rid of the first point in the tree
        }
      if (t.first == 0) {
        break;
      }
    }

    #ifdef HAVE_MPI
    for (auto ptr :  map_tree[time()]) {
       delete ptr.second;
    }
    #else
    delete map_tree[time()];  //get rid of the current point
    #endif

    #ifdef HAVE_MPI
    int32_t owner_rank;// TODO: Update rank
    #endif

    //get the tree before the current point
    if ((map_tree.find(((int)((local_time - 1) / tree_step) + 1) * tree_step) ==
        map_tree.end())
        ) {
       #ifdef HAVE_MPI
       for (int32_t rank = 0; rank < exp_manager->nb_rank(); rank++) {
       sprintf(tree_file_name, "tree/tree_" TIMESTEP_FORMAT "_%d.ae",
              ((int)((local_time - 1) / tree_step) + 1) * tree_step,rank);
       #else
       sprintf(tree_file_name, "tree/tree_" TIMESTEP_FORMAT ".ae",
              ((int)((local_time - 1) / tree_step) + 1) * tree_step);
       #endif
      //if this tree is the last one make a new tree map for the next
      //set of trees
      #ifdef HAVE_MPI
      map_tree[((int)((local_time - 1) / tree_step) + 1) * tree_step][rank] =
      #else
      map_tree[((int)((local_time - 1) / tree_step) + 1) * tree_step] =
      #endif
          new Tree(exp_manager, tree_file_name);
      tree = 
      #ifdef HAVE_MPI
      map_tree[((int)((local_time - 1) / tree_step) + 1) * tree_step][rank];
      #else
      map_tree[((int)((local_time - 1) / tree_step) + 1) * tree_step];
      #endif
      printf("Tree file name: %s \n", tree_file_name);
      printf("\n Loading tree %ld\n",
             ((int)((local_time - 1) / tree_step) + 1) * tree_step);
       #ifdef HAVE_MPI
       }
       #endif
    } else
    #ifdef HAVE_MPI
      tree = map_tree[((int)((local_time - 1) / tree_step) + 1) * tree_step][0];
    #else
      tree = map_tree[((int)((local_time - 1) / tree_step) + 1) * tree_step];
    #endif

    while (!coal_found) {
      //if we are at the origin, break
      if (local_time >= t_end) {  //local time is always 5000....
        printf("local time > t_end %ld \n", local_time);
        break;
      }

      //if the next point is a tree step
      if (Utils::mod(local_time - 1, tree_step) == 0) {

        //if we are at the end, generate the new tree map
        //it is doing this to fast.....
        if (map_tree.find(((int)((local_time - 1) / tree_step) + 1) *
                          tree_step) == map_tree.end()) {
          #ifdef HAVE_MPI
          for (int32_t rank = 0; rank < exp_manager->nb_rank(); rank++) {
          sprintf(tree_file_name, "tree/tree_" TIMESTEP_FORMAT "_%d.ae",
              ((int)((local_time - 1) / tree_step) + 1) * tree_step,rank);
          #else
          sprintf(tree_file_name, "tree/tree_" TIMESTEP_FORMAT ".ae",
                  ((int)((local_time - 1) / tree_step) + 1) * tree_step);
          #endif
       #ifdef HAVE_MPI
       map_tree[((int)((local_time - 1) / tree_step) + 1) * tree_step][rank] =
       #else
       map_tree[((int)((local_time - 1) / tree_step) + 1) * tree_step] =
       #endif
              new Tree(exp_manager, tree_file_name);
          tree =
          #ifdef HAVE_MPI
              map_tree[((int)((local_time - 1) / tree_step) + 1) * tree_step][rank];
          #else
              map_tree[((int)((local_time - 1) / tree_step) + 1) * tree_step];
          #endif
          printf("Tree file name: %s \n", tree_file_name);
          printf("Loading tree %ld\n",
                 ((int)((local_time - 1) / tree_step) + 1) * tree_step);
       #ifdef HAVE_MPI
       }
       #endif
        } else {
          tree =
          #ifdef HAVE_MPI
              map_tree[((int)((local_time - 1) / tree_step) + 1) * tree_step][0];
          #else
          map_tree[((int)((local_time - 1) / tree_step) + 1) * tree_step];
          #endif
        }
      }



#ifdef HAVE_MPI
      #pragma omp parallel for
      for (int32_t rank = 0; rank < exp_manager->nb_rank(); rank++) {
      ReplicationReport** reports =  map_tree[((int)((local_time - 1) / tree_step) + 1) * tree_step][rank]->reports(local_time);
#else
      //get the rep reports for the tree at the current time
      ReplicationReport** reports = tree->reports(local_time);
      #pragma omp parallel for
#endif
      for (int i = 0; i < pop_size; i++) {
        ReplicationReport* rep = new ReplicationReport(*(reports[i]));

        //find the first time parent_id occurs between begin and end, returns end if not there
        auto foundPrevious =
            std::find(previous.begin(), previous.end(),
            rep->parent_id()
            );

        //if the parent id is in previous add the current id to current vector
        if (foundPrevious != previous.end()) {
#pragma omp critical
          {
                 current.push_back(
                     rep->id()
                     );
          }
        }

        delete rep;
      }
      #ifdef HAVE_MPI
      }
      #endif

      //if all the ids in the population have been added to current then we have coalescence
      if (current.size() ==
      #ifdef HAVE_MPI
      global_pop_size
      #else
      pop_size
      #endif
      ) {
        coalescence_time[time()] = coal_time;
        coal_found               = true;
      } else {
        local_time++;
        coal_time++;
        previous.swap(current);
        current.clear();
      }
    }

    delete rep_f;
    // } else {
    //  coalescence_time[time()] = coalescence_time[AeTime::time() - 1] + 1;
    //}
    printf("Generation: %ld, coal_time: %d", AeTime::time(),
           coalescence_time[AeTime::time()]);

    coalescence_file << AeTime::time() << ","
                     << coalescence_time[AeTime::time()] << std::endl;

    aevol::AeTime::plusplus();  //take a time step
    if (verbose)
      printf(" OK\n");
  }

  //  for (int gen = 0; gen < t_end; gen++) {
  //
  //  }

  coalescence_file.flush();
  coalescence_file.close();

  //delete exp_manager;

  free(lineage_file_name);

  exit(EXIT_SUCCESS);
}

/*!
  \brief

*/
void print_help(char* prog_path) {
  // default values :
  // begin_gener = 0
  // indiv  = best individual at generation end_gener

  // there must be a genome backup file for begin_gener

  // not relevant if crossover

  printf("\n");
  printf("*********************** aevol - Artificial Evolution "
         "******************* \n");
  printf("*                                                                    "
         "  * \n");
  printf("*                      Lineage post-treatment program                "
         "  * \n");
  printf("*                                                                    "
         "  * \n");
  printf("*********************************************************************"
         "*** \n");
  printf("\n\n");
  printf("This program is Free Software. No Warranty.\n");
  printf("Copyright (C) 2009  LIRIS.\n");
  printf("\n");
#ifdef __REGUL
  printf("Usage : rlineage -h\n");
  printf("or :    rlineage [-vn] [-i index | -r rank] [-b gener1] -e end_gener "
         "\n");
#else
  printf("Usage : lineage -h\n");
  printf(
      "or :    lineage [-vn] [-i index | -r rank] [-b gener1] -e end_gener \n");
#endif
  printf("\n");
#ifdef __REGUL
  printf("This program retrieves the ancestral lineage of an individual and "
         "writes \n");
  printf("it in an output file called lineage.rae. Specifically, it retrieves "
         "the \n");
  printf(
      "lineage of the individual of end_gener whose index is index, going \n");
  printf("back in time up to gener1. This program requires at least one "
         "population backup\n");
  printf("file (for the generation gener1), one environment backup file (for "
         "the generation gener1)\n");
  printf("and all tree files for generations gener1 to end_gener.\n");
#else
  printf("This program retrieves the ancestral lineage of an individual and "
         "writes \n");
  printf("it in an output file called lineage.ae. Specifically, it retrieves "
         "the \n");
  printf(
      "lineage of the individual of end_gener whose index is index, going \n");
  printf("back in time up to gener1. This program requires at least one "
         "population backup\n");
  printf("file (for the generation gener1), one environment backup file (for "
         "the generation gener1)\n");
  printf("and all tree files for generations gener1 to end_gener.\n");
#endif
  printf("\n");
  printf("WARNING: This program should not be used for simulations run with "
         "lateral\n");
  printf("transfer. When an individual has more than one parent, the notion of "
         "lineage\n");
  printf("used here is not relevant.\n");
  printf("\n");
  printf("\t-h or --help    : Display this help.\n");
  printf("\n");
  printf("\t-v or --verbose : Be verbose, listing generations as they are \n");
  printf("\t                  treated.\n");
  printf("\n");
  printf(
      "\t-n or --nocheck    : Disable genome sequence checking. Makes the \n");
  printf(
      "\t                       program faster, but it is not recommended. \n");
  printf(
      "\t                       It is better to let the program check that \n");
  printf("\t                       when we rebuild the genomes of the "
         "ancestors\n");
  printf("\t                       from the lineage file, we get the same "
         "sequences\n");
  printf("\t                       as those stored in the backup files.\n");
  printf("\n");
  printf("\t-c or --fullcheck  : Will perform the genome checks every "
         "<BACKUP_STEP>\n");
  printf("\t                       generations. Default behaviour is lighter "
         "as it\n");
  printf("\t                       only performs these checks at the ending "
         "generation.\n");
  printf("\n");
  printf("\t-i index or --index index : \n");
  printf("\t                  Retrieve the lineage of the individual whose\n");
  printf("\t                  index is index. The index must be comprised \n");
  printf("\t                  between 0 and N-1, with N the size of the \n");
  printf(
      "\t                  population at the ending generation. If neither\n");
  printf("\t                  index nor rank are specified, the program "
         "computes \n");
  printf("\t                  the lineage of the best individual of the ending "
         "\n");
  printf("\t                  generation.\n");
  printf("\n");
  printf("\t-r rank or --rank rank : \n");
  printf("\t                  Retrieve the lineage of the individual whose\n");
  printf("\t                  rank is rank. The rank must be comprised \n");
  printf("\t                  between 1 and N, with N the size of the \n");
  printf(
      "\t                  population at the endind generation. If neither\n");
  printf("\t                  index nor rank are specified, the program "
         "computes \n");
  printf("\t                  the lineage of the best individual of the ending "
         "\n");
  printf("\t                  generation.\n");
  printf("\n");
  printf("\t-b gener1 or --begin gener1 : \n");
  printf("\t                  Retrieve the lineage up to generation gener1.\n");
  printf("\t                  There must be a genome backup file for this\n");
  printf("\t                  generation. If not specified, the program \n");
  printf("\t                  retrieves the lineage up to generation 0.\n");
  printf("\n");
  printf("\t-e end_gener or --end end_gener : \n");
  printf("\t                  Retrieve the lineage of the individual of "
         "end_gener \n");
  printf("\t                  (default: that contained in file last_gener.txt, "
         "if any)\n");
  printf("\n");
}
