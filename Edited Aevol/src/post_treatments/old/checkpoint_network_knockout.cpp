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
#include <errno.h>
#include <inttypes.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <sys/stat.h>

#include <list>

#include <cstdint>
#include <fstream>
#include <limits>
#include <string>
// =================================================================
//                            Project Files
// =================================================================
#include "aevol.h"

using namespace aevol;

enum check_type
{
  FULL_CHECK  = 0,
  LIGHT_CHECK = 1,
  NO_CHECK    = 2
};

// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);
void extract_network(int time, Individual_R* indiv, double* fabs_metaerror_loss, double* fabs_fitness_loss,
                     double* fabs_metaerror_loss_percent, double* fabs_fitness_loss_percent);
void dump_network(int time, Individual_R* indiv, double* fabs_metaerror_loss, double* fabs_fitness_loss,
                  double* fabs_metaerror_loss_percent, double* fabs_fitness_loss_percent);
void filter_network(int time, Individual_R* indiv, double* fabs_metaerror_loss, double* fabs_fitness_loss,
                    double* fabs_metaerror_loss_percent, double* fabs_fitness_loss_percent);
void extract_network_single_target_model(int time, Individual_R* best, int nb_phenotypic_target_models,
                                         double** ptm_fabs_metaerror_loss, double** ptm_fabs_fitness_loss,
                                         double** ptm_fabs_metaerror_loss_percent,
                                         double** ptm_fabs_fitness_loss_percent);

int main(int argc, char** argv)
{
  // The output file (lineage.ae or lineage.rae) contains the following information:
  //
  // - common data                                                (ae_common::write_to_backup)
  // - begin gener                                                (int32_t)
  // - end gener                                                  (int32_t)
  // - final individual index                                     (int32_t)
  // - initial genome size                                        (int32_t)
  // - initial ancestor (nb genetic units + sequences)            (Individual::write_to_backup)
  // - replication report of ancestor at generation begin_gener+1 (ae_replic_report::write_to_backup)
  // - replication report of ancestor at generation begin_gener+2 (ae_replic_report::write_to_backup)
  // - replication report of ancestor at generation begin_gener+3 (ae_replic_report::write_to_backup)
  // - ...
  // - replication report of ancestor at generation end_gener     (ae_replic_report::write_to_backup)


  printf("\n  WARNING : Parameters' change in the middle of a simulation is not managed.\n");


  // =====================
  //  Parse command line
  // =====================

  // Default values
  check_type  check_genome      = LIGHT_CHECK;
  bool verbose = false;
  int64_t t0 = 0;
  int64_t t_end = -1;
  int32_t final_indiv_index = -1;
  int32_t final_indiv_rank  = -1;
  char tree_file_name[50];

  const char * short_options = "hVvncb:i:r:e:";
  static struct option long_options[] = {
    {"help",      no_argument,       NULL,  'h'},
    {"version",   no_argument,       NULL,  'V'},
    {"verbose",   no_argument,       NULL,  'v'},
    {"nocheck",   no_argument,       NULL,  'n'},
    {"fullcheck", no_argument,       NULL,  'c'},
    {"begin",     required_argument, NULL,  'b'},
    {"index",     required_argument, NULL,  'i'},
    {"rank",      required_argument, NULL,  'r'},
    {"end",       required_argument,  NULL, 'e' },
    {0, 0, 0, 0}
  };

  int option;
  while((option = getopt_long(argc, argv, short_options, long_options, NULL)) != -1)
  {
    switch(option)
    {
      case 'h' :
      {
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      }
      case 'V' :
      {
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      }
      case 'v' : verbose = true;                    break;
      case 'n' : check_genome = NO_CHECK;           break;
      case 'c' : check_genome = FULL_CHECK;         break;
      case 'b' : t0  = atol(optarg);                break;
      case 'i' : final_indiv_index  = atol(optarg); break;
      case 'r' : final_indiv_rank  = atol(optarg);  break;
      case 'e' :
      {
        if (strcmp(optarg, "") == 0)
        {
          printf("%s: error: Option -e or --end : missing argument.\n", argv[0]);
          exit(EXIT_FAILURE);
        }

        t_end = atol(optarg);

        break;
      }
    }
  }

  // Set undefined command line parameters to default values
  if (t_end == -1) {
    // Set t_end to the content of the LAST_GENER file if it exists.
    // If it doesn't, print help and exit
    FILE* lg_file = fopen(LAST_GENER_FNAME, "r");
    if (lg_file != NULL) {
      if (fscanf(lg_file, "%" PRId64, &t_end) == EOF) {
        printf("ERROR: failed to read last generation from file %s\n",
               LAST_GENER_FNAME);
        exit(EXIT_FAILURE);
      }
      fclose(lg_file);
    }
    else {
      printf("%s: error: You must provide a generation number.\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  // Load the simulation
  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(t_end, true, false);

    World* world = exp_manager->world();
    int16_t grid_width  = world->width();
    int16_t grid_height = world->height();
  // Check that the tree was recorded
  if (not exp_manager->record_tree()) {
    Utils::ExitWithUsrMsg("The phylogenetic tree wasn't recorded during "
                              "evolution, could not reconstruct the lineage");
  }

  int64_t tree_step = exp_manager->tree_step();

  //delete exp_manager;


  // The tree
  Tree* tree = NULL;

  // Indices, ranks and replication reports of the individuals in the lineage
  int32_t* indices = new int32_t[t_end - t0 + 1];
  //~ int32_t *                 ranks   = new int32_t[end_gener - begin_gener + 1];
  ReplicationReport** reports = new ReplicationReport*[t_end - t0];
  // NB: we do not need the report of the ancestor at generation begin_gener
  // (it might be the generation 0, for which we have no reports)
  // reports[0] = how ancestor at generation begin_gener + 1 was created
  // reports[i] = how ancestor at generation begin_gener + i + 1 was created
  // reports[end_gener - begin_gener - 1] = how the final individual was created
  //
  //            -----------------------------------------------------------------------------------------
  //  reports  | gener_0 => gener_1 | gener_1 => gener_2 | ... | gener_n-1 => gener_n | //////////////// |
  //            -----------------------------------------------------------------------------------------
  //  indices  |  index at gener_0  |  index at gener_1  | ... |  index at gener_n-1  | index at gener_n |
  //            -----------------------------------------------------------------------------------------


    // ============================
    // Init files
    // ============================
    std::ofstream network;
    network.open("checkpoint_network_knockout.csv",std::ofstream::trunc);
    network<<"Generation,"<<"Enhancer_or_Inhibitor,"<<"Value,"<<"Metaerror_lost,"<<"Fitness_lost,Metaerror_lost_percent,Fitness_lost_percent"
           <<std::endl;
    network.flush();
    network.close();

    float filter_values[3] = {0.01, 0.001, 0.005};

    for (float filter_value : filter_values) {

        std::string str_filter_value = std::to_string(filter_value);
        std::string file_name = "checkpoint_network_filtered_" + str_filter_value + ".csv";
        network.open(file_name, std::ofstream::trunc);
        network << "Generation," << "Enhancer," << "Inhibitor," << "Both," << "Value"
                << "Metaerror_lost,Fitness_lost,Metaerror_lost_percent,Fitness_lost_percent" << std::endl;
        network.flush();
        network.close();

        file_name = "checkpoint_network_edges_" + str_filter_value + ".csv";

        network.open(file_name, std::ofstream::trunc);
        network << "Generation," << "nb_enhancing," << "nb_inhibitor," << "nb_both,nb_edges," << "filter_nb_enhancing,"
                << "filter_nb_inhibitor," << "filter_nb_both,filter_nb_edges" << std::endl;
        network.flush();
        network.close();
    }

    float filter_values_2[4] = {0.0, 0.01, 0.001, 0.005};

    for (float filter_value : filter_values_2) {

        std::string str_filter_value = std::to_string(filter_value);
        std::string file_name = "checkpoint_network_dump_" + str_filter_value + ".csv";

        network.open(file_name, std::ofstream::trunc);
        network << "Generation," << "Source," << "Destination," << "Enhancer_or_Inhibitor," <<
                "Value" << "Metaerror_lost,Fitness_lost,Metaerror_lost_percent,Fitness_lost_percent" << std::endl;
        network.flush();
        network.close();
    }

    network.open("checkpoint_network_knockout_single_env.csv",std::ofstream::trunc);
    network<<"Generation,"<<"Enhancer_or_Inhibitor,"<<"TargetModel,"<<"Value"<<"Metaerror_lost,Fitness_lost,Metaerror_lost_percent,Fitness_lost_percent"<<std::endl;
    network.flush();
    network.close();


  // =========================
  //  Load the last tree file
  // =========================

  if (verbose)
  {
    printf("\n\n");
    printf("====================================\n");
    printf(" Loading the last tree file ... ");
    fflush(stdout);
  }


  // Example for ae_common::rec_params->tree_step() == 100 :
  //
  // tree_000100.ae ==>  generations   1 to 100.
  // tree_000200.ae ==>  generations 101 to 200.
  // tree_000300.ae ==>  generations 201 to 300.
  // etc.
  //
  // Thus, the information for generation end_gener are located
  // in the file called (end_gener/ae_common::rec_params->tree_step() + 1) * ae_common::rec_params->tree_step(),
  // except if end_gener%ae_common::rec_params->tree_step()==0.

  #ifdef __REGUL
    sprintf(tree_file_name,"tree/tree_%06" PRId64 ".ae", t_end);
  #else
    sprintf(tree_file_name,"tree/tree_%06" PRId64 ".ae", t_end);
  #endif

  tree = new Tree(exp_manager, tree_file_name);

  if (verbose)
  {
    printf("OK\n");
    printf("====================================\n");
  }


  // ============================================================================
  //  Find the index of the final individual and retrieve its replication report
  // ============================================================================
  if (final_indiv_index != -1)
  {
    // The index was directly provided, get the replication report and update the indices and ranks tables
    reports[t_end - t0 - 1] =
        new ReplicationReport(*(tree->report_by_index(t_end,
                                                          final_indiv_index)));
    final_indiv_rank = reports[t_end - t0 - 1]->rank();

    indices[t_end - t0]  = final_indiv_index;
  }
  else
  {
    if (final_indiv_rank == -1)
    {
      // No index nor rank was given in the command line.
      // By default, we construct the lineage of the best individual, the rank of which
      // is simply the number of individuals in the population.
      final_indiv_rank = exp_manager->nb_indivs();
    }

    // Retrieve the replication report of the individual of interest (at t_end)
    reports[t_end - t0 - 1] = new ReplicationReport(*(tree->report_by_rank(t_end, final_indiv_rank)));
    final_indiv_index = reports[t_end - t0 - 1]->id();

    indices[t_end - t0]  = final_indiv_index;
    //~ ranks[end_gener - begin_gener]    = final_indiv_rank;
  }

  if (verbose) printf("The final individual has the index %" PRId32 " (rank %" PRId32 ")\n", final_indiv_index, final_indiv_rank);


  // =======================
  //  Open the output file
  // =======================
  char output_file_name[101];

  #ifdef __REGUL
    snprintf(output_file_name, 100,
        "lineage-b%06" PRId64 "-e%06" PRId64 "-i%" PRId32 "-r%" PRId32 ".ae",
        t0, t_end, final_indiv_index, final_indiv_rank);
  #else
    snprintf(output_file_name, 100,
        "lineage-b%06" PRId64 "-e%06" PRId64 "-i%" PRId32 "-r%" PRId32 ".ae",
        t0, t_end, final_indiv_index, final_indiv_rank);
  #endif

  gzFile lineage_file = gzopen(output_file_name, "w");
  if (lineage_file == NULL)
  {
    fprintf(stderr, "File %s could not be created, exiting.\n", output_file_name);
    fprintf(stderr, "Please check your permissions in this directory.\n");
    exit(EXIT_FAILURE);
  }




  // ===================================================
  //  Retrieve the replication reports of the ancestors
  // ===================================================

  if (verbose)
  {
    printf("\n\n\n");
    printf("======================================================================\n");
    printf(" Parsing tree files to retrieve the ancestors' replication reports... \n");
    printf("======================================================================\n");
  }


  // Retrieve the index of the first ancestor from the last replication report
  indices[t_end - t0 -1] = reports[t_end - t0 - 1]->parent_id();

  // For each generation (going backwards), retrieve the index of the parent and
  // the corresponding replication report
  for (int64_t i = t_end - t0 - 2 ; i >= 0 ; i--)
  {
    int64_t t = t0 + i + 1;

    // We want to fill reports[i], that is to say, how the ancestor
    // at generation begin_gener + i + 1  was created
    if (verbose)
      printf("Getting the replication report for the ancestor at generation %" PRId64 "\n", t);

    // If we've exhausted the current tree file, load the next one
    if (Utils::mod(t, tree_step) == 0)
    {
      // Change the tree file
      delete tree;

      #ifdef __REGUL
        sprintf(tree_file_name,"tree/tree_%06" PRId64 ".ae", t);
      #else
        sprintf(tree_file_name,"tree/tree_%06" PRId64 ".ae", t);
      #endif

      tree = new Tree(exp_manager, tree_file_name);
    }

    // Copy the replication report of the ancestor
      //printf("Looking for the report %d at %d\n",t,indices[i + 1]);
    //bool notfound = true;
      for (int16_t x = 0 ; x < grid_width ; x++)
          for (int16_t y = 0 ; y < grid_height ; y++) {
              ReplicationReport* rep =  new ReplicationReport(*(tree->report_by_index(t,
                                                                                      x * grid_height + y)));
              if (rep->id() == indices[i + 1]) {
                  reports[i] = rep;
                  //printf("FOUND !!\n");
                  //notfound=false;
                  break;
              } else
                  delete rep;
          }

          //if (notfound) printf("ERROR NOT FOUND\n");



    // Retreive the index and rank of the next ancestor from the report
    indices[i] = reports[i]->parent_id();
  }
  delete exp_manager;


  if (verbose)  printf("OK\n");


  // =============================================================================
  //  Get the initial genome from the backup file and write it in the output file
  // =============================================================================

  if (verbose)
  {
    printf("\n\n\n");
    printf("=============================================== \n");
    printf(" Getting the initial genome sequence... ");
    fflush(NULL);
  }

  printf("Computing Network Knockout for generation %d\n",t0);
  // Load the simulation
  exp_manager = new ExpManager();
  exp_manager->load(t0, true, false);

  // Copy the initial ancestor
  // NB : The list of individuals is sorted according to the index
  Individual* initial_ancestor;
  //= exp_manager->indiv_by_id(indices[0]);
    for (int16_t x = 0 ; x < grid_width ; x++)
        for (int16_t y = 0 ; y < grid_height ; y++) {
            if (exp_manager->world()->indiv_at(x,y)->id() == indices[0]) {
                initial_ancestor = exp_manager->world()->indiv_at(x,y);
            }
        }


    Individual_R* best = dynamic_cast<Individual_R*>(initial_ancestor);
    best->do_transcription_translation_folding();

    int nb_edges = 0;
    for (auto &rna: best->get_rna_list_coding()) {
        nb_edges+=((Rna_R *) rna)->nb_influences();
    }

    double* fabs_metaerror_loss = new double[nb_edges];
    double* fabs_fitness_loss = new double[nb_edges];
    double* fabs_metaerror_loss_percent = new double[nb_edges];
    double* fabs_fitness_loss_percent = new double[nb_edges];

    for (int i = 0; i < nb_edges; i++) {
        fabs_metaerror_loss[i] = 0;
        fabs_fitness_loss[i] = 0;
        fabs_metaerror_loss_percent[i] = 0;
        fabs_fitness_loss_percent[i] = 0;
    }


    int nb_iteration = 10;
    printf("Running %d evals for %d edges\n",nb_iteration,nb_edges);
    for (int i = 0; i < nb_iteration; i++) {
        printf("Iteration %d\n",i);
        exp_manager->world()->ApplyHabitatVariation();

        best->evaluated_ = false;
        best->Evaluate();

        double base_metaerror = best->dist_to_target_by_feature(METABOLISM);
        double base_fitness = best->fitness();

        int i_edges = 0;

        for (auto &rna: best->get_rna_list_coding()) {
            for (unsigned int i = 0; i < ((Rna_R *) rna)->nb_influences(); i++) {
                double enhance_backup = ((Rna_R *) rna)->_enhancing_coef_list[i];
                double operate_backup = ((Rna_R *) rna)->_operating_coef_list[i];
                ((Rna_R *) rna)->_enhancing_coef_list[i] = 0;
                ((Rna_R *) rna)->_operating_coef_list[i] = 0;

                best->evaluated_ = false;
                best->Evaluate();

                fabs_metaerror_loss[i_edges] += std::fabs(base_metaerror-best->dist_to_target_by_feature(METABOLISM));
                fabs_fitness_loss[i_edges] += std::fabs(base_fitness-best->fitness());

                fabs_metaerror_loss_percent[i_edges] += (std::fabs(base_metaerror-best->dist_to_target_by_feature(METABOLISM)))/best->dist_to_target_by_feature(METABOLISM);
                fabs_fitness_loss_percent[i_edges] += (std::fabs(base_fitness-best->fitness()))/best->fitness();

                ((Rna_R *) rna)->_enhancing_coef_list[i] = enhance_backup;
                ((Rna_R *) rna)->_operating_coef_list[i] = operate_backup;

                i_edges++;
            }
        }
    }

    for (int i = 0; i < nb_edges; i++) {
        fabs_metaerror_loss[i] /= nb_iteration;
        fabs_fitness_loss[i] /= nb_iteration;
        fabs_metaerror_loss_percent[i] /= nb_iteration;
        fabs_fitness_loss_percent[i] /= nb_iteration;
    }

    extract_network(t0,best,fabs_metaerror_loss,fabs_fitness_loss,
                    fabs_metaerror_loss_percent,fabs_fitness_loss_percent);
    filter_network(t0,best,fabs_metaerror_loss,fabs_fitness_loss,fabs_metaerror_loss_percent,fabs_fitness_loss_percent);
    dump_network(t0,best,fabs_metaerror_loss,fabs_fitness_loss,fabs_metaerror_loss_percent,fabs_fitness_loss_percent);

    int nb_phenotypic_target_models = dynamic_cast<PhenotypicTargetHandler_R*>(exp_manager->world()->
            phenotypic_target_handler())->phenotypic_target_models_.size();
    printf("Running with a single phenotypic target models : %d\n",nb_phenotypic_target_models);

    double** ptm_fabs_metaerror_loss = new double*[nb_phenotypic_target_models];
    double** ptm_fabs_fitness_loss = new double*[nb_phenotypic_target_models];
    double** ptm_fabs_metaerror_loss_percent = new double*[nb_phenotypic_target_models];
    double** ptm_fabs_fitness_loss_percent = new double*[nb_phenotypic_target_models];

    for (int i = 0; i < nb_phenotypic_target_models; i++) {
        ptm_fabs_metaerror_loss[i] = new double[nb_edges];
        ptm_fabs_fitness_loss[i] = new double[nb_edges];
        ptm_fabs_metaerror_loss_percent[i] = new double[nb_edges];
        ptm_fabs_fitness_loss_percent[i] = new double[nb_edges];
    }


    for (int target_id = 0; target_id < nb_phenotypic_target_models; target_id++) {

        for (int j = 0; j < nb_edges; j++) {
            ptm_fabs_metaerror_loss[target_id][j] = 0;
            ptm_fabs_fitness_loss[target_id][j] = 0;
            ptm_fabs_metaerror_loss_percent[target_id][j] = 0;
            ptm_fabs_fitness_loss_percent[target_id][j] = 0;
        }

        printf("Testing with phenotypic target model %d\n",target_id);
        dynamic_cast<PhenotypicTargetHandler_R*>(exp_manager->world()->phenotypic_target_handler())->set_single_env(target_id);

        best->evaluated_ = false;
        best->Evaluate();

        double base_metaerror = best->dist_to_target_by_feature(METABOLISM);
        double base_fitness = best->fitness();

        int i_edges = 0;

        for (auto &rna: best->get_rna_list_coding()) {
            for (unsigned int i = 0; i < ((Rna_R *) rna)->nb_influences(); i++) {
                double enhance_backup = ((Rna_R *) rna)->_enhancing_coef_list[i];
                double operate_backup = ((Rna_R *) rna)->_operating_coef_list[i];
                ((Rna_R *) rna)->_enhancing_coef_list[i] = 0;
                ((Rna_R *) rna)->_operating_coef_list[i] = 0;

                best->evaluated_ = false;
                best->Evaluate();

                //printf("Testing with phenotypic target model %d : %lf %lf\n",target_id,base_metaerror,best->dist_to_target_by_feature(METABOLISM));

                ptm_fabs_metaerror_loss[target_id][i_edges] += std::fabs(base_metaerror-best->dist_to_target_by_feature(METABOLISM));
                ptm_fabs_fitness_loss[target_id][i_edges] += std::fabs(base_fitness-best->fitness());

                ptm_fabs_metaerror_loss_percent[target_id][i_edges] += (std::fabs(base_metaerror-best->dist_to_target_by_feature(METABOLISM)))/best->dist_to_target_by_feature(METABOLISM);
                ptm_fabs_fitness_loss_percent[target_id][i_edges] += (std::fabs(base_fitness-best->fitness()))/best->fitness();

                ((Rna_R *) rna)->_enhancing_coef_list[i] = enhance_backup;
                ((Rna_R *) rna)->_operating_coef_list[i] = operate_backup;

                i_edges++;
            }
        }
    }

    extract_network_single_target_model(t0,best,nb_phenotypic_target_models,ptm_fabs_metaerror_loss,ptm_fabs_fitness_loss,ptm_fabs_metaerror_loss_percent,ptm_fabs_fitness_loss_percent);


    delete best;

    if (verbose)
  {
    printf("OK\n");
    printf("=============================================== \n");
  }


  // ===============================================================================
  //  Write the replication reports of the successive ancestors in the output file
  //  (and, optionally, check that the rebuilt genome is correct each time a backup
  //  is available)
  // ===============================================================================

  if (verbose)
  {
    printf("\n\n\n");
    printf("============================================================ \n");
    printf(" Write the replication reports in the output file... \n");
    printf("============================================================ \n");
  }

  std::list<GeneticUnit>::const_iterator unit;

  Individual* stored_indiv = nullptr;
  std::list<GeneticUnit>::const_iterator stored_gen_unit;

  ExpManager* exp_manager_backup = NULL;

  // NB: I must keep the genome encapsulated inside an Individual, because
  // replaying the mutations has side effects on the list of promoters,
  // which is stored in the individual
  bool check_genome_now = false;

  for (int64_t i = 0 ; i < t_end - t0 ; i++)
  {
    // Where are we in time...
    int64_t t = t0 + i + 1;

      if (Utils::mod(t, exp_manager->backup_step()) == 0) {

          printf("Computing Network Knockout for generation %d\n",t);

          // Load the simulation
          exp_manager_backup = new ExpManager();
          exp_manager_backup->load(t, true, false);

          // Copy the ancestor from the backup
          for (int16_t x = 0; x < grid_width; x++)
              for (int16_t y = 0; y < grid_height; y++) {
                  if (exp_manager_backup->world()->indiv_at(x, y)->id() == indices[i + 1]) {
                      stored_indiv = exp_manager_backup->world()->indiv_at(x, y);
                      break;
                  }
              }

          Individual_R *best = dynamic_cast<Individual_R *>(stored_indiv);
          best->do_transcription_translation_folding();

          nb_edges = 0;
          for (auto &rna: best->get_rna_list_coding()) {
              nb_edges += ((Rna_R *) rna)->nb_influences();
          }

          delete[] fabs_metaerror_loss;
          delete[] fabs_fitness_loss;
          delete[] fabs_metaerror_loss_percent;
          delete[] fabs_fitness_loss_percent;

          fabs_metaerror_loss = new double[nb_edges];
          fabs_fitness_loss = new double[nb_edges];
          fabs_metaerror_loss_percent = new double[nb_edges];
          fabs_fitness_loss_percent = new double[nb_edges];

          for (int i = 0; i < nb_edges; i++) {
              fabs_metaerror_loss[i] = 0;
              fabs_fitness_loss[i] = 0;
              fabs_metaerror_loss_percent[i] = 0;
              fabs_fitness_loss_percent[i] = 0;
          }

          printf("Running %d evals for %d edges\n", nb_iteration, nb_edges);
          for (int i = 0; i < nb_iteration; i++) {
              printf("Iteration %d\n", i);
              exp_manager->world()->ApplyHabitatVariation();

              best->evaluated_ = false;
              best->Evaluate();

              double base_metaerror = best->dist_to_target_by_feature(METABOLISM);
              double base_fitness = best->fitness();

              int i_edges = 0;

              for (auto &rna: best->get_rna_list_coding()) {
                  for (unsigned int i = 0; i < ((Rna_R *) rna)->nb_influences(); i++) {
                      double enhance_backup = ((Rna_R *) rna)->_enhancing_coef_list[i];
                      double operate_backup = ((Rna_R *) rna)->_operating_coef_list[i];
                      ((Rna_R *) rna)->_enhancing_coef_list[i] = 0;
                      ((Rna_R *) rna)->_operating_coef_list[i] = 0;

                      best->evaluated_ = false;
                      best->Evaluate();

                      fabs_metaerror_loss[i_edges] += std::fabs(
                              base_metaerror - best->dist_to_target_by_feature(METABOLISM));
                      fabs_fitness_loss[i_edges] += std::fabs(base_fitness - best->fitness());

                      fabs_metaerror_loss_percent[i_edges] +=
                              (std::fabs(base_metaerror - best->dist_to_target_by_feature(METABOLISM))) /
                              best->dist_to_target_by_feature(METABOLISM);
                      fabs_fitness_loss_percent[i_edges] +=
                              (std::fabs(base_fitness - best->fitness())) / best->fitness();

                      ((Rna_R *) rna)->_enhancing_coef_list[i] = enhance_backup;
                      ((Rna_R *) rna)->_operating_coef_list[i] = operate_backup;

                      i_edges++;
                  }
              }
          }

          for (int i = 0; i < nb_edges; i++) {
              fabs_metaerror_loss[i] /= nb_iteration;
              fabs_fitness_loss[i] /= nb_iteration;
              fabs_metaerror_loss_percent[i] /= nb_iteration;
              fabs_fitness_loss_percent[i] /= nb_iteration;
          }

          extract_network(t,best, fabs_metaerror_loss, fabs_fitness_loss,
                          fabs_metaerror_loss_percent, fabs_fitness_loss_percent);
          filter_network(t,best, fabs_metaerror_loss, fabs_fitness_loss, fabs_metaerror_loss_percent,
                         fabs_fitness_loss_percent);
          dump_network(t,best, fabs_metaerror_loss, fabs_fitness_loss, fabs_metaerror_loss_percent,
                       fabs_fitness_loss_percent);

          printf("Running with a single phenotypic target models : %d\n", nb_phenotypic_target_models);

          for (int i = 0; i < nb_phenotypic_target_models; i++) {
              delete[] ptm_fabs_metaerror_loss[i];
              delete[] ptm_fabs_fitness_loss[i];
              delete[] ptm_fabs_metaerror_loss_percent[i];
              delete[] ptm_fabs_fitness_loss_percent[i];
          }

          for (int i = 0; i < nb_phenotypic_target_models; i++) {
              ptm_fabs_metaerror_loss[i] = new double[nb_edges];
              ptm_fabs_fitness_loss[i] = new double[nb_edges];
              ptm_fabs_metaerror_loss_percent[i] = new double[nb_edges];
              ptm_fabs_fitness_loss_percent[i] = new double[nb_edges];
          }


          for (int target_id = 0; target_id < nb_phenotypic_target_models; target_id++) {

              for (int j = 0; j < nb_edges; j++) {
                  ptm_fabs_metaerror_loss[target_id][j] = 0;
                  ptm_fabs_fitness_loss[target_id][j] = 0;
                  ptm_fabs_metaerror_loss_percent[target_id][j] = 0;
                  ptm_fabs_fitness_loss_percent[target_id][j] = 0;
              }

              printf("Testing with phenotypic target model %d\n", target_id);
              dynamic_cast<PhenotypicTargetHandler_R *>(exp_manager->world()->phenotypic_target_handler())->set_single_env(
                      target_id);

              best->evaluated_ = false;
              best->Evaluate();

              double base_metaerror = best->dist_to_target_by_feature(METABOLISM);
              double base_fitness = best->fitness();

              int i_edges = 0;

              for (auto &rna: best->get_rna_list_coding()) {
                  for (unsigned int i = 0; i < ((Rna_R *) rna)->nb_influences(); i++) {
                      double enhance_backup = ((Rna_R *) rna)->_enhancing_coef_list[i];
                      double operate_backup = ((Rna_R *) rna)->_operating_coef_list[i];
                      ((Rna_R *) rna)->_enhancing_coef_list[i] = 0;
                      ((Rna_R *) rna)->_operating_coef_list[i] = 0;

                      best->evaluated_ = false;
                      best->Evaluate();

                      //printf("Testing with phenotypic target model %d : %lf %lf\n",target_id,base_metaerror,best->dist_to_target_by_feature(METABOLISM));

                      ptm_fabs_metaerror_loss[target_id][i_edges] += std::fabs(
                              base_metaerror - best->dist_to_target_by_feature(METABOLISM));
                      ptm_fabs_fitness_loss[target_id][i_edges] += std::fabs(base_fitness - best->fitness());

                      ptm_fabs_metaerror_loss_percent[target_id][i_edges] +=
                              (std::fabs(base_metaerror - best->dist_to_target_by_feature(METABOLISM))) /
                              best->dist_to_target_by_feature(METABOLISM);
                      ptm_fabs_fitness_loss_percent[target_id][i_edges] +=
                              (std::fabs(base_fitness - best->fitness())) / best->fitness();

                      ((Rna_R *) rna)->_enhancing_coef_list[i] = enhance_backup;
                      ((Rna_R *) rna)->_operating_coef_list[i] = operate_backup;

                      i_edges++;
                  }
              }
          }

          extract_network_single_target_model(t,best, nb_phenotypic_target_models, ptm_fabs_metaerror_loss,
                                              ptm_fabs_fitness_loss, ptm_fabs_metaerror_loss_percent,
                                              ptm_fabs_fitness_loss_percent);


          delete best;
          delete exp_manager_backup;
      }

  }


  gzclose(lineage_file);
  delete [] reports;
  delete exp_manager;

  exit(EXIT_SUCCESS);
}

/*!
  \brief

*/
void print_help(char* prog_path)
{
  // default values :
  // begin_gener = 0
  // indiv  = best individual at generation end_gener

  // there must be a genome backup file for begin_gener

  // not relevant if crossover

  printf("\n");
  printf("*********************** aevol - Artificial Evolution ******************* \n");
  printf("*                                                                      * \n");
  printf("*                      Lineage post-treatment program                  * \n");
  printf("*                                                                      * \n");
  printf("************************************************************************ \n");
  printf("\n\n");
  printf("This program is Free Software. No Warranty.\n");
  printf("Copyright (C) 2009  LIRIS.\n");
  printf("\n");
#ifdef __REGUL
  printf("Usage : rlineage -h\n");
  printf("or :    rlineage [-vn] [-i index | -r rank] [-b gener1] -e end_gener \n");
#else
  printf("Usage : lineage -h\n");
  printf("or :    lineage [-vn] [-i index | -r rank] [-b gener1] -e end_gener \n");
#endif
  printf("\n");
#ifdef __REGUL
  printf("This program retrieves the ancestral lineage of an individual and writes \n");
  printf("it in an output file called lineage.rae. Specifically, it retrieves the \n");
  printf("lineage of the individual of end_gener whose index is index, going \n");
  printf("back in time up to gener1. This program requires at least one population backup\n");
  printf("file (for the generation gener1), one environment backup file (for the generation gener1)\n");
  printf("and all tree files for generations gener1 to end_gener.\n");
#else
  printf("This program retrieves the ancestral lineage of an individual and writes \n");
  printf("it in an output file called lineage.ae. Specifically, it retrieves the \n");
  printf("lineage of the individual of end_gener whose index is index, going \n");
  printf("back in time up to gener1. This program requires at least one population backup\n");
  printf("file (for the generation gener1), one environment backup file (for the generation gener1)\n");
  printf("and all tree files for generations gener1 to end_gener.\n");
#endif
  printf("\n");
  printf("WARNING: This program should not be used for simulations run with lateral\n");
  printf("transfer. When an individual has more than one parent, the notion of lineage\n");
  printf("used here is not relevant.\n");
  printf("\n");
  printf("\t-h or --help    : Display this help.\n");
  printf("\n");
  printf("\t-v or --verbose : Be verbose, listing generations as they are \n");
  printf("\t                  treated.\n");
  printf("\n");
  printf("\t-n or --nocheck    : Disable genome sequence checking. Makes the \n");
  printf("\t                       program faster, but it is not recommended. \n");
  printf("\t                       It is better to let the program check that \n");
  printf("\t                       when we rebuild the genomes of the ancestors\n");
  printf("\t                       from the lineage file, we get the same sequences\n");
  printf("\t                       as those stored in the backup files.\n");
  printf("\n");
  printf("\t-c or --fullcheck  : Will perform the genome checks every <BACKUP_STEP>\n");
  printf("\t                       generations. Default behaviour is lighter as it\n");
  printf("\t                       only performs these checks at the ending generation.\n");
  printf("\n");
  printf("\t-i index or --index index : \n");
  printf("\t                  Retrieve the lineage of the individual whose\n");
  printf("\t                  index is index. The index must be comprised \n");
  printf("\t                  between 0 and N-1, with N the size of the \n");
  printf("\t                  population at the ending generation. If neither\n");
  printf("\t                  index nor rank are specified, the program computes \n");
  printf("\t                  the lineage of the best individual of the ending \n");
  printf("\t                  generation.\n");
  printf("\n");
  printf("\t-r rank or --rank rank : \n");
  printf("\t                  Retrieve the lineage of the individual whose\n");
  printf("\t                  rank is rank. The rank must be comprised \n");
  printf("\t                  between 1 and N, with N the size of the \n");
  printf("\t                  population at the endind generation. If neither\n");
  printf("\t                  index nor rank are specified, the program computes \n");
  printf("\t                  the lineage of the best individual of the ending \n");
  printf("\t                  generation.\n");
  printf("\n");
  printf("\t-b gener1 or --begin gener1 : \n");
  printf("\t                  Retrieve the lineage up to generation gener1.\n");
  printf("\t                  There must be a genome backup file for this\n");
  printf("\t                  generation. If not specified, the program \n");
  printf("\t                  retrieves the lineage up to generation 0.\n");
  printf("\n");
  printf("\t-e end_gener or --end end_gener : \n");
  printf("\t                  Retrieve the lineage of the individual of end_gener \n");
  printf("\t                  (default: that contained in file last_gener.txt, if any)\n");
  printf("\n");
}



void extract_network(int time, Individual_R* indiv, double* fabs_metaerror_loss, double* fabs_fitness_loss, double* fabs_metaerror_loss_percent, double* fabs_fitness_loss_percent) {
    std::ofstream network;
    network.open("checkpoint_network_knockout.csv",std::ofstream::app);

    int i_edges = 0;

    for (auto& rna: indiv->get_rna_list_coding()) {
        for (unsigned int i = 0; i < ((Rna_R*)rna)->nb_influences(); i++) {
            //std::cout<<"Influence "<<i<<" value is "<<((Rna_R*)rna)->_enhancing_coef_list[i]<<" "<<((Rna_R*)rna)->_operating_coef_list[i]<<std::endl;
            //compute the activity
            if (((Rna_R*)rna)->_enhancing_coef_list[i] > 0)
            {
                network<<time<<",1,"<<((Rna_R*)rna)->_enhancing_coef_list[i]<<","<<fabs_metaerror_loss[i_edges]<<","
                       <<fabs_fitness_loss[i_edges]<<","<<fabs_metaerror_loss_percent[i_edges]<<","<<fabs_fitness_loss_percent[i_edges]<<std::endl;
            }

            if (((Rna_R*)rna)->_operating_coef_list[i] > 0)
            {
                network<<time<<",0,"<<((Rna_R*)rna)->_operating_coef_list[i]<<","<<fabs_metaerror_loss[i_edges]<<","
                       <<fabs_fitness_loss[i_edges]<<","<<fabs_metaerror_loss_percent[i_edges]<<","<<fabs_fitness_loss_percent[i_edges]<<std::endl;
            }
            i_edges++;
        }
    }

    network.flush();
    network.close();
}

void filter_network(int time, Individual_R* indiv, double* fabs_metaerror_loss, double* fabs_fitness_loss, double* fabs_metaerror_loss_percent, double* fabs_fitness_loss_percent) {
    float filter_values[3] = {0.01, 0.001, 0.005};

    for (float filter_value : filter_values) {

        std::string str_filter_value = std::to_string(filter_value);
        std::string file_name = "checkpoint_network_filtered_" + str_filter_value + ".csv";
        std::ofstream network;
        network.open(file_name, std::ofstream::app);

        int i_edges = 0;

        int nb_edges_enhance = 0, nb_edges_operating = 0, nb_edges_both = 0, nb_edges = 0;
        int filter_nb_edges_enhance = 0, filter_nb_edges_operating = 0, filter_nb_edges_both = 0, filter_nb_edges = 0;

        for (auto &rna: indiv->get_rna_list_coding()) {
            for (unsigned int i = 0; i < ((Rna_R *) rna)->nb_influences(); i++) {
                int both = 0;
                if (fabs_metaerror_loss[i_edges] >= filter_value) {
                    if ((((Rna_R *) rna)->_enhancing_coef_list[i] > 0) &&
                        (((Rna_R *) rna)->_operating_coef_list[i] > 0)) {
                        network << time << ",1,1,1," << ((Rna_R *) rna)->_enhancing_coef_list[i] << ","
                                << fabs_metaerror_loss[i_edges] << "," << fabs_fitness_loss[i_edges]<<
                                ","<<fabs_metaerror_loss_percent[i_edges]<<","<<fabs_fitness_loss_percent[i_edges]<< std::endl;
                        filter_nb_edges_both++;
                    } else {
                        //std::cout<<"Influence "<<i<<" value is "<<((Rna_R*)rna)->_enhancing_coef_list[i]<<" "<<((Rna_R*)rna)->_operating_coef_list[i]<<std::endl;
                        //compute the activity
                        if (((Rna_R *) rna)->_enhancing_coef_list[i] > 0) {
                            network << time << ",1,0,0," << ((Rna_R *) rna)->_enhancing_coef_list[i] << ","
                                    << fabs_metaerror_loss[i_edges] << "," << fabs_fitness_loss[i_edges]<<
                                    ","<<fabs_metaerror_loss_percent[i_edges]<<","<<fabs_fitness_loss_percent[i_edges]<< std::endl;
                            filter_nb_edges_enhance++;
                        }

                        if (((Rna_R *) rna)->_operating_coef_list[i] > 0) {
                            network << time << ",0,1,0," << ((Rna_R *) rna)->_operating_coef_list[i] << ","
                                    << fabs_metaerror_loss[i_edges] << "," << fabs_fitness_loss[i_edges]<<
                                    ","<<fabs_metaerror_loss_percent[i_edges]<<","<<fabs_fitness_loss_percent[i_edges]<< std::endl;
                            filter_nb_edges_operating++;
                        }
                    }
                    filter_nb_edges++;
                }

                if ((((Rna_R *) rna)->_enhancing_coef_list[i] > 0) && (((Rna_R *) rna)->_operating_coef_list[i] > 0)) {
                    nb_edges_both++;
                } else {
                    //std::cout<<"Influence "<<i<<" value is "<<((Rna_R*)rna)->_enhancing_coef_list[i]<<" "<<((Rna_R*)rna)->_operating_coef_list[i]<<std::endl;
                    //compute the activity
                    if (((Rna_R *) rna)->_enhancing_coef_list[i] > 0) {
                        nb_edges_enhance++;
                    }

                    if (((Rna_R *) rna)->_operating_coef_list[i] > 0) {
                        nb_edges_operating++;
                    }
                }
                nb_edges++;

                i_edges++;
            }
        }

        network.flush();
        network.close();

        file_name = "checkpoint_network_edges_" + str_filter_value + ".csv";

        network.open(file_name, std::ofstream::app);
        network << time << "," << nb_edges_enhance << "," << nb_edges_operating << "," << nb_edges_both << ","
                << nb_edges << "," <<
                filter_nb_edges_enhance << "," << filter_nb_edges_operating << "," << filter_nb_edges_both << ","
                << filter_nb_edges << std::endl;
        network.close();

    }

}



void dump_network(int time, Individual_R* indiv, double* fabs_metaerror_loss, double* fabs_fitness_loss, double* fabs_metaerror_loss_percent, double* fabs_fitness_loss_percent) {

    float filter_values[4] = {0.0, 0.01, 0.001, 0.005};

    for (float filter_value : filter_values) {

        std::string str_filter_value = std::to_string(filter_value);
        std::string file_name = "checkpoint_network_dump_"+str_filter_value+".csv";
        std::ofstream network;
        network.open(file_name, std::ofstream::app);

        int i_edges = 0;

        int nb_edges_enhance = 0, nb_edges_operating = 0, nb_edges_both = 0, nb_edges = 0;
        int filter_nb_edges_enhance = 0, filter_nb_edges_operating = 0, filter_nb_edges_both = 0, filter_nb_edges = 0;

        for (auto &rna: indiv->get_rna_list_coding()) {
            for (unsigned int i = 0; i < ((Rna_R *) rna)->nb_influences(); i++) {
                for (auto& protein : rna->transcribed_proteins()) {
                    if (fabs_metaerror_loss[i_edges] >= filter_value) {
                        if (((Rna_R *) rna)->_enhancing_coef_list[i] > 0) {
                            network << time << "," << protein->shine_dal_pos() << ","
                                    << dynamic_cast<Rna_R*>(rna)->_protein_list[i]->shine_dal_pos()<<","
                                    << "1," << ((Rna_R *) rna)->_enhancing_coef_list[i] << ","
                                    << fabs_metaerror_loss[i_edges] << "," << fabs_fitness_loss[i_edges] << std::endl;
                        }

                        if (((Rna_R *) rna)->_operating_coef_list[i] > 0) {
                            network << time << "," << protein->shine_dal_pos() << ","
                                    << dynamic_cast<Rna_R*>(rna)->_protein_list[i]->shine_dal_pos()<<","
                                    << "0," << ((Rna_R *) rna)->_operating_coef_list[i] << ","
                                    << fabs_metaerror_loss[i_edges] << "," << fabs_fitness_loss[i_edges] << std::endl;
                        }
                    }
                }

                i_edges++;
            }
        }

        network.flush();
        network.close();

    }

}

void extract_network_single_target_model(int time, Individual_R* indiv, int nb_phenotypic_target_models,
                                         double** ptm_fabs_metaerror_loss, double** ptm_fabs_fitness_loss,
                                         double** ptm_fabs_metaerror_loss_percent,
                                         double** ptm_fabs_fitness_loss_percent) {
    std::ofstream network;
    network.open("checkpoint_network_knockout_single_env.csv",std::ofstream::trunc);


    for (int target_id = 0; target_id < nb_phenotypic_target_models; target_id++) {
        int i_edges = 0;
        for (auto &rna: indiv->get_rna_list_coding()) {
            for (unsigned int i = 0; i < ((Rna_R *) rna)->nb_influences(); i++) {
                //std::cout<<"Influence "<<i<<" value is "<<((Rna_R*)rna)->_enhancing_coef_list[i]<<" "<<((Rna_R*)rna)->_operating_coef_list[i]<<std::endl;
                //compute the activity
                if (((Rna_R *) rna)->_enhancing_coef_list[i] > 0) {
                    network << time << ",1," <<target_id<<","<< ((Rna_R *) rna)->_enhancing_coef_list[i] << ","
                            << ptm_fabs_metaerror_loss[target_id][i_edges]<<","
                            << ptm_fabs_fitness_loss[target_id][i_edges] <<","
                            << ptm_fabs_metaerror_loss_percent[target_id][i_edges]<<","
                            << ptm_fabs_fitness_loss_percent[target_id][i_edges] << std::endl;
                }

                if (((Rna_R *) rna)->_operating_coef_list[i] > 0) {
                    network << time << ",0," <<target_id<<","<< ((Rna_R *) rna)->_operating_coef_list[i] << ","
                            << ptm_fabs_metaerror_loss[target_id][i_edges]<<","
                            << ptm_fabs_fitness_loss[target_id][i_edges] << ","
                            << ptm_fabs_metaerror_loss_percent[target_id][i_edges] <<","
                            << ptm_fabs_fitness_loss_percent[target_id][i_edges] << std::endl;
                }
                i_edges++;
            }
        }
    }

    network.flush();
    network.close();

}
