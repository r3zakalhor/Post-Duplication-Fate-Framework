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

using namespace aevol;

// Helper functions
void interpret_cmd_line_options(int argc, char* argv[]);
void print_help(char* prog_path);

// Command-line option variables
static char* lineage_file_name = nullptr;
static bool verbose = false;
static bool full_check = false;
static bool trace_mutations = false;

static const double FILTER_PERCENT_VALUE[3] = {0.01, 0.005, 0.001};

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
    exit(EXIT_FAILURE);
  }

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



  // =============================
  //  Open the experiment manager
  // =============================
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
  std::string prefix = "stats/ancestor_stats/ancestor_network.csv";


  // ==================================================
  //  Prepare the initial ancestor and write its stats
  // ==================================================
  GridCell* grid_cell = new GridCell(lineage_file, exp_manager, nullptr);
#ifdef __REGUL
  auto* indiv = dynamic_cast<Individual_R*>(grid_cell->individual());
#else
  auto* indiv = grid_cell->individual();
#endif

#ifdef __REGUL
    std::ofstream network;
    network.open(prefix, std::ofstream::trunc);
    network << "Generation," << "nb_enhancing," << "nb_inhibitor," << "nb_both,nb_edges," << "filter_nb_enhancing,"
            << "filter_nb_inhibitor," << "filter_nb_both,filter_nb_edges" << std::endl;
    network.flush();


    Individual_R *best = dynamic_cast<Individual_R *>(indiv);
    best->clear_everything_except_dna_and_promoters();
    //best->do_transcription_translation_folding();
    //printf("---------> stoch env at %d\n",t);
    exp_manager->world()->ApplyHabitatVariation();

    //printf("---------> Evaluate indiv at %d\n",t);
    //best->evaluated_ = false;
    best->Evaluate();
    best->compute_statistical_data();
    printf("Initial fitness     = %e\n", best->fitness());
    printf("Initial genome size = %" PRId32 "\n", best->total_genome_size());
    
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

    double base_metaerror = 0;
    double base_fitness = 0;
    int nb_iteration = 1;

    for (int i = 0; i < nb_iteration; i++) {
        exp_manager->world()->ApplyHabitatVariation();
        Individual_R* cloned = new Individual_R(*best);
        cloned->set_grid_cell(exp_manager->world()->grid(0,0));

        cloned->clear_everything_except_dna_and_promoters();
        cloned->Evaluate();

// #pragma omp atomic
        base_metaerror = cloned->dist_to_target_by_feature(METABOLISM);

// #pragma omp atomic
        base_fitness = cloned->fitness();

        int i_edges = 0;

        for (auto &rna: cloned->get_rna_list_coding()) {
          for (unsigned int i = 0; i < ((Rna_R *) rna)->nb_influences(); i++) {
                double enhance_backup = ((Rna_R *) rna)->_enhancing_coef_list[i];
                double operate_backup = ((Rna_R *) rna)->_operating_coef_list[i];
                ((Rna_R *) rna)->_enhancing_coef_list[i] = 0;
                ((Rna_R *) rna)->_operating_coef_list[i] = 0;

                cloned->clear_everything_except_dna_and_promoters();
                cloned->Evaluate();

                fabs_metaerror_loss[i_edges] += std::fabs(base_metaerror-cloned->dist_to_target_by_feature(METABOLISM));
                fabs_fitness_loss[i_edges] += std::fabs(base_fitness-cloned->fitness());

                fabs_metaerror_loss_percent[i_edges] += (std::fabs(base_metaerror-cloned->dist_to_target_by_feature(METABOLISM)))
                                                          /cloned->dist_to_target_by_feature(METABOLISM);
                fabs_fitness_loss_percent[i_edges] += (std::fabs(base_fitness-cloned->fitness()))
                                                          /cloned->fitness();

                ((Rna_R *) rna)->_enhancing_coef_list[i] = enhance_backup;
                ((Rna_R *) rna)->_operating_coef_list[i] = operate_backup;

                i_edges++;
          }
        }
        delete cloned;
    }

    for (int i = 0; i < nb_edges; i++) {
      fabs_metaerror_loss[i] /= nb_iteration;
      fabs_fitness_loss[i] /= nb_iteration;
      fabs_metaerror_loss_percent[i] /= nb_iteration;
      fabs_fitness_loss_percent[i] /= nb_iteration;
    }


    for (auto filter : FILTER_PERCENT_VALUE) {
          // Counting edges
      int i_edges = 0;

      nb_edges = 0;
      int nb_edges_enhance = 0, nb_edges_operating = 0, nb_edges_both = 0;
      int filter_nb_edges_enhance = 0, filter_nb_edges_operating = 0, filter_nb_edges_both = 0, filter_nb_edges = 0;

      for (auto &rna: indiv->get_rna_list_coding()) {
        for (unsigned int i = 0; i < ((Rna_R *) rna)->nb_influences(); i++) {
          int both = 0;
          if (fabs_fitness_loss_percent[i_edges] >= filter) {
            if ((((Rna_R *) rna)->_enhancing_coef_list[i] > 0) &&
                (((Rna_R *) rna)->_operating_coef_list[i] > 0)) {
              filter_nb_edges_both++;
            } else {
              if (((Rna_R *) rna)->_enhancing_coef_list[i] > 0) {
                filter_nb_edges_enhance++;
              }

              if (((Rna_R *) rna)->_operating_coef_list[i] > 0) {
                filter_nb_edges_operating++;
              }
            }
            filter_nb_edges++;
          }

          if ((((Rna_R *) rna)->_enhancing_coef_list[i] > 0) && (((Rna_R *) rna)->_operating_coef_list[i] > 0)) {
            nb_edges_both++;
          } else {
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

      network << aevol::AeTime::time() << "," << filter<<","<<nb_edges_enhance << "," << nb_edges_operating << "," << nb_edges_both << ","
            << nb_edges << "," <<
            filter_nb_edges_enhance << "," << filter_nb_edges_operating << "," << filter_nb_edges_both << ","
            << filter_nb_edges << std::endl;
    }


    delete [] fabs_metaerror_loss;
    delete [] fabs_fitness_loss;
    delete [] fabs_metaerror_loss_percent;
    delete [] fabs_fitness_loss_percent;
#else
  indiv->Evaluate();
  indiv->compute_statistical_data();
  indiv->compute_non_coding();

  mystats->write_statistics_of_this_indiv(t0,indiv, nullptr);
#endif


  if (verbose) {
    printf("Initial fitness     = %f\n", indiv->fitness());
    printf("Initial genome size = %" PRId32 "\n", indiv->total_genome_size());
  }

  // ==========================================================================
  //  Replay the mutations to get the successive ancestors and analyze them
  // ==========================================================================
  ReplicationReport* rep = nullptr;
  int32_t index;
  ExpManager* exp_manager_backup = nullptr;
  int32_t unitlen_before;
  double metabolic_error_before;
  double impact_on_metabolic_error;
  double fitness_before;
  double impact_on_fitness;
  char mut_descr_string[255];

  bool check_now = false;

  aevol::AeTime::plusplus();
  while (time() <= t_end)
  {
    rep = new ReplicationReport(lineage_file, indiv);
    index = rep->id(); // who we are building...

    // Check now?
    check_now = time() == t_end ||
        (full_check && Utils::mod(time(), backup_step) == 0);

    if (verbose)
        printf("Rebuilding ancestor at generation %" PRId64
            " (index %" PRId32 ")...", time(), index);

    indiv->Reevaluate();

    // 2) Replay replication (create current individual's child)
    GeneticUnit& gen_unit = indiv->genetic_unit_nonconst(0);
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
    for (const auto& mut: dnarep.HT())
      gen_unit.dna()->undergo_this_mutation(*mut);

    for (const auto& mut: dnarep.rearrangements()) {
      if (trace_mutations) {
        // Store initial values before the mutation
        metabolic_error_before = indiv->dist_to_target_by_feature(METABOLISM);
	fitness_before = indiv->fitness();
        unitlen_before = gen_unit.dna()->length();
      }

      // Apply mutation
      gen_unit.dna()->undergo_this_mutation(*mut);

      if (trace_mutations) {
        indiv->Reevaluate();

        // Compute the metabolic impact of the mutation
        impact_on_metabolic_error =
            indiv->dist_to_target_by_feature(METABOLISM) -
            metabolic_error_before;

	impact_on_fitness =
	  indiv->fitness()-fitness_before;

        mut->generic_description_string(mut_descr_string);
      }
    }

    for (const auto& mut: dnarep.mutations()) {
      // Apply mutation
      gen_unit.dna()->undergo_this_mutation(*mut);
    }

    if (check_now) {
      if (verbose) {
        printf("Checking the sequence of the unit...");
        fflush(NULL);
      }

      char * str1 = new char[gen_unit.dna()->length() + 1];
      memcpy(str1, gen_unit.dna()->data(), \
             gen_unit.dna()->length()*sizeof(char));
      str1[gen_unit.dna()->length()] = '\0';

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
        delete indiv;
        delete stored_indiv;
        delete exp_manager_backup;
        delete exp_manager;
        exit(EXIT_FAILURE);
      }

      delete [] str1;
      delete [] str2;
    }

    // 3) All the mutations have been replayed, we can now evaluate the new individual

#ifdef __REGUL
    Individual_R *best = dynamic_cast<Individual_R *>(indiv);
    best->clear_everything_except_dna_and_promoters();
    
    exp_manager->world()->ApplyHabitatVariation();
    best->Evaluate();
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

    double base_metaerror = 0;
    double base_fitness = 0;
    int nb_iteration = 1;

    for (int i = 0; i < nb_iteration; i++) {
        exp_manager->world()->ApplyHabitatVariation();
        Individual_R* cloned = new Individual_R(*best);
        cloned->set_grid_cell(exp_manager->world()->grid(0,0));

        cloned->clear_everything_except_dna_and_promoters();
        cloned->Evaluate();

        base_metaerror = cloned->dist_to_target_by_feature(METABOLISM);

        base_fitness = cloned->fitness();

        int i_edges = 0;

        #pragma omp parallel for
        for (int j = 0; j < cloned->_rna_list_coding.size(); j++) {
          #pragma omp parallel for
          for (unsigned int i = 0; i < cloned->_rna_list_coding[j]->nb_influences(); i++) {
                Individual_R* l_cloned = new Individual_R(*cloned);
                l_cloned->set_grid_cell(exp_manager->world()->grid(0,0));
                l_cloned->clear_everything_except_dna_and_promoters();
                l_cloned->Evaluate();

                double enhance_backup = l_cloned->_rna_list_coding[j]->_enhancing_coef_list[i];
                double operate_backup = l_cloned->_rna_list_coding[j]->_operating_coef_list[i];
                l_cloned->_rna_list_coding[j]->_enhancing_coef_list[i] = 0;
                l_cloned->_rna_list_coding[j]->_operating_coef_list[i] = 0;

                l_cloned->clear_everything_except_dna_and_promoters();
                l_cloned->Evaluate();

                #pragma omp atomic
                fabs_metaerror_loss[i_edges] += std::fabs(base_metaerror-l_cloned->dist_to_target_by_feature(METABOLISM));

                #pragma omp atomic
                fabs_fitness_loss[i_edges] += std::fabs(base_fitness-l_cloned->fitness());

                #pragma omp atomic
                fabs_metaerror_loss_percent[i_edges] += (std::fabs(base_metaerror-l_cloned->dist_to_target_by_feature(METABOLISM)))
                                                          /l_cloned->dist_to_target_by_feature(METABOLISM);
                
                #pragma omp atomic
                fabs_fitness_loss_percent[i_edges] += (std::fabs(base_fitness-l_cloned->fitness()))
                                                          /l_cloned->fitness();

                delete l_cloned;

                i_edges++;
          }
        }
        delete cloned;
    }

    for (int i = 0; i < nb_edges; i++) {
      fabs_metaerror_loss[i] /= nb_iteration;
      fabs_fitness_loss[i] /= nb_iteration;
      fabs_metaerror_loss_percent[i] /= nb_iteration;
      fabs_fitness_loss_percent[i] /= nb_iteration;
    }
    

    for (auto filter : FILTER_PERCENT_VALUE) {
          // Counting edges
      int i_edges = 0;

      nb_edges = 0;
      int nb_edges_enhance = 0, nb_edges_operating = 0, nb_edges_both = 0;
      int filter_nb_edges_enhance = 0, filter_nb_edges_operating = 0, filter_nb_edges_both = 0, filter_nb_edges = 0;

      // int i_rna = 0;
      for (auto &rna: indiv->get_rna_list_coding()) {
        for (unsigned int i = 0; i < ((Rna_R *) rna)->nb_influences(); i++) {
          int both = 0;

          if (fabs_fitness_loss_percent[i_edges] >= filter) {
            if ((((Rna_R *) rna)->_enhancing_coef_list[i] > 0) &&
                (((Rna_R *) rna)->_operating_coef_list[i] > 0)) {
              filter_nb_edges_both++;
            } else {
              if (((Rna_R *) rna)->_enhancing_coef_list[i] > 0) {
                filter_nb_edges_enhance++;
              }

              if (((Rna_R *) rna)->_operating_coef_list[i] > 0) {
                filter_nb_edges_operating++;
              }
            }
            filter_nb_edges++;
          }
          // printf("%ld -- Rna %d Influence %d : Loss %lf Filter %lf :: %d || %d + %d + %d = %d\n",AeTime::time(),
          //         i_rna,i,fabs_fitness_loss_percent[i_edges],filter,(fabs_fitness_loss_percent[i_edges] >= filter),
          //         filter_nb_edges_enhance,filter_nb_edges_operating,filter_nb_edges_both,filter_nb_edges);
          // i_rna++;

          if ((((Rna_R *) rna)->_enhancing_coef_list[i] > 0) && (((Rna_R *) rna)->_operating_coef_list[i] > 0)) {
            nb_edges_both++;
          } else {
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

      network << aevol::AeTime::time() << "," << filter <<","<<nb_edges_enhance << "," << nb_edges_operating << "," << nb_edges_both << ","
              << nb_edges << "," <<
              filter_nb_edges_enhance << "," << filter_nb_edges_operating << "," << filter_nb_edges_both << ","
              << filter_nb_edges << std::endl;
    }

    delete [] fabs_metaerror_loss;
    delete [] fabs_fitness_loss;
    delete [] fabs_metaerror_loss_percent;
    delete [] fabs_fitness_loss_percent;
#else
    indiv->Reevaluate();
    indiv->compute_statistical_data();
    indiv->compute_non_coding();

    mystats->write_statistics_of_this_indiv(time(),indiv, rep);
#endif
    // Additional outputs
#ifndef __REGUL
    write_environment_stats(time(), phenotypicTargetHandler, env_output_file);
    write_terminators_stats(time(), indiv, term_output_file);
    if(phenotypicTargetHandler->phenotypic_target().nb_segments() > 1) {
      write_zones_stats(time(), indiv, phenotypicTargetHandler,
                        zones_output_file);
    }
    write_operons_stats(time(), indiv, operons_output_file);
#endif

    if (verbose) printf(" OK\n");

    delete rep;

    if (check_now) {
      delete stored_indiv;
      delete exp_manager_backup;
    }

    aevol::AeTime::plusplus();
  }

  network.close();

  gzclose(lineage_file);

  // Additional outputs
  // fclose(env_output_file);
  // fclose(term_output_file);
  // if(phenotypicTargetHandler->phenotypic_target().nb_segments() > 1) {
  //   fclose(zones_output_file);
  // }
  // fclose(operons_output_file);

  delete exp_manager;
  // delete mystats;
  delete indiv;

  return EXIT_SUCCESS;
}


void interpret_cmd_line_options(int argc, char* argv[]) {
  // =====================
  //  Parse command line
  // =====================
  const char * short_options = "hVFvM";
  static struct option long_options[] = {
    {"help",                no_argument, NULL, 'h'},
    {"version",             no_argument, NULL, 'V'},
    {"full-check",          no_argument, NULL, 'F'},
    {"trace-mutations",     no_argument, NULL, 'M'},
    {"verbose",             no_argument, NULL, 'v'},
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
      default:
        // An error message is printed in getopt_long, we just need to exit
        exit(EXIT_FAILURE);
    }
  }

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
}
