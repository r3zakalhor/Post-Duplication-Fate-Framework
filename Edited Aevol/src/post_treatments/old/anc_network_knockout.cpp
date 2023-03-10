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
#include <inttypes.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <err.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
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
  ENV_CHECK   = 2,
  NO_CHECK    = 3
};



// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);

void extract_network(Individual_R* indiv, double* fabs_metaerror_loss, double* fabs_fitness_loss,
                     double* fabs_metaerror_loss_percent, double* fabs_fitness_loss_percent);
void dump_network(Individual_R* indiv, double* fabs_metaerror_loss, double* fabs_fitness_loss,
                  double* fabs_metaerror_loss_percent, double* fabs_fitness_loss_percent);
void filter_network(Individual_R* indiv, double* fabs_metaerror_loss, double* fabs_fitness_loss,
                    double* fabs_metaerror_loss_percent, double* fabs_fitness_loss_percent);
void extract_network_single_target_model(Individual_R* best, int nb_phenotypic_target_models,
                                         double** ptm_fabs_metaerror_loss, double** ptm_fabs_fitness_loss,
                                         double** ptm_fabs_metaerror_loss_percent,
                                         double** ptm_fabs_fitness_loss_percent);


int main(int argc, char** argv)
{
  // The input file (lineage.ae or lineage.rae) must contain the following information:
  //
  // - common data                                                (ae_common::write_to_backup)
  // - begin gener                                                (int64_t)
  // - end gener                                                  (int64_t)
  // - final individual index                                     (int32_t)
  // - initial genome size                                        (int32_t)
  // - initial ancestor (nb genetic units + sequences)            (Individual::write_to_backup)
  // - replication report of ancestor at time t0+1  (ae_replic_report::write_to_backup)
  // - replication report of ancestor at time t0+2  (ae_replic_report::write_to_backup)
  // - replication report of ancestor at time t0+3  (ae_replic_report::write_to_backup)
  // - ...
  // - replication report of ancestor at time t_end_ (ae_replic_report::write_to_backup)




  // =====================
  //  Parse command line
  // =====================

  // Default values
  char*       lineage_file_name   = NULL;
  bool        verbose             = false;
  check_type  check               = LIGHT_CHECK;
  double      tolerance           = 0;

  const char * short_options = "hVvncf:lt:";
  static struct option long_options[] =
  {
    {"help",        no_argument,       NULL, 'h'},
    {"version",     no_argument,       NULL, 'V' },
    {"verbose",     no_argument,       NULL, 'v'},
    {"nocheck",     no_argument,       NULL, 'n'},
    {"fullcheck",   no_argument,       NULL, 'c'},
    {"file",        required_argument, NULL, 'f'},
    {"tolerance",   required_argument, NULL, 't'},
    {0, 0, 0, 0}
  };

  int option;
  while ((option = getopt_long(argc, argv, short_options, long_options, NULL)) != -1)
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
      case 'n' : check = NO_CHECK;                  break;
      case 'c' : check = FULL_CHECK;                break;
      case 'f' :
      {
        if (strcmp(optarg, "") == 0)
        {
          fprintf(stderr, "ERROR : Option -f or --file : missing argument.\n");
          exit(EXIT_FAILURE);
        }
        lineage_file_name = new char[strlen(optarg) + 1];
        sprintf(lineage_file_name, "%s", optarg);
        break;
      }
      case 't' :
      {
        if (strcmp(optarg, "") == 0)
        {
          fprintf(stderr, "ERROR : Option -t or --tolerance : missing argument.\n");
          exit(EXIT_FAILURE);
        }
        check = ENV_CHECK;
        tolerance = atof(optarg);
        break;
      }
      default :
      {
        fprintf(stderr, "ERROR : Unknown option, check your syntax.\n");
        print_help(argv[0]);
        exit(EXIT_FAILURE);
      }
    }
  }



  if (lineage_file_name == NULL)
  {
    fprintf(stderr, "ERROR : Option -f or --file missing. \n");
    exit(EXIT_FAILURE);
  }


  printf("\n");
  printf("WARNING : Parameter change during simulation is not managed in general.\n");
  printf("          Only changes in environmental target done with aevol_modify are handled.\n");
  printf("\n");

  // =======================
  //  Open the lineage file
  // =======================
  gzFile lineage_file = gzopen(lineage_file_name, "r");
  if (lineage_file == Z_NULL)
  {
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

  if (verbose)
  {
    printf("\n\n");
    printf("===============================================================================\n");
    printf(" Statistics of the ancestors of indiv. %" PRId32
           " (rank %" PRId32 ") from time %" PRId64 " to %" PRId64 "\n",
           final_indiv_index, final_indiv_rank, t0, t_end);
    printf("================================================================================\n");
  }

  // ============================
  // Init files
  // ============================
  std::ofstream network;
  network.open("anc_network_knockout.csv",std::ofstream::trunc);
  network<<"Generation,"<<"Enhancer_or_Inhibitor,"<<"Value,"<<"Metaerror_lost,"<<"Fitness_lost,Metaerror_lost_percent,Fitness_lost_percent"
         <<std::endl;
  network.flush();
  network.close();

  float filter_values[3] = {0.01, 0.001, 0.005};

  for (float filter_value : filter_values) {

    std::string str_filter_value = std::to_string(filter_value);
    std::string file_name = "anc_network_filtered_" + str_filter_value + ".csv";
    network.open(file_name, std::ofstream::trunc);
    network << "Generation," << "Enhancer," << "Inhibitor," << "Both," << "Value"
            << "Metaerror_lost,Fitness_lost,Metaerror_lost_percent,Fitness_lost_percent" << std::endl;
    network.flush();
    network.close();

    file_name = "anc_network_edges_" + str_filter_value + ".csv";

    network.open(file_name, std::ofstream::trunc);
    network << "Generation," << "nb_enhancing," << "nb_inhibitor," << "nb_both,nb_edges," << "filter_nb_enhancing,"
            << "filter_nb_inhibitor," << "filter_nb_both,filter_nb_edges" << std::endl;
    network.flush();
    network.close();
  }

  float filter_values_2[4] = {0.0, 0.01, 0.001, 0.005};

  for (float filter_value : filter_values_2) {

    std::string str_filter_value = std::to_string(filter_value);
    std::string file_name = "anc_network_dump_" + str_filter_value + ".csv";

    network.open(file_name, std::ofstream::trunc);
    network << "Generation," << "Source," << "Destination," << "Enhancer_or_Inhibitor," <<
            "Value" << "Metaerror_lost,Fitness_lost,Metaerror_lost_percent,Fitness_lost_percent" << std::endl;
    network.flush();
    network.close();
  }

  network.open("anc_network_knockout_single_env.csv",std::ofstream::trunc);
  network<<"Generation,"<<"Enhancer_or_Inhibitor,"<<"TargetModel,"<<"Value"<<"Metaerror_lost,Fitness_lost,Metaerror_lost_percent,Fitness_lost_percent"<<std::endl;
  network.flush();
  network.close();

  // =============================
  //  Open the experience manager
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
  if (not (phenotypicTargetHandler->var_method() == NO_VAR))
    Utils::ExitWithUsrMsg("sorry, ancestor stats has not yet been implemented "
                              "for variable phenotypic targets");

  int64_t backup_step = exp_manager->backup_step();

  // ==================================================
  //  Prepare the initial ancestor and write its stats
  // ==================================================
  GridCell* grid_cell = new GridCell(lineage_file, exp_manager, nullptr);
  // Individual*indiv = Individual::CreateIndividual(exp_manager, lineage_file);
  auto* indiv = grid_cell->individual();

  Individual_R* best = dynamic_cast<Individual_R*>(indiv);
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


  int nb_iteration = 100;
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

  extract_network(best,fabs_metaerror_loss,fabs_fitness_loss,
                  fabs_metaerror_loss_percent,fabs_fitness_loss_percent);
  filter_network(best,fabs_metaerror_loss,fabs_fitness_loss,fabs_metaerror_loss_percent,fabs_fitness_loss_percent);
  dump_network(best,fabs_metaerror_loss,fabs_fitness_loss,fabs_metaerror_loss_percent,fabs_fitness_loss_percent);

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

  extract_network_single_target_model(best,nb_phenotypic_target_models,ptm_fabs_metaerror_loss,ptm_fabs_fitness_loss,ptm_fabs_metaerror_loss_percent,ptm_fabs_fitness_loss_percent);


  // ==========================================================================
  //  Replay the mutations to get the successive ancestors and analyze them
  // ==========================================================================
  ReplicationReport* rep = nullptr;

  int32_t index;

  ExpManager* exp_manager_backup = nullptr;
  Habitat *backup_habitat = nullptr;

  bool check_now = false;

  aevol::AeTime::plusplus();
  while (time() <= t_end)
  {
    rep = new ReplicationReport(lineage_file, indiv);
    index = rep->id(); // who we are building...

    // Check now?
    check_now = ((check == FULL_CHECK && Utils::mod(time(), backup_step) == 0) ||
                 (check == ENV_CHECK && Utils::mod(time(), backup_step) == 0) ||
                 (check == LIGHT_CHECK && time() == t_end));

    if (verbose)
        printf("Rebuilding ancestor at generation %" PRId64
            " (index %" PRId32 ")...", time(), index);

    indiv->Reevaluate();

    // TODO <david.parsons@inria.fr> Check for phenotypic variation has to be
    // done for all the grid cells, disable checking until coded

//    // Check, and possibly update, the environment according to the backup files
//    // (update necessary if the env. was modified by aevol_modify at some point)
//    if (Utils::mod(time(), backup_step) == 0)
//    {
//      char world_file_name[255];
//      sprintf(world_file_name, "./" WORLD_FNAME_FORMAT, time());
//      gzFile world_file = gzopen(world_file_name, "r");
//      backup_habitat = new Habitat(world_file, pth); // TODO vld: fix pth
//
//      if (! env->is_identical_to(*backup_env, tolerance))
//      {
//        printf("Warning: At time()=%" PRId64 ", the replayed environment is not the same\n", time());
//        printf("         as the one saved at time()=%" PRId64 "... \n", time());
//        printf("         with tolerance of %lg\n", tolerance);
//        printf("Replacing the replayed environment by the one stored in the backup.\n");
//        delete env;
//        h = new Habitat(*backup_habitat);
//      }
//      delete backup_habitat;
//    }


    // Warning: this portion of code won'time() work if the number of units changes
    // during the evolution

    // 2) Replay replication (create current individual's child)
    GeneticUnit& gen_unit = indiv->genetic_unit_nonconst(0);
    GeneticUnit* stored_gen_unit = nullptr;
    Individual* stored_indiv = nullptr;

    if (check_now)
    {
      exp_manager_backup = new ExpManager();
      exp_manager_backup->load(time(), true, false);
      stored_indiv = new Individual(
          *(Individual*) exp_manager_backup->indiv_by_id(index));
      stored_gen_unit = &(stored_indiv->genetic_unit_nonconst(0));
    }

    // For each genetic unit, replay the replication (undergo all mutations)
    // TODO <david.parsons@inria.fr> disabled for multiple GUs
    const auto& dnarep = rep->dna_replic_report();

    dnarep.iter_muts([&](const auto& mut) {
      gen_unit.dna()->undergo_this_mutation(*mut);
    });

    if (check_now)
    {
      if (verbose)
      {
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
    Individual_R* best = dynamic_cast<Individual_R*>(indiv);
    best->do_transcription_translation_folding();

    nb_edges = 0;
    for (auto &rna: best->get_rna_list_coding()) {
      nb_edges+=((Rna_R *) rna)->nb_influences();
    }

    delete [] fabs_metaerror_loss;
    delete [] fabs_fitness_loss;
    delete [] fabs_metaerror_loss_percent;
    delete [] fabs_fitness_loss_percent;

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

    extract_network(best,fabs_metaerror_loss,fabs_fitness_loss,
                    fabs_metaerror_loss_percent,fabs_fitness_loss_percent);
    filter_network(best,fabs_metaerror_loss,fabs_fitness_loss,fabs_metaerror_loss_percent,fabs_fitness_loss_percent);
    dump_network(best,fabs_metaerror_loss,fabs_fitness_loss,fabs_metaerror_loss_percent,fabs_fitness_loss_percent);

    printf("Running with a single phenotypic target models : %d\n",nb_phenotypic_target_models);

    for (int i = 0; i < nb_phenotypic_target_models; i++) {
      delete [] ptm_fabs_metaerror_loss[i];
      delete [] ptm_fabs_fitness_loss[i];
      delete [] ptm_fabs_metaerror_loss_percent[i];
      delete [] ptm_fabs_fitness_loss_percent[i];
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

    extract_network_single_target_model(best,nb_phenotypic_target_models,ptm_fabs_metaerror_loss,ptm_fabs_fitness_loss,ptm_fabs_metaerror_loss_percent,ptm_fabs_fitness_loss_percent);


    if (verbose) printf(" OK\n");

    delete rep;

    if (check_now)
    {
      delete stored_indiv;
      delete exp_manager_backup;
    }

    aevol::AeTime::plusplus();
  }

  gzclose(lineage_file);

  delete exp_manager;
  delete indiv;

  exit(EXIT_SUCCESS);
}



void extract_network(Individual_R* indiv, double* fabs_metaerror_loss, double* fabs_fitness_loss, double* fabs_metaerror_loss_percent, double* fabs_fitness_loss_percent) {
  std::ofstream network;
  network.open("network_knockout.csv",std::ofstream::app);

  int i_edges = 0;

  for (auto& rna: indiv->get_rna_list_coding()) {
    for (unsigned int i = 0; i < ((Rna_R*)rna)->nb_influences(); i++) {
      //std::cout<<"Influence "<<i<<" value is "<<((Rna_R*)rna)->_enhancing_coef_list[i]<<" "<<((Rna_R*)rna)->_operating_coef_list[i]<<std::endl;
      //compute the activity
      if (((Rna_R*)rna)->_enhancing_coef_list[i] > 0)
      {
        network<<aevol::AeTime::time()<<",1,"<<((Rna_R*)rna)->_enhancing_coef_list[i]<<","<<fabs_metaerror_loss[i_edges]<<","
               <<fabs_fitness_loss[i_edges]<<","<<fabs_metaerror_loss_percent[i_edges]<<","<<fabs_fitness_loss_percent[i_edges]<<std::endl;
      }

      if (((Rna_R*)rna)->_operating_coef_list[i] > 0)
      {
        network<<aevol::AeTime::time()<<",0,"<<((Rna_R*)rna)->_operating_coef_list[i]<<","<<fabs_metaerror_loss[i_edges]<<","
               <<fabs_fitness_loss[i_edges]<<","<<fabs_metaerror_loss_percent[i_edges]<<","<<fabs_fitness_loss_percent[i_edges]<<std::endl;
      }
      i_edges++;
    }
  }

  network.flush();
  network.close();
}

void filter_network(Individual_R* indiv, double* fabs_metaerror_loss, double* fabs_fitness_loss, double* fabs_metaerror_loss_percent, double* fabs_fitness_loss_percent) {
  float filter_values[3] = {0.01, 0.001, 0.005};

  for (float filter_value : filter_values) {

    std::string str_filter_value = std::to_string(filter_value);
    std::string file_name = "network_filtered_" + str_filter_value + ".csv";
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
            network << aevol::AeTime::time() << ",1,1,1," << ((Rna_R *) rna)->_enhancing_coef_list[i] << ","
                    << fabs_metaerror_loss[i_edges] << "," << fabs_fitness_loss[i_edges]<<
                    ","<<fabs_metaerror_loss_percent[i_edges]<<","<<fabs_fitness_loss_percent[i_edges]<< std::endl;
            filter_nb_edges_both++;
          } else {
            //std::cout<<"Influence "<<i<<" value is "<<((Rna_R*)rna)->_enhancing_coef_list[i]<<" "<<((Rna_R*)rna)->_operating_coef_list[i]<<std::endl;
            //compute the activity
            if (((Rna_R *) rna)->_enhancing_coef_list[i] > 0) {
              network << aevol::AeTime::time() << ",1,0,0," << ((Rna_R *) rna)->_enhancing_coef_list[i] << ","
                      << fabs_metaerror_loss[i_edges] << "," << fabs_fitness_loss[i_edges]<<
                      ","<<fabs_metaerror_loss_percent[i_edges]<<","<<fabs_fitness_loss_percent[i_edges]<< std::endl;
              filter_nb_edges_enhance++;
            }

            if (((Rna_R *) rna)->_operating_coef_list[i] > 0) {
              network << aevol::AeTime::time() << ",0,1,0," << ((Rna_R *) rna)->_operating_coef_list[i] << ","
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

    file_name = "network_edges_" + str_filter_value + ".csv";

    network.open(file_name, std::ofstream::app);
    network << aevol::AeTime::time() << "," << nb_edges_enhance << "," << nb_edges_operating << "," << nb_edges_both << ","
            << nb_edges << "," <<
            filter_nb_edges_enhance << "," << filter_nb_edges_operating << "," << filter_nb_edges_both << ","
            << filter_nb_edges << std::endl;
    network.close();

  }

}



void dump_network(Individual_R* indiv, double* fabs_metaerror_loss, double* fabs_fitness_loss, double* fabs_metaerror_loss_percent, double* fabs_fitness_loss_percent) {

  float filter_values[4] = {0.0, 0.01, 0.001, 0.005};

  for (float filter_value : filter_values) {

    std::string str_filter_value = std::to_string(filter_value);
    std::string file_name = "network_dump_"+str_filter_value+".csv";
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
              network << aevol::AeTime::time() << "," << protein->shine_dal_pos() << ","
                      << dynamic_cast<Rna_R*>(rna)->_protein_list[i]->shine_dal_pos()<<","
                      << "1," << ((Rna_R *) rna)->_enhancing_coef_list[i] << ","
                      << fabs_metaerror_loss[i_edges] << "," << fabs_fitness_loss[i_edges] << std::endl;
            }

            if (((Rna_R *) rna)->_operating_coef_list[i] > 0) {
              network << aevol::AeTime::time() << "," << protein->shine_dal_pos() << ","
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

void extract_network_single_target_model(Individual_R* indiv, int nb_phenotypic_target_models,
                                         double** ptm_fabs_metaerror_loss, double** ptm_fabs_fitness_loss,
                                         double** ptm_fabs_metaerror_loss_percent,
                                         double** ptm_fabs_fitness_loss_percent) {
  std::ofstream network;
  network.open("network_knockout_single_env.csv",std::ofstream::trunc);


  for (int target_id = 0; target_id < nb_phenotypic_target_models; target_id++) {
    int i_edges = 0;
    for (auto &rna: indiv->get_rna_list_coding()) {
      for (unsigned int i = 0; i < ((Rna_R *) rna)->nb_influences(); i++) {
        //std::cout<<"Influence "<<i<<" value is "<<((Rna_R*)rna)->_enhancing_coef_list[i]<<" "<<((Rna_R*)rna)->_operating_coef_list[i]<<std::endl;
        //compute the activity
        if (((Rna_R *) rna)->_enhancing_coef_list[i] > 0) {
          network << aevol::AeTime::time() << ",1," <<target_id<<","<< ((Rna_R *) rna)->_enhancing_coef_list[i] << ","
                  << ptm_fabs_metaerror_loss[target_id][i_edges]<<","
                  << ptm_fabs_fitness_loss[target_id][i_edges] <<","
                  << ptm_fabs_metaerror_loss_percent[target_id][i_edges]<<","
                  << ptm_fabs_fitness_loss_percent[target_id][i_edges] << std::endl;
        }

        if (((Rna_R *) rna)->_operating_coef_list[i] > 0) {
          network << aevol::AeTime::time() << ",0," <<target_id<<","<< ((Rna_R *) rna)->_operating_coef_list[i] << ","
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

/*!
  \brief

*/
void print_help(char* prog_path)
{
  printf("\n");
  printf("*********************** aevol - Artificial Evolution ******************* \n");
  printf("*                                                                      * \n");
  printf("*                      Ancstats post-treatment program                 * \n");
  printf("*                                                                      * \n");
  printf("************************************************************************ \n");
  printf("\n\n");
  printf("This program is Free Software. No Warranty.\n");
  printf("Copyright (C) 2009  LIRIS.\n");
  printf("\n");
#ifdef __REGUL
  printf("Usage : rancstats -h\n");
  printf("or :    rancstats [-vn] -f lineage_file \n");
#else
  printf("Usage : ancstats -h\n");
  printf("or :    ancstats [-vn] -f lineage_file \n");
#endif
  printf("\n");
  printf("This program compute some statistics for the individuals within lineage_file.\n");
  printf("\n");
  printf("\n");
  printf("\t-h or --help       : Display this help.\n");
  printf("\n");
  printf("\t-v or --verbose    : Be verbose, listing generations as they are \n");
  printf("\t                       treated.\n");
  printf("\n");
  printf("\t-n or --nocheck    : Disable genome sequence checking. Makes the \n");
  printf("\t                       program faster, but it is not recommended. \n");
  printf("\t                       It is better to let the program check that \n");
  printf("\t                       when we rebuild the genomes of the ancestors\n");
  printf("\t                       from the lineage file, we get the same sequences\n");
  printf("\t                       as those stored in the backup files.\n");
  printf("\n");
  printf("\t-c or --fullcheck  : Will perform the genome and environment checks every\n");
  printf("\t                       <BACKUP_STEP> generations. Default behaviour is\n");
  printf("\t                       lighter as it only perform sthese checks at the\n");
  printf("\t                       ending generation.\n");
  printf("\n");
  printf("\t-f lineage_file or --file lineage_file : \n");
  printf("\t                       Compute the statistics for the individuals within lineage_file.\n");
  printf("\t-t tolerance or --tolerance tolerance : \n");
  printf("\t                       Tolerance used to compare the replayed environment to environment in backup\n");
  printf("\n");
}
