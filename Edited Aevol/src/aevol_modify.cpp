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


const char* DEFAULT_PARAM_FILE_NAME = "param.in";

// ============================================================================
//                                   Includes
// ============================================================================
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>

#include <getopt.h>

#include <list>

#if __cplusplus == 201103L
  #include "make_unique.h"
#endif

#include "ParameterLine.h"

#ifdef __X11
  #include "ExpManager_X11.h"
#else
  #include "ExpManager.h"
  #include "ParameterLine.h"
  #include "JumpingMT.h"
#endif

#include "ExpManager_7.h"

using namespace aevol;

enum population_change_type {
  SUBPOPULATIONS_BASED_ON_NON_CODING_BASES = 3,
  REMOVE_NON_CODING_BASES_BEST_IND = 4,
  REMOVE_NON_CODING_BASES_POPULATION = 5,
  DOUBLE_NON_CODING_BASES_BEST_IND = 6,
  DOUBLE_NON_CODING_BASES_POPULATION = 7
};

// Helper functions
void print_help(char* prog_path);
void interpret_cmd_line_options(int argc, char* argv[]);
ParameterLine* get_line(FILE* param_file);
void format_line(ParameterLine* formatted_line, char* line,
                 bool* line_is_interpretable);
// void change_by_cloning_best(ae_population* pop, ae_exp_manager* exp_m);
// void change_based_on_non_coding_bases_of_best_individual(ae_population* pop, ae_exp_manager* exp_m, population_change_type type);
// void change_based_on_non_coding_bases_in_population(ae_population* pop, ae_exp_manager* exp_m, population_change_type type);

// Command-line option variables
char* param_file_name = nullptr;
bool verbose = false;
int64_t timestep = -1;


int main(int argc, char* argv[]) {
  interpret_cmd_line_options(argc, argv);

  // Initialize the experiment manager
  #ifndef __NO_X
    ExpManager* exp_manager = new ExpManager_X11();
  #else
    ExpManager* exp_manager = new ExpManager();
  #endif
  exp_manager->load(timestep, false, verbose);

  // Define shorthands
  //Environment* env = exp_manager->env();
  Selection* sel = exp_manager->sel();
  World* world = exp_manager->world();
  ExpSetup* exp_s = exp_manager->exp_s();


  // If relevant, load the tree information
  char tree_file_name[50];
  Tree* tree = nullptr;
  bool take_care_of_the_tree = exp_manager->record_tree() &&
                               time() > 0;

  if (take_care_of_the_tree) {
    // If a tree is available, assign the replication reports to the individuals
    #ifdef __REGUL
      sprintf(tree_file_name,"tree/tree_" TIMESTEP_FORMAT ".ae", timestep);
    #else
      sprintf(tree_file_name, "tree/tree_" TIMESTEP_FORMAT ".ae", timestep);
    #endif

    gzFile tree_file = gzopen(tree_file_name, "r");

    if (tree_file != Z_NULL) {
      tree = new Tree(exp_manager, tree_file_name);
    } else
      take_care_of_the_tree = false;
  }






  // Interpret and apply changes
  printf("Interpret and apply changes\n");
  FILE* param_file = fopen(param_file_name, "r");
  if (param_file == nullptr) {
    printf("%s:%d: error: could not open parameter file %s\n", __FILE__,
           __LINE__, param_file_name);
    exit(EXIT_FAILURE);
  }

  std::list<Gaussian> new_gaussians;
  bool phen_target_change = false;
  bool start_to_record_tree = false;
  bool set_tree_step = false;
  int32_t tree_step = 100;

  int32_t grid_width_ = -1;
  int32_t grid_height_ = -1;


  ParameterLine* line;
  int32_t cur_line = 0;
  while ((line = get_line(param_file)) != nullptr) {
    cur_line++;
    if (strcmp(line->words[0], "ENV_AXIS_FEATURES") == 0) {
      // TODO <dpa> adapt to new organization
      printf(
          "%s:%d: error: ENV_AXIS_FEATURES has to be adapted to the new organization.\n",
          __FILE__, __LINE__);
      exit(EXIT_FAILURE);
//      int16_t env_axis_nb_segments = line->nb_words / 2;
//      double* env_axis_segment_boundaries = new double [env_axis_nb_segments + 1];
//      env_axis_segment_boundaries[0] = X_MIN;
//      for (int16_t i = 1 ; i < env_axis_nb_segments ; i++)
//      {
//        env_axis_segment_boundaries[i] = atof(line->words[2*i]);
//      }
//      env_axis_segment_boundaries[env_axis_nb_segments] = X_MAX;
//
//      // Set segment features
//      PhenotypicFeature* env_axis_features = new PhenotypicFeature[env_axis_nb_segments];
//      for (int16_t i = 0 ; i < env_axis_nb_segments ; i++)
//      {
//        if (strcmp(line->words[2*i+1], "NEUTRAL") == 0)
//        {
//          env_axis_features[i] = NEUTRAL;
//        }
//        else if (strcmp(line->words[2*i+1], "METABOLISM") == 0)
//        {
//          env_axis_features[i] = METABOLISM;
//        }
//        else if (strcmp(line->words[2*i+1], "SECRETION") == 0)
//        {
//          exp_manager->exp_s()->set_with_secretion(true);
//          env_axis_features[i] = SECRETION;
//        }
//        else if (strcmp(line->words[2*i+1], "DONOR") == 0)
//        {
//          env_axis_features[i] = DONOR;
//        }
//        else if (strcmp(line->words[2*i+1], "RECIPIENT") == 0)
//        {
//          env_axis_features[i] = RECIPIENT;
//        }
//        else
//        {
//          printf("ERROR in param file \"%s\" on line %" PRId32 " : unknown axis feature \"%s\".\n",
//                  param_file_name, cur_line, line->words[2*i+1]);
//          exit(EXIT_FAILURE);
//        }
//      }
//      env->set_segmentation(env_axis_nb_segments,
//                             env_axis_segment_boundaries,
//                             env_axis_features);
//      env_hasbeenmodified = true;
//      delete env_axis_segment_boundaries;
//      delete env_axis_features;
    }
    else if (strcmp(line->words[0], "RECORD_TREE") == 0) {
      if (strncmp(line->words[1], "true", 4) == 0) {
        start_to_record_tree = true;
      }
      else if (strncmp(line->words[1], "false", 5) == 0) {
        printf("ERROR stop recording tree is not implemented yet.\n");
        exit(EXIT_FAILURE);
      }
      else {
        printf("ERROR in param file \"%s\" on line %" PRId32
        " : unknown tree recording option (use true/false).\n",
            param_file_name, cur_line);
        exit(EXIT_FAILURE);
      }
      if (exp_manager->output_m()->record_tree()) {
        printf(
            "ERROR modification of already existing tree not impemented yet\n");
        exit(EXIT_FAILURE);
      }
    }
    else if (strcmp(line->words[0], "TREE_STEP") == 0) {
      tree_step = atol(line->words[1]);
      set_tree_step = true;
    }
    else if (strcmp(line->words[0], "TREE_MODE") == 0) {
      printf("ERROR : Tree mode management has been removed.\n");
      exit(EXIT_FAILURE);
    }
    else if (strcmp(line->words[0], "DUMP_STEP") == 0) {
      int step = atoi(line->words[1]);
      if (step > 0) {
        exp_manager->output_m()->set_dump_step(step);
      }
    }
    else if (strcmp(line->words[0], "BACKUP_STEP") == 0) {
      exp_manager->output_m()->set_backup_step(atol(line->words[1]));
    }
    else if (strcmp(line->words[0], "BIG_BACKUP_STEP") == 0) {
      exp_manager->output_m()->set_big_backup_step(atol(line->words[1]));
    }
    else if (strcmp(line->words[0], "POPULATION_SIZE") == 0) {
      printf("ERROR in param file \"%s\" on line %" PRId32
      ": the change of population size is not implemented yet\n"
      " for spatially structured populations",
          param_file_name, cur_line);
      exit(EXIT_FAILURE);
    }
    else if (strcmp(line->words[0], "SELECTION_SCHEME") == 0) {
      if (strncmp(line->words[1], "lin", 3) == 0) {
        if (line->nb_words != 3) {
          printf("ERROR in param file \"%s\" on line %" PRId32
          " : selection pressure parameter is missing.\n",
              param_file_name, cur_line);
          exit(EXIT_FAILURE);
        }
        sel->set_selection_scheme(RANK_LINEAR);
        sel->set_selection_pressure(atof(line->words[2]));
      }
      else if (strncmp(line->words[1], "exp", 3) == 0) {
        if (line->nb_words != 3) {
          printf("ERROR in param file \"%s\" on line %" PRId32
          " : selection pressure parameter is missing.\n",
              param_file_name, cur_line);
          exit(EXIT_FAILURE);
        }
        sel->set_selection_scheme(RANK_EXPONENTIAL);
        sel->set_selection_pressure(atof(line->words[2]));
      }
      else if (strncmp(line->words[1], "fitness", 7) == 0) {
        if (line->nb_words != 3) {
          printf("ERROR in param file \"%s\" on line %" PRId32
          " : selection pressure parameter is missing.\n",
              param_file_name, cur_line);
          exit(EXIT_FAILURE);
        }
        sel->set_selection_scheme(FITNESS_PROPORTIONATE);
        sel->set_selection_pressure(atof(line->words[2]));
      }
      else if (strcmp(line->words[1], "fittest") == 0) {
        sel->set_selection_scheme(FITTEST);
      }
      else {
        printf("ERROR in param file \"%s\" on line %" PRId32
        " : unknown selection scheme \"%s\".\n",
            param_file_name, cur_line, line->words[1]);
        exit(EXIT_FAILURE);
      }
    }
    else if (strcmp(line->words[0], "SELECTION_PRESSURE") == 0) {
      printf(
          "WARNING: SELECTION_PRESSURE keyword is outdated, you should\n"
          "specify a value for selection pressure using SELECTION_SCHEME\n");
      sel->set_selection_pressure(atof(line->words[1]));
      printf("\tChange of selection pressure to %f\n", atof(line->words[1]));
    }
    else if (strcmp(line->words[0], "POINT_MUTATION_RATE") == 0) {
      double point_mutation_rate = atof(line->words[1]);

      auto param_mut = exp_s->mut_params();
      param_mut->set_point_mutation_rate(point_mutation_rate);
      exp_s->set_mutation_parameters(param_mut);

      for (auto& indiv: exp_manager->indivs())
        indiv->set_point_mutation_rate(point_mutation_rate);
      printf("\tChange of overall point mutation rate to %f\n",
             point_mutation_rate);
    }
    else if (strcmp(line->words[0], "SMALL_INSERTION_RATE") == 0) {
      double small_insertion_rate = atof(line->words[1]);

      auto param_mut = exp_s->mut_params();
      param_mut->set_small_insertion_rate(small_insertion_rate);
      exp_s->set_mutation_parameters(param_mut);

      for (auto& indiv: exp_manager->indivs())
        indiv->set_small_insertion_rate(small_insertion_rate);
      printf("\tChange of overall small insertion rate to %f\n",
             small_insertion_rate);
    }
    else if (strcmp(line->words[0], "SMALL_DELETION_RATE") == 0) {
      double small_deletion_rate = atof(line->words[1]);

      auto param_mut = exp_s->mut_params();
      param_mut->set_small_deletion_rate(small_deletion_rate);
      exp_s->set_mutation_parameters(param_mut);

      for (auto& indiv: exp_manager->indivs())
        indiv->set_small_deletion_rate(small_deletion_rate);
      printf("\tChange of overall small deletion rate to %f\n",
             small_deletion_rate);
    }
    else if (strcmp(line->words[0], "MAX_INDEL_SIZE") == 0) {
      int16_t max_indel_size = atol(line->words[1]);
      for (auto& indiv: exp_manager->indivs())
        indiv->set_max_indel_size(max_indel_size);
      printf("\tChange of overall maximum indel size to %d\n", max_indel_size);
    }
    else if (strcmp(line->words[0], "DUPLICATION_RATE") == 0) {
      double duplication_rate = atof(line->words[1]);

      auto param_mut = exp_s->mut_params();
      param_mut->set_duplication_rate(duplication_rate);
      exp_s->set_mutation_parameters(param_mut);

      for (auto& indiv: exp_manager->indivs())
        indiv->set_duplication_rate(duplication_rate);
      printf("\tChange of overall duplication rate to %f\n", duplication_rate);
    }
    else if (strcmp(line->words[0], "DELETION_RATE") == 0) {
      double deletion_rate = atof(line->words[1]);

      auto param_mut = exp_s->mut_params();
      param_mut->set_deletion_rate(deletion_rate);
      exp_s->set_mutation_parameters(param_mut);

      for (auto& indiv: exp_manager->indivs())
        indiv->set_deletion_rate(deletion_rate);
      printf("\tChange of overall deletion rate to %f\n", deletion_rate);
    }
    else if (strcmp(line->words[0], "TRANSLOCATION_RATE") == 0) {
      double translocation_rate = atof(line->words[1]);

      auto param_mut = exp_s->mut_params();
      param_mut->set_translocation_rate(translocation_rate);
      exp_s->set_mutation_parameters(param_mut);

      for (auto& indiv: exp_manager->indivs())
        indiv->set_translocation_rate(translocation_rate);
      printf("\tChange of overall translocation rate to %f\n",
             translocation_rate);
    }
    else if (strcmp(line->words[0], "INVERSION_RATE") == 0) {
      double inversion_rate = atof(line->words[1]);
      
      auto param_mut = exp_s->mut_params();
      param_mut->set_inversion_rate(inversion_rate);
      exp_s->set_mutation_parameters(param_mut);

      for (auto& indiv: exp_manager->indivs())
        indiv->set_inversion_rate(inversion_rate);
      printf("\tChange of overall inversion to %f\n", inversion_rate);
    }
    else if (strcmp(line->words[0], "TRANSFER_INS_RATE") == 0) {
      double transfer_ins_rate = atof(line->words[1]);
      for (auto& indiv: exp_manager->indivs())
        indiv->set_HT_ins_rate(transfer_ins_rate);
      exp_manager->set_HT_ins_rate(transfer_ins_rate);
      printf("\tChange of overall transfer insertion rate to %f\n",
             transfer_ins_rate);
    }
    else if (strcmp(line->words[0], "TRANSFER_REPL_RATE") == 0) {
      double transfer_repl_rate = atof(line->words[1]);
      for (auto& indiv: exp_manager->indivs())
        indiv->set_HT_repl_rate(transfer_repl_rate);
      exp_manager->set_HT_repl_rate(transfer_repl_rate);
      printf("\tChange of overall transfer replacement rate to %f\n",
             transfer_repl_rate);
    }
    else if ((strcmp(line->words[0], "ENV_ADD_GAUSSIAN") == 0) ||
             (strcmp(line->words[0], "ENV_GAUSSIAN") == 0)) {
      // TODO <dpa> adapt to new organization
//      printf("%s:%d: error: ENV_ADD_GAUSSIAN has to be adapted to the new organization.\n", __FILE__, __LINE__);
//      exit(EXIT_FAILURE);

      new_gaussians.emplace_back(atof(line->words[1]),
                                 atof(line->words[2]),
                                 atof(line->words[3]));
      printf("\tAdding a gaussian with %f, %f, %f \n",
             atof(line->words[1]), atof(line->words[2]), atof(line->words[3]));
      phen_target_change = true;
    }
    else if (strcmp(line->words[0], "ENV_ADD_POINT") == 0) {
      // custom_points
      printf("%s:%d: error: Custom points_ management has been removed.\n",
             __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    else if (strcmp(line->words[0], "ENV_VARIATION") == 0) {
      // TODO <dpa> adapt to new organization
//      printf("%s:%d: error: ENV_VARIATION has to be adapted to the new organization.\n", __FILE__, __LINE__);
//      exit(EXIT_FAILURE);

      static bool env_var_already_set = false;
      if (env_var_already_set) {
        printf("%s:%d: ERROR in param file : duplicate entry for %s.\n",
               __FILE__, __LINE__, line->words[0]);
        exit(EXIT_FAILURE);
      }
      env_var_already_set = true;

      if (strcmp(line->words[1], "none") == 0) {
        assert(line->nb_words == 2);
        exp_manager->world()->phenotypic_target_handler()->set_var_method(
            NO_VAR);
        printf("\tNo more environmental variation\n");
      }
      else if (strcmp(line->words[1], "autoregressive_mean_variation") == 0) {
        assert(line->nb_words == 5);
        auto pt_handler = exp_manager->world()->phenotypic_target_handler();
        pt_handler->set_var_method(AUTOREGRESSIVE_MEAN_VAR);
        pt_handler->set_var_sigma_tau(atof(line->words[2]),
                                      atol(line->words[3]));
        pt_handler->set_var_prng(
            std::make_shared<JumpingMT>(atoi(line->words[4])));
        printf(
            "\tChange of environmental variation to a autoregressive mean variation with sigma=%f, tau=%ld and seed=%d\n",
            atof(line->words[2]), atol(line->words[3]), atoi(line->words[4]));
      }
      else if (strcmp(line->words[1], "autoregressive_height_variation") == 0) {
        assert(line->nb_words == 5);
        auto pt_handler = exp_manager->world()->phenotypic_target_handler();
        pt_handler->set_var_method(AUTOREGRESSIVE_HEIGHT_VAR);
        pt_handler->set_var_sigma_tau(atof(line->words[2]),
                                      atol(line->words[3]));
        pt_handler->set_var_prng(
            std::make_shared<JumpingMT>(atoi(line->words[4])));
        printf(
            "\tChange of environmental variation to a autoregressive height variation with sigma=%f, tau=%ld and seed=%d\n",
            atof(line->words[2]), atol(line->words[3]), atoi(line->words[4]));
      }
      else if (strcmp(line->words[1], "add_local_gaussians") == 0) {
        assert(line->nb_words == 3);
        auto pt_handler = exp_manager->world()->phenotypic_target_handler();
        pt_handler->set_var_method(LOCAL_GAUSSIANS_VAR);
        pt_handler->set_var_prng(
            std::make_shared<JumpingMT>(atoi(line->words[2])));
        printf(
            "\tChange of environmental variation to a local gaussians variation with seed=%d\n",
            atoi(line->words[2]));
      }
      else {
        Utils::ExitWithUsrMsg("unknown environment variation method");
      }
    }
    else if (strcmp(line->words[0], "SECRETION_CONTRIB_TO_FITNESS") == 0) {
      exp_manager->exp_s()->set_secretion_contrib_to_fitness(
          atof(line->words[1]));
    }
    else if (strcmp(line->words[0], "SECRETION_COST") == 0) {
      exp_manager->exp_s()->set_secretion_cost(atof(line->words[1]));
    }
    else if (strcmp(line->words[0], "PLASMID_MINIMAL_LENGTH") == 0) {
      if (not exp_manager->with_plasmids()) {
        printf(
            "ERROR: option PLASMID_MINIMAL_LENGTH has no sense because there are no plasmids in this population.\n");
        exit(EXIT_FAILURE);
      }
      int32_t plasmid_minimal_length = atoi(line->words[1]);
      for (const auto& indiv: exp_manager->indivs()) {
        if (indiv->genetic_unit(1).seq_length() < plasmid_minimal_length) {
          printf(
              "ERROR: there is one genetic unit with a smaller length than the new minimum.\n");
          exit(EXIT_FAILURE);
        }
        indiv->genetic_unit_nonconst(1).set_min_gu_length(
            plasmid_minimal_length);
      }
    }
    else if (strcmp(line->words[0], "PLASMID_MAXIMAL_LENGTH") == 0) {
      if (!exp_manager->with_plasmids()) {
        printf(
            "ERROR: option PLASMID_MAXIMAL_LENGTH has no sense because there are no plasmids in this population.\n");
        exit(EXIT_FAILURE);
      }
      int32_t plasmid_maximal_length = atoi(line->words[1]);
      for (const auto& indiv: exp_manager->indivs()) {
        if (indiv->genetic_unit_nonconst(1).seq_length() >
            plasmid_maximal_length) {
          printf(
              "ERROR: there is one genetic unit with a higher length than the new maximum.\n");
          exit(EXIT_FAILURE);
        }
        indiv->genetic_unit_nonconst(1).set_max_gu_length(
            plasmid_maximal_length);
      }
    }
    else if (strcmp(line->words[0], "CHROMOSOME_MINIMAL_LENGTH") == 0) {
      int32_t chromosome_minimal_length = atoi(line->words[1]);
      for (const auto& indiv: exp_manager->indivs()) {
        if (indiv->genetic_unit_nonconst(0).seq_length() <
            chromosome_minimal_length) {
          printf(
              "ERROR: there is one genetic unit with a smaller length than the new minimum.\n");
          exit(EXIT_FAILURE);
        }
        indiv->genetic_unit_nonconst(0).set_min_gu_length(
            chromosome_minimal_length);
      }
    }
    else if (strcmp(line->words[0], "CHROMOSOME_MAXIMAL_LENGTH") == 0) {
      int32_t chromosome_maximal_length = atoi(line->words[1]);
      for (const auto& indiv: exp_manager->indivs()) {
        if (indiv->genetic_unit_nonconst(0).seq_length() >
            chromosome_maximal_length) {
          printf(
              "ERROR: there is one genetic unit with a higher length than the new maximum.\n");
          exit(EXIT_FAILURE);
        }
        indiv->genetic_unit_nonconst(0).set_max_gu_length(
            chromosome_maximal_length);
      }
    }
    else if (strcmp(line->words[0], "SEED") == 0) {
      int32_t seed = atoi(line->words[1]);
      std::shared_ptr<JumpingMT> prng = std::make_shared<JumpingMT>(seed);

      for (int16_t x = 0; x < world->width(); x++) {
        for (int16_t y = 0; y < world->height(); y++) {
          int32_t seed = prng->random(1000000);
#if __cplusplus == 201103L
          world->grid(x,y)->set_reprod_prng(make_unique<JumpingMT>(seed));
          world->grid(x,y)->set_reprod_prng_simd(make_unique<JumpingMT>(seed));
          world->grid(x,y)->set_mut_prng(std::make_shared<JumpingMT>(seed));
          world->grid(x,y)->set_stoch_prng(std::make_shared<JumpingMT>(seed));
#else
          world->grid(x,y)->set_reprod_prng(std::make_unique<JumpingMT>(seed));
          world->grid(x,y)->set_reprod_prng_simd(std::make_unique<JumpingMT>(seed));
          world->grid(x,y)->set_mut_prng(std::make_shared<JumpingMT>(seed));
          world->grid(x,y)->set_stoch_prng(std::make_shared<JumpingMT>(seed));
#endif
        }
      }

      printf("\tChange of the seed to %d in selection and world \n",
             atoi(line->words[1]));
    }
    else if (strcmp(line->words[0], "MUT_SEED") == 0) {
      int32_t mut_seed = atoi(line->words[1]);

      // Change mutation prng
      world->set_mut_prng(std::make_shared<JumpingMT>(mut_seed));
      for (int16_t x = 0; x < world->width(); x++) {
        for (int16_t y = 0; y < world->height(); y++) {
          world->grid(x, y)->set_mut_prng(
              std::make_shared<JumpingMT>(world->mut_prng()->random(1000000)));
          world->grid(x, y)->individual()->set_mut_prng(
              world->grid(x, y)->mut_prng());
        }
      }
      printf("\tChange of the seed to %d in mutations \n",
             atoi(line->words[1]));
    }
    else if (strcmp(line->words[0], "STOCH_SEED") == 0) {
      int32_t stoch_seed = atoi(line->words[1]);

      // Change stochasticity prng
      world->set_stoch_prng(std::make_shared<JumpingMT>(stoch_seed));
      for (int16_t x = 0; x < world->width(); x++) {
        for (int16_t y = 0; y < world->height(); y++) {
          world->grid(x, y)->set_stoch_prng(std::make_shared<JumpingMT>(
              world->stoch_prng()->random(1000000)));
          world->grid(x, y)->individual()->set_stoch_prng(
              world->grid(x, y)->stoch_prng());
        }
      }
      printf("\tChange of the seed to %d in individuals' stochasticity \n",
             atoi(line->words[1]));
    }else if (strcmp(line->words[0], "SIMD_METADATA_FLAVOR") == 0)
    {
      if (strncmp(line->words[1], "stdmap", 6) == 0)
      {
        exp_manager->exp_s()->set_simd_metadata_flavor(STD_MAP);
      }
      else if (strncmp(line->words[1], "dyntab", 6) == 0)
      {
        exp_manager->exp_s()->set_simd_metadata_flavor(DYN_TAB);
      }
      else if (strncmp(line->words[1], "list", 6) == 0)
      {
        exp_manager->exp_s()->set_simd_metadata_flavor(STD_LIST);
      }
    }
    else if (strcmp(line->words[0], "CLONE_BEST") == 0) {
      exp_manager->FillGridWithClones(*(exp_manager->best_indiv()));
      printf("\tChange of the population for a population with %" PRId32
      " individuals, all clones of the best one\n",
          exp_manager->nb_indivs());
    }
    else if (strcmp(line->words[0], "WORLD_SIZE") == 0) {
      #ifdef HAVE_MPI  
        printf("Do not change world_size but GLOBAL_WORLD_SIZE for MPI simulation.\n");
        exit(-42);
      #endif
        // Init Factory (Fuzzy/Dna)
      DnaFactory* dna_factory_ = new DnaFactory(DnaFactory_Policy::FIRSTFIT,32,5000);
      FuzzyFactory_7* fuzzy_factory_ = new FuzzyFactory_7(exp_manager->exp_s()->get_fuzzy_flavor(),exp_manager->nb_indivs()*4,
                            exp_manager->world()->phenotypic_target_handler()->sampling());

      // Change population size
      grid_width_ = atoi(line->words[1]);
      grid_height_ = atoi(line->words[2]);
      int32_t pop_size = grid_width_ * grid_height_;

      World* old_world = exp_manager->non_const_world();
      int32_t old_pop_size = old_world->width() * old_world->height();
      // Init new world
      World* new_world = new World();
      new_world->set_prng(old_world->non_const_prng());
      new_world->set_mut_prng(old_world->non_const_mut_prng());
      new_world->set_stoch_prng(old_world->non_const_stoch_prng());
      new_world->InitGrid(grid_width_, grid_height_, old_world->grid(0,0)->habitat_non_const(), true);

      

      // Reinit PRNG
      #ifdef HAVE_MPI  
        printf("NOT YET AVAILABLE\n");
        exit(-42);
        // delete [] exp_manager->prng_seed_;
        // exp_manager->prng_seed_ = new int32_t[global_pop_size_];

        // for (int16_t x = 0; x < global_grid_width_; x++) {
        //   for (int16_t y = 0; y < global_grid_height_; y++) {
        //     exp_m->prng_seed_[x*global_grid_height_+y] = prng_->random(1000000);
        //   }
        // }
      #endif

        for (int16_t x = 0; x < grid_width_; x++) {
          for (int16_t y = 0; y < grid_height_; y++) {
            #ifdef HAVE_MPI
            int32_t seed = exp_m->prng_seed_[x*exp_m->grid_height()+y];
            #else
            int32_t seed = new_world->prng()->random(1000000);
            #endif
      #if __cplusplus == 201103L
            new_world->grid(x,y)->set_reprod_prng(make_unique<JumpingMT>(seed));
            new_world->grid(x,y)->set_reprod_prng_simd(make_unique<JumpingMT>(seed));
      #else
            new_world->grid(x,y)->set_reprod_prng(std::make_unique<JumpingMT>(seed));
            new_world->grid(x,y)->set_reprod_prng_simd(std::make_unique<JumpingMT>(seed));
      #endif
          }
        }
      // Select and Place individuals
      double *  local_fit_array   = new double[old_pop_size];
      double *  probs             = new double[old_pop_size];
      double    sum_local_fit     = 0.0;

      double selection_pressure = exp_manager->selection_pressure();

      for (int loc_indiv_id = 0; loc_indiv_id < old_pop_size; loc_indiv_id++) {
        int32_t x = loc_indiv_id / old_world->height();
        int32_t y = loc_indiv_id % old_world->height();

        Individual_7* indiv = new Individual_7(exp_manager, old_world->grid(x,y)->individual()->w_max(),dna_factory_,fuzzy_factory_);
        indiv->dna_ = dna_factory_->get_dna(old_world->grid(x,y)->individual()->genetic_unit_seq_length(0));
        indiv->dna_->set_indiv(old_world->grid(x,y)->individual()->genetic_unit(0).dna(),dna_factory_);
        indiv->dna_->set_indiv(indiv);
        indiv->indiv_id = 0;
        indiv->parent_id = 0;
      
        exp_manager->exp_m_7_->evaluate(indiv,old_world->grid(x,y)->individual()->w_max(),selection_pressure);

        local_fit_array[loc_indiv_id] =
                indiv->fitness;
        sum_local_fit += local_fit_array[loc_indiv_id];

        delete indiv;
      }

      for (int loc_indiv_id = 0; loc_indiv_id < old_pop_size; loc_indiv_id++) {
        probs[loc_indiv_id] = local_fit_array[loc_indiv_id]/sum_local_fit;
      }
      
      int32_t* nb_offsprings = new int32_t[old_pop_size];
      exp_manager->world()->grid(0,0)->reprod_prng_simd_->multinomial_drawing(nb_offsprings, probs, pop_size, old_pop_size);

      int64_t index = 0;

      std::ofstream match_pop_change;
      match_pop_change.open("match_pop_change.csv", std::ofstream::trunc);
      match_pop_change << "Generation,old_index,new_index"<<std::endl;

      for (int32_t loc_indiv = 0; loc_indiv < old_pop_size; loc_indiv++) {
        for (int32_t j = 0; j < nb_offsprings[loc_indiv]; j++) {
          int32_t x = index / grid_height_;
          int32_t y = index % grid_height_;
        
          int32_t old_x = loc_indiv / old_world->height();
          int32_t old_y = loc_indiv % old_world->height();

          match_pop_change << AeTime::time() <<","<<loc_indiv<<","<<index<<std::endl;

          #ifdef __NO_X
            #ifndef __REGUL
            Individual* new_indiv =
                new Individual(old_world->indiv_at(old_x,old_y), index, new_world->grid(x, y)->mut_prng(),
                              new_world->grid(x, y)->stoch_prng());
            #else
            Individual_R* new_indiv =
                new Individual_R(dynamic_cast<Individual_R*>(old_world->indiv_at(old_x,old_y)), index,
                                new_world->grid(x, y)->mut_prng(),
                                new_world->grid(x, y)->stoch_prng());
            #endif
          #elif defined __X11
            #ifndef __REGUL
            Individual_X11* new_indiv =
                new Individual_X11(dynamic_cast<Individual_X11*>(old_world->indiv_at(old_x,old_y)), index,
                                  new_world->grid(x, y)->mut_prng(),
                                  new_world->grid(x, y)->stoch_prng());
            #else
            Individual_R_X11* new_indiv =
                new Individual_R_X11(dynamic_cast<Individual_R_X11*>(old_world->indiv_at(old_x,old_y)), index,
                                    new_world->grid(x, y)->mut_prng(),
                                    new_world->grid(x, y)->stoch_prng());
            #endif
          #endif

          new_world->PlaceIndiv(new_indiv, x, y, true);

          index++;
        }
      }

      delete [] nb_offsprings;
      delete [] local_fit_array;
      delete [] probs;



      // Update ExpManager
      exp_manager->set_world(new_world);
      delete old_world;

            // Init new Tree
      exp_manager->output_m()->init_tree(exp_manager,exp_manager->output_m()->tree()->tree_step());
      tree = exp_manager->tree();
      // delete [] exp_manager->next_generation_reproducer_;
      // next_generation_reproducer_ = new int32_t[exp_mananb_indivs()];
    }

    // TODO: re-enable these options
    // else if (strcmp(line->words[0], "CREATE_3_SUBPOPULATIONS_BASED_ON_NON_CODING_BASES") == 0)
    // {
    //   change_based_on_non_coding_bases_of_best_individual(pop, exp_manager, SUBPOPULATIONS_BASED_ON_NON_CODING_BASES);
    //   printf("\tChange of the population for a population with %" PRId32 " individuals in 3 equal subpopulations (A: clones of the previous best individual, B: clones of the previous best individual without any non coding bases, C: clones of the previous best individual with twice non bases\n",pop->nb_indivs());
    //   printf("WARNING: lineage will not work properly if called with \n");
    //   printf("         a begin generation anterior to this modification \n");
    // }
    // else if (strcmp(line->words[0], "REMOVE_NON_CODING_BASES_BEST") == 0)
    // {
    //   change_based_on_non_coding_bases_of_best_individual(pop, exp_manager, REMOVE_NON_CODING_BASES_BEST_IND);
    //   printf("\tChange of the population for a population with %" PRId32 " clones of the best individual ancestor without any non coding bases\n",pop->nb_indivs());
    //   printf("WARNING: lineage will not work properly if called with \n");
    //   printf("         a begin generation anterior to this modification \n");
    // }
    // else if (strcmp(line->words[0], "REMOVE_NON_CODING_BASES_POP") == 0)
    // {
    //   change_based_on_non_coding_bases_in_population(pop, exp_manager,  REMOVE_NON_CODING_BASES_POPULATION);
    //   printf("\tChange of the population for a population with %" PRId32 " individuals without any non coding bases\n",pop->nb_indivs());
    //   printf("WARNING: lineage will not work properly if called with \n");
    //   printf("         a begin generation anterior to this modification \n");
    // }
    // else if (strcmp(line->words[0], "DOUBLE_NON_CODING_BASES_BEST") == 0)
    // {
    //   change_based_on_non_coding_bases_of_best_individual(pop, exp_manager,  DOUBLE_NON_CODING_BASES_BEST_IND);
    //   printf("\tChange of the population for a population with %" PRId32 " clones of the best individual ancestor with twice the non coding bases number \n",pop->nb_indivs());
    //   printf("WARNING: lineage will not work properly if called with \n");
    //   printf("         a begin generation anterior to this modification \n");
    // }
    // else if (strcmp(line->words[0], "DOUBLE_NON_CODING_BASES_POP") == 0)
    // {
    //   change_based_on_non_coding_bases_in_population(pop, exp_manager, DOUBLE_NON_CODING_BASES_POPULATION);
    //   printf("\tChange of the population for a population with %" PRId32 " individuals with twice the non coding bases number\n",pop->nb_indivs());
    //   printf("WARNING: lineage will not work properly if called with \n");
    //   printf("         a begin generation anterior to this modification \n");
    // }
    else {
      printf("%s:%d: error: the change %s is not implemented yet \n", __FILE__,
             __LINE__, line->words[0]);
      exit(EXIT_FAILURE);
    }

    delete line;
  }
  fclose(param_file);

  printf("OK\n");

  if (phen_target_change) {
    // The current version doesn't allow for phenotypic variation nor for
    // different phenotypic targets among the grid
    if (not exp_manager->world()->phenotypic_target_shared()) {
      Utils::ExitWithUsrMsg("sorry, aevol_modify has not yet been implemented "
                                "for per grid-cell phenotypic target");
    }
    auto phenotypicTargetHandler =
        exp_manager->world()->phenotypic_target_handler();
    phenotypicTargetHandler->set_gaussians(new_gaussians);
    phenotypicTargetHandler->BuildPhenotypicTarget();
  }

  // 9) Save the modified experiment
  if (start_to_record_tree) {
    if (!set_tree_step) {
      printf("WARNING: you modifed parameter RECORD_TREE without specifying "
                 "TREE_STEP in the same parameter modification file. TREE_STEP will "
                 "be set to its default value even if you previously gave another "
                 "value.\n");
    }
    exp_manager->output_m()->init_tree(exp_manager, tree_step);
  }

  if (take_care_of_the_tree) {
    printf("Save the modified replication reports into tree...\t");

    #ifdef __REGUL
      sprintf(tree_file_name,"tree/tree_" TIMESTEP_FORMAT ".ae", timestep);
    #else
      sprintf(tree_file_name, "tree/tree_" TIMESTEP_FORMAT ".ae", timestep);
    #endif
    //gzFile tree_file = gzopen(tree_file_name, "w");
    tree->write_to_tree_file(tree_file_name);
    //gzclose(tree_file);
    printf("OK\n");
  }
  printf("Save the modified experiment into backup...\t");
  exp_manager->WriteSetupFiles();
  exp_manager->WriteDynamicFiles();
  printf("OK\n");

  delete exp_manager;
}


/*!
  \brief Get a line in a file and format it

  \param param_file file with param in which a line is reading
  \return line (pointer)

  \see format_line(ParameterLine* formatted_line, char* line,
                   bool* line_is_interpretable)
*/
ParameterLine* get_line(FILE* param_file) {
  char line[255];
  ParameterLine* formatted_line = new ParameterLine();

  bool found_interpretable_line = false;

  while (!feof(param_file) && !found_interpretable_line) {
    if (!fgets(line, 255, param_file)) {
      delete formatted_line;
      return nullptr;
    }
    format_line(formatted_line, line, &found_interpretable_line);
  }

  if (found_interpretable_line) {
    return formatted_line;
  }
  else {
    delete formatted_line;
    return nullptr;
  }
}

/*!
  \brief Format a line by parsing it and the words inside

  \param formatted_line the resulted formatted line
  \param line original line in char*
  \param line_is_interpretable boolean with about the possible intrepretation of the line
*/
void format_line(ParameterLine* formatted_line, char* line,
                 bool* line_is_interpretable) {
  int16_t i = 0;
  int16_t j;

  // Parse line
  while (line[i] != '\n' && line[i] != '\0' && line[i] != '\r') {
    j = 0;

    // Flush white spaces and tabs
    while (line[i] == ' ' || line[i] == 0x09)
      i++; // 0x09 is the ASCII code for TAB

    // Check comments
    if (line[i] == '#') break;

    // If we got this far, there is content in the line
    *line_is_interpretable = true;

    // Parse word
    while (line[i] != ' ' && line[i] != '\n' && line[i] != '\0' &&
           line[i] != '\r') {
      formatted_line->words[formatted_line->nb_words][j++] = line[i++];
    }

    // Add '\0' at end of word if it's not empty (line ending with space or tab)
    if (j != 0) {
      formatted_line->words[formatted_line->nb_words++][j] = '\0';
    }
  }
}






// /*!
//   \brief Change in the population based on non coding bases on the best individual. 3 types of changes

//   SUBPOPULATIONS_BASED_ON_NON_CODING_BASES:
//   Create the 3 subpopulations in the population. The definition of 3 subpopulations is based on non coding bases.

//   The subpopulation are clonal and based on the ancestor of best individual of pop at begin.
//   The individuals in first subpopulation are clones of the best individual.
//   The individuals in second subpopulation are clones of the best individual without any bases that are not in coding RNA.
//   The individuals in third subpopulation are clones of the best individual with addition of bases that are not in coding RNA to double them.

//   pop is changed into the new population with the 3 subpopulations

//   REMOVE_NON_CODING_BASES_BEST_IND:
//   The individuals of the new population are clones of the best individual but without any bases that are not in coding RNA.

//   DOUBLE_NON_CODING_BASES_BEST_IND:
//   The individuals of the new population are clones of the best individual but with addition of bases that are not in coding RNA to double them.

//   \param pop population to change
//   \param exp_m global exp_manager
//   \param type type of change in the population
// */
// void change_based_on_non_coding_bases_of_best_individual(ae_population* pop, ExpManager* exp_m, population_change_type type)
// {
//   if(type == SUBPOPULATIONS_BASED_ON_NON_CODING_BASES || type == REMOVE_NON_CODING_BASES_BEST_IND || type == DOUBLE_NON_CODING_BASES_BEST_IND)
//     {
//       // 1) Compute the population size
//       int32_t subpopulation_size = (int)floor(pop->nb_indivs()/3);

//       // 2) Get the best individual
//       ae_individual* best_indiv = exp_m->best_indiv();


//       // 3) Create the new population


//       std::list<ae_individual*> new_generation;

//       ae_individual* indiv = create_clone(best_indiv, -1);

//       ae_individual* only_coding_indiv = create_clone(best_indiv, -1); //one individual being the clone of the chosen individual but without any non coding bases
//       only_coding_indiv->remove_non_coding_bases();

//       ae_individual* twice_non_coding_indiv = create_clone(best_indiv, -1); //one individual being the clone of the chosen individual but without any non coding bases
//       twice_non_coding_indiv->double_non_coding_bases();


//       int32_t* probe_A = new int32_t[5];
//       int32_t* probe_B = new int32_t[5];
//       int32_t* probe_C = new int32_t[5];
//       for(int32_t i = 0 ; i<5; i++)
//         {
//           probe_A[i] = 1;
//           probe_B[i] = 10;
//           probe_C[i] = 100;
//         }
//       indiv->set_int_probes(probe_A);
//       only_coding_indiv->set_int_probes(probe_B);
//       twice_non_coding_indiv->set_int_probes(probe_C);

//       double* probe_double_A = new double[5];
//       double* probe_double_B = new double[5];
//       double* probe_double_C = new double[5];
//       for(int32_t i = 0 ; i<5; i++)
//         {
//           probe_double_A[i] = 1;
//           probe_double_B[i] = 10;
//           probe_double_C[i] = 100;
//         }
//       indiv->set_double_probes(probe_double_A);
//       only_coding_indiv->set_double_probes(probe_double_B);
//       twice_non_coding_indiv->set_double_probes(probe_double_C);


//       switch(type)
//         {
//         case SUBPOPULATIONS_BASED_ON_NON_CODING_BASES:
//           {
//             int32_t  index_new_indiv = 0;
//             for (int32_t i = 0 ; i < subpopulation_size ; i++) // clones of the 3 individuals
//               {
//                 new_generation.push_back(create_clone(indiv, index_new_indiv++));
//                 new_generation.push_back(create_clone(only_coding_indiv, index_new_indiv++));
//                 new_generation.push_back(create_clone(twice_non_coding_indiv, index_new_indiv++));
//               }
//             break;
//           }
//         case REMOVE_NON_CODING_BASES_BEST_IND:
//           {
//             for (int32_t i = 0 ; i < pop->nb_indivs() ; i++)
//               {
//                 new_generation.push_back(create_clone(only_coding_indiv, i));
//               }
//             break;
//           }
//         case DOUBLE_NON_CODING_BASES_BEST_IND:
//           {
//             for (int32_t i = 0 ; i < pop->nb_indivs() ; i++)
//               {
//                 new_generation.push_back(create_clone(twice_non_coding_indiv, i));
//               }
//             break;
//           }
//         default:
//           {
//             fprintf(stderr, "%s:%d: error: wrong population_change_type %d\n", __FILE__, __LINE__, type);
//             exit(EXIT_FAILURE);
//             break;
//           }
//         }

//       //  4) Replace the current population by the new one
//       //     -> Useless since it is done by replace_population.
//       pop->replace_population(std::move(new_generation));



//       // TODO
//       // If the population is spatially structured, set each individual's position
//       // There will be a problem however for the "3 subpopulations" type of change,
//       // if the population size has changed (which is likely given that we do not
//       // generally used population size that are multiple of 3)

//       pop->evaluate_individuals(exp_m->env());
//       pop->sort_individuals();
//     }
//   else
//     {
//       printf("%s:%d: error: wrong population_change_type %d\n", __FILE__, __LINE__, type);
//       exit(EXIT_FAILURE);
//     }
// }

// /*!
//   \brief Change in the population based on non coding bases. 2 types of changes

//   REMOVE_NON_CODING_BASES_POPULATION:
//   The individual of the new population are the individuals without any bases that are not in coding RNA.

//   DOUBLE_NON_CODING_BASES_POPULATION:
//   The individual of the new population are the individuals with addition of bases that are not in coding RNA to double them.

//   \param pop population to change
//   \param exp_m global exp_manager
//   \param type type of change in the population
// */
// void change_based_on_non_coding_bases_in_population(ae_population* pop, ExpManager* exp_m, population_change_type type)
// {
//   if(type == REMOVE_NON_CODING_BASES_POPULATION || type == DOUBLE_NON_CODING_BASES_POPULATION)
//     {
//       for (auto& indiv: pop->indivs())
//         if (type == REMOVE_NON_CODING_BASES_POPULATION)
//           indiv->remove_non_coding_bases();
//         else
//           indiv->double_non_coding_bases();
//     }
//   else
//     {
//       printf("%s:%d: error: wrong population_change_type %d\n", __FILE__, __LINE__, type);
//       exit(EXIT_FAILURE);
//     }
// }







/*!
  \brief

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
  printf("%s: modify an experiment as specified in PARAM_FILE.\n", prog_name);
  printf("\n");
  printf("Usage : %s -h or --help\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("   or : %s [-t TIMESTEP] [-f PARAM_FILE]\n", prog_name);
  printf("\nOptions\n");
  printf("  -h, --help\n"
             "\tprint this help, then exit\n");
  printf("  -V, --version\n"
             "\tprint version number, then exit\n");
  printf("  -t, --timestep TIMESTEP\n"
             "\tspecify timestep\n");
  printf("  -f, --file PARAM_FILE\n"
             "\tspecify parameter file (default: param.in)\n");
}

void interpret_cmd_line_options(int argc, char* argv[]) {
  // Define allowed options
  const char* options_list = "hf:t:V";
  static struct option long_options_list[] = {
      {"help",      no_argument,       nullptr, 'h'},
      {"version",   no_argument,       nullptr, 'V'},
      {"file",      required_argument, nullptr, 'f'},
      {"timestep",  required_argument, nullptr, 't'},
      {0, 0, 0, 0}
  };

  // Get actual values of the CLI options
  int option;
  while ((option = getopt_long(argc, argv, options_list, long_options_list,
                               nullptr)) != -1) {
    switch (option) {
      case 'h' : {
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      }
      case 'V' : {
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      }
      case 'f' : {
        if (strcmp(optarg, "") == 0) {
          printf("%s: error: Option -f or --file : missing argument.\n",
                 argv[0]);
          exit(EXIT_FAILURE);
        }

        param_file_name = optarg;
        break;
      }
      case 't' : {
        timestep = atoi(optarg);
        break;
      }
      default : {
        // An error message is printed in getopt_long, we just need to exit
        exit(EXIT_FAILURE);
      }
    }
  }

  // If param file name wasn't provided, use default
  if (param_file_name == nullptr) {
    param_file_name = new char[strlen(DEFAULT_PARAM_FILE_NAME) + 1];
    sprintf(param_file_name, "%s", DEFAULT_PARAM_FILE_NAME);
  }

  // If timestep wasn't provided, use default
  if (timestep == -1) {
    timestep = OutputManager::last_gener();
  }
}
