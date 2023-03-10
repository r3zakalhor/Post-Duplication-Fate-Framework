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
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cerrno>
#include <climits>
#include <ctime>

#include <list>

// =================================================================
//                            Project Files
// =================================================================
#include "ParamLoader.h"

#include "FuzzyFactory.h"

#if __cplusplus == 201103L
#include "make_unique.h"
#endif

#include "ExpManager.h"
#include "ExpSetup.h"
#include "OutputManager.h"
#include "Individual.h"
#include "IndividualFactory.h"
#include "ExpManager.h"


#ifdef __REGUL
#include "raevol/Individual_R.h"
#endif



#include "JumpingMT.h"
#include "Gaussian.h"
#include "PhenotypicSegment.h"
#include "Point.h"
#include "Alignment.h"
#include "World.h"


namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================


//##############################################################################
//                                                                             #
//                             Class ParamLoader                              #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================
static const int8_t STRAIN_NAME_DEFAULT_SIZE  = 20;
static const int8_t STRAIN_NAME_LOGIN_SIZE    = 10;

const char kTabChar = 0x09;

// =================================================================
//                             Constructors
// =================================================================
ParamLoader::ParamLoader(const char* file_name)
{
  // Give default values to parameters

  // ----------------------------------------- PseudoRandom Number Generators
  seed_           = 0;
  mut_seed_       = 0;
  stoch_seed_     = 0;
  env_var_seed_   = 0;
  env_noise_seed_ = 0;

  // ------------------------------------------------------------ Constraints
  min_genome_length_  = 1;
  max_genome_length_  = 10000000;
  w_max_              = 0.033333333;

  // ----------------------------------------------------- Initial conditions
  chromosome_initial_length_  = 5000;
  init_method_            = ONE_GOOD_GENE | CLONE;
  init_pop_size_          = 1024;
  strain_name_ = new char[STRAIN_NAME_DEFAULT_SIZE+1];

  for (int i = 0; i < STRAIN_NAME_DEFAULT_SIZE+1; i++)
    strain_name_[i] = 0;

  // ------------------------------------------------------------- Strain name
  char* login_name = new char[LOGIN_NAME_MAX+1];
  // Try get user login. If fail, replace by default value
  if(getlogin_r(login_name, LOGIN_NAME_MAX) != 0)
    strcpy(login_name, "anon");

  // Copy login into strain name with at most STRAIN_NAME_LOGIN_SIZE characters
  strncpy(strain_name_, login_name, STRAIN_NAME_LOGIN_SIZE);
  delete [] login_name;

  // Null-terminate the c-string if the max number of characters were copied
  if (strain_name_[STRAIN_NAME_LOGIN_SIZE] != 0)
    strain_name_[STRAIN_NAME_LOGIN_SIZE + 1] = 0;

  // -------------------------------------------------------- Phenotypic target
  env_sampling_ = 300;

  // ------------------------------------ Phenotypic target x-axis segmentation
  env_axis_nb_segments_         = 1;
  env_axis_segment_boundaries_  = NULL;
  env_axis_features_            = NULL;
  env_axis_separate_segments_   = false;

  // ---------------------------------------------- Phenotypic target variation
  env_var_method_ = NO_VAR;
  env_var_sigma_  = 0;
  env_var_tau_    = 0;

  // -------------------------------------------------- Phenotypic target noise
  env_noise_method_       = NO_NOISE;
  env_noise_alpha_        = 0;
  env_noise_sigma_        = 0;
  env_noise_prob_         = 0;
  env_noise_sampling_log_ = 0;

  // --------------------------------------------------------- Mutation rates
  point_mutation_rate_  = 1e-5;
  small_insertion_rate_ = 1e-5;
  small_deletion_rate_  = 1e-5;
  max_indel_size_       = 6;

  // -------------------------------------------- Rearrangements and Transfer
  with_4pts_trans_            = true;
  with_alignments_            = false;
  with_HT_                    = false;
  repl_HT_with_close_points_  = false;
  HT_ins_rate_                = 0.0;
  HT_repl_rate_               = 0.0;
  repl_HT_detach_rate_        = 0.0;

  // ------------------------------ Rearrangement rates (without alignements)
  duplication_rate_   = 1e-5;
  deletion_rate_      = 1e-5;
  translocation_rate_ = 1e-5;
  inversion_rate_     = 1e-5;

  // --------------------------------- Rearrangement rates (with alignements)
  neighbourhood_rate_       = 5e-5;
  duplication_proportion_   = 0.3;
  deletion_proportion_      = 0.3;
  translocation_proportion_ = 0.3;
  inversion_proportion_     = 0.3;

  // ------------------------------------------------------------ Alignements
  align_fun_shape_    = SIGMOID;
  align_sigm_lambda_  = 4;
  align_sigm_mean_    = 50;
  align_lin_min_      = 0;
  align_lin_max_      = 100;

  align_max_shift_      = 20;
  align_w_zone_h_len_   = 50;
  align_match_bonus_    = 1;
  align_mismatch_cost_  = 2;

  // ----------------------------------------------- Phenotypic Stochasticity
  with_stochasticity_ = false;

  // -------------------------------------------------------------- Selection
  selection_scheme_   = RANK_EXPONENTIAL;
  selection_pressure_ = 0.998;

  selection_scope_   = SCOPE_LOCAL;
  selection_scope_x_ = 3;
  selection_scope_y_ = 3;

  fitness_function_ = FITNESS_EXP;
  fitness_function_x_ = 3;
  fitness_function_y_ = 3;
  // -------------------------------------------------------------- Secretion
  with_secretion_               = false;
  secretion_contrib_to_fitness_ = 0;
  secretion_diffusion_prop_     = 0;
  secretion_degradation_prop_   = 0;
  secretion_cost_               = 0;
  secretion_init_               = 0;

  // --------------------------------------------------------------- Plasmids
  allow_plasmids_             = false;
  plasmid_initial_length_     = -1;
  plasmid_initial_gene_       = 0;
  plasmid_minimal_length_     = -1;
  plasmid_maximal_length_     = -1;
  chromosome_minimal_length_  = -1;
  chromosome_maximal_length_  = -1;
  prob_plasmid_HT_            = 0;
  tune_donor_ability_         = 0;
  tune_recipient_ability_     = 0;
  donor_cost_                 = 0;
  recipient_cost_             = 0;
  compute_phen_contrib_by_GU_ = false;
  swap_GUs_         = false;

  // ------------------------------------------------------- Translation cost
  translation_cost_ = 0;

  #ifdef BASE_4
  // ------------------------------------------------------------ Terminators
  term_polya_seq_length_ = 0;
  #endif

  // ---------------------------------------------------------------- Outputs
  stats_            = 0;
  delete_old_stats_ = false;

  // Backups
  backup_step_      = 500;
  big_backup_step_  = 10000;

  // Tree
  record_tree_  = false;
  tree_step_    = 100;

  //LightTree
  record_light_tree_ = false;

  // Dumps
  make_dumps_ = false;
  dump_step_  = 1000;

  // Logs
  logs_ = 0;

  // Other
  more_stats_ = false;

  _fuzzy_flavor = 1;

  simd_metadata_flavor_ = MetadataFlavor::STD_LIST;

#ifdef __REGUL
    // ------------------------------------------------------- Binding matrix
    _binding_zeros_percentage = 75;

    _protein_presence_limit = 1e-2;
    _degradation_rate  = 1;
    _nb_degradation_step  = 10;
    _nb_indiv_age         = 20;
    _with_heredity          = false;

    _hill_shape_n      = 4;
    _hill_shape_theta  = 0.5;
    _hill_shape        = std::pow( _hill_shape_theta, _hill_shape_n );

    _list_eval_step.insert(_nb_indiv_age);
    _env_switch_probability = 0.1;
  #endif

  // Read parameter file
  param_file_name_ = strdup(file_name);
  param_file_  = fopen(param_file_name_,  "r");

  if (param_file_ == NULL)
  {
    printf("ERROR : couldn't open file %s\n", file_name);
    exit(EXIT_FAILURE);
  }

  assert(param_file_);

  read_file();
}

// =================================================================
//                             Destructor
// =================================================================
ParamLoader::~ParamLoader()
{
  free(param_file_name_);
  fclose(param_file_);

  delete [] env_axis_segment_boundaries_;
  delete [] env_axis_features_;
  delete [] strain_name_;
}

// =================================================================
//                            Public Methods
// =================================================================
void ParamLoader::interpret_line(ParameterLine * line, int32_t cur_line)
{
  if (strcmp(line->words[0], "STRAIN_NAME") == 0)
  {
    delete [] strain_name_;
    strain_name_ = new char[strlen(line->words[1])+1];
    strcpy(strain_name_, line->words[1]);
  }
  else if (strcmp(line->words[0], "MIN_TRIANGLE_WIDTH") == 0)
  {
    printf("ERROR in param file \"%s\" on line %" PRId32
               ": %s is no longer a valid option, "
               "its value is fixed to 0.\n",
           param_file_name_, cur_line, line->words[0]);
    exit(EXIT_FAILURE);
  }
  else if (strcmp(line->words[0], "MAX_TRIANGLE_WIDTH") == 0)
  {
    w_max_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "ENV_AXIS_FEATURES") == 0)
  {
    // Set general segmentation data
    env_axis_nb_segments_ = line->nb_words / 2;

    // Set segmentation boundaries
    env_axis_segment_boundaries_ = new double [env_axis_nb_segments_ + 1];
    env_axis_segment_boundaries_[0] = X_MIN;
    for (int16_t i = 1 ; i < env_axis_nb_segments_ ; i++)
    {
      env_axis_segment_boundaries_[i] = atof(line->words[2*i]);
    }
    env_axis_segment_boundaries_[env_axis_nb_segments_] = X_MAX;

    // Set segment features
    env_axis_features_ = new PhenotypicFeature[env_axis_nb_segments_];
    for (int16_t i = 0 ; i < env_axis_nb_segments_ ; i++)
    {
      if (strcmp(line->words[2*i+1], "NEUTRAL") == 0)
      {
        env_axis_features_[i] = NEUTRAL;
      }
      else if (strcmp(line->words[2*i+1], "METABOLISM") == 0)
      {
        env_axis_features_[i] = METABOLISM;
      }
      else if (strcmp(line->words[2*i+1], "SECRETION") == 0)
      {
        with_secretion_ = true;
        env_axis_features_[i] = SECRETION;
      }
      else if (strcmp(line->words[2*i+1], "DONOR") == 0)
      {
        env_axis_features_[i] = DONOR;
      }
      else if (strcmp(line->words[2*i+1], "RECIPIENT") == 0)
      {
        env_axis_features_[i] = RECIPIENT;
      }
      else
      {
        printf("ERROR in param file \"%s\" on line %" PRId32
                   ": unknown axis feature \"%s\".\n",
               param_file_name_, cur_line, line->words[2*i+1]);
        exit(EXIT_FAILURE);
      }
    }
  }
  else if (strcmp(line->words[0], "ENV_SEPARATE_SEGMENTS") == 0)
  {
    env_axis_separate_segments_ = true;
  }
  else if (strcmp(line->words[0], "RECORD_TREE") == 0)
  {
    if (strncmp(line->words[1], "true", 4) == 0)
    {
      record_tree_ = true;
    }
    else if (strncmp(line->words[1], "false", 5) == 0)
    {
      record_tree_ = false;
    }
    else
    {
      printf("ERROR in param file \"%s\" on line %" PRId32
                 ": unknown tree recording option (use true/false).\n",
             param_file_name_, cur_line);
      exit(EXIT_FAILURE);
    }
  }
  else if (strcmp(line->words[0], "TREE_MODE") == 0)
  {
    printf("ERROR in param file \"%s\" on line %" PRId32 ": "
           "Tree mode management has been removed.\n",
           param_file_name_, cur_line);
    exit(EXIT_FAILURE);
  }
  else if (strcmp(line->words[0], "RECORD_LIGHT_TREE") == 0)
  {
    if (strncmp(line->words[1], "true", 4) == 0)
    {
      record_light_tree_ = true;
    }
    else if (strncmp(line->words[1], "false", 5) == 0)
    {
      record_light_tree_ = false;
    }
    else
    {
      printf("ERROR in param file \"%s\" on line %" PRId32
             ": unknown light tree recording option (use true/false).\n",
             param_file_name_, cur_line);
      exit(EXIT_FAILURE);
    }
  }
  else if (strcmp(line->words[0], "MORE_STATS") == 0)
  {
    if (strncmp(line->words[1], "true", 4) == 0)
    {
      more_stats_ = true;
    }
    else if (strncmp(line->words[1], "false", 5) == 0)
    {
      more_stats_ = false;
    }
    else
    {
      printf("ERROR in param file \"%s\" on line %" PRId32
                 ": unknown more stats option (use true/false).\n",
             param_file_name_, cur_line);
      exit(EXIT_FAILURE);
    }
  }
  else if (strcmp(line->words[0], "DUMP_PERIOD") == 0)
  {
    printf("ERROR in param file \"%s\" on line %" PRId32
               ": %s is no longer a valid option, "
               "use DUMP_STEP instead.\n",
           param_file_name_, cur_line, line->words[0]);
    exit(EXIT_FAILURE);
  }
  else if (strcmp(line->words[0], "DUMP_STEP") == 0)
  {
    dump_step_ = atol(line->words[1]);
    if (dump_step_>0) make_dumps_ = true;
  }
  else if (strcmp(line->words[0], "BACKUP_STEP") == 0)
  {
    backup_step_ = atol(line->words[1]);
  }
  else if (strcmp(line->words[0], "BIG_BACKUP_STEP") == 0)
  {
    big_backup_step_ = atol(line->words[1]);
  }
  else if (strcmp(line->words[0], "TREE_STEP") == 0)
  {
    tree_step_ = atol(line->words[1]);
  }
  else if (strcmp(line->words[0], "NB_GENER") == 0)
  {
    printf("ERROR in param file \"%s\" on line %" PRId32
               ": %s is no longer a valid option, "
               "use command line arguments of aevol_run instead.\n",
           param_file_name_, cur_line, line->words[0]);
    exit(EXIT_FAILURE);
  }
  else if (strcmp(line->words[0], "INITIAL_GENOME_LENGTH") == 0)
  {
    printf("ERROR in param file \"%s\" on line %" PRId32
               ": %s is no longer a valid option, "
               "use CHROMOSOME_INITIAL_LENGTH (and optionally "
               "PLASMID_INITIAL_LENGTH) instead.\n",
           param_file_name_, cur_line, line->words[0]);
    exit(EXIT_FAILURE);
  }
  else if (strcmp(line->words[0], "CHROMOSOME_INITIAL_LENGTH") == 0)
  {
    chromosome_initial_length_ = atol(line->words[1]);
  }
  else if (strcmp(line->words[0], "MIN_GENOME_LENGTH") == 0)
  {
    if (strncmp(line->words[1], "NONE", 4) == 0)
    {
      min_genome_length_ = 1; // Must not be 0
    }
    else
    {
      min_genome_length_ = atol(line->words[1]);
      if (min_genome_length_ == 0)
      {
        printf("ERROR in param file \"%s\" on line %" PRId32 " : MIN_GENOME_LENGTH must be > 0.\n",
               param_file_name_, cur_line);
        exit(EXIT_FAILURE);
      }
    }
  }
  else if (strcmp(line->words[0], "MAX_GENOME_LENGTH") == 0)
  {
    if (strncmp(line->words[1], "NONE", 4) == 0)
    {
      max_genome_length_ = INT32_MAX;
    }
    else
    {
      max_genome_length_ = atol(line->words[1]);
    }
  }
  else if (strcmp(line->words[0], "INIT_POP_SIZE") == 0)
  {
    init_pop_size_ = atol(line->words[1]);
  }
  else if (strcmp(line->words[0], "WORLD_SIZE") == 0)
  {
    grid_width_ = atoi(line->words[1]);
    grid_height_ = atoi(line->words[2]);
  }
  else if (strcmp(line->words[0], "POP_STRUCTURE") == 0)
  {
    printf("ERROR in param file \"%s\" on line %" PRId32
               ": %s is no longer a valid option, "
               "use WORLD_SIZE <width> <height> instead.\n",
           param_file_name_, cur_line, line->words[0]);
    exit(EXIT_FAILURE);
  }
  else if (strcmp(line->words[0], "MIGRATION_NUMBER") == 0)
  {
    printf("ERROR in param file \"%s\" on line %" PRId32
               ": %s is no longer a valid option, "
               "use INDIV_MIXING instead.\n",
           param_file_name_, cur_line, line->words[0]);
    printf("usage: INDIV_MIXING WELL_MIXED|NONE|PARTIAL <n>\n");
    exit(EXIT_FAILURE);
  }
  else if (strcmp(line->words[0], "INDIV_MIXING") == 0)
  {
    if (strcmp(line->words[1], "WELL_MIXED") == 0)
      well_mixed = true;
    else if (strcmp(line->words[1], "NONE") == 0)
      well_mixed = false;
    else if (strcmp(line->words[1], "PARTIAL") == 0)
      partial_mix_nb_permutations = atol(line->words[2]);
    else {
      printf("ERROR in param file \"%s\" on line %" PRId32
                 ": unknown mixing option.\n", param_file_name_, cur_line);
      printf("usage: INDIV_MIXING WELL_MIXED|NONE|PARTIAL <n>\n");
      exit(EXIT_FAILURE);
    }
  }
  else if (strcmp(line->words[0], "INIT_METHOD") == 0)
  {
    for (int8_t i = 1 ; i < line->nb_words ; i++)
    {
      if (strcmp(line->words[i], "ONE_GOOD_GENE") == 0)
      {
        init_method_ |= ONE_GOOD_GENE;
      }
      else if (strcmp(line->words[i], "CLONE") == 0)
      {
        init_method_ |= CLONE;
      }
      else if (strcmp(line->words[i], "WITH_INS_SEQ") == 0)
      {
        init_method_ |= WITH_INS_SEQ;
      }
      else
      {
        printf("ERROR in param file \"%s\" on line %" PRId32
                   ": unknown initialization method %s.\n",
               param_file_name_, cur_line, line->words[1]);
        exit(EXIT_FAILURE);
      }
    }
  }
  else if (strcmp(line->words[0], "POINT_MUTATION_RATE") == 0)
  {
    point_mutation_rate_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "SMALL_INSERTION_RATE") == 0)
  {
    small_insertion_rate_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "SMALL_DELETION_RATE") == 0)
  {
    small_deletion_rate_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "MAX_INDEL_SIZE") == 0)
  {
    max_indel_size_ = atol(line->words[1]);
  }
  else if (strcmp(line->words[0], "DUPLICATION_RATE") == 0)
  {
    duplication_rate_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "DELETION_RATE") == 0)
  {
    deletion_rate_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "TRANSLOCATION_RATE") == 0)
  {
    translocation_rate_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "INVERSION_RATE") == 0)
  {
    inversion_rate_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "NEIGHBOURHOOD_RATE") == 0)
  {
    neighbourhood_rate_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "DUPLICATION_PROPORTION") == 0)
  {
    duplication_proportion_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "DELETION_PROPORTION") == 0)
  {
    deletion_proportion_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "TRANSLOCATION_PROPORTION") == 0)
  {
    translocation_proportion_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "INVERSION_PROPORTION") == 0)
  {
    inversion_proportion_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "ALIGN_FUNCTION") == 0)
  {
    if (line->nb_words != 2 && line->nb_words != 4)
    {
      printf("ERROR in param file \"%s\" on line %" PRId32
                 ": incorrect number of parameters for keyword \"%s\".\n",
             param_file_name_, cur_line, line->words[0]);
      exit(EXIT_FAILURE);
    }

    if (strcmp(line->words[1], "LINEAR") == 0)
    {
      align_fun_shape_ = LINEAR;

      if (line->nb_words == 4)
      {
        align_lin_min_ = atol(line->words[2]);
        align_lin_max_ = atol(line->words[3]);
      }
    }
    else if (strcmp(line->words[1], "SIGMOID") == 0)
    {
      align_fun_shape_ = SIGMOID;

      if (line->nb_words == 4)
      {
        align_sigm_lambda_ = atol(line->words[2]);
        align_sigm_mean_ = atol(line->words[3]);
      }
    }
    else
    {
      printf("ERROR in param file \"%s\" on line %" PRId32
                 ": unknown align function shape \"%s\".\n",
             param_file_name_, cur_line, line->words[1]);
      exit(EXIT_FAILURE);
    }
  }
  else if (strcmp(line->words[0], "ALIGN_MAX_SHIFT") == 0)
  {
    align_max_shift_ = atol(line->words[1]);
  }
  else if (strcmp(line->words[0], "ALIGN_W_ZONE_H_LEN") == 0)
  {
    align_w_zone_h_len_ = atol(line->words[1]);
  }
  else if (strcmp(line->words[0], "ALIGN_MATCH_BONUS") == 0)
  {
    align_match_bonus_ = atol(line->words[1]);
  }
  else if (strcmp(line->words[0], "ALIGN_MISMATCH_COST") == 0)
  {
    align_mismatch_cost_ = atol(line->words[1]);
  }
  else if (strcmp(line->words[0], "STOCHASTICITY") == 0)
  {
    if (strncmp(line->words[1], "true", 4) == 0)
    {
      with_stochasticity_ = true;
    }
  }
  else if (strcmp(line->words[0], "SELECTION_SCHEME") == 0)
  {
    if (strncmp(line->words[1], "lin", 3) == 0)
    {
      if (line->nb_words != 3)
      {
        printf("ERROR in param file \"%s\" on line %" PRId32
                   ": selection pressure parameter is missing.\n",
               param_file_name_, cur_line);
        exit(EXIT_FAILURE);
      }

      selection_scheme_ = RANK_LINEAR;
      selection_pressure_ = atof(line->words[2]);
    }
    else if (strncmp(line->words[1], "exp", 3) == 0)
    {
      if (line->nb_words != 3)
      {
        printf("ERROR in param file \"%s\" on line %" PRId32
                   ": selection pressure parameter is missing.\n",
               param_file_name_, cur_line);
        exit(EXIT_FAILURE);
      }

      selection_scheme_ = RANK_EXPONENTIAL;
      selection_pressure_ = atof(line->words[2]);
    }
    else if (strncmp(line->words[1], "fitness", 7) == 0)
    {
      if (line->nb_words != 3)
      {
        printf("ERROR in param file \"%s\" on line %" PRId32
                   ": selection pressure parameter is missing.\n",
               param_file_name_, cur_line);
        exit(EXIT_FAILURE);
      }

      selection_scheme_ = FITNESS_PROPORTIONATE;
      selection_pressure_ = atof(line->words[2]);
    }
    else if (strcmp(line->words[1], "fittest") == 0)
    {
      selection_scheme_ = FITTEST;
    }
    else
    {
      printf("ERROR in param file \"%s\" on line %" PRId32
                 ": unknown selection scheme \"%s\".\n",
             param_file_name_, cur_line, line->words[1]);
      exit(EXIT_FAILURE);
    }
  }
  else if (strcmp(line->words[0], "SELECTION_SCOPE") == 0)
  {
    if (strncmp(line->words[1], "global", 6) == 0)
    {
      selection_scope_ = SCOPE_GLOBAL;
    }
    else if (strncmp(line->words[1], "local", 5) == 0)
    {
      if (line->nb_words != 4)
      {
        printf("ERROR in param file \"%s\" on line %" PRId32
                   ": selection scope parameter for local selection is missing (x,y).\n",
               param_file_name_, cur_line);
        exit(EXIT_FAILURE);
      }

      selection_scope_ = SCOPE_LOCAL;
      selection_scope_x_ = atoi(line->words[2]);
      selection_scope_y_ = atoi(line->words[2]);
    }
    else
    {
      printf("ERROR in param file \"%s\" on line %" PRId32
                 ": unknown selection scope \"%s\".\n",
             param_file_name_, cur_line, line->words[1]);
      exit(EXIT_FAILURE);
    }
  }
  else if (strcmp(line->words[0], "SEED") == 0)
  {
    static bool seed_already_set = false;
    if (seed_already_set)
    {
      printf("ERROR in param file \"%s\" on line %" PRId32
                 ": duplicate entry for SEED.\n",
             param_file_name_, cur_line);
      exit(EXIT_FAILURE);
    }
    seed_ = atol(line->words[1]);
    seed_already_set = true;
  }
  else if (strcmp(line->words[0], "MUT_SEED") == 0)
  {
    static bool mut_seed_already_set = false;
    if (mut_seed_already_set)
    {
      printf("ERROR in param file \"%s\" on line %" PRId32
                 ": duplicate entry for MUT_SEED.\n",
             param_file_name_, cur_line);
      exit(EXIT_FAILURE);
    }
    mut_seed_ = atol(line->words[1]);
    mut_seed_already_set = true;
  }
  else if (strcmp(line->words[0], "STOCH_SEED") == 0)
  {
    static bool stoch_seed_already_set = false;
    if (stoch_seed_already_set)
    {
      printf("ERROR in param file \"%s\" on line %" PRId32
                 ": duplicate entry for STOCH_SEED.\n",
             param_file_name_, cur_line);
      exit(EXIT_FAILURE);
    }
    stoch_seed_ = atol(line->words[1]);
    stoch_seed_already_set = true;
  }
  else if (strcmp(line->words[0], "WITH_4PTS_TRANS") == 0)
  {
    if (strncmp(line->words[1], "true", 4) == 0)
    {
      with_4pts_trans_ = true;
    }
    else if (strncmp(line->words[1], "false", 5) == 0)
    {
      printf("ERROR: 3 points_ translocation hasn't been implemented yet\n");
      exit(EXIT_FAILURE);
    }
    else
    {
      printf("ERROR in param file \"%s\" on line %" PRId32
                 ": unknown 4pts_trans option (use true/false).\n",
             param_file_name_, cur_line);
      exit(EXIT_FAILURE);
    }
  }
  else if (strcmp(line->words[0], "WITH_ALIGNMENTS") == 0)
  {
    if (strncmp(line->words[1], "true", 4) == 0)
    {
      with_alignments_ = true;
    }
    else if (strncmp(line->words[1], "false", 5) == 0)
    {
      with_alignments_ = false;
    }
    else
    {
      printf("ERROR in param file \"%s\" on line %" PRId32
                 ": unknown alignement option (use true/false).\n",
             param_file_name_, cur_line);
      exit(EXIT_FAILURE);
    }
  }
  else if (strcmp(line->words[0], "WITH_TRANSFER") == 0)
  {
    if (strncmp(line->words[1], "true", 4) == 0)
    {
      with_HT_ = true;
    }
    else if (strncmp(line->words[1], "false", 5) == 0)
    {
      with_HT_ = false;
    }
    else
    {
      printf("ERROR in param file \"%s\" on line %" PRId32
                 ": unknown transfer option (use true/false).\n",
             param_file_name_, cur_line);
      exit(EXIT_FAILURE);
    }
  }
  else if (strcmp(line->words[0], "REPL_TRANSFER_WITH_CLOSE_POINTS") == 0)
  {
    if (strncmp(line->words[1], "true", 4) == 0)
    {
      repl_HT_with_close_points_ = true;
    }
    else if (strncmp(line->words[1], "false", 5) == 0)
    {
      repl_HT_with_close_points_ = false;
    }
    else
    {
      printf("ERROR in param file \"%s\" on line %" PRId32
                 ": unknown transfer option (use true/false).\n",
             param_file_name_, cur_line);
      exit(EXIT_FAILURE);
    }
  }
  else if (strcmp(line->words[0], "SWAP_GUS") == 0)
  {
    if (strncmp(line->words[1], "true", 4) == 0)
    {
      swap_GUs_ = true;
    }
    else if (strncmp(line->words[1], "false", 5) == 0)
    {
      swap_GUs_ = false;
    }
    else
    {
      printf("ERROR in param file \"%s\" on line %" PRId32
                 ": unknown swap option (use true/false).\n",
             param_file_name_, cur_line);
      exit(EXIT_FAILURE);
    }
  }
  else if (strcmp(line->words[0], "TRANSFER_INS_RATE") == 0)
  {
    HT_ins_rate_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "TRANSFER_REPL_RATE") == 0)
  {
    HT_repl_rate_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "REPL_TRANSFER_DETACH_RATE") == 0)
  {
    repl_HT_detach_rate_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "TRANSLATION_COST") == 0)
  {
    translation_cost_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "ENV_ADD_POINT") == 0)
  {
    // custom_points
    printf("ERROR in param file \"%s\" on line %" PRId32
               ": Custom points management has been removed.\n",
        param_file_name_, cur_line);
    exit(EXIT_FAILURE);
  }
  else if ((strcmp(line->words[0], "ENV_ADD_GAUSSIAN") == 0) ||
      (strcmp(line->words[0], "ENV_GAUSSIAN") == 0))
  {
    #ifdef __REGUL
      // le premier chiffre est l'indice d'environment en convention humaine ( le premier a 1)
      // On vérifie que cet indice n'est pas trop élevé ni négatif pour éviter les crash
      if ( atoi(line->words[1]) - 1 < (int)_env_gaussians_list.size() && atoi(line->words[1]) > 0)
      {
        (_env_gaussians_list.at( atoi(line->words[1]) - 1)).push_back
        ( Gaussian(  atof( line->words[2] ), atof( line->words[3] ), atof( line->words[4] ) ) );
      }
      else
      {
        printf( " ERROR in param file \"%s\" on line %" PRId32 " : There is only %ld environment.\n",
         param_file_name_, cur_line, _env_gaussians_list.size() );
        exit( EXIT_FAILURE );
      }

    #else
          std_env_gaussians.push_back(
        Gaussian(atof(line->words[1]), atof(line->words[2]), atof(line->words[3])));
    #endif
  }
  else if (strcmp(line->words[0], "ENV_SAMPLING") == 0)
  {
    env_sampling_ = atoi(line->words[1]);
  }
  else if (strcmp(line->words[0], "ENV_VARIATION") == 0)
  {
    static bool env_var_already_set = false;
    if (env_var_already_set)
    {
      printf("ERROR in param file \"%s\" on line %" PRId32 " : "
                 "duplicate entry for %s.\n",
             param_file_name_, cur_line, line->words[0]);
      exit(EXIT_FAILURE);
    }
    env_var_already_set = true;

    if (strcmp(line->words[1], "none") == 0)
    {
      if (line->nb_words != 2) {
        printf("ERROR in param file \"%s\" on line %" PRId32
                   ": wrong number of parameters.\n",
               param_file_name_, cur_line);
        printf("usage: %s %s\n", line->words[0], line->words[1]);
        exit(EXIT_FAILURE);
      }
      env_var_method_ = NO_VAR;
    }
    else if (strcmp(line->words[1], "autoregressive_mean_variation") == 0)
    {
      if (line->nb_words != 5) {
        printf("ERROR in param file \"%s\" on line %" PRId32
                   ": wrong number of parameters.\n",
               param_file_name_, cur_line);
        printf("usage: %s %s sigma tau prng_seed\n",
               line->words[0], line->words[1]);
        exit(EXIT_FAILURE);
      }
      env_var_method_ = AUTOREGRESSIVE_MEAN_VAR;
      env_var_sigma_ = atof(line->words[2]);
      env_var_tau_ = atol(line->words[3]);
      env_var_seed_ = atoi(line->words[4]);
    }
    else if (strcmp(line->words[1], "autoregressive_height_variation") == 0)
    {
      if (line->nb_words != 5) {
        printf("ERROR in param file \"%s\" on line %" PRId32
                   ": wrong number of parameters.\n",
               param_file_name_, cur_line);
        printf("usage: %s %s sigma tau prng_seed\n",
               line->words[0], line->words[1]);
        exit(EXIT_FAILURE);
      }
      env_var_method_ = AUTOREGRESSIVE_HEIGHT_VAR;
      env_var_sigma_ = atof(line->words[2]);
      env_var_tau_ = atol(line->words[3]);
      env_var_seed_ = atoi(line->words[4]);
    }
    else if (strcmp(line->words[1], "add_local_gaussians") == 0)
    {
      if (line->nb_words != 3) {
        printf("ERROR in param file \"%s\" on line %" PRId32
                   ": wrong number of parameters.\n",
               param_file_name_, cur_line);
        printf("usage: %s %s prng_seed\n",
               line->words[0], line->words[1]);
        exit(EXIT_FAILURE);
      }
      env_var_method_ = LOCAL_GAUSSIANS_VAR;
      env_var_seed_ = atoi(line->words[2]);
    }
      #ifdef __REGUL
    else if (strcmp(line->words[1], "switch_in_a_list") == 0)
    {
      if (line->nb_words != 3) {
        printf("ERROR in param file \"%s\" on line %" PRId32
                   ": wrong number of parameters.\n",
               param_file_name_, cur_line);
        printf("usage: %s %s probability to switch between different environments\n",
               line->words[0], line->words[1]);
        exit(EXIT_FAILURE);
      }
      env_var_method_ = SWITCH_IN_A_LIST;
      _env_switch_probability = atof(line->words[2]);
    } else if (strcmp(line->words[1], "one_after_another") == 0)
    {
      if (line->nb_words != 2) {
        printf("ERROR in param file \"%s\" on line %" PRId32
                   ": wrong number of parameters.\n",
               param_file_name_, cur_line);
        printf("usage: %s %s probability to switch between different environments\n",
               line->words[0], line->words[1]);
        exit(EXIT_FAILURE);
      }
      env_var_method_ = ONE_AFTER_ANOTHER;
    }
    #endif
    else
    {
      printf("ERROR in param file \"%s\" on line %" PRId32
                 " : unknown phenotypic target variation method.\n",
             param_file_name_, cur_line);
      exit(EXIT_FAILURE);
    }
  }
  else if (strcmp(line->words[0], "ENV_NOISE") == 0)
  {
    static bool env_noise_already_set = false;
    if (env_noise_already_set)
    {
      printf("ERROR in param file \"%s\" on line %" PRId32
                 " : duplicate entry for %s.\n",
             param_file_name_, cur_line, line->words[0]);
      exit(EXIT_FAILURE);
    }
    env_noise_already_set = true;

    if (strcmp(line->words[1], "none") == 0)
    {
      assert(line->nb_words == 2);
      env_noise_method_ = NO_NOISE;
    }
    else if (strcmp(line->words[1], "FRACTAL") == 0)
    {
      assert(line->nb_words == 6);
      env_noise_method_ = FRACTAL;
      env_noise_sampling_log_ = atoi(line->words[2]);
      env_noise_sigma_ = atof(line->words[3]);
      env_noise_alpha_ = atof(line->words[4]);
      env_noise_prob_ = atof(line->words[5]);
      env_noise_seed_ = atoi(line->words[6]);
    }
    else
    {
      printf("ERROR in param file \"%s\" on line %" PRId32
                 " : unknown phenotypic target noise method.\n",
             param_file_name_,
             cur_line);
      exit(EXIT_FAILURE);
    }
  }
  else if (strcmp(line->words[0], "SECRETION_FITNESS_CONTRIB") == 0)
  {
    printf("ERROR in param file \"%s\" on line %" PRId32
               ": %s is no longer a valid option, "
               "use SECRETION_CONTRIB_TO_FITNESS instead.\n",
           param_file_name_, cur_line, line->words[0]);
    exit(EXIT_FAILURE);
  }
  else if (strcmp(line->words[0], "SECRETION_CONTRIB_TO_FITNESS") == 0)
  {
    secretion_contrib_to_fitness_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "SECRETION_DIFFUSION_PROP") == 0)
  {
    secretion_diffusion_prop_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "SECRETION_DEGRADATION_PROP") == 0)
  {
    secretion_degradation_prop_ = atof(line->words[1]);
    if (secretion_degradation_prop_ > 1 || secretion_degradation_prop_ < 0)
    {
      printf("ERROR in param file \"%s\" on line %" PRId32
                 ": degradation must be in (0,1).\n",
             param_file_name_, cur_line);
      exit(EXIT_FAILURE);
    }
  }
  else if (strcmp(line->words[0], "SECRETION_INITIAL") == 0)
  {
    secretion_init_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "SECRETION_COST") == 0)
  {
    secretion_cost_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "ALLOW_PLASMIDS") == 0)
  {
    if (strncmp(line->words[1], "true", 4) == 0)
    {
      allow_plasmids_ = true;
    }
    else if (strncmp(line->words[1], "false", 5) == 0)
    {
      allow_plasmids_ = false;
    }
    else
    {
      printf("ERROR in param file \"%s\" on line %" PRId32
                 ": unknown allow_plasmids option (use true/false).\n",
             param_file_name_, cur_line);
      exit(EXIT_FAILURE);
    }
  }
  else if (strcmp(line->words[0], "PLASMID_INITIAL_LENGTH") == 0)
  {
    plasmid_initial_length_ = atoi(line->words[1]);
  }
  else if (strcmp(line->words[0], "PLASMID_INITIAL_GENE") == 0)
  {
    plasmid_initial_gene_ = atoi(line->words[1]);
  }
  else if (strcmp(line->words[0], "PLASMID_MINIMAL_LENGTH") == 0)
  {
    plasmid_minimal_length_ = atoi(line->words[1]);
  }
  else if (strcmp(line->words[0], "PLASMID_MAXIMAL_LENGTH") == 0)
  {
    plasmid_maximal_length_ = atoi(line->words[1]);
  }
  else if (strcmp(line->words[0], "CHROMOSOME_MINIMAL_LENGTH") == 0)
  {
    chromosome_minimal_length_ = atoi(line->words[1]);
  }
  else if (strcmp(line->words[0], "CHROMOSOME_MAXIMAL_LENGTH") == 0)
  {
    chromosome_maximal_length_ = atoi(line->words[1]);
  }
  else if (strcmp(line->words[0], "PROB_HORIZONTAL_TRANS") == 0)
  {
    printf("ERROR in param file \"%s\" on line %" PRId32
               ": %s is no longer a valid option, "
               "did you mean PROB_PLASMID_HT ?.\n",
           param_file_name_, cur_line, line->words[0]);
    exit(EXIT_FAILURE);
  }
  else if (strcmp(line->words[0], "PROB_PLASMID_HT") == 0)
  {
    prob_plasmid_HT_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "TUNE_DONOR_ABILITY") == 0)
  {
    tune_donor_ability_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "TUNE_RECIPIENT_ABILITY") == 0)
  {
    tune_recipient_ability_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "DONOR_COST") == 0)
  {
    donor_cost_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "RECIPIENT_COST") == 0)
  {
    recipient_cost_ = atof(line->words[1]);
  }
  else if (strcmp(line->words[0], "COMPUTE_PHEN_CONTRIB_BY_GU") == 0)
  {
    if (strncmp(line->words[1], "true", 4) == 0)
    {
      compute_phen_contrib_by_GU_ = true;
    }
    else if (strncmp(line->words[1], "false", 5) == 0)
    {
      compute_phen_contrib_by_GU_ = false;
    }
    else
    {
      printf("ERROR in param file \"%s\" on line %" PRId32 " : unknown compute_phen_contrib_by_GU option (use true/false).\n",
             param_file_name_, cur_line);
      exit(EXIT_FAILURE);
    }
  }
  else if (strcmp(line->words[0], "LOG") == 0)
  {
    printf("LOGGING ");
    for (int8_t i = 1 ; i < line->nb_words ; i++)
    {
      if (strcmp(line->words[i], "TRANSFER") == 0)
      {
	printf("TRANSFER ");
        logs_ |= LOG_TRANSFER;
      }
      else if (strcmp(line->words[i], "REAR") == 0)
      {
	printf("REAR ");
        logs_ |= LOG_REAR;
      }
      else if (strcmp(line->words[i], "BARRIER") == 0)
      {
	printf("BARRIER ");
        logs_ |= LOG_BARRIER;
      }
        /*else if (strcmp(line->words[i], "LOADS") == 0)
        {
          tmp_to_be_logged |= LOG_LOADS;
        }   */
      else
      {
        printf("ERROR in param file \"%s\" on line %" PRId32 " : unknown log option %s.\n",
               param_file_name_, cur_line, line->words[1]);
        exit(EXIT_FAILURE);
      }
    }
    printf("(%d)\n",logs_);
  }
  else if (strcmp(line->words[0], "FUZZY_FLAVOR") == 0)
  {
    _fuzzy_flavor = atoi(line->words[1]);
  }  else if (strcmp(line->words[0], "SIMD_METADATA_FLAVOR") == 0)
  {
    if (strncmp(line->words[1], "stdmap", 6) == 0)
    {
      simd_metadata_flavor_ = STD_MAP;
      printf("Not fully working ATM, choose list\n"); exit(-1);
    }
    else if (strncmp(line->words[1], "dyntab", 6) == 0)
    {
      simd_metadata_flavor_ = DYN_TAB;
      printf("Not fully working ATM, choose list\n"); exit(-1);
    }
    else if (strncmp(line->words[1], "list", 6) == 0)
    {
      simd_metadata_flavor_ = STD_LIST;
    }
  }
  #ifdef BASE_4
else if(strcmp(line->words[0], "AMINO_ACID") == 0)
  {
    if(!mwh_bases_redefined_)
    {
      for(auto i = 0; i < NB_AMINO_ACIDS; i++)
      {
        aa_base_m_[i] = aa_base_w_[i] = aa_base_h_[i] = (int8_t) -1;
      }

      mwh_bases_redefined_ = true;
    }

    if(line->nb_words < 3)
    {
      printf("ERROR in param file \"%s\" on line %" PRId32 " : amino-acid base definition should at least comprise one amino-acid and one base digit\n",
             param_file_name_, cur_line);
      exit(EXIT_FAILURE);
    }

    AminoAcid amino_acid;
    if(str_to_aminoacid_.find(std::string(line->words[1])) != str_to_aminoacid_.end())
    {
      amino_acid = str_to_aminoacid_.at(std::string(line->words[1]));
    }
    else
    {
      printf("ERROR in param file \"%s\" on line %" PRId32 " : unrecognized amino-acid string : %s\n",
             param_file_name_, cur_line, line->words[1]);
      exit(EXIT_FAILURE);
    }

    for(auto i = 2; i < line->nb_words; i++)
    {
      if(strlen(line->words[i]) < 2)
      {
        printf("ERROR in param file \"%s\" on line %" PRId32 " : malformed base digit : %s\n",
               param_file_name_, cur_line, line->words[i]);
        exit(EXIT_FAILURE);
      }

      int8_t digit;

      try {
        digit = (int8_t) std::stoi(std::string(
            &(line->words[i][1])
        ));
      } catch (std::exception const &e) {
        printf("ERROR in param file \"%s\" on line %" PRId32 " : malformed base digit (invalid value) : %s\n",
               param_file_name_, cur_line, line->words[i]);
        exit(EXIT_FAILURE);
      }

      if(digit < 0)
      {
        printf("ERROR in param file \"%s\" on line %" PRId32 " : negative values aren't supported for base digits : %s\n",
               param_file_name_, cur_line, line->words[i]);
        exit(EXIT_FAILURE);
      }

      switch(line->words[i][0])
      {
      case 'M':
        aa_base_m_[amino_acid] = digit;
        break;
      case 'W':
        aa_base_w_[amino_acid] = digit;
        break;
      case 'H':
        aa_base_h_[amino_acid] = digit;
        break;
      }
    }
  }
  #endif

#ifdef __REGUL
    else if (strcmp(line->words[0], "HILL_SHAPE_N") == 0)
    {
      _hill_shape_n = atof(line->words[1]);
    }
    else if (strcmp(line->words[0], "HILL_SHAPE_THETA") == 0)
    {
      _hill_shape_theta = atof(line->words[1]);
    }
    else if (strcmp(line->words[0], "DEGRADATION_RATE") == 0)
    {
      _degradation_rate = atof(line->words[1]);
    }
    else if (strcmp(line->words[0], "NB_DEGRADATION_STEP") == 0)
    {
      _nb_degradation_step = atoi(line->words[1]);
    }
    else if (strcmp(line->words[0], "NB_INDIV_AGE") == 0)
    {
      _nb_indiv_age = atoi(line->words[1]);
    }
    else if (strcmp(line->words[0], "RANDOM_BINDING_MATRIX") == 0)
    {
        if (strncmp(line->words[1], "true", 4) == 0)
        {
        	_random_binding_matrix = true;
        }
        else if (strncmp(line->words[1], "false", 5) == 0)
        {
        	_random_binding_matrix = false;
        }
        else
        {
          printf("ERROR in param file \"%s\" on line %" PRId32 " : unknown more random_binding_matrix option (use true/false).\n",
                 param_file_name_, cur_line);
          exit(EXIT_FAILURE);
        }
    }
    else if (strcmp(line->words[0], "BINDING_ZEROS_PERCENTAGE") == 0)
    {
      _binding_zeros_percentage = atof(line->words[1]);
    }
    else if (strcmp(line->words[0], "INDIVIDUAL_EVALUATION_AGES") == 0)
    {
      _list_eval_step.clear();
      for (int i = 1; i < line->nb_words; i++) _list_eval_step.insert(atoi(line->words[i]));
    }
    else if (strcmp(line->words[0], "WITH_HEREDITY") == 0)
    {
      if (strncmp(line->words[1], "true", 4) == 0)
      {
        _with_heredity = true;
      }
      else if (strncmp(line->words[1], "false", 5) == 0)
      {
        _with_heredity = false;
      }
      else
      {
        printf("ERROR in param file \"%s\" on line %" PRId32 " : unknown with_heredity option (use true/false).\n",
               param_file_name_, cur_line);
        exit(EXIT_FAILURE);
      }
    }
    else if (strcmp(line->words[0], "PROTEIN_PRESENCE_LIMIT") == 0)
    {
      _protein_presence_limit = atof(line->words[1]);
    }
    else if (strcmp(line->words[0], "NB_ENVIRONMENTS") == 0)
    {
      int16_t nb_env = atoi( line->words[1] );

      if( nb_env < 1 )
      {
        printf( "ERROR in param file \"%s\" on line %" PRId32 " : you must have at least one environment\n", param_file_name_, cur_line );
        printf("you put %" PRId16 "\n", nb_env);
        exit( EXIT_FAILURE );
      }

      // Utile uniquement en cas de reprise sur backup
      // Je ne sais pas comment ça va se passer avec cette version ...
      if( _env_gaussians_list.size() > 0 )
      {
        _env_gaussians_list.clear();
      }


      if( _env_signals_list.size() > 0 )
      {
        _env_signals_list.clear();
      }


      for( int16_t i = 0; i < nb_env; i++)
      {
        _env_gaussians_list.push_back(std::list<Gaussian>());
        _env_signals_list.push_back(std::list<int16_t>());
      }
    }
    else if (strcmp(line->words[0], "CREATE_SIGNAL") == 0)
    {
      int signal_lenght = line->nb_words - 1;

      std::list<Codon*> codon_list;
      Codon* codon = NULL;
      for (int8_t i = 0; i < signal_lenght; i++)
      {
        if(strcmp(line->words[i+1], "h0")==0)
        {
          codon = new Codon(CODON_H0);
        }
        else if(strcmp(line->words[i+1], "h1")==0)
        {
          codon = new Codon(CODON_H1);
        }
        else if(strcmp(line->words[i+1], "w0")==0)
        {
          codon = new Codon(CODON_W0);
        }
        else if(strcmp(line->words[i+1], "w1")==0)
        {
          codon = new Codon(CODON_W1);
        }
        else if(strcmp(line->words[i+1], "m0")==0)
        {
          codon = new Codon(CODON_M0);
        }
        else if(strcmp(line->words[i+1], "m1")==0)
        {
          codon = new Codon(CODON_M1);
        }
        else
        {
          printf("Error this codon doesn't exist\n");
          exit( EXIT_FAILURE );
        }
        codon_list.push_back(codon);
      }
      _signals_models.push_back(new Protein_R(codon_list, 0.5, w_max_));

      for (auto cod : codon_list) delete cod;

      codon_list.clear();
    }
    else if (strcmp(line->words[0], "ENV_ADD_SIGNAL") == 0)
    {
      // le premier chiffre est l'indice d'environment en convention humaine ( le premier a 1)
      // On vérifie que cet indice n'est pas trop élevé ni négatif pour éviter les crash
      if ( atoi(line->words[1]) - 1 < (int)_env_signals_list.size() && atoi(line->words[1]) > 0)
      {
        (_env_signals_list.at( atoi(line->words[1]) - 1)).push_back(atoi(line->words[2]) - 1);
      }
      else
      {
        printf( " ERROR in param file \"%s\" on line %" PRId32 " : There are only %ld environment.\n",
         param_file_name_, cur_line, _env_gaussians_list.size() );
        exit( EXIT_FAILURE );
      }
    }
  #endif

  else
  {
    printf("ERROR in param file \"%s\" on line %" PRId32 " : undefined key word \"%s\"\n", param_file_name_, cur_line, line->words[0]);
    exit(EXIT_FAILURE);
  }
}

void ParamLoader::read_file()
{
  // The rewind is only necessary when using multiple param files
  rewind(param_file_);

  int32_t cur_line = 0;
  ParameterLine* parameter_line;

  // TODO : write parameter_line = new ParameterLine(param_file_) => ParameterLine::ParameterLine(char*)
  while ((parameter_line = line(&cur_line)) != NULL)
  {
    interpret_line(parameter_line, cur_line);
    delete parameter_line;
  }
}

void ParamLoader::CheckConsistency() {
  if (allow_plasmids_) {
    if (plasmid_initial_gene_ != 1) { // the plasmid will be copied from the chromosome
      if (plasmid_initial_length_ != -1) {
        printf(
            "WARNING: PLASMID_INITIAL_LENGTH is not taken into account because PLASMID_INITIAL_GENE is set to 0 (copy from chromosome)\n");
        plasmid_initial_length_ = chromosome_initial_length_;
      }
    }
    else if (compute_phen_contrib_by_GU_ == false) {
      printf("ERROR: when using PLASMID_INITIAL_GENE==1, the paramater COMPUTE_PHEN_CONTRIB_BY_GU should be set to true.\n");
      exit(EXIT_FAILURE);
    }

    if (plasmid_maximal_length_ == -1)
      plasmid_maximal_length_ = max_genome_length_;
    if (plasmid_minimal_length_ == -1)
      plasmid_minimal_length_ = min_genome_length_;
    if(plasmid_minimal_length_ > plasmid_initial_length_) {
      printf("ERROR: PLASMID_INITIAL_LENGTH is lower than PLASMID_MINIMAL_LENGTH\n");
      exit(EXIT_FAILURE);
    }
    if (plasmid_maximal_length_ < plasmid_initial_length_) {
      printf("ERROR: PLASMID_INITIAL_LENGTH is higher than PLASMID_MAXIMAL_LENGTH\n");
      exit(EXIT_FAILURE);
    }
  }
  if (chromosome_maximal_length_ == -1)
    chromosome_maximal_length_ = max_genome_length_;
  if (chromosome_minimal_length_ == -1)
    chromosome_minimal_length_ = min_genome_length_;
  if (chromosome_minimal_length_ > chromosome_initial_length_) {
    printf("ERROR: CHROMOSOME_INITIAL_LENGTH is lower than CHROMOSOME_MINIMAL_LENGTH\n");
    exit(EXIT_FAILURE);
  }
  if (chromosome_maximal_length_ < chromosome_initial_length_) {
    printf("ERROR: CHROMOSOME_INITIAL_LENGTH is higher than PLASMID_MAXIMAL_LENGTH\n");
    exit(EXIT_FAILURE);
  }
  // Check that the population fits in the spatial structure
  if (init_pop_size_ != grid_width_ * grid_height_)
  {
    printf("ERROR: the number of individuals (%" PRId32
               ") does not match the size of the grid  (%" PRId16
               " * %" PRId16 ")\n",
           init_pop_size_,
           grid_width_,
           grid_height_);
    exit(EXIT_FAILURE);
  }
}

FuzzyFactory* FuzzyFactory::fuzzyFactory = NULL;

void ParamLoader::load(ExpManager * exp_m, bool verbose,
                       char* chromosome, int32_t lchromosome,
                       char* plasmid, int32_t lplasmid
                       #ifdef HAVE_MPI
                      , int32_t nb_rank
                      #endif
                      ) {
  // Check consistency of min, max and initial length of chromosome and plasmid
  // Default for by GU minimal or maximal size is -1.
  // If equal to -1, maximal sizes of each GU will be replaced by total maximal size for the whole genome
  CheckConsistency();

#ifdef HAVE_MPI
  global_pop_size_ = init_pop_size_;
  global_grid_width_ = grid_width_;
  global_grid_height_ = grid_height_;
  rank_width_ = 1;
  rank_height_ = 1;

  int min_dist = grid_width_;
  printf("Init width rank %d (global_pop_size %d)\n",nb_rank,global_pop_size_%1);

  for (int32_t i = 2; i <= nb_rank; i++) {
    printf("Testing width rank %d (global_pop_size %d)\n",i,global_pop_size_%i);
    if (global_pop_size_ % i == 0) {
      int32_t seg_rank_width_ = i;
      int32_t seg_rank_height_ = nb_rank / i;

      if (seg_rank_height_ * seg_rank_width_ > nb_rank) {
        printf("Incorrect Nb Rank %d %d (%d)\n",seg_rank_width_,seg_rank_height_,nb_rank);
        continue;
      }

      int32_t seg_grid_width = global_grid_width_ / seg_rank_width_;
      int32_t seg_grid_height = global_grid_height_ / seg_rank_height_;

      if (seg_grid_width*seg_grid_height*nb_rank != init_pop_size_) {
        printf("Invalid population size %d (%d %d) Nb Rank %d Pop Size %d\n",seg_grid_width*seg_grid_height,seg_grid_width,seg_grid_height,nb_rank,init_pop_size_);
        continue;
      }

      printf("Evaluating %d x %d : %d (%d x %d) dist %d (min dist %d)\n",seg_rank_width_,seg_rank_height_,
            seg_grid_width*seg_grid_height,
            seg_grid_width,seg_grid_height,
            std::abs(seg_grid_width - seg_grid_height),min_dist);
      if (std::abs(seg_grid_width - seg_grid_height) < min_dist) {
        grid_width_ = seg_grid_width;
        grid_height_ = seg_grid_height;
        min_dist = std::abs(seg_grid_width - seg_grid_height);
        init_pop_size_ = grid_width_*grid_height_;
        rank_width_ = seg_rank_width_;
        rank_height_ = seg_rank_height_;
      }
    }
  }

  printf("Global population is %d (%d x %d) : Patch population is %d (%d x %d)\n",
      global_pop_size_,global_grid_width_,global_grid_height_,
      init_pop_size_,grid_width_,grid_height_);
#endif
  // Initialize prng_
  // This one will be used to create the initial genome(s) and to generate seeds for other prng
  prng_ = std::make_shared<JumpingMT>(seed_);

  // Initialize mut_prng, stoch_prng, world_prng :
  // if mut_seed (respectively stoch_seed) not given in param.in, choose it at random
  if (mut_seed_ == 0) {
    mut_seed_ = prng_->random(1000000);
  }
  if (stoch_seed_ == 0) {
    stoch_seed_ = prng_->random(1000000);
  }
  auto mut_prng   = std::make_shared<JumpingMT>(mut_seed_);
  auto stoch_prng = std::make_shared<JumpingMT>(stoch_seed_);
  auto world_prng = std::make_shared<JumpingMT>(prng_->random(1000000));

  // Create aliases
  ExpSetup* exp_s = exp_m->exp_s();
  Selection* sel = exp_m->sel();
  OutputManager* output_m = exp_m->output_m();
  output_m->InitStats();


  // 1) ------------------------------------- Initialize the experimental setup


  // ---------------------------------------------------------------- Selection
  sel->set_selection_scheme(selection_scheme_);
  sel->set_selection_pressure(selection_pressure_);

  sel->set_selection_scope(selection_scope_);
  sel->set_selection_scope_x(selection_scope_x_);
  sel->set_selection_scope_y(selection_scope_y_);

  sel->set_fitness_function(fitness_function_);
  sel->set_fitness_function_scope_x(fitness_function_x_);
  sel->set_fitness_function_scope_y(fitness_function_y_);
  // ----------------------------------------------------------------- Transfer
  exp_s->set_with_HT(with_HT_);
  exp_s->set_repl_HT_with_close_points(repl_HT_with_close_points_);
  exp_s->set_HT_ins_rate(HT_ins_rate_);
  exp_s->set_HT_repl_rate(HT_repl_rate_);
  exp_s->set_repl_HT_detach_rate(repl_HT_detach_rate_);

  // ----------------------------------------------------------------- Plasmids
  exp_s->set_with_plasmids(allow_plasmids_);
  exp_s->set_prob_plasmid_HT(prob_plasmid_HT_);
  exp_s->set_tune_donor_ability(tune_donor_ability_);
  exp_s->set_tune_recipient_ability(tune_recipient_ability_);
  exp_s->set_donor_cost(donor_cost_);
  exp_s->set_recipient_cost(recipient_cost_);
  exp_s->set_swap_GUs(swap_GUs_);
  output_m->set_compute_phen_contrib_by_GU(compute_phen_contrib_by_GU_);

  // ---------------------------------------------------------------- Secretion
  exp_s->set_with_secretion(with_secretion_);
  exp_s->set_secretion_contrib_to_fitness(secretion_contrib_to_fitness_);
  exp_s->set_secretion_cost(secretion_cost_);

  exp_s->set_fuzzy_flavor(_fuzzy_flavor);
  exp_s->set_simd_metadata_flavor(simd_metadata_flavor_);


  //------------------------------------------------------------------ Parameter for SIMD
  exp_s->set_min_genome_length(min_genome_length_);
  exp_s->set_max_genome_length(max_genome_length_);

#ifdef BASE_4
 // ----------------------------------------------------------------- Terminators
  exp_s->set_terminator_polya_sequence_length(term_polya_seq_length_);

  // -------------------------------------------------- MWH bases configuration
  if(mwh_bases_redefined_) {
    exp_s->set_aa_base_m(aa_base_m_);
    exp_s->set_aa_base_w(aa_base_w_);
    exp_s->set_aa_base_h(aa_base_h_);
  }
#endif

#ifdef __REGUL
  if (env_var_method_ == ONE_AFTER_ANOTHER)
    exp_s->set_with_heredity(false);
  else
    exp_s->set_with_heredity(_with_heredity);

  exp_s->set_degradation_rate(_degradation_rate);
  exp_s->set_nb_degradation_step(_nb_degradation_step);
  exp_s->set_nb_indiv_age(_nb_indiv_age);
  exp_s->set_list_eval_step(_list_eval_step);
  exp_s->set_protein_presence_limit(_protein_presence_limit);
  exp_s->set_hill_shape(pow( _hill_shape_theta, _hill_shape_n ));
  exp_s->set_hill_shape_n( _hill_shape_n );

  // printf("Random %d : Sparsity %lf\n",_random_binding_matrix,_binding_zeros_percentage);
  exp_s->init_binding_matrix(_random_binding_matrix,_binding_zeros_percentage,prng_);
#endif

  #ifdef HAVE_MPI
  //------------------------------------------------------ Global Grid characteristics
  exp_s->set_global_pop_size(global_pop_size_);
  exp_s->set_global_grid_width(global_grid_width_);
  exp_s->set_global_grid_height(global_grid_height_);

  exp_s->set_rank_width(rank_width_);
  exp_s->set_rank_height(rank_height_);
  #endif

  if (FuzzyFactory::fuzzyFactory == NULL)
    FuzzyFactory::fuzzyFactory = new FuzzyFactory();

  // 2) --------------------------------------------- Create and init a Habitat
  #ifndef __REGUL
  Habitat habitat;
  #else
  Habitat_R habitat;
  #endif

  // Shorthand for phenotypic target handler
  #ifndef __REGUL
  PhenotypicTargetHandler& phenotypic_target_handler =
      habitat.phenotypic_target_handler_nonconst();
  #else
  PhenotypicTargetHandler_R& phenotypic_target_handler =
      habitat.phenotypic_target_handler_nonconst();
  #endif

  // Move the gaussian list from the parameters to the phen target handler
  #ifndef __REGUL
  phenotypic_target_handler.set_gaussians(std_env_gaussians);
  #else
  phenotypic_target_handler.set_gaussians(_env_gaussians_list);
  phenotypic_target_handler.set_signals_models(_signals_models);
  phenotypic_target_handler.set_signals(_env_signals_list);
  #endif

  // Copy the sampling
  phenotypic_target_handler.set_sampling(env_sampling_);

  // Set phenotypic target segmentation

  if((env_axis_features_ != NULL) && (env_axis_segment_boundaries_ != NULL)) {
    // if param.in contained a line starting with ENV_AXIS_FEATURES,
    // we use the values indicated on this line
    phenotypic_target_handler.set_segmentation(env_axis_nb_segments_,
                                               env_axis_segment_boundaries_,
                                               env_axis_features_,
                                               env_axis_separate_segments_);
  }
  // else we leave the segmentation as it is by default
  // (one "metabolic" segment from X_MIN to X_MAX)


  // Set phenotypic target variation
  if (env_var_method_ != NO_VAR)
  {
    phenotypic_target_handler.set_var_method(env_var_method_);
    phenotypic_target_handler.set_var_prng(std::make_shared<JumpingMT>(env_var_seed_));
    phenotypic_target_handler.set_var_sigma_tau(env_var_sigma_, env_var_tau_);
#ifdef __REGUL
    phenotypic_target_handler.set_switch_probability(_env_switch_probability);
#endif
  }

  // Set phenotypic target noise
  if (env_noise_method_ != NO_NOISE)
  {
    phenotypic_target_handler.set_noise_method(env_noise_method_);
    phenotypic_target_handler.set_noise_sampling_log(env_noise_sampling_log_);
    phenotypic_target_handler.set_noise_prng(std::make_shared<JumpingMT>(env_noise_seed_));
    phenotypic_target_handler.set_noise_alpha(env_noise_alpha_);
    phenotypic_target_handler.set_noise_sigma(env_noise_sigma_);
    phenotypic_target_handler.set_noise_prob(env_noise_prob_);
  }

  // Build the phenotypic target
  #ifndef __REGUL
  phenotypic_target_handler.BuildPhenotypicTarget();
  #else
  printf("Init phenotypic target with %d\n",_nb_indiv_age);
  phenotypic_target_handler.InitPhenotypicTargetsAndModels( _nb_indiv_age );
  #endif

  if (verbose) {
    #ifndef __REGUL
    printf("Entire geometric area of the phenotypic target : %f\n",
           phenotypic_target_handler.get_geometric_area());
    #else
    phenotypic_target_handler.print_geometric_areas();
    #endif
  }


  // 3) --------------------------------------------- Create the new population
  list<Individual *> indivs;
  // Generate a model ae_mut_param object
  auto param_mut = std::make_shared<MutationParams>();
  param_mut->set_point_mutation_rate(point_mutation_rate_);
  param_mut->set_small_insertion_rate(small_insertion_rate_);
  param_mut->set_small_deletion_rate(small_deletion_rate_);
  param_mut->set_max_indel_size(max_indel_size_);
  param_mut->set_with_4pts_trans(with_4pts_trans_);
  param_mut->set_with_alignments(with_alignments_);
  param_mut->set_with_HT(with_HT_);
  param_mut->set_repl_HT_with_close_points(repl_HT_with_close_points_);
  param_mut->set_HT_ins_rate(HT_ins_rate_);
  param_mut->set_HT_repl_rate(HT_repl_rate_);
  param_mut->set_repl_HT_detach_rate(repl_HT_detach_rate_);
  param_mut->set_duplication_rate(duplication_rate_);
  param_mut->set_deletion_rate(deletion_rate_);
  param_mut->set_translocation_rate(translocation_rate_);
  param_mut->set_inversion_rate(inversion_rate_);
  param_mut->set_neighbourhood_rate(neighbourhood_rate_);
  param_mut->set_duplication_proportion(duplication_proportion_);
  param_mut->set_deletion_proportion(deletion_proportion_);
  param_mut->set_translocation_proportion(translocation_proportion_);
  param_mut->set_inversion_proportion(inversion_proportion_);

  exp_s->set_mutation_parameters(param_mut);

  Individual * indiv = NULL;
  int32_t id_new_indiv = 0;

  if (chromosome != NULL)
  {
    printf("Option -c is used: chromosome will be loaded from a text file\n");
    #ifndef __REGUL
    Individual * indiv = new Individual(exp_m,
                                             mut_prng,
                                             stoch_prng,
                                             param_mut,
                                             w_max_,
                                             min_genome_length_,
                                             max_genome_length_,
                                             allow_plasmids_,
                                             id_new_indiv++,
                                             strain_name_,
                                             0);
    #else
    Individual_R * indiv = new Individual_R(exp_m,
                                         mut_prng,
                                         stoch_prng,
                                         param_mut,
                                         w_max_,
                                         min_genome_length_,
                                         max_genome_length_,
                                         allow_plasmids_,
                                         id_new_indiv++,
                                         strain_name_,
                                         0);

    #endif

    indiv->add_GU(chromosome, lchromosome);
    indiv->genetic_unit_nonconst(0).set_min_gu_length(chromosome_minimal_length_);
    indiv->genetic_unit_nonconst(0).set_max_gu_length(chromosome_maximal_length_);

    if (plasmid != NULL)
    {
      printf("Option -p is used: plasmid will be loaded from a text file\n");
      if (! allow_plasmids_)
      {
        printf("ERROR: option -p requires ALLOW_PLASMIDS set to true\n");
        exit(EXIT_FAILURE);
      }
      indiv->add_GU(plasmid, lplasmid);
      indiv->genetic_unit_nonconst(1).set_min_gu_length(plasmid_minimal_length_);
      indiv->genetic_unit_nonconst(1).set_max_gu_length(plasmid_maximal_length_);
    }
    else if (allow_plasmids_)
    {
      printf("ERROR: if you use option -c and ALLOW_PLASMIDS is set to true, you must also use option -p. \n For now loading a genetic unit from text file and generating the other is not supported.\n");
      exit(EXIT_FAILURE);
    }

    indiv->set_with_stochasticity(with_stochasticity_);
    indiv->compute_statistical_data();
    indiv->EvaluateInContext(habitat);
    printf("Starting with a clonal population of individual with metabolic error %f and secretion error %f \n",indiv->dist_to_target_by_feature(METABOLISM),indiv->dist_to_target_by_feature(SECRETION));
    indivs.push_back(indiv);

    // Make the clones and add them to the list of individuals
    for (int32_t i = 1 ; i < init_pop_size_ ; i++)
    {
      #ifndef __REGUL
      Individual * clone = Individual::CreateClone(indiv, id_new_indiv++);
      #else
      Individual_R * clone = Individual_R::CreateClone(indiv, id_new_indiv++);
      #endif
      clone->EvaluateInContext(habitat);
      indivs.push_back(clone);
    }
  }
  else if (plasmid != NULL)
  {
    printf("ERROR: option -p can only be used in combination with option -c for now\n");
    exit(EXIT_FAILURE);
  }
  else if (init_method_ & ONE_GOOD_GENE)
  {
    if (init_method_ & CLONE)
    {
      // Create an individual with a "good" gene (in fact, make an indiv whose
      // fitness is better than that corresponding to a flat phenotype)
      // and set its id
      indiv = IndividualFactory::create_random_individual(
          exp_m,
          id_new_indiv++,
          param_mut,
          mut_prng,
          stoch_prng,
          habitat,
          w_max_,
          min_genome_length_,
          max_genome_length_,
          chromosome_initial_length_,
          allow_plasmids_,
          plasmid_initial_gene_,
          plasmid_initial_length_,
          strain_name_,
          prng_,
          true);
      indiv->genetic_unit_nonconst(0).set_min_gu_length(chromosome_minimal_length_);
      indiv->genetic_unit_nonconst(0).set_max_gu_length(chromosome_maximal_length_);

      if (allow_plasmids_)
      {
        indiv->genetic_unit_nonconst(1).set_min_gu_length(plasmid_minimal_length_);
        indiv->genetic_unit_nonconst(1).set_max_gu_length(plasmid_maximal_length_);
      }

      indiv->set_with_stochasticity(with_stochasticity_);

      // Add it to the list
      indivs.push_back(indiv);

      // Make the clones and add them to the list of individuals
      for (int32_t i = 1 ; i < init_pop_size_ ; i++)
      {
        // Add new clone to the list
        #ifndef __REGUL
        Individual * clone = Individual::CreateClone(indiv, id_new_indiv++);
        #else
        Individual_R * clone = Individual_R::CreateClone(dynamic_cast<Individual_R*>(indiv), id_new_indiv++);
        #endif
        clone->EvaluateInContext(habitat);
        indivs.push_back(clone);
      }
    }
    else // if (! CLONE)
    {
      for (int32_t i = 0 ; i < init_pop_size_ ; i++)
      {
        // Create an individual and set its id
        indiv = IndividualFactory::create_random_individual(
            exp_m,
            id_new_indiv++,
            param_mut,
            mut_prng,
            stoch_prng,
            habitat,
            w_max_,
            min_genome_length_,
            max_genome_length_,
            chromosome_initial_length_,
            allow_plasmids_,
            plasmid_initial_gene_,
            plasmid_initial_length_,
            strain_name_,
            prng_,
            true);
        indiv->genetic_unit_nonconst(0).set_min_gu_length(chromosome_minimal_length_);
        indiv->genetic_unit_nonconst(0).set_max_gu_length(chromosome_maximal_length_);
        if (allow_plasmids_)
        {
          indiv->genetic_unit_nonconst(1).set_min_gu_length(plasmid_minimal_length_);
          indiv->genetic_unit_nonconst(1).set_max_gu_length(plasmid_maximal_length_);
        }

        // Add it to the list
        indivs.push_back(indiv);
      }
    }
  }
  else // if (! ONE_GOOD_GENE)
  {
    if (init_method_ & CLONE)
    {
      // Create a random individual and set its id
      indiv = IndividualFactory::create_random_individual(
          exp_m,
          id_new_indiv++,
          param_mut,
          mut_prng,
          stoch_prng,
          habitat,
          w_max_,
          min_genome_length_,
          max_genome_length_,
          chromosome_initial_length_,
          allow_plasmids_,
          plasmid_initial_gene_,
          plasmid_initial_length_,
          strain_name_,
          prng_,
          false);
      indiv->genetic_unit_nonconst(0).set_min_gu_length(chromosome_minimal_length_);
      indiv->genetic_unit_nonconst(0).set_max_gu_length(chromosome_maximal_length_);
      if (allow_plasmids_)
      {
        indiv->genetic_unit_nonconst(1).set_min_gu_length(plasmid_minimal_length_);
        indiv->genetic_unit_nonconst(1).set_max_gu_length(plasmid_maximal_length_);
      }

      indiv->clear_everything_except_dna_and_promoters();
      indiv->EvaluateInContext(habitat);

      // Add it to the list
      indivs.push_back(indiv);

      // Make the clones and add them to the list of individuals
      for (int32_t i = 1 ; i < init_pop_size_ ; i++)
      {
        // Add clone to the list
        #ifndef __REGUL
        Individual * clone = Individual::CreateClone(indiv, id_new_indiv++);
        #else
        Individual_R * clone = Individual_R::CreateClone(dynamic_cast<Individual_R*>(indiv), id_new_indiv++);
        #endif
        clone->EvaluateInContext(habitat);
        indivs.push_back(clone);
      }
    }
    else // if (! CLONE)
    {
      for (int32_t i = 0 ; i < init_pop_size_ ; i++)
      {
        // Create a random individual and set its id
        indiv = IndividualFactory::create_random_individual(
            exp_m,
            id_new_indiv++,
            param_mut,
            mut_prng,
            stoch_prng,
            habitat,
            w_max_,
            min_genome_length_,
            max_genome_length_,
            chromosome_initial_length_,
            allow_plasmids_,
            plasmid_initial_gene_,
            plasmid_initial_length_,
            strain_name_,
            prng_,
            false);
        indiv->genetic_unit_nonconst(0).set_min_gu_length(chromosome_minimal_length_);
        indiv->genetic_unit_nonconst(0).set_max_gu_length(chromosome_maximal_length_);
        if (allow_plasmids_)
        {
          indiv->genetic_unit_nonconst(1).set_min_gu_length(plasmid_minimal_length_);
          indiv->genetic_unit_nonconst(1).set_max_gu_length(plasmid_maximal_length_);
        }

        // Add it to the list
        indivs.push_back(indiv);
      }
    }
  }

  // -------------------------------------------------------- Spatial structure  
  #ifdef HAVE_MPI
  exp_m->mut_prng_seed_ = new int32_t[global_pop_size_];
  exp_m->stoch_prng_seed_ = new int32_t[global_pop_size_];

  for (int16_t x = 0; x < global_grid_width_; x++) {
    for (int16_t y = 0; y < global_grid_height_; y++) {
      exp_m->mut_prng_seed_[x*global_grid_height_+y]  = mut_prng->random(1000000);
      exp_m->stoch_prng_seed_[x*global_grid_height_+y] = stoch_prng->random(1000000);
    }
  }
  #endif

  exp_m->InitializeWorld(grid_width_, grid_height_,
                         world_prng,mut_prng,stoch_prng,
                         habitat,
                         true);
  World* world = exp_m->world();

#ifdef HAVE_MPI  
  exp_m->prng_seed_ = new int32_t[global_pop_size_];

  for (int16_t x = 0; x < global_grid_width_; x++) {
    for (int16_t y = 0; y < global_grid_height_; y++) {
      exp_m->prng_seed_[x*global_grid_height_+y] = prng_->random(1000000);
    }
  }
#endif

  for (int16_t x = 0; x < exp_m->grid_width(); x++) {
    for (int16_t y = 0; y < exp_m->grid_height(); y++) {
      #ifdef HAVE_MPI
      int32_t seed = exp_m->prng_seed_[x*exp_m->grid_height()+y];
      #else
      int32_t seed = prng_->random(1000000);
      #endif
#if __cplusplus == 201103L
      exp_m->world()->grid(x,y)->set_reprod_prng(make_unique<JumpingMT>(seed));
      exp_m->world()->grid(x,y)->set_reprod_prng_simd(make_unique<JumpingMT>(seed));
#else
      exp_m->world()->grid(x,y)->set_reprod_prng(std::make_unique<JumpingMT>(seed));
      exp_m->world()->grid(x,y)->set_reprod_prng_simd(std::make_unique<JumpingMT>(seed));
#endif
    }
  }

  world->set_secretion_degradation_prop(secretion_degradation_prop_);
  world->set_secretion_diffusion_prop(secretion_diffusion_prop_);
  world->set_is_well_mixed(well_mixed);
  world->set_partial_mix_nb_permutations(partial_mix_nb_permutations);

  // Set each individual's position on the grid
  int16_t x, y;
  int16_t x_max = exp_m->grid_width();
  int16_t y_max = exp_m->grid_height();

  for (const auto& indiv: indivs) {
    do {
      x = exp_m->world()->prng()->random(x_max);
      y = exp_m->world()->prng()->random(y_max);
    } while (world->indiv_at(x, y) != NULL);

    world->PlaceIndiv(indiv, x, y, true);
  }

  world->set_best(0, 0);



  // 4) ------------------------------------------ Set the recording parameters
  output_m->set_backup_step(backup_step_);
  output_m->set_big_backup_step(big_backup_step_);

  if (record_tree_)
  {
    output_m->init_tree(exp_m, tree_step_);
  }

  output_m->init_light_tree(record_light_tree_,exp_m, tree_step_);

  if (make_dumps_)
  {
    output_m->set_dump_step(dump_step_);
  }
  output_m->set_logs(logs_);
}



// =================================================================
//                           Protected Methods
// =================================================================
/*!
  \brief Format a line by parsing it and the words inside

  \param formated_line the resulted formated line
  \param line original line in char*
  \param line_is_interpretable boolean with about the possible intrepretation of the line
*/
void ParamLoader::format_line(ParameterLine * formated_line, char* line, bool* line_is_interpretable)
{
  int32_t i = 0;
  int32_t j;

  // Parse line
  while (line[i] != '\n' && line[i] != '\0' && line[i] != '\r')
  {
    j = 0;

    // Flush white spaces and tabs
    while (line[i] == ' ' || line[i] == kTabChar) i++;

    // Check comments
    if (line[i] == '#') break;

    // If we got this far, there is content in the line
    *line_is_interpretable = true;

    // Parse word
    while (line[i] != ' ' && line[i] != kTabChar && line[i] != '\n' &&
        line[i] != '\0' && line[i] != '\r') {
      formated_line->words[formated_line->nb_words][j++] = line[i++];
    }

    // Add '\0' at end of word if it's not empty (line ending with space or tab)
    if (j != 0)
    {
      formated_line->words[formated_line->nb_words++][j] = '\0';
    }
  }
}

/*!
  \brief Get a line in a file and format it

  \return line (pointer)

  \see format_line(ParameterLine* formated_line, char* line, bool* line_is_interpretable)
*/
ParameterLine *ParamLoader::line(int32_t* cur_line_ptr) // void
{
  char line[4096];
  ParameterLine * formated_line = new ParameterLine();

  bool found_interpretable_line = false; // Found line that is neither a comment nor empty

  while (!feof(param_file_) && !found_interpretable_line)
  {
    if (!fgets(line, 4096, param_file_))
    {
      delete formated_line;
      return NULL;
    }
    (*cur_line_ptr)++;
    format_line(formated_line, line, &found_interpretable_line);
  }

  if (found_interpretable_line)
  {
    return formated_line;
  }
  else
  {
    delete formated_line;
    return NULL;
  }
}

void ParamLoader::print_to_file(FILE* file)
{
  // ------------------------------------------------------------ Constraints
  fprintf(file, "\nConstraints ---------------------------------------------\n");
  fprintf(file, "min_genome_length :          %" PRId32 "\n", min_genome_length_);
  fprintf(file, "max_genome_length :          %" PRId32 "\n", max_genome_length_);
  fprintf(file, "W_MAX :                      %f\n",        w_max_);

  // --------------------------------------------------------- Mutation rates
  fprintf(file, "\nMutation rates ------------------------------------------\n");
  fprintf(file, "point_mutation_rate :        %e\n",  point_mutation_rate_);
  fprintf(file, "small_insertion_rate :       %e\n",  small_insertion_rate_);
  fprintf(file, "small_deletion_rate :        %e\n",  small_deletion_rate_);
  fprintf(file, "max_indel_size :             %" PRId16 "\n", max_indel_size_);

  // -------------------------------------------- Rearrangements and Transfer
  fprintf(file, "\nRearrangements and Transfer -----------------------------\n");
  fprintf(file, "with_4pts_trans :            %s\n",  with_4pts_trans_? "true" : "false");
  fprintf(file, "with_alignments :            %s\n",  with_alignments_? "true" : "false");
  fprintf(file, "with_HT :                    %s\n",  with_HT_? "true" : "false");
  fprintf(file, "repl_HT_with_close_points :  %s\n",  repl_HT_with_close_points_? "true" : "false");
  fprintf(file, "HT_ins_rate :                %e\n",  HT_ins_rate_);
  fprintf(file, "HT_repl_rate :               %e\n",  HT_repl_rate_);

  // ---------------------------------------------------- Rearrangement rates
  if (with_alignments_)
  {
    fprintf(file, "\nRearrangement rates (with alignements) ------------------\n");
    fprintf(file, "neighbourhood_rate :         %e\n",  neighbourhood_rate_);
    fprintf(file, "duplication_proportion :     %e\n",  duplication_proportion_);
    fprintf(file, "deletion_proportion :        %e\n",  deletion_proportion_);
    fprintf(file, "translocation_proportion :   %e\n",  translocation_proportion_);
    fprintf(file, "inversion_proportion :       %e\n",  inversion_proportion_);
  }
  else
  {
    fprintf(file, "\nRearrangement rates (without alignements) ----------------\n");
    fprintf(file, "duplication_rate :           %e\n",  duplication_rate_);
    fprintf(file, "deletion_rate :              %e\n",  deletion_rate_);
    fprintf(file, "translocation_rate :         %e\n",  translocation_rate_);
    fprintf(file, "inversion_rate :             %e\n",  inversion_rate_);
  }

  // ------------------------------------------------------------ Alignements
  fprintf(file, "\nAlignements ---------------------------------------------\n");
  fprintf(file, "align_fun_shape :            %" PRId16 "\n", (int16_t) align_fun_shape_);
  fprintf(file, "align_sigm_lambda :          %f\n",        align_sigm_lambda_);
  fprintf(file, "align_sigm_mean :            %" PRId16 "\n", align_sigm_mean_);
  fprintf(file, "align_lin_min :              %" PRId16 "\n", align_lin_min_);
  fprintf(file, "align_lin_max :              %" PRId16 "\n", align_lin_max_);
  fprintf(file, "align_max_shift :            %" PRId16 "\n", align_max_shift_);
  fprintf(file, "align_w_zone_h_len :         %" PRId16 "\n", align_w_zone_h_len_);
  fprintf(file, "align_match_bonus :          %" PRId16 "\n", align_match_bonus_);
  fprintf(file, "align_mismatch_cost :        %" PRId16 "\n", align_mismatch_cost_);

  // -------------------------------------------------------------- Selection
  fprintf(file, "\nSelection -----------------------------------------------\n");
  switch (selection_scheme_)
  {
    case RANK_LINEAR :
    {
      fprintf(file, "selection_scheme :           RANK_LINEAR\n");
      break;
    }
    case RANK_EXPONENTIAL :
    {
      fprintf(file, "selection_scheme :           RANK_EXPONENTIAL\n");
      break;
    }
    case FITNESS_PROPORTIONATE :
    {
      fprintf(file, "selection_scheme :           FITNESS_PROPORTIONATE\n");
      break;
    }
    case FITTEST :
    {
      fprintf(file, "selection_scheme :           FITTEST\n");
      break;
    }
    default :
    {
      fprintf(file, "selection_scheme :           UNKNOWN\n");
      break;
    }
  }
  fprintf(file, "selection_pressure :         %e\n",  selection_pressure_);

  switch (selection_scope_)
  {
    case SCOPE_GLOBAL :
    {
      fprintf(file, "selection_scope :           GLOBAL\n");
      break;
    }
    case SCOPE_LOCAL :
    {
      fprintf(file, "selection_scope :          LOCAL %d %d\n",selection_scope_x_,selection_scope_y_);
      break;
    }
    default :
    {
      fprintf(file, "selection_scope :           UNKNOWN\n");
      break;
    }
  }

  // -------------------------------------------------------------- Secretion
  fprintf(file, "\nSecretion -----------------------------------------------\n");
  fprintf(file, "with_secretion :                %s\n", with_secretion_? "true" : "false");
  fprintf(file, "secretion_contrib_to_fitness :  %e\n", secretion_contrib_to_fitness_);
  fprintf(file, "secretion_diffusion_prop :      %e\n", secretion_diffusion_prop_);
  fprintf(file, "secretion_degradation_prop :    %e\n", secretion_degradation_prop_);
  fprintf(file, "secretion_cost :                %e\n", secretion_cost_);

  // --------------------------------------------------------------- Plasmids
  fprintf(file, "\nPlasmids ------------------------------------------------\n");
  fprintf(file, "allow_plasmids :             %s\n", allow_plasmids_? "true" : "false");
  fprintf(file, "plasmid_minimal_length :     %" PRId32 "\n", plasmid_minimal_length_);
  fprintf(file, "plasmid_maximal_length :     %" PRId32 "\n", plasmid_maximal_length_);
  fprintf(file, "chromosome_minimal_length :  %" PRId32 "\n", chromosome_minimal_length_);
  fprintf(file, "chromosome_maximal_length :  %" PRId32 "\n", chromosome_maximal_length_);
  fprintf(file, "prob_plasmid_HT :            %e\n", prob_plasmid_HT_);
  fprintf(file, "tune_donor_ability :         %e\n", tune_donor_ability_);
  fprintf(file, "tune_recipient_ability :     %e\n", tune_recipient_ability_);
  fprintf(file, "donor_cost :                 %e\n", donor_cost_);
  fprintf(file, "recipient_cost :             %e\n", recipient_cost_);
  fprintf(file, "compute_phen_contrib_by_GU : %s\n", compute_phen_contrib_by_GU_? "true" : "false");
  fprintf(file, "swap_GUs :                   %s\n",  swap_GUs_? "true" : "false");

  // ------------------------------------------------------- Translation cost
  fprintf(file, "\nTranslation cost ----------------------------------------\n");
  fprintf(file, "translation_cost :           %e\n",  translation_cost_);
}

} // namespace aevol
