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


#ifndef AEVOL_PARAM_LOADER_H_
#define AEVOL_PARAM_LOADER_H_


// =================================================================
//                              Includes
// =================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <list>
#include <vector>
#include <set>
#include <memory>

#include <zlib.h>

#include "IOJson.h"
#include "ParameterLine.h"
#include "MutationParams.h"
#include "JumpingMT.h"
#include "macros.h"
#include "ae_enums.h"
#include "Gaussian.h"
#include "Point.h"
#include "Habitat.h"
#ifdef __REGUL
#include "raevol/Protein_R.h"
#endif

#ifndef LOGIN_NAME_MAX
#define LOGIN_NAME_MAX 256
#endif


namespace aevol {
// =================================================================
//                          Class declarations
// =================================================================
class ExpManager;
class Individual;


class ParamLoader {
 public :
  // =========================================================================
  //                          Constructors & Destructor
  // =========================================================================
    ParamLoader() = delete; //< Default ctor
    ParamLoader(const ParamLoader&) = delete; //< Copy ctor
    ParamLoader(ParamLoader&&) = delete; //< Move ctor
    ParamLoader(const char* file_name);
    virtual ~ParamLoader(); //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================
  /// Copy assignment
  ParamLoader& operator=(const ParamLoader& other) = delete;
  /// Move assignment
  ParamLoader& operator=(ParamLoader&& other) = delete;

  // =========================================================================
  //                             Public Methods
  // =========================================================================
  void load(ExpManager * exp_m, bool verbose = false,
            char* chromosome = nullptr, int32_t lchromosome = 0,
            char* plasmid = nullptr, int32_t lplasmid = 0
            #ifdef HAVE_MPI
            , int32_t nb_rank = 1
            #endif
            );
  void print_to_file(FILE* file);

  friend  IOJson::IOJson(const std::string &param_in, const std::string &chromosome);

  // =========================================================================
  //                                  Getters
  // =========================================================================

  // =========================================================================
  //                                  Setters
  // =========================================================================

  ParameterLine * get_line( int32_t* );

 protected :
  // =========================================================================
  //                            Protected Methods
  // =========================================================================
  void CheckConsistency();
  void read_file();
  ParameterLine* line(int32_t*);
  static void format_line(ParameterLine *, char*, bool*);
  void interpret_line(ParameterLine * line, int32_t cur_line);

  // =========================================================================
  //                               Attributes
  // =========================================================================
  std::shared_ptr<JumpingMT> prng_;

  char*   param_file_name_;
  FILE*   param_file_;

  // ----------------------------------------- PseudoRandom Number Generators
  // Seed for the selection random generator
  int32_t seed_;
  // Seed for the mutations random generator
  int32_t mut_seed_;
  // Seed for the stochasticity random generator
  int32_t stoch_seed_;
  // Seed for the phenotypic target variation random generator
  int32_t env_var_seed_;
  // Seed for the phenotypic target noise random generator
  int32_t env_noise_seed_;

  // ------------------------------------------------------------ Constraints
  int32_t min_genome_length_;
  int32_t max_genome_length_;
  double  w_max_;

  // ----------------------------------------------------- Initial conditions
  int32_t  chromosome_initial_length_;
  int8_t   init_method_;
  int32_t  init_pop_size_;
  char*    strain_name_;

  // -------------------------------------------------------- Phenotypic target
  std::list<Gaussian> std_env_gaussians;
  int16_t  env_sampling_;

  // ------------------------------------ Phenotypic target x-axis segmentation
  /// Number of x-axis segments
  int16_t env_axis_nb_segments_;
  /// x-axis segment boundaries (sorted -- including MIN_X and MAX_X)
  double* env_axis_segment_boundaries_;
  /// x-axis segment features
  PhenotypicFeature * env_axis_features_;
  /// Whether to automatically separate segments
  bool env_axis_separate_segments_;

  // ---------------------------------------------- Phenotypic target variation
  PhenotypicTargetVariationMethod env_var_method_;
  double      env_var_sigma_;
  int32_t     env_var_tau_;

  // -------------------------------------------------- Phenotypic target noise
  PhenotypicTargetNoiseMethod env_noise_method_;   // Method... TODO
  double  env_noise_alpha_;         // Alpha value (variance coefficient)
  double  env_noise_sigma_;         // Variance of the noise
  double  env_noise_prob_;          // Probability of variation.
  int32_t env_noise_sampling_log_;  // Log2 of the number of points in the noise fuzzy_set

  // --------------------------------------------------------- Mutation rates
  double  point_mutation_rate_;
  double  small_insertion_rate_;
  double  small_deletion_rate_;
  int16_t max_indel_size_;

  // -------------------------------------------- Rearrangements and Transfer
  bool    with_4pts_trans_;
  bool    with_alignments_;
  bool    with_HT_;
  bool    repl_HT_with_close_points_;
  double  HT_ins_rate_;
  double  HT_repl_rate_;
  double  repl_HT_detach_rate_;

  // ------------------------------ Rearrangement rates (without alignements)
  double duplication_rate_;
  double deletion_rate_;
  double translocation_rate_;
  double inversion_rate_;

  // --------------------------------- Rearrangement rates (with alignements)
  double neighbourhood_rate_;
  double duplication_proportion_;
  double deletion_proportion_;
  double translocation_proportion_;
  double inversion_proportion_;

  // ------------------------------------------------------------ Alignements
  AlignmentFunctionShape align_fun_shape_;
  double  align_sigm_lambda_;
  int16_t align_sigm_mean_;
  int16_t align_lin_min_;
  int16_t align_lin_max_;

  int16_t align_max_shift_;     // Maximum shift of one seq on the other
  int16_t align_w_zone_h_len_;  // Work zone half length
  int16_t align_match_bonus_;   // Corresponding residues match bonus
  int16_t align_mismatch_cost_; // Corresponding residues mismatch cost

  // ----------------------------------------------- Phenotypic Stochasticity
  bool with_stochasticity_;

  // -------------------------------------------------------------- Selection
  SelectionScheme selection_scheme_;
  double               selection_pressure_;

  SelectionScope  selection_scope_;
  int32_t               selection_scope_x_;
  int32_t               selection_scope_y_;

  FitnessFunction fitness_function_;
  int32_t               fitness_function_x_;
  int32_t               fitness_function_y_;
  // ------------------------------------------------------ Spatial structure
  int16_t grid_width_  = 32;
  int16_t grid_height_ = 32;
  bool    well_mixed = false;
  int32_t partial_mix_nb_permutations = 0;

  // -------------------------------------------------------------- Secretion
  bool   with_secretion_;
  // Proportion of the fitness contributed by secretion
  double secretion_contrib_to_fitness_;      // (0,1)
  // Proportion that diffuses into each cell, every generation
  // (0 for no diffusion)
  double secretion_diffusion_prop_;
  // Proportion of secreted substance that degrades every generation
  double secretion_degradation_prop_;
  // Cost of secreting the compound, as a proportion of the amount secreted
  double secretion_cost_;
  // Starting configuration of secretion grid
  // 0, all are 0; 1, point source of secreted compund
  double secretion_init_;

  // --------------------------------------------------------------- Plasmids
  bool      allow_plasmids_;
  int32_t   plasmid_initial_length_;
  int32_t   plasmid_initial_gene_;
  int32_t   plasmid_minimal_length_;
  int32_t   plasmid_maximal_length_;
  int32_t   chromosome_minimal_length_;
  int32_t   chromosome_maximal_length_;
  double    prob_plasmid_HT_;
  double    tune_donor_ability_;
  double    tune_recipient_ability_;
  double    donor_cost_;
  double    recipient_cost_;
  bool      compute_phen_contrib_by_GU_;
  bool      swap_GUs_;

  // ------------------------------------------------------- Translation cost
  double translation_cost_;

  #ifdef BASE_4

  // ------------------------------------------------------------ Terminators
  int8_t term_polya_seq_length_;

  // -------------------------------------------------- MWH bases configuration
  bool mwh_bases_redefined_ = false;
  int8_t aa_base_m_[NB_AMINO_ACIDS];
  int8_t aa_base_w_[NB_AMINO_ACIDS];
  int8_t aa_base_h_[NB_AMINO_ACIDS];

  std::map<std::string, AminoAcid> str_to_aminoacid_ = {
      {std::string("phe"), PHENYLALANINE},
      {std::string("leu"), LEUCINE},
      {std::string("iso"), ISOLEUCINE},
      {std::string("met"), METHIONINE},
      {std::string("val"), VALINE},
      {std::string("ser"), SERINE},
      {std::string("pro"), PROLINE},
      {std::string("thr"), THREONINE},
      {std::string("ala"), ALANINE},
      {std::string("tyr"), TYROSINE},
      {std::string("sto"), STOP},
      {std::string("his"), HISTIDINE},
      {std::string("glu"), GLUTAMINE},
      {std::string("asp"), ASPARAGINE},
      {std::string("lys"), LYSINE},
      {std::string("asa"), ASPARTIC_ACID},
      {std::string("gla"), GLUTAMIC_ACID},
      {std::string("cys"), CYSTEINE},
      {std::string("try"), TRYPTOPHAN},
      {std::string("arg"), ARGININE},
      {std::string("gly"), GLYCINE}
  };
  #endif

  // ---------------------------------------------------------------- Outputs
  // Stats
  int8_t stats_;
  // Whether to delete the existing statistics file
  // (otherwise kept with the suffix ".old")
  bool delete_old_stats_;

  // Backups
  int32_t backup_step_;
  int32_t big_backup_step_;

  // Tree
  bool record_tree_;
  int32_t tree_step_;


  // LightTree
  bool record_light_tree_;

  // Dumps // TODO : explain
  bool    make_dumps_;
  int32_t dump_step_;

  // Logs
  int8_t logs_;

  // Fuzzy set flavor
  int _fuzzy_flavor;

  // SIMD Metadata set flavor
  int simd_metadata_flavor_;

  // Other
  bool more_stats_;  // TODO : explain

  #ifdef __REGUL
    // Regulation factors
    double  _hill_shape_n;
    double  _hill_shape_theta;
    double  _hill_shape;

    // Degradation equation
    double  _degradation_rate;
    int     _nb_degradation_step;

    // Individual life
    int     _nb_indiv_age;

    // List of evaluation step
    std::set<int>  _list_eval_step;

    // Binding matrix
    double _binding_zeros_percentage;
    bool   _random_binding_matrix;

    // Heredity
    bool    _with_heredity;
    double  _protein_presence_limit;

    //Specific variatio method
    double _env_switch_probability;
    std::vector<std::list<Gaussian>> _env_gaussians_list;
    std::vector<std::list<int16_t>> _env_signals_list;
    std::vector<Protein_R*> _signals_models;
  #endif

  #ifdef HAVE_MPI
  // Define global grid size
  int32_t global_pop_size_;
  int32_t global_grid_width_;
  int32_t global_grid_height_;

  int32_t rank_width_;
  int32_t rank_height_;
  #endif
};

} // namespace aevol
#endif // AEVOL_PARAM_LOADER_H_
