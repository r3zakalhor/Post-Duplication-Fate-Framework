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


#ifndef AEVOL_INDIVIDUAL_H_
#define AEVOL_INDIVIDUAL_H_


// =================================================================
//                              Includes
// =================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <zlib.h>

#include <list>
#include <memory>

#include "Metrics.h"
#include "NonCodingMetrics.h"
#include "GeneticUnit.h"
#include "Dna.h"
#include "Rna.h"
#include "Protein.h"
#include "Phenotype.h"
#include "MutationParams.h"
#include "Fuzzy.h"
#include "PhenotypicTarget.h"
#include "Habitat.h"
#include "Observable.h"


namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================
class ExpManager;
class GridCell;
class Dna_7;

/// Models an individual cell.
///
/// Proteins and RNAs are shared with descent.
/// Genetic units are an individual's own.
class Individual : public Observable {
  friend class Dna;
    friend class Selection;
 public :

  // =================================================================
  //                             Constructors
  // =================================================================
  Individual() = delete;

  Individual(const Individual& other);

  Individual(ExpManager* exp_m,
             std::shared_ptr<JumpingMT> mut_prng,
             std::shared_ptr<JumpingMT> stoch_prng,
             std::shared_ptr<MutationParams> param_mut,
             double w_max,
             int32_t min_genome_length,
             int32_t max_genome_length,
             bool allow_plasmids,
             int32_t id,
             const char* strain_name,
             int32_t age);

  Individual(ExpManager* exp_m, gzFile backup_file);

  Individual(const Individual* parent, int32_t id,
             std::shared_ptr<JumpingMT> mut_prng,
             std::shared_ptr<JumpingMT> stoch_prng);

  static Individual* CreateIndividual(ExpManager* exp_m,
                                      gzFile backup_file);

  /**
   * \brief Create of clone of an Individual
   *
   * A clone is ... TODO
   *
   * \param dolly original individual to be cloned
   * \param id ID of the clone in the population
   * \return clone of dolly
   */
  static Individual* CreateClone(const Individual* dolly, int32_t id);

  // =================================================================
  //                             Destructors
  // =================================================================
  virtual ~Individual() noexcept(true);

  // =================================================================
  //                        Accessors: Getters
  // =================================================================
  const char* strain_name() const;

  int32_t age() const;

  ExpManager* exp_m() const;

  int16_t nb_genetic_units() const;

  int32_t nb_plasmids() const;

  int32_t amount_of_dna() const;

  const GeneticUnit& genetic_unit(int16_t num_unit) const;

  GeneticUnit& genetic_unit_nonconst(int16_t num_unit);

  double dist_to_target_by_feature(PhenotypicFeature feature) const;

  double fitness() const;

  double fitness_by_feature(PhenotypicFeature feature) const;

  GridCell* grid_cell() const;

  const Habitat& habitat() const;

  bool placed_in_population() const;

  std::shared_ptr<JumpingMT> mut_prng() const;

  std::shared_ptr<JumpingMT> stoch_prng() const;

  const std::shared_ptr<MutationParams>& mut_params() const;

  const std::list<GeneticUnit>& genetic_unit_list() const;

  std::list<GeneticUnit>& genetic_unit_list_nonconst();

  /// Keep only the first (main chromosome) and the last chromosomes from the GU.
  void drop_nested_genetic_units();

  const char* genetic_unit_sequence(int16_t num_unit) const;

  int32_t genetic_unit_seq_length(int16_t num_unit) const;

    unsigned long long id() const;

    unsigned long long long_id() const { return long_id_;};

  int32_t rank() const;


  Phenotype* phenotype() const;

  AbstractFuzzy* phenotype_activ() const;
  AbstractFuzzy* phenotype_inhib() const;

#ifndef __REGUL
  const PhenotypicTarget& phenotypic_target() const;
#endif

  double* dist_to_target_by_segment() const;

  const std::list<Protein*>& protein_list() const;

  const std::list<const Rna*>& rna_list() const;

  double w_max() const;

  // Genome size constraints
  int32_t min_genome_length() const;

  int32_t max_genome_length() const;

  // Plasmids settings
  bool allow_plasmids() const;

  // Mutation rates etc...
  double point_mutation_rate() const;

  double small_insertion_rate() const;

  double small_deletion_rate() const;

  int16_t max_indel_size() const;

  double duplication_rate() const;

  double deletion_rate() const;

  double translocation_rate() const;

  double inversion_rate() const;

  double neighbourhood_rate() const;

  double duplication_proportion() const;

  double deletion_proportion() const;

  double translocation_proportion() const;

  double inversion_proportion() const;

  // Transfer
  bool with_4pts_trans() const;

  bool with_HT() const;

  bool repl_HT_with_close_points() const;

  double HT_ins_rate() const;

  double HT_repl_rate() const;

  double repl_HT_detach_rate() const;

  // ------------------------------------------------------------ Alignements
  bool with_alignments() const;

  AlignmentFunctionShape align_fun_shape() const;

  double align_sigm_lambda() const;

  int16_t align_sigm_mean() const;

  int16_t align_lin_min() const;

  int16_t align_lin_max() const;

  // Maximum shift of one seq on the other
  int16_t align_max_shift() const;

  // Work zone half length
  int16_t align_w_zone_h_len() const;

  // Corresponding residues match bonus
  int16_t align_match_bonus() const;

  // Corresponding residues mismatch cost
  int16_t align_mismatch_cost() const;


  // ----------------------------------------------- Phenotypic stochasticity
  bool with_stochasticity() const;

  // Statistical data
  int32_t total_genome_size() const; // TODO: duplicate with amount_of_dna?
  int16_t nb_coding_RNAs() const;

  int16_t nb_non_coding_RNAs() const;

  int32_t overall_size_coding_RNAs() const;

  double av_size_coding_RNAs() const;

  int32_t overall_size_non_coding_RNAs() const;

  double av_size_non_coding_RNAs() const;

  int16_t nb_genes_activ() const;

  int16_t nb_genes_inhib() const;

  int16_t nb_functional_genes() const;

  int16_t nb_non_functional_genes() const;

  int32_t overall_size_functional_genes() const;

  double av_size_functional_genes() const;

  int32_t overall_size_non_functional_genes() const;

  double av_size_non_functional_genes() const;

  int32_t nb_bases_in_0_CDS() const;

  int32_t nb_bases_in_0_functional_CDS() const;

  int32_t nb_bases_in_0_non_functional_CDS() const;

  int32_t nb_bases_in_0_RNA() const;

  int32_t nb_bases_in_0_coding_RNA() const;

  int32_t nb_bases_in_0_non_coding_RNA() const;

  int32_t nb_bases_in_neutral_regions() const;

  int32_t nb_neutral_regions() const;

  double modularity(); // Not const

  int32_t* int_probes() const;

  double* double_probes() const;


  bool transcribed() const { return transcribed_; }

  // =================================================================
  //                        Accessors: Setters
  // =================================================================
  void set_strain_name(char* name);

  void set_exp_m(ExpManager* exp_m);

  void set_id(int32_t id);

  void set_rank(int32_t rank);

  void set_grid_cell(GridCell* grid_cell);

  void set_placed_in_population(bool placed_in_population);

  void reset_dist_to_target_by_segment(double* dist_to_target_by_segment);

  void set_w_max(double w_max);

  // Genome size constraints
  void set_min_genome_length(int32_t min_genome_length);

  void set_max_genome_length(int32_t max_genome_length);

  // Plasmids
  void set_allow_plasmids(bool allow_plasmids);

  // Mutation rates etc...
  void set_point_mutation_rate(double point_mutation_rate);

  void set_small_insertion_rate(double small_insertion_rate);

  void set_small_deletion_rate(double small_deletion_rate);

  void set_max_indel_size(int16_t max_indel_size);

  void set_duplication_rate(double duplication_rate);

  void set_deletion_rate(double deletion_rate);

  void set_translocation_rate(double translocation_rate);

  void set_inversion_rate(double inversion_rate);

  void set_neighbourhood_rate(double neighbourhood_rate);

  void set_duplication_proportion(double duplication_proportion);

  void set_deletion_proportion(double deletion_proportion);

  void set_translocation_proportion(double translocation_proportion);

  void set_inversion_proportion(double inversion_proportion);

  // Transfer
  void set_with_4pts_trans(bool with_4pts_trans);

  void set_with_alignments(bool with_alignments);

  void set_with_HT(bool with_HT);

  void set_repl_HT_with_close_points(bool repl_HT_with_close_points);

  void set_HT_ins_rate(double HT_ins_rate);

  void set_HT_repl_rate(double HT_repl_rate);

  void set_repl_HT_detach_rate(double repl_HT_detach_rate);

  // ----------------------------------------------- Phenotypic stochasticity
  void set_with_stochasticity(bool with_stoch);

  void set_mut_prng(std::shared_ptr<JumpingMT> prng);

  void set_stoch_prng(std::shared_ptr<JumpingMT> prng);

  //------------------------------------------------ Generic probes
  void set_int_probes(int32_t* int_probes);

  void set_double_probes(double* double_probes);

  void set_fitness(double fitness) { fitness_ = fitness;}
  void set_metaerror(double metaerror) { dist_to_target_by_feature_[METABOLISM] = metaerror; }

  // =================================================================
  //                            Public Methods
  // =================================================================
  // warning: the individual is left in a totally "cleared" state but not
  // reevaluated
  void add_GU(char*& sequence,  int32_t length);
  // void add_GU(GeneticUnit&& unit); // warning: the individual is left in
  // a totally "cleared" state but not reevaluated
  void add_GU(Individual* indiv, int32_t chromosome_length,
              std::shared_ptr<JumpingMT> prng);

  // void add_GU(Individual* indiv, int32_t length, JumpingMT* prng);
  // warning: the individual is left in a totally "cleared" state but not
  // reevaluated
  void remove_GU(int16_t num_unit);
  void renew_dist_to_target_by_feature();

  void renew_fitness_by_feature();

  void inject_GU(Individual* donor);

  void inject_2GUs(Individual* partner);

  /**
   * Main evaluation method
   */
  virtual void Evaluate(bool no_signal = true);

  /**
   * Evaluate within the provided context
   */
  virtual void EvaluateInContext(const Habitat& habitat);

  /**
   * Reset expression, proteins, phenotype and so on and evaluate
   */
  virtual void Reevaluate();

  /**
   * Reset expression, proteins, phenotype and so on and evaluate in context
   */
  virtual void ReevaluateInContext(const Habitat& habitat);

  virtual void clear_everything_except_dna_and_promoters();

  void do_transcription_translation_folding();

  void do_transcription();

  void do_translation();

  void do_folding();

  void compute_phenotype();

  void compute_distance_to_target(const PhenotypicTarget& target);

  // Computation of a "proper" fitness value (one that increases when the individual is fitter)
  // The behaviour of this function depends on many parameters and most notably on whether it is
  // a "composite" fitness or not, and on the selection scheme.
  void compute_fitness(const PhenotypicTarget& target);

  void compute_statistical_data();

  void compute_non_coding();

  virtual void save(gzFile backup_file) const;

  int32_t nb_terminators();

#ifdef DEBUG

  void assert_promoters();

  void assert_promoters_order();

#endif

  void compute_experimental_f_nu(int32_t nb_children,
                                 double* reproduction_statistics,
                                 double* offsprings_statistics = NULL,
                                 FILE* replication_file = NULL);

  double compute_theoritical_f_nu();
  // These functions compute the probability of neutral reproduction (F_nu).
  // The first method replicates the individual "nb_children" times and counts how often
  // a child has the same fitness as its parent (and if second argument not NULL, how often
  // a child has same or better fitness and stores the result at the adress contained in the
  // pointer). Results are proportions.
  // The second is an estimate based on genome structure as defined by Carole.
  // They have been implemented on the chromosome only !

  void remove_non_coding_bases();

  void double_non_coding_bases();

  void copy_parent(const Individual* parent, bool env_will_change);


  // =================================================================
  //                           Public Attributes
  // =================================================================



    GridCell* grid_cell_ = NULL;
    std::list<Protein*> protein_list_;
    int number_of_clones_ = 1;
    bool evaluated_;


        // The chromosome and plasmids (if allowed)
        std::list<GeneticUnit> genetic_unit_list_;

        int32_t parent_id_;
    protected :
  // =================================================================
  //                           Protected Methods
  // =================================================================
  virtual void make_protein_list();

  virtual void make_rna_list();

  // =================================================================
  //                          Protected Attributes
  // =================================================================
  ExpManager* exp_m_;

  // Name and "Age" of the strain
  char* strain_name_;
  int32_t age_;

  // Random number generators
  // These are shared pointers because depending on the configuration,
  // they can either be exclusive to the individual or grid cell, or they
  // can be mutualized (shared) for all or part of the population
  std::shared_ptr<JumpingMT> mut_prng_;
  std::shared_ptr<JumpingMT> stoch_prng_;

  // Individual ID and rank of the individual in the population
  // WARNING : The ID is no longer corresponding to the rank of the individual.
  //           The reason for this change is that we now need an identifier for the individuals
  //           as soon as they are created (the rank is only known when all the individuals have been evaluated).
  //           The rank will now be handled in a specific new attribute. (1 for the worst indiv, POP_SIZE for the best)
  unsigned long long id_ = -1;   // [0 ; POP_SIZE[
  unsigned long long long_id_ = -1;   // [0 ; POP_SIZE[

  int32_t rank_; // [1 ; POP_SIZE]

  // Total activation (resp. inhibition) of metabolic functions
  AbstractFuzzy* phenotype_activ_ = nullptr;
  AbstractFuzzy* phenotype_inhib_ = nullptr;

  // The phenotype, roughly corresponding to the sum of activ and inhib
  Phenotype* phenotype_ = nullptr;

  // Array containing the partial area of the "gap" (difference between the
  // phenotype and the target) for each phenotypic target segment.
  // Note: if the  phenotypic target is not segmented, there will be a single
  // cell in this array
  double* dist_to_target_by_segment_;

  // This array contains the aggregated area of the gap for each type of
  // segment (PhenotypicFeature).
  // When the  phenotypic target is not segmented, the only meaningful value
  // will be held in dist_to_target_[METABOLISM]
  double* dist_to_target_by_feature_;

  // This array contains the fitness contribution for each type of segment
  // (PhenotypicFeature).
  // For metabolic segment, this is the metabolic fitness contribution;
  // For secretion segment, this is the amount of compound that gets secreted
  // into the  phenotypic target;
  // Note: total fitness is the combination of metabolic fitness and the amount
  // of compound present in the habitat, not the amount of compound secreted by
  // the individual.
  double* fitness_by_feature_;

  // THE fitness
  double fitness_ = -1;

  // When using structured population, this is the cell the individual is in

  // int16_t x, y;


  // Access lists to all the proteins/RNAs of the individual.
  // Please note that these proteins/RNAs are actually managed (i.e. newed and deleted) via genetic units.

  std::list<const Rna*> rna_list_;

  // Generic probes
  int32_t* int_probes_ = nullptr; // Array of 5 int32_t values to be used as one wishes
  double* double_probes_ = nullptr; // Array of 5 double values to be used as one wishes

  // Mutation rates etc...
  std::shared_ptr<MutationParams> mut_params_;

  // ----------------------------------------------- Phenotypic stochasticity
  bool with_stochasticity_ = false;

  // Artificial chemistry
  double w_max_ = -1;

  // Genome size constraints
  int32_t min_genome_length_ = -1;
  int32_t max_genome_length_ = -1;

  // Plasmids settings
  bool allow_plasmids_ = false;


  // --------------------------------------------------
  // "State" of the individual
  // --------------------------------------------------
  // We keep trace of what we have already computed to avoid double computation (mainly in post-treaments)
  bool transcribed_ = false;
  bool translated_ = false;
  bool folded_ = false;
  bool phenotype_computed_ = false;
  bool distance_to_target_computed_ = false;
  bool fitness_computed_ = false;

  bool placed_in_population_ = true; // TODO: spatial ?



  // ----------------------------------------
  // Statistical data
  // ----------------------------------------
  // Genome, RNAs and genes size and stuff
  Metrics* metrics_ = nullptr;
  // Coding / non-coding
  NonCodingMetrics* nc_metrics_ = nullptr;

  // Mutation/Rearrangement statistics are managed in the replication report


  double modularity_ = -1; // Ratio between the pairwise distance between genes whose corresponding
  // phenotypic triangles overlap and the average intergenic distance
  // (ignoring non-functional genes)void compute_phenotype();



};

} // namespace aevol

#endif // AEVOL_INDIVIDUAL_H_
