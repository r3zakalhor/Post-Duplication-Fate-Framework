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

#include <cstdio>
#include <cmath>
#include <sys/stat.h>
#ifndef __OPENMP_GPU
#include <algorithm>
#else

#endif
#include <cassert>
#include <list>

#include<chrono>
using namespace std::chrono;
#include "ae_logger.h"

#include "Codon.h"
#include "ExpSetup.h"
#include "ExpManager.h"
#include "GridCell.h"
#include "GeneticUnit.h"
#include "VisAVis.h"
#include "Utils.h"
#include "HybridFuzzy.h"

#ifdef __NO_X
  #ifdef __REGUL
    #include "raevol/Individual_R.h"
  #else
    #include "Individual.h"
  #endif
#elif defined __X11
  #ifdef __REGUL
    #include "raevol/Individual_R_X11.h"
  #else
    #include "Individual_X11.h"
  #endif
#endif
namespace aevol {

// ===========================================================================
//                             Constructors
// ===========================================================================
/**
 * // TODO <david.parsons@inria.fr>
 */
Individual::Individual(ExpManager* exp_m,
                       std::shared_ptr<JumpingMT> mut_prng,
                       std::shared_ptr<JumpingMT> stoch_prng,
                       std::shared_ptr<MutationParams> param_mut,
                       double w_max,
                       int32_t min_genome_length,
                       int32_t max_genome_length,
                       bool allow_plasmids,
                       int32_t id,
                       const char* strain_name,
                       int32_t age) {
  // Experiment manager
  exp_m_ = exp_m;

  // PRNGs
  mut_prng_ = mut_prng;
  stoch_prng_ = stoch_prng;

  // ID and rank of the indiv ; name and "age" of the strain
  set_id(id);
  rank_ = -1; // TODO: UNRANKED
  age_ = age;
  strain_name_ = new char[strlen(strain_name) + 1];
  strcpy(strain_name_, strain_name);

  phenotype_activ_ = NULL;
  phenotype_inhib_ = NULL;
  phenotype_ = NULL;

  dist_to_target_by_segment_ = NULL;
  dist_to_target_by_feature_ = new double[NB_FEATURES];
  fitness_by_feature_ = new double[NB_FEATURES];
  for (int i = 0; i < NB_FEATURES; i++) {
    dist_to_target_by_feature_[i] = 0;
    fitness_by_feature_[i] = 0;
  }

  fitness_ = 0.0;

  // TODO <david.parsons@inria.fr> Deal with cell coordinates on the grid
  // When using structured population, this is the cell the individual is in
  // x = y = -1;

  // The chromosome and plasmids (if allowed)
  // TODO <david.parsons@inria.fr> ???

  // Generic probes
  int_probes_ = new int32_t[5];
  double_probes_ = new double[5];
  for (int8_t i = 0; i < 5; i++) {
    int_probes_[i] = 0;
    double_probes_[i] = 0.0;
  }

  // Mutation rates etc...
  mut_params_ = param_mut;

  // Artificial chemistry
  w_max_ = w_max;

  // Genome size constraints
  min_genome_length_ = min_genome_length;
  max_genome_length_ = max_genome_length;

  // Plasmids settings
  allow_plasmids_ = allow_plasmids;


  // --------------------------------------------------
  // "State" of the individual
  // --------------------------------------------------
  evaluated_ = false;
  transcribed_ = false;
  translated_ = false;
  folded_ = false;
  phenotype_computed_ = false;
  distance_to_target_computed_ = false;
  fitness_computed_ = false;
  placed_in_population_ = false;



  // ----------------------------------------
  // Statistical data
  // ----------------------------------------
  modularity_ = 0.0;
}

/**
 * This constructor retrieves an individual from a backup file.
 *
 * Since this generation has already been processed, no unnecessary calculation
 * (e.g. fitness) will be done.
 * No transcription, translation or other process of that kind is performed.
 */
Individual::Individual(ExpManager* exp_m, gzFile backup_file) {
  exp_m_ = exp_m;

  // Retrieve the name and "age" of the strain
  int8_t strain_string_len;
  gzread(backup_file, &strain_string_len, sizeof(strain_string_len));
  strain_name_ = new char[strain_string_len + 1];
  gzread(backup_file, strain_name_, strain_string_len + 1);
  gzread(backup_file, &age_, sizeof(age_));

  // Retrieve the PRNGs
  if (exp_m == NULL) {
    // Detached mode
    mut_prng_ = NULL;
    stoch_prng_ = NULL;
  }
  else {
    // TODO: => prngs as parameters
    mut_prng_ = exp_m->world()->mut_prng();
    stoch_prng_ = exp_m->world()->stoch_prng();
    assert(mut_prng_);
    assert(stoch_prng_);
  }

  // Retrieve id and rank
  gzread(backup_file, &id_, sizeof(id_));

      gzread(backup_file, &long_id_, sizeof(long_id_));
  gzread(backup_file, &rank_, sizeof(rank_));

  // Retrieve spatial coordinates
  // gzread(backup_file, &x, sizeof(x));
  // gzread(backup_file, &y, sizeof(y));
  placed_in_population_ = false;

  // Retrieve generic probes
  int_probes_ = new int32_t[5];
  double_probes_ = new double[5];
  gzread(backup_file, int_probes_, 5 * sizeof(*int_probes_));
  gzread(backup_file, double_probes_, 5 * sizeof(*double_probes_));

  // Retrieve mutational parameters
  mut_params_ = std::make_shared<MutationParams>(backup_file);

  // ------------------------------------------------- Phenotypic stochasticity
  gzread(backup_file, &with_stochasticity_, sizeof(with_stochasticity_));

  // Retrieve artificial chemistry parameters
  gzread(backup_file, &w_max_, sizeof(w_max_));

  // Retrieve genome size constraints
  gzread(backup_file, &min_genome_length_, sizeof(min_genome_length_));
  gzread(backup_file, &max_genome_length_, sizeof(max_genome_length_));

  // Retrieve plasmids settings
  int8_t tmp_allow_plasmids;
  gzread(backup_file, &tmp_allow_plasmids, sizeof(tmp_allow_plasmids));
  allow_plasmids_ = tmp_allow_plasmids ? 1 : 0;

  // Retrieve genetic units
  int16_t nb_gen_units;
  gzread(backup_file, &nb_gen_units, sizeof(nb_gen_units));

  int lid = 0;
  for (int16_t i = 0; i < nb_gen_units; i++) {
    genetic_unit_list_.emplace_back(this, backup_file);
    genetic_unit_list_.back().set_local_id(lid++);
  }

  // --------------------------------------------------------------------------
  // No more data to retrieve, the following are only structural
  // initializations (no data is set)
  // --------------------------------------------------------------------------
  // --------------------------------------------------------------------------

  // Create empty fuzzy sets for activation and inhibition
  phenotype_activ_ = NULL;
  phenotype_inhib_ = NULL;
  phenotype_ = NULL;

  dist_to_target_by_segment_ = NULL;
  dist_to_target_by_feature_ = new double[NB_FEATURES];
  fitness_by_feature_ = new double[NB_FEATURES];

  for (int8_t i = 0; i < NB_FEATURES; i++) {
    dist_to_target_by_feature_[i] = 0.0;
    fitness_by_feature_[i] = 0.0;
  }

  // Initialize the computational state of the individual
  evaluated_ = false;
  transcribed_ = false;
  translated_ = false;
  folded_ = false;
  phenotype_computed_ = false;
  distance_to_target_computed_ = false;
  fitness_computed_ = false;

  modularity_ = -1;
}

/**
 * Copy constructor
 */
Individual::Individual(const Individual& other) {
  exp_m_ = other.exp_m_;

  // PRNGs
  mut_prng_ = other.mut_prng_;
  stoch_prng_ = other.stoch_prng_;

  int strain_string_len = strlen(other.strain_name_);
  strain_name_ = new char[strain_string_len + 1];
  memcpy(strain_name_, other.strain_name_, strain_string_len + 1);
  age_ = other.age_;

  id_ = other.id_;
  long_id_ = other.long_id_;
  rank_ = other.rank_;

  evaluated_ = other.evaluated_;
  transcribed_ = other.transcribed_;
  translated_ = other.translated_;
  folded_ = other.folded_;
  phenotype_computed_ = other.phenotype_computed_;

  with_stochasticity_ = other.with_stochasticity_;

  // Artificial chemistry parameters
  w_max_ = other.w_max_;

  distance_to_target_computed_ = other.distance_to_target_computed_;
  fitness_computed_ = other.fitness_computed_;
  placed_in_population_ = other.placed_in_population_;

  // Copy genetic units from other
  // Should actually use GeneticUnit copy ctor which is disabled.
  for (const auto& gu: other.genetic_unit_list_)
    genetic_unit_list_.emplace_back(this, gu);

  // Copy phenotype
  if (phenotype_computed_) {
    phenotype_activ_  = FuzzyFactory::fuzzyFactory->create_fuzzy((*(other.phenotype_activ_)));
    phenotype_inhib_  = FuzzyFactory::fuzzyFactory->create_fuzzy((*(other.phenotype_inhib_)));
    phenotype_        = FuzzyFactory::fuzzyFactory->create_fuzzy((*(other.phenotype_)));
  }
  else {
    phenotype_activ_ = NULL;
    phenotype_inhib_ = NULL;
    phenotype_ = NULL;
  }


  // Copy fitness-related stuff
  dist_to_target_by_segment_ = NULL;
  dist_to_target_by_feature_ = new double[NB_FEATURES];
  fitness_by_feature_ = new double[NB_FEATURES];

  for (int8_t i = 0; i < NB_FEATURES; i++) {
    dist_to_target_by_feature_[i] = other.dist_to_target_by_feature_[i];
    fitness_by_feature_[i] = other.fitness_by_feature_[i];
  }

  fitness_ = other.fitness_;

  // Copy statistical data
  metrics_ = other.metrics_ ?
             new Metrics(*other.metrics_) :
             nullptr;
  nc_metrics_ = other.nc_metrics_ ?
                new NonCodingMetrics(*other.nc_metrics_) :
                nullptr;

  modularity_ = other.modularity_;

  // Generic probes
  int_probes_ = new int32_t[5];
  double_probes_ = new double[5];
  for (int8_t i = 0; i < 5; i++) {
    int_probes_[i] = other.int_probes_[i];
    double_probes_[i] = other.double_probes_[i];
  }

  // Mutation rates etc...
  mut_params_ = std::make_shared<MutationParams>(*(other.mut_params_));

  // Genome size constraints
  min_genome_length_ = other.min_genome_length_;
  max_genome_length_ = other.max_genome_length_;

  // Plasmids settings
  allow_plasmids_ = other.allow_plasmids_;


  grid_cell_ = other.grid_cell_;
}

/**
 * Reproduction constructor
 *
 * This constructor creates a new individual with the same genome as it's
 * parent. The location of promoters will be copied but no further process will
 * be performed.
 *
 * The phenotype and the fitness are not set, neither is the statistical data.
*/
Individual::Individual(const Individual* parent, int32_t id,
                       std::shared_ptr<JumpingMT> mut_prng,
                       std::shared_ptr<JumpingMT> stoch_prng) {
  exp_m_ = parent->exp_m_;

  // PRNGs
  mut_prng_ = mut_prng;
  stoch_prng_ = stoch_prng;

  int strain_string_len = strlen(parent->strain_name_);
  strain_name_ = new char[strain_string_len + 1];
  memcpy(strain_name_, parent->strain_name_, strain_string_len + 1);
  age_ = parent->age_ + 1;

  id_ = id;
  long_id_ = id_+AeTime::time()*exp_m_->nb_indivs();
  parent_id_ = parent->id();
  rank_ = -1;

  evaluated_ = false;
  transcribed_ = false;
  translated_ = false;
  folded_ = false;
  phenotype_computed_ = false;
  distance_to_target_computed_ = false;
  fitness_computed_ = false;

  placed_in_population_ = false;
  // x = y = -1;

  with_stochasticity_ = parent->with_stochasticity_;

  // Artificial chemistry parameters
  w_max_ = parent->w_max_;

  // Create new genetic units with their DNA copied from here
  // NOTE : The RNA lists (one per genetic unit) will also be copied so that we don't
  // need to look for promoters on the whole genome
  for (auto& gu: parent->genetic_unit_list_) {
    genetic_unit_list_.emplace_back(this, &gu);
    //genetic_unit_list_.back().set_local_id(gu.local_id());
  }

  phenotype_activ_ = NULL;
  phenotype_inhib_ = NULL;
  phenotype_ = NULL;

  // Initialize all the fitness-related stuff
  dist_to_target_by_segment_ = NULL;
  dist_to_target_by_feature_ = new double[NB_FEATURES];
  fitness_by_feature_ = new double[NB_FEATURES];

  for (int8_t i = 0; i < NB_FEATURES; i++) {
    dist_to_target_by_feature_[i] = 0.0;
    fitness_by_feature_[i] = 0.0;
  }

  // Generic probes
  int_probes_ = new int32_t[5];
  double_probes_ = new double[5];
  for (int8_t i = 0; i < 5; i++) {
    int_probes_[i] = parent->int_probes_[i];
    double_probes_[i] = parent->double_probes_[i];
  }

  // Mutation rates etc...
  mut_params_ = std::make_shared<MutationParams>(*(parent->mut_params_));

  // Genome size constraints
  min_genome_length_ = parent->min_genome_length_;
  max_genome_length_ = parent->max_genome_length_;

  // Plasmids settings
  allow_plasmids_ = parent->allow_plasmids_;

  // Initialize statistical data
  modularity_ = -1;

  grid_cell_ = parent->grid_cell_;
}

Individual* Individual::CreateIndividual(ExpManager* exp_m,
                                         gzFile backup_file) {
  Individual* indiv = NULL;
  #ifdef __NO_X
    #ifndef __REGUL
      indiv = new Individual(exp_m, backup_file);
    #else
      indiv = new Individual_R(exp_m, backup_file);
    #endif
  #elif defined __X11
    #ifndef __REGUL
      indiv = new Individual_X11(exp_m, backup_file);
    #else
      indiv = new Individual_R_X11(exp_m, backup_file);
    #endif
  #endif

  return indiv;
}

/*!
  \brief Create a clone

  \param dolly original individual to be cloned
  \param id ID of the clone
  \return clone of dolly (not evaluated)
*/
Individual* Individual::CreateClone(const Individual* dolly, int32_t id) {
  Individual* indiv = new Individual(*dolly);
  indiv->set_id(id);
  return indiv;
}


// =================================================================
//                             Destructor
// =================================================================
Individual::~Individual() noexcept {
  clear_everything_except_dna_and_promoters();
  
  delete[] strain_name_;

  // Proteins and RNAs are recycled, don't delete them.

  delete phenotype_activ_;
  delete phenotype_inhib_;
  delete phenotype_;

  if (dist_to_target_by_segment_ != NULL) delete[] dist_to_target_by_segment_;
  delete[] dist_to_target_by_feature_;

  delete[] fitness_by_feature_;

  // Generic probes
  delete[] int_probes_;
  delete[] double_probes_;

  delete metrics_;
  delete nc_metrics_;
}

// =================================================================
//                        Non-inline Accessors
// =================================================================
void Individual::set_exp_m(ExpManager* exp_m) {
  exp_m_ = exp_m;

  // Update pointer to exp_manager in each GU
  for (auto& gen_unit: genetic_unit_list_)
    gen_unit.set_exp_m(exp_m_);
}

/// TODO
void Individual::set_grid_cell(GridCell* grid_cell) {
  grid_cell_ = grid_cell;
  placed_in_population_ = true;
  // x = grid_cell->x();
  // y = grid_cell->y();
  if (grid_cell->individual() != this) {
    grid_cell->set_individual(this);
  }
}

/// TODO
const char* Individual::strain_name() const {
  return strain_name_;
}

/// TODO
int32_t Individual::age() const {
  return age_;
}

/// TODO
unsigned long long Individual::id() const {
  return id_;
}

/// TODO
double* Individual::dist_to_target_by_segment() const {
  return dist_to_target_by_segment_;
}

/*!
  Get the individual's rank in the population (1 for the worst indiv, POP_SIZE for the best)

  Warning: be sure you call sort_individuals() before using rank_in_population
*/
int32_t Individual::rank() const {
  return rank_;
}

/// TODO
ExpManager* Individual::exp_m() const {
  return exp_m_;
}

/// TODO
std::shared_ptr<JumpingMT> Individual::mut_prng() const {
  return mut_prng_;
}

/// TODO
std::shared_ptr<JumpingMT> Individual::stoch_prng() const {
  return stoch_prng_;
}

const std::shared_ptr<MutationParams>& Individual::mut_params() const {
  return mut_params_;
}

/*!
  Returns the number of genetic units
*/
int16_t Individual::nb_genetic_units() const {
  return genetic_unit_list_.size();
}

/// Get the number of plasmids. That is, the number of genetic units
/// minus one DNA-based genetic unit.
int32_t Individual::nb_plasmids() const {
  return genetic_unit_list_.size() - 1;
}

/// TODO
int32_t Individual::amount_of_dna() const {
  int32_t amount = 0;
  for (const auto& gen_unit: genetic_unit_list_)
    amount += gen_unit.dna()->length();
  return amount;
}

/// Return the list of genetic units.
const std::list<GeneticUnit>& Individual::genetic_unit_list() const {
  return genetic_unit_list_;
}

std::list<GeneticUnit>& Individual::genetic_unit_list_nonconst() {
  return genetic_unit_list_;
}

/// Remove all the elements from the GU list except the firt and the
/// last ones. If the GU list has less that 2 elements, do nothing.
void Individual::drop_nested_genetic_units() {
  if (genetic_unit_list_.size() <= 2) {
    return;
  }

  genetic_unit_list_.erase(std::next(genetic_unit_list_.begin()),
                           std::prev(genetic_unit_list_.end()));
  assert(genetic_unit_list_.size() <= 2);
}

/// Returns genetic unit number `num_unit` (0 for main chromosome)
const GeneticUnit& Individual::genetic_unit(int16_t num_unit) const {
  assert(num_unit < static_cast<int32_t>(genetic_unit_list_.size()));
  auto it = genetic_unit_list_.cbegin();
  std::advance(it, num_unit);
  return *it;
}

/// Returns genetic unit number `num_unit` (0 for main chromosome) as
/// a non-constant reference. To be used when the purpose is to alter
/// the individual.
GeneticUnit& Individual::genetic_unit_nonconst(int16_t num_unit) {
  assert(num_unit < static_cast<int32_t>(genetic_unit_list_.size()));
  auto it = genetic_unit_list_.begin();
  std::advance(it, num_unit);
  return *it;
}

/// TODO
double Individual::dist_to_target_by_feature(
    PhenotypicFeature feature) const {
  assert(distance_to_target_computed_);

  return dist_to_target_by_feature_[feature];
}

/// TODO
double Individual::fitness() const {
//  assert(fitness_computed_);

  return fitness_;
}

/// TODO
double Individual::fitness_by_feature(PhenotypicFeature feature) const {
  assert(fitness_computed_);

  return fitness_by_feature_[feature];
}

/// TODO
GridCell* Individual::grid_cell() const {
  return grid_cell_;
}

/// TODO
const Habitat& Individual::habitat() const {
  return grid_cell_->habitat();
}

/// TODO
bool Individual::placed_in_population() const {
  return placed_in_population_;
}

/*!
  Returns the sequence of genetic unit number <num_unit> (0 for main chromosome)
*/
const char* Individual::genetic_unit_sequence(int16_t num_unit) const {
  return genetic_unit(num_unit).sequence();
}

/*!
  Returns the sequence length of genetic unit number <num_unit> (0 for main chromosome)
*/
int32_t Individual::genetic_unit_seq_length(int16_t num_unit) const {
  return genetic_unit(num_unit).seq_length();
}

/// TODO
AbstractFuzzy* Individual::phenotype_activ() const {
  return phenotype_activ_;
}

/// TODO
AbstractFuzzy* Individual::phenotype_inhib() const {
  return phenotype_inhib_;
}

/// TODO
Phenotype* Individual::phenotype() const {
  return phenotype_;
}

#ifndef __REGUL
const PhenotypicTarget& Individual::phenotypic_target() const {
  return grid_cell_->phenotypic_target();
}
#endif

/// TODO
const std::list<Protein*>& Individual::protein_list() const {
  return protein_list_;
}

/// TODO
const std::list<const Rna*>& Individual::rna_list() const {
  return rna_list_;
}

/// TODO
int32_t Individual::total_genome_size() const {
  assert(metrics_ != nullptr);
  return metrics_->total_genome_size();
}

/// TODO
int16_t Individual::nb_coding_RNAs() const {
  assert(metrics_ != nullptr);
  return metrics_->nb_coding_RNAs();
}

/// TODO
int16_t Individual::nb_non_coding_RNAs() const {
  assert(metrics_ != nullptr);
  return metrics_->nb_non_coding_RNAs();
}

/// TODO
int32_t Individual::overall_size_coding_RNAs() const {
  assert(metrics_ != nullptr);
  return metrics_->overall_size_coding_RNAs();
}

/// TODO
double Individual::av_size_coding_RNAs() const {
  assert(metrics_ != nullptr);
  return metrics_->nb_coding_RNAs() == 0 ?
         0.0 :
         metrics_->overall_size_coding_RNAs() /
         metrics_->nb_coding_RNAs();
}

/// TODO
int32_t Individual::overall_size_non_coding_RNAs() const {
  assert(metrics_ != nullptr);
  return metrics_->overall_size_non_coding_RNAs();
}

/// TODO
double Individual::av_size_non_coding_RNAs() const {
  assert(metrics_ != nullptr);
  return metrics_->nb_non_coding_RNAs() == 0 ?
         0.0 :
         metrics_->overall_size_non_coding_RNAs() /
         metrics_->nb_non_coding_RNAs();
}

/// TODO
int16_t Individual::nb_genes_activ() const {
  assert(metrics_ != nullptr);
  return metrics_->nb_genes_activ();
}

/// TODO
int16_t Individual::nb_genes_inhib() const {
  assert(metrics_ != nullptr);
  return metrics_->nb_genes_inhib();
}

/// TODO
int16_t Individual::nb_functional_genes() const {
  assert(metrics_ != nullptr);
  return metrics_->nb_functional_genes();
}

/// TODO
int16_t Individual::nb_non_functional_genes() const {
  assert(metrics_ != nullptr);
  return metrics_->nb_non_functional_genes();
}

/// TODO
int32_t Individual::overall_size_functional_genes() const {
  assert(metrics_ != nullptr);
  return metrics_->overall_size_functional_genes();
}

/// TODO
double Individual::av_size_functional_genes() const {
  assert(metrics_ != nullptr);
  return metrics_->nb_functional_genes() == 0 ?
         0.0 :
         metrics_->overall_size_functional_genes() /
         metrics_->nb_functional_genes();
}

/// TODO
int32_t Individual::overall_size_non_functional_genes() const {
  assert(metrics_ != nullptr);
  return metrics_->overall_size_non_functional_genes();
}

/// TODO
double Individual::av_size_non_functional_genes() const {
  assert(metrics_ != nullptr);
  return metrics_->nb_non_functional_genes() == 0 ?
         0.0 :
         metrics_->overall_size_non_functional_genes() /
         metrics_->nb_non_functional_genes();
}

/// TODO
int32_t Individual::nb_bases_in_0_CDS() const {
  assert(nc_metrics_ != nullptr);
  return nc_metrics_->nb_bases_in_0_CDS();
}

/// TODO
int32_t Individual::nb_bases_in_0_functional_CDS() const {
  assert(nc_metrics_ != nullptr);
  return nc_metrics_->nb_bases_in_0_functional_CDS();
}

/// TODO
int32_t Individual::nb_bases_in_0_non_functional_CDS() const {
  assert(nc_metrics_ != nullptr);
  return nc_metrics_->nb_bases_in_0_non_functional_CDS();
}

/// TODO
int32_t Individual::nb_bases_in_0_RNA() const {
  assert(nc_metrics_ != nullptr);
  return nc_metrics_->nb_bases_in_0_RNA();
}

/// TODO
int32_t Individual::nb_bases_in_0_coding_RNA() const {
  assert(nc_metrics_ != nullptr);
  return nc_metrics_->nb_bases_in_0_coding_RNA();
}

/// TODO
int32_t Individual::nb_bases_in_0_non_coding_RNA() const {
  assert(nc_metrics_ != nullptr);
  return nc_metrics_->nb_bases_in_0_non_coding_RNA();
}

/// TODO
int32_t Individual::nb_bases_in_neutral_regions() const {
  assert(nc_metrics_ != nullptr);
  return nc_metrics_->nb_bases_in_neutral_regions();
}

/// TODO
int32_t Individual::nb_neutral_regions() const {
  assert(nc_metrics_ != nullptr);
  return nc_metrics_->nb_neutral_regions();
}

/// TODO
double Individual::modularity() {
  printf("\n  WARNING : modularity measure not yet implemented.\n");
  //~ if (modularity_ < 0) compute_modularity();
  //~ return modularity_;
  return 0;
}

/// TODO
double Individual::w_max() const {
  return w_max_;
}

// ------------------------------------------------------------- Mutation rates
/// TODO
double Individual::point_mutation_rate() const {
  return mut_params_->point_mutation_rate();
}

/// TODO
double Individual::small_insertion_rate() const {
  return mut_params_->small_insertion_rate();
}

/// TODO
double Individual::small_deletion_rate() const {
  return mut_params_->small_deletion_rate();
}

/// TODO
int16_t Individual::max_indel_size() const {
  return mut_params_->max_indel_size();
}

// ---------------------------------- Rearrangement rates (without alignements)
/// TODO
double Individual::duplication_rate() const {
  return mut_params_->duplication_rate();
}

/// TODO
double Individual::deletion_rate() const {
  return mut_params_->deletion_rate();
}

/// TODO
double Individual::translocation_rate() const {
  return mut_params_->translocation_rate();
}

/// TODO
double Individual::inversion_rate() const {
  return mut_params_->inversion_rate();
}

// ------------------------------------- Rearrangement rates (with alignements)
/// TODO
double Individual::neighbourhood_rate() const {
  return mut_params_->neighbourhood_rate();
}

/// TODO
double Individual::duplication_proportion() const {
  return mut_params_->duplication_proportion();
}

/// TODO
double Individual::deletion_proportion() const {
  return mut_params_->deletion_proportion();
}

/// TODO
double Individual::translocation_proportion() const {
  return mut_params_->translocation_proportion();
}

/// TODO
double Individual::inversion_proportion() const {
  return mut_params_->inversion_proportion();
}

// ---------------------------------------------------------------- Transfer
bool Individual::with_4pts_trans() const {
  return mut_params_->with_4pts_trans();
}

bool Individual::with_HT() const {
  return mut_params_->with_HT();
}

bool Individual::repl_HT_with_close_points() const {
  return mut_params_->repl_HT_with_close_points();
}

double Individual::HT_ins_rate() const {
  return mut_params_->HT_ins_rate();
}

double Individual::HT_repl_rate() const {
  return mut_params_->HT_repl_rate();
}

double Individual::repl_HT_detach_rate() const {
  return mut_params_->repl_HT_detach_rate();
}


// ---------------------------------------------------------------- Alignements
bool Individual::with_alignments() const {
  return mut_params_->with_alignments();
}

AlignmentFunctionShape Individual::align_fun_shape() const {
  return mut_params_->align_fun_shape();
}

double Individual::align_sigm_lambda() const {
  return mut_params_->align_sigm_lambda();
}

int16_t Individual::align_sigm_mean() const {
  return mut_params_->align_sigm_mean();
}

int16_t Individual::align_lin_min() const {
  return mut_params_->align_lin_min();
}

int16_t Individual::align_lin_max() const {
  return mut_params_->align_lin_max();
}

int16_t Individual::align_max_shift() const {
  return mut_params_->align_max_shift();
}

int16_t Individual::align_w_zone_h_len() const {
  return mut_params_->align_w_zone_h_len();
}

int16_t Individual::align_match_bonus() const {
  return mut_params_->align_match_bonus();
}

int16_t Individual::align_mismatch_cost() const {
  return mut_params_->align_mismatch_cost();
}

/// TODO
bool Individual::with_stochasticity() const {
  return with_stochasticity_;
}

void Individual::set_allow_plasmids(bool allow_plasmids) {
  allow_plasmids_ = allow_plasmids;
}

// Genome size constraints
/// TODO
int32_t Individual::min_genome_length() const {
  return min_genome_length_;
}

/// TODO
int32_t Individual::max_genome_length() const {
  return max_genome_length_;
}

// Plasmids settings
/// TODO
bool Individual::allow_plasmids() const {
  return allow_plasmids_;
}

/*!
  \brief Return the int_probes_

  \return int_probes_
*/
int32_t* Individual::int_probes() const {
  return int_probes_;
}

/*!
  \brief Return the double_probes_

  \return double_probes_
*/
double* Individual::double_probes() const {
  return double_probes_;
}


// =====================================================================
//                           Setters' definitions
// =====================================================================
void Individual::set_strain_name(char* name) {
  assert(name && strlen(name) < INT8_MAX); // Conservative, could be <=
  int8_t name_len = strlen(name);
  delete[] strain_name_;
  strain_name_ = new char[name_len + 1];
  memcpy(strain_name_, name, name_len + 1);
}

/// TODO
void Individual::set_id(int32_t id) {
  id_ = id;
}

/// TODO
void Individual::set_rank(int32_t rank) {
  rank_ = rank;
}

/// TODO
void Individual::set_placed_in_population(bool placed_in_population) {
  placed_in_population_ = placed_in_population;
}

/// TODO
void Individual::set_w_max(double w_max) {
  w_max_ = w_max;
}


// Genome size constraints
/// TODO
void Individual::set_min_genome_length(int32_t min_genome_length) {
  min_genome_length_ = min_genome_length;
}

/// TODO
void Individual::set_max_genome_length(int32_t max_genome_length) {
  max_genome_length_ = max_genome_length;
}


void Individual::set_point_mutation_rate(double point_mutation_rate) {
  mut_params_->set_point_mutation_rate(point_mutation_rate);
}

void Individual::set_small_insertion_rate(double small_insertion_rate) {
  mut_params_->set_small_insertion_rate(small_insertion_rate);
}

void Individual::set_small_deletion_rate(double small_deletion_rate) {
  mut_params_->set_small_deletion_rate(small_deletion_rate);
}

void Individual::set_max_indel_size(int16_t max_indel_size) {
  mut_params_->set_max_indel_size(max_indel_size);
}

void Individual::set_duplication_rate(double duplication_rate) {
  mut_params_->set_duplication_rate(duplication_rate);
}

void Individual::set_deletion_rate(double deletion_rate) {
  mut_params_->set_deletion_rate(deletion_rate);
}

void Individual::set_translocation_rate(double translocation_rate) {
  mut_params_->set_translocation_rate(translocation_rate);
}

void Individual::set_inversion_rate(double inversion_rate) {
  mut_params_->set_inversion_rate(inversion_rate);
}

void Individual::set_neighbourhood_rate(double neighbourhood_rate) {
  mut_params_->set_neighbourhood_rate(neighbourhood_rate);
}

void Individual::set_duplication_proportion(double duplication_proportion) {
  mut_params_->set_duplication_proportion(duplication_proportion);
}

void Individual::set_deletion_proportion(double deletion_proportion) {
  mut_params_->set_deletion_proportion(deletion_proportion);
}

void Individual::set_translocation_proportion(double translocation_proportion) {
  mut_params_->set_translocation_proportion(translocation_proportion);
}

void Individual::set_inversion_proportion(double inversion_proportion) {
  mut_params_->set_inversion_proportion(inversion_proportion);
}

void Individual::set_with_4pts_trans(bool with_4pts_trans) {
  mut_params_->set_with_4pts_trans(with_4pts_trans);
}

void Individual::set_with_alignments(bool with_alignments) {
  mut_params_->set_with_alignments(with_alignments);
}

void Individual::set_with_HT(bool with_HT) {
  mut_params_->set_with_HT(with_HT);
}

void Individual::set_repl_HT_with_close_points(bool repl_HT_with_close_points) {
  mut_params_->set_repl_HT_with_close_points(repl_HT_with_close_points);
}

void Individual::set_HT_ins_rate(double HT_ins_rate) {
  mut_params_->set_HT_ins_rate(HT_ins_rate);
}

void Individual::set_HT_repl_rate(double HT_repl_rate) {
  mut_params_->set_HT_repl_rate(HT_repl_rate);
}

void Individual::set_repl_HT_detach_rate(double repl_HT_detach_rate) {
  mut_params_->set_repl_HT_detach_rate(repl_HT_detach_rate);
}


void Individual::set_with_stochasticity(bool with_stoch) {
  with_stochasticity_ = with_stoch;
}

void Individual::set_stoch_prng(std::shared_ptr<JumpingMT> prng) {
  stoch_prng_ = prng;
}

void Individual::set_mut_prng(std::shared_ptr<JumpingMT> prng) {
  mut_prng_ = prng;
}

/*!
  \brief Change the int_probes_

  \param int_probes 5 int32_t* that constitute a probe
*/
void Individual::set_int_probes(int32_t* int_probes) {
  int_probes_ = int_probes;
}

/*!
  \brief Change the double_probes_

  \param double_probes 5 double* that constitute a probe
*/
void Individual::set_double_probes(double* double_probes) {
  double_probes_ = double_probes;
}

// =====================================================================
//                       functions' definition
// =====================================================================

void Individual::reset_dist_to_target_by_segment(
    double* dist_to_target_by_segment) {
  if (dist_to_target_by_segment_ != NULL) delete[] dist_to_target_by_segment_;
  dist_to_target_by_segment_ = dist_to_target_by_segment;
}

void Individual::renew_dist_to_target_by_feature() {
  if (dist_to_target_by_feature_ != NULL) delete[] dist_to_target_by_feature_;
  dist_to_target_by_feature_ = new double[NB_FEATURES];
}


void Individual::renew_fitness_by_feature() {
  if (fitness_by_feature_ != NULL) delete[] fitness_by_feature_;
  fitness_by_feature_ = new double[NB_FEATURES];
}

void Individual::do_transcription_translation_folding() {
  if (transcribed_ == true && translated_ == true && folded_ == true) return;

#ifdef __TRACING__
  auto t1 = high_resolution_clock::now();
#endif

  do_transcription();

#ifdef __TRACING__
  auto t2 = high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	ae_logger::addLog(DNA_TO_RNA,duration);

  t1 = high_resolution_clock::now();
#endif

  do_translation();

#ifdef __TRACING__
  t2 = high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	ae_logger::addLog(RNA_TO_PROTEIN,duration);

  t1 = high_resolution_clock::now();
#endif

  do_folding();

#ifdef __TRACING__
  t2 = high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	ae_logger::addLog(COMPUTE_PROTEIN,duration);
#endif

  make_protein_list();
}

#ifdef DEBUG

void Individual::assert_promoters() {
  // Perform assertion for each genetic unit
  // for (auto& gen_unit: genetic_unit_list_)
  //   gen_unit.assert_promoters();
}

void Individual::assert_promoters_order() {
  // Perform assertion for each genetic unit
  // for (auto& gen_unit: genetic_unit_list_)
  //   gen_unit.assert_promoters_order();
}

#endif

// =================================================================
//                            Public Methods
// =================================================================
void Individual::compute_phenotype() {
  if (phenotype_computed_) return; // Phenotype has already been computed, nothing to do.
  phenotype_computed_ = true;

  // Make sure the transcription, translation and folding stages have taken place
  do_transcription_translation_folding();


  // We will use two fuzzy sets :
  //   * phenotype_activ_ for the proteins realising a set of functions
  //   * phenotype_inhib_ for the proteins inhibiting a set of functions
  // The phenotype will then be given by the sum of these 2 fuzzy sets
  phenotype_activ_ = FuzzyFactory::fuzzyFactory->create_fuzzy();
  phenotype_inhib_ = FuzzyFactory::fuzzyFactory->create_fuzzy();

  for (const auto& gen_unit: genetic_unit_list_) {
    phenotype_activ_->add(*gen_unit.activ_contribution());
    phenotype_inhib_->add(*gen_unit.inhib_contribution());
  }

        phenotype_activ_->clip(AbstractFuzzy::max,   Y_MAX);
  phenotype_inhib_->clip(AbstractFuzzy::min, - Y_MAX);

  phenotype_ = FuzzyFactory::fuzzyFactory->create_fuzzy();
  phenotype_->add(*phenotype_activ_);
  phenotype_->add(*phenotype_inhib_);

  phenotype_->clip(AbstractFuzzy::min, Y_MIN);

  phenotype_->simplify();
}

void Individual::compute_distance_to_target(const PhenotypicTarget& target) {
// Compute the areas between the phenotype and the target for each segment
// If the target is not segmented, the total area is computed
  if (distance_to_target_computed_) {
    return;
  } // distance_to_target_ has already been computed, nothing to do.

  distance_to_target_computed_ = true;

  if (not phenotype_computed_)
    compute_phenotype();

  // Compute the difference between the (whole) phenotype and the target
  AbstractFuzzy* delta = FuzzyFactory::fuzzyFactory->create_fuzzy(*phenotype_);
  bool verbose = false;
  // if ((id_ == 966) && AeTime::time()==447) {
  //     printf("Target %lf :: I %lf\n",target.fuzzy()->get_geometric_area(),delta->get_geometric_area());
  //     verbose = true;
  // }
  
  delta->sub(*(target.fuzzy()),verbose);
  // if (id_ == 966 && AeTime::time()==447) {
  //     printf("Delta %lf\n",delta->get_geometric_area());
  //     // delta->print();
  // }
  PhenotypicSegment ** segments = target.segments();
  delete [] dist_to_target_by_segment_;
  dist_to_target_by_segment_ = new double [target.nb_segments()];

  for (size_t i = 0 ; i < static_cast<size_t>(target.nb_segments()) ; i++) {
    dist_to_target_by_segment_[i] = 0;
  }

  // TODO : We should take into account that we compute the areas in order (from the leftmost segment, rightwards)
  //   => We shouldn't parse the whole list of points on the left of the segment we are considering (we have
  //      already been through them!)

  for (size_t i = 0; i < static_cast<size_t>(target.nb_segments()); i++) {

    dist_to_target_by_segment_[i] = delta->get_geometric_area(
      segments[i]->start, segments[i]->stop);

    dist_to_target_by_feature_[segments[i]->feature] += dist_to_target_by_segment_[i];
  }

  delete delta;
}

/*!
  Computation of a "proper" fitness value (one that increases when the individual is fitter)

  Computation of a "proper" fitness value (one that increases when the individual is fitter)
  The behaviour of this function depends on many parameters and most notably on whether it is
  a "composite" fitness or not, and on the selection scheme.
*/
void Individual::compute_fitness(const PhenotypicTarget& target) {
  if (fitness_computed_) return; // Fitness has already been computed, nothing to do.
  fitness_computed_ = true;

#ifdef NORMALIZED_FITNESS
  for (int8_t i = 0 ; i < NB_FEATURES ; i++) {
    if (envir->area_by_feature(i)==0.) {
      fitness_by_feature_[i] = 0.;
    }
    else {
      fitness_by_feature_[i] =  (envir->area_by_feature(i) - dist_to_target_by_feature_[i]) / envir->area_by_feature(i);
      if ((fitness_by_feature_[i] < 0.) && (i != METABOLISM)) // non-metabolic fitness can NOT be lower than zero (we do not want individual to secrete a negative quantity of public good)
      {
        fitness_by_feature_[i] = 0.;
      }
    }
  }

  if ((! placed_in_population_) || (! exp_m_->with_secretion())) {
    fitness_ = fitness_by_feature_[METABOLISM];
  }
  else {
    fitness_ =  fitness_by_feature_[METABOLISM] * (1 + exp_m_->secretion_contrib_to_fitness() * (grid_cell_->compound_amount() - exp_m_->secretion_cost() * fitness_by_feature_[SECRETION]));
  }

  if (exp_m_->selection_scheme() == FITNESS_PROPORTIONATE) // Then the exponential selection is integrated inside the fitness value
  {
    fitness_ = exp(-exp_m_->selection_pressure() * (1 - fitness_));
  }
#else
  for (int8_t i = 0; i < NB_FEATURES; i++) {
    if (i == SECRETION) {
      fitness_by_feature_[SECRETION] = exp(-exp_m_->selection_pressure() *
                                           dist_to_target_by_feature_[SECRETION])
                                       - exp(-exp_m_->selection_pressure() *
                                             target.area_by_feature(SECRETION));

      if (fitness_by_feature_[i] < 0) {
        fitness_by_feature_[i] = 0;
      }
    }
    else {
      fitness_by_feature_[i] = exp(
          -exp_m_->selection_pressure() * dist_to_target_by_feature_[i]);

    }
  }

  // Calculate combined, total fitness here!
  // Multiply the contribution of metabolism and the amount of compound in the
  // habitat
  if ((!placed_in_population_) || (!exp_m_->with_secretion())) {
    fitness_ = fitness_by_feature_[METABOLISM];
  }
  else {
    fitness_ = fitness_by_feature_[METABOLISM]
               * (1 + exp_m_->secretion_contrib_to_fitness() *
                      grid_cell()->compound_amount()
                  - exp_m_->secretion_cost() *
                    fitness_by_feature_[SECRETION]);
  }
#endif
}


void Individual::clear_everything_except_dna_and_promoters() {
  protein_list_.clear();

  evaluated_ = false;
  transcribed_ = false;
  translated_ = false;
  folded_ = false;
  phenotype_computed_ = false;
  distance_to_target_computed_ = false;
  fitness_computed_ = false;

  for (auto& gen_unit: genetic_unit_list_)
    gen_unit.reset_expression();

  if (phenotype_activ_ != NULL) {
    delete phenotype_activ_;
    phenotype_activ_ = NULL;
  }

  if (phenotype_inhib_ != NULL) {
    delete phenotype_inhib_;
    phenotype_inhib_ = NULL;
  }

  if (phenotype_ != NULL) {
    delete phenotype_;
    phenotype_ = NULL;
  }

  // Initialize all the fitness-related stuff
  if (dist_to_target_by_segment_ != NULL) {
    delete[] dist_to_target_by_segment_;
    dist_to_target_by_segment_ = NULL;
  }

  for (int8_t i = 0; i < NB_FEATURES; i++) {
    dist_to_target_by_feature_[i] = 0.0;
    fitness_by_feature_[i] = 0.0;
  }


  // For each RNA / individual / genetic_unit delete proteins it knows
  // Deleting the protein itself is made only once

  for (auto& gen_unit: genetic_unit_list_)
    gen_unit.clear_transcribed_proteins();

  // Clear RNA and proteins
  rna_list_.clear();
  protein_list_.clear();

  // Reset statistical data
  delete metrics_;
  metrics_ = nullptr;
  delete nc_metrics_;
  nc_metrics_ = nullptr;

  modularity_ = -1;
}

void Individual::Reevaluate() {
  // useful for post-treatment programs that replay mutations
  // on a single individual playing the role of the successive
  // ancestors

  clear_everything_except_dna_and_promoters();
  Evaluate();
}

void Individual::ReevaluateInContext(const Habitat& habitat) {
  // useful for post-treatment programs that replay mutations
  // on a single individual playing the role of the successive
  // ancestors

  clear_everything_except_dna_and_promoters();
  EvaluateInContext(habitat);
}


void Individual::add_GU(char*& sequence, int32_t length) {
  clear_everything_except_dna_and_promoters();
  genetic_unit_list_.emplace_back(this, sequence, length);
  //genetic_unit_list_.back().set_local_id(cpt_local_id++);
}

/// Overloaded version to prevent the use of GeneticUnit disabled
/// copy ctor. Forwards arguments to GeneticUnit's ctor.
void Individual::add_GU(Individual* indiv,
                        int32_t chromosome_length,
                        std::shared_ptr<JumpingMT> prng) {
  clear_everything_except_dna_and_promoters();
  genetic_unit_list_.emplace_back(indiv, chromosome_length, prng);
  //genetic_unit_list_.back().set_local_id(cpt_local_id++);
}

void Individual::remove_GU(int16_t num_unit) {
  clear_everything_except_dna_and_promoters();
  auto it = genetic_unit_list_.begin();
  std::advance(it, num_unit);
  genetic_unit_list_.erase(it);
}


void Individual::do_transcription() {
  if (transcribed_) {
    return;
  } // Transcription has already been performed, nothing to do.
  transcribed_ = true;

  for (auto& gen_unit: genetic_unit_list_) {
    gen_unit.do_transcription();
    const auto& rna_list = gen_unit.rna_list();
    for (auto& strand: {LEADING, LAGGING}) {
      for (auto& rna: rna_list[strand])
        rna_list_.push_back(&rna);
    }
  }
}

void Individual::do_translation() {
  if (translated_) {
    return;
  } // ARNs have already been translated, nothing to do.
  translated_ = true;

  if (not transcribed_)
    do_transcription();

  for (auto& gen_unit: genetic_unit_list_) {
    gen_unit.do_translation();
    // append all proteins from `gen_unit` to `protein_list_`
    for (auto& strand_id: {LEADING, LAGGING}) {
      auto& strand = gen_unit.protein_list(strand_id);
      for (auto& p: strand)
        protein_list_.push_back(&p);
    }
  }
}

void Individual::do_folding() {
  if (folded_) {
    return;
  } // Proteins have already been folded, nothing to do.
  folded_ = true;

  if (not translated_)
    do_translation();

  for (auto& gen_unit: genetic_unit_list_)
    gen_unit.compute_phenotypic_contribution();
}

void Individual::Evaluate(bool no_signal) {
  EvaluateInContext(grid_cell_->habitat());
}

void Individual::EvaluateInContext(const Habitat& habitat) {
  if (evaluated_ == true) {
    return;
  } // Individual has already been evaluated, nothing to do.
  evaluated_ = true;

  // ----------------------------------------------------------------------
  // Transcription - Translation - Folding
  // ----------------------------------------------------------------------
  do_transcription_translation_folding();

  // ----------------------------------------------------------------------
  // Compute phenotype and compare it to the target => fitness
  // ----------------------------------------------------------------------
#ifdef __TRACING__
  auto t1 = high_resolution_clock::now();
#endif
  compute_phenotype();
  compute_distance_to_target(habitat.phenotypic_target());
  compute_fitness(habitat.phenotypic_target());
#ifdef __TRACING__
  auto t2 = high_resolution_clock::now();
	  	  auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	  	  ae_logger::addLog(PROTEIN_TO_PHENOTYPE,duration);
#endif

  if (exp_m_->output_m()->compute_phen_contrib_by_GU())
    for (auto& gen_unit: genetic_unit_list_) {
      gen_unit.compute_distance_to_target(habitat.phenotypic_target());
      gen_unit.compute_fitness(habitat.phenotypic_target());
    }
}


void Individual::inject_GU(Individual* donor) {
  // Add the GU at the end of the list
  genetic_unit_list_.emplace_back(this, donor->genetic_unit_list_.back());
}

void Individual::inject_2GUs(Individual* partner) {
  // We swap GUs from the end of the list.

  // TODO vld: As far as I understood the old code (47b27578), the
  // elements were not swapped but appended to the end of the other GU
  // lists. Error?  Original author (Dule, commit 47b27578), asked for
  // clarification by e-mail on 2015-02-23.

  const auto& gu_list_back_it = std::prev(
      genetic_unit_list_.end()); // initial last cell from genetic_unit_list_
  genetic_unit_list_.splice(genetic_unit_list_.end(),
                            partner->genetic_unit_list_,
                            std::prev(partner->genetic_unit_list_.end()));
  partner->genetic_unit_list_.splice(partner->genetic_unit_list_.end(),
                                     genetic_unit_list_, gu_list_back_it);
}

void Individual::compute_statistical_data() {
  if (metrics_ != nullptr) {
    return;
  } // Statistical data has already been computed,
  // nothing to do.

  metrics_ = new Metrics();

  if (not phenotype_computed_)
    compute_phenotype();

  for (const auto& gen_unit : genetic_unit_list_) {
    metrics_->Accumulate(gen_unit);
  }
}

void Individual::compute_non_coding() {
  if (nc_metrics_ != nullptr) return; // NC stats have already been computed,
  // nothing to do.
  nc_metrics_ = new NonCodingMetrics();

  for (auto& gen_unit: genetic_unit_list_) {
    gen_unit.compute_non_coding();
    nc_metrics_->Accumulate(gen_unit);
  }
}

void Individual::save(gzFile backup_file) const {
  //printf("Appel à la sauvegarde de Individual\n");
  // Write the name and "age" of the strain
  int8_t strain_string_len = strlen(strain_name_);
  gzwrite(backup_file, &strain_string_len, sizeof(strain_string_len));
  gzwrite(backup_file, strain_name_, strain_string_len + 1);
  gzwrite(backup_file, &age_, sizeof(age_));

  // Write id and rank
  gzwrite(backup_file, &id_, sizeof(id_));
  gzwrite(backup_file, &long_id_, sizeof(long_id_));
  gzwrite(backup_file, &rank_, sizeof(rank_));

  // Write the position of the individual
  // gzwrite(backup_file, &x, sizeof(x));
  // gzwrite(backup_file, &y, sizeof(y));

  // Write generic probes
  gzwrite(backup_file, int_probes_, 5 * sizeof(*int_probes_));
  gzwrite(backup_file, double_probes_, 5 * sizeof(*double_probes_));

  // Write mutational parameters
  mut_params_->save(backup_file);

  // ------------------------------------------------- Phenotypic stochasticity
  gzwrite(backup_file, &with_stochasticity_, sizeof(with_stochasticity_));

  // Write artificial chemistry parameters
  gzwrite(backup_file, &w_max_, sizeof(w_max_));

  // Write genome size constraints
  gzwrite(backup_file, &min_genome_length_, sizeof(min_genome_length_));
  gzwrite(backup_file, &max_genome_length_, sizeof(max_genome_length_));

  // Write plasmids settings
  int8_t tmp_allow_plasmids = allow_plasmids_;
 gzwrite(backup_file, &tmp_allow_plasmids, sizeof(tmp_allow_plasmids));

  // Write genetic units
  int16_t nb_gen_units = genetic_unit_list_.size();
  gzwrite(backup_file, &nb_gen_units, sizeof(nb_gen_units));

  for (const auto& gen_unit: genetic_unit_list_)
    gen_unit.save(backup_file);
}

int32_t Individual::nb_terminators() {
  int32_t nb_term = 0;
  for (auto& gen_unit: genetic_unit_list_)
    nb_term += gen_unit.nb_terminators();
  return nb_term;
}


/*!
  \brief Compute reproduction statistics and statistics about the offsprings of the current individual

  * Make nb_children replications of the current individual.
  * For each replication, determine if the offsprings is neutral, beneficial or deleterious by comparison of fitness with the current individual (the parent)
  * If statistics about offsprings are required (offsprings_statistics != NULL), fitness mean, fitness variance, size mean, size variance, functional gene number mean, functional gene number variance fo the nb_children offsprings are computed
  * If information about each children are required (replication_file != NULL), fitness, genome_size, nb of functional genes, number of coding bases, number of transcribed but not translated bases, number of non transcribed bases of each offsprings are written in replication_file

  \param nb_children              number of replications made to have the statistics
  \param reproduction_statistics  statistics about the replications (proportion of neutral offsprings, proportion of beneficial offsprings, proportion of deleterious offsprings)
  \param offsprings_statistics    statistics about the nb_children offsprings (fitness mean, fitness variance, size mean, size variance, functional gene number mean,
                                    functional gene number variance) compute if not null
  \param replication_file         file with information about each children of the current individual (fitness, genome_size, nb of functional genes, number of coding bases,
                                    number of transcribed but not translated bases, number of non transcribed bases) if not null
*/
void Individual::compute_experimental_f_nu(int32_t nb_children,
                                           double* reproduction_statistics,
                                           double* offsprings_statistics,
                                           FILE* replication_file) {
  double initial_fitness = fitness();

  if (reproduction_statistics != NULL) {
    reproduction_statistics[0] = 0; // proportion of neutral offsprings
    reproduction_statistics[1] = 0; // proportion of beneficial offsprings
    reproduction_statistics[2] = 0; // proportion of deleterious offsprings
  }
  else {
    printf("%s:%d: error: reproduction_statistics was not initialized\n",
           __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  if (offsprings_statistics != NULL) {
    offsprings_statistics[0] = 0; // offspring fitness mean
    offsprings_statistics[1] = 0; // offspring fitness variance
    offsprings_statistics[2] = 0; // offspring size mean
    offsprings_statistics[3] = 0; // offspring size variance
    offsprings_statistics[4] = 0; // offspring functional gene number mean
    offsprings_statistics[5] = 0; // offspring functional gene number variance
  }

  // ------------------------------------------
  //      Simulate fitness degradation
  // ------------------------------------------

  double fitness_child = 0.0;
  double metabolic_error_child = 0.0;

  // replicate this individual to create 'nb_children' children
  Individual* child = NULL;

  int32_t genome_size = 0;
  int32_t nb_functional_genes = 0;
  int32_t nb_bases_in_0_functional_CDS = 0;
  int32_t nb_bases_in_0_coding_RNA = 0;

  for (int i = 0; i < nb_children; i++) {
    int8_t tm;
    child = exp_m_->exp_s()->sel()->do_replication(this, id_,tm);
    fitness_child = child->fitness();
    metabolic_error_child = child->dist_to_target_by_feature(METABOLISM);

    if (fabs(initial_fitness - fitness_child) <
        1e-10 * std::max(initial_fitness, fitness_child)) {
      reproduction_statistics[0] += 1;
    }
    else if (fitness_child > initial_fitness) {
      reproduction_statistics[1] += 1;
    }
    else {
      reproduction_statistics[2] += 1;
    }

    genome_size = child->total_genome_size();
    nb_functional_genes = child->nb_functional_genes();
    nb_bases_in_0_functional_CDS = child->nb_bases_in_0_functional_CDS();
    nb_bases_in_0_coding_RNA = child->nb_bases_in_0_coding_RNA();

    if (offsprings_statistics != NULL) {
      offsprings_statistics[0] += fitness_child;
      offsprings_statistics[1] += pow(fitness_child, 2);
      offsprings_statistics[2] += (double) genome_size;
      offsprings_statistics[3] += pow((double) genome_size, 2);
      offsprings_statistics[4] += (double) nb_functional_genes;
      offsprings_statistics[5] += pow((double) nb_functional_genes, 2);
    }

    if (replication_file != NULL) {
      fprintf(replication_file,
              "%le %le %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 "\n",
              fitness_child, metabolic_error_child, genome_size,
              nb_functional_genes, genome_size - nb_bases_in_0_functional_CDS,
              nb_bases_in_0_functional_CDS - nb_bases_in_0_coding_RNA,
              nb_bases_in_0_coding_RNA);
    }

    delete child;
  }

  //compute Fv
  reproduction_statistics[0] /= (double) nb_children;
  reproduction_statistics[1] /= (double) nb_children;
  reproduction_statistics[2] /= (double) nb_children;

  if (offsprings_statistics != NULL) {
    offsprings_statistics[0] /= (double) nb_children;
    offsprings_statistics[1] /= (double) nb_children;
    offsprings_statistics[2] /= (double) nb_children;
    offsprings_statistics[3] /= (double) nb_children;
    offsprings_statistics[4] /= (double) nb_children;
    offsprings_statistics[5] /= (double) nb_children;

    offsprings_statistics[1] -= pow(offsprings_statistics[0], 2);
    offsprings_statistics[3] -= pow(offsprings_statistics[2], 2);
    offsprings_statistics[5] -= pow(offsprings_statistics[4], 2);
  }
}


/// Compute reproduction theoretical proportion of neutral offsprings.
///
/// Compute the theoretical proportion of neutral offsprings given the
/// Carole's formula, based on the mutations and rearrangement rates
/// and not on multiple replications.
///
/// \return theoretical proportion of neutral offsprings
double Individual::compute_theoritical_f_nu() {
  // We first have to collect information about genome structure.
  // Abbreviations are chosen according to Carole's formula.
  // Please notice that compared to the formula we have the beginning
  // and ends of neutral regions instead of 'functional regions'
  GeneticUnit& chromosome = genetic_unit_list_.front();
  int32_t L = chromosome.dna()->length();
  int32_t N_G = chromosome.nb_neutral_regions(); // which is not exactly Carole's original definition
  int32_t* b_i = chromosome.beginning_neutral_regions();
  int32_t* e_i = chromosome.end_neutral_regions();
  int32_t lambda = chromosome.nb_bases_in_neutral_regions();
  int32_t l = L - lambda; // nb bases in 'functional regions'

  int32_t* lambda_i = NULL;  // nb bases in ith neutral region
  if (N_G > 0) // all the chromosome may be functional
  {
    lambda_i = new int32_t[N_G];

    for (int32_t i = 0; i < N_G - 1; i++) {
      lambda_i[i] = e_i[i] - b_i[i] + 1;
    }
    if (b_i[N_G - 1] > e_i[N_G -
                           1]) // last neutral region is overlapping on the beginning of chromosome
    {
      lambda_i[N_G - 1] = (e_i[N_G - 1] + L) - b_i[N_G - 1] + 1;
    }
    else // no overlap
    {
      lambda_i[N_G - 1] = e_i[N_G - 1] - b_i[N_G - 1] + 1;
    }
  }

  // we now compute the probabilities of neutral reproduction for
  // each type of mutation and rearrangement and update Fv
  double Fv = 1;

  // mutation + insertion + deletion
  double nu_local_mutation = 1 - ((double) l) / L;
  Fv = pow(1 - point_mutation_rate() * (1 - nu_local_mutation), L);
  Fv *= pow(1 - small_insertion_rate() * (1 - nu_local_mutation), L);
  Fv *= pow(1 - small_deletion_rate() * (1 - nu_local_mutation), L);

  // inversion ~ two local mutations
  double nu_inversion = nu_local_mutation * nu_local_mutation;
  Fv *= pow(1 - inversion_rate() * (1 - nu_inversion), L);

  // translocation ~ inversion + insertion (mathematically)
  Fv *= pow(
      1 - translocation_rate() * (1 - nu_inversion * nu_local_mutation), L);

  // long deletion
  double nu_deletion = 0; // if N_G == 0, a deletion is always not neutral
  for (int32_t i = 0; i < N_G; i++) {
    nu_deletion += lambda_i[i] * (lambda_i[i] + 1);
  }
  nu_deletion /= ((double) 2 * L * L);
  Fv *= pow(1 - deletion_rate() * (1 - nu_deletion), L);

  // duplication ~ big deletion + insertion
  Fv *= pow(1 - duplication_rate() * (1 - nu_deletion * nu_local_mutation),
            L);

  if (lambda_i != NULL) delete[] lambda_i;

  return Fv;
}

/// Remove the bases that are not in coding RNA.
///
/// Remove the bases that are not in coding RNA and test at each loss
/// that fitness is not changed.
void Individual::remove_non_coding_bases() {
  for (auto& gen_unit: genetic_unit_list_)
    gen_unit.remove_non_coding_bases();

  // Delete the obsolete stats
  delete metrics_;
  metrics_ = NULL;
  delete nc_metrics_;
  nc_metrics_ = NULL;

#ifdef DEBUG
  compute_statistical_data();
  compute_non_coding();
  assert(nb_bases_in_0_coding_RNA() == 0);
#endif
}

/// Double the bases that are not in coding RNA.
///
/// Double the bases that are not in coding RNA by addition of random
/// bases and test at each addition that fitness is not changed.
void Individual::double_non_coding_bases() {
  metrics_->total_genome_size_ = 0;
  int32_t initial_non_coding_base_nb = nb_bases_in_0_coding_RNA();

  for (auto& gen_unit: genetic_unit_list_)
    gen_unit.double_non_coding_bases();

  // Delete the obsolete stats
  delete metrics_;
  metrics_ = NULL;
  delete nc_metrics_;
  nc_metrics_ = NULL;

#ifdef DEBUG
  compute_statistical_data();
  compute_non_coding();
  assert(nb_bases_in_0_coding_RNA() == 2 * initial_non_coding_base_nb);
#endif
}
/*
void Individual::copy_parent(const Individual* parent, bool env_will_change) {
  // copy fitness related data
  for (int8_t i = 0; i < NB_FEATURES; i++) {
    dist_to_target_by_feature_[i] = parent->dist_to_target_by_feature_[i];
    fitness_by_feature_[i] = parent->fitness_by_feature_[i];
  }

  // Copy genetic unit
  for (auto& gen_unit: parent->genetic_unit_list_) {
    for (auto& gunit: genetic_unit_list_) {
      if (gen_unit.local_id() == gunit.local_id()) {
        gunit->copy_parent(gen_unit,env_will_change);
      }
    }
  }

  make_protein_list();
  make_rna_list();


  transcribed_ = true;
  translated_  = true;
  folded_  = true;

  if (!env_will_change) {
    phenotype_computed_ = true;
    distance_to_target_computed_ = true;
    fitness_computed_ = true;
    evaluated_ = true;
  } else {
    phenotype_computed_ = false;
    distance_to_target_computed_ = false;
    fitness_computed_ = false;
    evaluated_ = false;
  }
}*/

// =================================================================
//                           Protected Methods
// =================================================================

// TODO vld: refactor make_protein_list and make_rna_list
void Individual::make_protein_list() {
  // Clear list
  protein_list_.clear();

  // Make a copy of each genetic unit's protein list
  for (auto& gen_unit: genetic_unit_list_) {
    // append all proteins from `gen_unit` to `protein_list_`
    for (auto& strand_id: {LEADING, LAGGING}) {
      auto& strand = gen_unit.protein_list(strand_id);
      for (auto& p: strand)
        protein_list_.push_back(&p);
    }
  }
}

void Individual::make_rna_list() {
  // Clear list
  rna_list_.clear();

  // Make a copy of each genetic unit's rna list
  for (const auto& gen_unit: genetic_unit_list_) {
    // Create proxies
    const auto& rna_list = gen_unit.rna_list();

    // append pointers to rna material to local rna_list_
    for (auto& strand: {LEADING, LAGGING})
      for (auto& rna: rna_list[strand])
        rna_list_.push_back(&rna);
  }
}
} // namespace aevol
