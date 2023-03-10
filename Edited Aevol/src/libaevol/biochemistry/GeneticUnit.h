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

#ifndef AEVOL_GENETIC_UNIT_H_
#define AEVOL_GENETIC_UNIT_H_

#include <cinttypes>
#include <cstdio>
#include <memory>
#include <array>

#include "Dna.h"
#include "Rna.h"
#include "Protein.h"
#include "AbstractFuzzy.h"
#include "FuzzyFactory.h"
#include "JumpingMT.h"
#include "Utils.h"
#include "PhenotypicTarget.h"

#ifdef __REGUL
#include "raevol/Protein_R.h"
#endif

using std::vector;
using std::list;

namespace aevol {

#ifndef __REGUL
// TODO vld: check types
using Promoters1Strand = std::list<Rna>;
#else
using Promoters1Strand = std::list<Rna_R>;
#endif

using Promoters2Strands = std::vector<Promoters1Strand>;

class ExpManager;
class MutationEvent;

class GeneticUnit {
 public :

  // =================================================================
  //                             Constructors
  // =================================================================

  GeneticUnit(Individual* indiv);

  GeneticUnit(Individual* indiv,
              int32_t length,
              std::shared_ptr<JumpingMT> prng);

  GeneticUnit(Individual* indiv,
              char* seq,
              int32_t length,
              const Promoters2Strands& prom_list = {{},
                                                    {}});

  GeneticUnit(Individual* indiv, const GeneticUnit& model);

  GeneticUnit(Individual* indiv, const GeneticUnit* parent);

  GeneticUnit(Individual* indiv, gzFile backup_file);

  GeneticUnit(Individual* indiv, char* organism_file_name);

  GeneticUnit() = delete;

  GeneticUnit(const GeneticUnit&) = delete;

  // =================================================================
  //                             Destructors
  // =================================================================
  virtual ~GeneticUnit();

  // =================================================================
  //                              Accessors
  // =================================================================

  ExpManager* exp_m() const;

  Individual* indiv() const;

  Dna* dna() const;

  AbstractFuzzy* activ_contribution() const;

  AbstractFuzzy* inhib_contribution() const;

  AbstractFuzzy* phenotypic_contribution() const;

  const Promoters2Strands& rna_list() const;

  // TODO return as reference
#ifndef __REGUL
  std::list<Protein>& protein_list(Strand strand);
#else
  std::list<Protein_R>& protein_list(Strand strand);
  #endif
  void clear_protein_list(Strand strand);
  void clear_transcribed_proteins();


  // Direct DNA access
  const char* sequence() const;

  int32_t seq_length() const;

  // Statistical data
  int32_t nb_coding_RNAs() const;

  int32_t nb_non_coding_RNAs() const;

  double overall_size_coding_RNAs() const;

  double av_size_coding_RNAs() const;

  double overall_size_non_coding_RNAs() const;

  double av_size_non_coding_RNAs() const;

  int32_t nb_genes_activ() const;

  int32_t nb_genes_inhib() const;

  int32_t nb_functional_genes() const;

  int32_t nb_non_functional_genes() const;

  double overall_size_functional_genes() const;

  double av_size_functional_genes() const;

  double overall_size_non_functional_genes() const;

  double av_size_non_functional_genes() const;

  int32_t nb_bases_in_0_CDS() const;

  int32_t nb_bases_in_0_functional_CDS() const;

  int32_t nb_bases_in_0_non_functional_CDS() const;

  int32_t nb_bases_in_0_RNA() const;

  int32_t nb_bases_in_0_coding_RNA() const;

  int32_t nb_bases_in_0_non_coding_RNA() const;

  int32_t nb_bases_non_essential() const;

  int32_t nb_bases_non_essential_including_nf_genes() const;

  int32_t nb_bases_in_neutral_regions() const;

  int32_t nb_neutral_regions() const;

  int32_t* beginning_neutral_regions() const;

  int32_t* end_neutral_regions() const;

  double modularity() const;

  double dist_to_target_by_feature(PhenotypicFeature feature) const;

  double fitness() const;

  double fitness_by_feature(PhenotypicFeature feature) const;

  int32_t min_gu_length() const;

  int32_t max_gu_length() const;

  void set_min_gu_length(int32_t min_gu_length);

  void set_max_gu_length(int32_t max_gu_length);

  void set_exp_m(ExpManager* exp_m);

  int local_id() { return local_id_;}

  void set_local_id(int local_id) { local_id_ = local_id;}

  // =================================================================
  //                            Public Methods
  // =================================================================
  void locate_promoters();

  void do_transcription();

  void do_translation();

  void compute_phenotypic_contribution(int indiv_id = -1);

  void take_ownership_of_all_rnas() { Dna::set_GU(rna_list(), this); };


  // DM: these two are identical to functions from Individual
  void compute_distance_to_target(const PhenotypicTarget& target);

  void compute_fitness(const PhenotypicTarget& target);

  void reset_expression(); // useful for post-treatment programs

  void print_rnas() const;

  void print_coding_rnas();

  static void print_rnas(const Promoters2Strands& rnas);

  static void print_rnas(const Promoters1Strand& rnas, Strand strand);

  void print_proteins() const;

  bool is_promoter(Strand strand, int32_t pos, int8_t& dist) const;
  //bool is_promoterX(Strand strand, int32_t pos, int8_t& dist) const;

  bool is_terminator(Strand strand, int32_t pos) const;

  bool is_shine_dalgarno(Strand strand, int32_t pos) const;

  bool is_start(Strand strand, int32_t pos) const;

  bool is_stop(Strand strand, int32_t pos) const;

  #ifdef BASE_2
  int8_t codon(Strand strand, int32_t pos) const;
  #endif
  
  void compute_non_coding();


  // these functions are called once, they should likely not be public methods
  void duplicate_promoters_included_in(int32_t pos_1,
                                       int32_t pos_2,
                                       Promoters2Strands& duplicated_promoters);

  void promoters_included_in(int32_t pos_1,
                                 int32_t pos_2,
                                 Promoters2Strands& promoters_list);

  void promoters(Strand strand_id,
                     Position start,
                     int32_t pos1,
                     int32_t pos2,
                     Promoters1Strand& promoters);

  void invert_promoters_included_in(int32_t pos_1, int32_t pos_2);

  static void invert_promoters(Promoters2Strands& promoter_lists,
                               int32_t seq_length);

  static void invert_promoters(Promoters2Strands& promoter_lists,
                               int32_t pos_1, int32_t pos_2); // See WARNING

  void extract_promoters_included_in(int32_t pos_1,
                                     int32_t pos_2,
                                     Promoters2Strands& extracted_promoters);

  void extract_promoters_starting_between(
      int32_t pos_1, int32_t pos_2, Promoters2Strands& extracted_promoters);

  void extract_leading_promoters_starting_between(
      int32_t pos_1, int32_t pos_2, Promoters1Strand& extracted_promoters);

  void extract_lagging_promoters_starting_between(
      int32_t pos_1, int32_t pos_2, Promoters1Strand& extracted_promoters);





  // end comment

  static void shift_promoters(Promoters2Strands& promoters_to_shift,
                              int32_t delta_pos,
                              int32_t seq_length);

  void insert_promoters(Promoters2Strands& promoters_to_insert);

  void insert_promoters_at(Promoters2Strands& promoters_to_insert,
                           int32_t pos);

  void remove_promoters_around(int32_t pos);

  void remove_promoters_around(int32_t pos_1, int32_t pos_2);

  void remove_all_promoters();

  void look_for_new_promoters_around(int32_t pos);

  void look_for_new_promoters_around(int32_t pos_1, int32_t pos_2);

  void move_all_promoters_after(int32_t pos, int32_t delta_pos);

  void move_all_leading_promoters_after(int32_t pos, int32_t delta_pos);

  void move_all_lagging_promoters_after(int32_t pos, int32_t delta_pos);

  void copy_promoters_included_in(
      int32_t pos_1, int32_t pos_2, Promoters2Strands& new_promoter_lists);

  void copy_promoters_starting_between(
      int32_t pos_1, int32_t pos_2, Promoters2Strands& new_promoter_lists);

  void copy_leading_promoters_starting_between(
      int32_t pos_1, int32_t pos_2, Promoters1Strand& new_promoter_list);

  void copy_lagging_promoters_starting_between(
      int32_t pos_1, int32_t pos_2, Promoters1Strand& new_promoter_list);

  void save(gzFile backup_file) const;

  int32_t nb_terminators();

  //~ // set the genetic unit of all promoters to the appropriate
  //~ void reasign_promoter_genetic_unit ();

  void assert_promoters();

  void assert_promoters_order();

  bool* is_belonging_to_coding_RNA();

  void remove_non_coding_bases();

  void double_non_coding_bases();

  // WARNING: The method below works properly only in the case of a single
  // genetic unit (no plasmid).
  // Translocations between different genetic units is not managed.
  void compute_nb_of_affected_genes(const Mutation* mut,
                                    int& nb_genes_at_breakpoints,
                                    int& nb_genes_in_segment,
                                    int& nb_genes_in_replaced_segment);

  // =================================================================
  //                           Public Attributes
  // =================================================================

#ifdef __REGUL
    Individual_R* indiv_r_;
#endif


        bool phenotypic_contributions_computed_;
  Promoters2Strands rna_list_ = {{},
                                 {}};
 protected :
  // =================================================================
  //                           Protected Methods
  // =================================================================
  void init_statistical_data();

  void remove_leading_promoters_starting_between(int32_t pos_1, int32_t pos_2);

  void remove_lagging_promoters_starting_between(int32_t pos_1, int32_t pos_2);

  void remove_leading_promoters_starting_before(int32_t pos);

  void remove_lagging_promoters_starting_before(int32_t pos);

  void remove_leading_promoters_starting_after(int32_t pos);

  void remove_lagging_promoters_starting_after(int32_t pos);

  void look_for_new_leading_promoters_starting_between(int32_t pos_1,
                                                       int32_t pos_2);

  void look_for_new_lagging_promoters_starting_between(int32_t pos_1,
                                                       int32_t pos_2);

  void look_for_new_leading_promoters_starting_after(int32_t pos);

  void look_for_new_lagging_promoters_starting_after(int32_t pos);

  void look_for_new_leading_promoters_starting_before(int32_t pos);

  void look_for_new_lagging_promoters_starting_before(int32_t pos);


  // =================================================================
  //                          Protected Attributes
  // =================================================================

  ExpManager* exp_m_ = nullptr;

  Individual* indiv_ = nullptr;
  Dna* dna_ = nullptr;
  AbstractFuzzy* activ_contribution_ = nullptr;
  AbstractFuzzy* inhib_contribution_ = nullptr;
  // NB : phenotypic_contribution_ is only an indicative value,
  // not used for the whole phenotype computation
  AbstractFuzzy* phenotypic_contribution_ = nullptr;

  // rna_list_ always has 2 elements: make it an std::array

  // protein_list_ always has 2 elements: make it an std::array
#ifndef __REGUL
  std::array<std::list<Protein>, 2> protein_list_; // = {{},{}};
#else
  std::array<std::list<Protein_R>, 2> protein_list_; // = {{},{}};
  #endif

  // DM: For plasmid work, we sometimes *need* all the data
  // (e.g. fitness, secretion) calculated for each GU
  double* dist_to_target_per_segment_ = nullptr;
  double* dist_to_target_by_feature_ = nullptr;
  double* fitness_by_feature_ = nullptr;
  double fitness_;

  // Statistical data (intrinsic for this genetic unit)
  /**
   * Number of coding RNAs (at least one gene on RNA)
   */
  int32_t nb_coding_RNAs_;

  /**
   * Number of non-coding-RNAs
   */
  int32_t nb_non_coding_RNAs_;

  /**
   * Average size of coding RNAs
   */
  double overall_size_coding_RNAs_;

  /**
   * Average size of non-coding RNAs
   */
  double overall_size_non_coding_RNAs_;

  /**
   * Number of genes realizing a function
   */
  int32_t nb_genes_activ_;

  /**
   * Number of genes inhibitting a function
   */
  int32_t nb_genes_inhib_;

  /**
   * Number of functional genes
   */
  int32_t nb_fun_genes_;

  /**
   * Number of non-functional genes
   */
  int32_t nb_non_fun_genes_;

  /**
   * Average size of functional genes
   */
  double overall_size_fun_genes_;

  /**
   * Average size of non-functional genes
   */
  double overall_size_non_fun_genes_;

  /**
   * Number of bases that are not included in any gene
   */
  int32_t nb_bases_in_0_CDS_;

  /**
   * Number of bases that are not included in any non-degenerated gene
   */
  int32_t nb_bases_in_0_functional_CDS_;

  /**
   * Number of bases that are not included in any degenerated gene
   */
  int32_t nb_bases_in_0_non_functional_CDS_;

  /**
   * Number of bases that are not included in any RNA
   */
  int32_t nb_bases_in_0_RNA_;

  /**
   * Number of bases that are not included in any coding RNA
   * (RNAs containing at least one CDS)
   */
  int32_t nb_bases_in_0_coding_RNA_;

  /**
   * Number of bases that are not included in any non coding RNA
   */
  int32_t nb_bases_in_0_non_coding_RNA_;

  /**
   * Number of bases that are not included in any non-degenerated gene or
   * involved in its expression
   */
  int32_t nb_bases_non_essential_;

  /**
   * Number of bases that are not included in any gene or involved in its
   * expression
   */
  int32_t nb_bases_non_essential_including_nf_genes_;

  /**
   * Number of bases that are not included in a neutral region
   */
  int32_t nb_bases_in_neutral_regions_;

  /**
   * A base is considered neutral when neither itself NOR its corresponding
   * base on the other strand belongs to a coding [promoter->terminator] region
   */
  /**
   * Number of neutral regions
   */
  int32_t nb_neutral_regions_;

  /**
   * Beginning of neutral regions
   */
  int32_t* beginning_neutral_regions_ = nullptr;

  /**
   * Corresponding ends of neutral regions
   */
  int32_t* end_neutral_regions_ = nullptr;

  /**
   * Ratio between the pairwise distance between genes whose corresponding
   * phenotypic triangles overlap and the average intergenic distance
   * (ignoring non-functional genes)
   */
  double modularity_;

  int32_t min_gu_length_;
  int32_t max_gu_length_;

  // We keep trace of what we have already computed to avoid double
  // computation (mainly in post-treaments)
  bool transcribed_;
  bool translated_;
  bool non_coding_computed_;
  bool distance_to_target_computed_;
  bool fitness_computed_;


  int local_id_;
  int cpt_local_id=0;
};

} // namespace aevol
#endif // AEVOL_GENETIC_UNIT_H_
