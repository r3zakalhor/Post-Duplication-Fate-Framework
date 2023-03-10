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


#ifndef AEVOL_STAT_RECORD_H_
#define AEVOL_STAT_RECORD_H_


// =================================================================
//                              Libraries
// =================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <list>



// =================================================================
//                            Project Files
// =================================================================


namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================
class ExpManager;
class Stats;
class Individual;
class ae_population;
class ExpSetup;
class ReplicationReport;

enum indiv_or_pop
{
  INDIV           = 0,
  POP             = 1,
  STDEVS          = 2,
  SKEWNESS        = 3,
  NB_INDIV_OR_POP = 4
};

enum chrom_or_gen_unit
{
  ALL_GU          = 0,
  CHROM           = 1,
  PLASMIDS        = 2,
  NB_CHROM_OR_GU  = 3
};

enum best_or_glob
{
  BEST            = 0,
  GLOB            = 1,
  SDEV            = 2,
  SKEW            = 3,
  NB_BEST_OR_GLOB = 4
};

enum stats_type
{
  FITNESS_STATS   = 0,
  MUTATION_STATS  = 1,
  GENES_STATS     = 2,
  BP_STATS        = 3,
  REAR_STATS      = 4,
  NB_STATS_TYPES  = 5
};






class StatRecord
{
  friend class Stats;
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    StatRecord() = delete;
    StatRecord(ExpManager * exp_m);
    StatRecord(const StatRecord &model);
    StatRecord(ExpSetup* exp_s,
               Individual * indiv,
               ReplicationReport* replic_report,
               chrom_or_gen_unit chrom_or_gu = CHROM,
               bool compute_non_coding = true);
    StatRecord(ExpSetup* exp_s,
               std::list<std::pair<Individual*,
                       ReplicationReport*>> annotated_indivs,
               chrom_or_gen_unit chrom_or_gu = CHROM);
    StatRecord(ExpSetup* exp_s,
               std::list<std::pair<Individual*,
                       ReplicationReport*>> annotated_indivs,
               const StatRecord * means,
               chrom_or_gen_unit chrom_or_gu = CHROM);
    StatRecord(ExpSetup* exp_s,
               std::list<std::pair<Individual*,
                       ReplicationReport*>> annotated_indivs,
               const StatRecord * means,
               const StatRecord* stdevs,
               chrom_or_gen_unit chrom_or_gu = CHROM);

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~StatRecord();

    // =================================================================
    //                              Accessors
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    void initialize_data();
    void write_to_file(int64_t time, FILE* stat_file, stats_type stat_type_to_print) const;

    void divide(double divisor);
    void divide_record(StatRecord const * means, double power);

    void add(StatRecord * to_add, int32_t index);
    void substract_power(StatRecord const * means, StatRecord const * to_substract, double power);

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :
    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ExpManager * exp_m_;

    // NB : All the attributes are doubles because they will be used to
    //      compute averages over the population.
    indiv_or_pop record_type_;

    int32_t pop_size_ = 0.0;

    double  fitness_ = 0.0;

    double  metabolic_error_ = 0.0;
    double  parent_metabolic_error_ = 0.0;
    double  metabolic_fitness_ = 0.0;

    double  secretion_error_ = 0.0;
    double  parent_secretion_error_ = 0.0;
    double  secretion_fitness_ = 0.0;

    double  compound_amount_ = 0.0;

    int32_t  amount_of_dna_ = 0.0;
    int32_t  nb_coding_rnas_ = 0.0;
    int32_t  nb_non_coding_rnas_ = 0.0;

    // NOT including promoter but including terminator
    double  av_size_coding_rnas_ = 0.0;
    double  av_size_non_coding_rnas_ = 0.0;

    int32_t  nb_functional_genes_ = 0.0;
    int32_t  nb_non_functional_genes_ = 0.0;

    // NOT including START and STOP codons
    double  av_size_functional_gene_ = 0.0;
    double  av_size_non_functional_gene_ = 0.0;

    int32_t  nb_mut_ = 0.0;
    int32_t  nb_rear_ = 0.0;
    int32_t  nb_switch_ = 0.0;
    int32_t  nb_indels_ = 0.0;
    int32_t  nb_dupl_ = 0.0;
    int32_t  nb_del_ = 0.0;
    int32_t  nb_trans_ = 0.0;
    int32_t  nb_inv_ = 0.0;

    double  dupl_rate_ = 0.0;
    double  del_rate_ = 0.0;
    double  trans_rate_ = 0.0;
    double  inv_rate_ = 0.0;
    double  mean_align_score_ = 0.0;

    int32_t  nb_bases_in_0_CDS_ = 0.0;
    int32_t  nb_bases_in_0_functional_CDS_ = 0.0;
    int32_t  nb_bases_in_0_non_functional_CDS_ = 0.0;
    int32_t  nb_bases_in_0_RNA_ = 0.0;
    int32_t  nb_bases_in_0_coding_RNA_ = 0.0;
    int32_t  nb_bases_in_0_non_coding_RNA_ = 0.0;

    int32_t  nb_bases_non_essential_ = 0.0;
    int32_t  nb_bases_non_essential_including_nf_genes_ = 0.0;

#ifdef __REGUL
        int32_t  nb_influences_ = 0.0;
      int32_t  nb_enhancing_influences_ = 0.0;
      int32_t  nb_operating_influences_ = 0.0;
      double  av_value_influences_ = 0.0;
      double  av_value_enhancing_influences_ = 0.0;
      double  av_value_operating_influences_ = 0.0;
      int32_t  nb_TF_;
      int32_t  nb_pure_TF_;
#endif
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol

#endif // AEVOL_STAT_RECORD_H_
