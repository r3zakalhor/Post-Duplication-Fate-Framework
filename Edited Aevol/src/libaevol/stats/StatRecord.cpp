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
#include <list>



// =================================================================
//                            Project Files
// =================================================================
#include "StatRecord.h"

#include "ExpManager.h"
#include "ExpSetup.h"
#include "Individual.h"
#include "GeneticUnit.h"
#include "ReplicationReport.h"
#include "DnaReplicationReport.h"

#ifdef __REGUL
#include "raevol/Rna_R.h"
#include "raevol/Protein_R.h"
#include "raevol/Individual_R.h"
#endif

using std::list;


namespace aevol {



//##############################################################################
//                                                                             #
//                             Class StatRecord                            #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
StatRecord::StatRecord(ExpManager * exp_m)
{
  exp_m_ = exp_m;
  initialize_data();
}

StatRecord::StatRecord(const StatRecord &model)
{
  exp_m_ = model.exp_m_;

  pop_size_  = model.pop_size_;

  metabolic_error_         = model.metabolic_error_;
  metabolic_fitness_       = model.metabolic_fitness_;
  parent_metabolic_error_  = model.parent_metabolic_error_;

  secretion_error_         = model.secretion_error_;
  secretion_fitness_       = model.secretion_fitness_;
  parent_secretion_error_  = model.parent_secretion_error_;

  compound_amount_   = model.compound_amount_;

  fitness_ = model.fitness_;

  amount_of_dna_                = model.amount_of_dna_;
  nb_coding_rnas_               = model.nb_coding_rnas_;
  nb_non_coding_rnas_           = model.nb_non_coding_rnas_;
  av_size_coding_rnas_          = model.av_size_coding_rnas_;
  av_size_non_coding_rnas_      = model.av_size_non_coding_rnas_;
  nb_functional_genes_          = model.nb_functional_genes_;
  nb_non_functional_genes_      = model.nb_non_functional_genes_;
  av_size_functional_gene_      = model.av_size_functional_gene_;
  av_size_non_functional_gene_  = model.av_size_non_functional_gene_;

  nb_mut_    = model.nb_mut_;
  nb_rear_   = model.nb_rear_;
  nb_switch_ = model.nb_switch_;
  nb_indels_ = model.nb_indels_;
  nb_dupl_   = model.nb_dupl_;
  nb_del_    = model.nb_del_;
  nb_trans_  = model.nb_trans_;
  nb_inv_    = model.nb_inv_;

  dupl_rate_        = model.dupl_rate_;
  del_rate_         = model.del_rate_;
  trans_rate_       = model.trans_rate_;
  inv_rate_         = model.inv_rate_;
  mean_align_score_ = model.mean_align_score_;

  nb_bases_in_0_CDS_                = model.nb_bases_in_0_CDS_;
  nb_bases_in_0_functional_CDS_     = model.nb_bases_in_0_functional_CDS_;
  nb_bases_in_0_non_functional_CDS_ = model.nb_bases_in_0_non_functional_CDS_;
  nb_bases_in_0_RNA_                = model.nb_bases_in_0_RNA_;
  nb_bases_in_0_coding_RNA_         = model.nb_bases_in_0_coding_RNA_;
  nb_bases_in_0_non_coding_RNA_     = model.nb_bases_in_0_non_coding_RNA_;

  nb_bases_non_essential_                     = model.nb_bases_non_essential_;
  nb_bases_non_essential_including_nf_genes_  = model.nb_bases_non_essential_including_nf_genes_;

  #ifdef __REGUL
    nb_TF_                         = model.nb_TF_;
    nb_pure_TF_                    = model.nb_pure_TF_;
    nb_influences_                 = model.nb_influences_;
    nb_enhancing_influences_       = model.nb_enhancing_influences_;
    nb_operating_influences_       = model.nb_operating_influences_;
    av_value_influences_           = model.av_value_influences_;
    av_value_enhancing_influences_ = model.av_value_enhancing_influences_;
    av_value_operating_influences_ = model.av_value_operating_influences_;
  #endif
}


    StatRecord::StatRecord(ExpSetup* exp_s,
                           Individual* indiv,
                           ReplicationReport* replic_report,
                           chrom_or_gen_unit chrom_or_gu,
                           bool compute_non_coding) {
        record_type_ = INDIV;

        // ---------------
        // Simulation data
        // ---------------
        pop_size_ = 0; // The pop_size value is irrelevant when dealing with a single
        // individual. It is present for column alignment.

#ifdef __REGUL
        int32_t nb_activators = 0;
  int32_t nb_operators = 0;
  double mean_activator_activity = 0.0;
  double mean_operator_activity = 0.0;

  Individual_R* indiv_r = dynamic_cast<Individual_R*>(indiv);

  for (auto& rna: indiv_r->get_rna_list_coding()) {
    for (int i = 0; i < ((Rna_R*)rna)->nb_influences(); i++) {
      //compute the activity
      if (((Rna_R*)rna)->_enhancing_coef_list[i] > 0)
      {
        nb_activators++;
	      mean_activator_activity += ((Rna_R*)rna)->_enhancing_coef_list[i];
      }

      if (((Rna_R*)rna)->_operating_coef_list[i] > 0)
      {
	      nb_operators++;
	      mean_operator_activity += ((Rna_R*)rna)->_operating_coef_list[i];
      }
    }
  }


  nb_enhancing_influences_       = nb_activators;
  nb_operating_influences_       = nb_operators;
  nb_influences_                 = nb_operating_influences_ + nb_enhancing_influences_;
  av_value_influences_           = ( mean_activator_activity + mean_operator_activity ) / double ( nb_activators + nb_operators);
  av_value_enhancing_influences_ = ( mean_activator_activity ) / double ( nb_activators );
  av_value_operating_influences_ = ( mean_operator_activity ) / double ( nb_operators);

  //ajout raevol_yo_2
  int32_t nb_TF = 0;
  int32_t nb_pure_TF = 0;

  for (auto& prot: indiv_r->protein_list()) {
    if(prot->is_functional())
    {
      if(!((Protein_R*)prot)->is_TF_)
      {
				nb_TF+=1;
      }
    }
    else
    {
      if(((Protein_R*)prot)->is_TF_)
      {
				nb_TF+=1;
				nb_pure_TF+=1;
      }
    }
  }

	nb_TF_ = nb_TF;
	nb_pure_TF_ = nb_pure_TF;
#endif

        // TODO : These conditions are not well managed!!!
        if (indiv->nb_genetic_units() == 1) {
            // -------------------------------------------------
            // Compute statistical data for the given individual
            // -------------------------------------------------
            if (compute_non_coding)
                indiv->compute_non_coding();

            const GeneticUnit& gen_unit = *indiv->genetic_unit_list().begin();

            // Metabolic error stats
            metabolic_error_ = indiv->dist_to_target_by_feature(METABOLISM);
            metabolic_fitness_ = indiv->fitness_by_feature(METABOLISM);
            parent_metabolic_error_ = (replic_report != NULL) ?
                                      replic_report->parent_metabolic_error() :
                                      0.0;

            // Fitness
            fitness_ = indiv->fitness();

            // Secretion stats
            if (exp_s->with_secretion()) {
                secretion_error_   = indiv->dist_to_target_by_feature(SECRETION);
                secretion_fitness_ = indiv->fitness_by_feature(SECRETION);
                compound_amount_   = indiv->grid_cell()->compound_amount();
                parent_secretion_error_ = 0.0;

                if (replic_report != NULL)
                {
                    parent_secretion_error_ = replic_report->parent_secretion_error();
                }
            }
            else
            {
                secretion_error_   = 0.0;
                secretion_fitness_ = 0.0;
                compound_amount_   = 0.0;
                parent_secretion_error_ = 0.0;
            }

            // Genes and RNA stats
            amount_of_dna_               = gen_unit.dna()->length();
            nb_coding_rnas_              = gen_unit.nb_coding_RNAs();
            nb_non_coding_rnas_          = gen_unit.nb_non_coding_RNAs();
            av_size_coding_rnas_         = gen_unit.av_size_coding_RNAs();
            av_size_non_coding_rnas_     = gen_unit.av_size_non_coding_RNAs();
            nb_functional_genes_         = gen_unit.nb_functional_genes();
            nb_non_functional_genes_     = gen_unit.nb_non_functional_genes();
            av_size_functional_gene_     = gen_unit.av_size_functional_genes();
            av_size_non_functional_gene_ = gen_unit.av_size_non_functional_genes();

            // Non coding stats
            if (compute_non_coding) {
                nb_bases_in_0_CDS_                = gen_unit.nb_bases_in_0_CDS();
                nb_bases_in_0_functional_CDS_     = gen_unit.nb_bases_in_0_functional_CDS();
                nb_bases_in_0_non_functional_CDS_ = gen_unit.nb_bases_in_0_non_functional_CDS();
                nb_bases_in_0_RNA_                = gen_unit.nb_bases_in_0_RNA();
                nb_bases_in_0_coding_RNA_         = gen_unit.nb_bases_in_0_coding_RNA();
                nb_bases_in_0_non_coding_RNA_     = gen_unit.nb_bases_in_0_non_coding_RNA();

                nb_bases_non_essential_                     = gen_unit.nb_bases_non_essential();
                nb_bases_non_essential_including_nf_genes_  = gen_unit.nb_bases_non_essential_including_nf_genes();
            }

            // Mutation stats
            if (replic_report != NULL)
            {
                nb_mut_    = replic_report->nb(S_MUT);
                nb_rear_   = replic_report->nb(REARR);
                nb_switch_ = replic_report->nb(SWITCH);
                nb_indels_ = replic_report->nb(INDEL);
                nb_dupl_   = replic_report->nb(DUPL);
                nb_del_    = replic_report->nb(DEL);
                nb_trans_  = replic_report->nb(TRANS);
                nb_inv_    = replic_report->nb(INV);

                // Rearrangement rate stats
                int32_t parent_genome_size = replic_report->parent_genome_size();
                dupl_rate_  = nb_dupl_  / parent_genome_size;
                del_rate_   = nb_del_   / parent_genome_size;
                trans_rate_ = nb_trans_ / parent_genome_size;
                inv_rate_   = nb_inv_   / parent_genome_size;

                //~ // <DEBUG>
                //~ if (nb_dupl_ + nb_del_ + nb_trans_ + nb_inv_ != 0)
                //~ {
                //~ printf("nb_dupl_ : %"PRId32"\n_nb_del : %"PRId32"\n_nb_trans : %"PRId32"\n_nb_inv : %"PRId32"\n",
                //~ (int32_t) nb_dupl_, (int32_t) nb_del_, (int32_t) nb_trans_, (int32_t) nb_inv_);
                //~ printf("parent genome size : %"PRId32"\n", parent_genome_size);
                //~ printf("dupl_rate_ : %f\n_del_rate : %f\n_trans_rate : %f\n_inv_rate : %f\n",
                //~ dupl_rate_, del_rate_, trans_rate_, inv_rate_);
                //~ getchar();
                //~ }
                //~ // </DEBUG>

                mean_align_score_ = replic_report->mean_align_score();
            }
        }
        else if (chrom_or_gu == ALL_GU)
        {
            // -------------------------------------------------
            // Compute statistical data for the given individual
            // -------------------------------------------------
            // Metabolic error stats
            metabolic_error_ = (double) indiv->dist_to_target_by_feature(METABOLISM);
            metabolic_fitness_ = (double) indiv->fitness_by_feature(METABOLISM);
            parent_metabolic_error_ = (replic_report != NULL) ? replic_report->parent_metabolic_error() : 0.0;

            // Fitness
            fitness_ = indiv->fitness();

            // Secretion stats
            if (exp_s->with_secretion()) {
                secretion_error_ = (double) indiv->dist_to_target_by_feature(SECRETION);
                secretion_fitness_ = (double) indiv->fitness_by_feature(SECRETION);
                compound_amount_   = (double) indiv->grid_cell()->compound_amount();
                parent_secretion_error_ = 0.0;

                if (replic_report != NULL)
                {
                    parent_secretion_error_ = replic_report->parent_secretion_error();
                }
            }
            else
            {
                secretion_error_   = 0.0;
                secretion_fitness_ = 0.0;
                compound_amount_   = 0.0;
                parent_secretion_error_ = 0.0;
            }

            for (auto& gen_unit: indiv->genetic_unit_list_nonconst()) {
                // Genes and RNA stats
                amount_of_dna_               += gen_unit.dna()->length();
                nb_coding_rnas_              += gen_unit.nb_coding_RNAs();
                nb_non_coding_rnas_          += gen_unit.nb_non_coding_RNAs();
                av_size_coding_rnas_         += gen_unit.av_size_coding_RNAs();
                av_size_non_coding_rnas_     += gen_unit.av_size_non_coding_RNAs();
                nb_functional_genes_         += gen_unit.nb_functional_genes();
                nb_non_functional_genes_     += gen_unit.nb_non_functional_genes();
                av_size_functional_gene_     += gen_unit.av_size_functional_genes();
                av_size_non_functional_gene_ += gen_unit.av_size_non_functional_genes();

                // Non coding stats
                if (compute_non_coding)
                {
                    nb_bases_in_0_CDS_                += gen_unit.nb_bases_in_0_CDS();
                    nb_bases_in_0_functional_CDS_     += gen_unit.nb_bases_in_0_functional_CDS();
                    nb_bases_in_0_non_functional_CDS_ += gen_unit.nb_bases_in_0_non_functional_CDS();
                    nb_bases_in_0_RNA_                += gen_unit.nb_bases_in_0_RNA();
                    nb_bases_in_0_coding_RNA_         += gen_unit.nb_bases_in_0_coding_RNA();
                    nb_bases_in_0_non_coding_RNA_     += gen_unit.nb_bases_in_0_non_coding_RNA();

                    nb_bases_non_essential_                     += gen_unit.nb_bases_non_essential();
                    nb_bases_non_essential_including_nf_genes_  += gen_unit.nb_bases_non_essential_including_nf_genes();
                }

                // Mutation stats
                if (replic_report != NULL)
                {
                    nb_mut_    += replic_report->nb(S_MUT);
                    nb_rear_   += replic_report->nb(REARR);
                    nb_switch_ += replic_report->nb(SWITCH);
                    nb_indels_ += replic_report->nb(INDEL);
                    nb_dupl_   += replic_report->nb(DUPL);
                    nb_del_    += replic_report->nb(DEL);
                    nb_trans_  += replic_report->nb(TRANS);
                    nb_inv_    += replic_report->nb(INV);
                }
            }

            // Rearrangement rate stats
            if (replic_report != NULL)
            {
                int32_t parent_genome_size = replic_report->parent_genome_size();
                dupl_rate_  = nb_dupl_  / parent_genome_size;
                del_rate_   = nb_del_   / parent_genome_size;
                trans_rate_ = nb_trans_ / parent_genome_size;
                inv_rate_   = nb_inv_   / parent_genome_size;
                mean_align_score_ = replic_report->mean_align_score();
            }
        }
        else // => We have a multi-GU individual and we want only the main chromosome or only the plasmids
            // WARNING (TODO) As it is coded, this will work only if there is ONE SINGLE PLASMID!
        {
            if (chrom_or_gu != PLASMIDS and chrom_or_gu != CHROM) {
                printf("%s: error: StatRecord called with inappropriate `chrom_or_gu`\n", __FILE__);
                exit(EXIT_FAILURE);
            }

            GeneticUnit& gen_unit = (chrom_or_gu == PLASMIDS) ?
                                    *std::next(indiv->genetic_unit_list_nonconst().begin()) :
                                    indiv->genetic_unit_list_nonconst().front();

            // -------------------------------------------------
            // Compute statistical data for the given individual
            // -------------------------------------------------
            // Metabolic error stats
            metabolic_error_ = (double) gen_unit.dist_to_target_by_feature(METABOLISM);
            metabolic_fitness_ = (double) gen_unit.fitness_by_feature(METABOLISM);
            parent_metabolic_error_ = (replic_report != NULL) ? replic_report->parent_metabolic_error() : 0.0;

            // Fitness
            fitness_ = indiv->fitness();

            // Secretion stats
            if (exp_s->with_secretion()) {
                secretion_error_ = (double) gen_unit.dist_to_target_by_feature(SECRETION);
                secretion_fitness_ = (double) gen_unit.fitness_by_feature(SECRETION);
                compound_amount_   = (double) indiv->grid_cell()->compound_amount();
                parent_secretion_error_ = 0.0;

                if (replic_report != NULL)
                {
                    parent_secretion_error_ = replic_report->parent_secretion_error();
                }
            }
            else
            {
                secretion_error_   = 0.0;
                secretion_fitness_ = 0.0;
                compound_amount_   = 0.0;
                parent_secretion_error_ = 0.0;
            }

            // Genes and RNA stats
            amount_of_dna_               = gen_unit.dna()->length();
            nb_coding_rnas_              = gen_unit.nb_coding_RNAs();
            nb_non_coding_rnas_          = gen_unit.nb_non_coding_RNAs();
            av_size_coding_rnas_         = gen_unit.av_size_coding_RNAs();
            av_size_non_coding_rnas_     = gen_unit.av_size_non_coding_RNAs();
            nb_functional_genes_         = gen_unit.nb_functional_genes();
            nb_non_functional_genes_     = gen_unit.nb_non_functional_genes();
            av_size_functional_gene_     = gen_unit.av_size_functional_genes();
            av_size_non_functional_gene_ = gen_unit.av_size_non_functional_genes();

            // Non coding stats
            if (compute_non_coding)
            {
                nb_bases_in_0_CDS_                  = gen_unit.nb_bases_in_0_CDS();
                nb_bases_in_0_functional_CDS_       = gen_unit.nb_bases_in_0_functional_CDS();
                nb_bases_in_0_non_functional_CDS_   = gen_unit.nb_bases_in_0_non_functional_CDS();
                nb_bases_in_0_RNA_                  = gen_unit.nb_bases_in_0_RNA();
                nb_bases_in_0_coding_RNA_           = gen_unit.nb_bases_in_0_coding_RNA();
                nb_bases_in_0_non_coding_RNA_       = gen_unit.nb_bases_in_0_non_coding_RNA();

                nb_bases_non_essential_                     = gen_unit.nb_bases_non_essential();
                nb_bases_non_essential_including_nf_genes_  = gen_unit.nb_bases_non_essential_including_nf_genes();
            }

            // Mutation stats
            // TODO <david.parsons@inria.fr> Disabled
//    if (gen_unit.dna()->replication_report() != NULL)
//    {
//      nb_mut_    = gen_unit.dna()->replication_report()->nb(S_MUT);
//      nb_rear_   = gen_unit.dna()->replication_report()->nb(REARR);
//      nb_switch_ = gen_unit.dna()->replication_report()->nb(SWITCH);
//      nb_indels_ = gen_unit.dna()->replication_report()->nb(INDEL);
//      nb_dupl_   = gen_unit.dna()->replication_report()->nb(DUPL);
//      nb_del_    = gen_unit.dna()->replication_report()->nb(DEL);
//      nb_trans_  = gen_unit.dna()->replication_report()->nb(TRANS);
//      nb_inv_    = gen_unit.dna()->replication_report()->nb(INV);
//    }

            // Rearrangement rate stats
            if (replic_report != NULL)
            {
                int32_t parent_genome_size = replic_report->parent_genome_size();
                dupl_rate_  = nb_dupl_  / parent_genome_size;
                del_rate_   = nb_del_   / parent_genome_size;
                trans_rate_ = nb_trans_ / parent_genome_size;
                inv_rate_   = nb_inv_   / parent_genome_size;
                mean_align_score_ = replic_report->mean_align_score();
            }
        }
    }
/*
StatRecord::StatRecord(ExpManager* exp_m,
                       Individual* indiv,
                       chrom_or_gen_unit chrom_or_gu,
                       bool compute_non_coding)
{
  exp_m_ = exp_m;
  initialize_data();
  record_type_ = INDIV;

  // ---------------
  // Simulation data
  // ---------------
  pop_size_ = 0; // The pop_size value is irrelevant when dealing with a single
                 // individual. It is present for column alignment.

  #ifdef __REGUL
  int32_t nb_activators = 0;
  int32_t nb_operators = 0;
  double mean_activator_activity = 0.0;
  double mean_operator_activity = 0.0;

  Individual_R* indiv_r = dynamic_cast<Individual_R*>(indiv);

  for (auto& rna: indiv_r->get_rna_list_coding()) {
    for (unsigned int i = 0; i < ((Rna_R*)rna)->nb_influences(); i++) {
      //compute the activity
      if (((Rna_R*)rna)->_enhancing_coef_list[i] > 0)
      {
        nb_activators++;
	      mean_activator_activity += ((Rna_R*)rna)->_enhancing_coef_list[i];
      }

      if (((Rna_R*)rna)->_operating_coef_list[i] > 0)
      {
	      nb_operators++;
	      mean_operator_activity += ((Rna_R*)rna)->_operating_coef_list[i];
      }
    }
  }


  nb_enhancing_influences_       = nb_activators;
  nb_operating_influences_       = nb_operators;
  nb_influences_                 = nb_operating_influences_ + nb_enhancing_influences_;
  av_value_influences_           = ( mean_activator_activity + mean_operator_activity ) / double ( nb_activators + nb_operators);
  av_value_enhancing_influences_ = ( mean_activator_activity ) / double ( nb_activators );
  av_value_operating_influences_ = ( mean_operator_activity ) / double ( nb_operators);

  //ajout raevol_yo_2
  int32_t nb_TF = 0;
  int32_t nb_pure_TF = 0;

  for (auto& prot: indiv_r->protein_list()) {
    if(prot->is_functional())
    {
      if(!((Protein_R*)prot)->is_TF_)
      {
				nb_TF+=1;
      }
    }
    else
    {
      if(((Protein_R*)prot)->is_TF_)
      {
				nb_TF+=1;
				nb_pure_TF+=1;
      }
    }
  }

	nb_TF_ = nb_TF;
	nb_pure_TF_ = nb_pure_TF;
  #endif  
    
  // TODO : These conditions are not well managed!!!
  if (indiv->nb_genetic_units() == 1)
  { // One single Genetic Unit
    // -------------------------------------------------
    // Compute statistical data for the given individual
    // -------------------------------------------------
    ReplicationReport* replic_report = nullptr;
    if (exp_m_->tree() != nullptr)
      replic_report = exp_m_->tree()->report_by_index(AeTime::time(),
                                                      indiv->grid_cell()->x() *
                                                      indiv->exp_m()->grid_height()
                                                      + indiv->grid_cell()->y());

    if (compute_non_coding)
      indiv->compute_non_coding();

    const GeneticUnit& gen_unit = *indiv->genetic_unit_list().begin();

    // Metabolic error stats
    metabolic_error_ = indiv->dist_to_target_by_feature(METABOLISM);
    metabolic_fitness_ = indiv->fitness_by_feature(METABOLISM);
    parent_metabolic_error_ = (replic_report != NULL) ?
                              replic_report->parent_metabolic_error() :
                              0.0;

    // Fitness
    fitness_ = indiv->fitness();

    // Secretion stats
    if (exp_m_->with_secretion())
    {
       secretion_error_   = indiv->dist_to_target_by_feature(SECRETION);
       secretion_fitness_ = indiv->fitness_by_feature(SECRETION);
       compound_amount_   = indiv->grid_cell()->compound_amount();
       parent_secretion_error_ = 0.0;

      if (replic_report != NULL)
      {
        parent_secretion_error_ = replic_report->parent_secretion_error();
      }
    }
    else
    {
      secretion_error_   = 0.0;
      secretion_fitness_ = 0.0;
      compound_amount_   = 0.0;
      parent_secretion_error_ = 0.0;
    }

    // Genes and RNA stats
    amount_of_dna_               = gen_unit.dna()->length();
    nb_coding_rnas_              = gen_unit.nb_coding_RNAs();
    nb_non_coding_rnas_          = gen_unit.nb_non_coding_RNAs();
    av_size_coding_rnas_         = gen_unit.av_size_coding_RNAs();
    av_size_non_coding_rnas_     = gen_unit.av_size_non_coding_RNAs();
    nb_functional_genes_         = gen_unit.nb_functional_genes();
    nb_non_functional_genes_     = gen_unit.nb_non_functional_genes();
    av_size_functional_gene_     = gen_unit.av_size_functional_genes();
    av_size_non_functional_gene_ = gen_unit.av_size_non_functional_genes();

    // Non coding stats
    if (compute_non_coding) {
      nb_bases_in_0_CDS_                = gen_unit.nb_bases_in_0_CDS();
      nb_bases_in_0_functional_CDS_     = gen_unit.nb_bases_in_0_functional_CDS();
      nb_bases_in_0_non_functional_CDS_ = gen_unit.nb_bases_in_0_non_functional_CDS();
      nb_bases_in_0_RNA_                = gen_unit.nb_bases_in_0_RNA();
      nb_bases_in_0_coding_RNA_         = gen_unit.nb_bases_in_0_coding_RNA();
      nb_bases_in_0_non_coding_RNA_     = gen_unit.nb_bases_in_0_non_coding_RNA();

      nb_bases_non_essential_                     = gen_unit.nb_bases_non_essential();
      nb_bases_non_essential_including_nf_genes_  = gen_unit.nb_bases_non_essential_including_nf_genes();
    }

    // Mutation stats
    if (replic_report != NULL)
    {
      nb_mut_    = replic_report->nb(S_MUT);
      nb_rear_   = replic_report->nb(REARR);
      nb_switch_ = replic_report->nb(SWITCH);
      nb_indels_ = replic_report->nb(INDEL);
      nb_dupl_   = replic_report->nb(DUPL);
      nb_del_    = replic_report->nb(DEL);
      nb_trans_  = replic_report->nb(TRANS);
      nb_inv_    = replic_report->nb(INV);

      // Rearrangement rate stats
      int32_t parent_genome_size = replic_report->parent_genome_size();
      dupl_rate_  = nb_dupl_  / parent_genome_size;
      del_rate_   = nb_del_   / parent_genome_size;
      trans_rate_ = nb_trans_ / parent_genome_size;
      inv_rate_   = nb_inv_   / parent_genome_size;

      //~ // <DEBUG>
      //~ if (nb_dupl_ + nb_del_ + nb_trans_ + nb_inv_ != 0)
      //~ {
        //~ printf("nb_dupl_ : %"PRId32"\n_nb_del : %"PRId32"\n_nb_trans : %"PRId32"\n_nb_inv : %"PRId32"\n",
                //~ (int32_t) nb_dupl_, (int32_t) nb_del_, (int32_t) nb_trans_, (int32_t) nb_inv_);
        //~ printf("parent genome size : %"PRId32"\n", parent_genome_size);
        //~ printf("dupl_rate_ : %f\n_del_rate : %f\n_trans_rate : %f\n_inv_rate : %f\n",
                //~ dupl_rate_, del_rate_, trans_rate_, inv_rate_);
        //~ getchar();
      //~ }
      //~ // </DEBUG>

      mean_align_score_ = replic_report->mean_align_score();
    }
  }
  else if (chrom_or_gu == ALL_GU)
  {
    // -------------------------------------------------
    // Compute statistical data for the given individual
    // -------------------------------------------------
    ReplicationReport* replic_report = nullptr;
    if (exp_m_->tree() != nullptr)
      printf("Size %d %d %d\n",indiv->grid_cell()->x(),indiv->grid_cell()->y(),indiv->id());

      replic_report = exp_m_->tree()->report_by_index(AeTime::time(),
                                                      indiv->grid_cell()->x() *
                                                      indiv->exp_m()->grid_height()
                                                      + indiv->grid_cell()->y());

    // Metabolic error stats
    metabolic_error_ = (double) indiv->dist_to_target_by_feature(METABOLISM);
    metabolic_fitness_ = (double) indiv->fitness_by_feature(METABOLISM);
    parent_metabolic_error_ = (replic_report != NULL) ? replic_report->parent_metabolic_error() : 0.0;

    // Fitness
    fitness_ = indiv->fitness();

    // Secretion stats
    if (exp_m_->with_secretion())
    {
       secretion_error_ = (double) indiv->dist_to_target_by_feature(SECRETION);
       secretion_fitness_ = (double) indiv->fitness_by_feature(SECRETION);
       compound_amount_   = (double) indiv->grid_cell()->compound_amount();
       parent_secretion_error_ = 0.0;

      if (replic_report != NULL)
      {
        parent_secretion_error_ = replic_report->parent_secretion_error();
      }
    }
    else
    {
      secretion_error_   = 0.0;
      secretion_fitness_ = 0.0;
      compound_amount_   = 0.0;
      parent_secretion_error_ = 0.0;
    }

    for (auto& gen_unit: indiv->genetic_unit_list_nonconst()) {
      // Genes and RNA stats
      amount_of_dna_               += gen_unit.dna()->length();
      nb_coding_rnas_              += gen_unit.nb_coding_RNAs();
      nb_non_coding_rnas_          += gen_unit.nb_non_coding_RNAs();
      av_size_coding_rnas_         += gen_unit.av_size_coding_RNAs();
      av_size_non_coding_rnas_     += gen_unit.av_size_non_coding_RNAs();
      nb_functional_genes_         += gen_unit.nb_functional_genes();
      nb_non_functional_genes_     += gen_unit.nb_non_functional_genes();
      av_size_functional_gene_     += gen_unit.av_size_functional_genes();
      av_size_non_functional_gene_ += gen_unit.av_size_non_functional_genes();

      // Non coding stats
      if (compute_non_coding)
      {
        nb_bases_in_0_CDS_                += gen_unit.nb_bases_in_0_CDS();
        nb_bases_in_0_functional_CDS_     += gen_unit.nb_bases_in_0_functional_CDS();
        nb_bases_in_0_non_functional_CDS_ += gen_unit.nb_bases_in_0_non_functional_CDS();
        nb_bases_in_0_RNA_                += gen_unit.nb_bases_in_0_RNA();
        nb_bases_in_0_coding_RNA_         += gen_unit.nb_bases_in_0_coding_RNA();
        nb_bases_in_0_non_coding_RNA_     += gen_unit.nb_bases_in_0_non_coding_RNA();

        nb_bases_non_essential_                     += gen_unit.nb_bases_non_essential();
        nb_bases_non_essential_including_nf_genes_  += gen_unit.nb_bases_non_essential_including_nf_genes();
      }

      // Mutation stats
      if (replic_report != NULL)
      {
        nb_mut_    += replic_report->nb(S_MUT);
        nb_rear_   += replic_report->nb(REARR);
        nb_switch_ += replic_report->nb(SWITCH);
        nb_indels_ += replic_report->nb(INDEL);
        nb_dupl_   += replic_report->nb(DUPL);
        nb_del_    += replic_report->nb(DEL);
        nb_trans_  += replic_report->nb(TRANS);
        nb_inv_    += replic_report->nb(INV);
      }
    }

    // Rearrangement rate stats
    if (replic_report != NULL)
    {
      int32_t parent_genome_size = replic_report->parent_genome_size();
      dupl_rate_  = nb_dupl_  / parent_genome_size;
      del_rate_   = nb_del_   / parent_genome_size;
      trans_rate_ = nb_trans_ / parent_genome_size;
      inv_rate_   = nb_inv_   / parent_genome_size;
      mean_align_score_ = replic_report->mean_align_score();
    }
  }
  else // => We have a multi-GU individual and we want only the main chromosome or only the plasmids
  // WARNING (TODO) As it is coded, this will work only if there is ONE SINGLE PLASMID!
  {
    if (chrom_or_gu != PLASMIDS and chrom_or_gu != CHROM) {
      printf("%s: error: StatRecord called with inappropriate `chrom_or_gu`\n", __FILE__);
      exit(EXIT_FAILURE);
    }

    GeneticUnit& gen_unit = (chrom_or_gu == PLASMIDS) ?
        *std::next(indiv->genetic_unit_list_nonconst().begin()) :
        indiv->genetic_unit_list_nonconst().front();

    // -------------------------------------------------
    // Compute statistical data for the given individual
    // -------------------------------------------------
    ReplicationReport* replic_report = nullptr;
    if (exp_m_->tree() != nullptr)
      replic_report = exp_m_->tree()->report_by_index(AeTime::time(),
                                                      indiv->grid_cell()->x() *
                                                      indiv->exp_m()->grid_height()
                                                      + indiv->grid_cell()->y());

    // Metabolic error stats
    metabolic_error_ = (double) gen_unit.dist_to_target_by_feature(METABOLISM);
    metabolic_fitness_ = (double) gen_unit.fitness_by_feature(METABOLISM);
    parent_metabolic_error_ = (replic_report != NULL) ? replic_report->parent_metabolic_error() : 0.0;

    // Fitness
    fitness_ = indiv->fitness();

    // Secretion stats
    if (exp_m_->with_secretion())
    {
       secretion_error_ = (double) gen_unit.dist_to_target_by_feature(SECRETION);
       secretion_fitness_ = (double) gen_unit.fitness_by_feature(SECRETION);
       compound_amount_   = (double) indiv->grid_cell()->compound_amount();
       parent_secretion_error_ = 0.0;

      if (replic_report != NULL)
      {
        parent_secretion_error_ = replic_report->parent_secretion_error();
      }
    }
    else
    {
      secretion_error_   = 0.0;
      secretion_fitness_ = 0.0;
      compound_amount_   = 0.0;
      parent_secretion_error_ = 0.0;
    }

      // Genes and RNA stats
    amount_of_dna_               = gen_unit.dna()->length();
    nb_coding_rnas_              = gen_unit.nb_coding_RNAs();
    nb_non_coding_rnas_          = gen_unit.nb_non_coding_RNAs();
    av_size_coding_rnas_         = gen_unit.av_size_coding_RNAs();
    av_size_non_coding_rnas_     = gen_unit.av_size_non_coding_RNAs();
    nb_functional_genes_         = gen_unit.nb_functional_genes();
    nb_non_functional_genes_     = gen_unit.nb_non_functional_genes();
    av_size_functional_gene_     = gen_unit.av_size_functional_genes();
    av_size_non_functional_gene_ = gen_unit.av_size_non_functional_genes();

      // Non coding stats
    if (compute_non_coding)
    {
      nb_bases_in_0_CDS_                  = gen_unit.nb_bases_in_0_CDS();
      nb_bases_in_0_functional_CDS_       = gen_unit.nb_bases_in_0_functional_CDS();
      nb_bases_in_0_non_functional_CDS_   = gen_unit.nb_bases_in_0_non_functional_CDS();
      nb_bases_in_0_RNA_                  = gen_unit.nb_bases_in_0_RNA();
      nb_bases_in_0_coding_RNA_           = gen_unit.nb_bases_in_0_coding_RNA();
      nb_bases_in_0_non_coding_RNA_       = gen_unit.nb_bases_in_0_non_coding_RNA();

      nb_bases_non_essential_                     = gen_unit.nb_bases_non_essential();
      nb_bases_non_essential_including_nf_genes_  = gen_unit.nb_bases_non_essential_including_nf_genes();
    }

    // Mutation stats
    // TODO <david.parsons@inria.fr> Disabled
//    if (gen_unit.dna()->replication_report() != NULL)
//    {
//      nb_mut_    = gen_unit.dna()->replication_report()->nb(S_MUT);
//      nb_rear_   = gen_unit.dna()->replication_report()->nb(REARR);
//      nb_switch_ = gen_unit.dna()->replication_report()->nb(SWITCH);
//      nb_indels_ = gen_unit.dna()->replication_report()->nb(INDEL);
//      nb_dupl_   = gen_unit.dna()->replication_report()->nb(DUPL);
//      nb_del_    = gen_unit.dna()->replication_report()->nb(DEL);
//      nb_trans_  = gen_unit.dna()->replication_report()->nb(TRANS);
//      nb_inv_    = gen_unit.dna()->replication_report()->nb(INV);
//    }

    // Rearrangement rate stats
    if (replic_report != NULL)
    {
      int32_t parent_genome_size = replic_report->parent_genome_size();
      dupl_rate_  = nb_dupl_  / parent_genome_size;
      del_rate_   = nb_del_   / parent_genome_size;
      trans_rate_ = nb_trans_ / parent_genome_size;
      inv_rate_   = nb_inv_   / parent_genome_size;
      mean_align_score_ = replic_report->mean_align_score();
    }
  }
}*/

// Calculate average statistics for all the recorded values
    StatRecord::StatRecord(ExpSetup* exp_s,
                           std::list<std::pair<Individual*,
                                   ReplicationReport*>> annotated_indivs,
                           chrom_or_gen_unit chrom_or_gu) {
        record_type_ = POP;

        // ---------------
        // Simulation data
        // ---------------
        pop_size_ = static_cast<int32_t>(annotated_indivs.size());

        // ------------------------------------------------------------------
        // Compute statistical data for the each individual in the population
        // ------------------------------------------------------------------
        for (const auto& annotated_indiv : annotated_indivs) {
            StatRecord indiv_stat_record(exp_s,
                                         annotated_indiv.first,
                                         annotated_indiv.second,
                                         chrom_or_gu, false);
            this->add(&indiv_stat_record, annotated_indiv.first->id());
        }


        // ------------------------------------------------------------------
        // Divide every accumulator by the number of indivs in the population
        // ------------------------------------------------------------------
        this->divide(pop_size_);
    }

// Calculate standard deviation for all the recorded values
    StatRecord::StatRecord(ExpSetup* exp_s,
                           std::list<std::pair<Individual*,
                                   ReplicationReport*>> annotated_indivs,
                           const StatRecord * means,
                           chrom_or_gen_unit chrom_or_gu) {
        record_type_ = STDEVS;

        // ---------------
        // Simulation data
        // ---------------
        pop_size_ = static_cast<int32_t>(annotated_indivs.size());

        // ------------------------------------------------------------------
        // Compute statistical data for the each individual in the population
        // ------------------------------------------------------------------
        for (const auto& annotated_indiv : annotated_indivs) {
            StatRecord indiv_stat_record(exp_s,
                                         annotated_indiv.first,
                                         annotated_indiv.second,
                                         chrom_or_gu, false);
            this->substract_power(means, &indiv_stat_record, 2);
        }

        // ---------------------------------------------------------------------------------
        // Divide every accumulator by the square root of number of indivs in the population
        // ---------------------------------------------------------------------------------
        this->divide(pow((pop_size_-1), 0.5));
    }

    // Calculate skewness for all the recorded values
    StatRecord::StatRecord(ExpSetup* exp_s,
                           std::list<std::pair<Individual*,
                                   ReplicationReport*>> annotated_indivs,
                           const StatRecord* means,
                           const StatRecord* stdevs,
                           chrom_or_gen_unit chrom_or_gu) {
        record_type_ = SKEWNESS;

        // ---------------
        // Simulation data
        // ---------------
        pop_size_ = static_cast<int32_t>(annotated_indivs.size());

        // ------------------------------------------------------------------
        // Compute statistical data for the each individual in the population
        // ------------------------------------------------------------------
        for (const auto& annotated_indiv : annotated_indivs) {
            StatRecord indiv_stat_record(exp_s,
                                         annotated_indiv.first,
                                         annotated_indiv.second,
                                         chrom_or_gu, false);
            this->substract_power(means, &indiv_stat_record, 3);
        }

        this->divide(-pop_size_);

        this->divide_record(stdevs, 3/2);
    }


// =================================================================
//                             Destructors
// =================================================================
StatRecord::~StatRecord()
{
}

// =================================================================
//                            Public Methods
// =================================================================
void StatRecord::initialize_data()
{
  pop_size_  = 0.0;

  metabolic_error_         = 0.0;
  metabolic_fitness_       = 0.0;
  parent_metabolic_error_  = 0.0;

  secretion_error_         = 0.0;
  parent_secretion_error_  = 0.0;

  secretion_fitness_ = 0.0;
  compound_amount_   = 0.0;

  fitness_ = 0.0;

  amount_of_dna_               = 0.0;
  nb_coding_rnas_              = 0.0;
  nb_non_coding_rnas_          = 0.0;
  av_size_coding_rnas_         = 0.0;
  av_size_non_coding_rnas_     = 0.0;
  nb_functional_genes_         = 0.0;
  nb_non_functional_genes_     = 0.0;
  av_size_functional_gene_     = 0.0;
  av_size_non_functional_gene_ = 0.0;

  nb_mut_    = 0.0;
  nb_rear_   = 0.0;
  nb_switch_ = 0.0;
  nb_indels_ = 0.0;
  nb_dupl_   = 0.0;
  nb_del_    = 0.0;
  nb_trans_  = 0.0;
  nb_inv_    = 0.0;

  dupl_rate_  = 0.0;
  del_rate_   = 0.0;
  trans_rate_ = 0.0;
  inv_rate_   = 0.0;
  mean_align_score_ = 0.0;

  nb_bases_in_0_CDS_                = 0.0;
  nb_bases_in_0_functional_CDS_     = 0.0;
  nb_bases_in_0_non_functional_CDS_ = 0.0;
  nb_bases_in_0_RNA_                = 0.0;
  nb_bases_in_0_coding_RNA_         = 0.0;
  nb_bases_in_0_non_coding_RNA_     = 0.0;

  nb_bases_non_essential_                     = 0.0;
  nb_bases_non_essential_including_nf_genes_  = 0.0;

  #ifdef __REGUL
    nb_TF_                         = 0.0;
    nb_pure_TF_                    = 0.0;
    nb_influences_                 = 0.0;
    nb_enhancing_influences_       = 0.0;
    nb_operating_influences_       = 0.0;
    av_value_influences_           = 0.0;
    av_value_enhancing_influences_ = 0.0;
    av_value_operating_influences_ = 0.0;
  #endif
}
    void StatRecord::write_to_file(int64_t time, FILE* stat_file, stats_type stat_type_to_print) const
    {
        if (stat_type_to_print == FITNESS_STATS)
        {
            fprintf(stat_file,
                    "%" PRId64 " %" PRId32 " %e %" PRId32 " %e %e %e %e %e %e %e",
                    time,
                    pop_size_,
                    fitness_,
                    amount_of_dna_,
                    metabolic_error_,
                    parent_metabolic_error_,
                    metabolic_fitness_,
                    secretion_error_,
                    parent_secretion_error_,
                    secretion_fitness_,
                    compound_amount_);
#ifdef __REGUL
            fprintf(stat_file,
          " %" PRId32 " %" PRId32 " %" PRId32 " %f %f %f",
          nb_influences_,
          nb_enhancing_influences_,
          nb_operating_influences_,
          av_value_influences_,
          av_value_enhancing_influences_,
          av_value_operating_influences_);
#endif
            fflush (stat_file);

        }
        if (stat_type_to_print == MUTATION_STATS)
        {
            fprintf(stat_file,
                    "%" PRId64 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 "",
                    time,
                    nb_mut_,
                    nb_rear_,
                    nb_switch_,
                    nb_indels_,
                    nb_dupl_,
                    nb_del_,
                    nb_trans_,
                    nb_inv_);
        }
        if (stat_type_to_print == GENES_STATS)
        {
            fprintf(stat_file,
                    "%" PRId64 " %" PRId32 " %" PRId32 " %f %f %" PRId32 " %" PRId32 " %f %f ",
                    time,
                    nb_coding_rnas_,
                    nb_non_coding_rnas_,
                    av_size_coding_rnas_,
                    av_size_non_coding_rnas_,
                    nb_functional_genes_,
                    nb_non_functional_genes_,
                    av_size_functional_gene_,
                    av_size_non_functional_gene_);
          #ifdef __REGUL
            fprintf(stat_file,
          " %" PRId32 " %" PRId32 ,
          nb_TF_,
          nb_pure_TF_);
          #endif
        }
        if (stat_type_to_print == BP_STATS) {
            if (record_type_ == INDIV) {
                fprintf(stat_file,
                        "%" PRId64 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 "",
                        time,
                        nb_bases_in_0_CDS_,
                        nb_bases_in_0_functional_CDS_,
                        nb_bases_in_0_non_functional_CDS_,
                        nb_bases_in_0_RNA_,
                        nb_bases_in_0_coding_RNA_,
                        nb_bases_in_0_non_coding_RNA_,
                        nb_bases_non_essential_,
                        nb_bases_non_essential_including_nf_genes_);
            }
            else // if record_type_ == POP
            {
                // TO DO (if needed) : base-pair stats for all individuals, not just for the best one.
                //
                // fprintf(stat_file, "%" PRId64 " %f %f %f %f %f %f %f %f",
                //         AeTime::time(),
                //         nb_bases_in_0_CDS_,
                //         nb_bases_in_0_functional_CDS_,
                //         nb_bases_in_0_non_functional_CDS_,
                //         nb_bases_in_0_RNA_,
                //         nb_bases_in_0_coding_RNA_,
                //         nb_bases_in_0_non_coding_RNA_,
                //         nb_bases_non_essential_,
                //         nb_bases_non_essential_including_nf_genes_);
            }
        }
        if (stat_type_to_print == REAR_STATS)
        {
            fprintf(stat_file,
                    "%" PRId64 " %e %e %e %e %f",
                    time,
                    dupl_rate_,
                    del_rate_,
                    trans_rate_,
                    inv_rate_,
                    mean_align_score_);
        }

        fprintf(stat_file, "\n");
    }

void StatRecord::divide(double divisor)
{
  // NB : pop_size is a "global" value and must not be divided.

  fitness_                 /= divisor;

  metabolic_error_         /= divisor;
  parent_metabolic_error_  /= divisor;
  metabolic_fitness_       /= divisor;

  secretion_error_         /= divisor;
  parent_secretion_error_  /= divisor;

  secretion_fitness_       /= divisor;
  compound_amount_         /= divisor;

  amount_of_dna_               /= divisor;
  nb_coding_rnas_              /= divisor;
  nb_non_coding_rnas_          /= divisor;
  av_size_coding_rnas_         /= divisor;
  av_size_non_coding_rnas_     /= divisor;
  nb_functional_genes_          /= divisor;
  nb_non_functional_genes_      /= divisor;
  av_size_functional_gene_      /= divisor;
  av_size_non_functional_gene_  /= divisor;

  nb_mut_    /= divisor;
  nb_rear_   /= divisor;
  nb_switch_ /= divisor;
  nb_indels_ /= divisor;
  nb_dupl_   /= divisor;
  nb_del_    /= divisor;
  nb_trans_  /= divisor;
  nb_inv_    /= divisor;

  //~ printf("PREFINAL %f %f %f %f\n", dupl_rate_, del_rate_, trans_rate_, inv_rate_);
  dupl_rate_  /= divisor;
  del_rate_   /= divisor;
  trans_rate_ /= divisor;
  inv_rate_   /= divisor;
  //~ printf("FINAL %f %f %f %f\n", dupl_rate_, del_rate_, trans_rate_, inv_rate_);
  //~ getchar();
  mean_align_score_ /= divisor;

  nb_bases_in_0_CDS_                /= divisor;
  nb_bases_in_0_functional_CDS_     /= divisor;
  nb_bases_in_0_non_functional_CDS_ /= divisor;
  nb_bases_in_0_RNA_                /= divisor;
  nb_bases_in_0_coding_RNA_         /= divisor;
  nb_bases_in_0_non_coding_RNA_     /= divisor;

  nb_bases_non_essential_                     /= divisor;
  nb_bases_non_essential_including_nf_genes_  /= divisor;

  #ifdef __REGUL
    nb_influences_                 /= divisor;
    nb_enhancing_influences_       /= divisor;
    nb_operating_influences_       /= divisor;
    av_value_influences_           /= divisor;
    av_value_enhancing_influences_ /= divisor;
    av_value_operating_influences_ /= divisor;

    nb_TF_ /= divisor;
    nb_pure_TF_ /= divisor;
  #endif
}


void StatRecord::divide_record(StatRecord const * to_divide, double power)
{
  // NB : pop_size is a "global" value and must not be divided.

  if (to_divide->fitness_ != 0) { fitness_    /= pow(to_divide->fitness_, power); }

  if (to_divide->metabolic_error_ != 0)        { metabolic_error_         /= pow(to_divide->metabolic_error_, power); }
  if (to_divide->parent_metabolic_error_ != 0) { parent_metabolic_error_  /= pow(to_divide->parent_metabolic_error_, power); }
  if (to_divide->metabolic_fitness_ != 0)        { metabolic_fitness_         /= pow(to_divide->metabolic_fitness_, power); }

  if (to_divide->secretion_error_ != 0)        { secretion_error_         /= pow(to_divide->secretion_error_, power); }
  if (to_divide->parent_secretion_error_ != 0) { parent_secretion_error_  /= pow(to_divide->parent_secretion_error_, power); }

  if (to_divide->secretion_fitness_ != 0)       { secretion_fitness_       /= pow(to_divide->secretion_fitness_, power); }
  if (to_divide->compound_amount_ != 0)        { compound_amount_         /= pow(to_divide->compound_amount_, power); }

  if (to_divide->amount_of_dna_ != 0)               { amount_of_dna_               /= pow(to_divide->amount_of_dna_, power); }
  if (to_divide->nb_coding_rnas_ != 0)              { nb_coding_rnas_              /= pow(to_divide->nb_coding_rnas_, power); }
  if (to_divide->nb_non_coding_rnas_ != 0)          { nb_non_coding_rnas_          /= pow(to_divide->nb_non_coding_rnas_, power); }
  if (to_divide->av_size_coding_rnas_ != 0)         { av_size_coding_rnas_         /= pow(to_divide->av_size_coding_rnas_, power); }
  if (to_divide->av_size_non_coding_rnas_ != 0)     { av_size_non_coding_rnas_     /= pow(to_divide->av_size_non_coding_rnas_, power); }
  if (to_divide->nb_functional_genes_ != 0)         { nb_functional_genes_         /= pow(to_divide->nb_functional_genes_, power); }
  if (to_divide->nb_non_functional_genes_ != 0)     { nb_non_functional_genes_     /= pow(to_divide->nb_non_functional_genes_, power); }
  if (to_divide->av_size_functional_gene_ != 0)     { av_size_functional_gene_     /= pow(to_divide->av_size_functional_gene_, power); }
  if (to_divide->av_size_non_functional_gene_ != 0) { av_size_non_functional_gene_ /= pow(to_divide->av_size_non_functional_gene_, power); }

  if (to_divide->nb_mut_ != 0)     { nb_mut_    /= pow(to_divide->nb_mut_, power); }
  if (to_divide->nb_rear_ != 0)    { nb_rear_   /= pow(to_divide->nb_rear_, power); }
  if (to_divide->nb_switch_ != 0)  { nb_switch_ /= pow(to_divide->nb_switch_, power); }
  if (to_divide->nb_indels_ != 0)  { nb_indels_ /= pow(to_divide->nb_indels_, power); }
  if (to_divide->nb_dupl_ != 0)    { nb_dupl_   /= pow(to_divide->nb_dupl_, power); }
  if (to_divide->nb_del_ != 0)     { nb_del_    /= pow(to_divide->nb_del_, power); }
  if (to_divide->nb_trans_ != 0)   { nb_trans_  /= pow(to_divide->nb_trans_, power); }
  if (to_divide->nb_inv_ != 0)     { nb_inv_    /= pow(to_divide->nb_inv_, power); }

  if (to_divide->dupl_rate_ != 0)        { dupl_rate_  /= pow(to_divide->dupl_rate_, power); }
  if (to_divide->del_rate_ != 0)         { del_rate_   /= pow(to_divide->del_rate_, power); }
  if (to_divide->trans_rate_ != 0)       { trans_rate_ /= pow(to_divide->trans_rate_, power); }
  if (to_divide->inv_rate_ != 0)         { inv_rate_   /= pow(to_divide->inv_rate_, power); }
  if (to_divide->mean_align_score_ != 0) { mean_align_score_ /= pow(to_divide->mean_align_score_, power); }

  if (to_divide->nb_bases_in_0_CDS_ != 0)                { nb_bases_in_0_CDS_                /= pow(to_divide->nb_bases_in_0_CDS_, power); }
  if (to_divide->nb_bases_in_0_functional_CDS_ != 0)     { nb_bases_in_0_functional_CDS_     /= pow(to_divide->nb_bases_in_0_functional_CDS_, power); }
  if (to_divide->nb_bases_in_0_non_functional_CDS_ != 0) { nb_bases_in_0_non_functional_CDS_ /= pow(to_divide->nb_bases_in_0_non_functional_CDS_, power); }
  if (to_divide->nb_bases_in_0_RNA_ != 0)                { nb_bases_in_0_RNA_                /= pow(to_divide->nb_bases_in_0_RNA_, power); }
  if (to_divide->nb_bases_in_0_coding_RNA_ != 0)         { nb_bases_in_0_coding_RNA_         /= pow(to_divide->nb_bases_in_0_coding_RNA_, power); }
  if (to_divide->nb_bases_in_0_non_coding_RNA_ != 0)     { nb_bases_in_0_non_coding_RNA_     /= pow(to_divide->nb_bases_in_0_non_coding_RNA_, power); }

  if (to_divide->nb_bases_non_essential_ != 0)                    { nb_bases_non_essential_                     /= pow(to_divide->nb_bases_non_essential_, power); }
  if (to_divide->nb_bases_non_essential_including_nf_genes_ != 0) { nb_bases_non_essential_including_nf_genes_  /= pow(to_divide->nb_bases_non_essential_including_nf_genes_, power); }

  #ifdef __REGUL
    if (to_divide->nb_influences_ != 0)                 { nb_influences_                 /= pow(to_divide->nb_influences_, power); }
    if (to_divide->nb_enhancing_influences_ != 0)       { nb_enhancing_influences_       /= pow(to_divide->nb_enhancing_influences_, power); }
    if (to_divide->nb_TF_ != 0)                 { nb_TF_                 /= pow(to_divide->nb_TF_, power); }
    if (to_divide->nb_pure_TF_ != 0)       { nb_pure_TF_       /= pow(to_divide->nb_pure_TF_, power); }
    if (to_divide->nb_operating_influences_ != 0)       { nb_operating_influences_       /= pow(to_divide->nb_operating_influences_, power); }
    if (to_divide->av_value_influences_ != 0)           { av_value_influences_           /= pow(to_divide->av_value_influences_, power); }
    if (to_divide->av_value_enhancing_influences_ != 0) { av_value_enhancing_influences_ /= pow(to_divide->av_value_enhancing_influences_, power); }
    if (to_divide->av_value_operating_influences_ != 0) { av_value_operating_influences_ /= pow(to_divide->av_value_operating_influences_, power); }
  #endif
}

void StatRecord::add(StatRecord * to_add, int32_t index)
{
  // NB : pop_size is a global values and must not be summed.

  fitness_                 += to_add->fitness_;

  metabolic_error_         += to_add->metabolic_error_;
  parent_metabolic_error_  += to_add->parent_metabolic_error_;
  metabolic_fitness_       += to_add->metabolic_fitness_;

  secretion_error_         += to_add->secretion_error_;
  parent_secretion_error_  += to_add->parent_secretion_error_;

  secretion_fitness_       += to_add->secretion_fitness_;
  compound_amount_         += to_add->compound_amount_;

  amount_of_dna_               += to_add->amount_of_dna_;
  nb_coding_rnas_              += to_add->nb_coding_rnas_;
  nb_non_coding_rnas_          += to_add->nb_non_coding_rnas_;
  av_size_coding_rnas_         += to_add->av_size_coding_rnas_;
  av_size_non_coding_rnas_     += to_add->av_size_non_coding_rnas_;
  nb_functional_genes_         += to_add->nb_functional_genes_;
  nb_non_functional_genes_     += to_add->nb_non_functional_genes_;
  av_size_functional_gene_     += to_add->av_size_functional_gene_;
  av_size_non_functional_gene_ += to_add->av_size_non_functional_gene_;

  nb_mut_    += to_add->nb_mut_;
  nb_rear_   += to_add->nb_rear_;
  nb_switch_ += to_add->nb_switch_;
  nb_indels_ += to_add->nb_indels_;
  nb_dupl_   += to_add->nb_dupl_;
  nb_del_    += to_add->nb_del_;
  nb_trans_  += to_add->nb_trans_;
  nb_inv_    += to_add->nb_inv_;

  dupl_rate_  += to_add->dupl_rate_;
  del_rate_   += to_add->del_rate_;
  trans_rate_ += to_add->trans_rate_;
  inv_rate_   += to_add->inv_rate_;
  //~ printf("%f %f %f %f\n", to_add->dupl_rate_, to_add->del_rate_, to_add->trans_rate_, to_add->inv_rate_);
  mean_align_score_ += to_add->mean_align_score_;

  nb_bases_in_0_CDS_                += to_add->nb_bases_in_0_CDS_;
  nb_bases_in_0_functional_CDS_     += to_add->nb_bases_in_0_functional_CDS_;
  nb_bases_in_0_non_functional_CDS_ += to_add->nb_bases_in_0_non_functional_CDS_;
  nb_bases_in_0_RNA_                += to_add->nb_bases_in_0_RNA_;
  nb_bases_in_0_coding_RNA_         += to_add->nb_bases_in_0_coding_RNA_;
  nb_bases_in_0_non_coding_RNA_     += to_add->nb_bases_in_0_non_coding_RNA_;

  nb_bases_non_essential_                     += to_add->nb_bases_non_essential_;
  nb_bases_non_essential_including_nf_genes_  += to_add->nb_bases_non_essential_including_nf_genes_;

  #ifdef __REGUL
    nb_TF_                         += to_add->nb_TF_;
    nb_pure_TF_                    += to_add->nb_pure_TF_;
    nb_influences_                 += to_add->nb_influences_;
    nb_enhancing_influences_       += to_add->nb_enhancing_influences_;
    nb_operating_influences_       += to_add->nb_operating_influences_;
    av_value_influences_           += to_add->av_value_influences_;
    av_value_enhancing_influences_ += to_add->av_value_enhancing_influences_;
    av_value_operating_influences_ += to_add->av_value_operating_influences_;
  #endif
}

void StatRecord::substract_power(const StatRecord * means,
                                     const StatRecord * to_substract,
                                     double power)
{
  // NB : pop_size is a "global" value and must not be summed.
  fitness_                 += pow(means->fitness_ - to_substract->fitness_, power);

  metabolic_error_         += pow(means->metabolic_error_ - to_substract->metabolic_error_, power);
  parent_metabolic_error_  += pow(means->parent_metabolic_error_ - to_substract->parent_metabolic_error_, power);
  metabolic_fitness_         += pow(means->metabolic_fitness_ - to_substract->metabolic_fitness_, power);

  secretion_error_         += pow(means->secretion_error_ - to_substract->secretion_error_, power);
  parent_secretion_error_  += pow(means->parent_secretion_error_ - to_substract->parent_secretion_error_, power);

  secretion_fitness_       += pow(means->secretion_fitness_ - to_substract->secretion_fitness_, power);
  compound_amount_         += pow(means->compound_amount_ - to_substract->compound_amount_, power);

  amount_of_dna_               += pow(means->amount_of_dna_ - to_substract->amount_of_dna_, power);
  nb_coding_rnas_              += pow(means->nb_coding_rnas_ - to_substract->nb_coding_rnas_, power);
  nb_non_coding_rnas_          += pow(means->nb_non_coding_rnas_ - to_substract->nb_non_coding_rnas_, power);
  av_size_coding_rnas_         += pow(means->av_size_coding_rnas_ - to_substract->av_size_coding_rnas_, power);
  av_size_non_coding_rnas_     += pow(means->av_size_non_coding_rnas_ - to_substract->av_size_non_coding_rnas_, power);
  nb_functional_genes_         += pow(means->nb_functional_genes_ - to_substract->nb_functional_genes_, power);
  nb_non_functional_genes_     += pow(means->nb_non_functional_genes_ - to_substract->nb_non_functional_genes_, power);
  av_size_functional_gene_     += pow(means->av_size_functional_gene_ - to_substract->av_size_functional_gene_, power);
  av_size_non_functional_gene_ += pow(means->av_size_non_functional_gene_ - to_substract->av_size_non_functional_gene_, power);

  nb_mut_    += pow(means->nb_mut_ - to_substract->nb_mut_, power);
  nb_rear_   += pow(means->nb_rear_ - to_substract->nb_rear_, power);
  nb_switch_ += pow(means->nb_switch_ - to_substract->nb_switch_, power);
  nb_indels_ += pow(means->nb_indels_ - to_substract->nb_indels_, power);
  nb_dupl_   += pow(means->nb_dupl_ - to_substract->nb_dupl_, power);
  nb_del_    += pow(means->nb_del_ - to_substract->nb_del_, power);
  nb_trans_  += pow(means->nb_trans_ - to_substract->nb_trans_, power);
  nb_inv_    += pow(means->nb_inv_ - to_substract->nb_inv_, power);

  dupl_rate_  += pow(means->dupl_rate_ - to_substract->dupl_rate_, power);
  del_rate_   += pow(means->del_rate_ - to_substract->del_rate_, power);
  trans_rate_ += pow(means->trans_rate_ - to_substract->trans_rate_, power);
  inv_rate_   += pow(means->inv_rate_ - to_substract->inv_rate_, power);

  mean_align_score_ += pow(means->mean_align_score_ - to_substract->mean_align_score_, power);

  nb_bases_in_0_CDS_                += pow(means->nb_bases_in_0_CDS_ - to_substract->nb_bases_in_0_CDS_, power);
  nb_bases_in_0_functional_CDS_     += pow(means->nb_bases_in_0_functional_CDS_ - to_substract->nb_bases_in_0_functional_CDS_, power);
  nb_bases_in_0_non_functional_CDS_ += pow(means->nb_bases_in_0_non_functional_CDS_ - to_substract->nb_bases_in_0_non_functional_CDS_, power);
  nb_bases_in_0_RNA_                += pow(means->nb_bases_in_0_RNA_ - to_substract->nb_bases_in_0_RNA_, power);
  nb_bases_in_0_coding_RNA_         += pow(means->nb_bases_in_0_coding_RNA_ - to_substract->nb_bases_in_0_coding_RNA_, power);
  nb_bases_in_0_non_coding_RNA_     += pow(means->nb_bases_in_0_non_coding_RNA_ - to_substract->nb_bases_in_0_non_coding_RNA_, power);

  nb_bases_non_essential_                     += pow(means->nb_bases_non_essential_ - to_substract->nb_bases_non_essential_, power);
  nb_bases_non_essential_including_nf_genes_  += pow(means->nb_bases_non_essential_including_nf_genes_ - to_substract->nb_bases_non_essential_including_nf_genes_, power);

  #ifdef __REGUL
  nb_TF_                 += pow(means->nb_TF_ - to_substract->nb_TF_, power);
  nb_pure_TF_                 += pow(means->nb_pure_TF_ - to_substract->nb_pure_TF_, power);
    nb_influences_                 += pow(means->nb_influences_ - to_substract->nb_influences_, power);
    nb_enhancing_influences_       += pow(means->nb_enhancing_influences_ - to_substract->nb_enhancing_influences_, power);
    nb_operating_influences_       += pow(means->nb_operating_influences_ - to_substract->nb_operating_influences_, power);
    av_value_influences_           += pow(means->av_value_influences_ - to_substract->av_value_influences_, power);
    av_value_enhancing_influences_ += pow(means->av_value_enhancing_influences_ - to_substract->av_value_enhancing_influences_, power);
    av_value_operating_influences_ += pow(means->av_value_operating_influences_ - to_substract->av_value_operating_influences_, power);
  #endif
}


// =================================================================
//                           Protected Methods
// =================================================================
} // namespace aevol
