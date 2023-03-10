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

#include "Stats_7.h"
#include "7/Individual_7.h"
#include "7/Protein_7.h"
#include "7/Abstract_Metadata.h"
#include <fstream>
#include <iostream>
#include <string>

namespace aevol {

Stats_7::Stats_7(ExpManager_7* simd_individual, int64_t generation, bool best_or_not, char* prefix
               #ifdef HAVE_MPI
               , int32_t rank
               #endif
               , bool non_coding
) {
  simd_individual_ = simd_individual;
  is_indiv_ = best_or_not;
  generation_ = generation;

  pop_size_ = 0;

  fitness_ = 0;
  metabolic_error_ = 0;

  amount_of_dna_ = 0;
  nb_coding_rnas_ = 0;
  nb_non_coding_rnas_ = 0;

  nb_functional_genes_ = 0;
  nb_non_functional_genes_ = 0;

  nb_mut_ = 0;
  nb_rear_ = 0;
  nb_switch_ = 0;
  nb_indels_ = 0;
  nb_dupl_ = 0;
  nb_del_ = 0;
  nb_trans_ = 0;
  nb_inv_ = 0;


  nb_bases_in_0_CDS_ = 0;
  nb_bases_in_0_functional_CDS_ = 0;
  nb_bases_in_0_non_functional_CDS_ = 0;
  nb_bases_in_0_RNA_ = 0;
  nb_bases_in_0_coding_RNA_ = 0;
  nb_bases_in_0_non_coding_RNA_ = 0;

  nb_bases_non_essential_ = 0;
  nb_bases_non_essential_including_nf_genes_ = 0;

  if (generation_==0) {
      #ifdef PROGENY_STATS
      std::string file = "stats/progeny.csv";
      std::string file_rep = "stats/reproducer_progeny.csv";

      std::ofstream progeny;
      progeny.open(file,std::ofstream::trunc);
      std::ofstream progeny_rep;
      progeny_rep.open(file_rep,std::ofstream::trunc);
      progeny_rep.close();
      progeny.close();
      #endif

      if (prefix != nullptr) {
        std::string iprefix = std::string(prefix)+"/stats_ancestor_best.csv";
        statfile_best_.open(iprefix,std::ofstream::trunc);
      } else if (is_indiv_) {
        std::string file = "stats/stats_simd_best";
        #ifdef HAVE_MPI
        file+="_"+std::to_string(rank);
        #endif
        file+=".csv";
        statfile_best_.open(file,std::ofstream::trunc);
      } else {
        std::string file = "stats/stats_simd_mean";
        #ifdef HAVE_MPI
        file+="_"+std::to_string(rank);
        #endif
        file+=".csv";
        statfile_mean_.open(file,std::ofstream::trunc);
      }
    if (is_indiv_) {
        statfile_best_ << "Generation" << "," <<"pop_size"<<","<<"nb_clones"<<"," << "fitness" << "," << "metabolic_error" << "," <<
                       "amount_of_dna" << "," << "nb_coding_rnas" << "," << "nb_non_coding_rnas" << "," <<
                       "nb_functional_genes" << "," << "nb_non_functional_genes" << "," << "nb_mut"
                       << "," << "nb_switch" << "," << "nb_indels" << "," << "nb_rear" << "," << "nb_dupl" << "," <<
                       "nb_del" << "," << "nb_trans" << "," << "nb_inv" << "," << "dupl_rate" << "," << "del_rate"
                       << "," << "trans_rate" << "," << "inv_rate"
                       #ifdef __REGUL
                       << "," << "nb_link" << "," <<"nb_pos_link" << "," << "nb_neg_link"
                       << "," << "avg_influence" << "," << "avg_pos_influence" << "," << "avg_neg_influence"
                       << "," << "nb_TF" << "," << "nb_pure_TF"
                       #endif
                       ;
                       
        if (non_coding) {
          statfile_best_<< ","<< "nb_bp_not_incl_cds" << "," <<
                           "nb_bp_not_incl_func_cds" << "," <<
                           "nb_bp_not_incl_non_func_cds" << "," <<
                           "nb_bp_not_incl_rna" << "," <<
                           "nb_bp_not_incl_coding_rna" << "," <<
                           "nb_bp_not_incl_non_coding_rna" << "," <<
                           "nb_bp_not_essential" << "," <<
                           "nb_bp_not_essential_incl_non_func_genes";
        }
        statfile_best_<< std::endl;
        statfile_best_.flush();
    } else {
        statfile_mean_ << "Generation" << "," <<"pop_size"<<","<<"nb_clones"<< "," << "fitness" << "," << "metabolic_error" << "," <<
                       "amount_of_dna" << "," << "nb_coding_rnas" << "," << "nb_non_coding_rnas" << "," <<
                       "nb_functional_genes" << "," << "nb_non_functional_genes" << "," << "nb_mut"
                       << "," << "nb_switch" << "," << "nb_indels" << "," << "nb_rear" << "," << "nb_dupl" << "," <<
                       "nb_del" << "," << "nb_trans" << "," << "nb_inv" << "," << "dupl_rate" << "," << "del_rate"
                       << "," << "trans_rate" << "," << "inv_rate"
                       #ifdef __REGUL
                       << "," << "nb_link" << "," <<"nb_pos_link" << "," << "nb_neg_link"
                       << "," << "avg_influence" << "," << "avg_pos_influence" << "," << "avg_neg_influence"
                       << "," << "nb_TF" << "," << "nb_pure_TF"
                       #endif
                       << std::endl;
        statfile_mean_.flush();
    }
  } else {
      printf("Resume without rheader\n");
    std::ifstream tmp_mean;
    std::ifstream tmp_best;

    if (prefix != nullptr) {
        std::string iprefix = std::string(prefix)+"/stats_ancestor_best.csv";
        statfile_best_.open(iprefix,std::ofstream::trunc);
    } else if (is_indiv_) {
      std::string file = "stats/stats_simd_best";
      #ifdef HAVE_MPI
      file+="_"+std::to_string(rank);
      #endif
      file+=".csv";
      std::string fileTmp = "stats/stats_simd_best";
      #ifdef HAVE_MPI
      fileTmp+="_"+std::to_string(rank);
      #endif
      fileTmp+=".csv.tmp";

      tmp_best.open(file,std::ifstream::in);
      statfile_best_.open(fileTmp, std::ofstream::trunc);
    } else {
      std::string file = "stats/stats_simd_mean";
      #ifdef HAVE_MPI
      file+="_"+std::to_string(rank);
      #endif
      file+=".csv";
      std::string fileTmp = "stats/stats_simd_mean";
      #ifdef HAVE_MPI
      fileTmp+="_"+std::to_string(rank);
      #endif
      fileTmp+=".csv.tmp";

      tmp_mean.open(file,std::ifstream::in);
      statfile_mean_.open(fileTmp, std::ofstream::trunc);
    }


    if (prefix == nullptr) {  
      std::string str;
      for (int i = 0; i <= generation_; i++) {
        if (is_indiv_) {
          std::getline(tmp_best, str);
          statfile_best_ << str << std::endl;
        } else {
          std::getline(tmp_mean, str);
          statfile_mean_ << str << std::endl;
        }
      }

      if (is_indiv_) {
        statfile_best_.flush();
        statfile_best_.close();
      } else {
        statfile_mean_.flush();
        statfile_mean_.close();
      }

      if (is_indiv_) {
        std::string file = "stats/stats_simd_best";
        #ifdef HAVE_MPI
        file+="_"+std::to_string(rank);
        #endif
        file+=".csv";

        std::string fileTmp = "stats/stats_simd_best";
        #ifdef HAVE_MPI
        fileTmp+="_"+std::to_string(rank);
        #endif
        fileTmp+=".csv.tmp";

        statfile_best_.open(file, std::ofstream::trunc);
          tmp_best.close();
        tmp_best.open(fileTmp, std::ifstream::in);
          tmp_best.seekg(0, std::ios::beg);
      } else {
        std::string file = "stats/stats_simd_mean";
        #ifdef HAVE_MPI
        file+="_"+std::to_string(rank);
        #endif
        file+=".csv";

        std::string fileTmp = "stats/stats_simd_mean";
        #ifdef HAVE_MPI
        fileTmp+="_"+std::to_string(rank);
        #endif
        fileTmp+=".csv.tmp";

        statfile_mean_.open(file, std::ofstream::trunc);
          tmp_mean.close();
        tmp_mean.open(fileTmp, std::ifstream::in);
          tmp_mean.seekg(0, std::ios::beg);
      }
      
      
      for (int i = 0; i <= generation_; i++) {
        if (is_indiv_) {
          std::getline(tmp_best, str);
          if (str!="")
            statfile_best_ << str << std::endl;
        } else {
          std::getline(tmp_mean, str);
          if (str!="")
            statfile_mean_ << str << std::endl;
        }
      }

        tmp_best.close();
        tmp_mean.close();

        if (is_indiv_) {
          statfile_best_.flush();
        } else {
          statfile_mean_.flush();
        }
    }
  }
}

void Stats_7::compute_best(bool non_coding) {
  compute(simd_individual_->best_indiv, non_coding);
}


void Stats_7::compute(Individual_7* indiv, bool non_coding) {

  if (non_coding) indiv->compute_non_coding();
  //  printf("Compute BEST\n");
  is_indiv_ = true;

  nb_clones_ = simd_individual_->nb_clones_;
  pop_size_ = simd_individual_->nb_indivs_;

  fitness_ = indiv->fitness;
  metabolic_error_  = indiv->metaerror;

  amount_of_dna_ = indiv->dna_->length();

  nb_coding_rnas_ = indiv->nb_coding_RNAs;
  nb_non_coding_rnas_ = indiv->nb_non_coding_RNAs;

  nb_functional_genes_ = indiv->nb_func_genes;
  nb_non_functional_genes_ = indiv->nb_non_func_genes;


  nb_mut_ = indiv->dna_->nb_mut_;
  nb_rear_ = indiv->dna_->nb_rear_;
  nb_switch_ = indiv->dna_->nb_swi_;
  nb_indels_ = indiv->dna_->nb_indels_;
  nb_dupl_ = indiv->dna_->nb_large_dupl_;
  nb_del_ = indiv->dna_->nb_large_del_;
  nb_trans_ = indiv->dna_->nb_large_trans_;
  nb_inv_ = indiv->dna_->nb_large_inv_;

  dupl_rate_  = nb_dupl_  / (double)( indiv->dna_->parent_length());
  del_rate_   = nb_del_   / (double)( indiv->dna_->parent_length());
  trans_rate_ = nb_trans_ / (double)( indiv->dna_->parent_length());
  inv_rate_   = nb_inv_   / (double)( indiv->dna_->parent_length());

#ifdef __REGUL
  int32_t nb_activators = 0;
  int32_t nb_operators = 0;
  double mean_activator_activity = 0.0;
  double mean_operator_activity = 0.0;

  indiv->metadata_->rna_begin();
  for (int i = 0; i < indiv->metadata_->rna_count(); i++) {
    Rna_7*rna = indiv->metadata_->rna_next();
    if (rna != nullptr) {
      if (rna->is_coding_) {
        for (auto affinity: rna->affinity_list) {
          if (affinity.enhancer_factor > 0.0)
          {
            nb_activators++;
            mean_activator_activity += affinity.enhancer_factor;
          }

          if (affinity.operator_factor > 0.0)
          {
            nb_operators++;
            mean_operator_activity += affinity.operator_factor;
          }
        }
      }
    }
  }


  nb_enhancing_influences_       = nb_activators;
  nb_operating_influences_       = nb_operators;
  nb_influences_                 = nb_operating_influences_ + nb_enhancing_influences_;
  av_value_influences_           = ( mean_activator_activity + mean_operator_activity ) / double ( nb_activators + nb_operators);
  av_value_enhancing_influences_ = ( mean_activator_activity ) / double ( nb_activators );
  av_value_operating_influences_ = ( mean_operator_activity ) / double ( nb_operators);

  int32_t nb_TF = 0;
  int32_t nb_pure_TF = 0;

  indiv->metadata_->protein_begin();
  for (int j = 0; j < indiv->metadata_->proteins_count(); j++) {
    Protein_7* prot = indiv->metadata_->protein_next();
    if (prot != nullptr) {
      if (prot->is_init_) {
        if(prot->is_functional)
        {
          if (prot->is_TF_)
          {
            nb_TF+=1;
          }
        }
        else
        {
          if (prot->is_TF_)
          {
            nb_TF+=1;
            nb_pure_TF+=1;
          }
        }
      }
    }
  }

	nb_TF_ = nb_TF;
	nb_pure_TF_ = nb_pure_TF;
#endif

  is_computed_ = true;
}


void Stats_7::compute_average(bool non_coding) {
  is_indiv_ = false;
  nb_clones_ = simd_individual_->nb_clones_;
  pop_size_ = simd_individual_->nb_indivs_;
  

  for (int indiv_id = 0; indiv_id < pop_size_; indiv_id++) {

    if (non_coding) simd_individual_->previous_individuals[indiv_id]->compute_non_coding();

    fitness_ += simd_individual_->previous_individuals[indiv_id]->fitness;
    metabolic_error_ += simd_individual_->previous_individuals[indiv_id]->metaerror;

    amount_of_dna_ += simd_individual_->previous_individuals[indiv_id]->dna_->length();

    nb_coding_rnas_ += simd_individual_->previous_individuals[indiv_id]->nb_coding_RNAs;
    nb_non_coding_rnas_ += simd_individual_->previous_individuals[indiv_id]->nb_non_coding_RNAs;

    nb_functional_genes_ += simd_individual_->previous_individuals[indiv_id]->nb_func_genes;
    nb_non_functional_genes_ += simd_individual_->previous_individuals[indiv_id]->nb_non_func_genes;

    nb_mut_ += simd_individual_->previous_individuals[indiv_id]->dna_->nb_mut_;
    nb_rear_ += simd_individual_->previous_individuals[indiv_id]->dna_->nb_rear_;
    nb_switch_ += simd_individual_->previous_individuals[indiv_id]->dna_->nb_swi_;
    nb_indels_ += simd_individual_->previous_individuals[indiv_id]->dna_->nb_indels_;
    nb_dupl_ += simd_individual_->previous_individuals[indiv_id]->dna_->nb_large_dupl_;
    nb_del_ += simd_individual_->previous_individuals[indiv_id]->dna_->nb_large_del_;
    nb_trans_ += simd_individual_->previous_individuals[indiv_id]->dna_->nb_large_trans_;
    nb_inv_ += simd_individual_->previous_individuals[indiv_id]->dna_->nb_large_inv_;

    dupl_rate_ += nb_dupl_ / simd_individual_->previous_individuals[indiv_id]->dna_->parent_length();
    del_rate_ += nb_del_ / simd_individual_->previous_individuals[indiv_id]->dna_->parent_length();
    trans_rate_ +=
        nb_trans_ / simd_individual_->previous_individuals[indiv_id]->dna_->parent_length();
    inv_rate_ += nb_inv_ / simd_individual_->previous_individuals[indiv_id]->dna_->parent_length();

  }

  fitness_ /= pop_size_;
  metabolic_error_ /= pop_size_;

  amount_of_dna_ /= pop_size_;
  nb_coding_rnas_ /= pop_size_;
  nb_non_coding_rnas_ /= pop_size_;

  nb_functional_genes_ /= pop_size_;
  nb_non_functional_genes_ /= pop_size_;

  nb_mut_ /= pop_size_;
  nb_rear_ /= pop_size_;
  nb_switch_ /= pop_size_;
  nb_indels_ /= pop_size_;
  nb_dupl_ /= pop_size_;
  nb_del_ /= pop_size_;
  nb_trans_ /= pop_size_;
  nb_inv_ /= pop_size_;

  dupl_rate_ /= pop_size_;
  del_rate_ /= pop_size_;
  trans_rate_ /= pop_size_;
  inv_rate_ /= pop_size_;

  is_computed_ = true;
}

void Stats_7::write_best(bool non_coding) {
  write(simd_individual_->best_indiv,non_coding);
}

void Stats_7::write(Individual_7* indiv,bool non_coding) {
  if (is_indiv_ && !is_computed_)
    compute(indiv,non_coding);



  if (is_indiv_ && is_computed_) {
    // Write best stats
    statfile_best_<<generation_;
    statfile_best_<<","<<pop_size_<<","<<nb_clones_;
    statfile_best_<<","<<fitness_;
    statfile_best_<<","<<metabolic_error_;
    statfile_best_<<","<<amount_of_dna_;
    statfile_best_<<","<<nb_coding_rnas_;
    statfile_best_<<","<<nb_non_coding_rnas_;
    statfile_best_<<","<<nb_functional_genes_;
    statfile_best_<<","<<nb_non_functional_genes_;
    statfile_best_<<","<<nb_mut_;
    statfile_best_<<","<<nb_switch_;
    statfile_best_<<","<<nb_indels_;
    statfile_best_<<","<<nb_rear_;
    statfile_best_<<","<<nb_dupl_;
    statfile_best_<<","<<nb_del_;
    statfile_best_<<","<<nb_trans_;
    statfile_best_<<","<<nb_inv_;
    statfile_best_<<","<<dupl_rate_;
    statfile_best_<<","<<del_rate_;
    statfile_best_<<","<<trans_rate_;
    statfile_best_<<","<<inv_rate_;
    #ifdef __REGUL
      statfile_best_<< "," << nb_influences_;
      statfile_best_<< "," << nb_enhancing_influences_;
      statfile_best_<< "," << nb_operating_influences_;
      statfile_best_<< "," << av_value_influences_;
      statfile_best_<< "," << av_value_enhancing_influences_;
      statfile_best_<< "," << av_value_operating_influences_;
      statfile_best_<< "," << nb_TF_;
      statfile_best_<< "," << nb_pure_TF_;
    #endif
    if (non_coding) {
      statfile_best_<< "," << indiv->nb_bases_in_0_CDS_;
      statfile_best_<< "," << indiv->nb_bases_in_0_functional_CDS_;
      statfile_best_<< "," << indiv->nb_bases_in_0_non_functional_CDS_;
      statfile_best_<< "," << indiv->nb_bases_in_0_RNA_;
      statfile_best_<< "," << indiv->nb_bases_in_0_coding_RNA_;
      statfile_best_<< "," << indiv->nb_bases_in_0_non_coding_RNA_;
      statfile_best_<< "," << indiv->nb_bases_non_essential_;
      statfile_best_<< "," << indiv->nb_bases_non_essential_including_nf_genes_;
    }
    statfile_best_<<std::endl;
    statfile_best_.flush();
  }
}

void Stats_7::write_average(bool non_coding) {
  if (!is_indiv_ && !is_computed_)
    compute_average();

  if (!is_indiv_ && is_computed_) {
    // Write average stats
    statfile_mean_<<generation_<<","<<pop_size_<<","<<nb_clones_<<","<<fitness_<<","<<metabolic_error_<<","<<
                  amount_of_dna_<<","<<nb_coding_rnas_<<","<<nb_non_coding_rnas_<<","<<
                  nb_functional_genes_<<","<<nb_non_functional_genes_<<","<<nb_mut_
                  <<","<<nb_switch_<<","<<nb_indels_<<","<<nb_rear_<<","<<nb_dupl_<<","<<
                  nb_del_<<","<<nb_trans_<<","<<nb_inv_<<","<<dupl_rate_<<","<<del_rate_
                  <<","<<trans_rate_<<","<<inv_rate_
                  <<std::endl;
    statfile_mean_.flush();
  }
}

    void Stats_7::reinit(int64_t generation) {
      generation_ = generation;

      pop_size_ = 0;

      fitness_ = 0;
      metabolic_error_ = 0;

      amount_of_dna_ = 0;
      nb_coding_rnas_ = 0;
      nb_non_coding_rnas_ = 0;

      nb_functional_genes_ = 0;
      nb_non_functional_genes_ = 0;

      nb_mut_ = 0;
      nb_rear_ = 0;
      nb_switch_ = 0;
      nb_indels_ = 0;
      nb_dupl_ = 0;
      nb_del_ = 0;
      nb_trans_ = 0;
      nb_inv_ = 0;


      nb_bases_in_0_CDS_ = 0;
      nb_bases_in_0_functional_CDS_ = 0;
      nb_bases_in_0_non_functional_CDS_ = 0;
      nb_bases_in_0_RNA_ = 0;
      nb_bases_in_0_coding_RNA_ = 0;
      nb_bases_in_0_non_coding_RNA_ = 0;

      nb_bases_non_essential_ = 0;
      nb_bases_non_essential_including_nf_genes_ = 0;

      is_computed_ = false;
    }

}
