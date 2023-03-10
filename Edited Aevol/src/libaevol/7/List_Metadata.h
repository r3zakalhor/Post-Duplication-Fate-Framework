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


#ifndef AEVOL_LIST_METADATA_H
#define AEVOL_LIST_METADATA_H


#include "Protein_7.h"
#include "Rna_7.h"
#include "Abstract_Metadata.h"
#include "Dna_7.h"
#include "Individual_7.h"


#include <algorithm>

using std::list;


namespace aevol {

    using Promoters1Strand_7     = std::list<PromoterStruct>;
    using Promoters2Strands_7    = std::vector<Promoters1Strand_7>;

    class List_Metadata : public Abstract_Metadata {
    public:
    // List_Metadata(Individual_7* indiv) { indiv_ = indiv; set_iterators(); };

     List_Metadata(Individual_7* indiv, List_Metadata* metadata) : Abstract_Metadata(indiv,metadata) {

            for (auto& strand: {LEADING, LAGGING}) {
                for (auto& rna: metadata->promoters_list_[strand]) {
                    promoters_list_[strand].emplace_back(rna);
                    promoters_list_[strand].back().to_compute_end_ = false;
                    #ifdef WITH_OPTIMIZE_DIFF_SEARCH
                          if (rna.rna != nullptr) {
                            int32_t glob_rna_idx = rna_count_++;
                            auto prom_it = promoters_list_[strand].end();
                            prom_it--;
                            PromoterStruct* promx = &*(prom_it);
                            Rna_7* rna_ad = new Rna_7(rna.rna,promx);
                            rna_add(glob_rna_idx,rna_ad,promx);
                            // rna_ad->prom = &(promoters_list_[strand].back());
                            // promoters_list_[strand].back().rna = rna_ad;
                            rna_ad->to_compute_start_pos_ = false;
                        }
                    #endif
                }
            }

            #ifdef __REGUL
            // printf("With heredity %d\n",indiv->exp_m_->exp_s()->get_with_heredity());

            if (indiv->exp_m_->exp_s()->get_with_heredity()) {
                inherited_proteins_.clear();

                // printf("---------- %d -------------- COPY PROTEINS -- begin\n",indiv->indiv_id);
                for (auto prot : metadata->proteins_) {
                    if (prot!=nullptr && prot->is_init_ && (!prot->signal_ && !prot->inherited_))
                        if (prot->e >
                            indiv->exp_m_->exp_s()->get_protein_presence_limit() && !(prot->signal_)) {

                                // if (prot->e > 1) {
                                //     printf("P_TO_CPY %f %f %f %f :: %d %d %d\n",prot->h,prot->w,prot->m,prot->e,prot->is_init_,prot->signal_,prot->inherited_);

                                //     printf("ERROR e is too HIGH\n");
                                //     exit(-22);
                                // }
                                    
                            Protein_7* inherited_prot = new Protein_7(prot,indiv->exp_m_);
                            inherited_prot->inherited_ = true;
                            inherited_proteins_.push_back(inherited_prot);
                        }
                }
                // printf("---------- %d -------------- COPY PROTEINS -- END\n",indiv->indiv_id);
                // printf("Size of inherited is %d (%d)\n",inherited_proteins_.size(), metadata->proteins_.size());
            }
            #endif
            // TODO: Add inherited protein if REGUL

            set_iterators();
        };

     List_Metadata(Individual_7* indiv,
                                int32_t* lead_prom_pos = nullptr, int8_t* lead_prom_error = nullptr, int32_t lead_prom_size = -1,
                                int32_t* lag_prom_pos = nullptr, int8_t* lag_prom_error = nullptr, int32_t lag_prom_size = -1) 
                                : Abstract_Metadata(indiv,lead_prom_pos,lead_prom_error,lead_prom_size,lag_prom_pos,lag_prom_error,lag_prom_size) {

            // printf("Leading/Lagging Length %d %d\n",lead_prom_size,lag_prom_size);
            for (int i = 0; i < lead_prom_size; i++) {
                
                PromoterStruct* nprom = new PromoterStruct(lead_prom_pos[i],lead_prom_error[i],
                                                   true);
                int prom_idx = -1;

                // prom_idx = indiv->metadata_->promoter_count();
                // indiv->metadata_->set_promoters_count(
                //     indiv->metadata_->promoter_count()+ 1);

                promoter_add(prom_idx, nprom);
                // printf("%d -- %d  -- %d -- Add promoter LEADING %d => %d\n",indiv->indiv_id,i,prom_idx,lead_prom_pos[i],lead_prom_error[i]);

                delete nprom;
            }

            for (int i = lag_prom_size-1; i >= 0; i--) {
                
                
                PromoterStruct* nprom = new PromoterStruct(lag_prom_pos[i],lag_prom_error[i],
                                                        false);
                int prom_idx = -1;
                // prom_idx = indiv->metadata_->promoter_count();
                // indiv->metadata_->set_promoters_count(
                //     indiv->metadata_->promoter_count() + 1);

                promoter_add(prom_idx, nprom);
                // printf("%d -- %d -- %d -- Add promoter LAGGING %d => %d\n",indiv->indiv_id,i,prom_idx,lag_prom_pos[i],lag_prom_error[i]);

                delete nprom;
            }

            set_iterators();
        };

        ~List_Metadata() override {

        //    printf("Delete NB RNA : %d\n",rnas_.size());
        //    int nb_del = 0;
           for (std::list<Rna_7*>::iterator it_rna = rnas_.begin(); it_rna != rnas_.end(); it_rna++) {
               delete (*(it_rna));
            //    nb_del++;
           }
           rnas_.clear();
        //    printf("Delete NB RNA : %d %d %d\n",rnas_.size(),nb_del,rna_count_add_);

            for (std::list<Protein_7*>::iterator it_protein = proteins_.begin(); it_protein != proteins_.end(); it_protein++) {
                // #ifdef __REGUL
                // if (! (*(it_protein))->inherited_)
                // #endif
                delete (*(it_protein));
            }

            // for (std::vector<Protein_7*>::iterator it_protein = signal_proteins_.begin(); it_protein != signal_proteins_.end(); it_protein++) {
            //     // #ifdef __REGUL
            //     // if (! (*(it_protein))->inherited_)
            //     // #endif
            //     delete (*(it_protein));
            // }

          

#ifdef __REGUL
        signal_proteins_.clear();

        if (indiv_->exp_m_->exp_s()->get_with_heredity()) {
          for (std::list<Protein_7*>::iterator it_protein = inherited_proteins_.begin(); it_protein != inherited_proteins_.end(); it_protein++) {
            delete (*(it_protein));
          }
        }

        //   for (std::vector<Protein_7*>::iterator it_protein = signal_proteins_.begin(); it_protein != signal_proteins_.end(); it_protein++) {
        //     delete (*(it_protein));
        //   }
#endif



            promoters_list_.clear();
           terminator_lead_.clear();
           terminator_lag_.clear();

        };

        void clean_remote() {
           terminator_lead_.clear();
           terminator_lag_.clear();

            for (std::list<Protein_7*>::iterator it_protein = proteins_.begin(); it_protein != proteins_.end(); it_protein++) {
                // #ifdef __REGUL
                // if (! (*(it_protein))->inherited_)
                // #endif
                delete (*(it_protein));
            }

            proteins_.clear();

#ifdef __REGUL
          for (std::list<Protein_7*>::iterator it_protein = inherited_proteins_.begin(); it_protein != inherited_proteins_.end(); it_protein++) {
            delete (*(it_protein));
          }

          inherited_proteins_.clear();

          signal_proteins_.clear();
#endif

           for (std::list<Rna_7*>::iterator it_rna = rnas_.begin(); it_rna != rnas_.end(); it_rna++) {
               delete (*(it_rna));
           }

           rnas_.clear();

           set_iterators();
        };

#ifdef __REGUL
        void add_inherited_proteins(Individual_7* indiv = nullptr);
#endif
        /** Getter **/

        void set_iterators() {
            it_promoter_ = promoters_list_[LEADING].begin();
            it_rna_ = rnas_.begin();
            it_protein_ = proteins_.begin();
        };

        /*** Promoters ***/
        PromoterStruct* promoters(int idx) override;
        void promoter_add(int idx, PromoterStruct* prom) override;

        PromoterStruct* promoter_next() override ;
        void promoter_begin() override ;
        bool promoter_end() override ;

        int promoter_count() override;
        void set_promoters_count(int pcount) override;

        /*** Terminators ***/
        int terminator_count(int LoL) override;
        void terminator_add(int LoL, int dna_pos) override;

        int next_terminator(int LoL, int dna_pos) override;

        void terminators_clear() override;

        /*** RNAs ***/
        Rna_7* rnas(int idx) override;
        void rna_add(int idx, Rna_7* rna, PromoterStruct* prom = nullptr) override; 
        void rna_add(int idx, int32_t t_begin, int32_t t_end,
                     int8_t t_leading_lagging, double t_e,
                     int32_t t_length, PromoterStruct* prom = nullptr) override ;

        Rna_7* rna_next() override ;
        void rna_begin() override ;
        bool rna_end() override ;

        int rna_count() override;
        void set_rna_count(int rcount) override;

        void rnas_resize(int resize) override;
        void rnas_clear() override;

        /*** Proteins ***/
        Protein_7* proteins(int idx) override;
        void protein_add(int idx, Protein_7* prot) override;
        void protein_add(int idx, int32_t t_protein_start,
            int32_t t_protein_end,
            int32_t t_protein_length,
            int8_t t_leading_lagging,
            double t_e, Rna_7* rna);

        Protein_7* protein_next() override ;
        void protein_begin() override ;
        bool protein_end() override ;

        int proteins_count() override;
        void set_proteins_count(int pcount) override;

        void proteins_resize(int resize) override;
        void proteins_clear() override;

        void proteins_remove_signal();

        void proteins_print(int step = -1, int indiv_id = -1);

        /*** Promoters ***/
        void lst_promoters(bool lorl,
                           Position before_after_btw, // with regard to the strand's reading direction
                           int32_t pos1,
                           int32_t pos2,
                           std::list<PromoterStruct*>&  motif_list) override;

        /** Search and update **/
        void remove_promoters_around(int32_t pos_1) override;
        void remove_promoters_around(int32_t pos_1, int32_t pos_2) override;
        void remove_all_promoters() override;

        void look_for_new_promoters_around(int32_t pos_1, int32_t pos_2) override;
        void look_for_new_promoters_around(int32_t pos) override;

        void locate_promoters() override;

        void move_all_promoters_after(int32_t pos, int32_t delta_pos) override;

        void duplicate_promoters_included_in(int32_t pos_1,
                                             int32_t pos_2,
                                             std::vector<std::list<PromoterStruct*>>& duplicated_promoters) override;
        void extract_promoters_included_in(int32_t pos_1,
                                           int32_t pos_2, std::vector<std::list<PromoterStruct*>>& extracted_promoters) override;
        void insert_promoters(std::vector<std::list<PromoterStruct*>>& promoters_to_insert) override;
        void insert_promoters_at(std::vector<std::list<PromoterStruct*>>& promoters_to_insert,
                                 int32_t pos) override;

        void invert_promoters_included_in(int32_t pos1,
                                          int32_t pos2) override;


        static void shift_promoters(
                std::vector<std::list<PromoterStruct*>>& promoters_to_shift,
                int32_t delta_pos,
                int32_t seq_length);
        static void invert_promoters(std::vector<std::list<PromoterStruct*>>& promoter_lists,
                                     int32_t pos1,
                                     int32_t pos2);

        void remove_leading_promoters_starting_between(int32_t pos_1,
                                                       int32_t pos_2) override;
        void remove_leading_promoters_starting_after(int32_t pos) override;
        void remove_leading_promoters_starting_before(int32_t pos) override;

        void remove_lagging_promoters_starting_between(int32_t pos_1,
                                                       int32_t pos_2) override;
        void remove_lagging_promoters_starting_after(int32_t pos) override;
        void remove_lagging_promoters_starting_before(int32_t pos) override;

        void move_all_leading_promoters_after(int32_t pos, int32_t delta_pos) override;
        void move_all_lagging_promoters_after(int32_t pos,int32_t delta_pos) override;

        void look_for_new_leading_promoters_starting_between(int32_t pos_1, int32_t pos_2) override;
        void look_for_new_leading_promoters_starting_after(int32_t pos) override;
        void look_for_new_leading_promoters_starting_before(int32_t pos) override;

        void look_for_new_lagging_promoters_starting_between(int32_t pos_1,int32_t pos_2) override;
        void look_for_new_lagging_promoters_starting_after(int32_t pos) override;
        void look_for_new_lagging_promoters_starting_before(int32_t pos) override;

        void promoters_included_in(int32_t pos_1,
                                   int32_t pos_2,
                                   std::vector<std::list<PromoterStruct*>>& promoters_list) override;

        void extract_leading_promoters_starting_between(int32_t pos_1,
                                                        int32_t pos_2, std::list<PromoterStruct*>& extracted_promoters) override;

        void extract_lagging_promoters_starting_between(int32_t pos_1,
                                                        int32_t pos_2,
                                                        std::list<PromoterStruct*>& extracted_promoters) override;

    #ifdef WITH_OPTIMIZE_DIFF_SEARCH

        void reinit_rnas(int32_t pos);

        void reinit_rnas(int32_t pos_1, int32_t pos_2) {
            // if (indiv_->indiv_id==6)
            //     printf("Init RNA %d %d\n",pos_1,pos_2);
            for (auto i = pos_1; i <= pos_2;i++) {
                reinit_rnas(i);
                reinit_rnas(i);
            }
        }
    #endif

    #ifdef __REGUL
        void reinit_proteins_inherited(Individual_7* indiv) {
            //proteins_remove_signal();
            
            if (indiv->exp_m_->exp_s()->get_with_heredity()) {
                for (std::list<Protein_7*>::iterator it_protein = inherited_proteins_.begin(); it_protein != inherited_proteins_.end(); it_protein++) {
                    delete (*(it_protein));
                }
                inherited_proteins_.clear();

                for (auto prot : proteins_) {
                    if (prot!=nullptr && prot->is_init_ && prot->translated_)
                        if (prot->e > indiv->exp_m_->exp_s()->get_protein_presence_limit()) {
                            Protein_7* inherited_prot = new Protein_7(prot,indiv->exp_m_);
                            inherited_prot->inherited_ = true;
                            inherited_proteins_.push_back(inherited_prot);
                        }
                }

                for (std::list<Protein_7*>::iterator it_protein = proteins_.begin(); it_protein != proteins_.end(); it_protein++) {
                    if (!(*(it_protein))->translated_)
                        delete (*(it_protein));
                }
                proteins_.clear();

                for (auto& prot : proteins_translated_) {
                    if (prot!=nullptr && prot->is_init_) {
		        for (auto&& rna: prot->rna_list_) {
			   rna->affinity_list.clear();
			}
                        proteins_.push_back(prot);
                    }
                }


                rna_begin();
                for (int i = 0; i < rna_count(); i++) {
                    Rna_7* rna = rna_next();
                    if (rna != nullptr) {
                        if (rna->is_coding_) {
                            rna->affinity_list.clear();
                        }
                    }
                }
            }
        }
    #endif

        Promoters2Strands_7 promoters_list_ = {{},
                                                  {}};

        std::list<Protein_7*> proteins_;
#ifdef __REGUL
  std::list<Protein_7*> inherited_proteins_;
  std::vector<Protein_7*> signal_proteins_;
  std::list<Protein_7*> proteins_translated_;
  bool inherited_added_ = false;
#endif
        std::list<Rna_7*> rnas_;

        // Individual_7* indiv_;
 protected:
     Promoters1Strand_7::iterator it_promoter_;
        int it_promoter_pos_;

        std::list<Rna_7*>::iterator it_rna_;
        std::list<Protein_7*>::iterator it_protein_;

        std::set<int> terminator_lag_;
        std::set<int> terminator_lead_;


        int32_t protein_count_ = 0;
        int32_t rna_count_add_ = 0;

        int cmp_rna = 0;

    };
}


#endif //AEVOL_LIST_METADATA_H
