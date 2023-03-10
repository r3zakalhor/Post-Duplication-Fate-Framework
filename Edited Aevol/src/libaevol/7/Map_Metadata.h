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


#ifndef AEVOL_MAP_METADATA_H
#define AEVOL_MAP_METADATA_H


#include "Rna_7.h"
#include "Protein_7.h"
#include "Abstract_Metadata.h"
#include "Dna_7.h"
#include "Individual_7.h"

namespace aevol {
    class Map_Metadata : public Abstract_Metadata {
    public:
     Map_Metadata(Individual_7* indiv) : Abstract_Metadata(indiv) { count_promoters_ = 0; };

     Map_Metadata(Individual_7* indiv, Map_Metadata* metadata) : Abstract_Metadata(indiv,metadata) {
            count_promoters_ = 0;

            for (const auto& prom : metadata->promoters_) {
                if (prom.second != nullptr) {
                    auto prom_copy = new PromoterStruct(prom.second->pos, prom.second->error,
                                                        prom.second->leading_or_lagging);
                    promoters_[count_promoters_] = prom_copy;


                    if (prom.second->leading_or_lagging) {
                        leading_prom_pos_[prom_copy->pos] = count_promoters_;
                    } else {
                        lagging_prom_pos_[prom_copy->pos] = count_promoters_;
                    }

                    count_promoters_++;
                }
            }
        };

        ~Map_Metadata() override {
            if (! promoters_.empty()) {
                for (auto element = promoters_.begin();
                     element != promoters_.end(); ++element) {
                    if (element->second != nullptr) delete element->second;
                }
            }

            promoters_.clear();

            leading_prom_pos_.clear();
            lagging_prom_pos_.clear();


            for (auto rn : rnas_) {
                delete rn;
            }

            for (auto prot : proteins_) {
                delete prot;
            }

            rnas_.clear();
            proteins_.clear();
        };

        /** Getter **/

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
        void rna_add(int idx, Rna_7* rna) override;
        void rna_add(int idx, int32_t t_begin, int32_t t_end,
                     int8_t t_leading_lagging, double t_e,
                      int32_t t_length) override {exit(1);};

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

        Protein_7* protein_next() override ;
        void protein_begin() override ;
        bool protein_end() override ;

        int proteins_count() override;
        void set_proteins_count(int pcount) override;

        void proteins_resize(int resize) override;
        void proteins_clear() override;

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


        void rebuild_index() {
            if (count_promoters_ > (int)promoters_.size()/2) {
                /**
                 * Do the reindexation process
                 */
                auto old_promoters = promoters_;
                promoters_.clear();
                leading_prom_pos_.clear();
                lagging_prom_pos_.clear();
                count_promoters_ = 0;
                for (auto prom : old_promoters) {
                    promoters_[count_promoters_] = prom.second;
                    if (prom.second->leading_or_lagging) {
                        leading_prom_pos_[prom.second->pos] = count_promoters_;
                    } else {
                        lagging_prom_pos_[prom.second->pos] = count_promoters_;
                    }
                    count_promoters_++;
                }
            }
        }

    protected:
        std::map<int32_t, PromoterStruct*> promoters_;
        std::map<int32_t,int32_t> leading_prom_pos_;
        std::map<int32_t,int32_t> lagging_prom_pos_;

        std::set<int> terminator_lag_;
        std::set<int> terminator_lead_;
        std::vector<Rna_7*> rnas_;
        std::vector<Protein_7*> proteins_;

        int it_promoter_ = 0;
        int it_rna_ = 0;
        int it_protein_ = 0;

        int32_t count_promoters_;
        int32_t protein_count_ = 0;
        Individual_7* indiv_;
    };
}


#endif
