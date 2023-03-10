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


#include "DynTab_Metadata.h"

#include "Rna_7.h"
#include "AeTime.h"

namespace aevol {




    void DynTab_Metadata::lst_promoters(bool lorl, Position before_after_btw, int32_t pos1, int32_t pos2,
                                          std::list<PromoterStruct*>& motif_list) {

        /*
         * pos_1    >=      LEADING
         * pos_1    <       LAGGING
         *
         * pos_2    >=      LEADING
         * pos_2    <       LAGGING
         */

        for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
            if (promoters_[prom_idx] != nullptr) {


                if (before_after_btw == BEFORE) {
                    //  begin -> pos_2
                    if (promoters_[prom_idx]->leading_or_lagging && lorl == LEADING) {
                        if (promoters_[prom_idx]->pos < pos2) {
                            motif_list.push_back(promoters_[prom_idx]);
                        }
                    } else if (!promoters_[prom_idx]->leading_or_lagging && lorl == LAGGING)
                        if (promoters_[prom_idx]->pos >= pos2) {
                            motif_list.push_back(promoters_[prom_idx]);
                        }
                } else if (before_after_btw == AFTER) {
                    // pos_1 -> end
                    if (promoters_[prom_idx]->leading_or_lagging && lorl == LEADING) {
                        if (promoters_[prom_idx]->pos >= pos1) {
                            motif_list.push_back(promoters_[prom_idx]);
                        }
                    } else if (!promoters_[prom_idx]->leading_or_lagging && lorl == LAGGING)
                        if (promoters_[prom_idx]->pos < pos1) {
                            motif_list.push_back(promoters_[prom_idx]);
                        }
                } else {
                    // pos_1 -> pos_2
                    if (promoters_[prom_idx]->leading_or_lagging && lorl == LEADING) {
                        if ((promoters_[prom_idx]->pos >= pos1) && (promoters_[prom_idx]->pos < pos2)) {
                            motif_list.push_back(promoters_[prom_idx]);
                        }
                    } else if (!promoters_[prom_idx]->leading_or_lagging && lorl == LAGGING)
                        if ((promoters_[prom_idx]->pos < pos1) && (promoters_[prom_idx]->pos >= pos2)) {
                            motif_list.push_back(promoters_[prom_idx]);
                        }
                }
            }
        }

    }

    void DynTab_Metadata::remove_promoters_around(int32_t pos) {
        if (length() >= PROM_SIZE) {
            remove_leading_promoters_starting_between(Utils::mod(pos - PROM_SIZE + 1,
                                                                 length()),
                                                      pos);
            remove_lagging_promoters_starting_between(pos,
                                                      Utils::mod(pos + PROM_SIZE - 1,
                                                                 length()));
        }
        else {
            remove_all_promoters();
        }
    }

    void DynTab_Metadata::remove_promoters_around(int32_t pos_1, int32_t pos_2) {
        if (Utils::mod(pos_1 - pos_2, length()) >= PROM_SIZE) {
            remove_leading_promoters_starting_between(Utils::mod(pos_1 - PROM_SIZE + 1,
                                                                 length()),
                                                      pos_2);
            remove_lagging_promoters_starting_between(pos_1,
                                                      Utils::mod(pos_2 + PROM_SIZE - 1,
                                                                 length()));
        }
        else {
            remove_all_promoters();
        }
    }

    void DynTab_Metadata::remove_all_promoters() {

        for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
            delete promoters_[prom_idx];
            promoters_[prom_idx] = nullptr;
        }

        count_promoters_ = 0;
    }

    void DynTab_Metadata::look_for_new_promoters_around(int32_t pos_1, int32_t pos_2) {
        if (length() >= PROM_SIZE) {
            look_for_new_leading_promoters_starting_between(
                    Utils::mod(pos_1 - PROM_SIZE + 1,
                               length()), pos_2);

            look_for_new_lagging_promoters_starting_between(pos_1, Utils::mod(
                    pos_2 + PROM_SIZE - 1,
                    length()));
        }
    }

    void DynTab_Metadata::look_for_new_promoters_around(int32_t pos) {
        if (length() >= PROM_SIZE) {
            look_for_new_leading_promoters_starting_between(
                    Utils::mod(pos - PROM_SIZE + 1, length()),
                    pos);
            look_for_new_lagging_promoters_starting_between(
                    pos,
                    Utils::mod(pos + PROM_SIZE - 1, length()));
        }
    }

    void DynTab_Metadata::locate_promoters() {
        look_for_new_leading_promoters_starting_between(0,length());
        look_for_new_lagging_promoters_starting_between(0,length());
    }

    void DynTab_Metadata::move_all_promoters_after(int32_t pos, int32_t delta_pos) {
        move_all_leading_promoters_after(pos, delta_pos);
        move_all_lagging_promoters_after(pos, delta_pos);
    }

    void DynTab_Metadata::duplicate_promoters_included_in(int32_t pos_1, int32_t pos_2,
                                                            std::vector<std::list<PromoterStruct*>> &duplicated_promoters) {
        // 1) Get promoters to be duplicated
        std::vector<std::list<PromoterStruct*>> retrieved_promoters = {{},
                                                                       {}};

        promoters_included_in(pos_1, pos_2, retrieved_promoters);

        // 2) Set RNAs' position as their position on the duplicated segment
        for (auto& strand: {LEADING, LAGGING}) {
            for (auto& prom : retrieved_promoters[strand]) {
                // Make a copy of current RNA inside container
                duplicated_promoters[strand].push_back(new PromoterStruct(prom));

                // Set RNA's position as it's position on the duplicated segment
                duplicated_promoters[strand].back()->pos = Utils::mod(duplicated_promoters[strand].back()->pos -pos_1,
                                                                      length());
            }
        }

    }

    void DynTab_Metadata::extract_promoters_included_in(int32_t pos_1, int32_t pos_2,
                                                          std::vector<std::list<PromoterStruct*>> &extracted_promoters) {
        if (pos_2 - pos_1 < PROM_SIZE) {
            return;
        }

        extract_leading_promoters_starting_between(pos_1, pos_2 - PROM_SIZE + 1,
                                                   extracted_promoters[LEADING]);
        extract_lagging_promoters_starting_between(pos_1 + PROM_SIZE - 1, pos_2,
                                                   extracted_promoters[LAGGING]);
    }

    void DynTab_Metadata::insert_promoters(std::vector<std::list<PromoterStruct*>> &promoters_to_insert) {
        if (count_promoters_ + promoters_to_insert[LEADING].size() + promoters_to_insert[LAGGING].size() >= dyntab_size_)
            reallocate_promoters();

        for (auto strand: {LEADING, LAGGING}) {
            if (promoters_to_insert[strand].empty()) {
                continue;
            }
            // Insert the promoters in the individual's RNA list
            for (auto& to_insert: promoters_to_insert[strand]) {

                bool to_add = true;
                if (strand == LEADING) {
                    for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                        if (promoters_[prom_idx] != nullptr) {
                            if ((promoters_[prom_idx]->pos == to_insert->pos) &&
                                (promoters_[prom_idx]->leading_or_lagging == to_insert->leading_or_lagging)) {
                                to_add = false;
                                break;
                            }
                        }
                    }
                } else {
                    for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                        if (promoters_[prom_idx] != nullptr) {
                            if ((promoters_[prom_idx]->pos == to_insert->pos) &&
                                (promoters_[prom_idx]->leading_or_lagging == to_insert->leading_or_lagging)) {
                                to_add = false;
                                break;
                            }
                        }
                    }
                }
                if (to_add) {
                    int prom_idx;
                    {
                        prom_idx = count_promoters_;
                        count_promoters_ = count_promoters_ + 1;
                    }
                    promoters_[prom_idx] = to_insert;
                }

            }
        }
    }

    void DynTab_Metadata::insert_promoters_at(std::vector<std::list<PromoterStruct*>>& promoters_to_insert,
                             int32_t pos) {
        if (count_promoters_ + promoters_to_insert[LEADING].size() + promoters_to_insert[LAGGING].size() >= dyntab_size_)
            reallocate_promoters();

        for (auto strand: {LEADING, LAGGING}) {
            if (promoters_to_insert[strand].size() <= 0) {
                continue;
            }

            // Insert the promoters in the individual's RNA list
            for (auto &to_insert: promoters_to_insert[strand]) {
                // Update promoter position
                to_insert->pos = Utils::mod(to_insert->pos + pos, length());
                bool to_add = true;
                if (strand == LEADING) {
                    for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                        if (promoters_[prom_idx] != nullptr) {
                            if ((promoters_[prom_idx]->pos == to_insert->pos) &&
                                (promoters_[prom_idx]->leading_or_lagging == to_insert->leading_or_lagging)) {
                                to_add = false;
                                break;
                            }
                        }
                    }
                } else {
                    for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                        if (promoters_[prom_idx] != nullptr) {
                            if ((promoters_[prom_idx]->pos == to_insert->pos) &&
                                (promoters_[prom_idx]->leading_or_lagging == to_insert->leading_or_lagging)) {
                                to_add = false;
                                break;
                            }
                        }
                    }
                }
                if (to_add) {
                    int prom_idx;
                    {
                        prom_idx = count_promoters_;
                        count_promoters_ = count_promoters_ + 1;
                    }
                    promoters_[prom_idx] = to_insert;
                }

            }
        }
    }

    void DynTab_Metadata::invert_promoters_included_in(int32_t pos1,
                                      int32_t pos2) {
        int32_t segment_length = pos2 - pos1;

        if (segment_length < PROM_SIZE) {
            return;
        }

        std::vector<std::list<PromoterStruct*>> inverted_promoters = {{},
                                                                      {}};

        // 1) Extract the promoters completely included on the segment to be inverted
        extract_promoters_included_in(pos1, pos2, inverted_promoters);

        // 2) Invert segment's promoters
        DynTab_Metadata::invert_promoters(inverted_promoters, pos1, pos2);

        // 3) Reinsert the inverted promoters
        insert_promoters(inverted_promoters);
    }

    void DynTab_Metadata::shift_promoters(
            std::vector<std::list<PromoterStruct*>>& promoters_to_shift,
            int32_t delta_pos,
            int32_t seq_length) {
        for (auto& strand: {LEADING, LAGGING})
            for (auto& prom: promoters_to_shift[strand])
                prom->pos = Utils::mod(prom->pos + delta_pos, seq_length);
    }

    void DynTab_Metadata::invert_promoters(std::vector<std::list<PromoterStruct*>>& promoter_lists,
                                         int32_t pos1,
                                         int32_t pos2) {
        // Exchange LEADING and LAGGING lists
        promoter_lists[LEADING].swap(promoter_lists[LAGGING]);

        // Update the position and strand of each promoter to be inverted...
        for (auto& strand: {LEADING, LAGGING})
            for (auto& prom: promoter_lists[strand]) {
                prom->pos = pos1 + pos2 - prom->pos - 1;
                prom->leading_or_lagging = !prom->leading_or_lagging;
            }
    }

    void DynTab_Metadata::remove_leading_promoters_starting_between(int32_t pos_1,
                                                   int32_t pos_2) {
        if (pos_1 > pos_2) {
            remove_leading_promoters_starting_after(pos_1);
            remove_leading_promoters_starting_before(pos_2);
        }
        else {
            for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                if (promoters_[prom_idx] != nullptr) {
                    if (promoters_[prom_idx]->leading_or_lagging)
                        if ((promoters_[prom_idx]->pos >= pos_1) && (promoters_[prom_idx]->pos < pos_2)) {
                            delete promoters_[prom_idx];
                            promoters_[prom_idx] = nullptr;
                        }
                }
            }
        }

    }
    void DynTab_Metadata::remove_leading_promoters_starting_after(int32_t pos) {
        for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
            if (promoters_[prom_idx] != nullptr) {
                if (promoters_[prom_idx]->leading_or_lagging)
                    if (promoters_[prom_idx]->pos >= pos) {
                        delete promoters_[prom_idx];
                        promoters_[prom_idx] = nullptr;
                    }
            }
        }
    }

    void
    DynTab_Metadata::remove_leading_promoters_starting_before(int32_t pos) {
        for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
            if (promoters_[prom_idx] != nullptr) {
                if (promoters_[prom_idx]->leading_or_lagging)
                    if (promoters_[prom_idx]->pos < pos) {
                        delete promoters_[prom_idx];
                        promoters_[prom_idx] = nullptr;
                    }
            }
        }
    }

    void DynTab_Metadata::remove_lagging_promoters_starting_between(int32_t pos_1,
                                                   int32_t pos_2) {
        if (pos_1 == length()) pos_1 = 0;
        if (pos_2 == 0) pos_2 = length();

        if (pos_1 >
            pos_2) { // vld: that's a weird case... really do this? used from remove_promoters_around()
            remove_lagging_promoters_starting_after(pos_1);
            remove_lagging_promoters_starting_before(pos_2);
        } else {

            // Delete RNAs until we pass pos_1 (or we reach the end of the list)
            for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                if (promoters_[prom_idx] != nullptr) {
                    if (!promoters_[prom_idx]->leading_or_lagging)
                        if ((promoters_[prom_idx]->pos >= pos_1) && (promoters_[prom_idx]->pos < pos_2)) {
                            delete promoters_[prom_idx];
                            promoters_[prom_idx] = nullptr;
                        }
                }
            }
        }
    }

    void DynTab_Metadata::remove_lagging_promoters_starting_after(int32_t pos) {
        for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
            if (promoters_[prom_idx] != nullptr) {
                if (!promoters_[prom_idx]->leading_or_lagging)
                    if (promoters_[prom_idx]->pos >= pos) {
                        delete promoters_[prom_idx];
                        promoters_[prom_idx] = nullptr;
                    }
            }
        }
    }

    void
    DynTab_Metadata::remove_lagging_promoters_starting_before(int32_t pos) {
        // Delete RNAs until we reach pos (or we reach the end of the list)
        // TODO: optimize by starting from the end (with reverse iterators)
        for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
            if (promoters_[prom_idx] != nullptr) {
                if (!promoters_[prom_idx]->leading_or_lagging)
                    if (promoters_[prom_idx]->pos < pos) {
                        delete promoters_[prom_idx];
                        promoters_[prom_idx] = nullptr;
                    }
            }
        }
    }

    void DynTab_Metadata::move_all_leading_promoters_after(int32_t pos, int32_t delta_pos) {
        for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
            if (promoters_[prom_idx] != nullptr) {
                if (promoters_[prom_idx]->leading_or_lagging)
                    if (promoters_[prom_idx]->pos >= pos) {
                        promoters_[prom_idx]->pos = Utils::mod(promoters_[prom_idx]->pos + delta_pos, length());
                    }
            }
        }
    }

    void DynTab_Metadata::move_all_lagging_promoters_after(int32_t pos,int32_t delta_pos) {
        for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
            if (promoters_[prom_idx] != nullptr) {
                if (!promoters_[prom_idx]->leading_or_lagging)
                    if (promoters_[prom_idx]->pos >= pos) {
                        promoters_[prom_idx]->pos = Utils::mod(promoters_[prom_idx]->pos + delta_pos, length());
                    }
            }
        }
    }

    void DynTab_Metadata::look_for_new_leading_promoters_starting_between(int32_t pos_1, int32_t pos_2) {
        // When pos_1 > pos_2, we will perform the search in 2 steps.
        // As positions  0 and dna_->length() are equivalent, it's preferable to
        // keep 0 for pos_1 and dna_->length() for pos_2.


        if (pos_1 >= pos_2) {
            if (count_promoters_ + (pos_1 + (length() - pos_2))/PROM_SIZE >= dyntab_size_)
                reallocate_promoters();

            look_for_new_leading_promoters_starting_after(pos_1);
            look_for_new_leading_promoters_starting_before(pos_2);
            return;
        }
        // Hamming distance of the sequence from the promoter consensus

        if (count_promoters_ + (pos_2 - pos_1) / PROM_SIZE >= dyntab_size_)
            reallocate_promoters();

        for (int32_t i = pos_1; i < pos_2; i++) {

            int8_t dist = is_promoter_leading(i);

            if (dist <= 4) {
                bool notFound = true;
                for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                    if (promoters_[prom_idx] != nullptr) {
                        if (promoters_[prom_idx]->leading_or_lagging)
                            if (promoters_[prom_idx]->pos == i) {
                                notFound = false;
                            }
                    }
                }

                if (notFound) {
                    int prom_idx;
                    {
                        prom_idx = count_promoters_;
                        count_promoters_ = count_promoters_ + 1;
                    }

                    promoters_[prom_idx] = new PromoterStruct(i, dist, true);
                }

            }
        }
    }

    void DynTab_Metadata::look_for_new_leading_promoters_starting_after(int32_t pos) {
        // Hamming distance of the sequence from the promoter consensus
        //int8_t dist = 8;

        if (count_promoters_ + (length() - pos)/PROM_SIZE >= dyntab_size_)
            reallocate_promoters();

        for (int32_t i = pos; i < length(); i++) {

            int8_t dist = is_promoter_leading(i);
            if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
                bool notFound = true;
                for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                    if (promoters_[prom_idx] != nullptr) {
                        if (promoters_[prom_idx]->leading_or_lagging)
                            if (promoters_[prom_idx]->pos == i) {
                                notFound = false;
                            }
                    }
                }

                if (notFound) {
                    int prom_idx;
                    {
                        prom_idx = count_promoters_;
                        count_promoters_ = count_promoters_ + 1;
                    }

                    promoters_[prom_idx] = new PromoterStruct(i, dist, true);;
                }
            }
        }
    }

    void DynTab_Metadata::look_for_new_leading_promoters_starting_before(int32_t pos) {
        // Hamming distance of the sequence from the promoter consensus
        if (count_promoters_ + (pos)/PROM_SIZE >= dyntab_size_)
            reallocate_promoters();

        for (int32_t i = 0; i < pos; i++) {

            int8_t dist = is_promoter_leading(i);
            if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
                bool notFound = true;
                for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                    if (promoters_[prom_idx] != nullptr) {
                        if (promoters_[prom_idx]->leading_or_lagging)
                            if (promoters_[prom_idx]->pos == i) {
                                notFound = false;
                            }
                    }
                }

                if (notFound) {
                    int prom_idx;
                    {
                        prom_idx = count_promoters_;
                        count_promoters_ = count_promoters_ + 1;
                    }

                    promoters_[prom_idx] = new PromoterStruct(i, dist, true);;
                }
            }
        }
    }

    void DynTab_Metadata::look_for_new_lagging_promoters_starting_between(int32_t pos_1,int32_t pos_2) {
        // When pos_1 > pos_2, we will perform the search in 2 steps.
        // As positions  0 and dna_->length() are equivalent, it's preferable to
        // keep 0 for pos_1 and dna_->length() for pos_2.

        if (pos_1 >= pos_2) {
            if (count_promoters_ + (pos_1 + (length() - pos_2))/PROM_SIZE >= dyntab_size_)
                reallocate_promoters();

            look_for_new_lagging_promoters_starting_after(pos_1);
            look_for_new_lagging_promoters_starting_before(pos_2);
            return;
        }

        if (count_promoters_ + (pos_2 - pos_1) / PROM_SIZE >= dyntab_size_)
            reallocate_promoters();

        // Hamming distance of the sequence from the promoter consensus
        for (int32_t i = pos_2 - 1; i >= pos_1; i--) {

            int8_t dist = is_promoter_lagging(i);
            if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
                bool notFound = true;
                for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                    if (promoters_[prom_idx] != nullptr) {
                        if (!promoters_[prom_idx]->leading_or_lagging)
                            if (promoters_[prom_idx]->pos == i) {
                                notFound = false;
                            }
                    }
                }

                if (notFound) {
                    int prom_idx;
                    {
                        prom_idx = count_promoters_;
                        count_promoters_ = count_promoters_ + 1;
                    }

                    promoters_[prom_idx] = new PromoterStruct(i, dist, false);;
                }
            }
        }
    }

    void DynTab_Metadata::look_for_new_lagging_promoters_starting_after(int32_t pos) {
        // Hamming distance of the sequence from the promoter consensus
        if (count_promoters_ + (length() - pos)/PROM_SIZE >= dyntab_size_)
            reallocate_promoters();

        for (int32_t i = length() - 1; i >= pos; i--) {

            int8_t dist = is_promoter_lagging(i);
            if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
                bool notFound = true;
                for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                    if (promoters_[prom_idx] != nullptr) {
                        if (!promoters_[prom_idx]->leading_or_lagging)
                            if (promoters_[prom_idx]->pos == i) {
                                notFound = false;
                            }
                    }
                }

                if (notFound) {
                    int prom_idx;
                    {
                        prom_idx = count_promoters_;
                        count_promoters_ = count_promoters_ + 1;
                    }

                    promoters_[prom_idx] = new PromoterStruct(i, dist, false);;
                }
            }
        }
    }

    void DynTab_Metadata::look_for_new_lagging_promoters_starting_before(int32_t pos) {
        // Hamming distance of the sequence from the promoter consensus
        if (count_promoters_ + (pos)/PROM_SIZE >= dyntab_size_)
            reallocate_promoters();

        for (int32_t i = pos - 1; i >= 0; i--) {

            int8_t dist = is_promoter_lagging(i);
            if (dist <= 4) { // dist takes the hamming distance of the sequence from the consensus
                bool notFound = true;
                for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
                    if (promoters_[prom_idx] != nullptr) {
                        if (!promoters_[prom_idx]->leading_or_lagging)
                            if (promoters_[prom_idx]->pos == i) {
                                notFound = false;
                            }
                    }
                }

                if (notFound) {
                    int prom_idx;
                    {
                        prom_idx = count_promoters_;
                        count_promoters_ = count_promoters_ + 1;
                    }

                    promoters_[prom_idx] = new PromoterStruct(i, dist, false);;
                }
            }
        }
    }

    void DynTab_Metadata::promoters_included_in(int32_t pos_1,
                               int32_t pos_2,
                               std::vector<std::list<PromoterStruct*>>& promoters_list) {
        if (pos_1 < pos_2) {
            int32_t seg_length = pos_2 - pos_1;

            if (seg_length >= PROM_SIZE) {
                lst_promoters(LEADING, BETWEEN, pos_1, pos_2 - PROM_SIZE + 1,
                              promoters_list[LEADING]);
                lst_promoters(LAGGING, BETWEEN, pos_2, pos_1 + PROM_SIZE - 1,
                              promoters_list[LAGGING]);
            }
        }
        else {
            int32_t seg_length = length() + pos_2 - pos_1;
            if (seg_length >= PROM_SIZE) {
                bool is_near_end_of_genome = (pos_1 + PROM_SIZE > length());
                bool is_near_beginning_of_genome = (pos_2 - PROM_SIZE < 0);

                if (!is_near_end_of_genome && !is_near_beginning_of_genome) {
                    lst_promoters(LEADING, AFTER, pos_1, -1, promoters_list[LEADING]);
                    lst_promoters(LEADING, BEFORE, -1, pos_2 - PROM_SIZE + 1,promoters_list[LEADING]);
                    lst_promoters(LAGGING, AFTER, pos_2, -1, promoters_list[LAGGING]);
                    lst_promoters(LAGGING, BEFORE, -1, pos_1 + PROM_SIZE - 1,
                                  promoters_list[LAGGING]);
                }
                else if (!is_near_end_of_genome) // => && is_near_beginning_of_genome
                {
                    lst_promoters(LEADING, BETWEEN, pos_1, pos_2 - PROM_SIZE + 1 +
                                                           length(),
                                  promoters_list[LEADING]);
                    lst_promoters(LAGGING, AFTER, pos_2, -1, promoters_list[LAGGING]);
                    lst_promoters(LAGGING, BEFORE, -1, pos_1 + PROM_SIZE - 1,
                                  promoters_list[LAGGING]);
                }
                else if (!is_near_beginning_of_genome) // => && is_near_end_of_genome
                {
                    lst_promoters(LEADING, AFTER, pos_1, -1, promoters_list[LEADING]);
                    lst_promoters(LEADING, BEFORE, -1, pos_2 - PROM_SIZE + 1,
                                  promoters_list[LEADING]);
                    lst_promoters(LAGGING, BETWEEN, pos_2, pos_1 + PROM_SIZE - 1 -
                                                           length(),
                                  promoters_list[LAGGING]);
                }
                else // is_near_end_of_genome && is_near_beginning_of_genome
                {
                    lst_promoters(LEADING, BETWEEN, pos_1, pos_2 - PROM_SIZE + 1 +
                                                           length(),
                                  promoters_list[LEADING]);
                    lst_promoters(LAGGING, BETWEEN, pos_2, pos_1 + PROM_SIZE - 1 -
                                                           length(),
                                  promoters_list[LAGGING]);
                }
            }
        }
    }

    void DynTab_Metadata::extract_leading_promoters_starting_between(int32_t pos_1,
                                                    int32_t pos_2, std::list<PromoterStruct*>& extracted_promoters) {
        for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
            if (promoters_[prom_idx] != nullptr) {
                if (promoters_[prom_idx]->leading_or_lagging)
                    if ((promoters_[prom_idx]->pos >= pos_1) && (promoters_[prom_idx]->pos < pos_2)) {
                        extracted_promoters.push_back(promoters_[prom_idx]);
                        promoters_[prom_idx] = nullptr;
                    }
            }
        }
    }


    void DynTab_Metadata::extract_lagging_promoters_starting_between(int32_t pos_1,
                                                    int32_t pos_2,
                                                    std::list<PromoterStruct*>& extracted_promoters) {

        for (int prom_idx = 0; prom_idx < count_promoters_; prom_idx++) {
            if (promoters_[prom_idx] != nullptr) {
                if (!promoters_[prom_idx]->leading_or_lagging)
                    if ((promoters_[prom_idx]->pos >= pos_1) && (promoters_[prom_idx]->pos < pos_2)) {
                        extracted_promoters.push_back(promoters_[prom_idx]);
                        promoters_[prom_idx] = nullptr;
                    }
            }
        }
    }

    PromoterStruct* DynTab_Metadata::promoters(int idx) {
        return promoters_[idx];
    }

    void DynTab_Metadata::promoter_add(int idx, PromoterStruct*prom) {
        promoters_[idx] = prom;
    }

    PromoterStruct* DynTab_Metadata::promoter_next() {

      PromoterStruct* prom = promoters_[it_promoter_];
        it_promoter_++;

        return prom;
    }

    void DynTab_Metadata::promoter_begin() {
        it_promoter_ = 0;
        it_promoter_count_ = 0;
    }

    bool DynTab_Metadata::promoter_end() {
        return it_promoter_count_ == promoter_count();
    }


    int DynTab_Metadata::promoter_count() {
        return count_promoters_;
    }

    void DynTab_Metadata::set_promoters_count(int pcount) {
        count_promoters_ = pcount;
    }

    int DynTab_Metadata::terminator_count(int LoL) {
        if (LoL == LEADING)
            return (int) terminator_lead_.size();
        else
            return (int) terminator_lag_.size();
    }

    void DynTab_Metadata::terminator_add(int LoL, int dna_pos) {
        if (LoL == LEADING)
            terminator_lead_.insert(dna_pos);
        else
            terminator_lag_.insert(dna_pos);
    }

    int DynTab_Metadata::next_terminator(int LoL, int dna_pos) {
        if (LoL == LEADING) {
            auto it_rna_end = terminator_lead_.lower_bound(dna_pos);

            if (it_rna_end == terminator_lead_.end()) {
                it_rna_end = terminator_lead_.begin();
            }

            return *it_rna_end;
        } else {
            auto it_rna_end = terminator_lag_.upper_bound(dna_pos);


            if (it_rna_end == terminator_lag_.begin()) {
                it_rna_end = terminator_lag_.end();
                it_rna_end--;
            } else if ((*it_rna_end) != dna_pos)
                it_rna_end--;

            return *it_rna_end;
        }

    }

    void DynTab_Metadata::terminators_clear() {
        terminator_lead_.clear();
        terminator_lag_.clear();
    }

    Rna_7* DynTab_Metadata::rnas(int idx) {
        return rnas_[idx];
    }

    void DynTab_Metadata::rna_add(int idx, Rna_7*rna) {
        rnas_[idx] = rna;
    }

    void DynTab_Metadata::rna_add(int idx, int32_t t_begin, int32_t t_end,
                 int8_t t_leading_lagging, double t_e,
                 int32_t t_length) {
        rnas_[idx] = new Rna_7(t_begin, t_end,
                t_leading_lagging, t_e,
                t_length);
    }

    Rna_7* DynTab_Metadata::rna_next() {
      Rna_7* rna = *it_rna_;
        it_rna_++;
        return rna;
    }

    void DynTab_Metadata::rna_begin() {
        it_rna_ = rnas_.begin();
    }

    bool DynTab_Metadata::rna_end() {
        return it_rna_ == rnas_.end();
    }

    int DynTab_Metadata::rna_count() {
        return rna_count_;
    }

    void DynTab_Metadata::set_rna_count(int rcount) {
        rna_count_ = rcount;
    }

    void DynTab_Metadata::rnas_resize(int resize) {
        rnas_.resize(resize);
    }

    void DynTab_Metadata::rnas_clear() {
        rnas_.clear();
    }

    Protein_7* DynTab_Metadata::proteins(int idx) {
        return proteins_[idx];
    }

    void DynTab_Metadata::protein_add(int idx, Protein_7*prot) {
        //printf("Add proteins at %d\n",idx);
        proteins_[idx] = prot;
        protein_count_++;
    }

    Protein_7* DynTab_Metadata::protein_next() {
      Protein_7* prot = *it_protein_;
        it_protein_++;
        return prot;
    }

    void DynTab_Metadata::protein_begin() {
        it_protein_ = proteins_.begin();
    }

    bool DynTab_Metadata::protein_end() {
        return it_protein_ == proteins_.end();
    }

    int DynTab_Metadata::proteins_count() {
        return protein_count_;
    }

    void DynTab_Metadata::set_proteins_count(int pcount) {
        //protein_count_ = pcount;
    }

    void DynTab_Metadata::proteins_resize(int resize) {
        proteins_.resize(resize);
    }

    void DynTab_Metadata::proteins_clear() {
        proteins_.clear();
    }
}
