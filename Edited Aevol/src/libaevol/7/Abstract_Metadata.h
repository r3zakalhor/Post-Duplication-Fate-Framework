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

#ifndef AEVOL_ABSTRACT_METADATA_H
#define AEVOL_ABSTRACT_METADATA_H

#include "Dna_7.h"
#include "Individual_7.h"
#include "ExpManager_7.h"
#include "ae_enums.h"
#include "Promoter.h"
#include "ExpManager.h"

#include <cstdint>
#include <list>
#include <vector>

namespace aevol {
    class Rna_7;
    class Protein_7;

    class Abstract_Metadata {
    public:
     Abstract_Metadata(Individual_7* indiv,
                       Abstract_Metadata* metadata) {
            indiv_ = indiv;
        }

        Abstract_Metadata(Individual_7* indiv,
                                int32_t* lead_prom_pos = nullptr, int8_t* lead_prom_error = nullptr, int32_t lead_prom_size = -1,
                                int32_t* lag_prom_pos = nullptr, int8_t* lag_prom_error = nullptr, int32_t lag_prom_size = -1) {
            indiv_ = indiv;
        }

        Abstract_Metadata(Individual_7* indiv) {
            indiv_ = indiv;
        }


        virtual ~Abstract_Metadata() {};
        /** Getter **/
        /*** Promoters ***/
        inline int8_t is_promoter_leading(int pos) {
            #ifdef BASE_2
                #ifdef VANILLA_SEARCH
                int8_t prom_dist_leading[26];
                int len = indiv_->dna_->length();

                for (int motif_id = 0; motif_id < 22; motif_id++) {
                        prom_dist_leading[motif_id] =
                                PROM_SEQ_LEAD[motif_id] ==
                                indiv_->dna_->data_[pos + motif_id >= len ? pos + motif_id - len : pos + motif_id]
                                ? 0 : 1;
                }

                return prom_dist_leading[0] +
                       prom_dist_leading[1] +
                       prom_dist_leading[2] +
                       prom_dist_leading[3] +
                       prom_dist_leading[4] +
                       prom_dist_leading[5] +
                       prom_dist_leading[6] +
                       prom_dist_leading[7] +
                       prom_dist_leading[8] +
                       prom_dist_leading[9] +
                       prom_dist_leading[10] +
                       prom_dist_leading[11] +
                       prom_dist_leading[12] +
                       prom_dist_leading[13] +
                       prom_dist_leading[14] +
                       prom_dist_leading[15] +
                       prom_dist_leading[16] +
                       prom_dist_leading[17] +
                       prom_dist_leading[18] +
                       prom_dist_leading[19] +
                       prom_dist_leading[20] +
                       prom_dist_leading[21];
#elif CMOD_SEARCH
            int8_t prom_dist_leading[26];
                int len = indiv_->dna_->length();

                for (int motif_id = 0; motif_id < 22; motif_id++) {
                        prom_dist_leading[motif_id] =
                                PROM_SEQ_LEAD[motif_id] ==
                                indiv_->dna_->data_[Utils::mod(pos + motif_id, len)]
                                ? 0 : 1;
                }

                return prom_dist_leading[0] +
                       prom_dist_leading[1] +
                       prom_dist_leading[2] +
                       prom_dist_leading[3] +
                       prom_dist_leading[4] +
                       prom_dist_leading[5] +
                       prom_dist_leading[6] +
                       prom_dist_leading[7] +
                       prom_dist_leading[8] +
                       prom_dist_leading[9] +
                       prom_dist_leading[10] +
                       prom_dist_leading[11] +
                       prom_dist_leading[12] +
                       prom_dist_leading[13] +
                       prom_dist_leading[14] +
                       prom_dist_leading[15] +
                       prom_dist_leading[16] +
                       prom_dist_leading[17] +
                       prom_dist_leading[18] +
                       prom_dist_leading[19] +
                       prom_dist_leading[20] +
                       prom_dist_leading[21];

#elif AMOD_SEARCH
            int8_t prom_dist_leading[26];
                int len = indiv_->dna_->length();

                if (pos+22 >= len) {
                    for (int motif_id = 0; motif_id < 22; motif_id++) {
                            prom_dist_leading[motif_id] =
                                    PROM_SEQ_LEAD[motif_id] ==
                                    indiv_->dna_->data_[Utils::mod(pos + motif_id, len)]
                                    ? 0 : 1;
                    }
                } else {
                    for (int motif_id = 0; motif_id < 22; motif_id++) {
                            prom_dist_leading[motif_id] =
                                    PROM_SEQ_LEAD[motif_id] ==
                                    indiv_->dna_->data_[pos + motif_id]
                                    ? 0 : 1;
                    }
                }

                return prom_dist_leading[0] +
                       prom_dist_leading[1] +
                       prom_dist_leading[2] +
                       prom_dist_leading[3] +
                       prom_dist_leading[4] +
                       prom_dist_leading[5] +
                       prom_dist_leading[6] +
                       prom_dist_leading[7] +
                       prom_dist_leading[8] +
                       prom_dist_leading[9] +
                       prom_dist_leading[10] +
                       prom_dist_leading[11] +
                       prom_dist_leading[12] +
                       prom_dist_leading[13] +
                       prom_dist_leading[14] +
                       prom_dist_leading[15] +
                       prom_dist_leading[16] +
                       prom_dist_leading[17] +
                       prom_dist_leading[18] +
                       prom_dist_leading[19] +
                       prom_dist_leading[20] +
                       prom_dist_leading[21];

#elif B1MOD_SEARCH
            int8_t prom_dist_leading[26];
                int32_t len = indiv_->dna_->length();

                for (int motif_id = 0; motif_id < 22; motif_id++) {
                        prom_dist_leading[motif_id] =
                                PROM_SEQ_LEAD[motif_id] ==
                                indiv_->dna_->data_[(pos + motif_id) - ((pos + motif_id) >= len) * len]
                                ? 0 : 1;
                }

                return prom_dist_leading[0] +
                       prom_dist_leading[1] +
                       prom_dist_leading[2] +
                       prom_dist_leading[3] +
                       prom_dist_leading[4] +
                       prom_dist_leading[5] +
                       prom_dist_leading[6] +
                       prom_dist_leading[7] +
                       prom_dist_leading[8] +
                       prom_dist_leading[9] +
                       prom_dist_leading[10] +
                       prom_dist_leading[11] +
                       prom_dist_leading[12] +
                       prom_dist_leading[13] +
                       prom_dist_leading[14] +
                       prom_dist_leading[15] +
                       prom_dist_leading[16] +
                       prom_dist_leading[17] +
                       prom_dist_leading[18] +
                       prom_dist_leading[19] +
                       prom_dist_leading[20] +
                       prom_dist_leading[21];
#elif B2MOD_SEARCH
            int8_t prom_dist_leading[26];
                int32_t len = indiv_->dna_->length();

                for (int motif_id = 0; motif_id < 22; motif_id++) {
                        prom_dist_leading[motif_id] =
                                PROM_SEQ_LEAD[motif_id] ==
                                indiv_->dna_->data_[(pos + motif_id) + ((((uint32_t)((pos + motif_id) - len)) >> 31) -1 )* len]
                                ? 0 : 1;
                }

                return prom_dist_leading[0] +
                       prom_dist_leading[1] +
                       prom_dist_leading[2] +
                       prom_dist_leading[3] +
                       prom_dist_leading[4] +
                       prom_dist_leading[5] +
                       prom_dist_leading[6] +
                       prom_dist_leading[7] +
                       prom_dist_leading[8] +
                       prom_dist_leading[9] +
                       prom_dist_leading[10] +
                       prom_dist_leading[11] +
                       prom_dist_leading[12] +
                       prom_dist_leading[13] +
                       prom_dist_leading[14] +
                       prom_dist_leading[15] +
                       prom_dist_leading[16] +
                       prom_dist_leading[17] +
                       prom_dist_leading[18] +
                       prom_dist_leading[19] +
                       prom_dist_leading[20] +
                       prom_dist_leading[21];
#else // OLD_SEARCH
            int8_t prom_dist_leading = 0;
            int len = indiv_->dna_->length();

            for (int motif_id = 0; motif_id < 22; motif_id++) {
                prom_dist_leading +=
                        PROM_SEQ_LEAD[motif_id] ==
                        indiv_->dna_->data_[pos + motif_id >= len ? pos + motif_id - len : pos + motif_id]
                        ? 0 : 1;
                if (prom_dist_leading>PROM_MAX_DIFF)
                    break;
            }

            return prom_dist_leading;
#endif
#elif BASE_4
            int8_t prom_dist_leading = 0;
            int len = indiv_->dna_->length();

            for (int motif_id = 0; motif_id < PROM_SIZE; motif_id++) {
                prom_dist_leading +=
                    PROM_SEQ[motif_id] != indiv_->dna_->data_[
                        pos + motif_id >= len ?
                        pos - len + motif_id : pos + motif_id
                    ];

                if (prom_dist_leading > PROM_MAX_DIFF)
                    break;
            }

            return prom_dist_leading;
#endif
        }

        inline int8_t is_promoter_lagging(int pos) {
#ifdef BASE_2
#ifdef VANILLA_SEARCH
            int8_t prom_dist[26];
                int len = indiv_->dna_->length();

                for (int motif_id = 0; motif_id < 22; motif_id++) {
                        prom_dist[motif_id] =
                                PROM_SEQ_LAG[motif_id] ==
                                indiv_->dna_->data_[pos - motif_id < 0 ? len + pos - motif_id : pos - motif_id]
                                ? 0 : 1;
                }

                return prom_dist[0] +
                       prom_dist[1] +
                       prom_dist[2] +
                       prom_dist[3] +
                       prom_dist[4] +
                       prom_dist[5] +
                       prom_dist[6] +
                       prom_dist[7] +
                       prom_dist[8] +
                       prom_dist[9] +
                       prom_dist[10] +
                       prom_dist[11] +
                       prom_dist[12] +
                       prom_dist[13] +
                       prom_dist[14] +
                       prom_dist[15] +
                       prom_dist[16] +
                       prom_dist[17] +
                       prom_dist[18] +
                       prom_dist[19] +
                       prom_dist[20] +
                       prom_dist[21];
#elif CMOD_SEARCH
            int8_t prom_dist[26];
            int len = indiv_->dna_->length();

            for (int motif_id = 0; motif_id < 22; motif_id++) {
                prom_dist[motif_id] =
                        PROM_SEQ_LAG[motif_id] ==
                        indiv_->dna_->data_[Utils::mod(pos - motif_id,len)]
                        ? 0 : 1;
            }

            return prom_dist[0] +
                   prom_dist[1] +
                   prom_dist[2] +
                   prom_dist[3] +
                   prom_dist[4] +
                   prom_dist[5] +
                   prom_dist[6] +
                   prom_dist[7] +
                   prom_dist[8] +
                   prom_dist[9] +
                   prom_dist[10] +
                   prom_dist[11] +
                   prom_dist[12] +
                   prom_dist[13] +
                   prom_dist[14] +
                   prom_dist[15] +
                   prom_dist[16] +
                   prom_dist[17] +
                   prom_dist[18] +
                   prom_dist[19] +
                   prom_dist[20] +
                   prom_dist[21];
#elif AMOD_SEARCH
            int8_t prom_dist[26];
            int len = indiv_->dna_->length();

            if (pos-22 < 0) {

                for (int motif_id = 0; motif_id < 22; motif_id++) {
                    prom_dist[motif_id] =
                            PROM_SEQ_LAG[motif_id] ==
                            indiv_->dna_->data_[Utils::mod(pos - motif_id,len)]
                            ? 0 : 1;
                }
            } else {
                for (int motif_id = 0; motif_id < 22; motif_id++) {
                    prom_dist[motif_id] =
                            PROM_SEQ_LAG[motif_id] ==
                            indiv_->dna_->data_[pos - motif_id]
                            ? 0 : 1;
                }
            }

            return prom_dist[0] +
                   prom_dist[1] +
                   prom_dist[2] +
                   prom_dist[3] +
                   prom_dist[4] +
                   prom_dist[5] +
                   prom_dist[6] +
                   prom_dist[7] +
                   prom_dist[8] +
                   prom_dist[9] +
                   prom_dist[10] +
                   prom_dist[11] +
                   prom_dist[12] +
                   prom_dist[13] +
                   prom_dist[14] +
                   prom_dist[15] +
                   prom_dist[16] +
                   prom_dist[17] +
                   prom_dist[18] +
                   prom_dist[19] +
                   prom_dist[20] +
                   prom_dist[21];
#elif B1MOD
            int8_t prom_dist[26];
                int32_t len = indiv_->dna_->length();

                for (int motif_id = 0; motif_id < 22; motif_id++) {
                        prom_dist[motif_id] =
                                PROM_SEQ_LAG[motif_id] ==
                                indiv_->dna_->data_[(pos - motif_id) + ((pos - motif_id) < 0) * len]
                                ? 0 : 1;
                }

                return prom_dist[0] +
                       prom_dist[1] +
                       prom_dist[2] +
                       prom_dist[3] +
                       prom_dist[4] +
                       prom_dist[5] +
                       prom_dist[6] +
                       prom_dist[7] +
                       prom_dist[8] +
                       prom_dist[9] +
                       prom_dist[10] +
                       prom_dist[11] +
                       prom_dist[12] +
                       prom_dist[13] +
                       prom_dist[14] +
                       prom_dist[15] +
                       prom_dist[16] +
                       prom_dist[17] +
                       prom_dist[18] +
                       prom_dist[19] +
                       prom_dist[20] +
                       prom_dist[21];

#elif B2MOD
            int8_t prom_dist[26];
                int32_t len = indiv_->dna_->length();

                for (int motif_id = 0; motif_id < 22; motif_id++) {
                        prom_dist[motif_id] =
                                PROM_SEQ_LAG[motif_id] ==
                                indiv_->dna_->data_[ (pos - motif_id) + (((uint32_t)((pos - motif_id))) >> 31) * len]
                                ? 0 : 1;
                }

                return prom_dist[0] +
                       prom_dist[1] +
                       prom_dist[2] +
                       prom_dist[3] +
                       prom_dist[4] +
                       prom_dist[5] +
                       prom_dist[6] +
                       prom_dist[7] +
                       prom_dist[8] +
                       prom_dist[9] +
                       prom_dist[10] +
                       prom_dist[11] +
                       prom_dist[12] +
                       prom_dist[13] +
                       prom_dist[14] +
                       prom_dist[15] +
                       prom_dist[16] +
                       prom_dist[17] +
                       prom_dist[18] +
                       prom_dist[19] +
                       prom_dist[20] +
                       prom_dist[21];
#else
            int8_t prom_dist = 0;
            int len = indiv_->dna_->length();

            for (int motif_id = 0; motif_id < 22; motif_id++) {
                prom_dist+=
                        PROM_SEQ_LAG[motif_id] ==
                        indiv_->dna_->data_[Utils::mod(pos - motif_id,len)]
                        ? 0 : 1;
                if (prom_dist>PROM_MAX_DIFF)
                    break;
            }

            return prom_dist;
#endif
#elif BASE_4
            int8_t prom_dist = 0;
            int len = indiv_->dna_->length();

            for (int motif_id = 0; motif_id < PROM_SIZE; motif_id++) {
                prom_dist +=
                    PROM_SEQ[motif_id]
                    != get_complementary_base(indiv_->dna_->data_[  // this is the same as comparing the two opposites (comparing opposite PROM with opposite strand)
                        pos - motif_id < 0 ?
                        pos + len - motif_id : pos - motif_id
                    ]);
                if (prom_dist > PROM_MAX_DIFF)
                    break;
            }

            return prom_dist;
#endif
        }

        virtual void lst_promoters(bool lorl,
                           Position before_after_btw, // with regard to the strand's reading direction
                           int32_t pos1,
                           int32_t pos2, std::list<PromoterStruct*>& motif_list) = 0;

        /*** Shine Dal + Start Codon ***/
        inline int8_t is_shine_dal_start_prot_leading(int pos) {
            #ifdef BASE_2
            int32_t len = indiv_->dna_->length();
            int8_t start[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

            for (int k = 0; k < 9; k++) {
                int k_t = k >= 6 ? k + 4 : k;

                if (indiv_->dna_->data_[(pos + k_t) + ((((uint32_t)((pos + k_t) - len)) >> 31) -1 )* len] ==
                    SHINE_DAL_SEQ_LEAD_7[k]) {
                    start[k] = true;
                } else {
                    start[k] = false;
                    break;
                }
            }

            return start[0] + start[1] + start[2] + start[3] + start[4] + start[5] + start[6]
                   + start[7] + start[8];
            #elif BASE_4
            int32_t len = indiv_->dna_->length();

            
            for (int k = 0; k < SHINE_DAL_SIZE; k++) {
                if (indiv_->dna_->get_lead(pos + k) != SHINE_DAL_SEQ[k]) {
                    // printf("FAILED Search SHine Dal at %d => %d\n",pos,k);
                    return false;
                }
            }

            int codon_pos = pos + SHINE_DAL_SIZE + SHINE_START_SPACER;

            // if (indiv_->indiv_id == 884) printf("%d -- %d -- Shine dal found %d, start %d : %c %c %c SHine %c %c %c %c %c %c\n",
            //         AeTime::time(),
            //         indiv_->indiv_id,pos,bases_to_codon_value(
            //         indiv_->dna_->get_lead(codon_pos    ),
            //         indiv_->dna_->get_lead(codon_pos + 1),
            //         indiv_->dna_->get_lead(codon_pos + 2)),
            //         indiv_->dna_->get_lead(codon_pos    ),
            //         indiv_->dna_->get_lead(codon_pos + 1),
            //         indiv_->dna_->get_lead(codon_pos + 2),
            //         indiv_->dna_->get_lead(pos    ),
            //         indiv_->dna_->get_lead(pos + 1),
            //         indiv_->dna_->get_lead(pos + 2),
            //         indiv_->dna_->get_lead(pos + 3),
            //         indiv_->dna_->get_lead(pos + 4),
            //         indiv_->dna_->get_lead(pos + 5));


            if(!is_start_codon(bases_to_codon_value(
                    indiv_->dna_->get_lead(codon_pos    ),
                    indiv_->dna_->get_lead(codon_pos + 1),
                    indiv_->dna_->get_lead(codon_pos + 2)
            ))) {
              return false;
            }

            return true;
            #endif
        }

        inline int8_t is_shine_dal_start_prot_lagging(int pos) {
            #ifdef BASE_2
            int32_t len = indiv_->dna_->length();
            int8_t start = false;

            for (int k = 0; k < 9; k++) {
                int k_t = k >= 6 ? k + 4 : k;

                if (indiv_->dna_->data_[(pos - k_t) + (((uint32_t)((pos - k_t))) >> 31) * len] ==
                    SHINE_DAL_SEQ_LAG_7[k]) {
                    start = true;
                } else {
                    start = false;
                    break;
                }
            }
            return start;
            #elif BASE_4
            int32_t len = indiv_->dna_->length();

            for (int k = 0; k < SHINE_DAL_SIZE; k++) {
              // call to lead instead of lag to simplify computation
              // => (complementary(a) != complementary(b) <=> a != b)
            //   if (pos == 2002) printf("%d :: %d --> %c != %c\n",pos-k,k,indiv_->dna_->get_lead(pos - k),SHINE_DAL_SEQ_LAG[k]);

              if (indiv_->dna_->get_lag(pos - k) != SHINE_DAL_SEQ[k]) {
                return false;
              }
            }

            int codon_pos = pos - SHINE_DAL_SIZE - SHINE_START_SPACER;

            // if (indiv_->indiv_id == 884) printf("%d -- %d -- Shine dal found %d, start %d : %c %c %c   %c %c %c %c %c %c\n",AeTime::time(),
            //         indiv_->indiv_id,pos,bases_to_codon_value(
            //         indiv_->dna_->get_lag(codon_pos    ),
            //         indiv_->dna_->get_lag(codon_pos - 1),
            //         indiv_->dna_->get_lag(codon_pos - 2)),
            //         indiv_->dna_->get_lag(codon_pos    ),
            //         indiv_->dna_->get_lag(codon_pos - 1),
            //         indiv_->dna_->get_lag(codon_pos - 2),
            //         indiv_->dna_->get_lag(pos    ),
            //         indiv_->dna_->get_lag(pos - 1),
            //         indiv_->dna_->get_lag(pos - 2),
            //         indiv_->dna_->get_lag(pos - 3),
            //         indiv_->dna_->get_lag(pos - 4),
            //         indiv_->dna_->get_lag(pos - 5));

            if(!is_start_codon(bases_to_codon_value(
                    indiv_->dna_->get_lag(codon_pos    ),
                    indiv_->dna_->get_lag(codon_pos - 1),
                    indiv_->dna_->get_lag(codon_pos - 2)
            ))) {
              return false;
            }

            return true;
            #endif
        }

        /*** Promoters ***/
        virtual PromoterStruct* promoters(int idx) = 0;
        virtual void promoter_add(int idx, PromoterStruct* prom) = 0;

        virtual int promoter_count() = 0;
        virtual void set_promoters_count(int rcount) = 0;

        virtual PromoterStruct* promoter_next() = 0;
        virtual void promoter_begin() = 0;
        virtual bool promoter_end() = 0;

        /*** Terminators ***/
        #ifdef BASE_4
        inline bool is_terminator_leading(int pos) {
            int32_t len = indiv_->dna_->length();

            // if (indiv_->indiv_id == 817 && pos == 653) {
            //     printf("SIMD/NOT FOUND a Terminator %d : \n",pos);
            //     for(auto i = 0; i < TERM_STEM_SIZE; i++) {
            //         printf("%c == %c\n",indiv_->dna_->data_[Utils::mod(pos + i,len)], 
            //                     indiv_->dna_->data_[Utils::mod(pos + (TERM_SIZE-1) - i,len)]);
            //     }
            //     printf("\n");
            // }

            // Palindrome
            // if (pos + TERM_SIZE -1 > len) {
            //     for(auto i = 0; i < TERM_STEM_SIZE; i++) {
            //         if (!is_complementary_base(indiv_->dna_->data_[pos + i], 
            //                     indiv_->dna_->data_[pos + (TERM_SIZE-1) - i])) {
            //         return false;
            //         }
            //     }
            // } else {
                for(auto i = 0; i < TERM_STEM_SIZE; i++) {
                    if (!is_complementary_base(indiv_->dna_->data_[Utils::mod(pos + i,len)], 
                                indiv_->dna_->data_[Utils::mod(pos + (TERM_SIZE-1) - i,len)])) {
                        return false;
                    }
                }
            // }


            // if (indiv_->indiv_id == 817) {
            //     printf("SIMD/Found a Terminator %d : \n",pos);
            //     for(auto i = 0; i < TERM_STEM_SIZE; i++) {
            //         printf("%c == %c\n",indiv_->dna_->data_[Utils::mod(pos + i,len)], 
            //                     indiv_->dna_->data_[Utils::mod(pos + (TERM_SIZE-1) - i,len)]);
            //     }
            //     printf("\n");
            // }

            // Poly-A sequence
            for (auto i = 0;
                 i < indiv_->exp_m_->exp_s()->terminator_polya_sequence_length();
                 i++) {
                if (indiv_->dna_->data_[Utils::mod(pos + TERM_SIZE + i, len)] != BASE_A) {
                  return false;
                }
            }

            return true;
        }

        inline bool is_terminator_lagging(int pos) {
            int32_t len = indiv_->dna_->length();

            // Palindrome
            if (pos - TERM_SIZE -1 >= 0) {
                for(auto i = 0; i < TERM_STEM_SIZE; i++)
                {
                    if (!is_complementary_base(indiv_->dna_->data_[pos - i], indiv_->dna_->data_[pos - (TERM_SIZE-1) + i])) {
                    return false;
                    }
                }
            } else {
                for(auto i = 0; i < TERM_STEM_SIZE; i++)
                {
                    if (!is_complementary_base(indiv_->dna_->data_[Utils::mod(pos - i,len)], indiv_->dna_->data_[Utils::mod(pos - (TERM_SIZE-1) + i,len)])) {
                    return false;
                    }
                }
            }


            // Poly-A sequence
            for (auto i = 0;
                 i < indiv_->exp_m_->exp_s()->terminator_polya_sequence_length();
                 i++) {
                if (indiv_->dna_->data_[Utils::mod(pos - TERM_SIZE - i, len)] != BASE_A) {
                  return false;
                }
            }

            return true;
        }
        #endif

        virtual int terminator_count(int LoL) = 0;
        virtual void terminator_add(int LoL, int dna_pos) = 0;

        virtual int next_terminator(int LoL, int dna_pos) = 0;

        virtual void terminators_clear() = 0;

        /*** RNAs ***/
        virtual Rna_7* rnas(int idx) = 0;
        virtual void rna_add(int idx, Rna_7* rna, PromoterStruct* prom = nullptr) = 0;
        virtual void rna_add(int idx, int32_t t_begin, int32_t t_end, int8_t t_leading_lagging, double t_e,
                             int32_t t_length, PromoterStruct* prom = nullptr) = 0;
        // virtual void rna_update(Rna_7* rna, int32_t t_begin, int32_t t_end, int8_t t_leading_lagging, double t_e,
        //                      int32_t t_length) = 0;

        virtual int rna_count() = 0;
        virtual void set_rna_count(int rcount) = 0;

        virtual void rnas_resize(int resize) = 0;
        virtual void rnas_clear() = 0;

        virtual Rna_7* rna_next() = 0;
        virtual void rna_begin() = 0;
        virtual bool rna_end() = 0;

        /*** Proteins ***/
        virtual Protein_7* proteins(int idx) = 0;
        virtual void protein_add(int idx, Protein_7* prot) = 0;

        virtual int proteins_count() = 0;
        virtual void set_proteins_count(int pcount) = 0;

        virtual void proteins_resize(int resize) = 0;
        virtual void proteins_clear() = 0;


        virtual Protein_7* protein_next() = 0;
        virtual void protein_begin() = 0;
        virtual bool protein_end() = 0;
        /** Search and update **/
        virtual void remove_promoters_around(int32_t pos_1) = 0;
        virtual void remove_promoters_around(int32_t pos_1, int32_t pos_2) = 0;
        virtual void remove_all_promoters() = 0;

        virtual void look_for_new_promoters_around(int32_t pos_1, int32_t pos_2) = 0;
        virtual void look_for_new_promoters_around(int32_t pos) = 0;

        virtual void locate_promoters() = 0;

        virtual void move_all_promoters_after(int32_t pos, int32_t delta_pos) = 0;

        virtual void duplicate_promoters_included_in(int32_t pos_1,
                                             int32_t pos_2,
                                             std::vector<std::list<PromoterStruct*>>& duplicated_promoters) = 0;
        virtual void extract_promoters_included_in(int32_t pos_1,
                                           int32_t pos_2, std::vector<std::list<PromoterStruct*>>& extracted_promoters) = 0;
        virtual void insert_promoters(std::vector<std::list<PromoterStruct*>>& promoters_to_insert) = 0;
        virtual void insert_promoters_at(std::vector<std::list<PromoterStruct*>>& promoters_to_insert,
                                 int32_t pos) = 0;

        virtual void invert_promoters_included_in(int32_t pos1,
                                          int32_t pos2) = 0;


        static void shift_promoters(
                std::vector<std::list<PromoterStruct*>>& promoters_to_shift,
                int32_t delta_pos,
                int32_t seq_length);
        static void invert_promoters(std::vector<std::list<PromoterStruct*>>& promoter_lists,
                                     int32_t pos1,
                                     int32_t pos2);

        virtual void remove_leading_promoters_starting_between(int32_t pos_1,
                                                       int32_t pos_2) = 0;
        virtual void remove_leading_promoters_starting_after(int32_t pos) = 0;
        virtual void remove_leading_promoters_starting_before(int32_t pos) = 0;

        virtual void remove_lagging_promoters_starting_between(int32_t pos_1,
                                                       int32_t pos_2) = 0;
        virtual void remove_lagging_promoters_starting_after(int32_t pos) = 0;
        virtual void remove_lagging_promoters_starting_before(int32_t pos) = 0;

        virtual void move_all_leading_promoters_after(int32_t pos, int32_t delta_pos) = 0;
        virtual void move_all_lagging_promoters_after(int32_t pos,int32_t delta_pos) = 0;

        virtual void look_for_new_leading_promoters_starting_between(int32_t pos_1, int32_t pos_2) = 0;
        virtual void look_for_new_leading_promoters_starting_after(int32_t pos) = 0;
        virtual void look_for_new_leading_promoters_starting_before(int32_t pos) = 0;

        virtual void look_for_new_lagging_promoters_starting_between(int32_t pos_1,int32_t pos_2) = 0;
        virtual void look_for_new_lagging_promoters_starting_after(int32_t pos) = 0;
        virtual void look_for_new_lagging_promoters_starting_before(int32_t pos) = 0;

        virtual void promoters_included_in(int32_t pos_1,
                                   int32_t pos_2,
                                   std::vector<std::list<PromoterStruct*>>& promoters_list) = 0;

        virtual void extract_leading_promoters_starting_between(int32_t pos_1,
                                                        int32_t pos_2, std::list<PromoterStruct*>& extracted_promoters) = 0;

        virtual void extract_lagging_promoters_starting_between(int32_t pos_1,
                                                        int32_t pos_2,
                                                        std::list<PromoterStruct*>& extracted_promoters) = 0;


        int32_t length() { return indiv_->dna_->length(); };


        int32_t rna_count_ = 0;
    // protected:
     Individual_7* indiv_;
    };
}

#endif
