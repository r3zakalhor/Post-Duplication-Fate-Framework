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
//                              Includes
// =================================================================
#include "Alignment.h"
#include "Utils.h"
#include "Individual.h"

namespace aevol {


// ############################################################################
//
//                                Class Alignment
//
// ############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================

// =================================================================
//                             Destructors
// =================================================================

// =================================================================
//                            Public Methods
// =================================================================
VisAVis* Alignment::search_alignment_direct(const Dna* chrom_1,
                                            const int32_t seed_1,
                                            const Dna* chrom_2,
                                            const int32_t seed_2,
                                            const int16_t needed_score) {
  VisAVis* best_alignment = NULL;

  int16_t nb_diags = 2 * chrom_1->indiv()->align_max_shift() + 1;
  VisAVis* cur_vav = NULL;

  // TODO <david.parsons@inria.fr>: As VisAVis now contains its score, we
  // should adapt the code to make it more integrated
  int16_t cur_score;

  // Zone 1 (Index on the chromosome)
  // First base in working zone 1
  int32_t w_zone_1_first =
      seed_1 - chrom_1->indiv()->align_w_zone_h_len();
  // Last  base in working zone 1
  int32_t w_zone_1_last =
      seed_1 + chrom_1->indiv()->align_w_zone_h_len();
  // First base in extended zone 1
  int32_t x_zone_1_first =
      w_zone_1_first - chrom_1->indiv()->align_max_shift();

  // Zone 2 (Index on the chromosome)
  // First base in working zone 2
  int32_t w_zone_2_first =
      seed_2 - chrom_2->indiv()->align_w_zone_h_len();
  // Last  base in working zone 2
  int32_t w_zone_2_last =
      seed_2 + chrom_2->indiv()->align_w_zone_h_len();

  // This doesn't represent any point of interest
  // in the sequence but will spare some calculation
  int32_t w_zone_2_shifted_first =
      w_zone_2_first + chrom_2->indiv()->align_max_shift();


  // Parse diagonals
  for (int16_t cur_diag = 0 ; cur_diag < nb_diags ; cur_diag++) {
    // Initialize cur_vav according to the diagonal we are on
    if (cur_diag < chrom_1->indiv()->align_max_shift()) {
      cur_vav = new VisAVis(chrom_1, chrom_2,
                            x_zone_1_first + cur_diag, w_zone_2_first, DIRECT);
    }
    else if (cur_diag > chrom_1->indiv()->align_max_shift()) {
      cur_vav = new VisAVis(chrom_1, chrom_2,
                            w_zone_1_first, w_zone_2_shifted_first - cur_diag,
                            DIRECT);
    }
    else { // Central diagonal
      cur_vav = new VisAVis(chrom_1, chrom_2, w_zone_1_first, w_zone_2_first,
                            DIRECT);
    }

    // A sequence against itself is not an alignment
    if (chrom_1 == chrom_2 &&
        Utils::mod(cur_vav->i_1_, chrom_1->length()) ==
            Utils::mod(cur_vav->i_2_, chrom_2->length())) {
      delete cur_vav;
      continue;
    }

    cur_score = 0;

    // Parse current diagonal
    while (cur_vav->i_1_ <= w_zone_1_last || cur_vav->i_2_ <= w_zone_2_last) {
      // Re-initialize score and potential alignment starting point if
      // score <= 0
      if (cur_score <= 0) {
        cur_score = 0;
        if (best_alignment != NULL) {
          best_alignment->copy(cur_vav);
        }
        else {
          best_alignment = new VisAVis(*cur_vav);
        }
      }

      // Update Score
      if (cur_vav->match()) {
        cur_score += chrom_1->indiv()->align_match_bonus();

        // Check whether score is high enough to rearrange
        if (cur_score >= needed_score) {
          delete cur_vav;
          best_alignment->check_indices();

          best_alignment->score_ = cur_score;
          return best_alignment;
        }
      }
      else {
        cur_score -= chrom_1->indiv()->align_mismatch_cost();
      }

      // Step forward
      cur_vav->step_fwd();
    }

    delete cur_vav;
  }


  if (best_alignment != NULL) {
    delete best_alignment;
  }

  return NULL; // Didn't find any alignment with sufficient score.
}


VisAVis* Alignment::search_alignment_indirect(const Dna* chrom_1,
                                              const int32_t seed_1,
                                              const Dna* chrom_2,
                                              const int32_t seed_2,
                                              const int16_t needed_score) {
  VisAVis * best_alignment = NULL;

  int16_t nb_diags = 2 * chrom_1->indiv()->align_max_shift() + 1;
  int16_t cur_score;
  VisAVis * cur_vav = NULL;

  // Zone 1 (Indice on the chromosome)
  // First base in working zone 1
  int32_t w_zone_1_first =
      seed_1 - chrom_1->indiv()->align_w_zone_h_len();
  // Last  base in working zone 1
  int32_t w_zone_1_last =
      seed_1 + chrom_1->indiv()->align_w_zone_h_len();
  // First base in extended zone 1
  int32_t x_zone_1_first =
      w_zone_1_first - chrom_1->indiv()->align_max_shift();

  // Zone 2 (Indice on the chromosome)
  // ********** WARNING **********
  // Because we are reading backwards, every read is shifted by -1.
  //
  //    0 1 2 3 4 5 6 7 8 9
  //    |_|_|_|_|_|_|_|_|_|_|   On this sequence, if you consider index 5,
  //                            the nucleotide you need to consider
  //                            for a backwards read is seq[4], NOT seq[5].
  //
  // This is managed by class VisAVis.
  //
  //    a b c d e f g h i j       Creating a new VisAVis(F, 5, LAGGING);
  //    |_|_|_|_|_|_|_|_|_|_|     and then asking: my_vis_a_vis->match() will
  //    | | | | | | | | | | |     compare seq1[F] and seq2[4]
  //      9 8 7 6 5 4 3 2 1 0     (i.e. its "real" vis_a_vis).
  //                              my_vis_a_vis->match() returns true
  //                              if (seq1[F] == complementary(seq2[4]))
  //
  // In this example and with a align_w_zone_h_len of 2 and a max_shift of 1,
  // the working_zone and extended_zone will be defined by
  // w_zone_1_first = d   // (f-2)
  // w_zone_1_last  = h   // (f+2)
  // x_zone_1_first = c   // (f-2-1)
  // w_zone_2_first = 7   // The "real" vis_a_vis for d (5+2)
  // w_zone_2_last  = 3   // The "real" vis_a_vis for h (5-2)
  // w_zone_2_shifted_first = 4 // (5-1).
  //
  // It is worth noting that the working zone is defined as follows:
  //
  //    d e f g h       Hence, while 'd' and 'h' do correspond to the zone's
  //    |_|_|_|_|_|     first and last nucleotides' indices,
  //    | | | | | |     '7' and '3' are indeed shifted and it seems they don't
  //      6 5 4 3 2     correspond to anything. I should have coded the limits
  //                    (i.e. 'd' and 'i' + '7' and '2') in the first place...
  //
  //

  // First base in working zone 2
  int32_t w_zone_2_first =
      seed_2 + chrom_2->indiv()->align_w_zone_h_len();
  // Last  base in working zone 2
  int32_t w_zone_2_last =
      seed_2 - chrom_2->indiv()->align_w_zone_h_len();

  // This doesn't represent any interesting point
  // in the sequence but will spare some calculation
  int32_t w_zone_2_shifted_first =
      w_zone_2_first - chrom_2->indiv()->align_max_shift();


  // Parse diagonals
  for (int16_t cur_diag = 0 ; cur_diag < nb_diags ; cur_diag++) {
    // Initialize cur_vav according to the diagonal we are on
    if (cur_diag < chrom_1->indiv()->align_max_shift()) {
      cur_vav = new VisAVis(chrom_1, chrom_2,
                            x_zone_1_first + cur_diag, w_zone_2_first,
                            INDIRECT);
    }
    else if (cur_diag > chrom_1->indiv()->align_max_shift()) {
      cur_vav = new VisAVis(chrom_1, chrom_2,
                            w_zone_1_first, w_zone_2_shifted_first + cur_diag,
                            INDIRECT);
    }
    else { // Central diagonal
      cur_vav = new VisAVis(chrom_1, chrom_2, w_zone_1_first, w_zone_2_first,
                            INDIRECT);
    }

    cur_score = 0;


    // Parse diagonal
    while (cur_vav->i_1_ <= w_zone_1_last || cur_vav->i_2_ >= w_zone_2_last) {
      // Re-initialize score and potential alignment starting point if
      // score <= 0
      if (cur_score <= 0) {
        cur_score = 0;
        if (best_alignment != NULL) {
          best_alignment->copy(cur_vav);
        }
        else {
          best_alignment = new VisAVis(*cur_vav);
        }
      }

      // Update Score
      if (cur_vav->match()) {
        cur_score += chrom_1->indiv()->align_match_bonus();

        // Check whether score is high enough to rearrange
        if (cur_score >= needed_score) {
          delete cur_vav;
          best_alignment->check_indices();

          best_alignment->score_ = cur_score;
          return best_alignment;
        }
      }
      else {
        cur_score -= chrom_1->indiv()->align_mismatch_cost();
      }

      // Step forward
      cur_vav->step_fwd();
    }

    delete cur_vav;
  }


  if (best_alignment != NULL) {
    delete best_alignment;
  }

  return NULL; // Didn't find any alignment with sufficient score.
}

} // namespace aevol
// =================================================================
//                           Protected Methods
// =================================================================
