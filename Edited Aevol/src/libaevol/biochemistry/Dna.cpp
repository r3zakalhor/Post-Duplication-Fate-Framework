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
#include "Dna.h"

#include <cinttypes>
#include <cstdio>
#include <cmath>

#include <list>
#include <vector>

#ifndef __OPENMP_GPU
#include <algorithm>
#endif

#include "ExpManager.h"
#include "ExpSetup.h"
#include "GeneticUnit.h"
#include "Individual.h"
#include "Rna.h"
#include "Utils.h"
#include "VisAVis.h"
#include "Alignment.h"
#include "Mutation.h"
#include "PointMutation.h"
#include "SmallInsertion.h"
#include "SmallDeletion.h"
#include "Duplication.h"
#include "Deletion.h"
#include "Translocation.h"
#include "Inversion.h"
#include "InsertionHT.h"
#include "ReplacementHT.h"
#include "MutationEvent.h"
#include "DnaMutator.h"

namespace aevol {

// ###########################################################################
//
//                                 Class Dna
//
// ###########################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
/**
 * Create a random dna sequence of length <length> belonging to <gen_unit>.
 */
Dna::Dna(GeneticUnit* gen_unit,
         int32_t length,
         std::shared_ptr<JumpingMT> prng) :
    ae_string(length, prng) {
  gen_unit_ = gen_unit;
  exp_m_ = gen_unit->exp_m();
  indiv_ = gen_unit->indiv();
}

/**
 * Create a new piece of dna identical to the model but belonging to <gen_unit>
 */
Dna::Dna(GeneticUnit* gen_unit, const Dna& model) :
    ae_string(model) {
  gen_unit_ = gen_unit;
  exp_m_ = gen_unit->exp_m();
  indiv_ = gen_unit->indiv();
}

/**
 * Create a new piece of dna identical to the parent's but belonging to
 * <gen_unit>
 */
Dna::Dna(GeneticUnit* gen_unit, Dna* const parent_dna) :
    ae_string(parent_dna->data_, parent_dna->length_) {
  gen_unit_ = gen_unit;
  exp_m_ = gen_unit->exp_m();
  indiv_ = gen_unit->indiv();
}

/**
 * Create a new piece of dna with sequence <seq> (of length <length>).
 * WARNING : <seq> will be used directly as the new dna sequence (it will not
 *           be copied), which means the caller must not delete it.
 */
// TODO(dpa) make seq a rvalue ref and set it to NULL ?
Dna::Dna(GeneticUnit* gen_unit, char* seq, int32_t length) :
    ae_string(seq, length, true) {
  gen_unit_ = gen_unit;
  exp_m_ = gen_unit->exp_m();
  indiv_ = gen_unit->indiv();
}

/**
 * Load a piece of dna from <backup_file>
 */
Dna::Dna(GeneticUnit* gen_unit, gzFile backup_file) :
    ae_string(backup_file) {
  gen_unit_ = gen_unit;
  exp_m_ = gen_unit->exp_m();
  indiv_ = gen_unit->indiv();
}

/**
 * Create a dna sequence from a text file
 */
Dna::Dna(GeneticUnit* gen_unit, char* organism_file_name) :
    ae_string(organism_file_name) {
  gen_unit_ = gen_unit;
  exp_m_ = gen_unit->exp_m();
  indiv_ = gen_unit->indiv();
}

// =================================================================
//                             Destructors
// =================================================================
Dna::~Dna() = default;

// =================================================================
//                         Non inline Accessors
// =================================================================
char* Dna::subsequence(int32_t from, int32_t to, Strand strand) const {
  char* subseq = NULL;

  from = Utils::mod(from, length_);
  to = Utils::mod(to, length_);

  if (strand == LEADING) {
    if (from < to) {
      subseq = new char[to - from + 1];
      subseq[to - from] = '\0';
      strncpy(subseq, &(data_[from]), to - from);
    }
    else {
      subseq = new char[length_ - from + to + 1];
      subseq[length_ - from + to] = '\0';
      strncpy(subseq, &(data_[from]), length_ - from);
      strncpy(&subseq[length_ - from], data_, to);
    }
  }
  else { // if (strand == LAGGING)
    if (from > to) {
      subseq = new char[from - to + 1];
      subseq[from - to] = '\0';

      for (int32_t i = 0; i < from - to; i++) {
        #ifdef BASE_2
        subseq[i] = (data_[from - 1 - i] == '1') ? '0' : '1';
        #elif BASE_4
        subseq[i] = get_complementary_base(data_[from - 1 - i]);
        #endif
      }
    }
    else {
      subseq = new char[from + length_ - to + 1];
      subseq[from + length_ - to] = '\0';

      for (int32_t i = 0; i < from; i++) {
        #ifdef BASE_2
        subseq[i] = (data_[from - 1 - i] == '1') ? '0' : '1';
        #elif BASE_4
        subseq[i] = get_complementary_base(data_[from - 1 - i]);
        #endif
      }
      for (int32_t i = 0; i < length_ - to; i++) {
        #ifdef BASE_2
        subseq[from + i] = (data_[length_ - 1 - i] == '1') ? '0' : '1';
        #elif BASE_4
        subseq[i] = get_complementary_base(data_[length_ - 1 - i]);
        #endif
      }
    }
  }

  return subseq;
}

// =================================================================
//                            Public Methods
// =================================================================
/// Perform mutations and record how many of them occurred
int32_t Dna::perform_mutations(int32_t parent_id) {
  int32_t nb_events = 0;

  if (indiv_->with_HT())
    nb_events += do_transfer(parent_id);

  if (indiv_->with_alignments())
    nb_events += do_rearrangements_with_align();
  else
    nb_events += do_rearrangements();

  nb_events += do_small_mutations();

  return nb_events;
}

int32_t Dna::do_small_mutations() {
  // ==============================================================
  //  1. Compute how many rearrangements this genome will undertake
  // ==============================================================
  //
  // Given the rate p (by nucl.) of insertion - for instance -, the number of
  // insertions we perform on the genome follows a binomial law B(n,p), with
  // n = genome length.

  int32_t nb_swi = indiv_->mut_prng_->
      binomial_random(length_, indiv_->point_mutation_rate());
  int32_t nb_ins = indiv_->mut_prng_->
      binomial_random(length_, indiv_->small_insertion_rate());
  int32_t nb_del = indiv_->mut_prng_->
      binomial_random(length_, indiv_->small_deletion_rate());
  int32_t nb_mut = nb_swi + nb_ins + nb_del;



  if (nb_mut > 0) {

/*    printf("%d -- Swi %d Ins %d Del %d\n",
           indiv()->id(),
           nb_swi,
           nb_ins,nb_del);*/
    if (!hasMutated) {
      hasMutated = true;
    }
  }


  // ====================================================
  //  2. Perform those small mutations in a random order
  // ====================================================
  //
  // We put the 'nb_small_mutations_' mutation events in an "urn".
  // Then we repeat a random drawing of one mutation event in this urn,
  // without replacement, until no mutation event is left in the urn.
  // Here is the "urn" we use at the beginning:
  //
  //     -----------------------------------------------------------
  //    | swi | swi | swi | ins | ins | ins | del | del | del | del |
  //     -----------------------------------------------------------
  //                      ^                 ^                       ^
  //                  nb_swi             nb_swi                   nb_swi
  //                                    +nb_ins                  +nb_ins
  //                                                             +nb_del
  //
  // Random draw of one mutation = random draw of one position in this "urn".
  // Given this position, we know what kind of mutation we have drawn.


  int32_t random_value;
  Mutation* mut = nullptr;

  for (int32_t i = nb_mut; i >= 1; i--) {
    random_value = indiv_->mut_prng_->random(i);

    if (random_value < nb_swi) {
      mut = do_switch();

      nb_swi--;  // updating the urn (no replacement!)...
    }
    else if (random_value < nb_swi + nb_ins) {
      mut = do_small_insertion();

      nb_ins--;
    }
    else { // (random_value >= nb_swi + nb_ins) => del
      mut = do_small_deletion();

      nb_del--;
    }

    // Record mutation in tree
    if (mut != NULL) {
      //indiv_->notifyObservers(MUTATION, mut);
        if (exp_m_->record_tree() || exp_m_->light_tree()) {
          indiv_->exp_m_->tree()->report_by_index(AeTime::time(),indiv_->grid_cell()->x() *
                                                                indiv_->exp_m()->grid_height()
                                                                + indiv_->grid_cell()->y())->dna_replic_report().add_mut(mut);
        }
      delete mut;
    }
  }

  return nb_mut;
}

int32_t Dna::do_rearrangements() {
  // ==============================================================
  //  1. Compute how many rearrangements this genome will undertake
  // ==============================================================
  //
  // Given the rate p (by nucl.) of duplication - for instance -, the number of
  // duplications we perform on the genome follows a binomial law B(n, p), with
  // n = genome length.
  int32_t nb_dupl = indiv_->mut_prng_->
      binomial_random(length_, indiv_->duplication_rate());
  int32_t nb_del = indiv_->mut_prng_->
      binomial_random(length_, indiv_->deletion_rate());
  int32_t nb_trans = indiv_->mut_prng_->
      binomial_random(length_, indiv_->translocation_rate());
  int32_t nb_inv = indiv_->mut_prng_->
      binomial_random(length_, indiv_->inversion_rate());
  int32_t nb_rear = nb_dupl + nb_del + nb_trans + nb_inv;




  if (nb_rear > 0) {
/*    printf("%d -- Dupli %d Large_Del %d Trans %d Inv %d\n",indiv()->id(),
           nb_dupl,nb_del,nb_trans,nb_inv);*/
    if (!hasMutated) {
      hasMutated = true;
    }
  }

  // ===================================================
  //  2. Perform those rearrangements in a random order
  // ===================================================
  //
  // We put the nb_rea rearrangements in an "urn". Then we repeat a random draw
  // of one rearrangement in this urn, without replacement, until no rearrange-
  // -ment is left in the urn. Here is the "urn" we use at the beginning:
  //
  //     ----------------------------------------------------------------
  //    | Dupl | Dupl | Del  Del | Del | Trans | Trans | Inv | Inv | Inv |
  //     ----------------------------------------------------------------
  //                  ^                ^               ^                 ^
  //               nb_dupl          nb_dupl         nb_dupl           nb_dupl
  //                               +nb_del         +nb_del           +nb_del
  //                                               +nb_trans         +nb_trans
  //                                                                 +nb_inv
  //
  // Random draw of one rearrangement = random draw of one position in this urn.
  // Given this position, we know what kind of rearrangement we have drawn.

  int32_t random_value;
  Mutation* mut = nullptr;

  for (int32_t i = nb_rear; i >= 1; i--) {
    random_value = indiv_->mut_prng_->random(i);

    if (random_value < nb_dupl) {
      mut = do_duplication();
      nb_dupl--;  // Updating the urn (no replacement!)...
    }
    else if (random_value < nb_dupl + nb_del) {
      mut = do_deletion();
      nb_del--;
    }
    else if (random_value < nb_dupl + nb_del + nb_trans) {
      mut = do_translocation();
      nb_trans--;
    }
    else {
      mut = do_inversion();
      nb_inv--;
    }

    // Record rearrangement in tree
    if (mut != NULL) {
//      indiv_->notifyObservers(MUTATION, mut);
        if (exp_m_->record_tree() || exp_m_->light_tree()) {
          indiv_->exp_m_->tree()->report_by_index(AeTime::time(),indiv_->grid_cell()->x() *
                                                                indiv_->exp_m()->grid_height()
                                                                + indiv_->grid_cell()->y())->dna_replic_report().add_mut(mut);
        }
        delete mut;
    }
  }
  return nb_rear;
}

int32_t Dna::do_rearrangements_with_align() {
  // Whether we look for a direct or indirect alignment
  bool direct_sense;
  // Determines the type of rearrangement that will be done if an alignment
  // is found
  double rand1 = 0.0;
  // Minimum alignment score needed to recombine (stochastic)
  int16_t needed_score;
  // Points defining the sequences between which we will look for an alignment
  int32_t seed1, seed2;

  // Indiv's Time To Live
  double ttl = 1.0;
  // Number of pairs of sequences we will try to align
  int32_t nb_pairs;
  // Keep trace of the original length of the genome
  int32_t genome_size = length_;
  int32_t nb_rearr = 0;

  Mutation* mut = nullptr;
  VisAVis* alignment = NULL;

  /////////////////////////////////////////////////////////////////////////////
  // For each pair of points to be tested
  // (i.e. while the organism is still "alive"),
  // 1) Draw a random sense (direct or indirect).
  // 2) Determine the minimum alignment score needed for a rear to occur.
  // 3) Test the existence of an alignment with a high enough score.
  // 4) If such an alignment was found, determine the type of rear to be
  //    performed and proceed (WARNING : translocations require another
  //    alignment to be found between the sequence to be translocated and
  //    the rest of the chromosome).
  // 5) If there was a change in the chromosome's length, update the
  //    individual's TTL and nb_pairs according to new genome size.
  // 6) If there was a rearrangement, we either save its record in the tree or
  //    delete it.
  /////////////////////////////////////////////////////////////////////////////
  nb_pairs = static_cast<int32_t>(
      ceil(ttl * length_ * indiv_->neighbourhood_rate()));
  for (; nb_pairs > 0; nb_pairs--) {
    /////////////////////////////////////////////////
    // 1) Draw a random sense (direct or indirect) //
    /////////////////////////////////////////////////
    // Determine whether we look for a direct or indirect alignment
    direct_sense = (indiv_->mut_prng_->random() < 0.5);
    // Determine the type of rearrangement to be done. This is an
    // anticipation on step 4) for optimization purpose (save computation
    // time if the type of rear is "none"
    rand1 = indiv_->mut_prng_->random();


    ///////////////////////////////////////////////////////////////////////////
    // 2) Determine the minimum alignment score needed for a rearrangement
    // to occur
    ///////////////////////////////////////////////////////////////////////////
    if (indiv_->align_fun_shape() == LINEAR) {
      needed_score = static_cast<int16_t>(
          ceil(indiv_->align_lin_min() +
               indiv_->mut_prng_->random() *
               (indiv_->align_lin_max() -
                indiv_->align_lin_min())));
    }
    else {
      // I want the probability of rearrangement for an alignment of score
      // <score> to be prob = 1 / (1 + exp(-(score-mean)/lambda))
      // The score needed for a rearrangement to take place with a given random
      // drawing is hence needed_score = ceil(-lambda * log(1/rand - 1) + mean)
      needed_score = static_cast<int16_t>(
          ceil(-indiv_->align_sigm_lambda() *
               log(1 / indiv_->mut_prng_->random() - 1) +
               indiv_->align_sigm_mean()));
      if (needed_score < 0) needed_score = 0;

      //~ <DEBUG>
      //~ FILE* tmp_file = fopen("scores.out", "a");
      //~ fprintf(tmp_file, "%"PRId16"\n", needed_score);
      //~ fclose(tmp_file);
      //~ </DEBUG>
    }

    // Determine where to look for an alignment (draw seeds)
    seed1 = indiv_->mut_prng_->random(length_);
    seed2 = indiv_->mut_prng_->random(length_);


    if (direct_sense) {
      if (rand1 >= indiv_->duplication_proportion() +
                   indiv_->deletion_proportion() +
                   indiv_->translocation_proportion()) {
        // rand1 corresponds to "no rearrangement" => Nothing to do
        continue;
      }

      ////////////////////////////////////////////////////////////////////
      // 3) Test the existence of an alignment with a high enough score //
      ////////////////////////////////////////////////////////////////////
      alignment = Alignment::search_alignment_direct(this, seed1,
                                                     this, seed2,
                                                     needed_score);

      if (alignment == NULL) {
        // No alignment found
        continue;
      }

      //~ printf("direct   needed_score : %"PRId32"\n", needed_score);

      ////////////////////////////////////////////////////////////////////////
      // 4) Determine the type of rearrangement to be performed and proceed //
      ////////////////////////////////////////////////////////////////////////
      if (rand1 < indiv_->duplication_proportion()) {
        // Remember the length of the segment to be duplicated and of the
        // genome before the duplication
        int32_t segment_length = Utils::mod(alignment->i_2() -
                                            alignment->i_1(),
                                            length_);
        int32_t gu_size_before = length_;
        int32_t gu_size_after = gu_size_before + segment_length;
        int32_t genome_size_before = indiv_->amount_of_dna();
        int32_t genome_size_after = genome_size_before + segment_length;

        if ((genome_size_after > indiv_->max_genome_length()) ||
            (gu_size_after > gen_unit_->max_gu_length())) {
          if (exp_m_->output_m()->is_logged(LOG_BARRIER)) {
            // Write an entry in the barrier log file
            fprintf(exp_m_->output_m()->log(LOG_BARRIER),
                    "%" PRId64 " %" PRId32 " DUPLICATION %" PRId32 " %" PRId32
                        " %" PRId32 " %" PRId32 "\n",
                    AeTime::time(), indiv_->grid_cell()->x() *
                                    indiv_->exp_m()->grid_height()
                                    + indiv_->grid_cell()->y(), segment_length, 0,
                    gu_size_before, genome_size_before);
          }
        }
        else {
          // Perform in situ (tandem) DUPLICATION
          do_duplication(alignment->i_1(),
                         alignment->i_2(),
                         alignment->i_2());

          // Report the duplication
          mut = new Duplication(alignment->i_1(),
                                alignment->i_2(),
                                alignment->i_2(),
                                segment_length, needed_score);
	  nb_rearr++;
          // Write a line in rearrangement logfile
          if (exp_m_->output_m()->is_logged(LOG_REAR)) {
            fprintf(exp_m_->output_m()->log(LOG_REAR),
                    "%" PRId64 " %" PRId32 " %" PRId8 " %" PRId32 " %" PRId32
                        " %" PRId16 "\n",
                    AeTime::time(), indiv_->grid_cell()->x() *
                                    indiv_->exp_m()->grid_height()
                                    + indiv_->grid_cell()->y(), int8_t(DUPL),
                    segment_length, genome_size_before, needed_score);
          }
        }
      }
      else if (rand1 < indiv_->duplication_proportion() +
                       indiv_->deletion_proportion()) {
        // Remember the length of the segment to be duplicated and of the
        // genome before the deletion
        int32_t segment_length = Utils::mod(alignment->i_2() -
                                            alignment->i_1() - 1,
                                            length_) + 1;
        int32_t gu_size_before = length_;
        int32_t gu_size_after = gu_size_before - segment_length;
        int32_t genome_size_before = indiv_->amount_of_dna();
        int32_t genome_size_after = genome_size_before - length_;

        if ((genome_size_after < indiv_->min_genome_length()) ||
            (gu_size_after < gen_unit_->min_gu_length())) {
          if (exp_m_->output_m()->is_logged(LOG_BARRIER)) {
            // Write an entry in the barrier log file
            fprintf(exp_m_->output_m()->log(LOG_BARRIER),
                    "%" PRId64 " %" PRId32 " DELETION %" PRId32 " %" PRId32
                        " %" PRId32 " %" PRId32 "\n",
                    AeTime::time(), indiv_->grid_cell()->x() *
                                    indiv_->exp_m()->grid_height()
                                    + indiv_->grid_cell()->y(), segment_length, 0,
                    gu_size_before, genome_size_before);
          }
        }
        else {
          // Perform DELETION
          do_deletion(alignment->i_1(), alignment->i_2());

          // Report the deletion
          mut = new Deletion(alignment->i_1(),
                             alignment->i_2(),
                             segment_length, needed_score);
	  nb_rearr++;

          // Write a line in rearrangement logfile
          if (exp_m_->output_m()->is_logged(LOG_REAR)) {
            fprintf(exp_m_->output_m()->log(LOG_REAR),
                    "%" PRId64 " %" PRId32 " %" PRId8 " %" PRId32 " %" PRId32
                        " %" PRId16 "\n",
                    AeTime::time(), indiv_->grid_cell()->x() *
                                    indiv_->exp_m()->grid_height()
                                    + indiv_->grid_cell()->y(), int8_t(DEL),
                    segment_length, genome_size_before, needed_score);
          }
        }
      }
      else {
        assert(rand1 < indiv_->duplication_proportion() +
                       indiv_->deletion_proportion() +
                       indiv_->translocation_proportion());

        // Perform TRANSLOCATION
        // Make sure the segment to be translocated doesn't contain OriC
        // TODO(dpa) is that still necessary?
        if (alignment->i_1() > alignment->i_2()) {
          alignment->swap();
        }

        // Remember the length of the segment to be translocated
        int32_t segment_length =
            Utils::mod(alignment->i_2() - alignment->i_1(), length_);

        // Extract the segment to be translocated
        GeneticUnit* translocated_segment =
            extract_into_new_GU(alignment->i_1(), alignment->i_2());

        // Look for an alignment between the segment to be translocated and
        // the rest of the genome
        bool direct_sense;
        int16_t needed_score_2;
        VisAVis* alignment_2 = NULL;
        int32_t seed1, seed2;
        nb_pairs = static_cast<int32_t>(
            ceil(ttl * length_ * indiv_->neighbourhood_rate()));
        for (; nb_pairs > 0; nb_pairs--) {
          direct_sense = (indiv_->mut_prng_->random() < 0.5);

          if (indiv_->align_fun_shape() == LINEAR) {
            needed_score_2 = static_cast<int16_t>(
                ceil(indiv_->align_lin_min() +
                     indiv_->mut_prng_->random() *
                     (indiv_->align_lin_max() -
                      indiv_->align_lin_min())));
          }
          else {
            needed_score_2 = static_cast<int16_t>(
                ceil(-indiv_->align_sigm_lambda() *
                     log(1 / indiv_->mut_prng_->random() - 1) +
                     indiv_->align_sigm_mean()));
            if (needed_score_2 < 0) needed_score_2 = 0;
          }

          seed1 = indiv_->mut_prng_->random(length_);
          seed2 = indiv_->mut_prng_->random(segment_length);

          if (direct_sense) {
            alignment_2 = Alignment::search_alignment_direct(
                this, seed1,
                translocated_segment->dna(), seed2,
                needed_score_2);
          }
          else { // if indirect
            alignment_2 = Alignment::search_alignment_indirect(
                this, seed1,
                translocated_segment->dna(), seed2,
                needed_score_2);
          }

          if (alignment_2 != NULL) {
            //~ printf("transloc needed_score : %"PRId32"\n", needed_score_2);
            break;
          }
        }


        // If an alignment was found between the segment to be translocated
        // and the rest of the genome, proceed to the translocation.
        // Else, replace the extracted segment at its former position
        // (cancel the translocation event).
        if (alignment_2 != NULL) {
          // Proceed to the translocation
          insert_GU(translocated_segment,
                    alignment_2->i_1(), alignment_2->i_2(),
                    (alignment_2->sense() == INDIRECT));

          // Report the translocation
          mut = new Translocation(alignment->i_1(),
                                  alignment->i_2(),
                                  alignment_2->i_1(),
                                  alignment_2->i_2(),
                                  segment_length,
                                  (alignment_2->sense() == INDIRECT),
                                  needed_score, needed_score_2);
	  nb_rearr++;

          // Write a line in rearrangement logfile
          if (exp_m_->output_m()->is_logged(LOG_REAR)) {
            fprintf(exp_m_->output_m()->log(LOG_REAR),
                    "%" PRId64 " %" PRId32 " %" PRId8 " %" PRId32 " %" PRId32
                        " %" PRId16 "\n",
                    AeTime::time(), indiv_->grid_cell()->x() *
                                    indiv_->exp_m()->grid_height()
                                    + indiv_->grid_cell()->y(), int8_t(TRANS),
                    segment_length, length_, needed_score_2);
          }
          delete alignment_2;
        }
        else {
          // Cancel the translocation
          // (re-place the extracted segment at its former position)
          insert_GU(translocated_segment, alignment->i_1(), 0, false);
        }

        delete translocated_segment;
      }

      delete alignment;
    }
    else { // if indirect
      if (rand1 >= indiv_->inversion_proportion()) {
        // rand1 corresponds to no rearrangement => Nothing to do
        continue;
      }

      //~ printf("indirect needed_score : %"PRId32"\n", needed_score);

      ////////////////////////////////////////////////////////////////////
      // 3) Test the existence of an alignment with a high enough score //
      ////////////////////////////////////////////////////////////////////
      alignment = Alignment::search_alignment_indirect(this, seed1,
                                                       this, seed2,
                                                       needed_score);

      if (alignment == NULL) {
        // No alignment found
        continue;
      }

      /////////////////////////////
      // 4) Proceed to inversion //
      /////////////////////////////
      // Make sure the segment to be inverted doesn't contain OriC
      if (alignment->i_1() > alignment->i_2()) {
        alignment->swap();
      }

      // Remember the length of the segment to be duplicated
      int32_t segment_length = Utils::mod(alignment->i_2() -
                                          alignment->i_1(),
                                          length_);

      // Proceed
      do_inversion(alignment->i_1(), alignment->i_2());

      // Report the inversion
      mut = new Inversion(alignment->i_1(),
                          alignment->i_2(),
                          segment_length, needed_score);
      nb_rearr++;
      // Write a line in rearrangement logfile
      if (exp_m_->output_m()->is_logged(LOG_REAR)) {
        fprintf(exp_m_->output_m()->log(LOG_REAR),
                "%" PRId64 " %" PRId32 " %" PRId8 " %" PRId32 " %" PRId32
                    " %" PRId16 "\n",
                AeTime::time(), indiv_->grid_cell()->x() *
                                    indiv_->exp_m()->grid_height()
                                    + indiv_->grid_cell()->y(), int8_t(INV),
                segment_length, length_, needed_score);
      }

      delete alignment;
    }


    ///////////////////////////////////////////////////////////////////////////
    // 5) If there was a change in the chromosome's length,
    //    update the individual's TTL and nb_pairs according to new genome size
    ///////////////////////////////////////////////////////////////////////////
    if (genome_size != length_) {
      ttl = (static_cast<double>(nb_pairs - 1)) /
            (static_cast<double>(genome_size)) /
            indiv_->neighbourhood_rate();
      genome_size = length_;
      nb_pairs = static_cast<int32_t>(
                     ceil(ttl * length_ * indiv_->neighbourhood_rate())) + 1;
    }

    ///////////////////////////////////////////////////////////////////////////
    // 6) If there was a rearrangement, we either save its record in the tree
    //    or delete it.
    ///////////////////////////////////////////////////////////////////////////
    if (mut != NULL) {
//      indiv_->notifyObservers(MUTATION, mut);
        if (exp_m_->record_tree() || exp_m_->light_tree()) {
          indiv_->exp_m_->tree()->report_by_index(AeTime::time(),indiv_->grid_cell()->x() *
                                                                indiv_->exp_m()->grid_height()
                                                                + indiv_->grid_cell()->y())->dna_replic_report().add_mut(mut);
        }
      delete mut;
    }
  }
  return nb_rearr;
}

int32_t Dna::do_transfer(int32_t parent_id) {
  Mutation* mut = nullptr;
  int32_t nb_transfer = 0;

  if (indiv_->mut_prng()->random() < indiv_->HT_ins_rate()) {
    mut = do_ins_HT(parent_id);
    if (mut != nullptr) {
//      indiv_->notifyObservers(MUTATION, mut);
        if (exp_m_->record_tree() || exp_m_->light_tree()) {
          indiv_->exp_m_->tree()->report_by_index(AeTime::time(),indiv_->grid_cell()->x() *
                                                                indiv_->exp_m()->grid_height()
                                                                + indiv_->grid_cell()->y())->dna_replic_report().add_mut(mut);
        }
      nb_transfer++;
      delete mut;
    }
  }

  if (indiv_->mut_prng()->random() < indiv_->HT_repl_rate()) {
    mut = do_repl_HT(parent_id);
    if (mut != nullptr) {
      //indiv_->notifyObservers(MUTATION, mut);
        if (exp_m_->record_tree() || exp_m_->light_tree()) {
          indiv_->exp_m_->tree()->report_by_index(AeTime::time(),indiv_->grid_cell()->x() *
                                                                indiv_->exp_m()->grid_height()
                                                                + indiv_->grid_cell()->y())->dna_replic_report().add_mut(mut);
        }
      nb_transfer++;
      delete mut;
    }
  }
  return nb_transfer;
}

PointMutation* Dna::do_switch() {
  PointMutation* mut = nullptr;

  int32_t pos = indiv_->mut_prng_->random(length_);
  
  #ifdef BASE_2
  if (do_switch(pos)) {
    // Report the mutation
    mut = new PointMutation(pos);
  }
  #elif BASE_4
  char choice = static_cast<char>('0' + indiv_->mut_prng_->random(NB_BASE - 1));
  
  if (do_switch(pos,choice)) {
    // Report the mutation
    mut = new PointMutation(pos,choice);
  }
  #endif
  
  return mut;
}

SmallInsertion* Dna::do_small_insertion() {
  SmallInsertion* mut = nullptr;

  // Determine the position and size of the small insertion
  int32_t pos = indiv_->mut_prng_->random(length_);
  int16_t nb_insert;
  if (indiv_->max_indel_size() == 1) {
    nb_insert = 1;
  }
  else {
    nb_insert = 1 + indiv_->mut_prng_->random(indiv_->max_indel_size());
    // <nb_insert> must be in [1 ; max_indel_size]
  }

  // Check that the insertion won't throw the genome size over the limit
  if ((indiv_->amount_of_dna() + nb_insert >
       indiv_->max_genome_length()) ||
      (length_ + nb_insert > gen_unit_->max_gu_length())) {
    if (exp_m_->output_m()->is_logged(LOG_BARRIER)) {
      // Write an entry in the barrier log file
      fprintf(exp_m_->output_m()->log(LOG_BARRIER),
              "%" PRId64 " %" PRId32 " S_INS %" PRId32 " %" PRId32 " %" PRId32
                  " %" PRId32 "\n",
              AeTime::time(), indiv_->grid_cell()->x() *
                                    indiv_->exp_m()->grid_height()
                                    + indiv_->grid_cell()->y(), nb_insert, 0, length_,
              indiv_->amount_of_dna());
    }

    return NULL;
  }

  // Prepare the sequence to be inserted
  char* inserted_seq = new char[nb_insert + 1];
  char inserted_char;
  for (int16_t j = 0; j < nb_insert; j++) {
    inserted_char = static_cast<char>('0' + indiv_->mut_prng_->random(NB_BASE));
    inserted_seq[j] = inserted_char;
  }
  inserted_seq[nb_insert] = '\0';

  // Proceed to the insertion and report it
  if (do_small_insertion(pos, nb_insert, inserted_seq)) {
    // Report the insertion
    mut = new SmallInsertion(pos, nb_insert, inserted_seq);
  }

  // Delete the sequence
  delete[] inserted_seq;

  return mut;
}

SmallDeletion* Dna::do_small_deletion() {
  SmallDeletion* mut = nullptr;

  // Determine the position and size of the small deletion
  int32_t pos = indiv_->mut_prng_->random(length_);
  int16_t nb_del;
  if (indiv_->max_indel_size() == 1) {
    nb_del = 1;
  }
  else {
    nb_del = 1 + indiv_->mut_prng_->random(indiv_->max_indel_size());
    // <nb_del> must be in [1 ; max_indel_size]
  }

  // Check that the insertion will shrink neither the genome nor the GU size
  // under their respective limit
  if ((indiv_->amount_of_dna() - nb_del <
       indiv_->min_genome_length()) ||
      (length_ - nb_del < gen_unit_->min_gu_length())) {
    if (exp_m_->output_m()->is_logged(LOG_BARRIER)) {
      // Write an entry in the barrier log file
      fprintf(exp_m_->output_m()->log(LOG_BARRIER),
              "%" PRId64 " %" PRId32 " S_DEL %" PRId32 " %" PRId32
                  " %" PRId32 " %" PRId32 "\n",
              AeTime::time(), indiv_->grid_cell()->x() *
                                    indiv_->exp_m()->grid_height()
                                    + indiv_->grid_cell()->y(), nb_del, 0, length_,
              indiv_->amount_of_dna());
    }

    return nullptr;
  }

  if (do_small_deletion(pos, nb_del)) {
    mut = new SmallDeletion(pos, nb_del);
  }

  return mut;
}

bool Dna::do_switch(int32_t pos
#ifdef BASE_4
, char choice
#endif
) {
  #ifdef BASE_2
  // Perform the mutation
  if (data_[pos] == '0') data_[pos] = '1';
  else data_[pos] = '0';
  #elif BASE_4
  // Perform the mutation
  char current_nucleotide = data_[pos];
  data_[pos] = choice; // generated in range [0; NB_BASE-1[
  // exclude current_nucleotide by shifting generated nucleotide by 1 when equal or above
  if (data_[pos] >= current_nucleotide) data_[pos] += 1;
  #endif

  // Remove promoters containing the switched base
  gen_unit_->remove_promoters_around(pos, Utils::mod(pos + 1, length_));

  // Look for potential new promoters containing the switched base
  if (length_ >= PROM_SIZE)
    gen_unit_->look_for_new_promoters_around(pos, Utils::mod(pos + 1, length_));

  return true;
}

bool Dna::do_small_insertion(int32_t pos, int32_t nb_insert, char* seq) {
  // Check genome size limit
  assert(length_ + nb_insert <= gen_unit_->max_gu_length());
  assert(indiv_->amount_of_dna() + nb_insert <=
         indiv_->max_genome_length());

  // Remove the promoters that will be broken
  gen_unit_->remove_promoters_around(pos);

  // Insert the sequence
  insert(pos, seq, nb_insert);

  // Look for new promoters
  if (length_ >= PROM_SIZE) {
    if (length_ - nb_insert < PROM_SIZE) {
      // Special case where the genome was smaller than a promoter before the
      // insertion and greater than (or as big as) a promoter after the
      // insertion.
      // In that case, we must look for new promoters thoroughly on the whole
      // genome using locate_promoters
      gen_unit_->locate_promoters();
    }
    else {
      gen_unit_->move_all_promoters_after(pos, nb_insert);
      gen_unit_->look_for_new_promoters_around(pos, Utils::mod(pos + nb_insert,
                                                               length_));
    }
  }

  return true;
}

bool Dna::do_small_deletion(int32_t pos, int16_t nb_del) {
  // Check genome size limit
  assert(length_ - nb_del >= gen_unit_->min_gu_length());
  assert(indiv_->amount_of_dna() - nb_del >=
         indiv_->min_genome_length());

  // Remove promoters containing at least one nucleotide from the sequence to
  // delete
//  printf("Remove %d %d (%d %d) %d\n",pos,Utils::mod(pos + nb_del, length()),pos,nb_del,length());
//  fflush(stdout);

  gen_unit_->remove_promoters_around(pos, Utils::mod(pos + nb_del, length_));

  // Do the deletion and update promoter list
  if (pos + nb_del <= length_) { // the deletion does not contain the origin of
    // replication
    // Do the deletion
    remove(pos, pos + nb_del);

    // Update promoter list
    if (length_ >= PROM_SIZE) {
      gen_unit_->move_all_promoters_after(pos, -nb_del);
      gen_unit_->look_for_new_promoters_around(Utils::mod(pos, length_));
    }
  }
  else { // the deletion contains the origin of replication
    // Do the deletion
    int32_t nb_del_at_pos_0 = nb_del - length_ + pos;
    remove(pos, length_);
    remove(0, nb_del_at_pos_0);
    pos -= nb_del_at_pos_0;

    // Update promoter list
    if (length_ >= PROM_SIZE) {
      gen_unit_->move_all_promoters_after(0, -nb_del_at_pos_0);
      gen_unit_->look_for_new_promoters_around(0);
    }
  }

  return true;
}

Duplication* Dna::do_duplication() {
  Duplication* mut = nullptr;
  if (length_ == 1)
    {
      printf("*** genome of size 1 ; duplication not done *** \n");
      return mut;
    }
  int32_t pos_1, pos_2, pos_3;
  pos_1 = indiv_->mut_prng_->random(length_);
  pos_2 = indiv_->mut_prng_->random(length_);

  // TODO: why are complete duplications forbidden ?
  while (pos_2 == pos_1)
    {
      pos_2 = indiv_->mut_prng_->random(length_);
    }
  pos_3 = indiv_->mut_prng_->random(length_);

  // Remember the length of the segment to be duplicated and of the former
  // genome
  int32_t segment_length = Utils::mod(pos_2 - pos_1 - 1, length_) + 1;
  int32_t gu_size_before = length_;
  int32_t gu_size_after = gu_size_before + segment_length;
  int32_t genome_size_before = indiv_->amount_of_dna();
  int32_t genome_size_after = genome_size_before + segment_length;
  if ((gu_size_after > gen_unit_->max_gu_length()) ||
      (genome_size_after > indiv_->max_genome_length())) {
    if (exp_m_->output_m()->is_logged(LOG_BARRIER)) {
      // Write an entry in the barrier log file
      fprintf(exp_m_->output_m()->log(LOG_BARRIER),
              "%" PRId64 " %" PRId32 " DUPLICATION %" PRId32 " %" PRId32
                  " %" PRId32 " %" PRId32 "\n",
              AeTime::time(), indiv_->grid_cell()->x() *
                                    indiv_->exp_m()->grid_height()
                                    + indiv_->grid_cell()->y(), segment_length, 0,
              gu_size_before, genome_size_before);
    }
  }
  else {
    // Perform the duplication
    do_duplication(pos_1, pos_2, pos_3);

    // Report the duplication
    mut = new Duplication(pos_1, pos_2, pos_3, segment_length);

    // Write a line in rearrangement logfile
    if (exp_m_->output_m()->is_logged(LOG_REAR)) {
      fprintf(exp_m_->output_m()->log(LOG_REAR),
              "%" PRId64 " %" PRId32 " %" PRId8 " %" PRId32 " %" PRId32 "\n",
              AeTime::time(), indiv_->grid_cell()->x() *
                                    indiv_->exp_m()->grid_height()
                                    + indiv_->grid_cell()->y(), int8_t(DUPL),
              segment_length, genome_size_before);
    }
  }

  return mut;
}

Deletion* Dna::do_deletion() {
  Deletion* mut = nullptr;

  int32_t pos_1, pos_2;
  if (length_ == 1)
    {
      printf("*** genome of size 1 ; deletion not done *** \n");
      return mut;
    }
  pos_1 = indiv_->mut_prng_->random(length_);
  pos_2 = indiv_->mut_prng_->random(length_);
  while (pos_2 == pos_1)
    {
      pos_2 = indiv_->mut_prng_->random(length_);
    }
  // Remember the length of the segment to be deleted and of the genome
  // before the deletion
  int32_t segment_length = Utils::mod(pos_2 - pos_1 - 1, length_) + 1;
  int32_t gu_size_before = length_;
  int32_t gu_size_after = gu_size_before - segment_length;
  int32_t genome_size_before = indiv_->amount_of_dna();
  int32_t genome_size_after = genome_size_before - segment_length;

  // TODO : please, never commit that again.
//  if(genome_size_after < indiv_->min_genome_length())
//    {
//      printf("the shit hits the fan in Dna.cpp (%d %d)\n",genome_size_after,indiv_->min_genome_length());
//      exit(0);
//    }

  if ((gu_size_after < gen_unit_->min_gu_length()) ||
      (genome_size_after < indiv_->min_genome_length())) {
    if (exp_m_->output_m()->is_logged(LOG_BARRIER)) {
      // Write an entry in the barrier log file
      fprintf(exp_m_->output_m()->log(LOG_BARRIER),
              "%" PRId64 " %" PRId32 " DELETION %" PRId32 " %" PRId32
                  " %" PRId32 " %" PRId32 "\n",
              AeTime::time(), indiv_->grid_cell()->x() *
                                    indiv_->exp_m()->grid_height()
                                    + indiv_->grid_cell()->y(), segment_length, 0,
              gu_size_before, genome_size_before);
    }
  }
  else {
    // Perform the deletion
    do_deletion(pos_1, pos_2);

    // Report the deletion
    mut = new Deletion(pos_1, pos_2, segment_length);

    // Write a line in rearrangement logfile
    if (exp_m_->output_m()->is_logged(LOG_REAR)) {
      fprintf(exp_m_->output_m()->log(LOG_REAR),
              "%" PRId64 " %" PRId32 " %" PRId8 " %" PRId32 " %" PRId32 "\n",
              AeTime::time(), indiv_->grid_cell()->x() *
                                    indiv_->exp_m()->grid_height()
                                    + indiv_->grid_cell()->y(), int8_t(DEL), segment_length,
              genome_size_before);
    }
  }

  return mut;
}

Translocation* Dna::do_translocation() {
  Translocation* mut = nullptr;

  int32_t pos_1, pos_2, pos_3, pos_4;
  int32_t segment_length;
  bool invert;

  if (indiv_->allow_plasmids()) {
    // -----------------------------------------------------------------
    // WARNING : This is only valid when there is only 1 plasmid allowed
    // -----------------------------------------------------------------
    int32_t pos_1_rel, pos_2_rel, pos_3_rel, pos_4_rel;

    Individual* indiv = indiv_;
    const GeneticUnit* chromosome = &indiv->genetic_unit_list().front();
    const GeneticUnit* plasmid =
        &*std::next(indiv->genetic_unit_list().begin());
    int32_t chrom_length = chromosome->dna()->length();
    int32_t total_amount_of_dna = indiv->amount_of_dna();

    // 1) What sequence are we translocating?
    pos_1_rel = indiv_->mut_prng_->random(length_);
    pos_2_rel = indiv_->mut_prng_->random(length_);

    int32_t segment_length = Utils::mod(pos_2_rel - pos_1_rel, length_);

    pos_3_rel =
        Utils::mod(pos_1_rel + indiv_->mut_prng_->random(segment_length),
                   length_);

    if (gen_unit_ == chromosome) {
      pos_1 = pos_1_rel;
      pos_2 = pos_2_rel;
      pos_3 = pos_3_rel;
    }
    else { // (gen_unit_ == plasmid)
      pos_1 = pos_1_rel + chrom_length;
      pos_2 = pos_2_rel + chrom_length;
      pos_3 = pos_3_rel + chrom_length;
    }


    // 2) Where are we translocating it?
    pos_4 = indiv_->mut_prng_->random(total_amount_of_dna - segment_length);

    if (gen_unit_ == chromosome) {
      if (pos_1 <= pos_2) {
        if (pos_4 >= pos_1) {
          pos_4 += segment_length;
        }
      }
      else {
        if (pos_4 >= chrom_length - segment_length) {
          pos_4 += segment_length;
        }
        else {
          pos_4 += pos_2;
        }
      }
      if (pos_4 >= chrom_length) {
        pos_4_rel = pos_4 - chrom_length;
      }
      else {
        pos_4_rel = pos_4;
      }
    }
    else { // (gen_unit_ == plasmid)
      if (pos_1 <= pos_2) {
        if (pos_4 >= pos_1) {
          pos_4 += segment_length;
        }
      }
      else {
        if (pos_4 >= chrom_length) {
          pos_4 += pos_2_rel;
        }
      }

      if (pos_4 >= chrom_length) {
        pos_4_rel = pos_4 - chrom_length;
      }
      else {
        pos_4_rel = pos_4;
      }
    }

    invert = (indiv_->mut_prng_->random(2) == 0);


    // If inter GU translocation
    if ((gen_unit_ == chromosome && pos_4 >= chrom_length) ||
        (gen_unit_ == plasmid && pos_4 < chrom_length)) {
      if (do_inter_GU_translocation(pos_1_rel, pos_2_rel,
                                    pos_3_rel, pos_4_rel, invert)) {
        // Report the translocation
        mut = new Translocation(pos_1_rel, pos_2_rel,
                                pos_3_rel, pos_4_rel,
                                segment_length, invert);

        // Write a line in rearrangement logfile
        if (exp_m_->output_m()->is_logged(LOG_REAR)) {
          fprintf(exp_m_->output_m()->log(LOG_REAR),
                  "%" PRId64 " %" PRId32 " %" PRId8 " %" PRId32 " %" PRId32
                      "\n",
                  AeTime::time(), indiv_->grid_cell()->x() *
                                    indiv_->exp_m()->grid_height()
                                    + indiv_->grid_cell()->y(), int8_t(TRANS),
                  segment_length, length_);
        }
      }
    }
    else {
      if (do_translocation(pos_1_rel, pos_2_rel, pos_3_rel, pos_4_rel,
                           invert)) { // NOLINT(whitespace/braces)
        mut = new Translocation(pos_1_rel, pos_2_rel,
                                pos_3_rel, pos_4_rel,
                                segment_length, invert);

        // Write a line in rearrangement logfile
        if (exp_m_->output_m()->is_logged(LOG_REAR)) {
          fprintf(exp_m_->output_m()->log(LOG_REAR),
                  "%" PRId64 " %" PRId32 " %" PRId8 " %" PRId32 " %" PRId32
                      "\n",
                  AeTime::time(), indiv_->grid_cell()->x() *
                                    indiv_->exp_m()->grid_height()
                                    + indiv_->grid_cell()->y(), int8_t(TRANS),
                  segment_length, length_);
        }
      }
    }
  }
  else { // (! ae_common::params->allow_plasmids())
    pos_1 = indiv_->mut_prng_->random(length_);
    pos_2 = indiv_->mut_prng_->random(length_);
    if (pos_1 == pos_2) return NULL;

    // As it is commented in do_translocation(int32_t pos_1, int32_t pos_2,
    // int32_t pos_3, int32_t pos_4, bool invert), translocating segment
    // [pos_1, pos_2] is the same as translocating segment [pos_2, pos_1]
    // Since OriC must be at position 0, we will always translocate segment
    // [pos_1, pos_2] with pos_1 < pos_2
    if (pos_1 > pos_2) Utils::exchange(pos_1, pos_2);

    segment_length = pos_2 - pos_1;

    // Generate a position between pos_1 and pos_2
    pos_3 = pos_1 + indiv_->mut_prng_->random(segment_length);

    // Generate a position that is NOT between pos_1 and pos_2
    pos_4 = indiv_->mut_prng_->random(length_ - segment_length);
    if (pos_4 >= pos_1) pos_4 += segment_length;

    invert = (indiv_->mut_prng_->random(2) == 0);

    if (do_translocation(pos_1, pos_2, pos_3, pos_4, invert)) {
      // Report the translocation
      mut = new Translocation(pos_1, pos_2, pos_3, pos_4,
                              segment_length, invert);

      // Write a line in rearrangement logfile
      if (exp_m_->output_m()->is_logged(LOG_REAR)) {
        fprintf(exp_m_->output_m()->log(LOG_REAR),
                "%" PRId64 " %" PRId32 " %" PRId8 " %" PRId32 " %" PRId32 "\n",
                AeTime::time(), indiv_->grid_cell()->x() *
                                    indiv_->exp_m()->grid_height()
                                    + indiv_->grid_cell()->y(), int8_t(TRANS),
                segment_length, length_);
      }
    }
  }

  return mut;
}

Inversion* Dna::do_inversion() {
  Inversion* mut = nullptr;

  int32_t segment_length;
  int32_t pos_1, pos_2;
  if (length_ == 1)
    {
      printf("*** genome of size 1 ; inversion not done *** \n");
      return mut;
    }
  pos_1 = indiv_->mut_prng_->random(length_);
  pos_2 = indiv_->mut_prng_->random(length_);
  while (pos_2 == pos_1)
    {
      pos_2 = indiv_->mut_prng_->random(length_);
    }

  if (pos_1 == pos_2) return NULL; // Invert everything <=> Invert nothing!

  // Invert the segment that don't contain OriC
  if (pos_1 > pos_2) Utils::exchange(pos_1, pos_2);

  segment_length = pos_2 - pos_1;

  if (do_inversion(pos_1, pos_2)) {
    // Report the inversion
    mut = new Inversion(pos_1, pos_2, segment_length);

    // Write a line in rearrangement logfile
    if (exp_m_->output_m()->is_logged(LOG_REAR)) {
      fprintf(exp_m_->output_m()->log(LOG_REAR),
              "%" PRId64 " %" PRId32 " %" PRId8 " %" PRId32 " %" PRId32 "\n",
              AeTime::time(), indiv_->grid_cell()->x() *                                     indiv_->exp_m()->grid_height()                                     + indiv_->grid_cell()->y(), int8_t(INV),
              segment_length, length_);
    }
  }

  return mut;
}

Mutation* Dna::do_insertion(const char* seq_to_insert,
                            int32_t seq_length /*= -1*/) {
  Mutation* mut = nullptr;

  Utils::ExitWithDevMsg("Not implemented yet", __FILE__, __LINE__);

//  // Compute seq_length if not known
//  if (seq_length == -1) {
//    seq_length = strlen(seq_to_insert);
//  }
//
//  // Where to insert the sequence
//  int32_t pos = indiv_->mut_prng_->random(length_);
//
//  if (do_insertion(pos, seq_to_insert, seq_length)) {
//    // Report the insertion
//    mut = new Mutation();
//    mut->report_insertion(pos, seq_length, seq_to_insert);
//  }

  return mut;
}


bool Dna::do_duplication(int32_t pos_1, int32_t pos_2, int32_t pos_3) {
  //printf("%d %d -- Do duplication is %d %d %d -- %d\n",indiv_->id(), indiv()->parent_id_, pos_1,pos_2,pos_3,length());
  // Duplicate segment [pos_1; pos_2[ and insert the duplicate before pos_3
  char* duplicate_segment = NULL;
  int32_t seg_length;

  if (pos_1 < pos_2) {
    //
    //       pos_1         pos_2                   -> 0-
    //         |             |                   -       -
    // 0--------------------------------->      -         -
    //         ===============                  -         - pos_1
    //           tmp (copy)                      -       -
    //                                             -----      |
    //                                             pos_2    <-'
    //

    seg_length = pos_2 - pos_1;
    duplicate_segment = new char[seg_length + 1];
    memcpy(duplicate_segment, &data_[pos_1], seg_length);
    duplicate_segment[seg_length] = '\0';
  }
  else { // if (pos_1 >= pos_2)
    // The segment to duplicate includes the origin of replication.
    // The copying process will be done in two steps.
    //
    //                                            ,->
    //    pos_2                 pos_1            |      -> 0-
    //      |                     |                   -       - pos_2
    // 0--------------------------------->     pos_1 -         -
    // ======                     =======            -         -
    //  tmp2                        tmp1              -       -
    //                                                  -----
    //
    //

    int32_t tmp1_len = length_ - pos_1;
    int32_t tmp2_len = pos_2;
    seg_length = tmp1_len + tmp2_len;
    duplicate_segment = new char[seg_length + 1];
    memcpy(duplicate_segment, &data_[pos_1], tmp1_len);     // Copy tmp1
    memcpy(&duplicate_segment[tmp1_len], data_, tmp2_len);  // Copy tmp2
    duplicate_segment[seg_length] = '\0';
  }
  // Create a copy of the promoters beared by the segment to be duplicated
  // (they will be inserted in the individual's RNA list later)
  Promoters2Strands duplicated_promoters = {{},
                                            {}};

  gen_unit_->duplicate_promoters_included_in(pos_1, pos_2,
                                             duplicated_promoters);

  gen_unit_->remove_promoters_around(pos_3);

  insert(pos_3, duplicate_segment, seg_length);

  if (length_ >= PROM_SIZE) {
    if (length_ - seg_length < PROM_SIZE) {
      // Special case where the genome was smaller than a promoter before
      // the insertion and greater than (or as big as) a promoter after the
      // insertion.
      // In that case, we must look for new promoters thoroughly on the whole
      // genome using locate_promoters
      gen_unit_->locate_promoters();
    }
    else {
      gen_unit_->move_all_promoters_after(pos_3, seg_length);

      gen_unit_->insert_promoters_at(duplicated_promoters, pos_3);

      gen_unit_->look_for_new_promoters_around(pos_3);
      gen_unit_->look_for_new_promoters_around(pos_3 + seg_length);
    }
  }

  delete[] duplicate_segment;

  return true;
}

bool Dna::do_deletion(int32_t pos_1, int32_t pos_2) {
  //printf("%d %d -- DO DELETION is %d %d -- %d\n",indiv_->id(), indiv()->parent_id_,pos_1,pos_2,length());
  // printf("LARGE DELETION %d %d (%d): %d %d\n",indiv()->grid_cell()->x(),indiv()->grid_cell()->y(),
  //        indiv_->grid_cell_->x()*indiv_->exp_m_->world()->height()+indiv_->grid_cell_->y(),
  //       pos_1,pos_2);


// Delete segment going from pos_1 (included) to pos_2 (excluded)
  if (pos_1 < pos_2) {
    //
    //       pos_1         pos_2                   -> 0-
    //         |             |                   -       -
    // 0--------------------------------->      -         -
    //         ===============                  -         - pos_1
    //           tmp (copy)                      -       -
    //                                             -----      |
    //                                             pos_2    <-'
    //

    int32_t segment_length = pos_2 - pos_1;

    // Remove promoters containing at least one nucleotide from the sequence
    // to delete
    gen_unit_->remove_promoters_around(pos_1, pos_2);

    // Delete the sequence between pos_1 and pos_2
    remove(pos_1, pos_2);

    // Update promoter list
    if (length_ >= PROM_SIZE) {
      gen_unit_->move_all_promoters_after(pos_1, -segment_length);

      gen_unit_->look_for_new_promoters_around(pos_1);
    }
  }
  else { // if (pos_1 >= pos_2)
    // The segment to delete includes the origin of replication.
    // The deletion process will be done in two steps.
    //
    //                                            ,->
    //    pos_2                 pos_1            |      -> 0-
    //      |                     |                   -       - pos_2
    // 0--------------------------------->     pos_1 -         -
    // =====                      =======            -         -
    //  tmp2                        tmp1              -       -
    //                                                  -----
    //
    //

    // int32_t segment_length = length_ + pos_2 - pos_1; //useless variable

    // Remove promoters containing at least one nucleotide from the sequence
    // to delete
    gen_unit_->remove_promoters_around(pos_1, pos_2);

    // Delete the sequence between pos_1 and pos_2
    remove(pos_1, length_); // delete tmp1 from genome
    remove(0, pos_2);       // delete tmp2 from genome

    // Update promoter list
    if (length_ >= PROM_SIZE) {
      gen_unit_->move_all_promoters_after(0, -pos_2);
      gen_unit_->look_for_new_promoters_around(0);
    }
  }

  return true;
}

bool Dna::do_translocation(int32_t pos_1, int32_t pos_2, int32_t pos_3,
                           int32_t pos_4, bool invert) {
    //printf("%d %d -- Do Translocation %d %d %d %d %d\n",indiv_->id(), indiv()->parent_id_,pos_1,pos_2,pos_3,pos_4,invert);
  // Provided that OriC must be at position 0
  //
  //    1) Note that in Case 1 (without inversion), whichever position
  //       comes first, translocating segment [pos_1->pos_2] to pos_4
  //       through pos_3 is always equivalent to rearrange the sequences from
  //       an ABCDE order to ADCBE
  //    2) In Case 2, depending on which position comes first, we may have the
  //       following rearrangements :
  //       (ABCDE => ADB'C'E) or (ABCDE => AC'D'BE)
  //       where X' stands for "inverted X"
  //
  //  Case 1 : Without inversion
  //
  //         A      B        C       D       E
  //      |----->=======[>=======>-------[>-------|
  //          pos_1   pos_3    pos_2   pos_4
  //                         |
  //                         V
  //         A      D        C       B        E
  //      |----->-------[>=======>=======[>-------|
  //
  //
  //         A      B        C       D       E
  //      |=====>-------[>------->=======[>=======|
  //          pos_2   pos_4    pos_1   pos_3
  //                         |
  //                         V
  //         A      D        C       B        E
  //      |=====>=======[>------->-------[>=======|
  //
  //
  //         A      B        C       D       E
  //      |====[>========>-------[>------->=======|
  //          pos_3    pos_2    pos_4   pos_1
  //                         |
  //                         V
  //         A       D       C        B        E
  //      |=====[>------->-------[>=======[>=======|
  //
  //
  //         A      B        C       D       E
  //      |----[>-------->=======[>=======>-------|
  //          pos_4    pos_1    pos_3   pos_2
  //                         |
  //                         V
  //         A       D       C        B        E
  //      |-----[>=======>=======[>-------[>-------|
  //
  //  Case 2 : With inversion
  //
  //    Case 2.A
  //
  //         A      B        C       D        E
  //      |----->=======[>=======>-------<]-------|
  //          pos_1   pos_3    pos_2   pos_4
  //                         |
  //                         V
  //         A      D        B'      C'       E
  //      |----->-------<]=======<=======<]-------|
  //
  //
  //         A      B        C       D       E
  //      |=====>-------[>------->=======<]=======|
  //          pos_2   pos_4    pos_1   pos_3
  //                         |
  //                         V
  //         A      D        B'      C'       E
  //      |=====>=======<]-------<-------<]=======|
  //
  //    Case 2.B
  //
  //         A      B        C       D       E
  //      |====[>========>-------<]------->=======|
  //          pos_3    pos_2    pos_4   pos_1
  //                         |
  //                         V
  //         A       C'      D'       B       E
  //      |=====[>-------<-------[>=======>=======|
  //
  //
  //         A      B        C       D       E
  //      |----<]-------->=======[>=======>-------|
  //          pos_4    pos_1    pos_3   pos_2
  //                         |
  //                         V
  //         A       C'      D'       B       E
  //      |-----<]=======>=======<]------->-------|



  // Determine which position comes first and do the corresponding rearrangement
  // TODO(dpa) use min from std
  int32_t pos_min = Utils::min(pos_1,
                               Utils::min(pos_2, Utils::min(pos_3, pos_4)));

  if (not invert) {
    if (pos_min == pos_1) {
      ABCDE_to_ADCBE(pos_1, pos_3, pos_2, pos_4);
    }
    else if (pos_min == pos_2) {
      ABCDE_to_ADCBE(pos_2, pos_4, pos_1, pos_3);
    }
    else if (pos_min == pos_3) {
      ABCDE_to_ADCBE(pos_3, pos_2, pos_4, pos_1);
    }
    else { // if (pos_min == pos_4)
      ABCDE_to_ADCBE(pos_4, pos_1, pos_3, pos_2);
    }
  }
  else { // invert
    if (pos_min == pos_1) {
      ABCDE_to_ADBpCpE(pos_1, pos_3, pos_2, pos_4);
    }
    else if (pos_min == pos_2) {
      ABCDE_to_ADBpCpE(pos_2, pos_4, pos_1, pos_3);
    }
    else if (pos_min == pos_3) {
      ABCDE_to_ACpDpBE(pos_3, pos_2, pos_4, pos_1);
    }
    else { // if (pos_min == pos_4)
      ABCDE_to_ACpDpBE(pos_4, pos_1, pos_3, pos_2);
    }
  }

  return true;
}

bool Dna::do_inter_GU_translocation(int32_t pos_1_rel, int32_t pos_2_rel,
                                    int32_t pos_3_rel, int32_t pos_4_rel,
                                    bool invert) {
  // TODO(???) check GU lengths according to positions and size limit
  int32_t segment_length = Utils::mod(pos_2_rel - pos_1_rel, length_);

  if (pos_1_rel == pos_2_rel) { // TODO(???) shouldn't that raise an error?
    return false;
  }

  // Do not allow translocation if it would decrease the size of the origin
  // GU below a given threshold
  if ((length_ - segment_length) < gen_unit_->min_gu_length()) {
    if (exp_m_->output_m()->is_logged(LOG_BARRIER)) {
      // Write an entry in the barrier log file
      fprintf(exp_m_->output_m()->log(LOG_BARRIER),
              "%" PRId64 " %" PRId32 " TRANS %" PRId32 " %" PRId32 " %" PRId32
                  " %" PRId32 "\n",
              AeTime::time(), indiv_->grid_cell()->x() *
                                    indiv_->exp_m()->grid_height()
                                    + indiv_->grid_cell()->y(), segment_length, 0, length_,
              indiv_->amount_of_dna());
    }
    return false;
  }

  //
  const GeneticUnit& chromosome = indiv_->genetic_unit(0);
  const GeneticUnit& plasmid = indiv_->genetic_unit(1);
  // TODO(vld) (2015-02-23): check if this == is sound
  const GeneticUnit& destination_GU = (gen_unit_ == &chromosome) ?
                                      plasmid :
                                      chromosome;

  int32_t dest_gu_size_before = destination_GU.seq_length();

  // Do not allow translocation if it would increase the size of the receiving
  // GU above a given threshold
  if (dest_gu_size_before + segment_length >
      destination_GU.max_gu_length()) {
    if (exp_m_->output_m()->is_logged(LOG_BARRIER)) {
      // Write an entry in the barrier log file
      fprintf(exp_m_->output_m()->log(LOG_BARRIER),
              "%" PRId64 " %" PRId32 " TRANS %" PRId32 " %" PRId32 " %" PRId32
                  " %" PRId32 "\n",
              AeTime::time(), indiv_->grid_cell()->x() *
                                    indiv_->exp_m()->grid_height()
                                    + indiv_->grid_cell()->y(), segment_length, 0,
              dest_gu_size_before, indiv_->amount_of_dna());
    }
    return false;
  }

  // Provided that OriC must be at position 0
  //
  //   Case 1: inter_GU_ABCDE_to_ACDBE
  //
  //         A    B     C        D     E
  //      |---->=====>----|    |----[>----|
  //          p1r   p2r            p4r
  //                      |
  //                      V
  //        A    C          D     B     E
  //      |---->----|     |----[>====[>----|
  //          p1r             p4r   p4r+(p2r-p1r)
  //
  //
  //   Case 2: inter_GU_ABCDE_to_BDCAE
  //
  //         A    B     C        D     E
  //      |====>----->====|    |----[>----|
  //          p2r   p1r             p4r
  //                      |
  //                      V
  //          B               D     C     A   E
  //       |-----|         |----[>====[>====>----|
  //                           p4r    |
  //                             p4r+(length_-p1r)
  //                                        |
  //                              p4r+(length_-(p1r-p2r))

  // Determine which position comes first and do the corresponding rearrangement
  // int32_t pos_min = Utils::min(pos_1, pos_2);

  if (not invert) {
    if (pos_1_rel < pos_2_rel) {
      segment_length = Utils::mod(pos_2_rel - pos_1_rel, length_);
      inter_GU_ABCDE_to_ACDBE(pos_1_rel, pos_2_rel, pos_4_rel);
    }
    else {
      segment_length = Utils::mod(pos_1_rel - pos_2_rel, length_);
      inter_GU_ABCDE_to_BDCAE(pos_2_rel, pos_1_rel, pos_4_rel);
    }
  }
  else { // invert
    if (pos_1_rel < pos_2_rel) {
      segment_length = Utils::mod(pos_2_rel - pos_1_rel, length_);
      do_inversion(pos_1_rel, pos_2_rel);
      inter_GU_ABCDE_to_ACDBE(pos_1_rel, pos_2_rel, pos_4_rel);
    }
    else { // pos_1_rel > pos_2_rel
      segment_length = Utils::mod(pos_1_rel - pos_2_rel, length_);
      if (pos_2_rel != 0) {
        do_inversion(0, pos_2_rel);
      }
      if (pos_1_rel != length_) {
        do_inversion(pos_1_rel, length_);
      }
      inter_GU_ABCDE_to_BDCAE(pos_2_rel, pos_1_rel, pos_4_rel);
    }
  }

  return true;
}

bool Dna::do_inversion(int32_t pos_1, int32_t pos_2) {
// Invert segment going from pos_1 (included) to pos_2 (excluded)
// Exemple : sequence 011101001100 => 110011010001
  if (pos_1 == pos_2) return false; // Invert everything <=> Invert nothing!
  assert(pos_1 < pos_2);
  //
  //       pos_1         pos_2                   -> 0-
  //         |             |                   -       -
  // 0--------------------------------->      -         -
  //         ===============                  -         - pos_1
  //           tmp (copy)                      -       -
  //                                             -----      |
  //                                             pos_2    <-'
  //

  int32_t seg_length = pos_2 - pos_1;

  // Create the inverted sequence
  char* inverted_segment = NULL;
  inverted_segment = new char[seg_length + 1];

//#pragma simd
//#pragma distribute_point
/*  for (int32_t i = 0, j = pos_2 - 1; i < seg_length; i = i + 4, j = j - 4) {
    if (data_[j] == '0') inverted_segment[i] = '1';
    else inverted_segment[i] = '0';
    if (data_[j - 1] == '0') inverted_segment[i + 1] = '1';
    else inverted_segment[i + 1] = '0';
    if (data_[j - 2] == '0') inverted_segment[i + 2] = '1';
    else inverted_segment[i + 2] = '0';
    if (data_[j - 3] == '0') inverted_segment[i + 3] = '1';
    else inverted_segment[i + 3] = '0';
  }
  inverted_segment[seg_length] = '\0';
*/
#ifdef __SIMD
  #pragma omp simd
#endif
  for (int32_t i = 0, j = pos_2 - 1; i < seg_length; i++, j--) {
    #ifdef BASE_2
    if (data_[j] == '0') inverted_segment[i] = '1';
    else inverted_segment[i] = '0';
    #elif BASE_4
    inverted_segment[i] = get_complementary_base(data_[j]);
    #endif
  }
  inverted_segment[seg_length] = '\0';

  // Remove promoters that included a breakpoint
  gen_unit_->remove_promoters_around(pos_1);
  gen_unit_->remove_promoters_around(pos_2);

  // Invert the sequence
  replace(pos_1, inverted_segment, seg_length);

  // Update promoter list
  if (length_ >= PROM_SIZE) {
    gen_unit_->invert_promoters_included_in(pos_1, pos_2);

    gen_unit_->look_for_new_promoters_around(pos_1);
    gen_unit_->look_for_new_promoters_around(pos_2);

  }

  delete[] inverted_segment;

  return true;
}

bool Dna::do_insertion(int32_t pos, const char* seq_to_insert,
                       int32_t seq_length) {
  // Remove the promoters that will be broken
  gen_unit_->remove_promoters_around(pos);

  // Insert the sequence
  insert(pos, seq_to_insert, seq_length);

  // Look for new promoters
  if (length_ >= PROM_SIZE) {
    gen_unit_->move_all_promoters_after(pos, seq_length);
    gen_unit_->look_for_new_promoters_around(pos, pos + seq_length);
  }

  return true;
}


Mutation* Dna::do_ins_HT(int32_t parent_id) {
  Mutation* mut = nullptr;

  // TODO(dpa) disabled
//  int32_t nb_indivs = exp_m_->pop()->nb_indivs();
//
//  // Insertion transfer
//  // Requirements:
//  //    * A circular exogenote => an alignment on the donor chromosome
//  //    * An alignment between the exogenote and the endogenote
//
//  // 1) Draw a random donor (uniform drawing).
//  // We use the rank because indivs are sorted by rank (1 for the worst,
//  // POP_SIZE for the best).
//  Individual * donor = NULL;
//  do {
//    donor = exp_m_->pop()->
//        indiv_by_rank(exp_m_->sel()->prng()->random(nb_indivs) +
//                              1);
//  }
//  while (donor->id() == parent_id);
//
//  // 2) Look for an alignment within the donor genome
//  VisAVis* alignment_1   = NULL;
//  Dna* donor_dna = donor->genetic_unit(0).dna();
//  int32_t nb_pairs_1 =
//      static_cast<int32_t>(
//          ceil(donor_dna->length() * indiv_->neighbourhood_rate()));
//
//  alignment_1 = donor_dna->search_alignment(donor_dna, nb_pairs_1, DIRECT);
//
//    if (alignment_1 != NULL) {
//      // 3) Make a copy of the sequence to be transferred (the exogenote)
//      GeneticUnit* exogenote =
//          donor_dna->copy_into_new_GU(alignment_1->i_1(),
//                                      alignment_1->i_2());
//
//      // 4) Look for an alignments between the exogenote and the endogenote
//      VisAVis* alignment_2 = NULL;
//      int32_t nb_pairs_2 = static_cast<int32_t>(
//          ceil(length() * indiv_->neighbourhood_rate()));
//
//      alignment_2 =
//          exogenote->dna()->search_alignment(this, nb_pairs_2,
//                                                 BOTH_SENSES);
//
//      if (alignment_2 != NULL) {
//        int32_t gu_length_before = length_;
//        int32_t gu_length_after =
//            gu_length_before + exogenote->dna()->length();
//        int32_t genome_length_before = indiv_->amount_of_dna();
//        int32_t genome_length_after =
//            genome_length_before + exogenote->dna()->length();
//
//        if ((genome_length_after > indiv_->max_genome_length()) ||
//            (gu_length_after > gen_unit_->max_gu_length())) {
//          if (exp_m_->output_m()->is_logged(LOG_BARRIER)) {
//            // Write an entry in the barrier log file
//            fprintf(exp_m_->output_m()->log(LOG_BARRIER),
//                "%" PRId64 " %" PRId32 " INS_TRANSFER %" PRId32 " %" PRId32
//                " %" PRId32 " %" PRId32 "\n",
//                AeTime::time(),
//                indiv_->grid_cell()->x() *                                     indiv_->exp_m()->grid_height()                                     + indiv_->grid_cell()->y(),
//                exogenote->dna()->length(),
//                0,
//                gu_length_before,
//                genome_length_before);
//          }
//        }
//        else {
//          insert_GU(exogenote,
//              alignment_2->i_2(), alignment_2->i_1(),
//              (alignment_2->sense() == INDIRECT));
//          //~ fprintf(logfile, "RESULT:\n%s\n\n\n", new_indiv_dna->data());
//          //~ fflush(logfile);
//
//          // Write a line in transfer logfile
//          if (exp_m_->output_m()->is_logged(LOG_TRANSFER)) {
//            fprintf(exp_m_->output_m()->log(LOG_TRANSFER),
//                    "%" PRId64 " %" PRId32 " %" PRId32 " 0 %" PRId32
//                    " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32
//                    " %" PRId32 " %" PRId16 " %" PRId32 " %" PRId32
//                    " %" PRId16 "\n",
//                    AeTime::time(),
//                    indiv_->grid_cell()->x() *                                     indiv_->exp_m()->grid_height()                                     + indiv_->grid_cell()->y(),
//                    donor->id(),
//                    exogenote->dna()->length(),
//                    0,
//                    genome_length_before,
//                    length(),
//                    alignment_1->i_1(),
//                    alignment_1->i_2(),
//                    alignment_1->score(),
//                    alignment_2->i_1(),
//                    alignment_2->i_2(),
//                    alignment_2->score());
//          }
//
//          #ifdef BIG_DEBUG
//            ae_common::sim->logs()->flush();
//            indiv_->assert_promoters();
//            indiv_->assert_promoters_order();
//          #endif
//
//          if (exp_m_->output_m()->record_tree() &&
//              exp_m_->output_m()->tree_mode() == NORMAL) {
//            char* donor_seq;
//            if (alignment_2->sense() == DIRECT) {
//              donor_seq = exogenote->dna()->
//                  subsequence(alignment_2->i_1(),
//                                  alignment_2->i_1(), LEADING);
//            }
//            else {
//              donor_seq = exogenote->dna()->
//                  subsequence(alignment_2->i_1(),
//                                  alignment_2->i_1(), LAGGING);
//            }
//            // Report the transfer
//            mut = new Mutation();
//            mut->report_ins_HT(alignment_1->i_1(), alignment_1->i_2(),
//                               alignment_2->i_1(), alignment_2->i_2(),
//                               exogenote->dna()->length(),
//                               alignment_1->score(),
//                               alignment_2->score(),
//                               donor->id(),
//                               alignment_2->sense(),
//                               donor_seq);
//            delete [] donor_seq;
//          }
//        }
//
//        delete alignment_2;
//      }
//
//      delete exogenote;
//      delete alignment_1;
//    }

  return mut;
}

Mutation* Dna::do_repl_HT(int32_t parent_id) {
  Mutation* mut = nullptr;

  // TODO(dpa) disabled
//  int32_t nb_indivs = exp_m_->pop()->nb_indivs();
//
//  // Replacement transfer
//  // Requirements:
//  //   * 2 distinct alignments between the (linear) exogenote and the
//  //     endogenote
//
//  // 1) Draw a random donor (uniform drawing).
//  // We use the rank because indivs are sorted by rank (1 for the worst,
//  // POP_SIZE for the best).
//  Individual * donor = NULL;
//  do {
//    donor = exp_m_->pop()->
//        indiv_by_rank(exp_m_->sel()->prng()->random(nb_indivs)
//            + 1);
//  }
//  while (donor->id() == parent_id);
//
//  // 2) Look for an alignment between the parent genome and the donor genome
//  VisAVis* alignment_1   = NULL;
//  VisAVis* alignment_2   = NULL;
//  Dna* donor_dna = donor->genetic_unit(0).dna();
//  AlignmentSense sense = (exp_m_->sel()->prng()->random() < 0.5) ?
//                         DIRECT : INDIRECT;
//  int32_t nb_pairs_1 = static_cast<int32_t>(
//      ceil(length() * indiv_->neighbourhood_rate()));
//  int32_t nb_pairs_2 = static_cast<int32_t>(
//      ceil(length() * indiv_->neighbourhood_rate()));
//  int8_t search_sense = 0;
//
//  alignment_1 = search_alignment(donor_dna, nb_pairs_1, sense);
//  if (alignment_1 != NULL) {
//    if (exp_m_->repl_HT_with_close_points()) {
//      alignment_2 =
//          search_alignment_around_positions(donor_dna,
//                                            alignment_1->i_1(),
//                                            alignment_1->i_2(),
//                                            alignment_1->sense(),
//                                            search_sense);
//      if (alignment_2 != NULL &&
//          alignment_2->i_1() == alignment_1->i_1()) {
//        delete alignment_2;
//        alignment_2 = NULL;
//      }
//      if (alignment_2 != NULL) {
//        // If the second alignment is found upstream of the first alignment,
//        // they are inverse to facilitate
//        if (search_sense == -1) {
//          VisAVis* tmp_alignment = new VisAVis(*alignment_1);
//          alignment_1->copy(alignment_2);
//          alignment_2->copy(tmp_alignment);
//          delete tmp_alignment;
//        }
//      }
//    }
//    else {
//      // Look for a second alignement between the parent and the donor
//      // (must be different from alignment_1)
//      while (alignment_2 == NULL && nb_pairs_2 > 0) {
//        alignment_2 = search_alignment(donor_dna, nb_pairs_2, sense);
//
//        // Forbid the replacement of the whole genome of the parent
//        if (alignment_2 != NULL &&
//            alignment_2->i_1() == alignment_1->i_1()) {
//          delete alignment_2;
//          alignment_2 = NULL;
//        }
//      }
//    }
//
//    // If both alignments were found, proceed to the transfer
//    if (alignment_2 != NULL) {
//      int32_t gu_length_before  = length_;
//      int32_t exogenote_length =
//          Utils::mod(alignment_2->i_2() - alignment_1->i_2() - 1,
//                     donor_dna->length()) + 1;
//      int32_t replaced_seq_length =
//          Utils::mod(alignment_2->i_1() - alignment_1->i_1() - 1,
//                     gu_length_before) + 1;
//      int32_t gu_length_after =
//          gu_length_before - replaced_seq_length + exogenote_length;
//
//      int32_t genome_length_before = indiv_->amount_of_dna();
//      int32_t genome_length_after =
//          genome_length_before - replaced_seq_length + exogenote_length;
//
//      if (genome_length_after < indiv_->min_genome_length() ||
//          genome_length_after > indiv_->max_genome_length() ||
//          gu_length_after < gen_unit_->min_gu_length() ||
//          gu_length_after > gen_unit_->max_gu_length()) {
//        if (exp_m_->output_m()->is_logged(LOG_BARRIER)) {
//          // Write an entry in the barrier log file
//          fprintf(exp_m_->output_m()->log(LOG_BARRIER),
//              "%" PRId64 " %" PRId32 " REPL_TRANSFER %" PRId32 " %" PRId32
//              " %" PRId32 " %" PRId32 "\n",
//              AeTime::time(),
//              indiv_->grid_cell()->x() *                                     indiv_->exp_m()->grid_height()                                     + indiv_->grid_cell()->y(),
//              exogenote_length,
//              replaced_seq_length,
//              gu_length_before,
//              genome_length_before);
//        }
//      }
//      else {
//        // 3) Make a copy of the sequence to be transferred (the exogenote)
//        GeneticUnit* exogenote = NULL;
//        if (sense == DIRECT) {
//          exogenote = donor_dna->copy_into_new_GU(alignment_1->i_2(),
//                                                  alignment_2->i_2());
//        }
//        else {
//          exogenote = donor_dna->copy_into_new_GU(alignment_2->i_2(),
//                                                  alignment_1->i_2());
//        }
//
//        char* alignment1_parent_dna = nullptr;
//        char* alignment2_parent_dna = nullptr;
//        char* alignment1_donor_dna  = nullptr;
//        char* alignment2_donor_dna  = nullptr;
//        if (exp_m_->output_m()->is_logged(LOG_TRANSFER) == true) {
//          if (sense  == DIRECT) {
//            alignment1_parent_dna =
//                subsequence(alignment_1->i_1(),
//                                alignment_1->i_1() +
//                                    2 * indiv_->align_w_zone_h_len() +
//                                    indiv_->align_max_shift(), LEADING);
//            alignment2_parent_dna =
//                subsequence(alignment_2->i_1(),
//                                alignment_2->i_1() +
//                                    2 * indiv_->align_w_zone_h_len() +
//                                    indiv_->align_max_shift(), LEADING);
//            alignment1_donor_dna =
//                donor_dna->subsequence(
//                    alignment_1->i_2(),
//                    alignment_1->i_2() +
//                        2 * indiv_->align_w_zone_h_len() +
//                        indiv_->align_max_shift(), LEADING);
//            alignment2_donor_dna =
//                donor_dna->subsequence(
//                    alignment_2->i_2(),
//                    alignment_2->i_2() +
//                        2 * indiv_->align_w_zone_h_len() +
//                        indiv_->align_max_shift(), LEADING);
//          }
//          else {
//            alignment1_parent_dna =
//                subsequence(alignment_1->i_1(),
//                                alignment_1->i_1() +
//                                    2 * indiv_->align_w_zone_h_len() +
//                                    indiv_->align_max_shift(), LEADING);
//            alignment2_parent_dna =
//                subsequence(alignment_2->i_1(),
//                                alignment_2->i_1() +
//                                    2 * indiv_->align_w_zone_h_len() +
//                                    indiv_->align_max_shift(), LEADING);
//            alignment1_donor_dna =
//                donor_dna->subsequence(
//                    alignment_1->i_2(),
//                    alignment_1->i_2() -
//                        2 * indiv_->align_w_zone_h_len() -
//                        indiv_->align_max_shift(), LAGGING);
//            alignment2_donor_dna =
//                donor_dna->subsequence(
//                    alignment_2->i_2(),
//                    alignment_2->i_2() -
//                        2 * indiv_->align_w_zone_h_len() -
//                        indiv_->align_max_shift(), LAGGING);
//          }
//        }
//
//        // Delete the sequence to be replaced
//        do_deletion(alignment_1->i_1(), alignment_2->i_1());
//        if (alignment_1->i_1() < alignment_2->i_1()) {
//          insert_GU(exogenote, alignment_1->i_1(), 0, sense == INDIRECT);
//        }
//        else {
//          insert_GU(exogenote, 0, 0, sense == INDIRECT);
//        }
//
//        // Write a line in transfer logfile
//        if (exp_m_->output_m()->is_logged(LOG_TRANSFER)) {
//          fprintf(exp_m_->output_m()->log(LOG_TRANSFER),
//              "%" PRId64 " %" PRId32 " %" PRId32 " 1 %" PRId32 " %" PRId32
//              " %" PRId32 " %" PRId32 " %" PRId16 " %" PRId32 " %" PRId32
//              " %" PRId16 " %" PRId32 " %" PRId32 " %" PRId16
//              " %" PRId16 "\n",
//              AeTime::time(),
//              indiv_->grid_cell()->x() *                                     indiv_->exp_m()->grid_height()                                     + indiv_->grid_cell()->y(),
//              donor->id(),
//              exogenote->dna()->length(),
//              replaced_seq_length,
//              genome_length_before,
//              length(),
//              (int16_t) alignment_1->sense(),
//              alignment_1->i_1(),
//              alignment_1->i_2(),
//              alignment_1->score(),
//              alignment_2->i_1(),
//              alignment_2->i_2(),
//              alignment_2->score() ,
//              (int16_t) search_sense);
//
//        fprintf(exp_m_->output_m()->log(LOG_TRANSFER),
//            "\tAlignment 1:\n\t\t%s\n\t\t%s\n"
//            "\tAlignment 2:\n\t\t%s\n\t\t%s\n",
//            alignment1_parent_dna, alignment1_donor_dna,
//            alignment2_parent_dna, alignment2_donor_dna);
//
//            delete [] alignment1_parent_dna;
//            delete [] alignment2_parent_dna;
//            delete [] alignment1_donor_dna;
//            delete [] alignment2_donor_dna;
//        }
//
//        if (exp_m_->output_m()->record_tree() &&
//            exp_m_->output_m()->tree_mode() == NORMAL) {
//          // Report the transfer
//          char* donor_seq;
//          if (alignment_2->sense() == DIRECT) {
//            donor_seq = exogenote->dna()->
//                subsequence(0, exogenote->dna()->length(), LEADING);
//          }
//          else {
//            donor_seq = exogenote->dna()->
//                subsequence(0, exogenote->dna()->length(), LAGGING);
//          }
//          mut = new Mutation();
//          mut->report_repl_HT(alignment_1->i_1(), alignment_1->i_2(),
//                              alignment_2->i_1(), alignment_2->i_2(),
//                              replaced_seq_length,
//                              exogenote->dna()->length(),
//                              alignment_1->score(),
//                              alignment_2->score(),
//                              donor->id(),
//                              alignment_2->sense(),
//                              donor_seq);
//          delete [] donor_seq;
//        }
//
//        delete exogenote;
//      }
//      delete alignment_2;
//    }
//    delete alignment_1;
//  }
  return mut;
}

bool Dna::do_ins_HT(int32_t pos, const char* seq_to_insert,
                    int32_t seq_length) { // NOLINT(whitespace/braces)
  // Remove the promoters that will be broken
  gen_unit_->remove_promoters_around(pos);

  // Insert the sequence
  insert(pos, seq_to_insert, seq_length);

  // Look for new promoters
  if (length_ >= PROM_SIZE) {
    gen_unit_->move_all_promoters_after(pos, seq_length);
    gen_unit_->look_for_new_promoters_around(pos, pos + seq_length);
  }

  return true;
}

bool Dna::do_repl_HT(int32_t pos1, int32_t pos2, const char* seq_to_insert,
                     int32_t seq_length) {
  // Remove the promoters that will be broken
  gen_unit_->remove_promoters_around(pos1);

  // Delete the replaced segment
  do_deletion(pos1, pos2);

  // Insert the sequence
  int32_t insertion_position;
  if (pos1 < pos2) {
    insertion_position = pos1;
  }
  else {
    insertion_position = 0;
  }
  insert(insertion_position, seq_to_insert, seq_length);

  // Look for new promoters
  if (length_ >= PROM_SIZE) {
    gen_unit_->move_all_promoters_after(insertion_position, seq_length);
    gen_unit_->look_for_new_promoters_around(insertion_position,
                                             insertion_position + seq_length);
  }

  return true;
}

std::tuple<int32_t, int32_t, int32_t, int32_t, bool> Dna::undergo_this_mutation(const Mutation& mut) {
  //******************************************************
  
  //    printf("UNDERGO MUT : \n");
  switch (mut.mut_type()) {
  case SWITCH :{
    //        printf("%d -- DO Switch %d\n",indiv()->id(),dynamic_cast<const PointMutation&>(mut).pos());
    //**************************************************
    const auto& s_pos = dynamic_cast<const PointMutation&>(mut).pos();
    //**************************************************
    do_switch(s_pos);
    return std::make_tuple(s_pos, 0, 0, 0, false);
    break;
  }
  case S_INS : {
    const auto& s_ins = dynamic_cast<const SmallInsertion&>(mut);
    //      printf("%d -- DO INSERT ",indiv()->id());
    //      fflush(stdout);
    //
    //             printf("%d ",s_ins.pos());
    //      fflush(stdout);
    //             printf("%d ",s_ins.length());
    //      fflush(stdout);
    //      for (int i = 0; i < s_ins.length(); i++)
    //          printf("%c",s_ins.seq()[i]);
    //      printf("\n");
    //
    //      printf("%s\n", s_ins.seq());
    //      fflush(stdout);
    
    do_small_insertion(s_ins.pos(), s_ins.length(), s_ins.seq());
    return std::make_tuple(s_ins.pos(), s_ins.length(), 0, 0, false);
    break;
  }
  case S_DEL : {
    const auto& s_del = dynamic_cast<const SmallDeletion&>(mut);
    //        printf("%d -- DO DELETE %d %d\n",indiv()->id(),s_del.pos(), s_del.length());
    do_small_deletion(s_del.pos(), s_del.length());
    return std::make_tuple(s_del.pos(), s_del.length(), 0, 0, false);
    break;
  }
  case DUPL : {
    const auto& dupl = dynamic_cast<const Duplication&>(mut);
    //      printf("%d -- DO DUPLICATION %d %d %d\n",indiv()->id(),dupl.pos1(), dupl.pos2(), dupl.pos3());
    do_duplication(dupl.pos1(), dupl.pos2(), dupl.pos3());
    //std::cout << "pos1 and pos2:" <<dupl.pos1()<<"&" <<dupl.pos2()<<"pos3:" << dupl.pos3();
    return std::make_tuple(dupl.pos1(), dupl.pos2(), dupl.pos3(), 0, false);
    break;
  }
  case DEL : {
    const auto& del = dynamic_cast<const Deletion&>(mut);
    //      printf("%d -- DO DELETE LARGE %d %d\n",indiv()->id(),del.pos1(), del.pos2());
    do_deletion(del.pos1(), del.pos2());
    return std::make_tuple(del.pos1(), del.pos2(), 0, 0, false);
    break;
  }
  case TRANS : {
    const auto& trans = dynamic_cast<const Translocation&>(mut);
    //      printf("%d -- DO TRANS %d %d %d %d %d\n",indiv()->id(),trans.pos1(), trans.pos2(), trans.pos3(), trans.pos4(),
    //             trans.invert());
    do_translocation(trans.pos1(), trans.pos2(), trans.pos3(), trans.pos4(),
                     trans.invert());
    return std::make_tuple(trans.pos1(), trans.pos2(), trans.pos3(), trans.pos4(),
                           trans.invert());
    break;
  }
  case INV : {
    const auto& inv = dynamic_cast<const Inversion&>(mut);
    //      printf("%d -- DO INVERSION %d %d\n",indiv()->id(),inv.pos1(), inv.pos2());
    do_inversion(inv.pos1(), inv.pos2());
    return std::make_tuple(inv.pos1(), inv.pos2(), 0, 0, false);
    break;
  }
  case INS_HT : {
    const auto& ins_ht = dynamic_cast<const InsertionHT&>(mut);
    //      printf("NON NON \n");
    do_ins_HT(ins_ht.receiver_pos(), ins_ht.seq(), ins_ht.length());
    return std::make_tuple(0, 0, 0, 0, false);
    break;
  }
  case REPL_HT : {
    const auto& repl_ht = dynamic_cast<const ReplacementHT&>(mut);
    //        printf("NON NON NOPEP \n");
    do_repl_HT(repl_ht.receiver_pos1(), repl_ht.receiver_pos2(),
               repl_ht.seq(), repl_ht.length());
    return std::make_tuple(0, 0, 0, 0, false);
    break;
  }
  default :
    fprintf(stderr, "ERROR, invalid mutation type in file %s:%d\n",
            __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    break;
  }
}

void Dna::compute_statistical_data() {
// ****************************************************************************
//                            WARNING     Deprecated
// ****************************************************************************
  assert(false);
}

#ifndef __REGUL

void Dna::set_GU(std::vector<std::list < Rna>>

rna_list,
const GeneticUnit* GU
) {
for (

auto& strand : { LEADING, LAGGING }

)
for (
auto& rna :
rna_list[strand])
rna.
set_genetic_unit(GU);
}
#else
void Dna::set_GU(std::vector<std::list<Rna_R>> rna_list, const GeneticUnit* GU) {
  for (auto& strand: {LEADING, LAGGING})
    for (auto& rna: rna_list[strand])
      rna.set_genetic_unit(GU);
}
#endif

GeneticUnit* Dna::extract_into_new_GU(int32_t pos_1, int32_t pos_2) {
  assert(pos_1 < pos_2);
  int32_t seq_length = pos_2 - pos_1;

  // =============== Remove/Extract promoters from old sequence ===============
  // Remove promoters around breakpoints
  gen_unit_->remove_promoters_around(pos_1);
  gen_unit_->remove_promoters_around(pos_2);

  // Remove promoters belonging to the sequence (to be extracted) from the
  // "old" GU and put them in a stand-alone promoter list (with indices
  // ranging from 0 to seq_length-1)
  Promoters2Strands proms_GU_1 = {{},
                                  {}};
  gen_unit_->extract_promoters_included_in(pos_1, pos_2, proms_GU_1);
  GeneticUnit::shift_promoters(proms_GU_1, -pos_1, length_);

  // ==================== Manage sequences ====================
  // Copy the sequence in a stand-alone char* (size must be multiple of
  // BLOCK_SIZE)
  int32_t length_GU_1 = seq_length;
  char* sequence_GU_1 = new char[nb_blocks(length_GU_1) * BLOCK_SIZE];
  memcpy(sequence_GU_1, &data_[pos_1], length_GU_1 * sizeof(char));
  sequence_GU_1[length_GU_1] = '\0';

  // Remove the sequence from the "old" GU
  int32_t length_GU_0 = length_ - seq_length;
  char* sequence_GU_0 = new char[nb_blocks(length_GU_0) * BLOCK_SIZE];
  memcpy(sequence_GU_0, data_, pos_1 * sizeof(char));
  memcpy(&sequence_GU_0[pos_1], &data_[pos_2],
         (length_ - pos_2) * sizeof(char));
  sequence_GU_0[length_GU_0] = '\0';

  set_data(sequence_GU_0, length_GU_0);

  // ==================== Create the new genetic unit ====================
  GeneticUnit* GU_1 =
      new GeneticUnit(indiv_, sequence_GU_1, length_GU_1, proms_GU_1);

  // ==================== Update promoter lists ====================
  // Shift the position of the promoters of the "old" GU
  gen_unit_->move_all_promoters_after(pos_1, -seq_length);

  // Look for new promoters around breakpoints
  gen_unit_->look_for_new_promoters_around(pos_1);
  GU_1->look_for_new_promoters_around(0);

  return GU_1;
}


/*!
  \brief Copy the sequence going from pos_1 (included) to pos_2 (excluded)
  into a new standalone genetic unit.

  The new genetic unit's list of promoter is up-to-date.
  if (pos_1 == pos_2), the whole genome is copied
*/
GeneticUnit* Dna::copy_into_new_GU(int32_t pos_1, int32_t pos_2) const {
  int32_t seq_length = Utils::mod(pos_2 - pos_1, length_);
  if (seq_length == 0) seq_length = length_;

  // ==================== Copy promoters from old sequence ====================
  // Copy the promoters belonging to the sequence to be copied from the "old"
  // GU into a stand-alone promoter list (with indices ranging from 0 to
  // seq_length - 1)
  Promoters2Strands proms_new_GU = {{},
                                    {}};
  gen_unit_->copy_promoters_included_in(pos_1, pos_2, proms_new_GU);
  GeneticUnit::shift_promoters(proms_new_GU, -pos_1, length_);


  // ==================== Manage sequences ====================
  // Copy the sequence in a stand-alone char* (size must be multiple of
  // BLOCK_SIZE)
  int32_t length_new_GU = seq_length;
  char* sequence_new_GU = new char[nb_blocks(length_new_GU) * BLOCK_SIZE];
  if (pos_1 < pos_2) {
    memcpy(sequence_new_GU, &data_[pos_1], length_new_GU * sizeof(char));
  }
  else {
    memcpy(sequence_new_GU, &data_[pos_1], (length_ - pos_1) * sizeof(char));
    memcpy(&sequence_new_GU[length_ - pos_1], &data_[0], pos_2 * sizeof(char));
  }
  sequence_new_GU[length_new_GU] = '\0';


  // ==================== Create the new genetic unit ====================
  GeneticUnit* new_GU =
      new GeneticUnit(indiv_, sequence_new_GU, length_new_GU, proms_new_GU);

  // ==================== Update new GU promoter list ====================
  // Look for new promoters around breakpoints
  new_GU->look_for_new_promoters_around(0);

  return new_GU;
}

/**
 * \brief Insert the genetic unit GU_to_insert at pos_B, through pos_D.
 *
 * Sequence is inverted if invert == true
 *
 *
 * GU to insert: segments C and D, breakpoint pos_D.
 * Current GU:   segments A and B, breakpoint pos_B.

  \verbatim
  If invert is false, the insertion will render ADCB
           A        B                  C        D
       |-------[>-------|     +     |=====[>========|
             pos_B                      pos_D
                              |
                              V
                   A       D       C        B
               |-------[>=====|========[>-------|

  If invert is true, the insertion will render ACpDpB
  with Cp (resp. Dp) = inverted C (resp. D).

           A        B                  C        D
       |-------[>-------|     +     |=====<]========|
             pos_B                      pos_D
                              |
                              V
                   A       Cp     Dp        B
               |-------[>=====|========[>-------|

  \endverbatim

  Sequence from GU_to_insert is untouched but its list of promoters is emptied
*/
void Dna::insert_GU(GeneticUnit* GU_to_insert, int32_t pos_B, int32_t pos_D,
                    bool invert) {
  // Compute segment lengths
  const char* GUti_data = GU_to_insert->dna()->data();
  int32_t len_A = pos_B;
  int32_t len_B = length_ - pos_B;
  int32_t len_C = pos_D;
  int32_t len_D = GU_to_insert->dna()->length() - pos_D;
  int32_t len_AB = length_;
  int32_t len_CD = GU_to_insert->dna()->length();
  int32_t len_ABCD = len_AB + len_CD;


  // ==================== Insert the sequence ====================
  // Create new genome
  char* new_seq = new char[BLOCK_SIZE * nb_blocks(len_ABCD)];

  // Insert A
  memcpy(new_seq, data_, len_A * sizeof(char));

  // Insert C and D (inverted?)
  if (invert) {
    // Build Cp and Dp
    char* seq_Cp = new char[pos_D + 1];
    for (int32_t i = 0, j = pos_D - 1; i < len_C; i++, j--) {
      #ifdef BASE_2
      if (GUti_data[j] == '0') seq_Cp[i] = '1';
      else seq_Cp[i] = '0';
      #elif BASE_4
      seq_Cp[i] = get_complementary_base(GUti_data[j]);
      #endif
    }
    seq_Cp[len_C] = '\0';

    char* seq_Dp = new char[len_D + 1];
    for (int32_t i = 0, j = len_CD - 1; i < len_D; i++, j--) {
      #ifdef BASE_2
      if (GUti_data[j] == '0') seq_Dp[i] = '1';
      else seq_Dp[i] = '0';
      #elif BASE_4
      seq_Dp[i] = get_complementary_base(GUti_data[j]);
      #endif
    }
    seq_Dp[len_D] = '\0';

    // Insert Cp and DP
    // TODO(???) Maybe we should construct Cp and Dp directly at their rightful
    // place...
    memcpy(&new_seq[len_A], seq_Cp, len_C * sizeof(char));
    memcpy(&new_seq[len_A + len_C], seq_Dp, len_D * sizeof(char));

    delete[] seq_Cp;
    delete[] seq_Dp;
  }
  else { // if (invert == false)
    // Insert D and C
    memcpy(&new_seq[len_A], &GUti_data[pos_D], len_D * sizeof(char));
    memcpy(&new_seq[len_A + len_D], GUti_data, len_C * sizeof(char));
  }

  // Insert B
  memcpy(&new_seq[len_A + len_C + len_D], &data_[pos_B], len_B * sizeof(char));
  new_seq[len_ABCD] = '\0';


  // Remove promoters that are astride segments A and B : breakpoint pos_B
  gen_unit_->remove_promoters_around(pos_B);


  set_data(new_seq, len_ABCD);




  // ==================== Manage promoters ====================
  // Remove promoters that are astride segments C and D : breakpoint pos_D
  GU_to_insert->remove_promoters_around(pos_D);

  // Shift the position of the promoters of segment B
  gen_unit_->move_all_promoters_after(pos_B, len_CD);

  // Extract promoters of segments C and D.
  // NOTE : We want ALL THE PROMOTERS to be transfered, not only those that
  //        are completely included in segment C or D. Hence, we will use
  //        extract_promoters_starting_between() instead of
  //        extract_promoters_included_in().
  // NOTE : Once removed the promoters starting on sequence D, the remaining
  //        is precisely the promoters starting on sequence C (and they are
  //        at their rightful position). We can hence directly use the list
  //        of promoters from GU_to_insert.
  Promoters2Strands proms_C = {{},
                               {}};
  Promoters2Strands proms_D = {{},
                               {}};

  if (pos_D != 0) {
    // TODO(???) : Manage this in the different functions? with a parameter
    // WholeGenomeEventHandling ?
    GU_to_insert->extract_promoters_starting_between(0, pos_D, proms_C);
  }
  GU_to_insert->extract_promoters_starting_between(pos_D, len_CD, proms_D);
  assert(GU_to_insert->rna_list()[LEADING].empty());
  assert(GU_to_insert->rna_list()[LAGGING].empty());
  GeneticUnit::shift_promoters(proms_D, -len_C, len_D);

  if (invert) {
    //~ printf("+++++++++++++++++++++++++++++++++++++\n");
    //~ GeneticUnit::print_rnas(proms_C);
    //~ GeneticUnit::print_rnas(proms_D);
    //~ printf("//////////////////////////////////////\n");

    GeneticUnit::invert_promoters(proms_C, len_C);
    GeneticUnit::invert_promoters(proms_D, len_D);

    //~ GeneticUnit::print_rnas(proms_C);
    //~ GeneticUnit::print_rnas(proms_D);
    //~ printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");

    gen_unit_->insert_promoters_at(proms_C, len_A);
    gen_unit_->insert_promoters_at(proms_D, len_A + len_C);
  }
  else { // if (invert == false)
    gen_unit_->insert_promoters_at(proms_D, len_A);
    gen_unit_->insert_promoters_at(proms_C, len_A + len_D);
  }

  // Look for new promoters around breakpoints
  gen_unit_->look_for_new_promoters_around(pos_B);
  gen_unit_->look_for_new_promoters_around(pos_B + len_CD);

  gen_unit_->take_ownership_of_all_rnas();
}


/**
 * \brief Looks for an alignment between this and chrom2 in the given sense
 * with max nb_pairs trials. nb_pairs is updated accordingly
 *
 * Performs local alignment searches between this and chrom2 around randomly
 * drawn pairs of points.
 * The minimum score will be generated according to align_fun_shape and
 * associated parameters for each pair of points.
 * The parameter nb_pairs will be updated according to how many trials
 * were necessary for an alignment to be found.
 *
 * The sense of the searched alignment can be either DIRECT, INDIRECT or
 * BOTH_SENSE.
 * In the latter case, the sense will be randomly drawn (uniformly between
 * DIRECT and INDIRECT) for each pair of points.
 */
VisAVis* Dna::search_alignment(Dna* chrom2, int32_t& nb_pairs,
                               AlignmentSense sense) {
  VisAVis* alignment = NULL;
  // Which sense (direct or indirect)
  AlignmentSense cur_sense = sense;
  // Minimum alignement score needed to recombine (stochastic)
  int16_t needed_score;

  for (; nb_pairs > 0; nb_pairs--) {
    ///////////////////////////////////////
    //  1) Draw random sense (if needed) //
    ///////////////////////////////////////
    if (sense == BOTH_SENSES) {
      cur_sense = (indiv_->mut_prng_->random() < 0.5) ? DIRECT : INDIRECT;
    }

    /////////////////////////////////////////////////////
    // 2) Determine the minimum alignment score needed //
    /////////////////////////////////////////////////////
    if (indiv_->align_fun_shape() == LINEAR) {
      needed_score = static_cast<int16_t>(
          ceil(indiv_->align_lin_min() +
               indiv_->mut_prng_->random() *
               (indiv_->align_lin_max() -
                indiv_->align_lin_min())));
    }
    else {
      // I want the probability of rearrangement for an alignment of score
      // <score> to be prob = 1 / (1 + exp(-(score-mean_score)/lambda))
      // The score needed for a rearrangement to take place with a given
      // random drawing is hence
      // needed_score = ceil(-lambda * log(1/rand - 1) + mean)
      needed_score = static_cast<int16_t>(
          ceil(-indiv_->align_sigm_lambda() *
               log(1 / indiv_->mut_prng_->random() - 1) +
               indiv_->align_sigm_mean()));
      if (needed_score < 0) needed_score = 0;
    }

    ///////////////////////////////////////////////////////////////
    // 3) Determine where to look for an alignement (draw seeds) //
    ///////////////////////////////////////////////////////////////
    int32_t seed1 = indiv_->mut_prng_->random(length_);
    int32_t seed2 = indiv_->mut_prng_->random(chrom2->length());

    ////////////////////////////////////////////////////////////////////
    // 3) Test the existence of an alignment with a high enough score //
    ////////////////////////////////////////////////////////////////////
    if (cur_sense == DIRECT) {
      alignment = Alignment::search_alignment_direct(this, seed1,
                                                     chrom2, seed2,
                                                     needed_score);
      if (alignment != NULL) {
        return alignment;
      }
    }
    else { // if (cur_sense = INDIRECT)
      alignment = Alignment::search_alignment_indirect(this, seed1,
                                                       chrom2, seed2,
                                                       needed_score);
      if (alignment != NULL) {
        return alignment;
      }
    }
  }

  return NULL;
}

/**
 * \brief Looks for an alignment between this and chrom2 in the given sense
 *  around the given positions
 *
 *  Performs local alignment searches between this and chrom2 around the given
 *  positions
 *  The minimum score will be generated according to align_fun_shape and
 *  associated parameters for each pair of points.
 *  The parameter nb_pairs will be updated according to how many trials were
 *  necessary for an alignment to be found.
 *
 *  The sense of the searched alignment can be either DIRECT, INDIRECT or
 *  BOTH_SENSE.
 *  In the latter case, the sense will be randomly drawn (uniformly between
 *  DIRECT and INDIRECT) for each pair of points.
 */
VisAVis* Dna::search_alignment_around_positions(Dna* chrom2,
                                                int32_t chrom1_pos_1,
                                                int32_t chrom2_pos_1,
                                                AlignmentSense sense,
                                                int8_t& search_sense) {
  VisAVis* alignment = NULL;
  VisAVis* tmp_alignment = NULL;
  // Which sense (direct or indirect)
  AlignmentSense cur_sense = sense;
  // Minimum alignement score needed to recombine (stochastic)
  int16_t needed_score;
  int32_t chrom1_pos_for_research;
  int32_t chrom2_pos_for_research;
  int32_t size_between_two_alignments = 3 * indiv_->align_w_zone_h_len();

  ///////////////////////////////////////
  //  1) Draw random sense (if needed) //
  ///////////////////////////////////////
  if (sense == BOTH_SENSES) {
    printf("WARNING : Alignment could not be searched in both senses "
               "in %s:%d\n", __FILE__, __LINE__);
    return (NULL);
  }

  /////////////////////////////////////////////////////
  // 2) Determine the minimum alignment score needed //
  /////////////////////////////////////////////////////
  if (indiv_->align_fun_shape() == LINEAR) {
    needed_score = static_cast<int16_t>(
        ceil(indiv_->align_lin_min() +
             indiv_->mut_prng_->random() *
             (indiv_->align_lin_max() -
              indiv_->align_lin_min())));
  }
  else {
    // I want the probability of rearrangement for an alignment of score
    // <score> to be prob = 1 / (1 + exp(-(score-mean_score)/lambda))
    // The score needed for a rearrangement to take place with a given
    // random drawing is hence
    // needed_score = ceil(-lambda * log(1/rand - 1) + mean)
    needed_score = static_cast<int16_t>(
        ceil(-indiv_->align_sigm_lambda() *
             log(1 / indiv_->mut_prng_->random() - 1) +
             indiv_->align_sigm_mean()));
    if (needed_score < 0) needed_score = 0;
  }

  /////////////////////////////////////////////////////////
  // 3) Determine the sense by which the research begins //
  /////////////////////////////////////////////////////////
  int16_t first_research_sense = (indiv_->mut_prng_->random() < 0.5) ? 1 : -1;
  int16_t second_research_sense = -1 * first_research_sense;

  /////////////////////////////////////////////////////////////////////////////
  // 4) Test the first sense for the existence of an alignment with a high
  //    enough score
  /////////////////////////////////////////////////////////////////////////////
  chrom1_pos_for_research = chrom1_pos_1;
  chrom2_pos_for_research = chrom2_pos_1;
  int16_t i = 0;

  while (indiv_->mut_prng_->random() < 1 - exp_m_->repl_HT_detach_rate()) {
    chrom1_pos_for_research =
        Utils::mod(chrom1_pos_for_research +
                   first_research_sense * size_between_two_alignments,
                   this->length());
    if (cur_sense == DIRECT) {
      chrom2_pos_for_research =
          Utils::mod(chrom2_pos_for_research +
                     first_research_sense * size_between_two_alignments,
                     chrom2->length());
      tmp_alignment =
          Alignment::search_alignment_direct(this, chrom1_pos_for_research,
                                             chrom2, chrom2_pos_for_research,
                                             needed_score);
    }
    else { // if (cur_sense = INDIRECT)
      chrom2_pos_for_research =
          Utils::mod(chrom2_pos_for_research -
                     first_research_sense * size_between_two_alignments,
                     chrom2->length());
      tmp_alignment =
          Alignment::search_alignment_indirect(this, chrom1_pos_for_research,
                                               chrom2, chrom2_pos_for_research,
                                               needed_score);
    }

    if (tmp_alignment == NULL) {
      if (alignment != NULL) {
        search_sense = first_research_sense;
        return alignment;
      }
      else {
        break;
      }
    }
    else {
      if (alignment != NULL) {
        alignment->copy(tmp_alignment);
      }
      else {
        alignment = new VisAVis(*tmp_alignment);
      }
      delete tmp_alignment;
      chrom1_pos_for_research = alignment->i_1();
      chrom2_pos_for_research = alignment->i_2();
    }
    i++;
  }

  if (alignment != NULL) {
    search_sense = first_research_sense;
    return alignment;
  }

  /////////////////////////////////////////////////////////////////////////////
  // 5) Test the second sense for the existence of an alignment with a high
  //    enough score
  /////////////////////////////////////////////////////////////////////////////
  alignment = NULL;
  chrom1_pos_for_research = chrom1_pos_1;
  chrom2_pos_for_research = chrom2_pos_1;
  i = 0;
  while (indiv_->mut_prng_->random() < 1 - exp_m_->repl_HT_detach_rate()) {
    chrom1_pos_for_research =
        Utils::mod(chrom1_pos_for_research +
                   second_research_sense * size_between_two_alignments,
                   this->length());
    if (cur_sense == DIRECT) {
      chrom2_pos_for_research =
          Utils::mod(chrom2_pos_for_research +
                     second_research_sense * size_between_two_alignments,
                     chrom2->length());
      tmp_alignment =
          Alignment::search_alignment_direct(this, chrom1_pos_for_research,
                                             chrom2, chrom2_pos_for_research,
                                             needed_score);
    }
    else { // if (cur_sense = INDIRECT)
      chrom2_pos_for_research =
          Utils::mod(chrom2_pos_for_research -
                     second_research_sense * size_between_two_alignments,
                     chrom2->length());
      tmp_alignment =
          Alignment::search_alignment_indirect(this, chrom1_pos_for_research,
                                               chrom2, chrom2_pos_for_research,
                                               needed_score);
    }

    if (tmp_alignment == NULL) {
      if (alignment != NULL) {
        search_sense = second_research_sense;
        return alignment;
      }
      else {
        break;
      }
    }
    else {
      if (alignment != NULL) {
        alignment->copy(tmp_alignment);
      }
      else {
        alignment = new VisAVis(*tmp_alignment);
      }
      delete tmp_alignment;
      chrom1_pos_for_research = alignment->i_1();
      chrom2_pos_for_research = alignment->i_2();
    }
    i++;
  }
  if (alignment != NULL) {
    search_sense = second_research_sense;
    return alignment;
  }

  return NULL;
}

// =================================================================
//                           Protected Methods
// =================================================================
/**
 * Extract the sequence going from pos_1 (included) to pos_2 (excluded)
 * into a new standalone genetic unit.
 * Promoter lists are created / updated accordingly
 */
void Dna::ABCDE_to_ADCBE(int32_t pos_B, int32_t pos_C, int32_t pos_D,
                         int32_t pos_E) {
  // Rearrange the sequence from ABCDE to ADCBE (complex translocation
  // of segment defined between positions pos_B and pos_D)
  //
  // Segments are identified by pos_x values as shown below.
  //
  // TODO(dpa) CHECK THIS !!!
  // WARNING : Segment C includes nucleotide at pos_D // NOTE : WTF???
  //
  //         A      B        C       D       E
  //      |----->=======[>=======>-------[>-------|        =>
  //          pos_B   pos_C    pos_D   pos_E
  //
  //                         |
  //                         V
  //
  //         A      D        C       B        E
  //      |----->-------[>=======>=======[>-------|
  // Check points' order and range
  assert(pos_B >= 0 && pos_C >= pos_B && pos_D >= pos_C && pos_E >= pos_D &&
         pos_E <= length_);

  // Compute segment lengths
  int32_t len_A = pos_B;
  int32_t len_B = pos_C - pos_B;
  int32_t len_C = pos_D - pos_C;
  int32_t len_D = pos_E - pos_D;
  int32_t len_E = length_ - pos_E;
  int32_t len_AD = len_A + len_D;
  int32_t len_ADC = len_AD + len_C;
  int32_t len_ADCB = len_ADC + len_B;

  // Create new sequence
  char* new_genome = new char[nb_blocks_ * BLOCK_SIZE];

  memcpy(new_genome, data_, len_A * sizeof(char));
  memcpy(&new_genome[len_A], &data_[pos_D], len_D * sizeof(char));
  memcpy(&new_genome[len_AD], &data_[pos_C], len_C * sizeof(char));
  memcpy(&new_genome[len_ADC], &data_[pos_B], len_B * sizeof(char));
  memcpy(&new_genome[len_ADCB], &data_[pos_E], len_E * sizeof(char));
  new_genome[length_] = '\0';

  // Replace sequence
  // NB : The size of the genome doesn't change. Therefore, we don't nee
  // to update length_ and nb_blocks_
  delete[] data_;
  data_ = new_genome;


  // ========== Update promoter list ==========
  if (length_ >= PROM_SIZE) {

    // Remove promoters that include a breakpoint
    gen_unit_->remove_promoters_around(pos_B);
    gen_unit_->remove_promoters_around(pos_C);
    gen_unit_->remove_promoters_around(pos_D);
    gen_unit_->remove_promoters_around(pos_E);

    // Create temporary lists for promoters to move and/or invert
    Promoters2Strands promoters_B = {{},
                                     {}};
    Promoters2Strands promoters_C = {{},
                                     {}};
    Promoters2Strands promoters_D = {{},
                                     {}};

    // Extract promoters that are totally included in each segment to be moved
    // and shift them to their new positions
    if (len_B >= PROM_SIZE) {
      gen_unit_->extract_promoters_included_in(pos_B, pos_C, promoters_B);

      GeneticUnit::shift_promoters(promoters_B, len_D + len_C,
                                   gen_unit_->dna()->length());
    }
    if (len_C >= PROM_SIZE) {
      gen_unit_->extract_promoters_included_in(pos_C, pos_D, promoters_C);

      GeneticUnit::shift_promoters(promoters_C, len_D - len_B,
                                   gen_unit_->dna()->length());
    }
    if (len_D >= PROM_SIZE) {
      gen_unit_->extract_promoters_included_in(pos_D, pos_E, promoters_D);

      GeneticUnit::shift_promoters(promoters_D, -len_B - len_C,
                                   gen_unit_->dna()->length());
    }

    // Reinsert the shifted promoters
    gen_unit_->insert_promoters(promoters_B);
    gen_unit_->insert_promoters(promoters_C);
    gen_unit_->insert_promoters(promoters_D);

    // 5) Look for new promoters including a breakpoint
    gen_unit_->look_for_new_promoters_around(len_A);
    gen_unit_->look_for_new_promoters_around(len_AD);
    gen_unit_->look_for_new_promoters_around(len_ADC);
    gen_unit_->look_for_new_promoters_around(len_ADCB);


      // printf("%d -- CPU_TRANSLOC\n",indiv_->id_);
      // for (auto& strand: {LEADING, LAGGING}) {
      //   for (auto& rna: gen_unit_->rna_list_[strand]) {
      //     printf("%d -- Promoter list %d\n",indiv_->id_,rna.promoter_pos());
      //   }
      // }
  }
}

void Dna::ABCDE_to_ADBpCpE(int32_t pos_B, int32_t pos_C, int32_t pos_D,
                           int32_t pos_E) {
  // Rearrange the sequence from ABCDE to ADBpCpE (complex translocation
  // with inversion of segment defined between positions pos_B and pos_D)
  // Bp (resp Cp) stands for inverted B (resp C)
  //
  // Segments are identified by pos_x values as shown below.
  //
  // TODO(dpa) CHECK THIS !!!
  // WARNING : Segment C includes nucleotide at pos_D // NOTE : WTF???
  //
  //         A      B        C       D        E
  //      |----->=======[>=======>-------<]-------|
  //          pos_B   pos_C    pos_D   pos_E
  //
  //                         |
  //                         V
  //
  //         A      D        Bp      Cp       E
  //      |----->-------<]=======<=======<]-------|


  // Check points' order and range
  assert(pos_B >= 0 && pos_C >= pos_B && pos_D >= pos_C && pos_E >= pos_D &&
                                                           pos_E <= length_);

  // Compute segment lengths
  int32_t len_A = pos_B;
  int32_t len_B = pos_C - pos_B;
  int32_t len_C = pos_D - pos_C;
  int32_t len_D = pos_E - pos_D;
  int32_t len_E = length_ - pos_E;
  int32_t len_AD = len_A + len_D;
  int32_t len_ADB = len_AD + len_B;
  int32_t len_ADBC = len_ADB + len_C;

  // Create new sequence
  char* new_genome = new char[nb_blocks_ * BLOCK_SIZE];

  // Copy segments A and D
  memcpy(new_genome, data_, len_A * sizeof(char));
  memcpy(&new_genome[len_A], &data_[pos_D], len_D * sizeof(char));


  // Build Bp and put it in the new genome
  char* inverted_segment = new char[len_B + 1];

//#pragma simd
//#pragma distribute_point
 #ifdef __SIMD
#pragma omp simd
 #endif
  for (int32_t i = 0, j = pos_C - 1; i < len_B; i++, j--) {
    #ifdef BASE_2
    if (data_[j] == '0') inverted_segment[i] = '1';
    else inverted_segment[i] = '0';
    #elif BASE_4
    inverted_segment[i] = get_complementary_base(data_[j]);
    #endif
  }
  inverted_segment[len_B] = '\0';

  memcpy(&new_genome[len_AD], inverted_segment, len_B * sizeof(char));

  delete[] inverted_segment;


  // Build Cp and put it in the new genome
  inverted_segment = new char[len_C + 1];

//#pragma simd
//#pragma distribute_point
#ifdef __SIMD
#pragma omp simd
#endif
  for (int32_t i = 0, j = pos_D - 1; i < len_C; i++, j--) {
    #ifdef BASE_2
    if (data_[j] == '0') inverted_segment[i] = '1';
    else inverted_segment[i] = '0';
    #elif BASE_4
    inverted_segment[i] = get_complementary_base(data_[j]);
    #endif
  }
  inverted_segment[len_C] = '\0';

  memcpy(&new_genome[len_ADB], inverted_segment, len_C * sizeof(char));

  delete[] inverted_segment;

  // Copy segment E into the new genome
  memcpy(&new_genome[len_ADBC], &data_[pos_E], len_E * sizeof(char));
  new_genome[length_] = '\0';


  // Replace sequence
  delete[] data_;
  data_ = new_genome;


  // ========== Update promoter list ==========
  if (length_ >= PROM_SIZE) {
    // Remove promoters that include a breakpoint
    gen_unit_->remove_promoters_around(pos_B);
    gen_unit_->remove_promoters_around(pos_C);
    gen_unit_->remove_promoters_around(pos_D);
    gen_unit_->remove_promoters_around(pos_E);

    // Create temporary lists for promoters to move and/or invert
    Promoters2Strands promoters_B = {{},
                                     {}};
    Promoters2Strands promoters_C = {{},
                                     {}};
    Promoters2Strands promoters_D = {{},
                                     {}};

    // 2) Extract promoters that are totally included in each segment to be
    //    moved (B, C and D)
    if (len_B >= PROM_SIZE) {
      gen_unit_->extract_promoters_included_in(pos_B, pos_C, promoters_B);
    }
    if (len_C >= PROM_SIZE) {
      gen_unit_->extract_promoters_included_in(pos_C, pos_D, promoters_C);
    }
    if (len_D >= PROM_SIZE) {
      gen_unit_->extract_promoters_included_in(pos_D, pos_E, promoters_D);
    }

    // 3a) Invert promoters of segments B and C
    GeneticUnit::invert_promoters(promoters_B, pos_B, pos_C);
    GeneticUnit::invert_promoters(promoters_C, pos_C, pos_D);

    // 3b) Shift these promoters positions
    GeneticUnit::shift_promoters(promoters_B, len_D,
                                 gen_unit_->dna()->length());
    GeneticUnit::shift_promoters(promoters_C, len_D,
                                 gen_unit_->dna()->length());
    GeneticUnit::shift_promoters(promoters_D, -len_B - len_C,
                                 gen_unit_->dna()->length());

    // 4) Reinsert the shifted promoters
    gen_unit_->insert_promoters(promoters_C);
    gen_unit_->insert_promoters(promoters_B);
    gen_unit_->insert_promoters(promoters_D);

    // 5) Look for new promoters including a breakpoint
    gen_unit_->look_for_new_promoters_around(len_A);
    gen_unit_->look_for_new_promoters_around(len_AD);
    gen_unit_->look_for_new_promoters_around(len_ADB);
    gen_unit_->look_for_new_promoters_around(len_ADBC);
  }
}

void Dna::ABCDE_to_ACpDpBE(int32_t pos_B, int32_t pos_C, int32_t pos_D,
                           int32_t pos_E) {
  // Rearrange the sequence from ABCDE to ACpDpBE (complex translocation with
  // inversion of segment defined between positions pos_C and pos_E)
  // Cp (resp Dp) stands for inverted C (resp D)
  //
  // Segments are identified by pos_x values as shown below.
  //
  // TODO(dpa) CHECK THIS !!!
  // WARNING : Segment D includes nucleotide at pos_E // NOTE : WTF???
  //
  //         A      B        C       D       E
  //      |----<]-------->=======[>=======>-------|
  //          pos_B    pos_C    pos_D   pos_E
  //
  //                         |
  //                         V
  //
  //          A       C'      D'       B       E
  //       |-----<]=======>=======<]------->-------|


  // Check points' order and range
  assert(pos_B >= 0 && pos_C >= pos_B && pos_D >= pos_C && pos_E >= pos_D &&
                                                           pos_E <= length_);

  // Compute segment lengths
  int32_t len_A = pos_B;
  int32_t len_B = pos_C - pos_B;
  int32_t len_C = pos_D - pos_C;
  int32_t len_D = pos_E - pos_D;
  int32_t len_E = length_ - pos_E;
  int32_t len_AC = len_A + len_C;
  int32_t len_ACD = len_AC + len_D;
  int32_t len_ACDB = len_ACD + len_B;

  // Create new sequence
  char* new_genome = new char[nb_blocks_ * BLOCK_SIZE];

  // Copy segment A
  memcpy(new_genome, data_, len_A * sizeof(char));


  // Build Cp and put it in the new genome
  char* inverted_segment = new char[len_C + 1];

#ifdef __SIMD
#pragma simd
#pragma distribute_point
#endif
  for (int32_t i = 0, j = pos_D - 1; i < len_C; i++, j--) {
    #ifdef BASE_2
    if (data_[j] == '0') inverted_segment[i] = '1';
    else inverted_segment[i] = '0';
    #elif BASE_4
    inverted_segment[i] = get_complementary_base(data_[j]);
    #endif
  }
  inverted_segment[len_C] = '\0';

  memcpy(&new_genome[len_A], inverted_segment, len_C * sizeof(char));

  delete[] inverted_segment;


  // Build Dp and put it in the new genome
  inverted_segment = new char[len_D + 1];

#ifdef __SIMD
#pragma simd
#pragma distribute_point
#endif
  for (int32_t i = 0, j = pos_E - 1; i < len_D; i++, j--) {
    #ifdef BASE_2
    if (data_[j] == '0') inverted_segment[i] = '1';
    else inverted_segment[i] = '0';
    #elif BASE_4
    inverted_segment[i] = get_complementary_base(data_[j]);
    #endif
  }
  inverted_segment[len_D] = '\0';

  memcpy(&new_genome[len_AC], inverted_segment, len_D * sizeof(char));

  delete[] inverted_segment;

  // Copy segments B and E
  memcpy(&new_genome[len_ACD], &data_[pos_B], len_B * sizeof(char));
  memcpy(&new_genome[len_ACDB], &data_[pos_E], len_E * sizeof(char));
  new_genome[length_] = '\0';


  // Replace sequence
  delete[] data_;
  data_ = new_genome;


  // ========== Update promoter list ==========
  // 1) Remove promoters that include a breakpoint
  // 2) Extract promoters that are totally included in each segment to be
  //    moved (B, C and D)
  // 3) Shift (and invert when needed) these promoters positions
  // 4) Reinsert the shifted promoters
  // 5) Look for new promoters including a breakpoint
  if (length_ >= PROM_SIZE) {
    // 1) Remove promoters that include a breakpoint
    gen_unit_->remove_promoters_around(pos_B);
    gen_unit_->remove_promoters_around(pos_C);
    gen_unit_->remove_promoters_around(pos_D);
    gen_unit_->remove_promoters_around(pos_E);

    // Create temporary lists for promoters to move and/or invert
    Promoters2Strands promoters_B = {{},
                                     {}};
    Promoters2Strands promoters_C = {{},
                                     {}};
    Promoters2Strands promoters_D = {{},
                                     {}};

    // 2) Extract promoters that are totally included in each segment to be
    //    moved (B, C and D)
    if (len_B >= PROM_SIZE) {
      gen_unit_->extract_promoters_included_in(pos_B, pos_C, promoters_B);
    }
    if (len_C >= PROM_SIZE) {
      gen_unit_->extract_promoters_included_in(pos_C, pos_D, promoters_C);
    }
    if (len_D >= PROM_SIZE) {
      gen_unit_->extract_promoters_included_in(pos_D, pos_E, promoters_D);
    }

    // 3a) Invert promoters of segments C and D
    GeneticUnit::invert_promoters(promoters_C, pos_C, pos_D);

    GeneticUnit::invert_promoters(promoters_D, pos_D, pos_E);

    // 3b) Shift these promoters positions
    GeneticUnit::shift_promoters(promoters_B, len_C + len_D,
                                 gen_unit_->dna()->length());

    GeneticUnit::shift_promoters(promoters_C, -len_B,
                                 gen_unit_->dna()->length());

    GeneticUnit::shift_promoters(promoters_D, -len_B,
                                 gen_unit_->dna()->length());

    // 4) Reinsert the shifted promoters
    gen_unit_->insert_promoters(promoters_B);
    gen_unit_->insert_promoters(promoters_D);
    gen_unit_->insert_promoters(promoters_C);

    // 5) Look for new promoters including a breakpoint
    gen_unit_->look_for_new_promoters_around(len_A);
    gen_unit_->look_for_new_promoters_around(len_AC);
    gen_unit_->look_for_new_promoters_around(len_ACD);
    gen_unit_->look_for_new_promoters_around(len_ACDB);
  }
}

void Dna::inter_GU_ABCDE_to_ACDBE(int32_t pos_B, int32_t pos_C, int32_t pos_E) {
  // Check points' order and range
  assert((pos_B >= 0 && pos_C >= pos_B) && (pos_E >= 0));

  if (pos_B != pos_C) {
    // Useful values
    Individual* indiv = indiv_;
    GeneticUnit& chromosome = indiv->genetic_unit_nonconst(0);
    GeneticUnit& plasmid = indiv->genetic_unit_nonconst(1);
    GeneticUnit& destination_GU = (gen_unit_ == &chromosome) ?
                                  plasmid : chromosome;

    // Compute segment lengths
    int32_t len_A = pos_B;
    int32_t len_B = pos_C - pos_B;
    int32_t len_C = length_ - pos_C;
    int32_t len_D = pos_E;
    int32_t len_E = destination_GU.dna()->length() - pos_E;
    int32_t len_AC = len_A + len_C;
    int32_t len_DB = len_D + len_B;
    int32_t len_DBE = len_DB + len_E;


    // Create the new sequence of this genetic unit
    int32_t tmp = ae_string::nb_blocks(len_AC);
    char* new_sequence_this = new char[tmp * BLOCK_SIZE];

    memcpy(new_sequence_this, data_, len_A * sizeof(char));
    memcpy(&new_sequence_this[len_A], &data_[pos_C], len_C * sizeof(char));
    new_sequence_this[len_AC] = '\0';

    // Create the new sequence of the destination genetic unit
    tmp = ae_string::nb_blocks(len_DBE) * BLOCK_SIZE;
    char* new_sequence_dest = new char[tmp];
    const char* dest_GU_former_seq = destination_GU.dna()->data();

    memcpy(new_sequence_dest, dest_GU_former_seq, len_D * sizeof(char));
    memcpy(&new_sequence_dest[len_D], &data_[pos_B], len_B * sizeof(char));
    memcpy(&new_sequence_dest[len_DB], &dest_GU_former_seq[pos_E],
           len_E * sizeof(char));
    new_sequence_dest[len_DBE] = '\0';



    // ========== Update promoter list ==========
    // 1) Remove promoters that include a breakpoint
    // 2) Extract promoters that are totally included in each segment to be
    //    moved (B, C and E)
    // NB : Sequences have to be updated at this stage in order to have the
    // correct lengths when managing new promoter positions
    // ........
    // 3) Shift these promoters positions
    // 4) Reinsert the shifted promoters
    // 5) Look for new promoters including a breakpoint


    // 1) Remove promoters that include a breakpoint
    gen_unit_->remove_promoters_around(pos_B);
    gen_unit_->remove_promoters_around(pos_C);
    destination_GU.remove_promoters_around(pos_E);

    // Create temporary lists for promoters to move and/or invert
    Promoters2Strands promoters_B = {{},
                                     {}};

    // 2) Extract promoters that are totally included in each segment to be
    // moved (B, C and E)
    if (len_B >= PROM_SIZE) {
      gen_unit_->extract_promoters_included_in(pos_B, pos_C, promoters_B);
    }

    // ========== Replace sequences ==========
    set_data(new_sequence_this, len_AC);
    destination_GU.dna()->set_data(new_sequence_dest, len_DBE);

    // 3) Shift these promoters positions
    GeneticUnit::shift_promoters(promoters_B, len_D - len_A,
                                 destination_GU.dna()->length());

    // Reassign promoters to their new genetic unit
    for (auto& strand : {LEADING, LAGGING})
      for (auto& rna : promoters_B[strand])
        rna.set_genetic_unit(&destination_GU);

    // Shift the promoters of sequences C and E
    gen_unit_->move_all_promoters_after(pos_C, -len_B);
    destination_GU.move_all_promoters_after(pos_E, len_B);

    // 4) Reinsert the shifted promoters
    destination_GU.insert_promoters(promoters_B);

    // 5) Look for new promoters including a breakpoint
    gen_unit_->look_for_new_promoters_around(0);
    gen_unit_->look_for_new_promoters_around(len_A);
    destination_GU.look_for_new_promoters_around(0);
    destination_GU.look_for_new_promoters_around(len_D);
    destination_GU.look_for_new_promoters_around(len_DB);
  }
}


void Dna::inter_GU_ABCDE_to_BDCAE(int32_t pos_B, int32_t pos_C, int32_t pos_E) {
  // Check points' order and range
  assert((pos_B >= 0 && pos_C >= pos_B) && (pos_E >= 0));

  // Compute segment lengths
  int32_t len_A = pos_B;
  int32_t len_B = pos_C - pos_B;
  int32_t len_C = length_ - pos_C;
  //~ int32_t len_ABC = length_;
  int32_t len_DA = len_A + pos_E;

  inter_GU_ABCDE_to_ACDBE(0, pos_B, pos_E);
  inter_GU_ABCDE_to_ACDBE(len_B, (len_B + len_C), len_DA);
}

void Dna::apply_mutations() {
  int32_t segment_length;
  Mutation* mut = nullptr;

  MutationEvent* repl = nullptr;
// int idx = 0;

//     if (indiv_->grid_cell()->x() *
//                                         exp_m_->world()->height() +
//                                     indiv_->grid_cell()->y() == 6) {
//         for (auto rna : indiv_->genetic_unit(0).rna_list()[LEADING]) {
//             printf("RNA CPU %d Start %d Stop %d Leading/Lagging %d Length %d Basal %lf\n", idx,
//                    rna.promoter_pos(), rna.last_transcribed_pos(), rna.strand(),
//                    rna.transcript_length(),rna.basal_level());
//           idx++;
//         }
         
//          for (auto rna : indiv_->genetic_unit(0).rna_list()[LAGGING]) {
//             printf("RNA CPU %d Start %d Stop %d Leading/Lagging %d Length %d Basal %lf\n", idx,
//                    rna.promoter_pos(), rna.last_transcribed_pos(), rna.strand(),
//                    rna.transcript_length(),rna.basal_level());
//           idx++;
//         }                           }
  do {
    repl = exp_m_
               ->dna_mutator_array_[indiv_->grid_cell()->x() *
                                        exp_m_->world()->height() +
                                    indiv_->grid_cell()->y()]
               ->generate_next_mutation(length_);


    if (repl != nullptr) {
        // if (indiv_->grid_cell()->x() *
        //                                 exp_m_->world()->height() +
        //                             indiv_->grid_cell()->y() == 6)
        //     printf("CPU_Mutation %d : POS %d %d\n",repl->type(),repl->pos_1(),repl->pos_2());

      switch (repl->type()) {
        case DO_SWITCH:
//          printf("%d -- %d -- Switch at %d\n",AeTime::time(),indiv()->id(),repl->pos_1());
          #ifdef BASE_2
          mut = new PointMutation(repl->pos_1());
          do_switch(repl->pos_1());
          #elif BASE_4
          do_switch(repl->pos_1(), repl->base());
          mut = new PointMutation(repl->pos_1(), repl->base());
          #endif
          break;
        case SMALL_INSERTION:
//          printf("%d -- %d -- Insertion at %d size %d\n",AeTime::time(),indiv()->id(),repl->pos_1(),repl->number());

              mut = new SmallInsertion(repl->pos_1(), repl->number(), repl->seq());
              do_small_insertion(repl->pos_1(), repl->number(), repl->seq());

          break;
        case SMALL_DELETION:
//          printf("%d -- %d -- Deletion at %d size %d\n",AeTime::time(),indiv()->id(),repl->pos_1(),repl->number());
          mut = new SmallDeletion(repl->pos_1(), repl->number());
          do_small_deletion(repl->pos_1(), repl->number());
          break;
        case DUPLICATION:
          segment_length =
              Utils::mod(repl->pos_2() - repl->pos_1() - 1, length_) + 1;

          mut = new Duplication(repl->pos_1(), repl->pos_2(), repl->pos_3(),
                                segment_length);
          do_duplication(repl->pos_1(), repl->pos_2(), repl->pos_3());


          break;
        case TRANSLOCATION:
          segment_length = repl->pos_2() - repl->pos_1();
        //  printf("%d -- Translocation pos_1 %d pos_2 %d pos_3 %d pos_4 %d seg_lengh %d\n",
        //         indiv()->id(),repl->pos_1(),repl->pos_2(),repl->pos_3(),
        //         repl->pos_4(),segment_length);

          mut = new Translocation(repl->pos_1(), repl->pos_2(), repl->pos_3(),
                                  repl->pos_4(), segment_length,
                                  repl->invert());
          do_translocation(repl->pos_1(), repl->pos_2(), repl->pos_3(),
                           repl->pos_4(), repl->invert());
          break;
        case INVERSION:
          segment_length = repl->pos_2() - repl->pos_1();
//          printf("%d -- %d -- Inversion pos_1 %d pos_2 %d seg_lengh %d\n",AeTime::time(),
//                 indiv()->id(),repl->pos_1(),repl->pos_2(),segment_length);
          mut = new Inversion(repl->pos_1(), repl->pos_2(), segment_length);
          do_inversion(repl->pos_1(), repl->pos_2());
          break;
        case DELETION:
          segment_length =
              Utils::mod(repl->pos_2() - repl->pos_1() - 1, length_) + 1;
//          printf("%d -- %d -- Deletion pos_1 %d pos_2 %d seg_lengh %d length %d\n",
//                  AeTime::time(),indiv()->id(),repl->pos_1(),repl->pos_2(),segment_length, length());
          mut = new Deletion(repl->pos_1(), repl->pos_2(), segment_length);
          do_deletion(repl->pos_1(), repl->pos_2());
          break;
      }
      if (mut != nullptr) {
#pragma omp critical
        {
          //indiv_->notifyObservers(MUTATION, mut);
          if (exp_m_->record_tree() || exp_m_->light_tree()) {
            indiv_->exp_m_->tree()->report_by_index(AeTime::time(),indiv_->grid_cell()->x() *
                                                                  indiv_->exp_m()->grid_height()
                                                                  + indiv_->grid_cell()->y())->dna_replic_report().add_mut(mut);
          }
          delete mut;
        }
      }
    }

  } while (exp_m_
               ->dna_mutator_array_[indiv_->grid_cell()->x() *
                                        exp_m_->world()->height() +
                                    indiv_->grid_cell()->y()]
               ->mutation_available() > 0);
    // if (indiv_->grid_cell()->x() *
    //                                     exp_m_->world()->height() +
    //                                 indiv_->grid_cell()->y() == 6) {
    //     for (auto rna : indiv_->genetic_unit(0).rna_list()[LEADING]) {
    //         printf("RNA CPU %d Start %d Stop %d Leading/Lagging %d Length %d Basal %lf\n", idx,
    //                rna.promoter_pos(), rna.last_transcribed_pos(), rna.strand(),
    //                rna.transcript_length(),rna.basal_level());
    //       idx++;
    //     }
         
    //      for (auto rna : indiv_->genetic_unit(0).rna_list()[LAGGING]) {
    //         printf("RNA CPU %d Start %d Stop %d Leading/Lagging %d Length %d Basal %lf\n", idx,
    //                rna.promoter_pos(), rna.last_transcribed_pos(), rna.strand(),
    //                rna.transcript_length(),rna.basal_level());
    //       idx++;
    //     }                           }
}

}// namespace aevol
