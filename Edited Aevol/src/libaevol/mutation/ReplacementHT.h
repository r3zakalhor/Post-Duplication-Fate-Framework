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

#ifndef AEVOL_REPLACEMENT_HT_H_
#define AEVOL_REPLACEMENT_HT_H_


// ============================================================================
//                                   Includes
// ============================================================================

#include "HorizontalTransfer.h"
#include "VisAVis.h"

namespace aevol {

/**
 *
 */
class ReplacementHT : public HorizontalTransfer {

 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  ReplacementHT() = default; //< Default ctor
  ReplacementHT(const ReplacementHT&) = default; //< Copy ctor
  // The move constructor is implicitly deleted.
  ReplacementHT(const VisAVis& align1, const VisAVis& align2,
                int32_t length, int32_t replaced_seq_length,
                char* seq, int32_t donor_id);

  virtual Mutation* Clone() const override { return new ReplacementHT(*this); };

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~ReplacementHT() noexcept; //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================
  // The copy and move operators are implicitly deleted.

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  virtual void save(gzFile backup_file) const override;
  virtual void load(gzFile backup_file) override;
  void generic_description_string(char* str) const override;

  // ==========================================================================
  //                                 Getters
  // ==========================================================================
  virtual MutationType mut_type() const override {
    return REPL_HT;
  };

  int32_t donor_pos1() const {
    return align1_.i_2();
  }

  int32_t donor_pos2() const {
    return align2_.i_2();
  }

  int32_t receiver_pos1() const {
    return align1_.i_1();
  }

  int32_t receiver_pos2() const {
    return align2_.i_1();
  }

  AlignmentSense sense() const {
    return align2_.sense();
  }

  char* seq() const {
    return seq_;
  }

  int32_t length() const {
    return length_;
  }

  // ==========================================================================
  //                                 Setters
  // ==========================================================================

 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================

  // ==========================================================================
  //                               Attributes
  // ==========================================================================
  /**
   * Alignments between the receiver and the donor
   *
   * The first point of each alignment corresponds to the receiver, the sequence
   * that will be replaced lies between align1_->i1 and align2_->i1. It will be
   * replaced by the seq. btw align1_->i2 and align2_->i2 on the donor.
   */
  VisAVis align1_, align2_;
  int32_t length_;
  int32_t replaced_seq_length_;
  char* seq_ = nullptr;
  int32_t donor_id_ = -1;
};

} // namespace aevol
#endif //AEVOL_REPLACEMENT_HT_H_
