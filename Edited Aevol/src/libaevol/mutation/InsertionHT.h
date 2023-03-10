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

#ifndef AEVOL_INSERTION_HT_H_
#define AEVOL_INSERTION_HT_H_


// ============================================================================
//                                   Includes
// ============================================================================

#include "HorizontalTransfer.h"
#include "VisAVis.h"

namespace aevol {

/**
 *
 */
class InsertionHT : public HorizontalTransfer {

 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  InsertionHT() = default; //< Default ctor
  InsertionHT(const InsertionHT&) = default; //< Copy ctor
  // The move constructor is implicitly deleted.
  InsertionHT(VisAVis& donor_donor_align, VisAVis& exo_recv_align,
              int32_t length, char* seq, int32_t donor_id);

  virtual Mutation* Clone() const override { return new InsertionHT(*this); };

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~InsertionHT() noexcept; //< Destructor

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
    return INS_HT;
  };

  int32_t donor_pos1() const {
    return donor_donor_align_.i_1();
  }

  int32_t donor_pos2() const {
    return donor_donor_align_.i_2();
  }

  int32_t exogenote_pos() const {
    return exo_recv_align_.i_1();
  }

  int32_t receiver_pos() const {
    return exo_recv_align_.i_2();
  }

  AlignmentSense sense() const {
    return exo_recv_align_.sense();
  }

  int32_t length() const {
    return length_;
  }

  char* seq() const {
    return seq_;
  }

  int32_t donor_id() const {
    return donor_id_;
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
  VisAVis donor_donor_align_, exo_recv_align_;
  int32_t length_;
  char* seq_ = nullptr;
  int32_t donor_id_ = -1;
};

} // namespace aevol
#endif //AEVOL_INSERTION_HT_H_
