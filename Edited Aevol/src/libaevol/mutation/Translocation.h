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

#ifndef AEVOL_TRANSLOCATION_H_
#define AEVOL_TRANSLOCATION_H_


// ============================================================================
//                                   Includes
// ============================================================================

#include "Rearrangement.h"

namespace aevol {

/**
 *
 */
class Translocation : public Rearrangement {

 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  Translocation() = default; //< Default ctor
  Translocation(const Translocation&) = default; //< Copy ctor
  // The move constructor is implicitly deleted.

  Translocation(int32_t pos1, int32_t pos2, int32_t pos3, int32_t pos4,
                int32_t length,
                bool invert,
                int16_t align_score_1 = -1, int16_t align_score_2 = -1);

  virtual Mutation* Clone() const override { return new Translocation(*this); };

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~Translocation() noexcept = default; //< Destructor

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
    return TRANS;
  };

  int32_t pos1() const {
    return pos1_;
  }

  int32_t pos2() const {
    return pos2_;
  }

  int32_t pos3() const {
    return pos3_;
  }

  int32_t pos4() const {
    return pos4_;
  }

  int32_t length() const {
    return length_;
  }

  bool invert() const {
    return invert_;
  }

  int16_t align_score_1() const {
    return align_score_1_;
  }

  int16_t align_score_2() const {
    return align_score_2_;
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
  int32_t pos1_, pos2_, pos3_, pos4_;
  int32_t length_;
  bool invert_;
  int16_t align_score_1_ = -1;
  int16_t align_score_2_ = -1;
};

} // namespace aevol
#endif //AEVOL_TRANSLOCATION_H_
