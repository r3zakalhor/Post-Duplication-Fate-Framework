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

#ifndef AEVOL_INVERSION_H_
#define AEVOL_INVERSION_H_


// ============================================================================
//                                   Includes
// ============================================================================

#include "Rearrangement.h"

namespace aevol {

/**
 *
 */
class Inversion : public Rearrangement{

 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  Inversion() = default; //< Default ctor
  Inversion(const Inversion&) = default; //< Copy ctor
  // The move constructor is implicitly deleted.
  Inversion(int32_t pos1, int32_t pos2,
            int32_t length, int16_t align_score = -1);

  virtual Mutation* Clone() const override { return new Inversion(*this); };

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~Inversion() noexcept = default; //< Destructor

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
    return INV;
  };

  int32_t pos1() const {
    return pos1_;
  }

  int32_t pos2() const {
    return pos2_;
  }

  int32_t length() const {
    return length_;
  }

  int16_t align_score() const {
    return align_score_;
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
  int32_t pos1_, pos2_;
  int32_t length_;
  int16_t align_score_ = -1;
};

} // namespace aevol
#endif //AEVOL_INVERSION_H_
