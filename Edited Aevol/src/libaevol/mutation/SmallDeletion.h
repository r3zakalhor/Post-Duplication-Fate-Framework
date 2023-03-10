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

#ifndef AEVOL_SMALLDELETION_H_
#define AEVOL_SMALLDELETION_H_


// ============================================================================
//                                   Includes
// ============================================================================
#include <cstdint>

#include "LocalMutation.h"

namespace aevol {

/**
 *
 */
class SmallDeletion : public LocalMutation {

 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  SmallDeletion() = default; //< Default ctor
  SmallDeletion(const SmallDeletion&) = default; //< Copy ctor
  // The move constructor is implicitly deleted.
  SmallDeletion(int32_t pos, int16_t length);

  virtual Mutation* Clone() const override { return new SmallDeletion(*this); };

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~SmallDeletion() noexcept = default; //< Destructor

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
    return S_DEL;
  };

  int32_t pos() const {
    return pos_;
  };

  int16_t length() const {
    return length_;
  };

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
  int32_t pos_;
  int16_t length_;
};

} // namespace aevol
#endif //AEVOL_SMALLDELETION_H_
