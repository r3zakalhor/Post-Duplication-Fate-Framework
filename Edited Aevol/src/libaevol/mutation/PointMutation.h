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

#ifndef AEVOL_POINTMUTATION_H_
#define AEVOL_POINTMUTATION_H_


// ============================================================================
//                                   Includes
// ============================================================================
#include <cstdint>

#include "LocalMutation.h"

namespace aevol {

/**
 *
 */
class PointMutation : public LocalMutation {

 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  PointMutation() = default; //< Default ctor
  PointMutation(const PointMutation&) = default; //< Copy ctor
  // The move constructor is implicitly deleted.
  #ifdef BASE_2
  PointMutation(int32_t pos);
  #elif BASE_4
  PointMutation(int32_t pos, char choice);
  #endif

  virtual Mutation* Clone() const override { return new PointMutation(*this); };

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~PointMutation() noexcept = default; //< Destructor

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
    return SWITCH;
  };

  int32_t pos() const {
    return pos_;
  }

  #ifdef BASE_4
  char base() const {
    return base_;
  }
  #endif
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

  #ifdef BASE_4
  char base_;
  #endif
};

} // namespace aevol
#endif //AEVOL_POINTMUTATION_H_
