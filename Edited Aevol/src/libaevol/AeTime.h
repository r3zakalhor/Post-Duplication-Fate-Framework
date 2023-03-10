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


#ifndef AEVOL_TIME_H_
#define AEVOL_TIME_H_


// =================================================================
//                              Includes
// =================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>


namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================





/**
 * The time_ value represents the step that is currently being computed
 *
 * e.g. when creating generation 1 from generation 0, time_ == 1
 */
class AeTime {
 public :
  // =================================================================
  //                             Constructors
  // =================================================================
  AeTime() = delete; //< Default ctor
  AeTime(const AeTime &) = delete; //< Copy ctor
  AeTime(AeTime &&) = delete; //< Move ctor

  // =================================================================
  //                             Destructors
  // =================================================================
  virtual ~AeTime() = delete;

  // =================================================================
  //                        Accessors: getters
  // =================================================================
  static inline int64_t time() {return time_;}

  // =================================================================
  //                        Accessors: setters
  // =================================================================
  static inline void set_time(int64_t t) { time_ = t;}

  // =================================================================
  //                              Operators
  // =================================================================

  // =================================================================
  //                            Public Methods
  // =================================================================
  static inline void plusplus() { time_++;}

  // =================================================================
  //                           Public Attributes
  // =================================================================





 protected :
  // =================================================================
  //                           Protected Methods
  // =================================================================

  // =================================================================
  //                          Protected Attributes
  // =================================================================
  static int64_t time_;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================

// =====================================================================
//                           Setters' definitions
// =====================================================================

// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

int64_t time();

} // namespace aevol

#endif // AEVOL_TIME_H_
