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


#ifndef AEVOL_MUTATION_H_
#define AEVOL_MUTATION_H_


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>




// =================================================================
//                            Project Files
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "ae_enums.h"
#include <zlib.h>

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================


enum MutationType {
  // Simple mutation types.
  SWITCH  = 0,
  S_INS,
  S_DEL,
  DUPL,
  DEL,
  TRANS,
  INV,
  INSERT,
  INS_HT,
  REPL_HT,

  // Composite mutation types follow. They represent categories of
  // several simple mutation types. Therefore, they should not be used
  // as array index for counters.
  //
  // The composite mutations should extend MutationType but
  // C++ enums can't be inherited directly.
  S_MUT, // SWITCH or S_INS or S_DEL
  REARR, // DUPL or DEL or TRANS or INV
  H_T,    // INS_HT or REPL_HT
  INDEL  // S_INS or S_DEL
};

class Mutation {
 public :
  // =================================================================
  //                             Constructors
  // =================================================================
  Mutation() = default;
  Mutation(const Mutation& model) = default;
  Mutation(Mutation&& model) = delete;

  virtual Mutation* Clone() const = 0;
  static Mutation* Load(gzFile backup);

  // =================================================================
  //                             Destructor
  // =================================================================
  virtual ~Mutation() noexcept = default;

  // ==========================================================================
  //                                Operators
  // ==========================================================================
  /// Copy assignment
  Mutation& operator=(const Mutation& other) = delete;

  /// Move assignment
  Mutation& operator=(Mutation&& other) = delete;

  // =================================================================
  //                            Public Methods
  // =================================================================
  virtual void save(gzFile backup_file) const = 0;
  virtual void load(gzFile backup_file) = 0;
  virtual void generic_description_string(char* str) const = 0;

  // =================================================================
  //                        Accessors: Getters
  // =================================================================
  virtual MutationType mut_type() const = 0;
  virtual bool is_local_mut() const { return false; };
  virtual bool is_rear() const { return false; };
  virtual bool is_ht() const { return false; };

  // =================================================================
  //                        Accessors: Setters
  // =================================================================

 protected :
  // =================================================================
  //                           Protected Methods
  // =================================================================

  // =================================================================
  //                          Protected Attributes
  // =================================================================
};

// =====================================================================
//                          Accessors' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol
#endif // AEVOL_MUTATION_H_
