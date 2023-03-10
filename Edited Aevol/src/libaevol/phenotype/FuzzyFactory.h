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
//*****************************************************************************


#ifndef AEVOL_FUZZY_FACTORY_H__
#define AEVOL_FUZZY_FACTORY_H__


// ============================================================================
//                                   Includes
// ============================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "AbstractFuzzy.h"

namespace aevol {

// ============================================================================
//                          Class declarations
// ============================================================================
class ExpSetup;

class FuzzyFactory
{
 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  FuzzyFactory() {}; //< Default ctor
  FuzzyFactory(ExpSetup* exp_s) {_exp_s = exp_s;};

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~FuzzyFactory(void) {}; //< Destructor

  // ==========================================================================
  //                                 Getters
  // ==========================================================================
  int get_fuzzy_flavor();
  // ==========================================================================
  //                                 Setters
  // ==========================================================================

  // ==========================================================================
  //                                Operators
  // ==========================================================================

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
/*  static AbstractFuzzy* create_fuzzy(ExpSetup* exp_s);
  static AbstractFuzzy* create_fuzzy(ExpSetup* exp_s, const AbstractFuzzy& copy);
  static AbstractFuzzy* create_fuzzy(ExpSetup* exp_s, const gzFile backup);*/

  AbstractFuzzy* create_fuzzy();
  AbstractFuzzy* create_fuzzy(const AbstractFuzzy& copy);
  AbstractFuzzy* create_fuzzy(const gzFile backup);

  static FuzzyFactory* fuzzyFactory;
 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================

  // ==========================================================================
  //                               Attributes
  // ==========================================================================
  ExpSetup* _exp_s;
};



// ============================================================================
//                           Getters' definitions
// ============================================================================

// ============================================================================
//                           Setters' definitions
// ============================================================================

// ============================================================================
//                          Operators' definitions
// ============================================================================

// ============================================================================
//                       Inline functions' definition
// ============================================================================

} // namespace aevol

#endif //AEVOL_FUZZYFACTORY_H
