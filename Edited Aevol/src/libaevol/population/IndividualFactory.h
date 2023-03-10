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


#ifndef AEVOL_INDIVIDUAL_FACTORY_H_
#define AEVOL_INDIVIDUAL_FACTORY_H_


// ============================================================================
//                                   Includes
// ============================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#ifdef __NO_X
  #ifndef __REGUL
    #include <Individual.h>
  #else
    #include <raevol/Individual_R.h>
  #endif
#elif defined __X11
  #ifndef __REGUL
    #include <Individual_X11.h>
  #else
    #include <raevol/Individual_R_X11.h>
  #endif
#endif

#include "Habitat.h"


namespace aevol {

// ============================================================================
//                          Class declarations
// ============================================================================






class IndividualFactory
{
 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  IndividualFactory() = default; //< Default ctor
  IndividualFactory(const IndividualFactory&) = delete; //< Copy ctor
  IndividualFactory(IndividualFactory&&) = delete; //< Move ctor

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~IndividualFactory() = default; //< Destructor

  // ==========================================================================
  //                                 Getters
  // ==========================================================================

  // ==========================================================================
  //                                 Setters
  // ==========================================================================

  // ==========================================================================
  //                                Operators
  // ==========================================================================

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  static Individual* create_random_individual(
      ExpManager* exp_m,
      int32_t id,
      std::shared_ptr<MutationParams> param_mut,
      std::shared_ptr<JumpingMT> mut_prng,
      std::shared_ptr<JumpingMT> stoch_prng,
      const Habitat& habitat,
      double w_max,
      int32_t min_genome_length,
      int32_t max_genome_length,
      int32_t chromosome_initial_length,
      bool allow_plasmids,
      bool plasmid_initial_gene,
      int32_t plasmid_initial_length,
      char* strain_name,
      std::shared_ptr<JumpingMT> local_prng,
      bool better_than_flat);





 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================

  // ==========================================================================
  //                               Attributes
  // ==========================================================================
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

#endif // AEVOL_INDIVIDUAL_FACTORY_H_
