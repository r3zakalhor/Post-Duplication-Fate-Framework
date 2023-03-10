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


#ifndef AEVOL_HABITAT_H_
#define AEVOL_HABITAT_H_


// ============================================================================
//                                   Includes
// ============================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <memory>

#include <zlib.h>

#include "PhenotypicTargetHandler.h"


namespace aevol {

// ============================================================================
//                          Class declarations
// ============================================================================






class Habitat
{
 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  Habitat(); //< Default ctor
  Habitat(const Habitat&) = delete; //< Copy ctor
  Habitat(Habitat&&) = delete; //< Move ctor
  Habitat(const Habitat&, bool share_phenotypic_target);
  Habitat(gzFile backup_file,
          std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler);

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~Habitat() = default; //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================
  Habitat& operator=(const Habitat&) = default;
  Habitat& operator=(Habitat&&) = default;

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  virtual void ApplyVariation();
  void save(gzFile backup_file,
            bool skip_phenotypic_target = false) const;
  void load(gzFile backup_file,
            std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler);

  // ==========================================================================
  //                                 Getters
  // ==========================================================================
  double compound_amount() const {return compound_amount_;};
  const PhenotypicTarget& phenotypic_target() const {
    return phenotypic_target_handler_->phenotypic_target();
  }
  virtual const PhenotypicTargetHandler& phenotypic_target_handler() const {
    return *phenotypic_target_handler_;
  }
  virtual PhenotypicTargetHandler& phenotypic_target_handler_nonconst() const {
    return *phenotypic_target_handler_;
  }

  double mean_environmental_area() const {
    return phenotypic_target_handler_->mean_environmental_area();
  }

  // ==========================================================================
  //                                 Setters
  // ==========================================================================
  void set_compound_amount(double compound_amount) {
    compound_amount_ = compound_amount;
  };

  void set_phenotypic_target_handler(std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler) {
    phenotypic_target_handler_ = phenotypic_target_handler;
  };

 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================

  // ==========================================================================
  //                               Attributes
  // ==========================================================================
  // Amount of secreted compound currently present in the grid cell
  double compound_amount_ = 0;

  /** Handler for the phenotypic target and its "evolution" over time */
  std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler_;
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

#endif // AEVOL_HABITAT_H_
