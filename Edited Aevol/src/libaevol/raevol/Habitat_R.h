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


#ifndef AEVOL_HABITAT_R_H__
#define AEVOL_HABITAT_R_H__


// ============================================================================
//                                   Includes
// ============================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <memory>

#include <zlib.h>

#include "PhenotypicTargetHandler_R.h"
#include "Habitat.h"


namespace aevol {

// ============================================================================
//                          Class declarations
// ============================================================================






class Habitat_R : public virtual Habitat
{
 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  Habitat_R(void); //< Default ctor
  Habitat_R(const Habitat_R&) = delete; //< Copy ctor
  Habitat_R(Habitat_R&&) = delete; //< Move ctor
  Habitat_R(const Habitat_R&, bool share_phenotypic_target);
  Habitat_R(gzFile backup_file,
            std::shared_ptr<PhenotypicTargetHandler_R> phenotypic_target_handler_);

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~Habitat_R(void) = default; //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  virtual void ApplyVariation();
  virtual void save(gzFile backup_file,
                    bool skip_phenotypic_target = false) const;
  virtual void load(gzFile backup_file,
                    std::shared_ptr<PhenotypicTargetHandler_R> phenotypic_target_handler);
  // ==========================================================================
  //                                 Getters
  // ==========================================================================
  /*virtual const PhenotypicTarget& phenotypic_target() const {
    Utils::ExitWithDevMsg("You should not call a phenotypic target without age id in RAevol", __FILE__, __LINE__);
    return ;
  }*/

  const PhenotypicTarget_R& phenotypic_target(  int16_t age ) const {
    return (std::dynamic_pointer_cast<PhenotypicTargetHandler_R>(phenotypic_target_handler_))->phenotypic_target( age );
  }

    const PhenotypicTarget_R& phenotypic_target_model(  int16_t env_id ) const {
        return (std::dynamic_pointer_cast<PhenotypicTargetHandler_R>(phenotypic_target_handler_))->phenotypic_target_model( env_id );
    }

  int16_t number_of_phenotypic_targets() const {
    return (std::dynamic_pointer_cast<PhenotypicTargetHandler_R>(phenotypic_target_handler_))->number_of_phenotypic_targets();
  }

    int16_t number_of_phenotypic_target_models() const {
      return (std::dynamic_pointer_cast<PhenotypicTargetHandler_R>(phenotypic_target_handler_))->number_of_phenotypic_target_models();
    }

  virtual const PhenotypicTargetHandler_R& phenotypic_target_handler() const {
    return *(std::dynamic_pointer_cast<PhenotypicTargetHandler_R> (phenotypic_target_handler_));
  }
  virtual PhenotypicTargetHandler_R& phenotypic_target_handler_nonconst() const {
    return *(std::dynamic_pointer_cast<PhenotypicTargetHandler_R> (phenotypic_target_handler_));;
  }

    virtual std::shared_ptr<PhenotypicTargetHandler_R> phenotypic_target_handler_ptr() const {
      return (std::dynamic_pointer_cast<PhenotypicTargetHandler_R> (phenotypic_target_handler_));;
    }

  const std::list<Protein_R*> signals() const {
    return (std::dynamic_pointer_cast<PhenotypicTargetHandler_R>(phenotypic_target_handler_))->signals();
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

#endif // AEVOL_HABITAT_R_H__
