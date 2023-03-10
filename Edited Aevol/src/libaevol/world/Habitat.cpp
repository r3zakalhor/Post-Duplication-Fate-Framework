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




// ============================================================================
//                                   Includes
// ============================================================================
#if __cplusplus == 201103L
#include "make_unique.h"
#endif

#include "Habitat.h"

#include <iostream>


using std::cout;
using std::endl;


namespace aevol {


//##############################################################################
//                                                                             #
//                                Class Habitat                                #
//                                                                             #
//##############################################################################

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================
Habitat::Habitat() {
  compound_amount_ = 0.0;
  phenotypic_target_handler_ = std::make_shared<PhenotypicTargetHandler>();
}

Habitat::Habitat(const Habitat& rhs, bool share_phenotypic_target) {
  assert(share_phenotypic_target);
  compound_amount_ = rhs.compound_amount_;
  phenotypic_target_handler_ = rhs.phenotypic_target_handler_;
}

Habitat::Habitat(gzFile backup_file,
                 std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler) {
  load(backup_file, phenotypic_target_handler);
}

// ============================================================================
//                                 Destructor
// ============================================================================

// ============================================================================
//                                   Methods
// ============================================================================
void Habitat::ApplyVariation() {
  //printf("Appel au apply_variation de habitat\n");
  phenotypic_target_handler_->ApplyVariation();
}

void Habitat::save(gzFile backup_file,
                   bool skip_phenotypic_target /*=false*/) const {
  //printf("Appel a la sauvegarde de Habitat\n");
  gzwrite(backup_file, &compound_amount_, sizeof(compound_amount_));
  if (not skip_phenotypic_target)
    phenotypic_target_handler_->save(backup_file);
}

void Habitat::load(gzFile backup_file,
                   std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler) {
  //printf("Appel au chargement de Habitat\n");
  gzread(backup_file, &compound_amount_, sizeof(compound_amount_));
  if (phenotypic_target_handler == nullptr)
    phenotypic_target_handler_ = std::make_shared<PhenotypicTargetHandler>(backup_file);
  else
    phenotypic_target_handler_ = phenotypic_target_handler;
}

// ============================================================================
//                            Non inline accessors
// ============================================================================

} // namespace aevol
