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
#include "PhenotypicTarget_R.h"
#include "FuzzyFactory.h"

#include <cstring>


namespace aevol {


//##############################################################################
//                                                                             #
//                           Class PhenotypicTarget                            #
//                                                                             #
//##############################################################################

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================
PhenotypicTarget_R::PhenotypicTarget_R() : PhenotypicTarget() {
  id_ = 0;
  //signals_ = NULL;
}

PhenotypicTarget_R::PhenotypicTarget_R( int16_t id) : PhenotypicTarget() {
  id_ = id;
  //signals_ = NULL;
}

PhenotypicTarget_R::PhenotypicTarget_R(const PhenotypicTarget_R& rhs) : PhenotypicTarget(rhs) {
  id_ = rhs.id_;
  signals_ = rhs.signals_;
}

// ============================================================================
//                                 Destructor
// ============================================================================
PhenotypicTarget_R::~PhenotypicTarget_R() {

}

// ============================================================================
//                                   Methods
// ============================================================================
  void PhenotypicTarget_R::save(gzFile backup_file) const {
    //printf("Appel a la sauvegarde de PhenotypicTarget_R\n");
    SaveSegmentation(backup_file);
    // Backup id
    gzwrite(backup_file, &id_, sizeof(id_));
  }

  void PhenotypicTarget_R::load(gzFile backup_file) {
    //printf("Appel au chargement de PhenotypicTarget_R\n");
    LoadSegmentation(backup_file);
    //  Retrieve id
    gzread(backup_file, &id_, sizeof(id_));
  }

// ============================================================================
//                            Non inline accessors
// ============================================================================
} // namespace aevol
