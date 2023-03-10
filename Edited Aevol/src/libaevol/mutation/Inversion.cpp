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
#include "Inversion.h"

namespace aevol {

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================
Inversion::Inversion(int32_t pos1, int32_t pos2, int32_t length,
                     int16_t align_score) :
    pos1_(pos1), pos2_(pos2),
    length_(length), align_score_(align_score) {
}

// ============================================================================
//                                 Destructor
// ============================================================================

// ============================================================================
//                                   Methods
// ============================================================================
void Inversion::save(gzFile backup_file) const {
  int8_t tmp_mut_type = INV;
  gzwrite(backup_file, &tmp_mut_type, sizeof(tmp_mut_type));
  gzwrite(backup_file, &pos1_, sizeof(pos1_));
  gzwrite(backup_file, &pos2_, sizeof(pos2_));
  gzwrite(backup_file, &length_, sizeof(length_));
  gzwrite(backup_file, &align_score_, sizeof(align_score_));
}

void Inversion::load(gzFile backup_file) {
  gzread(backup_file, &pos1_, sizeof(pos1_));
  gzread(backup_file, &pos2_, sizeof(pos2_));
  gzread(backup_file, &length_, sizeof(length_));
  gzread(backup_file, &align_score_, sizeof(align_score_));
}

void Inversion::generic_description_string(char* str) const {
  sprintf(str, "%" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32
      " %" PRId8 " %" PRId16 " %" PRId16 " %" PRId32 " %" PRId32,
          mut_type(), pos1(), pos2(), -1, -1,
          (int8_t) -1, align_score_, (int16_t) -1, length_, -1);
}

// ============================================================================
//                            Non inline accessors
// ============================================================================

} // namespace aevol
