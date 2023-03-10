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
#include "PointMutation.h"


namespace aevol {

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================
#ifdef BASE_2
PointMutation::PointMutation(int32_t pos) : pos_(pos) {}
#elif BASE_4
PointMutation::PointMutation(int32_t pos, char choice) : pos_(pos), base_(choice) {}
#endif
// ============================================================================
//                                 Destructor
// ============================================================================

// ============================================================================
//                                   Methods
// ============================================================================
void PointMutation::save(gzFile backup_file) const {
  int8_t tmp_mut_type = SWITCH;
  gzwrite(backup_file, &tmp_mut_type, sizeof(tmp_mut_type));
  gzwrite(backup_file, &pos_, sizeof(pos_));
  #ifdef BASE_4
  gzwrite(backup_file, &base_, sizeof(base_));
  #endif
}

void PointMutation::load(gzFile backup_file) {
  gzread(backup_file, &pos_, sizeof(pos_));
  #ifdef BASE_4
  gzread(backup_file, &base_, sizeof(base_));
  #endif
}

void PointMutation::generic_description_string(char* str) const {
  #ifdef BASE_2
  sprintf(str, "%" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32
      " %" PRId8 " %" PRId16 " %" PRId16 " %" PRId32 " %" PRId32,
          mut_type(), pos_, -1, -1, -1, (int8_t) -1, (int16_t) -1, (int16_t) -1,
          -1, -1);
  #elif BASE_4
  sprintf(str, "%" PRId32 " %" PRId32 " %" PRId8 " %" PRId32 " %" PRId32 " %" PRId32
      " %" PRId8 " %" PRId16 " %" PRId16 " %" PRId32 " %" PRId32,
          mut_type(), pos_, base_, -1, -1, -1, (int8_t) -1, (int16_t) -1, (int16_t) -1,
          -1, -1);
  #endif
}

// ============================================================================
//                            Non inline accessors
// ============================================================================

} // namespace aevol
