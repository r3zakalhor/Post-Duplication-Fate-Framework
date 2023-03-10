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
#include "Translocation.h"

namespace aevol {

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================
Translocation::Translocation(int32_t pos1, int32_t pos2, int32_t pos3,
                             int32_t pos4, int32_t length, bool invert,
                             int16_t align_score_1, int16_t align_score_2) :
    pos1_(pos1), pos2_(pos2), pos3_(pos3), pos4_(pos4),
    length_(length), invert_(invert),
    align_score_1_(align_score_1), align_score_2_(align_score_2) {
}

// ============================================================================
//                                 Destructor
// ============================================================================

// ============================================================================
//                                   Methods
// ============================================================================
void Translocation::save(gzFile backup_file) const {
  int8_t tmp_mut_type = TRANS;
  gzwrite(backup_file, &tmp_mut_type, sizeof(tmp_mut_type));
  gzwrite(backup_file, &pos1_, sizeof(pos1_));
  gzwrite(backup_file, &pos2_, sizeof(pos2_));
  gzwrite(backup_file, &pos3_, sizeof(pos3_));
  gzwrite(backup_file, &pos4_, sizeof(pos4_));
  gzwrite(backup_file, &length_, sizeof(length_));
  gzwrite(backup_file, &invert_, sizeof(invert_));
  gzwrite(backup_file, &align_score_1_, sizeof(align_score_1_));
  gzwrite(backup_file, &align_score_2_, sizeof(align_score_2_));
}

void Translocation::load(gzFile backup_file) {
  gzread(backup_file, &pos1_, sizeof(pos1_));
  gzread(backup_file, &pos2_, sizeof(pos2_));
  gzread(backup_file, &pos3_, sizeof(pos3_));
  gzread(backup_file, &pos4_, sizeof(pos4_));
  gzread(backup_file, &length_, sizeof(length_));
  gzread(backup_file, &invert_, sizeof(invert_));
  gzread(backup_file, &align_score_1_, sizeof(align_score_1_));
  gzread(backup_file, &align_score_2_, sizeof(align_score_2_));
}

void Translocation::generic_description_string(char* str) const {
  sprintf(str, "%" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32
      " %" PRId8 " %" PRId16 " %" PRId16 " %" PRId32 " %" PRId32,
          mut_type(), pos1(), pos2(), pos3(), pos4(),
          (int8_t) invert_, align_score_1(), align_score_2(), length_, -1);
}

// ============================================================================
//                            Non inline accessors
// ============================================================================

} // namespace aevol
