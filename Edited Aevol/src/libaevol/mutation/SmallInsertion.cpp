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
#include "SmallInsertion.h"

#include <cstring>


namespace aevol {

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================
SmallInsertion::SmallInsertion(const SmallInsertion& other) :
    pos_(other.pos_),
    length_(other.length_) {
  seq_ = new char[length_ + 1];
  memcpy(seq_, other.seq_, length_ + 1);
}

SmallInsertion::SmallInsertion(int32_t pos, int32_t length, const char* seq) :
    pos_(pos),
    length_(length) {
  seq_ = new char[length_ + 1];
  memcpy(seq_, seq, length_);
  seq_[length_] = '\0';
}

// ============================================================================
//                                 Destructor
// ============================================================================
SmallInsertion::~SmallInsertion() noexcept {
  delete [] seq_;
}

// ============================================================================
//                                 Operators
// ============================================================================


// ============================================================================
//                                   Methods
// ============================================================================
void SmallInsertion::save(gzFile backup_file) const {
  int8_t tmp_mut_type = S_INS;
  gzwrite(backup_file, &tmp_mut_type, sizeof(tmp_mut_type));
  gzwrite(backup_file, &pos_, sizeof(pos_));
  gzwrite(backup_file, &length_, sizeof(length_));
  gzwrite(backup_file, seq_, length_ * sizeof(seq_[0]));
}

void SmallInsertion::load(gzFile backup_file) {
  gzread(backup_file, &pos_, sizeof(pos_));
  gzread(backup_file, &length_, sizeof(length_));
  seq_ = new char[length_ + 1];
  gzread(backup_file, seq_, length_ * sizeof(seq_[0]));
  seq_[length_] = '\0';
}

void SmallInsertion::generic_description_string(char* str) const {
  sprintf(str, "%" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32
      " %" PRId8 " %" PRId16 " %" PRId16 " %" PRId32 " %" PRId32,
          mut_type(), pos_, -1, -1, -1, (int8_t) -1, (int16_t) -1, (int16_t) -1,
          length_, -1);
}

// ============================================================================
//                            Non inline accessors
// ============================================================================

} // namespace aevol
