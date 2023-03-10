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
#include "InsertionHT.h"

namespace aevol {

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================
InsertionHT::InsertionHT(VisAVis& donor_donor_align, VisAVis& exo_recv_align,
                         int32_t length, char* seq, int32_t donor_id)
    : donor_donor_align_(donor_donor_align), exo_recv_align_(exo_recv_align),
      length_(length), donor_id_(donor_id) {
  assert(donor_donor_align_.sense() == DIRECT);
  seq_ = new char[length_ + 1];
  memcpy(seq_, seq, length_ + 1);
}

// ============================================================================
//                                 Destructor
// ============================================================================
InsertionHT::~InsertionHT() noexcept {
  delete seq_;
}

// ============================================================================
//                                   Methods
// ============================================================================
void InsertionHT::save(gzFile backup_file) const {
  Utils::ExitWithDevMsg("Not implemented yet", __FILE__, __LINE__);
}

void InsertionHT::load(gzFile backup_file) {
  Utils::ExitWithDevMsg("Not implemented yet", __FILE__, __LINE__);
}

void InsertionHT::generic_description_string(char* str) const {
  sprintf(str, "%" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32
      " %" PRId32 " %" PRId16 " %" PRId16 " %" PRId32 " %" PRId32,
          mut_type(),
          donor_pos1(), donor_pos2(), exogenote_pos(), receiver_pos(),
          sense(),
          donor_donor_align_.score(), exo_recv_align_.score(),
          length_, -1);
}

// ============================================================================
//                            Non inline accessors
// ============================================================================

} // namespace aevol
