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
#include "ReplacementHT.h"

namespace aevol {

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================
ReplacementHT::ReplacementHT(const VisAVis& align1, const VisAVis& align2,
                             int32_t length, int32_t replaced_seq_length,
                             char* seq, int32_t donor_id)
    : align1_(align1), align2_(align2),
      length_(length), replaced_seq_length_(replaced_seq_length),
      donor_id_(donor_id) {
  seq_ = new char[length_ + 1];
  memcpy(seq_, seq, length_ + 1);
}

// ============================================================================
//                                 Destructor
// ============================================================================
ReplacementHT::~ReplacementHT() noexcept {
  delete seq_;
}

// ============================================================================
//                                   Methods
// ============================================================================
void ReplacementHT::save(gzFile backup_file) const {
  Utils::ExitWithDevMsg("Not implemented yet", __FILE__, __LINE__);
}

void ReplacementHT::load(gzFile backup_file) {
  Utils::ExitWithDevMsg("Not implemented yet", __FILE__, __LINE__);
}

void ReplacementHT::generic_description_string(char* str) const {
  sprintf(str, "%" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32
      " %" PRId32 " %" PRId16 " %" PRId16 " %" PRId32 " %" PRId32 " ",
          mut_type(),
          donor_pos1(), donor_pos2(),
          receiver_pos1(), receiver_pos2(),
          sense(),
          align1_.score(), align2_.score(),
          length_, replaced_seq_length_);
}

// ============================================================================
//                            Non inline accessors
// ============================================================================

} // namespace aevol
