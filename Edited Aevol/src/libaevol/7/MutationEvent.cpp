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


#include "MutationEvent.h"

namespace aevol {

void MutationEvent::switch_pos(int32_t pos
#ifdef BASE_4
, char new_base
#endif
) {
  type_ = MutationEventType::DO_SWITCH;
  pos_1_ = pos;
  #ifdef BASE_4
  base_ = new_base;
  #endif
}

void MutationEvent::small_insertion(int32_t pos, int32_t number, char* seq) {
  type_ = MutationEventType::SMALL_INSERTION;
  pos_1_ = pos;
  number_ = number;
  seq_ = seq;
}

void MutationEvent::small_deletion(int32_t pos, int32_t number) {
  type_ = MutationEventType::SMALL_DELETION;
  pos_1_ = pos;
  number_ = number;
}

void MutationEvent::duplication(int32_t pos1, int32_t pos2, int32_t pos3) {
  type_ = MutationEventType::DUPLICATION;
  pos_1_ = pos1;
  pos_2_ = pos2;
  pos_3_ = pos3;
}

void MutationEvent::translocation(int32_t pos1, int32_t pos2, int32_t pos3,
                                    int32_t pos4, bool invert) {
  type_ = MutationEventType::TRANSLOCATION;
  pos_1_ = pos1;
  pos_2_ = pos2;
  pos_3_ = pos3;
  pos_4_ = pos4;
  invert_ = invert;
}

void MutationEvent::inversion(int32_t pos1, int32_t pos2) {
  type_ = MutationEventType::INVERSION;
  pos_1_ = pos1;
  pos_2_ = pos2;
}

void MutationEvent::deletion(int32_t pos1, int32_t pos2) {
  type_ = MutationEventType::DELETION;
  pos_1_ = pos1;
  pos_2_ = pos2;
}

MutationEvent::~MutationEvent() {
  if (seq_ != nullptr)
    delete [] seq_;
}
}
