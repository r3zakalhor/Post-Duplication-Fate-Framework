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

#ifndef AEVOL_PROMOTER_H
#define AEVOL_PROMOTER_H

#include <cstdint>
#include <map>

namespace aevol {
  class Rna_7;

class PromoterStruct {
 public:
  PromoterStruct(int32_t t_pos, int8_t t_error, bool lead) {
    pos                = t_pos;
    error              = t_error;
    leading_or_lagging = lead;
  }

  PromoterStruct(const PromoterStruct& clone) {
    pos                = clone.pos;
    error              = clone.error;
    leading_or_lagging = clone.leading_or_lagging;
  }

  PromoterStruct(PromoterStruct* clone) {
    pos                = clone->pos;
    error              = clone->error;
    leading_or_lagging = clone->leading_or_lagging;
  }

  int32_t pos  = -1;
  int8_t error = -1;
  bool leading_or_lagging; // TRUE = leading / FALSE = lagging


  bool to_compute_end_ = true;
  Rna_7* rna = nullptr;
};

class PromoterList {
  PromoterList() = default;

  std::map<int32_t, int> leading_prom_pos;
  std::map<int32_t, int> lagging_prom_pos;
};
}

#endif //AEVOL_PROMOTER_H
