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

#ifndef AEVOL_POINT_H_
#define AEVOL_POINT_H_

#include <zlib.h>
#include <utility>

#include "Protein.h"

namespace aevol {

class Point {
public:
  ProteinConcentration x;
  mutable ProteinConcentration y;
  Point(ProteinConcentration x_, ProteinConcentration y_): x(x_), y(y_) {};
  Point() {};
  Point(Point* p): x(p->x), y(p->y) {};

  bool operator<(const Point & other) const {
      return (x <  other.x);
  }
};

Point readpoint(const gzFile backup_file);
void writepoint(const Point& p, gzFile backup_file);

} // namespace aevol
#endif // AEVOL_POINT_H_
