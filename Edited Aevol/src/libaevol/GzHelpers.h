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

#ifndef AEVOL_GZHELPERS_H_
#define AEVOL_GZHELPERS_H_

#include <zlib.h>

// ===========================================================================
//                                 Public Methods
// ===========================================================================
enum class GzAction {READ, WRITE};
// begin variadic template gzwrite(...)

// Base case for the next template.
// Useless by itself. Not intented to be called directly.
void gz(GzAction action, gzFile file) {
  return;
}

/// Read/write variables to gzip file
/// \param `action` tells whether to read or write
/// \param `file` an open gzip file
/// \param `field` the field to be written
/// \param `fields_list` the remaining fields
/// The function is simply called like: `gz(action, file, x, y, z, t)`.
/// Warning: as it is currently written, this template overrides any
/// call to gzwrite, which causes writing the length of the field to the file.
template<typename Field, typename... Args>
void gz(GzAction action, gzFile file, Field& field, Args&&... fields_list) {
  // This switch is unfortunate especially since the alternatives are pretty much the same.
  // But there is a subtle difference that prevented me from a trivial factorization:
  // gzwrite takes a _const_ void pointer as second argument.
  switch (action) {
    case GzAction::READ:
      ::gzread(file, &field, sizeof(field));
      break;
    case GzAction::WRITE:
      ::gzwrite(file, &field, sizeof(field));
      break;
  }
  gz(action, file, fields_list...);
}
// end variadic template gzwrite

template<typename... Args>
void gzwrite(Args... args) {
  gz(GzAction::WRITE, args...);
}

template<typename... Args>
void gzread(Args&&... args) {
  gz(GzAction::READ, args...);
}

#endif //AEVOL_GZHELPERS_H_
