// ****************************************************************************
//
//          Aevol - An in silico experimental evolution platform
//
// ****************************************************************************
//
// Copyright: See the AUTHORS file provided with the package or <www.aevol.fr>
// Web: http://www.aevol.fr/
// E-mail: See <http://www.aevol.fr/contact/>
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons, Jonathan Rouzaud-Cornabas
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


#ifndef AEVOL_ABSTRACTFUZZY_7_H
#define AEVOL_ABSTRACTFUZZY_7_H


#include <list>

#include "macros.h"
#include "Point.h"
#include "Fuzzy.h"

namespace aevol {

class AbstractFuzzy_7
{
 public:
  // ==========================================================================
  //                               Constructors
  // ==========================================================================

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~AbstractFuzzy_7() {};

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  virtual void save(gzFile backup) const = 0;
  virtual void load(gzFile backup) = 0;
  virtual void reset() = 0;
  virtual void simplify() = 0;
  virtual void add_triangle(ProteinConcentration mean, ProteinConcentration width, ProteinConcentration height, bool verbose = false)  = 0;
  
  virtual void add(const AbstractFuzzy_7* f, bool verbose = false)  = 0;
  virtual void sub(const AbstractFuzzy_7* f, bool verbose = false) = 0;
  virtual void copy(const AbstractFuzzy_7* f, bool verbose = false) = 0;
  virtual void copy(const Fuzzy* f, ExpManager* exp_m = nullptr, bool verbose = false) = 0;

  virtual void add_point(ProteinConcentration x, ProteinConcentration y) = 0;

  /// `clipping_direction` is only used for `clip` function's keyword.
  enum clipping_direction: bool {min, max};
  virtual void clip(clipping_direction direction, ProteinConcentration bound) = 0;

  // ==========================================================================
  //                                 Getters
  // ==========================================================================

  virtual ProteinConcentration get_geometric_area(bool verbose = false) const = 0;
  virtual ProteinConcentration get_geometric_area(ProteinConcentration start_segment, ProteinConcentration end_segment) const = 0;

  virtual bool is_identical_to(const AbstractFuzzy_7& fs, ProteinConcentration tolerance) const = 0;

  virtual void print() const = 0;

  virtual void clear() = 0;

  // ==========================================================================
  //                                 Setters
  // ==========================================================================

  // ==========================================================================
  //                                Operators
  // ==========================================================================

 protected:
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================


  // ==========================================================================
  //                               Attributes
  // ==========================================================================

};

} // namespace aevol


#endif //AEVOL_AbstractFuzzy_7_H
