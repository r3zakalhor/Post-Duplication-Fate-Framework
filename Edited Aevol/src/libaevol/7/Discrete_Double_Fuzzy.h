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

#ifndef AEVOL_DISCRETE_DOUBLE_FUZZY_H
#define AEVOL_DISCRETE_DOUBLE_FUZZY_H

#include <set>

#include "macros.h"
#include "Point.h"
#include "AbstractFuzzy_7.h"
#include "FuzzyFactory_7.h"
#include "Fuzzy.h"
namespace aevol {

/// Triangular fuzzy sets.
///
/// This class provides management tools for "fuzzy sets" abscissa are
/// bound between X_MIN and X_MAX (defined in macros.h) A "fuzzy
/// set" should always have at least two points of abscissa X_MIN and
/// X_MAX.
///
/// A fuzzy set holds elements in a set X together with a probability
/// function X → [0,1] which tells how likely element x ∈ X beholds to
/// the fuzzy set. With these triangular fuzzy sets, the probability
/// function is a finite sum of isosceles triangles.
///
/// The current class models fuzzy sets over range X = [X_MIN; X_MAX]
/// by representing the probability function as the list of singular
/// points on its graph.
///
/// \verbatim
/// \\code{.unparsed}
///           ^
///       y2  +...              X                                            ....
///           |...             / \                                           ....
///       y6  +...            /   \               X                          ....
///           |...           /     \             / \                         ....
///       y4  +...          /       \   X       /   \                        ....
///           |...         /         \ / \     /     \                       ....
///      y3,y9+...        /           X   \   /       \             X        ....
///           |...       /                 \ /         \           / \       ....
///       y5  +...      /                   X           \         /   \      ....
///           |...     /                                 \       /     \     ....
///       0   +--X----X---------|-----|-|---|-----|-------X-----X-------X----X--->
///            X_MIN x1        x2    x3 x4 x5    x6      x7    x8  x9  x10 X_MAX
/// \\endcode
/// \endverbatim
/// fs.points_ would hold the list {(X_MIN,0),(x1,y1),...,(x10,y10)(X_MAX,0)}
///
/// \invariant{`points_.size()` ≥ 2}
/// \invariant{`points_.begin()->x == X_MIN`}
/// \invariant{`prev(points_.end())->x == X_MAX`}
/// \invariant{`is_increasing()`}
class Discrete_Double_Fuzzy : public AbstractFuzzy_7
{
 public:
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  Discrete_Double_Fuzzy(int vector_size) {
      length_ = vector_size;
      d_length_ = (double) length_;
      posix_memalign((void**)&points_, 64, vector_size*sizeof(double));//new double[vector_size];
      for (int i = 0; i < vector_size; i++) {
        points_[i] = 0;
      }
  };
  Discrete_Double_Fuzzy(int vector_size,const Discrete_Double_Fuzzy& f){
      length_ = vector_size;
      d_length_ = (double) length_;
      //points_ = new double[vector_size];
      posix_memalign((void**)&points_, 64, vector_size*sizeof(double));
      for (int i = 0; i < vector_size; i++) {
        points_[i] = f.points_[i];
      }  
  };

  Discrete_Double_Fuzzy(const Fuzzy& f){
      assert(false);
  };

  Discrete_Double_Fuzzy(const gzFile backup) { load(backup); };

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~Discrete_Double_Fuzzy() { free(points_);};

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  void save(gzFile backup) const;
  void load(gzFile backup);
  void reset();
  void simplify();
  void add_triangle(ProteinConcentration mean, ProteinConcentration width, ProteinConcentration height, bool verbose = false);
  void add(const AbstractFuzzy_7* f, bool verbose = false);
  void sub(const AbstractFuzzy_7* f, bool verbose = false);
  void copy(const AbstractFuzzy_7* f, bool verbose = false);
  void copy(const Fuzzy* f, ExpManager* exp_m = nullptr, bool verbose = false);
  void add_point(ProteinConcentration x, ProteinConcentration y);
  void add_point(int32_t x, ProteinConcentration y);

  void clip(clipping_direction direction, ProteinConcentration bound);
  // ==========================================================================
  //                                 Getters
  // ==========================================================================
  ProteinConcentration get_geometric_area(bool verbose = false) const;
  ProteinConcentration get_geometric_area(ProteinConcentration start_segment, ProteinConcentration end_segment) const;
  bool is_identical_to(const AbstractFuzzy_7& fs, ProteinConcentration tolerance) const;
  void print() const;
  void clear();
  double* points() { return points_; };
  int length() { return length_; };
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
  int length_;
  double d_length_;
  double* points_;
};
} // namespace aevol
#endif // AEVOL_FUZZY_H_
