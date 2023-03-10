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

/// TODO: add unit tests
/// Why should there always be points (X_MIN, 0),(X_MAX, 0) ?
/// Many tests for ProteinConcentration-type equality ==. Should't we check mod ε?

#include "AeTime.h"
#include "ExpManager.h"
#include "Fuzzy.h"
#include "Point.h"
#include "Discrete_Double_Fuzzy.h"
#include "macros.h"
#include "FuzzyFactory_7.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>

using std::cout;
using std::endl;
using std::fabs;

namespace aevol {
void Discrete_Double_Fuzzy::simplify() {
  // NOT USEFUL
}

/// Add a triangle to the fuzzy set.
/// \param mean abscissa of its apex
/// \param width of the side opposite to the apex
/// \param height ordinate of the apex
void Discrete_Double_Fuzzy::add_triangle(ProteinConcentration mean,
                                         ProteinConcentration width,
                                         ProteinConcentration height,
                                         bool verbose) {
  assert(width > 0.0);
  assert(X_MIN <= mean and mean <= X_MAX);
  assert(W_MIN <= width);  // the maximum width depends on each individual
  // assert(MIN_H <= height and height <= MAX_H); Not necessarily because the concentration can be > 1
  #if defined(__INTEL_COMPILER)
  __declspec(align(points_));
  #elif defined(__INTEL_LLVM_COMPILER)
  void* vec_r = __builtin_assume_aligned(points_,64);
  #endif

  const ProteinConcentration threshold =
      1e-15;  // TODO: should it not be the machine epsilon?
              // if not, it should at least be a class constant

  if (fabs(width) < threshold or fabs(height) < threshold)
    return;

  double x0 = mean - width;
  double x1 = mean;
  double x2 = mean + width;

  int loop_A_start = (int)std::ceil(x0 * d_length_);
  loop_A_start     = loop_A_start < 0 ? 0 : loop_A_start;
  loop_A_start = loop_A_start > length_ ? length_
                                                      : loop_A_start;

  int loop_A_end = (int)std::ceil(x1 * d_length_);
  loop_A_end     = loop_A_end < 0 ? 0 : loop_A_end;
  loop_A_end =
      loop_A_end > length_ ? length_ : loop_A_end;

  for (int i = loop_A_start; i < loop_A_end; i++) {
    points_[i] += (((i / d_length_) - x0) / (x1 - x0)) * height;
  }

  // Compute the second equation of the triangle
  // Updating value between x1 and x2
  int loop_B_start = (int)std::ceil(x1 * d_length_);
  loop_B_start     = loop_B_start < 0 ? 0 : loop_B_start;
  loop_B_start = loop_B_start > length_ ? length_
                                                      : loop_B_start;

  int loop_B_end = (int)std::ceil(x2 * d_length_);
  if (loop_B_end > length_) {
    points_[length_ - 1] += height * ((x2 - 1.0) / (x2 - x1));
  }

  loop_B_end = loop_B_end < 0 ? 0 : loop_B_end;
  loop_B_end =
      loop_B_end > length_ ? length_ : loop_B_end;

  for (int i = loop_B_start; i < loop_B_end; i++) {
    points_[i] += height * ((x2 - (i / d_length_)) / (x2 - x1));
  }
}

/// Add a fuzzy set to the current one.
///
/// Should actually be called `operator+=()`.
///
/// Semantically speaking, we deal with fuzzy sets over the same
/// range. So adding two fuzzy sets sums up to adding the probability
/// functions.
void Discrete_Double_Fuzzy::add(const AbstractFuzzy_7* f,bool verbose) {
  Discrete_Double_Fuzzy* fs = (Discrete_Double_Fuzzy*)(f);
  #if defined(__INTEL_COMPILER)
  __declspec(align(fs->points_));
  __declspec(align(points_));
  #elif defined(__INTEL_LLVM_COMPILER)
  void* vec_r = __builtin_assume_aligned(fs->points_,64);
  vec_r = __builtin_assume_aligned(points_,64);
  #endif

  for (int i = 0; i < length_; i++) {
    points_[i] += fs->points_[i];
  }
}

/// Substract to the current fuzzy set.
///
/// TODO: Dumb version (?), to be completed.
void Discrete_Double_Fuzzy::sub(const AbstractFuzzy_7* f, bool verbose) {
  Discrete_Double_Fuzzy* fs = (Discrete_Double_Fuzzy*)(f);
  #if defined(__INTEL_COMPILER)
  __declspec(align(fs->points_));
  __declspec(align(points_));
    #elif defined(__INTEL_LLVM_COMPILER)
  void* vec_r = __builtin_assume_aligned(fs->points_,64);
  vec_r = __builtin_assume_aligned(points_,64);
  #endif

  for (int i = 0; i < length_; i++) {
    points_[i] -= fs->points_[i];
  }

}

void Discrete_Double_Fuzzy::copy(const AbstractFuzzy_7* f, bool verbose) {
  Discrete_Double_Fuzzy* fs = (Discrete_Double_Fuzzy*)(f);
  #if defined(__INTEL_COMPILER)
  __declspec(align(fs->points_));
  __declspec(align(points_));
    #elif defined(__INTEL_LLVM_COMPILER)
  void* vec_r = __builtin_assume_aligned(fs->points_,64);
  vec_r = __builtin_assume_aligned(points_,64);
  #endif

  for (int i = 0; i < length_; i++) {
    points_[i] = fs->points_[i];
  }
}

void Discrete_Double_Fuzzy::copy(const Fuzzy* f, ExpManager* exp_m, bool verbose) {
  for (int i = 0; i < length_; i++) {
    double tmp = ((Fuzzy*)f)
                     ->y(((double)i) / d_length_);
    points_[i] = tmp;
  }
}

ProteinConcentration
Discrete_Double_Fuzzy::get_geometric_area(bool verbose) const {
  return get_geometric_area(X_MIN, X_MAX);
}

ProteinConcentration Discrete_Double_Fuzzy::get_geometric_area(
    ProteinConcentration x_start, ProteinConcentration x_stop) const {
  // assert(invariant());
  // Precondition: X_MIN ≤ x_start < x_stop ≤ X_MAX
  assert(X_MIN <= x_start and x_start < x_stop and x_stop <= X_MAX);
  ProteinConcentration area = 0;
  #if defined(__INTEL_COMPILER)
  __declspec(align(points_));
  #elif defined(__INTEL_LLVM_COMPILER)
  void* vec_r = __builtin_assume_aligned(points_,64);
  #endif

  for (int fuzzy_idx = 0; fuzzy_idx < length_-1; fuzzy_idx++) {
    area += fabs(
        ((points_[fuzzy_idx] +
          points_[fuzzy_idx + 1]) /
         (d_length_ * 2)));
  }

  return area;
}

/// Probability function gets clipped either upwise ou downwise.
///
/// `pf` := min(`pf`, `upper_bound`)
///
///            X    above: removed               |
///           / \                                |
///          /   \               X      bound    |
/// --------o-----o-------------o-o--------      |
///        /       \   X       /   \             |
///       X         \ / \     /     \            |
///                  X   \   /       X           |
///                       \ /                    |
///      underneath: kept  X                     |

/// `pf` := max(`pf`, `lower_bound`)
void Discrete_Double_Fuzzy::clip(clipping_direction direction,
                                 ProteinConcentration bound) {
  // assert(invariant());
  #if defined(__INTEL_COMPILER)
  __declspec(align(points_));
  #elif defined(__INTEL_LLVM_COMPILER)
  void* vec_r = __builtin_assume_aligned(points_,64);
  #endif

  for (int i = 0; i < length_; i++) {
    if ((direction == AbstractFuzzy_7::clipping_direction::min and
         points_[i] < bound) or
        (direction == AbstractFuzzy_7::clipping_direction::max and
         points_[i] > bound))
      points_[i] = bound;
  }

  // assert(invariant());
}

bool Discrete_Double_Fuzzy::is_identical_to(
    const AbstractFuzzy_7& f, ProteinConcentration tolerance) const {
  assert(false);
}

void Discrete_Double_Fuzzy::save(gzFile backup_file) const { assert(false); }

void Discrete_Double_Fuzzy::load(gzFile backup_file) { assert(false); }

/// Set all points ordinate to 0
///
// TODO <david.parsons@inria.fr> Not sure if it's useful.
void Discrete_Double_Fuzzy::reset() {
      // set are ordered data structure and then only generate const iterator
    // But y attribute of Point is marked as mutable, i.e. still mutable even
    // through a const reference
  #if defined(__INTEL_COMPILER)
  __declspec(align(points_));
  #elif defined(__INTEL_LLVM_COMPILER)
  void* vec_r = __builtin_assume_aligned(points_,64);
  #endif

  for (int i = 0; i < length_; i++)
    points_[i] = 0.0;
}

void Discrete_Double_Fuzzy::clear() { reset(); }

void Discrete_Double_Fuzzy::add_point(ProteinConcentration x,
                                      ProteinConcentration y) {
  
  if ((int)std::ceil(x*d_length_) < d_length_) {
    points_[(int)std::ceil(x*d_length_)] = y;
    printf("Add point %lf (%d) :: %lf\n",x,(int)std::ceil(x*d_length_),y);
  }
}

void Discrete_Double_Fuzzy::add_point(int32_t x,
                                      ProteinConcentration y) {
    points_[x] = y;  
}

void Discrete_Double_Fuzzy::print() const {
  for (int i = 0; i < length_; i++)
    printf("[%f : %f]\n", i / d_length_, points_[i]);
  printf("\n");
}
}  // namespace aevol
