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

#ifndef AEVOL_GAUSSIAN_H_
#define AEVOL_GAUSSIAN_H_

// =================================================================
//                              Libraries
// =================================================================
#include <cstdlib>
#include <cmath>

#include <zlib.h>

namespace aevol {

// =================================================================
//                            Project Files
// =================================================================

// =================================================================
//                          Class declarations
// =================================================================

class Gaussian {
 public :

  // =================================================================
  //                             Constructors
  // =================================================================
  Gaussian(double height, double mean, double width) : height_{height}, mean_{mean}, width_{width} {}
  Gaussian(const Gaussian& model) : height_{model.height_}, mean_{model.mean_}, width_{model.width_} {}
  Gaussian(gzFile backup_file);

  // =================================================================
  //                             Destructor
  // =================================================================
  virtual ~Gaussian() {}

  // =================================================================
  //                              Accessors
  // =================================================================
  double height() const { return height_; }
  double mean() const { return mean_; }
  double width() const { return width_; }
  void   set_height(double height) { height_ = height; }
  void   set_mean(double mean) { mean_ = mean; }

  // =================================================================
  //                            Public Methods
  // =================================================================
#ifdef PHENOTYPIC_TARGET_TRIANGLE
double compute_y(double x) const {
    if (x < mean_ - width_) return 0;
    if (x > mean_ + width_) return 0;
    if (x >= mean_) return (height_ - (height_ * (x - mean_) / width_));
    return (height_ - (height_ * (mean_ - x) / width_));
  }
#else
    double compute_y(double x) const { return height_ * exp(-(x- mean_)*(x- mean_) / (2* width_ * width_)); }
#endif

  void save(gzFile backup_file) const;

  // =================================================================
  //                           Public Attributes
  // =================================================================

 protected :
  // =================================================================
  //                         Forbidden Constructors
  // =================================================================
  Gaussian() = delete;

  // =================================================================
  //                           Protected Methods
  // =================================================================

  // =================================================================
  //                          Protected Attributes
  // =================================================================
  double height_;
  double mean_;
  double width_; // In fact half-width to the inflexion points
};


// =====================================================================
//                               Constructors
// =====================================================================
inline Gaussian::Gaussian(gzFile backup_file) {
  gzread(backup_file, &height_,  sizeof(height_));
  gzread(backup_file, &mean_,    sizeof(mean_));
  gzread(backup_file, &width_,   sizeof(width_));
}

// =====================================================================
//                               Destructor
// =====================================================================

// =====================================================================
//                          Accessors' definitions
// =====================================================================

// =====================================================================
//                       functions' definition
// =====================================================================
inline void Gaussian::save(gzFile backup_file) const {
  gzwrite(backup_file, &height_, sizeof(height_));
  gzwrite(backup_file, &mean_, sizeof(mean_));
  gzwrite(backup_file, &width_, sizeof(width_));
}

} // namespace aevol

#endif // AEVOL_GAUSSIAN_H_
