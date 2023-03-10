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


#ifndef AEVOL_PHENOTYPIC_TARGET_H_
#define AEVOL_PHENOTYPIC_TARGET_H_


// ============================================================================
//                                   Includes
// ============================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "Fuzzy.h"
#include "PhenotypicSegment.h"
#include "ae_enums.h"


namespace aevol {

// ============================================================================
//                          Class declarations
// ============================================================================






class PhenotypicTarget
{
  friend class PhenotypicTargetHandler;

 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  PhenotypicTarget(); //< Default ctor
  PhenotypicTarget(const PhenotypicTarget&); //< Copy ctor
  PhenotypicTarget(PhenotypicTarget&&) = delete; //< Move ctor

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~PhenotypicTarget(); //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  void ComputeArea();

  void SaveSegmentation(gzFile backup_file) const;
  void LoadSegmentation(gzFile backup_file);

  // ==========================================================================
  //                                 Getters
  // ==========================================================================
  double area_by_feature(int8_t feature) const {
    return area_by_feature_[feature];
  }
  int8_t nb_segments() const {
    return nb_segments_;
  }
  PhenotypicSegment ** segments() const {
    return segments_;
  }

  AbstractFuzzy* fuzzy() const {
    return fuzzy_;
  }
// ==========================================================================
  //                                 Setters
  // ==========================================================================
  void set_segmentation(int16_t nb_segments,
                        double* boundaries,
                        PhenotypicFeature * features,
                        bool separate_segments);




 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================

  // ==========================================================================
  //                               Attributes
  // ==========================================================================
  /// Number segments.
  int8_t nb_segments_;
  /// Ordered array of segments.
  /// Each PhenotypicSegment knows its boundaries and corresponding feature.
  /// When the phenotypic target is not segmented, this array contains a single
  /// segment with feature METABOLIC and boundaries MIN_X and MAX_X
  PhenotypicSegment ** segments_;
  /// Geometric area of each feature
  double* area_by_feature_;

  AbstractFuzzy* fuzzy_;
};


// ============================================================================
//                           Getters' definitions
// ============================================================================

// ============================================================================
//                           Setters' definitions
// ============================================================================

// ============================================================================
//                          Operators' definitions
// ============================================================================

// ============================================================================
//                       Inline functions' definition
// ============================================================================

} // namespace aevol

#endif // AEVOL_PHENOTYPIC_TARGET_H_
