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
#include "PhenotypicTarget.h"
#include "FuzzyFactory.h"

#include <cstring>


namespace aevol {


//##############################################################################
//                                                                             #
//                           Class PhenotypicTarget                            #
//                                                                             #
//##############################################################################

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================
PhenotypicTarget::PhenotypicTarget() {
  fuzzy_ = FuzzyFactory::fuzzyFactory->create_fuzzy();
  nb_segments_     = 1;
  segments_        = new PhenotypicSegment * [1];
  segments_[0]     = new PhenotypicSegment(X_MIN, X_MAX, METABOLISM);
  area_by_feature_ = new double [NB_FEATURES];
}

PhenotypicTarget::PhenotypicTarget(const PhenotypicTarget& rhs) {
  fuzzy_ = FuzzyFactory::fuzzyFactory->create_fuzzy(*(rhs.fuzzy()));
  nb_segments_     = rhs.nb_segments_;
  segments_        = new PhenotypicSegment * [nb_segments_];
  for (int8_t i = 0 ; i < nb_segments_ ; i++)
    segments_[i] = new PhenotypicSegment(*(rhs.segments_[i]));
  area_by_feature_ = new double [NB_FEATURES];
  memcpy(area_by_feature_,
         rhs.area_by_feature_,
         NB_FEATURES * sizeof(*area_by_feature_));
}

// ============================================================================
//                                 Destructor
// ============================================================================
PhenotypicTarget::~PhenotypicTarget() {
  if (segments_ != NULL) {
    for (int8_t i = 0 ; i < nb_segments_ ; i++)
      delete segments_[i];
    delete [] segments_;
  }
  delete [] area_by_feature_;
  delete fuzzy_;
}

// ============================================================================
//                                   Methods
// ============================================================================
void PhenotypicTarget::set_segmentation(int16_t nb_segments,
                                        double *boundaries,
                                        PhenotypicFeature *features,
                                        bool separate_segments) {
  // Delete the data to be replaced
  for (int16_t i = 0 ; i < nb_segments_ ; i++)
    delete segments_[i];
  delete[] segments_;

  // Now replace with the new data
  nb_segments_  = nb_segments;
  segments_     = new PhenotypicSegment * [nb_segments_];

  for (int16_t i = 0 ; i < nb_segments_; i++)
    segments_[i] = new PhenotypicSegment(boundaries[i], boundaries[i+1], features[i]);

  // TODO <dpa>: Manage separate_segments
}

void PhenotypicTarget::ComputeArea() {
  for (int8_t i = 0 ; i < NB_FEATURES ; i++)
    area_by_feature_[i] = 0.0;

  // TODO <dpa>: We should take into account that we compute the areas in a specific order (from the leftmost segment, rightwards)
  //   => We shouldn't parse the whole list of points on the left of the segment we are considering (we have
  //      already been through them!)
  for (int8_t i = 0 ; i < nb_segments_ ; i++) {
    area_by_feature_[segments_[i]->feature] +=
      fuzzy_->get_geometric_area(segments_[i]->start, segments_[i]->stop);
  }
}

void PhenotypicTarget::SaveSegmentation(gzFile backup_file) const {
  // --------------------------------------------------------------------------
  //  Write x-axis segmentation
  gzwrite(backup_file, &nb_segments_, sizeof(nb_segments_));

  for (int8_t i = 0 ; i < nb_segments_; i++) // TODO <david.parsons@inria.fr> suppress warning
    segments_[i]->save(backup_file);
}

void PhenotypicTarget::LoadSegmentation(gzFile backup_file) {
  // Delete obsolete segmentation data
  for (int8_t i = 0 ; i < nb_segments_ ; i++)
    delete segments_[i];
  delete [] segments_;

  // Replace by data from the backup
  gzread(backup_file, &nb_segments_, sizeof(nb_segments_));
  segments_ = new PhenotypicSegment * [nb_segments_];
  for (int8_t i = 0 ; i < nb_segments_; i++)
    segments_[i] = new PhenotypicSegment(backup_file);
}

// ============================================================================
//                            Non inline accessors
// ============================================================================
} // namespace aevol
