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


#ifndef AEVOL_PHENOTYPIC_SEGMENT_H_
#define AEVOL_PHENOTYPIC_SEGMENT_H_


// =================================================================
//                              Includes
// =================================================================
#include <cinttypes>

#include <zlib.h>

#include "macros.h"
#include "ae_enums.h"


namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================






class PhenotypicSegment
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    PhenotypicSegment() = delete;
    inline PhenotypicSegment(double start, double stop, PhenotypicFeature feature);
    inline PhenotypicSegment(const PhenotypicSegment & source);
    inline PhenotypicSegment(gzFile backup_file);

    // =================================================================
    //                             Destructors
    // =================================================================
    inline virtual ~PhenotypicSegment();

    // =================================================================
    //                              Accessors
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    inline void save(gzFile backup_file) const;
    inline void load(gzFile backup_file);

    // =================================================================
    //                           Public Attributes
    // =================================================================
    double start;
    double stop;
    PhenotypicFeature feature;





  protected :
    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
};




//##############################################################################
//                                                                             #
//                             Class PhenotypicSegment                            #
//                                                                             #
//##############################################################################

// =================================================================
//                             Constructors
// =================================================================
//~ inline PhenotypicSegment::PhenotypicSegment()
//~ {
  //~ start   = X_MIN;
  //~ stop    = X_MAX;
  //~ feature = NEUTRAL;
//~ }

inline PhenotypicSegment::PhenotypicSegment(double start, double stop, PhenotypicFeature feature)
{
  this->start   = start;
  this->stop    = stop;
  this->feature = feature;
}

inline PhenotypicSegment::PhenotypicSegment(const PhenotypicSegment & source)
{
  this->start   = source.start;
  this->stop    = source.stop;
  this->feature = source.feature;
}

inline PhenotypicSegment::PhenotypicSegment(gzFile backup_file)
{
  load(backup_file);
}

// =================================================================
//                             Destructors
// =================================================================
inline PhenotypicSegment::~PhenotypicSegment()
{
}

// =====================================================================
//                          Accessors definitions
// =====================================================================

// =================================================================
//                            Public Methods
// =================================================================
inline void PhenotypicSegment::save(gzFile backup_file) const
{
  gzwrite(backup_file, &start, sizeof(start));
  gzwrite(backup_file, &stop,  sizeof(stop));
  int8_t tmp_feature = feature;
  gzwrite(backup_file, &tmp_feature, sizeof(tmp_feature));
}

inline void PhenotypicSegment::load(gzFile backup_file)
{
  gzread(backup_file, &start,  sizeof(start));
  gzread(backup_file, &stop,   sizeof(stop));
  int8_t tmp_feature;
  gzread(backup_file, &tmp_feature, sizeof(tmp_feature));
  feature = (PhenotypicFeature) tmp_feature;
}

// =================================================================
//                           Protected Methods
// =================================================================


} // namespace aevol
#endif // AEVOL_PHENOTYPIC_SEGMENT_H_
