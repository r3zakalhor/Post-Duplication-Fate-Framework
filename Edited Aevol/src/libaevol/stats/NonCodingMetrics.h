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


#ifndef AEVOL_NON_CODING_STATS_H_
#define AEVOL_NON_CODING_STATS_H_


// ============================================================================
//                                   Includes
// ============================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "GeneticUnit.h"


namespace aevol {

// ============================================================================
//                          Class declarations
// ============================================================================





// TODO <david.parsons@inria.fr> Not used ?
class NonCodingMetrics
{
  friend Individual;

 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  NonCodingMetrics() = default; //< Default ctor
  NonCodingMetrics(const NonCodingMetrics&) = default; //< Copy ctor
  NonCodingMetrics(NonCodingMetrics&&) = default; //< Move ctor

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~NonCodingMetrics() = default; //< Destructor

  // ==========================================================================
  //                                 Getters
  // ==========================================================================
  int32_t nb_bases_in_0_CDS() const {
    return nb_bases_in_0_CDS_;
  };
  int32_t nb_bases_in_0_functional_CDS() const {
    return nb_bases_in_0_functional_CDS_;
  };
  int32_t nb_bases_in_0_non_functional_CDS() const {
    return nb_bases_in_0_non_functional_CDS_;
  };
  int32_t nb_bases_in_0_RNA() const {
    return nb_bases_in_0_RNA_;
  };
  int32_t nb_bases_in_0_coding_RNA() const {
    return nb_bases_in_0_coding_RNA_;
  };
  int32_t nb_bases_in_0_non_coding_RNA() const {
    return nb_bases_in_0_non_coding_RNA_;
  };
  int32_t nb_bases_in_neutral_regions() const {
    return nb_bases_in_neutral_regions_;
  };
  int32_t nb_neutral_regions() const {
    return nb_neutral_regions_;
  };

  // ==========================================================================
  //                                 Setters
  // ==========================================================================

  // ==========================================================================
  //                                Operators
  // ==========================================================================

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  inline void Accumulate(const GeneticUnit& gen_unit);





 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================

  // ==========================================================================
  //                               Attributes
  // ==========================================================================
  /// Number of bases that are not included in any gene
  int32_t nb_bases_in_0_CDS_;
  /// Number of bases that are not included in any functional gene
  int32_t nb_bases_in_0_functional_CDS_;
  /// Number of bases that are not included in any degenerated gene
  int32_t nb_bases_in_0_non_functional_CDS_;
  /// Number of bases that are not included in any RNA
  int32_t nb_bases_in_0_RNA_;
  /// Number of bases that are not included in any coding RNA
  /// (RNAs containing at least one CDS)
  int32_t nb_bases_in_0_coding_RNA_;
  /// Number of bases that are not included in any non coding RNA
  int32_t nb_bases_in_0_non_coding_RNA_;
  /// Number of bases that are in a neutral region
  /// A base is considered neutral when neither itself NOR its corresponding base on the other
  /// strand belongs to a coding promoter->terminator region (both included)
  int32_t nb_bases_in_neutral_regions_;
  /// Number of neutral regions
  int32_t nb_neutral_regions_;
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
void NonCodingMetrics::Accumulate(const GeneticUnit& gen_unit) {
  nb_bases_in_0_CDS_ += gen_unit.nb_bases_in_0_CDS();
  nb_bases_in_0_functional_CDS_ += gen_unit.nb_bases_in_0_functional_CDS();
  nb_bases_in_0_non_functional_CDS_ +=
      gen_unit.nb_bases_in_0_non_functional_CDS();
  nb_bases_in_0_RNA_ += gen_unit.nb_bases_in_0_RNA();
  nb_bases_in_0_coding_RNA_ += gen_unit.nb_bases_in_0_coding_RNA();
  nb_bases_in_0_non_coding_RNA_ += gen_unit.nb_bases_in_0_non_coding_RNA();
  nb_bases_in_neutral_regions_ += gen_unit.nb_bases_in_neutral_regions();
  nb_neutral_regions_ += gen_unit.nb_neutral_regions();
}
} // namespace aevol

#endif // AEVOL_NON_CODING_STATS_H_
