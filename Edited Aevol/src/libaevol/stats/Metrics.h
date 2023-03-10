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


#ifndef AEVOL_INDIV_STATS_H_
#define AEVOL_INDIV_STATS_H_


// ============================================================================
//                                   Includes
// ============================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "Individual.h"
#include "GeneticUnit.h"


namespace aevol {

// ============================================================================
//                          Class declarations
// ============================================================================
class Individual;






class Metrics
{
  friend Individual;

 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  Metrics() = default; //< Default ctor
  Metrics(const Metrics&) = default; //< Copy ctor
  Metrics(Metrics&&) = default; //< Move ctor

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~Metrics() = default; //< Destructor

  // ==========================================================================
  //                                 Getters
  // ==========================================================================
  int32_t total_genome_size() const {
    return total_genome_size_;
  };
  int16_t nb_coding_RNAs() const {
    return nb_coding_RNAs_;
  };
  int16_t nb_non_coding_RNAs() const {
    return nb_non_coding_RNAs_;
  };
  int32_t overall_size_coding_RNAs() const {
    return overall_size_coding_RNAs_;
  };
  int32_t overall_size_non_coding_RNAs() const {
    return overall_size_non_coding_RNAs_;
  };
  int16_t nb_genes_activ() const {
    return nb_genes_activ_;
  };
  int16_t nb_genes_inhib() const {
    return nb_genes_inhib_;
  };
  int16_t nb_functional_genes() const {
    return nb_functional_genes_;
  };
  int16_t nb_non_functional_genes() const {
    return nb_non_functional_genes_;
  };
  int32_t overall_size_functional_genes() const {
    return overall_size_functional_genes_;
  };
  int32_t overall_size_non_functional_genes() const {
    return overall_size_non_functional_genes_;
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
  inline void Reset();
  inline void Accumulate(const GeneticUnit& gen_unit);





 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================

  // ==========================================================================
  //                               Attributes
  // ==========================================================================
  /// Sum of sizes of the genetic units
  int32_t total_genome_size_ = 0;
  /// Number of coding RNAs (at least one gene on RNA)
  int16_t nb_coding_RNAs_ = 0;
  /// Number of non-coding-RNAs
  int16_t nb_non_coding_RNAs_ = 0;
  /// Cumulated size of all coding RNAs
  int32_t overall_size_coding_RNAs_ = 0;
  /// Cumulated size of all non-coding RNAs
  int32_t overall_size_non_coding_RNAs_ = 0;
  /// Number of genes realizing a function
  int16_t nb_genes_activ_ = 0;
  /// Number of genes inhibitting a function
  int16_t nb_genes_inhib_ = 0;
  /// Number of functional genes
  int16_t nb_functional_genes_ = 0;
  /// Number of non-functional genes
  int16_t nb_non_functional_genes_ = 0;
  /// Cumulated size of all functional genes
  int32_t overall_size_functional_genes_ = 0;
  /// Cumulated size of all non-functional genes
  int32_t overall_size_non_functional_genes_ = 0;
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
void Metrics::Reset() {
  total_genome_size_                  = 0;
  nb_coding_RNAs_                     = 0;
  nb_non_coding_RNAs_                 = 0;
  overall_size_coding_RNAs_           = 0;
  overall_size_non_coding_RNAs_       = 0;
  nb_genes_activ_                     = 0;
  nb_genes_inhib_                     = 0;
  nb_functional_genes_                = 0;
  nb_non_functional_genes_            = 0;
  overall_size_functional_genes_      = 0;
  overall_size_non_functional_genes_  = 0;
}

void Metrics::Accumulate(const GeneticUnit& gen_unit) {
  total_genome_size_ += gen_unit.dna()->length();
  nb_coding_RNAs_ += gen_unit.nb_coding_RNAs();
  nb_non_coding_RNAs_ += gen_unit.nb_non_coding_RNAs();
  overall_size_coding_RNAs_ += gen_unit.overall_size_coding_RNAs();
  overall_size_non_coding_RNAs_ += gen_unit.overall_size_non_coding_RNAs();
  nb_genes_activ_ += gen_unit.nb_genes_activ();
  nb_genes_inhib_ += gen_unit.nb_genes_inhib();
  nb_functional_genes_ += gen_unit.nb_functional_genes();
  nb_non_functional_genes_ += gen_unit.nb_non_functional_genes();
  overall_size_functional_genes_ +=
      gen_unit.overall_size_functional_genes();
  overall_size_non_functional_genes_ +=
      gen_unit.overall_size_non_functional_genes();
};

} // namespace aevol

#endif // AEVOL_INDIV_STATS_H_
