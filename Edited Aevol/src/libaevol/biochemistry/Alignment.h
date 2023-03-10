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


#ifndef AEVOL_ALIGNMENT_H_
#define AEVOL_ALIGNMENT_H_


// =================================================================
//                              Includes
// =================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "VisAVis.h"

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================
class Dna;






class Alignment {
 public :

  // =================================================================
  //                             Constructors
  // =================================================================
  Alignment() = delete;
  Alignment(const Alignment&) = delete;

  // =================================================================
  //                             Destructors
  // =================================================================

  // =================================================================
  //                              Accessors
  // =================================================================

  // =================================================================
  //                            Public Methods
  // =================================================================
  static VisAVis* search_alignment_direct(const Dna * chrom1,
                                          const int32_t seed1,
                                          const Dna* chrom2,
                                          const int32_t seed2,
                                          int16_t needed_score);
  static VisAVis* search_alignment_indirect(const Dna* chrom1,
                                            const int32_t seed1,
                                            const Dna* chrom2,
                                            const int32_t seed2,
                                            int16_t needed_score);

  // =================================================================
  //                           Public Attributes
  // =================================================================

  static bool with_alignments;

  static AlignmentFunctionShape align_fun_shape;

  static double  align_sigm_lambda;
  static int16_t align_sigm_mean;
  static int16_t align_lin_min;
  static int16_t align_lin_max;

  // Maximum shift of one seq on the other
  static int16_t align_max_shift;
  // Work zone half length
  static int16_t align_w_zone_h_len;
  // Corresponding residues match bonus
  static int16_t align_match_bonus;
  // Corresponding residues mismatch cost
  static int16_t align_mismatch_cost;





 protected :

  // =================================================================
  //                           Protected Methods
  // =================================================================

  // =================================================================
  //                          Protected Attributes
  // =================================================================
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol

#endif // AEVOL_ALIGNMENT_H_
