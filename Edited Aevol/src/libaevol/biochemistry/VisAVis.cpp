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




// =================================================================
//                              Libraries
// =================================================================



// =================================================================
//                            Project Files
// =================================================================
#include "VisAVis.h"

namespace aevol {



// ############################################################################
//
//                              Class VisAVis
//
// ############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
VisAVis::VisAVis(const Dna* chrom_1, const Dna* chrom_2,
                 int32_t i_1, int32_t i_2, AlignmentSense sense /*= DIRECT*/) {
  chrom_1_  = chrom_1;
  chrom_2_  = chrom_2;
  i_1_      = i_1;
  i_2_      = i_2;
  sense_    = sense;
  score_    = 0;
}

VisAVis::VisAVis(const VisAVis & orig) {
  chrom_1_  = orig.chrom_1_;
  chrom_2_  = orig.chrom_2_;
  i_1_      = orig.i_1_;
  i_2_      = orig.i_2_;
  sense_    = orig.sense_;
  score_    = orig.score_;
}

//~ VisAVis::VisAVis(const VisAVis* orig)
//~ {
  //~ chrom_1_  = orig->chrom_1_;
  //~ chrom_2_  = orig->chrom_2_;
  //~ i_1_      = orig->i_1_;
  //~ i_2_      = orig->i_2_;
  //~ sense_    = orig->sense_;
//~ }

// =================================================================
//                             Destructors
// =================================================================
VisAVis::~VisAVis() {
}

// =================================================================
//                            Public Methods
// =================================================================

// =================================================================
//                           Protected Methods
// =================================================================
} // namespace aevol
