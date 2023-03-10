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

#include <list>

// =================================================================
//                            Project Files
// =================================================================
#include "Rna.h"
#include "GeneticUnit.h"
#include "Individual.h"

namespace aevol {



//##############################################################################
//                                                                             #
//                                 Class Rna                                #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
Rna::Rna(GeneticUnit* gen_unit, const Rna &model)
{
  // Copy "trivial" attributes
  gen_unit_  = gen_unit;

  strand_             = model.strand_;
  pos_                = model.pos_;
  transcript_length_  = model.transcript_length_;
  basal_level_        = model.basal_level_;

  // Copy transcribed proteins
  // WARNING : Since this list do not "own" the proteins (they will not be deleted)
  //            proteins must NOT be CREATED here.

  // TODO : Not needed for the moment...
  // for (const auto& protein: model.transcribed_proteins_)
  //   transcribed_proteins_.push_back(protein);
}

Rna::Rna(GeneticUnit* gen_unit, Strand strand, int32_t pos, int8_t diff)
{
  gen_unit_  = gen_unit;
  strand_ = strand;
  pos_    = pos;

  transcript_length_  = -1;
  #ifdef BASE_2
#ifdef __POW_BASAL
  basal_level_ = (double) 1.0 / (1<<diff);
#else
  basal_level_        = 1 - (double)diff / (PROM_MAX_DIFF + 1);
#endif
  #elif BASE_4
  basal_level_ = BASAL_LEVELS[diff];
  #endif
}

/*
Rna::Rna(Rna* parent)
{
  gen_unit_           = parent->gen_unit_;
  strand_             = parent->strand_;
  pos_                = parent->pos_;
  transcript_length_  = parent->transcript_length_;
  basal_level_        = parent->basal_level_;
}
*/

// =================================================================
//                             Destructors
// =================================================================
Rna::~Rna()
{
}

// =================================================================
//                            Public Methods
// =================================================================
int32_t Rna::first_transcribed_pos() const
{
  if (strand_ == LEADING)
  {
    return Utils::mod(pos_ + PROM_SIZE, gen_unit_->dna()->length());
  }
  else
  {
    return Utils::mod(pos_ - PROM_SIZE, gen_unit_->dna()->length());
  }
}

int32_t Rna::last_transcribed_pos() const
{
  if (strand_ == LEADING)
  {
    return Utils::mod(pos_ +  PROM_SIZE + transcript_length_ - 1,
                       gen_unit_->dna()->length());
  }
  else
  {
    return Utils::mod(pos_ - (PROM_SIZE + transcript_length_ - 1),
                       gen_unit_->dna()->length());
  }
}

void Rna::copy_parent(const Rna* parent, bool env_will_changed) {

}

// =================================================================
//                           Protected Methods
// =================================================================
} // namespace aevol
