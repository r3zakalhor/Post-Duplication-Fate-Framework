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


#ifndef AEVOL_RNA_H_
#define AEVOL_RNA_H_


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>




// =================================================================
//                            Project Files
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Dna.h"
#include "Utils.h"
#include "Protein.h"

namespace aevol {



// =================================================================
//                          Class declarations
// =================================================================
class Individual;
class GeneticUnit;

class Rna
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    Rna() = delete;
    Rna(const GeneticUnit&) = delete;
    Rna(GeneticUnit* gen_unit, const Rna &model);
    Rna(GeneticUnit* gen_unit, Strand strand, int32_t index, int8_t diff);
    //Rna(Rna* parent);

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~Rna();

    // =================================================================
    //                              Accessors
    // =================================================================

    // <DEBUG>
    void check(GeneticUnit* gen_unit) { assert(gen_unit == gen_unit_); };
    //~ void* indiv() const { return (void*)indiv_; };
    // </DEBUG>

    inline const GeneticUnit * genetic_unit() const;
    inline void set_genetic_unit(const GeneticUnit*  gen_unit);
    inline Strand strand() const;
    inline void       set_strand(Strand strand);
    inline int32_t    promoter_pos() const;
    inline void       set_promoter_pos(int32_t pos);
    inline ProteinConcentration     basal_level() const;
    inline int32_t    transcript_length() const; // The promoter is NOT transcribed.
    inline void       set_transcript_length(int32_t length);
    inline bool       is_coding() const;

    inline const std::list<Protein *>& transcribed_proteins() const;
    inline void clear_transcribed_proteins() { transcribed_proteins_.clear(); };

    // =================================================================
    //                            Public Methods
    // =================================================================
    int32_t first_transcribed_pos() const;   // The promoter is NOT transcribed.
    int32_t last_transcribed_pos() const;    // The terminator is transcribed.
    inline void add_transcribed_protein(Protein * prot);
    inline void shift_position(int32_t delta_pos, int32_t genome_length);

    void copy_parent(const Rna* parent, bool env_will_changed);

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================

    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    const GeneticUnit*  gen_unit_;
    Strand strand_;
    int32_t pos_; // Index of the promoter on the genome.
                  // The promoter itself is NOT transcribed
                  // The terminator is transcribed.
    int32_t transcript_length_;
    ProteinConcentration basal_level_;

    // Access list to the proteins transcribed by this rna
    std::list<Protein*> transcribed_proteins_;
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================
inline const GeneticUnit*Rna::genetic_unit() const
{
  return gen_unit_;
}

inline void Rna::set_genetic_unit(const GeneticUnit*  gen_unit)
{
  gen_unit_ = gen_unit;
}

inline Strand Rna::strand() const
{
  return strand_;
}

inline void Rna::set_strand(Strand strand)
{
  strand_ = strand;
}

void Rna::set_promoter_pos(int32_t pos)
{
  pos_ = pos;
}

inline int32_t Rna::promoter_pos() const
{
  return pos_;
}

inline ProteinConcentration Rna::basal_level() const
{
  return basal_level_;
}

inline int32_t Rna::transcript_length() const
{
  return transcript_length_;
}

inline void Rna::set_transcript_length(int32_t transcript_length)
{
  transcript_length_ = transcript_length;
}

inline const std::list<Protein *>&Rna::transcribed_proteins() const {
  return transcribed_proteins_;
}

inline bool Rna::is_coding() const
{
  return (not transcribed_proteins_.empty());
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
void Rna::add_transcribed_protein(Protein * prot)
{
  transcribed_proteins_.push_back(prot);
}

void Rna::shift_position(int32_t delta_pos, int32_t genome_length)
{
  pos_ = Utils::mod(pos_ + delta_pos, genome_length);
}

} // namespace aevol
#endif // AEVOL_RNA_H_
