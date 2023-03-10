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


#ifndef AEVOL_VIS_A_VIS_H_
#define AEVOL_VIS_A_VIS_H_


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


namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================







class VisAVis
{
  friend class Alignment;

  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    VisAVis() = default;
    VisAVis(const Dna * chrom_1, const Dna * chrom_2,
            int32_t i_1, int32_t i_2, AlignmentSense sense = DIRECT);
    VisAVis(const VisAVis & orig);

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~VisAVis();

    // =================================================================
    //                              Accessors
    // =================================================================
    inline const Dna*  chrom_1() const;
    inline const Dna*  chrom_2() const;
    inline int32_t        i_1() const;
    inline int32_t        i_2() const;
    inline int16_t        score() const;
    inline AlignmentSense sense() const;

    // =================================================================
    //                              Operators
    // =================================================================
    inline bool operator <  (VisAVis &cmp);
    inline bool operator <= (VisAVis &cmp);
    inline bool operator >  (VisAVis &cmp);
    inline bool operator >= (VisAVis &cmp);

    // =================================================================
    //                            Public Methods
    // =================================================================
    inline bool match();
    inline void step_fwd();
    inline void step_back();
    inline void add(int common_inc);
    inline void add(int inc_1, int inc_2);
    inline void sub(int common_inc);
    inline void sub(int inc_1, int inc_2);
    inline void swap();

    inline void copy(VisAVis * source);
    inline void check_indices();

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :
    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    const Dna* chrom_1_ = nullptr;
    const Dna* chrom_2_ = nullptr;
    int32_t i_1_ = -1; //< Index on chrom_1
    int32_t i_2_ = -1; //< Index on chrom_2
    int16_t score_ = -1;
    // Sense (DIRECT or INDIRECT) of the vis_a_vis (alignement)
    AlignmentSense sense_ = DIRECT;
    // Say we have the following sequences :
    //    0 1 2 3 4 5 6 7 8 9             0 1 2 3 4 5 6 7 8 9
    //    |a|b|c|d|e|f|g|h|i|j|           |a|b|c|d|e|f|g|h|i|j|
    //
    // The DIRECT vis_a_vis between i_1_ = 3 and i_2_ = 7 is 'd' with 'h' (caracters at indices 3 and 7 resp.).
    //   The corresponding alignment would be "defgh" with "hijkl"
    //
    // WARNING! The INDIRECT vis_a_vis between the same i_1_ = 3 and i_2 = 7 is 'd' with 'g' (and not 'h'!).
    // This is because we are reading backwards (towards the left). Directly left to index 7 is 'g' which corresponds to index 6.
    //   The corresponding alignment would hence be "defgh" with "!g!f!e!d!c" ("!x" means "complementary of x")
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================
inline const Dna *VisAVis::chrom_1() const
{
  return chrom_1_;
}

inline const Dna *VisAVis::chrom_2() const
{
  return chrom_2_;
}

inline int32_t VisAVis::i_1() const
{
  return i_1_;
}

inline int32_t VisAVis::i_2() const
{
  return i_2_;
}

inline int16_t VisAVis::score() const
{
  return score_;
}

inline AlignmentSense VisAVis::sense() const
{
  return sense_;
}



// =====================================================================
//                          Operators' definitions
// =====================================================================
inline bool VisAVis::operator < (VisAVis &cmp)
{
  return (i_1_ < cmp.i_1_);
}

inline bool VisAVis::operator <= (VisAVis &cmp)
{
  return (i_1_ <= cmp.i_1_);
}

inline bool VisAVis::operator > (VisAVis &cmp)
{
  return (i_1_ > cmp.i_1_);
}

inline bool VisAVis::operator >= (VisAVis &cmp)
{
  return (i_1_ >= cmp.i_1_);
}


// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline bool VisAVis::match()
{
  if (sense_ == DIRECT)
  {
    return (chrom_1_->data()[Utils::mod(i_1_, chrom_1_->length())] ==
        chrom_2_->data()[Utils::mod(i_2_, chrom_2_->length())]);
  }
  else // (sense_ == INDIRECT)
  {
    // Note that we are reading the sequence backwards, The nucleotide corresponding to a breakpoint at point <i>
    // is hence stored at index <i-1>
    //    a b c d e f g h i j
    //    |_|_|_|_|_|_|_|_|_|_|
    //    | | | | | | | | | | |
    //      9 8 7 6 5 4 3 2 1 0
    //
    // The breakpoint F-5 puts into a vis_a_vis the nucleotide at index F on seq1 and that at index 4 (not 5!!!) on seq2
    return (chrom_1_->data()[Utils::mod(i_1_, chrom_1_->length())] !=
        chrom_2_->data()[Utils::mod(i_2_-1, chrom_2_->length())]);
  }
}

inline void VisAVis::step_fwd()
{
  if (sense_ == DIRECT)
  {
    i_1_++;
    i_2_++;
  }
  else // (sense_ == INDIRECT)
  {
    i_1_++;
    i_2_--;
  }
}

inline void VisAVis::step_back()
{
  if (sense_ == DIRECT)
  {
    i_1_--;
    i_2_--;
  }
  else // (sense_ == INDIRECT)
  {
    i_1_--;
    i_2_++;
  }
}

inline void VisAVis::add(int common_inc)
{
  if (sense_ == DIRECT)
  {
    i_1_ += common_inc;
    i_2_ += common_inc;
  }
  else // (sense_ == INDIRECT)
  {
    i_1_ += common_inc;
    i_2_ -= common_inc;
  }
}

inline void VisAVis::add(int inc_1, int inc_2)
{
  if (sense_ == DIRECT)
  {
    i_1_ += inc_1;
    i_2_ += inc_2;
  }
  else // (sense_ == INDIRECT)
  {
    i_1_ += inc_1;
    i_2_ -= inc_2;
  }
}

inline void VisAVis::sub(int common_inc)
{
  if (sense_ == DIRECT)
  {
    i_1_ -= common_inc;
    i_2_ -= common_inc;
  }
  else // (sense_ == INDIRECT)
  {
    i_1_ -= common_inc;
    i_2_ += common_inc;
  }
}

inline void VisAVis::sub(int inc_1, int inc_2)
{
  if (sense_ == DIRECT)
  {
    i_1_ -= inc_1;
    i_2_ -= inc_2;
  }
  else // (sense_ == INDIRECT)
  {
    i_1_ -= inc_1;
    i_2_ += inc_2;
  }
}

inline void VisAVis::swap()
{
  const Dna *  tmp_chrom = chrom_1_;
  int32_t         tmp_i     = i_1_;

  chrom_1_  = chrom_2_;
  i_1_      = i_2_;

  chrom_2_  = tmp_chrom;
  i_2_      = tmp_i;
}

inline void VisAVis::copy(VisAVis * source)
{
  i_1_ = source->i_1_;
  i_2_ = source->i_2_;
  chrom_1_ = source->chrom_1_;
  chrom_2_ = source->chrom_2_;
  sense_ = source->sense_;
  score_ = source->score_;
}

inline void VisAVis::check_indices()
{
  i_1_ = Utils::mod(i_1_, chrom_1_->length());
  i_2_ = Utils::mod(i_2_, chrom_2_->length());
}
} // namespace aevol

#endif // AEVOL_VIS_A_VIS_H_
