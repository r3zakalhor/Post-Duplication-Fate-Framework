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

#ifndef AEVOL_CODON_H_
#define AEVOL_CODON_H_

// =================================================================
//                              Includes
// =================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "Dna.h"
#include "macros.h"

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================

class Codon {
 public :
  // =================================================================
  //                             Constructors
  // =================================================================
  Codon();
  explicit Codon(const Codon &model);
  explicit Codon(int8_t value);

  #ifdef BASE_4
  explicit Codon(char b1, char b2, char b3);
  #endif

  Codon(Dna* genome, Strand strand, int32_t index);
  explicit Codon(gzFile backup_file);

  // =================================================================
  //                             Destructors
  // =================================================================
  virtual ~Codon();

  // =================================================================
  //                              Accessors
  // =================================================================
  int8_t value() { return value_; }

  #ifdef BASE_4
  AminoAcid amino_acid() { return amino_acid_; }
  #endif

  // =================================================================
  //                            Public Methods
  // =================================================================
  #ifdef BASE_2
  bool is_start() { return value_ == CODON_START; }
  bool is_stop() { return value_ == CODON_STOP; }
  #elif BASE_4
  bool is_start() { return amino_acid_ == METHIONINE; }
  bool is_stop() { return amino_acid_ == STOP; }
  #endif

  Codon* copy() { return new Codon(value_); } // TODO(dpa) use copy ctor instead!
  void save(gzFile backup_file);

  #ifdef BASE_4
  inline AminoAcid to_aminoacid() const;
  #endif

 protected :
  // =================================================================
  //                           Protected Methods
  // =================================================================

  // =================================================================
  //                          Protected Attributes
  // =================================================================
  int8_t value_;

  #ifdef BASE_4
  AminoAcid amino_acid_ = -1;
  #endif
};

// =====================================================================
//                          Accessors' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

#ifdef BASE_4
inline int8_t bases_to_codon_value(char b1, char b2, char b3) {
  return ((b1-'0') << 4) | ((b2-'0') << 2) | (b3-'0');
}

inline AminoAcid codon_value_to_aminoacid(int8_t value) {
  char base_1 = '0' + ((value >> 4) & 3),
      base_2 = '0' + ((value >> 2) & 3),
      base_3 = '0' + (value & 3);

  if(base_2 == BASE_T)                                  // ** -T- **
  {
    if(base_1 == BASE_G) return VALINE;                   // GT-
    if(base_1 == BASE_C) return LEUCINE;                  // CT-

    if(base_1 == BASE_A) {
      if(base_3 == BASE_G) return METHIONINE;             // ATG
      else return ISOLEUCINE;                             // AT(T/C/A)
    }

    // base_1 == BASE_T
    if(base_3 == BASE_T || base_3 == BASE_C)
      return PHENYLALANINE;                               // TT(T/C)
    else
      return LEUCINE;                                     // TT(A/G)
  }
  else if(base_2 == BASE_C)                             // ** -C- **
  {
    switch(base_1)
    {
    case BASE_T:
      return SERINE;                                    // TC-

    case BASE_C:
      return PROLINE;                                   // CC-

    case BASE_A:
      return THREONINE;                                 // AC-

    case BASE_G:
      return ALANINE;                                   // GC-
    }
  }
  else if(base_2 == BASE_A)                             // ** -A- **
  {
    switch(base_1)
    {
    case BASE_T:
      if(base_3 == BASE_T || base_3 == BASE_C)
        return TYROSINE;                                // TA(T/C)
      else
        return STOP;                                    // TA(A/G)

    case BASE_C:
      if(base_3 == BASE_T || base_3 == BASE_C)
        return HISTIDINE;                               // CA(T/C)
      else
        return GLUTAMINE;                               // CA(A/G)

    case BASE_A:
      if(base_3 == BASE_T || base_3 == BASE_C)
        return ASPARAGINE;                              // AA(T/C)
      else
        return LYSINE;                                  // AA(A/G)

    case BASE_G:
      if(base_3 == BASE_T || base_3 == BASE_C)
        return ASPARTIC_ACID;                           // GA(T/C)
      else
        return GLUTAMIC_ACID;                           // GA(A/G)
    }
  }
  else if(base_2 == BASE_G)                             // ** -G- **
  {
    switch(base_1)
    {
    case BASE_T:
      if(base_3 == BASE_A) return STOP;                 // TGA
      if(base_3 == BASE_G) return TRYPTOPHAN;           // TGG

      return CYSTEINE;                                  // TG(T/C)

    case BASE_C:
      return ARGININE;                                  // CG-

    case BASE_A:
      if(base_3 == BASE_A || base_3 == BASE_G)
        return ARGININE;                                // CG(A/G)

      return SERINE;                                    // CG(T/C)

    case BASE_G:
      return GLYCINE;                                   // GG-
    }
  }

  std::cerr << "An error occurred: codon "
            << base_1 << "," << base_2 << "," << base_3
            << " cannot be translated to aminoacid. Aborting." << std::endl;
  exit(1);
}

inline bool is_start_codon(int8_t value) {
  return codon_value_to_aminoacid(value) == METHIONINE;
}

inline AminoAcid Codon::to_aminoacid() const {
  return codon_value_to_aminoacid(value_);
}
#endif

} // namespace aevol
#endif // AEVOL_CODON_H_
