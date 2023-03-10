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
//                              Includes
// =================================================================
#include "Codon.h"

#include <cmath>

#include "Utils.h"

namespace aevol {



// ############################################################################
//
//                                Class Codon
//
// ############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
Codon::Codon() {
  value_ = -1;
}

Codon::Codon(const Codon &model) {
  value_ = model.value_;
  #ifdef BASE_4
  amino_acid_ = model.amino_acid_;
  #endif
}

Codon::Codon(int8_t value) {
  value_ = value;
  #ifdef BASE_4
  amino_acid_ =  to_aminoacid();
  #endif
}

#ifdef BASE_4
Codon::Codon(char b1, char b2, char b3) {
  amino_acid_ = to_aminoacid();
  value_ = bases_to_codon_value(b1, b2, b3);
}
#endif

Codon::Codon(Dna * dna, Strand strand, int32_t index) {
  const char* gen = dna->data();
  int32_t     len = dna->length();

  #ifdef BASE_2
  value_ = 0;

  if (strand == LEADING) {
    for (int8_t i = 0 ; i < CODON_SIZE ; i++) {
      if (gen[Utils::mod((index+i), len)] == '1') {
        value_ += 1 << (CODON_SIZE - i - 1);
      }
    }
  }
  else { // (strand == LAGGING)
    for (int8_t i = 0 ; i < CODON_SIZE ; i++) {
      if (gen[Utils::mod((index-i), len)] != '1') {
        value_ += 1 << (CODON_SIZE - i - 1);
      }
    }
  }
  #elif BASE_4
    char base1,
       base2,
       base3;

  if(strand == LEADING) {
    base1 = gen[Utils::mod(index, len)];
    base2 = gen[Utils::mod(index + 1, len)];
    base3 = gen[Utils::mod(index + 2, len)];
  }
  else { // (strand == LAGGING)
    // lagging strain, need to get the complementary base
    base1 = get_complementary_base(gen[Utils::mod(index, len)]);
    base2 = get_complementary_base(gen[Utils::mod(index - 1, len)]);
    base3 = get_complementary_base(gen[Utils::mod(index - 2, len)]);
  }

  value_ = bases_to_codon_value(base1, base2, base3);
  amino_acid_ = to_aminoacid();
  #endif
}

Codon::Codon(gzFile backup_file) {
  gzread(backup_file, &value_, sizeof(value_));
  #ifdef BASE_4
  gzread(backup_file, &amino_acid_, sizeof(amino_acid_));
  #endif
}

// =================================================================
//                             Destructors
// =================================================================
Codon::~Codon() = default;

// =================================================================
//                            Public Methods
// =================================================================
void Codon::save(gzFile backup_file) {
  gzwrite(backup_file, &value_, sizeof(value_));
  #ifdef BASE_4
  gzwrite(backup_file, &amino_acid_, sizeof(amino_acid_));
  #endif
}
// =================================================================
//                           Protected Methods
// =================================================================
} // namespace aevol
