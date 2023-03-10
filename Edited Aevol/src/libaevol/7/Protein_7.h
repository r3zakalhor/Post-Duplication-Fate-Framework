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

#ifndef AEVOL_PROTEIN_7_H
#define AEVOL_PROTEIN_7_H


#include <cstdint>
#include "Rna_7.h"

namespace aevol {

class Protein_7 {
 public:
  Protein_7(){};
  Protein_7(int32_t t_protein_start,
            int32_t t_protein_end,
            int32_t t_protein_length,
            int8_t t_leading_lagging,
            double t_e, Rna_7* rna, bool translated = false) {

              
    protein_start   = t_protein_start;
    protein_end     = t_protein_end;
    protein_length  = t_protein_length;
    leading_lagging = t_leading_lagging;
    e               = t_e;
    is_init_        = true;

#ifdef __REGUL
    rna_list_.push_back(rna);
    rna->protein_list_.push_back(this);
    initial_e_ = e;
    translated_ = translated;
#endif

#ifdef BASE_4
    rna_list_.push_back(rna);
#endif
    
    // printf("P_CLASSIC ??? %f %f %f %f :: %d %d %d\n",h,w,m,e,is_init_,signal_,inherited_);
// if (t_e > 1) {
//                                     printf("P_TO_CREATE %f\n",t_e);

//                                     printf("ERROR CREATE e is too HIGH\n");
//                                     exit(-22);
//                                 }
  }

#ifdef __REGUL
  Protein_7(Protein_R* prot_sig);


  Protein_7(Protein_7* prot, ExpManager* exp_m);
#endif

  bool operator<(const Protein_7& other);

  ~Protein_7() {
    rna_list_.clear();
  }

  int32_t protein_start = -1;
  int32_t protein_end = -1;
  int32_t protein_length = -1;
  int8_t leading_lagging = -1; // 0 = leading, 1 = lagging
  double m = -1.0;
  double w = -1.0;
  double h = -1.0;
  double e = -1.0;
  bool is_functional = false;

  bool is_init_ = false;

  #ifdef BASE_2
  int8_t codon_list[64*4] = {};
  #elif BASE_4
  vector<int8_t> codon_list = {};
  #endif

  int16_t nb_codons_ = 0;

#ifdef __REGUL
  bool is_TF_ = false;

  double    delta_concentration_;
  bool      inherited_ = false;
  bool      signal_ = false;
  bool      translated_ = false;
#endif

  std::list<Rna_7*> rna_list_;
  double initial_e_ = -1;

  int32_t protein_id_ = -1; // Only set and used by post treatments.

  int32_t network_id = -2;
};
}

#endif //AEVOL_PROTEIN_7_H
