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


#include "Protein_7.h"
#ifdef __REGUL
#include "raevol/Protein_R.h"
#endif

namespace aevol {
#ifdef __REGUL
Protein_7::Protein_7(Protein_R* prot_sig) {
  
  protein_length = prot_sig->length();
  e              = prot_sig->concentration();

  signal_ = prot_sig->is_signal();

  nb_codons_ = prot_sig->AA_list().size();

  w = 0;
  m = 0;
  h = 0;
// printf("P_SIG ??? %f %f %f %f :: %d %d %d\n",h,w,m,e,is_init_,signal_,inherited_);
  for (int i = 0; i < nb_codons_; i++)
    codon_list[i] = prot_sig->_cod_tab[i];

}

Protein_7::Protein_7(Protein_7* prot, ExpManager* exp_m) {
  protein_start   = prot->protein_start;
  protein_end     = prot->protein_end;
  protein_length  = prot->protein_length;
  
  // printf("Protein Length %d :: %d\n",protein_length,prot->protein_length);

  leading_lagging = prot->leading_lagging;

  if (exp_m->exp_s()->get_with_heredity() && !prot->signal_) {
    e               = prot->e;
  } else {
    e               = prot->e;
  }

  w = prot->w;
  m = prot->m;
  h = prot->h;

  
  signal_         = prot->signal_;
  nb_codons_      = prot->nb_codons_;

  for (int i = 0; i < nb_codons_; i++)
    codon_list[i] = prot->codon_list[i];

  initial_e_ = prot->initial_e_;
  if (! signal_) {
    inherited_ = true;
    is_init_        = true;
  }
  
  translated_ = false;
  if (! exp_m->exp_s()->get_with_heredity()) {
    for (auto& rna : prot->rna_list_)
      rna_list_.push_back(rna);
  }
  // printf("P_INHERIT_CPY %f %f %f %f :: %d %d %d\n",prot->h,prot->w,prot->m,prot->e,prot->is_init_,prot->signal_,prot->inherited_);
  // printf("P_INHERIT_THIS %f %f %f %f :: %d %d %d\n",h,w,m,e,is_init_,signal_,inherited_);
}
#endif
    bool Protein_7::operator<(const Protein_7& other){
      // printf("P_THIS %f %f %f %f :: %d %d %d\n",h,w,m,e,is_init_,signal_,inherited_);
      // printf("P_OTHER %f %f %f %f :: %d %d %d\n",other.h,other.w,other.m,other.e,is_init_,signal_,inherited_);
        return (h <  other.h)
               || (h == other.h && m < other.m)
               || (h == other.h && m == other.m && w < other.w)
               || (h == other.h && m == other.m && w == other.w && e < other.e);
    }

    // bool Protein_7::operator<(const Protein_7& other){



    //     return (h <  other.h)
    //            || (h == other.h && m < other.m)
    //            || (h == other.h && m == other.m && w < other.w)
    //            || (h == other.h && m == other.m && w == other.w && e < other.e);
    // }

}