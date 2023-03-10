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


#ifndef AEVOL_RNA_7_H
#define AEVOL_RNA_7_H


#include <list>
#include <cstdint>

#include "ae_enums.h"
#include "Individual_7.h"
#include "ExpManager.h"
#include "Promoter.h"
namespace aevol {
class Protein_7;
class promoterStruct;

class AffinityFactor {
 public:
  AffinityFactor(Protein_7* prot, double efactor, double ofactor) {
    protein = prot;
    enhancer_factor = efactor;
    operator_factor = ofactor;
  }

  double concentration();

  Protein_7* protein;
  double enhancer_factor;
  double operator_factor;
};


class Rna_7 {
 public:
 /*
  RNA Constructor
  */
  Rna_7(){};
  Rna_7(int32_t t_begin,
        int32_t t_end,
        bool t_leading_lagging,
        double t_e,
        int32_t t_length,PromoterStruct* promo = nullptr);

  Rna_7(const Rna_7& clone) {
    pos                = clone.pos;
    error              = clone.error;
    leading_or_lagging = clone.leading_or_lagging;
    begin             = clone.begin;
    end               = clone.end;
    leading_lagging   = clone.leading_lagging;
    e                 = clone.e;
    length            = clone.length;
    is_coding_        = clone.is_coding_;
    is_init_          = clone.is_init_;
    start_prot_count_ = clone.start_prot_count_;
    
    // if (clone.start_prot != nullptr) {
    //   start_prot = new std::list<int32_t>();
      for (auto pos : (clone.start_prot))
        start_prot.push_back(pos);
    // }
    // start_prot        = clone.start_prot;
    //to_compute_end_ =   clone.to_compute_end_;
    to_compute_start_pos_ = clone.to_compute_start_pos_;
    prom = clone.prom;
  }

  Rna_7(Rna_7* clone,PromoterStruct* promo = nullptr);

  ~Rna_7() { 
  #ifdef __REGUL
  affinity_list.clear();
  #endif

  }

void init(PromoterStruct* prom_x, int32_t t_begin,
        int32_t t_end,
        int8_t t_leading_lagging,
        double t_e,
        int32_t t_length) {
    begin             = t_begin;
    end               = t_end;
    leading_lagging   = t_leading_lagging;
    e                 = t_e;
    length            = t_length;
    is_coding_        = false;
    is_init_          = true;
    start_prot_count_ = 0;
    start_prot.clear();
    // if (start_prot != nullptr)
    //   delete start_prot;
    // start_prot = new std::list<int32_t>();
    //to_compute_end_ = false;
    to_compute_start_pos_ = true;
    prom = prom_x;
    prom_x->rna = this;
  }

#ifdef __REGUL
  std::vector<AffinityFactor> affinity_list;

    int nb_influences_ = 0;

    int32_t enhancer_position(int32_t length) {
      if(leading_lagging == LEADING)
      {
        return Utils::mod(begin - 20, length);
      }
      else  // strand_ = LAGGING
      {
        return Utils::mod(begin + 20, length);
      }
    }

  int32_t operator_position(int32_t length) {
    if(leading_lagging == LEADING)
    {
      return Utils::mod(begin + PROM_SIZE, length);
    }
    else  // strand_ = LAGGING
    {
      return Utils::mod(begin - PROM_SIZE, length);
    }
  }
  double compute_synthesis_rate(Individual_7* indiv);

  double affinity_with_protein( int32_t index, Protein_7* protein,
                                     Individual_7* indiv,
                               ExpManager* exp_m);
#endif

  /*
    Promoter back compatibility
  */ 
  int32_t pos  = -1;
  int8_t error = -1;
  bool leading_or_lagging;

  /*
    RNA variables
  */
  int32_t begin = -1;
  int32_t end = -1;
  int8_t leading_lagging = -1; // 0 = leading, 1 = lagging
  double e = -1.0;
  std::list<int32_t> start_prot;
  std::list<Protein_7*> protein_list_;
  int32_t start_prot_count_ = 0;
  int32_t length = -1;
  bool is_coding_ = false;

  bool is_init_ = false;
  bool to_compute_end_ = true;
  bool to_compute_start_pos_ = true;


  PromoterStruct* prom = nullptr;

  int32_t network_id = -2;
};

}

#endif //AEVOL_RNA_7_H
