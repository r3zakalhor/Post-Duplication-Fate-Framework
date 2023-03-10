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


#include "Rna_7.h"
#include "ae_enums.h"
#include "Protein_7.h"
#include "Individual_7.h"
#include "ExpManager_7.h"
namespace aevol {

    Rna_7::Rna_7(int32_t t_begin,
        int32_t t_end,
        bool t_leading_lagging,
        double t_e,
        int32_t t_length,PromoterStruct* promo) {
    begin             = t_begin;
    end               = t_end;
    leading_lagging   = t_leading_lagging;
    e                 = t_e;
    length            = t_length;
    is_coding_        = false;
    is_init_          = true;
    start_prot_count_ = 0;
    to_compute_start_pos_ = true;
    prom = promo;
    promo->rna = this;
    //  start_prot = new std::list<int32_t>();
  }

Rna_7::Rna_7(Rna_7* clone,PromoterStruct* promo) {
    pos                = clone->pos;
    error              = clone->error;
    leading_or_lagging = clone->leading_or_lagging;
    begin             = clone->begin;
    end               = clone->end;
    leading_lagging   = clone->leading_lagging;
    e                 = clone->e;
    length            = clone->length;
    is_coding_        = clone->is_coding_;
    is_init_          = clone->is_init_;
    start_prot_count_ = clone->start_prot_count_;
    // start_prot = new std::list<int32_t>();
    // if (clone->start_prot != nullptr) {
      for (auto pos : (clone->start_prot))
        start_prot.push_back(pos);
    // }
    //to_compute_end_ =   clone->to_compute_end_;
    if (promo != nullptr)
      prom = promo;
    else
      prom = clone->prom;

    to_compute_start_pos_ = clone->to_compute_start_pos_;
    promo->rna = this;
  }
#ifdef __REGUL
double Rna_7::affinity_with_protein( int32_t index, Protein_7 *protein, Individual_7* indiv, ExpManager* exp_m ) {
  int32_t len = protein->protein_length;

  if (len > 5) {
    double max = 0;
    double temp = 1;

    int32_t quadon_tab[5];

    for (int32_t pos = index; pos < index+5; pos++) {

      int8_t quadon[4];

      if (leading_lagging == LEADING) {
        for (int8_t i = 0; i < QUADON_SIZE; i++) {
          quadon[i] = (indiv->dna_->get_lead(pos + i) == '1')
                      ? 1 << (QUADON_SIZE - i - 1)
                      : 0;
        }
      } else {
        for (int8_t i = 0; i < QUADON_SIZE; i++) {
          quadon[i] = (indiv->dna_->get_lag(pos - i) != '1')
                      ? 1 << (QUADON_SIZE - i - 1)
                      : 0;
        }
      }

      // printf("Quadon %d : %d / %d %d %d %d\n",index,pos,quadon[0],quadon[1],quadon[2],quadon[3]);
      quadon_tab[pos - index] = quadon[0] + quadon[1] + quadon[2] + quadon[3];
    }

    for (int32_t i = 0; i < len - 4; i++) {
      temp = 1;

      for (int8_t j = 0; j < 5; j++) {
        if (protein->codon_list[i+j]>=MAX_CODON || quadon_tab[j]>=MAX_QUADON) {
          printf("Individual %d Protein %d Init %d\n",indiv->indiv_id,protein->protein_start,protein->is_init_);
          printf("Codon[%d] (i %d j %d) %d out of %d\n",i+j,i,j,protein->codon_list[i+j],MAX_CODON);
          printf("Protein Length %d\n",protein->protein_length);
        } else {
          temp *= exp_m->exp_s()->get_binding_matrix(quadon_tab[j],
                                                   protein->codon_list[i + j]);
        }
      }
      // if (temp > 0.0)
        // printf("%d -- %d => %lf\n",indiv->indiv_id,index,
        //                                             temp);
      max = (max < temp) ? temp : max;
    }

    return max;
  } else {
    return 0.0;
  }
}

double Rna_7::compute_synthesis_rate(Individual_7* indiv) {
  
  if (is_init_) {

    double enhancer_activity = 0;
    double operator_activity = 0;

//for (auto affinity = affinity_list.begin(); affinity != affinity_list.end(); affinity++) {
    #pragma omp simd reduction(+:enhancer_activity,operator_activity) 
    for (int32_t i = 0; i < affinity_list.size(); i++) {
      enhancer_activity +=
          affinity_list[i].enhancer_factor * affinity_list[i].concentration();
      operator_activity +=
          affinity_list[i].operator_factor * affinity_list[i].concentration();

    //  if (AeTime::time() == 85303)
      //  printf("%d -- %d -- RNA %d Protein %d (%lf) :: Enhancer %lf Operator %lf\n",AeTime::time(),
      //                            #ifdef HAVE_MPI
      // indiv->exp_m_->exp_m_7_->localXtoGlobalX(indiv->indiv_id / indiv->exp_m_->exp_m_7_->grid_height_)*indiv->exp_m_->exp_m_7_->global_grid_height_+indiv->exp_m_->exp_m_7_->localYtoGlobalY(indiv->indiv_id % indiv->exp_m_->exp_m_7_->grid_height_)
      // indiv_id / grid_height_,indiv_id % grid_height_,
      // localXtoGlobalX(indiv_id / grid_height_),localYtoGlobalY(indiv_id % grid_height_),
      // #else
      // indiv->indiv_id
      // #endif
      // ,begin,affinity.protein->protein_start,
      //         affinity.concentration(), affinity.enhancer_factor, affinity.operator_factor);
    }

//    if (indiv->indiv_id==137)
    //  if (enhancer_activity !=0 ||  operator_activity != 0) 
    //   printf("%d -- %d -- RNA %d Enhancer %lf Operator %lf\n",AeTime::time(),
    //                              #ifdef HAVE_MPI
    //   indiv->exp_m_->exp_m_7_->localXtoGlobalX(indiv->indiv_id / indiv->exp_m_->exp_m_7_->grid_height_)*indiv->exp_m_->exp_m_7_->global_grid_height_+indiv->exp_m_->exp_m_7_->localYtoGlobalY(indiv->indiv_id % indiv->exp_m_->exp_m_7_->grid_height_)
    //   // indiv_id / grid_height_,indiv_id % grid_height_,
    //   // localXtoGlobalX(indiv_id / grid_height_),localYtoGlobalY(indiv_id % grid_height_),
    //   #else
    //   indiv->indiv_id
    //   #endif
    //       ,begin,enhancer_activity,operator_activity);

    ProteinConcentration enhancer_activity_pow_n =
        enhancer_activity == 0
        ? 0
        : pow(enhancer_activity,
              indiv->exp_m_->exp_s()->get_hill_shape_n());
    ProteinConcentration operator_activity_pow_n =
        operator_activity == 0
        ? 0
        : pow(operator_activity,
              indiv->exp_m_->exp_s()->get_hill_shape_n());

    // if ((indiv->exp_m_->exp_s()->get_hill_shape() /
    //      (operator_activity_pow_n + indiv->exp_m_->exp_s()->get_hill_shape())) *
    //     (1 + ((1 / e) - 1) * (enhancer_activity_pow_n /
    //                           (enhancer_activity_pow_n +
    //                            indiv->exp_m_->exp_s()->get_hill_shape()))) != 1.0)
    // printf("%d -- %d -- RNA %d : E %lf New_E %lf -- OP %lf OP_POW %lf :: EP %lf EP_POW %lf\n",
    // AeTime::time(),
    //                              #ifdef HAVE_MPI
    //   indiv->exp_m_->exp_m_7_->localXtoGlobalX(indiv->indiv_id / indiv->exp_m_->exp_m_7_->grid_height_)*indiv->exp_m_->exp_m_7_->global_grid_height_+indiv->exp_m_->exp_m_7_->localYtoGlobalY(indiv->indiv_id % indiv->exp_m_->exp_m_7_->grid_height_)
    //   // indiv_id / grid_height_,indiv_id % grid_height_,
    //   // localXtoGlobalX(indiv_id / grid_height_),localYtoGlobalY(indiv_id % grid_height_),
    //   #else
    //   indiv->indiv_id
    //   #endif
    //       ,
    // begin,e,
    //                            e *
    //     (indiv->exp_m_->exp_s()->get_hill_shape() /
    //      (operator_activity_pow_n + indiv->exp_m_->exp_s()->get_hill_shape())) *
    //     (1 + ((1 / e) - 1) * (enhancer_activity_pow_n /
    //                           (enhancer_activity_pow_n +
    //                            indiv->exp_m_->exp_s()->get_hill_shape()))), operator_activity,operator_activity_pow_n,enhancer_activity,enhancer_activity_pow_n);
    return
        e *
        (indiv->exp_m_->exp_s()->get_hill_shape() /
         (operator_activity_pow_n + indiv->exp_m_->exp_s()->get_hill_shape())) *
        (1 + ((1 / e) - 1) * (enhancer_activity_pow_n /
                              (enhancer_activity_pow_n +
                               indiv->exp_m_->exp_s()->get_hill_shape())));
  }
  return 0;
}

double AffinityFactor::concentration() { return protein->e; };
#endif
}
