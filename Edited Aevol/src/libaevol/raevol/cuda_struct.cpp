//
// Created by arrouan on 10/02/16.
//

#include "cuda_struct.h"

namespace aevol {
/*
void cuda_struct::init_struct(int max_protein, int max_rna, int max_influence,
                              int nb_signals, int life_time, int nb_eval,
                              double selection_pressure) {
  max_protein_ = max_protein;
  max_rna_ = max_rna;
  max_influence_ = max_influence;
  nb_signals_ = nb_signals;
  life_time_ = life_time;
  nb_eval_ = nb_eval;
  selection_pressure_ = selection_pressure;

  phenotype_inhib = new float[300 * 1024];
  phenotype_activ = new float[300 * 1024];
  phenotype = new float[300 * 1024];
  dist_sum = new float[1024];

  environment = new float[300 * life_time];
  delta = new float[300 * 1024];

  eval_step = new int[life_time];

  protein_concentration = new float[max_protein * 1024];

  signals_concentration = new float[nb_signals * life_time];

  enhance_coef = new float[max_influence * max_rna * 1024];
  operate_coef = new float[max_influence * max_rna * 1024];

  rna_synthesis = new float[max_rna * 1024];
  basal_level = new float[max_rna * 1024];

  protein_influence = new int[max_influence * max_rna * 1024];

  protein_influenced = new float[max_protein * max_rna * 1024];

  protein_triangle_ix0 = new float[max_protein * 1024];
  protein_triangle_ix1 = new float[max_protein * 1024];
  protein_triangle_ix2 = new float[max_protein * 1024];

  protein_triangle_height = new float[max_protein * 1024];
}

void cuda_struct::delete_struct() {
  delete[] phenotype_inhib;
  delete[] phenotype_activ;
  delete[] phenotype;
  delete[] dist_sum;

  delete[] environment;
  delete[] delta;

  delete[] eval_step;

  delete[] protein_concentration;

  delete[] signals_concentration;

  delete[] enhance_coef;
  delete[] operate_coef;
  delete[] rna_synthesis;
  delete[] basal_level;

  delete[] protein_influence;

  delete[] protein_influenced;

  delete[] protein_triangle_ix0;
  delete[] protein_triangle_ix1;
  delete[] protein_triangle_ix2;

  delete[] protein_triangle_height;
}


void cuda_struct::transfert_to_gpu(ExpManager* exp_m) {
  std::set<int>* eval = exp_m->exp_s()->get_list_eval_step();
  const Habitat_R& hab = dynamic_cast<const Habitat_R&>(exp_m->world()->grid(0,0)->habitat());

  for (int8_t i = 1; i <= exp_m->exp_s()->get_nb_indiv_age(); i++) {
    //Set the concentration of signals for this age
    for (Protein_R* prot1 : hab.signals()) {
      signals_concentration[(i-1) * nb_signals_ + prot1->get_local_id()] = 0.0;
      //printf("Index %d (%d) MAX %d\n",(i-1) * life_time_ + prot1->get_local_id(),
      //       prot1->get_local_id(),life_time_*nb_signals_);

 //     printf("A -- Signal %d at %d : %f\n",prot1->get_local_id(),i,signals_concentration[i*life_time_+prot1->get_local_id()]);
    }

    for (Protein_R* prot2 : hab.phenotypic_target(i).signals()) {
      signals_concentration[(i-1) * nb_signals_ + prot2->get_local_id()] = 0.9;
//      printf("B -- Signal %d at %d : %f\n",prot2->get_local_id(),i,signals_concentration[i*life_time_+prot2->get_local_id()]);
    }

//    for (int j = 0; j < nb_signals_; j++) {
//      printf("C -- Signal %d at %d : %f\n",j,i,signals_concentration[i*life_time_+j]);
//    }

    if (eval->find(i) != eval->end())
      eval_step[(i-1)] = 1;
    else
      eval_step[(i-1)] = 0;

    for (int j = 0; j < 300; j++)
      environment[(i-1) * 300 + j] = ((HybridFuzzy*) hab.phenotypic_target(
          i).fuzzy())->points()[j];
  }


  int i = 0;
  for (int ki = 0; ki < 32; ki++)
    for (int kj = 0; kj < 32; kj++) {
      Individual_R* indiv = dynamic_cast<Individual_R*>(exp_m->world()->indiv_at(ki, kj));

      i = ki * 32 + kj;

      for (int j = 0; j < 300; j++) {
        phenotype[i * 300 + j] = 0;
        phenotype_activ[i * 300 + j] = 0;
        phenotype_inhib[i * 300 + j] = 0;
        delta[i * 300 + j] = 0;
      }

      for (int je = 0; je < max_protein_; je++) {
        protein_concentration[i * max_protein_ + je] = 0;

        protein_triangle_ix0[i * max_protein_ + je] = 0;
        protein_triangle_ix1[i * max_protein_ + je] = 0;
        protein_triangle_ix2[i * max_protein_ + je] = 0;

        protein_triangle_height[i * max_protein_ + je] = 0;


        for (int k = 0; k < max_rna_; k++)
          protein_influenced[i * max_protein_ * max_rna_ + je * max_protein_ + k] = 0;
      }

      for (auto prot_a : indiv->protein_list()) {
        Protein_R* prot = dynamic_cast<Protein_R*>(prot_a);

        prot->reset_concentration();
        protein_concentration[i * max_protein_ +
                              prot->get_local_id()] = prot->concentration();

        float x0 = prot->mean() - prot->width();
        float x1 = prot->mean();
        float x2 = prot->mean() + prot->width();

        int ix0 = (int) (x0 * 300);
        int ix1 = (int) (x1 * 300);
        int ix2 = (int) (x2 * 300);

        if (ix0 < 0) ix0 = 0; else if (ix0 > (300 - 1)) ix0 = 300 - 1;
        if (ix1 < 0) ix1 = 0; else if (ix1 > (300 - 1)) ix1 = 300 - 1;
        if (ix2 < 0) ix2 = 0; else if (ix2 > (300 - 1)) ix2 = 300 - 1;

        protein_triangle_ix0[i * max_protein_ + prot->get_local_id()] = ix0;
        protein_triangle_ix1[i * max_protein_ + prot->get_local_id()] = ix1;
        protein_triangle_ix2[i * max_protein_ + prot->get_local_id()] = ix2;

        protein_triangle_height[i * max_protein_ +
                                prot->get_local_id()] = prot->height();



        for (auto rna : prot->_rna_R_list)
          protein_influenced[i * max_protein_ * max_rna_ +
                             prot->get_local_id() * max_protein_
                             + rna->get_local_id()] = 1;
      }

      for (int je = 0; je < max_rna_; je++) {
        rna_synthesis[i * max_rna_ +
                      je] = 0;
        basal_level[i * max_rna_ +
                      je] = 0;
        for (int ke = 0; ke < max_influence_; ke++) {
          enhance_coef[i * max_rna_ * max_influence_ +
                       je * max_influence_ + ke] = 0;
          operate_coef[i * max_rna_ * max_influence_ +
                       je * max_influence_ + ke] = 0;
          protein_influence[i * max_rna_ * max_influence_ +
                            je * max_influence_ + ke] = 0;
        }
      }


      for (auto rna : indiv->_rna_list_coding) {

        rna_synthesis[i * max_rna_ +
                      rna->get_local_id()] = rna->get_synthesis_rate();

        basal_level[i * max_rna_ +
                      rna->get_local_id()] = rna->basal_level();

        for (int k = 0; k < rna->_nb_influences; k++) {
          enhance_coef[i * max_rna_ * max_influence_ +
                       rna->get_local_id() * max_influence_ +
                       k] = rna->_enhancing_coef_list[k];
          operate_coef[i * max_rna_ * max_influence_ +
                       rna->get_local_id() * max_influence_ +
                       k] = rna->_operating_coef_list[k];
          if (rna->_protein_list[k]->is_signal())
            protein_influence[i * max_rna_ * max_influence_ +
                              rna->get_local_id() * max_influence_ +
                              k] =
                rna->_protein_list[k]->get_local_id() + max_protein_;
          else
            protein_influence[i * max_rna_ * max_influence_ +
                              rna->get_local_id() * max_influence_ +
                              k] = rna->_protein_list[k]->get_local_id();
        }
      }

    }


}
*/
void cuda_struct::compute_a_generation(ExpManager* exp_m) {

  int degradation_step = exp_m->exp_s()->get_nb_degradation_step();
  float degradation_rate = exp_m->exp_s()->get_degradation_rate();
  float hill_shape_n = exp_m->exp_s()->get_hill_shape_n();
  float hill_shape = exp_m->exp_s()->get_hill_shape();

#pragma omp parallel for collapse(2)
  for (int ki = 0; ki < 32; ki++)
    for (int kj = 0; kj < 32; kj++) {
      int indiv_id = ki * 32 + kj;

      for (int8_t i = 0; i < life_time_; i++) {
        for (int j = 0; j < degradation_step; j++) {
          // Compute synthesis rate of RNA
          for (int m = 0; m < max_rna_; m++) {
            float enhancer_activity = 0, operator_activity = 0;
#ifdef __SIMD
            #pragma omp simd
#endif
            for (int n = 0; n < max_influence_; n++) {
              int prot_id = protein_influence
              [ki * 32 + kj * max_rna_ * max_influence_ + m * max_influence_ +
               n];

              if (prot_id > max_protein_) {
                enhancer_activity +=
                    enhance_coef[indiv_id * max_rna_ * max_influence_ +
                                 m * max_influence_ + n] *
                    signals_concentration[i * nb_signals_ + prot_id -
                                          max_protein_];

                operator_activity +=
                    operate_coef[indiv_id * max_rna_ * max_influence_ +
                                 m * max_influence_ + n] *
                    signals_concentration[i * nb_signals_ + prot_id -
                                          max_protein_];
              } else {
                enhancer_activity +=
                    enhance_coef[indiv_id * max_rna_ * max_influence_ +
                                 m * max_influence_ + n] *
                    protein_concentration[indiv_id * max_protein_ + prot_id];

                operator_activity +=
                    operate_coef[indiv_id * max_rna_ * max_influence_ +
                                 m * max_influence_ + n] *
                    protein_concentration[indiv_id * max_protein_ + prot_id];
              }
            }

            float enhancer_activity_pow_n = pow(enhancer_activity,
                                                hill_shape_n);
            float operator_activity_pow_n = pow(operator_activity,
                                                hill_shape_n);

            rna_synthesis[indiv_id * max_rna_ + m] = basal_level[indiv_id * max_rna_ + m]
                                                         * (hill_shape
                                                            /
                                                            (operator_activity_pow_n +
                                                             hill_shape))
                                                         * (1 + ((1 /
                basal_level[indiv_id * max_rna_ + m]) -
                                                                 1)
                                                                *
                                                                (enhancer_activity_pow_n /
                                                                 (enhancer_activity_pow_n +
                                                                  hill_shape)));
          }

          // Compute concentration
          for (int l = 0; l < max_protein_; l++) {
            float _delta_concentration;

            for (int m = 0; m < max_rna_; m++) {
              if (protein_influenced[indiv_id * max_protein_ * max_rna_ +
                                     l * max_protein_
                                     + m])
                _delta_concentration += rna_synthesis[indiv_id * max_rna_ +
                                                      m];
            }

            _delta_concentration -= degradation_rate * protein_concentration[indiv_id * max_protein_ +
                                                                             l];
            _delta_concentration *= 1 / ((float) degradation_step);

            protein_concentration[indiv_id * max_protein_ +
                                  l] += _delta_concentration;
          }
        }

        // If we have to evaluate the individual at this age
        if (eval_step[i]) {
          /** update phenotype **/
#ifdef __SIMD
          #pragma omp simd
#endif
          for (int j = 0; j < 300; j++) {
            phenotype[indiv_id * 300 + j] = 0;
            phenotype_activ[indiv_id * 300 + j] = 0;
            phenotype_inhib[indiv_id * 300 + j] = 0;
          }

          for (int l = 0; l < max_protein_; l++) {
            float height = (protein_triangle_height[indiv_id * max_protein_ +
                                                    l] *
                            protein_concentration[indiv_id * max_protein_ + l]);
            float incY =
                protein_triangle_height[indiv_id * max_protein_ + l] *
                protein_concentration[indiv_id * max_protein_ + l] /
                (protein_triangle_ix1[indiv_id * max_protein_ + l] -
                 protein_triangle_ix0[indiv_id * max_protein_ + l]);

            if (protein_triangle_height[indiv_id * max_protein_ + l] > 0) {
              for (int j = 0; j < 300; j++) {
                float fj = (float) j;
                if (fj > protein_triangle_ix0[indiv_id * max_protein_ + l]
                    and
                    fj < protein_triangle_ix1[indiv_id * max_protein_ + l])
                  phenotype_activ[indiv_id * 300 + j] +=
                      incY * (fj -
                              protein_triangle_ix0[indiv_id * max_protein_ +
                                                   l]);
                else if (
                    fj > protein_triangle_ix1[indiv_id * max_protein_ + l]
                    and
                    fj < protein_triangle_ix2[indiv_id * max_protein_ + l])
                  phenotype_activ[indiv_id * 300 + j] += height
                                                         -
                                                         incY * (fj -
                                                                 protein_triangle_ix1[
                                                                     indiv_id *
                                                                     max_protein_ +
                                                                     l]);
                else if (fj ==
                         protein_triangle_ix1[indiv_id * max_protein_ + l])
                  phenotype_activ[indiv_id * 300 + j] +=
                      height;
              }
            } else {
              for (int j = 0; j < 300; j++) {
                float fj = (float) j;

                if (fj > protein_triangle_ix0[indiv_id * max_protein_ + l]
                    and
                   fj < protein_triangle_ix1[indiv_id * max_protein_ + l])
                  phenotype_inhib[indiv_id * 300 + j] += incY * (fj -
                                                                 protein_triangle_ix0[
                                                                     indiv_id *
                                                                     max_protein_ +
                                                                     l]);
                else if (
                    fj > protein_triangle_ix1[indiv_id * max_protein_ + l]
                    and
                    fj < protein_triangle_ix2[indiv_id * max_protein_ + l])
                  phenotype_inhib[indiv_id * 300 + j] +=
                      height -
                      incY * (fj -
                              protein_triangle_ix1[indiv_id * max_protein_ +
                                                   l]);
                else if (fj ==
                         protein_triangle_ix1[indiv_id * max_protein_ + l])
                  phenotype_inhib[indiv_id * 300 + j] +=
                      height;
              }
            }
          }

#ifdef __SIMD
          #pragma omp simd
#endif
          for (int j = 0; j < 300; j++) {
            phenotype_activ[indiv_id * 300 + j] =
                phenotype_activ[indiv_id * 300 + j] > Y_MAX ? Y_MAX :
                phenotype_activ[indiv_id * 300 + j];
            phenotype_inhib[indiv_id * 300 + j] =
                phenotype_inhib[indiv_id * 300 + j] < -Y_MAX ? -Y_MAX :
                phenotype_inhib[indiv_id * 300 + j];
            phenotype[indiv_id * 300 + j] =
                phenotype_activ[indiv_id * 300 + j] +
                phenotype_inhib[indiv_id * 300 + j];
            phenotype[indiv_id * 300 + j] =
                phenotype[indiv_id * 300 + j] > Y_MIN ?
                Y_MIN : phenotype[indiv_id * 300 + j];
            delta[indiv_id * 300 + j] =
                phenotype[indiv_id * 300 + j] -
                environment[i * 300 + j];
          }

          /** compute distance to target **/

          float area = 0;

#ifdef __SIMD
          #pragma omp simd
#endif
          for (int j = 0; j < 299; j++) {
            dist_sum[indiv_id] += ((delta[indiv_id * 300 + j]
                                    + delta[indiv_id * 300 + j + 1]) / 600.0);
          }
        }
      }

      dist_sum[indiv_id] = exp(
          selection_pressure_ * (dist_sum[indiv_id] / nb_eval_));
    }
}


void cuda_struct::compute_a_generation_v2(ExpManager* exp_m) {

  int degradation_step = exp_m->exp_s()->get_nb_degradation_step();
  float degradation_rate = exp_m->exp_s()->get_degradation_rate();
  float hill_shape_n = exp_m->exp_s()->get_hill_shape_n();
  float hill_shape = exp_m->exp_s()->get_hill_shape();


  for (int8_t i = 0; i < life_time_; i++) {
    for (int j = 0; j < degradation_step; j++) {
      // Compute synthesis rate of RNA
      for (int m = 0; m < max_rna_; m++) {
        float enhancer_activity = 0, operator_activity = 0;
        for (int n = 0; n < max_influence_; n++) {
          #pragma omp parallel for
          for (int indiv_id = 0; indiv_id < 1024; indiv_id++) {

            int prot_id = protein_influence
            [indiv_id * max_rna_ * max_influence_ + m * max_influence_ +
             n];

            if (prot_id > max_protein_) {
              enhancer_activity +=
                  enhance_coef[indiv_id * max_rna_ * max_influence_ +
                               m * max_influence_ + n] *
                  signals_concentration[i * nb_signals_ + prot_id -
                                        max_protein_];

              operator_activity +=
                  operate_coef[indiv_id * max_rna_ * max_influence_ +
                               m * max_influence_ + n] *
                  signals_concentration[i * nb_signals_ + prot_id -
                                        max_protein_];
            } else {
              enhancer_activity +=
                  enhance_coef[indiv_id * max_rna_ * max_influence_ +
                               m * max_influence_ + n] *
                  protein_concentration[indiv_id * max_protein_ + prot_id];

              operator_activity +=
                  operate_coef[indiv_id * max_rna_ * max_influence_ +
                               m * max_influence_ + n] *
                  protein_concentration[indiv_id * max_protein_ + prot_id];
            }
          }
        }

        #pragma omp parallel for
        for (int indiv_id = 0; indiv_id < 1024; indiv_id++) {

          float enhancer_activity_pow_n = pow(enhancer_activity,
                                              hill_shape_n);
          float operator_activity_pow_n = pow(operator_activity,
                                              hill_shape_n);

          rna_synthesis[indiv_id * max_rna_ + m] =
              basal_level[indiv_id * max_rna_ + m]
              * (hill_shape
                 /
                 (operator_activity_pow_n +
                  hill_shape))
              * (1 + ((1 /
                       basal_level[indiv_id * max_rna_ + m]) -
                      1)
                     *
                     (enhancer_activity_pow_n /
                      (enhancer_activity_pow_n +
                       hill_shape)));
        }
      }

      // Compute concentration
      for (int l = 0; l < max_protein_; l++) {
        float _delta_concentration;
        for (int m = 0; m < max_rna_; m++) {
          #pragma omp parallel for
          for (int indiv_id = 0; indiv_id < 1024; indiv_id++)

            if (protein_influenced[indiv_id * max_protein_ * max_rna_ +
                                 l * max_protein_
                                 + m])
              _delta_concentration += rna_synthesis[indiv_id * max_rna_ +
                                                  m];
        }

        #pragma omp parallel for
        for (int indiv_id = 0; indiv_id < 1024; indiv_id++) {

          _delta_concentration -=
              degradation_rate * protein_concentration[indiv_id * max_protein_ +
                                                       l];
          _delta_concentration *= 1 / ((float) degradation_step);

          protein_concentration[indiv_id * max_protein_ +
                                l] += _delta_concentration;
        }
      }
    }

    // If we have to evaluate the individual at this age
    if (eval_step[i]) {
      /** update phenotype **/
      for (int j = 0; j < 300; j++) {
        #pragma omp parallel for
        for (int indiv_id = 0; indiv_id < 1024; indiv_id++) {
          phenotype[indiv_id * 300 + j] = 0;
          phenotype_activ[indiv_id * 300 + j] = 0;
          phenotype_inhib[indiv_id * 300 + j] = 0;
        }
      }

      for (int l = 0; l < max_protein_; l++) {
#pragma omp parallel for
        for (int indiv_id = 0; indiv_id < 1024; indiv_id++) {

          float height = (protein_triangle_height[indiv_id * max_protein_ +
                                                  l] *
                          protein_concentration[indiv_id * max_protein_ + l]);
          float incY =
              protein_triangle_height[indiv_id * max_protein_ + l] *
              protein_concentration[indiv_id * max_protein_ + l] /
              (protein_triangle_ix1[indiv_id * max_protein_ + l] -
               protein_triangle_ix0[indiv_id * max_protein_ + l]);

          if (protein_triangle_height[indiv_id * max_protein_ + l] > 0) {
#ifdef __SIMD
            #pragma omp simd
#endif
            for (int j = 0; j < 300; j++) {
              if (j > protein_triangle_ix0[indiv_id * max_protein_ + l]
                  and
                  j < protein_triangle_ix1[indiv_id * max_protein_ + l])
                phenotype_activ[indiv_id * 300 + j] +=
                    incY * (j -
                            protein_triangle_ix0[indiv_id * max_protein_ +
                                                 l]);
              else if (
                  j > protein_triangle_ix1[indiv_id * max_protein_ + l]
                  and
                  j < protein_triangle_ix2[indiv_id * max_protein_ + l])
                phenotype_activ[indiv_id * 300 + j] += height
                                                       -
                                                       incY * (j -
                                                               protein_triangle_ix1[
                                                                   indiv_id *
                                                                   max_protein_ +
                                                                   l]);
              else if (j ==
                       protein_triangle_ix1[indiv_id * max_protein_ + l])
                phenotype_activ[indiv_id * 300 + j] +=
                    height;
            }
          } else {
#ifdef __SIMD
            #pragma omp simd
#endif
            for (int j = 0; j < 300; j++) {
              if (j > protein_triangle_ix0[indiv_id * max_protein_ + l]
                  and
                  j < protein_triangle_ix1[indiv_id * max_protein_ + l])
                phenotype_inhib[indiv_id * 300 + j] += incY * (j -
                                                               protein_triangle_ix0[
                                                                   indiv_id *
                                                                   max_protein_ +
                                                                   l]);
              else if (
                  j > protein_triangle_ix1[indiv_id * max_protein_ + l]
                  and
                  j < protein_triangle_ix2[indiv_id * max_protein_ + l])
                phenotype_inhib[indiv_id * 300 + j] +=
                    height -
                    incY * (j -
                            protein_triangle_ix1[indiv_id * max_protein_ +
                                                 l]);
              else if (j ==
                       protein_triangle_ix1[indiv_id * max_protein_ + l])
                phenotype_inhib[indiv_id * 300 + j] +=
                    height;
            }
          }
        }
      }

      #pragma omp parallel for
      for (int indiv_id = 0; indiv_id < 1024; indiv_id++) {
#ifdef __SIMD
        #pragma omp simd
#endif
        for (int j = 0; j < 300; j++) {
          phenotype_activ[indiv_id * 300 + j] =
              phenotype_activ[indiv_id * 300 + j] > Y_MAX ? Y_MAX :
              phenotype_activ[indiv_id * 300 + j];
          phenotype_inhib[indiv_id * 300 + j] =
              phenotype_inhib[indiv_id * 300 + j] < -Y_MAX ? -Y_MAX :
              phenotype_inhib[indiv_id * 300 + j];
          phenotype[indiv_id * 300 + j] =
              phenotype_activ[indiv_id * 300 + j] +
              phenotype_inhib[indiv_id * 300 + j];
          phenotype[indiv_id * 300 + j] =
              phenotype[indiv_id * 300 + j] > Y_MIN ?
              Y_MIN : phenotype[indiv_id * 300 + j];
          delta[indiv_id * 300 + j] =
              phenotype[indiv_id * 300 + j] -
              environment[i * 300 + j];
        }
      }

      /** compute distance to target **/

      #pragma omp parallel for
      for (int indiv_id = 0; indiv_id < 1024; indiv_id++) {

        float area = 0;
#ifdef __SIMD
        #pragma omp simd
#endif
        for (int j = 0; j < 299; j++) {
          dist_sum[indiv_id] += ((delta[indiv_id * 300 + j]
                                  + delta[indiv_id * 300 + j + 1]) / 600.0);
        }
      }
    }
  }

  #pragma omp parallel for
  for (int indiv_id = 0; indiv_id < 1024; indiv_id++)
    dist_sum[indiv_id] = exp(
          selection_pressure_ * (dist_sum[indiv_id] / nb_eval_));

}

/*
void cuda_struct::print_dist(ExpManager* exp_m) {
  std::set<int>* eval = exp_m->exp_s()->get_list_eval_step();
  const Habitat_R& habitat = dynamic_cast<const Habitat_R&>(exp_m->world()->grid(0,0)->habitat());

  int degradation_step = exp_m->exp_s()->get_nb_degradation_step();
  float degradation_rate = exp_m->exp_s()->get_degradation_rate();
  float hill_shape_n = exp_m->exp_s()->get_hill_shape_n();
  float hill_shape = exp_m->exp_s()->get_hill_shape();

  Individual_R* indiv = dynamic_cast<Individual_R*>(exp_m->best_indiv());
  int indiv_id = indiv->grid_cell_->x() * 32 +indiv->grid_cell_->y();

  indiv->_initial_protein_list = indiv->protein_list_;

  //printf("MAX %d %d (%d %d)\n",max_protein_,max_rna_, indiv->grid_cell_->x() ,indiv->grid_cell_->y());

  float dist_temp = 0;

  //_protein_list.insert(_protein_list.end(), habitat.signals().begin(), habitat.signals().end());
  for(Protein_R* prot : habitat.signals()) {
    indiv->protein_list_.push_back(prot);
  }
  for (const auto& prot : indiv->protein_list_) {
    ((Protein_R*)prot)->reset_concentration();
  }

  indiv->set_influences();

  for (int j = 0; j < 300; j++) {
    phenotype[indiv_id * 300 + j] = 0;
    phenotype_activ[indiv_id * 300 + j] = 0;
    phenotype_inhib[indiv_id * 300 + j] = 0;
    delta[indiv_id * 300 + j] = 0;
  }

  for (int je = 0; je < max_protein_; je++) {
    protein_concentration[indiv_id * max_protein_ + je] = 0;

    protein_triangle_ix0[indiv_id * max_protein_ + je] = 0;
    protein_triangle_ix1[indiv_id * max_protein_ + je] = 0;
    protein_triangle_ix2[indiv_id * max_protein_ + je] = 0;

    protein_triangle_height[indiv_id * max_protein_ + je] = 0;


    for (int k = 0; k < max_rna_; k++)
      protein_influenced[indiv_id * max_protein_ * max_rna_ + je * max_protein_ + k] = 0;
  }

  for (auto prot_a : indiv->protein_list_) {
    Protein_R* prot = dynamic_cast<Protein_R*>(prot_a);

    prot->reset_concentration();
    protein_concentration[indiv_id * max_protein_ +
                          prot->get_local_id()] = prot->concentration();

    float x0 = prot->mean() - prot->width();
    float x1 = prot->mean();
    float x2 = prot->mean() + prot->width();

    int ix0 = (int) (x0 * 300);
    int ix1 = (int) (x1 * 300);
    int ix2 = (int) (x2 * 300);

    if (ix0 < 0) ix0 = 0; else if (ix0 > (300 - 1)) ix0 = 300 - 1;
    if (ix1 < 0) ix1 = 0; else if (ix1 > (300 - 1)) ix1 = 300 - 1;
    if (ix2 < 0) ix2 = 0; else if (ix2 > (300 - 1)) ix2 = 300 - 1;

    protein_triangle_ix0[indiv_id * max_protein_ + prot->get_local_id()] = ix0;
    protein_triangle_ix1[indiv_id * max_protein_ + prot->get_local_id()] = ix1;
    protein_triangle_ix2[indiv_id * max_protein_ + prot->get_local_id()] = ix2;

    protein_triangle_height[indiv_id * max_protein_ +
                            prot->get_local_id()] = prot->height();


    for (auto rna : prot->_rna_R_list) {
      protein_influenced[indiv_id * max_protein_ * max_rna_ +
                         prot->get_local_id() * max_protein_
                         + rna->get_local_id()] = 1;
    }
  }

  for (int je = 0; je < max_rna_; je++) {
    for (int ke = 0; ke < max_influence_; ke++) {
      enhance_coef[indiv_id * max_rna_ * max_influence_ +
                   je * max_influence_ + ke] = 0;
      operate_coef[indiv_id * max_rna_ * max_influence_ +
                   je * max_influence_ + ke] = 0;
      protein_influence[indiv_id * max_rna_ * max_influence_ +
                        je * max_influence_ + ke] = 0;
    }
  }


  for (auto rna : indiv->_rna_list_coding) {

    rna_synthesis[indiv_id * max_rna_ +
                  rna->get_local_id()] = rna->get_synthesis_rate();

    basal_level[indiv_id * max_rna_ +
                rna->get_local_id()] = rna->basal_level();


    for (int k = 0; k < rna->_nb_influences; k++) {
      enhance_coef[indiv_id * max_rna_ * max_influence_ +
                   rna->get_local_id() * max_rna_ +
                   k] = rna->_enhancing_coef_list[k];
      operate_coef[indiv_id * max_rna_ * max_influence_ +
                   rna->get_local_id() * max_rna_ +
                   k] = rna->_operating_coef_list[k];
      if (rna->_protein_list[k]->is_signal())
        protein_influence[indiv_id * max_rna_ * max_influence_ +
                          rna->get_local_id() * max_rna_ +
                          k] =
            rna->_protein_list[k]->get_local_id() + max_protein_;
      else
        protein_influence[indiv_id * max_rna_ * max_influence_ +
                          rna->get_local_id() * max_rna_ +
                          k] = rna->_protein_list[k]->get_local_id();
    }
  }



  for (int8_t i = 1; i <= life_time_; i++) {
    for (int j = 0; j < 300; j++) {
      environment[(i - 1) * 300 +
                  j] = ((HybridFuzzy*) habitat.phenotypic_target(
          i).fuzzy())->points()[j];

    }

    for (int j = 0; j < degradation_step; j++) {
      // Compute synthesis rate of RNA
      for (int m = 0; m < max_rna_; m++) {
        float enhancer_activity = 0, operator_activity = 0;
        for (int n = 0; n < max_influence_; n++) {
          int prot_id = protein_influence
          [indiv_id * max_rna_ * max_influence_ + m * max_influence_ +
           n];

          if (prot_id > max_protein_) {
            enhancer_activity +=
                enhance_coef[indiv_id * max_rna_ * max_influence_ +
                             m * max_influence_ + n] *
                signals_concentration[i * nb_signals_ + prot_id -
                                      max_protein_];

            operator_activity +=
                operate_coef[indiv_id * max_rna_ * max_influence_ +
                             m * max_influence_ + n] *
                signals_concentration[i * nb_signals_ + prot_id -
                                      max_protein_];
          } else {
            enhancer_activity +=
                enhance_coef[indiv_id * max_rna_ * max_influence_ +
                             m * max_influence_ + n] *
                protein_concentration[indiv_id * max_protein_ + prot_id];

            operator_activity +=
                operate_coef[indiv_id * max_rna_ * max_influence_ +
                             m * max_influence_ + n] *
                protein_concentration[indiv_id * max_protein_ + prot_id];
          }
        }

        float enhancer_activity_pow_n = pow(enhancer_activity,
                                            hill_shape_n);
        float operator_activity_pow_n = pow(operator_activity,
                                            hill_shape_n);

        rna_synthesis[indiv_id * max_rna_ + m] =
            basal_level[indiv_id * max_rna_ + m]
            * (hill_shape
               /
               (operator_activity_pow_n +
                hill_shape))
            * (1 + ((1 /
                     basal_level[indiv_id * max_rna_ + m]) -
                    1)
                   *
                   (enhancer_activity_pow_n /
                    (enhancer_activity_pow_n +
                     hill_shape)));
      }

      // Compute concentration
      for (int l = 0; l < max_protein_; l++) {
        float _delta_concentration = 0;
        for (int m = 0; m < max_rna_; m++) {
          if (protein_influenced[indiv_id * max_protein_ * max_rna_ +
                                 l * max_protein_
                                 + m] == 1) {
             _delta_concentration += rna_synthesis[indiv_id * max_rna_ +
                                                   m];
          }
        }

        _delta_concentration -= degradation_rate *
                                protein_concentration[indiv_id * max_protein_ +
                                                      l];
      _delta_concentration *= 1 / ((float) degradation_step);

        protein_concentration[indiv_id * max_protein_ +
                              l] += _delta_concentration;
      }

      indiv->update_concentrations();
    }

    for (auto rna : indiv->_rna_list_coding) {

    }

    for (auto prot : indiv->protein_list_)

      indiv->update_phenotype();
// If we have to evaluate the individual at this age
      if (eval_step[i]) {
        for (int j = 0; j < 300; j++) {
          phenotype[indiv_id * 300 + j] = 0;
          phenotype_activ[indiv_id * 300 + j] = 0;
          phenotype_inhib[indiv_id * 300 + j] = 0;
        }

        for (int l = 0; l < max_protein_; l++) {
          float height = (protein_triangle_height[indiv_id * max_protein_ +
                                                  l] *
                          protein_concentration[indiv_id * max_protein_ + l]);
          float incY =
              protein_triangle_height[indiv_id * max_protein_ + l] *
              protein_concentration[indiv_id * max_protein_ + l] /
              (protein_triangle_ix1[indiv_id * max_protein_ + l] -
               protein_triangle_ix0[indiv_id * max_protein_ + l]);


          if (protein_triangle_height[indiv_id * max_protein_ + l] > 0) {
            for (int j = 0; j < 300; j++) {
              if (j > protein_triangle_ix0[indiv_id * max_protein_ + l]
                  and
                  j < protein_triangle_ix1[indiv_id * max_protein_ + l])
                phenotype_activ[indiv_id * 300 + j] +=
                    incY * (j -
                            protein_triangle_ix0[indiv_id * max_protein_ +
                                                 l]);
              else if (
                  j > protein_triangle_ix1[indiv_id * max_protein_ + l]
                  and
                  j < protein_triangle_ix2[indiv_id * max_protein_ + l])
                phenotype_activ[indiv_id * 300 + j] += height
                                                       -
                                                       incY * (j -
                                                               protein_triangle_ix1[
                                                                   indiv_id *
                                                                   max_protein_ +
                                                                   l]);
              else if (j ==
                       protein_triangle_ix1[indiv_id * max_protein_ + l])
                phenotype_activ[indiv_id * 300 + j] +=
                    height;
            }
          } else {
            for (int j = 0; j < 300; j++) {
              if (j > protein_triangle_ix0[indiv_id * max_protein_ + l]
                  and
                  j < protein_triangle_ix1[indiv_id * max_protein_ + l])
                phenotype_inhib[indiv_id * 300 + j] += incY * (j -
                                                               protein_triangle_ix0[
                                                                   indiv_id *
                                                                   max_protein_ +
                                                                   l]);
              else if (
                  j > protein_triangle_ix1[indiv_id * max_protein_ + l]
                  and
                  j < protein_triangle_ix2[indiv_id * max_protein_ + l])
                phenotype_inhib[indiv_id * 300 + j] +=
                    height -
                    incY * (j -
                            protein_triangle_ix1[indiv_id * max_protein_ +
                                                 l]);
              else if (j ==
                       protein_triangle_ix1[indiv_id * max_protein_ + l])
                phenotype_inhib[indiv_id * 300 + j] +=
                    height;
            }
          }
        }

        for (int j = 0; j < 300; j++) {
          phenotype_activ[indiv_id * 300 + j] =
              phenotype_activ[indiv_id * 300 + j] > Y_MAX ? Y_MAX :
              phenotype_activ[indiv_id * 300 + j];
          phenotype_inhib[indiv_id * 300 + j] =
              phenotype_inhib[indiv_id * 300 + j] < -Y_MAX ? -Y_MAX :
              phenotype_inhib[indiv_id * 300 + j];
          phenotype[indiv_id * 300 + j] =
              phenotype_activ[indiv_id * 300 + j] +
              phenotype_inhib[indiv_id * 300 + j];
          phenotype[indiv_id * 300 + j] =
              phenotype[indiv_id * 300 + j] < Y_MIN ?
              Y_MIN : phenotype[indiv_id * 300 + j];
          delta[indiv_id * 300 + j] =
              phenotype[indiv_id * 300 + j] -
              environment[(i-1) * 300 + j];


          //HybridFuzzy* phen = new HybridFuzzy(*((HybridFuzzy*)indiv->phenotype_));
          //AbstractFuzzy* delta = FuzzyFactory::fuzzyFactory->create_fuzzy();
          //phen->sub(*(((HybridFuzzy*) habitat.phenotypic_target(i).fuzzy())));
          //HybridFuzzy* del = (HybridFuzzy*)delta;

        }


        float area = 0;

        for (int j = 0; j < 299; j++) {
          dist_sum[indiv_id] += ((fabs(delta[indiv_id * 300 + j])
                                  + fabs(delta[indiv_id * 300 + j + 1])) / 600.0);
          area += ((fabs(delta[indiv_id * 300 + j])
                    + fabs(delta[indiv_id * 300 + j + 1])) / 600.0);
        }

        for (int a = 0; a < NB_FEATURES; a++) {
          indiv->dist_to_target_by_feature_[a] = 0;
        }

        indiv->distance_to_target_computed_ = false;
        indiv->phenotype_computed_ = true;

        indiv->compute_distance_to_target(habitat.phenotypic_target(i));

        dist_temp += indiv->dist_to_target_by_feature_[METABOLISM];

        //printf("Dist sum at %d : %f %f %f\n",i,dist_sum[indiv_id],indiv->dist_to_target_by_feature_[METABOLISM],area);
      }
  }

  indiv->dist_to_target_by_feature_[METABOLISM] = dist_temp / (float) (nb_eval_);

  printf("Mean dist_to_target =  %e %e\n",dist_temp/(float)nb_eval_,dist_sum[indiv_id]/(float)nb_eval_);
  //printf("fitness CUDA!!! %e %e\n",selection_pressure_, dist_sum[indiv_id]/(float)nb_eval_);

    dist_sum[indiv_id] = exp(
        -selection_pressure_ * (dist_sum[indiv_id] / (float)nb_eval_));
  indiv->fitness_computed_ = false;

  indiv->compute_fitness(habitat.phenotypic_target( life_time_ ));

  printf("Fitness =  %e %e\n",indiv->fitness(),dist_sum[indiv_id]);

  indiv->protein_list_.clear();
  indiv->protein_list_ = indiv->_initial_protein_list;

}
*/

void cuda_struct::compute_a_generation_v3(ExpManager* exp_m) {

  int degradation_step = exp_m->exp_s()->get_nb_degradation_step();
  float degradation_rate = exp_m->exp_s()->get_degradation_rate();
  float hill_shape_n = exp_m->exp_s()->get_hill_shape_n();
  float hill_shape = exp_m->exp_s()->get_hill_shape();

#pragma omp parallel num_threads(4)
  {
#pragma omp single nowait
    {
      for (int ki = 0; ki < 32; ki++)
        for (int kj = 0; kj < 32; kj++) {
          int indiv_id = ki * 32 + kj;

          for (int8_t i = 0; i < life_time_; i++) {
            for (int j = 0; j < degradation_step; j++) {
              // Compute synthesis rate of RNA
              for (int m = 0; m < max_rna_; m++) {

/*#pragma omp task \
                  depend(in: protein_influence[indiv_id * max_rna_ * max_influence_ + m * max_influence_:max_influence_]) \
                  depend(in: enhance_coef[indiv_id * max_rna_ * max_influence_ + m * max_influence_:max_influence_]) \
                  depend(in: enhance_coef[indiv_id * max_rna_ * max_influence_ + m * max_influence_:max_influence_]) \
                  depend(in: signals_concentration[i*nb_signals_:nb_signals_]) \
                  depend(in: protein_concentration[indiv_id * max_protein_:max_protein_]) \
                  depend(out: rna_synthesis[indiv_id * max_rna_ + m])*/
                {
                  float enhancer_activity = 0, operator_activity = 0;
                  for (int n = 0; n < max_influence_; n++) {
                    int prot_id = protein_influence
                    [ki * 32 + kj * max_rna_ * max_influence_ +
                     m * max_influence_ +
                     n];

                    if (prot_id > max_protein_) {
                      enhancer_activity +=
                          enhance_coef[indiv_id * max_rna_ * max_influence_ +
                                       m * max_influence_ + n] *
                          signals_concentration[i * nb_signals_ + prot_id -
                                                max_protein_];

                      operator_activity +=
                          operate_coef[indiv_id * max_rna_ * max_influence_ +
                                       m * max_influence_ + n] *
                          signals_concentration[i * nb_signals_ + prot_id -
                                                max_protein_];
                    } else {
                      enhancer_activity +=
                          enhance_coef[indiv_id * max_rna_ * max_influence_ +
                                       m * max_influence_ + n] *
                          protein_concentration[indiv_id * max_protein_ +
                                                prot_id];

                      operator_activity +=
                          operate_coef[indiv_id * max_rna_ * max_influence_ +
                                       m * max_influence_ + n] *
                          protein_concentration[indiv_id * max_protein_ +
                                                prot_id];
                    }
                  }

                  float enhancer_activity_pow_n = pow(enhancer_activity,
                                                      hill_shape_n);
                  float operator_activity_pow_n = pow(operator_activity,
                                                      hill_shape_n);

                  rna_synthesis[indiv_id * max_rna_ + m] =
                      basal_level[indiv_id * max_rna_ + m]
                      * (hill_shape
                         /
                         (operator_activity_pow_n +
                          hill_shape))
                      * (1 + ((1 /
                               basal_level[indiv_id * max_rna_ + m]) -
                              1)
                             *
                             (enhancer_activity_pow_n /
                              (enhancer_activity_pow_n +
                               hill_shape)));
                }
              }

              // Compute concentration
              for (int l = 0; l < max_protein_; l++) {
/*#pragma omp task \
                  depend(in: protein_influenced[indiv_id * max_protein_ * max_rna_ +l * max_protein_:max_rna_]) \
                  depend(in: rna_synthesis[indiv_id * max_rna_:max_rna_]) \
                  depend(inout: protein_concentration[indiv_id * max_protein_ +l])*/
                {

                  float _delta_concentration;

                  for (int m = 0; m < max_rna_; m++) {
                    if (protein_influenced[indiv_id * max_protein_ * max_rna_ +
                                           l * max_protein_
                                           + m])
                      _delta_concentration += rna_synthesis[
                          indiv_id * max_rna_ +
                          m];
                  }

                  _delta_concentration -=
                      degradation_rate * protein_concentration[
                          indiv_id * max_protein_ +
                          l];
                  _delta_concentration *= 1 / ((float) degradation_step);

                  protein_concentration[indiv_id * max_protein_ +
                                        l] += _delta_concentration;
                }
              }
            }

            // If we have to evaluate the individual at this age
            if (eval_step[i]) {
              /** update phenotype **/
#pragma omp taskwait

/*#pragma omp task depend(out: phenotype[indiv_id * 300:300]) \
                           depend(out: phenotype_activ[indiv_id * 300:300]) \
                           depend(out: phenotype_inhib[indiv_id * 300:300])*/
              {
#ifdef __SIMD
#pragma omp simd
#endif
                for (int j = 0; j < 300; j++) {
                  phenotype[indiv_id * 300 + j] = 0;
                  phenotype_activ[indiv_id * 300 + j] = 0;
                  phenotype_inhib[indiv_id * 300 + j] = 0;
                }
              }

              for (int l = 0; l < max_protein_; l++) {
/*#pragma omp task depend(in: protein_triangle_height[indiv_id * max_protein_:max_protein_]) \
                           depend(in: protein_concentration[indiv_id * max_protein_:max_protein_]) \
                           depend(in: protein_triangle_ix0[indiv_id * max_protein_:max_protein_]) \
                           depend(in: protein_triangle_ix1[indiv_id * max_protein_:max_protein_]) \
                           depend(in: protein_triangle_ix2[indiv_id * max_protein_:max_protein_]) \
                           depend(inout: phenotype_activ[indiv_id * 300:300]) \
                           depend(inout: phenotype_inhib[indiv_id * 300:300])*/
                {
                  float height = (
                      protein_triangle_height[indiv_id * max_protein_ +
                                              l] *
                      protein_concentration[indiv_id * max_protein_ +
                                            l]);
                  float incY =
                      protein_triangle_height[indiv_id * max_protein_ + l] *
                      protein_concentration[indiv_id * max_protein_ + l] /
                      (protein_triangle_ix1[indiv_id * max_protein_ + l] -
                       protein_triangle_ix0[indiv_id * max_protein_ + l]);

                  if (protein_triangle_height[indiv_id * max_protein_ + l] >
                      0) {
                    for (int j = 0; j < 300; j++) {
                      if (j > protein_triangle_ix0[indiv_id * max_protein_ + l]
                          and
                          j < protein_triangle_ix1[indiv_id * max_protein_ + l])
                        phenotype_activ[indiv_id * 300 + j] +=
                            incY * (j -
                                    protein_triangle_ix0[
                                        indiv_id * max_protein_ +
                                        l]);
                      else if (
                          j > protein_triangle_ix1[indiv_id * max_protein_ + l]
                          and
                          j < protein_triangle_ix2[indiv_id * max_protein_ + l])
                        phenotype_activ[indiv_id * 300 + j] += height
                                                               -
                                                               incY * (j -
                                                                       protein_triangle_ix1[
                                                                           indiv_id *
                                                                           max_protein_ +
                                                                           l]);
                      else if (j ==
                               protein_triangle_ix1[indiv_id * max_protein_ +
                                                    l])
                        phenotype_activ[indiv_id * 300 + j] +=
                            height;
                    }
                  } else {

                    for (int j = 0; j < 300; j++) {
                      if (j > protein_triangle_ix0[indiv_id * max_protein_ + l]
                          and
                          j < protein_triangle_ix1[indiv_id * max_protein_ + l])
                        phenotype_inhib[indiv_id * 300 + j] += incY * (j -
                                                                       protein_triangle_ix0[
                                                                           indiv_id *
                                                                           max_protein_ +
                                                                           l]);
                      else if (
                          j > protein_triangle_ix1[indiv_id * max_protein_ + l]
                          and
                          j < protein_triangle_ix2[indiv_id * max_protein_ + l])
                        phenotype_inhib[indiv_id * 300 + j] +=
                            height -
                            incY * (j -
                                    protein_triangle_ix1[
                                        indiv_id * max_protein_ +
                                        l]);
                      else if (j ==
                               protein_triangle_ix1[indiv_id * max_protein_ +
                                                    l])
                        phenotype_inhib[indiv_id * 300 + j] +=
                            height;
                    }
                  }
                }
              }


/*#pragma omp task depend(inout: phenotype[indiv_id * 300:300]) \
                           depend(inout: phenotype_activ[indiv_id * 300:300]) \
                           depend(inout: phenotype_inhib[indiv_id * 300:300]) \
                           depend(out: delta[indiv_id * 300:300]) \
                           depend(in: environment[indiv_id * 300:300])*/
              {
#ifdef __SIMD
#pragma omp simd
#endif
                for (int j = 0; j < 300; j++) {
                  phenotype_activ[indiv_id * 300 + j] =
                      phenotype_activ[indiv_id * 300 + j] > Y_MAX ? Y_MAX :
                      phenotype_activ[indiv_id * 300 + j];
                  phenotype_inhib[indiv_id * 300 + j] =
                      phenotype_inhib[indiv_id * 300 + j] < -Y_MAX ? -Y_MAX :
                      phenotype_inhib[indiv_id * 300 + j];
                  phenotype[indiv_id * 300 + j] =
                      phenotype_activ[indiv_id * 300 + j] +
                      phenotype_inhib[indiv_id * 300 + j];
                  phenotype[indiv_id * 300 + j] =
                      phenotype[indiv_id * 300 + j] > Y_MIN ?
                      Y_MIN : phenotype[indiv_id * 300 + j];
                  delta[indiv_id * 300 + j] =
                      phenotype[indiv_id * 300 + j] -
                      environment[i * 300 + j];
                }
              }
              /** compute distance to target **/


/*#pragma omp task depend(in: delta[indiv_id * 300:300]) \
                           depend(out: dist_sum[indiv_id])*/
              {
                for (int j = 0; j < 299; j++) {
                  dist_sum[indiv_id] += ((delta[indiv_id * 300 + j]
                                          + delta[indiv_id * 300 + j + 1]) /
                                         600.0);
                }
              }
            }
          }

//#pragma omp task depend(inout: dist_sum[indiv_id])
          {
            dist_sum[indiv_id] = exp(
                selection_pressure_ * (dist_sum[indiv_id] / nb_eval_));
          }
        }
    }
  }
}
}
