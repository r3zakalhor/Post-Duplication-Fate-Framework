// ****************************************************************************
//
//          Aevol - An in silico experimental evolution platform
//
// ****************************************************************************
//
// Copyright: See the AUTHORS file provided with the package or <www.aevol.fr>
// Web: http://www.aevol.fr/
// E-mail: See <http://www.aevol.fr/contact/>
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons, Jonathan Rouzaud-Cornabas
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

#include "SIMD_PhenotypicTargetHandler_R.h"
#include "Discrete_Double_Fuzzy.h"
#include "ExpSetup.h"
#include "Protein_7.h"

namespace aevol {

SIMD_PhenotypicTargetHandler_R::SIMD_PhenotypicTargetHandler_R(
    std::shared_ptr<PhenotypicTargetHandler_R> handler,
    ExpSetup* exp_s,
    FuzzyFactory_7* fuzzy_factory_,
    bool check_simd) {

  var_method_ = handler->var_method();

  if (!check_simd && var_method_ != NO_VAR)
    var_prng_ = std::make_shared<JumpingMT>(*(handler->var_prng_));

  if (var_method_ == NO_VAR)
    var_method_ = SWITCH_IN_A_LIST;

  env_gaussians_list_.resize(handler->env_gaussians_list_.size());
  nb_env_ = env_gaussians_list_.size();

  int i = 0;
  for (auto gaussian_env: handler->env_gaussians_list_) {
    for (auto gauss: gaussian_env) {
      env_gaussians_list_[i].emplace_back(gauss);
    }
    handler->env_signals_list_[i];
    i++;
  }

  env_signals_list_.resize(handler->env_signals_list_.size());
  i = 0;
  for (auto env_signal: handler->env_signals_list_) {
    for (auto e_signal: env_signal) {
      env_signals_list_[i].push_back(e_signal);
    }
    i++;
  }

  i = 0;
  signals_models_.resize(handler->signals_models_.size());
  for (auto signal_protein: handler->signals_models_) {
    Protein_7* prot    = new Protein_7(signal_protein);
    prot->protein_id_ = -i;
    signals_models_[i] = prot;
    i++;
  }

  env_switch_probability_ = handler->env_switch_probability_;
  nb_indiv_age_           = handler->_nb_indiv_age;

  //printf("VAR METHOD %d ! %d :: %d\n",check_simd, var_method_,handler->var_method());

  sampling_ = handler->sampling();

  targets_fuzzy_by_id_ = new AbstractFuzzy_7*[nb_env_];

  for (int env_id = 0; env_id < nb_env_; env_id++) {
    targets_fuzzy_by_id_[env_id] = fuzzy_factory_->get_fuzzy();

    if (not env_gaussians_list_.at(env_id).empty()) {
      for (int16_t i = 0; i < sampling_; i++) {
        Point new_point =
            Point(X_MIN + (double)i * (X_MAX - X_MIN) / (double)sampling_, 0.0);
        int gi = 0;
        for (const Gaussian& g: env_gaussians_list_.at(env_id)) {
          gi++;
          new_point.y += g.compute_y(new_point.x);
          // printf("SIMD -- %d -- Compute point %e %e\n",env_id,new_point.x,new_point.y);
        }
        // if (new_point.y == 0) printf("ERROR :: ");
        // printf("SIMD -- %d -- Add point %e %e\n",env_id,new_point.x,new_point.y);

        ((Discrete_Double_Fuzzy*)targets_fuzzy_by_id_[env_id])->add_point(i, new_point.y);
      }
    }
    // Add lower and upper bounds
    targets_fuzzy_by_id_[env_id]->clip(AbstractFuzzy_7::min, Y_MIN);
    targets_fuzzy_by_id_[env_id]->clip(AbstractFuzzy_7::max, Y_MAX);

    // Simplify (get rid of useless points)
    targets_fuzzy_by_id_[env_id]->simplify();

    // targets_fuzzy_by_id_[env_id]->print();
  }

  if (var_method_ == SWITCH_IN_A_LIST) {
    targets_fuzzy_ = new AbstractFuzzy_7*[nb_indiv_age_];
    list_env_id_   = new int16_t[nb_indiv_age_];

    if (nb_env_ <= 1) {
      for (int age = 0; age < nb_indiv_age_; age++) {
        targets_fuzzy_[age] = targets_fuzzy_by_id_[0];
        list_env_id_[age]   = 0;
      }
    } else {
      if (env_switch_probability_ <= 0.0) {
        int32_t duration = std::floor(nb_indiv_age_ / nb_env_);
        int32_t cur_env_id = 0, cur_duration = 0;
        hasChanged_ = false;
        
        printf("Duration per env %d\n",duration);
        for (int age = 0; age < nb_indiv_age_; age++) {
          targets_fuzzy_[age] =
              targets_fuzzy_by_id_[cur_env_id];
          list_env_id_[age] = cur_env_id;
          printf("Init Env at age %d is %d :: %lf // L %d\n",age,cur_env_id,
                targets_fuzzy_[age]->get_geometric_area(),list_env_id_[age]);
          cur_duration++;
          if ((cur_duration >= duration) && (cur_env_id + 1 < nb_env_)) {
            cur_duration = 0;
            cur_env_id += 1;
          }
        }
      } else {
        for (int age = 0; age < nb_indiv_age_; age++) {
          targets_fuzzy_[age] =
              targets_fuzzy_by_id_[handler->phenotypic_targets_[age]->get_id()];
          list_env_id_[age] = handler->phenotypic_targets_[age]->get_id();
          // printf("Init Env at age %d is %d :: %lf\n",age,handler->phenotypic_targets_[age]->get_id(),
          //       targets_fuzzy_[age]->get_geometric_area());
        }
        // targets_fuzzy_by_id_[1]->get_geometric_area(true);
      }
    }

    nb_eval_ = exp_s->get_list_eval_step()->size();
  }

  handler_    = handler;
  check_simd_ = check_simd;
}

SIMD_PhenotypicTargetHandler_R::SIMD_PhenotypicTargetHandler_R(
    SIMD_PhenotypicTargetHandler_R* handler, FuzzyFactory_7* fuzzy_factory_, ExpManager* exp_m) {

  var_method_ = handler->var_method_;

  //   if (!check_simd && var_method_ != NO_VAR)
  //     var_prng_ = handler->var_prng_;

  //   if (var_method_ == NO_VAR)
  //     var_method_ = SWITCH_IN_A_LIST;

  // env_gaussians_list_ = handler->env_gaussians_list_);
  nb_env_ = handler->nb_env_;

  //   int i = 0;
  //   for (auto gaussian_env : handler->env_gaussians_list_) {
  //     for (auto gauss : gaussian_env) {
  //       env_gaussians_list_[i].emplace_back(gauss);
  //     }
  //     handler->env_signals_list_[i];
  //     i++;
  //   }

  env_signals_list_ = handler->env_signals_list_;
  //   i=0;
  //   for (auto env_signal : handler->env_signals_list_) {
  //     for (auto e_signal : env_signal) {
  //       env_signals_list_[i].push_back(e_signal);
  //     }
  //     i++;

  //   }

  //   i = 0;
  
signals_models_.resize(handler->signals_models_.size());
int32_t i = 0;
    for (auto signal_protein : handler->signals_models_) {
      Protein_7* prot = new Protein_7(signal_protein,exp_m);
      prot->protein_id_ = -i;
      signals_models_[i] = prot;
      i++;
    }
  env_switch_probability_ = handler->env_switch_probability_;
  nb_indiv_age_           = handler->nb_indiv_age_;

  //   //printf("VAR METHOD %d ! %d :: %d\n",check_simd, var_method_,handler->var_method());

  sampling_ = handler->sampling_;

  targets_fuzzy_by_id_ = new AbstractFuzzy_7*[nb_env_];

  for (int env_id = 0; env_id < nb_env_; env_id++) {
    targets_fuzzy_by_id_[env_id] = fuzzy_factory_->get_fuzzy();
    targets_fuzzy_by_id_[env_id]->copy(handler->targets_fuzzy_by_id_[env_id]);
  }

  if (var_method_ == SWITCH_IN_A_LIST) {
    targets_fuzzy_ = new AbstractFuzzy_7*[nb_indiv_age_];
    list_env_id_   = new int16_t[nb_indiv_age_];


      if (env_switch_probability_ == 0.0) {
        int32_t duration = std::floor(nb_indiv_age_ / nb_env_);
        int32_t cur_env_id = 0, cur_duration = 0;
        hasChanged_ = false;
        
        printf("Duration per env %d\n",duration);
        for (int age = 0; age < nb_indiv_age_; age++) {
          targets_fuzzy_[age] =
              targets_fuzzy_by_id_[cur_env_id];
          list_env_id_[age] = cur_env_id;
          printf("Init Env at age %d is %d :: %lf // L %d\n",age,cur_env_id,
                targets_fuzzy_[age]->get_geometric_area(),list_env_id_[age]);
          cur_duration++;
          if ((cur_duration >= duration) && (cur_env_id + 1 < nb_env_)) {
            cur_duration = 0;
            cur_env_id += 1;
          }
        }
      } else {
        for (int age = 0; age < nb_indiv_age_; age++) {
          targets_fuzzy_[age] = targets_fuzzy_by_id_[handler->list_env_id_[age]];
          list_env_id_[age]   = handler->list_env_id_[age];
        }
      }

    nb_eval_ = handler->nb_eval_;
  }

  handler_    = handler->handler_;
  check_simd_ = handler->check_simd_;
}

SIMD_PhenotypicTargetHandler_R::~SIMD_PhenotypicTargetHandler_R() {
  
  if (var_method_ == SWITCH_IN_A_LIST) {
    delete [] targets_fuzzy_;
    delete [] list_env_id_ ;
  }
  
  for (int i = 0; i < nb_env_; i++)
    delete targets_fuzzy_by_id_[i];

  delete[] targets_fuzzy_by_id_;

  for (std::vector<Protein_7*>::iterator it_protein = signals_models_.begin();
       it_protein != signals_models_.end(); it_protein++) {
    delete (*(it_protein));
  }

  signals_models_.clear();

  int i = 0;
  for (auto &gaussian_env: env_gaussians_list_) {
    gaussian_env.clear();
    i++;
  }

  env_gaussians_list_.clear();
for (auto &env_signals: env_signals_list_) {
    env_signals.clear();
    i++;
  }
  env_signals_list_.clear();

}

void SIMD_PhenotypicTargetHandler_R::ApplyVariation() {

  // printf("ApplyVar\n");
  switch (var_method_) {
  case NO_VAR:
    return;
  case SWITCH_IN_A_LIST: {
    if (nb_env_ <= 1) {
      break;
    }

    int16_t* list_of_old_target_id;
    int16_t id_new_env;
    int16_t id_old_env;

    if (env_switch_probability_ > 0.0) {
      list_of_old_target_id = list_env_id_;
      list_env_id_                   = new int16_t[nb_indiv_age_];

      // Shortcuts used
      id_new_env = list_of_old_target_id[nb_indiv_age_ - 1];
      id_old_env = list_of_old_target_id[nb_indiv_age_ - 1];
    }
    hasChanged_ = false;

    if (check_simd_) {
      hasChanged_ = handler_->hasChanged();
      for (int age = 0; age < nb_indiv_age_; age++) {
        targets_fuzzy_[age] =
            targets_fuzzy_by_id_[handler_->phenotypic_targets_[age]->get_id()];
        list_env_id_[age] = handler_->phenotypic_targets_[age]->get_id();
      }
    } else {
      if (env_switch_probability_ > 0.0) {
        for (int16_t i = 0; i < nb_indiv_age_ ; i++) {
            // if we have to change of environment :
            double env_chang = var_prng_->random();

            if (env_chang < env_switch_probability_) {
              //we have to change to a new env that have an id different from the old one
              while (id_new_env == id_old_env) {
                id_new_env = var_prng_->random(nb_env_);
              }
                    //The environment has changed
              id_old_env = id_new_env;
            }

      
            list_env_id_[i] = id_new_env;
            targets_fuzzy_[i] = targets_fuzzy_by_id_[id_new_env];

    // printf("Env at %d : %d (OLD %d)\n",i,list_env_id_[i],list_of_old_target_id[i]);
  
            if (list_env_id_[i] != list_of_old_target_id[i])
              hasChanged_ = true;
        }
        delete [] list_of_old_target_id;
      }

    }
    
    

    break; }
  case ONE_AFTER_ANOTHER:
    break;
  default:
    Utils::ExitWithDevMsg("Unknown variation method", __FILE__, __LINE__);
    break;
  }
}

void SIMD_PhenotypicTargetHandler_R::set_single_env(int16_t id) {
	printf("Set %d for %d age\n",id,nb_indiv_age_);

  for (int16_t i = 0; i < nb_indiv_age_; i++) {
    list_env_id_[i]   = id;
    targets_fuzzy_[i] = targets_fuzzy_by_id_[id];
  }
}
}  // namespace aevol
