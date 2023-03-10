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

#ifndef AEVOL_SIMD_PHENOTYPICTARGETHANDLER_R_H
#define AEVOL_SIMD_PHENOTYPICTARGETHANDLER_R_H

#include "population/SIMD_Individual.h"
#include "PhenotypicTargetHandler_R.h"
#include "phenotype/Gaussian.h"
#include "raevol/Protein_R.h"
#include "AbstractFuzzy_7.h"
#include "FuzzyFactory_7.h"
#include <list>
#include <vector>

namespace aevol {
class Protein_7;

class SIMD_PhenotypicTargetHandler_R {
 public:
  SIMD_PhenotypicTargetHandler_R(std::shared_ptr<PhenotypicTargetHandler_R> handler, 
    ExpSetup* exp_s,  FuzzyFactory_7* fuzzy_factory, bool check_simd = false);

    SIMD_PhenotypicTargetHandler_R(SIMD_PhenotypicTargetHandler_R* handler, FuzzyFactory_7* fuzzy_factory_, ExpManager* exp_m);

    ~SIMD_PhenotypicTargetHandler_R();

  void ApplyVariation();
  
  // For post_treatment evaluate_regulation
  void set_single_env(int16_t id);

  std::vector<Protein_7*> signals_models_;
  std::vector<std::list<int16_t>> env_signals_list_;

  AbstractFuzzy_7** targets_fuzzy_;
  AbstractFuzzy_7** targets_fuzzy_by_id_;

  PhenotypicTargetVariationMethod var_method_;

  int16_t* list_env_id_;

  int16_t nb_indiv_age_;
  int16_t nb_eval_;
  int16_t nb_env_;
  bool hasChanged_ = false;

  std::shared_ptr<JumpingMT> var_prng_;

  double env_switch_probability_;

 protected:
  std::vector<std::list<Gaussian>> env_gaussians_list_;

  std::shared_ptr<PhenotypicTargetHandler_R> handler_;



  int16_t sampling_;



  bool check_simd_ = false;
};

}

#endif //AEVOL_SIMD_PHENOTYPICTARGETHANDLER_R_H
