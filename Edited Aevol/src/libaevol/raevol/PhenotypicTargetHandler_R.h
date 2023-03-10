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
//*****************************************************************************


#ifndef AEVOL_PHENOTYPIC_TARGET_HANDLER_R_H__
#define AEVOL_PHENOTYPIC_TARGET_HANDLER_R_H__


// ============================================================================
//                                   Includes
// ============================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <memory>
//#include <list>
#include <vector>

#include "PhenotypicTargetHandler.h"
#include "PhenotypicTarget_R.h"
#include "Utils.h"
//#include "Habitat_R.h"

//using std::list;


namespace aevol {

// ============================================================================
//                          Class declarations
// ============================================================================
class Habitat_R;
class ExpSetup;



/**
 * Manages a phenotypic target and its "evolution" over time
 *
 * Handles a phenotypic target, the variation and/or noise that may be applied
 * to it as well as the set of possible phenotypic targets and the rules that
 * define how and when we switch from one to another
 */
class PhenotypicTargetHandler_R : public virtual PhenotypicTargetHandler
{
 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  PhenotypicTargetHandler_R(void); //< Default ctor
  PhenotypicTargetHandler_R(const PhenotypicTargetHandler_R&); //< Copy ctor
  PhenotypicTargetHandler_R(PhenotypicTargetHandler_R&&) = delete; //< Move ctor
  PhenotypicTargetHandler_R(gzFile backup_file);

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~PhenotypicTargetHandler_R(void); //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  virtual void ApplyVariation();
  void InitPhenotypicTargetsAndModels(int16_t nb_indiv_age);
  void print_geometric_areas();
  virtual void save(gzFile backup_file) const;
  virtual void load(gzFile backup_file);

  // ==========================================================================
  //                                 Getters
  // ==========================================================================
  const std::vector<PhenotypicTarget_R*>& phenotypic_targets() const {
    return phenotypic_targets_;
  }

  const PhenotypicTarget_R& phenotypic_target_model(int16_t env_id) const {
    assert(env_id >= 0 && env_id <= (int16_t) phenotypic_target_models_.size());
    return *(phenotypic_target_models_.at(env_id));
  }

  const PhenotypicTarget_R& phenotypic_target(int16_t age) const {
    //printf("AX --- %d %d\n",age,phenotypic_targets_.size());

    //assert(age >= 0 && age <= (int16_t) phenotypic_targets_.size());
    return *(phenotypic_targets_.at(age-1));
  }

  int16_t number_of_phenotypic_targets() const {
    return phenotypic_targets_.size();
  }


  int16_t number_of_phenotypic_target_models() const {
    return phenotypic_target_models_.size();
  }

  virtual double mean_environmental_area(ExpSetup* exp_s) const;

  const std::list<Protein_R*> signals() const {
    return signals_models_list_;
  }

  // ==========================================================================
  //                                 Setters
  // ==========================================================================
  void set_gaussians(const std::vector<std::list<Gaussian>>& gaussians_list) {
    env_gaussians_list_ = gaussians_list;
  }

  void set_signals(const std::vector<std::list<int16_t>>& signals_list) {
    env_signals_list_ = signals_list;
  }

  void set_signals_models(const std::vector<Protein_R*>& signals_list) {
    signals_models_ = signals_list;
    std::list<Protein_R*> temp_list;
    int local_id = 0;
    for(Protein_R* prot : signals_list) {
      temp_list.push_back(prot);
      prot->set_local_id(local_id++);
    }
    signals_models_list_ = temp_list;
    for(Protein_R* prot : signals_models_list_) {
      printf("Signals prot id %d\n",prot->get_local_id());
    }
  }

  void set_switch_probability(double p) {
    env_switch_probability_ = p;
  }

  virtual void set_segmentation(int16_t nb_segments,
                        double* boundaries,
                        PhenotypicFeature * features,
                        bool separate_segments = false) {
    for(PhenotypicTarget_R* phenotypic_target : phenotypic_target_models_) {
      phenotypic_target->set_segmentation(nb_segments,
                                         boundaries,
                                         features,
                                         separate_segments);
    }
  }

  bool hasChanged() { return hasChanged_; }
  std::vector<PhenotypicTarget_R*> phenotypic_target_models_;

    // For post_treatment evaluate_regulation
    void set_single_env(int16_t id);
    void set_two_env(int16_t id_1, int16_t id_2);

    void ShuffleRandomlySignals();

  std::vector<PhenotypicTarget_R*> phenotypic_targets_;
  std::vector<std::list<Gaussian>> env_gaussians_list_;
  std::vector<std::list<int16_t>> env_signals_list_;
  std::vector<Protein_R*> signals_models_;
  std::list<Protein_R*> signals_models_list_;
  double env_switch_probability_;
  int16_t _nb_indiv_age;

 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================
  void InitPhenotypicTargetsModels();
  void BuildPhenotypicTargetsModels();
  void BuildPhenotypicTargetModel( int16_t id);
  // This function keep only the last element of the vector
  void ResetPhenotypicTargets();
  void InitPhenotypicTargets(int16_t nb_indiv_age);
  void addEnv( int time, int16_t env_id );
  void changeEnv( int16_t ind, int16_t env_id );


  // ==========================================================================
  //                               Attributes
  // ==========================================================================


  bool hasChanged_;

  bool init_2_env = false;

  bool deterministic = false;

  bool long_period = true;

};

// ============================================================================
//                       Inline functions' definition
// ============================================================================

} // namespace aevol

#endif // AEVOL_PHENOTYPIC_TARGET_HANDLER_R_H__
