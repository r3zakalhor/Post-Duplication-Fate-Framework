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


#ifndef AEVOL_PHENOTYPIC_TARGET_HANDLER_H_
#define AEVOL_PHENOTYPIC_TARGET_HANDLER_H_


// ============================================================================
//                                   Includes
// ============================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <memory>
#include <list>

#include "PhenotypicTarget.h"
#include "Gaussian.h"
#include "ae_enums.h"
#include "JumpingMT.h"
#include "AbstractFuzzy.h"

using std::list;


namespace aevol {

// ============================================================================
//                          Class declarations
// ============================================================================
class ExpSetup;



/**
 * Manages a phenotypic target and its "evolution" over time
 *
 * Handles a phenotypic target, the variation and/or noise that may be applied
 * to it as well as the set of possible phenotypic targets and the rules that
 * define how and when we switch from one to another
 */
class PhenotypicTargetHandler
{
 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  PhenotypicTargetHandler(); //< Default ctor
  PhenotypicTargetHandler(const PhenotypicTargetHandler&); //< Copy ctor
  PhenotypicTargetHandler(PhenotypicTargetHandler&&) = delete; //< Move ctor
  PhenotypicTargetHandler(gzFile backup_file);

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~PhenotypicTargetHandler(); //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  void BuildPhenotypicTarget();
  virtual void ApplyVariation();

  virtual void save(gzFile backup_file) const;
  virtual void load(gzFile backup_file);

  // ==========================================================================
  //                                 Getters
  // ==========================================================================
  const PhenotypicTarget& phenotypic_target() const {
    return *phenotypic_target_;
  }
  double get_geometric_area() const {
    return phenotypic_target_->fuzzy()->get_geometric_area();
  }
  double area_by_feature(int8_t feature) const {
    return phenotypic_target_->area_by_feature(feature);
  }
  const list<Gaussian>& gaussians() const {
    return initial_gaussians_;
  }

    const list<Gaussian>& current_gaussians() const {
      return current_gaussians_;
    }

  const PhenotypicTargetVariationMethod& var_method() const {
    return var_method_;
  }

  virtual double mean_environmental_area(ExpSetup* exp_s = nullptr) const {
    return phenotypic_target_->area_by_feature(METABOLISM);
  }

  int16_t sampling() const {
      return sampling_;
  }

  double noise_alpha() const {
      return noise_alpha_;
  }

  double noise_prob() const {
      return noise_prob_;
  }

  int8_t noise_sampling_log() const {
      return noise_sampling_log_;
  }

  double noise_sigma() const {
      return noise_sigma_;
  }

  double var_sigma() const {
      return var_sigma_;
  }

  double var_tau() const {
      return var_tau_;
  }

  std::shared_ptr<JumpingMT> noise_prng() const {
      return noise_prng_;
  }

  std::shared_ptr<JumpingMT> var_prng() const {
      return var_prng_;
  }

  // ==========================================================================
  //                                 Setters
  // ==========================================================================
  void set_gaussians(const list<Gaussian>& gaussians) {
    current_gaussians_ = initial_gaussians_ = gaussians;
  }
  void set_sampling(int16_t val){
    sampling_ = val;
  }
  virtual void set_segmentation(int16_t nb_segments,
                        double* boundaries,
                        PhenotypicFeature * features,
                        bool separate_segments = false) {
    phenotypic_target_->set_segmentation(nb_segments,
                                         boundaries,
                                         features,
                                         separate_segments);
  }
  void set_var_method(PhenotypicTargetVariationMethod var_method) {
    var_method_ = var_method;
  }
  void set_var_prng(std::shared_ptr<JumpingMT> prng) {
    var_prng_ = prng;
  }
  void set_var_sigma(double sigma) {
    var_sigma_ = sigma;
  }
  void set_var_tau(int32_t tau) {
    var_tau_ = tau;
  }
  void set_var_sigma_tau(double sigma, int32_t tau) {
    var_sigma_  = sigma;
    var_tau_    = tau;
  }
  void set_noise_method(PhenotypicTargetNoiseMethod noise_method) {
    noise_method_ = noise_method;
  }
  void set_noise_prng(std::shared_ptr<JumpingMT> prng) {
    noise_prng_ = prng;
  }
  void set_noise_sigma(double sigma) {
    noise_sigma_ = sigma;
  }
  void set_noise_alpha(double alpha) {
    noise_alpha_ = alpha;
  }
  void set_noise_prob(double prob) {
    noise_prob_ = prob;
  }
  void set_noise_sampling_log(int8_t sampling_log) {
    noise_sampling_log_ = sampling_log;
  }

  /// PRNG used for variation
  std::shared_ptr<JumpingMT> var_prng_;
 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================
  void ApplyAutoregressiveMeanVariation();
  void ApplyAutoregressiveHeightVariation();

  // ==========================================================================
  //                               Attributes
  // ==========================================================================
  // ------------------------------------------------ Current Phenotypic Target
  std::unique_ptr<PhenotypicTarget> phenotypic_target_;

  // ---------------------------------------------------------------- Gaussians
  /// Phenotypic target's constitutive Gaussians in their initial state
  std::list<Gaussian> initial_gaussians_;
    /// Phenotypic target's constitutive Gaussians in their current state
  std::list<Gaussian> current_gaussians_;

  // ----------------------------------------------------------------- Sampling
  /// Number of points to be generated from the gaussians.
  int16_t sampling_;

  // ---------------------------------------------------------------- Variation
  /// Variation method
  PhenotypicTargetVariationMethod var_method_;
  /// Autoregressive mean variation sigma parameter
  double var_sigma_;
  /// Autoregressive mean variation tau parameter
  int16_t var_tau_;

  // -------------------------------------------------------------------- Noise
  /// Current noise (pure noise that is added to the phenotypic target)
  AbstractFuzzy* cur_noise_ = NULL;
  /// PRNG used for noise
  std::shared_ptr<JumpingMT> noise_prng_;
  PhenotypicTargetNoiseMethod noise_method_;
  /// Alpha value (variance coefficient)
  double noise_alpha_;
  /// Variance of the noise
  double noise_sigma_;
  /// Probability of variation.
  double noise_prob_;
  /// Log2 of the number of points in the noise fuzzy_set
  int8_t noise_sampling_log_;
};

// ============================================================================
//                       Inline functions' definition
// ============================================================================

} // namespace aevol

#endif // AEVOL_PHENOTYPIC_TARGET_HANDLER_H_
