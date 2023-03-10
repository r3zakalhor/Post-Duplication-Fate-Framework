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




// ============================================================================
//                                   Includes
// ============================================================================
#if __cplusplus == 201103L
#include "make_unique.h"
#endif

#include "PhenotypicTargetHandler.h"
#include "ExpSetup.h"
#include "HybridFuzzy.h"
#include "Utils.h"

#include <iostream>


using std::cout;
using std::endl;


namespace aevol {


//##############################################################################
//                                                                             #
//                        Class PhenotypicTargetHandler                        #
//                                                                             #
//##############################################################################

// ============================================================================
//                       Definition of static attributes
// ============================================================================

// ============================================================================
//                                Constructors
// ============================================================================
PhenotypicTargetHandler::PhenotypicTargetHandler() {
  // The phenotypic target
#if __cplusplus == 201103L
  phenotypic_target_ = make_unique<PhenotypicTarget>();
#else
  phenotypic_target_ = std::make_unique<PhenotypicTarget>();
#endif

  // Sampling
  sampling_ = 0;

  // Variation
  var_prng_   = NULL;
  var_method_ = NO_VAR;
  var_sigma_  = 0.0;
  var_tau_    = 0;

  // Noise
  cur_noise_          = NULL;
  noise_method_       = NO_NOISE;
  noise_prng_         = NULL;
  noise_prob_         = 0.0;
  noise_alpha_        = 0.0;
  noise_sigma_        = 0.0;
  noise_sampling_log_ = 8;
}

PhenotypicTargetHandler::PhenotypicTargetHandler(
    const PhenotypicTargetHandler& rhs) {
  // ------------------------------------------------ Current Phenotypic Target
#if __cplusplus == 201103L
  phenotypic_target_ = make_unique<PhenotypicTarget>(*(rhs.phenotypic_target_));
#else
  phenotypic_target_ = std::make_unique<PhenotypicTarget>(*(rhs.phenotypic_target_));
#endif

  // ---------------------------------------------------------------- Gaussians
  initial_gaussians_ = rhs.initial_gaussians_;
  current_gaussians_ = rhs.current_gaussians_;

  // ----------------------------------------------------------------- Sampling
  sampling_ = rhs.sampling_;

  // ---------------------------------------------------------------- Variation
  var_method_ = rhs.var_method_;
  var_prng_ = rhs.var_prng_;
  var_sigma_ = rhs.var_sigma_;
  var_tau_ = rhs.var_tau_;

  // -------------------------------------------------------------------- Noise
  cur_noise_ = rhs.cur_noise_;
  noise_prng_ = rhs.noise_prng_;
  noise_method_ = rhs.noise_method_;
  noise_alpha_ = rhs.noise_alpha_;
  noise_sigma_ = rhs.noise_sigma_;
  noise_prob_ = rhs.noise_prob_;
  noise_sampling_log_ = rhs.noise_sampling_log_;
}

PhenotypicTargetHandler::PhenotypicTargetHandler(gzFile backup_file) {
  load(backup_file);
}

// ============================================================================
//                                 Destructor
// ============================================================================
PhenotypicTargetHandler::~PhenotypicTargetHandler() {
  delete cur_noise_;
}

// ============================================================================
//                                   Methods
// ============================================================================
void PhenotypicTargetHandler::BuildPhenotypicTarget() {
  // NB : Extreme points (at abscissa X_MIN and X_MAX) will be generated, we need to erase the list first
  phenotypic_target_->fuzzy()->clear();

  // Generate sample points from gaussians
  if (not current_gaussians_.empty()) {

      for (int16_t i = 0; i <= sampling_; i++) {
          Point new_point = Point(
                  X_MIN + (double) i * (X_MAX - X_MIN) / (double) sampling_, 0.0);
          for (const Gaussian &g: current_gaussians_)
              new_point.y += g.compute_y(new_point.x);


//          printf("Add point [%f %d] = %f\n",new_point.x,i,new_point.y);

//          if (FuzzyFactory::fuzzyFactory->get_fuzzy_flavor() == 0) {
              phenotypic_target_->fuzzy()->add_point(new_point.x, new_point.y);
//              printf("Point added [%d %f] = %f\n",i,new_point.x,((Fuzzy*)phenotypic_target_->fuzzy())->y(new_point.x));
//          } else {
//              ((HybridFuzzy *) phenotypic_target_->fuzzy())->points()[i] = new_point.y;
////              printf("Point added [%d %f] = %f\n",i,new_point.x,((HybridFuzzy*)phenotypic_target_->fuzzy())->points()[i]);
//          }
      }

  }



  // Add lower and upper bounds
  phenotypic_target_->fuzzy()->clip(AbstractFuzzy::min, Y_MIN);
  phenotypic_target_->fuzzy()->clip(AbstractFuzzy::max, Y_MAX);

//        if (FuzzyFactory::fuzzyFactory->get_fuzzy_flavor() == 0) {
//            printf("Size of points %d\n",((Fuzzy*)phenotypic_target_->fuzzy())->points().size());
//        }
//
//  printf("=========== BEFORE SIMPLIFY============\n");
  phenotypic_target_->ComputeArea();
  // Simplify (get rid of useless points)
  phenotypic_target_->fuzzy()->simplify();
//        printf("=========== AFTER SIMPLIFY ============\n");
  // Compute areas (total and by feature)
  phenotypic_target_->ComputeArea();
}

void PhenotypicTargetHandler::ApplyVariation() {
  switch (var_method_) {
    case NO_VAR :
      return;
    case AUTOREGRESSIVE_MEAN_VAR :
      ApplyAutoregressiveMeanVariation();
      break;
    case AUTOREGRESSIVE_HEIGHT_VAR :
      ApplyAutoregressiveHeightVariation();
      break;
    default :
      Utils::ExitWithDevMsg("Unknown variation method", __FILE__, __LINE__);
  }

  // Phenotypic target has changed, recompute its area
  phenotypic_target_->ComputeArea();
}

void PhenotypicTargetHandler::ApplyAutoregressiveMeanVariation() {
  // For each gaussian :
  // current_mean = ref_mean + delta_m, where
  // delta_m follows an autoregressive stochastic process
  // with the parameters var_sigma_ and var_tau_
  for (auto cur_gaussian = current_gaussians_.begin(),
           initial_gaussian = initial_gaussians_.begin() ;
       cur_gaussian != current_gaussians_.end() ;
       cur_gaussian++, initial_gaussian++) {
    // Find the current delta_mean = current_mean - ref_mean
    double delta_mean = cur_gaussian->mean() - initial_gaussian->mean();

    // Compute the next value :
    // Dm(t+1) = Dm(t)*(1-1/tau) + ssd/tau*sqrt(2*tau-1)*normal_random()
    Utils::ApplyAutoregressiveStochasticProcess(delta_mean,
                                                var_sigma_,
                                                var_tau_,
                                                *var_prng_);

    // Deduce the new value of the mean : ref_mean + delta_m
    cur_gaussian->set_mean(initial_gaussian->mean() + delta_mean);
  }

  BuildPhenotypicTarget();
}

void PhenotypicTargetHandler::ApplyAutoregressiveHeightVariation() {
  // For each gaussian :
  // current_height = ref_height + delta_h, where
  // delta_h follows an autoregressive stochastic process
  // with the parameters var_sigma_ and var_tau_
  for (auto cur_gaussian = current_gaussians_.begin(),
           initial_gaussian = initial_gaussians_.begin() ;
       cur_gaussian != current_gaussians_.end() ;
       cur_gaussian++, initial_gaussian++) {
    // Find the current delta_height = current_height - ref_height
    double delta_height = cur_gaussian->height() -
      initial_gaussian->height();

    // Compute the next value :
    // Dm(t+1) = Dm(t)*(1-1/tau) + ssd/tau*sqrt(2*tau-1)*normal_random()
    Utils::ApplyAutoregressiveStochasticProcess(delta_height,
                                                var_sigma_,
                                                var_tau_,
                                                *var_prng_);

    // Deduce the new value of the height : ref_height + delta_h
    cur_gaussian->set_height(initial_gaussian->height() + delta_height);
  }

  BuildPhenotypicTarget();
}

void PhenotypicTargetHandler::save(gzFile backup_file) const {
  //printf("Appel a la sauvegarde de PhenotypicTargetHandler\n");
  // --------------------------------------------------------------------------
  //  Write phenotypic target segmentation
  phenotypic_target_->SaveSegmentation(backup_file);

  // --------------------------------------------------------------------------
  //  Write current gaussians (initial gaussians will be stored later if
  // necessary)
  int8_t nb_gaussians = current_gaussians_.size();
  gzwrite(backup_file, &nb_gaussians, sizeof(nb_gaussians));
  for (const Gaussian & g: current_gaussians_)
    g.save(backup_file);

  // --------------------------------------------------------------------------
  //  Write sampling
  gzwrite(backup_file, &sampling_, sizeof(sampling_));

  // --------------------------------------------------------------------------
  //  Write variation data
  int8_t tmp_var_method = var_method_;
  gzwrite(backup_file, &tmp_var_method,  sizeof(tmp_var_method));

  if (var_method_ != NO_VAR) {
    var_prng_->save(backup_file);
    gzwrite(backup_file, &var_sigma_, sizeof(var_sigma_));
    gzwrite(backup_file, &var_tau_,   sizeof(var_tau_));
  }

  // --------------------------------------------------------------------------
  //  Write noise data
  int8_t tmp_noise_method = noise_method_;
  gzwrite(backup_file, &tmp_noise_method, sizeof(tmp_noise_method));

  if (noise_method_ != NO_NOISE) {
    int8_t tmp_save_cur_noise = (cur_noise_ != NULL);
    gzwrite(backup_file, &tmp_save_cur_noise,  sizeof(tmp_save_cur_noise));
    if (tmp_save_cur_noise) cur_noise_->save(backup_file);

    noise_prng_->save(backup_file);
    gzwrite(backup_file, &noise_alpha_,  sizeof(noise_alpha_));
    gzwrite(backup_file, &noise_sigma_,  sizeof(noise_sigma_));
    gzwrite(backup_file, &noise_prob_,   sizeof(noise_prob_));
    gzwrite(backup_file, &noise_sampling_log_, sizeof(noise_sampling_log_));
  }

  // ---------------------------------------------------------------
  //  If needed, keep a copy of the initial state of the gaussians
  // ---------------------------------------------------------------
  if (var_method_ != NO_VAR || noise_method_ != NO_NOISE) {
    size_t nb_gaussians = initial_gaussians_.size();
    gzwrite(backup_file, &nb_gaussians, sizeof(nb_gaussians));
    for (const Gaussian & g: initial_gaussians_)
      g.save(backup_file);
  }
}

void PhenotypicTargetHandler::load(gzFile backup_file) {
  //printf("Appel au chargement de PhenotypicTargetHandler\n");
  // --------------------------------------------------------------------------
  //  Retrieve phenotypic target segmentation
#if __cplusplus == 201103L
  phenotypic_target_ = make_unique<PhenotypicTarget>();
#else
  phenotypic_target_ = std::make_unique<PhenotypicTarget>();
#endif

  phenotypic_target_->LoadSegmentation(backup_file);

  // --------------------------------------------------------------------------
  //  Retrieve current gaussians
  int8_t nb_gaussians;
  gzread(backup_file, &nb_gaussians, sizeof(nb_gaussians));
  for (int8_t i = 0 ; i < nb_gaussians ; i++)
    current_gaussians_.push_back(Gaussian(backup_file));

  // --------------------------------------------------------------------------
  //  Retrieve sampling
  gzread(backup_file, &sampling_, sizeof(sampling_));

  // --------------------------------------------------------------------------
  //  Retrieve variation data
  int8_t tmp_var_method;
  gzread(backup_file, &tmp_var_method, sizeof(tmp_var_method));
  var_method_ = (PhenotypicTargetVariationMethod) tmp_var_method;

  if (var_method_ != NO_VAR) {
    var_prng_ = std::make_shared<JumpingMT>(backup_file);
    gzread(backup_file, &var_sigma_, sizeof(var_sigma_));
    gzread(backup_file, &var_tau_, sizeof(var_tau_));
  }

  // --------------------------------------------------------------------------
  //  Retrieve noise data
  int8_t tmp_noise_method;
  gzread(backup_file, &tmp_noise_method, sizeof(tmp_noise_method));
  noise_method_ = (PhenotypicTargetNoiseMethod) tmp_noise_method;

  if (noise_method_ != NO_NOISE) {
    int8_t tmp_cur_noise_saved;
    gzread(backup_file, &tmp_cur_noise_saved,  sizeof(tmp_cur_noise_saved));
    if (tmp_cur_noise_saved)
      cur_noise_ = FuzzyFactory::fuzzyFactory->create_fuzzy(backup_file);

    noise_prng_ = std::make_shared<JumpingMT>(backup_file);
    gzread(backup_file, &noise_alpha_, sizeof(noise_alpha_));
    gzread(backup_file, &noise_sigma_, sizeof(noise_sigma_));
    gzread(backup_file, &noise_prob_,  sizeof(noise_prob_));
    gzread(backup_file, &noise_sampling_log_, sizeof(noise_sampling_log_));
  }

  // --------------------------------------------------------------------------
  //  If needed, retrieve a copy of the initial state of the gaussians
  if (var_method_ != NO_VAR || noise_method_ != NO_NOISE) {
    size_t nb_gaussians;
    gzread(backup_file, &nb_gaussians, sizeof(nb_gaussians));
    
    for (size_t i = 0 ; i < nb_gaussians ; i++)
      initial_gaussians_.emplace_back(backup_file);
  }

  // --------------------------------------------------------------------------
  //  Build the phenotypic target
  BuildPhenotypicTarget();
}

// ============================================================================
//                            Non inline accessors
// ============================================================================

} // namespace aevol
