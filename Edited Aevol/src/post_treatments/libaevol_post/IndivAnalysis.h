//
// Created by dparsons on 31/05/16.
//

#ifndef AEVOL_INDIVANALYSIS_H__
#define AEVOL_INDIVANALYSIS_H__


// ============================================================================
//                                   Includes
// ============================================================================
#include "Individual.h"

namespace aevol {

/**
 *
 */
class IndivAnalysis : public Individual {
 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  // base constructor deleted in mother class
  IndivAnalysis(const IndivAnalysis&) = delete; //< Copy ctor
  IndivAnalysis(IndivAnalysis&&) = delete; //< Move ctor

  explicit IndivAnalysis(const Individual&);

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  ~IndivAnalysis() override = default; //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================
  /// Copy assignment
  IndivAnalysis& operator=(const IndivAnalysis& other) = delete;
  /// Move assignment
  IndivAnalysis& operator=(IndivAnalysis&& other) = delete;

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  double compute_theoritical_f_nu();
  void compute_experimental_f_nu(int32_t nb_indiv,
                                 std::shared_ptr<JumpingMT> prng,
                                 FILE* output_summary = nullptr,
                                 bool verbose = false,
				 bool full_output = false);

  void compute_experimental_mutagenesis(int32_t nb_indiv,
					int32_t mutation_type,
					std::shared_ptr<JumpingMT> prng,
					FILE* output_summary = nullptr,
					bool verbose = false,
					bool full_output = false);
  

  // ==========================================================================
  //                                Accessors
  // ==========================================================================

 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================

  // ==========================================================================
  //                               Attributes
  // ==========================================================================
};

} // namespace aevol
#endif //AEVOL_INDIVANALYSIS_H__
