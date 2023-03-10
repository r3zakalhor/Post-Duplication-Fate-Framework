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


#ifndef AEVOL_MUTATION_PARAMS_H_
#define AEVOL_MUTATION_PARAMS_H_


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <zlib.h>



// =================================================================
//                            Project Files
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "ae_enums.h"

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================






class MutationParams
{
  public :
    // =================================================================
    //                             Constructors
    // =================================================================
    MutationParams();
    MutationParams(const MutationParams & model);
    MutationParams(gzFile backup_file);


    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~MutationParams();

    // =================================================================
    //                        Accessors: getters
    // =================================================================

    // --------------------------------------------------------- Mutation rates
    inline double   point_mutation_rate() const;
    inline double   small_insertion_rate() const;
    inline double   small_deletion_rate() const;
    inline int16_t  max_indel_size() const;

    // -------------------------------------------- Rearrangements and Transfer
    inline bool   with_4pts_trans() const;
    inline bool   with_alignments() const;
    inline bool   with_HT() const;
    inline bool   repl_HT_with_close_points() const;
    inline double HT_ins_rate() const;
    inline double HT_repl_rate() const;
    inline double repl_HT_detach_rate() const;

    // ------------------------------ Rearrangement rates (without alignements)
    inline double duplication_rate() const;
    inline double deletion_rate() const;
    inline double translocation_rate() const;
    inline double inversion_rate() const;

    // --------------------------------- Rearrangement rates (with alignements)
    inline double neighbourhood_rate() const;
    inline double duplication_proportion() const;
    inline double deletion_proportion() const;
    inline double translocation_proportion() const;
    inline double inversion_proportion() const;

    // ------------------------------------------------------------ Alignements
    inline AlignmentFunctionShape align_fun_shape() const;
    inline double  align_sigm_lambda() const;
    inline int16_t align_sigm_mean() const;
    inline int16_t align_lin_min() const;
    inline int16_t align_lin_max() const;

    // Maximum shift of one seq on the other
    inline int16_t align_max_shift() const;
    // Work zone half length
    inline int16_t align_w_zone_h_len() const;
    // Corresponding residues match bonus
    inline int16_t align_match_bonus() const;
    // Corresponding residues mismatch cost
    inline int16_t align_mismatch_cost() const;

    // =================================================================
    //                        Accessors: setters
    // =================================================================

    // --------------------------------------------------------- Mutation rates
    inline void set_point_mutation_rate(double point_mutation_rate);
    inline void set_small_insertion_rate(double small_insertion_rate);
    inline void set_small_deletion_rate(double small_deletion_rate);
    inline void set_max_indel_size(int16_t max_indel_size);

    // -------------------------------------------- Rearrangements and Transfer
    inline void set_with_4pts_trans(bool with_4pts_trans);
    inline void set_with_alignments(bool with_alignments);
    inline void set_with_HT(bool with_HT);
    inline void set_repl_HT_with_close_points(bool repl_HT_with_close_points);
    inline void set_HT_ins_rate(double HT_ins_rate);
    inline void set_HT_repl_rate(double HT_repl_rate);
    inline void set_repl_HT_detach_rate(double repl_HT_detach_rate);

    // ------------------------------ Rearrangement rates (without alignements)
    inline void set_duplication_rate(double duplication_rate);
    inline void set_deletion_rate(double deletion_rate);
    inline void set_translocation_rate(double translocation_rate);
    inline void set_inversion_rate(double inversion_rate);

    // --------------------------------- Rearrangement rates (with alignements)
    inline void set_neighbourhood_rate(double neighbourhood_rate);
    inline void set_duplication_proportion(double duplication_proportion);
    inline void set_deletion_proportion(double deletion_proportion);
    inline void set_translocation_proportion(double translocation_proportion);
    inline void set_inversion_proportion(double inversion_proportion);

    // ------------------------------------------------------------ Alignements
    inline void set_align_fun_shape(AlignmentFunctionShape align_fun_shape);
    inline void set_align_sigm_lambda(double align_sigm_lambda);
    inline void set_align_sigm_mean(int16_t align_sigm_mean);
    inline void set_align_lin_min(int16_t align_lin_min);
    inline void set_align_lin_max(int16_t align_lin_max);

    inline void set_align_max_shift(int16_t align_max_shift);
    inline void set_align_w_zone_h_len(int16_t align_w_zone_h_len);
    inline void set_align_match_bonus(int16_t align_match_bonus);
    inline void set_align_mismatch_cost(int16_t align_mismatch_cost);

    // =================================================================
    //                              Operators
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    void save(gzFile backup_file) const;

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    /*MutationParams()
    {
      printf("%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    };
    MutationParams(const MutationParams &model)
    {
      printf("%s:%d: error: call to forbidden constructor.\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    };*/


    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================

    // --------------------------------------------------------- Mutation rates
    double  point_mutation_rate_;
    double  small_insertion_rate_;
    double  small_deletion_rate_;
    int16_t max_indel_size_;

    // -------------------------------------------- Rearrangements and Transfer
    bool    with_4pts_trans_;
    bool    with_alignments_;
    bool    with_HT_;
    bool    repl_HT_with_close_points_;
    double  HT_ins_rate_;
    double  HT_repl_rate_;
    double  repl_HT_detach_rate_;

    // ------------------------------ Rearrangement rates (without alignements)
    double duplication_rate_;
    double deletion_rate_;
    double translocation_rate_;
    double inversion_rate_;

    // --------------------------------- Rearrangement rates (with alignements)
    double neighbourhood_rate_;
    double duplication_proportion_;
    double deletion_proportion_;
    double translocation_proportion_;
    double inversion_proportion_;

    // ------------------------------------------------------------ Alignements
    AlignmentFunctionShape align_fun_shape_;
    double  align_sigm_lambda_;
    int16_t align_sigm_mean_;
    int16_t align_lin_min_;
    int16_t align_lin_max_;

    int16_t align_max_shift_;     // Maximum shift of one seq on the other
    int16_t align_w_zone_h_len_;  // Work zone half length
    int16_t align_match_bonus_;   // Corresponding residues match bonus
    int16_t align_mismatch_cost_; // Corresponding residues mismatch cost
};


// =====================================================================
//                           Getters' definitions
// =====================================================================

// ------------------------------------------------------------- Mutation rates
inline double MutationParams::point_mutation_rate() const
{
  return point_mutation_rate_;
}

inline double MutationParams::small_insertion_rate() const
{
  return small_insertion_rate_;
}

inline double MutationParams::small_deletion_rate() const
{
  return small_deletion_rate_;
}

inline int16_t MutationParams::max_indel_size() const
{
  return max_indel_size_;
}

// ------------------------------------------------ Rearrangements and Transfer
inline bool MutationParams::with_4pts_trans() const
{
  return with_4pts_trans_;
}

inline bool MutationParams::with_alignments() const
{
  return with_alignments_;
}

inline bool MutationParams::with_HT() const
{
  return with_HT_;
}

inline bool MutationParams::repl_HT_with_close_points() const
{
  return repl_HT_with_close_points_;
}

inline double MutationParams::HT_ins_rate() const
{
  return HT_ins_rate_;
}

inline double MutationParams::HT_repl_rate() const
{
  return HT_repl_rate_;
}

inline double MutationParams::repl_HT_detach_rate() const
{
  return repl_HT_detach_rate_;
}

// ---------------------------------- Rearrangement rates (without alignements)
inline double MutationParams::duplication_rate() const
{
  return duplication_rate_;
}

inline double MutationParams::deletion_rate() const
{
  return deletion_rate_;
}

inline double MutationParams::translocation_rate() const
{
  return translocation_rate_;
}

inline double MutationParams::inversion_rate() const
{
  return inversion_rate_;
}

// ------------------------------------- Rearrangement rates (with alignements)
inline double MutationParams::neighbourhood_rate() const
{
  return neighbourhood_rate_;
}

inline double MutationParams::duplication_proportion() const
{
  return duplication_proportion_;
}

inline double MutationParams::deletion_proportion() const
{
  return deletion_proportion_;
}

inline double MutationParams::translocation_proportion() const
{
  return translocation_proportion_;
}

inline double MutationParams::inversion_proportion() const
{
  return inversion_proportion_;
}

// ---------------------------------------------------------------- Alignements
inline AlignmentFunctionShape MutationParams::align_fun_shape() const
{
 return align_fun_shape_;
}

inline double MutationParams::align_sigm_lambda() const
{
 return align_sigm_lambda_;
}

inline int16_t MutationParams::align_sigm_mean() const
{
 return align_sigm_mean_;
}

inline int16_t MutationParams::align_lin_min() const
{
 return align_lin_min_;
}

inline int16_t MutationParams::align_lin_max() const
{
 return align_lin_max_;
}

inline int16_t MutationParams::align_max_shift() const
{
 return align_max_shift_;
}

inline int16_t MutationParams::align_w_zone_h_len() const
{
 return align_w_zone_h_len_;
}

inline int16_t MutationParams::align_match_bonus() const
{
 return align_match_bonus_;
}

inline int16_t MutationParams::align_mismatch_cost() const
{
 return align_mismatch_cost_;
}

// =====================================================================
//                           Setters' definitions
// =====================================================================

// ------------------------------------------------------------- Mutation rates
inline void MutationParams::set_point_mutation_rate(double point_mutation_rate)
{
  point_mutation_rate_ = point_mutation_rate;
}

inline void MutationParams::set_small_insertion_rate(double small_insertion_rate)
{
  small_insertion_rate_ = small_insertion_rate;
}

inline void MutationParams::set_small_deletion_rate(double small_deletion_rate)
{
  small_deletion_rate_ = small_deletion_rate;
}

inline void MutationParams::set_max_indel_size(int16_t max_indel_size)
{
  max_indel_size_ = max_indel_size;
}

// ------------------------------------------------ Rearrangements and Transfer
inline void MutationParams::set_with_4pts_trans(bool with_4pts_trans)
{
  with_4pts_trans_ = with_4pts_trans;
}

inline void MutationParams::set_with_alignments(bool with_alignments)
{
  with_alignments_ = with_alignments;
}

inline void MutationParams::set_with_HT(bool with_HT)
{
  with_HT_ = with_HT;
}

inline void MutationParams::set_repl_HT_with_close_points(bool repl_HT_with_close_points)
{
  repl_HT_with_close_points_ = repl_HT_with_close_points;
}

inline void MutationParams::set_HT_ins_rate(double HT_ins_rate)
{
  HT_ins_rate_ = HT_ins_rate;
}

inline void MutationParams::set_HT_repl_rate(double HT_repl_rate)
{
  HT_repl_rate_ = HT_repl_rate;
}

inline void MutationParams::set_repl_HT_detach_rate(double repl_HT_detach_rate)
{
  repl_HT_detach_rate_ = repl_HT_detach_rate;
}

// ---------------------------------- Rearrangement rates (without alignements)
inline void MutationParams::set_duplication_rate(double duplication_rate)
{
  duplication_rate_ = duplication_rate;
}

inline void MutationParams::set_deletion_rate(double deletion_rate)
{
  deletion_rate_ = deletion_rate;
}

inline void MutationParams::set_translocation_rate(double translocation_rate)
{
  translocation_rate_ = translocation_rate;
}

inline void MutationParams::set_inversion_rate(double inversion_rate)
{
  inversion_rate_ = inversion_rate;
}

// ------------------------------------- Rearrangement rates (with alignements)
inline void MutationParams::set_neighbourhood_rate(double neighbourhood_rate)
{
  neighbourhood_rate_ = neighbourhood_rate;
}

inline void MutationParams::set_duplication_proportion(double duplication_proportion)
{
  duplication_proportion_ = duplication_proportion;
}

inline void MutationParams::set_deletion_proportion(double deletion_proportion)
{
  deletion_proportion_ = deletion_proportion;
}

inline void MutationParams::set_translocation_proportion(double translocation_proportion)
{
  translocation_proportion_ = translocation_proportion;
}

inline void MutationParams::set_inversion_proportion(double inversion_proportion)
{
  inversion_proportion_ = inversion_proportion;
}

// ---------------------------------------------------------------- Alignements
inline void MutationParams::set_align_fun_shape(
    AlignmentFunctionShape align_fun_shape)
{
  align_fun_shape_ = align_fun_shape;
}

inline void MutationParams::set_align_sigm_lambda(double align_sigm_lambda)
{
  align_sigm_lambda_ = align_sigm_lambda;
}

inline void MutationParams::set_align_sigm_mean(int16_t align_sigm_mean)
{
  align_sigm_mean_ = align_sigm_mean;
}

inline void MutationParams::set_align_lin_min(int16_t align_lin_min)
{
  align_lin_min_ = align_lin_min;
}

inline void MutationParams::set_align_lin_max(int16_t align_lin_max)
{
  align_lin_max_ = align_lin_max;
}

inline void MutationParams::set_align_max_shift(int16_t align_max_shift)
{
  align_max_shift_ = align_max_shift;
}

inline void MutationParams::set_align_w_zone_h_len(int16_t align_w_zone_h_len)
{
  align_w_zone_h_len_ = align_w_zone_h_len;
}

inline void MutationParams::set_align_match_bonus(int16_t align_match_bonus)
{
  align_match_bonus_ = align_match_bonus;
}

inline void MutationParams::set_align_mismatch_cost(int16_t align_mismatch_cost)
{
  align_mismatch_cost_ = align_mismatch_cost;
}

// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol

#endif // AEVOL_MUTATION_PARAMS_H_
