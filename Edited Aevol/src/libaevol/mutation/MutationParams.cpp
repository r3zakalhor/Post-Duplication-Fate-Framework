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




// =================================================================
//                              Libraries
// =================================================================



// =================================================================
//                            Project Files
// =================================================================
#include "MutationParams.h"
#include "Alignment.h"

namespace aevol {

//##############################################################################
//                                                                             #
//                             Class MutationParams                             #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
MutationParams::MutationParams()
{
  // --------------------------------------------------------- Mutation rates
  point_mutation_rate_  = 0.0;
  small_insertion_rate_ = 0.0;
  small_deletion_rate_  = 0.0;
  max_indel_size_       = 0;

  // -------------------------------------------- Rearrangements and Transfer
  with_4pts_trans_            = false;
  with_alignments_            = false;
  with_HT_                    = false;
  repl_HT_with_close_points_  = false;
  HT_ins_rate_                = 0.0;
  HT_repl_rate_               = 0.0;
  repl_HT_detach_rate_        = 0.0;

  // ------------------------------ Rearrangement rates (without alignements)
  duplication_rate_   = 0.0;
  deletion_rate_      = 0.0;
  translocation_rate_ = 0.0;
  inversion_rate_     = 0.0;

  // --------------------------------- Rearrangement rates (with alignements)
  neighbourhood_rate_ = 0.0;

  duplication_proportion_   = 0.0;
  deletion_proportion_      = 0.0;
  translocation_proportion_ = 0.0;
  inversion_proportion_     = 0.0;

  // ------------------------------------------------------------ Alignements
  align_fun_shape_    = SIGMOID;
  align_sigm_lambda_  = 4;
  align_sigm_mean_    = 50;
  align_lin_min_      = 0;
  align_lin_max_      = 100;

  align_max_shift_      = 20;
  align_w_zone_h_len_   = 50;
  align_match_bonus_    = 1;
  align_mismatch_cost_  = 2;
}

MutationParams::MutationParams(const MutationParams & model)
{
  // --------------------------------------------------------- Mutation rates
  point_mutation_rate_  = model.point_mutation_rate_;
  small_insertion_rate_ = model.small_insertion_rate_;
  small_deletion_rate_  = model.small_deletion_rate_;
  max_indel_size_       = model.max_indel_size_;

  // -------------------------------------------- Rearrangements and Transfer
  with_4pts_trans_            = model.with_4pts_trans_;
  with_alignments_            = model.with_alignments_;
  with_HT_                    = model.with_HT_;
  repl_HT_with_close_points_  = model.repl_HT_with_close_points_;
  HT_ins_rate_                = model.HT_ins_rate_;
  HT_repl_rate_               = model.HT_repl_rate_;
  repl_HT_detach_rate_        = model.repl_HT_detach_rate_;

  // ------------------------------ Rearrangement rates (without alignements)
  duplication_rate_   = model.duplication_rate_;
  deletion_rate_      = model.deletion_rate_;
  translocation_rate_ = model.translocation_rate_;
  inversion_rate_     = model.inversion_rate_;

  // --------------------------------- Rearrangement rates (with alignements)
  neighbourhood_rate_ = model.neighbourhood_rate_;

  duplication_proportion_   = model.duplication_proportion_;
  deletion_proportion_      = model.deletion_proportion_;
  translocation_proportion_ = model.translocation_proportion_;
  inversion_proportion_     = model.inversion_proportion_;

  // ------------------------------------------------------------ Alignements
  align_fun_shape_    = model.align_fun_shape_;
  align_sigm_lambda_  = model.align_sigm_lambda_;
  align_sigm_mean_    = model.align_sigm_mean_;
  align_lin_min_      = model.align_lin_min_;
  align_lin_max_      = model.align_lin_max_;

  align_max_shift_      = model.align_max_shift_;
  align_w_zone_h_len_   = model.align_w_zone_h_len_;
  align_match_bonus_    = model.align_match_bonus_;
  align_mismatch_cost_  = model.align_mismatch_cost_;
}

MutationParams::MutationParams(gzFile backup_file)
{
  // --------------------------------------------------------- Mutation rates
  gzread(backup_file, &point_mutation_rate_,  sizeof(point_mutation_rate_));
  gzread(backup_file, &small_insertion_rate_, sizeof(small_insertion_rate_));
  gzread(backup_file, &small_deletion_rate_,  sizeof(small_deletion_rate_));
  gzread(backup_file, &max_indel_size_,       sizeof(max_indel_size_));

  // -------------------------------------------- Rearrangements and Transfer
  int8_t tmp;
  gzread(backup_file, &tmp, sizeof(tmp));
  with_4pts_trans_ = (tmp != 0);
  gzread(backup_file, &tmp, sizeof(tmp));
  with_alignments_ = (tmp != 0);
  gzread(backup_file, &tmp, sizeof(tmp));
  with_HT_ = (tmp != 0);
  gzread(backup_file, &tmp, sizeof(tmp));
  repl_HT_with_close_points_ = (tmp != 0);
  gzread(backup_file, &HT_ins_rate_,  sizeof(HT_ins_rate_));
  gzread(backup_file, &HT_repl_rate_, sizeof(HT_repl_rate_));
  gzread(backup_file, &repl_HT_detach_rate_, sizeof(repl_HT_detach_rate_));

  // ------------------------------ Rearrangement rates (without alignements)
  gzread(backup_file, &duplication_rate_,   sizeof(duplication_rate_));
  gzread(backup_file, &deletion_rate_,      sizeof(deletion_rate_));
  gzread(backup_file, &translocation_rate_, sizeof(translocation_rate_));
  gzread(backup_file, &inversion_rate_,     sizeof(inversion_rate_));

  // --------------------------------- Rearrangement rates (with alignements)
  gzread(backup_file, &neighbourhood_rate_,       sizeof(neighbourhood_rate_));
  gzread(backup_file, &duplication_proportion_,   sizeof(duplication_proportion_));
  gzread(backup_file, &deletion_proportion_,      sizeof(deletion_proportion_));
  gzread(backup_file, &translocation_proportion_, sizeof(translocation_proportion_));
  gzread(backup_file, &inversion_proportion_,     sizeof(inversion_proportion_));

  // ------------------------------------------------------------ Alignements
  gzread(backup_file, &align_fun_shape_,     sizeof(align_fun_shape_));
  gzread(backup_file, &align_sigm_lambda_,   sizeof(align_sigm_lambda_));
  gzread(backup_file, &align_sigm_mean_,     sizeof(align_sigm_mean_));
  gzread(backup_file, &align_lin_min_,       sizeof(align_lin_min_));
  gzread(backup_file, &align_lin_max_,       sizeof(align_lin_max_));

  gzread(backup_file, &align_max_shift_,     sizeof(align_max_shift_));
  gzread(backup_file, &align_w_zone_h_len_,  sizeof(align_w_zone_h_len_));
  gzread(backup_file, &align_match_bonus_,   sizeof(align_match_bonus_));
  gzread(backup_file, &align_mismatch_cost_, sizeof(align_mismatch_cost_));

  //Alignment::align_fun_shape     = align_fun_shape_;
  //Alignment::align_sigm_lambda   = align_sigm_lambda_;
  //Alignment::align_sigm_mean     = align_sigm_mean_;
  //Alignment::align_lin_min       = align_lin_min_;
  //Alignment::align_lin_max       = align_lin_max_;

  //Alignment::align_max_shift     = align_max_shift_;
  //Alignment::align_w_zone_h_len  = align_w_zone_h_len_;
  //Alignment::align_match_bonus   = align_match_bonus_;
  //Alignment::align_mismatch_cost = align_mismatch_cost_;
}

// =================================================================
//                             Destructors
// =================================================================
MutationParams::~MutationParams()
{
}

// =================================================================
//                            Public Methods
// =================================================================
void MutationParams::save(gzFile backup_file) const
{
  // --------------------------------------------------------- Mutation rates
  gzwrite(backup_file, &point_mutation_rate_,  sizeof(point_mutation_rate_));
  gzwrite(backup_file, &small_insertion_rate_, sizeof(small_insertion_rate_));
  gzwrite(backup_file, &small_deletion_rate_,  sizeof(small_deletion_rate_));
  gzwrite(backup_file, &max_indel_size_,       sizeof(max_indel_size_));

  // -------------------------------------------- Rearrangements and Transfer
  int8_t tmp = with_4pts_trans_ ? 1 : 0;
  gzwrite(backup_file, &tmp,  sizeof(tmp));
  tmp = with_alignments_ ? 1 : 0;
  gzwrite(backup_file, &tmp,  sizeof(tmp));
  tmp = with_HT_ ? 1 : 0;
  gzwrite(backup_file, &tmp,  sizeof(tmp));
  tmp = repl_HT_with_close_points_ ? 1 : 0;
  gzwrite(backup_file, &tmp,  sizeof(tmp));
  gzwrite(backup_file, &HT_ins_rate_,  sizeof(HT_ins_rate_));
  gzwrite(backup_file, &HT_repl_rate_, sizeof(HT_repl_rate_));
  gzwrite(backup_file, &repl_HT_detach_rate_, sizeof(repl_HT_detach_rate_));

  // ------------------------------ Rearrangement rates (without alignements)
  gzwrite(backup_file, &duplication_rate_,   sizeof(duplication_rate_));
  gzwrite(backup_file, &deletion_rate_,      sizeof(deletion_rate_));
  gzwrite(backup_file, &translocation_rate_, sizeof(translocation_rate_));
  gzwrite(backup_file, &inversion_rate_,     sizeof(inversion_rate_));

  // --------------------------------- Rearrangement rates (with alignements)
  gzwrite(backup_file, &neighbourhood_rate_,       sizeof(neighbourhood_rate_));
  gzwrite(backup_file, &duplication_proportion_,   sizeof(duplication_proportion_));
  gzwrite(backup_file, &deletion_proportion_,      sizeof(deletion_proportion_));
  gzwrite(backup_file, &translocation_proportion_, sizeof(translocation_proportion_));
  gzwrite(backup_file, &inversion_proportion_,     sizeof(inversion_proportion_));

  // ------------------------------------------------------------ Alignements
  gzwrite(backup_file, &align_fun_shape_,     sizeof(align_fun_shape_));
  gzwrite(backup_file, &align_sigm_lambda_,   sizeof(align_sigm_lambda_));
  gzwrite(backup_file, &align_sigm_mean_,     sizeof(align_sigm_mean_));
  gzwrite(backup_file, &align_lin_min_,       sizeof(align_lin_min_));
  gzwrite(backup_file, &align_lin_max_,       sizeof(align_lin_max_));

  gzwrite(backup_file, &align_max_shift_,     sizeof(align_max_shift_));
  gzwrite(backup_file, &align_w_zone_h_len_,  sizeof(align_w_zone_h_len_));
  gzwrite(backup_file, &align_match_bonus_,   sizeof(align_match_bonus_));
  gzwrite(backup_file, &align_mismatch_cost_, sizeof(align_mismatch_cost_));
}

// =================================================================
//                           Protected Methods
// =================================================================

// =================================================================
//                          Non inline accessors
// =================================================================
} // namespace aevol
