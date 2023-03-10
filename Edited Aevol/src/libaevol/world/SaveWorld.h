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


#ifndef AEVOL_SAVEWORLD_H_
#define AEVOL_SAVEWORLD_H_


// ============================================================================
//                                   Includes
// ============================================================================
#include "GridCell.h"
#include "Individual_7.h"

#include <cassert>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <zlib.h>

namespace aevol {



// =================================================================
//                          Class declarations
// =================================================================

class SaveGridCell {
 public :

  SaveGridCell(ExpManager* exp_m, int16_t x, int16_t y,
#ifdef __REGUL
               std::unique_ptr<Habitat_R>&& habitat,
#else
              std::unique_ptr<Habitat>&& habitat,
#endif
              Individual * indiv, std::shared_ptr<JumpingMT> mut_prng,
              std::shared_ptr<JumpingMT> stoch_prng,
              std::shared_ptr<JumpingMT> reprod_prng, std::shared_ptr<JumpingMT> reprod_prng_simd, double w_max)
  {
    x_ = x;
    y_ = y;

    mut_prng_ = mut_prng;
    stoch_prng_ = stoch_prng;
    reprod_prng_ = reprod_prng;
    reprod_prng_simd_ = reprod_prng_simd;

    individual_ = indiv;
    habitat_ = std::move(habitat);

    exp_m_ = exp_m;
    w_max_ = w_max;
  };

    SaveGridCell(ExpManager* exp_m, int16_t x, int16_t y,
#ifdef __REGUL
                 std::unique_ptr<Habitat_R>&& habitat,
#else
            std::unique_ptr<Habitat>&& habitat,
#endif
               Individual_7* indiv, std::shared_ptr<JumpingMT> mut_prng,
                 std::shared_ptr<JumpingMT> stoch_prng,
                 std::shared_ptr<JumpingMT> reprod_prng, std::shared_ptr<JumpingMT> reprod_prng_simd, double w_max)
    {
      x_ = x;
      y_ = y;

      mut_prng_ = mut_prng;
      stoch_prng_ = stoch_prng;
      reprod_prng_ = reprod_prng;
      reprod_prng_simd_ = reprod_prng_simd;

      simd_individual_ = indiv;
      habitat_ = std::move(habitat);

      exp_m_ = exp_m;
      w_max_ = w_max;
    };


    virtual ~SaveGridCell() = default;

  void save(gzFile backup_file,
                          bool skip_phenotypic_target /*=false*/, bool create = false);

  int16_t x_;
  int16_t y_;

  double w_max_;

  // Pointer to the individual in this cell
  Individual* individual_ = NULL;
  Individual_7* simd_individual_ = NULL;
#ifdef __REGUL
        std::unique_ptr<Habitat_R> habitat_ = nullptr;
#else
  std::unique_ptr<Habitat> habitat_ = nullptr;
#endif
  std::shared_ptr<JumpingMT> mut_prng_ = nullptr;
  std::shared_ptr<JumpingMT> stoch_prng_ = nullptr;
  std::shared_ptr<JumpingMT> reprod_prng_ = nullptr;
  std::shared_ptr<JumpingMT> reprod_prng_simd_ = nullptr;

  ExpManager* exp_m_;
};

class SaveWorld
{

 friend class World;
 public :
  // =================================================================
  //                             Constructors
  // =================================================================
  SaveWorld() = default;
  SaveWorld(const SaveWorld &) = delete;

  // =================================================================
  //                             Destructors
  // =================================================================
  virtual ~SaveWorld();

  // =================================================================
  //                            Public Methods
  // =================================================================
  void save(gzFile backup_file, bool create = false);

 protected :
  // =================================================================
  //                           Protected Methods
  // =================================================================
  void MallocGrid();

  // =================================================================
  //                          Protected Attributes
  // =================================================================
  std::shared_ptr<JumpingMT> prng_ = nullptr;

  std::shared_ptr<JumpingMT> mut_prng_ = nullptr;
  std::shared_ptr<JumpingMT> stoch_prng_ = nullptr;
  std::unique_ptr<JumpingMT> stoch_prng_bak_ = nullptr;

  int16_t width_  = -1;
  int16_t height_ = -1;

  int16_t x_best = -1;
  int16_t y_best = -1;

  SaveGridCell*** grid_ = nullptr;

  bool is_well_mixed_ = false;
  int32_t partial_mix_nb_permutations_ = 0;

  bool phenotypic_target_shared_ = true;
  std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler_ = nullptr;


  double  secretion_diffusion_prop_ = -1;
  double  secretion_degradation_prop_ = -1;


};


// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol
#endif // AEVOL_SAVEWORLD_H_
