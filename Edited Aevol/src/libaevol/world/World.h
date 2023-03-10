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


#ifndef AEVOL_WORLD_H_
#define AEVOL_WORLD_H_


// ============================================================================
//                                   Includes
// ============================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <memory>

#include <zlib.h>

#include "GridCell.h"

#include "SaveWorld.h"

namespace aevol {



// =================================================================
//                          Class declarations
// =================================================================
class ExpManager;






class World
{
    friend class GridCell;
 public :
  // =================================================================
  //                             Constructors
  // =================================================================
  World() = default;
  World(const World &) = delete;

  // =================================================================
  //                             Destructors
  // =================================================================
  virtual ~World();

  // =================================================================
  //                        Accessors: getters
  // =================================================================
  // PRNGs
  std::shared_ptr<JumpingMT> prng() const;
  std::shared_ptr<JumpingMT> mut_prng() const;
  std::shared_ptr<JumpingMT> stoch_prng() const;
  

  std::shared_ptr<JumpingMT> non_const_prng();
  std::shared_ptr<JumpingMT> non_const_mut_prng();
  std::shared_ptr<JumpingMT> non_const_stoch_prng();

  std::list<Individual *> indivs() const;
  inline int32_t          nb_indivs() const;
  inline Individual *   best_indiv() const;
  int16_t          width()  const {return width_;};
  int16_t          height() const {return height_;};
  inline int32_t          partial_mix_nb_permutations() const;
  GridCell ***  grid() const {return grid_;};
  inline GridCell*    grid(int16_t x, int16_t y) const;
  inline Individual*   indiv_at(int16_t x, int16_t y) const;

    inline Individual*   indiv_at_nonconst(int16_t x, int16_t y) {return grid_[x][y]->individual_nonconst();}

  Individual* indiv_by_id(int32_t id) const;
  Individual* indiv_by_rank(int32_t rank) const;

  double secretion_degradation_prop() const { return secretion_degradation_prop_; };
  double secretion_diffusion_prop() const { return secretion_diffusion_prop_; };

  inline double** secretion_present_grid() const;
  inline double** secreted_amount_grid() const;
  inline double** metabolic_fitness_grid() const;
  inline double** total_fitness_grid() const;

  bool phenotypic_target_shared() const {
    return phenotypic_target_shared_;
  }

  std::shared_ptr<PhenotypicTargetHandler>
  phenotypic_target_handler() const {
    return phenotypic_target_handler_;
  }

  // =================================================================
  //                        Accessors: setters
  // =================================================================
  // PRNGs
  inline void set_prng(std::shared_ptr<JumpingMT> prng);
  void set_mut_prng(std::shared_ptr<JumpingMT> prng);
  void set_stoch_prng(std::shared_ptr<JumpingMT> prng);

  inline void set_is_well_mixed(bool is_well_mixed);
  inline void set_partial_mix_nb_permutations(int32_t nb_permutations);
  inline void set_secretion_degradation_prop(double degradation_prop);
  inline void set_secretion_diffusion_prop(double diffusion_prop);
  inline void set_best(int16_t x, int16_t y);

  // =================================================================
  //                              Operators
  // =================================================================

  // =================================================================
  //                            Public Methods
  // =================================================================
  void InitGrid(int16_t width, int16_t height,
                Habitat& habitat,
                bool share_phenotypic_target);
  void PlaceIndiv(Individual * indiv, int16_t x, int16_t y, bool set_prng);
  void FillGridWithClones(Individual & dolly);
  void evaluate_individuals();
  void update_secretion_grid();
  void MixIndivs();
  void update_best();
  void ApplyHabitatVariation();

  void save(gzFile backup_file) const;
  void load(gzFile backup_file, ExpManager * exp_man);

  void set_phen_target_prngs(std::shared_ptr<JumpingMT> var_prng,
                             std::shared_ptr<JumpingMT> noise_prng);

  SaveWorld* make_save(ExpManager* exp_m, std::list<Individual*> indivs, bool share_phenotypic_target);
  SaveWorld* make_save(ExpManager* exp_m,
                       Individual_7** indivs,
                       Individual_7* best_indiv);

  int16_t x_best = -1;
  int16_t y_best = -1;
 protected :
  // =================================================================
  //                           Protected Methods
  // =================================================================
  void MallocGrid();
  void WellMixIndivs();
  void PartiallyMixIndivs();
  void backup_stoch_prng();

  // =================================================================
  //                          Protected Attributes
  // =================================================================
  std::shared_ptr<JumpingMT> prng_ = nullptr;

  std::shared_ptr<JumpingMT> mut_prng_ = nullptr;
  std::shared_ptr<JumpingMT> stoch_prng_ = nullptr;
  std::unique_ptr<JumpingMT> stoch_prng_bak_ = nullptr;

  int16_t width_  = -1;
  int16_t height_ = -1;


  GridCell*** grid_ = nullptr;
  GridCell** grid_1d_ = nullptr;

  bool is_well_mixed_ = false;
  int32_t partial_mix_nb_permutations_ = 0;

  bool phenotypic_target_shared_ = true;
  std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler_;

  double  secretion_diffusion_prop_ = -1;
  double  secretion_degradation_prop_ = -1;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================
inline int32_t World::nb_indivs() const
{
  return width_ * height_;
}

inline Individual *World::best_indiv() const
{
  return grid_[x_best][y_best]->individual();
}

inline int32_t World::partial_mix_nb_permutations() const
{
  return partial_mix_nb_permutations_;
}

inline GridCell *World::grid(int16_t x, int16_t y) const
{
  return grid_[x][y];
}

inline Individual *World::indiv_at(int16_t x, int16_t y) const
{
  return grid_[x][y]->individual();
}

inline double**World::secretion_present_grid() const
{
  double** ret = new double*[width_];

  for (int16_t x = 0; x < width_ ; x++)
  {
    ret[x] = new double[height_];
    for (int16_t y = 0; y < height_ ; y++)
    {
      ret[x][y] = grid_[x][y]->compound_amount();
    }
  }

  return ret;
}

inline double**World::secreted_amount_grid() const
{
  double** ret = new double*[width_];
  for (int16_t x = 0 ; x < width_ ; x++)
  {
    ret[x] = new double[height_];
    for (int16_t y = 0; y < height_ ; y++)
    {
      ret[x][y] = grid_[x][y]->secreted_amount();
    }
  }

  return ret;
}

inline double**World::metabolic_fitness_grid() const
{
  double** ret = new double*[width_];
  for (int16_t x = 0 ; x < width_ ; x++)
  {
    ret[x] = new double[height_];
    for (int16_t y = 0; y < height_ ; y++)
    {
      ret[x][y] = grid_[x][y]->metabolic_fitness();
    }
  }

  return ret;
}

inline double**World::total_fitness_grid() const
{
  double** ret = new double*[width_];
  for (int16_t x = 0 ; x < width_ ; x++)
  {
    ret[x] = new double[height_];
    for (int16_t y = 0; y < height_ ; y++)
    {
      ret[x][y] = grid_[x][y]->total_fitness();
    }
  }

  return ret;
}

// =====================================================================
//                           Setters' definitions
// =====================================================================
inline void World::set_prng(std::shared_ptr<JumpingMT> prng) {
  prng_ = prng;
}

inline void World::set_is_well_mixed(bool is_well_mixed) {
  is_well_mixed_ = is_well_mixed;
}

inline void World::set_partial_mix_nb_permutations(int32_t nb_permutations) {
  partial_mix_nb_permutations_ = nb_permutations;
}

inline void World::set_secretion_degradation_prop(double degradation_prop) {
  secretion_degradation_prop_=degradation_prop;
}

inline void World::set_secretion_diffusion_prop(double diffusion_prop) {
  secretion_diffusion_prop_=diffusion_prop;
}

inline void World::set_best(int16_t x, int16_t y) {
  x_best = x;
  y_best = y;
}


// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol
#endif // AEVOL_WORLD_H_
