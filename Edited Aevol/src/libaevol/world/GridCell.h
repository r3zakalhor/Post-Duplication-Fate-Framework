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


#ifndef AEVOL_GRID_CELL_H_
#define AEVOL_GRID_CELL_H_


// =================================================================
//                              Libraries
// =================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <zlib.h>



// =================================================================
//                            Project Files
// =================================================================
#include "Individual.h"
#include "HabitatFactory.h"
#ifndef __REGUL
#include "Habitat.h"
#else
#include "raevol/Habitat_R.h"
#endif

namespace aevol {



// =================================================================
//                          Class declarations
// =================================================================
class ExpManager;





class GridCell
{
    friend class World;
 public :
  // =================================================================
  //                             Constructors
  // =================================================================
  GridCell() = delete;
  GridCell(const GridCell &) = delete;
  GridCell(int16_t x, int16_t y,
#ifdef __REGUL
               std::unique_ptr<Habitat_R>&& habitat,
#else
           std::unique_ptr<Habitat> habitat,
#endif
               Individual * indiv, std::shared_ptr<JumpingMT> mut_prng,
              std::shared_ptr<JumpingMT> stoch_prng);
  GridCell(gzFile backup_file,
               ExpManager * exp_m,
               std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler);

  // =================================================================
  //                             Destructors
  // =================================================================
  virtual ~GridCell();


  // =================================================================
  //                        Accessors: getters
  // =================================================================
  inline int16_t x() const {return x_;};
  inline int16_t y() const {return y_;};
  inline double compound_amount() const;
  inline Individual * individual() const;
  inline Individual * individual_nonconst() {return individual_;}

  inline double secreted_amount() const;
  inline double metabolic_fitness() const;
  inline double total_fitness() const;

#ifdef __REGUL
        const Habitat_R& habitat() const {
            return *habitat_;
        }
#else
  const Habitat& habitat() const {
    return *habitat_;
  }
#endif

#ifdef __REGUL
        Habitat_R& habitat_non_const() {
            return *habitat_;
        }
#else
    Habitat& habitat_non_const() {
      return *habitat_;
    }
#endif
    Habitat* habitat_ptr() {
      return habitat_.get();
    }

#ifndef __REGUL
  const PhenotypicTarget& phenotypic_target() const {
    return habitat_->phenotypic_target();
  }
#endif


  std::shared_ptr<JumpingMT> mut_prng() const;
  std::shared_ptr<JumpingMT> stoch_prng() const;
  // =================================================================
  //                        Accessors: setters
  // =================================================================
  inline void set_compound_amount(double compound_amount);
  inline void set_individual(Individual * indiv);
  inline void set_mut_prng(std::shared_ptr<JumpingMT> prng);
  inline void set_stoch_prng(std::shared_ptr<JumpingMT> prng);
  inline void set_reprod_prng(std::unique_ptr<JumpingMT>&& prng);
    inline void set_reprod_prng_simd(std::unique_ptr<JumpingMT>&& prng);
  // =================================================================
  //                            Public Methods
  // =================================================================
  void ApplyHabitatVariation();
  void save(gzFile backup_file,
            bool skip_phenotypic_target = false) const;
  void load(gzFile backup_file,
            ExpManager * exp_m,
            std::shared_ptr<PhenotypicTargetHandler> phenotypic_target_handler);

  std::unique_ptr<JumpingMT> reprod_prng_ = nullptr;
    std::unique_ptr<JumpingMT> reprod_prng_simd_ = nullptr;
    double *  probs;
    double *  local_fit_array;
    double *  local_meta_array;
    double    sum_local_fit;
    double*    loc_phenotype;
    int *  indiv_index;
    Individual* old_one;
 protected :
  // =================================================================
  //                           Protected Methods
  // =================================================================

  // =================================================================
  //                          Protected Attributes
  // =================================================================
  // Position on the grid
  int16_t x_;
  int16_t y_;

  // Pointer to the individual in this cell
  Individual* individual_ = NULL;

#ifdef __REGUL
        std::unique_ptr<Habitat_R> habitat_ = nullptr;
#else
        std::unique_ptr<Habitat> habitat_ = nullptr;
#endif


  std::shared_ptr<JumpingMT> mut_prng_ = nullptr;
  std::shared_ptr<JumpingMT> stoch_prng_ = nullptr;


};


// =====================================================================
//                           Getters' definitions
// =====================================================================
inline double GridCell::compound_amount() const
{
  return habitat_->compound_amount();
}

inline Individual *GridCell::individual() const
{
  return individual_;
}

inline double GridCell::secreted_amount() const
{
  return individual_->fitness_by_feature(SECRETION);
}

inline double GridCell::metabolic_fitness() const
{
  return individual_->fitness_by_feature(METABOLISM);
}

inline double GridCell::total_fitness() const
{
  return individual_->fitness();
}

// =====================================================================
//                           Setters' definitions
// =====================================================================
inline void GridCell::set_compound_amount(double compound_amount)
{
  habitat_->set_compound_amount(compound_amount);
}

inline void GridCell::set_individual(Individual * indiv)
{
  individual_ = indiv;
  if (individual_->grid_cell() != this)
  {
    individual_->set_grid_cell(this);
  }
}

inline void GridCell::set_mut_prng(std::shared_ptr<JumpingMT> prng) {
  mut_prng_ = prng;
  // printf("%d %d -- Update Mut PRNG\n",x(),y());
}

inline void GridCell::set_stoch_prng(std::shared_ptr<JumpingMT> prng) {
  stoch_prng_ = prng;
  // printf("%d %d -- Update Stoch PRNG\n",x(),y());
}


inline void GridCell::set_reprod_prng(std::unique_ptr<JumpingMT>&& prng)
{
  reprod_prng_ = std::move(prng);
}

    inline void GridCell::set_reprod_prng_simd(std::unique_ptr<JumpingMT>&& prng)
    {
      reprod_prng_simd_ = std::move(prng);
    }
// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol
#endif // AEVOL_GRID_CELL_H_
