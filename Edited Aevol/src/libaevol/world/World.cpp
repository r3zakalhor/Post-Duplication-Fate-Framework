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
//                              Includes
// =================================================================
#include "World.h"
#include "HabitatFactory.h"
#include "ExpManager.h"

#include "ExpManager_7.h"

#if __cplusplus == 201103L
#include "make_unique.h"
#endif

#include <list>
#include <iostream>

#ifdef __REGUL
#include "raevol/Individual_R.h"
#endif

using std::cout;
using std::endl;
using std::list;


namespace aevol {


//##############################################################################
//                                                                             #
//                                Class World                                  #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================

// =================================================================
//                             Destructor
// =================================================================
World::~World()
{
/*  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++) {
      for (int16_t x2 = x ; x2 < width_ ; x2++)
        for (int16_t y2 = 0 ; y2 < height_ ; y2++) {
          if (grid_[x][y]->individual_ != nullptr && grid_[x2][y2]->individual_ != nullptr)
          if (grid_[x][y]->individual_->id() == grid_[x2][y2]->individual_->id())
            grid_[x2][y2]->individual_ = nullptr;
    }
  }*/

  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++)
      delete grid_[x][y];

  // grid_ is 2D accessible but 1D allocated, there were only 2 new
  // statements and these are the corresponding deletes
  delete [] grid_1d_;
  delete [] grid_;
}

// =================================================================
//                            Public Methods
// =================================================================
void World::InitGrid(int16_t width, int16_t height,
                     Habitat& habitat,
                     bool share_phenotypic_target) {
  assert(share_phenotypic_target);

  if (share_phenotypic_target) {
    #ifndef __REGUL
    phenotypic_target_handler_ = std::make_shared<PhenotypicTargetHandler>(habitat.phenotypic_target_handler());
    #else
    phenotypic_target_handler_ = std::make_shared<PhenotypicTargetHandler_R>((dynamic_cast<Habitat_R&>(habitat)).phenotypic_target_handler());
    #endif
  }

  width_  = width;
  height_ = height;

  MallocGrid();

  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++)
    {
      int mut_seed = mut_prng_->random(1000000);
      int stoch_seed = stoch_prng_->random(1000000);
      // printf("Seed Grid C %d : %d %d\n",x*height+y,mut_seed,stoch_seed);

      if (share_phenotypic_target)
        grid_[x][y] =
                new GridCell(x, y,
#ifdef __REGUL
                             HabitatFactory::create_unique_habitat(dynamic_cast<Habitat_R&>(habitat),share_phenotypic_target),
#else
                             HabitatFactory::create_unique_habitat(habitat,share_phenotypic_target),
#endif
                             NULL,std::make_shared<JumpingMT>(mut_seed),
                             std::make_shared<JumpingMT>(stoch_seed));
    }
}

void World::MallocGrid()
{
  // Although grid_ is a 2D array, we want all its cells to be contiguous
  // in memory. However, we also want it to be 2D-accessible i.e. we want to
  // be able to access a cell with grid_[x][y].
  // The following code does just this
  grid_1d_ = new GridCell * [width_ * height_];
  grid_ = new GridCell ** [width_];
  for (int16_t x = 0 ; x < width_ ; x++)
    grid_[x] = &(grid_1d_[x * height_]);
}

void World::PlaceIndiv(Individual * indiv, int16_t x, int16_t y, bool set_prng) {
  grid_[x][y]->set_individual(indiv);
  if (set_prng) {
    indiv->set_mut_prng(grid_[x][y]->mut_prng());
    indiv->set_stoch_prng(grid_[x][y]->stoch_prng());
  }
  indiv->set_id(x*height_+y);
}

void World::FillGridWithClones(Individual & dolly)
{
  int32_t id_new_indiv = 0;
  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++)
      #ifndef __REGUL
      PlaceIndiv(Individual::CreateClone(&dolly, id_new_indiv++), x, y, true);
      #else
      PlaceIndiv(Individual_R::CreateClone(dynamic_cast<Individual_R*>(&dolly), id_new_indiv++), x, y, true);
      #endif
}

void World::evaluate_individuals()
{
  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++) {
      #ifndef __REGUL
      Individual* indiv       = indiv_at_nonconst(x, y);
      #else
      Individual_R* indiv       = dynamic_cast <Individual_R*> (indiv_at(x, y));
      #endif

      indiv->Evaluate();
      indiv->compute_statistical_data();
    }
}

void World::update_secretion_grid()
{
  int16_t cur_x, cur_y;

  double ** new_secretion = new double*[width_];
  for (int16_t x = 0 ; x < width_ ; x++)
  {
    new_secretion[x] = new double[height_];
    for (int16_t y = 0 ; y < height_ ; y++)
    {
      new_secretion[x][y] = grid_[x][y]->compound_amount();
    }
  }

  for (int16_t x = 0 ; x < width_ ; x++)
  {
    for (int16_t y = 0 ; y < height_ ; y++)
    {
      // look at the entire neighborhood
      for (int8_t i = -1 ; i < 2 ; i++)
      {
        for (int8_t j = -1 ; j < 2 ; j ++)
        {
          cur_x = (x + i + width_)  % width_;
          cur_y = (y + j + height_) % height_;

          // add the diffusion from the neighboring cells
          new_secretion[x][y] +=
              grid_[cur_x][cur_y]->compound_amount() *
                  secretion_diffusion_prop_;
        }
      }
    }
  }

  // Substract what has diffused from each cell, and calculate the
  // compound degradation
  for (int16_t x = 0 ; x < width_ ; x++)
  {
    for (int16_t y = 0 ; y < height_ ; y++)
    {
      grid_[x][y]->set_compound_amount(new_secretion[x][y] -
          9 * grid_[x][y]->compound_amount() * secretion_diffusion_prop_);
      grid_[x][y]->set_compound_amount(grid_[x][y]->compound_amount() *
          (1 - secretion_degradation_prop_));
    }
  }
  for (int16_t x = 0 ; x < width_ ; x++)
  {
    delete [] new_secretion[x];
  }
  delete [] new_secretion;
}

/*
 * Perform mixing of individuals
 *
 * Depending on parameters this can either well-mix or partially-mix
 * the population
 */
void World::MixIndivs()
{
  if (is_well_mixed_)
    WellMixIndivs();
  else if (partial_mix_nb_permutations_ > 0)
    PartiallyMixIndivs();
}

/*
 * Suffle individuals randomly using Fisher-Yates shuffle
 */
void World::WellMixIndivs()
{
  for (int16_t i = width_ * height_ - 1 ; i > 0 ; i--) {
    int16_t j = prng_->random(i + 1); // random in [0, 1]

    // Swap individuals btw cells i and j
    Individual * tmp = grid_1d_[i]->individual();
    grid_1d_[i]->set_individual(grid_1d_[j]->individual());
    grid_1d_[j]->set_individual(tmp);
    tmp->set_mut_prng(grid_1d_[j]->mut_prng());
    tmp->set_stoch_prng(grid_1d_[j]->stoch_prng());
  }
}

/*
 * Perform permutations between individuals randomly
 *
 * The number of permutations is given by partial_mix_nb_permutations_
 */
void World::PartiallyMixIndivs()
{
  for (int32_t i = 0 ; i < partial_mix_nb_permutations_ ; i++)
  {
    int16_t old_x = prng_->random(width_);
    int16_t old_y = prng_->random(height_);
    int16_t new_x = prng_->random(width_);
    int16_t new_y = prng_->random(height_);

    // Swap the individuals in these grid cells...
    Individual * tmp_swap = grid_[old_x][old_y]->individual();
    grid_[old_x][old_y]->set_individual(grid_[new_x][new_y]->individual());
    grid_[new_x][new_y]->set_individual(tmp_swap);

    grid_[old_x][old_y]->individual()->set_mut_prng(grid_[old_x][old_y]->mut_prng());
    grid_[old_x][old_y]->individual()->set_stoch_prng(grid_[old_x][old_y]->stoch_prng());

    grid_[new_x][new_y]->individual()->set_mut_prng(grid_[new_x][new_y]->mut_prng());
    grid_[new_x][new_y]->individual()->set_stoch_prng(grid_[new_x][new_y]->stoch_prng());
  }
}

void World::update_best()
{
  x_best = y_best = 0;
  double fit_best = indiv_at(0, 0)->fitness();
  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++)
      if (indiv_at(x, y)->fitness() > fit_best) {
        x_best = x;
        y_best = y;
        fit_best = indiv_at(x, y)->fitness();
      }

  //printf("Position %d / %d :: %e / %d\n",x_best,y_best,fit_best, indiv_at(x_best, y_best)->genetic_unit_seq_length(0));
}

void World::ApplyHabitatVariation() {
  if (phenotypic_target_shared_)
    phenotypic_target_handler_->ApplyVariation();
  else
    for (int16_t x = 0 ; x < width_ ; x++)
      for (int16_t y = 0 ; y < height_ ; y++)
        grid_[x][y]->ApplyHabitatVariation();
}

void World::save(gzFile backup_file) const
{
  if (prng_ == nullptr)
  {
    printf("%s:%d: error: PRNG not initialized.\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  prng_->save(backup_file);

  mut_prng_->save(backup_file);

  int8_t tmp_with_stoch = static_cast<int8_t>(stoch_prng_ == nullptr ? 0 : 1);
  gzwrite(backup_file, &tmp_with_stoch, sizeof(tmp_with_stoch));
  if (tmp_with_stoch)
  {
    stoch_prng_->save(backup_file);
  }
  if (grid_ == nullptr)
  {
    printf("%s:%d: error: grid not initialized.\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  // Manage shared or private phenotypic targets
  int8_t tmp_phenotypic_target_shared =
      static_cast<int8_t>(phenotypic_target_shared_ ? 1 : 0);
  gzwrite(backup_file,
          &tmp_phenotypic_target_shared,
          sizeof(tmp_phenotypic_target_shared));
  if (phenotypic_target_shared_) {
    phenotypic_target_handler_->save(backup_file);
  }


  gzwrite(backup_file, &width_,   sizeof(width_));
  gzwrite(backup_file, &height_,  sizeof(height_));

  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++)
      grid_[x][y]->save(backup_file, phenotypic_target_shared_);

  gzwrite(backup_file, &x_best, sizeof(x_best));
  gzwrite(backup_file, &y_best, sizeof(y_best));

  gzwrite(backup_file, &is_well_mixed_,
          sizeof(is_well_mixed_));
  gzwrite(backup_file, &partial_mix_nb_permutations_,
          sizeof(partial_mix_nb_permutations_));
  gzwrite(backup_file, &secretion_diffusion_prop_,
          sizeof(secretion_diffusion_prop_));
  gzwrite(backup_file, &secretion_degradation_prop_,
          sizeof(secretion_degradation_prop_));
}

void World::load(gzFile backup_file, ExpManager * exp_man)
{
  // Retrieve PRNGs
  prng_ = std::make_shared<JumpingMT>(backup_file);
  mut_prng_ = std::make_shared<JumpingMT>(backup_file);

  int8_t tmp_with_stoch;
  gzread(backup_file, &tmp_with_stoch, sizeof(tmp_with_stoch));
  if (tmp_with_stoch)
  {
    stoch_prng_ = std::make_shared<JumpingMT>(backup_file);
  }

  // Manage shared or private phenotypic targets
  int8_t tmp_phenotypic_target_shared;
  gzread(backup_file,
          &tmp_phenotypic_target_shared,
          sizeof(tmp_phenotypic_target_shared));
  phenotypic_target_shared_ = tmp_phenotypic_target_shared;
  if (phenotypic_target_shared_) {
    phenotypic_target_handler_ =
    #ifndef __REGUL
        std::make_shared<PhenotypicTargetHandler>(backup_file);
    #else
        std::make_shared<PhenotypicTargetHandler_R>(backup_file);
    #endif
  }

  //printf("Phenotypic targets %ld\n", dynamic_cast<PhenotypicTargetHandler_R*>(phenotypic_target_handler_)->phenotypic_targets().size());
  // A priori useless car déjà fait dans le constructeur de reprise sur backup
  //phenotypic_target_handler_->BuildPhenotypicTarget();

  gzread(backup_file, &width_,  sizeof(width_));
  gzread(backup_file, &height_, sizeof(height_));


  MallocGrid();

  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++) {
        grid_[x][y] = new GridCell(backup_file,
                                   exp_man,
                                   phenotypic_target_handler_);
      //  printf("Individual %d %d %p %p %p\n",x,y,grid_[x][y],grid_[x][y]->individual_nonconst(),
      //         grid_[x][y]->individual_nonconst()->grid_cell());
    }

  gzread(backup_file, &x_best, sizeof(x_best));
  gzread(backup_file, &y_best, sizeof(y_best));

  gzread(backup_file, &is_well_mixed_,
         sizeof(is_well_mixed_));
  gzread(backup_file, &partial_mix_nb_permutations_,
         sizeof(partial_mix_nb_permutations_));
  gzread(backup_file, &secretion_diffusion_prop_,
         sizeof(secretion_diffusion_prop_));
  gzread(backup_file, &secretion_degradation_prop_,
         sizeof(secretion_degradation_prop_));
}


    SaveWorld* World::make_save(ExpManager* exp_m, std::list<Individual*> indivs, bool share_phenotypic_target) {
      SaveWorld* backup = new SaveWorld();
      //random generator
      backup->prng_       = std::make_shared<JumpingMT>(*prng_);
      backup->mut_prng_   = std::make_shared<JumpingMT>(*mut_prng_);
      backup->stoch_prng_ = std::make_shared<JumpingMT>(*stoch_prng_);

      //non Pointer
      backup->width_  = width_;
      backup->height_ = height_;
      backup->x_best  = x_best;
      backup->y_best  = y_best;
      backup->is_well_mixed_ = is_well_mixed_;
      backup->partial_mix_nb_permutations_ = partial_mix_nb_permutations_;
      backup->secretion_degradation_prop_ = secretion_degradation_prop_;
      backup->secretion_diffusion_prop_ = secretion_diffusion_prop_;

      //pointers
      backup->phenotypic_target_handler_ = phenotypic_target_handler_;

      //the grid
      backup->MallocGrid();
      for (int16_t x = 0 ; x < width_ ; x++)
        for (int16_t y = 0 ; y < height_ ; y++) {
          Individual* indiv = indivs.front();
          backup->grid_[x][y] = new SaveGridCell(exp_m,x, y,
#ifdef __REGUL
                                                 std::make_unique<Habitat_R>(grid_[x][y]->habitat(), true),
#else
                  std::make_unique<Habitat>(grid_[x][y]->habitat(), true),
#endif
                                                 indiv,
                                                 std::make_shared<JumpingMT>(*grid_[x][y]->mut_prng()),
                                                 std::make_shared<JumpingMT>(*grid_[x][y]->stoch_prng()),
                                                 std::make_shared<JumpingMT>(*grid_[x][y]->reprod_prng_),
                                                 std::make_shared<JumpingMT>(*grid_[x][y]->reprod_prng_simd_),
                                                         indiv->w_max());
          indivs.pop_front();
        }

      return backup;
    }

    SaveWorld* World::make_save(ExpManager* exp_m,
                                Individual_7** indivs,
                                Individual_7* best_indiv) {
      SaveWorld *backup = new SaveWorld();
      //random generator
      backup->prng_ = std::make_shared<JumpingMT>(*prng_);
      backup->mut_prng_ = std::make_shared<JumpingMT>(*mut_prng_);
      backup->stoch_prng_ = std::make_shared<JumpingMT>(*stoch_prng_);

      //non Pointer
      backup->width_ = width_;
      backup->height_ = height_;
      backup->x_best = best_indiv->indiv_id / height_;
      backup->y_best = best_indiv->indiv_id % height_;
      backup->is_well_mixed_ = is_well_mixed_;
      backup->partial_mix_nb_permutations_ = partial_mix_nb_permutations_;
      backup->secretion_degradation_prop_ = secretion_degradation_prop_;
      backup->secretion_diffusion_prop_ = secretion_diffusion_prop_;

      //pointers
      backup->phenotypic_target_handler_ = phenotypic_target_handler_;

      //the grid
      backup->MallocGrid();
        for (int16_t x = 0; x < width(); x++) {
            for (int16_t y = 0; y < height(); y++) {

                if (ExpManager_7::standalone_simd) {
                    backup->grid_[x][y] = new SaveGridCell(exp_m, x, y,
#ifdef __REGUL
                                                           std::make_unique<Habitat_R>(grid_[x][y]->habitat(), true),
#else
                                                           std::make_unique<Habitat>(grid_[x][y]->habitat(), true),
#endif
                                                           indivs[x*height()+y],
                                                           std::make_shared<JumpingMT>(*grid_[x][y]->mut_prng()),
                                                           std::make_shared<JumpingMT>(*grid_[x][y]->stoch_prng()),
                                                           std::make_shared<JumpingMT>(*grid_[x][y]->reprod_prng_),
                                                           std::make_shared<JumpingMT>(*grid_[x][y]->reprod_prng_simd_),
                                                           indivs[x*height()+y]->w_max_);
                } else {
                    backup->grid_[x][y] = new SaveGridCell(exp_m, x, y,
#ifdef __REGUL
                                                           std::make_unique<Habitat_R>(grid_[x][y]->habitat(), true),
#else
                            std::make_unique<Habitat>(grid_[x][y]->habitat(), true),
#endif
                                                           indivs[x*height()+y],
                                                           std::make_shared<JumpingMT>(*grid_[x][y]->mut_prng()),
                                                           std::make_shared<JumpingMT>(*grid_[x][y]->stoch_prng()),
                                                           std::make_shared<JumpingMT>(*grid_[x][y]->reprod_prng_),
                                                           std::make_shared<JumpingMT>(*grid_[x][y]->reprod_prng_simd_),
                                                           indivs[x*height()+y]->w_max_);
                }
            }
        }

/*      printf("Backup pointer : %p\n",backup);
        gzFile world_file = gzopen("toto.aex", "w");
        backup->save(world_file);
        gzclose(world_file);*/
      return backup;
    }

// =================================================================
//                           Protected Methods
// =================================================================
void World::backup_stoch_prng()
{
  // Store a copy of stoch_prng_ in stoch_prng_bak_
#if __cplusplus == 201103L
  stoch_prng_bak_ = make_unique<JumpingMT>(*stoch_prng_);
#else
  stoch_prng_bak_ = std::make_unique<JumpingMT>(*stoch_prng_);
#endif
}

// =================================================================
//                          Non inline accessors
// =================================================================
std::shared_ptr<JumpingMT> World::prng() const
{
  return prng_;
}

std::shared_ptr<JumpingMT> World::mut_prng() const
{
  return mut_prng_;
}

std::shared_ptr<JumpingMT> World::stoch_prng() const
{
  return stoch_prng_;
}

std::shared_ptr<JumpingMT> World::non_const_prng()
{
  return prng_;
}

std::shared_ptr<JumpingMT> World::non_const_mut_prng()
{
  return mut_prng_;
}

std::shared_ptr<JumpingMT> World::non_const_stoch_prng()
{
  return stoch_prng_;
}


list<Individual *> World::indivs() const
{
  list<Individual *> r;

  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++)
      r.push_back(indiv_at(x, y));

  return r;
}

void World::set_mut_prng(std::shared_ptr<JumpingMT> prng)
{
  mut_prng_ = prng;

  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++) {
      Individual * indiv;
      if ((indiv = indiv_at(x, y)))
        indiv->set_mut_prng(mut_prng_);
    }
}

void World::set_stoch_prng(std::shared_ptr<JumpingMT> prng)
{
  stoch_prng_ = prng;

  for (int16_t x = 0 ; x < width_ ; x++)
    for (int16_t y = 0 ; y < height_ ; y++) {
      Individual * indiv;
      if ((indiv = indiv_at(x, y)))
        indiv->set_stoch_prng(stoch_prng_);
    }
}

void World::set_phen_target_prngs(std::shared_ptr<JumpingMT> var_prng,
                                  std::shared_ptr<JumpingMT> noise_prng) {
  assert(phenotypic_target_shared_);
  phenotypic_target_handler_->set_var_prng(var_prng);
  phenotypic_target_handler_->set_noise_prng(noise_prng);
}

// This method is disabled because with clones the ID is not injective
// any longer. Use indiv_by_rank() or indiv_at() instead.
Individual* World::indiv_by_id(int32_t id) const {
  assert(false);
  return nullptr;
  Individual* indiv = grid_1d_[id]->individual();
  // When the population isn't mixed at all, the individual with id n is in
  // grid_1d_[n]. Try this first...
  if ((int32_t)indiv->id() == id)
    return indiv;
  // ... If it isn't, do a basic search
  int32_t nb_indivs = width_ * height_;
  for (int32_t i = 0 ; i < nb_indivs ; i++) {
    if ((int32_t ) grid_1d_[i]->individual()->id() == id)
      return grid_1d_[i]->individual();
  }
  return nullptr;
}

Individual* World::indiv_by_rank(int32_t rank) const {
  int32_t nb_indivs = width_ * height_;
  for (int32_t i = 0 ; i < nb_indivs ; i++) {
    if (grid_1d_[i]->individual()->rank() == rank)
      return grid_1d_[i]->individual();
  }
  return nullptr;
}
} // namespace aevol
