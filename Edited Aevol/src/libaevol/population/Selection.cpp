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
#include "Selection.h"

#include "7/Dna_7.h"
#include "7/ExpManager_7.h"

#include "DnaMutator.h"
//#include <math.h>

#ifdef _OPENMP
#ifndef __OPENMP_GPU
#include <omp.h>
#endif
#endif

#if __cplusplus == 201103L
#include "make_unique.h"
#endif

#include<chrono>

#include <iostream>
#include <unordered_map>
using namespace std;
using namespace std::chrono;
#include "ae_logger.h"
#include "ExpManager.h"
#include "VisAVis.h"
#include "HybridFuzzy.h"
#ifdef __NO_X
  #ifndef __REGUL
    #include "Individual.h"
  #else
    #include "raevol/Individual_R.h"
  #endif
#elif defined __X11
  #ifndef __REGUL
    #include "Individual_X11.h"
  #else
    #include "raevol/Individual_R_X11.h"
  #endif
#endif
#include "7/Individual_7.h"

namespace aevol {


//##############################################################################
//
//                              Class Selection
//
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
Selection::Selection(ExpManager* exp_m) {
  exp_m_ = exp_m;

  // ----------------------------------------- Pseudo-random number generator
  //prng_ = NULL;

  // -------------------------------------------------------------- Selection
  selection_scheme_   = RANK_EXPONENTIAL;
  selection_pressure_ = 0.998;


        fitness_function_ = FITNESS_EXP;
        fitness_function_scope_x_ = 3;
        fitness_function_scope_y_ = 3;

  // --------------------------- Probability of reproduction of each organism
  prob_reprod_ = NULL;

#ifdef WITH_PERF_TRACES
        std::ofstream perf_traces_file_;
        perf_traces_file_.open("vanilla_perf_traces.csv",std::ofstream::trunc);
        perf_traces_file_<<"Generation,Indiv_ID,Runtime"<<std::endl;
        perf_traces_file_.close();
#endif

}

// =================================================================
//                             Destructors
// =================================================================
Selection::~Selection() {
  delete [] prob_reprod_;
}


int mutator = 0;

// =================================================================
//                            Public Methods
// =================================================================
void Selection::step_to_next_generation() {
  int32_t* nb_offsprings;
  // Create proxies
  World* world = exp_m_->world();
  int16_t grid_width  = world->width();
  int16_t grid_height = world->height();
  int32_t pop_size = grid_width*grid_height;

  GridCell*** pop_grid = exp_m_->grid();

  int16_t x, y;
  int8_t what;
  high_resolution_clock::time_point t1,t2;

  std::list<Individual*> new_generation;

#pragma omp single
  {
  // To create the new generation, we must create nb_indivs new individuals
  // (offspring) and "kill" the existing ones.
  // The number of offspring on a given individual will be given by a stochastic
  // process biased on it's fitness value (the selection process).
  // There are 3 possible selection schemes :
  //    * Linear Ranking
  //    * Exponential Ranking
  //    * Fitness proportionate
  //
  // Whichever method is chosen, we will
  // 1) Compute the probability of reproduction of each individual in the population
  // 2) Simulate the stochastic process by a multinomial drawing (based upon the probabilities computed in 1)
  // 3) Make the selected individuals reproduce, thus creating the new generation
  // 4) Replace the current generation by the newly created one.
  // 5) Sort the newly created population*
/*
  if (prng_ == NULL) {
    printf("%s:%d: error: PRNG not initialized.\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
*/
  // -------------------------------------------------------------------------------
  // 1) Compute the probability of reproduction of each individual in the population
  // -------------------------------------------------------------------------------
 #ifndef FIXED_POPULATION_SIZE
    #error this method is not ready for variable population size
    compute_local_prob_reprod();
  #else
    // The function compute_local_prob_reprod creates and fills the array prob_reprod_, which is telling us the probability of being picked for reproduction according to the rank of an individual in its neighboorhood.
    // It is only usefull when selection is rank based. When selection scheme is FITNESS_PROPORTIONATE, we do not need to call it.
    // It shoud only be called once in the simulation and not at each generation. So if prob_reprod_ already exists we do not need to call it.
    if ((selection_scheme_ != FITNESS_PROPORTIONATE) && (prob_reprod_ == NULL)) {
      compute_local_prob_reprod();
    }
  #endif
/*
  if (prng_ == NULL) {
    printf("%s:%d: error: PRNG not initialized.\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
*/

  // create a temporary grid to store the reproducers
  reproducers = new Individual ** [grid_width];
  for (int16_t i = 0 ; i < grid_width ; i++) {
    reproducers[i] = new Individual* [grid_height];
  }

  if (selection_scope_ == SCOPE_GLOBAL) {
    delete[] prob_reprod_;

    prob_reprod_ = new double[pop_size];

    double* fitnesses = new double[pop_size];
    double sum = 0;

    size_t i = 0;
    for (const auto& indiv: exp_m_->indivs()) {
      fitnesses[i] = indiv->fitness();
      sum += fitnesses[i];
      ++i;
    }

    for (int32_t i = 0; i < pop_size; i++) {
      prob_reprod_[i] = fitnesses[i] / sum;
    }

    delete[] fitnesses;

    nb_offsprings = new int32_t[pop_size];
    exp_m_->world()->grid(0,0)->reprod_prng_->multinomial_drawing(nb_offsprings, prob_reprod_, pop_size, pop_size);

    int index = 0;
    i = 0;

    for (const auto& indiv: exp_m_->indivs()) {
      for (int32_t j = 0; j < nb_offsprings[i]; j++) {
        x = index / grid_height;
        y = index % grid_height;

        reproducers[x][y] = indiv;

        index++;
      }
      i++;
    }
  } else {
    if (fitness_function_ == FITNESS_GLOBAL_SUM) {
#ifdef __REGUL
        int number_of_phenotypic_target_models = dynamic_cast<const Habitat_R&> (world->grid(0,0)->habitat()).number_of_phenotypic_target_models();

        fitness_sum_tab_ = new double[number_of_phenotypic_target_models];
      for (int env_id = 0; env_id < number_of_phenotypic_target_models; env_id++) {
          fitness_sum_tab_[env_id] = 0;
          for (int i = 0; i < exp_m_->world()->width(); i++)
            for (int j = 0; j < exp_m_->world()->height(); j++) {
              fitness_sum_tab_[env_id] +=
                  dynamic_cast<Individual_R*>(world->indiv_at(i, j))
                      ->fitness(env_id);
            }
        }
#else
        printf("Fitness local sum is not supported for Aevol (only R-Aevol)\n");
        exit(-1);
#endif
    }
  }
  }

    // Do local competitions
    #pragma omp for schedule(dynamic) private(x,y)
    for (int32_t index = 0; index < grid_width * grid_height; index++) {
      x = index / grid_height;
      y = index % grid_height;
      reproducers[x][y] = do_local_competition(x, y);
/*
      exp_m_->exp_m_7_->current_individuals[x*exp_m_->world()->height()+y] =
          new Individual_7(exp_m_,exp_m_->exp_m_7_->previous_individuals
          [reproducers[x][y]->grid_cell()->x()*exp_m_->world()->height()+reproducers[x][y]->grid_cell()->y()],false);

      exp_m_->exp_m_7_->current_individuals[x*exp_m_->world()->height()+y]->indiv_id = x*exp_m_->world()->height()+y;
      exp_m_->exp_m_7_->current_individuals[x*exp_m_->world()->height()+y]->parent_id =
          reproducers[x][y]->grid_cell()->x()*exp_m_->world()->height()+reproducers[x][y]->grid_cell()->y();*/

    }

std::unordered_map<unsigned long long, Individual*> unique_individual;

#pragma omp single
  {
    // TODO : Why is that not *after* the creation of the new population ?
    // Add the compound secreted by the individuals
    if (exp_m_->with_secretion()) {
      for (int16_t x = 0; x < grid_width; x++) {
        for (int16_t y = 0; y < grid_height; y++) {
          pop_grid[x][y]->set_compound_amount(
              pop_grid[x][y]->compound_amount() +
              pop_grid[x][y]->individual()->fitness_by_feature(SECRETION));
        }
      }

      // Diffusion and degradation of compound in the habitat
      world->update_secretion_grid();
    }

    // Create the new generation
    std::list<Individual*> old_generation = exp_m_->indivs();

#ifdef __DETECT_CLONE
    mutator = 0;
    for (auto indiv: old_generation) {
      indiv->number_of_clones_ = 0;
      unique_individual[indiv->id()] = indiv;
      (&indiv->genetic_unit_list().front())->dna()->set_hasMutate(false);
    }
#endif

  to_evaluate.clear();
  }



#pragma omp for schedule(dynamic)  private(x,y,what)
  for (int32_t index = 0; index < grid_width * grid_height; index++) {
    x = index / grid_height;
    y = index % grid_height;

    do_replication(reproducers[x][y],
                   x * grid_height + y, what, x, y);

#ifdef __DETECT_CLONE
    if (what == 1 || what == 2) {
#endif

#pragma omp critical(updateindiv)
      {
        to_evaluate.push_back(index);
      }
#ifdef __DETECT_CLONE
    }
#endif

  }

#pragma omp single
  {
    t1 = high_resolution_clock::now();
  }

#pragma omp for schedule(dynamic) private(x,y)
  for (int i = 0; i < (int) to_evaluate.size(); i++) {
    x = to_evaluate[i] / grid_height;
    y = to_evaluate[i] % grid_height;

    Individual* l_indiv = world->indiv_at(x,y);
#ifdef __REGUL
    if (!l_indiv->evaluated_) {
      run_life(dynamic_cast<Individual_R*>(l_indiv));
    }
#else
    run_life(l_indiv);
#endif

  }

#pragma omp single
  {
    for (int32_t index = 0; index < grid_width * grid_height; index++) {
      x = index / grid_height;
      y = index % grid_height;
      if (exp_m_->record_tree() || exp_m_->light_tree()) {
        {
          //EndReplicationEvent *eindiv = new EndReplicationEvent(
          //                    world->indiv_at(x, y), x, y);
          // Tell observers the replication is finished
          //->notifyObservers(END_REPLICATION, eindiv);
          // world->indiv_at(x, y)->compute_statistical_data();
          world->indiv_at(x, y)->compute_non_coding();
          //delete eindiv;
        }
      }
    }

    for (int32_t index = 0; index < grid_width * grid_height; index++) {
      x = index / grid_height;
      y = index % grid_height;
      if (exp_m_->record_tree() || exp_m_->light_tree()) {
        {
          //EndReplicationEvent *eindiv = new EndReplicationEvent(
          //                    world->indiv_at(x, y), x, y);
          // Tell observers the replication is finished
          //->notifyObservers(END_REPLICATION, eindiv);
//          world->indiv_at(x, y)->compute_non_coding();
#ifdef __REGUL
          exp_m_->tree()
              ->report_by_index(AeTime::time(), x * grid_height + y)
          ->signal_end_of_replication(dynamic_cast<Individual_R*>(world->indiv_at(x, y)));
#else
          exp_m_->tree()
              ->report_by_index(AeTime::time(), x * grid_height + y)
              ->signal_end_of_replication(world->indiv_at(x, y));
#endif
          //delete eindiv;
        }
      }
    }

    for (int16_t x = 0; x < grid_width; x++)
      for (int16_t y = 0; y < grid_height; y++)
        new_generation.push_back(pop_grid[x][y]->individual());

    // delete the temporary grid and the parental generation
    for (int16_t x = 0; x < grid_width; x++) {
      for (int16_t y = 0; y < grid_height; y++) {
        reproducers[x][y] = nullptr;
      }
      delete[] reproducers[x];
    }
    delete[] reproducers;

#ifndef __DETECT_CLONE
#ifdef __OPENMP_TASK
  #pragma omp parallel
  #pragma omp single
  {
#endif
  for (auto element = old_generation.begin();
       element != old_generation.end(); ++element) {
#ifdef __OPENMP_TASK
    #pragma omp task
      {
#endif
    delete *element;
#ifdef __OPENMP_TASK
    }
#endif
  }
#ifdef __OPENMP_TASK
  }
#endif
#endif

    // Compute the rank of each individual

    new_generation.sort([](Individual* lhs, Individual* rhs) {
      return lhs->fitness() < rhs->fitness();
    });

    int rank = 1;

    for (Individual* indiv: new_generation) {
      indiv->set_rank(rank++);
    }
    // randomly migrate some organisms, if necessary
    world->MixIndivs();

    PerformPlasmidTransfers();

    // Update the best individual
    exp_m_->update_best();

    // Notify observers of the end of the generation
    if (exp_m_->record_tree() || exp_m_->light_tree()) {
      exp_m_->tree()->signal_end_of_generation();
    }

#ifdef WITH_PERF_TRACES
    std::ofstream perf_traces_file_;
    perf_traces_file_.open("vanilla_perf_traces.csv", std::ofstream::app);
    for (int indiv_id = 0; indiv_id < (int)exp_m_->nb_indivs(); indiv_id++) {
      perf_traces_file_ << AeTime::time() << "," << indiv_id << ","
                        << apply_mutation[indiv_id] << std::endl;
    }
    perf_traces_file_.close();
#endif
  }
#ifdef __DETECT_CLONE
//    int number_of_clones = 0;
#ifdef __OPENMP_TASK
  #pragma omp parallel
  #pragma omp single
  {
#endif
    for (auto element = unique_individual.begin();
         element != unique_individual.end(); ++element) {
#ifdef __OPENMP_TASK
    #pragma omp task
      {
#endif
        if (element->second->number_of_clones_ == 0) {
          delete element->second;
        }
#ifdef __OPENMP_TASK
      }
#endif
    }
#ifdef __OPENMP_TASK
  }
#endif
#endif


}

void Selection::PerformPlasmidTransfers() {
  if (exp_m_->with_plasmids() &&
      ((exp_m_->prob_plasmid_HT() != 0.0) ||
        (exp_m_->tune_donor_ability() != 0.0) ||
        (exp_m_->tune_recipient_ability() != 0.0))) {

    // Create proxies
    World* world = exp_m_->world();
    int16_t grid_width  = world->width();
    int16_t grid_height = world->height();

    int16_t x_offset, y_offset, new_x, new_y;

    // Shuffle the grid:
    int16_t total_size = ((grid_width)*(grid_height));
    int16_t** shuffled_table = new int16_t* [total_size];
    for (int16_t z = 0 ; z < total_size ; z++) {
      shuffled_table[z] = new int16_t[2];
      int16_t quotient = z / grid_width;
      int16_t remainder = z % grid_width;
      shuffled_table[z][0] = (int16_t) remainder;
      shuffled_table[z][1] = (int16_t) quotient;
    }

    for (int16_t z = 0 ;z < total_size - 1 ; z++) {

      int16_t x_prng = z / grid_width;
      int16_t y_prng = z % grid_width;

      int16_t rand_nb = world->grid(x_prng,y_prng)->reprod_prng_->random((int16_t) (total_size-z));
      int16_t* tmp=shuffled_table[z+rand_nb];
      shuffled_table[z+rand_nb]=shuffled_table[z];
      shuffled_table[z]=tmp;
    }


    // First transfer all the plasmids, but just add them at the end of the list of the GUs
    for (int16_t z = 0 ; z < total_size ; z++) { // for each individual x
      int16_t x=shuffled_table[z][0];
      int16_t y=shuffled_table[z][1];

      for (int16_t n = 0 ; n < 9 ; n++) { // for each neighbour n of x
        x_offset = (n / 3) - 1;
        y_offset = (n % 3) - 1;

        new_x = (x+x_offset+grid_width) % grid_width;
        new_y = (y+y_offset+grid_height) % grid_height;

        if ((new_x != x)||(new_y != y)) {
          double ptransfer = exp_m_->prob_plasmid_HT() + exp_m_->tune_donor_ability()
                            * world->indiv_at(x, y)->fitness_by_feature(DONOR)
                            +
            exp_m_->tune_recipient_ability() * world->indiv_at(new_x, new_y)->fitness_by_feature(RECIPIENT) ;
          if (world->grid(x,y)->reprod_prng_->random() < ptransfer) { // will x give a plasmid to n ?
            if (exp_m_->swap_GUs()) {
              world->indiv_at(new_x, new_y)->inject_2GUs(world->indiv_at(x, y));
            }
            else {
              world->indiv_at(new_x, new_y)->inject_GU(world->indiv_at(x, y));
            }
          }
        }
      }
    }

    for(int16_t z=0;z <total_size;z++) {
      delete [] shuffled_table[z];
    }
    delete [] shuffled_table;



    // If an individual has more than 2 GUs, we keep only the first (main chromosome) and the last one
    // and re-evaluate the individual
    for (int16_t x = 0 ; x < grid_width ; x++) {
      for (int16_t y = 0 ; y < grid_height ; y++) {
        bool reevaluate = (world->indiv_at(x, y)->nb_genetic_units() > 2);
        world->indiv_at(x, y)->drop_nested_genetic_units();
        if (reevaluate)
          world->indiv_at(x, y)->Reevaluate();
      }
    }
  }
}

/*!
*/
void Selection::write_setup_file(gzFile exp_setup_file) const {
  // ---------------------------------------------------- Selection Parameters
  int8_t tmp_sel_scheme = static_cast<int8_t>(selection_scheme_);
  gzwrite(exp_setup_file, &tmp_sel_scheme,      sizeof(tmp_sel_scheme));
  gzwrite(exp_setup_file, &selection_pressure_, sizeof(selection_pressure_));

  int8_t tmp_sel_scope = static_cast<int8_t>(selection_scope_);
  gzwrite(exp_setup_file, &tmp_sel_scope,      sizeof(tmp_sel_scope));
  gzwrite(exp_setup_file, &selection_scope_x_, sizeof(selection_scope_x_));
  gzwrite(exp_setup_file, &selection_scope_y_, sizeof(selection_scope_y_));


        int8_t tmp_fit_func= fitness_function_;
        gzwrite(exp_setup_file, &tmp_fit_func,      sizeof(tmp_fit_func));
        gzwrite(exp_setup_file, &fitness_function_scope_x_, sizeof(fitness_function_scope_x_));
        gzwrite(exp_setup_file, &fitness_function_scope_y_, sizeof(fitness_function_scope_y_));
}

/*!
*/
void Selection::save(gzFile& backup_file) const {
  /*if (prng_ == NULL) {
    printf("%s:%d: error: PRNG not initialized.\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }*/

  // ----------------------------------------- Pseudo-random number generator
  //prng_->save(backup_file);
}

void Selection::load(gzFile& exp_setup_file,
                     gzFile& backup_file,
                     bool verbose) {
  // ---------------------------------------------------- Selection parameters
  int8_t tmp_sel_scheme;
  gzread(exp_setup_file, &tmp_sel_scheme, sizeof(tmp_sel_scheme));
  selection_scheme_ = static_cast<SelectionScheme>(tmp_sel_scheme);
  gzread(exp_setup_file, &selection_pressure_, sizeof(selection_pressure_));

  int8_t tmp_sel_scope;
  gzread(exp_setup_file, &tmp_sel_scope, sizeof(tmp_sel_scope));
  selection_scope_ = static_cast<SelectionScope>(tmp_sel_scope);
  gzread(exp_setup_file, &selection_scope_x_, sizeof(selection_scope_x_));
  gzread(exp_setup_file, &selection_scope_y_, sizeof(selection_scope_y_));


    int8_t tmp_fit_func;
    gzread(exp_setup_file, &tmp_fit_func, sizeof(tmp_fit_func));
    fitness_function_ = (FitnessFunction) tmp_fit_func;
    gzread(exp_setup_file, &fitness_function_scope_x_, sizeof(fitness_function_scope_x_));
    gzread(exp_setup_file, &fitness_function_scope_y_, sizeof(fitness_function_scope_y_));
  // ----------------------------------------- Pseudo-random number generator
  /*
#if __cplusplus == 201103L
  prng_ = make_unique<JumpingMT>(backup_file);
#else
  prng_ = std::make_unique<JumpingMT>(backup_file);
#endif*/
}


// =================================================================
//                           Protected Methods
// =================================================================
void Selection::compute_prob_reprod() { // non spatially structured only
  if (prob_reprod_ != NULL) { // TODO <david.parsons@inria.fr> remove
    delete [] prob_reprod_;
  }

  int32_t nb_indivs = exp_m_->nb_indivs();
  prob_reprod_ = new double[nb_indivs];

  if (selection_scheme_ == RANK_LINEAR) {
    // The probability of reproduction for an individual is given by
    // (2-SP + 2 * (SP-1) * (R-1)/(N-1)) / N
    // With :
    //      SP : selective pressure. Linear ranking allows values of SP in [1.0, 2.0].
    //      R  : the rank of the individual in the population (1 for the worst individual)
    //      N  : the number of individuals in the population
    //
    // We can transform this expression into (2-SP)/N + ((2*(SP-1)) / (N*(N-1))) * (R-1)
    // Furthermore, (R-1) is given directly by <i> (the index of our probability table)
    //
    // probs[0] will hence be given by (2-SP)/N
    // probs[i+1] can then be expressed by probs[i] + (2*(SP-1)) / (N*(N-1))

    double increment = (2 * (selection_pressure_-1)) / (nb_indivs * (nb_indivs-1));
    prob_reprod_[0]  = (2 - selection_pressure_) / nb_indivs;

    for (int32_t i = 1 ; i < nb_indivs ; i++) {
      prob_reprod_[i] = prob_reprod_[i-1] + increment;
    }

    // No need to normalize: The sum is always 1 for linear ranking
  }
  else if (selection_scheme_ == RANK_EXPONENTIAL) {
    // The probability of reproduction for an individual is given by
    // ((SP-1) * SP^(N-R)) / (SP^N - 1)
    // Which is equivalent to
    // ((SP-1) * SP^N) / ((SP^N - 1) * SP^R)
    // With :
    //      SP : selective pressure. Exponential ranking allows values of SP in ]0.0, 1.0[
    //      R  : the rank of the individual in the population (1 for the worst individual)
    //      N  : the number of individuals in the population
    //
    // NB : The only rank-dependent term is SP^R
    //
    // Because we don't allow ex-aequo,
    // probs[i+1] can hence be expressed as (probs[i] / SP)
    // We will hence compute probs[0] with the original formula and infer the remaining values

    double SP_N = pow(selection_pressure_, nb_indivs); // SP^N
    prob_reprod_[0] = ((selection_pressure_ - 1) * SP_N) /
                      ((SP_N - 1) * selection_pressure_);

    for (int32_t i = 1 ; i < nb_indivs ; i++) {
      prob_reprod_[i] = prob_reprod_[i-1] / selection_pressure_;
    }

    // No need to normalize: We don't allow ex-aequo
  }
  else if (selection_scheme_ == FITNESS_PROPORTIONATE) {
    // The probability of reproduction for an individual is given by
    // exp(-SP * gap) / sum of this measure on all individuals
    //    SP : selective pressure. Fitness proportionate allows values of SP in ]0, +inf[
    //                             The closer SP to 0, the closer the selection to being linear.

    double* fitnesses = new double[nb_indivs];
    double  sum       = 0;

    size_t i = 0;
    for (const auto& indiv: exp_m_->indivs()) {
      fitnesses[i] = indiv->fitness();
      sum += fitnesses[i];
      ++i;
    }

    for (int32_t i = 0 ; i < nb_indivs ; i++) {
      prob_reprod_[i] = fitnesses[i] / sum;
    }

    delete [] fitnesses;
  }
  else if (selection_scheme_ == FITTEST) {
    printf("ERROR, fittest selection scheme is meant to be used for spatially structured populations %s:%d\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  else {
    printf("ERROR, invalid selection scheme in file %s:%d\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
}

void Selection::compute_local_prob_reprod() {
  int16_t neighborhood_size = selection_scope_x_*selection_scope_y_;

  if (prob_reprod_ != NULL) {
    printf ("Warning, already defined %s:%d\n", __FILE__, __LINE__);
    delete [] prob_reprod_;
  }

  prob_reprod_ = new double[neighborhood_size];

  if (selection_scheme_ == RANK_LINEAR) {
    double increment = (2 * (selection_pressure_-1)) / (neighborhood_size * (neighborhood_size-1));
    double init_prob = (2 - selection_pressure_) / neighborhood_size;

    for (int16_t i = 0 ; i < neighborhood_size ; i++) {
      prob_reprod_[i] = init_prob + increment * i;
    }
  }
  else if (selection_scheme_ == RANK_EXPONENTIAL) {
    double SP_N = pow(selection_pressure_, neighborhood_size);
    prob_reprod_[0] = ((selection_pressure_ - 1) * SP_N) /
    ((SP_N - 1) * selection_pressure_);

    for (int16_t i = 1 ; i < neighborhood_size ; i++) {
      prob_reprod_[i] =  prob_reprod_[i-1] /  selection_pressure_;
    }
  }
  else if (selection_scheme_ == FITTEST) {
    for (int16_t i = 0 ; i < neighborhood_size-1 ; i++) {
      prob_reprod_[i] = 0.;
    }
    prob_reprod_[neighborhood_size-1] = 1.;
  }
  else if (selection_scheme_ == FITNESS_PROPORTIONATE) {
    printf("ERROR, this function is not intented to be use with this selection scheme %s:%d\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  else {
    printf("ERROR, invalid selection scheme in file %s:%d\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
}

Individual* Selection::do_replication(Individual* parent,
                                      unsigned long long index,
                                      int8_t& type_mutate,
                                      int16_t x,
                                      int16_t y) {
// ===========================================================================
//  1) Copy parent
// ===========================================================================
#ifdef __NO_X
  #ifndef __REGUL
  Individual* new_indiv =
      new Individual(parent, index, exp_m_->world()->grid(x, y)->mut_prng(),
                     exp_m_->world()->grid(x, y)->stoch_prng());
  #else
  Individual_R* new_indiv =
      new Individual_R(dynamic_cast<Individual_R*>(parent), index,
                       exp_m_->world()->grid(x, y)->mut_prng(),
                       exp_m_->world()->grid(x, y)->stoch_prng());
  #endif
#elif defined __X11
  #ifndef __REGUL
  Individual_X11* new_indiv =
      new Individual_X11(dynamic_cast<Individual_X11*>(parent), index,
                         exp_m_->world()->grid(x, y)->mut_prng(),
                         exp_m_->world()->grid(x, y)->stoch_prng());
  #else
  Individual_R_X11* new_indiv =
      new Individual_R_X11(dynamic_cast<Individual_R_X11*>(parent), index,
                           exp_m_->world()->grid(x, y)->mut_prng(),
                           exp_m_->world()->grid(x, y)->stoch_prng());
  #endif
#endif

  // Set the new individual's location on the grid
#pragma omp critical(placeindiv)
  {
    exp_m_->world()->PlaceIndiv(new_indiv, x, y, true);
  }
  //NewIndivEvent *eindiv = new NewIndivEvent(new_indiv,parent, x, y);
  //notifyObservers(NEW_INDIV, eindiv);
  delete exp_m_->dna_mutator_array_[x * exp_m_->world()->height() + y];

  exp_m_->dna_mutator_array_[x * exp_m_->world()->height() + y] =
      new DnaMutator(exp_m_->world()->grid(x, y)->mut_prng(),
                     parent->amount_of_dna(),
                     exp_m_->exp_s()->mut_params()->duplication_rate(),
                     exp_m_->exp_s()->mut_params()->deletion_rate(),
                     exp_m_->exp_s()->mut_params()->translocation_rate(),
                     exp_m_->exp_s()->mut_params()->inversion_rate(),
                     exp_m_->exp_s()->mut_params()->point_mutation_rate(),
                     exp_m_->exp_s()->mut_params()->small_insertion_rate(),
                     exp_m_->exp_s()->mut_params()->small_deletion_rate(),
                     exp_m_->exp_s()->mut_params()->max_indel_size(),
                     exp_m_->exp_s()->min_genome_length(),
                     exp_m_->exp_s()->max_genome_length(),
                     x * exp_m_->world()->height() + y, x, y);
  exp_m_->dna_mutator_array_[x * exp_m_->world()->height() + y]
      ->generate_mutations();

  // printf("Generate Mutator %d : %d\n", x * exp_m_->world()->height() + y, 
  //     exp_m_->dna_mutator_array_[x * exp_m_->world()->height() + y]->hasMutate());

  bool mutate = true;

  // Perform transfer, rearrangements and mutations
  if (not new_indiv->allow_plasmids()) {
    const GeneticUnit* chromosome = &new_indiv->genetic_unit_list().front();

#ifdef __DETECT_CLONE
    if (!exp_m_->dna_mutator_array_[x * exp_m_->world()->height() + y]
             ->hasMutate()) {

      bool firstClone;

  #pragma omp critical(firstclone)
      {
        firstClone = parent->number_of_clones_ == 0;

        firstClone ? type_mutate = 2 : type_mutate = 0;

        parent->number_of_clones_++;
      }
      delete new_indiv;
      mutate = false;
  #ifdef __NO_X
    #ifndef __REGUL
      new_indiv = parent;
    #else
      new_indiv = dynamic_cast<Individual_R*>(parent);
    #endif
  #elif defined __X11
    #ifndef __REGUL
      new_indiv = dynamic_cast<Individual_X11*>(parent);
    #else
      new_indiv = dynamic_cast<Individual_R_X11*>(parent);
    #endif
#endif
if (exp_m_->record_tree() || exp_m_->light_tree()) {
  #pragma omp critical(placeindiv)
        {
#ifdef __REGUL
          NewIndivEvent* eindiv = new NewIndivEvent(
              dynamic_cast<Individual_R*>(new_indiv), dynamic_cast<Individual_R*>(parent), x, y, index,
              exp_m_->next_generation_reproducer_[index]);
#else
          NewIndivEvent* eindiv = new NewIndivEvent(
              new_indiv, parent, x, y, index,
              exp_m_->next_generation_reproducer_[index]);

#endif
          //notifyObservers(NEW_INDIV, eindiv);
          exp_m_->tree()->update_new_indiv(eindiv);
          delete eindiv;
        }
      }

      // Notify observers that a new individual was created from <parent>
#pragma omp critical(placeindiv)
      {
        exp_m_->world()->PlaceIndiv(new_indiv, x, y, false);
      }

//      printf("%d -- CPU -- Indiv %d (Parent %d) has NOT mutate ? DNA Size %d\n",time(),x*exp_m_->world()->height()+y,new_indiv->parent_id_,
//             new_indiv->genetic_unit_seq_length(0));

#ifdef WITH_PERF_TRACES
        apply_mutation[index] = -1;
#endif
    } else {
#endif
    if (exp_m_->record_tree() || exp_m_->light_tree()) {
        {
#ifdef __REGUL
          NewIndivEvent* eindiv = new NewIndivEvent(
              dynamic_cast<Individual_R*>(new_indiv), dynamic_cast<Individual_R*>(parent), x, y, index,
              exp_m_->next_generation_reproducer_[index]);
#else
        NewIndivEvent* eindiv = new NewIndivEvent(
            new_indiv, parent, x, y, index,
            exp_m_->next_generation_reproducer_[index]);
#endif
          //notifyObservers(NEW_INDIV, eindiv);
          exp_m_->tree()->update_new_indiv(eindiv);
          delete eindiv;
        }
      }
#ifdef WITH_PERF_TRACES
      auto t_start = std::chrono::steady_clock::now();
#endif
    chromosome->dna()->apply_mutations();
#ifdef WITH_PERF_TRACES
      auto t_end = std::chrono::steady_clock::now();
      apply_mutation[index] = t_end.time_since_epoch().count() - t_start.time_since_epoch().count();
#endif

#ifdef __DETECT_CLONE
    }
#endif

    //#pragma omp critical(newindivevent)

  } else { // For each GU, apply mutations
    // Randomly determine the order in which the GUs will undergo mutations
    bool inverse_order =
        (exp_m_->world()->grid(x, y)->reprod_prng_->random((int32_t)2) < 0.5);

    if (not inverse_order) { // Apply mutations in normal GU order
      for (const auto& gen_unit: new_indiv->genetic_unit_list()) {
        gen_unit.dna()->perform_mutations(parent->id());
      }
    } else { // Apply mutations in inverse GU order
      const auto& gul = new_indiv->genetic_unit_list();
      for (auto gen_unit = gul.crbegin(); gen_unit != gul.crend(); ++gen_unit) {
        gen_unit->dna()->perform_mutations(parent->id());
      }
    }
  }

#ifdef __DETECT_CLONE
  if (mutate) {
  #ifdef __REGUL
    new_indiv->init_indiv();
  #else
    new_indiv->Evaluate();
  #endif
    type_mutate = 1;
  }
#endif


  return new_indiv;
}

#ifdef __REGUL
void Selection::run_life(Individual_R* new_indiv) {

    if (dynamic_cast<PhenotypicTargetHandler_R*>(&new_indiv->grid_cell()->habitat().
        phenotypic_target_handler_nonconst())->hasChanged()) {
      new_indiv->evaluated_ = false;
    }

    // Evaluate new individual
    new_indiv->Evaluate();

    // Compute statistics
    new_indiv->compute_statistical_data();

}
#else
void Selection::run_life(Individual* new_indiv) {
  // Evaluate new individual
  new_indiv->Evaluate();

  // Compute statistics
  new_indiv->compute_statistical_data();

}
#endif

Individual *Selection::do_local_competition (int16_t x, int16_t y) {
  // This function uses the array prob_reprod_ when selection scheme is
  // RANK_LINEAR, RANK_EXPONENTIAL, or FITTEST. For these selection schemes,
  // the function compute_local_prob_reprod (creating the array prob_reprod_)
  // must have been called before.
  // When selection scheme is FITNESS_PROPORTIONATE, this function only uses
  // the fitness values

  World* world = exp_m_->world();

  int16_t neighborhood_size = selection_scope_x_*selection_scope_y_;
  int16_t grid_width  = world->width();
  int16_t grid_height = world->height();
  int16_t cur_x;
  int16_t cur_y;

  // Build a temporary local array of fitness values
  double *  local_fit_array   = new double[neighborhood_size];
  //double *  local_meta_array   = new double[neighborhood_size];
  double *  sort_fit_array    = new double[neighborhood_size];
  int16_t * initial_location  = new int16_t[neighborhood_size];
  double *  probs             = new double[neighborhood_size];
  int16_t   count             = 0;
  double    sum_local_fit     = 0.0;
  //double* loc_phenotype = new double[300];

#ifdef __REGUL
  double ** fitness_sum_local_tab_;
  int number_of_phenotypic_target_models = dynamic_cast<const Habitat_R&> (world->grid(x,y)->habitat()).number_of_phenotypic_target_models();
#endif

  if (fitness_function_ == FITNESS_LOCAL_SUM) {
#ifdef __REGUL
    fitness_sum_local_tab_ = new double*[fitness_function_scope_x_*fitness_function_scope_y_];
    for (int tab_id = 0; tab_id < fitness_function_scope_x_*fitness_function_scope_y_; tab_id++)
      fitness_sum_local_tab_[tab_id] = new double[number_of_phenotypic_target_models];

    for (int env_id = 0; env_id < number_of_phenotypic_target_models; env_id++) {

      int tab_id = 0;
      fitness_sum_local_tab_[tab_id][env_id] = 0;

      for (int8_t i = -(selection_scope_x_/2) ; i <= (selection_scope_x_/2) ; i++) {
        for (int8_t j = -(selection_scope_y_/2) ; j <= (selection_scope_y_/2) ; j++) {
          cur_x = (x + i + grid_width) % grid_width;
          cur_y = (y + j + grid_height) % grid_height;

            int16_t new_x,new_y;
          for (int8_t ii = -(fitness_function_scope_x_/2); ii <= (fitness_function_scope_x_/2); ii++) {
            for (int8_t jj = -(fitness_function_scope_y_/2); jj <= (fitness_function_scope_y_/2); jj++) {
                //TODO: Check values HERE !

              new_x = (cur_x + ii + grid_width) % grid_width;
              new_y = (cur_y + jj + grid_height) % grid_height;

              fitness_sum_local_tab_[tab_id][env_id] +=  dynamic_cast<Individual_R*>(world->indiv_at(new_x, new_y))->fitness(env_id);
            }
          }

          tab_id++;
        }
      }
    }
#else
    printf("Fitness local sum is not supported for Aevol (only R-Aevol)\n");
    exit(-1);
#endif
  }

  int tab_id = 0;
  for (int8_t i = -(selection_scope_x_/2) ; i <= (selection_scope_x_/2) ; i++) {
    for (int8_t j = -(selection_scope_y_/2) ; j <= (selection_scope_y_/2) ; j++) {
      cur_x = (x + i + grid_width)  % grid_width;
      cur_y = (y + j + grid_height) % grid_height;

      if (fitness_function_ == FITNESS_EXP)
        local_fit_array[count]  = world->indiv_at(cur_x, cur_y)->fitness();
      else if (fitness_function_ == FITNESS_GLOBAL_SUM) {
#ifdef __REGUL
        double composed_fitness = 0;
        for (int env_id = 0; env_id < number_of_phenotypic_target_models; env_id++) {
          composed_fitness +=  dynamic_cast<Individual_R*>(world->indiv_at(cur_x, cur_y))->fitness(env_id) / fitness_sum_tab_[env_id];
        }
        composed_fitness/=number_of_phenotypic_target_models;
        local_fit_array[count]  = composed_fitness;
#else
          printf("Fitness local sum is not supported for Aevol (only R-Aevol)\n");
          exit(-1);
#endif
      } else if (fitness_function_ == FITNESS_LOCAL_SUM) {
#ifdef __REGUL
        double composed_fitness = 0;
        for (int env_id = 0; env_id < number_of_phenotypic_target_models; env_id++) {
          composed_fitness +=  dynamic_cast<Individual_R*>(world->indiv_at(cur_x, cur_y))->fitness(env_id) / fitness_sum_local_tab_[tab_id][env_id];
        }
        composed_fitness/=number_of_phenotypic_target_models;
        local_fit_array[count]  = composed_fitness;
#else
          printf("Fitness local sum is not supported for Aevol (only R-Aevol)\n");
          exit(-1);
#endif
      }

      local_fit_array[count]  = world->indiv_at(cur_x, cur_y)->fitness();
      //local_meta_array[count]  = world->indiv_at(cur_x, cur_y)->dist_to_target_by_feature(METABOLISM);
      sort_fit_array[count]   = local_fit_array[count];
      initial_location[count] = count;
      sum_local_fit += local_fit_array[count];
      count++;
      tab_id++;
    }
  }

  if (fitness_function_ == FITNESS_LOCAL_SUM) {
#ifdef __REGUL
    for (int tab_id = 0; tab_id < fitness_function_scope_x_ * fitness_function_scope_y_; tab_id++)
      delete[] fitness_sum_local_tab_[tab_id];
    delete[] fitness_sum_local_tab_;
#else
      printf("Fitness local sum is not supported for Aevol (only R-Aevol)\n");
      exit(-1);
#endif
  }
  //printf("Competition 2\n");
  // Do the competitions between the individuals, based on one of the 4 methods:
  // 1. Rank linear
  // 2. Rank exponential
  // 3. Fitness proportionate
  // 4. Fittest individual

  // Any rank based selection
  switch (selection_scheme_) {
    case RANK_LINEAR :
    case RANK_EXPONENTIAL :
    case FITTEST : {
      assert(prob_reprod_);
      // First we sort the local fitness values using bubble sort :
      // we sort by increasing order, so the first element will have the worst fitness.
      bool swaped = true;
      int16_t loop_length = neighborhood_size-1;
      double  tmp_holder;
      int16_t tmp_holder2;
      while (swaped == true) {
        swaped = false;
        for (int16_t i = 0 ; i < loop_length ; i++) {
          //if the first is higher than the second,  exchange them
          if (sort_fit_array[i] > sort_fit_array[i+1]) {
            tmp_holder = sort_fit_array[i];
            sort_fit_array[i] = sort_fit_array[i+1];
            sort_fit_array[i+1] = tmp_holder;

            tmp_holder2 = initial_location[i];
            initial_location[i] = initial_location[i+1];
            initial_location[i+1] = tmp_holder2;

            swaped = true;
          }
        }

        loop_length = loop_length - 1;
      }


      // Then we use the already computed probabilities
      for (int16_t i = 0 ; i < neighborhood_size ; i++) {
        probs[initial_location[i]] = prob_reprod_[i];
      }

      break;
    }
    // Fitness proportionate selection
    case FITNESS_PROPORTIONATE : {
      for(int16_t i = 0 ; i < neighborhood_size ; i++) {
        probs[i] = local_fit_array[i]/sum_local_fit;
      }

      break;
    }
    default : {
      printf("ERROR, invalid selection scheme in file %s:%d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }

  // pick one organism to reproduce, based on probs[] calculated above, using roulette selection
    /*bool verbose = false;*/
  int16_t found_org = world->grid(x,y)->reprod_prng_->roulette_random(probs, neighborhood_size);

  int16_t x_offset = (found_org / selection_scope_x_) - 1;
  int16_t y_offset = (found_org % selection_scope_x_) - 1;


  delete [] local_fit_array;
  delete [] sort_fit_array;
  delete [] initial_location;
  delete [] probs;


  /*world->grid(x,y)->probs = probs;
    world->grid(x,y)->local_fit_array = local_fit_array;
    world->grid(x,y)->sum_local_fit = sum_local_fit;
    world->grid(x,y)->local_meta_array = local_meta_array;
    world->grid(x,y)->loc_phenotype = loc_phenotype;*/

    exp_m_->next_generation_reproducer_[x*grid_height+y] = ((x+x_offset+grid_width)  % grid_width)*grid_height+
                                            ((y+y_offset+grid_height) % grid_height);

  return world->indiv_at((x+x_offset+grid_width)  % grid_width,
                             (y+y_offset+grid_height) % grid_height);
}

// =================================================================
//                          Non inline accessors
// =================================================================
} // namespace aevol
