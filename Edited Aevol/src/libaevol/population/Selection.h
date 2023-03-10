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


#ifndef AEVOL_SELECTION_H_
#define AEVOL_SELECTION_H_


// =================================================================
//                              Libraries
// =================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <memory>



// =================================================================
//                            Project Files
// =================================================================
#include "World.h"
#include "Observable.h"

#ifdef __REGUL
#include "raevol/Individual_R.h"
#endif

namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================
class ExpManager;






class Selection : public Observable
{
  public :
    // =================================================================
    //                             Constructors
    // =================================================================
    Selection() = delete;
    Selection(const Selection&) = delete;
    Selection(ExpManager* exp_m);

    // =================================================================
    //                             Destructors
    // =================================================================

    virtual ~Selection();

    // =================================================================
    //                        Accessors: getters
    // =================================================================
    inline SelectionScheme selection_scheme() const;
    inline double               selection_pressure() const;
    inline double*              prob_reprod() const;

    inline SelectionScope selection_scope() const;
    inline int32_t               selection_scope_x() const;
    inline int32_t               selection_scope_y() const;
    // inline std::unique_ptr<JumpingMT> prng() const;

    inline FitnessFunction       fitness_func() const;
    inline int32_t               fitness_function_scope_x() const;
    inline int32_t               fitness_function_scope_y() const;

    // =================================================================
    //                        Accessors: setters
    // =================================================================
    // ----------------------------------------- Pseudo-random number generator
    //inline void set_prng(std::unique_ptr<JumpingMT>&& prng);

    // -------------------------------------------------------------- Selection
    inline void set_selection_scheme(SelectionScheme sel_scheme);
    inline void set_selection_pressure(double sel_pressure);

    inline void set_selection_scope(SelectionScope sel_scope);
    inline void set_selection_scope_x(int32_t sel_scope_x);
    inline void set_selection_scope_y(int32_t sel_scope_y);

    inline void set_fitness_function(FitnessFunction fit_func) {fitness_function_ = fit_func; };
    inline void set_fitness_function_scope_x(int32_t fit_scope_x) { fitness_function_scope_x_ = fit_scope_x; };
    inline void set_fitness_function_scope_y(int32_t fit_scope_y) { fitness_function_scope_y_ = fit_scope_y; };

    // =================================================================
    //                              Operators
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    void step_to_next_generation();
    void PerformPlasmidTransfers();
    void write_setup_file(gzFile setup_file) const;
    void save(gzFile& backup_file) const;
    void load(gzFile& exp_setup_file, gzFile& backup_file, bool verbose);
#ifndef __REGUL
    void run_life(Individual* new_indiv);
#else
    void run_life(Individual_R* new_indiv);
#endif

    Individual* do_replication(Individual* parent,
                               unsigned long long index,
                               int8_t &type_mutate,
                               int16_t x = -1,
                               int16_t y = -1);

    // =================================================================
    //                           Public Attributes
    // =================================================================


    void set_unique_id(unsigned long long uid) { unique_id = uid; }

    
    void compute_prob_reprod();
    void compute_local_prob_reprod();
    Individual* do_local_competition(int16_t x, int16_t y);
    double* prob_reprod_;


 protected :
    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =======================================================================
    //                             Protected Attributes
    // =======================================================================
    ExpManager* exp_m_;

    // ----------------------------------------- Pseudo-random number generator


    // -------------------------------------------------------------- Selection
    SelectionScheme selection_scheme_;
    double selection_pressure_;

    SelectionScope selection_scope_;
    int32_t selection_scope_x_;
    int32_t selection_scope_y_;

    FitnessFunction fitness_function_;
    int32_t fitness_function_scope_x_;
    int32_t fitness_function_scope_y_;

    double* fitness_sum_tab_;
    // --------------------------- Probability of reproduction of each organism

  Individual *** reproducers;


    unsigned long long unique_id = 16000;

    long apply_mutation[1024];

    std::vector<int32_t> to_evaluate;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================
// inline std::unique_ptr<JumpingMT> Selection::prng() const
// {
//   return prng_;
// }

inline SelectionScheme Selection::selection_scheme() const
{
  return selection_scheme_;
}

inline double Selection::selection_pressure() const
{
  return selection_pressure_;
}

inline SelectionScope Selection::selection_scope() const
{
  return selection_scope_;
}

inline int32_t Selection::selection_scope_x() const
{
  return selection_scope_x_;
}

inline int32_t Selection::selection_scope_y() const
{
  return selection_scope_y_;
}

    inline FitnessFunction Selection::fitness_func() const
    {
        return fitness_function_;
    }

    inline int32_t Selection::fitness_function_scope_x() const
    {
        return fitness_function_scope_x_;
    }

    inline int32_t Selection::fitness_function_scope_y() const
    {
        return fitness_function_scope_y_;
    }

inline double*Selection::prob_reprod() const
{
  if (prob_reprod_ == NULL)
  {
    printf("ERROR, prob_reprod_ has not been computed %s:%d\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  return prob_reprod_;
}

// =====================================================================
//                           Setters' definitions
// =====================================================================
// ----------------------------------------- Pseudo-random number generator
/*inline void Selection::set_prng(std::unique_ptr<JumpingMT>&& prng)
{
  prng_ = std::move(prng);
}*/

// -------------------------------------------------------------- Selection
inline void Selection::set_selection_scheme(SelectionScheme sel_scheme)
{
  selection_scheme_ = sel_scheme;
}

inline void Selection::set_selection_pressure(double sel_pressure)
{
  selection_pressure_ = sel_pressure;
}

inline void Selection::set_selection_scope(SelectionScope sel_scope)
{
  selection_scope_ = sel_scope;
}

inline void Selection::set_selection_scope_x(int32_t sel_scope_x)
{
  selection_scope_x_ = sel_scope_x;
}

inline void Selection::set_selection_scope_y(int32_t sel_scope_y)
{
  selection_scope_y_ = sel_scope_y;
}
// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol
#endif // AEVOL_SELECTION_H_
