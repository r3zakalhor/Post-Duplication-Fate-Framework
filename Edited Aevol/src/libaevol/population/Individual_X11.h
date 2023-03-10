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


#ifndef  AEVOL_INDIVIDUAL_X11_H_
#define  AEVOL_INDIVIDUAL_X11_H_


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>



// =================================================================
//                            Project Files
// =================================================================
#include "Individual.h"
#include "X11Window.h"

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================






class Individual_X11 : public virtual Individual
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    Individual_X11(ExpManager * exp_manager, gzFile backup_file);
    //Individual_X11(const Individual_X11 &model, bool replication_report_copy);
    Individual_X11(const Individual_X11& other);
    Individual_X11(ExpManager * exp_m,
      std::shared_ptr<JumpingMT> mut_prng,
      std::shared_ptr<JumpingMT> stoch_prng,
      std::shared_ptr<MutationParams> param_mut,
      double w_max,
      int32_t min_genome_length,
      int32_t max_genome_length,
      bool allow_plasmids,
      int32_t id,
      const char* strain_name,
      int32_t age);
    Individual_X11(Individual_X11 * const parent, int32_t id,
                      std::shared_ptr<JumpingMT> mut_prng,
                      std::shared_ptr<JumpingMT> stoch_prng);
    Individual_X11() = delete; // forbidden constructor

    // =================================================================
    //                             Destructors
    // =================================================================

    virtual ~Individual_X11();

    // =================================================================
    //                              Accessors
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    virtual void display();
    virtual void display_cdss(X11Window * win);
    virtual void display_rnas(X11Window * win);

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :
    // =================================================================
    //                           Protected Methods
    // =================================================================
    void reset_sectors();
    void add_layer();
    void init_occupied_sectors();

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    // These are used to manage overlapping CDS and RNA display
    int16_t outmost_layer_;
    bool*   occupied_sectors_[2][100];  // TODO : find a way to manage this table's size properly?
};


// =====================================================================
//                          Accessors definitions
// =====================================================================

} // namespace aevol

#endif // AEVOL_INDIVIDUAL_X11_H_
