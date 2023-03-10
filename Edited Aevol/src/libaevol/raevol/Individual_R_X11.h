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


 #ifndef AEVOL_INDIVIDUAL_R_X11_H_
#define  AEVOL_INDIVIDUAL_R_X11_H_


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>

// =================================================================
//                            Project Files
// =================================================================
#include "Individual_R.h"
#include "Individual_X11.h"

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================
class Individual_R_X11 : public Individual_R, Individual_X11
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
	Individual_R_X11( const Individual_R_X11 &model  );
  Individual_R_X11(ExpManager * exp_m,
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
	Individual_R_X11(  Individual_R_X11* parent, int32_t id,
                     std::shared_ptr<JumpingMT> mut_prng,
                     std::shared_ptr<JumpingMT> stoch_prng);
	Individual_R_X11( ExpManager* exp_m, gzFile backup_file );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~Individual_R_X11( void ) noexcept;
    // =================================================================
    //                              Accessors
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    virtual void display_regulation( X11Window* win );
    virtual void display_concentrations( X11Window* win );
    virtual void display_phenotype( X11Window* win, const Habitat_R& habitat );
    // =================================================================
    //                           Public Attributes
    // =================================================================

  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    /*    ae_individual_R_X11(const ae_individual &model)
    {
      printf("ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      };*/

    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================

};

// =====================================================================
//                          Accessors definitions
// =====================================================================
} // namespace aevol
#endif // AEVOL_INDIVIDUAL_R_X11_H_
