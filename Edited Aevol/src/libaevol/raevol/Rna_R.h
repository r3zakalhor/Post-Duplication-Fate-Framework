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


#ifndef  AEVOL_RNA_R_H_
#define  AEVOL_RNA_R_H_


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <unordered_map>
// =================================================================
//                            Project Files
// =================================================================
#include "Rna.h"
#include "Protein_R.h"
#include "macros.h"

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================

class Protein_R;
class Individual_R;

class Rna_R : public Rna
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    Rna_R( GeneticUnit* gen_unit, const Rna_R &model );
	  Rna_R( GeneticUnit* gen_unit, Strand strand, int32_t index, int8_t diff );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~Rna_R( void );

    // =================================================================
    //                              Accessors
    // =================================================================
    inline std::vector<Protein_R*> get_protein_list( void );
    void set_local_id(int local_id) { _local_id = local_id; }

    int get_local_id() { return _local_id; }
    // =================================================================
    //                            Public Methods
    // =================================================================
    void    set_influences( std::list<Protein*>& protein_list, int id );
    ProteinConcentration  get_synthesis_rate( void );
    ProteinConcentration  get_affinity_with_protein( int32_t index, Protein *protein );
    int32_t get_enhancer_position( void );
    int32_t get_operator_position( void );

    long get_id() { return _id; };

    int nb_influences() { return _nb_influences; }

#ifdef __PROXY_POW_APPROX
    static void load_lookup_table();
    static void create_lookup_table(double hill_shape_theta, double hill_shape_n);
#endif
    // =================================================================
    //                           Public Attributes
    // =================================================================
    std::vector<Protein_R*> _protein_list;
    std::vector<ProteinConcentration> _enhancing_coef_list;
    std::vector<ProteinConcentration> _operating_coef_list;
/*
 *     std::unordered_map<int,Protein_R*> _protein_list;
    std::unordered_map<int,double> _enhancing_coef_list;
    std::unordered_map<int,double> _operating_coef_list;
 */
    //std::vector<double> _protein_concentration_list;
    static long id;
    static int8_t quadon_map[4];

#ifdef __PROXY_POW_LOOKUP
    static std::map<double,double> pow_cache;
#elif __PROXY_POW_APPROX
#ifdef __REGUL
    static double lookup_table_pow[LOOKUP_TABLE_SIZE];
#endif
#endif


    int _nb_influences = 0;
  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    /*    ae_rna()
    {
      printf("ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      };*/

    // =================================================================
    //                           Protected Methods
    // =================================================================
    //inline  ae_rna_R* copy( void );
    ProteinConcentration    affinity_with_protein( int32_t index, Protein *protein );

    // =================================================================
    //                          Protected Attributes
    // ================================================================

    long _id;
    int       _local_id;


    Individual_R* indiv_;

    double hill_shape_n_;
    double hill_shape_;

};




// =====================================================================
//                          Accessors definitions
// =====================================================================
// =====================================================================
//                       Inline functions' definition
// =====================================================================


} // namespace aevol

#endif // AEVOL_RNA_R_H_
