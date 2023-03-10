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

#ifndef AE_PROTEIN_R_H
#define AE_PROTEIN_R_H

// =================================================================
//                              Libraries
// =================================================================

// =================================================================
//                            Project Files
// =================================================================
#include "Codon.h"
#include "Protein.h"
#include "Rna_R.h"

namespace aevol {

//##############################################################################
//                                                                             #
//                           Class ae_protein_R                                #
//                                                                             #
//##############################################################################
class GeneticUnit;

class Protein_R : public Protein
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
  Protein_R() = delete;
  Protein_R(const Protein_R &model) = delete;
	Protein_R( GeneticUnit* gen_unit, const Protein_R &model );
  Protein_R( GeneticUnit* gen_unit, const Protein_R &model, ExpManager* exp_m );
	Protein_R( GeneticUnit* gen_unit,
    		const std::list<Codon*> codon_list,
    		Strand strand,
    		int32_t shine_dal_pos,
    		Rna* rna,
        double w_max ); // TODO Rna_R?
	Protein_R( const std::list<Codon*> codon_list, ProteinConcentration concentration, double w_max);
    Protein_R( Protein_R* signal );
    Protein_R( gzFile backup_file );

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~Protein_R( void );

    // =================================================================
    //                              Accessors
    // =================================================================
    inline void     set_inherited( bool is_inherited );
    inline void     set_signal( bool is_signal);
    inline bool     is_inherited( void );
    inline bool     is_signal( void );

    void set_local_id(int local_id) { _local_id = local_id; }

    // =================================================================
    //                            Public Methods
    // =================================================================
    //inline ae_protein_R* copy( void );
    inline void    multiply_concentration( ProteinConcentration factor );
    inline void    set_concentration ( ProteinConcentration concentration);
    inline void    update_concentration( void );
    inline void    reset_concentration( void );
    inline void    set_initial_concentration( void );
           void    compute_delta_concentration( ExpManager* exp_m );
           int8_t  get_codon( int32_t index );
//           void    add_influence( ae_influence_R* influence );
	         void    save( gzFile backup_file );
//	         void    remove_influence( ae_influence_R* influence );
    inline int8_t  get_cod_tab(int32_t index) const;

    void  add_RNA( Rna * rna );

    long get_id() { return _id; };

    int get_local_id() { return _local_id; }

    // =================================================================
    //                           Public Attributes
    // =================================================================
    std::vector<Rna_R*>  _rna_R_list;
	bool is_TF_;

    static long id;

    int8_t*   _cod_tab;
    //bool _concentration_has_change = true;
    ProteinConcentration    _initial_concentration = -1.0; // concentration at cell birth

  ProteinConcentration    _delta_concentration;
  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================

    // =================================================================
    //                           Protected Methods
    // =================================================================
    void remove_influences( void );

    // =================================================================
    //                          Protected Attributes
    // =================================================================

    bool      _inherited;
    bool      _signal;

    long      _id;
    int       _local_id;
    ExpManager* exp_m_;

};

// =====================================================================
//                          Accessors definitions
// =====================================================================
//std::vector<ae_influence_R*> ae_protein_R::get_influence_list( void )
//{
//  return _influence_list;
//}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline void Protein_R::update_concentration( void )
{
 // _concentration_has_change = _delta_concentration != 0 ? true : false;

  concentration_ += _delta_concentration;
}

inline void Protein_R::set_inherited( bool is_inherited )
{
  _inherited = is_inherited;
}

inline void Protein_R::set_signal( bool is_signal )
{
  _signal = is_signal;
}

inline void Protein_R::reset_concentration( void )
{
  concentration_ = _initial_concentration;
}

inline void Protein_R::set_initial_concentration( void )
{
  _initial_concentration = concentration_;
}

inline bool Protein_R::is_inherited( void )
{
  return _inherited;
}

inline bool Protein_R::is_signal( void )
{
  return _signal;
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================

inline void Protein_R::multiply_concentration( ProteinConcentration factor )
{
  concentration_ *= factor;
}

inline void Protein_R::set_concentration( ProteinConcentration concentration )
{
  concentration_ = concentration;
}

inline int8_t Protein_R::get_cod_tab(int32_t index) const
{
  return _cod_tab[index];
}

} // namespace aevol


#endif // AE_PROTEIN_R_H
