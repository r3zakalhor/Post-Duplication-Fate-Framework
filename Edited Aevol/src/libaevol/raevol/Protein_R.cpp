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
//                              Libraries
// =================================================================

// =================================================================
//                            Project Files
// =================================================================
#include "Protein_R.h"
#include "Codon.h"
#include "GeneticUnit.h"
#include "ExpManager.h"

#ifndef __OPENMP_GPU
#include <algorithm>
#endif

namespace aevol {

//##############################################################################
//                                                                             #
//                           Class Protein_R                                #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================
long Protein_R::id = 0;

// =================================================================
//                             Constructors
// =================================================================
Protein_R::Protein_R( GeneticUnit* gen_unit, const Protein_R &model ) : Protein::Protein( gen_unit, model ) {
  concentration_ = model.concentration_;
  _initial_concentration = model.concentration_;
  _delta_concentration = model._delta_concentration;
  _signal = model._signal;
  _inherited = model._inherited;
  is_TF_ = model.is_TF_;
  _id = id++;
  _local_id = model._local_id;

  if (!AA_list_.empty()) {
    _cod_tab = new int8_t[AA_list_.size()];
    int i = 0;
    for (auto cod : AA_list_) {
      _cod_tab[i] = cod->value();
      i++;
    }
  }
}

// Constructor for heritable protein
Protein_R::Protein_R( GeneticUnit* gen_unit, const Protein_R &model, ExpManager* exp_m ) :
    Protein::Protein( gen_unit, model, exp_m ) {
  concentration_ = model.concentration_;
  _initial_concentration = model.concentration_;
  _delta_concentration = model._delta_concentration;
  _signal = model._signal;
  _inherited = model._inherited;
  is_TF_ = model.is_TF_;
  _id = id++;
  _local_id = model._local_id;


  if (exp_m->exp_s()->get_with_heredity()) {
    for (Codon* cod : model.AA_list_) {
      AA_list_.push_back(new Codon(*cod));
    }
    exp_m_ = exp_m;
  }

  if (!AA_list_.empty()) {
    _cod_tab = new int8_t[AA_list_.size()];
    int i = 0;
    for (auto cod : AA_list_) {
      _cod_tab[i] = cod->value();
      i++;
    }
  }
}


Protein_R::Protein_R( GeneticUnit* gen_unit, const std::list<Codon*> codon_list,
							Strand strand, int32_t shine_dal_pos,
                            Rna* rna, double w_max )  :
		Protein::Protein( gen_unit, codon_list, strand, shine_dal_pos, rna, w_max )
{
  _rna_R_list.push_back((Rna_R*)rna);

	_initial_concentration = concentration_;
  _delta_concentration   = 0;
  _inherited             = false;
  _signal                = false;
  is_TF_			       = false;
  _id = id++;

  if (!AA_list_.empty()) {
    _cod_tab = new int8_t[AA_list_.size()];
    int i = 0;
    for (auto cod : AA_list_) {
      _cod_tab[i] = cod->value();
      i++;
    }
  }
}

//used to build the signal protein
Protein_R::Protein_R( const std::list<Codon*> codon_list, ProteinConcentration concentration, double w_max)  :
		Protein::Protein( codon_list, concentration, w_max )
{
  _initial_concentration = 0;
  _delta_concentration  = 0;
  _inherited            = false;
  _signal               = true;

  is_TF_			 = false;
  _id = id++;

  if (!AA_list_.empty()) {
    _cod_tab = new int8_t[AA_list_.size()];
    int i = 0;
    for (auto cod : AA_list_) {
      _cod_tab[i] = cod->value();
      i++;
    }
  }
}

//used to clone the signal protein
Protein_R::Protein_R( Protein_R* signal )  :
    Protein::Protein( signal)
{
  _initial_concentration = 0;
  _delta_concentration  = 0;
  _inherited            = false;
  _signal               = true;
  _local_id = signal->_local_id;

  is_TF_			 = false;
  _id = signal->_id;

  if (!AA_list_.empty()) {
    _cod_tab = new int8_t[AA_list_.size()];
    int i = 0;
    for (auto cod : AA_list_) {
      _cod_tab[i] = cod->value();
      i++;
    }
  }
}

Protein_R::Protein_R( gzFile backup_file ) : Protein::Protein( backup_file )
{
  // the Influence list is re-calculate afterward, and then is not saved, nor use in this consctructor.
  double delta;
  gzread( backup_file, &delta,   	sizeof(delta) );
  _delta_concentration = (ProteinConcentration) delta;
  _initial_concentration = concentration_;

  gzread( backup_file, &_inherited,   			sizeof(_inherited) );
  gzread( backup_file, &_signal,   			sizeof(_signal) );

  gzread( backup_file, &_local_id,   			sizeof(_local_id) );

  is_TF_			 = false;
  _id = id++;

  if (!AA_list_.empty()) {
    _cod_tab = new int8_t[AA_list_.size()];
    int i = 0;
    for (auto cod : AA_list_) {
      _cod_tab[i] = cod->value();
      i++;
    }
  }
}

// =================================================================
//                             Destructors
// =================================================================
Protein_R::~Protein_R( void )
{
  for (unsigned int i = 0; i < _rna_R_list.size(); i++)
    _rna_R_list[i] = NULL;

	_rna_R_list.clear();

  delete [] _cod_tab;
}

// =================================================================
//                            Public Methods
// =================================================================
void Protein_R::compute_delta_concentration( ExpManager* exp_m )
{
  _delta_concentration = 0;

  if( _signal == false )
  {
    //printf("Protein %ld is generated by ",_id);
    int loop_size = _rna_R_list.size();
    for (int i = 0; i < loop_size; i++) {
	//for (auto& rna: _rna_R_list)
    //{
      //if (_id == 34483) printf("%ld (influenced by %ld) at %f  - ",rna->get_id(),rna->_operating_coef_list.size(),rna->get_synthesis_rate());
      //assert( _inherited == false);

      /*if (gen_unit_->indiv()->id() == 12885)
        printf("12608 RNA %d synthesis of %d is %f %f\n",rna->get_id(),get_id(),
               _delta_concentration,rna->get_synthesis_rate());*/

      double synthesis_rate = _rna_R_list[i]->get_synthesis_rate();
      // if (gen_unit_->indiv()->id() == 543 && AeTime::time() == 5895) printf("CPU -- Protein %d synthesis by RNA %d at rate %lf : DELTA BEFORE %f :: ",shine_dal_pos_,_rna_R_list[i]->first_transcribed_pos(),
      //                      synthesis_rate,_delta_concentration);

      _delta_concentration += synthesis_rate;

      // if (gen_unit_->indiv()->id() == 543 && AeTime::time() == 5895) printf("DELTA AFTER %lf : %lf\n",_delta_concentration,synthesis_rate);

//      if (indiv()->id() == 389) printf("UPDATE_NETWORK_SYN_UPDATE Protein CPU %d :: %lf DELTA %lf - %lf -- %d\n",
//                                                                   first_translated_pos(),concentration(),_delta_concentration,_rna_R_list[i]->get_synthesis_rate(),
//                                                                       _rna_R_list[i]->first_transcribed_pos()
//        );
    }
    //if (_id == 34483) printf("\n");
    //if (_id == 34483)  printf("Prot %ld BEFORE DEGRADATION concentration %f %f\n",_id,concentration_,_delta_concentration);
    /*if (gen_unit_->indiv()->id() == 12608)
      printf("12608 RNA synthesis of %d is %f %f %f\n",get_id(),
             _delta_concentration,concentration_,_initial_concentration);*/

    _delta_concentration -= exp_m->exp_s()->get_degradation_rate() * concentration_;
    _delta_concentration *= 1/((ProteinConcentration)exp_m->exp_s()->get_nb_degradation_step());

    //if (_id == 34483)  printf("Prot %ld AFTER degradation concentration %f %f\n",_id,concentration_,_delta_concentration);
  }
}

int8_t Protein_R::get_codon( int32_t index )
{
  return _cod_tab[index];
}

void Protein_R::save( gzFile backup_file )
{
  Protein::save( backup_file );

  // the Influence list is re-calculate afterward, and then is not saved.
  double delta = (double) _delta_concentration;
  gzwrite( backup_file, &delta,   	sizeof(delta) );

  gzwrite( backup_file, &_inherited,   			sizeof(_inherited) );
  gzwrite( backup_file, &_signal,   			sizeof(_signal) );

  gzwrite( backup_file, &_local_id,   			sizeof(_local_id) );
}
// =================================================================
//                           Protected Methods
// =================================================================
void Protein_R::remove_influences( void )
{
  printf("ALERTE la proteine veut détruire une influence !!!\n");

  _rna_R_list.clear();
}


void Protein_R::add_RNA( Rna * rna )
{
  Protein::add_RNA(rna);
  _initial_concentration += rna->basal_level();
  _rna_R_list.push_back((Rna_R*)rna);
}

} // namespace aevol
