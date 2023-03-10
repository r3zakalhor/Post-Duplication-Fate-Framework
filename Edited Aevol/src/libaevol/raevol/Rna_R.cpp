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
//*****************************************************************************




// =================================================================
//                              Libraries
// =================================================================
//#include <math.h>
#ifdef __BLAS__
#include <cblas.h>
#endif
// =================================================================
//                            Project Files
// =================================================================
#include "Rna_R.h"
#include "ExpManager.h"

namespace aevol {
//##############################################################################
//                                                                             #
//                           Class ae_rna_R                                    #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

long Rna_R::id = 0;
int8_t Rna_R::quadon_map[4] = {8, 4, 2, 1};

#ifdef __PROXY_POW_LOOKUP
std::map<double,double> Rna_R::pow_cache;
#elif __PROXY_POW_APPROX
double Rna_R::lookup_table_pow[LOOKUP_TABLE_SIZE];
#endif
// =================================================================
//                             Constructors
// =================================================================
Rna_R::Rna_R( GeneticUnit* gen_unit, const Rna_R &model ) : Rna( gen_unit, model )
{
  _protein_list = model._protein_list;
  _enhancing_coef_list = model._enhancing_coef_list;
  _operating_coef_list = model._operating_coef_list;
  _id = model._id;
  _nb_influences = model._nb_influences;

  hill_shape_n_= gen_unit_->exp_m()->exp_s()->get_hill_shape_n();
  hill_shape_ = gen_unit_->exp_m()->exp_s()->get_hill_shape();
}

Rna_R::Rna_R( GeneticUnit* gen_unit, Strand strand, int32_t index, int8_t diff ) :
		Rna( gen_unit, strand, index, diff )
{
  _id = id++;
  _nb_influences = 0;

  hill_shape_n_= gen_unit_->exp_m()->exp_s()->get_hill_shape_n();
  hill_shape_ = gen_unit_->exp_m()->exp_s()->get_hill_shape();
}

/*
ae_rna_R::ae_rna_R( ae_rna_R* parent ) :
ae_rna( parent )
{
  _influence_list = parent->_influence_list->copy();
}
*/

// =================================================================
//                             Destructors
// =================================================================
Rna_R::~Rna_R( void )
{
	_protein_list.clear();
	_enhancing_coef_list.clear();
	_operating_coef_list.clear();
}

// =================================================================
//                            Public Methods
// =================================================================
void Rna_R::set_influences( std::list<Protein*>& protein_list, int id )
{
	int32_t enhancer_position = get_enhancer_position();
	int32_t operator_position = get_operator_position();

	_protein_list.clear();
  std::vector<Protein_R*>::iterator itprot = _protein_list.begin();
  _protein_list.resize(0);
 _enhancing_coef_list.clear();
  std::vector<ProteinConcentration>::iterator itenh = _enhancing_coef_list.begin();
	_enhancing_coef_list.resize(0);
  _operating_coef_list.clear();
  std::vector<ProteinConcentration>::iterator itope = _operating_coef_list.begin();
	_operating_coef_list.resize(0);


  int i = 0;
  ProteinConcentration enhance=0,operate=0;

  //#pragma omp simd
	for (auto& prot : protein_list) {
    enhance = affinity_with_protein( enhancer_position, prot );
    operate = affinity_with_protein( operator_position, prot );

    if (enhance != 0.0 || operate != 0.0) {
//      if (id==389) printf("Add Affinity for RNA %d with Protein %d : E %lf O %lf\n",
//             first_transcribed_pos(),i,enhance,operate);

      //  if (id==68 && AeTime::time() == 4) 
      // if (id==0)
      //   printf("CPU -- Affinity between RNA %d and Protein %d : %lf %lf\n",
      //          first_transcribed_pos(),prot->shine_dal_pos(),enhance,operate);

      _protein_list.insert(itprot,(Protein_R*) prot);

      _enhancing_coef_list.insert(itenh, enhance);
      _operating_coef_list.insert(itope, operate);

  //    if (std::isnan(_protein_list[i]->concentration_)) {
  //      printf("Wierd concentration?\n");
  //    }

      _protein_list[i]->is_TF_ = true;

     // _protein_concentration_list[i] = prot->concentration();
      i++;
      itprot = _protein_list.begin()+i;
      itenh = _enhancing_coef_list.begin()+i;
      itope = _operating_coef_list.begin()+i;
    }
    //else
    // _protein_list[i] = nullptr;
	}
  _nb_influences = i;

  //TODO NOT USEFUL ??? _nb_influences = i==0 ? 0 : i-1;

  /*if (gen_unit_->indiv()->id() == 12608)
    printf("12608 RNA %d is influenced by %d proteins\n",_id,_nb_influences);*/

  /*if (protein_list.size() > 0)
    printf("Set Influences of RNA %ld with %ld %ld %ld\n",_id,_enhancing_coef_list.size(),_operating_coef_list.size(),
         _protein_list.size());*/

}

ProteinConcentration Rna_R::get_synthesis_rate( void )
{
  ProteinConcentration enhancer_activity  = 0;
  ProteinConcentration operator_activity  = 0;

//#ifndef __BLAS__
  for (int i = 0; i < _nb_influences; i++) {
    /*if (gen_unit_->indiv()->id() == 12608)
      printf("12608 RNA %d  due to Protein %d is %f %f %f\n",
             _id,_protein_list[i]->get_id(),
             _enhancing_coef_list[i],_operating_coef_list[i],_protein_list[i]->concentration_);*/
//    printf("E[%d] %f %f %f\n",i,enhancer_activity,_enhancing_coef_list[i],_protein_list[i]->concentration_);

    enhancer_activity  += _enhancing_coef_list[i] * _protein_list[i]->concentration_;

//    printf("O[%d] %f %f %f\n",i,operator_activity,_operating_coef_list[i],_protein_list[i]->concentration_);
    operator_activity  += _operating_coef_list[i] * _protein_list[i]->concentration_;
  //  if (gen_unit_->indiv()->id()==70 && AeTime::time() == 1595)
  //    printf("CPU -- RNA %d Protein %d (%lf) :: Enhancer %lf Operator %lf\n",first_transcribed_pos(),_protein_list[i]->first_translated_pos(),
  //           _protein_list[i]->concentration_,  _enhancing_coef_list[i], _operating_coef_list[i]);
  }


//  if (gen_unit_->indiv()->id()==137)
//    printf("CPU -- RNA %d Enhancer %lf Operator %lf\n",first_transcribed_pos(),enhancer_activity,operator_activity);
/*#else
  ProteinConcentration enhancer_tab[_nb_influences];
  ProteinConcentration operator_tab[_nb_influences];

#ifdef __SIMD
  #pragma omp simd
#endif
  for (int i = 0; i < _nb_influences; i++) {
  	enhancer_tab[i] = _enhancing_coef_list[i] * _protein_list[i]->concentration_;
  	}

#ifdef __SIMD
  #pragma omp simd
#endif
  for (int i = 0; i < _nb_influences; i++) {
    operator_tab[i] = _operating_coef_list[i] * _protein_list[i]->concentration_;
  }

#ifdef __FLOAT_CONCENTRATION
  ProteinConcentration enhancer_activity_blas = cblas_sasum(_nb_influences,enhancer_tab,1);
  ProteinConcentration operator_activity_blas = cblas_sasum(_nb_influences,operator_tab,1);
#else
  ProteinConcentration enhancer_activity_blas = cblas_dasum(_nb_influences,enhancer_tab,1);
  ProteinConcentration operator_activity_blas = cblas_dasum(_nb_influences,operator_tab,1);
#endif
#endif
*/

#ifndef __PROXY_POW


  ProteinConcentration enhancer_activity_pow_n = enhancer_activity == 0 ? 0 :
                                                 pow( enhancer_activity, hill_shape_n_ );
  ProteinConcentration operator_activity_pow_n = operator_activity == 0 ? 0 :
                                                 pow( operator_activity, hill_shape_n_ );

#elif __PROXY_POW_LOOKUP

  ProteinConcentration enhancer_activity_pow_n;
  ProteinConcentration operator_activity_pow_n;

  if (enhancer_activity==0) {
    enhancer_activity_pow_n=0;
  } else {
#pragma omp critical(readpowcache)
  {
    try {
      enhancer_activity_pow_n = pow_cache[enhancer_activity];
    } catch(std::out_of_range& e) {
      pow_cache[enhancer_activity] = pow( enhancer_activity, hill_shape_n_ );
      enhancer_activity_pow_n = pow_cache[enhancer_activity];
    }
  }
  }

  if (operator_activity==0) {
    operator_activity_pow_n=0;
  } else {
#pragma omp critical(readpowcache)
  {
    try {
      operator_activity_pow_n = pow_cache[operator_activity];
    } catch(std::out_of_range& e) {
      pow_cache[enhancer_activity] = pow( operator_activity, hill_shape_n_ );
      operator_activity_pow_n = pow_cache[operator_activity];
    }
  }
  }
#elif __PROXY_POW_APPROX
//  printf("E %f %d %d\n",enhancer_activity,(int)(enhancer_activity*LOOKUP_TABLE_SIZE),
//  ((int)(enhancer_activity*LOOKUP_TABLE_SIZE))+1);
//
//
//printf("O %f %d %d\n",operator_activity,(int)(operator_activity*LOOKUP_TABLE_SIZE),
//  ((int)(operator_activity*LOOKUP_TABLE_SIZE))+1);


  ProteinConcentration enhancer_activity_pow_n = enhancer_activity == 0 ? 0 :
      enhancer_activity <= 1 ? ( lookup_table_pow[(int)(enhancer_activity*LOOKUP_TABLE_SIZE)] +
      lookup_table_pow[((int)(enhancer_activity*LOOKUP_TABLE_SIZE))+1] )/2 :
      pow( enhancer_activity, hill_shape_n_ );
  ProteinConcentration operator_activity_pow_n = operator_activity == 0 ? 0 :
      operator_activity <= 1 ? ( lookup_table_pow[(int)(operator_activity*LOOKUP_TABLE_SIZE)] +
      lookup_table_pow[((int)(operator_activity*LOOKUP_TABLE_SIZE))+1] )/2 :
      pow( operator_activity, hill_shape_n_ );
#endif
  //if (enhancer_activity != 0.0 || operator_activity != 0.0)
  /*if (_id == 132073) printf("Synthesis of RNA %ld : E %f O %f EP %f OP %f SN %f S %f B %f\n",_id,enhancer_activity,operator_activity,enhancer_activity_pow_n,
                                                operator_activity_pow_n,gen_unit_->exp_m()->exp_s()->get_hill_shape_n(),
         gen_unit_->exp_m()->exp_s()->get_hill_shape(),basal_level_);*/
  /*if (gen_unit_->indiv()->id() == 12608)
    printf("12608 RNA %d is %f %f %f\n",_id,enhancer_activity,operator_activity,basal_level_);*/


  ProteinConcentration ret = basal_level_
           * (hill_shape_
              / (operator_activity_pow_n + hill_shape_))
           * (1 + ((1 / basal_level_) - 1)
                  * (enhancer_activity_pow_n /
                     (enhancer_activity_pow_n + hill_shape_)));

//  if (std::isnan(ret)) {
//    printf("syn is nan B %f O %f ON %f E %f EN %f\n",basal_level_,
//          operator_activity,operator_activity_pow_n,enhancer_activity,enhancer_activity_pow_n);
//    exit(88);
//  }
  return ret;
}

// =================================================================
//                           Protected Methods
// =================================================================
int32_t Rna_R::get_enhancer_position( void )
{
  int32_t length = gen_unit_->dna()->length();

  if(strand_ == LEADING)
  {
    return (pos_ - 20)  % ( length ) < 0 ?
           ((pos_ - 20)  % ( length )) + ( length ) :
           (pos_ - 20)  % ( length );
  }
  else  // strand_ = LAGGING
  {
    return (pos_ + 20)  % ( length ) < 0 ?
           ((pos_ + 20)  % ( length )) + ( length ) :
           (pos_ + 20)  % ( length );
  }
}

int32_t Rna_R::get_operator_position( void )
{
  int32_t length = gen_unit_->dna()->length();

  if(strand_ == LEADING)
  {
    return (pos_ + PROM_SIZE)  % ( length ) < 0 ?
           (pos_ + PROM_SIZE)  % ( length ) + (length) :
           (pos_ + PROM_SIZE)  % ( length );
  }
  else  // strand_ = LAGGING
  {
    return (pos_ - PROM_SIZE)  % ( length ) < 0 ?
           (pos_ - PROM_SIZE)  % ( length ) + (length) :
           (pos_ - PROM_SIZE)  % ( length );
  }
}

ProteinConcentration Rna_R::affinity_with_protein( int32_t index, Protein *protein )
{
  int32_t len = protein->length();

  if (len > 5) {

    ProteinConcentration max = 0;
    ProteinConcentration temp = 1;

    int32_t quadon_tab[5];
    Protein_R* prot = NULL;
    prot = (Protein_R*) (protein);

    const char* dna = gen_unit_->dna()->data();
    int32_t  len_dna    = gen_unit_->dna()->length();

    for (int32_t pos = index; pos < index+5; pos++) {

      int8_t quadon[4];

      if ( strand_ == LEADING )
      {
        for (int8_t i = 0; i < QUADON_SIZE; i++) {
          quadon[i] = (dna[Utils::mod(pos + i,len_dna)] == '1')
                      ? 1 << (QUADON_SIZE - i - 1)
                      : 0;
        }
      } else {
        for (int8_t i = 0; i < QUADON_SIZE; i++) {
          quadon[i] = (dna[Utils::mod(pos - i,len_dna)] != '1')
                      ? (QUADON_SIZE - i - 1)
                      : 0;
        }
      }/*
        if (pos+QUADON_SIZE < len_dna) {
          for (int8_t i = 0; i < QUADON_SIZE; i++) {
            quadon[i] = (dna[pos + i] == '1') ? 1 <<  (QUADON_SIZE - i - 1) : 0;
          }
        } else {
          for (int8_t i = 0; i < QUADON_SIZE; i++) {
            quadon[i] = (dna[((pos + i) % len_dna < 0 ? (pos + i) % len_dna +
                                                        len_dna : (pos + i) %
                                                                  len_dna)] ==
                         '1') ?  (QUADON_SIZE - i - 1) : 0;
          }
        }
      }
      else  // ( strand == LAGGING )
      {
        if (pos-QUADON_SIZE >= 0) {
          for (int8_t i = 0; i < QUADON_SIZE; i++) {
            quadon[i] = (dna[pos - i] != '1') ?  (QUADON_SIZE - i - 1) : 0;
          }
        } else {
          for (int8_t i = 0; i < QUADON_SIZE; i++) {
            quadon[i] = (dna[((pos - i) % len_dna < 0 ? (pos - i) % len_dna +
                                                        len_dna : (pos - i) %
                                                                  len_dna)] !=
                         '1') ?  (QUADON_SIZE - i - 1) : 0;
          }
        }
      }*/

      quadon_tab[pos-index] = quadon[0]+quadon[1]+quadon[2]+quadon[3];//gen_unit_->indiv_r_->get_quadon(gen_unit_, strand_, (index + i));
    }

    for (int32_t i = 0; i < len - 4; i++) {
      temp = 1;

      for (int8_t j = 0; j < 5; j++) {
        // if (gen_unit_->indiv()->id()==660) {
        //   printf("Individual %d Protein %d\n",gen_unit_->indiv()->id(),prot->shine_dal_pos());
        //   printf("Codon[%d] (i %d j %d) %d out of %d\n",i+j,i,j,prot->_cod_tab[i+j],MAX_CODON);
        //   printf("Protein Length %d\n",prot->length());
        // }

        temp *= gen_unit_->exp_m()->exp_s()->_binding_matrix[quadon_tab[j]][prot->_cod_tab[
            i + j]];
      }

      max = (max < temp) ? temp : max;

    }

    return max;
  } else {
    return 0.0;
  }
}

#ifdef __PROXY_POW_APPROX
void Rna_R::load_lookup_table() {
  char* lookup_table_file_name = new char[100];

  sprintf( lookup_table_file_name, "lookup_table.ae" );

  gzFile lookup_table_file = gzopen( lookup_table_file_name, "r" );

  if ( lookup_table_file == Z_NULL )
  {
    printf( "ERROR : Could not read lookup table file %s\n", lookup_table_file_name );
    exit( EXIT_FAILURE );
  }

  double value;
  for (int i=0; i < LOOKUP_TABLE_SIZE; i++) {
    gzread( lookup_table_file, &value, sizeof(double));
    lookup_table_pow[i] = (ProteinConcentration) value;
  }

  gzclose( lookup_table_file );

  delete[] lookup_table_file_name;
}

void Rna_R::create_lookup_table(double hill_shape_theta, double hill_shape_n) {
  char* lookup_table_file_name = new char[100];

  sprintf( lookup_table_file_name, "lookup_table.ae" );

  gzFile lookup_table_file = gzopen( lookup_table_file_name, "w" );

  if ( lookup_table_file == Z_NULL )
  {
    printf( "ERROR : Could not write lookup table file %s\n", lookup_table_file_name );
    exit( EXIT_FAILURE );
  }

  double value;
  double hill_n = pow( hill_shape_theta, hill_shape_n );

  for (int i=0; i < LOOKUP_TABLE_SIZE; i++) {
    value = pow(i / LOOKUP_TABLE_SIZE, hill_n);

    gzwrite(lookup_table_file, &value, sizeof(double));
  }

  gzclose( lookup_table_file );

  delete[] lookup_table_file_name;
}
#endif
} // namespace aevol
