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
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <set>

// =================================================================
//                            Project Files
// =================================================================
#include "ExpSetup.h"
#include "JumpingMT.h"


//#include "GzHelpers.h"

namespace aevol {
//##############################################################################
//                                                                             #
//                              Class ExpSetup                                 #
//                                                                             #
//##############################################################################

// ===========================================================================
//                         Definition of static attributes
// ===========================================================================

// ===========================================================================
//                                  Constructors
// ===========================================================================
ExpSetup::ExpSetup( ExpManager * exp_m )
{
  exp_m_  = exp_m;

  // -------------------------------------------------------------- Selection
  sel_ = new Selection( exp_m );

  // --------------------------------------------------------------- Transfer
  with_HT_                    = false;
  repl_HT_with_close_points_  = false;
  HT_ins_rate_                = 0.0;
  HT_repl_rate_               = 0.0;
  repl_HT_detach_rate_         = 0.0;

  // --------------------------------------------------------------- Plasmids
  with_plasmids_    = false;
  prob_plasmid_HT_  = 0.0;
  tune_donor_ability_     = 0.0;
  tune_recipient_ability_ = 0.0;
  donor_cost_       = 0.0;
  recipient_cost_   = 0.0;
  swap_GUs_         = false;

  // -------------------------------------------------------------- Secretion
  with_secretion_ = false;
  secretion_contrib_to_fitness_ = 0.0;
  secretion_cost_               = 0.0;

#ifdef BASE_4
  // -------------------------------------------------------------- Terminators
  terminator_polya_sequence_length_ = 0;

  // -------------------------------------------------------------- MWH bases configuration
  for(int8_t i = 0; i < NB_AMINO_ACIDS; i++) {
    aa_base_m_[i] = -1;
    aa_base_w_[i] = -1;
    aa_base_h_[i] = -1;
  }

   // randomly generated no epistasis code
  
  aa_base_m_[ARGININE] = 2;    
  aa_base_m_[ASPARAGINE] = 4;
  aa_base_m_[ASPARTIC_ACID] = 5;
  aa_base_m_[GLUTAMIC_ACID] = 6;
  aa_base_m_[GLUTAMINE] = 0; 
  aa_base_m_[HISTIDINE] = 1;   
  aa_base_m_[SERINE] = 3;       

  aa_base_w_[ALANINE] = 5;  
  aa_base_w_[CYSTEINE] = 0;   
  aa_base_w_[ISOLEUCINE] = 4;
  aa_base_w_[LEUCINE] = 3; 
  aa_base_w_[THREONINE] = 1;      
  aa_base_w_[PHENYLALANINE] = 2;   
  
  aa_base_h_[GLYCINE] = 3;    
  aa_base_h_[METHIONINE] = 2;
  aa_base_h_[LYSINE] = 5;          
  aa_base_h_[PROLINE] = 6;       
  aa_base_h_[TRYPTOPHAN] = 1;    
  aa_base_h_[TYROSINE] = 4;   
  aa_base_h_[VALINE] = 0; 

  /*          
  aa_base_x_[ALANINE] = ;  
  aa_base_x_[ARGININE] = ;    
  aa_base_x_[ASPARAGINE] = ;
  aa_base_x_[ASPARTIC_ACID] = ;
  aa_base_x_[CYSTEINE] = ;
  aa_base_x_[GLUTAMIC_ACID] = ;
  aa_base_x_[GLUTAMINE] = ; 
  aa_base_x_[GLYCINE] = ;
  aa_base_x_[HISTIDINE] = ;
  aa_base_x_[ISOLEUCINE] = ; 
  aa_base_x_[LEUCINE] = ; 
  aa_base_x_[LYSINE] = ;       
  aa_base_x_[METHIONINE] = ;
  aa_base_x_[PHENYLALANINE] = ;    
  aa_base_x_[PROLINE] = ;      
  aa_base_x_[SERINE] = ;     
  aa_base_x_[THREONINE] = ;   
  aa_base_x_[TRYPTOPHAN] = ;    
  aa_base_x_[TYROSINE] = ;   
  aa_base_x_[VALINE] = ;  
 */

  // GB 07/01/22 randomly generated full code

  // aa_base_m_[ALANINE] = 9;  
  // aa_base_m_[ARGININE] = 2;    
  // aa_base_m_[ASPARAGINE] = 4;
  // aa_base_m_[ASPARTIC_ACID] = 5;
  // aa_base_m_[CYSTEINE] = 14;
  // aa_base_m_[GLUTAMIC_ACID] = 6;
  // aa_base_m_[GLUTAMINE] = 0; 
  // aa_base_m_[GLYCINE] = 16;
  // aa_base_m_[HISTIDINE] = 1;
  // aa_base_m_[ISOLEUCINE] = 17; 
  // aa_base_m_[LEUCINE] = 19; 
  // aa_base_m_[LYSINE] = 11;       
  // aa_base_m_[METHIONINE] = 10;
  // aa_base_m_[PHENYLALANINE] = 12;    
  // aa_base_m_[PROLINE] = 15;      
  // aa_base_m_[SERINE] = 3;     
  // aa_base_m_[THREONINE] = 8;   
  // aa_base_m_[TRYPTOPHAN] = 7;    
  // aa_base_m_[TYROSINE] = 13;   
  // aa_base_m_[VALINE] = 18;

  // aa_base_w_[ALANINE] = 5;  
  // aa_base_w_[ARGININE] = 0;    
  // aa_base_w_[ASPARAGINE] = 7;
  // aa_base_w_[ASPARTIC_ACID] = 4;
  // aa_base_w_[CYSTEINE] = 13;
  // aa_base_w_[GLUTAMIC_ACID] = 15;
  // aa_base_w_[GLUTAMINE] = 10; 
  // aa_base_w_[GLYCINE] = 12;
  // aa_base_w_[HISTIDINE] = 19;
  // aa_base_w_[ISOLEUCINE] = 11; 
  // aa_base_w_[LEUCINE] = 3; 
  // aa_base_w_[LYSINE] = 1;       
  // aa_base_w_[METHIONINE] = 18;
  // aa_base_w_[PHENYLALANINE] = 14;    
  // aa_base_w_[PROLINE] = 8;      
  // aa_base_w_[SERINE] = 17;     
  // aa_base_w_[THREONINE] = 6;   
  // aa_base_w_[TRYPTOPHAN] = 2;    
  // aa_base_w_[TYROSINE] = 9;   
  // aa_base_w_[VALINE] = 16;  
  
  // aa_base_h_[ALANINE] = 13;  
  // aa_base_h_[ARGININE] = 3;    
  // aa_base_h_[ASPARAGINE] = 2;
  // aa_base_h_[ASPARTIC_ACID] = 15;
  // aa_base_h_[CYSTEINE] = 14;
  // aa_base_h_[GLUTAMIC_ACID] = 10;
  // aa_base_h_[GLUTAMINE] = 7; 
  // aa_base_h_[GLYCINE] = 9;
  // aa_base_h_[HISTIDINE] = 8;
  // aa_base_h_[ISOLEUCINE] = 17; 
  // aa_base_h_[LEUCINE] = 12; 
  // aa_base_h_[LYSINE] = 5;       
  // aa_base_h_[METHIONINE] = 19;
  // aa_base_h_[PHENYLALANINE] = 16;    
  // aa_base_h_[PROLINE] = 6;      
  // aa_base_h_[SERINE] = 18;     
  // aa_base_h_[THREONINE] = 11;   
  // aa_base_h_[TRYPTOPHAN] = 1;    
  // aa_base_h_[TYROSINE] = 4;   
  // aa_base_h_[VALINE] = 0;  
  
  
  /* Standard code
aa_base_m_[PHENYLALANINE] = 0;  // 2
  aa_base_m_[LEUCINE] = 1;        // 6
  aa_base_m_[ISOLEUCINE] = 2;     // 3+1
  aa_base_m_[METHIONINE] = 2;
  aa_base_m_[VALINE] = 3;         // 4
  aa_base_m_[SERINE] = 4;         // 6
                                // TOTAL = 2+6+4+4+6 = 22

  aa_base_w_[PROLINE] = 0;        // 4
  aa_base_w_[THREONINE] = 1;      // 4
  aa_base_w_[ALANINE] = 2;        // 4
  aa_base_w_[TYROSINE] = 3;       // 2+2
  aa_base_w_[HISTIDINE] = 3;
  aa_base_w_[GLUTAMINE] = 4;      // 2+2
  aa_base_w_[ASPARAGINE] = 4;
  aa_base_w_[LYSINE] = 5;         // 2
                                // TOTAL = 4+4+4+4+4+2 = 22

  aa_base_h_[ASPARAGINE] = 0;     // 2+2
  aa_base_h_[LYSINE] = 0;
  aa_base_h_[ASPARTIC_ACID] = 1;  // 2+2+2
  aa_base_h_[CYSTEINE] = 1;
  aa_base_h_[GLUTAMIC_ACID] = 1;
  aa_base_h_[TRYPTOPHAN] = 2;     // 1+4
  aa_base_h_[ARGININE] = 3;       // 4
  aa_base_h_[GLYCINE] = 4;        // 4
                                // TOTAL = 4+6+5+4+4 = 23
				*/
  // Compute bases sizes
  aa_base_m_size_ = 0;
  aa_base_w_size_ = 0;
  aa_base_h_size_ = 0;

  for(int8_t i = 0; i < NB_AMINO_ACIDS; i++) {
    if(aa_base_m_[i] > aa_base_m_size_)
      aa_base_m_size_ = aa_base_m_[i];

    if(aa_base_w_[i] > aa_base_w_size_)
      aa_base_w_size_ = aa_base_w_[i];

    if(aa_base_h_[i] > aa_base_h_size_)
      aa_base_h_size_ = aa_base_h_[i];
  }
  aa_base_m_size_++; // base_size = max_digit + 1
  aa_base_w_size_++;
  aa_base_h_size_++;
#endif

  fuzzy_flavor_                 = 0;
  simd_metadata_flavor_         = STD_MAP;

#ifdef __REGUL
  _protein_presence_limit = 1e-2;
  _degradation_rate  = 1;
  _nb_degradation_step  = 10;
  _with_heredity          = false;
  _nb_indiv_age      = 20;

  _hill_shape_n      = 4;
  _hill_shape_theta  = 0.5;
  _hill_shape        = std::pow( _hill_shape_theta, _hill_shape_n );

  _list_eval_step    = new std::set<int>();

  for (int i=0; i < MAX_QUADON; i++) {
		for (int j=0; j < MAX_CODON; j++) {
      _binding_matrix[i][j] = 0.0;
    }
  }
#endif

  #ifdef HAVE_MPI
  nb_rank_ = 1;
  #endif
}


// ===========================================================================
//                                 Destructor
// ===========================================================================
ExpSetup::~ExpSetup() noexcept {
  delete sel_;
#ifdef __REGUL
  delete _list_eval_step;
#endif
}

#ifdef BASE_4

// =====================================================================
//                          Accessors' definitions
// =====================================================================
void ExpSetup::set_aa_base_m(int8_t *base_m) {
  aa_base_m_size_ = 0;
  for(int8_t i = 0; i < NB_AMINO_ACIDS; i++) {
    aa_base_m_[i] = base_m[i];

    if(base_m[i] > aa_base_m_size_)
      aa_base_m_size_ = base_m[i];
  }
  aa_base_m_size_++; // base_size = base_max + 1
}

void ExpSetup::set_aa_base_w(int8_t *base_w) {
  aa_base_w_size_ = 0;
  for(int8_t i = 0; i < NB_AMINO_ACIDS; i++) {
    aa_base_w_[i] = base_w[i];

    if(base_w[i] > aa_base_w_size_)
      aa_base_w_size_ = base_w[i];
  }
  aa_base_w_size_++; // base_size = base_max + 1
}

void ExpSetup::set_aa_base_h(int8_t *base_h) {
  aa_base_h_size_ = 0;
  for(int8_t i = 0; i < NB_AMINO_ACIDS; i++) {
    aa_base_h_[i] = base_h[i];

    if(base_h[i] > aa_base_h_size_)
      aa_base_h_size_ = base_h[i];
  }
  aa_base_h_size_++; // base_size = base_max + 1
}
#endif

/*!
*/
void ExpSetup::write_setup_file(gzFile exp_setup_file) const {
  gzwrite(exp_setup_file,
          &fuzzy_flavor_,sizeof(fuzzy_flavor_));
      gzwrite(exp_setup_file,
              &simd_metadata_flavor_,sizeof(simd_metadata_flavor_));
  // --------------------------------------------------------------- Transfer
  int8_t tmp_with_HT = static_cast<int8_t>(with_HT_);
  int8_t tmp_repl_HT_with_close_points = static_cast<int8_t>(repl_HT_with_close_points_);

  gzwrite(exp_setup_file,
          &tmp_with_HT,
          sizeof(tmp_with_HT));

  gzwrite(exp_setup_file,&tmp_repl_HT_with_close_points,sizeof(tmp_repl_HT_with_close_points));

  if (with_HT_) {
    gzwrite(exp_setup_file,
            &HT_ins_rate_,sizeof(HT_ins_rate_));
    gzwrite(exp_setup_file,&HT_repl_rate_,sizeof(HT_repl_rate_));
  }


  if(repl_HT_with_close_points_)
    gzwrite(exp_setup_file,
            &repl_HT_detach_rate_,sizeof(repl_HT_detach_rate_));


  // --------------------------------------------------------------- Plasmids
  int8_t tmp_with_plasmids = static_cast<int8_t>(with_plasmids());
  gzwrite(exp_setup_file, &tmp_with_plasmids,sizeof(tmp_with_plasmids));

  if (tmp_with_plasmids) {
    gzwrite(exp_setup_file,
            &prob_plasmid_HT_, sizeof(prob_plasmid_HT_));
    gzwrite(exp_setup_file,&tune_donor_ability_,sizeof(tune_donor_ability_));
    gzwrite(exp_setup_file,&tune_recipient_ability_,sizeof(tune_recipient_ability_));
    gzwrite(exp_setup_file,&donor_cost_,sizeof(donor_cost_));
    gzwrite(exp_setup_file,&recipient_cost_,sizeof(recipient_cost_));
    int tmp_swap_GUs_ = static_cast<int8_t>(swap_GUs_);
    gzwrite(exp_setup_file,&tmp_swap_GUs_,sizeof(tmp_swap_GUs_));
  }

  // -------------------------------------------------------------- Secretion
  int8_t tmp_with_secretion_ = static_cast<int8_t>(with_secretion_);
  gzwrite(exp_setup_file,
          &tmp_with_secretion_,sizeof(tmp_with_secretion_));
  gzwrite(exp_setup_file,&secretion_contrib_to_fitness_,sizeof(secretion_contrib_to_fitness_));
  gzwrite(exp_setup_file,&secretion_cost_,sizeof(secretion_cost_));

  sel()->write_setup_file(exp_setup_file);

  mut_params_->save(exp_setup_file);

  #ifdef BASE_4
  gzwrite(exp_setup_file, &terminator_polya_sequence_length_, sizeof(terminator_polya_sequence_length_));
  #endif

  // ---------------------------------------------- Retrieve SIMD parameters
  gzwrite(exp_setup_file, &min_genome_length_,sizeof(min_genome_length_));
  gzwrite(exp_setup_file, &max_genome_length_,sizeof(max_genome_length_));


#ifdef __REGUL
  //printf("Nb indiv age : %d\n",_nb_indiv_age);
  gzwrite( exp_setup_file,&_hill_shape,sizeof(_hill_shape));
  gzwrite( exp_setup_file, &_hill_shape_n,sizeof(_hill_shape_n));
  gzwrite( exp_setup_file, &_hill_shape_theta,sizeof(_hill_shape_theta));
  gzwrite( exp_setup_file, &_degradation_rate,sizeof(_degradation_rate));
  gzwrite( exp_setup_file,  &_nb_degradation_step,sizeof(_nb_degradation_step));
  gzwrite( exp_setup_file, &_nb_indiv_age,sizeof(_nb_indiv_age));
  gzwrite( exp_setup_file, &_with_heredity,sizeof(_with_heredity));
  gzwrite( exp_setup_file, &_protein_presence_limit,sizeof(_protein_presence_limit));

  char* binding_matrix_file_name = new char[100];

  sprintf( binding_matrix_file_name, "binding_matrix.rae" );

  gzFile binding_matrix_file = gzopen( binding_matrix_file_name, "w" );

  if ( binding_matrix_file == Z_NULL )
  {
    printf( "ERROR : Could not write binding matrix file %s\n", binding_matrix_file_name );
    exit( EXIT_FAILURE );
  }

  write_binding_matrix_to_backup( binding_matrix_file );
  gzclose( binding_matrix_file );

  delete[] binding_matrix_file_name;

  unsigned int eval_step_size = _list_eval_step->size();
  gzwrite(exp_setup_file, &eval_step_size,sizeof(eval_step_size));

  for(auto eval_step : *_list_eval_step) {
    gzwrite(exp_setup_file, &eval_step,sizeof(eval_step));
  }
#endif

#ifdef HAVE_MPI
  gzwrite( exp_setup_file, &nb_rank_,sizeof(nb_rank_));

  gzwrite( exp_setup_file, &global_pop_size_,sizeof(global_pop_size_));
  gzwrite( exp_setup_file, &global_grid_width_,sizeof(global_grid_width_));
  gzwrite( exp_setup_file, &global_grid_height_,sizeof(global_grid_height_));

  gzwrite( exp_setup_file, &rank_width_,sizeof(rank_width_));
  gzwrite( exp_setup_file, &rank_height_,sizeof(rank_height_));
#endif
}

void ExpSetup::load(gzFile setup_file, gzFile backup_file, bool verbose) {
  gzread(setup_file,&fuzzy_flavor_,sizeof(fuzzy_flavor_));
  gzread(setup_file,&simd_metadata_flavor_,sizeof(simd_metadata_flavor_));
  // -------------------------------------------- Retrieve transfer parameters
  int8_t tmp_with_HT;
  int8_t tmp_repl_HT_with_close_points;


  gzread(setup_file,
         &tmp_with_HT,sizeof(tmp_with_HT)
  );

  gzread(setup_file,
         &tmp_repl_HT_with_close_points,sizeof(tmp_repl_HT_with_close_points)
  );

  with_HT_ = static_cast<bool>(tmp_with_HT);

  repl_HT_with_close_points_ = static_cast<bool>(tmp_repl_HT_with_close_points);



  if (with_HT_)
  {
    gzread(setup_file,
           &HT_ins_rate_,sizeof(HT_ins_rate_));
    gzread(setup_file,&HT_repl_rate_,sizeof(HT_repl_rate_));
  }
   if(repl_HT_with_close_points_)
    gzread(setup_file, &repl_HT_detach_rate_,sizeof(repl_HT_detach_rate_));


  // -------------------------------------------- Retrieve plasmid parameters
  int8_t tmp_with_plasmids;
  gzread(setup_file,&tmp_with_plasmids,sizeof(tmp_with_plasmids));
  with_plasmids_ = static_cast<bool>(tmp_with_plasmids);

  if (with_plasmids_)
  {
    int8_t tmp_swap_GUs;
    gzread(setup_file,
           &prob_plasmid_HT_,sizeof(prob_plasmid_HT_));
    gzread(setup_file,&tune_donor_ability_,sizeof(tune_donor_ability_));
    gzread(setup_file,&tune_recipient_ability_,sizeof(tune_recipient_ability_));
    gzread(setup_file,&donor_cost_,sizeof(donor_cost_));
    gzread(setup_file,&recipient_cost_,sizeof(recipient_cost_));
    gzread(setup_file,&tmp_swap_GUs,sizeof(tmp_swap_GUs));
    swap_GUs_ = static_cast<bool>(tmp_swap_GUs);
  }

  // ------------------------------------------ Retrieve secretion parameters
  int8_t tmp_with_secretion;
  gzread(setup_file, &tmp_with_secretion,sizeof(tmp_with_secretion));
  gzread(setup_file, &secretion_contrib_to_fitness_,sizeof(secretion_contrib_to_fitness_));
  gzread(setup_file, &secretion_cost_,sizeof(secretion_cost_));
  with_secretion_ = static_cast<bool>(tmp_with_secretion);
  // ---------------------------------------------- Retrieve selection context

  sel()->load(setup_file, backup_file, verbose);

  mut_params_ = std::make_shared<MutationParams>(setup_file);

  #ifdef BASE_4
  // ---------------------------------------------- Retrieve terminators configuration
  gzread(setup_file, &terminator_polya_sequence_length_, sizeof(terminator_polya_sequence_length_));
  #endif

  // ---------------------------------------------- Retrieve SIMD parameters
  gzread(setup_file, &min_genome_length_,sizeof(min_genome_length_));
  gzread(setup_file, &max_genome_length_,sizeof(max_genome_length_));

#ifdef __REGUL
  gzread( setup_file, &_hill_shape,sizeof(_hill_shape));
  gzread( setup_file, &_hill_shape_n,sizeof(_hill_shape_n));
  gzread( setup_file, &_hill_shape_theta,sizeof(_hill_shape_theta));
  gzread( setup_file, &_degradation_rate,sizeof(_degradation_rate));
  gzread( setup_file, &_nb_degradation_step,sizeof(_nb_degradation_step));
  gzread( setup_file, &_nb_indiv_age,sizeof(_nb_indiv_age));
  gzread( setup_file, &_with_heredity,sizeof(_with_heredity));
  gzread( setup_file, &_protein_presence_limit,sizeof(_protein_presence_limit));


  //printf("Nb indiv age : %d\n",_nb_indiv_age);

  char* binding_matrix_file_name = new char[100];
//    _binding_matrix = new double[MAX_QUADON][MAX_CODON];

  sprintf( binding_matrix_file_name, "binding_matrix.rae" );

  gzFile binding_matrix_file = gzopen( binding_matrix_file_name, "r" );

  if ( binding_matrix_file == Z_NULL )
  {
    printf( "ERROR : Could not read binding matrix file %s\n", binding_matrix_file_name );
    exit( EXIT_FAILURE );
  }

  read_binding_matrix_from_backup( binding_matrix_file );
  gzclose( binding_matrix_file );

  delete[] binding_matrix_file_name;

  unsigned int eval_step_size;
  gzread(setup_file, &eval_step_size,sizeof(eval_step_size));

  int eval_val;
  for(unsigned int i = 0; i < eval_step_size; i++) {
    gzread(setup_file,&eval_val,sizeof(eval_val));
    _list_eval_step->insert(eval_val);
  }
#endif
#ifdef HAVE_MPI
  gzread( setup_file, &nb_rank_,sizeof(nb_rank_));

  gzread( setup_file, &global_pop_size_,sizeof(global_pop_size_));
  gzread( setup_file, &global_grid_width_,sizeof(global_grid_width_));
  gzread( setup_file, &global_grid_height_,sizeof(global_grid_height_));

  gzread( setup_file, &rank_width_,sizeof(rank_width_));
  gzread( setup_file, &rank_height_,sizeof(rank_height_));
#endif
}

#ifdef __REGUL
void ExpSetup::init_binding_matrix( bool random_binding_matrix, double binding_zeros_percentage,
		std::shared_ptr<JumpingMT> prng)
{
  if(random_binding_matrix==1)
  {
    // printf("Generate the matrix %lf\n",binding_zeros_percentage);
    for( int8_t i = 0; i < MAX_QUADON; i++ )  // i for the quadons
    {
      for( int8_t j = 0; j < MAX_CODON; j++ )  // j for the codons
      {
        double rand = prng->random();
        // printf("PRNG %lf\n",rand);
        if( rand > binding_zeros_percentage)
        {
        	_binding_matrix[i][j] = (ProteinConcentration) prng->random();
           	  // printf("m[%d][%d] = %f\n",i,j, _binding_matrix[i][j]);

        }
        else
        {
        	_binding_matrix[i][j] = 0;
        }
      }
    }
  }
  else // random_binding_matrix == 0
  {
    char* binding_matrix_file_name = new char[100];
//    _binding_matrix = new double[MAX_QUADON][MAX_CODON];

    sprintf( binding_matrix_file_name, "binding_matrix.rae" );

    gzFile binding_matrix_file = gzopen( binding_matrix_file_name, "r" );

    if ( binding_matrix_file == Z_NULL )
    {
      printf( "ERROR : Could not read binding matrix file %s\n", binding_matrix_file_name );
      exit( EXIT_FAILURE );
    }

    read_binding_matrix_from_backup( binding_matrix_file );
    gzclose( binding_matrix_file );

    delete[] binding_matrix_file_name;
  }




}

void ExpSetup::read_binding_matrix_from_backup(gzFile binding_matrix_file) {
  double value;
	for (int i=0; i < MAX_QUADON; i++) {
		for (int j=0; j < MAX_CODON; j++) {
			gzread( binding_matrix_file, &value, sizeof(value));
			_binding_matrix[i][j] = (ProteinConcentration) value;
		}
	}
}

void ExpSetup::write_binding_matrix_to_backup(gzFile binding_matrix_file) const {
	double value;
	for (int i=0; i < MAX_QUADON; i++) {
		for (int j=0; j < MAX_CODON; j++) {
			value = (double)_binding_matrix[i][j];
			gzwrite( binding_matrix_file,&value,sizeof(value));
		}
	}
}

void ExpSetup::write_binding_matrix_to_file( FILE* file ) const
{
  for( int16_t row = 0 ; row < MAX_QUADON ; row++ )
  {
    for( int16_t column = 0 ; column < MAX_CODON ; column++ )
    {
      fprintf( file, "\t%e", _binding_matrix[row][column] );
    }
    fprintf( file, "\n");
  }
}

void ExpSetup::print_binding_matrix()
{
  for( int16_t row = 0 ; row < MAX_QUADON ; row++ )
  {
    for( int16_t column = 0 ; column < MAX_CODON ; column++ )
    {
      printf("M[%d][%d] = %e\n", row, column, _binding_matrix[row][column] );
    }
  }
}
#endif

// ===========================================================================
//                                Protected Methods
// ===========================================================================
} // namespace aevol
