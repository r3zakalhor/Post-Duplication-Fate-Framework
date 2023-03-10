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
//#include <math.h>

#include <list>

// =================================================================
//                            Project Files
// =================================================================
#include "Protein.h"

#include "ExpManager.h"
#include "Codon.h"
#include "Individual.h"
#include "GeneticUnit.h"
#include "Rna.h"
#include "Utils.h"

namespace aevol {



//##############################################################################
//                                                                             #
//                              Class Protein                               #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
/**
 * Copy constructor.
 *
 * Copies the protein but does nothing regarding the RNAs transcribing it (creates an empty list).
 */
Protein::Protein(GeneticUnit* gen_unit, const Protein &model)
{
  gen_unit_  = gen_unit;

  strand_                 = model.strand_;
  shine_dal_pos_          = model.shine_dal_pos_;
  first_translated_pos_   = model.first_translated_pos_;
  last_translated_pos_    = model.last_translated_pos_;
  length_                 = model.length_;
  concentration_          = model.concentration_;
  is_functional_          = model.is_functional_;

  // Copy the list of amino-acids

  // TODO vld: check if deep copy needed (at least, in case of heritage yes)
  //AA_list_ = model.AA_list_;
          for (auto aa_c : model.AA_list_) {
            #ifdef BASE_2
            AA_list_.push_back(new Codon(aa_c->value()));
            #elif BASE_4
            AA_list_.push_back(new Codon(*aa_c));
            #endif
          }

  // Copy triangle parameters
  mean_   = model.mean_;
  width_  = model.width_;
  height_ = model.height_;
}

Protein::Protein(GeneticUnit* gen_unit, const Protein &model, ExpManager* exp_m)
{
  gen_unit_  = gen_unit;

  strand_                 = model.strand_;
  shine_dal_pos_          = model.shine_dal_pos_;
  first_translated_pos_   = model.first_translated_pos_;
  last_translated_pos_    = model.last_translated_pos_;
  length_                 = model.length_;
  concentration_          = model.concentration_;
  is_functional_          = model.is_functional_;

  // Copy triangle parameters
  mean_   = model.mean_;
  width_  = model.width_;
  height_ = model.height_;
}

Protein::Protein(GeneticUnit* gen_unit,
                 const std::list<Codon*>& codon_list,
                 Strand strand,
                 int32_t shine_dal_pos,
                 Rna* rna,
                 double w_max)
{
  assert(shine_dal_pos >= 0);
  assert(shine_dal_pos < gen_unit->seq_length());

  gen_unit_       = gen_unit;
  strand_         = strand;
  shine_dal_pos_  = shine_dal_pos;
  length_         = codon_list.size();

  #ifndef __REGUL
    // In Aevol the concentration of a new protein is set at the basal level
    concentration_  = rna->basal_level();
  #else
    // In Raevol, there is two case, depending on the heredity
    if (gen_unit_->exp_m()->exp_s()->get_with_heredity() )
    {
      // With heredity the new protein has a concentration set at 0, because there are inherited proteins which allow the regulation
      concentration_ = 0;
    }
    else
    {
      // Without heredity, we use the same concentration as in Aevol (No inherited proteins)
      concentration_ = rna->basal_level();
    }
  #endif

  // TODO : make this cleaner...
  AA_list_ = codon_list;

  rna_list_.push_back(rna);

  if (strand_ == LEADING)
  {
    first_translated_pos_ = Utils::mod(shine_dal_pos_ + (SHINE_DAL_SIZE + SHINE_START_SPACER + CODON_SIZE),
                                        gen_unit_->dna()->length());
    last_translated_pos_  = Utils::mod(first_translated_pos_ + (length_ * CODON_SIZE - 1),
                                        gen_unit_->dna()->length());
  }
  else
  {
    first_translated_pos_ = Utils::mod(shine_dal_pos_ - (SHINE_DAL_SIZE + SHINE_START_SPACER + CODON_SIZE),
                                        gen_unit_->dna()->length());
    last_translated_pos_ = Utils::mod(first_translated_pos_ - (length_ * CODON_SIZE - 1),
                                       gen_unit_->dna()->length());
  }



  // ============================================================================
  // Folding process (compute phenotypic contribution parameters from codon list)
  // ============================================================================
  //  1) Compute values for M, W and H
  //  2) Normalize M, W and H values according to number of codons of each kind
  //  3) Normalize M, W and H values according to the allowed ranges (defined in macros.h)
  #ifdef BASE_4
  auto base_m = gen_unit_->exp_m()->exp_s()->aa_base_m();
  auto base_w = gen_unit_->exp_m()->exp_s()->aa_base_w();
  auto base_h = gen_unit_->exp_m()->exp_s()->aa_base_h();

  int8_t base_m_size = gen_unit_->exp_m()->exp_s()->aa_base_m_size();
  int8_t base_w_size = gen_unit_->exp_m()->exp_s()->aa_base_w_size();
  int8_t base_h_size = gen_unit_->exp_m()->exp_s()->aa_base_h_size();


  double m_digit_factor = 1.0;
  double w_digit_factor = 1.0;
  double h_digit_factor = 1.0;
  #endif

  //  --------------------------------
  //  1) Compute values for M, W and H
  //  --------------------------------
  long double M = 0.0;
  long double W = 0.0;
  long double H = 0.0;

  int32_t nb_m = 0;
  int32_t nb_w = 0;
  int32_t nb_h = 0;

  #ifdef BASE_2
  bool bin_m = false; // Initializing to false will yield a conservation of the high weight bit
  bool bin_w = false; // when applying the XOR operator for the Gray to standard conversion
  bool bin_h = false;

  for (const auto& codon: codon_list) {
    switch (codon->value())
    {
      case CODON_M0 :
      {
        // M codon found
        nb_m++;

        // Convert Gray code to "standard" binary code
        bin_m ^= false; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

        // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
        //~ M <<= 1;
        M *= 2;

        // Add this nucleotide's contribution to M
        if (bin_m) M += 1;

        break;
      }
      case CODON_M1 :
      {
        // M codon found
        nb_m++;

        // Convert Gray code to "standard" binary code
        bin_m ^= true; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

        // A lower-than-the-previous-lowest bit was found, make a left bitwise shift
        //~ M <<= 1;
        M *= 2;

        // Add this nucleotide's contribution to M
        if (bin_m) M += 1;

        break;
      }
      case CODON_W0 :
      {
        // W codon found
        nb_w++;

        // Convert Gray code to "standard" binary code
        bin_w ^= false; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

        // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
        //~ W <<= 1;
        W *= 2;

        // Add this nucleotide's contribution to W
        if (bin_w) W += 1;

        break;
      }
      case CODON_W1 :
      {
        // W codon found
        nb_w++;

        // Convert Gray code to "standard" binary code
        bin_w ^= true; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

        // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
        //~ W <<= 1;
        W *= 2;

        // Add this nucleotide's contribution to W
        if (bin_w) W += 1;

        break;
      }
      case CODON_H0 :
      case CODON_START : // Start codon codes for the same amino-acid as H0 codon
      {
        // H codon found
        nb_h++;

        // Convert Gray code to "standard" binary code
        bin_h ^= false; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

        // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
        //~ H <<= 1;
        H *= 2;

        // Add this nucleotide's contribution to H
        if (bin_h) H += 1;

        break;
      }
      case CODON_H1 :
      {
        // H codon found
        nb_h++;

        // Convert Gray code to "standard" binary code
        bin_h ^= true; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

        // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
        //~ H <<= 1;
        H *= 2;

        // Add this nucleotide's contribution to H
        if (bin_h) H += 1;

        break;
      }
    }
  }
  #elif BASE_4
  for(auto codon_it = codon_list.rbegin(); codon_it != codon_list.rend(); codon_it++) {
    AminoAcid amino_acid = (*codon_it)->to_aminoacid();

    if(base_m[amino_acid] != -1) {
      M += base_m[amino_acid] * m_digit_factor;
      m_digit_factor *= base_m_size;
      nb_m++;
    }

    if(base_w[amino_acid] != -1) {
      W += base_w[amino_acid] * w_digit_factor;
      w_digit_factor *= base_w_size;
      nb_w++;
    }

    if(base_h[amino_acid] != -1) {
      H += base_h[amino_acid] * h_digit_factor;
      h_digit_factor *= base_h_size;
      nb_h++;
    }
  }
  #endif


  //  ----------------------------------------------------------------------------------
  //  2) Normalize M, W and H values in [0;1] according to number of codons of each kind
  //  ----------------------------------------------------------------------------------
  #ifdef BASE_2
  if (nb_m != 0)  mean_ = M / (pow(2, nb_m) - 1);
  else              mean_ = 0.5;
  if (nb_w != 0)  width_ = W / (pow(2, nb_w) - 1);
  else              width_ = 0.0;
  if (nb_h != 0)  height_ = H / (pow(2, nb_h) - 1);
  else              height_ = 0.5;
  #elif BASE_4
  if(nb_m != 0) {
    mean_ = M / (pow(base_m_size, nb_m) - 1);
  } else {
    mean_ = 0.5;
  }

  if(nb_w != 0) {
    width_ = W / (pow(base_w_size, nb_w) - 1);
  }
  else {
    width_ = 0.0;
  }

  if(nb_h != 0) {
    height_ = H / (pow(base_h_size, nb_h) - 1);
  } else {
    height_ = 0.5;
  }
  #endif
  
  assert(mean_ >= 0.0 && mean_ <= 1.0);
  assert(width_ >= 0.0 && width_ <= 1.0);
  assert(height_ >= 0.0 && height_ <= 1.0);



  //  ------------------------------------------------------------------------------------
  //  3) Normalize M, W and H values according to the allowed ranges (defined in macros.h)
  //  ------------------------------------------------------------------------------------
  // x_min <= M <= x_max
  // w_min <= W <= w_max
  // h_min <= H <= h_max
  mean_   = (X_MAX - X_MIN) * mean_ + X_MIN;
  width_  = (w_max - W_MIN) * width_ + W_MIN;
  height_ = (H_MAX - H_MIN) * height_ + H_MIN;

  if (nb_m == 0 || nb_w == 0 || nb_h == 0 || width_ == 0.0 || height_ == 0.0)
  {
    is_functional_ = false;
  }
  else
  {
    is_functional_ = true;
  }

  assert(mean_ >= X_MIN && mean_ <= X_MAX);
  assert(width_ >= W_MIN && width_ <= indiv()->w_max());
  assert(height_ >= H_MIN && height_ <= H_MAX);
}

/*
Protein::Protein(Protein* parent)
{
  gen_unit_             = parent->gen_unit_;
  strand_               = parent->strand_;
  shine_dal_pos_        = parent->shine_dal_pos_;
  first_translated_pos_ = parent->first_translated_pos_;
  last_translated_pos_  = parent->last_translated_pos_;
  length_               = parent->length_;
  concentration_        = parent->concentration_;
  mean_                 = parent->mean_;
  width_                = parent->width_;
  height_               = parent->height_;
}
*/

//Constructor for the signal proteins
//modif raevol_yo_3 : now we really copy the codon list
Protein::Protein(const std::list<Codon*> codon_list, ProteinConcentration concentration, double w_max)
// WARNING this constructor should not being used for other purpose
// In particular codon list should be destroyer by the caller of this constructor
// since we will copy the list to be sure to own it
{
  gen_unit_             = NULL;
  strand_               = LEADING;
  shine_dal_pos_        = 0;
  length_               = codon_list.size();
  first_translated_pos_ = 0;
  last_translated_pos_  = 0;
  concentration_        = concentration;

  // Copy the list of amino-acids
  for(Codon* c : codon_list) {
    AA_list_.push_back(new Codon(*c));
  }



  // ============================================================================
  // Folding process (compute phenotypic contribution parameters from codon list)
  // ============================================================================
  //  1) Compute values for M, W and H
  //  2) Normalize M, W and H values according to number of codons of each kind
  //  3) Normalize M, W and H values according to the allowed ranges (defined in macros.h)
  #ifdef BASE_4
  auto base_m = gen_unit_->exp_m()->exp_s()->aa_base_m();
  auto base_w = gen_unit_->exp_m()->exp_s()->aa_base_w();
  auto base_h = gen_unit_->exp_m()->exp_s()->aa_base_h();

  int8_t base_m_size = gen_unit_->exp_m()->exp_s()->aa_base_m_size();
  int8_t base_w_size = gen_unit_->exp_m()->exp_s()->aa_base_w_size();
  int8_t base_h_size = gen_unit_->exp_m()->exp_s()->aa_base_h_size();
  #endif

  //  --------------------------------
   //  1) Compute values for M, W and H
   //  --------------------------------
   long double M = 0.0;
   long double W = 0.0;
   long double H = 0.0;

   int32_t nb_m = 0;
   int32_t nb_w = 0;
   int32_t nb_h = 0;

   #ifdef BASE_2
   bool bin_m = false; // Initializing to false will yield a conservation of the high weight bit
   bool bin_w = false; // when applying the XOR operator for the Gray to standard conversion
   bool bin_h = false;

   for (const auto& codon: codon_list) {
     switch ( codon->value() )
     {
       case CODON_M0 :
       {
         // M codon found
         nb_m++;

         // Convert Gray code to "standard" binary code
         bin_m ^= false; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

         // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
         //~ M <<= 1;
         M *= 2;

         // Add this nucleotide's contribution to M
         if ( bin_m ) M += 1;

         break;
       }
       case CODON_M1 :
       {
         // M codon found
         nb_m++;

         // Convert Gray code to "standard" binary code
         bin_m ^= true; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

         // A lower-than-the-previous-lowest bit was found, make a left bitwise shift
         //~ M <<= 1;
         M *= 2;

         // Add this nucleotide's contribution to M
         if ( bin_m ) M += 1;

         break;
       }
       case CODON_W0 :
       {
         // W codon found
         nb_w++;

         // Convert Gray code to "standard" binary code
         bin_w ^= false; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

         // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
         //~ W <<= 1;
         W *= 2;

         // Add this nucleotide's contribution to W
         if ( bin_w ) W += 1;

         break;
       }
       case CODON_W1 :
       {
         // W codon found
         nb_w++;

         // Convert Gray code to "standard" binary code
         bin_w ^= true; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

         // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
         //~ W <<= 1;
         W *= 2;

         // Add this nucleotide's contribution to W
         if ( bin_w ) W += 1;

         break;
       }
       case CODON_H0 :
       case CODON_START : // Start codon codes for the same amino-acid as H0 codon
       {
         // H codon found
         nb_h++;

         // Convert Gray code to "standard" binary code
         bin_h ^= false; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

         // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
         //~ H <<= 1;
         H *= 2;

         // Add this nucleotide's contribution to H
         if ( bin_h ) H += 1;

         break;
       }
       case CODON_H1 :
       {
         // H codon found
         nb_h++;

         // Convert Gray code to "standard" binary code
         bin_h ^= true; // as bin_m was initialized to false, the XOR will have no effect on the high weight bit

         // A lower-than-the-previous-lowest weight bit was found, make a left bitwise shift
         //~ H <<= 1;
         H *= 2;

         // Add this nucleotide's contribution to H
         if ( bin_h ) H += 1;

         break;
       }
     }
   }
   #elif BASE_4
   /*
  double m_digit_factor = 1.0;
  double w_digit_factor = 1.0;
  double h_digit_factor = 1.0;

  for(auto codon_it = codon_list.rbegin(); codon_it != codon_list.rend(); codon_it++) {
    AminoAcid amino_acid = (*codon_it)->to_aminoacid();

    if(base_m[amino_acid] != -1) {
      M += base_m[amino_acid] * m_digit_factor;
      m_digit_factor *= base_m_size;
      nb_m++;
    }

    if(base_w[amino_acid] != -1) {
      W += base_w[amino_acid] * w_digit_factor;
      w_digit_factor *= base_w_size;
      nb_w++;
    }

    if(base_h[amino_acid] != -1) {
      H += base_h[amino_acid] * h_digit_factor;
      h_digit_factor *= base_h_size;
      nb_h++;
    }
    }*/

    // new version with long genes

  
  double m_digit_factor = 1.0;
  double w_digit_factor = 1.0;
  double h_digit_factor = 1.0;

  
  double m_base_exponent=1.25;
  double w_base_exponent=1.25;
  double h_base_exponent=1.25;
  

  for(auto codon_it = codon_list.rbegin(); codon_it != codon_list.rend(); codon_it++) {
    AminoAcid amino_acid = (*codon_it)->to_aminoacid();

    if(base_m[amino_acid] != -1) {
      M += base_m[amino_acid] * m_digit_factor;
      m_digit_factor *= base_m_size;
      nb_m++;
    }

    if(base_w[amino_acid] != -1) {
      W += base_w[amino_acid] * w_digit_factor;
      w_digit_factor *= base_w_size;
      nb_w++;
    }

    if(base_h[amino_acid] != -1) {
      H += base_h[amino_acid] * h_digit_factor;
      h_digit_factor *= base_h_size;
      nb_h++;
    }
  }

   #endif



   //  ----------------------------------------------------------------------------------
   //  2) Normalize M, W and H values in [0;1] according to number of codons of each kind
   //  ----------------------------------------------------------------------------------
   #ifdef BASE_2
   if ( nb_m != 0 )  mean_ = M / (pow(2, nb_m) - 1);
   else              mean_ = 0.5;
   if ( nb_w != 0 )  width_ = W / (pow(2, nb_w) - 1);
   else              width_ = 0.0;
   if ( nb_h != 0 )  height_ = H / (pow(2, nb_h) - 1);
   else              height_ = 0.5;
   #elif BASE_4

   /*
  if(nb_m != 0) {
    mean_ = M / (pow(base_m_size, nb_m) - 1);
  } else {
    mean_ = 0.5;
  }

  if(nb_w != 0) {
    width_ = W / (pow(base_w_size, nb_w) - 1);
  }
  else {
    width_ = 0.0;
  }

  if(nb_h != 0) {
    height_ = H / (pow(base_h_size, nb_h) - 1);
  } else {
    height_ = 0.5;
    }*/

   // new version with long genes

   
  if(nb_m != 0) {
    mean_ = M / ((base_m_size)*(pow(m_base_exponent, nb_m) - 1)/(m_base_exponent - 1));
  } else {
    mean_ = 0.5;
  }

  if(nb_w != 0) {
    width_ = W / ((base_w_size)*(pow(w_base_exponent, nb_w) - 1)/(w_base_exponent - 1));
  }
  else {
    width_ = 0.0;
  }

  if(nb_h != 0) {
    height_ = H / ((base_h_size)*(pow(h_base_exponent, nb_h) - 1)/(h_base_exponent - 1));
  } else {
    height_ = 0.5;
  }
   
   #endif
   assert( mean_ >= 0.0 && mean_ <= 1.0 );
   assert( width_ >= 0.0 && width_ <= 1.0 );
   assert( height_ >= 0.0 && height_ <= 1.0 );



   //  ------------------------------------------------------------------------------------
   //  3) Normalize M, W and H values according to the allowed ranges (defined in macros.h)
   //  ------------------------------------------------------------------------------------
   // x_min <= M <= x_max
   // w_min <= W <= w_max
   // h_min <= H <= h_max
   mean_   = (X_MAX - X_MIN) * mean_ + X_MIN;
   width_  = (w_max - W_MIN) * width_ + W_MIN;
   height_ = (H_MAX - H_MIN) * height_ + H_MIN;

   if ( nb_m == 0 || nb_w == 0 || nb_h == 0 || width_ == 0.0 || height_ == 0.0 )
   {
     is_functional_ = false;
   }
   else
   {
     is_functional_ = true;
   }

   assert( mean_ >= X_MIN && mean_ <= X_MAX );
   //NOT WORKING if indiv is nullptr
  // assert( width_ >= W_MIN && width_ <= indiv()->w_max() );
   assert( height_ >= H_MIN && height_ <= H_MAX );
}

//Constructor for the cloning of signal proteins
//modif raevol_yo_3 : now we really copy the codon list
Protein::Protein(Protein* signal)
// WARNING this constructor should not being used for other purpose
// In particular codon list should be destroyer by the caller of this constructor
// since we will copy the list to be sure to own it
{
  gen_unit_             = NULL;
  strand_               = LEADING;
  shine_dal_pos_        = 0;
  length_               = signal->AA_list_.size();
  first_translated_pos_ = 0;
  last_translated_pos_  = 0;
  concentration_        = signal->concentration();

  // Copy the list of amino-acids
  for(Codon* c : signal->AA_list_) {
    AA_list_.push_back(new Codon(*c));
  }


  mean_   = signal->mean_;
  width_  = signal->width_;
  height_ = signal->height_;

 is_functional_ = signal->is_functional_;

}

Protein::Protein(gzFile backup_file)
{
  gen_unit_ = NULL;
  int8_t tmp_strand;
  gzread(backup_file, &tmp_strand, sizeof(tmp_strand));
  strand_ = (Strand) tmp_strand;
  gzread(backup_file, &shine_dal_pos_,			    sizeof(shine_dal_pos_));
  gzread(backup_file, &first_translated_pos_, 	sizeof(first_translated_pos_));
  gzread(backup_file, &last_translated_pos_,  	sizeof(last_translated_pos_));
  gzread(backup_file, &length_,     			      sizeof(length_));

  double tmp;
  gzread(backup_file, &tmp,     		sizeof(tmp));
  concentration_ = (ProteinConcentration) tmp;

  gzread(backup_file, &is_functional_,         sizeof(is_functional_));

  gzread(backup_file, &tmp,  			          sizeof(tmp));
  mean_ = (ProteinConcentration) tmp;

  gzread(backup_file, &tmp,    			        sizeof(tmp));
  width_ = (ProteinConcentration) tmp;

  gzread(backup_file, &tmp,                sizeof(tmp));
  height_ = (ProteinConcentration) tmp;

  // Retreive the AA
  int16_t nb_AA = 0;
  gzread(backup_file, &nb_AA,  sizeof(nb_AA));

  for (int16_t i = 0 ; i < nb_AA ; i++)
    AA_list_.push_back(new Codon(backup_file));

}

// =================================================================
//                             Destructors
// =================================================================
Protein::~Protein()
{
  for (const auto& AA: AA_list_)
    delete AA;
}

// =================================================================
//                            Public Methods
// =================================================================
int32_t Protein::last_STOP_base_pos() const
{
  if (strand_ == LEADING)
  {
    return Utils::mod(last_translated_pos_ + 3, gen_unit_->dna()->length());
  }
  else
  {
    return Utils::mod(last_translated_pos_ - 3, gen_unit_->dna()->length());
  }
}

/*
void Protein::copy_parent(Rna* rna) {
  rna_list_.push_back(rna);
}*/

void Protein::add_RNA(Rna * rna)
{

  rna_list_.push_back(rna);

#ifndef __REGUL
  concentration_ += rna->basal_level();
#else
	if ( gen_unit_->exp_m()->exp_s()->get_with_heredity() )
	{
	  // With heredity the new protein has a concentration set at 0, because there are inherited proteins which allow the regulation
	 concentration_ = 0;
	}
	else
	{
	  // Without heredity, we use the same concentration as in Aevol (No inherited proteins)
	  concentration_ += rna->basal_level();
	}
#endif
}

char* Protein::AA_sequence(char separator /*= ' '*/) const
{
  #ifdef BASE_2
  char* seq = new char[3*length_]; // + 1 (for the '\0')  - 1 (length_ - 1 spaces)

  int32_t i = 0;
  for (const auto& codon: AA_list_) {
    if (i != 0) seq[i++] = separator;
    switch (codon->value())
    {
      case CODON_START :
      {
        seq[i++] = 'S';
        seq[i++] = 'T';
        break;
      }
      case CODON_M0 :
      {
        seq[i++] = 'M';
        seq[i++] = '0';
        break;
      }
      case CODON_M1 :
      {
        seq[i++] = 'M';
        seq[i++] = '1';
        break;
      }
      case CODON_W0 :
      {
        seq[i++] = 'W';
        seq[i++] = '0';
        break;
      }
      case CODON_W1 :
      {
        seq[i++] = 'W';
        seq[i++] = '1';
        break;
      }
      case CODON_H0 :
      {
        seq[i++] = 'H';
        seq[i++] = '0';
        break;
      }
      case CODON_H1 :
      {
        seq[i++] = 'H';
        seq[i++] = '1';
        break;
      }
    }
  }

  seq[3*length_-1] = '\0';
  #elif BASE_4
  size_t seq_length = 0;

  auto base_m = gen_unit_->exp_m()->exp_s()->aa_base_m();
  auto base_w = gen_unit_->exp_m()->exp_s()->aa_base_w();
  auto base_h = gen_unit_->exp_m()->exp_s()->aa_base_h();


  AminoAcid amino_acid;
  for (const auto &codon : AA_list_) {
    amino_acid = codon->amino_acid();
    if (base_m[amino_acid] != -1) {
      seq_length += 3; // 2 chars + 1 separator  (1st separator counts for \0 terminator)
    }
    if (base_w[amino_acid] != -1) {
      seq_length += 3; // 2 chars + 1 separator  (1st separator counts for \0 terminator)
    }
    if (base_h[amino_acid] != -1) {
      seq_length += 3; // 2 chars + 1 separator  (1st separator counts for \0 terminator)
    }
  }

  char* seq = new char[seq_length];
  seq[seq_length - 1] = '\0';

  int32_t i = 0;
  for (const auto& codon: AA_list_) {

    amino_acid = codon->amino_acid();
    if (base_m[amino_acid] != -1) {
      if (i != 0) seq[i++] = separator;
      seq[i++] = 'M';
      seq[i++] = '0' + base_m[amino_acid];
    }
    if (base_w[amino_acid] != -1) {
      if (i != 0) seq[i++] = separator;
      seq[i++] = 'W';
      seq[i++] = '0' + base_w[amino_acid];
    }
    if (base_h[amino_acid] != -1) {
      if (i != 0) seq[i++] = separator;
      seq[i++] = 'H';
      seq[i++] = '0' + base_h[amino_acid];
    }
  }
  #endif
  return seq;
}

void Protein::save(gzFile backup_file)
{
  // The rna_list_ is not write because there is no need to, it is an empty list.
  int8_t tmp_strand = strand_;
  gzwrite(backup_file, &tmp_strand,            sizeof(tmp_strand));
  gzwrite(backup_file, &shine_dal_pos_,        sizeof(shine_dal_pos_));
  gzwrite(backup_file, &first_translated_pos_, sizeof(first_translated_pos_));
  gzwrite(backup_file, &last_translated_pos_,  sizeof(last_translated_pos_));
  gzwrite(backup_file, &length_,     			    sizeof(length_));
  double tmp = concentration_;

  gzwrite(backup_file, &tmp,     		sizeof(tmp));
  gzwrite(backup_file, &is_functional_,        sizeof(is_functional_));

  tmp = mean_;
  gzwrite(backup_file, &tmp,     		sizeof(tmp));

  tmp = width_;
  gzwrite(backup_file, &tmp,     		sizeof(tmp));

  tmp = height_;
  gzwrite(backup_file, &tmp,     		sizeof(tmp));

  // Write the Acide Amino in the backup file
  int16_t nb_AA = AA_list_.size();
  gzwrite(backup_file, &nb_AA,  sizeof(nb_AA));

  for (const auto& AA: AA_list_)
    AA->save(backup_file);
}

void Protein::recompute_concentration() {
    concentration_=0.0; for (auto rna : rna_list_) { concentration_+=rna->basal_level();}}
// =================================================================
//                        Overloaded Operators
// =================================================================
    bool Protein::operator<(const Protein & other){
      return (height_ <  other.height_)
             || (height_ == other.height_ && mean_ < other.mean_)
             || (height_ == other.height_ && mean_ == other.mean_ && width_ < other.width_)
             || (height_ == other.height_ && mean_ == other.mean_ && width_ == other.width_ && concentration_ < other.concentration_);
    }

// =================================================================
//                           Protected Methods
// =================================================================

// =================================================================
//                          Non inline accessors
// =================================================================
Individual *Protein::indiv() const
{
  return gen_unit_->indiv();
}


GeneticUnit* Protein::get_gen_unit( void ) const
{
  return gen_unit_;
}
} // namespace aevol
