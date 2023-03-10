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

/** \class
 *  \brief
 */


// =================================================================
//                              Libraries
// =================================================================
#include <assert.h>



// =================================================================
//                            Project Files
// =================================================================

#include "Mutation.h"
#include "GeneMutation.h"
#include "Rna.h"
#include "macros.h"

namespace aevol {

//##############################################################################
//                                                                             #
//                        Class GeneMutation                               #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================



// =================================================================
//                             Constructors
// =================================================================




// Creates a copy of the mutation mut, but enriched with the generation when it occured
// and the position where it occurred in the RNA, relative to the first bp of the promoter
GeneMutation::GeneMutation(Mutation const & mut, int32_t gener, int32_t cdsPosBefore, Strand strandBefore, ae_gene_mutation_region region) : Mutation(mut)
{
  generation_ = gener;
  impact_on_metabolic_error_ = 0.0; /* should be set to its real value when known */
  region_ = region;

  /* Compute position_relative_to_shine_dal_ */

  switch (mut_type_)
    {
    case SWITCH :
      position_relative_to_shine_dal_ = new int32_t;
      if (strandBefore == LEADING) {position_relative_to_shine_dal_[0] = pos_[0] - cdsPosBefore;}
      else                              {position_relative_to_shine_dal_[0] = cdsPosBefore - pos_[0];}
      break;
    case S_INS :
      position_relative_to_shine_dal_ = new int32_t;
      if (strandBefore == LEADING) {position_relative_to_shine_dal_[0] = pos_[0] - cdsPosBefore;}
      else                              {position_relative_to_shine_dal_[0] = cdsPosBefore - pos_[0];}
      break;
    case S_DEL :
      position_relative_to_shine_dal_ = new int32_t;
      if (strandBefore == LEADING) {position_relative_to_shine_dal_[0] = pos_[0] - cdsPosBefore;}
      else                              {position_relative_to_shine_dal_[0] = cdsPosBefore - pos_[0];}
      break;
    case DUPL :
      /* A duplication can affect a gene in two ways:
         1) The reinsertion point of the duplicated segment is located within the gene => stored
         2) The gene is partly or completely duplicated, but this does not change its sequence => nothing to store
      */
      /* We should enter here only in case (1). Note that in this case, the relative positions for beginseg and endseg may be outside the gene */
      position_relative_to_shine_dal_ = new int32_t[3];
      if (strandBefore == LEADING)
        {
          position_relative_to_shine_dal_[0] = pos_[0] - cdsPosBefore;
          position_relative_to_shine_dal_[1] = pos_[1] - cdsPosBefore;
          position_relative_to_shine_dal_[2] = pos_[2] - cdsPosBefore;
        }
      else
        {
          position_relative_to_shine_dal_[0] = cdsPosBefore - pos_[0];
          position_relative_to_shine_dal_[1] = cdsPosBefore - pos_[1];
          position_relative_to_shine_dal_[2] = cdsPosBefore - pos_[2];
        }
      break;
    case DEL :
      position_relative_to_shine_dal_ = new int32_t[2];
      if (strandBefore == LEADING)
        {
          position_relative_to_shine_dal_[0] = pos_[0] - cdsPosBefore;
          position_relative_to_shine_dal_[1] = pos_[1] - cdsPosBefore;
        }
      else
        {
          position_relative_to_shine_dal_[0] = cdsPosBefore - pos_[0];
          position_relative_to_shine_dal_[1] = cdsPosBefore - pos_[1];
        }
      break;
    case TRANS :
      position_relative_to_shine_dal_ = new int32_t[4];
      if (strandBefore == LEADING)
        {
          position_relative_to_shine_dal_[0] = pos_[0] - cdsPosBefore;
          position_relative_to_shine_dal_[1] = pos_[1] - cdsPosBefore;
          position_relative_to_shine_dal_[2] = pos_[2] - cdsPosBefore;
          position_relative_to_shine_dal_[3] = pos_[3] - cdsPosBefore;
        }
      else
        {
          position_relative_to_shine_dal_[0] = cdsPosBefore - pos_[0];
          position_relative_to_shine_dal_[1] = cdsPosBefore - pos_[1];
          position_relative_to_shine_dal_[2] = cdsPosBefore - pos_[2];
          position_relative_to_shine_dal_[3] = cdsPosBefore - pos_[3];
        }
      break;
    case INV :
      position_relative_to_shine_dal_ = new int32_t[2];
      if (strandBefore == LEADING)
        {
          position_relative_to_shine_dal_[0] = pos_[0] - cdsPosBefore;
          position_relative_to_shine_dal_[1] = pos_[1] - cdsPosBefore;
        }
      else
        {
          position_relative_to_shine_dal_[0] = cdsPosBefore - pos_[0];
          position_relative_to_shine_dal_[1] = cdsPosBefore - pos_[1];
        }
      break;
    case INSERT :
      position_relative_to_shine_dal_ = new int32_t;
      if (strandBefore == LEADING) position_relative_to_shine_dal_[0] = pos_[0] - cdsPosBefore;
      else                                 position_relative_to_shine_dal_[0] = cdsPosBefore - pos_[0];
      break;
    default :
      fprintf(stderr, "ERROR, invalid mutation type \"%d\" in file %s:%d.\n", mut_type_, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      break;
    }
}


// =================================================================
//                             Destructors
// =================================================================

GeneMutation::~GeneMutation() noexcept
{
  /* Mutation::Mutationwill be called automatically by the compiler for the other attributes */
  switch (mut_type_)
  {
  case SWITCH :
    delete position_relative_to_shine_dal_;
    break;
  case S_INS :
    delete position_relative_to_shine_dal_;
    break;
  case S_DEL :
    delete position_relative_to_shine_dal_;
    break;
  case DUPL :
    delete [] position_relative_to_shine_dal_;
    break;
  case DEL :
    delete [] position_relative_to_shine_dal_;
    break;
  case TRANS :
    delete [] position_relative_to_shine_dal_;
    break;
  case INV :
      delete [] position_relative_to_shine_dal_;
      break;
  case INSERT :
    delete position_relative_to_shine_dal_;
      break;
  default :
    fprintf(stderr, "ERROR, invalid mutation type \"%d\" in file %s:%d.\n" , mut_type_, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    break;
  }

} ;




// =================================================================
//                        Public methods
// =================================================================

// 0 if local mut, 1 if rearrangement, 2 if transfer
int8_t GeneMutation::type_of_event()
{
  if ((mut_type_ == SWITCH) || (mut_type_ == S_INS) || (mut_type_ == S_DEL)) return 0;
  else if ((mut_type_ == DUPL) || (mut_type_ == DEL) || (mut_type_ == TRANS) || (mut_type_ == INV)) return 1;
  else return 2;


}



// str must be at least of size 60
void GeneMutation::description_string_for_gene_mut(char * str)
{
   switch (mut_type_)
  {
    case SWITCH :
    {
      sprintf(str, "%" PRId32 " SWITCH %" PRId32 " %.10f " , generation_, position_relative_to_shine_dal_[0], impact_on_metabolic_error_);
      break;
    }
    case S_INS :
    {
      sprintf(str, "%" PRId32 " SMALL_INS %" PRId32 " %" PRId32 " %s %.10f " , generation_, position_relative_to_shine_dal_[0], length_[0], seq_, impact_on_metabolic_error_);
      break;
    }
    case S_DEL :
    {
      sprintf(str, "%" PRId32 " SMALL_DEL %" PRId32 " %" PRId32 " %.10f " , generation_, position_relative_to_shine_dal_[0], length_[0], impact_on_metabolic_error_);
      break;
    }
    case DUPL :
    {
      sprintf(str, "%" PRId32 " INSERTION_OF_DUPLICATED_DNA %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %.10f " , generation_, position_relative_to_shine_dal_[0], position_relative_to_shine_dal_[1], position_relative_to_shine_dal_[2], length_[0], impact_on_metabolic_error_);
      break;
    }
    case DEL :
    {
      sprintf(str, "%" PRId32 " LARGE_DEL %" PRId32 " %" PRId32 " %" PRId32 " %.10f " , generation_, position_relative_to_shine_dal_[0], position_relative_to_shine_dal_[1], length_[0], impact_on_metabolic_error_);
      break;
    }
    case TRANS :
    {
      sprintf(str, "%" PRId32 " TRANSLOC %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %.10f " , generation_, position_relative_to_shine_dal_[0], position_relative_to_shine_dal_[1],  position_relative_to_shine_dal_[2],  position_relative_to_shine_dal_[3], length_[0], impact_on_metabolic_error_);
      break;
    }
    case INV :
    {
      sprintf(str, "%" PRId32 " INV %" PRId32 " %" PRId32 " %" PRId32 " %.10f " , generation_, position_relative_to_shine_dal_[0], position_relative_to_shine_dal_[1], length_[0], impact_on_metabolic_error_);
      break;
    }
   case INSERT :
    {
      sprintf(str, "%" PRId32 " INSERTION_OF_FOREIGN_DNA %" PRId32 " %" PRId32 " %.10f " , generation_, position_relative_to_shine_dal_[0], length_[0], impact_on_metabolic_error_);
      break;
    }
    default :
    {
      fprintf(stderr, "ERROR, invalid mutation type \"%d\" in file %s:%d.\n" , mut_type_, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      break;
    }
  }

   if (region_ == CDS)            strcat(str, "CDS ");
   else if (region_ == UPSTREAM)  strcat(str, "UPSTREAM ");
   else if (region_ == BOTH)      strcat(str, "BOTH ");


}

} // namespace aevol
