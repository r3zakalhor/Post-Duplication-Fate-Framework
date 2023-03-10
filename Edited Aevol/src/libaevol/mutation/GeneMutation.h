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

/*! \class GeneMutation
    \brief Currently used only by post-treatments, on a specific lineage, to monitor the fate of paralogs.
         Each paralog maintains a list of the mutations it underwent. A gene mutation is a mutation, but
         enriched with the generation when it occurred and the position where it occurred in the coding RNA
         (relative to the first bp of the promoter).
*/


 #ifndef AEVOL_GENE_MUTATION_H_
#define  AEVOL_GENE_MUTATION_H_


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>




// =================================================================
//                            Project Files
// =================================================================
#include "Mutation.h"
#include "Rna.h"
#include "ae_enums.h"


namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================



enum ae_gene_loss_type
{
  NOT_LOST_YET = 0,
  LOST_BY_LOCAL_MUTATION  = 1,
  DELETED = 2,
  BROKEN_BY_REAR = 3,
  DUPLICATED = 4
};




enum ae_gene_mutation_region
{
  UPSTREAM = 0,
  CDS = 1,
  BOTH = 2,
};





class GeneMutation : public Mutation
{
  friend class GeneTreeNode;

 public :

  // =================================================================
  //                             Constructors
  // =================================================================

  // Creates a copy of the mutation mut, but enriched with the generation when it occured
  // and the position where it occurred in the RNA, relative to the first bp of the promoter
  GeneMutation(Mutation const & mut, int32_t gener, int32_t cdsPosBefore, Strand strandBefore, ae_gene_mutation_region region);

  GeneMutation(const GeneMutation &model);

  // =================================================================
  //                             Destructors
  // =================================================================


  virtual ~GeneMutation();
  // =================================================================
  //                              Accessors
  // =================================================================

  inline int32_t generation() const;
  inline double impact_on_metabolic_error() const;
  inline ae_gene_mutation_region region();
  inline void set_impact_on_metabolic_error(double impact);


  // =================================================================
  //                            Public Methods
  // =================================================================
  void description_string_for_gene_mut(char * str); // str must be at least of size 60
  int8_t type_of_event(); // 0 if local mut, 1 if rearrangement, 2 if transfer


 protected :

  // =================================================================
  //                         Forbidden Constructors
  // =================================================================

  GeneMutation()
    {
      printf("ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    };
  /* GeneMutation(const GeneMutation &model) */
  /*   { */
  /*     printf("ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__); */
  /*     exit(EXIT_FAILURE); */
  /*   }; */

  // =================================================================
  //                           Protected Methods
  // =================================================================

  // =================================================================
  //                          Protected Attributes
  // =================================================================

  int32_t*  position_relative_to_shine_dal_; /* array of positions similar to the pos_ array of the Mutation class (size 1 for the switch, 2 for an inversion, etc.) */
  int32_t   generation_;
  double    impact_on_metabolic_error_;
  ae_gene_mutation_region region_;

};





// =====================================================================
//                         Inline Accessors' definitions
// =====================================================================

inline int32_t GeneMutation::generation() const
{
  return generation_;
}

inline double GeneMutation::impact_on_metabolic_error() const
{
  return impact_on_metabolic_error_;

}


inline void GeneMutation::set_impact_on_metabolic_error(double impact)
{
  impact_on_metabolic_error_ = impact;
}


inline ae_gene_mutation_region GeneMutation::region()
{
  return region_;
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol
#endif // AEVOL_GENE_MUTATION_H_
