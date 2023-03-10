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


#ifndef AEVOL_INDIVIDUAL_7_H
#define AEVOL_INDIVIDUAL_7_H

#include <vector>
#include <map>
#include <set>

#include "DnaFactory.h"
#include "Observable.h"
#include "PhenotypicTargetHandler.h"
#include "Stats.h"
#include "Vector_Fuzzy.h"
#include "FuzzyFactory_7.h"
#include "ae_enums.h"

namespace aevol {

class ExpManager;
class Dna_7;
class Abstract_Metadata;
class Map_Metadata;
class Rna_7;

class Individual_7 : public Observable {
 public:
  Individual_7(ExpManager* exp_m, double w_max, DnaFactory* dna_factory, FuzzyFactory_7* fuzzy_factory);

  Individual_7(ExpManager* exp_m, Individual_7* clone, DnaFactory* dna_factory,  FuzzyFactory_7* fuzzy_factory, bool no_metadata = false);

  Individual_7(ExpManager* exp_m, double w_max,
                               char* dna_clone,
                               int32_t dna_length,
                               DnaFactory* dna_factory,
                                FuzzyFactory_7* fuzzy_factory,
                                int32_t* lead_prom_pos = nullptr, int8_t* lead_prom_error = nullptr, int32_t lead_prom_size = -1,
                                int32_t* lag_prom_pos = nullptr, int8_t* lag_prom_error = nullptr, int32_t lag_prom_size = -1);
  ~Individual_7();


  void search_start_protein(Rna_7* rna, int32_t pos_1, int32_t pos_2);

  void compute_non_coding();

#ifdef PHENOTYPE_VECTOR
  double phenotype[PHENOTYPE_VECTOR_SIZE];
  double delta[PHENOTYPE_VECTOR_SIZE];
#else
  AbstractFuzzy_7* phenotype = nullptr;
  AbstractFuzzy_7* delta = nullptr;
#endif
  double fitness;
  double metaerror;


  double* fitness_by_env_id_;
  double* metaerror_by_env_id_;

  Dna_7* dna_ = nullptr;

  int32_t indiv_id;
  int32_t parent_id;
  int32_t last_id;
  bool first_to_add = true;
  int32_t added_id = -1;

  int32_t usage_count_ = 1;

  ExpManager* exp_m_;
  DnaFactory* dna_factory_;
  FuzzyFactory_7* fuzzy_factory_;
  int global_id = -1;

  double w_max_;

  /** Variables for Tree mgmt **/
  int32_t nb_genes_activ     = 0;
  int32_t nb_genes_inhib     = 0;
  int32_t nb_func_genes      = 0;
  int32_t nb_non_func_genes  = 0;
  int32_t nb_coding_RNAs     = 0;
  int32_t nb_non_coding_RNAs = 0;
  /** END of Variables for Tree Mgmt **/

/** Non-coding variables **/
  double overall_size_coding_RNAs_;
  double overall_size_non_coding_RNAs_;

  double overall_size_fun_genes_;
  double overall_size_non_fun_genes_;

  int32_t nb_bases_in_0_CDS_;
  int32_t nb_bases_in_0_functional_CDS_;
  int32_t nb_bases_in_0_non_functional_CDS_;
  int32_t nb_bases_in_0_RNA_;
  int32_t nb_bases_in_0_coding_RNA_;
  int32_t nb_bases_in_0_non_coding_RNA_;
  int32_t nb_bases_non_essential_;
  int32_t nb_bases_non_essential_including_nf_genes_;
  int32_t nb_bases_in_neutral_regions_;
  int32_t nb_neutral_regions_;
  int32_t* beginning_neutral_regions_ = nullptr;
  int32_t* end_neutral_regions_ = nullptr;

  void reset_stats() {
    nb_genes_activ     = 0;
    nb_genes_inhib     = 0;
    nb_func_genes      = 0;
    nb_non_func_genes  = 0;
    nb_coding_RNAs     = 0;
    nb_non_coding_RNAs = 0;
  
    overall_size_coding_RNAs_ = 0;
    overall_size_non_coding_RNAs_ = 0;

    overall_size_fun_genes_ = 0;
    overall_size_non_fun_genes_ = 0;

    nb_bases_in_0_CDS_ = 0;
    nb_bases_in_0_functional_CDS_ = 0;
    nb_bases_in_0_non_functional_CDS_ = 0;
    nb_bases_in_0_RNA_ = 0;
    nb_bases_in_0_coding_RNA_ = 0;
    nb_bases_in_0_non_coding_RNA_ = 0;
    nb_bases_non_essential_ = 0;
    nb_bases_non_essential_including_nf_genes_ = 0;
    nb_bases_in_neutral_regions_ = 0;
    nb_neutral_regions_ = 0;
  }

  void rebuild_index();

  void reset_metadata();

  Abstract_Metadata* metadata_;

  bool is_at_border_ = false;

  bool non_coding_computed_ = false;
  
};



}
#endif
