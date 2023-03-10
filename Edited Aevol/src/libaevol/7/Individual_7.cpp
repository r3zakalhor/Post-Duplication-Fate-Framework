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


#include "Individual_7.h"
#include "DnaMutator.h"
#include "ExpManager.h"
#include "Stats_7.h"
// #include "7/DynTab_Metadata.h"
#include "7/List_Metadata.h"
// #include "7/Map_Metadata.h"
#include "7/Protein_7.h"

#include "Fuzzy.h"

#include <algorithm>

namespace aevol {

/** Individual_7 Constructor and Destructor **/
Individual_7::Individual_7(ExpManager* exp_m, double w_max,
                           DnaFactory* dna_factory, FuzzyFactory_7* fuzzy_factory) {
        exp_m_ = exp_m;
        w_max_ = w_max;
        usage_count_ = 1;

        // if (exp_m_->exp_s()->get_simd_metadata_flavor() ==
        //     MetadataFlavor::STD_MAP)
        //     metadata_ = new Map_Metadata(this);
        // else if (exp_m_->exp_s()->get_simd_metadata_flavor() ==
        //          MetadataFlavor::DYN_TAB)
        //     metadata_ = new DynTab_Metadata(this);
        // else if (exp_m_->exp_s()->get_simd_metadata_flavor() ==
        //          MetadataFlavor::STD_LIST)
            metadata_ = new List_Metadata(this);

        dna_factory_ = dna_factory;
        fuzzy_factory_ = fuzzy_factory;
    }

    Individual_7::Individual_7(ExpManager* exp_m,
                               Individual_7* clone,
                               DnaFactory* dna_factory,
                                FuzzyFactory_7* fuzzy_factory,
                               bool no_metadata) {
    w_max_ = clone->w_max_;

  exp_m_ = exp_m;

  usage_count_ = 1;
  dna_ = dna_factory->get_dna(clone->dna_->length());
  //printf("DNA Factory -- %p %p\n",dna_,dna_->data_);
  dna_->set_indiv(clone->dna_,this);

  dna_factory_ = dna_factory;
  fuzzy_factory_ = fuzzy_factory;


  if (no_metadata) {
    // if (exp_m_->exp_s()->get_simd_metadata_flavor() == MetadataFlavor::STD_MAP)
    //     metadata_ = new Map_Metadata(this);
    // else if (exp_m_->exp_s()->get_simd_metadata_flavor() ==
    //         MetadataFlavor::DYN_TAB)
    //     metadata_ = new DynTab_Metadata(this);
    // else if (exp_m_->exp_s()->get_simd_metadata_flavor() ==
    //         MetadataFlavor::STD_LIST)
        metadata_ = new List_Metadata(this);

  } else {
    // if (exp_m_->exp_s()->get_simd_metadata_flavor() == MetadataFlavor::STD_MAP)
    //     metadata_ = new Map_Metadata(this,dynamic_cast<Map_Metadata*>(clone->metadata_));
    // else if (exp_m_->exp_s()->get_simd_metadata_flavor() ==
    //         MetadataFlavor::DYN_TAB)
    //     metadata_ = new DynTab_Metadata(this,dynamic_cast<DynTab_Metadata*>(clone->metadata_));
    // else if (exp_m_->exp_s()->get_simd_metadata_flavor() ==
    //         MetadataFlavor::STD_LIST)
        metadata_ = new List_Metadata(this,dynamic_cast<List_Metadata*>(clone->metadata_));
  }
  fitness = clone->fitness;
  metaerror = clone->metaerror;

}

    Individual_7::Individual_7(ExpManager* exp_m, double w_max,
                               char* dna_clone,
                               int32_t dna_length,
                               DnaFactory* dna_factory,
                                FuzzyFactory_7* fuzzy_factory,
                                int32_t* lead_prom_pos, int8_t* lead_prom_error, int32_t lead_prom_size,
                                int32_t* lag_prom_pos, int8_t* lag_prom_error, int32_t lag_prom_size) {
    w_max_ = w_max;

  exp_m_ = exp_m;

  usage_count_ = 1;
  dna_ = dna_factory->get_dna(dna_length);
  //printf("DNA Factory -- %p %p\n",dna_,dna_->data_);
  dna_->set_indiv(dna_clone,dna_length,this);

  dna_factory_ = dna_factory;
  fuzzy_factory_ = fuzzy_factory;

  if (lead_prom_size == -1 && lag_prom_size == -1)
    metadata_ = new List_Metadata(this);
  else
    metadata_ = new List_Metadata(this,lead_prom_pos,lead_prom_error,lead_prom_size,lag_prom_pos,lag_prom_error,lag_prom_size);

}

Individual_7::~Individual_7() {
  dna_factory_->give_back(dna_);

  if (phenotype!=nullptr) {
    fuzzy_factory_->give_back(phenotype);
    phenotype = nullptr;
  }

  delete metadata_;

  clearAllObserver();
}

/**
 * We need some index for the promoter optimization
 */
void Individual_7::rebuild_index() {
        // if (exp_m_->exp_s()->get_simd_metadata_flavor() == MetadataFlavor::STD_MAP)
        //     dynamic_cast<Map_Metadata*>(metadata_)->rebuild_index();
}


void Individual_7::reset_metadata() {
    delete metadata_;
    // if (exp_m_->exp_s()->get_simd_metadata_flavor() == MetadataFlavor::STD_MAP)
    //     metadata_ = new Map_Metadata(this);
    // else if (exp_m_->exp_s()->get_simd_metadata_flavor() ==
    //         MetadataFlavor::DYN_TAB)
    //     metadata_ = new DynTab_Metadata(this);
    // else if (exp_m_->exp_s()->get_simd_metadata_flavor() ==
    //         MetadataFlavor::STD_LIST)
        metadata_ = new List_Metadata(this);

}


#ifdef BASE_2
void Individual_7::search_start_protein(Rna_7* rna, int32_t pos_1, int32_t pos_2) {
        int32_t dna_length =  dna_->length_;
      char* dna_data = dna_->data_;
      #if defined(__INTEL_COMPILER)
        __declspec(align(dna_data));
      #elif defined(__INTEL_LLVM_COMPILER)
         void* vec_r =  __builtin_assume_aligned(dna_data,64);
      #endif

      int32_t rna_length = rna->length;
      bool rna_leading_lagging = rna->leading_lagging;

  if (rna->is_init_) {
        int32_t s_pos = pos_1;
        if (rna_length >= 21) {


          for (int32_t loop_size = 0; loop_size < pos_2-pos_1; loop_size++) {
            #if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
            __declspec(align(64)) bool start[14] = {false,false,false,false,false,false,false,false,false,false,false,false,false,false};
            #else
            bool start[14] = {false,false,false,false,false,false,false,false,false,false,false,false,false,false};
            #endif

            int32_t c_pos = s_pos;
            if (rna_leading_lagging == 0) {
              c_pos+=loop_size;
              c_pos =
                  c_pos >= dna_length ? c_pos - dna_length
                                      : c_pos;
            } else {
              c_pos-=loop_size;
              c_pos = c_pos < 0 ? dna_length + c_pos : c_pos;
            }


            
            if (rna_leading_lagging == 0) {
              // Search for Shine Dalgarro + START codon on LEADING
              if (c_pos + 15 < dna_length && c_pos + 15 >= 0) {
                #pragma omp simd
                for (int32_t k = 0; k < 12; k++) {
                  // int32_t k_t = k >= 6 ? k + 4 : k;
                 if (dna_data[c_pos+k] == SHINE_DAL_SEQ_LEAD_7[k])
                  start[k] = true;
                }
                if (dna_data[c_pos+12] == SHINE_DAL_SEQ_LEAD_7[12]) start[12] = true;
              } else {
                for (int32_t k = 0; k < 9; k++) {
                  int32_t k_t = k >= 6 ? k + 4 : k;
                  int32_t pos_m = c_pos + k_t;

                  while (pos_m < 0)  pos_m += dna_length;
                  while (pos_m >= dna_length) pos_m -= dna_length;

                  start[k_t] = (dna_data[pos_m] == SHINE_DAL_SEQ_LEAD_7[k_t]) ? true: false;
                }
              }
            } else {
              // Search for Shine Dalgarro + START codon on LAGGING
              if (c_pos - 15 < dna_length && c_pos - 15 >= 0) {
                #pragma omp simd
                for (int32_t k = 0; k < 12; k++) {
                  // int32_t k_t = k >= 6 ? k + 4 : k;
                  if (dna_data[c_pos - k] == SHINE_DAL_SEQ_LAG_7[k])
                    start[k] = true;
                }
                if (dna_data[c_pos - 12] == SHINE_DAL_SEQ_LAG_7[12])
                    start[12] = true;
              } else {
                for (int32_t k = 0; k < 9; k++) {
                  int32_t k_t = k >= 6 ? k + 4 : k;
                  int32_t pos_m = c_pos - k_t;

                  while (pos_m < 0)  pos_m += dna_length;
                  while (pos_m >= dna_length) pos_m -= dna_length;

                  start[k_t] = (dna_data[pos_m] == SHINE_DAL_SEQ_LAG_7[k_t]) ? true : false;
                }
              }
            }

            if (start[0] && start[1] && start[2] && start[3] && start[4] && start[5] && start[10] && start[11] && start[12]) {
              rna->start_prot_count_++;
              rna->start_prot.push_back(c_pos);
            }
          }
        }
      }
      }
#endif

void Individual_7::compute_non_coding() {
  if (non_coding_computed_) return;
  non_coding_computed_ = true;

  // printf("Compute non coding\n");

  // Create a table of <genome_length> bools initialized to false (non-coding)
  int32_t genome_length = dna_->length_;

  // Including Shine-Dalgarno, spacer, START and STOP
  bool* belongs_to_CDS;
  bool* belongs_to_functional_CDS;
  bool* belongs_to_non_functional_CDS; // non-functional CDSs are those that have a null area or that lack a kind of codons (M, W or H)

  // Including Promoters and terminators
  bool* belongs_to_RNA;
  bool* belongs_to_coding_RNA;
  bool* belongs_to_non_coding_RNA;

  // Genes + prom + term (but not UTRs)
  bool* is_essential_DNA;
  bool* is_essential_DNA_including_nf_genes; // Adds non-functional genes + promoters & terminators

  bool* is_not_neutral;            // prom + term + everything in between (as opposed to neutral)


  belongs_to_CDS = new bool[genome_length];
  belongs_to_functional_CDS = new bool[genome_length];
  belongs_to_non_functional_CDS = new bool[genome_length];
  belongs_to_RNA = new bool[genome_length];
  belongs_to_coding_RNA = new bool[genome_length];
  belongs_to_non_coding_RNA = new bool[genome_length];
  is_essential_DNA = new bool[genome_length];
  is_essential_DNA_including_nf_genes = new bool[genome_length];
  is_not_neutral = new bool[genome_length];

  memset(belongs_to_CDS, 0, genome_length);
  memset(belongs_to_functional_CDS, 0, genome_length);
  memset(belongs_to_non_functional_CDS, 0, genome_length);
  memset(belongs_to_RNA, 0, genome_length);
  memset(belongs_to_coding_RNA, 0, genome_length);
  memset(belongs_to_non_coding_RNA, 0, genome_length);
  memset(is_essential_DNA, 0, genome_length);
  memset(is_essential_DNA_including_nf_genes, 0, genome_length);
  memset(is_not_neutral, 0, genome_length);


  // Parse protein lists and mark the corresponding bases as coding
  metadata_->protein_begin();
  for (int protein_idx = 0; protein_idx < metadata_->proteins_count(); protein_idx++) {
    Protein_7* prot = metadata_->protein_next();
    #ifdef __REGUL
    if (!prot->signal_ && prot->is_init_) {
    #else
    if (prot->is_init_) {
    #endif
      int32_t first;
      int32_t last;

      

      switch (prot->leading_lagging) {
        case 0:
          first = prot->protein_start;
          last = prot->protein_end;
          break;
        case 1:
          last = prot->protein_start;
          first = prot->protein_end;
          break;
        default:
          assert(false); // error: should never happen
      }

      if (first <= last) {
        for (int32_t i = first; i <= last; i++) {
          belongs_to_CDS[i] = true;
          if (prot->is_functional) is_essential_DNA[i] = true;
          is_essential_DNA_including_nf_genes[i] = true;
        }
      }
      else {
        for (int32_t i = first; i < genome_length; i++) {
          belongs_to_CDS[i] = true;
          if (prot->is_functional) is_essential_DNA[i] = true;
          is_essential_DNA_including_nf_genes[i] = true;
        }
        for (int32_t i = 0; i <= last; i++) {
          belongs_to_CDS[i] = true;
          if (prot->is_functional) is_essential_DNA[i] = true;
          is_essential_DNA_including_nf_genes[i] = true;
        }
      }

      // Include the promoter and terminator to essential DNA
      // Mark everything between promoter and terminator as not neutral
      for (Rna_7* rna : prot->rna_list_) {

        int32_t prom_first;
        int32_t prom_last;
        int32_t term_first;
        int32_t term_last;
        int32_t rna_first;
        int32_t rna_last;

        if (prot->leading_lagging == 0) { // LEADING
          prom_first = rna->begin;
          prom_last = Utils::mod(prom_first + PROM_SIZE - 1, genome_length);
          term_last = rna->end;
          term_first = Utils::mod(term_last - TERM_SIZE + 1, genome_length);
          rna_first = prom_first;
          rna_last = term_last;
        }
        else {
          prom_last = rna->begin;
          prom_first = Utils::mod(prom_last - PROM_SIZE + 1, genome_length);
          term_first = rna->end;
          term_last = Utils::mod(term_first + TERM_SIZE - 1, genome_length);
          rna_first = term_first;
          rna_last = prom_last;
        }

        // Let us begin with "non-neutral" regions...
        if (rna_first <= rna_last) {
          for (int32_t i = rna_first;
               i <= rna_last; i++) { is_not_neutral[i] = true; }
        }
        else {
          for (int32_t i = rna_first;
               i < genome_length; i++) { is_not_neutral[i] = true; }
          for (int32_t i = 0; i <= rna_last; i++) { is_not_neutral[i] = true; }
        }

        // ...and go on with essential DNA
        if (prom_first <= prom_last) {
          for (int32_t i = prom_first; i <= prom_last; i++) {
            //~ printf("%ld ", i);
            if (prot->is_functional) is_essential_DNA[i] = true;
            is_essential_DNA_including_nf_genes[i] = true;
          }
        }
        else {
          for (int32_t i = prom_first; i < genome_length; i++) {
            //~ printf("%ld ", i);
            if (prot->is_functional) is_essential_DNA[i] = true;
            is_essential_DNA_including_nf_genes[i] = true;
          }
          for (int32_t i = 0; i <= prom_last; i++) {
            //~ printf("%ld ", i);
            if (prot->is_functional) is_essential_DNA[i] = true;
            is_essential_DNA_including_nf_genes[i] = true;
          }
        }
        //~ printf("\n");

        //~ printf("term ");
        if (term_first <= term_last) {
          for (int32_t i = term_first; i <= term_last; i++) {
            //~ printf("%ld ", i);
            if (prot->is_functional) is_essential_DNA[i] = true;
            is_essential_DNA_including_nf_genes[i] = true;
          }
        }
        else {
          for (int32_t i = term_first; i < genome_length; i++) {
            //~ printf("%ld ", i);
            if (prot->is_functional) is_essential_DNA[i] = true;
            is_essential_DNA_including_nf_genes[i] = true;
          }
          for (int32_t i = 0; i <= term_last; i++) {
            //~ printf("%ld ", i);
            if (prot->is_functional) is_essential_DNA[i] = true;
            is_essential_DNA_including_nf_genes[i] = true;
          }
        }
        //~ printf("\n");
        //~ getchar();
      }


      if (prot->is_functional) {
        if (first <= last) {
          for (int32_t i = first; i <= last; i++) {
            belongs_to_functional_CDS[i] = true;
          }
        }
        else {
          for (int32_t i = first; i < genome_length; i++) {
            belongs_to_functional_CDS[i] = true;
          }
          for (int32_t i = 0; i <= last; i++) {
            belongs_to_functional_CDS[i] = true;
          }
        }
      }
      else // degenerated protein
      {
        if (first <= last) {
          for (int32_t i = first; i <= last; i++) {
            belongs_to_non_functional_CDS[i] = true;
          }
        }
        else {
          for (int32_t i = first; i < genome_length; i++) {
            belongs_to_non_functional_CDS[i] = true;
          }
          for (int32_t i = 0; i <= last; i++) {
            belongs_to_non_functional_CDS[i] = true;
          }
        }
      }
    }
  }


  // Parse RNA lists and mark the corresponding bases as coding (only for the coding RNAs)
  // TODO vld: this block cries for refactoring
  metadata_->rna_begin();

  for (int rna_idx = 0; rna_idx <
                        (int)metadata_->rna_count(); rna_idx++) {
      Rna_7* rna =
          metadata_->rna_next();
      int32_t first;
      int32_t last;

      if (rna->leading_lagging == 0) {
        first = rna->begin;
        last = rna->end;
      }
      else { // (strand == LAGGING)
        first = rna->end;
        last = rna->begin;
      }

      assert(first < genome_length);
      assert(last < genome_length);

      if (first <= last) {
        for (int32_t i = first; i <= last; i++) {
          belongs_to_RNA[i] = true;
        }
      }
      else {
        for (int32_t i = first; i < genome_length; i++)
          belongs_to_RNA[i] = true;
        for (int32_t i = 0; i <= last; i++)
          belongs_to_RNA[i] = true;
      }

      bool is_coding_rna = false;
      metadata_->protein_begin();
      for (int protein_idx = 0; protein_idx < metadata_->proteins_count(); protein_idx++) {
        Protein_7* prot = metadata_->protein_next();
        if (prot->is_init_) {
          for (auto i_rna : prot->rna_list_) {
            if (i_rna->begin == rna->begin && i_rna->end == rna->end && i_rna->leading_lagging == rna->leading_lagging) {
              is_coding_rna = true;
            }

          }
        }
      }


      if (is_coding_rna) { // coding RNA
        if (first <= last) {
          for (int32_t i = first; i <= last; i++)
            belongs_to_coding_RNA[i] = true;
        }
        else {
          for (int32_t i = first; i < genome_length; i++)
            belongs_to_coding_RNA[i] = true;
          for (int32_t i = 0; i <= last; i++)
            belongs_to_coding_RNA[i] = true;
        }
      }
      else // non coding RNA
      {
        if (first <= last) {
          for (int32_t i = first; i <= last; i++)
            belongs_to_non_coding_RNA[i] = true;
        }
        else {
          for (int32_t i = first; i < genome_length; i++)
            belongs_to_non_coding_RNA[i] = true;
          for (int32_t i = 0; i <= last; i++)
            belongs_to_non_coding_RNA[i] = true;
        }
      }
    }
  

  // Count non-coding bases
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

  // We do not know how many neutral regions there will be, but
  // there should be less than nb_coding_RNAs_ + 1
  // As we will see, there may be a shift in values so we take size nb_coding_RNAs_ + 2
  int32_t* tmp_beginning_neutral_regions = new int32_t[nb_coding_RNAs + 2];
  int32_t* tmp_end_neutral_regions = new int32_t[nb_coding_RNAs + 2];
  memset(tmp_beginning_neutral_regions, -1, nb_coding_RNAs + 2);
  memset(tmp_end_neutral_regions, -1, nb_coding_RNAs + 2);

  for (int32_t i = 0; i < genome_length; i++) {
    if (belongs_to_CDS[i] == false) {
      nb_bases_in_0_CDS_++;
    }
    if (belongs_to_functional_CDS[i] == false) {
      nb_bases_in_0_functional_CDS_++;
    }
    if (belongs_to_non_functional_CDS[i] == false) {
      nb_bases_in_0_non_functional_CDS_++;
    }
    if (belongs_to_RNA[i] == false) {
      nb_bases_in_0_RNA_++;
    }
    if (belongs_to_coding_RNA[i] == false) {
      nb_bases_in_0_coding_RNA_++;
    }
    if (belongs_to_non_coding_RNA[i] == false) {
      nb_bases_in_0_non_coding_RNA_++;
    }
    if (is_essential_DNA[i] == false) {
      nb_bases_non_essential_++;
    }
    if (is_essential_DNA_including_nf_genes[i] == false) {
      nb_bases_non_essential_including_nf_genes_++;
    }
    if (is_not_neutral[i] == false) {
      nb_bases_in_neutral_regions_++;
    }
    if (i != 0) {
      if (is_not_neutral[i] != is_not_neutral[i - 1]) {
        if (is_not_neutral[i - 1] == true) // beginning of a neutral region
        {
          tmp_beginning_neutral_regions[nb_neutral_regions_] = i;
        }
        else // end of a neutral region
        {
          tmp_end_neutral_regions[nb_neutral_regions_] = i - 1;
          nb_neutral_regions_++;
        }
      }
    }
    else // i = 0
    {
      // we arbitrarily set 0 as the beginning of a neutral region (linkage with end of genetic unit
      // will be done later)
      if (is_not_neutral[0] ==
          false) { tmp_beginning_neutral_regions[nb_neutral_regions_] = 0; }
    }
  }

  // we have to treat specifically the last base of the genetic unit in order to link neutral regions
  // at the end and the beginning of genetic unit
  int32_t shift = 0;
  if (is_not_neutral[genome_length - 1] == false) {
    if (is_not_neutral[0] == true) // end of a neutral region
    {
      tmp_end_neutral_regions[nb_neutral_regions_] = genome_length - 1;
      nb_neutral_regions_++;
    }
    else // neutral region goes on after base 0, linkage to be done
    {
      if (nb_neutral_regions_ != 0) {
        tmp_end_neutral_regions[nb_neutral_regions_] = tmp_end_neutral_regions[0];
        // the first neutral region is only a subpart of the last one, it should not be
        // taken into account. When we transfer values to the final array, we introduce a shift
        shift = 1;
        // we do not ++ nb_neutral_regions_ as it was already counted
      }
      else // no neutral region detected until now -> all the genetic unit is neutral
      {
        // as all the chromosome is neutral, we indicate 0 as the beginning of the region
        // and genome_length - 1 as its end
        tmp_end_neutral_regions[0] = genome_length - 1;
        nb_neutral_regions_++;
      }
    }
  }

  // now that we know how many neutral regions there are, we can transfer data to correctly sized arrays
  assert(nb_neutral_regions_ <= nb_coding_RNAs + 1);
  if (beginning_neutral_regions_ !=
      NULL) { delete[] beginning_neutral_regions_; }
  if (end_neutral_regions_ != NULL) { delete[] end_neutral_regions_; }

  if (nb_neutral_regions_ >
      0) // as unlikely as it seems, there may be no neutral region
  {
    beginning_neutral_regions_ = new int32_t[nb_neutral_regions_];
    end_neutral_regions_ = new int32_t[nb_neutral_regions_];
    // transfer from tmp to attributes
    for (int32_t i = 0; i < nb_neutral_regions_; i++) {
      beginning_neutral_regions_[i] = tmp_beginning_neutral_regions[i + shift];
      end_neutral_regions_[i] = tmp_end_neutral_regions[i + shift];
    }
  }
  else // nb_neutral_regions_ == 0
  {
    beginning_neutral_regions_ = NULL;
    end_neutral_regions_ = NULL;
  }

  // printf("nb_bases_in_0_CDS_ %d\n",nb_bases_in_0_CDS_);

  delete[] tmp_beginning_neutral_regions;
  delete[] tmp_end_neutral_regions;

  delete[] belongs_to_CDS;
  delete[] belongs_to_functional_CDS;
  delete[] belongs_to_non_functional_CDS;
  delete[] belongs_to_RNA;
  delete[] belongs_to_coding_RNA;
  delete[] belongs_to_non_coding_RNA;
  delete[] is_essential_DNA;
  delete[] is_essential_DNA_including_nf_genes;
  delete[] is_not_neutral;
}
}
