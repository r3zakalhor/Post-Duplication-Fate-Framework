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


#ifndef AEVOL_STATS_7_H
#define AEVOL_STATS_7_H


#include "ExpManager_7.h"
#include <cstdint>
#include <fstream>
#include <limits>

namespace aevol {

    class ExpManager_7;

class Stats_7 {
 public:
  Stats_7(ExpManager_7* simd_individual, int64_t generation,
               bool best_or_not, char* prefix = nullptr
               #ifdef HAVE_MPI
               , int32_t rank = 1
               #endif
               ,bool non_coding = false
  );

    ~Stats_7() {
      if (is_indiv_) {
        statfile_best_.flush();
        statfile_best_.close();
      } else {
        statfile_mean_.flush();
        statfile_mean_.close();
      }
    }

    void compute(Individual_7* indiv,bool non_coding = false);
    void compute_best(bool non_coding = false);
    void compute_average(bool non_coding = false);

    void write(Individual_7* indiv,bool non_coding = false);
    void write_best(bool non_coding = false);
    void write_average(bool non_coding = false);

    void reinit(int64_t generation);

    bool is_indiv() { return is_indiv_; }


 protected:
  ExpManager_7* simd_individual_;

    int64_t generation_;

    bool is_indiv_ = true;

    int32_t pop_size_ = 0;
    int32_t nb_clones_ = 0;

    double fitness_ = 0;
    double metabolic_error_ = 0;

    int32_t amount_of_dna_ = 0;
    int32_t nb_coding_rnas_ = 0;
    int32_t nb_non_coding_rnas_ = 0;

    int32_t nb_functional_genes_ = 0;
    int32_t nb_non_functional_genes_ = 0;

    int32_t nb_mut_ = 0;
    int32_t nb_rear_ = 0;
    int32_t nb_switch_ = 0;
    int32_t nb_indels_ = 0;
    int32_t nb_dupl_ = 0;
    int32_t nb_del_ = 0;
    int32_t nb_trans_ = 0;
    int32_t nb_inv_ = 0;

    double dupl_rate_ = 0;
    double del_rate_ = 0;
    double trans_rate_ = 0;
    double inv_rate_ = 0;

    int32_t nb_bases_in_0_CDS_ = 0;
    int32_t nb_bases_in_0_functional_CDS_ = 0;
    int32_t nb_bases_in_0_non_functional_CDS_ = 0;
    int32_t nb_bases_in_0_RNA_ = 0;
    int32_t nb_bases_in_0_coding_RNA_ = 0;
    int32_t nb_bases_in_0_non_coding_RNA_ = 0;

    int32_t nb_bases_non_essential_ = 0;
    int32_t nb_bases_non_essential_including_nf_genes_ = 0;

#ifdef __REGUL
    int32_t  nb_influences_ = 0.0;
    int32_t  nb_enhancing_influences_ = 0.0;
    int32_t  nb_operating_influences_ = 0.0;
    double  av_value_influences_ = 0.0;
    double  av_value_enhancing_influences_ = 0.0;
    double  av_value_operating_influences_ = 0.0;
    int32_t  nb_TF_ = 0;
    int32_t  nb_pure_TF_ = 0;
#endif

    bool is_computed_ = false;

    // Stats

    std::ofstream statfile_best_;
    std::ofstream statfile_mean_;

};

}
#endif //RAEVOL_CUDA_STATS_SIMD_H
