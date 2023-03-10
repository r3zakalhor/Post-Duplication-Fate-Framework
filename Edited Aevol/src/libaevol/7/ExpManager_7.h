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


#ifndef AEVOL_EXPMANAGER_7_H
#define AEVOL_EXPMANAGER_7_H

#include "Stats.h"
#include "DnaFactory.h"
#include "Observable.h"
#include "PhenotypicTargetHandler.h"
#include "raevol/SIMD_PhenotypicTargetHandler_R.h"
#include "FuzzyFactory_7.h"

namespace aevol {

#ifdef BASE_2
constexpr const char* PROM_SEQ_LEAD = "0101011001110010010110";
constexpr const char* PROM_SEQ_LAG  = "1010100110001101101001";

#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
__declspec(align(64)) constexpr const char* SHINE_DAL_SEQ_LEAD_7 = "0110111111000";
__declspec(align(64)) constexpr const char* SHINE_DAL_SEQ_LAG_7  = "1001001111111";
#else
constexpr const char* SHINE_DAL_SEQ_LEAD_7 = "0110111111000";
constexpr const char* SHINE_DAL_SEQ_LAG_7  = "1001001111111";
#endif

// constexpr const char* SHINE_DAL_SEQ_LEAD = "011011000";
// constexpr const char* SHINE_DAL_SEQ_LAG  = "100100111";

constexpr const char* PROTEIN_END_LEAD = "001";
constexpr const char* PROTEIN_END_LAG  = "110";
#endif

class ToFetchIndividual {
  public:
    ToFetchIndividual(int32_t id, int32_t x, int32_t y, int32_t rank) 
    {
      id_ = id; x_ = x; y_ = y; rank_ = rank; 
      // printf("---> Add FETCHIndividual %d (%d %d)\n",id_,x_,y_);
    
    }

    int32_t x_;
    int32_t y_;
    int32_t id_;
    int32_t rank_;


    int32_t local_x_;
    int32_t local_y_;
    
    int32_t mpi_id_;

    char* dna_ = nullptr;
    int32_t dna_length_ = -1;

    int32_t lead_prom_size;

    int32_t* lead_prom_pos = nullptr;
    int8_t* lead_prom_error = nullptr;
            
    int32_t lag_prom_size;
            
    int32_t* lag_prom_pos = nullptr;
    int8_t* lag_prom_error = nullptr;


    bool fetched_ = false;
};

class ExpManager;
class Stats_7;
// class FuzzyFactory_7;

class ExpManager_7 : public Observable{
 public:
  ExpManager_7(ExpManager* exp_m);

  ~ExpManager_7();

  void setup_individuals(double w_max, double selection_pressure);
  void run_a_step(double w_max, double selection_pressure);

// Function with ID
  void apply_mutations(int indiv_id, double w_max = -1.0, double selection_pressure = -1.0);

  void do_mutation(int indiv_id, double w_max = -1.0, double selection_pressure = -1.0);
  void start_stop_RNA(int indiv_id);
  void opt_prom_compute_RNA(int indiv_id);
  void compute_RNA(int indiv_id);
  void start_protein(int indiv_id);
  void compute_protein(int indiv_id);
  void translate_protein(int indiv_id, double w_max);
  void compute_phenotype(int indiv_id);
  void compute_fitness(int indiv_id, double selection_pressure, int env_id = -1);
#ifdef __REGUL
  void compute_network(int indiv_id, double selection_pressure);
  void update_network(int indiv_id, double selection_pressure);
  void evaluate_network(int indiv_id, double selection_pressure, int env_id);
  void finalize_network(int indiv_id, double selection_pressure);
  void solve_network(int indiv_id, double selection_pressure);
  void update_phenotype( int indiv_id );
#endif

  void write_stat(bool just_best = true, bool non_coding = false);
  void write_stat(Stats_7* stats, Individual_7* indiv, int32_t generation = -1, bool non_coding = false);
  void check_result();
  void check_dna();
  void check_struct();
  void check_individual(int indiv_id, int x, int y);
  static bool standalone() { return standalone_simd; }


// Function with Pointer
  void evaluate(Individual_7* indiv, double w_max, double selection_pressure, SIMD_PhenotypicTargetHandler_R* pth = nullptr);
  void do_mutation(Individual_7* indiv, int indiv_id);
  void start_stop_RNA(Individual_7* indiv);
  void opt_prom_compute_RNA(Individual_7* indiv);
  void compute_RNA(Individual_7* indiv);
  void start_protein(Individual_7* indiv);
  void compute_protein(Individual_7* indiv);
  void translate_protein(Individual_7* indiv, double w_max);
  void compute_phenotype(Individual_7* indiv);
  void compute_fitness(Individual_7* indiv, double selection_pressure, int env_id = -1
  #ifdef __REGUL
    , SIMD_PhenotypicTargetHandler_R* pth = nullptr
  #endif
  );
#ifdef __REGUL
  void compute_network(Individual_7* indiv, double selection_pressure, SIMD_PhenotypicTargetHandler_R* pth = nullptr);
  void update_network(Individual_7* indiv, double selection_pressure, int degradation_step = -1,int lifestep_step = -1);
  void evaluate_network(Individual_7* indiv, double selection_pressure, int env_id, SIMD_PhenotypicTargetHandler_R* pth = nullptr);
  void finalize_network(Individual_7* indiv, double selection_pressure, SIMD_PhenotypicTargetHandler_R* pth = nullptr);
  void solve_network(Individual_7* indiv, double selection_pressure, SIMD_PhenotypicTargetHandler_R* pth = nullptr, bool verbose = true);
  void recompute_network(Individual_7* indiv, double selection_pressure, SIMD_PhenotypicTargetHandler_R* pth);
  void update_phenotype(Individual_7* indiv);
#endif

#ifdef HAVE_MPI
void send_broadcast_border();
void recv_broadcast_border();
void Fetch_Remote_Individual();
#endif

  Individual_7** current_individuals;
  Individual_7** previous_individuals;
  Individual_7* best_indiv;

  std::vector<int> mutant_list_;

  #ifdef HAVE_MPI
  double** fitness_at_border_;
  std::list<ToFetchIndividual> individual_to_fetch_at_border_;
  std::list<ToFetchIndividual> individual_sent_at_border_;

  std::list<Individual_7*> to_delete_individuals_;
  #endif

  int32_t nb_indivs_;
  int32_t nb_clones_;

  static bool standalone_simd;//= true;
  int rna_grain_size = 32;
  int protein_grain_size = 32;

  // static bool compute_diff_rnas;

  DnaFactory* dna_factory_;
  FuzzyFactory_7* fuzzy_factory_;

  double* fitness_sum_tab_;
#ifdef __REGUL
  SIMD_PhenotypicTargetHandler_R* phenotypic_target_handler_;
#else
  AbstractFuzzy_7* target;
#endif

  int32_t grid_width_;
  int32_t grid_height_;

#ifdef HAVE_MPI
  int32_t global_grid_width_;
  int32_t global_grid_height_;

  int32_t local_rank_x_;
  int32_t local_rank_y_;

  int32_t x_offset_;
  int32_t y_offset_;

  int32_t rank_x_;
  int32_t rank_y_;
  
  int32_t localXtoGlobalX(int32_t lx);
  int32_t globalXtoLocalX(int32_t gx);

  int32_t localYtoGlobalY(int32_t ly);
  int32_t globalYtoLocalY(int32_t gy);

  int32_t xOffset();
  int32_t yOffset();
  int32_t rankOf(int32_t x, int32_t y);

  void setRank(int32_t rank);
#endif
 private:
  ExpManager* exp_m_;

  Stats_7* stats_best = nullptr;
  Stats_7* stats_mean = nullptr;

  void selection(int indiv_id);

  void check_selection(int indiv_id);

#ifdef WITH_PERF_TRACES_PER_INDIV
  long* apply_mutation;
  long* compute_rna_;
  long* start_protein_;
  long* compute_protein_;
  long* translate_protein_;
  long* compute_phenotype_;
  long* compute_fitness_;
  long* total_;

  long* allocate_individual_start_;
  long* apply_mutation_start_;
  long* compute_rna_start_;
  long* start_protein_start_;
  long* compute_protein_start_;
  long* translate_protein_start_;
  long* compute_phenotype_start_;
  long* compute_fitness_start_;
  
  long* allocate_individual_stop_;
  long* apply_mutation_stop_;
  long* compute_rna_stop_;
  long* start_protein_stop_;
  long* compute_protein_stop_;
  long* translate_protein_stop_;
  long* compute_phenotype_stop_;
  long* compute_fitness_stop_;

  long* total_start_;
  long* total_stop_;
  long* omp_tid_;
#endif


};
}

#endif //AEVOL_EXPMANAGER_7_H
