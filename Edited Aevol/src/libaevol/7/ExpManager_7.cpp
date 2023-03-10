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


#include "ExpManager_7.h"
#include "AbstractFuzzy_7.h"
#include "7/Stats_7.h"
#include "Abstract_Metadata.h"
#include "DnaMutator.h"
#include "ExpManager.h"
#include "Promoter.h"
#include "Protein_7.h"
#include "Rna_7.h"
#include "List_Metadata.h"
#include "FuzzyFactory_7.h"

#include <algorithm>
#include <err.h>
#include <fstream>
#include <sys/stat.h>
#include<chrono>
#include <omp.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef __OMP_LIST_SORT
  #if 1 == __OMP_LIST_SORT
    #include <random>       // std::default_random_engine
  #endif
#endif

//#include "/opt/intel/advisor_2020/include/advisor-annotate.h"
namespace aevol {

#ifndef WITH_STANDALONE_SIMD
bool ExpManager_7::standalone_simd = false;
#else
bool ExpManager_7::standalone_simd = true;
#endif

#define __VECTORIZE_STRCMP

// #ifdef WITH_OPTIMIZE_DIFF_SEARCH
// bool ExpManager_7::compute_diff_rnas = true;
// #else
// bool ExpManager_7::compute_diff_rnas = false;
// #endif

ExpManager_7::ExpManager_7(ExpManager* exp_m) {
  printf("  Loading SIMD Controller...");

  exp_m_ = exp_m;

  nb_indivs_ = exp_m_->nb_indivs();

  current_individuals       = new Individual_7*[exp_m_->nb_indivs()];
  previous_individuals      = new Individual_7*[exp_m_->nb_indivs()];

  // Allocate Dna_7
  int max_size_dna = -1;
  for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++) {
    int x = indiv_id / exp_m_->world()->height();
    int y = indiv_id % exp_m_->world()->height();

    max_size_dna = max_size_dna < exp_m->world()->grid(x, y)->individual()->genetic_unit_seq_length(0) ?
                   exp_m->world()->grid(x, y)->individual()->genetic_unit_seq_length(0) : max_size_dna;
  }

  // printf("NDA Factory\n");
  #if 0 == __DNA_FACTORY_ALG
  dna_factory_ = new DnaFactory(DnaFactory_Policy::FIRST,exp_m_->nb_indivs()*3,max_size_dna,nb_indivs_);
  #elif 1 == __DNA_FACTORY_ALG
  dna_factory_ = new DnaFactory(DnaFactory_Policy::FIRSTFIT,exp_m_->nb_indivs()*3,max_size_dna,nb_indivs_);
  #elif 2 == __DNA_FACTORY_ALG
  dna_factory_ = new DnaFactory(DnaFactory_Policy::LOCAL_GLOBAL_FIT,exp_m_->nb_indivs()*3,max_size_dna,nb_indivs_);
  #elif 3 == __DNA_FACTORY_ALG
  dna_factory_ = new DnaFactory(DnaFactory_Policy::ALLOCATE,exp_m_->nb_indivs()*3,max_size_dna,nb_indivs_);
  #endif
  
  fuzzy_factory_ = new FuzzyFactory_7(exp_m_->exp_s()->get_fuzzy_flavor(),exp_m_->nb_indivs()*4,
                        exp_m->world()->phenotypic_target_handler()->sampling(),nb_indivs_);

  for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++) {
    int x = indiv_id / exp_m_->world()->height();
    int y = indiv_id % exp_m_->world()->height();

    current_individuals[indiv_id] = new Individual_7(exp_m_, exp_m_->world()->indiv_at(0,0)->w_max(),dna_factory_,fuzzy_factory_);
    current_individuals[indiv_id]->dna_ = dna_factory_->get_dna(exp_m->world()->grid(x, y)->individual()->genetic_unit_seq_length(0));
    current_individuals[indiv_id]->dna_->set_indiv(exp_m->world()->grid(x, y)->individual()->genetic_unit(0).dna(),dna_factory_);
    current_individuals[indiv_id]->indiv_id = indiv_id;
    current_individuals[indiv_id]->parent_id = indiv_id;
    previous_individuals[indiv_id] = current_individuals[indiv_id];
    current_individuals[indiv_id]->global_id = AeTime::time() * 1024 + indiv_id;
  }

#ifdef __REGUL
  printf("VAR Method %d\n",exp_m->world()->phenotypic_target_handler()->var_method());
  phenotypic_target_handler_ = new SIMD_PhenotypicTargetHandler_R(
        std::dynamic_pointer_cast<PhenotypicTargetHandler_R>(exp_m->world()->phenotypic_target_handler()),
        exp_m->exp_s(),fuzzy_factory_,exp_m->check_simd());
#else
  target = fuzzy_factory_->get_fuzzy();
  target->copy((Fuzzy*)(exp_m->world()->phenotypic_target_handler()->phenotypic_target().fuzzy()));
#endif

#ifdef WITH_PERF_TRACES
  std::ofstream perf_traces_file_;

  #ifdef HAVE_MPI
    char filename[256];
    sprintf(filename,"mpi_perf_traces_%d.csv",exp_m_->rank());
    perf_traces_file_.open(filename, std::ofstream::trunc);
    perf_traces_file_ << "generation,runtime,io,exch_fitness,selection,exch_indiv,evaluate" << std::endl;
    perf_traces_file_.close();
  #else
        perf_traces_file_.open("sm_perf_traces.csv", std::ofstream::trunc);
        perf_traces_file_ << "generation,runtime,io" << std::endl;
        perf_traces_file_.close();
  #endif
#endif


  grid_width_ = exp_m_->grid_width();
  grid_height_ = exp_m_->grid_height();

#ifdef HAVE_MPI
fitness_at_border_ = new double*[4];
fitness_at_border_[TOP] = new double[exp_m_->world()->width()+2];
fitness_at_border_[DOWN] = new double[exp_m_->world()->width()+2];
fitness_at_border_[LEFT] = new double[exp_m_->world()->height()+2];
fitness_at_border_[RIGHT] = new double[exp_m_->world()->height()+2];


  global_grid_width_ = exp_m_->exp_s()->global_grid_width();
  global_grid_height_ = exp_m_->exp_s()->global_grid_height();


  rank_x_ = exp_m_->exp_s()->rank_width();
  rank_y_ = exp_m_->exp_s()->rank_height();

  local_rank_x_ = exp_m_->rank() / rank_y_;
  local_rank_y_ = exp_m_->rank() % rank_y_;

  x_offset_ = (local_rank_x_) * grid_width_; 
  y_offset_ = (local_rank_y_) * grid_height_;
#endif

  printf(" OK\n");
  #ifdef WITH_PERF_TRACES_PER_INDIV
  apply_mutation = new long[exp_m_->nb_indivs()];
  compute_rna_ = new long[exp_m_->nb_indivs()];
  start_protein_ = new long[exp_m_->nb_indivs()];
  compute_protein_ = new long[exp_m_->nb_indivs()];
  translate_protein_ = new long[exp_m_->nb_indivs()];
  compute_phenotype_ = new long[exp_m_->nb_indivs()];
  compute_fitness_ = new long[exp_m_->nb_indivs()];
  total_ = new long[exp_m_->nb_indivs()];


  allocate_individual_start_ = new long[exp_m_->nb_indivs()];
  apply_mutation_start_ = new long[exp_m_->nb_indivs()];
  compute_rna_start_ = new long[exp_m_->nb_indivs()];
  start_protein_start_ = new long[exp_m_->nb_indivs()];
  compute_protein_start_ = new long[exp_m_->nb_indivs()];
  translate_protein_start_ = new long[exp_m_->nb_indivs()];
  compute_phenotype_start_ = new long[exp_m_->nb_indivs()];
  compute_fitness_start_ = new long[exp_m_->nb_indivs()];
  
  allocate_individual_stop_ = new long[exp_m_->nb_indivs()];
  apply_mutation_stop_ = new long[exp_m_->nb_indivs()];
  compute_rna_stop_ = new long[exp_m_->nb_indivs()];
  start_protein_stop_ = new long[exp_m_->nb_indivs()];
  compute_protein_stop_ = new long[exp_m_->nb_indivs()];
  translate_protein_stop_ = new long[exp_m_->nb_indivs()];
  compute_phenotype_stop_ = new long[exp_m_->nb_indivs()];
  compute_fitness_stop_ = new long[exp_m_->nb_indivs()];

  total_start_ = new long[exp_m_->nb_indivs()];
  total_stop_ = new long[exp_m_->nb_indivs()];
  omp_tid_ = new long[exp_m_->nb_indivs()];
  #endif
}

#ifdef HAVE_MPI
int32_t ExpManager_7::xOffset() {
  return x_offset_;
}


int32_t ExpManager_7::yOffset() {
  return y_offset_;
}

int32_t ExpManager_7::rankOf(int32_t x, int32_t y) {
  // printf("Rank Of (%d %d) %d %d : RX %d RY %d R %d\n",x,y,(x / grid_width_),(y / grid_height_),rank_x_,rank_y_,(x / grid_width_) * rank_y_ + (y / grid_height_));
  return (x / grid_width_) * rank_y_ + (y / grid_height_);
}
  
int32_t ExpManager_7::localXtoGlobalX(int32_t lx) {
  return (lx + x_offset_) % global_grid_width_;
}

int32_t ExpManager_7::globalXtoLocalX(int32_t gx) {
  return (gx - local_rank_x_ * grid_width_ + grid_width_)%grid_width_;
}

int32_t ExpManager_7::localYtoGlobalY(int32_t ly) {
  return  (ly + y_offset_) % global_grid_height_;
}

int32_t ExpManager_7::globalYtoLocalY(int32_t gy) {
  return (gy - local_rank_y_ * grid_height_ + grid_height_)%grid_height_;
}

  void ExpManager_7::setRank(int32_t rank) {
    exp_m_->set_rank(rank);
    local_rank_x_ = exp_m_->rank() / rank_y_;
    local_rank_y_ = exp_m_->rank() % rank_y_;

    x_offset_ = (local_rank_x_) * grid_width_; 
    y_offset_ = (local_rank_y_) * grid_height_;
  }
#endif

#ifdef HAVE_MPI
void ExpManager_7::send_broadcast_border() {
  int32_t x_offset = xOffset();
  int32_t y_offset = yOffset();

  MPI_Request request = MPI_REQUEST_NULL;

#ifdef DEBUG_MPI
  printf("%d -- Send fitness\n",exp_m_->rank());
#endif
  // 8 SEND...
  // SEND 1 :
  // 1 value [0,0]

  int32_t i = (-1+x_offset+global_grid_width_) % global_grid_width_;
  int32_t j = (-1+y_offset+global_grid_height_) % global_grid_height_;

  int32_t to_send_rank = rankOf(i,j);

#ifdef DEBUG_MPI
  printf("%d -- (Grid Rank %d %d) Rank %d %d (Offset %d %d): Position %d %d => %d %d (global grid %d %d local grid %d %d): to send %d\n",
          exp_m_->rank(),rank_x,rank_y,local_rank_x,local_rank_y,
          x_offset,y_offset,i,j,
          (i / grid_width), (j / grid_height),
          global_grid_width,global_grid_height,grid_width,grid_height,to_send_rank);
  fflush(stdout);
#endif
    double fitness_to_send = previous_individuals[0]->fitness;

  if (exp_m_->rank() != to_send_rank) {


  #ifdef DEBUG_MPI
    printf("%d -- SF -- STEP 1 -- Send fitness (%d %d -- %d %d) [%d] (%e) to %d :: %d %d %d\n",exp_m_->rank(),i,j,
            x_offset,y_offset,0,fitness_to_send,to_send_rank,(i / grid_width_) , rank_x_ , (j / grid_height_));
    fflush(stdout);
  #endif

    MPI_Isend(&fitness_to_send, 1, MPI_DOUBLE, to_send_rank, 11, MPI_COMM_WORLD, &request);
  }

  // SEND 2 :
  i = (-1+x_offset+global_grid_width_) % global_grid_width_;
  j = ( 0+y_offset+global_grid_height_) % global_grid_height_;
  to_send_rank = rankOf(i,j);
  double fitness_array_to_send[grid_height_];
  if (exp_m_->rank() != to_send_rank) {


      for (int32_t jj = 0; jj < grid_height_; jj++) {
      fitness_array_to_send[jj] = previous_individuals[jj]->fitness;
      // printf("%d -- %d -- RIGHT SEND %d : %e\n",AeTime::time(),exp_m_->rank(),jj,fitness_array_to_send[jj]);
    }
  #ifdef DEBUG_MPI

    printf("%d -- SF -- STEP 2 -- Send fitness array (%d) to %d\n",exp_m_->rank(),grid_height_,to_send_rank);
  #endif
  
    MPI_Isend(&fitness_array_to_send, grid_height_, MPI_DOUBLE, to_send_rank, 12, MPI_COMM_WORLD, &request);
  }

  // SEND 3 :
  i = (-1+x_offset+global_grid_width_) % global_grid_width_;
  j = (grid_height_+y_offset+global_grid_height_) % global_grid_height_;

  to_send_rank = rankOf(i,j);
  if (exp_m_->rank() != to_send_rank) {

    fitness_to_send = previous_individuals[grid_height_-1]->fitness;

  #ifdef DEBUG_MPI
    printf("%d -- SF -- STEP 3 -- Send fitness [%d] (%e) to %d\n",exp_m_->rank(),0,fitness_to_send,to_send_rank);
  #endif

    MPI_Isend(&fitness_to_send, 1, MPI_DOUBLE, to_send_rank, 13, MPI_COMM_WORLD, &request);
  }

  // SEND 4:
  i = ( 0+x_offset+global_grid_width_) % global_grid_width_;
  j = (grid_height_+y_offset+global_grid_height_) % global_grid_height_;
  to_send_rank = rankOf(i,j);
  double fitness_array_to_send_2[grid_width_];
  if (exp_m_->rank() != to_send_rank) {

    for (int32_t jj = 0; jj < grid_width_; jj++) {
      fitness_array_to_send_2[jj] = previous_individuals[(grid_height_-1)+(jj*grid_height_)]->fitness;
    }
  #ifdef DEBUG_MPI
    printf("%d -- SF -- STEP 4 -- Send fitness array (%d) to %d :: %d %d\n",exp_m_->rank(),grid_width_,to_send_rank,i,j);
  #endif

    MPI_Isend(&fitness_array_to_send_2, grid_width_, MPI_DOUBLE, to_send_rank, 14, MPI_COMM_WORLD, &request);
  }

  // SEND 5:
  i = (grid_width_+x_offset+global_grid_width_) % global_grid_width_;
  j = (grid_height_+y_offset+global_grid_height_) % global_grid_height_;

  to_send_rank = rankOf(i,j);
  if (exp_m_->rank() != to_send_rank) {

    fitness_to_send = previous_individuals[grid_width_*grid_height_-1]->fitness;

  #ifdef DEBUG_MPI
    printf("%d -- SF -- STEP 5 -- Send fitness [%d] (%e) to %d\n",exp_m_->rank(),0,fitness_to_send,to_send_rank);
  #endif

    MPI_Isend(&fitness_to_send, 1, MPI_DOUBLE, to_send_rank, 15, MPI_COMM_WORLD, &request);
  }
  // SEND 6:
  i = (grid_width_+x_offset+global_grid_width_) % global_grid_width_;
  j = (grid_height_-1+y_offset+global_grid_height_) % global_grid_height_;
  to_send_rank = rankOf(i,j);
  if (exp_m_->rank() != to_send_rank) {

    for (int32_t jj = 0; jj < grid_height_; jj++) {
      fitness_array_to_send[jj] = previous_individuals[jj+((grid_width_-1)*grid_height_)]->fitness;
    }

  #ifdef DEBUG_MPI
    printf("%d -- SF -- STEP 6 -- Send fitness array (%d) to %d\n",exp_m_->rank(),grid_height_,to_send_rank);
  #endif

    MPI_Isend(&fitness_array_to_send, grid_height_, MPI_DOUBLE, to_send_rank, 16, MPI_COMM_WORLD, &request);
  }
  // SEND 7:
  i = (grid_width_+1+x_offset+global_grid_width_) % global_grid_width_;
  j = (-1+y_offset+global_grid_height_) % global_grid_height_;

  to_send_rank =  rankOf(i,j);
  if (exp_m_->rank() != to_send_rank) {

    fitness_to_send = previous_individuals[(grid_width_-1)*grid_height_]->fitness;

  #ifdef DEBUG_MPI
    printf("%d -- %d -- SF -- STEP 7 -- Send fitness [%d] (%e) to %d (%d %d)\n",AeTime::time(),exp_m_->rank(),0,fitness_to_send,to_send_rank,grid_width_,grid_height_);
  #endif

    MPI_Isend(&fitness_to_send, 1, MPI_DOUBLE, to_send_rank, 17, MPI_COMM_WORLD, &request);
  }

  // SEND 8:
  i = (grid_width_-1+x_offset+global_grid_width_) % global_grid_width_;
  j = (-1+y_offset+global_grid_height_) % global_grid_height_;

  to_send_rank =  rankOf(i,j);

  if (exp_m_->rank() != to_send_rank) {

    double fitness_array_to_send_3[grid_width_];
    for (int32_t jj = 0; jj < grid_width_; jj++) {
      fitness_array_to_send_3[jj] = previous_individuals[(jj*grid_height_)]->fitness;
      // printf("%d -- %d -- SEND DOWN Fitness %d to %d : %e\n",AeTime::time(),exp_m_->rank(),jj,to_send_rank,fitness_array_to_send_2[jj]);
    }

  #ifdef DEBUG_MPI
    printf("%d -- SF -- STEP 8 -- Send fitness array (%d) to %d\n",exp_m_->rank(),grid_width_,to_send_rank);
  #endif

    MPI_Isend(&fitness_array_to_send_3, grid_width_, MPI_DOUBLE, to_send_rank, 18, MPI_COMM_WORLD, &request);
    // printf("%d -- Finishing send fitness...\n",exp_m_->rank());
  }
  // MPI_Barrier(MPI_COMM_WORLD);
    // printf("%d -- Barrier send fitness...\n",exp_m_->rank());

}

void ExpManager_7::recv_broadcast_border() {
  int32_t x_offset = xOffset();
  int32_t y_offset = yOffset();

  MPI_Request request = MPI_REQUEST_NULL;
  MPI_Status status;

#ifdef DEBUG_MPI
  printf("%d -- Receiving fitness...\n",exp_m_->rank());
#endif

  // 8 RECV...
  // SEND 1 :
  int32_t i = (-1+x_offset+global_grid_width_) % global_grid_width_;
  int32_t j = (-1+y_offset+global_grid_height_) % global_grid_height_;
  // printf("Local Grid %d %d Global Grid %d %d Offset %d %d\n",grid_width_,grid_height_,global_grid_width_,global_grid_height_,x_offset_,y_offset_);
  int32_t to_recv_rank = rankOf(i,j);

  double fitness_to_recv;
  // printf("%d -- RF -- STEP 1 BEFORE -- Receive fitness %d (%lf) from %d\n",exp_m_->rank(),0,fitness_to_recv,to_recv_rank);
  // fflush(stdout);

  if (exp_m_->rank() != to_recv_rank) {
    MPI_Irecv(&fitness_to_recv, 1, MPI_DOUBLE, to_recv_rank, 15, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);

  #ifdef DEBUG_MPI
    printf("%d -- RF -- STEP 1 -- Receive fitness %d -- %d %d -- (%e) from %d\n",exp_m_->rank(),0,i,j,fitness_to_recv,to_recv_rank);
  #endif

    fitness_at_border_[TOP][0] = fitness_to_recv;
    fitness_at_border_[LEFT][0] = fitness_to_recv;
  }

  // SEND 2 :
  i = (-1+x_offset+global_grid_width_) % global_grid_width_;
  j = ( 0+y_offset+global_grid_height_) % global_grid_height_;
  to_recv_rank = rankOf(i,j);
  double fitness_array_to_recv[grid_height_];
  if (exp_m_->rank() != to_recv_rank) {
    MPI_Irecv(&fitness_array_to_recv, grid_height_, MPI_DOUBLE, to_recv_rank, 16, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);

    #ifdef DEBUG_MPI
    printf("%d -- RF -- STEP 2 -- Receive fitness array (%d) from %d : %e\n",exp_m_->rank(),grid_width_,to_recv_rank,fitness_array_to_recv[0]);
    #endif

    for (int32_t jj = 0; jj < grid_height_; jj++) {
      fitness_at_border_[LEFT][jj+1] = fitness_array_to_recv[jj];
    }
  }
  // SEND 3 :
  i = (-1+x_offset+global_grid_width_) % global_grid_width_;
  j = (grid_height_+y_offset+global_grid_height_) % global_grid_height_;

  to_recv_rank = rankOf(i,j);
  // printf("%d -- RF -- STEP 3BEFORE -- Receive fitness (%e) from %d\n",exp_m_->rank(),fitness_to_recv,to_recv_rank);

  if (exp_m_->rank() != to_recv_rank) {
    MPI_Irecv(&fitness_to_recv, 1, MPI_DOUBLE, to_recv_rank, 17, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);

  #ifdef DEBUG_MPI
    printf("%d -- %d -- RF -- STEP 3 -- Receive fitness (%e) from %d\n",AeTime::time(),exp_m_->rank(),fitness_to_recv,to_recv_rank);
  #endif

    fitness_at_border_[LEFT][grid_height_+1] = fitness_to_recv;
    fitness_at_border_[DOWN][0] = fitness_to_recv;
  }

  // SEND 4:
  i = ( 0+x_offset+global_grid_width_) % global_grid_width_;
  j = (grid_height_+y_offset+global_grid_height_) % global_grid_height_;
  to_recv_rank = rankOf(i,j);
  double fitness_array_to_recv_2[grid_width_];
  if (exp_m_->rank() != to_recv_rank) {
    // printf("%d -- RF -- STEP 4 BEFORE -- Receive fitness array (%d) from %d :: %d %d Offset %d %d\n",
    //         exp_m_->rank(),grid_width_,to_recv_rank,i,j,x_offset,y_offset);

    MPI_Irecv(&fitness_array_to_recv_2, grid_width_, MPI_DOUBLE, to_recv_rank, 18, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);
  

  #ifdef DEBUG_MPI
    printf("%d -- RF -- STEP 4 -- Receive fitness array (%d) from %d : %e\n",exp_m_->rank(),grid_width_,to_recv_rank,fitness_array_to_recv_2[0]);
  #endif

    for (int32_t jj = 0; jj < grid_width_; jj++) {
      fitness_at_border_[DOWN][jj+1] = fitness_array_to_recv_2[jj];
      // printf("%d -- %d -- Fitness %d from %d : %e\n",AeTime::time(),exp_m_->rank(),jj,to_recv_rank,fitness_array_to_recv_2[jj]);
    }
  }

  // SEND 5:
  i = (grid_width_+x_offset+global_grid_width_) % global_grid_width_;
  j = (grid_height_+y_offset+global_grid_height_) % global_grid_height_;

  to_recv_rank = rankOf(i,j);
  if (exp_m_->rank() != to_recv_rank) {

    MPI_Irecv(&fitness_to_recv, 1, MPI_DOUBLE, to_recv_rank, 11, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);

  #ifdef DEBUG_MPI
    printf("%d -- RF -- STEP 5 -- Receive fitness %d %d (%e) from %d\n",exp_m_->rank(),i,j,fitness_to_recv,to_recv_rank);
  #endif

    fitness_at_border_[RIGHT][grid_height_+1] = fitness_to_recv;
    fitness_at_border_[DOWN][grid_width_+1] = fitness_to_recv;
  }

  // SEND 6:
  i = (grid_width_+x_offset+global_grid_width_) % global_grid_width_;
  j = (grid_height_-1+y_offset+global_grid_height_) % global_grid_height_;
  to_recv_rank = rankOf(i,j);
  if (exp_m_->rank() != to_recv_rank) {

    MPI_Irecv(&fitness_array_to_recv, grid_height_, MPI_DOUBLE, to_recv_rank, 12, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);

  #ifdef DEBUG_MPI
    printf("%d -- RF -- STEP 6 -- Receive fitness array (%d) from %d\n",exp_m_->rank(),grid_width_,to_recv_rank);
  #endif

    for (int32_t jj = 0; jj < grid_height_; jj++) {
      fitness_at_border_[RIGHT][jj+1] = fitness_array_to_recv[jj];
    }
    
  }

  // SEND 7:
  i = (grid_width_+1+x_offset+global_grid_width_) % global_grid_width_;
  j = (-1+y_offset+global_grid_height_) % global_grid_height_;

  to_recv_rank = rankOf(i,j);
  if (exp_m_->rank() != to_recv_rank) {

    MPI_Irecv(&fitness_to_recv, 1, MPI_DOUBLE, to_recv_rank, 13, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);

  #ifdef DEBUG_MPI
    printf("%d -- RF -- STEP 7 -- Receive fitness (%e) from %d\n",exp_m_->rank(),fitness_to_recv,to_recv_rank);
  #endif

    fitness_at_border_[RIGHT][0] = fitness_to_recv;
    fitness_at_border_[TOP][grid_width_+1] = fitness_to_recv;
  }

  // SEND 8:
  i = (grid_width_-1+x_offset+global_grid_width_) % global_grid_width_;
  j = (-1+y_offset+global_grid_height_) % global_grid_height_;

  to_recv_rank = rankOf(i,j);
  if (exp_m_->rank() != to_recv_rank) {

    MPI_Irecv(&fitness_array_to_recv_2, grid_width_, MPI_DOUBLE, to_recv_rank, 14, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);

  #ifdef DEBUG_MPI
    printf("%d -- RF -- STEP 8 -- Receive fitness array (%d) from %d\n",exp_m_->rank(),grid_width_,to_recv_rank,fitness_array_to_recv_2[0]);
  #endif
  // 
    for (int32_t jj = 0; jj < grid_width_; jj++) {
      fitness_at_border_[TOP][jj+1] = fitness_array_to_recv_2[jj];
    }
  }
  // Verify the correct values have been received
}
#endif


void ExpManager_7::selection(int indiv_id) {
  if (exp_m_->sel()->selection_scope() == SCOPE_LOCAL) {
  // printf("Selection %d\n",indiv_id);
  int32_t selection_scope_x = exp_m_->sel()->selection_scope_x();
  int32_t selection_scope_y = exp_m_->sel()->selection_scope_y();
bool is_at_border = false;
#ifdef HAVE_MPI
  int32_t local_x = indiv_id / grid_height_;
  int32_t local_y = indiv_id % grid_height_;
#endif


  int16_t neighborhood_size = selection_scope_x * selection_scope_y;
  FitnessFunction fitness_function = exp_m_->sel()->fitness_func();

  double *  local_fit_array   = new double[(selection_scope_x) * (selection_scope_y)];
  double *  local_meta_array   = new double[(selection_scope_x) * (selection_scope_y)];
  double *  probs             = new double[(selection_scope_x) * (selection_scope_y)];
  int16_t   count             = 0;
  double    sum_local_fit     = 0.0;

#ifdef HAVE_MPI
  int32_t x = localXtoGlobalX(local_x);
  int32_t y = localYtoGlobalY(local_y);
#else
  int32_t x = indiv_id / grid_height_;
  int32_t y = indiv_id % grid_height_;
#endif
  int cur_x,cur_y;

#ifdef __REGUL
  double** fitness_sum_local_tab_;
  int32_t fitness_function_scope_x = exp_m_->sel()->fitness_function_scope_x();
  int32_t fitness_function_scope_y = exp_m_->sel()->fitness_function_scope_y();
#endif

  if (fitness_function == FITNESS_LOCAL_SUM) {
#ifdef __REGUL
    fitness_sum_local_tab_ = new double*[fitness_function_scope_x*fitness_function_scope_y];
        for (int tab_id = 0; tab_id < fitness_function_scope_x*fitness_function_scope_y; tab_id++)
          fitness_sum_local_tab_[tab_id] = new double[phenotypic_target_handler_->nb_env_];

        for (int env_id = 0; env_id < phenotypic_target_handler_->nb_env_; env_id++) {

          int tab_id = 0;
          fitness_sum_local_tab_[tab_id][env_id] = 0;

          for (int8_t i = -(selection_scope_x/2) ; i <= (selection_scope_x/2) ; i++) {
            for (int8_t j = -(selection_scope_y/2) ; j <= (selection_scope_y/2) ; j++) {
              #ifdef HAVE_MPI
              cur_x = (x + i + global_grid_width_) % global_grid_width_;
              cur_y = (y + j + global_grid_height_) % global_grid_height_;
              printf("NOT YET Available with MPI\n");
              exit(-1);
              #else
              cur_x = (x + i + grid_width_) % grid_width_;
              cur_y = (y + j + grid_height_) % grid_height_;
              #endif

              int16_t new_x,new_y;
              for (int8_t ii = -(fitness_function_scope_x/2); ii <= (fitness_function_scope_x/2); ii++) {
                for (int8_t jj = -(fitness_function_scope_y/2); jj <= (fitness_function_scope_y/2); jj++) {
                  //TODO: Check values HERE !

                  new_x = (cur_x + ii + grid_width_) % grid_width_;
                  new_y = (cur_y + jj + grid_height_) % grid_height_;

                  fitness_sum_local_tab_[tab_id][env_id] +=
                      previous_individuals[new_x * grid_height_ + new_y]->fitness_by_env_id_[env_id];
                }
              }
              tab_id++;
            }
          }
        }
#else
    printf("Fitness local sum is not supported for Aevol (only R-Aevol)\n");
    exit(-1);
#endif
  }

  int tab_id = 0;

  for (int8_t i = -(selection_scope_x/2) ; i <= (selection_scope_x/2) ; i++) {
    for (int8_t j = -(selection_scope_y/2) ; j <= (selection_scope_y/2) ; j++) {
      #ifdef HAVE_MPI
      cur_x = (x + i + global_grid_width_) % global_grid_width_;
      cur_y = (y + j + global_grid_height_) % global_grid_height_;

      int32_t owner_rank = rankOf(cur_x,cur_y);
      is_at_border = owner_rank != exp_m_->rank();
      #else
      cur_x = (x + i + grid_width_) % grid_width_;
      cur_y = (y + j + grid_height_) % grid_height_;
      #endif

      if (fitness_function == FITNESS_EXP) {
        #ifdef HAVE_MPI
        if (is_at_border) {
          // int map = -1;
          // printf("%d -- %d -- Stencil Global %d %d (Local %d %d):: Current %d %d :: Local %d %d (Grid %d %d)\n",
          // AeTime::time(),
          // // " :: X/Y (offset) %d (%d) %d (%d)\n",                
          //       x*global_grid_height_+y,
          //       x,y,
          //       local_x,local_y,
          //       cur_x,cur_y,
          //       local_x+i,local_y+j,
          //       grid_width_,grid_height_
          //       // ,
          //       // x,x+i,
          //       // y,y+j
          //       );

          int32_t local_cur_x = local_x + i;
          int32_t local_cur_y = local_y + j;

          // printf("%d -- %d -- Remote FITNESS CUR %d %d LOC CUR %d %d\n",
          // AeTime::time(),x *
          //       #ifdef HAVE_MPI
          //       global_grid_height_
          //       #else
          //       grid_height_
          //       #endif
          //        + y,cur_x,cur_y,local_cur_x,local_cur_y);

          if (local_x+i < 0) {
            local_fit_array[count] = fitness_at_border_[LEFT][local_y+j+1];
            // printf("%d -- %d -- Fitness at border LEFT %d (%d %d -- %d %d => %d %d): %e\n",AeTime::time(),x *
            //     #ifdef HAVE_MPI
            //     global_grid_height_
            //     #else
            //     grid_height_
            //     #endif
            //      + y,local_y+j+1,i,j,x,y,cur_x,cur_y,fitness_at_border_[LEFT][local_y+j+1]);
            // map = 0;
          } else if (local_x+i >= grid_width_) {
            local_fit_array[count] = fitness_at_border_[RIGHT][local_y+j+1];
            // printf("%d -- %d -- Fitness at border RIGHT %d (%d %d -- %d %d => %d %d): %e\n",AeTime::time(),x *
            //     #ifdef HAVE_MPI
            //     global_grid_height_
            //     #else
            //     grid_height_
            //     #endif
            //      + y,local_y+j+1,i,j,x,y,cur_x,cur_y,fitness_at_border_[RIGHT][local_y+j+1]);
            // map = 1;
          } else if (local_y+j < 0) {
            local_fit_array[count] = fitness_at_border_[TOP][local_x+i+1];
            // printf("%d -- %d -- Fitness at border TOP %d (%d %d -- %d %d => %d %d): %e\n",AeTime::time(),x *
            //     #ifdef HAVE_MPI
            //     global_grid_height_
            //     #else
            //     grid_height_
            //     #endif
            //      + y,local_x+i+1,i,j,x,y,cur_x,cur_y,fitness_at_border_[TOP][local_x+i+1]);
            // map = 2;
          } else if (local_y+j >= grid_height_) {
            local_fit_array[count] = fitness_at_border_[DOWN][local_x+i+1];
            // printf("%d -- %d -- Fitness at border DOWN %d (%d %d -- %d %d => %d %d): %e\n",AeTime::time(),x *
            //     #ifdef HAVE_MPI
            //     global_grid_height_
            //     #else
            //     grid_height_
            //     #endif
            //      + y,local_x+i+1,i,j,x,y,cur_x,cur_y,fitness_at_border_[LEFT][local_x+i+1]);
            // map = 3;
          } 
          // else {
          //   printf("SHOULD NOT BERE HERE CENTER %d %d LOOKING %d %d :%d %d: %d %d\n",x,y,cur_x,cur_y,i,j,local_x+i,local_y+i);exit(-1);
          // }
          // printf("%d -- %d -- Fitness for %d %d :: (%d %d) %d : %e\n",                x * 
          //       #ifdef HAVE_MPI
          //       global_grid_height_
          //       #else
          //       grid_height_
          //       #endif
          //        + y,count,x,y,cur_x,cur_y, is_at_border, local_fit_array[count]);

        } else {
          int32_t local_cur_x = (local_x + i + grid_width_) % grid_width_;
          int32_t local_cur_y = (local_y + j + grid_height_) % grid_height_;
          local_fit_array[count] =
            previous_individuals[local_cur_x * grid_height_ + local_cur_y]->fitness;
          // if (x * 
          //       #ifdef HAVE_MPI
          //       global_grid_height_
          //       #else
          //       grid_height_
          //       #endif
          //        + y == 382 && AeTime::time()==6)
          //   printf("%d -- %d -- %d -- Fitness for %d %d :: (%d %d) %d -- %d : %e\n",AeTime::time(),
          //                   x * 
          //       #ifdef HAVE_MPI
          //       global_grid_height_
          //       #else
          //       grid_height_
          //       #endif
          //        + y,count,x,y,cur_x,cur_y, cur_x * global_grid_width_ + cur_y, is_at_border, local_fit_array[count]);
        }

        #else
        local_fit_array[count] =
            previous_individuals[cur_x * grid_height_ + cur_y]->fitness;
        // if (x * 
        //         #ifdef HAVE_MPI
        //         global_grid_height_
        //         #else
        //         grid_height_
        //         #endif
        //          + y == 382 && AeTime::time()==6)
        //   printf("%d -- %d -- Fitness for %d %d :: (%d %d) [%d] : %e\n",                x * 
        //         #ifdef HAVE_MPI
        //         global_grid_height_
        //         #else
        //         grid_height_
        //         #endif
        //          + y,count,x,y,cur_x,cur_y, cur_x * grid_height_ + cur_y, local_fit_array[count]);
        #endif
      } else if (fitness_function == FITNESS_GLOBAL_SUM) {
        #ifdef HAVE_MPI
              printf("NOT YET Available with MPI\n");
              exit(-1);
        #endif
#ifdef __REGUL
        double composed_fitness = 0;
            for (int env_id = 0; env_id < phenotypic_target_handler_->nb_env_;
                 env_id++) {
              composed_fitness +=
                  previous_individuals[cur_x * grid_height_ + cur_y]->fitness_by_env_id_[env_id] /
                  fitness_sum_tab_[env_id];
            }
            composed_fitness /= phenotypic_target_handler_->nb_env_;
            local_fit_array[count] = composed_fitness;
#else
        printf("Fitness local sum is not supported for Aevol (only R-Aevol)\n");
        exit(-1);
#endif
      } else if (fitness_function == FITNESS_LOCAL_SUM) {
        #ifdef HAVE_MPI
              printf("NOT YET Available with MPI\n");
              exit(-1);
        #endif
#ifdef __REGUL
        double composed_fitness = 0;
            for (int env_id = 0; env_id < phenotypic_target_handler_->nb_env_;
                 env_id++) {
              composed_fitness +=
                  previous_individuals[cur_x * grid_height_ + cur_y]->fitness_by_env_id_[env_id] /
                  fitness_sum_local_tab_[tab_id][env_id];
            }
            composed_fitness /= phenotypic_target_handler_->nb_env_;
            local_fit_array[count] = composed_fitness;
#else
        printf("Fitness local sum is not supported for Aevol (only R-Aevol)\n");
        exit(-1);
#endif
      }

      sum_local_fit += local_fit_array[count];
      // printf("%d -- SUM LOCAL FIT %d :: %e (%e)\n",exp_m_->rank(),count,sum_local_fit,local_fit_array[count]);

      count++;
      tab_id++;
    }
  }

  if (fitness_function == FITNESS_LOCAL_SUM) {
#ifdef __REGUL
    for (int tab_id = 0; tab_id < fitness_function_scope_x * fitness_function_scope_y; tab_id++)
          delete[] fitness_sum_local_tab_[tab_id];
        delete[] fitness_sum_local_tab_;
#else
    printf("Fitness local sum is not supported for Aevol (only R-Aevol)\n");
    exit(-1);
#endif
  }

  for(int16_t i = 0 ; i < neighborhood_size ; i++) {
    probs[i] = local_fit_array[i]/sum_local_fit;
    // printf("%d -- %d -- I %d -- P %e (LFit %e SFit %e)\n",
    //       AeTime::time(),
    //             x * 
    //             #ifdef HAVE_MPI
    //             global_grid_height_
    //             #else
    //             grid_height_
    //             #endif
    //              + y,
    //       i,probs[i], local_fit_array[i], sum_local_fit);
  }
  // if (local_x < 0 || local_x >= 3 || local_y < 0 || local_y >= 3)  
    // printf("%d -- Reprod roulette for [ %d %d ] : %p (%d)\n",exp_m_->rank(),local_x,local_y,
    //         exp_m_->world()->grid(local_x,local_y)->reprod_prng_simd_.get(),neighborhood_size);
#ifdef HAVE_MPI
  int16_t found_org = exp_m_->world()->grid(local_x,local_y)->reprod_prng_simd_->roulette_random(probs, neighborhood_size);
#else
  int16_t found_org = exp_m_->world()->grid(x,y)->reprod_prng_simd_->roulette_random(probs, neighborhood_size);
#endif

#ifndef HAVE_MPI
int32_t x_offset,y_offset;

  x_offset = (found_org / selection_scope_x) - 1;
  y_offset = (found_org % selection_scope_x) - 1;
#else
int32_t g_x_offset=xOffset();
int32_t g_y_offset=yOffset();


  g_x_offset = (found_org / selection_scope_x) - 1;
  g_y_offset = (found_org % selection_scope_x) - 1;
#endif


  delete [] local_fit_array;
  delete [] local_meta_array;
  delete [] probs;

#ifdef HAVE_MPI
  
  int32_t owner_rank = rankOf( ((x+g_x_offset+global_grid_width_)  % global_grid_width_),
                               ((y+g_y_offset+global_grid_height_) % global_grid_height_) );

  is_at_border = owner_rank != exp_m_->rank();

  // printf("%d -- %d -- Global ID %d (%d %d) -- SID %d (%d %d) -- RX %d RY %d - GW %d GH %d\n",
  //                                         exp_m_->rank(),indiv_id,
  //                                         ((x+x_offset+global_grid_width)  % global_grid_width)*global_grid_height+
  //                                         ((y+y_offset+global_grid_height) % global_grid_height),
  //                                         ((x+x_offset+global_grid_width)  % global_grid_width),
  //                                         ((y+y_offset+global_grid_height) % global_grid_height),
                                          
  //                                         to_send_rank,
  //                                         ((x+x_offset+global_grid_width) % global_grid_width) / grid_width,
  //                                         ((y+y_offset+global_grid_height) % global_grid_height) / grid_height,
  //                                         grid_width,grid_height,global_grid_height,global_grid_width);
  //  printf("Indiv %d from %d (%d)\n",localXtoGlobalX(indiv_id%grid_height_) * global_grid_height_ + localYtoGlobalY(indiv_id/grid_height_),((x+g_x_offset+global_grid_width_)  % global_grid_width_)*global_grid_height_+
  //                                         ((y+g_y_offset+global_grid_height_) % global_grid_height_),is_at_border);

  if (!is_at_border) {
    // x = indiv_id / grid_height_;
    // y = indiv_id % grid_height_;
    exp_m_->next_generation_reproducer_[indiv_id] = ((x+g_x_offset+global_grid_width_)  % global_grid_width_)*global_grid_height_+
                                          ((y+g_y_offset+global_grid_height_) % global_grid_height_);
    // #ifdef DEBUG_MPI
    // printf("%d -- LOCAL -- R %d -- L_ID %d -- Cell %d (%d %d) -- Reproducer %d (%d %d) -- Local Rank %d\n",
    //       AeTime::time(),exp_m_->rank(),
    //       indiv_id,
    //       localXtoGlobalX(indiv_id%grid_height_) * grid_height_ + localYtoGlobalY(indiv_id/grid_height_),
    //       localXtoGlobalX(indiv_id%grid_height_),
    //       localYtoGlobalY(indiv_id/grid_height_),
    //         ((x+g_x_offset+global_grid_width_)  % global_grid_width_)*global_grid_height_+
    //                                       ((y+g_y_offset+global_grid_height_) % global_grid_height_),
    //        ((x+g_x_offset+global_grid_width_)  % global_grid_width_),
    //                                       ((y+g_y_offset+global_grid_height_) % global_grid_height_),
    //       owner_rank
    //       );
    // #endif
    // printf("%d -- Indiv %d from %d\n",AeTime::time(),
    //                                       ((x+global_grid_width_)  % global_grid_width_)*global_grid_height_+
    //                                       ((y+global_grid_height_) % global_grid_height_),
    //                                       ((x+g_x_offset+global_grid_width_)  % global_grid_width_)*global_grid_height_+
    //                                       ((y+g_y_offset+global_grid_height_) % global_grid_height_));
  } else {
    exp_m_->next_generation_reproducer_[indiv_id] = ((x+g_x_offset+global_grid_width_)  % global_grid_width_)*global_grid_height_+
                                          ((y+g_y_offset+global_grid_height_) % global_grid_height_);
    // printf("%d -- %d (%d %d => G %d %d) is at border fetch to %d (%d %d :=> %d) Found org (%d => %d %d)\n",exp_m_->rank(),indiv_id,
    //         local_x,local_y,x,y,
    //         owner_rank,
    //         ((x+x_offset+global_grid_width)  % global_grid_width),
    //                     ((y+y_offset+global_grid_height) % global_grid_height),
    //                     (((x+x_offset+global_grid_width)  % global_grid_width) ) * global_grid_height
    //                     + ( ((y+y_offset+global_grid_height) % global_grid_height) ) ,
    //                     found_org,x_offset,y_offset);
    // #ifdef DEBUG_MPI
    // printf("%d -- BORDER -- R %d -- L_ID %d -- Global %d (%d %d) -- Reproducer %d %d :: %d -- Remote Rank %d\n",
    //       AeTime::time(),exp_m_->rank(),
    //       indiv_id,
    //       x * grid_height_ +y,
    //       x,
    //       y,
    //       x+g_x_offset,y+g_y_offset,
    //       ((x+g_x_offset)%global_grid_width_)*global_grid_height_+((y+g_y_offset)%global_grid_height_),
    //       owner_rank
    //       );
    // #endif
    // printf("%d -- Indiv %d from %d\n",AeTime::time(),((x+global_grid_width_)  % global_grid_width_)*global_grid_height_+
    //                                       ((y+global_grid_height_) % global_grid_height_),
    //                                       ((x+g_x_offset+global_grid_width_)  % global_grid_width_)*global_grid_height_+
    //                                       ((y+g_y_offset+global_grid_height_) % global_grid_height_));
    #pragma omp critical(update_ifb)
    {
    individual_to_fetch_at_border_.push_back(
      ToFetchIndividual((((x+g_x_offset+global_grid_width_)  % global_grid_width_) ) * global_grid_height_
                        + ( ((y+g_y_offset+global_grid_height_) % global_grid_height_) ) ,
                        ((x+g_x_offset+global_grid_width_)  % global_grid_width_),
                        ((y+g_y_offset+global_grid_height_) % global_grid_height_),
                        owner_rank));
    }
  }

#else
  exp_m_->next_generation_reproducer_[indiv_id] = ((x+x_offset+grid_width_)  % grid_width_)*grid_height_+
                                          ((y+y_offset+grid_height_) % grid_height_);
//  printf("%d -- Indiv %d from %d\n",AeTime::time(),indiv_id,exp_m_->next_generation_reproducer_[indiv_id]);
#endif
      // fflush(stdout);
      // #ifdef HAVE_MPI
      // FILE *fp;
      // char filename[20];
      // sprintf(filename,"test_log_%d",exp_m_->rank());
      // fp = fopen(filename, "a");
      // fprintf(fp,"%d -- Indiv %d from %d (%d)\n",AeTime::time(),
      //           x * 
      //           #ifdef HAVE_MPI
      //           global_grid_height_
      //           #else
      //           grid_height_
      //           #endif
      //            + y,
      //           exp_m_->next_generation_reproducer_[indiv_id],found_org);
      // // fprintf(fp,"%d -- %d is remote ? %d\n",AeTime::time(),x * 
      // //           #ifdef HAVE_MPI
      // //           global_grid_height_
      // //           #else
      // //           grid_height_
      // //           #endif
      // //            + y,is_at_border);                
      // fclose(fp);
      // #else
      // printf("%d -- Indiv %d from %d (%d)\n",AeTime::time(),
      //           x * 
      //           #ifdef HAVE_MPI
      //           global_grid_height_
      //           #else
      //           grid_height_
      //           #endif
      //            + y,
      //           exp_m_->next_generation_reproducer_[indiv_id],found_org);
      // // printf("%d -- %d is remote ? %d\n",AeTime::time(),x * 
      // //           #ifdef HAVE_MPI
      // //           global_grid_height_
      // //           #else
      // //           grid_height_
      // //           #endif
      // //            + y,is_at_border);
      // #endif
                // fflush(stdout);
//     printf("%d -- Indiv %d from %d (%d %d :: %d) (%d %d) (%d %d)\n",AeTime::time(), 
//                 // indiv_id,
//                 x * 
//                 #ifdef HAVE_MPI
//                 global_grid_height_
//                 #else
//                 grid_height_
//                 #endif
//                  + y,
//                 exp_m_->next_generation_reproducer_[indiv_id],
//                 #ifdef HAVE_MPI
//                 x+g_x_offset,
//                 y+g_y_offset,
//                 #else
//                 x+x_offset,
//                 y+y_offset,
//                 #endif
//                 #ifdef HAVE_MPI
// (((x+g_x_offset+global_grid_width_)  % global_grid_width_) ) * global_grid_height_
//                         + ( ((y+g_y_offset+global_grid_height_) % global_grid_height_) ),
//                 #else
// ((x+x_offset+grid_width_)  % grid_width_)*grid_height_+
//                                           ((y+y_offset+grid_height_) % grid_height_),
//                 #endif
//                                 #ifdef HAVE_MPI
//                 x+g_x_offset,
//                 y+g_y_offset,
//                 #else
//                 x+x_offset,
//                 y+y_offset,
//                 #endif
//                 #ifdef HAVE_MPI
//                 g_x_offset,
//                 g_y_offset
//                 #else
//                 x_offset,
//                 y_offset
//                 #endif
//                 );
                // ,exp_m_->rank(),
                // local_x,local_y,x,y);

  }

}

#ifdef HAVE_MPI
void ExpManager_7::Fetch_Remote_Individual() {
  MPI_Request request = MPI_REQUEST_NULL;
  MPI_Status status;
  
  int32_t selection_scope_x = exp_m_->sel()->selection_scope_x();
  int32_t selection_scope_y = exp_m_->sel()->selection_scope_y();
  int16_t neighborhood_size = selection_scope_x * selection_scope_y - 1;


//       int32_t owner_rank = (cur_x / rank_x) + (cur_y / rank_y) * rank_y;
  int32_t neighbourhood_rank[8];
  int32_t neighbourhood_nb_request[8];
  for (int i = 0; i < neighborhood_size; i++) {
    int32_t nb_indiv = 0;
    int32_t rank_i;
    int32_t rank_i_x;
    int32_t rank_i_y;
    neighbourhood_nb_request[i] = 0;

    switch(i) {
      case 0:
        rank_i_x = (-1 + local_rank_x_ + rank_x_) % rank_x_;
        rank_i_y = (-1 + local_rank_y_ + rank_y_) % rank_y_;
        rank_i = rank_i_x * rank_y_ + rank_i_y ;
        break;
      case 1:
        rank_i_x = (-1 + local_rank_x_+ rank_x_) % rank_x_;
        rank_i_y = ( 0 + local_rank_y_+ rank_y_) % rank_y_;
        rank_i = rank_i_x * rank_y_ + rank_i_y ;
        break;
      case 2:
        rank_i_x = (-1 + local_rank_x_+ rank_x_) % rank_x_;
        rank_i_y = (+1 + local_rank_y_+ rank_y_) % rank_y_;
        rank_i = rank_i_x * rank_y_ + rank_i_y ;
        break;
      case 3:
        rank_i_x = ( 0 + local_rank_x_+ rank_x_) % rank_x_;
        rank_i_y = (+1 + local_rank_y_+ rank_y_) % rank_y_;
        rank_i = rank_i_x * rank_y_ + rank_i_y ;
        break;
      case 4:
        rank_i_x = (+1 + local_rank_x_+ rank_x_) % rank_x_;
        rank_i_y = (+1 + local_rank_y_+ rank_y_) % rank_y_;
        rank_i = rank_i_x * rank_y_ + rank_i_y ;
        break;
      case 5:
        rank_i_x = (+1 + local_rank_x_+ rank_x_) % rank_x_;
        rank_i_y = ( 0 + local_rank_y_+ rank_y_) % rank_y_;
        rank_i = rank_i_x * rank_y_ + rank_i_y ;
        break;
      case 6:
        rank_i_x = (+1 + local_rank_x_+ rank_x_) % rank_x_;
        rank_i_y = (-1 + local_rank_y_+ rank_y_) % rank_y_;
        rank_i = rank_i_x * rank_y_ + rank_i_y ;
        break;
      case 7:
        rank_i_x = (0 + local_rank_x_+ rank_x_) % rank_x_;
        rank_i_y = (-1 + local_rank_y_+ rank_y_) % rank_y_;
        rank_i = rank_i_x * rank_y_ + rank_i_y ;
        break;
    }

    // printf("%d [%d %d]-- N %d -- RANK [%d] = [%d %d]\n",exp_m_->rank(),local_rank_x,local_rank_y,i,rank_i,rank_i_x,rank_i_y);
    neighbourhood_rank[i] = rank_i;

      for (auto indiv_to_fetch = individual_to_fetch_at_border_.begin(); 
            indiv_to_fetch != individual_to_fetch_at_border_.end(); indiv_to_fetch++ ) {
              if ((*indiv_to_fetch).rank_ == rank_i) {
                bool toSend = true;

                if (indiv_to_fetch !=  individual_to_fetch_at_border_.begin()) {
                  for (auto ifetch = individual_to_fetch_at_border_.begin(); ifetch != indiv_to_fetch; ifetch++) {
                      // printf("%d -- Looking for %d AT %d\n",exp_m_->rank(),(*indiv_to_fetch).id_, (*ifetch).id_);
                      if ((*indiv_to_fetch).id_ == (*ifetch).id_)
                        toSend = false;
                  } 
                }
                  
                if (toSend) neighbourhood_nb_request[i]++;
            }
      
      }
    // for (auto indiv_to_fetch : individual_to_fetch_at_border_) {
    //   if (indiv_to_fetch.rank_ == rank_i) neighbourhood_nb_request[i]++;
    // }
  }
      // for (auto indiv_to_fetch = individual_to_fetch_at_border_.begin();
      //             indiv_to_fetch != individual_to_fetch_at_border_.end(); indiv_to_fetch++ ) {
      //               printf("%d -- TOFETCH %d from %d\n",exp_m_->rank(),(*indiv_to_fetch).id_,(*indiv_to_fetch).rank_);
      //             } 


  // for (int i = 0; i < neighborhood_size; i++)
    // printf("%d -- NB Request to %d : %d\n",exp_m_->rank(),neighbourhood_rank[i],neighbourhood_nb_request[i]);

  for (int i = 0; i < neighborhood_size; i++) {
    bool toSend = true;
    
    if (i-1 >= 0)
      for (int j = 0; j <= i-1; j++) {
        if (neighbourhood_rank[j] == neighbourhood_rank[i])
          toSend = false;
      }

    if (toSend && neighbourhood_rank[i] != exp_m_->rank()) {
      #ifdef DEBUG_MPI
      printf("%d -- Request %d indivs from %d\n",exp_m_->rank(),
             neighbourhood_nb_request[i],neighbourhood_rank[i]);
      #endif

      MPI_Isend(&neighbourhood_nb_request[i], 1, MPI_INT, neighbourhood_rank[i], 145, MPI_COMM_WORLD, &request);
    }
  }

  // MPI_Barrier(MPI_COMM_WORLD);

  int32_t to_send_dna_per_rank[8];
  for (int i = 0; i < neighborhood_size; i++) {
        bool toRecv = true;
    
    if (i-1 >= 0)
      for (int j = 0; j <= i-1; j++) {
        if (neighbourhood_rank[j] == neighbourhood_rank[i])
          toRecv = false;
      }

    if (toRecv && neighbourhood_rank[i] != exp_m_->rank()) {
      MPI_Recv(&to_send_dna_per_rank[i], 1, MPI_INT, neighbourhood_rank[i], 145, MPI_COMM_WORLD, &status);//, &request);
      // MPI_Wait(&request, &status);
      #ifdef DEBUG_MPI
      printf("%d -- %d is requesting %d indivs\n",exp_m_->rank(),neighbourhood_rank[i],to_send_dna_per_rank[i]);
      #endif
    } else {
      to_send_dna_per_rank[i] = 0;
    }
  }

  // Request DNA size
  int32_t index_per_rank[8];

  for (int i = 0; i < neighborhood_size; i++) {
      index_per_rank[i] = 0;
  }

  for (auto indiv_to_fetch = individual_to_fetch_at_border_.begin(); 
            indiv_to_fetch != individual_to_fetch_at_border_.end(); indiv_to_fetch++ ) {
          bool toSend = true;

          if (indiv_to_fetch !=  individual_to_fetch_at_border_.begin()) {

            for (auto ifetch = individual_to_fetch_at_border_.begin(); ifetch != indiv_to_fetch; ifetch++) {
  //               // printf("%d -- Looking for %d AT %d\n",exp_m_->rank(),(*indiv_to_fetch).id_, (*ifetch).id_);
                if ((*indiv_to_fetch).id_ == (*ifetch).id_)
                  toSend = false;
            } 
          }

    if (toSend) {
  //           // Request DNA
            #ifdef DEBUG_MPI
            printf("%d -- Request DNA size for %d from %d\n",exp_m_->rank(),
                    (*indiv_to_fetch).id_,(*indiv_to_fetch).rank_);
            #endif
            MPI_Isend(&(*indiv_to_fetch).id_, 1, MPI_INT, (*indiv_to_fetch).rank_, 567, MPI_COMM_WORLD, &request);
  //           (*indiv_to_fetch).mpi_id_ = index_per_rank[(*indiv_to_fetch).rank_];
  //           index_per_rank[(*indiv_to_fetch).rank_]++;
  //   } else {
  //     printf("%d -- Indiv %d is already requested to %d\n",exp_m_->rank(),
  //                   (*indiv_to_fetch).id_,(*indiv_to_fetch).rank_);
    }
        
  }
  
  // MPI_Barrier(MPI_COMM_WORLD);

  // printf("%d -- Receiving DNA requests\n ",exp_m_->rank());
  for (int i = 0; i < neighborhood_size; i++) {
    bool toSend = true;
    if (i-1 >=0 )
      for (int j = 0; j < i; j++) {
        if (neighbourhood_rank[i] == neighbourhood_rank[j]) toSend = false;
      }

    if (toSend) {
      // printf("%d -- %d --  Number of indiv requested %d by %d\n",exp_m_->rank(),i,to_send_dna_per_rank[i],neighbourhood_rank[i]);
      for (int j = 0; j < to_send_dna_per_rank[i]; j++) {
        // #ifdef DEBUG_MPI
              // printf("%d -- %d -- Which indiv %d is requesting ? (NB %d) \n",exp_m_->rank(),i, 
              //               neighbourhood_rank[i], to_send_dna_per_rank[i]);
        // #endif

        // Receive DNA size request
        int32_t indiv_id = -1;
        MPI_Recv(&indiv_id, 1, MPI_INT, neighbourhood_rank[i], 567, MPI_COMM_WORLD, &status);
        // MPI_Wait(&request, &status);
        #ifdef DEBUG_MPI
        printf("%d -- %d is requesting DNA SIZE for Individual %d\n",exp_m_->rank(),neighbourhood_rank[i],indiv_id);
        #endif
        // Send DNA size
        // for (auto indiv_to_fetch : individual_to_fetch_at_border_) {
        //   if (indiv_to_fetch.rank_ == neighbourhood_rank[i] && indiv_to_fetch.id_ == indiv_id) {
            int32_t x = (indiv_id / global_grid_height_);
            int32_t y = (indiv_id % global_grid_height_);

            int32_t local_x = globalXtoLocalX(x);
            int32_t local_y = globalYtoLocalY(y);
            int32_t local_id_ = local_x * grid_height_ + local_y;

            // printf("%d -- Get local IDX %d (%d %d -- %d) (offset %d %d) for global IDX %d (%d %d) rank %d %d (%d %d)\n",exp_m_->rank(),
            //         local_id_,local_x,local_y,global_grid_height,
            //         local_rank_x * rank_x,local_rank_y * rank_y,
            //         indiv_id,
            //         x,y,
            //         local_rank_x,local_rank_y,rank_x,rank_y);
            
            // printf("%d -- ==========> Reading info for %d (LOC %d %d) Global %d (GLOB %d %d) (BEFORE SENDING)\n",
            //         exp_m_->rank(),local_id_,local_x,local_y,indiv_id,x,y);
            int32_t dna_length = previous_individuals[local_id_]->dna_->length_;
            char* dna = previous_individuals[local_id_]->dna_->data_;
  // #ifdef DEBUG_MPI

            // printf("%d -- %d -- Sending DNA SIZE (%d) and DNA for Individual %d (L %d -- %d %d) to %d TAG %d\n",
            //         AeTime::time(),
            //         exp_m_->rank(),dna_length,indiv_id,
            //         local_id_,local_x,local_y,
            //         neighbourhood_rank[i],indiv_id*10+1);
            // fflush(stdout);
// #endif
            // auto pair = std::make_pair(indiv_id,dna_length);

            // MPI_Isend(&pair, sizeof(std::pair<int,int>), MPI_BYTE, neighbourhood_rank[i], indiv_id+1024, MPI_COMM_WORLD, &request);
            MPI_Isend(dna, dna_length+1, MPI_CHAR, neighbourhood_rank[i], indiv_id*10+1, MPI_COMM_WORLD, &request);

            int32_t lead_prom_size = 0;
            int32_t lag_prom_size = 0;
            

            int lead_idx = 0;
            int lag_idx = 0;
              
              for (int rna_idx = 0; rna_idx <
                                  (int)previous_individuals[local_id_]->metadata_->promoter_count(); rna_idx++) {
              
                if (previous_individuals[local_id_]->metadata_->promoters(rna_idx) != nullptr) {
                  if (previous_individuals[local_id_]->metadata_->promoters(rna_idx)->leading_or_lagging) {
                                lead_prom_size++;
                  } else {
                                lag_prom_size++;
                  }
                }
              }

                          // printf("TOSEND -- Leading/Lagging Length %d %d\n",lead_prom_size,lag_prom_size);


            int32_t* lead_prom_pos = new int32_t[lead_prom_size];
            int8_t* lead_prom_error = new int8_t[lead_prom_size];
            
            int32_t* lag_prom_pos = new int32_t[lag_prom_size];
            int8_t* lag_prom_error = new int8_t[lag_prom_size];

            for (int rna_idx = 0; rna_idx <
                                  (int)previous_individuals[local_id_]->metadata_->promoter_count(); rna_idx++) {
              
                if (previous_individuals[local_id_]->metadata_->promoters(rna_idx) != nullptr) {
                  if (previous_individuals[local_id_]->metadata_->promoters(rna_idx)->leading_or_lagging) {
                                lead_prom_pos[lead_idx] = previous_individuals[local_id_]->metadata_->promoters(rna_idx)->pos;
                                lead_prom_error[lead_idx] = previous_individuals[local_id_]->metadata_->promoters(rna_idx)->error;
                                // printf("%d :: LEAD :: %d %d\n",indiv_id,lead_prom_pos[lead_idx],lead_prom_error[lead_idx]);
                                lead_idx++;
                  } else {
                                lag_prom_pos[lag_idx] = previous_individuals[local_id_]->metadata_->promoters(rna_idx)->pos;
                                lag_prom_error[lag_idx] = previous_individuals[local_id_]->metadata_->promoters(rna_idx)->error;
                                // printf("%d :: LAG :: %d %d\n",indiv_id,lag_prom_pos[lag_idx],lag_prom_error[lag_idx]);
                                lag_idx++;
                  }
                }
              }
              #pragma omp critical(update_isb)
              {
                individual_sent_at_border_.push_back(
                  ToFetchIndividual(indiv_id,x,y,-1));
                individual_sent_at_border_.back().lead_prom_pos = lead_prom_pos;
                individual_sent_at_border_.back().lead_prom_error = lead_prom_error;

                individual_sent_at_border_.back().lag_prom_pos = lag_prom_pos;
                individual_sent_at_border_.back().lag_prom_error = lag_prom_error;
              }


      // if ((current_individuals[indiv_id]->is_at_border_)) {
        // printf("T %d -- %d -- Individual %d (Global %d) : SEND to %d (SIZE %d) from scratch brder %d (%d %d):: %e\n",AeTime::time(),
        // exp_m_->rank(),local_id_,indiv_id,
        // neighbourhood_rank[i],dna_length
        // ,previous_individuals[local_id_]->metadata_->promoter_count(),
        //                       lead_prom_size,lag_prom_size,previous_individuals[local_id_]->fitness);
      // }


            // // for (auto& prom: ((List_Metadata*)previous_individuals[local_id_]->metadata_)->promoters_list_[LEADING]) {
            // for (auto it = ((List_Metadata*)previous_individuals[local_id_]->metadata_)->promoters_list_[LEADING].begin(); 
            //         it != ((List_Metadata*)previous_individuals[local_id_]->metadata_)->promoters_list_[LEADING].end(); it++) {
            //           if ((it) != nullptr) {
            //             lead_prom_pos[idx] = it->pos;
            //             lead_prom_error[idx] = it->error;
            //             printf("%d :: LEAD :: %d %d\n",indiv_id,lead_prom_pos[idx],lead_prom_error[idx]);
            //             idx++;
            //           }
            // }


            // idx = 0;
            // // for (auto& prom: ((List_Metadata*)previous_individuals[local_id_]->metadata_)->promoters_list_[LAGGING]) {
            // for (auto it = ((List_Metadata*)previous_individuals[local_id_]->metadata_)->promoters_list_[LAGGING].begin(); 
            //         it != ((List_Metadata*)previous_individuals[local_id_]->metadata_)->promoters_list_[LAGGING].end(); it++) {
            //           if ((it) != nullptr) {
            //             lag_prom_pos[idx] = it->pos;
            //             lag_prom_error[idx] = it->error;
            //             printf("%d :: LAG :: %d %d\n",indiv_id,lead_prom_pos[idx],lead_prom_error[idx]);
            //             idx++;
            //           }
            // }

            // printf("%d -- Send to %d : 4 messages size (%d %d)\n",exp_m_->rank(),neighbourhood_rank[i],lead_prom_size,lag_prom_size);

            MPI_Isend(lead_prom_pos, lead_prom_size, MPI_INT32_T, neighbourhood_rank[i], indiv_id*10+2, MPI_COMM_WORLD, &request);
            MPI_Isend(lead_prom_error, lead_prom_size, MPI_INT8_T, neighbourhood_rank[i], indiv_id*10+3, MPI_COMM_WORLD, &request);
            MPI_Isend(lag_prom_pos, lag_prom_size, MPI_INT32_T, neighbourhood_rank[i], indiv_id*10+4, MPI_COMM_WORLD, &request);
            MPI_Isend(lag_prom_error, lag_prom_size, MPI_INT8_T, neighbourhood_rank[i], indiv_id*10+5, MPI_COMM_WORLD, &request);
// #ifdef DEBUG_MPI
            // printf("%d -- DNA %d sent to %d\n",
            //         exp_m_->rank(),dna_length,indiv_id,
            //         local_id_,local_x,local_y,
            //         neighbourhood_rank[i]);
// #endif
        //   }
        // }
      } 
    }
  }

  #ifdef DEBUG_MPI
  printf("%d -- Finish Sending DNA\n",exp_m_->rank());
  #endif

    // MPI_Barrier(MPI_COMM_WORLD);
// return;
  
#ifdef DEBUG_MPI
  printf("%d -- Receiving DNA\n",exp_m_->rank());
#endif

  // for (int i = 0; i < neighborhood_size; i++) {
      // for (auto indiv_to_fetch = individual_to_fetch_at_border_.begin(); 
      //       indiv_to_fetch != individual_to_fetch_at_border_.end(); indiv_to_fetch++ ) {
      //     bool toSend = true;

      //     if (indiv_to_fetch !=  individual_to_fetch_at_border_.begin()) {

      //       for (auto ifetch = individual_to_fetch_at_border_.begin(); ifetch != indiv_to_fetch; ifetch++) {
      //           if ((*indiv_to_fetch).id_ == (*ifetch).id_)
      //             toSend = false;
      //       } 
      //     }
  for (int i = 0; i < neighborhood_size; i++) {
    bool toSend = true;
    if (i-1 >=0 )
    for (int j = 0; j < i; j++) {
      if (neighbourhood_rank[i] == neighbourhood_rank[j]) toSend = false;
    }

    if (toSend && neighbourhood_rank[i] != exp_m_->rank())
      for (int j = 0; j < neighbourhood_nb_request[i]*5; j ++) {
          #ifdef DEBUG_MPI
        printf("%d -- Fetching %d out of %d from %d\n",exp_m_->rank(),j,neighbourhood_nb_request[i],neighbourhood_rank[j]);
        #endif
            // if (toSend) {
              // Fetch DNA size
              // int32_t dna_length = -1;
              // #ifdef DEBUG_MPI
              // printf("%d -- Fetching DNA size for %d from %d : %d\n",exp_m_->rank(),(*indiv_to_fetch).id_,(*indiv_to_fetch).rank_,
              //        (*indiv_to_fetch).mpi_id_);
              // #endif

              // std::pair<int,int> pair_i;
              // MPI_Recv(&pair_i, sizeof(std::pair<int,int>), MPI_BYTE, neighbourhood_rank[i],  MPI::ANY_TAG, MPI_COMM_WORLD,&status);
              // MPI_Wait(&request, &status);
              // #ifdef DEBUG_MPI
              // #endif

              // printf("%d -- Request message from %d : %d out of %d\n",exp_m_->rank(),neighbourhood_rank[i],j,neighbourhood_nb_request[i]);


              MPI_Probe(neighbourhood_rank[i], MPI_ANY_TAG, MPI_COMM_WORLD, &status);



              int32_t fetched_indiv_id = status.MPI_TAG / 10;
              int32_t fetched_type_id = status.MPI_TAG % 10;

                  // MPI_Wait(&request, &status);
              // printf("%d -- Received message (tag %d) type %d for individual %d\n",exp_m_->rank(),status.MPI_TAG,fetched_type_id,fetched_indiv_id);

              for (auto indiv_to_fetch = individual_to_fetch_at_border_.begin(); 
                indiv_to_fetch != individual_to_fetch_at_border_.end(); indiv_to_fetch++ ) {
                if ((neighbourhood_rank[i] == ((*indiv_to_fetch).rank_)) && (fetched_indiv_id == (*indiv_to_fetch).id_)) {
                  (*indiv_to_fetch).local_x_ = globalXtoLocalX((*indiv_to_fetch).id_ / global_grid_height_);
                  (*indiv_to_fetch).local_y_ = globalYtoLocalY((*indiv_to_fetch).id_ % global_grid_height_);
                  // int32_t local_id_ = local_x * grid_height + local_y;
                  if (fetched_type_id == 1) {
                    int32_t dna_length = -1;
                    MPI_Get_count(&status, MPI_CHAR, &dna_length);

                    char* dna = (char*)malloc((dna_length+1)*sizeof(char));
                    // printf("%d -- START -- Receive DNA %d for %d\n",exp_m_->rank(),dna_length,neighbourhood_rank[i]);
                    MPI_Recv(dna, dna_length, MPI_CHAR, neighbourhood_rank[i], fetched_indiv_id*10+fetched_type_id, MPI_COMM_WORLD,&status);
                    // printf("%d -- END -- Receive DNA %d for %d\n",exp_m_->rank(),dna_length,neighbourhood_rank[i]);
                    (*indiv_to_fetch).dna_length_ = dna_length-1;
                    (*indiv_to_fetch).dna_ = dna;
                    (*indiv_to_fetch).fetched_ = true;
                  }
                   else if (fetched_type_id == 2) {
                    int32_t lead_prom_size = -1;
                    MPI_Get_count(&status, MPI_INT32_T, &lead_prom_size);

                    int32_t* lead_prom_pos = (int32_t*)malloc((lead_prom_size)*sizeof(int32_t));
                    MPI_Recv(lead_prom_pos, lead_prom_size, MPI_INT32_T, neighbourhood_rank[i], fetched_indiv_id*10+fetched_type_id, MPI_COMM_WORLD,&status);
                    
                    (*indiv_to_fetch).lead_prom_size = lead_prom_size;
                    (*indiv_to_fetch).lead_prom_pos = lead_prom_pos;
                    (*indiv_to_fetch).fetched_ = true;

                    // if (fetched_indiv_id == 207) {
                      // printf("%d -- (%d %d :: %d -- %d) -- Update lead_prom_size %d == %d\n",exp_m_->rank(),fetched_indiv_id,fetched_type_id,(*indiv_to_fetch).id_,
                      //                 status.MPI_TAG,lead_prom_size,(*indiv_to_fetch).lead_prom_size);
                    // }
                  }
                   else if (fetched_type_id == 3) {
                    int32_t lead_prom_size = -1;
                    MPI_Get_count(&status, MPI_INT8_T, &lead_prom_size);

                    int8_t* lead_prom_error = (int8_t*)malloc((lead_prom_size)*sizeof(int8_t));
                    MPI_Recv(lead_prom_error, lead_prom_size, MPI_INT8_T, neighbourhood_rank[i], fetched_indiv_id*10+fetched_type_id, MPI_COMM_WORLD,&status);
                    
                    (*indiv_to_fetch).lead_prom_size = lead_prom_size;
                    (*indiv_to_fetch).lead_prom_error = lead_prom_error;
                    (*indiv_to_fetch).fetched_ = true;

                    // if (fetched_indiv_id == 207) {
                      // printf("%d -- (%d %d :: %d -- %d) -- Update lead_prom_size ERROR %d\n",exp_m_->rank(),fetched_indiv_id,
                      //       fetched_type_id,(*indiv_to_fetch).id_,status.MPI_TAG,lead_prom_size);
                    // }
                  } else if (fetched_type_id == 4) {
                    int32_t lag_prom_size = -1;
                    MPI_Get_count(&status, MPI_INT32_T, &lag_prom_size);

                    int32_t* lag_prom_pos = (int32_t*)malloc((lag_prom_size)*sizeof(int32_t));
                    MPI_Recv(lag_prom_pos, lag_prom_size, MPI_INT32_T, neighbourhood_rank[i], fetched_indiv_id*10+fetched_type_id, MPI_COMM_WORLD,&status);
                    
                    (*indiv_to_fetch).lag_prom_size = lag_prom_size;
                    (*indiv_to_fetch).lag_prom_pos = lag_prom_pos;
                    (*indiv_to_fetch).fetched_ = true;

                    // if (fetched_indiv_id == 207) {
                      // printf("%d -- (%d %d :: %d -- %d) -- Update lag_prom_size %d\n",exp_m_->rank(),fetched_indiv_id,fetched_type_id,
                      //           (*indiv_to_fetch).id_,status.MPI_TAG,lag_prom_size);
                    // }
                  } else if (fetched_type_id == 5) {
                    int32_t lag_prom_size = -1;
                    MPI_Get_count(&status, MPI_INT8_T, &lag_prom_size);

                    int8_t* lag_prom_error = (int8_t*)malloc((lag_prom_size)*sizeof(int8_t));
                    MPI_Recv(lag_prom_error, lag_prom_size, MPI_INT8_T, neighbourhood_rank[i], fetched_indiv_id*10+fetched_type_id, MPI_COMM_WORLD,&status);
                    
                    (*indiv_to_fetch).lag_prom_size = lag_prom_size;
                    (*indiv_to_fetch).lag_prom_error = lag_prom_error;
                    (*indiv_to_fetch).fetched_ = true;

                    // if (fetched_indiv_id == 207) {
                      // printf("%d -- (%d %d :: %d -- %d) -- Update lag_prom_size ERROR %d\n",exp_m_->rank(),fetched_indiv_id,fetched_type_id,
                      //           (*indiv_to_fetch).id_,status.MPI_TAG,lag_prom_size);
                    // }
                  }
      // // #ifdef DEBUG_MPI
                  // if (fetched_type_id == 1) printf("%d -- Fetched DNA for %d from %d : %d\n",exp_m_->rank(),(*indiv_to_fetch).id_,(*indiv_to_fetch).rank_,(*indiv_to_fetch).dna_length_);
      // // #endif
                  break;
                }

              }
                              // printf("%d -- Go to next message\n",exp_m_->rank());
        }
      }
  

  // MPI_Barrier(MPI_COMM_WORLD);
  // MPI_Finalize();
// #ifdef DEBUG_MPI
  // printf("%d -- Fetch DNA is OVER\n",exp_m_->rank());
// #endif
    // exit(1);
  // printf("Initializing MPI Environment\n");
  // MPI_Init(NULL, NULL);
}
#endif

void ExpManager_7::check_selection(int indiv_id) {
  int32_t selection_scope_x = exp_m_->sel()->selection_scope_x();
  int32_t selection_scope_y = exp_m_->sel()->selection_scope_y();

  int32_t grid_width = exp_m_->grid_width();
  int32_t grid_height = exp_m_->grid_height();

  int16_t neighborhood_size = selection_scope_x * selection_scope_y;

  double *  local_fit_array   = new double[neighborhood_size];
  double *  local_meta_array   = new double[neighborhood_size];
  double *  probs             = new double[neighborhood_size];
  int16_t   count             = 0;
  double    sum_local_fit     = 0.0;

  int * indiv_index = new int[neighborhood_size];

  int32_t x = indiv_id / grid_height;
  int32_t y = indiv_id % grid_height;

  int cur_x,cur_y;

  for (int8_t i = -(selection_scope_x/2) ; i <= (selection_scope_x/2) ; i++) {
    for (int8_t j = -(selection_scope_y/2) ; j <= (selection_scope_y/2) ; j++) {
      cur_x = (x + i + grid_width) % grid_width;
      cur_y = (y + j + grid_height) % grid_height;


      local_fit_array[count] =
          previous_individuals[cur_x * grid_height + cur_y]->fitness;

      local_meta_array[count] =
          previous_individuals[cur_x * grid_height + cur_y]->metaerror;

      if (local_meta_array[count] !=
          previous_individuals[cur_x * grid_height + cur_y]->metaerror) {
        printf("NONONNONONONN\n");
        exit(-1);
      }
      sum_local_fit += local_fit_array[count];


      indiv_index[count] = cur_x * grid_height + cur_y;
      count++;
    }
  }

  for(int16_t i = 0 ; i < neighborhood_size ; i++) {
    probs[i] = local_fit_array[i]/sum_local_fit;
  }

  int16_t found_org = exp_m_->world()->grid(x,y)->reprod_prng_simd_->roulette_random(probs, neighborhood_size);

  int16_t x_offset = (found_org / selection_scope_x) - 1;
  int16_t y_offset = (found_org % selection_scope_x) - 1;


  int found_id = ((x+x_offset+grid_width)  % grid_width)*grid_height+
                 ((y+y_offset+grid_height) % grid_height);


  if (found_id != exp_m_->next_generation_reproducer_[indiv_id]) {


    printf("For individual %d: Selection is diff SIMD %d CPU %d (Meta error %f -- %f || Fitness %e -- %e) \n",
           indiv_id, found_id, exp_m_->next_generation_reproducer_[indiv_id],
        previous_individuals[exp_m_->next_generation_reproducer_[indiv_id]]->metaerror,
        previous_individuals[found_id]->metaerror,
        previous_individuals[exp_m_->next_generation_reproducer_[indiv_id]]->fitness,
        previous_individuals[found_id]->fitness);



    for (int i = 0; i < neighborhood_size; i++) {
      if (i==8) {
        int v_x = indiv_index[i] / grid_height;
        int v_y = indiv_index[i] % grid_height;

        printf(
            "A-A-ERROR -- Individual %d (%d,%d): Metaerror (CPU/GPU) : %e/%e || Fitness (CPU/GPU) : %e/%e \n",
            indiv_index[i],v_x,v_y,
            exp_m_->world()->grid(v_x, v_y)->individual()->dist_to_target_by_feature(
                METABOLISM),
               previous_individuals[indiv_index[i]]->metaerror,
            exp_m_->world()->grid(v_x, v_y)->individual()->fitness(),
               previous_individuals[indiv_index[i]]->fitness);

        printf("ID CPU %lld SIMD %d -- PARENT ID CPU %d SIMD %d\n",
               exp_m_->world()->grid(v_x, v_y)->individual()->id(),
               previous_individuals[indiv_index[i]]->indiv_id,
               exp_m_->world()->grid(v_x, v_y)->individual()->parent_id_,
               previous_individuals[indiv_index[i]]->parent_id);

        // printf(
        //     "Nb RNA SIMD/CPU %ud/%ld Protein %ud/%ld\n",
        //     previous_individuals[indiv_index[i]]->metadata_->rna_count(),
        //     exp_m_->world()->grid(v_x, v_y)->individual()->rna_list().size(),
        //     previous_individuals[indiv_index[i]]->metadata_->proteins_count(),
        //     exp_m_->world()->grid(v_x, v_y)->individual()->protein_list().size());

        int idx = 0;

        // for (auto rna : exp_m_->world()->grid(v_x, v_y)->individual()->rna_list()) {
        //   printf("RNA CPU %d Start %d Stop %d Leading/Lagging %d Length %d\n", idx,
        //          rna->promoter_pos(), rna->last_transcribed_pos(), rna->strand(), rna->transcript_length());
        //   idx++;
        // }

        idx = 0;
        // for (idx = 0; idx < (int) (previous_individuals[indiv_index[i]]->metadata_->rna_count()); idx++) {
        //   printf("RNA SIMD %d Start %d Stop %d Leading/Lagging %d Length %d\n", idx,
        //       previous_individuals[indiv_index[i]]->metadata_->rnas(idx)->begin,
        //       previous_individuals[indiv_index[i]]->metadata_->rnas(idx)->end,
        //       previous_individuals[indiv_index[i]]->metadata_->rnas(idx)->leading_lagging,
        //       previous_individuals[indiv_index[i]]->metadata_->rnas(idx)->length);
        // }


      }

      printf("%d -- (Probs %e %e -- Fit Array %e %e -- Sum Fit %e %e\n",i,
             exp_m_->world()->grid(x,y)->probs[i],probs[i],
             exp_m_->world()->grid(x,y)->local_fit_array[i],local_fit_array[i],
             exp_m_->world()->grid(x,y)->sum_local_fit,sum_local_fit);
    }

    exit(-44);
  }



  delete [] local_fit_array;
  delete [] local_meta_array;
  delete [] probs;
}


void ExpManager_7::do_mutation(int indiv_id, double w_max, double selection_pressure) {
  // printf("%d -- DO MUTATION for %d\n",exp_m_->rank(),indiv_id);

  int32_t dna_length = -1;

  int32_t parent_id = -1;
  
  #ifdef HAVE_MPI
  bool is_at_border;
    int32_t local_x ;
      int32_t local_y ;
  char* remote_dna;
  int remote_parent_id;
  int32_t reprod_id_;
  int32_t remote_rank;

  // for (auto fetch_indiv : individual_to_fetch_at_border_) {
  //             dna_length = fetch_indiv.dna_length_;
  //             remote_dna = fetch_indiv.dna_;
  //             printf("%d -- %d -- Found remote DNA for %d (length %d)\n",exp_m_->rank(),
  //               indiv_id,
  //               fetch_indiv.id_,fetch_indiv.dna_length_);
            
  //         }
  #endif

  if (!exp_m_->check_simd()) {

    int32_t x = indiv_id / exp_m_->world()->height();
    int32_t y = indiv_id % exp_m_->world()->height();
    delete exp_m_->dna_mutator_array_[indiv_id];

    #ifdef HAVE_MPI
      int32_t reprod_global_x = exp_m_->next_generation_reproducer_[indiv_id] / global_grid_height_;
      int32_t reprod_global_y = exp_m_->next_generation_reproducer_[indiv_id] % global_grid_height_;

      int32_t reprod_pos_local_x = globalXtoLocalX(reprod_global_x);
      int32_t reprod_pos_local_y = globalYtoLocalY(reprod_global_y);
      reprod_id_ = reprod_pos_local_x * grid_height_ + reprod_pos_local_y;

      local_x = indiv_id / grid_height_;
      local_y = indiv_id % grid_height_;

      int32_t owner_rank = rankOf(reprod_global_x,reprod_global_y);
      
      is_at_border = owner_rank != exp_m_->rank();
    #endif

        #ifdef HAVE_MPI
        if (is_at_border) {
          bool found = false;
         
         #ifdef DEBUG_MPI
          printf("%d -- Looking for %d %d (local %d): Parent ID %d :: IFB %d\n",exp_m_->rank(),localXtoGlobalX(x),localYtoGlobalY(y),
                  indiv_id,exp_m_->next_generation_reproducer_[indiv_id], individual_to_fetch_at_border_.size());
         #endif

          // for (auto fetch_indiv : individual_to_fetch_at_border_) {
          //     printf("%d -- %d  (%d %d) -- Looking remote DNA for %d: Fetched Indiv %d (length %d) :: %d %d (%d %d)\n",exp_m_->rank(),
          //       indiv_id,local_x,local_y,exp_m_->next_generation_reproducer_[indiv_id],
          //       fetch_indiv.id_,fetch_indiv.dna_length_,
          //       fetch_indiv.x_,fetch_indiv.y_,fetch_indiv.local_x_,
          //       fetch_indiv.local_y_);
          // }

          for (auto fetch_indiv : individual_to_fetch_at_border_) {
              dna_length = fetch_indiv.dna_length_;
              remote_dna = fetch_indiv.dna_;
              remote_rank = fetch_indiv.rank_;
              remote_parent_id = fetch_indiv.id_;
              
              // printf("%d -- %d  (%d %d) -- Looking remote DNA for %d: Fetched Indiv %d (length %d)\n",exp_m_->rank(),
              //   indiv_id,local_x,local_y,exp_m_->next_generation_reproducer_[indiv_id],
              //   fetch_indiv.id_,fetch_indiv.dna_length_);

            if (exp_m_->next_generation_reproducer_[indiv_id] == fetch_indiv.id_ && fetch_indiv.fetched_) {
              // dna_length = fetch_indiv.dna_length_;
              // remote_dna = fetch_indiv.dna_;
              #ifdef DEBUG_MPI
              printf("%d -- %d -- FOUND (%d %d :: %d) remote DNA for %d (length %d) : %d %d (%d %d) :: %d \n",exp_m_->rank(),
                indiv_id,local_x,local_y,exp_m_->next_generation_reproducer_[indiv_id],
                fetch_indiv.id_,fetch_indiv.dna_length_,fetch_indiv.x_,fetch_indiv.y_,fetch_indiv.local_x_,
                fetch_indiv.local_y_,fetch_indiv.fetched_);
              // printf("%d -- Found remote DNA for %d (length %d)\n",exp_m_->rank(),indiv_id,fetch_indiv.dna_length_);
              #endif
              found = true;
              break;
            }
          }

          if (!found) {printf("%d -- Error indiv %d not found\n",exp_m_->rank(),indiv_id);exit(-1);}
        } else {
        dna_length = previous_individuals[reprod_id_]->dna_->length();
        }
        #else
        dna_length = previous_individuals[exp_m_->next_generation_reproducer_[indiv_id]]->dna_->length();
        #endif

    // printf("Create DNA mutator %d\n",indiv_id);
    exp_m_->dna_mutator_array_[indiv_id] = new DnaMutator(
      #ifdef HAVE_MPI
        exp_m_->world()->grid(local_x, local_y)->mut_prng(),
      #else
        exp_m_->world()->grid(x, y)->mut_prng(),
      #endif
        dna_length,
        exp_m_->exp_s()->mut_params()->duplication_rate(),
        exp_m_->exp_s()->mut_params()->deletion_rate(),
        exp_m_->exp_s()->mut_params()->translocation_rate(),
        exp_m_->exp_s()->mut_params()->inversion_rate(),
        exp_m_->exp_s()->mut_params()->point_mutation_rate(),
        exp_m_->exp_s()->mut_params()->small_insertion_rate(),
        exp_m_->exp_s()->mut_params()->small_deletion_rate(),
        exp_m_->exp_s()->mut_params()->max_indel_size(),
        exp_m_->exp_s()->min_genome_length(),
        exp_m_->exp_s()->max_genome_length(), indiv_id,
      #ifdef HAVE_MPI
        local_x,local_y
      #else
        x,y
      #endif  
        );
    exp_m_->dna_mutator_array_[indiv_id]->generate_mutations();
  }

  // {
  //     int32_t global_grid_width = exp_m_->exp_s()->global_grid_width();
  //     int32_t global_grid_height = exp_m_->exp_s()->global_grid_height();

  //     int32_t grid_width = exp_m_->grid_width();
  //     int32_t grid_height = exp_m_->grid_height();

  //     int32_t rank_x = exp_m_->exp_s()->rank_width();
  //     int32_t rank_y = exp_m_->exp_s()->rank_height();

  //     int32_t local_rank_x = exp_m_->rank() / rank_y;
  //     int32_t local_rank_y = exp_m_->rank() % rank_y;

  //     int32_t x_offset = (local_rank_x) * grid_width; 
  //     int32_t y_offset = (local_rank_y) * grid_height;

  //     int32_t i = (local_x+x_offset+global_grid_width) % global_grid_width;
  //     int32_t j = (local_y+y_offset+global_grid_height) % global_grid_height;
      
  //     #ifdef DEBUG_MPI
  //     printf("%d -- Parent %d => Child %d (Mutate %d)\n",exp_m_->rank(),i*global_grid_height+j,
  //         exp_m_->next_generation_reproducer_[indiv_id],
  //         exp_m_->dna_mutator_array_[indiv_id]->hasMutate());
  //     #endif

  // }

  // printf("%d -- Read Mutator %d : %d (%d %d) -- BORDER %d\n", exp_m_->rank(), indiv_id, 
  //     exp_m_->dna_mutator_array_[indiv_id]->hasMutate(),
  //     exp_m_->dna_mutator_array_[indiv_id]->nb_mut_,exp_m_->dna_mutator_array_[indiv_id]->nb_rear_,is_at_border);



  if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate()) {
    if (ExpManager_7::standalone() && !exp_m_->check_simd())
      exp_m_->dna_mutator_array_[indiv_id]->generate_all_mutations(exp_m_->dna_mutator_array_[indiv_id]->length_);
    // #ifdef HAVE_MPI
    // if (!is_at_border) {
    // #endif
    
    #pragma omp critical
    {
    mutant_list_.push_back(indiv_id);
    }
// printf("MUTATED %d\n",indiv_id);
// #ifdef HAVE_MPI
//     }
// #endif

  } else {

#ifdef HAVE_MPI
    int32_t parent_id = reprod_id_;
#else
    int32_t parent_id = exp_m_->next_generation_reproducer_[indiv_id];
#endif

    #ifdef __REGUL
    bool first_to_add = false;
    #pragma omp critical
    {
      first_to_add = previous_individuals[parent_id]->first_to_add;
      previous_individuals[parent_id]->first_to_add = false;
    }
    
    if (first_to_add)
        previous_individuals[parent_id]->added_id = indiv_id;

    
    if (first_to_add && phenotypic_target_handler_->hasChanged_) {
      #pragma omp critical
      {
      mutant_list_.push_back(indiv_id);
      }
    }
    #endif
      #ifdef HAVE_MPI
  if (is_at_border) {
    exp_m_->dna_mutator_array_[indiv_id]->setIsAtBorder(true);
    // printf("%d -- COPYING %d with parent %d\n",exp_m_->rank(),indiv_id,parent_id);
    #pragma omp critical
    {
    mutant_list_.push_back(indiv_id);
    }
  } else {
#endif
#pragma omp atomic
    nb_clones_++;


#ifdef DEBUG_MPI
    printf("%d -- Linking %d with parent %d\n",exp_m_->rank(),indiv_id,parent_id);
    #endif
    current_individuals[indiv_id] = previous_individuals[parent_id];


    #pragma omp atomic
    current_individuals[indiv_id]->usage_count_++;


    #pragma omp critical
    {
    current_individuals[indiv_id]->last_id = indiv_id;
    current_individuals[indiv_id]->indiv_id = indiv_id;
    }

// #ifdef __REGUL
//     // if (!first_to_add)
//     // if (indiv_id == 205)
//       printf("CLONED %d (p %d) : %d :: E %d :: %p == %p == %p\n",indiv_id,parent_id,previous_individuals[parent_id]->added_id,
//               phenotypic_target_handler_->hasChanged_,
//               current_individuals[indiv_id],previous_individuals[parent_id],current_individuals[current_individuals[indiv_id]->added_id]);
// #endif

    if (ExpManager_7::standalone() && exp_m_->record_tree()) {
      int x = indiv_id / exp_m_->world()->height();
      int y = indiv_id % exp_m_->world()->height();

      #ifdef HAVE_MPI
      int32_t i = localXtoGlobalX(x);
      int32_t j = localYtoGlobalY(y);
      #endif

      #ifdef DEBUG_MPI
      int32_t reprod_global_x = exp_m_->next_generation_reproducer_[indiv_id] / global_grid_height_;
      int32_t reprod_global_y = exp_m_->next_generation_reproducer_[indiv_id] % global_grid_height_;

      int32_t reprod_pos_local_x = globalXtoLocalX(reprod_global_x);
      int32_t reprod_pos_local_y = globalYtoLocalY(reprod_global_y);

      printf("%d -- C -- New Indiv from local clone : Local ID %d Parent ID %d :: [%d => %d] :: %d %d => %d %d\n",// Parent ID %d :: %d // PSize %d\n",
                exp_m_->rank(),indiv_id,//local_x,local_y,i*global_grid_height_+j,i,j,
                exp_m_->next_generation_reproducer_[indiv_id],
                current_individuals[indiv_id]->usage_count_,
                previous_individuals[reprod_pos_local_x*grid_height_+reprod_pos_local_y]->usage_count_,
                reprod_global_x,reprod_global_y,
                reprod_pos_local_x,reprod_pos_local_y);
      #endif
      NewIndivEvent *eindiv = new NewIndivEvent(
          #ifdef HAVE_MPI
            current_individuals[indiv_id],(Individual_7*)nullptr,
            local_x,local_y,i*global_grid_height_+j,
          #else
            current_individuals[indiv_id],
            previous_individuals[exp_m_->next_generation_reproducer_[indiv_id]],
            x, y,indiv_id,
          #endif       
            exp_m_->next_generation_reproducer_[indiv_id]);
      exp_m_->tree()->update_new_indiv(eindiv);
      delete eindiv;
    }

#ifdef WITH_PERF_TRACES_PER_INDIV
    apply_mutation[indiv_id] = -1;
#endif
#ifdef HAVE_MPI
    }
#endif
  }

  // printf("%d -- Mutation for %d IS DONE\n",exp_m_->rank(),indiv_id);
}

void ExpManager_7::apply_mutations(int indiv_id, double w_max, double selection_pressure) {

#ifdef WITH_PERF_TRACES_PER_INDIV
    allocate_individual_start_[indiv_id] = std::chrono::steady_clock::now().time_since_epoch().count();
#endif

  int32_t dna_length = -1;

  int32_t parent_id = -1;
      #ifdef HAVE_MPI
  bool is_at_border;
    int32_t local_x ;
      int32_t local_y ;
  char* remote_dna;
  int remote_parent_id;
  int32_t reprod_id_;
  int32_t remote_rank;

    int32_t lead_prom_size;

    int32_t* lead_prom_pos;
    int8_t* lead_prom_error;
            
    int32_t lag_prom_size;
            
    int32_t* lag_prom_pos;
    int8_t* lag_prom_error;

      int32_t reprod_global_x = exp_m_->next_generation_reproducer_[indiv_id] / global_grid_height_;
      int32_t reprod_global_y = exp_m_->next_generation_reproducer_[indiv_id] % global_grid_height_;

      int32_t reprod_pos_local_x = globalXtoLocalX(reprod_global_x);
      int32_t reprod_pos_local_y = globalYtoLocalY(reprod_global_y);
      reprod_id_ = reprod_pos_local_x * grid_height_ + reprod_pos_local_y;

      local_x = indiv_id / grid_height_;
      local_y = indiv_id % grid_height_;

      int32_t owner_rank = rankOf(reprod_global_x,reprod_global_y);
      
      is_at_border = owner_rank != exp_m_->rank();

          
      if (!is_at_border) {
    #endif

      current_individuals[indiv_id] =
        new Individual_7(exp_m_, 
        #ifdef HAVE_MPI
        previous_individuals[reprod_id_],
        #else
         previous_individuals[exp_m_->next_generation_reproducer_[indiv_id]],
        #endif
        dna_factory_,fuzzy_factory_);
    current_individuals[indiv_id]->global_id = AeTime::time()*1024+indiv_id;
    current_individuals[indiv_id]->indiv_id = indiv_id;
    current_individuals[indiv_id]->last_id = indiv_id;
    current_individuals[indiv_id]->parent_id =
        exp_m_->next_generation_reproducer_[indiv_id];
        current_individuals[indiv_id]->dna_->dna_factory_ = dna_factory_;

    if (ExpManager_7::standalone() && exp_m_->record_tree()) {
      int x = indiv_id / exp_m_->world()->height();
      int y = indiv_id % exp_m_->world()->height();

      #ifdef HAVE_MPI
      int32_t i = localXtoGlobalX(local_x);
      int32_t j = localYtoGlobalY(local_y);
      #endif

      #ifdef DEBUG_MPI
      //  (lx + x_offset_) % global_grid_width_;
      printf("%d -- B -- New Indiv from local mutate : Local ID %d (%d x %d) Global ID %d (%d x %d) Parent ID %d\n",// Parent ID %d :: %d // PSize %d\n",
                exp_m_->rank(),indiv_id,local_x,local_y,i*global_grid_height_+j,i,j,
                exp_m_->next_generation_reproducer_[indiv_id]); //exp_m_->next_generation_reproducer_[indiv_id],
                //current_individuals[indiv_id]->dna_->length_,previous_individuals[reprod_id_]->dna_->length_);
      #endif
      NewIndivEvent *eindiv = new NewIndivEvent(
          #ifdef HAVE_MPI
            current_individuals[indiv_id],(Individual_7*)nullptr,
            local_x,local_y,i*global_grid_height_+j,
          #else
            current_individuals[indiv_id],
            previous_individuals[exp_m_->next_generation_reproducer_[indiv_id]],
            x, y,indiv_id,
          #endif                                      
            exp_m_->next_generation_reproducer_[indiv_id]);

      exp_m_->tree()->update_new_indiv(eindiv);
      delete eindiv;
    }
    #ifdef HAVE_MPI
      } else {
          bool found = false;
      

          for (auto fetch_indiv : individual_to_fetch_at_border_) {
              
              // printf("%d -- %d  (%d %d) -- Looking remote DNA for %d: Fetched Indiv %d (length %d)\n",exp_m_->rank(),
              //   indiv_id,local_x,local_y,exp_m_->next_generation_reproducer_[indiv_id],
              //   fetch_indiv.id_,fetch_indiv.dna_length_);

            if (exp_m_->next_generation_reproducer_[indiv_id] == fetch_indiv.id_ && fetch_indiv.fetched_) {
              dna_length = fetch_indiv.dna_length_;
              remote_dna = fetch_indiv.dna_;
              remote_rank = fetch_indiv.rank_;
              remote_parent_id = fetch_indiv.id_;
              
              lead_prom_size = fetch_indiv.lead_prom_size;
              lead_prom_pos = fetch_indiv.lead_prom_pos;
              lead_prom_error = fetch_indiv.lead_prom_error;
            
              lag_prom_size = fetch_indiv.lag_prom_size;
              lag_prom_pos = fetch_indiv.lag_prom_pos;
              lag_prom_error = fetch_indiv.lag_prom_error;
              // dna_length = fetch_indiv.dna_length_;
              // remote_dna = fetch_indiv.dna_;
              #ifdef DEBUG_MPI

      int32_t i = localXtoGlobalX(local_x);
      int32_t j = localYtoGlobalY(local_y);
              printf("T %d -- %d -- Local %d Global %d -- FOUND (%d %d :: %d) remote DNA for %d (length %d) : %d %d (%d %d) :: %d :: %d %d\n",
              AeTime::time(),exp_m_->rank(),
                indiv_id,i*global_grid_width_+j,local_x,local_y,exp_m_->next_generation_reproducer_[indiv_id],
                fetch_indiv.id_,fetch_indiv.dna_length_,fetch_indiv.x_,fetch_indiv.y_,fetch_indiv.local_x_,
                fetch_indiv.local_y_,fetch_indiv.fetched_,fetch_indiv.lead_prom_size,fetch_indiv.lag_prom_size);
              // printf("%d -- Found remote DNA for %d (length %d)\n",exp_m_->rank(),indiv_id,fetch_indiv.dna_length_);
              #endif
              found = true;
              break;
            }
          }

          if (!found) {printf("%d -- Error indiv %d not found\n",exp_m_->rank(),indiv_id);exit(-1);}

  #ifdef DEBUG_MPI
      printf("%d -- %d => Creating from remote DNA %d\n",exp_m_->rank(),indiv_id,dna_length);
      #endif

      current_individuals[indiv_id] = new Individual_7(exp_m_, exp_m_->best_indiv()->w_max(), remote_dna, dna_length, 
                                                        dna_factory_,fuzzy_factory_,lead_prom_pos,lead_prom_error,lead_prom_size,
                                                        lag_prom_pos,lag_prom_error,lag_prom_size);
      // current_individuals[indiv_id]->dna_ = dna_factory_->get_dna(dna_length);
      // current_individuals[indiv_id]->dna_->set_indiv(remote_dna,dna_length,dna_factory_, current_individuals[indiv_id]);

  #ifdef DEBUG_MPI
      printf("%d -- %d => DNA Individual %p (%p) \n",exp_m_->rank(),indiv_id,current_individuals[indiv_id]->dna_->indiv_,current_individuals[indiv_id]);
  #endif

      current_individuals[indiv_id]->global_id = AeTime::time()*1024+indiv_id;
      current_individuals[indiv_id]->indiv_id = indiv_id;
      current_individuals[indiv_id]->last_id = indiv_id;
      current_individuals[indiv_id]->parent_id = remote_parent_id;
      current_individuals[indiv_id]->is_at_border_ = true;

#ifdef DEBUG_MPI
      opt_prom_compute_RNA(current_individuals[indiv_id]);
      start_protein(current_individuals[indiv_id]);
      compute_protein(current_individuals[indiv_id]);
      translate_protein(current_individuals[indiv_id], exp_m_->w_max_);
      compute_phenotype(current_individuals[indiv_id]);
      compute_fitness(current_individuals[indiv_id], exp_m_->selection_pressure());

      int32_t i = localXtoGlobalX(local_x);
      int32_t j = localYtoGlobalY(local_y);
      printf("T %d -- %d -- Individual %d (G %d): must be computed from scratch brder %d (%d %d):: %e\n",AeTime::time(),exp_m_->rank(),
                indiv_id,
                i*global_grid_width_+j,
                current_individuals[indiv_id]->metadata_->promoter_count(),
                lead_prom_size,lag_prom_size,
                current_individuals[indiv_id]->fitness);
      #ifdef HAVE_MPI

                FILE *fp;
      char filename[20];
      sprintf(filename,"test_log_%d",exp_m_->rank());
      fp = fopen(filename, "a");
#else

      char filename[20];
      sprintf(filename,"test_log_0");
      fp = fopen(filename, "a");
#endif
      fprintf(fp,"%d -- %d -- BEFORE_MUT -- Promoter list LEADING : ",AeTime::time(),indiv_id);
      for (auto prot : ((List_Metadata*)current_individuals[indiv_id]->metadata_)->promoters_list_[LEADING]) {
        fprintf(fp,"%d ",prot.pos);
      }
      fprintf(fp,"\n");

      fprintf(fp,"%d -- %d -- BEFORE_MUT -- Promoter list LAGGING : ",AeTime::time(),indiv_id);
      for (auto prot : ((List_Metadata*)current_individuals[indiv_id]->metadata_)->promoters_list_[LAGGING]) {
        fprintf(fp,"%d ",prot.pos);
      }
      fprintf(fp,"\n");
      fclose(fp);
      #endif
      // evaluate(current_individuals[indiv_id],w_max,selection_pressure);
      
      if (ExpManager_7::standalone() && exp_m_->record_tree()) {
        int32_t i = localXtoGlobalX(local_x);
        int32_t j = localYtoGlobalY(local_y);
        #ifdef DEBUG_MPI
        printf("%d -- A -- New Indiv from Remote : Local ID %d (%d x %d) Global ID %d (%d x %d) Parent ID %d\n",
                  exp_m_->rank(),indiv_id,local_x,local_y,i*global_grid_height_+j,i,j,remote_parent_id
                  );
        #endif
        NewIndivEvent *eindiv = new NewIndivEvent(
            current_individuals[indiv_id],
            nullptr,
                                                  local_x, local_y,
                                                  i*global_grid_height_+j,remote_parent_id,true,remote_rank);

        exp_m_->tree()->update_new_indiv(eindiv);
        delete eindiv;
      }
    ((List_Metadata*)current_individuals[indiv_id]->metadata_)->clean_remote();

      }
      #endif

#ifdef WITH_PERF_TRACES_PER_INDIV
    auto t_start = std::chrono::steady_clock::now();
    allocate_individual_stop_[indiv_id] = t_start.time_since_epoch().count();
    apply_mutation_start_[indiv_id] = t_start.time_since_epoch().count();
#endif
((List_Metadata*)current_individuals[indiv_id]->metadata_)->indiv_ = current_individuals[indiv_id];

// if (indiv_id == 6) {
//           printf("BM -- Promoters LEADING : ");
//         for (auto& prom : ((List_Metadata*)current_individuals[indiv_id]->metadata_)->promoters_list_[LEADING]) {
//           printf("%d ",prom.pos);
//         }
//         printf("\n");
//         printf("BM -- Promoters LAGGING : ");
//         for (auto& prom : ((List_Metadata*)current_individuals[indiv_id]->metadata_)->promoters_list_[LAGGING]) {
//           printf("%d ",prom.pos);
//         }
//         printf("\n");
// }
    if (ExpManager_7::standalone() && !exp_m_->check_simd())
      current_individuals[indiv_id]->dna_->apply_mutations_standalone();
    else
      current_individuals[indiv_id]->dna_->apply_mutations();

// if (indiv_id == 6) {
//           printf("BM -- Promoters LEADING : ");
//         for (auto& prom : ((List_Metadata*)current_individuals[indiv_id]->metadata_)->promoters_list_[LEADING]) {
//           printf("%d ",prom.pos);
//         }
//         printf("\n");
//         printf("BM -- Promoters LAGGING : ");
//         for (auto& prom : ((List_Metadata*)current_individuals[indiv_id]->metadata_)->promoters_list_[LAGGING]) {
//           printf("%d ",prom.pos);
//         }
//         printf("\n");
// }

#ifdef DEBUG_MPI
#ifdef HAVE_MPI
      int32_t i = localXtoGlobalX(local_x);
      int32_t j = localYtoGlobalY(local_y);
          printf("T %d -- %d -- Individual %d (G %d): has BEEN computed from scratch brder %d :: %e\n",AeTime::time(),exp_m_->rank(),
                indiv_id,
                i*global_grid_width_+j,
                current_individuals[indiv_id]->metadata_->promoter_count(),
                current_individuals[indiv_id]->fitness);
                                FILE *fp;

                      char filename[20];
      sprintf(filename,"test_log_%d",exp_m_->rank());
      fp = fopen(filename, "a");
#else
      int32_t i = indiv_id / exp_m_->grid_height();
      int32_t j = indiv_id % exp_m_->grid_height();

          printf("T %d -- Individual %d (G %d): has BEEN computed from scratch brder %d :: %e\n",AeTime::time(),
                indiv_id,
                i*grid_width_+j,
                current_individuals[indiv_id]->metadata_->promoter_count(),
                current_individuals[indiv_id]->fitness);
                FILE *fp;

                      char filename[20];
      sprintf(filename,"test_log_0");
      fp = fopen(filename, "a");
#endif

fprintf(fp,"%d -- %d -- AFTER_MUT -- Promoter list LEADING : ",AeTime::time(),indiv_id);
      for (auto prot : ((List_Metadata*)current_individuals[indiv_id]->metadata_)->promoters_list_[LEADING]) {
        fprintf(fp,"%d ",prot.pos);
      }
      fprintf(fp,"\n");

      fprintf(fp,"%d -- %d -- AFTER_MUT -- Promoter list LAGGING : ",AeTime::time(),indiv_id);
      for (auto prot : ((List_Metadata*)current_individuals[indiv_id]->metadata_)->promoters_list_[LAGGING]) {
        fprintf(fp,"%d ",prot.pos);
      }
      fprintf(fp,"\n");
      fclose(fp);
      #ifdef DEBUG_MPI
    {
      int32_t x = indiv_id / exp_m_->grid_height();
      int32_t y = indiv_id % exp_m_->grid_height();

      int32_t global_x = localXtoGlobalX(x);
      int32_t global_y = localYtoGlobalY(y);
      int32_t global_id = global_x * global_grid_height_ + global_y;
    printf("%d -- %d -- %d -- Dna length %d\n",AeTime::time(),exp_m_->rank(),global_id,current_individuals[indiv_id]->dna_->length_);
  }
    #endif
    #endif
#ifdef WITH_PERF_TRACES_PER_INDIV
    auto t_end = std::chrono::steady_clock::now();
    apply_mutation_stop_[indiv_id] = t_end.time_since_epoch().count();
    apply_mutation[indiv_id] = t_end.time_since_epoch().count() - t_start.time_since_epoch().count();
#endif

}

ExpManager_7::~ExpManager_7() {

  printf("Destroy SIMD Controller\n");

  delete stats_best;
  delete stats_mean;

  /* No need to delete current_individuals, at the end of a generation all the
   * element of the table are nullptr
   */
  #pragma omp single
  {
    int cpt = 0;
  for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
      if (previous_individuals[indiv_id]->usage_count_ > 1)
        previous_individuals[indiv_id]->usage_count_--;
      else {
          // for (int rn = 0; rn < previous_individuals[indiv_id]->metadata_->rna_count(); rn++) {
          //   delete previous_individuals[indiv_id]->metadata_->rnas(rn);
          // }

          // previous_individuals[indiv_id]->metadata_->rnas_clear();
          // for (int rn = 0; rn < previous_individuals[indiv_id]->metadata_->proteins_count(); rn++) {
          //   delete previous_individuals[indiv_id]->metadata_->proteins(rn);
          // }
          // previous_individuals[indiv_id]->metadata_->proteins_clear();
          previous_individuals[indiv_id]->metadata_->indiv_ = previous_individuals[indiv_id];
          delete previous_individuals[indiv_id];
          previous_individuals[indiv_id] = nullptr;
          cpt++;
      }

    delete exp_m_->dna_mutator_array_[indiv_id];
  }
    // printf("CPT %d\n",cpt);
  }

  // #ifdef __REGUL
  // #pragma omp single
  // {

  // int nb_not_clean = 0;
  //   for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
  //     if (current_individuals[indiv_id]!=nullptr) {
  //       nb_not_clean++;
  //       delete current_individuals[indiv_id];
  //     }
  //   }

  //   printf("Not clean %d\n",nb_not_clean);
  // }
  // #endif

  delete[] previous_individuals;
  delete[] current_individuals;


#ifdef __REGUL
  delete phenotypic_target_handler_;
#else
  fuzzy_factory_->give_back(target);
#endif

  delete dna_factory_;
  delete fuzzy_factory_;
}

void ExpManager_7::start_stop_RNA(int indiv_id) {
  start_stop_RNA(current_individuals[indiv_id]);
}

void ExpManager_7::start_stop_RNA(Individual_7* indiv) {
  for (int dna_pos = 0; dna_pos < indiv->dna_->length(); dna_pos++) {
    int len = indiv->dna_->length();
    if (len >= PROM_SIZE) {
      #ifdef BASE_2
      int prom_dist_leading[26];
      int prom_dist_lagging[26];

      int term_dist_leading[4];
      int term_dist_lagging[4];

      for (int motif_id = 0; motif_id < 52; motif_id++) {
        if (motif_id >= 26 && motif_id < 48) {
          // LAGGING
          int t_motif_id = motif_id - 26;
          prom_dist_lagging[t_motif_id] =
              PROM_SEQ_LAG[t_motif_id] ==
                      indiv->dna_->data_[
                  dna_pos - t_motif_id < 0 ? len +
                                             dna_pos -
                                             t_motif_id :
                  dna_pos - t_motif_id]
              ? 0 : 1;
        } else if (motif_id < 22) {
          // LEADING
          prom_dist_leading[motif_id] =
              PROM_SEQ_LEAD[motif_id] ==
                      indiv->dna_->data_[
                  dna_pos + motif_id >= len ? dna_pos +
                                              motif_id -
                                              len
                                            : dna_pos +
                                              motif_id]
              ? 0
              : 1;
        } else if (motif_id >= 22 && motif_id < 26) {
          int t_motif_id = motif_id - 22;
          // LEADING
          term_dist_leading[t_motif_id] =
              indiv->dna_->data_[
                  dna_pos + t_motif_id >= len ? dna_pos +
                                                t_motif_id -
                                                len :
                  dna_pos + t_motif_id] !=
                      indiv->dna_->data_[
                  dna_pos - t_motif_id + 10 >= len ?
                  dna_pos - t_motif_id + 10 - len :
                  dna_pos -
                  t_motif_id +
                  10] ? 1
                      : 0;
        } else {
          int t_motif_id = motif_id - 48;
          term_dist_lagging[t_motif_id] =
              indiv->dna_->data_[
                  dna_pos - t_motif_id < 0 ? dna_pos -
                                             t_motif_id +
                                             len
                                           : dna_pos -
                                             t_motif_id] !=
              indiv->dna_->data_[
                  dna_pos + t_motif_id - 10 < 0 ? dna_pos +
                                                  t_motif_id -
                                                  10 + len
                                                :
                  dna_pos + t_motif_id - 10] ? 1 : 0;
        }
      }


      int i_prom_dist_leading = prom_dist_leading[0] +
                      prom_dist_leading[1] +
                      prom_dist_leading[2] +
                      prom_dist_leading[3] +
                      prom_dist_leading[4] +
                      prom_dist_leading[5] +
                      prom_dist_leading[6] +
                      prom_dist_leading[7] +
                      prom_dist_leading[8] +
                      prom_dist_leading[9] +
                      prom_dist_leading[10] +
                      prom_dist_leading[11] +
                      prom_dist_leading[12] +
                      prom_dist_leading[13] +
                      prom_dist_leading[14] +
                      prom_dist_leading[15] +
                      prom_dist_leading[16] +
                      prom_dist_leading[17] +
                      prom_dist_leading[18] +
                      prom_dist_leading[19] +
                      prom_dist_leading[20] +
                      prom_dist_leading[21];

     int i_prom_dist_lagging = prom_dist_lagging[0] +
                     prom_dist_lagging[1] +
                     prom_dist_lagging[2] +
                     prom_dist_lagging[3] +
                     prom_dist_lagging[4] +
                     prom_dist_lagging[5] +
                     prom_dist_lagging[6] +
                     prom_dist_lagging[7] +
                     prom_dist_lagging[8] +
                     prom_dist_lagging[9] +
                     prom_dist_lagging[10] +
                     prom_dist_lagging[11] +
                     prom_dist_lagging[12] +
                     prom_dist_lagging[13] +
                     prom_dist_lagging[14] +
                     prom_dist_lagging[15] +
                     prom_dist_lagging[16] +
                     prom_dist_lagging[17] +
                     prom_dist_lagging[18] +
                     prom_dist_lagging[19] +
                     prom_dist_lagging[20] +
                     prom_dist_lagging[21];

      int i_term_leading = term_dist_leading[0] +
                           term_dist_leading[1] +
                           term_dist_leading[2] +
                           term_dist_leading[3];


      int i_term_lagging = term_dist_lagging[0] +
                          term_dist_lagging[1] +
                          term_dist_lagging[2] +
                          term_dist_lagging[3];

      #elif BASE_4
      int i_prom_dist_leading = indiv->metadata_->is_promoter_leading(dna_pos);
      int i_prom_dist_lagging = indiv->metadata_->is_promoter_lagging(dna_pos);
      int i_term_leading = indiv->metadata_->is_terminator_leading(dna_pos);
      int i_term_lagging = indiv->metadata_->is_terminator_lagging(dna_pos);
      #endif

      if (i_prom_dist_leading <= PROM_MAX_DIFF) {
        PromoterStruct* nprom = new PromoterStruct(dna_pos, i_prom_dist_leading,
                                                   true);
        int prom_idx;

        prom_idx = indiv->metadata_->promoter_count();
        indiv->metadata_->set_promoters_count(
            indiv->metadata_->promoter_count()+ 1);

        indiv->metadata_->promoter_add(prom_idx, nprom);

        delete nprom;
      }


      #ifdef BASE_2
      if (i_term_leading == 4) {
      #elif BASE_4
      if (i_term_leading) {
      #endif
        indiv->metadata_->terminator_add(LEADING, dna_pos);
      }


      if (i_prom_dist_lagging <= PROM_MAX_DIFF) {
        PromoterStruct* nprom = new PromoterStruct(dna_pos, i_prom_dist_lagging,
                                                   false);
        int prom_idx;
        prom_idx = indiv->metadata_->promoter_count();
        indiv->metadata_->set_promoters_count(
            indiv->metadata_->promoter_count() + 1);

        indiv->metadata_->promoter_add(prom_idx, nprom);
        delete nprom;
      }




      #ifdef BASE_2
      if (i_term_lagging == 4) {
      #elif BASE_4
      if (i_term_lagging) {
      #endif        
        indiv->metadata_->terminator_add(LAGGING, dna_pos);
      }
    }
  }
}

void ExpManager_7::opt_prom_compute_RNA(int indiv_id) {
  opt_prom_compute_RNA(current_individuals[indiv_id]);
}

void ExpManager_7::opt_prom_compute_RNA(Individual_7* indiv) {
    #ifdef WITH_PERF_TRACES_PER_INDIV
      auto runtime_start = std::chrono::steady_clock::now();
      compute_rna_start_[indiv->indiv_id] = runtime_start.time_since_epoch().count();
    #endif

    indiv->metadata_->proteins_clear();
    indiv->metadata_->rnas_clear();
    indiv->metadata_->terminators_clear();
 
    indiv->metadata_->rnas_resize(
        indiv->metadata_->promoter_count());

    indiv->metadata_->promoter_begin();

    #ifdef BASE_4
    int8_t term_size = TERM_SIZE + exp_m_->exp_s()->terminator_polya_sequence_length();
    #endif

    // ANNOTATE_SITE_BEGIN(mainloop);
    for (int prom_idx = 0;
         prom_idx < (int)indiv->metadata_->promoter_count(); prom_idx++) {
      PromoterStruct*prom =
          indiv->metadata_->promoter_next();

      if (prom != nullptr) {
        // if (indiv_id == 279 && AeTime::time() == 95)
          
          // printf("%d -- RNA %d => %p (%d %d)\n",indiv->indiv_id,prom->pos,prom->rna,prom->to_compute_end_ ,
          // ExpManager_7::compute_diff_rnas );
        // if (prom->to_compute_end_ || !ExpManager_7::compute_diff_rnas ) {
          #ifdef WITH_OPTIMIZE_DIFF_SEARCH
            if (prom->to_compute_end_ ) {
          #else
            if (true) {
          #endif
        Dna_7*dna = indiv->dna_;
        int32_t dna_length = dna->length();
        int prom_pos;
        bool lead_lag;
        double prom_error;

        prom_pos = prom->pos;
        lead_lag = prom->leading_or_lagging;
        prom_error = fabs(
            ((float) prom->error));

          int32_t length = dna->length_;
          char *dna_data = dna->data_;
          #if defined(__INTEL_COMPILER)
          __declspec(align(dna_data));
          #elif defined(__INTEL_LLVM_COMPILER)
          void* vec_r = __builtin_assume_aligned(dna_data,64);
          #endif

        if (lead_lag) {
          /* Search for terminators */
          int cur_pos =
              Utils::mod(prom_pos + PROM_SIZE,dna_length);
          int start_pos = cur_pos;

          bool terminator_found = false;
          bool no_terminator = false;
          int8_t term_dist_leading = 0;

          int loop_size = 0;


          while (!terminator_found) {
            #ifdef BASE_2
            // ANNOTATE_SITE_BEGIN(leading_motif_loop);
            if (cur_pos + 10 < length && cur_pos + 10 >= 0) {
              #pragma omp simd
              for (int32_t t_motif_id = 0; t_motif_id < 4; t_motif_id++) {
                int32_t pos_begin = cur_pos + t_motif_id;
                int32_t pos_end   = cur_pos - t_motif_id + 10;

                // ANNOTATE_ITERATION_TASK(searchLeading);
                term_dist_leading +=
                    dna_data[pos_begin] !=
                    dna_data[pos_end]
                    ? 1 : 0;
              }
            } else {
              for (int32_t t_motif_id = 0; t_motif_id < 4; t_motif_id++) {
                int32_t pos_begin = cur_pos + t_motif_id;
                int32_t pos_end   = cur_pos - t_motif_id + 10;
                
                while (pos_begin < 0)       pos_begin += length;
                while (pos_begin >= length) pos_begin -= length;

                while (pos_end < 0)       pos_end += length;
                while (pos_end >= length) pos_end -= length;

                // ANNOTATE_ITERATION_TASK(searchLeading);
                term_dist_leading +=
                    dna_data[pos_begin] !=
                    dna_data[pos_end]
                    ? 1 : 0;
              }
              // for (int t_motif_id = 0; t_motif_id < 4; t_motif_id++) {

              //   term_dist_leading +=
              //       dna->get_lead(cur_pos + t_motif_id) !=
              //       dna->get_lead(cur_pos - t_motif_id + 10)
              //       ? 1 : 0;
              // }
            }
            
            // ANNOTATE_SITE_END(leading_motif_loop);

            if (term_dist_leading == 4) {
              terminator_found = true;
            }
            else {
              cur_pos = Utils::mod(cur_pos + 1,dna_length);

              term_dist_leading = 0;
              if (cur_pos == start_pos) {
                no_terminator = true;
                terminator_found = true;
              }
            }


            #elif BASE_4
            terminator_found = indiv->metadata_->is_terminator_leading(cur_pos);

            // if (indiv->indiv_id == 817 && terminator_found) printf("SIMD/Found a term at %d\n",cur_pos);

            if(!terminator_found) {
              cur_pos = cur_pos + 1 >= dna_length ?
                        cur_pos + 1 - dna_length : cur_pos + 1;

              if (cur_pos == start_pos) {
                no_terminator = true;
                terminator_found = true;
              }
            }
            #endif

            loop_size++;

          }

          if (!no_terminator) {
            #ifdef BASE_2  
            int32_t rna_end =
                cur_pos;
            int32_t rna_length = 0;

            rna_length = loop_size + TERM_SIZE -1;

            if (rna_length > 0) {


              int glob_rna_idx = -1;
              glob_rna_idx =
                  indiv->metadata_->rna_count_++;
              indiv->metadata_->rna_add(glob_rna_idx,
                                                                 prom_pos,
                                                                 rna_end,
                                                                 !lead_lag,
                                                                 1.0 -
                                                                 prom_error /
                                                                 5.0, rna_length,prom);

            }
            #elif BASE_4
            int32_t rna_end = Utils::mod(cur_pos + term_size, dna_length);
            int32_t rna_length = loop_size-1 + term_size;

            if (rna_length >= term_size) {  // check that terminator is not overlapping with promoter
              int glob_rna_idx = indiv->metadata_->rna_count_++;
              indiv->metadata_->rna_add(glob_rna_idx,
                                         prom_pos,
                                         rna_end,
                                         !lead_lag,
                                         BASAL_LEVELS[prom_error],
                                         rna_length, prom);
            }
            #endif
          }
        } else {
          /* Search for terminator */
          int cur_pos =
              prom_pos - 22;
          cur_pos =
              cur_pos < 0 ? dna_length + (cur_pos) : cur_pos;
          int start_pos = cur_pos;
          bool terminator_found = false;
          bool no_terminator = false;
          int term_dist_lagging = 0;

          int loop_size = 0;

          while (!terminator_found) {
            #ifdef BASE_2
            // ANNOTATE_SITE_BEGIN(lagging_motif_loop);
            if (cur_pos - 10 < length && cur_pos - 10 >= 0) {
              #pragma omp simd
              for (int32_t t_motif_id = 0; t_motif_id < 4; t_motif_id++) {
                int32_t pos_begin = cur_pos - t_motif_id;
                int32_t pos_end   = cur_pos + t_motif_id - 10;

                // ANNOTATE_ITERATION_TASK(searchLeading);
                term_dist_lagging +=
                    dna_data[pos_begin] !=
                    dna_data[pos_end]
                    ? 1 : 0;
              }
            } else {
              for (int32_t t_motif_id = 0; t_motif_id < 4; t_motif_id++) {
                int32_t pos_begin = cur_pos - t_motif_id;
                int32_t pos_end   = cur_pos + t_motif_id - 10;
                
                while (pos_begin < 0)       pos_begin += length;
                while (pos_begin >= length) pos_begin -= length;

                while (pos_end < 0)       pos_end += length;
                while (pos_end >= length) pos_end -= length;

                // ANNOTATE_ITERATION_TASK(searchLeading);
                term_dist_lagging +=
                    dna_data[pos_begin] !=
                    dna_data[pos_end]
                    ? 1 : 0;
              }
              // for (int32_t t_motif_id = 0; t_motif_id < 4; t_motif_id++) {
              //   term_dist_lagging +=
              //       dna->get_lag(cur_pos - t_motif_id) !=
              //       dna->get_lag(cur_pos + t_motif_id - 10)
              //       ? 1 : 0;
              // }
            }


            // ANNOTATE_SITE_END(lagging_motif_loop);

            if (term_dist_lagging == 4) {
              terminator_found = true;
            }
            else {
              cur_pos =
                  cur_pos - 1 < 0 ? dna_length + (cur_pos - 1)
                                  : cur_pos - 1;
              term_dist_lagging = 0;
              if (cur_pos == start_pos) {
                no_terminator = true;
                terminator_found = true;
              }
            }
            #elif BASE_4
            terminator_found = indiv->metadata_->is_terminator_lagging(cur_pos);

            if(!terminator_found) {
              cur_pos = cur_pos - 1 < 0 ?
                        dna_length + (cur_pos - 1) : cur_pos - 1;

              if (cur_pos == start_pos) {
                no_terminator = true;
                terminator_found = true;
              }
            }
            #endif
            loop_size++;
          }

          if (!no_terminator) {
            #ifdef BASE_2

            int32_t rna_end = Utils::mod(cur_pos - 10,dna_length);

            int32_t rna_length = 0;
            rna_length = loop_size + TERM_SIZE -1;

            if (rna_length >= 0) {
              int glob_rna_idx = -1;
              glob_rna_idx =
                  indiv->metadata_->rna_count_++;

              indiv->metadata_->rna_add(glob_rna_idx,
                                                                 prom_pos,
                                                                 rna_end,
                                                                 !lead_lag,
                                                                 1.0 -
                                                                 prom_error /
                                                                 5.0, rna_length,prom);
            }

            #elif BASE_4
            int32_t rna_end = Utils::mod(cur_pos - term_size, dna_length);
            int32_t rna_length = loop_size-1 + term_size;

            if (rna_length >= term_size) {
              int glob_rna_idx = -1;
              glob_rna_idx =
                  indiv->metadata_->rna_count_++;

              indiv->metadata_->rna_add(glob_rna_idx,
                                         prom_pos,
                                         rna_end,
                                         !lead_lag,
                                         BASAL_LEVELS[prom_error],
                                         rna_length,prom);
            }
            #endif
          }
        }
        } 
      }
    }
    // ANNOTATE_SITE_END(mainloop);

    #ifdef WITH_PERF_TRACES_PER_INDIV
      auto runtime_end = std::chrono::steady_clock::now();
      compute_rna_stop_[indiv->indiv_id] = runtime_end.time_since_epoch().count();
      compute_rna_[indiv->indiv_id] = runtime_end.time_since_epoch().count() - runtime_start.time_since_epoch().count();
    #endif
  
}

void ExpManager_7::compute_RNA(int indiv_id) {
    compute_RNA(current_individuals[indiv_id]);
}

void ExpManager_7::compute_RNA(Individual_7* indiv) {
  {
    indiv->metadata_->rnas_resize(
        indiv->metadata_->promoter_count());
    int32_t dna_length = indiv->dna_->length();

    #ifdef BASE_2
    int8_t term_size = 10;
    #elif BASE_4
    int8_t term_size = TERM_SIZE + exp_m_->exp_s()->terminator_polya_sequence_length();
    #endif

#ifdef WITH_FINETASKLOOP
#pragma omp taskloop grainsize(rna_grain_size)
#endif
    for (int rna_idx = 0; rna_idx <
                          (int)indiv->metadata_->promoter_count(); rna_idx++) {
      {
        if (indiv->metadata_->promoters(rna_idx) != nullptr) {
          if (indiv->metadata_->promoters(rna_idx)->leading_or_lagging) {
            if (indiv->metadata_->terminator_count(LEADING) != 0) {


              int k = indiv->metadata_->promoters(rna_idx)->pos + PROM_SIZE;

              k = k >= dna_length ? k - dna_length : k;

              int32_t next_rna_end =
                  indiv->metadata_->next_terminator(LEADING,k);

              int32_t rna_end =
                  next_rna_end + term_size >= dna_length ?
                  next_rna_end + term_size - dna_length :
                  next_rna_end + term_size;

              int32_t rna_length = 0;

#ifdef BASE_2
              if (k > next_rna_end)
                rna_length = dna_length + next_rna_end - k;
              else
                rna_length = next_rna_end - k;

              rna_length += TERM_SIZE;

              if (rna_length >= 0) {
#elif BASE_4
              if (indiv->metadata_->promoters(rna_idx)->pos > rna_end) {
                rna_length = dna_length -
                             indiv->metadata_->promoters(rna_idx)->pos +
                             rna_end + 1;
              } else {
                rna_length =
                    rna_end - indiv->metadata_->promoters(rna_idx)->pos + 1;
              }
              rna_length -= PROM_SIZE;

              if (rna_length >= term_size) { // check that terminator is not overlapping with promoter
#endif

                int glob_rna_idx = -1;
// #pragma omp critical
//                 {
                  glob_rna_idx =
                      indiv->metadata_->rna_count();
                  indiv->metadata_->set_rna_count(
                      indiv->metadata_->rna_count() + 1);
                

                indiv->metadata_->rna_add(glob_rna_idx, new Rna_7(indiv->metadata_->promoters(rna_idx)->pos,
                    rna_end,
                    !indiv->metadata_->promoters(rna_idx)->leading_or_lagging,
#ifdef BASE_2
                    1.0 -
                    fabs(
                        ((float)indiv->metadata_->promoters(rna_idx)->error)) /
                    5.0, rna_length,indiv->metadata_->promoters(rna_idx)));
#elif BASE_4
                    BASAL_LEVELS[abs(indiv->metadata_->promoters(rna_idx)->error)],
                    rna_length,indiv->metadata_->promoters(rna_idx)));
#endif
              }
            }
          } else {
            // LAGGING
            if (indiv->metadata_->terminator_count(LAGGING) != 0) {




              // Search for terminator
              int k = indiv->metadata_->promoters(rna_idx)->pos - PROM_SIZE;
              k = k < 0 ? dna_length + k : k;

              int32_t next_rna_end =
                  indiv->metadata_->next_terminator(LAGGING,k);


              int32_t rna_end =
                  next_rna_end - term_size < 0 ? dna_length + (next_rna_end - term_size)
                                        :
                  next_rna_end - term_size;

              int32_t rna_length = 0;
              
              #ifdef BASE_2
              if (k < next_rna_end)
                rna_length = k + dna_length - next_rna_end;
              else
                rna_length = k - next_rna_end;

              rna_length += TERM_SIZE;


              if (rna_length >= 0) {
              #elif BASE_4
              if (indiv->metadata_->promoters(rna_idx)->pos < rna_end) {
                rna_length = indiv->metadata_->promoters(rna_idx)->pos +
                             dna_length - rna_end + 1;
              }
              else {
                rna_length =
                    indiv->metadata_->promoters(rna_idx)->pos - rna_end + 1;
              }

              rna_length -= PROM_SIZE;

              if (rna_length >= term_size) { // check that terminator is not overlapping with promoter
              #endif

                int glob_rna_idx;
// #pragma omp critical
//                 {
                  glob_rna_idx =
                      indiv->metadata_->rna_count();
                  indiv->metadata_->set_rna_count(
                      indiv->metadata_->rna_count() + 1);
                

                indiv->metadata_->rna_add(glob_rna_idx, new Rna_7(indiv->metadata_->promoters(rna_idx)->pos,
                    rna_end,
                    !indiv->metadata_->promoters(rna_idx)->leading_or_lagging,
                    #ifdef BASE_2
                    1.0 -
                    fabs(
                        ((float)indiv->metadata_->promoters(rna_idx)->error)) /
                    5.0, rna_length,indiv->metadata_->promoters(rna_idx)));
                    #elif BASE_4
                    BASAL_LEVELS[abs(indiv->metadata_->promoters(rna_idx)->error)],
                    rna_length,indiv->metadata_->promoters(rna_idx)));
                    #endif

              }
            }
          }
        }
      }
    }
  }
}

void ExpManager_7::start_protein(int indiv_id) {
  start_protein(current_individuals[indiv_id]);
}

void ExpManager_7::start_protein(Individual_7* indiv) {
        #ifdef WITH_PERF_TRACES_PER_INDIV
      auto runtime_start = std::chrono::steady_clock::now();
      start_protein_start_[indiv->indiv_id] = runtime_start.time_since_epoch().count();
    #endif

  indiv->metadata_->rna_begin();
  #ifdef BASE_4
  int8_t term_size = TERM_SIZE + exp_m_->exp_s()->terminator_polya_sequence_length();
  #endif

  for (int rna_idx = 0; rna_idx <
                        (int)indiv->metadata_->rna_count(); rna_idx++) {
      Rna_7* rna =
          indiv->metadata_->rna_next();

      int32_t dna_length = indiv->dna_->length_;
      
      #ifdef BASE_2
      #ifdef WITH_OPTIMIZE_DIFF_SEARCH
      if (true) {
      #else
      if (rna->to_compute_start_pos_) { //
      #endif
      //if (true) {
          
          rna->start_prot.clear();
          rna->start_prot_count_ = 0;
          
      char* dna_data = indiv->dna_->data_;
      #if defined(__INTEL_COMPILER)
        __declspec(align(dna_data));
      #elif defined(__INTEL_LLVM_COMPILER)
          void* vec_r = __builtin_assume_aligned(dna_data,64);
      #endif

      int32_t rna_length = rna->length;
      bool rna_leading_lagging = rna->leading_lagging;


          // if (indiv_id==0)
          //   printf("[1/3] Looking for start in %d %d %d\n",rna->begin,rna->is_init_,rna_length);
      if (rna->is_init_) {
        int32_t s_pos = rna->begin;

          // if (indiv_id==373)
          //   printf("[2/3] Looking for start in %d %d %d\n",rna->begin,rna->is_init_,rna_length);

        if (rna_length >= 21) {
          //bool a_start[rna_length];
          //memset(a_start, 0, rna_length);
          
          // if (indiv_id==373)
          //   printf("[3/3] Looking for start in %d %d %d\n",rna->begin,rna->is_init_,rna_length);
          if (rna_leading_lagging == 0) {
            s_pos += 22;
            s_pos =
                s_pos >= dna_length ? s_pos - dna_length
                                    : s_pos;
          } else {
            s_pos -= 22;
            s_pos =
                s_pos < 0 ? ((int) dna_length) + s_pos : s_pos;
          }


          int32_t loop_size = 0;
          //#pragma omp simd

          #ifdef __VECTORIZE_STRCMP
          char buffer_str[13];
          #endif
    
          for (int32_t loop_size = 0; loop_size < rna_length - DO_TRANSLATION_LOOP; loop_size++) {


          // int32_t c_pos = s_pos;
          #ifdef __VECTORIZE_STRCMP
            bool b_start = false;
            bool start[14] = {false,false,false,false,false,false,false,false,false,false,false,false,false,false};
          #else
            #if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
              __declspec(align(64)) bool start[14] = {false,false,false,false,false,false,false,false,false,false,false,false,false,false};
            #else
              bool start[14] = {false,false,false,false,false,false,false,false,false,false,false,false,false,false};
            #endif
          #endif
          // while (loop_size+DO_TRANSLATION_LOOP < rna->length) {
            // bool start = false;
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
                #ifdef __VECTORIZE_STRCMP
                bool compare_1 = strncmp(&(dna_data[c_pos]),SHINE_DAL_SEQ_LEAD_7,6) == 0 ? true : false;
                bool compare_2 = strncmp(&(dna_data[c_pos+10]),&(SHINE_DAL_SEQ_LEAD_7[10]),3) == 0 ? true : false;
                b_start = compare_1 && compare_2;
                #else
                // char c_sub[13] = {};
                #if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
                #pragma omp simd aligned(dna_data:64) aligned(SHINE_DAL_SEQ_LEAD_7:64) aligned(start:64)
                #else
                #pragma omp simd
                #endif
                for (int32_t k = 0; k < 12; k++) {
                  // int32_t k_t = k >= 6 ? k + 4 : k;
                 start[k] = dna_data[c_pos+k] == SHINE_DAL_SEQ_LEAD_7[k];
                //  c_sub[k] = dna_data[c_pos+k];
                  // start[0] = (dna_data[c_pos    ] == SHINE_DAL_SEQ_LEAD[0]) ? true: false;
                  // start[1] = (dna_data[c_pos + 1] == SHINE_DAL_SEQ_LEAD[1]) ? true: false;
                  // start[2] = (dna_data[c_pos + 2] == SHINE_DAL_SEQ_LEAD[2]) ? true: false;
                  // start[3] = (dna_data[c_pos + 3] == SHINE_DAL_SEQ_LEAD[3]) ? true: false;
                  // start[4] = (dna_data[c_pos + 4] == SHINE_DAL_SEQ_LEAD[4]) ? true: false;
                  // start[5] = (dna_data[c_pos + 5] == SHINE_DAL_SEQ_LEAD[5]) ? true: false;
                  // start[6] = (dna_data[c_pos + 10] == SHINE_DAL_SEQ_LEAD[6]) ? true: false;
                  // start[7] = (dna_data[c_pos + 11] == SHINE_DAL_SEQ_LEAD[7]) ? true: false;
                  // start[8] = (dna_data[c_pos + 12] == SHINE_DAL_SEQ_LEAD[8]) ? true: false;
                }
                start[12] = dna_data[c_pos+12] == SHINE_DAL_SEQ_LEAD_7[12];
                // c_sub[12] = dna_data[c_pos+12];

                // #ifdef __VECTORIZE_STRCMP
                // b_start = start[0] && start[1] && start[2] && start[3] && start[4] && start[5] && start[10] && start[11] && start[12];
                // bool compare_1 = strncmp(&(dna_data[c_pos]),SHINE_DAL_SEQ_LEAD_7,6) == 0 ? true : false;
                // bool compare_2 = strncmp(&(dna_data[c_pos+10]),&(SHINE_DAL_SEQ_LEAD_7[10]),3) == 0 ? true : false;
                // if ((compare_1 && compare_2) != b_start) {
                //   for (int32_t k = 0; k < 13; k++) {
                //     if (k < 6 or k > 9) printf("%c",dna_data[c_pos+k]);
                //     else if (k==6) printf(" ");
                //   }
                //   printf(" -- ");
                //   for (int32_t k = 0; k < 13; k++) {
                //     printf("%c",SHINE_DAL_SEQ_LEAD_7[k]);
                //     if (k==5) { printf(" "); k+=4; }
                //   }
                //   printf(" -- Incorrect CMP %.6s != %.6s (%.6s)",
                //     &(dna_data[c_pos]),SHINE_DAL_SEQ_LEAD_7,c_sub);
                //   printf(" AND %.3s != %.3s : %d != (%d && %d) -- %d %d %d %d %d %d %d %d %d\n",
                //     &(dna_data[c_pos+10]),&(SHINE_DAL_SEQ_LEAD_7[6]),
                //     b_start,compare_1,compare_2,start[0],start[1],start[2],start[3],start[4],start[5],start[10],start[11],start[12]);
                // }
                // #endif
                #endif
              } else {
                
                for (int32_t k = 0; k < 9; k++) {
                  int32_t k_t = k >= 6 ? k + 4 : k;
                  int32_t pos_m = c_pos + k_t;

                  while (pos_m < 0)  pos_m += dna_length;
                  while (pos_m >= dna_length) pos_m -= dna_length;

                  // #ifdef __VECTORIZE_STRCMP
                  // start = (dna_data[pos_m] == SHINE_DAL_SEQ_LEAD_7[k_t]) ? true: false;
                  // if (!start)
                  //   break;
                  // #else
                  start[k_t] = (dna_data[pos_m] == SHINE_DAL_SEQ_LEAD_7[k_t]) ? true: false;
                  // #endif
                }
                #ifdef __VECTORIZE_STRCMP
                b_start = start[0] && start[1] && start[2] && start[3] && start[4] && start[5] && start[10] && start[11] && start[12];
                #endif
              }
              // for (int k = 0; k < 9; k++) {
              //   int32_t k_t = k >= 6 ? k + 4 : k;

              //   if (current_individuals[indiv_id]->dna_->data_[Utils::mod((c_pos + k_t),dna_length)] ==
              //       SHINE_DAL_SEQ_LEAD[k]) {
              //     start = true;
              //   } else {
              //     start = false;
              //     break;
              //   }

              // }

            } else {
              // Search for Shine Dalgarro + START codon on LAGGING
              if (c_pos - 15 < dna_length && c_pos - 15 >= 0) {
                #ifdef __VECTORIZE_STRCMP
                for (int i = 0; i < 13; i++) buffer_str[i] = dna_data[c_pos-i];
                bool compare_1 = strncmp(buffer_str,SHINE_DAL_SEQ_LAG_7,6) == 0 ? true : false;
                bool compare_2 = strncmp(&(buffer_str[10]),&(SHINE_DAL_SEQ_LAG_7[10]),3) == 0 ? true : false;

                b_start = compare_1 && compare_2;
                #else
                #pragma omp simd
                for (int32_t k = 0; k < 12; k++) {
                  // int32_t k_t = k >= 6 ? k + 4 : k;
                  start[k] = dna_data[c_pos - k] == SHINE_DAL_SEQ_LAG_7[k];
                  // start[0] = (dna_data[c_pos    ] == SHINE_DAL_SEQ_LAG[0]) ? true: false;
                  // start[1] = (dna_data[c_pos - 1] == SHINE_DAL_SEQ_LAG[1]) ? true: false;
                  // start[2] = (dna_data[c_pos - 2] == SHINE_DAL_SEQ_LAG[2]) ? true: false;
                  // start[3] = (dna_data[c_pos - 3] == SHINE_DAL_SEQ_LAG[3]) ? true: false;
                  // start[4] = (dna_data[c_pos - 4] == SHINE_DAL_SEQ_LAG[4]) ? true: false;
                  // start[5] = (dna_data[c_pos - 5] == SHINE_DAL_SEQ_LAG[5]) ? true: false;
                  // start[6] = (dna_data[c_pos - 10] == SHINE_DAL_SEQ_LAG[6]) ? true: false;
                  // start[7] = (dna_data[c_pos - 11] == SHINE_DAL_SEQ_LAG[7]) ? true: false;
                  // start[8] = (dna_data[c_pos - 12] == SHINE_DAL_SEQ_LAG[8]) ? true: false;
                }
                start[12] = dna_data[c_pos - 12] == SHINE_DAL_SEQ_LAG_7[12];
                // #ifdef __VECTORIZE_STRCMP
                // b_start = start[0] && start[1] && start[2] && start[3] && start[4] && start[5] && start[10] && start[11] && start[12];
                // #endif
                #endif    
              } else {
                for (int32_t k = 0; k < 9; k++) {
                  int32_t k_t = k >= 6 ? k + 4 : k;
                  int32_t pos_m = c_pos - k_t;

                  while (pos_m < 0)  pos_m += dna_length;
                  while (pos_m >= dna_length) pos_m -= dna_length;

                  // #ifdef __VECTORIZE_STRCMP
                  // start = (dna_data[pos_m] == SHINE_DAL_SEQ_LAG_7[k_t]) ? true: false;
                  // if (!start)
                  //   break;
                  // #else
                  start[k_t] = (dna_data[pos_m] == SHINE_DAL_SEQ_LAG_7[k_t]) ? true : false;
                  // #endif
                }
                #ifdef __VECTORIZE_STRCMP
                b_start = start[0] && start[1] && start[2] && start[3] && start[4] && start[5] && start[10] && start[11] && start[12];
                #endif
              }
              // for (int k = 0; k < 9; k++) {
              //   int32_t k_t = k >= 6 ? k + 4 : k;

              //   if (current_individuals[indiv_id]->dna_->data_[Utils::mod((c_pos - k_t),dna_length)] ==
              //       SHINE_DAL_SEQ_LAG[k]) {
              //     start = true;
              //   } else {
              //     start = false;
              //     break;
              //   }
              // }

            }

            #ifdef __VECTORIZE_STRCMP
            if (b_start) {
            #else
            if (start[0] && start[1] && start[2] && start[3] && start[4] && start[5] && start[10] && start[11] && start[12]) {
            #endif
              rna->start_prot_count_++;
                          
              // if (indiv_id==373)
                // if (indiv->indiv_id == 0) printf("[5/3] Adding start %d\n",c_pos);

              // a_start[loop_size] = start[0] && start[1] && start[2] && start[3] && start[4] && start[5] && start[6] && start[7] && start[8];
              rna->start_prot.push_back(c_pos);

          // printf("RNA %d : Protein %d\n",prom_idx,c_pos);
            }
            // if (rna->leading_lagging == 0) {
            //   c_pos++;
            //   c_pos =
            //       c_pos >= dna_length ? c_pos - dna_length
            //                           : c_pos;
            // } else {
            //   c_pos--;
            //   c_pos = c_pos < 0 ? dna_length + c_pos : c_pos;
            // }
            // loop_size++;


            //loop_size++;
          }

          // for (int i = 0; i < rna_length; i++) {
          //   if (a_start[i]) {
          //     rna->start_prot_count_++;
            
          //     int32_t c_pos = s_pos;
          //     if (rna_leading_lagging == 0) {
          //       c_pos+=loop_size;
          //       c_pos =
          //           c_pos >= dna_length ? c_pos - dna_length
          //                               : c_pos;
          //     } else {
          //       c_pos-=loop_size;
          //       c_pos = c_pos < 0 ? dna_length + c_pos : c_pos;
          //     }

          //     rna->start_prot.push_back(c_pos);
          //   }
          // }
        }
      }

      // if (indiv_id==373)
      //   for (auto posi : rna->start_prot)
      //     printf("RNA %d : Protein %d\n",prom_idx,posi);

      
      }
  }
  #elif BASE_4

      if (rna->is_init_ && rna->length >= DO_TRANSLATION_LOOP) {
        // if (indiv->indiv_id == 6) printf("%d -- %d -- Process RNA %d => %d\n",AeTime::time(),indiv->indiv_id,rna->begin,rna->end);

        int c_pos = rna->begin;

        if(rna->leading_lagging == LEADING) {
          c_pos += PROM_SIZE;
          c_pos =
              c_pos >= dna_length ? c_pos - dna_length : c_pos;
        } else {
          c_pos -= PROM_SIZE;
          c_pos =
              c_pos < 0 ? dna_length + c_pos : c_pos;
        }

        int loop_size = 0;
        while (loop_size+DO_TRANSLATION_LOOP < rna->length) {
          bool start = false;

          if (rna->leading_lagging == LEADING) {
            // Search for Shine Dalgarro + START codon on LEADING
            start = indiv->metadata_->is_shine_dal_start_prot_leading(c_pos);
          } else {
            // Search for Shine Dalgarro + START codon on LAGGING
            start = indiv->metadata_->is_shine_dal_start_prot_lagging(c_pos);
          }

            // if (indiv->indiv_id == 229 && c_pos == 238) printf("%d -- %d -- Found a start at %d for RNA %d (%s) : %c %c %c %c %c %c\n",AeTime::time(),
            //         indiv->indiv_id,c_pos,rna->begin,
            //         rna->leading_lagging == LEADING ? "LEADING" : "LAGGING",
            //         indiv->dna_->data_[c_pos],
            //         indiv->dna_->data_[c_pos - 1],
            //         indiv->dna_->data_[c_pos - 2],
            //         indiv->dna_->data_[c_pos - 3],
            //         indiv->dna_->data_[c_pos - 4],
            //         indiv->dna_->data_[c_pos - 5],
            //         indiv->dna_->data_[c_pos - SHINE_DAL_SIZE - SHINE_START_SPACER],
            //         indiv->dna_->data_[c_pos - SHINE_DAL_SIZE - SHINE_START_SPACER - 1],
            //         indiv->dna_->data_[c_pos - SHINE_DAL_SIZE - SHINE_START_SPACER - 2]);

          if (start) {
            rna->start_prot_count_++;
            rna->start_prot.push_back(c_pos);

            // if (indiv->indiv_id == 229) printf("%d -- %d -- Found a start at %d for RNA %d (%s) : %c %c %c %c %c %c\n",AeTime::time(),
            //         indiv->indiv_id,c_pos,rna->begin,
            //         rna->leading_lagging == LEADING ? "LEADING" : "LAGGING",
            //         indiv->dna_->data_[c_pos],
            //         indiv->dna_->data_[c_pos - 1],
            //         indiv->dna_->data_[c_pos - 2],
            //         indiv->dna_->data_[c_pos - 3],
            //         indiv->dna_->data_[c_pos - 4],
            //         indiv->dna_->data_[c_pos - 5],
            //         indiv->dna_->data_[c_pos - SHINE_DAL_SIZE - SHINE_START_SPACER],
            //         indiv->dna_->data_[c_pos - SHINE_DAL_SIZE - SHINE_START_SPACER - 1],
            //         indiv->dna_->data_[c_pos - SHINE_DAL_SIZE - SHINE_START_SPACER - 2]);
          }

          if (rna->leading_lagging == LEADING) {
            c_pos++;
            c_pos =
                c_pos >= dna_length ? c_pos - dna_length
                                    : c_pos;
          } else {
            c_pos--;
            c_pos = c_pos < 0 ? dna_length + c_pos : c_pos;
          }
          loop_size++;
        }
      }
    }
    #endif

      #ifdef WITH_PERF_TRACES_PER_INDIV
      auto runtime_end = std::chrono::steady_clock::now();
      start_protein_stop_[indiv->indiv_id] = runtime_end.time_since_epoch().count();
      start_protein_[indiv->indiv_id] = runtime_end.time_since_epoch().count() - runtime_start.time_since_epoch().count();
    #endif
}

void ExpManager_7::compute_protein(int indiv_id) {
  compute_protein(current_individuals[indiv_id]);
}

void ExpManager_7::compute_protein(Individual_7* indiv) {
      #ifdef WITH_PERF_TRACES_PER_INDIV
      auto runtime_start = std::chrono::steady_clock::now();
      compute_protein_start_[indiv->indiv_id] = runtime_start.time_since_epoch().count();
    #endif

  int resize_to = 0;

    
    indiv->metadata_->rna_begin();
    
    // ANNOTATE_SITE_BEGIN(mainloop);
    for (int prom_idx = 0;
         prom_idx < (int)indiv->metadata_->rna_count(); prom_idx++) {
          //  ANNOTATE_ITERATION_TASK(searchMain);
      Rna_7* rna =
         indiv->metadata_->rna_next();

    if (rna->is_init_)
      resize_to += rna->start_prot_count_;
  }

  indiv->
      metadata_->proteins_resize(resize_to);

  Dna_7* dna = indiv->dna_;
  int32_t dna_length = dna->length();

  #ifdef BASE_4
  int8_t term_size = TERM_SIZE + exp_m_->exp_s()->terminator_polya_sequence_length();
  #endif

  indiv->metadata_->rna_begin();
  for (int rna_idx = 0; rna_idx <
                        (int)indiv->metadata_->rna_count(); rna_idx++) {

    Rna_7* rna = indiv->metadata_->rna_next();

    if (rna->is_init_) {
      for (auto it_start_pos = rna->start_prot.begin(); it_start_pos != rna->start_prot.end(); it_start_pos++) {

        int start_prot = *it_start_pos;
        int start_protein_pos = rna->leading_lagging == 0 ?
                                start_prot +
                                PROT_START_SIZE :
                                start_prot -
                                PROT_START_SIZE;
        int length;

        if (rna->leading_lagging == 0) {
          start_protein_pos = start_protein_pos >= dna_length ?
                              start_protein_pos - dna_length
                                                              : start_protein_pos;

          if (start_prot < rna->end) {
            length = rna->end - start_prot;
          } else {
            length = dna_length -
                     start_prot +
                     rna->end;

          }

          length -= PROT_START_SIZE;
        } else {


          start_protein_pos = start_protein_pos < 0 ?
                              dna_length + start_protein_pos
                                                    : start_protein_pos;

          if (start_prot > rna->end) {
            length = (*it_start_pos) - rna->end;
          } else {
            length = *it_start_pos +
                     dna_length -
                     rna->end;
          }

          length -= PROT_START_SIZE;
        }

        bool is_protein = false;

        length += 1;
        length = length - (length % CODON_SIZE);

        int j = 0;
        int32_t transcribed_start = 0;

        if (rna->leading_lagging == 0) {
          transcribed_start = rna->begin + PROM_SIZE;
          transcribed_start = transcribed_start >= dna_length ?
                              transcribed_start - dna_length
                                                              : transcribed_start;

          if (transcribed_start <= start_prot) {
            j = start_prot - transcribed_start;
          } else {
            j = dna_length -
                transcribed_start +
                start_prot;

          }
        } else {
          transcribed_start = rna->begin - PROM_SIZE;
          transcribed_start = transcribed_start < 0 ?
                              dna_length + transcribed_start
                                                    : transcribed_start;

          if (transcribed_start >=
              start_prot) {
            j = transcribed_start -
                start_prot;
          } else {
            j = transcribed_start +
                dna_length - start_prot;
          }
        }

        j += PROT_START_SIZE;

        while (rna->length - j >= CODON_SIZE) {
          int t_k;

          if (rna->leading_lagging == 0) {
            start_protein_pos =
                start_protein_pos >= dna_length ?
                start_protein_pos - dna_length
                                                : start_protein_pos;

            #ifdef BASE_2
            is_protein = false;

            for (int k = 0; k < 3; k++) {
              t_k = start_protein_pos + k;

              if (dna->get_lead(t_k) ==
                  PROTEIN_END_LEAD[k]) {
                is_protein = true;
              } else {
                is_protein = false;
                break;
              }
            }
            #elif BASE_4
            is_protein = codon_value_to_aminoacid(bases_to_codon_value(
                  dna->get_lead(start_protein_pos),
                  dna->get_lead(start_protein_pos + 1),
                  dna->get_lead(start_protein_pos + 2)
                )) == STOP;

            // if (indiv->indiv_id == 884)
            //   printf("Search for STOP %d :: %c %c %c :: %d\n",start_protein_pos,
            //       dna->get_lead(start_protein_pos),
            //       dna->get_lead(start_protein_pos + 1),
            //       dna->get_lead(start_protein_pos + 2),is_protein);

                  t_k = start_protein_pos + CODON_SIZE;
            #endif

            if (is_protein) {
              int prot_length = -1;

              if (start_prot + PROT_START_SIZE < t_k) {
                prot_length = t_k -
                              (start_prot +
                               PROT_START_SIZE);
              } else {
                prot_length = dna_length -
                              (start_prot +
                               PROT_START_SIZE) + t_k;
              }

              // printf("Protein Length %d\n",prot_length);

              if (prot_length >= CODON_SIZE) {
                int32_t glob_prot_idx = -1;
                glob_prot_idx =
                    indiv->metadata_->proteins_count();
                indiv->metadata_->set_proteins_count(
                    indiv->metadata_->proteins_count() +
                    1);

                #ifdef BASE_4
                int32_t prot_len = prot_length/CODON_SIZE;
                prot_len--;
                if (prot_len > 0) {
                #endif

                // printf("%d -- COMPUTE PROT -- Indiv %d Add prot %d\n",AeTime::time(),indiv->indiv_id,glob_prot_idx);

                indiv->
                    metadata_->protein_add(glob_prot_idx, new Protein_7(
                    Utils::mod(start_prot+PROT_START_SIZE,dna_length), Utils::mod(t_k,dna_length),
                        #ifdef BASE_4
                        prot_len,
                        #else
                        prot_length/3,
                        #endif                    
                    rna->leading_lagging,
                    rna->e,rna,true
                ));

                rna->is_coding_ = true;
                #ifdef BASE_4
                }
                #endif
              }

              break;
            }

            start_protein_pos += CODON_SIZE;
            start_protein_pos =
                start_protein_pos >= dna_length ?
                start_protein_pos - dna_length
                                                : start_protein_pos;
          } else {


            start_protein_pos = start_protein_pos < 0 ?
                                dna_length + start_protein_pos
                                                      : start_protein_pos;

            #ifdef BASE_2
            is_protein = false;
            for (int k = 0; k < 3; k++) {
              t_k = start_protein_pos - k;

              if (dna->get_lag(t_k) ==
                  PROTEIN_END_LAG[k]) {
                is_protein = true;
              } else {
                is_protein = false;
                break;
              }
            }
            #elif BASE_4
            is_protein = codon_value_to_aminoacid(bases_to_codon_value(
                dna->get_lag(start_protein_pos),
                dna->get_lag(start_protein_pos - 1),
                dna->get_lag(start_protein_pos - 2)
            )) == STOP;


            // if (indiv->indiv_id == 884 && start_protein_pos == 2408)
            //   printf("Search for STOP :: %c %c %c\n",
            //       dna->get_lead(start_protein_pos),
            //       dna->get_lead(start_protein_pos + 1),
            //       dna->get_lead(start_protein_pos + 2));

                   t_k = start_protein_pos - CODON_SIZE;
            #endif

            if (is_protein) {
              int prot_length = -1;
              if (start_prot - PROT_START_SIZE > t_k) {
                prot_length =
                    (start_prot - PROT_START_SIZE) -
                    t_k;
              } else {
                prot_length =
                    (start_prot - PROT_START_SIZE) +
                    dna_length - t_k;
              }
              if (prot_length >= CODON_SIZE) {
                int32_t glob_prot_idx = -1;

                glob_prot_idx =
                    indiv->metadata_->proteins_count();
                indiv->metadata_->set_proteins_count(
                    indiv->metadata_->proteins_count() +
                    1);


                #ifdef BASE_4
                int32_t prot_len = prot_length/CODON_SIZE;
                prot_len--;

                if (prot_len > 0) {
                #endif

                // printf("%d -- COMPUTE PROT -- Indiv %d Add prot %d\n",AeTime::time(),indiv->indiv_id,glob_prot_idx);

                ((List_Metadata*)indiv->metadata_)->protein_add(
                    glob_prot_idx,
                        new Protein_7(Utils::mod(start_prot-PROT_START_SIZE,dna_length), Utils::mod(t_k,dna_length),
                        #ifdef BASE_4
                        prot_len,
                        #else
                        prot_length/3,
                        #endif
                        rna->leading_lagging,
                        rna->e,rna,true
                    ));
                rna->is_coding_ = true;
                #ifdef BASE_4
                }
                #endif
              }
              break;
            }
            start_protein_pos = start_protein_pos - CODON_SIZE;
            start_protein_pos = start_protein_pos < 0 ?
                                dna_length + start_protein_pos
                                                      : start_protein_pos;
          }
          j += CODON_SIZE;
        }
      }
    }
  }


    #ifdef WITH_PERF_TRACES_PER_INDIV
      auto runtime_end = std::chrono::steady_clock::now();
      compute_protein_stop_[indiv->indiv_id] = runtime_end.time_since_epoch().count();
      compute_protein_[indiv->indiv_id] = runtime_end.time_since_epoch().count() - runtime_start.time_since_epoch().count();
    #endif
}

void ExpManager_7::translate_protein(int indiv_id, double w_max) {
  translate_protein(current_individuals[indiv_id],w_max);
}

void ExpManager_7::translate_protein(Individual_7* indiv, double w_max) {
      #ifdef WITH_PERF_TRACES_PER_INDIV
      auto runtime_start = std::chrono::steady_clock::now();
      translate_protein_start_[indiv->indiv_id] = runtime_start.time_since_epoch().count();
    #endif

  int32_t dna_length = indiv->dna_->length();

  indiv->metadata_->protein_begin();
  for (int protein_idx = 0; protein_idx <
                            (int)indiv->metadata_->proteins_count(); protein_idx++) {
    Protein_7* prot = indiv->metadata_->protein_next();

    if (prot->is_init_) {
      int c_pos = prot->protein_start, t_pos;
      int end_pos = prot->protein_end;
      if (prot->leading_lagging ==
          0) {
        end_pos -= 3;

        c_pos =
            c_pos >= dna_length ? c_pos - dna_length
                                        : c_pos;
        end_pos = end_pos < 0 ? dna_length + end_pos : end_pos;
      } else {
        end_pos += 3;

        end_pos =
            end_pos >= dna_length ? end_pos - dna_length
                                          : end_pos;
        c_pos = c_pos < 0 ? dna_length + c_pos : c_pos;
      }

      int8_t value = 0;
      int16_t codon_idx = 0;
      int32_t count_loop = 0;

      #ifdef BASE_4
      char base_1 = -1;
      char base_2 = -1;
      char base_3 = -1;
      #endif

      if (prot->leading_lagging ==
          0) {
        // LEADING
        #ifdef BASE_2
        while (count_loop <
               prot->protein_length &&
               codon_idx < 64*3) {
          value = 0;
          for (int8_t i = 0; i < 3; i++) {
            t_pos = c_pos + i;
            if (indiv->dna_->get_lead(t_pos) ==
                '1')
              value += 1 << (CODON_SIZE - i - 1);
          }

          prot->codon_list[codon_idx] = value;
        #elif BASE_4
        while (count_loop < prot->protein_length) {
          base_1 = indiv->dna_->get_lead(c_pos);
          base_2 = indiv->dna_->get_lead(c_pos+1);
          base_3 = indiv->dna_->get_lead(c_pos+2);

          value = bases_to_codon_value(base_1, base_2, base_3);
          prot->codon_list.push_back(value);
        #endif

          codon_idx++;

          count_loop++;
          c_pos += CODON_SIZE;
          c_pos =
              c_pos >= dna_length ? c_pos - dna_length
                                          : c_pos;
        }
      } else {
        // LAGGING
        #ifdef BASE_2
        while (count_loop <
               prot->protein_length &&
               codon_idx < 64*3) {
          value = 0;
          for (int8_t i = 0; i < 3; i++) {
            t_pos = c_pos - i;
            if (indiv->dna_->get_lag(t_pos) !=
                '1')
              value += 1 << (CODON_SIZE - i - 1);
          }

          // if (indiv_id == 660) printf("Protein %d :: Add codon %d : %d\n",prot->protein_start,codon_idx,value);
          prot->codon_list[codon_idx] = value;

        #elif BASE_4
        while (count_loop < prot->protein_length) {

          // get complementary base on lagging strain
          base_1 = get_complementary_base(indiv->dna_->get_lead(c_pos));
          base_2 = get_complementary_base(indiv->dna_->get_lead(c_pos - 1));
          base_3 = get_complementary_base(indiv->dna_->get_lead(c_pos - 2));

          value = bases_to_codon_value(base_1, base_2, base_3);
          prot->codon_list.push_back(value);
        #endif

          codon_idx++;

          count_loop++;

          c_pos -= CODON_SIZE;
          c_pos = c_pos < 0 ? c_pos + dna_length : c_pos;
        }
      }

      #ifdef BASE_2
      if (codon_idx>=64*3) {
        std::ofstream last_gener_file("aevol_run.log",
                                    std::ofstream::out);
        last_gener_file << "Stop translating protein before end (length is greater than 196" << std::endl;
        std::cout << "Stop translating protein before end (length is greater than 196" << std::endl;
        last_gener_file.close();
      }

      double M = 0.0;
      double W = 0.0;
      double H = 0.0;
      #elif BASE_4
      //  --------------------------------
      //  1) Compute values for M, W and H
      //  --------------------------------
      auto base_m = exp_m_->exp_s()->aa_base_m();
      auto base_w = exp_m_->exp_s()->aa_base_w();
      auto base_h = exp_m_->exp_s()->aa_base_h();

      int8_t base_m_size = exp_m_->exp_s()->aa_base_m_size();
      int8_t base_w_size = exp_m_->exp_s()->aa_base_w_size();
      int8_t base_h_size = exp_m_->exp_s()->aa_base_h_size();

      long double M = 0.0;
      long double W = 0.0;
      long double H = 0.0;

      prot->nb_codons_ = codon_idx;
      #endif

      int32_t nb_m = 0;
      int32_t nb_w = 0;
      int32_t nb_h = 0;

      #ifdef BASE_2
      bool bin_m = false; // Initializing to false will yield a conservation of the high weight bit
      bool bin_w = false; // when applying the XOR operator for the Gray to standard conversion
      bool bin_h = false;

      prot->nb_codons_ = codon_idx-1;


      for (int i = 0; i < codon_idx; i++) {
        switch (prot->codon_list[i]) {
        case CODON_M0 : {
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
        case CODON_M1 : {
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
        case CODON_W0 : {
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
        case CODON_W1 : {
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
        case CODON_H1 : {
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

      
      double m_base_exponent=1.25;
      double w_base_exponent=1.25;
      double h_base_exponent=1.25;

      double m_digit_factor = 1.0;
      double w_digit_factor = 1.0;
      double h_digit_factor = 1.0;

      for(int i = 0; i < prot->nb_codons_; i++) {
        AminoAcid amino_acid = codon_value_to_aminoacid(prot->codon_list.at(prot->nb_codons_ - i - 1));

        if(base_m[amino_acid] != -1) {
          M += base_m[amino_acid] * m_digit_factor;
          m_digit_factor *= m_base_exponent;//base_m_size;
          nb_m++;
        }

        if(base_w[amino_acid] != -1) {
          W += base_w[amino_acid] * w_digit_factor;
          w_digit_factor *= w_base_exponent;//base_w_size;
          nb_w++;
        }

        if(base_h[amino_acid] != -1) {
          H += base_h[amino_acid] * h_digit_factor;
          h_digit_factor *= h_base_exponent;//base_h_size;
          nb_h++;
        }
      }
      #endif


      //  ----------------------------------------------------------------------------------
      //  2) Normalize M, W and H values in [0;1] according to number of codons of each kind
      //  ----------------------------------------------------------------------------------
      #ifdef BASE_2
      prot->m =
          nb_m != 0 ? M / (pow(2, nb_m) - 1) : 0.5;
      prot->w =
          nb_w != 0 ? W / (pow(2, nb_w) - 1) : 0.0;
      prot->h =
          nb_h != 0 ? H / (pow(2, nb_h) - 1) : 0.5;
      #elif BASE_4
      /*
      if(nb_m != 0) {
        prot->m = M / (pow(base_m_size, nb_m) - 1);
      } else {
        prot->m = 0.5;
      }

      if(nb_w != 0) {
        prot->w = W / (pow(base_w_size, nb_w) - 1);
      }
      else {
        prot->w = 0.0;
      }

      if(nb_h != 0) {
        prot->h = H / (pow(base_h_size, nb_h) - 1);
      } else {
        prot->h = 0.5;
	}*/

      
      if(nb_m != 0) {
        prot->m = M / ((base_m_size)*(pow(m_base_exponent, nb_m) - 1)/(m_base_exponent - 1));
      } else {
        prot->m = 0.5;
      }

      if(nb_w != 0) {
      	prot->w = W / ((base_w_size)*(pow(w_base_exponent, nb_w) - 1)/(w_base_exponent - 1));
      }
      else {
        prot->w = 0.0;
      }

      if(nb_h != 0) {
        prot->h = H / ((base_h_size)*(pow(h_base_exponent, nb_h) - 1)/(h_base_exponent - 1));
      } else {
        prot->h = 0.5;
      }
 
      #endif
      //  ------------------------------------------------------------------------------------
      //  3) Normalize M, W and H values according to the allowed ranges (defined in macros.h)
      //  ------------------------------------------------------------------------------------
      // x_min <= M <= x_max
      // w_min <= W <= w_max
      // h_min <= H <= h_max
      prot->m =
          (X_MAX - X_MIN) *
          prot->m +
          X_MIN;
      prot->w =
          (w_max - W_MIN) *
          prot->w +
          W_MIN;
      prot->h =
          (H_MAX - H_MIN) *
          prot->h +
          H_MIN;

      if (nb_m == 0 || nb_w == 0 || nb_h == 0 ||
          prot->w ==
          0.0 ||
          prot->h ==
          0.0) {
        prot->is_functional = false;
      } else {
        prot->is_functional = true;
      }
    }
  }

std::map<int32_t, Protein_7*> lookup;

  indiv->metadata_->protein_begin();

  for (int protein_idx = 0; protein_idx <
                            (int)indiv->metadata_->proteins_count(); protein_idx++) {
    {

      Protein_7* prot  =
          indiv->metadata_->protein_next();
      if (prot->is_init_ && prot->leading_lagging==0) {
        if (lookup.find(prot->protein_start) ==
            lookup.end()) {
          lookup[prot->protein_start] = prot;
        } else {
          #ifdef __MULTI_PROMOTERS
            #if 0 == __MULTI_PROMOTERS
              if(prot->e > lookup[prot->protein_start]->e)
                lookup[prot->protein_start]->e = prot->e;
              if(prot->initial_e_ > lookup[prot->protein_start]->initial_e_)
                lookup[prot->protein_start]->initial_e_ = prot->initial_e_;
            #endif
          #else
          lookup[prot->protein_start]->e += prot->e;
          lookup[prot->protein_start]->initial_e_ += prot->initial_e_;
          #endif
          lookup[prot->protein_start]->rna_list_.insert(
              lookup[prot->protein_start]->rna_list_.end(),
              prot->rna_list_.begin(),prot->rna_list_.end());
          prot->is_init_ = false;
        }
      }
    }
  }

  lookup.clear();

  indiv->metadata_->protein_begin();

  for (int protein_idx = 0; protein_idx <
                            (int)indiv->metadata_->proteins_count(); protein_idx++) {
    {

      Protein_7* prot  =
          indiv->metadata_->protein_next();
      if (prot->is_init_ && prot->leading_lagging==1) {
        if (lookup.find(prot->protein_start) ==
            lookup.end()) {
          lookup[prot->protein_start] = prot;
        } else {
          #ifdef __MULTI_PROMOTERS
            #if 0 == __MULTI_PROMOTERS
              if(prot->e > lookup[prot->protein_start]->e)
                lookup[prot->protein_start]->e = prot->e;
              if(prot->initial_e_ > lookup[prot->protein_start]->initial_e_)
                lookup[prot->protein_start]->initial_e_ = prot->initial_e_;
            #endif
          #else
          lookup[prot->protein_start]->e += prot->e;
          lookup[prot->protein_start]->initial_e_ += prot->initial_e_;
          #endif
          lookup[prot->protein_start]->rna_list_.insert(
              lookup[prot->protein_start]->rna_list_.end(),
              prot->rna_list_.begin(),prot->rna_list_.end());
          prot->is_init_ = false;
        }
      }
    }
  }


  #ifdef __MULTI_PROMOTERS
  #if 1 == __MULTI_PROMOTERS
  indiv->metadata_->protein_begin();

  for (int protein_idx = 0; protein_idx <
                            (int)indiv->metadata_->proteins_count(); protein_idx++) {
    Protein_7* prot  =
          indiv->metadata_->protein_next();
    
    if (prot->is_init_) {
      Rna_7* first_rna = *(prot->rna_list_.begin());
      int32_t max_length = -1;

      if (prot->leading_lagging == 1) {
        max_length = prot->protein_start > first_rna->begin ? 
                              prot->protein_start - first_rna->begin : prot->protein_start + (indiv->dna_->length_ - first_rna->begin);
      } else {
        max_length = prot->protein_start < first_rna->begin ? 
                              first_rna->begin - prot->protein_start : first_rna->begin + (indiv->dna_->length_ - prot->protein_start);
      }

      for (auto rna : prot->rna_list_) {
        int32_t rna_length = -1;

        if (prot->leading_lagging == 1) {
          max_length = prot->protein_start > rna->begin ? 
                                prot->protein_start - rna->begin : prot->protein_start + (indiv->dna_->length_ - rna->begin);
        } else {
          max_length = prot->protein_start < rna->begin ? 
                                rna->begin - prot->protein_start : rna->begin + (indiv->dna_->length_ - prot->protein_start);
        }
      
        if (rna_length > max_length) {
          first_rna = rna;
          max_length = rna_length;
        }
      }

      prot->e = first_rna->e;
      prot->initial_e_ = first_rna->e;
    }

  }
  #endif
  #endif
    

      #ifdef WITH_PERF_TRACES_PER_INDIV
      auto runtime_end = std::chrono::steady_clock::now();
      translate_protein_stop_[indiv->indiv_id] = runtime_end.time_since_epoch().count();
      translate_protein_[indiv->indiv_id] = runtime_end.time_since_epoch().count() - runtime_start.time_since_epoch().count();
    #endif
}

void ExpManager_7::compute_phenotype(int indiv_id) {
  compute_phenotype(current_individuals[indiv_id]);
}

void ExpManager_7::compute_phenotype(Individual_7* indiv) {
      #ifdef WITH_PERF_TRACES_PER_INDIV
      auto runtime_start = std::chrono::steady_clock::now();
      compute_protein_start_[indiv->indiv_id]= runtime_start.time_since_epoch().count();
    #endif

  AbstractFuzzy_7* activ_phenotype = fuzzy_factory_->get_fuzzy();
  AbstractFuzzy_7* inhib_phenotype = fuzzy_factory_->get_fuzzy();


  std::vector<Protein_7*> protein_vector;
  indiv->metadata_->protein_begin();
  for (int protein_idx = 0; protein_idx < indiv->metadata_->proteins_count(); protein_idx++) {
    Protein_7* prot = indiv->metadata_->protein_next();
    if (prot->is_init_) {
      if (prot->leading_lagging==0)
        protein_vector.push_back(prot);
    }
  }

  indiv->metadata_->protein_begin();
  for (int protein_idx = 0; protein_idx < indiv->metadata_->proteins_count(); protein_idx++) {
    Protein_7* prot = indiv->metadata_->protein_next();
    if (prot->is_init_) {
      if (prot->leading_lagging==1)
        protein_vector.push_back(prot);
    }
  }


  // for (auto prot : protein_vector) {
  //           printf("%d -- %d -- P_THIS %f %f %f %d (S %d E %d)\n",AeTime::time(),indiv->indiv_id, prot->h,prot->m,prot->w,prot->protein_start,prot->signal_,prot->inherited_);
  // }
  std::sort(protein_vector.begin(), protein_vector.end(),
            [](Protein_7*a, Protein_7*b) { return *a < *b;});
            // if (AeTime::time() >= 3) printf("%d -- %d -- Compute phenotype %ld\n",AeTime::time(),indiv_id,protein_vector.size());
  for (auto prot : protein_vector) {
    if (fabs(
        prot->w) >=
        1e-15 &&
        fabs(
            prot->h) >=
        1e-15) {

      if (prot->is_functional) {

      // if (current_individuals[indiv_id]->indiv_id==99 && AeTime::time() == 4)
      //  printf("SIMD -- Add triangle %lf %lf %lf (%lf %lf)\n",prot->m,
      //         prot->w,
      //         prot->h * prot->e,
      //         prot->h, prot->e );
        // Compute triangle points' coordinates
        
        bool verbose = false;
        // if (AeTime::time() == 447 && indiv_id==966) {
        //   verbose = true;
        // }


        if (prot->h > 0)
          activ_phenotype->add_triangle(prot->m, prot->w, prot->h *
                                                          prot->e, verbose);
        else
          inhib_phenotype->add_triangle(prot->m, prot->w, prot->h *
                                                          prot->e, verbose);

      // if (indiv_id==41 && AeTime::time() == 1)
        // printf("SIMD -- Phenotype : %lf %lf\n",activ_phenotype->get_geometric_area(),inhib_phenotype->get_geometric_area());
      }
    }
  }



  // if (AeTime::time()==3 && indiv_id == 781) {
  //   activ_phenotype->print();
  //   inhib_phenotype->print();
  // }

  activ_phenotype->clip(AbstractFuzzy_7::max,   Y_MAX);
  inhib_phenotype->clip(AbstractFuzzy_7::min, - Y_MAX);

  activ_phenotype->simplify();
  inhib_phenotype->simplify();

  if (indiv->phenotype != nullptr) {
    fuzzy_factory_->give_back(indiv->phenotype);
    indiv->phenotype = nullptr;
  }

      
  indiv->phenotype = fuzzy_factory_->get_fuzzy();
  indiv->phenotype->copy(activ_phenotype);
  indiv->phenotype->add(inhib_phenotype);
  indiv->phenotype->clip(AbstractFuzzy_7::min, Y_MIN);
  indiv->phenotype->simplify();

      // if (indiv->indiv_id==205) {
      //   printf("%p == %p\n",indiv,current_individuals[206]);
        // printf("SIMD -- %d -- Phenotype : %lf %lf :: %lf\n",indiv->indiv_id,activ_phenotype->get_geometric_area(),inhib_phenotype->get_geometric_area(),indiv->phenotype->get_geometric_area());
      // }
  //  if (indiv_id==0) printf("Meta error A %lf I %lf P %lf\n",activ_phenotype->get_geometric_area(),
  //         inhib_phenotype->get_geometric_area(), current_individuals[indiv_id]->phenotype->get_geometric_area());
  fuzzy_factory_->give_back(activ_phenotype);
  fuzzy_factory_->give_back(inhib_phenotype);

    
    #ifdef WITH_PERF_TRACES_PER_INDIV
      auto runtime_end = std::chrono::steady_clock::now();
      compute_phenotype_stop_[indiv->indiv_id] = runtime_end.time_since_epoch().count();
      compute_phenotype_[indiv->indiv_id] = runtime_end.time_since_epoch().count() - runtime_start.time_since_epoch().count();
    #endif
}

void ExpManager_7::compute_fitness(int indiv_id, double selection_pressure, int env_id) {
  compute_fitness(current_individuals[indiv_id],selection_pressure,env_id);
}

void ExpManager_7::compute_fitness(Individual_7* indiv, double selection_pressure, int env_id
#ifdef __REGUL
, SIMD_PhenotypicTargetHandler_R* pth
#endif
) {
#ifdef __REGUL
  if (pth == nullptr)
    pth = phenotypic_target_handler_;
#endif

      #ifdef WITH_PERF_TRACES_PER_INDIV
      auto runtime_start = std::chrono::steady_clock::now();
      compute_fitness_start_[indiv->indiv_id] = runtime_start.time_since_epoch().count();
    #endif
#ifdef __REGUL
  AbstractFuzzy_7* delta = fuzzy_factory_->get_fuzzy();
  delta->copy(indiv->phenotype);
  // if (indiv->indiv_id==0) printf("SIMD -- Step A -- Delta %lf\n",delta->get_geometric_area());
  
      if (pth->var_method_ == SWITCH_IN_A_LIST)
        delta->sub(pth->targets_fuzzy_by_id_[env_id]);
      else if (pth->var_method_ == ONE_AFTER_ANOTHER)
        delta->sub(pth->targets_fuzzy_by_id_[env_id]);
      

      // if (
    //     #ifdef HAVE_MPI
    //     (exp_m_->rank() == 0) && 
    //     #endif
    //     (indiv->indiv_id == 176)) {//81
    //   #ifdef HAVE_MPI
    //   FILE *fp;
    //   char filename[20];
    //   sprintf(filename,"test_log_meta_%d",exp_m_->rank());
    //   fp = fopen(filename, "a");
    //   #endif
    // // printf("SIMD -- %d -- Step B -- %d -- Delta %lf (SUB %lf)\n",indiv->indiv_id,env_id,delta->get_geometric_area(), 
    // // phenotypic_target_handler_->targets_fuzzy_by_id_[env_id]->get_geometric_area());
    //   #ifdef HAVE_MPI
    //   fprintf(fp,
    //   #else
    //   printf(
    //   #endif
        
    //     "%d -- Metaerror G %d (E %d PC %d SIG %ld RNA %d): Phenotype %e Target %e SUB %e\n",AeTime::time(),
    //   #ifdef HAVE_MPI
    //   localXtoGlobalX(indiv->indiv_id / grid_height_)*global_grid_height_+localYtoGlobalY(indiv->indiv_id % grid_height_),
    //   // indiv_id / grid_height_,indiv_id % grid_height_,
    //   // localXtoGlobalX(indiv_id / grid_height_),localYtoGlobalY(indiv_id % grid_height_),
    //   #else
    //   indiv->indiv_id,
    //   #endif
    //   env_id,
    //   indiv->metadata_->proteins_count(),((List_Metadata*)indiv->metadata_)->signal_proteins_.size(),indiv->metadata_->rna_count(),
    //   indiv->phenotype->get_geometric_area(),phenotypic_target_handler_->targets_fuzzy_by_id_[env_id]->get_geometric_area(),
    //   delta->get_geometric_area());

    //   #ifdef HAVE_MPI
    //   fclose(fp);
    //   #endif
    //   }

      if (pth->var_method_ == SWITCH_IN_A_LIST)
        indiv->metaerror_by_env_id_[env_id] = delta->get_geometric_area();
      else if (pth->var_method_ == ONE_AFTER_ANOTHER)
        indiv->metaerror_by_env_id_[env_id] += delta->get_geometric_area();

//  if (indiv_id==68 && AeTime::time() == 4) {
    // if (indiv->indiv_id==0) {
    //   printf("Delta : %lf (TARGET %d : %lf // %lf ||  PHEN %lf)\n",delta->get_geometric_area(),
    //     env_id,
    //     phenotypic_target_handler_->targets_fuzzy_[env_id]->get_geometric_area(),
    //     phenotypic_target_handler_->targets_fuzzy_by_id_[env_id]->get_geometric_area(),
    //     indiv->phenotype->get_geometric_area());
    //   delta->print();
    //   phenotypic_target_handler_->targets_fuzzy_[env_id]->print();
    //   current_individuals[indiv_id]->phenotype->print();
    // }
//     // delta->print();

//     printf("Phenotype \n");
//     current_individuals[indiv_id]->phenotype->get_geometric_area(true);
//       // current_individuals[indiv_id]->phenotype->print();

//     printf("Target %d :: ID %d\n",env_id,phenotypic_target_handler_->list_env_id_[env_id]);
//     phenotypic_target_handler_->targets_fuzzy_[env_id]->get_geometric_area(true);
//       // phenotypic_target_handler_->targets_fuzzy_[env_id]->print();
//  }

      delete delta;

      if (pth->var_method_ == SWITCH_IN_A_LIST)
        indiv->fitness_by_env_id_[env_id] = exp(
            -selection_pressure *
            ((double) indiv->metaerror_by_env_id_[env_id]));
#else
  AbstractFuzzy_7* delta = fuzzy_factory_->get_fuzzy();
  // printf("%d -- Phenotype %lf : %p (%p)\n",indiv_id,current_individuals[indiv_id]->phenotype->get_geometric_area(),
  //                                       current_individuals[indiv_id]->phenotype,delta);

  delta->copy(indiv->phenotype);
  bool verbose = false;
  // if (indiv_id==966 && AeTime::time() == 447) {
  //   verbose = true;
    
  // }
  delta->sub(target,verbose);
  
  // printf("SIMD -- %d -- Target %lf I %lf :: %lf\n",indiv_id,target->get_geometric_area(),delta->get_geometric_area(),
  //         current_individuals[indiv_id]->phenotype->get_geometric_area());

  // if (indiv_id==966 && AeTime::time() == 447) {

  //   printf("SIMD -- Delta %lf\n",delta->get_geometric_area());
  //   // delta->print();
  // }

//  if (indiv_id==68 && AeTime::time() == 4) {
//   printf("Delta SIMD : \n");
//     delta->print();
//  }


  indiv->metaerror = delta->get_geometric_area();
    // printf("Metadata %f\n",current_individuals[indiv_id]->metaerror);

  fuzzy_factory_->give_back(delta);

  indiv->fitness = exp(
      -selection_pressure *
      ((double)indiv->metaerror));
#endif

    #ifdef WITH_PERF_TRACES_PER_INDIV
      auto runtime_end = std::chrono::steady_clock::now();
      compute_fitness_stop_[indiv->indiv_id] = runtime_end.time_since_epoch().count();
      compute_fitness_[indiv->indiv_id] = runtime_end.time_since_epoch().count() - runtime_start.time_since_epoch().count();
    #endif
}


#ifdef __REGUL
void ExpManager_7::compute_network(int indiv_id, double selection_pressure) {
  compute_network(current_individuals[indiv_id],selection_pressure);
}

void ExpManager_7::compute_network(Individual_7* indiv, double selection_pressure, SIMD_PhenotypicTargetHandler_R* pth) {
  if (pth == nullptr)
    pth = phenotypic_target_handler_;

  // Allocate fitness and metarror array
  // if (pth->var_method_ == SWITCH_IN_A_LIST) {
  //   indiv->fitness_by_env_id_ = new double[pth->nb_indiv_age_];
  //   indiv->metaerror_by_env_id_ = new double[pth->nb_indiv_age_];
  // } else if (pth->var_method_ == ONE_AFTER_ANOTHER) {
  //   indiv->fitness_by_env_id_ = new double[pth->nb_env_];
  //   indiv->metaerror_by_env_id_ = new double[pth->nb_env_];

  //   for (int env_id = 0; env_id < pth->nb_env_; env_id++) {
  //     indiv->fitness_by_env_id_[env_id] = 0;
  //     indiv->metaerror_by_env_id_[env_id] = 0;
  //   }
  // }
// printf("%d -- Protein Size ==> BEFORE SIG SOLVE %d\n",indiv->indiv_id,indiv->metadata_->proteins_count());
  // Add signals
  //((List_Metadata*)indiv->metadata_)->proteins_remove_signal();
  ((List_Metadata*)indiv->metadata_)->signal_proteins_.resize(pth->signals_models_.size());
  int i = 0;
  for (auto signal : pth->signals_models_) {
    int glob_prot_idx = indiv->metadata_->proteins_count();
    indiv->metadata_->set_proteins_count(
        indiv->metadata_->proteins_count() +
        1);

                // printf("%d -- ADD SIGNAL -- Indiv %d Add prot %d\n",AeTime::time(),indiv->indiv_id,glob_prot_idx);

    Protein_7* prot = new Protein_7(signal, exp_m_);
    indiv->metadata_->protein_add(glob_prot_idx, prot);
    ((List_Metadata*)indiv->metadata_)->signal_proteins_[i] = (prot);
    i++;
  }
          // printf("%d -- Protein Size ==> AFTER SIG SOLVE %d\n",indiv->indiv_id,indiv->metadata_->proteins_count());

  indiv->metadata_->rna_begin();
  for (int i = 0; i < indiv->metadata_->rna_count(); i++) {
    Rna_7* rna = indiv->metadata_->rna_next();
    if (rna != nullptr) {
      if (rna->is_coding_) {
        rna->affinity_list.clear();
      }
    }
  }
                
        indiv->metadata_->protein_begin();
        // printf("%d -- Protein Size ==> SOLVE %d\n",indiv->indiv_id,indiv->metadata_->proteins_count());
        for (int j = 0; j < indiv->metadata_->proteins_count(); j++) {
          Protein_7* prot =
              indiv->metadata_->protein_next();
          if (prot != nullptr) {
 
  
		        for (auto&& rna: prot->rna_list_) {
			   rna->affinity_list.clear();
			}
                    }
                }


  // Set influences
  indiv->metadata_->rna_begin();
  for (int i = 0; i < indiv->metadata_->rna_count(); i++) {
    Rna_7* rna = indiv->metadata_->rna_next();
    if (rna != nullptr) {
      if (rna->is_coding_) {
        rna->nb_influences_ = 0;

        int32_t enhancer_position = rna->enhancer_position(indiv->dna_->length());
        int32_t operator_position = rna->operator_position(indiv->dna_->length());

        double enhance,operate;
        indiv->metadata_->protein_begin();
        // printf("%d -- Protein Size ==> SOLVE %d\n",indiv->indiv_id,indiv->metadata_->proteins_count());
        for (int j = 0; j < indiv->metadata_->proteins_count(); j++) {
          Protein_7* prot =
              indiv->metadata_->protein_next();
          if (prot != nullptr) {
            if (prot->is_init_ || prot->signal_) {
              enhance = rna->affinity_with_protein(
                  enhancer_position, prot, indiv,
                  exp_m_);
              operate = rna->affinity_with_protein(
                  operator_position, prot, indiv,
                  exp_m_);

              if (enhance != 0.0 || operate != 0.0) {

                rna->affinity_list.push_back(
                    AffinityFactor(prot, enhance, operate));
                prot->is_TF_ = true;

                rna->nb_influences_++;

                //  if (indiv_id==68 && AeTime::time() == 4)
      //           if (indiv->indiv_id == 0)
      //             printf("%d -- %d -- Affinity between RNA %d and Protein %d : %lf %lf\n",
      //                    AeTime::time(),
      //                            #ifdef HAVE_MPI
      // localXtoGlobalX(indiv->indiv_id / grid_height_)*global_grid_height_+localYtoGlobalY(indiv->indiv_id % grid_height_)
      // // indiv_id / grid_height_,indiv_id % grid_height_,
      // // localXtoGlobalX(indiv_id / grid_height_),localYtoGlobalY(indiv_id % grid_height_),
      // #else
      // indiv->indiv_id
      // #endif
      //                    ,rna->begin, prot->protein_start, enhance, operate);
              }
            }
          }
        }
      }
    }
  }
}

void ExpManager_7::update_network(int indiv_id, double selection_pressure) {
  update_network(current_individuals[indiv_id],selection_pressure);
}

void ExpManager_7::update_network(Individual_7* indiv, double selection_pressure, int degradation_step, int lifestep_step) {

  indiv->metadata_->protein_begin();
  for (int j = 0; j < indiv->metadata_->proteins_count(); j++) {
    Protein_7* prot =
        indiv->metadata_->protein_next();
    if (!prot->signal_) {
      if (prot->is_init_) {

        prot->delta_concentration_ = 0;

        for (auto rna: prot->rna_list_) {
          double synthesis_rate = rna->compute_synthesis_rate(indiv);

          prot->delta_concentration_ += synthesis_rate;
        }
        
        // printf("SIMD Life %d Deg %d -- Protein BEFORE DEGRADATION %d (H %d) : %lf = E %lf + D %lf (%lf) : RNAList %zu\n", lifestep_step, degradation_step,
        //          prot->protein_start,prot->inherited_, prot->e + prot->delta_concentration_,
        //          prot->e, prot->delta_concentration_, exp_m_->exp_s()->get_degradation_rate(), prot->rna_list_.size());

        prot->delta_concentration_ -=
            exp_m_->exp_s()->get_degradation_rate() * prot->e;


        // printf("SIMD Life %d Deg %d -- Protein BETWEEN DEGRADATION %d (H %d) : %lf = %lf + %lf (%lf)\n", lifestep_step, degradation_step,
        //          prot->protein_start,prot->inherited_, prot->e + prot->delta_concentration_,
        //          prot->e, prot->delta_concentration_, exp_m_->exp_s()->get_degradation_rate());

        prot->delta_concentration_ *=
            1.0 / exp_m_->exp_s()->get_nb_degradation_step();

        // printf("SIMD Life %d Deg %d -- Protein AFTER DEGRADATION %d (H %d) : %lf = %lf + %lf\n", lifestep_step, degradation_step,
        //          prot->protein_start,prot->inherited_,prot->e + prot->delta_concentration_,
        //          prot->e, prot->delta_concentration_);
      }
    }
  }


  // Apply the changes in concentrations we have just computed


//  if (indiv_id==137)
//    ((List_Metadata*)current_individuals[indiv_id]->metadata_)->proteins_print();

  indiv->metadata_->protein_begin();
  for (int j = 0; j < indiv->metadata_->proteins_count(); j++) {
    Protein_7* prot =
        indiv->metadata_->protein_next();
//    printf("SIMD -- Protein %d : %d %d\n",prot->protein_start,prot->signal_,prot->is_init_);
    if (!prot->signal_) {
      if (prot->is_init_) {
        // if (indiv_id == 70 && AeTime::time() == 1595)
        // if (prot->delta_concentration_ >= 0.000001)
        //   printf("SIMD -- Protein %d : %e + %e\n", prot->protein_start,
        //          prot->e, prot->delta_concentration_);
      //   if (indiv_id==188 && AeTime::time() == 524)
      //     printf("%d -- %d -- Protein %d : %lf + %lf\n",
      //                              AeTime::time(),
      //                            #ifdef HAVE_MPI
      // localXtoGlobalX(indiv->indiv_id / grid_height_)*global_grid_height_+localYtoGlobalY(indiv->indiv_id % grid_height_)
      // indiv_id / grid_height_,indiv_id % grid_height_,
      // localXtoGlobalX(indiv_id / grid_height_),localYtoGlobalY(indiv_id % grid_height_),
      // #else
      // indiv->indiv_id
      // #endif
      //     ,prot->protein_start,
      //       prot->e,prot->delta_concentration_);

        prot->e += prot->delta_concentration_;
      }
    }
  }
}

void ExpManager_7::evaluate_network(int indiv_id, double selection_pressure, int env_id) {
  evaluate_network(current_individuals[indiv_id],selection_pressure,env_id);
}

void ExpManager_7::evaluate_network(Individual_7* indiv, double selection_pressure, int env_id, SIMD_PhenotypicTargetHandler_R* pth) {
  if (pth == nullptr)
    pth = phenotypic_target_handler_;
    
  update_phenotype(indiv);

  //   current_individuals[indiv_id]->metadata_->protein_begin();
  // for (int j = 0; j < current_individuals[indiv_id]->metadata_->proteins_count(); j++) {
  //   Protein_7* prot =
  //       current_individuals[indiv_id]->metadata_->protein_next();

  //   if (!prot->signal_)
  //     if (prot->is_init_) {
  //       double old_e = prot->e;
  //       prot->e += prot->delta_concentration_;
  //       printf("SIMD -- Evaluate/Update %d : %lf + %lf\n", prot->protein_start,
  //                prot->e, prot->delta_concentration_);
  //     }
  // }

  // printf("Compute fitness %d : ID %d\n",env_id,phenotypic_target_handler_->list_env_id_[env_id]);
  compute_fitness(indiv,selection_pressure,env_id,pth);

  if (pth->var_method_ == SWITCH_IN_A_LIST)
    indiv->metaerror += indiv->metaerror_by_env_id_[env_id];
}

void ExpManager_7::finalize_network(int indiv_id, double selection_pressure) {
  finalize_network(current_individuals[indiv_id],selection_pressure);
}

void ExpManager_7::finalize_network(Individual_7* indiv, double selection_pressure, SIMD_PhenotypicTargetHandler_R* pth) {
   //double sum_meta = current_individuals[indiv_id]->metaerror;
   if (pth == nullptr)
    pth = phenotypic_target_handler_;

  if (pth->var_method_ == SWITCH_IN_A_LIST) {
    indiv->metaerror =
        indiv->metaerror /
        (double)(pth->nb_eval_);
  }
  else if (pth->var_method_ == ONE_AFTER_ANOTHER) {
    for (int env_id = 0; env_id < pth->nb_eval_; env_id++) {
      indiv->metaerror_by_env_id_[env_id] =
          indiv->metaerror_by_env_id_[env_id] / 10.0;
      indiv->metaerror +=
          indiv->metaerror_by_env_id_[env_id];

      indiv->fitness_by_env_id_[env_id] = exp(
          -selection_pressure *
          ((double) indiv->metaerror_by_env_id_[env_id]));
    }
    indiv->metaerror =
        indiv->metaerror / (double) (pth->nb_env_);
  }

  indiv->fitness = exp(
      -selection_pressure *
      ((double) indiv->metaerror));


  // if (indiv->indiv_id == 205) {
  //     printf("%d -- SIMD -- Finalize Network :: %lf %lf (%ld)\n",indiv->indiv_id,indiv->metaerror, indiv->fitness,
  //             pth->nb_eval_);
  //     if (current_individuals[206] != nullptr) printf("%d (L %d) -- SIMD -- Finalize Network :: %lf %lf (%ld)\n",206,current_individuals[206]->indiv_id,
  //             current_individuals[206]->metaerror, current_individuals[206]->fitness,
  //             pth->nb_eval_);

  // }
}

void ExpManager_7::solve_network(int indiv_id, double selection_pressure) {
  // printf("%d -- %d -- SOLVE FUNCTION NETWORK\n",AeTime::time(),indiv_id);
  solve_network(current_individuals[indiv_id],selection_pressure);
}

void ExpManager_7::solve_network(Individual_7* indiv, double selection_pressure, SIMD_PhenotypicTargetHandler_R* pth, bool verbose) {

  
  if (pth == nullptr)
    pth = phenotypic_target_handler_;

  // Allocate fitness and metarror array
  if (pth->var_method_ == SWITCH_IN_A_LIST) {
    indiv->fitness_by_env_id_ = new double[pth->nb_indiv_age_];
    indiv->metaerror_by_env_id_ = new double[pth->nb_indiv_age_];
  } else if (pth->var_method_ == ONE_AFTER_ANOTHER) {
    indiv->fitness_by_env_id_ = new double[pth->nb_env_];
    indiv->metaerror_by_env_id_ = new double[pth->nb_env_];

    for (int env_id = 0; env_id < pth->nb_env_; env_id++) {
      indiv->fitness_by_env_id_[env_id] = 0;
      indiv->metaerror_by_env_id_[env_id] = 0;
    }
  }

  indiv->metadata_->protein_begin();
  for (int j = 0;
       j < indiv->metadata_->proteins_count(); j++) {
    Protein_7* prot =
        indiv->metadata_->protein_next();
    if (!prot->signal_) {
      if (prot->is_init_) {
        prot->e = prot->initial_e_;
      }
    }
  }

  indiv->metaerror = 0;

  // printf("%d -- %d -- SOLVE COMPUTE_NETWORK NETWORK\n",AeTime::time(),indiv->indiv_id);

  compute_network(indiv, selection_pressure,pth);

  indiv->metaerror = 0;

      //   ((List_Metadata*)indiv->metadata_)->proteins_print(
      //   AeTime::time(),      
      //   #ifdef HAVE_MPI
      // localXtoGlobalX(indiv->indiv_id / grid_height_)*global_grid_height_+localYtoGlobalY(indiv->indiv_id % grid_height_)
      // // indiv_id / grid_height_,indiv_id % grid_height_,
      // // localXtoGlobalX(indiv_id / grid_height_),localYtoGlobalY(indiv_id % grid_height_),
      // #else
      // indiv->indiv_id
      // #endif
      // );

  // if (indiv_id==190 && AeTime::time() == 1936)
  //  ((List_Metadata*)current_individuals[indiv_id]->metadata_)->proteins_print(0);
  // printf("%d -- %d -- SOLVE BEFORE_SOLVE NETWORK\n",AeTime::time(),indiv->indiv_id);

  if (pth->var_method_ == ONE_AFTER_ANOTHER) {
    for (int16_t env_i = 0; env_i < pth->nb_env_; env_i++) {

      //Set the concentration of signals for this age
      for (auto prot1: ((List_Metadata*)indiv->metadata_)->signal_proteins_) {
        prot1->e = 0;
      }

      for (auto prot_id : pth->env_signals_list_[env_i]) {
        ((List_Metadata*)indiv->metadata_)->signal_proteins_[prot_id]->e = 0.9;
      }

      // TODO REINIT PROT ?

      for (int16_t i = 0; i < 10; i++) {

        for (int j = 0; j < 10; j++) {
          update_network(indiv,selection_pressure);
        }

        // If we have to evaluate the individual at this age
        evaluate_network(indiv,selection_pressure,env_i,pth);

                // }
      }
    }

    finalize_network(indiv,selection_pressure);
  } else {
    std::set<int>* eval = exp_m_->exp_s()->get_list_eval_step();
    // i is thus the age of the individual
      // printf("%d -- %d -- SOLVE PHEN NETWORK\n",AeTime::time(),indiv->indiv_id);
    int32_t prev_env_id = pth->list_env_id_[0];
    for (int16_t i = 0; i < pth->nb_indiv_age_; i++) {


      //Set the concentration of signals for this age
      for (auto prot1: ((List_Metadata*)indiv->metadata_)->signal_proteins_) {
        prot1->e = 0;
      }

      for (auto prot_id : pth->env_signals_list_[pth->list_env_id_[i]]) {
        ((List_Metadata*)indiv->metadata_)->signal_proteins_[prot_id]->e = 0.9;
      }

      if (pth->env_switch_probability_ < 0) {
        if (prev_env_id != pth->list_env_id_[i]) {        
          indiv->metadata_->protein_begin();
          for (int j = 0;
              j < indiv->metadata_->proteins_count(); j++) {
            Protein_7* prot =
                indiv->metadata_->protein_next();
            if (!prot->signal_) {
              if (prot->is_init_) {
                prot->e = prot->initial_e_;
              }
            }
          }
        }

        prev_env_id = pth->list_env_id_[i];
      }

      // printf("%d -- %d -- SOLVE BEFORE_UPDATE NETWORK\n",AeTime::time(),indiv->indiv_id);

      //         ((List_Metadata*)indiv->metadata_)->proteins_print(
      //   AeTime::time(),      
      //   #ifdef HAVE_MPI
      // localXtoGlobalX(indiv->indiv_id / grid_height_)*global_grid_height_+localYtoGlobalY(indiv->indiv_id % grid_height_)
      // // indiv_id / grid_height_,indiv_id % grid_height_,
      // // localXtoGlobalX(indiv_id / grid_height_),localYtoGlobalY(indiv_id % grid_height_),
      // #else
      // indiv->indiv_id
      // #endifn
      // );


      for (int j = 0; j < exp_m_->exp_s()->get_nb_degradation_step(); j++) {
        update_network(indiv,selection_pressure,j,i);
      }


      // if (indiv->indiv_id == 0) {
      // std::ofstream concentrationfile;
      // concentrationfile.open("stats/online_e7_concentration_visu.csv",std::ofstream::app);

      //   // Save concentration
      //   indiv->metadata_->protein_begin();
      //   int32_t id = 0;
      //   for (int j = 0;
      //       j < indiv->metadata_->proteins_count(); j++) {
      //     Protein_7* prot =
      //         indiv->metadata_->protein_next();
      //     if (prot->is_init_ || prot->signal_) {
      //         concentrationfile<<AeTime::time()<<","<<indiv->indiv_id<<","<<i<<","<<id<<","<<prot->e<<","<<(prot->signal_? 1 : 0)<<std::endl;
      //         id++;
      //     }
      //   }

      //   concentrationfile.close();
      // }
      // printf("%d -- %d -- SOLVE AFTER_UPDATE NETWORK\n",AeTime::time(),indiv->indiv_id);

      // if (indiv->indiv_id == 294)
      // ((List_Metadata*)indiv->metadata_)->proteins_print(
      //   i,      
      //   #ifdef HAVE_MPI
      // localXtoGlobalX(indiv->indiv_id / grid_height_)*global_grid_height_+localYtoGlobalY(indiv->indiv_id % grid_height_)
      // // indiv_id / grid_height_,indiv_id % grid_height_,
      // // localXtoGlobalX(indiv_id / grid_height_),localYtoGlobalY(indiv_id % grid_height_),
      // #else
      // indiv->indiv_id
      // #endif
      // );

  // if (indiv_id==68 && AeTime::time() == 4)
  //  ((List_Metadata*)current_individuals[indiv_id]->metadata_)->proteins_print(i+1);

      // If we have to evaluate the individual at this age
      if (eval->find(i+1) != eval->end()) {//} || (indiv_id==543 && AeTime::time() == 5895)) {// ||( (indiv_id == 70) && (AeTime::time()>=1570))) {
        evaluate_network(indiv,selection_pressure, pth->list_env_id_[i]);

/*
          printf(
                  "%d -- SIMD --  Metaerror -- Evaluate Network at %d (env id %d):: %lf %lf -- %lf\n",

      indiv->indiv_id,

                  i+1,pth->list_env_id_[i],
                         indiv->metaerror,
               indiv->metaerror_by_env_id_[pth->list_env_id_[i]],
                         pth->targets_fuzzy_by_id_[pth->list_env_id_[i]]->get_geometric_area());*/
      //   #pragma omp critical
      //   {
      //   printf("%d -- %d -- Evaluate Network %d (E %d) -- %d : \n",
      //   AeTime::time(),      
      //   #ifdef HAVE_MPI
      //   localXtoGlobalX(indiv->indiv_id / grid_height_)*global_grid_height_+localYtoGlobalY(indiv->indiv_id % grid_height_)
      //   #else
      //   indiv->indiv_id
      //   #endif
      //   ,i+1,phenotypic_target_handler_->list_env_id_[i],indiv->last_id
      //   );

      //   ((List_Metadata*)indiv->metadata_)->proteins_print(
      //   AeTime::time(),      
      //   #ifdef HAVE_MPI
      // localXtoGlobalX(indiv->indiv_id / grid_height_)*global_grid_height_+localYtoGlobalY(indiv->indiv_id % grid_height_)
      // // indiv_id / grid_height_,indiv_id % grid_height_,
      // // localXtoGlobalX(indiv_id / grid_height_),localYtoGlobalY(indiv_id % grid_height_),
      // #else
      // indiv->indiv_id
      // #endif
      // );
        
        // }
        //        if (
        // #ifdef HAVE_MPI
        // (exp_m_->rank() == 0) && 
        // #endif
        // (indiv->indiv_id == 112)) {
      //     if (AeTime::time() > 0 && verbose && (phenotypic_target_handler_->hasChanged_ || exp_m_->dna_mutator_array_[indiv->indiv_id]->hasMutate())) {
      //     #pragma omp critical
      //     {
      //     #ifdef HAVE_MPI
      // FILE *fp;
      // char filename[20];
      // sprintf(filename,"test_log_meta_%d",exp_m_->rank());
      // fp = fopen(filename, "a");
      // fprintf(fp,
      //     #else
      //     printf(
      //     #endif
      //         //   if (indiv_id==41 && AeTime::time() == 1)  {
      //           "%d -- Metaerror -- %d : %d\n",AeTime::time(),
      // #ifdef HAVE_MPI
      // localXtoGlobalX(indiv->indiv_id / grid_height_)*global_grid_height_+localYtoGlobalY(indiv->indiv_id % grid_height_),
      // // indiv_id / grid_height_,indiv_id % grid_height_,
      // // localXtoGlobalX(indiv_id / grid_height_),localYtoGlobalY(indiv_id % grid_height_),
      // #else
      // indiv->indiv_id,
      // #endif
                
      //           ((List_Metadata*)indiv->metadata_)->proteins_count());

      // if (indiv->indiv_id == 205)
      //     #ifdef HAVE_MPI
      //     fprintf(fp,
      //     #else
      //     printf(
      //     #endif
      //             "%d -- SIMD --  Metaerror -- Evaluate Network at %d :: %lf %lf -- %lf\n",
                  
      // #ifdef HAVE_MPI
      // localXtoGlobalX(indiv->indiv_id / grid_height_)*global_grid_height_+localYtoGlobalY(indiv->indiv_id % grid_height_),
      // // indiv_id / grid_height_,indiv_id % grid_height_,
      // // localXtoGlobalX(indiv_id / grid_height_),localYtoGlobalY(indiv_id % grid_height_),
      // #else
      // indiv->indiv_id,
      // #endif
                  
      //             i+1,
      //                    indiv->metaerror,
      //          indiv->metaerror_by_env_id_[phenotypic_target_handler_->list_env_id_[i]],
      //                    phenotypic_target_handler_->targets_fuzzy_by_id_[phenotypic_target_handler_->list_env_id_[i]]->get_geometric_area());


      //   #ifdef HAVE_MPI
      //   fclose(fp);
      //   #endif
      //   }

      //     }
      }
    }


    finalize_network(indiv,selection_pressure);
  }
      // printf("%d -- %d -- SOLVE END NETWORK\n",AeTime::time(),indiv->indiv_id);

  if (pth->var_method_ == SWITCH_IN_A_LIST) {
    delete [] indiv->fitness_by_env_id_;
    indiv->fitness_by_env_id_ = nullptr;
    delete [] indiv->metaerror_by_env_id_;
    indiv->metaerror_by_env_id_ = nullptr;
  } else if (pth->var_method_ == ONE_AFTER_ANOTHER) {
    delete [] indiv->fitness_by_env_id_;
    delete [] indiv->metaerror_by_env_id_;
  }
}


void ExpManager_7::recompute_network(Individual_7* indiv, double selection_pressure, SIMD_PhenotypicTargetHandler_R* pth) {

  // Allocate fitness and metarror array
  if (pth->var_method_ == SWITCH_IN_A_LIST) {
    indiv->fitness_by_env_id_ = new double[pth->nb_indiv_age_];
    indiv->metaerror_by_env_id_ = new double[pth->nb_indiv_age_];
  } else if (pth->var_method_ == ONE_AFTER_ANOTHER) {
    indiv->fitness_by_env_id_ = new double[pth->nb_env_];
    indiv->metaerror_by_env_id_ = new double[pth->nb_env_];

    for (int env_id = 0; env_id < pth->nb_env_; env_id++) {
      indiv->fitness_by_env_id_[env_id] = 0;
      indiv->metaerror_by_env_id_[env_id] = 0;
    }
  }

  indiv->metadata_->protein_begin();
  for (int j = 0;
       j < indiv->metadata_->proteins_count(); j++) {
    Protein_7* prot =
        indiv->metadata_->protein_next();
    if (!prot->signal_) {
      if (prot->is_init_) {
        prot->e = prot->initial_e_;
      }
    }
  }

  indiv->metaerror = 0;

  if (pth->var_method_ == ONE_AFTER_ANOTHER) {
    for (int16_t env_i = 0; env_i < pth->nb_env_; env_i++) {

      //Set the concentration of signals for this age
      for (auto prot1: ((List_Metadata*)indiv->metadata_)->signal_proteins_) {
        prot1->e = 0;
      }

      for (auto prot_id : pth->env_signals_list_[env_i]) {
        ((List_Metadata*)indiv->metadata_)->signal_proteins_[prot_id]->e = 0.9;
      }

      for (int16_t i = 0; i < 10; i++) {

        for (int j = 0; j < 10; j++) {
          update_network(indiv,selection_pressure);
        }

        // If we have to evaluate the individual at this age
        evaluate_network(indiv,selection_pressure,env_i,pth);
      }
    }

    finalize_network(indiv,selection_pressure,pth);
  } else {
    std::set<int>* eval = exp_m_->exp_s()->get_list_eval_step();
    // i is thus the age of the individual

    for (int16_t i = 0; i < pth->nb_indiv_age_; i++) {


      //Set the concentration of signals for this age
      for (auto prot1: ((List_Metadata*)indiv->metadata_)->signal_proteins_) {
        prot1->e = 0;
      }

      for (auto prot_id : pth->env_signals_list_[pth->list_env_id_[i]]) {
        ((List_Metadata*)indiv->metadata_)->signal_proteins_[prot_id]->e = 0.9;
      }

      for (int j = 0; j < exp_m_->exp_s()->get_nb_degradation_step(); j++) {
        update_network(indiv,selection_pressure);
      }

      if (eval->find(i+1) != eval->end()) {
        evaluate_network(indiv,selection_pressure, pth->list_env_id_[i],pth);
        // printf("Eval at %d using env %d :: %lf\n",i,pth->list_env_id_[i],indiv->metaerror);
      }
    }


    finalize_network(indiv,selection_pressure,pth);
  }

  if (pth->var_method_ == SWITCH_IN_A_LIST) {
    delete [] indiv->fitness_by_env_id_;
    delete [] indiv->metaerror_by_env_id_;
  } else if (pth->var_method_ == ONE_AFTER_ANOTHER) {
    delete [] indiv->fitness_by_env_id_;
    delete [] indiv->metaerror_by_env_id_;
  }
}

void ExpManager_7::update_phenotype( int indiv_id ) {
  update_phenotype(current_individuals[indiv_id]);
}

void ExpManager_7::update_phenotype( Individual_7* indiv ) {
  if (indiv->phenotype!=nullptr) {
    fuzzy_factory_->give_back(indiv->phenotype);
    indiv->phenotype = nullptr;
  }

  compute_phenotype(indiv);
}
#endif

void ExpManager_7::setup_individuals(double w_max, double selection_pressure) {
  if (exp_m_->sel()->fitness_func() == FITNESS_GLOBAL_SUM) {
#ifdef __REGUL
    fitness_sum_tab_ = new double[phenotypic_target_handler_->nb_env_];
    for (int env_id = 0; env_id < phenotypic_target_handler_->nb_env_; env_id++) {
      fitness_sum_tab_[env_id] = 0;
      for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++)
        fitness_sum_tab_[env_id] += previous_individuals[indiv_id]->fitness_by_env_id_[env_id];
    }
#else
    printf("Fitness local sum is not supported for Aevol (only R-Aevol)\n");
      exit(-1);
#endif
  }

  for (int indiv_id = 0; indiv_id < nb_indivs_; ++indiv_id) {
    start_stop_RNA(indiv_id);
    compute_RNA(indiv_id);
    start_protein(indiv_id);
    compute_protein(indiv_id);
    translate_protein(indiv_id, w_max);
            // ((List_Metadata*)current_individuals[indiv_id]->metadata_)->proteins_print();

#ifdef __REGUL
    solve_network(indiv_id,selection_pressure);
#else
    compute_phenotype(indiv_id);
    compute_fitness(indiv_id, selection_pressure);
#endif
  }

  for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
    previous_individuals[indiv_id] = current_individuals[indiv_id];
  }

    // for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
    //   printf("%d -- %d -- To Delete %d : Usage %d\n",AeTime::time(),exp_m_->rank(),indiv_id,previous_individuals[indiv_id]->usage_count_);
    // }
  stats_best = new Stats_7(this, AeTime::time(), true);
  write_stat();

  // Rewritting tree to avoid index mismatch

//   if (exp_m_->record_tree() && AeTime::time() %  exp_m_->tree_step() == 0) {
//         int status;
//         status = mkdir(TREE_DIR, 0755);
//         if ((status == -1) && (errno != EEXIST))
//         {
//           err(EXIT_FAILURE, "Impossible to create the directory %s", TREE_DIR);
//         }

//         char tree_file_name[50];

//         #ifdef HAVE_MPI
//         sprintf(tree_file_name, "tree/tree_" TIMESTEP_FORMAT "_%d.ae", AeTime::time(),exp_m_->rank());
//         #else
//         sprintf(tree_file_name, "tree/tree_" TIMESTEP_FORMAT ".ae", AeTime::time());
//         #endif
// std::cout << "writing 7777 tree for gen : " << AeTime::time()
//                     << '\n';
//         exp_m_->output_m()->tree()->write_to_tree_file(tree_file_name);
//   }


//   if (!exp_m_->check_simd() && AeTime::time() % exp_m_->backup_step() == 0) {
// #pragma omp single
//     {
//       for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
//         int x = indiv_id / exp_m_->world()->height();
//         int y = indiv_id % exp_m_->world()->height();
        

//         exp_m_->world()->grid(x, y)->individual()->clear_everything_except_dna_and_promoters();
//         exp_m_->world()->grid(x, y)->individual()->genetic_unit_list_nonconst().clear();
//         delete exp_m_->world()->grid(x, y)->individual();

// #ifdef __REGUL
//         Individual_R *indiv = new Individual_R(exp_m_,
//                                                          exp_m_->world()->grid(x, y)->mut_prng(),
//                                                          exp_m_->world()->grid(x, y)->stoch_prng(),
//                                                          exp_m_->exp_s()->mut_params(),
//                                                          w_max,
//                                                          exp_m_->exp_s()->min_genome_length(),
//                                                          exp_m_->exp_s()->max_genome_length(),
//                                                          false,
//                                                          indiv_id,
//                                                          "",
//                                                          0);
// #else
//         Individual *indiv = new Individual(exp_m_,
//                                            exp_m_->world()->grid(x, y)->mut_prng(),
//                                            exp_m_->world()->grid(x, y)->stoch_prng(),
//                                            exp_m_->exp_s()->mut_params(),
//                                            w_max,
//                                            exp_m_->exp_s()->min_genome_length(),
//                                            exp_m_->exp_s()->max_genome_length(),
//                                            false,
//                                            indiv_id,
//                                            "",
//                                            0);
// #endif

//         int32_t nb_blocks_ =
//             previous_individuals[indiv_id]->dna_->nb_block();
//         char *dna_string = new char[nb_blocks_ * BLOCK_SIZE];
//         memset(dna_string, 0,
//                (previous_individuals[indiv_id]->dna_->length() + 1) * sizeof(char));


//         char *to_copy = previous_individuals[indiv_id]->dna_->to_char();


//         memcpy(dna_string, to_copy,
//                (previous_individuals[indiv_id]->dna_->length() + 1) * sizeof(char));

//         dna_string[previous_individuals[indiv_id]->dna_->length()] = '\0';

//         indiv->add_GU(dna_string,
//                       previous_individuals[indiv_id]->dna_->length());
//         indiv->genetic_unit_nonconst(0).set_min_gu_length(exp_m_->exp_s()->min_genome_length());
//         indiv->genetic_unit_nonconst(0).set_max_gu_length(exp_m_->exp_s()->max_genome_length());
//         indiv->compute_statistical_data();
//         indiv->EvaluateInContext(exp_m_->world()->grid(x, y)->habitat());


//         exp_m_->world()->grid(x, y)->set_individual(indiv);
//       //   {
//       // int32_t global_x = localXtoGlobalX(x);
//       // int32_t global_y = localYtoGlobalY(y);
//       // int32_t global_id = global_x * exp_m_->exp_s()->global_grid_height() + global_y;
        
//       //   printf("%d -- R %d -- Indiv %d DNa Length %d\n",
//       //   AeTime::time(),exp_m_->rank(),
//       //   global_id,
//       //   exp_m_->world()->indiv_at(x,y)->genetic_unit_seq_length(0));
//       //   }
//       }

//       // Create missing directories
//       exp_m_->WriteDynamicFiles();

//       std::ofstream last_gener_file(LAST_GENER_FNAME,
//                                     std::ofstream::out);

//       last_gener_file << AeTime::time() << std::endl;
//       last_gener_file.close();
//     }
//   }

}

void ExpManager_7::evaluate(Individual_7* indiv, double w_max, double selection_pressure, SIMD_PhenotypicTargetHandler_R* pth) {
  // printf("Evaluate %d : %lf %lf\n",indiv->indiv_id,w_max,selection_pressure);
    start_stop_RNA(indiv);
    compute_RNA(indiv);
    start_protein(indiv);
    compute_protein(indiv);
    translate_protein(indiv, w_max);

  
#ifdef __REGUL
          // printf("%d -- %d -- SOLVE EVALUATE NETWORK\n",AeTime::time(),indiv->indiv_id);
    solve_network(indiv,selection_pressure, pth, false);
#else
    compute_phenotype(indiv);
    compute_fitness(indiv, selection_pressure);
#endif
}

void ExpManager_7::run_a_step(double w_max, double selection_pressure) {
  #ifdef WITH_PERF_TRACES
      auto runtime_start = std::chrono::steady_clock::now();
      auto runtime_end = std::chrono::steady_clock::now();
      auto runtime_counter = runtime_end.time_since_epoch().count() - runtime_start.time_since_epoch().count();
      runtime_counter = -1;
      auto stats_tree_backup_start = std::chrono::steady_clock::now();
      auto stats_tree_backup_end = std::chrono::steady_clock::now();
      auto stats_tree_backup_counter = stats_tree_backup_end.time_since_epoch().count() - stats_tree_backup_start.time_since_epoch().count();

      auto selection_start = std::chrono::steady_clock::now();
      auto selection_end = std::chrono::steady_clock::now();
      auto selection_counter = selection_end.time_since_epoch().count() - selection_start.time_since_epoch().count();

      auto evaluate_start = std::chrono::steady_clock::now();
      auto evaluate_end = std::chrono::steady_clock::now();
      auto evalute_counter = evaluate_end.time_since_epoch().count() - evaluate_start.time_since_epoch().count();

      auto sort_start = std::chrono::steady_clock::now();
      auto sort_end = std::chrono::steady_clock::now();
      auto sort_counter = sort_end.time_since_epoch().count() - sort_start.time_since_epoch().count();

      auto tree_start = std::chrono::steady_clock::now();
      auto tree_end = std::chrono::steady_clock::now();
      auto tree_counter = tree_end.time_since_epoch().count() - tree_start.time_since_epoch().count();

      auto clean_start = std::chrono::steady_clock::now();
      auto clean_end = std::chrono::steady_clock::now();
      auto clean_counter = clean_end.time_since_epoch().count() - clean_start.time_since_epoch().count();

      auto stats_start = std::chrono::steady_clock::now();
      auto stats_end = std::chrono::steady_clock::now();
      auto stats_counter = clean_end.time_since_epoch().count() - clean_start.time_since_epoch().count();


      auto tree_backup_start = std::chrono::steady_clock::now();
      auto tree_backup_end = std::chrono::steady_clock::now();
      auto tree_backup_counter = clean_end.time_since_epoch().count() - clean_start.time_since_epoch().count();

      #ifdef HAVE_MPI
      auto exch_fitness_start = std::chrono::steady_clock::now();
      auto exch_fitness_end = std::chrono::steady_clock::now();
      auto exch_indiv_start = std::chrono::steady_clock::now();
      auto exch_indiv_end = std::chrono::steady_clock::now();
      auto exch_selection_start = std::chrono::steady_clock::now();
      auto exch_selection_end = std::chrono::steady_clock::now();
      auto eval_start = std::chrono::steady_clock::now();
      auto eval_end = std::chrono::steady_clock::now();
      auto eval_counter = eval_end.time_since_epoch().count() - eval_start.time_since_epoch().count();
      auto exch_indiv_counter = exch_indiv_end.time_since_epoch().count() - exch_indiv_start.time_since_epoch().count();
      auto exch_selection_counter = selection_end.time_since_epoch().count() - selection_start.time_since_epoch().count();
      auto exch_fitness_counter = exch_fitness_end.time_since_epoch().count() - exch_fitness_start.time_since_epoch().count();
      #endif
  #endif

  #ifdef WITH_PERF_TRACES_PER_INDIV
    for (int i = 0; i < nb_indivs_; i++) {
        apply_mutation[i] = -1;
        compute_rna_[i] = -1;
        start_protein_[i] = -1;
        compute_protein_[i] = -1;
        translate_protein_[i] = -1;
        compute_phenotype_[i] = -1;
        compute_fitness_[i] = -1;
        total_[i] = -1;


        allocate_individual_start_[i] = -1;
        apply_mutation_start_[i] = -1;
        compute_rna_start_[i] = -1;
        start_protein_start_[i] = -1;
        compute_protein_start_[i] = -1;
        translate_protein_start_[i] = -1;
        compute_phenotype_start_[i] = -1;
        compute_fitness_start_[i] = -1;
        allocate_individual_stop_[i] = -1;
        apply_mutation_stop_[i] = -1;
        compute_rna_stop_[i] = -1;
        start_protein_stop_[i] = -1;
        compute_protein_stop_[i] = -1;
        translate_protein_stop_[i] = -1;
        compute_phenotype_stop_[i] = -1;
        compute_fitness_stop_[i] = -1;
        total_start_[i] = -1;
        total_stop_[i] = -1;
        
        
        omp_tid_[i] = -1;
    }
  #endif


  #ifdef WITH_PERF_TRACES
  #ifdef HAVE_MPI
#pragma omp single copyprivate(runtime_start,exch_fitness_start,selection_start)
  #else
#pragma omp single copyprivate(runtime_start,selection_start)
  #endif
#else
#pragma omp single
#endif
  {
    nb_clones_ = 0;
    mutant_list_.clear();
    mutant_list_.reserve(nb_indivs_);

    if (exp_m_->sel()->fitness_func() == FITNESS_GLOBAL_SUM) {
#ifdef __REGUL
      fitness_sum_tab_ = new double[phenotypic_target_handler_->nb_env_];
        for (int env_id = 0; env_id < phenotypic_target_handler_->nb_env_; env_id++) {
          fitness_sum_tab_[env_id] = 0;
          for (int indiv_id = 0; indiv_id < exp_m_->nb_indivs(); indiv_id++)
            fitness_sum_tab_[env_id] += previous_individuals[indiv_id]->fitness_by_env_id_[env_id];
        }
#else
      printf("Fitness local sum is not supported for Aevol (only R-Aevol)\n");
      exit(-1);
#endif
    }
    #ifdef WITH_PERF_TRACES
      runtime_start = std::chrono::steady_clock::now();
      selection_start = std::chrono::steady_clock::now();
      #ifdef HAVE_MPI
      exch_fitness_start = std::chrono::steady_clock::now();
      #endif
    #endif
  }

#ifdef HAVE_MPI
  #ifdef WITH_PERF_TRACES
#pragma omp single copyprivate(exch_fitness_end,exch_selection_start)
#else
#pragma omp single
#endif
{
  send_broadcast_border();
  // MPI_Barrier(MPI_COMM_WORLD);
  recv_broadcast_border();


  #ifdef WITH_PERF_TRACES
    exch_fitness_end = std::chrono::steady_clock::now();
    exch_selection_start =  std::chrono::steady_clock::now();
  #endif
}
#endif

#pragma omp single
{
  if (exp_m_->sel()->selection_scope() == SCOPE_GLOBAL) {
    #ifdef HAVE_MPI
    printf("SELECTION_SCOPE GLOBAL is not available for MPI yet\n");
    exit(-1);
    #endif

    double *  local_fit_array   = new double[nb_indivs_];
    double *  probs             = new double[nb_indivs_];
    double    sum_local_fit     = 0.0;
    
    for (int loc_indiv_id = 0; loc_indiv_id < nb_indivs_; loc_indiv_id++) {
      local_fit_array[loc_indiv_id] =
              previous_individuals[loc_indiv_id]->fitness;
      sum_local_fit += local_fit_array[loc_indiv_id];
    }

    for (int loc_indiv_id = 0; loc_indiv_id < nb_indivs_; loc_indiv_id++) {
      probs[loc_indiv_id] = local_fit_array[loc_indiv_id]/sum_local_fit;
    }
    
    int32_t* nb_offsprings = new int32_t[nb_indivs_];
    exp_m_->world()->grid(0,0)->reprod_prng_simd_->multinomial_drawing(nb_offsprings, probs, nb_indivs_, nb_indivs_);

    int index = 0;
    
    int32_t x = 0;
    int32_t y = 0;

    for (int32_t loc_indiv = 0; loc_indiv < nb_indivs_; loc_indiv++) {
      for (int32_t j = 0; j < nb_offsprings[loc_indiv]; j++) {
        x = index / grid_height_;
        y = index % grid_height_;
        exp_m_->next_generation_reproducer_[x*grid_height_+y] = loc_indiv;
        index++;
      }
    }

    delete [] nb_offsprings;
    delete [] local_fit_array;
    delete [] probs;
  }
}
  
#pragma omp for schedule(dynamic)
  for (int indiv_id = 0; indiv_id < nb_indivs_; indiv_id++) {

    if (!exp_m_->check_simd()) {
      // printf("BEFORE Selection %d\n",indiv_id);
      // if (indiv_id == 1017) {
      //         for (int32_t jj = 0; jj < grid_height_; jj++) {
      //           printf("1017 -- FITNESS AT BORDER %d : %e\n",jj,previous_individuals[jj]->fitness);
      //         }
      // }
      selection(indiv_id);

    }

#ifndef HAVE_MPI
    do_mutation(indiv_id,w_max,selection_pressure);
#endif

  // #if defined(__REGUL) or defined(HAVE_MPI)
}

#ifdef HAVE_MPI
  #ifdef WITH_PERF_TRACES
#pragma omp single copyprivate(exch_selection_end,exch_indiv_start,exch_indiv_end,eval_start)
#else
#pragma omp single
#endif
{
    #ifdef WITH_PERF_TRACES
      exch_selection_end = std::chrono::steady_clock::now();
      exch_indiv_start = std::chrono::steady_clock::now();
    #endif
  // MPI_Barrier(MPI_COMM_WORLD);
  Fetch_Remote_Individual();
  // MPI_Barrier(MPI_COMM_WORLD);

    #ifdef WITH_PERF_TRACES
      exch_indiv_end = std::chrono::steady_clock::now();
      eval_start = std::chrono::steady_clock::now();
    #endif
}

  #pragma omp for schedule(dynamic)
  for (int indiv_id = 0; indiv_id < nb_indivs_; indiv_id++) {
    do_mutation(indiv_id,w_max,selection_pressure);
}

  for (auto ite = to_delete_individuals_.begin(); ite != to_delete_individuals_.end(); ite++) {
    delete *ite;
  }
  to_delete_individuals_.clear();
#endif

// for (int index = 0; index < mutant_list_.size(); index++) {
// printf("Mutant %d\n",mutant_list_[index]);
// }
#ifdef WITH_PERF_TRACES
#ifdef HAVE_MPI
#pragma omp single copyprivate(selection_end,sort_start,eval_end)
#else
#pragma omp single copyprivate(selection_end,sort_start)
#endif
{
  #ifdef HAVE_MPI
        eval_end = std::chrono::steady_clock::now();
  #endif
  selection_end = std::chrono::steady_clock::now();
  sort_start = std::chrono::steady_clock::now();
}
#endif

#ifdef __OMP_LIST_SORT
#pragma omp single
{
  #if 0 == __OMP_LIST_SORT
  sort(mutant_list_.begin(), mutant_list_.end(), [this](int a, int b) {
        return exp_m_->dna_mutator_array_[a]->length_ > exp_m_->dna_mutator_array_[b]->length_;
    });
  #elif 1 == __OMP_LIST_SORT
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

  std::shuffle (mutant_list_.begin(), mutant_list_.end(), std::default_random_engine(seed));
  #endif
}
#endif


#ifdef WITH_PERF_TRACES
#pragma omp single copyprivate(sort_end,evaluate_start)
{
  sort_end = std::chrono::steady_clock::now();
  evaluate_start = std::chrono::steady_clock::now();
}
#endif

#ifdef __GNUC__
#pragma omp for schedule(dynamic,1)
#else
#pragma omp for schedule(monotonic:dynamic,1)
#endif
for (int index = 0; index < mutant_list_.size(); index++) {
  int32_t indiv_id = mutant_list_[index];

    #ifdef WITH_PERF_TRACES_PER_INDIV
      auto indiv_runtime_start = std::chrono::steady_clock::now();
      omp_tid_[indiv_id] = omp_get_thread_num();
    #endif

      // if (indiv_id == 205)
      //   printf("ALREADY V00 ADDED  %p == %p\n",
      //           current_individuals[indiv_id],current_individuals[current_individuals[indiv_id]->added_id]);

  // for (int indiv_id = 0; indiv_id < nb_indivs_; indiv_id++) {
  // #endif

      // printf("%d -- %d -- EVALUATE %d (%d) %d :: %d (%d) (%p)\n",AeTime::time(),indiv_id,current_individuals[indiv_id]->last_id, 
      //                 current_individuals[indiv_id]->last_id == indiv_id,
      //                 current_individuals[indiv_id]->indiv_id, exp_m_->dna_mutator_array_[indiv_id]->hasMutate(),
      //                 (phenotypic_target_handler_->hasChanged_ && current_individuals[indiv_id]->last_id == indiv_id),
      //                 current_individuals[indiv_id] );
    if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate() 
    #ifdef HAVE_MPI
    || exp_m_->dna_mutator_array_[indiv_id]->isAtBorder()
    #endif
    )
      apply_mutations(indiv_id,w_max,selection_pressure);
}

      // #pragma omp single
      // {
      //   std::map<int64_t,int8_t> count;
      //   for (int32_t id = 0; id < nb_indivs_; id++) {
      //     count[id] = 0;
      //     // printf("%d\n",id);
      //   }

      //   int32_t nb_it = 0;

      //   for (auto id : mutant_list_) {
      //     // printf("Count %d : %d\n",mutant_list_[id],current_individuals[mutant_list_[id]]->indiv_id);
      //     count[current_individuals[id]->indiv_id] = count[current_individuals[id]->indiv_id] + 1;
      //     nb_it++;
      //   }
      //   printf("Mutant list %d\n",mutant_list_.size(),nb_it);


      //   for (auto id : mutant_list_) {
      //     if (count[current_individuals[id]->indiv_id] > 1) {
      //       printf("Error multiple %d : %d\n",current_individuals[id]->indiv_id,count[current_individuals[id]->indiv_id]);
      //       exit(-1);
      //     }
      //   }
      // }


#ifdef __GNUC__
#pragma omp for schedule(dynamic,1)
#else
#pragma omp for schedule(monotonic:dynamic,1)
#endif
for (int index = 0; index < mutant_list_.size(); index++) {
  int32_t indiv_id = mutant_list_[index];

    #ifdef WITH_PERF_TRACES_PER_INDIV
      auto indiv_runtime_start = std::chrono::steady_clock::now();
      omp_tid_[indiv_id] = omp_get_thread_num();
    #endif

    // printf("Indiv %d :: %d\n",indiv_id,exp_m_->dna_mutator_array_[indiv_id]->hasMutate());
    if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate() 
    #ifdef __REGUL
    || phenotypic_target_handler_->hasChanged_//|| (current_individuals[indiv_id]->last_id == indiv_id) // phenotypic_target_handler_->hasChanged_ && 
    #endif
    #ifdef HAVE_MPI
    || (current_individuals[indiv_id]->is_at_border_)
    #endif
    ) {

      #ifdef __REGUL
      if (exp_m_->dna_mutator_array_[indiv_id]->hasMutate()) {
      #endif

      // if (indiv_id == 205)
      //   printf("ALREADY ADDED  %p == %p\n",
      //           current_individuals[indiv_id],current_individuals[current_individuals[indiv_id]->added_id]);
        // printf("NoCLone %d\n",indiv_id);

      opt_prom_compute_RNA(indiv_id);
      start_protein(indiv_id);
      compute_protein(indiv_id);
      translate_protein(indiv_id, w_max);
      #ifdef __REGUL
      //((List_Metadata*)current_individuals[indiv_id]->metadata_)->proteins_remove_signal();
      if (exp_m_->exp_s()->get_with_heredity())
        ((List_Metadata*)current_individuals[indiv_id]->metadata_)->add_inherited_proteins(current_individuals[indiv_id]);

      } else {
        // printf("Inherited %d\n",indiv_id);
        // ((List_Metadata*)current_individuals[indiv_id]->metadata_)->proteins_print();
        ((List_Metadata*)current_individuals[indiv_id]->metadata_)->reinit_proteins_inherited(current_individuals[indiv_id]);
        // ((List_Metadata*)current_individuals[indiv_id]->metadata_)->proteins_print();

      }
      #endif
#ifdef __REGUL
      solve_network(indiv_id,selection_pressure);
#else
      compute_phenotype(indiv_id);
      compute_fitness(indiv_id, selection_pressure);
#endif


#ifdef DEBUG_MPI
      #ifdef HAVE_MPI
      int32_t i = localXtoGlobalX(indiv_id / exp_m_->grid_height());
      int32_t j = localYtoGlobalY(indiv_id % exp_m_->grid_height());
          printf("T %d -- %d -- Individual %d (G %d): has FINALLY computed from scratch brder %d :: %e\n",AeTime::time(),exp_m_->rank(),
                indiv_id,
                i*global_grid_width_+j,
                current_individuals[indiv_id]->metadata_->promoter_count(),
                current_individuals[indiv_id]->fitness);

                FILE *fp;
      char filename[20];
      sprintf(filename,"test_log_%d",exp_m_->rank());
      fp = fopen(filename, "a");
#else
      int32_t i = indiv_id / exp_m_->grid_height();
      int32_t j = indiv_id % exp_m_->grid_height();

          printf("T %d -- Individual %d (G %d): has FINALLY computed from scratch brder %d :: %e\n",AeTime::time(),
                indiv_id,
                i*grid_width_+j,
                current_individuals[indiv_id]->metadata_->promoter_count(),
                current_individuals[indiv_id]->fitness);
                FILE *fp;

      char filename[20];
      sprintf(filename,"test_log_0");
      fp = fopen(filename, "a");
#endif
      fprintf(fp,"%d -- %d -- Promoter list LEADING : ",AeTime::time(),indiv_id);
      for (auto prot : ((List_Metadata*)current_individuals[indiv_id]->metadata_)->promoters_list_[LEADING]) {
        fprintf(fp,"%d ",prot.pos);
      }
      fprintf(fp,"\n");

      fprintf(fp,"%d -- %d -- Promoter list LAGGING : ",AeTime::time(),indiv_id);
      for (auto prot : ((List_Metadata*)current_individuals[indiv_id]->metadata_)->promoters_list_[LAGGING]) {
        fprintf(fp,"%d ",prot.pos);
      }
      fprintf(fp,"\n");
      fclose(fp);
      #endif
        }

    #ifdef WITH_PERF_TRACES_PER_INDIV
      auto indiv_runtime_end = std::chrono::steady_clock::now();
      total_[indiv_id] = indiv_runtime_end.time_since_epoch().count() - indiv_runtime_start.time_since_epoch().count();
      total_start_[indiv_id] = indiv_runtime_start.time_since_epoch().count();
      total_stop_[indiv_id] = indiv_runtime_end.time_since_epoch().count();
    #endif
}



// Traces
#ifdef WITH_PERF_TRACES_PER_INDIV
#pragma omp single
{
    std::ofstream perf_traces_file_;
            perf_traces_file_.open("simd_perf_traces.csv", std::ofstream::app);

                            std::ofstream gantt_traces_file_;
            gantt_traces_file_.open("gantt_perf_traces.csv", std::ofstream::app);
            
            for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
                perf_traces_file_ << AeTime::time() << "," << indiv_id << "," << apply_mutation[indiv_id] 
                      << "," << start_protein_[indiv_id] << "," << compute_protein_[indiv_id] 
                      << "," << translate_protein_[indiv_id] << "," << compute_phenotype_[indiv_id] 
                      << "," << compute_fitness_[indiv_id] << "," << total_[indiv_id]  << std::endl;
  

                gantt_traces_file_ << AeTime::time() << "," << indiv_id << "," << omp_tid_[indiv_id] << "," 
                      << ",TOTAL," << total_start_[indiv_id] 
                      << "," << total_stop_[indiv_id] << "," << exp_m_->dna_mutator_array_[indiv_id]->length_
                      << "," << current_individuals[indiv_id]->dna_->length_ 
                      << "," << current_individuals[indiv_id]->metadata_->rna_count()
                      << "," << current_individuals[indiv_id]->metadata_->proteins_count() << std::endl;

                gantt_traces_file_ << AeTime::time() << "," << indiv_id << "," << omp_tid_[indiv_id] << "," 
                      << ",ALLOCATE," << allocate_individual_start_[indiv_id] 
                      << "," << allocate_individual_stop_[indiv_id] << "," << exp_m_->dna_mutator_array_[indiv_id]->length_
                      << "," << current_individuals[indiv_id]->dna_->length_ 
                      << "," << current_individuals[indiv_id]->metadata_->rna_count()
                      << "," << current_individuals[indiv_id]->metadata_->proteins_count() << std::endl;

                gantt_traces_file_ << AeTime::time() << "," << indiv_id << "," << omp_tid_[indiv_id] << "," 
                      << ",APPLY_MUTATION," << apply_mutation_start_[indiv_id] 
                      << "," << apply_mutation_stop_[indiv_id] << "," << exp_m_->dna_mutator_array_[indiv_id]->length_
                      << "," << current_individuals[indiv_id]->dna_->length_ 
                      << "," << current_individuals[indiv_id]->metadata_->rna_count()
                      << "," << current_individuals[indiv_id]->metadata_->proteins_count() << std::endl;

                gantt_traces_file_ << AeTime::time() << "," << indiv_id << "," << omp_tid_[indiv_id] << "," 
                      << ",COMPUTE_RNA," << compute_rna_start_[indiv_id] 
                      << "," << compute_rna_stop_[indiv_id] << "," << exp_m_->dna_mutator_array_[indiv_id]->length_
                      << "," << current_individuals[indiv_id]->dna_->length_ 
                      << "," << current_individuals[indiv_id]->metadata_->rna_count()
                      << "," << current_individuals[indiv_id]->metadata_->proteins_count() << std::endl;

                gantt_traces_file_ << AeTime::time() << "," << indiv_id << "," << omp_tid_[indiv_id] << "," 
                      << ",START_PROTEIN," << start_protein_start_[indiv_id] 
                      << "," << start_protein_stop_[indiv_id] << "," << exp_m_->dna_mutator_array_[indiv_id]->length_
                      << "," << current_individuals[indiv_id]->dna_->length_ 
                      << "," << current_individuals[indiv_id]->metadata_->rna_count()
                      << "," << current_individuals[indiv_id]->metadata_->proteins_count() << std::endl;

                gantt_traces_file_ << AeTime::time() << "," << indiv_id << "," << omp_tid_[indiv_id] << "," 
                      << ",COMPUTE_PROTEIN," << compute_protein_start_[indiv_id] 
                      << "," << compute_protein_stop_[indiv_id] << "," << exp_m_->dna_mutator_array_[indiv_id]->length_
                      << "," << current_individuals[indiv_id]->dna_->length_ 
                      << "," << current_individuals[indiv_id]->metadata_->rna_count()
                      << "," << current_individuals[indiv_id]->metadata_->proteins_count() << std::endl;

                gantt_traces_file_ << AeTime::time() << "," << indiv_id << "," << omp_tid_[indiv_id] << "," 
                      << ",TRANSLATE_PROTEIN," << translate_protein_start_[indiv_id] 
                      << "," << translate_protein_stop_[indiv_id] << "," << exp_m_->dna_mutator_array_[indiv_id]->length_
                      << "," << current_individuals[indiv_id]->dna_->length_ 
                      << "," << current_individuals[indiv_id]->metadata_->rna_count()
                      << "," << current_individuals[indiv_id]->metadata_->proteins_count() << std::endl;

                gantt_traces_file_ << AeTime::time() << "," << indiv_id << "," << omp_tid_[indiv_id] << "," 
                      << ",COMPUTE_PHENOTYPE," << compute_phenotype_start_[indiv_id] 
                      << "," << compute_phenotype_stop_[indiv_id] << "," << exp_m_->dna_mutator_array_[indiv_id]->length_
                      << "," << current_individuals[indiv_id]->dna_->length_ 
                      << "," << current_individuals[indiv_id]->metadata_->rna_count()
                      << "," << current_individuals[indiv_id]->metadata_->proteins_count() << std::endl;

                gantt_traces_file_ << AeTime::time() << "," << indiv_id << "," << omp_tid_[indiv_id] << "," 
                      << ",COMPUTE_FITNESS," << compute_fitness_start_[indiv_id] 
                      << "," << compute_fitness_stop_[indiv_id] << "," << exp_m_->dna_mutator_array_[indiv_id]->length_
                      << "," << current_individuals[indiv_id]->dna_->length_ 
                      << "," << current_individuals[indiv_id]->metadata_->rna_count()
                      << "," << current_individuals[indiv_id]->metadata_->proteins_count() << std::endl;

                apply_mutation_start_[indiv_id] = -1;
                compute_rna_start_[indiv_id] = -1;
                start_protein_start_[indiv_id] = -1;
                compute_protein_start_[indiv_id] = -1;
                translate_protein_start_[indiv_id] = -1;
                compute_phenotype_start_[indiv_id] = -1;
                compute_fitness_start_[indiv_id] = -1;
                
                apply_mutation_stop_[indiv_id] = -1;
                compute_rna_stop_[indiv_id] = -1;
                start_protein_stop_[indiv_id] = -1;
                compute_protein_stop_[indiv_id] = -1;
                translate_protein_stop_[indiv_id] = -1;
                compute_phenotype_stop_[indiv_id] = -1;
                compute_fitness_stop_[indiv_id] = -1;
                
                apply_mutation[indiv_id] = -1;
                compute_rna_[indiv_id] = -1;
                start_protein_[indiv_id] = -1;
                compute_protein_[indiv_id] = -1;
                translate_protein_[indiv_id] = -1;
                compute_phenotype_[indiv_id] = -1;
                compute_fitness_[indiv_id] = -1;
                total_[indiv_id] = -1;
                total_start_[indiv_id] = -1;
                total_stop_[indiv_id] = -1;
                omp_tid_[indiv_id] = -1;
            }
            perf_traces_file_.close();
            gantt_traces_file_.close();
}
#endif

#ifdef WITH_PERF_TRACES
#pragma omp single copyprivate(evaluate_end,tree_start)
{
  evaluate_end = std::chrono::steady_clock::now();
  tree_start = std::chrono::steady_clock::now();
}
#endif

#pragma omp for schedule(dynamic)
  for (int indiv_id = 0; indiv_id < nb_indivs_; indiv_id++) {
    current_individuals[indiv_id]->is_at_border_ = false;
    if (exp_m_->record_tree()) {
      int x = indiv_id / exp_m_->world()->height();
      int y = indiv_id % exp_m_->world()->height();

      // printf("End replication %d : %p == Prev %p\n",indiv_id,current_individuals[indiv_id],previous_individuals[indiv_id]);
      auto* eindiv =
          new EndReplicationEvent(current_individuals[indiv_id], x, y);
      // Tell observers the replication is finished
      exp_m_->tree()->update_end_replication(eindiv);
      delete eindiv;
    }
    
    current_individuals[indiv_id]->first_to_add = true;
  }

#ifdef WITH_PERF_TRACES
#ifdef HAVE_MPI
#pragma omp single copyprivate(runtime_counter,eval_counter,exch_indiv_counter,exch_selection_counter,exch_fitness_counter,stats_tree_backup_start,tree_end,clean_start)
#else
#pragma omp single copyprivate(runtime_counter,tree_end,clean_start)
#endif
#else
#pragma omp single
#endif
  {
    #ifdef WITH_PERF_TRACES
      runtime_end = std::chrono::steady_clock::now();
      runtime_counter = runtime_end.time_since_epoch().count() - runtime_start.time_since_epoch().count();
      
      #ifdef HAVE_MPI
      eval_counter = eval_end.time_since_epoch().count() - eval_start.time_since_epoch().count();
      exch_indiv_counter = exch_indiv_end.time_since_epoch().count() - exch_indiv_start.time_since_epoch().count();
      exch_selection_counter = exch_selection_end.time_since_epoch().count() - exch_selection_start.time_since_epoch().count();
      exch_fitness_counter = exch_fitness_end.time_since_epoch().count() - exch_fitness_start.time_since_epoch().count();
      #endif

      stats_tree_backup_start = std::chrono::steady_clock::now();
    #endif

    if (exp_m_->record_tree())
      exp_m_->tree()->update_end_generation();
    
    #ifdef WITH_PERF_TRACES
      tree_end =  std::chrono::steady_clock::now();
      clean_start = std::chrono::steady_clock::now();
    #endif

    #ifdef HAVE_MPI
    for (auto ite = individual_to_fetch_at_border_.begin(); ite != individual_to_fetch_at_border_.end(); ite++) {
      free( (*ite).dna_ );
      free( (*ite).lead_prom_pos );
      free( (*ite).lead_prom_error );
      free( (*ite).lag_prom_pos );
      free( (*ite).lag_prom_error );
    }

    individual_to_fetch_at_border_.clear();

    for (auto ite = individual_sent_at_border_.begin(); ite != individual_sent_at_border_.end(); ite++) {
      delete [] (*ite).lead_prom_pos;
      delete [] (*ite).lead_prom_error;

      delete [] (*ite).lag_prom_pos;
      delete [] (*ite).lag_prom_error;
    }

    individual_sent_at_border_.clear();
    #endif


  }

#pragma omp for schedule(static)
    for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {

      bool toDelete = false;

#pragma omp critical(indiv_list)
      {
        // printf("%d -- %d -- To Delete %d : Usage %p\n",AeTime::time(),exp_m_->rank(),indiv_id,previous_individuals[indiv_id]);

        previous_individuals[indiv_id]->usage_count_--;
        if (previous_individuals[indiv_id]->usage_count_ == 0) {
          toDelete = true;
        }
      }


      if (toDelete) {
        #ifdef HAVE_MPI
        to_delete_individuals_.push_back(previous_individuals[indiv_id]);
        #else
        delete previous_individuals[indiv_id];
        #endif
      }
      

      previous_individuals[indiv_id] = current_individuals[indiv_id];
      current_individuals[indiv_id]  = nullptr;


    }


    #pragma omp single
    {
      if (AeTime::time() % DnaFactory::cleanup_step == 0) dna_factory_->reduce_space(this);
    }
#ifdef WITH_PERF_TRACES
#pragma omp single copyprivate(clean_end,stats_start)
  {
      clean_end = std::chrono::steady_clock::now();
      stats_start = std::chrono::steady_clock::now();
  }
#endif

#ifndef AEVOL_NO_STATS
#pragma omp single
  {
      double best_fitness = previous_individuals[0]->fitness;
  int idx_best = 0;

  for (int indiv_id = 1; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
    if (previous_individuals[indiv_id]->fitness > best_fitness) {
      idx_best = indiv_id;
      best_fitness = previous_individuals[indiv_id]->fitness;
  
    }
    
  }
    write_stat();



  }
#endif


    #ifdef WITH_PERF_TRACES
#pragma omp single copyprivate(stats_end,tree_backup_start)
  {
      stats_end = std::chrono::steady_clock::now();
      tree_backup_start = std::chrono::steady_clock::now();
  }
#endif

#pragma omp single
    {
      if (exp_m_->record_light_tree()) {
        exp_m_->output_m()->light_tree()->update_tree(AeTime::time(),
                                                      previous_individuals);

        if (AeTime::time() % exp_m_->backup_step() == 0) {
          std::cout << "writing light tree for gen : " << AeTime::time()
                    << '\n';
          exp_m_->output_m()->write_light_tree(AeTime::time());
        }
      }

      if (exp_m_->record_tree() && AeTime::time() %  exp_m_->tree_step() == 0) {
        int status;
        status = mkdir(TREE_DIR, 0755);
        if ((status == -1) && (errno != EEXIST))
        {
          err(EXIT_FAILURE, "Impossible to create the directory %s", TREE_DIR);
        }

        char tree_file_name[50];

        #ifdef HAVE_MPI
        sprintf(tree_file_name, "tree/tree_" TIMESTEP_FORMAT "_%d.ae", AeTime::time(),exp_m_->rank());
        #else
        sprintf(tree_file_name, "tree/tree_" TIMESTEP_FORMAT ".ae", AeTime::time());
        #endif
std::cout << "writing 7777 tree for gen : " << AeTime::time()
                    << '\n';
        exp_m_->output_m()->tree()->write_to_tree_file(tree_file_name);
      }

      #ifdef HAVE_MPI
      individual_to_fetch_at_border_.clear();
      #endif
    }


  if (!exp_m_->check_simd() && AeTime::time() % exp_m_->backup_step() == 0) {
#pragma omp single
    {
      for (int indiv_id = 0; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
        int x = indiv_id / exp_m_->world()->height();
        int y = indiv_id % exp_m_->world()->height();
        

        exp_m_->world()->grid(x, y)->individual()->clear_everything_except_dna_and_promoters();
        exp_m_->world()->grid(x, y)->individual()->genetic_unit_list_nonconst().clear();
        delete exp_m_->world()->grid(x, y)->individual();

#ifdef __REGUL
        Individual_R *indiv = new Individual_R(exp_m_,
                                                         exp_m_->world()->grid(x, y)->mut_prng(),
                                                         exp_m_->world()->grid(x, y)->stoch_prng(),
                                                         exp_m_->exp_s()->mut_params(),
                                                         w_max,
                                                         exp_m_->exp_s()->min_genome_length(),
                                                         exp_m_->exp_s()->max_genome_length(),
                                                         false,
                                                         indiv_id,
                                                         "",
                                                         0);
#else
        Individual *indiv = new Individual(exp_m_,
                                           exp_m_->world()->grid(x, y)->mut_prng(),
                                           exp_m_->world()->grid(x, y)->stoch_prng(),
                                           exp_m_->exp_s()->mut_params(),
                                           w_max,
                                           exp_m_->exp_s()->min_genome_length(),
                                           exp_m_->exp_s()->max_genome_length(),
                                           false,
                                           indiv_id,
                                           "",
                                           0);
#endif

        int32_t nb_blocks_ =
            previous_individuals[indiv_id]->dna_->nb_block();
        char *dna_string = new char[nb_blocks_ * BLOCK_SIZE];
        memset(dna_string, 0,
               (previous_individuals[indiv_id]->dna_->length() + 1) * sizeof(char));


        char *to_copy = previous_individuals[indiv_id]->dna_->to_char();


        memcpy(dna_string, to_copy,
               (previous_individuals[indiv_id]->dna_->length() + 1) * sizeof(char));

        dna_string[previous_individuals[indiv_id]->dna_->length()] = '\0';

        indiv->add_GU(dna_string,
                      previous_individuals[indiv_id]->dna_->length());
        indiv->genetic_unit_nonconst(0).set_min_gu_length(exp_m_->exp_s()->min_genome_length());
        indiv->genetic_unit_nonconst(0).set_max_gu_length(exp_m_->exp_s()->max_genome_length());
        indiv->compute_statistical_data();
        indiv->EvaluateInContext(exp_m_->world()->grid(x, y)->habitat());


        exp_m_->world()->grid(x, y)->set_individual(indiv);
      //   {
      // int32_t global_x = localXtoGlobalX(x);
      // int32_t global_y = localYtoGlobalY(y);
      // int32_t global_id = global_x * exp_m_->exp_s()->global_grid_height() + global_y;
        
      //   printf("%d -- R %d -- Indiv %d DNa Length %d\n",
      //   AeTime::time(),exp_m_->rank(),
      //   global_id,
      //   exp_m_->world()->indiv_at(x,y)->genetic_unit_seq_length(0));
      //   }
      }
        
      //printf("Position_7 BUF %d / %d :: %e\n",best_indiv->indiv_id/exp_m_->world()->height(),best_indiv->indiv_id%exp_m_->world()->height(),best_indiv->fitness);
      //exp_m_->update_best();
      //int32_t idx_bb = exp_m_->world()->x_best*exp_m_->world()->height()+ exp_m_->world()->y_best;
      double best_fitness = previous_individuals[0]->fitness;
      int idx_best = 0;
      for (int indiv_id = 1; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
        if (previous_individuals[indiv_id]->fitness > best_fitness) {
          idx_best = indiv_id;
          best_fitness = previous_individuals[indiv_id]->fitness;
      
        }
      }
      best_indiv = previous_individuals[idx_best];

      exp_m_->world()->x_best = best_indiv->indiv_id / exp_m_->world()->height();
      exp_m_->world()->y_best = best_indiv->indiv_id % exp_m_->world()->height();
      // printf("Position_7 BUF %d / %d :: %e // %d\n",previous_individuals[idx_bb]->indiv_id/exp_m_->world()->height(),
      //                           previous_individuals[idx_bb]->indiv_id%exp_m_->world()->height(),previous_individuals[idx_bb]->fitness,
      //                           previous_individuals[idx_bb]->dna_->length_);

      // Create missing directories
      exp_m_->WriteDynamicFiles();

      std::ofstream last_gener_file(LAST_GENER_FNAME,
                                    std::ofstream::out);

      last_gener_file << AeTime::time() << std::endl;
      last_gener_file.close();
    }
  }

#ifdef WITH_PERF_TRACES
#pragma omp single copyprivate(tree_backup_end)
  {
      tree_backup_end = std::chrono::steady_clock::now();
  }
#endif

  #pragma omp single
  {
#ifdef WITH_PERF_TRACES
    clean_end = std::chrono::steady_clock::now();

    stats_tree_backup_end = std::chrono::steady_clock::now();
    stats_tree_backup_counter = stats_tree_backup_end.time_since_epoch().count() - stats_tree_backup_start.time_since_epoch().count();

      selection_counter = selection_end.time_since_epoch().count() - selection_start.time_since_epoch().count();

      evalute_counter = evaluate_end.time_since_epoch().count() - evaluate_start.time_since_epoch().count();
      sort_counter = sort_end.time_since_epoch().count() - sort_start.time_since_epoch().count();
      tree_counter = tree_end.time_since_epoch().count() - tree_start.time_since_epoch().count();

      clean_counter = clean_end.time_since_epoch().count() - clean_start.time_since_epoch().count();

      stats_counter = stats_end.time_since_epoch().count() - stats_start.time_since_epoch().count();

      clean_counter = tree_backup_end.time_since_epoch().count() - tree_backup_start.time_since_epoch().count();
    std::ofstream perf_traces_file_;

#ifdef HAVE_MPI

    char filename[256];
    sprintf(filename,"mpi_perf_traces_%d.csv",exp_m_->rank());

    perf_traces_file_.open(filename, std::ofstream::app);
    perf_traces_file_<< AeTime::time() << "," << runtime_counter << "," << stats_tree_backup_counter << "," <<
                     exch_fitness_counter << "," << selection_counter << "," <<
                     exch_indiv_counter << "," << eval_counter << "," << sort_counter <<"," << evalute_counter << "," << 
                     tree_counter << "," <<clean_counter<<","<<stats_counter<<","<<tree_backup_counter<<"," <<
                     nb_clones_ <<std::endl;
    perf_traces_file_.close();
#else

    perf_traces_file_.open("sm_perf_traces.csv", std::ofstream::app);
    perf_traces_file_<< AeTime::time() << "," << runtime_counter << "," << stats_tree_backup_counter << "," <<
                       selection_counter << "," << sort_counter <<"," << evalute_counter << "," << 
                       tree_counter << "," <<clean_counter<<","<<stats_counter<<","<<tree_backup_counter<<"," <<
                       nb_clones_ <<std::endl;
    perf_traces_file_.close();
#endif
#endif
  }
}

void ExpManager_7::write_stat(bool just_best, bool non_coding) {
  #ifdef PROGENY_STATS
  std::map<int32_t,int32_t> progeny_per_indiv;

  for (int indiv_id = 1; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
    int32_t parent_id = exp_m_->next_generation_reproducer_[indiv_id];

    if ( progeny_per_indiv.find(parent_id) == progeny_per_indiv.end() ) {
      progeny_per_indiv[parent_id]=1;
    } else {
      progeny_per_indiv[parent_id]++;
    }
  }

  std::string file = "stats/progeny.csv";
  std::string file_rep = "stats/reproducer_progeny.csv";

  std::ofstream progeny;
  progeny.open(file,std::ofstream::app);
  std::ofstream progeny_rep;
  progeny_rep.open(file_rep,std::ofstream::app);
  progeny_rep<<AeTime::time()<<","<<progeny_per_indiv.size()<<std::endl;
  progeny_rep.close();

  for (auto prog : progeny_per_indiv) {
    progeny<<AeTime::time()<<","<<prog.first<<","<<prog.second<<std::endl;
  }

  progeny.close();
  #endif


  stats_best->reinit(AeTime::time());

  // Search for the best
  double best_fitness = previous_individuals[0]->fitness;
  int idx_best = 0;
  for (int indiv_id = 1; indiv_id < (int) exp_m_->nb_indivs(); indiv_id++) {
    if (previous_individuals[indiv_id]->fitness > best_fitness) {
      idx_best = indiv_id;
      best_fitness = previous_individuals[indiv_id]->fitness;
  
    }

  
    #ifdef __REGUL
int32_t nb_activators = 0;
  int32_t nb_operators = 0;
  double mean_activator_activity = 0.0;
  double mean_operator_activity = 0.0;
 previous_individuals[indiv_id]->metadata_->rna_begin();
      for (int i = 0; i < previous_individuals[indiv_id]->metadata_->rna_count(); i++) {
        Rna_7*rna = previous_individuals[indiv_id]->metadata_->rna_next();
        if (rna != nullptr) {
          if (rna->is_coding_) {
            for (auto affinity: rna->affinity_list) {
              if (affinity.enhancer_factor > 0.0)
              {
                nb_activators++;
                mean_activator_activity += affinity.enhancer_factor;
              }

              if (affinity.operator_factor > 0.0)
              {
                nb_operators++;
                mean_operator_activity += affinity.operator_factor;
              }
            }
          }
        }
      }

      #endif
    // if (nb_activators+nb_operators > 0)
    //   printf("%d -- Nb Link %d (%d %d) (%lf %lf)\n",indiv_id,nb_activators+nb_operators,nb_activators,nb_operators,mean_activator_activity,mean_operator_activity);
  }
  best_indiv = previous_individuals[idx_best];

// if (!just_best)  {
  best_indiv->reset_stats();
    best_indiv->metadata_->rna_begin();
    
    // ANNOTATE_SITE_BEGIN(mainloop);
    for (int prom_idx = 0;
         prom_idx < (int)best_indiv->metadata_->rna_count(); prom_idx++) {
          //  ANNOTATE_ITERATION_TASK(searchMain);
      Rna_7* rna =
          best_indiv->metadata_->rna_next();
    if (rna != nullptr) {
      if (rna->is_coding_)
        best_indiv->nb_coding_RNAs++;
      else
        best_indiv->nb_non_coding_RNAs++;
    }
  }


  for (int i = 0; i < best_indiv->metadata_->proteins_count(); i++) {
    Protein_7*prot = best_indiv->metadata_->proteins(i);
    if (prot != nullptr) {
      if (prot->is_functional) {
        best_indiv->nb_func_genes++;
      } else {
        best_indiv->nb_non_func_genes++;
      }
      if (prot->h > 0) {
        best_indiv->nb_genes_activ++;
      } else {
        best_indiv->nb_genes_inhib++;
      }
    }
  }

  stats_best->write(best_indiv,non_coding);
  // }
}



void ExpManager_7::write_stat(Stats_7* stats, Individual_7* indiv, int32_t generation, bool non_coding) {
  if (generation == -1) stats->reinit(AeTime::time());
  else stats->reinit(generation);

  indiv->reset_stats();
  indiv->metadata_->rna_begin();
  for (int i = 0; i < indiv->metadata_->rna_count(); i++) {
    Rna_7*rna = indiv->metadata_->rna_next();
    if (rna != nullptr) {
      if (rna->is_coding_)
        indiv->nb_coding_RNAs++;
      else
        indiv->nb_non_coding_RNAs++;
    }
  }


  for (int i = 0; i < indiv->metadata_->proteins_count(); i++) {
    Protein_7*prot = indiv->metadata_->proteins(i);
    if (prot != nullptr) {
      if (prot->is_functional) {
        indiv->nb_func_genes++;
      } else {
        indiv->nb_non_func_genes++;
      }
      if (prot->h > 0) {
        indiv->nb_genes_activ++;
      } else {
        indiv->nb_genes_inhib++;
      }
    }
  }

  stats->write(indiv,non_coding);
}

void ExpManager_7::check_dna() {

  int x, y;
  for (int i = 0; i < (int) exp_m_->nb_indivs(); i++) {
    x = i / exp_m_->world()->height();
    y = i % exp_m_->world()->height();

    for (int dna_pos = 0; dna_pos < current_individuals[i]->dna_->length(); dna_pos++) {
      if (exp_m_->world()->grid(x, y)->individual()->genetic_unit(
          0).dna()->data()[dna_pos] !=
          current_individuals[i]->dna_->data_[dna_pos]) {

        printf("Check DNA indiv %d %d %d --- NB Mutation %ld\n",i,current_individuals[i]->dna_->length(),
            exp_m_->world()->grid(x, y)->individual()->genetic_unit(
            0).dna()->length(),exp_m_->dna_mutator_array_[i]->mutation_list_.size());

        printf("Divergence between classic DNA and SIMD DNA %d %d at pos %d\n",
               exp_m_->world()->grid(x, y)->individual()->genetic_unit(
                   0).dna()->data()[dna_pos],
               current_individuals[i]->dna_->data_[dna_pos],
               dna_pos);
        for (auto mute : exp_m_->dna_mutator_array_[i]->mutation_list_) {
          printf("Mutation type %d\n",mute->type());
        }

        break;
      }
    }
  }
}

void ExpManager_7::check_individual(int i, int x, int y) {
  exp_m_->world()->grid(x, y)->set_individual(exp_m_->world()->grid(x, y)->old_one);
  exp_m_->world()->grid(x, y)->old_one->Reevaluate();

  printf("%d %d %d -- ",i,x,y);

  // printf(
  //     "Nb RNA SIMD/CPU %ud/%ld Protein %ud/%ld Metaerror %f/%f Fitness %e/%e DNA Size %d/%d\n",
  //     previous_individuals[i]->metadata_->rna_count(),
  //     exp_m_->world()->grid(x, y)->individual()->rna_list().size(),
  //     previous_individuals[i]->metadata_->proteins_count(),
  //     exp_m_->world()->grid(x, y)->individual()->protein_list().size(),
  //     previous_individuals[i]->metaerror,
  //     exp_m_->world()->grid(x, y)->individual()->dist_to_target_by_feature(
  //         METABOLISM),
  //     previous_individuals[i]->fitness,
  //     exp_m_->world()->grid(x, y)->individual()->fitness(),
  //     previous_individuals[i]->dna_->length(),
  //     exp_m_->world()->grid(x, y)->individual()->genetic_unit(
  //         0).seq_length());

  int idx = 0;

  // for (auto rna : exp_m_->world()->grid(x, y)->old_one->rna_list()) {
  //   printf("RNA CPU %d Start %d Stop %d Leading/Lagging %d Length %d\n", idx,
  //          rna->promoter_pos(), rna->last_transcribed_pos(), rna->strand(), rna->transcript_length());
  //   idx++;
  // }
  idx = 0;
  for (idx = 0; idx < (previous_individuals[i]->metadata_->promoter_count()); idx++) {
    if (previous_individuals[i]->metadata_->promoters(idx) != nullptr)
      printf("Promoters found at %d\n",
             previous_individuals[i]->metadata_->promoters(idx)->pos);
  }

  idx = 0;
  // for (idx = 0; idx < (previous_individuals[i]->metadata_->rna_count()); idx++) {
  //   if (previous_individuals[i]->metadata_->rnas(idx) != nullptr) {
  //     printf("RNA SIMD %d Start %d Stop %d Leading/Lagging %d Length %d\n", idx,
  //            previous_individuals[i]->metadata_->rnas(idx)->begin,
  //            previous_individuals[i]->metadata_->rnas(idx)->end,
  //            previous_individuals[i]->metadata_->rnas(idx)->leading_lagging,
  //            previous_individuals[i]->metadata_->rnas(idx)->length);
  //   }
  // }

  int prot_cpt_b=0;
  idx = 0;
  for (auto prot : exp_m_->world()->grid(x, y)->old_one->protein_list()) {
    bool found = false;

    for (int pidx = 0; pidx <
                       (int)previous_individuals[i]->metadata_->proteins_count(); pidx++) {
      if (previous_individuals[i]->metadata_->proteins(pidx)->is_init_) {
        if ((previous_individuals[i]->metadata_->proteins(pidx)->e == prot->concentration()) &&
            (previous_individuals[i]->metadata_->proteins(pidx)->protein_end == prot->last_STOP_base_pos())) {
          found = true;
          break;
        }
      }
    }

    if (!found) {
      printf("Proteins CPU %d Start %d (end %d stop %d) Length %d Leading/Lagging %d M/W/H %f/%f/%f Func %d -- Concentration %f\n",
             idx,
             prot->first_translated_pos(), prot->last_translated_pos(), prot->last_STOP_base_pos(),
             prot->length(), prot->strand(),
             prot->mean(), prot->width(), prot->height(), prot->is_functional(), prot->concentration());
    }
    idx++;
  }

  for (int idx = 0; idx <
                    (int)previous_individuals[i]->metadata_->proteins_count(); idx++) {
    if (previous_individuals[i]->metadata_->proteins(idx)->is_init_) {


      bool found = false;

      for (auto prot : exp_m_->world()->grid(x, y)->old_one->protein_list()) {
        if ((previous_individuals[i]->metadata_->proteins(idx)->e ==  prot->concentration()) &&
            (previous_individuals[i]->metadata_->proteins(idx)->protein_end ==  prot->last_STOP_base_pos())) {
          found = true;
          break;
        }
      }

      if (!found)
        printf("Proteins SIMD %d Start %d (end %d) Length %d Leading/Lagging %d M/W/H %f/%f/%f Func %d -- Concentration %f\n", idx,
            previous_individuals[i]->metadata_->proteins(idx)->protein_start,
            previous_individuals[i]->metadata_->proteins(idx)->protein_end,
            previous_individuals[i]->metadata_->proteins(idx)->protein_length,
            previous_individuals[i]->metadata_->proteins(idx)->leading_lagging,
            previous_individuals[i]->metadata_->proteins(idx)->m,
            previous_individuals[i]->metadata_->proteins(idx)->w,
            previous_individuals[i]->metadata_->proteins(idx)->h,
            previous_individuals[i]->metadata_->proteins(idx)->is_functional,
            previous_individuals[i]->metadata_->proteins(idx)->e
        );
      prot_cpt_b++;
    }
  }

}

void ExpManager_7::check_struct() {
  for (int i = 0; i < (int) exp_m_->nb_indivs(); i++) {

  }
}

void ExpManager_7::check_result() {


#pragma omp single
  {
    printf("Check results !!!!!\n");
    int x, y;

    for (int i = 0; i < (int) exp_m_->nb_indivs(); i++) {
      //if (i != 905) continue;

      x = i / exp_m_->world()->height();
      y = i % exp_m_->world()->height();

      auto* ae6_indiv          = exp_m_->world()->grid(x, y)->individual();
      auto* ae7_indiv          = previous_individuals[i];
      auto* ae7_indiv_metadata = ae7_indiv->metadata_;

      double fit_1 = ae6_indiv->dist_to_target_by_feature(
          METABOLISM);
      double fit_2 = ae7_indiv->metaerror;
      // ((List_Metadata*)ae7_indiv->metadata_)->proteins_remove_signal();

      
      float i_fit_1 = roundf(fit_1 * 100);
      float i_fit_2 = roundf(fit_2 * 100);


      int count_prot = 0;

      for (int pidx = 0; pidx < ae7_indiv_metadata->proteins_count(); pidx++) {
        if (ae7_indiv_metadata->proteins(pidx)->is_init_) {
          count_prot++;
        }
      }

      int count_rna_cpu = 0;
            int count_rna_7 = 0;

    ae7_indiv->metadata_->rna_begin();
    
    // ANNOTATE_SITE_BEGIN(mainloop);
    for (int prom_idx = 0;
         prom_idx < (int)ae7_indiv->metadata_->rna_count(); prom_idx++) {
                  Rna_7*prom =
          ae7_indiv->metadata_->rna_next();

          if (prom->is_init_) count_rna_7++;
         }

      for (auto rna : ae6_indiv->rna_list()) {
        if (rna->transcript_length() >= 0) {
          count_rna_cpu++;
        }
      }


//      int idx = 0, fidx = 0;
//      for (auto prot : ae6_indiv->protein_list()) {
//        bool found = false;
//        fidx = 0;
//
//        for (int pidx = 0; pidx < ae7_indiv_metadata->proteins_count(); pidx++) {
//          if (ae7_indiv_metadata->proteins(pidx)->is_init_) {
//            if ((ae7_indiv_metadata->proteins(pidx)->e ==
//                 prot->concentration()) &&
//                (ae7_indiv_metadata->proteins(pidx)->protein_end ==
//                 prot->last_STOP_base_pos())) {
//              if ((ae7_indiv_metadata->proteins(pidx)
//                       ->protein_length == prot->length()) &&
//                  (ae7_indiv_metadata->proteins(pidx)
//                       ->protein_start == prot->first_translated_pos())) {
//                found = true;
//                fidx  = pidx;
//                break;
//              } else {
//                fidx = pidx;
//              }
//            }
//          }
//        }
//
//        if (!found) {
//          printf("==================-------------------------======================\n");
//          printf(
//              "Proteins CPU %d Start %d (end %d stop %d) Length %d "
//              "Leading/Lagging %d M/W/H %f/%f/%f Func %d -- Concentration %f\n",
//              idx, prot->first_translated_pos(), prot->last_translated_pos(),
//              prot->last_STOP_base_pos(), prot->length(), prot->strand(),
//              prot->mean(), prot->width(), prot->height(),
//              prot->is_functional(), prot->concentration());
//
//          if (fidx < ae7_indiv_metadata->proteins_count())
//            printf("Proteins SIMD %d Start %d (end %d) Length %d Leading/Lagging %d M/W/H %f/%f/%f Func %d -- Concentration %f\n",
//                fidx, ae7_indiv_metadata->proteins(fidx)->protein_start,
//                ae7_indiv_metadata->proteins(fidx)->protein_end,
//                ae7_indiv_metadata->proteins(fidx)->protein_length,
//                ae7_indiv_metadata->proteins(fidx)->leading_lagging,
//                ae7_indiv_metadata->proteins(fidx)->m,
//                ae7_indiv_metadata->proteins(fidx)->w,
//                ae7_indiv_metadata->proteins(fidx)->h,
//                ae7_indiv_metadata->proteins(fidx)->is_functional,
//                ae7_indiv_metadata->proteins(fidx)->e);
//          printf("==================-------------------------======================\n");
//        }
//        idx++;
//      }
//
//
//      for (int j = 0; j < ae7_indiv->dna_->length(); j++) {
//        if (ae7_indiv->dna_->data_[j] !=
//            ae6_indiv->genetic_unit(0).dna()->data()[j]) {
//          printf("%ld -- %d -- DNA is different at %d !!!\n", AeTime::time(), i, j);
//
//          exit(-1);
//        }
//      }


      int prot_size = (int)ae6_indiv->protein_list().size();
      if ((((count_rna_7 != count_rna_cpu) ||
           (count_prot != prot_size)) ||
           (i_fit_1 != i_fit_2)) && ae7_indiv->dna_->length() > 300) {


        printf(
            "X-X-ERROR -- %ld -- Individual %d  (%llu %d) -- %d %d --(P %d / %d): Metaerror (CPU/GPU) : %e/%e || Fitness (CPU/GPU) : %e/%e \n",
            AeTime::time(), i, ae6_indiv->id(), ae7_indiv->indiv_id,
            x,y, ae6_indiv->parent_id_,
               ae7_indiv->parent_id,
               ae6_indiv->dist_to_target_by_feature(METABOLISM),
               ae7_indiv->metaerror, ae6_indiv->fitness(), ae7_indiv->fitness);

        printf(
            "Nb RNA SIMD/CPU %ud/%ld Protein %ud/%ld Metaerror %f/%f Fitness %e/%e DNA Size %d/%d\n",
               count_rna_7, ae6_indiv->rna_list().size(),
            count_prot, ae6_indiv->protein_list().size(),
               ae7_indiv->metaerror,
               ae6_indiv->dist_to_target_by_feature(METABOLISM),
               ae7_indiv->fitness, ae6_indiv->fitness(), ae7_indiv->dna_->length(), ae6_indiv->genetic_unit(
                0).seq_length());

        printf("Promoters LEADING : ");
        for (auto& prom : ((List_Metadata*)ae7_indiv_metadata)->promoters_list_[LEADING]) {
          printf("%d ",prom.pos);
        }
        printf("\n");
        printf("Promoters LAGGING : ");
        for (auto& prom : ((List_Metadata*)ae7_indiv_metadata)->promoters_list_[LAGGING]) {
          printf("%d ",prom.pos);
        }
        printf("\n");

        int idx = 0;

        for (auto rna : ae6_indiv->rna_list()) {
          bool found = false;
              ae7_indiv->metadata_->rna_begin();

             for (int prom_idx = 0;
         prom_idx < (int)ae7_indiv->metadata_->rna_count(); prom_idx++) {
                    Rna_7*prom =
            ae7_indiv->metadata_->rna_next();

            if (prom->is_init_) {
            if ((rna->promoter_pos() ==
                  prom->begin) &&
                (rna->transcript_length() ==
                  prom->length)){
              found = true;
              break;
            }
          }
         }

//          if (i == 392) found = false;

          if (!found)
            printf("RNA CPU %d Start %d Stop %d Leading/Lagging %d Length %d Basal %lf\n", idx,
                   rna->promoter_pos(), rna->last_transcribed_pos(), rna->strand(),
                   rna->transcript_length(),rna->basal_level());
          idx++;
        }
    ae7_indiv->metadata_->rna_begin();

         for (int prom_idx = 0;
         prom_idx < (int)ae7_indiv->metadata_->rna_count(); prom_idx++) {
                    Rna_7*prom =
            ae7_indiv->metadata_->rna_next();
          bool found = false;
         for (auto rna : ae6_indiv->rna_list()) {
           if ((rna->promoter_pos() == prom->begin) &&
               (rna->transcript_length() ==
                prom->length)) {
             found = true;
             break;
           }
         }

        //  if (i == 392) found = false;

          if (!found)
            printf("RNA SIMD %d Start %d Stop %d Leading/Lagging %d Length %d  Basal %lf\n", idx, prom->begin,
                   prom->end,
                   prom->leading_lagging,
                   prom->length,
                   prom->e);

        }


        idx = 0;

        for (auto prot : ae6_indiv->protein_list()) {
          bool found = false;

          ae7_indiv->metadata_->protein_begin();
          for (int protein_idx = 0; protein_idx <
                                    (int)ae7_indiv->metadata_->proteins_count(); protein_idx++) {
           Protein_7* protb = ae7_indiv->metadata_->protein_next();
           if ((protb->protein_start == prot->first_translated_pos()) &&
               (protb->protein_length ==
                prot->length()) &&
               (protb->e ==
                prot->concentration())) {
             found = true;
             break;
           }
         }

          // if (i == 392) found = false;
  #ifdef __REGUL
            if (!found && !((Protein_R*)prot)->is_signal()) {
  #else
            if (!found) {
  #endif
              printf("Proteins CPU %d Start %d (end %d stop %d) Length %d Leading/Lagging %d M/W/H %f/%f/%f Func %d -- Concentration % RNA : \n",
                    idx,
                    prot->first_translated_pos(), prot->last_translated_pos(),
                    prot->last_STOP_base_pos(),
                    prot->length(), prot->strand(),
                    prot->mean(), prot->width(), prot->height(), prot->is_functional(),
                    prot->concentration());

              // for (auto rna : prot->rna_list()) {
              //   printf("[%d => %d]\n",rna->promoter_pos(),rna->last_transcribed_pos());
              // }

            //  ae7_indiv->metadata_->rna_begin();

            //  for (int prom_idx = 0;
            //     prom_idx < (int)ae7_indiv->metadata_->rna_count(); prom_idx++) {
            //     Rna_7*prom =
            //     ae7_indiv->metadata_->rna_next();

            //     if (prom->is_init_ && prom->start_prot.size() > 0) {
            //       printf("RNA %d => %d-- Protein Start : ",prom->pos,prom->end);
            //       for (auto start : (prom->start_prot))
            //         printf("%d ",start);
            //       printf("\n");
            //     }
            //   }

            }

            idx++;
          }

        int prot_cpt_b = 0;
        for (int idx = 0; idx < ae7_indiv_metadata->proteins_count(); idx++) {
          if (ae7_indiv_metadata->proteins(idx)->is_init_) {


            bool found = false;
            for (auto prot : ae6_indiv->protein_list()) {
              if ((ae7_indiv_metadata->proteins(idx)->protein_start == prot->first_translated_pos()) &&
                  (ae7_indiv_metadata->proteins(idx)->protein_length ==
                    prot->length()) &&
                  (ae7_indiv_metadata->proteins(idx)->e ==
                    prot->concentration())) {
                found = true;
                break;
              }
            }
            // if (i == 392) found = false;

            //for (idx = 0; idx < (int) (current_individuals[i]->proteins.size()); idx++) {
            if (!found 
            #ifdef __REGUL
            && !ae7_indiv_metadata->proteins(idx)->signal_
            #endif
            ) {
              printf("Proteins SIMD %d Start %d (end %d) Length %d Leading/Lagging %d M/W/H %f/%f/%f Func %d -- Concentration %e\n",
                     idx, ae7_indiv_metadata->proteins(idx)->protein_start,
                     ae7_indiv_metadata->proteins(idx)->protein_end,
                     ae7_indiv_metadata->proteins(idx)->protein_length,
                     ae7_indiv_metadata->proteins(idx)->leading_lagging,
                     ae7_indiv_metadata->proteins(idx)->m,
                     ae7_indiv_metadata->proteins(idx)->w,
                     ae7_indiv_metadata->proteins(idx)->h,
                     ae7_indiv_metadata->proteins(idx)->is_functional,
                     ae7_indiv_metadata->proteins(idx)->e
              );

              // for (auto rna: previous_individuals[i]
              //                    ->metadata_->proteins(idx)
              //                    ->rna_list_) {
              //   printf("[%d => %d]\n", rna->begin, rna->end);
              // }
            }
            prot_cpt_b++;

          }
        }

//        printf("Start prot LEAD : ");
//        for (int pidx = 0; pidx < (int) (ae7_indiv_metadata->rna_count());
//             pidx++) {
//          if (ae7_indiv_metadata->rnas(pidx)->leading_lagging == 0) {
//            for (int pos : ae7_indiv_metadata->rnas(pidx)->start_prot) {
//              printf("%d ",pos);
//            }
//          }
//        }
//        printf("\n");

//        printf("Start prot LAG : ");
//        for (int pidx = 0; pidx < (int) (ae7_indiv_metadata->rna_count());
//             pidx++) {
//          if (ae7_indiv_metadata->rnas(pidx)->leading_lagging == 1) {
//            for (int pos : ae7_indiv_metadata->rnas(pidx)->start_prot) {
//              printf("%d ",pos);
//            }
//          }
//        }
//        printf("\n");

        // exp_m_->world()->grid(x, y)->individual()->phenotype()->print();

        // previous_individuals[i]->phenotype->print();

        exit(-1);
      }

    }

    printf("Generation %ld is replicated with SIMD without diff\n", AeTime::time());
  }

}

}
