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
//                              Includes
// =================================================================
#ifndef AEVOL_EXP_MANAGER_H_
#define AEVOL_EXP_MANAGER_H_

#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <list>

#include "AeTime.h"
#include "JumpingMT.h"
#include "ExpSetup.h"
#include "OutputManager.h"
#include "World.h"
#include "Observer.h"
#include "ObservableEvent.h"
#include "SaveWorld.h"

namespace aevol {
// =================================================================
//                          Class declarations
// =================================================================
class ExpManager_7;
class DnaMutator;
/// Allows for high-level experiment management. (This is Aevol's top-level class.)
///
/// An experiment manager allows one to... manage an experiment.
/// It owns a population and an experimental_setup that can be loaded from a
/// pair of aevol binary files (pop and exp_setup)
class ExpManager : public Observer {
 public:
  // =======================================================================
  //                                Constructors
  // =======================================================================
  ExpManager();

  // =======================================================================
  //                                Destructors
  // =======================================================================
  virtual ~ExpManager() noexcept;

  // =======================================================================
  //                                 Algorithms
  // =======================================================================

  // =======================================================================
  //                           Accessors: getters
  // =======================================================================
  ExpSetup* exp_s() const { return exp_s_; }
  Selection* sel() const { return exp_s()->sel(); }
  OutputManager* output_m() const { return output_m_; }
  bool quit_signal_received() const { return quit_signal_received_; }
  double selection_pressure() const { return sel()->selection_pressure(); }
  #ifdef HAVE_MPI
  int nb_rank() const { return exp_s()->nb_rank(); }
  int32_t rank() const { return rank_; }
  void set_rank(int32_t rank ) { rank_ = rank; }
  #endif

  // Spatial structure
  World* world() const { return world_; }
  World* non_const_world() { return world_; }
  int16_t grid_width() const { return world()->width(); }
  int16_t grid_height() const { return world()->height(); }
  GridCell*** grid() const { return world()->grid(); }

  // Global settings
  double repl_HT_detach_rate() const { return exp_s()->repl_HT_detach_rate();}

  // The ability to own a plasmid is a property of the individuals (allow_plasmids_) because it is used during mutations.
  // However, the experimental setup's member variable with_plasmids_ indicates whether plasmids are used
  // because the replication and loading/writting processes need this information.
  // For now when plasmids are used each individual has one and only one plasmid (so these variables should always be
  // equal), this may change in the future, though.
  // Member variable with_plasmids_HT_ has been removed because the ability to transfer is evolvable and may thus depend
  // on the plasmid itself

  bool with_plasmids() const { return exp_s()->with_plasmids(); }
  double prob_plasmid_HT() const { return exp_s()->prob_plasmid_HT(); }
  double tune_donor_ability() const { return exp_s()->tune_donor_ability(); }
  double tune_recipient_ability() const { return exp_s()->tune_recipient_ability(); }
  bool swap_GUs() const { return exp_s()->swap_GUs(); }
  bool with_secretion() const { return exp_s()->with_secretion(); }
  double secretion_contrib_to_fitness() const { return exp_s()->secretion_contrib_to_fitness(); }
  double secretion_cost() const { return exp_s()->secretion_cost(); }

  // Accessors to population stuff
  std::list<Individual*> indivs() const { return world()->indivs(); }
  std::list<std::pair<Individual*, ReplicationReport*>> indivs_annotated()
      const;
  int32_t nb_indivs() const { return world()->nb_indivs(); }
  Individual* best_indiv() const { return world()->best_indiv(); }
  Individual* indiv_by_id(int32_t id) const;
  Individual* indiv_by_rank(int32_t rank) const;
  Individual* indiv_by_position(int16_t x, int16_t y) const;

  // Accessors to output manager stuff
  int64_t	backup_step() const { return output_m()->backup_step(); }
  int64_t	end_step() const { return t_end_; }
  bool record_tree() const { return output_m()->record_tree(); }
  bool record_light_tree() const { return output_m()->record_light_tree(); }
  int32_t tree_step() const { return static_cast<int32_t>(output_m()->tree_step()); }
  Tree* tree() const { return output_m()->tree(); }
  LightTree* light_tree() const { return output_m()->light_tree(); }
  bool check_simd() const { return check_simd_; }

  // =======================================================================
  //                          Accessors: setters
  // =======================================================================
  void set_world(World* new_world) { world_ = new_world; }
  void set_t_end(int64_t t_end) { t_end_ = t_end; }
  void set_HT_ins_rate(double HT_ins_rate) { exp_s_->set_HT_ins_rate(HT_ins_rate); }
  void set_HT_repl_rate(double HT_repl_rate) { exp_s_->set_HT_repl_rate(HT_repl_rate); }
  void set_with_mrca(bool w_mrca) { with_mrca_ = w_mrca; }
  void set_anc_stat(bool w_anc_stat) { anc_stat_ = w_anc_stat; }
  #ifdef HAVE_MPI
  void set_nb_rank(int nb_rank) { exp_s()->set_nb_rank(nb_rank); }
  #endif
  // =======================================================================
  //                                 Operators
  // =======================================================================

  // =======================================================================
  //                               Public Methods
  // =======================================================================
  void InitializeWorld(int16_t grid_width,
                       int16_t grid_height,
                       std::shared_ptr<JumpingMT> prng,
                       std::shared_ptr<JumpingMT> mut_prng,
                       std::shared_ptr<JumpingMT> stoch_prng,
                       Habitat& habitat,
                       bool share_phenotypic_target);
  void Save(bool create = false) const;
  void WriteSetupFiles(bool first = false) const;
  void WriteDynamicFiles(bool first = false) const;
  void WriteDynamicFiles(int64_t gen, SaveWorld* world, bool create = false) const;
  void save_copy(char* dir, int64_t time = 0) const;
  void load(int64_t first_gener, bool verbose = false, bool to_be_run = true
  #ifdef HAVE_MPI
  , int32_t rank = -1
  #endif
  ) { load(".", first_gener, verbose, to_be_run
    #ifdef HAVE_MPI
    ,rank
    #endif
    ); 
  }

  void load(const char* dir, int64_t t0,
            bool verbose = false, bool to_be_run = true
            #ifdef HAVE_MPI
            , int32_t rank = -1
            #endif
            );
  /// Load an experiment with default files from the current directory
  void load(int64_t t0, char* exp_s_file, char* exp_backup_file, char* sp_struct_file, char* out_p_file,
            bool verbose = false, bool to_be_run = true
            #ifdef HAVE_MPI
            , int32_t rank = -1
            #endif
            );
  void run_evolution();
  virtual void display() {};
  void update_best();
  void FillGridWithClones(Individual& dolly) { world_->FillGridWithClones(dolly); }
  void update(Observable& o, ObservableEvent e, void* arg) override {}

  // DNA Mutator tables
    DnaMutator** dna_mutator_array_ = nullptr;

    // SIMD Stuff
    ExpManager_7* exp_m_7_ = nullptr;

  /// Which individual will be placed in each GridCell
  int32_t* next_generation_reproducer_ = nullptr;

    int grain_size = 1;

      #ifdef HAVE_MPI
  int32_t rank_;
  int32_t* prng_seed_;
  int32_t* mut_prng_seed_;
  int32_t* stoch_prng_seed_;
  #endif
  
        double w_max_ = 0;
 protected:
  // =======================================================================
  //                              Protected Methods
  // =======================================================================
  void step_to_next_generation();

  void load(gzFile& exp_s_file,
            gzFile& exp_backup_file,
            gzFile& sp_struct_file,
            gzFile& out_p_file,
            bool verbose = false,
            bool to_be_run = true
            #ifdef HAVE_MPI
            , int32_t rank = -1
            #endif
            );
  void create_missing_directories(const char* dir = ".") const;
  void open_backup_files(gzFile& sel_file,
                         gzFile& sp_struct_file,
                         int64_t t,
                         const char mode[3],
                         const char* dir = "."
                         #ifdef HAVE_MPI
                         , int rank = -1
                         #endif
                         ) const;
  void close_backup_files(gzFile& sel_file,
                          gzFile& sp_struct_file) const;
  void open_setup_files(gzFile& exp_s_gzfile,
                        gzFile& out_p_gzfile,
                        int64_t t,
                        const char mode[3],
                        const char* dir = "."
                        #ifdef HAVE_MPI
                        , int32_t rank = -1
                        #endif
                        ) const;
  void close_setup_files(gzFile& exp_s_gzfile, gzFile& out_p_gzfile) const;

  // =======================================================================
  //                             Protected Attributes
  // =======================================================================
  /// Experimental setup
  ExpSetup* exp_s_;

  /// Spatial structure
  World* world_;

  /// Output manager
  OutputManager*output_m_;
  /// Time step up to which we want to simulate
  int64_t t_end_;

  /// Should the simulation be stopped? Set to true when ctrl-Q is received.
  /// Will cause the simulation to be ended after the current time step is completed.
  bool quit_signal_received_;

    bool first_gen = true;

    bool check_simd_ = false;

        //wether you stop the simulation based on the age of mrca or not
        bool with_mrca_ = false;

        //do you record the stat of the lineage ?
        bool anc_stat_ = false;






};

// ===========================================================================
//                             Getters' definitions
// ===========================================================================

// ===========================================================================
//                             Setters' definitions
// ===========================================================================

// ===========================================================================
//                            Operators' definitions
// ===========================================================================

// ===========================================================================
//                         Inline methods' definition
// ===========================================================================

} // namespace aevol
#endif // AEVOL_EXP_MANAGER_H_
