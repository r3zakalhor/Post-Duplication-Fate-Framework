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



// =================================================================
//                            Project Files
// =================================================================
#include "Tree.h"
#include "7/ExpManager_7.h"

#include "7/Individual_7.h"
#include "ExpManager.h"
#include "ExpSetup.h"
#include "Individual.h"
#include "Utils.h"
#include "macros.h"

#include <algorithm>

namespace aevol {


//##############################################################################
//                                                                             #
//                                Class Tree                                #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================
const int32_t Tree::NO_PARENT = -1;

// =================================================================
//                             Constructors
// =================================================================
Tree::Tree(ExpManager* exp_m, int64_t tree_step) {
  exp_m_ = exp_m;
  tree_step_ = tree_step;

  replics_ = new ReplicationReport** [tree_step_];

  for (int32_t time = 0 ; time < tree_step ; time++) {
    replics_[time] = new ReplicationReport* [exp_m->nb_indivs()];
    for (int32_t num_indiv = 0 ;
         num_indiv < exp_m->nb_indivs() ;
         num_indiv++) {
      replics_[time][num_indiv] = new ReplicationReport();
    }
  }
}


/**
 *
 */
Tree::Tree(ExpManager* exp_m, char* tree_file_name) {
  exp_m_ = exp_m;
  tree_step_ = exp_m_->tree_step();

  gzFile tree_file = gzopen(tree_file_name, "r");
  if (tree_file == Z_NULL) {
    printf("ERROR : Could not read tree file %s\n", tree_file_name);
    exit(EXIT_FAILURE);
  }

  replics_ = new ReplicationReport** [tree_step_];

  #ifdef HAVE_MPI
  int32_t global_grid_width = exp_m->exp_s()->global_grid_width();
  int32_t global_grid_height = exp_m->exp_s()->global_grid_height();

  int32_t rank_x = exp_m->exp_s()->rank_width();
  int32_t rank_y = exp_m->exp_s()->rank_height();

  int32_t grid_width = exp_m->grid_width();
  int32_t grid_height = exp_m->grid_height();

  int32_t local_rank_x = exp_m->rank() / rank_y;
  int32_t local_rank_y = exp_m->rank() % rank_y;
  #endif
  //for (int64_t t = AeTime::time()-tree_step_+1 ; t <= AeTime::time() ; t++) {
  for (int64_t t = 0 ; t < tree_step_ ; t++) {
    replics_[t] = new ReplicationReport* [exp_m_->nb_indivs()];

    for (int32_t indiv_i = 0 ;
         indiv_i < exp_m_->nb_indivs() ;
         indiv_i++) {
      // Retrieve a replication report
      ReplicationReport* replic_report = new ReplicationReport(tree_file,
                                                               nullptr);

      // Put it at its rightful position
      //printf("Replication Report %d\n",replic_report->id());


        #ifdef HAVE_MPI
        int32_t local_x = exp_m_->exp_m_7_->globalXtoLocalX(replic_report->id() / global_grid_height);
        int32_t local_y = exp_m_->exp_m_7_->globalYtoLocalY(replic_report->id() % global_grid_height);
        int32_t local_id_ = local_x * grid_height + local_y;
        // printf("Replication Report %d (%d %d) Local ID %d (%d %d) R %d (%d %d) GH %d %d OFFS %d %d\n",replic_report->id(),
        //         (replic_report->id() / global_grid_height),
        //         (replic_report->id() % global_grid_height),
        //         local_id_,local_x,local_y,
        //         exp_m->rank(),local_rank_x,local_rank_y,global_grid_height,grid_height,
        //         local_rank_x * rank_x,
        //         local_rank_y * rank_y);
        // printf("%d -- Loading Replication Report %d : GID %d PID %d\n",t,local_id_,replic_report->id(),replic_report->parent_id());
        replics_[t][local_id_] = replic_report;
        #else
        replics_[t][replic_report->id()] = replic_report;
        #endif
    }
  }
  gzclose(tree_file);
}




// =================================================================
//                             Destructors
// =================================================================
Tree::~Tree() {
  if (replics_ != NULL)  {
    for (int32_t i = 0 ; i < tree_step_ ; i++)
      if (replics_[i] != NULL) {
        for (int32_t j = 0 ; j < exp_m_->nb_indivs() ; j++)
          delete replics_[i][j];
        delete [] replics_[i];
      }
    delete [] replics_;
  }
}

// =================================================================
//                            Public Methods
// =================================================================
    ReplicationReport** Tree::reports(int64_t t) const {
        return replics_[Utils::mod(t - 1, tree_step_)];
    }

    ReplicationReport* Tree::report_by_index(int64_t t, int32_t index) const {
        // printf("Loading from %d Indiv %d (%d) : %p\n",Utils::mod(t - 1, tree_step_),
        //         index,replics_[Utils::mod(t - 1, tree_step_)][index]->id(),replics_[Utils::mod(t - 1, tree_step_)][index]);
        return replics_[Utils::mod(t - 1, tree_step_)][index];
    }


    ReplicationReport* Tree::report_by_rank(int64_t t, int32_t rank) const {
         printf("error!!! Rank is not supported anymore\n");exit(-1);
      }


    void Tree::signal_end_of_generation() {
        auto cur_reports = reports(AeTime::time());
        for (int32_t i = 0; i < exp_m_->nb_indivs(); i++) {
            if (cur_reports[i] == nullptr){
                printf("error!!!\n");exit(-1);
            }

            cur_reports[i]->signal_end_of_generation();
        }
    }

    void Tree::write_to_tree_file(char* tree_file_name) {
      //printf("%d -- Saving Tree\n",exp_m_->rank());
        gzFile tree_file = gzopen( tree_file_name, "w" );
        // Write the tree in the backup
        for (int64_t t = 0 ; t < tree_step_ ; t++)
            for (int32_t indiv_i = 0 ; indiv_i < exp_m_->nb_indivs() ; indiv_i++) {
          //              if (exp_m_->rank() == 0) printf("%ld -- %d -- Saving to File -- Local ID %d GID %d-- Parent ID : %d\n",t,exp_m_->rank(),
          //     indiv_i,
          // replics_[t][indiv_i]->id_,
          // replics_[t][indiv_i]->parent_id_);
        // printf("%d -- %d -- Saving Replication Report %d : GID %d PID %d\n",t,exp_m_->rank(),indiv_i,replics_[t][indiv_i]->id(),replics_[t][indiv_i]->parent_id());

                assert(replics_[t][indiv_i] != NULL);
                replics_[t][indiv_i]->write_to_tree_file(tree_file);

            }

        gzclose(tree_file);


        // Reinitialize the tree
/*        for (int32_t time = 0 ; time < tree_step_ ; time++) {
            for (int32_t num_indiv = 0 ;
                 num_indiv < exp_m_->nb_indivs() ;
                 num_indiv++) {
                delete replics_[time][num_indiv];
                replics_[time][num_indiv] = nullptr;
            }
        }*/

    }

void Tree::update_new_indiv(NewIndivEvent* evt) {
  //  printf("%ld -- Tree Update -- Local ID %d (%d %d) -- Global ID : %d (Parent %d)\n",AeTime::time(),evt->x *
  //                                                                                                evt->simd_child->exp_m_->grid_height()
  //                                                                                                + evt->y,
  //         evt->x,
  //         evt->y,evt->indiv_id_,evt->parent_id_);

    if (ExpManager_7::standalone_simd) {
      // printf("%d -- ",Utils::mod(AeTime::time() - 1, tree_step_));
        replics_[Utils::mod(AeTime::time() - 1, tree_step_)][evt->x *
                                                             evt->simd_child->exp_m_->grid_height()
                                                             + evt->y]->
                init(this, evt->simd_child, evt->simd_parent, evt->indiv_id_, evt->parent_id_, evt->remote_);
    } else {
        replics_[Utils::mod(AeTime::time() - 1, tree_step_)][evt->x *
                                                             evt->child->exp_m()->grid_height()
                                                             + evt->y]->
                init(this, evt->child, evt->parent, evt->indiv_id_, evt->parent_id_);
    }

        //  if (exp_m_->rank() == 0) printf("%ld -- %d -- Local ID %d GID %d-- Parent ID : %d\n",AeTime::time(),exp_m_->rank(),
        //   evt->x *
        //                                                      evt->simd_child->exp_m_->grid_height()
        //                                                      + evt->y,
        //   replics_[Utils::mod(AeTime::time() - 1, tree_step_)][evt->x *
        //                                                      evt->simd_child->exp_m_->grid_height()
        //                                                      + evt->y]->id_,
        //   replics_[Utils::mod(AeTime::time() - 1, tree_step_)][evt->x *
        //                                                      evt->simd_child->exp_m_->grid_height()
        //                                                      + evt->y]->parent_id_);
}

void Tree::update_end_replication(EndReplicationEvent* evt) {
    if (ExpManager_7::standalone_simd) {
        replics_[Utils::mod(AeTime::time() - 1, tree_step_)][evt->x *

                                                         evt->simd_child->exp_m_->grid_height()
                                                         + evt->y]->signal_end_of_replication(evt->simd_child);
    } else {
        replics_[Utils::mod(AeTime::time() - 1, tree_step_)][evt->x *
                                                             evt->child->exp_m()->grid_height()
                                                             + evt->y]->signal_end_of_replication(evt->child);
    }

          //    if (exp_m_->rank() == 0) printf("%ld -- %d -- ER -- Local ID %d GID %d-- Parent ID : %d\n",AeTime::time(),exp_m_->rank(),
          // evt->x *
          //                                                    evt->simd_child->exp_m_->grid_height()
          //                                                    + evt->y,
          // replics_[Utils::mod(AeTime::time() - 1, tree_step_)][evt->x *
          //                                                    evt->simd_child->exp_m_->grid_height()
          //                                                    + evt->y]->id_,
          // replics_[Utils::mod(AeTime::time() - 1, tree_step_)][evt->x *
          //                                                    evt->simd_child->exp_m_->grid_height()
          //                                                    + evt->y]->parent_id_);
}


void Tree::update_end_generation() {
    signal_end_of_generation();
}

void Tree::update(Observable& o, ObservableEvent e, void* arg) {

  switch (e) {
    case NEW_INDIV : {
      // Initialize the replication report corresponding to the new individual

      auto ievent = reinterpret_cast<NewIndivEvent*>(arg);
       if (ExpManager_7::standalone_simd) {
           replics_[Utils::mod(AeTime::time() - 1, tree_step_)][ievent->x *
                                              ievent->simd_child->exp_m_->grid_height()
                                              + ievent->y]->
                      init(this, ievent->simd_child, ievent->simd_parent, ievent->indiv_id_, ievent->parent_id_);
          } else {
           replics_[Utils::mod(AeTime::time() - 1, tree_step_)][ievent->x *
                                              ievent->child->exp_m()->grid_height()
                                              + ievent->y]->
                      init(this, ievent->child, ievent->parent, ievent->indiv_id_, ievent->parent_id_);
          }

      break;
    }
    case END_GENERATION : {
      signal_end_of_generation();
      break;
    }
    case END_REPLICATION : {
      auto ievent = reinterpret_cast<EndReplicationEvent*>(arg);

            if (ExpManager_7::standalone_simd) {
                //printf("EoR %d : %p -- %p\n",ievent->simd_child->indiv_id,ievent->simd_child, replics_[AeTime::time()][ievent->simd_child->indiv_id]);
                replics_[Utils::mod(AeTime::time() - 1, tree_step_)][ievent->x *
                                                ievent->simd_child->exp_m_->grid_height()
                                                + ievent->y]->signal_end_of_replication(
                        ievent->simd_child);
            } else {
                replics_[Utils::mod(AeTime::time() - 1, tree_step_)][ievent->x *
                                                ievent->child->exp_m()->grid_height()
                                                + ievent->y]->signal_end_of_replication(
                        ievent->child);
            }

      break;
    }
    default : {
      Utils::ExitWithDevMsg("Event not handled", __FILE__, __LINE__);
    }
  }
}



// =================================================================
//                  Non-inline accessors' definition
// =================================================================

// =================================================================
//                           Protected Methods
// =================================================================
} // namespace aevol
