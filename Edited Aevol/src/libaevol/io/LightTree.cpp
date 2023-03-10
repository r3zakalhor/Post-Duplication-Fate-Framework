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

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <fstream>

// =================================================================
//                            Project Files
// =================================================================
#include "7/Dna_7.h"
#include "7/Individual_7.h"
#include "7/ExpManager_7.h"
#include "AeTime.h"
#include "AncestorStats.h"
#include "LightTree.h"

namespace aevol {


//##############################################################################
//                                                                             #
//                                Class LightTree                              #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================

LightTree::LightTree(ExpManager* exp_m, int64_t tree_step) {
  anc_stat_          = nullptr;
  saved_indivs_time_ = -1;
  exp_m_ = exp_m;
  tree_step_ = tree_step;
}

// =================================================================
//                             Destructors
// =================================================================

LightTree::~LightTree() noexcept {
  for(auto node_time : allNodes_) {
    for(auto node : node_time.second) {
      allNodes_[node_time.first].erase(node.first);
      delete node.second;
    }
  }

  for(auto p_indiv : saved_indivs_) {
    delete p_indiv.second;
  }

  for(auto p_indiv : saved_simd_indivs_) {
    delete p_indiv.second;
  }

  if(anc_stat_ != nullptr) {
    delete anc_stat_;
  }
}

// =================================================================
//                            Public Methods
// =================================================================


  void LightTree::init_tree(int64_t time, std::list<Individual*> root_indiv) {
  if(time > 0) {
    read_from_tree_file();
  }
  else {
    mrca_time_ = 0;
    saved_mrca_time_ = 0;
    for(auto indiv : root_indiv) {
      Node* leaf = new Node(AeTime::time(), indiv->id());
      allNodes_[leaf->tid_][leaf->id_] = leaf;
    }
  }
}

void LightTree::update_tree(int64_t gen,
                            Individual_7** prev_internal_simd_struct){
  parentsNodes_.clear();
  // debug std::cout << "New generation : " << gen << '\n';
  std::unordered_map<int32_t, Node*> previous = allNodes_[gen-1];

  link_nodes(gen);
  //parentsNodes_ is filled

  //check for dead-end in the tree and cut the branch
  for(auto oldNode : previous) {
    if(std::find(parentsNodes_.begin(), parentsNodes_.end(), oldNode.first) == parentsNodes_.end()) {
      // debug std::cout << "Death of : " << oldNode.second->id << '\n';
      prune(oldNode.second);
    }
  }

  int64_t prev_mrca = mrca_time_;
  update_mrca(gen);
  // debug std::cout << "Generation of the MRCA : " << mrca_time_ << '\n';
  if(anc_stat_ != nullptr) {
    if(prev_mrca != mrca_time_) {
      if(prev_mrca == 0) {
          Individual* to_set;
          if (ExpManager_7::standalone_simd) {
            Individual_7* simd_indiv = prev_internal_simd_struct[allNodes_[0].begin()->second->id_];
              int x = simd_indiv->indiv_id / exp_m_->world()->height();
              int y = simd_indiv->indiv_id % exp_m_->world()->height();

            to_set = new Individual(exp_m_,
                                    exp_m_->world()->grid(x,y)->mut_prng(),
                                    exp_m_->world()->grid(x,y)->stoch_prng(),
                                    exp_m_->exp_s()->mut_params(),
                                    exp_m_->best_indiv()->w_max(),
                                    exp_m_->exp_s()->min_genome_length(),
                                    exp_m_->exp_s()->max_genome_length(),
                                    false,
                                    simd_indiv->indiv_id,
                                    "",
                                    0);
              int32_t nb_blocks_ = simd_indiv->dna_->length()/BLOCK_SIZE + 1;
              char* dna_string = new char[nb_blocks_ * BLOCK_SIZE];
              memset(dna_string,0,
                     (simd_indiv->dna_->length()+1) * sizeof(char));


              char* to_copy = simd_indiv->dna_->to_char();


              //printf("Copy DNA for indiv %d size %d (%d x %d)\n",indiv_id,previous_individuals[indiv_id]->dna_->length(),nb_blocks_,BLOCK_SIZE);
              memcpy(dna_string, to_copy,
                     (simd_indiv->dna_->length()+1) * sizeof(char));


              to_set->genetic_unit_list_.clear();
              to_set->add_GU(dna_string, simd_indiv->dna_->length());
              to_set->genetic_unit_nonconst(0).set_min_gu_length(exp_m_->exp_s()->min_genome_length());
              to_set->genetic_unit_nonconst(0).set_max_gu_length(exp_m_->exp_s()->max_genome_length());
              to_set->EvaluateInContext(exp_m_->world()->grid(x,y)->habitat());
              to_set->compute_statistical_data();
          } else
            to_set = exp_m_->indiv_by_id(allNodes_[0].begin()->second->id_);

          anc_stat_->setup_anc_indiv(to_set);
      }
      for(int64_t t = prev_mrca+1 ; t <= mrca_time_ ; t++) {
        anc_stat_->process_evolution(allNodes_[t].begin()->second->replics_, t);
      }
      anc_stat_->Flush();
    }
  }
}

void LightTree::write_tree(int64_t t /* = -1 */) {
  if(t == -1) {
    t = AeTime::time();
  }

  if(t == mrca_time_) {
    std::unordered_map<int32_t, Node*> Nodes = allNodes_[t];

    std::ofstream light_tree_file;
    light_tree_file.open (LIGHTTREE_FILE_NAME);
    if(!light_tree_file) {
      std::cerr << "ERROR : could not open the file" << std::endl;
      exit(EXIT_FAILURE);
    }

    light_tree_file << mrca_time_<< std::endl;

    for(auto root_node : Nodes) {
      light_tree_file << root_node.second->nhx_ << ":" << mrca_time_ << ";" << std::endl;
    }
    light_tree_file.close();
  }
  else {
    std::unordered_map<int32_t, Node*> Nodes = allNodes_[t-1];
    for(auto pair_node : Nodes) {
      std::string local = "(";
      bool first = true;
      for (auto child_time : pair_node.second->childs_) {
        for (auto child : child_time.second){
          if (!first)
            local=local+",";
          else
            first = false;

          if (child.second->nhx_ == "")
            local=local+std::to_string(child.second->id_);
          else {
            local = local + child.second->nhx_;
          }
        }
      }
      local=local+")";
      pair_node.second->nhx_ = local;
    }
    write_tree(t-1);
  }
}

void LightTree::write_to_tree_file(int64_t gen, gzFile trunc_file, gzFile branches_file) {
  int64_t saved_time = saved_mrca_time_;
  //in case MRCA dose not change between two saves
  if(saved_time < mrca_time_) {
    for( ; saved_time < mrca_time_ ; saved_time++) {
      if(saved_time > 0) {
        for(auto node : allNodes_[saved_time]) {
          (node.second->replics_)->write_to_tree_file(trunc_file);
        }
      }
    }
    //after saving the trunc, you erase it from the RAM
    // debug std::cout << "deleting report from " << saved_mrca_time_ << " to " << saved_time-1 << '\n';
    ((allNodes_[saved_time].begin())->second)->parent_ = nullptr;
    prune((allNodes_[saved_time-1].begin())->second);
  }

  for( ; saved_time <= gen ; saved_time++) {
    if(saved_time > 0) {
      for(auto node : allNodes_[saved_time]) {
        (node.second->replics_)->write_to_tree_file(branches_file);
        int64_t generation = node.second->tid_;
        gzwrite(branches_file, &generation, sizeof(generation));
      }
    }
  }
  saved_mrca_time_ = mrca_time_;
}

void LightTree::read_from_tree_file() {
  std::ifstream light_tree_file;
  light_tree_file.open(LIGHTTREE_FILE_NAME);
  if(light_tree_file and light_tree_file.peek() != EOF) {
    light_tree_file >> mrca_time_;
    saved_mrca_time_ = mrca_time_;
  }
  else {
    exit(EXIT_FAILURE);
  }

  char branches_file_name[50];
  sprintf(branches_file_name, "lightTree/tree_branches.ae");
  gzFile branches_file = gzopen( branches_file_name, "r" );

  int c = gzgetc(branches_file);

  while(c > -1) {
    gzungetc(c, branches_file);
    ReplicationReport* rep = new ReplicationReport(branches_file, nullptr);
    int64_t generation;
    gzread(branches_file, &generation, sizeof(generation));
    Node* new_node = new Node(generation, rep->id());
    new_node->replics_ = rep;
    allNodes_[new_node->tid_][new_node->id_] = new_node;

    c = gzgetc(branches_file);
  }

  gzclose(branches_file);

  for(int64_t i = mrca_time_+1 ; i <= AeTime::time() ; i++) {
    link_nodes(i);
  }
}
#ifdef __REGUL
void LightTree::keep_indivs(std::list<Individual_R*> indivs) {
#else
  void LightTree::keep_indivs(std::list<Individual*> indivs) {
#endif
  #ifdef _OPENMP
  #pragma omp taskgroup
  {
  #endif
#ifdef __REGUL
  for(Individual_R* indiv : indivs)
#else
    for(Individual* indiv : indivs)
#endif
    #ifdef _OPENMP
    #pragma omp task firstprivate(indiv)
    #endif
#ifdef __REGUL
    saved_indivs_[indiv->id()] = new Individual_R(*indiv);
#else
  saved_indivs_[indiv->id()] = new Individual(*indiv);
#endif
  #ifdef _OPENMP
  }
  #endif


  saved_indivs_time_ = AeTime::time();//(*indivs.begin())->age();
}

void LightTree::keep_indivs(std::list<Individual_7*> indivs,
                            DnaFactory* dna_factory, FuzzyFactory_7* fuzzy_factory) {
#ifdef _OPENMP
      #pragma omp taskgroup
  {
#endif
      for(Individual_7* indiv : indivs)
#ifdef _OPENMP
#pragma omp task firstprivate(indiv)
#endif
        saved_simd_indivs_[indiv->indiv_id] = new Individual_7(exp_m_,indiv,dna_factory,fuzzy_factory);
#ifdef _OPENMP
      }
#endif
      saved_indivs_time_ = AeTime::time();//(*indivs.begin())->age();
}

void LightTree::save_mrca_indiv() {
    printf("  Printing the mrca individual into " "mrca_dna.txt" "\n");
    FILE* org_file = fopen("mrca_dna.txt", "w");
    gzFile indiv_file = gzopen("mrca_indiv.ae", "w");

    if (ExpManager_7::standalone_simd) {
      Individual_7*mrca = saved_simd_indivs_.begin()->second;
      fputs(mrca->dna_->to_char(), org_file);

      int x=0,y=0;


      Individual * indiv = new Individual(exp_m_,
                                          exp_m_->world()->grid(x,y)->mut_prng(),
                                          exp_m_->world()->grid(x,y)->stoch_prng(),
                                          exp_m_->exp_s()->mut_params(),
                                          exp_m_->best_indiv()->w_max(),
                                          exp_m_->exp_s()->min_genome_length(),
                                          exp_m_->exp_s()->max_genome_length(),
                                          false,
                                          mrca->indiv_id,
                                          "",
                                          0);
      int32_t nb_blocks_ = mrca->dna_->length()/BLOCK_SIZE + 1;
      char* dna_string = new char[nb_blocks_ * BLOCK_SIZE];
      memset(dna_string,0,
             (mrca->dna_->length()+1) * sizeof(char));


      char* to_copy = mrca->dna_->to_char();


      //printf("Copy DNA for indiv %d size %d (%d x %d)\n",indiv_id,previous_individuals[indiv_id]->dna_->length(),nb_blocks_,BLOCK_SIZE);
      memcpy(dna_string, to_copy,
             (mrca->dna_->length()+1) * sizeof(char));


      indiv->genetic_unit_list_.clear();
      indiv->add_GU(dna_string, mrca->dna_->length());
      indiv->genetic_unit_nonconst(0).set_min_gu_length(exp_m_->exp_s()->min_genome_length());
      indiv->genetic_unit_nonconst(0).set_max_gu_length(exp_m_->exp_s()->max_genome_length());
      indiv->EvaluateInContext(exp_m_->world()->grid(x,y)->habitat());
      indiv->compute_statistical_data();

      indiv->save(indiv_file);
  } else {
    Individual *mrca = saved_indivs_.begin()->second;
      fputs(mrca->genetic_unit_sequence(0), org_file);
        mrca->save(indiv_file);
  }
  fclose(org_file);

  gzclose(indiv_file);
}

void LightTree::signal_end_of_generation() {
  for(auto leaf : allNodes_[AeTime::time()]) {
    leaf.second->replics_->signal_end_of_generation();
  }
}

void LightTree::update(Observable& o, ObservableEvent e, void* arg) {
  switch (e) {
    case NEW_INDIV : {
      // Initialize the replication report corresponding to the new individual
      auto ievent = reinterpret_cast<NewIndivEvent*>(arg);

      // debug std::cout << "new indiv id : " << new_indiv->id() << " from parent : " << parent->id() << '\n';
      ReplicationReport* new_rep = new ReplicationReport();
      if (ExpManager_7::standalone_simd) {
        //printf("Update with %d %d\n",ievent->indiv_id_,ievent->parent_id_);
        new_rep->
                init(this, ievent->simd_child, ievent->simd_parent,ievent->indiv_id_,ievent->parent_id_);
      } else {
        new_rep->
                init(this, ievent->child, ievent->parent,ievent->indiv_id_,ievent->parent_id_);
      }

      Node *new_node  = new Node(AeTime::time(), ievent->indiv_id_);
      new_node->replics_ = new_rep;
      #ifdef __OPENMP_TASK
      #pragma omp critical(lighttree)
      {
      #endif
      allNodes_[AeTime::time()][ievent->indiv_id_] = new_node;
      #ifdef __OPENMP_TASK
      }
      #endif
      break;
    }
    case END_GENERATION :
      signal_end_of_generation();
      break;
    case END_REPLICATION :
      break;
    default :
      Utils::ExitWithDevMsg("Event not handled", __FILE__, __LINE__);
  }
}

void LightTree::setup_anc_stat() {
  anc_stat_ = new AncestorStats();
  anc_stat_->setup_anc_stat(mrca_time_);
}

void LightTree::close_anc_stat() {
  anc_stat_->Close();
}

// =================================================================
//                  Non-inline accessors' definition
// =================================================================

// =================================================================
//                           Protected Methods
// =================================================================

void LightTree::prune(Node * dead_node) {
  Node * father = dead_node->parent_;
  int64_t tid_dead = dead_node->tid_;
  int32_t id_dead = dead_node->id_;
  delete dead_node;
  if (ExpManager_7::standalone_simd) {
    if (tid_dead == saved_indivs_time_ && saved_simd_indivs_.size() > 1) {
      delete saved_simd_indivs_[id_dead];
      saved_simd_indivs_.erase(id_dead);
    }
  } else {
    if (tid_dead == saved_indivs_time_ && saved_indivs_.size() > 1) {
      delete saved_indivs_[id_dead];
      saved_indivs_.erase(id_dead);
    }
  }
  // debug std::cout << "Killing : " << id_dead << " from time : " << tid_dead << '\n';
  #ifdef __OPENMP_TASK
  #pragma omp critical(lighttree)
  {
  #endif
  allNodes_[tid_dead].erase(id_dead);
  if(allNodes_[tid_dead].empty()) {
    allNodes_.erase(tid_dead);
  }
  #ifdef __OPENMP_TASK
  }
  #endif
  if (father != nullptr) {

    (father->childs_)[tid_dead].erase(id_dead);

    int size = 0;
    for(auto child_time : father->childs_) {
      size += (child_time.second).size();
    }
    if (size < 1) {
      prune(father);
    }
  }
}

void LightTree::update_mrca(int64_t gen) {
  for(int64_t i = mrca_time_ ; i < gen ; i++) {
    if(allNodes_[i].size() > 1) {
      mrca_time_ = std::max(i-1, (int64_t)0);
      break;
    }
  }
}

  void LightTree::link_nodes(int64_t t) {
  std::unordered_map<int32_t, Node*> previous = allNodes_[t-1];
  for(auto leaf_pair : allNodes_[t]) {

    Node* leaf = leaf_pair.second;
    auto foundFather = previous.find(leaf->replics_->parent_id());

    //we schould always find a father (no orphan, it is crual)
    if (foundFather != previous.end()) {
      Node* father = foundFather->second;
      leaf->parent_ = father;
      father->childs_[leaf->tid_][leaf->id_] = leaf;

      if(std::find(parentsNodes_.begin(), parentsNodes_.end(), father->id_) == parentsNodes_.end()) {
        parentsNodes_.push_back(father->id_);
      }
    }
  }
}

} // namespace aevol
