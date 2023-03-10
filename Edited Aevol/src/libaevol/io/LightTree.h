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


#ifndef AEVOL_LIGHTTREE_H_
#define AEVOL_LIGHTTREE_H_


// =================================================================
//                              Includes
// =================================================================
#include "DnaFactory.h"
#include "ObservableEvent.h"
#include "Observer.h"
#include "ReplicationReport.h"
#include "FuzzyFactory_7.h"
#include <cstdio>
#include <cstdlib>
#include <inttypes.h>
#include <unordered_map>
#include <vector>

namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================
class AncestorStats;
class Individual;
class Individual_7;
class ExpManager;

class Node {
public:
  Node(int64_t t, int32_t lid) { tid_ = t; id_ = lid; };

  ~Node() { delete replics_; };
  ReplicationReport* replics() const { return replics_; };

  int64_t tid_;
  int32_t id_;
  std::unordered_map<int64_t, std::unordered_map<int32_t, Node*>> childs_;
  Node* parent_ = nullptr;
  ReplicationReport * replics_ = nullptr;
  std::string nhx_ = "";
};

class LightTree : public Observer
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    LightTree(ExpManager* exp_m, int64_t tree_step);

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~LightTree() noexcept;

    // =================================================================
    //                        Accessors: getters
    // =================================================================

    int64_t mrca_time() const { return mrca_time_; };
    int64_t saved_mrca_time() const { return saved_mrca_time_; };
    int64_t saved_indivs_time() const { return saved_indivs_time_; };
    ReplicationReport* get_first_replics(int64_t t) const { return allNodes_.at(t).begin()->second->replics(); };

    // =================================================================
    //                        Accessors: setters
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================

    void init_tree(int64_t time, std::list<Individual*> root_indiv);
    // for the generation gen, the methode link the nodes with the generation gen-1
    // by creating the link parent/children
    // it also prune the tree
    // if ask it perform ancestor_stat
    void update_tree(int64_t gen,
                     Individual_7** prev_internal_simd_struct = nullptr);

    // write the Newick format tree
    void write_tree(int64_t t = -1);

    void write_to_tree_file(int64_t gen, gzFile trunc_file, gzFile branches_file);

    void read_from_tree_file();

#ifdef __REGUL
    void keep_indivs(std::list<Individual_R*> indivs);
#else
  void keep_indivs(std::list<Individual*> indivs);
#endif
    void keep_indivs(std::list<Individual_7*> indivs, DnaFactory* dna_factory, FuzzyFactory_7* fuzzy_factory);

    void save_mrca_indiv();

    void setup_anc_stat();

    void close_anc_stat();

    void signal_end_of_generation();

  void update(Observable& o, ObservableEvent e, void* arg) override;

    inline int64_t  tree_step() const {
        return tree_step_;
    };

    // =================================================================
    //                           Public Attributes
    // =================================================================

  protected :
    // =================================================================
    //                           Protected Methods
    // =================================================================

    void prune(Node * oldNode);
    void update_mrca(int64_t gen);

    // attribute the relation child/parent for the nodes of the t generation
    // to the t-1 generation
    // it fills the parentsNodes_
    void link_nodes(int64_t t);

    // =================================================================
    //                          Protected Attributes
    // =================================================================

    // the actual tree
    std::unordered_map<int64_t, std::unordered_map<int32_t, Node*>> allNodes_;
    // all the nodes that have children in the current generation
    // if a node is not in this list, we prune its branche
    std::vector<int32_t> parentsNodes_;
    // generation of the actual mrca
    int64_t mrca_time_;
    // generation of the last saved mrca
    int64_t saved_mrca_time_;

    // list of all the Individual at the generation wented to stop the simulation
    std::unordered_map<int32_t, Individual*> saved_indivs_;
    std::unordered_map<int32_t, Individual_7*> saved_simd_indivs_;
    int64_t saved_indivs_time_;

    AncestorStats* anc_stat_;

    ExpManager* exp_m_;


    int64_t tree_step_;
};


// =====================================================================
//                           Getters' definitions
// =====================================================================

// =====================================================================
//                           Setters' definitions
// =====================================================================


// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol

#endif // AEVOL_TREE_H_
