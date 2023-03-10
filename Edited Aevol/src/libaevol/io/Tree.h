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


#ifndef AEVOL_TREE_H_
#define AEVOL_TREE_H_


// =================================================================
//                              Includes
// =================================================================
#include <inttypes.h>

#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <map>

#include "ReplicationReport.h"
#include "Observer.h"
#include "ObservableEvent.h"
#include "ae_enums.h"


namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================
class ExpManager;

class RepTree
{
public:
    RepTree() { a = 1; }
    ~RepTree() = default;

    int a = -1;
};

class Tree : public Observer
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    Tree() = delete;
    Tree(const Tree &model) = delete;
    // To be used when we want to run a simulation.
    Tree(ExpManager* exp_m, int64_t tree_step);
    // To be used when we want to INSPECT a tree,
    // not when we want to run a simulation.
    Tree(ExpManager* exp_m, char* tree_file_name);

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~Tree() noexcept;

    // =================================================================
    //                        Accessors: getters
    // =================================================================
    inline int64_t  tree_step() const {
      return tree_step_;
    };

    // Precondition for the following methods:
    // the tree was emptied every TREE_STEP generations ==> it contains
    // only the last generations since the last emptying ==> do not ask
    // something about an older generation
    ReplicationReport** reports(int64_t t) const;
    ReplicationReport* report_by_index(int64_t t, int32_t index) const;
    ReplicationReport* report_by_rank(int64_t t, int32_t rank) const;


    // =================================================================
    //                        Accessors: setters
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    void signal_end_of_generation();
    void write_to_tree_file(char* tree_file_name);

  void update(Observable& o, ObservableEvent e, void* arg) override;

    void update_new_indiv(NewIndivEvent* evt);

    void update_end_replication(EndReplicationEvent* evt);

    void update_end_generation();

    // =================================================================
    //                           Public Attributes
    // =================================================================
    static const int32_t NO_PARENT;


    ExpManager* exp_m_;


  protected :
    // =================================================================
    //                           Protected Methods
    // =================================================================

    // =================================================================
    //                          Protected Attributes
    // =================================================================


    int64_t tree_step_;

    ReplicationReport*** replics_;
    // Two-dimensional table of ReplicationReport*
    //    dimension 1 (lines)   : generation
    //    dimension 2 (columns) : individual
    //
    // !!!!! WARNING !!!!!
    // The report at line l, column c is for the
    // replication that created the indiv with index c of generation l+1

    // light tree representation
    int32_t** parent_;
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
