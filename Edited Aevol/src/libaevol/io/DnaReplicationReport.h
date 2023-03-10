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

#ifndef AEVOL_DNA_REPLICATION_REPORT_H_
#define AEVOL_DNA_REPLICATION_REPORT_H_

#include <list>
#include <memory>

#include "Mutation.h"
#include "LocalMutation.h"
#include "HorizontalTransfer.h"
#include "Rearrangement.h"

namespace aevol {

class Dna;

class DnaReplicationReport {
  friend class Dna;

 public :

  // =================================================================
  //                             Constructors
  // =================================================================
  DnaReplicationReport() = default;
  DnaReplicationReport(const DnaReplicationReport&);
  DnaReplicationReport(DnaReplicationReport&&) = delete;

  // =================================================================
  //                             Destructor
  // =================================================================
  ~DnaReplicationReport() = default;

  // ==========================================================================
  //                                Operators
  // ==========================================================================
  /// Copy assignment
  DnaReplicationReport& operator=(const DnaReplicationReport& other) = delete;
  /// Move assignment
  DnaReplicationReport& operator=(DnaReplicationReport&& other) = delete;

  // Accessors
  const std::list<std::unique_ptr<const LocalMutation>>& mutations() const {
    return mutations_;
  };
  const std::list<std::unique_ptr<const Rearrangement>>& rearrangements() const {
    return rearrangements_;
  };
  const std::list<std::unique_ptr<const HorizontalTransfer>>& HT() const {
    return ht_;
  };
  int32_t nb(MutationType t) const;

  // Public Methods
  void compute_stats();  // useful when we inspect a tree file
  void add_mut(Mutation* mut);
  void add_local_mut(Mutation* mut);
  void add_rear(Mutation* mut);
  void add_HT(Mutation* mut);

  void clear();

  void write_to_tree_file(gzFile tree_file) const;
  void read_from_tree_file(gzFile tree_file);

  template <typename F> void iter_muts(const F& f) const {
    for (const auto& mut: ht_) {
      f(mut);
    }
    for (const auto& mut: rearrangements_) {
      f(mut);
    }
    for (const auto& mut: mutations_) {
      f(mut);
    }
  }

 protected :
  /// Lists of mutations, rearrangements and undergone
  std::list<std::unique_ptr<const LocalMutation>> mutations_;
  std::list<std::unique_ptr<const Rearrangement>> rearrangements_;
  std::list<std::unique_ptr<const HorizontalTransfer>> ht_;
  // Number of mutations/rearrangements/HT of each (simple) type undergone
  int32_t nb_mut_[10] = {0,0,0,0,0,0,0,0,0,0};
};

} // namespace aevol
#endif // AEVOL_DNA_REPLICATION_REPORT_H_
