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

#include <cinttypes>
#include <list>
#include <cstdlib>
#include <cassert>

#if __cplusplus == 201103L
#include "make_unique.h"
#endif

#include "DnaReplicationReport.h"
#include "InsertionHT.h"
#include "ReplacementHT.h"

#include "Mutation.h"
#include "Duplication.h"
#include "Translocation.h"
#include "Inversion.h"
#include "Deletion.h"
#include "SmallDeletion.h"
#include "PointMutation.h"
#include "SmallInsertion.h"

namespace aevol {

void DnaReplicationReport::clear() {
/*  for (auto it = ht_.begin(); it < ht_.end(); it++) {
    delete (*it);
  }*/
        ht_.clear();
        rearrangements_.clear();
        mutations_.clear();

    nb_mut_[SWITCH] = 0;
    nb_mut_[S_INS]  = 0;
    nb_mut_[S_DEL]  = 0;
    nb_mut_[DUPL]   = 0;
    nb_mut_[DEL]    = 0;
    nb_mut_[TRANS]  = 0;
    nb_mut_[INV]    = 0;
    nb_mut_[INS_HT] = 0;
    nb_mut_[REPL_HT]= 0;
    }

DnaReplicationReport::DnaReplicationReport(const DnaReplicationReport& other) {
  Mutation* mut = nullptr;

  for (auto& ht : other.ht_) {
      mut = ht->Clone();
      add_HT(mut);
      delete mut;
  }

  for (auto& rear : other.rearrangements_) {
      mut = rear->Clone();
      add_rear(mut);
      delete mut;
  }

  for (auto& pmut : other.mutations_) {
      mut = pmut->Clone();
      add_local_mut(mut);
      delete mut;
  }
}

int32_t DnaReplicationReport::nb(MutationType t)  const {
  switch (t) {
    case S_MUT:
      assert(mutations_.size() ==
             static_cast<size_t>(nb_mut_[SWITCH] +
                                 nb_mut_[S_INS] +
                                 nb_mut_[S_DEL]));
      return mutations_.size();
    case REARR:
      assert(rearrangements_.size() ==
             static_cast<size_t>(nb_mut_[DUPL] +
                                 nb_mut_[DEL] +
                                 nb_mut_[TRANS] +
                                 nb_mut_[INV]));
      return rearrangements_.size();
    case H_T:
      assert(ht_.size() ==
             static_cast<size_t>(nb_mut_[INS_HT] +
                                 nb_mut_[REPL_HT]));
      return ht_.size();
    case INDEL:
      return nb_mut_[S_INS] + nb_mut_[S_DEL];
    default: // Simple mutation type.
      return nb_mut_[t];
  };
}

void DnaReplicationReport::add_mut(Mutation* mut) {
  if (mut->is_local_mut()) {
    add_local_mut(mut);
  }
  else if (mut->is_rear()) {
    add_rear(mut);
  }
  else if (mut->is_ht()) {
    add_HT(mut);
  }
}

void DnaReplicationReport::add_local_mut(Mutation* mut) {
  assert(mut->is_local_mut());
  std::unique_ptr<const LocalMutation> cmut = nullptr;
  switch(mut->mut_type()) {
    case SWITCH:
#if __cplusplus == 201103L
      cmut = make_unique<const PointMutation>(static_cast<PointMutation&>(*mut));
#else
      cmut = std::make_unique<const PointMutation>(static_cast<PointMutation&>(*mut));
#endif
      break;
    case S_DEL:
#if __cplusplus == 201103L
      cmut = make_unique<const SmallDeletion>(static_cast<SmallDeletion&>(*mut));
#else
      cmut = std::make_unique<const SmallDeletion>(static_cast<SmallDeletion&>(*mut));
#endif
      break;
    case S_INS:
#if __cplusplus == 201103L
      cmut = make_unique<const SmallInsertion>(static_cast<SmallInsertion&>(*mut));
#else
      cmut = std::make_unique<const SmallInsertion>(static_cast<SmallInsertion&>(*mut));
#endif
      break;
    default:
#if __cplusplus == 201103L
      cmut = make_unique<const SmallInsertion>(static_cast<SmallInsertion&>(*mut));
#else
      cmut = std::make_unique<const SmallInsertion>(static_cast<SmallInsertion&>(*mut));
#endif
      break;
  }
  mutations_.push_back(std::move(cmut));
  nb_mut_[mut->mut_type()]++;
}

void DnaReplicationReport::add_rear(Mutation* mut) {
  assert(mut->is_rear());

  std::unique_ptr<const Rearrangement> cmut = nullptr;
  switch(mut->mut_type()) {
    case DUPL:
#if __cplusplus == 201103L
      cmut = make_unique<const Duplication>(static_cast<Duplication&>(*mut));
#else
      cmut = std::make_unique<const Duplication>(static_cast<Duplication&>(*mut));
#endif
      break;
    case DEL:
#if __cplusplus == 201103L
      cmut = make_unique<const Deletion>(static_cast<Deletion&>(*mut));
#else
      cmut = std::make_unique<const Deletion>(static_cast<Deletion&>(*mut));
#endif
      break;
    case TRANS:
#if __cplusplus == 201103L
      cmut = make_unique<const Translocation>(static_cast<Translocation&>(*mut));
#else
      cmut = std::make_unique<const Translocation>(static_cast<Translocation&>(*mut));
#endif
      break;
    case INV:
#if __cplusplus == 201103L
      cmut = make_unique<const Inversion>(static_cast<Inversion&>(*mut));
#else
      cmut = std::make_unique<const Inversion>(static_cast<Inversion&>(*mut));
#endif
      break;
    default:
#if __cplusplus == 201103L
      cmut = make_unique<const Inversion>(static_cast<Inversion&>(*mut));
#else
      cmut = std::make_unique<const Inversion>(static_cast<Inversion&>(*mut));
#endif
      break;
  }
  rearrangements_.push_back(std::move(cmut));
  nb_mut_[mut->mut_type()]++;
}

void DnaReplicationReport::add_HT(Mutation* mut) {
  assert(mut->is_ht());

  std::unique_ptr<const HorizontalTransfer> cmut = nullptr;
  switch(mut->mut_type()) {
    case INS_HT:
#if __cplusplus == 201103L
      cmut = make_unique<const InsertionHT>(static_cast<InsertionHT&>(*mut));
#else
      cmut = std::make_unique<const InsertionHT>(static_cast<InsertionHT&>(*mut));
#endif
      break;
    case REPL_HT:
#if __cplusplus == 201103L
      cmut = make_unique<const ReplacementHT>(static_cast<ReplacementHT&>(*mut));
#else
      cmut = std::make_unique<const ReplacementHT>(static_cast<ReplacementHT&>(*mut));
#endif
      break;
    default:
#if __cplusplus == 201103L
      cmut = make_unique<const ReplacementHT>(static_cast<ReplacementHT&>(*mut));
#else
      cmut = std::make_unique<const ReplacementHT>(static_cast<ReplacementHT&>(*mut));
#endif
      break;
  }
  ht_.push_back(std::move(cmut));
  nb_mut_[mut->mut_type()]++;
}


/// Useful when we inspect a tree file
/// because stats are not saved in the file.
void DnaReplicationReport::compute_stats()
{
  nb_mut_[SWITCH] = 0;
  nb_mut_[S_INS]  = 0;
  nb_mut_[S_DEL]  = 0;
  nb_mut_[DUPL]   = 0;
  nb_mut_[DEL]    = 0;
  nb_mut_[TRANS]  = 0;
  nb_mut_[INV]    = 0;
  nb_mut_[INS_HT] = 0;
  nb_mut_[REPL_HT]= 0;

  for (const auto& ht : ht_) {
    assert(ht->mut_type() == INS_HT or
           ht->mut_type() == REPL_HT);
    nb_mut_[ht->mut_type()]++;
  }

  for (const auto& rear : rearrangements_) {
    assert(rear->mut_type() == DUPL or
           rear->mut_type() == DEL or
           rear->mut_type() == TRANS or
           rear->mut_type() == INV);
    nb_mut_[rear->mut_type()]++;
  }

  for (const auto& mut : mutations_) {
    assert(mut->mut_type() == SWITCH or
           mut->mut_type() == S_INS or
           mut->mut_type() == S_DEL);
    nb_mut_[mut->mut_type()]++;
  }
}

void DnaReplicationReport::write_to_tree_file(gzFile tree_file) const {
  // Write the mutations and rearrangements undergone during replication
  // Store HT
  int32_t nb_HT = nb(H_T);
  gzwrite(tree_file, &nb_HT, sizeof(nb_HT));
  for (const auto& ht : ht_) {
    switch(ht->mut_type()) {
      case INS_HT:
        ht->save(tree_file);
        break;
      case REPL_HT:
        ht->save(tree_file);
        break;
      default:
        ht->save(tree_file);
        break;
    }
  }


  // Store rearrangements
  int32_t nb_rears = nb(REARR);
  gzwrite(tree_file, &nb_rears, sizeof(nb_rears));
  for (const auto& rear : rearrangements_) {
    switch(rear->mut_type()) {
      case DUPL:
        rear->save(tree_file);
        break;
      case DEL:
        rear->save(tree_file);
        break;
      case TRANS:
        rear->save(tree_file);
        break;
      case INV:
        rear->save(tree_file);
        break;
      default:
        rear->save(tree_file);
        break;
    }
  }

  // Store mutations
  int32_t nb_muts = nb(S_MUT);
  gzwrite(tree_file, &nb_muts, sizeof(nb_muts));
  for (const auto& mut : mutations_)
    switch(mut->mut_type()) {
      case SWITCH:
        mut->save(tree_file);
        break;
      case S_DEL:
        mut->save(tree_file);
        break;
      case S_INS:
        mut->save(tree_file);
        break;
      default:
        mut->save(tree_file);
        break;
    }
}

void DnaReplicationReport::read_from_tree_file(gzFile tree_file) {
  int32_t nb_rears, nb_muts, nb_HT;
    Mutation* mut = nullptr;

  gzread(tree_file, &nb_HT, sizeof(nb_HT));
  for (int i = 0 ; i < nb_HT ; i++) {
      mut = Mutation::Load(tree_file);
      add_HT(mut);
      delete mut;
  }


  gzread(tree_file, &nb_rears, sizeof(nb_rears));
  for (int i = 0 ; i < nb_rears ; i++) {
      mut = Mutation::Load(tree_file);
      add_rear(mut);
      delete mut;
  }

  gzread(tree_file, &nb_muts, sizeof(nb_muts));
  for(int i = 0 ; i < nb_muts ; i++) {
      mut = Mutation::Load(tree_file);
      add_mut(mut);
      delete mut;
  }
}
} // namespace aevol
