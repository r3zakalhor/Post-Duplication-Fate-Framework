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

#ifndef AEVOL_DNAMUTATOR_H
#define AEVOL_DNAMUTATOR_H

#include "Individual.h"
#include "JumpingMT.h"
#include "MutationEvent.h"

namespace aevol {
class DnaMutator {
 public:
    DnaMutator(std::shared_ptr<JumpingMT> mut_prng, int32_t length,
               double duplication_rate,
               double deletion_rate, double translocation_rate,
               double inversion_rate,
               double point_mutation_rate, double small_insertion_rate,
               double small_deletion_rate, int16_t max_indel_size,
               int32_t min_genome_length, int32_t max_genome_length, int indiv_id, int x, int y);

    DnaMutator(Individual* indiv, int x, int y);

    ~DnaMutator() {
      int cpt = 0;
      for (auto repl : mutation_list_) {
        delete repl;
        cpt++;
      }

      mutation_list_.clear();
    }

    void generate_mutations();
    void generate_rearrangements();
    void generate_small_mutations();

    MutationEvent* generate_next_mutation(int32_t length = -1);
    void generate_all_mutations(int32_t length);
    bool mutation_available() { return ((cpt_rear_ + cpt_mut_) > 0); }

    std::list<MutationEvent*> mutation_list_;

    bool hasMutate() {return hasMutate_;}

    void setMutate(bool mutate) {hasMutate_ = mutate;}

#ifdef HAVE_MPI
    bool isAtBorder() {return is_at_border_;}

    void setIsAtBorder(bool mutate) {is_at_border_ = mutate;}
#endif

    int x_,y_;

 //private:
    std::shared_ptr<JumpingMT> mut_prng_;
    int32_t length_;

    // ------------------------------ Rearrangement rates (without alignements)
    double duplication_rate_;
    double deletion_rate_;
    double translocation_rate_;
    double inversion_rate_;

    // --------------------------------------------------------- Mutation rates
    double point_mutation_rate_;
    double small_insertion_rate_;
    double small_deletion_rate_;
    int16_t max_indel_size_;

    //--------------------------- Mutation counters
    int32_t nb_swi_;
    int32_t nb_ins_;
    int32_t nb_del_;
    int32_t nb_mut_;

    int32_t nb_large_dupl_;
    int32_t nb_large_del_;
    int32_t nb_large_trans_;
    int32_t nb_large_inv_;
    int32_t nb_rear_;

    int32_t cpt_rear_;
    int32_t cpt_mut_;

    int32_t min_genome_length_;
    int32_t max_genome_length_;

    bool hasMutate_ = false;
    unsigned long long int id_;

    #ifdef HAVE_MPI
    bool is_at_border_ = false;
    #endif
};
}


#endif