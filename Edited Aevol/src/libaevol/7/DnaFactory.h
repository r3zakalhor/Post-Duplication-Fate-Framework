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

#ifndef AEVOL_DNAFACTORY_H
#define AEVOL_DNAFACTORY_H

#include "Dna_7.h"

#include <vector>

namespace aevol {
    enum DnaFactory_Policy {
        FIRST = 0,
        FIRSTFIT = 1,
        BESTFIT = 2,
        LOCAL_GLOBAL_FIT = 3,
        ALLOCATE = 4
    };

    enum DnaFactory_Garbage_Policy {
        MAXSIZE = 0,
        MAXMEM = 1
    };

    class ExpManager_7;

    class DnaFactory {
    public:
     DnaFactory(DnaFactory_Policy policy, int pool_size, int init_size, int pop_size = -1) {
            policy_ = policy;
            pool_size_ = pool_size;
            // printf("Pool init \n");
            init(init_size,pop_size);
        }

        ~DnaFactory() {
            for (auto&& it_dna : list_unused_dna_) {
                delete it_dna;
            }
            list_unused_dna_.clear();

            for (size_t j = 0; j < local_list_unused_dna_.size(); j++) { 
                for (auto&& it_dna : local_list_unused_dna_[j]) {
                    delete it_dna;
                }   
                local_list_unused_dna_[j].clear();
            }
            local_list_unused_dna_.clear();

        }

        void init(int init_size, int pop_size = -1);

        void stats() {
            size_t total_length_ = 0;
            size_t* a_total_length_ = new size_t[local_list_unused_dna_.size()];
            int32_t nb_dna = 0;
            for (std::list<Dna_7 *>::iterator it_dna = list_unused_dna_.begin();
                 it_dna != list_unused_dna_.end(); it_dna++) {
                total_length_ += (*it_dna)->nb_block() * BLOCK_SIZE * sizeof(char);
                nb_dna++;
            }

            for (size_t j = 0; j < local_list_unused_dna_.size(); j++) { 
                a_total_length_[j] = 0;
                for (std::list<Dna_7 *>::iterator it_dna = local_list_unused_dna_[j].begin();
                    it_dna != local_list_unused_dna_[j].end(); it_dna++) {
                    total_length_ += (*it_dna)->nb_block() * BLOCK_SIZE * sizeof(char);
                    a_total_length_[j]+= (*it_dna)->nb_block() * BLOCK_SIZE * sizeof(char);
                    nb_dna++;
                }
                
                a_total_length_[j] /= 1000000;
            }

            total_length_ /= 1000000;
            printf("DNA_FACTORY_STATS --  Number of DNAs %d - Combined size %zu Mb (%zu):: ",
                //    local_empty,empty_global,
                //    free_local,free_space_local,free_global,
                   nb_dna, total_length_,local_list_unused_dna_.size());
            for (size_t j = 0; j < local_list_unused_dna_.size(); j++) { 
                printf("%zu ",a_total_length_[j]);
            }

            printf("\n");
            // empty_global = 0;
            // local_empty = 0;

            // free_space_local = 0;
            // free_local = 0;
            // free_global = 0;

            delete [] a_total_length_;
        }

        Dna_7 *get_dna(int request_size);

        void give_back(Dna_7 *dna);

        void reduce_space(ExpManager_7* exp_m);

        const static int32_t cleanup_step = 100;

        const static DnaFactory_Garbage_Policy garbage_policy = DnaFactory_Garbage_Policy::MAXMEM;

    private:
        std::list<Dna_7 *> list_unused_dna_;
        std::vector<std::list<Dna_7 *>> local_list_unused_dna_;
        DnaFactory_Policy policy_;
        int pool_size_;

        int global_pool_size_ = 256;
        int local_pool_size_ = 256;

        size_t max_pool_size = 16000;

        // int empty_global = 0;
        // int local_empty = 0;

        // int free_space_local = 0;
        // int free_local = 0;
        // int free_global = 0;
    };

}
#endif //AEVOL_DNAFACTORY_H
