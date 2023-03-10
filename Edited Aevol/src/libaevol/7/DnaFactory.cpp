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

#include "DnaFactory.h"
#include "ExpManager_7.h"
#include "Individual_7.h"

#include <list>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0;}
inline omp_int_t omp_get_max_threads() { return 1;}
#endif
namespace aevol {
    void DnaFactory::init(int init_size, int pop_size) {
        if (policy_ != DnaFactory_Policy::LOCAL_GLOBAL_FIT) {
            if (policy_ != DnaFactory_Policy::ALLOCATE) {
                for (int i = 0; i < pool_size_; i++) {
                    list_unused_dna_.push_back(new Dna_7(init_size,this));
                }
            }
        } else {
            int n = 0;
            #pragma omp parallel reduction(+:n)
            {
                n += 1;
            }

            // printf("NB OpenMP : %d\n",n);

            local_list_unused_dna_.resize(n);
            if (pop_size != -1) {
                local_pool_size_ = (pop_size / n) * 4;
                global_pool_size_ = pop_size * 2;
            }

            for (int j = 0; j < n; j++) { 
                for (int i = 0; i < local_pool_size_; i++) {
                    local_list_unused_dna_[j].push_back(new Dna_7(init_size,this));
                }
            }

            for (int i = 0; i < global_pool_size_; i++) {
                list_unused_dna_.push_back(new Dna_7(init_size,this));
            }
        }
    }

    Dna_7 *DnaFactory::get_dna(int request_size) {
        request_size++; // Count the end \0
        int req_block = Dna_7::nb_blocks(request_size);
        if (policy_ == DnaFactory_Policy::FIRST) {
          Dna_7 *pop = nullptr;
          bool allocate = false;
            #pragma omp critical(pop_dna)
            {
                if (list_unused_dna_.empty()) {
                    allocate = true;
                } else {
                    pop = list_unused_dna_.front();
                    list_unused_dna_.pop_front();
                }
            }
            if (allocate)
                pop = new Dna_7(request_size,this);
            
            return pop;
            // return new Dna_7(request_size,this);
        } else if (policy_ == DnaFactory_Policy::FIRSTFIT) {
          Dna_7 *pop = nullptr;

#pragma omp critical(pop_dna)
            {
                if (list_unused_dna_.empty()) {
                    pop = new Dna_7(request_size,this);
                } else {
                    std::list<Dna_7 *>::iterator found_it;
                    for (auto it = list_unused_dna_.begin(); it != list_unused_dna_.end(); it++) {
                        if ((*it)->nb_block() >= req_block) {
                            found_it = it;
                            pop = (*it);
                            break;
                        }
                    }

                    if (pop == nullptr) {
                        pop = list_unused_dna_.front();
                        list_unused_dna_.pop_front();
                    } else {
                        list_unused_dna_.erase(found_it);
                    }
                }
            }
            return pop;
        } else if (policy_ == DnaFactory_Policy::LOCAL_GLOBAL_FIT) {
            Dna_7 *pop = nullptr;

            if (local_list_unused_dna_[omp_get_thread_num()].empty()) {
                // #pragma omp atomic
                // empty_global+=1;

                if (list_unused_dna_.empty()) {
                    pop = new Dna_7(request_size,this);
                } else {
                    // Go to global (and lock)
                    #pragma omp critical(pop_dna)
                    {
                        // printf("Search global array of DNA\n");
                        std::list<Dna_7 *>::iterator found_it;
                        for (auto it = list_unused_dna_.begin(); it != list_unused_dna_.end(); it++) {
                            if ((*it)->nb_block() >= req_block) {
                                found_it = it;
                                pop = (*it);
                                break;
                            }
                        }

                        if (pop == nullptr) {
                            pop = list_unused_dna_.front();
                            list_unused_dna_.pop_front();
                        } else {
                            list_unused_dna_.erase(found_it);
                        }
                    }    
                }
            } else {
                // Go to local (and lock free)
                std::list<Dna_7 *>::iterator found_it;

                // #pragma omp atomic
                // local_empty++;

                for (auto it = local_list_unused_dna_[omp_get_thread_num()].begin(); 
                            it != local_list_unused_dna_[omp_get_thread_num()].end(); it++) {
                    if ((*it)->nb_block() >= req_block) {
                        found_it = it;
                        pop = (*it);
                        break;
                    }
                }

                if (pop == nullptr) {
                    pop = local_list_unused_dna_[omp_get_thread_num()].front();
                    local_list_unused_dna_[omp_get_thread_num()].pop_front();
                } else {
                    local_list_unused_dna_[omp_get_thread_num()].erase(found_it);
                }
            }
            return pop;
        } else if (policy_ == DnaFactory_Policy::ALLOCATE) {
            return new Dna_7(request_size,this);
        }
        return nullptr;
    }

void DnaFactory::give_back(Dna_7 *dna) {
  dna->reset_stat();
//   delete dna;
    if (policy_ == DnaFactory_Policy::LOCAL_GLOBAL_FIT) {
        if (local_list_unused_dna_[omp_get_thread_num()].size() == local_pool_size_) {
            
            // #pragma omp atomic
            // free_space_local++;

            Dna_7* pop = local_list_unused_dna_[omp_get_thread_num()].back();
            local_list_unused_dna_[omp_get_thread_num()].pop_back();

            #pragma omp critical(pop_dna)
            {
                list_unused_dna_.push_back(pop);
            }  
        } 
        
        // #pragma omp atomic
        // free_local++;
            
        local_list_unused_dna_[omp_get_thread_num()].push_front(dna);
    } else if (policy_ == DnaFactory_Policy::ALLOCATE) {
        delete dna;    
    } else {
        // #pragma omp atomic
        // free_global++;
            
        #pragma omp critical(pop_dna)
        {
            list_unused_dna_.push_back(dna);
        }
    }
}

void DnaFactory::reduce_space(ExpManager_7* exp_m) {
    if (garbage_policy == DnaFactory_Garbage_Policy::MAXSIZE) {
        int32_t max_size = 0;

        for (int32_t indiv_id = 0; indiv_id < exp_m->nb_indivs_; indiv_id++) {
            max_size = exp_m->previous_individuals[indiv_id]->dna_->length_ > max_size ? 
                            exp_m->previous_individuals[indiv_id]->dna_->length_ : max_size;
        }

        for (int j = 0; j < local_list_unused_dna_.size(); j++) { 
            for (std::list<Dna_7 *>::iterator it_dna = local_list_unused_dna_[j].begin();
                    it_dna != local_list_unused_dna_[j].end(); ) {
                bool toRemove = false;
                auto it_remove = it_dna;
                if ((*it_dna)->length_ > max_size) {
                    toRemove = true;
                }
                it_dna++;
                if (toRemove) {
                    printf("Delete %d (max %d) from local pool %d\n",(*it_remove)->length_,max_size,j);
                    delete (*it_remove);
                    local_list_unused_dna_[j].erase(it_remove);
                }
            }
        }

        for (std::list<Dna_7 *>::iterator it_dna = list_unused_dna_.begin();
                    it_dna != list_unused_dna_.end();) {
                bool toRemove = false;
                auto it_remove = it_dna;
                if ((*it_dna)->length_ > max_size) {
                    toRemove = true;
                    printf("MARK -- Delete %d (max %d) from global pool\n",(*it_remove)->length_,max_size);

                }
                
                it_dna++;

                if (toRemove) {
                    printf("Delete %d (max %d) from global pool\n",(*it_remove)->length_,max_size);
                    delete (*it_remove);
                    list_unused_dna_.erase(it_remove);

                }
        }
    } else if (garbage_policy == DnaFactory_Garbage_Policy::MAXMEM) {
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

        printf("total length %zu ou of %zu\n",total_length_,max_pool_size);

        if (total_length_ > max_pool_size) {
            size_t to_remove_size = total_length_ - max_pool_size;

            list_unused_dna_.sort([](Dna_7 *a, Dna_7 *b) {
                    return a->length_ > b->length_;
                });

            for (std::list<Dna_7 *>::iterator it_dna = list_unused_dna_.begin();
                    it_dna != list_unused_dna_.end();) {
                auto it_remove = it_dna;
                it_dna++;

                to_remove_size -= ((*it_remove)->length_ / 1000000);

                printf("Delete %d (remaining to remove %zu / max pool size %zu tt %zu) from global pool\n",
                                            ((*it_remove)->length_ / 1000000),to_remove_size,max_pool_size,total_length_);                    
                delete (*it_remove);
                list_unused_dna_.erase(it_remove);
                
                if (to_remove_size < 0) break;
            }

            if (to_remove_size > 0) {
                std::list<Dna_7 *>::iterator list_iter[local_list_unused_dna_.size()];
                
                for (int j = 0; j < local_list_unused_dna_.size(); j++) {   
                    local_list_unused_dna_[j].sort([](Dna_7 *a, Dna_7 *b) {
                        return a->length_ > b->length_;
                    });
                    list_iter[j] = local_list_unused_dna_[j].begin();
                }

                while (to_remove_size > 0) {
                
                    for (int j = 0; j < local_list_unused_dna_.size(); j++) { 
                        if (list_iter[j] != local_list_unused_dna_[j].end()) {
                            auto it_remove = list_iter[j];
                            list_iter[j]++;
                                        
                            to_remove_size -= ((*it_remove)->length_ / 1000000);

                            printf("Delete %d (remaining to remove %zu / max pool size %zu tt %zu) from global pool\n",
                                            ((*it_remove)->length_ / 1000000),to_remove_size,max_pool_size,total_length_);
                            
                            delete (*it_remove);
                            local_list_unused_dna_[j].erase(it_remove);
                        }

                        if (to_remove_size < 0) break;
                    }

                }
            }
        }

        delete [] a_total_length_;
    }

}

}