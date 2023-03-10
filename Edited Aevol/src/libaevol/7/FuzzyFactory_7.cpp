// ****************************************************************************
//
//          Aevol - An in silico experimental evolution platform
//
// ****************************************************************************
//
// Copyright: See the AUTHORS file provided with the package or <www.aevol.fr>
// Web: http://www.aevol.fr/
// E-mail: See <http://www.aevol.fr/contact/>
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons, Jonathan Rouzaud-Cornabas
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

#include "FuzzyFactory_7.h"
#include "AbstractFuzzy_7.h"
#include <list>
#include <omp.h>
#include "Vector_Fuzzy.h"
#include "List_Fuzzy.h"
#include "Discrete_Double_Fuzzy.h"
namespace aevol {
    AbstractFuzzy_7* FuzzyFactory_7::createFuzzy() {
        AbstractFuzzy_7* fuzz = nullptr;
        switch (flavor_)
        {
            case FuzzyFactory_Flavor::LIST:
                //fuzz = new List_Fuzzy();
                printf("List Fuzzy are currently not available\n");
                exit(-1);
                break;
            case FuzzyFactory_Flavor::VECTOR:
                fuzz = new Vector_Fuzzy();
                break;
            case FuzzyFactory_Flavor::DISCRETE_DOUBLE_TABLE:
                // printf("Size %d \n",PHENOTYPE_VECTOR_SIZE);
                fuzz = new Discrete_Double_Fuzzy(PHENOTYPE_VECTOR_SIZE);
                break;
            default:
                fuzz = new Vector_Fuzzy();
                break;
        }
        return fuzz;
    }
    void FuzzyFactory_7::init(int init_size, int pop_size) {
        int n = 0;
            #pragma omp parallel reduction(+:n)
            {
                n += 1;
            }

            nb_local_pool = n;

            if (pop_size != -1) {
                local_pool_size_ = (pop_size / n) * 4;
                global_pool_size_ = pop_size * 2;
            }

            local_list_unused_fuzzy_.resize(n);

            for (int j = 0; j < n; j++) { 
                for (int i = 0; i < local_pool_size_; i++) {
                    local_list_unused_fuzzy_[j].push_back(createFuzzy());
                }
            }

        // printf("SAMPLING %d :: %d :: %d\n",pool_size_,sampling_,PHENOTYPE_VECTOR_SIZE);
        for (int i = 0; i < global_pool_size_; i++) {
            list_unused_fuzzy_.push_back(createFuzzy());
        }
    }

    AbstractFuzzy_7 *FuzzyFactory_7::get_fuzzy() {
            AbstractFuzzy_7 *pop = nullptr;

            #ifdef _OPENMP
            if (local_list_unused_fuzzy_[omp_get_thread_num()].empty()) {
            #endif
                // #pragma omp atomic
                // empty_global+=1;

                if (list_unused_fuzzy_.empty()) {
                    pop = createFuzzy();
                } else {
                    // Go to global (and lock)
                    #pragma omp critical(pop_fuzzy)
                    {
                        if (list_unused_fuzzy_.empty()) {
                            pop = createFuzzy();
                        } else {    
                            pop = list_unused_fuzzy_.front();
                            list_unused_fuzzy_.pop_front();
                        }
                    }    
                }
            #ifdef _OPENMP
            } else {
                #pragma omp critical(pop_fuzzy)
                    {
                        pop = local_list_unused_fuzzy_[omp_get_thread_num()].front();
                        local_list_unused_fuzzy_[omp_get_thread_num()].pop_front();
                    }
            }
            #endif

            return pop;
    }

void FuzzyFactory_7::give_back(AbstractFuzzy_7 *fuzz) {
  fuzz->clear();
   #ifdef _OPENMP
          if (local_list_unused_fuzzy_[omp_get_thread_num()].size() == local_pool_size_) {
    
            // #pragma omp atomic
            // free_space_local++;

            AbstractFuzzy_7* pop = local_list_unused_fuzzy_[omp_get_thread_num()].back();
            local_list_unused_fuzzy_[omp_get_thread_num()].pop_back();

            #pragma omp critical(pop_fuzzy)
            {
                if (list_unused_fuzzy_.size() == global_pool_size_)
                    delete pop;
                else
                    list_unused_fuzzy_.push_back(pop);
            }  
         } 
        
        // #pragma omp atomic
        // free_local++;
        local_list_unused_fuzzy_[omp_get_thread_num()].push_front(fuzz);
    #else
            #pragma omp critical(pop_fuzzy)
            {
                if (list_unused_fuzzy_.size() == global_pool_size_)
                    delete fuzz;
                else
                    list_unused_fuzzy_.push_back(fuzz);
            }
    #endif
}
    }