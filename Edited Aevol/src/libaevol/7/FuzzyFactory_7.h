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

#ifndef AEVOL_FUZZYFACTORY_H
#define AEVOL_FUZZYFACTORY_H

#include "AbstractFuzzy_7.h"

namespace aevol {
    enum FuzzyFactory_Flavor {
        LIST = 0,
        VECTOR = 1,
        DISCRETE_FLOAT_TABLE = 2,
        DISCRETE_DOUBLE_TABLE = 3
    };



    class FuzzyFactory_7 {
    public:
     FuzzyFactory_7(int flavor, int pool_size, int sampling, int pop_size = -1) {
         
            flavor_ = (FuzzyFactory_Flavor)flavor;
        //  printf("FuzzyFactory %d\  n",flavor_);
            pool_size_ = pool_size;
            if (flavor_>=2) {
                PHENOTYPE_VECTOR_SIZE=sampling;
                D_PHENOTYPE_VECTOR_SIZE=((double)sampling);
                // printf("SAMPLING %d :: %d\n",sampling,PHENOTYPE_VECTOR_SIZE);

            }
            init(pool_size_,pop_size);
        }

        ~FuzzyFactory_7() {
            for (auto&& it_fuzzy : list_unused_fuzzy_) {
                delete it_fuzzy;
            }

            for (int i = 0; i < nb_local_pool; i++)
                for (auto && it_fuzzy : local_list_unused_fuzzy_[i]) {
                    delete it_fuzzy;
                }
            list_unused_fuzzy_.clear();
        }

        void init(int init_size, int pop_size = -1);


        void stats() {
            int total_length_ = ((flavor_>=2)? PHENOTYPE_VECTOR_SIZE: 1)*list_unused_fuzzy_.size();
            printf("FUZZY_FACTORY_STATS -- Number of Fuzzys %ld - Combined size %d elements\n",list_unused_fuzzy_.size(),
            total_length_);
        }

        AbstractFuzzy_7 *get_fuzzy();

        void give_back(AbstractFuzzy_7 *fuzzy);

    private:
        AbstractFuzzy_7* createFuzzy();

        std::list<AbstractFuzzy_7 *> list_unused_fuzzy_;
        std::vector<std::list<AbstractFuzzy_7 *>> local_list_unused_fuzzy_;

        FuzzyFactory_Flavor flavor_;
        int pool_size_;

        int global_pool_size_ = 256;
        int local_pool_size_ = 256;

        int32_t nb_local_pool = 0;

    int32_t PHENOTYPE_VECTOR_SIZE            = -1;
    double D_PHENOTYPE_VECTOR_SIZE = ((double)PHENOTYPE_VECTOR_SIZE);      
    };

}
#endif //AEVOL_DNAFACTORY_H
