//
// Created by duazel on 11/03/2020.
//

#ifndef AEVOL_NEUTRAL_MUTATION_EXP_H
#define AEVOL_NEUTRAL_MUTATION_EXP_H

#include <string>

#include "ExpManager.h"
#include "Individual.h"
#include "MutationParams.h"

using namespace aevol;

Individual *new_child(Individual *individual,
                      std::vector<MutationEvent *> &mutations,
                      std::vector<int32_t> &lengths);
Individual *new_child(Individual *individual, std::vector<MutationEvent *> &mutations,
                      std::vector<bool> &is_neutral);

Individual *neutral_mutation(Individual *individual);

Individual * run_generations(unsigned int nb_generations, Individual* indiv);
Individual * run_to_size(int32_t wanted_size, Individual* indiv);
Individual * run_to_homogenize(int32_t nb_gen, Individual* indiv);

#endif //AEVOL_NEUTRAL_MUTATION_EXP_H
