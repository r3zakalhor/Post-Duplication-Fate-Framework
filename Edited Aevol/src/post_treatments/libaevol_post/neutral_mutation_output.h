//
// Created by duazel on 12/03/2020.
//

#ifndef AEVOL_NEUTRAL_MUTATION_OUTPUT_H
#define AEVOL_NEUTRAL_MUTATION_OUTPUT_H

#include <fstream>
#include "MutationEvent.h"

using namespace aevol;

namespace out {
  extern std::ofstream result_file;
  extern std::ofstream mutation_file;

  void init(const char *name_result,
            const char *name_mutation);

  void close();

  void record_mutation_event(MutationEvent &mutationEvent, unsigned int number_generation, unsigned int size);
}

#endif //AEVOL_NEUTRAL_MUTATION_OUTPUT_H
