//
// Created by duazel on 11/03/2020.
//

#include "neutral_mutation_exp.h"

#include "Dna.h"
#include "DnaMutator.h"
#include "ExpSetup.h"
#include "FuzzyFactory.h"
#include "Individual.h"
#include "JumpingMT.h"
#include "MutationEvent.h"
#include "MutationParams.h"
#include "neutral_mutation_output.h"

#include <cmath>
#include <cstring>
#include <iostream>
#include <utility>
#include <vector>

using namespace aevol;

void perform_mutation(MutationEvent *mutEvent, Dna *dna) {
  switch (mutEvent->type()) {
  case MutationEventType::DELETION :
    dna->do_deletion(mutEvent->pos_1(), mutEvent->pos_2());
    break;

  case MutationEventType::DO_SWITCH :
  #ifdef BASE_2
    dna->do_switch(mutEvent->pos_1());
  #elif BASE_4
    dna->do_switch(mutEvent->pos_1(), mutEvent->base());
  #endif
    break;

  case MutationEventType::DUPLICATION :
    dna->do_duplication(mutEvent->pos_1(), mutEvent->pos_2(),
                        mutEvent->pos_3());
    break;

  case MutationEventType::INVERSION :
    dna->do_inversion(mutEvent->pos_1(), mutEvent->pos_2());
    break;

  case MutationEventType::SMALL_DELETION :
    dna->do_small_deletion(mutEvent->pos_1(), mutEvent->number());
    break;

  case MutationEventType::SMALL_INSERTION :
    dna->do_small_insertion(mutEvent->pos_1(), mutEvent->number(),
                            mutEvent->seq());
    break;

  case MutationEventType::TRANSLOCATION :
    dna->do_translocation(mutEvent->pos_1(), mutEvent->pos_2(),
                          mutEvent->pos_3(), mutEvent->pos_4(),
                          mutEvent->invert());
  }
}

Individual *new_child(Individual *individual, std::vector<MutationEvent *> &mutations,
                      std::vector<int32_t> &lengths) {
  auto child = new Individual(*individual);
  DnaMutator dnaMutator(child, 0, 0);
  dnaMutator.generate_mutations();
  MutationEvent *mutEvent;
  mutEvent = dnaMutator.generate_next_mutation(child->amount_of_dna());
  Dna *dna = child->genetic_unit(0).dna();
  while (mutEvent) {
    lengths.emplace_back(child->amount_of_dna());
    perform_mutation(mutEvent, dna);
    mutations.emplace_back(new MutationEvent(*mutEvent));
    mutEvent = dnaMutator.generate_next_mutation(child->amount_of_dna());
  }
  return child;
}

Individual *new_child(Individual *individual, std::vector<MutationEvent *> &mutations,
                      std::vector<bool> &is_neutral) {
  auto child = new Individual(*individual);
  DnaMutator dnaMutator(child, 0, 0);
  dnaMutator.generate_mutations();
  MutationEvent *mutEvent;
  mutEvent = dnaMutator.generate_next_mutation(child->amount_of_dna());
  Dna *dna = child->genetic_unit(0).dna();
  while (mutEvent) {
    perform_mutation(mutEvent, dna);
    child->clear_everything_except_dna_and_promoters();
    child->compute_phenotype();
    is_neutral.emplace_back(individual->phenotype()->is_identical_to(*child->phenotype(),0));
    mutations.emplace_back(new MutationEvent(*mutEvent));
    mutEvent = dnaMutator.generate_next_mutation(child->amount_of_dna());
  }
  return child;
}

Individual *neutral_mutation(Individual *individual, unsigned int current_generation) {
  Individual *child;
  std::vector<MutationEvent *> mutations;
  std::vector<int32_t> lengths;
  while (true) {
    child = new_child(individual, mutations, lengths);
    if (mutations.empty()) {
      delete child;
      return individual;
    } else {
      child->clear_everything_except_dna_and_promoters();
      child->compute_phenotype();
      if (individual->phenotype()->is_identical_to(*child->phenotype(), 0)) {
        break;
      } else {
        mutations.clear();
        delete child;
      }
    }
  }
  typedef std::pair<std::vector<MutationEvent *>::iterator, std::vector<int32_t>::iterator> pairIter;
  for (pairIter p(mutations.begin(), lengths.begin()); p.first != mutations.end();
       p.first++, p.second++) {
    out::record_mutation_event(**p.first, current_generation, *(p.second));
  }
  delete individual;
  return child;
}

Individual *neutral_mutation_closer_to(Individual *individual, int32_t size_wanted, unsigned int current_generation) {
  Individual *child;
  std::vector<MutationEvent *> mutations;
  std::vector<int32_t> lengths;
  while (true) {
    child = new_child(individual, mutations, lengths);
    if (mutations.empty()) {
      delete child;
      return individual;
    } else {
      child->clear_everything_except_dna_and_promoters();
      child->compute_phenotype();
      if (individual->phenotype()->is_identical_to(*child->phenotype(), 0) &&
          std::abs(size_wanted - individual->amount_of_dna()) >=
              std::abs(size_wanted - child->amount_of_dna())) {
        break;
      } else {
        mutations.clear();
        lengths.clear();
        delete child;
      }
    }
  }
  typedef std::pair<std::vector<MutationEvent *>::iterator, std::vector<int32_t>::iterator> pairIter;
  for (pairIter p(mutations.begin(), lengths.begin()); p.first != mutations.end();
       p.first++, p.second++) {
    out::record_mutation_event(**p.first, current_generation, *(p.second));
  }
  delete individual;
  return child;
}

Individual *neutral_mutation_constant_size(Individual *individual, unsigned int current_generation) {
  Individual *child;
  std::vector<MutationEvent *> mutations;
  std::vector<int32_t> lengths;
  while (true) {
    child = new_child(individual, mutations, lengths);
    if (mutations.empty()) {
      delete child;
      return individual;
    } else {
      child->clear_everything_except_dna_and_promoters();
      child->compute_phenotype();
      if (individual->phenotype()->is_identical_to(*child->phenotype(), 0) && (individual->amount_of_dna() == child->amount_of_dna())) {
        break;
      } else {
        mutations.clear();
        lengths.clear();
        delete child;
      }
    }
  }
  typedef std::pair<std::vector<MutationEvent *>::iterator, std::vector<int32_t>::iterator> pairIter;
  for (pairIter p(mutations.begin(), lengths.begin()); p.first != mutations.end();
       p.first++, p.second++) {
    out::record_mutation_event(**p.first, current_generation, *(p.second));
  }
  delete individual;
  return child;
}


Individual * run_to_size(int32_t wanted_size, Individual* indiv) {
  auto *individual = new Individual(*indiv);
  auto *expSetup = new ExpSetup(nullptr);
  FuzzyFactory::fuzzyFactory = new FuzzyFactory(expSetup);
  individual->clear_everything_except_dna_and_promoters();
  individual->compute_phenotype();
  unsigned int current_generation = 0;
  while (individual->amount_of_dna() != wanted_size) {
    if (current_generation % 1000 == 0) {
      std::cout << current_generation << "," << individual->genetic_unit(0).dna()->length() << std::endl;
    }
    out::result_file << current_generation << "," << individual->genetic_unit(0).dna()->length() << std::endl;
    individual = neutral_mutation_closer_to(individual, wanted_size, current_generation);
    current_generation++;
  }
  out::result_file << current_generation << "," << individual->genetic_unit(0).dna()->length() << std::endl;
  std::cout << current_generation << "," << individual->genetic_unit(0).dna()->length() << std::endl;
  return  individual;
}

Individual * run_to_homogenize(int32_t nb_gen, Individual* indiv) {
  auto *individual = new Individual(*indiv);
  auto *expSetup = new ExpSetup(nullptr);
  FuzzyFactory::fuzzyFactory = new FuzzyFactory(expSetup);
  individual->clear_everything_except_dna_and_promoters();
  individual->compute_phenotype();
  unsigned int current_generation = 0;
  while (current_generation<=nb_gen) {
    if (current_generation % 1000 == 0) {
      std::cout << current_generation << "," << individual->genetic_unit(0).dna()->length() << std::endl;
    }
    out::result_file << current_generation << "," << individual->genetic_unit(0).dna()->length() << std::endl;
    individual = neutral_mutation_constant_size(individual, current_generation);
    current_generation++;
  }
  out::result_file << current_generation << "," << individual->genetic_unit(0).dna()->length() << std::endl;
  std::cout << current_generation << "," << individual->genetic_unit(0).dna()->length() << std::endl;
  return  individual;
}

Individual * run_generations(unsigned int nb_generations, Individual* indiv) {
  auto *individual = new Individual(*indiv);
  auto *expSetup = new ExpSetup(nullptr);
  FuzzyFactory::fuzzyFactory = new FuzzyFactory(expSetup);
  individual->clear_everything_except_dna_and_promoters();
  individual->compute_phenotype();

  for (unsigned int current_generation = 0; current_generation < nb_generations; ++current_generation) {
    if (current_generation % 100 == 0) {
      std::cout << current_generation << "," << individual->genetic_unit(0).dna()->length() << std::endl;
    }
    out::result_file << current_generation << "," << individual->genetic_unit(0).dna()->length() << std::endl;
    individual = neutral_mutation(individual, current_generation);
  }
  out::result_file << nb_generations << "," << individual->genetic_unit(0).dna()->length() << std::endl;
  std::cout << nb_generations << "," << individual->genetic_unit(0).dna()->length() << std::endl;
  return  individual;
}
