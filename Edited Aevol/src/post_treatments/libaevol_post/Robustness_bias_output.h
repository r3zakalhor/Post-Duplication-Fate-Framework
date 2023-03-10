//
// Created by duazel on 27/04/2020.
//

#ifndef AEVOL_ROBUSTNESS_BIAS_OUTPUT_H
#define AEVOL_ROBUSTNESS_BIAS_OUTPUT_H

#include "Individual.h"
#include "MutationEvent.h"
#include <fstream>

using namespace aevol;

class Robustness_bias_output {
 public:
  Robustness_bias_output(const Individual & indiv, const char * filename_indiv, const char * filename_mutation);
  virtual ~Robustness_bias_output();

  //void add_indiv(int id, int seq_length, bool has_mutate, bool is_neutral);
  void add_mutation(int id, int seq_length, bool has_mutate, bool is_neutral, int type,
      int impact_on_length, bool mutation_neutral, int16_t nb_non_coding_RNA);
  //void add_duplication(int id, bool is_neutral, int16_t nb_non_coding_RNA, int nb_touched_genes, double percent_touched);

  void record_replication(unsigned int id, Individual& child, Individual& parent, std::vector<MutationEvent*>& mutations, std::vector<bool>& mut_is_neutral);
  void print_summary(std::ostream & os);

 private:
  std::ofstream file_indiv_;
  std::ofstream file_mutations_;
  std::ofstream file_duplications_;
  Individual individual_;

  unsigned int nb_replication=0;
  unsigned int nb_neutre=0;
  unsigned int nb_mutant=0;
  unsigned int nb_mutations=0;
  unsigned int nb_neutral_mutations=0;
  unsigned int nb_mutations_type[7]={0,0,0,0,0,0,0};
  unsigned int nb_neutral_mutations_type[7]={0,0,0,0,0,0,0};
  unsigned int sum_size_neutral_duplication=0;
  unsigned int sum_size_neutral_deletion=0;
  unsigned int multiple_mutations=0;
  unsigned int mult_unneutral_mut_into_neutral_replication = 0;
  long sum_size_variation_neutral_mutation=0;
  double percent_touched(int d1, int f1, int d2, int f2, int size);
  int nb_touched_gene(const Individual& individual, int d1, int f1, int size);
  bool is_touched_gene(Protein* prot, int d1, int f1, int size);
  double
  mean_percent_touched(const Individual& individual, int d1, int f1, int size);
};

#endif //AEVOL_ROBUSTNESS_BIAS_OUTPUT_H
