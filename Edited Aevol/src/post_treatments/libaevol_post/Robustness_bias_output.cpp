//
// Created by duazel on 27/04/2020.
//

#include "Robustness_bias_output.h"
#include <iostream>
#include <iomanip>

Robustness_bias_output::Robustness_bias_output(const Individual & indiv, const char* filename_indiv,
                                               const char* filename_mutation)
    : file_indiv_(filename_indiv), file_mutations_(filename_mutation),
      file_duplications_("duplications.csv"), individual_(indiv){
  /*file_indiv_ << "id,";
  file_indiv_ << "seq_length,";
  file_indiv_ << "has_mutate,";
  file_indiv_ << "is_neutral" << std::endl;*/

  file_mutations_ << "id,";
  file_mutations_ << "seq_length,";
  file_mutations_ << "has_mutate,";
  file_mutations_ << "is_neutral,";
  file_mutations_ << "type,";
  file_mutations_ << "impact_on_length,";
  file_mutations_ << "mutation_neutral" << ",";
  file_mutations_ << "nb_non_coding_RNA" << std::endl;

  /*file_duplications_ << "id,";
  file_duplications_ << "is_neutral,";
  file_duplications_ << "nb_non_coding_RNA,";
  file_duplications_ << "nb_touched_genes," << ",";
  file_duplications_ << "percent_touched" << std::endl;*/
}

Robustness_bias_output::~Robustness_bias_output() {
  file_mutations_.close();
  file_indiv_.close();
  file_duplications_.close();
}

int Robustness_bias_output::nb_touched_gene(const Individual & individual, int d1, int f1, int size) {
  int count = 0;

  for (auto prot : individual.protein_list()) {
    if (is_touched_gene(prot, d1, f1, size)) count++;
  }

  return count;
}

bool Robustness_bias_output::is_touched_gene(Protein *prot, int d1, int f1, int size) {
  int d2,f2;

  if (prot->strand() == Strand::LEADING) {
    d2 = prot->shine_dal_pos();
    f2 = prot->last_translated_pos();
  } else {
    f2 = prot->shine_dal_pos();
    d2 = prot->last_translated_pos();
  }

  if (f1<d1) f1 += size;
  if (d2<d1) d2 += size;
  if (f2<d1) f2 += size;

  return !(f1 < d2 && d2 < f2);
}

double Robustness_bias_output::mean_percent_touched(const Individual & individual, int d1, int f1, int size) {
  double sum = 0;
  int count = 0;

  for (auto prot : individual.protein_list()) {
    int d2,f2;
    if (prot->strand() == Strand::LEADING) {
      d2 = prot->shine_dal_pos();
      f2 = prot->last_translated_pos();
    } else {
      f2 = prot->shine_dal_pos();
      d2 = prot->last_translated_pos();
    }
    double value = percent_touched(d1,f1,d2,f2,size);
    if(value < -100 || value > 100) {
    }
    sum += value;
    if (value>0) count++;
  }

  return sum/count;
}

double Robustness_bias_output::percent_touched(int d1, int f1, int d2, int f2, int size) {
  if (f1<d1) f1 += size;
  if (d2<d1) d2 += size;
  if (f2<d1) f2 += size;

  if (f1<d2 && d2<f2) return 0;
  if (f1<f2 && f2<d2) return 100.0*(f1-d1)/(f2+size-d2);
  if (d2<f1 && f1<f2) return 100.0*(f1-d2)/(f2-d2);
  if (d2<f2 && f2<f1) return 100;
  if (f2<f1 && f1<d2) return 100.0*(f2-d1)/(f2+size-d2);
  if (f2<d2 && d2<f1) return 100.0*(f2-d1+f1-d2)/(f2+size-d2);
  return -1;
}

/*void Robustness_bias_output::add_duplication(int id,
                                       bool is_neutral,
                                       int16_t nb_non_coding_RNA,
                                       int nb_touched_genes,
                                       double percent_touched) {
  file_duplications_ << id << ",";
  file_duplications_ << is_neutral << ",";
  file_duplications_ << nb_non_coding_RNA << ",";
  file_duplications_ << nb_touched_genes << ",";
  file_duplications_ << percent_touched << std::endl;
}

void Robustness_bias_output::add_indiv(int id,
                                       int seq_length,
                                       bool has_mutate,
                                       bool is_neutral) {
  file_indiv_ << id << ",";
  file_indiv_ << seq_length << ",";
  file_indiv_ << has_mutate << ",";
  file_indiv_ << is_neutral << std::endl;

  nb_replication++;
  if (has_mutate) nb_mutant++;
  if (is_neutral) nb_neutre++;
}*/

void Robustness_bias_output::add_mutation(int id,
                                          int seq_length,
                                          bool has_mutate,
                                          bool is_neutral,
                                          int type,
                                          int impact_on_length,
                                          bool mutation_neutral,
                                          int16_t nb_non_coding_RNA) {
  std::string type_name[] = {"switch", "small_insertion", "small_deletion", "duplication",
                       "deletion", "translocation", "inversion"};
  nb_mutations++;
  nb_mutations_type[type]++;
  file_mutations_ << id << ";";
  file_mutations_ << seq_length << ";";
  file_mutations_ << has_mutate << ";";
  file_mutations_ << is_neutral << ";";
  file_mutations_ << type << ";";
  file_mutations_ << impact_on_length << ";";
  file_mutations_ << mutation_neutral << ";";
  file_mutations_ << nb_non_coding_RNA << std::endl;

  if(mutation_neutral) {
    nb_neutral_mutations++;
    sum_size_variation_neutral_mutation += impact_on_length;
    nb_neutral_mutations_type[type]++;
  }
}


int32_t impact_on_length(MutationEvent &mutationEvent,
                                       unsigned int size) {
  int32_t ret = 0;

  if (mutationEvent.type() == MutationEventType::DELETION
      || mutationEvent.type() == MutationEventType::DUPLICATION) {
    if (mutationEvent.pos_2() > mutationEvent.pos_1()) {
      ret = mutationEvent.pos_2() - mutationEvent.pos_1();
    } else {
      ret = mutationEvent.pos_2() + size - mutationEvent.pos_1(); //bug correction 17/05/2022 : switch pos_1 and pos_2
    }
  } else {
    if (mutationEvent.type() == MutationEventType::SMALL_INSERTION
        || mutationEvent.type() == MutationEventType::SMALL_DELETION) {
      ret = mutationEvent.number();
    }
  }
  if (mutationEvent.type() == MutationEventType::DELETION
      || mutationEvent.type() == MutationEventType::SMALL_DELETION){
    return -ret;
  } else {
    return ret;
  }
}

void Robustness_bias_output::record_replication(unsigned int id,
    Individual& child, Individual& parent, std::vector<MutationEvent*>& mutations, std::vector<bool>& mut_is_neutral) {
  bool has_mutate;
  bool is_neutral;
  if(mutations.empty()){
    has_mutate = false;
    is_neutral = true;
  } else {
    //std::cout << "Mutations Not Empty" << std::endl;
    has_mutate = true;
    child.clear_everything_except_dna_and_promoters();
    child.compute_phenotype();
    is_neutral = individual_.phenotype()->is_identical_to(*child.phenotype(),0);
  }

  //add_indiv(id, child.amount_of_dna(), has_mutate, is_neutral);
  nb_replication++;
   if(is_neutral) nb_neutre++;
   if(has_mutate) nb_mutant++;
  auto mut_neutral = mut_is_neutral.begin();

  bool unneutral_mut = false;
  for (auto mut = mutations.begin(); mut != mutations.end(); mut++, mut_neutral++) {
    int length_mut = impact_on_length(**mut, parent.amount_of_dna()); // bug correction 16/05/2022 : replace "child." by "parent."
    child.compute_statistical_data();
    add_mutation(id, child.amount_of_dna(), has_mutate, is_neutral, (*mut)->type(), length_mut, *mut_neutral, child.nb_non_coding_RNAs());

    
    /*if ((*mut)->type() == MutationEventType::DUPLICATION) {
      int d1 = (*mut)->pos_1();
      int f1 = (*mut)->pos_2();
      int size = child.amount_of_dna();
      add_duplication(id,is_neutral, child.nb_non_coding_RNAs(),
                      nb_touched_gene(parent,d1,f1,size),
                      mean_percent_touched(parent,d1,f1,size));
    }*/

    if (!(*mut_neutral)) {
      unneutral_mut = true;
    }
  }
  if (mutations.size()>1) {
    multiple_mutations++;
    if (unneutral_mut && is_neutral) {
      mult_unneutral_mut_into_neutral_replication++;
    }
  }
}

void Robustness_bias_output::print_summary(std::ostream& os) {
  std::string type_name[] = {"switch", "small_insertion", "small_deletion", "duplication",
                             "deletion","translocation", "inversion"};
  /*  os << "------------- Summary after " << nb_replication << " replications -------------\n";
  os << "nb neutral replications = " << nb_neutre << "  |  " << ((double) nb_neutre)/nb_replication*100 << "%\n";
  os << "nb mutants              = " << nb_mutant << "  |  " << ((double) nb_mutant)/nb_replication*100 << "%\n";
  os << "\n";
  os << "+---------------+---------------+---------------+---------------+\n";
  os << "| mutation type |  nb  neutral  | nb  mutations |   % neutral   |\n";
  os << "+---------------+---------------+---------------+---------------+\n";
  for (int i=0; i<7; i++) {
    os << "|" << std::setw(15) << type_name[i] << "|" << std::setw(15) << nb_neutral_mutations_type[i] << "|" << std::setw(15) << nb_mutations_type[i];
    os << "|" << std::setw(15) << std::setprecision(2) << ((double)nb_neutral_mutations_type[i])*100 / nb_mutations_type[i] << "|\n";
    os << "+---------------+---------------+---------------+---------------+\n";
  }
  os << "|" << std::setw(15) << "Total" << "|" << std::setw(15) << nb_neutral_mutations << "|" << std::setw(15) << nb_mutations;
  os << "|" << std::setw(15) << std::setprecision(2) << std::fixed << ((double)nb_neutral_mutations)*100 / nb_mutations << "|\n";
  os << "+---------------+---------------+---------------+---------------+\n";
  os << std::setw(1);
  os << "nb mutations            = " << nb_mutations << "\n";
  os << "nb neutral mutations    = " << nb_neutral_mutations << " | " << ((double) nb_neutral_mutations)/nb_mutations*100 << "%\n";
  os << "mean size variation on neutral mutation = " << ((double) sum_size_variation_neutral_mutation)/nb_neutral_mutations << "\n";
  os << "mean size variation on neutral replication = " << ((double) sum_size_variation_neutral_mutation)/nb_neutre << "\n";
  os << "\n";
  os << "nb multiple mutations   = " << multiple_mutations << "\n";
  os << "nb multiple unneutral mutations creating neutral replication = " << mult_unneutral_mut_into_neutral_replication <<"\n";
  os << "------------------------------------------------------------------------" << std::endl;*/

os << nb_replication << ";";
  os <<  nb_neutre << ";";
  os <<  nb_mutant << ";";
  os <<  nb_mutations << ";" << nb_neutral_mutations << ";";
  os <<  sum_size_variation_neutral_mutation << ";" << multiple_mutations << ";" << mult_unneutral_mut_into_neutral_replication;
  for (int i=0; i<7; i++) {
    os << ";" << nb_mutations_type[i] << ";" << nb_neutral_mutations_type[i] ;
  }
  os << std::endl;
  
}
