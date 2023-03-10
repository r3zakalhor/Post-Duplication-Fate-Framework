//
// Created by duazel on 12/03/2020.
//

#include "neutral_mutation_output.h"

unsigned int indent_level = 0;

namespace out {
  std::ofstream log;
  std::ofstream result_file;
  std::ofstream mutation_file;
}

void out::init(const char *name_result,
               const char *name_mutation) {
  result_file.open(name_result);
  mutation_file.open(name_mutation);
  printf("opening %s\n",name_mutation);
  indent_level = 0;

  mutation_file << "# ##############################################################################\n";
  mutation_file << "#  Mutations during the accumulations of neutral mutations\n";
  mutation_file << "# ##############################################################################\n";
  mutation_file << "#  1.  Generation       (mut. occurred when producing the indiv. of this\n"
                   "#                        generation)\n";
  mutation_file << "#  2.  Genetic unit     (which underwent the mutation, 0 = chromosome) \n";
  mutation_file << "#  3.  Mutation type    (0: switch, 1: smallins, 2: smalldel, 3:dupl, 4: del,\n"
                   "#                       5:trans, 6:inv, 7:insert)\n";
  mutation_file << "#  4.  pos_0            (position for the small events,\n"
                   "#                       begin_segment for the rearrangements\n";
  mutation_file << "#  5.  pos_1            (-1 for the small events, \n"
                   "#                        end_segment for the rearrangements\n";
  mutation_file << "#  6.  pos_2            (reinsertion point for duplic.,\n"
                   "#                        cutting point in segment for transloc.,\n"
                   "#                        -1 for other events)\n";
  mutation_file << "#  7.  pos_3            (reinsertion point for transloc.,\n"
                   "#                        -1 for other events)\n";
  mutation_file << "#  8.  invert           (transloc, was the segment inverted (0/1)?)\n";
  mutation_file << "#  9.  align_score      (score that was needed for the rearrangement to occur,\n"
                   "#       char                 score of the first alignment for ins_HT and repl_HT)\n";
  mutation_file << "#  10. align_score2     (score for the reinsertion for transloc,\n"
                   "#                        score of the second alignment for ins_HT and repl_HT)\n";
  mutation_file << "#  11. seg_len          (segment length for rearrangement)\n";
  mutation_file << "#  12. repl_seg_len     (replaced segment length for repl_HT,\n"
                   "#                        -1 for the others)\n";
  mutation_file << "#  13. GU_length        (before the event)\n";
  mutation_file << "#  14. Impact of the mutation on the metabolic error\n"
                   "#                       (negative value = smaller gap after = beneficial mutation)\n";
  mutation_file << "#  14. Impact of the mutation on the fitness\n"
                   "#                       (negative value = smaller gap after = beneficial mutation)\n";
  mutation_file << "#  15. Number of coding RNAs possibly disrupted by the breakpoints\n";
  mutation_file << "#  16. Number of coding RNAs completely included in the segment\n"
                   "#                       (donor segment in the case of a transfer)\n";
  mutation_file << "#  17. Number of coding RNAs that were completely included in the replaced segment\n"
                   "#                       (meaningful only for repl_HT) \n";
  mutation_file << "################################################################################\n";
  mutation_file << "#\n";
  mutation_file << "# Header for R\n";
  mutation_file << "gener,gen_unit,mut_type,pos_0,pos_1,pos_2,pos_3,invert,align_score,align_score_2,seg_len"
                   ",repl_seg_len,GU_len,impact,nbgenesatbreak,nbgenesinseg,nbgenesinreplseg"
                << std::endl;
}

void out::close() {
  result_file.close();
  mutation_file.close();
}

int32_t mutation_event_pos_2(MutationEvent &mutationEvent) {
  if (mutationEvent.type() == MutationEventType::DELETION
      || mutationEvent.type() == MutationEventType::DUPLICATION
      || mutationEvent.type() == MutationEventType::INVERSION
      || mutationEvent.type() == MutationEventType::TRANSLOCATION) {
    return mutationEvent.pos_2();
  } else {
    return -1;
  }
}

int32_t mutation_event_pos_3(MutationEvent &mutationEvent) {
  if (mutationEvent.type() == MutationEventType::DUPLICATION
      || mutationEvent.type() == MutationEventType::TRANSLOCATION) {
    return mutationEvent.pos_3();
  } else {
    return -1;
  }
}

int32_t mutation_event_pos_4(MutationEvent &mutationEvent) {
  if (mutationEvent.type() == MutationEventType::TRANSLOCATION) {
    return mutationEvent.pos_4();
  } else {
    return -1;
  }
}

int32_t mutation_event_invert(MutationEvent &mutationEvent) {
  if (mutationEvent.type() == MutationEventType::TRANSLOCATION) {
    return mutationEvent.invert();
  } else {
    return 0;
  }
}

int32_t mutation_event_length_mutation(MutationEvent &mutationEvent,
                                       unsigned int size) {
  if (mutationEvent.type() == MutationEventType::DELETION
      || mutationEvent.type() == MutationEventType::DUPLICATION
      || mutationEvent.type() == MutationEventType::INVERSION
      || mutationEvent.type() == MutationEventType::TRANSLOCATION){
    int intern_mut_size = std::max(mutationEvent.pos_1()-mutationEvent.pos_2() , mutationEvent.pos_2()-mutationEvent.pos_1()); //length of the segment that do not cover the origin
    int extern_mut_size = size - std::max(mutationEvent.pos_1(),mutationEvent.pos_2()) + std::min(mutationEvent.pos_1(),mutationEvent.pos_2()); //length of the segment covering the origin
    return std::min(intern_mut_size, extern_mut_size);
  } else {
    if (mutationEvent.type() == MutationEventType::SMALL_INSERTION
        || mutationEvent.type() == MutationEventType::SMALL_DELETION) {
      return mutationEvent.number();
    } else {
      return 0;
    }
  }
}

void out::record_mutation_event(MutationEvent &mutationEvent,
                                unsigned int number_generation,
                                unsigned int size) {
  mutation_file << number_generation << ",";
  mutation_file << 0 << ",";  // GeneticUnit
  mutation_file << mutationEvent.type() << ",";
  mutation_file << mutationEvent.pos_1() << ",";
  mutation_file << mutation_event_pos_2(mutationEvent) << ",";
  mutation_file << mutation_event_pos_3(mutationEvent) << ",";
  mutation_file << mutation_event_pos_4(mutationEvent) << ",";
  mutation_file << mutation_event_invert(mutationEvent) << ",";
  mutation_file << 0 << ","; // align_score
  mutation_file << 0 << ","; // align_score2
  mutation_file << mutation_event_length_mutation(mutationEvent, size) << ",";
  mutation_file << -1 << ","; // rep_seg_length
  mutation_file << size << ",";
  mutation_file << 0 << ","; // Impact of the mutation on the metabolic error
  mutation_file << 0 << ","; // Nb of coding RNAs possibly disrupted by the breakpoints
  mutation_file << 0 << ","; // Nb of coding RNAs included in the segment
  mutation_file << 0 << std::endl; // Nb of coding RNAs that were completely included in the replaced segment
}
