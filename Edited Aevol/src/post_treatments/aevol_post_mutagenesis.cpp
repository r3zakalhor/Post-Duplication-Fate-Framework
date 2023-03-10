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




// =================================================================
//                              Includes
// =================================================================
#include <cerrno>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cinttypes>
#include <cmath>
#include <cassert>
#include <getopt.h>
#include <sys/stat.h>

#include <list>

#include <zlib.h>

// =================================================================
//                            Project Files
// =================================================================
#include "aevol.h"


using std::list;
using namespace aevol;

// =================================================================
//                     Command line option variables
// =================================================================
int32_t wanted_rank = -1;
int32_t wanted_index = -1;
int64_t num_gener = 0;
int32_t mutation_type = 0;
int32_t nb_mutants = -1;

// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);
void interpret_cmd_line_options(int argc, char* argv[]);

int main(int argc, char* argv[]) {
  interpret_cmd_line_options(argc, argv);

  // ------------------------------------------------------
  //  Load the backup and get the individual to be mutated
  // ------------------------------------------------------

  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(num_gener, true, false);

  if (exp_manager->output_m()->record_tree() == false) {
    // The following instruction is needed to ensure that methods
    // like ae_dna::do_deletion, ae_dna::do_inversion, etc
    // will create ae_mutation objects (otherwise they return NULL)
    exp_manager->output_m()->init_tree(exp_manager, 100);
  }


  if ((wanted_rank == -1) && (wanted_index == -1)) {
    wanted_rank = exp_manager->nb_indivs();  // the best one has rank N
  }

  // TODO: factor with duplicated code in robustness.cpp
  Individual* initial_indiv = nullptr;

  { // (local scope for variable `indivs` used as a shorthand)
    bool found = false;
    int32_t current_rank = -1;
    int32_t current_index = -1;
    list<Individual*> indivs = exp_manager->indivs();
    for (auto indiv = indivs.rbegin();
         not found and indiv != indivs.rend(); ++indiv) {
      current_index = (*indiv)->id();
      current_rank = (*indiv)->rank();

      if (wanted_index != -1 and current_index == wanted_index) {
        found = true;
        initial_indiv = (*indiv);
        wanted_rank = current_rank;
      }
      else if (current_rank == wanted_rank) {
        // no index was specified, we use the desired rank
        found = true;
        initial_indiv = (*indiv);
        wanted_index = current_index;
      }
    }

    if (not found) {
      Utils::ExitWithUsrMsg("sorry, the individual you have requested has not "
                            "been found");
    }
  }

  initial_indiv->Evaluate();
  initial_indiv->compute_statistical_data();
  initial_indiv->compute_non_coding();



  // ---------------------
  //  Prepare the output
  // ---------------------


  char mutation_type_name[24];
  switch (mutation_type) {
    case SWITCH: {
      snprintf(mutation_type_name, 23, "point-mutation");
      break;
    }
    case S_INS: {
      snprintf(mutation_type_name, 23, "small-insertion");
      break;
    }
    case S_DEL: {
      snprintf(mutation_type_name, 23, "small-deletion");
      break;
    }
    case DUPL: {
      snprintf(mutation_type_name, 23, "duplication");
      break;
    }
    case DEL: {
      snprintf(mutation_type_name, 23, "large-deletion");
      break;
    }
    case TRANS: {
      snprintf(mutation_type_name, 23, "translocation");
      break;
    }
    case INV: {
      snprintf(mutation_type_name, 23, "inversion");
      break;
    }
    default: {
      fprintf(stderr, "Error, unexpected mutation type.\n");
      exit(EXIT_FAILURE);
    }
  }

  char directory_name[64];
  snprintf(directory_name, 63, "analysis-generation_" TIMESTEP_FORMAT,
           num_gener);

  // Check whether the directory already exists and is writable
  if (access(directory_name, F_OK) == 0) {
    if (access(directory_name, X_OK | W_OK) != 0) {
      fprintf(stderr, "Error: cannot enter or write in directory %s.\n",
              directory_name);
      exit(EXIT_FAILURE);
    }
  }
  else {
    // Create the directory with permissions : rwx r-x r-x
    if (mkdir(directory_name,
              S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) != 0) {
      fprintf(stderr, "Error: cannot create directory %s.\n", directory_name);
      exit(EXIT_FAILURE);
    }
  }


  char output_file_name[256];
  snprintf(output_file_name, 255,
           "%s/mutagenesis-t" TIMESTEP_FORMAT "-i%" PRId32 "-r%" PRId32 "-%s.out", \
            directory_name, num_gener, wanted_index, wanted_rank,
           mutation_type_name);

  FILE* output = fopen(output_file_name, "w");
  if (output == NULL) {
    fprintf(stderr, "ERROR : Could not create the output file %s\n",
            output_file_name);
    exit(EXIT_FAILURE);
  }


  // Write the header

  int16_t col = 1;

  fprintf(output,
          "# ####################################################################################################\n");
  fprintf(output,
          "#   Single %s mutants of individual %" PRId32 " (rank %" PRId32 ") at generation %" PRId64 "\n",
          mutation_type_name, wanted_index, wanted_rank, num_gener);
  fprintf(output,
          "# ####################################################################################################\n");
  fprintf(output,
          "#  %" PRId16 ".  Mutation type    (0: switch, 1: smallins, 2: smalldel, 3:dupl, 4: del, 5:trans, 6:inv) \n",
          col);
  col++;
  fprintf(output,
          "#  %" PRId16 ".  Genetic unit which underwent the mutation (0 = chromosome) \n",
          col);
  col++;
  fprintf(output,
          "#  %" PRId16 ".  Length of this genetic unit before the event \n",
          col);
  col++;

  switch (mutation_type) {
    case SWITCH: {
      // Even though not all five columns are relevant for point mutations, we still write them all
      // to make the statistical analysis easier if other types of mutants are to be generated.
      // This way, for a given experiment, the number of columns will be the same for all types of mutants.
      fprintf(output,
              "#  %" PRId16 ".  pos0            (position of the point mutation on the genetic unit) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  pos1            (irrelevant for point mutations) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  pos2            (irrelevant for point mutations) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  pos3            (irrelevant for point mutations) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  invert          (irrelevant for point mutations) \n",
              col);
      col++;
      break;
    }
    case S_INS: {
      fprintf(output,
              "#  %" PRId16 ".  pos0            (position of the small insertion on the genetic unit) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  pos1            (irrelevant for small insertions) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  pos2            (irrelevant for small insertions) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  pos3            (irrelevant for small insertions) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  invert          (irrelevant for small insertions) \n",
              col);
      col++;
      break;
    }
    case S_DEL: {
      fprintf(output,
              "#  %" PRId16 ".  pos0            (position of the small deletion on the genetic unit) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  pos1            (irrelevant for small deletions) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  pos2            (irrelevant for small deletions) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  pos3            (irrelevant for small deletions) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  invert          (irrelevant for small deletions) \n",
              col);
      col++;
      break;
    }
    case DUPL: {
      fprintf(output,
              "#  %" PRId16 ".  pos0            (begin of the duplicated segment) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  pos1            (end of the duplicated segment) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  pos2            (reinsertion point of the duplicate in the genetic unit) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  pos3            (irrelevant for duplications) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  invert          (irrelevant for duplications) \n",
              col);
      col++;
      break;
    }
    case DEL: {
      fprintf(output,
              "#  %" PRId16 ".  pos0            (begin of the deleted segment) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  pos1            (end of the deleted segment) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  pos2            (irrelevant for large deletions) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  pos3            (irrelevant for large deletions) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  invert          (irrelevant for large deletions) \n",
              col);
      col++;
      break;
    }
    case TRANS: {
      fprintf(output,
              "#  %" PRId16 ".  pos0            (begin of the excised segment) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  pos1            (end of the excised segment) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  pos2            (cutting point in the excised segment when reinserted) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  pos3            (reinsertion point of the segment in the genetic unit) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  invert          (was the segment inverted when reinserted (0: no, 1: yes)) \n",
              col);
      col++;
      break;
    }
    case INV: {
      fprintf(output,
              "#  %" PRId16 ".  pos0            (begin of the inverted segment) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  pos1            (end of the inverted segment) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  pos2            (irrelevant for inversions) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  pos3            (irrelevant for inversions) \n",
              col);
      col++;
      fprintf(output,
              "#  %" PRId16 ".  invert          (irrelevant for inversions) \n",
              col);
      col++;
      break;
    }
    default: {
      fprintf(stderr, "Error: unexpected mutation type.\n");
      exit(EXIT_FAILURE);
    }
  }


  if (initial_indiv->with_alignments()) {
    fprintf(output,
            "#  %" PRId16 ".  align_score1    (score that was needed for the rearrangement to occur)\n",
            col);
    col++;
    fprintf(output,
            "#  %" PRId16 ".  align_score2    (score for the reinsertion for translocations)\n",
            col);
    col++;
  }


  fprintf(output,
          "#  %" PRId16 ".  Length of the {inserted, deleted, duplicated, translocated, inverted} segment \n",
          col);
  col++;
  fprintf(output,
          "#  %" PRId16 ".  Number of coding RNAs possibly disrupted by the breakpoint(s) \n",
          col);
  col++;
  fprintf(output,
          "#  %" PRId16 ".  Number of coding RNAs completely included in the segment \n",
          col);
  col++;
  fprintf(output, "#  %" PRId16 ".  Metabolic error after the mutation \n",
          col);
  col++;

  if (exp_manager->with_secretion()) {
    fprintf(output, "#  %" PRId16 ".  Secretion error after the mutation \n",
            col);
    col++;
  }

  if ((exp_manager->with_plasmids()) &&
      (exp_manager->tune_donor_ability() != 0.0)) {
    fprintf(output,
            "#  %" PRId16 ".  Error on the donor ability after the mutation \n",
            col);
    col++;
  }

  if ((exp_manager->with_plasmids()) &&
      (exp_manager->tune_recipient_ability() != 0.0)) {
    fprintf(output,
            "#  %" PRId16 ".  Error on the recipient ability after the mutation \n",
            col);
    col++;
  }

  fprintf(output, "#  %" PRId16 ".  Total genome size after the mutation \n",
          col);
  col++;
  fprintf(output,
          "#  %" PRId16 ".  Number of coding RNAs after the mutation \n", col);
  col++;
  fprintf(output,
          "#  %" PRId16 ".  Number of non coding RNAs after the mutation \n",
          col);
  col++;
  fprintf(output,
          "#  %" PRId16 ".  Number of functional coding sequences after the mutation \n",
          col);
  col++;
  fprintf(output,
          "#  %" PRId16 ".  Number of non functional coding sequences after the mutation \n",
          col);
  col++;
  fprintf(output,
          "#  %" PRId16 ".  Number of coding bases after the mutation \n", col);
  col++;
  fprintf(output,
          "#  %" PRId16 ".  Number of transcribed but not translated bases after the mutation\n",
          col);
  col++;
  fprintf(output,
          "#  %" PRId16 ".  Number of non transcribed bases after the mutation \n",
          col);
  col++;
  fprintf(output,
          "#  %" PRId16 ".  Number of bases belonging to at least one coding RNA, after the mutation \n",
          col);
  col++;
  fprintf(output,
          "#  %" PRId16 ".  Number of bases not belonging to any coding RNA, after the mutation \n",
          col);
  col++;
  fprintf(output,
          "####################################################################################################################\n");
  fprintf(output,
          "#  Values for the initial individual [irr = irrelevant]: \n");
  fprintf(output,
          "####################################################################################################################\n");
  fprintf(output, "#  ");
  fprintf(output, "irr ");
  fprintf(output, "irr "); // genetic unit number (0 for the chromosome)
  fprintf(output, "irr ");
  fprintf(output, "irr ");
  fprintf(output, "irr ");
  fprintf(output, "irr ");
  fprintf(output, "irr ");
  fprintf(output, "irr ");
  if (initial_indiv->with_alignments()) {
    fprintf(output, "irr ");
    fprintf(output, "irr ");
  }
  fprintf(output, "irr ");
  fprintf(output, "irr ");
  fprintf(output, "irr ");
  fprintf(output, "%e ",
          initial_indiv->dist_to_target_by_feature(METABOLISM));
  if (exp_manager->with_secretion()) {
    fprintf(output, "%e ",
            initial_indiv->dist_to_target_by_feature(SECRETION));
  }
  if ((exp_manager->with_plasmids()) &&
      (exp_manager->tune_donor_ability() != 0.0)) {
    fprintf(output, "%e ", initial_indiv->dist_to_target_by_feature(DONOR));
  }
  if ((exp_manager->with_plasmids()) &&
      (exp_manager->tune_recipient_ability() != 0.0)) {
    fprintf(output, "%e ",
            initial_indiv->dist_to_target_by_feature(RECIPIENT));
  }
  fprintf(output, "%" PRId32 " ", initial_indiv->total_genome_size());
  fprintf(output, "%" PRId32 " ", initial_indiv->nb_coding_RNAs());
  fprintf(output, "%" PRId32 " ", initial_indiv->nb_non_coding_RNAs());
  fprintf(output, "%" PRId32 " ", initial_indiv->nb_functional_genes());
  fprintf(output, "%" PRId32 " ", initial_indiv->nb_non_functional_genes());
  fprintf(output, "%" PRId32 " ", initial_indiv->total_genome_size() -
                                  initial_indiv->nb_bases_in_0_CDS()); // coding bp
  fprintf(output, "%" PRId32 " ", initial_indiv->total_genome_size() -
                                  initial_indiv->nb_bases_in_0_RNA() -
                                  (initial_indiv->total_genome_size() -
                                   initial_indiv->nb_bases_in_0_CDS())); // transcribed but not translated bp
  fprintf(output, "%" PRId32 " ",
          initial_indiv->nb_bases_in_0_RNA()); // not transcribed bp
  fprintf(output, "%" PRId32 " ", initial_indiv->total_genome_size() -
                                  initial_indiv->nb_bases_in_0_coding_RNA());
  fprintf(output, "%" PRId32 " ",
          initial_indiv->nb_bases_in_0_coding_RNA());
  fprintf(output, "\n");
  fprintf(output,
          "####################################################################################################################\n");

  // ---------------------------------------
  //  Create the mutants and evaluate them
  // ---------------------------------------

  Individual* mutant = NULL;
  Mutation* mut = NULL;
  int32_t nb_genetic_units = initial_indiv->nb_genetic_units();
  double* relative_lengths_genetic_units = NULL;
  int32_t u = 0;
  double alea, cumul;
  int32_t pos, pos0, pos1, pos2, pos3;
  int32_t mut_length;
  int16_t align_score1, align_score2;
  bool invert;
  VisAVis* alignment_1 = NULL;
  VisAVis* alignment_2 = NULL;
  int32_t nb_pairs;
  Dna* initial_dna = NULL;
  int32_t initial_len;
  bool rear_done;
  int32_t nb_genes_at_breakpoints;
  int32_t nb_genes_in_segment;
  int32_t nb_genes_in_replaced_segment;
  double metabolic_error_after = -1.0, secretion_error_after = -1.0;

  if (mutation_type == SWITCH) {
    // *********************  Exhaustive mutagenesis  *************************
    pos0 = pos1 = pos2 = pos3 = -1;
    mut_length = -1;
    align_score1 = align_score2 = -1;
    invert = false;
    metabolic_error_after = -1.0;
    secretion_error_after = -1.0;

    for (const auto& gu: initial_indiv->genetic_unit_list()) {
      initial_len = gu.dna()->length();

      for (pos = 0; pos < initial_len; pos++) {
        mutant = new Individual(*initial_indiv);
        
        #ifdef BASE_2
        mutant->genetic_unit(u).dna()->do_switch(pos);
        mut = new PointMutation(pos);
        #elif BASE_4
        char base = mutant->genetic_unit(u).dna()->data()[pos];
        base = get_complementary_base(base);


        mutant->genetic_unit(u).dna()->do_switch(pos,base);
        mut = new PointMutation(pos,base);
        #endif

        mut_length = 1;
        pos0 = pos;

        initial_indiv->genetic_unit_nonconst(
            u).compute_nb_of_affected_genes(mut, nb_genes_at_breakpoints,
                                            nb_genes_in_segment,
                                            nb_genes_in_replaced_segment);

        // Evaluate the mutant, compute its statistics
        mutant->ReevaluateInContext(initial_indiv->habitat());
        mutant->compute_statistical_data();
        mutant->compute_non_coding();

        metabolic_error_after = mutant->dist_to_target_by_feature(
            METABOLISM);
        if (exp_manager->with_secretion()) {
          secretion_error_after = mutant->dist_to_target_by_feature(
              SECRETION);
        }

        // Write the description of the mutant in the output file
        fprintf(output, "%" PRId32 " ", mutation_type);
        fprintf(output, "%" PRId32 " ",
                u); // genetic unit number (0 for the chromosome)
        fprintf(output, "%" PRId32 " ", initial_indiv->genetic_unit(
            u).dna()->length()); // Length of GU before the event
        fprintf(output, "%" PRId32 " ", pos0);
        fprintf(output, "%" PRId32 " ", pos1);
        fprintf(output, "%" PRId32 " ", pos2);
        fprintf(output, "%" PRId32 " ", pos3);
        if (invert) fprintf(output, "1 "); else fprintf(output, "0 ");
        if (mutant->with_alignments()) {
          fprintf(output, "%" PRId16 " ", align_score1);
          fprintf(output, "%" PRId16 " ", align_score2);
        }
        fprintf(output, "%" PRId32 " ", mut_length);
        fprintf(output, "%" PRId32 " ", nb_genes_at_breakpoints);
        fprintf(output, "%" PRId32 " ", nb_genes_in_segment);
        fprintf(output, "%e ", metabolic_error_after);
        if (exp_manager->with_secretion()) {
          fprintf(output, "%e ", secretion_error_after);
        }
        if ((exp_manager->with_plasmids()) &&
            (exp_manager->tune_donor_ability() != 0.0)) {
          fprintf(output, "%e ", mutant->dist_to_target_by_feature(DONOR));
        }
        if ((exp_manager->with_plasmids()) &&
            (exp_manager->tune_recipient_ability() != 0.0)) {
          fprintf(output, "%e ",
                  mutant->dist_to_target_by_feature(RECIPIENT));
        }
        fprintf(output, "%" PRId32 " ", mutant->total_genome_size());
        fprintf(output, "%" PRId32 " ", mutant->nb_coding_RNAs());
        fprintf(output, "%" PRId32 " ", mutant->nb_non_coding_RNAs());
        fprintf(output, "%" PRId32 " ", mutant->nb_functional_genes());
        fprintf(output, "%" PRId32 " ", mutant->nb_non_functional_genes());
        fprintf(output, "%" PRId32 " ", mutant->total_genome_size() -
                                        mutant->nb_bases_in_0_CDS()); // coding bp
        fprintf(output, "%" PRId32 " ", mutant->total_genome_size() -
                                        mutant->nb_bases_in_0_RNA() -
                                        (mutant->total_genome_size() -
                                         mutant->nb_bases_in_0_CDS())); // transcribed but not translated bp
        fprintf(output, "%" PRId32 " ",
                mutant->nb_bases_in_0_RNA()); // not transcribed bp
        fprintf(output, "%" PRId32 " ", mutant->total_genome_size() -
                                        mutant->nb_bases_in_0_coding_RNA());
        fprintf(output, "%" PRId32 " ", mutant->nb_bases_in_0_coding_RNA());
        fprintf(output, "\n");

        delete mutant;
        delete mut;
      }

      u++;
    }

  }
  else {
    // *******************************  Sampling nb_mutants mutants  **********************************

    relative_lengths_genetic_units = new double[nb_genetic_units];

    for (const auto& gu: initial_indiv->genetic_unit_list())
      relative_lengths_genetic_units[u++] =
          gu.dna()->length() /
          static_cast<double>(initial_indiv->total_genome_size());

    for (int32_t i = 0; i < nb_mutants; i++) {
      mutant = new Individual(*initial_indiv);

      // Pick the genetic unit which will undergo the mutation
      alea = mutant->mut_prng()->random();
      u = 0;
      cumul = relative_lengths_genetic_units[0];
      while (alea > cumul) {
        u++;
        cumul += relative_lengths_genetic_units[u];
      }


      // Ask the genetic unit to perform the mutation and store it
      pos0 = pos1 = pos2 = pos3 = -1;
      mut_length = -1;
      align_score1 = align_score2 = -1;
      invert = false;

      alignment_1 = NULL;
      alignment_2 = NULL;
      initial_dna = initial_indiv->genetic_unit(u).dna();
      initial_len = initial_dna->length();
      metabolic_error_after = -1.0;
      secretion_error_after = -1.0;

      switch (mutation_type) {
        // Locally we need the precise type of the mutation. Outside of the
        // switch, we will need the generic 'mut'. We could have static_casted
        // 'mut' here but it felt better declaring a specific variable
        // instead...
        case S_INS: {
          // cf. comment at top of switch statement
          SmallInsertion* small_ins;
          do {
            mut = small_ins = mutant->genetic_unit(u).dna()->
                do_small_insertion();
          } while (mut == NULL);
          pos0 = small_ins->pos();
          mut_length = small_ins->length();
          break;
        }
        case S_DEL: {
          // cf. comment at top of switch statement
          SmallDeletion* small_del;
          do {
            mut = small_del = mutant->genetic_unit(u).dna()->
                do_small_deletion();
          } while (mut == NULL);
          pos0 = small_del->pos();
          mut_length = small_del->length();
          break;
        }
        case DUPL: {
          // cf. comment at top of switch statement
          Duplication* duplication;
          if (mutant->with_alignments()) {
            // TODO(dpa) Encapsulate in method do_duplication_align()
            rear_done = false;
            do {
              do {
                nb_pairs = initial_len;
                alignment_1 = initial_dna->search_alignment(initial_dna,
                                                            nb_pairs, DIRECT);
              } while (alignment_1 == NULL);
              mut_length = Utils::mod(
                  alignment_1->i_2() - alignment_1->i_1(), initial_len);
              rear_done = mutant->genetic_unit(u).dna()->do_duplication(
                  alignment_1->i_1(), alignment_1->i_2(),
                  alignment_1->i_2());
            } while (!rear_done);

            mut = duplication = new Duplication(alignment_1->i_1(),
                                                alignment_1->i_2(),
                                                alignment_1->i_2(),
                                                mut_length,
                                                alignment_1->score());
          }
          else {
            do {
              mut = duplication = mutant->genetic_unit(u).dna()->
                  do_duplication();
            } while (mut == NULL);
          }

          pos0 = duplication->pos1();
          pos1 = duplication->pos2();
          pos2 = duplication->pos3();
          align_score1 = duplication->align_score();
          mut_length = duplication->length();
          break;
        }
        case DEL: {
          // cf. comment at top of switch statement
          Deletion* deletion;
          if (mutant->with_alignments()) {
            rear_done = false;
            do {
              do {
                nb_pairs = initial_len;
                alignment_1 = initial_dna->search_alignment(initial_dna,
                                                            nb_pairs, DIRECT);
              } while (alignment_1 == NULL);
              mut_length = Utils::mod(
                  alignment_1->i_2() - alignment_1->i_1(), initial_len);
              rear_done = mutant->genetic_unit(u).dna()->do_deletion(
                  alignment_1->i_1(), alignment_1->i_2());
            } while (!rear_done);

            mut = deletion = new Deletion(alignment_1->i_1(),
                                          alignment_1->i_2(),
                                          mut_length,
                                          alignment_1->score());
          }
          else {
            do {
              mut = deletion = mutant->genetic_unit(u).dna()->
                  do_deletion();
            } while (mut == NULL);
          }

          pos0 = deletion->pos1();
          pos1 = deletion->pos2();
          align_score1 = deletion->align_score();
          mut_length = deletion->length();
          break;
        }
        case TRANS: {
          // cf. comment at top of switch statement
          Translocation* translocation;
          // TO DO: problems might arise because the ae_mutation does not
          //        record whether it was an intra- or interGU translocation

          if (mutant->with_alignments()) {
            // TODO(dpa) Encapsulate in method do_duplication_align()
            rear_done = false;
            do {
              do {
                nb_pairs = initial_len;
                alignment_1 = initial_dna->search_alignment(initial_dna,
                                                            nb_pairs, DIRECT);
              } while (alignment_1 == NULL);
              // Make sure the segment to be translocated doesn't contain OriC
              // TODO(dpa) is that still necessary?
              if (alignment_1->i_1() > alignment_1->i_2()) {
                alignment_1->swap();
              }
              mut_length = Utils::mod(
                  alignment_1->i_2() - alignment_1->i_1(), initial_len);

              // Extract the segment to be translocated
              GeneticUnit* tmp_segment = mutant->genetic_unit(
                  u).dna()->extract_into_new_GU(alignment_1->i_1(),
                                                    alignment_1->i_2());
              // Look for a "new" alignment between this segment and the
              // remaining of the chromosome
              do {
                nb_pairs = initial_len;
                alignment_2 = tmp_segment->dna()->search_alignment(
                    mutant->genetic_unit(u).dna(), nb_pairs,
                    BOTH_SENSES);
              } while (alignment_2 == NULL);
              invert = (alignment_2->sense() == INDIRECT);
              // Reinsert the segment into the genetic unit
              mutant->genetic_unit(u).dna()->
                  insert_GU(tmp_segment, alignment_2->i_2(),
                            alignment_2->i_1(), invert);
              rear_done = true;
              delete tmp_segment;
            } while (!rear_done);

            mut = translocation = new Translocation(alignment_1->i_1(),
                                                    alignment_1->i_2(),
                                                    alignment_2->i_1(),
                                                    alignment_2->i_2(),
                                                    mut_length, invert,
                                                    alignment_1->score(),
                                                    alignment_2->score());
          }
          else {
            do {
              mut = translocation = mutant->genetic_unit(u).dna()->
                  do_translocation();
            }
            while (mut == NULL);
          }

          pos0 = translocation->pos1();
          pos1 = translocation->pos2();
          pos2 = translocation->pos3();
          pos3 = translocation->pos4();
          align_score1 = translocation->align_score_1();
          align_score2 = translocation->align_score_2();
          invert = translocation->invert();
          mut_length = translocation->length();
          break;
        }
        case INV: {
          // cf. comment at top of switch statement
          Inversion* inversion;
          if (mutant->with_alignments()) {
            // TODO(dpa) Encapsulate in method do_duplication_align()
            rear_done = false;
            do {
              do {
                nb_pairs = initial_len;
                alignment_1 = initial_dna->search_alignment(initial_dna,
                                                            nb_pairs, INDIRECT);
              } while (alignment_1 == NULL);
              // Make sure the segment to be inverted doesn't contain OriC
              if (alignment_1->i_1() > alignment_1->i_2()) {
                alignment_1->swap();
              }
              mut_length = Utils::mod(
                  alignment_1->i_2() - alignment_1->i_1(), initial_len);
              rear_done = mutant->genetic_unit(u).dna()->do_inversion(
                  alignment_1->i_1(), alignment_1->i_2());
            } while (!rear_done);
            mut = inversion = new Inversion(alignment_1->i_1(),
                                            alignment_1->i_2(),
                                            mut_length,
                                            alignment_1->score());
          }
          else {
            do {
              mut = inversion = mutant->genetic_unit(u).dna()->
                  do_inversion();
            }
            while (mut == NULL);
          }

          pos0 = inversion->pos1();
          pos1 = inversion->pos2();
          align_score1 = inversion->align_score();
          mut_length = inversion->length();
          break;
        }
        default: {
          fprintf(stderr, "Error, unexpected mutation type\n");
          break;
        }
      }

      // TO DO: improve this method to make it work also with
      // interGU translocations
      initial_indiv->genetic_unit_nonconst(u).compute_nb_of_affected_genes(
          mut, nb_genes_at_breakpoints, nb_genes_in_segment,
          nb_genes_in_replaced_segment);


      // Evaluate the mutant, compute its statistics
      mutant->ReevaluateInContext(initial_indiv->habitat()); //mutant->Reevaluate();
      mutant->compute_statistical_data();
      mutant->compute_non_coding();

      metabolic_error_after = mutant->dist_to_target_by_feature(METABOLISM);
      if (exp_manager->with_secretion()) {
        secretion_error_after = mutant->dist_to_target_by_feature(
            SECRETION);
      }

      // Write the description of the mutant in the output file
      fprintf(output, "%" PRId32 " ", mutation_type);
      fprintf(output, "%" PRId32 " ",
              u); // genetic unit number (0 for the chromosome)
      fprintf(output, "%" PRId32 " ", initial_indiv->genetic_unit(
          u).dna()->length()); // Length of GU before the event
      fprintf(output, "%" PRId32 " ", pos0);
      fprintf(output, "%" PRId32 " ", pos1);
      fprintf(output, "%" PRId32 " ", pos2);
      fprintf(output, "%" PRId32 " ", pos3);
      if (invert) fprintf(output, "1 "); else fprintf(output, "0 ");
      if (mutant->with_alignments()) {
        fprintf(output, "%" PRId16 " ", align_score1);
        fprintf(output, "%" PRId16 " ", align_score2);
      }
      fprintf(output, "%" PRId32 " ", mut_length);
      fprintf(output, "%" PRId32 " ", nb_genes_at_breakpoints);
      fprintf(output, "%" PRId32 " ", nb_genes_in_segment);
      fprintf(output, "%e ", metabolic_error_after);
      if (exp_manager->with_secretion()) {
        fprintf(output, "%e ", secretion_error_after);
      }
      if ((exp_manager->with_plasmids()) &&
          (exp_manager->tune_donor_ability() != 0.0)) {
        fprintf(output, "%e ", mutant->dist_to_target_by_feature(DONOR));
      }
      if ((exp_manager->with_plasmids()) &&
          (exp_manager->tune_recipient_ability() != 0.0)) {
        fprintf(output, "%e ",
                mutant->dist_to_target_by_feature(RECIPIENT));
      }
      fprintf(output, "%" PRId32 " ", mutant->total_genome_size());
      fprintf(output, "%" PRId32 " ", mutant->nb_coding_RNAs());
      fprintf(output, "%" PRId32 " ", mutant->nb_non_coding_RNAs());
      fprintf(output, "%" PRId32 " ", mutant->nb_functional_genes());
      fprintf(output, "%" PRId32 " ", mutant->nb_non_functional_genes());
      fprintf(output, "%" PRId32 " ", mutant->total_genome_size() -
                                      mutant->nb_bases_in_0_CDS()); // coding bp
      fprintf(output, "%" PRId32 " ", mutant->total_genome_size() -
                                      mutant->nb_bases_in_0_RNA() -
                                      (mutant->total_genome_size() -
                                       mutant->nb_bases_in_0_CDS())); // transcribed but not translated bp
      fprintf(output, "%" PRId32 " ",
              mutant->nb_bases_in_0_RNA()); // not transcribed bp
      fprintf(output, "%" PRId32 " ", mutant->total_genome_size() -
                                      mutant->nb_bases_in_0_coding_RNA());
      fprintf(output, "%" PRId32 " ", mutant->nb_bases_in_0_coding_RNA());
      fprintf(output, "\n");


      delete mutant;
      delete mut;
    }

    delete[] relative_lengths_genetic_units;
  }


  delete exp_manager;

  return EXIT_SUCCESS;
}

void interpret_cmd_line_options(int argc, char* argv[]) {
  const char* options_list = "hVg:r:i:m:n:";
  static struct option long_options_list[] = {
      {"help",          no_argument,       NULL, 'h'},
      {"version",       no_argument,       NULL, 'V'},
      {"gener",         required_argument, NULL, 'g'},
      {"rank",          required_argument, NULL, 'r'},
      {"index",         required_argument, NULL, 'i'},
      {"mutation-type", required_argument, NULL, 'm'},
      {"nb-mutants",    required_argument, NULL, 'n'},
      {0, 0, 0,                                  0}
  };

  int option = -1;
  bool rank_already_set = false;
  bool index_already_set = false;
  while ((option = getopt_long(argc, argv, options_list, long_options_list,
                               NULL)) != -1) {
    switch (option) {
      case 'h':
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      case 'V':
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      case 'g':
        if (strcmp(optarg, "") == 0) {
          fprintf(stderr,
                  "%s: error: Option -g or --gener : missing argument.\n",
                  argv[0]);
          exit(EXIT_FAILURE);
        }
        num_gener = atol(optarg);
        break;
      case 'r':
        if (index_already_set) {
          fprintf(stderr,
                  "%s: error: Options -r and -i are incompatible. Please choose one of them only.\n",
                  argv[0]);
          exit(EXIT_FAILURE);
        }
        wanted_rank = atol(optarg);
        rank_already_set = true;
        break;
      case 'i':
        if (rank_already_set) {
          fprintf(stderr,
                  "%s: error: Options -r and -i are incompatible. Please choose one of them only.\n",
                  argv[0]);
          fprintf(stderr, "           Use %s --help for more information.\n",
                  argv[0]);
          exit(EXIT_FAILURE);
        }
        wanted_index = atol(optarg);
        index_already_set = true;
        break;
      case 'm':
        mutation_type = (MutationType) atol(optarg);
        if (mutation_type == SWITCH) {
        }
        else if ((mutation_type == S_INS) || (mutation_type == S_DEL) ||
                 (mutation_type == DUPL) || (mutation_type == DEL) ||
                 (mutation_type == TRANS) || (mutation_type == INV)) {
        }
        else {
          fprintf(stderr,
                  "%s: error: So far, mutagenesis is implemented only for "
                      "point mutations, small insertions, \n"
                      "           small deletions, duplications, deletions, "
                      "translocations or inversions.\n"
                      "           It is not available yet for lateral transfer.\n",
                  argv[0]);
          exit(EXIT_FAILURE);
        }
        break;
      case 'n':
        nb_mutants = atol(optarg);
        if (nb_mutants <= 0) {
          fprintf(stderr,
                  "%s: error: The number of mutants (option -n) must be "
                      "positive.\n",
                  argv[0]);
          exit(EXIT_FAILURE);
        }
        break;
    }
  }


  if ((mutation_type == SWITCH) && (nb_mutants != -1)) {
    printf("For point mutations, the mutagenesis will be exhaustive, "
               "the number of \n"
               "mutants will be ignored.\n");
  }
  else if ((mutation_type != SWITCH) && (nb_mutants == -1)) {
    nb_mutants = 1000;
    printf("Mutagenesis cannot be exhaustive in a reasonable time for "
               "mutations \n");
    printf("other than point mutations. A sample of %" PRId32 " mutants will "
               "be generated.\n",
           nb_mutants);
  }
}

/*!
  \brief

*/
void print_help(char* prog_path) {
  // Get the program file-name in prog_name (strip prog_path of the path)
  char* prog_name; // No new, it will point to somewhere inside prog_path
  if ((prog_name = strrchr(prog_path, '/'))) {
    prog_name++;
  }
  else { prog_name = prog_path; }

  printf("\n");
  printf(
      "*********************** aevol - Artificial Evolution ******************* \n");
  printf(
      "*                                                                      * \n");
  printf(
      "*                    Mutagenesis post-treatment program                * \n");
  printf(
      "*                                                                      * \n");
  printf(
      "************************************************************************ \n");
  printf("\n\n");
  printf("This program is Free Software. No Warranty.\n");
  printf("\n");
  printf("Usage : %s -h\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf(
      "   or : %s -g NUMGENER [-r RANK | -i INDEX] [-m MUTATIONTYPE] [-n NBMUTANTS]\n",
      prog_name);
  printf("\n");
  printf(
      "This program creates and evaluates single mutants of an individual saved in a backup, \n");
  printf(
      "by default the best of its generation. Use either the -r or the -i option to\n");
  printf(
      "select another individual than the best one: with -i, you have to provide the\n");
  printf(
      "ID of the individual, and with -r the rank (1 for the individual with the lowest\n");
  printf("fitness, N for the fittest one). \n\n");
  printf(
      "The type of mutations to perform must be specified with the -m option. \n");
  printf(
      "Choose 0 to create mutants with a point mutation, 1 for a small insertion, \n");
  printf(
      "2 for a small deletion, 3 for a duplication, 4 for a large deletion, \n");
  printf("5 for a translocation or 6 for an inversion. \n\n");
  printf(
      "For the point mutations, all single mutants will be created and evaluated. For the\n");
  printf(
      "other mutation types, an exhaustive mutagenesis would be too long, hence only a\n");
  printf(
      "sample of mutants (1000 by default) will be generated. Use option -n to specify\n");
  printf("another sample size.\n\n");
  printf(
      "The output file will be placed in a subdirectory called analysis-generationNUMGENER.\n");
  printf("\n");
  printf("\n");
  printf("\t-h or --help    : Display this help, then exit\n");
  printf("\n");
  printf("\t-V or --version : Print version number, then exit\n");
  printf("\n");
  printf("\t-g NUMGENER or --gener NUMGENER : \n");
  printf(
      "\t                  Generation of the backup containing the individual of interest\n");
  printf("\n");
  printf("\t-i INDEX or --index INDEX : \n");
  printf(
      "\t                  Index of individual of interest. Should be comprised between 0 and N-1, where\n");
  printf("\t                  N is the size of the population.\n");
  printf("\n");
  printf("\t-r RANK or --rank RANK : \n");
  printf(
      "\t                  Rank of individual of interest. Should be comprised between 1 and N, where\n");
  printf(
      "\t                  N is the size of the population. Default = N (fittest individual).\n");
  printf("\n");
  printf("\t-m MUTATIONTYPE or --mutation-type MUTATIONTYPE : \n");
  printf(
      "\t                  Integer type of the mutation carried by each mutant: 0 for a point mutation, 1 for a \n");
  printf(
      "\t                  small insertion, 2 for small deletions, 3 for a duplication, 4 for a large deletion, \n");
  printf("\t                  5 for a translocation or 6 for an inversion. \n");
  printf("\n");
  printf("\t-n NBMUTANTS or --nb-mutants NBMUTANTS : \n");
  printf(
      "\t                  Number of single mutants to create and evaluate. Default = 1000. Note that this option\n");
  printf(
      "\t                  is ignored in the case of point mutations, where all single mutants are created \n");
  printf("\t                  (exhaustive mutagenesis). \n");
  printf("\n");

  printf("\n");
}
