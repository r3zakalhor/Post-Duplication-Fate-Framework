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
#include <cinttypes>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cerrno>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

#include "AncestorStats.h"

namespace aevol {


//#############################################################################
//
//                                Class AncestorStats
//
//#############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================

AncestorStats::AncestorStats() {
  stat_files_ = new FILE*[NB_TYPES];
  stat_files_names_ = new char*[NB_TYPES];
  for(int8_t s_type = 0 ; s_type < NB_TYPES ; s_type++) {
    stat_files_[s_type] = nullptr;
    stat_files_names_[s_type] = nullptr;
  }
  exp_m_ = nullptr;
  phenotypicTH_ = nullptr;
  mystats_ = nullptr;
  indiv_ = nullptr;
  grid_cell_ = nullptr;
}

// =================================================================
//                             Destructors
// =================================================================

AncestorStats::~AncestorStats() {
  delete exp_m_;
  delete mystats_;
  if(grid_cell_ != nullptr) {
    delete &(grid_cell_->habitat().phenotypic_target_handler_nonconst());
    delete grid_cell_;
  }

  delete[] stat_files_;

  for (int8_t s_type = 0 ; s_type < NB_TYPES ; s_type++) {
    if(stat_files_names_[s_type] != nullptr) {
      delete[] stat_files_names_[s_type];
    }
  }
  delete[] stat_files_names_;
}

// =================================================================
//                            Public Methods
// =================================================================

void AncestorStats::setup_anc_stat(int64_t t_mrca) {
  exp_m_ = new ExpManager();
  exp_m_->load(-1, false, false);
  phenotypicTH_ = exp_m_->world()->phenotypic_target_handler();

  b_zones_ = (phenotypicTH_->phenotypic_target().nb_segments() > 1);

  Set_file_names();
  Open_files(t_mrca);
  auto prefix = "ancestor_stats/ancestor_stats";
  bool best_indiv_only = true;
  bool addition_old_stats = (t_mrca > 0);
  bool delete_old_stats = false;

  mystats_ = new Stats(exp_m_, t_mrca, best_indiv_only, prefix, "",
                             addition_old_stats, delete_old_stats);

  Flush();

  if(t_mrca > 0) {
    load_grid_cell(exp_m_);
  }
}

void AncestorStats::setup_anc_indiv(Individual* indiv) {
  indiv_ = indiv;
  indiv_->Evaluate();
  indiv_->compute_statistical_data();
  indiv_->compute_non_coding();

  mystats_->write_statistics_of_this_indiv(0, indiv_, nullptr);
  write_environment_stats(0, phenotypicTH_, stat_files_[ENVIRONMENT_STAT]);
  write_terminators_stats(0, indiv_, stat_files_[TERMINATOR_STAT]);
  if(b_zones_)
  {
    write_zones_stats(0, indiv_, phenotypicTH_, stat_files_[ZONES_STAT]);
  }
  write_operons_stats(0, indiv_, stat_files_[OPERONS_STAT]);

  Flush();
}

void AncestorStats::process_evolution(ReplicationReport* rep, int64_t t) {
  int32_t unitlen_before;
  double metabolic_error_before;
  double fitness_before;
  double impact_on_metabolic_error;
  double impact_on_fitness;
  char mut_descr_string[255];

  indiv_->Reevaluate();

  GeneticUnit& gen_unit = indiv_->genetic_unit_nonconst(0);

  const auto& dnarep = rep->dna_replic_report();

  dnarep.iter_muts([&](const auto& mut) {
    metabolic_error_before = indiv_->dist_to_target_by_feature(METABOLISM);
    fitness_before = indiv_->fitness();
    unitlen_before = gen_unit.dna()->length();

    gen_unit.dna()->undergo_this_mutation(*mut);

    indiv_->Reevaluate();
    impact_on_metabolic_error =
        indiv_->dist_to_target_by_feature(METABOLISM) -
        metabolic_error_before;
    impact_on_fitness = indiv_->fitness() - fitness_before;
    mut->generic_description_string(mut_descr_string);

    fprintf(stat_files_[FIXED_MUTATION_STAT],
            "%" PRId64 " %" PRId32 " %s %" PRId32 " %e %e \n",
            t, 0, mut_descr_string, unitlen_before,
            impact_on_metabolic_error, impact_on_fitness);
  });

  indiv_->Reevaluate();
  indiv_->compute_statistical_data();
  indiv_->compute_non_coding();

  mystats_->write_statistics_of_this_indiv(t, indiv_, rep);
  write_environment_stats(t, phenotypicTH_, stat_files_[ENVIRONMENT_STAT]);
  write_terminators_stats(t, indiv_, stat_files_[TERMINATOR_STAT]);
  if(b_zones_)
  {
    write_zones_stats(t, indiv_, phenotypicTH_, stat_files_[ZONES_STAT]);
  }
  write_operons_stats(t, indiv_, stat_files_[OPERONS_STAT]);
}

void AncestorStats::Flush() {
  mystats_->flush();
  for (int8_t s_type = 0 ; s_type < NB_TYPES ; s_type++) {
    if(stat_files_[s_type] != nullptr) {
      fflush(stat_files_[s_type]);
    }
  }
}

void AncestorStats::Close() {
  if(indiv_ != nullptr)
    save_grid_cell();

  for (int8_t s_type = 0 ; s_type < NB_TYPES ; s_type++) {
    if(stat_files_[s_type] != nullptr) {
      fclose(stat_files_[s_type]);
    }
  }
}

// =================================================================
//                           Protected Methods
// =================================================================

void AncestorStats::Set_file_names() {

  mkdir("stats/ancestor_stats/", 0755);

  const char* stat_type[NB_TYPES] =
        {"_envir", "_nb_term", "_zones", "_operons", "_fixed_mutations"};

  auto prefix = "ancestor_stats/ancestor_stats";
  for(int8_t s_type = 0 ; s_type < NB_TYPES ; s_type++) {
    if(s_type == ZONES_STAT && !b_zones_)
      continue;
    stat_files_names_[s_type] = new char[255];
    sprintf(stat_files_names_[s_type], STATS_DIR"/%s%s.out", prefix, stat_type[s_type]);
  }
}

void AncestorStats::Open_files(int64_t time) {
  if(time > 0) {
    Resume_stats(time);
  }
  else {
    for(int8_t s_type = 0 ; s_type < NB_TYPES ; s_type++) {
      if(stat_files_names_[s_type] != nullptr) {
        stat_files_[s_type] = fopen(stat_files_names_[s_type], "w");
      }
    }
    Write_headers();
  }
}

void AncestorStats::Write_headers() {
  if(stat_files_[ENVIRONMENT_STAT] != nullptr) {
    fprintf(stat_files_[ENVIRONMENT_STAT],
            "# Each line contains: Generation, and then, for each gaussian: M W H.\n");
    fprintf(stat_files_[ENVIRONMENT_STAT], "#\n");
  }

  if(stat_files_[TERMINATOR_STAT] != nullptr) {
    fprintf(stat_files_[TERMINATOR_STAT], "# Each line contains : \n");
    fprintf(stat_files_[TERMINATOR_STAT], "#   * Generation\n");
    fprintf(stat_files_[TERMINATOR_STAT], "#   * Genome size\n");
    fprintf(stat_files_[TERMINATOR_STAT], "#   * Terminator number\n");
    fprintf(stat_files_[TERMINATOR_STAT], "#\n");
  }

  if(stat_files_[ZONES_STAT] != nullptr) {
    fprintf(stat_files_[ZONES_STAT], "# Each line contains : Generation, and then, for each zone:\n");
    fprintf(stat_files_[ZONES_STAT], "#   * Number of activation genes\n");
    fprintf(stat_files_[ZONES_STAT], "#   * Number of inhibition genes\n");
    fprintf(stat_files_[ZONES_STAT], "#   * Geometric area of the activation genes\n");
    fprintf(stat_files_[ZONES_STAT], "#   * Geometric area of the inhibition genes\n");
    fprintf(stat_files_[ZONES_STAT], "#   * Geometric area of the resulting phenotype\n");
    fprintf(stat_files_[ZONES_STAT], "#\n");
  }

  if(stat_files_[OPERONS_STAT] != nullptr) {
    fprintf(stat_files_[OPERONS_STAT], "# Each line contains : Generation, and then, for 20 RNA, the number of genes inside the RNA\n");
  }

  if(stat_files_[FIXED_MUTATION_STAT] != nullptr) {
    fprintf(stat_files_[FIXED_MUTATION_STAT], "# #################################################################\n");
    fprintf(stat_files_[FIXED_MUTATION_STAT], "#                      Mutations in the lineage                   #\n");
    fprintf(stat_files_[FIXED_MUTATION_STAT], "# #################################################################\n");
    fprintf(stat_files_[FIXED_MUTATION_STAT], "#  1.  Generation       (mut. occurred when producing the indiv. of this generation)\n");
    fprintf(stat_files_[FIXED_MUTATION_STAT], "#  2.  Genetic unit     (which underwent the mutation, 0 = chromosome) \n");
    fprintf(stat_files_[FIXED_MUTATION_STAT], "#  3.  Mutation type    (0: switch, 1: smallins, 2: smalldel, 3:dupl, 4: del, 5:trans, 6:inv, 7:insert, 8:ins_HT, 9:repl_HT) \n");
    fprintf(stat_files_[FIXED_MUTATION_STAT], "#  4.  pos_0            (position for the small events, begin_segment for the rearrangements, begin_segment of the inserted segment for ins_HT, begin_segment of replaced segment for repl_HT) \n");
    fprintf(stat_files_[FIXED_MUTATION_STAT], "#  5.  pos_1            (-1 for the small events, end_segment for the rearrangements, end_segment of the inserted segment for ins_HT, begin_segment of donor segment for repl_HT) \n");
    fprintf(stat_files_[FIXED_MUTATION_STAT], "#  6.  pos_2            (reinsertion point for duplic., cutting point in segment for transloc., insertion point in the receiver for ins_HT, end_segment of the replaced segment for repl_HT, -1 for other events)\n");
    fprintf(stat_files_[FIXED_MUTATION_STAT], "#  7.  pos_3            (reinsertion point for transloc., breakpoint in the donor for ins_HT, end_segment of the donor segment for repl_HT, -1 for other events)\n");
    fprintf(stat_files_[FIXED_MUTATION_STAT], "#  8.  invert           (transloc, was the segment inverted (0/1)?, sense of insertion for ins_HT (0=DIRECT, 1=INDIRECT), sense of the donor segment for repl_HT (0=DIRECT, 1=INDIRECT),-1 for other events)\n");
    fprintf(stat_files_[FIXED_MUTATION_STAT], "#  9.  align_score      (score that was needed for the rearrangement to occur, score of the first alignment for ins_HT and repl_HT)\n");
    fprintf(stat_files_[FIXED_MUTATION_STAT], "#  10. align_score2     (score for the reinsertion for transloc, score of the second alignment for ins_HT and repl_HT)\n");
    fprintf(stat_files_[FIXED_MUTATION_STAT], "#  11. seg_len          (segment length for rearrangement, donor segment length for ins_HT and repl_HT)\n");
    fprintf(stat_files_[FIXED_MUTATION_STAT], "#  12. repl_seg_len     (replaced segment length for repl_HT, -1 for the others)\n");
    fprintf(stat_files_[FIXED_MUTATION_STAT], "#  13. GU_length        (before the event)\n");
    fprintf(stat_files_[FIXED_MUTATION_STAT], "#  14. Impact of the mutation on the metabolic error (negative value = smaller gap after = beneficial mutation) \n");
    fprintf(stat_files_[FIXED_MUTATION_STAT], "#  15. Impact of the mutation on the fitness \n");
    fprintf(stat_files_[FIXED_MUTATION_STAT], "####################################################################################################################\n");
    fprintf(stat_files_[FIXED_MUTATION_STAT], "#\n");
    fprintf(stat_files_[FIXED_MUTATION_STAT], "# Header for R\n");
    fprintf(stat_files_[FIXED_MUTATION_STAT], "# gener gen_unit mut_type pos_0 pos_1 pos_2 pos_3 invert align_score align_score_2 seg_len repl_seg_len GU_len met_err_impact fit_impact\n");
  }
}

void AncestorStats::Resume_stats(int64_t time) {
  char* cur_file_name;  // Syntaxic sugar for stat_files_names_[]
  FILE* cur_file;
  char* tmp_file_name = new char[100];
  FILE* tmp_file;
  char  line[500];

  for(int8_t s_type = 0 ; s_type < NB_TYPES ; s_type++) {
    cur_file_name = stat_files_names_[s_type];
    if(cur_file_name != nullptr) {
      sprintf(tmp_file_name, "%s.tmp", cur_file_name);

      cur_file = fopen(cur_file_name, "r");
      tmp_file = fopen(tmp_file_name, "w");

      if (fgets(line, 500, cur_file) == nullptr) {
        // TODO check for error
      }
      while (!feof(cur_file) && line[0] == '#') {
        fputs(line, tmp_file);
        if (fgets(line, 500, cur_file) == nullptr) {
          // TODO check for error
        }
      }
      // Copy stats until time (included)
      while ((int64_t)atol(line) <= time && !feof(cur_file)) {
        fputs(line, tmp_file);
        if (fgets(line, 500, cur_file)) {
          // TODO check for error
        }
      }
      fclose(cur_file);
      fclose(tmp_file);

      remove(cur_file_name);
      rename(tmp_file_name, cur_file_name);

      stat_files_[s_type] = fopen(cur_file_name, "a");
    }
  }
  delete[] tmp_file_name;
}

void AncestorStats::write_environment_stats(int64_t t, const std::shared_ptr<PhenotypicTargetHandler> pth,
                             FILE* env_output_file) {
  // Num gener
  fprintf(env_output_file, "%" PRId64, t);

  for (const Gaussian& g: pth->gaussians())
    fprintf(env_output_file,
            "     %.16f %.16f %.16f",
            g.mean(), g.width(), g.height());

  fprintf(env_output_file, "\n");
}

void AncestorStats::write_terminators_stats(int64_t t,  Individual* indiv,
                             FILE* term_output_file) {
  fprintf(term_output_file, "%" PRId64 " %" PRId32 " %" PRId32 "\n",
            t,
            indiv->total_genome_size(),
            indiv->nb_terminators());
}

void AncestorStats::write_zones_stats(int64_t t,
                       Individual* indiv,
                       const std::shared_ptr<PhenotypicTargetHandler> phenotypicTargetHandler,
                       FILE* zones_output_file)
{
  assert(phenotypicTargetHandler->phenotypic_target().nb_segments() > 1);

  int16_t nb_segments = phenotypicTargetHandler->phenotypic_target().nb_segments();
  int16_t num_segment = 0;
  PhenotypicSegment** segments =
      phenotypicTargetHandler->phenotypic_target().segments();

  // Tables : index 0 for the 0 segment
  //                1 for the neutral segment
  int32_t nb_genes_activ[nb_segments];
  int32_t nb_genes_inhib[nb_segments];
  double  geom_area_activ[nb_segments];
  double  geom_area_inhib[nb_segments];
  double  geom_area_phen[nb_segments];

  for (num_segment = 0 ; num_segment < nb_segments ; num_segment++)
  {
    nb_genes_activ[num_segment]   = 0;
    nb_genes_inhib[num_segment]   = 0;
    geom_area_activ[num_segment]  = 0.0;
    geom_area_inhib[num_segment]  = 0.0;
    geom_area_phen[num_segment]   = 0.0;
  }


  AbstractFuzzy* activ = NULL;
  AbstractFuzzy* inhib = NULL;
  Phenotype* phen  = NULL;



  // Compute number of genes in each segment
  for (const auto& prot: indiv->protein_list()) {
    // Go to the corresponding segment
    num_segment = 0;
    while (prot->mean() > segments[num_segment]->stop)
    {
      num_segment++;
    }

    // Add a genes (activ or inhib)
    if (prot->is_functional())
    {
      if (prot->height() > 0)
      {
        nb_genes_activ[num_segment]++;
      }
      else if (prot->height() < 0)
      {
        nb_genes_inhib[num_segment]++;
      }

      // It the gene is exactly at the frontier between 2 zones, mark it in both
      if (prot->mean() == segments[num_segment]->stop && num_segment < nb_segments - 1)
      {
        if (prot->height() > 0)
        {
          nb_genes_activ[num_segment+1]++;
        }
        else if (prot->height() < 0)
        {
          nb_genes_inhib[num_segment+1]++;
        }
      }
    }
  }

  // Compute the geometric areas
  activ = indiv->phenotype_activ();
  inhib = indiv->phenotype_inhib();
  phen  = indiv->phenotype();

  for (num_segment = 0 ; num_segment < nb_segments ; num_segment++)
  {
    geom_area_activ[num_segment]  = activ->get_geometric_area(segments[num_segment]->start, segments[num_segment]->stop);
    geom_area_inhib[num_segment]  = inhib->get_geometric_area(segments[num_segment]->start, segments[num_segment]->stop);
    geom_area_phen[num_segment]   = phen->get_geometric_area(segments[num_segment]->start, segments[num_segment]->stop);
  }


  // Print stats to file
  fprintf(zones_output_file, "%" PRId64, t);

  for (num_segment = 0 ; num_segment < nb_segments ; num_segment++)
  {
    fprintf(zones_output_file, "     %" PRId32 " %" PRId32 " %lf %lf %lf",
              nb_genes_activ[num_segment],
              nb_genes_inhib[num_segment],
              geom_area_activ[num_segment],
              geom_area_inhib[num_segment],
              geom_area_phen[num_segment]);
  }

  fprintf(zones_output_file, "\n");
}

void AncestorStats::write_operons_stats(int64_t t, Individual* indiv, FILE*  operons_output_file)
{
  int32_t nb_genes_per_rna[20];
  for (int i = 0 ; i < 20 ; i++)
  {
    nb_genes_per_rna[i] = 0;
  }

  for (const auto& rna: indiv->rna_list()) {
    if (rna->transcribed_proteins().size() >= 20)
    {
      printf("Found operon with 20 genes or more : %zu\n", rna->transcribed_proteins().size());
    }

    nb_genes_per_rna[rna->transcribed_proteins().size()]++;
  }

  fprintf(operons_output_file, "%" PRId64 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 " %" PRId32 "\n",
            t,
            nb_genes_per_rna[0],
            nb_genes_per_rna[1],
            nb_genes_per_rna[2],
            nb_genes_per_rna[3],
            nb_genes_per_rna[4],
            nb_genes_per_rna[5],
            nb_genes_per_rna[6],
            nb_genes_per_rna[7],
            nb_genes_per_rna[8],
            nb_genes_per_rna[9],
            nb_genes_per_rna[10],
            nb_genes_per_rna[11],
            nb_genes_per_rna[12],
            nb_genes_per_rna[13],
            nb_genes_per_rna[14],
            nb_genes_per_rna[15],
            nb_genes_per_rna[16],
            nb_genes_per_rna[17],
            nb_genes_per_rna[18],
            nb_genes_per_rna[19]);
}

void AncestorStats::save_grid_cell(std::string grid_file_name) const {

  gzFile grid_file = gzopen(grid_file_name.c_str(), "w");
  indiv_->grid_cell()->save(grid_file);

  gzclose(grid_file);
}

void AncestorStats::load_grid_cell(ExpManager* exp_m) {
  char grid_file_name[255];
  sprintf(grid_file_name, STATS_DIR"/ancestor_stats/saved_gridcell.ae");

  gzFile grid_file = gzopen(grid_file_name, "r");

  grid_cell_ = new GridCell(grid_file, exp_m, nullptr);
  indiv_ = grid_cell_->individual();

  gzclose(grid_file);
}

} // namespace aevol
