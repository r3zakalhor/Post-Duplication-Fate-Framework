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
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <err.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <list>




// =================================================================
//                            Project Files
// =================================================================
#include "aevol.h"


using namespace aevol;


enum check_type
{
  FULL_CHECK  = 0,
  LIGHT_CHECK = 1,
  ENV_CHECK   = 2,
  NO_CHECK    = 3
};



// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);

FILE* open_environment_stat_file(const char * prefix);
void write_environment_stats(int64_t t,
                             const PhenotypicTargetHandler* pth,
                             FILE* env_file);

FILE* open_terminators_stat_file(const char * prefix);
void write_terminators_stats(int64_t t,  Individual* indiv, FILE* terminator_file);

FILE* open_zones_stat_file(const char * prefix);
void write_zones_stats(int64_t t,
                       Individual* indiv,
                       const PhenotypicTargetHandler* phenotypicTargetHandler,
                       FILE* zone_file);

FILE* open_operons_stat_file(const char * prefix);
void write_operons_stats(int64_t t, Individual* indiv, FILE* operon_file);


double* dist_to_target_segment;


int main(int argc, char** argv)
{
  // The input file (lineage.ae or lineage.rae) must contain the following information:
  //
  // - common data                                                (ae_common::write_to_backup)
  // - begin gener                                                (int64_t)
  // - end gener                                                  (int64_t)
  // - final individual index                                     (int32_t)
  // - initial genome size                                        (int32_t)
  // - initial ancestor (nb genetic units + sequences)            (Individual::write_to_backup)
  // - replication report of ancestor at time t0+1  (ae_replic_report::write_to_backup)
  // - replication report of ancestor at time t0+2  (ae_replic_report::write_to_backup)
  // - replication report of ancestor at time t0+3  (ae_replic_report::write_to_backup)
  // - ...
  // - replication report of ancestor at time t_end_ (ae_replic_report::write_to_backup)




  // =====================
  //  Parse command line
  // =====================

  // Default values
  char*       lineage_file_name   = NULL;
  bool        verbose             = false;
  check_type  check               = LIGHT_CHECK;
  double      tolerance           = 0;

  const char * short_options = "hVvncf:lt:";
  static struct option long_options[] =
  {
    {"help",        no_argument,       NULL, 'h'},
    {"version",     no_argument,       NULL, 'V' },
    {"verbose",     no_argument,       NULL, 'v'},
    {"nocheck",     no_argument,       NULL, 'n'},
    {"fullcheck",   no_argument,       NULL, 'c'},
    {"file",        required_argument, NULL, 'f'},
    {"tolerance",   required_argument, NULL, 't'},
    {0, 0, 0, 0}
  };

  int option;
  while ((option = getopt_long(argc, argv, short_options, long_options, NULL)) != -1)
  {
    switch(option)
    {
      case 'h' :
      {
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      }
      case 'V' :
      {
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      }
      case 'v' : verbose = true;                    break;
      case 'n' : check = NO_CHECK;                  break;
      case 'c' : check = FULL_CHECK;                break;
      case 'f' :
      {
        if (strcmp(optarg, "") == 0)
        {
          fprintf(stderr, "ERROR : Option -f or --file : missing argument.\n");
          exit(EXIT_FAILURE);
        }
        lineage_file_name = new char[strlen(optarg) + 1];
        sprintf(lineage_file_name, "%s", optarg);
        break;
      }
      case 't' :
      {
        if (strcmp(optarg, "") == 0)
        {
          fprintf(stderr, "ERROR : Option -t or --tolerance : missing argument.\n");
          exit(EXIT_FAILURE);
        }
        check = ENV_CHECK;
        tolerance = atof(optarg);
        break;
      }
      default :
      {
        fprintf(stderr, "ERROR : Unknown option, check your syntax.\n");
        print_help(argv[0]);
        exit(EXIT_FAILURE);
      }
    }
  }



  if (lineage_file_name == NULL)
  {
    fprintf(stderr, "ERROR : Option -f or --file missing. \n");
    exit(EXIT_FAILURE);
  }


  printf("\n");
  printf("WARNING : Parameter change during simulation is not managed in general.\n");
  printf("          Only changes in environmental target done with aevol_modify are handled.\n");
  printf("\n");

  // =======================
  //  Open the lineage file
  // =======================
  gzFile lineage_file = gzopen(lineage_file_name, "r");
  if (lineage_file == Z_NULL)
  {
    fprintf(stderr, "ERROR : Could not read the lineage file %s\n", lineage_file_name);
    exit(EXIT_FAILURE);
  }

  int64_t t0 = 0;
  int64_t t_end = 0;
  int32_t final_indiv_index = 0;
  int32_t final_indiv_rank  = 0;


  gzread(lineage_file, &t0, sizeof(t0));
  gzread(lineage_file, &t_end, sizeof(t_end));
  gzread(lineage_file, &final_indiv_index, sizeof(final_indiv_index));
  gzread(lineage_file, &final_indiv_rank,  sizeof(final_indiv_rank));

  if (verbose)
  {
    printf("\n\n");
    printf("===============================================================================\n");
    printf(" Statistics of the ancestors of indiv. %" PRId32
           " (rank %" PRId32 ") from time %" PRId64 " to %" PRId64 "\n",
           final_indiv_index, final_indiv_rank, t0, t_end);
    printf("================================================================================\n");
  }



  // =============================
  //  Open the experience manager
  // =============================
  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(t0, true, false);

  // The current version doesn't allow for phenotypic variation nor for
  // different phenotypic targets among the grid
  if (not exp_manager->world()->phenotypic_target_shared())
    Utils::ExitWithUsrMsg("sorry, ancestor stats has not yet been implemented "
                              "for per grid-cell phenotypic target");
  auto phenotypicTargetHandler =
      exp_manager->world()->phenotypic_target_handler();
  if (not (phenotypicTargetHandler->var_method() == NO_VAR))
    Utils::ExitWithUsrMsg("sorry, ancestor stats has not yet been implemented "
                              "for variable phenotypic targets");

  int64_t backup_step = exp_manager->backup_step();


  // =========================
  //  Open the output file(s)
  // =========================
  // Create missing directories
  int status;
  status = mkdir("stats/ancstats/", 0755);
  if ((status == -1) && (errno != EEXIST))
  {
    err(EXIT_FAILURE, "stats/ancstats/");
  }

  auto prefix = "ancestor_stats/ancestor_stats";
    char postfix[255];
    snprintf(postfix, 255,
             "-b" TIMESTEP_FORMAT "-e" TIMESTEP_FORMAT "-i%" PRId32 "-r%" PRId32,
            t0, t_end, final_indiv_index, final_indiv_rank);
  bool best_indiv_only = true;
  bool addition_old_stats = false;
  bool delete_old_stats = true;
    Stats* mystats = new Stats(exp_manager, t0, best_indiv_only, prefix, postfix,
                               addition_old_stats, delete_old_stats);
  //mystats->write_headers();

  // Optional outputs
  FILE* env_output_file = open_environment_stat_file(prefix);
  FILE* term_output_file = open_terminators_stat_file(prefix);
  FILE* zones_output_file = NULL;

  // Next line patchy (specific for the constraints mentioned earlier, i.e.
  // works only for shared and unvarying phenotypic target)
  if (phenotypicTargetHandler->phenotypic_target().nb_segments() > 1)
  {
    zones_output_file = open_zones_stat_file(prefix);
  }
  FILE* operons_output_file = open_operons_stat_file(prefix);


  // ==================================================
  //  Prepare the initial ancestor and write its stats
  // ==================================================
  GridCell* grid_cell = new GridCell(lineage_file, exp_manager, nullptr);

  // Individual*indiv = Individual::CreateIndividual(exp_manager, lineage_file);
  auto* indiv = grid_cell->individual();
  indiv->Evaluate();
  indiv->compute_statistical_data();
  indiv->compute_non_coding();

  mystats->write_statistics_of_this_indiv(t0,indiv,nullptr);

  // Optional outputs
  write_environment_stats(t0, phenotypicTargetHandler, env_output_file);
  write_terminators_stats(t0, indiv, term_output_file);
  if(phenotypicTargetHandler->phenotypic_target().nb_segments() > 1)
  {
    write_zones_stats(t0, indiv, phenotypicTargetHandler, zones_output_file);
  }
  write_operons_stats(t0, indiv, operons_output_file);


  if (verbose)
  {
    printf("Initial fitness     = %f\n", indiv->fitness());
    printf("Initial genome size = %" PRId32 "\n", indiv->total_genome_size());
  }

  //delete exp_manager;

  // ==========================================================================
  //  Replay the mutations to get the successive ancestors and analyze them
  // ==========================================================================
  ReplicationReport* rep = nullptr;

  int32_t index;

  ExpManager* exp_manager_backup = nullptr;
  Habitat *backup_habitat = nullptr;

  bool check_now = false;

  aevol::AeTime::plusplus();
  while (time() <= t_end)
  {
    rep = new ReplicationReport(lineage_file, indiv);
    index = rep->id(); // who we are building...

    // Check now?
    check_now = ((check == FULL_CHECK && Utils::mod(time(), backup_step) == 0) ||
                 (check == ENV_CHECK && Utils::mod(time(), backup_step) == 0) ||
                 (check == LIGHT_CHECK && time() == t_end));

    if (verbose)
        printf("Rebuilding ancestor at generation %" PRId64
            " (index %" PRId32 ")...", time(), index);

    indiv->Reevaluate();

    // TODO <david.parsons@inria.fr> Check for phenotypic variation has to be
    // done for all the grid cells, disable checking until coded

//    // Check, and possibly update, the environment according to the backup files
//    // (update necessary if the env. was modified by aevol_modify at some point)
//    if (Utils::mod(time(), backup_step) == 0)
//    {
//      char world_file_name[255];
//      sprintf(world_file_name, "./" WORLD_FNAME_FORMAT, time());
//      gzFile world_file = gzopen(world_file_name, "r");
//      backup_habitat = new Habitat(world_file, pth); // TODO vld: fix pth
//
//      if (! env->is_identical_to(*backup_env, tolerance))
//      {
//        printf("Warning: At time()=%" PRId64 ", the replayed environment is not the same\n", time());
//        printf("         as the one saved at time()=%" PRId64 "... \n", time());
//        printf("         with tolerance of %lg\n", tolerance);
//        printf("Replacing the replayed environment by the one stored in the backup.\n");
//        delete env;
//        h = new Habitat(*backup_habitat);
//      }
//      delete backup_habitat;
//    }


    // Warning: this portion of code won'time() work if the number of units changes
    // during the evolution

    // 2) Replay replication (create current individual's child)
    GeneticUnit& gen_unit = indiv->genetic_unit_nonconst(0);
    GeneticUnit* stored_gen_unit = nullptr;
    Individual* stored_indiv = nullptr;

    if (check_now)
    {
      exp_manager_backup = new ExpManager();
      exp_manager_backup->load(time(), true, false);
      stored_indiv = new Individual(
          *(Individual*) exp_manager_backup->indiv_by_id(index));
      stored_gen_unit = &(stored_indiv->genetic_unit_nonconst(0));
    }

    // For each genetic unit, replay the replication (undergo all mutations)
    // TODO <david.parsons@inria.fr> disabled for multiple GUs
    const auto& dnarep = rep->dna_replic_report();

    dnarep.iter_muts([&](const auto& mut) {
      gen_unit.dna()->undergo_this_mutation(*mut);
    });

    if (check_now)
    {
      if (verbose)
      {
        printf("Checking the sequence of the unit...");
        fflush(NULL);
      }

      char * str1 = new char[gen_unit.dna()->length() + 1];
      memcpy(str1, gen_unit.dna()->data(), \
             gen_unit.dna()->length()*sizeof(char));
      str1[gen_unit.dna()->length()] = '\0';

      char * str2 = new char[(stored_gen_unit->dna())->length() + 1];
      memcpy(str2, (stored_gen_unit->dna())->data(),
             (stored_gen_unit->dna())->length()*sizeof(char));
      str2[(stored_gen_unit->dna())->length()] = '\0';

      if (strncmp(str1, str2, stored_gen_unit->dna()->length()) == 0) {
        if (verbose)
          printf(" OK\n");
      }
      else {
        if (verbose) printf(" ERROR !\n");
        fprintf(stderr, "Error: the rebuilt genetic unit is not the same as \n");
        fprintf(stderr, "the one saved at generation %" PRId64 "... ", time());
        fprintf(stderr, "Rebuilt unit : %" PRId32 " bp\n %s\n", (int32_t)strlen(str1), str1);
        fprintf(stderr, "Stored unit  : %" PRId32 " bp\n %s\n", (int32_t)strlen(str2), str2);

        delete [] str1;
        delete [] str2;
        gzclose(lineage_file);
        delete indiv;
        delete stored_indiv;
        delete exp_manager_backup;
        delete exp_manager;
        exit(EXIT_FAILURE);
      }

      delete [] str1;
      delete [] str2;
    }

    // 3) All the mutations have been replayed, we can now evaluate the new individual
    indiv->Reevaluate();
    indiv->compute_statistical_data();
    indiv->compute_non_coding();

    mystats->write_statistics_of_this_indiv(time(),indiv,rep);

    // Optional outputs
    write_environment_stats(time(), phenotypicTargetHandler, env_output_file);
    write_terminators_stats(time(), indiv, term_output_file);
    if(phenotypicTargetHandler->phenotypic_target().nb_segments() > 1)
    {
      write_zones_stats(time(), indiv, phenotypicTargetHandler, zones_output_file);
    }
    write_operons_stats(time(), indiv, operons_output_file);

    if (verbose) printf(" OK\n");

    delete rep;

    if (check_now)
    {
      delete stored_indiv;
      delete exp_manager_backup;
    }

    aevol::AeTime::plusplus();
  }

  gzclose(lineage_file);

  // Optional outputs
  fclose(env_output_file);
  fclose(term_output_file);
  if(phenotypicTargetHandler->phenotypic_target().nb_segments() > 1)
  {
    fclose(zones_output_file);
  }
  fclose(operons_output_file);

  delete exp_manager;
  delete mystats;
  delete indiv;

  exit(EXIT_SUCCESS);
}




FILE* open_environment_stat_file(const char * prefix)
{
  // Open file
  char* env_output_file_name = new char[80];
  sprintf(env_output_file_name, "stats/%s_envir.out",prefix);
  FILE* env_output_file = fopen(env_output_file_name, "w");
  delete env_output_file_name;

  // Write headers
  // TODO vld: was limited to "if environment->gaussians_provided"
  // are gaussians always available now?
  fprintf(env_output_file, "# Each line contains : Generation, and then, for each gaussian: M W H.\n");
  fprintf(env_output_file, "#\n");

  return env_output_file;
}


void write_environment_stats(int64_t t,
                             const PhenotypicTargetHandler* pth,
                             FILE* env_output_file)
{
  // Num gener
  fprintf(env_output_file, "%" PRId64, t);

  for (const Gaussian& g: pth->gaussians())
    fprintf(env_output_file,
            "     %.16f %.16f %.16f",
            g.mean(), g.width(), g.height());

  fprintf(env_output_file, "\n");
}



FILE* open_terminators_stat_file(const char * prefix)
{
  char* term_output_file_name = new char[80];
  sprintf(term_output_file_name, "stats/%s_nb_term.out",prefix);
  FILE* term_output_file = fopen(term_output_file_name, "w");
  delete [] term_output_file_name;

  // Write headers
  fprintf(term_output_file, "# Each line contains : \n");
  fprintf(term_output_file, "#   * Generation\n");
  fprintf(term_output_file, "#   * Genome size\n");
  fprintf(term_output_file, "#   * Terminator number\n");
  fprintf(term_output_file, "#\n");

  return term_output_file;
}

void write_terminators_stats(int64_t t,  Individual* indiv, FILE* term_output_file)
{
  fprintf(term_output_file, "%" PRId64 " %" PRId32 " %" PRId32 "\n",
            t,
            indiv->total_genome_size(),
            indiv->nb_terminators());
}



FILE* open_zones_stat_file(const char * prefix)
{
  // Open file
  char* zones_output_file_name = new char[80];
  sprintf(zones_output_file_name, "stats/%s_zones.out",prefix);
  FILE* zones_output_file = fopen(zones_output_file_name, "w");
  delete [] zones_output_file_name;

  // Write headers
  fprintf(zones_output_file, "# Each line contains : Generation, and then, for each zone:\n");
  fprintf(zones_output_file, "#   * Number of activation genes\n");
  fprintf(zones_output_file, "#   * Number of inhibition genes\n");
  fprintf(zones_output_file, "#   * Geometric area of the activation genes\n");
  fprintf(zones_output_file, "#   * Geometric area of the inhibition genes\n");
  fprintf(zones_output_file, "#   * Geometric area of the resulting phenotype\n");
  fprintf(zones_output_file, "#\n");

  return zones_output_file;
}

void write_zones_stats(int64_t t,
                       Individual* indiv,
                       const PhenotypicTargetHandler* phenotypicTargetHandler,
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



FILE* open_operons_stat_file(const char * prefix)
{
  char* operons_output_file_name = new char[80];
  sprintf(operons_output_file_name, "stats/%s_operons.out",prefix);
  FILE* operons_output_file = fopen(operons_output_file_name, "w");
  delete [] operons_output_file_name,

  // Write headers
  fprintf(operons_output_file, "# Each line contains : Generation, and then, for 20 RNA, the number of genes inside the RNA\n");
  return operons_output_file;
}

void write_operons_stats(int64_t t, Individual* indiv, FILE*  operons_output_file)
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




/*!
  \brief

*/
void print_help(char* prog_path)
{
  printf("\n");
  printf("*********************** aevol - Artificial Evolution ******************* \n");
  printf("*                                                                      * \n");
  printf("*                      Ancstats post-treatment program                 * \n");
  printf("*                                                                      * \n");
  printf("************************************************************************ \n");
  printf("\n\n");
  printf("This program is Free Software. No Warranty.\n");
  printf("Copyright (C) 2009  LIRIS.\n");
  printf("\n");
#ifdef __REGUL
  printf("Usage : rancstats -h\n");
  printf("or :    rancstats [-vn] -f lineage_file \n");
#else
  printf("Usage : ancstats -h\n");
  printf("or :    ancstats [-vn] -f lineage_file \n");
#endif
  printf("\n");
  printf("This program compute some statistics for the individuals within lineage_file.\n");
  printf("\n");
  printf("\n");
  printf("\t-h or --help       : Display this help.\n");
  printf("\n");
  printf("\t-v or --verbose    : Be verbose, listing generations as they are \n");
  printf("\t                       treated.\n");
  printf("\n");
  printf("\t-n or --nocheck    : Disable genome sequence checking. Makes the \n");
  printf("\t                       program faster, but it is not recommended. \n");
  printf("\t                       It is better to let the program check that \n");
  printf("\t                       when we rebuild the genomes of the ancestors\n");
  printf("\t                       from the lineage file, we get the same sequences\n");
  printf("\t                       as those stored in the backup files.\n");
  printf("\n");
  printf("\t-c or --fullcheck  : Will perform the genome and environment checks every\n");
  printf("\t                       <BACKUP_STEP> generations. Default behaviour is\n");
  printf("\t                       lighter as it only perform sthese checks at the\n");
  printf("\t                       ending generation.\n");
  printf("\n");
  printf("\t-f lineage_file or --file lineage_file : \n");
  printf("\t                       Compute the statistics for the individuals within lineage_file.\n");
  printf("\t-t tolerance or --tolerance tolerance : \n");
  printf("\t                       Tolerance used to compare the replayed environment to environment in backup\n");
  printf("\n");
}
