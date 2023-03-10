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



//
// This program extracts data about the individuals and write
// them into text files easy to parse with matlab.
//
// Options are :
//     * -f / --file   backup_file
//            backup file where individuals should be extracted
//     * -o / --output   main_output_file
//            output directory where statistics should be printed
//     * -n / --neutral  output_file
//            print position of neutral regions in 'ouput_file'
//     * -b / --best
//            process only the best individual
//
// Data concerning genomes is printed in 'main_output_file'. A space delimits two pieces
// of information, a new line two individuals. Format is
//       "nc1 nc2 nc3 nc4 nc5 nc6 nc7 total\n"
// where:
//     * nc1: number of bases in neutral regions
//     * nc2: number of bases outside CDS
//     * nc3: number of bases outside functional CDS
//     * nc4: number of bases outside non functional CDS
//     * nc5: number of bases outside RNAs
//     * nc6: number of bases outside functional RNAs
//     * nc7: number of bases outside non functional RNAs
//     * total: total size of genome
//
// It is also possible to print neutral regions with the '-n output_file' option.
// In this case, "output_file" contains:
//       "# chromosome length: cl, nb neutral bases, nr neutral regions\n
//       "bnr_1 bnr_2 ... bnr_n\n"
//       "enr_1 enr_2 ... enr_n\n"
// where:
//     * cl:   number of bases in chromosome
//     * nb:   number of bases in neutral regions
//     * nr:   number of neutral regions
//     * bnr_i: beginning of ith neutral region
//     * enr_i: end of ith neutral region
//
// Examples :
//
// For generation 20000, write info about the genomes of all the
// individuals in out_020000 :
//
//    genome_stats -f backup/gen_020000.ae -o out_020000
//
// For generation 20000, write the best individual's info in
// out_020000_best :
//
//    genome_stats -b -f backup/gen_020000.ae -o out_020000_best
//
// For generation 20000, write beginning and ends of neutral regions
// in neutral.out for all the individuals
//
//    genome_stats -f backup/gen_020000.ae -n neutral.out
//




// =================================================================
//                              Libraries
// =================================================================
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>

#include <list>

// =================================================================
//                            Project Files
// =================================================================
#include "aevol.h"

using namespace aevol;

// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);

// print information about the indivdual's genome to file
void print_genome_info(ae_individual* indiv, FILE* output_file);

// print information about the indivdual's neutral regions to file
void print_neutral_regions(ae_individual* indiv, FILE* output_file);


// =================================================================
//                         Main Function
// =================================================================

int main(int argc, char* argv[])
{
  // ---------------------------------------
  //      command-line option parsing
  // ---------------------------------------
  // Initialize command-line option variables with default values
  char* backup_file_name           = NULL;
  char* main_output_name           = NULL;
  char* neutral_region_output_name = NULL;
  bool best_only = false;

  // Define allowed options
  const char * options_list = "hVf:o:bn:";
  static struct option long_options_list[] =
  {
    {"help",    no_argument,        NULL, 'h'},
    {"version", no_argument,        NULL, 'V'},
    {"file",    required_argument,  NULL, 'f' },
    {"output",  required_argument,  NULL, 'o' },
    {"best",    no_argument,        NULL, 'b' },
    {"neutral", no_argument,        NULL, 'n' },
    { 0, 0, 0, 0 }
  };

  // Get actual values of the command-line options
  int option;
  while ((option = getopt_long(argc, argv, options_list, long_options_list, NULL)) != -1)
  {
    switch (option)
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
      case 'f' :
        backup_file_name = new char[strlen(optarg) + 1];
        sprintf(backup_file_name, "%s", optarg);
        break;
      case 'o' :
        main_output_name = new char[strlen(optarg) + 1];
        sprintf(main_output_name, "%s", optarg);
        break;
      case 'n' :
        neutral_region_output_name = new char[strlen(optarg) + 1];
        sprintf(neutral_region_output_name, "%s", optarg);
        break;
      case 'b' :
        best_only = true;
        break;
    }
  }

  // -------------------------------
  //          Initialize
  // -------------------------------
  FILE* main_output           = NULL;
  FILE* neutral_region_output = NULL;

  if (backup_file_name == NULL)
  {
    printf("You must specify a backup file. Please use the option -f or --file.\n");
    exit(EXIT_FAILURE);
  }
  if (main_output_name != NULL)
  {
    main_output = fopen(main_output_name,"w");
    if (main_output == NULL)
    {
      fprintf(stderr, "Warning: Could not open file %s.\n", main_output_name);
    }
  }
  if (neutral_region_output_name != NULL)
  {
    neutral_region_output = fopen(neutral_region_output_name,"w");
    if (neutral_region_output == NULL)
    {
      fprintf(stderr, "Warning: Could not open file %s.\n", neutral_region_output_name);
    }
  }
  fflush(stderr);

  printf("Reading backup file <%s>... \n", backup_file_name);
  fflush(stdout);

  // Load the simulation from backup
  ae_common::sim = new ae_experiment();
  ae_common::sim->load_backup(backup_file_name, false, NULL);
  printf("done\n");
  delete [] backup_file_name;

  printf("Computing phenotypes... \n");
  fflush(stdout);

  // Evaluate the individuals
  ae_common::pop->evaluate_individuals(ae_common::sim->env());

  int i = 0;
  int nb_indiv = ae_common::pop->nb_indivs();

  // --------------------------------
  //         Parse individuals
  // --------------------------------
  if (best_only)
  {
    ae_individual* best = ae_common::pop->best();
    if (main_output != NULL)           { print_genome_info(best, main_output); }
    if (neutral_region_output != NULL) { print_neutral_regions(best, neutral_region_output); }
  }
  else
  {
    if (ae_common::pop_structure)
    {
      ae_grid_cell*** pop_grid_ = ae_common::pop->pop_grid();
      for (int16_t x = 0 ; x < ae_common::grid_x ; x++)
      {
        for (int16_t y = 0 ; y < ae_common::grid_y ; y++)
        {
          ae_individual* indiv = (pop_grid_[x][y]->individual());
	  if (main_output != NULL)           { print_genome_info(indiv, main_output); }
	  if (neutral_region_output != NULL) { print_neutral_regions(indiv, neutral_region_output); }
          i++;
        }
      }
    }
    else
    {
      for (const auto& indiv: pop->indivs()) {
        if (main_output != NULL)           { print_genome_info(indiv, main_output); }
	if (neutral_region_output != NULL) { print_neutral_regions(indiv, neutral_region_output); }
      }
    }
  }

  if (main_output != NULL) { fclose(main_output); }
  if (neutral_region_output != NULL) { fclose(neutral_region_output); }

  if (main_output_name != NULL)           { delete [] main_output_name; }
  if (neutral_region_output_name != NULL) { delete [] neutral_region_output_name; }

  ae_common::clean();

  return EXIT_SUCCESS;
}


// =================================================================
//              Implementation of Secondary Functions
// =================================================================


// The export fonction
inline void print_genome_info(ae_individual* indiv, FILE* output_file)
{
  int32_t nb_bases = indiv->total_genome_size();
  int32_t nb_bases_in_neutral_regions = indiv->nb_bases_in_neutral_regions();
  int32_t nb_bases_in_0_CDS = indiv->nb_bases_in_0_CDS();
  int32_t nb_bases_in_0_functional_CDS = indiv->nb_bases_in_0_functional_CDS();
  int32_t nb_bases_in_0_non_functional_CDS = indiv->nb_bases_in_0_non_functional_CDS();
  int32_t nb_bases_in_0_RNA = indiv->nb_bases_in_0_RNA();
  int32_t nb_bases_in_0_coding_RNA = indiv->nb_bases_in_0_coding_RNA();
  int32_t nb_bases_in_0_non_coding_RNA = indiv->nb_bases_in_0_non_coding_RNA();
  if (output_file != NULL)
  {
    fprintf(output_file, "%"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32" %"PRId32"\n",
	     nb_bases_in_neutral_regions, nb_bases_in_0_CDS, nb_bases_in_0_functional_CDS,
	     nb_bases_in_0_non_functional_CDS, nb_bases_in_0_RNA, nb_bases_in_0_coding_RNA,
	     nb_bases_in_0_non_coding_RNA, nb_bases);
  }
}

inline void print_neutral_regions(ae_individual* indiv, FILE* output_file)
{
  if (output_file == NULL) return;

  //header
  GeneticUnit* chromosome = *indiv->genetic_unit_list().begin();
  int32_t nb_neutral_regions  = chromosome->nb_neutral_regions();
  fprintf(output_file, "# chromosome length: %"PRId32", %"PRId32" neutral bases, %"PRId32" neutral regions\n",
	   chromosome->dna()->length(), chromosome->nb_bases_in_neutral_regions(), nb_neutral_regions);

  //neutral regions
  if (nb_neutral_regions > 0)
  {
    int32_t* beginning_nr = chromosome->beginning_neutral_regions();
    int32_t* end_nr = chromosome->end_neutral_regions();

    for (int32_t i=0; i<nb_neutral_regions; i++) fprintf(output_file, "%"PRId32"\t", beginning_nr[i]);
    fprintf(output_file, "\n");
    for (int32_t i=0; i<nb_neutral_regions; i++) fprintf(output_file, "%"PRId32"\t", end_nr[i]);
    fprintf(output_file, "\n");
  }
}

// TODO: update
void print_help(char* prog_name)
{
  printf("\n\
Usage : genome_stats -h\n\
   or : genome_stats -f source [-o output_file] [-b] \n\n\
\t-h : display this screen\n\
\t--file source : read from the backup file source\n\
\t--output of : extract and save some infos about the genomes of the individuals to file of\
\t--best : only treat the best individual\n");
}
