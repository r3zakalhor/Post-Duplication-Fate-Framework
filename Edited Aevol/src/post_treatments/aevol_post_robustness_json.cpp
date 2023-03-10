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

// ============================================================================
//                                   Includes
// ============================================================================
#include <list>
#include <getopt.h>

#include "aevol.h"
#include "IndivAnalysis.h"

using namespace aevol;

// Helper functions
void interpret_cmd_line_options(int argc, char* argv[]);
void print_help(char* prog_path);
int write_headers(FILE* output_file,bool full_output);

// Command-line option variables
static char* json_file_name = nullptr;
static int32_t nb_mutants = 1000; //< Number of mutants per individual
static int32_t begin = 0; //< First generation to analyse
static int32_t end = -1; //< Last generation to analyse
static int32_t period = 1; //< Period of analysis
static char* output_file_name = nullptr;
static bool verbose = false;
static bool full_output = false;

int main(int argc, char* argv[]) {
  output_file_name = new char[256];
  strcpy(output_file_name, "robustness_summary.txt");
  interpret_cmd_line_options(argc, argv);

  /* Open output file */
  FILE* output_summary = fopen(output_file_name, "w");
  if (output_summary == nullptr) {
    Utils::ExitWithUsrMsg(std::string("Could not create ") + output_file_name);
  }
  write_headers(output_summary,full_output);


  IOJson* iojson = new IOJson(json_file_name);
  iojson->setRecordTree(false);

  /* Analyse each individuals of the json file */
   for(Individual* indiv: iojson->getIndividuals()) {
    IndivAnalysis wanted_indiv(*indiv);
    wanted_indiv.set_grid_cell(indiv->grid_cell());
    wanted_indiv.grid_cell()->set_individual((Individual*)&wanted_indiv);
    wanted_indiv.Evaluate();
    wanted_indiv.compute_statistical_data();
    wanted_indiv.compute_non_coding();

     if ((time() >= begin) && ((time() < end) || (end == -1))) {

      wanted_indiv.compute_experimental_f_nu(
          nb_mutants, std::make_shared<JumpingMT>(time(nullptr)),
          output_summary, verbose, full_output);
    }
  }

  fclose(output_summary);

  return EXIT_SUCCESS;
}


void interpret_cmd_line_options(int argc, char* argv[]) {
  const char* short_options = "hVvfn:b:e:o:";
  static struct option long_options[] = {
      {"help",        no_argument,       nullptr, 'h'},
      {"version",     no_argument,       nullptr, 'V'},
      {"verbose",     no_argument,       nullptr, 'v'},
      {"full",        no_argument,       nullptr, 'f'},
      {"nb-mutants",  required_argument, nullptr, 'n'},
      {"begin",       required_argument, nullptr, 'b'},
      {"end",         required_argument, nullptr, 'e'},
      {"output",      required_argument, nullptr, 'o'},
      {0, 0, 0, 0}
  };

  int option;
  while ((option = getopt_long(argc, argv, short_options, long_options,
                               nullptr)) != -1) {
    switch (option) {
    case 'h' :
      print_help(argv[0]);
      exit(EXIT_SUCCESS);
    case 'V' :
      Utils::PrintAevolVersion();
      exit(EXIT_SUCCESS);
    case 'v' :
      verbose = true;
      break;
    case 'f' :
      full_output = true;
      break;
    case 'b' :
      begin = atol(optarg);
      break;
    case 'e' :
      end = atol(optarg);
      break;
    case 'n' :
      nb_mutants = atol(optarg);
      break;
    case 'o' :
      delete [] output_file_name;
      output_file_name = new char[strlen(optarg) + 1];
      sprintf(output_file_name, "%s", optarg);
      break;
    default:
      // An error message is printed in getopt_long, we just need to exit
      exit(EXIT_FAILURE);
    }
  }

  // There should be only one remaining arg: the lineage file
  if (optind != argc - 1) {
    Utils::ExitWithUsrMsg("please specify a json file");
  }

  json_file_name = new char[strlen(argv[optind]) + 1];
  sprintf(json_file_name, "%s", argv[optind]);
}

void print_help(char* prog_path) {
  // Get the program file-name in prog_name (strip prog_path of the path)
  char* prog_name; // No new, it will point to somewhere inside prog_path
  if ((prog_name = strrchr(prog_path, '/'))) {
    prog_name++;
  }
  else {
    prog_name = prog_path;
  }

  printf("******************************************************************************\n");
  printf("*                                                                            *\n");
  printf("*                        aevol - Artificial Evolution                        *\n");
  printf("*                                                                            *\n");
  printf("* Aevol is a simulation platform that allows one to let populations of       *\n");
  printf("* digital organisms evolve in different conditions and study experimentally  *\n");
  printf("* the mechanisms responsible for the structuration of the genome and the     *\n");
  printf("* transcriptome.                                                             *\n");
  printf("*                                                                            *\n");
  printf("******************************************************************************\n");
  printf("\n");
  printf("%s: generate and analyse mutants for the provided json file.\n",
         prog_name);
  printf("\n");
  printf("Usage : %s -h or --help\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("   or : %s JSON_FILE [-b TIMESTEP] [-e TIMESTEP] [-n NB_MUTANTS] [-o output] [-v]\n",
         prog_name);
  printf("\nOptions\n");
  printf("  -h, --help\n\tprint this help, then exit\n");
  printf("  -V, --version\n\tprint version number, then exit\n");
  printf("  -b, --begin TIMESTEP\n");
  printf("\ttimestep at which to start the analysis\n");
  printf("  -e, --end TIMESTEP\n");
  printf("\ttimestep at which to stop the analysis\n");
  printf("  -n, --nb-mutants NB_MUTANTS\n");
  printf("\tnumber of mutants to be generated\n");
  printf("  -f, --full\n");
  printf("\tfull output (otherwize synthetic output will be produced). -f option must be used with care as it may produce very large output files\n");
  printf("  -o, --output\n");
  printf("\toutput file name\n");
  printf("  -v, --verbose\n\tbe verbose\n");
}


int write_headers(FILE* output_file,bool full_output) {
  // --------------------------------------
  //  Write headers in robustness files
  // --------------------------------------
  if (!full_output)
  {
    fprintf(output_file,"# ------------------------------------------------------\n");
    fprintf(output_file,"# Evolvability, Robustness and Antirobustness statistics\n");
    fprintf(output_file,"# ------------------------------------------------------\n");
    fprintf(output_file,"# \n");
    fprintf(output_file,"# 1. Generation \n");
    fprintf(output_file,"# 2. Fraction of positive offspring  (2*10 is Evolvability) \n");
    fprintf(output_file,"# 3. Fraction of neutral offspring (aka reproductive robustness) \n");
    fprintf(output_file,"# 4. Fraction of neutral mutants (aka mutational robustness) \n");
    fprintf(output_file,"# 5. Fraction of negative offspring \n");
    fprintf(output_file,"# 8. Cumul of delta-gaps of positive offspring\n");
    fprintf(output_file,"# 9. Cumul of delta-gaps of negative offspring\n");
    fprintf(output_file,"# 6. Delta-gap for the best offspring \n");
    fprintf(output_file,"# 7. Delta-gap for the worst offspring \n");
    fprintf(output_file,"# 10. Cumul of delta-fitness of positive offspring (2*10 is Evolvability)\n");
    fprintf(output_file,"# 11. Cumum of delta-fitness of negative offspring\n");
    fprintf(output_file,"# 12. Delta-fitness for the best offspring\n");
    fprintf(output_file,"# 13. Delta-fitness for the worst offspring\n\n\n");
  }
  else
  {
    fprintf(output_file,"# ------------------------------------------------------\n");
    fprintf(output_file,"# Evolvability, Robustness and Antirobustness statistics\n");
    fprintf(output_file,"# ------------------------------------------------------\n");
    fprintf(output_file,"# \n");
    fprintf(output_file,"# 1. Generation \n");
    fprintf(output_file,"# 2 to n+2. delta-fitness of each tested offspring \n\n\n");

  }
  return 0;
}