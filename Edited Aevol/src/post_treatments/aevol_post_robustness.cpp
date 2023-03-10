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
#include <cerrno>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cinttypes>
#include <cmath>

#include <getopt.h>
#include <zlib.h>
#include <inttypes.h>
#include <getopt.h>
#include <sys/stat.h>

#include <list>

// =================================================================
//                            Project Files
// =================================================================
#include "aevol.h"
#include "IndivAnalysis.h"

using std::list;
using namespace aevol;

// Helper functions
void print_help(char* prog_path);
void interpret_cmd_line_options(int argc, char* argv[]);

// Command-line option variables
static int32_t nb_children = 1000;
static int32_t wanted_rank = -1;
static int32_t wanted_index = -1;
static int64_t timestep = -1;

int main(int argc, char* argv[]) {
  interpret_cmd_line_options(argc, argv);

  // ----------------------
  //  Prepare the outputs
  // ----------------------
  char directory_name[255];
  snprintf(directory_name, 255, "analysis-generation_" TIMESTEP_FORMAT,
           timestep);
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

  // -----------------------------------------------------------------------------------
  //  Load the backup and get the individual for which detailed information is desired
  // -----------------------------------------------------------------------------------
  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(timestep, true, false);

  Individual* indiv_tmp = nullptr;

  // If neither a rank nor an index was provided, consider the best indiv
  if ((wanted_rank == -1) && (wanted_index == -1)) {
    indiv_tmp = exp_manager->best_indiv();
  }
  else if (wanted_index != -1) {
    indiv_tmp = exp_manager->indiv_by_id(wanted_index);
  }
  else {
    indiv_tmp = exp_manager->indiv_by_rank(wanted_rank);
  }
  // Update id and rank
  wanted_index = indiv_tmp->id();
  wanted_rank = indiv_tmp->rank();

  IndivAnalysis wanted_indiv(*indiv_tmp);
  wanted_indiv.set_grid_cell(indiv_tmp->grid_cell());
  wanted_indiv.grid_cell()->set_individual((Individual*)&wanted_indiv);;

  // Now that we have the index and rank of the indiv of interest, we can
  // generate the output file name and hence open that file
  char filename[255];
  snprintf(filename, 255, "%s/robustness-summary-%" PRId64 "-i%" PRId32 "-r%" PRId32,
           directory_name, timestep, wanted_index, wanted_rank);
  FILE* output_summary = fopen(filename, "w");
  if (output_summary == nullptr) {
    Utils::ExitWithUsrMsg(std::string("Could not open file ") + filename);
  }
  snprintf(filename, 255, "%s/robustness-detailed-%" PRId64 "-i%" PRId32 "-r%" PRId32,
           directory_name, timestep, wanted_index, wanted_rank);
  FILE* output_detailed = fopen(filename, "w");
  if (output_detailed == nullptr) {
    Utils::ExitWithUsrMsg(std::string("Could not open file ") + filename);
  }

  fprintf(output_summary, "###############################################################################\n");
  fprintf(output_summary, "#  Summary of the mutants generated from individual with rank %" PRId32
                          " and index %" PRId32 " at timestep %" PRId64 " \n",
          wanted_rank, wanted_index, timestep);
  fprintf(output_summary, "###############################################################################\n");
  fprintf(output_summary, "#  1.  Timestep\n");
  fprintf(output_summary, "#  2.  Proportion of mutants that are better than their parent\n");
  fprintf(output_summary, "#  3.  Proportion of mutants that are as good as their parent was\n");
  fprintf(output_summary, "#  4.  Proportion of mutants that are worse than their parent\n");
  fprintf(output_summary, "#  5.  Average difference in metabolic error btw good mutants and their parent\n");
  fprintf(output_summary, "#  6.  Average difference in metabolic error btw bad mutants and their parent\n");
  fprintf(output_summary, "#  7.  Maximum gain in metabolic error among good mutants\n");
  fprintf(output_summary, "#  8.  Maximum loss of metabolic error among bad mutants\n");
  fprintf(output_summary, "###############################################################################\n");

  fprintf(output_detailed, "###############################################################################\n");
  fprintf(output_detailed, "#  Mutants generated from individual with rank %" PRId32
                           " and index %" PRId32 " at timestep %" PRId64 " \n",
          wanted_rank, wanted_index, timestep);
  fprintf(output_detailed, "###############################################################################\n");
  fprintf(output_detailed, "#  1.  Parent id\n");
  fprintf(output_detailed, "#  2.  Parent metabolic error\n");
  fprintf(output_detailed, "#  3.  Parent secretion\n");
  fprintf(output_detailed, "#  4.  Mutant metabolic error\n");
  fprintf(output_detailed, "#  5.  Mutant secretion\n");
  fprintf(output_detailed, "#  6.  Mutant genome size\n");
  fprintf(output_detailed, "#  7.  Mutant number of functional genes\n");
  fprintf(output_detailed, "###############################################################################\n");

  wanted_indiv.Evaluate();
  wanted_indiv.compute_statistical_data();
  wanted_indiv.compute_non_coding();

  wanted_indiv.compute_experimental_f_nu(
      nb_children,
      std::make_shared<JumpingMT>(time(nullptr)),
      output_summary,
      output_detailed);

  indiv_tmp->set_grid_cell(wanted_indiv.grid_cell());
  indiv_tmp->grid_cell()->set_individual(indiv_tmp);

  fclose(output_summary);
  fclose(output_detailed);
  delete exp_manager;


  return EXIT_SUCCESS;
}


/**
 * \brief print help and exist
 */
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
  printf("%s:\n", prog_name);
  printf("\tComputes replication statistics of a given individual\n");
  printf("\tThis is achieved by performing NB_MUTANTS replications\n");
  printf("\tof the individual of interest and evaluating them\n");
  printf("\n");
  printf("Usage : %s -h or --help\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("   or : %s [-t TIMESTEP] [-I INDEX | -R RANK] [-n NB_MUTANTS]\n",
         prog_name);
  printf("\nOptions\n");
  printf("  -h, --help\n\tprint this help, then exit\n");
  printf("  -V, --version\n\tprint version number, then exit\n");
  printf("  -t, --timestep TIMESTEP\n");
  printf("\tspecify timestep of the individual of interest\n");
  printf("  -I, --index INDEX\n");
  printf("\tspecify the index of the individual of interest\n");
  printf("  -R, --rank RANK\n");
  printf("\tspecify the rank of the individual of interest\n");
  printf("  -n, --nb-mutants NB_MUTANTS\n");
  printf("\tspecify the number of mutants to be generated\n");
}

void interpret_cmd_line_options(int argc, char* argv[]) {
  const char* options_list = "hVt:n:R:I:";
  static struct option long_options_list[] = {
      {"help",       no_argument,       NULL, 'h'},
      {"version",    no_argument,       NULL, 'V'},
      {"timestep",   required_argument, NULL, 't'},
      {"nb-mutants", required_argument, NULL, 'n'},
      {"rank",       required_argument, NULL, 'R'},
      {"index",      required_argument, NULL, 'I'},
      {0, 0, 0, 0}
  };

  int option;
  while ((option = getopt_long(argc, argv, options_list, long_options_list,
                               NULL)) != -1) {
    switch (option) {
    case 'h' :
      print_help(argv[0]);
      exit(EXIT_SUCCESS);
    case 'V' :
      Utils::PrintAevolVersion();
      exit(EXIT_SUCCESS);
    case 't' :
      timestep = atol(optarg);
      break;
    case 'n' :
      nb_children = atol(optarg);
      break;
    case 'R' :
      if (wanted_index != -1) {
        Utils::ExitWithUsrMsg("Options -R and -I are incompatible");
      }
      wanted_rank = atol(optarg);
      break;
    case 'I' :
      if (wanted_rank != -1) {
        Utils::ExitWithUsrMsg("Options -R and -I are incompatible");
      }
      wanted_index = atol(optarg);
      break;
    default :
      // An error message is printed in getopt_long, we just need to exit
      exit(EXIT_FAILURE);
    }
  }

  // If timestep wasn't provided, use default
  if (timestep < 0) {
    timestep = OutputManager::last_gener();
  }
}
