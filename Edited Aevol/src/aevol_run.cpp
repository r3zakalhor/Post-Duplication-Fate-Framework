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
#include <cinttypes>
#include <cstdlib>
#include <cstdio>
#include <csignal>

#include <sys/stat.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <sys/stat.h>

#ifdef _OPENMP
#ifndef __OPENMP_GPU
#include <omp.h>
#endif
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num() { return 0;}
inline omp_int_t omp_get_max_threads() { return 1;}
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

// =================================================================
//                            Project Files
// =================================================================
#ifdef __X11
  #include "ExpManager_X11.h"
#else
  #include "ExpManager.h"
#endif

#include "ExpManager_7.h"
#include "macros.h"

#include "aevol_version.h"

using namespace aevol;

#ifndef __NO_X
  void catch_usr1(int sig_num);
#endif

// Helper functions
void print_help(char* prog_path);
void interpret_cmd_line_options(int argc, char* argv[]);

// Command-line option variables
// static bool pause_on_startup = false;
static bool verbose = false;
static int64_t t0 = -1;
static int64_t t_end = -1;
static int64_t nb_steps = -1;
static int grain_size = 1;
static bool w_mrca = false;

#ifndef __NO_X
  static bool show_display_on_startup = true;
#endif

#ifdef _OPENMP
static bool run_in_parallel = false;
#endif

// Other file-scope variables
static ExpManager* exp_manager = nullptr;


int main(int argc, char* argv[]) {
  // Set handlers for signals to be caught
  #ifndef __NO_X
    signal(SIGUSR1, catch_usr1);
  #endif

  ExpManager_7::standalone_simd = true;
  
  interpret_cmd_line_options(argc, argv);

  // Print the hash and date of the commit used to compile this executable
  printf("Running aevol version %s\n", version_string);

  #ifdef HAVE_MPI
  printf("Initializing MPI Environment\n");
  MPI_Init(NULL, NULL);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  #endif

  // =================================================================
  //                          Load the simulation
  // =================================================================
  #ifndef __NO_X
    exp_manager = new ExpManager_X11();
  #else
    exp_manager = new ExpManager();
  #endif

  exp_manager->load(t0, verbose, true
  #ifdef HAVE_MPI
  ,world_rank
  #endif
  );
  exp_manager->set_t_end(t_end);

  #ifdef HAVE_MPI
  if (world_size != exp_manager->exp_s()->nb_rank()) {
    printf("MANDATORY : the nb_rank specify through aevol_create should be the same than the one for aevol_run\n");
    exit(-1);
  }
  #endif


  // Make a numbered copy of each static input file
  // (dynamic files are saved elsewhere)
  // TODO (?)

  #ifndef __NO_X
    if (show_display_on_startup) {
      ((ExpManager_X11*) exp_manager)->toggle_display_on_off();
    }
  #endif

  // =================================================================
  //                         Run the simulation
  // =================================================================
  exp_manager->run_evolution();
  
  #ifdef HAVE_MPI
  MPI_Finalize();
  #endif

  delete exp_manager;
  //~ return EXIT_SUCCESS;
}


#ifndef __NO_X
  void catch_usr1(int sig_num) {
    signal(SIGUSR1, catch_usr1);

    printf("display on/off\n");
    if (exp_manager != nullptr) {
      ((ExpManager_X11*) exp_manager)->toggle_display_on_off();
    }
  }
#endif


/*!
  \brief

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
	printf("%s: run an aevol simulation.\n", prog_name);
  printf("\n");
	printf("Usage : %s -h or --help\n", prog_name);
	printf("   or : %s -V or --version\n", prog_name);
	printf("   or : %s [-b TIMESTEP] [-e TIMESTEP|-n NB_TIMESTEPS] [-p NB_THREADS] [-vwx]\n",
         prog_name);
	printf("\nOptions\n");
	printf("  -h, --help\n\tprint this help, then exit\n");
	printf("  -V, --version\n\tprint version number, then exit\n");
  printf("  -b, --begin TIMESTEP\n");
  printf("\tspecify time t0 to resume simulation at (default read in last_gener.txt)\n");
  printf("  -e, --end TIMESTEP\n");
  printf("\tspecify time of the end of the simulation\n");
  printf("  -n, --nb-timesteps NB_TIMESTEPS\n");
  printf("\tspecify number of timesteps to be simulated (default 1000)\n");
  printf("  -p, --parallel NB_THREADS\n");
  printf("\trun on NB_THREADS threads (use -1 for system default)\n");
  printf("  -v, --verbose\n\tbe verbose\n");
  printf("  -w, --wait\n\tpause after loading\n");
  printf("  -x, --noX\n\tdon't display X outputs upon start\n");
  printf("\tsend SIGUSR1 to switch X output on/off\n");
}

void interpret_cmd_line_options(int argc, char* argv[]) {
  // Define allowed options
  const char* options_list = "hVb:e:n:vwxp:6:";
  static struct option long_options_list[] = {
      {"help",          no_argument,       nullptr, 'h'},
      {"version",       no_argument,       nullptr, 'V'},
      {"begin",         required_argument, nullptr, 'b'},
      {"end",           required_argument, nullptr, 'e'},
      {"nb-timesteps",  required_argument, nullptr, 'n'},
      {"verbose",       no_argument,       nullptr, 'v'},
      {"wait",          no_argument,       nullptr, 'w'},
      {"noX",           no_argument,       nullptr, 'x'},
      {"parallel",      required_argument, nullptr, 'p'},
      {"aevol_6",       no_argument,       nullptr, '6'},
      {0, 0, 0, 0}
  };

  // Get actual values of the CLI options
  int option;
  while ((option =
              getopt_long(argc, argv, options_list, long_options_list, nullptr))
         != -1) {
    switch (option) {
      case 'h' : {
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      }
      case 'V' : {
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      }
      case 'b' : {
        t0 = atol(optarg);
        break;
      }
      case 'e' : {
        if (nb_steps != -1) {
          Utils::ExitWithUsrMsg("use either option -n or -e, not both");
        }

        t_end = atol(optarg);
        w_mrca = true;
        break;
      }
      case 'n' : {
        if (t_end != -1) {
          Utils::ExitWithUsrMsg("use either option -n or -e, not both");
        }

        nb_steps = atol(optarg);
        break;
      }
      case 'v' : {
        verbose = true;
        break;
      }
      case 'w' : {
        // pause_on_startup = true;
        break;
      }
    case '6' : {
      ExpManager_7::standalone_simd = false;
      break;
    }
      case 'x' : {
        #ifdef __NO_X
        printf("%s: error: Program was compiled with __NO_X option, "
                 "no visualisation available.\n", argv[0]);
          exit(EXIT_FAILURE);
        #else
        show_display_on_startup = false;
        #endif

        break;
      }
      case 'p' : {
        #ifdef _OPENMP
        run_in_parallel = true;
          if (atoi(optarg) > 0) {
            omp_set_num_threads(atoi(optarg));
          }
        #endif
        break;
      }
      case 'g' : {
        grain_size = atoi(optarg);
        break;
      }
      default : {
        // An error message is printed in getopt_long, we just need to exit
        exit(EXIT_FAILURE);
      }
    }
  }

  // If t0 wasn't provided, use default
  if (t0 < 0) {
    t0 = OutputManager::last_gener();
  }

  // If t_end_ wasn't provided, set it according to nb_steps or use default
  // (run for 1000 timesteps)
  if (t_end < 0) {
    if (nb_steps >= 0) {
      t_end = t0 + nb_steps;
    }
    else {
      t_end = t0 + 1000;
    }
  }

  // It the user didn't ask for a parallel run, set number of threads to 1
  #ifdef _OPENMP
  if (not run_in_parallel) {
      omp_set_num_threads(1);
    }
  #endif
}
