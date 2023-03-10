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
#include <cstdlib>
#include <cstdio>

#include <getopt.h>

#ifdef __NO_X
  #error This program requires graphics libraries
#else
  #include <X11/Xlib.h>
#endif

#include "aevol.h"

using namespace aevol;

// Helper functions
void print_help(char* prog_path);
void interpret_cmd_line_options(int argc, char* argv[]);

// Command-line option variables
static int64_t timestep = -1;


int main(int argc, char* argv[]) {
  interpret_cmd_line_options(argc, argv);

  // If timestep wasn't provided, use default
  if (timestep < 0) {
    timestep = OutputManager::last_gener();
  }

  printf("Displaying timestep %" PRId64 "...\n", timestep);

  // =================================================================
  //                       Read the backup file
  // =================================================================
  // Load simulation from backup
  ExpManager_X11* exp_manager = new ExpManager_X11();
  exp_manager->load(timestep, false, false);

  // =================================================================
  //                       Draw the windows
  // =================================================================
  // Display is off by default, switch it on
  exp_manager->toggle_display_on_off();

  // Display is usually triggered in ExpManager::run_evolution(), here we want
  // to call it manually
  exp_manager->display();

  // Handle user events until he quits
  while (not exp_manager->quit_signal_received()) {
    exp_manager->handle_events();
  }

  delete exp_manager;
  return EXIT_SUCCESS;
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
  printf("%s: view the simulation at the provided timestep\n",
         prog_name);
  printf("\n");
  printf("Usage : %s -h or --help\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("   or : %s [-t TIMESTEP]\n",
         prog_name);
  printf("\nOptions\n");
  printf("  -h, --help\n\tprint this help, then exit\n\n");
  printf("  -V, --version\n\tprint version number, then exit\n\n");
  printf("  -t, --timestep TIMESTEP\n");
  printf("\tspecify timestep to display (default value read in last_gener.txt)\n");
}

void interpret_cmd_line_options(int argc, char* argv[]) {
  // Define allowed options
  const char* options_list = "hVt:";
  static struct option long_options_list[] = {
      {"help",       no_argument,       nullptr, 'h'},
      {"version",    no_argument,       nullptr, 'V'},
      {"timestep",   required_argument, nullptr, 't'},
      {0, 0, 0, 0}
  };

  // Get actual values of the CLI options
  int option;
  while ((option = getopt_long(argc, argv, options_list, long_options_list,
                               nullptr)) != -1) {
    switch (option) {
      case 'h' : {
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      }
      case 'V' : {
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      }
      case 't' : {
        timestep = atol(optarg);
        break;
      }
      default : {
        // An error message is printed in getopt_long, we just need to exit
        exit(EXIT_FAILURE);
      }
    }
  }
}
