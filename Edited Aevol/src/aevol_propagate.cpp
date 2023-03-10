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


const char* DEFAULT_PARAM_FILE_NAME = "param.in";

// ============================================================================
//                                   Includes
// ============================================================================
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <err.h>
#include <errno.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#include <limits>

// TODO(dpa) will be merged into C++17, remove dependency when possible
#include <boost/filesystem.hpp>


// =================================================================
//                            Project Files
// =================================================================
#if __cplusplus == 201103L
  #include "make_unique.h"
#endif

#ifdef __X11
  #include "ExpManager_X11.h"
#else
  #include "ExpManager.h"
#endif

#include "ParamLoader.h"

using namespace aevol;

// Helper functions
void print_help(char* prog_path);
void interpret_cmd_line_options(int argc, char* argv[]);

// Command-line option variables
static int64_t timestep = -1;
static char* input_dir = nullptr;
static char* output_dir = nullptr;
static bool keep_prng_states = false;
static bool verbose = false;

static int32_t selseed        = -1;

int main(int argc, char* argv[]) {
  interpret_cmd_line_options(argc, argv);

  // If output directory doesn't exist, create it
  struct stat stat_buf;
  if ((stat(output_dir, &stat_buf) == -1) && (errno == ENOENT)) {
    if (not boost::filesystem::create_directories(output_dir)) {
      err(EXIT_FAILURE, output_dir, errno);
    }
  }



  // =================================================================
  //                    Load the model experiment
  // =================================================================
#ifndef __NO_X
  ExpManager *exp_manager = new ExpManager_X11();
#else
  ExpManager* exp_manager = new ExpManager();
#endif

  exp_manager->load(input_dir, timestep, verbose, true);

  if (not keep_prng_states) {

      auto max = std::numeric_limits<int32_t>::max();

#if __cplusplus == 201103L
    auto prng = make_unique<JumpingMT>(time(nullptr));

    exp_manager->sel()->set_prng(
        make_unique<JumpingMT>(prng->random(max)));
    exp_manager->world()->set_prng(
        make_unique<JumpingMT>(prng->random(max)));
#else
    auto prng = std::make_unique<JumpingMT>(time(nullptr));

    for (int16_t x = 0; x < exp_manager->world()->width(); x++) {
      for (int16_t y = 0; y < exp_manager->world()->height(); y++) {
        int32_t seed = prng->random(1000000);
#if __cplusplus == 201103L
        exp_manager->world()->grid(x,y)->set_reprod_prng(make_unique<JumpingMT>(seed));
          exp_manager->world()->grid(x,y)->set_reprod_prng_simd(make_unique<JumpingMT>(seed));
          exp_manager->world()->grid(x,y)->set_mut_prng(std::make_shared<JumpingMT>(seed));
          exp_manager->world()->grid(x,y)->set_stoch_prng(std::make_shared<JumpingMT>(seed));
#else
        exp_manager->world()->grid(x, y)->set_reprod_prng(std::make_unique<JumpingMT>(seed));
        exp_manager->world()->grid(x, y)->set_reprod_prng_simd(std::make_unique<JumpingMT>(seed));
        exp_manager->world()->grid(x, y)->set_mut_prng(std::make_shared<JumpingMT>(seed));
        exp_manager->world()->grid(x, y)->set_stoch_prng(std::make_shared<JumpingMT>(seed));
#endif
      }
    }
#endif


    exp_manager->world()->set_phen_target_prngs(
            std::make_shared<JumpingMT>(prng->random(max)),
            std::make_shared<JumpingMT>(prng->random(max)));

      exp_manager->save_copy(output_dir, timestep);

  } else {

    if (selseed != -1) {
#if __cplusplus == 201103L
      auto prng = make_unique<JumpingMT>(selseed);
#else
      auto prng = std::make_unique<JumpingMT>(selseed);
#endif

      for (int16_t x = 0; x < exp_manager->world()->width(); x++) {
        for (int16_t y = 0; y < exp_manager->world()->height(); y++) {
          int32_t seed = prng->random(1000000);
#if __cplusplus == 201103L
          exp_manager->world()->grid(x,y)->set_reprod_prng(make_unique<JumpingMT>(seed));
          exp_manager->world()->grid(x,y)->set_reprod_prng_simd(make_unique<JumpingMT>(seed));
#else
          exp_manager->world()->grid(x, y)->set_reprod_prng(std::make_unique<JumpingMT>(seed));
          exp_manager->world()->grid(x, y)->set_reprod_prng_simd(std::make_unique<JumpingMT>(seed));
#endif
        }
      }

    }


      exp_manager->save_copy(output_dir, timestep);
  }
}







/*!
  \brief print help and exist

*/
    void print_help(char *prog_path) {
      // Get the program file-name in prog_name (strip prog_path of the path)
      char *prog_name; // No new, it will point to somewhere inside prog_path
      if ((prog_name = strrchr(prog_path, '/'))) {
        prog_name++;
      } else {
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
      printf("\tCreate a fresh copy of the experiment as it was at the given timestep.\n");
      printf("\tThe timestep number of the copy will be reset to 0.\n");
      printf("\n");
      printf("Usage : %s -h or --help\n", prog_name);
      printf("   or : %s -V or --version\n", prog_name);
      printf("   or : %s [-t TIMESTEP] [-K] [-o OUTDIR] [-v]\n",
             prog_name);
      printf("\nOptions\n");
      printf("  -h, --help\n\tprint this help, then exit\n");
      printf("  -V, --version\n\tprint version number, then exit\n");
      printf("  -t, --timestep TIMESTEP\n");
      printf("\tspecify timestep to propagate\n");
      printf("  -K, --keep-prng-st\n\tdo not alter prng states\n");
//  printf("  -i, --in INDIR\n"
//             "\tspecify input directory (default \".\")\n");
      printf("  -o, --out OUTDIR\n"
             "\tspecify output directory (default \"./output\")\n");
      printf("  -v, --verbose\n\tbe verbose\n");
    }

    void interpret_cmd_line_options(int argc, char *argv[]) {
      // Define allowed options
      const char *options_list = "hVt:o:v";
      static struct option long_options_list[] = {
              {"help",         no_argument,       nullptr, 'h'},
              {"version",      no_argument,       nullptr, 'V'},
              {"timestep",     required_argument, nullptr, 't'},
              {"keep-prng-st", no_argument,       nullptr, 'K'},
//      {"in",             required_argument, nullptr, 'i'},
              {"out",          required_argument, nullptr, 'o'},
              {"verbose",      no_argument,       nullptr, 'v'},
              {0, 0,                              0,       0}
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
            timestep = atoi(optarg);
            break;
          }
          case 'K' : {
            keep_prng_states = true;
            break;
          }
//      case 'i' : {
//        input_dir = new char[strlen(optarg) + 1];
//        strcpy(input_dir, optarg);
//        break;
//      }
          case 'o' : {
            output_dir = new char[strlen(optarg) + 1];
            strcpy(output_dir, optarg);
            break;
          }
          case 'v' : {
            verbose = true;
            break;
          }
          default : {
            // An error message is printed in getopt_long, we just need to exit
            exit(EXIT_FAILURE);
          }
        }
      }

      // If input directory wasn't provided, use default
      if (input_dir == nullptr) {
        input_dir = new char[255];
        sprintf(input_dir, "%s", ".");
      }
      // If output directory wasn't provided, use default
      if (output_dir == nullptr) {
        output_dir = new char[255];
        sprintf(output_dir, "%s", "output");
      }

      // If timestep wasn't provided, use default
      if (timestep == -1) {
        timestep = OutputManager::last_gener();
      }
    }

