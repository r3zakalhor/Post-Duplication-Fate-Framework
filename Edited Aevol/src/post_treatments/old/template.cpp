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
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <getopt.h>



// =================================================================
//                            Project Files
// =================================================================
#include "aevol.h"

using namespace aevol;

// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);





int main(int argc, char* argv[])
{
  // 1) Initialize command-line option variables with default values
  // bool  opt = false;
  char* input_dir = NULL;
  int32_t num_gener = -1;


  // 2) Define allowed options
  const char * options_list = "hVoi:g:";
  static struct option long_options_list[] = {
    { "help",     no_argument,        NULL, 'h' },
    { "version",  no_argument,        NULL, 'V' },
    { "opt",      no_argument,        NULL, 'o' },
    { "in",       required_argument,  NULL, 'i' },
    { "gener",    required_argument,  NULL, 'g' },
    { 0, 0, 0, 0 }
  };


  // 3) Get actual values of the command-line options
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
      case 'o' :
      {
        // opt = true;
        break;
      }
      case 'i' :
      {
        input_dir = new char[strlen(optarg)+1];
        strcpy(input_dir, optarg);
        break;
      }
      case 'g' :
      {
        num_gener = atol(optarg);
        break;
      }
      default :
      {
        // An error message is printed in getopt_long, we just need to exit
        exit(EXIT_FAILURE);
      }
    }
  }

  // 4) Check for missing mandatory arguments
  if (input_dir == NULL)
  {
    printf("%s: error: You must provide an input directory.\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  if (num_gener == -1)
  {
    printf("%s: error: You must provide a generation number.\n", argv[0]);
    exit(EXIT_FAILURE);
  }


  // Load the experiment
  ae_exp_manager* exp = new ae_exp_manager();
  exp->load(input_dir, num_gener, false, false, false);


  // *********************************************************************** //
  // Do something here
  // *********************************************************************** //

  delete exp;
}


/*!
  \brief

*/
void print_help(char* prog_path)
{
  // Get the program file-name in prog_name (strip prog_path of the path)
  char* prog_name; // No new, it will point to somewhere inside prog_path
  if ((prog_name = strrchr(prog_path, '/'))) prog_name++;
  else prog_name = prog_path;

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
  printf("%s: does nothing (template file).\n", prog_name);
  printf("\n");
  printf("Usage : %s -h or --help\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("   or : %s [-o] -a arg\n", prog_name);
  printf("\nOptions\n");
  printf("  -h, --help\n\tprint this help, then exit\n\n");
  printf("  -V, --version\n\tprint version number, then exit\n\n");
  printf("  -o, --opt\n\toption with no argument\n");
  printf("  -a, --arg argument\n\toption with an argument\n");
}
