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
#include <signal.h>
#include <zlib.h>


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
  // Initialize command-line option variables with default values
  char* init_file_name  = NULL;

  // Define allowed options
  const char * options_list = "hVf:";
  static struct option long_options_list[] = {
    {"help",    no_argument,        NULL, 'h'},
    {"version", no_argument,        NULL, 'V'},
    { "from",   required_argument,  NULL, 'f' },
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
      {
        init_file_name = new char[strlen(optarg) + 1];
        sprintf(init_file_name, "%s", optarg);
        break;
      }
    }
  }


  if (init_file_name == NULL)
  {
    printf("You must specify an indiv file. Please use the option -f or --from.\n");
    exit(EXIT_FAILURE);
  }

  int64_t num_gener_ = 0;

  //~ // We create a new simulation
  //~ printf("Creating the new simulation.\n");
  ae_experiment* sim1 = new ae_experiment();
  //~ ae_param_loader* param_loader = new ae_param_loader("param.in");
  //~ sim1->load_params(param_loader, NULL);
  //~ delete param_loader;


  // We modify its individuals
  printf("Loading the individual.\n");
  gzFile init_file = gzopen(init_file_name, "r");
  ae_individual* indiv = new ae_individual(init_file);
  printf("Closing the file.\n");
  gzclose(init_file);

  // ae_list* list_indivs = NULL;//sim1->pop()->indivs();
  // int32_t nb_indivs = NULL;//sim1->pop()->nb_indivs();

  // ae_list_node*  indiv_node = list_indivs->first();
  // int i;
  // for (i=0;i<nb_indivs;i++)
  // {
  //   ae_individual* new_indiv = new ae_individual(*indiv);
  //   ae_individual* old_indiv = (ae_individual*) indiv_node->obj();
  //   delete old_indiv;
  //   indiv_node->set_obj(new_indiv);
  //   indiv_node = indiv_node->next();
  // }

  // if (ae_common::pop_structure == true)
  // {
  //   ae_list_node*  indiv_node = list_indivs->first();
  //   ae_grid_cell*** pop_grid_ = NULL;//sim1->pop()->pop_grid();
  //   for (int16_t x = 0 ; x < ae_common::grid_x ; x++)
  //   {
  //     for (int16_t y = 0 ; y < ae_common::grid_y ; y++)
  //     {
  //       pop_grid_[x][y]->set_individual((ae_individual*) indiv_node->obj());
  //       pop_grid_[x][y]->individual()->set_grid_cell(pop_grid_[x][y]);
  //       indiv_node = indiv_node->next();
  //     }
  //   }
  // }

  // Evaluate the new individuals
  ae_population* pop = NULL;
  /*sim1->pop()*/pop->evaluate_individuals(sim1->env());

  // Delete the backup files created by the ae_experiment constructor
  char backup_file_name[50];
  char best_indiv_file_name[50];

#ifdef __REGUL
  sprintf(backup_file_name, "backup/gen_" TIMESTEP_FORMAT ".rae", num_gener_);
  sprintf(best_indiv_file_name, "backup/best_" TIMESTEP_FORMAT ".rae", num_gener_);
#else
  sprintf(backup_file_name, "backup/gen_" TIMESTEP_FORMAT ".ae", num_gener_);
  sprintf(best_indiv_file_name, "backup/best_" TIMESTEP_FORMAT ".ae", num_gener_);
#endif

  remove(backup_file_name);
  remove(best_indiv_file_name);

  // And create new ones
  printf("Writting the backup file\n");
  sim1->write_backup();
  return EXIT_SUCCESS;

}

void print_help(char* prog_name)
{
  printf("\n************* aevol - Artificial Evolution ************* \n\n");
  printf("This program is Free Software. No Warranty.\n\n\
Usage : set_genome_from_indiv -h\n\
   or : set_genome_from_indiv -f source\n\n\
\t-h : Display this screen\n\
\t--from source    : Create a new simulation from param.in then load the individual from the file 'source'\n");
}
