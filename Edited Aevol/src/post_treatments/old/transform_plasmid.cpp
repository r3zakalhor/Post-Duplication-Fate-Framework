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

// This program transforms a plasmid (read as a text file) into each individual.

// =================================================================
//                              Libraries
// =================================================================
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>

// =================================================================
//                            Project Files
// =================================================================
#include "aevol.h"


using namespace aevol;

void print_help(char* prog_name);

int main(int argc, char* argv[])
{
  // Initialize command-line option variables with default values
  char* plasmid_file_name  = NULL;
  int32_t num_gener = -1;
  int32_t plasmid_maximal_length = 10000000;
  int32_t plasmid_minimal_length = 40;

  // Define allowed options
  // TODO <david.parsons@inria.fr> version
  const char * options_list = "hr:p:m:M:";
  static struct option long_options_list[] = {
    { "help", 1, NULL, 'h' },
    { "resume", 1, NULL, 'r' },
    { "plasmid", 1, NULL, 'p' },
    { "min", 1, NULL, 'm' },
    { "max", 1, NULL, 'M' },
    { 0, 0, 0, 0 }
  };

  // Get actual values of the command-line options
  int option;
  while ((option = getopt_long(argc, argv, options_list, long_options_list, NULL)) != -1)
  {
    switch (option)
    {
      case 'h' :
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
        break;
      case 'r':
        num_gener = atol(optarg);
        break;
      case 'p':
        plasmid_file_name = new char [strlen(optarg)+1];
        strcpy(plasmid_file_name, optarg);
        break;
      case 'm':
        plasmid_minimal_length = atol(optarg);
        break;
      case 'M':
        plasmid_maximal_length = atol(optarg);
        break;
    }
  }

  // Open the files
  ae_exp_manager* exp_manager = new ae_exp_manager();

  // We need a full backup
  if (num_gener == -1)
  {
    printf("You must specify a generation number");
    exit(EXIT_FAILURE);
  }
  exp_manager->load(num_gener, false, false, false);

  // And a plasmid
  if (plasmid_file_name==NULL)
  {
    printf("You must specify a plasmid");
    exit(EXIT_FAILURE);
  }

  char rawplasmid[1000000];
  int32_t lplasmid=-1;
  FILE* plasmid_file = fopen(plasmid_file_name,"r");
  if (fgets(rawplasmid, 1000000, plasmid_file) == nullptr) {
    printf("Failed to read from %s\n", plasmid_file_name);
    exit(EXIT_FAILURE);
  }
  lplasmid = strlen(rawplasmid)-1;
  fclose(plasmid_file);

  // We know need to make several changes in experiment so it 'accepts' this new plasmid
  exp_manager->exp_s()->set_with_plasmids(true); // We need to change the allow_plasmid parameter, otherwise aevol_run will crash
  exp_manager->output_m()->set_compute_phen_contrib_by_GU(true);

  // We parse the individuals and transform them
  for (const auto& indiv: exp_manager->world()->indivs()) {
    char* plasmid=new char[lplasmid+1]; // Warning: will become the DNA of the first individual created -> no not delete, will be freed in ~ae_dna.
    strncpy(plasmid, rawplasmid, lplasmid);
    plasmid[lplasmid]='\0';
    indiv->add_GU(plasmid,lplasmid);
    indiv->set_allow_plasmids(true);
    indiv->genetic_unit_nonconst(1).set_min_gu_length(plasmid_minimal_length);
    indiv->genetic_unit_nonconst(1).set_max_gu_length(plasmid_maximal_length);
    indiv->set_replication_report(NULL); // plasmid's DNA should not have replic reports otherwise stat_record will try to access it.
    indiv->genetic_unit_nonconst(1).dna()->set_replic_report(NULL);
    indiv->genetic_unit_nonconst(1).compute_phenotypic_contribution();
    indiv->reevaluate();
  }

  // We save and exit
  exp_manager->write_setup_files();
  exp_manager->save();
  delete exp_manager;
  delete [] plasmid_file_name;

  return EXIT_SUCCESS;
}

void print_help(char* prog_name)
{
  printf("\n");
  printf("Add a plasmid to each individual \n");
  printf("Usage : transform_indiv -h\n");
  printf("or :    transform_indiv -r ng -p plasmid_file \n");
  printf("\t-h : display this screen\n");
  printf("\t-r ng  : modify generation ng (need a full backup), \n");
  printf("\t-p plasmid_file : read plasmid sequence from file plasmid_file\n");
}
