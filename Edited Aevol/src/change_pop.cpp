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
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>



// =================================================================
//                            Project Files
// =================================================================
#include <f_line.h>
#include "ae_population.h"
#ifdef __X11
  #include "ae_exp_manager_X11.h"
#else
  #include "ae_exp_manager.h"
#endif

using namespace aevol;

// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_name);
f_line* get_line(FILE* param_file);
void format_line(f_line* formated_line, char* line, bool* line_is_interpretable);

int main(int argc, char* argv[])
{
  // 1) Initialize command-line option variables with default values
  char* param_file_name = NULL;
  char* pop_file_name   = NULL;
  bool verbose          = false;

  // 2) Define allowed options
  const char * options_list = "hp:f:";
  static struct option long_options_list[] = {
    { "help", no_argument,        NULL, 'h' },
    { "pop",  required_argument,  NULL, 'p' }, // Provide population file,
    { "file", required_argument,  NULL, 'f' }, // Provide file with population parameters to change
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
      case 'f' :
      {
        if (strcmp(optarg, "") == 0)
        {
          printf("%s: error: Option -f or --file : missing argument.\n", argv[0]);
          exit(EXIT_FAILURE);
        }

        param_file_name = optarg;
        break;
      }
      case 'p' :
      {
        if (strcmp(optarg, "") == 0)
        {
          printf("%s: error: Option -p or --pop : missing argument.\n", argv[0]);
          exit(EXIT_FAILURE);
        }

        pop_file_name = optarg;/*new char[strlen(optarg)+1];
        memcpy(pop_file_name, optarg, strlen(optarg)+1);*/

        break;
      }
      default :
      {
        // An error message is printed in getopt_long, we just need to exit
        exit(EXIT_FAILURE);
      }
    }
  }

  // 4) Check the consistancy of the command-line options
  if (param_file_name == NULL || pop_file_name == NULL)
  {
    printf("%s: error: You must provide both a parameter file and a population backup.\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  // 5) Initialize an empty experiment manager
  #ifndef __NO_X
    ae_exp_manager* exp_manager = new ae_exp_manager_X11();
  #else
    ae_exp_manager* exp_manager = new ae_exp_manager();
  #endif

  // 6) Initialize an empty population then load the backup population
  printf("Loading the backup population\t");
  ae_population* pop = new ae_population(exp_manager);
  gzFile pop_file = gzopen(pop_file_name, "r");
  if (pop_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, pop_file_name);
    exit(EXIT_FAILURE);
  }
  pop->load(pop_file, verbose);
  printf("Ok\n");

  // 7) Interpret and apply changes
  printf("Interpret and apply changes\n");
  FILE* param_file  = fopen(param_file_name,  "r");
  if (param_file == NULL)
  {
    printf("%s:%d: error: could not open parameter file %s\n", __FILE__, __LINE__, param_file_name);
    exit(EXIT_FAILURE);
  }

  f_line* line;
  while ((line = line(param_file)) != NULL)
  {
    if (strcmp(line->words[0], "POINT_MUTATION_RATE") == 0)
    {
      pop->set_overall_point_mutation_rate(atof(line->words[1]));
      printf("\tChange of overall point mutation rate to %f\n",atof(line->words[1]));
    }
    else if (strcmp(line->words[0], "SMALL_INSERTION_RATE") == 0)
    {
      pop->set_overall_small_insertion_rate(atof(line->words[1]));
      printf("\tChange of overall small insertion rate to %f\n",atof(line->words[1]));
    }
    else if (strcmp(line->words[0], "SMALL_DELETION_RATE") == 0)
    {
      pop->set_overall_small_deletion_rate(atof(line->words[1]));
      printf("\tChange of overall small deletion rate to %f\n",atof(line->words[1]));
    }
    else if (strcmp(line->words[0], "MAX_INDEL_SIZE") == 0)
    {
      pop->set_overall_max_indel_size(atol(line->words[1]));
      printf("\tChange of overall maximum indel size to %f\n",atof(line->words[1]));
    }
    else if (strcmp(line->words[0], "DUPLICATION_RATE") == 0)
    {
      pop->set_overall_duplication_rate(atof(line->words[1]));
      printf("\tChange of overall duplication rate to %f\n",atof(line->words[1]));
    }
    else if (strcmp(line->words[0], "DELETION_RATE") == 0)
    {
      pop->set_overall_deletion_rate(atof(line->words[1]));
      printf("\tChange of overall deletion rate to %f\n",atof(line->words[1]));
    }
    else if (strcmp(line->words[0], "TRANSLOCATION_RATE") == 0)
    {
      pop->set_overall_translocation_rate(atof(line->words[1]));
      printf("\tChange of overall translocation rate to %f\n",atof(line->words[1]));
    }
    else if (strcmp(line->words[0], "INVERSION_RATE") == 0)
    {
      pop->set_overall_inversion_rate(atof(line->words[1]));
      printf("\tChange of overall inversion to %f\n",atof(line->words[1]));
    }

    delete line;
  }
  fclose(param_file);
  printf("Ok\n");

  // 8) Save the changed population in a new population file (similar name with _changed)
  printf("Save the changed population\t");
  char* new_pop_file_name   = NULL;
  new_pop_file_name = new char[strlen(pop_file_name)-3+strlen("_changed.ae")+1];
  sprintf(new_pop_file_name, "%s_changed.ae",strtok(pop_file_name, "."));
  gzFile new_pop_file = gzopen(new_pop_file_name, "w");
  if (new_pop_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, new_pop_file_name);
    exit(EXIT_FAILURE);
  }
  pop->save(new_pop_file);
  gzclose(new_pop_file);
  printf("Ok\n");
}



f_line* get_line(FILE* param_file)
{
  char line[255];
  f_line* formated_line = new f_line();

  bool found_interpretable_line = false;

  while (!feof(param_file) && !found_interpretable_line)
  {
    if (!fgets(line, 255, param_file))
    {
      delete formated_line;
      return NULL;
    }
    format_line(formated_line, line, &found_interpretable_line);
  }

  if (found_interpretable_line)
  {
    return formated_line;
  }
  else
  {
    delete formated_line;
    return NULL;
  }
}

void format_line(f_line* formated_line, char* line, bool* line_is_interpretable)
{
  int16_t i = 0;
  int16_t j;

  // Parse line
  while (line[i] != '\n' && line[i] != '\0' && line[i] != '\r')
  {
    j = 0;

    // Flush white spaces and tabs
    while (line[i] == ' ' || line[i] == 0x09) i++; // 0x09 is the ASCII code for TAB

    // Check comments
    if (line[i] == '#') break;

    // If we got this far, there is content in the line
    *line_is_interpretable = true;

    // Parse word
    while (line[i] != ' '  && line[i] != '\n' && line[i] != '\0' && line[i] != '\r')
    {
      formated_line->words[formated_line->nb_words][j++] = line[i++];
    }

    // Add '\0' at end of word if it's not empty (line ending with space or tab)
    if (j != 0)
    {
      formated_line->words[formated_line->nb_words++][j] = '\0';
    }
  }
}


void print_help(char* prog_name)
{
	printf("******************************************************************************\n");
	printf("*                        aevol - Artificial Evolution                        *\n");
	printf("******************************************************************************\n");
	printf("Usage : change_pop -h\n");
	printf("   or : change_pop [-f param_file -p pop_file]\n");
	printf("  -h, --help  Display this screen\n");
	printf("  -p, --pop  Specify population to change\n");
	printf("  -f, --file  Specify file with parameters to change\n");
	printf("Change a population file as specified in the parameter file.\n");
    printf("(default: param_to_change.in)\n\n");
}
