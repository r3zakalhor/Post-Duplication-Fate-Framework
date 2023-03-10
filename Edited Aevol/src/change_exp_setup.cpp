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
#include <sys/stat.h>
#include <errno.h>


// =================================================================
//                            Project Files
// =================================================================
#include <f_line.h>
#include "ae_population.h"
#include "ae_jumping_mt.h"
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
  char* exp_setup_file_name = new char[63];
  char* out_prof_file_name  = new char[63];
  strcpy(exp_setup_file_name,  "exp_setup.ae");
  strcpy(out_prof_file_name,   "output_profile.ae");
  char* env_file_name       = NULL;
  char* pop_file_name       = NULL;
  char* sp_struct_file_name = NULL;
  bool verbose          = false;

  int32_t num_gener = 0;

  // 2) Define allowed options
  const char * options_list = "hf:n:";
  static struct option long_options_list[] = {
    { "help", no_argument,        NULL, 'h' },
    { "file", required_argument,  NULL, 'f' }, // Provide file with parameters to change
    { "num_gener",    required_argument,  NULL, 'n' }, // Provide generation number corresponding to population, environment, exp_setup and output_profile files
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
      case 'n' :
      {
        if (strcmp(optarg, "") == 0)
        {
          printf("%s: error: Option -r or --resume : missing argument.\n", argv[0]);
          exit(EXIT_FAILURE);
        }

        num_gener = atol(optarg);

        env_file_name       = new char[255];
        pop_file_name       = new char[255];
        sp_struct_file_name = new char[255];
        exp_setup_file_name = new char[255];
        out_prof_file_name  = new char[255];

        sprintf(env_file_name,       ENV_FNAME_FORMAT,       num_gener);
        sprintf(pop_file_name,       POP_FNAME_FORMAT,       num_gener);
        sprintf(sp_struct_file_name, SP_STRUCT_FNAME_FORMAT, num_gener);
        sprintf(exp_setup_file_name, EXP_S_FNAME_FORMAT,     num_gener);
        sprintf(out_prof_file_name,  OUT_P_FNAME_FORMAT,     num_gener);


        // Check existence of optional files in file system.
        // Missing files will cause the corresponding file_name variable to be nullified
        struct stat stat_buf;
        if (stat(sp_struct_file_name, &stat_buf) == -1)
        {
          if (errno == ENOENT)
          {
            delete [] sp_struct_file_name;
            sp_struct_file_name = NULL;
          }
          else
          {
            printf("%s:%d: error: unknown error.\n", __FILE__, __LINE__);
            exit(EXIT_FAILURE);
          }
        }

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
  if (param_file_name == NULL)
  {
    printf("%s: error: You must provide both a parameter and a exp_setup file. \n", argv[0]);
    exit(EXIT_FAILURE);
  }

  if (env_file_name == NULL || pop_file_name == NULL)
  {
    printf("%s: error: You must provide both an environment backup and a population backup.\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  // 5) Initialize an empty experiment manager
  printf("Load previous experiment\n");
  #ifndef __NO_X
    ae_exp_manager* exp_manager = new ae_exp_manager_X11();
  #else
    ae_exp_manager* exp_manager = new ae_exp_manager();
  #endif

  if (num_gener > 0)
  {
    exp_manager->set_first_gener(num_gener);
  }
  exp_manager->load_experiment(exp_setup_file_name, out_prof_file_name, env_file_name, pop_file_name, sp_struct_file_name, true);
  printf("Generation : %"PRId32"\n", exp_manager->num_gener());

  // 6) Interpret and apply changes
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
    if (strcmp(line->words[0], "SEED") == 0)
    {
      int32_t seed = atoi(line->words[1]) ;
      printf("\tChange of the seed to %d\n",atoi(line->words[1]));

      ae_jumping_mt* prng = new ae_jumping_mt(seed);

      // Change prng in ae_exp_manager, ae_selection and ae_spatial_structure
      printf("Change of the seed in ae_exp_nanager, ae_selection and ae_spatial_structure\t");
      ae_selection* sel = exp_manager->exp_s()->sel();
      //printf("\n%.5f ", sel->prng()->random());
      sel->set_prng(new ae_jumping_mt(*prng));
      //printf("%.5f\n", sel->prng()->random());
      if (sel->is_spatially_structured())
      {
        sel->spatial_structure()->set_prng(new ae_jumping_mt(*prng));
      }
      printf("Ok\n");

      // Change prng of each individual
      if (pop_file_name == NULL)
      {
        printf("%s: error: You must provide a backup file if any change in seed. \n", argv[0]);
        exit(EXIT_FAILURE);
      }
      printf("Loading the backup population\t");
      ae_population* pop = exp_manager->pop();
      gzFile pop_file = gzopen(pop_file_name, "r");
      if (pop_file == Z_NULL)
      {
        printf("%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, pop_file_name);
        exit(EXIT_FAILURE);
      }
      pop->load(pop_file, verbose);
      printf("Ok\n");

      #warning PRNG change disabled
    }

    delete line;
  }
  fclose(param_file);

  // 7) Save the change in the static setup files (experimental setup and output profile)
  exp_manager->write_setup_files();

  // 8) Save the change in backups
  exp_manager->save_experiment();

  // 7) Save the changed population in a new population file (similar name with _changed)
  /*printf("Save the changed population\t");
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
  printf("Ok\n");*/
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
