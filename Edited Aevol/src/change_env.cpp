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
#include <environment.h>

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
  char* env_file_name   = NULL;
  bool verbose          = false;

  // 2) Define allowed options
  const char * options_list = "he:f:";
  static struct option long_options_list[] = {
    { "help", no_argument,        NULL, 'h' },
    { "env",  required_argument,  NULL, 'e' }, // Provide environment file,
    { "file", required_argument,  NULL, 'f' }, // Provide file with environment parameters to change
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
      case 'e' :
      {
        if (strcmp(optarg, "") == 0)
        {
          printf("%s: error: Option -e or --env : missing argument.\n", argv[0]);
          exit(EXIT_FAILURE);
        }

        env_file_name = optarg;

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
  if (param_file_name == NULL || env_file_name == NULL)
  {
    printf("%s: error: You must provide both a parameter file and a environment backup.\n", argv[0]);
    exit(EXIT_FAILURE);
  }


  // 5) Initialize an empty environment then load the backup environment
  printf("Loading the backup environment\t");
  Environment* env = new Environment();
  gzFile env_file = gzopen(env_file_name, "r");
  if (env_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, env_file_name);
    exit(EXIT_FAILURE);
  }
  env->load(env_file);
  printf("Ok\n");

  // 6) Interpret and apply changes
  printf("Interpret and apply changes\n");
  FILE* param_file = fopen(param_file_name,  "r");
  if (param_file == NULL)
  {
    printf("%s:%d: error: could not open parameter file %s\n", __FILE__, __LINE__, param_file_name);
    exit(EXIT_FAILURE);
  }

  f_line* line;
  while ((line = line(param_file)) != NULL)
  {
    if (strcmp(line->words[0], "ENV_ADD_GAUSSIAN") == 0)
    {
      env->add_gaussian(atof(line->words[1]), atof(line->words[2]), atof(line->words[3]));
      env->add_initial_gaussian(atof(line->words[1]), atof(line->words[2]), atof(line->words[3])); //usefull in case of autoregressive mean variation to compute delta_m
      printf("\tAddition of a gaussian with %f, %f, %f \n",atof(line->words[1]), atof(line->words[2]), atof(line->words[3]));
    }
    else if (strcmp(line->words[0], "ENV_VARIATION") == 0)
    {
      static bool env_var_already_set = false;
      if (env_var_already_set)
      {
        printf("%s:%d: ERROR in param file : duplicate entry for %s.\n", __FILE__, __LINE__, line->words[0]);
        exit(EXIT_FAILURE);
      }
      env_var_already_set = true;

      if (strcmp(line->words[1], "none") == 0)
      {
        assert(line->nb_words == 2);
        env->set_var_method(NO_VAR);
        printf("\tNo more environmental variation\n");
      }
      else if (strcmp(line->words[1], "autoregressive_mean_variation") == 0)
      {
        assert(line->nb_words == 5);
        env->set_var_method(AUTOREGRESSIVE_MEAN_VAR);
        env->set_var_sigma(atof(line->words[2]));
        env->set_var_tau(atol(line->words[3]));
        env->set_var_prng(new ae_jumping_mt(atoi(line->words[4])));
        printf("\tChange of environmental variation to a autoregressive mean variation with sigma=%f, tau=%ld and seed=%d\n", atof(line->words[2]),atol(line->words[3]),atoi(line->words[4]));
      }
      else if (strcmp(line->words[1], "add_local_gaussians") == 0)
      {
        assert(line->nb_words == 3);
        env->set_var_method(LOCAL_GAUSSIANS_VAR);
        env->set_var_prng(new ae_jumping_mt(atoi(line->words[2])));
        printf("\tChange of environmental variation to a local gaussians variation with seed=%d\n", atoi(line->words[2]));
      }
      else
      {
        printf("%s:%d: ERROR in param file : unknown environment variation method.\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
      }
    }

    delete line;
  }
  fclose(param_file);
  printf("Ok\n");

  // 7) Save the changed environment in a new environment file (similar name with _changed)
  printf("Save the changed environment\t");
  char* new_env_file_name   = NULL;
  new_env_file_name = new char[strlen(env_file_name)-3+strlen("_changed.ae")+1];
  sprintf(new_env_file_name, "%s_changed.ae",strtok(env_file_name, "."));
  gzFile new_env_file = gzopen(new_env_file_name, "w");
  if (new_env_file == Z_NULL)
  {
    printf("%s:%d: error: could not open backup file %s\n", __FILE__, __LINE__, env_file_name);
    exit(EXIT_FAILURE);
  }
  env->save(new_env_file);
  gzclose(new_env_file);
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
	printf("Usage : change_env -h\n");
	printf("   or : change_env [-f param_file -p env_file]\n");
	printf("  -h, --help  Display this screen\n");
	printf("  -e, --env  Specify environment to change\n");
	printf("  -f, --file  Specify file with parameters to change\n");
	printf("Change a environment file as specified in the parameter file.\n");
    printf("(default: param_to_change.in)\n\n");
}
