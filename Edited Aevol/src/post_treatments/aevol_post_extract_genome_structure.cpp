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
//                              Includes
// =================================================================
#include <cerrno>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cinttypes>
#include <cmath>
#include <cassert>
#include <getopt.h>
#include <sys/stat.h>
#include<iostream>
#include<fstream>
using namespace std;

#include <list>

#include <zlib.h>

// =================================================================
//                            Project Files
// =================================================================
#include "aevol.h"


using std::list;
using namespace aevol;

// =================================================================
//                     Command line option variables
// =================================================================
int32_t wanted_rank = -1;
int32_t wanted_index = -1;
int64_t num_gener = 0;
int32_t mutation_type = 0;
int32_t nb_mutants = -1;

// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);
void interpret_cmd_line_options(int argc, char* argv[]);

int main(int argc, char* argv[]) {
#ifdef BASE_2
  interpret_cmd_line_options(argc, argv);



  // ------------------------------------------------------
  //  Load the backup and get the individual to be mutated
  // ------------------------------------------------------

  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(num_gener, true, false);

  if (exp_manager->output_m()->record_tree() == false) {
    // The following instruction is needed to ensure that methods
    // like ae_dna::do_deletion, ae_dna::do_inversion, etc
    // will create ae_mutation objects (otherwise they return NULL)
    exp_manager->output_m()->init_tree(exp_manager, 100);
  }


  // TODO: factor with duplicated code in robustness.cpp
  Individual* initial_indiv = nullptr;

  { // (local scope for variable `indivs` used as a shorthand)
    bool found = false;
    int32_t best_index = -1;
    double best_fitness = 0;
    int32_t current_index = -1;
    list<Individual*> indivs = exp_manager->indivs();
    for (auto indiv = indivs.rbegin();
         not found and indiv != indivs.rend(); ++indiv) {
      current_index = (*indiv)->id();
      if (wanted_index != -1 and current_index == wanted_index) {
	printf("  Searching for individual of index %d (%d)\r",wanted_index,current_index);
        found = true;
        initial_indiv = (*indiv);
      }
      else
	{
	  
	printf("  No index specified, searching for the best individual (%d)\r",current_index++);
	 (*indiv)->Evaluate();
      if ((*indiv)->dist_to_target_by_feature(METABOLISM) > best_fitness)
	{
	  initial_indiv = (*indiv);
	  best_fitness = (*indiv)->dist_to_target_by_feature(METABOLISM);
	  found = true;
	}
  }
    }

    if (not found) {
      Utils::ExitWithUsrMsg("sorry, the individual you have requested has not "
                            "been found");
    }
  }

  printf("\n Processing individual %llu",initial_indiv->id());
  
  initial_indiv->Evaluate();
  //  initial_indiv->compute_statistical_data();
  //  initial_indiv->compute_non_coding();



  // ---------------------
  //  Prepare the output
  // ---------------------


  char mutation_type_name[24];
 

  Individual* mutant = NULL;
  int32_t nb_genetic_units = initial_indiv->nb_genetic_units();
  double* relative_lengths_genetic_units = NULL;
  int32_t u = 0;
  double alea, cumul;
  int32_t pos1;
  int32_t mut_length;
  int16_t align_score1, align_score2;
  bool invert;
  VisAVis* alignment_1 = NULL;
  VisAVis* alignment_2 = NULL;
  int32_t nb_pairs;
  Dna* initial_dna = NULL;
  int32_t initial_len;
  bool rear_done;
  int32_t nb_genes_at_breakpoints;
  int32_t nb_genes_in_segment;
  int32_t nb_genes_in_replaced_segment;
  double metabolic_error_after_inversion = -1.0;
  double metabolic_error_after_switch = -1.0;




  int32_t non_coding_started = 0;
  int32_t non_coding_length = 0;
  int32_t total_non_coding_length = 0;
  int32_t ncid = 0;
  int32_t init_pos = 0;
  int32_t end_pos = 0;

    // *********************  Exhaustive insertion of terminators at each position of the genome  *************************
    pos1 = -1;
    mut_length = -1;
    align_score1 = align_score2 = -1;
    invert = false;


    
    //   int32_t result_tab[30000];
    //   int32_t neutral_pos[30000];

    for (const auto& gu: initial_indiv->genetic_unit_list()) {
      initial_len = gu.dna()->length();

      printf (" (genome length: %ld)\n",initial_len);

      // if (initial_len > 30000)
      // 	{
      // 	  printf("Genomes of more than 30000 bp not handled ; abort\n");
      // 	  return(0);
      // 	}

      //     for (pos1 = 0; pos1 < initial_len; pos1++)
      //    {
      //	result_tab[pos1] = 0;
      //	neutral_pos[pos1] = 0;
      //  }


       // initialize neutral_pos table
       int32_t nb_neutral_pos = 0;
       ofstream file;

       int start_neutral = 0;
       int in_first_zone = 0;
       unsigned long int end_of_ORI_zone = -1;
       unsigned long int last_transition;
       int previous = -1;
       int nb_zone = 0;
       int nc_length = 0;

       file.open("genome_structure.txt");
       file << "position;non_neutral_to_switch;non_neutral_to_terminator_insertion\n";
       
       for (pos1 = 0;pos1 < initial_len;pos1++)
	{

	  // testing point mutations
	  mutant = new Individual(*initial_indiv);
	  mutant->genetic_unit(u).dna()->do_switch(pos1);
	  // Evaluate the mutant, compute its statistics
	  mutant->ReevaluateInContext(initial_indiv->habitat());
	  metabolic_error_after_switch = mutant->dist_to_target_by_feature(METABOLISM);
	  delete mutant;	  
	  
	  if (metabolic_error_after_switch == initial_indiv->dist_to_target_by_feature(METABOLISM))
	    {
	      file << pos1 << ",0,";
	    }
	  else
	    {

	      file << pos1 << ",1,";
	    }

	  
	  // testing insertion of a terminator
	  mutant = new Individual(*initial_indiv);
	  mutant->genetic_unit(u).dna()->do_small_insertion(pos1,30,"000000000000000111111111111111");
	  // Evaluate the mutant, compute its statistics
	  mutant->ReevaluateInContext(initial_indiv->habitat());
	  metabolic_error_after_switch = mutant->dist_to_target_by_feature(METABOLISM);
	  delete mutant;

	  if (pos1 == initial_len -1)
	    {
	      if (start_neutral == 1)
	      {
		printf("%d %d %d %d\n",++nb_zone,end_of_ORI_zone+pos1-last_transition,last_transition,end_of_ORI_zone);
		nc_length +=end_of_ORI_zone+pos1-last_transition;
	      }
	    }
	  
	  if (metabolic_error_after_switch == initial_indiv->dist_to_target_by_feature(METABOLISM))
	    {

	      file <<"0\n";
	      if (pos1 == 0)
		{
		  start_neutral = 1;
		  in_first_zone = 1;
		}
	      if (previous == 1)
		{ // start of a neutral zone
			      previous = 0;
			      last_transition = pos1;
		}

	      
	    }
	  else
	    {
	      file <<"1\n";
	      if (previous == 0)
		{ //end of a neutral zone
		  if (in_first_zone)
		    {
		      end_of_ORI_zone = pos1-1;
		      in_first_zone = 0;
		    }
		  else
		    {
		      nc_length += pos1 - (last_transition + 1);
		      if (pos1 - (last_transition + 1) > 2)
		      printf("%d %d %d %d\n",++nb_zone,pos1 - (last_transition + 1),last_transition,pos1);
		    }
		}
	      previous = 1;
	    }
	}

       printf("\n\ntotal nb of neutral zones : %d --- total neutral size : %d\n",nb_zone,nc_length);
      
    file.close();
    }
    u++;

    


    delete exp_manager;
    
    return EXIT_SUCCESS;
#elif BASE_4
    printf("Post-Treatment not available in 4 bases version of Aevol\n");
#endif
}

void interpret_cmd_line_options(int argc, char* argv[]) {
  const char* options_list = "hVg:r:i:m:n:";
  static struct option long_options_list[] = {
      {"help",          no_argument,       NULL, 'h'},
      {"version",       no_argument,       NULL, 'V'},
      {"gener",         required_argument, NULL, 'g'},
      {"rank",          required_argument, NULL, 'r'},
      {"index",         required_argument, NULL, 'i'},
      {"mutation-type", required_argument, NULL, 'm'},
      {"nb-mutants",    required_argument, NULL, 'n'},
      {0, 0, 0,                                  0}
  };

  int option = -1;
  bool rank_already_set = false;
  bool index_already_set = false;
  while ((option = getopt_long(argc, argv, options_list, long_options_list,
                               NULL)) != -1) {
    switch (option) {
      case 'h':
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      case 'V':
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      case 'g':
        if (strcmp(optarg, "") == 0) {
          fprintf(stderr,
                  "%s: error: Option -g or --gener : missing argument.\n",
                  argv[0]);
          exit(EXIT_FAILURE);
        }
        num_gener = atol(optarg);
        break;
      case 'i':
        wanted_index = atol(optarg);
        index_already_set = true;
        break;
    }
  }


}

/*!
  \brief

*/
void print_help(char* prog_path) {
  // Get the program file-name in prog_name (strip prog_path of the path)
  char* prog_name; // No new, it will point to somewhere inside prog_path
  if ((prog_name = strrchr(prog_path, '/'))) {
    prog_name++;
  }
  else { prog_name = prog_path; }

  printf("\n");
  printf(
      "*********************** aevol - Artificial Evolution ******************* \n");
  printf(
      "*                                                                      * \n");
  printf(
      "*                  Genome Structure post-treatment program             * \n");
  printf(
      "*                                                                      * \n");
  printf(
      "************************************************************************ \n");
  printf("\n\n");
  printf("This program is Free Software. No Warranty.\n");
  printf("\n");
  printf("Usage : %s -h\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf(
      "   or : %s -g NUMGENER [-i INDEX]\n",
      prog_name);
  printf("\n");
  printf(
      "This program insert terminators at each position of an initial individual, \n");
  printf(
      "and output the neutral vs non-neutral result. For sake of comparison, it also \n");
  printf(
      "performs base switch at the same positions.\n");
  printf(
	 "The processed individual is extracted from an aevol backup file at generation -g.\n");
  printf("the index -i of the individual can be provided, otherwise the best individual of the population is processed \n\n");
  printf(
      "detailed output is recorded in file genome_structure.txt.\n");
  printf("output summary is displayed in the console (id of the neutral zone, length of the neutral zone, starting position, ending position)\n");
  printf("\n");
  printf("\t-h or --help    : Display this help, then exit\n");
  printf("\n");
  printf("\t-V or --version : Print version number, then exit\n");
  printf("\n");
  printf("\t-g NUMGENER or --gener NUMGENER : \n");
  printf(
      "\t                  Generation of the backup containing the individual of interest\n");
  printf("\n");
  printf("\t-i INDEX or --index INDEX : \n");
  printf(
      "\t                  Index of individual of interest. Should be comprised between 0 and N-1, where\n");
  printf("\t                  N is the size of the population. if no index is specified, the program uses\n");
  printf("\t                  the best individual in the population\n");
  printf("\n");
  printf("\n");

  printf("\n");
}
