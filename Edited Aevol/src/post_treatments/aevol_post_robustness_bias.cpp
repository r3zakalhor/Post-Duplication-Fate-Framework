//
// Created by duazel on 27/04/2020.
//

//
// Created by duazel on 16/04/2020.
//
#include "IOJson.h"
#include "Individual.h"
#include "neutral_mutation_exp.h"
#include "Robustness_bias_output.h"

#include <iostream>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <getopt.h>
#include <string>
#include <vector>

#include <sys/stat.h>

unsigned int print_every = 0;
uint32_t seed_prng = 0;
std::string inputFile = "input.json";
std::string outputFile = "none";
unsigned int number_replications;

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
  printf("%s: Generate replications of an individual \n",
         prog_name);
  printf("\n");
  printf("Usage : %s -h or --help\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("\nOptions\n");
  printf("  -n, --nb_replications\n\tnumber of replications\n");
  printf("  -f, --file\n\tinput file\n");
  printf("  -h, --help\n\tprint this help, then exit\n");
  printf("  -V, --version\n\tprint version number, then exit\n");
  printf("  -V, --version\n\tprint version number, then exit\n");
}

void interpret_cmd_line_options(int argc, char **argv) {
  static struct option long_options[] = {{"help", no_argument, nullptr, 'h'},
                                         {"version", no_argument, nullptr, 'V'},
                                         {"print", required_argument, nullptr, 'p'},
                                         {"nb_replications", required_argument, nullptr, 'n'},
                                         {"seed", required_argument, nullptr, 's'},
                                         {"file", required_argument, nullptr, 'f'},
                                         {"output", required_argument, nullptr, 'o'},
                                         {0, 0, 0, 0}};

  int option;
  while ((option = getopt_long(argc, argv, "hVvl:n:s:f:o:", long_options,
                               nullptr)) != -1) {
    switch (option) {
    case 'h' :
      print_help(argv[0]);
      exit(EXIT_SUCCESS);
    case 'V' :
      Utils::PrintAevolVersion();
      exit(EXIT_SUCCESS);
    case 'p':
      print_every = atol(optarg);
      break;
    case 'n':
      number_replications = atol(optarg);
      break;
    case 's':
      seed_prng = atol(optarg);
      break;
    case 'f':
      inputFile.assign(optarg);
      break;
    case 'o':
      outputFile.assign(optarg);
      break;
    default:
      // An error message is printed in getopt_long, we just need to exit0
      exit(EXIT_FAILURE);
    }
  }

  if (number_replications == 0) {
    Utils::ExitWithUsrMsg("The parameter number_replications is needed");
  }
  if (seed_prng == 0) {
    Utils::ExitWithUsrMsg("The parameter seed_prng is needed");
  }
}

int main(int argc, char ** argv) {
  interpret_cmd_line_options(argc, argv);

  print_every = 10000;

  std::string folder_name = "robustness_bias";
  int folder = mkdir(folder_name.c_str(), 0755);
  std::string summary = folder_name + "/" + inputFile.substr(0, inputFile.size()-5) + "_summary_seed_" + std::to_string(seed_prng) + ".txt";
  std::ofstream summary_file(summary.c_str());

  std::cout << argv[0] << " started with :" << std::endl;
  std::cout << " - number_replications = " << number_replications << std::endl;
  std::cout << " - seed = " << seed_prng << std::endl;
  std::cout << " - file = " << inputFile << std::endl;

  IOJson inputJson(inputFile);

  auto mut_prng   = std::make_shared<JumpingMT>(seed_prng);
  auto stoch_prng = std::make_shared<JumpingMT>(seed_prng);
  Individual ancestor = Individual(inputJson.getIndividuals().front(), 0, mut_prng, stoch_prng);

  Individual* individual = new Individual(ancestor);
  individual->clear_everything_except_dna_and_promoters();
  individual->compute_phenotype();

  std::string mutation_file = folder_name + "/" + inputFile.substr(0, inputFile.size()-5) + "_mutation_seed_" + std::to_string(seed_prng) + ".csv";
  Robustness_bias_output out(*individual, "indiv.csv", mutation_file.c_str());

  summary_file << "nb_replication;nb_neutre;nb_mutant;nb_mutations;nb_neutral_mutations;sum_size_variation_neutral_mutation;multiple_mutations;mult_unneutral_mut_into_neutral_replication;nb_switch;nb_neutral_switch;nb_small_dup;nb_neutral_small_dup;nb_small_del;nb_neutral_small_del;nb_large_dup;nb_neutral_large_dup;nb_large_del;nb_neutral_large_del;nb_transloc;nb_neutral_transloc;nb_invers;nb_neutral_invers"<< std::endl;
  
  for (unsigned int id = 1; id <= number_replications; ++id) {
    std::vector<MutationEvent *> mutations;
    std::vector<bool> is_neutral;
    Individual* child = new_child(individual, mutations, is_neutral);
    out.record_replication(id, *child, *individual, mutations, is_neutral);
    delete child;
    if (id%print_every == 0){
      out.print_summary(std::cout);
      out.print_summary(summary_file);
    }// else if (id%100 == 0) {
     // std::cout << id << std::endl;
    //}
  }
  out.print_summary(summary_file);

  return 0;
}
