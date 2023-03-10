//
// Created by duazel on 16/04/2020.
//
#include "IOJson.h"
#include "Individual.h"
#include "neutral_mutation_exp.h"
#include "neutral_mutation_output.h"

#include <iostream>
#include <cstdio>
#include <cstring>
#include <getopt.h>
#include <string>
#include <vector>

#include <sys/stat.h>

bool verbose = false;
uint32_t seed_prng = 0;
std::string inputFile = "input.json";
std::string outputFile = "none";
unsigned int wanted_size = 0;
int delta = 0;

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
  printf("%s: Accumulate neutral mutations on an individual to change its size without\n"
         "changing its phenotype.\n", prog_name);
  printf("\n");
  printf("Usage : %s -h or --help\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("\nOptions\n");
  printf("  -l, --length\n\twanted sequence length\n");
  printf("  -d, --delta\n\twanted modification of the size\n");
  printf("  -s, --seed\n\trandom generator seed\n");
  printf("  -h, --help\n\tprint this help, then exit\n");
  printf("  -V, --version\n\tprint version number, then exit\n");
  printf("  -f, --input_file_name\n\tdefault is input.json\n");
  printf("  -o, --output_file_name\n");
}

void interpret_cmd_line_options(int argc, char **argv) {
  static struct option long_options[] = {{"help", no_argument, nullptr, 'h'},
                                         {"version", no_argument, nullptr, 'V'},
                                         {"verbose", no_argument, nullptr, 'v'},
                                         {"length", required_argument, nullptr, 'l'},
                                         {"delta", required_argument, nullptr, 'd'},
                                         {"seed", required_argument, nullptr, 's'},
                                         {"file", required_argument, nullptr, 'f'},
                                         {"output", required_argument, nullptr, 'o'},
                                         {0, 0, 0, 0}};

  int option;
  while ((option = getopt_long(argc, argv, "hVvl:d:s:f:o:", long_options,
                               nullptr)) != -1) {
    switch (option) {
    case 'h' :
      print_help(argv[0]);
      exit(EXIT_SUCCESS);
    case 'V' :
      Utils::PrintAevolVersion();
      exit(EXIT_SUCCESS);
    case 'v':
      verbose = true;
      break;
    case 'l':
      wanted_size = atol(optarg);
      break;
    case 'd':
      delta = atol(optarg);
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
      // An error message is printed in getopt_long, we just need to exit
      exit(EXIT_FAILURE);
    }
  }

  if (wanted_size == 0 && delta == 0) {
    Utils::ExitWithUsrMsg("One of the parameters wanted_size or delta is needed");
  } else if (wanted_size != 0 && delta != 0) {
    Utils::ExitWithUsrMsg("Only one of the parameters wanted_size or delta has to be defined.\n"
                          "Otherwise its ambiguous");
  }
  if (seed_prng == 0) {
    Utils::ExitWithUsrMsg("The parameter seed_prng is needed");
  }

}

int main(int argc, char ** argv) {
  interpret_cmd_line_options(argc, argv);

  std::cout << argv[0] << " started with :" << std::endl;
  if (wanted_size != 0) {
    std::cout << " - wanted_sequence_length = " << wanted_size << std::endl;
  } else {
    std::cout << " - delta = " << delta << std::endl;
  }
  std::cout << " - seed = " << seed_prng << std::endl;
  std::cout << " - input file = " << inputFile << std::endl;
  std::cout << " - output file = " << outputFile << std::endl;
  IOJson inputJson(inputFile);

  std::string folder_name = "directional_neutral_mutation_accumulation";
  int folder = mkdir(folder_name.c_str(), 0755);
  std::string summary = folder_name + "/" + inputFile.substr(0, inputFile.size()-5) + "_summary_seed_" + std::to_string(seed_prng) + ".csv";
  std::string mutation = folder_name + "/" + inputFile.substr(0, inputFile.size()-5) + "_mutation_seed_" + std::to_string(seed_prng) + ".csv";
  out::init(summary.c_str(), mutation.c_str());

  if (wanted_size == 0) {
    wanted_size = inputJson.getIndividuals()[0]->amount_of_dna() + delta;
  }

  auto mut_prng   = std::make_shared<JumpingMT>(seed_prng);
  auto stoch_prng = std::make_shared<JumpingMT>(seed_prng);
  Individual ancestor = Individual(inputJson.getIndividuals()[0], 0, mut_prng, stoch_prng);

  Individual* indiv = run_to_size(wanted_size, &ancestor);
  
  std::string genetic_seq = indiv->genetic_unit_sequence(0);
  inputJson.setIndividualSequence(0,0, genetic_seq.c_str());

  if (outputFile == "none"){
    if (delta != 0){
      if (delta > 0){
        outputFile = inputFile.substr(0, inputFile.size()-5) + "_+" + std::to_string(delta) + "bp_seed_" + std::to_string(seed_prng) + ".json";
      }
      else{
        outputFile = inputFile.substr(0, inputFile.size()-5) + "_" + std::to_string(delta) + "bp_seed_" + std::to_string(seed_prng) + ".json";
      }      
    }
    else{
      outputFile = inputFile.substr(0, inputFile.size()-5) + "_" + std::to_string(genetic_seq.length()) + "bp_seed_" + std::to_string(seed_prng) + ".json";
    }
    
  }
  outputFile = folder_name + "/" + outputFile;
  std::cout << outputFile;
  inputJson.write(outputFile);

  return 0;
}
