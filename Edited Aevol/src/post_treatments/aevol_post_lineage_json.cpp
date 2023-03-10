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

// ============================================================================
//                                   Includes
// ============================================================================
#include <list>
#include <getopt.h>
#include <assert.h>

#include "aevol.h"
#include "IndivAnalysis.h"

using namespace aevol;

// Helper functions
void interpret_cmd_line_options(int argc, char* argv[]);
void print_help(char* prog_path);
json extract_gu(IndivAnalysis* indiv);

// Command-line option variables
static char* json_file_name = nullptr;
static char* lineage_file_name = nullptr;
static int32_t begin = 0; //< First generation to analyse
static int32_t end = -1; //< Last generation to analyse
static int32_t period = 1; //< Period of analysis
static int32_t final_indiv_index = -1;
static int32_t final_indiv_rank = -1;

int main(int argc, char* argv[]) {
  lineage_file_name = new char[256];
  strcpy(lineage_file_name,"lineage.json");
  interpret_cmd_line_options(argc, argv);

  // =======================
  //  Open the lineage file
  // =======================
  gzFile lineage_file = gzopen(lineage_file_name, "r");
  if (lineage_file == Z_NULL) {
    Utils::ExitWithUsrMsg(std::string("Could not read lineage file ") +
                          lineage_file_name + "\n");
  }

  int64_t t0 = 0;
  int64_t t_end = 0;

  gzread(lineage_file, &t0, sizeof(t0));
  gzread(lineage_file, &t_end, sizeof(t_end));
  gzread(lineage_file,&final_indiv_index,sizeof(final_indiv_index));
  gzread(lineage_file,&final_indiv_rank,sizeof(final_indiv_rank));

  // =============================
  //  Open the experience manager
  // =============================
  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(t0, true, false);

  IOJson* iojson = new IOJson(exp_manager);

  FILE* json_file = nullptr;
  json_file = fopen(json_file_name,"w");

  // The current version doesn't allow for phenotypic variation nor for
  // different phenotypic targets among the grid
  if (not exp_manager->world()->phenotypic_target_shared()) {
    Utils::ExitWithUsrMsg("sorry, ancestor stats has not yet been implemented "
                          "for per grid-cell phenotypic target\n");
  }
  auto phenotypicTargetHandler =
      exp_manager->world()->phenotypic_target_handler();
  if (phenotypicTargetHandler->var_method() != NO_VAR) {
    Utils::ExitWithUsrMsg("sorry, ancestor stats has not yet been implemented "
                          "for variable phenotypic targets\n");
  }


  std::shared_ptr <JumpingMT> prng = std::make_shared<JumpingMT>(9695);

  // ==============================
  //  Prepare the initial ancestor
  // ============================
  // ==
  GridCell* grid_cell = new GridCell(lineage_file, exp_manager, nullptr);
  IndivAnalysis indiv(*(grid_cell->individual()));

  if (begin == 0) {
    iojson->addIndividual(&indiv,extract_gu(&indiv));
  }

  // ==========================================================================
  //  Replay the mutations to get the successive ancestors and analyze them
  // ==========================================================================
  ReplicationReport* rep = nullptr;

  aevol::AeTime::plusplus();
  while ((time() <= t_end) && (((time() < end) || (end == -1)))) {
#ifdef __REGUL
    rep = nullptr;
    printf("Not supported with RAEvol\n");
    exit(-11);
#else
    rep = new ReplicationReport(lineage_file, &indiv);
#endif

    // 2) Replay replication (create current individual's child)
    GeneticUnit& gen_unit = indiv.genetic_unit_nonconst(0);


    // For each genetic unit, replay the replication (undergo all mutations)
    // TODO <david.parsons@inria.fr> disabled for multiple GUs
    const auto& dnarep = rep->dna_replic_report();


    dnarep.iter_muts([&](const auto& mut) {
      gen_unit.dna()->undergo_this_mutation(*mut);

      if (period == -1) {
        iojson->addIndividual(&indiv,extract_gu(&indiv));
      }
    });

    if (period != -1) {
      if ((time() >= begin) && ((time() < end) || (end == -1)) &&
          (((time() - begin) % period) == 0)) {
        iojson->addIndividual(&indiv, extract_gu(&indiv));
      }
    }

    delete rep;

    aevol::AeTime::plusplus();
  }

  iojson->write(json_file_name);

  fclose(json_file);
  gzclose(lineage_file);
  delete [] lineage_file_name;
  delete exp_manager;
  return EXIT_SUCCESS;
}

void interpret_cmd_line_options(int argc, char* argv[]) {
  const char* short_options = "hVb:e:P:j:";
  static struct option long_options[] = {
      {"help",        no_argument,       nullptr, 'h'},
      {"version",     no_argument,       nullptr, 'V'},
      {"begin",       required_argument, nullptr, 'b'},
      {"end",         required_argument, nullptr, 'e'},
      {"period",      required_argument, nullptr, 'P'},
      {"json",        required_argument, nullptr, 'j'},
      {0, 0, 0, 0}
  };

  int option;
  while ((option = getopt_long(argc, argv, short_options, long_options,
                               nullptr)) != -1) {
    switch (option) {
    case 'h' :
      print_help(argv[0]);
      exit(EXIT_SUCCESS);
    case 'V' :
      Utils::PrintAevolVersion();
      exit(EXIT_SUCCESS);
    case 'b' :
      begin = atol(optarg);
      break;
    case 'e' :
      end = atol(optarg);
      break;
    case 'P' :
      period = atol(optarg);
      break;
    case 'j' :
      delete [] json_file_name;
      json_file_name = new char[strlen(optarg) + 1];
      sprintf(json_file_name, "%s", optarg);
      break;
    default:
      // An error message is printed in getopt_long, we just need to exit
      exit(EXIT_FAILURE);
    }
  }

  // There should be only one remaining arg: the lineage file
  if (optind != argc - 1) {
    Utils::ExitWithUsrMsg("please specify a lineage file");
  }

  assert(period > 0 ||period == -1);

  lineage_file_name = new char[strlen(argv[optind]) + 1];
  sprintf(lineage_file_name, "%s", argv[optind]);
}

void print_help(char* prog_path) {
  // Get the program file-name in prog_name (strip prog_path of the path)
  char* prog_name; // No new, it will point to somewhere inside prog_path
  if ((prog_name = strrchr(prog_path, '/'))) {
    prog_name++;
  }
  else {
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
  printf("%s: generate a json file for the provided lineage.\n",
         prog_name);
  printf("\n");
  printf("Usage : %s -h or --help\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("   or : %s LINEAGE_FILE [-b TIMESTEP] [-e TIMESTEP] [-P PERIOD] [-j json]\n",
         prog_name);
  printf("\nOptions\n");
  printf("  -h, --help\n\tprint this help, then exit\n");
  printf("  -V, --version\n\tprint version number, then exit\n");
  printf("  -b, --begin TIMESTEP\n");
  printf("\ttimestep at which to start the extraction\n");
  printf("  -e, --end TIMESTEP\n");
  printf("\ttimestep at which to stop the extraction\n");
  printf("  -P, --period\n");
  printf("\tperiod with which to perform the extraction\n");
  printf("  -j, --output\n");
  printf("\tjson file name\n");
}

json extract_gu(IndivAnalysis* indiv){
  json gu_list = json::array();
  for (auto& gen_unit: indiv->genetic_unit_list_nonconst()) {
    std::string dna = gen_unit.dna()->data();
    int32_t length = gen_unit.dna()->length();
    dna.resize(length);
    json a_gu;
    a_gu["seq"] = dna;
    gu_list.emplace_back(a_gu);
  }
  return gu_list;
}