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
#include <getopt.h>
#include <cstdlib>
#include <cstdio>
#include "IOJson.h"
#include <string>
#include <assert.h>

// =================================================================
//                            Project Files
// =================================================================
#include "aevol.h"

using namespace aevol;

// Command-line option variables
static char* triangles_file_name  = nullptr;
static char* sequence_file_name  = nullptr;
static char* json_file_name = nullptr;
static bool all_indiv = false;
static bool by_index = false;
static bool x_axis = false;
static bool y_axis = false;
static int16_t gu = -1;
static int32_t timestep = -1;
static int32_t ind = -1;
static int16_t x_pos = -1;
static int16_t y_pos = -1;

// Helper functions
void print_help(char* prog_path);
void interpret_cmd_line_options(int argc, char* argv[]);

json analyse_indiv(Individual* indiv, FILE* triangles_file, FILE* sequence_file,
                   FILE* json_file, int16_t gu, const PhenotypicTarget& phenotypicTarget);
void analyse_gu(GeneticUnit* gen_unit, int32_t gen_unit_number, FILE* triangles_file,
                const PhenotypicTarget& phenotypicTarget);




int main(int argc, char* argv[]) {
  interpret_cmd_line_options(argc, argv);

  assert(x_axis == y_axis);
  //make sure that no more than one option is selected
  int nb_option = 0;
  if(all_indiv){nb_option++;}
  if(by_index){nb_option++;}
  if(x_axis){nb_option++;}
  assert(nb_option <= 1);

  // Open the files
  FILE* triangles_file = nullptr;
  FILE* sequence_file = nullptr;
  FILE* json_file = nullptr;

  if (triangles_file_name != nullptr) {
    triangles_file = fopen(triangles_file_name, "w");

    // Write file headers
    int key = 1;
    fprintf(triangles_file, "# %2.d individual's identifier (id)\n", key++);
    fprintf(triangles_file, "# %2.d chromosome or plasmid (c_or_p)\n", key++);
    fprintf(triangles_file, "# %2.d strand\n", key++);
    fprintf(triangles_file, "# %2.d protein position (pos)\n", key++);
    fprintf(triangles_file, "# %2.d length (len)\n", key++);
    fprintf(triangles_file, "# %2.d position of last translated nucleotide (lpos)\n", key++);
    fprintf(triangles_file, "# %2.d primary sequence (sequence)\n", key++);
    fprintf(triangles_file, "# %2.d mean (m)\n", key++);
    fprintf(triangles_file, "# %2.d width (w)\n", key++);
    fprintf(triangles_file, "# %2.d height (h)\n", key++);
    fprintf(triangles_file, "# %2.d concentration (c)\n", key++);
    fprintf(triangles_file, "# %2.d feature (f)\n", key++);
    fprintf(triangles_file, "# %2.d promoter position (prom_pos)\n", key++);
    fprintf(triangles_file, "# %2.d RNA length (rna_len)\n", key++);
    fprintf(triangles_file, "# %2.d basal level (basal_level)\n", key++);
    fprintf(triangles_file, "\n");
    fprintf(triangles_file,
            "id c_or_p strand pos len lpos sequence m w h c f "
                "prom_pos rna_len basal_level\n");
  }
  if (sequence_file_name != nullptr) {
    sequence_file = fopen(sequence_file_name,"w");
  }
  if (json_file_name != nullptr) {
    json_file = fopen(json_file_name,"w");
  }

  auto exp_manager = new ExpManager();
  exp_manager->load(timestep, false, false);

  IOJson* io_json = new IOJson(exp_manager);
  json gu_list = json::array();

  if(by_index){
      assert(ind >= 0 && ind < exp_manager->grid_width() * exp_manager->grid_height());
  }
  if(x_axis){
      //we already know that x_axis = y_axis
      assert(x_pos >= 0 && x_pos < exp_manager->grid_width());
      assert(y_pos >= 0 && y_pos < exp_manager->grid_height());
  }

  // The best individual is already known because it is the last in the list
  // Thus we do not need to know anything about the environment and to evaluate
  // the individuals

  // Parse the individuals
  if(by_index)
  {
      int32_t orig_x = ind / exp_manager->world()->height();
      int32_t orig_y = ind % exp_manager->world()->height();
      Individual* indiv = exp_manager->indiv_by_position(orig_x,orig_y);
      
      gu_list = analyse_indiv(indiv, triangles_file, sequence_file, json_file, gu, 
      #ifdef __REGUL
      dynamic_cast<PhenotypicTargetHandler_R*>(&indiv->habitat().
        phenotypic_target_handler_nonconst())->phenotypic_target_model(0)
      #else
      indiv->habitat().phenotypic_target()
      #endif
      );
      io_json->addIndividual(indiv, gu_list);
  }
  else if(x_axis){
      Individual* indiv = exp_manager->indiv_by_position(x_pos, y_pos);

      gu_list = analyse_indiv(indiv, triangles_file, sequence_file, json_file, gu, 
      #ifdef __REGUL
      dynamic_cast<PhenotypicTargetHandler_R*>(&indiv->habitat().
        phenotypic_target_handler_nonconst())->phenotypic_target_model(0)      
      #else
      indiv->habitat().phenotypic_target()
      #endif
      );
      io_json->addIndividual(indiv, gu_list);
  }
  else if (all_indiv) {
    for (const auto& indiv: exp_manager->indivs()) {
      indiv->do_transcription_translation_folding(); // We need to recompute proteins if not already done (ie if using a population file and not a full backup)

      gu_list = analyse_indiv(indiv, triangles_file, sequence_file, json_file, gu, 
      #ifdef __REGUL
      dynamic_cast<PhenotypicTargetHandler_R*>(&indiv->habitat().
        phenotypic_target_handler_nonconst())->phenotypic_target_model(0)
      #else
      indiv->habitat().phenotypic_target()
      #endif
      );
      io_json->addIndividual(indiv, gu_list);
    }
  }
  else
  {
    Individual* best = exp_manager->best_indiv();
    best->do_transcription_translation_folding(); // We need to recompute proteins if not already done (ie if using a population file and not a full backup)
    gu_list = analyse_indiv(best, triangles_file, sequence_file, json_file, gu, best->habitat().phenotypic_target()); // list of GU of the individual
    io_json->addIndividual(best, gu_list);
  }

  if (sequence_file_name != nullptr) {
    fclose(sequence_file);
  }
  if (triangles_file_name != nullptr) {
    fclose(triangles_file);
  }
  if (json_file_name != nullptr) {
    io_json->write(json_file_name);
    fclose(json_file);
  }

  delete [] triangles_file_name;
  delete [] sequence_file_name;
  delete [] json_file_name;

  delete exp_manager;

  return EXIT_SUCCESS;
}

// Parsing an individual
inline json analyse_indiv(Individual* indiv, FILE* triangles_file,
                          FILE* sequence_file, FILE* json_file, int16_t gu,
                          const PhenotypicTarget & phenotypicTarget) {

  json gu_list = json::array();

  if (gu == -1) { // We want to treat all genetic units
    int32_t gen_unit_number = 0;
    for (auto& gen_unit: indiv->genetic_unit_list_nonconst()) {

      // Get the sequence of the GU
      std::string dna = gen_unit.dna()->data();
      int32_t length = gen_unit.dna()->length();
      dna.resize(length);
      json a_gu;
      a_gu["seq"] = dna;
      gu_list.emplace_back(a_gu);

      if(triangles_file != nullptr) {
        analyse_gu(&gen_unit, gen_unit_number, triangles_file,
                   phenotypicTarget);
      }
      if (sequence_file != nullptr) {
        const char* dna = gen_unit.dna()->data();
        // The sequences of different GUs are separated by a space
        if (gen_unit_number > 0) fprintf(sequence_file, " ");
 #ifdef BASE_2
        fprintf(sequence_file, "%.*s", length, dna);
#else
	int i;
	for (i=0;i<length;i++)
	  {
	    if (dna[i] == '0') fprintf(sequence_file, "A");
	    if (dna[i] == '1') fprintf(sequence_file, "T");
	    if (dna[i] == '2') fprintf(sequence_file, "C");
	    if (dna[i] == '3') fprintf(sequence_file, "G");
	  }
	fprintf(sequence_file, "\n");
#endif
      }

      gen_unit_number++;
    }
  }
  else { // User has specified a genetic unit
    GeneticUnit* gen_unit = &indiv->genetic_unit_nonconst(gu);

    // Get the sequence of the GU
    std::string dna = gen_unit->dna()->data();
    int32_t length = gen_unit->dna()->length();
    dna.resize(length);
    json a_gu;
    a_gu["seq"] = dna;
    gu_list.emplace_back(a_gu);

    if(triangles_file != nullptr) {
        analyse_gu(gen_unit, gu, triangles_file, phenotypicTarget);
    }
    if (sequence_file != nullptr) {
      const char* dna = gen_unit->dna()->data();
      fprintf(sequence_file, "%.*s", length, dna);
    }
  }
  // We go to next line in each file
  if (triangles_file != nullptr) {
    fprintf(triangles_file, "\n");
  }
  if (sequence_file != nullptr) {
    fprintf(sequence_file, "\n");
  }
  if (json_file != nullptr) {
  }

  return gu_list;
}

// Parsing a GU
inline void analyse_gu(GeneticUnit* gen_unit, int32_t gen_unit_number,
                       FILE* triangles_file,
                       const PhenotypicTarget& phenotypicTarget) {
  // Construct the list of all rnas
  auto llrnas = gen_unit->rna_list();
  auto lrnas = llrnas[LEADING];
  lrnas.splice(lrnas.end(), llrnas[LAGGING]);
  // Parse this list
  int rna_nb = 0;
  int i = 0;

  for (const auto& rna: lrnas) {

    for (const auto& protein: rna.transcribed_proteins()) {
      double mean = protein->mean();
      int nfeat = -1;
      for (size_t i = 0 ;
           i <= static_cast<size_t>(phenotypicTarget.nb_segments()) - 1 ;
           ++i) {
        if ((mean > phenotypicTarget.segments()[i]->start) and
            (mean < phenotypicTarget.segments()[i]->stop)) {
          nfeat = phenotypicTarget.segments()[i]->feature;
          break;
        }
      }

      char *dummy;
      fprintf(triangles_file,
        "%llu %s %s %" PRId32 " %" PRId32 " %" PRId32
        " %s %f %f %f %f %d %" PRId32 " %" PRId32 " %f\n",
        gen_unit->indiv()->id(),
        gen_unit_number != 0 ? "PLASMID" :
        "CHROM  ",
        protein->strand() == LEADING ? "LEADING" :
        "LAGGING",
        protein->first_translated_pos(),
        protein->length(),
        protein->last_translated_pos(),
        dummy = protein->AA_sequence('_'),
        mean,
        protein->width(),
        protein->height(),
        protein->concentration(),
        nfeat,
        rna.promoter_pos(),
        rna.transcript_length(),
        rna.basal_level());
    }
  }
}


void print_help(char* prog_path) {
    // Get the program file-name in prog_name (strip prog_path of the path)
    char *prog_name; // No new, it will point to somewhere inside prog_path
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
    printf("%s:\n", prog_name);
    printf("\tExtracts the genotype and/or data about the phenotype of individuals\n");
    printf("\tin the provided population and write them into text files easy to parse\n");
    printf("\twith e.g. matlab.\n");
    printf("\n");
    printf("Usage : %s -h\n", prog_name);
    printf("   or : %s -V or --version\n", prog_name);
    printf("   or : %s [-t TIMESTEP] [-S SEQ_FILE] [-T TRIANGLE_FILE] [-J JSON_FILE] [-U NUM_GU] [-a] [-x X -y Y]\n",
           prog_name);
    printf("\nOptions\n");
    printf("  -h, --help\n\tprint this help, then exit\n");
    printf("  -V, --version\n\tprint version number, then exit\n");
    printf("  -t TIMESTEP\n");
    printf("\tspecify timestep of the individual(s) of interest\n");
    printf("  -S SEQ_FILE\n");
    printf("\textract sequences into file SEQ_FILE\n");
    printf("  -T TRIANGLE_FILE\n");
    printf("\textract phenotypic data into file TRIANGLE_FILE\n");
    printf("  -J JSON_FILE\n");
    printf("\textract phenotypic data into file JSON_FILE\n");
    printf("  -U NUM_GU\n");
    printf("\tonly treat genetic unit #NUM_GU (default: treat all genetic units)\n");
    printf("  -a\n");
    printf("\ttreat all the individuals (default: treat only the best)\n");
    printf("  -i IND\n");
    printf("\tonly treat individual #IND (default: treat only the best)\n");
    printf("  -x X -y Y\n");
    printf("\tonly treat individual at position X, Y on the grid\n");


    printf("\n\
This program extracts some data about the individuals and write\n\
them into text files easy to parse with e.g. matlab.\n\
\n\
Two kinds of data can be extracted :\n\
\n\
 * data about the phenotype (option -t) : write information about\n\
   the proteins in a text file. A space delimits two proteins, a\n\
   new line delimits two individuals. For each protein, the output\n\
   is \"m_h_w_c_r_s_f_l_z_g\" where :\n\
       * m, h, w and c are the mean, height, width and concentration of the protein\n\
       * r is an identifier of the rna it belongs (useful to\n\
           know if several proteins are on the same rna)\n\
       * s indicates the strand (LEADING/LAGGING)\n\
       * f and l are the first and last translated base\n\
       * z indicates the feature (at the center of the protein)\n\
       * g indicates the genetic unit to which the protein belongs (0=chromosome, 1=plasmid)\n\
\n\
 * sequences of the individuals (option -s) : write the sequences\n\
   in a text file. A new line delimits two individuals. In case\n\
   there are several GU, they are separated by whitespaces.\n\
\n\
With option -b, only the best individual is treated.\n\
\n\
The input can be either a generation number, in which case we\n\
will attempt to load a full backup tree, or a population file,\n\
in which case features of the proteins won't be outputed as we\n\
need to know the environment to infer them.\n\
\n\
Examples :\n\
\n\
For generation 20000, write infos about the phenotypes of all the\n\
individuals in phe_020000 and the sequences of all the\n\
individuals in seq_020000 :\n\
\n\
   extract -r 20000 -t phe_020000 -s seq_020000\n\
\n\
For generation 20000, write the best individual's sequence in\n\
seq_020000_best :\n\
\n\
   extract -b -r 20000 -s seq_020000_best\n\
or extract -b -p populations/pop_020000.ae -s seq_020000_best\n");
}

void interpret_cmd_line_options(int argc, char* argv[]) {
  // Define allowed options
  const char * options_list = "hVt:aU:S:T:J:i:x:y:";
  static struct option long_options_list[] = {
      {"help",      no_argument,        nullptr, 'h'},
      {"version",   no_argument,        nullptr, 'V'},
      {"timestep",  required_argument,  nullptr, 't'},
      {"all",       no_argument,        nullptr, 'a'},
      {"gu",        required_argument,  nullptr, 'U'},
      {"sequence",  required_argument,  nullptr, 'S'},
      {"triangles", required_argument,  nullptr, 'T'},
      {"json",      required_argument,  nullptr, 'J'},
      {"index",     required_argument,  nullptr, 'i'},
      {"xaxis",     required_argument,  nullptr, 'x'},
      {"yaxis",     required_argument,  nullptr, 'y'},
      {0, 0, 0, 0}
  };

  // Get actual values of the command-line options
  int option;
  while ((option = getopt_long(argc, argv, options_list,
                               long_options_list, nullptr)) != -1) {
    switch (option) {
      case 'h' : {
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      }
      case 'V' : {
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      }
      case 't' : {
        timestep = atol(optarg);
        break;
      }
      case 'a' : {
        all_indiv = true;
        break;
      }
      case 'U' : {
        gu = atoi(optarg);
        break;
      }
      case 'S' : {
        sequence_file_name = new char[strlen(optarg) + 1];
        sprintf(sequence_file_name, "%s", optarg);
        break;
      }
      case 'T' : {
        triangles_file_name = new char[strlen(optarg) + 1];
        sprintf(triangles_file_name, "%s", optarg);
        break;
      }
      case 'J' : {
        json_file_name = new char[strlen(optarg) + 1];
        sprintf(json_file_name, "%s", optarg);
        break;
      }
      case 'i' : {
        ind = atoi(optarg);
        by_index = true;
        break;
      }
      case 'x' : {
        x_pos = atoi(optarg);
        x_axis = true;
        break;
      }
      case 'y' : {
        y_pos = atol(optarg);
        y_axis = true;
        break;
      }
    }
  }

  // If timestep wasn't provided, use default
  if (timestep < 0) {
    timestep = OutputManager::last_gener();
  }

  // If neither the sequence_file_name nor the triangles_file_name was provided,
  // we will output only the sequence in a default-named file
  if (sequence_file_name == nullptr && triangles_file_name == nullptr && json_file_name == nullptr) {
    sequence_file_name = new char[255];
    strcpy(sequence_file_name, "sequence");
  }
}
