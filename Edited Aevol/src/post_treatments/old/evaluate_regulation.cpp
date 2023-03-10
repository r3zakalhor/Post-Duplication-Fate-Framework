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

#include <list>
#include <getopt.h>
#include <libgen.h>

#include <iostream>
#include <fstream>

#include "aevol.h"

using namespace aevol;


enum check_type {
    FULL_CHECK  = 0,
    LIGHT_CHECK = 1,
    ENV_CHECK   = 2,
    NO_CHECK    = 3
};

void print_help(char* prog_name);

int main(int argc, char* argv[]) {
  int32_t generation = 0;        // Default starting generation

  char*   output_file_name    = NULL;
  char*   lineage_file_name   = NULL;
  bool    verbose             = false;

  const char* short_options = "hVv:g:b";
  static struct option long_options[] =
      {
          {"help",        no_argument,       NULL, 'h'},
          {"version",     no_argument,       NULL, 'V'},
          {"verbose",     no_argument,       NULL, 'v'},
          {"bestonly",    no_argument,       NULL, 'b'},
          {"generation",  required_argument, NULL, 'g'},
          {0, 0, 0, 0}
      };

  bool best_only = false;
  int option;
  while ((option = getopt_long(argc, argv, short_options, long_options, NULL)) != -1) {
    switch (option) {
      case 'h' :
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      case 'V' :
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      case 'v' :
        verbose = true;
        break;
      case 'g' :
        generation = atol(optarg);
        break;
      case 'b' : best_only = true; break;
    }
  }

  // =============================
  //  Open the experience manager
  // =============================
  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(generation, true, true);

  int nb_model_env = dynamic_cast<PhenotypicTargetHandler_R*>(exp_manager->world()->phenotypic_target_handler())->
      phenotypic_target_models_.size();

  double* meta_error_with_signal;
  double* meta_error_without_signal;
  double* diff_meta_error;

  double* dual_meta_error_with_signal;
  double* dual_meta_error_without_signal;
  double* dual_diff_meta_error;

  int* env_a = new int[nb_model_env*nb_model_env];
  int* env_b = new int[nb_model_env*nb_model_env];

  double* var_env_meta_error_with_signal;
  double* var_env_meta_error_without_signal;
  double* var_env_diff_meta_error;

  int total_dual_env = 0;
  int nb_tested_env = 0;

  if (best_only) {
    meta_error_with_signal = new double[1*nb_model_env];
    meta_error_without_signal = new double[1*nb_model_env];
    diff_meta_error = new double[1*nb_model_env];

    dual_meta_error_with_signal = new double[1*nb_model_env*nb_model_env];
    dual_meta_error_without_signal = new double[1*nb_model_env*nb_model_env];
    dual_diff_meta_error = new double[1*nb_model_env*nb_model_env];

    var_env_meta_error_with_signal = new double[50];
    var_env_meta_error_without_signal = new double[50];
    var_env_diff_meta_error = new double[50];
  } else {
    meta_error_with_signal = new double[exp_manager->indivs().size()*nb_model_env
                                                    +nb_model_env*nb_model_env];
    meta_error_without_signal = new double[exp_manager->indivs().size()*(nb_model_env
                                        +nb_model_env*nb_model_env)];
    diff_meta_error = new double[exp_manager->indivs().size()*nb_model_env
                                        +nb_model_env*nb_model_env];

    dual_meta_error_with_signal = new double[exp_manager->indivs().size()*nb_model_env*nb_model_env];
    dual_meta_error_without_signal = new double[exp_manager->indivs().size()*nb_model_env*nb_model_env];
    dual_diff_meta_error = new double[exp_manager->indivs().size()*nb_model_env*nb_model_env];

    var_env_meta_error_with_signal = new double[exp_manager->indivs().size()*50];
    var_env_meta_error_without_signal = new double[exp_manager->indivs().size()*50];
    var_env_diff_meta_error = new double[exp_manager->indivs().size()*50];
  }


  for (int i = 0; i < nb_model_env; i++) {
    for (int j = 0; j < nb_model_env; j++) {
      if (i != j) {
        total_dual_env++;
      }
    }
  }

  /** Test each individual in each env model (ie constant env) **/
  for (int i = 0; i < nb_model_env; i++) {


    /** Without signal **/
    if (best_only) {
      Individual_R* indiv = dynamic_cast<Individual_R*>(exp_manager->best_indiv());

      dynamic_cast<PhenotypicTargetHandler_R*>(dynamic_cast<const Habitat_R&>(
          dynamic_cast<Individual_R*>(indiv)->grid_cell()->habitat()).phenotypic_target_handler_ptr())->set_single_env(i);

      indiv->evaluated_ = false;

      indiv->
          Evaluate(true);
      meta_error_without_signal[i] =
          indiv->
          dist_to_target_by_feature(METABOLISM);
    } else {
      int indiv_id = 0;
      for (auto indiv : exp_manager->indivs()) {
        indiv->evaluated_ = false;

        dynamic_cast<PhenotypicTargetHandler_R*>(dynamic_cast<const Habitat_R&>(
            dynamic_cast<Individual_R*>(indiv)->grid_cell()->habitat()).phenotypic_target_handler_ptr())->set_single_env(i);

        dynamic_cast<Individual_R*>(indiv)->
            EvaluateInContext(dynamic_cast<const Habitat_R&>(
                                  dynamic_cast<Individual_R*>(indiv)->grid_cell()->habitat()),true);
        meta_error_without_signal[i+indiv_id*(nb_model_env)] =
            dynamic_cast<Individual_R*>(indiv)->
                dist_to_target_by_feature(METABOLISM);
        indiv_id++;
      }
    }

    /** With signal **/
    if (best_only) {
      Individual_R* indiv = dynamic_cast<Individual_R*>(exp_manager->best_indiv());
      indiv->evaluated_ = false;

      dynamic_cast<PhenotypicTargetHandler_R*>(dynamic_cast<const Habitat_R&>(
          dynamic_cast<Individual_R*>(indiv)->grid_cell()->habitat()).phenotypic_target_handler_ptr())->set_single_env(i);

      indiv->
          EvaluateInContext(dynamic_cast<const Habitat_R&>(indiv->grid_cell()->habitat()),false);
      meta_error_with_signal[i] =
          dynamic_cast<Individual_R*>(exp_manager->best_indiv())->
              dist_to_target_by_feature(METABOLISM);
    } else {
      int indiv_id = 0;
      for (auto indiv : exp_manager->indivs()) {
        indiv->evaluated_ = false;

        dynamic_cast<PhenotypicTargetHandler_R*>(dynamic_cast<const Habitat_R&>(
            dynamic_cast<Individual_R*>(indiv)->grid_cell()->habitat()).phenotypic_target_handler_ptr())->set_single_env(i);

        dynamic_cast<Individual_R*>(indiv)->
            EvaluateInContext(dynamic_cast<const Habitat_R&>(
                                  dynamic_cast<Individual_R*>(indiv)->grid_cell()->habitat()),false);
        meta_error_with_signal[i+indiv_id*(nb_model_env)] =
            dynamic_cast<Individual_R*>(exp_manager->best_indiv())->
                dist_to_target_by_feature(METABOLISM);
        indiv_id++;
      }
    }


    /** Test each individual in each combination of envs **/
    for (int j = 0; j < nb_model_env; j++) {
      if (i!=j) {
        env_a[nb_tested_env] = i;
        env_b[nb_tested_env] = j;

        /** Without signal **/
        if (best_only) {
          Individual_R* indiv = dynamic_cast<Individual_R*>(exp_manager->best_indiv());
          indiv->evaluated_ = false;

          dynamic_cast<PhenotypicTargetHandler_R*>(dynamic_cast<const Habitat_R&>(
              dynamic_cast<Individual_R*>(indiv)->grid_cell()->habitat()).phenotypic_target_handler_ptr())->set_two_env(i,j);

          indiv->
              EvaluateInContext(dynamic_cast<const Habitat_R&>(indiv->grid_cell()->habitat()),true);
          dual_meta_error_without_signal[nb_tested_env] =
              dynamic_cast<Individual_R*>(exp_manager->best_indiv())->
                  dist_to_target_by_feature(METABOLISM);
        } else {
          int indiv_id = 0;
          for (auto indiv : exp_manager->indivs()) {
            indiv->evaluated_ = false;

            dynamic_cast<PhenotypicTargetHandler_R*>(dynamic_cast<const Habitat_R&>(
                dynamic_cast<Individual_R*>(indiv)->grid_cell()->habitat()).phenotypic_target_handler_ptr())->set_two_env(i,j);

            dynamic_cast<Individual_R*>(indiv)->
                EvaluateInContext(dynamic_cast<const Habitat_R&>(
                                      dynamic_cast<Individual_R*>(indiv)->grid_cell()->habitat()),true);;
            dual_meta_error_without_signal[nb_tested_env + indiv_id * (total_dual_env)] =
                dynamic_cast<Individual_R*>(exp_manager->best_indiv())->
                    dist_to_target_by_feature(METABOLISM);
            indiv_id++;
          }
        }

        /** With signal **/
        if (best_only) {
          Individual_R* indiv = dynamic_cast<Individual_R*>(exp_manager->best_indiv());
          indiv->evaluated_ = false;

          dynamic_cast<PhenotypicTargetHandler_R*>(dynamic_cast<const Habitat_R&>(
              dynamic_cast<Individual_R*>(indiv)->grid_cell()->habitat()).phenotypic_target_handler_ptr())->set_two_env(i,j);

          indiv->
              EvaluateInContext(dynamic_cast<const Habitat_R&>(indiv->grid_cell()->habitat()),false);
          dual_meta_error_with_signal[nb_tested_env] =
              dynamic_cast<Individual_R*>(exp_manager->best_indiv())->
                  dist_to_target_by_feature(METABOLISM);
        } else {
          int indiv_id = 0;
          for (auto indiv : exp_manager->indivs()) {
            indiv->evaluated_ = false;

            dynamic_cast<PhenotypicTargetHandler_R*>(dynamic_cast<const Habitat_R&>(
                dynamic_cast<Individual_R*>(indiv)->grid_cell()->habitat()).phenotypic_target_handler_ptr())->set_two_env(i,j);

            dynamic_cast<Individual_R*>(indiv)->
                EvaluateInContext(dynamic_cast<const Habitat_R&>(
                                      dynamic_cast<Individual_R*>(indiv)->grid_cell()->habitat()),false);;
            dual_meta_error_with_signal[nb_tested_env + indiv_id * (total_dual_env)] =
                dynamic_cast<Individual_R*>(exp_manager->best_indiv())->
                    dist_to_target_by_feature(METABOLISM);
            indiv_id++;
          }
        }

        nb_tested_env++;
      }
    }
  }

  for (int i = 0; i < 50; i++) {
    exp_manager->world()->ApplyHabitatVariation();

    if (best_only) {
      /** Without signal **/
      Individual_R* indiv = dynamic_cast<Individual_R*>(exp_manager->best_indiv());

      indiv->evaluated_ = false;


      indiv->
          Evaluate(true);
      var_env_meta_error_without_signal[i] =
          indiv->
              dist_to_target_by_feature(METABOLISM);

      /** With signal **/
      indiv = dynamic_cast<Individual_R*>(exp_manager->best_indiv());
      indiv->evaluated_ = false;
      indiv->clear_dist_sum();

      indiv->
          EvaluateInContext(dynamic_cast<const Habitat_R&>(indiv->grid_cell()->habitat()),false);
      var_env_meta_error_with_signal[i] =
          dynamic_cast<Individual_R*>(exp_manager->best_indiv())->
              dist_to_target_by_feature(METABOLISM);
    } else {
      /** Without signal **/
      int indiv_id = 0;
      for (auto indiv : exp_manager->indivs()) {
        indiv->evaluated_ = false;

        dynamic_cast<Individual_R*>(indiv)->
            EvaluateInContext(dynamic_cast<const Habitat_R&>(
                                  dynamic_cast<Individual_R*>(indiv)->grid_cell()->habitat()),
                              true);
        var_env_meta_error_without_signal[i + indiv_id * 50] =
            dynamic_cast<Individual_R*>(indiv)->
                dist_to_target_by_feature(METABOLISM);
        indiv_id++;
      }

      /** With signal **/
      indiv_id = 0;
      for (auto indiv : exp_manager->indivs()) {
        indiv->evaluated_ = false;

        dynamic_cast<Individual_R*>(indiv)->
            EvaluateInContext(dynamic_cast<const Habitat_R&>(
                                  dynamic_cast<Individual_R*>(indiv)->grid_cell()->habitat()),false);;
        var_env_meta_error_with_signal[i + indiv_id * 50] =
            dynamic_cast<Individual_R*>(exp_manager->best_indiv())->
                dist_to_target_by_feature(METABOLISM);
        indiv_id++;
      }
    }
  }


  if (best_only) {
    for (int i = 0; i < nb_model_env; i++) {
      diff_meta_error[i] = meta_error_without_signal[i]
                           - meta_error_with_signal[i];
    }

    for (int i = 0; i < total_dual_env; i++) {
      dual_diff_meta_error[i] = dual_meta_error_without_signal[i]
                                - dual_meta_error_with_signal[i];
    }

    for (int i = 0; i < 50; i++) {
      var_env_diff_meta_error[i] = var_env_meta_error_without_signal[i]
                           - var_env_meta_error_with_signal[i];
    }
  } else {
    for (int i = 0; i < exp_manager->indivs().size()*
                        (nb_model_env); i++) {
      diff_meta_error[i] = meta_error_without_signal[i]
                           - meta_error_with_signal[i];
    }

    for (int i = 0; i < exp_manager->indivs().size()*
                        (total_dual_env); i++) {
      dual_diff_meta_error[i] = dual_meta_error_without_signal[i]
                                - dual_meta_error_with_signal[i];
    }

    for (int i = 0; i < exp_manager->indivs().size()*50; i++) {
      var_env_diff_meta_error[i] = var_env_meta_error_without_signal[i] -
          var_env_meta_error_with_signal[i];
    }
  }


  std::ofstream myfile;
  myfile.open ("evaluate_regulation.csv");
  myfile<<"Generation,"<<"indiv_id,"<<"single_dual_multi,"<<"iteration,"<<"Env_Number_1,"<<"Env_Number_2,"<<"diff_meta_error"<<std::endl;

  if (best_only) {
    for (int i = 0; i < nb_model_env; i++) {
      myfile<<generation<<",0,"<<"single,0,"<<std::to_string(i)<<",-1,"<<std::to_string(diff_meta_error[i])<<std::endl;
    }
    for (int i = 0; i < total_dual_env; i++) {
      myfile<<generation<<",0,"<<"dual,0,"<<std::to_string(env_a[i])<<","<<
          std::to_string(env_b[i])<<","<<std::to_string(dual_diff_meta_error[i])<<std::endl;
    }
    for (int i = 0; i < 50; i++) {
      myfile<<generation<<",0,"<<"multi,"<<i<<","<<"-1"<<","<<
            "-1"<<","<<std::to_string(var_env_diff_meta_error[i])<<std::endl;
    }
  } else {
    for (int indiv_id=0; indiv_id < exp_manager->indivs().size(); indiv_id++) {
      for (int i = 0; i < nb_model_env; i++) {
        myfile << generation <<","<< indiv_id <<","<< "single,0" <<","<< std::to_string(i) <<","<< "-1" <<","<<
        std::to_string(diff_meta_error[i + indiv_id * (nb_model_env)]) <<
        std::endl;
      }
      for (int i = 0; i < total_dual_env; i++) {
        myfile << generation << "," << indiv_id <<","<< "dual,0" <<","<< std::to_string(env_a[i]) <<","<<
        std::to_string(env_b[i]) <<","<<
        std::to_string(dual_diff_meta_error[i + indiv_id * (total_dual_env)]) <<
        std::endl;
      }
      for (int i = 0; i < 50; i++) {
        myfile<<generation<<","<<indiv_id<<","<<"multi,"<<i<<","<<"-1"<<","<<
              "-1"<<","<<std::to_string(var_env_diff_meta_error[i+indiv_id*50])<<std::endl;
      }
    }
  }

  myfile.close();

  delete exp_manager;
  // delete mystats;
  return EXIT_SUCCESS;
}


// Print help
void print_help(char* prog_name) {
  printf("\n\
%s is a post-treatment that generates and analyses a large quantity of mutants for a lineage of ancestors.\
For each mutant we record the phenotypic effect on metabolism.\n\n\
Usage: %s [-h] -i input_file_name -o output_file_name [-b start_at_generation] [-e end_at_generation] [-p period] [-n num_mutants] [-r] [-h bin_size] [-v verbose] [-s mutation_seed]\n\
\t-h: display this screen\n\
\t-f input_file_name: lineage file to be analyzed\n\
\t-o output_file_name: name of the output file (to be written in ./stats/ancstats). In case of histogram output (-h) one file will be produced for each histogram and output_file_name will be postfixed with the generation number\n\
\t-b start_at_generation: first generation of the lineage to be analyzed (default: 0)\n\
\t-e end_at_generation: last generation of the lineage to be analyzed (default: last generation stored in the input file)\n\
\t-p period: temporal resolution of the analyze (default: 1)\n\
\t-n nb_mutants : generate and analyse nb_mutants per individual (default: 1000)\n\
Example:\n\t%s -i lineage_file -o toto.out -b 4000 -e 5000 -p 10 -n 100000 -s 19769\n", prog_name, prog_name, prog_name);

// \t-r: raw output; store the difference of metabolic error for each mutant generated (warning: the output file may quickly grow)\n
// \t-h bin_size: store the histogram with a bin_size resolution. One output file is generated for each histogram (postfixed with the generation number)\n
// \t-s mutation_seed: specify the seed to be used for the mutation random generator\n\n

}
