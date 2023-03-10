// ****************************************************************************
//
//          Aevol - An in silico experimental evolution platform
//
// ****************************************************************************
//
// Copyright: See the AUTHORS file provided with the package or <www.aevol.fr>
// Web: http://www.aevol.fr/
// E-mail: See <http://www.aevol.fr/contact/>
// Original Authors : Guillaume Beslon, Carole Knibbe, David Parsons, Jonathan Rouzaud-Cornabas
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
#include <signal.h>

#include <cstdint>
#include <fstream>
#include <limits>

// =================================================================
//                            Project Files
// =================================================================
#include "aevol.h"

using namespace aevol;

// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);
void extract_network(Individual_R* indiv);

int main(int argc, char* argv[])
{
    // Initialize command-line option variables with default values
    bool best_only = true;
    int32_t num_gener = -1;

    // 2) Define allowed options
    const char * options_list = ":hVv:b:g:";
    static struct option long_options_list[] = {
            {"help",      no_argument,       NULL, 'h'},
            {"version",   no_argument,       NULL, 'V' },
            {"best",   no_argument,       NULL, 'b' },
            {"gener",     required_argument, NULL, 'g' },
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
                //print_help(argv[0]);
                exit(EXIT_SUCCESS);
            }
            case 'V' :
            {
                Utils::PrintAevolVersion();
                exit(EXIT_SUCCESS);
            }
            case 'b' :
            {
                best_only = true;
                break;
            }
            case 'g' :
            {
                if (strcmp(optarg, "") == 0)
                {
                    printf("%s: error: Option -g or --gener : missing argument.\n", argv[0]);
                    exit(EXIT_FAILURE);
                }

                num_gener = atol(optarg);

                break;
            }
        }
    }

    // If num_gener is not provided, assume last gener
    if (num_gener == -1) {
        num_gener = OutputManager::last_gener();
    }

    // Open the files
    auto exp_manager = new ExpManager();
    exp_manager->load(num_gener, false, false);


    // The best individual is already known because it is the last in the list
    // Thus we do not need to know anything about the environment and to evaluate the individuals

    // Parse the individuals
    if (best_only)
    {
        Individual_R* best = dynamic_cast<Individual_R*>(exp_manager->best_indiv());
        best->do_transcription_translation_folding(); // We need to recompute proteins if not already done (ie if using a population file and not a full backup)
        extract_network(best);
    }
    else
    {
        for (const auto& indiv: exp_manager->indivs()) {
            indiv->do_transcription_translation_folding();
            extract_network(dynamic_cast<Individual_R*>(indiv));
        }
    }

    delete exp_manager;

    return EXIT_SUCCESS;
}

void extract_network(Individual_R* indiv) {
    std::ofstream network;
    network.open("network.csv",std::ofstream::trunc);
    network<<"Individual,"<<"Enhancer_or_Inhibitor,"<<"Value"<<std::endl;

    for (auto& rna: indiv->get_rna_list_coding()) {
        for (unsigned int i = 0; i < ((Rna_R*)rna)->nb_influences(); i++) {
            //compute the activity
            if (((Rna_R*)rna)->_enhancing_coef_list[i] > 0)
            {
                network<<indiv->id()<<",1,"<<((Rna_R*)rna)->_enhancing_coef_list[i]<<std::endl;
            }

            if (((Rna_R*)rna)->_operating_coef_list[i] > 0)
            {
                network<<indiv->id()<<",0,"<<((Rna_R*)rna)->_enhancing_coef_list[i]<<std::endl;
            }
        }
    }

    network.flush();
    network.close();
}