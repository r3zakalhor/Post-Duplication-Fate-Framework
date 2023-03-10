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
#include <string>
// =================================================================
//                            Project Files
// =================================================================
#include "aevol.h"

using namespace aevol;

// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);
void extract_network(Individual_R* indiv, double* fabs_metaerror_loss, double* fabs_fitness_loss,
                     double* fabs_metaerror_loss_percent, double* fabs_fitness_loss_percent);
void dump_network(Individual_R* indiv, double* fabs_metaerror_loss, double* fabs_fitness_loss,
                  double* fabs_metaerror_loss_percent, double* fabs_fitness_loss_percent);
void filter_network(Individual_R* indiv, double* fabs_metaerror_loss, double* fabs_fitness_loss,
                    double* fabs_metaerror_loss_percent, double* fabs_fitness_loss_percent);
void extract_network_single_target_model(Individual_R* best, int nb_phenotypic_target_models,
                                         double** ptm_fabs_metaerror_loss, double** ptm_fabs_fitness_loss,
                                         double** ptm_fabs_metaerror_loss_percent,
                                         double** ptm_fabs_fitness_loss_percent);

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
        best->do_transcription_translation_folding();

        int nb_edges = 0;
        for (auto &rna: best->get_rna_list_coding()) {
            nb_edges+=((Rna_R *) rna)->nb_influences();
        }

        double* fabs_metaerror_loss = new double[nb_edges];
        double* fabs_fitness_loss = new double[nb_edges];
        double* fabs_metaerror_loss_percent = new double[nb_edges];
        double* fabs_fitness_loss_percent = new double[nb_edges];

        for (int i = 0; i < nb_edges; i++) {
            fabs_metaerror_loss[i] = 0;
            fabs_fitness_loss[i] = 0;
            fabs_metaerror_loss_percent[i] = 0;
            fabs_fitness_loss_percent[i] = 0;
        }


        int nb_iteration = 100;
        printf("Running %d evals for %d edges\n",nb_iteration,nb_edges);
        for (int i = 0; i < nb_iteration; i++) {
            printf("Iteration %d\n",i);
            exp_manager->world()->ApplyHabitatVariation();

            best->evaluated_ = false;
            best->Evaluate();

            double base_metaerror = best->dist_to_target_by_feature(METABOLISM);
            double base_fitness = best->fitness();

            int i_edges = 0;

            for (auto &rna: best->get_rna_list_coding()) {
                for (unsigned int i = 0; i < ((Rna_R *) rna)->nb_influences(); i++) {
                    double enhance_backup = ((Rna_R *) rna)->_enhancing_coef_list[i];
                    double operate_backup = ((Rna_R *) rna)->_operating_coef_list[i];
                    ((Rna_R *) rna)->_enhancing_coef_list[i] = 0;
                    ((Rna_R *) rna)->_operating_coef_list[i] = 0;

                    best->evaluated_ = false;
                    best->Evaluate();

                    fabs_metaerror_loss[i_edges] += std::fabs(base_metaerror-best->dist_to_target_by_feature(METABOLISM));
                    fabs_fitness_loss[i_edges] += std::fabs(base_fitness-best->fitness());

                    fabs_metaerror_loss_percent[i_edges] += (std::fabs(base_metaerror-best->dist_to_target_by_feature(METABOLISM)))/best->dist_to_target_by_feature(METABOLISM);
                    fabs_fitness_loss_percent[i_edges] += (std::fabs(base_fitness-best->fitness()))/best->fitness();

                    ((Rna_R *) rna)->_enhancing_coef_list[i] = enhance_backup;
                    ((Rna_R *) rna)->_operating_coef_list[i] = operate_backup;

                    i_edges++;
                }
            }
        }

        for (int i = 0; i < nb_edges; i++) {
            fabs_metaerror_loss[i] /= nb_iteration;
            fabs_fitness_loss[i] /= nb_iteration;
            fabs_metaerror_loss_percent[i] /= nb_iteration;
            fabs_fitness_loss_percent[i] /= nb_iteration;
        }

        extract_network(best,fabs_metaerror_loss,fabs_fitness_loss,
                        fabs_metaerror_loss_percent,fabs_fitness_loss_percent);
        filter_network(best,fabs_metaerror_loss,fabs_fitness_loss,fabs_metaerror_loss_percent,fabs_fitness_loss_percent);
        dump_network(best,fabs_metaerror_loss,fabs_fitness_loss,fabs_metaerror_loss_percent,fabs_fitness_loss_percent);

        int nb_phenotypic_target_models = dynamic_cast<PhenotypicTargetHandler_R*>(exp_manager->world()->
                phenotypic_target_handler())->phenotypic_target_models_.size();
        printf("Running with a single phenotypic target models : %d\n",nb_phenotypic_target_models);

        double** ptm_fabs_metaerror_loss = new double*[nb_phenotypic_target_models];
        double** ptm_fabs_fitness_loss = new double*[nb_phenotypic_target_models];
        double** ptm_fabs_metaerror_loss_percent = new double*[nb_phenotypic_target_models];
        double** ptm_fabs_fitness_loss_percent = new double*[nb_phenotypic_target_models];

        for (int i = 0; i < nb_phenotypic_target_models; i++) {
            ptm_fabs_metaerror_loss[i] = new double[nb_edges];
            ptm_fabs_fitness_loss[i] = new double[nb_edges];
            ptm_fabs_metaerror_loss_percent[i] = new double[nb_edges];
            ptm_fabs_fitness_loss_percent[i] = new double[nb_edges];
        }


        for (int target_id = 0; target_id < nb_phenotypic_target_models; target_id++) {

            for (int j = 0; j < nb_edges; j++) {
                ptm_fabs_metaerror_loss[target_id][j] = 0;
                ptm_fabs_fitness_loss[target_id][j] = 0;
                ptm_fabs_metaerror_loss_percent[target_id][j] = 0;
                ptm_fabs_fitness_loss_percent[target_id][j] = 0;
            }

            printf("Testing with phenotypic target model %d\n",target_id);
            dynamic_cast<PhenotypicTargetHandler_R*>(exp_manager->world()->phenotypic_target_handler())->set_single_env(target_id);

            best->evaluated_ = false;
            best->Evaluate();

            double base_metaerror = best->dist_to_target_by_feature(METABOLISM);
            double base_fitness = best->fitness();

            int i_edges = 0;

            for (auto &rna: best->get_rna_list_coding()) {
                for (unsigned int i = 0; i < ((Rna_R *) rna)->nb_influences(); i++) {
                    double enhance_backup = ((Rna_R *) rna)->_enhancing_coef_list[i];
                    double operate_backup = ((Rna_R *) rna)->_operating_coef_list[i];
                    ((Rna_R *) rna)->_enhancing_coef_list[i] = 0;
                    ((Rna_R *) rna)->_operating_coef_list[i] = 0;

                    best->evaluated_ = false;
                    best->Evaluate();

                    //printf("Testing with phenotypic target model %d : %lf %lf\n",target_id,base_metaerror,best->dist_to_target_by_feature(METABOLISM));

                    ptm_fabs_metaerror_loss[target_id][i_edges] += std::fabs(base_metaerror-best->dist_to_target_by_feature(METABOLISM));
                    ptm_fabs_fitness_loss[target_id][i_edges] += std::fabs(base_fitness-best->fitness());

                    ptm_fabs_metaerror_loss_percent[target_id][i_edges] += (std::fabs(base_metaerror-best->dist_to_target_by_feature(METABOLISM)))/best->dist_to_target_by_feature(METABOLISM);
                    ptm_fabs_fitness_loss_percent[target_id][i_edges] += (std::fabs(base_fitness-best->fitness()))/best->fitness();

                    ((Rna_R *) rna)->_enhancing_coef_list[i] = enhance_backup;
                    ((Rna_R *) rna)->_operating_coef_list[i] = operate_backup;

                    i_edges++;
                }
            }
        }

        extract_network_single_target_model(best,nb_phenotypic_target_models,ptm_fabs_metaerror_loss,ptm_fabs_fitness_loss,ptm_fabs_metaerror_loss_percent,ptm_fabs_fitness_loss_percent);

    }
    else
    {
/*        for (const auto& indiv: exp_manager->indivs()) {
            indiv->do_transcription_translation_folding();
            //indiv->Evalutate();
            extract_network(dynamic_cast<Individual_R*>(indiv));
        }*/
    }

    delete exp_manager;

    return EXIT_SUCCESS;
}

void extract_network(Individual_R* indiv, double* fabs_metaerror_loss, double* fabs_fitness_loss, double* fabs_metaerror_loss_percent, double* fabs_fitness_loss_percent) {
    std::ofstream network;
    network.open("network_knockout.csv",std::ofstream::trunc);
    network<<"Individual,"<<"Enhancer_or_Inhibitor,"<<"Value,"<<"Metaerror_lost,"<<"Fitness_lost,Metaerror_lost_percent,Fitness_lost_percent" <<std::endl;

    int nb_edges_enhance = 0, nb_edges_operating = 0, nb_edges_both = 0, nb_edges = 0;
    int i_edges = 0;

    for (auto& rna: indiv->get_rna_list_coding()) {
        for (unsigned int i = 0; i < ((Rna_R*)rna)->nb_influences(); i++) {
            //std::cout<<"Influence "<<i<<" value is "<<((Rna_R*)rna)->_enhancing_coef_list[i]<<" "<<((Rna_R*)rna)->_operating_coef_list[i]<<std::endl;
            //compute the activity

            if (((Rna_R*)rna)->_enhancing_coef_list[i] > 0)
            {
                network<<indiv->id()<<",1,"<<((Rna_R*)rna)->_enhancing_coef_list[i]<<","<<fabs_metaerror_loss[i_edges]<<","
                       <<fabs_fitness_loss[i_edges]<<","<<fabs_metaerror_loss_percent[i_edges]<<","<<fabs_fitness_loss_percent[i_edges]<<std::endl;
            }

            if (((Rna_R*)rna)->_operating_coef_list[i] > 0)
            {
                network<<indiv->id()<<",0,"<<((Rna_R*)rna)->_operating_coef_list[i]<<","<<fabs_metaerror_loss[i_edges]<<","
                        <<fabs_fitness_loss[i_edges]<<","<<fabs_metaerror_loss_percent[i_edges]<<","<<fabs_fitness_loss_percent[i_edges]<<std::endl;
            }

            if ((((Rna_R *) rna)->_enhancing_coef_list[i] > 0) && (((Rna_R *) rna)->_operating_coef_list[i] > 0)) {
                nb_edges_both++;
            } else {
                //std::cout<<"Influence "<<i<<" value is "<<((Rna_R*)rna)->_enhancing_coef_list[i]<<" "<<((Rna_R*)rna)->_operating_coef_list[i]<<std::endl;
                //compute the activity
                if (((Rna_R *) rna)->_enhancing_coef_list[i] > 0) {
                    nb_edges_enhance++;
                }

                if (((Rna_R *) rna)->_operating_coef_list[i] > 0) {
                    nb_edges_operating++;
                }
            }

            i_edges++;
            nb_edges++;
        }
    }

    network.flush();
    network.close();

    float filter_value = 0.0;
    std::string str_filter_value = std::to_string(filter_value);

    std::string file_name = "network_edges_" + str_filter_value + ".csv";

    network.open(file_name, std::ofstream::trunc);
    network << "Individual," << "nb_enhancing," << "nb_inhibitor," << "nb_both,nb_edges"<< std::endl;
    network << indiv->id() << "," << nb_edges_enhance << "," << nb_edges_operating << "," << nb_edges_both << ","
            << nb_edges << std::endl;
    network.close();
}

void filter_network(Individual_R* indiv, double* fabs_metaerror_loss, double* fabs_fitness_loss, double* fabs_metaerror_loss_percent, double* fabs_fitness_loss_percent) {
    float filter_values[3] = {0.00001, 0.0001, 0.001};

    for (float filter_value : filter_values) {

        std::string str_filter_value = std::to_string(filter_value);
        std::string file_name = "network_filtered_" + str_filter_value + ".csv";
        std::ofstream network;
        network.open(file_name, std::ofstream::trunc);
        network << "Individual," << "Enhancer," << "Inhibitor," << "Both," << "Value,"
                << "Metaerror_lost,Fitness_lost,Metaerror_lost_percent,Fitness_lost_percent" << std::endl;

        int i_edges = 0;

        int nb_edges_enhance = 0, nb_edges_operating = 0, nb_edges_both = 0, nb_edges = 0;
        int filter_nb_edges_enhance = 0, filter_nb_edges_operating = 0, filter_nb_edges_both = 0, filter_nb_edges = 0;

        for (auto &rna: indiv->get_rna_list_coding()) {
            for (unsigned int i = 0; i < ((Rna_R *) rna)->nb_influences(); i++) {
                int both = 0;
                if (fabs_fitness_loss_percent[i_edges] >= filter_value) {
                    if ((((Rna_R *) rna)->_enhancing_coef_list[i] > 0) &&
                        (((Rna_R *) rna)->_operating_coef_list[i] > 0)) {
                        network << indiv->id() << ",1,1,1," << ((Rna_R *) rna)->_enhancing_coef_list[i] << ","
                                << fabs_metaerror_loss[i_edges] << "," << fabs_fitness_loss[i_edges]<<
                        ","<<fabs_metaerror_loss_percent[i_edges]<<","<<fabs_fitness_loss_percent[i_edges]<< std::endl;
                        filter_nb_edges_both++;
                    } else {
                        //std::cout<<"Influence "<<i<<" value is "<<((Rna_R*)rna)->_enhancing_coef_list[i]<<" "<<((Rna_R*)rna)->_operating_coef_list[i]<<std::endl;
                        //compute the activity
                        if (((Rna_R *) rna)->_enhancing_coef_list[i] > 0) {
                            network << indiv->id() << ",1,0,0," << ((Rna_R *) rna)->_enhancing_coef_list[i] << ","
                                    << fabs_metaerror_loss[i_edges] << "," << fabs_fitness_loss[i_edges]<<
                            ","<<fabs_metaerror_loss_percent[i_edges]<<","<<fabs_fitness_loss_percent[i_edges]<< std::endl;
                            filter_nb_edges_enhance++;
                        }

                        if (((Rna_R *) rna)->_operating_coef_list[i] > 0) {
                            network << indiv->id() << ",0,1,0," << ((Rna_R *) rna)->_operating_coef_list[i] << ","
                                    << fabs_metaerror_loss[i_edges] << "," << fabs_fitness_loss[i_edges]<<
                            ","<<fabs_metaerror_loss_percent[i_edges]<<","<<fabs_fitness_loss_percent[i_edges]<< std::endl;
                            filter_nb_edges_operating++;
                        }
                    }
                    filter_nb_edges++;
                }

                if ((((Rna_R *) rna)->_enhancing_coef_list[i] > 0) && (((Rna_R *) rna)->_operating_coef_list[i] > 0)) {
                    nb_edges_both++;
                } else {
                    //std::cout<<"Influence "<<i<<" value is "<<((Rna_R*)rna)->_enhancing_coef_list[i]<<" "<<((Rna_R*)rna)->_operating_coef_list[i]<<std::endl;
                    //compute the activity
                    if (((Rna_R *) rna)->_enhancing_coef_list[i] > 0) {
                        nb_edges_enhance++;
                    }

                    if (((Rna_R *) rna)->_operating_coef_list[i] > 0) {
                        nb_edges_operating++;
                    }
                }
                nb_edges++;

                i_edges++;
            }
        }

        network.flush();
        network.close();

        file_name = "network_edges_" + str_filter_value + ".csv";

        network.open(file_name, std::ofstream::trunc);
        network << "Individual," << "nb_enhancing," << "nb_inhibitor," << "nb_both,nb_edges," << "filter_nb_enhancing,"
                << "filter_nb_inhibitor," << "filter_nb_both,filter_nb_edges" << std::endl;
        network << indiv->id() << "," << nb_edges_enhance << "," << nb_edges_operating << "," << nb_edges_both << ","
                << nb_edges << "," <<
                filter_nb_edges_enhance << "," << filter_nb_edges_operating << "," << filter_nb_edges_both << ","
                << filter_nb_edges << std::endl;
        network.close();

    }

}



void dump_network(Individual_R* indiv, double* fabs_metaerror_loss, double* fabs_fitness_loss, double* fabs_metaerror_loss_percent, double* fabs_fitness_loss_percent) {

    float filter_values[4] = {0.0, 0.00001, 0.0001, 0.001};

    for (float filter_value : filter_values) {

        std::string str_filter_value = std::to_string(filter_value);
        std::string file_name = "network_dump_"+str_filter_value+".csv";
        std::ofstream network;
        network.open(file_name, std::ofstream::trunc);
        network << "Individual," << "Source,"<<"Destination,"<<"Enhancer_or_Inhibitor," <<
                "Value" << "Metaerror_lost,Fitness_lost,Metaerror_lost_percent,Fitness_lost_percent" << std::endl;

        int i_edges = 0;

        int nb_edges_enhance = 0, nb_edges_operating = 0, nb_edges_both = 0, nb_edges = 0;
        int filter_nb_edges_enhance = 0, filter_nb_edges_operating = 0, filter_nb_edges_both = 0, filter_nb_edges = 0;

        for (auto &rna: indiv->get_rna_list_coding()) {
            for (unsigned int i = 0; i < ((Rna_R *) rna)->nb_influences(); i++) {
                for (auto& protein : rna->transcribed_proteins()) {
                    if (fabs_fitness_loss_percent[i_edges] >= filter_value) {
                        if (((Rna_R *) rna)->_enhancing_coef_list[i] > 0) {
                            network << indiv->id() << "," << protein->shine_dal_pos() << ","
                                   << dynamic_cast<Rna_R*>(rna)->_protein_list[i]->shine_dal_pos()<<","
                                   << "1," << ((Rna_R *) rna)->_enhancing_coef_list[i] << ","
                                   << fabs_metaerror_loss[i_edges] << "," << fabs_fitness_loss[i_edges]<<","
                                    << fabs_metaerror_loss_percent[i_edges] << "," << fabs_fitness_loss_percent[i_edges] << std::endl;
                        }

                        if (((Rna_R *) rna)->_operating_coef_list[i] > 0) {
                            network << indiv->id() << "," << protein->shine_dal_pos() << ","
                                    << dynamic_cast<Rna_R*>(rna)->_protein_list[i]->shine_dal_pos()<<","
                                    << "0," << ((Rna_R *) rna)->_operating_coef_list[i] << ","
                                    << fabs_metaerror_loss[i_edges] << "," << fabs_fitness_loss[i_edges] <<","
                                    << fabs_metaerror_loss_percent[i_edges] << "," << fabs_fitness_loss_percent[i_edges]<< std::endl;
                        }
                    }
                }

                i_edges++;
            }
        }

        network.flush();
        network.close();

    }

}

void extract_network_single_target_model(Individual_R* indiv, int nb_phenotypic_target_models,
                                         double** ptm_fabs_metaerror_loss, double** ptm_fabs_fitness_loss,
                                         double** ptm_fabs_metaerror_loss_percent,
                                         double** ptm_fabs_fitness_loss_percent) {
    std::ofstream network;
    network.open("network_knockout_single_env.csv",std::ofstream::trunc);
    network<<"Individual,"<<"Enhancer_or_Inhibitor,"<<"TargetModel,"<<"Value"<<"Metaerror_lost,Fitness_lost,Metaerror_lost_percent,Fitness_lost_percent"<<std::endl;



    for (int target_id = 0; target_id < nb_phenotypic_target_models; target_id++) {
        int i_edges = 0;
        for (auto &rna: indiv->get_rna_list_coding()) {
            for (unsigned int i = 0; i < ((Rna_R *) rna)->nb_influences(); i++) {
                //std::cout<<"Influence "<<i<<" value is "<<((Rna_R*)rna)->_enhancing_coef_list[i]<<" "<<((Rna_R*)rna)->_operating_coef_list[i]<<std::endl;
                //compute the activity
                if (((Rna_R *) rna)->_enhancing_coef_list[i] > 0) {
                    network << indiv->id() << ",1," <<target_id<<","<< ((Rna_R *) rna)->_enhancing_coef_list[i] << ","
                            << ptm_fabs_metaerror_loss[target_id][i_edges]<<","
                            << ptm_fabs_fitness_loss[target_id][i_edges] <<","
                            << ptm_fabs_metaerror_loss_percent[target_id][i_edges]<<","
                            << ptm_fabs_fitness_loss_percent[target_id][i_edges] << std::endl;
                }

                if (((Rna_R *) rna)->_operating_coef_list[i] > 0) {
                    network << indiv->id() << ",0," <<target_id<<","<< ((Rna_R *) rna)->_operating_coef_list[i] << ","
                            << ptm_fabs_metaerror_loss[target_id][i_edges]<<","
                            << ptm_fabs_fitness_loss[target_id][i_edges] << ","
                            << ptm_fabs_metaerror_loss_percent[target_id][i_edges] <<","
                            << ptm_fabs_fitness_loss_percent[target_id][i_edges] << std::endl;
                }
                i_edges++;
            }
        }
    }

    network.flush();
    network.close();

}