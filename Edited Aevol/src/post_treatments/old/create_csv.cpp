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
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/stat.h>  // for the permission symbols used with mkdir
#include <errno.h>
#include <HybridFuzzy.h>


// =================================================================
//                            Project Files
// =================================================================
#include "aevol.h"

using namespace aevol;

// ========
//  TO DO
// ========
//
//  * option --color ?
//  * Raevol-specific output (EPS file with the network) ?



// =================================================================
//                         Function declarations
// =================================================================
void print_help(char* prog_path);


// The height of each triangle is proportional to the product c*m, where c is the
// concentration of the protein and m its intrinsic efficiency (depending on its
// aminoacid sequence). In the case of Raevol, the concentration used here is the
// final one, i.e. the one reached after all the time steps of the lifetime.
// If a coding sequence has several promoters, only one triangle is drawn.
void draw_triangles(Individual* indiv, const PhenotypicTarget& target, char * directoryName);



// In the case of Raevol, the profile is drawn using the final concentrations
// of the proteins, i.e. the ones reached after all the time steps of the lifetime.
void draw_pos_neg_profiles(Individual * indiv, const PhenotypicTarget& target, char * directoryName);



// In the case of Raevol, the phenotype is drawn using the final concentrations
// of the proteins, i.e. the ones reached after all the time steps of the lifetime.
void draw_phenotype(Individual * indiv, const PhenotypicTarget& target, char * directoryName, int generation = -1);




// The chromosome is drawn as a circle. Coding sequences are represented by arcs starting
// at the start codon and ending at the stop codon. Coding sequences that are on the
// leading (resp. lagging) strand are drawn outside (resp. inside) the circle. If a coding
// sequence has several promoters, only one arc is drawn.
void draw_genetic_unit_with_CDS(GeneticUnit* gen_unit, char * directoryName);


// The chromosome is drawn as a circle. Transcribed sequences are represented by arcs starting
// at the first trasncribed position codon and ending at the last transcribed position. mRNAs
// that are on the leading (resp. lagging) strand are drawn outside (resp. inside) the circle.
// mRNAs which include at least one coding sequence are black, the others are gray.
 void draw_genetic_unit_with_mRNAs(GeneticUnit* gen_unit, char * directoryName);


int main(int argc, char* argv[])
{
  // =================================================================
  //                      Get command-line options
  // =================================================================
  //
  // 1) Initialize command-line option variables with default values

  //~ bool  verbose           = false;
  //~ char* backup_file_name  = NULL;
  int64_t num_gener       = -1;
  int32_t indiv_index     = -1;
  int32_t indiv_rank      = -1;

  // 2) Define allowed options
  const char * options_list = "hVvi:r:g:";
  static struct option long_options_list[] = {
    {"help",      no_argument,       NULL, 'h'},
    {"version",   no_argument,       NULL, 'V' },
    {"index",     required_argument, NULL, 'i'},
    {"rank",      required_argument, NULL, 'r'},
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
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      }
      case 'V' :
      {
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      }
      case 'i' :
        indiv_index  = atol(optarg);
        break;
      case 'r' :
        indiv_rank  = atol(optarg);
        break;
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

  // Check mandatory arguments
  if (num_gener == -1)
  {
    printf("%s: error: You must provide a generation number.\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  if (indiv_index != -1 && indiv_rank != -1)
  {
    printf("%s: error: You must provide either the index or the rank of the individual to plot.\n", argv[0]);
    exit(EXIT_FAILURE);
  }



  printf("Creating eps files for generation %" PRId64 "...\n", num_gener);

  // =================================================================
  //                       Read the backup file
  // =================================================================
  Individual*  indiv;

  // Load the simulation
  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(num_gener, true, false);

  if (indiv_index == -1 && indiv_rank == -1)
  {
    indiv = exp_manager->best_indiv();
  }
  else
  {
    // TODO <david.parsons@inria.fr> tmp disabled
//     if (indiv_rank != -1) {
//       indiv = new Individual(*exp_manager->indiv_by_rank(indiv_rank), false);
//     }
//     else {
//       indiv = new Individual(*exp_manager->indiv_by_id(indiv_index), false);
//     }
  }

  // The constructor of the exp_manager has read the genomes of the individuals
  // and located their promoters, but has not performed the translation nor the
  // phenotype computation. We must do it now.
  // However, as the individuals in the backups are sorted, we don't need to
  // evaluate all the individuals, only the one we are interested in
#ifdef __REGUL
  Individual_R* rindiv = dynamic_cast<Individual_R*>(indiv);
#else
  indiv->Evaluate();
#endif



  // =================================================================
  //                      Create the EPS files
  // =================================================================
  char directory_name[64];
  snprintf(directory_name, 63, "analysis-generation%06" PRId64, num_gener);

  // Check whether the directory already exists and is writable
  if (access(directory_name, F_OK) == 0)
  {
//       struct stat status;
//       stat(directory_name, &status);
//       if (status.st_mode & S_IFDIR) cout << "The directory exists." << endl;
//       else cout << "This path is a file." << endl;

    if (access(directory_name, X_OK | W_OK) != 0)
    {
      fprintf(stderr, "Error: cannot enter or write in directory %s.\n", directory_name);
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    // Create the directory with permissions : rwx r-x r-x
    if (mkdir(directory_name, S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) != 0)
    {
      fprintf(stderr, "Error: cannot create directory %s.\n", directory_name);
      exit(EXIT_FAILURE);
    }
  }



  // =================================================================
  //                  Write the data in the EPS files
  // =================================================================

  GeneticUnit*  indiv_main_genome = &indiv->genetic_unit_nonconst(0);

  printf("Creating the EPS file with the triangles of the chosen individual... ");
  fflush(stdout);
  draw_triangles(indiv, indiv->phenotypic_target(), directory_name);
  printf("OK\n");

  printf("Creating the EPS file with the positive and negatives profiles of the chosen individual... ");
  fflush(stdout);
  draw_pos_neg_profiles(indiv, indiv->phenotypic_target(), directory_name);
  printf("OK\n");

#ifndef __REGUL
  printf("Creating the EPS file with the phenotype of the chosen individual... ");
  fflush(stdout);
  draw_phenotype(indiv, indiv->phenotypic_target(), directory_name);
  printf("OK\n");
#endif

  printf("Creating the EPS file with the CDS of the chosen individual... ");
  fflush(stdout);
  draw_genetic_unit_with_CDS(indiv_main_genome, directory_name);
  printf("OK\n");

  printf("Creating the EPS file with the mRNAs of the chosen individual... ");
  fflush(stdout);
  draw_genetic_unit_with_mRNAs(indiv_main_genome, directory_name);
  printf("OK\n");

#ifdef __REGUL

  for (const auto& prot : rindiv->protein_list_) {
      ((Protein_R*) prot)->reset_concentration();
    }

    rindiv->init_indiv();


  char filename[128];

  snprintf(filename, 127, "%s/best_concentration.csv", directory_name);

  FILE * drawingfile = fopen(filename, "w");

  if (drawingfile == NULL)
    {
      fprintf(stderr, "Error: could not create the file %s\n", filename);
    }


  fprintf(drawingfile, "Generation,Protein,Concentration\n");

  for (const auto& prot : rindiv->protein_list_) {
    ((Protein_R*)prot)->reset_concentration();
  }

    for (const auto& prot : rindiv->protein_list_) {
      printf("Protein List %d -> %f (signal: %d)\n",((Protein_R*)prot)->get_id(),prot->concentration(),((Protein_R*)prot)->is_signal() );
    }


  std::set<int>* eval = exp_manager->exp_s()->get_list_eval_step();
  // i is thus the age of the individual
  for (int16_t i = 1; i <= exp_manager->exp_s()->get_nb_indiv_age(); i++) {
    //Set the concentration of signals for this age
    printf("Prot id from signal_list is ");
    for (auto prot1 : rindiv->signal_list) {
      printf("%d ",prot1.second->get_id());
      prot1.second->set_concentration(0.0);
    }
    printf("\n");

    printf("Prot id from phenotypic target is ");
    for (Protein_R* prot2 : dynamic_cast<const Habitat_R&>(rindiv->habitat()).phenotypic_target(i).signals()) {
      printf("%d ",prot2->get_id());
      rindiv->signal_list[prot2->get_id()]->set_concentration(0.9);
    }
    printf("\n");

    for (auto prot : rindiv->signal_list) {
      printf("before %d prot %d : %f (%p)\n",i,prot.second->get_id(),prot.second->concentration(),prot.second);
    }


    for (const auto& prot : rindiv->protein_list_) {
      if (((Protein_R*)prot)->is_signal()) printf("Protein List %d -> %f (signal: %d) (%p)\n",((Protein_R*)prot)->get_id(),prot->concentration(),((Protein_R*)prot)->is_signal(),prot);
    }

    for (int j = 0; j < exp_manager->exp_s()->get_nb_degradation_step(); j++) {
      rindiv->one_step();
    }


    for (auto prot : rindiv->signal_list) {
      printf("after %d prot %d : %f\n",i,prot.second->get_id(),prot.second->concentration());
    }


    for (const auto& prot : rindiv->protein_list_) {
      if (((Protein_R*)prot)->is_signal()) printf("Protein List %d -> %f (signal: %d)\n",((Protein_R*)prot)->get_id(),prot->concentration(),((Protein_R*)prot)->is_signal() );
    }

    int prot_idx = 0;
    for (const auto& prot : rindiv->protein_list_) {
      if (((Protein_R*)prot)->is_signal()) printf("%d -- %d is signal %f\n",i,((Protein_R*)prot)->get_id(),prot->concentration());
      fprintf(drawingfile, "%d,%d,%lf\n",i,prot_idx,prot->concentration());
      prot_idx++;
    }

    rindiv->update_phenotype();


    printf("Creating the EPS file with the phenotype of the chosen individual at step %d... ",i);
    fflush(stdout);
    draw_phenotype(indiv, dynamic_cast<const Habitat_R&>(rindiv->habitat()).phenotypic_target( i ), directory_name,i);
    printf("OK\n");

  }

  fclose(drawingfile);

#endif



  delete exp_manager;

  return EXIT_SUCCESS;
}







/*!
  \brief

*/
void print_help(char* prog_path)
{
  printf("\n");
  printf("*********************** aevol - Artificial Evolution ******************* \n");
  printf("*                                                                      * \n");
  printf("*                      EPS creation post-treatment program             * \n");
  printf("*                                                                      * \n");
  printf("************************************************************************ \n");
  printf("\n\n");
  printf("This program is Free Software. No Warranty.\n");
  printf("Copyright (C) 2009  LIRIS.\n");
  printf("\n");
#ifdef __REGUL
  printf("Usage : rcreate_eps -h\n");
  printf("or :    rcreate_eps [-i INDEX | -r RANK] -g GENER \n");
#else
  printf("Usage : create_eps -h\n");
  printf("or :    create_eps  [-i INDEX | -r RANK] -g GENER \n");
#endif
  printf("\n");
  printf("This program creates 5 EPS files with the triangles, the positive and negative \n");
  printf("profiles, the phenotype, the CDS, the mRNAs of the chosen individual (the best by default)\n");
  printf("of generation GENER. This program requires at least one population backup file and one\n");
  printf("environment backup file at generation GENER. The EPS files are placed in a \n");
  printf("subdirectory called analysis-generationGENER.\n");
  printf("\n");
  printf("\n");
  printf("\t-h or --help    : Display this help.\n");
  printf("\n");
  printf("\t-i INDEX or --index INDEX : \n");
  printf("\t                  Creates the EPS files for the individual whose\n");
  printf("\t                  index is INDEX. The index must be comprised \n");
  printf("\t                  between 0 and N-1, with N the size of the \n");
  printf("\t                  population at the ending generation. If neither\n");
  printf("\t                  index nor rank are specified, the program creates \n");
  printf("\t                  the EPS files of the best individual. \n");
  printf("\n");
  printf("\t-r RANK or --rank RANK : \n");
  printf("\t                  Creates the EPS files for the individual whose\n");
  printf("\t                  rank is rank. The rank must be comprised \n");
  printf("\t                  between 1 and N, with N the size of the \n");
  printf("\t                  population at the ending generation. If neither\n");
  printf("\t                  index nor rank are specified, the program creates \n");
  printf("\t                  the EPS files of the best individual. \n");
  printf("\n");
  printf("\t-g GENER or --gener GENER : \n");
  printf("\t                  Create the EPS files for the chosen individual of \n");
  printf("\t                  generation GENER. \n");
  printf("\n");
}




void draw_triangles(Individual* indiv, const PhenotypicTarget& target, char * directoryName)
{
  const uint8_t bbsize = 200;  // a4 paper: 595*842
  double margin = 0.1;
  double scalex = 0.8*(1 - 2*margin);
  double scaley = 0.4*(1 - 2*margin);

  char filename[128];
  snprintf(filename, 127, "%s/best_triangles.eps", directoryName);
  FILE * drawingfile = fopen(filename, "w");


  fprintf(drawingfile, "%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(drawingfile, "%%%%BoundingBox: 0 0 %d %d\n", bbsize, bbsize);
  fprintf(drawingfile, "%d %d scale\n", bbsize, bbsize);
  fprintf(drawingfile, "%d %d 8 [100 0 0 -100 0 100]\n",bbsize, bbsize);
  fprintf(drawingfile, "{currentfile 3 100 mul string readhexstring pop} bind\n");


  // -----------------------------
  //  paint neutral zones in grey
  // -----------------------------
  if (target.nb_segments() > 1)
  {
    int16_t nb_segments = target.nb_segments();
    PhenotypicSegment** segments = target.segments();

    for (int16_t i = 0 ; i < nb_segments ; i++)
    {
      if (segments[i]->feature == NEUTRAL)
      {
        fprintf(drawingfile, "%lf 0 moveto\n", margin + scalex * segments[i]->start);
        fprintf(drawingfile, "%lf 1 lineto\n", margin + scalex * segments[i]->start);
        fprintf(drawingfile, "%lf 1 lineto\n", margin + scalex * segments[i]->stop);
        fprintf(drawingfile, "%lf 0 lineto\n", margin + scalex * segments[i]->stop);
        fprintf(drawingfile, "closepath\n");
        fprintf(drawingfile, "0.8 setgray\n");
        fprintf(drawingfile, "fill\n");
        fprintf(drawingfile, "0 setgray\n");
      }
    }
  }

  // -----------
  //  draw axis
  // -----------

  double arrowsize = 0.03;
  double arrowangle = 3.14/6;

  fprintf(drawingfile, "0.001 setlinewidth\n");
  fprintf(drawingfile, "0 0 0 setrgbcolor\n");

  // axis X + arrow
  fprintf(drawingfile, "%lf %lf moveto\n", margin/2, 0.5);
  fprintf(drawingfile, "%lf %lf lineto\n", 1-margin, 0.5);
  fprintf(drawingfile, "stroke\n");
  fprintf(drawingfile, "%lf %lf moveto\n", 1-margin, 0.5);
  fprintf(drawingfile, "%lf %lf lineto\n", 1-margin-arrowsize*cos(arrowangle), 0.5 + arrowsize*sin(arrowangle));
  fprintf(drawingfile, "stroke\n");
  fprintf(drawingfile, "%lf %lf moveto\n", 1-margin, 0.5);
  fprintf(drawingfile, "%lf %lf lineto\n", 1-margin-arrowsize*cos(arrowangle), 0.5 - arrowsize*sin(arrowangle));
  fprintf(drawingfile, "stroke\n");

  // axis Y + arrow
  fprintf(drawingfile, "%lf %lf moveto\n", margin, margin/2);
  fprintf(drawingfile, "%lf %lf lineto\n", margin, 1-margin);
  fprintf(drawingfile, "%lf %lf moveto\n", margin, 1-margin);
  fprintf(drawingfile, "%lf %lf lineto\n", margin-arrowsize*sin(arrowangle), 1 - margin - arrowsize*cos(arrowangle));
  fprintf(drawingfile, "stroke\n");
  fprintf(drawingfile, "%lf %lf moveto\n", margin, 1-margin);
  fprintf(drawingfile, "%lf %lf lineto\n", margin+arrowsize*sin(arrowangle), 1 - margin - arrowsize*cos(arrowangle));
  fprintf(drawingfile, "stroke\n");

  // max degree = 1
  fprintf(drawingfile, "[0.02 0.02] 0 setdash\n");
  fprintf(drawingfile, "%lf %lf moveto\n", margin, 0.5 + 1.0*scaley);
  fprintf(drawingfile, "%lf %lf lineto\n", 1-margin, 0.5 + 1.0*scaley);
  fprintf(drawingfile, "stroke\n");

  // ----------------
  //  draw triangles
  // ----------------

  fprintf(drawingfile,"[ ] 0 setdash\n");

  double h;

  for (auto& gu: indiv->genetic_unit_list_nonconst()) { // should use const version
    for (const auto& prot: gu.protein_list(LEADING)) {
      h = prot.height() * prot.concentration();
      fprintf(drawingfile, "%lf %lf moveto\n", margin, 0.5);
      fprintf(drawingfile, "%lf %lf lineto\n", margin + scalex*(prot.mean() - prot.width()), 0.5);
      fprintf(drawingfile, "%lf %lf lineto\n", margin + scalex*(prot.mean()), 0.5 + scaley*(h));
      fprintf(drawingfile, "%lf %lf lineto\n", margin + scalex*(prot.mean() + prot.width()), 0.5);
      fprintf(drawingfile, "%lf %lf moveto\n", margin + scalex*(1), 0.5);
      fprintf(drawingfile, "stroke\n");
    }

    for (const auto& prot: gu.protein_list(LAGGING)) {
      h = prot.height() * prot.concentration();
      fprintf(drawingfile, "%lf %lf moveto\n", margin, 0.5);
      fprintf(drawingfile, "%lf %lf lineto\n", margin + scalex*(prot.mean() - prot.width()), 0.5);
      fprintf(drawingfile, "%lf %lf lineto\n", margin + scalex*(prot.mean()), 0.5 + scaley*(h));
      fprintf(drawingfile, "%lf %lf lineto\n", margin + scalex*(prot.mean() + prot.width()), 0.5);
      fprintf(drawingfile, "%lf %lf moveto\n", margin + scalex*(1), 0.5);
      fprintf(drawingfile, "stroke\n");
    }
  }



  fprintf(drawingfile,"%%%%EOF\n");
  fclose(drawingfile);

}



void draw_pos_neg_profiles(Individual * indiv, const PhenotypicTarget& target, char * directoryName)
{
  const uint8_t bbsize = 200;  // a4 paper: 595*842
  double margin = 0.1;
  double scale = 0.8*(1 - 2*margin);

  char filename[128];
  snprintf(filename, 127, "%s/best_pos_neg_profiles.eps", directoryName);
  FILE * drawingfile = fopen(filename, "w");


  fprintf(drawingfile, "%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(drawingfile, "%%%%BoundingBox: 0 0 %d %d\n", bbsize, bbsize);
  fprintf(drawingfile, "%d %d scale\n", bbsize, bbsize);
  fprintf(drawingfile, "%d %d 8 [100 0 0 -100 0 100]\n",bbsize, bbsize);
  fprintf(drawingfile, "{currentfile 3 100 mul string readhexstring pop} bind\n");


  // -----------------------------
  //  paint neutral zones in grey
  // -----------------------------
  if (target.nb_segments() > 1)
  {
    int16_t nb_segments = target.nb_segments();
    PhenotypicSegment** segments = target.segments();

    for (int16_t i = 0 ; i < nb_segments ; i++)
    {
      if (segments[i]->feature == NEUTRAL)
      {
        fprintf(drawingfile, "%lf 0 moveto\n", margin + scale * segments[i]->start);
        fprintf(drawingfile, "%lf 1 lineto\n", margin + scale * segments[i]->start);
        fprintf(drawingfile, "%lf 1 lineto\n", margin + scale * segments[i]->stop);
        fprintf(drawingfile, "%lf 0 lineto\n", margin + scale * segments[i]->stop);
        fprintf(drawingfile, "closepath\n");
        fprintf(drawingfile, "0.8 setgray\n");
        fprintf(drawingfile, "fill\n");
        fprintf(drawingfile, "0 setgray\n");
      }
    }
  }


  // -----------
  //  draw axis
  // -----------

  double arrowsize = 0.03;
  double arrowangle = 3.14/6;

  fprintf(drawingfile, "0.001 setlinewidth\n");
  fprintf(drawingfile, "0 0 0 setrgbcolor\n");

  // axis X + arrow
  fprintf(drawingfile, "%lf %lf moveto\n", margin/2, 0.5);
  fprintf(drawingfile, "%lf %lf lineto\n", 1-margin, 0.5);
  fprintf(drawingfile, "stroke\n");
  fprintf(drawingfile, "%lf %lf moveto\n", 1-margin, 0.5);
  fprintf(drawingfile, "%lf %lf lineto\n", 1-margin-arrowsize*cos(arrowangle), 0.5 + arrowsize*sin(arrowangle));
  fprintf(drawingfile, "stroke\n");
  fprintf(drawingfile, "%lf %lf moveto\n", 1-margin, 0.5);
  fprintf(drawingfile, "%lf %lf lineto\n", 1-margin-arrowsize*cos(arrowangle), 0.5 - arrowsize*sin(arrowangle));
  fprintf(drawingfile, "stroke\n");

  // axis Y + arrow
  fprintf(drawingfile, "%lf %lf moveto\n", margin, margin/2);
  fprintf(drawingfile, "%lf %lf lineto\n", margin, 1-margin);
  fprintf(drawingfile, "%lf %lf moveto\n", margin, 1-margin);
  fprintf(drawingfile, "%lf %lf lineto\n", margin-arrowsize*sin(arrowangle), 1 - margin - arrowsize*cos(arrowangle));
  fprintf(drawingfile, "stroke\n");
  fprintf(drawingfile, "%lf %lf moveto\n", margin, 1-margin);
  fprintf(drawingfile, "%lf %lf lineto\n", margin+arrowsize*sin(arrowangle), 1 - margin - arrowsize*cos(arrowangle));
  fprintf(drawingfile, "stroke\n");


  // -----------------------
  //  draw positive profile
  // -----------------------

  fprintf(drawingfile,"[ ] 0 setdash\n");
  fprintf(drawingfile, "0.002 setlinewidth\n");
  fprintf(drawingfile, "%lf %lf moveto\n", margin, 0.5);

  if (indiv->exp_m()->exp_s()->get_fuzzy_flavor() == 0)
    for (const auto& p: ((Fuzzy*)indiv->phenotype_activ())->points())
      fprintf(drawingfile, "%lf %lf lineto\n", margin + scale * p.x, 0.5 + scale * p.y);
  else
    for (int i=0; i < ((HybridFuzzy*)indiv->phenotype_activ())->get_pheno_size(); i++) {
      int xi = (int) ( i / ((HybridFuzzy*)indiv->phenotype_activ())->get_pheno_size());
      fprintf(drawingfile, "%lf %lf lineto\n", margin +
                                               scale * xi, 0.5 + scale *
                                               ((HybridFuzzy*) indiv->phenotype_activ())->points()[i]);
    }
  fprintf(drawingfile, "stroke\n" );

  // -----------------------
  //  draw negative profile
  // -----------------------

  fprintf(drawingfile,"[ ] 0 setdash\n");
  fprintf(drawingfile, "0.002 setlinewidth\n");
  fprintf(drawingfile, "%lf %lf moveto\n", margin, 0.5);

  if (indiv->exp_m()->exp_s()->get_fuzzy_flavor() == 0)
    for (const auto& p: ((Fuzzy*)indiv->phenotype_inhib())->points())
      fprintf( drawingfile, "%lf %lf lineto\n", margin + scale * p.x, 0.5 + scale * p.y);
  else
    for (int i=0; i < ((HybridFuzzy*)indiv->phenotype_inhib())->get_pheno_size(); i++) {
      int xi = (int) ( i / ((HybridFuzzy*)indiv->phenotype_inhib())->get_pheno_size());
      fprintf(drawingfile, "%lf %lf lineto\n", margin +
                                               scale * xi, 0.5 + scale *
                                                                 ((HybridFuzzy*) indiv->phenotype_inhib())->points()[i]);
    }
  fprintf( drawingfile, "stroke\n" );

  fprintf(drawingfile,"%%%%EOF\n");
  fclose(drawingfile);
}




void draw_phenotype(Individual* indiv, const PhenotypicTarget& target, char* directoryName, int generation)
{
  const uint8_t bbsize = 200;  // a4 paper: 595*842
  double margin = 0.1;
  double scale = 0.8*(1 - 2*margin);


  char filename[128];
#ifndef __REGUL
  snprintf(filename, 127, "%s/best_phenotype.csv", directoryName);
#else
  snprintf(filename, 127, "%s/best_phenotype_%d.csv", directoryName,generation);
#endif
  FILE * drawingfile = fopen(filename, "w");

  if (drawingfile == NULL)
    {
      fprintf(stderr, "Error: could not create the file %s\n", filename);
      return;
    }


  fprintf(drawingfile, "Generation,IndividualOrPhenotype,X,Y\n");

  if (indiv->exp_m()->exp_s()->get_fuzzy_flavor() == 0)
    for (const auto& p: ((Fuzzy*)indiv->phenotype_activ())->points())
      fprintf(drawingfile, "%d,0,%lf,%lf\n",generation,p.x,p.y);
  else {
    for (int i=0; i < ((HybridFuzzy*)indiv->phenotype())->get_pheno_size(); i++) {

      float xi =  (i / (float)
                      ((HybridFuzzy*) indiv->phenotype())->get_pheno_size());
      fprintf(drawingfile, "%d,0,%lf,%lf\n",generation,xi,((HybridFuzzy*) indiv->phenotype())->points()[i]);
    }
  }


  // ------------------
  //  draw environment
  // ------------------
  if (indiv->exp_m()->exp_s()->get_fuzzy_flavor() == 0)
    for (const auto& p: ((Fuzzy*)target.fuzzy())->points())
      fprintf(drawingfile, "%d,1,%lf,%lf\n", generation,p.x, p.y);
  else
  for (int i=0; i < ((HybridFuzzy*)target.fuzzy())->get_pheno_size(); i++) {
    float xi =  (i / (float)
                         ((HybridFuzzy*)target.fuzzy())->get_pheno_size());
    fprintf(drawingfile, "%d,1,%lf,%lf\n",generation,xi, ((HybridFuzzy*) target.fuzzy())->points()[i]);
  }

  fclose(drawingfile);

}







void draw_genetic_unit_with_CDS(GeneticUnit* gen_unit, char * directoryName)
{
  const uint8_t bbsize = 200;  // a4 paper: 595*842
  int32_t gen_length = (gen_unit->dna())->length();
  double r = 0.35;
  double scale = 2*M_PI*r/gen_length;

  char filename[128];
  snprintf(filename, 127, "%s/best_genome_with_CDS.eps", directoryName);
  FILE * drawingfile = fopen(filename, "w");

  fprintf(drawingfile, "%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(drawingfile, "%%%%BoundingBox: 0 0 %d %d\n", bbsize, bbsize);
  fprintf(drawingfile, "%d %d scale\n", bbsize, bbsize);
  fprintf(drawingfile, "%d %d 8 [100 0 0 -100 0 100]\n",bbsize, bbsize);
  fprintf(drawingfile, "{currentfile 3 100 mul string readhexstring pop} bind\n");

  // -----------
  //  chromosome
  // -----------

  fprintf(drawingfile, "0.001 setlinewidth\n");
  fprintf(drawingfile, "0 0 0 setrgbcolor\n");
  fprintf(drawingfile, "%lf %lf %lf 0 360 arc\n", 0.5, 0.5, r); // arcn = clockwise arc
  fprintf(drawingfile, "stroke\n");

  // -----------
  //  scale
  // -----------

  double scalesize = 0.15;
  fprintf(drawingfile, "%lf %lf moveto\n", 0.5-scalesize/2, 0.5);
  fprintf(drawingfile, "%lf %lf lineto\n", 0.5+scalesize/2, 0.5);
  fprintf(drawingfile, "stroke\n");
  fprintf(drawingfile, "/Helvetica findfont\n");
  fprintf(drawingfile, "0.035 scalefont\n");
  fprintf(drawingfile, "setfont\n");
  fprintf(drawingfile, "newpath\n");
  fprintf(drawingfile, "%lf %lf moveto\n", 0.5-scalesize/3, 0.52);
  fprintf(drawingfile, "(scale : %.0lf bp) show\n", scalesize/scale);

  // -----------
  //  genes
  // -----------

  int32_t first;
  int32_t last;
  int8_t  layer = 0;
  int8_t  outmost_layer = 1;
  int16_t nb_sect;
  int16_t rho;
  int16_t angle;
  bool    sectors_free;
  int16_t alpha_first;
  int16_t alpha_last;
  int16_t theta_first;
  int16_t theta_last;


  // To handle gene overlaps, we remember where we have aldready drawn
  // something, with a precision of 1 degree. We handle up to 100 layers:
  //  - 50 layers inside the circle (lagging strand),
  //  - 50 layers outside the cricle (leading strand).
  // At first, only one layer is created, we create new layers when we
  // need them.
  bool* occupied_sectors[2][50];
  occupied_sectors[LEADING][0] = new bool[360];
  occupied_sectors[LAGGING][0] = new bool[360];
  for (int16_t angle = 0 ; angle < 360 ; angle++)
  {
    occupied_sectors[LEADING][0][angle] = false;
    occupied_sectors[LAGGING][0][angle] = false;
  }


  // printf("LEADING\n");
  for (const auto& prot: gen_unit->protein_list(LEADING)) {
    first = prot.first_translated_pos();
    last = prot.last_translated_pos();
    // h = prot.height() * prot.concentration();

    alpha_first   = (int16_t) round((double)(360 * first) / (double)gen_length);  //  == sect1 == alphaB
    alpha_last    = (int16_t) round((double)(360 * last)  / (double)gen_length);  //  == sect2 == alphaA
    theta_first   = Utils::mod(90 - alpha_first, 360);  //  == tetaB
    theta_last    = Utils::mod(90 - alpha_last, 360);  //   == tetaA
    if (theta_first == theta_last) theta_first = Utils::mod(theta_first + 1, 360);

    nb_sect = Utils::mod(theta_first - theta_last + 1, 360);


    // Outside the circle, look for the inmost layer that has all the sectors between
    // theta_first and theta_last free
    layer = 0;
    sectors_free = false;
    while (! sectors_free)
    {
      sectors_free = true;
      for (rho = 0 ; rho < nb_sect ; rho++)
      {
        if (occupied_sectors[LEADING][layer][Utils::mod(theta_first - rho, 360)])
        {
          sectors_free = false;
          break;
        }
      }

      if (sectors_free)
      {
        break; // All the needed sectors are free on the current layer
      }
      else
      {
        layer++;

        if ((layer >= outmost_layer) && (layer < 49))
        {
          // Add a new layer (actually, 2 layers, to maintain the symmetry)
          occupied_sectors[LEADING][outmost_layer] = new bool[360];
          occupied_sectors[LAGGING][outmost_layer] = new bool[360];
          for (int16_t angle = 0 ; angle < 360 ; angle++)
          {
            occupied_sectors[LEADING][outmost_layer][angle] = false;
            occupied_sectors[LAGGING][outmost_layer][angle] = false;
          }

          outmost_layer++;
          break; // A new layer is necessarily free, no need to look further
        }
        if (layer == 49)
        {
          // We shall not create a 51th layer, the CDS will be drawn on the
          // layer, probably over another CDS
          break;
        }
      }
    }

    // printf("f %d, l %d, af %d, al %d, tf %d, tl %d, nbsect %d, layer %d\n", first, last, alpha_first, alpha_last, theta_first, theta_last, nb_sect, layer);

    // Mark sectors to be drawn as occupied
    for (rho = 0 ; rho < nb_sect ; rho++)
    {
      occupied_sectors[LEADING][layer][Utils::mod(theta_first - rho, 360)] = true;
    }
    // Mark flanking sectors as occupied
    occupied_sectors[LEADING][layer][Utils::mod(theta_first + 1, 360)] = true;
    occupied_sectors[LEADING][layer][Utils::mod(theta_first - nb_sect, 360)] = true;


    // draw !
    fprintf(drawingfile, "0.018 setlinewidth\n");
    // fprintf(drawingfile, "%lf %lf %lf setrgbcolor\n",  1-(0.8*h/max_height + 0.2), 1-(0.8*h/max_height + 0.2),1-(0.8*h/max_height + 0.2));
    layer++; // index starting at 0 but needed to start at 1

    if (theta_last > theta_first)
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r + layer*0.02, theta_last, 360);
      fprintf(drawingfile, "stroke\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r + layer*0.02, 0, theta_first);
      fprintf(drawingfile, "stroke\n");
    }
    else
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r + layer*0.02, theta_last, theta_first);
      fprintf(drawingfile, "stroke\n");
    }
  }


  // printf("LAGGING\n");
  for (const auto& prot: gen_unit->protein_list(LAGGING)) {
    first = prot.first_translated_pos();
    last = prot.last_translated_pos();
    // h = prot.height() * prot.concentration();

    alpha_first   = (int16_t) round((double)(360 * first) / (double)gen_length);
    alpha_last    = (int16_t) round((double)(360 * last)  / (double)gen_length);
    theta_first   = Utils::mod(90 - alpha_first, 360);
    theta_last    = Utils::mod(90 - alpha_last, 360);
    if (theta_first == theta_last) theta_last = Utils::mod(theta_last + 1, 360);

    nb_sect = Utils::mod(theta_last - theta_first + 1, 360);


    // Inside the circle, look for the inmost layer that has all the sectors between
    // theta_first and theta_last free
    layer = 0;
    sectors_free = false;
    while (! sectors_free)
    {
      sectors_free = true;
      for (rho = 0 ; rho < nb_sect ; rho++)
      {
        if (occupied_sectors[LAGGING][layer][Utils::mod(theta_first + rho, 360)])
        {
          sectors_free = false;
          break;
        }
      }

      if (sectors_free)
      {
        break; // All the needed sectors are free on the current layer
      }
      else
      {
        layer++;

        if ((layer >= outmost_layer) && (layer < 49))
        {
          // Add a new layer (actually, 2 layers, to maintain the symmetry)
          occupied_sectors[LEADING][outmost_layer] = new bool[360];
          occupied_sectors[LAGGING][outmost_layer] = new bool[360];
          for (angle = 0 ; angle < 360 ; angle++)
          {
            occupied_sectors[LEADING][outmost_layer][angle] = false;
            occupied_sectors[LAGGING][outmost_layer][angle] = false;
          }

          outmost_layer++;
          break; // A new layer is necessarily free, no need to look further
        }
        if (layer == 49)
        {
          // We shall not create a 51th layer, the CDS will be drawn on the
          // layer, probably over another CDS
          break;
        }
      }
    }

    // printf("f %d, l %d, af %d, al %d, tf %d, tl %d, nbsect %d, layer %d\n", first, last, alpha_first, alpha_last, theta_first, theta_last, nb_sect, layer);

    // Mark sectors to be drawn as occupied
    for (rho = 0 ; rho < nb_sect ; rho++)
    {
      occupied_sectors[LAGGING][layer][Utils::mod(theta_first + rho, 360)] = true;
    }
    // Mark flanking sectors as occupied
    occupied_sectors[LAGGING][layer][Utils::mod(theta_first - 1, 360)] = true;
    occupied_sectors[LAGGING][layer][Utils::mod(theta_first + nb_sect, 360)] = true;


    // draw !
    fprintf(drawingfile, "0.018 setlinewidth\n");
    // fprintf(drawingfile, "%lf %lf %lf setrgbcolor\n",  1-(0.8*h/max_height + 0.2), 1-(0.8*h/max_height + 0.2),1-(0.8*h/max_height + 0.2));
    layer++; // index starting at 0 but needed to start at 1

    if (theta_first > theta_last)
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r - layer*0.02, theta_first, 360);
      fprintf(drawingfile, "stroke\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r - layer*0.02, 0, theta_last);
      fprintf(drawingfile, "stroke\n");
    }
    else
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r - layer*0.02, theta_first, theta_last);
      fprintf(drawingfile, "stroke\n");
    }
  }


  fprintf(drawingfile,"showpage\n");
  fprintf(drawingfile,"%%%%EOF\n");
  fclose(drawingfile);

  for (layer = 0 ; layer < outmost_layer ; layer++)
  {
    delete occupied_sectors[LEADING][layer];
    delete occupied_sectors[LAGGING][layer];
  }

}




void draw_genetic_unit_with_mRNAs(GeneticUnit* gen_unit, char * directoryName)
{
  const uint8_t bbsize = 200;  // a4 paper: 595*842
  int32_t gen_length = (gen_unit->dna())->length();
  double r = 0.35;
  double scale = 2*M_PI*r/gen_length;

  char filename[128];
  snprintf(filename, 127, "%s/best_genome_with_mRNAs.eps", directoryName);
  FILE * drawingfile = fopen(filename, "w");

  fprintf(drawingfile, "%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(drawingfile, "%%%%BoundingBox: 0 0 %d %d\n", bbsize, bbsize);
  fprintf(drawingfile, "%d %d scale\n", bbsize, bbsize);
  fprintf(drawingfile, "%d %d 8 [100 0 0 -100 0 100]\n",bbsize, bbsize);
  fprintf(drawingfile, "{currentfile 3 100 mul string readhexstring pop} bind\n");

  // -----------
  //  chromosome
  // -----------

  fprintf(drawingfile, "0.001 setlinewidth\n");
  fprintf(drawingfile, "0 0 0 setrgbcolor\n");
  fprintf(drawingfile, "%lf %lf %lf 0 360 arc\n", 0.5, 0.5, r); // arcn = clockwise arc
  fprintf(drawingfile, "stroke\n");

  // -----------
  //  scale
  // -----------

  double scalesize = 0.15;
  fprintf(drawingfile, "%lf %lf moveto\n", 0.5-scalesize/2, 0.5);
  fprintf(drawingfile, "%lf %lf lineto\n", 0.5+scalesize/2, 0.5);
  fprintf(drawingfile, "stroke\n");
  fprintf(drawingfile, "/Helvetica findfont\n");
  fprintf(drawingfile, "0.035 scalefont\n");
  fprintf(drawingfile, "setfont\n");
  fprintf(drawingfile, "newpath\n");
  fprintf(drawingfile, "%lf %lf moveto\n", 0.5-scalesize/3, 0.52);
  fprintf(drawingfile, "(scale : %.0lf bp) show\n", scalesize/scale);

  // -----------
  //  mRNAs
  // -----------

  int32_t first;
  int32_t last;
  int8_t  layer = 0;
  int8_t  outmost_layer = 1;
  int16_t nb_sect;
  int16_t rho;
  int16_t angle;
  bool    sectors_free;
  int16_t alpha_first;
  int16_t alpha_last;
  int16_t theta_first;
  int16_t theta_last;


  // To handle gene overlaps, we remember where we have aldready drawn
  // something, with a precision of 1 degree. We handle up to 100 layers:
  //  - 50 layers inside the circle (lagging strand),
  //  - 50 layers outside the cricle (leading strand).
  // At first, only one layer is created, we create new layers when we
  // need them.
  bool* occupied_sectors[2][50];
  occupied_sectors[LEADING][0] = new bool[360];
  occupied_sectors[LAGGING][0] = new bool[360];
  for (int16_t angle = 0 ; angle < 360 ; angle++)
  {
    occupied_sectors[LEADING][0][angle] = false;
    occupied_sectors[LAGGING][0][angle] = false;
  }



  for (const auto& rna: gen_unit->rna_list()[LEADING]) {
    first = rna.first_transcribed_pos();
    last = rna.last_transcribed_pos();


    alpha_first   = (int16_t) round((double)(360 * first) / (double)gen_length);  //  == sect1 == alphaB
    alpha_last    = (int16_t) round((double)(360 * last)  / (double)gen_length);  //  == sect2 == alphaA
    theta_first   = Utils::mod(90 - alpha_first, 360);  //  == tetaB
    theta_last    = Utils::mod(90 - alpha_last, 360);  //   == tetaA
    if (theta_first == theta_last) theta_first = Utils::mod(theta_first + 1, 360);

    nb_sect = Utils::mod(theta_first - theta_last + 1, 360);


    // Outside the circle, look for the inmost layer that has all the sectors between
    // theta_first and theta_last free
    layer = 0;
    sectors_free = false;
    while (! sectors_free)
    {
      sectors_free = true;
      for (rho = 0 ; rho < nb_sect ; rho++)
      {
        if (occupied_sectors[LEADING][layer][Utils::mod(theta_first - rho, 360)])
        {
          sectors_free = false;
          break;
        }
      }

      if (sectors_free)
      {
        break; // All the needed sectors are free on the current layer
      }
      else
      {
        layer++;

        if ((layer >= outmost_layer) && (layer < 49))
        {
          // Add a new layer (actually, 2 layers, to maintain the symmetry)
          occupied_sectors[LEADING][outmost_layer] = new bool[360];
          occupied_sectors[LAGGING][outmost_layer] = new bool[360];
          for (int16_t angle = 0 ; angle < 360 ; angle++)
          {
            occupied_sectors[LEADING][outmost_layer][angle] = false;
            occupied_sectors[LAGGING][outmost_layer][angle] = false;
          }

          outmost_layer++;
          break; // A new layer is necessarily free, no need to look further
        }
        if (layer == 49)
        {
          // We shall not create a 51th layer, the CDS will be drawn on the
          // layer, probably over another CDS
          break;
        }
      }
    }

    // Mark sectors to be drawn as occupied
    for (rho = 0 ; rho < nb_sect ; rho++)
    {
      occupied_sectors[LEADING][layer][Utils::mod(theta_first - rho, 360)] = true;
    }

    // Mark flanking sectors as occupied
    occupied_sectors[LEADING][layer][Utils::mod(theta_first + 1, 360)] = true;
    occupied_sectors[LEADING][layer][Utils::mod(theta_first - nb_sect, 360)] = true;


    // draw !
    fprintf(drawingfile, "0.018 setlinewidth\n");
    if (rna.is_coding()) fprintf(drawingfile, "0 0 0 setrgbcolor\n");
    else fprintf(drawingfile, "0.7 0.7 0.7 setrgbcolor\n");
    layer++; // index starting at 0 but needed to start at 1

    if (theta_last > theta_first)
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r + layer*0.02, theta_last, 360);
      fprintf(drawingfile, "stroke\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r + layer*0.02, 0, theta_first);
      fprintf(drawingfile, "stroke\n");
    }
    else
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r + layer*0.02, theta_last, theta_first);
      fprintf(drawingfile, "stroke\n");
    }
  }



  for (const auto& rna: gen_unit->rna_list()[LAGGING]) {
    first = rna.first_transcribed_pos();
    last = rna.last_transcribed_pos();


    alpha_first   = (int16_t) round((double)(360 * first) / (double)gen_length);
    alpha_last    = (int16_t) round((double)(360 * last)  / (double)gen_length);
    theta_first   = Utils::mod(90 - alpha_first, 360);
    theta_last    = Utils::mod(90 - alpha_last, 360);
    nb_sect = Utils::mod(alpha_first - alpha_last + 1,  360);


    // Inside the circle, look for the inmost layer that has all the sectors between
    // theta_first and theta_last free
    layer = 0;
    sectors_free = false;
    while (! sectors_free)
    {
      sectors_free = true;
      for (rho = 0 ; rho < nb_sect ; rho++)
      {
        if (occupied_sectors[LAGGING][layer][Utils::mod(theta_first + rho, 360)])
        {
          sectors_free = false;
          break;
        }
      }

      if (sectors_free)
      {
        break; // All the needed sectors are free on the current layer
      }
      else
      {
        layer++;

        if ((layer >= outmost_layer) && (layer < 49))
        {
          // Add a new layer (actually, 2 layers, to maintain the symmetry)
          occupied_sectors[LEADING][outmost_layer] = new bool[360];
          occupied_sectors[LAGGING][outmost_layer] = new bool[360];
          for (angle = 0 ; angle < 360 ; angle++)
          {
            occupied_sectors[LEADING][outmost_layer][angle] = false;
            occupied_sectors[LAGGING][outmost_layer][angle] = false;
          }

          outmost_layer++;
          break; // A new layer is necessarily free, no need to look further
        }
        if (layer == 49)
        {
          // We shall not create a 51th layer, the CDS will be drawn on the
          // layer, probably over another CDS
          break;
        }
      }
    }

    // Mark sectors to be drawn as occupied
    for (rho = 0 ; rho < nb_sect ; rho++)
    {
      occupied_sectors[LAGGING][layer][Utils::mod(theta_first + rho, 360)] = true;
    }

    // Mark flanking sectors as occupied
    occupied_sectors[LAGGING][layer][Utils::mod(theta_first - 1, 360)] = true;
    occupied_sectors[LAGGING][layer][Utils::mod(theta_first + nb_sect, 360)] = true;


    // draw !
    fprintf(drawingfile, "0.018 setlinewidth\n");
    if (rna.is_coding()) fprintf(drawingfile, "0 0 0 setrgbcolor\n");
    else fprintf(drawingfile, "0.7 0.7 0.7 setrgbcolor\n");
    layer++; // index starting at 0 but needed to start at 1

    if (theta_first > theta_last)
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r - layer*0.02, theta_first, 360);
      fprintf(drawingfile, "stroke\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r - layer*0.02, 0, theta_last);
      fprintf(drawingfile, "stroke\n");
    }
    else
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r - layer*0.02, theta_first, theta_last);
      fprintf(drawingfile, "stroke\n");
    }
  }


  fprintf(drawingfile,"showpage\n");
  fprintf(drawingfile,"%%%%EOF\n");
  fclose(drawingfile);

  for (layer = 0 ; layer < outmost_layer ; layer++)
  {
    delete occupied_sectors[LEADING][layer];
    delete occupied_sectors[LAGGING][layer];
  }

}
