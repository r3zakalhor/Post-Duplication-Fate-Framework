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
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cerrno>
#include <getopt.h>
#include <sys/stat.h>  // for the permission symbols used with mkdir
#include <errno.h>
#include <HybridFuzzy.h>
#include <fstream>

// =================================================================
//                            Project Files
// =================================================================
#include "aevol.h"
#include "Rna_7.h"
#include "Protein_7.h"
#include "Abstract_Metadata.h"
#include "AbstractFuzzy_7.h"
#include "Discrete_Double_Fuzzy.h"
#include "List_Metadata.h"
#include "Vector_Fuzzy.h"
using namespace aevol;

// ========
//  TO DO
// ========
//
//  * option --color ?
//  * Raevol-specific output (EPS file with the network) ?

// Command-line option variables
static int64_t timestep = -1;
static int32_t indiv_index = -1;
static int32_t indiv_rank = -1;

// Helper functions
void print_help(char* prog_path);
void interpret_cmd_line_options(int argc, char* argv[]);

// The height of each triangle is proportional to the product c*m,
// where c is the concentration of the protein and m its intrinsic efficacy
// (depending on its aminoacid sequence).
// In the case of Raevol, the concentration used here is the final one,
// i.e. the one reached after all the time steps of the lifetime.
// If a coding sequence has several promoters, only one triangle is drawn.
void draw_triangles(Individual_7* indiv, AbstractFuzzy_7* target,
                    char* directoryName);


// In the case of Raevol, the profile is drawn using the final concentrations
// of the proteins, i.e. the ones reached after all the time steps of the
// lifetime.
void draw_pos_neg_profiles(Individual_7* indiv, AbstractFuzzy_7* target,
                           char* directoryName);


// In the case of Raevol, the phenotype is drawn using the final concentrations
// of the proteins, i.e. the ones reached after all the time steps of the
// lifetime.
void draw_phenotype(Individual_7* indiv, AbstractFuzzy_7* target,
                    char* directoryName, int32_t timestep, 
                    int32_t env_id = -1, bool single_env = true, SIMD_PhenotypicTargetHandler_R* handler = nullptr,
                    FuzzyFactory_7* fuzzy_factory = nullptr, ExpManager* exp_m = nullptr,
                    DnaFactory* dna_factory = nullptr);


// The chromosome is drawn as a circle. Coding sequences are represented by
// arcs starting at the start codon and ending at the stop codon.
// Coding sequences that are on the leading (resp. lagging) strand are
// drawn outside (resp. inside) the circle. If a coding sequence has several
// promoters, only one arc is drawn.
void draw_genetic_unit_with_CDS(Individual_7* indiv, char* directoryName, int64_t timestep);


// The chromosome is drawn as a circle. Transcribed sequences are represented
// by arcs starting at the first transcribed position and ending at the last
// transcribed position. mRNAs that are on the leading (resp. lagging) strand
// are drawn outside (resp. inside) the circle.
// mRNAs that include at least one coding sequence are black,
// the others are gray.
void draw_genetic_unit_with_mRNAs(Individual_7* indiv, char* directoryName, int64_t timestep);


int main(int argc, char* argv[]) {
  interpret_cmd_line_options(argc, argv);

  printf("Creating eps files for generation %" PRId64 "...\n", timestep);

  // =================================================================
  //                       Read the backup file
  // =================================================================
  Individual*  indiv;

  // Load the simulation
  ExpManager_7::standalone_simd = true;
  ExpManager* exp_manager = new ExpManager();
  exp_manager->load(timestep, true, false);

  if (indiv_index == -1 && indiv_rank == -1) {
    indiv = exp_manager->best_indiv();
  }
  else {
    // TODO <david.parsons@inria.fr> tmp disabled
//     if (indiv_rank != -1) {
//       indiv = new Individual(*exp_manager->indiv_by_rank(indiv_rank), false);
//     }
//     else {
//       indiv = new Individual(*exp_manager->indiv_by_id(indiv_index), false);
//     }
  }

  // Init Factory (Fuzzy/Dna)
  DnaFactory* dna_factory_ = new DnaFactory(DnaFactory_Policy::FIRSTFIT,32,5000);
  FuzzyFactory_7* fuzzy_factory_ = new FuzzyFactory_7(exp_manager->exp_s()->get_fuzzy_flavor(),exp_manager->nb_indivs()*4,
                        exp_manager->world()->phenotypic_target_handler()->sampling());

  // The constructor of the exp_manager has read the genomes of the individuals
  // and located their promoters, but has not performed the translation nor the
  // phenotype computation. We must do it now.
  // However, as the individuals in the backups are sorted, we don't need to
  // evaluate all the individuals, only the one we are interested in
  Individual_7* indiv_7 = new Individual_7(exp_manager, indiv->w_max(),dna_factory_,fuzzy_factory_);
  indiv_7->dna_ = dna_factory_->get_dna(indiv->genetic_unit_seq_length(0));
  indiv_7->dna_->set_indiv(indiv->genetic_unit(0).dna(),dna_factory_);
  indiv_7->dna_->set_indiv(indiv_7);
  indiv_7->indiv_id = 0;
  indiv_7->parent_id = 0;

  SIMD_PhenotypicTargetHandler_R* pth_target = nullptr;
  
  #ifdef __REGUL
  pth_target = new SIMD_PhenotypicTargetHandler_R(exp_manager->exp_m_7_->phenotypic_target_handler_,
                fuzzy_factory_,exp_manager);

  pth_target->var_prng_ = std::make_shared<JumpingMT>(
    exp_manager->exp_m_7_->phenotypic_target_handler_->var_prng_->random(100000000));
  pth_target->ApplyVariation();

  #endif
  exp_manager->exp_m_7_->evaluate(indiv_7,indiv->w_max(),exp_manager->selection_pressure(),pth_target);

  // =================================================================
  //                      Create the EPS files
  // =================================================================
  char directory_name[64];
  snprintf(directory_name, 63, "analysis-generation_" TIMESTEP_FORMAT,
           timestep);

  // Check whether the directory already exists and is writable
  if (access(directory_name, F_OK) == 0) {
//       struct stat status;
//       stat(directory_name, &status);
//       if (status.st_mode & S_IFDIR) cout << "The directory exists." << endl;
//       else cout << "This path is a file." << endl;

    if (access(directory_name, X_OK | W_OK) != 0) {
      fprintf(stderr, "Error: cannot enter or write in directory %s.\n", directory_name);
      exit(EXIT_FAILURE);
    }
  }
  else {
    // Create the directory with permissions : rwx r-x r-x
    if (mkdir(directory_name,
              S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH) != 0) {
      fprintf(stderr, "Error: cannot create directory %s.\n", directory_name);
      exit(EXIT_FAILURE);
    }
  }



  // =================================================================
  //                  Write the data in the EPS files
  // =================================================================

#ifndef __REGUL
  printf("Creating the EPS file with the triangles of the chosen individual... ");
  fflush(stdout);
  draw_triangles(indiv_7, exp_manager->exp_m_7_->target, directory_name);
  printf("OK\n");

  //TODO : Should be ported

  // printf("Creating the EPS file with the positive and negatives profiles of the chosen individual... ");
  // fflush(stdout);
  // draw_pos_neg_profiles(indiv_7, exp_manager->exp_m_7_->target, directory_name);
  // printf("OK\n");

  printf("Creating the EPS file with the phenotype of the chosen individual... ");
  fflush(stdout);
  draw_phenotype(indiv_7, exp_manager->exp_m_7_->target, directory_name,timestep,-1,true,nullptr,fuzzy_factory_,exp_manager,dna_factory_);
  printf("OK\n");

#else
  printf("Creating the EPS file with the phenotype of the chosen individual... ");
  fflush(stdout);

  std::ofstream mullerfile;
  mullerfile.open("stats/phenotype_visu.csv",std::ofstream::trunc);
  mullerfile<<"generation,iteration,single_env,type,x_start,y_start,fill_color,env_id"<<std::endl;
  mullerfile.close();
  
  // for (int32_t env_id = 0; env_id < exp_manager->exp_m_7_->phenotypic_target_handler_->nb_env_; env_id++) {
  //   draw_phenotype(indiv_7, exp_manager->exp_m_7_->phenotypic_target_handler_->targets_fuzzy_by_id_[env_id], directory_name,timestep,env_id,true,
  //                   exp_manager->exp_m_7_->phenotypic_target_handler_,fuzzy_factory_,exp_manager,dna_factory_);
  // }


  std::ofstream concentrationfile;
  concentrationfile.open("stats/concentration_visu.csv",std::ofstream::trunc);
  concentrationfile<<"generation,iteration,age,prot_id,concentration,signal"<<std::endl;
  concentrationfile.close();

  draw_phenotype(indiv_7, exp_manager->exp_m_7_->phenotypic_target_handler_->targets_fuzzy_by_id_[0], directory_name,timestep,0,false,
                    exp_manager->exp_m_7_->phenotypic_target_handler_,fuzzy_factory_,exp_manager,dna_factory_);

  printf("OK\n");
#endif

  printf("Creating the EPS file with the CDS of the chosen individual... ");
  fflush(stdout);
  draw_genetic_unit_with_CDS(indiv_7, directory_name,timestep);
  printf("OK\n");

  printf("Creating the EPS file with the mRNAs of the chosen individual... ");
  fflush(stdout);
  draw_genetic_unit_with_mRNAs(indiv_7, directory_name,timestep);
  printf("OK\n");



// #ifdef __REGUL


//   std::set<int>* eval = exp_manager->exp_s()->get_list_eval_step();
//   // i is thus the age of the individual
//   for (int16_t i = 1; i <= exp_manager->exp_s()->get_nb_indiv_age(); i++) {
//     //Set the concentration of signals for this age
//     for (auto prot1 : rindiv->signal_list) {
//       prot1.second->set_concentration(0.0);
//     }

//     for (Protein_R* prot2 : dynamic_cast<const Habitat_R&>(rindiv->habitat()).phenotypic_target(i).signals()) {
//       rindiv->signal_list[prot2->get_id()]->set_concentration(0.9);
//     }


//     for (int j = 0; j < exp_manager->exp_s()->get_nb_degradation_step(); j++) {
//       rindiv->one_step();
//     }

//     rindiv->update_phenotype();


//     printf("Creating the EPS file with the phenotype of the chosen individual at step %d... ",i);
//     fflush(stdout);
//     draw_phenotype(indiv, dynamic_cast<const Habitat_R&>(rindiv->habitat()).phenotypic_target( i ), directory_name);
//     printf("OK\n");

//   }



// #endif



  delete exp_manager;

  return EXIT_SUCCESS;
}


void draw_triangles(Individual_7* indiv, AbstractFuzzy_7* target,
                    char* directoryName) {
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
  indiv->metadata_->protein_begin();
  for (int j = 0; j < indiv->metadata_->proteins_count(); j++) {
    Protein_7* prot = indiv->metadata_->protein_next();
    if (prot != nullptr) {
      if (prot->is_init_) {
        if (prot->h > 0) {
          h = prot->h * prot->e;
          fprintf(drawingfile, "%lf %lf moveto\n", margin, 0.5);
          fprintf(drawingfile, "%lf %lf lineto\n", margin + scalex*(prot->m - prot->w), 0.5);
          fprintf(drawingfile, "%lf %lf lineto\n", margin + scalex*(prot->m), 0.5 + scaley*(h));
          fprintf(drawingfile, "%lf %lf lineto\n", margin + scalex*(prot->m + prot->w), 0.5);
          fprintf(drawingfile, "%lf %lf moveto\n", margin + scalex*(1), 0.5);
          fprintf(drawingfile, "stroke\n");
        } else {
          h = prot->h * prot->e;
          fprintf(drawingfile, "%lf %lf moveto\n", margin, 0.5);
          fprintf(drawingfile, "%lf %lf lineto\n", margin + scalex*(prot->m - prot->w), 0.5);
          fprintf(drawingfile, "%lf %lf lineto\n", margin + scalex*(prot->m), 0.5 + scaley*(h));
          fprintf(drawingfile, "%lf %lf lineto\n", margin + scalex*(prot->m + prot->w), 0.5);
          fprintf(drawingfile, "%lf %lf moveto\n", margin + scalex*(1), 0.5);
          fprintf(drawingfile, "stroke\n");
        }
      }
    }
  }

  fprintf(drawingfile,"%%%%EOF\n");
  fclose(drawingfile);

}


void draw_pos_neg_profiles(Individual_7* indiv, AbstractFuzzy_7* target,
                           char* directoryName) {
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

  // if (indiv->exp_m()->exp_s()->get_fuzzy_flavor() == 1)
  //   for (const auto& p: ((Vector_Fuzzy*)indiv->phenotype_activ())->points())
  //     fprintf(drawingfile, "%lf %lf lineto\n", margin + scale * p.x, 0.5 + scale * p.y);
  // else if (indiv->exp_m()->exp_s()->get_fuzzy_flavor() == 3)
  //   for (int i=0; i < ((HybridFuzzy*)indiv->phenotype_activ())->get_pheno_size(); i++) {
  //     int xi = (int) ( i / ((HybridFuzzy*)indiv->phenotype_activ())->get_pheno_size());
    
  //     fprintf(drawingfile, "%lf %lf lineto\n", margin +
  //                                              scale * xi, 0.5 + scale *
  //                                              ((HybridFuzzy*) indiv->phenotype_activ())->points()[i]);
  //   }
  fprintf(drawingfile, "stroke\n" );

  // -----------------------
  //  draw negative profile
  // -----------------------

  fprintf(drawingfile,"[ ] 0 setdash\n");
  fprintf(drawingfile, "0.002 setlinewidth\n");
  fprintf(drawingfile, "%lf %lf moveto\n", margin, 0.5);

  // if (indiv->exp_m()->exp_s()->get_fuzzy_flavor() == 0)
  //   for (const auto& p: ((Fuzzy*)indiv->phenotype_inhib())->points())
  //     fprintf( drawingfile, "%lf %lf lineto\n", margin + scale * p.x, 0.5 + scale * p.y);
  // else
  //   for (int i=0; i < ((HybridFuzzy*)indiv->phenotype_inhib())->get_pheno_size(); i++) {
  //     int xi = (int) ( i / ((HybridFuzzy*)indiv->phenotype_inhib())->get_pheno_size());
  //     fprintf(drawingfile, "%lf %lf lineto\n", margin +
  //                                              scale * xi, 0.5 + scale *
  //                                                                ((HybridFuzzy*) indiv->phenotype_inhib())->points()[i]);
  //   }
  fprintf( drawingfile, "stroke\n" );

  fprintf(drawingfile,"%%%%EOF\n");
  fclose(drawingfile);
}


void draw_phenotype(Individual_7* indiv, AbstractFuzzy_7* target,
                    char* directoryName, int32_t timestep, int32_t env_id, bool single_env, SIMD_PhenotypicTargetHandler_R* handler,
                    FuzzyFactory_7* fuzzy_factory, ExpManager* exp_m, DnaFactory* dna_factory) {
  const uint8_t bbsize = 200;  // a4 paper: 595*842
  double margin = 0.1;
  double scale = 0.8*(1 - 2*margin);

  double w_max = exp_m->best_indiv()->w_max();
  double selection_pressure = exp_m->selection_pressure();

  char filename[128];
#ifndef __REGUL
  snprintf(filename, 127, "%s/best_phenotype.eps", directoryName);
#else
  snprintf(filename, 127, "%s/best_phenotype.eps", directoryName);
#endif
  FILE * drawingfile = fopen(filename, "w");

  if (drawingfile == NULL)
    {
      fprintf(stderr, "Error: could not create the file %s\n", filename);
      return;
    }


#ifdef __REGUL
  std::ofstream mullerfile;
  mullerfile.open("stats/phenotype_visu.csv",std::ofstream::app);
#else
  std::ofstream mullerfile;
  mullerfile.open("stats/phenotype_visu.csv",std::ofstream::trunc);
  mullerfile<<"generation,type,x_start,y_start,fill_color,env_id"<<std::endl;
#endif

  #ifdef __REGUL
  std::ofstream concentrationfile;
  concentrationfile.open("stats/concentration_visu.csv",std::ofstream::app);
  #endif
  
  fprintf(drawingfile, "%%!PS-Adobe-3.0 EPSF-3.0\n");
  fprintf(drawingfile, "%%%%BoundingBox: 0 0 %d %d\n", bbsize, bbsize);
  fprintf(drawingfile, "%d %d scale\n", bbsize, bbsize);
  fprintf(drawingfile, "%d %d 8 [100 0 0 -100 0 100]\n",bbsize, bbsize);
  fprintf(drawingfile, "{currentfile 3 100 mul string readhexstring pop} bind\n");
  // -----------
  //  draw axis
  // -----------

  double arrowsize = 0.03;
  double arrowangle = 3.14/6;

  fprintf(drawingfile, "0.001 setlinewidth\n");
  fprintf(drawingfile, "0 0 0 setrgbcolor\n");

  // axis X + arrow
  fprintf(drawingfile, "%lf %lf moveto\n", margin/2, margin);
  fprintf(drawingfile, "%lf %lf lineto\n", 1-margin, margin);
  fprintf(drawingfile, "stroke\n");
  fprintf(drawingfile, "%lf %lf moveto\n", 1-margin, margin);
  fprintf(drawingfile, "%lf %lf lineto\n", 1-margin-arrowsize*cos(arrowangle), margin + arrowsize*sin(arrowangle));
  fprintf(drawingfile, "stroke\n");
  fprintf(drawingfile, "%lf %lf moveto\n", 1-margin, margin);
  fprintf(drawingfile, "%lf %lf lineto\n", 1-margin-arrowsize*cos(arrowangle), margin - arrowsize*sin(arrowangle));
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
  fprintf(drawingfile, "%lf %lf moveto\n", margin, margin + 1*scale);
  fprintf(drawingfile, "%lf %lf lineto\n", 1-margin, margin + 1*scale);
  fprintf(drawingfile, "stroke\n");


  // ----------------
  //  draw phenotype
  // ----------------

  fprintf( drawingfile,"[ ] 0 setdash\n" );
  fprintf( drawingfile, "0.002 setlinewidth\n" );
  fprintf( drawingfile, "%lf %lf moveto\n", margin, margin);

#ifdef __REGUL
  SIMD_PhenotypicTargetHandler_R* pth_target = new SIMD_PhenotypicTargetHandler_R(handler,fuzzy_factory,exp_m);

  pth_target->var_prng_ = std::make_shared<JumpingMT>(handler->var_prng_->random(100000000));
  
  if (single_env) {
    pth_target->set_single_env(env_id);
    Individual_7* cloned = new Individual_7(exp_m,indiv,dna_factory,fuzzy_factory,true);
          
    exp_m->exp_m_7_->evaluate(cloned,w_max,selection_pressure,pth_target);
    

    for (int i = 0; i <  ((Discrete_Double_Fuzzy*) cloned->phenotype)->length() - 1; i++) {
        double xi_start =  (i / (double)
                        ((Discrete_Double_Fuzzy*) cloned->phenotype)->length());
        double xi_end =  ((i+1) / (double)
                        ((Discrete_Double_Fuzzy*) cloned->phenotype)->length());
        // fprintf(drawingfile, "%lf %lf lineto\n", margin +
        //                                          scale * xi, margin + scale *
        //                                                               ((Discrete_Double_Fuzzy*) cloned->phenotype)->points()[i]);
        mullerfile<<timestep<<","<<"0"<<","<<"1,phenotype"<<","<<xi_start<<","<<((Discrete_Double_Fuzzy*) cloned->phenotype)->points()[i]<<","<<i<<","<<env_id<<std::endl;
        mullerfile<<timestep<<","<<"0"<<","<<"1,phenotype"<<","<<xi_end<<","<<((Discrete_Double_Fuzzy*) cloned->phenotype)->points()[i+1]<<","<<i<<","<<env_id<<std::endl;
    }
    fprintf( drawingfile, "stroke\n" );


    // ------------------
    //  draw environment
    // ------------------
    fprintf( drawingfile,"[ ] 0 setdash\n" );
    fprintf( drawingfile, "0.001 setlinewidth\n" );
    fprintf( drawingfile, "%lf %lf moveto\n", margin, margin);
    
    for (int i=0; i < ((Discrete_Double_Fuzzy*)target)->length() - 1; i++) {

        double xi_start =  (i / (double)
                        ((Discrete_Double_Fuzzy*)target)->length());
        double xi_end =  ((i+1) / (double)
                        ((Discrete_Double_Fuzzy*)target)->length());
        // fprintf(drawingfile, "%lf %lf lineto\n", margin +
        //                                          scale * xi, margin + scale *
        //                                                               ((Discrete_Double_Fuzzy*) target)->points()[i]);


        mullerfile<<timestep<<","<<"0"<<","<<"1,env"<<","<<xi_start<<","<<((Discrete_Double_Fuzzy*) target)->points()[i]<<","<<i<<","<<env_id<<std::endl;
        mullerfile<<timestep<<","<<"0"<<","<<"1,env"<<","<<xi_end<<","<<((Discrete_Double_Fuzzy*) target)->points()[i+1]<<","<<i<<","<<env_id<<std::endl;
        //((Discrete_Double_Fuzzy*) target)->print();
    }
    fprintf( drawingfile, "stroke\n" );



    fprintf(drawingfile,"%%%%EOF\n");
    fclose(drawingfile);
  } else {

    for (int repeat = 0; repeat < 10; repeat++) {
      pth_target->ApplyVariation();
      Individual_7* cloned = new Individual_7(exp_m,indiv,dna_factory,fuzzy_factory,false);
      exp_m->exp_m_7_->start_stop_RNA(cloned);
      exp_m->exp_m_7_->compute_RNA(cloned);
      exp_m->exp_m_7_->start_protein(cloned);
      exp_m->exp_m_7_->compute_protein(cloned);
      exp_m->exp_m_7_->translate_protein(cloned, w_max);

      cloned->fitness_by_env_id_ = new double[pth_target->nb_indiv_age_];
      cloned->metaerror_by_env_id_ = new double[pth_target->nb_indiv_age_];

      cloned->metadata_->protein_begin();
      for (int j = 0;
          j < cloned->metadata_->proteins_count(); j++) {
        Protein_7* prot =
            cloned->metadata_->protein_next();
        if (!prot->signal_) {
          if (prot->is_init_) {
            prot->e = prot->initial_e_;
          }
        }
      }

      // printf("Source %d\n",indiv->metadata_->proteins_count());

      ((List_Metadata*)cloned->metadata_)->add_inherited_proteins(indiv);

      cloned->metaerror = 0;
      exp_m->exp_m_7_->compute_network(cloned, selection_pressure,pth_target);
      cloned->metaerror = 0;

      // Save network
      int32_t network_id = 0;

      std::ofstream nodefile;
      nodefile.open("stats/node_graph.csv",std::ofstream::trunc);
      nodefile<<"generation"<<","<<"id"<<","<<"type"<<std::endl;

      cloned->metadata_->rna_begin();
      for (int i = 0; i < cloned->metadata_->rna_count(); i++) {
        Rna_7* rna = cloned->metadata_->rna_next();
        if (rna != nullptr) {
          if (rna->is_coding_) {
            rna->network_id = network_id++;
            nodefile<<timestep<<","<<rna->network_id<<",rna"<<std::endl;
          }
        }
      }

      cloned->metadata_->protein_begin();
      for (int j = 0; j < cloned->metadata_->proteins_count(); j++) {
        Protein_7* prot = cloned->metadata_->protein_next();
        if (prot != nullptr) {
          if (prot->is_init_ || prot->signal_) {
            prot->network_id = network_id++;
            if (prot->signal_)
              nodefile<<timestep<<","<<prot->network_id<<",signal"<<std::endl;
            else
              nodefile<<timestep<<","<<prot->network_id<<",prot"<<std::endl;
          }
        }
      }

      nodefile.close();

      std::ofstream edgefile;
      edgefile.open("stats/edge_graph.csv",std::ofstream::trunc);
      edgefile<<"generation"<<","<<"src"<<","<<"destination"<<","<<"value"<<","<<"type"<<std::endl;

      cloned->metadata_->rna_begin();
      for (int i = 0; i < cloned->metadata_->rna_count(); i++) {
        Rna_7* rna = cloned->metadata_->rna_next();
        if (rna != nullptr) {
          if (rna->is_coding_) {

            int32_t enhancer_position = rna->enhancer_position(cloned->dna_->length());
            int32_t operator_position = rna->operator_position(cloned->dna_->length());

            double enhance,operate;
            cloned->metadata_->protein_begin();
            for (int j = 0; j < cloned->metadata_->proteins_count(); j++) {
              Protein_7* prot =
                  cloned->metadata_->protein_next();
              if (prot != nullptr) {
                if (prot->is_init_ || prot->signal_) {
                  enhance = rna->affinity_with_protein(
                      enhancer_position, prot, cloned,
                      exp_m);
                  operate = rna->affinity_with_protein(
                      operator_position, prot, cloned,
                      exp_m);

                  if (enhance != 0.0) {
                    edgefile<<timestep<<","<<prot->network_id<<","<<rna->network_id<<","<<enhance<<","<<"enhance"<<std::endl;
                  } else if (operate != 0.0) {
                    edgefile<<timestep<<","<<prot->network_id<<","<<rna->network_id<<","<<operate<<","<<"operate"<<std::endl;
                  }
                }
              }
            }
          }
        }
      }

      cloned->metadata_->rna_begin();
      for (int i = 0; i < cloned->metadata_->rna_count(); i++) {
        Rna_7* rna = cloned->metadata_->rna_next();
        if (rna != nullptr) {
          if (rna->is_coding_) {
            for (auto&& prot : rna->protein_list_) {
              if (prot != nullptr) {
                if (prot->is_init_) {
                  edgefile<<timestep<<","<<rna->network_id<<","<<prot->network_id<<","<<"1"<<","<<"produce"<<std::endl;
                }
              }
            }
          }
        }
      }

      // indiv->metadata_->protein_begin();
      // for (int j = 0; j < indiv->metadata_->proteins_count(); j++) {
      //   Protein_7* prot =
      //       indiv->metadata_->protein_next();
      //   if (!prot->signal_) {
      //     if (prot->is_init_) {
      //       for (auto rna: prot->rna_list_) {
      //         edgefile<<timestep<<","<<rna->network_id<<","<<prot->network_id<<","<<"1"<<","<<"produce"<<std::endl;
      //       }
      //     }
      //   }
      // }
      edgefile.close();

      if (pth_target->var_method_ == ONE_AFTER_ANOTHER) {
        for (int16_t env_i = 0; env_i < pth_target->nb_env_; env_i++) {

          //Set the concentration of signals for this age
          for (auto prot1: ((List_Metadata*)indiv->metadata_)->signal_proteins_) {
            prot1->e = 0;
          }

          for (auto prot_id : pth_target->env_signals_list_[env_i]) {
            ((List_Metadata*)indiv->metadata_)->signal_proteins_[prot_id]->e = 0.9;
          }

          for (int16_t i = 0; i < 10; i++) {

            for (int j = 0; j < 10; j++) {
              exp_m->exp_m_7_->update_network(indiv,selection_pressure);
            }

            // If we have to evaluate the individual at this age
            exp_m->exp_m_7_->evaluate_network(indiv,selection_pressure,env_i,pth_target);

            for (int j = 0; j <  ((Discrete_Double_Fuzzy*) cloned->phenotype)->length() - 1; j++) {
              double xi_start =  (j / (double)
                              ((Discrete_Double_Fuzzy*) cloned->phenotype)->length());
              double xi_end =  ((j+1) / (double)
                              ((Discrete_Double_Fuzzy*) cloned->phenotype)->length());
              mullerfile<<timestep<<","<<repeat<<","<<"0,phenotype"<<","<<xi_start<<","<<((Discrete_Double_Fuzzy*) cloned->phenotype)->points()[j]<<","<<pth_target->list_env_id_[i]<<","<<i<<std::endl;
              mullerfile<<timestep<<","<<repeat<<","<<"0,phenotype"<<","<<xi_end<<","<<((Discrete_Double_Fuzzy*) cloned->phenotype)->points()[j+1]<<","<<pth_target->list_env_id_[i]<<","<<i<<std::endl;
            }

            for (int j=0; j < ((Discrete_Double_Fuzzy*)pth_target->targets_fuzzy_[i])->length() - 1; j++) {
              double xi_start =  (j / (double)
                              ((Discrete_Double_Fuzzy*)pth_target->targets_fuzzy_[i])->length());
              double xi_end =  ((j+1) / (double)
                              ((Discrete_Double_Fuzzy*)pth_target->targets_fuzzy_[i])->length());
              mullerfile<<timestep<<","<<repeat<<","<<"0,env"<<","<<xi_start<<","<<((Discrete_Double_Fuzzy*) pth_target->targets_fuzzy_[i])->points()[j]<<","<<pth_target->list_env_id_[i]<<","<<i<<std::endl;
              mullerfile<<timestep<<","<<repeat<<","<<"0,env"<<","<<xi_end<<","<<((Discrete_Double_Fuzzy*) pth_target->targets_fuzzy_[i])->points()[j+1]<<","<<pth_target->list_env_id_[i]<<","<<i<<std::endl;
            }   
// Save concentration
          cloned->metadata_->protein_begin();
          int32_t id = 0;
          for (int j = 0;
              j < cloned->metadata_->proteins_count(); j++) {
            Protein_7* prot =
                cloned->metadata_->protein_next();
            if (prot->is_init_ || prot->signal_) {
                concentrationfile<<timestep<<","<<repeat<<","<<i<<","<<id<<","<<prot->e<<","<<(prot->signal_? 1 : 0)<<std::endl;
                id++;
            }
          }
          }
        }

        exp_m->exp_m_7_->finalize_network(indiv,selection_pressure);
      } else {
        std::set<int>* eval = exp_m->exp_s()->get_list_eval_step();
        int32_t prev_env_id = pth_target->list_env_id_[0];

        for (int16_t i = 0; i < pth_target->nb_indiv_age_; i++) {
          //Set the concentration of signals for this age
          for (auto prot1: ((List_Metadata*)cloned->metadata_)->signal_proteins_) {
            prot1->e = 0;
          }

          for (auto prot_id : pth_target->env_signals_list_[pth_target->list_env_id_[i]]) {
            ((List_Metadata*)cloned->metadata_)->signal_proteins_[prot_id]->e = 0.9;
          }

          for (int j = 0; j < exp_m->exp_s()->get_nb_degradation_step(); j++) {
            exp_m->exp_m_7_->update_network(cloned,selection_pressure,j,i);
          }

          if (pth_target->env_switch_probability_ < 0) {
            if (prev_env_id != pth_target->list_env_id_[i]) {        
              indiv->metadata_->protein_begin();
              for (int j = 0;
                  j < indiv->metadata_->proteins_count(); j++) {
                Protein_7* prot =
                    indiv->metadata_->protein_next();
                if (!prot->signal_) {
                  if (prot->is_init_) {
                    prot->e = prot->initial_e_;
                  }
                }
              }
            }

            prev_env_id = pth_target->list_env_id_[i];
          }


          if (eval->find(i+1) != eval->end()) {
            exp_m->exp_m_7_->evaluate_network(cloned,selection_pressure, pth_target->list_env_id_[i]);
            for (int j = 0; j <  ((Discrete_Double_Fuzzy*) cloned->phenotype)->length() - 1; j++) {
              double xi_start =  (j / (double)
                              ((Discrete_Double_Fuzzy*) cloned->phenotype)->length());
              double xi_end =  ((j+1) / (double)
                              ((Discrete_Double_Fuzzy*) cloned->phenotype)->length());
              mullerfile<<timestep<<","<<repeat<<","<<"0,phenotype"<<","<<xi_start<<","<<((Discrete_Double_Fuzzy*) cloned->phenotype)->points()[j]<<","<<pth_target->list_env_id_[i]<<","<<i<<std::endl;
              mullerfile<<timestep<<","<<repeat<<","<<"0,phenotype"<<","<<xi_end<<","<<((Discrete_Double_Fuzzy*) cloned->phenotype)->points()[j+1]<<","<<pth_target->list_env_id_[i]<<","<<i<<std::endl;
            }

            for (int j=0; j < ((Discrete_Double_Fuzzy*)pth_target->targets_fuzzy_[i])->length() - 1; j++) {
              double xi_start =  (j / (double)
                              ((Discrete_Double_Fuzzy*)pth_target->targets_fuzzy_[i])->length());
              double xi_end =  ((j+1) / (double)
                              ((Discrete_Double_Fuzzy*)pth_target->targets_fuzzy_[i])->length());
              mullerfile<<timestep<<","<<repeat<<","<<"0,env"<<","<<xi_start<<","<<((Discrete_Double_Fuzzy*) pth_target->targets_fuzzy_[i])->points()[j]<<","<<pth_target->list_env_id_[i]<<","<<i<<std::endl;
              mullerfile<<timestep<<","<<repeat<<","<<"0,env"<<","<<xi_end<<","<<((Discrete_Double_Fuzzy*) pth_target->targets_fuzzy_[i])->points()[j+1]<<","<<pth_target->list_env_id_[i]<<","<<i<<std::endl;
            }   
          }

          // Save concentration
          cloned->metadata_->protein_begin();
          int32_t id = 0;
          for (int j = 0;
              j < cloned->metadata_->proteins_count(); j++) {
            Protein_7* prot =
                cloned->metadata_->protein_next();
            if (prot->is_init_ || prot->signal_) {
                concentrationfile<<timestep<<","<<repeat<<","<<i<<","<<id<<","<<prot->e<<","<<(prot->signal_? 1 : 0)<<std::endl;
                id++;
            }
          }
        }
      }

      delete [] cloned->fitness_by_env_id_;
      delete [] cloned->metaerror_by_env_id_;
    }
    
  }
    #else
    Individual_7* cloned = new Individual_7(exp_m,indiv,dna_factory,fuzzy_factory,true);
          
    exp_m->exp_m_7_->evaluate(cloned,w_max,selection_pressure);
    

    for (int i = 0; i <  ((Discrete_Double_Fuzzy*) cloned->phenotype)->length() - 1; i++) {
        double xi_start =  (i / (double)
                        ((Discrete_Double_Fuzzy*) cloned->phenotype)->length());
        double xi_end =  ((i+1) / (double)
                        ((Discrete_Double_Fuzzy*) cloned->phenotype)->length());
        // fprintf(drawingfile, "%lf %lf lineto\n", margin +
        //                                          scale * xi, margin + scale *
        //                                                               ((Discrete_Double_Fuzzy*) cloned->phenotype)->points()[i]);
        mullerfile<<timestep<<","<<"phenotype"<<","<<xi_start<<","<<((Discrete_Double_Fuzzy*) cloned->phenotype)->points()[i]<<","<<i<<","<<env_id<<std::endl;
        mullerfile<<timestep<<","<<"phenotype"<<","<<xi_end<<","<<((Discrete_Double_Fuzzy*) cloned->phenotype)->points()[i+1]<<","<<i<<","<<env_id<<std::endl;
    }
    fprintf( drawingfile, "stroke\n" );


    // ------------------
    //  draw environment
    // ------------------
    fprintf( drawingfile,"[ ] 0 setdash\n" );
    fprintf( drawingfile, "0.001 setlinewidth\n" );
    fprintf( drawingfile, "%lf %lf moveto\n", margin, margin);
    
    for (int i=0; i < ((Discrete_Double_Fuzzy*)target)->length() - 1; i++) {

        double xi_start =  (i / (double)
                        ((Discrete_Double_Fuzzy*)target)->length());
        double xi_end =  ((i+1) / (double)
                        ((Discrete_Double_Fuzzy*)target)->length());
        // fprintf(drawingfile, "%lf %lf lineto\n", margin +
        //                                          scale * xi, margin + scale *
        //                                                               ((Discrete_Double_Fuzzy*) target)->points()[i]);


        mullerfile<<timestep<<","<<"env"<<","<<xi_start<<","<<((Discrete_Double_Fuzzy*) target)->points()[i]<<","<<i<<","<<env_id<<std::endl;
        mullerfile<<timestep<<","<<"env"<<","<<xi_end<<","<<((Discrete_Double_Fuzzy*) target)->points()[i+1]<<","<<i<<","<<env_id<<std::endl;
        //((Discrete_Double_Fuzzy*) target)->print();
    }
    fprintf( drawingfile, "stroke\n" );



    fprintf(drawingfile,"%%%%EOF\n");
    fclose(drawingfile);
    #endif
}


void draw_genetic_unit_with_CDS(Individual_7* indiv, char* directoryName, int64_t timestep) {
  const uint8_t bbsize = 200;  // a4 paper: 595*842
  int32_t gen_length = (indiv->dna_)->length();
  double r = 0.35;
  double scale = 2*M_PI*r/gen_length;

  char filename[128];
  snprintf(filename, 127, "%s/best_genome_with_CDS.eps", directoryName);
  FILE * drawingfile = fopen(filename, "w");


  std::ofstream mullerfile;
  mullerfile.open("stats/dna_cds_visu.csv",std::ofstream::trunc);
  mullerfile<<"generation,type,radius,start_position,end_position,gid"<<std::endl;

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

  mullerfile<<timestep<<","<<"dna"<<","<<r<<","<<"0"<<","<<"360,-1"<<std::endl;

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

  int32_t gid = 0;


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
    indiv->metadata_->protein_begin();
  for (int j = 0; j < indiv->metadata_->proteins_count(); j++) {
    Protein_7* prot = indiv->metadata_->protein_next();
    if (prot != nullptr) {
      if (prot->is_init_) {
        if (prot->leading_lagging == 0) {

    first = prot->protein_start;
    last = prot->protein_end;
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


      mullerfile<<timestep<<","<<"cds"<<","<<r + layer*0.02<<","<<theta_last<<","<<"360"<<","<<gid<<std::endl;
      mullerfile<<timestep<<","<<"cds"<<","<<r + layer*0.02<<","<<"0"<<","<<theta_first<<","<<gid<<std::endl;
      gid++;
    }
    else
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r + layer*0.02, theta_last, theta_first);
      fprintf(drawingfile, "stroke\n");

            mullerfile<<timestep<<","<<"cds"<<","<<r + layer*0.02<<","<<theta_last<<","<<theta_first<<","<<gid<<std::endl;
            gid++;
    }
        }
        }
        }
  }


  // printf("LAGGING\n");
    indiv->metadata_->protein_begin();
  for (int j = 0; j < indiv->metadata_->proteins_count(); j++) {
    Protein_7* prot = indiv->metadata_->protein_next();
    if (prot != nullptr) {
      if (prot->is_init_) {
        if (prot->leading_lagging == 1) {
              first = prot->protein_start;
    last = prot->protein_end;
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


      mullerfile<<timestep<<","<<"cds"<<","<<r - layer*0.02<<","<<theta_last<<","<<"360"<<","<<gid<<std::endl;
      mullerfile<<timestep<<","<<"cds"<<","<<r - layer*0.02<<","<<"0"<<","<<theta_first<<","<<gid<<std::endl;
      gid++;
    }
    else
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r - layer*0.02, theta_first, theta_last);
      fprintf(drawingfile, "stroke\n");

            mullerfile<<timestep<<","<<"cds"<<","<<r - layer*0.02<<","<<theta_first<<","<<theta_last<<","<<gid<<std::endl;
            gid++;

    }
        }
      }
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


void draw_genetic_unit_with_mRNAs(Individual_7* indiv, char* directoryName, int64_t timestep) {
  const uint8_t bbsize = 200;  // a4 paper: 595*842
  int32_t gen_length = (indiv->dna_)->length();
  double r = 0.35;
  double scale = 2*M_PI*r/gen_length;
  int32_t gid = 0;


  char filename[128];
  snprintf(filename, 127, "%s/best_genome_with_mRNAs.eps", directoryName);
  FILE * drawingfile = fopen(filename, "w");


  std::ofstream mullerfile;
  mullerfile.open("stats/dna_mrna_visu.csv",std::ofstream::trunc);
  mullerfile<<"generation,type,radius,start_position,end_position,is_coding,gid"<<std::endl;

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
  
  
  mullerfile<<timestep<<","<<"dna"<<","<<r<<","<<"0"<<","<<"360,1,"<<gid<<std::endl;

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



  indiv->metadata_->rna_begin();
  for (int rna_idx = 0; rna_idx <
                        (int)indiv->metadata_->rna_count(); rna_idx++) {

    Rna_7* rna = indiv->metadata_->rna_next();

    if (rna->is_init_) {
        if (rna->leading_lagging == 0) {
    first = rna->begin;
    last = rna->end;


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
    int32_t is_coding = 0;
    fprintf(drawingfile, "0.018 setlinewidth\n");
    if (rna->is_coding_) { fprintf(drawingfile, "0 0 0 setrgbcolor\n"); is_coding = 1; }
    else fprintf(drawingfile, "0.7 0.7 0.7 setrgbcolor\n");
    layer++; // index starting at 0 but needed to start at 1

    if (theta_last > theta_first)
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r + layer*0.02, theta_last, 360);
      fprintf(drawingfile, "stroke\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r + layer*0.02, 0, theta_first);
      fprintf(drawingfile, "stroke\n");

      mullerfile<<timestep<<","<<"rna"<<","<<r + layer*0.02<<","<<theta_last<<","<<"360"<<","<<is_coding<<","<<gid<<std::endl;
      mullerfile<<timestep<<","<<"rna"<<","<<r + layer*0.02<<","<<"0"<<","<<theta_first<<","<<is_coding<<","<<gid<<std::endl;
      
    }
    else
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r + layer*0.02, theta_last, theta_first);
      fprintf(drawingfile, "stroke\n");

      mullerfile<<timestep<<","<<"rna"<<","<<r + layer*0.02<<","<<theta_last<<","<<theta_first<<","<<is_coding<<","<<gid<<std::endl;
      gid++;
    }
        }
    }
  }



indiv->metadata_->rna_begin();
  for (int rna_idx = 0; rna_idx <
                        (int)indiv->metadata_->rna_count(); rna_idx++) {

    Rna_7* rna = indiv->metadata_->rna_next();

    if (rna->is_init_) {
        if (rna->leading_lagging == 1) {
    first = rna->begin;
    last = rna->end;


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
    int32_t is_coding = 0;
    fprintf(drawingfile, "0.018 setlinewidth\n");
    if (rna->is_coding_) { fprintf(drawingfile, "0 0 0 setrgbcolor\n"); is_coding = 1; }
    else fprintf(drawingfile, "0.7 0.7 0.7 setrgbcolor\n");
    layer++; // index starting at 0 but needed to start at 1

    if (theta_first > theta_last)
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r - layer*0.02, theta_first, 360);
      fprintf(drawingfile, "stroke\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r - layer*0.02, 0, theta_last);
      fprintf(drawingfile, "stroke\n");

      mullerfile<<timestep<<","<<"rna"<<","<<r - layer*0.02<<","<<theta_last<<","<<"360"<<","<<is_coding<<std::endl;
      mullerfile<<timestep<<","<<"rna"<<","<<r - layer*0.02<<","<<"0"<<","<<theta_first<<","<<is_coding<<std::endl;
      gid++;
    }
    else
    {
      fprintf(drawingfile, "newpath\n");
      fprintf(drawingfile, "%lf %lf %lf %d %d arc\n", 0.5, 0.5, r - layer*0.02, theta_first, theta_last);
      fprintf(drawingfile, "stroke\n");

      mullerfile<<timestep<<","<<"rna"<<","<<r - layer*0.02<<","<<theta_first<<","<<theta_last<<","<<is_coding<<","<<gid<<std::endl;
      gid++;
    }
  }
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


/**
 * \brief print help and exist
 */
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
  printf("%s:\n", prog_name);
  printf("\tCreates EPS files with the triangles, the positive and negative\n");
  printf("\tprofiles, the phenotype, the CDS and the mRNAs of the\n");
  printf("\tindividual of interest.\n");
  printf("\n");
  printf("Usage : %s -h or --help\n", prog_name);
  printf("   or : %s -V or --version\n", prog_name);
  printf("   or : %s [-t TIMESTEP] [-I INDEX | -R RANK]\n",
         prog_name);
  printf("\nOptions\n");
  printf("  -h, --help\n\tprint this help, then exit\n");
  printf("  -V, --version\n\tprint version number, then exit\n");
  printf("  -t, --timestep TIMESTEP\n");
  printf("\tspecify timestep of the individual of interest\n");
  printf("  -I, --index INDEX\n");
  printf("\tspecify the index of the individual of interest\n");
  printf("  -R, --rank RANK\n");
  printf("\tspecify the rank of the individual of interest\n");
}


void interpret_cmd_line_options(int argc, char* argv[]) {
  // Define allowed options
  const char * options_list = "hVI:R:t:";
  static struct option long_options_list[] = {
      {"help",      no_argument,       NULL, 'h'},
      {"version",   no_argument,       NULL, 'V' },
      {"index",     required_argument, NULL, 'I'},
      {"rank",      required_argument, NULL, 'R'},
      {"timestep",  required_argument, NULL, 't' },
      { 0, 0, 0, 0 }
  };

  // Get actual values of the command-line options
  int option;
  while ((option = getopt_long(argc, argv, options_list,
                               long_options_list, NULL)) != -1) {
    switch (option) {
      case 'h' : {
        print_help(argv[0]);
        exit(EXIT_SUCCESS);
      }
      case 'V' : {
        Utils::PrintAevolVersion();
        exit(EXIT_SUCCESS);
      }
      case 'I' : {
        indiv_index  = atol(optarg);
        break;
      }
      case 'R' : {
        indiv_rank  = atol(optarg);
        break;
      }
      case 't' : {
        if (strcmp(optarg, "") == 0) {
          printf("%s: error: Option -t or --timestep: missing argument.\n",
                 argv[0]);
          exit(EXIT_FAILURE);
        }
        timestep = atol(optarg);
        break;
      }
    }
  }

  // If timestep wasn't provided, use default
  if (timestep < 0) {
    timestep = OutputManager::last_gener();
  }

  // If neither the rank nor the index were provided, the individual of interest
  // will be the best individual at the provided timestep
}
