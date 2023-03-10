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
#include "Stats.h"
#include "ExpManager_7.h"
#include <string>

#include <err.h>
#include <errno.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>

#include "StatRecord.h"
#include "ExpManager.h"
#include "ExpSetup.h"
#include "Individual.h"
#include "GeneticUnit.h"
#ifdef __REGUL
  #include "raevol/Protein_R.h"
#endif

using std::string;

namespace aevol {

//##############################################################################
//                                                                             #
//                                Class Stats                               #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
/**
 * This is a temporary patch for experiment propagation, it shall become
 * obsolete or need to be adapted when in/out dirs are managed properly
 */
    Stats::Stats(const string prefix, const string postfix /*= ""*/,
                 bool best_indiv_only /*= false*/) {
      init_data();
      set_file_names(prefix, postfix, best_indiv_only);
    }

/**
 * Create a NEW stat manager
 */
    Stats::Stats(ExpManager * exp_m,
                 bool best_indiv_only /*= false*/,
                 const string prefix /*= "stat"*/,
                 const string postfix /*= ""*/,
                 bool with_plasmids /*= false*/,
                 bool compute_phen_contrib_by_GU /*= false*/) {
      exp_m_ = exp_m;
      init_data();
      set_file_names(prefix, postfix, best_indiv_only, with_plasmids,
                     compute_phen_contrib_by_GU);
      open_files();
      write_headers();
    }

/**
 * Create a stat manager to append to existing stats
 */
    Stats::Stats(ExpManager * exp_m,
                 int64_t time,
                 bool best_indiv_only/* = false */,
                 const string prefix /*= "stat"*/,
                 const string postfix /*= ""*/,
                 bool addition_old_stats /* = true */,
                 bool delete_old_stats /* = true */) {
      exp_m_ = exp_m;
      init_data();
      set_file_names(prefix, postfix, best_indiv_only);

      if (addition_old_stats) {
        CreateTmpFiles(time);
        PromoteTmpFiles();
      }
      else { // ancstat case
        open_files();
        write_headers(true);
      }

      // Flush the new stat files
      flush();
    }

// =================================================================
//                             Destructors
// =================================================================
Stats::~Stats() {
  for (int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++) {
    for (int8_t best_or_glob = 0 ;
         best_or_glob < NB_BEST_OR_GLOB ;
         best_or_glob++) {
      for (int8_t stat_type = 0 ; stat_type < NB_STATS_TYPES ; stat_type++) {
        if(stat_files_names_[chrom_or_GU][best_or_glob][stat_type] != nullptr) {
          if (stat_files_[chrom_or_GU][best_or_glob][stat_type] != nullptr) {
            fclose(stat_files_[chrom_or_GU][best_or_glob][stat_type]);
            stat_files_[chrom_or_GU][best_or_glob][stat_type] = nullptr;
          }
          delete [] stat_files_names_[chrom_or_GU][best_or_glob][stat_type];
          stat_files_names_[chrom_or_GU][best_or_glob][stat_type] = nullptr;
        }
      }

      delete [] stat_files_[chrom_or_GU][best_or_glob];
      stat_files_[chrom_or_GU][best_or_glob] = nullptr;

      delete [] stat_files_names_[chrom_or_GU][best_or_glob];
      stat_files_names_[chrom_or_GU][best_or_glob] = nullptr;
    }

    delete [] stat_files_[chrom_or_GU];
    stat_files_[chrom_or_GU] = nullptr;

    delete [] stat_files_names_[chrom_or_GU];
    stat_files_names_[chrom_or_GU] = nullptr;
  }

  delete [] stat_files_;
  stat_files_ = nullptr;

  delete [] stat_files_names_;
  stat_files_names_ = nullptr;

  if (!ExpManager_7::standalone_simd) {
    std::map<long long int, Individual *> unique_individual;

    for (auto g_indivs = indivs_.begin(); g_indivs != indivs_.end(); ++g_indivs) {
      for (auto &indiv: g_indivs->second) {
        unique_individual[indiv->long_id()] = indiv;
        indiv->number_of_clones_--;
      }
    }

    for (auto indiv_it = unique_individual.begin(); indiv_it != unique_individual.end(); ++indiv_it) {
      delete indiv_it->second;
    }
  }
}


// =================================================================
//                            Public Methods
// =================================================================

inline double sqr(double x)
{
  return x*x;
}

inline double rsqr(double x)
{
  return (x < 0.00000001) ? 0. : sqrt(x);
}


void Stats::write_headers(bool ancstats_stats /* = false */)
{
  // Column key in the stat files
  int8_t key;

  // --------------------------------------
  //  Write headers in FITNESS_STATS files
  // --------------------------------------
  for (int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++)
  {
    if(ancstats_stats)
    {
      write_header(stat_files_[chrom_or_GU][BEST][FITNESS_STATS], "----------------------------------------");
      write_header(stat_files_[chrom_or_GU][BEST][FITNESS_STATS], " Lineage individuals fitness statistics ");
      write_header(stat_files_[chrom_or_GU][BEST][FITNESS_STATS], "----------------------------------------");
      write_header(stat_files_[chrom_or_GU][BEST][FITNESS_STATS], "");
    }
    else
    {
      if (stat_files_names_[chrom_or_GU][BEST][FITNESS_STATS] != nullptr)
      {
        write_header(stat_files_[chrom_or_GU][BEST][FITNESS_STATS], "---------------------------------------");
        write_header(stat_files_[chrom_or_GU][BEST][FITNESS_STATS], " Fittest individual fitness statistics ");
        write_header(stat_files_[chrom_or_GU][BEST][FITNESS_STATS], "---------------------------------------");
        write_header(stat_files_[chrom_or_GU][BEST][FITNESS_STATS], "");
      }
      if (stat_files_names_[chrom_or_GU][GLOB][FITNESS_STATS] != nullptr)
      {
        write_header(stat_files_[chrom_or_GU][GLOB][FITNESS_STATS], "------------------------------------------------");
        write_header(stat_files_[chrom_or_GU][GLOB][FITNESS_STATS], " Average fitness statistics over the population ");
        write_header(stat_files_[chrom_or_GU][GLOB][FITNESS_STATS], "------------------------------------------------");
        write_header(stat_files_[chrom_or_GU][GLOB][FITNESS_STATS], "");
      }

      if (stat_files_names_[chrom_or_GU][SDEV][FITNESS_STATS] != nullptr)
      {
        write_header(stat_files_[chrom_or_GU][SDEV][FITNESS_STATS], "------------------------------------------------------------");
        write_header(stat_files_[chrom_or_GU][SDEV][FITNESS_STATS], " Standard deviation, fitness statistics over the population ");
        write_header(stat_files_[chrom_or_GU][SDEV][FITNESS_STATS], "------------------------------------------------------------");
        write_header(stat_files_[chrom_or_GU][SDEV][FITNESS_STATS], "");
      }

      if (stat_files_names_[chrom_or_GU][SKEW][FITNESS_STATS] != nullptr)
      {
        write_header(stat_files_[chrom_or_GU][SKEW][FITNESS_STATS], "--------------------------------------------------");
        write_header(stat_files_[chrom_or_GU][SKEW][FITNESS_STATS], " Skewness statistics, fitness over the population ");
        write_header(stat_files_[chrom_or_GU][SKEW][FITNESS_STATS], "--------------------------------------------------");
        write_header(stat_files_[chrom_or_GU][SKEW][FITNESS_STATS], "");
      }
    }


    for (int8_t best_or_glob = 0 ; best_or_glob < NB_BEST_OR_GLOB ; best_or_glob++)
    {
      if (stat_files_names_[chrom_or_GU][best_or_glob][FITNESS_STATS] != nullptr)
      {
        assert(stat_files_[chrom_or_GU][best_or_glob][FITNESS_STATS] != nullptr);
        key = 1;

        write_header(stat_files_[chrom_or_GU][best_or_glob][FITNESS_STATS], "Generation", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][FITNESS_STATS], "Population size", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][FITNESS_STATS], "Fitness", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][FITNESS_STATS], "Genome size (amount of DNA)", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][FITNESS_STATS], "Metabolic error", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][FITNESS_STATS], "Parent's metabolic error", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][FITNESS_STATS], "Metabolic fitness", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][FITNESS_STATS], "Secretion error", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][FITNESS_STATS], "Parent's secretion error", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][FITNESS_STATS], "Secretion fitness", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][FITNESS_STATS], "Amount of compound present in the grid-cell", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][FITNESS_STATS], "Int probe", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][FITNESS_STATS], "Double probe", key++);

        #ifdef __REGUL
          write_header(stat_files_[chrom_or_GU][best_or_glob][FITNESS_STATS], "Number of links in the regulation graph", key++);
          write_header(stat_files_[chrom_or_GU][best_or_glob][FITNESS_STATS], "Number of positive links in the regulation graph", key++);
          write_header(stat_files_[chrom_or_GU][best_or_glob][FITNESS_STATS], "Number of negative links in the regulation graph", key++);
          write_header(stat_files_[chrom_or_GU][best_or_glob][FITNESS_STATS], "Average value of links in the regulation graph", key++);
          write_header(stat_files_[chrom_or_GU][best_or_glob][FITNESS_STATS], "Average value of positive links in the regulation graph", key++);
          write_header(stat_files_[chrom_or_GU][best_or_glob][FITNESS_STATS], "Average value of negative links in the regulation graph", key++);
        #endif

        write_header(stat_files_[chrom_or_GU][best_or_glob][FITNESS_STATS], "");
      }
    }
  }

  // ---------------------------------------
  //  Write headers in MUTATION_STATS files
  // ---------------------------------------
  for (int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++)
  {
    if(ancstats_stats)
    {
      write_header(stat_files_[chrom_or_GU][BEST][MUTATION_STATS], "-----------------------------------------");
      write_header(stat_files_[chrom_or_GU][BEST][MUTATION_STATS], " Lineage individuals mutation statistics ");
      write_header(stat_files_[chrom_or_GU][BEST][MUTATION_STATS], "-----------------------------------------");
      write_header(stat_files_[chrom_or_GU][BEST][MUTATION_STATS], "");
    }
    else
    {
      if (stat_files_names_[chrom_or_GU][BEST][MUTATION_STATS] != nullptr)
      {
        write_header(stat_files_[chrom_or_GU][BEST][MUTATION_STATS], "----------------------------------------");
        write_header(stat_files_[chrom_or_GU][BEST][MUTATION_STATS], " Fittest individual mutation statistics ");
        write_header(stat_files_[chrom_or_GU][BEST][MUTATION_STATS], "----------------------------------------");
        write_header(stat_files_[chrom_or_GU][BEST][MUTATION_STATS], "");
      }
      if (stat_files_names_[chrom_or_GU][GLOB][MUTATION_STATS] != nullptr)
      {
        write_header(stat_files_[chrom_or_GU][GLOB][MUTATION_STATS], "-------------------------------------------------");
        write_header(stat_files_[chrom_or_GU][GLOB][MUTATION_STATS], " Average mutation statistics over the population ");
        write_header(stat_files_[chrom_or_GU][GLOB][MUTATION_STATS], "-------------------------------------------------");
        write_header(stat_files_[chrom_or_GU][GLOB][MUTATION_STATS], "");
      }
    }

    for (int8_t best_or_glob = 0 ; best_or_glob < NB_BEST_OR_GLOB ; best_or_glob++)
    {
      if (stat_files_names_[chrom_or_GU][best_or_glob][MUTATION_STATS] != nullptr)
      {
        assert(stat_files_[chrom_or_GU][best_or_glob][MUTATION_STATS] != nullptr);
        key = 1;

        write_header(stat_files_[chrom_or_GU][best_or_glob][MUTATION_STATS], "Generation", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][MUTATION_STATS], "Number of local mutations undergone", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][MUTATION_STATS], "Number of chromosomic rearrangements undergone", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][MUTATION_STATS], "Number of switch undergone", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][MUTATION_STATS], "Number of indels undergone", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][MUTATION_STATS], "Number of duplications undergone", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][MUTATION_STATS], "Number of deletions undergone", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][MUTATION_STATS], "Number of translocations undergone", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][MUTATION_STATS], "Number of inversions undergone", key++);

        write_header(stat_files_[chrom_or_GU][best_or_glob][MUTATION_STATS], "");
      }
    }
  }

  // ---------------------------------------
  //  Write headers in GENES_STATS files
  // ---------------------------------------
  for (int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++)
  {
    if(ancstats_stats)
    {
      write_header(stat_files_[chrom_or_GU][BEST][GENES_STATS], "-------------------------------------");
      write_header(stat_files_[chrom_or_GU][BEST][GENES_STATS], " Lineage individuals gene statistics ");
      write_header(stat_files_[chrom_or_GU][BEST][GENES_STATS], "-------------------------------------");
      write_header(stat_files_[chrom_or_GU][BEST][GENES_STATS], "");
    }
    else
    {
      if (stat_files_names_[chrom_or_GU][BEST][GENES_STATS] != nullptr)
      {
        write_header(stat_files_[chrom_or_GU][BEST][GENES_STATS], "------------------------------------");
        write_header(stat_files_[chrom_or_GU][BEST][GENES_STATS], " Fittest individual gene statistics ");
        write_header(stat_files_[chrom_or_GU][BEST][GENES_STATS], "------------------------------------");
        write_header(stat_files_[chrom_or_GU][BEST][GENES_STATS], "");
      }
      if (stat_files_names_[chrom_or_GU][GLOB][GENES_STATS] != nullptr)
      {
        write_header(stat_files_[chrom_or_GU][GLOB][GENES_STATS], "---------------------------------------------");
        write_header(stat_files_[chrom_or_GU][GLOB][GENES_STATS], " Average gene statistics over the population ");
        write_header(stat_files_[chrom_or_GU][GLOB][GENES_STATS], "---------------------------------------------");
        write_header(stat_files_[chrom_or_GU][GLOB][GENES_STATS], "");
      }
    }

    for (int8_t best_or_glob = 0 ; best_or_glob < NB_BEST_OR_GLOB ; best_or_glob++)
    {
      if (stat_files_names_[chrom_or_GU][best_or_glob][GENES_STATS] != nullptr)
      {
        assert(stat_files_[chrom_or_GU][best_or_glob][GENES_STATS] != nullptr);
        key = 1;

        write_header(stat_files_[chrom_or_GU][best_or_glob][GENES_STATS], "Generation", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][GENES_STATS], "Number of coding RNAs (at least one gene on RNA)", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][GENES_STATS], "Number of non-coding RNAs", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][GENES_STATS], "Average size of coding RNAs (at least one gene on RNA)", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][GENES_STATS], "Average size of non-coding RNAs", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][GENES_STATS], "Number of functional genes", key++);
        // Non functional genes are those with width_ == 0 or height_ == 0 or those that lack one kind of codons (M, W or H)
        write_header(stat_files_[chrom_or_GU][best_or_glob][GENES_STATS], "Nb of non functional genes", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][GENES_STATS], "Average size of functional genes", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][GENES_STATS], "Average size of non functional genes (WARNING : bias towards 0)", key++);
#ifdef __REGUL
        write_header(stat_files_[chrom_or_GU][best_or_glob][GENES_STATS], "Number of Transcription Factors", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][GENES_STATS], "Number of pure Transcription Factors", key++);
#endif
        write_header(stat_files_[chrom_or_GU][best_or_glob][GENES_STATS], "");
      }
    }
  }

  // ---------------------------------------
  //  Write headers in BP_STATS files
  // ---------------------------------------
  for (int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++)
  {
    if(ancstats_stats)
    {
      write_header(stat_files_[chrom_or_GU][BEST][BP_STATS], "-------------------------------------------");
      write_header(stat_files_[chrom_or_GU][BEST][BP_STATS], " Lineage individuals non-coding statistics ");
      write_header(stat_files_[chrom_or_GU][BEST][BP_STATS], "-------------------------------------------");
      write_header(stat_files_[chrom_or_GU][BEST][BP_STATS], "");
    }
    else
    {
      if (stat_files_names_[chrom_or_GU][BEST][BP_STATS] != nullptr)
      {
        write_header(stat_files_[chrom_or_GU][BEST][BP_STATS], "------------------------------------------");
        write_header(stat_files_[chrom_or_GU][BEST][BP_STATS], " Fittest individual non-coding statistics ");
        write_header(stat_files_[chrom_or_GU][BEST][BP_STATS], "------------------------------------------");
        write_header(stat_files_[chrom_or_GU][BEST][BP_STATS], "");
      }
      if (stat_files_names_[chrom_or_GU][GLOB][BP_STATS] != nullptr)
      {
        write_header(stat_files_[chrom_or_GU][GLOB][BP_STATS], "---------------------------------------------------");
        write_header(stat_files_[chrom_or_GU][GLOB][BP_STATS], " Average non-coding statistics over the population ");
        write_header(stat_files_[chrom_or_GU][GLOB][BP_STATS], "---------------------------------------------------");
        write_header(stat_files_[chrom_or_GU][GLOB][BP_STATS], "");
        write_header(stat_files_[chrom_or_GU][GLOB][BP_STATS], " This data is not available");
        write_header(stat_files_[chrom_or_GU][GLOB][BP_STATS], " Computing bp stats for all individuals is extremely costly computationaly");
        write_header(stat_files_[chrom_or_GU][GLOB][BP_STATS], "");

        // Mark file as "not to be written into" and close it
        delete [] stat_files_names_[chrom_or_GU][GLOB][BP_STATS];
        stat_files_names_[chrom_or_GU][GLOB][BP_STATS] = nullptr;
        fclose(stat_files_[chrom_or_GU][GLOB][BP_STATS]);
        stat_files_[chrom_or_GU][GLOB][BP_STATS] = nullptr;
      }
    }

    for (int8_t best_or_glob = 0 ; best_or_glob < NB_BEST_OR_GLOB ; best_or_glob++)
    {
      if (stat_files_names_[chrom_or_GU][best_or_glob][BP_STATS] != nullptr)
      {
        assert(stat_files_[chrom_or_GU][best_or_glob][BP_STATS] != nullptr);
        key = 1;

        write_header(stat_files_[chrom_or_GU][best_or_glob][BP_STATS], "Generation", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][BP_STATS], "Number of bp not included in any CDS", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][BP_STATS], "Number of bp not included in any functional CDS", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][BP_STATS], "Number of bp not included in any non functional CDS", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][BP_STATS], "Number of bp not included in any RNA", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][BP_STATS], "Number of bp not included in any coding RNA", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][BP_STATS], "Number of bp not included in any non coding RNA", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][BP_STATS], "Number of non essential bp", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][BP_STATS], "Number of non essential bp including non fonctional genes", key++);

        write_header(stat_files_[chrom_or_GU][best_or_glob][BP_STATS], "");
        write_header(stat_files_[chrom_or_GU][best_or_glob][BP_STATS], "NOTE: a bp is considered \"essential\" when it is part of any [functional] CDS");
        write_header(stat_files_[chrom_or_GU][best_or_glob][BP_STATS], "  or any promoter or terminator corresponding to an RNA transcribing a [functional] CDS.");
        write_header(stat_files_[chrom_or_GU][best_or_glob][BP_STATS], "");
      }
    }
  }

  // ---------------------------------------
  //  Write headers in REAR_STATS files
  // ---------------------------------------
  for (int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++)
  {
    if(ancstats_stats)
    {
      write_header(stat_files_[chrom_or_GU][BEST][REAR_STATS], "----------------------------------------------");
      write_header(stat_files_[chrom_or_GU][BEST][REAR_STATS], " Lineage individuals rearrangement statistics ");
      write_header(stat_files_[chrom_or_GU][BEST][REAR_STATS], "----------------------------------------------");
      write_header(stat_files_[chrom_or_GU][BEST][REAR_STATS], "");
    }
    else
    {
      if (stat_files_names_[chrom_or_GU][BEST][REAR_STATS] != nullptr)
      {
        write_header(stat_files_[chrom_or_GU][BEST][REAR_STATS], "---------------------------------------------");
        write_header(stat_files_[chrom_or_GU][BEST][REAR_STATS], " Fittest individual rearrangement statistics ");
        write_header(stat_files_[chrom_or_GU][BEST][REAR_STATS], "---------------------------------------------");
        write_header(stat_files_[chrom_or_GU][BEST][REAR_STATS], "");
      }
      if (stat_files_names_[chrom_or_GU][GLOB][REAR_STATS] != nullptr)
      {
        write_header(stat_files_[chrom_or_GU][GLOB][REAR_STATS], "------------------------------------------------------");
        write_header(stat_files_[chrom_or_GU][GLOB][REAR_STATS], " Average rearrangement statistics over the population ");
        write_header(stat_files_[chrom_or_GU][GLOB][REAR_STATS], "------------------------------------------------------");
        write_header(stat_files_[chrom_or_GU][GLOB][REAR_STATS], "");
      }
    }

    for (int8_t best_or_glob = 0 ; best_or_glob < NB_BEST_OR_GLOB ; best_or_glob++)
    {
      if (stat_files_names_[chrom_or_GU][best_or_glob][REAR_STATS] != nullptr)
      {
        assert(stat_files_[chrom_or_GU][best_or_glob][REAR_STATS] != nullptr);
        key = 1;

        write_header(stat_files_[chrom_or_GU][best_or_glob][REAR_STATS], "Generation", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][REAR_STATS], "Actual duplication rate", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][REAR_STATS], "Actual deletion rate", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][REAR_STATS], "Actual translocation rate", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][REAR_STATS], "Actual inversion rate", key++);
        write_header(stat_files_[chrom_or_GU][best_or_glob][REAR_STATS], "Average alignment score (needed score)", key++);

        write_header(stat_files_[chrom_or_GU][best_or_glob][REAR_STATS], "");
      }
    }
  }

  flush();
}

    void Stats::write_current_generation_statistics(int64_t gen) {
      // debug std::cout << "writing stats for gen : " << gen << '\n';
      StatRecord** stat_records;

#ifdef __OPENMP_TASK
      bool stat_depend[4][NB_CHROM_OR_GU];
  bool stat_file_depend[4][NB_CHROM_OR_GU][NB_STATS_TYPES];
  ExpManager* exp_m = exp_m_;
#endif

      Individual* best_indiv = b_indiv(gen);
      ReplicationReport* best_replic_report = exp_m_->tree() ?
                                              exp_m_->tree()->report_by_index(gen, best_indiv->id()) :
                                              nullptr;
      std::list<std::pair<Individual*, ReplicationReport*>> annotated_indivs =
              indivs_annotated(gen);

      for (int8_t chrom_or_GU = 0; chrom_or_GU < NB_CHROM_OR_GU; chrom_or_GU++) {
        if ((not exp_m_->output_m()->compute_phen_contrib_by_GU()) &&
            chrom_or_GU > ALL_GU)
          continue;

        stat_records = new StatRecord* [NB_BEST_OR_GLOB];
#ifdef __OPENMP_TASK
        #pragma omp taskgroup
{
#pragma omp task shared(exp_m) depend(out: stat_depend[BEST][chrom_or_GU])
#endif
        stat_records[BEST] = new StatRecord(exp_m_->exp_s(),
                                            best_indiv,
                                            best_replic_report,
                                            (chrom_or_gen_unit) chrom_or_GU);

#ifdef __OPENMP_TASK
#pragma omp task shared(exp_m) depend(out: stat_depend[GLOB][chrom_or_GU])
#endif
        stat_records[GLOB] = new StatRecord(exp_m_->exp_s(),
                                            annotated_indivs,
                                            (chrom_or_gen_unit) chrom_or_GU);

#ifdef __OPENMP_TASK
        #pragma omp task shared(exp_m) depend(inout: stat_depend[GLOB][chrom_or_GU])\
                               depend(out: stat_depend[SDEV][chrom_or_GU])
#endif
        stat_records[SDEV] = new StatRecord(exp_m_->exp_s(),
                                            annotated_indivs,
                                            stat_records[GLOB],
                                            (chrom_or_gen_unit) chrom_or_GU);

#ifdef __OPENMP_TASK
        #pragma omp task shared(exp_m) depend(inout: stat_depend[GLOB][chrom_or_GU])\
                  depend(inout: stat_depend[SDEV][chrom_or_GU])\
                  depend(out: stat_depend[SKEW][chrom_or_GU])
#endif
        stat_records[SKEW] = new StatRecord(exp_m_->exp_s(),
                                            annotated_indivs,
                                            stat_records[GLOB],
                                            stat_records[SDEV],
                                            (chrom_or_gen_unit) chrom_or_GU);
//#pragma omp taskwait
        for (int8_t best_or_glob = 0; best_or_glob < NB_BEST_OR_GLOB; best_or_glob++) {
          for (int8_t stat_type = 0; stat_type < NB_STATS_TYPES; stat_type++) {
            if (stat_files_names_[chrom_or_GU][best_or_glob][stat_type] != NULL) {
#ifdef __OPENMP_TASK
              #pragma omp task \
                  depend(inout: stat_depend[best_or_glob][chrom_or_GU])\
                  depend(inout: stat_file_depend[best_or_glob][chrom_or_GU][stat_type])
#endif
              stat_records[best_or_glob]->write_to_file(gen,
                                                        stat_files_[chrom_or_GU][best_or_glob][stat_type],(stats_type) stat_type);
            }
          }

#ifdef __OPENMP_TASK
#pragma omp task depend(in: stat_depend[best_or_glob][chrom_or_GU])
#endif
          delete stat_records[best_or_glob];
        }
#ifdef __OPENMP_TASK
        }
#endif
        delete[] stat_records;
      }
    }

    void Stats::write_statistics_of_this_indiv(int64_t time, Individual * indiv,
                                               ReplicationReport* replic_report)
    {
      StatRecord * stat_record;

      for (int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++)
      {
        stat_record = new StatRecord(exp_m_->exp_s(), indiv, replic_report,
                                     (chrom_or_gen_unit) chrom_or_GU, true);

        for (int8_t stat_type = 0 ; stat_type < NB_STATS_TYPES ; stat_type++)
        {
          if (stat_files_names_[chrom_or_GU][BEST][stat_type] != nullptr)
          {
            assert(stat_files_[chrom_or_GU][BEST][stat_type] != nullptr);
            stat_record->write_to_file(time, stat_files_[chrom_or_GU][BEST][stat_type],
                                       (stats_type) stat_type);
          }
        }
        delete stat_record;
      }
    }



void Stats::flush() {
  for (int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++) {
    for (int8_t best_or_glob = 0 ;
         best_or_glob < NB_BEST_OR_GLOB ;
         best_or_glob++) {
      for (int8_t stat_type = 0 ; stat_type < NB_STATS_TYPES ; stat_type++) {
        if (stat_files_names_[chrom_or_GU][best_or_glob][stat_type] != nullptr) {
          assert(stat_files_[chrom_or_GU][best_or_glob][stat_type] != nullptr);
          fflush(stat_files_[chrom_or_GU][best_or_glob][stat_type]);
        }
      }
    }
  }
}

void Stats::add_indivs(int64_t gen, const std::list<Individual*> indivs) {
  indivs_[gen] = indivs;
}


    void Stats::add_indivs(int64_t gen, Individual_7** indivs) {
      std::list<Individual*> gen_indivs;
      //printf("At gen %ld AddIndivs %d\n",gen,exp_m_->nb_indivs());
      for (int i = 0; i < exp_m_->nb_indivs(); i++) {
          int x = indivs[i]->indiv_id / exp_m_->world()->height();
          int y = indivs[i]->indiv_id % exp_m_->world()->height();



        Individual * indiv = new Individual(exp_m_,
                                            exp_m_->world()->grid(x,y)->mut_prng(),
                                            exp_m_->world()->grid(x,y)->stoch_prng(),
                                            exp_m_->exp_s()->mut_params(),
                                            indivs[i]->w_max_,
                                            exp_m_->exp_s()->min_genome_length(),
                                            exp_m_->exp_s()->max_genome_length(),
                                            false,
                                            indivs[i]->indiv_id,
                                            "",
                                            0);
        gen_indivs.push_back(indiv);

      }
      indivs_[gen] = gen_indivs;
    }

// =================================================================
//                           Protected Methods
// =================================================================
/**
 * Allocate memory and initialize file handlers and file names to NULL
 */
void Stats::init_data() {
  stat_files_       = new FILE***[NB_CHROM_OR_GU];
  stat_files_names_ = new char***[NB_CHROM_OR_GU];

  for (int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++) {
    stat_files_[chrom_or_GU]        = new FILE**[NB_BEST_OR_GLOB];
    stat_files_names_[chrom_or_GU]  = new char**[NB_BEST_OR_GLOB];

    for (int8_t best_or_glob = 0 ;
         best_or_glob < NB_BEST_OR_GLOB ;
         best_or_glob++) {
      stat_files_[chrom_or_GU][best_or_glob] = new FILE*[NB_STATS_TYPES];
      stat_files_names_[chrom_or_GU][best_or_glob] = new char*[NB_STATS_TYPES];

      for (int8_t stat_type = 0 ; stat_type < NB_STATS_TYPES ; stat_type++) {
        stat_files_[chrom_or_GU][best_or_glob][stat_type]       = nullptr;
        stat_files_names_[chrom_or_GU][best_or_glob][stat_type] = nullptr;
      }
    }
  }
}

/**
 * Construct file names
 */
// NB: Here is where we "choose" which files we will write.
//     Files that are not wanted must be left with a nullptr name
//     (don't new char[] them)
//     There is an exception though: for the non-coding file for the population,
//     we will give it a name temporarily so that we can write the warning
//     headers. Once this is done, the name will be deleted to mark the file as
//     "not to be written into"
//
    void Stats::set_file_names(const string prefix,
                               const string postfix,
                               bool best_indiv_only,
                               bool with_plasmids /*= false*/,
                               bool compute_phen_contrib_by_GU /*= false*/) {
        // 1) Create stats directory
        int status;
        status = mkdir(STATS_DIR, 0755);
        if ((status == -1) && (errno != EEXIST)) {
            err(EXIT_FAILURE, STATS_DIR);
        }

        const char* chrom_or_gu_name[NB_CHROM_OR_GU] =
                {"", "_chromosome", "_plasmids"};
        const char* best_or_glob_name[NB_BEST_OR_GLOB] =
                {"_best", "_glob", "_sdev", "_skew"};
        const char* stat_type_name[NB_STATS_TYPES] =
                {"_fitness" ,"_mutation", "_genes", "_bp", "_rear"};

        for (int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++) {
            // If plasmids are not allowed, don't issue "chromosome" and
            // "plasmids" files
            if (not with_plasmids && chrom_or_GU > 0)
                continue;

            // Idem if COMPUTE_PHEN_CONTRIB_BY_GU not set
            if ((not compute_phen_contrib_by_GU && chrom_or_GU > ALL_GU))
                continue;


            for (int8_t best_or_glob = 0 ;
                 best_or_glob < NB_BEST_OR_GLOB ;
                 best_or_glob++) {
                if (best_indiv_only && best_or_glob != BEST) continue;

                for (int8_t stat_type = 0 ; stat_type < NB_STATS_TYPES ; stat_type++) {
                    // // We don't want REAR_STATS when rearrangements are done without
                    // alignments
                    // if (stat_type == REAR_STATS && ! exp_m_->with_alignments())
                    //   continue;

                    // For now, we only want sdev and skew for fitness data
                    if (best_or_glob > GLOB && stat_type > FITNESS_STATS) continue;
                    if ((chrom_or_GU != ALL_GU || best_or_glob != GLOB) &&
                        stat_type > REAR_STATS) continue;

                    stat_files_names_[chrom_or_GU][best_or_glob][stat_type] = new char[255];

                    // Construct the correct name
                    if (best_indiv_only) {
                        sprintf(stat_files_names_[chrom_or_GU][best_or_glob][stat_type],
                                STATS_DIR"/%s%s%s%s.out",
                                prefix.c_str(),
                                stat_type_name[stat_type],
                                chrom_or_gu_name[chrom_or_GU],
                                postfix.c_str());
                    }
                    else
                    {
                        sprintf(stat_files_names_[chrom_or_GU][best_or_glob][stat_type],
                                STATS_DIR"/%s%s%s%s%s.out",
                                prefix.c_str(),
                                stat_type_name[stat_type],
                                chrom_or_gu_name[chrom_or_GU],
                                best_or_glob_name[best_or_glob],
                                postfix.c_str());
                    }

                }
            }

        }
    }

/**
 * Open files that have a non nullptr name
 */
void Stats::open_files()
{
  for (int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++) {
    for (int8_t best_or_glob = 0 ;
         best_or_glob < NB_BEST_OR_GLOB ;
         best_or_glob++) {
      for (int8_t stat_type = 0 ; stat_type < NB_STATS_TYPES ; stat_type++) {
        if (stat_files_names_[chrom_or_GU][best_or_glob][stat_type] != nullptr) {
          stat_files_[chrom_or_GU][best_or_glob][stat_type] =
              fopen(stat_files_names_[chrom_or_GU][best_or_glob][stat_type],
                    "w");
        }
      }
    }
  }
}

/**
 * Create partial copies (up to a given timestamp) of all stat files
 */
void Stats::CreateTmpFiles(int64_t time) {
  char* old_file_name;  // Syntaxic sugar for stat_files_names_[][][]
  FILE* old_file;
  char* new_file_name = new char[500];
  FILE* new_file;
  char  line[500];

  for (int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++) {
    for (int8_t best_or_glob = 0 ;
         best_or_glob < NB_BEST_OR_GLOB ;
         best_or_glob++) {
      for (int8_t stat_type = 0 ; stat_type < NB_STATS_TYPES ; stat_type++) {
        old_file_name = stat_files_names_[chrom_or_GU][best_or_glob][stat_type];
        if (old_file_name != nullptr) {
          sprintf(new_file_name, "%s.tmp", old_file_name);

          old_file = fopen(old_file_name, "r");
          new_file = fopen(new_file_name, "w");

          if (old_file == nullptr)
                  continue;

          // Copy file header
          if (fgets(line, 500, old_file) == nullptr) {
            // TODO check for error
          }

          while (!feof(old_file) && line[0] == '#') {
            fputs(line, new_file);
            if (fgets(line, 500, old_file) == nullptr) {
              // TODO check for error
              break;
            }
          }

          while ((int64_t)atol(line) <= time && !feof(old_file)) {
            fputs(line, new_file);
            if (fgets(line, 500, old_file) == nullptr) {
              // TODO check for error
            }
          }

          fclose(old_file);
          fclose(new_file);
        }
      }
    }
  }

  delete [] new_file_name;
}

/**
 * Replace all the stat files by their tmp counterpart
 */
void Stats::PromoteTmpFiles() {
  char* cur_file_name;  // Syntaxic sugar for stat_files_names_[][][]
  char* tmp_file_name = new char[100];

  for (int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++) {
    for (int8_t best_or_glob = 0 ;
         best_or_glob < NB_BEST_OR_GLOB ;
         best_or_glob++) {
      for (int8_t stat_type = 0 ; stat_type < NB_STATS_TYPES ; stat_type++) {
        cur_file_name = stat_files_names_[chrom_or_GU][best_or_glob][stat_type];
        if (cur_file_name != NULL) {
          sprintf(tmp_file_name, "%s.tmp", cur_file_name);

          remove(cur_file_name);
          int renameOK = rename(tmp_file_name, cur_file_name);
          if (renameOK != 0)
            Utils::ExitWithUsrMsg(string("could not rename file ") +
                                      tmp_file_name + " into " +
                                      cur_file_name);

          // Reopen file
          if (stat_files_[chrom_or_GU][best_or_glob][stat_type] != NULL)
            fclose(stat_files_[chrom_or_GU][best_or_glob][stat_type]);
          stat_files_[chrom_or_GU][best_or_glob][stat_type] =
              fopen(cur_file_name, "a");
        }
      }
    }
  }

  delete [] tmp_file_name;
}

void Stats::MoveTmpFiles(const std::string& destdir) {
  char* cur_file_name;  // Syntaxic sugar for stat_files_names_[][][]
  string tmp_file_name;
  string dest_file_name;

  for (int8_t chrom_or_GU = 0 ; chrom_or_GU < NB_CHROM_OR_GU ; chrom_or_GU++) {
    for (int8_t best_or_glob = 0 ;
         best_or_glob < NB_BEST_OR_GLOB ;
         best_or_glob++) {
      for (int8_t stat_type = 0 ; stat_type < NB_STATS_TYPES ; stat_type++) {
        cur_file_name = stat_files_names_[chrom_or_GU][best_or_glob][stat_type];
        if (cur_file_name != NULL) {
          tmp_file_name = string(cur_file_name) + ".tmp";
          dest_file_name = destdir + "/" + cur_file_name;
          int renameOK = rename(tmp_file_name.c_str(),
                                dest_file_name.c_str());
          if (renameOK != 0)
            Utils::ExitWithUsrMsg(string("could not rename file ") +
                                      tmp_file_name + " into " +
                                      dest_file_name);
        }
      }
    }
  }
}


    void Stats::delete_indivs(int64_t gen) {
      //std::cout << "Deleting gen : " << gen << '\n';
      /*  int i = 0;
      for(auto indiv_it = indivs_[gen].begin() ; indiv_it != indivs_[gen].end() ; ++indiv_it) {
          if ((*indiv_it) != nullptr) {
              if ((*indiv_it)->number_of_clones_ == 0) {
                  printf("Indiv %d is %d -- %p\n",i,(*indiv_it)->number_of_clones_,*indiv_it);
                  delete *indiv_it;
                  (*indiv_it) = nullptr;
              }
          }
          i++;
      }*/

      if ((!ExpManager_7::standalone_simd) || (ExpManager_7::standalone_simd && exp_m_->check_simd())){
          std::unordered_map<unsigned long long, Individual *> unique_individual;

          for (auto indiv : indivs_[gen]) {
              unique_individual[indiv->long_id()] = indiv;
              indiv->number_of_clones_--;
          }

          for (auto indiv_it = unique_individual.begin(); indiv_it != unique_individual.end(); ++indiv_it) {
              if ((indiv_it->second)->number_of_clones_ == 0) {
                  delete indiv_it->second;
              }
          }
      } else {
          for (auto indiv_it = indivs_[gen].begin(); indiv_it != indivs_[gen].end(); ++indiv_it) {

              delete *indiv_it;
          }
      }


      indivs_.erase(gen);
    }

    Individual* Stats::b_indiv(int64_t gen) const {
      if (gen != AeTime::time()) {
        // printf("gen: %lld, cur gen: %lld\n", gen, AeTime::time());
        return nullptr;
      }
      return exp_m_->best_indiv();
      // if(gen < 1) {
      //   return exp_m_->best_indiv();
      // }
      // for(const auto& indiv : indivs_.at(gen)) {
      //   if(indiv->rank() == exp_m_->nb_indivs())
      //     return indiv;
      // }
      // return nullptr;
    }

/**
 * Returns a list of all the individuals with their replication report
 */
    std::list<std::pair<Individual*, ReplicationReport*>>
    Stats::indivs_annotated(int64_t gen) const {
      std::list<std::pair<Individual*, ReplicationReport*>> annotated_list;
      for (const auto& indiv : exp_m_->indivs()) {
        annotated_list.emplace_back(indiv, exp_m_->tree() ?
                                           exp_m_->tree()->report_by_index(gen, indiv->id()) : nullptr);
      }
      return annotated_list;
    }

} // namespace aevol
