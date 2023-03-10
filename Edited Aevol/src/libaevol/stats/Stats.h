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


/*! \class Stats
    \brief Manage statistics files
*/


#ifndef AEVOL_STATS_H_
#define AEVOL_STATS_H_


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdio.h>



// =================================================================
//                            Project Files
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Individual.h"
#include "ReplicationReport.h"


namespace aevol {


// =================================================================
//                          Class declarations
// =================================================================
class ExpManager;






class Stats
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    Stats() = delete;
    Stats(const Stats &) = delete;
    Stats(const std::string prefix, const std::string postfix = "",
          bool best_indiv_only = false);
    Stats(ExpManager* exp_m,
          bool best_indiv_only = false,
          const std::string prefix = "stat",
          const std::string postfix = "",
          bool with_plasmids = false,
          bool compute_phen_contrib_by_GU = false);
    Stats(ExpManager* exp_m,
          int64_t time,
          bool best_indiv_only = false,
          const std::string prefix = "stat",
          const std::string postfix = "",
          bool addition_old_stats = true,
          bool delete_old_stats = true);


    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~Stats();

    // =================================================================
    //                        Accessors: getters
    // =================================================================

    // =================================================================
    //                        Accessors: setters
    // =================================================================
    inline void set_exp_m(ExpManager * exp_m);

    // =================================================================
    //                            Public Methods
    // =================================================================
    void write_current_generation_statistics(int64_t gen);
    void write_statistics_of_this_indiv(int64_t time, Individual* indiv,
                                        ReplicationReport* replic_report);

    void flush();

    void write_headers(bool ancstats_stats = false);

    void CreateTmpFiles(int64_t time);
    void MoveTmpFiles(const std::string& destdir);
    void PromoteTmpFiles();

    void add_indivs(int64_t gen, const std::list<Individual*> indivs);
    void add_indivs(int64_t gen, Individual_7** indivs);
    void delete_indivs(int64_t gen);


 protected :
    // =================================================================
    //                           Protected Methods
    // =================================================================
    void init_data();
    void set_file_names(const std::string prefix,
                            const std::string postfix,
                            bool best_indiv_only,
                            bool with_plasmids = false,
                            bool compute_phen_contrib_by_GU = false);
    void open_files();

    Individual* b_indiv(int64_t gen) const;
    std::list<std::pair<Individual*, ReplicationReport*>> indivs_annotated(int64_t gen) const;

    inline void write_header(FILE* file_name, const char* header);
    inline void write_header(FILE* file_name, const char* header, int8_t key);

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    ExpManager * exp_m_;

        // 3D tables of stat files (FILE*) and their names (char*)
        // Dimensions are given by:
        //    * genetic unit (ALL_GU, CHROM or PLASMIDS)
        //    * BEST or GLOB
        //    * stat type (FITNESS_STATS, MUTATION_STATS, GENES_STATS, BP_STATS or REAR_STATS)
        // Files that are not wanted MUST have their name set to NULL.
        // The files themselves are also NULL because we don't fopen() them.
    FILE**** stat_files_;
    char**** stat_files_names_;

    //indivs
    std::map<int64_t, std::list<Individual*>> indivs_;
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================
inline void Stats::set_exp_m(ExpManager * exp_m)
{
  exp_m_ = exp_m;
}

inline void Stats::write_header(FILE* file_name, const char* header)
{
  if (file_name != NULL) fprintf(file_name, "# %s\n", header);
}

inline void Stats::write_header(FILE* file_name, const char* header, int8_t key)
{
  if (file_name != NULL) fprintf(file_name, "# %2d. %s\n", key, header);
}


} // namespace aevol
#endif // AEVOL_STATS_H_
