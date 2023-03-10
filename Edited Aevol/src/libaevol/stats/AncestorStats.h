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


#ifndef AEVOL_ANCESTORSTATS_H_
#define AEVOL_ANCESTORSTATS_H_


// =================================================================
//                              Libraries
// =================================================================
#include <string>


// =================================================================
//                            Project Files
// =================================================================
#include "ExpManager.h"
//#include "Individual.h"
//#include "PhenotypicTargetHandler.h"
//#include "Stats.h"
//#include "ReplicationReport.h"

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================
class Individual;
class PhenotypicTargetHandler;
class ReplicationReport;


enum types {
  ENVIRONMENT_STAT    = 0,
  TERMINATOR_STAT     = 1,
  ZONES_STAT          = 2,
  OPERONS_STAT        = 3,
  FIXED_MUTATION_STAT = 4,
  NB_TYPES            = 5
};


class AncestorStats {
 public :
  // =================================================================
  //                             Constructors
  // =================================================================
  AncestorStats();


  // =================================================================
  //                             Destructor
  // =================================================================
  virtual ~AncestorStats();


  // =================================================================
  //                            Public Methods
  // =================================================================

  void setup_anc_stat(int64_t t_mrca);

  void setup_anc_indiv(Individual* indiv);

  void process_evolution(ReplicationReport* rep, int64_t t);

  void Flush();
  void Close();

 protected :
  // =================================================================
  //                           Protected Methods
  // =================================================================

  void Set_file_names();

  void Open_files(int64_t time);

  void Write_headers();

  void Resume_stats(int64_t time);


  static FILE* open_environment_stat_file(const char* prefix, const char* postfix);
  static void write_environment_stats(int64_t t,
                               const std::shared_ptr<PhenotypicTargetHandler> pth,
                               FILE* env_file);
  static FILE* open_terminators_stat_file(const char* prefix, const char* postfix);
  static void write_terminators_stats(int64_t t, Individual* indiv,
                               FILE* terminator_file);
  static FILE* open_zones_stat_file(const char* prefix, const char* postfix);
  static void write_zones_stats(int64_t t,
                         Individual* indiv,
                         const std::shared_ptr<PhenotypicTargetHandler> phenotypicTargetHandler,
                         FILE* zone_file);
  static FILE* open_operons_stat_file(const char* prefix, const char* postfix);
  static void write_operons_stats(int64_t t, Individual* indiv, FILE* operon_file);
  static FILE* open_fixed_mutations_file(const char* prefix, const char* postfix);

  void save_grid_cell(std::string grid_file_name = STATS_DIR"/ancestor_stats/saved_gridcell.ae") const;
  void load_grid_cell(ExpManager* exp_m);

  // =================================================================
  //                          Protected Attributes
  // =================================================================

  ExpManager* exp_m_;
  std::shared_ptr<PhenotypicTargetHandler> phenotypicTH_;

  bool b_zones_;

  FILE ** stat_files_;
  char ** stat_files_names_;
  Stats* mystats_;

  Individual* indiv_;
  GridCell* grid_cell_;
};

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol
#endif // AEVOL_ANCESTORSTATS_H_
