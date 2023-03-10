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



// =================================================================
//                            Project Files
// =================================================================
#include "Logging.h"
#include "ExpSetup.h"

namespace aevol {

//##############################################################################
//                                                                             #
//                                Class Logging                                #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================

Logging::Logging()
{
  logs_ = 0;

  transfer_log_         = NULL;
  rear_log_             = NULL;
  barrier_log_          = NULL;
  //param_modification_log_ = NULL;
}

// =================================================================
//                             Destructors
// =================================================================
Logging::~Logging()
{
  if (logs_ & LOG_TRANSFER)
  {
    fclose(transfer_log_);
  }
  if (logs_ & LOG_REAR)
  {
    fclose(rear_log_);
  }
  if (logs_ & LOG_BARRIER)
  {
    fclose(barrier_log_);
  }
  /*if (logs_ & LOG_LOADS)
  {
    fclose(param_modification_log_);
  }*/
}

// =================================================================
//                            Public Methods
// =================================================================
/*void Logging::save(gzFile setup_file) const
{
  gzwrite(backup_file, &logs_, sizeof(logs_));
}*/

void Logging::load(int8_t logs, int32_t num_gener)
{
  char* line = new char[500];

  logs_ = logs;

  // Prepare required log files
  if (logs_ & LOG_TRANSFER)
  {
    rename ("log_transfer.out", "log_transfer.out.old");
    FILE* old_transfer_log = fopen("log_transfer.out.old", "r");
    if (old_transfer_log == NULL)
    {
      printf("Error: Failed to open file \"log_transfer.out.old\"\n");
      exit(EXIT_FAILURE);
    }

    transfer_log_ = fopen("log_transfer.out", "w");
    if (transfer_log_ == NULL)
    {
      printf("Error: Failed to open file \"log_transfer.out\"\n");
      exit(EXIT_FAILURE);
    }

    // Copy file headers
    if (fgets(line, 500, old_transfer_log) == NULL) {
      // TODO check for error
    }
    while (!feof(old_transfer_log) && line[0] == '#')
    {
      fputs(line, transfer_log_);
      if (fgets(line, 500, old_transfer_log) == NULL) {
	// TODO check for error
      }
    }
    // This is the empty line between the header and the values
    //fputs(line, transfer_log_);

    // Copy log entries until num_gener (excluded)
    while ((int32_t)atol(line) < num_gener && !feof(old_transfer_log))
    {
      fputs(line, transfer_log_);
      if (fgets(line, 500, old_transfer_log) == NULL) {
	// TODO check for error
      }
      while(!feof(old_transfer_log) & (line[0] == '\t' || line[0] == ' '))
      {
        fputs(line, transfer_log_);
        if (fgets(line, 500, old_transfer_log) == NULL) {
	  // TODO check for error
	}
      }
    }

    fclose(old_transfer_log);
    remove("log_transfer.out.old");
  }
  if (logs_ & LOG_REAR)
  {
    rename ("log_rear.out", "log_rear.out.old");
    FILE* old_rear_log = fopen("log_rear.out.old", "r");
    if (old_rear_log == NULL)
    {
      printf("Error: Failed to open file \"log_rear.out.old\"\n");
      exit(EXIT_FAILURE);
    }

    rear_log_ = fopen("log_rear.out", "w");
    if (rear_log_ == NULL)
    {
      printf("Error: Failed to open file \"log_rear.out\"\n");
      exit(EXIT_FAILURE);
    }

    // Copy file headers
    if (fgets(line, 500, old_rear_log) == NULL) {
      // TODO check for error
    }
    while (!feof(old_rear_log) && line[0] == '#')
    {
      fputs(line, rear_log_);
      if (fgets(line, 500, old_rear_log) == NULL) {
	// TODO check for error
      }
    }
    // This is the empty line between the header and the values
    //fputs(line, rear_log_);

    // Copy log entries until num_gener (excluded)
    while ((int32_t)atol(line) < num_gener && !feof(old_rear_log))
    {
      fputs(line, rear_log_);
      if (fgets(line, 500, old_rear_log) == NULL) {
	// TODO check for error
      }
    }

    fclose(old_rear_log);
    remove("log_rear.out.old");
  }
  if (logs_ & LOG_BARRIER)
  {
    rename ("log_barrier.out", "log_barrier.out.old");
    FILE* old_barrier_log = fopen("log_barrier.out.old", "r");
    if (old_barrier_log == NULL)
    {
      printf("Error: Failed to open file \"log_barrier.out.old\"\n");
      exit(EXIT_FAILURE);
    }

    barrier_log_ = fopen("log_barrier.out", "w");
    if (barrier_log_ == NULL)
    {
      printf("Error: Failed to open file \"log_barrier.out\"\n");
      exit(EXIT_FAILURE);
    }

    // Copy file headers
    if (fgets(line, 500, old_barrier_log) == NULL) {
      // TODO check for error
    }
    while (!feof(old_barrier_log) && line[0] == '#')
    {
      fputs(line, barrier_log_);
      if (fgets(line, 500, old_barrier_log) == NULL) {
	// TODO check for error
      }
    }
    // This is the empty line between the header and the values
    //fputs(line, barrier_log_);

    // Copy log entries until num_gener (excluded)
    while ((int32_t)atol(line) < num_gener && !feof(old_barrier_log))
    {
      fputs(line, barrier_log_);
      if (fgets(line, 500, old_barrier_log) == NULL) {
	// TODO check for error
      }
    }

    fclose(old_barrier_log);
    remove("log_barrier.out.old");
  }
  /*if (logs_ & LOG_LOADS)
  {
    rename ("log_param_modification.out", "log_param_modification.out.old");
    FILE* old_param_modification_log = fopen("log_param_modification.out.old", "r");
    if (old_param_modification_log == NULL)
    {
      printf("Error: Failed to open file \"log_param_modification.out.old\"\n");
      exit(EXIT_FAILURE);
    }

    param_modification_log_ = fopen("log_param_modification.out", "w");
    if (param_modification_log_ == NULL)
    {
      printf("Error: Failed to open file \"log_param_modification.out\"\n");
      exit(EXIT_FAILURE);
    }

    //Copy file headers
    ret = fgets(line, 500, old_param_modification_log);
    while (!feof(old_param_modification_log) && line[0] == '#')
    {
      fputs(line, param_modification_log_);
      ret = fgets(line, 500, old_param_modification_log);
    }
    // This is the empty line between the header and the values
    //fputs(line, param_modification_log_);

    // Copy log entries until num_gener (excluded)
    while ((int32_t)atol(line) < num_gener && !feof(old_param_modification_log))
    {
      fputs(line, param_modification_log_);
      ret = fgets(line, 500, old_param_modification_log);
    }

    fclose(old_param_modification_log);
    remove("log_param_modification.out.old");
  }*/

  delete [] line;
}

void Logging::print_to_file(FILE* file) const
{
  fprintf(file, "logs        :                %" PRId8 "\n", logs_);
}

void Logging::set_logs(int8_t logs)
{
  logs_ = logs;

  // Open required log files
  if (logs_ & LOG_TRANSFER)
  {
    transfer_log_ = fopen("log_transfer.out", "w");
    if (transfer_log_ == NULL)
    {
      printf("Error: Failed to open file \"log_transfer.out\"\n");
      exit(EXIT_FAILURE);
    }
  }
  if (logs_ & LOG_REAR)
  {
    rear_log_ = fopen("log_rear.out", "w");
    if (rear_log_ == NULL)
    {
      printf("Error: Failed to open file \"log_rear.out\"\n");
      exit(EXIT_FAILURE);
    }
  }
  if (logs_ & LOG_BARRIER)
  {
    barrier_log_ = fopen("log_barrier.out", "w");

    if (barrier_log_ == NULL)
    {
      printf("Error: Failed to open file \"log_barrier.out\"\n");
      exit(EXIT_FAILURE);
    }
  }
  /*if (logs_ & LOG_LOADS)
  {
    param_modification_log_ = fopen("log_param_modification.out", "w");

    if (param_modification_log_ == NULL)
    {
      printf("Error: Failed to open file \"log_param_modification.out\"\n");
      exit(EXIT_FAILURE);
    }
  }*/

  this->write_headers();
}

void Logging::flush()
{
  if (logs_ & LOG_TRANSFER)
  {
    fflush(transfer_log_);
  }
  if (logs_ & LOG_REAR)
  {
    fflush(rear_log_);
  }
  if (logs_ & LOG_BARRIER)
  {
    fflush(barrier_log_);
  }
  /*if (logs_ & LOG_LOADS)
  {
    fflush(param_modification_log_);
  }*/
}

// =================================================================
//                           Protected Methods
// =================================================================
void Logging::write_headers() const
{
  // ========== TRANSFER LOG ==========
  if (logs_ & LOG_TRANSFER)
  {
    fprintf(transfer_log_, "######################################################################\n");
    fprintf(transfer_log_, "#                 Horizontal transfer log\n");
    fprintf(transfer_log_, "#\n");
    fprintf(transfer_log_, "# Log of every horizontal transfer that occured during the simulation\n");
    fprintf(transfer_log_, "#\n");
    fprintf(transfer_log_, "# 1.  Generation\n");
    fprintf(transfer_log_, "# 2.  Index of the recepient\n");
    fprintf(transfer_log_, "# 3.  Index of the donor (generation n-1)\n");
    fprintf(transfer_log_, "# 4.  Type of transfer\n");
    fprintf(transfer_log_, "# 5.  Length of the transferred segment\n");
    fprintf(transfer_log_, "# 6.  Length of the replaced segment (if any)\n");
    fprintf(transfer_log_, "# 7.  Size of the genome before the transfer\n");
    fprintf(transfer_log_, "# 8.  Size of the genome after the transfer\n");
    fprintf(transfer_log_, "# 9.  Alignment 1 point 1\n");
    fprintf(transfer_log_, "# 10. Alignment 1 point 2\n");
    fprintf(transfer_log_, "# 11. Alignment 1 score\n");
    fprintf(transfer_log_, "# 12. Alignment 2 point 1\n");
    fprintf(transfer_log_, "# 13. Alignment 2 point 2\n");
    fprintf(transfer_log_, "# 14. Alignment 2 score\n");
    fprintf(transfer_log_, "#\n");
    fprintf(transfer_log_, "######################################################################\n");
    fprintf(transfer_log_, "#\n");
    fprintf(transfer_log_, "# Header for R\n");
    fprintf(transfer_log_, "gener recepient donor t_type seg_len replaced_len size_before size_after align1_pt1 align1_pt2 score1 align2_pt1 align2_pt2 score2\n");
    fprintf(transfer_log_, "#\n");
  }

  // ========== REAR LOG ==========
  if (logs_ & LOG_REAR)
  {
    fprintf(rear_log_, "######################################################################\n");
    fprintf(rear_log_, "#                 Chromosomal rearrangement log\n");
    fprintf(rear_log_, "#\n");
    fprintf(rear_log_, "# Log of every rearrangement that occured during the simulation\n");
    fprintf(rear_log_, "# (not just one lineage)\n");
    fprintf(rear_log_, "#\n");
    fprintf(rear_log_, "# 1. Generation\n");
    fprintf(rear_log_, "# 2. Index of the individual that has undergone the rearrangement\n");
    fprintf(rear_log_, "# 3. Type of rearrangement\n");
    fprintf(rear_log_, "# 4. Length of the rearranged segment\n");
    fprintf(rear_log_, "# 5. Size of the genome before the rearrangement\n");
    fprintf(rear_log_, "# 6. Alignment score that was needed for this rearrangement to occur\n");
    fprintf(rear_log_, "# 7. Second alignment score (translocation only)\n");
    fprintf(rear_log_, "#\n");
    fprintf(rear_log_, "######################################################################\n");
    fprintf(rear_log_, "#\n");
    fprintf(rear_log_, "# Header for R\n");
    fprintf(rear_log_, "gener indiv r_type seg_len genome_size score1 score2\n");
    fprintf(rear_log_, "#\n");
  }

  // ========== BARRIER LOG ==========
  if (logs_ & LOG_BARRIER)
  {
    fprintf(barrier_log_, "######################################################################\n");
    fprintf(barrier_log_, "#                     Genome size limits log\n");
    fprintf(barrier_log_, "#\n");
    fprintf(barrier_log_, "# An entry is written whenever a mutation would have produced a\n");
    fprintf(barrier_log_, "# genome whose size wouldn't lie in [min, max].\n");
    fprintf(barrier_log_, "# The corresponding mutation is \"cancelled\"\n");
    fprintf(barrier_log_, "#\n");
    fprintf(barrier_log_, "# 1. Generation\n");
    fprintf(barrier_log_, "# 2. Index of the individual\n");
    fprintf(barrier_log_, "# 3. Type of event\n");
    fprintf(barrier_log_, "# 4. Segment length\n");
    fprintf(barrier_log_, "# 5. Replaced segment length\n");
    fprintf(barrier_log_, "# 6. GU size (before the event)\n");
    fprintf(barrier_log_, "# 7. Genome size (before the event)\n");
    fprintf(barrier_log_, "#\n");
    fprintf(barrier_log_, "######################################################################\n");

  }

  // ========== LOADS LOG ==========
  /*if (logs_ & LOG_LOADS)
  {
    fprintf(param_modification_log_, "######################################################################\n");
    fprintf(param_modification_log_, "#                     Parameter modification log\n");
    fprintf(param_modification_log_, "#\n");
    fprintf(param_modification_log_, "# An entry is written whenever a parameter is modified by aevol_modify.\n");
    fprintf(param_modification_log_, "######################################################################\n");
  }*/
}
} // namespace aevol
