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


#ifndef AEVOL_LOGS_H_
#define AEVOL_LOGS_H_


// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <zlib.h>



// =================================================================
//                            Project Files
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "ae_enums.h"

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================



class Logging
{
  public :

    // =================================================================
    //                             Constructors
    // =================================================================
    Logging();

    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~Logging();

    // =================================================================
    //                              Accessors
    // =================================================================
    inline FILE* log(LogType log_type)   const;
    inline int8_t logs() const;
    inline bool  is_logged(LogType log_type) const;

    // =================================================================
    //                              Operators
    // =================================================================

    // =================================================================
    //                            Public Methods
    // =================================================================
    //void save(gzFile backup_file) const;
    void load(int8_t logs, int32_t num_gener);
    void print_to_file(FILE* file) const;

    void set_logs(int8_t logs);
    void flush();

    // =================================================================
    //                           Public Attributes
    // =================================================================





  protected :

    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    /*    Logging()
    {
      printf("ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    };
    Logging(const Logging &model)
    {
      printf("ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    };*/


    // =================================================================
    //                           Protected Methods
    // =================================================================
    void write_headers() const;

    // =================================================================
    //                          Protected Attributes
    // =================================================================
    int8_t  logs_; // Which logs are "turned on" (bitmap)
    FILE*   transfer_log_;
    FILE*   rear_log_;
    FILE*   barrier_log_;
    //FILE*   param_modification_log_;
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================
inline FILE*Logging::log(LogType log_type) const
{
  switch (log_type)
  {
    case LOG_TRANSFER :
    {
      return transfer_log_;
    }
    case LOG_REAR :
    {
      return rear_log_;
    }
    case LOG_BARRIER :
    {
      return barrier_log_;
    }
    /*case LOG_LOADS :
    {
      return param_modification_log_;
    }*/
    default:
    {
      printf("ERROR: unknown log_type in file %s : l%d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }
}

inline int8_t Logging::logs() const
{
  return logs_;
}

inline bool Logging::is_logged(LogType log_type) const
{
  switch (log_type)
  {
    case LOG_TRANSFER :
    {
      return (logs_ & LOG_TRANSFER);
    }
    case LOG_REAR :
    {
      return (logs_ & LOG_REAR);
    }
    case LOG_BARRIER :
    {
      return (logs_ & LOG_BARRIER);
    }
    /*case LOG_LOADS :
    {
      return (logs_ & LOG_LOADS);
    }*/
    default:
    {
      printf("ERROR: unknown log_type in file %s : l%d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
  }
}

// =====================================================================
//                          Operators' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol

#endif // AEVOL_LOGS_H_
