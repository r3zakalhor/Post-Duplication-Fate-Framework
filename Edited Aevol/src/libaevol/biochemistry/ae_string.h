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


#ifndef AEVOL_STRING_H_
#define AEVOL_STRING_H_


// =================================================================
//                              Libraries
// =================================================================
#include <cinttypes>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <zlib.h>

#include <memory>



// =================================================================
//                            Project Files
// =================================================================
#include "Utils.h"
#include "JumpingMT.h"

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================







#define BLOCK_SIZE INT32_C(1024)

class ae_string {
 public :

  // =================================================================
  //                             Constructors
  // =================================================================
  ae_string();
  ae_string(const ae_string &model);
  ae_string(int32_t length, std::shared_ptr<JumpingMT> prng);
  ae_string(const char* seq, int32_t length);
  ae_string(char* seq, int32_t length, bool use_seq);
  explicit ae_string(gzFile backup_file);
  explicit ae_string(char* organism_file_name);

  // =================================================================
  //                             Destructors
  // =================================================================
  virtual ~ae_string();

  // =================================================================
  //                              Accessors
  // =================================================================
  const char* data() const {return data_;}
  int32_t length() const {return length_;}
  inline void set_data(char* data, int32_t length = -1);

  // =================================================================
  //                            Public Methods
  // =================================================================
  void remove(int32_t first, int32_t last);
  void insert(int32_t pos, const char* seq, int32_t seq_length = -1);
  void replace(int32_t pos, char* seq, int32_t seq_length = -1);

  void save(gzFile backup_file);





 protected :
  // =================================================================
  //                           Protected Methods
  // =================================================================
  static inline int32_t nb_blocks(int32_t length);

  // =================================================================
  //                          Protected Attributes
  // =================================================================

  char* data_ = nullptr;
  int32_t length_ = -1;
  int32_t nb_blocks_ = -1;
};


// =====================================================================
//                          Accessors' definitions
// =====================================================================
void ae_string::set_data(char* data, int32_t length /* = -1 */) {
  if (data_ != NULL) {
    delete [] data_;
    data_ = NULL;
  }

  data_ = data;
  length_ = (length != -1) ? length : strlen(data_);
  nb_blocks_ = nb_blocks(length_);
}

// =====================================================================
//                       Inline functions' definition
// =====================================================================
int32_t ae_string::nb_blocks(int32_t length) {
  return length/BLOCK_SIZE + 1;
}

} // namespace aevol

#endif // AEVOL_STRING_H_
