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
#include <cassert>
#include <cstring>



// =================================================================
//                            Project Files
// =================================================================
#include "ae_string.h"
#include "ExpSetup.h"

namespace aevol {

//##############################################################################
//                                                                             #
//                               Class ae_string                               #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
ae_string::ae_string() {
  nb_blocks_ = 1;
  length_ = 0;
  data_ = new char[nb_blocks_ * BLOCK_SIZE * sizeof(char)];
  data_[length_] = '\0';
}

ae_string::ae_string(const ae_string &model) {
  nb_blocks_ = model.nb_blocks_;
  length_ = model.length_;
  data_ = new char[nb_blocks_ * BLOCK_SIZE * sizeof(char)];

  memcpy(data_, model.data_, (length_ +1) * sizeof(char));
}

/*!
  Creates a new random string with the given length.
*/
ae_string::ae_string(int32_t length, std::shared_ptr<JumpingMT> prng) {
  nb_blocks_ = nb_blocks(length);
  length_ = length;
  data_ = new char[nb_blocks_ * BLOCK_SIZE];

  // Generate a random genome
  for (int32_t i = 0 ; i < length_; i++) {
    data_[i] = '0' + prng->random(NB_BASE);
  }
  data_[length_] = '\0';
}

/**
 * Creates a new ae_string with sequence <seq> (having length <length>)
 */
ae_string::ae_string(const char* seq, int32_t length) {
  length_ = length;
  nb_blocks_ = nb_blocks(length);
  data_ = new char[nb_blocks_ * BLOCK_SIZE];
  memcpy(data_, seq, (length+1) * sizeof(char));
}

/**
 * Creates a new ae_string with sequence <seq> (having length <length>).
 * WARNING : <seq> is used directly which means the caller must not delete it.
 */
ae_string::ae_string(char* seq, int32_t length, bool use_seq) {
  assert(use_seq);

  length_ = length;
  nb_blocks_ = nb_blocks(length);
  data_ = seq;
}

ae_string::ae_string(gzFile backup_file) {
  gzread(backup_file, &nb_blocks_,  sizeof(nb_blocks_));
  data_ = new char[nb_blocks_ * BLOCK_SIZE];
  gzread(backup_file, &length_, sizeof(length_));
  gzread(backup_file, data_,
         static_cast<unsigned int>((length_ + 1) * sizeof(*data_)));
}

ae_string::ae_string(char* organism_file_name) {
  FILE* org_file = fopen(organism_file_name, "r");
  int32_t length;

  if (org_file == NULL) {
    printf("Error opening organism file\n");
  }
  else {
    fseek(org_file , 0 , SEEK_END);
    length = ftell(org_file);
    rewind(org_file);

    nb_blocks_ = nb_blocks(length);
    length_ = length;
    data_ = new char[nb_blocks_ * BLOCK_SIZE];
    for (int32_t i = 0 ; i < length_ -1 ; i++) {
      data_[i] = static_cast<char>(fgetc(org_file));
    }

    data_[length_] = '\0';
  }

  fclose(org_file);
}

// =================================================================
//                             Destructors
// =================================================================
ae_string::~ae_string() {
  delete [] data_;
}

// =================================================================
//                            Public Methods
// =================================================================
/**
 * Remove the sequence between positions 'pos_1' and 'pos_2'
 *
 * Character at pos_1 is removed
 * Character at pos_2 is not removed unless pos_1 == pos_2, in which case the
 * string will become empty
 */
void ae_string::remove(int32_t pos_1, int32_t pos_2) {
  assert(pos_1 >= 0 && pos_2 >= pos_1 && pos_2 <= length_);

  // Compute size of new genome
  int32_t new_length    = length_ - (pos_2 - pos_1);
  int32_t new_nb_blocks = nb_blocks(new_length);
  char*   new_genome    = new char[new_nb_blocks * BLOCK_SIZE * sizeof(char)];

  // Copy the remaining of the genome in tmp (preceeding and following parts)
  memcpy(new_genome, data_, pos_1 * sizeof(char));
  memcpy(&new_genome[pos_1], &data_[pos_2],
         (new_length - pos_1) * sizeof(char));
  new_genome[new_length] = '\0';

  // Replace previous genome with the new one
  delete [] data_;
  data_ = new_genome;

  // Update length data
  length_ = new_length;
  nb_blocks_ = new_nb_blocks;
}

void ae_string::insert(int32_t pos, const char* seq, int32_t seq_length) {
// Insert sequence 'seq' at position 'pos'
  assert(pos >= 0 && pos < length_);

  // If the sequence's length was not provided, compute it
  if (seq_length == -1) {
    seq_length = strlen(seq);
  }

  // Compute size of new genome
  int32_t new_length    = length_ + seq_length;
  int32_t new_nb_blocks = nb_blocks(new_length);
  char*   new_genome    = new char[new_nb_blocks * BLOCK_SIZE * sizeof(char)];

  // Build new genome from previous genome and sequence to insert
  memcpy(new_genome, data_, pos * sizeof(char));
  memcpy(&new_genome[pos], seq, seq_length * sizeof(char));
  memcpy(&new_genome[pos+seq_length], &data_[pos],
         (length_ - pos) * sizeof(char));
  new_genome[new_length] = '\0';

  // Replace the previous genome with the new one
  delete [] data_;
  data_ = new_genome;

  // Update length-related data
  length_ = new_length;
  nb_blocks_ = new_nb_blocks;
}

void ae_string::replace(int32_t pos, char* seq, int32_t seq_length) {
// Invert the sequence between positions 'first' and 'last'
  // Check pos value
  assert(pos >= 0 && pos < length_);

  // If the sequence's length was not provided, compute it
  if (seq_length == -1) {
    seq_length = strlen(seq);
  }

  // Check that the sequence is contiguous
  assert(pos + seq_length <= length_);

  // Perform the replacement
  memcpy(&data_[pos], seq, seq_length * sizeof(char));
}

void ae_string::save(gzFile backup_file) {
  gzwrite(backup_file, &nb_blocks_, sizeof(nb_blocks_));
  gzwrite(backup_file, &length_, sizeof(length_));
  gzwrite(backup_file, data_,
         static_cast<unsigned int>((length_ + 1) * sizeof(data_[0])));
}


// =================================================================
//                           Protected Methods
// =================================================================
} // namespace aevol
