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
#include <cstdlib>

#include "Mutation.h"
#include "PointMutation.h"
#include "SmallInsertion.h"
#include "SmallDeletion.h"
#include "Duplication.h"
#include "Deletion.h"
#include "Translocation.h"
#include "Inversion.h"
#include "InsertionHT.h"
#include "ReplacementHT.h"

namespace aevol {


//#############################################################################
//
//                                Class Mutation
//
//#############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
Mutation* Mutation::Load(gzFile backup_file) {
  // Retrieve mutation type
  int8_t tmp_mut_type;
  gzread(backup_file, &tmp_mut_type,  sizeof(tmp_mut_type));
  MutationType mut_type = (MutationType) tmp_mut_type;

  // Call the appropriate constructor accordingly
  Mutation* mut;


      // printf("Mut type %d\n",mut_type);

  switch (mut_type) {
    case SWITCH :
      mut = new PointMutation();
      break;
    case S_INS :
      mut = new SmallInsertion();
      break;
    case S_DEL :
      mut = new SmallDeletion();
      break;
    case DUPL :
      mut = new Duplication();
      break;
    case DEL :
      mut = new Deletion();
      break;
    case TRANS :
      mut = new Translocation();
      break;
    case INV :
      mut = new Inversion();
      break;
    case INS_HT :
      mut = new InsertionHT();
      break;
    case REPL_HT :
      mut = new ReplacementHT();
      break;
    default :
      Utils::ExitWithDevMsg("invalid mutation type ", __FILE__, __LINE__);
      exit(-1); // Superfluous but suppresses a warning
  }

  // Load from backup file
  mut->load(backup_file);
  return mut;
}




//Mutation::Mutation(gzFile backup_file)
//{
//  pos_ = NULL;
//  length_ = NULL;
//  seq_ = NULL;
//  align_score_ = NULL;
//
//  int8_t tmp_mut_type;
//  gzread(backup_file, &tmp_mut_type,  sizeof(tmp_mut_type));
//  mut_type_ = (MutationType) tmp_mut_type;
//  //~ printf("mut type %d\n", mut_type_);
//
//  switch (mut_type_)
//  {
//    case SWITCH :
//    {
//      pos_ = new int32_t;
//      gzread(backup_file, pos_,  sizeof(*pos_));
//      break;
//    }
//    case S_INS :
//    {
//      pos_ = new int32_t;
//      gzread(backup_file, pos_,      sizeof(*pos_));
//      length_ = new int32_t;
//      gzread(backup_file, length_,  sizeof(*length_));
//
//      seq_ = new char[length_[0] + 1];
//      gzread(backup_file, seq_,  length_[0] * sizeof(seq_[0]));
//      seq_[length_[0]] = '\0';
//      break;
//    }
//    case S_DEL :
//    {
//      pos_ = new int32_t;
//      gzread(backup_file, pos_,      sizeof(*pos_));
//      length_ = new int32_t;
//      gzread(backup_file, length_,  sizeof(*length_));
//      break;
//    }
//    case DUPL :
//    {
//      pos_ = new int32_t[3];
//      gzread(backup_file, pos_,  3 * sizeof(pos_[0]));
//      length_ = new int32_t;
//      gzread(backup_file, length_,  sizeof(*length_));
//      align_score_ = new int16_t;
//      gzread(backup_file, align_score_, sizeof(*align_score_));
//
//      break;
//    }
//    case DEL :
//    {
//      pos_ = new int32_t[2];
//      gzread(backup_file, pos_,  2 * sizeof(pos_[0]));
//      length_ = new int32_t;
//      gzread(backup_file, length_,  sizeof(*length_));
//      align_score_ = new int16_t;
//      gzread(backup_file, align_score_, sizeof(*align_score_));
//
//      break;
//    }
//    case TRANS :
//    {
//      pos_ = new int32_t[4];
//      gzread(backup_file, pos_,  4 * sizeof(pos_[0]));
//      int8_t tmp_invert;
//      gzread(backup_file, &tmp_invert,  sizeof(tmp_invert));
//      invert_ = (tmp_invert != 0);
//      length_ = new int32_t;
//      gzread(backup_file, length_,  sizeof(*length_));
//      align_score_ = new int16_t[2];
//      gzread(backup_file, align_score_, 2 * sizeof(align_score_[0]));
//
//      break;
//    }
//    case INV :
//    {
//      pos_ = new int32_t[2];
//      gzread(backup_file, pos_,  2 * sizeof(pos_[0]));
//      length_ = new int32_t;
//      gzread(backup_file, length_,  sizeof(*length_));
//      align_score_ = new int16_t;
//      gzread(backup_file, align_score_, sizeof(*align_score_));
//
//      break;
//    }
//    case INSERT:
//    {
//      pos_ = new int32_t;
//      gzread(backup_file, pos_,  sizeof(*pos_));
//     length_ = new int32_t;
//      gzread(backup_file, length_,  sizeof(*length_));
//      seq_ = new char[length_[0] + 1];
//      gzread(backup_file, seq_,  length_[0] * sizeof(seq_[0]));
//      seq_[length_[0]] = '\0';
//    }
//    case INS_HT:
//    {
//      pos_ = new int32_t[4];
//      gzread(backup_file, pos_,  4 * sizeof(pos_[0]));
//      length_ = new int32_t;
//      gzread(backup_file, length_,  sizeof(*length_));
//      seq_ = new char[length_[0] + 1];
//      gzread(backup_file, seq_,  length_[0] * sizeof(seq_[0]));
//      seq_[length_[0]] = '\0';
//      align_score_ = new int16_t[2];
//      gzread(backup_file, align_score_, 2 * sizeof(align_score_[0]));
//      gzread(backup_file, &donor_id_,  sizeof(donor_id_));
//      gzread(backup_file, &sense_,  sizeof(sense_));
//      break;
//    }
//    case REPL_HT:
//    {
//      pos_ = new int32_t[4];
//      gzread(backup_file, pos_,  4 * sizeof(pos_[0]));
//      length_ = new int32_t[2];
//      gzread(backup_file, length_, 2 *sizeof(length_[0]));
//      seq_ = new char[length_[1] + 1];
//      gzread(backup_file, seq_,  length_[1] * sizeof(seq_[0]));
//      seq_[length_[1]] = '\0';
//      align_score_ = new int16_t[2];
//      gzread(backup_file, align_score_, 2 * sizeof(align_score_[0]));
//      gzread(backup_file, &donor_id_,  sizeof(donor_id_));
//      gzread(backup_file, &sense_,  sizeof(sense_));
//      break;
//    }
//    default :
//    {
//      fprintf(stderr, "ERROR, invalid mutation type \"%d\" in file %s:%d.\n", mut_type_, __FILE__, __LINE__);
//      exit(EXIT_FAILURE);
//      break;
//    }
//  }
//
//
//}





// =================================================================
//                             Destructor
// =================================================================

// =================================================================
//                            Public Methods
// =================================================================

// =================================================================
//                           Protected Methods
// =================================================================
} // namespace aevol
