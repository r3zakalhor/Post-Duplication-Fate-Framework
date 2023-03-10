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
//*****************************************************************************


#ifndef AEVOL_PHENOTYPIC_TARGET_R_H__
#define AEVOL_PHENOTYPIC_TARGET_R_H__


// ============================================================================
//                                   Includes
// ============================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "PhenotypicTarget.h"
#include "Protein_R.h"


namespace aevol {

// ============================================================================
//                          Class declarations
// ============================================================================
//class Protein_R;





class PhenotypicTarget_R : public PhenotypicTarget
{
  // Faut il le laisser ?
  //friend class PhenotypicTargetHandler;

 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  PhenotypicTarget_R(void); //< Default ctor
  PhenotypicTarget_R(int16_t id);
  PhenotypicTarget_R(const PhenotypicTarget_R&); //< Copy ctor
  PhenotypicTarget_R(PhenotypicTarget_R&&) = delete; //< Move ctor

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~PhenotypicTarget_R(void); //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  void save(gzFile backup_file) const;
  void load(gzFile backup_file);
  // ==========================================================================
  //                                 Getters
  // ==========================================================================
  int16_t get_id() const {
    return id_;
  }
  std::list<Protein_R*> signals() const {
    return signals_;
  }

// ==========================================================================
  //                                 Setters
  // ==========================================================================
  void set_signals(const std::list<Protein_R*>& signals_list) {
      signals_ = signals_list;
  }

  void clear_signals() {
    for (auto prot : signals_) {
      delete prot;
    }
    signals_.clear();
  }

 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================

  // ==========================================================================
  //                               Attributes
  // ==========================================================================
  // An identifier used to know the position of this Phenotypic target in PhenotypicTargetHandler
  int16_t id_;
  std::list<Protein_R*> signals_;
};


// ============================================================================
//                           Getters' definitions
// ============================================================================

// ============================================================================
//                           Setters' definitions
// ============================================================================

// ============================================================================
//                          Operators' definitions
// ============================================================================

// ============================================================================
//                       Inline functions' definition
// ============================================================================

} // namespace aevol

#endif // AEVOL_PHENOTYPIC_TARGET_R_H__
