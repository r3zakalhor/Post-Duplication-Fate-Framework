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

#ifndef AEVOL_OBSERVER_H_
#define AEVOL_OBSERVER_H_


// ============================================================================
//                                   Includes
// ============================================================================
#include "ObservableEvent.h"

class Observable;


/**
 *
 */
class Observer {

 public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  Observer() = default; //< Default ctor
  Observer(const Observer&) = delete; //< Copy ctor
  Observer(Observer&&) = delete; //< Move ctor

  // ==========================================================================
  //                                Destructor
  // ==========================================================================

  virtual ~Observer() = default; //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================
  Observer& operator=(const Observer& other); //< Copy assignment
  Observer& operator=(const Observer&& other); //< Move assignment

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  /// This method is called whenever the observed object is changed.
  virtual void update(Observable& o, ObservableEvent e, void* arg) = 0;

  // ==========================================================================
  //                                 Getters
  // ==========================================================================

  // ==========================================================================
  //                                 Setters
  // ==========================================================================

 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================

  // ==========================================================================
  //                               Attributes
  // ==========================================================================
};


#endif // AEVOL_OBSERVER_H_
