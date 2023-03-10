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

#ifndef AEVOL_OBSERVABLE_H_
#define AEVOL_OBSERVABLE_H_


// ============================================================================
//                                   Includes
// ============================================================================
#include "ObservableEvent.h"
#include "Observer.h"

#include <list>
#include <map>

using std::list;
using std::map;




/**
 *
 */
  class Observable {
   public :
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  Observable() = default; //< Default ctor
  Observable(const Observable&) = delete; //< Copy ctor
  Observable(Observable&&) = delete; //< Move ctor

  // ==========================================================================
  //                                Destructor
  // ==========================================================================

  virtual ~Observable() = default; //< Destructor

  // ==========================================================================
  //                                Operators
  // ==========================================================================
  Observable& operator=(const Observable& other); //< Copy assignment
  Observable& operator=(const Observable&& other); //< Move assignment

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  void addObserver(Observer* o, ObservableEvent e) {
    observers_[e].emplace_back(o);
/*    if (e == MUTATION)
        printf("Size ADD-OBS %ld\n",observers_[e].size());*/
  };

  void deleteObserver(Observer* o, ObservableEvent e) {
    observers_[e].remove(o);
  };

  void clearAllObserver() {
      observers_[NEW_INDIV].clear();
      observers_[MUTATION].clear();
      observers_[END_REPLICATION].clear();
      observers_[END_GENERATION].clear();
  };

  void notifyObservers(ObservableEvent e, void* arg = nullptr);

  // ==========================================================================
  //                                 Getters
  // ==========================================================================

  // ==========================================================================
  //                                 Setters
  // ==========================================================================
  map<ObservableEvent, list<Observer*> > observers_;
 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================

  // ==========================================================================
  //                               Attributes
  // ==========================================================================

};


#endif // AEVOL_OBSERVABLE_H_
