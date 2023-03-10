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

#ifndef AEVOL_EXP_SETUP_X11_H_
#define AEVOL_EXP_SETUP_X11_H_

// =================================================================
//                              Libraries
// =================================================================
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <X11/Xlib.h>
#include <X11/keysym.h>


// =================================================================
//                            Project Files
// =================================================================
#include "ExpManager.h"
#include "X11Window.h"

namespace aevol {

// =================================================================
//                          Class declarations
// =================================================================
class ExpSetup;
class X11Window;

enum key_map {
  KEY_ESCAPE = 0,
  KEY_F1 = 1,
  KEY_F2 = 2,
  KEY_F3 = 3,
  KEY_F4 = 4,
  KEY_F5 = 5,
  KEY_F6 = 6,
  KEY_F7 = 7,
  KEY_F8 = 8,
  KEY_F9 = 9,
  KEY_F10 = 10,
  KEY_F11 = 11,
  KEY_F12 = 12,
  KEY_A = 13,
  KEY_Q = 14,
  KEY_W = 15,
  KEY_Z = 16,
  KEY_S = 17,
  KEY_X = 18,
  KEY_E = 19,
  KEY_D = 20,
  KEY_C = 21,
  KEY_R = 22,
  KEY_F = 23,
  KEY_V = 24,
  KEY_T = 25,
  KEY_G = 26,
  KEY_B = 27,
  KEY_Y = 28,
  KEY_H = 29,
  KEY_N = 30,
  KEY_U = 31,
  KEY_J = 32,
  KEY_I = 33,
  KEY_K = 34,
  KEY_O = 35,
  KEY_L = 36,
  KEY_P = 37,
  KEY_M = 38,
  KEY_1 = 41,
  KEY_2 = 42,
  KEY_3 = 43,
  KEY_4 = 44,
  KEY_5 = 45,
  KEY_6 = 46,
  KEY_7 = 47,
  KEY_8 = 48,
  KEY_9 = 49
};


class ExpManager_X11 : public ExpManager
{
  friend class ExpSetup;
  
  public :
    
    // =================================================================
    //                             Constructors
    // =================================================================
    ExpManager_X11(void);
  
    // =================================================================
    //                             Destructors
    // =================================================================
    virtual ~ExpManager_X11(void) noexcept;
  
    // =================================================================
    //                              Accessors
    // =================================================================
    inline bool display_on();
    inline Display* get_display();
    inline int8_t screen();
    inline Atom* atoms();
    bool show_window(int8_t win) { return static_cast<bool>((show_window_ >> win) & 1); }
    bool new_show_window(int8_t win) { return static_cast<bool>((new_show_window_ >> win) & 1); }
    inline X11Window * window(int8_t win);
  
    // =================================================================
    //                            Public Methods
    // =================================================================
    KeyCode* key_codes() { return key_codes_;  };
    virtual void display(void);
    void toggle_display_on_off(void);
    void handle_events(void);
    void display(X11Window * win, const AbstractFuzzy& fuzzy, color_map color,
        bool fill = false, bool bold = false);
    void display_3D(X11Window * win,
                                  const AbstractFuzzy& fuzz, color_map color,
                    int x0 , int y0, bool fill /*= false*/ );
    void display_grid(X11Window * win, double** cell_grid);

    // =================================================================
    //                           Public Attributes
    // =================================================================
  
  
  
  
  
    protected :
  
    // =================================================================
    //                         Forbidden Constructors
    // =================================================================
    //~ ExpManager_X11(void)
    //~ {
      //~ printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      //~ exit( EXIT_FAILURE );
    //~ };
    ExpManager_X11( const ExpManager_X11 &model )
    {
      printf( "ERROR : Call to forbidden constructor in file %s : l%d\n", __FILE__, __LINE__ );
      exit( EXIT_FAILURE );
    };
  
  // =================================================================
  //                           Protected Methods
  // =================================================================
  void initialize(bool with_grid = false, bool with_plasmids = false);
  void compute_colormap();
  void set_codes();
  int8_t identify_window(Window winID);
  void draw_window(int8_t win_number);
  void refresh_window(int8_t win_number);

  // =================================================================
  //                          Protected Attributes
  // =================================================================
  bool     display_on_;
  bool     handle_display_on_off_;
  uint32_t show_window_;     // (bitmap) windows that have to be displayed (user switches value pressing F1, F2, ...)
  uint32_t new_show_window_; // (bitmap) windows that have to be displayed but were not displayed at the last refresh
  Display* display_;
  int8_t   screen_;
  Atom*    atoms_;
  KeyCode* key_codes_;

  X11Window**    win_;      // Table containing the <nb_windows> windows
  char**         win_name_; // window names
  unsigned int** win_size_; // window sizes
  int**          win_pos_;  // window positions

  std::vector<char*> col_map_;
};

// =====================================================================
//                          Accessors' definitions
// =====================================================================

// =====================================================================
//                       Inline functions' definition
// =====================================================================

} // namespace aevol

#endif // AEVOL_EXP_SETUP_X11_H_
