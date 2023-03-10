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
#include <stdio.h>
#include <assert.h>
#include <string>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>
#include <unordered_map>
#include <string>

// =================================================================
//                            Project Files
// =================================================================
#include "X11Window.h"

#include "ExpSetup.h"
namespace aevol {

// =================================================================
//                       Basic X11/Xlib notions
// =================================================================
//
// THE DISPLAY
//
// The major notion of using Xlib is the X display. This is a structure
// representing the connection we have open with a given X server. It
// hides a queue of messages coming from the server, and a queue of
// pending requests that our client intends to send to the server.
// In Xlib, this structure is named 'Display'. When we open a connection
// to an X server, the library returns a pointer to such a structure.
// Later, we supply this pointer to any Xlib function that should send
// messages to the X server or receive messages from this server.
//
//
// THE WINDOWS
//
// X11 relies on a hierarchical model of rectangular areas called "Windows".
//
// 1. Each Window can be included in another Window (its parent) and may include
//    other Windows (its children). Windows sharing the same owner are called
//    siblings.
// 2. The screen itself is a Window (the Root Window) that contains all Windows.
// 3. A window can be above or behind a sibling Window. The Window which is above
//    hides partly or completely the other one.
// 4. Any drawing made in a Window is automatically "cut", meaning that only the
//    part of the drawing which is inside the Window is drawn.
// 5. A Window can be hidden or displayed ("mapped"). The drawing instructions
//    made on an unmapped Window are ignored. By default, newly created windows
//    are not mapped on the screen - they are invisible. In order to make a
//    window visible, we must use the XMapWindow() function.
// 6. Each event (keyboard, mouse) is aimed at a specific Window.
// 7. A Window does not memorize its content. Each time it must be re-displayed,
//    it gets an Expose event, and the content must be redrawn as a response to
//    this event.
//
//
// THE GC (GRAPHICS CONTEXT)
//
// When we perform various drawing operations (graphics, text, etc), we may
// specify various options for controlling how the data will be drawn - what
// foreground and background colors to use, how line edges will be connected,
// what font to use when drawing some text, etc). In order to avoid the need
// to supply zillions of parameters to each drawing function, a graphical context
// structure, of type 'GC' is used. We set the various drawing options in this
// structure, and then pass a pointer to this structure to any drawing routines.
// This is rather handy, as we often needs to perform several drawing requests
// with the same options. Thus, we would initialize a graphical context, set the
// desired options, and pass this GC structure to all drawing functions.
// Allocating a new GC is done using the XCreateGC() function.
//    GC XCreateGC(Display *display, Drawable d, uint32_t valuemask,
//                  XGCValues *values)
// Since a graphics context has zillions of attributes, and since often we want
// to define only few of them, we need to be able to tell the XCreateGC() which
// attributes we want to set. This is what the "valuemask" variable is for.
// We then use the "values" variable to specify actual values for the attributes
// we defined in the "valuesmask". The rest of the attributes of this GC will
// be set to their default values. Once we created a graphics context, we can
// use it in drawing functions. We can also modify its parameters using various
// functions (e.g. XSetForeground to change the foreground color of the GC).
//
//
// THE EVENTS
//
// A structure of type 'XEvent' is used to pass events received from the X server.
// Xlib supports a large amount of event types. The XEvent structure contains the
// type of event received, as well as the data associated with the event (e.g.
// position on the screen where the event was generated, mouse button associated
// with the event, region of screen associated with a 'redraw' event, etc). The way
// to read the event's data depends on the event type. Thus, an XEvent structure
// contains a C language union of all possible event types (if you're not sure what
// C unions are, it is time to check your favourite C language manual...). Thus,
// we could have an XExpose event, an XButton event, an XMotion event, etc.
// After a program creates a window (or several windows), it should tell the X
// server what types of events it wishes to receive for this window. By default,
// no events are sent to the program. This is done for optimizing the server-to-client
// connection (i.e. why send a program (that might even be running at the other
// side of the globe) an event it is not interested in?). It may register for
// various mouse (also called "pointer") events, keyboard events, expose events, etc.
// In Xlib, we use the XSelectInput() function to register for events. This function
// accepts 3 parameters - the display structure, an ID of a window, and a mask of
// the event types it wishes to get.




//##############################################################################
//                                                                             #
//                             Class X11Window                             #
//                                                                             #
//##############################################################################

// =================================================================
//                    Definition of static attributes
// =================================================================

// =================================================================
//                             Constructors
// =================================================================
X11Window::X11Window()
{
}


X11Window::X11Window(Display* display, int8_t screen, Atom* atoms,
                              uint16_t x, uint16_t y, uint16_t width, uint16_t height,
                              const char* caption)
{
  width_    = width;
  height_   = height;
  display_  = display;
  screen_   = screen;

  XSetWindowAttributes win_attributes;
  win_attributes.event_mask = StructureNotifyMask | ExposureMask | KeyPressMask;
  win_attributes.background_pixel = XBlackPixel(display_, screen_);

  window_ = XCreateWindow(display_, DefaultRootWindow(display_), x, y, width_, height_, 0,
                            CopyFromParent, CopyFromParent, CopyFromParent,
                            CWBackPixel|CWEventMask, &win_attributes);
  // NB: the 7th parameter is the width of the window's border, it has nothing to do with
  // the border appended by the window manager, so this is most often set to zero.

  // Define the title & iconname of the window
  XSetStandardProperties(display_, window_, caption, caption, None, NULL, 0, NULL);


  // We want to get MapNotify events, KeyPress events...
  XSelectInput(display_, window_, StructureNotifyMask | ExposureMask | KeyPressMask);


  // Create graphical contexts
  uint32_t whiteColor = WhitePixel(display_, screen_);
  XGCValues values;
  values.line_width = 1;

  values.foreground = pixel(display_, screen_, (char*)"white", whiteColor);
  values.background = pixel(display_, screen_, (char*)"black", whiteColor);
  gcWhite_ = XCreateGC(display_, window_, GCForeground|GCBackground|GCLineWidth, &values);


  values.foreground = pixel(display_, screen_, (char*)"black",whiteColor);
  values.background = pixel(display_, screen_, (char*)"white",whiteColor);
  gcBlack_ = XCreateGC(display_, window_, GCForeground|GCBackground|GCLineWidth, &values);

  values.foreground = pixel(display_, screen_, (char*)"red",  whiteColor);
  values.background = pixel(display_, screen_, (char*)"black",whiteColor);
  gcRed_ = XCreateGC(display_, window_, GCForeground|GCBackground|GCLineWidth, &values);

  values.foreground = pixel(display_, screen_, (char*)"green",whiteColor);
  values.background = pixel(display_, screen_, (char*)"green",whiteColor);
  gcGreen_ = XCreateGC(display_, window_, GCForeground|GCBackground|GCLineWidth, &values);

  values.foreground = pixel(display_, screen_, (char*)"blue", whiteColor);
  values.background = pixel(display_, screen_, (char*)"black",whiteColor);
  gcBlue_ = XCreateGC(display_, window_, GCForeground|GCBackground|GCLineWidth, &values);

  values.foreground = pixel(display_, screen_, (char*)"orange",whiteColor);
  values.background = pixel(display_, screen_, (char*)"orange",whiteColor);
  gcOrange_ = XCreateGC(display_, window_, GCForeground|GCBackground|GCLineWidth, &values);

  values.foreground = pixel(display_, screen_, (char*)"yellow",whiteColor);
  values.background = pixel(display_, screen_, (char*)"yellow",whiteColor);
  gcYellow_ = XCreateGC(display_, window_, GCForeground|GCBackground|GCLineWidth, &values);

  values.foreground = pixel(display_, screen_, (char*)"lightgrey",whiteColor);
  values.background = pixel(display_, screen_, (char*)"lightgrey",whiteColor);
  gcLightGrey_ = XCreateGC(display_, window_, GCForeground|GCBackground|GCLineWidth, &values);

  values.foreground = pixel(display_, screen_, (char*)"darkgrey", whiteColor);
  values.background = pixel(display_, screen_, (char*)"darkgrey", whiteColor);
  gcDarkGrey_ = XCreateGC(display_, window_, GCForeground|GCBackground|GCLineWidth, &values);

  values.foreground = pixel(display_, screen_, (char*)"grey15", whiteColor);
  values.background = pixel(display_, screen_, (char*)"grey15", whiteColor);
  gcDarkerGrey_ = XCreateGC(display_, window_, GCForeground|GCBackground|GCLineWidth, &values);

  values.foreground = pixel(display_, screen_, (char*)"grey",whiteColor);
  values.background = pixel(display_, screen_, (char*)"grey",whiteColor);
  gcGrey_ = XCreateGC(display_, window_, GCForeground|GCBackground|GCLineWidth, &values);


  XMapWindow(display_, window_);
  XMoveWindow(display_, window_, x, y);
  XFlush(display_);


  // Necessary to handle window closing
  XSetWMProtocols(display_, window_, atoms, 2);
}


// =================================================================
//                             Destructors
// =================================================================

X11Window::~X11Window()
{
  XFreeGC(display_, gcWhite_);
  XFreeGC(display_, gcBlack_);
  XFreeGC(display_, gcRed_);
  XFreeGC(display_, gcGreen_);
  XFreeGC(display_, gcBlue_);
  XFreeGC(display_, gcOrange_);
  XFreeGC(display_, gcYellow_);
  XFreeGC(display_, gcGrey_);
  XFreeGC(display_, gcLightGrey_);
  XFreeGC(display_, gcDarkGrey_);
  XFreeGC(display_, gcDarkerGrey_);

  XDestroyWindow(display_, window_);
}

// =================================================================
//                            Public Methods
// =================================================================



void X11Window::resize(unsigned int width, unsigned int height)
{
  width_  = width;
  height_ = height;
}


void X11Window::draw_string(int16_t x, int16_t y, char * str)
{
  XDrawImageString(display_, window_, gcWhite_, x, y, str, strlen(str));
}

void X11Window::draw_line(int16_t x1, int16_t y1, int16_t x2, int16_t y2, color_map color, bool bold /*= false*/)
{
  GC* gc = NULL;

  // Determine which GC to use
  switch (color)
  {
    case WHITE :
      gc = &gcWhite_;
      break;
    case BLACK :
      gc = & gcBlack_;
      break;
    case RED :
      gc = & gcRed_;
      break;
    case GREEN :
      gc = & gcGreen_;
      break;
    case BLUE :
      gc = & gcBlue_;
      break;
    case ORANGE :
      gc = & gcOrange_;
      break;
    case YELLOW :
      gc = & gcYellow_;
      break;
    case GREY :
      gc = & gcGrey_;
      break;
    case LIGHT_GREY :
      gc = & gcLightGrey_;
      break;
    case DARK_GREY :
      gc = & gcDarkGrey_;
      break;
    case DARKER_GREY :
      gc = & gcDarkerGrey_;
      break;
  }


  // Draw line (lines if bold)
  XDrawLine(display_, window_, *gc, x1, y1, x2, y2);
  if (bold)
  {
    XDrawLine(display_, window_, *gc, x1-1, y1, x2-1, y2);
    XDrawLine(display_, window_, *gc, x1+1, y1, x2+1, y2);
  }
}

void X11Window::draw_line(int16_t x1, int16_t y1, int16_t x2, int16_t y2, char* color, bool bold /*= false*/)
{
  // Create custom GC
  XGCValues values;
  values.foreground = pixel(display_, screen_, color, WhitePixel(display_,screen_));
  values.background = pixel(display_, screen_, color, WhitePixel(display_,screen_));
  GC tmp_gc = XCreateGC(display_, window_, GCForeground|GCBackground, &values);

  // Draw line (lines if bold)
  XDrawLine(display_, window_, tmp_gc, x1, y1, x2, y2);
  if (bold)
  {
    XDrawLine(display_, window_, tmp_gc, x1-1, y1, x2-1, y2);
    XDrawLine(display_, window_, tmp_gc, x1+1, y1, x2+1, y2);
  }

  XFreeGC(display_, tmp_gc);
}

void X11Window::draw_circle(int16_t x, int16_t y, int16_t diam)
{
  XDrawArc(display_, window_, gcWhite_, x, y, diam, diam, 0, 64*360);
}

void X11Window::draw_arc(int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2)
{
  XDrawArc(display_, window_, gcWhite_, x, y, diam, diam, 64*angle1, 64*angle2);
}

void X11Window::draw_arc(int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2, char* color)
{
  XGCValues values;
  values.line_width = 2;
  values.foreground = pixel(display_, screen_, color, WhitePixel(display_,screen_));
  values.background = pixel(display_, screen_, color, WhitePixel(display_,screen_));
  GC tmp_gc = XCreateGC(display_, window_, GCForeground|GCBackground|GCLineWidth, &values);

  XDrawArc(display_, window_, tmp_gc, x, y, diam, diam, 64*angle1, 64*angle2);

  XFreeGC(display_, tmp_gc);
}

void X11Window::draw_arc_64(int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2)
{
  XDrawArc(display_, window_, gcWhite_, x, y, diam, diam, angle1, angle2);
}

void X11Window::draw_arc_64(int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2, char* color)
{
  XGCValues values;
  values.line_width = 2;
  values.foreground = pixel(display_, screen_, color, WhitePixel(display_,screen_));
  values.background = pixel(display_, screen_, color, WhitePixel(display_,screen_));
  GC tmp_gc = XCreateGC(display_, window_, GCForeground|GCBackground|GCLineWidth, &values);

  XDrawArc(display_, window_, tmp_gc, x, y, diam, diam, angle1, angle2);

  XFreeGC(display_, tmp_gc);
}

void X11Window::fill_arc(int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2)
{
  XFillArc(display_, window_, gcWhite_, x, y, diam, diam, 64*angle1, 64*angle2);
}

void X11Window::fill_arc(int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2, char* color)
{
  XGCValues values;
  values.line_width = 2;
  values.foreground = pixel(display_, screen_, color, WhitePixel(display_,screen_));
  values.background = pixel(display_, screen_, color, WhitePixel(display_,screen_));
  GC tmp_gc = XCreateGC(display_, window_, GCForeground|GCBackground|GCLineWidth, &values);

  XFillArc(display_, window_, tmp_gc, x, y, diam, diam, 64*angle1, 64*angle2);

  XFreeGC(display_, tmp_gc);
}

void X11Window::fill_arc_64(int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2)
{
  XFillArc(display_, window_, gcWhite_, x, y, diam, diam, angle1, angle2);
}

void X11Window::fill_arc_64(int16_t x, int16_t y, int16_t diam, int16_t angle1, int16_t angle2, char* color)
{
  XGCValues values;
  values.line_width = 2;
  values.foreground = pixel(display_, screen_, color, WhitePixel(display_,screen_));
  values.background = pixel(display_, screen_, color, WhitePixel(display_,screen_));
  GC tmp_gc = XCreateGC(display_, window_, GCForeground|GCBackground|GCLineWidth, &values);

  XFillArc(display_, window_, tmp_gc, x, y, diam, diam, angle1, angle2);

  XFreeGC(display_, tmp_gc);
}

void X11Window::fill_rectangle(int16_t x, int16_t y, int16_t width, int16_t height, color_map color)
{
  switch (color)
  {
    case WHITE :
      XFillRectangle(display_, window_, gcWhite_, x, y, width, height);
      break;
    case BLACK :
      XFillRectangle(display_, window_, gcBlack_, x, y, width, height);
      break;
    case RED :
      XFillRectangle(display_, window_, gcRed_, x, y, width, height);
      break;
    case GREEN :
      XFillRectangle(display_, window_, gcGreen_, x, y, width, height);
      break;
    case BLUE :
      XFillRectangle(display_, window_, gcBlue_, x, y, width, height);
      break;
    case ORANGE :
      XFillRectangle(display_, window_, gcOrange_, x, y, width, height);
      break;
    case YELLOW :
      XFillRectangle(display_, window_, gcYellow_, x, y, width, height);
      break;
    case GREY :
      XFillRectangle(display_, window_, gcGrey_, x, y, width, height);
      break;
    case LIGHT_GREY :
      XFillRectangle(display_, window_, gcLightGrey_, x, y, width, height);
      break;
    case DARK_GREY :
      XFillRectangle(display_, window_, gcDarkGrey_, x, y, width, height);
      break;
    case DARKER_GREY :
      XFillRectangle(display_, window_, gcDarkerGrey_, x, y, width, height);
      break;
  }
}

void X11Window::fill_rectangle(int16_t x, int16_t y, int16_t width, int16_t height, char* color)
{
  XGCValues values;
  values.foreground = pixel(display_, screen_, color, WhitePixel(display_,screen_));
  values.background = pixel(display_, screen_, color, WhitePixel(display_,screen_));
  GC tmp_gc = XCreateGC(display_, window_, GCForeground|GCBackground, &values);

  XFillRectangle(display_, window_, tmp_gc, x, y, width, height);

  XFreeGC(display_, tmp_gc);
}

char*X11Window::color(double mean)
{
  int16_t red, green, blue;

  double  mean_range     = X_MAX - X_MIN;
  double  mean_range_5   = X_MIN + mean_range / 5;
  double  mean_range_2_5 = X_MIN + 2 * mean_range / 5;
  double  mean_range_3_5 = X_MIN + 3 * mean_range / 5;
  double  mean_range_4_5 = X_MIN + 4 * mean_range / 5;

  if (mean < mean_range_5)
  {
    red   = 0;
    green = 255 * (1.0 - ((mean_range_5 - mean) / mean_range_5));
    blue  = 255;
  }
  else if (mean < mean_range_2_5)
  {
    red   = 0;
    green = 255;
    blue  = 255 * ((mean_range_2_5 - mean) / mean_range_5);
  }
  else if (mean < mean_range_3_5)
  {
    red   = 255 * (1.0 - ((mean_range_3_5 - mean) / mean_range_5));
    green = 255;
    blue  = 0;
  }
  else if (mean < mean_range_4_5)
  {
    red   = 255;
    green = 255 * ((mean_range_4_5 - mean) / mean_range_5);
    blue  = 0;
  }
  else
  {
    red   = 255;
    green = 0;
    blue  = 255 * (1.0 - ((mean_range - mean) / mean_range_5));
  }

  char* color = new char[8];
  sprintf(color, "#%02x%02x%02x", red, green, blue);
  return color;
}


// =================================================================
//                           Protected Methods
// =================================================================

uint32_t X11Window::pixel(Display *display, int8_t screen, char *color_name, uint32_t default_color)
{
  // hacked on 2014-12-05 because XQuarz and Yosemite make XAllocColor veeery slow
  // dirty memoization workaround: display and screen are assumed to be constant
  // costs 21Kb of RAM on basic example (according to the measure given before the return)
  static std::unordered_map <std::string, unsigned long> color_memo;

  // if color_name is already recorded, compute and record it
  if (color_memo.find(color_name) == color_memo.end()) {
    XColor color;
    if (XParseColor(display, DefaultColormap(display,screen), color_name, &color) == 0) {
      fprintf(stderr, "Invalid color: %s\n", color_name);
      return default_color;
    }
    if (XAllocColor(display, DefaultColormap(display,screen), &color) == 0) {
      fprintf(stderr, "Could not allocate color %s\n", color_name);
      return default_color;
    }
    color_memo[color_name] = color.pixel;
  }
  // printf("ram used: %lu bytes\n", sizeof(color_memo) + color_memo.size() * (sizeof(std::string) + sizeof(unsigned long)));
  return color_memo[color_name];
}
} // namespace aevol
