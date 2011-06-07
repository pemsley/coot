
%module coot
%include "std_string.i"
%include "std_vector.i"
%include "std_pair.i"

%template(vector_string) std::vector<std::string>;


%{

#include <cstdio>
/* #ifdef USE_PYTHON */
/* #define SWIG_init SWIG_python_init */
/* #endif */

/* #ifdef USE_GUILE */
/* #define SWIG_init SWIG_guile_init */
/* #endif */

#include "globjects.h"  //includes gtk/gtk.h
#include "coot-coord-utils.hh"
#include "c-interface.h"
#include "cc-interface.hh"
#include "c-interface-database.hh"
#include "c-interface-python.hh"
#include "manipulation-modes.hh"
#include "rotamer-search-modes.hh"
#include "lbg-interface.hh"
%}


#include "globjects.h"  //includes gtk/gtk.h
#include "coot-coord-utils.hh"
/* actually we should ignore all GtkWidgets or similar.... */
%ignore main_window();
%ignore main_menubar();
%ignore main_statusbar();
%ignore main_toolbar();
%include "c-interface.h"
%include "cc-interface.hh"
%include "c-interface-database.hh"
%include "c-interface-python.hh"
%include "manipulation-modes.hh"
%include "rotamer-search-modes.hh"
%include "lbg-interface.hh"
