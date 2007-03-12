
%module coot

%{

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
%}


#include "globjects.h"  //includes gtk/gtk.h
#include "coot-coord-utils.hh"
%include "c-interface.h"
%include "cc-interface.hh"

