
%module coot
%include "std_string.i"
%include "std_vector.i"
%include "std_pair.i"



%{

#include <cstdio>
#include "globjects.h"  //includes gtk/gtk.h
#include "coot-coord-utils.hh"
#include "c-interface.h"
#include "cc-interface.hh"
#include "c-interface-database.hh"
#include "c-interface-python.hh"
#include "c-interface-ligands-swig.hh"
#include "c-interface-mogul.hh"
#include "c-interface-sequence.hh"
#include "manipulation-modes.hh"
#include "rotamer-search-modes.hh"
#include "lbg-interface.hh"
#include "sdf-interface.hh"
#include "probe-clash-score.hh"
%}


%template(vector_string) std::vector<std::string>;
%template(vector_atom_spec) std::vector<coot::atom_spec_t>;
%template(pairbf) std::pair<bool, float>;


#include "globjects.h"  //includes gtk/gtk.h
#include "coot-coord-utils.hh"
/* actually we should ignore all GtkWidgets or similar.... */
%ignore main_window();
%ignore main_menubar();
%ignore main_statusbar();
%ignore main_toolbar();
/* conflicts with redefinition */
%ignore list_nomenclature_errors(int);

%include "c-interface.h"
%include "cc-interface.hh"
%include "c-interface-database.hh"
%include "c-interface-python.hh"
%include "c-interface-ligands-swig.hh"
%include "c-interface-mogul.hh"
%include "c-interface-sequence.hh"
%include "manipulation-modes.hh"
%include "rotamer-search-modes.hh"
%include "lbg-interface.hh"
%include "sdf-interface.hh"
// %include "coot-coord-utils.hh"
%include "probe-clash-score.hh"
