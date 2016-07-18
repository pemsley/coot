
#ifdef USE_PYTHON
#include "Python.h"
#endif

#include<graphics-info.h>

/*  ----------------------------------------------------------------------- */
/*                  user-defined atom colours                               */
/*  ----------------------------------------------------------------------- */

#ifdef USE_PYTHON
void set_user_defined_atom_colour_by_residue_py(int imol, PyObject *residue_specs_colour_index_tuple_list_py);
void set_user_defined_atom_colour_py(int imol, PyObject *atom_specs_colour_index_tuple_list_py);
#endif // USE_PYTHON

#ifdef USE_GUILE
void set_user_defined_atom_colour_scm(int imol, SCM atom_specs_colour_index_tuple_list_scm);
void set_user_defined_atom_colour_by_residue_scm(int imol, SCM residue_specs_colour_index_tuple_list_scm);
#endif // USE_GUILE


