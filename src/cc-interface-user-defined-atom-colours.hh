
#ifdef USE_PYTHON
#include "Python.h"
#endif

#ifdef USE_GUILE
#include <libguile.h>
#endif

// #include<graphics-info.h>

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

//! \brief by default user-defined colour indices will use a colour wheel
//  but here we can set them ourselves by passing something like: [1, [0.2, 0.3, 0.4]], [2, [0.6, 0.4, 0.2]]]
//  i.e. a colour for each colour index. If the colour index goes beyond the limit, or there is a
//  missing value, then colour will be [0.5, 0.5, 0.5]
//
#ifdef USE_PYTHON
void set_user_defined_colours_py(PyObject *colour_list_py);
#endif

void clear_user_defined_atom_colours(int imol);

