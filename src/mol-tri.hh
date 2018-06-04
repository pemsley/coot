
#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <string>

void
graphics_to_ribbon_representation(int imol);

void
graphics_to_balls_representation(int imol);

void
graphics_to_sticks_representation(int imol);

// maybe we can make some settings so that the surface is electrostatic surface
void
graphics_to_surface_representation(int imol);


#ifdef USE_PYTHON
// think of a better name for this
void
mol_tri_settings_py(int imol, const std::string &param, PyObject *o);

#endif


#ifdef USE_GUILE
void
mol_tri_settings_scm(int imol, const std::string &param, SCM o);

#endif
