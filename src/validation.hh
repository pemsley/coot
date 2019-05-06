
// Many more things should go or be moved to here

#ifdef USE_PYTHON
#include "Python.h"
#endif // USE_PYTHON

#ifdef USE_GUILE
#include <libguile.h>
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *c_beta_deviations_py(int imol);
#endif


#ifdef USE_GUILE
SCM c_beta_deviations_scm(int imol);
#endif


