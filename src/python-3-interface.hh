
#ifdef USE_PYTHON

#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems

PyObject *myPyString_FromString(const char *str);

char *myPyString_AsString(PyObject *r);

#endif // USE_PYTHON
