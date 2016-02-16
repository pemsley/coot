
#ifdef USE_PYTHON
#include "Python.h"  // before guile includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "boot-python.hh"

void start_command_line_python_maybe(char **argv) {

#ifdef USE_PYTHON
   
#if PY_MAJOR_VERSION > 2 
   Py_Main(0, argv);
#else
#if PY_MINOR_VERSION > 2 
   Py_Main(0, argv);
#endif     // PY_MINOR_VERSION
#endif     // PY_MAJOR_VERSION

     //  Skip initialization registration of signal handlers, useful
     //  when Python is embedded. Version 2.4 or later. Thanks Stuart
     //  McNicholas for letting me know about this.
     //
     // Question: Do we need to check that we are not using command
     // line python no graphics before we use this?
     // 
#if PY_MAJOR_VERSION > 2
   Py_InitializeEx(0);
#endif     
#if PY_MAJOR_VERSION == 2
#if PY_MINOR_VERSION > 3
   Py_InitializeEx(0);
#endif
#endif
#endif
}
