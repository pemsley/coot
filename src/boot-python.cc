
#ifdef USE_PYTHON
#include "Python.h"  // before guile includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "boot-python.hh"

void start_command_line_python_maybe(char **argv) {

#ifdef USE_PYTHON
   
   // Py_Main(0, argv);
   Py_Main(0, 0); // fixme - how do I convert to argv to wchar_t - if at all?

   //  Skip initialization registration of signal handlers, useful
   //  when Python is embedded. Version 2.4 or later. Thanks Stuart
   //  McNicholas for letting me know about this.
   //
   // Question: Do we need to check that we are not using command
   // line python no graphics before we use this?
   //
   Py_InitializeEx(0);

#endif

}
