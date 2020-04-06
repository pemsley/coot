
#ifdef USE_PYTHON
#include "Python.h"  // before guile includes to stop "POSIX_C_SOURCE" redefined problems
#endif
#include <iostream>

#include "boot-python.hh"

// There is no maybe - rename.
//
void start_command_line_python_maybe(bool command_line_mode_flag, int argc, char **argv) {

#ifdef USE_PYTHON

   // we should only come here if --no-graphics --python was used
   //
   // for normal Python setup, we want to Py_InitializeEx(0) (so that Ctrl-C works)

   //  Skip initialization registration of signal handlers, useful
   //  when Python is embedded. Version 2.4 or later. Thanks Stuart
   //  McNicholas for letting me know about this.
   //
   // Question: Do we need to check that we are not using command
   // line python no graphics before we use this?
   //
   // I want Ctrl-C to exit the program, unless I am command line mode

   if (command_line_mode_flag)
      Py_InitializeEx(1);

   wchar_t** _argv = static_cast<wchar_t **>(PyMem_Malloc(sizeof(wchar_t*)));
   _argv[0] = Py_DecodeLocale(argv[0], NULL);

   Py_Main(1, _argv); // do I want command line argumesents? If so, see coot-setup-python.cc
                  // to see how to make wchar_** args

#endif

}
