
#include "Python.h"
#include <string>
#include <iostream>
#include "utils/coot-utils.hh"

void layla_setup_python_basic(int argc, char **argv) {

#ifdef USE_PYTHON
#ifdef USE_PYMAC_INIT

  //  (on Mac OS, call PyMac_Initialize() instead)
  //http://www.python.org/doc/current/ext/embedding.html
  //
  PyMac_Initialize();

#else

   wchar_t** _argv = static_cast<wchar_t **>(PyMem_Malloc(sizeof(wchar_t*)*argc));
   for (int i=0; i<argc; i++) {
      wchar_t* arg = Py_DecodeLocale(argv[i], NULL);
      _argv[i] = arg;
   }
   Py_InitializeEx(0);
   PySys_SetArgv(argc, _argv);

   // We expect these to be null because we are outside a python script.
   PyObject *globals = PyEval_GetGlobals();
   // std::cout << "in setup_python_basic() globals " << globals << std::endl;
   PyObject *locals  = PyEval_GetLocals();
   // std::cout << "in setup_python_basic() locals " << locals << std::endl;
   

#endif // USE_PYMAC_INIT

   auto get_pythondir = [] () {
                           std::string p = coot::prefix_dir();
                           std::string dp   = coot::util::append_dir_dir(p,   "lib");
                           std::string python_version = "python";
                           python_version += coot::util::int_to_string(PY_MAJOR_VERSION);
                           python_version += ".";
                           python_version += coot::util::int_to_string(PY_MINOR_VERSION);
                           std::string ddp  = coot::util::append_dir_dir(dp,  python_version);
                           std::string dddp = coot::util::append_dir_dir(ddp, "site-packages");
                           return dddp;
                        };
   auto get_pkgpythondir = [get_pythondir] () {
                              std::string d = get_pythondir();
                              std::string dp   = coot::util::append_dir_dir(d, "coot");
                              return dp;
                           };

   // std::string pkgpydirectory = PKGPYTHONDIR;
   // std::string pydirectory = PYTHONDIR;
   // use ${prefix}/lib/python3.9/site-package for PYTHONDIR
   // use ${pythondir}/coot' for PKGPYTHONDIR (i.e. PYTHONDIR + "/coot")

   std::string pkgpydirectory = get_pkgpythondir();
   std::string    pydirectory = get_pythondir();

   if (true) {
      std::cout << "debug:: in setup_python()    pydirectory is " << pydirectory << std::endl;
      std::cout << "debug:: in setup_python() pkgpydirectory is " << pkgpydirectory << std::endl;
   }

   PyObject *sys_path = PySys_GetObject("path");
   PyList_Append(sys_path, PyUnicode_FromString(pydirectory.c_str()));
   PyList_Append(sys_path, PyUnicode_FromString(pkgpydirectory.c_str()));

   // int err = PyRun_SimpleString("import coot");

#endif // USE_PYTHON

}

void layla_setup_python_coot_module() {

   PyObject *coot = PyImport_ImportModule("coot");
   if (! coot) {
      std::cout << "ERROR:: setup_python_coot_module() Null coot" << std::endl;
   } else {
      if (false)
         std::cout << "INFO:: setup_python_coot_module() good coot module" << std::endl;
   }
}


void layla_setup_python_module(const std::string &module_name) {

   PyObject *coot = PyImport_ImportModule(module_name.c_str());
   if (! coot) {
      std::cout << "ERROR:: setup_python_coot_module() Null module for " << module_name << std::endl;
   } else {
      if (false)
         std::cout << "INFO:: setup_python_coot_module() good module " << module_name << std::endl;
   }
}
