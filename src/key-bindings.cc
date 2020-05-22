
#ifdef USE_PYTHON
#include "Python.h"
#endif // USE_PYTHON

#include <iostream>
#include <string>

#include <gtk/gtk.h>

#include "key-bindings.hh"

void
key_bindings_t::run() const {
     if (type == BUILT_IN)
        func();
#ifdef USE_PYTHON
     if (type == PYTHON) {
        if (! scripting_function_text.empty())
           PyRun_SimpleString(scripting_function_text.c_str());

        // This is when we get a propr function - but (currently) if I pass that, then
        // the error messages disappear

        if (function_py) {
           PyObject *arg_list = PyTuple_New(0);
           PyObject *result_py = PyEval_CallObject(function_py, arg_list);
           if (result_py == NULL) {
              std::cout << "result_py was null" << std::endl;
              if (PyErr_Occurred())
                 PyErr_PrintEx(0);
           } else {
              const char *mess = "object: %s\n";
              PyObject *dest = PyUnicode_FromString(mess);
              PyObject *d_py = PyUnicode_Format(dest, result_py);
              if (PyUnicode_Check(d_py)) {
                 std::string s = PyBytes_AS_STRING(PyUnicode_AsUTF8String(d_py));
                 std::cout << s << std::endl;
              } else {
                 std::cout << "d_py was not unicode\n";
              }
           }
        }
     }

#endif
#ifdef USE_GUILE
     if (type == SCHEME) {
        SCM handler = scm_c_eval_string("(lambda (key . args) (display (list \"(key_bindings_t run()) Error in proc: key: \" key \" args: \" args)) (newline))");
        SCM v = scm_catch(SCM_BOOL_T, scm_thunk, handler);
     }
#endif

}

