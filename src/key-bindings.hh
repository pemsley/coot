
#ifndef KEY_BINDINGS_HH
#define KEY_BINDINGS_HH

// Python header has already been loaded by now.

#ifdef USE_GUILE
#include <cstdio> /* for std::FILE in gmp.h for libguile.h */
#include <libguile.h>
#endif

class keyboard_key_t {
public:
   int gdk_key;
   bool ctrl_is_pressed;
   keyboard_key_t(int g) {
      gdk_key = g;
      ctrl_is_pressed = false;
   }
   keyboard_key_t(int g, bool state) {
      gdk_key = g;
      ctrl_is_pressed = state;
   }
   bool operator<(const keyboard_key_t &other) const {
      if (other.gdk_key < gdk_key) {
         return true;
      } else {
         if (other.gdk_key == gdk_key) {
            if (other.ctrl_is_pressed) {
               if (ctrl_is_pressed) {
                  return false;
               } else {
                  return true;
               }
            } else {
               if (ctrl_is_pressed) {
                  return false; // hmm.
               } else {
                  return false;
               }
            }
         } else {
            return false;
         }
      }
   }
};

class key_bindings_t {

public:
  enum binding_type { NONE, SCHEME, PYTHON, BUILT_IN };
  binding_type type;
  std::string scripting_function_text; // used?
#ifdef USE_GUILE
  SCM thunk;
#endif
#ifdef USE_PYTHON
  PyObject *function_py;
#endif
  std::string description;
  void (*func)();
  key_bindings_t() { type = NONE; }
  // external
#ifdef USE_GUILE
  key_bindings_t(SCM thunk_in, const std::string &descr) {
     type = SCHEME;
     thunk = thunk_in;
     description = descr;
  }
#endif
#ifdef USE_PYTHON
  key_bindings_t(PyObject *func, const std::string &descr) {
     type = PYTHON;
     function_py = func;
     description = descr;
  }
#endif
  key_bindings_t(binding_type bt, const std::string &fn, const std::string &description_in) :
    type(bt), scripting_function_text(fn), description(description_in) {}
  // internal
  key_bindings_t(void (*func_in) (), const std::string &description_in) {
    description = description_in;
    type = BUILT_IN;
    func = func_in;
  }
  void run() const {
     if (type == BUILT_IN)
        func();
#ifdef USE_PYTHON
     if (type == PYTHON) {
        // PyRun_SimpleString(scripting_function_text.c_str()); // no
        PyObject *arg_list = PyTuple_New(0);
        PyObject *result_py = PyEval_CallObject(function_py, arg_list);
     }
#endif
#ifdef USE_GUILE
     if (type == SCHEME) {
        SCM handler = scm_c_eval_string("(lambda (key . args) (display (list \"(key_bindings_t run()) Error in proc: key: \" key \" args: \" args)) (newline))");
        SCM v = scm_catch(SCM_BOOL_T, scm_thunk, handler);
     }
#endif

  }

};

#endif // KEY_BINDINGS_HH

