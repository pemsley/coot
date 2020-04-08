
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
   std::string scripting_function_text;
#ifdef USE_GUILE
   SCM thunk;
#endif
#ifdef USE_PYTHON
   PyObject *function_py;
#endif
   std::string description;
   void (*func)();
   key_bindings_t() { type = NONE; function_py = 0; }
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
   key_bindings_t(const std::string function_as_string, const std::string &descr) {
      type = PYTHON;
      function_py = 0;
      scripting_function_text = function_as_string;
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
   void run() const;
};

#endif // KEY_BINDINGS_HH

