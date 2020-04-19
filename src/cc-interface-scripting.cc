
#ifdef USE_PYTHON
#include <Python.h>
#include "c-interface-python.hh"
#endif

#include <iostream>

// these are needed to comiple the cc-interface-scripting.hh
#include <string>
#include <vector>

#ifdef USE_GUILE
#include <cstdio> /* for std::FILE in gmp.h for libguile.h */
#include <libguile.h>
#endif

#include "graphics-info.h"
#include "cc-interface-scripting.hh"
#include "key-bindings.hh"

// put other scripting functions in here

#ifdef USE_PYTHON
// Prefered:
void add_key_binding_gtk3_py(PyObject *key_py, int ctrl_key, PyObject *func, const std::string &description) {

   int key = 0;

   if (PyLong_Check(key_py)) {
      key = PyLong_AsLong(key_py);
   } else {
      if (PyUnicode_Check(key_py)) {
         std::string key_str = PyBytes_AS_STRING(PyUnicode_AsUTF8String(key_py));
         // now convert the string to a key code.
         if (!key_str.empty()) {
            guint kk = gdk_unicode_to_keyval(key_str[0]);
            std::cout << "debug here with kk " << kk << std::endl;
            key = kk;
         } else {
            std::cout << "WARNING:: empty key " << std::endl;
         }
      }
   }
   if (key != 0) {
      keyboard_key_t k(key, ctrl_key);
      key_bindings_t kb(func, description);
      graphics_info_t::add_key_binding(k, kb);
   } else {
      // use display on this
      std::cout << "WARNGING:: add_key_binding_gtk3_py() failed to interpet "
                << PyBytes_AS_STRING(PyUnicode_AsUTF8String(display_python(key_py)))
                << std::endl;
   }
}
#endif


#ifdef USE_GUILE
void add_key_binding_gtk3_scm(int key, int ctrl_key, SCM thunk, const std::string &description) {

}
#endif


void reload_shaders() {

   std::cout << "Here in reload_shaders() " << std::endl;
   graphics_info_t g;
   g.screen_framebuffer.tear_down();
   g.blur_framebuffer.tear_down();
   g.init_shaders(); // regenerates screen and blur framebuffers also
   g.graphics_draw();

}

#ifdef USE_PYTHON
void set_light_position_py(int light_id, float x, float y, float z) {

   glm::vec4 pos(x,y,z,1.0);
   graphics_info_t::lights[light_id].position = pos;
   graphics_info_t g;
   g.graphics_draw();
   glFlush();
}
#endif
