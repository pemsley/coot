/* src/coot_pythonmodule.cc
 * 
 * Copyright 2007 by The University of York
 * Author: Bernhard Lohkamp
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif
#include <iostream>

#if 1 // was #ifdef USE_PYGTK, we want a test for PyGObject, I suppose.

#include <gtk/gtk.h>

#include <pygobject-3.0/pygobject.h>

#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "coot-glue.hh"

PyObject * _wrap_main_menubar(PyObject *self);
PyObject * _wrap_main_statusbar(PyObject *self);
PyObject * _wrap_main_toolbar(PyObject *self);
PyObject * _wrap_main_hbox(PyObject *self);

// static/const?
PyMethodDef coot_python_functions[] = {

   { "main_menubar", (PyCFunction)_wrap_main_menubar, METH_NOARGS,
     NULL },
   
   { "main_statusbar", (PyCFunction)_wrap_main_statusbar, METH_NOARGS,
     NULL },

   { "main_toolbar", (PyCFunction)_wrap_main_toolbar, METH_NOARGS,
     NULL },

   { "main_hbox", (PyCFunction)_wrap_main_hbox, METH_NOARGS,
     NULL },

   { NULL, NULL, 0, NULL }
};

// void coot_python_register_classes (PyObject *d);
// extern PyMethodDef coot_python_functions[];



/* ---------- types from other modules ---------- */
// TMP
// static
PyTypeObject *_PyGObject_Type;

// #define PyGObject_Type (*_PyGObject_Type)

#include "coot-glue.hh"
// for glue:


/* ---------- forward type declarations ---------- */

/* #line 22 "coot_python.c" */



/* ----------- functions ----------- */

// try not static
// TMP
PyObject *
_wrap_main_menubar(PyObject *self)
{
   GtkWidget *ret;

   ret = main_menubar();
   
   /* pygobject_new handles NULL checking */
   return pygobject_new((GObject *)ret);
}

PyObject *
_wrap_main_statusbar(PyObject *self)
{
   GtkWidget *ret;
   
   ret = main_statusbar();
   
   /* pygobject_new handles NULL checking */
   return pygobject_new((GObject *)ret);
}

PyObject *
_wrap_main_toolbar(PyObject *self)
{
   GtkWidget *ret;
   
   
   ret = main_toolbar();
   
   /* pygobject_new handles NULL checking */
   return pygobject_new((GObject *)ret);
}

PyObject *
_wrap_main_hbox(PyObject *self)
{
   GtkWidget *ret;
   ret = main_hbox();
   /* pygobject_new handles NULL checking */
   return pygobject_new((GObject *)ret);
}


/* initialise stuff extension classes */
void
coot_python_register_classes(PyObject *d) {

    PyObject *module;
       
    if ((module = PyImport_ImportModule("gobject")) != NULL) {
        _PyGObject_Type = (PyTypeObject *)PyObject_GetAttrString(module, "GObject");
        if (_PyGObject_Type == NULL) {
            PyErr_SetString(PyExc_ImportError,
                "cannot import name GObject from gobject");
            return ;
        }
    } else {
        PyErr_SetString(PyExc_ImportError,
            "could not import gobject");
        return ;
    }
}


/* not sure what this is...
DL_EXPORT(void)
#ifdef WIN32
__declspec(dllexport)
#endif
*/

PyObject *some_test_function_py(PyObject *a, PyObject *b) {
   PyObject *o = Py_None;
   return o;
}

void
initcoot_python_gobject() {

   return;
   
   std::cout << "-------------------- in initcoot_python() " << std::endl;

   int req_major = -1, req_minor = -1, req_micro = -1;
   pygobject_init(req_major, req_minor, req_micro);

   PyMethodDef XYZMethods[] = {
	 { "some_test_function", some_test_function_py, METH_O, NULL},
	 { NULL, NULL, 0, NULL }
   };

   /* Python2
   PyObject *m = Py_InitModule3("coot_python", coot_python_functions);
   PyObject *d = PyModule_GetDict(m);
   coot_python_register_classes (d);
   */

   static struct PyModuleDef XYZ_module_e = {
    PyModuleDef_HEAD_INIT,
    "_xyz_coot",
    "Some Doc String Here",
    -1,
    XYZMethods,
    NULL, // reload
    NULL, // traverse
    NULL, // clear
    NULL  // free
   };

   struct PyModuleDef *XYZ_module = new struct PyModuleDef;
   *XYZ_module = XYZ_module_e;

   PyObject *m = PyModule_Create(XYZ_module);

   if (m) {

      PyObject *pMap = NULL;
      pMap = PyDict_New();
      Py_INCREF(pMap);
      int mao = PyModule_AddObject(m, "spam", pMap);
      if (mao < 0) {
         // Py_DECREF(pMap);
      } else {
         std::cout << "WARNING:: PyModule_AddObject() returned bad " << std::endl;
      }
   } else {
      std::cout << "WARNING:: PyModule_Create() returned null " << std::endl;
   }
   if (PyErr_Occurred()) {
      Py_FatalError ("can't initialise module _xyz_coot");
   }
}


#endif // USE_PYGTK
