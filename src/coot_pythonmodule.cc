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

#include <gtk/gtk.h>

#ifdef USE_PYTHON
#include <pygobject-3.0/pygobject.h>
#endif

#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "coot-glue.hh"

PyObject * _wrap_main_menumodel(PyObject *self);
PyObject * _wrap_main_statusbar(PyObject *self);
PyObject * _wrap_main_toolbar(PyObject *self);
PyObject * _wrap_main_hbox(PyObject *self);
PyObject * _wrap_application(PyObject *self);

// static/const?
PyMethodDef coot_python_functions[] = {

   { "main_menubar",   (PyCFunction)_wrap_main_menumodel, METH_NOARGS, NULL },
   { "main_statusbar", (PyCFunction)_wrap_main_statusbar, METH_NOARGS, NULL },
   { "main_toolbar",   (PyCFunction)_wrap_main_toolbar,   METH_NOARGS, NULL },
   { "main_hbox",      (PyCFunction)_wrap_main_hbox,      METH_NOARGS, NULL },
   { "application",    (PyCFunction)_wrap_application,    METH_NOARGS, NULL },
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
_wrap_main_menumodel(PyObject *self)
{
   GMenuModel *mm = main_menumodel();
   /* pygobject_new handles NULL checking */
   GObject *o = G_OBJECT(mm);
   PyObject *pyo = pygobject_new(o);
   return pyo;

}

PyObject *
_wrap_main_statusbar(PyObject *self)
{
   GtkWidget *ret = main_statusbar();
   /* pygobject_new handles NULL checking */
   GObject *o = G_OBJECT(ret);
   PyObject *pyo = pygobject_new(o);
   return pyo;
}

PyObject *
_wrap_main_toolbar(PyObject *self)
{
   GtkWidget *ret = main_toolbar();
   /* pygobject_new handles NULL checking */
   // std::cout << "main_toolbar(): ret " << ret << std::endl;
   GObject *o = G_OBJECT(ret);
   // std::cout << "main_toolbar():   o " << o << std::endl;
   PyObject *pyo = pygobject_new(o);
   // std::cout << "main_toolbar(): pyo " << pyo << std::endl;
   return pyo;
}

PyObject *
_wrap_main_hbox(PyObject *self) {
   GtkWidget *ret = main_hbox();
   /* pygobject_new handles NULL checking */
   GObject *o = G_OBJECT(ret);
   PyObject *pyo = pygobject_new(o);
   return pyo;
}

#include "graphics-info.h"

PyObject *
_wrap_application(PyObject *self) {
   GtkApplication *app = graphics_info_t::application;
   if (! app) {
      std::cout << "WARNING: ion _wrap_application(), application is not initialized" << std::endl;
   }
   /* pygobject_new handles NULL checking */
   GObject *o = G_OBJECT(app);
   PyObject *pyo = pygobject_new(o);
   return pyo;
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


struct module_state {
    PyObject *error;
};

#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))

static PyObject *
error_out(PyObject *m) {
    struct module_state *st = GETSTATE(m);
    PyErr_SetString(st->error, "something bad happened");
    return NULL;
}

static PyMethodDef coot_gui_api_methods[] = {
    {"main_menumodel", (PyCFunction)_wrap_main_menumodel, METH_NOARGS, NULL},
    {"main_statusbar", (PyCFunction)_wrap_main_statusbar, METH_NOARGS, NULL},
    {"main_toolbar",   (PyCFunction)_wrap_main_toolbar,   METH_NOARGS, NULL},
    {"main_hbox",      (PyCFunction)_wrap_main_hbox,      METH_NOARGS, NULL},
    {"application",    (PyCFunction)_wrap_application,    METH_NOARGS, NULL},
    {"error_out",      (PyCFunction)error_out,            METH_NOARGS, NULL},
    {NULL, NULL, 0, NULL}
};

static int coot_gui_api_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int coot_gui_api_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "coot_gui_api",
        NULL,
        sizeof(struct module_state),
        coot_gui_api_methods,
        NULL,
        coot_gui_api_traverse,
        coot_gui_api_clear,
        NULL
};


PyObject *
PyInit_coot_gui_api(void) {

   PyObject *module = PyModule_Create(&moduledef);

   if (! module) {
      std::cout << "in PyInit_coot_gui_api() module null" << std::endl;
   }

   if (module == NULL)
      return NULL;
   struct module_state *st = GETSTATE(module);

   st->error = PyErr_NewException("coot_gui_api.Error", NULL, NULL);
   if (st->error == NULL) {
      Py_DECREF(module);
      std::cout << "in PyInit_coot_gui_api() st->error null" << std::endl;
      return NULL;
   }

   if (PyErr_Occurred())
      PyErr_PrintEx(0);

   // std::cout << "--------------------- in PyInit_coot_gui_api() returning module " << module << std::endl;
   return module;
}


void
initcoot_python_gobject() {

   int req_major = -1, req_minor = -1, req_micro = -1;
// #ifdef HAVE_PYGOBJECT
   pygobject_init(req_major, req_minor, req_micro);
// #endif

   if (true) {
      PyObject *o = PyInit_coot_gui_api();
      if (! o) {
         std::cout << "Null o" << std::endl;
         return;
      }

      // Insert this into sys.modules directly - thanks Nick!
      PyObject *sys = PyImport_ImportModule("sys");
      if (! sys) {
         std::cout << "Null sys" << std::endl;
         return;
      }
      PyObject *modules = PyObject_GetAttrString(sys, "modules");
      if (! modules) {
         std::cout << "Null modules" << std::endl;
         return;
      }
      PyDict_SetItemString(modules, "coot_gui_api", o);
      Py_DECREF(modules);
      Py_DECREF(sys);

   }
}


// #endif // USE_PYGTK
