
#ifdef USE_PYGTK

#if defined COOT_USE_GTK2_INTERFACE || (GTK_MAJOR_VERSION == 2)

#include <pygobject.h>
#include <pygtk/pygtk.h>
#include <c-interface.h>

void coot_python_register_classes (PyObject *d);
extern PyMethodDef coot_python_functions[];

DL_EXPORT(void)
#ifdef WIN32
__declspec(dllexport)
#endif
initcoot_python(void)
{
    PyObject *m, *d;

    init_pygobject ();
    init_pygtk ();

    m = Py_InitModule ("coot_python", coot_python_functions);
    d = PyModule_GetDict (m);

    coot_python_register_classes (d);

    if (PyErr_Occurred ()) {
        Py_FatalError ("can't initialise module coot_python");
    }
}

#else   // use gtk-1.2

// Actually dont know if it will work like this yet!

#include <pygtk/pygtk.h>
#include <c-interface.h>

void coot_python_register_classes (PyObject *d);
extern PyMethodDef coot_python_functions[];

DL_EXPORT(void)
#ifdef WIN32
__declspec(dllexport)
#endif
initcoot_python(void)
{
    PyObject *m, *d;

    init_pygobject ();
    init_pygtk ();

    m = Py_InitModule ("coot_python", coot_python_functions);
    d = PyModule_GetDict (m);

    coot_python_register_classes (d);

    if (PyErr_Occurred ()) {
        Py_FatalError ("can't initialise module coot_python");
    }
}

#endif // GTK2

#endif // USE_PYGTK
