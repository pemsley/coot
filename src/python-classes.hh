
#ifdef USE_PYTHON

void setup_python_classes();

PyMODINIT_FUNC init_pathology_data();

// static
PyObject *PathologyData_new(PyTypeObject *type, PyObject *args, PyObject *kwds);

#endif
