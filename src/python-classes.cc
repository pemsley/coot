
#ifdef USE_PYTHON

#include <iostream>

#include <Python.h>
#include "structmember.h"

typedef struct {
   PyObject_HEAD
   /* Type-specific fields go here. */
   PyObject *invresolsq_max; // as in clipper
   PyObject *fp_list;
    
} PathologyData;


static void
PathologyData_dealloc(PathologyData* self)
{
    Py_XDECREF(self->invresolsq_max);
    Py_XDECREF(self->fp_list);
    self->ob_type->tp_free((PyObject*) self);
}



static PyObject *
PathologyData_get_invresolsq_max(PathologyData *self, void *closure)
{
    Py_INCREF(self->invresolsq_max);
    return self->invresolsq_max;
}



static int
PathologyData_set_invresolsq_max(PathologyData *self, PyObject *value, void *closure)
{
  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the invresolsq_max attribute");
    return -1;
  }
  
  Py_DECREF(self->invresolsq_max);
  Py_INCREF(value);
  self->invresolsq_max = value;    

  return 0;
}


static PyObject *
PathologyData_get_fp_list(PathologyData *self, void *closure)
{
    Py_INCREF(self->fp_list);
    return self->fp_list;
}



static int
PathologyData_set_fp_list(PathologyData *self, PyObject *value, void *closure)
{
  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the fp_list attribute");
    return -1;
  }
  
  Py_DECREF(self->fp_list);
  Py_INCREF(value);
  self->fp_list = value;    

  return 0;
}



static PyGetSetDef PathologyData_getseters[] = {
   { (char *) "invresolsq_max", 
     (getter)PathologyData_get_invresolsq_max,
     (setter)PathologyData_set_invresolsq_max,
     (char *) "invresolsq_max",
     NULL},
   {(char *)"fp_list", 
    (getter)PathologyData_get_fp_list,
    (setter)PathologyData_set_fp_list,
    (char *)"fp_list",
    NULL},
   {NULL}  /* Sentinel */
};



static PyMemberDef PathologyData_members[] = {
   {(char *)"invresolsq_max", T_OBJECT_EX, offsetof(PathologyData, invresolsq_max), 0,
    (char *)"invresolsq_max resolution"},
   {(char *)"fp_list", T_OBJECT_EX, offsetof(PathologyData, fp_list), 0,
    (char *)"fp_list"},
    {NULL}  /* Sentinel */
};


static PyMethodDef PathologyData_methods[] = {
    {NULL}  /* Sentinel */
};

// static - Python API documentatation has static here
PyObject *
PathologyData_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PathologyData *self;

    self = (PathologyData *)type->tp_alloc(type, 0);
    if (self != NULL) {
       self->invresolsq_max = PyString_FromString("");
       if (self->invresolsq_max == NULL) {
	     Py_DECREF(self);
	     return NULL;
       }
        
       self->fp_list = PyString_FromString("");
       if (self->fp_list == NULL) {
	     Py_DECREF(self);
	     return NULL;
       }
    }
    return (PyObject *)self;
}



static int
PathologyData_init(PathologyData *self, PyObject *args, PyObject *kwds)
{
    PyObject *invresolsq_max=NULL, *fp_list=NULL, *tmp;

    static char *kwlist[] = {(char *)"invresolsq_max", (char *) "fp_list", NULL};

    if (! PyArg_ParseTupleAndKeywords(args, kwds, "|OO", kwlist, 
                                      &invresolsq_max, &fp_list))
        return -1; 

    if (invresolsq_max) {
        tmp = self->invresolsq_max;
        Py_INCREF(invresolsq_max);
        self->invresolsq_max = invresolsq_max;
        Py_XDECREF(tmp);
    }

    if (fp_list) {
       tmp = self->fp_list;
       Py_INCREF(fp_list);
       self->fp_list = fp_list;
       Py_XDECREF(tmp);
    }

    return 0;
}


static PyTypeObject pathology_data_PathologyDataType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "pathology_data.PathologyData", /*tp_name*/
    sizeof(PathologyData),     /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)PathologyData_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/
    "pathology data objects",  /* tp_doc */
    0,		               /* tp_traverse */
    0,		               /* tp_clear */
    0,		               /* tp_richcompare */
    0,		               /* tp_weaklistoffset */
    0,		               /* tp_iter */
    0,		               /* tp_iternext */
    PathologyData_methods,     /* tp_methods */
    PathologyData_members,     /* tp_members */
    PathologyData_getseters,   /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                            /* tp_descr_set */
    0,                            /* tp_dictoffset */
    (initproc)PathologyData_init, /* tp_init */
    0,                            /* tp_alloc */
    PathologyData_new,            /* tp_new */    
};

static PyMethodDef pathology_data_methods[] = {
    {NULL}  /* Sentinel */
};



#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
init_pathology_data() {
   
    pathology_data_PathologyDataType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&pathology_data_PathologyDataType) < 0)
        return;

    PyObject *m = Py_InitModule3("pathology_data", pathology_data_methods,
				 "Example module.");

    Py_INCREF(&pathology_data_PathologyDataType);
    PyModule_AddObject(m, "PathologyData", (PyObject *)&pathology_data_PathologyDataType);
}

#endif // USE_PYTHON

