
#include<iostream>
#include "specs.hh"

coot::py_atom_spec_t::py_atom_spec_t(PyObject *obj) {

   if (PyList_Check(obj)) {
      int len_view = PyList_Size(obj);
      if (len_view == 5) {
         char *errors = 0;
      	 PyObject *chain_id_py = PyList_GetItem(obj, 0);
         PyObject *s = PyUnicode_AsEncodedString(chain_id_py, "utf-8", errors);
         chain_id = PyBytes_AS_STRING(s);
	 PyObject *resno_python = PyList_GetItem(obj, 1);
	 res_no = PyLong_AsLong(resno_python);
	 PyObject *ins_code_py = PyList_GetItem(obj, 2);
	 ins_code = PyBytes_AS_STRING(PyUnicode_AsEncodedString(ins_code_py, "utf-8", errors));
	 PyObject *atom_name_py = PyList_GetItem(obj, 3);
	 atom_name = PyBytes_AS_STRING(PyUnicode_AsEncodedString(atom_name_py, "utf-8", errors));
	 PyObject *alt_conf_py = PyList_GetItem(obj, 4);
	 alt_conf = PyBytes_AS_STRING(PyUnicode_AsEncodedString(alt_conf_py, "utf-8", errors));
	 string_user_data = "OK";
      }
   }
}


// 5-item atom specs (atom_spec_to_py will give you a 6-item spec,
// that should be renamed to atom_spec_with_user_data_to_py)
//
PyObject *
coot::py_atom_spec_t::pyobject() const {

   PyObject *r = PyList_New(5);
   PyList_SetItem(r, 0, PyUnicode_FromString(chain_id.c_str()));
   PyList_SetItem(r, 1, PyLong_FromLong(res_no));
   PyList_SetItem(r, 2, PyUnicode_FromString(ins_code.c_str()));
   PyList_SetItem(r, 3, PyUnicode_FromString(atom_name.c_str()));
   PyList_SetItem(r, 4, PyUnicode_FromString(alt_conf.c_str()));
   return r;
}


coot::py_residue_spec_t::py_residue_spec_t(PyObject *obj) {

   if (PyList_Check(obj)) {
      int len = PyList_Size(obj);
      int offset = 0;
      if (len == 4)
	 int offset = 1;

      if (len > 2) {
         char *errors = 0;
	 PyObject *chain_id_py = PyList_GetItem(obj,  0+offset);
	 chain_id = PyBytes_AS_STRING(PyUnicode_AsEncodedString(chain_id_py, "utf-8", errors));
	 PyObject *resno_py = PyList_GetItem(obj, 1+offset);
	 res_no = PyLong_AsLong(resno_py);
	 PyObject *ins_code_py = PyList_GetItem(obj,  2+offset);
	 ins_code = PyBytes_AS_STRING(PyUnicode_AsEncodedString(ins_code_py, "utf-8", errors));
      }
   } else {
      std::cout << "WARNING:: oops in py_residue_spec_t() obj is not a list " << std::endl;
   }
}


// 5-item atom specs (atom_spec_to_py will give you a 6-item spec,
// that should be renamed to atom_spec_with_user_data_to_py)
//
PyObject *
coot::py_residue_spec_t::pyobject() const {
   PyObject *r = PyList_New(3);
   PyList_SetItem(r, 0, PyUnicode_FromString(chain_id.c_str()));
   PyList_SetItem(r, 1, PyLong_FromLong(res_no));
   PyList_SetItem(r, 2, PyUnicode_FromString(ins_code.c_str()));
#if 0
#endif
   return r;
}

