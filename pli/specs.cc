
#ifdef USE_PYTHON

#include "specs.hh"
#include<iostream>

coot::py_atom_spec_t::py_atom_spec_t(PyObject *obj) {

   if (PyList_Check(obj)) {
      int len_view = PyList_Size(obj);
      if (len_view == 5) {
      	 PyObject *chain_id_py = PyList_GetItem(obj, 0);
	 chain_id = PyString_AsString(chain_id_py);
	 PyObject *resno_python = PyList_GetItem(obj, 1);
	 res_no = PyInt_AsLong(resno_python);
	 PyObject *ins_code_py = PyList_GetItem(obj, 2);
	 ins_code = PyString_AsString(ins_code_py);
	 PyObject *atom_name_py = PyList_GetItem(obj, 3);
	 atom_name = PyString_AsString(atom_name_py);
	 PyObject *alt_conf_py = PyList_GetItem(obj, 4);
	 alt_conf = PyString_AsString(alt_conf_py);
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
   PyList_SetItem(r, 0, PyString_FromString(chain_id.c_str()));
   PyList_SetItem(r, 1, PyInt_FromLong(res_no));
   PyList_SetItem(r, 2, PyString_FromString(ins_code.c_str()));
   PyList_SetItem(r, 3, PyString_FromString(atom_name.c_str()));
   PyList_SetItem(r, 4, PyString_FromString(alt_conf.c_str()));
   return r;
}


coot::py_residue_spec_t::py_residue_spec_t(PyObject *obj) {

   if (PyList_Check(obj)) {
      int len = PyList_Size(obj);
      int offset = 0;
      if (len == 4)
	 int offset = 1;

      if (len > 2) {
	 PyObject *chain_id_py = PyList_GetItem(obj,  0+offset);
	 chain_id = PyString_AsString(chain_id_py);
	 PyObject *resno_python = PyList_GetItem(obj, 1+offset);
	 res_no = PyInt_AsLong(resno_python);
	 PyObject *ins_code_py = PyList_GetItem(obj,  2+offset);
	 ins_code = PyString_AsString(ins_code_py);
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
   PyList_SetItem(r, 0, PyString_FromString(chain_id.c_str()));
   PyList_SetItem(r, 1, PyInt_FromLong(res_no));
   PyList_SetItem(r, 2, PyString_FromString(ins_code.c_str()));
   return r;
}

#endif //USE_PYTHON
