/* src/c-interface-python.cc
 * 
 * Copyright 2008 by The University of York
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
#define PYTHONH
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#include "python-3-interface.hh"
#endif

#include "compat/coot-sysdep.h"


#include <string>
#include <vector>

#include "utils/coot-utils.hh"

#ifdef USE_PYTHON
#include "c-interface-python.hh"

#include "graphics-info.h"

PyObject *myPyString_FromString(const char *str) {

   PyObject *r = PyUnicode_FromString(str);
   return r;
}

char *myPyString_AsString(PyObject *r) {

   // maybe I should check the type before doing this coversion?
   return PyBytes_AS_STRING(PyUnicode_AsUTF8String(r));

}


// This is a common denominator really.  It does not depend on mmdb,
// but it can't be declared in c-interface.h because then we'd have to
// include c-interface.h which would cause (resolvable, I think, not
// checked) problems.
// BL says:: not sure about this, guess/hope it's ok as is!?
// 
// return a python string, decode to c++ using PyUnicode_AsUTF8String

PyObject * display_python(PyObject *o) {

   const char *mess = "object: %s\n";
   PyObject *dest = myPyString_FromString(mess);
   return PyUnicode_Format(dest, o);
}


// e.g. ["B", 41, "", " CA ", ""]
std::pair<bool, coot::atom_spec_t>
make_atom_spec_py(PyObject *spec) {

   bool good_spec = 0;
   coot::atom_spec_t as;
   int spec_length = PyObject_Length(spec);

   if (spec_length == 5) {
      PyObject  *chain_id_py = PyList_GetItem(spec, 0);
      PyObject     *resno_py = PyList_GetItem(spec, 1);
      PyObject  *ins_code_py = PyList_GetItem(spec, 2);
      PyObject *atom_name_py = PyList_GetItem(spec, 3);
      PyObject  *alt_conf_py = PyList_GetItem(spec, 4);
      if (PyUnicode_Check(chain_id_py)  &&
	  PyUnicode_Check(ins_code_py)  &&
	  PyUnicode_Check(atom_name_py) &&
	  PyUnicode_Check(alt_conf_py)  &&
	  PyLong_Check(resno_py)) {
         std::string chain_id = PyBytes_AS_STRING(PyUnicode_AsUTF8String(chain_id_py));
	 int resno = PyLong_AsLong(resno_py);
	 std::string ins_code  = PyBytes_AS_STRING(PyUnicode_AsUTF8String(ins_code_py));
	 std::string atom_name = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_name_py));
	 std::string alt_conf  = PyBytes_AS_STRING(PyUnicode_AsUTF8String(alt_conf_py));
	 as = coot::atom_spec_t(chain_id, resno, ins_code, atom_name, alt_conf);
	 good_spec = 1;
      } else {
	 std::cout << "WARNING:: badly formated atom spec: "
		   << PyBytes_AS_STRING(PyUnicode_AsUTF8String(display_python(spec)))
		   << std::endl;
      } 
   }
   return std::pair<bool, coot::atom_spec_t> (good_spec, as);
}

std::pair<bool, coot::residue_spec_t>
make_residue_spec_py(PyObject *spec) {

   bool good_spec = 0;
   coot::residue_spec_t rs("A", 1);
   int spec_length = PyObject_Length(spec);
   // we can now allow specs that are of length 4.  specs of length
   // are created by het-groups (amongst other things) and have a
   // state in the first position, which we skip (using offset = 1).
   int offset = 0;
   if (spec_length == 4) offset = 1;
   if ((spec_length == 3) || (spec_length == 4)) {
      PyObject  *chain_id_py = PyList_GetItem(spec, 0+offset);
      PyObject     *resno_py = PyList_GetItem(spec, 1+offset);
      PyObject  *ins_code_py = PyList_GetItem(spec, 2+offset);
      std::string chain_id = PyBytes_AS_STRING(PyUnicode_AsUTF8String(chain_id_py));
      int resno = PyLong_AsLong(resno_py);
      std::string ins_code  = PyBytes_AS_STRING(PyUnicode_AsUTF8String(ins_code_py));
      rs = coot::residue_spec_t(chain_id, resno, ins_code);
      good_spec = 1;
   }
   return std::pair<bool, coot::residue_spec_t> (good_spec, rs);
}

// return -1 on sting/symbol not found
int key_sym_code_py(PyObject *po) {

   int r = -1;
   if (PyUnicode_Check(po)) { 
      std::string s = PyBytes_AS_STRING(PyUnicode_AsUTF8String(po));
      r = coot::util::decode_keysym(s);
   }
   return r;
}


clipper::Spacegroup
py_symop_strings_to_space_group(PyObject *symop_string_list) {

   clipper::Spacegroup sg;
   if (PyList_Check(symop_string_list)) {
      int n = PyObject_Length(symop_string_list);
      std::string sgo;
      for (int i=0; i<n; i++) {
	 std::string se = PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyList_GetItem(symop_string_list, i)));
	 sgo += se;
	 sgo += " ; ";
      }
      if (sgo.length() > 0) {
	 try {
	    sg.init(clipper::Spgr_descr(sgo, clipper::Spgr_descr::Symops));
	 } catch ( clipper::Message_base exc ) {
	    std::string mess = "Can't make spacegroup from ";
	    mess += sgo;
	    std::cout << "WARNING:: " << mess << std::endl;
	 }
      }
   }
   return sg;
} 

PyObject *
atom_spec_to_py(const coot::atom_spec_t &spec) {

   graphics_info_t g;
   return g.atom_spec_to_py(spec);
}


#endif // USE_PYTHON
