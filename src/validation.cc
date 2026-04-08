/*
 * src/validation.cc
 *
 * Copyright 2018 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#include "validation.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "graphics-info.h"

#include "cc-interface.hh"  // for residue_spec_py
#include "coot-utils/c-beta-deviations.hh"


#ifdef USE_PYTHON
PyObject *c_beta_deviations_py(int imol) {

   PyObject *o = Py_False;

   if (graphics_info_t::is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      if (mol) {
	 std::map<mmdb::Residue *, std::map<std::string, coot::c_beta_deviation_t> >
	    residue_c_beta_map = coot::get_c_beta_deviations(mol);
	 std::map<mmdb::Residue *, std::map<std::string, coot::c_beta_deviation_t> >::const_iterator it;
	 // PyObject *dict_py = PyDict_New();
	 PyObject *outer_list_py = PyList_New(residue_c_beta_map.size());
	 unsigned int counter = 0;
	 for (it=residue_c_beta_map.begin(); it!=residue_c_beta_map.end(); ++it) {
	    // multiple alt confs
	    mmdb::Residue *res_key = it->first;
	    const std::map<std::string, coot::c_beta_deviation_t> &value_map = it->second;
	    std::map<std::string, coot::c_beta_deviation_t>::const_iterator it_inner;
	    PyObject *map_dict_py = PyDict_New();
	    for (it_inner=value_map.begin(); it_inner!=value_map.end(); ++it_inner) {
	       const std::string alt_conf_key = it_inner->first;
	       const coot::c_beta_deviation_t &cbd = it_inner->second;
	       PyDict_SetItemString(map_dict_py,
				    alt_conf_key.c_str(),
				    PyFloat_FromDouble(cbd.dist));
	    }
	    PyObject *residue_spec_py = residue_spec_to_py(coot::residue_spec_t(res_key));
	    PyObject *residue_pair_py = PyList_New(2);
	    PyList_SetItem(residue_pair_py, 0, residue_spec_py);
	    PyList_SetItem(residue_pair_py, 1, map_dict_py);
	    PyList_SetItem(outer_list_py, counter, residue_pair_py);
	    counter++;
	 }
	 o = outer_list_py;
      }
   }

   if (PyBool_Check(o))
      Py_INCREF(o);
   return o;
}
#endif // USE_PYTHON

#ifdef USE_GUILE
SCM c_beta_deviations_scm(int imol) {

   SCM r = SCM_BOOL_F;
   if (graphics_info_t::is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      if (mol) {
	 r = SCM_EOL;
	 std::map<mmdb::Residue *, std::map<std::string, coot::c_beta_deviation_t> >
	    residue_c_beta_map = coot::get_c_beta_deviations(mol);
	 std::map<mmdb::Residue *, std::map<std::string, coot::c_beta_deviation_t> >::const_iterator it;
	 SCM outer_list_scm = SCM_EOL;
	 unsigned int counter = 0;
	 for (it=residue_c_beta_map.begin(); it!=residue_c_beta_map.end(); it++) {
	    // multiple alt confs
	    mmdb::Residue *res_key = it->first;
	    const std::map<std::string, coot::c_beta_deviation_t> &value_map = it->second;
	    std::map<std::string, coot::c_beta_deviation_t>::const_iterator it_inner;
	    SCM map_dict_scm = SCM_EOL;
	    SCM residue_spec_scm = residue_spec_to_scm(coot::residue_spec_t(res_key));
	    for (it_inner=value_map.begin(); it_inner!=value_map.end(); it_inner++) {
	       const std::string alt_conf_key = it_inner->first;
	       const coot::c_beta_deviation_t &cbd = it_inner->second;
	       SCM item_scm = scm_list_2(scm_from_locale_string(alt_conf_key.c_str()),
					 scm_from_double(cbd.dist));
	       map_dict_scm = scm_cons(item_scm, map_dict_scm);
	    }
	    SCM l_scm = scm_list_2(residue_spec_scm, scm_reverse(map_dict_scm));
	    r = scm_cons(l_scm, r);
	 }
      }
   }
   return r;

}
#endif // USE_GUILE

void delete_unhappy_atom_markers() {

   for (unsigned int i=0; i<graphics_info_t::molecules.size(); i++) {
      if (! graphics_info_t::molecules[i].unhappy_atom_marker_positions.empty()) {
         graphics_info_t::molecules[i].unhappy_atom_marker_positions.clear();
         graphics_info_t::attach_buffers();
         const auto &positions = graphics_info_t::molecules[i].unhappy_atom_marker_positions;
         graphics_info_t::tmesh_for_unhappy_atom_markers.update_instancing_buffer_data(positions);
      }
   }
}

void remove_unhappy_atom_marker_py(int imol, PyObject *atom_spec_py) {
   // tricky lookup? Code needs reworking. Is it worth it?
   std::cout << "remove marker here for the atom spec " << std::endl;
}

void remove_all_unhappy_atom_markers() {
   graphics_info_t::remove_all_unhappy_atom_markers();
}

void add_unhappy_atom_marker(int imol, const coot::atom_spec_t &atom_spec) {

   if (graphics_info_t::is_valid_model_molecule(imol)) {
      graphics_info_t g;
      g.add_unhappy_atom_marker(imol, atom_spec);
   }
}

// atom_spec_list is a 5-member atom spec
void add_unhappy_atom_marker_py(int imol, PyObject *atom_spec_list_py) {

   if (graphics_info_t::is_valid_model_molecule(imol)) {
      if (PyList_Check(atom_spec_list_py)) {
         int n = PyObject_Length(atom_spec_list_py);
         if (n == 5) {
            PyObject *chain_id_py  = PyList_GetItem(atom_spec_list_py, 0);
            PyObject *res_no_py    = PyList_GetItem(atom_spec_list_py, 1);
            PyObject *ins_code_py  = PyList_GetItem(atom_spec_list_py, 2);
            PyObject *atom_name_py = PyList_GetItem(atom_spec_list_py, 3);
            PyObject *alt_conf_py  = PyList_GetItem(atom_spec_list_py, 4);
            if (PyLong_Check(res_no_py)) {
               if (PyUnicode_Check(chain_id_py)) {
                  if (PyUnicode_Check(ins_code_py)) {
                     if (PyUnicode_Check(atom_name_py)) {
                        if (PyUnicode_Check(alt_conf_py)) {
		           std::string chain_id  = PyBytes_AS_STRING(PyUnicode_AsEncodedString(chain_id_py, "UTF-8", "strict"));
		           std::string ins_code  = PyBytes_AS_STRING(PyUnicode_AsEncodedString(ins_code_py, "UTF-8", "strict"));
		           std::string atom_name = PyBytes_AS_STRING(PyUnicode_AsEncodedString(atom_name_py, "UTF-8", "strict"));
		           std::string alt_conf  = PyBytes_AS_STRING(PyUnicode_AsEncodedString(alt_conf_py, "UTF-8", "strict"));
                           long res_no = PyLong_AsLong(res_no_py);
                           coot::atom_spec_t atom_spec(chain_id, res_no, ins_code, atom_name, alt_conf);
                           add_unhappy_atom_marker(imol, atom_spec);
                        }
                     }
                  }
               }
            }
         }
      }
   }

}

