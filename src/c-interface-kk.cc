/* src/c-interface-kk.cc
 * 
 * Copyright 2009, 2010, 2011 by Kevin Keating
 * Copyright 2009 The University of Oxford
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
#include <Python.h>  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <clipper/core/coords.h>
#include <clipper/core/xmap.h>
#include "peak-search.hh"


#include "graphics-info.h"
#include "c-interface.h"
#include "c-interface-python.hh"
  
#ifdef USE_PYTHON
PyObject *map_peaks_near_point_from_list_py(int imol_map, PyObject *peak_list, float x, float y, float z,
					    float radius) {
   PyObject *r = Py_False;
     
   if (is_valid_map_molecule(imol_map)) {
      // convert *peak_list into a C++ vector
      std::vector<clipper::Coord_orth> peaks;
      int l = PyObject_Length(peak_list);
      for (unsigned int peak_num=0; peak_num<l; peak_num++) {
	 PyObject *cur_peak_py = PyList_GetItem(peak_list, peak_num);
             
	 double peak_x = PyFloat_AsDouble(PyList_GetItem(cur_peak_py, 0));
	 double peak_y = PyFloat_AsDouble(PyList_GetItem(cur_peak_py, 1));
	 double peak_z = PyFloat_AsDouble(PyList_GetItem(cur_peak_py, 2));
	 clipper::Coord_orth cur_peak (peak_x, peak_y, peak_z);
            
	 peaks.push_back(cur_peak);
      }
         
      // find all peaks close to (x,y,z)
      const clipper::Xmap<float> &xmap = graphics_info_t::molecules[imol_map].xmap_list[0];
      clipper::Coord_orth center (x,y,z);
      coot::peak_search ps(xmap);
      std::vector<clipper::Coord_orth> close_peaks = 
	 ps.get_peaks_from_list(xmap, center, radius, peaks);
         
      //convert close_peaks to a Python list
      r = PyList_New(close_peaks.size());
      for (unsigned int i=0; i<close_peaks.size(); i++) {
	 PyObject *coords = PyList_New(3);
	 PyList_SetItem(coords, 0, PyFloat_FromDouble(close_peaks[i].x()));
	 PyList_SetItem(coords, 1, PyFloat_FromDouble(close_peaks[i].y()));
	 PyList_SetItem(coords, 2, PyFloat_FromDouble(close_peaks[i].z()));
	 PyList_SetItem(r, i, coords);
      }
   }
     
   if (PyBool_Check(r)) {
      Py_INCREF(r);
   }
 
   return r;
}
#endif 
 

#ifdef USE_PYTHON
PyObject *screen_vectors_py() {
  PyObject *vecs = PyList_New(3);
  
  PyObject *x_vecs = PyList_New(3);
  PyObject *y_vecs = PyList_New(3);
  PyObject *z_vecs = PyList_New(3);
  
  //create a new ScreenVectors object
  coot::ScreenVectors screen_vec_object;

  //copy the ScreenVectors object into a Python list
  PyList_SetItem(x_vecs, 0, PyFloat_FromDouble(screen_vec_object.screen_x.x()));
  PyList_SetItem(x_vecs, 1, PyFloat_FromDouble(screen_vec_object.screen_x.y()));
  PyList_SetItem(x_vecs, 2, PyFloat_FromDouble(screen_vec_object.screen_x.z()));

  PyList_SetItem(y_vecs, 0, PyFloat_FromDouble(screen_vec_object.screen_y.x()));
  PyList_SetItem(y_vecs, 1, PyFloat_FromDouble(screen_vec_object.screen_y.y()));
  PyList_SetItem(y_vecs, 2, PyFloat_FromDouble(screen_vec_object.screen_y.z()));

  PyList_SetItem(z_vecs, 0, PyFloat_FromDouble(screen_vec_object.screen_z.x()));
  PyList_SetItem(z_vecs, 1, PyFloat_FromDouble(screen_vec_object.screen_z.y()));
  PyList_SetItem(z_vecs, 2, PyFloat_FromDouble(screen_vec_object.screen_z.z()));

  PyList_SetItem(vecs, 0, x_vecs);
  PyList_SetItem(vecs, 1, y_vecs);
  PyList_SetItem(vecs, 2, z_vecs);
  
  return vecs;
}
#endif

void
clear_extra_restraints(int imol) {
    if (is_valid_model_molecule(imol)) {
        graphics_info_t::molecules[imol].extra_restraints.bond_restraints.clear();
        graphics_info_t::molecules[imol].extra_restraints.torsion_restraints.clear();
        graphics_info_t::molecules[imol].extra_restraints.start_pos_restraints.clear();
    }
}

#ifdef USE_PYTHON
int mark_multiple_atoms_as_fixed_py(int imol, PyObject *atom_spec_list, int state) {
   int r = 0;
   int list_length = PyObject_Length(atom_spec_list);
   PyObject *atom_spec;
   
   for (int ispec = 0; ispec<list_length; ispec++) {
      atom_spec = PyList_GetItem(atom_spec_list, ispec);
      std::pair<bool, coot::atom_spec_t> p = make_atom_spec_py(atom_spec);
      if (p.first) {
	 graphics_info_t::mark_atom_as_fixed(imol, p.second, state);
	 r++; //keep track of how many atoms we successfully marked
      }
   }
   
   if (r > 0) {
      graphics_draw();
   }
   
   return r; //return a count of how many atoms we successfully marked
}
#endif // USE_PYTHON

int add_extra_start_pos_restraint(int imol, const char *chain_id_1, int res_no_1, const char *ins_code_1, const char *atom_name_1, const char *alt_conf_1, double esd) {

   int r = -1;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t as_1(chain_id_1, res_no_1, ins_code_1, atom_name_1, alt_conf_1);
      r = graphics_info_t::molecules[imol].add_extra_start_pos_restraint(as_1, esd);
      //graphics_draw();
   }
   return r;
}

#ifdef USE_PYTHON
PyObject *refine_zone_with_score_py(int imol, const char *chain_id,
		 int resno1,
		 int resno2,
		 const char *altconf) {

   graphics_info_t g;
   PyObject *rv = Py_False;
   
   if (is_valid_model_molecule(imol)) {
      CResidue *res_1 = g.molecules[imol].get_residue(chain_id, resno1, "");
      CResidue *res_2 = g.molecules[imol].get_residue(chain_id, resno2, "");
      if (res_1 && res_2) { 
	 std::string resname_1(res_1->GetResName());
	 std::string resname_2(res_2->GetResName());
	 bool is_water_like_flag = g.check_for_no_restraints_object(resname_1, resname_2);
	 coot::refinement_results_t rr =
            g.refine_residue_range(imol, chain_id, chain_id, resno1, "", resno2, "", altconf,
				is_water_like_flag);
         rv = g.refinement_results_to_py(rr);
      }
   }
   
   if (PyBool_Check(rv)) {
     Py_INCREF(rv);
   }

   return rv;
}
#endif	/* USE_PYTHON */


#ifdef USE_GUILE
SCM refine_zone_with_score_scm(int imol, const char *chain_id,
		 int resno1,
		 int resno2,
		 const char *altconf) {

   graphics_info_t g;
   SCM rv = SCM_BOOL_F;
   
   if (is_valid_model_molecule(imol)) {
      CResidue *res_1 = g.molecules[imol].get_residue(chain_id, resno1, "");
      CResidue *res_2 = g.molecules[imol].get_residue(chain_id, resno2, "");
      if (res_1 && res_2) { 
	 std::string resname_1(res_1->GetResName());
	 std::string resname_2(res_2->GetResName());
	 bool is_water_like_flag = g.check_for_no_restraints_object(resname_1, resname_2);
         coot::refinement_results_t rr =
            g.refine_residue_range(imol, chain_id, chain_id, resno1, "", resno2, "", altconf,
				is_water_like_flag);
         rv = g.refinement_results_to_scm(rr);
      }
   }
   
   return rv;
}
#endif // USE_GUILE
