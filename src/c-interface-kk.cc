/* src/c-interface-kk.cc
 * 
 * Copyright 2009 by Kevin Keating
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


#include <clipper/core/coords.h>
#include <clipper/core/xmap.h>
#include "peak-search.hh"


#include "graphics-info.h"
#include "c-interface.h"
  
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
    }
}
