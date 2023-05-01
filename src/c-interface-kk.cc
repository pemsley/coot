/* src/c-interface-kk.cc
 * 
 * Copyright 2009, 2010, 2011, 2012 by Kevin Keating
 * Copyright 2009 The University of Oxford
 * Copyright 2013 by Medical Research Council
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
#include "python-3-interface.hh"
#endif


#include <clipper/core/coords.h>
#include <clipper/core/xmap.h>
#include "coot-utils/peak-search.hh"


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
      for (int peak_num=0; peak_num<l; peak_num++) {
	 PyObject *cur_peak_py = PyList_GetItem(peak_list, peak_num);
             
	 double peak_x = PyFloat_AsDouble(PyList_GetItem(cur_peak_py, 0));
	 double peak_y = PyFloat_AsDouble(PyList_GetItem(cur_peak_py, 1));
	 double peak_z = PyFloat_AsDouble(PyList_GetItem(cur_peak_py, 2));
	 clipper::Coord_orth cur_peak (peak_x, peak_y, peak_z);
            
	 peaks.push_back(cur_peak);
      }
         
      // find all peaks close to (x,y,z)
      const clipper::Xmap<float> &xmap = graphics_info_t::molecules[imol_map].xmap;
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

   PyObject *vecs = Py_None;

   // we need graphics to create a ScreenVectors
   // 
   if (graphics_info_t::use_graphics_interface_flag) { 
   
      vecs = PyList_New(3);
  
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
   } else {
      Py_INCREF(vecs);
   }
  return vecs;
}
#endif

void
clear_extra_restraints(int imol) {
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].clear_extra_restraints(); 
      graphics_info_t::molecules[imol].set_display_extra_restraints(0);
   }
   graphics_draw();
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

int add_extra_angle_restraint(int imol, 
				const char *chain_id_1, int res_no_1, const char *ins_code_1, const char *atom_name_1, const char *alt_conf_1, 
				const char *chain_id_2, int res_no_2, const char *ins_code_2, const char *atom_name_2, const char *alt_conf_2, 
				const char *chain_id_3, int res_no_3, const char *ins_code_3, const char *atom_name_3, const char *alt_conf_3, 
				double angle, double esd) {

   int r = -1;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t as_1(chain_id_1, res_no_1, ins_code_1, atom_name_1, alt_conf_1);
      coot::atom_spec_t as_2(chain_id_2, res_no_2, ins_code_2, atom_name_2, alt_conf_2);
      coot::atom_spec_t as_3(chain_id_3, res_no_3, ins_code_3, atom_name_3, alt_conf_3);
      r = graphics_info_t::molecules[imol].add_extra_angle_restraint(as_1, as_2, as_3, angle, esd);
      graphics_draw();
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
      mmdb::Residue *res_1 = g.molecules[imol].get_residue(chain_id, resno1, "");
      mmdb::Residue *res_2 = g.molecules[imol].get_residue(chain_id, resno2, "");
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
      mmdb::Residue *res_1 = g.molecules[imol].get_residue(chain_id, resno1, "");
      mmdb::Residue *res_2 = g.molecules[imol].get_residue(chain_id, resno2, "");
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

#ifdef USE_PYTHON
PyObject *regularize_zone_with_score_py(int imol, const char *chain_id, int resno1, int resno2, const char *altconf) {
   int status = 0;
   
   PyObject *rv = Py_False;
   
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      // the "" is the insertion code (not passed to this function (yet)
      int index1 = graphics_info_t::molecules[imol].atom_index_first_atom_in_residue(chain_id, resno1, ""); 
      int index2 = graphics_info_t::molecules[imol].atom_index_first_atom_in_residue(chain_id, resno2, "");
      short int auto_range = 0;
      if (index1 >= 0) {
	 if (index2 >= 0) { 
	    coot::refinement_results_t rr = g.regularize(imol, auto_range, index1, index2);
	    std::cout << "debug:: restraints results " << rr.found_restraints_flag << " "
		      << rr.lights.size() << " " << rr.info_text << std::endl;
	    if ((rr.lights.size() > 0) || (rr.found_restraints_flag)) {
	       rv = g.refinement_results_to_py(rr);
	    }
	    
	 } else {
	    std::cout << "WARNING:: regularize_zone: Can't get index for resno2: "
		      << resno2 << std::endl;
	 } 
      } else {
	 std::cout << "WARNING:: regularize_zone: Can't get index for resno1: "
		   << resno1 << std::endl;
      }
   } else {
      std::cout << "Not a valid model molecule" << std::endl;
   }
   
   if (PyBool_Check(rv)) {
     Py_INCREF(rv);
   }

   return rv;
}
#endif	/* USE_PYTHON */

#ifdef USE_GUILE
SCM regularize_zone_with_score_scm(int imol, const char *chain_id, int resno1, int resno2, const char *altconf) {
   int status = 0;
   
   SCM rv = SCM_BOOL_F;
   
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      // the "" is the insertion code (not passed to this function (yet)
      int index1 = graphics_info_t::molecules[imol].atom_index_first_atom_in_residue(chain_id, resno1, ""); 
      int index2 = graphics_info_t::molecules[imol].atom_index_first_atom_in_residue(chain_id, resno2, "");
      short int auto_range = 0;
      if (index1 >= 0) {
	 if (index2 >= 0) { 
	    coot::refinement_results_t rr = g.regularize(imol, auto_range, index1, index2);
	    std::cout << "debug:: restraints results " << rr.found_restraints_flag << " "
		      << rr.lights.size() << " " << rr.info_text << std::endl;
	    if ((rr.lights.size() > 0) || (rr.found_restraints_flag)) {
	       rv = g.refinement_results_to_scm(rr);
	    }

	 } else {
	    std::cout << "WARNING:: regularize_zone: Can't get index for resno2: "
		      << resno2 << std::endl;
	 } 
      } else {
	 std::cout << "WARNING:: regularize_zone: Can't get index for resno1: "
		   << resno1 << std::endl;
      }
   } else {
      std::cout << "Not a valid model molecule" << std::endl;
   }

   return rv;
}
#endif // USE_GUILE

int
use_only_extra_torsion_restraints_for_torsions_state() {
   return graphics_info_t:: use_only_extra_torsion_restraints_for_torsions_flag;

}

#ifdef USE_PYTHON
PyObject *python_representation_kk(int imol) {
   
   PyObject *r = Py_False;
   
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      int nchains = mol->GetNumberOfChains(1);
      
      PyObject *chain_list = PyList_New(nchains);
      
      for (int ichain=0; ichain<nchains; ichain++) {
	 mmdb::Chain *chain_p = mol->GetChain(1,ichain);
	 PyObject *chain_id = myPyString_FromString(chain_p->GetChainID());
         int nres;
         mmdb::PResidue *residues;
         chain_p->GetResidueTable(residues, nres);
         
         PyObject *chain_info = PyList_New(2);
         PyList_SetItem(chain_info, 0, chain_id);
         PyObject *res_list = PyList_New(nres);
         
         for (int ires=0; ires<nres; ires++) {
            mmdb::Residue *this_res = residues[ires];
            
            //get the residue name, number, and insertion code
            PyObject *res_name, *res_num, *res_inscode;
            res_name    = myPyString_FromString(this_res->GetResName());
            res_num     = PyLong_FromLong(this_res->GetSeqNum());
            res_inscode = myPyString_FromString(this_res->GetInsCode());
            
            //store the residue name, number, and insertion code
            PyObject *res_info = PyList_New(4);
            PyList_SetItem(res_info, 0, res_num);
            PyList_SetItem(res_info, 1, res_inscode);
            PyList_SetItem(res_info, 2, res_name);
            
            //get the list of atoms (taken from residue_info_py)
            int n_atoms = this_res->GetNumberOfAtoms();
            PyObject *at_info = Py_False;
            PyObject *at_pos;
            PyObject *at_occ, *at_b, *at_biso, *at_ele, *at_name, *at_altconf;
            PyObject *at_segid;
            PyObject *at_x, *at_y, *at_z;
            PyObject *compound_name;
            PyObject *compound_attrib;
            PyObject *all_atoms = PyList_New(0);
            for (int iat=0; iat<n_atoms; iat++) {

               mmdb::Atom *at = this_res->GetAtom(iat);
               if (at->Ter) continue; //ignore TER records
               
               at_x  = PyFloat_FromDouble(at->x);
               at_y  = PyFloat_FromDouble(at->y);
               at_z  = PyFloat_FromDouble(at->z);
               at_pos = PyList_New(3);
               PyList_SetItem(at_pos, 0, at_x);
               PyList_SetItem(at_pos, 1, at_y);
               PyList_SetItem(at_pos, 2, at_z);

               at_occ = PyFloat_FromDouble(at->occupancy);
               at_biso= PyFloat_FromDouble(at->tempFactor);
               at_ele = myPyString_FromString(at->element);
               at_name = myPyString_FromString(at->name);
               at_segid = myPyString_FromString(at->segID);
               at_altconf = myPyString_FromString(at->altLoc);

               at_b = at_biso;
               if (at->WhatIsSet & mmdb::ASET_Anis_tFac) {
                  at_b = PyList_New(7);
                  PyList_SetItem(at_b, 0, at_biso);
                  PyList_SetItem(at_b, 1, PyFloat_FromDouble(at->u11));
                  PyList_SetItem(at_b, 2, PyFloat_FromDouble(at->u22));
                  PyList_SetItem(at_b, 3, PyFloat_FromDouble(at->u33));
                  PyList_SetItem(at_b, 4, PyFloat_FromDouble(at->u12));
                  PyList_SetItem(at_b, 5, PyFloat_FromDouble(at->u13));
                  PyList_SetItem(at_b, 6, PyFloat_FromDouble(at->u23));
               }

               compound_name = PyList_New(2);
               PyList_SetItem(compound_name, 0 ,at_name);
               PyList_SetItem(compound_name, 1 ,at_altconf);

               compound_attrib = PyList_New(4);
               PyList_SetItem(compound_attrib, 0, at_occ);
               PyList_SetItem(compound_attrib, 1, at_b);
               PyList_SetItem(compound_attrib, 2, at_ele);
               PyList_SetItem(compound_attrib, 3, at_segid);

               at_info = PyList_New(3);
               PyList_SetItem(at_info, 0, compound_name);
               PyList_SetItem(at_info, 1, compound_attrib);
               PyList_SetItem(at_info, 2, at_pos);

               PyList_Append(all_atoms, at_info);
            }
            
            //store the information about the current residue
            PyList_SetItem(res_info, 3, all_atoms);
            PyList_SetItem(res_list, ires, res_info);
         }
         
         //store the information about the current chain
         PyList_SetItem(chain_info, 1, res_list);
         PyList_SetItem(chain_list, ichain, chain_info);
      }
      
      //store the chain_list in r
      r = PyList_New(1);
      PyList_SetItem(r, 0, chain_list);
      //if we were retrieving multiple models for the current molecule, we would store chain_list at position imodel
   }
   
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
   
}
#endif // USE_PYTHON
