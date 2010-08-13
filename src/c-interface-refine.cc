/* src/c-interface-refine.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007, 2008 The University of York
 * Author: Paul Emsley
 * Copyright 2007 by Paul Emsley
 * Copyright 2007 by Bernhard Lohkamp
 * Copyright 2008 by Kevin Cowtan
 * Copyright 2007, 2008, 2009 The University of Oxford
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
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <stdlib.h>
#include <iostream>


#include "guile-fixups.h"

#include "graphics-info.h"

// Including python needs to come after graphics-info.h, because
// something in Python.h (2.4 - chihiro) is redefining FF1 (in
// ssm_superpose.h) to be 0x00004000 (Grrr).
//
// 20100813: Python.h needs to come before to stop"_POSIX_C_SOURCE" redefined problems 
//
// #ifdef USE_PYTHON
// #include "Python.h"
// #endif // USE_PYTHON

#include "c-interface.h"
#include "cc-interface.hh"
#include "c-interface-scm.hh" // for display_scm



/*  ----------------------------------------------------------------------- */
//                                 refinement
/*  ----------------------------------------------------------------------- */

void set_matrix(float f) {

   graphics_info_t::geometry_vs_map_weight = f;
}

float matrix_state() {
   return graphics_info_t::geometry_vs_map_weight;
}

void set_refine_auto_range_step(int i) { 
   graphics_info_t::refine_auto_range_step = i;
} 

void set_refine_max_residues(int n) { 
   graphics_info_t::refine_regularize_max_residues = n;
}

// (refine-zone-atom-index-define 0 688 688)
// 
void refine_zone_atom_index_define(int imol, int ind1, int ind2) {
   
   graphics_info_t g;

   if (is_valid_model_molecule(imol)) {
      if (g.molecules[imol].has_model()) {
	 if (g.molecules[imol].atom_sel.n_selected_atoms > ind1 &&
	     g.molecules[imol].atom_sel.n_selected_atoms > ind2) {
	    g.refine(imol, 0, ind1, ind2);
	 } else {
	    std::cout << "WARNING: atom index error in "
		      << "refine_zone_atom_index_define\n";
	 }
      } else {
	 std::cout << "WARNING: no model for molecule " << imol << " in "
		   << "refine_zone_atom_index_define\n";
      }
   } else {
      std::cout << "WARNING: no molecule " << imol << " in "
		<< "refine_zone_atom_index_define\n";
   }
}

void refine_zone(int imol, const char *chain_id,
		 int resno1,
		 int resno2,
		 const char *altconf) {

   graphics_info_t g;
   if (is_valid_model_molecule(imol)) {
      CResidue *res_1 = g.molecules[imol].get_residue(chain_id, resno1, "");
      CResidue *res_2 = g.molecules[imol].get_residue(chain_id, resno2, "");
      if (res_1 && res_2) { 
	 std::string resname_1(res_1->GetResName());
	 std::string resname_2(res_2->GetResName());
	 bool is_water_like_flag = g.check_for_no_restraints_object(resname_1, resname_2);
	 g.refine_residue_range(imol, chain_id, chain_id, resno1, "", resno2, "", altconf,
				is_water_like_flag);
      }
   }
}


void refine_auto_range(int imol, const char *chain_id, int resno1, const char *altconf) {


   if (is_valid_model_molecule(imol)) { 
      graphics_info_t g;
      int index1 = atom_index(imol, chain_id, resno1, " CA ");
      short int auto_range = 1;
      if (index1 >= 0) { 
	 g.refine(imol, auto_range, index1, index1);
      } else {
	 std::cout << "WARNING:: refine_auto_range: Can't get index for resno1: "
		   << resno1 << std::endl;
      }
   }
}

/*! \brief regularize a zone
  */
int regularize_zone(int imol, const char *chain_id, int resno1, int resno2, const char *altconf) {
   int status = 0;
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
		      << rr.lights.size() << " " << rr.info << std::endl;
	    if (rr.lights.size() > 0)
	       status = 1;
	    if (rr.found_restraints_flag)
	       status = 1;
	    
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
   return status;
} 


// This does not control if the atoms are accepted immediately, just
// whether the Accept Refinemnt gui is shown.
// 
void set_refinement_immediate_replacement(int istate) {
   graphics_info_t::refinement_immediate_replacement_flag = istate;
}

int  refinement_immediate_replacement_state() {
   return graphics_info_t::refinement_immediate_replacement_flag; 
} 


int imol_refinement_map() {

   graphics_info_t g;
   return g.Imol_Refinement_Map();

}

int set_imol_refinement_map(int imol) {

   int r = -1; 
   if (is_valid_map_molecule(imol)) { 
      graphics_info_t g;
      r = g.set_imol_refinement_map(imol);
   }
   return r;
}

/*  ----------------------------------------------------------------------- */
/*                  regularize/refine                                       */
/*  ----------------------------------------------------------------------- */

void do_regularize(short int state) { 

   //
   graphics_info_t g; 

   g.set_in_range_define_for_regularize(state);  // TRUE or FALSE
   if (state) { 
      g.untoggle_model_fit_refine_buttons_except("model_refine_dialog_regularize_zone_togglebutton");
      // and kill the delete dialog if it is there
      do_regularize_kill_delete_dialog();
      std::cout << "click on 2 atoms (in the same molecule)" << std::endl; 
      g.pick_cursor_maybe();
      g.pick_pending_flag = 1;
   } else { 
      g.normal_cursor();
   }
}

void do_refine(short int state) { 

   //
   graphics_info_t g; 

   g.set_in_range_define_for_refine(state);  // TRUE or Not...

   // std::cout << "DEBUG:: in do_refine" << std::endl;
   
   if (state) { 
      g.untoggle_model_fit_refine_buttons_except("model_refine_dialog_refine_togglebutton");
      // and kill the delete dialog if it is there
      do_regularize_kill_delete_dialog();
      
      int imol_map = g.Imol_Refinement_Map();
      // std::cout << "DEBUG:: in do_refine, imol_map: " << imol_map << std::endl;
      if (imol_map >= 0) {
	 if (g.molecules[imol_map].has_map()) { 
	    std::cout << "click on 2 atoms (in the same molecule)" << std::endl; 
	    g.pick_cursor_maybe();
	    g.pick_pending_flag = 1;
	    std::string s = "Pick 2 atoms or Autozone (pick 1 atom the press the A key)";
	    s += " [Ctrl Left-mouse rotates the view]";
	    s += "...";
	    g.statusbar_text(s);
	 } else {
	    g.show_select_map_dialog();
	    g.in_range_define_for_refine = 0;
	    g.model_fit_refine_unactive_togglebutton("model_refine_dialog_refine_togglebutton");
	 }
      } else {
	 // map chooser dialog
	 g.show_select_map_dialog();
	 g.in_range_define_for_refine = 0;
	 g.model_fit_refine_unactive_togglebutton("model_refine_dialog_refine_togglebutton");
      }
   } else { 
      g.normal_cursor();
      g.in_range_define_for_refine = 0;
      // g.pick_pending_flag = 0;
   }

}

void set_residue_selection_flash_frames_number(int i) {

   graphics_info_t::residue_selection_flash_frames_number = i;
}

void set_refinement_refine_per_frame(int i) {

   graphics_info_t::dragged_refinement_refine_per_frame_flag = i; 
}

int  refinement_refine_per_frame_state() {

   return graphics_info_t::dragged_refinement_refine_per_frame_flag;
}

#ifdef USE_GUILE
SCM refine_zone_with_full_residue_spec_scm(int imol, const char *chain_id,
					   int resno1,
					   const char *inscode_1,
					   int resno2,
					   const char *inscode_2,
					   const char *altconf) {
   SCM r = SCM_BOOL_F;
   graphics_info_t g;
   if (is_valid_model_molecule(imol)) {
      CResidue *res_1 = g.molecules[imol].get_residue(chain_id, resno1, inscode_1);
      CResidue *res_2 = g.molecules[imol].get_residue(chain_id, resno2, inscode_2);
      if (res_1 && res_2) { 
	 std::string resname_1(res_1->GetResName());
	 std::string resname_2(res_2->GetResName());
	 bool is_water_like_flag = g.check_for_no_restraints_object(resname_1, resname_2);
	 coot::refinement_results_t rr = 
	    g.refine_residue_range(imol, chain_id, chain_id, resno1, "", resno2, "", altconf,
				   is_water_like_flag);
	 r = g.refinement_results_to_scm(rr);
      }
   }

   return r;

}
#endif // USE_GUILE


#ifdef USE_PYTHON
PyObject *refine_zone_with_full_residue_spec_py(int imol, const char *chain_id,
						int resno1,
						const char*inscode_1,
						int resno2,
						const char*inscode_2,
						const char *altconf) {
   PyObject *r = Py_False;
   graphics_info_t g;
   if (is_valid_model_molecule(imol)) {
      CResidue *res_1 = g.molecules[imol].get_residue(chain_id, resno1, inscode_1);
      CResidue *res_2 = g.molecules[imol].get_residue(chain_id, resno2, inscode_2);
      if (res_1 && res_2) { 
	 std::string resname_1(res_1->GetResName());
	 std::string resname_2(res_2->GetResName());
	 bool is_water_like_flag = g.check_for_no_restraints_object(resname_1, resname_2);
	 coot::refinement_results_t rr = 
	    g.refine_residue_range(imol, chain_id, chain_id, resno1, "", resno2, "", altconf,
				   is_water_like_flag);
	 r = g.refinement_results_to_py(rr);
      }
   }

   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif // USE_PYTHON

/*! \brief add a user-define bond restraint

   to be used when the given atoms are selected.  */
int add_extra_bond_restraint(int imol, const char *chain_id_1, int res_no_1, const char *ins_code_1, const char *atom_name_1, const char *alt_conf_1, const char *chain_id_2, int res_no_2, const char *ins_code_2, const char *atom_name_2, const char *alt_conf_2, double bond_dist, double esd) {

   int r = -1;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t as_1(chain_id_1, res_no_1, ins_code_1, atom_name_1, alt_conf_1);
      coot::atom_spec_t as_2(chain_id_2, res_no_2, ins_code_2, atom_name_2, alt_conf_2);
      r = graphics_info_t::molecules[imol].add_extra_bond_restraint(as_1, as_2, bond_dist, esd);
      graphics_draw();
   }
   return r;

} 


int add_extra_torsion_restraint(int imol, 
				const char *chain_id_1, int res_no_1, const char *ins_code_1, const char *atom_name_1, const char *alt_conf_1, 
				const char *chain_id_2, int res_no_2, const char *ins_code_2, const char *atom_name_2, const char *alt_conf_2, 
				const char *chain_id_3, int res_no_3, const char *ins_code_3, const char *atom_name_3, const char *alt_conf_3, 
				const char *chain_id_4, int res_no_4, const char *ins_code_4, const char *atom_name_4, const char *alt_conf_4, 
				double torsion_angle, double esd, int period) {

   int r = -1;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t as_1(chain_id_1, res_no_1, ins_code_1, atom_name_1, alt_conf_1);
      coot::atom_spec_t as_2(chain_id_2, res_no_2, ins_code_2, atom_name_2, alt_conf_2);
      coot::atom_spec_t as_3(chain_id_3, res_no_3, ins_code_3, atom_name_3, alt_conf_3);
      coot::atom_spec_t as_4(chain_id_4, res_no_4, ins_code_4, atom_name_4, alt_conf_4);
      r = graphics_info_t::molecules[imol].add_extra_torsion_restraint(as_1, as_2, as_3, as_4, torsion_angle, esd, period);
      graphics_draw();
   }
   

   return r;
}


#ifdef USE_GUILE
SCM list_extra_restraints_scm(int imol) {

   SCM r = SCM_BOOL_F;

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g; // just because it's shorter
      if (graphics_info_t::molecules[imol].extra_restraints.has_restraints()) {
	 // reverse loop, put them in backwards (schemey thing)
	 r = SCM_EOL;
	 for (int ib=g.molecules[imol].extra_restraints.bond_restraints.size()-1; ib>=0; ib--) {
	    coot::atom_spec_t spec_1 = g.molecules[imol].extra_restraints.bond_restraints[ib].atom_1;
	    coot::atom_spec_t spec_2 = g.molecules[imol].extra_restraints.bond_restraints[ib].atom_2;
	    double d = g.molecules[imol].extra_restraints.bond_restraints[ib].bond_dist;
	    double esd = g.molecules[imol].extra_restraints.bond_restraints[ib].esd;
	    SCM spec_1_scm = atom_spec_to_scm(spec_1);
	    SCM spec_2_scm = atom_spec_to_scm(spec_2);
	    SCM l = scm_list_4(spec_1_scm, spec_2_scm, scm_double2num(d), scm_double2num(esd));
	    l = scm_cons(scm_str2symbol("bond"), l);
	    r = scm_cons(l, r);
	 }

	 for (int it=g.molecules[imol].extra_restraints.torsion_restraints.size()-1; it>=0; it--) {
	    coot::atom_spec_t spec_1 = g.molecules[imol].extra_restraints.torsion_restraints[it].atom_1;
	    coot::atom_spec_t spec_2 = g.molecules[imol].extra_restraints.torsion_restraints[it].atom_2;
	    coot::atom_spec_t spec_4 = g.molecules[imol].extra_restraints.torsion_restraints[it].atom_3;
	    coot::atom_spec_t spec_3 = g.molecules[imol].extra_restraints.torsion_restraints[it].atom_4;
	    double t = g.molecules[imol].extra_restraints.torsion_restraints[it].torsion_angle;
	    double e = g.molecules[imol].extra_restraints.torsion_restraints[it].esd;
	    int    p = g.molecules[imol].extra_restraints.torsion_restraints[it].period;
	    SCM l = scm_list_3(scm_double2num(t), scm_double2num(e), scm_int2num(p));
	    l = scm_cons(atom_spec_to_scm(spec_4), l);
	    l = scm_cons(atom_spec_to_scm(spec_3), l);
	    l = scm_cons(atom_spec_to_scm(spec_2), l);
	    l = scm_cons(atom_spec_to_scm(spec_1), l);
	    l = scm_cons(scm_str2symbol("torsion"), l);
	    r = scm_cons(l, r);
	 } 
	 
      }
   }
   return r;
} 
#endif	/* USE_GUILE */


#ifdef USE_PYTHON
PyObject *list_extra_restraints_py(int imol) {

   PyObject *r = Py_False;

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      if (graphics_info_t::molecules[imol].extra_restraints.has_restraints()) {
	 r = PyList_New(0);
	 for (int ib=0; ib<g.molecules[imol].extra_restraints.bond_restraints.size(); ib++) {
	    coot::atom_spec_t spec_1 = g.molecules[imol].extra_restraints.bond_restraints[ib].atom_1;
	    coot::atom_spec_t spec_2 = g.molecules[imol].extra_restraints.bond_restraints[ib].atom_2;
	    double d = g.molecules[imol].extra_restraints.bond_restraints[ib].bond_dist;
	    double esd = g.molecules[imol].extra_restraints.bond_restraints[ib].esd;
	    PyObject *spec_1_py = atom_spec_to_py(spec_1);
	    PyObject *spec_2_py = atom_spec_to_py(spec_2);
	    PyObject *l = PyList_New(5);
	    PyList_SetItem(l, 0, PyString_FromString("bond"));
	    PyList_SetItem(l, 1, spec_1_py);
	    PyList_SetItem(l, 2, spec_2_py);
	    PyList_SetItem(l, 3, PyFloat_FromDouble(d));
	    PyList_SetItem(l, 4, PyFloat_FromDouble(esd));
	    PyList_Append(r, l);
	 }

	 for (int it=0; it<g.molecules[imol].extra_restraints.torsion_restraints.size(); it++) {
	    coot::atom_spec_t spec_1 = g.molecules[imol].extra_restraints.torsion_restraints[it].atom_1;
	    coot::atom_spec_t spec_2 = g.molecules[imol].extra_restraints.torsion_restraints[it].atom_2;
	    coot::atom_spec_t spec_3 = g.molecules[imol].extra_restraints.torsion_restraints[it].atom_3;
	    coot::atom_spec_t spec_4 = g.molecules[imol].extra_restraints.torsion_restraints[it].atom_4;
	    PyObject *spec_1_py = atom_spec_to_py(spec_1);
	    PyObject *spec_2_py = atom_spec_to_py(spec_2);
	    PyObject *spec_3_py = atom_spec_to_py(spec_3);
	    PyObject *spec_4_py = atom_spec_to_py(spec_4);
	    double t = g.molecules[imol].extra_restraints.torsion_restraints[it].torsion_angle;
	    double e = g.molecules[imol].extra_restraints.torsion_restraints[it].esd;
	    int    p = g.molecules[imol].extra_restraints.torsion_restraints[it].period;
	    PyObject *l = PyList_New(8);
	    PyList_SetItem(l, 0, PyString_FromString("torsion"));
	    PyList_SetItem(l, 1, spec_1_py);
	    PyList_SetItem(l, 2, spec_2_py);
	    PyList_SetItem(l, 3, spec_3_py);
	    PyList_SetItem(l, 4, spec_4_py);
	    PyList_SetItem(l, 5, PyFloat_FromDouble(t));
	    PyList_SetItem(l, 6, PyFloat_FromDouble(e));
	    PyList_SetItem(l, 7, PyInt_FromLong(p));
	    PyList_Append(r, l);
	 }
      }
      if (PyBool_Check(r)) {
	 Py_INCREF(r);
      }
   }
   return r;
} 
#endif	/* USE_PYTHON */


#ifdef USE_GUILE
void
delete_extra_restraint_scm(int imol, SCM restraint_spec) {

   // for a bond restraint, the restraint_spec is something like:
   // (list restraint-type spec-1 spec-2)
   //
   // where restraint-type is a symbol, in the case of a bond
   // restraint is 'bond
   //
   //
   if (scm_is_true(scm_list_p(restraint_spec))) { 
      SCM restraint_spec_length_scm = scm_length(restraint_spec);
      int restraint_spec_length = scm_to_int(restraint_spec_length_scm);
      if (restraint_spec_length == 3) {
	 SCM restraint_type_scm = SCM_CAR(restraint_spec);
	 SCM spec_1_scm = scm_list_ref(restraint_spec, SCM_MAKINUM(1));
	 SCM spec_2_scm = scm_list_ref(restraint_spec, SCM_MAKINUM(2));
	 if (scm_is_true(scm_eq_p(restraint_type_scm, scm_str2symbol("bond")))) {
	    coot::atom_spec_t spec_1 = atom_spec_from_scm_expression(spec_1_scm);
	    coot::atom_spec_t spec_2 = atom_spec_from_scm_expression(spec_2_scm);
	    graphics_info_t::molecules[imol].remove_extra_bond_restraint(spec_1, spec_2);
	    graphics_draw();
	 }
      }
   }
   
} 
#endif // USE_GUILE

#ifdef USE_PYTHON
void
delete_extra_restraint_py(int imol, PyObject *restraint_spec) {

   // for a bond restraint, the restraint_spec is something like:
   // [restraint-type, spec-1 spec-2]
   //
   // where restraint-type is a symbol, in the case of a bond
   // restraint is "bond"
   //
   //
  if (PyList_Check(restraint_spec)) { 
    int restraint_spec_length = PyObject_Length(restraint_spec);
    if (restraint_spec_length == 3) {
      PyObject *restraint_type_py = PyList_GetItem(restraint_spec, 0);
      PyObject *spec_1_py = PyList_GetItem(restraint_spec, 1);
      PyObject *spec_2_py = PyList_GetItem(restraint_spec, 2);
      if (strcmp(PyString_AsString(restraint_type_py), "bond") == 0 ) {
	    coot::atom_spec_t spec_1 = atom_spec_from_python_expression(spec_1_py);
	    coot::atom_spec_t spec_2 = atom_spec_from_python_expression(spec_2_py);
	    graphics_info_t::molecules[imol].remove_extra_bond_restraint(spec_1, spec_2);
	    graphics_draw();
      }
    }
  }
   
} 
#endif // USE_PYTHON

void
delete_all_extra_restraints(int imol) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].extra_restraints.clear();
   } 
} 

void
set_use_only_extra_torsion_restraints_for_torsions(short int state) {

   graphics_info_t:: use_only_extra_torsion_restraints_for_torsions_flag = state;

} 
