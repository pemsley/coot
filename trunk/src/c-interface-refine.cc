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

#include <stdlib.h>
#include <iostream>


#include "guile-fixups.h"

#include "graphics-info.h"

// Including python needs to come after graphics-info.h, because
// something in Python.h (2.4 - chihiro) is redefining FF1 (in
// ssm_superpose.h) to be 0x00004000 (Grrr).
//
#ifdef USE_PYTHON
#include "Python.h"
#endif // USE_PYTHON

#include "c-interface.h"
#include "cc-interface.hh"


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
      CResidue *res_1 = g.molecules[imol].get_residue(resno1, "", chain_id);
      CResidue *res_2 = g.molecules[imol].get_residue(resno2, "", chain_id);
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
void regularize_zone(int imol, const char *chain_id, int resno1, int resno2, const char *altconf) {
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      // the "" is the insertion code (not passed to this function (yet)
      int index1 = graphics_info_t::molecules[imol].atom_index_first_atom_in_residue(chain_id, resno1, ""); 
      int index2 = graphics_info_t::molecules[imol].atom_index_first_atom_in_residue(chain_id, resno2, "");
      short int auto_range = 0;
      if (index1 >= 0) {
	 if (index2 >= 0) { 
	    g.regularize(imol, auto_range, index1, index2);
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
      CResidue *res_1 = g.molecules[imol].get_residue(resno1, inscode_1, chain_id);
      CResidue *res_2 = g.molecules[imol].get_residue(resno2, inscode_2, chain_id);
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
      CResidue *res_1 = g.molecules[imol].get_residue(resno1, inscode_1, chain_id);
      CResidue *res_2 = g.molecules[imol].get_residue(resno2, inscode_2, chain_id);
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

   return r;
}
#endif // USE_PYTHON

