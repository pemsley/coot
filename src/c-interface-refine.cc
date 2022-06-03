/* src/c-interface-refine.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007, 2008 The University of York
 * Author: Paul Emsley
 * Copyright 2007 by Paul Emsley
 * Copyright 2007 by Bernhard Lohkamp
 * Copyright 2008 by Kevin Cowtan
 * Copyright 2007, 2008, 2009 The University of Oxford
 * Copyright 2015 by Medical Research Council
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
#include "python-3-interface.hh"
#endif

#include <stdlib.h>
#include <iostream>

#include "compat/coot-sysdep.h"


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
#include "c-interface-gtk-widgets.h"
#include "cc-interface.hh"
#include "c-interface-scm.hh" // for display_scm
#include "c-interface-refine.h" // for crankshaft_peptide_rotation_optimization_intermediate_atoms()


/*  ----------------------------------------------------------------------- */
//                                 refinement
/*  ----------------------------------------------------------------------- */

void set_matrix(float f) {

   graphics_info_t::geometry_vs_map_weight = f;
}

float matrix_state() {
   return graphics_info_t::geometry_vs_map_weight;
}

float get_map_weight() {
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
   g.conditionally_wait_for_refinement_to_finish();
}

void refine_zone(int imol, const char *chain_id,
		 int resno1,
		 int resno2,
		 const char *altconf) {

   graphics_info_t g;
   g.residue_type_selection_was_user_picked_residue_range = false;
   if (is_valid_model_molecule(imol)) {
      mmdb::Residue *res_1 = g.molecules[imol].get_residue(chain_id, resno1, "");
      mmdb::Residue *res_2 = g.molecules[imol].get_residue(chain_id, resno2, "");
      if (res_1 && res_2) { 
	 std::string resname_1(res_1->GetResName());
	 std::string resname_2(res_2->GetResName());
	 bool is_water_like_flag = g.check_for_no_restraints_object(resname_1, resname_2);
	 // g.refine_residue_range(imol, chain_id, chain_id, resno1, "", resno2, "", altconf,
         //                        is_water_like_flag);
         mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
         std::vector<mmdb::Residue *> residues = coot::util::get_residues_in_range(mol, chain_id, resno1, resno2);

         std::string alt_conf(altconf);
         if (! residues.empty())
            coot::refinement_results_t rr = g.refine_residues_vec(imol, residues, alt_conf, mol);
      }
   }
   g.conditionally_wait_for_refinement_to_finish();
}

/* use stored atom indices to re-run the refinement using the same atoms as previous */
void repeat_refine_zone() {

   graphics_info_t g;
   g.repeat_refine_zone();
   g.conditionally_wait_for_refinement_to_finish();
}



void refine_auto_range(int imol, const char *chain_id, int resno1, const char *altconf) {


   if (is_valid_model_molecule(imol)) { 
      graphics_info_t g;
      int index1 = atom_index_full(imol, chain_id, resno1, "", " CA ", altconf);
      short int auto_range = 1;
      if (index1 >= 0) { 
	 g.refine(imol, auto_range, index1, index1);
      } else {
	 std::cout << "WARNING:: refine_auto_range: Can't get index for resno1: "
		   << resno1 << std::endl;
      }
      g.conditionally_wait_for_refinement_to_finish();
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
		      << rr.lights.size() << " " << rr.info_text << std::endl;
	    if (rr.lights.size() > 0)
	       status = 1;
	    if (rr.found_restraints_flag)
	       status = 1;
            g.conditionally_wait_for_refinement_to_finish();
	    
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
// whether the Accept Refinement gui is shown.
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
      // g.untoggle_model_fit_refine_buttons_except("model_refine_dialog_regularize_zone_togglebutton");
      // and kill the delete dialog if it is there
      // do_regularize_kill_delete_dialog(); // no longer
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
      // g.untoggle_model_fit_refine_buttons_except("model_refine_dialog_refine_togglebutton");
      // and kill the delete dialog if it is there
      // do_regularize_kill_delete_dialog(); // 20220602-PE no longer
      
      int imol_map = g.Imol_Refinement_Map();
      // std::cout << "DEBUG:: in do_refine, imol_map: " << imol_map << std::endl;
      if (imol_map < 0) {
          g.show_select_map_dialog();
          imol_map = g.Imol_Refinement_Map();
      }
      if (imol_map >= 0) {
	 if (g.molecules[imol_map].has_xmap()) { 
	    std::cout << "click on 2 atoms (in the same molecule)" << std::endl; 
	    g.pick_cursor_maybe();
	    g.pick_pending_flag = 1;
	    std::string s = "Pick 2 atoms or Autozone (pick 1 atom then press the A key)";
	    s += " [Ctrl Left-mouse rotates the view]";
	    s += "...";
	    g.add_status_bar_text(s);
	 } else {
	    g.show_select_map_dialog();
	    g.in_range_define_for_refine = 0;
	    g.model_fit_refine_unactive_togglebutton("model_refine_dialog_refine_togglebutton");
	 }
      } else {
	 // map chooser dialog
	 //g.show_select_map_dialog();
         // shouldnt get here any more?!? Only if we destroy the dialog above!
          g.in_range_define_for_refine = 0;
          g.model_fit_refine_unactive_togglebutton("model_refine_dialog_refine_togglebutton");
          info_dialog("WARNING:: Still, no refinement map has been set!");
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

/*! \brief - the elasticity of the dragged atom in refinement mode.

Default 0.1

 Bigger numbers mean bigger movement of the other atoms.*/
void set_refinement_drag_elasticity(float e) {
   graphics_info_t::refinement_drag_elasticity = e;
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
      mmdb::Residue *res_1 = g.molecules[imol].get_residue(chain_id, resno1, inscode_1);
      mmdb::Residue *res_2 = g.molecules[imol].get_residue(chain_id, resno2, inscode_2);
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

std::string mtz_file_name(int imol) {

   std::string r;
   if (is_valid_map_molecule(imol)) {
      r = graphics_info_t::molecules[imol].Refmac_mtz_filename();
   }
   return r;
}



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
      mmdb::Residue *res_1 = g.molecules[imol].get_residue(chain_id, resno1, inscode_1);
      mmdb::Residue *res_2 = g.molecules[imol].get_residue(chain_id, resno2, inscode_2);
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



#ifdef USE_PYTHON
PyObject *residues_distortions_py(int imol, PyObject *residue_specs_list_py) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      std::vector<coot::residue_spec_t> residue_specs = py_to_residue_specs(residue_specs_list_py);
      if (residue_specs.size() > 0) {
	 std::vector<mmdb::Residue *> residues;
	 for (unsigned int i=0; i<residue_specs.size(); i++) {
	    coot::residue_spec_t rs = residue_specs[i];
	    mmdb::Residue *r = graphics_info_t::molecules[imol].get_residue(rs);
	    if (r) {
	       residues.push_back(r);
	    }
	 }

	 if (residues.size() > 0) {
	    graphics_info_t g;
	    int imol_map = g.Imol_Refinement_Map();
	    if (! is_valid_map_molecule(imol_map)) { 
	       add_status_bar_text("Refinement map not set");
	    } else {
	       // normal
	       mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
	       graphics_info_t g;
	       std::vector<std::pair<bool,mmdb::Residue *> > local_residues;  // not fixed.
	       for (unsigned int i=0; i<residues.size(); i++)
		  local_residues.push_back(std::pair<bool, mmdb::Residue *>(false, residues[i]));
	       const coot::protein_geometry &geom = *g.Geom_p();
	       bool do_residue_internal_torsions = false;
	       bool do_trans_peptide_restraints = false;
	       bool do_rama_restraints = false;
	       float rama_plot_restraint_weight = 1.0;
	       coot::pseudo_restraint_bond_type pseudo_bonds_type = coot::NO_PSEUDO_BONDS;
	       const clipper::Xmap<float> &xmap = graphics_info_t::molecules[imol_map].xmap;
	       std::vector<coot::atom_spec_t> fixed_atom_specs;
	       std::vector<mmdb::Link> links;
	       coot::restraint_usage_Flags flags = coot::TYPICAL_RESTRAINTS;

	       coot::restraints_container_t restraints(local_residues, links, geom, mol, fixed_atom_specs, &xmap);
               unsigned int n_threads = coot::get_max_number_of_threads();
               restraints.thread_pool(&g.static_thread_pool, n_threads);
	       restraints.make_restraints(imol, geom, flags,
                                          do_residue_internal_torsions,
                                          do_trans_peptide_restraints,
                                          rama_plot_restraint_weight,
                                          do_rama_restraints,
                                          false, false, false,
                                          pseudo_bonds_type);
	       coot::geometry_distortion_info_container_t gd = restraints.geometric_distortions();
	       // std::cout << "Found " << gd.size() << " geometry distortions" << std::endl;
	       if (gd.size() > 0) {

		  r = PyList_New(gd.size());
		  for (std::size_t i=0; i<gd.geometry_distortion.size(); i++) {

		     PyList_SetItem(r, i, g.geometry_distortion_to_py(gd.geometry_distortion[i]));

		     if (false) { // debug to terminal
			std::cout << "   " << i << " ";
			for (std::size_t j=0; j<gd.geometry_distortion[i].atom_indices.size(); j++)
			   std::cout << " " << gd.geometry_distortion[i].atom_indices[j];
			std::cout << ": ";
			std::cout << " " << gd.geometry_distortion[i] << " ";
			std::cout << gd.geometry_distortion[i].distortion_score << std::endl;
		     }
		  }
	       }
	    }
	 }
      }
   }
   
   if (PyBool_Check(r))
      Py_INCREF(r);
   
   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *get_intermediate_atoms_distortions_py() {

   graphics_info_t g;
   PyObject *r = g.get_intermediate_atoms_distortions_py();
   if (PyBool_Check(r))
      Py_INCREF(r);
   return r;
}
#endif

#ifdef USE_GUILE
SCM residues_distortions_scm(int imol, SCM residue_specs_scm) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      std::vector<coot::residue_spec_t> residue_specs = scm_to_residue_specs(residue_specs_scm);
      if (! residue_specs.empty()) {
	 std::vector<mmdb::Residue *> residues;
	 for (std::size_t i=0; i<residue_specs.size(); i++) {
	    const coot::residue_spec_t &rs = residue_specs[i];
	    mmdb::Residue *resdiue_p = graphics_info_t::molecules[imol].get_residue(rs);
	    if (r) {
	       residues.push_back(resdiue_p);
	    }
	 }
	 if (! residues.empty()) {
	    graphics_info_t g;
	    int imol_map = g.Imol_Refinement_Map();
	    if (! is_valid_map_molecule(imol_map)) {
	       add_status_bar_text("Refinement map not set");
	    } else {
	       // Happy Path
	       mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;

	       std::vector<std::pair<bool,mmdb::Residue *> > local_residues;  // not fixed.
	       for (unsigned int i=0; i<residues.size(); i++)
		  local_residues.push_back(std::pair<bool, mmdb::Residue *>(false, residues[i]));
	       const coot::protein_geometry &geom = *g.Geom_p();
	       bool do_residue_internal_torsions = false;
	       bool do_trans_peptide_restraints = false;
	       bool do_rama_restraints = false;
	       float rama_plot_restraint_weight = 1.0;
	       coot::pseudo_restraint_bond_type pseudo_bonds_type = coot::NO_PSEUDO_BONDS;
	       const clipper::Xmap<float> &xmap = graphics_info_t::molecules[imol_map].xmap;
	       std::vector<coot::atom_spec_t> fixed_atom_specs;
	       std::vector<mmdb::Link> links;
	       coot::restraint_usage_Flags flags = coot::TYPICAL_RESTRAINTS;

	       coot::restraints_container_t restraints(local_residues, links, geom, mol,
						       fixed_atom_specs, &xmap);
               unsigned int n_threads = coot::get_max_number_of_threads();
               restraints.thread_pool(&g.static_thread_pool, n_threads);
	       // don't do auto-secondary structure restraints
               restraints.make_restraints(imol, geom, flags,
                                          do_residue_internal_torsions,
                                          do_trans_peptide_restraints,
                                          rama_plot_restraint_weight,
                                          do_rama_restraints, false, false, false,
                                          pseudo_bonds_type);
	       coot::geometry_distortion_info_container_t gd = restraints.geometric_distortions();
	       if (gd.size() > 0) {
		  r = SCM_EOL;
		  for (std::size_t i=0; i<gd.geometry_distortion.size(); i++) {
		     SCM gd_scm = g.geometry_distortion_to_scm(gd.geometry_distortion[i]);
		     r = scm_cons(gd_scm, r);
		  }
	       }
	    }
	 }
      }
   }
   return r;
}

#endif // USE_GUILE


/*! read in prosmart (typically) extra restraints */
void add_refmac_extra_restraints(int imol, const char *file_name) {
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].add_refmac_extra_restraints(file_name);
      graphics_draw();
   }
}

void set_cryo_em_refinement(bool mode) {
   graphics_info_t::cryo_EM_refinement_flag = mode;
}

bool get_cryo_em_refinement() {
   return graphics_info_t::cryo_EM_refinement_flag;
}

void write_interpolated_extra_restraints(int imol_1, int imol_2, int n_steps, char *file_name_stub) {

   if (is_valid_model_molecule(imol_1)) {
      if (is_valid_model_molecule(imol_2)) {
	 if (n_steps > 2) {
	    if (n_steps < 5000) { 
	       const coot::extra_restraints_t &e_1 = graphics_info_t::molecules[imol_1].extra_restraints;
	       const coot::extra_restraints_t &e_2 = graphics_info_t::molecules[imol_2].extra_restraints;
	       e_1.write_interpolated_restraints(e_2, n_steps, file_name_stub);
	    } else {
	       std::cout << "too many steps" << std::endl;
	    } 
	 } else {
	    std::cout << "too few steps" << std::endl;
	 }
      }
   }
}

/*! \brief proSMART interpolated restraints for model morphing and write interpolated model

interpolation_mode is currently dummy - in due course I will addd torion angle interpolation.
*/
void write_interpolated_models_and_extra_restraints(int imol_1, int imol_2, int n_steps, char *file_name_stub,
						    int interpolation_mode) {

   if (is_valid_model_molecule(imol_1)) {
      if (is_valid_model_molecule(imol_2)) {
	 if (n_steps > 2) {
	    if (n_steps < 5000) {
	       mmdb::Manager *mol_1 = graphics_info_t::molecules[imol_1].atom_sel.mol;
	       mmdb::Manager *mol_2 = graphics_info_t::molecules[imol_2].atom_sel.mol;

	       if (mol_1 && mol_2) { 
		  const coot::extra_restraints_t &e_1 = graphics_info_t::molecules[imol_1].extra_restraints;
		  const coot::extra_restraints_t &e_2 = graphics_info_t::molecules[imol_2].extra_restraints;
		  e_1.write_interpolated_models_and_restraints(e_2, mol_1, mol_2, n_steps, file_name_stub);
	       }
	    } else {
	       std::cout << "too many steps" << std::endl;
	    } 
	 } else {
	    std::cout << "too few steps" << std::endl;
	 }
      } else {
	 std::cout << "WARNING:: " << imol_2 << " is not a valid model molecule " << std::endl;
      }
   } else {
      std::cout << "WARNING:: " << imol_1 << " is not a valid model molecule " << std::endl;
   }
}

void clear_atom_pull_restraint_on_accept_reject_destroy() {

   graphics_info_t g;
   g.clear_all_atom_pull_restraints(false);
   graphics_draw();
}

void clear_all_atom_pull_restraints() {

   graphics_info_t g;
   g.clear_all_atom_pull_restraints(true);
   graphics_draw();
}

void set_auto_clear_atom_pull_restraint(int state) {
   graphics_info_t g;
   g.auto_clear_atom_pull_restraint_flag = state;
   std::cout << "------------------ set_auto_clear_atom_pull_restraint_state ";
   std::cout << "------------------ update the accept_reject_refinement_dialog here. " << std::endl;
} 

int  get_auto_clear_atom_pull_restraint_state() {
   graphics_info_t g;
   return g.auto_clear_atom_pull_restraint_flag;
}




void set_show_extra_restraints(int imol, int state) {
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_display_extra_restraints(state);
   }
   graphics_draw();
}

void set_show_parallel_plane_restraints(int imol, int state) {
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_display_parallel_plane_restraints(state);
   }
   graphics_draw();
}

int extra_restraints_are_shown(int imol) {
   int r = 0;
   if (is_valid_model_molecule(imol))
      r = graphics_info_t::molecules[imol].draw_it_for_extra_restraints;
   return r;
}

int parallel_plane_restraints_are_shown(int imol) {
   int r = 0;
   if (is_valid_model_molecule(imol))
      r = graphics_info_t::molecules[imol].draw_it_for_parallel_plane_restraints;
   return r;
}

void add_parallel_plane_restraint(int imol,
				  const char *chain_id_1, int res_no_1, const char *ins_code_1,
				  const char *chain_id_2, int res_no_2, const char *ins_code_2) {

   coot::residue_spec_t spec_1(chain_id_1, res_no_1, ins_code_1);
   coot::residue_spec_t spec_2(chain_id_2, res_no_2, ins_code_2);
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].add_parallel_plane_restraint(spec_1, spec_2);
   }
   graphics_draw();

} 


void set_extra_restraints_representation_for_bonds_go_to_CA(int imol, short int state) {

   if (is_valid_model_molecule(imol))
      graphics_info_t::molecules[imol].set_extra_restraints_representation_for_bonds_go_to_CA(state);
   graphics_draw();
} 


/*! \brief often we don't want to see all prosmart restraints, just the (big) violations */
void set_extra_restraints_prosmart_sigma_limits(int imol, double limit_low, double limit_high) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_extra_restraints_prosmart_sigma_limits(limit_low, limit_high);
   }
   graphics_draw();
} 


void generate_local_self_restraints(int imol, const char *chain_id, float local_dist_max) {

   if (is_valid_model_molecule(imol)) {
      // like prosmart self restraints
      graphics_info_t::molecules[imol].generate_local_self_restraints(local_dist_max, chain_id,
								      *graphics_info_t::Geom_p());
   }
   graphics_draw();
}

/*! \brief generate external distance all-molecule self restraints */
void generate_self_restraints(int imol, float local_dist_max) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].generate_self_restraints(local_dist_max,
								*graphics_info_t::Geom_p());
   } 
   graphics_draw();
}





#ifdef USE_GUILE
void generate_local_self_restraints_by_residues_scm(int imol, SCM residue_specs_scm, float local_dist_max) {

   std::vector<coot::residue_spec_t> residue_specs = scm_to_residue_specs(residue_specs_scm);
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].generate_local_self_restraints(local_dist_max, residue_specs,
								      *graphics_info_t::Geom_p());
      graphics_draw();
   }
}
#endif // USE_GUILE
#ifdef USE_PYTHON
void generate_local_self_restraints_by_residues_py(int imol, PyObject *residue_specs_py,
						   float local_dist_max) {

   std::vector<coot::residue_spec_t> residue_specs = py_to_residue_specs(residue_specs_py);
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].generate_local_self_restraints(local_dist_max, residue_specs,
								      *graphics_info_t::Geom_p());
      graphics_draw();
   }
}
#endif // USE_PYTHON


/* ! \brief delete the restraints for the given comp_id (i.e. residue name)  */
// return 0 or 1
// 
int delete_restraints(const char *comp_id) {

   int imol = coot::protein_geometry::IMOL_ENC_ANY; // perhaps this should be passed

   graphics_info_t g;
   return g.Geom_p()->delete_mon_lib(comp_id, imol);
   
}


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

int add_extra_geman_mcclure_restraint(int imol, const char *chain_id_1, int res_no_1, const char *ins_code_1, const char *atom_name_1, const char *alt_conf_1, const char *chain_id_2, int res_no_2, const char *ins_code_2, const char *atom_name_2, const char *alt_conf_2, double bond_dist, double esd) {

   int r = -1;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t as_1(chain_id_1, res_no_1, ins_code_1, atom_name_1, alt_conf_1);
      coot::atom_spec_t as_2(chain_id_2, res_no_2, ins_code_2, atom_name_2, alt_conf_2);
      r = graphics_info_t::molecules[imol].add_extra_geman_mcclure_restraint(as_1, as_2, bond_dist, esd);
      graphics_draw();
   }
   return r;
}


#ifdef USE_GUILE
int add_extra_bond_restraints_scm(int imol, SCM extra_bond_restraints_scm) {

   std::vector<coot::extra_restraints_t::extra_bond_restraint_t> ebr_vec;
   if (is_valid_model_molecule(imol)) {
      if (scm_is_true(scm_list_p(extra_bond_restraints_scm))) {
	 SCM l_scm = scm_length(extra_bond_restraints_scm);
	 int l = scm_to_int(l_scm);
	 for (int i=0; i<l; i++) {
	    SCM restr_descr_scm = scm_list_ref(extra_bond_restraints_scm, scm_from_int(i));
	    if (scm_is_true(scm_list_p(restr_descr_scm))) {
	       SCM r_l_scm = scm_length(restr_descr_scm);
	       int r_l = scm_to_int(r_l_scm);
	       if (r_l == 4) {
		  coot::atom_spec_t atom_1_spec = atom_spec_from_scm_expression(scm_list_ref(restr_descr_scm, scm_from_int(0)));
		  coot::atom_spec_t atom_2_spec = atom_spec_from_scm_expression(scm_list_ref(restr_descr_scm, scm_from_int(1)));
		  SCM target_dist_scm = scm_list_ref(restr_descr_scm, scm_from_int(2));
		  SCM dist_esd_scm    = scm_list_ref(restr_descr_scm, scm_from_int(3));
		  double target_dist = scm_to_double(target_dist_scm);
		  double dist_esd = scm_to_double(dist_esd_scm);
		  coot::extra_restraints_t::extra_bond_restraint_t ebr(atom_1_spec, atom_2_spec, target_dist, dist_esd);
		  ebr_vec.push_back(ebr);
	       }
	    }
	 }
	 int r = graphics_info_t::molecules[imol].add_extra_bond_restraints(ebr_vec);
	 graphics_draw();
      }
   }
   return ebr_vec.size();
}
#endif // USE_GUILE

#ifdef USE_PYTHON
int add_extra_bond_restraints_py(int imol, PyObject *extra_bond_restraints_py) {

   std::vector<coot::extra_restraints_t::extra_bond_restraint_t> ebr_vec;
   if (is_valid_model_molecule(imol)) {
      if (PyList_Check(extra_bond_restraints_py)) {
	 int l = PyObject_Length(extra_bond_restraints_py);
	 for (int i=0; i<l; i++) {
	    PyObject *item_py = PyList_GetItem(extra_bond_restraints_py, i);
	    int item_l = PyObject_Length(item_py);
	    if (item_l == 4) {
	       coot::atom_spec_t atom_1_spec = atom_spec_from_python_expression(PyList_GetItem(item_py, 0));
	       coot::atom_spec_t atom_2_spec = atom_spec_from_python_expression(PyList_GetItem(item_py, 1));
	       double d = PyFloat_AsDouble(PyList_GetItem(item_py, 2));
	       double e = PyFloat_AsDouble(PyList_GetItem(item_py, 3));
	       coot::extra_restraints_t::extra_bond_restraint_t ebr(atom_1_spec, atom_2_spec, d, e);
	       ebr_vec.push_back(ebr);
	    }
	 }
	 int r = graphics_info_t::molecules[imol].add_extra_bond_restraints(ebr_vec);
	 graphics_draw();
      }
   }
   return ebr_vec.size();
}
#endif // USE_GUILE


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

int add_extra_target_position_restraint(int imol,
					const char *chain_id,
					int res_no,
					const char *ins_code,
					const char *atom_name,
 					const char *alt_conf, float x, float y, float z,
					float weight) {

   int r = -1;
   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t as(chain_id, res_no, ins_code, atom_name, alt_conf);
      clipper::Coord_orth pos(x,y,z);
      graphics_info_t g;
      r = g.molecules[imol].add_extra_target_position_restraint(as, pos, weight);
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
	    SCM l = scm_list_4(spec_1_scm, spec_2_scm, scm_from_double(d), scm_from_double(esd));
	    l = scm_cons(scm_string_to_symbol(scm_from_locale_string("bond")), l);
	    r = scm_cons(l, r);
	 }
         
         for (int it=g.molecules[imol].extra_restraints.angle_restraints.size()-1; it>=0; it--) {
	    coot::atom_spec_t spec_1 = g.molecules[imol].extra_restraints.angle_restraints[it].atom_1;
	    coot::atom_spec_t spec_2 = g.molecules[imol].extra_restraints.angle_restraints[it].atom_2;
	    coot::atom_spec_t spec_3 = g.molecules[imol].extra_restraints.angle_restraints[it].atom_3;
	    double a = g.molecules[imol].extra_restraints.angle_restraints[it].angle;
	    double e = g.molecules[imol].extra_restraints.angle_restraints[it].esd;
	    SCM l = scm_list_2(scm_from_double(a), scm_from_double(e));
	    l = scm_cons(atom_spec_to_scm(spec_3), l);
	    l = scm_cons(atom_spec_to_scm(spec_2), l);
	    l = scm_cons(atom_spec_to_scm(spec_1), l);
	    l = scm_cons(scm_string_to_symbol(scm_from_locale_string("angle")), l);
	    r = scm_cons(l, r);
	 }
         
	 for (int it=g.molecules[imol].extra_restraints.torsion_restraints.size()-1; it>=0; it--) {
	    coot::atom_spec_t spec_1 = g.molecules[imol].extra_restraints.torsion_restraints[it].atom_1;
	    coot::atom_spec_t spec_2 = g.molecules[imol].extra_restraints.torsion_restraints[it].atom_2;
	    coot::atom_spec_t spec_3 = g.molecules[imol].extra_restraints.torsion_restraints[it].atom_3;
	    coot::atom_spec_t spec_4 = g.molecules[imol].extra_restraints.torsion_restraints[it].atom_4;
	    double t = g.molecules[imol].extra_restraints.torsion_restraints[it].torsion_angle;
	    double e = g.molecules[imol].extra_restraints.torsion_restraints[it].esd;
	    int    p = g.molecules[imol].extra_restraints.torsion_restraints[it].period;
	    SCM l = scm_list_3(scm_from_double(t), scm_from_double(e), scm_from_int(p));
	    l = scm_cons(atom_spec_to_scm(spec_4), l);
	    l = scm_cons(atom_spec_to_scm(spec_3), l);
	    l = scm_cons(atom_spec_to_scm(spec_2), l);
	    l = scm_cons(atom_spec_to_scm(spec_1), l);
	    l = scm_cons(scm_string_to_symbol(scm_from_locale_string("torsion")), l);
	    r = scm_cons(l, r);
	 }
         
         for (int ib=g.molecules[imol].extra_restraints.start_pos_restraints.size()-1; ib>=0; ib--) {
	    coot::atom_spec_t spec_1 = g.molecules[imol].extra_restraints.start_pos_restraints[ib].atom_1;
	    double esd = g.molecules[imol].extra_restraints.start_pos_restraints[ib].esd;
	    SCM spec_1_scm = atom_spec_to_scm(spec_1);
	    SCM l = scm_list_2(spec_1_scm, scm_from_double(esd));
	    l = scm_cons(scm_string_to_symbol(scm_from_locale_string("start-pos")), l);
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
	 for (unsigned int ib=0; ib<g.molecules[imol].extra_restraints.bond_restraints.size(); ib++) {
	    coot::atom_spec_t spec_1 = g.molecules[imol].extra_restraints.bond_restraints[ib].atom_1;
	    coot::atom_spec_t spec_2 = g.molecules[imol].extra_restraints.bond_restraints[ib].atom_2;
	    double d = g.molecules[imol].extra_restraints.bond_restraints[ib].bond_dist;
	    double esd = g.molecules[imol].extra_restraints.bond_restraints[ib].esd;
	    PyObject *spec_1_py = atom_spec_to_py(spec_1);
	    PyObject *spec_2_py = atom_spec_to_py(spec_2);
	    PyObject *l = PyList_New(5);
	    PyList_SetItem(l, 0, myPyString_FromString("bond"));
	    PyList_SetItem(l, 1, spec_1_py);
	    PyList_SetItem(l, 2, spec_2_py);
	    PyList_SetItem(l, 3, PyFloat_FromDouble(d));
	    PyList_SetItem(l, 4, PyFloat_FromDouble(esd));
	    PyList_Append(r, l);
	 }
	 
	 for (unsigned int it=0; it<g.molecules[imol].extra_restraints.angle_restraints.size(); it++) {
	    coot::atom_spec_t spec_1 = g.molecules[imol].extra_restraints.angle_restraints[it].atom_1;
	    coot::atom_spec_t spec_2 = g.molecules[imol].extra_restraints.angle_restraints[it].atom_2;
	    coot::atom_spec_t spec_3 = g.molecules[imol].extra_restraints.angle_restraints[it].atom_3;
	    PyObject *spec_1_py = atom_spec_to_py(spec_1);
	    PyObject *spec_2_py = atom_spec_to_py(spec_2);
	    PyObject *spec_3_py = atom_spec_to_py(spec_3);
	    double a = g.molecules[imol].extra_restraints.angle_restraints[it].angle;
	    double e = g.molecules[imol].extra_restraints.angle_restraints[it].esd;
	    PyObject *l = PyList_New(6);
	    PyList_SetItem(l, 0, myPyString_FromString("angle"));
	    PyList_SetItem(l, 1, spec_1_py);
	    PyList_SetItem(l, 2, spec_2_py);
	    PyList_SetItem(l, 3, spec_3_py);
	    PyList_SetItem(l, 4, PyFloat_FromDouble(a));
	    PyList_SetItem(l, 5, PyFloat_FromDouble(e));
	    PyList_Append(r, l);
	 }
	 
	 for (unsigned int it=0; it<g.molecules[imol].extra_restraints.torsion_restraints.size(); it++) {
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
	    PyList_SetItem(l, 0, myPyString_FromString("torsion"));
	    PyList_SetItem(l, 1, spec_1_py);
	    PyList_SetItem(l, 2, spec_2_py);
	    PyList_SetItem(l, 3, spec_3_py);
	    PyList_SetItem(l, 4, spec_4_py);
	    PyList_SetItem(l, 5, PyFloat_FromDouble(t));
	    PyList_SetItem(l, 6, PyFloat_FromDouble(e));
	    PyList_SetItem(l, 7, PyLong_FromLong(p));
	    PyList_Append(r, l);
	 }
	 
	 for (unsigned int is=0; is<g.molecules[imol].extra_restraints.start_pos_restraints.size(); is++) {
	    coot::atom_spec_t spec_1 = g.molecules[imol].extra_restraints.start_pos_restraints[is].atom_1;
	    double esd = g.molecules[imol].extra_restraints.start_pos_restraints[is].esd;
	    PyObject *spec_1_py = atom_spec_to_py(spec_1);
	    PyObject *l = PyList_New(3);
	    PyList_SetItem(l, 0, myPyString_FromString("start pos"));
	    PyList_SetItem(l, 1, spec_1_py);
	    PyList_SetItem(l, 2, PyFloat_FromDouble(esd));
	    PyList_Append(r, l);
	 }
      }
   }
   if (PyBool_Check(r)) {
	 Py_INCREF(r);
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
      if (restraint_spec_length == 2) {
	 SCM restraint_type_scm = SCM_CAR(restraint_spec);
	 SCM spec_1_scm = scm_list_ref(restraint_spec, scm_from_int(1));
	 if (scm_is_true(scm_eq_p(restraint_type_scm, scm_string_to_symbol(scm_from_locale_string("start-pos"))))) {
	    coot::atom_spec_t spec_1 = atom_spec_from_scm_expression(spec_1_scm);
	    graphics_info_t::molecules[imol].remove_extra_start_pos_restraint(spec_1);
	    //graphics_draw(); //there is currently no graphical representation for start_pos restraints
	 }
         
      } else if (restraint_spec_length == 3) {
	 SCM restraint_type_scm = SCM_CAR(restraint_spec);
	 SCM spec_1_scm = scm_list_ref(restraint_spec, scm_from_int(1));
	 SCM spec_2_scm = scm_list_ref(restraint_spec, scm_from_int(2));
	 if (scm_is_true(scm_eq_p(restraint_type_scm, scm_string_to_symbol(scm_from_locale_string("bond"))))) {
	    coot::atom_spec_t spec_1 = atom_spec_from_scm_expression(spec_1_scm);
	    coot::atom_spec_t spec_2 = atom_spec_from_scm_expression(spec_2_scm);
	    graphics_info_t::molecules[imol].remove_extra_bond_restraint(spec_1, spec_2);
	    graphics_draw();
	 }
         
      } else if (restraint_spec_length == 4) {
	 SCM restraint_type_scm = SCM_CAR(restraint_spec);
	 SCM spec_1_scm = scm_list_ref(restraint_spec, scm_from_int(1));
	 SCM spec_2_scm = scm_list_ref(restraint_spec, scm_from_int(2));
	 SCM spec_3_scm = scm_list_ref(restraint_spec, scm_from_int(3));
	 if (scm_is_true(scm_eq_p(restraint_type_scm, scm_string_to_symbol(scm_from_locale_string("angle"))))) {
	    coot::atom_spec_t spec_1 = atom_spec_from_scm_expression(spec_1_scm);
	    coot::atom_spec_t spec_2 = atom_spec_from_scm_expression(spec_2_scm);
	    coot::atom_spec_t spec_3 = atom_spec_from_scm_expression(spec_3_scm);
	    graphics_info_t::molecules[imol].remove_extra_angle_restraint(spec_1, spec_2, spec_3);
	    //graphics_draw(); //there is currently no graphical representation for torsion restraints
	 }
      
      } else if (restraint_spec_length == 5) {
	 SCM restraint_type_scm = SCM_CAR(restraint_spec);
	 SCM spec_1_scm = scm_list_ref(restraint_spec, scm_from_int(1));
	 SCM spec_2_scm = scm_list_ref(restraint_spec, scm_from_int(2));
	 SCM spec_3_scm = scm_list_ref(restraint_spec, scm_from_int(3));
	 SCM spec_4_scm = scm_list_ref(restraint_spec, scm_from_int(4));
	 if (scm_is_true(scm_eq_p(restraint_type_scm, scm_string_to_symbol(scm_from_locale_string("torsion"))))) {
	    coot::atom_spec_t spec_1 = atom_spec_from_scm_expression(spec_1_scm);
	    coot::atom_spec_t spec_2 = atom_spec_from_scm_expression(spec_2_scm);
	    coot::atom_spec_t spec_3 = atom_spec_from_scm_expression(spec_3_scm);
	    coot::atom_spec_t spec_4 = atom_spec_from_scm_expression(spec_4_scm);
	    graphics_info_t::molecules[imol].remove_extra_torsion_restraint(spec_1, spec_2, spec_3, spec_4);
	    //graphics_draw(); //there is currently no graphical representation for torsion restraints
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
      if (restraint_spec_length == 2) {
         PyObject *restraint_type_py = PyList_GetItem(restraint_spec, 0);
         PyObject *spec_1_py = PyList_GetItem(restraint_spec, 1);
         if ((strcmp(myPyString_AsString(restraint_type_py), "start pos") == 0) ||
             (strcmp(myPyString_AsString(restraint_type_py), "start_pos") == 0) ||
             (strcmp(myPyString_AsString(restraint_type_py), "start-pos") == 0)) {
            coot::atom_spec_t spec_1 = atom_spec_from_python_expression(spec_1_py);
            graphics_info_t::molecules[imol].remove_extra_start_pos_restraint(spec_1);
            //graphics_draw(); //there is currently no graphical representation for start_pos restraints
         }
      
      } else if (restraint_spec_length == 3) {
         PyObject *restraint_type_py = PyList_GetItem(restraint_spec, 0);
         PyObject *spec_1_py = PyList_GetItem(restraint_spec, 1);
         PyObject *spec_2_py = PyList_GetItem(restraint_spec, 2);
         if (strcmp(myPyString_AsString(restraint_type_py), "bond") == 0 ) {
            coot::atom_spec_t spec_1 = atom_spec_from_python_expression(spec_1_py);
            coot::atom_spec_t spec_2 = atom_spec_from_python_expression(spec_2_py);
            graphics_info_t::molecules[imol].remove_extra_bond_restraint(spec_1, spec_2);
            graphics_draw();
         }
      
      } else if (restraint_spec_length == 4) {
         PyObject *restraint_type_py = PyList_GetItem(restraint_spec, 0);
         PyObject *spec_1_py = PyList_GetItem(restraint_spec, 1);
         PyObject *spec_2_py = PyList_GetItem(restraint_spec, 2);
         PyObject *spec_3_py = PyList_GetItem(restraint_spec, 3);
         
         if (strcmp(myPyString_AsString(restraint_type_py), "angle") == 0 ) {
            coot::atom_spec_t spec_1 = atom_spec_from_python_expression(spec_1_py);
            coot::atom_spec_t spec_2 = atom_spec_from_python_expression(spec_2_py);
            coot::atom_spec_t spec_3 = atom_spec_from_python_expression(spec_3_py);
            graphics_info_t::molecules[imol].remove_extra_angle_restraint(spec_1, spec_2, spec_3);
            //graphics_draw(); //there is currently no graphical representation for torsion restraints
         }
      
      } else if (restraint_spec_length == 5) {
         PyObject *restraint_type_py = PyList_GetItem(restraint_spec, 0);
         PyObject *spec_1_py = PyList_GetItem(restraint_spec, 1);
         PyObject *spec_2_py = PyList_GetItem(restraint_spec, 2);
         PyObject *spec_3_py = PyList_GetItem(restraint_spec, 3);
         PyObject *spec_4_py = PyList_GetItem(restraint_spec, 4);
         
         if (strcmp(myPyString_AsString(restraint_type_py), "torsion") == 0 ) {
            coot::atom_spec_t spec_1 = atom_spec_from_python_expression(spec_1_py);
            coot::atom_spec_t spec_2 = atom_spec_from_python_expression(spec_2_py);
            coot::atom_spec_t spec_3 = atom_spec_from_python_expression(spec_3_py);
            coot::atom_spec_t spec_4 = atom_spec_from_python_expression(spec_4_py);
            graphics_info_t::molecules[imol].remove_extra_torsion_restraint(spec_1, spec_2, spec_3, spec_4);
            //graphics_draw(); //there is currently no graphical representation for torsion restraints
         }
      }
   }
   
} 
#endif // USE_PYTHON

void
delete_all_extra_restraints(int imol) {

   // c.f. clear_extra_restraints()
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].clear_extra_restraints();
   }
   graphics_draw();
}


void delete_extra_restraints_for_residue(int imol, const char *chain_id, int res_no, const char *ins_code) {

   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t rs(chain_id, res_no, ins_code);
      graphics_info_t::molecules[imol].delete_extra_restraints_for_residue(rs);
   } 
   graphics_draw();
}

#ifdef USE_GUILE
void delete_extra_restraints_for_residue_spec_scm(int imol, SCM residue_spec_in) {

   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t spec = residue_spec_from_scm(residue_spec_in);
      graphics_info_t::molecules[imol].delete_extra_restraints_for_residue(spec);
   }

}
#endif // USE_GUILE

#ifdef USE_PYTHON
void delete_extra_restraints_for_residue_spec_py(int imol, PyObject *residue_spec_in_py) {

   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t spec = residue_spec_from_py(residue_spec_in_py);
      graphics_info_t::molecules[imol].delete_extra_restraints_for_residue(spec);
   }
}
#endif // USE_PYTHON


void delete_extra_restraints_worse_than(int imol, float n_sigma) { 

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].delete_extra_restraints_worse_than(n_sigma);
   }
   graphics_draw();
} 


void
set_use_only_extra_torsion_restraints_for_torsions(short int state) {
   graphics_info_t:: use_only_extra_torsion_restraints_for_torsions_flag = state;
} 

void
add_initial_position_restraints(int imol, const std::vector<coot::residue_spec_t> &residue_specs,
				double weight) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      for (unsigned int i=0; i<residue_specs.size(); i++) {
	 mmdb::Residue *residue_p = g.molecules[imol].get_residue(residue_specs[i]);
	 if (residue_p) {
	    mmdb::PPAtom residue_atoms = 0;
	    int n_residue_atoms;
	    residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
	    for (int iat=0; iat<n_residue_atoms; iat++) {
	       mmdb::Atom *at = residue_atoms[iat];
	       add_extra_start_pos_restraint(imol,
					     at->GetChainID(), 
					     at->GetSeqNum(), 
					     at->GetInsCode(), 
					     at->GetAtomName(), 
					     at->altLoc,
					     weight);
	    }
	 } 
      }
   }
}

void
remove_initial_position_restraints(int imol, const std::vector<coot::residue_spec_t> &residue_specs) {
   delete_all_extra_restraints(imol);
}



void set_show_intermediate_atoms_rota_markup(short int state) {
   graphics_info_t::do_intermediate_atoms_rota_markup = state;
}


void set_show_intermediate_atoms_rama_markup(short int state) {
   graphics_info_t::do_intermediate_atoms_rama_markup = state;

}

#ifdef USE_PYTHON
void register_post_intermediate_atoms_moved_hook(PyObject *function) {

   graphics_info_t g;
   g.register_post_intermediate_atoms_moved_hook(function);

}
#endif


// trash the multimodal (sp3) ring torsions and use
// only unimodal restraints
void use_unimodal_ring_torsion_restraints(const std::string &res_name) {

   // uses auto-load if not already present in the store

   bool minimal = false; // don't allow minimal
   int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
   graphics_info_t::Geom_p()->use_unimodal_ring_torsion_restraints(imol_enc, res_name, minimal);

}


#ifdef USE_PYTHON
//! \brief use unimodal ring torsion restraints (e.g. for carbohydrate pyranose)
//
//         allow user definition of torsions for given residue
//  @var{torsions_info_list} is a list of item that are of the form
//  @var{[atom_name_1, atom_name_2, atom_name_3, atom_name_4, double torsion_1234]}
void use_unimodal_ring_torsion_restraints_for_residue(const std::string &res_name, PyObject *torsions_info_list) {

   if (PyList_Check(torsions_info_list)) {
      unsigned int n_torsions = PyObject_Length(torsions_info_list);
      std::vector<coot::atom_name_torsion_quad> tors_info_vec;
      for (unsigned int i=0; i<n_torsions; i++) {
	 PyObject *tors_info = PyList_GetItem(torsions_info_list, i);
	 if (PyList_Check(tors_info)) {
	    unsigned int n_eles = PyObject_Length(tors_info);
	    if (n_eles == 5) {
	       PyObject *at_1_py = PyList_GetItem(tors_info, 0);
	       PyObject *at_2_py = PyList_GetItem(tors_info, 1);
	       PyObject *at_3_py = PyList_GetItem(tors_info, 2);
	       PyObject *at_4_py = PyList_GetItem(tors_info, 3);
	       PyObject *tors_py = PyList_GetItem(tors_info, 4);
	       if (PyUnicode_Check(at_1_py)) {
		  if (PyUnicode_Check(at_2_py)) {
		     if (PyUnicode_Check(at_3_py)) {
			if (PyUnicode_Check(at_4_py)) {
			   if (PyFloat_Check(tors_py)) {
                              std::string at_name_1 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(at_1_py));
                              std::string at_name_2 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(at_2_py));
                              std::string at_name_3 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(at_3_py));
                              std::string at_name_4 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(at_4_py));
			      double tors = PyFloat_AsDouble(tors_py);
			      std::string id = "ring-torsion-";
			      id += coot::util::int_to_string(tors_info_vec.size()+1);
			      coot::atom_name_torsion_quad tors_info(id, at_name_1, at_name_2, at_name_3, at_name_4, tors);
			      tors_info_vec.push_back(tors_info);
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }

      if (tors_info_vec.size() > 0) {
	 int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
	 graphics_info_t::Geom_p()->use_unimodal_ring_torsion_restraints(imol_enc, res_name, tors_info_vec, graphics_info_t::cif_dictionary_read_number);
      }
   }
}
#endif


void set_refinement_geman_mcclure_alpha(float alpha) {

   graphics_info_t g;
   g.set_geman_mcclure_alpha(alpha);

}

//! \brief get the Geman-McClure distance alpha value (weight)
float get_refinement_geman_mcclure_alpha() {

   return graphics_info_t::geman_mcclure_alpha;
}


//! \brief set the Lennard Jones epsilon parameter
void set_refinement_lennard_jones_epsilon(float epsilon) {

   graphics_info_t g;
   g.set_lennard_jones_epsilon(epsilon);

}

//! \brief set the log cosh scale factor for target position restraints
void set_log_cosh_target_distance_scale_factor(float sf) {
   graphics_info_t::log_cosh_target_distance_scale_factor = sf;
}



#ifdef USE_GUILE
void crankshaft_peptide_rotation_optimization_scm(int imol, SCM residue_spec_scm) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      coot::residue_spec_t rs = residue_spec_from_scm(residue_spec_scm);
      unsigned int n_peptides = 3;
      int n_samples = -1; // auto

      int imol_map = g.Imol_Refinement_Map();
      if (is_valid_map_molecule(imol_map)) {
	 const clipper::Xmap<float> &xmap = g.molecules[imol_map].xmap;
	 float w = g.geometry_vs_map_weight;
	 int n_threads = coot::get_max_number_of_threads() - 1;
	 if (n_threads < 1) n_threads = 1;

	 g.molecules[imol].crankshaft_peptide_rotation_optimization(rs, n_peptides, xmap, w, n_samples,
								    &g.static_thread_pool, n_threads);
	 g.update_validation_graphs(imol);
	 graphics_draw();
      }
   }

}
#endif

#ifdef USE_PYTHON
void crankshaft_peptide_rotation_optimization_py(int imol, PyObject *residue_spec_py) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      coot::residue_spec_t rs = residue_spec_from_py(residue_spec_py);
      unsigned int n_peptides = 3;
      int n_samples = -1; // auto

      int imol_map = g.Imol_Refinement_Map();
      if (is_valid_map_molecule(imol_map)) {
	 const clipper::Xmap<float> &xmap = g.molecules[imol_map].xmap;
	 float w = g.geometry_vs_map_weight;
	 int n_threads = coot::get_max_number_of_threads() - 1;
	 if (n_threads < 1) n_threads = 1;
	 g.molecules[imol].crankshaft_peptide_rotation_optimization(rs, n_peptides, xmap, w, n_samples,
								    &g.static_thread_pool, n_threads);
	 g.update_validation_graphs(imol);
	 graphics_draw();
      }
   }
}
#endif


int crankshaft_peptide_rotation_optimization_intermediate_atoms() {

   graphics_info_t g;
   return g.crankshaft_peptide_rotation_optimization_intermediate_atoms();

}


void set_regenerate_bonds_needs_make_bonds_type_checked(bool state) {
   graphics_info_t g;
   g.set_regenerate_bonds_needs_make_bonds_type_checked(state);
}

bool get_regenerate_bonds_needs_make_bonds_type_checked_state() {
   graphics_info_t g;
   return g.get_regenerate_bonds_needs_make_bonds_type_checked_state();
}

/*! \brief shiftfield B-factor refinement */
void shiftfield_b_factor_refinement(int imol) {
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      g.shiftfield_b_factor_refinement(imol);
   }
}

/*! \brief shiftfield xyz refinement */
void shiftfield_xyz_factor_refinement(int imol) {
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      g.shiftfield_xyz_factor_refinement(imol);
   }
}

void convert_dictionary_planes_to_improper_dihedrals() {
   graphics_info_t g;
   g.Geom_p()->all_plane_restraints_to_improper_dihedrals();
   g.Geom_p()->delete_plane_restraints();
   g.set_convert_dictionary_planes_to_improper_dihedrals(true);

}

