/* src/graphics-info-defines.cc
 * 
 * Copyright 2004, 2005 by The University of York
 * Copyright 2007 by The University of York
 * Copyright 2008, 2009 by The University of Oxford
 * Author: Paul Emsley
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1335, USA
 */

#if defined (USE_PYTHON)
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"


#if defined _MSC_VER
#include <windows.h>
#endif

#include <algorithm>

#include "graphics-info.h"
#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "interface.h"

#include "rotate-translate-modes.hh"
#include "manipulation-modes.hh"

#include "coot-fileselections.h"

void
graphics_info_t::clear_pending_picks() {
   a_is_pressed = 0;
   in_range_define = 0;
   in_range_define_for_refine = 0;
   in_pepflip_define = 0;
   in_rigid_body_define = 0;
   in_terminal_residue_define = 0;
   in_rot_trans_object_define = 0;
   in_residue_info_define = 0;
   in_distance_define = 0;
   in_angle_define = 0;
   in_torsion_define = 0;
   in_rotamer_define = 0;
   in_mutate_define = 0;
   in_mutate_auto_fit_define = 0;
   in_auto_fit_define = 0;
   in_db_main_define           = 0;
   in_edit_phi_psi_define      = 0;
   in_add_alt_conf_define      = 0;
   in_save_symmetry_define     = 0;
   in_cis_trans_convert_define = 0;
   in_180_degree_flip_define   = 0;
   in_reverse_direction_define = 0;
   in_dynamic_distance_define  = 0;
   in_torsion_general_define   = 0;
   pick_pending_flag           = 0;
   in_user_defined_define      = 0;
   in_multi_residue_torsion_define = 0; 
   in_edit_chi_angles_define   = 0; // Added 20100815: is there a reason why this was missing?

   std::vector<std::string> button_name_vec =
      model_fit_refine_toggle_button_name_list();
   for (unsigned int i=0; i<button_name_vec.size(); i++)
      model_fit_refine_unactive_togglebutton(button_name_vec[i]);

   button_name_vec = other_modelling_tools_toggle_button_name_list();
   for (unsigned int i=0; i<button_name_vec.size(); i++)
      other_modelling_tools_unactive_togglebutton(button_name_vec[i]);

   normal_cursor();
   std::cout << "Pending Picks Cleared" << std::endl;
}

// 20051004 now use use this interface rather direct access to pick_pending_flag
//          because a overal pick_pending_flag doesn't work when we consider turning
//          off various different model/fit/refine toggle-buttons.  So
//          static_graphics_pick_pending() checks each of the pending picks which
//          would make the graphics not rotate (ctrl key issues).
short int
graphics_info_t::static_graphics_pick_pending() const {

   // Maybe more stuff needs to go here.  I can't think at the moment.
   // in_range_define is regularization
   return in_range_define_for_refine || in_range_define;

}


//static
std::vector<std::string>
graphics_info_t::model_fit_refine_toggle_button_name_list() {


   std::vector<std::string> names;
   names.push_back("model_refine_dialog_regularize_zone_togglebutton");
   names.push_back("model_refine_dialog_refine_togglebutton");
   names.push_back("model_refine_dialog_pepflip_togglebutton");
   names.push_back("model_refine_dialog_rigid_body_togglebutton");
   names.push_back("model_refine_dialog_fit_terminal_residue_togglebutton");
   names.push_back("model_refine_dialog_rot_trans_togglebutton");
   names.push_back("model_refine_dialog_rotamer_togglebutton");
   names.push_back("model_refine_dialog_mutate_togglebutton");
   names.push_back("model_refine_dialog_mutate_auto_fit_togglebutton");
   names.push_back("model_refine_dialog_auto_fit_rotamer_togglebutton");
   names.push_back("model_refine_dialog_edit_phi_psi_togglebutton");
   names.push_back("model_refine_dialog_edit_chi_angles_togglebutton");
   names.push_back("model_refine_dialog_torsion_general_togglebutton");
   names.push_back("model_refine_dialog_do_180_degree_sidechain_flip_togglebutton");
   names.push_back("model_refine_dialog_edit_backbone_torsions_togglebutton");
   return names;
}

//static
std::vector<std::string>
graphics_info_t::other_modelling_tools_toggle_button_name_list() {

   std::vector<std::string> names;
   names.push_back("cis_trans_conversion_toggle_button");
   names.push_back("model_refine_dialog_db_main_togglebutton");

   return names;
}




// static
std::vector<std::string>
graphics_info_t::model_fit_refine_button_name_list() {

   std::vector<std::string> names;
   names.push_back("model_refine_dialog_refine_params_button");
   names.push_back("model_refine_dialog_map_select_button");
#if (GTK_MAJOR_VERSION >1)
   names.push_back("model_refine_dialog_fixed_atoms_button");
#else
   names.push_back("model_refine_dialog_fix_atoms_button");
#endif
//   names.push_back("model_refine_dialog_find_waters_button");
   names.push_back("model_refine_dialog_add_alt_conf_button");
   names.push_back("model_refine_dialog_pointer_atom_button");
   names.push_back("model_refine_dialog_clear_pending_button");
   names.push_back("model_refine_dialog_delete_button");
   names.push_back("model_refine_dialog_undo_button");
   names.push_back("model_refine_dialog_refmac_button");
   return names;
}

//static
std::vector<std::string>
graphics_info_t::other_modelling_tools_button_name_list() {

   std::vector<std::string> names;
   names.push_back("model_refine_dialog_find_waters_button"); 
   names.push_back("model_refine_dialog_find_ligands_button");
   names.push_back("model_refine_dialog_fast_sss_button"); 
   names.push_back("model_refine_dialog_baton_button"); 
   names.push_back("model_refine_dialog_add_OXT_button");
   names.push_back("place_helix_here_button");
   return names;
}


// gtk_widget_set_name, so that we can find the widget name in the
// resources file.
// 
// static
void
graphics_info_t::set_model_fit_refine_button_names(GtkWidget *widget) { 

   std::vector<std::string> toggle_button_names =
      model_fit_refine_toggle_button_name_list();

   std::vector<std::string> normal_button_names =
      model_fit_refine_button_name_list();

   std::vector<std::string> button_names = toggle_button_names;
   for (unsigned int i=0; i<normal_button_names.size(); i++)
      button_names.push_back(normal_button_names[i]);

   for (unsigned int i=0; i<button_names.size(); i++) {
      GtkWidget *w = lookup_widget(widget, button_names[i].c_str());
      if (w) {
	 gtk_widget_set_name(w, button_names[i].c_str());
      }
   }
}

// and for the other modelling toolbar
void
graphics_info_t::set_other_modelling_tools_button_names(GtkWidget *widget) { 

   std::vector<std::string> other_button_names =
     other_modelling_tools_button_name_list();

   std::vector<std::string> button_names = other_button_names;
   // dont need extra ones here yet
   //for (unsigned int i=0; i<other_button_names.size(); i++)
   //   button_names.push_back(other_button_names[i]);

   for (unsigned int i=0; i<button_names.size(); i++) {
      GtkWidget *w = lookup_widget(widget, button_names[i].c_str());
      if (w) {
        gtk_widget_set_name(w, button_names[i].c_str());
      }
   }
}

   


// static
void
graphics_info_t::untoggle_model_fit_refine_buttons_except(const std::string &button_name) {

   std::vector<std::string> button_name_vec =
      graphics_info_t::model_fit_refine_toggle_button_name_list();
   for (unsigned int i=0; i<button_name_vec.size(); i++)
      if (button_name_vec[i] != button_name)
	 model_fit_refine_unactive_togglebutton(button_name_vec[i]);
}


int
graphics_info_t::check_if_in_range_defines(GdkEventButton *event,
					   const GdkModifierType &state) {

   int iv = 0;

   check_if_in_user_defined_define(event);
   check_if_in_residue_info_define(event);
   iv += check_if_in_regularize_define(event);
   iv += check_if_in_refine_define(event);
   check_if_in_rot_trans_define(event); // iv += ? not currently.
   check_if_in_rigid_body_define(event);

   check_if_in_180_degree_flip_define(event);
   check_if_in_geometry_range_defines(event);
   check_if_in_pepflip_define(event);
   check_if_in_terminal_residue_define(event);
   check_if_in_delete_item_define(event, state);
   check_if_in_rotamer_define(event);
   check_if_in_mutate_define(event); 
   check_if_in_mutate_auto_fit_define(event);
   check_if_in_auto_fit_define(event);
   check_if_in_add_alt_conf_define(event);
   check_if_in_edit_phi_psi_define(event);
   check_if_in_edit_chi_angles_define(event);
   check_if_in_edit_backbone_torsion_define(event);
   check_if_in_save_symmetry_define(event);
   check_if_in_cis_trans_convertion_define(event);
   check_if_in_db_main_define(event);
   check_if_in_reverse_direction_define(event);
   check_if_in_lsq_plane_define(event);
   check_if_in_lsq_plane_deviant_atom_define(event);
   check_if_in_torsion_general_define(event);
   check_if_in_residue_partial_alt_locs(event);
   check_if_in_fixed_atom_define(event, state);
   check_if_in_base_pairing_define(event);
   check_if_in_multi_residue_torsion_define(event);

   return iv;
}

void
graphics_info_t::check_if_in_user_defined_define(GdkEventButton *event) {

   graphics_info_t g;
   if (g.in_user_defined_define) {
      pick_info nearest_atom_index_info = atom_pick(event);
      if (nearest_atom_index_info.success == GL_TRUE) {
	 in_user_defined_define--;
	 int im = nearest_atom_index_info.imol;
	 molecules[im].add_to_labelled_atom_list(nearest_atom_index_info.atom_index);
	 mmdb::Atom *at = molecules[im].atom_sel.atom_selection[nearest_atom_index_info.atom_index];
	 if (at) { 
	    coot::atom_spec_t spec(at);
	    spec.int_user_data = im;
	    user_defined_atom_pick_specs.push_back(spec);
	    graphics_draw(); // let's see the label
	    if (in_user_defined_define == 0) {
	       run_user_defined_click_func(); // uses user_defined_atom_pick_specs
	    }
	 }
	 normal_cursor();
      }
   } 
} 


void
graphics_info_t::check_if_in_residue_info_define(GdkEventButton *event) {

   // Eleanor's Residue Info: 
   //
   graphics_info_t info;
   if (info.in_residue_info_define == 1) { 
      pick_info nearest_atom_index_info = atom_pick(event);
      if (nearest_atom_index_info.success == GL_TRUE) { 
	 int im = nearest_atom_index_info.imol; 
	 // std::cout << "info: clicked on imol: " << im << std::endl;
	 // a c-interface-info function...
	 output_residue_info_dialog(im, nearest_atom_index_info.atom_index); 
	 info.in_residue_info_define = 0;
	 normal_cursor();
	 pick_pending_flag = 0;
      }
   }
}


int
graphics_info_t::check_if_in_refine_define(GdkEventButton *event) { 

   int iv = 0;

   if (in_range_define_for_refine) { 

      iv = 1;

      //       int auto_range_flag = 0;
      pick_info naii = atom_pick(event);
      if (naii.success == GL_TRUE) { 
	 molecules[naii.imol].add_to_labelled_atom_list(naii.atom_index);

	 if (in_range_define_for_refine == 1) { 
	    
	    residue_range_atom_index_1 = naii.atom_index;
	    residue_range_mol_no = naii.imol;
	    in_range_define_for_refine = 2;
	    // removed autorange code, on suggestion of Kevin.  Put
	    // into A key callback now.

	 } else { 

	    // (in_range_define_for_refine == 2)

	    if (naii.imol == residue_range_mol_no) {
	       watch_cursor();
	       residue_range_atom_index_2 = naii.atom_index;
	       int auto_range_flag = 0;
	       rot_trans_rotation_origin_atom = 0; // flag for Ctrl left
						   // mouse behaviour (we
						   // don't want to rotate
						   // the atoms)
	       refine(residue_range_mol_no,
		      auto_range_flag,
		      residue_range_atom_index_1,
		      residue_range_atom_index_2);
	    }
	    
	    in_range_define_for_refine = 0;
	    normal_cursor();
	    pick_pending_flag = 0;
	    model_fit_refine_unactive_togglebutton("model_refine_dialog_refine_togglebutton");
	 }
	 graphics_draw(); // let's see the label
      }
   }
   return iv;
}

int 
graphics_info_t::check_if_in_regularize_define(GdkEventButton *event) { 

   int iv = 0;

   if (in_range_define) { // regularization

      iv = 1;
      
      pick_info naii = atom_pick(event);
      int auto_range_flag = 0;
      if (naii.success == GL_TRUE) { 
	 int imol = naii.imol;
	 molecules[naii.imol].add_to_labelled_atom_list(naii.atom_index);
	 if (in_range_define == 1) { 
	    residue_range_atom_index_1 = naii.atom_index;
	    residue_range_mol_no = imol;
	    in_range_define = 2;
	    if (a_is_pressed) { 
	       auto_range_flag = 1;
	       rot_trans_rotation_origin_atom = 0; // flag for Ctrl left
						   // mouse behaviour (we
						   // don't want to rotate
						   // the atoms)
	       regularize(residue_range_mol_no,
			  auto_range_flag,
			  residue_range_atom_index_1,
			  residue_range_atom_index_1);
	       in_range_define = 0;
	       normal_cursor();
	       pick_pending_flag = 0;
	       model_fit_refine_unactive_togglebutton("model_refine_dialog_regularize_zone_togglebutton");

	    }
	 } else { 
	    // 
	    if (naii.imol == residue_range_mol_no) { 
	       watch_cursor();
	       residue_range_atom_index_2 = naii.atom_index;
	       auto_range_flag = 0;
	       rot_trans_rotation_origin_atom = 0; // flag for Ctrl left
						   // mouse behaviour (we
						   // don't want to rotate
						   // the atoms)
	       regularize(residue_range_mol_no,
			  auto_range_flag,
			  residue_range_atom_index_1,
			  residue_range_atom_index_2);
	    }
	    in_range_define = 0;
	    normal_cursor();
	    pick_pending_flag = 0;
	    model_fit_refine_unactive_togglebutton("model_refine_dialog_regularize_zone_togglebutton");
	 }
	 graphics_draw(); // let's see the label

      } 
   }
   return iv;
}

// distances/angles
void
graphics_info_t::check_if_in_geometry_range_defines(GdkEventButton *event) {

   if (in_distance_define) {
      
      pick_info nearest_atom_index_info; 
      nearest_atom_index_info = atom_pick(event);
	    
      if ( nearest_atom_index_info.success == GL_TRUE ) {

	 int im = nearest_atom_index_info.imol; 
	 std::cout << "geometry: on molecule number: " << im << std::endl;
	 // some visual feedback, label the atom:
	 molecules[im].add_to_labelled_atom_list(nearest_atom_index_info.atom_index);

	 if (in_distance_define == 1) {
	    in_distance_define = 2; // flag for next atom pick
	    geometry_atom_index_1 = nearest_atom_index_info.atom_index;
	    geometry_atom_index_1_mol_no = nearest_atom_index_info.imol;
	    mmdb::Atom *atom1 = molecules[im].atom_sel.atom_selection[geometry_atom_index_1];
	    distance_pos_1 = coot::Cartesian(atom1->x, atom1->y, atom1->z);
	    std::cout << "click on a second atom" << std::endl;
	    graphics_draw();
	 } else {

	    // in_distance_define == 2
	    geometry_atom_index_2 =
	       nearest_atom_index_info.atom_index;
	    geometry_atom_index_2_mol_no =
	       nearest_atom_index_info.imol;

	    mmdb::Atom *atom2 = molecules[im].atom_sel.atom_selection[geometry_atom_index_2];
	    coot::Cartesian pos2 = coot::Cartesian(atom2->x, atom2->y, atom2->z);

	    // 20190104-PE Why were we using the symmetry function?
// 	    display_geometry_distance_symm(geometry_atom_index_1_mol_no, distance_pos_1,
// 					   geometry_atom_index_2_mol_no, pos2);

	    display_geometry_distance(geometry_atom_index_1_mol_no, distance_pos_1,
				      geometry_atom_index_2_mol_no, pos2); // calls graphics_draw()

	    unset_geometry_dialog_distance_togglebutton();
	    in_distance_define = 0;  // clear flag
	    pick_pending_flag = 0;
	    normal_cursor();
	 }
      } else { 
	 // let's try symmetry pick:
	 // 
	 coot::Symm_Atom_Pick_Info_t symm_nearest_atom_index_info = symmetry_atom_pick(); 

	 if (symm_nearest_atom_index_info.success == GL_TRUE) {

	    int im = symm_nearest_atom_index_info.imol; 
	    // some visual feedback, label the atom:
	    molecules[im].add_atom_to_labelled_symm_atom_list(symm_nearest_atom_index_info.atom_index,
							      symm_nearest_atom_index_info.symm_trans,
							      symm_nearest_atom_index_info.pre_shift_to_origin);


	    if (in_distance_define == 1) {
	       distance_pos_1 = symm_nearest_atom_index_info.hybrid_atom.pos;
	       geometry_atom_index_1_mol_no = symm_nearest_atom_index_info.imol;
	       in_distance_define = 2;
	       std::cout << "click on a second atom" << std::endl;
	       graphics_draw();

	    } else { 

	       // in_distance_define == 2
	       coot::Cartesian pos2 = symm_nearest_atom_index_info.hybrid_atom.pos;
	       geometry_atom_index_2_mol_no = symm_nearest_atom_index_info.imol;
	       display_geometry_distance(geometry_atom_index_1_mol_no, distance_pos_1,
					 geometry_atom_index_2_mol_no, pos2);
	       unset_geometry_dialog_distance_togglebutton();
	       in_distance_define = 0;
	       pick_pending_flag = 0;
	       normal_cursor();
	    }
	 }
      }
   }

   // Angle define
   if (in_angle_define) {
      
      // We need a cleaner way to know if this was an atom pick or a symm atom pick.
      // Let's sort it out in the beginning:
      short int picked = 0; 
      pick_info nearest_atom_index_info = atom_pick(event);
      mmdb::Atom *atom = 0;
      coot::Cartesian pos;
      
      if (nearest_atom_index_info.success == GL_TRUE) {
	 picked = 1;
	 int im = nearest_atom_index_info.imol;
	 atom = molecules[im].atom_sel.atom_selection[nearest_atom_index_info.atom_index];
	 molecules[im].add_to_labelled_atom_list(nearest_atom_index_info.atom_index);
	 pos = coot::Cartesian(atom->x, atom->y, atom->z);
      } else {
	 coot::Symm_Atom_Pick_Info_t symm_nearest_atom_index_info = symmetry_atom_pick(); 
	 if (symm_nearest_atom_index_info.success == GL_TRUE) {
	    picked = 1;
	    int im = symm_nearest_atom_index_info.imol; 
	    // some visual feedback, label the atom:
	    molecules[im].add_atom_to_labelled_symm_atom_list(symm_nearest_atom_index_info.atom_index,
							      symm_nearest_atom_index_info.symm_trans,
							      symm_nearest_atom_index_info.pre_shift_to_origin);
	    pos = symm_nearest_atom_index_info.hybrid_atom.pos;
	 }
      }

      if (picked) {
	       
	 if (in_angle_define == 1) {
	    in_angle_define = 2; // flag for next atom pick
	    angle_tor_pos_1 = pos;
	    graphics_draw();
	    
	 } else {

	    if (in_angle_define == 2) {
	       in_angle_define = 3; // flag for next atom pick
	       angle_tor_pos_2 = pos;
	       graphics_draw();

	    } else { 
	       // in_angle_define == 3
	       angle_tor_pos_3 = pos;
	       graphics_draw();
	       
	       display_geometry_angle(); // uses class members
	                                 // that we have just set

	       in_angle_define = 0;  // clear flag
	       pick_pending_flag = 0;
	       normal_cursor();
	       unset_geometry_dialog_angle_togglebutton();
	    }
	 }	 
	 graphics_draw();
      }
   }

   //
   if (in_torsion_define) {

      // We need a cleaner way to know if this was an atom pick or a symm atom pick.
      // Let's sort it out in the beginning:
      short int picked = 0; 
      pick_info nearest_atom_index_info = atom_pick(event);
      mmdb::Atom *atom = 0;
      coot::Cartesian pos;
      
      if (nearest_atom_index_info.success == GL_TRUE) {
	 picked = 1;
	 int im = nearest_atom_index_info.imol;
	 atom = molecules[im].atom_sel.atom_selection[nearest_atom_index_info.atom_index];
	 molecules[im].add_to_labelled_atom_list(nearest_atom_index_info.atom_index);
	 pos = coot::Cartesian(atom->x, atom->y, atom->z);
      } else {
	 coot::Symm_Atom_Pick_Info_t symm_nearest_atom_index_info = symmetry_atom_pick(); 
	 if (symm_nearest_atom_index_info.success == GL_TRUE) {
	    picked = 1;
	    int im = symm_nearest_atom_index_info.imol; 
	    // some visual feedback, label the atom:
	    molecules[im].add_atom_to_labelled_symm_atom_list(symm_nearest_atom_index_info.atom_index,
							      symm_nearest_atom_index_info.symm_trans,
							      symm_nearest_atom_index_info.pre_shift_to_origin);
	    pos = symm_nearest_atom_index_info.hybrid_atom.pos;
	 }
      }

      if (picked) {

	 if (in_torsion_define == 1) {
	    angle_tor_pos_1 = pos;
	    in_torsion_define = 2; // flag for next atom pick
	    graphics_draw();
	 } else {
	    if (in_torsion_define == 2) {
	       angle_tor_pos_2 = pos;
	       in_torsion_define = 3; // flag for next atom pick
	       graphics_draw();
	    } else {
	       if (in_torsion_define == 3) {
		  angle_tor_pos_3 = pos;
		  in_torsion_define = 4; // flag for next atom pick
		  graphics_draw();
	       } else {
		  // in_torsion_define == 4
		  angle_tor_pos_4 = pos;
		  display_geometry_torsion(); // does a draw
		  in_torsion_define = 0; // clear up.
		  pick_pending_flag = 0;
		  normal_cursor();
		  unset_geometry_dialog_torsion_togglebutton(); 
	       }
	    }
	 }
      }
   }

   if (in_dynamic_distance_define) {
      if (! moving_atoms_asc) {
	 std::cout << "No intermediate atoms available" << std::endl;
	 add_status_bar_text("No intermediate atoms available");
      } else {
	 std::cout << "in_dynamic_distance_define check "
		   << in_dynamic_distance_define << std::endl;
	 if (in_dynamic_distance_define == 1) { 
	    pick_info nearest_atom_index_info = atom_pick(event);
	    pick_info nearest_intermediate_atom_info = pick_intermediate_atom(*moving_atoms_asc);
	    bool do_static_atom = 0;
	    bool do_intermediate_atom = 0;
	    if (nearest_atom_index_info.success)
	       if (nearest_intermediate_atom_info.success)
		  if (nearest_atom_index_info.min_dist < nearest_intermediate_atom_info.min_dist)
		     do_static_atom = 1;
	       else
		  do_intermediate_atom = 1;
	    else
	       do_static_atom = 1;
	    else
	       if (nearest_intermediate_atom_info.success)
		  do_intermediate_atom = 1;
	 
	    mmdb::Atom *atom = 0;
	    coot::Cartesian pos;
	    if (do_static_atom) {
	       std::cout << "in_dynamic_distance_define pick 1 static atom "
			 << in_dynamic_distance_define << std::endl;
	       int im = nearest_atom_index_info.imol;
	       atom = molecules[im].atom_sel.atom_selection[nearest_atom_index_info.atom_index];
	       molecules[im].add_to_labelled_atom_list(nearest_atom_index_info.atom_index);
	       pos = coot::Cartesian(atom->x, atom->y, atom->z);
	       running_dynamic_distance = coot::intermediate_atom_distance_t(pos);
	       in_dynamic_distance_define = 2;
	       graphics_draw();
	    } else {
	       if (do_intermediate_atom) { 
		  std::cout << "in_dynamic_distance_define pick 1 dynamic atom "
			    << in_dynamic_distance_define << std::endl;
		  int ai = nearest_intermediate_atom_info.atom_index;
		  mmdb::Atom *at = moving_atoms_asc->atom_selection[ai];
		  running_dynamic_distance = coot::intermediate_atom_distance_t(at);
		  add_status_bar_text("Now click on a molecule atom");
		  in_dynamic_distance_define = 2; 
		  graphics_draw();
	       }
	    }
	 } else {

	    if (in_dynamic_distance_define == 2) {
	       if (running_dynamic_distance.atom_is_filled()) {
		  pick_info naii = atom_pick(event);
		  if (naii.success == GL_TRUE) {
		     std::cout << "in_dynamic_distance_define pick 2 static atom "
			       << in_dynamic_distance_define << std::endl;
		     int im = naii.imol;
		     mmdb::Atom *at = molecules[im].atom_sel.atom_selection[naii.atom_index];
		     coot::Cartesian pt(at->x, at->y, at->z);
		     running_dynamic_distance.add_static_point(pt);
		  }
	       }
	       
	       if (running_dynamic_distance.static_position_is_filled()) {
		  pick_info niaii = pick_intermediate_atom(*moving_atoms_asc);
		  if (niaii.success == GL_TRUE) {
		     std::cout << "in_dynamic_distance_define pick 2 dynamic atom "
			       << in_dynamic_distance_define << std::endl;
		     mmdb::Atom *at = moving_atoms_asc->atom_selection[niaii.atom_index];
		     running_dynamic_distance.add_atom(at);
		  }
	       }
	    }
	 }

	 if (running_dynamic_distance.filled()) {
	    dynamic_distances.push_back(running_dynamic_distance);
	    graphics_draw();
	    normal_cursor();
	    in_dynamic_distance_define = 0;
	    running_dynamic_distance = coot::intermediate_atom_distance_t();
	    unset_geometry_dialog_dynamic_distance_togglebutton();
	 }
      }
   }
}


void
graphics_info_t::check_if_in_pepflip_define(GdkEventButton *event) {

   if (in_pepflip_define == 1) {
      pick_info nearest_atom_index_info = atom_pick(event);
      if (nearest_atom_index_info.success == GL_TRUE) {
	 imol_pepflip = nearest_atom_index_info.imol;
	 atom_index_pepflip = nearest_atom_index_info.atom_index;

	 pepflip();
	 in_pepflip_define = 0;
	 pick_pending_flag = 0;
	 normal_cursor();
	 model_fit_refine_unactive_togglebutton("model_refine_dialog_pepflip_togglebutton");
	 graphics_draw();
      }
   } 
}


void
graphics_info_t::check_if_in_delete_item_define(GdkEventButton *event,
						const GdkModifierType &state) {

   if (false)
      std::cout << "DEBUG:: delete_item atom " << delete_item_atom
		<< " residue " << delete_item_residue
		<< " water " << delete_item_water
		<< " sidechain " << delete_item_sidechain
		<< " sidechain-range " << delete_item_sidechain_range
		<< " chain " << delete_item_chain
		<< " range " << delete_item_residue_zone
		<< " hydrogens " << delete_item_residue_hydrogens
		<< std::endl;

   bool item_deleted = false;
   int imol_delete = -1;
   graphics_info_t g;
   short int destroy_delete_dialog_flag_by_ctrl_press = 1;
   if (state & GDK_CONTROL_MASK)
      destroy_delete_dialog_flag_by_ctrl_press = 0;
   
   if (g.delete_item_widget) { 
      // atom
      if (g.delete_item_atom) {
	 pick_info naii = atom_pick(event);
	 if (naii.success == GL_TRUE) {
	    mmdb::Residue *res = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index]->residue;
	    std::string resname(res->name);
	    if (resname == "WAT" || resname == "HOH") {
	       std::string w = "Protection measure: ";
	       w += "Waters can only be deleted using Delete Water";
	       add_status_bar_text(w);
	    } else {
	       normal_cursor();
	       delete_atom_by_atom_index(naii.imol, naii.atom_index,
					 destroy_delete_dialog_flag_by_ctrl_press);
	       run_post_manipulation_hook(naii.imol, DELETED);
	       pick_pending_flag = 0;
	       item_deleted = true;
	       imol_delete = naii.imol;
	    }
	 } else { 

	    if (show_symmetry) {
	       coot::Symm_Atom_Pick_Info_t symm_nearest_atom_index_info = symmetry_atom_pick(); 

	       if (1) {

		  if (symm_nearest_atom_index_info.success == GL_TRUE) {
		     int im = symm_nearest_atom_index_info.imol;
		     int index = symm_nearest_atom_index_info.atom_index;
		     delete_atom_by_atom_index(im, index,
					       destroy_delete_dialog_flag_by_ctrl_press);
		     normal_cursor();
		     run_post_manipulation_hook(naii.imol, DELETED);
		     pick_pending_flag = 0;
		     item_deleted = true;
		     imol_delete = im;
		  }
		  
	       } else {

		  // old code - WARNING, not executed
		  if (symm_nearest_atom_index_info.success == GL_TRUE) {

		     std::string s = "That was a symmetry atom\n";
		     s += "Coot currently doesn't delete symmetry atoms";
		     GtkWidget *w = wrapped_nothing_bad_dialog(s);
		     gtk_widget_show(w);
		  }
	       }
	    }
	 }
      }

      // water
      if (g.delete_item_water) {
	 pick_info naii = atom_pick(event);
	 if (naii.success == GL_TRUE) {
	    mmdb::Atom *at = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index];
	    mmdb::Residue *res = at->residue;
	    std::string resname(res->name);
	    if (resname == "WAT" || resname == "HOH") { 
	       normal_cursor();

	       // Delete (any) hydrogens in the residue, then delete
	       // the atom.
	       // What if they click on a hydrogen in the residue?
	       // The oxygen gets left.
	       //
	       //   So, we should check the element and if it is a
	       //      hydrogen, we should delete the residue.
	       //   else
	       //      do what we do now
	       //        (i.e. delete residue hydrogens,
	       //         delete_atom())
	       // OK.

	       std::string chain_id(res->GetChainID());
	       int resno = res->GetSeqNum();
	       std::string altloc(at->altLoc);
	       std::string inscode(at->GetInsCode());
	       std::string atom_name(at->name);
	       std::string ele = at->element;
	       if (ele == " H") {
		  delete_residue_with_full_spec(naii.imol, naii.model_number, chain_id.c_str(),
						resno, inscode.c_str(), altloc.c_str());
	       } else { 
		  molecules[naii.imol].delete_residue_hydrogens(chain_id, resno, inscode, altloc);
		  delete_atom(naii.imol, chain_id.c_str(), resno, inscode.c_str(),
			      atom_name.c_str(), altloc.c_str());
		  delete_object_handle_delete_dialog(destroy_delete_dialog_flag_by_ctrl_press);
		  pick_pending_flag = 0;
		  run_post_manipulation_hook(naii.imol, DELETED);
		  item_deleted = true;
		  imol_delete = naii.imol;
	       }
	    }
	 } else {

	    if (show_symmetry) {
	       coot::Symm_Atom_Pick_Info_t symm_nearest_atom_index_info = symmetry_atom_pick(); 

	       if (symm_nearest_atom_index_info.success == GL_TRUE) {
		  int im = symm_nearest_atom_index_info.imol;
		  int index = symm_nearest_atom_index_info.atom_index;
		  mmdb::Residue *res = molecules[im].atom_sel.atom_selection[index]->residue;
		  std::string resname(res->name);
		  if (resname == "WAT" || resname == "HOH") {
		     // Note of course if we don't delete residue
		     // hydrogens as we do above, then the atom_index
		     // doesn't go out of date, so delete_atom_by_atom_index() is fine.
		     delete_atom_by_atom_index(im, index,
					       destroy_delete_dialog_flag_by_ctrl_press);
		     normal_cursor();
		     run_post_manipulation_hook(im, DELETED);
		     pick_pending_flag = 0;
		     item_deleted = true;
		     imol_delete = im;
		  }
	       }
	    }
	 }
      }

      // side chain
      if (g.delete_item_sidechain) {
	 pick_info naii = atom_pick(event);
	 if (naii.success == GL_TRUE) {
	    mmdb::Residue *res = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index]->residue;
	    std::string resname(res->name);
	    if (resname != "WAT" && resname != "HOH") {
	       normal_cursor();
	       int resno = res->GetSeqNum();
	       const char *ins_code = res->GetInsCode();
	       const char *chain_id = res->GetChainID();
	       delete_residue_sidechain(naii.imol, chain_id, resno, ins_code,
					destroy_delete_dialog_flag_by_ctrl_press);
	       g.update_environment_distances_maybe(naii.atom_index, naii.imol);
	       run_post_manipulation_hook(naii.imol, DELETED);
	       item_deleted = true;
	       imol_delete = naii.imol;
	    }
	 }
      }

      // side chain range
      if (g.delete_item_sidechain_range) {
	 pick_info naii = atom_pick(event);
	 if (naii.success == GL_TRUE) {
	    if (g.delete_item_sidechain_range == 1) {
	       // This was the first click:
	       molecules[naii.imol].add_to_labelled_atom_list(naii.atom_index);
	       mmdb::Atom *at = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index];
	       g.delete_item_sidechain_range_1 = coot::residue_spec_t(at->residue);
	       g.delete_item_sidechain_range_1_imol = naii.imol;
	       // so set up to pick another atom
	       g.delete_item_sidechain_range = 2;
	    } else {
	       if (g.delete_item_sidechain_range == 2) {
		  // This was the second click:
		  mmdb::Atom *at = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index];
		  coot::residue_spec_t res2(at->residue);
		  if (naii.imol == g.delete_item_sidechain_range_1_imol) {
		     if (res2.model_number == g.delete_item_sidechain_range_1.model_number) {
			g.delete_sidechain_range(naii.imol, g.delete_item_sidechain_range_1, res2);
			pick_pending_flag = 0;
			g.delete_item_sidechain_range = 0; //reset for next time, or 1?
			run_post_manipulation_hook(naii.imol, DELETED);
			item_deleted = true;
			imol_delete = naii.imol;
		     } else {
			pick_pending_flag = 0;
			normal_cursor();
			std::string s = "Picked atoms not in same model.";
			add_status_bar_text(s);
		     }
		  } else {
		     pick_pending_flag = 0;
		     normal_cursor();
		     std::string s = "Picked atoms not in same molecule.";
		     add_status_bar_text(s);
		  }
	       }
	    }
	    graphics_draw();
	 }
      }

      // chain
      if (g.delete_item_chain) {
	 pick_info naii = atom_pick(event);
	 if (naii.success == GL_TRUE) {
	    normal_cursor();
	    mmdb::Atom *at = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index];
	    std::string chain_id = at->residue->chain->GetChainID();
	    delete_chain(naii.imol, chain_id.c_str()); // handles dialog
	    graphics_draw();
	    run_post_manipulation_hook(naii.imol, DELETED);
	    item_deleted = true;
	    imol_delete = naii.imol;
	 }
      }

      // residue
      if (g.delete_item_residue) {
	 pick_info naii = atom_pick(event);
	 if (naii.success == GL_TRUE) {
	    mmdb::Residue *res = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index]->residue;
	    std::string resname(res->name);
	    if (resname != "WAT" && resname != "HOH") {
	       normal_cursor();
	       delete_residue_by_atom_index(naii.imol, naii.atom_index,
					    destroy_delete_dialog_flag_by_ctrl_press);
	       g.update_environment_distances_maybe(naii.atom_index, naii.imol);
	       run_post_manipulation_hook(naii.imol, DELETED);
	       pick_pending_flag = 0;
	       item_deleted = true;
	       imol_delete = naii.imol;
	    }
	 } else { 

	    if (show_symmetry) { 
	       coot::Symm_Atom_Pick_Info_t symm_nearest_atom_index_info = symmetry_atom_pick();

	       if (1) {
		  if (symm_nearest_atom_index_info.success == GL_TRUE) {
		     int im = symm_nearest_atom_index_info.imol;
		     int index = symm_nearest_atom_index_info.atom_index;
		     mmdb::Residue *res = molecules[im].atom_sel.atom_selection[index]->residue;
		     std::string resname(res->name);
		     if (resname != "WAT" && resname != "HOH") {
			pick_pending_flag = 0;
			normal_cursor();
			int im = symm_nearest_atom_index_info.imol;
			delete_residue_by_atom_index(im, index,
						     destroy_delete_dialog_flag_by_ctrl_press);
			run_post_manipulation_hook(im, DELETED);
			item_deleted = true;
			imol_delete = im;
		     }
		  }
	       } else { // not used
	       
		  if (symm_nearest_atom_index_info.success == GL_TRUE) {
		  
		     std::string s = "That was a symmetry atom\n";
		     s += "Coot currently doesn't delete symmetry items";
		     GtkWidget *w = wrapped_nothing_bad_dialog(s);
		     gtk_widget_show(w);
		  }
	       } 
	    }
	 }
      }

      // residue's hydrogens
      if (g.delete_item_residue_hydrogens) {
	 pick_info naii = atom_pick(event);
	 if (naii.success == GL_TRUE) {
	    delete_residue_hydrogens_by_atom_index(naii.imol, naii.atom_index,
						   destroy_delete_dialog_flag_by_ctrl_press);
	    normal_cursor();
	    pick_pending_flag = 0;
	 } else { 
	    // Let's face it, this is pretty unlikely ever to happen....

	    if (show_symmetry) { 
	       coot::Symm_Atom_Pick_Info_t symm_nearest_atom_index_info = symmetry_atom_pick();

	       if (1) {
		  if (symm_nearest_atom_index_info.success == GL_TRUE) {
		     int index = symm_nearest_atom_index_info.atom_index;
		     int im = symm_nearest_atom_index_info.imol;
		     delete_residue_hydrogens_by_atom_index(im, index,
							    destroy_delete_dialog_flag_by_ctrl_press);
		     run_post_manipulation_hook(im, DELETED);
		     pick_pending_flag = 0;
		     normal_cursor();
		     item_deleted = true;
		     imol_delete = im;
		  }

	       } else { // old unused code

		  if (symm_nearest_atom_index_info.success == GL_TRUE) {

		     // int im = symm_nearest_atom_index_info.imol;
		     std::string s = "That was a symmetry atom\n";
		     s += "Coot currently doesn't delete symmetry atoms";
		     GtkWidget *w = wrapped_nothing_bad_dialog(s);
		     gtk_widget_show(w);
		  }
	       }
	    }
	 }
      }

      // residue zone
      if (g.delete_item_residue_zone) {
	 pick_info naii = atom_pick(event);
	 if (naii.success == GL_TRUE) {
	    if (g.delete_item_residue_zone == 1) {
	       // This was the first click:
	       molecules[naii.imol].add_to_labelled_atom_list(naii.atom_index);

	       // save a residue spec and the molecule number for the first residue
	       mmdb::Atom *at = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index];
	       g.delete_item_residue_zone_1 = coot::residue_spec_t(at->residue);
	       g.delete_item_residue_zone_1_imol = naii.imol;
	       // so set up to pick another atom
	       g.delete_item_residue_zone = 2;
	    } else {

	       // This was the second atom picked
	       mmdb::Atom *at = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index];
	       coot::residue_spec_t res2(at->residue);
	       // This is a graphics_info_t function, unlike the other
	       // delete functions.  This is because we use
	       // coot::residue_spec_t which can't be introduced int
	       // c-interface.h
	       if (naii.imol == g.delete_item_residue_zone_1_imol) {
		  if (res2.model_number == g.delete_item_residue_zone_1.model_number) {

		     delete_residue_range(naii.imol, g.delete_item_residue_zone_1, res2);
		     pick_pending_flag = 0;
		     g.delete_item_residue_zone = 1; //reset for next time
		     run_post_manipulation_hook(naii.imol, DELETED);
		     item_deleted = true;
		     imol_delete = naii.imol;
		  } else {
		     pick_pending_flag = 0;
		     normal_cursor();
		     std::string s = "Picked atoms not in same model.";
		     add_status_bar_text(s);
		  } 
	       } else {
		  pick_pending_flag = 0;
		  normal_cursor();
		  std::string s = "Picked atoms not in same molecule.";
		  add_status_bar_text(s);
	       }
	    }
	    graphics_draw();
	 }
      }
   }
}

int
graphics_info_t::check_if_in_rigid_body_define(GdkEventButton *event) {

   int state = 0;
   graphics_info_t g;
   if (g.in_rigid_body_define) {
      
      pick_info naii = atom_pick(event);
      if ( naii.success == GL_TRUE ) {
	 molecules[naii.imol].add_to_labelled_atom_list(naii.atom_index);

	 if (g.in_rigid_body_define == 1) {

	    g.residue_range_atom_index_1 = naii.atom_index;
	    g.in_rigid_body_define = 2;
	    g.imol_rigid_body_refine = naii.imol;

	    if (a_is_pressed) { 
	       short int auto_range_flag = 1;
	       execute_rigid_body_refine(auto_range_flag);
	       in_range_define = 0;
	       pick_pending_flag = 0;
	       normal_cursor();
	       model_fit_refine_unactive_togglebutton("model_refine_dialog_rigid_body_togglebutton");
	    }

	 } else {

	    if (g.in_rigid_body_define == 2) {
	       if (naii.imol == g.imol_rigid_body_refine) {
		  g.residue_range_atom_index_2 = naii.atom_index;
		  
		  // now do it
		  execute_rigid_body_refine(0);

	       } else {
		  cout << "Rigid Body: That atom was not in the "; 
		  cout << "same molecule as the previous atom" << endl; 
		  cout << "Cancelling selection" << endl;
	       }
	    }
	    g.in_rigid_body_define = 0;
	    pick_pending_flag = 0;
	    normal_cursor();
	    model_fit_refine_unactive_togglebutton("model_refine_dialog_rigid_body_togglebutton");
	 }
	 graphics_draw(); // let's see the label
      }
   }
   return state; // ignored currently.
}

void
graphics_info_t::check_if_in_terminal_residue_define(GdkEventButton *event) { 
   graphics_info_t g;
   if (g.in_terminal_residue_define) { 
      pick_info naii = atom_pick(event);
      if (naii.success == GL_TRUE) { 
	 
	 // need molecule number, chain_id, terminus type and
	 // residue type (string).

	 std::string term_type = g.molecules[naii.imol].get_term_type(naii.atom_index);
	 mmdb::Residue *res_p = g.molecules[naii.imol].atom_sel.atom_selection[naii.atom_index]->GetResidue();
	 std::string chain_id(g.molecules[naii.imol].atom_sel.atom_selection[naii.atom_index]->GetChainID());
	 mmdb::Atom *at = g.molecules[naii.imol].atom_sel.atom_selection[naii.atom_index];

	 // we check on term_type inside the function
	 //
	 watch_cursor();
	 short int add_it_now_flag = g.add_terminal_residue_immediate_addition_flag;
	 if (add_terminal_residue_do_post_refine) {
	    add_it_now_flag = 1;
	 }
	 if (!coot::util::is_nucleotide_by_dict_dynamic_add(res_p, Geom_p())) {
	    g.execute_add_terminal_residue(naii.imol,
					   term_type,
					   res_p,
					   chain_id,
					   g.add_terminal_residue_type,  // eg. "ALA" or "UNK"
					   add_it_now_flag);
	 } else {
	    g.execute_simple_nucleotide_addition(naii.imol, term_type, res_p, chain_id);
	 }
	 g.in_terminal_residue_define = 0;
	 pick_pending_flag = 0;
	 
	 if (add_terminal_residue_do_post_refine) {
	    // Run refine zone with autoaccept, autorange on
	    // the "clicked" atom:
	    refine_auto_range(naii.imol, chain_id.c_str(), at->GetSeqNum(),
			      at->altLoc);
	 }

	 normal_cursor();
	 model_fit_refine_unactive_togglebutton("model_refine_dialog_fit_terminal_residue_togglebutton");
      } 
   } 
} 

int
graphics_info_t::check_if_in_rot_trans_define(GdkEventButton *event) { 
   
   int state = 0;
   graphics_info_t g;
   if (g.in_rot_trans_object_define) {
      if (g.rot_trans_object_type == ROT_TRANS_TYPE_ZONE) { 
	 pick_info naii = atom_pick(event); 
	 if (naii.success == GL_TRUE) { 
	    molecules[naii.imol].add_to_labelled_atom_list(naii.atom_index);
	    if (g.in_rot_trans_object_define == 1) { 
	       g.rot_trans_atom_index_1 = naii.atom_index;
	       g.in_rot_trans_object_define = 2;
	       g.imol_rot_trans_object = naii.imol;
	    } else { 
	       
	       if (g.in_rot_trans_object_define == 2) { 
		  if (naii.imol == g.imol_rot_trans_object) { 
		     g.rot_trans_atom_index_2 = naii.atom_index;
		     
		     // now we are setup to move:
		     g.execute_rotate_translate_ready();
		     fleur_cursor();
		  } else { 
		     std::cout << "Rotation/Translations: that atom was ";
		     std::cout << "not in the same molecule as the "
			       << "previous atom" << std::endl;
		     std::cout << "Cancelling selection" << std::endl;
		  }
	       }
	       g.in_rot_trans_object_define = 0;
	       pick_pending_flag = 0;
	       normal_cursor();
	       model_fit_refine_unactive_togglebutton("model_refine_dialog_rot_trans_togglebutton");
	    }
	    graphics_draw(); // let's see the label
	 }
      }

      if (g.rot_trans_object_type == ROT_TRANS_TYPE_CHAIN) {
	 // One click
	 pick_info naii = atom_pick(event); 
	 if (naii.success == GL_TRUE) { 
	    molecules[naii.imol].add_to_labelled_atom_list(naii.atom_index);
	    g.imol_rot_trans_object = naii.imol;
	    g.rot_trans_atom_index_1 = naii.atom_index;
	    g.execute_rotate_translate_ready();
	    fleur_cursor();
	    g.in_rot_trans_object_define = 0;
	    model_fit_refine_unactive_togglebutton("model_refine_dialog_rot_trans_togglebutton");
	 }
      }

      
      if (g.rot_trans_object_type == ROT_TRANS_TYPE_MOLECULE) {
	 // One click
	 pick_info naii = atom_pick(event); 
	 if (naii.success == GL_TRUE) { 
	    molecules[naii.imol].add_to_labelled_atom_list(naii.atom_index);
	    g.imol_rot_trans_object = naii.imol;
	    g.rot_trans_atom_index_1 = naii.atom_index;
	    g.execute_rotate_translate_ready();
	    fleur_cursor();
	    g.in_rot_trans_object_define = 0;
	    model_fit_refine_unactive_togglebutton("model_refine_dialog_rot_trans_togglebutton");
	 }
      }
      
   }
   return state; // ignored, currently
} 


void
graphics_info_t::check_if_in_db_main_define(GdkEventButton *event) {

   graphics_info_t g;

   // 20180721 change this so that it needs only a single click.

   if (g.in_db_main_define) { 
      pick_info naii = atom_pick(event);
      if (naii.success == GL_TRUE) { 
	 molecules[naii.imol].add_to_labelled_atom_list(naii.atom_index);
	 if (g.in_db_main_define == 1) {
	    g.db_main_atom_index_1 = naii.atom_index;
	    g.in_db_main_define = 2;
	    g.db_main_imol = naii.imol;
	    g.execute_db_main();
	    g.in_db_main_define = 0;
	    pick_pending_flag = 0;
	    normal_cursor();
	    other_modelling_tools_unactive_togglebutton("model_refine_dialog_db_main_togglebutton");
	 }
	 graphics_draw(); // let's see the label
      }
   }
}

void
graphics_info_t::check_if_in_rotamer_define(GdkEventButton *event) {

   graphics_info_t g;
   if (g.in_rotamer_define) {
      pick_info naii = atom_pick(event);
      if (naii.success == GL_TRUE) {
	 g.do_rotamers(naii.atom_index, naii.imol);
	 g.in_rotamer_define = 0;
	 pick_pending_flag = 0;
	 normal_cursor();
	 model_fit_refine_unactive_togglebutton("model_refine_dialog_rotamer_togglebutton");
	 if (0) { 
	    if (moving_atoms_asc) {
	       std::cout << "debug moving atoms to moving-atoms.pdb" << std::endl;
	       moving_atoms_asc->mol->WritePDBASCII("moving-atoms.pdb");
	    } else {
	       std::cout << "debug no moving atoms object" << std::endl;
	    }
	 }
      }
   }
}

void
graphics_info_t::check_if_in_mutate_define(GdkEventButton *event) {

   graphics_info_t g;
   if (g.in_mutate_define) {
      pick_info naii = atom_pick(event);
      if (naii.success == GL_TRUE) {
	 g.mutate_residue_imol = naii.imol;
	 g.mutate_residue_atom_index = naii.atom_index;
	 mmdb::Residue *r = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index]->residue;

	 // is it sensible to do it two ways like this?  I'm not sure.
	 // Perhaps there should be one function that wraps both
	 // methods?
	 bool is_nuc = 0;
	 is_nuc = coot::util::is_nucleotide_by_dict_dynamic_add(r, Geom_p());
	 if (! is_nuc)
	    is_nuc = coot::util::is_nucleotide(r);
	 
	 if (is_nuc) {
	    GtkWidget *w = create_nucleic_acid_base_chooser_dialog();
	    gtk_widget_show(w);
	 } else { 
	    GtkWidget *widget = wrapped_create_residue_type_chooser_window(1);
	    gtk_widget_show(widget);
	 }
	 g.in_mutate_define = 0;
	 g.residue_type_chooser_auto_fit_flag = 0;
	 pick_pending_flag = 0;
	 normal_cursor();
	 model_fit_refine_unactive_togglebutton("model_refine_dialog_mutate_togglebutton");
      }
   }
} 

void
graphics_info_t::check_if_in_mutate_auto_fit_define(GdkEventButton *event) {

   graphics_info_t g;
   if (g.in_mutate_auto_fit_define) {
      pick_info naii = atom_pick(event);
      if (naii.success == GL_TRUE) {
	 g.mutate_auto_fit_residue_imol = naii.imol;
	 g.mutate_auto_fit_residue_atom_index = naii.atom_index;
	 GtkWidget *widget = wrapped_create_residue_type_chooser_window(0);
	 gtk_widget_show(widget);
	 g.in_mutate_auto_fit_define = 0;
	 g.residue_type_chooser_auto_fit_flag = 1;
	 pick_pending_flag = 0;
	 normal_cursor();
	 model_fit_refine_unactive_togglebutton("model_refine_dialog_mutate_auto_fit_togglebutton");
      }
   } 
}

GtkWidget *
graphics_info_t::wrapped_create_residue_type_chooser_window(bool show_stub_option_flag) const {

   GtkWidget *w = create_residue_type_chooser_window();
   GtkWidget *b = lookup_widget(w, "residue_type_chooser_stub_checkbutton");

   if (show_stub_option_flag == 0) 
      gtk_widget_hide(b);

   return w;
} 

void
graphics_info_t::check_if_in_auto_fit_define(GdkEventButton *event) { 

   if (in_auto_fit_define) { 

      pick_info naii = atom_pick(event);
      if (naii.success == GL_TRUE) { 
	 
	 int imol = naii.imol;
	 int imol_map = Imol_Refinement_Map();  // -1 is a magic number

	 if (imol_map == -1 ) { // magic number
	    //
	    std::cout << "Please set a map against which the refinement should occur"
		      << std::endl;
	    show_select_map_dialog();
	 } else {

	    mmdb::Atom *atom_p = molecules[imol].atom_sel.atom_selection[naii.atom_index];
	    int resno = atom_p->GetSeqNum();
	    std::string inscode(atom_p->GetInsCode());
	    std::string chain(atom_p->GetChainID());
	    std::string altloc(atom_p->altLoc);
	    
	    float f = molecules[imol].auto_fit_best_rotamer(rotamer_search_mode,
							    resno, altloc, inscode, chain,
							    imol_map, rotamer_fit_clash_flag,
							    rotamer_lowest_probability, *Geom_p());
	    if (rotamer_auto_fit_do_post_refine_flag) {
	      // Run refine zone with autoaccept, autorange on the "clicked" atom:
	      refine_auto_range(naii.imol, chain.c_str(), resno, altloc.c_str());
	    }
	    if (graphics_info_t::reset_b_factor_moved_atoms_flag) {
	      reset_b_factor_residue_range(naii.imol, chain.c_str(), resno, resno);
	    }
	    update_geometry_graphs(&atom_p->residue, 1, imol, imol_map);
	    std::cout << "Fitting score for best rotamer: " << f << std::endl;
	    graphics_draw();
	 }
	 in_auto_fit_define = 0;
	 pick_pending_flag = 0;
	 normal_cursor();
	 model_fit_refine_unactive_togglebutton("model_refine_dialog_auto_fit_rotamer_togglebutton");
      }
   } 
} 

void
graphics_info_t::check_if_in_add_alt_conf_define(GdkEventButton *event) {
   
   if (in_add_alt_conf_define) { 

      pick_info naii = atom_pick(event);
      
      if (naii.success == GL_TRUE) { 

	 int imol = naii.imol;

	 // bla bla... more stuff here. [20080618 Eh? what do I mean?]

	 if (alt_conf_split_type < 2) { 
	    // not a range, do it now on a residue

	    split_residue(imol, naii.atom_index);
	    in_add_alt_conf_define = 0;
	    if (add_alt_conf_dialog)
	       gtk_widget_destroy(add_alt_conf_dialog);
	    pick_pending_flag = 0;
	    normal_cursor();
	 } else { 

	    // was a range

	    // first point?
	    if (in_add_alt_conf_define == 1) { 
	       add_alt_conf_atom_index = naii.atom_index;
	       add_alt_conf_imol = imol;
	       // setup for next point:
	       std::cout << "Now pick another atom\n";
	       in_add_alt_conf_define = 2;
	    } else {
	       // this was the second point
	       if (imol == add_alt_conf_imol) { 
		  split_residue_range(imol, add_alt_conf_atom_index, naii.atom_index);
	       } else { 
		  std::cout << "Sorry those 2 atoms were not in the same molecule" 
			    << std::endl;
	       }
	       in_add_alt_conf_define = 0;
	       pick_pending_flag = 0;
	       normal_cursor();
	       if (add_alt_conf_dialog)
		  gtk_widget_destroy(add_alt_conf_dialog);
	    }
	 }
      } 
   }
}



   // I think we want a button, not a togglebutton so this does not apply:
	 // 	 model_fit_refine_unactive_togglebutton("model_refine_dialog_add_alt_conf_togglebutton");

void 
graphics_info_t::check_if_in_edit_phi_psi_define(GdkEventButton *event) {

   if (in_edit_phi_psi_define) { 

      pick_info naii = atom_pick(event);
      
      if (naii.success == GL_TRUE) { 

	 // generate phi/psi moving atoms 
	 
	 edit_phi_psi_atom_index = naii.atom_index;
	 execute_edit_phi_psi(naii.atom_index, naii.imol);

	 in_edit_phi_psi_define = 0;
	 pick_pending_flag = 0;
	 normal_cursor();
 	 model_fit_refine_unactive_togglebutton("model_refine_dialog_edit_phi_psi_togglebutton");
      }
   }
}

void
graphics_info_t::check_if_in_edit_chi_angles_define(GdkEventButton *event) {

   if (in_edit_chi_angles_define) {

      pick_info naii = atom_pick(event);
      if (naii.success == GL_TRUE) {

	 execute_edit_chi_angles(naii.atom_index, naii.imol); // pop up the dialog (and set
	                                                      // atom spec of clicked atom)
	 in_edit_chi_angles_define = 0;
	 pick_pending_flag = 0;
	 normal_cursor();
 	 model_fit_refine_unactive_togglebutton("model_refine_dialog_edit_chi_angles_togglebutton");
	 moving_atoms_move_chis_flag = 1;
      }
   }
}


void
graphics_info_t::check_if_in_180_degree_flip_define(GdkEventButton *event) {

   if (in_180_degree_flip_define) {
      pick_info naii = atom_pick(event);
      if (naii.success == GL_TRUE) {
	 mmdb::Atom *at = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index];
	 int resno = at->GetSeqNum();
	 mmdb::Residue *residue = at->residue;
	 std::string chain_id = at->GetChainID();
	 std::string inscode  = at->GetInsCode();
	 std::string alt_conf = at->altLoc;
	 std::string resname  = at->GetResName();
	 int istatus =
	    molecules[naii.imol].do_180_degree_side_chain_flip(chain_id, resno,
								 inscode, alt_conf, Geom_p());
	 std::string s;
	 if (istatus) { 
	    s = "Chi angle on residue ";
	    s += chain_id;
	    s += graphics_info_t::int_to_string(resno);
	    s += " ";
	    s += resname;
	    s += " successfully flipped.";
	    // update graphs here
	    // make a molecule from the residue so that we can call
	    // update_geometry_graphs(*moving_atoms_asc, imol_moving_atoms);
	    std::pair<mmdb::Manager *, int> mp = 
	       coot::util::create_mmdbmanager_from_res_selection(molecules[naii.imol].atom_sel.mol,
								 &residue, 1, 0, 0,
								 alt_conf, chain_id, 0);
	    atom_selection_container_t asc = make_asc(mp.first);
	    asc.UDDOldAtomIndexHandle = mp.second;
	    update_geometry_graphs(asc, naii.imol);
	    graphics_draw();
	 } else {
	    s = "Problem flipping chi angle on residue ";
	    s += chain_id;
	    s += graphics_info_t::int_to_string(resno);
	    s += " ";
	    s += resname;
	    s += ". Not done.";
	 }
	 add_status_bar_text(s);
	 
	 in_180_degree_flip_define = 0;
	 pick_pending_flag = 0;
	 normal_cursor();
 	 model_fit_refine_unactive_togglebutton("model_refine_dialog_do_180_degree_sidechain_flip_togglebutton");
      }
   }
}

void
graphics_info_t::check_if_in_torsion_general_define(GdkEventButton *event) {

   if (in_torsion_general_define) {
      pick_info nearest_atom_index_info = atom_pick(event);
      if (nearest_atom_index_info.success == GL_TRUE) {
	 std::cout << "in_torsion_general_define was " << in_torsion_general_define
		   << std::endl;
	 int im = nearest_atom_index_info.imol;
	 molecules[im].add_to_labelled_atom_list(nearest_atom_index_info.atom_index);
	 if (in_torsion_general_define == 1) {
	    std::cout << " 1" << std::endl;
	    torsion_general_atom_index_1 = nearest_atom_index_info.atom_index;
	    torsion_general_atom_index_1_mol_no = im;
	    in_torsion_general_define = 2;
	 } else { 
	    if (in_torsion_general_define == 2) {
	       std::cout << " 2" << std::endl;
	       torsion_general_atom_index_2 = nearest_atom_index_info.atom_index;
	       torsion_general_atom_index_2_mol_no = im;
	       in_torsion_general_define = 3;
	    } else { 
	       if (in_torsion_general_define == 3) {
		  std::cout << " 3" << std::endl;
		  torsion_general_atom_index_3 = nearest_atom_index_info.atom_index;
		  torsion_general_atom_index_3_mol_no = im;
		  in_torsion_general_define = 4;
	       } else { 
		  if (in_torsion_general_define == 4) {
		     std::cout << " 4" << std::endl;
		     torsion_general_atom_index_4 = nearest_atom_index_info.atom_index;
		     torsion_general_atom_index_4_mol_no = im;
		     // act on the torsion general setup
		     execute_torsion_general();
		     in_torsion_general_define = 0;
		     model_fit_refine_unactive_togglebutton("model_refine_dialog_torsion_general_togglebutton");
		     normal_cursor();
		  }
	       }
	    }
	 }
	 graphics_draw();
      }
   }
}


void 
graphics_info_t::check_if_in_edit_backbone_torsion_define(GdkEventButton *event) {

   if (in_backbone_torsion_define) { 
      
      pick_info naii = atom_pick(event);
      if (naii.success == GL_TRUE) { 
	 
	 execute_setup_backbone_torsion_edit(naii.imol, naii.atom_index);
	 in_backbone_torsion_define = 0;
	 pick_pending_flag = 0;
	 normal_cursor();
	 model_fit_refine_unactive_togglebutton("model_refine_dialog_edit_backbone_torsions_togglebutton");
      }
   } 
}

void
graphics_info_t::check_if_in_residue_partial_alt_locs(GdkEventButton *event) {

   // like edit chi angles 
   if (in_residue_partial_alt_locs_define) {

      pick_info naii = atom_pick(event);
      if (naii.success == GL_TRUE) {

	 int imol = naii.imol;
	 in_residue_partial_alt_locs_define = 0;
	 pick_pending_flag = 0;
	 normal_cursor();
 	 model_fit_refine_unactive_togglebutton("model_refine_dialog_residue_partial_alt_locs_togglebutton");
	 moving_atoms_move_chis_flag = 1;

	 imol_residue_partial_alt_locs = imol; // used in setup_flash_bond()
	 mmdb::Residue *residue_p = molecules[imol].atom_sel.atom_selection[naii.atom_index]->GetResidue();
	 residue_partial_alt_locs_spec = coot::residue_spec_t(residue_p);
	 std::string res_type(residue_p->GetResName());
	 //
	 edit_chi_edit_type mode = RESIDUE_PARTIAL_ALT_LOCS;
	 int ires = wrapped_create_edit_chi_angles_dialog(res_type, mode);
	 
      }
   }
}



void
graphics_info_t::check_if_in_save_symmetry_define(GdkEventButton *event) {

   if (in_save_symmetry_define) {

      coot::Symm_Atom_Pick_Info_t naii = symmetry_atom_pick();

      if (naii.success == GL_TRUE) { 

	 in_save_symmetry_define = 0;
	 GtkWidget *w = coot_save_symmetry_chooser();

	 // attach symmetry info to the widget for use in the OK
	 // button callback:
	 coot::Symm_Atom_Pick_Info_t *save_pick_info = new coot::Symm_Atom_Pick_Info_t;
	 *save_pick_info = naii;

	 g_object_set_data(G_OBJECT(w), "save_pick_info", save_pick_info);

// 	 std::string filename = "molecule-";
// 	 filename += int_to_string(naii.imol);
	 int imol = naii.imol;
	 std::string filename = molecules[imol].name_sans_extension(0);
	 filename += "_symmetry_";
	 // 20111117 - previously we had just the symmetry number.
	 // Now let's add the translation too.
	 std::string fill_char = "0";
	 filename += int_to_string(naii.symm_trans.isym()+1); // +1 because of
             	                                              // zero indexing -> real-world/dictionary indexing.
	 // x
	 if (naii.symm_trans.x() < 0)
	    fill_char = "";
	 else
	    fill_char = "0";
	 filename += fill_char;
	 filename += int_to_string(naii.symm_trans.x());
	 // y 
	 if (naii.symm_trans.y() < 0)
	    fill_char = "";
	 else
	    fill_char = "0";
	 filename += fill_char;
	 filename += int_to_string(naii.symm_trans.y());
	 // z
	 if (naii.symm_trans.z() < 0)
	    fill_char = "";
	 else
	    fill_char = "0";
	 filename += fill_char;
	 filename += int_to_string(naii.symm_trans.z());
	 
	 // 
	 filename += ".pdb";

	 gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(w), filename.c_str());

	 add_ccp4i_project_optionmenu(w, COOT_COORDS_FILE_SELECTION);
         
         add_filename_filter_button(w, COOT_COORDS_FILE_SELECTION);

	 normal_cursor();
	 gtk_widget_show(w);
	 pick_pending_flag = 0;
      }
   } 

}


void
graphics_info_t::check_if_in_cis_trans_convertion_define(GdkEventButton *event) {

   if (in_cis_trans_convert_define) {
      pick_info naii = atom_pick(event);
      if (naii.success == GL_TRUE) { 
	 in_cis_trans_convert_define = 0;
	 pick_pending_flag = 0;
	 short int is_N_flag = 0;
	 mmdb::Atom *at = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index];
	 std::string atom_name = at->name;
	 if (atom_name == " N  ")
	    is_N_flag = 1;
	 cis_trans_conversion(at, naii.imol, is_N_flag);
	 normal_cursor();
	 other_modelling_tools_unactive_togglebutton("cis_trans_conversion_toggle_button");
      }
   }
}


void
graphics_info_t::check_if_in_reverse_direction_define(GdkEventButton *event) {

   if (in_reverse_direction_define) {
      pick_info naii = atom_pick(event);
      if (naii.success == GL_TRUE) {
	 in_reverse_direction_define = 0;
	 pick_pending_flag = 0;
	 int imol = naii.imol;
	 mmdb::Atom *at = molecules[imol].atom_sel.atom_selection[naii.atom_index];
	 const char *chain_id = at->GetChainID();
	 int resno = at->GetSeqNum();
	 reverse_direction_of_fragment(imol, chain_id, resno);
	 other_modelling_tools_unactive_togglebutton("reverse_fragment_direction_togglebutton"); 
      }
   }
}

void
graphics_info_t::check_if_in_lsq_plane_define(GdkEventButton *event) {

   if (in_lsq_plane_define) {
      pick_info naii = atom_pick(event);
      if (naii.success == GL_TRUE) {
	 add_lsq_plane_atom(naii.imol, naii.atom_index);
      }
   }
}



void
graphics_info_t::check_if_in_lsq_plane_deviant_atom_define(GdkEventButton *event) {

   if (in_lsq_plane_deviation) {
      pick_info naii = atom_pick(event);
      if (naii.success == GL_TRUE) {
	 measure_lsq_plane_deviant_atom(naii.imol, naii.atom_index);
      }
   }
}


void
graphics_info_t::check_if_in_fixed_atom_define(GdkEventButton *event,
					       const GdkModifierType &state) {

   
   if (in_fixed_atom_define != coot::FIXED_ATOM_NO_PICK) {
      // we were listening for a pick then...
      bool pick_state = 0;
      if (in_fixed_atom_define == coot::FIXED_ATOM_FIX) {
	 pick_state = 1;
      }
      if (in_fixed_atom_define == coot::FIXED_ATOM_UNFIX) {
	 pick_state = 0;
      }

      // pick and fix are interchanged here. I mean on
      // UNFIX/pick_state=0 that an atom that is marked as FIXED (it
      // should have a dot) becomes unmarked as fixed.
      // 
      pick_info naii = atom_pick(event);
      if (naii.success == GL_TRUE) {
	 coot::atom_spec_t as(molecules[naii.imol].atom_sel.atom_selection[naii.atom_index]);
	 mark_atom_as_fixed(naii.imol, as, pick_state);
	 std::cout << "   " << as << " is a marked as fixed " << pick_state << std::endl;
	 graphics_draw();

	 // Sadly, Ctrl + left mouse click is intercepted upstream of
	 // this and we don't get to see it here.  Currently (20080212).
	 
	 if (! (state & GDK_CONTROL_MASK)) { 
	    // Ctrl key is not pressed.
	    if (!fixed_atom_dialog) {
	       std::cout << "Ooops fixed atom dialog has gone!" << std::endl;
	    } else { 
	       GtkWidget *button1 = lookup_widget(fixed_atom_dialog,   "fix_atom_togglebutton");
	       GtkWidget *button2 = lookup_widget(fixed_atom_dialog, "unfix_atom_togglebutton");
	       if (button1)
		  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button1), FALSE);
	       if (button2)
		  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button2), FALSE);
	       in_fixed_atom_define = coot::FIXED_ATOM_NO_PICK;
	       normal_cursor();
	    }
	 }
      }
   }
}


void
graphics_info_t::check_if_in_base_pairing_define(GdkEventButton *event) {

   if (in_base_paring_define) {
      pick_info naii = atom_pick(event);
      if (naii.success == TRUE) {
	 in_base_paring_define = 0;
	 mmdb::Atom *at = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index];
	 int res_no = at->GetSeqNum();
	 std::string chain_id = at->GetChainID();
	 watson_crick_pair(naii.imol, chain_id.c_str(), res_no);
	 if (other_modelling_tools_dialog) {
	    GtkWidget *w = lookup_widget(other_modelling_tools_dialog,
					 "other_tools_base_pair_toggle_button");
	    if (w) {
	       gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), FALSE);
	    } 
	 } 
      } 
   } 
} 

void
graphics_info_t::check_if_in_multi_residue_torsion_define(GdkEventButton *event) {

   // Note: when we set in_mutate_auto_fit_define, we need to clear
   // the picked atom spec list
   
   if (in_multi_residue_torsion_define) {
      pick_info naii = atom_pick(event);
      if (naii.success == TRUE) {
	 int im = naii.imol;
	 mmdb::Atom *at = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index];
	 coot::residue_spec_t residue_spec(at);
	 if (std::find(multi_residue_torsion_picked_residue_specs.begin(),
		       multi_residue_torsion_picked_residue_specs.end(),
		       residue_spec) ==
	     multi_residue_torsion_picked_residue_specs.end()) 
	    multi_residue_torsion_picked_residue_specs.push_back(residue_spec);
	 
	 multi_residue_torsion_picked_residues_imol = naii.imol;
	 molecules[im].add_to_labelled_atom_list(naii.atom_index);
	 graphics_draw();
      }
   }
}

