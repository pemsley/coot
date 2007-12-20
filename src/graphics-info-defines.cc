/* src/graphics-info-defines.cc
 * 
 * Copyright 2004, 2005 by The University of York
 * Copyright 2007 by The University of York
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
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#if defined _MSC_VER
#include <windows.h>
#endif

#include "graphics-info.h"
#include "c-interface.h"
#include "interface.h"

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
   pick_pending_flag           = 0;

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
   names.push_back("model_refine_dialog_find_waters_button");
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
   names.push_back("model_refine_dialog_baton_button"); 
   names.push_back("model_refine_dialog_find_ligands_button");
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

   return iv;
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
	 std::cout << "info: clicked on imol: " << im << std::endl;
	 // a c-interface-info function...
	 output_residue_info_dialog(nearest_atom_index_info.atom_index, im); 
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
	    CAtom *atom1 = molecules[im].atom_sel.atom_selection[geometry_atom_index_1];
	    distance_pos_1 = coot::Cartesian(atom1->x, atom1->y, atom1->z);
	    std::cout << "click on a second atom" << std::endl;
	    graphics_draw();
	 } else {

	    // in_distance_define == 2
	    geometry_atom_index_2 =
	       nearest_atom_index_info.atom_index;
	    geometry_atom_index_2_mol_no =
	       nearest_atom_index_info.imol;

	    CAtom *atom1 = molecules[im].atom_sel.atom_selection[geometry_atom_index_2];
	    coot::Cartesian pos2 = coot::Cartesian(atom1->x, atom1->y, atom1->z);

	    display_geometry_distance_symm(distance_pos_1, pos2);

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
	       in_distance_define = 2;
	       std::cout << "click on a second atom" << std::endl;
	       graphics_draw();

	    } else { 

	       // in_distance_define == 2
	       coot::Cartesian pos2 = symm_nearest_atom_index_info.hybrid_atom.pos;
	       display_geometry_distance_symm(distance_pos_1, pos2);
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
      CAtom *atom = 0;
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
      CAtom *atom = 0;
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

//    std::cout << "DEBUG:: delete_item atom " << delete_item_atom
// 	     << " residue " << delete_item_residue
// 	     << " water " << delete_item_water
// 	     << " range " << delete_item_residue_zone
// 	     << " hydrogens " << delete_item_residue_hydrogens
// 	     << std::endl;

   graphics_info_t g;
   short int destroy_delete_dialog_flag_by_ctrl_press = 1;
   if (state & GDK_CONTROL_MASK)
      destroy_delete_dialog_flag_by_ctrl_press = 0;
   
   if (g.delete_item_widget) { 
      // atom
      if (g.delete_item_atom) {
	 pick_info naii = atom_pick(event);
	 if (naii.success == GL_TRUE) {
	    CResidue *res = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index]->residue;
	    std::string resname(res->name);
	    if (resname == "WAT" || resname == "HOH") {
	       std::string w = "Protection measure: ";
	       w += "Waters can only be deleted using Delete Water";
	       statusbar_text(w);
	    } else {
	       normal_cursor();
	       delete_atom_by_atom_index(naii.imol, naii.atom_index,
					 destroy_delete_dialog_flag_by_ctrl_press);
	       pick_pending_flag = 0;
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
		     pick_pending_flag = 0;
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
	    CAtom *at = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index];
	    CResidue *res = at->residue;
	    std::string resname(res->name);
	    if (resname == "WAT" || resname == "HOH") { 
	       normal_cursor();

	       // Delete (any) hydrogens in the residue, then delete
	       // the atom.
	       // What if they click on a hydrogen in the residue?
	       // The oxygen gets left. 
	       std::string chain_id(res->GetChainID());
	       int resno = res->GetSeqNum();
	       std::string altloc(at->altLoc);
	       std::string inscode(at->GetInsCode());
	       graphics_info_t::molecules[naii.imol].delete_residue_hydrogens(chain_id, resno, inscode, altloc);
	       delete_atom_by_atom_index(naii.imol, naii.atom_index,
					 destroy_delete_dialog_flag_by_ctrl_press);
	       pick_pending_flag = 0;
	    }
	 } else { 

	    if (show_symmetry) {
	       coot::Symm_Atom_Pick_Info_t symm_nearest_atom_index_info = symmetry_atom_pick(); 

	       if (symm_nearest_atom_index_info.success == GL_TRUE) {
		  int im = symm_nearest_atom_index_info.imol;
		  int index = symm_nearest_atom_index_info.atom_index;
		  CResidue *res = molecules[im].atom_sel.atom_selection[index]->residue;
		  std::string resname(res->name);
		  if (resname == "WAT" || resname == "HOH") { 
		     delete_atom_by_atom_index(im, index,
					       destroy_delete_dialog_flag_by_ctrl_press);
		     normal_cursor();
		     pick_pending_flag = 0;
		  }
	       }
	    }
	 }
      }

      // side chain
      if (g.delete_item_sidechain) {
	 pick_info naii = atom_pick(event);
	 if (naii.success == GL_TRUE) {
	    CResidue *res = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index]->residue;
	    std::string resname(res->name);
	    if (resname != "WAT" && resname != "HOH") {
	       normal_cursor();
	       int resno = res->GetSeqNum();
	       const char *ins_code = res->GetInsCode();
	       const char *chain_id = res->GetChainID();
	       delete_residue_sidechain(naii.imol, chain_id, resno, ins_code,
					destroy_delete_dialog_flag_by_ctrl_press);
	    }
	 }
      } 

      // residue
      if (g.delete_item_residue) {
	 pick_info naii = atom_pick(event);
	 if (naii.success == GL_TRUE) {
	    CResidue *res = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index]->residue;
	    std::string resname(res->name);
	    if (resname != "WAT" && resname != "HOH") {
	       normal_cursor();
	       delete_residue_by_atom_index(naii.imol, naii.atom_index,
					    destroy_delete_dialog_flag_by_ctrl_press);
	       pick_pending_flag = 0;
	    }
	 } else { 

	    if (show_symmetry) { 
	       coot::Symm_Atom_Pick_Info_t symm_nearest_atom_index_info = symmetry_atom_pick();

	       if (1) {
		  if (symm_nearest_atom_index_info.success == GL_TRUE) {
		     int im = symm_nearest_atom_index_info.imol;
		     int index = symm_nearest_atom_index_info.atom_index;
		     CResidue *res = molecules[im].atom_sel.atom_selection[index]->residue;
		     std::string resname(res->name);
		     if (resname != "WAT" && resname != "HOH") {
			pick_pending_flag = 0;
			normal_cursor();
			int im = symm_nearest_atom_index_info.imol;
			delete_residue_by_atom_index(im, index,
						     destroy_delete_dialog_flag_by_ctrl_press);
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
		     pick_pending_flag = 0;
		     normal_cursor();
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

	       // save a residue spec for the first residue
	       CAtom *at = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index];
	       g.delete_item_residue_zone_1 = coot::residue_spec_t(at->residue);
	       // now pick another atom
	       g.delete_item_residue_zone = 2;
	    } else {

	       // This was the second atom picked
	       CAtom *at = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index];
	       coot::residue_spec_t res2(at->residue);
	       // This is a graphics_info_t function, unlike the other
	       // delete functions.  This is because we use
	       // coot::residue_spec_t which can't be introduced int
	       // c-interface.h
	       delete_residue_range(naii.imol, g.delete_item_residue_zone_1, res2);
	       pick_pending_flag = 0;

	       // normal_cursor(); not necessarily, the dialog may not close
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
	 const CResidue *res_p = g.molecules[naii.imol].atom_sel.atom_selection[naii.atom_index]->GetResidue();
	 std::string chain_id(g.molecules[naii.imol].atom_sel.atom_selection[naii.atom_index]->GetChainID());
	 CAtom *at = g.molecules[naii.imol].atom_sel.atom_selection[naii.atom_index];

	 // we check on term_type inside the function
	 //
	 watch_cursor();
	 short int add_it_now_flag = g.add_terminal_residue_immediate_addition_flag;
	 if (add_terminal_residue_do_post_refine) {
	    add_it_now_flag = 1;
	 }
	 CResidue *r = (CResidue *) res_p;
	 if (!coot::util::is_nucleotide(r)) {
	    // const CResidue * is pain:
	    g.execute_add_terminal_residue(naii.imol,
					   term_type,
					   res_p,
					   chain_id,
					   g.add_terminal_residue_type,  // eg. "ALA" or "UNK"
					   add_it_now_flag);
	 } else {
	    g.execute_simple_nucleotide_addition(naii.imol, term_type, r, chain_id);
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
   return state; // ignored, currently
} 


void
graphics_info_t::check_if_in_db_main_define(GdkEventButton *event) { 

   graphics_info_t g;
   if (g.in_db_main_define) { 
      pick_info naii = atom_pick(event);
      if (naii.success == GL_TRUE) { 
	 molecules[naii.imol].add_to_labelled_atom_list(naii.atom_index);
	 if (g.in_db_main_define == 1) {
	    g.db_main_atom_index_1 = naii.atom_index;
	    g.in_db_main_define = 2;
	    g.db_main_imol = naii.imol;
	 } else {
	    
	    if (g.in_db_main_define == 2) {
	       if (naii.imol == g.db_main_imol) { 
		  g.db_main_atom_index_2 = naii.atom_index;
		  
		  // now do it
		  //
		  watch_cursor();
		  g.execute_db_main();
	       } else {
		  std::cout << "db-main: that atom was ";
		  std::cout << "not in the same molecule as the "
			    << "previous atom" << std::endl;
		  std::cout << "Cancelling selection" << std::endl;
	       }
	    }
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
	 CResidue *r = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index]->residue;
	 if (coot::util::is_nucleotide(r)) {
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

	    CAtom *atom_p = molecules[imol].atom_sel.atom_selection[naii.atom_index];
	    
	    int resno = atom_p->GetSeqNum();
	    std::string inscode(atom_p->GetInsCode());
	    std::string chain(atom_p->GetChainID());
	    std::string altloc(atom_p->altLoc);
	    
	    float f = molecules[imol].auto_fit_best_rotamer(resno, altloc, inscode, chain, imol_map, graphics_info_t::rotamer_fit_clash_flag, graphics_info_t::rotamer_lowest_probability);
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

	 // bla bla... more stuff here.

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
	 CAtom *at = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index];
	 int resno = at->GetSeqNum();
	 CResidue *residue = at->residue;
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
	    std::pair<CMMDBManager *, int> mp = 
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
	 statusbar_text(s);
	 
	 in_180_degree_flip_define = 0;
	 pick_pending_flag = 0;
	 normal_cursor();
 	 model_fit_refine_unactive_togglebutton("model_refine_dialog_do_180_degree_sidechain_flip_togglebutton");
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
graphics_info_t::check_if_in_save_symmetry_define(GdkEventButton *event) {

   if (in_save_symmetry_define) {

      coot::Symm_Atom_Pick_Info_t naii = symmetry_atom_pick();

      if (naii.success == GL_TRUE) { 

	 in_save_symmetry_define = 0;
	 GtkWidget *w = create_save_symmetry_coords_fileselection();

	 // attach symmetry info to the widget for use in the OK
	 // button callback:
	 coot::Symm_Atom_Pick_Info_t *save_pick_info = new coot::Symm_Atom_Pick_Info_t;
	 *save_pick_info = naii;

	 gtk_object_set_user_data(GTK_OBJECT(w), save_pick_info);

// 	 std::string filename = "molecule-";
// 	 filename += int_to_string(naii.imol);
	 int imol = naii.imol;
	 std::string filename = molecules[imol].name_for_display_manager();
	 filename += "-symmetry-";
	 filename += int_to_string(naii.symm_trans.isym());
	 filename += ".pdb";

	 gtk_file_selection_set_filename(GTK_FILE_SELECTION(w),
					 filename.c_str());
	 
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
	 CAtom *at = molecules[naii.imol].atom_sel.atom_selection[naii.atom_index];
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
	 CAtom *at = molecules[imol].atom_sel.atom_selection[naii.atom_index];
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


