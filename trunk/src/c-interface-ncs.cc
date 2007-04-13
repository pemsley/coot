/* src/c-interface-ncs.cc
 * 
 * Copyright 2004, 2005 The University of York
 * Author: Paul Emsley 
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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



#include <stdlib.h>
#include <iostream>

#if defined _MSC_VER
#include <windows.h>
#endif

#ifdef USE_GUILE
#include <guile/gh.h>
#endif // USE_GUILE

#ifdef USE_PYTHON
#include "Python.h"
#endif // USE_PYTHON

 
#include "globjects.h" //includes gtk/gtk.h

#include "callbacks.h"
#include "interface.h" // now that we are moving callback
		       // functionality to the file, we need this
		       // header since some of the callbacks call
		       // functions are built by glade.

#include <vector>
#include <string>

#include "mmdb_manager.h"
#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb-crystal.h"

#include "graphics-info.h"
#include "atom-utils.h" // asc_to_graphics

#include "c-interface.h"
#include "cc-interface.hh"

int add_strict_ncs_matrix(int imol,
			  const char *this_chain_id,
			  const char *target_chain_id,
			  float m11, float m12, float m13, 
			  float m21, float m22, float m23, 
			  float m31, float m32, float m33, 
			  float t1,  float t2,  float t3)
{

   int istat = 0;
   if (is_valid_model_molecule(imol)) {
      clipper::Mat33<double> m(m11, m12, m13,
			       m21, m22, m23,
			       m31, m32, m33);
      clipper::Coord_orth c(t1, t2, t3);
      clipper::RTop_orth rtop(m, c);

      coot::coot_mat44 cm44;
      // The orientation is a guess: FIXME or testme.
      cm44.m[0].v4[0] = m11;
      cm44.m[0].v4[1] = m12;
      cm44.m[0].v4[2] = m13;
      cm44.m[1].v4[0] = m21;
      cm44.m[1].v4[1] = m22;
      cm44.m[1].v4[2] = m23;
      cm44.m[2].v4[0] = m31;
      cm44.m[2].v4[1] = m32;
      cm44.m[2].v4[2] = m33;
      // translation
      cm44.m[0].v4[3] = t1;
      cm44.m[1].v4[3] = t2;
      cm44.m[2].v4[3] = t3;
      // sensibles
      cm44.m[3].v4[0] = 0.0;
      cm44.m[3].v4[1] = 0.0;
      cm44.m[3].v4[2] = 0.0;
      cm44.m[3].v4[3] = 1.0;
      
      istat = 1;
      std::string tch = target_chain_id;
      std::string chain_id = this_chain_id;

      graphics_info_t::molecules[imol].add_strict_ncs_matrix(chain_id,
							     tch,
							     cm44);

      graphics_draw();
   }
   return istat;
}


GtkWidget *wrapped_create_ncs_maps_dialog() {

   GtkWidget *dialog = create_ncs_maps_dialog();
   short int diff_maps_only_flag = 0;
   int ifound;

   // Maps: 
   ifound = fill_ligands_dialog_map_bits_by_dialog_name(dialog,
							"ncs_maps_maps",
							diff_maps_only_flag);
   if (ifound == 0) {
      std::cout << "Error: you must have a difference map to analyse!" << std::endl;
      GtkWidget *none_frame = lookup_widget(dialog, "no_maps_frame");
      gtk_widget_show(none_frame);
   }

   // Models:
   short int have_ncs_flag = 1;
   ifound = fill_vbox_with_coords_options_by_dialog_name(dialog, "ncs_maps_models",
							 have_ncs_flag);
   if (ifound == 0) {
      std::cout << "You must have molecules with NCS to use this function\n";
      GtkWidget *none_frame = lookup_widget(dialog, "no_models_frame");
      gtk_widget_show(none_frame);
   }


   return dialog;
}

void add_ncs_matrix(int imol,
		    const char *chain_id, 
		    const char *target_chain_id, 
		    float m11, float m12, float m13, 
		    float m21, float m22, float m23, 
		    float m31, float m32, float m33, 
		    float t1,  float t2,  float t3) {

   if (is_valid_model_molecule(imol)) {
      clipper::Mat33<double> m(m11, m12, m13,
			       m21, m22, m23,
			       m31, m32, m33);
      clipper::Coord_orth c(t1, t2, t3);
      clipper::RTop_orth rtop(m, c);

      graphics_info_t::molecules[imol].add_ncs_ghost(chain_id, target_chain_id, rtop);
      graphics_draw();
   }
}

int show_strict_ncs_state(int imol) {

   if (is_valid_model_molecule(imol)) {
      return graphics_info_t::molecules[imol].show_strict_ncs_flag;
   }
   return 0;
}

void set_show_strict_ncs(int imol, int state) {

   short int istate = 0;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].show_strict_ncs_flag = state;
      istate = 1;
      graphics_draw();
   }
}


/* At what level of homology should we say that we can't see homology
   for NCS calculation? (default 0.8) */
void set_ncs_homology_level(float flev) {

   graphics_info_t::ncs_homology_level = flev;

}

/*  ----------------------------------------------------------------------- */
/*                  NCS                                                     */
/*  ----------------------------------------------------------------------- */
/* section NCS */
void set_draw_ncs_ghosts(int imol, int istate) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_show_ghosts(istate);
      graphics_draw();
   }
}

void set_ncs_ghost_bond_thickness(int imol, float f) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_ghost_bond_thickness(f);
      graphics_draw();
   }
}


void ncs_update_ghosts(int imol) {

   if (is_valid_model_molecule(imol)) {
      int incs = graphics_info_t::molecules[imol].update_ncs_ghosts();
      if (incs)
	 graphics_draw();
   }
}

GtkWidget *wrapped_create_ncs_control_dialog() {

   GtkWidget *w = create_ncs_control_dialog();

   for (int imol=0; imol<graphics_info_t::n_molecules; imol++)
      if (is_valid_model_molecule(imol))
	 graphics_info_t::molecules[imol].fill_ncs_control_frame(w);
   return w; 
}


void ncs_control_change_ncs_master_to_chain(int imol, int ichain) {

   if (is_valid_model_molecule(imol)) {
      std::vector<std::string> chain_ids =
	 coot::util::chains_in_molecule(graphics_info_t::molecules[imol].atom_sel.mol);
      if (ichain < chain_ids.size())
	 graphics_info_t::molecules[imol].set_ncs_master_chain(chain_ids[ichain]);
      graphics_draw();
   } 
}

void ncs_control_change_ncs_master_to_chain_update_widget(GtkWidget *w, int imol, int ichain) {

   if (is_valid_model_molecule(imol)) {
      ncs_control_change_ncs_master_to_chain(imol, ichain); 

      // Now we want to update the widget.  We need to change the sensitivity of
      // all the Chain check boxes in the dispaly ncs chain vbox.
      // 
      // We need to change to desensitve the chain that matches ichain.
      //
      // This is a -widget-work function
      graphics_info_t::molecules[imol].ncs_control_change_ncs_master_to_chain_update_widget(w, ichain);
   } 
}


void ncs_control_display_chain(int imol, int ichain, int state) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_display_ncs_ghost_chain(ichain, state);
      graphics_draw();
   }
}

// The NCS "Yes" was pressed, do we need to make the rtops for this
// molecule?  We do if there ncs_ghosts.size() > 0.
//
// If ncs_ghosts.size() == 0, then this widget should not be active.
// 
// 
void
make_ncs_ghosts_maybe(int imol) {

   // c.f. fill_bond_parameters_internals()
   
   // int imol = graphics_info_t::bond_parameters_molecule; old way.  broken.
   
   if (is_valid_model_molecule(imol)) {  // it should be!
      if (graphics_info_t::molecules[imol].has_ncs_p()) {
	 if (graphics_info_t::molecules[imol].ncs_ghosts_have_rtops_p() == 0) {
	    // This popup causes a crash on mac...
	    // GtkWidget *w = popup_window("Calculating NCS pair transformations...");
	    graphics_info_t::molecules[imol].fill_ghost_info(1, graphics_info_t::ncs_homology_level);
	    // gtk_widget_destroy(w);
	 }
      }
   } 
}

// 
int make_dynamically_transformed_ncs_maps(int imol_model, int imol_map) {

   int nmaps=0;

   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {

	 nmaps = graphics_info_t::molecules[imol_model].make_dynamically_transformed_maps(imol_map, graphics_info_t::ncs_maps_do_average_flag, graphics_info_t::ncs_homology_level);
	 
      } else {
	 std::cout << "WARNING:: molecule number " << imol_map
		   << " is not a valid map molecule\n";
      }
   } else {
      std::cout << "WARNING:: molecule number " << imol_model
		<< " is not a valid model molecule\n";
   }

   if (nmaps > 0)
      graphics_draw();
   
   return nmaps;
}

// Also this should go elsewhere
// 
int make_dynamically_transformed_ncs_maps_by_widget(GtkWidget *dialog) {

   // Decode the options and call the above function:

   int imol_map = -1;
   int imol_coords = -1;
 
   GtkWidget *map_button;
   short int found_active_button_for_map = 0;
   for (int imol=0; imol<graphics_info_t::n_molecules; imol++) {
      if (graphics_info_t::molecules[imol].has_map()) {
	 std::string map_str = "ncs_maps_maps_radiobutton_";
	 map_str += graphics_info_t::int_to_string(imol);
	 map_button = lookup_widget(dialog, map_str.c_str());
	 if (map_button) {
	    if (GTK_TOGGLE_BUTTON(map_button)->active) { 
	       imol_map = imol;
	       found_active_button_for_map = 1;
	       break;
	    }
	 } else {
	    std::cout << "WARNING:: (error) " << map_str << " map button not found in "
		      << "make_dynamically_transformed_ncs_maps_by_widget" << std::endl;
	 }
      }
   }

   GtkWidget *coords_button;
   short int found_active_button_for_coords = 0;
   for (int imol=0; imol<graphics_info_t::n_molecules; imol++) {
      if (graphics_info_t::molecules[imol].has_model()) {
	 if (graphics_info_t::molecules[imol].has_ncs_p()) {
	    std::string coords_str = "ncs_maps_models_radiobutton_";
	    coords_str += graphics_info_t::int_to_string(imol);
	    coords_button = lookup_widget(dialog, coords_str.c_str());
	    if (coords_button) {
	       if (GTK_TOGGLE_BUTTON(coords_button)->active) {
		  imol_coords = imol;
		  found_active_button_for_coords = 1;
	       }
	    } else {
	       std::cout << "WARNING:: (error) " << coords_str
			 << " coords button not found in "
			 << "make_dynamically_transformed_ncs_maps_by_widget"
			 << std::endl;
	    }
	 }
      }
   }

   if (!found_active_button_for_coords) {
      std::cout << "You need to define a set of coordinates for NCS maping\n";
   } else {
      if (!found_active_button_for_map) {
	 std::cout << "You need to define a map for NCS maping\n";
      } else {
	 make_dynamically_transformed_ncs_maps(imol_coords, imol_map);
      }
   }

   return 0;
}


