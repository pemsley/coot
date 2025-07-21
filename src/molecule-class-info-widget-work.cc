/* src/molecule-class-info-widget-work.cc
 *
 * Copyright 2005, 2006 by The University of York
 * Author: Paul Emsley
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
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street, 02110-1335, USA.
 */


#ifdef _MSC_VER
#include <windows.h>
#endif

#include <vector>
#include "molecule-class-info.h"

#include <mmdb2/mmdb_manager.h>
#include "coords/Cartesian.hh"
#include "coords/mmdb-extras.h"
#include "coords/mmdb-crystal.h"

#include "utils/coot-utils.hh"
// #include "globjects.h" // for rotate_rgb.  That should be a utility
  			// function, not in globjects.hh

#include "widget-from-builder.hh"
#include "c-interface.h"
#include "c-interface-gtk-widgets.h"

void
molecule_class_info_t::update_map_colour_menu_maybe(int imol) {
   // or maybe not.
}




// extern "C" G_MODULE_EXPORT
void
on_ncs_controller_molecule_n_display_ncs_checkbutton_toggled_gtkbuilder_callback
                                        (GtkCheckButton *checkbutton,
                                         gpointer         user_data) {

  int imol = GPOINTER_TO_INT(user_data);
  int state = 0;
  if (gtk_check_button_get_active(checkbutton)) {
    state = 1;
    make_ncs_ghosts_maybe(imol);
  }
  /*    printf("NCS_controller Display NCS ghosts for imol %d %d\n", imol, state); */
  set_draw_ncs_ghosts(imol, state);
}


void
on_ncs_controller_molecule_n_display_chain_ich_checkbutton_toggled_gtkbuilder_callback
                                        (GtkCheckButton *checkbutton,
                                         gpointer         user_data) {
   int imol_chain = GPOINTER_TO_INT(user_data);
   int imol = imol_chain/1000;
   int ich = imol_chain - imol*1000;
   int state = 0;
   if (gtk_check_button_get_active(checkbutton)) {
     state = 1;
   }
   printf("\nNCS_controller display chain toggled for imol %d chain %d state %d\n",
	  imol, ich, state);
   ncs_control_display_chain(imol, ich, state);
}


void
on_ncs_controller_ncs_master_chain_ich_radiobutton_toggled_gtkbuilder_callback(GtkCheckButton *checkbutton,
                                                                               gpointer        user_data)
{
   GtkWidget *w = widget_from_builder("ncs_control_dialog");
   int imol_chain = GPOINTER_TO_INT(user_data);
   int imol = imol_chain/1000;
   int ich = imol_chain - imol*1000;
   /*    printf("==== DEBUG:: chain raiobutton toggled: imol %d ich %d active-state: %d \n",  */
   /* 	  imol, gtk_toggle_button_get_active(ich, togglebutton)); */
   if (gtk_check_button_get_active(checkbutton)) {
      /*      printf("NCS_controller_ncs_master_chain_ich_radiobutton_toggled on for imol %d %d %d\n",  */
      /* 	    imol_chain, imol, ich); */

      /*      ncs_control_change_ncs_master_to_chain(imol, ich); (done in the following function) */

      ncs_control_change_ncs_master_to_chain_update_widget(w, imol, ich);
   }
}


void
on_molecule_0_checkbutton_toggled_gtkbuilder_callback(GtkCheckButton *checkbutton,
                                                      gpointer        user_data) {

  int imol = GPOINTER_TO_INT(user_data);
  if (gtk_check_button_get_active(checkbutton))
    set_show_symmetry_molecule(imol, 1);
  else
    set_show_symmetry_molecule(imol, 0);

}

void
on_colour_symm_std_molecule_0_toggled_gtkbuilder_callback(GtkCheckButton *checkbutton,
                                                          gpointer        user_data) {

  int imol = GPOINTER_TO_INT(user_data);
  if (gtk_check_button_get_active(checkbutton)) {
    set_symmetry_colour_by_symop(imol, 0);
    set_symmetry_molecule_rotate_colour_map(imol, 0);
  }
}

void
on_display_sphere_radiobutton_molecule_0_toggled_gtkbuilder_callback(GtkCheckButton *checkbutton,
                                                                     gpointer        user_data) {

  int imol = GPOINTER_TO_INT(user_data);
  if (gtk_check_button_get_active(checkbutton)) {
     set_symmetry_whole_chain(imol, 0);
     symmetry_as_calphas(imol, 0); /* does an update_symmetry() */
  }
}


void
on_display_all_radiobutton_molecule_0_toggled_gtkbuilder_callback
                                        (GtkCheckButton *checkbutton,
                                        gpointer          user_data) {

  int imol = GPOINTER_TO_INT(user_data);
  if (gtk_check_button_get_active(checkbutton)) {
     symmetry_as_calphas(imol, 0);
     set_symmetry_whole_chain(imol, 1);
/*   } else { */
/*     symmetry_as_calphas(imol, 1); */
/*     printf("DEBUG:: all for molecule %d CA state 1\n", imol); */
   }
}

void
on_colour_symm_by_symop_molecule_0_toggled_gtkbuilder_callback
                                        (GtkCheckButton *togglebutton,
                                         gpointer         user_data) {

  int imol = GPOINTER_TO_INT(user_data);
  if (gtk_check_button_get_active(togglebutton)) {
    set_symmetry_molecule_rotate_colour_map(imol, 1); /* yes, I mean this */
    set_symmetry_colour_by_symop(imol, 1);
  }
}

void
on_display_CA_radiobutton_molecule_0_toggled_gtkbuilder_callback
                                        (GtkCheckButton *checkbutton,
                                        gpointer         user_data)
{

  int imol = GPOINTER_TO_INT(user_data);
  if (gtk_check_button_get_active(checkbutton)) {
     symmetry_as_calphas(imol, 1);
  }
}

void
on_colour_symm_by_molecule_molecule_0_toggled_gtkbuilder_callback(GtkCheckButton *checkbutton,
                                                                  gpointer         user_data) {

  int imol = GPOINTER_TO_INT(user_data);
  if (gtk_check_button_get_active(checkbutton)) {
    set_symmetry_colour_by_symop(imol, 0);
    set_symmetry_molecule_rotate_colour_map(imol, 1);
  }
}









void
molecule_class_info_t::handle_map_colour_change_rotate_difference_map(bool swap_difference_map_colours_flag) {

   // Usually (by default) the colours for the difference map are  green and red. Some people like red and
   // green. set_last_map_colour() calls this function and it is here that we decide on the second (negative
   // level) colour.
   float rotation_size = rotate_colour_map_for_difference_map/360.0;
   if (swap_difference_map_colours_flag)
      rotation_size = (360.0 - rotate_colour_map_for_difference_map)/360.0;

   // std::vector<float> rgb_new = rotate_rgb(orig_colours, rotation_size);
   // map_colour_negative_level.red   = rgb_new[0];
   // map_colour_negative_level.green = rgb_new[1];
   // map_colour_negative_level.blue  = rgb_new[2];

   coot::colour_holder ch(map_colour.red, map_colour.green, map_colour.blue);
   ch.rotate_by(rotation_size);
   map_colour_negative_level.red   = ch.red;
   map_colour_negative_level.green = ch.green;
   map_colour_negative_level.blue  = ch.blue;
}

void
molecule_class_info_t::handle_map_colour_change(GdkRGBA map_col_in,
                                                bool swap_difference_map_colours_flag,
                                                bool main_or_secondary,
                                                clipper::Coord_orth centre,
                                                float radius) {

   if (false)
      std::cout << "debug:: handle_map_colour_change() handle change to map colour "
                << map_col_in.red << " "
                << map_col_in.green << " "
                << map_col_in.blue << std::endl;

   map_colour = map_col_in;

   // input map colours are now in the range 0 to 1
   // map_colour.red   = map_col_in.red  /65535.0;
   // map_colour.green = map_col_in.green/65535.0;
   // map_colour.blue  = map_col_in.blue /65535.0;

   if (xmap_is_diff_map) {
      // std::cout << "calling handle_map_colour_change_rotate_difference_map() map colour for " << imol_no
      // << " with swap flag " << swap_difference_map_colours_flag << std::endl;
      handle_map_colour_change_rotate_difference_map(swap_difference_map_colours_flag);
   }

   // ideally, just change the colour buffer, but this will do for now.

   // clipper::Coord_orth centre = graphics_info_t::get_rotation_centre_co();
   // float radius = graphics_info_t::box_radius_xray;
   setup_glsl_map_rendering(centre, radius);

   // main 0: secondary: 1
   // compile_density_map_display_list(main_or_secondary);
}

// symmetry control
//
// We create a frame and add it to the viewport that's passed.  It is used to fill the
// symmetry control widget (requested by Frank von Delft)
//
// Perhaps this should be in a file molecule-class-info-widget-work.cc
void
molecule_class_info_t::fill_symmetry_control_frame(GtkWidget *symmetry_controller_dialog) const {

   std::string s = "Molecule ";
   std::string imol_str = coot::util::int_to_string(imol_no);
   std::string molecule_n = "molecule_";
   molecule_n += imol_str;
   s          += imol_str;
   s          += " ";
   s          += name_for_display_manager();

   GtkWidget *molecule_0_frame;
   GtkWidget *vbox168;
   GtkWidget *molecule_0_checkbutton;
   GtkWidget *frame162;
   GtkWidget *table4;
   GSList *symm_display_mol_0_gr_group = NULL;
   GtkWidget *display_sphere_radiobutton_molecule_0 = nullptr;
   GtkWidget *display_all_radiobutton_molecule_0 = nullptr;
   GtkWidget *display_CA_radiobutton_molecule_0 = nullptr;
   GSList *symm_colour_mol_0_gr_group = NULL;
   GtkWidget *colour_symm_std_molecule_0 = nullptr;
   GtkWidget *colour_symm_by_symop_molecule_0 = nullptr;
   GtkWidget *colour_symm_by_molecule_molecule_0 = nullptr;
   GtkWidget *symmetry_control_vbox = nullptr;

   // symmetry_control_vbox = lookup_widget(symmetry_controller_dialog, "symmetry_controller_vbox");
   symmetry_control_vbox = widget_from_builder("symmetry_controller_vbox");

   molecule_0_frame = gtk_frame_new (s.c_str());
   // gtk_widget_ref (molecule_0_frame);
   std::string t = molecule_n + "_frame";
   g_object_set_data_full (G_OBJECT (symmetry_controller_dialog),
			     t.c_str(),
			     molecule_0_frame, NULL);
   gtk_box_append(GTK_BOX(symmetry_control_vbox), molecule_0_frame);

   vbox168 = gtk_box_new (GTK_ORIENTATION_VERTICAL, 0);
   // gtk_widget_ref (vbox168);
   g_object_set_data_full (G_OBJECT (symmetry_controller_dialog), "vbox168", vbox168, NULL);
   gtk_widget_set_visible (vbox168, TRUE);
   gtk_frame_set_child(GTK_FRAME(molecule_0_frame), vbox168);

   molecule_0_checkbutton = gtk_check_button_new_with_label (" Show Symmetry?");
   // gtk_widget_ref (molecule_0_checkbutton);
   std::string molecule_n_checkbutton = molecule_n + "_checkbutton";
   g_object_set_data_full (G_OBJECT (symmetry_controller_dialog),
			     molecule_n_checkbutton.c_str(),
			     molecule_0_checkbutton, NULL);
   gtk_widget_set_visible (molecule_0_checkbutton, TRUE);

   gtk_box_append (GTK_BOX (vbox168), molecule_0_checkbutton);
   if (show_symmetry)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(molecule_0_checkbutton), TRUE);

   frame162 = gtk_frame_new ("Display Options");
   // gtk_widget_ref (frame162);
   g_object_set_data_full (G_OBJECT (symmetry_controller_dialog), "frame162", frame162, NULL);
   gtk_widget_set_visible (frame162, TRUE);
   gtk_box_append (GTK_BOX (vbox168), frame162);

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
      // 20220528-PE FIXME table to grid
#else
   // table4 = gtk_table_new (3, 2, FALSE);
   table4 = gtk_grid_new ();
   // gtk_widget_ref (table4);
   g_object_set_data_full (G_OBJECT (symmetry_controller_dialog), "table4", table4, NULL);
   gtk_widget_set_visible (table4, TRUE);
   gtk_container_add (GTK_CONTAINER (frame162), table4);

   display_sphere_radiobutton_molecule_0 = gtk_radio_button_new_with_label (symm_display_mol_0_gr_group, "Display Sphere");
   symm_display_mol_0_gr_group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (display_sphere_radiobutton_molecule_0));
   // gtk_widget_ref (display_sphere_radiobutton_molecule_0);
   s = "display_sphere_radiobutton_" + molecule_n;
   g_object_set_data_full (G_OBJECT (symmetry_controller_dialog),
			     s.c_str(),
			     display_sphere_radiobutton_molecule_0, NULL);
   gtk_widget_set_visible (display_sphere_radiobutton_molecule_0, TRUE);

   gtk_table_attach (GTK_TABLE (table4), display_sphere_radiobutton_molecule_0, 0, 1, 0, 1,
		     (GtkAttachOptions) (GTK_FILL),
		     (GtkAttachOptions) (0), 0, 0);
#endif

#if (GTK_MAJOR_VERSION >= 4)
   // 20220602-PE FIXME radio buttons
#else
   display_all_radiobutton_molecule_0 = gtk_radio_button_new_with_label (symm_display_mol_0_gr_group,
									 "Display Near Chains");
   symm_display_mol_0_gr_group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (display_all_radiobutton_molecule_0));
#endif
   // gtk_widget_ref (display_all_radiobutton_molecule_0);

   // set display_all_radiobutton_
   if (symmetry_whole_chain_flag)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(display_all_radiobutton_molecule_0), TRUE);

   s = "display_all_radiobutton_" + molecule_n;
   g_object_set_data_full (G_OBJECT (symmetry_controller_dialog),
			     s.c_str(),
			     display_all_radiobutton_molecule_0, NULL);
   gtk_widget_set_visible (display_all_radiobutton_molecule_0, TRUE);

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)

#else
   gtk_table_attach (GTK_TABLE (table4), display_all_radiobutton_molecule_0, 0, 1, 1, 2,
		     (GtkAttachOptions) (GTK_FILL),
		     (GtkAttachOptions) (0), 0, 0);
#endif

#if (GTK_MAJOR_VERSION >= 4)
   // 20220602-PE FIXME radio buttons
#else
   display_CA_radiobutton_molecule_0 = gtk_radio_button_new_with_label (symm_display_mol_0_gr_group, "Display as CAs");
   symm_display_mol_0_gr_group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (display_CA_radiobutton_molecule_0));
#endif

   // gtk_widget_ref (display_CA_radiobutton_molecule_0);

   // set display_CA_radiobutton_
   if (symmetry_as_calphas)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(display_CA_radiobutton_molecule_0), TRUE);

   s = "display_CA_radiobutton_" + molecule_n;
   g_object_set_data_full (G_OBJECT (symmetry_controller_dialog),
			     s.c_str(),
			     display_CA_radiobutton_molecule_0, NULL);
   gtk_widget_set_visible (display_CA_radiobutton_molecule_0, TRUE);

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)

#else
   gtk_table_attach (GTK_TABLE (table4), display_CA_radiobutton_molecule_0, 0, 1, 2, 3,
		     (GtkAttachOptions) (GTK_FILL),
		     (GtkAttachOptions) (0), 0, 0);
#endif

#if (GTK_MAJOR_VERSION >= 4)
   // 20220602-PE FIXME radio buttons
#else
   colour_symm_std_molecule_0 = gtk_radio_button_new_with_label (symm_colour_mol_0_gr_group, "Standard Colouring");
   symm_colour_mol_0_gr_group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (colour_symm_std_molecule_0));
#endif
   // gtk_widget_ref (colour_symm_std_molecule_0);

   s = "colour_symm_std_" + molecule_n;
   g_object_set_data_full (G_OBJECT (symmetry_controller_dialog),
			     s.c_str(),
			     colour_symm_std_molecule_0, NULL);
   gtk_widget_set_visible (colour_symm_std_molecule_0, TRUE);

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)

#else
   gtk_table_attach (GTK_TABLE (table4), colour_symm_std_molecule_0, 1, 2, 0, 1,
		     (GtkAttachOptions) (GTK_FILL),
		     (GtkAttachOptions) (0), 0, 0);
#endif

   // set the colour radiobutton
   if (symmetry_colour_by_symop_flag == 0 && symmetry_rotate_colour_map_flag == 0)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(colour_symm_std_molecule_0), TRUE);

#if (GTK_MAJOR_VERSION >= 4)
   // 20220602-PE FIXME radio buttons
#else
   colour_symm_by_symop_molecule_0 = gtk_radio_button_new_with_label (symm_colour_mol_0_gr_group,
								      "Colour by Symop");
   symm_colour_mol_0_gr_group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (colour_symm_by_symop_molecule_0));
#endif
   // gtk_widget_ref (colour_symm_by_symop_molecule_0);

   s = "colour_symm_by_symop_" + molecule_n;
   g_object_set_data_full (G_OBJECT (symmetry_controller_dialog),
			     s.c_str(),
			     colour_symm_by_symop_molecule_0, NULL);
   gtk_widget_set_visible (colour_symm_by_symop_molecule_0, TRUE);

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)

#else
   gtk_table_attach (GTK_TABLE (table4), colour_symm_by_symop_molecule_0, 1, 2, 1, 2,
		     (GtkAttachOptions) (GTK_FILL),
		     (GtkAttachOptions) (0), 0, 0);
#endif

   // set the colour radiobutton
   if (symmetry_colour_by_symop_flag)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(colour_symm_by_symop_molecule_0), TRUE);


#if (GTK_MAJOR_VERSION >= 4)
   // 20220602-PE FIXME radio buttons
#else
   colour_symm_by_molecule_molecule_0 = gtk_radio_button_new_with_label (symm_colour_mol_0_gr_group, "Colour by Molecule");
   symm_colour_mol_0_gr_group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (colour_symm_by_molecule_molecule_0));
#endif

   // gtk_widget_ref (colour_symm_by_molecule_molecule_0);
   s = "colour_symm_by_molecule_" + molecule_n;
   g_object_set_data_full (G_OBJECT (symmetry_controller_dialog),
			     s.c_str(),
			     colour_symm_by_molecule_molecule_0, NULL);
   gtk_widget_set_visible (colour_symm_by_molecule_molecule_0, TRUE);

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)

#else
   gtk_table_attach (GTK_TABLE (table4), colour_symm_by_molecule_molecule_0, 1, 2, 2, 3,
		     (GtkAttachOptions) (GTK_FILL),
		     (GtkAttachOptions) (0), 0, 0);
#endif

   // set the colour radiobutton
   if (symmetry_rotate_colour_map_flag)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(colour_symm_by_molecule_molecule_0), TRUE);

   // 20220405-PE Porting to GTK4:
   // these callbacks are in glade-callbacks.cc - that's unusual and required glade-callbacks.hh

   g_signal_connect (G_OBJECT (molecule_0_checkbutton), "toggled",
		       G_CALLBACK (on_molecule_0_checkbutton_toggled_gtkbuilder_callback),
		       GINT_TO_POINTER(imol_no));
   g_signal_connect (G_OBJECT (display_sphere_radiobutton_molecule_0), "toggled",
		       G_CALLBACK (on_display_sphere_radiobutton_molecule_0_toggled_gtkbuilder_callback),
		       GINT_TO_POINTER(imol_no));
   g_signal_connect (G_OBJECT (display_all_radiobutton_molecule_0), "toggled",
		       G_CALLBACK (on_display_all_radiobutton_molecule_0_toggled_gtkbuilder_callback),
		       GINT_TO_POINTER(imol_no));
   g_signal_connect (G_OBJECT (display_CA_radiobutton_molecule_0), "toggled",
		       G_CALLBACK (on_display_CA_radiobutton_molecule_0_toggled_gtkbuilder_callback),
		       GINT_TO_POINTER(imol_no));
   g_signal_connect (G_OBJECT (colour_symm_std_molecule_0), "toggled",
		       G_CALLBACK (on_colour_symm_std_molecule_0_toggled_gtkbuilder_callback),
		       GINT_TO_POINTER(imol_no));
   g_signal_connect (G_OBJECT (colour_symm_by_symop_molecule_0), "toggled",
		       G_CALLBACK (on_colour_symm_by_symop_molecule_0_toggled_gtkbuilder_callback),
		       GINT_TO_POINTER(imol_no));
   g_signal_connect (G_OBJECT (colour_symm_by_molecule_molecule_0), "toggled",
		       G_CALLBACK (on_colour_symm_by_molecule_molecule_0_toggled_gtkbuilder_callback),
		       GINT_TO_POINTER(imol_no));

   gtk_widget_set_visible(molecule_0_frame, TRUE);
}

// NCS control
void
molecule_class_info_t::fill_ncs_control_frame(GtkWidget *ncs_control_dialog) const {

   if (atom_sel.n_selected_atoms > 0) {
      if (ncs_ghosts.size() > 0) {
         fill_ncs_control_frame_internal(ncs_control_dialog);
      }
   }
}

#include "graphics-info.h" // 20220315-PE becaause we want clear_out_container().
                           // (no, it's not a good arangement)

// NCS control
void
molecule_class_info_t::fill_ncs_control_frame_internal(GtkWidget *ncs_control_dialog) const {

   if (! ncs_control_dialog) return;

   GtkWidget *ncs_control_vbox = widget_from_builder("ncs_control_vbox");
   graphics_info_t::clear_out_container(ncs_control_vbox);

   std::string m("Molecule ");
   m += dotted_chopped_name();

   GtkWidget *frame = gtk_frame_new(m.c_str());
   GtkWidget *box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 6);
   GtkWidget *display_checkbutton = gtk_check_button_new_with_label("Display Non-crystallographic Ghosts");
   GtkWidget *sep = gtk_separator_new(GTK_ORIENTATION_HORIZONTAL);
   GtkWidget *grid = gtk_grid_new();
   GtkWidget *ld      = gtk_label_new("Displayed Chains");
   GtkWidget *lmaster = gtk_label_new("NCS Master Chain");

   if (show_ghosts_flag)
      gtk_check_button_set_active(GTK_CHECK_BUTTON(display_checkbutton), TRUE);

   gtk_box_append(GTK_BOX(ncs_control_vbox), frame);
   gtk_frame_set_child(GTK_FRAME(frame), box);
   gtk_box_append(GTK_BOX(box), display_checkbutton);
   gtk_box_append(GTK_BOX(box), sep);
   gtk_box_append(GTK_BOX(box), grid);
   gtk_grid_attach(GTK_GRID(grid), ld,      0, 0, 1, 1);
   gtk_grid_attach(GTK_GRID(grid), lmaster, 1, 0, 1, 1);

   gtk_widget_set_sensitive(grid, FALSE); // insensitive when Display button is off

   gtk_widget_set_margin_start(frame, 6);
   gtk_widget_set_margin_end(frame, 6);
   gtk_widget_set_margin_top(frame, 2);
   gtk_widget_set_margin_bottom(frame, 2);
   gtk_widget_set_margin_start(ld, 6);
   gtk_widget_set_margin_end(ld, 10);
   gtk_widget_set_margin_start(lmaster, 6);
   gtk_widget_set_margin_end(lmaster, 6);

   auto display_checkbutton_callback = +[] (GtkCheckButton *cb, gpointer data) {

      GtkWidget *grid = static_cast<GtkWidget *>(data);
      int n_chains = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(data), "n_chains"));
      int imol     = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(data), "imol"));
      if (gtk_check_button_get_active(cb)) {
         gtk_widget_set_sensitive(grid, TRUE);
         if (n_chains > 1) {
            for (int ii=0; ii<n_chains; ii++) {
               GtkWidget *cbd = gtk_grid_get_child_at(GTK_GRID(grid), 0, ii+1);
               GtkWidget *cbm = gtk_grid_get_child_at(GTK_GRID(grid), 1, ii+1);
               if (gtk_check_button_get_active(GTK_CHECK_BUTTON(cbd))) {
                  // show NCS for chain ii
                  ncs_control_display_chain(imol, ii, true);
               } else {
                  ncs_control_display_chain(imol, ii, false);
               }
            }
         }
      } else {
         graphics_info_t::molecules[imol].set_show_ghosts(0);
         gtk_widget_set_sensitive(grid, FALSE);
         graphics_info_t::graphics_draw();
      }
   };
   g_signal_connect(G_OBJECT(display_checkbutton), "toggled",
                    G_CALLBACK(display_checkbutton_callback), grid);

   std::vector<std::string> v = coot::util::chains_in_molecule(atom_sel.mol);
   int n_chains = v.size();
   g_object_set_data(G_OBJECT(grid), "n_chains", GINT_TO_POINTER(n_chains));
   g_object_set_data(G_OBJECT(grid), "imol",     GINT_TO_POINTER(imol_no));
   GtkWidget *master_chain_group = nullptr;
   auto cbd_toggled = +[] (GtkCheckButton *checkbutton, gpointer data) {

      int imol    = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(checkbutton), "imol"));
      int i_chain = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(checkbutton), "i_chain"));
      GtkWidget *display_checkbutton = GTK_WIDGET(g_object_get_data(G_OBJECT(checkbutton), "display_checkbutton"));
      std::cout << "debug:: checkbutton " << checkbutton << std::endl;
      std::cout << "debug:: display_checkbutton " << display_checkbutton << std::endl;
      if (gtk_check_button_get_active(checkbutton)) {
         if (gtk_check_button_get_active(GTK_CHECK_BUTTON(display_checkbutton)))
            ncs_control_display_chain(imol, i_chain, true);
      } else {
         ncs_control_display_chain(imol, i_chain, false);
      }
   };
   auto cbm_toggled = +[] (GtkCheckButton *checkbutton, gpointer data) {

      std::cout << "do something with the master chain change " << checkbutton << std::endl;
      if (gtk_check_button_get_active(checkbutton)) {
         int i_chain = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(checkbutton), "i_chain"));
         int imol    = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(checkbutton), "imol"));
         std::cout << "set chain " << i_chain << " to master chain" << std::endl;
         GtkWidget *w = nullptr; // not used
         ::ncs_control_change_ncs_master_to_chain_update_widget(w, imol, i_chain);
      }
   };
   for (int ich=0; ich<n_chains; ich++) {
      std::string label = "Chain ";
      label += v[ich];
      GtkWidget *cbd = gtk_check_button_new_with_label(label.c_str());
      GtkWidget *cbm = gtk_check_button_new_with_label(label.c_str());
      g_object_set_data(G_OBJECT(cbd), "imol",    GINT_TO_POINTER(imol_no));
      g_object_set_data(G_OBJECT(cbd), "i_chain", GINT_TO_POINTER(ich));
      g_object_set_data(G_OBJECT(cbm), "imol",    GINT_TO_POINTER(imol_no));
      g_object_set_data(G_OBJECT(cbm), "i_chain", GINT_TO_POINTER(ich));
      g_object_set_data(G_OBJECT(cbd), "display_checkbutton", display_checkbutton);
      gtk_grid_attach(GTK_GRID(grid), cbd, 0, ich+1, 1, 1);
      gtk_grid_attach(GTK_GRID(grid), cbm, 1, ich+1, 1, 1);
      g_signal_connect(G_OBJECT(cbd), "toggled", G_CALLBACK(cbd_toggled), nullptr);
      g_signal_connect(G_OBJECT(cbm), "toggled", G_CALLBACK(cbm_toggled), nullptr);
      if (ich > 0)
         gtk_check_button_set_active(GTK_CHECK_BUTTON(cbd), TRUE);
      if (master_chain_group) {
         gtk_check_button_set_group(GTK_CHECK_BUTTON(cbm), GTK_CHECK_BUTTON(master_chain_group));
      } else {
         // the first chain is the master by default
         master_chain_group = cbm;
         gtk_check_button_set_active(GTK_CHECK_BUTTON(cbm), TRUE);
      }
   }
}


// widget is no longer used (it used to be needed for lookups)
void
molecule_class_info_t::ncs_control_change_ncs_master_to_chain_update_widget(GtkWidget *widget, int imaster) const {

   // Now we want to update the widget.  We need to change the sensitivity of
   // all the Chain check boxes in the dispaly ncs chain vbox.
   //
   // We need to change to desensitve the chain that matches ichain.
   //

   // First find imaster
   std::vector<std::string> chain_ids = coot::util::chains_in_molecule(atom_sel.mol);

   if (imaster != -1) {
      // GtkWidget *vbox = lookup_widget(w, "ncs_controller_molecule_n_display_chain_vbox");
      GtkWidget *vbox = widget_from_builder("ncs_controller_molecule_n_display_chain_vbox");
      std::string imol_str = coot::util::int_to_string(imol_no);
      for (unsigned int i=0; i<chain_ids.size(); i++) {
         std::string name = "ncs_controller_molecule_";
         name += imol_str;
         name += "_display_chain_";
         name += coot::util::int_to_string(i);
         name += "_checkbutton";
         // GtkWidget *checkbutton = lookup_widget(vbox, name.c_str());
         GtkWidget *checkbutton = 0;
         std::cout << "in ncs_control_change_ncs_master_to_chain_update_widget() set the checkbutton correctly" << std::endl;
         if (checkbutton) {
            if (int(i) == imaster) {
               gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), FALSE);
               gtk_widget_set_sensitive(checkbutton, FALSE);
            } else {
               gtk_widget_set_sensitive(checkbutton, TRUE);
               // ncs control turns on all chains when we change the master
               gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
            }
         }
      }
   }
}
