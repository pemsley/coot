/*
 * src/single-map-properties-dialog.cc
 *
 * Copyright 2021 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
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
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#include <string>
#include <iostream>

#include "utils/coot-utils.hh"
#include "single-map-properties-dialog.hh"
#include "graphics-info.h"
#include "c-interface-gtk-widgets.h"

std::pair<GtkWidget *, GtkBuilder *> create_single_map_properties_dialog_gtk3() {

   GtkWidget *single_map_properties_dialog = 0;
   GtkBuilder *builder = gtk_builder_new();

   std::string dir = coot::package_data_dir();
   // glade to ui
   std::string dir_ui = coot::util::append_dir_dir(dir, "ui");
   // std::string glade_file_name = "single-map-properties-dialog.glade";
   std::string ui_file_name = "single-map-properties-dialog-gtk4.ui";
   std::string ui_file_full = coot::util::append_dir_file(dir_ui, ui_file_name);
   if (coot::file_exists(ui_file_name))
      ui_file_full = ui_file_name;

   GError *error = NULL;
   guint add_from_file_status = gtk_builder_add_from_file(builder, ui_file_full.c_str(), &error);

   if (add_from_file_status) {
      single_map_properties_dialog = GTK_WIDGET(gtk_builder_get_object(builder, "single_map_properties_dialog"));
      // gtk_builder_connect_signals(builder, single_map_properties_dialog); // automatic now (gtk4)
   } else {
      std::cout << "ERROR:: create_single_map_properties_dialog_gtk3() failed to get builder file for single-map-properties dialog"
                << std::endl;
      std::cout << "ERROR:: " << ui_file_full << std::endl;
      std::cout << "ERROR::" << error->message << std::endl;
   }

   return std::make_pair(single_map_properties_dialog, builder);

}

void fill_single_map_properties_dialog_gtk3(std::pair<GtkWidget *, GtkBuilder *> w_and_b, int imol) {

   graphics_info_t g;
   if (! g.is_valid_map_molecule(imol)) return;

   GtkWidget *dialog   = w_and_b.first;
   GtkBuilder *builder = w_and_b.second;
   g_object_set_data(G_OBJECT(dialog), "imol", GINT_TO_POINTER(imol));

   auto widget_from_builder = [builder] (const std::string &wid) {
                                 return GTK_WIDGET(gtk_builder_get_object(builder, wid.c_str()));
                              };

   auto setup_contour_buttons_and_entry = [widget_from_builder, imol] (float contour_level, bool contour_by_rmsd_flag) {

                                                GtkWidget *cl_apply_button = widget_from_builder("single_map_properties_contour_level_apply_button");
                                                GtkWidget *cl_entry        = widget_from_builder("single_map_properties_contour_level_entry");
                                                GtkWidget *level_type_absolute_radiobutton = widget_from_builder("single_map_properties_absolute_radiobutton");
                                                GtkWidget *level_type_rmsd_radiobutton     = widget_from_builder("single_map_properties_rmsd_radiobutton");


                                                // we can't set radiobuttons to active false! (makes sense)
                                                if (contour_by_rmsd_flag)
                                                   gtk_check_button_set_active(GTK_CHECK_BUTTON(level_type_rmsd_radiobutton), TRUE);

                                                g_object_set_data(G_OBJECT(level_type_absolute_radiobutton), "contour_level_entry", cl_entry);
                                                std::string entry_text = coot::util::float_to_string_using_dec_pl(contour_level, 3);
                                                gtk_editable_set_text(GTK_EDITABLE(cl_entry), entry_text.c_str());
                                                g_object_set_data(G_OBJECT(cl_apply_button), "imol", GINT_TO_POINTER(imol));
                                                g_object_set_data(G_OBJECT(cl_apply_button), "contour_level_entry", cl_entry);
                                                g_object_set_data(G_OBJECT(cl_apply_button), "single_map_properties_absolute_radiobutton", level_type_absolute_radiobutton);

                                             };

   // 20230429-PE my brain is melting with boredom in this dialog - just get rid of sigma step size for now.
   //
   auto setup_step_size_widget = [widget_from_builder] (int imol) {

                                                GtkWidget *step_size_is_in_rmsd_checkbutton = widget_from_builder("single_map_properties_step_in_rmsd_checkbutton");
                                                GtkWidget *step_size_entry = widget_from_builder("single_map_properties_step_size_entry");

                                                g_object_set_data(G_OBJECT(step_size_entry),                  "imol", GINT_TO_POINTER(imol));
                                                g_object_set_data(G_OBJECT(step_size_is_in_rmsd_checkbutton), "imol", GINT_TO_POINTER(imol));

                                                g_object_set_data(G_OBJECT(step_size_entry), "step_size_checkbutton", step_size_is_in_rmsd_checkbutton);
                                                g_object_set_data(G_OBJECT(step_size_is_in_rmsd_checkbutton), "step_size_checkbutton", step_size_entry);

   };

   GtkWidget *cell_text = widget_from_builder("single_map_properties_cell_label");
   GtkWidget *spgr_text = widget_from_builder("single_map_properties_symmetry_label");
   GtkWidget *reso_text = widget_from_builder("single_map_properties_resolution_label");

   std::string cell_text_string;
   std::string spgr_text_string;
   std::string reso_text_string;

   std::string title = "Coot: Properties for Map " + std::to_string(imol);
   gtk_window_set_title(GTK_WINDOW(dialog), title.c_str());

   const clipper::Xmap<float> &xmap = graphics_info_t::molecules[imol].xmap;
   cell_text_string = graphics_info_t::molecules[imol].cell_text_with_embeded_newline();
   spgr_text_string = "   ";
   spgr_text_string += xmap.spacegroup().descr().symbol_hm();
   spgr_text_string += "  [";
   spgr_text_string += xmap.spacegroup().descr().symbol_hall();
   spgr_text_string += "]";
   float r = graphics_info_t::molecules[imol].data_resolution();
   if (r < 0) {
      r = 2.0 * xmap.cell().descr().a()/static_cast<float>(xmap.grid_sampling().nu());
      reso_text_string = " ";
      reso_text_string += coot::util::float_to_string(r);
   } else {
      reso_text_string = coot::util::float_to_string(r);
   }
   // now add the grid info to the reso text
   reso_text_string += " Å (by grid) ";
   clipper::Grid_sampling gs = xmap.grid_sampling();
   reso_text_string += coot::util::int_to_string(gs.nu()) + " ";
   reso_text_string += coot::util::int_to_string(gs.nv()) + " ";
   reso_text_string += coot::util::int_to_string(gs.nw());

   gtk_label_set_text(GTK_LABEL(cell_text), cell_text_string.c_str());
   gtk_label_set_text(GTK_LABEL(spgr_text), spgr_text_string.c_str());
   gtk_label_set_text(GTK_LABEL(reso_text), reso_text_string.c_str());

   // return;  OK

   // And now the map rendering style: transparent surface or standard lines:
   GtkWidget *rb_1  = widget_from_builder("display_map_style_as_lines_radiobutton");
   GtkWidget *rb_2  = widget_from_builder("display_map_style_surface_radiobutton");
   GtkWidget *map_opacity_scale = widget_from_builder("map_opacity_hscale");

   g_object_set_data(G_OBJECT(rb_1), "imol", GINT_TO_POINTER(imol));

   GtkWidget *level_type_radiobutton = widget_from_builder("single_map_properties_absolute_radiobutton");
   g_object_set_data(G_OBJECT(level_type_radiobutton), "imol", GINT_TO_POINTER(imol));

   const molecule_class_info_t &m = g.molecules[imol];

   if (! m.draw_it_for_map_standard_lines) {
      gtk_check_button_set_active(GTK_CHECK_BUTTON(rb_2), TRUE);
   } else {
      gtk_check_button_set_active(GTK_CHECK_BUTTON(rb_1), TRUE);
   }

   // return; OK

   g_object_set_data(G_OBJECT(map_opacity_scale), "imol", GINT_TO_POINTER(imol));
   GtkAdjustment *adjustment = gtk_range_get_adjustment(GTK_RANGE(map_opacity_scale));
   float op = m.density_surface_opacity;
   gtk_adjustment_set_value(adjustment, 100.0*op);


   // resurect this one day
   // GtkWidget *map_contour_frame = widget_from_builder("single_map_properties_map_histogram_frame");
   // GtkWidget *alignment = widget_from_builder("alignment_for_map_density_histogram");
   // fill_map_histogram_widget(imol, alignment);

   // return; OK

   float contour_level = m.contour_level;
   if (m.contour_by_sigma_flag)
      contour_level = contour_level / m.map_sigma();

   setup_contour_buttons_and_entry(contour_level, m.contour_by_sigma_flag);

   // resurect this one day
   /*  and now the skeleton buttons */
   // GtkWidget *frame = lookup_widget(dialog, "single_map_skeleton_frame");
   // set_on_off_single_map_skeleton_radio_buttons(frame, imol);
   /* contour by sigma step */

   // -------------------------------------------------------------------------------
   // colour
   // -------------------------------------------------------------------------------

   GtkWidget *colour_button = widget_from_builder("single_map_properties_colour_button");
   if (colour_button) {
      g_object_set_data(G_OBJECT(colour_button), "imol", GINT_TO_POINTER(imol));
      GdkRGBA map_colour = get_map_colour(imol);
      if (false)
         std::cout << "DEBUG:: ---------------- got colour "
                   << map_colour.red << " " << map_colour.green << " " << map_colour.blue << " " << map_colour.alpha << std::endl;
      if (true) { // 2026-02-15-PE get_map_colour multiplies - so now we divide
                  // (maybe just don't multiply?)
         map_colour.red   /= 65535.0; 
         map_colour.green /= 65535.0; 
         map_colour.blue  /= 65535.0; 
         map_colour.alpha /= 65535.0; 
      }
      gtk_color_chooser_set_rgba(GTK_COLOR_CHOOSER(colour_button), &map_colour);
   } else {
      std::cout << "ERROR:: --------------- no colour_button found!" << std::endl;
   }

   // -------------------------------------------------------------------------------
   // Specularity
   // -------------------------------------------------------------------------------

   GtkWidget *specularity_checkbutton = widget_from_builder("map_properties_dialog_specularity_state_checkbutton");
   if (specularity_checkbutton) {
      
      GtkWidget *strength_entry  = widget_from_builder("map_properties_dialog_specularity_strength_entry");
      GtkWidget *shininess_entry = widget_from_builder("map_properties_dialog_specularity_shininess_entry");

      float specular_strength = m.material_for_maps.specular_strength;
      float shininess = m.material_for_maps.shininess;

      g_object_set_data(G_OBJECT(specularity_checkbutton), "imol", GINT_TO_POINTER(imol));
      g_object_set_data(G_OBJECT(specularity_checkbutton), "strength_entry",   strength_entry);
      g_object_set_data(G_OBJECT(specularity_checkbutton), "shininess_entry", shininess_entry);

      g_object_set_data(G_OBJECT(strength_entry),  "specularity_checkbutton", specularity_checkbutton);
      g_object_set_data(G_OBJECT(shininess_entry), "specularity_checkbutton", specularity_checkbutton);

      g_object_set_data(G_OBJECT(strength_entry),  "imol", GINT_TO_POINTER(imol));
      g_object_set_data(G_OBJECT(shininess_entry), "imol", GINT_TO_POINTER(imol));

      gtk_editable_set_text(GTK_EDITABLE(strength_entry),  coot::util::float_to_string_using_dec_pl(specular_strength, 1).c_str());
      gtk_editable_set_text(GTK_EDITABLE(shininess_entry), coot::util::float_to_string_using_dec_pl(shininess, 1).c_str());

      std::cout << "debug:: fill_single_map_properties_dialog_gtk3() imol " << imol
                << " m.material_for_maps.do_specularity " << m.material_for_maps.do_specularity << std::endl;

      if (m.material_for_maps.do_specularity) {
         gtk_check_button_set_active(GTK_CHECK_BUTTON(specularity_checkbutton), TRUE);
      }
   }

   // -------------------------------------------------------------------------------
   // Fresnel
   // -------------------------------------------------------------------------------

   GtkWidget *fresnel_checkbutton = widget_from_builder("map_properties_dialog_fresnel_state_checkbutton");
   if (fresnel_checkbutton) {
      GtkWidget *bias_entry  = widget_from_builder("map_properties_dialog_fresnel_bias_entry");
      GtkWidget *scale_entry = widget_from_builder("map_properties_dialog_fresnel_scale_entry");
      GtkWidget *power_entry = widget_from_builder("map_properties_dialog_fresnel_power_entry");

      float bias  = m.fresnel_settings.bias;
      float scale = m.fresnel_settings.scale;
      float power = m.fresnel_settings.power;

      g_object_set_data(G_OBJECT(fresnel_checkbutton), "imol", GINT_TO_POINTER(imol));
      g_object_set_data(G_OBJECT(fresnel_checkbutton), "bias_entry",   bias_entry);
      g_object_set_data(G_OBJECT(fresnel_checkbutton), "scale_entry", scale_entry);
      g_object_set_data(G_OBJECT(fresnel_checkbutton), "power_entry", power_entry);

      g_object_set_data(G_OBJECT(bias_entry),  "fresnel_checkbutton", fresnel_checkbutton);
      g_object_set_data(G_OBJECT(scale_entry), "fresnel_checkbutton", fresnel_checkbutton);
      g_object_set_data(G_OBJECT(power_entry), "fresnel_checkbutton", fresnel_checkbutton);

      g_object_set_data(G_OBJECT(bias_entry),  "imol", GINT_TO_POINTER(imol));
      g_object_set_data(G_OBJECT(scale_entry), "imol", GINT_TO_POINTER(imol));
      g_object_set_data(G_OBJECT(power_entry), "imol", GINT_TO_POINTER(imol));

      gtk_editable_set_text(GTK_EDITABLE(bias_entry),  coot::util::float_to_string_using_dec_pl(bias,  2).c_str());
      gtk_editable_set_text(GTK_EDITABLE(scale_entry), coot::util::float_to_string_using_dec_pl(scale, 1).c_str());
      gtk_editable_set_text(GTK_EDITABLE(power_entry), coot::util::float_to_string_using_dec_pl(power, 1).c_str());

      if (m.fresnel_settings.state)
         gtk_check_button_set_active(GTK_CHECK_BUTTON(fresnel_checkbutton), TRUE);
   }

}

GtkWidget *wrapped_create_single_map_properties_dialog_gtk3(int imol) {

   auto w_and_b = create_single_map_properties_dialog_gtk3();
   if (!w_and_b.first) return 0;
   fill_single_map_properties_dialog_gtk3(w_and_b, imol);
   return w_and_b.first;
}


extern "C" G_MODULE_EXPORT
void
on_single_map_properties_dialog_close(GtkDialog *dialog,
                                      gpointer   user_data) {

   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
}

extern "C" G_MODULE_EXPORT
void
on_single_map_properties_dialog_response(GtkDialog *dialog,
                                         gint       response_id,
                                         gpointer   user_data) {

   gtk_widget_set_visible(GTK_WIDGET(dialog), FALSE);
}

#include "gtk-manual.h"
#include "c-interface.h"
#include "c-interface-gtk-widgets.h"


// where is the widget for this?
extern "C" G_MODULE_EXPORT
void
on_single_map_properties_ok_button_clicked(GtkButton       *button,
                                           gpointer         user_data) {

   GtkWidget *window = widget_from_builder("single_map_properties_dialog");

   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(window), "imol"));

  if (is_valid_map_molecule(imol)) {
    set_contour_by_sigma_step_maybe(window, imol);
    skeletonize_map_single_map_maybe(window, imol);
  }
  gtk_widget_set_visible(window, FALSE);

}

/* Not sure that this exists any more... */
extern "C" G_MODULE_EXPORT
void
on_single_map_properties_cancel_button_clicked(GtkButton       *button,
                                               gpointer         user_data) {

   // gtk_widget_set_visible(window, FALSE);

}

extern "C" G_MODULE_EXPORT
void
on_single_map_properties_colour_button_color_set(GtkColorButton *colorbutton,
                                                 gpointer        user_data) {

   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(colorbutton), "imol"));
   if (is_valid_map_molecule(imol)) {
      GdkRGBA rgba;
      gtk_color_chooser_get_rgba(GTK_COLOR_CHOOSER(colorbutton), &rgba);
      graphics_info_t::molecules[imol].set_map_colour(rgba);
      graphics_draw();
   }
}




void on_colour_chooser_dialog_response(GtkDialog *dialog,
                                       int response) {

   if (response == GTK_RESPONSE_OK) {
      GdkRGBA color;
      gtk_color_chooser_get_rgba(GTK_COLOR_CHOOSER(dialog), &color);
      int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(dialog), "imol"));
      if (is_valid_map_molecule(imol)) {
         graphics_info_t g;
         g.molecules[imol].set_map_colour(color);
         // std::cout << "set map colour to " << color.red << " " << color.green << " " << color.blue << std::endl;
         graphics_draw();
      }
   }
   gtk_window_destroy(GTK_WINDOW(dialog));
}

void show_map_colour_selector_with_parent(int imol, GtkWidget *parent_window) {

   if (is_valid_map_molecule(imol)) {

      std::string label = std::string("Coot: Map ") + std::to_string(imol) + std::string(" Colour Selection");
      std::cout << "label: " << label << std::endl;

      // GtkWidget *dialog = gtk_color_dialog_new(); 20230429-PE in 4.10, but not in 4.4 - sadge
      // gtk_color_chooser_dialog_new is/will be deprecated.
      GtkWidget *colour_chooser_dialog = gtk_color_chooser_dialog_new("Test", GTK_WINDOW(parent_window));

      g_object_set_data(G_OBJECT(colour_chooser_dialog), "imol", GINT_TO_POINTER(imol));

      // I don't see this having and effect
      GdkRGBA map_colour = get_map_colour(imol);
      GdkRGBA *map_colour_p = new GdkRGBA(map_colour);
      gtk_color_chooser_set_rgba(GTK_COLOR_CHOOSER(colour_chooser_dialog), map_colour_p);

      gtk_widget_set_visible(colour_chooser_dialog, TRUE);
      GCallback callback = G_CALLBACK(on_colour_chooser_dialog_response);
      g_signal_connect(G_OBJECT(colour_chooser_dialog), "response", callback, GINT_TO_POINTER(imol));
   }
}



void handle_map_properties_specularity_change(int imol, GtkWidget *checkbutton) {

   if (! is_valid_map_molecule(imol)) return;

   molecule_class_info_t &m = graphics_info_t::molecules[imol];

   if (gtk_check_button_get_active(GTK_CHECK_BUTTON(checkbutton))) {
      std::cout << "Turn on specularity " << std::endl;
      GtkWidget *strength_entry  = GTK_WIDGET(g_object_get_data(G_OBJECT(checkbutton),  "strength_entry"));
      GtkWidget *shininess_entry = GTK_WIDGET(g_object_get_data(G_OBJECT(checkbutton), "shininess_entry"));
      std::string strength_entry_text  = gtk_editable_get_text(GTK_EDITABLE(strength_entry));
      std::string shininess_entry_text = gtk_editable_get_text(GTK_EDITABLE(shininess_entry));
      float f1 = coot::util::string_to_float(strength_entry_text);
      float f2 = coot::util::string_to_float(shininess_entry_text);
      m.material_for_maps.specular_strength = f1;
      m.material_for_maps.shininess         = f2;
      m.material_for_maps.turn_specularity_on(true);
      std::cout << "in handle_map_properties_specularity_change() imol: " << imol << " do: " <<  m.material_for_maps.do_specularity
                << " strength " << m.material_for_maps.specular_strength << " shiny " << m.material_for_maps.shininess << std::endl;
   } else {
      std::cout << "Turn off specularity " << std::endl;
      m.material_for_maps.turn_specularity_on(false);
   }
   graphics_draw();
}

