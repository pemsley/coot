
#include <string>
#include <iostream>

#include "utils/coot-utils.hh"
#include "single-map-properties-dialog.hh"
#include "graphics-info.h"

std::pair<GtkWidget *, GtkBuilder *> create_single_map_properties_dialog_gtk3() {

   auto expand_file = [] (const std::string &file_name) {
                         // FIXME
                         std::string ff = file_name;
                         return ff;
                   };

   GtkWidget *single_map_properties_dialog = 0;
   GtkBuilder *builder = gtk_builder_new();
   std::string glade_file_name = "single-map-properties-dialog.glade";
   std::string single_map_properties_dialog_glade_file_name = expand_file(glade_file_name);
   guint add_from_file_status = gtk_builder_add_from_file(builder, single_map_properties_dialog_glade_file_name.c_str(), NULL);

   if (add_from_file_status) {
      std::cout << "::::::::::::::::::::::::::: single_map_properties_dialog from gtkbuilder" << std::endl;
      single_map_properties_dialog = GTK_WIDGET(gtk_builder_get_object(builder, "single_map_properties_dialog"));
      gtk_builder_connect_signals(builder, single_map_properties_dialog); // if "nothing happens" then you've missed this call
   } else {
      std::cout << "ERROR:: failed to get builder file for single-map-properties dialog" << std::endl;
   }

   return std::make_pair(single_map_properties_dialog, builder);

}

void fill_single_map_properties_dialog_gtk3(std::pair<GtkWidget *, GtkBuilder *> w_and_b, int imol) {

   std::cout << "::::::::::::::::::::::::::: fill_single_map_properties_dialog_gtk3() for imol " << imol << std::endl;

   graphics_info_t g;
   if (! g.is_valid_map_molecule(imol)) return;

   GtkWidget *dialog   = w_and_b.first;
   GtkBuilder *builder = w_and_b.second;
   g_object_set_data(G_OBJECT(dialog), "imol", GINT_TO_POINTER(imol));

   auto widget_from_builder = [builder] (const std::string &wid) {
                                 return GTK_WIDGET(gtk_builder_get_object(builder, wid.c_str()));
                              };

   auto set_contour_sigma_button_and_entry = [widget_from_builder, imol] (float contour_rmsd_step,
                                                                          bool contour_by_rmsd_flag) {

                                                GtkWidget *cl_apply_button = widget_from_builder("single_map_properties_contour_level_apply_button");
                                                GtkWidget *cl_entry        = widget_from_builder("single_map_properties_contour_level_entry");
                                                GtkWidget *level_type_radiobutton =
                                                   widget_from_builder("single_map_properties_absolute_radiobutton");
                                                if (contour_by_rmsd_flag)
                                                   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(level_type_radiobutton), FALSE);
                                                g_object_set_data(G_OBJECT(level_type_radiobutton), "contour_level_entry", cl_entry);
                                                std::string entry_text = coot::util::float_to_string_using_dec_pl(contour_rmsd_step, 2);
                                                gtk_entry_set_text(GTK_ENTRY(cl_entry), entry_text.c_str());
                                                g_object_set_data(G_OBJECT(cl_apply_button), "imol", GINT_TO_POINTER(imol));
                                                g_object_set_data(G_OBJECT(cl_apply_button), "contour_level_entry", cl_entry);
                                                g_object_set_data(G_OBJECT(cl_apply_button), "single_map_properties_absolute_radiobutton", level_type_radiobutton);

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

   std::string title = "Properties for Map " + std::to_string(imol);
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

   // And now the map rendering style: transparent surface or standard lines:
   GtkWidget *rb_1  = widget_from_builder("display_map_style_as_lines_radiobutton");
   GtkWidget *rb_2  = widget_from_builder("display_map_style_surface_radiobutton");
   GtkWidget *scale = widget_from_builder("map_opacity_hscale");

   g_object_set_data(G_OBJECT(rb_1), "imol", GINT_TO_POINTER(imol));

   GtkWidget *level_type_radiobutton = widget_from_builder("single_map_properties_absolute_radiobutton");
   g_object_set_data(G_OBJECT(level_type_radiobutton), "imol", GINT_TO_POINTER(imol));

   const molecule_class_info_t &m = g.molecules[imol];

   if (! m.draw_it_for_map_standard_lines) {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(rb_2), TRUE);
   } else {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(rb_1), TRUE);
   }

   g_object_set_data(G_OBJECT(scale), "imol", GINT_TO_POINTER(imol));
   GtkAdjustment *adjustment = gtk_range_get_adjustment(GTK_RANGE(scale));
   float op = m.density_surface_opacity;
   gtk_adjustment_set_value(adjustment, 100.0*op);

   GtkWidget *map_contour_frame = widget_from_builder("single_map_properties_map_histogram_frame");
   GtkWidget *alignment = widget_from_builder("alignment_for_map_density_histogram");
   // fill_map_histogram_widget(imol, alignment);
   
   set_contour_sigma_button_and_entry(m.contour_sigma_step, m.contour_by_sigma_flag);


   /*  and now the skeleton buttons */
   // GtkWidget *frame = lookup_widget(dialog, "single_map_skeleton_frame");
   // set_on_off_single_map_skeleton_radio_buttons(frame, imol);
   /* contour by sigma step */

   // -------------------------------------------------------------------------------
   // colour
   // -------------------------------------------------------------------------------

   GtkWidget *colour_button = widget_from_builder("single_map_properties_colour_button");
   if (colour_button) {
      g_object_set_data(G_OBJECT(colour_button), "imol",   GINT_TO_POINTER(imol));
      g_object_set_data(G_OBJECT(colour_button), "parent", dialog);
      
   }

   // -------------------------------------------------------------------------------
   // Specularity
   // -------------------------------------------------------------------------------

   GtkWidget *specularity_checkbutton = widget_from_builder("map_properties_dialog_specularity_state_checkbutton");
   if (specularity_checkbutton) {
      
      GtkWidget *strength_entry  = widget_from_builder("map_properties_dialog_specularity_strength_entry");
      GtkWidget *shininess_entry = widget_from_builder("map_properties_dialog_specularity_shininess_entry");

      molecule_class_info_t &m = graphics_info_t::molecules[imol];
      float specular_strength = m.material_for_maps.specular_strength;
      float shininess = m.material_for_maps.shininess;

      g_object_set_data(G_OBJECT(specularity_checkbutton), "imol", GINT_TO_POINTER(imol));
      g_object_set_data(G_OBJECT(specularity_checkbutton), "strength_entry",   strength_entry);
      g_object_set_data(G_OBJECT(specularity_checkbutton), "shininess_entry", shininess_entry);

      g_object_set_data(G_OBJECT(strength_entry),  "specularity_checkbutton", specularity_checkbutton);
      g_object_set_data(G_OBJECT(shininess_entry), "specularity_checkbutton", specularity_checkbutton);

      g_object_set_data(G_OBJECT(strength_entry),  "imol", GINT_TO_POINTER(imol));
      g_object_set_data(G_OBJECT(shininess_entry), "imol", GINT_TO_POINTER(imol));

      gtk_entry_set_text(GTK_ENTRY(strength_entry),  coot::util::float_to_string_using_dec_pl(specular_strength, 1).c_str());
      gtk_entry_set_text(GTK_ENTRY(shininess_entry), coot::util::float_to_string_using_dec_pl(shininess, 1).c_str());

      std::cout << "debug:: fill_single_map_properties_dialog_gtk3() imol " << imol
                << " m.material_for_maps.do_specularity " << m.material_for_maps.do_specularity << std::endl;

      if (m.material_for_maps.do_specularity) {
         gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(specularity_checkbutton), TRUE);
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

      molecule_class_info_t &m = graphics_info_t::molecules[imol];
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
      g_object_set_data(G_OBJECT(scale_entry), "imol", GINT_TO_POINTER(scale));
      g_object_set_data(G_OBJECT(power_entry), "imol", GINT_TO_POINTER(power));

      gtk_entry_set_text(GTK_ENTRY(bias_entry),  coot::util::float_to_string_using_dec_pl(bias,  2).c_str());
      gtk_entry_set_text(GTK_ENTRY(scale_entry), coot::util::float_to_string_using_dec_pl(scale, 1).c_str());
      gtk_entry_set_text(GTK_ENTRY(power_entry), coot::util::float_to_string_using_dec_pl(power, 1).c_str());

      if (m.fresnel_settings.state)
         gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(fresnel_checkbutton), TRUE);
   }

   printf("!!!!!!!!! fill_single_map_properties_dialog_gtk3() done\n");
}

GtkWidget *wrapped_create_single_map_properties_dialog_gtk3(int imol) {

   auto w_and_b = create_single_map_properties_dialog_gtk3();
   fill_single_map_properties_dialog_gtk3(w_and_b, imol);
   return w_and_b.first;
}


extern "C" G_MODULE_EXPORT
void
on_single_map_properties_dialog_close_gtkbuilder_callback (GtkDialog *dialog,
                                                           gpointer   user_data) {

}

extern "C" G_MODULE_EXPORT
void
on_single_map_properties_dialog_response_gtkbuilder_callback (GtkDialog *dialog,
                                                              gint       response_id,
                                                              gpointer   user_data) {
   gtk_widget_destroy(GTK_WIDGET(dialog));
}

