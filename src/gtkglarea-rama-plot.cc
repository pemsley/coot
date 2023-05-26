
// this is the Gtk Wrapper of the OpengGL Rama plot

#include "graphics-info.h"

#include "gtkglarea-rama-plot.hh"
#include "gl-rama-plot.hh"
#include "widget-from-builder.hh"

void
gtkgl_rama_realize(GtkWidget *gtk_gl_area) {

   if (!gtk_gl_area) return;

   // ----- expand the Rama pane if needed

   graphics_info_t g;
   if (!g.rama_plot_boxes.empty()) {
      GtkWidget *paned = widget_from_builder("main_window_graphics_rama_vs_graphics_pane");
      int position = gtk_paned_get_position(GTK_PANED(paned));
      std::cout << ":::::::::::: in gtkgl_rama_realize() the paned had position " << position << std::endl;
      if (position < 10) {
         gtk_paned_set_position(GTK_PANED(paned), 400);
         std::cout << ":::::::::::: gtk_paned_set_position 400 here " << std::endl;
      }
   }

   // ----- find the right rama plot box for this gtk_gl_area and set up the rama
   //       in it with the molecule from the imol in the box.

   bool done = false;
   for (unsigned int i=0; i<g.rama_plot_boxes.size(); i++) {
      auto &rama_box = g.rama_plot_boxes[i];
      if (rama_box.gtk_gl_area == gtk_gl_area) {

         // Dangerous for main window code
         // GdkGLContext *context = gtk_gl_area_get_context(GTK_GL_AREA(gtk_gl_area)); // needed?
         gtk_gl_area_make_current(GTK_GL_AREA (gtk_gl_area));

         rama_box.rama.setup_buffers(0.9);
         int imol = g.rama_plot_boxes[i].imol;
         const std::string residue_selection = rama_box.residue_selection;
         auto &m = graphics_info_t::molecules[imol];
         g.rama_plot_boxes[i].rama.setup_from(imol, m.atom_sel.mol, residue_selection);
         done = true;
      }
   }

   if (!done) {
      std::cout << "WARNING:: oops - failed to setup in gtkgl_rama_realize() " << gtk_gl_area
                << " with " << g.rama_plot_boxes.size() << " rama-boxs " << std::endl;
   }

}

void
gtkgl_rama_unrealize(GtkWidget *gl_area) {


}

void
gtkgl_rama_on_glarea_render(GtkWidget *gtk_gl_area) {

   graphics_info_t g;
   for (unsigned int i=0; i<g.rama_plot_boxes.size(); i++) {
      if (g.rama_plot_boxes[i].gtk_gl_area == gtk_gl_area) {

         // Dangerous for main window code
         // GdkGLContext *context = gtk_gl_area_get_context(GTK_GL_AREA(gtk_gl_area)); // needed?
         gtk_gl_area_make_current(GTK_GL_AREA (gtk_gl_area));

         GtkAllocation allocation;
         gtk_widget_get_allocation(GTK_WIDGET(gtk_gl_area), &allocation);
         int w = allocation.width;
         int h = allocation.height;

         g.rama_plot_boxes[i].rama.draw(&g.shader_for_rama_plot_axes_and_ticks,
                                        &g.shader_for_rama_plot_phi_phis_markers, // instanced
                                        &g.shader_for_hud_image_texture,
                                        w, h, w, h); // background texture (not text!), uses window_resize_position_correc
      }
   }
}

void
gtkgl_rama_on_glarea_resize(GtkWidget *gl_area, gint width, gint height) {

   std::cout << "resize gl rama to " << width << " " << height << std::endl;
}


void show_opengl_ramachandran_plot(int imol, const std::string &residue_selection) {

   // find a better name for this function?

   auto on_rama_glarea_click = +[] (GtkGestureClick* click_gesture,
                                   gint n_press,
                                   gdouble x,
                                   gdouble y,
                                   gpointer user_data) {

      GtkWidget *gl_area = GTK_WIDGET(user_data);
      
      GtkAllocation allocation;
      gtk_widget_get_allocation(gl_area, &allocation);
      int w = allocation.width;
      int h = allocation.height;

      // std::cout << "Rama click! " << x << " " << y << " width " << w << " height " << h << std::endl;

      // Find the right rama plot and get the mouse-over hit, and if it was a residue, go there.
      graphics_info_t g;
      for (unsigned int i=0; i<g.rama_plot_boxes.size(); i++) {
         const auto &rama_box = g.rama_plot_boxes[i];
         if (rama_box.matches_gl_area(gl_area)) {
            auto rama_plot_hit = rama_box.rama.get_mouse_over_hit(x, y, w, h);
            if (rama_plot_hit.residue_was_clicked) {
               int imol = rama_box.imol;
               g.go_to_residue(imol, rama_plot_hit.residue_spec);
            }
         }
      }
   };

   graphics_info_t g;
   if (g.is_valid_model_molecule(imol)) {

      GtkWidget *box_for_all_plots = widget_from_builder("ramachandran_plots_vbox");

      GtkWidget *box_for_this_plot = gtk_box_new(GTK_ORIENTATION_VERTICAL, 2);
      GtkWidget *gl_area = gtk_gl_area_new();
      GtkWidget *close_button = gtk_button_new_with_label("Close");
      GtkWidget *box_for_selection = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 2);
      GtkWidget *selection_label = gtk_label_new("Selection: ");
      GtkWidget *selection_entry = gtk_entry_new();
      gtk_editable_set_text(GTK_EDITABLE(selection_entry), residue_selection.c_str());
      gtk_widget_set_margin_start(selection_label, 6);
      gtk_widget_set_margin_start(box_for_this_plot, 6);
      gtk_widget_set_margin_start(close_button, 6);
      gtk_widget_set_margin_end(close_button, 6);

      gl_rama_plot_t rama;
      graphics_info_t::widgeted_rama_plot_t wr(imol, residue_selection, rama, gl_area, close_button, box_for_this_plot);
      g.rama_plot_boxes.push_back(wr);

      gtk_widget_set_size_request(gl_area, 400, 400);
      g_signal_connect(gl_area, "realize",   G_CALLBACK(gtkgl_rama_realize),   NULL);
      g_signal_connect(gl_area, "unrealize", G_CALLBACK(gtkgl_rama_unrealize), NULL);
      g_signal_connect(gl_area, "render",    G_CALLBACK(gtkgl_rama_on_glarea_render),  NULL);
      g_signal_connect(gl_area, "resize",    G_CALLBACK(gtkgl_rama_on_glarea_resize),  NULL);

      gtk_widget_set_can_focus(gl_area, TRUE);
      gtk_widget_set_focusable(gl_area, TRUE);

      gtk_widget_set_hexpand(gl_area, FALSE);
      gtk_widget_set_vexpand(gl_area, FALSE);

      GtkGesture *click_controller          = gtk_gesture_click_new();
      gtk_widget_add_controller(GTK_WIDGET(gl_area), GTK_EVENT_CONTROLLER(click_controller));
      g_signal_connect(click_controller, "pressed",  G_CALLBACK(on_rama_glarea_click), gl_area);

      auto selection_entry_activate_callback = +[] (GtkWidget *entry, gpointer user_data) {
         std::string entry_string = gtk_editable_get_text(GTK_EDITABLE(entry));
         std::cout << "Now do something with " << entry_string << std::endl;
      };

      g_signal_connect(G_OBJECT(selection_entry), "activate", G_CALLBACK(selection_entry_activate_callback), selection_entry);

      auto close_callback = +[] (GtkWidget *close_button, gpointer user_data) {
         GtkWidget *box_for_all_plots = widget_from_builder("ramachandran_plots_vbox");
         GtkWidget *box_for_this_plot = GTK_WIDGET(user_data);
         graphics_info_t g;
         g.remove_plot_from_rama_plots(box_for_this_plot); // hides the rama pane in main_window_graphics_rama_vs_graphics_pane 
                                                           // if there are no plots left.
         gtk_box_remove(GTK_BOX(box_for_all_plots), box_for_this_plot);
      };

      g_signal_connect(G_OBJECT(close_button), "clicked", G_CALLBACK(close_callback), box_for_this_plot);

      gtk_box_append(GTK_BOX(box_for_selection), selection_label);
      gtk_box_append(GTK_BOX(box_for_selection), selection_entry);

      gtk_box_append(GTK_BOX(box_for_this_plot), gl_area);
      gtk_box_append(GTK_BOX(box_for_this_plot), box_for_selection);
      gtk_box_append(GTK_BOX(box_for_this_plot), close_button);
      gtk_box_append(GTK_BOX(box_for_all_plots), box_for_this_plot);

      gtk_widget_show(gl_area);

   }
}

// static
void
graphics_info_t::remove_plot_from_rama_plots(GtkWidget *plot_box) {

   std::vector<widgeted_rama_plot_t>::const_iterator it;
   for (it=rama_plot_boxes.begin(); it!=rama_plot_boxes.end(); ++it) {
      GtkWidget *this_plot_box = it->box;
      if (this_plot_box == plot_box) {
         rama_plot_boxes.erase(it);
         break;
      }
   }

   if (rama_plot_boxes.empty()) {
      GtkWidget *scrolled = widget_from_builder("ramachandran_plots_scrolled_window");
      if (GTK_IS_SCROLLED_WINDOW(scrolled))
         gtk_widget_set_visible(scrolled, FALSE);
      else
         std::cout << "Not a scrolled window " << scrolled << std::endl;
   }
}



//
//  static
void
graphics_info_t::draw_rama_plots() {

   for (unsigned int i=0; i<rama_plot_boxes.size(); i++) {
      GtkGLArea *gl_area = GTK_GL_AREA(rama_plot_boxes[i].gtk_gl_area);
      if (GTK_IS_GL_AREA(gl_area)) {
         GdkGLContext *context = gtk_gl_area_get_context(gl_area); // needed?
         gtk_gl_area_make_current(GTK_GL_AREA (gl_area));

         GtkAllocation allocation;
         gtk_widget_get_allocation(GTK_WIDGET(gl_area), &allocation);
         int w = allocation.width;
         int h = allocation.height;
         rama_plot_boxes[i].rama.draw(&shader_for_rama_plot_axes_and_ticks,
                                      &shader_for_rama_plot_phi_phis_markers, // instanced
                                      &shader_for_hud_image_texture,
                                      w, h, w, h);
      } else {
         std::cout << "ERROR:: ploting rama plot " << i << " which hash gl_area that has gone out of scope!"
                   << std::endl;
      }
   }
}
