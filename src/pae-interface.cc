
#include "graphics-info.h"
#include "c-interface.h"

#include "MoleculesToTriangles/CXXClasses/MyMolecule.h"
#include "cc-interface-molecular-representation.hh"

void
grey_ribbons_for_active_molecule() {

   auto add_colour = [] (const coot::colour_holder &col,
			 unsigned int idx,
			 std::vector<coot::colour_holder> *colours_p) {

      if (idx < colours_p->size()) {
	 colours_p->at(idx) = col;
      } else {
	 unsigned int cc = colours_p->capacity();
	 if (cc < (idx+1))
	    colours_p->reserve(2 * idx);
	 colours_p->resize(idx+1);
	 colours_p->at(idx) = col;
      }
   };

   graphics_info_t g;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = g.active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      std::string atom_selection = "/";
      std::string ColorScheme = "Chain";
      std::string style = "Ribbon";

      std::vector<std::pair<std::string, unsigned int> > cis;
      std::vector<coot::colour_holder> colours;
      coot::colour_holder base(0.7, 0.7, 0.7);
      cis.push_back(std::make_pair("//", 60));
      add_colour(base, 60, &colours);
      g.set_user_defined_colours(colours);
      g.molecules[imol].set_user_defined_colour_indices_by_selections(cis);

      add_ribbon_representation_with_user_defined_colours(imol, "AlphaFold-PAE");
      g.graphics_draw();
      
      set_mol_displayed(imol, 0);
   }

}


#include "coot-utils/pae.hh"

void
display_pae_from_file_in_a_dialog(const std::string &file_name) {

   auto draw_test = +[] (GtkDrawingArea *da, cairo_t *cr, int w, int h, gpointer data) {
      // Set the color to red
      cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);

      // Calculate the center of the window
      double center_x = w / 2.0;
      double center_y = h / 2.0;

      // Draw the circle
      double radius = (w < h ? w : h) / 4.0;
      cairo_arc(cr, center_x, center_y, radius, 0, 2 * G_PI);
      cairo_fill(cr);
   };

   auto draw_pae = +[] (GtkDrawingArea *da, cairo_t *cr, int w, int h, gpointer data) {

      // Set the background to white
      cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
      cairo_rectangle(cr, 0, 0, w, h);
      cairo_fill(cr);

      cairo_set_source_rgb(cr, 1.0, 0.0, 0.0);
      double center_x = w / 2.0;
      double center_y = h / 2.0;
      double radius = 5.0;
      cairo_arc(cr, center_x, center_y, radius, 0, 2 * G_PI);
      cairo_fill(cr);

      pae_t *pae_p = static_cast<pae_t *>(data);
      int n_pixels_for_pae_image = pae_p->n_pixels - 100;

      // move the PAE graph in the base image
      int offset_x = 80;
      int offset_y = 20;
      unsigned int n_residues = pae_p->pae_vecs.size();

      // let's make the actual data plot as image data:
      unsigned char *image_data = new unsigned char[n_pixels_for_pae_image*n_pixels_for_pae_image*4];
      float max_value = pae_p->get_max_value(pae_p->pae_vecs); // silly
      for (int i=0; i<n_pixels_for_pae_image; i++) {
         for (int j=0; j<n_pixels_for_pae_image; j++) {
            int idx = 4*(i*n_pixels_for_pae_image+j);
            float f_x = static_cast<float>(i)/static_cast<float>(n_pixels_for_pae_image);
            float f_y = static_cast<float>(j)/static_cast<float>(n_pixels_for_pae_image);
            int i_res_index = static_cast<int>(f_x * static_cast<float>(n_residues));
            int j_res_index = static_cast<int>(f_y * static_cast<float>(n_residues));
            float v = static_cast<float>(pae_p->pae_vecs[i_res_index][j_res_index]);
            auto col = pae_p->value_to_colour(v, max_value);
            image_data[idx]   = col[0];
            image_data[idx+1] = col[1];
            image_data[idx+2] = col[2];
         }
      }

      cairo_surface_t *pae_img_surface =
         cairo_image_surface_create_for_data(image_data, CAIRO_FORMAT_RGB24,
                                             n_pixels_for_pae_image, n_pixels_for_pae_image,
                                             n_pixels_for_pae_image*4);
      cairo_set_source_surface(cr, pae_img_surface, 80, 20);
      cairo_paint(cr);

      // draw a box around the pae plot
      cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
      cairo_rectangle(cr, offset_x, offset_y, n_pixels_for_pae_image, n_pixels_for_pae_image); // x, y, width, height
      cairo_stroke(cr);

      cairo_set_font_size(cr, offset_y);

      // "scored residue" label
      cairo_move_to(cr, 260, 573);
      cairo_show_text(cr, "Scored Residue");

      // "aligned residue" label
      cairo_move_to(cr, 26, 360);
      cairo_save(cr);
      cairo_rotate(cr, - M_PI / 2.0); //
      cairo_show_text(cr, "Aligned Residue");
      cairo_restore(cr);

      // "Expected position error" label
      cairo_move_to(cr, 150, 680);
      cairo_show_text(cr, "Expected position error (Ångströms)");

      // actual tick marks would be good here.

      // tick labels - x axis
      unsigned int tick_res_no = 0;
      while (tick_res_no < n_residues) {
         float f = static_cast<float>(tick_res_no) / static_cast<float>(n_residues);
         float pixel_for_tick_res_no = static_cast<float>(n_pixels_for_pae_image) * f;
         int x = offset_x + static_cast<int>(pixel_for_tick_res_no);
         int x_tweak = -20;
         if (tick_res_no == 0) x_tweak = -5;
         x += x_tweak;
         cairo_move_to(cr, x, 545);
         std::string text = std::to_string(tick_res_no);
         cairo_show_text(cr, text.c_str());
         tick_res_no += 100;
      }
      // tick labels - y axis
      tick_res_no = 0;
      while (tick_res_no < n_residues) {
         float f = static_cast<float>(tick_res_no) / static_cast<float>(n_residues);
         float pixel_for_tick_res_no = static_cast<float>(n_pixels_for_pae_image) * f;
         int y = offset_y + static_cast<int>(pixel_for_tick_res_no);
         y += 8; // so that the middle of the label text is at the tick position (rather than the bottom)
         double text_width = 42.0;
         if (tick_res_no == 0) text_width = 18.0;
         cairo_move_to(cr, offset_x-text_width, y);
         std::string text = std::to_string(tick_res_no);
         cairo_show_text(cr, text.c_str());
         tick_res_no += 100;
      }

      // legend colour ramp
      //
      for(double epe=0.0; epe<max_value; epe += 0.5) {
         auto col = pae_p->value_to_colour(epe, max_value);
         double r = static_cast<double>(col[0])/255.0;
         double g = static_cast<double>(col[1])/255.0;
         double b = static_cast<double>(col[2])/255.0;
         cairo_set_source_rgb(cr, r, g, b);
         double f = epe/max_value;
         double x = offset_x + f * static_cast<double>(n_pixels_for_pae_image);
         double y = offset_y + n_pixels_for_pae_image + 80.0;
         cairo_rectangle(cr, x, y, 9.0, 20.0); // x, y, width, height
         cairo_fill(cr);
      }

      // draw a box around the legend
      cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
      float offset_y_legend = offset_y + n_pixels_for_pae_image + 80.0;
      cairo_rectangle(cr, offset_x, offset_y_legend, n_pixels_for_pae_image, 20.0); // x, y, width, height
      cairo_stroke(cr);

      // tick labels for the legend
      int tick_pos_error = 0;
      while (tick_pos_error < max_value) {
         float f = static_cast<float>(tick_pos_error) / static_cast<float>(max_value);
         float pixel_for_tick = static_cast<float>(n_pixels_for_pae_image) * f; // the legend is the same length (width)
         int x = offset_x + static_cast<int>(pixel_for_tick);
         int x_tweak = -20;
         if (tick_pos_error < 10) x_tweak = -7;
         x += x_tweak;
         cairo_move_to(cr, x, 645);
         std::string text = std::to_string(tick_pos_error);
         cairo_show_text(cr, text.c_str());
         tick_pos_error += 5;
      }
   };


   if (coot::file_exists(file_name)) {

      pae_t pae(file_name, 600);
      GtkWidget *window = gtk_window_new();
      std::string title = "Coot: " + coot::util::file_name_non_directory(file_name);
      gtk_window_set_title(GTK_WINDOW(window), title.c_str());
      GtkWidget *vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
      GtkWidget *drawing_area = gtk_drawing_area_new();
      gtk_widget_set_hexpand(drawing_area, TRUE); // Allow expansion in x direction
      gtk_widget_set_vexpand(drawing_area, TRUE); // Allow expansion in y direction
      gtk_window_set_child(GTK_WINDOW(window), vbox);
      gtk_box_append(GTK_BOX(vbox), drawing_area);

      GtkWidget *buttons_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
      GtkWidget *button_close = gtk_button_new_with_label("Close");
      GtkWidget *ribbons_button = gtk_button_new_with_label("Ribbons");
      gtk_widget_set_halign(buttons_box, GTK_ALIGN_END);
      gtk_box_append(GTK_BOX(buttons_box), ribbons_button);
      gtk_box_append(GTK_BOX(buttons_box), button_close);
      gtk_box_append(GTK_BOX(vbox), buttons_box);

      for (auto b : {button_close, ribbons_button}) {
	 gtk_widget_set_margin_start(b, 10);
	 gtk_widget_set_margin_end(b, 10);
	 gtk_widget_set_margin_top(b, 14);
	 gtk_widget_set_margin_bottom(b, 10);
      }

      auto close_button_callback = +[] (GtkButton *button, gpointer data) {
         GtkWindow *window = static_cast<GtkWindow *>(data);
         gtk_window_destroy(window);
      };
      g_signal_connect(G_OBJECT(button_close), "clicked", G_CALLBACK(close_button_callback), window);
      auto ribbon_button_callback = +[] (GtkButton *button, gpointer data) {
	 std::cout << "clicked!" << std::endl;
	 grey_ribbons_for_active_molecule();
      };
      g_signal_connect(G_OBJECT(ribbons_button), "clicked", G_CALLBACK(ribbon_button_callback), window);

      gtk_window_set_default_size(GTK_WINDOW(window), 700, 780);
      pae_t *pae_p = new pae_t(pae);
      gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(drawing_area), draw_pae, (gpointer) (pae_p), NULL);

      // attach the pae_p (set "pae_p") to the pae_p so that it can be delete on closing the window.

      gtk_widget_set_visible(window, TRUE);
      graphics_info_t g;
      g.set_transient_for_main_window(window);

   }

}
