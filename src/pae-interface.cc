
#include "graphics-info.h"
#include "c-interface.h"

#include "MoleculesToTriangles/CXXClasses/MyMolecule.h"
#include "cc-interface-molecular-representation.hh"

#include "coot-utils/pae.hh"

// static bool pae_dragging = false;
// static double pae_drag_start_x = 0, pae_drag_start_y = 0;
// static double pae_drag_end_x = 0, pae_drag_end_y = 0;

// typedef struct {
//     bool dragging;
//     double start_x, start_y;
//     double end_x, end_y;
// } SelectionData;


class augmented_pae_info_t {
public:

   class box_info_t {
   public:
      double start_x;
      double start_y;
      double end_x;
      double end_y;
      glm::vec4 col;
      box_info_t(double sx, double sy, double ex, double ey, const glm::vec4 &c) :
         start_x(sx), start_y(sy), end_x(ex), end_y(ey), col(c) {}
   };
   int imol;
   pae_t pae_info;
   std::vector<box_info_t> boxes;
   bool dragging_mode;
   cairo_t *cr;
   augmented_pae_info_t(int imol, const pae_t &p) : imol(imol), pae_info(p), dragging_mode(false), cr(nullptr) {}

   // fill colours_p as needed
   // 20260122-PE I don't think this function is now needed.
   static void add_colour(const coot::colour_holder &col,
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

};

void
grey_ribbons_for_active_molecule() {

   graphics_info_t g;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = g.active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      std::string atom_selection = "/";
      std::string ColorScheme = "Chain";
      std::string style = "Ribbon";

      std::vector<std::pair<std::string, unsigned int> > cis;
      std::vector<std::pair<unsigned int, coot::colour_holder> > colours;
      coot::colour_holder base(0.7, 0.7, 0.7);
      cis.push_back(std::make_pair("//", 60));
      colours.push_back(std::make_pair(60, base));
      g.set_user_defined_colours(colours);
      g.molecules[imol].set_user_defined_colour_indices_by_selections(cis);

      // add_ribbon_representation_with_user_defined_colours(imol, "AlphaFold-PAE");
      g.graphics_draw();

      set_mol_displayed(imol, 0);
   }
}


void
draw_drag_box(augmented_pae_info_t *data) {

   const augmented_pae_info_t api = *data;
   for (unsigned int ibox=0; ibox<api.boxes.size(); ibox++) {

      const augmented_pae_info_t::box_info_t &box = api.boxes[ibox];
      double x = std::min(box.start_x, box.end_x);
      double y = std::min(box.start_y, box.end_y);
      double w = fabs(box.end_x - box.start_x);
      double h = fabs(box.end_y - box.start_y);

      cairo_t *cr = data->cr;
      // cairo_set_source_rgba(cr, 0.2, 0.4, 0.8, 0.3);
      cairo_set_source_rgba(cr, box.col.r, box.col.g, box.col.b, box.col.a);
      cairo_rectangle(cr, x, y, w, h);
      cairo_fill(cr);
      cairo_set_source_rgba(cr, box.col.r, box.col.g, box.col.b, 0.7);
      cairo_rectangle(cr, x, y, w, h);
      cairo_stroke(cr);
   }

}

static void
on_drag_begin(GtkGestureDrag *gesture, double start_x, double start_y, gpointer user_data) {

   // std::cout << "drag begin" << std::endl;

   augmented_pae_info_t *data = static_cast<augmented_pae_info_t*>(user_data);
   data->dragging_mode = true;
   unsigned int n_boxes = data->boxes.size();
   float f =static_cast<float>(n_boxes) / 10.0f;
   float r = 0.2 + 2.0 * f;
   float g = 0.4 + f * 0.5;
   float b = 0.8 - 2.0 * f;
   augmented_pae_info_t::box_info_t box(start_x, start_y, start_x, start_y, glm::vec4(r,g,b,0.3));
   data->boxes.push_back(box);
   GtkWidget *widget = gtk_event_controller_get_widget(GTK_EVENT_CONTROLLER(gesture));
   gtk_widget_queue_draw(widget);
}

void
pae_colour_ribbon_model(int imol,
                        unsigned int box_index,
                        const augmented_pae_info_t::box_info_t &box,
                        const std::vector<std::vector<int> > &pae_vecs) {

   std::string chain_id = "A"; // this is a hack, we should get the chain id from the pae_vecs
   unsigned int n_residues = pae_vecs.size();

   std::cout << "pae_colour_ribbon_model" << std::endl;
   std::cout << "start_x: " << box.start_x << " start_y: " << box.start_y
             << " end_x: " << box.end_x << " end_y: " << box.end_y << std::endl;

   graphics_info_t g;
   float f1 = static_cast<float>(box.start_x) / static_cast<float>(n_residues);
   float f2 = static_cast<float>(box.end_x)   / static_cast<float>(n_residues);
   int start_residue_index = static_cast<int>(f1 * static_cast<float>(n_residues));
   int   end_residue_index = static_cast<int>(f2 * static_cast<float>(n_residues));

   int max_index = n_residues - 1;
   if (start_residue_index < 0) start_residue_index = 0;
   if (  end_residue_index < 0)   end_residue_index = 0;
   if (start_residue_index > max_index) start_residue_index = max_index;
   if (  end_residue_index > max_index)   end_residue_index = max_index;

   int start_residue_number = pae_vecs[0][start_residue_index];
   int end_residue_number = pae_vecs[0][end_residue_index];

   std::string selection = "//" + chain_id + "/" +
      std::to_string(start_residue_number) + "-" + std::to_string(end_residue_number);
   coot::colour_holder ch(box.col.r, box.col.g, box.col.b, 1.0);

   std::vector<std::pair<unsigned int, coot::colour_holder> > colours = g.user_defined_colours;
   unsigned int colour_index = 60 + box_index + 1;
   // augmented_pae_info_t::add_colour(ch, colour_index, &colours);
   colours.push_back(std::make_pair(colour_index, ch));
   g.set_user_defined_colours(colours);
   std::vector<std::pair<std::string, unsigned int> > cis;
   std::cout << "selection: " << selection << " colour: " << ch << std::endl;
   std::pair<std::string, unsigned int> cc(selection, colour_index);
   cis.push_back(cc);
   g.molecules[imol].set_user_defined_colour_indices_by_selections(cis);
   g.undisplay_all_molecule_meshes(imol);
   add_ribbon_representation_with_user_defined_colours(imol, "AlphaFold-PAE");

   GtkWidget *w = widget_from_builder("molecular_representations_dialog");
   gtk_widget_set_visible(w, TRUE);
   g.set_transient_for_main_window(w);
   g.update_molecular_representation_widgets();

}

static void
on_drag_update(GtkGestureDrag *gesture, double offset_x, double offset_y, gpointer user_data) {

   // std::cout << "drag update" << std::endl;
   augmented_pae_info_t *data = static_cast<augmented_pae_info_t*>(user_data);

   // Drag offset is relative to drag begin
   double start_x, start_y;
   gtk_gesture_drag_get_start_point(gesture, &start_x, &start_y);
   augmented_pae_info_t::box_info_t &box = data->boxes.back();
   // data->end_x = start_x + offset_x;
   // data->end_y = start_y + offset_y;
   box.end_x = start_x + offset_x;
   box.end_y = start_y + offset_y;
   GtkWidget *widget = gtk_event_controller_get_widget(GTK_EVENT_CONTROLLER(gesture));
   gtk_widget_queue_draw(widget);
}

static void
on_drag_end(GtkGestureDrag *gesture, double offset_x, double offset_y, gpointer user_data) {

   // std::cout << "drag end" << std::endl;
   augmented_pae_info_t *data = static_cast<augmented_pae_info_t*>(user_data);
   double start_x, start_y;
   gtk_gesture_drag_get_start_point(gesture, &start_x, &start_y);
   // data->end_x = start_x + offset_x;
   // data->end_y = start_y + offset_y;
   data->dragging_mode = false;
   GtkWidget *widget = gtk_event_controller_get_widget(GTK_EVENT_CONTROLLER(gesture));
   gtk_widget_queue_draw(widget);
   // You can now use data->start_x/Y and data->end_x/Y as selection
   augmented_pae_info_t::box_info_t latest_box = data->boxes.back();
   // data->pae_info.pae_vecs.size();
   unsigned int box_index = data->boxes.size() - 1;
   int imol = data->imol;
   pae_colour_ribbon_model(imol, box_index, latest_box, data->pae_info.pae_vecs);
}


void
display_pae_from_file_in_a_dialog(int imol, const std::string &file_name) {


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

      augmented_pae_info_t *augmented_pae_p = static_cast<augmented_pae_info_t *>(data);
      const pae_t &pae_info = augmented_pae_p->pae_info;
      int n_pixels_for_pae_image = augmented_pae_p->pae_info.n_pixels - 100;

      // move the PAE graph in the base image
      int offset_x = 80;
      int offset_y = 20;
      unsigned int n_residues = augmented_pae_p->pae_info.pae_vecs.size();

      // let's make the actual data plot as image data:
      unsigned char *image_data = new unsigned char[n_pixels_for_pae_image*n_pixels_for_pae_image*4];
      float max_value = augmented_pae_p->pae_info.get_max_value(augmented_pae_p->pae_info.pae_vecs); // silly
      for (int i=0; i<n_pixels_for_pae_image; i++) {
         for (int j=0; j<n_pixels_for_pae_image; j++) {
            int idx = 4*(i*n_pixels_for_pae_image+j);
            float f_x = static_cast<float>(i)/static_cast<float>(n_pixels_for_pae_image);
            float f_y = static_cast<float>(j)/static_cast<float>(n_pixels_for_pae_image);
            int i_res_index = static_cast<int>(f_x * static_cast<float>(n_residues));
            int j_res_index = static_cast<int>(f_y * static_cast<float>(n_residues));
            float v = static_cast<float>(augmented_pae_p->pae_info.pae_vecs[i_res_index][j_res_index]);
            auto col = augmented_pae_p->pae_info.value_to_colour(v, max_value);
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
         auto col = pae_info.value_to_colour(epe, max_value);
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

      // if (augmented_pae_p->dragging_mode) {
      augmented_pae_p->cr = cr;
      draw_drag_box(augmented_pae_p);
      // }
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

      graphics_info_t g;
      std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = g.active_atom_spec();
      if (pp.first) {
	 int imol = pp.second.first;
	 augmented_pae_info_t *augmented_pae_p = new augmented_pae_info_t(imol, pae);
	 gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(drawing_area), draw_pae,
					static_cast<gpointer>(augmented_pae_p), NULL);
	 augmented_pae_p->dragging_mode = false;

	 GtkGesture *drag = gtk_gesture_drag_new();
	 gtk_widget_add_controller(drawing_area, GTK_EVENT_CONTROLLER(drag));
	 g_signal_connect(drag, "drag-begin",  G_CALLBACK(on_drag_begin),  augmented_pae_p);
	 g_signal_connect(drag, "drag-update", G_CALLBACK(on_drag_update), augmented_pae_p);
	 g_signal_connect(drag, "drag-end",    G_CALLBACK(on_drag_end),    augmented_pae_p);

	 // attach the pae_p (set "pae_p") to the pae_p so that it can be delete on closing the window.

	 gtk_widget_set_visible(window, TRUE);
	 graphics_info_t g;
	 g.set_transient_for_main_window(window);

      } else {
	 std::cout << "WARNING:: display_pae_from_file_in_a_dialog(): no model " << std::endl;
      }
   }

}
