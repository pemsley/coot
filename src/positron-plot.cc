/*
 * src/positron-plot.cc
 *
 * Copyright 2024 by Medical Research Council
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

#include <iostream>
#include <fstream>
#include <gtk/gtk.h>
#include <glib-2.0/glib.h>
#include <cairo/cairo.h>
#include "coot-utils/positron.hh"


class plot_data_t {
   void set_cm_tropical() {
      cm_tropical.resize(11);
      cm_tropical[0].red = 1.0;   cm_tropical[0].green = 1.0;   cm_tropical[0].blue = 1.0;   cm_tropical[0].alpha = 1.0;
      cm_tropical[1].red = 0.267; cm_tropical[1].green = 0.987; cm_tropical[1].blue = 0.988; cm_tropical[1].alpha = 1.0;
      cm_tropical[2].red = 0.154; cm_tropical[2].green = 0.934; cm_tropical[2].blue = 0.722; cm_tropical[2].alpha = 1.0;
      cm_tropical[3].red = 0.429; cm_tropical[3].green = 0.843; cm_tropical[3].blue = 0.431; cm_tropical[3].alpha = 1.0;
      cm_tropical[4].red = 0.647; cm_tropical[4].green = 0.719; cm_tropical[4].blue = 0.203; cm_tropical[4].alpha = 1.0;
      cm_tropical[5].red = 0.772; cm_tropical[5].green = 0.580; cm_tropical[5].blue = 0.031; cm_tropical[5].alpha = 1.0;
      cm_tropical[6].red = 0.837; cm_tropical[6].green = 0.429; cm_tropical[6].blue = 0.067; cm_tropical[6].alpha = 1.0;
      cm_tropical[7].red = 0.850; cm_tropical[7].green = 0.273; cm_tropical[7].blue = 0.195; cm_tropical[7].alpha = 1.0;
      cm_tropical[8].red = 0.808; cm_tropical[8].green = 0.111; cm_tropical[8].blue = 0.354; cm_tropical[8].alpha = 1.0;
      cm_tropical[9].red = 0.699; cm_tropical[9].green = 0.022; cm_tropical[9].blue = 0.528; cm_tropical[9].alpha = 1.0;
      cm_tropical[10].red = 0.415; cm_tropical[10].green = 0.024; cm_tropical[10].blue = 0.500; cm_tropical[10].alpha = 1.0;
      // cm_tropical[10].red = 0.565; cm_tropical[10].green = 0.054; cm_tropical[10].blue = 0.646; cm_tropical[10].alpha = 1.0;
   }
   GdkRGBA get_colour_mix(const GdkRGBA &s, const GdkRGBA &e, double f) {
      double r = s.red   * (1.0 - f) + e.red * f;
      double g = s.green * (1.0 - f) + e.green * f;
      double b = s.blue  * (1.0 - f) + e.blue * f;
      GdkRGBA c;
      c.red = r; c.green = g; c.blue = b; c.alpha = 1.0;
      return c;
   }
   public:
   plot_data_t() { init(); }
   GtkWidget *drawing_area;
   cairo_t *cairo;
   cairo_surface_t *image_surface;
   double data[512][512];
   double min_v, max_v;
   int window_size_x, window_size_y;
   coot::positron_metadata_container_t mdc;
   std::vector<std::pair<double, double> > user_clicks;
   double colour_scale_factor; // this can be user-setable.
   std::vector<GdkRGBA> cm_expanded;
   std::vector<GdkRGBA> cm_tropical;
   std::vector<GdkRGBA> expand_colour_map(const std::vector<GdkRGBA> &cm_in, unsigned int n_bins) {
      std::vector<GdkRGBA> v(n_bins);
      unsigned int n_ranges = cm_in.size() -1;
      // std::cout << "n_ranges: " << n_ranges << std::endl;
      for (unsigned int i=0; i<n_bins; i++) {
         float fr = static_cast<float>(i)/static_cast<float>(n_bins);
         unsigned int range_index = n_ranges * i/n_bins;
         // std::cout << "Range index for i " << i << " " << range_index << std::endl;
         const GdkRGBA &s = cm_in[range_index];
         const GdkRGBA &e = cm_in[range_index+1];
         float range_index_fr_start = static_cast<float>(range_index)/static_cast<float>(n_ranges);
         float f_delta = fr - range_index_fr_start;
         GdkRGBA c_mix = get_colour_mix(s,e,f_delta);
         v[i] = c_mix;
      }
      return v;
   }
   void init() {
      for (int i=0; i<512; i++)
         for (int j=0; j<512; j++)
            data[i][j] = 0.0;

      min_v = 0;
      max_v = 0;
      window_size_x = 0;
      window_size_y = 0;
      set_cm_tropical();
      cm_expanded = expand_colour_map(cm_tropical, 100);
      colour_scale_factor = 0.08;
      drawing_area = nullptr;
      cairo = nullptr;
      image_surface = nullptr;
   }
   void cleanup() {
      if (image_surface)
         cairo_surface_destroy(image_surface);
   }
   void set_window_size(int x, int y) {
      window_size_x = x;
      window_size_y = y;
   }
   void set_drawing_area(GtkWidget *da) { drawing_area = da; }
   void fill_plot_data(const std::string &plot_data_file_name) {
      if (g_file_test(plot_data_file_name.c_str(), G_FILE_TEST_IS_REGULAR)) {
         std::ifstream f(plot_data_file_name);
         int x, y;
         double fv;
         while (f) {
            f >> x >> y >> fv;
            data[x][y] = fv;
         }
      }
   }
   // this should be a member function of positron_metadata_container_t
   std::pair<double, double> get_sane_min_v_max_v(const std::vector<coot::positron_metadata_t> &d) const {
      // return mdc.get_sane_min_v_max_v();
      return std::make_pair(-7.0, 6.0);
   }

   void add_user_clicked_point(double x, double y) {
      user_clicks.push_back(std::make_pair(x,y));
   }

   void add_user_clicked_point(cairo_t *cr, double x, double y) {
      if (drawing_area) {
         std::cout << "add_user_clicked_point() B " << x << " " << y << " drawing_area: " << drawing_area << std::endl;
         GdkRGBA color;
         color.red   = 0.8;
         color.green = 0.2;
         color.blue  = 0.2;
         color.alpha = 1.0;
         gdk_cairo_set_source_rgba(cr, &color);
         cairo_arc(cr, x, y, 23.0 , 0.0, 2.0 * G_PI);
         cairo_fill(cr);
         gtk_widget_queue_draw(drawing_area);
      }
   }

   void fill_plot_data_from_positron_metadata(const std::string &fn1, const std::string &fn2) {

      coot::read_positron_metadata(&mdc.metadata, fn1, fn2);
      // range of data is -8 to 5
      std::pair<double, double> sane_min_max = get_sane_min_v_max_v(mdc.metadata);
      min_v = sane_min_max.first;
      max_v = sane_min_max.second;

      double min_max_range = max_v - min_v;
      for (unsigned int i=0; i<mdc.metadata.size(); i++) {
         const auto &dd = mdc.metadata[i];
         float x = ((dd.x - min_v)/min_max_range) * 512.0;
         float y = ((dd.y - min_v)/min_max_range) * 512.0;
         int ix = x;
         int iy = y;
         // a few points may lie outside the desired plotting min and max. So they have out of bounds
         // indices, in which case just move on.
         if (ix < 0)   continue;
         if (iy < 0)   continue;
         if (ix > 511) continue;
         if (iy > 511) continue;
         data[ix][iy] += 1;
      }
   }

   // caller handles the memory
   cairo_surface_t *make_image_from_plot_data(unsigned char *image_data) {

      int str = 4;
      for (int ix=0; ix<512; ix++) {
         for (int iy=0; iy<512; iy++) {
            int idx = ix * 512 * str + iy * str;
            const double &f = data[ix][iy];
            GdkRGBA c = get_colour_for_value(f * colour_scale_factor);

            image_data[idx]   = 255;
            image_data[idx+1] = 255;
            image_data[idx+2] = 255;

            int c_r = static_cast<int>(255.0 * c.red);
            int c_g = static_cast<int>(255.0 * c.green);
            int c_b = static_cast<int>(255.0 * c.blue);

            image_data[idx]   = c_b; // blue
            image_data[idx+1] = c_g; // green
            image_data[idx+2] = c_r; // red
            image_data[idx+3] = 255; // ignored?

            if (true)
               if (f > 0.2)
                  std::cout << f << " " << ix << " " << iy << " col "
                            << c_r << " " << c_g << " " << c_b << " image-data: "
                            << static_cast<int>(image_data[idx]) << " "
                            << static_cast<int>(image_data[idx+1]) << " "
                            << static_cast<int>(image_data[idx+2]) << std::endl;
         }
      }

      int stride = cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, 512);
      std::cout << "stride: " << stride << std::endl;
      cairo_surface_t *surface = cairo_image_surface_create_for_data(image_data, CAIRO_FORMAT_RGB24, 512, 512, 512*4);

      if (cairo_surface_status(surface) == CAIRO_STATUS_SUCCESS) {
         std::cout << "########### cairo_surface_status() success " << std::endl;
      } else {
         std::cout << "########### cairo_surface_status() fail " << std::endl;
      }

      if((cairo_surface_status(surface) != CAIRO_STATUS_SUCCESS)) {
         std::cout << "ERROR:: ", cairo_status_to_string(cairo_surface_status(surface));
         surface = nullptr;
      } else {
         int w = cairo_image_surface_get_width(surface);
         int h = cairo_image_surface_get_height(surface);
         std::cout << "make_image_from_plot_data(): image surface w h " << w << " " << h << std::endl;
      }
      return surface;
   };

   GdkRGBA get_colour_for_value(double f) const {
      if (f < 0.0) f = 0.0;
      if (f > 0.99999) f = 0.99999;
      int f_idx = f * static_cast<float>(cm_expanded.size());
      int idx = f_idx;
      GdkRGBA r = cm_expanded[idx];
      return r;
   }
};

void on_draw(GtkDrawingArea *area,
             cairo_t        *cr,
             int             width,
             int             height,
             gpointer        user_data) {

   auto function_value_to_color = [] (double fv, const plot_data_t &plot_data, GdkRGBA *c) {
         c->red   = 1.0 - fv;
         c->green = 1.0 - fv;
         c->blue  = 1.0 - fv;
         c->alpha = 1.0;
         *c = plot_data.get_colour_for_value(fv);
   };

   GdkRGBA color;

   cairo_arc (cr,
             width / 2.0, height / 2.0,
             MIN (width, height) / 10.0,
             0, 2 * G_PI);

   std::cout << "------------- on_draw(): test cairo_arc " << width/2.0 << " " << height/2.0
             << " " << MIN(width, height) / 10.0 << std::endl;

   gtk_widget_get_color(GTK_WIDGET(area), &color);
   gdk_cairo_set_source_rgba(cr, &color);
   cairo_fill(cr);

   plot_data_t *pdp = static_cast<plot_data_t *>(user_data);
   pdp->cairo = cr;

   const plot_data_t &plot_data(*pdp);
   float colour_scale_factor = plot_data.colour_scale_factor;

   if (false) {
      for(int i=0; i<512; i++) {
         for(int j=0; j<512; j++) {
            // continue;
            double di = width  * static_cast<double>(i)/512.0 + 5.0;
            double dj = height * static_cast<double>(j)/512.0 + 5.0;
            // cairo_rectangle(cr, di, dj, 0.01, 0.02);
            cairo_arc(cr, di, dj, 1.75, 0.0, 2.0 * G_PI);
            double fv = plot_data.data[i][j] * colour_scale_factor;
            if (false)
               if (i%10 == 0)
                  if (j%10 == 0)
                     std::cout << "Sample " << di << " " << dj << " " << plot_data.data[i][j] << std::endl;
            function_value_to_color(fv, plot_data, &color);
            gdk_cairo_set_source_rgba(cr, &color);
            cairo_fill(cr);
         }
      }
   }
   if (true) {
      std::cout << "on_draw() here with image_surface " << pdp->image_surface << std::endl;
      if (pdp->image_surface) {
         int w = cairo_image_surface_get_width(pdp->image_surface);
         int h = cairo_image_surface_get_height(pdp->image_surface);
         std::cout << "on_draw() image_surface w h " << w << " " << h << std::endl;
         cairo_set_source_surface(cr, pdp->image_surface, 0,0);
         cairo_paint(cr);
      } else {
         std::cout << "on_draw(): null surface " << std::endl;
      }
   }
   color.red   = 0.2;
   color.green = 0.2;
   color.blue  = 0.6;
   color.alpha = 1.0;
   for (unsigned int i=0; i<pdp->user_clicks.size(); i++) {
      double di = pdp->user_clicks[i].first;
      double dj = pdp->user_clicks[i].second;
      std::cout << "drawing user_click " << i << " " << di << " " << dj << std::endl;
      cairo_arc(cr, di, dj, 6.75, 0.0, 2.0 * G_PI);
      gdk_cairo_set_source_rgba(cr, &color);
      cairo_fill(cr);
   }
}


void quit(GtkWidget *w, gpointer user_data) {
   std::cout << "bye" << std::endl;
}

struct application_activate_data {
    int argc;
    char** argv;
    GtkWidget* splash_screen;
    GtkApplication* application;
    GtkWidget* app_window;
    // command_line_data cld;

    application_activate_data(int _argc, char** _argv) {
       argc = _argc;
       argv = _argv;
       splash_screen = nullptr;
       application = nullptr;
       app_window = nullptr;
    }
};


void
positron_plot_application_activate(GtkApplication *application,
                                   gpointer user_data) {

   application_activate_data* activate_data = static_cast<application_activate_data *>(user_data);
   activate_data->application = application;

   std::string window_name("Positron Plot");
   GtkWidget *app_window = gtk_application_window_new(application);
   gtk_window_set_application(GTK_WINDOW(app_window), application);
   gtk_window_set_title(GTK_WINDOW(app_window), window_name.c_str());

}

void on_positron_plot_drag_begin_primary(GtkGestureDrag *gesture,
                                         double          x,
                                         double          y,
                                         GtkWidget      *area) {

   // std::cout << "drag-begin " << x << " " << y << " " << std::endl;

}

void on_positron_plot_drag_update_primary(GtkGestureDrag *gesture,
                                          double          delta_x,
                                          double          delta_y,
                                          GtkWidget      *area) {

   // std::cout << "drag-update " << delta_x << " " << delta_y << " " << std::endl;

}

void on_positron_plot_drag_end_primary(GtkGestureDrag *gesture,
                                       double          x,
                                       double          y,
                                       GtkWidget      *area) {

   // std::cout << "drag-end " << x << " " << y << " " << std::endl;
}



void
on_positron_plot_click(GtkGestureClick* click_gesture,
                       gint n_press,
                       gdouble x,
                       gdouble y,
                       gpointer user_data) {

   std::cout << "click! " << x << " " << y << " " << n_press << std::endl;
   // undo this:
   // gpointer user_data = static_cast<void *>(plot_data_p);

   plot_data_t *plot_data_p = static_cast<plot_data_t *>(user_data);
   const auto &plot_data(*plot_data_p);

   // make this a member function
   double f_x = static_cast<float>(x)/static_cast<float>(plot_data.window_size_x);
   double f_y = static_cast<float>(x)/static_cast<float>(plot_data.window_size_y);
   std::pair<float, float> z(plot_data.min_v + f_x * (plot_data.max_v-plot_data.min_v),
                             plot_data.min_v + f_y * (plot_data.max_v-plot_data.min_v));
   int idx = plot_data.mdc.get_closest_positron_metadata_point(z); // -1 on failure
   std::cout << "idx: " << idx << std::endl;

   plot_data_p->add_user_clicked_point(x, y);

   if (plot_data_p->cairo)
      plot_data_p->add_user_clicked_point(plot_data.cairo, x, y);
}




int main(int argc, char *argv[]) {

   gtk_init();

   GtkApplication *app = gtk_application_new ("org.emsley.coot", (GApplicationFlags) (G_APPLICATION_NON_UNIQUE));

   GtkWidget *window = gtk_window_new();
   gtk_window_set_title(GTK_WINDOW(window), "Coot-Positron 2D Array Plot");
   int plot_window_x_size = 512;
   int plot_window_y_size = 700;
   gtk_window_set_default_size(GTK_WINDOW(window), plot_window_x_size, plot_window_y_size);

   GtkWidget *drawing_area = gtk_drawing_area_new();
   gtk_drawing_area_set_content_width( GTK_DRAWING_AREA(drawing_area), 512);
   gtk_drawing_area_set_content_height(GTK_DRAWING_AREA(drawing_area), 512);

   plot_data_t *plot_data_p = new plot_data_t;
   plot_data_p->set_window_size(plot_window_x_size, plot_window_y_size);
   plot_data_p->set_drawing_area(drawing_area);

   // std::string plot_data_file_name = "plot.data";
   //      if (argc > 1)
   // plot_data_file_name = argv[1];
   // plot_data_p->fill_plot_data(plot_data_file_name);

   std::string fn1 = "metadata_z.csv";
   std::string fn2 = "metadata_s.csv";
   plot_data_p->fill_plot_data_from_positron_metadata(fn1, fn2);
   unsigned char *image_data = new unsigned char[512*512*4]; // something something stride.
   cairo_surface_t *surface = plot_data_p->make_image_from_plot_data(image_data);
   plot_data_p->image_surface = surface;

   gpointer user_data = static_cast<void *>(plot_data_p);
   gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(drawing_area), on_draw, user_data, NULL);

   gtk_window_set_child(GTK_WINDOW(window), drawing_area);
   gtk_widget_set_visible(window, TRUE);

   GtkGesture *drag_controller_primary = gtk_gesture_drag_new();
   gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(drag_controller_primary), GDK_BUTTON_PRIMARY);
   gtk_widget_add_controller(GTK_WIDGET(drawing_area), GTK_EVENT_CONTROLLER(drag_controller_primary));
   g_signal_connect(drag_controller_primary, "drag-begin",  G_CALLBACK(on_positron_plot_drag_begin_primary),  user_data);
   g_signal_connect(drag_controller_primary, "drag-update", G_CALLBACK(on_positron_plot_drag_update_primary), user_data);
   g_signal_connect(drag_controller_primary, "drag-end",    G_CALLBACK(on_positron_plot_drag_end_primary),    user_data);

   GtkGesture *click_controller = gtk_gesture_click_new();
   gtk_widget_add_controller(GTK_WIDGET(drawing_area), GTK_EVENT_CONTROLLER(click_controller));
   g_signal_connect(click_controller, "pressed",  G_CALLBACK(on_positron_plot_click), user_data);

   application_activate_data *activate_data = new application_activate_data(argc, argv);
   g_signal_connect(app, "activate", G_CALLBACK(positron_plot_application_activate), activate_data);

   int status = g_application_run(G_APPLICATION(app), 1, argv);
   std::cout << "g_application_run returns with status " << status << std::endl;

   return status;
}
