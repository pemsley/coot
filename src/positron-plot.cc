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
      cm_tropical[10].red = 0.515; cm_tropical[10].green = 0.024; cm_tropical[10].blue = 0.600; cm_tropical[10].alpha = 1.0;
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
   double data[512][512];
   std::vector<GdkRGBA> cm_expanded;
   std::vector<GdkRGBA> cm_tropical;
   std::vector<GdkRGBA> expand_colour_map(const std::vector<GdkRGBA> &cm_in, unsigned int n_bins) {
      std::vector<GdkRGBA> v(n_bins);
      unsigned int n_ranges = cm_in.size() -1;
      std::cout << "n_ranges: " << n_ranges << std::endl;
      unsigned int n_bins_per_range = n_bins/n_ranges;
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
      set_cm_tropical();
      cm_expanded = expand_colour_map(cm_tropical, 100);
   }
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

   void fill_plot_data_from_positron_metadata(const std::string &fn1, const std::string &fn2) {
      std::cout << "fill_plot_data_from_positron_metadata() --- start ---" << std::endl;
      std::vector<coot::positron_metadata_t> d;
      // range of data is -8 to 5
      coot::read_positron_metadata(&d, fn1, fn2);
      std::cout << "fill_plot_data_from_positron_metadata() A d.size() " << d.size() << std::endl;

      for (int i=0; i<512; i++)
         for (int j=0; j<512; j++)
            data[i][j] = 0.0;

      std::cout << "fill_plot_data_from_positron_metadata() B d.size() " << d.size() << std::endl;
      double min_v = -8.0;
      double max_v =  5.0;
      double min_max_range = max_v - min_v;
      for (unsigned int i=0; i<d.size(); i++) {
         const auto &dd = d[i];
         float x = ((dd.x - min_v)/min_max_range) * 512.0;
         float y = ((dd.y - min_v)/min_max_range) * 512.0;
         int ix = x;
         int iy = y;
         // a few points may lie outside the desired plotting min and max. So they have out of bounds
         // indices, in which case just move on.
         if (ix < 0) continue;
         if (iy < 0) continue;
         if (iy > 511) continue;
         if (iy > 511) continue;
         // std::cout << "adding obs in fill-plot-data " << dd.x << " " << dd.y << " " << ix << " " << iy << "\n";
         data[ix][iy] += 1;
      }
   }

   GdkRGBA get_colour_for_value(double f) const {
      if (f < 0.0) f = 0.0;
      if (f > 0.99999) f = 0.99999;
      int f_idx = f * static_cast<float>(cm_expanded.size());
      int idx = f_idx;
      GdkRGBA r = cm_expanded[idx];
      // std::cout << "f " << f << " idx " << idx << "\n";
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

   gtk_widget_get_color(GTK_WIDGET(area), &color);
   gdk_cairo_set_source_rgba(cr, &color);
   cairo_fill(cr);

   plot_data_t *pdp = static_cast<plot_data_t *>(user_data);
   const plot_data_t &plot_data(*pdp);
   for(int i=0; i<512; i++) {
      for(int j=0; j<512; j++) {
         double di = width  * static_cast<double>(i)/512.0 + 5.0;
         double dj = height * static_cast<double>(j)/512.0 + 5.0;
         // cairo_rectangle(cr, di, dj, 0.01, 0.02);
         cairo_arc(cr, di, dj, 1.75, 0.0, 2.0 * G_PI);
         double fv = plot_data.data[i][j] * 0.1;
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
   application_activate_data* activate_data = (application_activate_data*) user_data;
   activate_data->application = application;

   std::string window_name("Positron Plot");
   GtkWidget *app_window = gtk_application_window_new(application);
   gtk_window_set_application(GTK_WINDOW(app_window), application);
   gtk_window_set_title(GTK_WINDOW(app_window), window_name.c_str());

}


int main(int argc, char *argv[]) {

   gtk_init();

   GtkApplication *app = gtk_application_new ("org.emsley.coot",
      (GApplicationFlags) (G_APPLICATION_NON_UNIQUE));

   GtkWidget *window = gtk_window_new();
   gtk_window_set_title(GTK_WINDOW(window), "2D Array Plot");
   gtk_window_set_default_size(GTK_WINDOW(window), 800, 600);

   GtkWidget *drawing_area = gtk_drawing_area_new();
   gtk_drawing_area_set_content_width(GTK_DRAWING_AREA(drawing_area), 512);
   gtk_drawing_area_set_content_height(GTK_DRAWING_AREA(drawing_area), 512);

   plot_data_t *plot_data_p = new plot_data_t;
   std::string plot_data_file_name = "plot.data";
   if (argc > 1)
      plot_data_file_name = argv[1];
   // plot_data_p->fill_plot_data(plot_data_file_name);
   std::string fn1 = "metadata_z.csv";
   std::string fn2 = "metadata_s.csv";
   plot_data_p->fill_plot_data_from_positron_metadata(fn1, fn2);

   gpointer user_data = static_cast<void *>(plot_data_p);
   gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(drawing_area), on_draw, user_data, NULL);

   gtk_window_set_child(GTK_WINDOW(window), drawing_area);
   gtk_widget_set_visible(window, TRUE);

   application_activate_data *activate_data = new application_activate_data(argc, argv);
   g_signal_connect(app, "activate", G_CALLBACK(positron_plot_application_activate), activate_data);

   int status = g_application_run(G_APPLICATION(app), 1, argv);
   std::cout << "g_application_run returns with status " << status << std::endl;

   return status;
}
