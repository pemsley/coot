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
#include "cc-interface.hh"
#include "c-interface.h" // set contour level

class positron_plot_user_click_info_t {
public:
   positron_plot_user_click_info_t() : x(-1),  y(-1), imol_map(-1) {}
   positron_plot_user_click_info_t(double xx, double yy, int idx) : x(xx),  y(yy), imol_map(idx) {}
   double x, y; // in canvas coordinates (so 0-512)
   int imol_map; // the index of the resulting map
};

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
   GtkWidget *interpolation_entry;
   GtkWidget *animate_switch;
   GtkWidget *animate_spin_button;
   GtkWidget *animate_reverse_button;
   enum class animate_direction_t { FORWARDS, BACKWARDS};
   animate_direction_t animate_direction;
   cairo_t *cairo;
   cairo_surface_t *image_surface;
   double data[512][512];
   double min_v, max_v;
   int window_size_x, window_size_y;
   std::vector<float> map_colour;
   coot::positron_metadata_container_t mdc;
   std::vector<positron_plot_user_click_info_t> user_clicks;
   int animate_function_user_clicks_index;
   double colour_scale_factor; // this can be user-setable.
   float default_contour_level;
   int timeout_function_handle;
   bool stop_animation_flag;
   std::vector<int> basis_set_map_list;
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
      interpolation_entry = nullptr;
      animate_switch = nullptr;
      animate_spin_button = nullptr;
      animate_reverse_button = nullptr;
      cairo = nullptr;
      image_surface = nullptr;
      default_contour_level = 0.02;
      animate_function_user_clicks_index = -1;
      map_colour = {0.7, 0.7, 0.8};
      animate_direction = animate_direction_t::FORWARDS;
      stop_animation_flag = false;
      timeout_function_handle = -1;
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
   void set_spin_button(GtkWidget *sb) { animate_spin_button = sb; }
   void set_reverse_button(GtkWidget *rb) { animate_reverse_button = rb; }
   void set_animate_switch(GtkWidget *s) { animate_switch = s; }
   void set_interpolation_entry(GtkWidget *e) { interpolation_entry = e; }

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

   // positron data is stored as floats
   std::pair<float, float>
   canvas_coords_to_positron_coords(double cx, double cy) const {

      double f_x = cx/static_cast<double>(window_size_x);
      double f_y = cy/static_cast<double>(window_size_y);
      std::pair<float, float> z(min_v + f_x * (max_v-min_v),
                                min_v + f_y * (max_v-min_v));
      return z;
   }

   // not every user-click results in a map.
   unsigned int add_user_clicked_point_and_make_map(double x, double y) {
      std::cout << "user_click " << x << " " << y << std::endl;
      int imol_map_new = make_map(x,y);
      positron_plot_user_click_info_t puci(x,y,imol_map_new);
      user_clicks.push_back(puci);
      return user_clicks.size() - 1; // the index of the just-added point
   }

   void remove_last_user_click() {
      if (! user_clicks.empty()) {
         positron_plot_user_click_info_t puci = user_clicks.back();
         close_molecule(puci.imol_map);
         user_clicks.pop_back();
      }
      if (user_clicks.empty())
         animate_function_user_clicks_index = -1;
   }

   void clear() {
      for (unsigned int i=0; i<user_clicks.size(); i++) {
         int imol_map = user_clicks[i].imol_map;
         if (is_valid_map_molecule(imol_map))
            close_molecule(imol_map);
      }
      user_clicks.clear();
      animate_function_user_clicks_index = -1;
      if (animate_switch)
         gtk_switch_set_active(GTK_SWITCH(animate_switch), FALSE);
   }

   std::vector<std::pair<int, float> > make_weighted_map_indices(const coot::positron_metadata_t &pmdi) const {
      std::vector<std::pair<int, float> > v;
      if (pmdi.params.size() == basis_set_map_list.size()) {
         for (unsigned int i=0; i<pmdi.params.size(); i++) {
            std::pair<int, float> p(basis_set_map_list[i], pmdi.params[i]);
            v.push_back(p);
         }
      }
      return v;
   };

   // x and y are canvas coords
   int make_map(double x, double y) {
      int imol_map_new = -1;
      std::pair<float, float> z = canvas_coords_to_positron_coords(x, y);
      int idx_close = mdc.get_closest_positron_metadata_point(z);
      if (false) {
         std::cout << "in make_map() we have canvas coords: x and y " << x << " " << y << std::endl;
         std::cout << "in make_map() we have z " << z.first << " " << z.second << std::endl;
         std::cout << "in make_map() we have meta data mdc size " << mdc.size() << std::endl;
         std::cout << "in make_map() we have idx_close: " << idx_close << std::endl;
      }
      if (idx_close != -1) {
         coot::positron_metadata_t pmdi = mdc.metadata[idx_close];
         std::vector<std::pair<int, float> > weighted_map_indices = make_weighted_map_indices(pmdi);
         imol_map_new = make_weighted_map_simple_internal(weighted_map_indices);
         if (imol_map_new != -1) {
            std::string name = "Particle " + std::to_string(idx_close);
            set_map_colour(imol_map_new, map_colour[0], map_colour[1], map_colour[2]);
            set_contour_level_absolute(imol_map_new, default_contour_level);
            set_molecule_name(imol_map_new, name.c_str());
         }
      } else {
         std::cout << "No map created for this click" << std::endl;
      }
      return imol_map_new;
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

            if (false)
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
         // std::cout << "########### cairo_surface_status() success " << std::endl;
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

   int get_time_step_from_spinbutton() {
      if (animate_spin_button) // animate_time_step_spin_button
         return gtk_spin_button_get_value(GTK_SPIN_BUTTON(animate_spin_button));
      return 50;
   }

   void stop_animation() {
      stop_animation_flag = true;
   }

   int animate_timeout_func_inner(int continuous_mode) {

      // std::cout << "animate_timeout_func_inner() --- start ---" << std::endl;

      if (stop_animation_flag) {
         stop_animation_flag = false; // reset
         animate_function_user_clicks_index = -1; // ready to run again.
         return 0;
      }

      if (user_clicks.size() < 2) return 0;

      const auto &uc = user_clicks[animate_function_user_clicks_index];
      int imol_map = uc.imol_map;
      int final_ending_value = 0;
      if (continuous_mode) final_ending_value = 1; // there is no end

      if (false)
         std::cout << "--------- animate_function_user_clicks_index " << animate_function_user_clicks_index
                   << " user-clicks size(): " << user_clicks.size()
                   << " dir " << static_cast<int>(animate_direction) << " display-only-map " << imol_map
                   << std::endl;
      undisplay_all_maps_except(imol_map);
      int uc_size = user_clicks.size();
      if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(animate_reverse_button))) {
         if (animate_direction == animate_direction_t::FORWARDS) {
            animate_function_user_clicks_index++;
            if (animate_function_user_clicks_index == uc_size) {
               animate_function_user_clicks_index--;
               animate_direction = animate_direction_t::BACKWARDS;
               return 1;
            } else {
               return 1;
            }
         } else {
            animate_function_user_clicks_index--;
            if (final_ending_value == 0) {
               if (animate_function_user_clicks_index == -1) {
                  animate_direction = animate_direction_t::FORWARDS;
                  return final_ending_value;
               } else {
                  return 1;
               }
            } else {
               if (animate_function_user_clicks_index == -1) {
                  animate_direction = animate_direction_t::FORWARDS;
                  animate_function_user_clicks_index = 0; // correct out of bounds
               }
               return final_ending_value;
            }
         }
      } else {
         animate_function_user_clicks_index++;
         if (animate_function_user_clicks_index == uc_size) {
            animate_function_user_clicks_index = -1;  // so that we can start another animation
            if (final_ending_value == 1)
               animate_function_user_clicks_index = 0; // correct out of bouunds
            animate_direction = animate_direction_t::FORWARDS;
            return final_ending_value;
         } else {
            return 1; // keep going
         }
      }
      return final_ending_value;
   }

   int single_pass_animate_timeout_func_inner() {
      return animate_timeout_func_inner(0);
   }

   static int single_pass_animate_timeout_func(gpointer user_data) {

      plot_data_t *plot_data_p = static_cast<plot_data_t *>(user_data);
      return plot_data_p->single_pass_animate_timeout_func_inner();
   }

   int continuous_animation_timeout_func_inner() {
      return animate_timeout_func_inner(1);
   }

   static int continuous_animation_timeout_func(gpointer user_data) {

      plot_data_t *plot_data_p = static_cast<plot_data_t *>(user_data);
      return plot_data_p->continuous_animation_timeout_func_inner();
   }

   void single_pass_animate() {
      int time_step = get_time_step_from_spinbutton();
      if (animate_function_user_clicks_index == -1) {
         animate_function_user_clicks_index = 0;
         int tof = g_timeout_add(time_step, GSourceFunc(single_pass_animate_timeout_func), this);
         timeout_function_handle = tof;
      } else {
         std::cout << "active animation trap single-pass " << animate_function_user_clicks_index << std::endl;
      }
   }

   void continuous_animation() {
      int time_step = get_time_step_from_spinbutton();
      if (animate_function_user_clicks_index == -1) {
         animate_function_user_clicks_index = 0;
         int tof = g_timeout_add(time_step, GSourceFunc(continuous_animation_timeout_func), this);
         timeout_function_handle = tof;
      } else {
         std::cout << "active animation trap continuous " << animate_function_user_clicks_index << std::endl;
      }
   }

   int get_n_interpolation_points() const {
      const char *t = gtk_editable_get_text(GTK_EDITABLE(interpolation_entry));
      if (t) {
         std::string s(t);
         try {
            int rv = coot::util::string_to_int(s);
            return rv;
         }
         catch (const std::runtime_error &e) {
            std::cout << "WARNING::" << e.what() << std::endl;
         }
      }
      return 10;
   }

   void interpolate_clicks() {

      int n_interpolation_points = get_n_interpolation_points();

      if (user_clicks.size() * n_interpolation_points > 200) {
         std::cout << "Too many maps " << std::endl;
         return;
      }

      if (user_clicks.size() > 1) {
         std::vector<positron_plot_user_click_info_t> new_user_clicks;
         unsigned int lim = user_clicks.size() - 1;
         for (unsigned int i=0; i<lim; i++) {
            unsigned int idx_this = i;
            unsigned int idx_next = i+1;
            const auto &user_click_this = user_clicks[idx_this];
            const auto &user_click_next = user_clicks[idx_next];
            for (int j=0; j<n_interpolation_points; j++) {
               float f_in_range = static_cast<float>(j)/static_cast<float>(n_interpolation_points);
               float xx = user_click_this.x + (user_click_next.x - user_click_this.x) * f_in_range;
               float yy = user_click_this.y + (user_click_next.y - user_click_this.y) * f_in_range;
               std::pair<float, float> z = canvas_coords_to_positron_coords(xx, yy);
               int idx_close = mdc.get_closest_positron_metadata_point(z);
               // std::cout << "xx " << xx  << " yy " << yy << " idx_close " << idx_close << std::endl;
               if (idx_close != -1) {
                  int imol_map = make_map(xx, yy); // in canvas coords
                  if (imol_map != -1) {
                     positron_plot_user_click_info_t nuc(xx, yy, imol_map);
                     new_user_clicks.push_back(nuc);
                  }
               }
            }
         }
         new_user_clicks.push_back(user_clicks.back()); // add the very end point
         user_clicks = new_user_clicks;
         gtk_widget_queue_draw(drawing_area);
      }
   }

};

void on_draw_positron_plot(GtkDrawingArea *area,
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

   if (false)
      std::cout << "------------- on_draw_positron_plot(): test cairo_arc " << width/2.0 << " " << height/2.0
                << " " << MIN(width, height) / 10.0 << std::endl;

#if GTK_MINOR_VERSION >= 10
   gtk_widget_get_color(GTK_WIDGET(area), &color);
#endif
   gdk_cairo_set_source_rgba(cr, &color);
   cairo_fill(cr);

   plot_data_t *pdp = static_cast<plot_data_t *>(user_data);
   pdp->cairo = cr;
   std::cout << "---------------- starting with plot_data_p " << pdp << std::endl;

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
      if (pdp->image_surface) {
         int w = cairo_image_surface_get_width(pdp->image_surface);
         int h = cairo_image_surface_get_height(pdp->image_surface);
         // std::cout << "on_draw_positron_plot() image_surface w h " << w << " " << h << std::endl;
         cairo_set_source_surface(cr, pdp->image_surface, 0,0);
         cairo_paint(cr);
      } else {
         std::cout << "on_draw_positron_plot(): null surface " << std::endl;
      }
   }
   if (true) {
      color.alpha = 1.0;
      for (unsigned int i=0; i<pdp->user_clicks.size(); i++) {
         double di = pdp->user_clicks[i].x;
         double dj = pdp->user_clicks[i].y;
         int imol = pdp->user_clicks[i].imol_map;
         // std::cout << "drawing user_click " << i << " " << di << " " << dj << std::endl;
         if (imol == -1) {
            color.red   = 0.8;
            color.green = 0.8;
            color.blue  = 0.8;
         } else {
            color.red   = 0.2;
            color.green = 0.3;
            color.blue  = 0.5;
         }
         cairo_arc(cr, di, dj, 6.75, 0.0, 2.0 * G_PI);
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

   // std::cout << "click! " << x << " " << y << " " << n_press << std::endl;

   plot_data_t *plot_data_p = static_cast<plot_data_t *>(user_data);
   const auto &plot_data(*plot_data_p);

   // make this a member function
   // double f_x = static_cast<float>(x)/static_cast<float>(plot_data.window_size_x);
   // ndouble f_y = static_cast<float>(y)/static_cast<float>(plot_data.window_size_y);
   // std::pair<float, float> z(plot_data.min_v + f_x * (plot_data.max_v-plot_data.min_v),
   //                           plot_data.min_v + f_y * (plot_data.max_v-plot_data.min_v));

   unsigned int user_click_idx = plot_data_p->add_user_clicked_point_and_make_map(x, y);
   std::cout << "user-click idx " << user_click_idx
             << " new map index: " << plot_data_p->user_clicks[user_click_idx].imol_map << std::endl;

   gtk_widget_queue_draw(plot_data.drawing_area);
}

#include "utils/coot-utils.hh"

extern "C" G_MODULE_EXPORT
void
on_positron_map_undo_button_clicked(GtkButton *button,
                                    gpointer   user_data) {

   void *obj = g_object_get_data(G_OBJECT(button), "plot-data");
   plot_data_t *plot_data_p = static_cast<plot_data_t *>(obj);

   plot_data_p->remove_last_user_click();
   gtk_widget_queue_draw(plot_data_p->drawing_area);
}

extern "C" G_MODULE_EXPORT
void
on_positron_map_clear_button_clicked(GtkButton *button,
                                     gpointer   user_data) {

   void *obj = g_object_get_data(G_OBJECT(button), "plot-data");
   plot_data_t *plot_data_p = static_cast<plot_data_t *>(obj);

   plot_data_p->clear();
   gtk_widget_queue_draw(plot_data_p->drawing_area);
}

extern "C" G_MODULE_EXPORT
void
on_positron_interpolate_button_clicked(GtkButton *button,
                                       gpointer   user_data) {

   std::cout << "------------- button clicked " << std::endl;

   void *obj = g_object_get_data(G_OBJECT(button), "plot-data");
   plot_data_t *plot_data_p = static_cast<plot_data_t *>(obj);
   if (plot_data_p)
      plot_data_p->interpolate_clicks();

}

extern "C" G_MODULE_EXPORT
void
on_positron_animate_single_pass_button_clicked(GtkButton *button,
                                               gpointer   user_data) {

   void *obj = g_object_get_data(G_OBJECT(button), "plot-data");
   plot_data_t *plot_data_p = static_cast<plot_data_t *>(obj);
   if (plot_data_p)
      plot_data_p->single_pass_animate();
}

extern "C" G_MODULE_EXPORT
void
on_positron_animate_switch_activate(GtkSwitch *sw, gpointer user_data) {

   void *obj = g_object_get_data(G_OBJECT(sw), "plot-data");
   plot_data_t *plot_data_p = static_cast<plot_data_t *>(obj);
   if (plot_data_p) {
      if (gtk_switch_get_active(sw)) {
         std::cout << ".... start continuous animation" << std::endl;
         plot_data_p->continuous_animation();
      } else {
         plot_data_p->stop_animation();
         std::cout << ".... stop animation" << std::endl;
      }
   }
}


#ifdef STANDALONE_POSITRON_PLOT
GtkWidget *widget_from_builder(const std::string &w_name, GtkBuilder *builder) {
   GtkWidget *w = GTK_WIDGET(gtk_builder_get_object(GTK_BUILDER(builder), w_name.c_str()));
   return  w;
}
#else
GtkWidget *widget_from_builder(const std::string &w_name, GtkBuilder *builder);
#endif

#ifdef STANDALONE_POSITRON_PLOT

int main(int argc, char *argv[]) {

   gtk_init();

   GtkApplication *app = gtk_application_new ("org.emsley.coot", (GApplicationFlags) (G_APPLICATION_NON_UNIQUE));

   // GtkWidget *window = gtk_window_new();
   // gtk_window_set_title(GTK_WINDOW(window), "Coot-Positron 2D Array Plot");

   const std::vector<int> base_map_index_list; // fill this
   positron_plot_internal("metadata_z.csv", "metadata_s.csv", base_map_index_list);

   application_activate_data *activate_data = new application_activate_data(argc, argv);
   g_signal_connect(app, "activate", G_CALLBACK(positron_plot_application_activate), activate_data);

   int status_run = g_application_run(G_APPLICATION(app), 1, argv);
   std::cout << "g_application_run returns with status " << status_run << std::endl;

   return status;
}

#endif

void positron_plot_internal(const std::string &fn_z_csv, const std::string &fn_s_cvs,
                            const std::vector<int> &base_map_index_list) {

   GtkBuilder *builder = gtk_builder_new();
   std::string ui_file_name = "positron.ui";
   std::string dir = coot::package_data_dir();
   std::string dir_ui = coot::util::append_dir_dir(dir, "ui");
   std::string ui_file_full = coot::util::append_dir_file(dir_ui, ui_file_name);
   if (coot::file_exists(ui_file_name))
      ui_file_full = ui_file_name;
   GError* error = NULL;
   gboolean status = gtk_builder_add_from_file(builder, ui_file_full.c_str(), &error);
   if (status == FALSE) {
      std::cout << "ERROR:: Failure to read or parse " << ui_file_full << std::endl;
      std::cout << error->message << std::endl;
      exit(0);
   }

   GtkWidget *dialog = widget_from_builder("positron-dialog", builder);

   int plot_window_x_size = 512;
   int plot_window_y_size = 660;
   gtk_window_set_default_size(GTK_WINDOW(dialog), plot_window_x_size, plot_window_y_size);

   // GtkWidget *drawing_area = gtk_drawing_area_new();
   GtkWidget *spin_button    = widget_from_builder("positron_animate_time_step_spinbutton",      builder);
   GtkWidget *drawing_area   = widget_from_builder("positron_drawing_area",                      builder);
   GtkWidget *animate_switch = widget_from_builder("positron_animate_switch",                    builder);
   GtkWidget *reverse_button = widget_from_builder("positron_animate_with_reverse_togglebutton", builder);
   GtkWidget *interpolation_entry = widget_from_builder("positron_interpolate_steps_per_click_entry", builder);
   gtk_drawing_area_set_content_width( GTK_DRAWING_AREA(drawing_area), 512);
   gtk_drawing_area_set_content_height(GTK_DRAWING_AREA(drawing_area), 512);

   GtkAdjustment *adjustment = gtk_adjustment_new (50.0, 0.0, 250.0, 1.0, 5.0, 0.0);
   gtk_spin_button_set_adjustment(GTK_SPIN_BUTTON(spin_button), adjustment);

   plot_data_t *plot_data_p = new plot_data_t;
   plot_data_p->set_window_size(plot_window_x_size, plot_window_y_size);
   plot_data_p->set_drawing_area(drawing_area);
   plot_data_p->set_spin_button(spin_button);
   plot_data_p->set_reverse_button(reverse_button);
   plot_data_p->set_animate_switch(animate_switch);
   plot_data_p->set_interpolation_entry(interpolation_entry);
   plot_data_p->basis_set_map_list = base_map_index_list;
   plot_data_p->default_contour_level = 0.03;

   plot_data_p->fill_plot_data_from_positron_metadata(fn_z_csv, fn_s_cvs);
   unsigned char *image_data = new unsigned char[512*512*4]; // something something stride.
   cairo_surface_t *surface = plot_data_p->make_image_from_plot_data(image_data);
   plot_data_p->image_surface = surface;

   GtkWidget *undo_button        = widget_from_builder("positron_map_undo_button",  builder);
   GtkWidget *clear_button       = widget_from_builder("positron_map_clear_button", builder);
   GtkWidget *single_pass_button = widget_from_builder("positron_animate_single_pass_button", builder);
   GtkWidget *interpolate_button = widget_from_builder("positron_interpolate_button", builder);
   g_object_set_data(G_OBJECT(undo_button),        "plot-data", plot_data_p);
   g_object_set_data(G_OBJECT(clear_button),       "plot-data", plot_data_p);
   g_object_set_data(G_OBJECT(animate_switch),     "plot-data", plot_data_p);
   g_object_set_data(G_OBJECT(single_pass_button), "plot-data", plot_data_p);
   g_object_set_data(G_OBJECT(interpolate_button), "plot-data", plot_data_p);

   gpointer user_data = static_cast<void *>(plot_data_p);
   gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(drawing_area), on_draw_positron_plot, user_data, NULL);

   // gtk_window_set_child(GTK_WINDOW(window), drawing_area);
   gtk_widget_set_visible(dialog, TRUE);

   GtkGesture *drag_controller_primary = gtk_gesture_drag_new();
   gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(drag_controller_primary), GDK_BUTTON_PRIMARY);
   gtk_widget_add_controller(GTK_WIDGET(drawing_area), GTK_EVENT_CONTROLLER(drag_controller_primary));
   g_signal_connect(drag_controller_primary, "drag-begin",  G_CALLBACK(on_positron_plot_drag_begin_primary),  user_data);
   g_signal_connect(drag_controller_primary, "drag-update", G_CALLBACK(on_positron_plot_drag_update_primary), user_data);
   g_signal_connect(drag_controller_primary, "drag-end",    G_CALLBACK(on_positron_plot_drag_end_primary),    user_data);

   GtkGesture *click_controller = gtk_gesture_click_new();
   gtk_widget_add_controller(GTK_WIDGET(drawing_area), GTK_EVENT_CONTROLLER(click_controller));
   g_signal_connect(click_controller, "pressed",  G_CALLBACK(on_positron_plot_click), user_data);

}

#ifdef USE_PYTHON
void positron_plot_py(const std::string &fn_z_csv, const std::string &fn_s_csv,
                      PyObject *base_map_index_list) {

   std::vector<int> v;
   if (PyList_Check(base_map_index_list)) {
      long l = PyObject_Length(base_map_index_list);
      for (long i=0; i<l; i++) {
         PyObject *o = PyList_GetItem(base_map_index_list, i);
         long idx = PyLong_AsLong(o);
         v.push_back(idx);
      }
   }
   positron_plot_internal(fn_z_csv, fn_s_csv, v);
}
#endif
