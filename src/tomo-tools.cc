/*
 * src/tomo-tools.cc
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

#include <chrono>
#include <iomanip>

#include <clipper/ccp4/ccp4_map_io.h>
#include "utils/coot-utils.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/coot-map-heavy.hh"
#include "coot-utils/texture-as-floats.hh"
#include "coot-utils/mini-texture.hh"
#include "graphics-info.h"

//! negative becomes positive and positive becomes negative.
//! Apply an offset so that most of the map is above zero.
//!
void reverse_map(int imol_map) {

   graphics_info_t g;
   if (g.is_valid_map_molecule(imol_map)) {
      g.molecules[imol_map].reverse_map();
      g.graphics_draw();
   }
}


void tomo_section_view(int imol, int axis_id);

void set_tomo_section_view_section(int imol, int section_index) {

   GtkWidget *scale = widget_from_builder("tomo_scale");
   GtkAdjustment *adjustment = gtk_range_get_adjustment(GTK_RANGE(scale));
   if (adjustment) {
      gtk_adjustment_set_value(adjustment, section_index);
   }
}


void tomo_scale_adjustment_changed(GtkAdjustment *adj, gpointer user_data) {

   double v = gtk_adjustment_get_value(adj);
   int section = static_cast<int>(v);
   int imol = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(adj), "imol"));
   graphics_info_t g;
   g.set_tomo_section_view_section(imol, section);
   g.graphics_grab_focus();

}

int setup_tomo_slider(int imol) {

   int section = -1;
   graphics_info_t g;
   if (g.is_valid_map_molecule(imol)) {

      const auto &xmap = g.molecules[imol].xmap;

      // move this into graphics_info_t
      clipper::Grid_sampling gs = xmap.grid_sampling();
      int n_sections = gs.nw();
      int mid_section = n_sections/2;
      section = mid_section;

      GtkWidget *toolbar_vbox = widget_from_builder("main_window_vbox_inner");
      gtk_widget_set_visible(toolbar_vbox, FALSE);

      GtkWidget *tomo_box = widget_from_builder("tomo_slider_vbox");
      GtkWidget *scale = widget_from_builder("tomo_scale");
      gtk_widget_set_visible(tomo_box, TRUE);
      GtkAdjustment *adjustment = gtk_range_get_adjustment(GTK_RANGE(scale));
      if (adjustment) {
         g_object_set_data(G_OBJECT(adjustment), "imol", GINT_TO_POINTER(imol));
         std::string end_str = std::to_string(n_sections - 1);
         double dns_1 = static_cast<double>(n_sections-1);
         gtk_adjustment_set_lower(adjustment, 0);
         gtk_adjustment_set_upper(adjustment, n_sections-1);
         gtk_adjustment_set_step_increment(adjustment, 1.0);
         gtk_adjustment_set_page_increment(adjustment, 5.0);
         gtk_adjustment_set_page_size(adjustment, 0.0);
         gtk_adjustment_set_value(adjustment, mid_section);
         gtk_scale_set_draw_value(GTK_SCALE(scale), TRUE);
         gtk_scale_set_digits(GTK_SCALE(scale), 0);
         gtk_scale_add_mark(GTK_SCALE(scale), 0, GTK_POS_BOTTOM, "0");
         gtk_scale_add_mark(GTK_SCALE(scale), dns_1, GTK_POS_BOTTOM, end_str.c_str());
      }
   }
   return section;
}

void tomo_section_view(int imol, int section_index) {

   // auto tp_start = std::chrono::high_resolution_clock::now();
   graphics_info_t g;
   if (g.is_valid_map_molecule(imol)) {
      g.zoom = 20000.0f; // or should be in setup_tomo_slider()?
      g.set_tomo_section_view_section(imol, section_index);
      g.clipping_front = 4.5;       // c.f. scene_preset_figure_making_action()
      g.clipping_back  = 1.3;

      clipper::Cell cell = g.molecules[imol].xmap.cell();
      clipper::Coord_orth pt(cell.a() * 0.5f, cell.b() * 0.5, cell.c() * 0.5f);
      g.set_rotation_centre(pt);
   }
   // auto tp_now = std::chrono::high_resolution_clock::now();
   // auto delta = std::chrono::duration_cast<std::chrono::milliseconds>(tp_now - tp_start);

}

void set_tomo_picker_mode_is_active(short int state) {

   graphics_info_t::tomo_picker_flag = state;

}

#include "analysis/stats.hh"

#ifdef USE_PYTHON
void tomo_map_analysis(int imol_map, PyObject *spot_positions) {

   auto do_stats = [] (unsigned int point_index,
                       const clipper::Coord_orth &pt,
                       const clipper::Xmap<float> &xmap) {
      double delta_size = 10.0;
      for (int ix=-50; ix<50; ix++) {
         for (int iy=-5; iy<5; iy++) {
            double dix = static_cast<double>(ix);
            double diy = static_cast<double>(iy);
            clipper::Coord_orth p = pt + clipper::Coord_orth(dix * delta_size, diy * delta_size, 0.0);
            float f = coot::util::density_at_point(xmap, p);
            if (false)
               std::cout << "map-point " << point_index << " " << ix << " "
                         << dix * delta_size << " " << diy * delta_size << " "
                         << p.x() << " " << p.y() << " " << p.z() << " value: " << f << std::endl;
         }
      }
   };

   class profile_t {
      public:
      profile_t() { n_pixels_per_half_edge = -1; width = -1; height = -1; }
      profile_t(int half_width_n_pixels_in, int w, int h) : n_pixels_per_half_edge(half_width_n_pixels_in), width(w), height(h) {
         mean_and_var.resize(w*h);}
      int n_pixels_per_half_edge;
      int width;
      int height;
      std::vector<std::pair<float, float> > mean_and_var;
      void set_value(int w, int h, const float &f) { mean_and_var[width * w + h].first = f; }
      float get(int iw, int ih) const { return mean_and_var[width * (iw + n_pixels_per_half_edge) + ih + n_pixels_per_half_edge].first; }
      void add_profile(const profile_t &other) {
         int size = width * height;
         for (int i=0; i<size; i++) {
            // std::cout << "adding " << other.mean_and_var[i].first << " to " << mean_and_var[i].first << "\n";
            mean_and_var[i].first += other.mean_and_var[i].first;
         }
      }
   };

   // make this a member function of the above class
   auto print_profile = [] (const profile_t &prof, bool as_single_digits) {

      int n = prof.n_pixels_per_half_edge;
      if (as_single_digits) {

         float sf = 22.0; // scale factor
         float offset = +0.15;

         for (int iw=-n; iw<=n; iw++) {
            for (int ih=-n; ih<=n; ih++) {
               float ff = sf * (prof.get(iw, ih) + offset);
               int ii = static_cast<int>(ff);
               if (ii<0) {
                  std::cout << "- ";
               } else {
                  if (ii>9) ii = 9;
                  std::cout << ii << " ";
               }
            }
            std::cout << "\n";
         }
      } else {
         for (int iw=-n; iw<=n; iw++) {
            for (int ih=-n; ih<=n; ih++) {
               float ff = prof.get(iw, ih);
               std::cout << std::fixed << std::setw(6) << std::setprecision(3) << ff << " ";
            }
            std::cout << "\n";
         }
      }
   };

   auto map_scan = [] (const std::map<int, double> &offset_values_map, const clipper::Xmap<float> &xmap, int section_index) {

      clipper::Grid_sampling gs = xmap.grid_sampling();
      clipper::Coord_grid cg_0(0,0,section_index);
      clipper::Coord_grid cg_1(gs.nu()-1, gs.nv()-1, section_index);
      clipper::Grid_map grid(cg_0, cg_1);
      clipper::Xmap_base::Map_reference_coord ix( xmap, grid.min()), iu, iv, iw;
      int nv = gs.nv();
      int nu = gs.nu();
      clipper::Cell cell = xmap.cell();
      float x_size = cell.a();
      float y_size = cell.b();

      // now make the data for a mini-texture - this is the thing we return
      unsigned int nunv = nu * nv;
      std::vector<float> tafd(nunv, 0);
      int c_u = 0;
      for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
         int c_v = 0;
         for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
            int c_w = 0;
            for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
               // for Z-sections (c.f. mini-texture.cc)
               int img_x_coord = c_u; int img_y_coord = c_v; int img_n_rows = nv;
               int idx = (img_y_coord + img_n_rows * img_x_coord);
               clipper::Coord_grid cg = iw.coord();
               clipper::Coord_frac cf = cg.coord_frac(xmap.grid_sampling());
               clipper::Coord_orth co = cf.coord_orth(xmap.cell());
               int delta_size_i = 10;
               double delta_size_d = 10.0; // same
               float s = 0.0f;
               float sum_d = 0.0f;
               for (int iix=-50; iix<50; iix++) {
                  for (int iiy=-5; iiy<5; iiy++) {
                     double dix = static_cast<double>(iix);
                     double diy = static_cast<double>(iiy);
                     clipper::Coord_orth p = co + clipper::Coord_orth(dix * delta_size_d, diy * delta_size_d, 0.0);
                     float f_pt = coot::util::density_at_point(xmap, p);
                     float f_av_profile = 10.0f * static_cast<float>(offset_values_map.at(iix*delta_size_i));
                     float d = f_pt - f_av_profile;
                     sum_d += d;
                     float ad = fabsf(d);
                     s += ad;
                  }
               }
               std::cout << "pos " << co.x() << " " << co.y() << " " << co.z() << " score " << s << " sum_d: " << sum_d << std::endl;
               tafd[idx] = s;
            }
         }
      }
      float z_pos = cell.c() * 2.0;
      texture_as_floats_t texture_as_floats(nu, nv, tafd, x_size, y_size, z_pos);
      std::vector<double> tafdd(tafd.size());
      for (unsigned int i=0; i< tafd.size(); i++) {tafdd[i] = tafd[i]; }
      coot::stats::single s(tafdd);
      float value_for_bottom = s.mean() - 2.0 * std::sqrt(s.variance());
      float value_for_top    = s.mean() + 3.0 * std::sqrt(s.variance());
      mini_texture_t mini_texture(texture_as_floats, value_for_bottom, value_for_top);
      return mini_texture;
   };

   auto orth_to_grid = [] (const clipper::Coord_orth &pos, const clipper::Xmap<float> &xmap) {
      clipper::Coord_frac cf = pos.coord_frac(xmap.cell());
      clipper::Coord_grid cg = cf.coord_grid(xmap.grid_sampling());
      return cg;
   };

   auto get_profile = [orth_to_grid] (const clipper::Coord_orth &co, int section_index, const clipper::Xmap<float> &xmap,
                                     int n_pixels_per_half_edge) {

      clipper::Coord_frac a_cf = co.coord_frac(xmap.cell());
      clipper::Coord_map  a_cm = a_cf.coord_map(xmap.grid_sampling());
      // clipper::Interp_nearest::interp(xmap, a_cm, dv);

      int n_pixels = 2 * n_pixels_per_half_edge + 1;
      profile_t profile(n_pixels_per_half_edge, n_pixels, n_pixels);
      clipper::Coord_grid cg_centre = orth_to_grid(co, xmap);
      clipper::Coord_grid cg_0 = cg_centre + clipper::Coord_grid(-n_pixels_per_half_edge, -n_pixels_per_half_edge, 0);
      clipper::Coord_grid cg_1 = cg_centre + clipper::Coord_grid( n_pixels_per_half_edge,  n_pixels_per_half_edge, 0);
      clipper::Grid_map grid(cg_0, cg_1);
      clipper::Xmap_base::Map_reference_coord ix( xmap, grid.min()), iu, iv, iw;
      int c_u = 0;
      for (iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u()) {
         int c_v = 0;
         for (iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v()) {
            int c_w = 0;
            for (iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w()) {
               const float &f = xmap[iw];
               profile.set_value(c_u, c_v, f);
               c_w++;
            }
            c_v++;
         }
         c_u++;
      }
      return profile;
   };

   auto generate_average_profile = [get_profile, print_profile] (const std::vector<clipper::Coord_orth> &positions,
                                                 const clipper::Xmap<float> &xmap,
                                                 int section_index,
                                                 int n_pixels_per_half_edge) {
      int n_pixels = 2 * n_pixels_per_half_edge + 1;
      profile_t av_prof(n_pixels_per_half_edge, n_pixels, n_pixels);
      profile_t sum_profile(n_pixels_per_half_edge, n_pixels, n_pixels);
      for (unsigned int i=0; i<positions.size(); i++) {
         const auto &p = positions[i];
         auto patch_profile = get_profile(p, section_index, xmap, n_pixels_per_half_edge);
         sum_profile.add_profile(patch_profile);
         if (false) {
            std::cout << "----------- pixed point profile " << i << std::endl;
            print_profile(patch_profile, false);
         }
      }

      // std::cout << "---------------- sum profile " << std::endl;
      // print_profile(sum_profile, false);

      // now fill the values in av_profile
      int n_pixels_sq = av_prof.mean_and_var.size();
      int n_position = positions.size();
      // std::cout << "####### compare n_pixels: " << n_pixels << " " << n_pixels_sq << std::endl;
      for (int i=0; i<n_pixels_sq; i++) {
         av_prof.mean_and_var[i].first = sum_profile.mean_and_var[i].first / static_cast<float>(n_position);
      }
      return av_prof;
   };

   auto score_each_pixel = [] (const profile_t &av_profile, const clipper::Xmap<float> &xmap, int section_index) {

      int hw = av_profile.n_pixels_per_half_edge;
      int n_pixels = hw * 2 + 1;
      clipper::Grid_sampling gs = xmap.grid_sampling();
      texture_as_floats_t results_texture(gs.nu(), gs.nv(), xmap.cell().a(), xmap.cell().b());
      clipper::Coord_grid cg_0(0,0,section_index);
      clipper::Coord_grid cg_1(gs.nu()-1, gs.nv()-1, section_index);
      clipper::Grid_map grid(cg_0, cg_1);
      clipper::Xmap_base::Map_reference_coord ix(xmap, grid.min()), iu, iv, iw;

      int c_u = 0;
      for (iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
         // if (c_u >= 300) continue;
         int c_v = 0;
         for (iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
            for (iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {

               clipper::Coord_grid cgi_0(iw.coord().u()-hw, iw.coord().v()-hw, iw.coord().w());
               clipper::Coord_grid cgi_1(iw.coord().u()+hw, iw.coord().v()+hw, iw.coord().w());
               clipper::Grid_map grid_inner(cgi_0, cgi_1);

               float delta_sum = 0.0;
               clipper::Xmap_base::Map_reference_coord iix(xmap, grid_inner.min()), iiu, iiv, iiw;
               int x_profile = -av_profile.n_pixels_per_half_edge;
               for (iiu = iix; iiu.coord().u() <= grid_inner.max().u(); iiu.next_u() ) {
                  int y_profile = -av_profile.n_pixels_per_half_edge;
                  for (iiv = iiu; iiv.coord().v() <= grid_inner.max().v(); iiv.next_v() ) {
                     for (iiw = iiv; iiw.coord().w() <= grid_inner.max().w(); iiw.next_w() ) {
                        float f = xmap[iiw];
                        if (false)
                           std::cout << "inner f " << f
                                     << " x_profile " << x_profile << " y_profile " << y_profile
                                     << " profile_value " << av_profile.get(x_profile, y_profile) << std::endl;
                        float delta = f-av_profile.get(x_profile, y_profile);
                        delta_sum += fabsf(delta);
                     }
                     y_profile++;
                  }
                  x_profile++;
               }
               int idx = c_v * gs.nu() + c_u;
               results_texture.image_data[idx] = delta_sum;
               if (false)
                  std::cout << "set-image-data "  << " idx " << idx << " size: " << results_texture.image_data.size()
                            << " c_u " << c_u << " c_v " << c_v << " delta_sum " << delta_sum << std::endl;

            }
            c_v++;
         }
         std::cout << "done strip c_u " << c_u << std::endl;
         c_u++;
      }
      return results_texture;
   };

   // ------------------------------------------------------------------------------------------------------
   // ------------------------------------------------------------------------------------------------------

   std::vector<clipper::Coord_orth> positions;
   if (PyList_Check(spot_positions)) {
      long l = PyObject_Length(spot_positions);
      for (unsigned int i=0; i<l; i++) {
         PyObject *item = PyList_GetItem(spot_positions, i);
         if (PyDict_Check(item)) {
            // std::cout << "found a dict " << item << std::endl;
            PyObject *pos_item = PyDict_GetItemString(item, "position"); // a borrowed reference
            if (PyList_Check(pos_item)) {
               long lpi = PyObject_Length(pos_item);
               if (lpi == 3) {
                  double x = PyFloat_AsDouble(PyList_GetItem(pos_item, 0));
                  double y = PyFloat_AsDouble(PyList_GetItem(pos_item, 1));
                  double z = PyFloat_AsDouble(PyList_GetItem(pos_item, 2));
                  clipper::Coord_orth p(x,y,z);
                  positions.push_back(p);
               }
            }
         }
      }
   }

   if (positions.empty()) {
      std::cout << "no positions " << std::endl;
   } else {

      if (false) { // checking that the values of the averaged histogram match those of the image
                   // they do no.
         int section_index = 62;
         graphics_info_t g;
         if (g.is_valid_map_molecule(imol_map)) {
            const auto &xmap = g.molecules[imol_map].xmap;
            clipper::Grid_sampling gs = xmap.grid_sampling();
            clipper::Coord_grid cg_0(0,0,section_index);
            clipper::Coord_grid cg_1(gs.nu()-1, gs.nv()-1, section_index);
            clipper::Grid_map grid(cg_0, cg_1);
            clipper::Xmap_base::Map_reference_coord ix( xmap, grid.min()), iu, iv, iw;

            int count = 0;
            for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
               for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
                  for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
                     clipper::Coord_grid cg = iw.coord();
                     clipper::Coord_frac cf = cg.coord_frac(xmap.grid_sampling());
                     clipper::Coord_orth co = cf.coord_orth(xmap.cell());
                     const float &f = xmap[iw];
                     count++;
                     if (count == 100) {
                        std::cout << "selected xyz " << co.x() << " " << co.y() << " " << co.z()
                                  << " " << f << std::endl;
                        count = 0;
                     }
                  }
               }
            }
         }
      }

      if (false) {
         graphics_info_t g;
         if (g.is_valid_map_molecule(imol_map)) {
            const auto &xmap = g.molecules[imol_map].xmap;
            for (unsigned int i=0; i<positions.size(); i++) {
               const auto &position = positions[i];
               std::cout << "pos " << i << " " << position.format() << std::endl;
               do_stats(i, position, xmap);
            }

            std::string averaged_spot_tab = "averaged-spots.tab";
            if (coot::file_exists(averaged_spot_tab)) {
               std::ifstream f(averaged_spot_tab);
               if (f) {
                  std::map<int, double> offset_values_map;
                  std::string line;
                  while(std::getline(f, line)) {
                     std::vector<float> row;
                     std::stringstream ss(line);
                     int x_offset = -9999;
                     double val = -9999;
                     ss >> x_offset;
                     ss >> val;
                     offset_values_map[x_offset] = val;
                  }
                  int section_index = 62;
                  mini_texture_t mt = map_scan(offset_values_map, xmap, section_index);
                  Texture t(mt, "mini-texture scores");
                  TextureMesh tm("min-texture scores texture-mesh");
                  TextureInfoType ti(t, "mini-texture scores info-type");
                  ti.unit = 0;
                  tm.add_texture(ti);
                  g.texture_meshes.push_back(tm);
               }
            } else {
               std::cout << "WARNING:: missing files " << averaged_spot_tab << std::endl;
            }
         }
      }

      if (true) {
         graphics_info_t g;
         const auto &xmap = g.molecules[imol_map].xmap;
         g.attach_buffers();
         // generate average profile from positions
         int section_index = 62;
         int n_pixels_per_half_edge = 8;
         profile_t av_profile = generate_average_profile(positions, xmap, section_index, n_pixels_per_half_edge);
         std::cout << "------------------ average profile ----- " << std::endl;
         print_profile(av_profile, false);
         std::cout << "------------------ average profile ----- " << std::endl;
         print_profile(av_profile, true);
         texture_as_floats_t results_texture = score_each_pixel(av_profile, xmap, section_index);
         // these 7 lines, converting from floats texture to mini-texture are common
         std::vector<double> tafdd(results_texture.image_data.size());
         for (unsigned int i=0; i< results_texture.image_data.size(); i++) {tafdd[i] = results_texture.image_data[i]; }
         coot::stats::single s(tafdd);
         float value_for_bottom = 0.0; // minimum possible value
         float value_for_top    = s.mean() + 1.0 * std::sqrt(s.variance());
         mini_texture_t mt(results_texture, value_for_bottom, value_for_top);
         std::cout << "mini-texture mt " << mt.x_size << " " << mt.y_size
                   << " mean " << s.mean() << " sd: " << std::sqrt(s.variance())
                   << " data range: " << value_for_bottom << " " << value_for_top << std::endl;
         //
         Texture t(mt, "mini-texture pixel scores");
         TextureMesh tm("min-texture pixel scores texture-mesh");
         TextureInfoType ti(t, "mini-texture pixel scores info-type");
         ti.unit = 0;
         tm.add_texture(ti);
         float z_pos = 5000.0;
         tm.setup_tomo_quad(mt.x_size, mt.y_size, 0, 0, z_pos, true);
         g.texture_meshes.push_back(tm);
         mt.clear();
         g.graphics_draw();
         std::cout << "tomo map analysis done!" << std::endl;
      }
   }
   std::cout << "returning from tomo_map_analysis() " << std::endl;

}
#endif


#ifdef STANDALONE
int main(int argc, char **argv) {
   int status = 0;
   if (argc > 1) {
      std::string map_file_name(argv[1]);
      clipper::CCP4MAPfile file;
      clipper::Xmap<float> xmap;
      if (coot::file_exists(map_file_name)) {
         file.open_read(map_file_name);
         file.import_xmap(xmap);
         file.close_read();
         mini_texture_t(xmap, 128, 0, 1);
      }
   }
   return status;
}
#endif
