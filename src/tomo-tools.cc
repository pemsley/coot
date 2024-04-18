
#include <chrono>

#include <clipper/ccp4/ccp4_map_io.h>
#include "utils/coot-utils.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/coot-map-heavy.hh"
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
            std::cout << "map-point " << point_index << " " << ix << " "
                      << dix * delta_size << " " << diy * delta_size << " "
                      << p.x() << " " << p.y() << " " << p.z() << " value: " << f << std::endl;
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

      for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
         for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
            for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
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
            }
         }
      }
   };

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
            int nv = gs.nv();
            int nu = gs.nu();

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
                        std::cout << "selected xyz " << co.x() << " " << co.y() << " " << co.z() << " " << f << std::endl;
                        count = 0;
                     }
                  }
               }
            }
         }
      }

      if (true) {
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
                  map_scan(offset_values_map, xmap, section_index);
               }
            } else {
               std::cout << "WARNING:: missing files " << averaged_spot_tab << std::endl;
            }
         }
      }
   }

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
