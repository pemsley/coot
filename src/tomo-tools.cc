

#include <clipper/ccp4/ccp4_map_io.h>
#include "utils/coot-utils.hh"
#include "coot-utils/coot-map-heavy.hh"
#include "mini-texture.hh"
#include "graphics-info.h"

void tomo_section(int imol, int section_index, int axis_id);

void tomo_scale_adjustment_changed(GtkAdjustment *adj, gpointer user_data) {

   double v = gtk_adjustment_get_value(adj);
   int section = static_cast<int>(v);
   tomo_section(0, section, 0);

}

void setup_tomo_slider(int imol) {

   graphics_info_t g;
   if (g.is_valid_map_molecule(imol)) {

      const auto &xmap = g.molecules[imol].xmap;

      // move this into graphics_info_t
      clipper::Grid_sampling gs = xmap.grid_sampling();
      int n_sections = gs.nw();
      int mid_section = n_sections/2;

      GtkWidget *scale = widget_from_builder("tomo_scale");
      GtkAdjustment *adjustment_current = gtk_range_get_adjustment(GTK_RANGE(scale));
      if (adjustment_current) {
         gtk_adjustment_set_lower(adjustment_current, 0);
         gtk_adjustment_set_upper(adjustment_current, n_sections-1);
         gtk_adjustment_set_step_increment(adjustment_current, 1.0);
         gtk_adjustment_set_page_increment(adjustment_current, 5.0);
         gtk_adjustment_set_page_size(adjustment_current, 0.0);
         gtk_adjustment_set_value(adjustment_current, mid_section);
         gtk_scale_set_draw_value(GTK_SCALE(scale), TRUE);
         gtk_scale_set_digits(GTK_SCALE(scale), 0);
      }
   }
}

void tomo_section(int imol, int section_index, int axis_id) {

   graphics_info_t g;
   if (g.is_valid_map_molecule(imol)) {
      g.tomo_section(imol, section_index);
   }

}

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
         mini_texture_t(xmap, 128);
      }
   }
   return status;
}
#endif
