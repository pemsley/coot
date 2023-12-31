
#include "diff-diff-map-peaks.hh"
#include "clipper/ccp4/ccp4_map_io.h"

int main(int argc, char **argv) {

   clipper::Xmap<float> m1;
   clipper::Xmap<float> m2;
   if (argc > 2) {
      std::string map_file_name_1(argv[1]);
      std::string map_file_name_2(argv[2]);
      try {
         clipper::CCP4MAPfile file;
         file.open_read(map_file_name_1);
         file.import_xmap(m1);
      }
      catch (const clipper::Message_base &exc) {
         std::cout << "WARNING:: failed to open " << map_file_name_1 << std::endl;
      }
      try {
         clipper::CCP4MAPfile file;
         file.open_read(map_file_name_2);
         file.import_xmap(m2);
      }
      catch (const clipper::Message_base &exc) {
         std::cout << "WARNING:: failed to open " << map_file_name_2 << std::endl;
      }

      clipper::Coord_orth screen_centre(10,11,12);
      float base_level = 0.2; //
      float radius = 20;
      std::vector<std::pair<clipper::Coord_orth, float> > map_peaks =
         coot::diff_diff_map_peaks(m1, m2, base_level, screen_centre, radius);

      std::cout << "Found " << map_peaks.size() << " peaks" << std::endl;
      for (size_t i = 0; i < map_peaks.size(); i++) {
         const auto &v   = map_peaks[i].second;
         const auto &pos = map_peaks[i].first;
         std::cout << "  " << i << " " << pos.format() << " " << v << std::endl;
      }
      std::vector<std::pair<clipper::Coord_orth, float> > peaks =
      coot::move_peaks_to_around_position(screen_centre, m1.spacegroup(), m1.cell(), map_peaks);
      for (size_t i = 0; i < peaks.size(); i++) {
         const auto &v   = peaks[i].second;
         const auto &pos = peaks[i].first;
         std::cout << "  " << i << " " << pos.format() << " " << v << std::endl;
      }
   }

   return 0;
}
