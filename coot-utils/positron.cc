
#include <fstream>
#include <iostream>

#include "positron.hh"
#include "utils/coot-utils.hh"

void
coot::read_positron_metadata(std::vector<positron_metadata_t> *data, const std::string &z_data, const std::string &s_data) {

   std::vector<positron_metadata_t> &vec(*data);
   vec.reserve(330000);

   if (coot::file_exists(z_data)) {
      if (coot::file_exists(s_data)) {
         std::ifstream file_z(z_data);
         std::ifstream file_s(s_data);
         float x, y;
         char c; // to eat the commas
         while (file_z >> x >> c >> y) {
            positron_metadata_t pm(x,y);
            file_s >> pm.params[0] >> c >> pm.params[1] >> c >> pm.params[2] >> c>> pm.params[3] >> c>> pm.params[4] >> c >> pm.params[5];
            vec.push_back(pm);
         }
      }
   }

   // std::cout << "vec.size() " << vec.size() << std::endl;

}

#include <cmath>

// 20240323-PE Should we also pass a "must-be-closer-than-limit?"
// Currently we use 0.1 in both x and y.
// If we fail to find a close point, then return -1.
int
coot::get_closest_positron_metadata_point(const std::vector<positron_metadata_t> &positron_metadata, float x, float y) {

   int idx = -1;
   float lim = 0.2;
   float closest_sqrd = lim;

   // std::cout << "HERE we are in get_closest_positron_metadata_point() with metadata size " << positron_metadata.size()
   // << " and test position " << x << " "  << y << std::endl;

   for (unsigned int i=0; i<positron_metadata.size(); i++) {
      const auto &md(positron_metadata[i]);
      if (std::fabs(md.x - x) < lim) {
         if (std::fabs(md.y - y) < lim) {
            float cs = (md.x - x) * (md.x - x) + (md.y - y) * (md.y - y);
            if (cs < closest_sqrd) {
               closest_sqrd = cs;
               idx = i;
               // std::cout << "new closest to x " << x << " " << y << " " << " at " << md.x << " " << md.y << std::endl;
            }
         }
      }
   }
   return idx;

}
