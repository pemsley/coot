/*
 * coot-utils/positron.cc
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

#include <fstream>
#include <iostream>
#include <sstream>

#include "positron.hh"
#include "utils/coot-utils.hh"

void
coot::read_positron_metadata(std::vector<positron_metadata_t> *data, const std::string &z_data, const std::string &s_data) {

   std::vector<positron_metadata_t> &vec(*data);
   vec.reserve(330000);

   // std::cout << "coot::read_positron_metadata() A " << z_data << std::endl;
   if (coot::file_exists(z_data)) {
      // std::cout << "coot::read_positron_metadata() B " << s_data << std::endl;
      if (coot::file_exists(s_data)) {
         // std::cout << "coot::read_positron_metadata() C" << std::endl;
         std::ifstream file_z(z_data);
         std::ifstream file_s(s_data);
         float x, y;
         char c; // to eat the commas
#if 1
         while (file_z >> x >> c >> y) {
            // std::cout << "coot::read_positron_metadata() D " << x << " " << y << std::endl;
            positron_metadata_t pm(x,y);
            file_s >> pm.params[0] >> c >> pm.params[1] >> c >> pm.params[2] >> c
                   >> pm.params[3] >> c >> pm.params[4] >> c >> pm.params[5];
            vec.push_back(pm);
         }
#endif
#if 0
         std::string z_line;
         while (std::getline(file_z, z_line)) {
            std::stringstream ss(z_line);
            float x;
            ss >> x;
            std::cout << "z_line: " << z_line << " x " << x << std::endl;
         }
#endif

      } else {
         std::cout << "Error: File does not exist: " << s_data << std::endl;
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

   if (false)
      std::cout << "HERE we are in get_closest_positron_metadata_point() with metadata size " << positron_metadata.size()
                << " and test position " << x << " "  << y << std::endl;

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
