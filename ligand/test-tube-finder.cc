/*
 * ligand/test-tube-finder.cc
 *
 * Copyright 2021 by Medical Research Council
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

#include <vector>
#include <string>
#include <clipper/core/xmap.h>
#include <clipper/ccp4/ccp4_map_io.h>
#include "tube-finder.hh"

int main(int argc, char **argv) {

   int status = 0;

   if (argc> 1) {
      std::string map_file_name = argv[1];
      clipper::CCP4MAPfile file;
      clipper::Xmap<float> xmap;
      file.open_read(map_file_name);
      file.import_xmap(xmap);
      coot::tube_finder_t tf(xmap);
      std::vector<clipper::Coord_orth> positions = tf.get_positions();
      for (unsigned int i=0; i<positions.size(); i++) {
         std::cout << i << " " << positions[i].format() << std::endl;
      }
   }

   return status;
}
