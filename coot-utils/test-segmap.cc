/*
 * coot-utils/test-segmap.cc
 *
 * Copyright 2019 by Medical Research Council
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
#include <string>
#include <iostream>
#include "segmap.hh"

int main(int argc, char **argv) {

   bool done = false;

   if (argc == 3) {
      // generate the stats from the sampled maps.
      std::string a1(argv[1]);
      if (a1 == "proc") {
         std::string map_file_name(argv[2]);
         clipper::CCP4MAPfile file;
         try {
            file.open_read(map_file_name.c_str());
            clipper::Xmap<float> xmap;
            file.import_xmap(xmap);
            coot::segmap s(xmap);
            bool do_write_flag = true;
            s.proc(do_write_flag, "segmap-out.map");
            done = true;
         }
         catch (const clipper::Message_base &exc) {
            std::cout << "WARNING:: failed to open " << map_file_name << std::endl;
         }
      }

      if (a1 == "dedust") {
         std::string map_file_name(argv[2]);
         clipper::CCP4MAPfile file;
         try {
            file.open_read(map_file_name.c_str());
            clipper::Xmap<float> xmap;
            file.import_xmap(xmap);
            coot::segmap s(xmap);
            s.dedust();
            done = true;
         }
         catch (const clipper::Message_base &exc) {
            std::cout << "WARNING:: failed to open " << map_file_name << std::endl;
         }
      }
   }

   return 0;
}
