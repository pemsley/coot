/*
 * geometry/srs-interface.cc
 *
 * Copyright 2016 by Medical Research Council
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
#include "srs-interface.hh"
#include <stdlib.h>

#include "utils/coot-utils.hh"

// use environment variables or ccp4/prefix-dir fall-back
std::string
coot::get_srs_dir() {

   std::string dir;
   const char *d1 = getenv(MONOMER_DIR_STR); // "COOT_CCP4SRS_DIR"
   const char *d2 = getenv("CCP4");

   if (d1) {
      if (file_exists(d1))
	 dir = d1;
   } else {
      if (d2) {
	 std::string dir_a = util::append_dir_dir(d2, "share");
	 std::string dir_b = util::append_dir_dir(dir_a, "ccp4srs");
	 if (file_exists(dir_b))
	    dir = dir_b;
      }
   }
   
   if (dir.length()) {
      std::cout << "INFO:: CCP4SRS::loadIndex from dir: " << dir << std::endl;
   }

   return dir;
}
