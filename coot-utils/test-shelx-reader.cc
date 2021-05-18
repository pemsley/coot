/* coot-utils/test-shelx-reader.cc
 * 
 * Copyright 2005, 2006 by Paul Emsley, The University of York
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


#include <map>

#include "utils/coot-utils.hh"
#include "coot-shelx.hh"
#include <iostream>

int
main(int argc, char **argv) {

   if (argc > 2) {
      coot::ShelxIns sh;
      coot::shelx_read_file_info_t p = sh.read_file(argv[1]);
      sh.write_ins_file(p.mol, std::string(argv[2]));
   } else {
      std::cout << "Usage: " << argv[0] << " shelx-ins-file-name out-file-name"
                << std::endl;
   }
   return 0;
}
