/* ideal/test-indexing.cc
 * 
 * Copyright 2004 The University of York
 * Author: Paul Emsley
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


#include <iostream>
#include <string>
#include <vector>

#include <sys/types.h> // for stating
#include <sys/stat.h>
#include <unistd.h>

#include <string.h>
#include <math.h>

#ifndef  __MMDB_MMCIF__
#include "mmdb_mmcif.h"
#endif
  

#include "mmdb_manager.h"
#include "mmdb-extras.h"
#include "mmdb.h"

#include "test-indexing.hh"


int
main(int argc, char **argv) {

   std::map <std::string, int > atom_map;
   atom_map[" CA "] = 1;
   atom_map[" N  "] = 2;
   
   std::map<std::string, int>::iterator cur = atom_map.find(" CA ");
   std::cout << (*cur).first << std::endl;

   // std::cout << "find: " << atom_map.find(" CA ") << " "
   // 	     << atom_map.find("not")
   //      << std::endl;
   
   coot::testclass c;
   int iresno = 68;
   
   c.add_residue_atom_map(iresno,atom_map);
   c.add_atom(67, " O  ", 50);

//    std::cout << c.atom_name_resno_to_index[iresno][" CA "] << std::endl;
//    std::cout << c.atom_name_resno_to_index[iresno][" N  "] << std::endl;

   // std::cout << c(iresno," CA ") << std::endl;
   // std::cout << c(iresno," N  ") << std::endl;

   coot::testclass d = c;

   std::cout << d.atom_name_resno_to_index[iresno  ][" CA "] << std::endl;
   std::cout << d.atom_name_resno_to_index[iresno-1][" CA "] << std::endl; 
   std::cout << d.atom_name_resno_to_index[iresno-1][" O  "] << std::endl;

   // d.set_big_index("A", iresno, " CB ",43);
   // std:: cout << d.big_index["A"][iresno][" CB "] << std::endl;

   return 0; 
}

