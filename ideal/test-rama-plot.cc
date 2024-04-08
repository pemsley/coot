/*
 * ideal/test-rama-plot.cc
 * 
 * Copyright 2018 by Medical Research Council
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
 * General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 */


#include <iostream>
#include "zo-rama.hh"

int main(int argc, char **argv) {

   zo::rama_table_set rts;
   std::map<std::string, zo::rama_table>::const_iterator it;
   for (it=rts.table_map.begin(); it!=rts.table_map.end(); it++) {
      const zo::rama_table &rt = it->second;
      std::string zo_residue_type = it->first;
      std::cout << zo_residue_type << " " << rt.rama_vec.size() << std::endl;
      std::string file_name = zo_residue_type + ".png";
      rt.make_a_png(600, file_name);
   }
   return 0;

}
