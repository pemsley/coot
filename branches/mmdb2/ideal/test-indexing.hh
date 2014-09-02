/* ideal/test-indexing.hh
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

#include <string>
#include <vector>
#include <map>


namespace coot {
   
   class testclass {
   public:
      
      // std::map <std::vector<std::map <std::string, int> > > big_index;

      std::vector<std::map <std::string, int> > atom_name_resno_to_index;
      void add_residue_atom_map(int iresno, const std::map<std::string, int> &atom_map) {
	 if (iresno > atom_name_resno_to_index.size() ) {
	    atom_name_resno_to_index.resize(iresno + 1);
	 }
	 atom_name_resno_to_index[iresno] = atom_map;
      }
      void add_atom(int iresno, const std::string &at_name, int atom_index) {
	 if (iresno > atom_name_resno_to_index.size() ) {
	    atom_name_resno_to_index.resize(iresno + 1);
	 }
	 atom_name_resno_to_index[iresno][at_name] = atom_index;
      }

      void set_big_index(const std::string &chain,
			 int iresno,
			 const std::string &at_name, int atom_index) {
	 
      }
   };

}
