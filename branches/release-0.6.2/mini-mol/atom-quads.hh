/* mini-mol/atom-quads.hh
 * 
 * Copyright  2009 The University of Oxford
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


#ifndef LIGAND_ATOM_QUAD_HH
#define LIGAND_ATOM_QUAD_HH

#ifndef HAVE_STRING
#include <string>
#define HAVE_STRING
#endif

#include <mmdb_manager.h>

namespace coot { 

   // 4 atom names for a torsion:
   class atom_name_quad {
   public:
      std::string atom1;
      std::string atom2;
      std::string atom3;
      std::string atom4;
      atom_name_quad() {} // sigh.  needed for torsioned_atoms_info_t,
			  // don't know why.
      atom_name_quad(const std::string &atom_name_1,
                     const std::string &atom_name_2,
                     const std::string &atom_name_3,
                     const std::string &atom_name_4) {
         atom1 = atom_name_1;
         atom2 = atom_name_2;
         atom3 = atom_name_3;
         atom4 = atom_name_4;
      }
      bool all_non_blank() const {
	 if ((atom1 == "") || (atom2 == "") || (atom3 == "") || (atom4 == ""))
	    return 0;
	 else
	    return 1;
      } 
      friend std::ostream& operator<<(std::ostream &o, const atom_name_quad &q);
   };
   std::ostream& operator<<(std::ostream &o, const atom_name_quad &q);

   class atom_index_quad {
   public:
      int index1;
      int index2;
      int index3;
      int index4;
      atom_index_quad() { // unassigned pair
         index1 = -1;
         index2 = -1;
         index3 = -1;
         index4 = -1;
      }
      atom_index_quad(int i1, int i2, int i3, int i4) {
         index1 = i1;
         index2 = i2;
         index3 = i3;
         index4 = i4;
      }

      // Return the torsion in degrees.  Use the indices to index into
      // residue res to find the atoms for the torsion.
      // 
      // Can throw an exception on not all indices found
      // 
      double torsion(CResidue *res) const;
      friend std::ostream& operator<<(std::ostream &o, const atom_index_quad &q);
   };
   std::ostream& operator<<(std::ostream &o, const atom_index_quad &q);

}

#endif // LIGAND_ATOM_QUAD_HH

