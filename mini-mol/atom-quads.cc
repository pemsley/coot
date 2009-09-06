/* mini-mol/atom-quads.cc
 * 
 * Copyright  2009 The University of York
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


#include <stdexcept>
#include "clipper/core/coords.h"

#include "atom-quads.hh"

std::ostream&
coot::operator<<(std::ostream &o, const coot::atom_name_quad &q) {

   o << "(" << q.atom1 << " " << q.atom2 << " " << q.atom3 << " " << q.atom4 << ")";
   return o;
}


std::ostream&
coot::operator<<(std::ostream &o, const coot::atom_index_quad &q) {

   o << "(" << q.index1 << " " << q.index2 << " " << q.index3 << " " << q.index4 << ")";
   return o;
}

// can throw an exception
double 
coot::atom_index_quad::torsion(CResidue *residue_p) const {

   double angle = 0;
   int n_residues_atoms;
   PPCAtom residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residues_atoms);
   for (unsigned int i=0; i<n_residues_atoms; i++) {
      bool good_indexing = 0;
      if ((index1 >= 0) && (index1 < n_residues_atoms)) { 
	 if ((index2 >= 0) && (index2 < n_residues_atoms)) { 
	    if ((index3 >= 0) && (index3 < n_residues_atoms)) { 
	       if ((index4 >= 0) && (index4 < n_residues_atoms)) {
		  good_indexing = 1;
	       }
	    }
	 }
      }
      if (! good_indexing) {
	 std::string mess = "bad atom indexing in atom_index_quad::torsion()";
	 throw std::runtime_error(mess);
      }
      clipper::Coord_orth pt_1(residue_atoms[index1]->x,
			       residue_atoms[index1]->y,
			       residue_atoms[index1]->z);
      clipper::Coord_orth pt_2(residue_atoms[index2]->x,
			       residue_atoms[index2]->y,
			       residue_atoms[index2]->z);
      clipper::Coord_orth pt_3(residue_atoms[index3]->x,
			       residue_atoms[index3]->y,
			       residue_atoms[index3]->z);
      clipper::Coord_orth pt_4(residue_atoms[index4]->x,
			       residue_atoms[index4]->y,
			       residue_atoms[index4]->z);

      angle = clipper::Util::rad2d(clipper::Coord_orth::torsion(pt_1, pt_2, pt_3, pt_4));
   } 
   return angle;
}
