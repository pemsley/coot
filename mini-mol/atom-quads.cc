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

// return the cente atom as atom_4, the others as numbered.
// If not all atoms found, then return all nulls.
//
// Just because I am lazy, setup wth the first atom names that match -
// don't consider alt confs.
// 
coot::atom_quad::atom_quad(CResidue *first, CResidue *second, const std::string &link) {

   // Fixup needed for PDBv3.
   //

   atom_1 = NULL;
   atom_2 = NULL;
   atom_3 = NULL;
   atom_4 = NULL; // atom_c for a chiral quad.
   
   std::string O_name = "";
   if (link == "ALPHA1-2" || link == "BETA1-2" )
      O_name = " O2 "; 
   if (link == "ALPHA1-3" || link == "BETA1-3" )
      O_name = " O3 "; 
   if (link == "ALPHA1-4" || link == "BETA1-4" )
      O_name = " O4 "; 
   if (link == "ALPHA1-6" || link == "BETA1-6" )
      O_name = " O6 ";
   if (O_name != "") { 
      PPCAtom residue_atoms = 0;
      int n_residue_atoms;
      first->GetAtomTable(residue_atoms, n_residue_atoms);
      for (unsigned int iat=0; iat<n_residue_atoms; iat++) {
	 std::string atom_name(residue_atoms[iat]->name);
	 if (atom_name == O_name) {
	    // the O atom name comes from the first residue
	    if (!atom_1) { 
	       atom_1 = residue_atoms[iat];
	       break;
	    }
	 }
      }
      // atoms 2 and 3 and the chiral centre atom come from the second residue.
      residue_atoms = NULL;
      second->GetAtomTable(residue_atoms, n_residue_atoms);
      for (unsigned int iat=0; iat<n_residue_atoms; iat++) {
	 CAtom *at = residue_atoms[iat];
	 std::string atom_name(at->name);
	 if (atom_name == " C1 ") {
	    if (! atom_4)
	       atom_4 = at; // chiral centre atom
	 }
	 if (atom_name == " O5 ") {
	    if (! atom_2)
	       atom_2 = at;
	 }
	 if (atom_name == " C2 ") {
	    if (! atom_3)
	       atom_3 = at;
	 }
      }
   }

   // all or nothing
   // 
   if (atom_1 == NULL || atom_2 == NULL || atom_3 == NULL || atom_4 == NULL) {
      atom_1 = NULL;
      atom_2 = NULL;
      atom_3 = NULL;
      atom_4 = NULL;
   } 
} 


std::ostream&
coot::operator<<(std::ostream &o, const coot::atom_name_quad &q) {

   o << "(" << q.atom_name(0) << " " << q.atom_name(1) << " " << q.atom_name(2) << " " << q.atom_name(3) << ")";
   return o;
}


std::ostream&
coot::operator<<(std::ostream &o, const coot::atom_index_quad &q) {

   o << "(" << q.index1 << " " << q.index2 << " " << q.index3 << " " << q.index4 << ")";
   return o;
}



std::ostream&
coot::operator<<(std::ostream &o, const coot::atom_quad &q) {

   // can't use specs at this level.  Maybe consider moving atom spec definition then...
   o << "("
     << q.atom_1 << " " << q.atom_2 << " "
     << q.atom_3 << " " << q.atom_4 << ")";
   return o;
}


// can throw an exception
double 
coot::atom_index_quad::torsion(CResidue *residue_p) const {

   double angle = 0;
   int n_residues_atoms;
   PPCAtom residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residues_atoms);
   return torsion(residue_atoms, n_residues_atoms);
}


// as above, but we reference atoms in a selection rather than a
// residue.
// 
double
coot::atom_index_quad::torsion(PPCAtom atom_selection, int n_selected_atoms) const {

   double angle = 0;
   for (unsigned int i=0; i<n_selected_atoms; i++) {
      bool good_indexing = 0;
      if ((index1 >= 0) && (index1 < n_selected_atoms)) { 
	 if ((index2 >= 0) && (index2 < n_selected_atoms)) { 
	    if ((index3 >= 0) && (index3 < n_selected_atoms)) { 
	       if ((index4 >= 0) && (index4 < n_selected_atoms)) {
		  good_indexing = 1;
	       }
	    }
	 }
      }
      if (! good_indexing) {
	 std::string mess = "bad atom indexing in atom_index_quad::torsion()";
	 throw std::runtime_error(mess);
      }
      clipper::Coord_orth pt_1(atom_selection[index1]->x,
			       atom_selection[index1]->y,
			       atom_selection[index1]->z);
      clipper::Coord_orth pt_2(atom_selection[index2]->x,
			       atom_selection[index2]->y,
			       atom_selection[index2]->z);
      clipper::Coord_orth pt_3(atom_selection[index3]->x,
			       atom_selection[index3]->y,
			       atom_selection[index3]->z);
      clipper::Coord_orth pt_4(atom_selection[index4]->x,
			       atom_selection[index4]->y,
			       atom_selection[index4]->z);

      angle = clipper::Util::rad2d(clipper::Coord_orth::torsion(pt_1, pt_2, pt_3, pt_4));
   } 
   return angle;
}

// Can throw a std::runtime_error if any of the atoms are null.
double
coot::atom_quad::torsion() const {

   if (atom_1 && atom_2 && atom_3 && atom_4) { 
      clipper::Coord_orth pt_1(atom_1->x, atom_1->y, atom_1->z);
      clipper::Coord_orth pt_2(atom_2->x, atom_2->y, atom_2->z);
      clipper::Coord_orth pt_3(atom_3->x, atom_3->y, atom_3->z);
      clipper::Coord_orth pt_4(atom_4->x, atom_4->y, atom_4->z);
      double angle = clipper::Util::rad2d(clipper::Coord_orth::torsion(pt_1, pt_2, pt_3, pt_4));
      return angle;
   } else {
      throw std::runtime_error("quad::torsion() Null atom(s)");
   }
} 

// can throw a std::runtime_error exception
// 
double
coot::atom_quad::chiral_volume() const {

   if (atom_1 == NULL || atom_2 == NULL || atom_3 == NULL || atom_4 == NULL) {
      throw std::runtime_error("Null atoms in quad for chiral volume");
   } else { 

      clipper::Coord_orth centre(atom_4->x, atom_4->y, atom_4->z);
      clipper::Coord_orth   at_1(atom_1->x, atom_1->y, atom_1->z);
      clipper::Coord_orth   at_2(atom_2->x, atom_2->y, atom_2->z);
      clipper::Coord_orth   at_3(atom_3->x, atom_3->y, atom_3->z);
      
      clipper::Coord_orth a = at_1 - centre;
      clipper::Coord_orth b = at_2 - centre;
      clipper::Coord_orth c = at_3 - centre;
      double cv = clipper::Coord_orth::dot(a, clipper::Coord_orth::cross(b,c));
      return cv;
   }
} 
