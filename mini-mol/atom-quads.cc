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
coot::atom_quad::atom_quad(mmdb::Residue *first, mmdb::Residue *second, const std::string &link) {

   // Fixup needed for PDBv3.
   //
   coot::atom_quad quad;
   
   // standard (1-linked) atom names for chiral centre:
   std::string atom_2_name = " O5 "; 
   std::string atom_3_name = " C2 "; 
   std::string atom_4_name = " C1 "; // chiral centre (on second residue)
   
   std::string O_name;
   
   if (link == "ALPHA1-2" || link == "BETA1-2" )
      O_name = " O2 "; 
   if (link == "ALPHA1-3" || link == "BETA1-3" )
      O_name = " O3 "; 
   if (link == "ALPHA1-4" || link == "BETA1-4" )
      O_name = " O4 "; 
   if (link == "ALPHA1-6" || link == "BETA1-6" )
      O_name = " O6 ";

   if (! O_name.empty()) {
      // is standard (non-SIA) link
      std::vector<std::string> v;
      v.push_back(atom_2_name);
      v.push_back(atom_3_name);
      v.push_back(atom_4_name); // chiral centre
      quad = setup_chiral_quad(first, second, O_name, v);

   } else { 

      if (link == "ALPHA2-3" || link == "BETA2-3" ) {
	 O_name = " O3 ";      // from first residue
	 atom_2_name = " O6 ";
	 atom_3_name = " C3 ";
	 atom_4_name = " C2 "; // chiral centre (on 1st residue)

	 std::vector<std::string> v;
	 v.push_back(atom_2_name);
	 v.push_back(atom_3_name);
	 v.push_back(atom_4_name); // chiral centre
	 quad = setup_chiral_quad(second, first, O_name, v);
      }

      if (link == "ALPHA2-6") {
	 O_name = " O6 ";
	 atom_2_name = " O6 ";
	 atom_3_name = " C3 ";
	 atom_4_name = " C2 "; // chiral centre (on 1st residue)
	 std::vector<std::string> v;
	 v.push_back(atom_2_name);
	 v.push_back(atom_3_name);
	 v.push_back(atom_4_name); // chiral centre
	 quad = setup_chiral_quad(first, second, O_name, v);
      }
   }

   // all or nothing
   // 
   if (quad.atom_1 == NULL || quad.atom_2 == NULL || quad.atom_3 == NULL || quad.atom_4 == NULL) {

      if (false) {
	 std::cout << "WARNING:: atom_quad() problems given "
		   << first->GetChainID() << " " << first->GetSeqNum() << " " << first->GetInsCode()
		   << "  and "
		   << second->GetChainID() << " " << second->GetSeqNum() << " " << second->GetInsCode()
		   << " and link type " << link << std::endl;
	 if (!quad.atom_1)
	    std::cout << "WARNING:: atom_quad() [for chiral] missing 1 \"" << O_name << "\"" << std::endl;
	 if (!quad.atom_2)
	    std::cout << "WARNING:: atom_quad() [for chiral] missing 2 \"" << atom_2_name << "\"" << std::endl;
	 if (!quad.atom_3)
	    std::cout << "WARNING:: atom_quad() [for chiral] missing 3 \"" << atom_3_name << "\"" << std::endl;
	 if (!quad.atom_4)
	    std::cout << "WARNING:: atom_quad() [for chiral] missing 4 \"" << atom_4_name << "\"" << std::endl;
      }
      atom_1 = NULL;
      atom_2 = NULL;
      atom_3 = NULL;
      atom_4 = NULL;
   } else {
      *this = quad;
      name = link;
   } 
}

bool
coot::atom_quad::filled_p() const { // ! were there any nulls?

   bool filled = true;
   if (atom_1 == NULL) filled = false;
   if (atom_2 == NULL) filled = false;
   if (atom_3 == NULL) filled = false;
   if (atom_4 == NULL) filled = false;
   return filled;
} 



// a bit strange :-)
// 
coot::atom_quad
coot::atom_quad::setup_chiral_quad(mmdb::Residue *residue_with_O, mmdb::Residue *residue_with_chiral_centre,
				   const std::string &O_name,
				   const std::vector<std::string> &chiral_atom_names) const {

   coot::atom_quad quad;

   // yes.
   std::string atom_2_name = chiral_atom_names[0];
   std::string atom_3_name = chiral_atom_names[1];
   std::string atom_4_name = chiral_atom_names[2]; // chiral atom

   if (O_name != "") { 
      mmdb::PPAtom residue_atoms = NULL;
      int n_residue_atoms;
      residue_with_O->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
	 std::string atom_name(residue_atoms[iat]->name);
	 if (atom_name == O_name) {
	    // the O atom name comes from the residue_with_O
	    if (!quad.atom_1) { 
	       quad.atom_1 = residue_atoms[iat];
	       break;
	    }
	 }
      }
      // atoms 2 and 3 and the chiral centre atom come from the residue_with_chiral_centre residue.
      residue_atoms = NULL;
      residue_with_chiral_centre->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
	 mmdb::Atom *at = residue_atoms[iat];
	 std::string atom_name(at->name);
	 if (atom_name == atom_4_name) {
	    if (! quad.atom_4)
	       quad.atom_4 = at; // chiral centre atom
	 }
	 if (atom_name == atom_2_name) {
	    if (! quad.atom_2)
	       quad.atom_2 = at;
	 }
	 if (atom_name == atom_3_name) {
	    if (! quad.atom_3)
	       quad.atom_3 = at;
	 }
      }
   }
   return quad;
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
coot::atom_index_quad::torsion(mmdb::Residue *residue_p) const {

   double angle = 0;
   int n_residues_atoms;
   mmdb::PPAtom residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residues_atoms);
   return torsion(residue_atoms, n_residues_atoms);
}


// as above, but we reference atoms in a selection rather than a
// residue.
// 
double
coot::atom_index_quad::torsion(mmdb::PPAtom atom_selection, int n_selected_atoms) const {

   double angle = 0;
   for (int i=0; i<n_selected_atoms; i++) {
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
coot::atom_quad::angle_2() const {  // angle 1-2-3 in degrees

   if (atom_1 && atom_2 && atom_3) { 
      clipper::Coord_orth pt_1(atom_1->x, atom_1->y, atom_1->z);
      clipper::Coord_orth pt_2(atom_2->x, atom_2->y, atom_2->z);
      clipper::Coord_orth pt_3(atom_3->x, atom_3->y, atom_3->z);
      double angle = clipper::Util::rad2d(clipper::Coord_orth::angle(pt_1, pt_2, pt_3));
      return angle;
   } else {
      throw std::runtime_error("quad::torsion() Null atom(s)");
   }
} 


double
coot::atom_quad::angle_3() const { // angle 2-3-4 in degrees

   if (atom_2 && atom_3 && atom_4) { 
      clipper::Coord_orth pt_2(atom_2->x, atom_2->y, atom_2->z);
      clipper::Coord_orth pt_3(atom_3->x, atom_3->y, atom_3->z);
      clipper::Coord_orth pt_4(atom_4->x, atom_4->y, atom_4->z);
      double angle = clipper::Util::rad2d(clipper::Coord_orth::angle(pt_2, pt_3, pt_4));
      return angle;
   } else {
      throw std::runtime_error("quad::torsion() Null atom(s)");
   }
}

coot::atom_name_quad
coot::atom_quad::get_atom_name_quad() const {

   if (atom_1 && atom_2 && atom_3 && atom_4) {
      return coot::atom_name_quad(atom_1->name,
				  atom_2->name,
				  atom_3->name,
				  atom_4->name);
   } else {
      throw std::runtime_error("atom_quad::atom_name_quad() Null atom(s)");
   }
}


double
coot::atom_name_quad::torsion(mmdb::Residue *residue) const {

   double r = -999.9;

   mmdb::Atom *at_0 = residue->GetAtom(atom_name_[0].c_str());
   mmdb::Atom *at_1 = residue->GetAtom(atom_name_[1].c_str());
   mmdb::Atom *at_2 = residue->GetAtom(atom_name_[2].c_str());
   mmdb::Atom *at_3 = residue->GetAtom(atom_name_[3].c_str());

   if (at_0 && at_1 && at_2 && at_3) {
      clipper::Coord_orth pt_0(at_0->x, at_0->y, at_0->z);
      clipper::Coord_orth pt_1(at_1->x, at_1->y, at_1->z);
      clipper::Coord_orth pt_2(at_2->x, at_2->y, at_2->z);
      clipper::Coord_orth pt_3(at_3->x, at_3->y, at_3->z);
      double angle = clipper::Util::rad2d(clipper::Coord_orth::torsion(pt_0, pt_1, pt_2, pt_3));
      return angle;
   } 
   return r;
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
