/* src/rotamer.cc
 * 
 * Copyright 2001, 2002, 2003, 2004, 2006 The University of York
 * Author: Paul Emsley
 * Copyright 2007 The University of Oxford
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
#include <fstream>
#include <stdexcept>

#include "compat/coot-sysdep.h"

#include "rotamer.hh"
#include "ccp4mg-utils/mgtree.h"

#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-extras.hh"


std::ostream &
coot::operator<<(std::ostream &s, const coot::rotamer_probability_info_t &rpi) {

   s << "[state: " << rpi.state << " name: \"" << rpi.rotamer_name << "\" prob: "
     << rpi.probability << "%]";

   return s;
} 


// Return a flag saying whether we did this or not.
//
// LEU, VAL, THR have "nomenclature" or real chiral centres - they are
// not dealt with here.
// 
// We deal with bifurcated symmetric non-chiral side chains (PHE, ASP,
// GLU, TYR)
// 
int 
coot::rotamer::optimize_rotamer_by_atom_names(bool apply_swap_if_found) {

   // So, the plan is to compare the rotamer probabilies of
   // 1: the input residue
   // 2: copy the input residue and swap the "rotamer" atoms of it.
   //
   // Now compare the probabilites of each rotamer
   //
   // If the second probability is greater than the first, then swap
   // the atom names.
   

   // Chi2 of PHE, TYR should be between -90 -> +90

   int imoved = 0;
   // p.first        // flag for assigned,                    1 
                     // unassigned due to missing atoms,      0 
                     // unassigned due to rotamer not found. -1
                     // unassigned due to GLY/ALA            -2
   coot::rotamer_probability_info_t p_init = probability_of_this_rotamer(); 

   std::string residue_name(residue->GetResName());
   if ((residue_name == "PHE") ||
       (residue_name == "ASP") ||
       (residue_name == "GLU") ||
       (residue_name == "TYR") ) {

      std::vector<std::pair<std::string, std::string> > swapper_atoms;

      if (residue_name == "PHE") {
	 swapper_atoms.push_back(std::pair<std::string, std::string>(" CD1", " CD2"));
	 swapper_atoms.push_back(std::pair<std::string, std::string>(" CE1", " CE2"));
      }
      if (residue_name == "ASP") {
	 swapper_atoms.push_back(std::pair<std::string, std::string>(" OD1", " OD2"));
      }
      if (residue_name == "GLU") {
	 swapper_atoms.push_back(std::pair<std::string, std::string>(" OE1", " OE2"));
      }
      if (residue_name == "TYR") {
	 swapper_atoms.push_back(std::pair<std::string, std::string>(" CD1", " CD2"));
	 swapper_atoms.push_back(std::pair<std::string, std::string>(" CE1", " CE2"));
      }

      // set up a copy of this residue
      // 
      mmdb::Residue *residue_copy = deep_copy_residue(residue); // Yes, gets deleted.
      mmdb::PPAtom residue_atoms;
      int n_residue_atoms;
      residue_copy->GetAtomTable(residue_atoms, n_residue_atoms);
      int nfound = 0; 
      for (unsigned int iswap=0; iswap<swapper_atoms.size(); iswap++) {
	 short int ifound = 0;
	 for (int ifirst=0; ifirst<n_residue_atoms; ifirst++) {
	    std::string first_atom_name(residue_atoms[ifirst]->GetAtomName());
	    if (first_atom_name == swapper_atoms[iswap].first) {
	       for (int isec=0; isec<n_residue_atoms; isec++) {
		  std::string sec_atom_name(residue_atoms[isec]->GetAtomName());
		  if (sec_atom_name == swapper_atoms[iswap].second) {
		     residue_atoms[ifirst]->SetAtomName(  sec_atom_name.c_str());
		     residue_atoms[isec  ]->SetAtomName(first_atom_name.c_str());
		     ifound = 1;
		     nfound++;
		  }
		  if (ifound)
		     break;
	       }
	       if (ifound)
		  break;
	    }
	 }
      }
      coot::rotamer dc(residue_copy);
      coot::rotamer_probability_info_t p_swapped = dc.probability_of_this_rotamer();

      // Delete residue_copy
      for (int ifirst=0; ifirst<n_residue_atoms; ifirst++) {
	 residue_copy->DeleteAtom(ifirst);
      }
      delete residue_copy;
      residue_copy = 0;

      // ==============================================
      
      if (p_init.state && p_swapped.state) {
	 if (0)
	    std::cout << "compare: intial orientation score " << p_init.probability
		      << "  vs swapped orientation score: "
		      << p_swapped.probability << "\n";
	 if (p_swapped.probability > p_init.probability) {
// 	    std::cout << "Fixing residue " << residue->GetSeqNum()
// 		      << " " << residue->GetChainID() << " " << residue->GetResName()
// 		      << "\n";
	    // make the changes to residue then
	    mmdb::PPAtom residue_atoms;
	    int n_residue_atoms;
	    residue->GetAtomTable(residue_atoms, n_residue_atoms);
	    int nfound = 0; 
	    for (unsigned int iswap=0; iswap<swapper_atoms.size(); iswap++) {
	       short int ifound = 0;
	       for (int ifirst=0; ifirst<n_residue_atoms; ifirst++) {
		  std::string first_atom_name(residue_atoms[ifirst]->GetAtomName());
		  if (first_atom_name == swapper_atoms[iswap].first) {
		     for (int isec=0; isec<n_residue_atoms; isec++) {
			std::string sec_atom_name(residue_atoms[isec]->GetAtomName());
			if (sec_atom_name == swapper_atoms[iswap].second) {
			   if (apply_swap_if_found) { 
			      residue_atoms[ifirst]->SetAtomName(  sec_atom_name.c_str());
			      residue_atoms[isec  ]->SetAtomName(first_atom_name.c_str());
			   }
			   ifound = 1;
			   nfound++;
			   imoved = 1;
			}
			if (ifound)
			   break;
		     }
		     if (ifound)
			break;
		  }
	       }
	    }
	    if (0)
	       std::cout << "   DEBUG:: Fixed " << nfound << " atom pairs of residue "
			 << residue->GetResName() << " " << residue->GetSeqNum()
			 << "in optimize_rotamer_by_atom_names()" << " of residue pointer "
			 << residue <<  std::endl;
	 }
      } else {
	 std::cout << "badness in atom selection " << std::endl;
      }
   }

   // Does this bit of code do anything!?  (This function is not
   // called from coot-nomenclature for LEU, VAL or THR).

   // These have Chiral Centres, we can't swap the atom names and
   // check the rotamer.  We must make sure that the chiral centre has
   // the right sign.
   //
   // Question:  What is the sign on these chiral centres?
   // What does procheck say?
   //
   // Note for LEU and VAL that these aren't a formal chiral centre
   // (because the swapper atoms are the same type) - they are just
   // "nomenclature volume" atoms (marked as "both" in the library).
   // 
   // CB in THR *is* a chiral centre and is marked as negative.
   // 
   if ((residue_name == "LEU") ||
       (residue_name == "VAL") ||
       (residue_name == "THR")) {

      std::vector<std::pair<std::string, std::string> > swapper_atoms;
      std::string chiral_centre;
      std::string other_atom; 
      if (residue_name == "LEU") {
	 chiral_centre = " CG ";
	 other_atom    = " CB ";
	 swapper_atoms.push_back(std::pair<std::string, std::string>(" CD1", " CD2"));
      }
      if (residue_name == "VAL") {
	 chiral_centre = " CB ";
	 other_atom    = " CA ";
	 swapper_atoms.push_back(std::pair<std::string, std::string>(" CG1", " CG2"));
      }
      if (residue_name == "THR") {
	 swapper_atoms.push_back(std::pair<std::string, std::string>(" OG1", " CG2"));
      }
   }
   return imoved;
}

coot::rotamer_probability_info_t
coot::rotamer::probability_of_this_rotamer() {

   short int state;  // flag for assigned,                    1 
                     // unassigned due to missing atoms,      0 
                     // unassigned due to rotamer not found. -1
                     // unassigned due to GLY/ALA            -2
   state = 0;

   // First (and mainly), what are the actual chi angles for (the
   // atoms of) this residue?

   std::string residue_name(residue->GetResName());

//    std::cout << "ROTAMER::           checking residue " << residue << " "
// 	     << residue->GetSeqNum() << " " << residue->GetChainID() << " "
// 	     << residue_name << std::endl;

   mmdb::PAtom *residue_atoms;
   int n_residue_atoms;
   residue->GetAtomTable(residue_atoms, n_residue_atoms);
   std::vector<std::vector<std::string> > rotamer_ats = rotamer_atoms(residue_name);

   if (rotamer_ats.size() == 0)
      return coot::rotamer_probability_info_t(-2, 0.0, "none"); // no chi-squared value for 
                                                                // residues with no side chain
   
   std::vector<std::vector<int> > atom_indices =
      rotamer_atom_names_to_indices(rotamer_ats, residue_atoms, n_residue_atoms);
   
   if (rotamer_ats.size() != atom_indices.size()) {
      // missing atom(s)
      return coot::rotamer_probability_info_t(0, 0.0, "none");
   }
   
   // Let's get the chi angles for this residue
   //
   std::vector<double> chi_angles;
   for (unsigned int ich=0; ich<rotamer_ats.size(); ich++) { 
      chi_angles.push_back(chi_torsion(atom_indices[ich], residue_atoms));
   }
   
   // chi_angles contains the chi angles for this residue
   // So, which rotamer are we in?
   
   std::vector<coot::simple_rotamer> rots = get_all_rotamers(residue_name);

   if (0) { // debug
      for (unsigned int ichi=0; ichi<chi_angles.size(); ichi++)
	 std::cout << "   " << chi_angles[ichi];
      std::cout << std::endl;
      std::cout << "calling probability_of_this_rotamer(chi_angles, rots)" << std::endl;
   } 
   coot::rotamer_probability_info_t d = probability_of_this_rotamer(chi_angles, rots);

   return d;
}

coot::rotamer_probability_info_t
coot::rotamer::probability_of_this_rotamer(const std::vector<double> &chi_angle_values,
					   const std::vector<coot::simple_rotamer> &rots) const {

   // std::cout << "probability_of_this_rotamer() " << rots.size() << std::endl;

   short int state;  // flag for assigned,                    1 
                     // unassigned due to missing atoms,      0 
                     // unassigned due to rotamer not found. -1
                     // unassigned due to GLY/ALA            -2
   double d = 0; // initially unassigned.
   state = -1; // rotamer not found initially
   std::string name; // rotamer name

//    std::cout << "DEBUG:: This residue's real chis: ";
//    for (int i=0; i<chi_angle_values.size(); i++)
//       std::cout << chi_angle_values[i] << " ";
//    std::cout << std::endl;
   
   if (rots.size() == 0) {
      state = -2;
   } else { 
      for (unsigned int irot=0; irot<rots.size(); irot++) {
	 int isimilar = 0;
	 for (unsigned int ich=0; ich<chi_angle_values.size(); ich++) {
	    double target = rots[irot].get_chi(ich+1);
	    double model = chi_angle_values[ich];
	    if (0) 
	       std::cout  << "DEBUG::    comparing chi " << ich << " target: " << target
			  <<  " model: " << model << std::endl;
	    if (similar_rotamer_chi(target, model))
	       isimilar++;
	    else
	       break;
	 }
	 if (0)
	    std::cout << "DEBUG::    rotamer " << irot << " isimilar = "
		      << isimilar << " out of " << chi_angle_values.size() << std::endl;
	 if (isimilar == int(chi_angle_values.size())) {
	    d = rots[irot].P_r1234();
	    name = rots[irot].rotamer_name();
	    if (d > 0) 
	       state = 1;
	 }
	 if (state == 1)
	    break;
      }
   }
   if (0)
      std::cout << "   returning from probability_of_this_rotamer() with "
		<< state << " " << d << " :" << name << ":" << std::endl;
   return coot::rotamer_probability_info_t(state, d, name);
}

std::vector<std::vector<std::string> >
coot::rotamer::rotamer_atoms(const std::string &residue_name) const {
   
   std::vector<std::vector<std::string> > r;
   
   std::vector<std::string> chi_1;
   std::vector<std::string> chi_2;
   std::vector<std::string> chi_3;
   std::vector<std::string> chi_4;
   std::vector<std::string> chi_5;

   
   chi_1.push_back(" N  "); 
   chi_1.push_back(" CA ");
   chi_1.push_back(" CB ");

   chi_2.push_back(" CA ");
   chi_2.push_back(" CB ");

   chi_3.push_back(" CB ");

   if (residue_name == "VAL") {
      chi_1.push_back(" CG1");
      // done
   }
   if (residue_name == "ASP") {
      chi_1.push_back(" CG ");
      chi_2.push_back(" CG ");
      chi_2.push_back(" OD1");
      // done
   }
   if (residue_name == "ASN") {
      chi_1.push_back(" CG ");
      chi_2.push_back(" CG ");
      chi_2.push_back(" OD1");
      // done
   }
   if (residue_name == "GLN") {
      chi_1.push_back(" CG ");
      chi_2.push_back(" CG ");
      chi_2.push_back(" CD ");
      chi_3.push_back(" CG ");
      chi_3.push_back(" CD ");
      chi_3.push_back(" OE1");
   }
   if (residue_name == "GLU") {
      chi_1.push_back(" CG ");
      chi_2.push_back(" CG ");
      chi_2.push_back(" CD ");
      chi_3.push_back(" CG ");
      chi_3.push_back(" CD ");
      chi_3.push_back(" OE1");
   }
   if (residue_name == "PHE") {
      chi_1.push_back(" CG ");
      chi_2.push_back(" CG ");
      chi_2.push_back(" CD1");
   }
   if (residue_name == "TYR") {
      chi_1.push_back(" CG ");
      chi_2.push_back(" CG ");
      chi_2.push_back(" CD1");
   }
   if (residue_name == "TRP") {
      chi_1.push_back(" CG ");
      chi_2.push_back(" CG ");
      chi_2.push_back(" CD1");
   }
   if (residue_name == "THR") {
      chi_1.push_back(" OG1");
   }
   if (residue_name == "LEU") {
      chi_1.push_back(" CG ");
      chi_2.push_back(" CG ");
      chi_2.push_back(" CD1");
   }
   if (residue_name == "ILE") {
      chi_1.push_back(" CG1");
      chi_2.push_back(" CG1");
      chi_2.push_back(" CD1");
   }

   if (residue_name == "ARG") {
      chi_1.push_back(" CG ");
      chi_2.push_back(" CG ");
      chi_2.push_back(" CD ");
      chi_3.push_back(" CG ");
      chi_3.push_back(" CD ");
      chi_3.push_back(" NE ");
      chi_4.push_back(" CG ");
      chi_4.push_back(" CD ");
      chi_4.push_back(" NE ");
      chi_4.push_back(" CZ ");
   }
   if (residue_name == "CYS") {
      chi_1.push_back(" SG ");
   }
   if (residue_name == "HIS") {
      chi_1.push_back(" CG ");
      chi_2.push_back(" CG ");
      chi_2.push_back(" ND1");
   }
   if (residue_name == "LYS") {
      chi_1.push_back(" CG ");
      chi_2.push_back(" CG ");
      chi_2.push_back(" CD ");
      chi_3.push_back(" CG ");
      chi_3.push_back(" CD ");
      chi_3.push_back(" CE ");
      chi_4.push_back(" CG ");
      chi_4.push_back(" CD ");
      chi_4.push_back(" CE ");
      chi_4.push_back(" NZ ");
   }
   if (residue_name == "MET") {
      chi_1.push_back(" CG ");
      chi_2.push_back(" CG ");
      chi_2.push_back(" SD ");
      chi_3.push_back(" CG ");
      chi_3.push_back(" SD ");
      chi_3.push_back(" CE ");
   }
   if (residue_name == "MSE") {
      chi_1.push_back(" CG ");
      chi_2.push_back(" CG ");
      chi_2.push_back("SE  ");
      chi_3.push_back(" CG ");
      chi_3.push_back("SE  ");
      chi_3.push_back(" CE ");
   }
   if (residue_name == "PRO") {
      chi_1.push_back(" CG ");
      chi_2.push_back(" CG ");
      chi_2.push_back(" CD ");
   }
   if (residue_name == "SER") {
      chi_1.push_back(" OG ");
      // done
   }

   if (chi_1.size() == 4)
      r.push_back(chi_1);
   if (chi_2.size() == 4)
      r.push_back(chi_2);
   if (chi_3.size() == 4)
      r.push_back(chi_3);
   if (chi_4.size() == 4)
      r.push_back(chi_4);

   return r;

} 

std::vector<std::vector<int> >
coot::rotamer::rotamer_atom_names_to_indices(const std::vector<std::vector<std::string> > &residue_rotamer_atoms,
					      mmdb::PAtom *residue_atoms,
					      int n_residue_atoms) const {

   bool verbose = false;
   std::vector<std::string> atom_indices(n_residue_atoms);
   std::vector<std::vector<int> > r;

   if (n_residue_atoms > 0) { 
      for (int iat=0; iat<n_residue_atoms; iat++) {
	 atom_indices[iat] = std::string(residue_atoms[iat]->name);
      }

      int ithis_atom;
      // for each chi atom name set:
      for (unsigned int ich=0; ich<residue_rotamer_atoms.size(); ich++) {
	 std::vector<int> single;
	 // for each atom name (in the chi set):
	 for (unsigned int iat=0; iat<residue_rotamer_atoms[ich].size(); iat++) {
	    ithis_atom = -1;
	    // Can we find in atom_indices the atom with the same name as
	    // the (dictionary) chi atom name?
	    for (int irat=0; irat<n_residue_atoms; irat++) {
	       if (atom_indices[irat] == residue_rotamer_atoms[ich][iat]) {
		  ithis_atom = irat;
		  break;
	       }
	    }
	    if (ithis_atom != -1) {
	       single.push_back(ithis_atom);
	    }
	 }
	 if (single.size() == 4) 
	    r.push_back(single);
	 else
	    if (verbose)
	       std::cout << "PROBLEM in coordinates file? failed to find all atoms in "
			 << "ich number " << ich << " "
			 << residue_atoms[0]->residue->name << " "
			 << residue_atoms[0]->residue->GetSeqNum() << " "
			 << residue_atoms[0]->residue->GetChainID() << std::endl;
      }
   }
   return r;
}





// return in degrees
double
coot::rotamer::chi_torsion(const std::vector<int> &chi_angle_atoms_indices,
			    mmdb::PAtom *residue_atoms) {

   double tors = 0.0;

   std::vector<clipper::Coord_orth> a;
   mmdb::Atom *at;

   for (unsigned int ich_at=0; ich_at<chi_angle_atoms_indices.size(); ich_at++) {
      at = residue_atoms[chi_angle_atoms_indices[ich_at]];
      a.push_back(clipper::Coord_orth(at->x, at->y, at->z));
   }

   double ctorsion = clipper::Coord_orth::torsion(a[0], a[1], a[2], a[3]);
   tors = clipper::Util::rad2d(ctorsion);

   return tors;
} 


std::vector<coot::simple_rotamer>
coot::rotamer::get_all_rotamers(const std::string &res_type) const {

   
   for(unsigned int i=0; i< typed_rotamers.size(); i++) {
      if (typed_rotamers[i].Type() == res_type) {
	 return typed_rotamers[i].simple_rotamers();
      }
   }
   std::vector<coot::simple_rotamer> r;
   return r;

}

// Return NULL if no residues available for this residue type
// 
// caller needs to delete returned residue and its chain
// 
mmdb::Residue *
coot::rotamer::GetResidue(const coot::dictionary_residue_restraints_t &rest,
			   int i_rot) const {

   unsigned int ui_rot = i_rot;
   bool debug = 0;
   // std::cout << "rotamer::GetResidue alt_conf is :" << alt_conf << ":" << std::endl;
   mmdb::Residue *rres = deep_copy_residue(Residue()); 
   std::string rt = Residue_Type();

   std::vector<coot::simple_rotamer> rots = rotamers(rt, probability_limit);

   if (debug) { // debug
      mmdb::PPAtom residue_atoms;
      int n_residue_atoms;
      Residue()->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
	 std::cout << "debug in rotamer::GetResidue(): " << iat << " "
		   << residue_atoms[iat] << std::endl;
      }
   } 

   if ((rots.size() == 0) || (ui_rot >= rots.size())) { 
      return rres; // or should this be null? 
   }
   coot::simple_rotamer this_rot = rots[ui_rot];
   set_dihedrals(rres, rest, this_rot);
   
   return rres;
}

mmdb::Residue *
coot::rotamer::GetResidue(const coot::dictionary_residue_restraints_t &rest,
			  const std::string &rotamer_name) const {

   mmdb::Residue *r = NULL; // returned value
   mmdb::Residue *rres = deep_copy_residue(Residue());
   if (rres) { 
      std::string rt = Residue_Type();
      std::vector<coot::simple_rotamer> rots = rotamers(rt, probability_limit);
      for (unsigned int i=0; i<rots.size(); i++) { 
	 if (rots[i].rotamer_name() == rotamer_name) {
	    const simple_rotamer &this_rot = rots[i];
	    set_dihedrals(rres, rest, this_rot);
	    r = rres;
	    break;
	 } 
      }
   }
   return r;
}

// move the atoms of rres
void
coot::rotamer::set_dihedrals(mmdb::Residue *rres,
			     const coot::dictionary_residue_restraints_t &rest,
			     const simple_rotamer &this_rot) const {

   bool debug = false;
   std::vector<coot::atom_name_quad> atom_name_quads = atom_name_quad_list(Residue_Type());

   if (debug) { 
      for (unsigned int ichi=0; ichi<atom_name_quads.size(); ichi++) {
	 std::cout << "Quad " << ichi << " "
		   << atom_name_quads[ichi].atom_name(0) << " "
		   << atom_name_quads[ichi].atom_name(1) << " "
		   << atom_name_quads[ichi].atom_name(2) << " "
		   << atom_name_quads[ichi].atom_name(3) << " "
		   << std::endl;
      }
   }

   for (unsigned int ichi=0; ichi<atom_name_quads.size(); ichi++) {
      double tors = this_rot[ichi];
      try {
	 coot::atom_tree_t t(rest, rres, alt_conf);
	 double new_angle = t.set_dihedral(atom_name_quads[ichi].atom_name(0),
					   atom_name_quads[ichi].atom_name(1),
					   atom_name_quads[ichi].atom_name(2),
					   atom_name_quads[ichi].atom_name(3),
					   tors);

	 if (debug) 
	    std::cout << "DEBUG:: in rotamer::GetResidue() setting torsion " 
		      << tors << " of " << atom_name_quads[ichi] << " results in angle "
		      << new_angle << std::endl;

      }
      catch (const std::runtime_error &rte) {
	 std::cout << "oops! in rotamer::GetResidue() " << rte.what() << std::endl;
      } 
   }
} 


// Return NULL if no residues available for this residue type
// 
mmdb::Residue *
coot::rotamer::GetResidue_old(int i_rot) const {

   mmdb::Residue *rres = deep_copy_residue(Residue());
   unsigned int ui_rot = i_rot;

//    std::cout << "debug:: deep_copy_residue of "
// 	     << Residue()->GetSeqNum() << Residue()->GetInsCode()
// 	     << " -> "
// 	     << rres->GetSeqNum() << rres->GetInsCode()
// 	     << "\n";

   std::string rt = Residue_Type();
   if (rt == "MSE")
      rt = "MET";
   std::vector<coot::simple_rotamer> rots = rotamers(rt, probability_limit);

   if ((rots.size() == 0) || (ui_rot >= rots.size())) { 
      return rres; // or should this be null? 
   }

   coot::simple_rotamer this_rot = rots[ui_rot];
   
   mmdb::PPAtom residue_atoms;
   int nResidueAtoms;
   rres->GetAtomTable(residue_atoms, nResidueAtoms); 
   if (nResidueAtoms == 0) {
      std::cout << " something broken in atom residue selection in ";
      std::cout << "GetResidue, got 0 atoms" << std::endl;
   } else {

      // Let's make the coordinates:
      //
      std::vector<mmdb::Atom *> ordered_atoms = ordered_residue_atoms(rres);
//       std::cout << "DEBUG:: nResidueAtoms " << nResidueAtoms
// 		<< " and ordered_atoms has size " << ordered_atoms.size() << std::endl;


      // debugging
      for(int i=0; i<nResidueAtoms; i++) {
// 	 smn_Cartesian c(residue_atoms[i]->x,
// 			 residue_atoms[i]->y,
// 			 residue_atoms[i]->z);
	 ::Cartesian d(ordered_atoms[i]->x,
		       ordered_atoms[i]->y,
		       ordered_atoms[i]->z);
// 	 std::cout << residue_atoms[i]->name << " " << c
// 		   << ordered_atoms[i]->name << " " << d << std::endl;
//	 std::cout << ordered_atoms[i]->name << " " << d << std::endl;
      }

      int nres_atoms = nResidueAtoms;
      mmdb::PAtom *ordered_residue_atoms_ppcatom = new mmdb::PAtom[nres_atoms];
      for(int i=0; i<nResidueAtoms; i++)
	 ordered_residue_atoms_ppcatom[i] = ordered_atoms[i];
      
      std::vector< ::Cartesian > coords;
      for(int i=0; i<nResidueAtoms; i++) {
	 ::Cartesian c(ordered_atoms[i]->x,
		       ordered_atoms[i]->y,
		       ordered_atoms[i]->z);
	 coords.push_back(c);
      }

      // Rotatable bonds -> coord indices
      // 
      std::vector<coot::atom_name_pair> atom_name_pairs =
	 atom_name_pair_list(Residue_Type());
      std::vector<coot::atom_index_pair> atom_index_pairs =
	 get_atom_index_pairs(atom_name_pairs, ordered_residue_atoms_ppcatom, nResidueAtoms);

      std::vector<coot::atom_name_quad> atom_name_quads =
	 atom_name_quad_list(Residue_Type());
      std::vector<coot::atom_index_quad> atom_index_quads =
	 get_atom_index_quads(atom_name_quads, ordered_residue_atoms_ppcatom, nResidueAtoms);
	 
	 
      
//       //
//       for (int i=0; i<atom_name_pairs.size(); i++) { 
// 	 std::cout << "DEBUG:: atom name pair " << i << " " 
// 		   << atom_name_pairs[i].atom1 << " " 
// 		   << atom_name_pairs[i].atom2 << std::endl;
//       }
//       // atom index pairs
//       for (int i=0; i<atom_index_pairs.size(); i++) { 
// 	 std::cout << "DEBUG:: atom index pair " << i << " " 
//  		   << atom_index_pairs[i].index1 << " "
//  		   << atom_index_pairs[i].index2 << std::endl;
//       } 

      // Contact indices:
      //
      atom_selection_container_t res_asc;
      res_asc.mol = stored_mol;
      res_asc.n_selected_atoms = nResidueAtoms;
      res_asc.atom_selection = ordered_residue_atoms_ppcatom;


      if (! res_asc.mol) {
	 std::cout << "WARNING NULL stored_mol in GetResidue() " << std::endl;

      } else { 

	 // mmdb_extras function
	 contact_info contact = getcontacts(res_asc);

	 std::vector<std::vector<int> > contact_indices(nResidueAtoms);
	 // std::cout << "Found " << contact_indices.size() << " atoms" << std::endl;

	 int in;
	 for (int i=0; i<contact.n_contacts(); i++) {
	    in = contact.contacts[i].id1;
	    contact_indices[contact.contacts[i].id2].push_back(in);
	 }
	 // std::cout << "filled " << contact_indices.size() << " contact_indices" << std::endl;

	 // debugging
	 if (0) { 
	    std::cout << " -----------  coords ---------------- " << std::endl;
	    for(unsigned int i=0; i<coords.size(); i++)
	       std::cout << i << " " << ordered_residue_atoms_ppcatom[i]->GetAtomName()
			 << " " << coords[i] << std::endl;
	 
	    // display contact indices here:
	    //
	    std::cout << " ---------------- contact indices ----------------------\n" ;
	    for (unsigned int ic=0; ic<contact_indices.size(); ic++) {
	       for (unsigned int jc=0; jc<contact_indices[ic].size(); jc++) {
		  std::cout << " contact " << ic << " "
			    << contact_indices[ic][jc] << std::endl;
	       }
	    }
	 }

	 //
	 float tors;
	 Tree tree;
	 tree.SetCoords(coords, 0, contact_indices);
	 
	 // test that the atom indices are sensible?
	 //
	 for(unsigned int ibond=0; ibond<atom_index_quads.size(); ibond++) { 
	    
	    tors = clipper::Util::d2rad(this_rot[ibond]);

	    if (tree.GetCoord(atom_index_pairs[ibond].index2)->GetNumberOfChildren() == 0) { 
	       std::cout << "WARNING: Can't torsion this bond - atom number " 
			 << atom_index_pairs[ibond].index2 << " " 
			 << ordered_residue_atoms_ppcatom[atom_index_pairs[ibond].index2] 
			 << " has 0 children" << std::endl;
	    } else { 
	       double ctors = tree.GetCoord(atom_index_pairs[ibond].index2)->GetChild(0)->GetParentDihedralAngle();
	       ::Cartesian child_0_coord = tree.GetCoord(atom_index_pairs[ibond].index2)->GetChild(0)->GetCoord();

	       // 	    std::cout << "DEBUG:: current " << ibond << " tors is: " << ctors 
	       // 		      << " (rads) " << clipper::Util::rad2d(ctors) << "  (degrees)" 
	       // 		      << " for child at "  << child_0_coord;

	       // 	    double dih1 = DihedralAngle(coords[0],
	       // 					coords[3],
	       // 					coords[2],
	       // 					coords[1]);
	       // 	    double dih2 = DihedralAngle(coords[6],
	       // 					coords[3],
	       // 					coords[2],
	       // 					coords[6]);

	       // 	    std::cout << "   Hand calculation of dihedrals: (pre-rotation) CG1: "
	       // 		      << clipper::Util::rad2d(dih1) << " CG2: " 
	       // 		      << clipper::Util::rad2d(dih2) << std::endl;
	    
	       // 	    std::cout << ":::: rotamer " << i_rot << " rotating "
	       // 		      << this_rot[ibond] << " for Chi_" << ibond+1
	       // 		      << " atom numbers: "
	       // 		      << atom_index_pairs[ibond].index1 << " "
	       // 		      << atom_index_pairs[ibond].index2 << std::endl;


	       // old style	    
	       // 	    tree.RotateAboutBond(atom_index_pairs[ibond].index1,
	       // 				 atom_index_pairs[ibond].index2, tors);

	       // 2-atom style
	       // 	    tree.SetDihedralAngle(atom_index_pairs[ibond].index1,
	       //				  atom_index_pairs[ibond].index2, tors);

	       // 20090213
	       tree.SetDihedralAngle(atom_index_quads[ibond].index1,
				     atom_index_quads[ibond].index2,
				     atom_index_quads[ibond].index3,
				     atom_index_quads[ibond].index4,
				     tors);

	       ctors = tree.GetCoord(atom_index_pairs[ibond].index2)->GetChild(0)->GetParentDihedralAngle();
	       child_0_coord = tree.GetCoord(atom_index_pairs[ibond].index2)->GetChild(0)->GetCoord();
	       // 	    std::cout << "DEBUG:: post rot " << ibond << " tors is: " << ctors 
	       // 		      << " " << clipper::Util::rad2d(ctors) << std::endl
	       // 		      << "     for child at " << child_0_coord;

	    }
	 }
	 
	 std::vector< ::Cartesian > coords_rotated =
	    tree.GetAllCartesians();

	 //       double dih1 = DihedralAngle(coords_rotated[0],
	 // 				  coords_rotated[3],
	 // 				  coords_rotated[2],
	 // 				  coords_rotated[1]);
	 //       double dih2 = DihedralAngle(coords_rotated[0],
	 // 				  coords_rotated[3],
	 // 				  coords_rotated[2],
	 // 				  coords_rotated[6]);
      
	 //       std::cout << "   Hand calculation of dihedrals: (post-rotation) CG1: "
	 // 		<< clipper::Util::rad2d(dih1) << " CG2: " 
	 // 		<< clipper::Util::rad2d(dih2) << std::endl;
      

	 
	 // debugging
	 //        for(int i=0; i<coords_rotated.size(); i++)
	 // 	  std::cout << i << " " << ordered_residue_atoms_ppcatom[i]->GetAtomName()
	 // 		    << " " << coords_rotated[i];
	 
	 if (int(coords_rotated.size()) != nResidueAtoms) {
	    std::cout << "disaster in atom selection, trees, dunbrack\n";
	 } else {
	    for (int iat=0; iat<nResidueAtoms; iat++) {
	       ordered_residue_atoms_ppcatom[iat]->x = coords_rotated[iat].get_x();
	       ordered_residue_atoms_ppcatom[iat]->y = coords_rotated[iat].get_y();
	       ordered_residue_atoms_ppcatom[iat]->z = coords_rotated[iat].get_z();
	    }
	 }
	 delete [] ordered_residue_atoms_ppcatom;
      }
   }
   return rres;
}

// The dunbrack rotamers (and the interaction with mgtree) depend on
// the atom ordering, e.g. for VAL CG1 should come before CG2.
// 
std::vector<mmdb::Atom *>
coot::rotamer::ordered_residue_atoms(mmdb::Residue *residue_p) const {

   std::vector<mmdb::Atom *> atom_vec;
   std::vector<mmdb::Atom *> store_vec;


   std::string residue_name(residue_p->GetResName());
   std::vector <std::string> atom_order;

   if (residue_name == "VAL") {
      atom_order.push_back(" CG1");
      atom_order.push_back(" CG2");
   }
   if (residue_name == "ASP") {
      atom_order.push_back(" OD1");
      atom_order.push_back(" OD2");
   }
   if (residue_name == "ASN") {
      atom_order.push_back(" OD1");
      atom_order.push_back(" ND2");
   }
   if (residue_name == "GLN") {
      atom_order.push_back(" OE1");
      atom_order.push_back(" NE2");
   }
   if (residue_name == "GLU") {
      atom_order.push_back(" OE1");
      atom_order.push_back(" OE2");
   }
   if (residue_name == "PHE") {
      atom_order.push_back(" CD1");
      atom_order.push_back(" CD2");
   }
   if (residue_name == "TYR") {
      atom_order.push_back(" CD1");
      atom_order.push_back(" CD2");
   }
   if (residue_name == "TRP") {
      atom_order.push_back(" CD1");
      atom_order.push_back(" CD2");
   }
   if (residue_name == "THR") {
      atom_order.push_back(" OG1");
      atom_order.push_back(" CG2");
   }
   if (residue_name == "LEU") {
      atom_order.push_back(" CD1");
      atom_order.push_back(" CD2");
   }
   if (residue_name == "ILE") {
      atom_order.push_back(" CG1");
      atom_order.push_back(" CG2");
   }
   if (residue_name == "PRO") {
      atom_order.push_back(" N  ");
      atom_order.push_back(" CA ");
      atom_order.push_back(" CB ");
      atom_order.push_back(" CG ");
      atom_order.push_back(" CD ");
      atom_order.push_back(" C  ");
      atom_order.push_back(" O  ");
   }

   if (atom_order.size() == 0) {
      // not an atom with possible ordering problems (e.g SER)
      mmdb::PPAtom residue_atoms;
      int nResidueAtoms;
      residue_p->GetAtomTable(residue_atoms, nResidueAtoms);
      for(int i=0; i<nResidueAtoms; i++)
	 atom_vec.push_back(residue_atoms[i]);

   } else { 

   
      mmdb::PPAtom residue_atoms;
      int nResidueAtoms;
      residue_p->GetAtomTable(residue_atoms, nResidueAtoms);
      if (nResidueAtoms == 0) {
	 std::cout << " something broken in atom residue selection in ";
	 std::cout << "ordered_residue_atoms, got 0 atoms" << std::endl;
      } else {
	 int sensitive_count = 0;
	 short int needs_reorder = 0;
	 std::string atom_name;
	 if (atom_order.size() == 2) { 
	    for(int i=0; i<nResidueAtoms; i++) {

	       // We are concerned about there being atom_order[0] atoms
	       // below (i.e. have a higher index than) atom_order[1]
	       // atoms.  This is the problem we need to fix.
	    
	       // So when we see an atom_order[1] atom, we set the
	       // sensitive_count and then if sensitive_count is set and we
	       // see an atom_order[0] atom, then we know we have to
	       // reorder.
	       //
	       // sensitive_count also counts the number of atoms that have
	       // name atom_order[0] (in the case of altlocs) that should be
	       // before atom_order[1] atoms.
	       // 

	       atom_name = std::string(residue_atoms[i]->GetAtomName());

	       if (atom_name == atom_order[1]) {
		  sensitive_count++;
	       } else {
		  if (sensitive_count && atom_name == atom_order[0]) {
		     needs_reorder = 1;
		  }
	       }
	    }

// 	    std::cout << "DEBUG:: needs_reorder: " << needs_reorder
// 		      << std::endl;
	    if (needs_reorder) {
	       // 	       std::cout << "DEBUG:: residue needs reordering\n";
	       int passed_sensitive_atoms = 0;
	       std::string atom_name;
	       for(int i=0; i<nResidueAtoms; i++) {
		  atom_name = std::string(residue_atoms[i]->GetAtomName());
		  if (passed_sensitive_atoms < sensitive_count) {
		     // we are the top of the list still
		     if (atom_name == atom_order[0]) {
			passed_sensitive_atoms++;
			atom_vec.push_back(residue_atoms[i]);
		     } else {
			if (atom_name == atom_order[1]) {
			   // need to store this one...
			   store_vec.push_back(residue_atoms[i]);
			} else {
			   atom_vec.push_back(residue_atoms[i]);
			}
		     }
		  } else {
		     atom_vec.push_back(residue_atoms[i]);
		  }
		  if (passed_sensitive_atoms == sensitive_count) {
		     for (unsigned int j=0; j<store_vec.size(); j++) {
			atom_vec.push_back(store_vec[j]);
		     }
		     passed_sensitive_atoms++; // so this test fails next time
		  }
	       }
	    } else {
	       for(int i=0; i<nResidueAtoms; i++)
		  atom_vec.push_back(residue_atoms[i]);
	    }
	 } else {

	    // more than 2 atoms to be ordered.  Run through the
	    // residue atoms twice.  Once to pick up the atoms and the
	    // second to order them.
	    std::vector<mmdb::Atom *> other_atoms;
	    for (unsigned int ian=0; ian<atom_order.size(); ian++) {
	       bool marked_atom = false;
	       for(int i=0; i<nResidueAtoms; i++) {
		  atom_name = std::string(residue_atoms[i]->GetAtomName());
		  if (atom_name == atom_order[ian]) {
		     atom_vec.push_back(residue_atoms[i]);
		     marked_atom = true;
		  }
	       }
	    }
	    
	    for(int i=0; i<nResidueAtoms; i++) {
	       bool found = false; 
	       for (unsigned int in_atom=0; in_atom<atom_vec.size(); in_atom++) {
		  if (atom_vec[in_atom] == residue_atoms[i]) {
		     found = true;
		     break;
		  }
	       }
	       if (found == false)
		  other_atoms.push_back(residue_atoms[i]);
	    }
	    // now merge other_atoms onto the end of atom_vec
	    for (unsigned int i=0; i<other_atoms.size(); i++)
	       atom_vec.push_back(other_atoms[i]);
	 } 
      }
   }
//    std::cout << "ordered residue atoms" << std::endl;
//    for (int i=0; i< atom_vec.size(); i++)
//       std::cout << i << " " << atom_vec[i] << std::endl;
   
   return atom_vec;
} 

std::vector<coot::simple_rotamer>
coot::rotamer::rotamers(const std::string &res_type, float prob_cut) const {

   std::vector<coot::simple_rotamer> rots;

   short int found_res = 0;
   for(unsigned int i=0; i< typed_rotamers.size(); i++) {
      if (typed_rotamers[i].Type() == res_type) {
	 found_res = 1;
	 rots = typed_rotamers[i].get_sorted_rotamers(prob_cut);
	 break;
      }
   }
   return rots;
}


// if we cared about efficiency, we'd tidy this up a bit.
float
coot::rotamer::Chi1(int irot) const {

   float v = -999;

   for (unsigned int i=0; i<typed_rotamers.size(); i++) {
      std::string rt = Residue_Type();
      if (typed_rotamers[i].Type() == rt) {
	 if (irot<int(rotamers(rt, Probability_limit()).size())) { 
	    v = rotamers(rt, Probability_limit())[irot].Chi1();
	    return v;
	 } else {
	    std::cout << "ERROR: asked for index " << irot << " but max rotamers was "
		      << rotamers(rt, Probability_limit()).size() << std::endl;
	 }
      }
   }
   return v;
}


std::string
coot::rotamer::rotamer_name(int irot) {

   std::string n = "";

   for (unsigned int i=0; i<typed_rotamers.size(); i++) {
      std::string rt = Residue_Type();
      if (typed_rotamers[i].Type() == rt) {
	 if (irot<int(rotamers(rt, Probability_limit()).size())) { 
	    n = rotamers(rt, Probability_limit())[irot].rotamer_name();
	    break;
	 } else {
	    std::cout << "ERROR: asked for index " << irot << " but max rotamers was "
		      << rotamers(rt, Probability_limit()).size() << std::endl;
	 }
      }
   }
   return n;
} 



