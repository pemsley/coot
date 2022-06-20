/* src/chi-angles.cc
 * 
 * Copyright 2001, 2002, 2003, 2004, 2006 The University of York
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

#include "compat/coot-sysdep.h"
#include "coot-utils/coot-coord-extras.hh"

#include "chi-angles.hh"
#include "ccp4mg-utils/mgtree.h"

 
coot::simple_rotamer::simple_rotamer(int rot1_in,  
				     int rot2_in,  
				     int rot3_in,  
				     int rot4_in,  
				     int n_r1_in,
				     int nr1234_in,
				     float p_r1234_in,
				     float sig_p_r1234_in,
				     float pr234_given_r1_in,
				     float sig_pr234_given_r1_in,
				     float chi1_in,
				     float sig_chi1_in,
				     float chi2_in,
				     float sig_chi2_in,
				     float chi3_in,
				     float sig_chi3_in,
				     float chi4_in,
				     float sig_chi4_in) {
   
   rotamer_type = coot::simple_rotamer::DUNBRACK_ROTAMER;
   name = ""; // no names for dunbrack rotamers
   rot1 =       	          rot1_in;  		      
   rot2 =       	          rot2_in;  		      
   rot3 =       	          rot3_in;  		      
   rot4 =       	          rot4_in;  		      
   n_r1 =       	          n_r1_in;		      
   nr1234 =     	          nr1234_in;		      
   p_r1234 =    	          p_r1234_in;		      
   sig_p_r1234 =	          sig_p_r1234_in;	      
   pr234_given_r1 =	          pr234_given_r1_in;	      
   sig_pr234_given_r1 =           sig_pr234_given_r1_in;    
   chi1 =       	          chi1_in;		      
   sig_chi1 =                     sig_chi1_in;	      
   chi2 =       	          chi2_in;		      
   sig_chi2 =                     sig_chi2_in;	      
   chi3 =       	          chi3_in;		      
   sig_chi3 =                     sig_chi3_in;	      
   chi4 =       	          chi4_in;		      
   sig_chi4 =                     sig_chi4_in;
   minus_one = -1;

}

// Constructor for richardson rotamer
coot::simple_rotamer::simple_rotamer(std::string rotamer_name,  
				     float percent_overall,
				     float percent_alpha,
				     float percent_beta,
				     float percent_other,
				     float chi_1_mode,
				     float chi_1_com,
				     float chi_2_mode,
				     float chi_2_com,
				     float chi_3_mode,
				     float chi_3_com,
				     float chi_4_mode,
				     float chi_4_com) {

   rotamer_type = coot::simple_rotamer::RICHARDSON_ROTAMER;
   name = rotamer_name;
   rot1 = 0;
   rot2 = 0; 
   rot3 = 0;
   rot4 = 0;
   p_r1234 = percent_overall;
   sig_p_r1234 = 0.0;
   chi1 = chi_1_mode;
   chi2 = chi_2_mode;
   chi3 = chi_3_mode;
   chi4 = chi_4_mode;

   nr1234  = -1; // not used

   // the awk program put in dummy -5555 for chi_mode values that have
   // not been assigned.  This is because they are "additionally
   // allowed" - not real rotamers.
   // 
   if (chi_1_mode < -555)
      chi1 = chi_1_com;
   if (chi_2_mode < -555)
      chi2 = chi_2_com;
   if (chi_3_mode < -555)
      chi3 = chi_3_com;
   if (chi_4_mode < -555)
      chi4 = chi_4_com;
   
   sig_chi1 = 40.0;
   sig_chi2 = 40.0;
   sig_chi3 = 40.0;
   sig_chi4 = 40.0;
}

coot::simple_rotamer
coot::simple_rotamer::rotate_chi2_180() const {

   float new_chi2 = chi2 + 180.0;
   if (chi2 > 0)
      new_chi2 = chi2 - 180.0;

   simple_rotamer r( rot1,  
		     rot2,  
		     rot3,  
		     rot4,  
		     n_r1,
		     nr1234,
		     p_r1234,
		     sig_p_r1234,
		     pr234_given_r1,
		     sig_pr234_given_r1,
		     chi1,
		     sig_chi1,
		     new_chi2, // here
		     sig_chi2,
		     chi3,
		     sig_chi3,
		     chi4,
		     sig_chi4);

   return r;

}

void
coot::chi_angles::add_IUPAC_extras_PHE_and_TYR_rotamers() {

   std::vector<std::string> restype;
   restype.push_back("PHE");
   restype.push_back("TYR");
   
   for (int irestype=0; irestype<2; irestype++) {
      for (unsigned int i=0; i<typed_rotamers.size(); i++) {
	 if (typed_rotamers[i].Type() == restype[irestype]) {

	    // look for rotamers with fabs(chi2) > 70:
	    std::vector<simple_rotamer> simple_rotamers = typed_rotamers[i].get_simple_rotamers();
	    for (unsigned int irot=0; irot<simple_rotamers.size(); irot++) {
	       if (fabs(simple_rotamers[irot].Chi2()) > 70.0) {
		  // create a simple rotamer and add it to typed_rotamers[i]
		  coot::simple_rotamer r = simple_rotamers[irot].rotate_chi2_180();
		  typed_rotamers[i].add_simple_rotamer(r);
	       }
	    }
	    break;
	 }
      }
   }
}


void
coot::chi_angles::setup_chi_atom_quads() {

   add_chi_quad("VAL", " N  ", " CA ", " CB ", " CG1");

   add_chi_quad("TYR", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("TYR", " CA ", " CB ", " CG ", " CD1");

   add_chi_quad("TRP", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("TRP", " CA ", " CB ", " CG ", " CD1");

   add_chi_quad("THR", " N  ", " CA ", " CB ", " OG1");

   add_chi_quad("SER", " N  ", " CA ", " CB ", " OG ");

   add_chi_quad("PRO", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("PRO", " CA ", " CB ", " CG ", " CD ");

   add_chi_quad("PHE", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("PHE", " CA ", " CB ", " CG ", " CD1");
   
   add_chi_quad("MET", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("MET", " CA ", " CB ", " CG ", " SD ");
   add_chi_quad("MET", " CB ", " CG ", " SD ", " CE ");

   // Use the monomer dictionary bonding instead   
   add_chi_quad("MSE", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("MSE", " CA ", " CB ", " CG ", "SE  ");
   add_chi_quad("MSE", " CB ", " CG ", "SE  ", " CE ");

   add_chi_quad("LYS", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("LYS", " CA ", " CB ", " CG ", " CD ");
   add_chi_quad("LYS", " CB ", " CG ", " CD ", " CE ");
   add_chi_quad("LYS", " CG ", " CD ", " CE ", " NZ ");

   add_chi_quad("LEU", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("LEU", " CA ", " CB ", " CG ", " CD1");

   add_chi_quad("ILE", " N  ", " CA ", " CB ", " CG1");
   add_chi_quad("ILE", " CA ", " CB ", " CG1", " CD1");

   add_chi_quad("HIS", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("HIS", " CA ", " CB ", " CG ", " ND1");

   add_chi_quad("GLU", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("GLU", " CA ", " CB ", " CG ", " CD ");
   add_chi_quad("GLU", " CB ", " CG ", " CD ", " OE1");

   add_chi_quad("GLN", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("GLN", " CA ", " CB ", " CG ", " CD ");
   add_chi_quad("GLN", " CB ", " CG ", " CD ", " OE1");

   add_chi_quad("CYS", " N  ", " CA ", " CB ", " SG ");
   
   add_chi_quad("ASP", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("ASP", " CA ", " CB ", " CG ", " OD1");

   add_chi_quad("ASN", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("ASN", " CA ", " CB ", " CG ", " OD1");

   add_chi_quad("ARG", " N  ", " CA ", " CB ", " CG ");
   add_chi_quad("ARG", " CA ", " CB ", " CG ", " CD ");
   add_chi_quad("ARG", " CB ", " CG ", " CD ", " NE ");
   add_chi_quad("ARG", " CG ", " CD ", " NE ", " CZ ");

}

// Return chi angle pair for the chi angles in setup_chi_atom_pairs
// (above).  Not used by generic ligands (as far as I can see).
void
coot::chi_angles::add_chi_quad(const std::string &residue_type,
			       const std::string &atom_name_1,
			       const std::string &atom_name_2,
			       const std::string &atom_name_3,
			       const std::string &atom_name_4) {

   short int found_res = 0;
   for(unsigned int i=0; i< typed_rotamers.size(); i++) {
      if (typed_rotamers[i].Type() == residue_type) {
	 found_res = 1;
	 // chi pairs:
	 typed_rotamers[i].add_torsion_bond_by_name(atom_name_2, atom_name_3);
	 typed_rotamers[i].add_torsion_bond_by_name(atom_name_1, atom_name_2,
						    atom_name_3, atom_name_4);
      } 
   }

   if (found_res == 0) {
      std::cout << "Oops, " << residue_type << " not found in typed_rotamers"
		<< std::endl;
   }
}


std::ostream&
coot::operator<<(std::ostream &s, coot::simple_rotamer rot) {

   if (rot.rotamer_type == simple_rotamer::RICHARDSON_ROTAMER) {

      s << " chi1: " << rot.chi1
        << " chi2: " << rot.chi2
        << " chi3: " << rot.chi3
        << " chi4: " << rot.chi4;

   } else {

      s <<
      rot.rot1 << " " <<     
      rot.rot2 << " " <<      
      rot.rot3 << " " <<      
      rot.rot4 << " " <<      
      rot.n_r1 << " " <<      
      rot.nr1234 << " " <<    
      rot.p_r1234 << " " <<   
      rot.sig_p_r1234 << " " <<   
      rot.pr234_given_r1 << " " <<   
      rot.sig_pr234_given_r1 << " " <<   
      rot.chi1 << " " <<     
      rot.sig_chi1 << " " <<  
      rot.chi2 << " " <<     
      rot.sig_chi2 << " " <<  
      rot.chi3 << " " <<      
      rot.sig_chi3 << " " <<  
      rot.chi4 << " " <<      
      rot.sig_chi4;
   }

   return s;
} 




void
coot::chi_angles::add_rotamer(std::string restype,
			    int rot1,  
			    int rot2,  
			    int rot3,  
			    int rot4,  
			    int n_r1,
			    int nr1234,
			    float p_r1234,
			    float sig_p_r1234,
			    float pr234_given_r1,
			    float sig_pr234_given_r1,
			    float chi1,
			    float sig_chi1,
			    float chi2,
			    float sig_chi2,
			    float chi3,
			    float sig_chi3,
			    float chi4,
			    float sig_chi4) {

   coot::simple_rotamer rot(rot1,  
			    rot2,  
			    rot3,  
			    rot4,  
			    n_r1,
			    nr1234,
			    p_r1234,
			    sig_p_r1234,
			    pr234_given_r1,
			    sig_pr234_given_r1,
			    chi1,
			    sig_chi1,
			    chi2,
			    sig_chi2,
			    chi3,
			    sig_chi3,
			    chi4,
			    sig_chi4);

   short int added=0;
   for (unsigned int i=0; i<typed_rotamers.size(); i++) {
      if (typed_rotamers[i].Type() == restype) {
	 typed_rotamers[i].add_simple_rotamer(rot);
	 added=1;
	 break;
      }
   }
   if (! added) {
      // std::cout  << "adding new restype " << restype << "to typed_rotamers" << std::endl;
      typed_rotamers.push_back(dunbrack_rotamer(restype, rot));
   }
}

// uses above function and wipes out Dunbrack rotamers
// 
void coot::chi_angles::use_richardson_rotamers() {

   typed_rotamers.clear();
   add_richardson_rotamers();
}

void
coot::chi_angles::add_richardson_rotamer(std::string restype,
					 std::string rotamer_name,
					 float anumber,
					 float percent_overall,
					 float percent_alpha,
					 float percent_beta,
					 float percent_other,
					 float chi_1_mode,
					 float chi_1_com, 
					 float chi_2_mode,
					 float chi_2_com, 
					 float chi_3_mode,
					 float chi_3_com, 
					 float chi_4_mode,
					 float chi_4_com) {

   coot::simple_rotamer rot(rotamer_name,
			    percent_overall,
			    percent_alpha,
			    percent_beta,
			    percent_other,
			    chi_1_mode,
			    chi_1_com,
			    chi_2_mode,
			    chi_2_com,
			    chi_3_mode,
			    chi_3_com,
			    chi_4_mode,
			    chi_4_com);

   bool added = 0;
   for (unsigned int i=0; i<typed_rotamers.size(); i++) {
      if (typed_rotamers[i].Type() == restype) {
	 typed_rotamers[i].add_simple_rotamer(rot);
	 added = 1;
	 break;
      }
   }
   if (! added) {
      typed_rotamers.push_back(dunbrack_rotamer(restype, rot));
   }
} 



// Return success status, 
// 0 means success
// 
// 1: we failed because we didn't find the residue type in
// typed_rotamers
// 
// 2: ichi was wrong [ichi is not zero indexed].
// 
// 3: ... 
//
//
// New style function in which the contact_indices are not passed, the
// are generated using geom_p.
// 
std::pair<short int, float>
coot::chi_angles::change_by(int ichi, double diff, coot::protein_geometry* geom_p) {


   bool add_reverse_contacts = 0;
   std::vector<std::vector<int> > contact_indices =
      util::get_contact_indices_from_restraints(residue, geom_p, 1, add_reverse_contacts);

   std::string resname(residue->GetResName());

   // PRO variables, used on the outside.
   mmdb::PPAtom ordered_residue_atoms = 0; 
   int CA_index = -1; // unfound initially (tree base)
   
   if (resname == "PRO") {
      mmdb::PPAtom residue_atoms;
      int nResidueAtoms;
      residue->GetAtomTable(residue_atoms, nResidueAtoms);

      mmdb::Atom *n = 0;
      mmdb::Atom *ca = 0;
      mmdb::Atom *cb = 0;
      mmdb::Atom *cg = 0;
      mmdb::Atom *cd = 0;
      if (nResidueAtoms > 2) {
	 for (int i=0; i<nResidueAtoms; i++) {
	    if (std::string(residue_atoms[i]->name) == " N  ")
	       n = residue_atoms[i];
	    if (std::string(residue_atoms[i]->name) == " CA ")
	       ca = residue_atoms[i];
	    if (std::string(residue_atoms[i]->name) == " CB ")
	       cb = residue_atoms[i];
	    if (std::string(residue_atoms[i]->name) == " CG ")
	       cg = residue_atoms[i];
	    if (std::string(residue_atoms[i]->name) == " CD ")
	       cd = residue_atoms[i];
	 }
	 
	 if (n && ca && cb && cg && cd) {
	    
	    ordered_residue_atoms = new mmdb::PAtom[nResidueAtoms];
	    ordered_residue_atoms[0] = n;
	    ordered_residue_atoms[1] = ca;
	    ordered_residue_atoms[2] = cb;
	    ordered_residue_atoms[3] = cg;
	    ordered_residue_atoms[4] = cd;
	    int atom_count = 5;
	    for (int i=0; i<nResidueAtoms; i++) {
	       if (std::string(residue_atoms[i]->name) == " N  " ||
		   std::string(residue_atoms[i]->name) == " CA " ||
		   std::string(residue_atoms[i]->name) == " CB " ||
		   std::string(residue_atoms[i]->name) == " CG " ||
		   std::string(residue_atoms[i]->name) == " CD " ) {
	       } else {
		  ordered_residue_atoms[atom_count] = residue_atoms[i];
		  atom_count++;
	       }
	    }
	 }
      }

      if (ordered_residue_atoms)
	 contact_indices =
	    util::get_contact_indices_for_PRO_residue(ordered_residue_atoms,
							    nResidueAtoms, geom_p);
   } // end of specific PRO-logic

   
   // Change the coordinates of the data member residue
   //
   // c.f. dunbrack::GetResidue but here we don't create and return a
   // residue, we just tinker with the one we have.
   //

//    if (ichi < 1 || ichi > 4) {
//       std::cout << "unacceptable chi bond: " << ichi << std::endl;
//       return 2;
//    } 

   // short int isuccess = 0;
   std::pair<short int, float> p(0, 0.0);
   mmdb::PPAtom residue_atoms;
   int nResidueAtoms;
   residue->GetAtomTable(residue_atoms, nResidueAtoms);
   if (resname == "PRO")
      residue_atoms = ordered_residue_atoms;
   if (nResidueAtoms == 0) {
      std::cout << " something broken in atom residue selection in ";
      std::cout << "change_by, got 0 atoms" << std::endl;
   } else {

      // Rotatable bonds -> coord indices
      // 
      std::vector<coot::atom_name_pair> atom_name_pairs
	 = atom_name_pair_list(Residue_Type());

      if (atom_name_pairs.size() == 1) { 
	 if (atom_name_pairs[0].atom1 == "empty") {
	    p.first = 1;
	    return p;
	 }
      } 

      coot::atom_spec_t dummy_tree_base_atom;
      dummy_tree_base_atom.int_user_data = -999; // special "don't use" value

      CA_index = 1;
      if (CA_index != -1) { 
	 // i.e. not unset
	 dummy_tree_base_atom.int_user_data = CA_index; // ignored!
	 dummy_tree_base_atom.atom_name = " N  ";
      }
      
      p = change_by_internal(ichi, diff, atom_name_pairs,
			     contact_indices, residue_atoms, nResidueAtoms,
			     dummy_tree_base_atom);
   }
   delete [] ordered_residue_atoms;
   return p;
}


// mmdb::Residue *
// coot::chi_angles::deep_copy_residue(const mmdb::Residue *residue) const {

//    // Horrible casting to mmdb::Residue because GetSeqNum and GetAtomTable
//    // are not const functions.
//    // 
//    mmdb::Residue *rres = new mmdb::Residue;
//    mmdb::Chain   *chain_p = new mmdb::Chain;
//    chain_p->SetChainID(((mmdb::Residue *)residue)->GetChainID());
//    rres->seqNum = ((mmdb::Residue *)residue)->GetSeqNum();
//    strcpy(rres->name, residue->name);

//    mmdb::PPAtom residue_atoms;
//    int nResidueAtoms;
//    ((mmdb::Residue *)residue)->GetAtomTable(residue_atoms, nResidueAtoms);
//    mmdb::Atom *atom_p;
   
//    for(int iat=0; iat<nResidueAtoms; iat++) { 
//       atom_p = new mmdb::Atom;
// //       atom_p->SetCoordinates(residue_atoms[iat]->x,
// // 			     residue_atoms[iat]->y,
// // 			     residue_atoms[iat]->z,
// // 			     residue_atoms[iat]->occupancy,
// // 			     residue_atoms[iat]->tempFactor);
// //       atom_p->SetAtomName(residue_atoms[iat]->name);
// //       strcpy(atom_p->element, residue_atoms[iat]->element);
//       atom_p->Copy(residue_atoms[iat]);
//       int i_add = rres->AddAtom(atom_p);
//    }
//    chain_p->AddResidue(rres);
//    return rres;
// }

// Return success status, 
// 0 means success
// 
// 1: we failed because we didn't find the residue type in
// typed_rotamers
// 
// 2: ichi was wrong [ichi is not zero indexed].
// 
// 3: ... 
//
//
// Old style function in which the contact_indices are passed.  Not used now?
// 
std::pair<short int, float>
coot::chi_angles::change_by(int ichi, double diff,
			    const std::vector<std::vector<int> > &contact_indices) {

   // Change the coordinates of the data member residue
   //
   // c.f. dunbrack::GetResidue but here we don't create and return a
   // residue, we just tinker with the one we have.
   //

//    if (ichi < 1 || ichi > 4) {
//       std::cout << "unacceptable chi bond: " << ichi << std::endl;
//       return 2;
//    } 

   // short int isuccess = 0;
   std::pair<short int, float> p(0, 0.0);
   mmdb::PPAtom residue_atoms;
   int nResidueAtoms;
   residue->GetAtomTable(residue_atoms, nResidueAtoms);
   if (nResidueAtoms == 0) {
      std::cout << " something broken in atom residue selection in ";
      std::cout << "change_by, got 0 atoms" << std::endl;
   } else {

      // Rotatable bonds -> coord indices
      // 
      std::vector<coot::atom_name_pair> atom_name_pairs
	 = atom_name_pair_list(Residue_Type());

      if (atom_name_pairs.size() == 1) { 
	 if (atom_name_pairs[0].atom1 == "empty") {
	    p.first = 1;
	    return p;
	 }
      } 

      coot::atom_spec_t dummy_tree_base_atom;
      dummy_tree_base_atom.int_user_data = -999; // special "don't use" value
      
      p = change_by_internal(ichi, diff, atom_name_pairs,
			     contact_indices, residue_atoms, nResidueAtoms,
			     dummy_tree_base_atom);
   }
   return p;
}

// Use the protein_geometry list to find the rotatable bonds:
// c.f. get_torsion_bonds_atom_pairs in wiggly-ligand
// 
std::pair<short int, float>
coot::chi_angles::change_by(int imol,
			    int ichi, double diff,
			    const std::vector<std::vector<int> > &contact_indices,
			    coot::protein_geometry *pg_p,
			    const coot::atom_spec_t &tree_base_atom,
			    short int include_hydrogen_torsions_flag) {

   std::pair<short int, float> p(1, 0.0);
   
   mmdb::PPAtom residue_atoms;
   int nResidueAtoms;
   residue->GetAtomTable(residue_atoms, nResidueAtoms);
   std::string residue_name = residue->name;
   // filter out CONST torsions when making atom_name_pairs
   std::vector<coot::atom_name_pair> atom_name_pairs = 
      get_torsion_bonds_atom_pairs(residue_name, imol, pg_p, include_hydrogen_torsions_flag);

   if (atom_name_pairs.size() == 0) {
      std::cout << " Sorry, can't find atom rotatable bonds for residue type ";
      std::cout << residue->name << "\n";
   } else { 
      if (nResidueAtoms == 0) {
	 std::cout << " something broken in atom residue selection in ";
	 std::cout << "change_by, got 0 atoms" << std::endl;
      } else {

	p  = change_by_internal(ichi, diff, atom_name_pairs,
				contact_indices, residue_atoms, nResidueAtoms,
				tree_base_atom);
      }
   }
   return p;
}

// Return 0 on success
// Other values for various kinds of failure
//
// the tree base atom spec contains a special "user data" specifying
// (basically) where it came from and whether or not to use the
// tree_base_atom spec;
std::pair<short int, float>
coot::chi_angles::change_by_internal(int ichi,
				     double diff,
				     const std::vector<coot::atom_name_pair> &atom_name_pairs,
				     const std::vector<std::vector<int> > &contact_indices,
				     mmdb::PPAtom residue_atoms_in,
				     int nResidueAtoms,
				     const coot::atom_spec_t &tree_base_atom) {

   std::pair<short int, float> p(0, 0.0);
   
   // Let's make the coordinates:
   // 
   std::vector< ::Cartesian > coords;
   for(int i=0; i<nResidueAtoms; i++) {
      ::Cartesian c(residue_atoms_in[i]->x,
		    residue_atoms_in[i]->y,
		    residue_atoms_in[i]->z);
      coords.push_back(c);
   }
   mmdb::PPAtom residue_atoms = residue_atoms_in;
   std::vector<coot::atom_index_pair> atom_index_pairs
      = get_atom_index_pairs(atom_name_pairs, residue_atoms_in, nResidueAtoms);


   // debugging
   if (false) {
      std::cout << " -----------   pairs ---------------- " << std::endl;
      for(unsigned int i=0; i<atom_index_pairs.size(); i++)
	 std::cout << "pair " << i << ": " << atom_index_pairs[i].index1
		   << " " << atom_name_pairs[i].atom1 << "       "
		   << atom_index_pairs[i].index2 << " "
		   << atom_name_pairs[i].atom2 << "\n";
      std::cout << " -----------  coords ---------------- " << std::endl;
      for(unsigned int i=0; i<coords.size(); i++)
	 std::cout << i << " " << residue_atoms[i]->GetAtomName()
		   << " " << coords[i] << "\n";;
   
      // display contact indices here:
      //
      std::cout << " ------------------ contact indices --------------------------\n" ;
      for (unsigned int ic=0; ic<contact_indices.size(); ic++) {
	 for (unsigned int jc=0; jc<contact_indices[ic].size(); jc++) {
	    std::cout << " contact " << ic << " "
		      << contact_indices[ic][jc] << std::endl;
	 }
      }
   }
   
   //
   int base_index = 0;
   
   if (tree_base_atom.int_user_data != -999) { // i.e. let's use the
					       // atom spec to find
					       // the base atom
      for(int i=0; i<nResidueAtoms; i++) {
	 if (tree_base_atom.atom_name == residue_atoms[i]->name) {
	    base_index = i;
	    // 	    std::cout << "DEBUG:: Using tree based on atom: "
	    // << tree_base_atom.atom_name << " index: " << base_index << std::endl;
	    break;
	 }
      }
   }
   
   Tree tree;
   tree.SetCoords(coords, base_index, contact_indices);
    // std::cout << tree << "\n";
    // std::cout << base_index << "\n";
   // test that the atom indices are sensible?
   //
   // std::cout << "DEBUG:: tree has " << tree.GetNumberOfVertices() << " vertices\n";

   int ibond = ichi - 1;
   if (ibond >= int(atom_index_pairs.size()) ) {
      std::cout << "ERROR: impossible ibond " << ibond << " (wanted rotamer index "
		<< ichi << ") in " << atom_index_pairs.size()
		<< " rotamer atom pairs" << std::endl;
      p.first = 2;
      return p;
   }

   float tors_orig;
   float tors;
   TreeVertex *tv = tree.GetCoord(atom_index_pairs[ibond].index2);
   if (tv->GetNumberOfChildren() > 0) {
      TreeVertex *tvc0 = tv->GetChild(0);
      tors_orig = tvc0->GetParentDihedralAngle();

      // tors += clipper::Util::d2rad(diff);
      tors = clipper::Util::d2rad(diff);
      p.second = tors_orig + tors; // in radians

//       std::cout << "moving coordinates by " << diff
// 		<< " degrees (tors = " << tors <<  ")." << std::endl;

//        std::cout << "rotating about atom indices: "
//  		<< atom_index_pairs[ibond].index2 << " "
//  		<< atom_index_pairs[ibond].index1 << "\n";


// 20090724
      tree.RotateAboutBond(atom_index_pairs[ibond].index2,
 			   atom_index_pairs[ibond].index1, tors);

      std::vector< ::Cartesian > coords_rotatated =
	 tree.GetAllCartesians();
      
      if (int(coords_rotatated.size()) != nResidueAtoms) {
	 std::cout << "disaster in atom selection, trees, dunbrack\n";
	 p.first = 3;
      } else {
	 for (int iat=0; iat<nResidueAtoms; iat++) {
	    if (0) { 
	       std::cout  << " From ("
			  << residue_atoms[iat]->x << ", "
			  << residue_atoms[iat]->y << ", "
			  << residue_atoms[iat]->z << ") to ("
			  << coords_rotatated[iat].get_x() << ", "
			  << coords_rotatated[iat].get_y() << ", "
			  << coords_rotatated[iat].get_z() << ")"
			  << std::endl;
	    }
	    residue_atoms[iat]->x = coords_rotatated[iat].get_x();
	    residue_atoms[iat]->y = coords_rotatated[iat].get_y();
	    residue_atoms[iat]->z = coords_rotatated[iat].get_z();
	 }
      }
   } else {
      std::cout << "WARNING: this vertex " << atom_index_pairs[ibond].index2
		<< " has no children (strangely)\n";
      std::cout << "         bond: " << ibond << " index2: "
		<< atom_index_pairs[ibond].index2 << "\n";
      TreeVertex *tv1 = tree.GetCoord(atom_index_pairs[ibond].index1);
      std::cout << "         tv1 (index: " << atom_index_pairs[ibond].index1
		<< ") has " << tv1->GetNumberOfChildren()
		<< " children \n";
      p.first = 1;
   } 

   if (0) {
      std::cout << "change_by_internal()... returns " << p.first << " "
		<< p.second << std::endl;
   }
   return p;
}

// On failure, return a vector of size one with pair with pair.first = "empty"
// 
std::vector<coot::atom_name_pair>
coot::chi_angles::atom_name_pair_list(const std::string &res_type) const {

   std::vector<atom_name_pair> pairs;
   // set up the not-filled value first:
   pairs.push_back(coot::atom_name_pair("empty", "empty"));
   for(unsigned int i=0; i<typed_rotamers.size(); i++) {
      if (typed_rotamers[i].Type() == res_type) {
	 pairs = typed_rotamers[i].AtomPairs();
	 break;
      } 
   }
   return pairs;
}


std::vector<coot::atom_name_quad>
coot::chi_angles::atom_name_quad_list(const std::string &residue_type) const {

   std::vector<coot::atom_name_quad> quads;
   quads.push_back(coot::atom_name_quad("empty", "empty", "empty", "empty"));
   for(unsigned int i=0; i<typed_rotamers.size(); i++) {
      if (typed_rotamers[i].Type() == residue_type) {
	 quads = typed_rotamers[i].AtomQuads();
	 break;
      } 
   }
   return quads;
} 

std::vector<coot::atom_name_pair>
coot::chi_angles::get_torsion_bonds_atom_pairs(const std::string &monomer_type,
					       int imol,
					       coot::protein_geometry *pg,
					       short int include_hydrogen_torsions_flag) const {

   std::vector<coot::atom_name_pair> atom_pairs;
   std::vector <coot::dict_torsion_restraint_t> monomer_torsions = 
      pg->get_monomer_torsions_from_geometry(monomer_type, imol);

   std::pair<short int, coot::dictionary_residue_restraints_t> r =
      pg->get_monomer_restraints(monomer_type, imol);

   std::string atom1, atom2;
   if (monomer_torsions.size() > 0) { 
      for(unsigned int i=0; i<monomer_torsions.size(); i++) {

	 if (!monomer_torsions[i].is_const()) {

	    atom1 = monomer_torsions[i].atom_id_1();
	    atom2 = monomer_torsions[i].atom_id_4();

	    if ( (!r.second.is_hydrogen(atom1) && !r.second.is_hydrogen(atom2))
		 || include_hydrogen_torsions_flag) {
	       coot::atom_name_pair pair(monomer_torsions[i].atom_id_2_4c(),
					 monomer_torsions[i].atom_id_3_4c());
	       atom_pairs.push_back(pair);
	    }
	 }
      }
   } else {

      std::cout << "WARNING: residue type " << monomer_type << " not found in "
		<< " restraints dictionary list!\n";
      std::cout << "         Are you sure you read the restraints file correctly!?\n";
   }
   return atom_pairs;
} 
//	 if ((monomer_torsions[i].atom_id_1_4c()[0] != 'H') &&
// 	     (monomer_torsions[i].atom_id_4_4c()[0] != 'H')) {

std::vector<coot::atom_index_pair> 
coot::chi_angles::get_atom_index_pairs(const std::vector<coot::atom_name_pair> &atom_name_pairs,
				     const mmdb::PPAtom atoms, int nresatoms) const {

   int i_store_index;
   std::vector<coot::atom_index_pair> index_pairs;

   for (unsigned int ipair=0; ipair<atom_name_pairs.size(); ipair++) {
      i_store_index = -1;
      for(int i=0; i<nresatoms; i++) {
	 std::string atomname = atoms[i]->name;
	 if (atomname == atom_name_pairs[ipair].atom1) {
	    i_store_index = i;
	 }
      }
      if (i_store_index > -1) { // i.e. we found the first atom
	 for(int i2=0; i2<nresatoms; i2++) {
	    std::string atomname = atoms[i2]->name;
	    if (atomname == atom_name_pairs[ipair].atom2) {
	       index_pairs.push_back(coot::atom_index_pair(i_store_index, i2));
	    }
	 }
      } else {
	 std::cout << "first atom :" << atom_name_pairs[ipair].atom1
		   << ": not found in residue\n";
      }
   }

   if (index_pairs.size() != atom_name_pairs.size()) {
      std::cout << "WARNING:: failure to find all atom pair in residue atoms!\n" ;
   } 
   return index_pairs;
} 


std::vector<coot::atom_index_quad>
coot::chi_angles::get_atom_index_quads(const std::vector<coot::atom_name_quad> &atom_name_quads,
				       const mmdb::PPAtom atoms, int nresatoms) const {

   std::vector<coot::atom_index_quad> index_quads;
   for (unsigned int iquad=0; iquad<atom_name_quads.size(); iquad++) {
      int index_1 = -1;
      int index_2 = -1;
      int index_3 = -1;
      int index_4 = -1;
      for (int iat=0; iat<nresatoms; iat++) {
	 std::string atomname = atoms[iat]->name;
	 if (atomname == atom_name_quads[iquad].atom_name(0))
	    index_1 = iat;
	 if (atomname == atom_name_quads[iquad].atom_name(1))
	    index_2 = iat;
	 if (atomname == atom_name_quads[iquad].atom_name(2))
	    index_3 = iat;
	 if (atomname == atom_name_quads[iquad].atom_name(3))
	    index_4 = iat;
      }
      if ((index_1 != -1) && (index_2 != -1) && (index_3 != -1) && (index_4 != -1)) {
	 coot::atom_index_quad iq(index_1, index_2, index_3, index_4);
	 index_quads.push_back(iq);
      } 
   } 

   return index_quads;
}



std::pair<std::string, std::string>
coot::chi_angles::atom_names_of_bond(int i) const {

   std::vector<atom_name_pair> pl = atom_name_pair_list(Residue_Type());
   std::string a1 = "";
   std::string a2 = "";

   if (pl.size() > 1) { // empty is not the first pair.

      if (i < int(pl.size())) {
	 if (i >= 0) { 
	    a1 = pl[i].atom1;
	    a2 = pl[i].atom2;
	 }
      }
   }
   return std::pair<std::string, std::string>(a1, a2);
}


std::vector<std::pair<int,float> > 
coot::chi_angles::get_chi_angles() const {

   std::vector<std::pair <int, float> > v;
   for (unsigned int ir=0; ir<typed_rotamers.size(); ir++) {
      if (typed_rotamers[ir].Type() == Residue_Type()) {
	 v = typed_rotamers[ir].get_chi_angles(residue);
	 break;
      }
   }
   return v;
} 
