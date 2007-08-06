
#include "rotamer.hh"


// Return a flag saying whether we did this or not.
//
// LEU, VAL, THR have "nomenclature" or real chiral centres - they are
// not dealt with here.
// 
// We deal with bifurcated symmetric non-chiral side chains (PHE, ASP,
// GLU, THR)
// 
int 
coot::rotamer::optimize_rotamer_by_atom_names() {

   // Chi2 of PHE, TYR should be between -90 -> +90

   int imoved = 0;
   // p.first        // flag for assigned,                    1 
                     // unassigned due to missing atoms,      0 
                     // unassigned due to rotamer not found. -1
                     // unassigned due to GLY/ALA            -2
   std::pair<short int, double> p_init = probability_of_this_rotamer(); 

   std::string residue_name(residue->GetResName());
   if ((residue_name == "PHE") ||
       (residue_name == "ASP") ||
       (residue_name == "GLU") ||
       (residue_name == "TYR")
       // 
       // These have Chiral Centres and should be treated with more
       // sophistication (or differently)
       // (residue_name == "LEU") ||
       // (residue_name == "VAL")
       ) {

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

      CResidue *residue_copy = deep_copy_residue(residue);
      PPCAtom residue_atoms;
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
//       std::cout << "swapped " << nfound << " atom pairs of residue "
// 		<< residue->GetResName() << " " << residue->GetSeqNum()
// 		<< std::endl;

      coot::rotamer dc(residue_copy);
      std::pair<short int, double> p_swapped = dc.probability_of_this_rotamer();

      // Delete residue_copy
      for (int ifirst=0; ifirst<n_residue_atoms; ifirst++) {
	 residue_copy->DeleteAtom(ifirst);
      }
      delete residue_copy;
      residue_copy = 0;
      
      if (p_init.first && p_swapped.first) { 
// 	 std::cout << "compare: " << p_init.second << "  " << p_swapped.second << "\n";
	 if (p_swapped.second > p_init.second) {
// 	    std::cout << "Fixing residue " << residue->GetSeqNum()
// 		      << " " << residue->GetChainID() << " " << residue->GetResName()
// 		      << "\n";
	    // make the changes to residue then
	    PPCAtom residue_atoms;
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
	    std::cout << "Fixed " << nfound << " atom pairs of residue "
		      << residue->GetResName() << " " << residue->GetSeqNum()
		      << std::endl;
	 }
      } else {
	 std::cout << "badness in atom selection " << std::endl;
      }
   }

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

std::pair<short int, double>
coot::rotamer::probability_of_this_rotamer() {

   short int state;  // flag for assigned,                    1 
                     // unassigned due to missing atoms,      0 
                     // unassigned due to rotamer not found. -1
                     // unassigned due to GLY/ALA            -2
   state = 0;

   // First (and mainly), what are the actual chi angles for (the
   // atoms of) this residue?

   std::string residue_name(residue->GetResName());

//    std::cout << "ROTAMER:: checking residue " << residue->GetSeqNum()
// 	     << residue->GetChainID() << " " << residue_name << std::endl;

   PCAtom *residue_atoms;
   int n_residue_atoms;
   residue->GetAtomTable(residue_atoms, n_residue_atoms);
   std::vector<std::vector<std::string> > rotamer_ats = rotamer_atoms(residue_name);
   if (residue_name == "MSE")
      residue_name = "MET";

   if (rotamer_ats.size() == 0)
      return std::pair<short int, double>(-2, 0.0); // no chi-squared value for residues with no side chain
   
   std::vector<std::vector<int> > atom_indices =
      rotamer_atom_names_to_indices(rotamer_ats, residue_atoms, n_residue_atoms);
   
   if (rotamer_ats.size() != atom_indices.size()) {
      // missing atom(s)
      return std::pair<short int, double>(0, 0.0);
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
   std::pair<short int, double> d = probability_of_this_rotamer(chi_angles, rots);
   
//    std::cout << "DEBUG:: probability_of_this_rotamer returns "
// 	     << d.first << " " << d.second << std::endl;
   
   return d;
}

std::pair<short int, double>
coot::rotamer::probability_of_this_rotamer(const std::vector<double> &chi_angle_values,
					    const std::vector<coot::simple_rotamer> &rots) const {

   short int state;  // flag for assigned,                    1 
                     // unassigned due to missing atoms,      0 
                     // unassigned due to rotamer not found. -1
                     // unassigned due to GLY/ALA            -2
   double d = 0; // initially unassigned.
   state = -1; // rotamer not found initially

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
	    // std::cout  << "DEBUG:: comparing chi " << ich << " " << target <<  " "
	    // << model << std::endl;
	    if (similar_rotamer_chi(target, model))
	       isimilar++;
	    else
	       break;
	 }
//  	 std::cout << "DEBUG:: rotamer " << irot << " isimilar = "
//  		   << isimilar << " out of " << chi_angle_values.size() << std::endl;
	 if (isimilar == int(chi_angle_values.size())) {
	    d = rots[irot].P_r1234();
	    if (d > 0) 
	       state = 1;
	 }
	 if (state == 1)
	    break;
      }
   }
      
   return std::pair<short int, double>(state, d);
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
					      PCAtom *residue_atoms,
					      int n_residue_atoms) const {

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
			    PCAtom *residue_atoms) {

   double tors = 0.0;

   std::vector<clipper::Coord_orth> a;
   CAtom *at;

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

