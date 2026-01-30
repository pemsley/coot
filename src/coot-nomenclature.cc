/* src/main.cc
 * 
 * Copyright 2005 by The University of York
 * Copyright 2015 by Medical Research Council
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */


#include <iostream>
#include "coot-nomenclature.hh"
#include "ideal/simple-restraint.hh"
#include "mmdb2/mmdb_mattype.h"


#ifdef USE_DUNBRACK_ROTAMERS
#include "ligand/dunbrack.hh"
#else 
#include "ligand/richardson-rotamer.hh"
#endif 


// Here we rename atoms to fix nomeclature errors. Note ILEs are not fixed
// by renaming atoms.
// 
std::vector<mmdb::Residue *>
coot::nomenclature::fix(const coot::protein_geometry *Geom_p) {

   std::vector<mmdb::Residue *> vr = fix_and_swap_maybe(Geom_p, 1);
   return vr;
}

// just list the swaps needed - don't apply them.
// 
std::vector<mmdb::Residue *>
coot::nomenclature::list(const coot::protein_geometry *Geom_p) {

   std::vector<mmdb::Residue *> vr = fix_and_swap_maybe(Geom_p, 0);
   return vr;
}


// Here we rename atoms to fix nomeclature errors. Note ILEs are not fixed
// by renaming atoms.
// 
std::vector<mmdb::Residue *>
coot::nomenclature::fix_and_swap_maybe(const coot::protein_geometry *Geom_p, bool apply_swaps) {

   std::vector<mmdb::Residue *> vr;
   if (mol_) { 

      int n_models = mol_->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) { 
      
	 mmdb::Model *model_p = mol_->GetModel(imod);
	 if (model_p) { 
	    mmdb::Chain *chain_p;
	    // run over chains of the existing mol
	    int nchains = model_p->GetNumberOfChains();
	    if (nchains <= 0) { 
	       std::cout << "bad nchains in molecule " << nchains
			 << std::endl;
	    } else { 
	       for (int ichain=0; ichain<nchains; ichain++) {
		  chain_p = model_p->GetChain(ichain);
		  if (chain_p == NULL) {  
		     // This should not be necessary. It seem to be a
		     // result of mmdb corruption elsewhere - possibly
		     // DeleteChain in update_molecule_to().
		     std::cout << "NULL chain in fix_nomenclature_errors "
			       << std::endl;
		  } else { 
		     int nres = chain_p->GetNumberOfResidues();
		     mmdb::PResidue residue_p;
		     for (int ires=0; ires<nres; ires++) { 
			residue_p = chain_p->GetResidue(ires);
			std::string residue_name(residue_p->GetResName());

			if ((residue_name == "PHE") ||
			    (residue_name == "TYR")) {

			   // if apply_swaps is true, apply swaps if they are found - otherwise just
			   // return a flag saying that they should be swapped.
			   int isw = test_and_fix_PHE_TYR_nomenclature_errors(residue_p, apply_swaps);
			   if (isw) {
			      if (apply_swaps) // don't tell us before swap is made
				 std::cout << "INFO:: (result) " << residue_name << " swapped atoms in "
					   << coot::residue_spec_t(residue_p)
					   << " " << residue_p->GetResName() << std::endl;
			      vr.push_back(residue_p);
			   }
			}
		     
			if ((residue_name == "ASP") ||
			    (residue_name == "GLU")) {

			   int isw = test_and_fix_ASP_GLU_nomenclature_errors(residue_p, apply_swaps);
			   if (isw) {
			      if (apply_swaps)
				 std::cout << "INFO:: (result) " << residue_name << " swapped atoms in "
					   << coot::residue_spec_t(residue_p)
					   << " " << residue_p->GetResName() << std::endl;
			      vr.push_back(residue_p);
			   }
			}

			if (residue_name == "THR") { 

			   // This is assigned a sign in the refmac
			   // dictionary, (unlike LEU and VAL).

			   int imol = 0; // dummy
			   std::vector<coot::dict_chiral_restraint_t> chiral_restraints = 
			      Geom_p->get_monomer_chiral_volumes(std::string(residue_p->name), imol);
			   coot::dict_chiral_restraint_t chiral_restraint;
			   for (unsigned int irestr=0; irestr<chiral_restraints.size(); irestr++) { 
			      chiral_restraint = chiral_restraints[irestr];
			   
			      // now what is the chiral centre for that restraint?
			    
			      if (chiral_restraint.atom_id_c_4c() == " CB ") {
			      
				 // make a synthetic restraint from
				 // chiral_restraint because we want to
				 // impose a sign on the CB chiral
				 // centre - in the chiral_restraint the 

				 // 			      std::cout << "For residue " << residue_p->GetSeqNum() << " "
				 // 					<< residue_p->GetResName() << " chiral centre: "
				 // 					<< " for restraint: " << irestr << " :" 
				 // 					<<  chiral_restraint.atom_id_c_4c() << ":\n";

				 std::vector<std::pair<short int, coot::atom_spec_t> > c = 
				    coot::is_inverted_chiral_atom_p(chiral_restraint, residue_p);
				 for (unsigned int ibad=0; ibad<c.size(); ibad++) { 
				    if (c[ibad].first) {
				       std::cout << "INFO:: found bad THR chiral atom: " 
						 << chain_p->GetChainID() << " " 
						 << residue_p->GetSeqNum() << " "
						 << residue_p->GetInsCode() << " "
						 << chiral_restraint.atom_id_c_4c() << " "
						 << c[ibad].second.alt_conf << std::endl;

				       // swap the OG1 and CG2 atoms of
				       // the residue for the given
				       // alt_conf
				       std::string alt_conf_bad = c[ibad].second.alt_conf;
				       mmdb::PPAtom residue_atoms;
				       int n_residue_atoms;
				       residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
				       mmdb::Atom *og1 = 0;
				       mmdb::Atom *cg2 = 0;
				       for (int iat=0; iat<n_residue_atoms; iat++) {
					  std::string alt_conf = residue_atoms[iat]->altLoc;
					  std::string atom_name = residue_atoms[iat]->name;
					  if (atom_name == " OG1" )
					     if (alt_conf == alt_conf_bad)
						og1 = residue_atoms[iat];
					  if (atom_name == " CG2" )
					     if (alt_conf == alt_conf_bad)
						cg2 = residue_atoms[iat];
				       }
				       if (og1 && cg2) {

					  if (apply_swaps) { 
					     og1->SetAtomName(" CG2");
					     cg2->SetAtomName(" OG1");
					     std::cout << "        CG2 and OG1 atoms swapped\n";
					     std::cout << "INFO:: swapped atoms in "
						       << coot::residue_spec_t(residue_p)
						       << " " << residue_p->GetResName() << std::endl;
					  }
					  vr.push_back(residue_p);
				       } else {
					  // This can't happen:
					  std::cout << "ERROR:: Bizarre missing atom scenario "
						    << "in fix_nomenclature_errors\n";
				       }
				    }
				 }
			      }
			   }
			}
		     
			if ((residue_name == "LEU") ||
			    (residue_name == "VAL")) {
			
			   int volume_sign = -1;
			   coot::dict_chiral_restraint_t synthetic_restraint;
			
			   if (residue_name == "VAL")
			      synthetic_restraint =
				 coot::dict_chiral_restraint_t(residue_name,
							       " CB ", " CA ", " CG1", " CG2",
							       volume_sign);
			   if (residue_name == "LEU")
			      synthetic_restraint =
				 coot::dict_chiral_restraint_t(residue_name,
							       " CG ", " CB ", " CD1", " CD2",
							       volume_sign);
			   
			   std::vector<std::pair<short int, coot::atom_spec_t> > c = 
			      coot::is_inverted_chiral_atom_p(synthetic_restraint,
							      residue_p);
			
			   for (unsigned int ibad=0; ibad<c.size(); ibad++) {
			      if (c[ibad].first) {
				 std::cout << "INFO:: found bad " << residue_name
					   << " chiral atom: " 
					   << chain_p->GetChainID() << " " 
					   << residue_p->GetSeqNum() << " "
					   << residue_p->GetInsCode() << " "
					   << synthetic_restraint.atom_id_c_4c() << " "
					   << c[ibad].second.alt_conf << std::endl;

				 // swap the CG1 and CG2 atoms of
				 // the residue for the given
				 // alt_conf
				 std::string alt_conf_bad = c[ibad].second.alt_conf;
				 std::string target_atom_1 = " CG1";
				 std::string target_atom_2 = " CG2";
				 if (residue_name == "LEU") {
				    target_atom_1 = " CD1";
				    target_atom_2 = " CD2";
				 }
				 mmdb::PPAtom residue_atoms;
				 int n_residue_atoms;
				 residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
				 mmdb::Atom *cg1 = 0; // cd1 and cd2 for LEU of course
				 mmdb::Atom *cg2 = 0;
				 for (int iat=0; iat<n_residue_atoms; iat++) {
				    std::string alt_conf = residue_atoms[iat]->altLoc;
				    std::string atom_name = residue_atoms[iat]->name;
				    if (atom_name == target_atom_1 )
				       if (alt_conf == alt_conf_bad)
					  cg1 = residue_atoms[iat];
				    if (atom_name == target_atom_2 )
				       if (alt_conf == alt_conf_bad)
					  cg2 = residue_atoms[iat];
				 }
				 if (cg1 && cg2) {
				    if (apply_swaps) { 
				       cg1->SetAtomName(target_atom_2.c_str());
				       cg2->SetAtomName(target_atom_1.c_str());
				       std::cout << "        " << target_atom_1 << " and "
						 << target_atom_2 << " atoms swapped\n";
				       std::cout << "INFO:: swapped atoms in "
						 << coot::residue_spec_t(residue_p)
						 << " " << residue_p->GetResName() << std::endl;
				    }
				    vr.push_back(residue_p);
				 } else {
				    // This can't happen:
				    std::cout << "ERROR:: Bizarre missing atom scenario "
					      << "in fix_nomenclature_errors, residue type:"
					      << residue_name << "\n";
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   if (apply_swaps)
      std::cout << "INFO:: " << vr.size() << " residues had their atoms swapped\n";
   return vr;
}


// test chi2 (and potentially fix) so that chi2 is the range -90 -> +90.
int
coot::nomenclature::test_and_fix_PHE_TYR_nomenclature_errors(mmdb::Residue *residue_p,
							     bool apply_swap_if_found) {

   // return 0; // testing memory corruption

   // PDBv3 FIXME - all this function.

   int iswapped = 0; // number of alt confs swapped in this residue
   mmdb::Atom **residue_atoms = 0;
   int n_residue_atoms = 0;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);

   // We have to test all the altconfs, so what are the altconfs of the CD1 atoms?
   
   std::vector<std::string> alt_conf_list;
   // first get the altconfs in the residue:
   for (int i=0; i<n_residue_atoms; i++) {
      std::string atom_name = residue_atoms[i]->name;
      if(atom_name == " CD1") {
	 alt_conf_list.push_back(residue_atoms[i]->altLoc);
      }
   }

   // so we need to test (and potentially fix) each altconf
   //
   // We first try to find all chi2 atoms with the same altconf as the
   // CD1 atom, failing that we find the atom with alt conf "".  And
   // failing that, there are missing atoms and we can't determine the
   // torsion
   for (unsigned int ialtconf=0; ialtconf<alt_conf_list.size(); ialtconf++) {
      mmdb::Atom *CA  = 0;
      mmdb::Atom *CB  = 0;
      mmdb::Atom *CG  = 0;
      mmdb::Atom *CD1 = 0;
      mmdb::Atom *CD2 = 0;
      for (int i=0; i<n_residue_atoms; i++) {
	 std::string atom_name = residue_atoms[i]->name;
	 std::string atom_altconf = residue_atoms[i]->altLoc;
	 if (atom_altconf == alt_conf_list[ialtconf]) {
	    if (atom_name == " CA ")
	       CA = residue_atoms[i];
	    if (atom_name == " CB ")
	       CB = residue_atoms[i];
	    if (atom_name == " CG ")
	       CG = residue_atoms[i];
	    if (atom_name == " CD1")
	       CD1 = residue_atoms[i];
	    if (atom_name == " CD2")
	       CD2 = residue_atoms[i];
	 }
      }
      if (CA==0 || CB==0 || CG==0) { // no need for CD1, it will be set
	 for (int i=0; i<n_residue_atoms; i++) {
	    std::string atom_name    = residue_atoms[i]->name;
	    std::string atom_altconf = residue_atoms[i]->altLoc;
	    if (atom_altconf == "") {
	       if (atom_name == " CA ")
		  CA = residue_atoms[i];
	       if (atom_name == " CB ")
		  CB = residue_atoms[i];
	       if (atom_name == " CG ")
		  CG = residue_atoms[i];
	       if (atom_name == " CD2")
		  CD2 = residue_atoms[i];
	    }
	 }
      }
//       std::cout << "DEBUG:: Pointers: " << CA << " " << CB << " " << CG
// 		<< " " << CD1 << std::endl;
      
      if (CA && CB && CG && CD1) {
	 
	 clipper::Coord_orth a1(CA->x,  CA->y,  CA->z);
	 clipper::Coord_orth a2(CB->x,  CB->y,  CB->z);
	 clipper::Coord_orth a3(CG->x,  CG->y,  CG->z);
	 clipper::Coord_orth a4(CD1->x, CD1->y, CD1->z);

	 double tors = clipper::Util::rad2d(clipper::Coord_orth::torsion(a1, a2, a3, a4));

	 // std::cout << "DEBUG:: CD1 torsion: " << tors << std::endl;
	 if (tors < -90.0 || tors > 90.0) {
	    // ooops, there was a problem with this torsion.

	    if (CD2) {
	       clipper::Coord_orth a4_o(CD2->x, CD2->y, CD2->z);
	       double cg2_tors =
		  clipper::Util::rad2d(clipper::Coord_orth::torsion(a1, a2, a3, a4_o));
	       // if cg2_tors is in range, then we swap atom names
	       if (cg2_tors > -90.0 && cg2_tors < 90.0) {
		  // find CE1 and CE2 and swap if both sets exists.
		  mmdb::Atom *CE1 = 0;
		  mmdb::Atom *CE2 = 0;
		  mmdb::Atom *HD1 = 0;
		  mmdb::Atom *HD2 = 0;
		  mmdb::Atom *HE1 = 0;
		  mmdb::Atom *HE2 = 0;
		  for (int ie=0; ie<n_residue_atoms; ie++) {
		     std::string e_atom_name = residue_atoms[ie]->name;
		     std::string e_atom_altconf = residue_atoms[ie]->altLoc;
		     if (e_atom_altconf == alt_conf_list[ialtconf]) {
			if (e_atom_name == " CE1") 
			   CE1 = residue_atoms[ie];
			if (e_atom_name == " CE2") 
			   CE2 = residue_atoms[ie];
			
			if (e_atom_name == " HD1")
			   HD1 = residue_atoms[ie];
			if (e_atom_name == " HD2") 
			   HD2 = residue_atoms[ie];
			if (e_atom_name == " HE1") 
			   HE1 = residue_atoms[ie];
			if (e_atom_name == " HE2") 
			   HE2 = residue_atoms[ie];
		     }
		  }
		  if (CE1 && CE2) {
		     if (apply_swap_if_found) { 
#if 0
			CD1->SetAtomName(" CD2");
			CD2->SetAtomName(" CD1");
			CE1->SetAtomName(" CE2");
			CE2->SetAtomName(" CE1");
#endif
                        mmdb::realtype pos_cd1[3] = {CD1->x, CD1->y, CD1->z};
                        mmdb::realtype pos_cd2[3] = {CD2->x, CD2->y, CD2->z};
                        mmdb::realtype pos_ce1[3] = {CE1->x, CE1->y, CE1->z};
                        mmdb::realtype pos_ce2[3] = {CE2->x, CE2->y, CE2->z};

                        CD1->x = pos_cd2[0]; CD1->y = pos_cd2[1]; CD1->z = pos_cd2[2];
                        CD2->x = pos_cd1[0]; CD2->y = pos_cd1[1]; CD2->z = pos_cd1[2];
                        CE1->x = pos_ce2[0]; CE1->y = pos_ce2[1]; CE1->z = pos_ce2[2];
                        CE2->x = pos_ce1[0]; CE2->y = pos_ce1[1]; CE2->z = pos_ce1[2];
		     }
		     if (false)
			std::cout << "DEBUG:: swapped in test_and_fix_PHE_TYR_nomenclature_errors()"
				  << std::endl;

		     iswapped++;
		  }

		  if (HD1 && HD2) {
		     if (apply_swap_if_found) { 
			HD1->SetAtomName(" HD2");
			HD2->SetAtomName(" HD1");
		     }
		  } 
		  if (HE1 && HE2) {
		     if (apply_swap_if_found) { 
			HE1->SetAtomName(" HE2");
			HE2->SetAtomName(" HE1");
		     }
		  }
	       }
	    }
	 }
						    
      } else {
	 std::cout << "WARNING:: PHE/TYR nomenc: missing atoms in " << residue_p->GetChainID() << " "
		   << residue_p->GetSeqNum() << " " << residue_p->GetResName()
		   << std::endl;
      }
   }
   return iswapped;
}

// ASP:: test chi2 (and potentially fix) so that chi2 is the range -90 -> +90.
// GLU:: test chi3 (and potentially fix) so that chi3 is the range -90 -> +90.
// 
int
coot::nomenclature::test_and_fix_ASP_GLU_nomenclature_errors(mmdb::Residue *residue_p,
							     bool apply_swap_if_found) {

   // return 0;  // testing memory corruption

   int iswapped = 0;
   mmdb::PPAtom residue_atoms;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);

   std::string residue_name = residue_p->GetResName();

   if (residue_name == "ASP" || residue_name == "GLU") {

      std::string test_atom_name = " OD1";
      if (residue_name == "GLU")
	 test_atom_name = " OE1";

      std::vector<std::string> alt_conf_list;
      // first get the altconfs in the residue:
      for (int i=0; i<n_residue_atoms; i++) {
	 std::string atom_name = residue_atoms[i]->name;
	 if(atom_name == test_atom_name) {
	    alt_conf_list.push_back(residue_atoms[i]->altLoc);
	 }
      }

      // OK, so now we have a list of all alt confs.

      for (unsigned int ialtconf=0; ialtconf<alt_conf_list.size(); ialtconf++) {

	 coot::atom_index_quad quad;
	 
	 for (int i=0; i<n_residue_atoms; i++) {
	    std::string atom_name = residue_atoms[i]->name;
	    std::string atom_altconf = residue_atoms[i]->altLoc;
	    if (atom_altconf == alt_conf_list[ialtconf]) {

	       if (residue_name == "ASP") {
		  if (atom_name == " CA ")
		     quad.index1 = i;
		  if (atom_name == " CB ")
		     quad.index2 = i;
		  if (atom_name == " CG ")
		     quad.index3 = i;
		  if (atom_name == " OD1")
		     quad.index4 = i;
	       }

	       if (residue_name == "GLU") {
		  if (atom_name == " CB ")
		     quad.index1 = i;
		  if (atom_name == " CG ")
		     quad.index2 = i;
		  if (atom_name == " CD ")
		     quad.index3 = i;
		  if (atom_name == " OE1")
		     quad.index4 = i;
	       }
	    }
	 }

	 try {
	    double torsion = quad.torsion(residue_p);
	    bool torsion_is_good = false;
	    if (torsion >= -90.0) {
	       if (torsion <= 90.0) {
		  torsion_is_good = true;
	       }
	    }
	    if (! torsion_is_good) {

	       // SWAP names

	       std::string swap_name_1 = " OD1";
	       std::string swap_name_2 = " OD2";
	       if (residue_name == "GLU") {
		  swap_name_1 = " OE1";
		  swap_name_2 = " OE2";
	       }

	       mmdb::Atom *at_1 = 0;
	       mmdb::Atom *at_2 = 0;
	       for (int i=0; i<n_residue_atoms; i++) {
		  std::string atom_name = residue_atoms[i]->name;
		  std::string atom_altconf = residue_atoms[i]->altLoc;
		  if (atom_altconf == alt_conf_list[ialtconf]) {
		     if (atom_name == swap_name_1)
			at_1 = residue_atoms[i];
		     if (atom_name == swap_name_2)
			at_2 = residue_atoms[i];
		  }
	       }
	       if (at_1 && at_2) {
		  if (apply_swap_if_found) {
		     at_1->SetAtomName(swap_name_2.c_str());
		     at_2->SetAtomName(swap_name_1.c_str());
		  }
		  iswapped = 1; // either they have been swapped or they need to be.
	       }
	    } 
	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << "WARNING:: missing atoms " << rte.what() << std::endl;
	 } 
      }
   } 
   return iswapped;
} 
