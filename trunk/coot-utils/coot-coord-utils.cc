/* coot-utils/coot-coord-utils.cc
 * 
 * Copyright 2006 by The University of York
 * Copyright 2007 by Paul Emsley
 * Copyright 2007, 2008, 2009 by The University of Oxford
 * Author: Paul Emsley
 * Author: Bernhard Lohkamp
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

#include <algorithm>
#include <stdexcept>

#include <string.h> // for strcpy
#include "coot-utils.hh"
#include "coot-coord-utils.hh"
#include "mmdb_tables.h"  // for Get1LetterCode()
#include "mmdb_graph.h" // for graph matching

#include "coot-sysdep.h"

#include "clipper/mmdb/clipper_mmdb.h"

std::vector<std::string>
coot::util::residue_types_in_molecule(CMMDBManager *mol) { 

   std::vector<std::string> v;

   if (mol) { 

      int n_models = mol->GetNumberOfModels();
      
      for (int imod=1; imod<=n_models; imod++) { 
      
	 CModel *model_p = mol->GetModel(imod);
   
	 CChain *chain;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 if (nchains <= 0) { 
	    std::cout << "bad nchains in trim molecule " << nchains
		      << std::endl;
	 } else { 
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain = model_p->GetChain(ichain);
	       if (chain == NULL) {  
		  // This should not be necessary. It seem to be a
		  // result of mmdb corruption elsewhere - possibly
		  // DeleteChain in update_molecule_to().
		  std::cout << "NULL chain in residues_in_molecule: "
			    << std::endl;
	       } else { 
		  int nres = chain->GetNumberOfResidues();
		  CResidue *residue_p;
		  for (int ires=0; ires<nres; ires++) { 
		     residue_p = chain->GetResidue(ires);

//		     int n_atoms = residue_p->GetNumberOfAtoms();
// 		     for (int iat=0; iat<n_atoms; iat++) {
// 			CAtom *atom_p = residue_p->GetAtom(iat);
// 		     }

		     std::string resname = residue_p->name;

		     if (! is_member_p(v, resname)) { 
			v.push_back(resname);
		     } 
		  }
	       }
	    }
	 }
      }
   }
   return v;
}

std::vector<std::string>
coot::util::non_standard_residue_types_in_molecule(CMMDBManager *mol) {

   std::vector<std::string> r;
   std::vector<std::string> v = residue_types_in_molecule(mol);
   std::vector<std::string> standards = coot::util::standard_residue_types();
   
   for (unsigned int i=0; i<v.size(); i++)
      if (! is_member_p(standards, v[i]))
	 r.push_back(v[i]);
	    
   return r; 
}

std::vector<std::string>
coot::util::standard_residue_types() {

   std::vector<std::string> v;
   v.push_back("ALA"); v.push_back("ARG"); v.push_back("ASP");
   v.push_back("ASN"); v.push_back("CYS"); v.push_back("SER");
   v.push_back("PRO"); v.push_back("PHE"); v.push_back("GLY");
   v.push_back("GLU"); v.push_back("GLN"); v.push_back("ILE");
   v.push_back("LEU"); v.push_back("TYR"); v.push_back("TRP");
   v.push_back("HIS"); v.push_back("LYS"); v.push_back("MET");
   v.push_back("VAL"); v.push_back("THR"); v.push_back("MSE"); 
   return v;
}

std::vector<std::string>
coot::util::PDB_standard_residue_types() {

   std::vector<std::string> v = coot::util::standard_residue_types();
   v.push_back("Td"); v.push_back("Tr"); v.push_back("T");
   v.push_back("Gd"); v.push_back("Gr"); v.push_back("G");
   v.push_back("Ad"); v.push_back("Ar"); v.push_back("A");

   v.push_back("DG"); v.push_back("DC"); v.push_back("DA");
   v.push_back("DU"); v.push_back("DT"); v.push_back("DI");

   v.push_back("UNK"); v.push_back("N");

   return v;
}




std::vector<std::pair<int, int> >
coot::util::pair_residue_atoms(CResidue *a_residue_p,
			       CResidue *b_residue_p) {

   std::vector<std::pair<int, int> > pv;

   PPCAtom residue_atoms_1 = NULL;
   PPCAtom residue_atoms_2 = NULL;
   int n_residue_atoms_1, n_residue_atoms_2;
   a_residue_p->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
   b_residue_p->GetAtomTable(residue_atoms_2, n_residue_atoms_2);

   for (int i=0; i<n_residue_atoms_1; i++) {
      std::string atn1(residue_atoms_1[i]->name);
      std::string alt1(residue_atoms_1[i]->altLoc);
      for (int j=0; j<n_residue_atoms_2; j++) {
	 std::string atn2(residue_atoms_2[j]->name);
	 std::string alt2(residue_atoms_2[j]->altLoc);
	 if (atn1 == atn2) {
	    if (alt1 == alt2) {
	       std::pair<int, int> p(i,j);
	       pv.push_back(p);
	       break;
	    }
	 }
      }
   }
   return pv;
}

CMMDBManager *
coot::util::copy_molecule(CMMDBManager *mol) {
   CMMDBManager *n = new CMMDBManager;
   n->Copy(mol, MMDBFCM_All);
   return n;
}

void
coot::util::translate_close_to_origin(CMMDBManager *mol) {

   try {
      std::pair<clipper::Cell, clipper::Spacegroup> csp = get_cell_symm(mol);
      clipper::Coord_frac cf = coot::util::shift_to_origin(mol);
      clipper::Coord_orth co = cf.coord_orth(csp.first);
      for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
	 int imod = 1;
	 CModel *model_p = mol->GetModel(imod);
	 CChain *chain_p;
	 int nchains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<nchains; ichain++) {
	    chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    CResidue *residue_p;
	    CAtom *at;
	    for (int ires=0; ires<nres; ires++) { 
	       residue_p = chain_p->GetResidue(ires);
	       int n_atoms = residue_p->GetNumberOfAtoms();
	       for (int iat=0; iat<n_atoms; iat++) {
		  at = residue_p->GetAtom(iat);
		  at->x += co.x();
		  at->y += co.y();
		  at->z += co.z();
	       }
	    }
	 }
      }
   } 
   catch (std::runtime_error rte) {
      std::cout << rte.what() << std::endl;
   } 
}


void
coot::sort_chains(CMMDBManager *mol) {

#ifdef MMDB_MAJOR_VERSION
#if ((MMDB_MAJOR_VERSION > 1) || ((MMDB_MAJOR_VERSION == 1) && (MMDB_MINOR_VERSION >= 22)))

   for (int imod=1; imod<=mol->GetNumberOfModels(); imod++) {
      CModel *model_p = mol->GetModel(imod);
      model_p->SortChains(SORT_CHAIN_ChainID_Asc); // "B" comes after "A"
   }
   mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
   mol->FinishStructEdit();
#endif   
#else
   
   // for (int imod=1; imod<=mol->GetNumberOfModels(); imod++) {
   int imod = 1;

   {
      CModel *model_p = mol->GetModel(imod);
      CChain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      std::vector<std::pair<CChain *, std::string> > chain_ids(nchains);
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 std::string chain_id = chain_p->GetChainID();
	 chain_ids[ichain] = std::pair<CChain *, std::string> (chain_p, chain_id);
      }
      // now chain_ids is full
      std::sort(chain_ids.begin(), chain_ids.end(), sort_chains_util);

      if (0) 
	 for (int ichain=0; ichain<nchains; ichain++)
	    std::cout << " Sorted chain order " << ichain << " "
		      << chain_ids[ichain].second << std::endl;
      
      CModel *new_model_p = new CModel;
      mol->AddModel(new_model_p);
      int new_model_number = new_model_p->GetSerNum();
      std::cout << "new model number : " << new_model_number << std::endl;
      for (int ichain=0; ichain<nchains; ichain++) {
	 CChain *new_chain_p = new CChain;
	 new_chain_p->Copy(chain_ids[ichain].first);
	 new_chain_p->SetChainID(chain_ids[ichain].second.c_str());
	 new_model_p->AddChain(new_chain_p);
	 std::cout << " adding new chain " << new_chain_p->GetChainID() << std::endl;
      }


      mol->SwapModels(imod, new_model_number);
      // Now delete the old model.
      mol->DeleteModel(2);
      
   }
   mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
   mol->FinishStructEdit();
#endif
   
}
      

bool
coot::sort_chains_util(const std::pair<CChain *, std::string> &a,
		       const std::pair<CChain *, std::string> &b) {

   return (a.second < b.second); 
}

// return residue specs for residues that have atoms that are
// closer than radius Angstroems to any atom in the residue
// specified by res_in.
// 
std::vector<coot::residue_spec_t>
coot::residues_near_residue(const coot::residue_spec_t &rs,
			    CMMDBManager *mol,
			    float radius) {

   std::vector<coot::residue_spec_t> r;

   CResidue *res_p = coot::util::get_residue(rs.chain, rs.resno, rs.insertion_code, mol);
   if (!res_p) {
      std::cout << "OOps failed to find " << rs << " in molecule\n";
   } else {

      // std::cout << "  Finding contacts of " << rs << " in molecule\n";

      std::vector<CResidue *> close_residues = residues_near_residue(res_p, mol, radius);
      
      for (unsigned int i=0; i<close_residues.size(); i++)
	 r.push_back(close_residues[i]);

   }
   return r;
}

std::vector<CResidue *>
coot::residues_near_residue(CResidue *res_ref,
			    CMMDBManager *mol,
			    float radius) {

   PPCAtom atom_selection;
   int n_selected_atoms;
   std::vector<CResidue *> close_residues;

   int SelectionHandle = mol->NewSelection();
   mol->SelectAtoms (SelectionHandle, 0, "*",
		     ANY_RES, // starting resno, an int
		     "*", // any insertion code
		     ANY_RES, // ending resno
		     "*", // ending insertion code
		     "*", // any residue name
		     "*", // atom name
		     "*", // elements
		     "*"  // alt loc.
		     );
      
   mol->GetSelIndex(SelectionHandle, atom_selection, n_selected_atoms);

   if (! res_ref)
      return close_residues;
   
	 
   PPCAtom res_atom_selection = 0;
   int n_res_atoms = 0;
   res_ref->GetAtomTable(res_atom_selection, n_res_atoms);


   if (n_res_atoms > 0) {

      PSContact pscontact = NULL;
      int n_contacts;
      float min_dist = 0.01;
      long i_contact_group = 1;
      mat44 my_matt;
      CSymOps symm;
      for (int i=0; i<4; i++) 
	 for (int j=0; j<4; j++) 
	    my_matt[i][j] = 0.0;      
      for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

      // by-pass bug noted below (?)
      if (min_dist < radius)
	 min_dist = 0.0;
      // bypass a bug when the n_contacts is 11m or so, but
      // pscontact[0] is NULL.  hence crash when we try to get the
      // residue from that atom index.
      // 
      if (radius > 0.0) {

	 mol->SeekContacts(res_atom_selection, n_res_atoms, 
			   atom_selection, n_selected_atoms,
			   min_dist, radius, // min, max distances
			   0,        // seqDist 0 -> in same res also
			   pscontact, n_contacts,
			   0, &my_matt, i_contact_group);
	 
	 //       std::cout << " Found " << n_contacts  << " contacts " << std::endl;
	 
	 if (n_contacts > 0) {
	    if (pscontact) { 
	       int n_cont_diff = 0; 
	       int n_cont_same = 0; 
	       for (int i=0; i<n_contacts; i++) {
		  // 	    std::cout << "   comparing " << atom_selection[pscontact[i].id2]
		  // 		      << " " << coot::atom_spec_t(atom_selection[pscontact[i].id2])
		  // 		      << " " << " to " << rs << " " << res_p << std::endl;
		  if (atom_selection[pscontact[i].id2]->residue != res_ref) {
		     n_cont_diff++;
		     std::vector<CResidue *>::iterator result =
			std::find(close_residues.begin(),
				  close_residues.end(),
				  atom_selection[pscontact[i].id2]->residue);
		  
		     if (result == close_residues.end()) { 
			close_residues.push_back(atom_selection[pscontact[i].id2]->residue);
		     }
		  } else {
		     n_cont_same++;
		  }
	       }
	    } else {
	       std::cout << "ERROR:: trapped null pscontact in residues_near_residue"
			 << std::endl;
	    } 
	 }
      }
      mol->DeleteSelection(SelectionHandle);
   }
   return close_residues;
}

std::vector<CResidue *>
coot::residues_near_position(const clipper::Coord_orth &pt,
			     CMMDBManager *mol,
			     double radius) {

   std::vector<CResidue *> v;
   int imod = 1;
   CModel *model_p = mol->GetModel(imod);
   CChain *chain_p;
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      PCResidue residue_p;
      CAtom *at;
      for (int ires=0; ires<nres; ires++) { 
	 residue_p = chain_p->GetResidue(ires);
	 int n_atoms = residue_p->GetNumberOfAtoms();
	 for (int iat=0; iat<n_atoms; iat++) {
	    at = residue_p->GetAtom(iat);
	    clipper::Coord_orth at_pt(at->x, at->y, at->z);
	    double d = clipper::Coord_orth::length(pt, at_pt);
	    if (d < radius) { 
	       v.push_back(residue_p);
	       break;
	    }
	 }
      }
   }
   return v;
}

// Don't include residues that are HOH residues that are not bonded to
// the protein (if bonding to res_ref but not protein, then reject.
// Reject waters that are not within water_dist_max to any atom in
// res_ref.
//
std::vector<CResidue *>
coot::filter_residues_by_solvent_contact(CResidue *res_ref,
					 CMMDBManager *mol,
					 const std::vector<CResidue *> &residues,
					 const double &water_dist_max) {
   std::vector<CResidue *> v;
   PPCAtom lig_residue_atoms = 0;
   int n_lig_residue_atoms;
   res_ref->GetAtomTable(lig_residue_atoms, n_lig_residue_atoms);
   for (unsigned int i=0; i<residues.size(); i++) { 
      std::string res_name = residues[i]->GetResName();
      if (res_name != "HOH") { 
	 v.push_back(residues[i]);
      } else {
	 PPCAtom residue_atoms = 0;
	 int n_residue_atoms;
	 residues[i]->GetAtomTable(residue_atoms, n_residue_atoms);
	 bool i_added = 0;
	 for (unsigned int jat=0; jat<n_lig_residue_atoms; jat++) {
	    clipper::Coord_orth lig_pt(lig_residue_atoms[jat]->x,
				       lig_residue_atoms[jat]->y,
				       lig_residue_atoms[jat]->z);
	    for (unsigned int iat=0; iat<n_residue_atoms; iat++) {
	       clipper::Coord_orth pt(residue_atoms[iat]->x,
				      residue_atoms[iat]->y,
				      residue_atoms[iat]->z);
	       if ((lig_pt-pt).lengthsq() < (water_dist_max*water_dist_max)) {
		  i_added = 1;
		  v.push_back(residues[i]);
		  break;
	       }
	    }
	    if (i_added)
	       break;
	 }
      } 
   }
   return v;
}


std::pair<bool,float>
coot::closest_approach(CMMDBManager *mol,
		       CResidue *r1, CResidue *r2) {

    int n_res_1_atoms; 
    int n_res_2_atoms;
    PPCAtom residue_atoms_1 = 0, residue_atoms_2 = 0;
    double dist_best = 9999999.9;
    bool good_d = 0;
       
    r1->GetAtomTable( residue_atoms_1, n_res_1_atoms);
    r2->GetAtomTable( residue_atoms_2, n_res_2_atoms);
    for (int i=0; i<n_res_1_atoms; i++) { 
       clipper::Coord_orth pt1(residue_atoms_1[i]->x,
			       residue_atoms_1[i]->y,
			       residue_atoms_1[i]->z);
       for (int j=0; j<n_res_2_atoms; j++) { 
	  clipper::Coord_orth pt2(residue_atoms_2[j]->x,
				  residue_atoms_2[j]->y,
				  residue_atoms_2[j]->z);
	  double d = clipper::Coord_orth::length(pt1, pt2);

	  if (d < dist_best) {
	     dist_best = d;
	     good_d = 1;
	  }
       }
    }

    return std::pair<bool, float> (good_d, dist_best);
    
	  
   // Older, faster method.
   // It doesn't work though
   // 
   // I don't know why this fail, but it always gives 0 ncontacts.
   
//    std::pair<bool,realtype> r(0, 0.0);
   
//    int n_res_1_atoms; 
//    int n_res_2_atoms;
//    PPCAtom residue_atoms_1 = 0, residue_atoms_2 = 0;
//    int ncontacts = 0;
//    realtype dist_closest = 9999999999.9;
//    PSContact contact = NULL;
//    long i_contact_group = 0;
//    mat44 my_matt;
//    for (int i=0; i<4; i++) 
//       for (int j=0; j<4; j++) 
// 	 my_matt[i][j] = 0.0;      
//    for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
   
//    r1->GetAtomTable( residue_atoms_1, n_res_1_atoms);
//    r2->GetAtomTable( residue_atoms_2, n_res_2_atoms);
//    mol->SeekContacts(residue_atoms_1, n_res_1_atoms,
// 		     residue_atoms_2, n_res_2_atoms,
// 		     0.01, 10.0, 1, // not in the same residue
// 		     contact, ncontacts,
// 		     0, &my_matt, i_contact_group);


//    // debug mol
//    if (0) { 
//       int imod = 1;
//       CModel *model_p = mol->GetModel(imod);
//       CChain *chain_p;
//       // run over chains of the existing mol
//       int nchains = model_p->GetNumberOfChains();
//       for (int ichain=0; ichain<nchains; ichain++) {
// 	 chain_p = model_p->GetChain(ichain);
// 	 int nres = chain_p->GetNumberOfResidues();
// 	 PCResidue residue_p;
// 	 CAtom *at;
// 	 for (int ires=0; ires<nres; ires++) { 
// 	    residue_p = chain_p->GetResidue(ires);
// 	    int n_atoms = residue_p->GetNumberOfAtoms();
	 
// 	    for (int iat=0; iat<n_atoms; iat++) {
// 	       at = residue_p->GetAtom(iat);
// 	       std::cout << "mol atom "
// 			 << at ->name << " " 
// 			 << at ->GetSeqNum() << " " 
// 			 << at ->GetInsCode() << " " 
// 			 << at ->GetChainID() << " "
// 			 << at ->z << " "
// 			 << at ->y << " "
// 			 << at ->z << std::endl;
// 	    }
// 	 }
//       }
//    }

   

//    for (int i=0; i<n_res_1_atoms; i++)
//       std::cout << "closest_approach res_1 " << i << " "
// 		<< residue_atoms_1[i]->name << " ("
// 		<< residue_atoms_1[i]->x << " "
// 		<< residue_atoms_1[i]->y << " "
// 		<< residue_atoms_1[i]->z << ") "
// 		<< std::endl;
//    for (int i=0; i<n_res_2_atoms; i++)
//       std::cout << "closest_approach res_2 " << i << " " 
// 		<< residue_atoms_2[i]->name << " ("
// 		<< residue_atoms_2[i]->x << " "
// 		<< residue_atoms_2[i]->y << " "
// 		<< residue_atoms_2[i]->z << ") "
// 		<< std::endl;


//    for (int i=0; i<ncontacts; i++) {
      
//       clipper::Coord_orth atom_1(residue_atoms_1[ contact[i].id1 ]->x,
// 				 residue_atoms_1[ contact[i].id1 ]->y,
// 				 residue_atoms_1[ contact[i].id1 ]->z);
//       clipper::Coord_orth atom_2(residue_atoms_2[ contact[i].id2 ]->x,
// 				 residue_atoms_2[ contact[i].id2 ]->y,
// 				 residue_atoms_2[ contact[i].id2 ]->z);
//       float d = clipper::Coord_orth::length(atom_1, atom_2);
//       if (d < dist_closest) {
// 	 dist_closest = d;
// 	 r = std::pair<bool,realtype>(1, d);
//       } 
//    }
//    std::cout << "DEbug:: closest_approach() d " << r.first << " " << r.second
// 	     << " using " << ncontacts << " contacts with n_res_atoms "
// 	     << n_res_1_atoms << " " << n_res_2_atoms << std::endl;
//    return r;
}


// Return dist in Angstroms, can throw an exception if any of the
// atoms is null.
// 
double
coot::distance(CAtom *at_1, CAtom *at_2) {

   double d = -1;
   if (at_1 && at_2) {
      clipper::Coord_orth pt_1(at_1->x, at_1->y, at_1->z);
      clipper::Coord_orth pt_2(at_2->x, at_2->y, at_2->z);
      d = clipper::Coord_orth::length(pt_1, pt_2);
   }

   return d;
}

// Return angle in degrees, can throw an exception if any of the
// atoms is null.
// 
double
coot::angle(CAtom *at_1, CAtom *at_2, CAtom *at_3) {

   double ang = -1;

   if (at_1 && at_2 && at_3) {

      clipper::Coord_orth pt_1(at_1->x, at_1->y, at_1->z);
      clipper::Coord_orth pt_2(at_2->x, at_2->y, at_2->z);
      clipper::Coord_orth pt_3(at_3->x, at_3->y, at_3->z);

      ang = clipper::Util::rad2d(clipper::Coord_orth::angle(pt_1, pt_2, pt_3));

   } 
   return ang; 
} 



clipper::RTop_orth
coot::util::matrix_convert(mat44 mat) {
   
   clipper::Mat33<double> clipper_mat(mat[0][0], mat[0][1], mat[0][2],
				      mat[1][0], mat[1][1], mat[1][2],
				      mat[2][0], mat[2][1], mat[2][2]);
   clipper::Coord_orth cco(mat[0][3], mat[1][3], mat[2][3]);
   return clipper::RTop_orth(clipper_mat, cco);
}

// -------------------------------------------------------------
//                       quaternions
// -------------------------------------------------------------

clipper::Mat33<double>
coot::util::quaternion::matrix() const {

   clipper::Mat33<double> mat;

   mat(0,0) = 1.0 - 2.0 * (q1 * q1 + q2 * q2);
   mat(0,1) = 2.0 * (q0 * q1 - q2 * q3);
   mat(0,2) = 2.0 * (q2 * q0 + q1 * q3);
   
   mat(1,0) = 2.0 * (q0 * q1 + q2 * q3);
   mat(1,1)= 1.0 - 2.0 * (q2 * q2 + q0 * q0);
   mat(1,2) = 2.0 * (q1 * q2 - q0 * q3);
   
   mat(2,0) = 2.0 * (q2 * q0 - q1 * q3);
   mat(2,1) = 2.0 * (q1 * q2 + q0 * q3);
   mat(2,2) = 1.0 - 2.0 * (q1 * q1 + q0 * q0);
   
   return mat;
}

coot::util::quaternion::quaternion(const clipper::Mat33<double> &m) {

   float pw = 1 + m(0,0) + m(1,1) + m(2,2);
   float px = 1 + m(0,0) - m(1,1) - m(2,2);
   float py = 1 - m(0,0) + m(1,1) - m(2,2); 
   float pz = 1 - m(0,0) - m(1,1) + m(2,2);

   float pr1 = sqrt( (pw>0) ? pw : 0) / 2.0;
   float pr2 = sqrt( (px>0) ? px : 0) / 2.0;
   float pr3 = sqrt( (py>0) ? py : 0) / 2.0;
   float pr4 = sqrt( (pz>0) ? pz : 0) / 2.0;

   q0 = convert_sign(pr2, m(2,1) - m(1,2));
   q1 = convert_sign(pr3, m(0,2) - m(2,0));
   q2 = convert_sign(pr4, m(1,0) - m(0,1));
   q3 = pr1;
   
}

// Return x with the sign of y.
float 
coot::util::quaternion::convert_sign(const float &x, const float &y) const {

   if ((x > 0) && (y > 0)) return  x;
   if ((x < 0) && (y > 0)) return -x; 
   if ((x > 0) && (y < 0)) return -x; 
   return  x; 
}

// static 
void
coot::util::quaternion::test_quaternion() {

   // currently quaternions are tested in testcootutils

} 

std::ostream&  coot::util::operator<<(std::ostream& s, const coot::util::quaternion &q) {

   s << "(" << q.q0 << ", " << q.q1 << ", " << q.q2 << ", " << q.q3 << ")";
   return s;
}

std::ofstream& coot::util::operator<<(std::ofstream &s, const coot::util::quaternion &q) {

   // s << "(" << q.q0 << ", " << q.q1 << ", " << q.q2 << ", " << q.q3 << ")";
   return s;
} 


// -------------------------------------------------------------


// Urgh.  Should use a template...
bool
coot::is_member_p(const std::vector<CResidue *> &v, CResidue *a) {
   
   bool ir = 0;
   unsigned int vsize = v.size();

   for (unsigned int i=0; i<vsize; i++) { 
      if (v[i] == a) { 
	 ir = 1;
	 break;
      } 
   }
   return ir;
} 

bool
coot::is_main_chain_p(CAtom *at) { 

   std::string mol_atom_name(at->name);
   if (mol_atom_name == " N  " ||
       mol_atom_name == " C  " ||
       mol_atom_name == " H  " ||
       mol_atom_name == " CA " ||
       mol_atom_name == " HA " || // CA hydrogen
       mol_atom_name == " O  ") {
      return 1;
   } else {
      return 0;
   } 
}

bool
coot::is_main_chain_or_cb_p(CAtom *at) { 

   std::string mol_atom_name(at->name);
   if (mol_atom_name == " N  " ||
       mol_atom_name == " C  " ||
       mol_atom_name == " H  " ||
       mol_atom_name == " CA " ||
       mol_atom_name == " CB " ||
       mol_atom_name == " HA " || // CA hydrogen
       mol_atom_name == " O  ") {
      return 1;
   } else {
      return 0;
   } 
}

// return 0 or 1
bool coot::is_hydrogen_p(CAtom *at) {

   std::string mol_atom_ele(at->element);
   if (mol_atom_ele == " H" ||
       mol_atom_ele == " D") {
      return 1;
   } else {
      return 0;
   }
}

bool
coot::residues_in_order_p(CChain *chain_p) {

   bool ordered_flag = 1;

   if (chain_p) {

      int n_residues = chain_p->GetNumberOfResidues();
      int current_resno = -9999999;
      for (int ires=0; ires<n_residues; ires++) {
	 CResidue *res_p = chain_p->GetResidue(ires);
	 int seqnum = res_p->GetSeqNum();
	 if (seqnum < current_resno) {
	    ordered_flag = 0;
	    break;
	 } else {
	    current_resno = seqnum;
	 }
      }
   }
   return ordered_flag;
}

// Throw an exception if there is no consistent seg id for the
// atoms in the given residue.
std::string
coot::residue_atoms_segid(CResidue *residue_p) {

   int n_residue_atoms;
   PPCAtom residue_atoms;

   std::vector<std::string> seg_ids;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      CAtom *at = residue_atoms[iat];
      std::string seg_id = at->segID;
      if (seg_ids.size() == 0) {
	 seg_ids.push_back(seg_id);
      } else {
	 if (!coot::is_member_p(seg_ids, seg_id)) {
	    std::string mess = "No consistent segids for residue ";
	    mess += coot::util::int_to_string(residue_p->GetSeqNum());
	    throw std::runtime_error(mess);
	 }
      }
   }

   if (seg_ids.size() == 0) {
      std::string mess = "No segids for residue ";
      mess += coot::util::int_to_string(residue_p->GetSeqNum());
      throw std::runtime_error(mess);
   }
      
   return seg_ids[0]; 
}


// Use the above function the get the segid and insert it into all
// the atoms of receiver.
bool
coot::copy_segid(CResidue *provider, CResidue *receiver) {

   try {
      std::string s = coot::residue_atoms_segid(provider);
      int n_residue_atoms;
      PPCAtom residue_atoms;
      
      receiver->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
	 CAtom *at = residue_atoms[iat];
	 // BL says:: to set the SegID we just strcpy it, we shouldnt fiddle
	 // with GetIndex().
	 strcpy(at->segID, s.c_str());
	 // we want to set just the segid, but there is no function to
	 // do that (there is for others, e.g. element, atom-name etc.).
	 //	 at->SetAtomName(at->GetIndex(),
	 //		 at->serNum,
	 //		 at->GetAtomName(),
	 //		 at->altLoc,
	 //		 s.c_str(),
	 //		 at->GetElementName());

      }
   }

   catch (std::runtime_error mess) {
      // maybe do this.. not sure.
      std::cout << "   INFO:: " << mess.what() << std::endl;
   }

   return 1;

}

// Throw an exception if there is no consistent seg id for the
// atoms in the given chain.
std::string
coot::chain_atoms_segid(CChain *chain_p) {

 
   int n_residue_atoms;
   PPCAtom residue_atoms;

   std::vector<std::string> seg_ids;

   int n_residues = chain_p->GetNumberOfResidues();

   for (int ires=0; ires<n_residues; ires++) {
      CResidue *residue_p = chain_p->GetResidue(ires);
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
	 CAtom *at = residue_atoms[iat];
	 std::string seg_id = at->segID;
	 if (seg_ids.size() == 0) {
	    seg_ids.push_back(seg_id);
	 } else {
	    if (!coot::is_member_p(seg_ids, seg_id)) {
	       std::string mess = "No consistent segids for chain ";
	       mess += chain_p->GetChainID();
	       throw std::runtime_error(mess);
	    }
	 }
      }
   }
   

   if (seg_ids.size() == 0) {
      std::string mess = "No segids for chain ";
      mess += chain_p->GetChainID();
      throw std::runtime_error(mess);
   }
   return seg_ids[0]; 
} 



std::vector<std::string>
coot::util::residue_types_in_chain(CChain *chain_p) {

   std::vector<std::string> v;
   
   PCResidue residue_p;
   int nres = chain_p->GetNumberOfResidues();
   for (int ires=0; ires<nres; ires++) { 
      residue_p = chain_p->GetResidue(ires);
      if (residue_p) {
	 std::string n(residue_p->name);
	 if (! is_member_p(v, n))
	    v.push_back(n);
      }
   }
   return v;
}

std::vector<std::string>
coot::util::residue_types_in_residue_vec(const std::vector<CResidue *> &residues) {

   std::vector<std::string> v;
   for (unsigned int ires=0; ires<residues.size(); ires++) { 
      if (residues[ires]) {
	 std::string n(residues[ires]->name);
	 if (! is_member_p(v, n))
	    v.push_back(n);
      }
   }
   return v;
} 


// Return -1 on badness 
int
coot::util::max_number_of_residues_in_chain(CMMDBManager *mol) {

   int max_number_of_residues_in_chain = -1;
   if (mol) { 

      int n_models = mol->GetNumberOfModels();
      
      for (int imod=1; imod<=n_models; imod++) { 
      
	 CModel *model_p = mol->GetModel(imod);
   
	 CChain *chain;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 if (nchains <= 0) { 
	    std::cout << "bad nchains in max_number_of_residues_in_chain "
		      << nchains << std::endl;
	 } else { 
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain = model_p->GetChain(ichain);
	       if (chain == NULL) {  
		  // This should not be necessary. It seem to be a
		  // result of mmdb corruption elsewhere - possibly
		  // DeleteChain in update_molecule_to().
		  std::cout << "NULL chain in residues_in_molecule: "
			    << std::endl;
	       } else { 
		  int nres = chain->GetNumberOfResidues();
		  if (nres > max_number_of_residues_in_chain) {
		     max_number_of_residues_in_chain = nres;
		  }
	       }
	    }
	 }
      }
   }
   return max_number_of_residues_in_chain;
}

// Return -1 on badness.
// 
// So that we can calculate the lenght of the graph x axis - there may
// be gaps, which is why max_number_of_residues_in_chain would fail.
// 
int
coot::util::max_min_max_residue_range(CMMDBManager *mol) {

   int max_min_max = -1;
   if (mol) { 

      int n_models = mol->GetNumberOfModels();
      
      for (int imod=1; imod<=n_models; imod++) { 
      
	 CModel *model_p = mol->GetModel(imod);
   
	 CChain *chain;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 if (nchains <= 0) { 
	    std::cout << "bad nchains in max_min_max_residue_range "
		      << nchains << std::endl;
	 } else { 
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain = model_p->GetChain(ichain);
	       if (chain == NULL) {  
		  // This should not be necessary. It seem to be a
		  // result of mmdb corruption elsewhere - possibly
		  // DeleteChain in update_molecule_to().
		  std::cout << "NULL chain in residues_in_molecule: "
			    << std::endl;
	       } else { 
		  int nres = chain->GetNumberOfResidues();
		  if (nres > 0) { 
		     int min_resno = 99999;
		     int max_resno = -99999;
		     int this_resno;
		     for (int i=0; i<nres; i++) {
			this_resno = chain->GetResidue(i)->GetSeqNum();
			if (this_resno > max_resno)
			   max_resno = this_resno;
			if (this_resno < min_resno)
			   min_resno = this_resno;
		     }
		     int range = max_resno - min_resno + 1; // residues 1, 2 -> 2
		     if (range > max_min_max)
			max_min_max = range;
		  }
	       }
	    }
	 }
      }
   }
   return max_min_max;
}

// Return a vector of residue that are in this fragment.
// Fragments are marked by consecutively numbered residues.  A
// gap in the sequence numbers marks the end/beginning of a
// fragment.
std::vector<PCResidue>
coot::util::get_residues_in_fragment(CChain *chain_p,
				     coot::residue_spec_t clicked_residue) {

   std::vector<PCResidue> r;
   int nSelResidues;
   PCResidue *SelResidues;

   chain_p->GetResidueTable(SelResidues, nSelResidues);
   for (int i=0; i<nSelResidues; i++)
      r.push_back(SelResidues[i]);

   return r; 
}

std::pair<bool, clipper::Coord_orth>
coot::util::get_residue_centre(CResidue *residue_p) {

   bool status = 0;
   clipper::Coord_orth centre(0,0,0);
   
   PPCAtom residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   if (n_residue_atoms>0) {
      status = 1;
      for (unsigned int i=0; i<n_residue_atoms; i++) {
	 clipper::Coord_orth pt(residue_atoms[i]->x,
				residue_atoms[i]->y,
				residue_atoms[i]->z);
	 centre += pt;
      }
      double scale = 1.0/double(n_residue_atoms);
      centre = scale * centre;
   }
   return std::pair<bool, clipper::Coord_orth> (status, centre);
 
}




// Return -1 on badness (actually, number of chains in the last model)
int
coot::util::number_of_chains(CMMDBManager *mol) {

   int nchains = -1;

   if (mol) { 

      int n_models = mol->GetNumberOfModels();
      
      for (int imod=1; imod<=n_models; imod++) { 
      
	 CModel *model_p = mol->GetModel(imod);
   
	 // run over chains of the existing mol
	 nchains = model_p->GetNumberOfChains();
	 if (nchains <= 0) { 
	    std::cout << "bad nchains in number_of_chains "
		      << nchains << std::endl;
	 }
      }
   }
   return nchains;
}

// Return -1 on badness:
int
coot::get_selection_handle(CMMDBManager *mol, const coot::atom_spec_t &at) {

   int SelHnd = -1;
   if (mol) { 
      SelHnd = mol->NewSelection();
      char *chain   = (char *) at.chain.c_str();
      char *inscode = (char *) at.insertion_code.c_str();
      char *atname  = (char *) at.atom_name.c_str(); // atom name
      char *altconf = (char *) at.alt_conf.c_str();
      mol->SelectAtoms (SelHnd, 0, chain,
			at.resno, // starting resno, an int
			inscode, // any insertion code
			at.resno, // ending resno
			inscode, // ending insertion code
			"*", // any residue name
			atname,
			"*", // elements
			altconf  // alt loc.
			);
   }
   return SelHnd;
}

std::ostream& coot::operator<< (std::ostream& s, const coot::atom_spec_t &spec) {

   s << "[spec: ";
   s << "\"";
   s << spec.chain;
   s << "\" ";
   s << spec.resno;
   s << " ";
   s << "\"";
   s << spec.insertion_code;
   s << "\"";
   s << " ";
   s << "\"";
   s  << spec.atom_name;
   s << "\"";
   s << " ";
   s << "\"";
   s << spec.alt_conf;
   s << "\"]";

   return s;

}

std::ostream& coot::operator<< (std::ostream& s, const coot::residue_spec_t &spec) {

   if (!spec.unset_p()) { 

      s << "[spec: ";
      s << "\"";
      s << spec.chain;
      s << "\" ";
      s << spec.resno;
      s << " ";
      s << "\"";
      s << spec.insertion_code;
      s << "\"]";
   } else {
      s << "{residue-spec-not-set}";
   } 
   return s;

}


// deleted by calling process
std::pair<CMMDBManager *, std::vector<coot::residue_spec_t> >
coot::util::get_fragment_from_atom_spec(const coot::atom_spec_t &atom_spec,
					CMMDBManager *mol_in) {

   CMMDBManager *mol = 0;
   std::vector<coot::residue_spec_t> v;
   
   // Plan:
   //
   // We only want one model, so take the first.
   // 
   // First we find the residue in mol_in.
   //
   // Then find the top residue that is attached by direct attachment
   // (continuously increasing residue number).
   //
   // Then find the bottom residue that is attached by direct attachment
   // (continuously decreasing residue number).
   //
   // Then use and return create_mmdbmanager_from_res_selection()


   CResidue *residue_bot = 0;
   CResidue *residue_top = 0;
   CAtom *search_atom = 0;

   int resno_top = -999;
   int resno_bot = -999;

   int imod = 1;

   CModel *model_p = mol_in->GetModel(imod);
   CChain *chain_p;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      std::string mol_chain_id = chain_p->GetChainID();
      if (mol_chain_id == atom_spec.chain) {
	 int nres = chain_p->GetNumberOfResidues();
	 PCResidue residue_p;
	 CAtom *at;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    if (residue_p->GetSeqNum() == atom_spec.resno) {
	       int n_atoms = residue_p->GetNumberOfAtoms();
	    
	       for (int iat=0; iat<n_atoms; iat++) {
		  at = residue_p->GetAtom(iat);
		  std::string mol_atom_name = at->name;
		  if (mol_atom_name == atom_spec.atom_name) {
		     std::string alt_conf = at->altLoc;
		     if (alt_conf == atom_spec.alt_conf) {
			search_atom = at;
		     }
		  }
		  if (search_atom)
		     break;
	       }
	    }
	    if (search_atom)
	       break;
	 }
      }
      if (search_atom)
	 break;
   }

   bool found_search = 0;
   if (search_atom) {
      // now try to set res_bot and res_top and their resnos.
      CChain *chain_p = search_atom->GetChain();
      int nres = chain_p->GetNumberOfResidues();
      PCResidue residue_p;
      for (int ires=0; ires<nres; ires++) { 
	 residue_p = chain_p->GetResidue(ires);
	 if (residue_p == search_atom->GetResidue()) {

	    residue_top = residue_p;
	    residue_bot = residue_p;
	    int resno_this = residue_p->GetSeqNum();
	    resno_bot = resno_this;
	    resno_top = resno_this;

	    // search forwards on this chain
	    for (int ires_search=ires+1; ires_search<nres; ires_search++) {
	       int ioff = ires_search - ires; 
	       CResidue *residue_search_p = chain_p->GetResidue(ires_search);
	       if (residue_search_p->GetSeqNum() == resno_this + ioff) {
		  residue_top = residue_search_p;
		  resno_top = residue_search_p->GetSeqNum();
	       } else {
		  break;
	       }
	    }

	    // search backwards on this chain
	    for (int ires_search=ires-1; ires_search>=0; ires_search--) {
	       int ioff = ires_search - ires; 
	       CResidue *residue_search_p = chain_p->GetResidue(ires_search);
	       if (residue_search_p->GetSeqNum() == (resno_this + ioff)) {
		  residue_bot = residue_search_p;
		  resno_bot = residue_search_p->GetSeqNum();
	       } else {
		  break;
	       }
	       
	    }
	    
	    found_search = 1;
	 }
	 
	 if (found_search)
	    break;
      }
   }

   if (residue_bot && residue_top) {

      PCResidue *SelResidues = 0; 
      int nSelResidues;
      int selHnd = mol_in->NewSelection();

      mol_in->Select(selHnd, STYPE_RESIDUE, 1,
		     (char *) atom_spec.chain.c_str(),
		     resno_bot, "",
		     resno_top, "",
		     "*",  // residue name
		     "*",  // Residue must contain this atom name?
		     "*",  // Residue must contain this Element?
		     "*",  // altLocs
		     SKEY_NEW // selection key
		     );
      mol_in->GetSelIndex(selHnd, SelResidues, nSelResidues);

      std::pair<CMMDBManager *, int> mol_info = 
	 create_mmdbmanager_from_res_selection(mol_in, 
					       SelResidues, 
					       nSelResidues, 0, 0, "", atom_spec.chain, 0);
      if (mol_info.second) { 
	 mol = mol_info.first;
	 for (int ires=0; ires<nSelResidues; ires++) {
	    v.push_back(coot::residue_spec_t(SelResidues[ires]->GetChainID(),
					     SelResidues[ires]->GetSeqNum(),
					     SelResidues[ires]->GetInsCode()));
	 }
      }

      mol_in->DeleteSelection(selHnd);		     

   } else {
      if (! residue_top) 
	 std::cout << "ERROR:: missing top residue in fragment" << std::endl;
      if (! residue_bot) 
	 std::cout << "ERROR:: missing bot residue in fragment" << std::endl;
   }

   return std::pair<CMMDBManager *, std::vector<coot::residue_spec_t> > (mol, v);
} 

bool
coot::atom_spec_t::matches_spec(CAtom *atom) const {

   if (atom_name == std::string(atom->name)) {

      if (alt_conf == std::string(atom->altLoc)) {

	 CResidue *residue_p = atom->residue;
	 
	 if (residue_p) { 
	    
	    if (resno == atom->GetSeqNum()) {
	       
	       if (insertion_code == std::string(atom->GetInsCode())) { 
		  
		  CChain *chain_p= atom->GetChain();
		  if (chain_p) {
		     if (chain == chain_p->GetChainID()) {
			// std::cout << atom_name << "a complete match " << std::endl;
			return 1;
		     } else {
			// std::cout << atom_name << "a chain mismatch " << std::endl;
			return 0;
		     }
		  } else {
		     // std::cout << atom_name << "a no chain match " << std::endl;
		     // no chain
		     return 1;
		  }
	       } else {
		  // std::cout << atom_name << "an inscode mismatch " << std::endl;
		  return 0;
	       }
	    } else {
	       // std::cout << atom_name << "a resno mismatch " << std::endl;
	       return 0;
	    }
	    
	 } else {
	    // no residue
	    // std::cout << atom_name << "a no chain match " << std::endl;
	    return 1;
	 }
      } else {
	 // std::cout << atom_name << "an altloc mismatch " << std::endl;
	 return 0;
      } 
   } else {
      // std::cout << atom_name << "an atom name mismatch :" << atom->name << ":" << std::endl;
      return 0;
   }
   std::cout << atom_name << " should not happen (matches_spec()) " << atom->name << ":" << std::endl;
   return 0;
}

std::vector<CAtom * >
coot::torsion::matching_atoms(CResidue *residue) {

   std::vector<CAtom *> v;
   
   CAtom *catom_1 = 0;  
   CAtom *catom_2 = 0; 
   CAtom *catom_3 = 0; 
   CAtom *catom_4 = 0;

   PPCAtom residue_atoms;
   int nResidueAtoms;
   residue->GetAtomTable(residue_atoms, nResidueAtoms);
   for (int iat=0; iat<nResidueAtoms; iat++) {
      if (atom_1.second.matches_spec(residue_atoms[iat]))
	 catom_1 = residue_atoms[iat];
      if (atom_2.second.matches_spec(residue_atoms[iat]))
	 catom_2 = residue_atoms[iat];
      if (atom_3.second.matches_spec(residue_atoms[iat]))
	 catom_3 = residue_atoms[iat];
      if (atom_4.second.matches_spec(residue_atoms[iat]))
	 catom_4 = residue_atoms[iat];
   }

   if (! (catom_1 && catom_2 && catom_3 && catom_4)) {
      if (!catom_1)
	 std::cout << " atom_1 is null for " << atom_1.second.atom_name << std::endl;
      if (!catom_2)
	 std::cout << " atom_2 is null for " << atom_2.second.atom_name << std::endl;
      if (!catom_3)
	 std::cout << " atom_3 is null for " << atom_3.second.atom_name << std::endl;
      if (!catom_4)
	 std::cout << " atom_4 is null for " << atom_4.second.atom_name << std::endl;
   } else { 
      v.push_back(catom_1);
      v.push_back(catom_2);
      v.push_back(catom_3);
      v.push_back(catom_4);
   }

   return v;
}


std::vector<std::string>
coot::util::chains_in_molecule(CMMDBManager *mol) { 

   std::vector<std::string> v;

   if (mol) { 

      int n_models = mol->GetNumberOfModels();
      
      for (int imod=1; imod<=n_models; imod++) { 
      
	 CModel *model_p = mol->GetModel(imod);
   
	 CChain *chain;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 if (nchains <= 0) { 
	    std::cout << "bad nchains in trim molecule " << nchains
		      << std::endl;
	 } else { 
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain = model_p->GetChain(ichain);
	       if (chain == NULL) {  
		  // This should not be necessary. It seem to be a
		  // result of mmdb corruption elsewhere - possibly
		  // DeleteChain in update_molecule_to().
		  std::cout << "NULL chain in residues_in_molecule: "
			    << std::endl;
	       } else {
		  v.push_back(chain->GetChainID());
	       }
	    }
	 }
      }
   }
   return v;
}

std::pair<bool, int>
coot::util::min_resno_in_chain(CChain *chain_p) {

   bool found_residues = 0;
   int min_resno = 99999999;
   
   if (chain_p == NULL) {  
      // This should not be necessary. It seem to be a
      // result of mmdb corruption elsewhere - possibly
      // DeleteChain in update_molecule_to().
      std::cout << "NULL chain in residues_in_molecule: "
		<< std::endl;
   } else { 
      int nres = chain_p->GetNumberOfResidues();
      CResidue *residue_p;
      int resno;
      for (int ires=0; ires<nres; ires++) {
	 residue_p = chain_p->GetResidue(ires);
	 resno = residue_p->seqNum;
	 if (resno < min_resno) {
	    min_resno = resno;
	    found_residues = 1;
	 }
      }
   }
   return std::pair<bool, int>(found_residues, min_resno);
}


std::pair<bool, int>
coot::util::max_resno_in_chain(CChain *chain_p) {

   bool found_residues = 0;
   int max_resno = -31999;
   
   if (chain_p == NULL) {  
      // This should not be necessary. It seem to be a
      // result of mmdb corruption elsewhere - possibly
      // DeleteChain in update_molecule_to().
      std::cout << "NULL chain in residues_in_molecule: "
		<< std::endl;
   } else { 
      int nres = chain_p->GetNumberOfResidues();
      CResidue *residue_p;
      int resno;
      if (nres > 0) { 
	 for (int ires=0; ires<nres; ires++) {
	    residue_p = chain_p->GetResidue(ires);
	    resno = residue_p->seqNum;
	    if (resno > max_resno) {
	       max_resno = resno;
	       found_residues = 1;
	    }
	 }
      } else {
	 // there was a chain, but no residues in it. which is what
	 // happens in pointer_atom_molecule() when we create a
	 // Pointer Atom molecule.
	 //
	 // In that case, we want the added atom to have residue
	 // number 1.
	 max_resno = 0;
      } 
   }
//    std::cout << "DEBUG:: max_resno_in_chain returning " << found_residues
// 	     << " " << max_resno << std::endl;
   return std::pair<bool, int>(found_residues, max_resno);
}

std::pair<bool, int>
coot::util::max_resno_in_molecule(CMMDBManager *mol) {

   bool found_residues = 0;
   int current_high = -31999;
   
   int n_models = mol->GetNumberOfModels();
   for (int imod=1; imod<=n_models; imod++) { 
      
      CModel *model_p = mol->GetModel(imod);
      CChain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 std::pair<bool, int> p = coot::util::max_resno_in_chain(chain_p);
	 if (p.first) { 
	    if (p.second > current_high) { 
	       current_high = p.second;
	       found_residues = 1;
	    }
	 }
      }
   }
   return std::pair<bool, int> (found_residues, current_high);
} 


int
coot::util::number_of_residues_in_molecule(CMMDBManager *mol) {

   int nres = 0;
   if (mol) { 

      int n_models = mol->GetNumberOfModels();
      
      for (int imod=1; imod<=n_models; imod++) { 
      
	 CModel *model_p = mol->GetModel(imod);
   
	 CChain *chain_p;
	 CResidue *residue_p;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 if (nchains <= 0) { 
	    std::cout << "bad nchains in trim molecule " << nchains
		      << std::endl;
	 } else { 
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain_p = model_p->GetChain(ichain);
	       if (chain_p == NULL) {  
		  // This should not be necessary. It seem to be a
		  // result of mmdb corruption elsewhere - possibly
		  // DeleteChain in update_molecule_to().
		  std::cout << "NULL chain in number residues_in_molecule: "
			    << std::endl;
	       } else { 
		  int nres = chain_p->GetNumberOfResidues();
		  for (int ires=0; ires<nres; ires++) { 
		     residue_p = chain_p->GetResidue(ires);
		     if (residue_p != NULL) {
			nres++;
		     }
		  }
	       }
	    }
	 }
      }
   }
   return nres;
}


// Enough kludgey code.  Utilities.
//
// On failure, return "";
std::string
coot::util::single_letter_to_3_letter_code(char code) {

   if (code == 'G') return std::string("GLY");
   if (code == 'A') return std::string("ALA");
   if (code == 'V') return std::string("VAL");
   if (code == 'S') return std::string("SER");
   if (code == 'N') return std::string("ASN");
   if (code == 'P') return std::string("PRO");
   if (code == 'D') return std::string("ASP");
   if (code == 'C') return std::string("CYS");
   if (code == 'Q') return std::string("GLN");
   if (code == 'E') return std::string("GLU");
   if (code == 'H') return std::string("HIS");
   if (code == 'I') return std::string("ILE");
   if (code == 'L') return std::string("LEU");
   if (code == 'K') return std::string("LYS");
   if (code == 'M') return std::string("MET");
   if (code == 'F') return std::string("PHE");
   if (code == 'T') return std::string("THR");
   if (code == 'W') return std::string("TRP");
   if (code == 'Y') return std::string("TYR");
   if (code == 'R') return std::string("ARG");

   return std::string("");
}


// Match on graph
// 
// Return the orientation matrix moving res_moving to res_reference
// and a flag letting us know that the match worked OK.
// 
coot::graph_match_info_t
coot::graph_match(CResidue *res_moving,
		  CResidue *res_reference,
		  bool apply_rtop_flag) {

  clipper::Mat33<double> m_dum(1,0,0,0,1,0,0,0,1);
  clipper::Coord_orth pt_dum(0,0,0);
  clipper::RTop_orth rtop(m_dum, pt_dum);
   bool success = 0;
   std::vector<std::pair<std::pair<std::string, std::string>, std::pair<std::string, std::string> > > best_matching_atoms;

   CGraph graph1;
   CGraph graph2;

   // These are deleted at the end
   CResidue *cleaned_res_moving    = coot::util::copy_and_delete_hydrogens(res_moving);
   CResidue *cleaned_res_reference = coot::util::copy_and_delete_hydrogens(res_reference);

   // debug 
   if (0) {
      int n_residue_atoms_1;
      PPCAtom residue_atoms_1;
      cleaned_res_moving->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
      int n_residue_atoms_2;
      PPCAtom residue_atoms_2;
      cleaned_res_moving->GetAtomTable(residue_atoms_2, n_residue_atoms_2);
      // are these the same atoms?
      for (int i=0; i<4; i++) {
	 std::cout << "moving and ref atoms: " << residue_atoms_1[i] <<  " " << residue_atoms_2[i]
		   << std::endl;
      }
   }

   graph1.MakeGraph(cleaned_res_moving);
   graph2.MakeGraph(cleaned_res_reference);

   // graph1.MakeSymmetryRelief ( False );
   // graph2.MakeSymmetryRelief ( False );

   int build_status1 = graph1.Build(1);
   int build_status2 = graph2.Build(1);
   double best_match_sum = 1e20;
   int best_n_match = -99;

   if (build_status1 != 0) {
      std::cout << "ERROR:: build_status1: " << build_status1 << std::endl;
   } else { 
      if (build_status2 != 0) {
	 std::cout << "ERROR:: build_status2: " << build_status2 << std::endl;
      } else {

	 int n_atoms_ref = cleaned_res_reference->GetNumberOfAtoms();
	 int n_atoms_mov = cleaned_res_moving->GetNumberOfAtoms();

	 int minMatch = 4;
	 int n_ref_frac = int(0.75*float(n_atoms_ref));
	 int n_mov_frac = int(0.75*float(n_atoms_mov));

	 int min_n = (n_ref_frac < n_mov_frac) ? n_ref_frac : n_mov_frac;
	 if (min_n > minMatch)
	    minMatch = min_n;
	 
	 CGraphMatch match;

	 std::cout << "match.MatchGraphs must match at least " << minMatch << " atoms."
		   << std::endl;
	 match.MatchGraphs(&graph1, &graph2, minMatch, 1);
	 int n_match = match.GetNofMatches();
	 std::cout << "INFO:: match NumberofMatches (potentially similar graphs) "
		   << n_match << std::endl;
	 // match.PrintMatches();

	 int best_match = -1;
	 clipper::Mat33<double> m_dum(1,0,0,0,1,0,0,0,1);
	 clipper::Coord_orth pt_dum(0,0,0);
	 clipper::RTop_orth best_rtop(m_dum, pt_dum);
	 for (int imatch=0; imatch<n_match; imatch++) {
	    std::vector<std::pair<std::pair<std::string, std::string>, std::pair<std::string, std::string> > > matching_atoms; 
	    int n;
	    realtype p1, p2;
	    ivector FV1, FV2;
	    match.GetMatch(imatch, FV1, FV2, n, p1, p2); // n p1 p2 set
// 	    For understanding only.  
// 	    std::cout << "Match number: " << imatch << "  " << p1*100 << "% "
// 		      << p2*100 << "% "<< std::endl;
	    std::vector<clipper::Coord_orth> coords_1_local;
	    std::vector<clipper::Coord_orth> coords_2_local;
	    for (int ipair=1; ipair<=n; ipair++) {
	       PCVertex V1 = graph1.GetVertex ( FV1[ipair] );
	       PCVertex V2 = graph2.GetVertex ( FV2[ipair] );
	       if ((!V1) || (!V2))  {
		  std::cout << "Can't get vertices for match "
			    << ipair << std::endl;
	       } else  {
// 		  printf(" %4i.  [%4s] <-> [%4s]\n",
// 			 ipair, V1->GetName(), V2->GetName());
		  CAtom *at1 = cleaned_res_moving->atom[V1->GetUserID()];
		  CAtom *at2 = cleaned_res_reference->atom[V2->GetUserID()];
		  coords_1_local.push_back(clipper::Coord_orth(at1->x, at1->y, at1->z));
		  coords_2_local.push_back(clipper::Coord_orth(at2->x, at2->y, at2->z));
		  std::pair<std::string, std::string> atom_info_1(at1->name, at1->altLoc);
		  std::pair<std::string, std::string> atom_info_2(at2->name, at2->altLoc);
		  std::pair<std::pair<std::string, std::string>, std::pair<std::string, std::string> > atom_pair(atom_info_1, atom_info_2);
		  matching_atoms.push_back(atom_pair);
	       }
	    }
	    
	    double dist_sum = 0.0;
	    clipper::RTop_orth rtop_local(clipper::Mat33<double>(0,0,0,0,0,0,0,0,0),
					  clipper::Coord_orth(0,0,0)); // unset
	    if (apply_rtop_flag) {

// 	       for (unsigned int iat=0; iat<4; iat++) {
// 		  std::cout << "debug:: getting rtop " << coords_1_local[iat].format() << " vs "
// 			    << coords_2_local[iat].format() << std::endl;
// 	       }
	       
	       rtop_local = clipper::RTop_orth(coords_1_local, coords_2_local);
	       for (unsigned int i=0; i<coords_1_local.size(); i++) {
		  dist_sum += clipper::Coord_orth::length(coords_2_local[i],
							  coords_1_local[i].transform(rtop_local));
	       }
	    } else {
	       for (unsigned int i=0; i<coords_1_local.size(); i++) {
		  dist_sum += clipper::Coord_orth::length(coords_2_local[i], coords_1_local[i]);
	       }
	    } 
	    if (dist_sum < best_match_sum) {
	       // Debugging
	       std::cout << "DEBUG:: better dist_sum: " << dist_sum << std::endl;
	       best_rtop = rtop_local;
	       best_match_sum = dist_sum;
	       best_match = imatch;
	       best_matching_atoms = matching_atoms;
	       best_n_match = coords_1_local.size();
	    }
	 } // imatch loop

	 if (best_match != -1) {
	    rtop = best_rtop;
	    success = 1;
	 }
      }
   }
   delete cleaned_res_reference;
   delete cleaned_res_moving;
   coot::graph_match_info_t gmi;
   gmi.success = success;
   gmi.rtop = rtop;
   gmi.dist_score = best_match_sum;
   gmi.matching_atom_names = best_matching_atoms;
   gmi.n_match = best_n_match;
   return gmi;
}

// Change the names in res_moving_names to match those in
// res_reference as much as possible.  When there is a name collision
// (i.e. the name maped from the res_reference is already in the
// res_moving_names (and that is not to be replace by anything)),
// invent a new name for the atom.  Use the internal
// matching_atom_names.
void
coot::graph_match_info_t::match_names(CResidue *res_with_moving_names) {

   if (!success) {
      std::cout << "Can't do name remapping, graph match failed" << std::endl;
   } else { 

      // non-mapped and not the same, that is.
      std::vector<std::string> orig_moving_atom_names_non_mapped_non_same;

      // not in the set of graph matched pairs, not changing (presumably).
      std::vector<std::string> orig_moving_atom_names_non_mapped_same;
   
      std::vector<std::string> residue_atom_names;
   

      // fill orig_moving_atom_names_non_mapped_non_same
      PPCAtom residue_atoms = 0;
      int n_residue_atoms;
      res_with_moving_names->GetAtomTable(residue_atoms, n_residue_atoms);
      for (unsigned int iat=0; iat<n_residue_atoms; iat++) {
	 std::string atom_name(residue_atoms[iat]->name);
	 // add the name if is not already there.
	 if (std::find(residue_atom_names.begin(), residue_atom_names.end(), atom_name)
	     == residue_atom_names.end())
	    residue_atom_names.push_back(atom_name);
	 bool found_match = 0;
	 for (unsigned int j_pair=0; j_pair<matching_atom_names.size(); j_pair++) {
	    // first is working atom spec
	    if (matching_atom_names[j_pair].first.first == atom_name) {
	       found_match = 1;
	       break;
	    }
	 }

	 std::cout << ".... atom name: " << atom_name << ": found_match " << found_match << std::endl;
	 if (! found_match) {
	    // this atom name was not in the list of working atoms that were matched.

	    // now, was this atom name in the list of matched target atoms?
	    // If so, it is a collision-atom
	    bool found_match_2 = 0;
	    for (unsigned int j_pair=0; j_pair<matching_atom_names.size(); j_pair++) {
	       // second is reference atom spec
	       if (matching_atom_names[j_pair].second.first == atom_name) {
		  found_match_2 = 1;
		  break;
	       }
	    }
	    std::cout << "      found_match_2 " << found_match_2 << std::endl;
	    if (found_match_2)
	       orig_moving_atom_names_non_mapped_non_same.push_back(atom_name);
	    else
	       orig_moving_atom_names_non_mapped_same.push_back(atom_name);
	    
	 } 
      }

      if (1) {
	 std::cout << "Mapped atom names: " << matching_atom_names.size() << std::endl;
	 // according to header ref is first, work is second.
	 for (unsigned int i=0; i<matching_atom_names.size(); i++) { 
	    std::cout << "   " << i << " :" << matching_atom_names[i].first.first 
		      << ": -> :" << matching_atom_names[i].second.first << ":"
		      << std::endl;
	 }
	 std::cout << "Non-mapped non-same atom names: " << orig_moving_atom_names_non_mapped_non_same.size()
		   << std::endl;
	 for (unsigned int i=0; i<orig_moving_atom_names_non_mapped_non_same.size(); i++) { 
	    std::cout << "   " << i << " :" << orig_moving_atom_names_non_mapped_non_same[i] << ":" << std::endl;
	 }
	 std::cout << "Non mapped same atom names: " << orig_moving_atom_names_non_mapped_same.size()
		   << std::endl;
	 for (unsigned int i=0; i<orig_moving_atom_names_non_mapped_same.size(); i++) { 
	    std::cout << "   " << i << " :" << orig_moving_atom_names_non_mapped_same[i] << ":" << std::endl;
	 }
      }

      for (unsigned int iat=0; iat<n_residue_atoms; iat++) {

	 // check for a collision.  Is the reference name in existing
	 // atom names that are not due to be replaced?
	 std::string this_atom_name(residue_atoms[iat]->name);

	 bool replace_name = 0;
	 std::string new_atom_name = "";
      
      
	 if (std::find(orig_moving_atom_names_non_mapped_non_same.begin(),
		       orig_moving_atom_names_non_mapped_non_same.end(),
		       this_atom_name) !=
	     orig_moving_atom_names_non_mapped_non_same.end()) {
	    // OK, this atom name is in the list of atoms non-mapped needing a name change
	    std::string ele = residue_atoms[iat]->element;
	    new_atom_name = invent_new_name(this_atom_name, ele, residue_atom_names);
	    residue_atom_names.push_back(new_atom_name);
	    replace_name = 1;
	    // std::cout << ":" << this_atom_name << ": non-mapped..." << std::endl;
	 } else {

	    // OK, this moving residue atom is not non-mapped and needing a name change.

	    // Is it non-mapped and not needing a name change?
	    //
	    //   If it is, don't set replace name,
	    ///  If is not, then it should be in the (graph) matching atom names, so
	    //      check there if the atom names are the same and if they are, do not
	    //      set the replace_name flag

	    if (std::find(orig_moving_atom_names_non_mapped_same.begin(),
			  orig_moving_atom_names_non_mapped_same.end(),
			  this_atom_name) !=
		orig_moving_atom_names_non_mapped_same.end()) {

	       // no change to replace_name.

	    } else { 

	       // std::cout << ":" << this_atom_name << ": mapped" << std::endl;
	       replace_name = 1;
	       for (unsigned int j_pair=0; j_pair<matching_atom_names.size(); j_pair++) {
		  if (matching_atom_names[j_pair].first.first == this_atom_name) { 
		     // if the atom name is the same the reference atom name, then
		     // nothing to do.
		     if (matching_atom_names[j_pair].second.first ==
			 matching_atom_names[j_pair].first.first) {
			replace_name = 0;
			break;
		     } else {
			new_atom_name = matching_atom_names[j_pair].second.first;
		     } 
		  }
	       }
	    }
	 }

	 std::cout << "debug atom name :" << this_atom_name
		   << ": replace status: " << replace_name
		   << " new name :" << new_atom_name << ":" << std::endl;
	 if (replace_name) {
	    residue_atoms[iat]->SetAtomName(new_atom_name.c_str());
	 } 
      
      }
   }
} 

std::string
coot::graph_match_info_t::invent_new_name(const std::string &name_in,
					  const std::string &ele,
					  const std::vector<std::string> &residue_atom_name) const {
   std::string a("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
   bool found = 0;
   std::string new_name("XXXX");
   for (unsigned int i=0; i<a.size(); i++) { 
      for (unsigned int j=0; j<a.size(); j++) {
	 std::string test_atom_name = "";
	 if (ele.length() == 1) { 
	    test_atom_name = " ";
	    test_atom_name += ele;
	 } else {
	    test_atom_name = ele;
	 }
	 test_atom_name += a[i];
	 test_atom_name += a[j];
	 if (std::find(residue_atom_name.begin(), residue_atom_name.end(), test_atom_name)
	     == residue_atom_name.end()) {
	    found = 1;
	    new_name = test_atom_name;
	 }
	 if (found)
	    break;
      }
      if (found)
	 break;
   }
   return new_name;
} 



float
coot::util::median_temperature_factor(PPCAtom atom_selection,
				      int n_atoms,
				      float low_cutoff,
				      float high_cutoff,
				      short int apply_low_cutoff,
				      short int apply_high_cutoff) {

   float this_b;
   float median = 0;
   std::vector<float> b;
   for (int i=0; i<n_atoms; i++) {
      this_b = atom_selection[i]->tempFactor;
      if ((apply_low_cutoff && (this_b > low_cutoff)) ||
	  !apply_low_cutoff) {
	 if ((apply_high_cutoff && (this_b > high_cutoff)) ||
	     !apply_high_cutoff) {
	    b.push_back(this_b);
	 }
      }
   }

   if (b.size() > 0) { 
      std::sort(b.begin(), b.end());
      int mid_pos = b.size()/2;
      median = b[mid_pos];
   }
   return median;
} 

float
coot::util::average_temperature_factor(PPCAtom atom_selection,
				       int n_atoms,
				       float low_cutoff,
				       float high_cutoff,
				       short int apply_low_cutoff,
				       short int apply_high_cutoff) {

   float this_b = 0.0;
   float b_sum = 0.0;
   int n_sum = 0;

   for (int i=0; i<n_atoms; i++) {
      this_b = atom_selection[i]->tempFactor;
      if ((apply_low_cutoff && (this_b > low_cutoff)) ||
	  !apply_low_cutoff) {
	 if ((apply_high_cutoff && (this_b > high_cutoff)) ||
	     !apply_high_cutoff) {
	    b_sum += this_b;
	    n_sum++;
	 }
      }
   }

   float mean = 0.0;
   if (n_atoms > 0)
      mean = b_sum/float(n_atoms);
   return mean;
} 

float
coot::util::standard_deviation_temperature_factor(PPCAtom atom_selection,
						  int n_atoms,
						  float low_cutoff,
						  float high_cutoff,
						  short int apply_low_cutoff,
						  short int apply_high_cutoff) {

   double this_b = 0.0;
   double b_sum = 0.0;
   double b_sum_sqs = 0.0;
   int n_sum = 0;

   for (int i=0; i<n_atoms; i++) {
      this_b = atom_selection[i]->tempFactor;
      if ((apply_low_cutoff && (this_b > low_cutoff)) ||
	  !apply_low_cutoff) {
	 if ((apply_high_cutoff && (this_b > high_cutoff)) ||
	     !apply_high_cutoff) {
	    b_sum += this_b;
	    b_sum_sqs += this_b * this_b;
	    n_sum++;
	 }
      }
   }

   double mean = 0.0;
   double var  = 0.0;
   float sd = 0.0;
   if (n_atoms > 0) { 
      mean = b_sum/double(n_atoms);
      var = b_sum_sqs/double(n_atoms) - mean * mean;
      if (var < 0.0)
	 var = 0.0;
      sd = sqrt(var);
   } 
   return sd;
} 

// Return NULL on residue not found in this molecule.
// 
CResidue *
coot::util::get_residue(const std::string &chain_id,
			int reso, const std::string &insertion_code,
			CMMDBManager *mol) {

   CResidue *res = NULL;
   bool found_res = 0;

   if (mol) {
      CModel *model_p = mol->GetModel(1);
      CChain *chain_p;
      int n_chains = model_p->GetNumberOfChains(); 
      for (int i_chain=0; i_chain<n_chains; i_chain++) {
	 chain_p = model_p->GetChain(i_chain);
	 std::string mol_chain(chain_p->GetChainID());
	 if (mol_chain == chain_id) {
	    int nres = chain_p->GetNumberOfResidues();
	    CResidue *residue_p;
	    for (int ires=0; ires<nres; ires++) { // ires is a serial number
	       residue_p = chain_p->GetResidue(ires);
	       if (residue_p->GetSeqNum() == reso) {
		  std::string ins_code(residue_p->GetInsCode());
		  if (insertion_code == ins_code) {
		     res = residue_p;
		     found_res = 1;
		     break;
		  }
	       }
	    }
	 }
	 if (found_res) break;
      }
   }
   return res;
}

CResidue *
coot::util::get_residue(const residue_spec_t &rs, CMMDBManager *mol) {
   return get_residue(rs.chain, rs.resno, rs.insertion_code, mol);
}
  

// Return NULL on residue not found in this molecule.
// 
CResidue *
coot::util::get_following_residue(const residue_spec_t &rs, 
				  CMMDBManager *mol) {

   CResidue *res = NULL;
   if (mol) {
      CModel *model_p = mol->GetModel(1);
      CChain *chain_p;
      CChain *chain_this_res = NULL;
      bool found_this_res = 0;
      
      int n_chains = model_p->GetNumberOfChains();
      for (int i_chain=0; i_chain<n_chains; i_chain++) {
	 chain_p = model_p->GetChain(i_chain);
	 std::string mol_chain_id(chain_p->GetChainID());
	 if (mol_chain_id == rs.chain) {
	    int nres = chain_p->GetNumberOfResidues();
	    CResidue *residue_p;
	    for (int ires=0; ires<nres; ires++) {
	       residue_p = chain_p->GetResidue(ires);
	       if (found_this_res == 0) { 
		  if (rs.resno == residue_p->GetSeqNum()) {
		     std::string ins_code = residue_p->GetInsCode();
		     if (ins_code == rs.insertion_code) {
			found_this_res = 1;
			chain_this_res = chain_p;
		     }
		  }
	       } else {
		  // previous was our residue
		  if (chain_p == chain_this_res) {
		     // next residue in the same chain
		     res = residue_p;
		     break;
		  } 
	       }
	    }
	 }
	 if (res) break;
      }
   }
   return res;
}

CResidue *
coot::util::get_first_residue(CMMDBManager *mol) {

   CResidue *res = NULL;
   if (mol) {
      CModel *model_p = mol->GetModel(1);
      CChain *chain_p;
      
      int n_chains = model_p->GetNumberOfChains();
      for (int i_chain=0; i_chain<n_chains; i_chain++) {
	 chain_p = model_p->GetChain(i_chain);
	 int nres = chain_p->GetNumberOfResidues();
	 CResidue *residue_p;
	 for (int ires=0; ires<nres; ires++) {
	    residue_p = chain_p->GetResidue(ires);
	    if (residue_p) {
	       res = residue_p;
	       break;
	    }
	 }
	 if (res)
	    break;
      }
   }
   return res;
}


// return success status, 1 is good, 0 is fail.  Use clipper::Coord_orth constructor
// 
bool
coot::util::add_atom(CResidue *res,
		     const std::string &atom_name_1,
		     const std::string &atom_name_2,
		     const std::string &atom_name_3,
		     const std::string &alt_conf,
		     double length,
		     double angle, // degrees
		     double torsion, // degrees
		     const std::string &new_atom_name,
		     const std::string &new_atom_ele,
		     float new_atom_occ,
		     float new_atom_b_factor) {

   bool added_status = 0;
   int n_found = 0; 
   CAtom *a = NULL;
   CAtom *b = NULL;
   CAtom *c = NULL;

   int nResidueAtoms;
   PPCAtom residue_atoms;
   if (res) { 
      res->GetAtomTable(residue_atoms, nResidueAtoms);
      for (int i=0; i<nResidueAtoms; i++) {
	 std::string atom_name(residue_atoms[i]->name);
	 std::string atom_alt_conf(residue_atoms[i]->altLoc);
	 if (atom_alt_conf == alt_conf) { 
	    if (atom_name == atom_name_1) {
	       a = residue_atoms[i];
	    } 
	    if (atom_name == atom_name_2) {
	       b = residue_atoms[i];
	    } 
	    if (atom_name == atom_name_3) {
	       c = residue_atoms[i];
	    }
	 }
      }

      if (a && b && c) {
	 clipper::Coord_orth ac(a->x, a->y, a->z);
	 clipper::Coord_orth bc(b->x, b->y, b->z);
	 clipper::Coord_orth cc(c->x, c->y, c->z);
	 double ang  = clipper::Util::d2rad(angle);
	 double tors = clipper::Util::d2rad(torsion);
	 clipper::Coord_orth pos(ac, bc, cc, length, ang, tors);
	 CAtom *new_atom = new CAtom();
	 new_atom->SetCoordinates(pos.x(), pos.y(), pos.z(), new_atom_occ, new_atom_b_factor);
	 new_atom->SetAtomName(new_atom_name.c_str());
	 new_atom->SetElementName(new_atom_ele.c_str());
	 res->AddAtom(new_atom);
	 added_status = 1;
      
      } else {
	 std::cout << "Failed to find all reference atoms : "
		   << atom_name_1 << " " 
		   << atom_name_2 << " " 
		   << atom_name_3 << ". Found " << n_found << " out of 3" << std::endl;
      }
   }
   return added_status;
} 



std::vector<std::string>
coot::util::get_residue_alt_confs(CResidue *res) {

   std::vector<std::string> v;
   PPCAtom residue_atoms;
   int nResidueAtoms;
   res->GetAtomTable(residue_atoms, nResidueAtoms);
   bool ifound = 0;
   for (int iat=0; iat<nResidueAtoms; iat++) {
      ifound = 0;
      for(unsigned int i=0; i<v.size(); i++) { 
	 if (std::string(residue_atoms[iat]->altLoc) == v[i]) {
	    ifound = 1;
	    break;
	 }
      }
      if (! ifound) 
	 v.push_back(std::string(residue_atoms[iat]->altLoc));
   }
   return v;
} 

      




// The flanking residues (if any) are in the residue selection (SelResidues).
// The flags are not needed now we have made adjustments in the calling
// function.
// 
// create_mmdbmanager_from_res_selection must make adjustments
// 
// Note: there is now a molecule-class-info version of this - perhaps
// we should call it?  Next bug fix here: move over to the function call.
// 
// We need to pass orig_mol because of atom index transfer
std::pair<CMMDBManager *, int>
coot::util::create_mmdbmanager_from_res_selection(CMMDBManager *orig_mol,
						  PCResidue *SelResidues, 
						  int nSelResidues, 
						  int have_flanking_residue_at_start,
						  int have_flanking_residue_at_end, 
						  const std::string &altconf,
						  const std::string &chain_id_1,
						  short int residue_from_alt_conf_split_flag) {

   int start_offset = 0;
   int end_offset = 0;
   
//    if (have_flanking_residue_at_start)
//       start_offset = -1;
//    if (have_flanking_residue_at_end)
//       end_offset = +1; 

   CMMDBManager *residues_mol = new CMMDBManager;
   CModel *model = new CModel;
   CChain *chain = new CChain;
   short int whole_res_flag = 0; // not all alt confs, only this one ("A") and "".

   // For the active residue range (i.e. not the flanking residues) we only want
   // to refine the atoms that have the alt conf the same as the picked atom
   // (and that is altconf, passed here).
   // 
   // However, for *flanking residues* it's different.  Say we are refining a
   // non-split residue with alt conf "": Say that residue has a flanking
   // residue that is completely split, into A and B.  In that case we want
   // either "" or "A" for the flanking atoms.
   // 
   // And say we want to refine the A alt conf of a completely split residue
   // that has a flanking neighbour that is completely unsplit (""), we want
   // atoms that are either "A" or "".
   // 
   // So let's try setting whole_res_flag to 1 for flanking residues.

   CResidue *r;
   int atom_index_handle = residues_mol->RegisterUDInteger(UDR_ATOM, "mol's atom index");
   int afix_handle_orig = orig_mol->GetUDDHandle(UDR_ATOM, "shelx afix");
   int afix_handle_new_mol = -1;
   if (afix_handle_orig >= 0) { 
      afix_handle_new_mol = residues_mol->RegisterUDInteger(UDR_ATOM, "shelx afix");
      // int udd_afix_handle_test = residues_mol->GetUDDHandle(UDR_ATOM, "shelx afix");
//       std::cout << "DEBUG:: in create_mmdbmanager_from_res_selection, "
// 		<< "test : " << afix_handle_new_mol << " vs " << udd_afix_handle_test << std::endl;
   }

   for (int ires=start_offset; ires<(nSelResidues + end_offset); ires++) { 

      if ( (ires == 0) || (ires == nSelResidues -1) ) { 
	 if (! residue_from_alt_conf_split_flag)
	    whole_res_flag = 1;
      } else { 
	 whole_res_flag = 0;
      }

      r = coot::util::deep_copy_this_residue_with_atom_index_and_afix_transfer(orig_mol, SelResidues[ires], altconf, whole_res_flag, atom_index_handle, afix_handle_new_mol);
      
      chain->AddResidue(r);
      r->seqNum = SelResidues[ires]->GetSeqNum();
      r->SetResName(SelResidues[ires]->GetResName());
   }
   chain->SetChainID(chain_id_1.c_str());
   model->AddChain(chain);
   residues_mol->AddModel(model);
   residues_mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
   residues_mol->FinishStructEdit();

   if (afix_handle_orig >= 0) { 
      afix_handle_new_mol = residues_mol->GetUDDHandle(UDR_ATOM, "shelx afix");
      int imod = 1;
      
      CModel *model_p = residues_mol->GetModel(imod);
      CChain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 PCResidue residue_p;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    // 	    int n_atoms = residue_p->GetNumberOfAtoms();
	    
// 	    for (int iat=0; iat<n_atoms; iat++) {
// 	       at = residue_p->GetAtom(iat);
// 	       int check_afix_number;
// 	       if (at->GetUDData(afix_handle_new_mol, check_afix_number) == UDDATA_Ok) {
// 		  std::cout << "atom " << at << " has afix handle " << check_afix_number
// 			    << std::endl;
// 	       } else {
// 		  std::cout << "Failed to get afix number right after set! "
// 			    << at << std::endl;
//	       }
//          }
	 }
      }
   }

//    int udd_afix_handle_inter = residues_mol->GetUDDHandle(UDR_ATOM, "shelx afix");
//    std::cout << "DEBUG:: about to return from create_mmdbmanager_from_res_selection, "
// 	     << "udd_afix_handle_inter : " << udd_afix_handle_inter << std::endl;

   return std::pair<CMMDBManager *, int>(residues_mol, atom_index_handle);
}

// More complex than above because res_vec is not sorted on chain or residue number.
//
// The residues are added in order, the chains are added in order.
// The returned mol has the same chain ids as do the input residues.
// 
std::pair<bool, CMMDBManager *>
coot::util::create_mmdbmanager_from_residue_vector(const std::vector<CResidue *> &res_vec) { 
  
   // So, first make a vector of residue sets, one residue set for each chain.
   std::vector<coot::util::chain_id_residue_vec_helper_t> residues_of_chain;
   
   for (unsigned int i=0; i<res_vec.size(); i++) { 
      std::string chain_id = res_vec[i]->GetChainID();
      
      // is chain_id already in residues_of_chain?  Do it in line here
      //
      bool found = 0;
      for (unsigned int ich=0; ich<residues_of_chain.size(); ich++) { 
	 if (residues_of_chain[ich].chain_id == chain_id) { 
	    found = 1;
	    break;
	 }
      }
      if (! found) { 
	 coot::util::chain_id_residue_vec_helper_t chir;
	 chir.chain_id = chain_id;
	 residues_of_chain.push_back(chir);
      } 
   }

   // now residues_of_chain is full of containers that have the chain_id specified.
   for (unsigned int i=0; i<res_vec.size(); i++) { 
      std::string chain_id = res_vec[i]->GetChainID();
      for (unsigned int ich=0; ich<residues_of_chain.size(); ich++) { 
	 if (residues_of_chain[ich].chain_id == chain_id) { 
	    residues_of_chain[ich].residues.push_back(res_vec[i]);
	    break;
	 } 
      }
   }
   
   for (unsigned int ich=0; ich<residues_of_chain.size(); ich++) { 
      residues_of_chain[ich].sort_residues();
   }

   std::sort(residues_of_chain.begin(), residues_of_chain.end());

   CModel *model_p = new CModel;
   for (unsigned int ich=0; ich<residues_of_chain.size(); ich++) { 
      CChain *chain_p = new CChain;
      chain_p->SetChainID(residues_of_chain[ich].chain_id.c_str());
      for (unsigned int ires=0; ires<residues_of_chain[ich].residues.size(); ires++) { 
	 CResidue *residue_p = 
            coot::util::deep_copy_this_residue(residues_of_chain[ich].residues[ires]);
	 chain_p->AddResidue(residue_p);
      }
      model_p->AddChain(chain_p);
   }

   CMMDBManager *mol = new CMMDBManager;
   mol->AddModel(model_p);

   return std::pair<bool, CMMDBManager *> (1, mol);

}


CMMDBManager *
coot::util::create_mmdbmanager_from_residue_specs(std::vector<coot::residue_spec_t> &r1,
						  CMMDBManager *mol) {
   std::vector<CResidue *> residues;

   for (unsigned int ires=0; ires<r1.size(); ires++) {
      CResidue *res = coot::util::get_residue(r1[ires], mol);
      if (res) {
	 residues.push_back(res);
      }
   }
   CMMDBManager *new_mol = coot::util::create_mmdbmanager_from_residue_vector(residues).second;
   return new_mol;
}



// a new residue for each point
CMMDBManager *coot::util::create_mmdbmanager_from_points(const std::vector<clipper::Coord_orth> &pts) {

   CMMDBManager *new_mol = new CMMDBManager;
   CModel *model_p = new CModel;
   CChain *chain_p = new CChain;
   chain_p->SetChainID("A");

   for (unsigned int i=0; i<pts.size(); i++) {
      CAtom *at = new CAtom;
      at->SetCoordinates(pts[i].x(), pts[i].y(), pts[i].z(), 1.0, 20.0);
      at->SetAtomName(" CA ");
      at->SetElementName(" C");
      CResidue *residue_p = new CResidue;
      residue_p->SetResName("ALA");
      residue_p->seqNum = i;
      residue_p->AddAtom(at);
      chain_p->AddResidue(residue_p);
   }
   model_p->AddChain(chain_p);
   new_mol->AddModel(model_p);
   return new_mol;
}

void 
coot::util::chain_id_residue_vec_helper_t::sort_residues() { 

   std::sort(residues.begin(), residues.end(), 
	     coot::util::chain_id_residue_vec_helper_t::residues_sort_func);
}


bool
coot::util::chain_id_residue_vec_helper_t::operator<(const chain_id_residue_vec_helper_t &c) const { 

   return (chain_id < c.chain_id);
}

// static 
bool
coot::util::chain_id_residue_vec_helper_t::residues_sort_func(CResidue *first, CResidue *second) { 
   
   if (first->GetSeqNum() < second->GetSeqNum()) { 
      return 1;
   } else { 
      if (first->GetSeqNum() > second->GetSeqNum()) { 
	 return 0; 
      } else { 
	 std::string inscode_1 = first->GetInsCode();
	 std::string inscode_2 = second->GetInsCode(); 
	 if (inscode_1 < inscode_2) 
	    return 1; 
	 else 
	    return 0;
      }
   }
   return 0; // not reached.
}


bool
coot::mol_has_symmetry(CMMDBManager *mol) {
   mat44 test_mat;
   int i_symm_err = mol->GetTMatrix(test_mat, 0, 0, 0, 0);
   return i_symm_err;
} 



// We don't mess with the chain ids, give as we get, but also
// return the handle for the atom index transfer.
std::pair<CMMDBManager *, int>
coot::util::create_mmdbmanager_from_mmdbmanager(CMMDBManager *mol_in) { 

   CMMDBManager *residues_mol = new CMMDBManager;
   int atom_index_handle = residues_mol->RegisterUDInteger(UDR_ATOM, "mol's atom index");
   int afix_handle_orig = mol_in->GetUDDHandle(UDR_ATOM, "shelx afix");
   int afix_handle_new_mol = -1;
   if (afix_handle_orig >= 0)
      afix_handle_new_mol = residues_mol->RegisterUDInteger(UDR_ATOM, "shelx afix");

   std::string altconf = ""; // dummy
   short int whole_res_flag = 1;

   for(int imod = 1; imod<=mol_in->GetNumberOfModels(); imod++) {
     int imod = 1;
     CModel *model_p = mol_in->GetModel(imod);
     CModel *new_model_p = new CModel;
     int nchains = model_p->GetNumberOfChains();
     for (int ichain=0; ichain<nchains; ichain++) {
       CChain *chain_p = model_p->GetChain(ichain);
       CChain *new_chain_p = new CChain;
       new_chain_p->SetChainID(chain_p->GetChainID());
       int nres = chain_p->GetNumberOfResidues();
       for (int ires=0; ires<nres; ires++) { 
	 CResidue *residue_p = chain_p->GetResidue(ires);
	 CResidue *r = coot::util::deep_copy_this_residue_with_atom_index_and_afix_transfer(mol_in, residue_p, altconf, whole_res_flag, atom_index_handle, afix_handle_new_mol);
	 new_chain_p->AddResidue(r);
       }
       new_model_p->AddChain(new_chain_p);
     }
     residues_mol->AddModel(new_model_p);
   }

   return std::pair<CMMDBManager *, int> (residues_mol, atom_index_handle);
} 


// ignore atom index transfer
CMMDBManager *
coot::util::create_mmdbmanager_from_atom_selection(CMMDBManager *orig_mol,
						   int SelectionHandle,
						   bool invert_seletion) {

   if (invert_seletion)
      return coot::util::create_mmdbmanager_from_inverted_atom_selection(orig_mol,
									 SelectionHandle);
   else
      return coot::util::create_mmdbmanager_from_atom_selection_straight(orig_mol,
									 SelectionHandle);
}
   
// ignore atom index transfer
CMMDBManager *
coot::util::create_mmdbmanager_from_atom_selection_straight(CMMDBManager *orig_mol,
							    int SelectionHandle) { 
   CMMDBManager *atoms_mol = new CMMDBManager;
   CModel *model = new CModel;
   atoms_mol->AddModel(model);

   PCAtom *atoms = 0;
   int n_selected_atoms;
   orig_mol->GetSelIndex(SelectionHandle, atoms, n_selected_atoms);
   
   for (int iatom=0; iatom<n_selected_atoms; iatom++) {

      int atom_residue_selection_handle = atoms_mol->NewSelection();
      CAtom *at = atoms[iatom];

      // Does this new atoms residue exist in new atoms_mol?
      // if yes,
      //   add a copy of this atom to that residue
      // else
      //   does this new atom's chain exist in new atoms_mol?
      //   if yes
      //      create a new residue in that chain
      //      add a copy of this atom to that new residue
      //   else
      //      create a new chain
      //      create a new residue in that chain
      //      add a copy of this atom to that new residue
      // 
      atoms_mol->Select (atom_residue_selection_handle,
			STYPE_RESIDUE, 0, // .. TYPE, iModel
			at->GetChainID(), // Chain(s)
			at->GetSeqNum(), at->GetInsCode(),  // starting res
			at->GetSeqNum(), at->GetInsCode(),  // ending res
			at->GetResName(),  // residue name
			"*",  // Residue must contain this atom name?
			"*",  // Residue must contain this Element?
			"*",  // altLocs
			SKEY_NEW // selection key
			);
      PCResidue *sel_residues;
      int n_sel_residues;
      atoms_mol->GetSelIndex(atom_residue_selection_handle,
			     sel_residues, n_sel_residues);

      if (n_sel_residues > 0) {
	 CAtom *new_atom = new CAtom;
	 new_atom->Copy(at);
	 sel_residues[0]->AddAtom(new_atom);
      } else {
	 // residue was not found!
	 // Now we select on the chain
	 int atom_chain_selection_handle = atoms_mol->NewSelection();
	 atoms_mol->Select (atom_chain_selection_handle,
			   STYPE_CHAIN, 0, // .. TYPE, iModel
			    at->GetChainID(), // Chain(s)
			    ANY_RES, "*",  // starting res
			    ANY_RES, "*",  // ending res
			    "*",  // residue name
			    "*",  // Residue must contain this atom name?
			    "*",  // Residue must contain this Element?
			    "*",  // altLocs
			    SKEY_NEW // selection key
			    );
	 PCChain *sel_chains;
	 int n_sel_chains;
	 atoms_mol->GetSelIndex(atom_chain_selection_handle, sel_chains,
				n_sel_chains);

	 if (n_sel_chains > 0) {
	    CResidue *residue = new CResidue(sel_chains[0],
					     at->GetResName(),
					     at->GetSeqNum(),
					     at->GetInsCode());
	    CAtom *new_atom = new CAtom; 
	    new_atom->Copy(at);
	    residue->AddAtom(new_atom);
	 } else {
	    // There was not even a chain...
	    CChain *chain = new CChain(model, at->GetChainID());
	    CResidue *residue = new CResidue(chain,
					     at->GetResName(),
					     at->GetSeqNum(),
					     at->GetInsCode());
	    CAtom *new_atom = new CAtom; 
	    new_atom->Copy(at);
	    residue->AddAtom(new_atom);
// 	    std::cout << "DEBUG:: straight: added atom (and res, chains) :"
// 		      << new_atom->GetAtomName()
// 		      << ": to "
// 		      << residue->GetChainID() << " " 
// 		      << residue->GetSeqNum()  << " " 
// 		      << residue->GetResName() << " " 
// 		      << std::endl;
// 	    std::cout << "DEBUG:: straight: added    residue " << residue->GetSeqNum()
// 		      << std::endl;
// 	    std::cout << "DEBUG:: straight: added       chain " << chain->GetChainID()
// 		      << std::endl;
	 } 
      }
      atoms_mol->DeleteSelection(atom_residue_selection_handle);
   }

   realtype a[6];
   realtype vol;
   int orthcode;
   orig_mol->GetCell(a[0], a[1], a[2], a[3], a[4], a[5], vol, orthcode);
   atoms_mol->SetCell(a[0], a[1], a[2], a[3], a[4], a[5]);
   cpstr sg = orig_mol->GetSpaceGroup();
   if (sg) { 
     atoms_mol->SetSpaceGroup(sg);
   }
   atoms_mol->FinishStructEdit();
   return atoms_mol;
}

// Beware This destroys (inverts) the atom selection as passed.
CMMDBManager *
coot::util::create_mmdbmanager_from_inverted_atom_selection(CMMDBManager *orig_mol,
							    int SelectionHandle) {

   // The idea here is that we want to have a selection that is
   // logical NOT of the SelectionHandle selection.
   // 
   // So we need to select everything in orig_mol and then use
   // SKEY_XOR to get a selection that is the NOT ofthe
   // SelectionHandle selection.

   orig_mol->Select(SelectionHandle, STYPE_ATOM, 0, "*", ANY_RES, "*", ANY_RES, "*",
		    "*", "*", "*", "*", SKEY_XOR);
   CMMDBManager *new_mol =
      coot::util::create_mmdbmanager_from_atom_selection(orig_mol, SelectionHandle);
   return new_mol;
}


// ignore atom index transfer, return NULL on error.
CMMDBManager *
coot::util::create_mmdbmanager_from_residue(CMMDBManager *orig_mol,
					    CResidue *res) {

   CResidue *r = coot::util::deep_copy_this_residue(res);
   CMMDBManager *mol = new CMMDBManager;
   CModel *model_p = new CModel;
   CChain *chain_p = new CChain;
   chain_p->AddResidue(r);
   model_p->AddChain(chain_p);
   mol->AddModel(model_p);
   if (mol) {
      chain_p->SetChainID(res->GetChainID());
   } else {
      chain_p->SetChainID("");
   } 
   return mol;
}

CMMDBManager *
coot::util::create_mmdbmanager_from_atom(CAtom *at) {
   CResidue *res = new CResidue;
   res->AddAtom(at);
   CChain *chain_p = new CChain;
   chain_p->AddResidue(res);
   chain_p->SetChainID("A");
   CModel *model_p = new CModel;
   model_p->AddChain(chain_p);
   CMMDBManager *mol = new CMMDBManager;
   mol->AddModel(model_p);
   return mol;
} 


// Note, we also create a chain and add this residue to that chain.
// We do this so that we have a holder for the segid.
// 
// whole_residue_flag: only copy atoms that are either in this altLoc,
// or has an altLoc of "".
// 
CResidue *
coot::util::deep_copy_this_residue_add_chain(CResidue *residue,
					     const std::string &altconf,
					     bool whole_residue_flag,
					     bool attach_to_new_chain_flag) {

   // Horrible casting to CResidue because GetSeqNum and GetAtomTable
   // are not const functions.
   // 
   CResidue *rres = new CResidue;
   CChain   *chain_p = NULL;
   if (attach_to_new_chain_flag) { 
      chain_p = new CChain;
      chain_p->SetChainID(residue->GetChainID());
   }
   rres->seqNum = residue->GetSeqNum();
   strcpy(rres->name, residue->name);

   PPCAtom residue_atoms;
   int nResidueAtoms;
   ((CResidue *)residue)->GetAtomTable(residue_atoms, nResidueAtoms);
   CAtom *atom_p;
   
   for(int iat=0; iat<nResidueAtoms; iat++) {
      if (! residue_atoms[iat]->isTer()) { 
	 std::string this_atom_alt_loc(residue_atoms[iat]->altLoc);
	 if (whole_residue_flag ||
	     this_atom_alt_loc  == altconf || this_atom_alt_loc == "") { 
	    atom_p = new CAtom;
	    atom_p->Copy(residue_atoms[iat]);
	    rres->AddAtom(atom_p);
	 }
      }
   }
   if (attach_to_new_chain_flag)
      chain_p->AddResidue(rres);
   return rres;
}

CResidue *
coot::util::deep_copy_this_residue(CResidue *residue) { 

   CResidue *rres = new CResidue;
   rres->seqNum = residue->GetSeqNum();
   strcpy(rres->name, residue->name);

   PPCAtom residue_atoms;
   int nResidueAtoms;
   ((CResidue *)residue)->GetAtomTable(residue_atoms, nResidueAtoms);
   CAtom *atom_p;
   
   for(int iat=0; iat<nResidueAtoms; iat++) {
      if (! residue_atoms[iat]->isTer()) { 
	 atom_p = new CAtom;
	 atom_p->Copy(residue_atoms[iat]);
	 rres->AddAtom(atom_p);
      }
   }
   return rres;
}


// Note, we also create a chain and add this residue to that chain.
// We do this so that we have a holder for the segid.
// 
// whole_residue_flag: only copy atoms that are either in this altLoc,
// or has an altLoc of "".
//
//
// Don't transfer afix number if afix_udd_handle (for the new atoms) is -1;
CResidue *
coot::util::deep_copy_this_residue_with_atom_index_and_afix_transfer(CMMDBManager *std_mol, 
								     const CResidue *residue,
								     const std::string &altconf,
								     short int whole_residue_flag,
								     int atom_index_udd_handle,
								     int afix_udd_handle) {

   // Note we pass the atom_index_udd_handle from the new residues
   // molecule and get the old one by using GetUDDHandle on the
   // old/std mol.

   // Horrible casting to CResidue because GetSeqNum and GetAtomTable
   // are not const functions.
   // 
   CResidue *rres = new CResidue;
   CChain   *chain_p = new CChain;
   chain_p->SetChainID(((CResidue *)residue)->GetChainID());
   rres->seqNum = ((CResidue *)residue)->GetSeqNum();
   strcpy(rres->name, residue->name);

   PPCAtom residue_atoms;
   int nResidueAtoms;
   ((CResidue *)residue)->GetAtomTable(residue_atoms, nResidueAtoms);
   CAtom *atom_p;

   // We are not passed the handle, we have to look it up.
   int mol_atom_index_handle = std_mol->GetUDDHandle(UDR_ATOM, "atom index");
   int mol_afix_handle = -1;
   int afix_number;
   if (afix_udd_handle >= 0) {
      mol_afix_handle = std_mol->GetUDDHandle(UDR_ATOM, "shelx afix");
   }
      
   for(int iat=0; iat<nResidueAtoms; iat++) {
      std::string this_atom_alt_loc(residue_atoms[iat]->altLoc);
      if (whole_residue_flag ||
	  this_atom_alt_loc  == altconf || this_atom_alt_loc == "") { 
	 atom_p = new CAtom;
	 atom_p->Copy(residue_atoms[iat]);
	 int mol_atom_index = -1;
	 residue_atoms[iat]->GetUDData(mol_atom_index_handle, mol_atom_index);
	 atom_p->PutUDData(atom_index_udd_handle, mol_atom_index);
	 // and shelx afix data:
	 if (mol_afix_handle >= 0) {
	    if (residue_atoms[iat]->GetUDData(mol_afix_handle, afix_number) == UDDATA_Ok) { 
// 	       std::cout << "DEBUG:: Transfering udd afix: " << afix_number
// 			 << " using std mol_udd_handle " << mol_afix_handle
// 			 << " to new atom with afix_udd_handle: " << afix_udd_handle
// 			 << std::endl;
	       atom_p->PutUDData(afix_udd_handle, afix_number);
// 	    } else {
// 	       std::cout << "Ooops! Can get afix number from original atom"
// 			 << std::endl;
	    } 
	 } 
	 rres->AddAtom(atom_p);
      }
   }
   chain_p->AddResidue(rres);
   return rres;
}

CResidue *coot::util::copy_and_delete_hydrogens(CResidue *residue_in) {

   CResidue *copy = coot::util::deep_copy_this_residue(residue_in);
   PPCAtom residue_atoms;
   int nResidueAtoms;
   copy->GetAtomTable(residue_atoms, nResidueAtoms);

   for(int i=0; i<nResidueAtoms; i++) {
      std::string element(residue_atoms[i]->element);
      if (element == " H" || element == " D") {
	 copy->DeleteAtom(i);
      }
   }
   copy->TrimAtomTable();
   return copy;
} 



// transform atoms in residue
void
coot::util::transform_atoms(CResidue *res, const clipper::RTop_orth &rtop) {

   PPCAtom residue_atoms;
   int natoms;
   clipper::Coord_orth co;
   clipper::Coord_orth trans_pos; 
   res->GetAtomTable(residue_atoms, natoms);
   for (int iatom=0; iatom<natoms; iatom++) {
      co = clipper::Coord_orth(residue_atoms[iatom]->x, 
			       residue_atoms[iatom]->y, 
			       residue_atoms[iatom]->z);
      trans_pos = co.transform(rtop);
      residue_atoms[iatom]->x = trans_pos.x();
      residue_atoms[iatom]->y = trans_pos.y();
      residue_atoms[iatom]->z = trans_pos.z();
   }
}

// transform all the atom in mol
void
coot::util::transform_mol(CMMDBManager *mol, const clipper::RTop_orth &rtop) {

   int n_models = mol->GetNumberOfModels();
   for (int imod=1; imod<=n_models; imod++) { 
      
      CModel *model_p = mol->GetModel(imod);
      CChain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 PCResidue residue_p;
	 CAtom *at;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    int n_atoms = residue_p->GetNumberOfAtoms();
	    for (int iat=0; iat<n_atoms; iat++) {
	       at = residue_p->GetAtom(iat);
	       clipper::Coord_orth co(at->x, at->y, at->z);
	       clipper::Coord_orth trans_pos = co.transform(rtop);
	       at->x = trans_pos.x();
	       at->y = trans_pos.y();
	       at->z = trans_pos.z();
	    }
	 }
      }
   }
} 

// a function to find the previous (next) residue.  Find residue
// by previous (next) serial number.
// Return NULL if prevous (next) resiude not found.
// 
CResidue *
coot::util::previous_residue(CResidue *this_residue) {

   CResidue *prev_res = NULL;
   CChain *chain_p = this_residue->GetChain();
   int nres = chain_p->GetNumberOfResidues();
   CResidue *residue_p;
   for (int ires=0; ires<nres; ires++) { 
      residue_p = chain_p->GetResidue(ires);
      if (this_residue == residue_p) {
	 if (ires>0)
	    prev_res = chain_p->GetResidue(ires-1);
	 break;
      } 
   }
   return prev_res;
}

// a function to find the previous (next) residue.  Find residue
// by previous (next) serial number.
// Return NULL if prevous (next) resiude not found.
// 
CResidue *
coot::util::next_residue(CResidue *this_residue) {

   CResidue *prev_res = NULL;
   CChain *chain_p = this_residue->GetChain();
   int nres = chain_p->GetNumberOfResidues();
   CResidue *residue_p;
   for (int ires=0; ires<nres; ires++) { 
      residue_p = chain_p->GetResidue(ires);
      if (this_residue == residue_p) {
	 if (ires < (nres-1))
	    prev_res = chain_p->GetResidue(ires+1);
	 break;
      } 
   }
   return prev_res;
}


CAtom *
coot::util::intelligent_this_residue_mmdb_atom(CResidue *res_p) {

   PCAtom *residue_atoms;
   int nResidueAtoms;
   
   res_p->GetAtomTable(residue_atoms, nResidueAtoms);
   for (int i=0; i<nResidueAtoms; i++) {
      std::string atom_name(residue_atoms[i]->name);
      if (atom_name == " CA ") {
	 return residue_atoms[i];
      }
   }

   if (nResidueAtoms > 0) {
      return residue_atoms[0];
   }

   // failure
   return NULL;

}

float
coot::util::occupancy_sum(PCAtom *atoms, int n_atoms) {

   float os = 0.0;

   for (int i=0; i<n_atoms; i++) {
      os += atoms[i]->occupancy;
   }
   return os;
} 

short int
coot::util::is_nucleotide(CResidue *residue_p) {
   
   short int nuc = 0;

   if (residue_p) { 
      std::string type(residue_p->name); // all spaces cut

      if (type == "Ad") { 
	 nuc = 1;
      } else { 
	 if (type == "Cd") { 
	    nuc = 1;
	 } else { 
	    if (type == "Gd") { 
	       nuc = 1;
	    } else { 
	       if (type == "Td") { 
		  nuc = 1;
	       } else { 
		  if (type == "Ar") { 
		     nuc = 1;
		  } else { 
		     if (type == "Cr") { 
			nuc = 1;
		     } else { 
			if (type == "Gr") { 
			   nuc = 1;
			} else { 
			   if (type == "Ur") { 
			      nuc = 1;
			   } else { 
			      if (type == "DG") { 
				 nuc = 1;
			      } else { 
				 if (type == "DC") { 
				    nuc = 1;
				 } else { 
				    if (type == "DA") { 
				       nuc = 1;
				    } else { 
				       if (type == "DU") { 
					  nuc = 1;
				       } else { 
					  if (type == "DT") { 
					     nuc = 1;
					  } else { 
					     if (type == "DI") { 
						nuc = 1;
					     } else {
						if (type == "Ud") {  // happens sometimes, e.g. brna.pdb
						   nuc = 1;
						} else {
						   if (type == "Tr") {  // not very likely
						      nuc = 1;
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
		  }
	       }
	    }
	 }
      }
   }
   return nuc;
}

// test for presence of O2'
bool
coot::util::nucleotide_is_DNA(CResidue *r) {

   bool has_o2_prime = 0;
   PPCAtom residue_atoms = NULL;
   int n_residue_atoms;
   r->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int i=0; i<n_residue_atoms; i++) {
      std::string atom_name = residue_atoms[i]->name;
      if (atom_name == " O2'") {
	 has_o2_prime = 1;
	 break;
      }
      if (atom_name == " O2*") {
	 has_o2_prime = 1;
	 break;
      }
   }

   if (has_o2_prime)
      return 0;
   else
      return 1;
}



// Return NULL on no such chain:
// 
CChain *
coot::util::chain_only_of_type(CMMDBManager *mol, const std::string &residue_type) {

   CChain *chain_p = NULL;
   CChain *single_type_chain_p = NULL;
   CModel *model_p = mol->GetModel(1);
   CResidue *residue_p;
   
   if (model_p) { 
      int nchains = model_p->GetNumberOfChains();
      for (int ich=0; ich<nchains; ich++) {
	 chain_p = model_p->GetChain(ich);
	 int nres = chain_p->GetNumberOfResidues();
	 short int all_same_type_flag = 1; 
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    std::string resname(residue_p->name);
	    if (! (resname == residue_type)) {
	       all_same_type_flag = 0;
	       break;
	    }
	 }
	 if (all_same_type_flag) { 
	    single_type_chain_p = chain_p;
	    break;
	 }
      }
   } 
   return single_type_chain_p;
} 


std::string
coot::util::three_letter_to_one_letter(const std::string &resname) {

   std::string n;
   char r[10];
   //    std::cout << "DEBUG:: resname length: " << resname << " " << resname.length() << std::endl;

   short int done_locally = 0;
   if (resname.length() == 2) {
      if ((resname[1] == 'd') ||
	  (resname[1] == 'r')) {
	 n = resname.substr(0, 1);
	 done_locally = 1;
      }
   }
   
   if (! done_locally) { 
      Get1LetterCode(resname.c_str(), r);
      n = r[0];
   }
   return n;
}

std::string
coot::util::three_letter_to_one_letter_with_specials(const std::string &resname) {

   std::string n;
   if (resname == "HOH") {
      n = "~";
   } else {
      n = coot::util::three_letter_to_one_letter(resname);
   } 
   return n;
}



std::pair<clipper::Coord_orth, clipper::Coord_orth>
coot::util::extents(CMMDBManager *mol) {


   int selHnd = mol->NewSelection();
   mol->SelectAtoms(selHnd, 0, "*", ANY_RES, "*", ANY_RES, "*",
		    "*", "*", "*", "*");
   std::pair<clipper::Coord_orth, clipper::Coord_orth> p =
      coot::util::extents(mol, selHnd);
   mol->DeleteSelection(selHnd);
   return p;
}

std::pair<clipper::Coord_orth, clipper::Coord_orth>
coot::util::extents(CMMDBManager *mol,
		     int SelectionHandle) {

   std::pair<clipper::Coord_orth, clipper::Coord_orth> p;

   PCAtom *atoms = NULL;
   int n_selected_atoms;
   mol->GetSelIndex(SelectionHandle, atoms, n_selected_atoms);
   float most_x = -99999;
   float most_y = -99999;
   float most_z = -99999;
   float least_x = 99999;
   float least_y = 99999;
   float least_z = 99999;
   for (int i=0; i<n_selected_atoms; i++) {
      if (atoms[i]->x > most_x) most_x = atoms[i]->x;
      if (atoms[i]->y > most_y) most_y = atoms[i]->y;
      if (atoms[i]->z > most_z) most_z = atoms[i]->z;

      if (atoms[i]->x < least_x) least_x = atoms[i]->x;
      if (atoms[i]->y < least_y) least_y = atoms[i]->y;
      if (atoms[i]->z < least_z) least_z = atoms[i]->z;
   }

   clipper::Coord_orth p1( most_x,  most_y,  most_z);
   clipper::Coord_orth p2(least_x, least_y, least_z);

   return std::pair<clipper::Coord_orth, clipper::Coord_orth> (p2, p1);
}


// pair.second = 0 for failure
// pair.first  = 1 for success
// 
std::pair<clipper::RTop_orth, short int>
coot::util::get_ori_to_this_res(CResidue *res) {

   std::pair<clipper::RTop_orth, short int>  pair;
   // clipper::RTop_orth op;

   PPCAtom residue_atoms;
   int nResidueAtoms;
   res->GetAtomTable(residue_atoms, nResidueAtoms);
   if (nResidueAtoms == 0) {
      std::cout << "ERROR:: something broken in atom residue selection before ";
      std::cout << "get_ori_to_this_res, 0 atoms in given residue" << std::endl;
   } else {

      clipper::Coord_orth ca(0,0,0), c(0,0,0), n(0,0,0);
      int n_found = 0;
      bool found_ca = 0;
      bool found_c  = 0;
      bool found_n  = 0;
      
      for(int iat=0; iat<nResidueAtoms; iat++) {
	 std::string atom_name = residue_atoms[iat]->name;
	 if (atom_name == " CA ") {
	    n_found++;
	    found_ca = 1;
	    ca = clipper::Coord_orth(residue_atoms[iat]->x,
				     residue_atoms[iat]->y,
				     residue_atoms[iat]->z);
	 }
	 if (atom_name == " C  ") {
	    n_found++;
	    found_c = 1;
	    c  = clipper::Coord_orth(residue_atoms[iat]->x,
				     residue_atoms[iat]->y,
				     residue_atoms[iat]->z);
	 }
	 if (atom_name == " N  ") {
	    n_found++;
	    found_n = 1;
	    n  = clipper::Coord_orth(residue_atoms[iat]->x,
				     residue_atoms[iat]->y,
				     residue_atoms[iat]->z);
	 }
      }

      if (n_found != 3) {
	 std::cout << "DISASTER! Not all necessary atoms found in residue ";
	 std::cout << res->GetChainID() << " " << res->GetSeqNum() << std::endl;
	 if (found_ca == 0)
	    std::cout << "    CA is missing " << std::endl;
	 if (found_c == 0)
	    std::cout << "    C is missing " << std::endl;
	 if (found_n == 0)
	    std::cout << "    N is missing " << std::endl;
	 pair.second = 0; // failure
	 return pair;
      }

      // testing
//       ca = clipper::Coord_orth(0,0,0);
//       n  = clipper::Coord_orth(0.87, 0, 1.23);
//       c  = clipper::Coord_orth(0.83, 0, -1.18);
//       clipper::Coord_orth cb(-1.03, -1.11, 0);

      
      // now get the vectors of the orientation:
      //
      clipper::Coord_orth can_unit = clipper::Coord_orth((n - ca).unit());
      clipper::Coord_orth cac_unit = clipper::Coord_orth((c - ca).unit());

      clipper::Coord_orth bisector ((can_unit + cac_unit).unit());
      clipper::Coord_orth diff_unit((can_unit - cac_unit).unit());

      clipper::Coord_orth cross_prod(clipper::Coord_orth::cross(diff_unit,bisector));
      clipper::Coord_orth cpu = clipper::Coord_orth(cross_prod.unit());

//       std::cout << "bisector   " << bisector.format() << std::endl;
//       std::cout << "diff_unit  " << diff_unit.format() << std::endl;
//       std::cout << "cross prod " << cross_prod.format() << std::endl;

      // bisector   -> new x axis
      // diff_unit  -> new z axis
      // cross_prod -> new y axis

      clipper::Mat33<double> m(  bisector.x(),   bisector.y(),   bisector.z(),
			              cpu.x(),        cpu.y(),        cpu.z(),
			        diff_unit.x(),  diff_unit.y(),  diff_unit.z());

      pair.first = clipper::RTop_orth(m.transpose(), ca);
      pair.second = 1;

   }
   return pair;
} 


// For sequence/sequence alignment


// Take into account the insertion code too:
std::vector<std::pair<CResidue *, int> >
coot::util::sort_residues_by_seqno(PCResidue *residues,
				   int nResidues) {

   std::vector<std::pair<CResidue *, int> > v;

   // construct the vector
   for (int i=0; i<nResidues; i++)
      v.push_back(std::pair<CResidue *, int> (residues[i], i));

//    // test
//    std::vector<int> vi;
//    vi.push_back(10);
//    vi.push_back(1);
//    vi.push_back(0);

//    std::vector<int>::iterator start = vi.begin();
//    std::vector<int>::iterator end   = vi.end();

//    sort(start, end);

   // and now with our lovely data (not int)
   // 
   std::vector<std::pair<CResidue *, int> >::iterator start = v.begin();
   std::vector<std::pair<CResidue *, int> >::iterator end   = v.end();
   std::sort(start, end, coot::util::compare_residues);
   // sort(start, end) for things that have implicit comparison function.
      
   return v;
}

bool
coot::util::compare_residues(const std::pair<CResidue *, int> &a,
			     const std::pair<CResidue *, int> &b) {

   int r1 = a.first->GetSeqNum();
   int r2 = b.first->GetSeqNum();

   if (r1 < r2) {
      return 1;
   } else {
      if (r1 > r2) {
	 return 0;
      } else {
	 std::string ins1(a.first->GetInsCode());
	 std::string ins2(b.first->GetInsCode());
	 if (ins1 > ins2) {
	    return 0;
	 } else {
	    return 1; // check
	 }
      }
   }

   return 1;
} 


// Use the results of the above to give us a sequence string:
std::string
coot::util::model_sequence(const std::vector<std::pair<CResidue *, int> > &sa) {

   std::string s;
   char r[10];
   for (unsigned int i=0; i<sa.size(); i++) {
      std::string this_residue = "X";
      std::string res_name = sa[i].first->GetResName();
      if (res_name != "HOH") { 
	 Get1LetterCode(res_name.c_str(), r);
	 this_residue = r[0];
	 s += this_residue;
      }
   }
   return s;
}


// residues with insertion codes
std::vector<CResidue *>
coot::util::residues_with_insertion_codes(CMMDBManager *mol) {

   std::vector<CResidue *> v;
   
   int imod = 1;
      
   CModel *model_p = mol->GetModel(imod);
   CChain *chain_p;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      PCResidue residue_p;
      for (int ires=0; ires<nres; ires++) { 
	 residue_p = chain_p->GetResidue(ires);
	 std::string inscode(residue_p->GetInsCode());
	 if (inscode != "") {
	    v.push_back(residue_p);
	 }
      }
   }
   return v;
}

// return the status too, torsion in radians.
std::pair<short int, double>
coot::util::omega_torsion(CResidue *C_residue, CResidue *N_residue, const std::string &altconf) {

   double omega = 0;
   short int istatus = 0; // initially unset
   PPCAtom C_residue_atoms = NULL;
   int nCResidueAtoms;
   C_residue->GetAtomTable(C_residue_atoms, nCResidueAtoms);
   CAtom *C_residue_CA_atom_p = NULL;
   CAtom *C_residue_C_atom_p = NULL; 

   PPCAtom N_residue_atoms = NULL;
   int nNResidueAtoms;
   N_residue->GetAtomTable(N_residue_atoms, nNResidueAtoms);
   CAtom *N_residue_CA_atom_p = NULL;
   CAtom *N_residue_N_atom_p = NULL;

   for (int i=0; i<nCResidueAtoms; i++) {
      std::string atom_name = C_residue_atoms[i]->name;
      std::string altconf_atom = C_residue_atoms[i]->altLoc;
      if (atom_name == " CA ")
	 if (altconf_atom == altconf)
	    C_residue_CA_atom_p = C_residue_atoms[i];
      if (atom_name == " C  ")
	 if (altconf_atom == altconf)
	    C_residue_C_atom_p = C_residue_atoms[i];
   }

   for (int i=0; i<nNResidueAtoms; i++) {
      std::string atom_name = N_residue_atoms[i]->name;
      std::string altconf_atom = N_residue_atoms[i]->altLoc;
      if (atom_name == " CA ")
	 if (altconf_atom == altconf)
	    N_residue_CA_atom_p = N_residue_atoms[i];
      if (atom_name == " N  ")
	 if (altconf_atom == altconf)
	    N_residue_N_atom_p = N_residue_atoms[i];
   }

   if (C_residue_CA_atom_p && C_residue_C_atom_p && N_residue_N_atom_p && N_residue_CA_atom_p) {
      clipper::Coord_orth ca1(C_residue_CA_atom_p->x,
			      C_residue_CA_atom_p->y,
			      C_residue_CA_atom_p->z);
      clipper::Coord_orth c1(C_residue_C_atom_p->x,
			     C_residue_C_atom_p->y,
			     C_residue_C_atom_p->z);
      clipper::Coord_orth ca2(N_residue_CA_atom_p->x,
			      N_residue_CA_atom_p->y,
			      N_residue_CA_atom_p->z);
      clipper::Coord_orth n2(N_residue_N_atom_p->x,
			     N_residue_N_atom_p->y,
			     N_residue_N_atom_p->z);

      omega = clipper::Coord_orth::torsion(ca1, c1, n2, ca2);
      istatus = 1;
   }
   return std::pair<short int, double> (istatus, omega);
}


// This function is based on mmdb's CATom::MakePDBAtomName(), so we share
// the copyright (with Eugene/EBI).
//
// 
// 
std::string
coot::pad_atom_name(const std::string &atom_id, const std::string &element) { 

   std::string name = atom_id;
   std::string new_name; // = atom_id;

   // std::cout << "DEBUG:: :" << atom_id << ": :" << element << ":" << std::endl;

   if (name.length() == 4) {

      new_name = name;

   } else {
   
      if (element == "") {
	 if (name.length() == 1) {
	    new_name = " ";
	    new_name += name;
	    new_name += "  ";
	 }
      
      } else { // element was defined

	 // consider the following cases of 
	 // atom_id and element:
	 //   "CA"  "C"  -> " CA "
	 //   "AC8" "C"  -> "AC8 " : in this case we post apply
	 //                          the space because the second char of
	 //                          the atom_id is equal to the
	 //                          element
	 //   "NN1" "N"  -> "NN1 " :
	 //   "NP"  "P"  -> "NP  "
	 
	 if (element.length() == 1) {
	    int k = atom_id.length();
	    if (k == 3) {
	       // std::cout << "comparing " << atom_id.substr(1,1) << " and " << element << std::endl;
	       if (atom_id.substr(1,1) == element) {
		  new_name = atom_id + " ";
	       } else {
		  new_name = " " + atom_id;
	       }
	    } else { 
	       // promote the characters one space
	       if (k==2) {
		  // e.g "NP" "P", or "CA" "C"
		  if ((atom_id.substr(1,1) == element) && (element != "H")) {
		     new_name = atom_id + "  ";
		  } else {
		     new_name = " " + atom_id;
		     new_name += " ";
		  }
	       } else {
		  new_name = " ";
		  new_name += atom_id;
		  // fill the rest with spaces
		  new_name += "  ";
	       }
	    }

	 } else {
	    // element was 2 (or more) characters
	    if (element[0] == ' ') { // unusual from dict, but usual for SHELXL reading
	       if (element[1] != name[1]) {
		  // promote the characters one space
		  new_name = " ";
		  new_name += atom_id;
		  int k = atom_id.length();
		  // fill the rest with spaces
		  if (k == 1) 
		     new_name += "  ";
		  if (k == 2)
		     new_name += " ";
	       } else {
		  int k = atom_id.length();
		  new_name = " ";
		  new_name += atom_id;
		  if (k == 1) 
		     new_name += "  ";
		  if (k == 2)
		     new_name += " ";
	       } 
	    } else {
	       // This is the usual case, there is no leading space
	       //
	       // left justify the name and pad with spaces
	       new_name = atom_id;
	       if (atom_id.size() == 1) // never happens?
		  new_name += "   ";
	       if (atom_id.size() == 2)
		  new_name += "  ";
	       if (atom_id.size() == 3)
		  new_name += " ";
	       
	    }
	 }
      }
   }

   // debug
//    CAtom at;
//    at.SetAtomName(atom_id.c_str());
//    at.SetElementName(type_symbol.c_str());
//    at.MakePDBAtomName();
//    if (new_name != std::string(at.name)) {
//       std::cout << "name pad failure, mmdb coot :" << at.name << ": :"
// 		<< new_name << ": for in_atom :" << atom_id.c_str()
// 		<< ": element :" << type_symbol.c_str() << ":" << std::endl;
//    } else {
//       std::cout << "name pad match " << new_name << std::endl;
//    }

//    std::cout << "new atom name :" << new_name << ": from :"
// 	     << atom_id << ": :" << element << ":" << std::endl;
   
   return new_name;
}

std::pair<double, double>
coot::lsq_plane_deviation(const std::vector<clipper::Coord_orth> &v,
			  const clipper::Coord_orth &pt) {

   coot::lsq_plane_info_t lpd(v);
   double val = lpd.plane_deviation(pt);
   double rms = lpd.plane_atoms_rms();
   return std::pair<double, double> (val, rms);
}

coot::lsq_plane_info_t::lsq_plane_info_t(const std::vector<clipper::Coord_orth> &v) {

   int n_atoms = v.size();
   clipper::Coord_orth sum(0,0,0);
   for (int i=0; i<n_atoms; i++)
      sum += v[i];
   double factor = 1/double(n_atoms);
   clipper::Coord_orth midpoint(sum.x()*factor, sum.y()*factor, sum.z()*factor);
   centre_ = midpoint;

   clipper::Matrix<double> mat(3,3);
   for (int i=0; i<n_atoms; i++) {
      mat(0,0) += (v[i].x() - midpoint.x()) * (v[i].x() - midpoint.x());
      mat(1,1) += (v[i].y() - midpoint.y()) * (v[i].y() - midpoint.y());
      mat(2,2) += (v[i].z() - midpoint.z()) * (v[i].z() - midpoint.z());
      mat(0,1) += (v[i].x() - midpoint.x()) * (v[i].y() - midpoint.y());
      mat(0,2) += (v[i].x() - midpoint.x()) * (v[i].z() - midpoint.z());
      mat(1,2) += (v[i].y() - midpoint.y()) * (v[i].z() - midpoint.z());
   }
   mat(1,0) = mat(0,1);
   mat(2,0) = mat(0,2);
   mat(2,1) = mat(1,2);

   if (0) { 
      std::cout << "  mat for eigens: " << std::endl;
      std::cout << "     " << mat(0,0) << "   " << mat(0,1) << "   " << mat(0,2) << std::endl;
      std::cout << "     " << mat(1,0) << "   " << mat(1,1) << "   " << mat(1,2) << std::endl;
      std::cout << "     " << mat(2,0) << "   " << mat(2,1) << "   " << mat(2,2) << std::endl;
   }
   std::vector<double> eigens = mat.eigen(true);
   // Let's now extract the values of a,b,c normalize them
   abcd.resize(4);
   
   abcd[0] = mat(0,0);
   abcd[1] = mat(1,0);
   abcd[2] = mat(2,0);

   if (0) 
      std::cout << " abcd - pre-values "
		<< abcd[0] << " "
		<< abcd[1] << " "
		<< abcd[2] << " "
		<< std::endl;
   
   double sqsum = 1e-20;
   
   for (int i=0; i<3; i++)
      sqsum += abcd[i] * abcd[i];
   for (int i=0; i<3; i++)
      abcd[i] /= sqsum;
   
   // set D, recall di = Axi+Byi+Czi-D, so when
   // xi = x_cen, yi = y_cen, zi = z_cen, d is 0,
   // so we can set D.
   // 
   abcd[3] = abcd[0]*midpoint.x() + abcd[1]*midpoint.y() + abcd[2]*midpoint.z();

   if (0) 
      std::cout << " abcd "
		<< abcd[0] << " "
		<< abcd[1] << " "
		<< abcd[2] << " "
		<< abcd[3] << std::endl;

   double var = 0;
   for (unsigned int i_plane_at=0; i_plane_at<v.size(); i_plane_at++) {
      double d =
	 abcd[0]*v[i_plane_at].x() +
	 abcd[1]*v[i_plane_at].y() +
	 abcd[2]*v[i_plane_at].z() - abcd[3];
      var += d*d;
   }
   rms = 0;
   if (v.size() > 0)
      rms = sqrt(var/double(v.size()));

} 

bool
coot::compare_atom_specs_user_float(const coot::atom_spec_t &a1, const coot::atom_spec_t &a2) {

   return a1.float_user_data < a2.float_user_data ? 1 : 0;

} 


// For use with interesting-things-gui, make the list argument
// from a vector of atom specs.
// 
// use the user data in the atom spec to give us the molecule number
// and the button label
// 
std::string
coot::util::interesting_things_list(const std::vector<atom_spec_t> &v) {

#ifdef USE_GUILE
   // e.g. (list) for empty v
   // (list (list "button label" imol-no chain-id resno atom-name)
   //       (list "button label" imol-no chain-id resno atom-name)
   // )

   std::string r = " (list ";

   for (unsigned int i=0; i<v.size(); i++) {

      std::string atom_str("\"");
      atom_str += v[i].chain;
      atom_str += "\" ";
      atom_str += int_to_string(v[i].resno);
      atom_str += " \"";
      atom_str += v[i].insertion_code;
      atom_str += "\" \"";
      atom_str += v[i].atom_name;
      atom_str += "\" \"";
      atom_str += v[i].alt_conf;
      atom_str += " \"";

      std::string button_label("Clash gap: ");
      button_label += float_to_string(v[i].float_user_data);
      button_label += " : ";
      button_label += v[i].chain;
      button_label += " ";
      button_label += int_to_string(v[i].resno);
      button_label += " ";
      if (v[i].insertion_code != "") {
	 button_label += v[i].insertion_code;
 	 button_label += " ";
      }
      button_label += v[i].atom_name;
      if (v[i].alt_conf != "") {
	 button_label += ",";
	 button_label += v[i].alt_conf;
	 button_label += " ";
      }

      std::string s = "(list ";
      s += single_quote(button_label);
      s += " ";
      s += int_to_string(v[i].int_user_data);
      s += " ";
      s += atom_str;
      s += ")\n";
      
      r += s;
   }

   r += ")";
   return r;
#else
#ifdef USE_PYTHON
// BL says:: we want to have [] lists in python, separated by commas (,)
   // e.g. [] for empty v
   // [["button label",imol-no,chain-id,resno,atom-name],
   //  ["button label",imol-no,chain-id,resno,atom-name]
   // ]

   std::string r = "[";

   for (unsigned int i=0; i<v.size(); i++) {

      std::string atom_str("\"");
      atom_str += v[i].chain;
      atom_str += "\",";
      atom_str += int_to_string(v[i].resno);
      atom_str += ",\"";
      atom_str += v[i].insertion_code;
      atom_str += "\",\"";
      atom_str += v[i].atom_name;
      atom_str += "\",\"";
      atom_str += v[i].alt_conf;
      atom_str += " \"";

      std::string button_label("Clash gap: ");
      button_label += float_to_string(v[i].float_user_data);
      button_label += " : ";
      button_label += v[i].chain;
      button_label += " ";
      button_label += int_to_string(v[i].resno);
      button_label += " ";
      if (v[i].insertion_code != "") {
         button_label += v[i].insertion_code;
         button_label += " ";
      }
      button_label += v[i].atom_name;
      if (v[i].alt_conf != "") {
         button_label += ",";
         button_label += v[i].alt_conf;
         button_label += " ";
      }

      std::string s = "[";
      s += single_quote(button_label);
      s += ",";
      s += int_to_string(v[i].int_user_data);
      s += ",";
      s += atom_str;
      s += "],\n";

      r += s;
   }

   r += "]";
   return r;
#else    
   return "";
#endif // PYTHON
   
#endif // GUILE
}

// we shall have an extra python version (as well)
std::string
coot::util::interesting_things_list_py(const std::vector<atom_spec_t> &v) {

#ifdef USE_PYTHON
   // BL says:: we want to have [] lists in python, separated by commas (,)
   // e.g. [] for empty v
   // [["button label",imol-no,chain-id,resno,atom-name],
   //  ["button label",imol-no,chain-id,resno,atom-name]
   // ]

   std::string r = "[";

   for (unsigned int i=0; i<v.size(); i++) {

      std::string atom_str("\"");
      atom_str += v[i].chain;
      atom_str += "\",";
      atom_str += int_to_string(v[i].resno);
      atom_str += ",\"";
      atom_str += v[i].insertion_code;
      atom_str += "\",\"";
      atom_str += v[i].atom_name;
      atom_str += "\",\"";
      atom_str += v[i].alt_conf;
      atom_str += " \"";

      std::string button_label("Clash gap: ");
      button_label += float_to_string(v[i].float_user_data);
      button_label += " : ";
      button_label += v[i].chain;
      button_label += " ";
      button_label += int_to_string(v[i].resno);
      button_label += " ";
      if (v[i].insertion_code != "") {
         button_label += v[i].insertion_code;
         button_label += " ";
      }
      button_label += v[i].atom_name;
      if (v[i].alt_conf != "") {
         button_label += ",";
         button_label += v[i].alt_conf;
         button_label += " ";
      }

      std::string s = "[";
      s += single_quote(button_label);
      s += ",";
      s += int_to_string(v[i].int_user_data);
      s += ",";
      s += atom_str;
      s += "],\n";

      r += s;
   }

   r += "]";
   return r;
#else
   return "";
#endif // PYTHON
}

// error_type is e.g. "Z score", "Clash gap"
// 
std::string
coot::util::interesting_things_list_with_fix(const std::vector<coot::util::atom_spec_and_button_info_t> &v,
					     const std::string error_type) {

#ifdef USE_GUILE
   // e.g. (list) for empty v
   // (list (list "button label" imol-no chain-id resno atom-name)
   //       (list "button label" imol-no chain-id resno atom-name)
   // )
   //
   // if we have a fix, the callback function is not "" and then the
   // returned thing becomes:
   // 
   // (list (list "button label" imol-no chain-id resno atom-name callback-func)
   //       (list "button label" imol-no chain-id resno atom-name callback-func)
   // )
   //
   // where callback-func is e.g. (lambda() (do-180-degree-side-chain-flip 0 "A" 45 "" ""))

   std::string r = " (list ";

   for (unsigned int i=0; i<v.size(); i++) {

      std::string atom_str("\"");
      atom_str += v[i].as.chain;
      atom_str += "\" ";
      atom_str += int_to_string(v[i].as.resno);
      atom_str += " \"";
      atom_str += v[i].as.insertion_code;
      atom_str += "\" \"";
      atom_str += v[i].as.atom_name;
      atom_str += "\" \"";
      atom_str += v[i].as.alt_conf;
      atom_str += " \"";

      std::string button_label = v[i].button_label;

      std::string s = "(list ";
      s += single_quote(button_label);
      s += " ";
      s += int_to_string(v[i].as.int_user_data);
      s += " ";
      s += atom_str;

      if (v[i].callback_func != "") {
	 s += " ";
	 s +=  v[i].callback_func;
      }
      
      s += ")\n";
      
      r += s;
   }

   r += ")";
   return r;
#else
#ifdef USE_PYTHON
// BL says:: here again we need a [] list in python 
   std::string r = "[";

   for (unsigned int i=0; i<v.size(); i++) {

      std::string atom_str("\"");
      atom_str += v[i].as.chain;
      atom_str += "\",";
      atom_str += int_to_string(v[i].as.resno);
      atom_str += ",\"";
      atom_str += v[i].as.insertion_code;
      atom_str += "\",\"";
      atom_str += v[i].as.atom_name;
      atom_str += "\",\"";
      atom_str += v[i].as.alt_conf;
      atom_str += " \"";

      std::string button_label = v[i].button_label;

      std::string s = "[";
      s += single_quote(button_label);
      s += ",";
      s += int_to_string(v[i].as.int_user_data);
      s += ",";
      s += atom_str;

      if (v[i].callback_func != "") {
         s += ",";
         s +=  v[i].callback_func;
      }

      if (i<(v.size()-1)) {
         s += "],\n";
      } else {
         s += "]\n";
      }

      r += s;
   }

   r += "]";
   return r;
#else
   return "";
#endif // PYTHON
#endif // GUILE
}

// python version
std::string
coot::util::interesting_things_list_with_fix_py(const std::vector<coot::util::atom_spec_and_button_info_t> &v,
						const std::string error_type) {
#ifdef USE_PYTHON
// BL says:: here again we need a [] list in python 
   std::string r = "[";

   for (unsigned int i=0; i<v.size(); i++) {

      std::string atom_str("\"");
      atom_str += v[i].as.chain;
      atom_str += "\",";
      atom_str += int_to_string(v[i].as.resno);
      atom_str += ",\"";
      atom_str += v[i].as.insertion_code;
      atom_str += "\",\"";
      atom_str += v[i].as.atom_name;
      atom_str += "\",\"";
      atom_str += v[i].as.alt_conf;
      atom_str += " \"";

      std::string button_label = v[i].button_label;

      std::string s = "[";
      s += single_quote(button_label);
      s += ",";
      s += int_to_string(v[i].as.int_user_data);
      s += ",";
      s += atom_str;

      if (v[i].callback_func != "") {
         s += ",";
         s +=  v[i].callback_func;
      }

      if (i<(v.size()-1)) {
         s += "],\n";
      } else {
         s += "]\n";
      }

      r += s;
   }

   r += "]";
   return r;
#else
   return "";
#endif // PYTHON
}


// Return the RTop that matches moving to reference.  Don't move
// moving though.
//
// reference is residue
// moving is std_base
// 
std::pair<bool, clipper::RTop_orth> coot::util::base_to_base(CResidue *residue,
							     CResidue *std_base) {

   bool good_rtop_flag = 0;

   clipper::Mat33<double> m_dum(1,0,0,0,1,0,0,0,1);
   clipper::Coord_orth pt_dum(0,0,0);
   clipper::RTop_orth rtop(m_dum, pt_dum);
   
   std::vector<std::string> adenine;  // Pirimidine
   adenine.push_back(" N9 ");
   adenine.push_back(" C8 ");
   adenine.push_back(" N7 ");
   adenine.push_back(" C5 ");
   adenine.push_back(" C4 ");
   // 
   adenine.push_back(" N1 ");
   adenine.push_back(" C2 ");
   adenine.push_back(" N3 ");
   adenine.push_back(" C6 ");
   adenine.push_back(" N6 ");

   std::vector<std::string> guanine; // Pirimidine
   guanine.push_back(" N9 ");
   guanine.push_back(" C8 ");
   guanine.push_back(" N7 ");
   guanine.push_back(" C5 ");
   guanine.push_back(" C4 ");
   //
   guanine.push_back(" N1 ");
   guanine.push_back(" C2 ");
   guanine.push_back(" N3 ");
   guanine.push_back(" C6 ");
   guanine.push_back(" O6 ");
   guanine.push_back(" N2 "); // No matcher for this in adenine

   std::vector<std::string> thymine;  // Purine
   thymine.push_back(" N1 ");
   thymine.push_back(" C2 ");
   thymine.push_back(" N3 ");
   thymine.push_back(" C4 ");
   thymine.push_back(" C5 ");
   thymine.push_back(" C6 ");
   // 
   thymine.push_back(" O2 ");
   thymine.push_back(" O4 ");
   thymine.push_back(" C5M");
   
   std::vector<std::string> cytosine;  // Purine
   cytosine.push_back(" N1 ");
   cytosine.push_back(" C2 ");
   cytosine.push_back(" N3 ");
   cytosine.push_back(" C4 ");
   cytosine.push_back(" C5 ");
   cytosine.push_back(" C6 ");
   // 
   cytosine.push_back(" O2 ");
   cytosine.push_back(" N4 ");
   
   std::vector<std::string> uracil;  // Purine
   uracil.push_back(" N1 ");
   uracil.push_back(" C2 ");
   uracil.push_back(" N3 ");
   uracil.push_back(" C4 ");
   uracil.push_back(" C5 ");
   uracil.push_back(" C6 ");
   // 
   uracil.push_back(" O2 ");
   uracil.push_back(" O4 ");
   

   // These next 2 are in match order, don't change it.
   std::vector<std::string> purine; // A and G
   purine.push_back(" N9 ");
   purine.push_back(" C4 ");
   purine.push_back(" C5 ");
   purine.push_back(" N7 ");
   purine.push_back(" C8 ");

   std::vector<std::string> pyrimidine; // T, U and C
   pyrimidine.push_back(" N1 ");
   pyrimidine.push_back(" C2 ");
   pyrimidine.push_back(" N3 ");
   pyrimidine.push_back(" C5 ");
   pyrimidine.push_back(" C6 ");
   pyrimidine.push_back(" C4 ");

   // We need to know whether we have purine or pyrimidine for both
   // the molecule base and the std_base.
   //
   // We need to get all (5 for pyrimidine, 6 for purine) the
   // coordinates for both bases.
   // 
   // If they are either or both are pyrimidine we match 5 atoms,
   // If they are both purine we match 6 atoms.


   // So what are the input base types?
   //
   // These for flags should be set to something after our test
   short int mol_base_is_pyrimidine = -1;
   short int mol_base_is_purine     = -1;
   short int std_base_is_pyrimidine = -1;
   short int std_base_is_purine     = -1;

   std::string mol_base_name = residue->GetResName();
   std::string std_base_name = std_base->GetResName();

   if (mol_base_name == "Ar" || mol_base_name == "Ad" ||
       mol_base_name == "Gr" || mol_base_name == "Gd") {
      mol_base_is_purine = 1;
      mol_base_is_pyrimidine = 0;
   }

   if (mol_base_name == "Cr" || mol_base_name == "Cd" ||
       mol_base_name == "Ur" || mol_base_name == "Ud" ||
       mol_base_name == "Tr" || mol_base_name == "Td") {
      mol_base_is_pyrimidine = 1;
      mol_base_is_purine = 0;
   }

   if (std_base_name == "Ar" || std_base_name == "Ad" ||
       std_base_name == "Gr" || std_base_name == "Gd") {
      std_base_is_purine = 1;
      std_base_is_pyrimidine = 0;
   }

   if (std_base_name == "Cr" || std_base_name == "Cd" ||
       std_base_name == "Tr" || std_base_name == "Td" ||
       std_base_name == "Ur" || std_base_name == "Ud") {
      std_base_is_pyrimidine = 1;
      std_base_is_purine = 0;
   }

   if ((mol_base_is_pyrimidine == -1) || (mol_base_is_purine == -1) || 
       (std_base_is_pyrimidine == -1) || (std_base_is_purine == -1) ) {

      std::cout << "ERROR:: unassigned type "
		<< "mol_base_is_pyrimidine:" << " "
		<< mol_base_is_pyrimidine << " "
		<< "mol_base_is_purine: " << " "
		<< mol_base_is_purine << " "
		<< "std_base_is_pyrimidine: " << " "
		<< std_base_is_pyrimidine << " "
		<< "std_base_is_purine: " << " "
		<< std_base_is_purine << " "
		<< residue->GetResName() << " " << std_base->GetResName()
		<< std::endl;

      
   } else {
   
//       std::cout << "DEBUG:: assigned types "
// 		<< "mol_base_is_pyrimidine:" << " "
// 		<< mol_base_is_pyrimidine << " "
// 		<< "mol_base_is_purine: " << " "
// 		<< mol_base_is_purine << " "
// 		<< "std_base_is_pyrimidine: " << " "
// 		<< std_base_is_pyrimidine << " "
// 		<< "std_base_is_purine: " << " "
// 		<< std_base_is_purine << " "
// 		<< residue->GetResName() << " " << std_base->GetResName()
// 		<< std::endl;

      int n_match_atoms = 5;
      if (mol_base_is_pyrimidine && std_base_is_pyrimidine)
	 n_match_atoms = 6;

      std::vector<std::string> moving_name_vector;
      std::vector<std::string> refrce_name_vector;

      if (std_base_is_purine)
	 moving_name_vector = purine;
      else
	 moving_name_vector = pyrimidine;

      if (mol_base_is_purine)
	 refrce_name_vector = purine;
      else
	 refrce_name_vector = pyrimidine;
      
      PCAtom *std_base_atoms;
      int n_std_base_atoms;

      PCAtom *mol_base_atoms;
      int n_mol_base_atoms;
      
      residue->GetAtomTable( mol_base_atoms, n_mol_base_atoms);
      std_base->GetAtomTable(std_base_atoms, n_std_base_atoms);

      std::vector<clipper::Coord_orth> refrce_atom_positions;
      std::vector<clipper::Coord_orth> moving_atom_positions;

//       for (unsigned int i=0; i<refrce_name_vector.size(); i++)
// 	 std::cout << "ref base search atom :" << refrce_name_vector[i]
// 		   << ":" << std::endl;
//       for (unsigned int i=0; i<moving_name_vector.size(); i++)
// 	 std::cout << "mov base search atom :" << moving_name_vector[i]
// 		   << ":" << std::endl;
      
      for (int j=0; j<n_match_atoms; j++) {
	 for (int i=0; i<n_mol_base_atoms; i++) {
	    std::string atom_name = mol_base_atoms[i]->name;
	    if (refrce_name_vector[j] == atom_name) {
	       refrce_atom_positions.push_back(clipper::Coord_orth(mol_base_atoms[i]->x,
								   mol_base_atoms[i]->y,
								   mol_base_atoms[i]->z));
	       // std::cout << "Found " << atom_name << " in reference " << std::endl;
	    }
	 }
      }

      for (int j=0; j<n_match_atoms; j++) {
	 for (int i=0; i<n_std_base_atoms; i++) {
	 std::string atom_name = std_base_atoms[i]->name;
	    if (moving_name_vector[j] == atom_name) {
	       moving_atom_positions.push_back(clipper::Coord_orth(std_base_atoms[i]->x,
								   std_base_atoms[i]->y,
								   std_base_atoms[i]->z));
	       // std::cout << "Found " << atom_name << " in moving (std) base " << std::endl;
	    }
	 }
      }

      if (int(refrce_atom_positions.size()) != n_match_atoms) {
	 std::cout << "ERROR:: wrong number of reference atoms found! "
		   << refrce_atom_positions.size() << std::endl;
      } else {

	 if (int(moving_atom_positions.size()) != n_match_atoms) {
	    std::cout << "ERROR:: wrong number of moving atoms found! "
		   << moving_atom_positions.size() << std::endl;

	 } else {

	    rtop = clipper::RTop_orth (moving_atom_positions, refrce_atom_positions);
	    good_rtop_flag = 1;
	 }
      }
   }
   return std::pair<bool, clipper::RTop_orth> (good_rtop_flag, rtop);
}

// Return the RTop that matches moving to reference.  Include the base
// atoms, *AND* the furanose and phosphate atoms.
//
// reference is residue
// moving is std_base
// 
std::pair<bool, clipper::RTop_orth> coot::util::nucleotide_to_nucleotide(CResidue *residue,
									 CResidue *std_base) {

   bool good_rtop_flag = 0;
   clipper::Mat33<double> m_dum(1,0,0,0,1,0,0,0,1);
   clipper::Coord_orth pt_dum(0,0,0);
   clipper::RTop_orth rtop(m_dum, pt_dum);
   
   std::vector<std::string> adenine;  // Pirimidine
   adenine.push_back(" N9 ");
   adenine.push_back(" C8 ");
   adenine.push_back(" N7 ");
   adenine.push_back(" C5 ");
   adenine.push_back(" C4 ");
   // 
   adenine.push_back(" N1 ");
   adenine.push_back(" C2 ");
   adenine.push_back(" N3 ");
   adenine.push_back(" C6 ");
   adenine.push_back(" N6 ");

   std::vector<std::string> guanine; // Pirimidine
   guanine.push_back(" N9 ");
   guanine.push_back(" C8 ");
   guanine.push_back(" N7 ");
   guanine.push_back(" C5 ");
   guanine.push_back(" C4 ");
   //
   guanine.push_back(" N1 ");
   guanine.push_back(" C2 ");
   guanine.push_back(" N3 ");
   guanine.push_back(" C6 ");
   guanine.push_back(" O6 ");
   guanine.push_back(" N2 "); // No matcher for this in adenine

   std::vector<std::string> thymine;  // Purine
   thymine.push_back(" N1 ");
   thymine.push_back(" C2 ");
   thymine.push_back(" N3 ");
   thymine.push_back(" C4 ");
   thymine.push_back(" C5 ");
   thymine.push_back(" C6 ");
   // 
   thymine.push_back(" O2 ");
   thymine.push_back(" O4 ");
   thymine.push_back(" C5M");
   
   std::vector<std::string> cytosine;  // Purine
   cytosine.push_back(" N1 ");
   cytosine.push_back(" C2 ");
   cytosine.push_back(" N3 ");
   cytosine.push_back(" C4 ");
   cytosine.push_back(" C5 ");
   cytosine.push_back(" C6 ");
   // 
   cytosine.push_back(" O2 ");
   cytosine.push_back(" N4 ");
   
   std::vector<std::string> uracil;  // Purine
   uracil.push_back(" N1 ");
   uracil.push_back(" C2 ");
   uracil.push_back(" N3 ");
   uracil.push_back(" C4 ");
   uracil.push_back(" C5 ");
   uracil.push_back(" C6 ");
   // 
   uracil.push_back(" O2 ");
   uracil.push_back(" O4 ");
   

   // These next 2 are in match order, don't change it.
   std::vector<std::string> purine; // A and G
   purine.push_back(" N9 ");
   purine.push_back(" C4 ");
   purine.push_back(" C5 ");
   purine.push_back(" N7 ");
   purine.push_back(" C8 ");

   std::vector<std::string> pyrimidine; // T, U and C
   pyrimidine.push_back(" N1 ");
   pyrimidine.push_back(" C2 ");
   pyrimidine.push_back(" N3 ");
   pyrimidine.push_back(" C5 ");
   pyrimidine.push_back(" C6 ");
   pyrimidine.push_back(" C4 ");

   // We need to know whether we have purine or pyrimidine for both
   // the molecule base and the std_base.
   //
   // We need to get all (5 for pyrimidine, 6 for purine) the
   // coordinates for both bases.
   // 
   // If they are either or both are pyrimidine we match 5 atoms,
   // If they are both purine we match 6 atoms.


   // So what are the input base types?
   //
   // These for flags should be set to something after our test
   short int mol_base_is_pyrimidine = -1;
   short int mol_base_is_purine     = -1;
   short int std_base_is_pyrimidine = -1;
   short int std_base_is_purine     = -1;

   std::string mol_base_name = residue->GetResName();
   std::string std_base_name = std_base->GetResName();

   if (mol_base_name == "Ar" || mol_base_name == "Ad" ||
       mol_base_name == "Gr" || mol_base_name == "Gd") {
      mol_base_is_purine = 1;
      mol_base_is_pyrimidine = 0;
   }

   if (mol_base_name == "Cr" || mol_base_name == "Cd" ||
       mol_base_name == "Ur" || mol_base_name == "Ud" ||
       mol_base_name == "Tr" || mol_base_name == "Td") {
      mol_base_is_pyrimidine = 1;
      mol_base_is_purine = 0;
   }

   if (std_base_name == "Ar" || std_base_name == "Ad" ||
       std_base_name == "Gr" || std_base_name == "Gd") {
      std_base_is_purine = 1;
      std_base_is_pyrimidine = 0;
   }

   if (std_base_name == "Cr" || std_base_name == "Cd" ||
       std_base_name == "Tr" || std_base_name == "Td" ||
       std_base_name == "Ur" || std_base_name == "Ud") {
      std_base_is_pyrimidine = 1;
      std_base_is_purine = 0;
   }

   if ((mol_base_is_pyrimidine == -1) || (mol_base_is_purine == -1) || 
       (std_base_is_pyrimidine == -1) || (std_base_is_purine == -1) ) {

      std::cout << "ERROR:: unassigned type "
		<< "mol_base_is_pyrimidine:" << " "
		<< mol_base_is_pyrimidine << " "
		<< "mol_base_is_purine: " << " "
		<< mol_base_is_purine << " "
		<< "std_base_is_pyrimidine: " << " "
		<< std_base_is_pyrimidine << " "
		<< "std_base_is_purine: " << " "
		<< std_base_is_purine << " "
		<< residue->GetResName() << " " << std_base->GetResName()
		<< std::endl;

      
   } else {
   
//       std::cout << "DEBUG:: assigned types "
// 		<< "mol_base_is_pyrimidine:" << " "
// 		<< mol_base_is_pyrimidine << " "
// 		<< "mol_base_is_purine: " << " "
// 		<< mol_base_is_purine << " "
// 		<< "std_base_is_pyrimidine: " << " "
// 		<< std_base_is_pyrimidine << " "
// 		<< "std_base_is_purine: " << " "
// 		<< std_base_is_purine << " "
// 		<< residue->GetResName() << " " << std_base->GetResName()
// 		<< std::endl;

      int n_match_atoms = 5;
      if (mol_base_is_pyrimidine && std_base_is_pyrimidine)
	 n_match_atoms = 6;

      std::vector<std::string> moving_name_vector;
      std::vector<std::string> refrce_name_vector;

      if (std_base_is_purine)
	 moving_name_vector = purine;
      else
	 moving_name_vector = pyrimidine;

      if (mol_base_is_purine)
	 refrce_name_vector = purine;
      else
	 refrce_name_vector = pyrimidine;
      
      PCAtom *std_base_atoms;
      int n_std_base_atoms;

      PCAtom *mol_base_atoms;
      int n_mol_base_atoms;
      
      residue->GetAtomTable( mol_base_atoms, n_mol_base_atoms);
      std_base->GetAtomTable(std_base_atoms, n_std_base_atoms);

      std::vector<clipper::Coord_orth> refrce_atom_positions;
      std::vector<clipper::Coord_orth> moving_atom_positions;

//       for (unsigned int i=0; i<refrce_name_vector.size(); i++)
// 	 std::cout << "ref base search atom :" << refrce_name_vector[i]
// 		   << ":" << std::endl;
//       for (unsigned int i=0; i<moving_name_vector.size(); i++)
// 	 std::cout << "mov base search atom :" << moving_name_vector[i]
// 		   << ":" << std::endl;
      
      for (int j=0; j<n_match_atoms; j++) {
	 for (int i=0; i<n_mol_base_atoms; i++) {
	    std::string atom_name = mol_base_atoms[i]->name;
	    if (refrce_name_vector[j] == atom_name) {
	       refrce_atom_positions.push_back(clipper::Coord_orth(mol_base_atoms[i]->x,
								   mol_base_atoms[i]->y,
								   mol_base_atoms[i]->z));
	       // std::cout << "Found " << atom_name << " in reference " << std::endl;
	    }
	 }
      }

      for (int j=0; j<n_match_atoms; j++) {
	 for (int i=0; i<n_std_base_atoms; i++) {
	 std::string atom_name = std_base_atoms[i]->name;
	    if (moving_name_vector[j] == atom_name) {
	       moving_atom_positions.push_back(clipper::Coord_orth(std_base_atoms[i]->x,
								   std_base_atoms[i]->y,
								   std_base_atoms[i]->z));
	       // std::cout << "Found " << atom_name << " in moving (std) base " << std::endl;
	    }
	 }
      }

      if (int(refrce_atom_positions.size()) != n_match_atoms) {
	 std::cout << "ERROR:: wrong number of reference atoms found! "
		   << refrce_atom_positions.size() << std::endl;
      } else {

	 if (int(moving_atom_positions.size()) != n_match_atoms) {
	    std::cout << "ERROR:: wrong number of moving atoms found! "
		   << moving_atom_positions.size() << std::endl;

	 } else {

	    // all nucleodites have these atoms, use them to do a match:
	    // 
	    std::vector<std::string> const_nuc_atoms;
	    const_nuc_atoms.push_back(" C1*");
	    const_nuc_atoms.push_back(" C2*");
	    const_nuc_atoms.push_back(" C3*");
	    const_nuc_atoms.push_back(" C4*");
	    const_nuc_atoms.push_back(" C5*");
	    const_nuc_atoms.push_back(" O3*");
	    const_nuc_atoms.push_back(" O4*");
	    const_nuc_atoms.push_back(" O5*");
	    const_nuc_atoms.push_back(" P  ");

	    // We want to match the bases too, don't we?
// 	    moving_atom_positions.clear();
// 	    refrce_atom_positions.clear(); 

	    for (unsigned int inuc=0; inuc<const_nuc_atoms.size(); inuc++) {
	       for (int istd=0; istd<n_std_base_atoms; istd++) {
		  std::string std_base_atom_name = std_base_atoms[istd]->name;
		  if (std_base_atom_name == const_nuc_atoms[inuc]) { 
		     for (int imol=0; imol<n_mol_base_atoms; imol++) {
			std::string mol_base_atom_name = mol_base_atoms[imol]->name;
			if (mol_base_atom_name == std_base_atom_name) {
			   std::string altconf1 = std_base_atoms[istd]->altLoc;
			   std::string altconf2 = mol_base_atoms[imol]->altLoc;
			   if (altconf1 == altconf2) {
			      clipper::Coord_orth s(std_base_atoms[istd]->x,
						    std_base_atoms[istd]->y,
						    std_base_atoms[istd]->z);
			      clipper::Coord_orth m(mol_base_atoms[imol]->x,
						    mol_base_atoms[imol]->y,
						    mol_base_atoms[imol]->z);
			   
// 			      std::cout << "---" << std::endl;
// 			      std::cout << std_base_atoms[istd]->GetSeqNum() << " "
// 					<< std_base_atoms[istd]->name << " ("
// 					<< std_base_atoms[istd]->x << ","
// 					<< std_base_atoms[istd]->y << ","
// 					<< std_base_atoms[istd]->z << ")" << std::endl;

// 			      std::cout << mol_base_atoms[imol]->GetSeqNum() << " "
// 					<< mol_base_atoms[imol]->name << " ("
// 					<< mol_base_atoms[imol]->x << ","
// 					<< mol_base_atoms[imol]->y << ","
// 					<< mol_base_atoms[imol]->z << ")" << std::endl;

			      refrce_atom_positions.push_back(m);
			      moving_atom_positions.push_back(s);
			   }
			}
		     }
		  }
	       }
	    }

// 	    std::cout << "debug:: matching "
// 		      << moving_atom_positions.size() << " atoms" << std::endl;
	    
	    rtop = clipper::RTop_orth (moving_atom_positions, refrce_atom_positions);
	    good_rtop_flag = 1;
	 }
      }
   }
   return std::pair<bool, clipper::RTop_orth> (good_rtop_flag, rtop);
}

void
coot::util::mutate_internal(CResidue *residue, CResidue *std_residue,
			    short int is_from_shelx_ins_flag) {

   PPCAtom residue_atoms;
   int nResidueAtoms;
   residue->GetAtomTable(residue_atoms, nResidueAtoms);
   std::string res_name(residue->GetResName());

   PPCAtom std_residue_atoms;
   int n_std_ResidueAtoms;
   std_residue->GetAtomTable(std_residue_atoms, n_std_ResidueAtoms);

   //
   std::string old_seg_id_for_residue_atoms;
   bool use_old_seg_id = 0;
   try {
      old_seg_id_for_residue_atoms = coot::residue_atoms_segid(residue);
      use_old_seg_id = 1;
   }
   catch (std::runtime_error mess) {
   } 

   bool verb = 0;
   if (verb) { 
      std::cout << "Mutate Atom Tables" << std::endl;
      std::cout << "Before" << std::endl;
      for(int i=0; i<nResidueAtoms; i++)
	 std::cout << residue_atoms[i]->name << std::endl;
      std::cout << "To be replaced by:" << std::endl;
      for(int i=0; i<n_std_ResidueAtoms; i++)
	 std::cout << std_residue_atoms[i]->name << std::endl;
   }

   // std::vector<CAtom *> keep_atoms;
   for(int i=0; i<nResidueAtoms; i++) {
      std::string residue_this_atom (residue_atoms[i]->name);
      if (coot::is_main_chain_p(residue_atoms[i])) {
	 // copy the atom, and add the copied atom to the keep_atoms
	 // vector
// 	 CAtom *c_copy = new CAtom;
// 	 c_copy->Copy(residue_atoms[i]);
// 	 keep_atoms.push_back(c_copy);
      } else { 
	 residue->DeleteAtom(i);
      }
   };

//    // now add the mainchain atoms back
//    for (unsigned int i=0; i<keep_atoms.size(); i++)
//       residue->AddAtom(keep_atoms[i]);
   
   // add all atoms of std residue to 
   for(int i=0; i<n_std_ResidueAtoms; i++) {
      std::string std_residue_this_atom (std_residue_atoms[i]->name);
      if (! coot::is_main_chain_p(std_residue_atoms[i])) {
	 if (is_from_shelx_ins_flag)
	    std_residue_atoms[i]->occupancy = 11.0;
	 CAtom *copy_at = new CAtom;
	 copy_at->Copy(std_residue_atoms[i]);
	 residue->AddAtom(copy_at);
	 if (use_old_seg_id) {
	    strcpy(copy_at->segID, old_seg_id_for_residue_atoms.c_str());
	 }
      }
   };

   residue->SetResName(std_residue->GetResName());
   residue->TrimAtomTable();
}

// return state, 0 bad, 1 good
// 
int
coot::util::mutate(CResidue *res, CResidue *std_res_unoriented, short int shelx_flag) {

   int istate = 0; 
   std::pair<clipper::RTop_orth, short int> rtop_pair =
      coot::util::get_ori_to_this_res(res);

   PPCAtom residue_atoms;
   int nResidueAtoms;
   std_res_unoriented->GetAtomTable(residue_atoms, nResidueAtoms);
   if (nResidueAtoms == 0) {
      std::cout << " something broken in atom residue selection in ";
      std::cout << "mutate, got 0 atoms" << std::endl;
   } else {
      for(int iat=0; iat<nResidueAtoms; iat++) {
	 clipper::Coord_orth co(residue_atoms[iat]->x,
				residue_atoms[iat]->y,
				residue_atoms[iat]->z);
	 clipper::Coord_orth rotted = co.transform(rtop_pair.first);
	 residue_atoms[iat]->x = rotted.x();
	 residue_atoms[iat]->y = rotted.y();
	 residue_atoms[iat]->z = rotted.z();
      }
      coot::util::mutate_internal(res, std_res_unoriented, shelx_flag); // it's oriented now.
      istate = 1;
   }
   return istate;
}



// Here std_base is at some arbitary position when passed.
// 
void
coot::util::mutate_base(CResidue *residue, CResidue *std_base) {

   std::vector<std::string> adenine;  // Pyrimidine
   adenine.push_back(" N9 ");
   adenine.push_back(" C8 ");
   adenine.push_back(" N7 ");
   adenine.push_back(" C5 ");
   adenine.push_back(" C4 ");
   // 
   adenine.push_back(" N1 ");
   adenine.push_back(" C2 ");
   adenine.push_back(" N3 ");
   adenine.push_back(" C6 ");
   adenine.push_back(" N6 ");

   std::vector<std::string> guanine; // Pyrimidine
   guanine.push_back(" N9 ");
   guanine.push_back(" C8 ");
   guanine.push_back(" N7 ");
   guanine.push_back(" C5 ");
   guanine.push_back(" C4 ");
   //
   guanine.push_back(" N1 ");
   guanine.push_back(" C2 ");
   guanine.push_back(" N3 ");
   guanine.push_back(" C6 ");
   guanine.push_back(" O6 ");
   guanine.push_back(" N2 "); // No matcher for this in adenine

   std::vector<std::string> thymine;  // Purine
   thymine.push_back(" N1 ");
   thymine.push_back(" C2 ");
   thymine.push_back(" N3 ");
   thymine.push_back(" C4 ");
   thymine.push_back(" C5 ");
   thymine.push_back(" C6 ");
   // 
   thymine.push_back(" O2 ");
   thymine.push_back(" O4 ");
   thymine.push_back(" C5M");
   
   std::vector<std::string> cytosine;  // Purine
   cytosine.push_back(" N1 ");
   cytosine.push_back(" C2 ");
   cytosine.push_back(" N3 ");
   cytosine.push_back(" C4 ");
   cytosine.push_back(" C5 ");
   cytosine.push_back(" C6 ");
   // 
   cytosine.push_back(" O2 ");
   cytosine.push_back(" N4 ");
   
   std::vector<std::string> uracil;  // Purine
   uracil.push_back(" N1 ");
   uracil.push_back(" C2 ");
   uracil.push_back(" N3 ");
   uracil.push_back(" C4 ");
   uracil.push_back(" C5 ");
   uracil.push_back(" C6 ");
   // 
   uracil.push_back(" O2 ");
   uracil.push_back(" O4 ");
   

   // These next 2 are in match order, don't change it.
   std::vector<std::string> purine; // A and G
   purine.push_back(" N9 ");
   purine.push_back(" C4 ");
   purine.push_back(" C5 ");
   purine.push_back(" N7 ");
   purine.push_back(" C8 ");

   std::vector<std::string> pyrimidine; // T, U and C
   pyrimidine.push_back(" N1 ");
   pyrimidine.push_back(" C2 ");
   pyrimidine.push_back(" N3 ");
   pyrimidine.push_back(" C5 ");
   pyrimidine.push_back(" C6 ");
   pyrimidine.push_back(" C4 ");


   std::string old_seg_id_for_residue_atoms;
   bool use_old_seg_id = 0;
   try {
      old_seg_id_for_residue_atoms = coot::residue_atoms_segid(residue);
      use_old_seg_id = 1;
   }
   catch (std::runtime_error mess) {
   } 


   // We need to know whether we have purine or pyrimidine for both
   // the molecule base and the std_base.
   //
   // We need to get all (5 for pyrimidine, 6 for purine) the
   // coordinates for both bases.
   // 
   // If they are either or both are pyrimidine we match 5 atoms,
   // If they are both purine we match 6 atoms.


   // So what are the input base types?
   //
   // These for flags should be set to something after our test
   short int mol_base_is_pyrimidine = -1;
   short int mol_base_is_purine     = -1;
   short int std_base_is_pyrimidine = -1;
   short int std_base_is_purine     = -1;

   std::string mol_base_name = residue->GetResName();
   std::string std_base_name = std_base->GetResName();

   if (mol_base_name == "Ar" || mol_base_name == "Ad" ||
       mol_base_name == "Gr" || mol_base_name == "Gd") {
      mol_base_is_purine = 1;
      mol_base_is_pyrimidine = 0;
   }

   if (mol_base_name == "Cr" || mol_base_name == "Cd" ||
       mol_base_name == "Ur" || mol_base_name == "Ud" ||
       mol_base_name == "Tr" || mol_base_name == "Td") {
      mol_base_is_pyrimidine = 1;
      mol_base_is_purine = 0;
   }

   if (std_base_name == "Ar" || std_base_name == "Ad" ||
       std_base_name == "Gr" || std_base_name == "Gd") {
      std_base_is_purine = 1;
      std_base_is_pyrimidine = 0;
   }

   if (std_base_name == "Cr" || std_base_name == "Cd" ||
       std_base_name == "Tr" || std_base_name == "Td" ||
       std_base_name == "Ur" || std_base_name == "Ud") {
      std_base_is_pyrimidine = 1;
      std_base_is_purine = 0;
   }

   if ((mol_base_is_pyrimidine == -1) || (mol_base_is_purine == -1) || 
       (std_base_is_pyrimidine == -1) || (std_base_is_purine == -1) ) {

      std::cout << "ERROR:: unassigned type "
		<< "mol_base_is_pyrimidine:" << " "
		<< mol_base_is_pyrimidine << " "
		<< "mol_base_is_purine: " << " "
		<< mol_base_is_purine << " "
		<< "std_base_is_pyrimidine: " << " "
		<< std_base_is_pyrimidine << " "
		<< "std_base_is_purine: " << " "
		<< std_base_is_purine << " "
		<< residue->GetResName() << " " << std_base->GetResName()
		<< std::endl;

      
   } else {
   
//       std::cout << "DEBUG:: assigned types "
// 		<< "mol_base_is_pyrimidine:" << " "
// 		<< mol_base_is_pyrimidine << " "
// 		<< "mol_base_is_purine: " << " "
// 		<< mol_base_is_purine << " "
// 		<< "std_base_is_pyrimidine: " << " "
// 		<< std_base_is_pyrimidine << " "
// 		<< "std_base_is_purine: " << " "
// 		<< std_base_is_purine << " "
// 		<< residue->GetResName() << " " << std_base->GetResName()
// 		<< std::endl;

      int n_match_atoms = 5;
      if (mol_base_is_pyrimidine && std_base_is_pyrimidine)
	 n_match_atoms = 6;

      std::vector<std::string> moving_name_vector;
      std::vector<std::string> refrce_name_vector;

      if (std_base_is_purine)
	 moving_name_vector = purine;
      else
	 moving_name_vector = pyrimidine;

      if (mol_base_is_purine)
	 refrce_name_vector = purine;
      else
	 refrce_name_vector = pyrimidine;
      
      PCAtom *std_base_atoms;
      int n_std_base_atoms;

      PCAtom *mol_base_atoms;
      int n_mol_base_atoms;
      
      residue->GetAtomTable( mol_base_atoms, n_mol_base_atoms);
      std_base->GetAtomTable(std_base_atoms, n_std_base_atoms);

      std::vector<clipper::Coord_orth> refrce_atom_positions;
      std::vector<clipper::Coord_orth> moving_atom_positions;

//       for (unsigned int i=0; i<refrce_name_vector.size(); i++)
// 	 std::cout << "ref base search atom :" << refrce_name_vector[i]
// 		   << ":" << std::endl;
//       for (unsigned int i=0; i<moving_name_vector.size(); i++)
// 	 std::cout << "mov base search atom :" << moving_name_vector[i]
// 		   << ":" << std::endl;
      
      for (int j=0; j<n_match_atoms; j++) {
	 for (int i=0; i<n_mol_base_atoms; i++) {
	    std::string atom_name = mol_base_atoms[i]->name;
	    if (refrce_name_vector[j] == atom_name) {
	       refrce_atom_positions.push_back(clipper::Coord_orth(mol_base_atoms[i]->x,
								   mol_base_atoms[i]->y,
								   mol_base_atoms[i]->z));
	       // std::cout << "Found " << atom_name << " in reference " << std::endl;
	    }
	 }
      }

      for (int j=0; j<n_match_atoms; j++) {
	 for (int i=0; i<n_std_base_atoms; i++) {
	 std::string atom_name = std_base_atoms[i]->name;
	    if (moving_name_vector[j] == atom_name) {
	       moving_atom_positions.push_back(clipper::Coord_orth(std_base_atoms[i]->x,
								   std_base_atoms[i]->y,
								   std_base_atoms[i]->z));
	       // std::cout << "Found " << atom_name << " in moving (std) base " << std::endl;
	    }
	 }
      }

      if (int(refrce_atom_positions.size()) != n_match_atoms) {
	 std::cout << "ERROR:: wrong number of reference atoms found! "
		   << refrce_atom_positions.size() << std::endl;
      } else {

	 if (int(moving_atom_positions.size()) != n_match_atoms) {
	    std::cout << "ERROR:: wrong number of moving atoms found! "
		   << moving_atom_positions.size() << std::endl;

	 } else {

	    clipper::RTop_orth rtop (moving_atom_positions, refrce_atom_positions);
      
	    double sum_dist = 0.0;
	    double sum_dist2 = 0.0;
	    double mind =  999999999.9;
	    double maxd = -999999999.9;
	    double d;
	    for (unsigned int i=0; i<refrce_atom_positions.size(); i++) {
	       d = clipper::Coord_orth::length(refrce_atom_positions[i],
					       clipper::Coord_orth(moving_atom_positions[i].transform(rtop)));
	       sum_dist  += d;
	       sum_dist2 += d*d;
	       if (d>maxd)
		  maxd = d;
	       if (d<mind)
		  mind = d;
	    }
	    double mean = sum_dist/double(moving_atom_positions.size());
	    double var  = sum_dist2/double(moving_atom_positions.size()); // no mean*mean
	    std::cout << "INFO:: " << moving_atom_positions.size() << " matched atoms had: \n"
		      << "   mean devi: " << mean << "\n"
		      << "    rms devi: " << sqrt(var) << "\n"
		      << "    max devi: " << maxd << "\n"
		      << "    min devi: " << mind << std::endl;

	    CAtom *at;
	    // We are going to delete the current atoms of the residue
	    // and add the std_base ones.  First what *are* the atom
	    // names we what to add or delete?
	    // 
	    std::vector<std::string> mol_base_atom_names;
	    if (mol_base_name == "Ar" || mol_base_name == "Ad")
	       mol_base_atom_names = adenine;
	    if (mol_base_name == "Gr" || mol_base_name == "Gd")
	       mol_base_atom_names = guanine;
	    if (mol_base_name == "Cr" || mol_base_name == "Cd")
	       mol_base_atom_names = cytosine;
	    if (mol_base_name == "Tr" || mol_base_name == "Td")
	       mol_base_atom_names = thymine;
	    if (mol_base_name == "Ur" || mol_base_name == "Ud")
	       mol_base_atom_names = uracil;

	    if (mol_base_atom_names.size() == 0) {
	       std::cout << "ERROR:: failed to find mol_base_name for mol_base_atom_names\n";
	    } else {
	       
	       std::vector<std::string> std_base_atom_names;
	       if (std_base_name == "Ar" || std_base_name == "Ad")
		  std_base_atom_names = adenine;
	       if (std_base_name == "Gr" || std_base_name == "Gd")
		  std_base_atom_names = guanine;
	       if (std_base_name == "Cr" || std_base_name == "Cd")
		  std_base_atom_names = cytosine;
	       if (std_base_name == "Tr" || std_base_name == "Td")
		  std_base_atom_names = thymine;
	       if (std_base_name == "Ur" || std_base_name == "Ud")
		  std_base_atom_names = uracil;
	    
	       if (std_base_atom_names.size() == 0) {
		  std::cout << "ERROR:: failed to find std_base_name for std_base_atom_names\n";
	       } else {

		  // now find the atoms of the given residue and apply
		  // the transformation and add them to the residue;
		  short int have_deleted = 0;
		  for (unsigned int iat=0; iat<mol_base_atom_names.size(); iat++) {
		     for (int i=0; i<n_mol_base_atoms; i++) {
			if (mol_base_atoms[i]) { 
			   if (mol_base_atom_names[iat] == mol_base_atoms[i]->name) {
// 			      std::cout << ".... Deleting Atom " << mol_base_atoms[i]->name
// 					<< std::endl;
			      residue->DeleteAtom(i);
			      mol_base_atoms[i] = NULL;
			      have_deleted = 1;
			      break;
			   }
			}
		     }
		  }
 		  if (have_deleted)
 		     residue->TrimAtomTable();
		  
		  for (unsigned int iat=0; iat<std_base_atom_names.size(); iat++) {
		     for (int i=0; i<n_std_base_atoms; i++) {
			if (std_base_atom_names[iat] == std_base_atoms[i]->name) {
			   clipper::Coord_orth p(std_base_atoms[i]->x,
						 std_base_atoms[i]->y,
						 std_base_atoms[i]->z);
			   clipper::Coord_orth pt = p.transform(rtop);
			   std::string ele = std_base_atom_names[iat].substr(0,2);
			   at = new CAtom;
// 			   std::cout << ".... Adding Atom " << std_base_atoms[i]->name
// 				     << std::endl;
			   at->SetCoordinates(pt.x(), pt.y(), pt.z(), 1.0, 20.0);
			   at->SetAtomName(std_base_atoms[i]->name);
			   at->SetElementName(ele.c_str());
			   std::string new_alt_conf("");
			   // force it down the atom's throat :) [is there a better way?]
			   strncpy(at->altLoc, new_alt_conf.c_str(), 2);
			   residue->AddAtom(at);
			   if (use_old_seg_id)
			      strcpy(at->segID, old_seg_id_for_residue_atoms.c_str());
			}
		     }
		  }
	    
	    
		  //       
		  std::string new_base_name = std_base_name;

		  if (mol_base_name.length() != 2) {
		     new_base_name = mol_base_name;
		  } else {
		     // the normal case
		     std::string debug_s = new_base_name;
// 		     new_base_name = std_base_name.substr(0,1) +
// 			mol_base_name.substr(1,2);
// 		     std::cout << "Changing base name from " << residue->name
// 			       << " to " << new_base_name << " not "
// 			       << debug_s << std::endl;
		  }

		  residue->SetResName(new_base_name.c_str());
		  residue->TrimAtomTable();

	       }
	    }
	 }
      }
   }
} 

std::vector<std::pair<coot::atom_spec_t, std::string> >
coot::util::gln_asn_b_factor_outliers(CMMDBManager *mol) {

   std::vector<std::pair<coot::atom_spec_t, std::string> > v;
   int imod = 1;
      
   CModel *model_p = mol->GetModel(imod);
   CChain *chain_p;

   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      PCResidue residue_p;
      CAtom *at;
      for (int ires=0; ires<nres; ires++) { 
	 residue_p = chain_p->GetResidue(ires);
	 int n_atoms = residue_p->GetNumberOfAtoms();
	 std::string residue_name(residue_p->GetResName());
	 if ((residue_name == "ASN") ||
	     (residue_name == "GLN")) {
	    float b_sum = 0.0;
	    float b_sum_sq = 0.0;
	    CAtom *oatom = NULL;
	    CAtom *natom = NULL;
	    CAtom *go_to_atom = NULL;
	    int n_residue_atoms = 0;
	    for (int iat=0; iat<n_atoms; iat++) {
	       at = residue_p->GetAtom(iat);
	       std::string altloc(at->altLoc);
	       if (altloc == "") {
		  std::string atom_name(at->GetAtomName());
		  
		  if (((residue_name == "GLN") &&
		       ((atom_name == " OE1") ||
			(atom_name == " NE2"))) ||
		      ((residue_name == "ASN") &&
		       ((atom_name == "OD1") ||
			(atom_name == "ND2")))) {
		     
		     if (residue_name == "GLN") { 
			if (atom_name == " OE1") {
			   oatom = at;
			}
			if (atom_name == " NE2") {
			   natom = at;
			}
		     }
		     if (residue_name == "ASN") { 
			if (atom_name == " OD1") {
			   oatom = at;
			}
			if (atom_name == " ND2") {
			   natom = at;
			}
		     }
		  } else {
		     // is a normal atom of the residue:
		     b_sum += at->tempFactor;
		     b_sum_sq += at->tempFactor * at->tempFactor;
		     n_residue_atoms++;
		  }
		  // find the atom to centre on when the button is
		  // clicked (that's why we need an atom spec vector,
		  // not a residue spec.
		  if (residue_name == "GLN") {
		     if (atom_name == " CD ")
			go_to_atom = at;
		  }
		  if (residue_name == "ASN") {
		     if (atom_name == " CG ")
			go_to_atom = at;
		  }
	       }
	    }
	    // end of residue atoms loop
	    if (oatom) {
	       if (natom) {
		  if (n_residue_atoms > 2) {
		     float mean = b_sum/float(n_residue_atoms);
		     float var = b_sum_sq/float(n_residue_atoms) - mean*mean;
		     float std_dev = sqrt(var);
		     float diff = (oatom->tempFactor - natom->tempFactor)/2.0;
		     // we are only interested in cases that have the
		     // O atom B-factor greater than the N atom
		     // B-factor because only they can be fixed by
		     // flipping.
		     if (diff > 0.0) { 
			float z = fabs(diff/std_dev);
			// std::cout << z << std::endl;
			// 		     std::cout << diff << "/sqrt(" << var << ") = " << z
			// 			       << "           " << b_sum << " " << b_sum_sq << " "
			// 			       << n_residue_atoms << std::endl;
			
			if (go_to_atom) { 
			   if (z > 2.25) {
			      coot::atom_spec_t as(go_to_atom);
			      std::string button_label = "Z score = ";
			      button_label += coot::util::float_to_string(z);
			      button_label += "   ";
			      button_label += go_to_atom->GetChainID();
			      button_label += " ";
			      button_label += int_to_string(go_to_atom->GetSeqNum());
			      button_label += " ";
			      button_label += go_to_atom->GetResName();
			      as.float_user_data = z;
			      std::pair<coot::atom_spec_t, std::string> p(as, button_label);
			      v.push_back(p);
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   if (v.size() > 1) {
      // sort v;
      std::sort(v.begin(), v.end(), coot::compare_atom_specs_user_float_in_pair);
   }
   return v;
}

bool
coot::compare_atom_specs_user_float_in_pair(const std::pair<coot::atom_spec_t, std::string> &a,
					    const std::pair<coot::atom_spec_t, std::string> &b) {

   return b.first.float_user_data < a.first.float_user_data ? 1 : 0;
}

bool
coot::util::residue_has_hydrogens_p(CResidue *res) {

   bool result = 0;

   if (res) { 
      PPCAtom residue_atoms;
      int natoms;
      res->GetAtomTable(residue_atoms, natoms);
      for (int iat=0; iat<natoms; iat++) {
	 CAtom *at = residue_atoms[iat];
	 std::string ele(at->element);
	 if ((ele == " H") || (ele == " D")) {
	    result = 1;
	    break;
	 }
      }
   }
   return result;
}

// angle in radians.
clipper::Coord_orth
coot::util::rotate_round_vector(const clipper::Coord_orth &direction,
				const clipper::Coord_orth &position,
				const clipper::Coord_orth &origin_shift,
				double angle) {
   
   clipper::Coord_orth unit_vec = clipper::Coord_orth(direction.unit());
   
   double l = unit_vec[0];
   double m = unit_vec[1];
   double n = unit_vec[2];

   double ll = l*l;
   double mm = m*m;
   double nn = n*n;
   double cosk = cos(angle);
   double sink = sin(angle);
   double I_cosk = 1.0 - cosk;
   
   // The Rotation matrix angle w about vector with direction cosines l,m,n.
   // 
   // ( l**2+(m**2+n**2)cos k     lm(1-cos k)-nsin k        nl(1-cos k)+msin k   )
   // ( lm(1-cos k)+nsin k        m**2+(l**2+n**2)cos k     mn(1-cos k)-lsin k   )
   // ( nl(1-cos k)-msin k        mn(1-cos k)+lsin k        n*2+(l**2+m**2)cos k )
   //
   // (Amore documentation) Thanks for that pointer EJD :).
   
   clipper::Mat33<double> r( ll+(mm+nn)*cosk,    l*m*I_cosk-n*sink,  n*l*I_cosk+m*sink,
			     l*m*I_cosk+n*sink,  mm+(ll+nn)*cosk,    m*n*I_cosk-l*sink,
			     n*l*I_cosk-m*sink,  m*n*I_cosk+l*sink,  nn+(ll+mm)*cosk );
   
   clipper::RTop_orth rtop(r, clipper::Coord_orth(0,0,0));
   return origin_shift + (position-origin_shift).transform(rtop);
}

// We ignore the issue of alt confs because from refinement/reg we
// will be only looking at a single stretch of amino acids, of a given
// alt conf (or blank).
// 
int
coot::util::count_cis_peptides(CMMDBManager *mol) {

   return cis_peptides_info_from_coords(mol).size();
}

std::vector<coot::util::cis_peptide_info_t>
coot::util::cis_peptides_info_from_coords(CMMDBManager *mol) {

   std::vector<coot::util::cis_peptide_info_t> v;
   
   int n_models = mol->GetNumberOfModels();
   if (n_models== 0) {
      return v;
   }
   
   int imod = 1;
   CModel *model_p = mol->GetModel(imod);
   CChain *chain_p;
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      //       std::cout << "DEBUGG:: in cis_peptides_info_from_coords ichain " << ichain
      // 		<< " of " << nchains << " chains " << std::endl;
      chain_p = model_p->GetChain(ichain);
      //       std::cout << "DEBUGG:: in cis_peptides_info_from_coords ichain " << ichain
      // << " chain_p " << chain_p << std::endl;
      if (chain_p) { 
	 int nres = chain_p->GetNumberOfResidues();
	 PCResidue residue_p_1;
	 PCResidue residue_p_2;
	 CAtom *at_1;
	 CAtom *at_2;
	 for (int ires=0; ires<(nres-1); ires++) { 
   
	    CAtom *ca_first = NULL, *c_first = NULL, *n_next = NULL, *ca_next = NULL;
	    residue_p_1 = chain_p->GetResidue(ires);
	    int n_atoms_1 = residue_p_1->GetNumberOfAtoms();
	    residue_p_2 = chain_p->GetResidue(ires+1);
	    int n_atoms_2 = residue_p_2->GetNumberOfAtoms();

	    if (residue_p_2->GetSeqNum() == (residue_p_1->GetSeqNum() + 1)) { 
	 
	       for (int iat=0; iat<n_atoms_1; iat++) {
		  at_1 = residue_p_1->GetAtom(iat);
		  if (std::string(at_1->GetAtomName()) == " CA ")
		     ca_first = at_1;
		  if (std::string(at_1->GetAtomName()) == " C  ")
		     c_first = at_1;
	       }

	       for (int iat=0; iat<n_atoms_2; iat++) {
		  at_2 = residue_p_2->GetAtom(iat);
		  if (std::string(at_2->GetAtomName()) == " CA ")
		     ca_next = at_2;
		  if (std::string(at_2->GetAtomName()) == " N  ")
		     n_next = at_2;
	       }
	    }
	 
	    if (ca_first && c_first && n_next && ca_next) {
	       clipper::Coord_orth caf(ca_first->x, ca_first->y, ca_first->z);
	       clipper::Coord_orth  cf( c_first->x,  c_first->y,  c_first->z);
	       clipper::Coord_orth can( ca_next->x,  ca_next->y,  ca_next->z);
	       clipper::Coord_orth  nn(  n_next->x,   n_next->y,   n_next->z);
	       double tors = clipper::Coord_orth::torsion(caf, cf, nn, can);
	       double torsion = clipper::Util::rad2d(tors);
	       double pos_torsion = (torsion > 0.0) ? torsion : 360.0 + torsion;
	       double distortion = fabs(180.0 - pos_torsion);
	       if (distortion > 90.0) {
		  coot::residue_spec_t rs1(residue_p_1);
		  coot::residue_spec_t rs2(residue_p_2); 
		  v.push_back(coot::util::cis_peptide_info_t(chain_p->GetChainID(),
							     rs1, rs2, imod, torsion));
	       }
	    }
	 } 
      }
   }
   return v;
} 


// remove wrong cis_peptides
void
coot::util::remove_wrong_cis_peptides(CMMDBManager *mol) {

#ifdef HAVE_MMDB_WITH_CISPEP
   
   std::vector<coot::util::cis_peptide_info_t> v_coords = 
      coot::util::cis_peptides_info_from_coords(mol);
//    std::cout << "INFO:: There were " << v_coords.size() << " CISPEPs from the coordinates"
// 	     << std::endl;


   PCCisPep       CisPep;
   int n_models = mol->GetNumberOfModels();
   for (int imod=1; imod<=n_models; imod++) { 
      std::vector<CCisPep> bad_cis_peptides;
      std::vector<CCisPep> good_cis_peptides;
      CModel *model_p = mol->GetModel(imod);
      int ncp = model_p->GetNumberOfCisPeps();
      for (int icp=1; icp<=ncp; icp++) {
	 CisPep = model_p->GetCisPep(icp);
	 if (CisPep)  {
// 	    std::cout << "mmdb:: " << " :" << CisPep->chainID1 << ": "<< CisPep->seqNum1 << " :" 
// 		      << CisPep->chainID2 << ": " << CisPep->seqNum2 << std::endl;
	    coot::util::cis_peptide_info_t cph(CisPep);

	    // Does that match any of the coordinates cispeps?
	    short int ifound = 0;
	    for (unsigned int iccp=0; iccp<v_coords.size(); iccp++) {
	       if (cph == v_coords[iccp]) {
		  // std::cout << " ......header matches" << std::endl;
		  ifound = 1;
		  break;
	       } else {
		  // std::cout << "       header not the same" << std::endl;
	       }
	    }
	    if (ifound == 0) {
	       // needs to be removed
	       std::cout << "INFO:: Removing CIS peptide from PDB header: " 
			 << cph.chain_id_1 << " "
			 << cph.resno_1 << " "
			 << cph.chain_id_2 << " "
			 << cph.resno_2 << " "
			 << std::endl;
	       bad_cis_peptides.push_back(*CisPep);
	    } else {
	       good_cis_peptides.push_back(*CisPep);
// 	       std::cout << "This CIS peptide was real: " 
// 			 << cph.chain_id_1 << " "
// 			 << cph.resno_1 << " "
// 			 << cph.chain_id_2 << " "
// 			 << cph.resno_2 << " "
// 			 << std::endl;
	    } 
	 }
      }
      if (bad_cis_peptides.size() > 0) {
	 // delete all CISPEPs and add back the good ones
	 model_p->RemoveCisPeps();
	 for (unsigned int igood=0; igood<good_cis_peptides.size(); igood++) {
	    PCCisPep good = new CCisPep;
	    *good = good_cis_peptides[igood];
	    model_p->AddCisPep(good);
	 }
      } 
   }
#endif // HAVE_MMDB_WITH_CISPEP
} 


CMMDBManager *
coot::mol_by_symmetry(CMMDBManager *mol, 
		      clipper::Cell cell, 
		      clipper::RTop_frac rtop_frac,
		      std::vector<int> pre_shift_to_origin_abc) {

   bool verbose_output = 0; // should be a passed param?
   
   CMMDBManager *mol2 = new CMMDBManager;
   mol2->Copy(mol, MMDBFCM_All);

   // Usually gets filled by GetTMatrix().
   mat44 mat_origin_shift; // shift needed to get close to the origin.
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 mat_origin_shift[i][j] = 0.0;
   for (int i=0; i<4; i++) mat_origin_shift[i][i] = 1.0;
   
   if (pre_shift_to_origin_abc.size() == 3) { 
      mol2->GetTMatrix(mat_origin_shift, 0,
		       pre_shift_to_origin_abc[0],
		       pre_shift_to_origin_abc[1],
		       pre_shift_to_origin_abc[2]);
   }
      
   // mol2->GetTMatrix(mat, symop_no, shift_a, shift_b, shift_c);

   clipper::Coord_orth origin_shift_orth(mat_origin_shift[0][3],
					 mat_origin_shift[1][3],
					 mat_origin_shift[2][3]);
   clipper::RTop_orth to_origin_rtop(clipper::Mat33<double>(1,0,0,0,1,0,0,0,1),
				     origin_shift_orth);

   clipper::RTop_orth rtop = rtop_frac.rtop_orth(cell);
   for(int imod = 1; imod<=mol2->GetNumberOfModels(); imod++) {
      CModel *model_p = mol2->GetModel(imod);
      CChain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 CResidue *residue_p;
	 CAtom *at;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    int n_atoms = residue_p->GetNumberOfAtoms();
	    for (int iat=0; iat<n_atoms; iat++) {
	       at = residue_p->GetAtom(iat);
	       clipper::Coord_orth co(at->x, at->y, at->z);
	       co -= origin_shift_orth;
	       clipper::Coord_orth to = co.transform(rtop);
	       to += origin_shift_orth;
// 	       std::cout << " atom from " 
// 			 << at->x << " " << at->y << " " << at->z << " " 
// 			 << " to " << to.format() << std::endl;
	       at->x = to.x(); at->y = to.y(); at->z = to.z(); 
	    }
	 }
      }
   }

   if (verbose_output) { 
      std::cout << "symmetry rtop_orth:\n" << rtop.format() << std::endl;
      std::cout << "symmetry rtop_frac:\n" << rtop_frac.format() << std::endl;
   }
   
   return mol2;
} 


// Perhaps this should be a class function of a class derived from CMMDBManager?
int
coot::write_coords_pdb(CMMDBManager *mol, const std::string &file_name) {


   coot::util::remove_wrong_cis_peptides(mol);
   int r = mol->WritePDBASCII(file_name.c_str());

   return r;
}


// convert atoms in residue to HETATMs
// 
int
coot::hetify_residue_atoms(CResidue *res) {

   int n = 0;
   if (res) { 
      PPCAtom residue_atoms;
      int natoms;
      res->GetAtomTable(residue_atoms, natoms);
      for (int iat=0; iat<natoms; iat++) {
	 CAtom *at = residue_atoms[iat];
	 at->Het = 1;
	 n++;
      } 
   }
   return n;
}

// convert atoms in residue to HETATMs
// 
int
coot::hetify_residue_atoms_as_needed(CResidue *res) {

   int r = 0;
   if (res) { 
      std::string resname = res->GetResName();
      if (! is_member_p(coot::util::PDB_standard_residue_types(), resname))
	 r = hetify_residue_atoms(res);
   }
   return r;
} 


// Interacting Residues: Return all residues of mol1, mol2 that
// have atoms that are closer that dist to atoms of mol2/mol1.
//
// Slightly horrifically, we have to create a new molecule from the 2
// given molecules so that we can run SeekContacts().  So we choose
// the first (non-blank) model of each of the molecules to go into
// MODEL 1 and MODEL 2 of the combined/new molecule.
//
// The returned CResidues are not in mol1 and mol2 - they are in combined_mol.
// 
std::pair<std::vector<CResidue *>, std::vector<CResidue *> >
coot::close_residues_from_different_molecules_t::close_residues(CMMDBManager *mol1,
								CMMDBManager *mol2, float dist) {

   std::vector<CResidue *> v1;
   std::vector<CResidue *> v2;

   if (mol1 && mol2) {

      combined_mol = new CMMDBManager;
      
      // combined_mol MODEL number 1 
      int n_models_mol_1 = mol1->GetNumberOfModels();
      for (int imod=1; imod<=n_models_mol_1; imod++) {
	 CModel *model_p = mol1->GetModel(imod);
	 if (model_p) {
	    CModel *new_model = new CModel;
	    new_model->Copy(model_p);
	    combined_mol->AddModel(new_model);
	    break;
	 } 
      } 

      // combined_mol MODEL number 2
      int n_models_mol_2 = mol2->GetNumberOfModels();
      for (int imod=1; imod<=n_models_mol_2; imod++) {
	 CModel *model_p = mol2->GetModel(imod);
	 if (model_p) {
	    CModel *new_model = new CModel;
	    new_model->Copy(model_p);
	    combined_mol->AddModel(new_model);
	    break;
	 } 
      } 

      int SelectionHandle_1 = combined_mol->NewSelection();
      PPCAtom atom_selection_1;
      combined_mol->SelectAtoms (SelectionHandle_1, 1, "*",
				 ANY_RES, // starting resno, an int
				 "*", // any insertion code
				 ANY_RES, // ending resno
				 "*", // ending insertion code
				 "*", // any residue name
				 "*", // atom name
				 "*", // elements
				 "*"  // alt loc.
				 );
      int n_selected_atoms_1;
      combined_mol->GetSelIndex(SelectionHandle_1, atom_selection_1, n_selected_atoms_1);
      
      int SelectionHandle_2 = combined_mol->NewSelection();
      PPCAtom atom_selection_2;
      combined_mol->SelectAtoms (SelectionHandle_2, 2, "*",
				 ANY_RES, // starting resno, an int
				 "*", // any insertion code
				 ANY_RES, // ending resno
				 "*", // ending insertion code
				 "*", // any residue name
				 "*", // atom name
				 "*", // elements
				 "*"  // alt loc.
				 );
      int n_selected_atoms_2;
      combined_mol->GetSelIndex(SelectionHandle_2, atom_selection_2, n_selected_atoms_2);

      std::cout << "INFO:: selected " << n_selected_atoms_1
		<< " from (copy of) 1st interaction molecule\n";
      std::cout << "INFO:: selected " << n_selected_atoms_2
		<< " from (copy of) 2nd interaction molecule\n";
      
      
      // (Sigh (of relief))...
      //
      // OK, now we can run SeekContacts();

      PSContact pscontact = NULL;
      int n_contacts;
      long i_contact_group = 1;
      mat44 my_matt;
      CSymOps symm;
      for (int i=0; i<4; i++) 
	 for (int j=0; j<4; j++) 
	    my_matt[i][j] = 0.0;      
      for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

      std::cout << "... SeekContacts() on " << n_selected_atoms_1
		<< " and " << n_selected_atoms_2 << " atoms" << std::endl;
      combined_mol->SeekContacts(atom_selection_1, n_selected_atoms_1,
				 atom_selection_2, n_selected_atoms_2,
				 0.0, dist,
				 1, pscontact, n_contacts,
				 0, &my_matt, i_contact_group);

      std::cout << "INFO:: Contacts between 2 molecules: found "
		<< n_contacts << " contacts" << std::endl;
      
      if (n_contacts > 0) {
	 if (pscontact) {
	    for (int i_contact=0; i_contact<n_contacts; i_contact++) {

	       CResidue *r1 = atom_selection_1[pscontact[i_contact].id1]->GetResidue();
	       CResidue *r2 = atom_selection_2[pscontact[i_contact].id2]->GetResidue();

	       CModel *model_1 = r1->GetModel();
	       CModel *model_2 = r2->GetModel();

	       if (model_1 != model_2) { 
	       
		  if (! coot::is_member_p(v1, r1))
		     v1.push_back(r1);
		  if (! coot::is_member_p(v2, r2))
		     v2.push_back(r2);
	       }
	    }
	 }
      } 
   }
   std::cout << "INFO:: interacting residues from molecules: "
	     << v1.size() << " and " << v2.size() << std::endl;
   return std::pair<std::vector<CResidue *>, std::vector<CResidue *> > (v1, v2);
}


// move waters round protein, fiddle with mol.
// return the number of moved waters.
int
coot::util::move_waters_around_protein(CMMDBManager *mol) {

   int n_moved = 0;
   std::vector<clipper::Coord_orth> protein_coords;
   std::vector<std::pair<CAtom*, clipper::Coord_orth> > water_atoms;

   // First we fill protein_atoms and water_atoms (water atoms are not
   // part of protein atoms)

   if (mol) { 
      int imod = 1;
      CModel *model_p = mol->GetModel(imod);
      CChain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 CResidue *residue_p;
	 CAtom *at = 0;
	 for (int ires=0; ires<nres; ires++) {
	    residue_p = chain_p->GetResidue(ires);
	    int n_atoms = residue_p->GetNumberOfAtoms();
	    std::string residue_name(residue_p->name);
	    if (residue_name == "WAT" ||
		residue_name == "HOH") {

	       for (int iat=0; iat<n_atoms; iat++) {
		  at = residue_p->GetAtom(iat);
		  if (! at->isTer()) { 
		     at = residue_p->GetAtom(iat);
		     clipper::Coord_orth c(at->x, at->y, at->z);
		     std::pair <CAtom *, clipper::Coord_orth> pair(at, c);
		     water_atoms.push_back(pair);
		  }
	       }
	    } else {
	       for (int iat=0; iat<n_atoms; iat++) {
		  at = residue_p->GetAtom(iat);
		  if (! at->isTer()) {
		     std::string ele(at->element);
		     if (ele  != " C") { 
			clipper::Coord_orth pt(at->x, at->y, at->z);
			protein_coords.push_back(pt);
		     }
		  }
	       }
	    } 
	 }
      }
   }

   // OK, so waters_atoms and protein atoms are filled.

   try { 
      // Now clipperize the variables.
      std::pair<clipper::Cell,clipper::Spacegroup> csp = get_cell_symm(mol);
      clipper::Cell cell = csp.first;
      clipper::Spacegroup spacegroup = csp.second;

      if (cell.is_null()) {
	 std::cout << "WARNING:: null cell in move_waters_around_protein" << std::endl;
      } else {
	 if (spacegroup.is_null()) { 
	    std::cout << "WARNING:: null spgr in move_waters_around_protein" << std::endl;
	 } else {

	    std::vector<std::pair<CAtom*, clipper::Coord_orth> > water_atoms_moved =
	       symmetry_move_atoms(protein_coords, water_atoms, cell, spacegroup);

	    for (unsigned int iw=0; iw<water_atoms_moved.size(); iw++) {
	       if (water_atoms_moved[iw].first) {
		  water_atoms_moved[iw].first->x = water_atoms_moved[iw].second.x();
		  water_atoms_moved[iw].first->y = water_atoms_moved[iw].second.y();
		  water_atoms_moved[iw].first->z = water_atoms_moved[iw].second.z();
		  n_moved++;
	       }
	    }
	 }
      }
   }
   catch (std::runtime_error rte) {
      std::cout << rte.what() << std::endl;
   }

   return n_moved;
}

// Return waters atoms of the same size as the input, except if the
// first is NULL, then the atom need not move, if it is not null, then
// the water atom (first) should be moved to the second position.
// 
std::vector<std::pair<CAtom *, clipper::Coord_orth> >
coot::util::symmetry_move_atoms(const std::vector<clipper::Coord_orth> &protein_coords,
				const std::vector<std::pair<CAtom*, clipper::Coord_orth> > &water_atoms_in,
				clipper::Cell cell,
				clipper::Spacegroup spacegroup) {


   clipper::Coord_frac pre_shift_frac = coot::util::shift_to_origin(protein_coords, cell, spacegroup);
   clipper::Coord_orth pre_shift_orth = pre_shift_frac.coord_orth(cell);

   std::vector<std::pair<CAtom*, clipper::Coord_orth> > water_atoms = water_atoms_in;
   
   if (0) 
      std::cout << "DEBUG:: pre_shift_frac " << pre_shift_frac.format()
		<< " pre_shift_orth " << pre_shift_orth.format()
		<< std::endl;

   // create shifted protein coords
   std::vector<clipper::Coord_orth> protein_coords_origin_shifted(protein_coords.size());
   for (unsigned int ip=0; ip<protein_coords.size(); ip++) { 
      protein_coords_origin_shifted[ip] =
	 protein_coords[ip] + pre_shift_orth;
      // 	 if (ip < 20)
      // 	    std::cout << "  shifting "
      // 		      << protein_coords[ip].format() << " by "
      // 		      << pre_shift_orth.format() << " gives "
      // 		      << protein_coords_origin_shifted[ip].format()
      // 		      << std::endl;
   }

   // Do the cell shift search
   int n = spacegroup.num_symops();
   clipper::Coord_frac cell_shift; 
   for (unsigned int iw=0; iw<water_atoms.size(); iw++) {
      clipper::Coord_orth water_pos_pre(water_atoms[iw].second);
      clipper::Coord_orth water_pos = translate_close_to_origin(water_pos_pre, cell);

      // std::cout << " water_pos " << water_pos.format() << std::endl;
      double d_best = 99999999.9;
      // The compiler doesn't like rtop_best being used below
      // without being initialized properly here.
      // clipper::RTop_orth rtop_best; // old
      clipper::RTop_orth rtop_best(clipper::Mat33<double>(1,0,0,0,1,0,0,0,1), clipper::Coord_orth(0,0,0));
      bool improved = 0;
      // 
      for (int isym=0; isym<n; isym++) {
	 for (int x_shift = -1; x_shift<2; x_shift++) { 
	    for (int y_shift = -1; y_shift<2; y_shift++) { 
	       for (int z_shift = -1; z_shift<2; z_shift++) {
		  cell_shift = clipper::Coord_frac(x_shift, y_shift, z_shift); 
		  clipper::RTop_orth orthop = clipper::RTop_frac(spacegroup.symop(isym).rot(), spacegroup.symop(isym).trn() + cell_shift).rtop_orth(cell);
		  clipper::Coord_orth t_point = water_pos.transform(orthop);
		  double t_dist = coot::util::min_dist_to_points(t_point, protein_coords_origin_shifted);
		  if (t_dist < d_best) {
		     // std::cout << " better dist " << t_dist << std::endl;
		     d_best = t_dist;
		     rtop_best = orthop;
		     improved = 1;
		  }
	       }
	    }
	 }
      }

      if (improved) { 
	 // Apply the transformation then.
	 clipper::Coord_orth t_point = water_pos.transform(rtop_best);
	 water_atoms[iw].second = t_point - pre_shift_orth;
      } else {
	 water_atoms[iw].first = 0; // NULL, don't move it.
      }
   }
   return water_atoms;
}



// Throw an std::runtime_error exception on
// not-able-to-extract-cell/symm-info.  (In such a case, we convert a
// clipper::Message_base to a std::runtime_error).
// 
std::pair<clipper::Cell, clipper::Spacegroup>
coot::util::get_cell_symm(CMMDBManager *mol) {

   // Now clipperize the variables.

   mat44 my_matt;
   int err = mol->GetTMatrix(my_matt, 0, 0, 0, 0);
   if (err != 0) {
      std::string mess = "No symmetry available";
      throw std::runtime_error(mess);
   } else { 
      try { 
	 const clipper::MMDBManager* pcmmdb =
	    static_cast<const clipper::MMDBManager*>( mol );

	 clipper::Spacegroup spacegroup(pcmmdb->spacegroup());
	 clipper::Cell cell(pcmmdb->cell());
	 if (spacegroup.is_null())
	    std::cout << "Null clipper spacegroup from " << mol->GetSpaceGroup()
		      << std::endl;
	 if (cell.is_null())
	    std::cout << "Null clipper cell  " << std::endl;
	 return std::pair<clipper::Cell, clipper::Spacegroup> (cell, spacegroup);
      }
      catch (clipper::Message_base except) {
	 std::string message = "Fail to make clipper::Spacegroup from ";
	 message += mol->GetSpaceGroup();
	 throw std::runtime_error(message);
      }
   }
   std::cout << "got to here - bad! in get_cell_symm()"
	     << std::endl;
}


// shove a cell from a clipper cell into the passed mol.
bool
coot::util::set_mol_cell(CMMDBManager *mol, clipper::Cell cell_local) {

   bool status = 0; 
   mol->SetCell(cell_local.a(), cell_local.b(), cell_local.c(),
		clipper::Util::rad2d(cell_local.alpha()),
		clipper::Util::rad2d(cell_local.beta()),
		clipper::Util::rad2d(cell_local.gamma()));

   realtype cell[6], vol;
   int orthog;
   
   mol->GetCell(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5], vol, orthog);
   if (fabs(cell[0] - cell_local.a()) < 0.1)
      if (fabs(cell[1] - cell_local.b()) < 0.1)
	 if (fabs(cell[2] - cell_local.c()) < 0.1)
	    if (fabs(clipper::Util::d2rad(cell[3]) - cell_local.alpha()) < 0.1)
	       if (fabs(clipper::Util::d2rad(cell[4]) - cell_local.beta()) < 0.1)
		  if (fabs(clipper::Util::d2rad(cell[5]) - cell_local.gamma()) < 0.1)
		     status = 1;

   return status;
}

//
clipper::Mat33<double>
coot::util::residue_orientation(CResidue *residue_p, const clipper::Mat33<double> &orientation_in) {

   std::vector<clipper::Coord_orth> pts;
   PPCAtom residue_atoms = 0;
   int n_residue_atoms;
   clipper::Mat33<double> r = orientation_in;
   clipper::Coord_orth n_vec(0,0,1);
   CAtom *ca = 0;
   CAtom *n = 0;
   
   if (residue_p) { 
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (unsigned int i=0; i<n_residue_atoms; i++) {
	 if (!is_main_chain_p(residue_atoms[i])) {
	    pts.push_back(clipper::Coord_orth(residue_atoms[i]->x,
					      residue_atoms[i]->y,
					      residue_atoms[i]->z));
	 } else {
	    std::string atom_name(residue_atoms[i]->name);
	    if (atom_name == " CA ")
	       ca = residue_atoms[i];
	    if (atom_name == " N  ")
	       n = residue_atoms[i];
	 } 
      }
      
      if (pts.size() > 0) {

	 if (ca) { 
	    clipper::Coord_orth ca_pos(ca->x, ca->y, ca->z);
	    clipper::Coord_orth average_pos = coot::util::average_position(pts);
	    clipper::Coord_orth u((average_pos-ca_pos).unit());

	    // reset n_vect to something sensible, if we have the CA and N.
	    if (ca && n) {
	       clipper::Coord_orth  n_pos( n->x,  n->y,  n->z);
	       n_vec = n_pos - ca_pos;
	    }
	 
	    // now make a mat...
	    //
	    clipper::Coord_orth n_vec_unit(n_vec.unit());
	    
	    clipper::Coord_orth p1(clipper::Coord_orth::cross(n_vec_unit, u).unit());
	    clipper::Coord_orth p2(clipper::Coord_orth::cross( p1, u).unit());
	    clipper::Coord_orth p3 = u;

	    std::cout << "primary sidechain direction: " << u.format() << std::endl;
	    r = clipper::Mat33<double> (p1.x(), p1.y(), p1.z(),
					p2.x(), p2.y(), p2.z(),
					p3.x(), p3.y(), p3.z());
	    std::cout << r.format() << std::endl;
	    std::cout << "determinant: " << r.det() << std::endl;
	    
	 }
      }
   }
   return r;
} 


// 
clipper::Coord_orth
coot::util::average_position(std::vector<clipper::Coord_orth> &pts) {
   
   if (pts.size() > 0) {
      double xsum = 0.0;
      double ysum = 0.0;
      double zsum = 0.0;
      for (unsigned int i=0; i<pts.size(); i++) { 
	 xsum += pts[i].x();
	 ysum += pts[i].y();
	 zsum += pts[i].z();
      }
      double denom=1.0/double(pts.size());
      return clipper::Coord_orth(denom*xsum, denom*ysum, denom*zsum);
   } else {
      return clipper::Coord_orth(0,0,0);
   }
}



// caller must check that others has some points in it.
// 
double
coot::util::min_dist_to_points(const clipper::Coord_orth &pt,
			       const std::vector<clipper::Coord_orth> &others) {

   double best_dist = 9999999.9;
   for (unsigned int i=0; i<others.size(); i++) {
      double d = (pt - others[i]).lengthsq();
      if (d<best_dist) {
	 best_dist = d;
      }
   }
   return sqrt(best_dist);
}


// Return the fractional shift needed to translate the protein
// as close as possible to the origin (do not apply the shift).
// 
//
// Throw a clipper::Message_base exception (e.g. no cell or symmetry).
//
// Can throw a std::runtime_error other times.
// 
clipper::Coord_frac
coot::util::shift_to_origin(CMMDBManager *mol) {

   // Throw a clipper::Message_base exception on no cell or symmetry.
   std::pair<clipper::Cell, clipper::Spacegroup> csp = get_cell_symm(mol);
   clipper::Cell cell = csp.first;
   clipper::Spacegroup spacegroup = csp.second;

   // Throws an exception
   clipper::Coord_orth median_pos = median_position(mol);

   clipper::Coord_frac mpf = median_pos.coord_frac(cell);
   clipper::Coord_frac rf (round(-mpf.u()), round(-mpf.v()), round(-mpf.w()));
   return rf;
}

clipper::Coord_frac
coot::util::shift_to_origin(const std::vector<clipper::Coord_orth> &protein_coords,
			    clipper::Cell cell,
			    clipper::Spacegroup spacegroup) {

   clipper::Coord_orth median_pos = median_position(protein_coords);
   clipper::Coord_frac mpf = median_pos.coord_frac(cell);
   clipper::Coord_frac rf (round(-mpf.u()), round(-mpf.v()), round(-mpf.w()));
   return rf;
}


// Can throw a std::runtime_error
// 
clipper::Coord_orth
coot::util::median_position(const std::vector<clipper::Coord_orth> &pts) {

   if (pts.size() == 0 ) {
      std::string message = "No atoms in molecule - no mediain position";
      throw std::runtime_error(message);
   }
   
   std::vector<float> pts_x;
   std::vector<float> pts_y;
   std::vector<float> pts_z;
   for (unsigned int i=0; i<pts.size(); i++) {
      pts_x.push_back(pts[i].x());
      pts_y.push_back(pts[i].y());
      pts_z.push_back(pts[i].z());
   }
   std::sort(pts_x.begin(), pts_x.end());
   std::sort(pts_y.begin(), pts_y.end());
   std::sort(pts_z.begin(), pts_z.end());
   unsigned int mid_index = pts_x.size()/2;
   return clipper::Coord_orth(pts_x[mid_index], pts_y[mid_index], pts_z[mid_index]);
}



// Can throw a std::runtime_error
// 
clipper::Coord_orth
coot::util::median_position(CMMDBManager *mol) {

   std::vector<float> pts_x;
   std::vector<float> pts_y;
   std::vector<float> pts_z;

   // for(int imod = 1; imod<=asc.mol->GetNumberOfModels(); imod++) {
   int imod = 1;
   CModel *model_p = mol->GetModel(imod);
   CChain *chain_p;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      CResidue *residue_p;
      CAtom *at;
      for (int ires=0; ires<nres; ires++) { 
	 residue_p = chain_p->GetResidue(ires);
	 int n_atoms = residue_p->GetNumberOfAtoms();
	 
	 for (int iat=0; iat<n_atoms; iat++) {
	    at = residue_p->GetAtom(iat);
	    if (! at->isTer()) {
	       pts_x.push_back(at->x);
	       pts_y.push_back(at->y);
	       pts_z.push_back(at->z);
	    }
	 }
      }
   }

   if (pts_x.size() == 0) {
      std::string message = "No atoms in molecule - no mediain position";
      throw std::runtime_error(message);
   }

   std::sort(pts_x.begin(), pts_x.end());
   std::sort(pts_y.begin(), pts_y.end());
   std::sort(pts_z.begin(), pts_z.end());
   unsigned int mid_index = pts_x.size()/2;

   return clipper::Coord_orth(pts_x[mid_index], pts_y[mid_index], pts_z[mid_index]);
} 



//
clipper::Coord_orth
coot::util::translate_close_to_origin(const clipper::Coord_orth pos,
				      const clipper::Cell &cell) {

   clipper::Coord_frac cf = pos.coord_frac(cell);
   clipper::Coord_frac cfi(round(-cf.u()), round(-cf.v()), round(-cf.w()));
   return pos + cfi.coord_orth(cell);
} 

void
coot::util::print_secondary_structure_info(CModel *model_p) {
   
   // secondary structure information
   //
   int nhelix = model_p->GetNumberOfHelices();
   int nsheet = model_p->GetNumberOfSheets();
   std::cout << "INFO:: There are " << nhelix << " helices and "
	     << nsheet << " sheets\n";
   PCHelix helix_p;
   PCSheet sheet_p;
   PCStrand strand_p;

   std::cout << "               Helix info: " << std::endl;
   std::cout << "------------------------------------------------\n";
   for (int ih=1; ih<=nhelix; ih++) {
      helix_p = model_p->GetHelix(ih);
      if (helix_p) { 
      std::cout << helix_p->serNum << " " << helix_p->helixID << " "
		<< helix_p->initChainID << " " << helix_p->initSeqNum
		<< " " << helix_p->endChainID << " " << helix_p->endSeqNum
		<< helix_p->length << " " << helix_p->comment << std::endl;
      } else {
	 std::cout << "ERROR: no helix!?" << std::endl;
      }
   }
   std::cout << "               Sheet info: " << std::endl;
   std::cout << "------------------------------------------------\n";
   for (int is=1; is<=nsheet; is++) {
      sheet_p = model_p->GetSheet(is);

      int nstrand = sheet_p->nStrands;
      for (int istrand=0; istrand<nstrand; istrand++) {
	 strand_p = sheet_p->Strand[istrand];
	 if (strand_p) { 
	    std::cout << strand_p->sheetID << " " << strand_p->strandNo << " "
		      << strand_p->initChainID << " " << strand_p->initSeqNum
		      << " " << strand_p->endChainID << " " << strand_p->endSeqNum
		      << std::endl;
	 }
      }
   }
   std::cout << "------------------------------------------------\n";
}
   

// return success status as first element
// 
std::pair<bool, clipper::Coord_orth>
coot::centre_of_molecule(CMMDBManager *mol) {

   bool status = 0;
   clipper::Coord_orth centre(0,0,0);

   if (mol) {

      int n_atoms = 0;
      double xs=0, ys=0, zs=0;

      for(int imod=1; imod<=mol->GetNumberOfModels(); imod++) {
	 CModel *model_p = mol->GetModel(imod);
	 CChain *chain_p;
	 int nchains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<nchains; ichain++) {
	    chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    CResidue *residue_p;
	    CAtom *at;
	    for (int ires=0; ires<nres; ires++) { 
	       residue_p = chain_p->GetResidue(ires);
	       int n_residue_atoms = residue_p->GetNumberOfAtoms();
	 
	       for (int iat=0; iat<n_residue_atoms; iat++) {
		  at = residue_p->GetAtom(iat);
		  if (! at->isTer()) { 
		     xs += at->x;
		     ys += at->y;
		     zs += at->z;
		     n_atoms++;
		  }
	       }
	    }
	 }
      }
      
      if (n_atoms > 0) {
	 status = 1;
	 double dna = static_cast<double> (n_atoms);
	 centre = clipper::Coord_orth(xs/dna, ys/dna, zs/dna);
      }
   }

   return std::pair<bool, clipper::Coord_orth> (status, centre);
} 

CResidue *
coot::nearest_residue_by_sequence(CMMDBManager *mol,
				  const residue_spec_t &spec) {

   CResidue *r = NULL;
   int resno_closest_high = -9999;
      
   if (mol) {
      int imod = 1;
      CModel *model_p = mol->GetModel(imod);
      CChain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 std::string chain_id = chain_p->GetChainID();
	 if (chain_id == spec.chain) { 
	    int nres = chain_p->GetNumberOfResidues();
	    CResidue *residue_p;
	    for (int ires=0; ires<nres; ires++) { 
	       residue_p = chain_p->GetResidue(ires);
	       int this_resno = residue_p->GetSeqNum();
	       if (labs(spec.resno - this_resno)
		   < (labs(spec.resno - resno_closest_high))) {
		  resno_closest_high = this_resno;
		  r = residue_p;
	       }
	    }
	 }
      }
   }
   return r;
}




// copy cell, symm, origin and scale cards from m1 to m2 (if possible)
bool
coot::util::copy_cell_and_symm_headers(CMMDBManager *m1, CMMDBManager *m2) {

   bool r = 0;

   if (m1 && m2) { 

      //       realtype a[6];
      //       realtype vol;
      //       int orthcode;
      
      //       m1->GetCell(a[0], a[1], a[2], a[3], a[4], a[5], vol, orthcode);
      //       char *sg = m1->GetSpaceGroup();
      //       m2->SetSpaceGroup(sg);
      //       m2->SetCell(a[0], a[1], a[2], a[3], a[4], a[5]);
      
      m2->Copy(m1, MMDBFCM_Cryst);
      r = 1;
   }
   return r;
}

