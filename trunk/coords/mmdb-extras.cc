/* coords/mmdb-extras.cc
 * 
 * Copyright 2006 by The University of York
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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

#include <stdlib.h>  // for getenv()
#include <string.h>  // for memcpy()
#include <string>
#include <vector> // for Cartesian


#include "coot-coord-utils.hh"

#include "mmdb-extras.h" 
#include "mmdb.h"

#include "coot-sysdep.h"


float max_bond_length(const std::string &element) { 

   if (element == " S") { 
      return 2.3;
   } else { 
      return 1.8;
   }
};

int
atom_colour(const std::string &element) {

   if (element == " C") {
      return yellow;
   } else {
      if (element == " N") {
	 return blue;
      } else {
	 if (element == " O") {
	    return red;
	 } else {
	    if (element == " S") {
	       return green;
	    } else {
	       if (element == " H") {
		  return grey;
	       }
	    }
	 }
      }
   }
   return 5;
}


void
debug_atom_selection_container(atom_selection_container_t asc) {

   //
   PCAtom ap;
   
   cout << "DEBUG: asc " << "mol=" << asc.mol << endl;
   cout << "DEBUG: asc " << "n_selected_atoms=" << asc.n_selected_atoms << endl;
   cout << "DEBUG: asc " << "atom_selection=" << asc.atom_selection << endl;
   cout << "DEBUG: asc " << "read_error_message=" << asc.read_error_message << endl;
   cout << "DEBUG: asc " << "read_success=" << asc.read_success << endl;

   cout << "DEBUG: asc " << "cell="
	<< asc.mol->get_cell_p()->a << " "
	<< asc.mol->get_cell_p()->b << " "
	<< asc.mol->get_cell_p()->c << " "
	<< asc.mol->get_cell_p()->alpha << " "
	<< asc.mol->get_cell_p()->beta << " "
	<< asc.mol->get_cell_p()->gamma << endl;
   
   cout << "DEBUG: asc " << "spacegroup=" << asc.mol->get_cell_p()->spaceGroup
	<< endl;

   if (asc.n_selected_atoms > 10) {
      cout << "DEBUG start 10 atoms: " << endl;
      for (int ii = 0; ii< 10; ii++) { 
	 cout << ii << " " << asc.atom_selection[ii] << " " ; 
	 ap = asc.atom_selection[ii];
	 cout << *ap << endl;
      }
      cout << "DEBUG end 10 atoms: " << endl;
      for (int ii = asc.n_selected_atoms - 10; ii< asc.n_selected_atoms; ii++) { 
	 cout << ii << " " << asc.atom_selection[ii] << " " ;
	 ap = asc.atom_selection[ii];
	 cout << *ap << endl;
      }
   }
}

// needs <stdlib.h>

// CCP4 symmetry library checking:
// 
// Return 0 on failure
// 
int check_ccp4_symm() {

   int i = 1;
   char *s = getenv("CLIBD");
   if ( s == NULL ) {
      std::cout << "WARNING: Unable to find CLIBD, symmetry will fail!\n";
      i = 0;
   }
   return i;
} 

atom_selection_container_t
make_asc(CMMDBManager *mol) {

   atom_selection_container_t asc;
   asc.mol = (MyCMMDBManager *) mol;
   
   asc.SelectionHandle = mol->NewSelection();
   asc.mol->SelectAtoms (asc.SelectionHandle, 0, "*",
		     ANY_RES, // starting resno, an int
		     "*", // any insertion code
		     ANY_RES, // ending resno
		     "*", // ending insertion code
		     "*", // any residue name
		     "*", // atom name
		     "*", // elements
		     "*"  // alt loc.
		     );
   int nSelAtoms;
   asc.mol->GetSelIndex(asc.SelectionHandle, asc.atom_selection, asc.n_selected_atoms);

   int uddHnd = mol->RegisterUDInteger(UDR_ATOM , "atom index");
   if (uddHnd < 0) {
      std::cout << " atom index registration failed.\n";
   } else {
      asc.UDDAtomIndexHandle = uddHnd; 
      for (int i=0; i<asc.n_selected_atoms; i++)
	 asc.atom_selection[i]->PutUDData(uddHnd,i);
   }
   asc.read_error_message = "No error";
   asc.read_success = 1;
   asc.UDDOldAtomIndexHandle = -1;

   return asc;
} 

// Note, we also create a chain and add this residue to that chain.
// We do this so that we have a holder for the segid.
// 
// whole_residue_flag: only copy atoms that are either in this altLoc,
// or has an altLoc of "".
// 
CResidue *
coot::deep_copy_this_residue(CResidue *residue,
			     const std::string &altconf,
			     short int whole_residue_flag,
			     int atom_index_handle) {

   // Horrible casting to CResidue because GetSeqNum and GetAtomTable
   // are not const functions.
   // 
   CResidue *rres = new CResidue;
   CChain   *chain_p = new CChain;
   chain_p->SetChainID(((CResidue *)residue)->GetChainID());
   rres->SetResID(residue->GetResName(),
		  residue->GetSeqNum(),
		  residue->GetInsCode());

   PPCAtom residue_atoms;
   int nResidueAtoms;
   ((CResidue *)residue)->GetAtomTable(residue_atoms, nResidueAtoms);
   CAtom *atom_p;
   
   for(int iat=0; iat<nResidueAtoms; iat++) {
      std::string this_atom_alt_loc(residue_atoms[iat]->altLoc);
      if (whole_residue_flag ||
	  this_atom_alt_loc  == altconf || this_atom_alt_loc == "") { 
	 atom_p = new CAtom;
	 atom_p->Copy(residue_atoms[iat]);
	 int i_add = rres->AddAtom(atom_p);
      }
   }
   chain_p->AddResidue(rres);
   return rres;
}

// Note, we also create a chain and add this residue to that chain.
// We do this so that we have a holder for the segid.
// 
// whole_residue_flag: only copy atoms that are either in this altLoc,
// or has an altLoc of "".
// 
std::pair<CResidue *, atom_selection_container_t>
coot::deep_copy_this_residue_and_make_asc(CMMDBManager *orig_mol,
					  CResidue *residue,
					  const std::string &altconf,
					  short int whole_residue_flag,
					  int atom_index_handle,
					  int udd_afix_handle) {

//    std::cout << "DEbbug:: in deep_copy_this_residue_and_make_asc is "
// 	     << udd_afix_handle << std::endl;

   // Horrible casting to CResidue because GetSeqNum and GetAtomTable
   // are not const functions.
   // 
   CResidue *rres = new CResidue;
   CChain   *chain_p = new CChain;
   std::string chain_id1 = ((CResidue *)residue)->GetChainID();
   chain_p->SetChainID(chain_id1.c_str());
   rres->seqNum = ((CResidue *)residue)->GetSeqNum();
   /* Copy insertion code - a char[10] */
   memcpy(rres->insCode, residue->insCode, sizeof(InsCode));
   memcpy(rres->name, residue->name, sizeof(ResName));
   
   PPCAtom residue_atoms;
   int nResidueAtoms;
   ((CResidue *)residue)->GetAtomTable(residue_atoms, nResidueAtoms);
   CAtom *atom_p;
   short int do_shelx_afix_data_flag = 0;
   if (udd_afix_handle >= 0)
      do_shelx_afix_data_flag = 1;
   
   for(int iat=0; iat<nResidueAtoms; iat++) {
      std::string this_atom_alt_loc(residue_atoms[iat]->altLoc);
      if (whole_residue_flag ||
	  this_atom_alt_loc  == altconf || this_atom_alt_loc == "") { 
	 atom_p = new CAtom;
	 atom_p->Copy(residue_atoms[iat]);
	 int i_add = rres->AddAtom(atom_p);
	 if (atom_index_handle > -1 ) { 
	    // copy across the atom indices:
	    
	 }
// 	 if (do_shelx_afix_data_flag) {
// 	    int ic;
// 	    if (residue_atoms[iat]->GetUDData(udd_afix_handle, ic) == UDDATA_Ok) {
// 	       std::cout << "Copying afix handle " << ic << " from " << residue_atoms[iat]
// 			 << "\n                        to " << atom_p << std::endl;
// 	       int istat = atom_p->PutUDData(udd_afix_handle, ic);
// 	       if (istat != UDDATA_Ok) { 
// 		  std::cout << "Opps! putuddata returns status " << istat << " vs "
// 			    << UDDATA_Ok << std::endl;
// 	       }
// 	    }
// 	 }
      }
   }
   chain_p->AddResidue(rres);

   // urgh.  This consting business makes a mess.
   CResidue *r = (CResidue *)residue;
   PCResidue *SelResidues = &r;
   std::pair<CMMDBManager *, int> mol_i =
      coot::util::create_mmdbmanager_from_res_selection(orig_mol, SelResidues,
							1, 0, 0, altconf, chain_id1, 1);

   atom_selection_container_t asc = make_asc(mol_i.first);
   int udd_afix_handle_inter = mol_i.first->GetUDDHandle(UDR_ATOM, "shelx afix");
//    std::cout << "DEBUG:: deep_copy_this_residue_and_make_asc got udd_afix_handle_inter : "
// 	     << udd_afix_handle_inter << std::endl;
   for (int i=0; i<asc.n_selected_atoms; i++) {
      int afix_number = -1;
      if (asc.atom_selection[i]->GetUDData(udd_afix_handle_inter, afix_number) == UDDATA_Ok)
	 std::cout << asc.atom_selection[i] << " has afix number " << afix_number
		   << std::endl;
//       else
// 	 std::cout << asc.atom_selection[i]
// 		   << " Failed get udd afix number in deep_copy_this_residue_and_make_asc"
// 		   << std::endl;
   }
   return std::pair<CResidue *, atom_selection_container_t>(rres, asc);
}


short int
coot::progressive_residues_in_chain_check(const CChain *chain_p) {

   int nres = ((CChain*)chain_p)->GetNumberOfResidues();
   CResidue *res_p;
   int previous_seq_num = -9999;
   int this_seq_no;

   for (int ires=0; ires<nres; ires++) { 
      res_p = ((CChain*)chain_p)->GetResidue(ires);
      if (res_p) {
	 this_seq_no = res_p->GetSeqNum();
	 if ( ! (this_seq_no >= (previous_seq_num+1)) )
	    return 0; 
      } else {
	 std::cout << "ERROR: null residue in progressive_residues_in_chain_check\n";
	 return 0;
      }
      previous_seq_num = this_seq_no;
   }
   return 1;
} 


coot::contact_info
coot::getcontacts(const atom_selection_container_t &asc) {

   PSContact pscontact = NULL;
   int n_contacts;
   float min_dist = 0.1;
   float max_dist = 1.9; // CB->SG CYS 1.8A
   long i_contact_group = 1;
   mat44 my_matt;
   CSymOps symm;
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

   asc.mol->SeekContacts(asc.atom_selection, asc.n_selected_atoms,
			 asc.atom_selection, asc.n_selected_atoms,
			 min_dist, max_dist, // min, max distances
			 0,        // seqDist 0 -> in same res also
			 pscontact, n_contacts,
			 0, &my_matt, i_contact_group);

   coot::contact_info ci(pscontact, n_contacts);

   // Now, do we need to handle MSE extra bonds?
   // 
   if (std::string(asc.atom_selection[0]->GetResName()) == "MSE") { 
      ci.add_MSE_Se_bonds(asc);
   }
   
   delete [] pscontact;
   return ci;
}

void
coot::contact_info::add_MSE_Se_bonds(const atom_selection_container_t &asc) {

   int SE_index = -1;
   int CE_index = -1;
   int CG_index = -1;
   for (int i=0; i<asc.n_selected_atoms; i++) {
      std::string atom_name = asc.atom_selection[i]->name;
      if (atom_name == "SE  ") SE_index = i;
      if (atom_name == " CE ") CE_index = i;
      if (atom_name == " CG ") CG_index = i;
   }
   if (SE_index != -1) { 
      if (CE_index != -1) { 
	 if (CG_index != -1) {
	    contacts.push_back(coot::contact_info::contacts_pair(CG_index, SE_index));
	    contacts.push_back(coot::contact_info::contacts_pair(SE_index, CE_index));
	 }
      }
   }
} 



// Typically this is used on an asc (moving atoms) to get the N of a
// peptide (say).  Return NULL on atom not found.
// 
CAtom *
coot::get_first_atom_with_atom_name(const std::string &atomname, 
				    const atom_selection_container_t &asc) { 

   CAtom *atom = NULL;

   for (int i=0; i<asc.n_selected_atoms; i++) { 
      std::string name(asc.atom_selection[i]->name);
      if (name == atomname) { 
	 atom = asc.atom_selection[i];
	 break;
      } 
   } 

   return atom;
} 

// tinker with asc
void 
coot::add_atom_index_udd_as_old(atom_selection_container_t asc) { 

   int old_atom_index;
   if (asc.n_selected_atoms > 0) {
      if (asc.UDDAtomIndexHandle >= 0) { 
	 int uddHnd = asc.mol->RegisterUDInteger(UDR_ATOM , "old atom index");
	 if (uddHnd >= 0) { 
	    asc.UDDOldAtomIndexHandle = uddHnd;
	    for (int i=0; i<asc.n_selected_atoms; i++) { 
	       if (asc.atom_selection[i]->GetUDData(asc.UDDAtomIndexHandle, old_atom_index) == UDDATA_Ok) { 
		  asc.atom_selection[i]->PutUDData(uddHnd, old_atom_index);
	       }
	    }
	 }
      }
   }
}
