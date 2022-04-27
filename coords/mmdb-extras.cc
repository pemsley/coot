/* coords/mmdb-extras.cc
 * 
 * Copyright 2006 by The University of York
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

//#include <stdlib.h>  // for getenv()
//#include <string.h>  // for memcpy()
#include <string>
#include <vector> // for Cartesian
//#include <map>
//#include <algorithm>

#include "compat/coot-sysdep.h"

#ifdef COOT_ENABLE_WINAPI_SUSPENSION
# undef AddAtom
#endif // COOT_ENABLE_WINAPI_SUSPENSION

#include "coot-utils/coot-coord-utils.hh"

#include "mmdb-extras.h" 
#include "mmdb.h"

#include <cstring>

float max_bond_length(const std::string &element) { 

   if (element == " S") { 
      return 2.3;
   } else { 
      return 1.8;
   }
};

int
get_atom_colour_from_element(const std::string &element) {

   if (element == " C") {
      return YELLOW_BOND;
   } else {
      if (element == " N") {
	 return BLUE_BOND;
      } else {
	 if (element == " O") {
	    return RED_BOND;
	 } else {
	    if (element == " S") {
	       return GREEN_BOND;
	    } else {
	       if (coot::is_hydrogen(element)) {
                  if (coot::is_deuterium(element))
                     return DEUTERIUM_PINK;
                  else
                     return HYDROGEN_GREY_BOND;
	       }
	    }
	 }
      }
   }
   return GREY_BOND;
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



// Note, we also create a chain and add this residue to that chain.
// We do this so that we have a holder for the segid.
// 
// whole_residue_flag: only copy atoms that are either in this altLoc,
// or has an altLoc of "".
//
// This always returns a mmdb::Residue pointer pointing to a allocated
// residue.  It may not have any atoms in it though.
// 
mmdb::Residue *
coot::deep_copy_this_residue_old_style(mmdb::Residue *residue,
			     const std::string &altconf,
			     short int whole_residue_flag,
			     int atom_index_handle,
			     bool embed_in_chain_flag) {

   // 20090622 altconf is "", whole_residue_flag = 0, and residue is
   // completely split into A and B.
   //
   // This would return a residue that is valid, but has no atoms in
   // it.  However when we come to do a GetAtomTable() on such as
   // residue (e.g. init_from_residue_vec in simple_restraint.cc) we
   // get a crash because for some reason I cannot understand, nAtoms
   // is not 0 and we try to read atoms in the residue_atoms that do
   // not exist...  So now, when there are no atoms added to the
   // returned residue, I will delete the residue and return a NULL.
   // I have changed all uses of this function to now check for NULL.
   //
   // 

   mmdb::Residue *rres = 0; 
   mmdb::PPAtom residue_atoms;
   int nResidueAtoms;
   residue->GetAtomTable(residue_atoms, nResidueAtoms);
   mmdb::Atom *atom_p;

   if (nResidueAtoms > 0) { 
      // 
      mmdb::Chain   *chain_p = NULL;
      rres = new mmdb::Residue();
      rres->SetResID(residue->GetResName(),
		     residue->GetSeqNum(),
		     residue->GetInsCode());
      
      int n_added_atoms = 0; 
      for(int iat=0; iat<nResidueAtoms; iat++) {
	 mmdb::Atom *at = residue_atoms[iat];
	 if (at) { // can have been deleted
	    if (! at->isTer()) { 
	       std::string this_atom_alt_loc(residue_atoms[iat]->altLoc);
	       if (whole_residue_flag ||
		   this_atom_alt_loc  == altconf || this_atom_alt_loc == "") { 
		  atom_p = new mmdb::Atom;
		  atom_p->Copy(residue_atoms[iat]);
		  int i_add = rres->AddAtom(atom_p);
		  n_added_atoms++;
	       }
	    }
	 }
      }
      if (n_added_atoms == 0) {
	 // reset the returned residue
	 delete rres;
	 rres = NULL;
      } else {
	 // As normal
	 
	 if (embed_in_chain_flag) { 
	    chain_p = new mmdb::Chain;
	    chain_p->SetChainID(residue->GetChainID());
	    chain_p->AddResidue(rres);
	 }

	 // debug
	 if (0) { 
	    std::cout << "debug:: coot::deep_copy_this_residue returns a residue copy of "
		      << residue ->GetChainID() << " " << residue->GetSeqNum() <<  " "
		      << residue->GetInsCode() << " " << "with alt confs \"" << altconf << "\""
		      << " which is " << rres <<  " with "  
		      << rres->nAtoms << " atoms " << " and atoms array " << rres->atom
		      << std::endl;
	 }
      }
   } else {
      if (0) {  // debug
	 std::cout << "Debug:: deep_copy_this_residue returns NULL" << std::endl;
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
std::pair<mmdb::Residue *, atom_selection_container_t>
coot::deep_copy_this_residue_and_make_asc(mmdb::Manager *orig_mol,
					  mmdb::Residue *residue,
					  const std::string &altconf,
					  short int whole_residue_flag,
					  int atom_index_handle,
					  int udd_afix_handle) {

//    std::cout << "DEbbug:: in deep_copy_this_residue_and_make_asc is "
// 	     << udd_afix_handle << std::endl;

   // Horrible casting to mmdb::Residue because GetSeqNum and GetAtomTable
   // are not const functions.
   // 
   mmdb::Residue *rres = new mmdb::Residue;
   mmdb::Chain   *chain_p = new mmdb::Chain;
   std::string chain_id1 = ((mmdb::Residue *)residue)->GetChainID();
   chain_p->SetChainID(chain_id1.c_str());
   rres->seqNum = ((mmdb::Residue *)residue)->GetSeqNum();
   /* Copy insertion code - a char[10] */
   std::memcpy(rres->insCode, residue->insCode, sizeof(mmdb::InsCode));
   std::memcpy(rres->name, residue->name, sizeof(mmdb::ResName));
   
   mmdb::PPAtom residue_atoms;
   int nResidueAtoms;
   residue->GetAtomTable(residue_atoms, nResidueAtoms);
   mmdb::Atom *atom_p;
   short int do_shelx_afix_data_flag = 0;
   if (udd_afix_handle >= 0)
      do_shelx_afix_data_flag = 1;
   
   for(int iat=0; iat<nResidueAtoms; iat++) {
      if (! residue_atoms[iat]->isTer()) { 
	 std::string this_atom_alt_loc(residue_atoms[iat]->altLoc);
	 if (whole_residue_flag ||
	     this_atom_alt_loc  == altconf || this_atom_alt_loc == "") { 
	    atom_p = new mmdb::Atom;
	    atom_p->Copy(residue_atoms[iat]);
	    int i_add = rres->AddAtom(atom_p);
	    if (atom_index_handle > -1 ) { 
	       // copy across the atom indices:
	    
	    }
// 	 if (do_shelx_afix_data_flag) {
// 	    int ic;
// 	    if (residue_atoms[iat]->GetUDData(udd_afix_handle, ic) == mmdb::UDDATA_Ok) {
// 	       std::cout << "Copying afix handle " << ic << " from " << residue_atoms[iat]
// 			 << "\n                        to " << atom_p << std::endl;
// 	       int istat = atom_p->PutUDData(udd_afix_handle, ic);
// 	       if (istat != mmdb::UDDATA_Ok) { 
// 		  std::cout << "Opps! putuddata returns status " << istat << " vs "
// 			    << mmdb::UDDATA_Ok << std::endl;
// 	       }
// 	    }
// 	 }
	 }
      }
   }
   chain_p->AddResidue(rres);

   // urgh.  This consting business makes a mess.
   mmdb::Residue *r = (mmdb::Residue *)residue;
   mmdb::PResidue *SelResidues = &r;
   std::pair<mmdb::Manager *, int> mol_i =
      coot::util::create_mmdbmanager_from_res_selection(orig_mol, SelResidues,
							1, 0, 0, altconf, chain_id1, 1);

   atom_selection_container_t asc = make_asc(mol_i.first);
   int udd_afix_handle_inter = mol_i.first->GetUDDHandle(mmdb::UDR_ATOM, "shelx afix");
//    std::cout << "DEBUG:: deep_copy_this_residue_and_make_asc got udd_afix_handle_inter : "
// 	     << udd_afix_handle_inter << std::endl;
   for (int i=0; i<asc.n_selected_atoms; i++) {
      int afix_number = -1;
      if (asc.atom_selection[i]->GetUDData(udd_afix_handle_inter, afix_number) == mmdb::UDDATA_Ok)
	 std::cout << asc.atom_selection[i] << " has afix number " << afix_number
		   << std::endl;
//       else
// 	 std::cout << asc.atom_selection[i]
// 		   << " Failed get udd afix number in deep_copy_this_residue_and_make_asc"
// 		   << std::endl;
   }
   return std::pair<mmdb::Residue *, atom_selection_container_t>(rres, asc);
}


short int
coot::progressive_residues_in_chain_check(const mmdb::Chain *chain_p) {

   int nres = ((mmdb::Chain*)chain_p)->GetNumberOfResidues();
   mmdb::Residue *res_p;
   int previous_seq_num = -9999;
   int this_seq_no;

   for (int ires=0; ires<nres; ires++) { 
      res_p = ((mmdb::Chain*)chain_p)->GetResidue(ires);
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


// Typically this is used on an asc (moving atoms) to get the N of a
// peptide (say).  Return NULL on atom not found.
// 
mmdb::Atom *
coot::get_first_atom_with_atom_name(const std::string &atomname, 
				    const atom_selection_container_t &asc) { 

   mmdb::Atom *atom = NULL;

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
	 int uddHnd = asc.mol->RegisterUDInteger(mmdb::UDR_ATOM , "old atom index");
	 if (uddHnd >= 0) { 
	    asc.UDDOldAtomIndexHandle = uddHnd;
	    for (int i=0; i<asc.n_selected_atoms; i++) { 
	       if (asc.atom_selection[i]->GetUDData(asc.UDDAtomIndexHandle, old_atom_index) == mmdb::UDDATA_Ok) { 
		  asc.atom_selection[i]->PutUDData(uddHnd, old_atom_index);
	       }
	    }
	 }
      }
   }
}
