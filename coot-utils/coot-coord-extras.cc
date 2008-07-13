/* coot-utils/coot-coord-extras.cc
 * 
 * Copyright 2004, 2005, 2006, 2007 by The University of York
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

// #include <algorithm>
#include "coot-utils.hh"
#include "coot-coord-utils.hh"
#include "coot-coord-extras.hh"


// Return 0 if any of the residues don't have a dictionary entry
// geom_p gets updated to include the residue restraints if necessary
// 
std::pair<int, std::vector<std::string> >
coot::util::check_dictionary_for_residues(PCResidue *SelResidues, int nSelResidues,
					  coot::protein_geometry *geom_p,
					  int read_number) {

   std::pair<int, std::vector<std::string> > r;

   int status;
   int fail = 0; // not fail initially.

   for (int ires=0; ires<nSelResidues; ires++) { 
      std::string resname(SelResidues[ires]->name);
      status = geom_p->have_dictionary_for_residue_type(resname, read_number);
      // This bit is redundant now that try_dynamic_add has been added
      // to have_dictionary_for_residue_type():
      if (status == 0) { 
	 status = geom_p->try_dynamic_add(resname, read_number);
	 if (status == 0) {
	    fail = 1; // we failed to find it then.
	    r.second.push_back(resname);
	 }
      }
   }

   if (fail)
      r.first = 0;
   return r;
} 





// We also now pass regular_residue_flag so that the indexing of the
// contacts is inverted in the case of not regular residue.  I don't
// know why this is necessary, but I have stared at it for hours, this
// is a quick (ugly hack) fix that works.  I suspect that there is
// some atom order dependency in mgtree that I don't understand.
// Please fix (remove the necessity of depending on
// regular_residue_flag) if you know how.
// 
std::vector<std::vector<int> >
coot::util::get_contact_indices_from_restraints(CResidue *residue,
						coot::protein_geometry *geom_p,
						short int regular_residue_flag) {

   int nResidueAtoms = residue->GetNumberOfAtoms(); 
   std::vector<std::vector<int> > contact_indices(nResidueAtoms);
   std::string restype(residue->name);
   CAtom *atom_p;

   int n_restr = geom_p->size();

   for (int icomp=0; icomp<n_restr; icomp++) {
      if ((*geom_p)[icomp].comp_id == restype) {
// 	 std::cout << "There are " << (*geom_p)[icomp].bond_restraint.size()
// 		   << " bond restraints " << "for " << restype << std::endl;
	 for (unsigned int ibr=0; ibr< (*geom_p)[icomp].bond_restraint.size(); ibr++) {
	    for (int iat=0; iat<nResidueAtoms; iat++) {
	       atom_p = residue->GetAtom(iat);
	       std::string at_name(atom_p->GetAtomName());
	       if ( (*geom_p)[icomp].bond_restraint[ibr].atom_id_1_4c() == at_name ) {
// 		  std::cout << "found a bond match "
// 			    << (*geom_p)[icomp].bond_restraint[ibr].atom_id_1_4c()
// 			    << " to "
// 			    << (*geom_p)[icomp].bond_restraint[ibr].atom_id_2_4c()
// 			    << std::endl;
		  int ibond_to = -1;  // initially unassigned.
		  std::string at_name_2;
		  for (int iat2=0; iat2<nResidueAtoms; iat2++) {
		     atom_p = residue->GetAtom(iat2);
		     at_name_2 = atom_p->GetAtomName();
		     if ( (*geom_p)[icomp].bond_restraint[ibr].atom_id_2_4c() == at_name_2 ) {
			ibond_to = iat2;
			break;
		     }
		  }
		  if (ibond_to > -1 ) { 
		     if (regular_residue_flag) {
			contact_indices[iat].push_back(ibond_to);  // for ALA etc
		     } else {
			contact_indices[ibond_to].push_back(iat);  // ligands
			// contact_indices[iat].push_back(ibond_to);  // ALA etc
		     }
		  } 
//		  else
		     // This spits out the names of Hydrogens, often.
//  		     std::cout << "failed to find bonded atom "
//  			       << (*geom_p)[icomp].bond_restraint[ibr].atom_id_2_4c()
//  			       << std::endl;
	       }
	    }
	 }
      }
   }
   return contact_indices;
}

// The atoms of residue_atoms are in the "right" order for not making
// a tree along the main chain.
// 
std::vector<std::vector<int> >
coot::util::get_contact_indices_for_PRO_residue(PPCAtom residue_atoms,
						int nResidueAtoms, 
						coot::protein_geometry *geom_p) { 

   std::vector<std::vector<int> > contact_indices(nResidueAtoms);
   CAtom *atom_p;
   int n_restr = geom_p->size();
   for (int icomp=0; icomp<n_restr; icomp++) {
      if ((*geom_p)[icomp].comp_id == "PRO") {
	 for (unsigned int ibr=0; ibr< (*geom_p)[icomp].bond_restraint.size(); ibr++) {
	    for (int iat=0; iat<nResidueAtoms; iat++) {
	       atom_p = residue_atoms[iat];
	       std::string at_name(atom_p->GetAtomName());
	       if ( (*geom_p)[icomp].bond_restraint[ibr].atom_id_1_4c() == at_name ) {
		  int ibond_to = -1;  // initially unassigned.
		  std::string at_name_2;
		  for (int iat2=0; iat2<nResidueAtoms; iat2++) {
		     atom_p = residue_atoms[iat2];
		     at_name_2 = atom_p->GetAtomName();
		     if ( (*geom_p)[icomp].bond_restraint[ibr].atom_id_2_4c() == at_name_2 ) {
			ibond_to = iat2;
			break;
		     }
		  }
		  if (ibond_to != -1)
		     contact_indices[iat].push_back(ibond_to);
	       }
	    }
	 }
      }
   }
   return contact_indices;
}


coot::util::dict_residue_atom_info_t::dict_residue_atom_info_t(const std::string &residue_name_in,
							       coot::protein_geometry *geom_p) {

   residue_name = residue_name_in;

   std::pair<short int, dictionary_residue_restraints_t> p = 
      geom_p->get_monomer_restraints(residue_name);

   if (p.first) {
      for (unsigned int iat=0; iat<p.second.atom_info.size(); iat++) {
	 std::string atom_name = p.second.atom_info[iat].atom_id_4c;
	 short int isHydrogen = 0;
	 if (p.second.atom_info[iat].type_symbol == "H" ||
	     p.second.atom_info[iat].type_symbol == "D") {
	    isHydrogen = 1;
	 }
	 atom_info.push_back(coot::util::dict_atom_info_t(atom_name, isHydrogen));
      }
   }

}

// This one we can do a dynamic add.
// 
short int
coot::util::is_nucleotide_by_dict_dynamic_add(CResidue *residue_p, coot::protein_geometry *geom_p) {

   short int is_nuc = 0;
   short int ifound = 0;
   std::string residue_name = residue_p->GetResName();

   int n_restr = geom_p->size();
   for (int icomp=0; icomp<n_restr; icomp++) {
      if ((*geom_p)[icomp].comp_id == residue_name) {
	 ifound = 1;
	 if ((*geom_p)[icomp].residue_info.group == "RNA" ||
	     (*geom_p)[icomp].residue_info.group == "DNA" ) {
	    is_nuc = 1;
	 }
	 break;
      }
   }

   int read_number = 40;
   if (ifound == 0) {
      int status = geom_p->try_dynamic_add(residue_name, read_number);
      if (status != 0) {
	 // we successfully added it, let's try to run this function
	 // again.  Or we could just test the last entry in
	 // geom_p->dict_res_restraints(), but it is not public, so
	 // it's messy.
	 // 
	 is_nuc = is_nucleotide_by_dict_dynamic_add(residue_p, geom_p);
      } 
   } 
   return is_nuc;
}


// This one we can NOT do a dynamic add.
//
short int
coot::util::is_nucleotide_by_dict(CResidue *residue_p, const coot::protein_geometry &geom) {

   short int is_nuc = 0;
   std::string residue_name = residue_p->GetResName();

   int n_restr = geom.size();
   for (int icomp=0; icomp<n_restr; icomp++) {
      if (geom[icomp].comp_id == residue_name) {
	 if (geom[icomp].residue_info.group == "RNA" ||
	     geom[icomp].residue_info.group == "DNA" ) {
	    is_nuc = 1;
	 }
	 break;
      }
   }

   return is_nuc;
}

