/* geometry/protein-geometry.cc
 * 
 * Copyright 2010 The University of Oxford
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
#include <stdlib.h>

#include "coot-utils.hh"
#include "protein-geometry.hh"

void
coot::protein_geometry::read_sbase_residues() {

   if (SBase) { 
      PCSBStructure SBS;
      std::vector<std::string> local_residue_codes;
      local_residue_codes.push_back("ASN");
      local_residue_codes.push_back("TRP");

      int nStructures = SBase->GetNofStructures();
      printf ( "  Total %i structures in SBASE\n",nStructures );
      std::cout << "INFO:: Total of " << nStructures << " in SBase\n";
      
      for (int i=0; i<local_residue_codes.size(); i++) {
	 SBS = SBase->GetStructure(local_residue_codes[i].c_str());
	 if (SBS) { 
	    std::cout << SBS->compoundID << " " << SBS->Name << std::endl;
	    for (int iat=0; iat<SBS->nAtoms; iat++) {
	       CSBAtom *at = SBS->Atom[iat];
	       std::cout << "    " << at->sca_name << " ("
			 << at->x << "," << at->y << ","
			 << at->z <<")\n";
	    } 
	 } else {
	    std::cout << "WARNING:: structure " << local_residue_codes[i]
		      << " not found in SBase " << std::endl;
	 } 
      }
   } else {
      std::cout << "WARNING:: SBase not initialised"  << std::endl;
   } 
}

// Here res_name is the tlc/comp_id.
// 
CResidue *
coot::protein_geometry::get_sbase_residue(const std::string &res_name) const {

   CResidue *residue_p = NULL;
   std::cout << "SBASE: " << SBase << std::endl;
   if (SBase) {
      CSBStructure *SBS = SBase->GetStructure(res_name.c_str());
      if (SBS) {
	 residue_p = new CResidue;
	 for (int iat=0; iat<SBS->nAtoms; iat++) {
	    CSBAtom *at = SBS->Atom[iat];

	    std::string new_atom_name = coot::atom_id_mmdb_expand(at->sca_name, at->element);

	    // we need to add a sanity check on the position of the
	    // atoms (e.g. S in SO4 is at (-1.7976e+308 -1.7976e+308
	    // -1.7976e+308)).

	    bool add_atom_ok = 0;
	    if (fabs(at->x) < 2000) { 
	       if (fabs(at->y) < 2000) { 
		  if (fabs(at->z) < 2000) {
		     add_atom_ok = 1;
		  }
	       }
	    }

	    if (add_atom_ok) { 
	       CAtom *new_atom = new CAtom;
	       new_atom->SetCoordinates(at->x, at->y, at->z, 1.0, 30.0);
	       new_atom->SetAtomName(new_atom_name.c_str());
	       new_atom->SetElementName(at->element);
	       new_atom->Het = true; 
	       residue_p->AddAtom(new_atom);
	    } else {
	       std::cout << "WARNING:: rejecting " << res_name << " SBase atom :" << new_atom_name
			 << ": at (" << at->x << " " << at->y << " " << at->z << ")" << std::endl;
	    }
	 }
      }
   }
   std::cout << "get_sbase_residue() returns " << residue_p << std::endl;
   return residue_p;
}


std::vector<std::pair<std::string, std::string> >
coot::protein_geometry::matching_sbase_residues_names(const std::string &compound_name_frag) const {

   std::string compound_name = coot::util::upcase(compound_name_frag);
   std::vector<std::pair<std::string,std::string> > v;
   if (SBase) {
      int nStructures = SBase->GetNofStructures();
      CFile *sf = SBase->GetStructFile();
      for (int is=0; is<nStructures; is++) {
	 CSBStructure *SBS = SBase->GetStructure(is, sf);
	 if (SBS) {
	    if (SBS->Name) { 
	       std::string sbase_residue_name(SBS->Name);
	       std::string::size_type imatch = sbase_residue_name.find(compound_name);
	       if (imatch != std::string::npos) {
		  if (0) 
		     std::cout << "DEBUG:: " << imatch << " :" << compound_name
			       << ": found in " << sbase_residue_name << std::endl;
		  std::pair<std::string,std::string> p(SBS->compoundID, SBS->Name);
		  v.push_back(p);
	       }
	    }
	 }
      }
   }
   return v;
}


// return mmdb sbase return codes
//
// Try to use the MONOMER_DIR_STR, ie. COOT_SBASE_DIR first, if that
// fails then use the fallback directory sbase_monomer_dir_in
// 
int
coot::protein_geometry::init_sbase(const std::string &sbase_monomer_dir_in) {

   // in the protein_geometry constructor SBase is set to NULL.
   
   int RC = SBASE_FileNotFound; // initial status.
   const char *monomer_dir = getenv(MONOMER_DIR_STR);
   
   if (!monomer_dir) {
      if (coot::is_directory_p(sbase_monomer_dir_in)) {
	 monomer_dir = sbase_monomer_dir_in.c_str();
      } else { 
	 RC = SBASE_FileNotFound; // fail
      }
   }

   if (monomer_dir) { 

      std::cout << "sbase monomer dir: " << monomer_dir << std::endl;
      SBase = new CSBase;
      RC = SBase->LoadIndex(monomer_dir);

      if (RC != SBASE_Ok) {
         std::cout << "sbase files not found in " << monomer_dir << std::endl;
	 delete SBase;
	 SBase = NULL;
      } else { 
         // std::cout << "sbase files found" << std::endl; 
      }
   }
   return RC; 
}


// used to created data from sbase to put into protein_geometry
// object.
//
// return a flag to signify success.
//
// This relies on SBase being setup before we get to make this
// call.
// 
bool
coot::protein_geometry::fill_using_sbase(const std::string &monomer_type) {

   bool success = 0;
   coot::dictionary_residue_restraints_t rest(true); // constructor for SBase
   rest.residue_info.comp_id = monomer_type;
   
   if (SBase) {
      CSBStructure *SBS = SBase->GetStructure(monomer_type.c_str());
      if (SBS) {
	 for (int iat=0; iat<SBS->nAtoms; iat++) {
	    CSBAtom *at = SBS->Atom[iat];
	    std::pair<bool, float> pc(0,0); // partial charge.
	    std::string atom_id = at->pdb_name;
	    std::string atom_id_4c = atom_id_mmdb_expand(atom_id);
	    std::string type_symbol = at->element;
	    std::string type_energy = at->energyType;
	    coot::dict_atom dict_at(atom_id,
				    atom_id_4c,
				    type_symbol,
				    type_energy,
				    pc);
	    rest.atom_info.push_back(dict_at);
	 }
      }

      if (rest.atom_info.size() > 1) {
	 for (int ib=0; ib<SBS->nBonds; ib++) {
	    CSBBond *bond = SBS->Bond[ib];
	    int ind_1 = bond->atom1 -1 ; // Yikes!
	    int ind_2 = bond->atom2 -1 ;
	    int order = bond->order;
	    double dist = bond->length;
	    double esd = bond->length_esd;
	    if (0) 
	       std::cout << "... "<< monomer_type << " atom index " << ind_1 << " " << ind_2
			 << " order " << order << std::endl;
	    std::string type = "single";
	    if (order == 2)
	       type = "double";
	    if (order == 3)
	       type = "triple";
	    std::string atom_name_1 = SBS->Atom[ind_1]->pdb_name;
	    std::string atom_name_2 = SBS->Atom[ind_2]->pdb_name;
	    coot::dict_bond_restraint_t dict_bond(atom_name_1, atom_name_2, type, dist, esd);
	    success = 1;
	    rest.bond_restraint.push_back(dict_bond);
	 }
      }
   }
   if (success)
      add(rest);
   
   return success;
}


// A new pdb file has been read in (say).  The residue types
// have been compared to the dictionary.  These (comp_ids) are
// the types that are not in the dictionary.  Try to load an
// sbase description at least so that we can draw their bonds
// correctly.  Use fill_using_sbase().
// 
bool
coot::protein_geometry::try_load_sbase_description(const std::vector<std::string> &comp_ids_with_duplicates) {

   bool status = 0; // none added initially.

   // If this is INH, DRG etc, don't try to auto-add
   // 
   if (is_non_auto_load_ligand(resname)) {
      std::cout << "INFO:: comp-id: " << resname << " is marked for non-autoloading - stopping now "
		<< std::endl;
      return false;
   }
   

   std::vector<std::string> uniques;
   for (unsigned int ic=0; ic<comp_ids_with_duplicates.size(); ic++) { 
      if (std::find(uniques.begin(), uniques.end(), comp_ids_with_duplicates[ic]) ==
	  uniques.end())
	 uniques.push_back(comp_ids_with_duplicates[ic]);
	  
   }

   
   if (SBase) {
      for (unsigned int i=0; i<uniques.size(); i++) {
	 std::cout << i << " " << uniques[i] << std::endl;
	 const std::string comp_id = uniques[i];
	 if (is_non_auto_load_ligand(comp_id)) {
	    std::cout << "INFO:: sbase-descriptions: comp-id: " << comp_id
		      << " is marked for non-autoloading - ignore " << std::endl;
	 } else { 
	    bool s = fill_using_sbase(comp_id);
	    if (s)
	       std::cout << "DEBUG:: sbase dictionary for " << comp_id << " successfully loaded "
			 << std::endl;
	    if (s)
	       status = 1;
	 }
      }
   }
   return status; // return true of something was added.
}
