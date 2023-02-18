/* coords/mmdb.cc
 * 
 * Copyright 2005 by The University of York
 * Copyright 2009 by the University of Oxford
 * Copyright 2013, 2015 by Medical Research Council
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

#include <iostream>
#include <string>
#include <vector> // for Cartesian.h

#include "string.h"


#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include "lidia-core/rdkit-interface.hh"
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/FileParsers/FileParsers.h>
#endif

#include <mmdb2/mmdb_manager.h>
#include "mmdb-extras.h"
#include "Cartesian.h"
#include "mmdb.hh"
#include "coot-utils/lidia-core-functions.hh"

#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-shelx.hh"

#include "coot-utils/read-sm-cif.hh"




// Dump this.  Use instead GetCentreOfMass(selhand)
//
coot::Cartesian
centre_of_molecule(atom_selection_container_t SelAtom) {

   coot::Cartesian centre; // defaults construction at (0,0,0).
   coot::Cartesian rs;     // running sum
   mmdb::PAtom atom;

   if (SelAtom.n_selected_atoms > 0) { 
      for (int i=0; i< SelAtom.n_selected_atoms; i++) {
	 
	 atom = SelAtom.atom_selection[i];
	 coot::Cartesian atom_pos (atom->x, atom->y, atom->z);
	 
	 rs += atom_pos;
      }
      
      float scale = 1/float (SelAtom.n_selected_atoms);
      centre = coot::Cartesian(rs.get_x()*scale, rs.get_y()*scale, rs.get_z()*scale);
   }
   
   return centre;
   
}


// should be a const reference in an ideal world.
//
std::ostream& operator<<(std::ostream& s, mmdb::Atom &atom) {

   //
   s << atom.GetModelNum() << "/" << atom.GetChainID() << "/"
     << atom.GetSeqNum() << atom.GetInsCode() << "/" << atom.GetResName()
     << "/" << atom.name << " altLoc :" << atom.altLoc << ": pos: ("
     << atom.x << "," << atom.y << "," << atom.z
     << ") B-factor: " << atom.tempFactor;

   return s;

}
  
// needs <iostream.h>
// 
std::ostream& operator<<(std::ostream& s, mmdb::PAtom atom) {

   //
   if (atom) {
      s << atom->GetModelNum() << "/" << atom->GetChainID() << "/"
	<< atom->GetSeqNum()   << atom->GetInsCode() << " {"
	<< atom->GetResName() << "}/"
	<< atom->name << " altLoc :" << atom->altLoc << ": segid :"
	<< atom->segID << ":" << " pos: ("
	<< atom->x << "," << atom->y << "," << atom->z
	<< ") B-factor: " << atom->tempFactor;
   } else {
      s << "NULL";
   } 

   return s;

}
  

int
write_atom_selection_file(atom_selection_container_t asc,
			  const std::string &filename,
			  bool write_as_cif_flag,
			  mmdb::byte gz,
			  bool write_hydrogens,     // optional arg
			  bool write_aniso_records, // optional arg
			  bool write_conect_records // optional arg
			  ) {

   int ierr = 0; 
   coot::util::remove_wrong_cis_peptides(asc.mol);
   mmdb::Manager *mol = asc.mol;
   bool mol_needs_deleting = false; // unless mol is reassigned...

   if (write_as_cif_flag) {

      // WriteCIFASCII() seems to duplicate the atoms (maybe related to aniso?)
      // So let's copy the molecule and throw away the copy, that way we don't
      // duplicate teh atoms in the original molecule.

      mmdb::Manager *mol_copy  = new mmdb::Manager;
      mol_copy->Copy(mol, mmdb::MMDBFCM_All);
      ierr = mol_copy->WriteCIFASCII(filename.c_str());
      delete mol_copy;

   } else {

      if (! write_hydrogens) {
	 mmdb::Manager *n = new mmdb::Manager;
	 n->Copy(mol, mmdb::MMDBFCM_All);
	 coot::delete_hydrogens_from_mol(n);
	 mol = n;
	 mol_needs_deleting = true;
      }

      if (! write_aniso_records) {
	 mmdb::Manager *n = new mmdb::Manager;
	 n->Copy(mol, mmdb::MMDBFCM_All);
	 coot::delete_aniso_records_from_atoms(n);
	 mol = n;
	 mol_needs_deleting = true;
      }

      if (! write_conect_records) {
	 mmdb::Manager *n = new mmdb::Manager;
	 n->Copy(mol, mmdb::MMDBFCM_All);
	 // Eugene's magic code
	 n->Delete ( mmdb::MMDBFCM_SC );
	 mol = n;
	 mol_needs_deleting = true;
      }

      coot::util::remove_long_links(mol, 2.1);
      
      // we need to put the hydrogen names back to how they used to be
      // when we read in the pdb file and then put them put them back
      // to how they currently are!

      int udd_old = mol->GetUDDHandle(mmdb::UDR_ATOM, "initial hydrogen name");
      int udd_new = mol->GetUDDHandle(mmdb::UDR_ATOM, "new hydrogen name");
//       std::cout << "udd_old: " << udd_old << std::endl;
//       std::cout << "udd_new: " << udd_new << std::endl;
      char *str;
      if (udd_old > 0 && udd_new > 0) { 
	 for (int i=0; i<asc.n_selected_atoms; i++) {
	    str = 0; 
	    if (asc.atom_selection[i]->GetUDData(udd_old, str) == mmdb::UDDATA_Ok) {
// 	       std::cout << "pre  UDD: " << asc.atom_selection[i]
// 			  << " gave udd str: " << str << std::endl;
	       asc.atom_selection[i]->SetAtomName(str);
	    } else {
	       // not a hydrogen, probably
	    }
	 }
      }
      ierr = mol->WritePDBASCII(filename.c_str());
      // now put the names back
      if (udd_old > 0 && udd_new > 0) { 
      	 for (int i=0; i<asc.n_selected_atoms; i++) {
	    str = 0; 
      	    if (asc.atom_selection[i]->GetUDData(udd_new, str) == mmdb::UDDATA_Ok) { 
       	       asc.atom_selection[i]->SetAtomName(str);
//        	       std::cout << "post UDD: " << asc.atom_selection[i]
//        			 << " gave udd str: " << str << std::endl;
	    } else {
	       // not a hydrogen, probably
       	    }
       	 }
      }
   }

   if (mol_needs_deleting)
      delete mol;
   
   return ierr; 
}


void
coot::delete_hydrogens_from_mol(mmdb::Manager *mol) {

   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      mmdb::Chain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 mmdb::Residue *residue_p;
	 mmdb::Atom *at;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    int n_atoms = residue_p->GetNumberOfAtoms();
	    bool deleted = 0;
	    for (int iat=0; iat<n_atoms; iat++) {
	       at = residue_p->GetAtom(iat);
	       std::string ele(at->element);
	       if (is_hydrogen(ele)) {
		  // delete this atom
		  deleted = 1;
		  residue_p->DeleteAtom(iat);
	       }
	    }
	    if (deleted)
	       residue_p->TrimAtomTable();
	 }
      }
   }
}


void
coot::delete_aniso_records_from_atoms(mmdb::Manager *mol) {

   std::cout << "ASET_Anis_tFac " << mmdb::ASET_Anis_tFac << " " << ~mmdb::ASET_Anis_tFac << std::endl;
   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      mmdb::Chain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 mmdb::Residue *residue_p;
	 mmdb::Atom *at;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    int n_atoms = residue_p->GetNumberOfAtoms();
	    for (int iat=0; iat<n_atoms; iat++) {
	       at = residue_p->GetAtom(iat);
	       at->WhatIsSet &=  ~mmdb::ASET_Anis_tFac;
	    }
	 }
      }
   }
}


// This is here because it's an mmdb function. There is nothing
// coordinates related here. Perhaps this should be a coot-utils
// function?

std::vector<std::string>
coot::get_compound_lines(mmdb::Manager *mol) {

   std::vector<std::string> compound_lines;

   access_mol *am = static_cast<access_mol *>(mol); // causes indent problem

   const mmdb::Title *tt = am->GetTitle();
   mmdb::Title *ttmp = const_cast<mmdb::Title *>(tt);
   access_title *at = static_cast<access_title *> (ttmp);
   mmdb::TitleContainer *compound = at->GetCompound();
   int cl = compound->Length();

   if (cl > 0) {
      for(int i=0; i<cl; i++)  {
	 mmdb::Compound *CLine = mmdb::PCompound(compound->GetContainerClass(i));
	 if (CLine) {
	    std::string line(CLine->Line);
	    compound_lines.push_back(line);
	 }
      }
   }
   return compound_lines;
}


// This is here because it's an mmdb function. There is nothing
// coordinates related here. Perhaps this should be a coot-utils
// function?
std::string coot::get_title(mmdb::Manager *mol) {

   std::string tt;

   char *title = new char[10240];

   char *t = mol->GetStructureTitle(title);

   if (t) {
      tt = std::string(t);
   }
   delete [] title;

   return tt;
}
