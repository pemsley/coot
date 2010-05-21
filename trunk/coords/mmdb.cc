/* coords/mmdb.cc
 * 
 * Copyright 2005 by The University of York
 * Copyright 2009 by the University of Oxford
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

#include "mmdb_manager.h"

#include "mmdb-extras.h"

#include "Cartesian.h"

#include "mmdb.h"

#include "coot-utils.hh"
#include "coot-coord-utils.hh"
#include "coot-shelx.hh"


// This is used for pick_test  (a function that returns
// an atom selection from a pdb_file name string is not
// generally so useful).
//
atom_selection_container_t
get_atom_selection(std::string pdb_name) {

   int err;
   MyCMMDBManager* MMDBManager;

   // Needed for the error message printing: 
   // MMDBManager->GetInputBuffer(S, lcount);
   // Used by reference and as a pointer.  Grimness indeed.
   int  error_count;
   char error_buf[500];

   //   Make routine initializations
   //
   InitMatType();

   atom_selection_container_t asc;

    std::string extention = coot::util::file_name_extension(pdb_name);

    // returns e.g. ".ins"


    if (coot::util::extension_is_for_shelx_coords(extention)) {


       coot::ShelxIns s;
       coot::shelx_read_file_info_t srf = s.read_file(pdb_name);
       // atom_selection_container_t.mol is of type MyCMMDBManager *
       // currently.
       asc = make_asc(srf.mol);
       MMDBManager = asc.mol;
       if (asc.mol)
	  asc.read_success = 1;  // a good idea?

    } else { 
   
       MMDBManager = new MyCMMDBManager();

       // For mmdb version 1.0.3:
       //    MMDBManager->SetFlag ( MMDBF_IgnoreBlankLines |
       // 			  MMDBF_IgnoreDuplSeqNum |
       // 			  MMDBF_IgnoreNonCoorPDBErrors);
       //
       // From mmdb versions 1.0.4 to 1.0.7:
       // 
       //    MMDBManager->SetFlag ( MMDBF_IgnoreBlankLines |
       // 			      MMDBF_IgnoreDuplSeqNum |
       // 			      MMDBF_IgnoreNonCoorPDBErrors |
       // 			      MMDBF_IgnoreRemarks);
       // 
       // For mmdb version 1.0.8 and beyond:

#ifdef HAVE_MMDB_IGNORE_HASH
       
       MMDBManager->SetFlag ( MMDBF_IgnoreBlankLines |
			      MMDBF_IgnoreDuplSeqNum |
			      MMDBF_IgnoreNonCoorPDBErrors |
			      MMDBF_IgnoreHash |
			      MMDBF_IgnoreRemarks);
#else
       MMDBManager->SetFlag ( MMDBF_IgnoreBlankLines |
			      MMDBF_IgnoreDuplSeqNum |
			      MMDBF_IgnoreNonCoorPDBErrors |
			      MMDBF_IgnoreRemarks);
#endif // HAVE_MMDB_IGNORE_HASH       
       
       std::cout << "Reading coordinate file: " << pdb_name.c_str() << "\n";
       err = MMDBManager->ReadCoorFile(pdb_name.c_str());
   
       if (err) {
	  // does_file_exist(pdb_name.c_str());
	  cout << "There was an error reading " << pdb_name.c_str() << ". \n";
	  cout << "ERROR " << err << " READ: "
	       << GetErrorDescription(err) << endl;
	  //
	  // This makes my stomach churn too. Sorry.
	  // 
	  MMDBManager->GetInputBuffer(error_buf, error_count);
	  if (error_count >= 0) { 
	     cout << "         LINE #" << error_count << "\n     "
		  << error_buf << endl << endl;
	  } else {
	     if (error_count == -1) { 
		cout << "       CIF ITEM: " << error_buf << endl << endl;
	     }
	  }
	  asc.read_success = 0; // FAIL
	  asc.read_error_message = error_buf; 
	  //
       } else {
	  // we read the coordinate file OK.
	  //
	  switch (MMDBManager->GetFileType())  {
	  case MMDB_FILE_PDB    :  cout << " PDB"         ;
	     break;
	  case MMDB_FILE_CIF    :  cout << " mmCIF"       ; 
	     break;
	  case MMDB_FILE_Binary :  cout << " MMDB binary" ;
	     break;
	  default:
	     cout << " Unknown (report as a bug!)\n";
	  }

	  MMDBManager->PDBCleanup(PDBCLEAN_ELEMENT);
	  
	  cout << " file " << pdb_name.c_str() << " has been read.\n";
	  asc.read_success = 1; // TRUE

	  // atom_selection_container.read_error_message = NULL; // its a string
	  asc.mol = MMDBManager;
       }
    }

    char *str = MMDBManager->GetSpaceGroup();
    if (str) { 
       std::string sgrp(str);
       std::cout << "Spacegroup: " << sgrp << "\n";
    } else {
       std::cout << "No Spacegroup found for this PDB file\n";
    } 
    
    std::cout << "Cell: "
	      << MMDBManager->get_cell().a << " "
	      << MMDBManager->get_cell().b << " "
	      << MMDBManager->get_cell().c << " "
	      << MMDBManager->get_cell().alpha << " "
	      << MMDBManager->get_cell().beta  << " "
	      << MMDBManager->get_cell().gamma << "\n";


    // 
    
    // Make handle_read_draw_molecule use make_asc which add the 
    // UDD "atom index".
    //
    if (asc.read_success) {
       asc = make_asc(asc.mol);

       // debug atom names
       if (0) { 
	  for (int i=0; i<asc.n_selected_atoms; i++) {
	     std::cout << i << " "
		       << asc.atom_selection[i]->GetChainID() << " "
		       << asc.atom_selection[i]->GetSeqNum() << " :"
		       << asc.atom_selection[i]->name << ":" <<std::endl;
	  }
       }

       
       fix_nucleic_acid_residue_names(asc);
       fix_away_atoms(asc);
       fix_wrapped_names(asc);
    }
    return asc; 
}

// Return the number of residue names changed.
//
// Tinker with asc as necessary.
int 
fix_nucleic_acid_residue_names(atom_selection_container_t asc) {

   int istat = 0;
   
   if (asc.n_selected_atoms > 0) { 

      int n_models = asc.mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) { 

	 CModel *model_p = asc.mol->GetModel(imod);
	 // model can legitimately be null if that particular model
	 // number was not in the PDB file.
	 if (model_p) { 
	    CChain *chain_p;
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
		     std::cout << "NULL chain in ... " << std::endl;
		  } else { 
		     int nres = chain_p->GetNumberOfResidues();
		     PCResidue residue_p;
		     for (int ires=0; ires<nres; ires++) { 
			residue_p = chain_p->GetResidue(ires);
			std::string residue_name(residue_p->name);

			if (residue_name == "T" ||
			    residue_name == "U" ||
			    residue_name == "A" ||
			    residue_name == "C" ||
			    residue_name == "G" ||
			    residue_name == "DA" ||
			    residue_name == "DG" ||
			    residue_name == "DT" ||
			    residue_name == "DC") {

			   istat += fix_nucleic_acid_residue_name(residue_p);
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   return istat;
}

int fix_nucleic_acid_residue_name(CResidue *r) {

   int istat=0;

   PCAtom *residue_atoms;
   int n_residue_atoms;
   bool found_o2_star = 0;

   r->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int i=0; i<n_residue_atoms; i++) {
      std::string atom_name(residue_atoms[i]->name);
      if (atom_name == " O2*") {
	 found_o2_star = 1;
	 break;
      }
      // I suppose we should call bases with O2' RNA too (not that it
      // does much good because the dictionary will not match the atom
      // names).
      if (atom_name == " O2'") {
	 found_o2_star = 1;
	 break;
      }
   }

   convert_to_old_nucleotide_atom_names(r);

   std::string res_name = r->name;
   std::string new_name_stub = res_name.substr(0,1);
   if (res_name == "DA" || res_name == "DT" ||
       res_name == "DC" || res_name == "DG")
      new_name_stub = res_name.substr(1,1);
       
   if (n_residue_atoms > 0)
      istat = 1;
   
   if (istat == 1) {
      if (found_o2_star) {
	 new_name_stub += "r";
      } else {
	 new_name_stub += "d";
      }
      r->SetResName(new_name_stub.c_str());
   }
   return istat;
}

// " H5'" -> "H5*1"
// " H5'" -> "H5*1"
// "H5''" -> "H5*2"
void
convert_to_old_nucleotide_atom_names(CResidue *r) {

   PCAtom *residue_atoms;
   int n_residue_atoms;
   r->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int i=0; i<n_residue_atoms; i++) {
      std::string atom_name(residue_atoms[i]->name);
      std::string name_orig = atom_name;
      std::string ele(residue_atoms[i]->element);
      char c3 = atom_name[2]; // 3rd char
      char c4 = atom_name[3]; // 4th char
      if (ele == " H") {
	 if (c3 == '\'') {
	    atom_name[2] = '*';
	    if (c4 == '\'')
	       atom_name[3] = '2';
	    else 
	       atom_name[3] = '1';
	 } else {
	    if (c4 == '\'') {
	       if (atom_name == " H5'")
		  atom_name = "H5*1";
	       else 
		  atom_name[3] = '*';
	    }
	 }
	 strncpy(residue_atoms[i]->name, atom_name.c_str(),5);
      } else {
	 // if it is not a hydrogen, simply change the prime to a star
	 if (c4 == '\'') {
	    atom_name[3] = '*';
	    strncpy(residue_atoms[i]->name, atom_name.c_str(),5);
	 }

	 if (atom_name == " OP1") {
	    atom_name = " O1P";
	    strncpy(residue_atoms[i]->name, atom_name.c_str(),5);
	 }
	 if (atom_name == " OP2") {
	    atom_name = " O2P";
	    strncpy(residue_atoms[i]->name, atom_name.c_str(),5);
	 }
      }
      // debug
      // std::cout << "from :" << name_orig << ": to :"
      // << atom_name << ":" << std::endl;
   }
}


int
fix_away_atoms(atom_selection_container_t asc) {

   int nat = 0;
   for (int i=0; i<asc.n_selected_atoms; i++) {
      if ((asc.atom_selection[i]->x > 9998.0) &&
	  (asc.atom_selection[i]->y > 9998.0) &&
	  (asc.atom_selection[i]->z > 9998.0)) {
	 asc.atom_selection[i]->x =  0.0;
	 asc.atom_selection[i]->y =  0.0;
	 asc.atom_selection[i]->z =  0.0;
	 nat++;
      }
   }
   return nat;
}

// Return the number of residue names changed.
//
// Tinker with asc as necessary.
int 
fix_wrapped_names(atom_selection_container_t asc) {

   int n_changed = 0;
   int uddHnd_old =
      asc.mol->RegisterUDString(UDR_ATOM , "initial hydrogen name");
   int uddHnd_new =
      asc.mol->RegisterUDString(UDR_ATOM , "new hydrogen name");
//    std::cout << "udd_old: create time " << uddHnd_old << std::endl;
//    std::cout << "udd_new: create time " << uddHnd_new << std::endl;

   // e.g. "3HB " -> " HB3", and "2HG2" -> "HG22"
   for (int i=0; i<asc.n_selected_atoms; i++) {
      // std::string ele(asc.atom_selection[i]->element);
      // if (ele == " H") {
      if (1) { 
	 std::string atom_name(asc.atom_selection[i]->name);
	 if (atom_name[0] == '1' ||
	     atom_name[0] == '2' ||
	     atom_name[0] == '3' ||
	     atom_name[0] == '4' ||
	     atom_name[0] == '*') {
	    // switch it.
	    std::string new_atom_name = atom_name.substr(1,3) + atom_name[0];
	    if (atom_name[3] != ' ') { 
	       if (atom_name[3] == ' ') {
		  new_atom_name = atom_name.substr(1,2) + atom_name[0];
		  new_atom_name += ' ';
	       }
	       if (atom_name[2] == ' ') {
		  new_atom_name = atom_name.substr(1,1) + atom_name[0];
		  new_atom_name += ' ';
		  new_atom_name += ' ';
	       }
	    } else {
	       // atom_name length is 3 presumably
	       new_atom_name = ' ';
	       new_atom_name += atom_name.substr(1,2) + atom_name[0];
	    } 
//   	    std::cout << "DEBUG:: atom_name switch :" <<  atom_name << ": -> :"
//   		      << new_atom_name << ":\n";
	    if (uddHnd_old >= 0)
	       asc.atom_selection[i]->PutUDData(uddHnd_old,
						asc.atom_selection[i]->name);
	    if (uddHnd_new >= 0)
	       asc.atom_selection[i]->PutUDData(uddHnd_new,
						new_atom_name.c_str());
	    asc.atom_selection[i]->SetAtomName(new_atom_name.c_str());
	    n_changed++;;
 	 } else {
	    // refmac calls it " H "
	    if (atom_name == " H0 ") {
	       std::string new_atom_name = " H  ";
	       if (uddHnd_old >= 0)
		  asc.atom_selection[i]->PutUDData(uddHnd_old,
						   asc.atom_selection[i]->name);
	       if (uddHnd_new >= 0)
		  asc.atom_selection[i]->PutUDData(uddHnd_new,
						   (char *) new_atom_name.c_str());
	       asc.atom_selection[i]->SetAtomName(new_atom_name.c_str());
	       n_changed++;
	    }
	 } 
      }
   }
   // std::cout << "done hydrogen names " << n_changed << std::endl;
   return n_changed; 
}



// Dump this.  Use instead GetCentreOfMass(selhand)
//
coot::Cartesian
centre_of_molecule(atom_selection_container_t SelAtom) {

   coot::Cartesian centre; // defaults construction at (0,0,0).
   coot::Cartesian rs;     // running sum
   PCAtom atom;

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
ostream& operator<<(ostream& s, CAtom &atom) {

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
ostream& operator<<(ostream& s, PCAtom atom) {

   //
   s << atom->GetModelNum() << "/" << atom->GetChainID() << "/"
     << atom->GetSeqNum()   << atom->GetInsCode() << "/"
     << atom->GetResName() << "/"
     << atom->name << " altLoc :" << atom->altLoc << ": segid :"
     << atom->segID << ":" << " pos: ("
     << atom->x << "," << atom->y << "," << atom->z
     << ") B-factor: " << atom->tempFactor;

   return s;

}
  

int
write_atom_selection_file(atom_selection_container_t asc,
			  const std::string &filename,
			  byte gz,
			  bool write_hydrogens,     // optional arg
			  bool write_aniso_records // optional arg
			  ) {

   int ierr = 0; 
   coot::util::remove_wrong_cis_peptides(asc.mol);
   CMMDBManager *mol = asc.mol;
   
   if (coot::is_mmcif_filename(filename)) {
      ierr = mol->WriteCIFASCII(filename.c_str());
      
   } else {

      if (! write_hydrogens) {
	 CMMDBManager *n = new CMMDBManager;
	 n->Copy(mol, MMDBFCM_All);
	 coot::delete_hydrogens_from_mol(n);
	 mol = n;
      }

      if (! write_aniso_records) {
	 CMMDBManager *n = new CMMDBManager;
	 n->Copy(mol, MMDBFCM_All);
	 coot::delete_aniso_records_from_atoms(n);
	 mol = n;
      }

      
      // we need to put the hydrogen names back to how they used to be
      // when we read in the pdb file and then put them put them back
      // to how they currently are!

      int udd_old = mol->GetUDDHandle(UDR_ATOM, "initial hydrogen name");
      int udd_new = mol->GetUDDHandle(UDR_ATOM, "new hydrogen name");
//       std::cout << "udd_old: " << udd_old << std::endl;
//       std::cout << "udd_new: " << udd_new << std::endl;
      char *str = 0; 
      if (udd_old > 0 && udd_new > 0) { 
	 for (int i=0; i<asc.n_selected_atoms; i++) {
	    str = 0; 
	    if (asc.atom_selection[i]->GetUDData(udd_old, str) == UDDATA_Ok) {
// 	       std::cout << "pre  UDD: " << asc.atom_selection[i]
// 			  << " gave udd str: " << str << std::endl;
	       asc.atom_selection[i]->SetAtomName(str);
	    } else {
	       // not a hydrogen, probably
	    }
	 }
      }
      ierr = mol->WritePDBASCII((char *)filename.c_str());
      // now put the names back
      if (udd_old > 0 && udd_new > 0) { 
      	 for (int i=0; i<asc.n_selected_atoms; i++) {
	    str = 0; 
      	    if (asc.atom_selection[i]->GetUDData(udd_new, str) == UDDATA_Ok) { 
       	       asc.atom_selection[i]->SetAtomName(str);
//        	       std::cout << "post UDD: " << asc.atom_selection[i]
//        			 << " gave udd str: " << str << std::endl;
	    } else {
	       // not a hydrogen, probably
       	    }
       	 }
      }
   }
   return ierr; 
}


void
coot::delete_hydrogens_from_mol(CMMDBManager *mol) {

   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
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
	    bool deleted = 0;
	    for (int iat=0; iat<n_atoms; iat++) {
	       at = residue_p->GetAtom(iat);
	       std::string ele(at->element);
	       if (ele == " H") {
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
coot::delete_aniso_records_from_atoms(CMMDBManager *mol) {

   std::cout << "ASET_Anis_tFac " << ASET_Anis_tFac << " " << ~ASET_Anis_tFac << std::endl;
   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
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
	       at->WhatIsSet &=  ~ASET_Anis_tFac;
	    }
	 }
      }
   }
}
