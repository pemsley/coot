/* ideal/extra-restraints.cc
 * 
 * Copyright 2010  by The University of Oxford
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


// we don't want to compile anything if we don't have gsl
#ifdef HAVE_GSL
#include <stdexcept>

#include "simple-restraint.hh"

void
coot::restraints_container_t::add_extra_restraints(const extra_restraints_t &extra_restraints) {

   add_extra_bond_restraints(extra_restraints);
   add_extra_angle_restraints(extra_restraints);
   add_extra_torsion_restraints(extra_restraints);
   add_extra_start_pos_restraints(extra_restraints);
}

void
coot::restraints_container_t::add_extra_bond_restraints(const extra_restraints_t &extra_restraints) {

   // don't add the restraint if both the residues are fixed.
   // 
   for (unsigned int i=0; i<extra_restraints.bond_restraints.size(); i++) {
      CResidue *r_1 = NULL;
      CResidue *r_2 = NULL;
      CAtom *at_1 = 0;
      CAtom *at_2 = 0;
      bool fixed_1 = 0;
      bool fixed_2 = 0;
      if (from_residue_vector) {
	 for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
	    if (coot::residue_spec_t(extra_restraints.bond_restraints[i].atom_1) ==
		coot::residue_spec_t(residues_vec[ir].second)) {
	       r_1 = residues_vec[ir].second;
	       fixed_1 = residues_vec[ir].first;
	    }
	    if (coot::residue_spec_t(extra_restraints.bond_restraints[i].atom_2) ==
		coot::residue_spec_t(residues_vec[ir].second)) {
	       r_2 = residues_vec[ir].second;
	       fixed_2 = residues_vec[ir].first;
	    }
	 }
      } else {

	 // bleugh.
	 int selHnd = mol->NewSelection();
	 mol->Select (selHnd, STYPE_RESIDUE, 1,       // .. TYPE, iModel
		      chain_id_save.c_str(), // Chain(s)
		      istart_res, "*", // starting res
		      iend_res,   "*", // ending   res
		      "*",  // residue name
		      "*",  // Residue must contain this atom name?
		      "*",  // Residue must contain this Element?
		      "*",  // altLocs
		      SKEY_NEW // selection key
		      );
	 int nSelResidues_local = 0;
	 PPCResidue SelResidue_local= 0;
	 mol->GetSelIndex (selHnd, SelResidue_local, nSelResidues_local);
	 for (int ir=0; ir<nSelResidues_local; ir++) {
	    if (coot::residue_spec_t(extra_restraints.bond_restraints[i].atom_1) ==
		coot::residue_spec_t(SelResidue_local[ir])) {
	       r_1 = SelResidue_local[ir];
               // BL says:: this is certainly not working, since we dont have the residues_vec
               // not sure how to get the fixed flag. Too lazy to look up. Paul. Fixme!?
	       //fixed_1 = residues_vec[ir].first;
	       //
	       // 20100807 corrected now.
	       // 
	       fixed_1 = fixed_check(ir);
	    }
	    if (coot::residue_spec_t(extra_restraints.bond_restraints[i].atom_2) ==
		coot::residue_spec_t(SelResidue_local[ir])) {
	       r_2 = SelResidue_local[ir];
	       fixed_1 = fixed_check(ir);
	    }
	 } 
	 mol->DeleteSelection(selHnd);

      }
      
      if (r_1 && r_2) {
	 if (! (fixed_1 && fixed_2)) {
	    PPCAtom residue_atoms_1 = 0;
	    PPCAtom residue_atoms_2 = 0;
	    int n_residue_atoms_1;
	    int n_residue_atoms_2;
	    r_1->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
	    r_2->GetAtomTable(residue_atoms_2, n_residue_atoms_2);

	    for (unsigned int iat=0; iat<n_residue_atoms_1; iat++) { 
	       std::string atom_name_1(residue_atoms_1[iat]->name);
	       if (atom_name_1 == extra_restraints.bond_restraints[i].atom_1.atom_name) {
		  std::string alt_loc_1(residue_atoms_1[iat]->altLoc);
		  if (alt_loc_1 == extra_restraints.bond_restraints[i].atom_1.alt_conf) {
		     at_1 = residue_atoms_1[iat];
		     break;
		  }
	       }
	    }
	    for (unsigned int iat=0; iat<n_residue_atoms_2; iat++) { 
	       std::string atom_name_2(residue_atoms_2[iat]->name);
	       if (atom_name_2 == extra_restraints.bond_restraints[i].atom_2.atom_name) {
		  std::string alt_loc_2(residue_atoms_2[iat]->altLoc);
		  if (alt_loc_2 == extra_restraints.bond_restraints[i].atom_2.alt_conf) {
		     at_2 = residue_atoms_2[iat];
		     break;
		  }
	       }
	    }

	    if (at_1 && at_2) {
	       int index_1 = -1; 
	       int index_2 = -1;
	       at_1->GetUDData(udd_atom_index_handle, index_1);
	       at_2->GetUDData(udd_atom_index_handle, index_2);
	       if ((index_1 != -1) && (index_2 != -1)) { 
		  std::vector<bool> fixed_flags = make_fixed_flags(index_1, index_2);
		  add(BOND_RESTRAINT, index_1, index_2, fixed_flags,
		      extra_restraints.bond_restraints[i].bond_dist,
		      extra_restraints.bond_restraints[i].esd,
		      1.2 /* dummy value */);
	       }
	    } 
	 } 
      } 
   }
}

void
coot::restraints_container_t::add_extra_torsion_restraints(const extra_restraints_t &extra_restraints) {

   for (unsigned int i=0; i<extra_restraints.torsion_restraints.size(); i++) {
      CResidue *r_1 = NULL;
      CResidue *r_2 = NULL;
      CResidue *r_3 = NULL;
      CResidue *r_4 = NULL;
      CAtom *at_1 = 0;
      CAtom *at_2 = 0;
      CAtom *at_3 = 0;
      CAtom *at_4 = 0;
      bool fixed_1 = 0;
      bool fixed_2 = 0;
      bool fixed_3 = 0;
      bool fixed_4 = 0;
      if (from_residue_vector) {
	 for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
	    if (coot::residue_spec_t(extra_restraints.torsion_restraints[i].atom_1) ==
		coot::residue_spec_t(residues_vec[ir].second)) {
	       r_1 = residues_vec[ir].second;
	       fixed_1 = residues_vec[ir].first;
	    }
	    if (coot::residue_spec_t(extra_restraints.torsion_restraints[i].atom_2) ==
		coot::residue_spec_t(residues_vec[ir].second)) {
	       r_2 = residues_vec[ir].second;
	       fixed_2 = residues_vec[ir].first;
	    }
	    if (coot::residue_spec_t(extra_restraints.torsion_restraints[i].atom_3) ==
		coot::residue_spec_t(residues_vec[ir].second)) {
	       r_3 = residues_vec[ir].second;
	       fixed_3 = residues_vec[ir].first;
	    }
	    if (coot::residue_spec_t(extra_restraints.torsion_restraints[i].atom_4) ==
		coot::residue_spec_t(residues_vec[ir].second)) {
	       r_4 = residues_vec[ir].second;
	       fixed_4 = residues_vec[ir].first;
	    }
	 }
      } else {
	 
	 // bleugh.
	 int selHnd = mol->NewSelection();
	 mol->Select (selHnd, STYPE_RESIDUE, 1,       // .. TYPE, iModel
		      chain_id_save.c_str(), // Chain(s)
		      istart_res, "*", // starting res
		      iend_res,   "*", // ending   res
		      "*",  // residue name
		      "*",  // Residue must contain this atom name?
		      "*",  // Residue must contain this Element?
		      "*",  // altLocs
		      SKEY_NEW // selection key
		      );
	 int nSelResidues_local = 0;
	 PPCResidue SelResidue_local= 0;
	 mol->GetSelIndex (selHnd, SelResidue_local, nSelResidues_local);
	 for (int ir=0; ir<nSelResidues_local; ir++) {
	    if (coot::residue_spec_t(extra_restraints.torsion_restraints[i].atom_1) ==
		coot::residue_spec_t(SelResidue_local[ir]))
	       r_1 = SelResidue_local[ir];
	    if (coot::residue_spec_t(extra_restraints.torsion_restraints[i].atom_2) ==
		coot::residue_spec_t(SelResidue_local[ir]))
	       r_2 = SelResidue_local[ir];
	    if (coot::residue_spec_t(extra_restraints.torsion_restraints[i].atom_3) ==
		coot::residue_spec_t(SelResidue_local[ir]))
	       r_3 = SelResidue_local[ir];
	    if (coot::residue_spec_t(extra_restraints.torsion_restraints[i].atom_4) ==
		coot::residue_spec_t(SelResidue_local[ir]))
	       r_4 = SelResidue_local[ir];
	 }
	 mol->DeleteSelection(selHnd);
      }

      if (r_1 && r_2 && r_3 && r_4) {
	 PPCAtom residue_atoms = 0;
	 int n_residue_atoms;
	 r_1->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (int iat=0; iat<n_residue_atoms; iat++) {
	    if (coot::atom_spec_t(residue_atoms[iat]) == extra_restraints.torsion_restraints[i].atom_1) {
	       at_1 = residue_atoms[iat];
	       break;
	    } 
	 }
	 residue_atoms = 0; // just to be safe
	 r_2->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (int iat=0; iat<n_residue_atoms; iat++) {
	    if (coot::atom_spec_t(residue_atoms[iat]) == extra_restraints.torsion_restraints[i].atom_2) {
	       at_2 = residue_atoms[iat];
	       break;
	    } 
	 }
	 residue_atoms = 0; // just to be safe
	 r_3->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (int iat=0; iat<n_residue_atoms; iat++) {
	    if (coot::atom_spec_t(residue_atoms[iat]) == extra_restraints.torsion_restraints[i].atom_3) {
	       at_3 = residue_atoms[iat];
	       break;
	    } 
	 }
	 residue_atoms = 0; // just to be safe
	 r_4->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (int iat=0; iat<n_residue_atoms; iat++) {
	    if (coot::atom_spec_t(residue_atoms[iat]) == extra_restraints.torsion_restraints[i].atom_4) {
	       at_4 = residue_atoms[iat];
	       break;
	    } 
	 }

	 if (at_1 && at_2 && at_3 && at_4) {
	    int index_1 = -1; 
	    int index_2 = -1;
	    int index_3 = -1; 
	    int index_4 = -1;
	    at_1->GetUDData(udd_atom_index_handle, index_1);
	    at_2->GetUDData(udd_atom_index_handle, index_2);
	    at_3->GetUDData(udd_atom_index_handle, index_3);
	    at_4->GetUDData(udd_atom_index_handle, index_4);
	    if ((index_1 != -1) && (index_2 != -1) && (index_3 != -1) && (index_4 != -1)) { 
	       std::vector<bool> fixed_flags = make_fixed_flags(index_1, index_2, index_3, index_4);
	       if (fixed_1) fixed_flags[0] = 1;
	       if (fixed_2) fixed_flags[1] = 1;
	       if (fixed_3) fixed_flags[2] = 1;
	       if (fixed_4) fixed_flags[3] = 1;

	       std::cout << "DEBUG:: adding user-defined torsion restraint with fixed flags: "
			 << "[" << index_1 << " " << coot::atom_spec_t(atom[index_1]) << " " << fixed_flags[0] << "]  " 
			 << "[" << index_2 << " " << coot::atom_spec_t(atom[index_2]) << " " << fixed_flags[1] << "]  " 
			 << "[" << index_3 << " " << coot::atom_spec_t(atom[index_3]) << " " << fixed_flags[2] << "]  " 
			 << "[" << index_4 << " " << coot::atom_spec_t(atom[index_4]) << " " << fixed_flags[3] << "]  " 
			 << std::endl;
	       
	       add_user_defined_torsion_restraint(TORSION_RESTRAINT,
						  index_1, index_2, index_3, index_4,
						  fixed_flags,
						  extra_restraints.torsion_restraints[i].torsion_angle,
						  extra_restraints.torsion_restraints[i].esd,
						  1.2, // dummy value
						  extra_restraints.torsion_restraints[i].period);
	    }
	 } 
      }
   } 
}


#endif // HAVE_GSL
