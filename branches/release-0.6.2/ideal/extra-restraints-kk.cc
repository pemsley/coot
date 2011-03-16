/* ideal/extra-restraints-kk.cc
 * 
 * Copyright 2011 by Kevin Keating
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
coot::restraints_container_t::add_extra_start_pos_restraints(const extra_restraints_t &extra_restraints) {

   for (unsigned int i=0; i<extra_restraints.start_pos_restraints.size(); i++) {
      
      CResidue *r_1 = NULL;
      CAtom *at_1 = 0;
      bool fixed_1 = 0;
      if (from_residue_vector) {
	 for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
	    if (coot::residue_spec_t(extra_restraints.start_pos_restraints[i].atom_1) ==
		coot::residue_spec_t(residues_vec[ir].second)) {
	       r_1 = residues_vec[ir].second;
	       fixed_1 = residues_vec[ir].first;
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
	    if (coot::residue_spec_t(extra_restraints.start_pos_restraints[i].atom_1) ==
		coot::residue_spec_t(SelResidue_local[ir])) {
	       r_1 = SelResidue_local[ir];
	    }
	 } 
	 mol->DeleteSelection(selHnd);

      }
      
      if (r_1) {
         PPCAtom residue_atoms_1 = 0;
         int n_residue_atoms_1;
         r_1->GetAtomTable(residue_atoms_1, n_residue_atoms_1);

         for (unsigned int iat=0; iat<n_residue_atoms_1; iat++) { 
            std::string atom_name_1(residue_atoms_1[iat]->name);
            if (atom_name_1 == extra_restraints.start_pos_restraints[i].atom_1.atom_name) {
               std::string alt_loc_1(residue_atoms_1[iat]->altLoc);
               if (alt_loc_1 == extra_restraints.start_pos_restraints[i].atom_1.alt_conf) {
                  at_1 = residue_atoms_1[iat];
                  break;
               }
            }
         }


         if (at_1) {
            int index_1 = -1; 
            at_1->GetUDData(udd_atom_index_handle, index_1);
            if (index_1 != -1) { 
               bool fixed_flag = fixed_check(index_1);
               if (! (fixed_flag)) {
                  //if the atom is fixed, then a start position retraint will have no effect, so we don't add it
                  add_user_defined_start_pos_restraint(START_POS_RESTRAINT, index_1, fixed_flag,
                      extra_restraints.start_pos_restraints[i].esd, 0.0);
               }
            }
         }  
      } 
   }
}

#endif // HAVE_GSL
