/* src/molecule-class-info-residues.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2007 by Paul Emsley
 * Copyright 2007, 2008, 2009 by The University of Oxford
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

// Must include these headers to get molecule_class_info_t.h to parse.
//

#include <string>
#include <stdexcept>

#include "mmdb_manager.h"
#include "mmdb-extras.h"
#include "Cartesian.h"
#include "mmdb-crystal.h"
// 
#include "molecule-class-info.h"
#include "coot-utils.hh"


// 1: success
// 0: failure
// 
bool
molecule_class_info_t::progressive_residues_in_chain_check_by_chain(const char *chain_id) const {

   bool r = 0;
   
   if (atom_sel.n_selected_atoms > 0) { 
      int nchains = atom_sel.mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) {
	 CChain *chain_p = atom_sel.mol->GetChain(1,ichain);
	 std::string mol_chain_id(chain_p->GetChainID());
	 if (mol_chain_id == std::string(chain_id)) { 
	    r = coot::progressive_residues_in_chain_check(chain_p);
	    break;
	 }
      }
   }

   return r; 
} 


// only apply charges if the molecule contains lots of hydrogens.
void
molecule_class_info_t::apply_charges(const coot::protein_geometry &geom) {

   // More than 20% of the atoms have to be hydrogens for use to set
   // the charges on all the atoms (to something other than
   // CXX_UNSET_CHARGE).
   
   float fraction_hydrogens = 0.2;

   if (atom_sel.n_selected_atoms > 0) { 
      CMMDBManager *mol = atom_sel.mol;

      int n_H = 0;
      int n_all = 0; 
      for (int i=0; i<atom_sel.n_selected_atoms; i++) {
	 std::string ele(atom_sel.atom_selection[i]->element);
	 if (ele == " H" || ele == " D") {
	    n_H++; 
	 }
	 n_all++;
      }
      
      if ( (float(n_H)/float(n_all) > fraction_hydrogens) || n_all < 40) {

	 // first set all atom charges to unset:
	 for (int i=0; i<atom_sel.n_selected_atoms; i++)
	    atom_sel.atom_selection[i]->charge = CXX_UNSET_CHARGE;


	 // Now add real charges from the dictionary
	 // 
	 int imod = 1;
	 CModel *model_p = mol->GetModel(imod);
	 CChain *chain_p;
	 int nchains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<nchains; ichain++) {
	    chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    PCResidue residue_p;
	    for (int ires=0; ires<nres; ires++) { 
	       residue_p = chain_p->GetResidue(ires);
	       std::string res_type = residue_p->GetResName();
	       std::pair<short int, coot::dictionary_residue_restraints_t> rp = 
		  geom.get_monomer_restraints(res_type);
	       if (rp.first) {
		  try { 
		     coot::dipole p(rp.second, residue_p);
		     p.fill_charged_atoms(residue_p, rp.second);
		  }
		  catch (std::runtime_error mess) {
		     std::cout << mess.what() << std::endl;
		  }
	       }
	    }
	 }
      }
   } 
} 

int
molecule_class_info_t::assign_hetatms() {

   int r = 0;
   for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
      int imod = 1;
      CModel *model_p = atom_sel.mol->GetModel(imod);
      CChain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 PCResidue residue_p;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    r += coot::hetify_residue_atoms_as_needed(residue_p);
	 }
      }
   }
   return r;
}
