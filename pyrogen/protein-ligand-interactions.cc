/* pyrogen/protein-ligand-interactions.cc
 * 
 * Copyright 2017 by Medical Research Council
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


#include "coot-utils/coot-h-bonds.hh"
#include "protein-ligand-interactions.hh"

// consider where a peptide is the ligand
void
coot::protein_ligand_interactions(mmdb::Residue *residue_p, mmdb::Manager *mol,
				  coot::protein_geometry *geom_p,
				  float h_bond_dist_max) {

   residue_spec_t spec(residue_p);

   int SelHnd_all = mol->NewSelection(); // d
   int SelHnd_lig = mol->NewSelection(); // d
   mol->SelectAtoms(SelHnd_all, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES,
		    "*", "*", "*", "*", "*");
   mol->SelectAtoms(SelHnd_lig, 0, spec.chain_id.c_str(),
		    spec.res_no, spec.ins_code.c_str(),
		    spec.res_no, spec.ins_code.c_str(),
		    "*", "*", "*", "*");

   coot::h_bonds hb;
   std::pair<bool, int> status = hb.check_hb_status(SelHnd_lig, mol, *geom_p);
   if (! status.first)
      std::cout << "WARNING:: no HB status on atoms of ligand\n";
   std::vector<h_bond> hbonds = hb.get_mcdonald_and_thornton(SelHnd_lig,
							     SelHnd_all,
							     mol, *geom_p, h_bond_dist_max);

   for (unsigned int i=0; i<hbonds.size(); i++) {
      if (true)
	 std::cout << "DEBUG:: in process_ligand() hbond [" << i << "] donor "
		   << coot::atom_spec_t(hbonds[i].donor) << "...to... "
		   << coot::atom_spec_t(hbonds[i].acceptor) << " with ligand donor flag "
		   << hbonds[i].ligand_atom_is_donor << std::endl;

      // override these 2 if ligand atom is donor
      //
      mmdb::Atom      *ligand_atom = hbonds[i].acceptor;
      mmdb::Atom *env_residue_atom = hbonds[i].donor;
      if (hbonds[i].ligand_atom_is_donor) {
	 ligand_atom = hbonds[i].donor;
	 env_residue_atom = hbonds[i].acceptor;
      }

      // ... other stuff

   }
   mol->DeleteSelection(SelHnd_all);
   mol->DeleteSelection(SelHnd_lig);

}
