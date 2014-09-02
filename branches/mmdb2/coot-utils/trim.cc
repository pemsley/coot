/* coot-utils/trim.cc
 * 
 * Copyright 2005, 2006 by Paul Emsley, The University of York
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

#include "coot-map-utils.hh"
#include "coot-trim.hh"

// return the number of trimmed atoms
int
coot::util::trim_molecule_by_map(CMMDBManager *mol,
				 const clipper::Xmap<float> &xmap,
				 float map_level,
				 short int remove_or_zero_occ_flag,
				 short int waters_only_flag) {

   int n_changed = 0;

   // Note mol->DeleteAtom(CAtom *at) doesn't exist
   // So we have to know the atom number in the residue
   // which means that we move away from simple atom_selection
   // looping to coordinate hierachy looping.
   //

   if (1) { 

      CModel *model_p = mol->GetModel(1);
   
      CChain *chain;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      if (nchains <= 0) { 
	 std::cout << "bad nchains in trim molecule " << nchains
		   << std::endl;
      } else { 
	 for (int ichain=0; ichain<nchains; ichain++) {
	    chain = model_p->GetChain(ichain);
	    if (chain == NULL) {  
	       // This should not be necessary. It seem to be a
	       // result of mmdb corruption elsewhere - possibly
	       // DeleteChain in update_molecule_to().
	       std::cout << "NULL chain in model_view_residue_button_info_t: "
			 << std::endl;
	    } else { 
	       int nres = chain->GetNumberOfResidues();
	       for (int ires=0; ires<nres; ires++) { 
		  PCResidue residue_p = chain->GetResidue(ires);
		  std::string resname = residue_p->name;
		  if (((resname == "WAT" || resname == "HOH") && waters_only_flag)
		      || !waters_only_flag) {

		     int n_atoms = residue_p->GetNumberOfAtoms();
		     for (int iat=0; iat<n_atoms; iat++) {

			CAtom *at = residue_p->GetAtom(iat);
			clipper::Coord_orth co(at->x, at->y, at->z);
			if (density_at_point(xmap, co) < map_level) {
			
			   // A baddie.  What do we do with it?  Set its
			   // occupancy to zero or delete it?

			   if (remove_or_zero_occ_flag == coot::util::TRIM_BY_MAP_DELETE) {
			      residue_p->DeleteAtom(iat);
			      n_changed++;
			   }
			
			   if (remove_or_zero_occ_flag == coot::util::TRIM_BY_MAP_ZERO_OCC) {
			      at->occupancy = 0.0;
			      n_changed++;
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   if ((n_changed > 0) && (remove_or_zero_occ_flag == coot::util::TRIM_BY_MAP_DELETE)) {
      mol->FinishStructEdit();
      mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
   }
   return n_changed;
}
	    


