/* ligand/test-ligand.cc
 * 
 * Copyright 2007 by The University of Oxford
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms ofn the GNU General Public License as published by
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
 * 02110-1301, USA.
*/

#include "mmdb-extras.h" 
#include "mmdb.h" 
#include "torsion-general.hh"

int main(int argc, char **argv) {

   int r=0;
   int istat = 1;

   if (argc > 1) {
      std::string pdb_filename = argv[1];
      atom_selection_container_t asc = get_atom_selection(pdb_filename);
      std::string chain_id = "B";
      int resno = 81;
      bool reverse_torsion = 0;
      CModel *model_p = asc.mol->GetModel(1);
      if (model_p) {
	 CChain *chain_p = model_p->GetChain(0);
	 if (chain_p) {
	    CResidue *residue_p = chain_p->GetResidue(0);
	    if (residue_p) {
	       std::vector<coot::atom_spec_t> torsion_general_atom_specs;
	       coot::atom_spec_t a1(chain_id, resno, "", " C  ", "");
	       coot::atom_spec_t a2(chain_id, resno, "", " CA ", "");
	       coot::atom_spec_t a3(chain_id, resno, "", " CB ", "");
	       coot::atom_spec_t a4(chain_id, resno, "", " CG ", "");
	       if (reverse_torsion) { 
		  torsion_general_atom_specs.push_back(a4);
		  torsion_general_atom_specs.push_back(a3);
		  torsion_general_atom_specs.push_back(a2);
		  torsion_general_atom_specs.push_back(a1);
	       } else {
		  torsion_general_atom_specs.push_back(a1);
		  torsion_general_atom_specs.push_back(a2);
		  torsion_general_atom_specs.push_back(a3);
		  torsion_general_atom_specs.push_back(a4);
	       }
		  
	       coot::torsion_general tg(residue_p, asc.mol,
					torsion_general_atom_specs);
	       double diff = 20;
	       istat = tg.change_by(diff);

	       asc.mol->WritePDBASCII("rotated.pdb");
	    }
	 }
      }
   }
   return r;
}
