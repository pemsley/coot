/* coords/test-phenix-geo.cc
 * 
 * Copyright 2014 by Medical Research Council
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#include "Bond_lines.hh"
#include "phenix-geo.hh"

int main(int argc, char **argv) {

   std::string file_name = "1yjp.geo";
   if (argc > 1 )
      file_name = argv[1];
   coot::phenix_geo::phenix_geometry pgb;
   pgb.parse(file_name);

   bool test_coords = false;
   if (argc > 2) {
      std::string t(argv[2]);
      if (t == "test-coords") {
	 test_coords = true;
      }
   }


   if (test_coords) {
      mmdb::Manager *mol = new mmdb::Manager;

      int ierr = mol->ReadCoorFile("1yjp.pdb");

      std::cout << "print ierr: " << ierr << std::endl;

      if (ierr != mmdb::Error_NoError) {
	 std::cout << "Problem reading pdb file "  << std::endl;
      } else {
	 if (false) {
	    int imod = 1;
	    mmdb::Model *model_p = mol->GetModel(imod);
	    if (! model_p) {
	       std::cout << "Null model" << std::endl;
	    } else {
	       mmdb::Chain *chain_p;
	       int n_chains = model_p->GetNumberOfChains();
	       for (int ichain=0; ichain<n_chains; ichain++) {
		  chain_p = model_p->GetChain(ichain);
		  if (! chain_p) {
		     std::cout << "Null chain" << std::endl;
		  } else {
		     int nres = chain_p->GetNumberOfResidues();
		     mmdb::Residue *residue_p;
		     mmdb::Atom *at;
		     for (int ires=0; ires<nres; ires++) {
			residue_p = chain_p->GetResidue(ires);
			if (! residue_p) {
			   std::cout << "Null residue " << std::endl;
			} else {
			   int n_atoms = residue_p->GetNumberOfAtoms();
			   // std::cout << "residue has " << n_atoms << " atoms " << std::endl;
			   for (int iat=0; iat<n_atoms; iat++) {
			      at = residue_p->GetAtom(iat);
			      std::cout << "   " << iat << " " << coot::atom_spec_t(at) << std::endl;
			   }
			}
		     }
		  }
	       }
	    }
	 }
	 Bond_lines_container blc(mol, pgb);
      }
   }
   return 0;
}
