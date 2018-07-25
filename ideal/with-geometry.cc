/* ideal/with-geometry.cc 
 * 
 * Copyright 2002, 2003 The University of York
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


#include <sys/types.h> // for stating
#include <sys/stat.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif

#include <iostream>
#include <string.h>
#include <math.h>

#ifndef  __MMDB_MMCIF__
#include <mmdb2/mmdb_utils.h>
#endif
  
#include <iostream>
#include <string>
#include <vector>

#include <mmdb2/mmdb_manager.h>
#include "coords/mmdb-extras.h"
#include "coords/mmdb.h"

#include "simple-restraint.hh"

#include "model-bond-deltas.hh"

int
main(int argc, char **argv) {

#ifndef HAVE_GSL
   std::cout << "We don't have GSL, this program does nothing" << std::endl;
#else 
   std::string dict_filename;
   coot::protein_geometry geom;

   if (argc < 2) {
      std::cout << "usage: " << argv[0] << " pdb_filename " << std::endl;

   } else {

      std::string argv_1 = argv[1];

      if (argv_1 == "--resolve-deltas") {

	 if (argc >= 3) {

	    geom.init_standard();
	    std::string pdb_file_name(argv[2]);
	    atom_selection_container_t asc = get_atom_selection(pdb_file_name, 1, 0); 
	    int imol = 0; // dummy
	    coot::model_bond_deltas deltas(asc.mol, imol, &geom);
	    deltas.resolve();

	 }

      } else {
      
	 geom.init_standard();
      
	 std::string pdb_file_name(argv[1]);

	 // if pdb_file_name does not exist -> crash?
	 atom_selection_container_t asc = get_atom_selection(pdb_file_name, 1, 0); 
	 //coot::restraints_container_t restraints(asc);

	 // So, we provide easy(?) access to the atoms of next and
	 // previous residues (to those in the atom selection
	 // moving_residue_atoms).  It is also possible to select "fixed"
	 // atoms in the graphics (so that they don't move).  Let's
	 // provide a vector of indices in the the moving_residue_atoms
	 // array to define those (lovely mixture of styles - heh).
	 //
	 std::vector<coot::atom_spec_t> fixed_atom_specs;

	 // This interface has been withdrawn because we need the whole
	 // molecule (acutally, a pointer to it) to do some atom selection.
	 // 
	 //       coot::restraints_container_t restraints(asc.atom_selection, // moving_residue_atoms,
	 // 					      asc.n_selected_atoms,
	 // 					      previous_residue,
	 // 					      next_residue,
	 // 					      fixed_atoms);
      
	 int istart_res = 72;
	 int iend_res   = 73;  // ropey.pdb

	 istart_res = 41;
	 iend_res = 41;   // bad-chiral.pdb

	 istart_res = 289;  // link torsion restraints frag.pdb
	 iend_res   = 290;

	 //      int istart_res = 1031;
	 //      int iend_res   = 1033;

	 char *chain_id   = asc.atom_selection[0]->residue->GetChainID();
	 std::string altloc("");

	 short int have_flanking_residue_at_start = 0;
	 short int have_flanking_residue_at_end   = 0;
	 short int have_disulfide_residues = 0;

	 clipper::Xmap<float> dummy;
      
	 coot::restraints_container_t restraints(istart_res,
						 iend_res,
						 have_flanking_residue_at_start,
						 have_flanking_residue_at_end,
						 have_disulfide_residues,
						 altloc,
						 chain_id,
						 asc.mol,
						 fixed_atom_specs, dummy);


	 // coot::restraint_usage_Flags flags = coot::BONDS;
	 // coot::restraint_usage_Flags flags = coot::BONDS_AND_ANGLES;
	 // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_AND_TORSIONS;
	 // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_AND_PLANES;
	 // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_TORSIONS_AND_PLANES; 
	 // coot::restraint_usage_Flags flags = coot::TORSIONS; ok
	 coot::restraint_usage_Flags flags = coot::BONDS_AND_PLANES;
	 // flags = coot::NON_BONDED;
	 flags = coot::CHIRAL_VOLUMES;
	 flags = coot::TORSIONS;
	 // flags = coot::BONDS_ANGLES_TORSIONS_AND_PLANES;
	 // flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRAL;

	 coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
	 bool do_trans_peptide_restraints = false;
	 bool do_link_restraints = true;
	 bool do_flank_restraints = true;
	 int imol = 0; // dummy
	 restraints.make_restraints(imol, geom, flags, 1, do_trans_peptide_restraints,
				    0.0, 0, pseudos, do_link_restraints, do_flank_restraints);

	 // restraints.set_do_numerical_gradients();
	 restraints.minimize(flags);
	 restraints.write_new_atoms("new.pdb");

      }
   }

#endif // HAVE_GSL

   return 0; 
} 
