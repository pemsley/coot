/* ideal/regularize-minimol.cc
 * 
 * Copyright 2004, 2005 The University of York
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

#ifdef HAVE_GSL

#include "regularize-minimol.hh"
#include "simple-restraint.hh"


coot::minimol::molecule
coot::regularize_minimol_molecule(const coot::minimol::molecule &molin,
				  const coot::protein_geometry &geom) {

   // By far the most often we get here with molin as a molecule made
   // from a single residue that we are regularizing as a wiggly
   // ligand residue.  So we will minimize the molecule by minimizing
   // a set of fragments.
   // 

   coot::minimol::molecule m;
   mmdb::Manager *mol = molin.pcmmdbmanager();

   // get resno_1 and resno_2 and chain_id
   // For now we presume that we have just the one chain.
   if (molin.fragments.size() > 0) {
      int ifrag = 0; // can make this a for loop variable if adventurous.
      int resno_1 = molin[ifrag].min_res_no();
      int resno_2 = molin[ifrag].max_residue_number();

      clipper::Xmap<float> dummy_xmap;
   
      // coot::restraints_container_t restraints(resno_1,
      //   				      resno_2,
      //   				      have_flanking_residue_at_start,
      //   				      have_flanking_residue_at_end,
      //   				      have_disulfide_residues,
      //   				      altconf,
      //   				      chn,
      //   				      mol,
      //   				      fixed_atom_specs,
      //   				      &dummy_xmap);

      std::string chain_id = molin[ifrag].fragment_id;
      std::vector<std::pair<bool,mmdb::Residue *> > residues;
      for (int ires=resno_1; ires<=resno_2; ires++) {
         mmdb::Residue *r = coot::util::get_residue(chain_id, ires, "", mol);
         if (r) {
            std::pair<bool, mmdb::Residue *> p(false, r);
            residues.push_back(p);
         }
      }
      coot::restraints_container_t restraints(residues, geom, mol, &dummy_xmap);

      int n_threads_max = get_max_number_of_threads();
      int n_threads = n_threads_max -1;
      if (n_threads < 1) n_threads = 1;
      ctpl::thread_pool tp(n_threads);
      // std::cout << "set thread pool " << n_threads << std::endl;
      restraints.thread_pool(&tp, n_threads);

      coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_AND_NON_BONDED;
      bool do_residue_internal_torsions = false;
      bool do_trans_peptide_restraints = true;
      coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
      int imol = 0; // dummy
      int nrestraints = restraints.make_restraints(imol,
						   geom, flags,
						   do_residue_internal_torsions,
						   do_trans_peptide_restraints,
						   0.0, 0, true, true, false,
						   pseudos);

      if (nrestraints > 0) { 
	 restraints.minimize(flags);
      }

      m = coot::minimol::molecule(mol);

      // set m for return
//       std::cout << "         DEBUG:: output minimol: " << std::endl;
//       m.check();
//       std::cout << "========================================="
// 		<< std::endl << std::endl;
      

   }
   return m;

} 

#endif // HAVE_GSL
