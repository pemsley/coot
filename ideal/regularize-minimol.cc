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


void
coot::refine_minimol_fragment(coot::minimol::fragment &frag,
                              const coot::protein_geometry &geom,
                              const clipper::Xmap<float> &xmap,
                              float weight,          // = 60.0,
                              bool  do_GM,           //   = false,
                              float GM_distance_lim, // = 4.0f,
                              float GM_alpha         // = 0.02f
                              ) {

   int n_steps_max = 3000;

   auto update_frag_atoms = [] (minimol::fragment &frag,
                                mmdb::Manager *mol) {
                               std::vector<minimol::atom *> frag_atoms = frag.select_atoms_serial();

                               int SelHnd = mol->NewSelection();
                               mmdb::Atom **atom_selection = NULL;
                               int n_selected_atoms = 0;
                               mol->SelectAtoms(SelHnd, 0,
                                                frag.fragment_id.c_str(),
                                                mmdb::ANY_RES, // starting resno, an int
                                                "*", // any insertion code
                                                mmdb::ANY_RES, // ending resno
                                                "*", // ending insertion code
                                                "*", // any residue name
                                                "*", // atom name
                                                "*", // elements
                                                "*"  // alt loc.
                                                );
                               mol->GetSelIndex(SelHnd, atom_selection, n_selected_atoms);

                               unsigned int n_selected_atoms_ui = n_selected_atoms;
                               if (n_selected_atoms_ui == frag_atoms.size()) {
                                  for (unsigned int iat=0; iat<frag_atoms.size(); iat++) {
                                     minimol::atom *at  = frag_atoms[iat];
                                     clipper::Coord_orth old_pos = at->pos;
                                     clipper::Coord_orth new_pos = co(atom_selection[iat]);
                                     std::cout << "updating atom " << iat << " " << old_pos.format() << " " << new_pos.format() << std::endl;
                                     at->pos = new_pos;
                                  }
                               }
                            };

   std::vector<std::pair<bool,mmdb::Residue *> > residues;
   minimol::molecule m_mol(frag);
   mmdb::Manager *mol = m_mol.pcmmdbmanager();

   mmdb::Model *model_p = mol->GetModel(1);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int n_res = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<n_res; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            if (residue_p) {
               if (residue_p->GetNumberOfAtoms() > 0) {
                  residues.push_back(std::make_pair(true, residue_p));
               }
            }
         }
      }
   }

   if (! residues.empty()) {
      int n_threads = 1; // needs optimizing
      ctpl::thread_pool tp(n_threads);
      coot::restraints_container_t restraints(residues, geom, mol, &xmap);
      restraints.thread_pool(&tp, n_threads);

      coot::restraint_usage_Flags flags = coot::TYPICAL_RESTRAINTS; // include GM
      bool do_residue_internal_torsions = false;
      bool do_trans_peptide_restraints = true;
      coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
      int imol = 0; // dummy
      restraints.make_restraints(imol,
                                 geom, flags,
                                 do_residue_internal_torsions,
                                 do_trans_peptide_restraints,
                                 0.0, 0, true, true, false,
                                 pseudos);

      restraints.add_map(weight);

      if (do_GM) {
         // how do I make and add GM distance restraints?
         restraints.set_geman_mcclure_alpha(GM_alpha);
      }

      restraints.minimize(flags, n_steps_max, false);

      update_frag_atoms(frag, mol);

   }
}


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
      //                                         resno_2,
      //                                         have_flanking_residue_at_start,
      //                                         have_flanking_residue_at_end,
      //                                         have_disulfide_residues,
      //                                         altconf,
      //                                         chn,
      //                                         mol,
      //                                         fixed_atom_specs,
      //                                         &dummy_xmap);

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
      restraints.set_quiet_reporting();

      int n_threads_max = get_max_number_of_threads();
      int n_threads = n_threads_max -1;
      if (n_threads < 1) n_threads = 1;
      // might it be a good idea, if we are optimizing several ligands at the same time
      // that each ligand uses just one thread?
      // n_threads = 1; // apparently not. I don't understand why. I don't understand
      // what gdb info threads is telling me. Perhaps that thread_pool.push() is slow?
      n_threads = 1;
      ctpl::thread_pool tp(n_threads);
      // std::cout << "set thread pool " << n_threads << std::endl;
      restraints.thread_pool(&tp, n_threads);

      coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_TORSIONS_NON_BONDED_CHIRALS_AND_PLANES;
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
//       << std::endl << std::endl;


   }
   return m;

}



#endif // HAVE_GSL
