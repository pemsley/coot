
#include "coot_molecule.hh"

// merge molecules
// Return +1 as status of pair if we did indeed do a merge
// Return a list of new chain ids as the second.
//
// Note, very often add_molecules will only be of size 1.
//
// Recall that we will often be merging ligands into (otherwise quite
// complete) proteins.  In that case, in which we add a single residue
// it would be The Right Thing to Do if we could find a chain that
// consisted only of the same residue type as is the new ligand and
// add it to that chain. If that is the case, we should return a spec
// for the residue, not just the chain id.
//
// If we can't find a chain like that - or the new molecule contains
// more than one residue, we simply add new chains (with new chain
// ids) to this molecule.
//
// Question to self: how do I deal with different models?
//
std::pair<int, std::vector<merge_molecule_results_info_t> >
coot::molecule_t::merge_molecules(const std::vector<atom_selection_container_t> &add_molecules) {

   // pass this?
   coot::residue_spec_t merge_molecules_ligand_spec;
   const coot::residue_spec_t &spec = merge_molecules_ligand_spec;

   int istat = 0;
   make_backup(); // could be more clever, by doing this only when needed.
   std::vector<merge_molecule_results_info_t> resulting_merge_info;
   std::pair<bool, coot::residue_spec_t> done_merge_ligand_to_near_chain;
   done_merge_ligand_to_near_chain.first = false;

   std::vector<std::string> this_model_chains = coot::util::chains_in_molecule(atom_sel.mol);

   atom_sel.delete_atom_selection();

   for (unsigned int imol=0; imol<add_molecules.size(); imol++) {
      mmdb::Manager *adding_mol = add_molecules[imol].mol;
      if (add_molecules[imol].n_selected_atoms > 0) {
         int nresidues = coot::util::number_of_residues_in_molecule(add_molecules[imol].mol);

         // We need to set multi_residue_add_flag appropriately.  We
         // unset it if the molecule to be added has only one residue
         //
         bool multi_residue_add_flag = true;

         std::vector<std::string> adding_model_chains
            = coot::util::chains_in_molecule(add_molecules[imol].mol);

         if (nresidues == 1) {

            bool done_add_specific = merge_molecules_just_one_residue_at_given_spec(add_molecules[imol], spec);
            bool done_homogeneous_addition_flag = false;

            // by "homogeneous" I mean is there a chain of residues of the same type as that we are adding
            // e.g. an SO4 to a chain of SO4s? This don't happen (much?) these days.
            // These days, with one residue, we expect to call merge_ligand_to_near_chain()

            if (! done_add_specific)
               done_homogeneous_addition_flag = merge_molecules_just_one_residue_homogeneous(add_molecules[imol]);

            if (done_add_specific)
               multi_residue_add_flag = false;
            else
               multi_residue_add_flag = ! done_homogeneous_addition_flag;

            if (! done_homogeneous_addition_flag) {

               if (! done_add_specific)
                  done_merge_ligand_to_near_chain = merge_ligand_to_near_chain(adding_mol);

               if (done_merge_ligand_to_near_chain.first) {
                  merge_molecule_results_info_t mmr;
                  mmr.is_chain = false;
                  mmr.spec = done_merge_ligand_to_near_chain.second;
                  resulting_merge_info.push_back(mmr);
                  istat = 1;
               } else {

                  if (done_add_specific) {
                     // JED ligand addition
                     merge_molecule_results_info_t mmr;
                     mmr.is_chain = false;
                     mmr.spec = spec;
                     // std::cout << "---- JED case pushing back mmr " << mmr.spec << std::endl;
                     resulting_merge_info.push_back(mmr);
                     istat = 1;
                  }
               }

               // set multi_residue_add_flag if that was not a successful merge
               if (! done_add_specific)
                  multi_residue_add_flag = ! done_merge_ligand_to_near_chain.first;

            }
         }

         // Now that multi_residue_add_flag has been set properly, we use it...

         if (multi_residue_add_flag) {
            // return state
            std::pair<bool, std::vector<std::string> > add_state = try_add_by_consolidation(adding_mol);

            // some mild hacking, but we need to return a proper state and the added chains
            istat = 0;
            for (unsigned int i=0; i<add_state.second.size(); i++) {
               merge_molecule_results_info_t mmr;
               mmr.is_chain = true;
               mmr.chain_id = add_state.second[i];
               resulting_merge_info.push_back(mmr);
            }
            if (add_state.first) {
               // update_molecule_after_additions();
               // if (graphics_info_t::show_symmetry == 1)
               // update_symmetry();
               multi_residue_add_flag = false; // we've added everything for this mol.
               istat = add_state.first;
            }
         }

         // this should happen rarely these days...
         //
         if (multi_residue_add_flag) {

            std::vector<std::string> mapped_chains = map_chains_to_new_chains(adding_model_chains, this_model_chains);

            std::cout << "INFO:: Merge From chains: " << std::endl;
            for (unsigned int ich=0; ich<adding_model_chains.size(); ich++)
               std::cout << " :" << adding_model_chains[ich] << ":";
            std::cout << std::endl;
            std::cout << "INFO:: Merge To chains: " << std::endl;
            for (unsigned int ich=0; ich<mapped_chains.size(); ich++)
               std::cout << " :" << mapped_chains[ich] << ":";
            std::cout << std::endl;

            if (mapped_chains.size() != adding_model_chains.size()) {
               // can't continue with merging - no chains available.
               std::cout << "can't continue with merging - no chains available." << std::endl;
            } else {
               // fine, continue

               // Add the chains of the new molecule to this atom_sel, chain by chain.
               int i_add_model = 1;
               int i_this_model = 1;

               mmdb::Model *model_p = add_molecules[imol].mol->GetModel(i_add_model);
               mmdb::Model *this_model_p = atom_sel.mol->GetModel(i_this_model);

               int n_add_chains = model_p->GetNumberOfChains();

               for (int iaddchain=0; iaddchain<n_add_chains; iaddchain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(iaddchain);
                  mmdb::Chain *copy_chain_p = new mmdb::Chain;
                  copy_chain_p->Copy(chain_p);
                  copy_chain_p->SetChainID(mapped_chains[iaddchain].c_str());
                  this_model_p->AddChain(copy_chain_p);
                  this_model_chains.push_back(mapped_chains[iaddchain].c_str());
                  merge_molecule_results_info_t mmr;
                  mmr.is_chain = true;
                  mmr.chain_id = mapped_chains[iaddchain];
                  resulting_merge_info.push_back(mmr);
               }

               if (n_add_chains > 0) {
                  atom_sel.mol->FinishStructEdit();
                  // update_molecule_after_additions();
                  // if (graphics_info_t::show_symmetry == 1)
                  // update_symmetry();
                  coot::util::pdbcleanup_serial_residue_numbers(atom_sel.mol);
               }
               istat = 1;
            }
         }
      }
   }

   atom_sel = make_asc(atom_sel.mol);

   // fill_ghost_info(true, 0.7);

   // std::cout << "------- resulting_merge_info has size " << resulting_merge_info.size() << std::endl;
   if (false)
      if (!resulting_merge_info.empty())
         std::cout << "-------- resulting_merge_info[0] " << resulting_merge_info[0].spec << std::endl;

   return std::pair<int, std::vector<merge_molecule_results_info_t> > (istat, resulting_merge_info);
}

// return the multi_residue_add_flag
// done_homogeneous_addition_flag
//
bool
coot::molecule_t::merge_molecules_just_one_residue_homogeneous(atom_selection_container_t molecule_to_add) {

   // If there is a chain that has only residues of the same
   // type as is the (single) residue in the new adding
   // molecule then we add by residue addition to chain
   // rather than add by chain (to molecule).

   // Are there chains in this model that only consist of
   // residues of type adding_model_chains[0]?
   //
   bool done_homogeneous_addition_flag = false;

   bool has_single_residue_type_chain_flag = false;

   int i_this_model = 1;

   mmdb::Chain *add_residue_to_this_chain = NULL;

   mmdb::Model *this_model_p = atom_sel.mol->GetModel(i_this_model);

   int n_this_mol_chains = this_model_p->GetNumberOfChains();

   for (int ithischain=0; ithischain<n_this_mol_chains; ithischain++) {
      mmdb::Chain *this_chain_p = this_model_p->GetChain(ithischain);
      std::vector<std::string> r = coot::util::residue_types_in_chain(this_chain_p);

      if (r.size() == 1) {
         std::string adding_model_resname(molecule_to_add.atom_selection[0]->residue->GetResName());
         if (r[0] == adding_model_resname) {
            // poly-ala helices (say) should not go into concatenated residues in the same chain
            if (adding_model_resname != "ALA") {
               add_residue_to_this_chain = this_chain_p;
               has_single_residue_type_chain_flag = true;
               break;
            }
         }
      }
   }

   if (has_single_residue_type_chain_flag) {
      if (molecule_to_add.n_selected_atoms > 0) {
         mmdb::Residue *add_model_residue = molecule_to_add.atom_selection[0]->residue;
         copy_and_add_residue_to_chain(add_residue_to_this_chain, add_model_residue);
         done_homogeneous_addition_flag = true;
         atom_sel.mol->FinishStructEdit();
         // update_molecule_after_additions();
         // if (graphics_info_t::show_symmetry == 1)
         // update_symmetry();
      }
   }
   return done_homogeneous_addition_flag;
}

bool
coot::molecule_t::merge_molecules_just_one_residue_at_given_spec(atom_selection_container_t molecule_to_add,
                                                                 coot::residue_spec_t target_spec) {

   bool status = false;

   if (! target_spec.empty()) {
      mmdb::Residue *residue_p = get_residue(target_spec);
      if (! residue_p) {
         // more checks: does molecule_to_add have only one residue?
         int i_model = 1;

         int n_res = coot::util::number_of_residues_in_molecule(molecule_to_add.mol);

         if (n_res == 1) {
            mmdb::Model *this_model_p = atom_sel.mol->GetModel(i_model);
            mmdb::Chain *this_chain_p = this_model_p->GetChain(target_spec.chain_id.c_str());
            if (! this_chain_p) {
               this_chain_p = new mmdb::Chain;
               this_chain_p->SetChainID(target_spec.chain_id.c_str());
               this_model_p->AddChain(this_chain_p);
            } else {
               std::cout << "INFO:: merge_molecules_just_one_residue_at_given_spec() "
                         << " this chain not found in molecule (good)" << std::endl;
            }
            mmdb::Residue *r = coot::util::get_first_residue(molecule_to_add.mol);
            if (r) {
               make_backup();
               mmdb::Residue *new_residue_p = copy_and_add_residue_to_chain(this_chain_p, r);
               new_residue_p->seqNum = target_spec.res_no;
               status = true;
            }
         } else {
            if (true) // debug
               std::cout << "debug:: merge_molecules_just_one_residue_at_given_spec() oops "
                         << " n_res is " << n_res << std::endl;
         }
      } else {
         std::cout << "WARNING:: merge_molecules_just_one_residue_at_given_spec() residue already exists "
                   << "in molecule " << target_spec << std::endl;
      }
   } else {
      if (false) // debug
         std::cout << "merge_molecules_just_one_residue_at_given_spec() null residue spec" << std::endl;
   }

   if (status) {
      atom_sel.mol->FinishStructEdit();
      // update_molecule_after_additions();
      // if (graphics_info_t::show_symmetry == 1)
      // update_symmetry();
   }

   if (false) // debug
      std::cout << "merge_molecules_just_one_residue_at_given_spec() returns " << status << std::endl;

   return status;
}

// This doesn't do a backup or finalise model.
mmdb::Residue *
coot::molecule_t::copy_and_add_residue_to_chain(mmdb::Chain *this_model_chain,
                                                mmdb::Residue *add_model_residue,
                                                bool new_resno_by_hundreds_flag) {

   mmdb::Residue *res_copied = NULL;
   if (add_model_residue) {
      bool whole_res_flag = true;
      int udd_atom_index_handle = 1; // does this matter?
      bool add_this = true;
      // check for overlapping water (could be generalised for same residue type?!
      std::vector<mmdb::Residue *> close_residues;
      close_residues = coot::residues_near_residue(add_model_residue, atom_sel.mol, 0.05);
      for (unsigned int i=0; i<close_residues.size(); i++) {
         if (close_residues[i]->isSolvent() && add_model_residue->isSolvent()) {
            add_this = false;
            std::cout<<"INFO:: not adding water because of overlap\n"<<std::endl;
            break;
         }
      }
      if (add_this) {

         /* No - this does an implicit embed-in-chain - that is not what we want
         mmdb::Residue *residue_copy = coot::deep_copy_this_residue(add_model_residue,
                                                                    "",
                                                                    whole_res_flag,
                                                                    udd_atom_index_handle);
         */
         mmdb::Residue *residue_copy = coot::util::deep_copy_this_residue(add_model_residue);

         if (residue_copy) {
            std::pair<short int, int> res_info =
               next_residue_number_in_chain(this_model_chain, new_resno_by_hundreds_flag);
            int new_res_resno = 9999;
            if (res_info.first)
               new_res_resno = res_info.second;
            residue_copy->seqNum = new_res_resno; // try changing the seqNum before AddResidue().
            this_model_chain->AddResidue(residue_copy);
            res_copied = residue_copy;
         }
      }
   }
   return res_copied;
}

// return "" on failure
std::string
coot::molecule_t::suggest_new_chain_id(const std::string &current_chain_id) const {

   // current_chain_id is the chain_id in the molecule that we are adding to this one.

   std::string new_chain_id;

   // 20200110-PE no more specials at all
   // std::string r("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz#$%^&@?/~|-+=(){}:;.,'");
   std::string r("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");

   std::vector<std::string> existing;
   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   mmdb::Chain *chain_p;
   int n_chains = model_p->GetNumberOfChains();
   // do we have current_chain_id in our list of chains already?  If not, then simply
   // choose current_chain_id.
   bool found_it = false;
   for (int ichain=0; ichain<n_chains; ichain++) {
      mmdb::Chain *chain_p = model_p->GetChain(ichain);
      std::string chid = chain_p->GetChainID();
      if (chid == current_chain_id) {
         found_it = true;
         break;
      }
   }
   if (! found_it)
      new_chain_id = current_chain_id; // all done!

   // how about a multichar post-fix? We only want to do that if current_chain_id
   // what multichar to begin with (from a pdbx file):
   // (Wolfram Tempel)
   //
   if (new_chain_id.empty()) {
      if (current_chain_id.length() > 1) {
         std::string trial_chain_id = current_chain_id + "2";
         if (trial_chain_id.length() < (10-1)) { // magic mmdb chain id max length (mmdb_defs.h)
                                               // (we need a space for the terminal NULL too).
            bool found_it = false;
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               std::string chid = chain_p->GetChainID();
               if (chid == trial_chain_id) {
                  found_it = true;
                  break;
               }
            }
            if (! found_it)
               new_chain_id = trial_chain_id; // all done!
         }
      }
   }

   if (new_chain_id.empty()) { // not set yet
      for (int ichain=0; ichain<n_chains; ichain++) {
         chain_p = model_p->GetChain(ichain);
         existing.push_back(chain_p->GetChainID());
      }
      unsigned int l = r.length();
      std::vector<std::string> candidates(l);
      for (unsigned int i=0; i<l; i++)
         candidates[i] = r[i];

      for (unsigned int i=0; i<existing.size(); i++)
         candidates.erase(std::remove(candidates.begin(), candidates.end(), existing[i]), candidates.end());

      if (candidates.size())
         new_chain_id = candidates[0];
   }
   return new_chain_id;
}


bool
coot::molecule_t::is_het_residue(mmdb::Residue *residue_p) const {

   bool status = false;

   if (residue_p) {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for(int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            if (at->Het) {
               status =  true;
               break;
            }
         }
      }
   }
   return status;
}



// return state, max_resno + 1, or 0, 1 of no residues in chain.
//
// new_res_no_by_hundreds is default false
std::pair<short int, int>
coot::molecule_t::next_residue_number_in_chain(mmdb::Chain *w,
                                               bool new_res_no_by_hundreds) const {

   std::pair<short int, int> p(0,1);
   int max_res_no = -9999;

   if (w) {
      int nres = w->GetNumberOfResidues();
      mmdb::Residue *residue_p;
      if (nres > 0) {
         for (int ires=nres-1; ires>=0; ires--) {
            residue_p = w->GetResidue(ires);
            if (residue_p->seqNum > max_res_no) {
               max_res_no = residue_p->seqNum;
               bool is_het_residue_flag = is_het_residue(residue_p);
               if (is_het_residue_flag) {
                  p = std::pair<short int, int>(1, residue_p->seqNum+1);
               } else {
                  if (new_res_no_by_hundreds) {
                     if (max_res_no < 9999) {
                        int res_no = coot::util::round_up_by_hundreds(max_res_no+1);
                        p = std::pair<short int, int>(1, res_no+1);
                     }
                  } else {
                     if (max_res_no < 9999) {
                        p = std::pair<short int, int>(1, max_res_no+1);
                     }
                  }
               }
            }
         }
         if (! p.first) {
            //  first the first space starting from the front
            int test_resno_start = 1001;
            bool is_clear = false;
            while (! is_clear) {
               is_clear = true;
               for (int iser=0; iser<nres; iser++) {
                  int resno_res = w->GetResidue(iser)->seqNum;
                  if (resno_res >= test_resno_start) {
                     if (resno_res <= (test_resno_start+10)) {
                        is_clear = false;
                     }
                  }
                  if (! is_clear)
                     break;
               }
               test_resno_start += 100;
            }
            p = std::pair<short int, int> (1, test_resno_start);
         }
      }
   }
   return p;
}


// There is (or should be) only one residue in mol
std::pair<bool, coot::residue_spec_t>
coot::molecule_t::merge_ligand_to_near_chain(mmdb::Manager *mol) {

   bool done_merge = false;
   coot::residue_spec_t res_spec;

   mmdb::Residue *adding_residue_p = 0;

   { // set adding_residue
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         mmdb::Chain *chain_p;
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            chain_p = model_p->GetChain(ichain);
            if (chain_p) {
               int nres = chain_p->GetNumberOfResidues();
               if (nres > 0)
                  adding_residue_p = chain_p->GetResidue(0);
            }
         }
      }
   }

   if (adding_residue_p) {

      // OK, what atoms in this molecule are close to adding_residue?
      // First, we need a vector of atoms in adding_residue
      // Then a set of atoms that are close to those positions
      // Then find the chain that has most atoms in that set
      // Then find a suitable residue number
      // Then add the residue

      std::vector<clipper::Coord_orth> ligand_atom_positions;
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      adding_residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++) {
         mmdb::Atom *atom = residue_atoms[i];
         if (!atom->isTer()) {
            clipper::Coord_orth pt = coot::co(atom);
            ligand_atom_positions.push_back(pt);
         }
      }
      std::set<mmdb::Atom *> near_atoms; // atoms in the protein (this molecule)
      double dist_crit = 4.2;
      double dist_crit_sqrd = dist_crit * dist_crit;

      // this is a slow position by position check. It can be speeded up, if needed,
      // by fiddling with residue and molecules (i.e. copying and merging) and using
      // SelectAtoms().
      //
      int imod = 1;
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            if (chain_p) {
               int nres = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<nres; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     int n_atoms = residue_p->GetNumberOfAtoms();
                     for (int iat=0; iat<n_atoms; iat++) {
                        mmdb::Atom *at = residue_p->GetAtom(iat);
                        if (at) {
                           if (! at->isTer()) {
                              clipper::Coord_orth this_at_co = coot::co(at);
                              for (std::size_t j=0; j<ligand_atom_positions.size(); j++) {
                                 const clipper::Coord_orth &pos = ligand_atom_positions[j];
                                 double this_dist_sqrd = (this_at_co-pos).lengthsq();
                                 if (this_dist_sqrd < dist_crit_sqrd) {
                                    near_atoms.insert(at);
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
         if (near_atoms.size() > 1) {
            // make a map of the number of residue in each chain that are close to the ligand atoms
            std::map<mmdb::Chain *, int> chain_map;
            std::set<mmdb::Atom *>::const_iterator it;
            for (it=near_atoms.begin(); it!=near_atoms.end(); it++) {
               mmdb::Chain *this_chain = (*it)->GetChain();
               if (chain_map.find(this_chain) == chain_map.end()) {
                  chain_map[this_chain] = 1;
               } else {
                  chain_map[this_chain]++;
               }
            }

            // find the "best" chain for this added residue/ligand
            mmdb::Chain *max_atom_chain = 0;
            int n_atoms_max = 0;
            std::map<mmdb::Chain *, int>::const_iterator chain_it;
            for (chain_it=chain_map.begin(); chain_it!=chain_map.end(); chain_it++) {
               int n = chain_it->second;
               if (n > n_atoms_max) {
                  n_atoms_max = n;
                  max_atom_chain = chain_it->first;
               }
            }

            if (n_atoms_max > 0) {
               if (max_atom_chain) {
                   // does not backup or finalize or update
                   bool new_res_no_by_hundreds = true;
                   mmdb::Residue *new_residue_p =
                      copy_and_add_residue_to_chain(max_atom_chain, adding_residue_p, new_res_no_by_hundreds);
                   done_merge = true;
                   res_spec = coot::residue_spec_t(new_residue_p);
               }
            }
         }
      }
   }

   if (done_merge) {
      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      // update_molecule_after_additions();
      // if (graphics_info_t::show_symmetry == 1)
      // update_symmetry();
   }
   return std::pair<bool, coot::residue_spec_t> (done_merge, res_spec);
}

// return status and vector of resulting chain ids.
//
std::pair<bool, std::vector<std::string> >
coot::molecule_t::try_add_by_consolidation(mmdb::Manager *adding_mol) {

   bool status = false;
   std::vector<std::string> chain_ids;

   // 20180104 Don't merge molecules made of ALA.

   // for this molecule molecule, make a map of chains that have one
   // residue type.  Could well be empty (or perhaps consist of just a
   // water chain)
   //
   std::map<std::string, std::pair<int, mmdb::Chain *> > single_res_type_map;
   for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         std::vector<std::string> residue_types;
         mmdb::Residue *residue_p;
         for (int ires=0; ires<nres; ires++) {
            residue_p = chain_p->GetResidue(ires);
            std::string res_name(residue_p->GetResName());
            if (std::find(residue_types.begin(), residue_types.end(), res_name) == residue_types.end())
               residue_types.push_back(res_name);
            if (residue_types.size() > 1)
               break;
         }
         if (residue_types.size() == 1)
            single_res_type_map[residue_types[0]] = std::pair<int, mmdb::Chain *> (imod, chain_p);
      }
   }

   for(int imod = 1; imod<=adding_mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = adding_mol->GetModel(imod);
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         bool done_this_chain = false;
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         std::vector<std::string> residue_types;
         mmdb::Residue *residue_p;
         for (int ires=0; ires<nres; ires++) {
            residue_p = chain_p->GetResidue(ires);
            std::string res_name(residue_p->GetResName());
            if (std::find(residue_types.begin(), residue_types.end(), res_name) == residue_types.end())
               residue_types.push_back(res_name);
         }

         if (residue_types.size() == 1) {
            if (residue_types[0] != "ALA") {
               std::map<std::string, std::pair<int, mmdb::Chain *> >::const_iterator it =
                  single_res_type_map.find(residue_types[0]);
               if (it != single_res_type_map.end()) {
                  if (it->second.first == imod) {

                     // We got a match, now add all of adding_mol chain_p
                     // to this molecule's chain

                     // BL says:: we check in copy_and_add_chain_residues_to_chain if there
                     // are overlapping waters. Alternativley we could do it here already.

                     copy_and_add_chain_residues_to_chain(chain_p, it->second.second);
                     done_this_chain = true;
                     std::string cid = it->second.second->GetChainID();
                     if (std::find(chain_ids.begin(), chain_ids.end(), cid) == chain_ids.end())
                        chain_ids.push_back(cid);
                  }
               }
            }
         }

         if (! done_this_chain) {
            // copy whole chain to a new chain
            mmdb::Model *this_model_p = atom_sel.mol->GetModel(imod);
            if (this_model_p) {
               std::string current_chain_id = chain_p->GetChainID();
               std::string new_chain_id = suggest_new_chain_id(current_chain_id);
               mmdb::Chain *copy_chain_p = new mmdb::Chain;
               copy_chain_p->Copy(chain_p);
               copy_chain_p->SetChainID(new_chain_id.c_str());
               this_model_p->AddChain(copy_chain_p);
               if (std::find(chain_ids.begin(), chain_ids.end(), new_chain_id) == chain_ids.end())
                  chain_ids.push_back(new_chain_id);
            }
         }
         atom_sel.mol->FinishStructEdit();
         status = true;
      }
   }
   return std::pair<bool, std::vector<std::string> > (status, chain_ids);
}
