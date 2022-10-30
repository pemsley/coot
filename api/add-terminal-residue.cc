
#include <iostream>
#include "add-terminal-residue.hh"

#include "coords/mmdb-extras.h"
#include "coot-utils/coot-coord-utils.hh"
#include "ligand/residue_by_phi_psi.hh"

// We need to find the serial number of the residue after the residue
// we want to insert (i.e. the new residue will be inserted just
// before the residue whose serial number we return).
//
// return -1 on error.
std::pair<int, mmdb::Residue *>
find_serial_number_for_insert(mmdb::Manager *mol,
                              int seqnum_for_new,
                              const std::string &ins_code_for_new,
                              const std::string &chain_id) {

   int iserial_no = -1;
   std::pair<int, std::string> current_diff(999999, "");
   int n_chains = mol->GetNumberOfChains(1);
   mmdb::Residue *res = NULL;

   for (int i_chain=0; i_chain<n_chains; i_chain++) {

      mmdb::Chain *chain_p = mol->GetChain(1,i_chain);

      if (chain_p) {

         std::string mol_chain(chain_p->GetChainID());

         if (chain_id == mol_chain) {

            // Find the first residue that has either the residue number or insertion code
            // greater than the passed parameters

            int nres = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<nres; ires++) { // ires is a serial number
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);

               // we need to consider insertion codes here

               int diff = residue_p->GetSeqNum() - seqnum_for_new;

               if (diff > 0) {
                  res = residue_p;
                  iserial_no = ires;
                  break;
               } else {
                  if (diff == 0) {
                     std::string ins_code_this = residue_p->GetInsCode();
                     if (ins_code_this > ins_code_for_new) {
                        res = residue_p;
                        iserial_no = ires;
                        break;
                     }
                  }
               }
            }
         }
      }
   }
   return std::pair<int, mmdb::Residue *> (iserial_no, res);
}

void
coot::remove_TER_internal(mmdb::Manager *mol, mmdb::Residue *residue_p) {

   int n_residue_atoms;
   mmdb::PPAtom residue_atoms;
   bool deleted = 0;
   if (residue_p) {
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++) {
	 if (residue_atoms[i]->isTer()) {
	    residue_p->DeleteAtom(i);
	    deleted = 1;
	 }
      }
   }
   if (deleted) {
      mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      mol->FinishStructEdit();
   }
}

// The atom_index is the atom index of the clicked atom.
//
// Initially, this routine tested for real terminii.
//
// The one day EJD asked me how I would build a few missing residues
// (in a loop or so) that arp-warp missed.  I said "add terminal
// residue".  Then I tried it and it failed, of course because that
// residue was not a real terminus.  On reflection, I don't think that
// it should fail in this situation.
//
// So, I don't want to test for a real terminus.
//
// Let's see if this residue has a residue on one side of it, but not
// the other.  If so, return N or C depending on whether the other
// residue is upstream or not.
//
// Return not-terminal-residue if 0 neighbours, return "M" (for mid)
// for both neighbours present (used by
// graphics_info_t::execute_add_terminal_residue()).  Realise that
// this is a bit of a kludge, because usually, the terminal type
// refers to the residue that we clicked on not the missing residue.
//
// In the "M" case, we refer to the missing residue.
//
// "singleton" is a possibile terminal type - for cases where this
// residue does not have neighbours.
//
// Note that this ignores altlocs and insertion codes. It should do
// altlocs at least.
//
std::string
coot::get_term_type(mmdb::Residue *residue_p, mmdb::Manager *mol) {

   std::string term_type = "not-terminal-residue"; // returned thing

   if (! residue_p) return term_type;
   if (! mol) return term_type;

   int ires_atom = residue_p->GetSeqNum();
   mmdb::Chain *chain = residue_p->GetChain();
   int nres = chain->GetNumberOfResidues();

   // including tests needed for single missing residue:
   short int has_up_neighb = 0;
   short int has_down_neighb = 0;
   short int has_up_up_neighb = 0;
   short int has_down_down_neighb = 0;

   // Check for neighbouring residues to the clicked atom. Don't count
   // waters as neighbours.
   //
   for (int ires=0; ires<nres; ires++) {
      mmdb::PResidue res = chain->GetResidue(ires);
      if (res) { // could have been deleted (NULL)
         if (res->GetSeqNum() == (ires_atom + 1))
            has_up_neighb = 1;
         if (res->GetSeqNum() == (ires_atom + 2))
            has_up_up_neighb = 1;
         if (res->GetSeqNum() == (ires_atom - 1))
            has_down_neighb = 1;
         if (res->GetSeqNum() == (ires_atom - 2))
            has_down_down_neighb = 1;
      }
   }

   if ( (has_up_neighb + has_down_neighb) == 1 ) {
      if (has_up_neighb)
         term_type = "N";
      if (has_down_neighb)
         term_type = "C";
   }

   if ((has_up_neighb == 0) && (has_down_neighb == 0))
      term_type = "singleton";

   // Now test for missing single residue, "M" (mid):
   //
   if ( (!has_up_neighb) && has_up_up_neighb)
      term_type = "MC"; // missing middle res, treat as C term

   if ( (!has_down_neighb) && has_down_down_neighb)
      term_type = "MN"; // missing middle res, treat as N term

   // std::cout << "DEBUG:: get_term_type Returning residue type " << term_type << std::endl;

   return term_type;
}


// this function should have its own file and header in coot-utils

// put this in the coot namespace
void insert_coords(mmdb::Manager *mol, int udd_atom_index,
                   atom_selection_container_t asc_for_atoms_to_be_added) {

   // run over each chain, residue of the asc (if terminal residue
   // fit only one chain, one residue, of course).

   // std::cout << "######################### insert_coords()! " << mol << std::endl;

   auto debug_asc_for_atoms_to_be_inserted = [] (atom_selection_container_t asc) {
      int imod = 1;
      mmdb::Model *model_p = asc.mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            std::cout << "chain " << chain_p << " of " << n_chains << " chains " << std::endl;
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               std::cout << "residue " << residue_p << " of "  << n_res << " residues" << std::endl;
               if (residue_p) {
                  int n_atoms = residue_p->GetNumberOfAtoms();
                  for (int iat=0; iat<n_atoms; iat++) {
                     mmdb::Atom *at = residue_p->GetAtom(iat);
                     std::cout << "  asc_for_atoms_to_be_added: " << coot::atom_spec_t(at) << std::endl;
                  }
               }
            }
         }
      }
   };

   asc_for_atoms_to_be_added.mol->WritePDBASCII("asc_for_atoms_to_be_added.pdb");
   debug_asc_for_atoms_to_be_inserted(asc_for_atoms_to_be_added);

   bool inserted = false; // not inserted yet
   int imod = 1;
   mmdb::Model *adding_frag_model_p = asc_for_atoms_to_be_added.mol->GetModel(imod);
   int adding_frag_n_chains = adding_frag_model_p->GetNumberOfChains();
   for (int i_adding_frag_chain=0; i_adding_frag_chain<adding_frag_n_chains; i_adding_frag_chain++) {
      mmdb::Chain *adding_frag_chain = asc_for_atoms_to_be_added.mol->GetChain(imod, i_adding_frag_chain);
      int nres_adding_frag = adding_frag_chain->GetNumberOfResidues();

      std::cout << "DEBUG:: There are " << nres_adding_frag << " residues in "
                << "protein_chain (chain id: " << adding_frag_chain->GetChainID() << ")." << std::endl;

      for (int ires_adding_frag=0; ires_adding_frag<nres_adding_frag; ires_adding_frag++) {
         std::cout << "------------------- New adding frag residue " << ires_adding_frag << std::endl;
         mmdb::Residue *adding_frag_residue = adding_frag_chain->GetResidue(ires_adding_frag);

         if (! adding_frag_residue) continue;

         // Now find the corresponding chain in our (protein) mol:

         int imodel = 1;
         int n_chains = mol->GetNumberOfChains(imodel);
         for (int i_chain=0; i_chain<n_chains; i_chain++) {

            std::cout << "------------------- New protein chain " << i_chain << std::endl;

            mmdb::Chain *chain = mol->GetChain(imod, i_chain);

            // test chains
            std::string protein_chain_str(chain->GetChainID());
            std::string adding_frag_chain_str(adding_frag_chain->GetChainID());

            if (true)
               std::cout << "comparing chain ids :" << protein_chain_str << ": :"
                         << adding_frag_chain_str << ":" << std::endl;

            if (protein_chain_str == adding_frag_chain_str) {

               // insert that residue!

               bool embed_in_chain_flag = false;
               // use the coot-utils function, not this one
               mmdb::Residue *res = coot::deep_copy_this_residue_old_style(adding_frag_residue, "", 1, udd_atom_index, embed_in_chain_flag);

               if (false) {
                  if (res->GetNumberOfAtoms() > 0) {
                     std::cout  << "DEBUG:: inserting residue in chain " << protein_chain_str
                                << " internal-index " << ires_adding_frag << " of " << nres_adding_frag
                                << " residue number " << adding_frag_residue->GetSeqNum()
                                << " with " << adding_frag_residue->GetNumberOfAtoms() << " atoms "
                                << " i_chain = " << i_chain << " ires_asc = " << ires_adding_frag << std::endl;
                  }
               }

               std::pair<int, mmdb::Residue *> serial_number =
                  find_serial_number_for_insert(mol,
                                                adding_frag_residue->GetSeqNum(),
                                                std::string(adding_frag_residue->GetInsCode()),
                                                protein_chain_str);

               std::cout << "DEBUG:: returned serial_number: " << serial_number.first << " " << serial_number.second << std::endl;

               if (res) {

                  if (serial_number.first != -1) {
                     // insert at this position (other residues are
                     // shifted up).
                     // std::cout << "########### InsResidue() " << std::endl;
                     chain->InsResidue(res, serial_number.first);
                     coot::copy_segid(serial_number.second, res);
                     inserted = 1;

                  } else {

                     // std::cout << "DEBUG:: insert_coords_internal() add residue\n";
                     mmdb::Residue *last_residue = coot::util::get_last_residue_in_chain(chain);
                     if (last_residue) {

                        if (false) { // debug::
                           int nat = last_residue->GetNumberOfAtoms();
                           for (int iat=0; iat<nat; iat++) {
                              std::cout << iat << " of " << nat << " "
                                        << last_residue->GetAtom(iat) << std::endl;
                           }
                        }

                        chain->AddResidue(res);
                        coot::copy_segid(last_residue, res);
                        inserted = 1;
                     }
                  }
               }
            }
            //if (inserted) break;
         }


         if (! inserted) {
            // OK, there was no chain in the current mol that matches
            // the chain of the asc.
            // Let's copy the asc chain and add it to atom_sel.mol
            mmdb::Chain *new_chain = new mmdb::Chain;
            mmdb::Model *this_model = mol->GetModel(imodel);
            this_model->AddChain(new_chain);
            new_chain->SetChainID(adding_frag_chain->GetChainID());

            std::cout << "DEBUG:: insert_coords() Creating a new chain " << adding_frag_chain->GetChainID() << std::endl;
            bool embed_in_chain_flag = false;
            // use the coot-utils function, not this one
            mmdb::Residue *res = coot::deep_copy_this_residue_old_style(adding_frag_residue, "", 1, udd_atom_index, embed_in_chain_flag);
            if (res) {
               new_chain->AddResidue(res);
               mol->FinishStructEdit(); // so that we don't keep adding a
                                        // new Chain to atom_sel.mol
            }

         }
         //if (inserted) break;
      }
      //if (inserted) break;
   }
   mol->FinishStructEdit();
   // update_molecule_after_additions();
   mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
   // do we need the residue index function here (I don't recall its name!)

   coot::util::pdbcleanup_serial_residue_numbers(mol);

}


// place the O (because we have added a new residue)
bool
move_atom(const std::string &atom_name_in, mmdb::Residue *res_p, const clipper::Coord_orth &new_O_pos) {

   // just change the position of the first atom that matches atom_name_in
   bool done = false;

   mmdb::Atom **residue_atoms = 0;
   int n_residue_atoms;
   res_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int i=0; i<n_residue_atoms; i++) {
      mmdb::Atom *at = residue_atoms[i];
      std::string atom_name(at->name);
      if (atom_name == atom_name_in) {
	 at->x = new_O_pos.x();
	 at->y = new_O_pos.y();
	 at->z = new_O_pos.z();
	 done = true;
	 break;
      }
   }
   return done;
}



// move to coot utils
//
// Return 0 on failure.
int
move_std_residue(mmdb::Residue *moving_residue, const mmdb::Residue *reference_residue) {

   std::map<std::string, clipper::RTop_orth> rtops =
      coot::util::get_ori_to_this_res((mmdb::Residue *)reference_residue);

   int istat = 0; // success

   if (!reference_residue) {
      std::cout << "This should not happen!" << std::endl;
      std::cout << "null reference residue in move_std_residue" << std::endl;
   } else {

      if (rtops.size()) { // successful attempt to get the matrix
         mmdb::PPAtom residue_atoms = NULL;
         int nResidueAtoms;
         moving_residue->GetAtomTable(residue_atoms, nResidueAtoms);
         if (nResidueAtoms == 0) {
            std::cout << " something broken in atom residue selection in ";
            std::cout << "mutate, got 0 atoms" << std::endl;
            istat = 0;
         } else {
            istat = 1;

            if (false)
               std::cout << "DEBUG:: move_std_residue: " << nResidueAtoms
                         << " atoms in residue "
                         << moving_residue << " " << moving_residue->seqNum << " "
                         << moving_residue->GetChainID() << std::endl;

            for (int iat=0; iat<nResidueAtoms; iat++) {
               if (residue_atoms[iat]) {
                  clipper::Coord_orth co(residue_atoms[iat]->x,
                                         residue_atoms[iat]->y,
                                         residue_atoms[iat]->z);
                  std::string alt_conf = residue_atoms[iat]->altLoc;

                  std::map<std::string, clipper::RTop_orth>::const_iterator it = rtops.find(alt_conf);

                  if (it != rtops.end()) {
                     clipper::Coord_orth rotted = co.transform(it->second); // an rtop
                     residue_atoms[iat]->x = rotted.x();
                     residue_atoms[iat]->y = rotted.y();
                     residue_atoms[iat]->z = rotted.z();
                  }
               } else {
                  istat = 0;
                  std::cout << "ERROR:: null residue atom in moving residue in move_std_residue: iat: "
                            << iat << std::endl;
               }
            }
         }
      } else {
         istat = 0; // failure
         std::cout << "DISASTER - failed to generate RTop for move_std_residue\n";
         if (reference_residue) {
            //          molecule-class-info.cc:4184: passing `const mmdb::Residue' as `this'
            // argument of `int mmdb::Residue::GetSeqNum ()' discards qualifiers
            mmdb::Residue *tmp = (mmdb::Residue *) reference_residue;
            std::cout << "mainchain atoms missing from residue "
                      << tmp->GetSeqNum()
                      << tmp->GetChainID() << std::endl;
         } else {
            std::cout << "This should not happen!" << std::endl;
            std::cout << "null residue in move_std_residue" << std::endl;
         }
      }
   }
   return istat;
}


// This should be in coot-coord-utils.
//
atom_selection_container_t
coot::add_side_chain_to_terminal_res(atom_selection_container_t asc,
                                     const std::string &res_type,
                                     const std::string &terminus_type,
                                     float b_factor_for_new_atoms,
                                     const coot::protein_geometry &geom) {

   // the add_other_residue_flag is now passed to this function so that, when building
   // forwards, we can add to the penultimate residue rather than the last residue
   // (get_last_residue_in_chain()).

   atom_selection_container_t rasc = asc;
   int istat;
   // molecule_class_info_t molci;

   // molci.install_model(0, asc, geom, "terminal residue", display_in_display_control_widget_status);

   // every (usually 1, occasionally 2) residue in the molecule
   mmdb::Model *model_p = asc.mol->GetModel(1);

   mmdb::Chain *chain;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   if (nchains <= 0) {
      std::cout << "bad nchains in add_cb_to_terminal_res: "
                << nchains << std::endl;
   } else {

      //std::string target_res_type("ALA");
      std::string target_res_type(res_type);
      mmdb::Residue *std_res = util::get_standard_residue_instance(target_res_type);

      if (std_res == NULL) {
         std::cout << "WARNING:: Can't find standard residue for "
                   << target_res_type << "\n";
      } else {

         for (int ichain=0; ichain<nchains; ichain++) {
            chain = model_p->GetChain(ichain);
            if (chain == NULL) {
               // This should not be happen.
               std::cout << "NULL chain in add_cb_to_terminal_res" << std::endl;
            } else {
               bool embed_in_chain_flag = false;
               mmdb::Residue *std_res_copy = coot::deep_copy_this_residue_old_style(std_res, "", 1, -1, embed_in_chain_flag);
               if (std_res_copy) {
                  mmdb::Residue *residue_p = 0;

                  if (terminus_type == "N" || terminus_type == "MN") {
                     residue_p = coot::util::get_first_residue_in_chain(chain);
                  }

                  if (terminus_type == "C" || terminus_type == "MC") {
                     residue_p = coot::util::get_last_residue_in_chain(chain);
                  }

                  if (terminus_type == "singleton")
                     residue_p = coot::util::get_last_residue_in_chain(chain);

                  if (false)
                     std::cout << "------- add_side_chain_to_terminal_res() here with residue_p "
                               << coot::residue_spec_t(residue_p) << std::endl;

                  if (residue_p) {

                     // int istat = molci.move_std_residue(std_res_copy, residue_p);
                     int istat = move_std_residue(std_res_copy, residue_p);

                     if (istat) {

                        mmdb::PPAtom residue_atoms;
                        int nResidueAtoms;
                        residue_p->GetAtomTable(residue_atoms, nResidueAtoms);

                        mmdb::PPAtom std_residue_atoms;
                        int n_std_ResidueAtoms;
                        std_res_copy->GetAtomTable(std_residue_atoms, n_std_ResidueAtoms);

                        // set the b factor for the new atoms.
                        for(int i=0; i<n_std_ResidueAtoms; i++) {
                           std_residue_atoms[i]->tempFactor = b_factor_for_new_atoms;
                        };

                        bool verb = false;
                        if (verb) {
                           std::cout << "Mutate Atom Tables" << std::endl;
                           std::cout << "Before" << std::endl;
                           for(int i=0; i<nResidueAtoms; i++) {
                              std::cout << residue_atoms[i] << std::endl;
                           };

                           std::cout << "To be replaced by:" << std::endl;
                           for(int i=0; i<n_std_ResidueAtoms; i++) {
                              std::cout << std_residue_atoms[i] << std::endl;
                           }
                        }

                        for(int i=0; i<nResidueAtoms; i++) {
                           std::string residue_this_atom (residue_atoms[i]->name);
                           if (residue_this_atom != " O  ")
                              residue_p->DeleteAtom(i);
                        };

                        for(int i=0; i<n_std_ResidueAtoms; i++) {
                           std::string std_residue_this_atom (std_residue_atoms[i]->name);
                           if (std_residue_this_atom != " O  ") {
                              // std::cout << "Adding atom " << std_residue_atoms[i] << std::endl;
                              residue_p->AddAtom(std_residue_atoms[i]);
                           }
                        };
                        // strcpy(residue_p->name, std_res->name);
                        residue_p->TrimAtomTable();

                     }
                     if (true)
                        std::cout << "INFO:: done mutating residue " << coot::residue_spec_t(residue_p)
                                  << " in add_cb_to_terminal_res\n";
                  }
               }

               // this looks a bit leaky!

               asc.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
               asc.mol->FinishStructEdit();

               rasc = make_asc(asc.mol);

//                for (int ii=0; ii<rasc.n_selected_atoms; ii++)
//                   std::cout << "New res with cb: " << rasc.atom_selection[ii] << "\n";
            }
         }
      }
   }
   return rasc;
}

std::pair<int, std::string>
coot::add_terminal_residue(int imol_no, const std::string &terminus_type, mmdb::Residue *residue_p,
                           mmdb::Manager *mol, int udd_atom_index,
                           const std::string &chain_id, const std::string &res_type, float b_factor_for_new_atoms,
                           const clipper::Xmap<float> &xmap, const coot::protein_geometry &geom) {

   int state = 0;
   std::string message;

   bool add_terminal_residue_debug_trials = false;
   int add_terminal_residue_n_phi_psi_trials = 3000;
   bool add_terminal_residue_add_other_residue_flag = false;
   bool is_from_shelx_ins = false;

   if (true) {
      int resno_added = -1; // was unset

      if (terminus_type == "not-terminal-residue") {
         std::string s = "That residue was not at a terminus";
         std::cout << s << std::endl;
         message = s;
      } else {

         mmdb::Residue *upstream_neighbour_residue_p = 0;
         mmdb::Residue *downstream_neighbour_residue_p = 0;
         std::string residue_type_string = res_type;
         int residue_number = residue_p->GetSeqNum();  // bleugh.
         if (residue_type_string == "auto") {
            if (terminus_type == "C" || terminus_type == "MC")
               resno_added = residue_number + 1;
            if (terminus_type == "N" || terminus_type == "MN")
               resno_added = residue_number - 1;
         }

         if (terminus_type == "C") {
            // we have upstream residue
            upstream_neighbour_residue_p = coot::util::previous_residue(residue_p);
         }
         if (terminus_type == "N") {
            // we have downstream residue
            downstream_neighbour_residue_p = coot::util::next_residue(residue_p);
         }

         float bf = b_factor_for_new_atoms;
         coot::residue_by_phi_psi addres(terminus_type, residue_p, chain_id, res_type, bf);

         if (upstream_neighbour_residue_p)
            addres.set_upstream_neighbour(upstream_neighbour_residue_p);
         if (downstream_neighbour_residue_p)
            addres.set_downstream_neighbour(downstream_neighbour_residue_p);

#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
         unsigned int n_threads = coot::get_max_number_of_threads();
         if (n_threads >= 1)
            addres.thread_pool(&static_thread_pool, n_threads);
#endif

         // std::cout << "DEBUG:: term_type: " << terminus_type << std::endl;

         // This was for debugging, so that we get *some* solutions at least.
         //
         addres.set_acceptable_fit_fraction(0.0); //  the default is 0.5, I think

         // do we want to output the trial solutions as pdbs?

         // std::cout << "--------------------- here with add_terminal_residue_debug_trials "
         // << add_terminal_residue_debug_trials << std::endl;

         if (add_terminal_residue_debug_trials)
            addres.write_trial_pdbs();

         // map value over protein, stops rigid body refinement down
         // into previous residue?  Yes, but only after I altered the
         // scoring in ligand to use this masked map (and not the
         // pristine_map as it has previously done).
         //
         // The mask value used to be -2.0, but what happened is that
         // in low resolution maps, the N atom of the fragment "felt"
         // some of the -2.0 because of interpolation/cubic-spline.
         // So -2.0 seems too much for low res.  Let's try -1.0.
         //
         float masked_map_val = -1.0;
         addres.set_map_atom_mask_radius(1.2);
         if (terminus_type == "MC" || terminus_type == "MN" || terminus_type == "singleton")
            masked_map_val = 0.0;
         addres.set_masked_map_value(masked_map_val);

         // addres.import_map_from(molecules[imol_map].xmap,
         // molecules[imol_map].map_sigma());

         // This masked map will be the one that is used for rigid
         // body refinement, unlike normal ligand class usage which
         // uses the xmap_pristine.
         //

         // old style mask by all the protein atoms - slow? (especially on sgis?)
         //          short int mask_waters_flag = 1;
         //      addres.mask_map(molecules[imol].atom_sel.mol, mask_waters_flag);
         //
         mmdb::PPAtom atom_sel = NULL;
         int n_selected_atoms = 0;
         mmdb::realtype radius = 8.0;  // more than enough for 2 residue mainchains.
         int SelHndSphere = mol->NewSelection();
         mmdb::Atom *terminal_at = NULL;
         std::string atom_name = "Unassigned";
         if (terminus_type == "MC" || terminus_type == "C" || terminus_type == "singleton")
            atom_name = " C  ";
         if (terminus_type == "MN" || terminus_type == "N")
            atom_name = " N  ";
         if (atom_name != "Unassigned") {
            mmdb::PPAtom residue_atoms;
            int nResidueAtoms;
            mmdb::Residue *res_tmp_p = residue_p;
            res_tmp_p->GetAtomTable(residue_atoms, nResidueAtoms);
            for (int i=0; i<nResidueAtoms; i++)
               if (atom_name == residue_atoms[i]->name) {
                  terminal_at = residue_atoms[i];
                  break;
               }

            if (terminal_at) {
               mol->SelectSphere(SelHndSphere, mmdb::STYPE_ATOM,
                                 terminal_at->x,
                                 terminal_at->y,
                                 terminal_at->z,
                                 radius, mmdb::SKEY_NEW);
               mol->GetSelIndex(SelHndSphere, atom_sel, n_selected_atoms);

               // Why is this commented out?
               // Because the map is added to addres after this point. Hmm. Masking is a good idea.
               // Maybe copy the molecules[imol_map].xmap and mask it making xmap_masked
               // and pass xmap_masked to the best_fit_phi_psi()?
               //
               // int invert_flag = 0;
               // addres.mask_map(molecules[imol].atom_sel.mol, SelHndSphere, invert_flag);

               mol->DeleteSelection(SelHndSphere);
            }
         } else {
            std::cout << "WARNING:: terminal atom not assigned - no masking!" << std::endl;
         }

         std::cout << "INFO:: fitting terminal residue with "
                   << add_terminal_residue_n_phi_psi_trials << " random trials"
                   << std::endl;

           coot::minimol::molecule mmol =
            addres.best_fit_phi_psi(add_terminal_residue_n_phi_psi_trials, false,
                                    add_terminal_residue_add_other_residue_flag,
                                    xmap);

         std::vector<coot::minimol::atom *> mmatoms = mmol.select_atoms_serial();

         if (mmol.is_empty()) {

            // this should not happen:
            std::cout <<  "WARNING: ------------- empty molecule: "
                      << "failed to find a fit for terminal residue"
                      << std::endl;

         } else {

            // check that we are adding some atoms:
            //
            if (mmatoms.size() == 0) {
               std::cout << "WARNING: failed to find a fit for terminal residue"
                         << std::endl;
               // if (use_graphics_interface_flag) {
               // // GtkWidget *w = create_add_terminal_residue_finds_none_dialog();
               // GtkWidget *w = widget_from_builder("add_terminal_residue_finds_none_dialog");
               // gtk_widget_show(w);
               // }

            } else {

               state = 1;

               atom_selection_container_t terminal_res_asc;

               // if this is begin added to a shelx molecule, then we
               // need to set the occs to 11.0
               //
               if (is_from_shelx_ins)
                  bf = 11.0;

               if (add_terminal_residue_add_other_residue_flag) {

                  // check that the other residue is not in the molecule before adding
                  // all of mmol. If it is already in the molecule, remove it from mmol
                  if (terminus_type == "C" || terminus_type == "MC") {
                     coot::residue_spec_t other_residue_spec(chain_id, resno_added+1, "");
                     mmdb::Residue *res_other = coot::util::get_residue(other_residue_spec, mol);
                     if (res_other) {
                        mmol[0][other_residue_spec.res_no].atoms.clear();
                     }
                  } else {
                     if (terminus_type == "N" || terminus_type == "MN") {
                        coot::residue_spec_t other_residue_spec(chain_id, resno_added-1, "");
                        mmdb::Residue *res_other = coot::util::get_residue(other_residue_spec, mol);
                        if (res_other) {
                           mmol[0][other_residue_spec.res_no].atoms.clear();
                        }
                     }
                  }
               }
               terminal_res_asc.mol = mmol.pcmmdbmanager();
               // terminal_res_asc.mol->WritePDBASCII("terminal_res_asc.pdb");

               int SelHnd = terminal_res_asc.mol->NewSelection();
               terminal_res_asc.mol->SelectAtoms(SelHnd, 0, "*",
                                                 mmdb::ANY_RES, // starting resno, an int
                                                 "*", // any insertion code
                                                 mmdb::ANY_RES, // ending resno
                                                 "*", // ending insertion code
                                                 "*", // any residue name
                                                 "*", // atom name
                                                 "*", // elements
                                                 "*"  // alt loc.
                                                 );
               terminal_res_asc.mol->GetSelIndex(SelHnd,
                                                 terminal_res_asc.atom_selection,
                                                 terminal_res_asc.n_selected_atoms);

               // Now we add in the cb of this residue (currently it
               // only has main chain atoms). This is somewhat
               // involved - the methods to manipulate the standard
               // residues are part of molecule_class_info_t - so we
               // need to make an instance of that class.
               //

                atom_selection_container_t tmp_asc =
                   add_side_chain_to_terminal_res(terminal_res_asc, res_type, terminus_type, b_factor_for_new_atoms, geom);


                // std::cout << "-------------- tmp_asc --------" << std::endl;
                // debug_atom_selection_container(tmp_asc);

               // If this is wrong also consider fixing execute_rigid_body_refine()
               //
//                std::cout << "debug: add_residue asc has n_selected_atoms = "
//                          << terminal_res_asc.n_selected_atoms << " "
//                          << terminal_res_asc.atom_selection << std::endl;

//                for (int i=0; i< terminal_res_asc.n_selected_atoms; i++) {
//                   std::cout << "debug: add_residue asc has chain_id: "
//                             << terminal_res_asc.atom_selection[i]->GetChainID()
//                             << " for " << terminal_res_asc.atom_selection[i]
//                             << std::endl;
//                }

               coot::residue_spec_t rs(residue_p);
               remove_TER_internal(mol, residue_p);

               if (false) {
#if 0
                  make_moving_atoms_graphics_object(imol, tmp_asc);
                  moving_atoms_asc_type = coot::NEW_COORDS_INSERT;
                  graphics_draw();
                  coot::refinement_results_t dummy;
                  if (use_graphics_interface_flag) {
                     do_accept_reject_dialog("Terminal Residue", dummy);
                  }
#endif
               } else {

                  if (is_from_shelx_ins) {
                     for (int i=0; i<tmp_asc.n_selected_atoms; i++) {
                        tmp_asc.atom_selection[i]->occupancy = 11.0;
                     }
                  }

                  insert_coords(mol, udd_atom_index, tmp_asc);

                  // when we place a new residue at the C-terminus, the oxygen
                  // position of this resiude (which is ignored in the selection
                  // of the position of the next residue) can well be in the
                  // wrong place!
                  // add_res knows the position of the residue being added here
                  // so we can use that to tell us where to place the O oxygen
                  // of the current residue.
                  //
                  // 201805014-PE merging - Oh, I've done it twice (forgetten first)
                  // by different methods - heyho
                  if (terminus_type == "C" || terminus_type == "MC") {
                     clipper::Coord_orth new_o_pos = addres.best_fit_phi_psi_attaching_oxygen_position_update(mmol, residue_p);
                     move_atom(" O  ", residue_p, new_o_pos);
                  }
                  // method from master:
                  // if (terminus_type == "C" || terminus_type == "MC")
                  //    molecules[imol_moving_atoms].move_O_atom_of_added_to_residue(res_p, chain_id);
                  // graphics_draw();
               }
            }
         }
      }
   }
   return std::pair<int, std::string>(state, message);
}

