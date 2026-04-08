/*
 * api/coot-molecule-replace-fragment.cc
 * 
 * Copyright 2020 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
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


#include "geometry/residue-and-atom-specs.hh"
#include "coot-molecule.hh"

// We need to find the serial number of the residue after the residue
// we want to insert (i.e. the new residue will be inserted just
// before the residue whose serial number we return).
//
// return -1 on error.
std::pair<int, mmdb::Residue *>
coot::molecule_t::find_serial_number_for_insert(int seqnum_for_new,
                                                const std::string &ins_code_for_new,
                                                const std::string &chain_id) const {

   int iserial_no = -1;
   std::pair<int, std::string> current_diff(999999, "");
   int n_chains = atom_sel.mol->GetNumberOfChains(1);
   mmdb::Residue *res = NULL;

   for (int i_chain=0; i_chain<n_chains; i_chain++) {

      mmdb::Chain *chain_p = atom_sel.mol->GetChain(1,i_chain);

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


// Replace the atoms in this molecule by those in the given atom selection.
int
coot::molecule_t::replace_fragment(atom_selection_container_t asc) {

   auto get_chain =  [] (const std::string& chain_id) {
      mmdb::Chain *chain_p = nullptr;
      return chain_p;
   };

   if (! asc.mol) return 0;

   bool move_zero_occ = true;

   // replace an atom if you can, otherwise create a new atom (and a new residue and chain if needed)

   make_backup("replace_fragment");

   for (int i=0; i<asc.n_selected_atoms; i++) {
      int idx = -1;
      mmdb::Atom *at = asc.atom_selection[i];

      // std::cout << "debug:: in replace_fragment() with asc.UDDOldAtomIndexHandle "
      //           << asc.UDDOldAtomIndexHandle << " " << coot::atom_spec_t(at) << std::endl;

      if (! at->isTer()) {
         // can we find the atom with fast indexing?
         if (asc.UDDOldAtomIndexHandle >= 0) {
            // OK for fast atom indexing
            int ref_index = -1;
            if (at->GetUDData(asc.UDDOldAtomIndexHandle, ref_index) == mmdb::UDDATA_Ok) {
               if (ref_index >= 0) {
                  if (moving_atom_matches(at, ref_index)) {
                     idx = ref_index; // yay.
                     std::cout << "yay " << idx << std::endl;
                  }
               }
            }
         }

         // std::cout << "here in replace_fragment() with idx " << idx << " for " << coot::atom_spec_t(at)
         //           << std::endl;

         if (idx == -1) {
            idx = full_atom_spec_to_atom_index(coot::atom_spec_t(at));
         }

         if (idx != -1) {
            mmdb::Atom *ref_atom = atom_sel.atom_selection[idx];
            ref_atom->x = at->x;
            ref_atom->y = at->y;
            ref_atom->z = at->z;

         } else {

            // add the atom
            mmdb::Chain *chain_p = get_chain(at->GetChainID());
            mmdb::Residue *residue_p = get_residue(coot::residue_spec_t(coot::atom_spec_t(at)));

            if (! chain_p) {
               int imod = 1;
               mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
               if (model_p) {
                  mmdb::Chain *chain_p = new mmdb::Chain;
                  chain_p->SetChainID(at->GetChainID());
                  residue_p = new mmdb::Residue;
                  residue_p->seqNum = at->GetSeqNum();
                  residue_p->SetResName(at->residue->GetResName());
                  chain_p->AddResidue(residue_p);
                  model_p->AddChain(chain_p);
                  atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
                  atom_sel.mol->FinishStructEdit();
               }
            } else {
               if (residue_p) {
                  // std::cout << "   ======= found the residue " << std::endl;
               } else {
                  residue_p = new mmdb::Residue;
                  residue_p->SetResID(at->residue->GetResName(), at->residue->seqNum, at->residue->insCode);
                  int res_no = at->GetSeqNum();
                  std::string ins_code(at->GetInsCode());
                  std::pair<int, mmdb::Residue *> sn =
                     find_serial_number_for_insert(res_no, ins_code, chain_p->GetChainID());

                  if (sn.first != -1) { // normal insert

                     int n_residues_before = chain_p->GetNumberOfResidues();
                     int n_chain_residues = chain_p->InsResidue(residue_p, sn.first);
                     mmdb::Residue *res_after_p = get_residue(coot::residue_spec_t(coot::atom_spec_t(at)));

                  } else {
                     chain_p->AddResidue(residue_p);
                     atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);// needed?
                     atom_sel.mol->FinishStructEdit();
                  }
               }
            }

            if (residue_p) {
               mmdb::Atom *at_copy(at);
               residue_p->AddAtom(at_copy);
               // residue_p->TrimAtomTable();
            }
         }
      }
   }

   atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);
   atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);// needed?
   coot::util::pdbcleanup_serial_residue_numbers(atom_sel.mol);
   atom_sel.mol->FinishStructEdit();

   atom_sel = make_asc(atom_sel.mol);

   // save_info.new_modification("replace_fragment");
   // have_unsaved_changes_flag = 1;
   // if (show_symmetry)
   //    update_symmetry();
   // make_bonds_type_checked();
   return 1;
}


// replace the atoms of SelHnd, which is a selection of mol_ref into this molecule.
// Use old_atom_index_handle for fast indexing.
int
coot::molecule_t::replace_fragment(mmdb::Manager *mol_ref, int old_atom_index_handle, int SelHnd) {

   auto tp_0 = std::chrono::high_resolution_clock::now();
   int status = 0;
   mmdb::Atom **selection_atoms = 0;
   int n_selection_atoms = 0;
   mol_ref->GetSelIndex(SelHnd, selection_atoms, n_selection_atoms);

   for (int i=0; i<n_selection_atoms; i++) {
      mmdb:: Atom *at_frag = selection_atoms[i];
      // std::cout << "### replace this: " << atom_spec_t(at_frag) << std::endl;
      int idx = -1;
      int ierr = at_frag->GetUDData(old_atom_index_handle, idx);
      if (ierr == mmdb::UDDATA_Ok) {
         mmdb::Atom *at = atom_sel.atom_selection[idx];
         // std::cout << "replacing position of " << atom_spec_t(at) << " from " << atom_spec_t(at_frag) << std::endl;
         at->x = at_frag->x;
         at->y = at_frag->y;
         at->z = at_frag->z;
         status = 1; // at least one atom was found
      } else {
         std::cout << "GetUDData() failed for " << coot::atom_spec_t(at_frag) << " " << old_atom_index_handle << std::endl;
      }
   }
   auto tp_1 = std::chrono::high_resolution_clock::now();
   auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_1 - tp_0).count();

   // std::cout << "debug (new) replace_fragment() took " << d10 << " microseconds" << std::endl;

   return status;
}
