/*
 * api/coot-molecule-change-chain-id.cc
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


#include "coot-molecule.hh"

// Change chain id
// return -1 on a conflict
// 1 on good.
// 0 on did nothing
//
std::pair<int, std::string>
coot::molecule_t::change_chain_id(const std::string &from_chain_id,
                                  const std::string &to_chain_id,
                                  bool use_resno_range,
                                  int start_resno, int end_resno) {



   // Return a copy of the pointer to the chain (only).  Return NULL
   // on chain with given chain ID not found.
   //
   auto get_chain = [] (const std::string &chain_id, mmdb::Manager *mol) {

      mmdb::Chain *r = NULL;
      if (mol) {
         int imod = 1;
         mmdb::Model *model_p = mol->GetModel(imod);
         mmdb::Chain *chain_p;
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            chain_p = model_p->GetChain(ichain);
            std::string mol_chain_id = chain_p->GetChainID();
            if (chain_id == mol_chain_id) {
               r = chain_p;
               break;
            }
         }
      }
      return r;
   };



   int istat = 0; // done nothing to start with
   std::string message("Nothing to say");

   //    std::cout << "DEBUG:: use_resno_range: " << use_resno_range << std::endl;

   if (atom_sel.n_selected_atoms > 0) {

      if (use_resno_range) {

         std::pair<int, std::string> r =
            change_chain_id_with_residue_range(from_chain_id, to_chain_id, start_resno, end_resno);
         istat = r.first;
         message = r.second;

      } else {
      // The usual case, I imagine

         bool target_chain_id_exists = false;

         int n_models = atom_sel.mol->GetNumberOfModels();
         for (int imod=1; imod<=n_models; imod++) {

            mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
            if (! model_p) continue;
            // run over chains of the existing mol
            int nchains = model_p->GetNumberOfChains();
            if (nchains <= 0) {
               std::cout << "bad nchains in molecule " << nchains
                         << std::endl;
            } else {
               for (int ichain=0; ichain<nchains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  if (chain_p == NULL) {
                     // This should not be necessary. It seem to be a
                     // result of mmdb corruption elsewhere - possibly
                     // DeleteChain in update_molecule_to().
                     std::cout << "NULL chain in change chain id" << std::endl;
                  } else {
                     std::string chain_id = chain_p->GetChainID();
                     if (to_chain_id == chain_id) {
                        target_chain_id_exists = true;
                        break;
                     }
                  }
               }
            }
         }

         if (!target_chain_id_exists) {

            n_models = atom_sel.mol->GetNumberOfModels();
            for (int imod=1; imod<=n_models; imod++) {

               mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
               if (! model_p) continue;
               // run over chains of the existing mol
               int nchains = model_p->GetNumberOfChains();
               if (nchains <= 0) {
                  std::cout << "bad nchains in molecule " << nchains
                            << std::endl;
               } else {
                  for (int ichain=0; ichain<nchains; ichain++) {
                     mmdb::Chain *chain_p = model_p->GetChain(ichain);
                     if (chain_p) {
                        std::string chain_id = chain_p->GetChainID();
                        if (from_chain_id == chain_id) {
                           make_backup("change chain id");
                           chain_p->SetChainID(to_chain_id.c_str());
                           // change the links here

                           int n_changed = coot::util::change_chain_in_links(model_p, from_chain_id, to_chain_id);
                           istat = 1;
                           // have_unsaved_changes_flag = 1;
                           atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
                           atom_sel.mol->FinishStructEdit();
                           atom_sel = make_asc(atom_sel.mol);
                           // make_bonds_type_checked(__FUNCTION__);
                        }
                     }
                  }
               }
            }

         } else {

            // OK, can we do a merge? (do we have non-overlapping residue ranges?)
            // If so, so a change_chain_id_with_residue_range()
            //
            mmdb::Chain *chain_p_from = get_chain(from_chain_id, atom_sel.mol);
            mmdb::Chain *chain_p_to   = get_chain(to_chain_id, atom_sel.mol);
            bool done_merge = false;
            if (chain_p_from) {
               if (chain_p_to) {
                  std::pair<bool, int> min_r_1 = coot::util::min_resno_in_chain(chain_p_from);
                  std::pair<bool, int> max_r_1 = coot::util::max_resno_in_chain(chain_p_from);

                  if (false) {
                     std::cout << "--------- here with min_r_1  " << min_r_1.first << " " << min_r_1.second << std::endl;
                     std::cout << "--------- here with max_r_1  " << max_r_1.first << " " << max_r_1.second << std::endl;
                     std::cout << "--------- here with from_chain_id " << from_chain_id << std::endl;
                     std::cout << "--------- here with to_chain_id " << to_chain_id << std::endl;
                  }

                  if (min_r_1.first) {
                     if (max_r_1.first) {
                        start_resno = min_r_1.second;
                        end_resno = max_r_1.second;
                        std::pair<int, std::string> r =
                           change_chain_id_with_residue_range(from_chain_id, to_chain_id, start_resno, end_resno);
                        istat = r.first;
                        message = r.second;
                     }
                  }
               }
            }

            if (! done_merge) {
               std::cout << "WARNING:: CONFLICT: target chain id " << to_chain_id << " already exists "
                         << "in this molecule" << std::endl;
               message = "WARNING:: CONFLICT: target chain id (";
               message += to_chain_id;
               message += ") already \nexists in this molecule!";
            }
         }
      } // residue range
   } // no atoms

   return std::pair<int, std::string> (istat, message);
}


std::pair<int, std::string>
coot::molecule_t::change_chain_id_with_residue_range(const std::string &from_chain_id,
                                                     const std::string &to_chain_id,
                                                     int start_resno,
                                                     int end_resno) {

   if (false) {
      std::cout << "-------------------- change_chain_id_with_residue_range ---- " << std::endl;
      std::cout << "-------------------- change_chain_id_with_residue_range ---- " << start_resno << std::endl;
      std::cout << "-------------------- change_chain_id_with_residue_range ---- " <<   end_resno << std::endl;
   }

   std::string message;
   int istat = 0;

   short int target_chain_id_exists = 0;
   int n_models = atom_sel.mol->GetNumberOfModels();
   for (int imod=1; imod<=n_models; imod++) {

      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      if (! model_p) continue;

      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      if (nchains <= 0) {
         std::cout << "bad nchains in molecule " << nchains
                   << std::endl;
      } else {
         for (int ichain=0; ichain<nchains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            if (chain_p == NULL) {
               // This should not be necessary. It seem to be a
               // result of mmdb corruption elsewhere - possibly
               // DeleteChain in update_molecule_to().
               std::cout << "NULL chain in change chain id" << std::endl;
            } else {
               std::string chain_id = chain_p->GetChainID();
               if (to_chain_id == chain_id) {
                  target_chain_id_exists = 1;
                  break;
               }
            }
         }
      }
   }

   if (target_chain_id_exists == 0) {

      // So we are moving residues 12->24 of Chain A to (new) Chain
      // C.  Not very frequent, I suspect.

      // make sure start and end are a sensible way round
      if (end_resno < start_resno) {
         int tmp = end_resno;
         end_resno = start_resno;
         start_resno = tmp;
      }

      for (int imod=1; imod<=n_models; imod++) {

         mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
         if (! model_p) continue;
         // run over chains of the existing mol
         int nchains = model_p->GetNumberOfChains();
         if (nchains <= 0) {
            std::cout << "bad nchains in molecule " << nchains << std::endl;
         } else {
            for (int ichain=0; ichain<nchains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               if (chain_p) {
                  std::string chain_id = chain_p->GetChainID();
                  if (from_chain_id == chain_id) {

                     std::cout << "matched from_chain_id " << from_chain_id << std::endl;

                     // So we have the chain from which we wish to move residues
                     make_backup("change chain id with residue range");
                     atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);
                     // Create a new chain and add it to the molecule
                     mmdb::Chain *new_chain_p = new mmdb::Chain;
                     new_chain_p->SetChainID(to_chain_id.c_str());
                     model_p->AddChain(new_chain_p);

                     int nresidues = chain_p->GetNumberOfResidues();
                     unsigned int n_changed = 0;
                     for (int ires=0; ires<nresidues; ires++) {
                        mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                        if (residue_p->GetSeqNum() >= start_resno &&
                            residue_p->GetSeqNum() <= end_resno) {
                           int iseqnum  = residue_p->GetSeqNum();
                           mmdb::pstr inscode = residue_p->GetInsCode();
                           mmdb::Residue *residue_copy = coot::util::deep_copy_this_residue(residue_p);
                           chain_p->DeleteResidue(iseqnum, inscode);
                           new_chain_p->AddResidue(residue_copy);
                           n_changed++;
                        }
                     }

                     if (n_changed > 0) {
                        istat = 1;
                        // have_unsaved_changes_flag = 1;
                        atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
                        atom_sel.mol->FinishStructEdit();
                        atom_sel = make_asc(atom_sel.mol);
                        // make_bonds_type_checked(__FUNCTION__);
                     }
                  }
               }
            }
         }
      }
   } else {

      // target chain alread exists.   Here is where we merge...

      // We need to check that we are not reproducing residues that
      // already exist in the chain.  If we are doing that, we stop
      // and give an error message back telling user that that
      // residues exists already.
      //
      for (int imod=1; imod<=n_models; imod++) {

         mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
         if (! model_p) continue;

         // run over chains of the existing mol
         int nchains = model_p->GetNumberOfChains();
         short int residue_already_exists_flag = 0;
         if (nchains <= 0) {
            std::cout << "bad nchains in molecule " << nchains
                      << std::endl;
         } else {
            for (int ichain=0; ichain<nchains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               if (chain_p) {
                  std::string chain_id = chain_p->GetChainID();
                  if (to_chain_id == chain_id) {

                     int nresidues = chain_p->GetNumberOfResidues();
                     int existing_residue_number = 0; // set later
                     for (int ires=0; ires<nresidues; ires++) {
                        mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                        if (residue_p->GetSeqNum() >= start_resno &&
                            residue_p->GetSeqNum() <= end_resno) {
                           residue_already_exists_flag = 1;
                           existing_residue_number = residue_p->GetSeqNum();
                           break;
                        }

                     }
                     if (residue_already_exists_flag) {
                        message += "CONFLICT!  Residue ";
                        message += coot::util::int_to_string(existing_residue_number);
                        message += " already exists in chain ";
                        message += chain_id;
                        message += "\nChange chain ID of residue range failed.\n";
                     } else {

                        // We are OK to move the residue into the existing chain
                        // (move is done by copy and delete)

                        mmdb::Chain *to_chain = NULL;
                        for (int ichain=0; ichain<nchains; ichain++) {
                           chain_p = model_p->GetChain(ichain);
                           if (chain_p) {
                              std::string chain_id = chain_p->GetChainID();
                              if (to_chain_id == chain_id) {
                                 to_chain = chain_p;
                              }
                           }
                        }


                        if (to_chain) {
                           make_backup("change chain id with residue-range 2");
                           for (int ichain=0; ichain<nchains; ichain++) {
                              chain_p = model_p->GetChain(ichain);
                              int to_chain_nresidues = chain_p->GetNumberOfResidues();
                              if (chain_p) {
                                 std::string chain_id = chain_p->GetChainID();
                                 if (from_chain_id == chain_id) {
                                    for (int ires=0; ires<to_chain_nresidues; ires++) {
                                       mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                                       if (residue_p->GetSeqNum() >= start_resno &&
                                           residue_p->GetSeqNum() <= end_resno) {

                                          int iseqnum  = residue_p->GetSeqNum();
                                          mmdb::pstr inscode = residue_p->GetInsCode();
                                          mmdb::Residue *residue_copy =
                                             coot::util::deep_copy_this_residue_add_chain(residue_p, "", 1, 1);
                                          // delete the residue in the "from" chain:
                                          chain_p->DeleteResidue(iseqnum, inscode);

                                          //
                                          change_chain_id_with_residue_range_helper_insert_or_add(to_chain, residue_copy);
                                          // this is done by the deep_copy
                                          // to_chain->AddResidue(residue_copy);
                                       }
                                    }
                                 }
                              }
                           }
                           istat = 1;
                           // have_unsaved_changes_flag = 1;
                           atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
                           atom_sel.mol->FinishStructEdit();
                           atom_sel = make_asc(atom_sel.mol);
                           // make_bonds_type_checked(__FUNCTION__);
                        }
                     }
                  }
               }
               if (residue_already_exists_flag)
                  break;
            }
         }
         if (residue_already_exists_flag)
            break;
      }
   }

   return std::pair<int, std::string> (istat, message);

}

void
coot::molecule_t::change_chain_id_with_residue_range_helper_insert_or_add(mmdb::Chain *to_chain_p, mmdb::Residue *new_residue) {

   // OK, if we can, let's try to *insert* the new_residue into the
   // right place in the to_chain_p.  If we don't manage to insert it,
   // let's fall back and simply add it.  The new residue is inserted
   // *before* the residue specified by the given serial number.

   // Let's use the serial number interface to InsResidue()

   int resno_new_residue = new_residue->GetSeqNum();
   std::string ins_code = new_residue->GetInsCode();
   int target_res_serial_number = coot::RESIDUE_NUMBER_UNSET;
   int target_res_seq_num = resno_new_residue; // simply ignore ins codes :)
   std::string target_res_ins_code = ""; // ignore this.  Fix later.
   int best_seq_num_diff = 99999999;

   int n_chain_residues;
   mmdb::PResidue *chain_residues;
   to_chain_p->GetResidueTable(chain_residues, n_chain_residues);
   for (int iserial=0; iserial<n_chain_residues; iserial++) {
      int chain_residue_seq_num = chain_residues[iserial]->GetSeqNum();
      int this_seq_num_diff = chain_residue_seq_num - resno_new_residue;
      if (this_seq_num_diff > 0) {
         if (this_seq_num_diff < best_seq_num_diff) {
            best_seq_num_diff = this_seq_num_diff;
            target_res_serial_number = iserial;
         }
      }
   }

   if (target_res_serial_number != coot::RESIDUE_NUMBER_UNSET) {
      // Good stuff
      // std::cout << "Debugging, inserting residue here...." << std::endl;
      to_chain_p->InsResidue(new_residue, target_res_serial_number);
   } else {
      // std::cout << "Debugging, adding residue here...." << std::endl;
      to_chain_p->AddResidue(new_residue);
   }
}
