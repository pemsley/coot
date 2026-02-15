/* src/molecule-class-info.cc
 *
 * Copyright 2002, 2003, 2004, 2005, 2006 by The University of York
 * Author: Paul Emsley
 * Copyright 2007 by Paul Emsley
 * Copyright 2013, 2014 by Medical Research Council
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
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"

#ifdef _MSC_VER
#include <windows.h>
#endif

#include <stdlib.h>
#include <string.h>  // strncpy

#include <stdexcept>

#include "clipper/core/xmap.h"
#include "density-contour/CIsoSurface.h"

#include "clipper/core/hkl_compute.h"
#include "clipper/clipper-phs.h"
#include "clipper/core/map_utils.h" // Map_stats


#include <mmdb2/mmdb_manager.h>
// #include "mmdb-extras.h"
// #include "mmdb.h"
// #include "mmdb-crystal.h"

#include "graphics-info.h"

// #include "xmap-utils.h"
// #include "coot-coord-utils.hh"

#include "geometry/mol-utils.hh"
#include "coot-utils/coot-map-utils.hh"
#include "rotamer-search-modes.hh"

#include "molecule-class-info.h"
#include <mmdb2/mmdb_math_align.h>
#include <mmdb2/mmdb_tables.h>

#include "utils/logging.hh"
extern logging logger;

int
molecule_class_info_t::mutate(int resno, const std::string &insertion_code,
			      const std::string &chain_id,
			      const std::string &residue_type) {

   int istat = -1;
   mmdb::Residue *res;

   if (atom_sel.mol) {
      mmdb::Residue *res;

      int nSelResidues;
      mmdb::PPResidue SelResidues;
      int SelHnd = atom_sel.mol->NewSelection();
      atom_sel.mol->Select(SelHnd, mmdb::STYPE_RESIDUE, 1,
			   chain_id.c_str(),
			   resno, insertion_code.c_str(),
			   resno, insertion_code.c_str(),
			   "*", "*", "*", "*",
			   mmdb::SKEY_NEW);

      atom_sel.mol->GetSelIndex ( SelHnd, SelResidues,nSelResidues );

      if (nSelResidues < 1) {
	 std::cout << "WARNING:: Can't find residue (mutate) "
		   << resno << " " << insertion_code << " chain-id \""
		   << chain_id << "\"\n";
      } else {
	 res = SelResidues[0];
	 istat = mutate(res, residue_type); // 1 is good.
      }
      // atom_sel.mol->DeleteSelection(SelHnd);
   } else {
      std::cout << "ERROR:: Null mol" << std::endl;
   }
   return istat;
}

// this is the interface from the GUI.
int
molecule_class_info_t::mutate(int atom_index, const std::string &residue_type,
			      short int make_stub_flag) {

   int r = -1;

   if (false) {
      // 20260119-PE this function is crashing. What's happend to the molecule?
      //
      // the problem was mixing old thinking about transfer of atom and residue info
      // into the new thinking.
      std::cout << "::::::::::::::: atom_index " << atom_index << std::endl;
      for (int iat=0; iat<atom_sel.n_selected_atoms; iat++) {
         mmdb::Atom *at = atom_sel.atom_selection[iat];
         std::cout << "atom " << at->residue->GetSeqNum() << " " << at->GetAtomName() << std::endl;
      }
   }

   if (atom_index < 0) return r;

   if (atom_index < atom_sel.n_selected_atoms) {
      mmdb:: Atom *at = atom_sel.atom_selection[atom_index];
      if (at) {
         mmdb::Residue *res = at->residue;
         if (res) {
            r = mutate(res, residue_type);
            if (atom_sel.mol) {
               if (make_stub_flag) {
                  int resno = res->GetSeqNum();
                  std::string chain_id(res->GetChainID());
                  std::string inscode(res->GetInsCode());
                  // BL says:: we made a new CB with a new B-factor
                  // we should either not delete it (pass do stub_flag) or
                  // save the b_factor and apply again. Not sure whats preferred,
                  // so leave it for now (shouldnt do stubbing anyway...?!)
                  // FIXME
                  delete_residue_sidechain(chain_id, resno, inscode);
               }
            }
         }
      }
   }
   return r;
}

// return 0 on fail.
//
// verbose_mode is an optional argument - default true
int
molecule_class_info_t::mutate(mmdb::Residue *res, const std::string &residue_type, bool verbose_mode) {

   graphics_info_t g;
   int istate = 0;

   if (verbose_mode)
      std::cout << "INFO:: mutate " << res->GetSeqNum() << " "
                << res->GetChainID() << " to a " << residue_type
                << std::endl;

   // get the standard orientation residue for this residue type
   mmdb::PPResidue SelResidue;
   int nSelResidues;

   if (g.standard_residues_asc.n_selected_atoms == 0) {
      std::cout << "WARNING:: 0 standard atoms selected in mutate" << std::endl
                << "WARNING:: did you fail to read the standard residues "
                << "correctly?" << std::endl;
      return 0;
   } else {
      if (g.standard_residues_asc.mol == NULL) {
         std::cout << "WARNING:: null standard_residues_asc in mutate" << std::endl
                   << "WARNING:: did you fail to read the standard residues "
                   << "correctly   ?" << std::endl;
         return 0;
      }
   }

   // std::cout << "--------------------- searching for type " << residue_type << " in standard residues mol "
   // << g.standard_residues_asc.mol << std::endl;

   int selHnd = g.standard_residues_asc.mol->NewSelection(); // deleted at end.
   g.standard_residues_asc.mol->Select (selHnd,mmdb::STYPE_RESIDUE, 1, // .. TYPE, iModel
					 "*", // Chain(s) it's "A" in this case.
					 mmdb::ANY_RES,"*",  // starting res
					 mmdb::ANY_RES,"*",  // ending res
					 residue_type.c_str(),  // residue name
					 "*",  // Residue must contain this atom name?
					 "*",  // Residue must contain this Element?
					 "*",  // altLocs
					 mmdb::SKEY_NEW // selection key
                                        );

   g.standard_residues_asc.mol->GetSelIndex ( selHnd, SelResidue,nSelResidues );

   if (nSelResidues != 1) {

      std::cout << "ERROR:: This should never happen - ";
      std::cout << "badness in mci:mutate() standard residue selection\n";

   } else {

      // get_ori_to_this_res() is for orienting a peptide based on the main-chain atoms.

      std::map<std::string, clipper::RTop_orth> rtops =
         coot::util::get_ori_to_this_res(res); // the passed res.

      if (rtops.size() == 0) {

         std::cout << "ERROR::: failure to get orientation matrix" << std::endl;

      } else {

         // we need to generate a std_residue for each alt conf,
         // (usually, just one of course).
         std::map<std::string, clipper::RTop_orth>::const_iterator it;

         for (it=rtops.begin(); it!=rtops.end(); ++it) {
            bool whole_res = true;
            bool embed_in_new_chain = false; // I think
            mmdb::Residue *std_residue = coot::deep_copy_this_residue_old_style(SelResidue[0], "", whole_res,
            atom_sel.UDDAtomIndexHandle, embed_in_new_chain);
            if (! std_residue) {

               std::cout << "ERROR:: failure to get std_residue in mutate()" << std::endl;

            } else {

               make_backup(__FUNCTION__);

               mmdb::PPAtom residue_atoms;
               int nResidueAtoms;
               std_residue->GetAtomTable(residue_atoms, nResidueAtoms);
               if (nResidueAtoms == 0) {
                  std::cout << "ERROR:: something broken in atom residue selection in ";
                  std::cout << "mutate, got 0 atoms" << std::endl;
               } else {
                  for(int iat=0; iat<nResidueAtoms; iat++) {

                     clipper::Coord_orth co(residue_atoms[iat]->x,
                        residue_atoms[iat]->y,
                        residue_atoms[iat]->z);
                        clipper::Coord_orth rotted = co.transform(it->second);
                        residue_atoms[iat]->x = rotted.x();
                        residue_atoms[iat]->y = rotted.y();
                        residue_atoms[iat]->z = rotted.z();
                     }

                     // 	 std::cout << " standard residue has moved to these positions: "
                     //             << std::endl;
                     // 	 for(int iat=0; iat<nResidueAtoms; iat++) {
                     // 	    std::cout << residue_atoms[iat]->name << " "
                     // 		      << residue_atoms[iat]->x << " "
                     // 		      << residue_atoms[iat]->y << " "
                     // 		      << residue_atoms[iat]->z << "\n";
                     // 	 }
                     // add the atom of std_res to res, deleting excess.
                     mutate_internal(res, std_residue, it->first);
                     istate = 1;
                  }
               }
            }
         }
      }
   g.standard_residues_asc.mol->DeleteSelection(selHnd);
   return istate;
}



// Delete all atoms of residue, and copy over all atoms of std_residue
// except the carbonyl O, which is "correct" in the initial
// coordinates and likely to be wrong in the standard residue.
//
void
molecule_class_info_t::mutate_internal(mmdb::Residue *residue, mmdb::Residue *std_residue,
				       const std::string &alt_conf) {

   atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);
   delete_ghost_selections();

   bool residue_had_hydrogen_atoms = coot::util::residue_has_hydrogens_p(residue);
   coot::residue_spec_t rs(residue);

   coot::util::mutate_internal(residue, std_residue, alt_conf, is_from_shelx_ins_flag,
                               graphics_info_t::default_new_atoms_b_factor);

   atom_sel.mol->FinishStructEdit(); // not sure if this is needed here.

   if (residue_had_hydrogen_atoms)
     add_hydrogen_atoms_to_residue(rs);

   atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
   atom_sel.mol->FinishStructEdit();
   coot::util::pdbcleanup_serial_residue_numbers(atom_sel.mol);

   // regenerate atom selection
   atom_sel = make_asc(atom_sel.mol);
   have_unsaved_changes_flag = 1;
   make_bonds_type_checked(); // calls update_ghosts();

}

void
molecule_class_info_t::mutate_chain(const std::string &chain_id,
				    const coot::chain_mutation_info_container_t &mutation_info,
				    mmdb::PResidue *SelResidues,
				    int nSelResidues,
				    bool renumber_residues_flag) {

   if (mutation_info.insertions.size() > 0 ||
       mutation_info.deletions.size() > 0 ||
       mutation_info.mutations.size() > 0) {
      make_backup(__FUNCTION__);

      // Don't backup each mutation, insertion etc - just do it before
      // and after.
      bool save_backup_state = backup_this_molecule;
      backup_this_molecule = false;

//       std::cout << "mutate chain " << mutation_info.insertions.size()
// 		<< " insertions" << std::endl;
//       std::cout << "mutate chain " << mutation_info.deletions.size()
// 		<< " deletions" << std::endl;
//       std::cout << "mutate chain " << mutation_info.mutations.size()
// 		<< " mutations" << std::endl;

      mutation_info.print();

      if (renumber_residues_flag) {
	 int imod=1;
	 mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
	 int n_chains = model_p->GetNumberOfChains();
	 for (int ich=0; ich<n_chains; ich++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ich);
	    std::string chain_chain_id = chain_p->GetChainID();
	    if (chain_chain_id == chain_id) {
	       simplify_numbering_internal(chain_p);
	    }
	 }
      }


      // do the operations in this order:
      // mutations
      // deletions
      // insertions
      int n_mutations  = 0;
      int n_mutations_failed = 0;
      int n_deletions  = 0;
      int n_insertions = 0;

      // But first let's save the original ordering
      std::vector<int> original_seqnums(nSelResidues);
      for (int i=0; i<nSelResidues; i++)
	 original_seqnums[i] = SelResidues[i]->GetSeqNum();


      // --------------- mutations ----------------
      for (unsigned int i=0; i<mutation_info.mutations.size(); i++) {

	 std::string target_type = mutation_info.mutations[i].second;
	 coot::residue_spec_t rs = mutation_info.mutations[i].first;
	 // std::cout << "DEBUG:: chains: " << chain_id << " " << rs.chain << std::endl;
	 int SelectionHandle = atom_sel.mol->NewSelection();
	 mmdb::PResidue *local_residues;
	 int local_n_selected_residues;
	 atom_sel.mol->Select(SelectionHandle, mmdb::STYPE_RESIDUE, 0,
			      chain_id.c_str(),
			      rs.res_no, rs.ins_code.c_str(),
			      rs.res_no, rs.ins_code.c_str(),
			      "*", "*", "*", "*",
			      mmdb::SKEY_NEW
			      );
	 atom_sel.mol->GetSelIndex(SelectionHandle, local_residues,
				   local_n_selected_residues);
	 // We save the local residues so that we can do a deletion of
	 // the atom selection before making the mutation. We can't do
	 // a DeleteSelection after a mutation (crash) and if we do it
	 // before, then we lose the selected residues!  So make a
	 // vector to keep them in.  Lovely mmdb.
	 std::vector<mmdb::Residue *> residues_vec(local_n_selected_residues);
	 for (int i=0; i<local_n_selected_residues; i++)
	    residues_vec[i] = local_residues[i];
	 atom_sel.mol->DeleteSelection(SelectionHandle);
	 if (local_n_selected_residues > 0) {
	    mutate(residues_vec[0], target_type);
	    n_mutations++;
	 } else {
	    std::cout << "ERROR:: No residue selected for  mutation " << chain_id
		      << " " << rs.res_no << " " << rs.ins_code << " -> "
		      << target_type << std::endl;
	    n_mutations_failed++;
	 }

	 // Nope.... Can't DeleteSelection after mods.
	 // atom_sel.mol->DeleteSelection(SelectionHandle);
      }

      if (n_mutations_failed > 0) {
	 // Why did the mutations go wrong?  Let's see the residues fo atom_sel.mol:
	 int SelectionHandle = atom_sel.mol->NewSelection();
	 mmdb::PResidue *local_residues;
	 int local_n_selected_residues;
	 atom_sel.mol->Select(SelectionHandle, mmdb::STYPE_RESIDUE, 0,
			      chain_id.c_str(),
			      mmdb::ANY_RES, "*",
			      mmdb::ANY_RES, "*",
			      "*", "*", "*", "*",
			      mmdb::SKEY_NEW
			      );
	 atom_sel.mol->GetSelIndex(SelectionHandle, local_residues,
				   local_n_selected_residues);
	 std::cout << "=== Residues currently in molecule === " << std::endl;
	 for (int ires=0; ires<local_n_selected_residues; ires++)
	    std::cout << local_residues[ires]->GetChainID() << " "
		      << local_residues[ires]->GetSeqNum() << " "
		      << local_residues[ires]->GetInsCode() << " "
		      << local_residues[ires]->GetResName() << "\n";
      }


      // --------------- insertions ----------------
      std::cout << "apply resno updates... " << std::endl;
      for (unsigned int i=0; i<mutation_info.insertions.size(); i++) {
	 coot::mutate_insertion_range_info_t r = mutation_info.insertions[i];
	 int offset = r.types.size();
	 n_insertions++;
	 for (int ires=0; ires<nSelResidues; ires++) {
	    if (SelResidues[ires]) {
	       if (SelResidues[ires]->seqNum >= r.start_resno) {
// 		  std::cout << "DEBUG:: for mut_ins[" << i << "] " << r << " incrementing "
// 			    << SelResidues[ires]->seqNum
// 			    << " to " << SelResidues[ires]->seqNum + offset << std::endl;
		  SelResidues[ires]->seqNum += offset;
	       }
	    }
	 }
      }
      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();

      // --------------- deletions ----------------
      std::vector<std::pair<mmdb::Residue *, int> > residues_for_deletion;
      for (unsigned int i=0; i<mutation_info.deletions.size(); i++) {
	 coot::residue_spec_t rs = mutation_info.deletions[i];
	 int SelectionHandle = atom_sel.mol->NewSelection();
	 mmdb::PResidue *local_residues;
	 int local_n_selected_residues;
	 atom_sel.mol->Select(SelectionHandle, mmdb::STYPE_RESIDUE, 0,
			      chain_id.c_str(),
			      rs.res_no, rs.ins_code.c_str(),
			      rs.res_no, rs.ins_code.c_str(),
			      "*", "*", "*", "*",
			      mmdb::SKEY_NEW
			      );
	 atom_sel.mol->GetSelIndex(SelectionHandle, local_residues,
				   local_n_selected_residues);
// 	 std::cout << "DEBUG:: delete section: residue spec :" << chain_id << ": "
// 		   << rs.resno << " :" << rs.insertion_code
// 		   << ": selected " << local_n_selected_residues
// 		   << " residues\n";
	 if (local_n_selected_residues > 0) {
// 	    std::cout << "DEBUG:: marking for deleting "
// 		      << local_residues[0]->GetChainID() << " "
// 		      << local_residues[0]->GetSeqNum()  << " "
// 		      << local_residues[0]->GetResName() << "\n";
	    n_deletions++;
	    residues_for_deletion.push_back(std::pair<mmdb::Residue *, int> (local_residues[0], rs.res_no));
	 }
	 atom_sel.mol->DeleteSelection(SelectionHandle);
      }

      for (unsigned int ird=0; ird<residues_for_deletion.size(); ird++) {
	 // delete residues_for_deletion[ird].first;
	 mmdb::Residue *residue_p = residues_for_deletion[ird].first;
	 if (residue_p) {
	    int seqnum = residue_p->GetSeqNum();
	    mmdb::pstr inscode = residue_p->GetInsCode();
	    residue_p->chain->DeleteResidue(seqnum, inscode);
	    residues_for_deletion[ird].first = NULL;
	 }
	 // Now renumber the followig residues because they should
	 // have a seqNum of 1 less than they had before the deletion.
	 for (int ires=0; ires<nSelResidues; ires++) {
	    if (SelResidues[ires]) {
	       if (original_seqnums[ires] > residues_for_deletion[ird].second) {
		  SelResidues[ires]->seqNum--;
	       }
	    } else {
	       std::cout << "INFO:: found a null residue at " << ires
			 << std::endl;
	    }
	 }
      }

      std::cout << "INFO:: Applied " << n_insertions << " insertions " << std::endl;
      std::cout << "INFO:: Applied " << n_mutations << " mutations " << std::endl;
      std::cout << "INFO:: Applied " << n_deletions << " deletions " << std::endl;

      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      update_molecule_after_additions();
      backup_this_molecule = save_backup_state;
   }
}



coot::chain_mutation_info_container_t
molecule_class_info_t::align_and_mutate(const std::string chain_id,
					const coot::fasta &fasta_seq,
					bool renumber_residues_flag,
					mmdb::realtype wgap, mmdb::realtype wspace) {

   coot::chain_mutation_info_container_t mutation_info;
   std::string target = fasta_seq.sequence;

   mmdb::Manager *mol = atom_sel.mol;
   if (mol) {
      int selHnd = mol->NewSelection();
      mmdb::PResidue *SelResidues = NULL;
      int nSelResidues;

      mol->Select(selHnd, mmdb::STYPE_RESIDUE, 0,
		  chain_id.c_str(),
		  mmdb::ANY_RES, "*",
		  mmdb::ANY_RES, "*",
		  "*",  // residue name
		  "*",  // Residue must contain this atom name?
		  "*",  // Residue must contain this Element?
		  "*",  // altLocs
		  mmdb::SKEY_NEW // selection key
		  );
      mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
      if (nSelResidues > 0) {

	 // I don't know if we can do this here, but I do know we
	 // mol->DeleteSelection(selHnd); // can't DeleteSelection after mods.
	 mutation_info = align_on_chain(chain_id, SelResidues, nSelResidues, target,
					wgap, wspace);
	 mutate_chain(chain_id,
		      mutation_info,
		      SelResidues, nSelResidues,
		      renumber_residues_flag);
      }
   } else {
      std::cout << "ERROR:: null mol in align_and_mutate" << std::endl;
   }
   return mutation_info;
}

#include "utils/align-utils.hh"

coot::chain_mutation_info_container_t
molecule_class_info_t::align_on_chain(const std::string &chain_id,
				      mmdb::PResidue *SelResidues,
				      int nSelResidues,
				      const std::string &target,
				      mmdb::realtype wgap,
				      mmdb::realtype wspace,
				      bool is_nucleic_acid_flag,
				      bool console_output) const {

   coot::chain_mutation_info_container_t ch_info(chain_id);
   bool debug = false;

   if (console_output || debug) {
      std::cout << "\n----- input to align_on_chain() -----------------" << std::endl;
      std::cout << "        chain " << chain_id << std::endl;
      if (0)
	 for (int i=0; i<nSelResidues; i++)
	    std::cout << "        sel-residue: " << coot::residue_spec_t(SelResidues[i])
		      << std::endl;

      std::cout << "        target " <<  target << std::endl;
      std::cout << "        wgap  " << wgap << std::endl;
      std::cout << "        is_nucleic_acid_flag " << is_nucleic_acid_flag << std::endl;
      std::cout << "        console_output " << console_output << std::endl;
      std::cout << "---------------------------------------------------" << std::endl;
   }

   std::vector<std::pair<mmdb::Residue *, int> > vseq =
      coot::util::sort_residues_by_seqno(SelResidues, nSelResidues);

   bool allow_ligands = false;
   std::string model = coot::util::model_sequence(vseq, allow_ligands);
   if (console_output || debug) {
      std::cout << "INFO:: input model  sequence: " << model  << std::endl;
      std::cout << "INFO:: input target sequence: " << target  << std::endl;
   }

   mmdb::math::Alignment align;

   // 20080601.  Don't monkey with these.  If I uncomment the next
   // line (as it used to be) then the I get an N-term alignment error
   // when there is a deletion at the N-term of my model (rnase).
   //
   // align.SetScores(0.5, -0.2);; // 2.0, -1 are the defaults.

   // It seems to me now that it is the gap (and space) penalty that
   // is the important issue.
   //
   // default values (it seems)  now passed parameters
   // mmdb::realtype wgap = 0.0;
   // mmdb::realtype wspace = -1.0;
   // wgap = -3.0;
   // wspace = -0.4;

   std::string stripped_target = coot::util::remove_whitespace(target);

   if (debug) {
      std::cout << "debug:::: Align() with gap penalty: " << wgap
		<< " and extension penalty: " << wspace << std::endl;
      std::cout << "calling Align with 1 " << model << std::endl;
      std::cout << "calling Align with 2 " << stripped_target << std::endl;
   }

   align.SetAffineModel(wgap, wspace);
   align.Align(model.c_str(), stripped_target.c_str());

   ch_info.alignedS = align.GetAlignedS();
   ch_info.alignedT = align.GetAlignedT();

   if (debug) {
      std::cout << "debug:::: Align() on model  " << model << std::endl;
      std::cout << "debug:::: Align() on target " << stripped_target << std::endl;
      std::cout << "debug:::: Align() GetAlignedS:  " << ch_info.alignedS << std::endl;
      std::cout << "debug:::: Align() GetAlignedT:  " << ch_info.alignedT << std::endl;
   }

   ch_info.alignedS_label = name_;
   ch_info.alignedT_label = "target sequence:";
   ch_info.alignment_score = std::pair<bool, float> (true, align.GetScore());

   { // make a string for the GUI - displayed as the result
      std::string as = "<tt>";
      as += ".  " + name_ + "\n";
      as += ". target sequence:\n";
      std::string aligned = align.GetAlignedS();
      std::string target  = align.GetAlignedT();
      std::string matches = coot::alignment_matches(aligned, target);
      as += "   " + aligned + "\n";
      as += "   " + matches + "\n";
      as += "   " + target  + "\n";
      as += "</tt>"; // is this pango at work?
      ch_info.alignment_string = as;
   }

   if (console_output) {
      bool new_style_output = true;
      if (! new_style_output) {
	 // old - alignment no so clear
	 std::cout << ">  ";
	 std::cout << name_ << std::endl;
	 std::cout << align.GetAlignedS() << std::endl;
	 std::cout << "> target seq: \n" << align.GetAlignedT() << std::endl;
	 std::cout << "INFO:: alignment score " << align.GetScore() << std::endl;
      } else {

	 std::string aligned = align.GetAlignedS();
	 std::string target  = align.GetAlignedT();
	 std::string matches = coot::alignment_matches(aligned, target);

	 if (false) { // old style
	    std::cout << ">  " << name_ << "\n";
	    std::cout << "> target sequence:\n";
	    std::cout << "   " << aligned << "\n";
	    std::cout << "   " << matches << "\n";
	    std::cout << "   " << target  << "\n";
	 }

	 // this doesn't make tt format on my mac, maybe it does on PC?

	 std::string s;
	 s += "<tt>Alignment model vs. target\n\n";
	 s += output_alignment_in_blocks(aligned, target, matches);
	 s += "INFO:: alignment score ";
	 s += coot::util::int_to_string(align.GetScore());
	 s += "</tt>\n";

	 // debug
	 // s = "<tt>some text here</tt>\nxyz\n";

	 std::cout << s;
	 ch_info.alignment_string = s;

      }
   }

   // Before the use of SetAffineModel()
   // we got something like:
   // DVSGTVCLSALPPEATDTLNLIASDGPFPYSQD----F--------Q-------NRESVLPTQSYGYYHEYTVITP
   // DTSGTVCLS-LPPEA------IASDGPFPYSQDGTERFDSCVNCAWQRTGVVFQNRESVLPTQSYGYYHEYTVITP
   // (bleugh).
   // 20091028: Now I don't know what alignment we will get.

   // Consider
   // pdb seq: "AVSGTVCLSA"
   // target: "TAVSGTVCLSAT"
   // s ->    "-AVSGTVCLSA-"
   // indexing 001234567899

   // So GetAlignedS() and GetAlignedT() return strings of the same length
   //
   // We run across the length of the returned string, looking for differences:

   // model  has residue             :  Simple mutate.
   // target has different residue   :

   // model  has residue             :  There is an insertion in the model.
   // target has "-"                 :  Delete that residue, offset residues numbers
   //                                :  to the right by -1.
   //
   // model  has "-"                 :  The is a deletion in the model.
   // target has residue             :  Insert a gap for the residue by offsetting
   //                                :  residue numbers to the right by +1.
   //
   // s is from source (our pdb file)
   // t is the target sequence

   std::string s=align.GetAlignedS();
   std::string t=align.GetAlignedT();

   // we need to match the position in SelResidues to the position
   // after alignment (it has had "-"s inserted into it).
   std::vector<int> selindex(s.length());
   int sel_offset = 0;
   for (unsigned int iseq_indx=0; iseq_indx<s.length(); iseq_indx++) {
      int trial_idx = iseq_indx - sel_offset;
      if (trial_idx < nSelResidues)
	 selindex[iseq_indx] = trial_idx;
      else
	 selindex[iseq_indx] = -1; // something bad happened in the setting of selindex values.

      if (debug)
	 std::cout << "assigned: selindex[" << iseq_indx << "]=" << selindex[iseq_indx]
		   << std::endl;

      // for next round
      if (s[iseq_indx] == '-')
	 sel_offset++;
   }

   if (s.length() == t.length()) {

      if (debug) {
	 for (unsigned int iseq_indx=0; iseq_indx<s.length(); iseq_indx++)
	    std::cout << "   s array: " << iseq_indx << " " << s[iseq_indx] << std::endl;

	 if (debug)
	    for (unsigned int i=0; i<selindex.size(); i++)
	       std::cout << "    selindex [" << i << "] is " << selindex[i] << "\n";


	 for (unsigned int i=0; i<s.length(); i++) {
	    if (selindex[i] > -1)
	       std::cout << "   sequence-check: " << i << " " << t[i] << " " << s[i] << " "
			 << coot::util::three_letter_to_one_letter(SelResidues[selindex[i]]->GetResName())
			 << std::endl;
	    else
	       std::cout << "   sequence-check: " << i << " target " << t[i] << " seq " << s[i]
                         << " -1" << std::endl;
	 }
      }

      // std::cout << "DEBUG:: s.length() " << s.length() << std::endl;
      // std::cout << "DEBUG:: nSelResidues " << nSelResidues << std::endl;

      std::string inscode("");
      int ires = 0;
      int res_no_running = 0;

      for (unsigned int iseq_indx=0; iseq_indx<s.length(); iseq_indx++) {

	 int res_idx = selindex[iseq_indx];

	 if (res_idx == -1) {

	    // the model is missing residues at the C-terminus

	    if (s[iseq_indx] == '-') {
	       if (t[iseq_indx] != '-') {

		  coot::residue_spec_t res_spec(res_no_running + 1);
		  res_no_running++; // for next insertion

		  std::string target_type(1, t[iseq_indx]);

		  if (! is_nucleic_acid_flag)
		     target_type = coot::util::single_letter_to_3_letter_code(t[iseq_indx]);

		  ch_info.add_insertion(res_spec, target_type);
		  if (false)
		     std::cout << "DEBUG:: Insert residue  " << res_spec << " " << target_type
			       << " " << iseq_indx << " " << t[iseq_indx]
			       << " with ires " << ires << std::endl;
	       }
	    }

	 } else {

	    // sane res_index

	    if (s[iseq_indx] != '-') {
	       if (debug)
		  std::cout << "DEBUG:: just set res_idx to " << res_idx << " because "
			    << "selindex[" << iseq_indx << "] was " << selindex[iseq_indx]
			    << " of " << nSelResidues << " selected residues" << std::endl;
	       res_no_running = SelResidues[res_idx]->GetSeqNum();
	    }

	    if (false)
	       std::cout << "iseq_indx " << iseq_indx << " "
			 << s[iseq_indx] << " vs " << t[iseq_indx] << std::endl;

	    if (s[iseq_indx] != t[iseq_indx]) {

	       // These only make sense when the aligned residue (in s) was not "-"
	       if (s[iseq_indx] != '-') {
		  ires            = SelResidues[selindex[iseq_indx]]->GetSeqNum();
		  inscode = SelResidues[selindex[iseq_indx]]->GetInsCode();
	       }

	       //	    std::cout << "DEBUG:: ires: " << ires << std::endl;

	       // Case 1: (simple mutate)
	       if ((s[iseq_indx] != '-') && t[iseq_indx] != '-') {
		  // 	       std::cout << "mutate res number " << ires << " "
		  // 			 << s[iseq_indx] << " to " << t[iseq_indx] << std::endl;
		  std::string target_type =
		     coot::util::single_letter_to_3_letter_code(t[iseq_indx]);
		  coot::residue_spec_t res_spec(ires);
		  ch_info.add_mutation(res_spec, target_type);
	       }

	       // Case 2: model had insertion
	       if ((s[iseq_indx] != '-') && t[iseq_indx] == '-') {

		  // 	       for (unsigned int i=iseq_indx+1; i<s.length(); i++)
		  // 		  resno_offsets[i] -= 1;

		  coot::residue_spec_t res_spec(ires);
		  if (0)
		     std::cout << "DEBUG:: Delete residue number " << iseq_indx << " "
			       << s[iseq_indx] << " " << res_spec << std::endl;

		  ch_info.add_deletion(res_spec);
	       }

	       // Case 3: model has a deletion
	       if ((s[iseq_indx] == '-') && t[iseq_indx] != '-') {
		  // 	       for (unsigned int i=iseq_indx+1; i<s.length(); i++)
		  // 		  resno_offsets[i] += 1;
		  // ires will be for the previous residue.  It was not
		  // set for this one.

		  // 20090902
		  coot::residue_spec_t res_spec(res_no_running + 1);
		  res_no_running++; // for next insertion

		  std::string target_type(1, t[iseq_indx]);

		  if (! is_nucleic_acid_flag)
		     target_type = coot::util::single_letter_to_3_letter_code(t[iseq_indx]);

		  ch_info.add_insertion(res_spec, target_type);
		  if (0)
		     std::cout << "DEBUG:: Insert residue  " << res_spec << " " << target_type
			       << " " << iseq_indx << " " << t[iseq_indx]
			       << " with ires " << ires << std::endl;
	       }
	    }
	 }
      }
   }
   ch_info.rationalize_insertions();
   return ch_info;
}



std::string
molecule_class_info_t::output_alignment_in_blocks(const std::string &aligned,
                                                  const std::string &target,
                                                  const std::string &matches) const {

   std::string o;

    // they should all be the same length, don't bother trying to handle the case
    // when they are not.

   std::size_t la = aligned.length();
   std::size_t lt = target.length();
   std::size_t lm = matches.length();

   if (la != lt) return o;
   if (la != lm) return o;

   std::size_t block_size = 80;
   bool do_blocks = true;
   std::string rem_ali = aligned;
   std::string rem_tar = target;
   std::string rem_mat = matches;
   while (do_blocks) {
      std::size_t la = rem_ali.length();
      std::size_t lt = rem_tar.length();
      std::size_t lm = rem_mat.length();
      // std::cout << "\n";
      // std::cout << " aligned: " << rem_ali.substr(0, block_size) << "\n";
      // std::cout << "          " << rem_mat.substr(0, block_size) << "\n";
      // std::cout << "  target: " << rem_tar.substr(0, block_size) << "\n";

      o += " aligned: ";
      o += rem_ali.substr(0, block_size);
      o += "\n";

      o += "          ";
      o += rem_mat.substr(0, block_size);
      o += "\n";

      o += "  target: ";
      o += rem_tar.substr(0, block_size);
      o += "\n\n";

      if (la < block_size) {
	 do_blocks = false;
      } else {
	 rem_ali = rem_ali.substr(80, std::string::npos);
	 rem_tar = rem_tar.substr(80, std::string::npos);
	 rem_mat = rem_mat.substr(80, std::string::npos);
      }
   }
   return o;
}



// Try to align on all chains - pick the best one and return it in the
// second.  If there is no chain that matches within match_frag
// (e.g. 0.95) then return 0 as first and a blank in second.  Also
// return the chain_id.
//
std::pair<bool, std::pair<std::string, coot::chain_mutation_info_container_t> >
molecule_class_info_t::try_align_on_all_chains(const std::string &target, float match_fragment_crit, mmdb::realtype wgap, mmdb::realtype wspace) const {

   coot::chain_mutation_info_container_t cmi;
   bool success = 0;
   float match_frag_best = 2;
   std::string chain_id_best;

   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   mmdb::Chain *chain_p;
   int n_chains = model_p->GetNumberOfChains();
   if (target.length() > 0) {
      for (int ichain=0; ichain<n_chains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres;
	 mmdb::PResidue *residue_table = 0;
	 chain_p->GetResidueTable(residue_table, nres);
	 std::string chain_id = chain_p->GetChainID();

	 // Only try to align if this chain does not have an assigned
	 // sequence already.
	 //
	 bool already_assigned = false;
	 for (unsigned int ii=0; ii<input_sequence.size(); ii++) {
	    if (input_sequence[ii].first == chain_id) {
	       already_assigned = true;
	       break;
	    }
	 }

	 if (! already_assigned) {
	    coot::chain_mutation_info_container_t mic =
	       align_on_chain(chain_id, residue_table, nres, target, wgap, wspace);
	    float sum_changes = mic.single_insertions.size() + mic.deletions.size() +  mic.mutations.size();
	    // was it close? (small number of differences)
	    float match_frag = sum_changes/float(target.length());
	    if (match_frag < (1 - match_fragment_crit)) {
	       if (match_frag < match_frag_best) {
		  match_frag_best = match_frag;
		  cmi = mic;
		  success = 1;
		  chain_id_best = chain_id;
	       }
	    }
	 }
      }
   }
   std::pair<std::string, coot::chain_mutation_info_container_t> p(chain_id_best, cmi);
   return std::pair<bool, std::pair<std::string, coot::chain_mutation_info_container_t> > (success, p);
}



// redundant now that we have coot-util functions.
//
std::string
molecule_class_info_t::make_model_string_for_alignment(mmdb::PResidue *SelResidues,
						       int nSelResidues) const {

   std::vector<std::pair<mmdb::Residue *, int> > vseq =
      coot::util::sort_residues_by_seqno(SelResidues, nSelResidues);
   return coot::util::model_sequence(vseq);
}


std::pair<bool, std::vector<coot::chain_mutation_info_container_t> >
molecule_class_info_t::residue_mismatches(mmdb::realtype alignment_wgap, mmdb::realtype alignment_wspace) const {

   std::vector<coot::chain_mutation_info_container_t> ar;
   bool status = 0;
   std::vector<coot::residue_spec_t> v;
   if (input_sequence.size() > 0) {

      int imod = 1;
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      mmdb::Chain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 std::string chain_id = chain_p->GetChainID();
	 for (unsigned int ich=0; ich<input_sequence.size(); ich++) {
	    if (input_sequence[ich].first == chain_id) {
	       status = 1;
	       mmdb::PPResidue SelResidues = 0;
	       int n_residues = 0;
	       chain_p->GetResidueTable(SelResidues, n_residues);
	       coot::chain_mutation_info_container_t ali =
		  align_on_chain(chain_id, SelResidues, n_residues,
				 input_sequence[ich].second,
				 alignment_wgap,
				 alignment_wspace);
	       ali.print();
	       ar.push_back(ali);
	    }
	 }
      }
   }
   return std::pair<bool, std::vector<coot::chain_mutation_info_container_t> > (status, ar);
}



std::pair<bool, std::string>
molecule_class_info_t::find_terminal_residue_type(const std::string &chain_id, int resno,
						  mmdb::realtype alignment_wgap,
						  mmdb::realtype alignment_wspace,
						  bool is_nucleic_acid_flag) const {


   bool found = false;
   std::string type = "None";
   std::string target = "";

   for (unsigned int iseq=0; iseq<input_sequence.size(); iseq++) {
      if (input_sequence[iseq].first == chain_id) {
	 target = input_sequence[iseq].second;
	 break;
      }
   }

   if (target != "") {

      mmdb::Manager *mol = atom_sel.mol;
      if (mol) {
	 int selHnd = mol->NewSelection(); // d
	 mmdb::PResidue *SelResidues = NULL;
	 int nSelResidues;

	 mol->Select(selHnd, mmdb::STYPE_RESIDUE, 0,
		     chain_id.c_str(),
		     mmdb::ANY_RES, "*",
		     mmdb::ANY_RES, "*",
		     "*",  // residue name
		     "*",  // Residue must contain this atom name?
		     "*",  // Residue must contain this Element?
		     "*",  // altLocs
		     mmdb::SKEY_NEW // selection key
		     );
	 mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
	 if (nSelResidues > 0) {

	    coot::chain_mutation_info_container_t mi =
	       align_on_chain(chain_id, SelResidues, nSelResidues, target,
			      alignment_wgap, alignment_wspace, is_nucleic_acid_flag);
	    // mi.print();

	    coot::residue_spec_t search_spec(chain_id, resno);
	    try {
	       type = mi.get_residue_type(search_spec);
	       found = true;
	    }
	    catch (const std::runtime_error &mess) {
	       std::cout << "WARNING:: on catch() failed to find " << search_spec
			 << " for an insertion " << mess.what() << std::endl;
	    }
	 }
	 mol->DeleteSelection(selHnd);
      }
   }

   return std::pair<bool, std::string> (found, type);
}

int
molecule_class_info_t::mutate_by_overlap(const std::string &chain_id, int res_no, const std::string &new_type) {

   int status = 0;
   graphics_info_t g;
   mmdb::Residue *residue_p = get_residue(coot::residue_spec_t(chain_id, res_no, ""));

   if (residue_p) {
      std::string current_residue_type = residue_p->GetResName();
      g.Geom_p()->check_and_try_dynamic_add(current_residue_type, imol_no, g.cif_dictionary_read_number);
      g.cif_dictionary_read_number++;
      std::pair<bool, coot::dictionary_residue_restraints_t> rp_current =
         g.Geom_p()->get_monomer_restraints(current_residue_type, imol_no);
      if (rp_current.first) {
         mmdb::Manager *mol = g.molecules[imol_no].atom_sel.mol;
         const auto &restraints_current_type = rp_current.second;

         g.Geom_p()->check_and_try_dynamic_add(new_type, imol_no, g.cif_dictionary_read_number);
         g.cif_dictionary_read_number++;
         std::pair<bool, coot::dictionary_residue_restraints_t> rp_new_type =
            g.Geom_p()->get_monomer_restraints(new_type, imol_no);
         if (rp_new_type.first) {
            const auto &restraints_new_type = rp_new_type.second;

            mmdb::Residue *restraints_new_type_residue_p = restraints_new_type.GetResidue(false, 10.0f);

            if (restraints_new_type_residue_p) {

               status = coot::util::mutate_by_overlap(residue_p, mol, restraints_current_type, restraints_new_type);
               atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
               atom_sel.mol->FinishStructEdit();
               atom_sel.regen_atom_selection();
               make_bonds_type_checked();

               if (status == 0)
                  logger.log(log_t::WARNING, "mutate_by_overlap() failed");
            } else {
               std::string m = "mutate_by_overlap() restraints_new_type_residue_p was null";
               logger.log(log_t::WARNING, m);
            }
         } else {
            std::string m = "failed to get restraints for " + new_type;
            logger.log(log_t::WARNING, logging::function_name_t(__FUNCTION__), m);
         }
      }
   }
   return status;
}



// Here is something that does DNA/RNA
int
molecule_class_info_t::mutate_base(const coot::residue_spec_t &res_spec, std::string type,
				   bool use_old_style_naming) {

   int istat=0;

   bool this_residue_is_DNA = false;
   std::string this_residue_res_name = get_residue_name(res_spec);
   if (this_residue_res_name == "DA" || this_residue_res_name == "DG" ||
       this_residue_res_name == "DC" || this_residue_res_name == "DT" ||
       this_residue_res_name == "D3" || this_residue_res_name == "DI" ||
       this_residue_res_name == "DN" || this_residue_res_name == "DU")
      this_residue_is_DNA = true;

   // refmac_nuc_type is the type of the residue that we extract from
   // the standard residues molecule.
   //
   std::string refmac_nuc_type = type;
   // we match the requested residue type to the residue type that is
   // in the standard residues file.

   {
      // modern names input,

      if (this_residue_is_DNA) {
         if (type == "A") refmac_nuc_type = "DA";
         if (type == "G") refmac_nuc_type = "DG";
         if (type == "T") refmac_nuc_type = "DT";
         if (type == "U") refmac_nuc_type = "DU";
         if (type == "C") refmac_nuc_type = "DC";
      } else {
         if (type == "A") refmac_nuc_type = "A";
         if (type == "G") refmac_nuc_type = "G";
         if (type == "T") refmac_nuc_type = "T";
         if (type == "U") refmac_nuc_type = "U";
         if (type == "C") refmac_nuc_type = "C";
      }
      if (type == "DA")	 refmac_nuc_type = "DA";
      if (type == "DG")	 refmac_nuc_type = "DG";
      if (type == "DT")	 refmac_nuc_type = "DT";
      if (type == "DC")	 refmac_nuc_type = "DC";
   }

   if (atom_sel.n_selected_atoms > 0) {

      int n_models = atom_sel.mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) {

	 mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
	 mmdb::Chain *chain_p;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 if (nchains <= 0) {
	    std::cout << "bad nchains in molecule " << nchains
		      << std::endl;
	 } else {
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain_p = model_p->GetChain(ichain);
	       if (chain_p != NULL) {
		  std::string mol_chain_id = chain_p->GetChainID();
		  if (mol_chain_id == res_spec.chain_id) {
		     int nres = chain_p->GetNumberOfResidues();
		     mmdb::PResidue residue_p;
		     for (int ires=0; ires<nres; ires++) {
			residue_p = chain_p->GetResidue(ires);
			if (residue_p->GetSeqNum() == res_spec.res_no) {
			   if (res_spec.ins_code == residue_p->GetInsCode()) {

			      // Found the residue (nucleotide in this case):

			      mmdb::Residue *std_base = get_standard_residue_instance(refmac_nuc_type);
			      if (std_base) {
				 mutate_base_internal(residue_p, std_base, use_old_style_naming);
				 istat = 1;
			      } else {
				 std::cout << "WARNING:: Oops - can't find standard residue for type "
					   << type << std::endl;
			      }
			   }
			}
			if (istat)
			   break;
		     }
		  }
	       }
	       if (istat)
		  break;
	    }
	 }
      }
   }
   if (istat) {
      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      make_bonds_type_checked();
   }
   // std::cout << "--------------------------------- done mci::mutate_base() --------------------------" << std::endl;
   return istat;
}



// Here std_base is at some arbitary position when passed.
//
void
molecule_class_info_t::mutate_base_internal(mmdb::Residue *residue, mmdb::Residue *std_base,
					    bool use_old_names) {

   make_backup(__FUNCTION__);

   if (0)
      std::cout << "DEBUG:: mutate_base_internal():: residue name: "
		<< residue->GetResName()
		<< " using std_base " << std_base->GetResName()
		<< " use_old_names flag " << use_old_names<< std::endl;
   coot::util::mutate_base(residue, std_base, use_old_names);
   have_unsaved_changes_flag = 1;
}

int
molecule_class_info_t::exchange_chain_ids_for_seg_ids() {

   short int changed_flag = 0;
   if (atom_sel.n_selected_atoms > 0) {

      // Make a list/vector of Chains so that we can delete them after
      // we have created the new chains.
      int n_models = atom_sel.mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) {

	 mmdb::Model *model_p = atom_sel.mol->GetModel(imod);

	 std::vector<int> chain_vec;
	 int nchains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<nchains; ichain++) {
	    chain_vec.push_back(ichain);
	 }

	 //
	 mmdb::Atom *at = 0;
	 std::vector<std::pair<std::vector<mmdb::Atom *>, std::string> > atom_chain_vec;
	 std::vector<mmdb::Atom *> running_atom_vec;
	 std::string current_seg_id = "unassigned-camel-burger";
	 std::string current_chain_id = "unassigned-camel-burger";
	 for (int iat=0; iat<atom_sel.n_selected_atoms; iat++) {
	    at = atom_sel.atom_selection[iat];
	    std::string this_atom_segid = at->segID;
	    std::string this_atom_chain_id = at->GetChainID();
	    if (current_seg_id != this_atom_segid) {
	       // add the running chain to the list then.  We construct a chain-id
	       // for it.
	       if (running_atom_vec.size() > 0) {
		  std::string constructed_seg_id = current_chain_id;  // (that of the previous atom)
		  constructed_seg_id += ":";
		  constructed_seg_id += current_seg_id; // also previous atom
		  atom_chain_vec.push_back(std::pair<std::vector<mmdb::Atom *>, std::string>(running_atom_vec, constructed_seg_id));
	       }
	       // start a new chain
	       running_atom_vec.clear();
	       running_atom_vec.push_back(at);
	       // set for next push
	       current_seg_id = this_atom_segid;
	       current_chain_id = this_atom_chain_id;
	    } else {
	       // all in all it's just a-nother atom in the chain...
	       running_atom_vec.push_back(at);
	    }
	 }
	 if (running_atom_vec.size() > 0) {
	    atom_chain_vec.push_back(std::pair<std::vector<mmdb::Atom *>, std::string>(running_atom_vec, current_seg_id));
	 }

	 // OK, so we have vector of vectors of atoms.  We need to make
	 // new atoms and residues to put them in.
	 std::cout << "INFO:: Creating " << atom_chain_vec.size() << " new chains\n";
	 for (unsigned int inch=0; inch<atom_chain_vec.size(); inch++) {
	    mmdb::Chain *chain_p = new mmdb::Chain;
	    const char *chid = atom_chain_vec[inch].second.c_str();
	    chain_p->SetChainID(chid);
	    // Add the chain to the model:
	    model_p->AddChain(chain_p);
	    mmdb::Residue *residue_p = 0;
	    int prev_resno = -999999;
	    const char *prev_ins = "";
	    for (unsigned int iat=0; iat<atom_chain_vec[inch].first.size(); iat++) {
	       at = atom_chain_vec[inch].first[iat];
	       char *resname = at->GetResName();
	       char *ins     = at->GetInsCode();
	       int seqnum      = at->GetSeqNum();
	       if ((seqnum != prev_resno) || strcmp(ins,prev_ins)) {
		  // we need to make a new residue attached to chain_p
// 		  std::cout << "debug:: Making new residue " << resname << " "
// 			    << seqnum << " :" << ins << ": prev_resno "
// 			    << prev_resno << " prev_ins :" << prev_ins << ":"
// 			    << std::endl;
		  residue_p = new mmdb::Residue(chain_p, resname, seqnum, ins);
		  prev_resno = seqnum;
		  prev_ins = ins;
	       }
	       mmdb::Atom *atom_p = new mmdb::Atom(residue_p); // does an AddAtom
	       atom_p->Copy(at); // doesn't touch atom's res
	    }
	 }
	 // now (finally) delete the chains of the model:
	 for (unsigned int ich=0; ich<chain_vec.size(); ich++) {
	    model_p->DeleteChain(chain_vec[ich]);
	 }
      }
   }

   // Let's just imagine that it always happens...
   have_unsaved_changes_flag = 1;
   atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
   atom_sel.mol->FinishStructEdit();
   atom_sel = make_asc(atom_sel.mol);
   make_bonds_type_checked();
   return changed_flag;
}

//
void
molecule_class_info_t::spin_search(clipper::Xmap<float> &xmap,
				   const std::string &chain_id,
				   int resno,
				   const std::string &ins_code,
				   const std::pair<std::string, std::string> &direction_atoms,
				   const std::vector<std::string> &moving_atoms_list) {

   mmdb::Residue *res = get_residue(chain_id, resno, ins_code);

   if (res) {
      // the first atom spec is not used in spin_search, so just bodge one in.
      coot::atom_spec_t atom_spec_tor_1(chain_id, resno, ins_code, direction_atoms.first, "");
      coot::atom_spec_t atom_spec_tor_2(chain_id, resno, ins_code, direction_atoms.first, "");
      coot::atom_spec_t atom_spec_tor_3(chain_id, resno, ins_code, direction_atoms.second, "");
      coot::atom_spec_t atom_spec_tor_4(chain_id, resno, ins_code, moving_atoms_list[0], "");

      coot::torsion tors(0, // not used
			 atom_spec_tor_1, atom_spec_tor_2,
			 atom_spec_tor_3, atom_spec_tor_4);

      // What is dir?
      mmdb::Atom *dir_atom_1 = 0;
      mmdb::Atom *dir_atom_2 = 0;
      mmdb::PPAtom residue_atoms;
      int nResidueAtoms;
      res->GetAtomTable(residue_atoms, nResidueAtoms);
      for (int iat=0; iat<nResidueAtoms; iat++) {
	 if (direction_atoms.first == residue_atoms[iat]->name)
	    dir_atom_1 = residue_atoms[iat];
	 if (direction_atoms.second == residue_atoms[iat]->name)
	    dir_atom_2 = residue_atoms[iat];
      }

      if (dir_atom_1 && dir_atom_2) {

	 float angle = coot::util::spin_search(xmap, res, tors).first;

	 if (angle < -1000) { // an error occured
	    std::cout << "ERROR:: something bad in spin_search" << std::endl;
	 } else {

	    make_backup(__FUNCTION__);

	    // Now rotate the atoms in moving_atoms_list about dir;
	    //
	    clipper::Coord_orth orig(dir_atom_2->x, dir_atom_2->y, dir_atom_2->z);
	    clipper::Coord_orth  dir(dir_atom_2->x - dir_atom_1->x,
				     dir_atom_2->y - dir_atom_1->y,
				     dir_atom_2->z - dir_atom_1->z);

	    for (unsigned int i_mov_atom = 0; i_mov_atom<moving_atoms_list.size(); i_mov_atom++) {


	       mmdb::PPAtom residue_atoms;
	       int nResidueAtoms;
	       res->GetAtomTable(residue_atoms, nResidueAtoms);
	       for (int iat=0; iat<nResidueAtoms; iat++) {
		  if (moving_atoms_list[i_mov_atom] == residue_atoms[iat]->name) {

		     clipper::Coord_orth pt(residue_atoms[iat]->x,
					    residue_atoms[iat]->y,
					    residue_atoms[iat]->z);

		     clipper::Coord_orth co = coot::util::rotate_around_vector(dir, pt, orig, angle);
		     residue_atoms[iat]->x = co.x();
		     residue_atoms[iat]->y = co.y();
		     residue_atoms[iat]->z = co.z();
		  }
	       }
	    }
	    //
	    have_unsaved_changes_flag = 1;
	    make_bonds_type_checked(); // calls update_ghosts();
	 }

      } else {
	 std::cout << "direction atoms not found" << std::endl;
      }

   } else {
      std::cout << "residue not found in coordinates molecule" << std::endl;
   }
}

std::vector<std::pair<coot::residue_spec_t, float> >
molecule_class_info_t::em_ringer(const clipper::Xmap<float> &xmap) const {

   std::vector<std::pair<coot::residue_spec_t, float> > vr;

   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<nres; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	    std::string residue_name = residue_p->GetResName();
	    // PRO has a CG but we don't want to ringer it. Other residues with CG we do.
	    // and spin_atom doesn't know about residue types or atom names, it's generic.
	    if (residue_name != "PRO") {
	       coot::residue_spec_t spec(residue_p);
	       std::vector<std::string> v;
	       v.push_back(" CG ");
	       density_results_container_t drc = spin_atom(xmap, spec, " N  ", " CA ", " CB ", v);
	       if (drc.scored_points.size() == 2) {
		  float a1 = drc.scored_points[0].angle;
		  float a2 = drc.scored_points[1].angle;
		  float delta = a2 - a1;
		  if (delta < -180.0) delta += 360.0;
		  if (delta >  180.0) delta -= 360.0;
		  // std::cout << "spin_atom " << spec << " " << a1 << " " << a2 << " " << delta << "\n";
		  std::pair<coot::residue_spec_t, float> p(spec, delta);
		  vr.push_back(p);
	       }
	    }
         }
      }
   }

   return vr;
}

// #include "density-results-container-t.hh" now in the molecule_class_info_t constructor.

// c.f. spin search above, here we need 3 atom names for the spin - so that we
// can use the torsion of the moving atom.
// tors_1 ref             e.g. N
// tors_2 base            e.g. CA
// tors_3 tip             e.g. CB
// tors_4 spnning atom    e.g. CG
//
density_results_container_t
molecule_class_info_t::spin_atom(const clipper::Xmap<float> &xmap,
                                 const coot::residue_spec_t &spec,
				 const std::string &direction_atoms_ref,
				 const std::string &direction_atoms_base,
				 const std::string &direction_atoms_tip,
				 const std::vector<std::string> &moving_atoms_list) const {

   // this is the wrong container. I (at the moment) want just current torsion and best torsion
   //
   density_results_container_t drc;

   mmdb::Residue *residue_p = get_residue(spec);
   if (residue_p) {
      if (moving_atoms_list.size() > 0) {
         // it would be nice if we could make atom specs from residue specs
         const std::string &chain_id = spec.chain_id;
         const std::string &ins_code = spec.ins_code;
         int res_no = spec.res_no;
         std::string alt_conf = "";
         coot::atom_spec_t atom_spec_tor_1(chain_id, res_no, ins_code, direction_atoms_ref,  alt_conf);
         coot::atom_spec_t atom_spec_tor_2(chain_id, res_no, ins_code, direction_atoms_base, alt_conf);
         coot::atom_spec_t atom_spec_tor_3(chain_id, res_no, ins_code, direction_atoms_tip,  alt_conf);
         coot::atom_spec_t atom_spec_tor_4(chain_id, res_no, ins_code, moving_atoms_list[0], alt_conf);

         coot::torsion tors(0, // imol, not used
			    atom_spec_tor_1, atom_spec_tor_2,
			    atom_spec_tor_3, atom_spec_tor_4);

         std::vector<mmdb::Atom *> ma = tors.matching_atoms(residue_p);
	 if (ma.size() == 4) {
	    float best_tors_angle = coot::util::spin_search(xmap, residue_p, tors).second; // degrees, relative to N
	    coot::atom_quad q(ma[0], ma[1], ma[2], ma[3]);
	    float tors_atoms = q.torsion(); // degrees
	    double delta = best_tors_angle - tors_atoms;
	    float dv = 0;
	    clipper::Coord_orth pos(0,0,0);
	    density_results_t sp_1(pos, tors_atoms, dv); // position angle density
	    density_results_t sp_2(pos, best_tors_angle, dv);
	    drc.scored_points.push_back(sp_1);
	    drc.scored_points.push_back(sp_2);
	 }
      }
   }
   return drc;
}



// maybe this can go down to coot coord utils?  Needs mutatation stuff though.
//
// Maybe this should be part of the the imol_coords molecule_class_info_t?
// so we can do mutations?
//
// Return a list residues that should be deleted from the original molecule.
//
// Manipulate mol (adding atoms).
//
// OK, poly_ala_mol is a fragment of structure that looks like a part
// of this mol (atom_sel.mol).  The residues of poly_ala_mol that
// overlaps atom_sel.mol have atoms specs in mmdb_residues - so that
// we can delete them afterwards (that is, if the can be mutated
// (i.e. they have a sensible letter (not ?))) in the best_seq.
//
// I need to fiddle with poly_ala_mol then merge it into atom_sel.mol
// and delete the residues *from atom_sel.mol* that got mutated in
// poly_ala_mol.
//
// We pass imol_map because we use auto_fit_best_rotamer which uses
// imol_map (eeugk!)
//
int
molecule_class_info_t::apply_sequence(int imol_map, mmdb::Manager *poly_ala_mol,
				      std::vector<coot::residue_spec_t> mmdb_residues,
				      std::string best_seq, std::string chain_id,
				      int resno_offset,
				      const coot::protein_geometry &pg) {

   std::cout << "--------------------------------- apply_sequence() --------------------------" << std::endl;

   int istat = 0;
   short int have_changes = 0;
   std::vector<coot::residue_spec_t> r_del;

//    std::cout << "DEBUG:: residue vector len " <<  mmdb_residues.size() << std::endl;
//    std::cout << "DEBUG:: best sequence  len " <<  best_seq.length() << std::endl;
   make_backup(__FUNCTION__);

   int selHnd = poly_ala_mol->NewSelection();
   poly_ala_mol->Select(selHnd, mmdb::STYPE_RESIDUE, 1,
			"*",
			mmdb::ANY_RES, "*",
			mmdb::ANY_RES, "*",
			"*",  // residue name
			"*",  // Residue must contain this atom name?
			"*",  // Residue must contain this Element?
			"*",  // altLocs
			mmdb::SKEY_NEW // selection key
			);
   mmdb::PResidue *SelResidues = 0;
   int nSelResidues;
   poly_ala_mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
   if (nSelResidues != int(best_seq.length())) {
      std::cout << "oops residue mismatch " << best_seq.length() << " " << nSelResidues
		<< std::endl;
   } else {

      for (unsigned int ichar=0; ichar<best_seq.length(); ichar++) {

	 char seq_char = best_seq[ichar];
	 if (seq_char == '?') {
	    std::cout << "bypassing ? at " << ichar << std::endl;

	    // but we still need to set the sequence offset
	    mmdb::Residue *poly_ala_res = SelResidues[ichar];

	    poly_ala_res->seqNum = resno_offset + ichar;
	    if (ichar < mmdb_residues.size())
		    r_del.push_back(mmdb_residues[ichar]);
	 } else {
	    std::string res_type = coot::util:: single_letter_to_3_letter_code(best_seq[ichar]);
	    if (res_type != "") {
	       have_changes = 1;
	       std::cout << "Mutating to type " << res_type << " at " << ichar << std::endl;
	       mmdb::Residue *std_res = get_standard_residue_instance(res_type);
	       if (std_res) {
		  // Get res in poly_ala_mol

		  mmdb::Residue *poly_ala_res = SelResidues[ichar];

		  std::cout << "Mutating poly_ala residue number " << poly_ala_res->GetSeqNum()
			    << std::endl;

		  std::string alt_conf = ""; // presumably.
		  coot::util::mutate(poly_ala_res, std_res, alt_conf, 0); // not shelx
		  poly_ala_res->GetChain()->SetChainID(chain_id.c_str());
		  poly_ala_res->seqNum = resno_offset + ichar;
		  if (ichar < mmdb_residues.size())
		     r_del.push_back(mmdb_residues[ichar]);
	       }
	    }
	 }
      }
   }
   poly_ala_mol->DeleteSelection(selHnd);

   if (have_changes) {
      poly_ala_mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      poly_ala_mol->FinishStructEdit();

      // now fiddle with atom_sel.mol, deleting r_del residues and
      // adding poly_ala_mol.

      // don't backup each residue deletion:
      int bs = backups_state();
      turn_off_backup();
      for (unsigned int ird=0; ird<r_del.size(); ird++) {
// 	 std::cout << "deleting :" << r_del[ird].chain << ": " << r_del[ird].resno
// 		   << std::endl;

	 int model_number_ANY = mmdb::MinInt4;
	 delete_residue(model_number_ANY,
			r_del[ird].chain_id,
			r_del[ird].res_no,
			r_del[ird].ins_code);
      }

      insert_coords_internal(make_asc(poly_ala_mol));

      // Now auto fit.  Get the residue list from poly_ala_mol, but
      // apply auto fits to atom_sel.mol:
      int imod = 1;
      mmdb::Model *poly_ala_model_p = poly_ala_mol->GetModel(imod);
      mmdb::Chain *poly_ala_chain_p;
      int nchains = poly_ala_model_p->GetNumberOfChains();
      graphics_info_t g;
      for (int ichain=0; ichain<nchains; ichain++) {
	 poly_ala_chain_p = poly_ala_model_p->GetChain(ichain);
	 int nres = poly_ala_chain_p->GetNumberOfResidues();
	 mmdb::Residue *poly_ala_residue_p = 0;
	 for (int ires=0; ires<nres; ires++) {
	    istat = 1;
	    poly_ala_residue_p = poly_ala_chain_p->GetResidue(ires);
	    int rotamer_mode = ROTAMERSEARCHLOWRES;
	    auto_fit_best_rotamer(rotamer_mode,
				  poly_ala_residue_p->GetSeqNum(), "",
 				  poly_ala_residue_p->GetInsCode(),
 				  poly_ala_residue_p->GetChainID(),
 				  imol_map,
 				  g.rotamer_fit_clash_flag,
 				  g.rotamer_lowest_probability,
				  *g.Geom_p());
	 }
      }
      if (bs)
	 turn_on_backup();
   }


   if (have_changes) {
      atom_sel = make_asc(atom_sel.mol);
      make_bonds_type_checked();
   }

   have_unsaved_changes_flag = 1;

   return istat;
}


// return success on residue type match
// success: 1, failure: 0.
int
molecule_class_info_t::mutate_single_multipart(int ires_serial,
					       const std::string &chain_id,
					       const std::string &target_res_type) {

   int istat = 0;
   if (atom_sel.n_selected_atoms > 0) {
      mmdb::Chain *chain_p;
      int nres;
      int nchains = atom_sel.mol->GetNumberOfChains(1) ;
      for (int ichain =0; ichain<nchains; ichain++) {
	 chain_p = atom_sel.mol->GetChain(1,ichain);
	 if (std::string(chain_id) == std::string(chain_p->GetChainID())) {
	    nres = chain_p->GetNumberOfResidues();
	    if (ires_serial < nres) {
	       mmdb::Residue *res_p = chain_p->GetResidue(ires_serial);
	       if (res_p) {

		  if (std::string(res_p->name) == target_res_type) {

                     if (false) // debug
                        std::cout << "residue type match for ires = " << ires_serial << std::endl;
		     istat = 1; // success

		  } else {

		     // OK, do the mutation:

		     // get an instance of a standard residue of type target_res_type
		     mmdb::Residue *std_res = get_standard_residue_instance(target_res_type); // a deep copy
		     // move the standard res to position of res_p
		     // move_std_residue(moving_residue, (const) reference_residue);
		     if (std_res) {
			istat = move_std_residue(std_res, res_p);

			if (istat) {
			   std::vector<std::string> alt_confs = coot::util::get_residue_alt_confs(res_p);
			   for (unsigned int i=0; i<alt_confs.size(); i++)
			      mutate_internal(res_p, std_res, alt_confs[i]);
			} else {
			   std::cout << "WARNING:  Not mutating residue due to missing atoms!\n";
			}
			// atom_selection and bonds regenerated in mutate_internal
		     } else {
			std::cout << "ERROR failed to get residue of type :" << target_res_type
				  << ":" << std::endl;
		     }
		  }
	       } else {
		  std::cout << "ERROR:: in mutate_single_multipart oops - can't get residue"
			    << " with ires_serial: " << ires_serial << std::endl;
	       }
	    } else {
	       std::cout << "PROGRAMMER ERROR: out of range residue indexing" << std::endl;
	    }
	 }
      }
   }
   return 0 + istat;
}

// trash all other residues in imol_ligand:
int
molecule_class_info_t::delete_all_except_res(mmdb::Residue *res) {

   int state = 0;
   if (atom_sel.n_selected_atoms > 0) {
      make_backup(__FUNCTION__);
//       std::cout << "DEBUG:: molecule number " << imol_no << " contains "
// 		<< atom_sel.mol->GetNumberOfModels() << " models"
// 		<< std::endl;
      for (int imod=1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
	 mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
	 mmdb::Chain *chain_p;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<nchains; ichain++) {
	    chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    mmdb::PResidue residue_p;
	    for (int ires=0; ires<nres; ires++) {
	       residue_p = chain_p->GetResidue(ires);
	       if (residue_p != res) {
		  // std::cout << "Deleting residue " << residue_p << std::endl;
		  chain_p->DeleteResidue(ires);
		  residue_p = NULL;
		  state = 1;
	       }
	    }
	 }
      }
      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      have_unsaved_changes_flag = 1;
      atom_sel = make_asc(atom_sel.mol);
      make_bonds_type_checked();
   }
   return state;
}

bool
molecule_class_info_t::residue_has_TER_atom(mmdb::Residue *res_p) const {

   int n_residue_atoms;
   mmdb::PPAtom residue_atoms;
   bool has_ter = 0;
   if (res_p) {
      res_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++) {
	 if (residue_atoms[i]->isTer()) {
	    has_ter = 1;
	    break;
	 }
      }
   }
   return has_ter;
}


void
molecule_class_info_t::remove_TER_on_last_residue(mmdb::Chain *chain_p) {

   int n_residues = chain_p->GetNumberOfResidues();
   if (n_residues > 0) {
      mmdb::Residue *r = chain_p->GetResidue(n_residues-1); // last residue
      if (r)
	 remove_TER_internal(r);
   }
}

void
molecule_class_info_t::remove_TER_on_residue_if_last_residue(mmdb::Chain *chain_p, mmdb::Residue *residue_p) {

   int n_residues = chain_p->GetNumberOfResidues();
   if (n_residues > 0) {
      mmdb::Residue *r = chain_p->GetResidue(n_residues-1); // last residue
      if (r == residue_p)
	 remove_TER_internal(r);
   }
}

// remove TER record from residue
//
void
molecule_class_info_t::remove_TER_internal(mmdb::Residue *res_p) {

   int n_residue_atoms;
   mmdb::PPAtom residue_atoms;
   bool deleted = 0;
   if (res_p) {
      res_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++) {
	 if (residue_atoms[i]->isTer()) {
	    res_p->DeleteAtom(i);
	    deleted = 1;
	 }
      }
   }
   if (deleted) {
      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
   }
}

void
molecule_class_info_t::remove_ter_atoms(const coot::residue_spec_t &spec) {  // from all models

   // do a backup only if this chain has a residue with a TER atom.
   bool has_ter = false;
   mmdb::Residue *residue_p = get_residue(spec);
   if (residue_p)
      has_ter = residue_has_TER_atom(residue_p);
   if (has_ter) {
      make_backup(__FUNCTION__);
      remove_TER_internal(residue_p);
   }
}



// Not a mutate function, (merely) a function that sees how close the
// sequence match is to each of the chains.
std::vector<coot::chain_mutation_info_container_t>
molecule_class_info_t::sequence_comparison_to_chains(const std::string &target_sequence) const {

   std::vector<coot::chain_mutation_info_container_t> mutation_info_vec;
   if (atom_sel.mol) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(1);
      if (model_p) {
	 int n_chains = model_p->GetNumberOfChains();
	 for (int ich=0; ich<n_chains; ich++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ich);
	    std::string chain_id = chain_p->GetChainID();
	    mmdb::realtype wgap   =  0.0;    //defaults
	    mmdb::realtype wspace = -1.0;
	    int nSelResidues;
	    mmdb::PPResidue SelResidues = 0;
	    chain_p->GetResidueTable(SelResidues, nSelResidues);
	    bool console_output = false;
	    coot::chain_mutation_info_container_t mutation_info =
	       align_on_chain(chain_id, SelResidues, nSelResidues, target_sequence, wgap, wspace,
			      console_output);
	    mutation_info_vec.push_back(mutation_info);
	 }
      }
   }
   return mutation_info_vec;
}


// mutate and autofit the residues
//
int
molecule_class_info_t::nudge_residue_sequence(const std::string &chain_id,
					      int res_no_range_start,
					      int res_no_range_end,
					      int nudge_by,
					      short int nudge_residue_numbers_also) {

   int status = 0;

   // we want to call this multiple times:
   // mutate_single_multipart(ires_serial, chain_id, const std::string &target_res_type)


   if (res_no_range_start < res_no_range_end) {
      // given 20, 22: we want a range of 2 and to mutate 20,21,22 (3 residues)
      int range = res_no_range_end - res_no_range_start;


      // first get the sequence of the residues as they currently are.
      bool seq_status = true;
      std::vector<std::string> current_types;
      for (int i=0; i<=range; i++) {
	 mmdb::Residue *r = get_residue(chain_id, res_no_range_start+i, "");
	 if (! r) {
	    status = false;
	    break;
	 } else {
	    current_types.push_back(r->GetResName());
	 }
      }

      if (seq_status && int(current_types.size()) == (range+1)) {

	 make_backup(__FUNCTION__);

	 status = 1; // we did something

	 for (int i_offset=0; i_offset<=range; i_offset++) {
	    int i_serial_no = residue_serial_number(chain_id, res_no_range_start+i_offset, "");
	    if (i_serial_no != -1) {
	       int new_type_idx = i_offset - nudge_by;
	       if (new_type_idx >= 0 && new_type_idx < int(current_types.size())) {
		  std::string new_res_type = current_types[new_type_idx];
		  mutate_single_multipart(i_serial_no, chain_id, new_res_type);
	       }
	    }
	 }

	 // if nudge_residue_numbers_also (as is often the case,
	 // because the residue number is correct, but the residues
	 // are out of position by one along the chain):

	 if (nudge_residue_numbers_also) {

	    std::vector<mmdb::Residue *> delete_these; // due to overlapping numbers
	    for (int i_offset=0; i_offset<=range; i_offset++) {

	    }

	    for (int i_offset=0; i_offset<=range; i_offset++) {
	       mmdb::Residue *r = get_residue(chain_id, res_no_range_start+i_offset, "");
	       if (r) {
		  r->seqNum -= nudge_by;
	       }
	    }
	 }

	 atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
	 atom_sel.mol->FinishStructEdit();
	 have_unsaved_changes_flag = 1;
	 make_bonds_type_checked();

      } else {
	 std::cout << "WARNING:: Null residue in nudge range " << std::endl;
      }
   } else {
      std::cout << "WARNING:: bad sequence numbering" << std::endl;
   }
   return status;
}
