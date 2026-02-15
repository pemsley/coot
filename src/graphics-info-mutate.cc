/* src/graphics-info-mutate.cc
 *
 * Copyright 2004, 2005 by The University of York
 * Copyright 2008, 2009 by The University of Oxford.
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
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 02110-1335, USA
 */

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"

// For stat, mkdir:
#include <sys/types.h>
#include <sys/stat.h>

#if !defined _MSC_VER
#include <unistd.h>
#else
#define PKGDATADIR "C:/coot/share"
#endif

#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif


#include <iostream>

#include <mmdb2/mmdb_manager.h>
#include "coords/mmdb-extras.hh"
#include "coords/mmdb.hh"
#include "coords/mmdb-crystal.hh"
#include "coords/Cartesian.hh"
#include "coords/Bond_lines.hh"

#include <gtk/gtk.h>  // must come after mmdb_manager on MacOS X Darwin
// #include <GL/glut.h>  // for some reason...  // Eh?

#include "graphics-info.h"
#include "interface.h"

#include "coot-utils/coot-coord-utils.hh"
#include "molecule-class-info.h"

#include "manipulation-modes.hh"

#include "utils/logging.hh"
extern logging logger;


// #include "coot-utils.hh"


// This should be in coot-coord-utils.
//
atom_selection_container_t
graphics_info_t::add_side_chain_to_terminal_res(atom_selection_container_t asc,
						const std::string &res_type,
						const std::string &terminus_type,
						bool add_other_residue_flag) {

   // the add_other_residue_flag is now passed to this function so that, when building
   // forwards, we can add to the penultimate residue rather than the last residue
   // (get_last_residue_in_chain()).

   if (false)
      std::cout << "here we are in add_side_chain_to_terminal_res " << res_type
		<< " " << terminus_type << " with add_other_residue_flag " << add_other_residue_flag
		<< std::endl;

   atom_selection_container_t rasc = asc;
   int istat;
   molecule_class_info_t molci;
   short int display_in_display_control_widget_status = 0;
   molci.install_model(0, asc, Geom_p(), "terminal residue", display_in_display_control_widget_status);

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
      mmdb::Residue *std_res = molci.get_standard_residue_instance(target_res_type);

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
		  int nres = chain->GetNumberOfResidues();
		  mmdb::Residue *residue_p = 0;

		  if (terminus_type == "N" || terminus_type == "MN") {
		     if (add_other_residue_flag)
			residue_p = coot::util::get_second_residue_in_chain(chain);
		     else
			residue_p = coot::util::get_first_residue_in_chain(chain);
		  }

		  if (terminus_type == "C" || terminus_type == "MC") {
		     if (add_other_residue_flag)
			residue_p = coot::util::get_penultimate_residue_in_chain(chain);
		     else
			residue_p = coot::util::get_last_residue_in_chain(chain);
		  }

                  if (terminus_type == "singleton")
                     residue_p = coot::util::get_last_residue_in_chain(chain);

                  if (false)
                     std::cout << "------- add_side_chain_to_terminal_res() here with residue_p "
                               << coot::residue_spec_t(residue_p) << std::endl;

		  if (residue_p) {

		     //

		     int istat = molci.move_std_residue(std_res_copy, residue_p);

		     if (istat) {

			mmdb::PPAtom residue_atoms;
			int nResidueAtoms;
			residue_p->GetAtomTable(residue_atoms, nResidueAtoms);

			mmdb::PPAtom std_residue_atoms;
			int n_std_ResidueAtoms;
			std_res_copy->GetAtomTable(std_residue_atoms, n_std_ResidueAtoms);

			// set the b factor for the new atoms.
			for(int i=0; i<n_std_ResidueAtoms; i++) {
			   std_residue_atoms[i]->tempFactor = default_new_atoms_b_factor;
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

	       molci.atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
	       molci.atom_sel.mol->FinishStructEdit();

	       rasc = make_asc(molci.atom_sel.mol);

// 	       for (int ii=0; ii<rasc.n_selected_atoms; ii++)
// 		  std::cout << "New res with cb: " << rasc.atom_selection[ii] << "\n";
	    }
	 }
      }
   }
   return rasc;
}

// We need the action of Add Alt Conf... to provide a new alt conf in
// the the most likely (rotamer 0) position.
//
// The problem is that the moving residues from an "Add Alt Conf" are
// not necessarily a whole residue.  We may need to dig out the main
// chain atoms too.
//
// Also, we have a moving atoms, so we don't need to generate a new
// one as we do above in generate_moving_atoms_from_rotamer
//
// Or do we want to act as if "Rotamers" was pressed at the end of an
// "Add Alt Conf" - yes, I think that we do... How does that work then?
// We simply need the an atom index of an atom in the new alt conf.
//
// atom_spec_to_atom_index?
// So I need to know the spec of a newly created atom
//

void
graphics_info_t::do_mutation(int imol, const coot::residue_spec_t &res_spec,
                              const std::string &residue_type, short int do_stub_flag) {

   // this is called from both simple mutate (in_mutate_define)
   // and auto_fit_mutate (in_mutate_auto_fit_define).
   // mutate_auto_fit_residue_imol will not be set in the simple mutate case

   if (false)
      std::cout << "::::::::::::::::::::::: in do_mutation() with imol " << imol
                << " and residue_type_chooser_auto_fit_flag "
                << residue_type_chooser_auto_fit_flag << std::endl;

   if (residue_type_chooser_auto_fit_flag) {

      if (! is_valid_model_molecule(imol)) {
         std::cout << "ERROR:: invalid model molecule number in do_mutation() "
                   << mutate_auto_fit_residue_imol << std::endl;
         return;
      }

      std::cout << "do_mutation() here with imol " << imol << std::endl;

      int local_atom_index = -1;

      // set an sane local_atom_index:
      mmdb::Residue *residue_p = molecules[imol].get_residue(res_spec);
      if (residue_p) {
         for (int iat=0; iat<molecules[imol].atom_sel.n_selected_atoms; iat++) {
            mmdb::Atom *at = molecules[imol].atom_sel.atom_selection[iat];
            std::string at_chain_id = at->GetChainID();
            if (res_spec.chain_id == at_chain_id) {
               int res_no = at->GetSeqNum();
               if (res_spec.res_no == res_no) {
                  std::string ins_code = at->GetInsCode();
                  if (res_spec.ins_code == ins_code) {
                     local_atom_index = iat;
                  }
               }
            }
            if (local_atom_index != -1)
               break;
         }
      }

      if (local_atom_index == -1) {
         std::cout << "WARNING:: bad local atom index " << local_atom_index << std::endl;
         return;
      }

      molecules[mutate_auto_fit_residue_imol].mutate(local_atom_index, residue_type, do_stub_flag);

      // 20071005 No longer check for stub state.  It doesn't make
      // sense to autofit a stub.  So ignore the stub state and just
      // autofit as normal.

      int imol_map = Imol_Refinement_Map();
      // std::cout << "DEBUG:: do_mutation(): here with imol_map = " << imol_map << std::endl;
      if (imol_map >= 0) {

         // float f =
         int mode = graphics_info_t::rotamer_search_mode;
         // this signature needs to be changed to a residue spec.
         float result = molecules[imol].auto_fit_best_rotamer(mode, local_atom_index,
                                                              imol_map, rotamer_fit_clash_flag,
                                                              rotamer_lowest_probability, *Geom_p());
         // std::cout << "DEBUG:: auto_fit_best_rotamer() returned " << result << std::endl;
         logger.log(log_t::INFO, logging::function_name_t("do_mutation"),
                    "auto_fit_best_rotamer returned", result);

         if (mutate_auto_fit_do_post_refine_flag) {
            // Run refine zone with autoaccept, autorange on
            // the "clicked" atom:
            short int auto_range = 1;
            refine(mutate_auto_fit_residue_imol, auto_range,
                   mutate_auto_fit_residue_atom_index,
                   mutate_auto_fit_residue_atom_index);
         }

         // This is the wrong function, isn't it?
         update_go_to_atom_window_on_changed_mol(imol);

         update_validation(imol);

         run_post_manipulation_hook(imol, MUTATED);

      } else {

         // imol map chooser
         show_select_map_frame();
      }

   } else {

      std::cout << "DEBUG:: do_mutation() here with mutate_residue_imol "
                << mutate_residue_imol << std::endl;
      if (is_valid_model_molecule(mutate_residue_imol)) {
         // simple mutation
         // molecules[mutate_residue_imol].mutate(mutate_residue_atom_index, residue_type, do_stub_flag);
         auto &m = molecules[mutate_residue_imol];
         mmdb::Residue *residue_p = m.get_residue(res_spec);
         if (residue_p) {
            m.mutate(residue_p, residue_type, true);
            update_go_to_atom_window_on_changed_mol(mutate_residue_imol);
            update_validation(mutate_auto_fit_residue_imol);
            run_post_manipulation_hook(mutate_residue_imol, MUTATED);
         }
      } else {
         std::cout << "WARNING:: mutate_residue_imol is " << mutate_residue_imol << std::endl;
      }
   }
   graphics_draw();
}

void
graphics_info_t::do_mutation_auto_fit(int imol, const coot::residue_spec_t &res_spec,
                  const std::string &residue_type,
				      short int do_stub_flag) {

   molecules[mutate_residue_imol].mutate(mutate_residue_atom_index, residue_type,
					 do_stub_flag);
   graphics_draw();
   run_post_manipulation_hook(mutate_residue_imol, MUTATED);
}


void
graphics_info_t::read_standard_residues() {

   std::string standard_env_dir = "COOT_STANDARD_RESIDUES";

   const char *filename = getenv(standard_env_dir.c_str());
   if (! filename) {

      std::string standard_file_name = coot::package_data_dir();
      standard_file_name += "/";
      standard_file_name += "standard-residues.pdb";

      struct stat buf;
      int status = stat(standard_file_name.c_str(), &buf);
      if (status != 0) { // standard-residues file was not found in
	                 // default location either...
	 std::cout << "WARNING: Can't find standard residues file in the "
		   << "default location \n";
	 std::cout << "         and environment variable for standard residues ";
	 std::cout << standard_env_dir << "\n";
	 std::cout << "         is not set.";
	 std::cout << " Mutations will not be possible\n";
	 // mark as not read then:
	 standard_residues_asc.read_success = 0;
	 standard_residues_asc.n_selected_atoms = 0;
	 // std::cout << "DEBUG:: standard_residues_asc marked as empty" << std::endl;
      } else {
	 // stat success:
	 standard_residues_asc = get_atom_selection(standard_file_name, false, true, false);
      }
   } else { 
      standard_residues_asc = get_atom_selection(filename, false, true, false);
   }
}



void
graphics_info_t::mutate_chain(int imol, const std::string &chain_id,
			      const std::string &seq,
			      bool do_auto_fit_flag,
			      bool renumber_residues_flag) {

   if (imol < n_molecules()) {
      if (imol >= 0) {
	 if (molecules[imol].has_model()) {
	    std::cout << "INFO:: aligning to mol number " << imol << " chain: "
		      << chain_id << std::endl;
	    coot::chain_mutation_info_container_t mutation_info =
	       molecules[imol].align_and_mutate(chain_id, coot::fasta(seq), renumber_residues_flag,
						alignment_wgap, alignment_wspace);

	    info_dialog_alignment(mutation_info); // protected

	    if (do_auto_fit_flag) {
	       int imol_map = Imol_Refinement_Map();
	       if (is_valid_map_molecule(imol_map)) {

		  // For now let's refine the whole chain
		  //
                  std::vector<std::string> s;
                  s.push_back("fit-chain");
                  s.push_back(coot::util::int_to_string(imol));
                  s.push_back(coot::util::single_quote(chain_id));

                  std::cout << ":::::::::::::::: here 1 with command string s " << schemize_command_strings(s) << std::endl;
#ifdef USE_GUILE
                  std::cout << ":::::::::::::::: here 2 with command string s " << schemize_command_strings(s) << std::endl;
                  safe_scheme_command(schemize_command_strings(s));
#endif // USE_GUILE


		  // This doesn't work (yet) because after alignment,
		  // the residue numbers of the mutations do not
		  // correspond the structure (because of insertions
		  // and deletions).
		  //
		  // =============================================
		  if (0) {
		     if (mutation_info.mutations.size()) {
			int replacement_state = refinement_immediate_replacement_flag;
			refinement_immediate_replacement_flag = 1;
			for (unsigned int i=0; i<mutation_info.mutations.size(); i++) {
			   coot::residue_spec_t res_spec = mutation_info.mutations[i].first;
			   // now autofit rotamer on res_spec.

			   int mode = graphics_info_t::rotamer_search_mode;
			   int clash_flag = 1;
			   float lowest_probability = 0.01;
			   std::string altloc = "";
			   molecules[imol].auto_fit_best_rotamer(mode,
								 res_spec.res_no,
								 altloc,
								 res_spec.ins_code,
								 res_spec.chain_id,
								 imol_map,
								 clash_flag,
								 lowest_probability,
								 *Geom_p());
			   short int auto_range = 1;
			   refine_residue_range(imol,
						res_spec.chain_id, res_spec.chain_id,
						res_spec.res_no, res_spec.ins_code,
						res_spec.res_no, res_spec.ins_code,
						altloc, 0);

			}
			refinement_immediate_replacement_flag = replacement_state;
		     }
		  } // end non-working block. =====================

	       } else {
		  std::cout << "WARNING:: refinement map set " << std::endl;
	       }
	    }
	 }
      }
   }
}



// --------------------------------------
// not really mutation, but hey ho...
// ------- cis trans conversion ---------
void
graphics_info_t::cis_trans_conversion(mmdb::Atom *at, int imol, short int is_N_flag) {

   if (molecules[imol].has_model()) {
      int istatus = molecules[imol].cis_trans_conversion(at, is_N_flag, standard_residues_asc.mol);
      if (istatus > 0)
	 graphics_draw();
   }
}
