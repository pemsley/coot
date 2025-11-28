/* src/graphics-info.cc
 *
 * Copyright 2002, 2003, 2004, 2005, 2006 by The University of York
 * Copyright 2007, 2008, 2009 by the University of Oxford
 * Copyright 2013, 2015, 2016 by Medical Research Council
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
#include "python-3-interface.hh"
#endif

#include "compat/coot-sysdep.h"

#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif

// was used for debugging - not needed now.
// #include <iomanip> // extern char *dcgettext (const char *__domainname, if placed later
#include <algorithm>

#include <iostream>
#include <stdexcept>

#include <gtk/gtk.h>  // must come after mmdb_manager on MacOS X Darwin
// #include <GL/glut.h>  // for some reason...  // Eh?

#include <gdk/gdkkeysyms.h> // for keyboarding (in this case nudge_active_residue) added 20091101

#include <sys/types.h> // for stating
#include <sys/stat.h>

#if !defined _MSC_VER
#include <unistd.h>
#endif

#include "guile-fixups.h"


#include <clipper/core/map_utils.h> // Map_stats
#include <mmdb2/mmdb_manager.h>

#include "coords/mmdb-extras.hh"
#include "coords/mmdb.hh"
#include "coords/mmdb-crystal.hh"
#include "coords/Cartesian.hh"
#include "coords/Bond_lines.hh"

#ifdef USE_DUNBRACK_ROTAMERS
#include "ligand/dunbrack.hh"
#else
#include "ligand/richardson-rotamer.hh"
#endif

#include "skeleton/graphical_skel.h"

#include "graphics-info.h"
#include "interface.h"

#include "molecule-class-info.h"
#include "coot-utils/coot-coord-extras.hh"


#include "globjects.h"
#include "ligand/torsion-general.hh"
#include "ligand/ligand.hh"
#include "ligand/ideal-rna.hh"
#include "ligand/residue_by_phi_psi.hh"

#include "rotate-translate-modes.hh"
#include "ideal/torsion-bonds.hh"

// Including python needs to come after graphics-info.h, because
// something in Python.h (2.4 - chihiro) is redefining FF1 (in
// ssm_superpose.h) to be 0x00004000 (Grrr).
//
// 20100813: Python.h needs to come before to stop"_POSIX_C_SOURCE" redefined problems
//
// #ifdef USE_PYTHON
// #include "Python.h"
// #endif // USE_PYTHON


#include "utils/coot-utils.hh"

#include "widget-from-builder.hh"


void
graphics_info_t::get_restraints_lock(const std::string &calling_function_name) {

   bool unlocked = false;

   // 20191127-PE not this formulation:
   // while (! graphics_info_t::restraints_lock.compare_exchange_weak(unlocked, true) && !unlocked) {

   while (! restraints_lock.compare_exchange_weak(unlocked, true)) {
      std::cout << "WARNING:: calling function: " << calling_function_name
		<< " restraints locked by " << restraints_locking_function_name
		<< std::endl;
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
      unlocked = false;
   }
   // std::cout << "debug:: Got the lock for " << calling_function_name << std::endl;
   restraints_locking_function_name = calling_function_name;
}

void
graphics_info_t::release_restraints_lock(const std::string &calling_function_name) {

   // std::cout << "debug:: release the restraints lock: " << calling_function_name << std::endl;
   restraints_lock = false;
   restraints_locking_function_name = "";

}

// similar for moving atoms:
void
graphics_info_t::get_moving_atoms_lock(const std::string &calling_function_name) {

   bool unlocked = false;
   while (! moving_atoms_lock.compare_exchange_weak(unlocked, true)) {
      std::cout << "WARNING:: calling function: " << calling_function_name
		<< " moving atoms locked by " << moving_atoms_locking_function_name
		<< std::endl;
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
      unlocked = false;
   }

   moving_atoms_locking_function_name = calling_function_name;
}

void
graphics_info_t::release_moving_atoms_lock(const std::string &calling_function_name) {

   moving_atoms_lock = false;
   moving_atoms_locking_function_name = "";
}


void
graphics_info_t::stop_refinement_internal() {

   // c.f. GDK_Escape key press:

   if (continue_threaded_refinement_loop) {
      continue_threaded_refinement_loop = false;
      threaded_refinement_needs_to_clear_up = true;
      clear_hud_buttons(); // if a refinement was running and Esc was pressed, here is where we catch it.
   }
   // now wait until refinement stops... (before we call clear_up_moving_atoms (from the
   // function that calls this function))
   // std::cout << "debug:: stop_refinement_internal() waiting for refinement to stop" << std::endl;
   get_restraints_lock(__FUNCTION__);
   release_restraints_lock(__FUNCTION__);
   // std::cout << "debug:: stop_refinement_internal() refinement stopped" << std::endl;

}



// Idealize the geometry without considering the map.
//
coot::refinement_results_t
graphics_info_t::copy_mol_and_regularize(int imol,
					 int resno_1,
					 std::string inscode_1,
					 int resno_2,
					 std::string inscode_2,
					 std::string altconf,// use this (e.g. "A") or "".
					 std::string chain_id_1) {

   return copy_mol_and_refine(imol, -1, resno_1, inscode_1, resno_2, inscode_2, altconf, chain_id_1);
}




// Regularize *and* fit to density.
//
// Cut and pasted from above.  You might ask why I didn't factor out
// the common stuff.
//
//
// Note that this refinement routine uses moving_atoms_asc.
//
//
coot::refinement_results_t
graphics_info_t::copy_mol_and_refine(int imol_for_atoms,
				     int imol_for_map,
				     int resno_1,
				     std::string inscode_1,
				     int resno_2,
				     std::string inscode_2,
				     std::string altconf,// use this (e.g. "A") or "".
				     std::string chain_id_1) {

   // This now wraps refine_residues_vec

   // can this function be deleted now?

   if (false)
      std::cout << "DEBUG:: In copy_mol_and_refine() refine range: "
		<< "chain  :" << chain_id_1 << ": "
		<< resno_1 << " :" << inscode_1 << ": "
		<< resno_2 << " :" << inscode_2 << ": "
		<< " altconf \"" << altconf << "\" "
		<< "coords mol: " << imol_for_atoms << " map mol: " << imol_for_map
		<< std::endl;

   short int irest = 0; // make 1 when restraints found.

   coot::refinement_results_t rr(0, GSL_CONTINUE, "");

   int imol = imol_for_atoms;
   imol_moving_atoms = imol_for_atoms;  // for use when we accept the
			      // regularization and want to copy the
			      // coordinates back.

   // make the selection and build a new molecule inside restraints.

   short int have_flanking_residue_at_start = 0;
   short int have_flanking_residue_at_end = 0;
   short int have_disulfide_residues = 0;  // other residues are included in the
                                           // residues_mol for disulfide restraints.

   // 9 Sept 2003: The atom selection goes mad if residue with seqnum
   // iend_res+1 does not exist, but is not at the end of the chain.

   // Therefore we will set 2 flags, which tell us if istart_res-1 and
   // iend_res+1 exist.  And we do that by trying to select atoms from
   // them - if they exist, the number of selected atoms will be more
   // than 0.

//    istart_minus_flag = 0;  // from simple restraint code
//    iend_plus_flag    = 0;

   mmdb::Manager *mol = molecules[imol].atom_sel.mol; // short-hand usage

   // We want to check for flanking atoms if the dictionary "group"
   // entry is not non-polymer.  So let's do a quick residue selection
   // of the first residue and find its residue type, look it up and
   // get the group.  If it is "non-polymer", then we can tinker with
   // the have_flanking_residue_at_* flags.

   int SelHnd_first = mol->NewSelection(); // d
   int n_residue_first;
   mmdb::PResidue *residue_first = NULL;
   mol->Select(SelHnd_first, mmdb::STYPE_RESIDUE, 0,
	       chain_id_1.c_str(),
	       resno_1, inscode_1.c_str(),
	       resno_1, inscode_1.c_str(),
	       "*",  // residue name
	       "*",  // Residue must contain this atom name?
	       "*",  // Residue must contain this Element?
	       "*",  // altLocs
	       mmdb::SKEY_NEW); // selection key
   mol->GetSelIndex(SelHnd_first, residue_first, n_residue_first);
   std::string group = "L-peptide";
   if (n_residue_first > 0) {
      std::string residue_type_first = residue_first[0]->name;
      // does a dynamic add if needed.

//       not used:
//       int status =
// 	 geom_p->have_dictionary_for_residue_type(residue_type_first,
// 						  cif_dictionary_read_number);

      std::pair<short int, coot::dictionary_residue_restraints_t> p =
	 geom_p->get_monomer_restraints(residue_type_first, imol_for_atoms);
      if (p.first) {
	 group = p.second.residue_info.group;
      }
      cif_dictionary_read_number++;
   }
   mol->DeleteSelection(SelHnd_first);

   if (group != "non-polymer") { // i.e. it is (or can be) a polymer
      int SelHnd_ends = mol->NewSelection();
      int n_atoms_ends;
      mmdb::PPAtom atoms_end = 0;
      mol->SelectAtoms(SelHnd_ends, 0, chain_id_1.c_str(),
		       resno_1-1, "*", resno_1-1, "*","*","*","*","*");
      mol->GetSelIndex(SelHnd_ends, atoms_end, n_atoms_ends);
      if (n_atoms_ends > 0)
	 have_flanking_residue_at_start = 1; // we have residue istart_res-1
      mol->DeleteSelection(SelHnd_ends);

      SelHnd_ends = mol->NewSelection();
      mol->SelectAtoms(SelHnd_ends, 0, chain_id_1.c_str(),
		       resno_2+1, "*", resno_2+1, "*","*","*","*","*");
      mol->GetSelIndex(SelHnd_ends, atoms_end, n_atoms_ends);
      if (n_atoms_ends > 0)
	 have_flanking_residue_at_end = 1; // we have residue iend_res+1
      mol->DeleteSelection(SelHnd_ends);
   }

   int selHnd = mol->NewSelection();  // d
   int nSelResidues;
   mmdb::PResidue *SelResidues = NULL;
   mol->Select(selHnd, mmdb::STYPE_RESIDUE, 0,
	       chain_id_1.c_str(),
	       resno_1, inscode_1.c_str(),
	       resno_2, inscode_2.c_str(),
	       "*",  // residue name
	       "*",  // Residue must contain this atom name?
	       "*",  // Residue must contain this Element?
	       "*",  // altLocs
	       mmdb::SKEY_NEW // selection key
	       );
   molecules[imol].atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

   // 20100201 - Happy Path

   bool check_hydrogens_too_flag = false;
   // convert to mmdb::Residues vector
   std::vector<mmdb::Residue *> residues;
   for (int ires=0; ires<nSelResidues; ires++) {
      residues.push_back(SelResidues[ires]);
   }

   mol->DeleteSelection(selHnd);

   std::pair<bool, std::vector<std::pair<mmdb::Residue *, std::vector<std::string> > > >
      icheck_atoms = Geom_p()->atoms_match_dictionary(imol_for_atoms, residues, check_hydrogens_too_flag, false);

   if (! icheck_atoms.first) {
      std::cout << "WARNING:: Fail atom check" << std::endl;
      info_dialog_refinement_non_matching_atoms(icheck_atoms.second);
      return rr; // fail
   } else {

      // 	 return copy_mol_and_refine_inner(imol_for_atoms,
      // 					  resno_1, resno_2,
      // 					  nSelResidues, SelResidues,
      // 					  chain_id_1, altconf,
      // 					  have_flanking_residue_at_start,
      // 					  have_flanking_residue_at_end,
      // 					  imol_for_map);

      if (imol_for_map == -1)
	 rr = regularize_residues_vec(imol_for_atoms, residues, altconf, mol);
      else
	 rr = refine_residues_vec(imol_for_atoms, residues, altconf, mol);
   }

   return rr;
}

void
graphics_info_t::show_missing_refinement_residues_dialog(const std::vector<std::string> &res_names,
							 bool run_get_monomer_post_fetch_flag) {

   if (use_graphics_interface_flag) {

      GtkWidget *dialog = widget_from_builder("download_monomers_dialog");
      gtk_widget_set_visible(dialog, TRUE);
      GtkWidget *vbox = widget_from_builder("download_monomers_dialog_vbox_inner");
      gtk_widget_set_visible(vbox, TRUE);
      clear_out_container(vbox);
      for (const auto &rn : res_names) {
         GtkWidget *label = gtk_label_new(rn.c_str());
         char* ccx = new char[rn.size() + 1];
         for (unsigned int ii=0; ii<=rn.size(); ii++) ccx[ii] = 0;
         std::copy(rn.begin(), rn.end(), ccx);
         g_object_set_data(G_OBJECT(label), "comp_id", ccx); // read in on_download_monomers_ok_button_clicked()
         gtk_box_append(GTK_BOX(vbox), label);
      }
      g_object_set_data(G_OBJECT(dialog), "run_get_monomer_post_fetch_flag",
                        GINT_TO_POINTER(run_get_monomer_post_fetch_flag));
   }
}


// static
void
graphics_info_t::info_dialog_missing_refinement_residues(const std::vector<std::string> &res_names) {

   std::set<std::string> res_set;
   for (unsigned int icheck_res=0; icheck_res<res_names.size(); icheck_res++)
      res_set.insert(res_names[icheck_res]);

   std::string problem_residues = "WARNING: Refinement setup failure.\nFailed to find restraints for ";
   problem_residues += "these ";
   problem_residues += std::to_string(res_set.size());
   problem_residues += " residue types:\n";

   std::set<std::string>::const_iterator it;
   unsigned int count = 0;
   for (it=res_set.begin(); it!=res_set.end(); ++it) {
      problem_residues += " ";
      problem_residues += *it;
      count++;
      if (count == 10) {
         problem_residues += "\n";
         count = 0;
      }
   }

   info_dialog(problem_residues);
}


int
graphics_info_t::copy_active_atom_molecule() {

   int imol = -1;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > aa = active_atom_spec();
   if (aa.first) {
      imol = copy_model_molecule(aa.second.first);
   }
   return imol;
}

int
graphics_info_t::copy_model_molecule(int imol) {
   int iret = -1;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      int new_mol_number = graphics_info_t::create_molecule();
      mmdb::Manager *m = graphics_info_t::molecules[imol].atom_sel.mol;
      mmdb::Manager *n = new mmdb::Manager;
      n->Copy(m, mmdb::MMDBFCM_All);
      atom_selection_container_t asc = make_asc(n);
      std::string label = "Copy_of_";
      label += graphics_info_t::molecules[imol].name_;
      const std::vector<drawn_ghost_molecule_display_t> &ghosts = g.molecules[imol].NCS_ghosts();
      std::vector<coot::ghost_molecule_display_t> empty_ghosts; // type mismatch - needed fixing?
      bool shelx_flag = g.molecules[imol].is_from_shelx_ins();
      g.molecules[new_mol_number].install_model_with_ghosts(new_mol_number, asc, g.Geom_p(), label, 1, empty_ghosts,
                                                            shelx_flag, false, false);
      update_go_to_atom_window_on_new_mol();
      iret = new_mol_number;
   }
   return iret;
}


void
graphics_info_t::save_accept_reject_dialog_window_position(GtkWidget *acc_rej_dialog) {

   // 20070801 crash reported by "Gajiwala, Ketan"

   // OK, we can reproduce a problem
   // Refine something
   // Close the window using WM delete window
   // Press return in Graphics window (globjects:key_press_event() GDK_Return case)
   //
   // So, we need to set graphics_info_t::accept_reject_dialog to NULL
   // when we get a WM delete event on the Accept/Reject box

   if (acc_rej_dialog) {
      gint upositionx, upositiony;
      // if (acc_rej_dialog->window) {
      if (true) { // no access to window
	    std::cout << "GTK-FIXME no root origin B" << std::endl;
	 // gdk_window_get_root_origin (acc_rej_dialog->window, &upositionx, &upositiony);
	 // graphics_info_t::accept_reject_dialog_x_position = upositionx;
	 //	 graphics_info_t::accept_reject_dialog_y_position = upositiony;
      } else {
	 std::cout << "ERROR:: Trapped an error in save_accept_reject_dialog_window_position\n"
		   << "        Report to Central Control!\n"
		   << "        (What did you do to make this happen?)\n";
      }
   }
}

void
graphics_info_t::clear_up_glsl_buffers_for_moving_atoms() {

}

void
graphics_info_t::clear_up_moving_atoms_wrapper() {

   hide_atom_pull_toolbar_buttons();

   rebond_molecule_corresponding_to_moving_atoms(); // in clear_up_moving_atoms_wrapper()

   // poke a value into the threaded refinement loop, to stop
   if (continue_threaded_refinement_loop) {
      // and tell it to clear up the moving atoms
      threaded_refinement_needs_to_clear_up = true;
      std::cout << ".... Esc key tells refinement to clean up" << std::endl;
      continue_threaded_refinement_loop = false;

      clear_hud_buttons(); // if a refinement was running and Esc was pressed, here is where we catch it.

   } else {

      // refinement was not running. we can clear up the atoms ourselves
      clear_up_moving_atoms();
      clear_up_glsl_buffers_for_moving_atoms();
      clear_moving_atoms_object();

      clear_hud_buttons();

      draw_bad_nbc_atom_pair_markers_flag = false;

      // remove this at some stage

      if (accept_reject_dialog) {
         if (accept_reject_dialog_docked_flag == coot::DIALOG) {
            save_accept_reject_dialog_window_position(accept_reject_dialog);
            // this calls clear_up_moving_atoms() and clears atom pull restraint.
            // gtk_widget_destroy(accept_reject_dialog);
            gtk_widget_set_visible(accept_reject_dialog, FALSE);
            accept_reject_dialog = 0;
         } else {
            gtk_widget_set_sensitive(graphics_info_t::accept_reject_dialog, FALSE);
         }
      }
   }
}

std::atomic<unsigned int> graphics_info_t::moving_atoms_bonds_lock(0);
std::atomic<bool> graphics_info_t::restraints_lock(false);
std::atomic<bool> graphics_info_t::moving_atoms_lock(false); // not locked
std::string graphics_info_t::restraints_locking_function_name = "unset";
std::string graphics_info_t::moving_atoms_locking_function_name = "unset";

int  graphics_info_t::threaded_refinement_loop_counter = 0;
int  graphics_info_t::threaded_refinement_loop_counter_bonds_gen = -1; // initial value is "less than" so that
                                                                       // the regeneration is activated.
bool graphics_info_t::threaded_refinement_needs_to_clear_up = false; // for Esc usage
bool graphics_info_t::threaded_refinement_needs_to_accept_moving_atoms = false; // for Return usage
bool graphics_info_t::continue_threaded_refinement_loop = false; // also for Esc usage
int  graphics_info_t::threaded_refinement_redraw_timeout_fn_id = -1;
bool graphics_info_t::refinement_of_last_restraints_needs_reset_flag = false;

// put this in graphics-info-intermediate-atoms?
//
// static
void
graphics_info_t::refinement_loop_threaded() {

   // ---- Don't touch the graphics or the gui! ---------

   if (graphics_info_t::restraints_lock) {
      return;
   }

   if (!graphics_info_t::last_restraints) {
      // nothing to refine - this should not happen
      return;
   }

   get_restraints_lock(__FUNCTION__);

   graphics_info_t::threaded_refinement_needs_to_clear_up = false; // set on Esc press
                                                                   // from the main loop
   graphics_info_t::threaded_refinement_needs_to_accept_moving_atoms = false;

   graphics_info_t g;

   coot::restraint_usage_Flags flags = g.set_refinement_flags(); // flags should not be needed for minimize()

   // continue_threaded_refinement_loop = true; not here - set it in the calling function
   while (continue_threaded_refinement_loop) {

      g.update_restraints_with_atom_pull_restraints();

      bool pr_chi_sqds = false; // print inital chi squareds
      int spf = dragged_refinement_steps_per_frame;

      if (refinement_of_last_restraints_needs_reset_flag) {
         g.last_restraints->set_needs_reset();
         refinement_of_last_restraints_needs_reset_flag = false;
      }

      // coot::refinement_results_t rr = g.last_restraints->minimize(flags, spf, pr_chi_sqds);
      coot::refinement_results_t rr = g.last_restraints->minimize(imol_moving_atoms, flags,
								  spf, pr_chi_sqds, *Geom_p());

      // 20231202-PE previous_round_nbc_baddies_atom_index_map is set her so that it
      // can be used in gone_contacts_from_nbc_baddies() in update_bad_nbc_atom_pair_marker_positions()
      //
      // previous_round_nbc_baddies_atom_index_map = rr.nbc_baddies_atom_index_map;
      saved_dragged_refinement_results = rr;

      if (false) {
         if (rr.refinement_results_contain_overall_nbc_score) {
            std::cout << "-------------- nbc baddies " << std::endl;
            for (unsigned int i=0; i<rr.sorted_nbc_baddies.size(); i++)
               std::cout << "       nbc number " << i
                         << ":  " << rr.sorted_nbc_baddies[i].atom_spec_1
                         << ":  " << rr.sorted_nbc_baddies[i].atom_spec_2
                         << " "  << rr.sorted_nbc_baddies[i].score << std::endl;
         }
      }
      if (false) {
         if (rr.refinement_results_contain_overall_rama_plot_score) {
            std::cout << "-------------- rama baddies " << std::endl;
            for (unsigned int i=0; i<rr.sorted_rama_baddies.size(); i++)
               std::cout << "       rama number " << i
                         << ":  " << rr.sorted_rama_baddies[i].first
                         << " "  << rr.sorted_rama_baddies[i].second << std::endl;
         }
      }

      // std::cout << "Here in refinement_loop_threaded() with rr.progress " << rr.progress << std::endl;

      if (rr.progress == GSL_SUCCESS) {
         continue_update_refinement_atoms_flag = false; // not sure what this does
         rr = g.saved_dragged_refinement_results;
         continue_threaded_refinement_loop = false;
         refinement_has_finished_moving_atoms_representation_update_needed_flag = true;

         // The hooray() function goes off too frequently.
         // Maybe It shouldn't go off if there is no intervention.
         // Unless it's really good.
         //
         std::pair<bool, std::string> hooray = rr.hooray();

         if (hooray.first) {
            // we can't touch Gtk or OpenGL because this we are in a thread
            // (not the main thread)
            // g.setup_draw_for_particles();
            g.setup_draw_for_particles_semaphore = true;
         } else {
            // set up a semaphore or something to put some text
            // (hooray.second) into the status bar
         }

      } else {
         if (rr.progress == GSL_FAILURE) {
            graphics_info_t::continue_update_refinement_atoms_flag = false;
            refinement_has_finished_moving_atoms_representation_update_needed_flag = true;
            continue_threaded_refinement_loop = false;
         } else {
            if (rr.progress == GSL_ENOPROG) {
               graphics_info_t::continue_update_refinement_atoms_flag = false;
               refinement_has_finished_moving_atoms_representation_update_needed_flag = true;
               continue_threaded_refinement_loop = false;
            }
         }
      }

      // 20220503-PE we do this every time the refinement updates - is that right?
      //             No. Becuase this is not the same thread as the graphics thread.
      // update_bad_nbc_atom_pair_marker_positions();

      graphics_info_t::threaded_refinement_loop_counter++;

      if (false)
	 std::cout << "threaded_refinement_loop_counter "
		   << threaded_refinement_loop_counter << std::endl;
   }

   // std::cout << "DEBUG:: refinement_loop_threaded() unlocking restraints_lock" << std::endl;
   release_restraints_lock(__FUNCTION__);
   // std::cout << "debug:: refinement_loop_threaded() goodbye" << std::endl;

   // when this function exits, the (detached) thread in which it's running ends
}

// static
void graphics_info_t::thread_for_refinement_loop_threaded() {

   // I think that there is a race condition here
   // check_and_warn_inverted_chirals_and_cis_peptides()
   // get called several times when the refine loop ends
   // (with success?).


   if (restraints_lock) {
      if (false)
         std::cout << "debug:: thread_for_refinement_loop_threaded() restraints locked by "
                   << restraints_locking_function_name << std::endl;
      return;
   } else {

      if (use_graphics_interface_flag) {

         if (!refinement_immediate_replacement_flag) {

            // if there's not a refinement redraw function already running start up a new one.
            if (threaded_refinement_redraw_timeout_fn_id == -1) {
               GSourceFunc cb = GSourceFunc(regenerate_intermediate_atoms_bonds_timeout_function_and_draw);
	       // int id = gtk_timeout_add(15, cb, NULL);

               int timeout_ms = 15;
               timeout_ms =  30; // 20220503-PE try this value
               timeout_ms = 120; // 20221227-PE try this value - Lucrezia crash
               timeout_ms = 30;
	       int id = g_timeout_add(timeout_ms, cb, NULL);
               threaded_refinement_redraw_timeout_fn_id = id;
            }
         }
      }

      continue_threaded_refinement_loop = true;
      std::thread r(refinement_loop_threaded);
      r.detach();
   }

}

void
graphics_info_t::poke_the_refinement() {

   if (moving_atoms_asc) {
      continue_threaded_refinement_loop = false;
      while (restraints_lock) {
         std::this_thread::sleep_for(std::chrono::milliseconds(2)); // not sure about the delay
      }
      if (last_restraints) {
         double tw = last_restraints->get_torsion_restraints_weight();
         last_restraints->set_map_weight(geometry_vs_map_weight);
         last_restraints->set_torsion_restraints_weight(torsion_restraints_weight);
         last_restraints->set_lennard_jones_epsilon(lennard_jones_epsilon);
         last_restraints->set_geman_mcclure_alpha(geman_mcclure_alpha);
         last_restraints->set_rama_plot_weight(rama_plot_restraints_weight);
         thread_for_refinement_loop_threaded(); // restart refinement if it's not running
      }
   }
}

// static
void graphics_info_t::refinement_of_last_restraints_needs_reset() {

   graphics_info_t::refinement_of_last_restraints_needs_reset_flag = true;

}


void
graphics_info_t::conditionally_wait_for_refinement_to_finish() {

   if (refinement_immediate_replacement_flag || !use_graphics_interface_flag) {
      while (restraints_lock) {
         // this is the main thread - it better be! :-)
	 // std::cout << "conditionally_wait_for_refinement_to_finish() "
         //           << restraints_locking_function_name << std::endl;
         std::this_thread::sleep_for(std::chrono::milliseconds(30));
      }
   }
}




coot::restraint_usage_Flags
graphics_info_t::set_refinement_flags() const {

   coot::restraint_usage_Flags flags = coot::TYPICAL_RESTRAINTS;
   flags = coot::TYPICAL_RESTRAINTS_WITH_IMPROPERS;

   // Oh, these will interact badly.

   if (do_torsion_restraints) {
      flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS;
   }

   if (do_rama_restraints)
      flags = coot::ALL_RESTRAINTS;

   return flags;
}

// reruns refinement if we have restraints
//
// static
void
graphics_info_t::set_geman_mcclure_alpha(float alpha) {

   graphics_info_t g;
   geman_mcclure_alpha = alpha;
   if (g.last_restraints_size() > 0) {
      thread_for_refinement_loop_threaded();
   }
}

// reruns refinement if we have restraints
//
// static
void
graphics_info_t::set_lennard_jones_epsilon(float epsilon) {

   graphics_info_t g;
   lennard_jones_epsilon = epsilon;
   if (g.last_restraints_size() > 0)
      thread_for_refinement_loop_threaded();
}


void
graphics_info_t::update_restraints_with_atom_pull_restraints() {

   // not sure if this should be static or not.

   // maybe we don't want to add... because it's already there at the given position
   //
   for (std::size_t j=0; j<atom_pulls.size(); j++) {
      const atom_pull_info_t &atom_pull = atom_pulls[j];
      if (atom_pull.get_status()) {
         // noise
         // 	 std::cout << "update_refinement_atoms() adding atom_pull_restraint "
         // 		   << atom_pull.spec << std::endl;
	 last_restraints->add_atom_pull_restraint(atom_pull.spec, atom_pull.pos); // mouse target position
      }
   }


   // we don't want to make reference to the moving atoms if they have been deleted.

   if (continue_threaded_refinement_loop) {

      if (auto_clear_atom_pull_restraint_flag) {
         // returns true when the restraint was turned off.
         // turn_off_atom_pull_restraints_when_close_to_target_position() should not
         // include the atom that the use is actively dragging
         //
         // i.e. except this one:
         mmdb::Atom *at_except = 0;
         coot::atom_spec_t except_dragged_atom(at_except);
         if (moving_atoms_currently_dragged_atom_index != -1) {
            if (moving_atoms_asc) {

               // I need the atoms lock here, because we don't want to access the moving atoms
               // if they have been deleted.
               bool unlocked = false;
               while (! moving_atoms_lock.compare_exchange_weak(unlocked, true) && !unlocked) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(1));
                    unlocked = 0;
               }
               if (moving_atoms_asc->atom_selection) {
		  // check that moving_atoms_currently_dragged_atom_index is set before using it
		  if (moving_atoms_currently_dragged_atom_index > -1) {
                     if (moving_atoms_currently_dragged_atom_index < moving_atoms_asc->n_selected_atoms) {
                        at_except = moving_atoms_asc->atom_selection[moving_atoms_currently_dragged_atom_index];
                        except_dragged_atom = coot::atom_spec_t(at_except);
                     }
		  }
               } else {
                  std::cout << "WARNING:: attempted use moving_atoms_asc->atom_selection, but NULL"
                            << std::endl;
               }
               moving_atoms_lock = false;
            }
         }

         std::vector<coot::atom_spec_t> specs_for_removed_restraints =
	    last_restraints->turn_off_atom_pull_restraints_when_close_to_target_position(except_dragged_atom);
         if (specs_for_removed_restraints.size()) {
            if (false) {
               for (unsigned int i=0; i<specs_for_removed_restraints.size(); i++) {
                  std::cout << "INFO:: Clear pull restraint on atom: " << specs_for_removed_restraints[i]
                            << std::endl;
               }
            }
            unsigned int unlocked = false;
            while (! moving_atoms_bonds_lock.compare_exchange_weak(unlocked, 1) && !unlocked) {
                 std::this_thread::sleep_for(std::chrono::milliseconds(1));
                 unlocked = false;
            }

	    atom_pulls_off(specs_for_removed_restraints);
	    clear_atom_pull_restraints(specs_for_removed_restraints, true);
            moving_atoms_bonds_lock = false;
         }
      }
   }
}

// static
gint
graphics_info_t::regenerate_intermediate_atoms_bonds_timeout_function_and_draw(gpointer data) {

   int continue_status = regenerate_intermediate_atoms_bonds_timeout_function();
   graphics_draw();

   if (continue_status == 0) { // if refinement was stopped, do we need to clean up
                               // or accept atoms?

      // ---- Here we can touch the graphics! ---------

      graphics_info_t g; // 37 nanoseconds

      if (threaded_refinement_needs_to_accept_moving_atoms) {
         g.accept_moving_atoms(); // calls clear_up_moving_atoms() which deletes last_restraints
      }

      if (threaded_refinement_needs_to_clear_up) {
         if (false)
            std::cout << "---------- in regenerate_intermediate_atoms_bonds_timeout_function() clear up moving atoms! "
                      << std::endl;
         g.clear_up_moving_atoms(); // get the restraints lock, deletes last_restraints
         g.clear_moving_atoms_object();

         if (glareas[0])
            gtk_widget_remove_tick_callback(glareas[0], wait_for_hooray_refinement_tick_id);

      }

      // no need to do this if Esc is pressed.
      if (! refinement_immediate_replacement_flag)
         g.check_and_warn_inverted_chirals_and_cis_peptides();
   }

   return continue_status;

}

// static
int
graphics_info_t::regenerate_intermediate_atoms_bonds_timeout_function() {

   // Sometimes this thread doesn't stop at the end of a refinement because
   // another mouse drag starts another refinement (and sets restraints_lock
   // to true before we get here). Hmm.

   int continue_status = 1;

   if (! restraints_lock) {
      threaded_refinement_redraw_timeout_fn_id = -1; // we've finished
      continue_status = 0; // stop the timeout function after this redraw
   }

   if (! use_graphics_interface_flag) {
      continue_status = 0;
      return continue_status;
   }

   if (moving_atoms_asc == NULL) {
      // std::cout << "DEBUG:: regenerate_intermediate_atoms_bonds_timeout_function() no moving_atoms_asc"
      // << std::endl;
      continue_status = 0;
      threaded_refinement_redraw_timeout_fn_id = -1; // we've finished
      return continue_status;
   }

   if (moving_atoms_asc->atom_selection == NULL) {
      // std::cout << "DEBUG:: regenerate_intermediate_atoms_bonds_timeout_function() no
      // oving_atoms_asc->atom_selection" << std::endl;
      threaded_refinement_redraw_timeout_fn_id = -1; // we've finished
      continue_status = 0;
      return continue_status;
   }

   // OK, we are interactive

   if (false) {
      if (threaded_refinement_loop_counter_bonds_gen < threaded_refinement_loop_counter) {
         std::cout << "DEBUG:: regenerate_intermediate_atoms_bonds_timeout_function signals make new bonds\n";
      } else {
         std::cout << "DEBUG:: regenerate_intermediate_atoms_bonds_timeout_function bonds - already done\n";
      }
   }

   bool do_the_redraw = false;
   if (threaded_refinement_loop_counter_bonds_gen < threaded_refinement_loop_counter)
      do_the_redraw = true;

   if (refinement_has_finished_moving_atoms_representation_update_needed_flag) {
      do_the_redraw = true;
      refinement_has_finished_moving_atoms_representation_update_needed_flag = false; // reset
   }

   if (do_the_redraw) {

      bool do_rota_markup = graphics_info_t::do_intermediate_atoms_rota_markup;

      // wrap the filling of the rotamer probability tables
      //
      coot::rotamer_probability_tables *tables_pointer = NULL;

      if (do_rota_markup) {
	 if (! rot_prob_tables.tried_and_failed()) {
	    if (rot_prob_tables.is_well_formatted()) {
	       tables_pointer = &rot_prob_tables;
	    } else {
	       rot_prob_tables.fill_tables();
	       if (rot_prob_tables.is_well_formatted()) {
		  tables_pointer = &rot_prob_tables;
	       }
	    }
	 } else {
	    do_rota_markup = false;
	 }
      }

      // we want to see the bonds if
      // threaded_refinement_loop_counter_bonds_gen < threaded_refinement_loop_counter
      // and that does not dependent on if the refinement was running.


      unsigned int unlocked = false;
      while (! moving_atoms_bonds_lock.compare_exchange_weak(unlocked, 1) && !unlocked) {
         std::this_thread::sleep_for(std::chrono::milliseconds(1));
         unlocked = 0;
      }

      // To stop the moving atoms begin deleted as we draw them, we lock them here
      //
      bool unlocked_atoms = false;
      while (! moving_atoms_lock.compare_exchange_weak(unlocked_atoms, 1) && !unlocked_atoms) {
         std::this_thread::sleep_for(std::chrono::milliseconds(1));
         unlocked_atoms = 0;
      }

      // tell the next round that we have (already) drawn the bonds for this set
      // of atom positions
      threaded_refinement_loop_counter_bonds_gen = threaded_refinement_loop_counter;
      graphics_info_t g;
      g.make_moving_atoms_graphics_object(imol_moving_atoms, *moving_atoms_asc);

      g.debug_refinement(); // checks COOT_DEBUG_REFINEMENT

      // Dots on then off but dots remain? Just undisplay them in the Generic Object manager
      if (do_coot_probe_dots_during_refine_flag) {
         g.do_interactive_coot_probe();
         graphics_draw(); // 20220316-PE is this what I want?
      }

      update_bad_nbc_atom_pair_marker_positions();
      update_bad_nbc_atom_pair_dashed_lines();
      update_hydrogen_bond_positions(); // if the intermediate atoms had hydrogen bond restraints

      moving_atoms_bonds_lock = 0;
      moving_atoms_lock = false;

   }

   if (! restraints_lock)
      continue_status = 0; // stop the timeout function

   if (continue_status == 0) {
      // unset the timeout function token
      graphics_info_t::threaded_refinement_redraw_timeout_fn_id = -1; // we've finished
   }

   return continue_status; // return 0 to stop
}

void
graphics_info_t::debug_refinement() {

   // we get the restraints restraints_lock here if needed

   // 20180923 Hmm... I guess that I need to lock the restraints? So that
   // they are not deleted as I try to write them out.
   //
   // Hideous race condition somehow? Yes, there was - when flags was passed to
   // the refinement and it updated the restraints_container_t's internal flags
   //

   bool do_tabulate_geometric_distortions_flag = false;
   char *env = getenv("COOT_DEBUG_REFINEMENT");
   if (env) {
      if (last_restraints) {
         do_tabulate_geometric_distortions_flag = true;
      }
   }

   if (do_debug_refinement)
      do_tabulate_geometric_distortions_flag = true;


   if (do_tabulate_geometric_distortions_flag) {
      if (last_restraints) {
         get_restraints_lock(__FUNCTION__);
         tabulate_geometric_distortions(*last_restraints);
         release_restraints_lock(__FUNCTION__);
      }
   }
}


coot::refinement_results_t
graphics_info_t::refine_molecule(int imol, mmdb::Manager *mol) {

   // This path is needed/is necessary when interactively refining ribosomes

   bool use_map_flag = true;
   coot::refinement_results_t rr = generate_molecule_from_molecule_and_refine(imol, mol, use_map_flag);

   if (rr.found_restraints_flag) {
      graphics_draw();
      if (! refinement_immediate_replacement_flag) {
	 if (use_graphics_interface_flag) {
	    do_accept_reject_dialog("Refinement", rr);
	    check_and_warn_inverted_chirals_and_cis_peptides();
	 }
      }
   }
   return rr;

}

coot::refinement_results_t
graphics_info_t::refine_chain(int imol, const std::string &chain_id, mmdb::Manager *mol) {

   // fill me in at some stage, something along the lines of refine_molecule()
   // Note that make_moving_atoms_asc() is slow if it's done residue by residue
   // So write another version of make_moving_atoms_asc() that takes the whole chain.

   coot::refinement_results_t rr(0, GSL_CONTINUE, "");

   return rr;
}


coot::refinement_results_t
graphics_info_t::refine_residues_vec(int imol,
				     const std::vector<mmdb::Residue *> &residues,
				     const std::string &alt_conf,
				     mmdb::Manager *mol) {

   bool use_map_flag = true;

   if (false)
      std::cout << "INFO:: refine_residues_vec() with altconf \"" << alt_conf << "\"" << std::endl;

   coot::refinement_results_t rr = generate_molecule_and_refine(imol, residues, alt_conf, mol, use_map_flag);

   short int istat = rr.found_restraints_flag;
   if (istat) {
      graphics_draw();
      if (! refinement_immediate_replacement_flag) {
	 if (use_graphics_interface_flag) {
	    do_accept_reject_dialog("Refinement", rr);
	    check_and_warn_inverted_chirals_and_cis_peptides();
	 }
      }
   }
   return rr;
}

coot::refinement_results_t
graphics_info_t::triple_refine_auto_accept() {

   coot::refinement_results_t rr(0, GSL_CONTINUE, "");
   std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = graphics_info_t::active_atom_spec();
   if (active_atom.first) {
      int imol = active_atom.second.first;
      coot::residue_spec_t res_spec(active_atom.second.second);
      molecule_class_info_t &m = molecules[imol];
      mmdb::Residue *r_this = m.get_residue(res_spec);
      if (r_this) {
         std::string alt_conf(active_atom.second.second.alt_conf);
         mmdb::Manager *mol = m.atom_sel.mol;
         float dist_crit = 2.2;
         std::vector<coot::residue_spec_t> rnr = molecules[imol].residues_near_residue(res_spec, dist_crit);
         std::vector<mmdb::Residue *> residues;
         residues.push_back(r_this);
         for (unsigned int i=0; i<rnr.size(); i++) {
            mmdb::Residue *r = m.get_residue(rnr[i]);
            if (r)
               residues.push_back(r);
         }
         int saved_state = refinement_immediate_replacement_flag;
         refinement_immediate_replacement_flag = true;
         refine_residues_vec(imol, residues, alt_conf, mol);
         if (last_restraints) {
#if 0 // 20220423-PE - this crashes. Although get_refinement_results() calls setup_minimize()
      // there is a failure to set *something* up correctly when the restraints and density fit
      // is being calculated.
            coot::refinement_results_t rr_start = last_restraints->get_refinement_results(); // calls setup_minimize()
            conditionally_wait_for_refinement_to_finish();
            std::cout << "---------- pre ------" << std::endl;
            rr_start.show();
            std::cout << "---------- post ------" << std::endl;
            rr.show();
#endif
            accept_moving_atoms();
         }
         refinement_immediate_replacement_flag = saved_state;
      }
   }
   return rr;  // nothing interesting now.
}


coot::refinement_results_t
graphics_info_t::regularize_residues_vec(int imol,
					 const std::vector<mmdb::Residue *> &residues,
					 const std::string &alt_conf,
					 mmdb::Manager *mol) {

   bool use_map_flag = 0;
   coot::refinement_results_t rr = generate_molecule_and_refine(imol, residues, alt_conf, mol, use_map_flag);
   short int istat = rr.found_restraints_flag;
   if (istat) {
      graphics_draw();
      if (! refinement_immediate_replacement_flag) {
	 if (use_graphics_interface_flag) {
	    do_accept_reject_dialog("Regularization", rr);
	    check_and_warn_inverted_chirals_and_cis_peptides();
	 }
      }
   }
   return rr;
}


// copy the functionality of generate_molecule_and_refine() but don't mess about with residues (until you have to)
// Factor out the construction of last_restraints and make_restraints(). Don't repeat that.
//
coot::refinement_results_t
graphics_info_t::generate_molecule_from_molecule_and_refine(int imol, mmdb::Manager *mol, bool use_map_flag) {

   // Note that make_moving_atoms_asc() is slow if it's done residue by residue.
   // So write another version of make_moving_atoms_asc() that takes the whole mol (trivial).

   coot::refinement_results_t rr(0, GSL_CONTINUE, "");

   std::cout << "fill me in" << std::endl;

   // check that if we are using a map that it is valid

   // Need this for the molecule:
//       std::vector<std::string> residue_types = coot::util::residue_types_in_residue_vec(residues);
//       // use try_dynamic_add()
//       bool have_restraints = geom_p->have_dictionary_for_residue_types(residue_types, imol,
// 								       cif_dictionary_read_number);
//       cif_dictionary_read_number += residue_types.size();

   // check that the names of the atoms match the dictionary:
// 	    bool check_hydrogens_too_flag = false;
// 	    std::pair<bool, std::vector<std::pair<std::string, std::vector<std::string> > > >
// 	       icheck_atoms = Geom_p()->atoms_match_dictionary(imol, residues, check_hydrogens_too_flag, false);

   // then create an mmdb::Manager that is the moving molecule
   // make a vector of all the residues in moving molecule (to be used in refinement).

   // set imol_moving_atoms = imol;
   // set moving_atoms_asc_type = coot::NEW_COORDS_REPLACE;

   // save moving_atoms_n_cis_peptides for use with the dialog later:
   // int n_cis = coot::util::count_cis_peptides(local_moving_atoms_asc.mol);
   // moving_atoms_n_cis_peptides = n_cis;

   // std::vector<coot::atom_spec_t> fixed_atom_specs = molecules[imol].get_fixed_atoms();

   // then call the (factored out in the future): last_restraints = new restraints_container_t(...)

   return rr;
}

#include "ligand/rotamer.hh"

std::vector<std::pair<mmdb::Residue *, std::vector<coot::dict_torsion_restraint_t> > >
graphics_info_t::make_rotamer_torsions(const std::vector<std::pair<bool, mmdb::Residue *> > &local_residues) const {

   std::vector<std::pair<mmdb::Residue *, std::vector<coot::dict_torsion_restraint_t> > > v;
   for (unsigned int i=0; i<local_residues.size(); i++) {
      if (! local_residues[i].first) {
         mmdb::Residue *residue_p = local_residues[i].second;
         std::string rn(residue_p->GetResName());
         if (coot::util::is_standard_amino_acid_name(rn)) {
            std::string alt_conf; // run through them all, ideally.
            coot::rotamer rot(residue_p, alt_conf, 1);
            coot::closest_rotamer_info_t cri = rot.get_closest_rotamer(rn);
            if (cri.residue_chi_angles.size() > 0) {
               std::vector<coot::dict_torsion_restraint_t> dictionary_vec;
               std::vector<std::vector<std::string> > rotamer_atom_names = rot.rotamer_atoms(rn);

               if (cri.residue_chi_angles.size() != rotamer_atom_names.size()) {

                  std::cout << "-------------- mismatch for " << coot::residue_spec_t(residue_p) << " "
                            << cri.residue_chi_angles.size() << " "  << rotamer_atom_names.size()
                            << " ---------------" << std::endl;
               } else {

                  for (unsigned int ichi=0; ichi<cri.residue_chi_angles.size(); ichi++) {
                     // we have to convert chi angles to atom names
                     double esd = 3.0; // 20210315-PE was 10.0. I want them tighter than that.
                     int per = 1;
                     std::string id = "chi " + coot::util::int_to_string(cri.residue_chi_angles[ichi].first);
                     const std::string &atom_name_1 = rotamer_atom_names[ichi][0];
                     const std::string &atom_name_2 = rotamer_atom_names[ichi][1];
                     const std::string &atom_name_3 = rotamer_atom_names[ichi][2];
                     const std::string &atom_name_4 = rotamer_atom_names[ichi][3];
                     double torsion = cri.residue_chi_angles[ichi].second;
                     coot::dict_torsion_restraint_t dr(id, atom_name_1, atom_name_2, atom_name_3, atom_name_4,
                                                       torsion, esd, per);
                     dictionary_vec.push_back(dr);
                  }

                  if (dictionary_vec.size() > 0) {
                     std::pair<mmdb::Residue *, std::vector<coot::dict_torsion_restraint_t> > p(residue_p, dictionary_vec);
                     v.push_back(p);
                  }
               }
            }
         }
      }
   }
   return v;
}


// return the state of having found restraints.
bool
graphics_info_t::make_last_restraints(const std::vector<std::pair<bool,mmdb::Residue *> > &local_residues,
				      const std::vector<mmdb::Link> &links,
				      const coot::protein_geometry &geom,
				      mmdb::Manager *mol_for_residue_selection,
				      const std::vector<coot::atom_spec_t> &fixed_atom_specs,
				      coot::restraint_usage_Flags flags,
				      bool use_map_flag,
				      const clipper::Xmap<float> *xmap_p) {

   if (last_restraints) {
      std::cout << "----------------------------------------------" << std::endl;
      std::cout << "----------------------------------------------" << std::endl;
      std::cout << "    ERROR:: A: last_restraints not cleared up " << std::endl;
      std::cout << "----------------------------------------------" << std::endl;
      std::cout << "----------------------------------------------" << std::endl;
   }

   if (false) { // these are the passed residues, nothing more.
      std::cout << "debug:: on construction of restraints_container_t local_residues: "
		<< std::endl;
      for (std::size_t jj=0; jj<local_residues.size(); jj++) {
	 std::cout << "   " << coot::residue_spec_t(local_residues[jj].second)
		   << " is fixed: " << local_residues[jj].first << std::endl;
      }
   }

   moving_atoms_extra_restraints_representation.clear();
   continue_threaded_refinement_loop = true; // no longer set in refinement_loop_threaded()

   // the refinment of torsion seems a bit confused? If it's in flags, why does it need an flag
   // of its own? I suspect that it doesn't. For now I will keep it (as it was).
   //
   bool do_residue_internal_torsions = false;
   if (do_torsion_restraints) {
      do_residue_internal_torsions = 1;
   }

   last_restraints = new
      coot::restraints_container_t(local_residues,
				   links,
				   *Geom_p(),
				   mol_for_residue_selection,
				   fixed_atom_specs, xmap_p);

   // std::cout << "debug:: on creation last_restraints is " << last_restraints << std::endl;

   bool verbose_refinement_geometry_reporting = false; // make this user-setable
   if (! verbose_refinement_geometry_reporting)
      last_restraints->set_quiet_reporting();

   last_restraints->set_torsion_restraints_weight(torsion_restraints_weight);

   if (convert_dictionary_planes_to_improper_dihedrals_flag) {
      last_restraints->set_convert_plane_restraints_to_improper_dihedral_restraints(true);
   }

   // This seems not to work yet.
   // last_restraints->set_dist_crit_for_bonded_pairs(9.0);

   if (use_map_flag)
      last_restraints->add_map(geometry_vs_map_weight);

   unsigned int n_threads = coot::get_max_number_of_threads();
   if (n_threads > 0)
      last_restraints->thread_pool(&static_thread_pool, n_threads);

   if (false)
      std::cout << "---------- debug:: in generate_molecule_and_refine() "
		<< " calling restraints.make_restraints() with imol "
		<< imol_moving_atoms << " "
		<< molecules[imol_moving_atoms].extra_restraints.has_restraints()
		<<  " " << molecules[imol_moving_atoms].extra_restraints.bond_restraints.size()
		<< std::endl;

   all_atom_pulls_off();
   particles_have_been_shown_already_for_this_round_flag = false;
   if (use_graphics_interface_flag)
      if (glareas[0])
         wait_for_hooray_refinement_tick_id =
            gtk_widget_add_tick_callback(glareas[0], wait_for_hooray_refinement_tick_func, 0, 0);
   //
   // elsewhere do this:
   // gtk_widget_remove_tick_callback(glareas[0], wait_for_hooray_refinement_tick_id);

   moving_atoms_visited_residues.clear(); // this is used for HUD label colour

   int n_restraints = last_restraints->make_restraints(imol_moving_atoms,
						       *Geom_p(), flags,
						       do_residue_internal_torsions,
						       do_trans_peptide_restraints,
						       rama_plot_restraints_weight,
						       do_rama_restraints,
						       true, true, make_auto_h_bond_restraints_flag,
						       pseudo_bonds_type);
                                                       // link and flank args default true

   if (use_harmonic_approximation_for_NBCs) {
      std::cout << "INFO:: using soft harmonic restraints for NBC" << std::endl;
      last_restraints->set_use_harmonic_approximations_for_nbcs(true);
   }

   if (pull_restraint_neighbour_displacement_max_radius > 1.99) {
      last_restraints->set_use_proportional_editing(true);
      last_restraints->pull_restraint_neighbour_displacement_max_radius =
         pull_restraint_neighbour_displacement_max_radius;
   }

   last_restraints->set_geman_mcclure_alpha(geman_mcclure_alpha);
   last_restraints->set_lennard_jones_epsilon(graphics_info_t::lennard_jones_epsilon);
   last_restraints->set_rama_type(restraints_rama_type);
   last_restraints->set_rama_plot_weight(rama_plot_restraints_weight); // >2? danger of non-convergence
                                                                       // if planar peptide restraints are used
   // Oh, I see... it's not just the non-Bonded contacts of the hydrogens.
   // It's the planes, chiral and angles too. Possibly bonds too.
   // How about marking non-H atoms in restraints that contain H atoms as
   // "invisible"? i.e. non-H atoms are not influenced by the positions of the
   // Hydrogen atoms (but Hydrogen atoms *are* influenced by the positions of the
   // non-Hydrogen atoms). This seems like a lot of work. Might be easier to turn
   // off angle restraints for H-X-X (but not H-X-H) in the first instance, that
   // should go most of the way to what "invisible" atoms would do, I imagine.
   // is_H_non_bonded_contact should be renamed to is_H_turn_offable_restraint
   // or something.
   //
   // last_restraints->set_apply_H_non_bonded_contacts(false);

   if (do_rotamer_restraints) {
      std::vector<std::pair<mmdb::Residue *, std::vector<coot::dict_torsion_restraint_t> > > rotamer_torsions = make_rotamer_torsions(local_residues);
      std::cout << "debug:: calling add_or_replace_torsion_restraints_with_closest_rotamer_restraints() from make_last_restraints() " << std::endl;
      last_restraints->add_or_replace_torsion_restraints_with_closest_rotamer_restraints(rotamer_torsions);
   }

   if (molecules[imol_moving_atoms].extra_restraints.has_restraints()) {
      std::cout << "debug:: calling add_extra_restraints() from make_last_restraints() " << std::endl;
      last_restraints->add_extra_restraints(imol_moving_atoms, "user-defined from make_last_restraints()",
                                            molecules[imol_moving_atoms].extra_restraints, *Geom_p());
   }

   if (do_numerical_gradients)
      last_restraints->set_do_numerical_gradients();

   bool found_restraints_flag = false;

   if (last_restraints->size() > 0) {

      last_restraints->analyze_for_bad_restraints();
      thread_for_refinement_loop_threaded();
      found_restraints_flag = true;
      // rr.found_restraints_flag = true;
      draw_gl_ramachandran_plot_flag = true;

      // are you looking for conditionally_wait_for_refinement_to_finish() ?

      if (refinement_immediate_replacement_flag) {
         // wait until refinement finishes
         while (restraints_lock) {
            std::this_thread::sleep_for(std::chrono::milliseconds(7));
            std::cout << "INFO:: make_last_restraints() [immediate] restraints locked by "
                      << restraints_locking_function_name << std::endl;
         }
      }

   } else {
      continue_threaded_refinement_loop = false;
      if (use_graphics_interface_flag) {
         // GtkWidget *widget = create_no_restraints_info_dialog();
         GtkWidget *widget = widget_from_builder("no_restraints_info_dialog");
         gtk_widget_set_visible(widget, TRUE);
      }
   }

   return found_restraints_flag;
}


// simple mmdb::Residue * interface to refinement.  20081216
//
coot::refinement_results_t
graphics_info_t::generate_molecule_and_refine(int imol,
					      const std::vector<mmdb::Residue *> &residues_in,
					      const std::string &alt_conf,
					      mmdb::Manager *mol,
					      bool use_map_flag) {

   if (false) {
      std::cout << "debug:: ----- start generate_molecule_and_refine() " << std::endl;
      for (unsigned int i=0; i<residues_in.size(); i++) {
	 std::cout << "   residue " << i << " of residues_in is " << residues_in[i] << std::endl;
      }
   }


   auto tp_0 = std::chrono::high_resolution_clock::now();

   coot::refinement_results_t rr(0, GSL_CONTINUE, "");

   if (is_valid_map_molecule(Imol_Refinement_Map()) || (! use_map_flag)) {
      // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
      coot::restraint_usage_Flags flags = set_refinement_flags();
      bool do_residue_internal_torsions = false;
      if (do_torsion_restraints) {
	 do_residue_internal_torsions = 1;
	 flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS;
      }

      if (do_rama_restraints)
	 // flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_RAMA;
	 flags = coot::ALL_RESTRAINTS;

      std::vector<coot::atom_spec_t> fixed_atom_specs = molecules[imol].get_fixed_atoms();

      // refinement goes a bit wonky if there are multiple occurrances of the same residue
      // in input residue vector, so let's filter out duplicates here
      //
      std::vector<mmdb::Residue *> residues;
      std::set<mmdb::Residue *> residues_set;
      std::set<mmdb::Residue *>::const_iterator it;
      for (std::size_t i=0; i<residues_in.size(); i++)
	 residues_set.insert(residues_in[i]);
      residues.reserve(residues_set.size());
      for(it=residues_set.begin(); it!=residues_set.end(); ++it)
	 residues.push_back(*it);

      // OK, so the passed residues are the residues in the graphics_info_t::molecules[imol]
      // molecule.  We need to do 2 things:
      //
      // convert the mmdb::Residue *s of the passed residues to the mmdb::Residue *s of residues mol
      //
      // and
      //
      // in create_mmdbmanager_from_res_vector() make sure that that contains the flanking atoms.
      // The create_mmdbmanager_from_res_vector() from this class is used, not coot::util
      //
      // The flanking atoms are fixed, the passed residues are not fixed.
      // Keep a clear head.

      std::vector<std::string> residue_types = coot::util::residue_types_in_residue_vec(residues);
      // use try_dynamic_add()
      bool have_restraints = geom_p->have_restraints_dictionary_for_residue_types(residue_types, imol,
                                                                                  cif_dictionary_read_number);

      cif_dictionary_read_number += residue_types.size();

      if (have_restraints) {

	 std::string residues_alt_conf = alt_conf;
	 imol_moving_atoms = imol;
	 std::pair<mmdb::Manager *, std::vector<mmdb::Residue *> > residues_mol_and_res_vec =
	    create_mmdbmanager_from_res_vector(residues, imol, mol, residues_alt_conf);

	 if (false) { // debug
	    mmdb::Manager *residues_mol = residues_mol_and_res_vec.first;
	    int imod = 1;
	    mmdb::Model *model_p = residues_mol->GetModel(imod);
	    if (model_p) {
	       int n_chains = model_p->GetNumberOfChains();
	       for (int ichain=0; ichain<n_chains; ichain++) {
		  mmdb::Chain *chain_p = model_p->GetChain(ichain);
		  std::cout << "residues_mol_and_res_vec mol: chain: " << chain_p->GetChainID() << std::endl;
		  int nres = chain_p->GetNumberOfResidues();
		  for (int ires=0; ires<nres; ires++) {
		     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
		     if (residue_p) {
			std::cout << "residues_mol_and_res_vec mol:   residue "
				  << coot::residue_spec_t(residue_p) << " residue "
				  << residue_p << " chain " << residue_p->chain << " index "
				  << residue_p->index << std::endl;
		     } else {
			std::cout << "residues_mol_and_res_vec: residue " << ires << " was null"
				  << std::endl;
		     }
		  }
	       }
	    }
	 }

	 // We only want to act on these new residues and molecule, if
	 // there is something there.
	 //
	 if (residues_mol_and_res_vec.first) {

	    // Now we want to do an atom name check.  This stops exploding residues.
	    //
	    bool check_hydrogens_too_flag = false;
	    std::pair<bool, std::vector<std::pair<mmdb::Residue *, std::vector<std::string> > > >
	       icheck_atoms = Geom_p()->atoms_match_dictionary(imol, residues, check_hydrogens_too_flag, false);

	    if (! icheck_atoms.first) {

	       // Oops. Just give us a dialog and don't start the refinement
	       info_dialog_refinement_non_matching_atoms(icheck_atoms.second); // protected

	    } else {

	       moving_atoms_have_hydrogens_displayed = true;
	       if (! molecules[imol].draw_hydrogens())
		  moving_atoms_have_hydrogens_displayed = false;

	       atom_selection_container_t local_moving_atoms_asc =
		  make_moving_atoms_asc(residues_mol_and_res_vec.first, residues);

               make_moving_atoms_graphics_object(imol, local_moving_atoms_asc);

               int n_cis = coot::util::count_cis_peptides(local_moving_atoms_asc.mol);
               moving_atoms_n_cis_peptides = n_cis;

	       std::vector<std::pair<bool,mmdb::Residue *> > local_residues;  // not fixed.
	       for (unsigned int i=0; i<residues_mol_and_res_vec.second.size(); i++)
		  local_residues.push_back(std::pair<bool, mmdb::Residue *>(0, residues_mol_and_res_vec.second[i]));

	       moving_atoms_asc_type = coot::NEW_COORDS_REPLACE;

	       int imol_for_map = Imol_Refinement_Map();
	       clipper::Xmap<float> *xmap_p = dummy_xmap;

	       if (is_valid_map_molecule(imol_for_map))
		  xmap_p = &molecules[imol_for_map].xmap;

	       bool found_restraints_flag = make_last_restraints(local_residues,
								 local_moving_atoms_asc.links,
								 *Geom_p(),
								 residues_mol_and_res_vec.first,
								 fixed_atom_specs,
								 flags, use_map_flag, xmap_p);

               if (last_restraints) {
                  // 20220423-PE I can't do this here because setup_minimize() has not been called yet
                  // rr = last_restraints->get_refinement_results();
               }
	       rr.found_restraints_flag = found_restraints_flag;

	    }
	 }
      } else {
	 // we didn't have restraints for everything.
	 //
	 std::pair<int, std::vector<std::string> > icheck =
	    check_dictionary_for_residue_restraints(imol, residues);
	 if (icheck.first == 0) {
	    // info_dialog_missing_refinement_residues(icheck.second);
	    show_missing_refinement_residues_dialog(icheck.second, false); // now had use_graphics_interface_flag protection
	 }
      }
   }

   if (false) {
      auto tp_1 = std::chrono::high_resolution_clock::now();
      auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
      std::cout << "INFO:: ---------- Timing for refinement " << d10 << " milliseconds" << std::endl;
   }

   return rr;

}

#include "coot-utils/atom-tools.hh"

// mol is new (not from molecules[imol]) molecule for the moving atoms.
//
// resno_1 and resno_2 need to be passed because in the
// conventional/linear case, we don't want all the residues selected
// (sigh). That's because in that case the residues_mol that is passed
// also has the flanking residues.
//
// Consider setting those to be RES_ANY for other cases.
//
atom_selection_container_t
graphics_info_t::make_moving_atoms_asc(mmdb::Manager *residues_mol,
				       int resno_1,
				       int resno_2) const {

   atom_selection_container_t local_moving_atoms_asc;
   if (residues_mol->GetNumberOfModels() > 0) {
      local_moving_atoms_asc.mol = residues_mol;
      local_moving_atoms_asc.UDDOldAtomIndexHandle = -1;  // true?
      local_moving_atoms_asc.UDDAtomIndexHandle = -1;
      if (residues_mol)
         local_moving_atoms_asc.read_success = 1;

      local_moving_atoms_asc.SelectionHandle = residues_mol->NewSelection();
      residues_mol->SelectAtoms (local_moving_atoms_asc.SelectionHandle, 0, "*",
				 resno_1, // starting resno, an int
				 "*", // any insertion code
				 resno_2, // ending resno
				 "*", // ending insertion code
				 "*", // any residue name
				 "*", // atom name
				 "*", // elements
				 "*"  // alt loc.
				 );

      residues_mol->GetSelIndex(local_moving_atoms_asc.SelectionHandle,
                                local_moving_atoms_asc.atom_selection,
                                local_moving_atoms_asc.n_selected_atoms);

      local_moving_atoms_asc.fill_links_using_mol(residues_mol);

   }
   return local_moving_atoms_asc;
}


atom_selection_container_t
graphics_info_t::make_moving_atoms_asc(mmdb::Manager *residues_mol,
                                       const std::vector<mmdb::Residue *> &residues) const {

   // This also rebonds the imol_moving_atoms molecule

   GLenum err;
   if (use_graphics_interface_flag) {
      err = glGetError();
      if (err) std::cout << "GL ERROR:: in make_moving_atoms_asc()-- start --\n"; // hmm were was thig from?
   }

   atom_selection_container_t local_moving_atoms_asc;
   local_moving_atoms_asc.UDDAtomIndexHandle = -1;
   local_moving_atoms_asc.UDDOldAtomIndexHandle = residues_mol->GetUDDHandle(mmdb::UDR_ATOM, "old atom index");

   int SelHnd = residues_mol->NewSelection();

   for (unsigned int ir=0; ir<residues.size(); ir++) {
      const char *chain_id = residues[ir]->GetChainID();
      const char *inscode = residues[ir]->GetInsCode();
      int resno = residues[ir]->GetSeqNum();
      residues_mol->Select(SelHnd, mmdb::STYPE_ATOM,
			   0, chain_id,
			   resno, // starting resno, an int
			   inscode, // any insertion code
			   resno, // ending resno
			   inscode, // ending insertion code
			   "*", // any residue name
			   "*", // atom name
			   "*", // elements
			   "*",  // alt loc.
			   mmdb::SKEY_OR);
   }

   local_moving_atoms_asc.mol = residues_mol;
   local_moving_atoms_asc.SelectionHandle = SelHnd;
   residues_mol->GetSelIndex(local_moving_atoms_asc.SelectionHandle,
			     local_moving_atoms_asc.atom_selection,
			     local_moving_atoms_asc.n_selected_atoms);


   if (false) {
      std::cout << "returning an atom selection for all moving atoms "
                << local_moving_atoms_asc.n_selected_atoms << " atoms "
                << std::endl;
   }

   // This new block added so that we don't draw atoms in the "static" molecule when we have the
   // corresponding atoms in the moving atoms.
   //
   const atom_selection_container_t &imol_asc = molecules[imol_moving_atoms].atom_sel;
   std::set<int> atom_set = coot::atom_indices_in_other_molecule(imol_asc, local_moving_atoms_asc);

   if (false) { // debug atoms in other molecule
      std::set<int>::const_iterator it;
      for (it=atom_set.begin(); it!=atom_set.end(); ++it) {
         int idx = *it;
         mmdb::Atom *at = imol_asc.atom_selection[idx];
         coot::atom_spec_t as(at);
         std::cout << " this is a moving atom: " << idx << " " << as << std::endl;
      }
   }

   if (false) { // debug old atom index
      for (int i=0; i<local_moving_atoms_asc.n_selected_atoms; i++) {
          mmdb::Atom *at = local_moving_atoms_asc.atom_selection[i];
          coot::atom_spec_t as(at);
          int idx = -1;
          at->GetUDData(local_moving_atoms_asc.UDDOldAtomIndexHandle, idx);
          std::cout << "DEBUG:: in make_moving_atoms_asc " << as << " idx " << idx << std::endl;
      }
   }

   if (use_graphics_interface_flag) {
      err = glGetError();
      if (err) std::cout << "GL ERROR:: in make_moving_atoms_asc()-- before end --\n";
   }

   // now rebond molecule imol without bonds to atoms in atom_set
   if (atom_set.size() > 0) {
      if (regenerate_bonds_needs_make_bonds_type_checked_flag) {
         molecules[imol_moving_atoms].make_bonds_type_checked(atom_set, __FUNCTION__);
      }
   }

   return local_moving_atoms_asc;
}

// surely we need to have control of one of the locks before we can do this?
// 20200220-PE - Sensible comment - I just discovered a crash here.
//
void
graphics_info_t::make_moving_atoms_restraints_graphics_object() {

   // 20230823-PE do we need any control over this (extra bond restraints) other than
   // the user_control?
   // Let's see...
   draw_it_for_moving_atoms_restraints_graphics_object = true;

   if (moving_atoms_asc) {
      if (last_restraints) {
         if (draw_it_for_moving_atoms_restraints_graphics_object) {
            if (draw_it_for_moving_atoms_restraints_graphics_object_user_control) {
               moving_atoms_extra_restraints_representation.clear();

               for (int i=0; i<last_restraints->size(); i++) {
                  // std::cout << "-------------------- in make_moving_atoms_restraints_graphics_object() D " << i << std::endl;
                  const coot::simple_restraint &rest = last_restraints->at(i);
                  if (rest.restraint_type == coot::BOND_RESTRAINT ||
                      rest.restraint_type == coot::GEMAN_MCCLURE_DISTANCE_RESTRAINT) {

                     if (rest.target_value > 2.15) {  // no real bond restraints
                        int idx_1 = rest.atom_index_1;
                        int idx_2 = rest.atom_index_2;
                        // we can't display bonds to non-moving atoms
                        if (idx_1 < moving_atoms_asc->n_selected_atoms) {
                           if (idx_2 < moving_atoms_asc->n_selected_atoms) {
                              mmdb::Atom *at_1 = moving_atoms_asc->atom_selection[idx_1];
                              mmdb::Atom *at_2 = moving_atoms_asc->atom_selection[idx_2];
                              if (at_1 && at_2) {
                                 clipper::Coord_orth p1 = coot::co(at_1);
                                 clipper::Coord_orth p2 = coot::co(at_2);
                                 // dd is the target value
                                 // de is the actual value
                                 const double &dd = rest.target_value;
                                 float def = sqrtf(clipper::Coord_orth(p1-p2).lengthsq());
                                 double de = static_cast<double>(def);
                                 bool do_it = true;
                                 std::string atom_name_1 = at_1->GetAtomName();
                                 std::string atom_name_2 = at_2->GetAtomName();
                                 if (atom_name_1 == " CA ")
                                    if (atom_name_2 == " CA ")
                                       do_it = false;
                                 if (do_it)
                                    moving_atoms_extra_restraints_representation.add_bond(p1, p2, dd, de);
                              }
                           }
                        }
                     }
                  }
               }

               if (false)
                  std::cout << "in make_moving_atoms_restraints_graphics_object calling make_extra_distance_restraints_objects() "
                            << std::endl;
               // now call the new graphics function:
               make_extra_distance_restraints_objects(); // make graphics objecs from moving_atoms_extra_restraints_representation()
            }
         }
      }
   }

}

// static
void
graphics_info_t::draw_moving_atoms_restraints_graphics_object() {

   std::cout << "FIXME in draw_moving_atoms_restraints_graphics_object() " << std::endl;

   // 20221221-PE now I have draw_extra_distance_restraints();

#if 0

   if (draw_it_for_moving_atoms_restraints_graphics_object) {
      if (moving_atoms_asc) {
         if (last_restraints) {
            if (moving_atoms_extra_restraints_representation.bonds.size() > 0) {
               glLineWidth(1.0); // 2 is smoother and fatter

               // maybe these can be overriddn by the molecules:
               // extra_restraints_representation.prosmart_restraint_display_limit_low
               // extra_restraints_representation.prosmart_restraint_display_limit_high
               //
               double z_lim_low  = -1.0;
               double z_lim_high =  1.0;

               glBegin(GL_LINES);
               for (unsigned int ib=0; ib<moving_atoms_extra_restraints_representation.bonds.size(); ib++) {

                  const coot::extra_restraints_representation_t::extra_bond_restraints_respresentation_t &res =
                     moving_atoms_extra_restraints_representation.bonds[ib];

                  // purple if actual distance is greater than target
                  //
                  double d_sqd = (res.second - res.first).clipper::Coord_orth::lengthsq();
                  double d = sqrt(d_sqd);
                  double esd = 0.05;
                  double z = (res.target_dist - d)/esd;
                  if (z > z_lim_low) {
                     if (z > z_lim_high) {
                        double b = 0.02 * z;
                        if (b >  0.4999) b =  0.4999;
                        if (b < -0.4999) b = -0.4999;
                        double b_green = b;
                        glColor3d(0.5-b, 0.5+b_green*0.9, 0.5-b);

                        glVertex3f(res.first.x(), res.first.y(), res.first.z());
                        glVertex3f(res.second.x(), res.second.y(), res.second.z());
                     }
                  }
               }
               glEnd();
            }
         }
      }
   }
#endif
}


// Return 0 (first) if any of the residues don't have a dictionary
// entry and a list of the residue type that don't have restraints.
//
std::pair<int, std::vector<std::string> >
graphics_info_t::check_dictionary_for_residue_restraints(int imol, mmdb::PResidue *SelResidues, int nSelResidues) {

   // Are you sure you want this function?

   bool status_OK = 1; // pass, by default
   std::vector<std::string> res_name_vec;

   for (int ires=0; ires<nSelResidues; ires++) {
      std::string resn(SelResidues[ires]->GetResName());
      std::string resname = adjust_refinement_residue_name(resn);
      bool status = geom_p->have_dictionary_for_residue_type(resname, imol, cif_dictionary_read_number);
      cif_dictionary_read_number++;
      if (! status) {
	 status_OK = false;
	 res_name_vec.push_back(resname);
      }

      if (false)
	 std::cout << "DEBUG:: have_dictionary_for_residues() on residue "
		   << ires << " of " << nSelResidues << ", "
		   << resname << " returns "
		   << status << std::endl;
      cif_dictionary_read_number++;
   }
   return std::pair<int, std::vector<std::string> > (status_OK, res_name_vec);
}

// Return 0 (first) if any of the residues don't have a dictionary
// entry and a list of the residue type that don't have restraints.
//
std::pair<int, std::vector<std::string> >
graphics_info_t::check_dictionary_for_residue_restraints(int imol, const std::vector<mmdb::Residue *> &residues) {

   int status = 1;
   std::vector<std::string> missings;
   std::vector<std::string> types;
   for (unsigned int i=0; i<residues.size(); i++)
      types.push_back(residues[i]->GetResName());

   if (false)
     for (unsigned int i=0; i<types.size(); i++)
       std::cout << " type :" << types[i] << ":" << std::endl;

   bool r_status = geom_p->have_restraints_dictionary_for_residue_types(types, imol, cif_dictionary_read_number++);
   if (! r_status) {
      // OK, compile a list of the missings
      for (const auto &rn : types) {
         std::vector<std::string> v = {rn};
         bool status_inner = geom_p->have_restraints_dictionary_for_residue_types(v, imol, cif_dictionary_read_number++);
         if (! status_inner) {
            missings.push_back(rn);
            status = 0;
         }
      }
   }
   return std::make_pair(status, missings);
}


std::string
graphics_info_t::adjust_refinement_residue_name(const std::string &resname) const {

   std::string r = resname;
   if (resname == "UNK") r = "ALA"; // hack for KC/buccaneer.
   if (resname.length() > 2)
      if (resname[2] == ' ')
	 r = resname.substr(0,2);
   return r;
}





// The flanking residues (if any) are in the residue selection (SelResidues).
// The flags are not needed now we have made adjustments in the calling
// function.
//
// create_mmdbmanager_from_res_selection must make adjustments
//
// Note: there is now a molecule-class-info version of this - perhaps
// we should call it?  Next bug fix here: move over to the function call.
//
// deep copy the passed residues
//
mmdb::Manager *
graphics_info_t::create_mmdbmanager_from_res_selection(mmdb::PResidue *SelResidues,
						       int nSelResidues,
						       int have_flanking_residue_at_start,
						       int have_flanking_residue_at_end,
						       const std::string &altconf,
						       const std::string &chain_id_1,
						       short int residue_from_alt_conf_split_flag,
						       int imol) {

   int start_offset = 0;
   int end_offset = 0;

//    if (have_flanking_residue_at_start)
//       start_offset = -1;
//    if (have_flanking_residue_at_end)
//       end_offset = +1;

   mmdb::Manager *residues_mol = new mmdb::Manager;
   mmdb::Model *model = new mmdb::Model;
   mmdb::Chain *chain = new mmdb::Chain;
   short int whole_res_flag = 0; // not all alt confs, only this one ("A") and "".

   // For the active residue range (i.e. not the flanking residues) we only want
   // to refine the atoms that have the alt conf the same as the picked atom
   // (and that is altconf, passed here).
   //
   // However, for *flanking residues* it's different.  Say we are refining a
   // non-split residue with alt conf "": Say that residue has a flanking
   // residue that is completely split, into A and B.  In that case we want
   // either "" or "A" for the flanking atoms.
   //
   // And say we want to refine the A alt conf of a completely split residue
   // that has a flanking neighbour that is completely unsplit (""), we want
   // atoms that are either "A" or "".
   //
   // So let's try setting whole_res_flag to 1 for flanking residues.

   mmdb::Residue *r;
   int atom_index_udd = molecules[imol].atom_sel.UDDAtomIndexHandle;
   for (int ires=start_offset; ires<(nSelResidues + end_offset); ires++) {

      if ( (ires == 0) || (ires == nSelResidues -1) ) {
	 if (! residue_from_alt_conf_split_flag)
	    whole_res_flag = 1;
      } else {
	 whole_res_flag = 0;
      }

//       std::cout << "DEBUG in create_mmdbmanager_from_res_selection, whole_res_flag is "
// 		<< whole_res_flag << " for altconf " << altconf
// 		<< "\n       residue_from_alt_conf_split_flag "
// 		<< residue_from_alt_conf_split_flag << std::endl;

      bool embed_in_chain_flag = false; // don't put r in a chain in deep_copy_this_residue()
                                        // because we put r in a chain here.
      r = coot::deep_copy_this_residue_old_style(SelResidues[ires], altconf, whole_res_flag,
                                                 atom_index_udd, embed_in_chain_flag);
      if (r) {
	 chain->AddResidue(r);
	 r->seqNum = SelResidues[ires]->GetSeqNum();
	 r->SetResName(SelResidues[ires]->GetResName());
      }
   }
   chain->SetChainID(chain_id_1.c_str());
   model->AddChain(chain);
   residues_mol->AddModel(model);
   residues_mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
   residues_mol->FinishStructEdit();

   return residues_mol;
}

// called by simple_refine_residues (a refinement from a vector of mmdb::Residues).
//
// The returned mol should have flanking residues too.
//
// return a NULL in the first of the pair if the past residue vector is of size 0.
//
std::pair<mmdb::Manager *, std::vector<mmdb::Residue *> >
graphics_info_t::create_mmdbmanager_from_res_vector(const std::vector<mmdb::Residue *> &residues,
						    int imol,
						    mmdb::Manager *mol_in,
						    std::string alt_conf) {

   // returned entities
   mmdb::Manager *new_mol = 0;
   std::vector<mmdb::Residue *> rv; // gets checked

   float dist_crit = 5.0;
   bool debug = false;

   if (debug) {
      std::cout << "############ starting create_mmdbmanager_from_res_vector() with these "
		<< " residues " << std::endl;
      for (std::size_t ii=0; ii<residues.size(); ii++)
	 std::cout << "   " << coot::residue_spec_t(residues[ii])  << std::endl;
      int udd_atom_index_handle = mol_in->GetUDDHandle(mmdb::UDR_ATOM, "atom index");
      std::cout << "############ udd for atom index from seeding molecule " << udd_atom_index_handle
		<< std::endl;
      for (std::size_t ii=0; ii<residues.size(); ii++) {
	 mmdb::Residue *residue_p = residues[ii];
	 mmdb::Atom **residue_atoms = 0;
	 int n_residue_atoms;
	 residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (int iat=0; iat<n_residue_atoms; iat++) {
	    mmdb::Atom *at = residue_atoms[iat];
	    int idx = -1;
	    at->GetUDData(udd_atom_index_handle, idx);
	    std::cout << "#### input residue atom " << coot::atom_spec_t(at) << " had udd index "
		      << idx << std::endl;
	 }
      }
   }

   int n_flanker = 0; // a info/debugging counter

   if (residues.size() > 0) {

      // Also add the index of the reference residue (the one in molecules[imol].atom_selection.mol)
      // to the molecule that we are construction here. So that we can properly link
      // the residues in restraints_container (there we rather need to know the references indices,
      // not the indices from the fragment molecule)
      //

      std::pair<bool,std::string> use_alt_conf(false, "");
      if (! alt_conf.empty())
	 use_alt_conf = std::pair<bool, std::string> (true, alt_conf);

      if (false) { // debug
         std::cout << "----------------- in create_mmdbmanager_from_res_vector() alt_conf is "
                   << "\"" << alt_conf << "\"" << std::endl;
         std::cout << "----------------- in create_mmdbmanager_from_res_vector() use_alt_conf is "
                   << use_alt_conf.first << "\"" << use_alt_conf.second << "\"" << std::endl;
      }

      std::pair<bool, mmdb::Manager *> n_mol_1 =
	 coot::util::create_mmdbmanager_from_residue_vector(residues, mol_in, use_alt_conf);

      // check that first is sane, so indent all this lot (when it works)

      if (n_mol_1.first) {

	 int index_from_reference_residue_handle =
	    n_mol_1.second->GetUDDHandle(mmdb::UDR_RESIDUE, "index from reference residue");

	 if (false) { // debug
	    int imod = 1;
	    mmdb::Model *model_p = n_mol_1.second->GetModel(imod);
	    if (model_p) {
	       int n_chains = model_p->GetNumberOfChains();
	       for (int ichain=0; ichain<n_chains; ichain++) {
		  mmdb::Chain *chain_p = model_p->GetChain(ichain);
		  int nres = chain_p->GetNumberOfResidues();
		  for (int ires=0; ires<nres; ires++) {
		     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
		     int n_atoms = residue_p->GetNumberOfAtoms();
		     for (int iat=0; iat<n_atoms; iat++) {
			mmdb::Atom *at = residue_p->GetAtom(iat);
			int idx = -1;
			at->GetUDData(index_from_reference_residue_handle, idx);
			std::cout << "   create_mmdbmanager_from_residue_vector() returns this mol atom "
				  << iat << " " << coot::atom_spec_t(at) << " with idx " << idx << std::endl;
		     }
		  }
	       }
	    }
	 }

	 new_mol = n_mol_1.second;
	 mmdb::Model *model_p = new_mol->GetModel(1);

	 // how many (free) residues were added to that model? (add them to rv)
	 //
	 int n_chains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<n_chains; ichain++) {
	    mmdb::Chain *chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    for (int ires=0; ires<nres; ires++) {
	       mmdb::Residue *residue_p = chain_p->GetResidue(ires);
	       rv.push_back(residue_p);
	    }
	 }

	 if (false) {
	    for (std::size_t ir=0; ir<rv.size(); ir++) {
	       mmdb::Residue *r = rv[ir];
	       std::cout << "Moving Residue " << coot::residue_spec_t(r) << std::endl;
	       mmdb::Atom **residue_atoms = 0;
	       int n_residue_atoms;
	       r->GetAtomTable(residue_atoms, n_residue_atoms);
	       for (int iat=0; iat<n_residue_atoms; iat++) {
		  mmdb::Atom *at = residue_atoms[iat];
		  std::cout << "    " << coot::atom_spec_t(at) << std::endl;
	       }
	    }
	 }

	 short int whole_res_flag = 0;
	 int atom_index_udd_handle = molecules[imol].atom_sel.UDDAtomIndexHandle;

	 // Now the flanking residues:
	 //
	 std::vector<mmdb::Residue *> flankers_in_reference_mol;
	 flankers_in_reference_mol.reserve(32); // say

	 // find the residues that are close to the residues of
	 // residues that are not part of residues
	 //
	 // We don't have quite the function that we need in coot-utils,
	 // so we need to munge residues in to local_residues:
	 std::vector<std::pair<bool, mmdb::Residue *> > local_residues;
	 local_residues.resize(residues.size());
	 for (std::size_t ires=0; ires<residues.size(); ires++)
	    local_residues[ires] = std::pair<bool, mmdb::Residue *>(false, residues[ires]);
	 std::map<mmdb::Residue *, std::set<mmdb::Residue *> > rnr =
	    coot::residues_near_residues(local_residues, mol_in, dist_crit);
	 // now fill @var{flankers_in_reference_mol} from rnr, avoiding residues
	 // already in @var{residues}.
	 std::map<mmdb::Residue *, std::set<mmdb::Residue *> >::const_iterator it;
	 for (it=rnr.begin(); it!=rnr.end(); ++it) {
	    const std::set<mmdb::Residue *> &s = it->second;
	    std::set<mmdb::Residue *>::const_iterator its;
	    for (its=s.begin(); its!=s.end(); ++its) {
	       mmdb::Residue *tres = *its;
	       if (std::find(residues.begin(), residues.end(), tres) == residues.end())
		  if (std::find(flankers_in_reference_mol.begin(), flankers_in_reference_mol.end(), tres) == flankers_in_reference_mol.end())
		     flankers_in_reference_mol.push_back(tres);
	    }
	 }

	 // So we have a vector of residues that were flankers in the
	 // reference molecule, we need to add copies of those to
	 // new_mol (making sure that they go into the correct chain).
	 //
	 if (false) { // debug
	    std::cout << "debug:: ############ Found " << flankers_in_reference_mol.size()
		      << " flanking residues" << std::endl;

	    for (unsigned int ires=0; ires<flankers_in_reference_mol.size(); ires++)
	       std::cout << "     #### flankers_in_reference_mol: " << ires << " "
			 << coot::residue_spec_t(flankers_in_reference_mol[ires]) << std::endl;
	 }


	 for (unsigned int ires=0; ires<flankers_in_reference_mol.size(); ires++) {
	    mmdb::Residue *r;

	    std::string ref_res_chain_id = flankers_in_reference_mol[ires]->GetChainID();

	    mmdb::Chain *chain_p = NULL;
	    int n_new_mol_chains = model_p->GetNumberOfChains();
	    for (int ich=0; ich<n_new_mol_chains; ich++) {
	       if (ref_res_chain_id == model_p->GetChain(ich)->GetChainID()) {
		  chain_p = model_p->GetChain(ich);
		  break;
	       }
	    }

	    if (! chain_p) {
	       // Add a new one then.
	       chain_p = new mmdb::Chain;
	       chain_p->SetChainID(ref_res_chain_id.c_str());
	       model_p->AddChain(chain_p);
	    }

	    if (false)
	       std::cout << "debug:: flankers_in_reference_mol " << ires << " "
			 << coot::residue_spec_t(flankers_in_reference_mol[ires]) << " "
			 << "had index " << flankers_in_reference_mol[ires]->index
			 << std::endl;

            // get rid of this function at some stage
            bool embed_in_chain = false;
	    r = coot::deep_copy_this_residue_old_style(flankers_in_reference_mol[ires],
					               alt_conf, whole_res_flag,
					               atom_index_udd_handle, embed_in_chain);

	    if (r) {

	       r->PutUDData(index_from_reference_residue_handle, flankers_in_reference_mol[ires]->index);

	       // copy over the atom indices. UDDAtomIndexHandle in mol_n becomes UDDOldAtomIndexHandle
	       // indices in the returned molecule

	       int sni = find_serial_number_for_insert(r->GetSeqNum(), r->GetInsCode(), chain_p);

	       if (false) { // debug
		  mmdb::Atom **residue_atoms = 0;
		  int n_residue_atoms;
		  std::cout << "Flanker Residue " << coot::residue_spec_t(r) << std::endl;
		  r->GetAtomTable(residue_atoms, n_residue_atoms);
		  for (int iat=0; iat<n_residue_atoms; iat++) {
		     mmdb::Atom *at = residue_atoms[iat];
		     std::cout << "    " << coot::atom_spec_t(at) << std::endl;
		  }
	       }

	       if (sni == -1)
		  chain_p->AddResidue(r); // at the end
	       else
		  chain_p->InsResidue(r, sni);
	       r->seqNum = flankers_in_reference_mol[ires]->GetSeqNum();
	       r->SetResName(flankers_in_reference_mol[ires]->GetResName());
	       n_flanker++;

	       if (false)
		  std::cout << "debug:: create_mmdbmanager_from_residue_vector() inserted/added flanker "
			    << coot::residue_spec_t(r) << std::endl;

	    }
	 }

	 // super-critical for correct peptide bonding in refinement!
	 //
	 coot::util::pdbcleanup_serial_residue_numbers(new_mol);

	 if (debug) {
	    int imod = 1;
	    mmdb::Model *model_p = new_mol->GetModel(imod);
	    if (model_p) {
	       int n_chains = model_p->GetNumberOfChains();
	       for (int ichain=0; ichain<n_chains; ichain++) {
		  mmdb::Chain *chain_p = model_p->GetChain(ichain);
		  int nres = chain_p->GetNumberOfResidues();
		  for (int ires=0; ires<nres; ires++) {
		     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
		     std::cout << "create_mmdb..  ^^^ " << coot::residue_spec_t(residue_p) << " "
			       << residue_p << " index " << residue_p->index
			       << std::endl;
		  }
	       }
	    }
	 }

	 if (debug)
	    std::cout << "DEBUG:: in create_mmdbmanager_from_res_vector: " << rv.size()
		      << " free residues and " << n_flanker << " flankers" << std::endl;
      }
   }

   return std::pair <mmdb::Manager *, std::vector<mmdb::Residue *> > (new_mol, rv);
}

// return -1 on failure to find a residue for insertion index
//
int
graphics_info_t::find_serial_number_for_insert(int seqnum_new,
					       const std::string &ins_code_for_new,
					       mmdb::Chain *chain_p) const {

   int iserial_no = -1;
   int current_diff = 999999;
   if (chain_p) {
      int nres = chain_p->GetNumberOfResidues();
      for (int ires=0; ires<nres; ires++) { // ires is a serial number
	 mmdb::Residue *residue = chain_p->GetResidue(ires);

	 // we are looking for the smallest negative diff:
	 //
	 int diff = residue->GetSeqNum() - seqnum_new;
	 if ( (diff > 0) && (diff < current_diff) ) {
	    iserial_no = ires;
	    current_diff = diff;
	 } else {
	    if (diff == 0) {
	       std::string ins_code_this = residue->GetInsCode();
	       if (ins_code_this > ins_code_for_new) {
		  iserial_no = ires;
		  break;
	       }
	    }
	 }
      }
   }
   return iserial_no;
}


// on reading a pdb file, we get a list of residues, use these to
// load monomers from the dictionary, to be used in refinement.
//
// This could be substantially speeded up, we currently do a most
// simple-minded search...  If both vectors were sorted on the
// monomer/residue name we could speed up a lot.
//
int
graphics_info_t::load_needed_monomers(const std::vector<std::string> &pdb_residue_types) {

   int iloaded=0;

   for (unsigned int ipdb=0; ipdb<pdb_residue_types.size(); ipdb++) {
      bool ifound = 0;

      int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
      if (geom_p->have_dictionary_for_residue_type_no_dynamic_add(pdb_residue_types[ipdb], imol_enc))
	 ifound = true;

      if (! ifound) {
	 // read in monomer for type pdb_residue_types[ipdb]
	 geom_p->try_dynamic_add(pdb_residue_types[ipdb],
				 cif_dictionary_read_number++);
	 iloaded++;
      }
   }
   return iloaded;
}



coot::refinement_results_t
graphics_info_t::regularize(int imol, short int auto_range_flag, int i_atom_no_1, int i_atom_no_2) {

   // 20230111-PE is this function used?

   // What are we going to do here:
   //
   // How do we get the atom selection (the set of atoms that will be
   // refined)?  Just do a SelectAtoms() (see notes on atom selection).
   //
   // Next flash the selection (to be refined).
   // How do we get the bonds for that?
   //
   // Next, convert the atom selection to the output of refmac.  The
   // output of refmac has a different atom indexing.  We need to
   // construct a ppcatom with the same indexing as the refmac file.
   // We will do that as part of the idealization class (currently
   // simple-restraints)
   //

   coot::refinement_results_t rr;
   int tmp;
   if (i_atom_no_1 > i_atom_no_2) {
      tmp = i_atom_no_1;
      i_atom_no_1 = i_atom_no_2;
      i_atom_no_2 = tmp;
   }
   // now i_atom_no_2 is greater than i_atom_no_1.

   // cout << "regularize: molecule " << imol << " atom index " << i_atom_no_1
   // << " to " << i_atom_no_2 << endl;

   int resno_1, resno_2;

   mmdb::PPAtom SelAtom = molecules[imol].atom_sel.atom_selection;

   resno_1 = SelAtom[i_atom_no_1]->residue->seqNum;
   resno_2 = SelAtom[i_atom_no_2]->residue->seqNum;

   std::string inscode_1 = SelAtom[i_atom_no_1]->residue->GetInsCode();
   std::string inscode_2 = SelAtom[i_atom_no_2]->residue->GetInsCode();

   if (resno_1 > resno_2) {
      tmp = resno_1;
      resno_1 = resno_2;
      resno_2 = tmp;
      std::string tmp_ins = inscode_2;
      inscode_2 = inscode_1;
      inscode_1 = tmp_ins;
   }

   std::string chain_id_1(SelAtom[i_atom_no_1]->residue->GetChainID());
   std::string chain_id_2(SelAtom[i_atom_no_2]->residue->GetChainID());
   std::string altconf(SelAtom[i_atom_no_2]->altLoc);

   if (auto_range_flag) {
      // we want the residues that are +/- 1 [typically] from the residues that
      // contains the atom with the index i_atom_no_1.
      std::pair<int, int> p = auto_range_residues(i_atom_no_1, imol);
      resno_1 = p.first;
      resno_2 = p.second;
      // std::cout << "DEBUG:: auto_range_residues: " << resno_1 << " " << resno_2 << std::endl;
   }

   // if ( chain_id_1 != chain_id_2 ) {
      // pointer comparison:
   if (SelAtom[i_atom_no_1]->GetChain() != SelAtom[i_atom_no_1]->GetChain()) {
      std::cout << "Picked atoms are not in the same chain.  Failure" << std::endl;
      std::cout << "FYI: chain ids are: \"" << chain_id_1
		<< "\" and \"" << chain_id_2 << "\"" << std::endl;
      std::cout << "Picked atoms are not in the same chain.  Failure" << std::endl;
   } else {
      // flash_selection(imol, resno_1, inscode_1, resno_2, inscode_2, altconf, chain_id_1);
      rr = copy_mol_and_regularize(imol, resno_1, inscode_1, resno_2, inscode_2, altconf, chain_id_1);
      short int istat = rr.found_restraints_flag;
      if (istat) {
	 graphics_draw();
	 if (! refinement_immediate_replacement_flag) {
	    // std::cout << "DEBUG:: Regularize: rr.info is " << rr.info << std::endl;
	    if (use_graphics_interface_flag) {
	       do_accept_reject_dialog("Regularization", rr);
	       check_and_warn_inverted_chirals_and_cis_peptides();
	    }
	 }
      } else {
	 std::cout << "No restraints: regularize()\n";
      }
   } // same chains test
   return rr;
}

std::pair<int, int>
graphics_info_t::auto_range_residues(int atom_index, int imol) const {
   std::pair<int, int> r;

   mmdb::Atom *this_atom =  molecules[imol].atom_sel.atom_selection[atom_index];
   mmdb::Residue *this_res = this_atom->residue;
   mmdb::Chain *this_chain = this_res->chain;
   int resno = this_res->GetSeqNum();
   char *inscode = this_res->GetInsCode();

   mmdb::Residue *prev_res = this_chain->GetResidue(resno-refine_auto_range_step, inscode);
   mmdb::Residue *next_res = this_chain->GetResidue(resno+refine_auto_range_step, inscode);

   // Warning: Enabling this code will cause a crash if prev_res or next_res are NULL.
//    std::cout << " debug:: in auto_range_residues() returns residues "
// 	     << prev_res->GetSeqNum() << " and " << next_res->GetSeqNum()
// 	     << " given refine_auto_range_step " << refine_auto_range_step
// 	     << std::endl;

   if (prev_res) {
      r.first = resno-refine_auto_range_step;
   } else {
      r.first = resno;
   }

   if (next_res) {
      r.second = resno+refine_auto_range_step;
   } else {
      r.second = resno;
   }

   return r;
}

#ifdef USE_GUILE
SCM
graphics_info_t::refinement_results_to_scm(const coot::refinement_results_t &rr) const {

   SCM r = SCM_BOOL_F;

   if (rr.found_restraints_flag) {
      SCM lights_scm = SCM_EOL;
      SCM progress_scm = scm_from_int(rr.progress);
      SCM info_scm = scm_from_locale_string(rr.info_text.c_str());
      for (int il=rr.lights.size()-1; il>=0; il--) {
	 SCM light_scm = SCM_EOL;
	 SCM value_scm = scm_from_double(rr.lights[il].value);
	 SCM label_scm = scm_from_locale_string(rr.lights[il].label.c_str());
 	 SCM  name_scm = scm_from_locale_string(rr.lights[il].name.c_str());

	 light_scm = scm_cons(value_scm, light_scm);
	 light_scm = scm_cons(label_scm, light_scm);
	 light_scm = scm_cons( name_scm, light_scm);

	 lights_scm = scm_cons(light_scm, lights_scm);
      }
      r = SCM_EOL;
      r = scm_cons(lights_scm,r);
      r = scm_cons(progress_scm,r);
      r = scm_cons(info_scm,r);
   }
   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *
graphics_info_t::refinement_results_to_py(const coot::refinement_results_t &rr) const {

   PyObject *r = Py_False;

   if (rr.found_restraints_flag) {
      PyObject *lights_py = Py_False;
      PyObject *progress_py = PyLong_FromLong(rr.progress);
      PyObject *info_py = myPyString_FromString(rr.info_text.c_str());
      if (rr.lights.size())
	lights_py = PyList_New(rr.lights.size());
      for (unsigned int il=0; il<rr.lights.size(); il++) {
	PyObject *light_py = PyList_New(3);
	PyObject *value_py = PyFloat_FromDouble(rr.lights[il].value);
	PyObject *label_py = myPyString_FromString(rr.lights[il].label.c_str());
	PyObject *name_py  = myPyString_FromString(rr.lights[il].name.c_str());

	PyList_SetItem(light_py, 0, name_py);
	PyList_SetItem(light_py, 1, label_py);
	PyList_SetItem(light_py, 2, value_py);

	PyList_SetItem(lights_py, il, light_py);
      }
      r = PyList_New(3);
      PyList_SetItem(r, 0, info_py);
      PyList_SetItem(r, 1, progress_py);
      PyList_SetItem(r, 2, lights_py);
   }

   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif



// static
void
graphics_info_t::flash_position(const clipper::Coord_orth &pos) {

   if (glareas[0]) {
      int n_flash = residue_selection_flash_frames_number; // default 3
      flash_intermediate_atom_pick_flag = 1;
      intermediate_flash_point = pos;
      for (int iflash=0; iflash<n_flash; iflash++) {
	 graphics_draw();
      }
      flash_intermediate_atom_pick_flag = 0;
   }
}


coot::refinement_results_t
graphics_info_t::refine(int imol, short int auto_range_flag, int i_atom_no_1, int i_atom_no_2) {

   coot::refinement_results_t rr;

   int tmp;
   if (i_atom_no_1 > i_atom_no_2) {
      tmp = i_atom_no_1;
      i_atom_no_1 = i_atom_no_2;
      i_atom_no_2 = tmp;
   }
   // now i_atom_no_2 is greater than i_atom_no_1.

   if (! is_valid_model_molecule(imol)) {
      std::cout << "ERROR:: bad molecule number " << imol << std::endl;
      return rr;
   }
   if (i_atom_no_1 < 0) {
      std::cout << "ERROR:: bad atom index 1 " << i_atom_no_1 << std::endl;
      return rr;
   }
   if (i_atom_no_2 < 0) {
      std::cout << "ERROR:: bad atom index 2 " << i_atom_no_2 << std::endl;
      return rr;
   }
   if (i_atom_no_2 >= molecules[imol].atom_sel.n_selected_atoms) {
      std::cout << "out of range atom 2 " << i_atom_no_2 << " vs "
		<< molecules[imol].atom_sel.n_selected_atoms<< std::endl;
      return rr;
   }

   int resno_1, resno_2;

   int imol_map = Imol_Refinement_Map();
   if (imol_map == -1) { // magic number check,
      // if not -1, then it has been set by user

      show_select_map_frame();

   } else {

      mmdb::PPAtom SelAtom = molecules[imol].atom_sel.atom_selection;

      resno_1 = SelAtom[i_atom_no_1]->GetSeqNum();
      resno_2 = SelAtom[i_atom_no_2]->GetSeqNum();

//       std::cout << "DEBUG:: refine: atom1 " << SelAtom[i_atom_no_1] << std::endl;
//       std::cout << "DEBUG:: refine: atom2 " << SelAtom[i_atom_no_2] << std::endl;

      if (auto_range_flag) {
	 // we want the residues that are +/- 1 [typically] from the residues that
	 // contains the atom with the index i_atom_no_1.
	 std::pair<int, int> p = auto_range_residues(i_atom_no_1, imol);
	 resno_1 = p.first;
	 resno_2 = p.second;
      }


      std::string chain_id_1(SelAtom[i_atom_no_1]->residue->GetChainID());
      std::string chain_id_2(SelAtom[i_atom_no_2]->residue->GetChainID());
      std::string altconf(SelAtom[i_atom_no_2]->altLoc);
      short int is_water_like_flag = 0;
      std::string resname_1(SelAtom[i_atom_no_1]->GetResName());
      std::string resname_2(SelAtom[i_atom_no_2]->GetResName());
      std::string inscode_1(SelAtom[i_atom_no_1]->GetInsCode());
      std::string inscode_2(SelAtom[i_atom_no_2]->GetInsCode());

      if (resno_1 > resno_2) {
	 tmp = resno_1;
	 resno_1 = resno_2;
	 resno_2 = tmp;
	 std::string tmp_ins = inscode_2;
	 inscode_2 = inscode_1;
	 inscode_1 = tmp_ins;
      }

      is_water_like_flag = check_for_no_restraints_object(resname_1, resname_2);
      if (! is_water_like_flag)
	 if (SelAtom[i_atom_no_1]->GetResidue() == SelAtom[i_atom_no_2]->GetResidue())
	    is_water_like_flag = check_for_single_hetatom(SelAtom[i_atom_no_2]->GetResidue());
      rr = refine_residue_range(imol, chain_id_1, chain_id_2, resno_1, inscode_1,
				resno_2, inscode_2, altconf, is_water_like_flag);
   }
   return rr;
}

coot::refinement_results_t
graphics_info_t::get_refinement_results() const {

   coot::refinement_results_t rr;

   std::this_thread::sleep_for(std::chrono::milliseconds(20)); // 20220423-PE why is this sleeping at all?
                                                               // was 200 ms.

   if (last_restraints)
      rr = last_restraints->get_refinement_results();

   return rr;
}


// I mean things like HOH, CL, BR etc
bool
graphics_info_t::check_for_no_restraints_object(std::string &resname_1, std::string &resname_2) const {

   // a better check would be to check in the geom for a dictionary of
   // that name and see if there are bonds for that residue type.

   bool r = 0;
   if (resname_1 == "WAT" || resname_1 == "HOH" ||
       resname_2 == "WAT" || resname_2 == "HOH")
      r = 1;
   if (resname_1 == "BR" || resname_1 == "CL" ||
       resname_2 == "BR" || resname_2 == "CL")
      r = 1;
   if (resname_1 == "NA" || resname_1 == "CA" ||
       resname_2 == "NA" || resname_2 == "CA")
      r = 1;
   if (resname_1 == "K" || resname_1 == "MG" ||
       resname_2 == "K" || resname_2 == "MG")
      r = 1;
   return r;

}

// I suppose that if check_for_no_restraints_object() was fully
// featured (i.e. it checked the restraints), we wouldn't need this
// function.
//
// We also check for Metal atom.
//
bool
graphics_info_t::check_for_single_hetatom(mmdb::Residue *res_p) const {

   bool r = 0;

   int n_atoms = res_p->GetNumberOfAtoms();
   if (n_atoms == 1) {
      mmdb::PPAtom residue_atoms;
      int nResidueAtoms;
      res_p->GetAtomTable(residue_atoms, nResidueAtoms);
      if (residue_atoms[0]->Het)
	 r = 1;
      if (residue_atoms[0]->isMetal())
	 r = 1;
   }
   return r;
}


// The calling function need to check that if chain_id_1 and
// chain_id_2 are not the same chain (mmdb::Chain *), then we don't call
// this function.  We don't want to do an atom selection here (we can
// do that in copy_mol_and_refine), so we need to pass is_water_like_flag
// and the auto_range is determined by the calling function.  Here we
// are passed the results of any auto_range calculation.
//
coot::refinement_results_t
graphics_info_t::refine_residue_range(int imol,
				      const std::string &chain_id_1,
				      const std::string &chain_id_2,
				      int resno_1,
				      const std::string &ins_code_1,
				      int resno_2,
				      const std::string &ins_code_2,
				      const std::string &altconf,
				      short int is_water_like_flag) {

   // 20230111-PE is this function used? Try deleting it one rainy day.

   if (false)
      std::cout << "DEBUG:: ================ refine_residue_range: "
		<< imol << " " << chain_id_1
		<< " " <<  resno_1 << ":" << ins_code_1 << ":"
		<< " " <<  resno_2 << ":" << ins_code_2 << ":"
		<< " " << ":" << altconf << ": " << is_water_like_flag << std::endl;

   coot::refinement_results_t rr;

   int imol_map = Imol_Refinement_Map();

   if (imol_map == -1) { // magic number check,
      // if not -1, then it has been set by user
      show_select_map_frame();

   } else {

      // if ( chain_id_1 != chain_id_2 ) {
      // Used to be pointer comparison, let that be done in the calling function.
      if (chain_id_1 != chain_id_2) {

	 // for now we will bug out.  In future, we will want to be
	 // able to refine glycosylation.
	 //
	 std::cout << "Picked atoms are not in the same chain.  Failure" << std::endl;
	 std::cout << "FYI: chain ids are: \"" << chain_id_1
		   << "\" and \"" << chain_id_2 << "\"" << std::endl;
      } else {
	 if (molecules[imol_map].has_xmap()) {  // it may have been
					        // closed after it was
					        // selected.
	    short int simple_water = 0;
	    if (resno_1 == resno_2) {
	       if (is_water_like_flag) {
		  simple_water = 1;
// 		  std::string s = "That water molecule doesn't have restraints.\n";
// 		  s += "Using Rigid Body Fit Zone";
// 		  GtkWidget *w = wrapped_nothing_bad_dialog(s);
// 		  gtk_widget_set_visible(w, TRUE);

		  // rigid body refine uses residue_range_atom_index_1
		  // and residue_range_atom_index_2, which should be
		  // defined in graphics-info-defines refine section.
		  //
		  imol_rigid_body_refine = imol;  // Fix GSK water refine crash.

		  // OK, in the simple water case, we do an atom selection

// 		  std::cout << "DEBUG:: before set_residue_range_refine_atoms altconf is :"
// 			    << altconf << ":" << std::endl;

		  set_residue_range_refine_atoms(chain_id_1, resno_1, resno_2, altconf, imol);
		  // There are now set by that function:
		  // residue_range_atom_index_1 = i_atom_no_1;
		  // residue_range_atom_index_2 = i_atom_no_1; // refining just one atom.

		  execute_rigid_body_refine(0); // no autorange for waters.
	       }
	    }
	    if (!simple_water) {
	       // flash_selection(imol, resno_1, ins_code_1, resno_2, ins_code_2, altconf, chain_id_1);
	       long t0 = 0; // glutGet(GLUT_ELAPSED_TIME);
               if (resno_2 < resno_1) std::swap(resno_1, resno_2);
	       rr = copy_mol_and_refine(imol, imol_map, resno_1, ins_code_1, resno_2, ins_code_2,
					altconf, chain_id_1);
	       short int istat = rr.found_restraints_flag;
	       // long t1 = 0; // glutGet(GLUT_ELAPSED_TIME);
               // std::cout << "Fix glut timer: Refinement elapsed time: " << float(t1-t0)/1000.0 << std::endl;
	       if (istat) {
		  graphics_draw();
		  if (! refinement_immediate_replacement_flag) {
		     if (use_graphics_interface_flag) {
			do_accept_reject_dialog("Refinement", rr);
			check_and_warn_inverted_chirals_and_cis_peptides();
		     }
		  }
	       }
	    }
	 } else {
	    std::cout << "Can't refine to a closed map.  Choose another map"
		      << std::endl;
	    show_select_map_frame();
	 }
      }
   } // same chains test
   return rr;
}

void
graphics_info_t::repeat_refine_zone() {

   if (is_valid_model_molecule(residue_range_mol_no)) {
      refine(residue_range_mol_no, false, residue_range_atom_index_1, residue_range_atom_index_2);
   }

}



// Question to self: Are you sure that imol_rigid_body_refine (the
// coordinates molecule) is set when we get here?
//
// Also: are residue_range_atom_index_1 and residue_range_atom_index_2 set?
// They should be.
//
void
graphics_info_t::execute_rigid_body_refine(short int auto_range_flag) {

   /* Atom picking has happened. Actually do it */

   mmdb::Atom *atom1 = NULL;
   mmdb::Atom *atom2 = NULL;

   int ires1;  // set according to auto_range_flag
   int ires2;
   char *chain_id_1;
   char *chain_id_2 = 0;
   bool mask_water_flag = 0; // don't mask waters

   if (auto_range_flag) {
      std::pair<int, int> p = auto_range_residues(residue_range_atom_index_1,
						  imol_rigid_body_refine);
      ires1 = p.first;
      ires2 = p.second;
      atom1 = molecules[imol_rigid_body_refine].atom_sel.atom_selection[residue_range_atom_index_1];
      chain_id_1 =  atom1->residue->GetChainID();
   } else {

      // make sure that the atom indices are in the right order:
      //
      if (residue_range_atom_index_1 > residue_range_atom_index_2) {
	 int tmp;
	 tmp = residue_range_atom_index_2;
	 residue_range_atom_index_2 = residue_range_atom_index_1;
	 residue_range_atom_index_1 = tmp;
      }

//       std::cout << "imol_rigid_body_refine " << imol_rigid_body_refine << std::endl;
//       std::cout << "residue_range_atom_index_1 "
// 		<< residue_range_atom_index_1 << std::endl;
//       std::cout << "residue_range_atom_index_2 "
// 		<< residue_range_atom_index_2 << std::endl;

      atom1 = molecules[imol_rigid_body_refine].atom_sel.atom_selection[residue_range_atom_index_1];
      atom2 = molecules[imol_rigid_body_refine].atom_sel.atom_selection[residue_range_atom_index_2];
      ires1 = atom1->residue->seqNum;
      ires2 = atom2->residue->seqNum;
      chain_id_1 =  atom1->residue->GetChainID();
      chain_id_2 =  atom2->residue->GetChainID();
      std::string resname_1(atom1->GetResName());
      std::string resname_2(atom2->GetResName());
      if (resname_1 == "WAT" || resname_1 == "HOH" ||
	  resname_2 == "WAT" || resname_2 == "HOH")
	 mask_water_flag = 1; // if the zone is (or contains) a water, then mask waters as if they
                              // were protein atoms.
   }

   // duck out now if the chains were not the same!
   if (chain_id_1 != chain_id_2) {
      std::string info_string("Atoms must be in the same chain");
      add_status_bar_text(info_string);
      return;
   }

   std::string altconf = atom1->altLoc;
   bool select_altconf = true; // only when refining a single ligand/residue or both atom
                               // have the same alt conf that is non-blank.

   if (atom1 && atom2) {
      if (ires1 != ires2) {
	 std::string alt_conf_1 = atom1->altLoc;
	 std::string alt_conf_2 = atom2->altLoc;
	 if (alt_conf_1 != alt_conf_2) {
	    select_altconf = false;
	 } else {
	    if (alt_conf_1.empty())
	       select_altconf = false;
	 }
      }
   }


   std::string chain(chain_id_1);

//    std::cout << "-----------------------------------------------------" << std::endl;
//    std::cout << "-----------------------------------------------------" << std::endl;
//    std::cout << " Rigid Body Refinement "
// 	     << " imol: " << imol_rigid_body_refine << " residue "
// 	     << ires1 << " to " << ires2 << " chain " << chain << std::endl;

   int imol_ref_map = Imol_Refinement_Map();  // -1 is a magic number

   if (Imol_Refinement_Map() == -1 ) { // magic number
      //
      std::cout << "Please set a map against which the refinement should occur"
		<< std::endl;
      show_select_map_frame();  // protected
   } else {

      coot::minimol::molecule mol(molecules[imol_rigid_body_refine].atom_sel.mol);

      coot::minimol::molecule range_mol;
      int ir = range_mol.fragment_for_chain(chain);

      // Fill range_mol and manipulate mol so that it has a blank (it
      // will get copied and used as to mask the map).
      //
      for (unsigned int ifrag=0; ifrag<mol.fragments.size(); ifrag++) {
	 if (mol[ifrag].fragment_id == chain) {
	    for (int ires=mol.fragments[ifrag].min_res_no();
		 ires<=mol.fragments[ifrag].max_residue_number();
		 ires++) {
	       if (ires>=ires1 && ires<=ires2) {

		  // a vector of atoms that should be deleted from the
		  // reference mol so that the density where the
		  // moving atoms are is not deleted.
		  std::vector<unsigned int> from_ref_delete_atom_indices;

		  // a vector of atoms that should be deleted from the
		  // moving mol, because they don't match the altconf
		  std::vector<unsigned int> from_mov_delete_atom_indices;

		  try { // (yes, belt and braces in this case)
		     if (! mol[ifrag][ires].is_undefined()) {
			// std::cout << "adding residue mol[" << ifrag<< "][" << ires << "]: "
			// << mol[ifrag][ires] << std::endl;
			range_mol[ir].addresidue(mol[ifrag][ires], 1);
		     }
		  }
		  catch (const std::runtime_error &rte) {
		     std::cout << "ERROR:: execute_rigid_body_refine() " << rte.what() << std::endl;
		  }


		  for (unsigned int iat=0; iat<mol[ifrag][ires].atoms.size(); iat++) {
		     if ((mol[ifrag][ires][iat].altLoc == altconf) || !select_altconf) {
// 			std::cout << "From ref res delete atom "
// 				  << mol[ifrag][ires][iat] << std::endl;
			from_ref_delete_atom_indices.push_back(iat);
		     } else {
// 			std::cout << "From mov res delete atom "
// 				  << range_mol[ifrag][ires][iat] << std::endl;
			from_mov_delete_atom_indices.push_back(iat);
		     }
		  }
// 		  std::cout << "--------------------------------" << std::endl;
// 		  mol.check();
// 		  mol[ifrag].check();
// 		  std::cout << "--------------------------------" << std::endl;
// 		  range_mol.check();
// 		  range_mol[ir].check();
// 		  std::cout << "--------------------------------" << std::endl;
          	     mol[ifrag][ires].delete_atom_indices(from_ref_delete_atom_indices);
		  range_mol[ir][ires].delete_atom_indices(from_mov_delete_atom_indices);
	       }
	    }
	 }
      }
      coot::minimol::molecule mol_without_moving_zone = mol;
      rigid_body_fit(mol_without_moving_zone, range_mol, imol_ref_map, mask_water_flag);
   } // valid map test
}

// replacing atom positions in imol_rigid_body_refine, so make sure
// that you set that correctly before calling this function.
bool
graphics_info_t::rigid_body_fit(const coot::minimol::molecule &mol_without_moving_zone,
				const coot::minimol::molecule &range_mol,
				int imol_ref_map,
				bool mask_water_flag) {

   bool success = false; // fail initially
   bool debug = false;

   if (! is_valid_map_molecule(imol_ref_map)) {
      std::cout << "WARNING:: not a valid map: " << imol_ref_map << std::endl;
      return success;// false
   }

   std::vector<coot::minimol::atom *> range_atoms = range_mol.select_atoms_serial();

   // debugging
   if (debug) {
      range_mol.write_file("rigid-body-range-mol.pdb", 44);
      mol_without_moving_zone.write_file("rigid-body-without-moving-zone.pdb", 44);
   }

   coot::ligand lig;
   lig.import_map_from(molecules[imol_ref_map].xmap,
		       molecules[imol_ref_map].map_sigma());

   lig.install_ligand(range_mol);
   lig.find_centre_by_ligand(0); // don't test ligand size
   // lig.set_map_atom_mask_radius(0.5);  // why did I have this?
   lig.mask_map(mol_without_moving_zone, mask_water_flag);
   if (debug)
      lig.output_map("rigid-body.map");
   lig.set_dont_write_solutions();
   lig.set_dont_test_rotations();
   lig.set_acceptable_fit_fraction(graphics_info_t::rigid_body_fit_acceptable_fit_fraction);
   lig.fit_ligands_to_clusters(1);
   unsigned int iclust = 0;
   unsigned int isol   = 0;
   coot::minimol::molecule moved_mol = lig.get_solution(iclust, isol);

   std::vector<coot::minimol::atom *> atoms = moved_mol.select_atoms_serial();
//       std::cout << "DEBUG:: There are " << atoms.size() << " atoms from fitted zone."
// 		<< std::endl;


   // lig.make_pseudo_atoms(); uncomment for a clipper mmdb crash (sigh)

   // range_mol.write_file("range_mol.pdb");
   // mol_without_moving_zone.write_file("mol_without_moving_zone.pdb");

   // Fine.  Now we have to go back to using MMDB to interface with
   // the rest of the program.  So let's create an asc that has
   // this atom_sel.mol and the moving atoms as the
   // atom_selection. (c.f. accepting refinement or
   // regularization).

   if (atoms.size() > 0) {

      atom_selection_container_t rigid_body_asc;
      // 	 rigid_body_asc.mol = (Mymmdb::Manager *) moved_mol.pcmmdbmanager();

      // 	 int SelHnd = rigid_body_asc.mol->NewSelection();
      // 	 rigid_body_asc.mol->SelectAtoms(SelHnd, 0, "*",
      // 					 mmdb::ANY_RES, // starting resno, an int
      // 					 "*", // any insertion code
      // 					 mmdb::ANY_RES, // ending resno
      // 					 "*", // ending insertion code
      // 					 "*", // any residue name
      // 					 "*", // atom name
      // 					 "*", // elements
      // 					 "*"  // alt loc.
      // 					 );
      // 	 rigid_body_asc.mol->GetSelIndex(SelHnd,
      // 					 rigid_body_asc.atom_selection,
      // 					 rigid_body_asc.n_selected_atoms);

      success = 1;
      rigid_body_asc = make_asc(moved_mol.pcmmdbmanager(), true);

      if (false)
	 std::cout << "debug in rigid_fit() post-fit: here UDDOldAtomIndexHandle is "
		   << rigid_body_asc.UDDOldAtomIndexHandle << std::endl;

      // this seems fine.
      if (debug)
	 rigid_body_asc.mol->WritePDBASCII("post-rigid-body-refine.pdb");

      moving_atoms_asc_type = coot::NEW_COORDS_REPLACE;
      imol_moving_atoms = imol_rigid_body_refine;
      int imol = 0; // dummy (we don't need dictionary for rigid)
      make_moving_atoms_graphics_object(imol, rigid_body_asc);
      // 	 std::cout << "DEBUG:: execute_rigid_body_refine "
      // 		   << " make_moving_atoms_graphics_object UDDOldAtomIndexHandle "
      // 		   << moving_atoms_asc->UDDOldAtomIndexHandle << std::endl;
      graphics_draw();
      if (! refinement_immediate_replacement_flag) {
	 if (use_graphics_interface_flag) {
            coot::refinement_results_t dummy;
	    do_accept_reject_dialog("Rigid Body Fit", dummy); // constructed ref res
	 }
      }
      //
   } else {
      if (use_graphics_interface_flag) {
	 // GtkWidget *w = create_rigid_body_refinement_failed_dialog();
	 GtkWidget *w = widget_from_builder("rigid_body_refinement_failed_dialog");
	 gtk_widget_set_visible(w, TRUE);
      }
   }
   return success;
}

// set residue_range_atom_index_1 and residue_range_atom_index_2.
//
void
graphics_info_t::set_residue_range_refine_atoms(const std::string &chain_id,
						int resno_start, int resno_end,
						const std::string &altconf,
						int imol) {

   // do 2 atom selections to find the atom indices
   if (imol < n_molecules()) {
      if (molecules[imol].has_model()) {

	 // recall that we can't use the
	 // full_atom_spec_to_atom_index() of the class because we
	 // don't have an atom name here.

	 int ind_1 = -1, ind_2 = -1; // flags, for having found atoms

	 int SelHnd = molecules[imol].atom_sel.mol->NewSelection();
	 mmdb::PPAtom selatoms;
	 int nselatoms;

         // 	 std::cout << "DEBUG:: in set_residue_range_refine_atoms altconf :"
         // 		   << altconf << ":" << std::endl;

	 molecules[imol].atom_sel.mol->SelectAtoms(SelHnd, 0, chain_id.c_str(),
						   resno_start, // starting resno, an int
						   "*", // any insertion code
						   resno_start, // ending resno
						   "*", // ending insertion code
						   "*", // any residue name
						   "*", // atom name
						   "*", // elements
						    altconf.c_str()  // alt loc.
						   );
	 molecules[imol].atom_sel.mol->GetSelIndex(SelHnd, selatoms, nselatoms);
         // 	 std::cout << "DEBUG:: in set_residue_range_refine_atoms nselatoms (1) "
         // 		   << nselatoms << std::endl;
	 if (nselatoms > 0) {
	    if (selatoms[0]->GetUDData(molecules[imol].atom_sel.UDDAtomIndexHandle, ind_1) == mmdb::UDDATA_Ok) {
	       residue_range_atom_index_1 = ind_1;
	    }
	 }
	 molecules[imol].atom_sel.mol->DeleteSelection(SelHnd);

	 // and again for the second atom seletion:
	 //
	 SelHnd = molecules[imol].atom_sel.mol->NewSelection();
	 molecules[imol].atom_sel.mol->SelectAtoms(SelHnd, 0, chain_id.c_str(),
						   resno_end, // starting resno, an int
						   "*", // any insertion code
						   resno_end, // ending resno
						   "*", // ending insertion code
						   "*", // any residue name
						   "*", // atom name
						   "*", // elements
						   altconf.c_str()  // alt loc.
						   );
	 molecules[imol].atom_sel.mol->GetSelIndex(SelHnd, selatoms, nselatoms);
// 	 std::cout << "DEBUG:: in set_residue_range_refine_atoms nselatoms (2) "
// 		   << nselatoms << std::endl;
	 if (nselatoms > 0) {
	    if (selatoms[0]->GetUDData(molecules[imol].atom_sel.UDDAtomIndexHandle, ind_2) == mmdb::UDDATA_Ok) {
	       residue_range_atom_index_2 = ind_2;
	    }
	 }
	 molecules[imol].atom_sel.mol->DeleteSelection(SelHnd);

	 //   if (ind_1 >= 0 && ind_2 >= 0)
	 //     execute_rigid_body_refine(0); // not autorange
      }
   }
}


// The passed residue type is either N, C or (now [20031222]) M.
//
int
graphics_info_t::execute_add_terminal_residue(int imol,
					      const std::string &terminus_type,
					      mmdb::Residue *res_p,
					      const std::string &chain_id,
					      const std::string &res_type_in,
					      bool immediate_addition_flag) {

   int state = 0;

   // Calling function also does a check for a map, I think.

   std::string res_type = res_type_in; // const
   int imol_map = Imol_Refinement_Map();
   if (imol_map == -1) {
      // just shove it on without a map
      if (molecules[imol].has_model()) {
	 // float phi = graphics_info_t::terminal_residue_addition_direct_phi;
	 // float psi = graphics_info_t::terminal_residue_addition_direct_psi;
	 mmdb::Manager *orig_mol = graphics_info_t::molecules[imol].atom_sel.mol;
	 //	    mmdb::Residue *res_new = add_terminal_residue_directly(terminus_type, res_p,
	 // chain_id, res_type, phi, psi);
	 mmdb::Residue *res_new = 0;
	 mmdb::Manager *new_mol = coot::util::create_mmdbmanager_from_residue(res_new);
	 if (new_mol) {
	    atom_selection_container_t extra_residue_asc = make_asc(new_mol);
	    graphics_info_t::molecules[imol].add_coords(extra_residue_asc);
	 }
      }
   } else {

      int resno_added = -1; // was unset

      if (terminus_type == "not-terminal-residue") {
	      std::string s = "That residue was not at a terminus";
	      std::cout << s << std::endl;
	      add_status_bar_text(s);
      } else {

	 imol_moving_atoms = imol;

	 mmdb::Residue *upstream_neighbour_residue_p = 0;
	 mmdb::Residue *downstream_neighbour_residue_p = 0;
	 std::string residue_type_string = res_type;
	 int residue_number = res_p->GetSeqNum();  // bleugh.
	 if (residue_type_string == "auto") {
	    if (terminus_type == "C" || terminus_type == "MC")
	       resno_added = residue_number + 1;
	    if (terminus_type == "N" || terminus_type == "MN")
	       resno_added = residue_number - 1;
	    std::pair<bool, std::string> p =
	       molecules[imol].find_terminal_residue_type(chain_id, resno_added,
							  alignment_wgap,
							  alignment_wspace);
	    if (p.first) {
	       res_type = p.second;
	    } else {
	       res_type = "ALA";
	    }
	 }

	 if (terminus_type == "C") {
	    // we have upstream residue
	    upstream_neighbour_residue_p = coot::util::previous_residue(res_p);
	 }
	 if (terminus_type == "N") {
	    // we have downstream residue
	    downstream_neighbour_residue_p = coot::util::next_residue(res_p);
	 }


	 float bf = default_new_atoms_b_factor;
	 coot::residue_by_phi_psi addres(terminus_type, res_p, chain_id, res_type, bf);

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
	 // 	 short int mask_waters_flag = 1;
	 //      addres.mask_map(molecules[imol].atom_sel.mol, mask_waters_flag);
	 //
	 mmdb::PPAtom atom_sel = NULL;
	 int n_selected_atoms = 0;
	 mmdb::realtype radius = 8.0;  // more than enough for 2 residue mainchains.
	 int SelHndSphere = molecules[imol].atom_sel.mol->NewSelection();
	 mmdb::Atom *terminal_at = NULL;
	 std::string atom_name = "Unassigned";
	 if (terminus_type == "MC" || terminus_type == "C" || terminus_type == "singleton")
	    atom_name = " C  ";
	 if (terminus_type == "MN" || terminus_type == "N")
	    atom_name = " N  ";
	 if (atom_name != "Unassigned") {
	    mmdb::PPAtom residue_atoms;
	    int nResidueAtoms;
	    mmdb::Residue *res_tmp_p = (mmdb::Residue *) res_p;
	    res_tmp_p->GetAtomTable(residue_atoms, nResidueAtoms);
	    for (int i=0; i<nResidueAtoms; i++)
	       if (atom_name == residue_atoms[i]->name) {
		  terminal_at = residue_atoms[i];
		  break;
	       }

	    if (terminal_at) {
	       molecules[imol].atom_sel.mol->SelectSphere(SelHndSphere, mmdb::STYPE_ATOM,
							  terminal_at->x,
							  terminal_at->y,
							  terminal_at->z,
							  radius, mmdb::SKEY_NEW);
	       molecules[imol].atom_sel.mol->GetSelIndex(SelHndSphere, atom_sel, n_selected_atoms);
	       int invert_flag = 0;

	       // Why is this commented out?
	       // Because the map is added to addres after this point. Hmm. Masking is a good idea.
	       // Maybe copy the molecules[imol_map].xmap and mask it making xmap_masked
	       // and pass xmap_masked to the best_fit_phi_psi()?
	       //
	       // addres.mask_map(molecules[imol].atom_sel.mol, SelHndSphere, invert_flag);

	       molecules[imol].atom_sel.mol->DeleteSelection(SelHndSphere);
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
				    molecules[imol_map].xmap);

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
	       if (use_graphics_interface_flag) {
		  // GtkWidget *w = create_add_terminal_residue_finds_none_dialog();
		  GtkWidget *w = widget_from_builder("add_terminal_residue_finds_none_dialog");
		  gtk_widget_set_visible(w, TRUE);
	       }

	    } else {

	       state = 1;

	       atom_selection_container_t terminal_res_asc;

	       // if this is begin added to a shelx molecule, then we
	       // need to set the occs to 11.0
	       //
	       if (graphics_info_t::molecules[imol].is_from_shelx_ins()) {
		  bf = 11.0;
	       }

	       if (add_terminal_residue_add_other_residue_flag) {

		  // check that the other residue is not in the molecule before adding
		  // all of mmol. If it is already in the molecule, remove it from mmol
		  if (terminus_type == "C" || terminus_type == "MC") {
		     coot::residue_spec_t other_residue_spec(chain_id, resno_added+1, "");
		     mmdb::Residue *res_other = graphics_info_t::molecules[imol].get_residue(other_residue_spec);
		     if (res_other) {
			mmol[0][other_residue_spec.res_no].atoms.clear();
		     }
		  } else {
		     if (terminus_type == "N" || terminus_type == "MN") {
			coot::residue_spec_t other_residue_spec(chain_id, resno_added-1, "");
			mmdb::Residue *res_other = graphics_info_t::molecules[imol].get_residue(other_residue_spec);
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
 		  add_side_chain_to_terminal_res(terminal_res_asc, res_type, terminus_type,
						 add_terminal_residue_add_other_residue_flag);


 	       // std::cout << "-------------- tmp_asc --------" << std::endl;
 	       // debug_atom_selection_container(tmp_asc);

	       // If this is wrong also consider fixing execute_rigid_body_refine()
	       //
// 	       std::cout << "debug: add_residue asc has n_selected_atoms = "
// 			 << terminal_res_asc.n_selected_atoms << " "
// 			 << terminal_res_asc.atom_selection << std::endl;

// 	       for (int i=0; i< terminal_res_asc.n_selected_atoms; i++) {
// 		  std::cout << "debug: add_residue asc has chain_id: "
// 			    << terminal_res_asc.atom_selection[i]->GetChainID()
// 			    << " for " << terminal_res_asc.atom_selection[i]
// 			    << std::endl;
// 	       }

	       coot::residue_spec_t rs(res_p);
	       graphics_info_t::molecules[imol].remove_ter_atoms(rs);

	       if (! immediate_addition_flag) {
		  make_moving_atoms_graphics_object(imol, tmp_asc);
		  moving_atoms_asc_type = coot::NEW_COORDS_INSERT;
		  graphics_draw();
		  coot::refinement_results_t dummy;
		  if (use_graphics_interface_flag) {
		     do_accept_reject_dialog("Terminal Residue", dummy);
		  }
	       } else {

		  if (molecules[imol].is_from_shelx_ins()) {
		     for (int i=0; i<tmp_asc.n_selected_atoms; i++) {
			tmp_asc.atom_selection[i]->occupancy = 11.0;
		     }
		  }

		  molecules[imol_moving_atoms].insert_coords(tmp_asc);

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
		     clipper::Coord_orth new_o_pos =
			addres.best_fit_phi_psi_attaching_oxygen_position_update(mmol, res_p);
		     molecules[imol_moving_atoms].move_atom(" O  ", res_p, new_o_pos);
		  }
                  // method from master:
		  // if (terminus_type == "C" || terminus_type == "MC")
		  //    molecules[imol_moving_atoms].move_O_atom_of_added_to_residue(res_p, chain_id);

                  graphics_draw();
               }
            }
         }
      }
   }
   return state;
}

void
graphics_info_t::execute_simple_nucleotide_addition(int imol, const std::string &chain_id, int res_no) {

   if (! is_valid_model_molecule(imol)) {
      std::cout << "WARNING:: wrong model " << imol << std::endl;
      return;
   }

   mmdb::Residue *residue_p = molecules[imol].get_residue(chain_id, res_no, "");
   if (! residue_p) {
      std::cout << "WARNING:: missing-residue" << chain_id << " " << res_no << std::endl;
   } else {
      std::string term_type = "";
      mmdb::Residue *r_p = molecules[imol].get_residue(chain_id, res_no-1, "");
      mmdb::Residue *r_n = molecules[imol].get_residue(chain_id, res_no+1, "");
      if (r_p  && ! r_n) term_type = "C";
      if (r_n  && ! r_p) term_type = "N";
      if (!r_n && ! r_p) term_type = "MC";
      execute_simple_nucleotide_addition(imol, term_type, residue_p, chain_id);
   }
}


void
graphics_info_t::execute_simple_nucleotide_addition(int imol, const std::string &term_type, 
                                                    mmdb::Residue *res_p, const std::string &chain_id) {


   // If it's RNA beam it in in ideal A form,
   // If it's DNA beam it in in ideal B form

   // What's the plan?
   //
   // OK the plan is to generate a 2 residue molecule of
   // single-stranded RNA (or DNA).
   //
   // Depending on if this is N or C terminal type, we define the
   // sequence, adding a "base" residue (that we'll use to match to
   // res_p);

   if (term_type == "not-terminal-residue") {
      std::cout << "That was not a terminal residue (check for neighbour solvent residues maybe) "
                << coot::residue_spec_t(res_p) << std::endl;
      add_status_bar_text("That was not a terminal residue.");
   } else {

      std::string seq = "aa";
      std::string RNA_or_DNA_str = "RNA";
      std::string form_str = "A";
      short int single_stranded_flag = 1;


      if (is_valid_model_molecule(imol)) {
	 int residue_number = res_p->GetSeqNum();
	 int resno_added = -1; // was unset
	 if (term_type == "C" || term_type == "MC")
	    resno_added = residue_number + 1;
	 if (term_type == "N" || term_type == "MN")
	    resno_added = residue_number - 1;
	 if (resno_added != -1) {
	    bool is_nucleic_acid = true;
	    std::pair<bool, std::string> p =
	       molecules[imol].find_terminal_residue_type(chain_id, resno_added,
							  alignment_wgap,
							  alignment_wspace, is_nucleic_acid);
	    if (p.first) {
	       seq = "a" + coot::util::downcase(p.second);
	    }
	 }
      }



      if (coot::util::nucleotide_is_DNA(res_p)) {
	 RNA_or_DNA_str = "DNA";
	 form_str = "B";
      }

      coot::ideal_rna ir(RNA_or_DNA_str, form_str, single_stranded_flag,
			 seq, graphics_info_t::standard_residues_asc.mol);
      ir.use_v3_names();
      mmdb::Manager *mol = ir.make_molecule();

      int match_resno;
      int interesting_resno;
      if (term_type == "C" || term_type == "MC") {
	 match_resno = 1;
	 interesting_resno = 2;
      } else {
	 interesting_resno = 1;
	 match_resno = 2;
      }

      mmdb::Residue *moving_residue_p = NULL;
      mmdb::Residue *interesting_residue_p = NULL;
      int imod = 1;
      // now set moving_residue_p and interesting_residue_p:
      mmdb::Model *model_p = mol->GetModel(imod);
      mmdb::Chain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 mmdb::PResidue residue_p;
	 for (int ires=0; ires<nres; ires++) {
	    residue_p = chain_p->GetResidue(ires);
	    // 	 std::cout << "testing vs resno " << residue_p->GetSeqNum()
	    // 		   << std::endl;
	    if (residue_p->GetSeqNum() == match_resno) {
	       moving_residue_p = residue_p;
	    }
	    if (residue_p->GetSeqNum() == interesting_resno) {
	       interesting_residue_p = residue_p;
	    }
	    if (moving_residue_p && interesting_residue_p)
	       break;
	 }
	 if (moving_residue_p && interesting_residue_p)
	    break;
      }

      if (interesting_residue_p) {
	 if (moving_residue_p) {
	    bool use_old_names = convert_to_v2_atom_names_flag;

	    std::pair<bool, clipper::RTop_orth> rtop_pair =
	       coot::util::nucleotide_to_nucleotide(res_p, moving_residue_p,
						    use_old_names);

	    // now apply rtop to mol:
	    if (rtop_pair.first) {
	       // fix up the residue number and chain id to match the clicked atom
	       int new_resno = res_p->GetSeqNum() + interesting_resno - match_resno;
	       interesting_residue_p->seqNum = new_resno;

               // we always want to remove OP3 from the residue to which a new residue
               // is added when we add to the "N-terminus"
               //
               if (term_type == "N" || term_type == "MN") {
                  mmdb::Atom **residue_atoms = 0;
                  int n_residue_atoms = 0;
                  bool deleted = false;
                  res_p->GetAtomTable(residue_atoms, n_residue_atoms);
                  for (int iat=0; iat<n_residue_atoms; iat++) {
                     mmdb::Atom *at = residue_atoms[iat];
                     if (at) {
                        std::string at_name(at->name);
                        if (at_name == " OP3") {  // PDBv3 FIXME
                           delete at;
                           at = NULL;
                           deleted = true;
                           break;
                        }
                     }
                  }
                  if (deleted)
                     res_p->TrimAtomTable();
               }

	       coot::util::transform_mol(mol, rtop_pair.second);
	       // byte gz = GZM_NONE;
	       // mol->WritePDBASCII("overlapped.pdb", gz);
	       mmdb::Manager *residue_mol =
		  coot::util::create_mmdbmanager_from_residue(interesting_residue_p);

	       atom_selection_container_t asc = make_asc(residue_mol);
	       // set the chain id of the chain that contains interesting_residue_p:
	       model_p = residue_mol->GetModel(imod);
	       // run over chains of the existing mol
	       nchains = model_p->GetNumberOfChains();
	       for (int ichain=0; ichain<nchains; ichain++) {
		  chain_p = model_p->GetChain(ichain);
		  int nres = chain_p->GetNumberOfResidues();
		  mmdb::PResidue residue_p;
		  for (int ires=0; ires<nres; ires++) {
		     residue_p = chain_p->GetResidue(ires);
		     if (residue_p->GetSeqNum() == interesting_residue_p->GetSeqNum()) {
			chain_p->SetChainID(chain_id.c_str());
		     }
		  }
	       }
	       graphics_info_t::molecules[imol].insert_coords(asc);

	       if (add_terminal_residue_do_post_refine) {
		  // shall we refine it?  If there is a map, yes.
		  int imol_map = Imol_Refinement_Map();
		  if (imol_map >= 0) {
		     refine_residue_range(imol, chain_id, chain_id, new_resno, "",
					  new_resno, "", "", 0);
		  }
	       }
	    }
	 }
      } else {
	 std::cout << "Failed to find interesting residue (with resno " << interesting_resno
		   << ")" << std::endl;
      }
      delete mol;
      graphics_draw();
   }
}


void
graphics_info_t::execute_rotate_translate_ready() { // manual movement

   // std::cout << "execute_rotate_translate_ready() --- start ---" << std::endl;

   // now we are called by chain and molecule pick (as well as the old
   // zone pick).

   // We use the rot_trans_object_type to distinguish.

   const char *chain_id = "*";
   int resno_1 = mmdb::ANY_RES;
   int resno_2 = mmdb::ANY_RES;
   std::string insertion_code_selection = "*"; // reset on start and stop residue in range being the same
   bool good_settings = 0; // fail initially
   mmdb::Atom *atom1 = molecules[imol_rot_trans_object].atom_sel.atom_selection[rot_trans_atom_index_1];
   const char *altLoc = atom1->altLoc;
   // This uses moving_atoms_asc internally, we don't need to pass it:
   coot::atom_spec_t origin_atom_spec(atom1);

//    std::cout << "debug:: rot_trans_atom_index_1 " << rot_trans_atom_index_1 << " gives atom "
// 	     << atom1 << std::endl;

   if (rot_trans_object_type == ROT_TRANS_TYPE_CHAIN) {
      chain_id = atom1->GetChainID();
      altLoc = "*";
      good_settings = 1;
   }

   if (rot_trans_object_type == ROT_TRANS_TYPE_MOLECULE) {
      // use default settings
      altLoc = "*";
      good_settings = 1;
   }

   // with ROT_TRANS_TYPE_RESIDUE, the atom_2 is valid and set to be the same as atom_1
   //
   if (rot_trans_object_type == ROT_TRANS_TYPE_ZONE || rot_trans_object_type == ROT_TRANS_TYPE_RESIDUE) {
      
      std::cout << "execute_rotate_translate_ready() idx_1 " << rot_trans_atom_index_1 << std::endl;
      std::cout << "execute_rotate_translate_ready() idx_2 " << rot_trans_atom_index_2 << std::endl;
      
      mmdb::Atom *atom2 = molecules[imol_rot_trans_object].atom_sel.atom_selection[rot_trans_atom_index_2];

      std::cout << "execute_rotate_translate_ready atom_1 " << atom1 << std::endl;
      std::cout << "execute_rotate_translate_ready atom_2 " << atom2 << std::endl;

      char *chain_id_1 = atom1->GetChainID();
      char *chain_id_2 = atom2->GetChainID();

      if (chain_id_1 != chain_id_2) {
	 std::string info_string("Atoms must be in the same chain");
	 add_status_bar_text(info_string);
      } else {
	 chain_id = chain_id_1;
	 resno_1 = atom1->GetSeqNum();
	 resno_2 = atom2->GetSeqNum();
	 if (resno_1 > resno_2) {
	    int tmp = resno_1;
	    resno_1 = resno_2;
	    resno_2 = tmp;
	 }
	 if (atom1->residue == atom2->residue) {
	    insertion_code_selection = atom1->GetInsCode();
	 }

	 origin_atom_spec = coot::atom_spec_t(atom2); // as it used to be

	 good_settings = 1;
      }
   }

   if (good_settings) {

      // GtkWidget *widget = create_rotate_translate_obj_dialog();
      GtkWidget *widget = widget_from_builder("rotate_translate_obj_dialog");
      GtkWindow *main_window = GTK_WINDOW(get_main_window());
      gtk_window_set_transient_for(GTK_WINDOW(widget), main_window);

      do_rot_trans_adjustments(widget);

      // set its position if it was shown before
      if (rotate_translate_x_position > -100) {
         /*
	 gtk_widget_set_uposition(widget,
				  rotate_translate_x_position,
				  rotate_translate_y_position);
         */
      }
      gtk_widget_set_visible(widget, TRUE);

      atom_selection_container_t rt_asc;
      // No! It cannot point to the same mmdb::Atoms.
      // rt_asc.mol = molecules[imol_rot_trans_object].atom_sel.mol;
      // Mymmdb::Manager *mol = new Mymmdb::Manager;
      // mol->Copy(molecules[imol_rot_trans_object].atom_sel.mol, mmdb::MMDBFCM_All);
      // how about we instead use:
      // mmdb::Manager *mol = create_mmdbmanager_from_res_selection();
      //
      mmdb::PResidue *sel_residues = NULL;
      int n_sel_residues;
      int selHnd = molecules[imol_rot_trans_object].atom_sel.mol->NewSelection();
      molecules[imol_rot_trans_object].atom_sel.mol->Select(selHnd, mmdb::STYPE_RESIDUE, 0,
							    chain_id,
							    resno_1, insertion_code_selection.c_str(),
							    resno_2, insertion_code_selection.c_str(),
							    "*",  // residue name
							    "*",  // Residue must contain this atom name?
							    "*",  // Residue must contain this Element?
							    "*",  // altLocs
							    mmdb::SKEY_NEW // selection key
							    );
      molecules[imol_rot_trans_object].atom_sel.mol->GetSelIndex(selHnd, sel_residues, n_sel_residues);

      short int alt_conf_split_flag = 0;
      std::string altloc_string(altLoc);
      if (altloc_string != "")
	 alt_conf_split_flag = 1;

      // create a complete new clean copy of chains/residues/atoms
      std::pair<mmdb::Manager *, int> mp(0, 0);


      if (rot_trans_object_type == ROT_TRANS_TYPE_ZONE || rot_trans_object_type == ROT_TRANS_TYPE_RESIDUE)
	mp =
	 coot::util::create_mmdbmanager_from_res_selection(molecules[imol_rot_trans_object].atom_sel.mol,
							   sel_residues, n_sel_residues,
							   0, 0, altloc_string, chain_id,
							   alt_conf_split_flag);


      if (rot_trans_object_type == ROT_TRANS_TYPE_CHAIN)
	mp =
	 coot::util::create_mmdbmanager_from_res_selection(molecules[imol_rot_trans_object].atom_sel.mol,
							   sel_residues, n_sel_residues,
							   0, 0, altloc_string, chain_id,
							   alt_conf_split_flag);

      if (rot_trans_object_type == ROT_TRANS_TYPE_MOLECULE)
	mp =
	  coot::util::create_mmdbmanager_from_mmdbmanager(molecules[imol_rot_trans_object].atom_sel.mol);

      rt_asc = make_asc(mp.first);
      rt_asc.UDDOldAtomIndexHandle = mp.second;

      molecules[imol_rot_trans_object].atom_sel.mol->DeleteSelection(selHnd);

      //    std::cout << "DEBUG:: rt_asc: has n_selected_atoms: " << rt_asc.n_selected_atoms
      // 	     << std::endl;
      imol_moving_atoms = imol_rot_trans_object;
      moving_atoms_asc_type = coot::NEW_COORDS_REPLACE;
      make_moving_atoms_graphics_object(imol_rot_trans_object, rt_asc, false, false); // shallow copy rt_asc to moving_atoms_asc

      if (false) {
         std::cout << "in execute_rotate_translate_ready() with moving_atoms_asc " << moving_atoms_asc << std::endl;
         std::cout << "in execute_rotate_translate_ready() with moving_atoms_asc mol " << moving_atoms_asc->mol << std::endl;
      }

      // set the rotation centre atom index:
      //   rot_trans_atom_index_rotation_origin_atom =
      //       find_atom_index_in_moving_atoms(chain_id,
      // 				      atom2->GetSeqNum(),
      // 				      atom2->name);  // uses moving_atoms_asc


      rot_trans_rotation_origin_atom = find_atom_in_moving_atoms(origin_atom_spec);

      if (false) {
	 if (rot_trans_rotation_origin_atom) {
	    std::cout << "DEBUG:: atom spec in moving atom " << origin_atom_spec << " returns "
		      << rot_trans_rotation_origin_atom << std::endl;
	 } else {
	    std::cout << "DEBUG:: atom spec in moving atom " << origin_atom_spec << " returns NULL "
		      << std::endl;
	 }
      }

      //    std::cout << "DEBUG:: in execute_rotate_translate_read, found rotation atom: "
      // 	     << rot_trans_rotation_origin_atom << std::endl;

      if (rot_trans_rotation_origin_atom == NULL) {
	 std::cout << "ERROR:: rot_trans_atom_rotation_origin not found" << std::endl;
      }
      graphics_draw();

      std::string info_string("Drag on an atom to translate residue, Ctrl Drag off atoms to rotate residue");
      add_status_bar_text(info_string);
   }

}


void
graphics_info_t::execute_torsion_general() {

   if (torsion_general_atom_index_1_mol_no == torsion_general_atom_index_2_mol_no) {
      if (torsion_general_atom_index_1_mol_no == torsion_general_atom_index_3_mol_no) {
	 if (torsion_general_atom_index_1_mol_no == torsion_general_atom_index_4_mol_no) {
	    if (torsion_general_atom_index_4_mol_no < n_molecules()) {

	       mmdb::Atom *atom_1 = 0;
	       mmdb::Atom *atom_2 = 0;
	       mmdb::Atom *atom_3 = 0;
	       mmdb::Atom *atom_4 = 0;
	       int im = torsion_general_atom_index_1_mol_no;

	       if (torsion_general_atom_index_1 < molecules[im].atom_sel.n_selected_atoms) {
		  if (torsion_general_atom_index_2 < molecules[im].atom_sel.n_selected_atoms) {
		     if (torsion_general_atom_index_3 < molecules[im].atom_sel.n_selected_atoms) {
			if (torsion_general_atom_index_4 < molecules[im].atom_sel.n_selected_atoms) {

			   atom_1 = molecules[im].atom_sel.atom_selection[torsion_general_atom_index_1];
			   atom_2 = molecules[im].atom_sel.atom_selection[torsion_general_atom_index_2];
			   atom_3 = molecules[im].atom_sel.atom_selection[torsion_general_atom_index_3];
			   atom_4 = molecules[im].atom_sel.atom_selection[torsion_general_atom_index_4];

			   mmdb::Residue *r1 = atom_1->GetResidue();
			   mmdb::Residue *r2 = atom_2->GetResidue();
			   mmdb::Residue *r3 = atom_3->GetResidue();
			   mmdb::Residue *r4 = atom_4->GetResidue();

			   // pointer comparison:
			   if (r1 == r2) {
			      if (r1 == r3) {
				 if (r1 == r4) {

				    moving_atoms_asc_type = coot::NEW_COORDS_REPLACE;
				    in_edit_torsion_general_flag = 1;
				    imol_moving_atoms = im;
				    int ai  = torsion_general_atom_index_1;
				    short int whole_res_flag = 0;
				    atom_selection_container_t residue_asc =
				       graphics_info_t::molecules[im].edit_residue_pull_residue(ai, whole_res_flag);
				    regularize_object_bonds_box.clear_up();
				    make_moving_atoms_graphics_object(im, residue_asc);

				    std::vector<coot::atom_spec_t> as;
				    as.push_back(coot::atom_spec_t(atom_1));
				    as.push_back(coot::atom_spec_t(atom_2));
				    as.push_back(coot::atom_spec_t(atom_3));
				    as.push_back(coot::atom_spec_t(atom_4));
				    torsion_general_atom_specs = as;
				    graphics_draw();
				    torsion_general_reverse_flag = 0;
				    mmdb::Residue *res_local = get_first_res_of_moving_atoms();
				    if (res_local) {

				       // save them for later usage (when the mouse is moved)
				       coot::contact_info contact = coot::getcontacts(*moving_atoms_asc);
				       // contact.print(); // debug
				       torsion_general_contact_indices = contact.get_contact_indices();
				       chi_angle_alt_conf = atom_4->altLoc;

				       coot::refinement_results_t dummy;
				       if (use_graphics_interface_flag)
					  do_accept_reject_dialog("Torsion General", dummy);
				    }
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
}

mmdb::Residue *
graphics_info_t::get_first_res_of_moving_atoms() {

   mmdb::Residue *r = 0;
   mmdb::Model *model_p = moving_atoms_asc->mol->GetModel(1);
   if (model_p) {
      mmdb::Chain *chain_p = model_p->GetChain(0);
      if (chain_p) {
	 mmdb::Residue *residue_p = chain_p->GetResidue(0);
	 if (residue_p) {
	    r = residue_p;
	 }
      }
   }
   return r;
}


void
graphics_info_t::do_rot_trans_adjustments(GtkWidget *dialog) {

   std::vector<std::string> hscale_lab;

   hscale_lab.push_back("rotate_translate_obj_xtrans_hscale");
   hscale_lab.push_back("rotate_translate_obj_ytrans_hscale");
   hscale_lab.push_back("rotate_translate_obj_ztrans_hscale");
   hscale_lab.push_back("rotate_translate_obj_xrot_hscale");
   hscale_lab.push_back("rotate_translate_obj_yrot_hscale");
   hscale_lab.push_back("rotate_translate_obj_zrot_hscale");

// GtkObject *gtk_adjustment_new( gfloat value,
//                                gfloat lower,
//                                gfloat upper,
//                                gfloat step_increment,
//                                gfloat page_increment,
//                                gfloat page_size );

   for (unsigned int i=0; i<hscale_lab.size(); i++) {
      // GtkWidget *hscale = lookup_widget(dialog, hscale_lab[i].c_str());
      GtkWidget *hscale = widget_from_builder(hscale_lab[i]); // 20220311-PE is this right?
      GtkAdjustment *adj = GTK_ADJUSTMENT(gtk_adjustment_new(0.0, -180.0, 360.0, 0.1, 1.0, 0));
      gtk_range_set_adjustment(GTK_RANGE(hscale), GTK_ADJUSTMENT(adj));
      g_signal_connect(G_OBJECT(adj),
			 "value_changed",
			 G_CALLBACK(graphics_info_t::rot_trans_adjustment_changed),
			 GINT_TO_POINTER(i));
   }
}


coot::ScreenVectors::ScreenVectors() {

   graphics_info_t g;
   glm::vec4 glm_centre = g.unproject(0, 0, 0.5);
   glm::vec4 glm_front  = g.unproject(0, 0, 0.0);
   glm::vec4 glm_right  = g.unproject(1, 0, 0.5);
   glm::vec4 glm_top    = g.unproject(0, 1, 0.5);

   coot::Cartesian centre(glm_centre.x, glm_centre.y, glm_centre.z);
   coot::Cartesian front(glm_front.x, glm_front.y, glm_front.z);
   coot::Cartesian right(glm_right.x, glm_right.y, glm_right.z);
   coot::Cartesian top(glm_top.x, glm_top.y, glm_top.z);

   screen_x = (right - centre);
   screen_y = (top   - centre);
   screen_z = (front - centre);

   screen_x.unit_vector_yourself();
   screen_y.unit_vector_yourself();
   screen_z.unit_vector_yourself();

}


// static
void
graphics_info_t::rot_trans_adjustment_changed(GtkAdjustment *adj, gpointer user_data) {

   // std::cout << "rot_trans_adjustment_changed() ... start ... " << std::endl;

   graphics_info_t g;  // because rotate_round_vector is not static - it should be.
                       // FIXME at some stage.
   double v = gtk_adjustment_get_value(adj);

   int i_hscale = GPOINTER_TO_INT(user_data);
   short int do_rotation;
   if (i_hscale < 3)
      do_rotation = 0;
   else
      do_rotation = 1;

   // std::cout << "rot_trans_adjustment_changed: user_data: " << i_hscale << std::endl;
   double x_diff = v;
   if ( previous_rot_trans_adjustment[i_hscale] > -9999) {
      x_diff = v - previous_rot_trans_adjustment[i_hscale];
   }
   previous_rot_trans_adjustment[i_hscale] = v;

   // std::cout << "using  " << x_diff << "  " << v << "  "
   //      << previous_rot_trans_adjustment[i_hscale] << std::endl;

   coot::ScreenVectors screen_vectors;

   float x_add = 0.0;
   float y_add = 0.0;
   float z_add = 0.0;

   if (i_hscale == 0) {
      x_add = screen_vectors.screen_x.x() * x_diff * 0.002 * zoom;
      y_add = screen_vectors.screen_x.y() * x_diff * 0.002 * zoom;
      z_add = screen_vectors.screen_x.z() * x_diff * 0.002 * zoom;
   }
   if (i_hscale == 1) {
      x_add = screen_vectors.screen_y.x() * x_diff * -0.002 * zoom;
      y_add = screen_vectors.screen_y.y() * x_diff * -0.002 * zoom;
      z_add = screen_vectors.screen_y.z() * x_diff * -0.002 * zoom;
   }
   if (i_hscale == 2) {
      x_add = screen_vectors.screen_z.x() * x_diff * 0.002 * zoom;
      y_add = screen_vectors.screen_z.y() * x_diff * 0.002 * zoom;
      z_add = screen_vectors.screen_z.z() * x_diff * 0.002 * zoom;
   }
   // std::cout << "Here 1 with x_add, y_add z_add " << x_add << " " << y_add << " " << z_add << std::endl;

   if (do_rotation) {

      clipper::Coord_orth screen_vector; // the vector to rotate about
      if (i_hscale == 3) {
	 do_rotation = 1;
	 screen_vector = clipper::Coord_orth(screen_vectors.screen_x.x(),
					     screen_vectors.screen_x.y(),
					     screen_vectors.screen_x.z());
      }
      if (i_hscale == 4) {
	 do_rotation = 1;
	 screen_vector = clipper::Coord_orth(screen_vectors.screen_y.x(),
					     screen_vectors.screen_y.y(),
					     screen_vectors.screen_y.z());
      }
      if (i_hscale == 5) {
	 do_rotation = 1;
	 screen_vector = clipper::Coord_orth(screen_vectors.screen_z.x(),
					     screen_vectors.screen_z.y(),
					     screen_vectors.screen_z.z());
      }


      // int indx = rot_trans_atom_index_rotation_origin_atom;
      mmdb::Atom *rot_centre = rot_trans_rotation_origin_atom;
      clipper::Coord_orth rotation_centre(0,0,0); // updated.

      // But! maybe we have a different rotation centre
      if (rot_trans_zone_rotates_about_zone_centre) {
	 if (moving_atoms_asc->n_selected_atoms  > 0) {
	    rotation_centre = g.moving_atoms_centre();
	 }
      } else {
	 if (rot_centre) {
	   rotation_centre = clipper::Coord_orth(rot_centre->x,
						 rot_centre->y,
						 rot_centre->z);
	 } else {
	   std::cout << "WARNING:: rot_centre atom not found" << std::endl;
	 }
      }


      for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
	 clipper::Coord_orth co(moving_atoms_asc->atom_selection[i]->x,
				moving_atoms_asc->atom_selection[i]->y,
				moving_atoms_asc->atom_selection[i]->z);

         // std::cout << "rotating atom by " << x_diff * 0.018 << std::endl;

	 clipper::Coord_orth new_pos =
	    coot::util::rotate_around_vector(screen_vector, co, rotation_centre, x_diff * 0.018);
	 moving_atoms_asc->atom_selection[i]->x = new_pos.x();
	 moving_atoms_asc->atom_selection[i]->y = new_pos.y();
	 moving_atoms_asc->atom_selection[i]->z = new_pos.z();
      }

   } else {

      for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
         // std::cout << "translating atom " << i << " by " << x_add << " " << y_add << " " << z_add << std::endl;
	 moving_atoms_asc->atom_selection[i]->x += x_add;
	 moving_atoms_asc->atom_selection[i]->y += y_add;
	 moving_atoms_asc->atom_selection[i]->z += z_add;
      }
   }

   int originating_molecule_bonds_box_type = molecules[imol_moving_atoms].Bonds_box_type();

   if (originating_molecule_bonds_box_type == coot::CA_BONDS ||
       originating_molecule_bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS ||
       originating_molecule_bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS_AND_SIDECHAINS ||
       originating_molecule_bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR ||
       originating_molecule_bonds_box_type == coot::COLOUR_BY_RAINBOW_BONDS) {

      Bond_lines_container bonds;
      bool draw_hydrogens_flag = false;
      if (molecules[imol_moving_atoms].draw_hydrogens())
	 draw_hydrogens_flag = true;

      bonds.do_Ca_plus_ligands_bonds(*moving_atoms_asc, imol_moving_atoms, Geom_p(), 1.0, 4.7,
                                     draw_missing_loops_flag, draw_hydrogens_flag);

      moving_atoms_molecule.atom_sel = *moving_atoms_asc;
      moving_atoms_molecule.bonds_box = bonds.make_graphical_bonds();
      attach_buffers();
      moving_atoms_molecule.make_meshes_from_bonds_box_instanced_version();

   } else {

      int do_disulphide_flag = 0;
      int do_bonds_to_hydrogens = 1;
      bool do_rama_markup = false;
      bool do_rota_markup = false;
      Bond_lines_container bonds(*moving_atoms_asc, imol_moving_atoms, do_disulphide_flag, do_bonds_to_hydrogens,
                                 do_rama_markup, do_rota_markup);

      moving_atoms_molecule.atom_sel = *moving_atoms_asc;
      moving_atoms_molecule.bonds_box = bonds.make_graphical_bonds();
      attach_buffers();
      moving_atoms_molecule.make_meshes_from_bonds_box_instanced_version();
   }

   graphics_draw();
}


// --- nudge active residue
// static
void
graphics_info_t::nudge_active_residue(guint direction) {

   std::cout << "nudge_active_residue() " << std::endl;

   std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = graphics_info_t::active_atom_spec();
   if (active_atom.first) {
      clipper::Coord_orth shift(0,0,0);
      clipper::Mat33<double> mat(1,0,0,0,1,0,0,0,1);
      double shift_scale_factor = 0.02 * zoom; // needs to be 0.04 for funny mode?
      coot::ScreenVectors screen_vectors;

      if (direction == GDK_KEY_Left) {
	 // std::cout << "Left nudge residue" << std::endl;
	 shift = clipper::Coord_orth(-shift_scale_factor * screen_vectors.screen_x.x(),
				     -shift_scale_factor * screen_vectors.screen_x.y(),
				     -shift_scale_factor * screen_vectors.screen_x.z());
      }
      if (direction == GDK_KEY_Right) {
	 // std::cout << "Right nudge residue" << std::endl;
	 shift = clipper::Coord_orth(shift_scale_factor * screen_vectors.screen_x.x(),
				     shift_scale_factor * screen_vectors.screen_x.y(),
				     shift_scale_factor * screen_vectors.screen_x.z());
      }
      if (direction == GDK_KEY_Up) {
	 // std::cout << "Up nudge residue" << std::endl;
	 shift = clipper::Coord_orth(-shift_scale_factor * screen_vectors.screen_y.x(),
				     -shift_scale_factor * screen_vectors.screen_y.y(),
				     -shift_scale_factor * screen_vectors.screen_y.z());
      }
      if (direction == GDK_KEY_Down) {
	 // std::cout << "Down nudge residue" << std::endl;
	 shift = clipper::Coord_orth(shift_scale_factor * screen_vectors.screen_y.x(),
				     shift_scale_factor * screen_vectors.screen_y.y(),
				     shift_scale_factor * screen_vectors.screen_y.z());
      }

      // all constructed.  Apply it
      clipper::RTop_orth rtop(mat, shift);
      int imol = active_atom.second.first;
      graphics_info_t::molecules[imol].transform_zone_by(active_atom.second.second.chain_id,
							 active_atom.second.second.res_no,
							 active_atom.second.second.res_no,
							 active_atom.second.second.ins_code,
							 rtop, 1);
      graphics_info_t g;

      // If this shift is not added to the rotation centre, we get
      // amusing action when this keypress is repeated.  That should
      // be exported to the scripting layer as an easter egg.
      //
      coot::Cartesian shift_cart(shift.x(), shift.y(), shift.z());
      g.add_vector_to_RotationCentre(shift_cart);
      graphics_draw();
   }
}


// --- nudge (rotate) active residue
// static
void
graphics_info_t::nudge_active_residue_by_rotate(guint direction) {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = graphics_info_t::active_atom_spec();
   if (active_atom.first) {
      graphics_info_t g;
      int imol = active_atom.second.first;
      double angle = M_PI/20;
      if (direction == GDK_KEY_Left)
	 angle = -angle;
      if (direction == GDK_KEY_Up)
	 angle *= 3.0;
      if (direction == GDK_KEY_Down)
	 angle *= -3.0;
      coot::Cartesian rc = g.RotationCentre();
      clipper::Coord_orth origin_offset(rc.x(), rc.y(), rc.z());
      glm::vec3 front_centre = unproject_to_world_coordinates(glm::vec3(0.0f, 0.0f,  1.0f));
      glm::vec3  back_centre = unproject_to_world_coordinates(glm::vec3(0.0f, 0.0f, -1.0f));
      glm::vec3 ftb = back_centre - front_centre;
      ftb = -ftb; // seems more sensible this way round
      clipper::Coord_orth around_vec(ftb.x, ftb.y, ftb.z);
      std::cout << "nudge_active_residue_by_rotate() around_vec " << around_vec.format() << std::endl;
      coot::residue_spec_t res_spec(coot::atom_spec_t(active_atom.second.second));
      g.molecules[imol].rotate_residue(res_spec, around_vec, origin_offset, angle);
      graphics_draw();
   }
}



// this is for the graphics
void
graphics_info_t::execute_db_main() {

   int imol = db_main_imol;
   mmdb::Atom *at_1 = molecules[imol].atom_sel.atom_selection[db_main_atom_index_1];

   // Replace this by the single click, double direction version
   //
   // mmdb::Atom *at2 = molecules[imol].atom_sel.atom_selection[db_main_atom_index_2];
   // std::string chain_id = at1->GetChainID();
   // int iresno_start = at1->GetSeqNum();
   // int iresno_end   = at2->GetSeqNum();

   // Replace this by the single click, double direction version
   //
   // std::string direction_string("forwards"); // forwards
   // execute_db_main(imol, chain_id, iresno_start, iresno_end, direction_string);

   coot::residue_spec_t residue_spec(at_1->GetResidue());
   std::pair<int, int> r = execute_db_main_fragment(imol, residue_spec);
}

// this is called by interactive function and scripting function
//
// return the new molecule number.
int
graphics_info_t::execute_db_main(int imol,
   std::string chain_id,
   int iresno_start,
   int iresno_end,
   std::string direction_string) {

      int imol_new = -1;

      int ilength = 6;
      int idbfrags = 0;

      if (main_chain.is_empty()) {
         idbfrags = main_chain.fill_with_fragments(ilength);
      }

      // should be filled now
      //
      if (main_chain.is_empty()) {
         std::cout << "Sorry cannot do a db fitting without reference structures"
         << std::endl;
         std::string s("Sorry cannot do a main-chain fitting without reference structures");
         wrapped_nothing_bad_dialog(s);
      } else {

         if (iresno_start > iresno_end) {
            int tmp = iresno_end;
            iresno_end = iresno_start;
            iresno_start = tmp;
         }

         mmdb::Manager *mol = molecules[imol].atom_sel.mol;
         if (!mol) return imol_new; // -1

         // mt is a minimol of the Baton Atoms:
         coot::minimol::molecule mt(molecules[imol].atom_sel.mol);
         coot::minimol::molecule target_ca_coords;

         if (direction_string != "backwards") {
            for (unsigned int i=0; i<mt.fragments.size(); i++)
            if (mt.fragments[i].fragment_id == chain_id) {
               std::cout << "not backwards " << mt.fragments[i] << std::endl;
               target_ca_coords.fragments.push_back(mt.fragments[i]);
            }
         } else { // backwards code.

            // Did this ever work!? (It seems to now)

            // std::cout << "---- backwards build" << std::endl;

            for (unsigned int i=0; i<mt.fragments.size(); i++) {
               if (mt[i].fragment_id == chain_id) {

                  // put in the residues of mt.fragments[i] backwards:

                  // The seqnum of the residues is ignored, the only
                  // important thing is the ires.

                  int ifrag = target_ca_coords.fragment_for_chain(chain_id);

                  if (false) {
                     std::cout << "here with ifrag " << ifrag << std::endl;
                     std::cout << "here with max_residue_number " << i << " " << mt[i].max_residue_number() << std::endl;
                  }

                  if (mt[i].max_residue_number() > 1) {
                     int mnr = mt[i].max_residue_number();
                     for (int ires=mnr; ires>=mt[i].min_res_no(); ires--) {

                        int ires_target = mnr-ires+1;

                        if (mt[i][ires].n_atoms() > 0) {
                           coot::minimol::atom ca = mt[i][ires][0];
                           coot::minimol::residue residue(ires_target);
                           residue.addatom(ca);

                           target_ca_coords[ifrag].addresidue(residue, false);
                        }

                     }
                     break;
                  }
               }
            }
         }

         if (direction_string == "backwards") {
            if (target_ca_coords.fragments.size() > 0) {
               iresno_start = target_ca_coords[0].min_res_no();
               iresno_end   = target_ca_coords[0].max_residue_number();
            }
         }

         if (false) {
            std::cout << "Here is target_ca_coords: " << std::endl;
            for(unsigned int ifrag=0; ifrag<target_ca_coords.fragments.size(); ifrag++) {
               for(int ires=target_ca_coords[ifrag].min_res_no(); ires<=target_ca_coords[ifrag].max_residue_number(); ires++) {
                  for (unsigned int iat=0; iat<target_ca_coords[ifrag][ires].atoms.size(); iat++) {
                     std::cout << " " << target_ca_coords[ifrag].fragment_id << " " << ires << " " << target_ca_coords[ifrag][ires]
                     << " " << target_ca_coords[ifrag][ires][iat].name
                     << " " << target_ca_coords[ifrag][ires][iat].pos.format() << std::endl;
                  }
               }
            }
         }

         // now target_ca_coords has only one chain, the chain of the zone.
         // Note that match_target_fragment selects CAs from target_ca_coords
         // so we don't need to filter them out here.


         // write out target_ca_coords:
         // if (direction_string == "backwards")
         // 	 target_ca_coords.write_file("target_ca_coords.pdb", 20);

         main_chain.match_target_fragment(target_ca_coords,
            iresno_start,
            iresno_end,
            ilength);

            float bf = default_new_atoms_b_factor;
            main_chain.merge_fragments();
            coot::minimol::molecule mmol;
            mmol.fragments.push_back(main_chain.mainchain_fragment());

            // if (direction_string == "backwards")
            // 	 mol.write_file("db-mainchain-backwards.pdb", bf);

            // std::cout << "DEBUG:: mol.is_empty() returns " << mol.is_empty() << std::endl;
            std::vector<coot::minimol::atom *> serial_atoms = mmol.select_atoms_serial();
            // std::cout << "DEBUG:: serial_atoms.size() returns " << serial_atoms.size() << std::endl;

            if (serial_atoms.size() > 0) {
               std::pair<std::vector<float>, std::string> cell_spgr =
                  molecules[imol].get_cell_and_symm();
               atom_selection_container_t asc = make_asc(mmol.pcmmdbmanager());
               set_mmdb_cell_and_symm(asc, cell_spgr); // tinker with asc.
               // Consider asc as an object.
               imol_new = create_molecule();
               std::string mol_name = "mainchain-";
               mol_name += direction_string;
               molecules[imol_new].install_model(imol_new, asc, Geom_p(), mol_name, 1);
               graphics_draw();
            } else {
               if (use_graphics_interface_flag) {
                  std::string s("Sorry, failed to convert that residue range.\nToo short, perhaps?");
                  GtkWidget *w = wrapped_nothing_bad_dialog(s);
                  gtk_widget_set_visible(w, TRUE);
               }
            }
            main_chain.clear_results();
         }

         return imol_new;
      }

// build both directions.
std::pair<int, int>
graphics_info_t::execute_db_main_fragment(int imol, coot::residue_spec_t spec) {

   std::pair<int, int> new_mols = std::pair<int, int> (-1, -1);

   if (is_valid_model_molecule(imol)) {

      mmdb::Manager *mol = molecules[imol].atom_sel.mol;
      float dist_max = 4.5; // CA-CA
      mmdb::Residue *residue_start_p = molecules[imol].get_residue(spec);
      if (residue_start_p) {
	 std::vector<mmdb::Residue *> residues =
	    coot::simple_residue_tree(residue_start_p, mol, dist_max);
	 if (residues.size() > 0) {
	    int found_max = -9999;
	    int found_min =  9999;
	    for (std::size_t i=0; i<residues.size(); i++) {
	       int resno_this = residues[i]->GetSeqNum();
	       if (resno_this < found_min) found_min = resno_this;
	       if (resno_this > found_max) found_max = resno_this;
	    }

	    std::cout << "-------------------------------------------------------------" << std::endl;
	    std::cout << "Here with " << found_min << " " << found_max << std::endl;
	    std::cout << "-------------------------------------------------------------" << std::endl;
	    int imol_new_1 = execute_db_main(imol, spec.chain_id, found_min, found_max, "forwards");
	    int imol_new_2 = execute_db_main(imol, spec.chain_id, found_min, found_max, "backwards");

	    std::pair<int, int> nm(imol_new_1, imol_new_2);
	    return nm;
	 }
      }
   }
   return new_mols;
}

// --------------------------------------------------------------------------------
//                 Rotamer stuff
// --------------------------------------------------------------------------------

void
graphics_info_t::do_rotamers(int atom_index, int imol) {

   // Note to self: this function is not called when clicking on the "Rotamers" button in the interface.

   if (use_graphics_interface_flag) {
      // display the buttons for the rotamer options and display
      // the most likely in the graphics as a
      // moving_atoms_asc.

      rotamer_residue_atom_index = atom_index;  // save for button
      // callbacks, so that we
      // can get the residue.
      rotamer_residue_imol = imol;
      mmdb:: Atom *active_atom = molecules[imol].atom_sel.atom_selection[atom_index];
      std::string altconf = active_atom->altLoc;
      std::cout << "debug:: altconf " << altconf << " with length " << altconf.length() << std::endl;
      bool is_alt_conf_dialog = false;
      if (altconf.length() > 0)
         is_alt_conf_dialog = true;

//       std::cout << "DEBUG:: in do_rotamers() atom_index is " << atom_index
//                 << " and alconf is :" <<  altconf << ":" << std::endl;

      // GtkWidget *dialog = create_rotamer_selection_dialog();
      GtkWidget *dialog = widget_from_builder("rotamer_selection_dialog");
      set_transient_and_position(COOT_ROTAMER_SELECTION_DIALOG, dialog);
      g_object_set_data(G_OBJECT(dialog), "imol", GINT_TO_POINTER(imol));

      // Test if this was an alt confed atom.
      // If it was, then we should set up the hscale.
      // It it was not, then we should hide the hscale
      //
      if (is_alt_conf_dialog) {
         // GtkWidget *hscale = lookup_widget(dialog, "new_alt_conf_occ_hscale");
         GtkWidget *hscale = widget_from_builder("new_alt_conf_occ_hscale");
         float v = add_alt_conf_new_atoms_occupancy;
         // The max value is 3rd arg - 6th arg (here 2 and 1 is the same as 1 and 0)
         GtkAdjustment *adj = GTK_ADJUSTMENT(gtk_adjustment_new(v, 0.0, 2.0, 0.01, 0.1, 1.0));
         gtk_range_set_adjustment(GTK_RANGE(hscale), GTK_ADJUSTMENT(adj));
         g_signal_connect(G_OBJECT(adj),
                          "value_changed",
                          G_CALLBACK(graphics_info_t::new_alt_conf_occ_adjustment_changed),
                          NULL);
         g_object_set_data(G_OBJECT(dialog), "type", GINT_TO_POINTER(1));

      } else {

         // GtkWidget *frame = lookup_widget(dialog, "new_alt_conf_occ_frame");
         GtkWidget *frame = widget_from_builder("new_alt_conf_occ_frame");
         gtk_widget_set_visible(frame, FALSE);
         g_object_set_data(G_OBJECT(dialog), "type", GINT_TO_POINTER(0));
         fill_rotamer_selection_buttons(dialog, active_atom, imol);
      }

      /* Capture keypress events */
      //    rotamer_key_press_event is not defined (yet)
      //    gtk_signal_connect(GTK_OBJECT(window), "key_press_event",
      //                       GTK_SIGNAL_FUNC(rotamer_key_press_event), NULL);
      /* set focus to glarea widget - we need this to get key presses. */
      // std::cout << "Focus on the table " << std::endl;
      // GTK_WIDGET_SET_FLAGS(dialog, GTK_CAN_FOCUS);
      gtk_widget_grab_focus(GTK_WIDGET(glareas[0])); // but set focus to the graphics.

      // old function no loger compiles
      // fill_rotamer_selection_buttons(dialog, atom_index, imol);

      // act as if the button for the first rotamer was pressed
      // short int stat = generate_moving_atoms_from_rotamer(0); // 20220812-PE no longer compiles

      if (true) // was stat != 0
         gtk_widget_set_visible(dialog, TRUE);
   }
}

// 20220812-PE gtk4 version
void
graphics_info_t::do_rotamers(int imol, mmdb::Atom *active_atom) {

   if (! use_graphics_interface_flag) return;
   if (! active_atom) {
      std::cout << "ERROR:: in do_rotamers() active_atom is null" << std::endl;
      return;
   }

   rotamer_residue_imol = imol;
   rotamer_residue_atom_index = 0; // we don't know it.
   std::string altconf(active_atom->altLoc);
   bool is_alt_conf_dialog = false; // does the clicked atom have an alt conf?
   if (! altconf.empty()) is_alt_conf_dialog = true;
   rotamer_residue_atom_spec = coot::atom_spec_t(active_atom);

   GtkWidget *dialog = widget_from_builder("rotamer_selection_dialog");
   set_transient_and_position(COOT_ROTAMER_SELECTION_DIALOG, dialog);
   // rotamer_dialog = dialog; // Hmm... this doesn't look good. // 20240304-PE removed.
   g_object_set_data(G_OBJECT(dialog), "imol", GINT_TO_POINTER(imol));
   if (is_alt_conf_dialog) {

   } else {
      GtkWidget *frame = widget_from_builder("new_alt_conf_occ_frame");
      gtk_widget_set_visible(frame, FALSE);
      g_object_set_data(G_OBJECT(dialog), "type", GINT_TO_POINTER(0));
   }
   fill_rotamer_selection_buttons(dialog, active_atom, imol);
   generate_moving_atoms_from_rotamer(imol, rotamer_residue_atom_spec, 0); // passed and data member - not good design
   gtk_widget_set_visible(dialog, TRUE);
}



// static
void graphics_info_t::new_alt_conf_occ_adjustment_changed(GtkAdjustment *adj,
                                                          gpointer user_data) {

   graphics_info_t g;
   g.add_alt_conf_new_atoms_occupancy = gtk_adjustment_get_value(adj);

   // Change the occupancies of the intermediate atoms:
   //
   if (moving_atoms_asc) {
      for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
         // this if test is a kludge!
         // Don't change the alt conf for fully occupied atoms.
         if (moving_atoms_asc->atom_selection[i]->occupancy < 0.99)
            moving_atoms_asc->atom_selection[i]->occupancy = gtk_adjustment_get_value(adj);
      }
   }
}

// static
void
graphics_info_t::drag_intermediate_atom(const coot::atom_spec_t &atom_spec, const clipper::Coord_orth &pt) {

   // std::cout << "DEBUG:: " << atom_spec << " to " << pt.format() << std::endl;

   if (! moving_atoms_asc) {
      std::cout << "WARNING:: No intermediate atoms - fail" << std::endl;
   } else {
      if (! moving_atoms_asc->mol) {
	 std::cout << "WARNING:: No intermediate atoms mol - fail" << std::endl;
      } else {
	 int imod = 1;
	 mmdb::Model *model_p = moving_atoms_asc->mol->GetModel(imod);
	 mmdb::Chain *chain_p;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<nchains; ichain++) {
	    chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    mmdb::PResidue residue_p;
	    mmdb::Atom *at;
	    for (int ires=0; ires<nres; ires++) {
	       residue_p = chain_p->GetResidue(ires);
	       int n_atoms = residue_p->GetNumberOfAtoms();
	       for (int iat=0; iat<n_atoms; iat++) {
		  at = residue_p->GetAtom(iat);
		  if (atom_spec.matches_spec(at)) {
		     at->x = pt.x();
		     at->y = pt.y();
		     at->z = pt.z();
		  }
	       }
	    }
	 }

	 // pre-threaded refinement
// 	 std::set<int> dummy;
// 	 Bond_lines_container bonds(*moving_atoms_asc, imol_moving_atoms, dummy, geom_p, 0, 1, 0);
// 	 regularize_object_bonds_box.clear_up();
// 	 regularize_object_bonds_box = bonds.make_graphical_bonds();
// 	 graphics_draw();

// 	 // now refine (again) with the new atom position
// 	 graphics_info_t g;
// 	 g.add_drag_refine_idle_function();

	 thread_for_refinement_loop_threaded();

      }
   }
}


// static
void
graphics_info_t::while_moving_atoms_active_mark_atom_as_fixed(int imol, const coot::atom_spec_t &atom_spec, bool state) {

   std::cout << "Here in while_moving_atoms_active_mark_atom_as_fixed() --- start --- maa: " << moving_atoms_asc << std::endl;

   // 20211202-PE I don't understand how this test works or why it is here.
   if (moving_atoms_asc) {
      std::cout << "Here in while_moving_atoms_active_mark_atom_as_fixed() 2" << std::endl;
      if ((imol >=0) && (imol < n_molecules())) {
         std::cout << "Here in while_moving_atoms_active_mark_atom_as_fixed() 3" << std::endl;
         if (graphics_info_t::molecules[imol].has_model()) {
            std::cout << "Here in while_moving_atoms_active_mark_atom_as_fixed() 4" << std::endl;
            graphics_info_t::molecules[imol].mark_atom_as_fixed(atom_spec, state);
            graphics_info_t g;
            g.setup_draw_for_anchored_atom_markers();
         }
      }
   } else {
      std::cout << "WARNING:: in while_moving_atoms_active_mark_atom_as_fixed() No intermediate atoms - fail" << std::endl;
   }
}

// static
void
graphics_info_t::mark_atom_as_fixed(int imol, const coot::atom_spec_t &atom_spec, bool state) {

   std::cout << "debug:: mark_atom_as_fixed() --- start --- maa: " << moving_atoms_asc << std::endl;

   if (true) {
      std::cout << "Here in mark_atom_as_fixed() 2" << std::endl;
      if ((imol >=0) && (imol < n_molecules())) {
         std::cout << "Here in mark_atom_as_fixed() 3" << std::endl;
         if (graphics_info_t::molecules[imol].has_model()) {
            std::cout << "Here in mark_atom_as_fixed() 4" << std::endl;
            graphics_info_t::molecules[imol].mark_atom_as_fixed(atom_spec, state);
            graphics_info_t g;
            g.setup_draw_for_anchored_atom_markers();
         }
      }
   }
}

void
graphics_info_t::fill_rotamer_selection_buttons(GtkWidget *dialog, mmdb::Atom *active_atom, int imol) const {

   std::cout << "in fill_rotamer_selection_buttons() with active_atom " << active_atom << std::endl;
   if (! active_atom) return; // error

   // for each rotamer do this:

   GSList *gr_group = NULL;
   GtkWidget *rotamer_selection_dialog = dialog;
   // GtkWidget *rotamer_selection_button_vbox = lookup_widget(window, "rotamer_selection_button_vbox");
   GtkWidget *rotamer_selection_button_vbox = widget_from_builder("rotamer_selection_button_vbox");
   std::string alt_conf = active_atom->altLoc;
   mmdb::Residue *residue = active_atom->residue;

   coot::richardson_rotamer d(residue, alt_conf,
                              molecules[imol].atom_sel.mol, rotamer_lowest_probability, 0);

   std::vector<float> probabilities = d.probabilities();

   if (false)
      std::cout << "debug:: in fill_rotamer_selection_buttons():: There are " << probabilities.size() << " probabilities" << std::endl;

   // Attach the number of residues to the dialog so that we can get
   // that data item when we make a synthetic key press due to
   // keyboard (arrow?) key press:
   g_object_set_data(G_OBJECT(dialog), "probabilities_size", GINT_TO_POINTER(probabilities.size()));

   GtkWidget *radio_group = nullptr;

   for (unsigned int i=0; i<probabilities.size(); i++) {
      std::string button_label = int_to_string(i+1);
      button_label += ":  ";
      button_label += d.rotamer_name(i);
      button_label += "  ";
      button_label += float_to_string(probabilities[i]);
      button_label += "% Chi_1 = ";
      button_label += float_to_string(d.Chi1(i));
      std::string button_name = "rotamer_selection_button_rot_";
      button_name += int_to_string(i);

      GtkWidget *rotamer_selection_radio_button =  nullptr;

      rotamer_selection_radio_button = gtk_check_button_new_with_label(button_label.c_str());
      if (radio_group)
         gtk_check_button_set_group(GTK_CHECK_BUTTON(rotamer_selection_radio_button), GTK_CHECK_BUTTON(radio_group));
      else
         radio_group = rotamer_selection_radio_button;
      // gtk_widget_ref (rotamer_selection_radio_button);
      g_object_set_data_full(G_OBJECT (rotamer_selection_dialog),
                             button_name.c_str(), rotamer_selection_radio_button,
                             NULL);

      int *iuser_data = new int;
      *iuser_data = i;
      g_signal_connect (G_OBJECT(rotamer_selection_radio_button), "toggled",
                        G_CALLBACK(on_rotamer_selection_button_toggled),
                        iuser_data);

       gtk_widget_set_visible (rotamer_selection_radio_button, TRUE);
       GtkWidget *frame = gtk_frame_new(NULL);
       gtk_frame_set_child(GTK_FRAME(frame), rotamer_selection_radio_button);
       gtk_widget_set_margin_start(GTK_WIDGET(frame), 6);
       gtk_widget_set_margin_end(  GTK_WIDGET(frame), 6);
       gtk_widget_set_margin_top(  GTK_WIDGET(frame), 6);
       gtk_widget_set_margin_bottom(GTK_WIDGET(frame), 6);

       gtk_box_append(GTK_BOX(rotamer_selection_button_vbox), frame);
       gtk_widget_set_visible(frame, TRUE);
   }
}


void
graphics_info_t::on_rotamer_selection_button_toggled(GtkCheckButton  *button,
                                                     gpointer         user_data) {

   int *i_tmp = (int *) user_data;
   int i = *i_tmp;

   graphics_info_t g;
   coot::atom_spec_t atom_spec = rotamer_residue_atom_spec;
   g.generate_moving_atoms_from_rotamer(rotamer_residue_imol, atom_spec, i);

}


// Return 1 for valid (i.e. non-GLY, non-ALA) residue, 0 otherwise
// (including residue type not found).
//
short int
graphics_info_t::generate_moving_atoms_from_rotamer(int imol, coot::atom_spec_t &atom_spec, int irot) {

   // int atom_index = rotamer_residue_atom_index;

   mmdb::Atom *at_rot = get_atom(imol, atom_spec);
   if (! at_rot) return 0;
   mmdb::Residue *residue  = at_rot->residue;
   int atom_index_udd = molecules[imol].atom_sel.UDDAtomIndexHandle;
   std::string altconf = at_rot->altLoc;

   if (std::string(residue->name) == "GLY" ||
       std::string(residue->name) == "ALA") {
      std::cout << "INFO:: This residue type ("<< residue->name << ") doesn't have rotamers\n";
      return 0;
   }

   // We need to filter out atoms that are not (either the same
   // altconf as atom_index or "")
   //
   // get rid of this function (needs a test)
   bool embed_in_chain_flag = false;
   mmdb::Residue *tres = coot::deep_copy_this_residue_old_style(residue,
						 std::string(at_rot->altLoc),
						 0, atom_index_udd, embed_in_chain_flag);
   if (!tres) {
      return 0;
   } else {
      mmdb::PPAtom residue_atoms;
      int nResidueAtoms;
      std::string mol_atom_altloc;
      std::string atom_altloc = at_rot->altLoc;
      tres->GetAtomTable(residue_atoms, nResidueAtoms);
      for (int iat=0; iat<nResidueAtoms; iat++) {
	 mol_atom_altloc = std::string(residue_atoms[iat]->altLoc);
	 if (! ((mol_atom_altloc ==  atom_altloc) || (mol_atom_altloc == ""))) {
	    tres->DeleteAtom(iat);
	 }
      }
      tres->TrimAtomTable();

      std::string monomer_type = tres->GetResName();
      std::pair<short int, coot::dictionary_residue_restraints_t> p =
	 Geom_p()->get_monomer_restraints(monomer_type, imol);

      coot::richardson_rotamer d(tres, altconf, molecules[imol].atom_sel.mol,
				 rotamer_lowest_probability, 0);

      if (p.first) {
	 // std::cout << "generate_moving_atoms_from_rotamer " << irot << std::endl;
	 // The magic happens here:
         const coot::dictionary_residue_restraints_t &dict = p.second;
	 mmdb::Residue *moving_res = d.GetResidue(dict, irot);
	 //
	 if (moving_res == NULL) {
	    std::cout << "Failure to find rotamer for residue type: "
		      << residue->name << std::endl;
	    return 0;
	 } else {

	    mmdb::Manager *mol = new mmdb::Manager;
	    mmdb::Model *model_p = new mmdb::Model;
	    mmdb::Chain *chain_p = new mmdb::Chain;
	    mmdb::Residue *res_p = new mmdb::Residue;
	    res_p->SetResID(residue->GetResName(),
			    residue->GetSeqNum(),
			    residue->GetInsCode());

	    mmdb::PPAtom residue_atoms_2 = NULL;
	    int nResidueAtoms_2 = 0;
	    moving_res->GetAtomTable(residue_atoms_2, nResidueAtoms_2);
	    mmdb::Atom *atom_p;
	    int i_add;
	    for(int iat=0; iat<nResidueAtoms_2; iat++) {
	       atom_p = new mmdb::Atom;
	       atom_p->Copy(residue_atoms_2[iat]);
	       i_add = res_p->AddAtom(atom_p);
	    }
	    chain_p->AddResidue(res_p);
	    chain_p->SetChainID(residue->GetChainID());
	    model_p->AddChain(chain_p);
	    mol->AddModel(model_p);
	    mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
	    mol->FinishStructEdit();

	    imol_moving_atoms = imol;

            if (! moving_atoms_asc)
               moving_atoms_asc = new atom_selection_container_t;

	    *moving_atoms_asc = make_asc(mol);
	    //    std::cout << "there are " << moving_atoms_asc->n_selected_atoms
	    // 	     << " selected atoms in the moving_atoms_asc" << std::endl;

	    moving_atoms_asc_type = coot::NEW_COORDS_REPLACE_CHANGE_ALTCONF;
	    make_moving_atoms_graphics_object(imol, *moving_atoms_asc);
	    if (do_probe_dots_on_rotamers_and_chis_flag) {
	       setup_for_probe_dots_on_chis_molprobity(imol);
	    }
	    graphics_draw();
	    return 1;
	 }
      }
   }
   return 0;
}

coot::rotamer_probability_info_t
graphics_info_t::get_rotamer_probability(mmdb::Residue *res,
					 const std::string &altconf,
					 mmdb::Manager *mol,
					 float lowest_probability,
					 short int add_extra_PHE_and_TYR_rotamers_flag) {

   bool debug = false;
   coot::rotamer_probability_info_t r(coot::rotamer_probability_info_t::MISSING_ATOMS,0,"");
   if (!rot_prob_tables.is_well_formatted()) {
      rot_prob_tables.fill_tables();
   }
   if (rot_prob_tables.is_well_formatted()) {
      try {
         std::string res_name(res->GetResName());
         if (coot::util::is_standard_amino_acid_name(res_name)) {
            std::vector<coot::rotamer_probability_info_t> v = rot_prob_tables.probability_this_rotamer(res);
            if (v.size() > 0) {
               r = v[0];
               if (debug)
                  std::cout << "  residue " << coot::residue_spec_t(res) << " " << v[0] << std::endl;
            }
         }

      }
      catch (const std::runtime_error &e) {
	 std::cout << "get_rotamer_probability: caught: " << e.what() << std::endl;
      }
   } else {
      coot::richardson_rotamer d(res, altconf, mol, rotamer_lowest_probability, 1);
      r  = d.probability_of_this_rotamer();
   }

   // flag for assigned,                    1
   // unassigned due to missing atoms,      0
   // unassigned due to rotamer not found. -1
   // unassigned due to GLY/ALA            -2

   return r;
}



// all molecule rotamer score, (depends on private rotamer probability tables)
coot::rotamer_score_t
graphics_info_t::all_molecule_rotamer_score(int imol) const {

   coot::rotamer_score_t rs;

   if (!rot_prob_tables.is_well_formatted()) {
      rot_prob_tables.fill_tables();
   }
   if (rot_prob_tables.is_well_formatted()) {
      if (is_valid_model_molecule(imol)) {
	 rs = graphics_info_t::molecules[imol].get_all_molecule_rotamer_score(rot_prob_tables);
      }
   }
   return rs;
}



// wiggly ligands support
//
std::vector <coot::dict_torsion_restraint_t>
graphics_info_t::get_monomer_torsions_from_geometry(const std::string &monomer_type) const {
   return geom_p->get_monomer_torsions_from_geometry(monomer_type, find_hydrogen_torsions_flag);
}

// Each new atom goes in its own residue.
// Residue type is HOH.
//
void
graphics_info_t::place_dummy_atom_at_pointer() {

   int imol = create_pointer_atom_molecule_maybe();
   molecules[imol].add_pointer_atom(RotationCentre()); // update bonds
   graphics_draw();

}

void
graphics_info_t::place_typed_atom_at_pointer(const std::string &type) {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      if (is_valid_model_molecule(imol)) {
         if (molecules[imol].is_displayed_p()) {
            std::pair<bool, std::string > status_mess =
               molecules[imol].add_typed_pointer_atom(RotationCentre(), type); // update bonds
            update_environment_distances_by_rotation_centre_maybe(imol);
            graphics_draw();
            if (status_mess.first == false) {
               std::string m = "WARNING:: disallowed ";
               m += status_mess.second;
               info_dialog(m);
            }
         } else {
            std::string message = "WARNING:: disallowed addition of ";
            message += type;
            message += "\nas the target molecule is not displayed";
            info_dialog(message);
         }
      }
   }
   graphics_grab_focus(); // 20250615-PE maybe a good idea
}

// Tinker with the atom positions of residue
// Return 1 on success.
// We need to pass the asc for the mol because we need it for seekcontacts()
// Of course the asc that is passed is the moving atoms asc.
//
// nth_chi is 1-based (i.e. rotating about CA-CB, nth_chi is 1).
//
short int
graphics_info_t::update_residue_by_chi_change(int imol, mmdb::Residue *residue,
					      atom_selection_container_t &asc,
					      int nth_chi, double diff) {
   short int istat = 0;
   double angle = diff/60.0;
   bool reverse = edit_chi_angles_reverse_fragment;

   std::string monomer_type = residue->GetResName();
   // this can throw an exception
   std::pair<short int, coot::dictionary_residue_restraints_t> p =
      Geom_p()->get_monomer_restraints(monomer_type, imol);

   if (p.first) {
      try {
	 std::pair<std::string, std::string> atom_names = get_chi_atom_names(residue, p.second, nth_chi);
	 std::string alt_conf = chi_angle_alt_conf;
	 try {
	    coot::atom_tree_t tree(p.second, residue, alt_conf);
	    // this can throw an exception
	    double new_torsion = tree.rotate_about(atom_names.first, atom_names.second, angle, reverse);
	    display_density_level_this_image = 1;
	    display_density_level_screen_string = "  Chi ";
	    display_density_level_screen_string += int_to_string(nth_chi);
	    display_density_level_screen_string += "  =  ";
	    display_density_level_screen_string += float_to_string(new_torsion);
	    add_status_bar_text(display_density_level_screen_string);
	 }
	 catch (const std::runtime_error &rte) {
	    // std::cout << rte.what() << std::endl;
	    int base_atom_index = 0;

	    // tmp hack for testing.
	    // coot::contact_info contact = coot::getcontacts(*moving_atoms_asc);

	    // c.f. get_contact_indices_from_restraints().  urgh. Same
	    // functionality: "written twice".
	    //
	    coot::contact_info contact = coot::getcontacts(*moving_atoms_asc, monomer_type, imol, Geom_p());
	    std::vector<std::vector<int> > contact_indices = contact.get_contact_indices_with_reverse_contacts();

	    try {
	       coot::atom_tree_t tree(contact_indices, base_atom_index, residue, alt_conf);
	       // this can throw an exception
	       double new_torsion = tree.rotate_about(atom_names.first, atom_names.second, angle, reverse);
	       display_density_level_this_image = 1;
	       display_density_level_screen_string = "  Chi ";
	       display_density_level_screen_string += int_to_string(nth_chi);
	       display_density_level_screen_string += "  =  ";
	       display_density_level_screen_string += float_to_string(new_torsion);
	       add_status_bar_text(display_density_level_screen_string);
	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << "Update chi - contact fall-back fails - " << rte.what() << std::endl;
	    }
	 }
      }
      catch (const std::runtime_error &rte) {
	 // atoms of the torsion not found.
	 std::cout << rte.what() << std::endl;
      }
   } else {

      // chi angles with no dictionary torsions.  No thanks.
      //
      // But hmmm... maybe I should...  :)

   }

   return istat;
}

// this can throw an std::runtime_error exception.
//
std::pair<std::string, std::string>
graphics_info_t::get_chi_atom_names(mmdb::Residue *residue,
				    const coot::dictionary_residue_restraints_t &rest,
				    int nth_chi) const {

   std::pair<std::string, std::string> r(" CA ", " CB "); // PDBv3 FIXME
   int torsion_index = nth_chi -1;
   std::vector <coot::dict_torsion_restraint_t> torsion_restraints =
      rest.get_non_const_torsions(find_hydrogen_torsions_flag);

   if ((torsion_index >=0) && (torsion_index < int(torsion_restraints.size()))) {
      r = std::pair<std::string, std::string> (torsion_restraints[torsion_index].atom_id_2(),
					       torsion_restraints[torsion_index].atom_id_3());
   } else {
      std::string mess = "No torsion found with index ";
      mess += coot::util::int_to_string(torsion_index);
      mess += " in ";
      mess += rest.residue_info.three_letter_code;
      std::runtime_error rte(mess);
      throw rte;
   }
   return r;
}




#include "pulse-data.hh"

// Called by mouse motion callback (in_edit_chi_mode_flag)
//
void
graphics_info_t::rotate_chi(double x, double y) {

   auto atom_to_glm = [] (mmdb::Atom *at) { return glm::vec3(at->x, at->y, at->z); };

   auto get_residue_from_moving_mol = [] () {

      mmdb::Residue *residue_p = nullptr;
      if (! moving_atoms_asc) {
         std::cout << "ERROR: moving_atoms_asc is NULL" << std::endl;
      } else {
         if (moving_atoms_asc->n_selected_atoms == 0) {
            std::cout << "ERROR: no atoms in moving_atoms_asc" << std::endl;
         } else {
            mmdb::Model *model_p = moving_atoms_asc->mol->GetModel(1);
            if (model_p) {
               mmdb::Chain *chain_p = model_p->GetChain(0);
               if (chain_p) {
                  if (chain_p->GetNumberOfResidues() > 0)
                     residue_p = chain_p->GetResidue(0);
               }
            }
         }
      }
      return residue_p;
   };

   auto setup_invalid_chi_angle_pulse = [atom_to_glm] (mmdb::Residue *residue_p) {

      bool broken_line_mode = false;
      lines_mesh_for_identification_pulse.setup_red_pulse(broken_line_mode);
      pulse_data_t *pulse_data = new pulse_data_t(0, 6);
      gpointer user_data = reinterpret_cast<void *>(pulse_data);
      std::vector<glm::vec3> positions;
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            positions.push_back(atom_to_glm(at));
         }
      }
      delete_item_pulse_centres = positions;
      gtk_widget_add_tick_callback(glareas[0], generic_pulse_function, user_data, NULL);
   };

      // the displacement of the mouse is the change in speed of the rotation
   // it's fun. Maybe tricky and conter-intuitive.

   // real values start at 1:
   if (edit_chi_current_chi <= 0) {

      mmdb::Residue *residue_p = get_residue_from_moving_mol();
      if (residue_p) {
         setup_invalid_chi_angle_pulse(residue_p);
         graphics_info_t::ephemeral_overlay_label("select_a_chi_angle_label");
      }
      return;
   }

   mouse_current_x = x;
   mouse_current_y = y;
   double diff = x + y;


   // diff *= 15;
   diff *= 10; // 20230519-PE slow it down, because it's acceleration

   int chi = edit_chi_current_chi;
   // std::cout << "graphics_info_t::rotate_chi " << chi << " by " << diff << std::endl;

   // c.f. generate_moving_atoms_from_rotamer(i), except here we will
   // not be changing our moving_atoms_asc, just updating the atom
   // positions.
   //


   short int istat = 1; // failure
   if (! moving_atoms_asc) {
      std::cout << "ERROR: moving_atoms_asc is NULL" << std::endl;
   } else {
      if (moving_atoms_asc->n_selected_atoms == 0) {
	 std::cout << "ERROR: no atoms in moving_atoms_asc" << std::endl;
      } else {
	 mmdb::Model *model_p = moving_atoms_asc->mol->GetModel(1);
	 if (model_p) {
	    mmdb::Chain *chain_p = model_p->GetChain(0);
	    if (chain_p) {
	       mmdb::Residue *residue_p = chain_p->GetResidue(0);
	       if (residue_p) {
		  istat = update_residue_by_chi_change(imol_moving_atoms, residue_p, *moving_atoms_asc, chi, diff);
	       }
	    }
	 }
      }
   }

   if (istat == 0) {
      // std::cout << "regenerating object" << std::endl;
      regularize_object_bonds_box.clear_up();
      make_moving_atoms_graphics_object(imol_moving_atoms, *moving_atoms_asc); // make new bonds
      graphics_draw();

      //    } else {
      // std::cout << "chi rotate failed  - not regenerating object" << std::endl;

   }
}

// Called by mouse motion callback (in_edit_chi_mode_flag)
//
void
graphics_info_t::rotate_chi_torsion_general(double x, double y) {

   mouse_current_x = x;
   mouse_current_y = y;
   double diff = mouse_current_x - GetMouseBeginX();
   diff += mouse_current_y - GetMouseBeginY();
   diff *= 0.5;

   std::vector<coot::atom_spec_t> specs_local = graphics_info_t::torsion_general_atom_specs;

   short int istat = 1; // failure
   if (! moving_atoms_asc) {
      std::cout << "ERROR:: No moving atoms in rotate_chi_torsion_general" << std::endl;
   } else {
      mmdb::Residue *residue_p = get_first_res_of_moving_atoms();
      if (residue_p) {

         if (specs_local.size() >= 3) {

            std::string altconf = chi_angle_alt_conf;
            try {
               // use class variable (previous saved)
               int base_atom_index = 0;
               coot::atom_tree_t tree(torsion_general_contact_indices, base_atom_index, residue_p, altconf);
               tree.rotate_about(specs_local[1].atom_name, specs_local[2].atom_name,
                                 diff, torsion_general_reverse_flag);
               regularize_object_bonds_box.clear_up();
               make_moving_atoms_graphics_object(imol_moving_atoms, *moving_atoms_asc);
               graphics_draw();
            }
            catch (const std::runtime_error &rte) {
               std::cout << "INFO:: tree by contacts failed " << rte.what() << std::endl;
            }
         } else {
            std::cout << "ERROR:: specs_local size is " << specs_local.size() << std::endl;
         }
      }
   }
}

void
graphics_info_t::rotate_multi_residue_torsion(double x, double y) {

   mouse_current_x = x;
   mouse_current_y = y;
   double diff = mouse_current_x - GetMouseBeginX();
   diff += mouse_current_y - GetMouseBeginY();
   diff *= 0.5; // angle (in degrees).

   if (! moving_atoms_asc->mol) {
      std::cout << "ERROR:: called rotate_multi_residue_torsion() but no moving mol"
		<< std::endl;
   } else {

      std::vector<mmdb::Residue *> residues;
      for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
	 mmdb::Residue *r = moving_atoms_asc->atom_selection[i]->residue;
	 if (std::find(residues.begin(), residues.end(), r) == residues.end())
	    residues.push_back(r);
      }

      std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > link_bond_atom_pairs =
	 coot::torsionable_link_bonds(residues, moving_atoms_asc->mol, Geom_p());
      coot::contact_info contacts(*moving_atoms_asc, imol_moving_atoms, geom_p, link_bond_atom_pairs);
      std::vector<std::vector<int> > contact_indices =
	 contacts.get_contact_indices_with_reverse_contacts();
      try {
	 coot::atom_tree_t tree(contact_indices, 0,
				moving_atoms_asc->mol,
				moving_atoms_asc->SelectionHandle);

	 int index_1 = multi_residue_torsion_rotating_atom_index_pair.first;
	 int index_2 = multi_residue_torsion_rotating_atom_index_pair.second;
	 tree.rotate_about(index_1, index_2, diff, multi_residue_torsion_reverse_fragment_mode);
	 make_moving_atoms_graphics_object(imol_moving_atoms, *moving_atoms_asc);
	 graphics_draw();
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "WARNING:: " << rte.what() << std::endl;
      }
   }
}



// 	 //  debug:
// 	 std::cout << "DEBUG:: residue_mol: ----------------- " << std::endl;
// 	 mmdb::Model *model_p = residues_mol->GetModel(1);
// 	 mmdb::Chain *chain_p;
// 	 int nchains = model_p->GetNumberOfChains();
// 	 std::cout << "DEBUG:: residue_mol: nchains " << nchains << std::endl;
// 	 for (int ichain=0; ichain<nchains; ichain++) {
// 	    chain_p = model_p->GetChain(ichain);
// 	    int nres = chain_p->GetNumberOfResidues();
// 	    for (int ires=0; ires<nres; ires++) {
// 	       mmdb::PResidue residue_p = chain_p->GetResidue(ires);
// 	       std::cout << "DEBUG:: residue " << residue_p->GetChainID()
// 			 << " " << residue_p->GetSeqNum()
// 			 << " " << residue_p->name
// 			 << std::endl;
// 	    }
// 	 }



// If we are splitting a residue, we may need move the altLoc of the
// existing residue from "" to "A". Let's create a new enumerated
// constant NEW_COORDS_INSERT_CHANGE_ALTCONF to flag that.
//
// What's in the residue     What we clicked   Old Coordinates   New Coordinates
//      ""                        ""                "" -> "A"       "B"
//    "A" "B"                     "A"              no change        "C"
//    "A" "B"                     "B"              no change        "C"
//    "" "A" "B"                  ""               [1]              "C"
//    "" "A" "B"                  "A"              [1]              "C"
//    "" "A" "B"                  "B"              [1]              "C"
//
// [1] depends on the split:
//     whole residue split: "" -> "A" , "A" and "B" remain the same
//     partial split:       no change
//
// Now that I think about it, it doesn't matter which atom we click.
//
// ....Oh but it does if we want to follow Stefano's suggestion and
// split at clicked residue...
//
// This is a molecule-class-info function.  What is it doing here?
// It's not here any more.  This is just a wrapper.
//
std::pair<bool,std::string>
graphics_info_t::split_residue(int imol, int atom_index) {

   std::pair<bool, std::string> p(0,"");
   // do moving molecule atoms:
   // short int do_intermediate_atoms = 0;

   // Actually, we don't want intermediate atoms in the usual case.
   //
   // We *do* want intermediate atoms if the user has set the flag so,
   // or if there are not all the necessary (mainchain) atoms to do a
   // rotamer.
   //
   // What are the issues for split position?  None.  We do however
   // need to know the what the alt conf (and atom spec) of a newly
   // created alt conf atom.

   if (imol<n_molecules()) {
      if (molecules[imol].has_model()) {
	 p = graphics_info_t::molecules[imol].split_residue(atom_index, alt_conf_split_type);
	 graphics_draw();
      } else {
	 std::cout << "WARNING:: split_residue: molecule has no model.\n";
      }
   } else {
      std::cout << "WARNING:: split_residue: no such molecule.\n";
   }
   return p;
}



// a wrapper to the lower-level split_residue() that uses an atom index
std::pair<bool,std::string>
graphics_info_t::split_residue(int imol, const std::string &chain_id,
			       int resno,
			       const std::string &ins_code,
			       const std::string &altconf) {

   std::pair<bool, std::string> p(0, "");

   mmdb::Residue *r = molecules[imol].get_residue(chain_id, resno, ins_code);
   if (!r) {
      std::cout << "WARNING:: Residue " << " chain-id :" << chain_id << ":  resno: " << resno
		<< " inscode :" << ins_code << ": not found" << std::endl;
   } else {
      mmdb::PPAtom residue_atoms;
      int n_residue_atoms;
      int at_index = -1;
      r->GetAtomTable(residue_atoms, n_residue_atoms);
      std::cout << "DEBUG:: split_residue table " << std::endl;
      for (int i=0; i<n_residue_atoms; i++) {
	 std::string atom_name(residue_atoms[i]->name);
	 std::string atom_alt_conf(residue_atoms[i]->altLoc);
	 std::cout << "   " << i << " " << atom_name << " :" << atom_alt_conf << ":" << std::endl;
	 if (atom_alt_conf == altconf) {
	    mmdb::Atom *at = residue_atoms[i];
	    int atom_index_udd = molecules[imol].atom_sel.UDDAtomIndexHandle;
	    int n_atoms = molecules[imol].atom_sel.n_selected_atoms;
	    at->GetUDData(atom_index_udd, at_index);
	    if (at_index >= 0 && at_index < n_atoms) {
	       break;
	    }
	 }
      }

      if (at_index != -1) {
	 p = split_residue(imol, at_index);
      } else {
	 std::cout << "WARNING:: atom without atom index in molecule: "
		   << imol << " chain-id :" << chain_id << ":  resno: " << resno << " inscode :"
		   << ins_code << ": altconf :" << altconf << ":"
		   << " split_residue() abandoned."
		   << std::endl;
      }
   }
   return p;
}

void
graphics_info_t::split_residue_range(int imol, int index_1, int index2) {

}


// delete zone
void
graphics_info_t::delete_residue_range(int imol,
				      const coot::residue_spec_t &res1_in,
				      const coot::residue_spec_t &res2_in) {

   if (is_valid_model_molecule(imol)) {

      coot::residue_spec_t res1 = res1_in;
      coot::residue_spec_t res2 = res2_in;

      if (res1.res_no > res2.res_no)
	 std::swap(res1, res2);

      molecules[imol].delete_zone(res1, res2);

      // cheap! I should find the residues with insertion codes in this range too.
      // How to do that? Hmm... Needs a class function. This will do for now
      //
      std::vector<coot::residue_spec_t> res_specs;
      for (int i=res1.res_no; i<=res2.res_no; i++) {
	 coot::residue_spec_t r(res1_in.chain_id, i, "");
	 res_specs.push_back(r);
      }
      delete_residues_from_geometry_graphs(imol, res_specs);

      // 20230426-PE GTK4: delete_molecule_from_from_display_manager() no longer works
      // if (! is_valid_model_molecule(imol))
      // delete_molecule_from_from_display_manager(imol, false);

      if (delete_item_widget) { // what is this?
	 // GtkWidget *checkbutton = lookup_widget(graphics_info_t::delete_item_widget, "delete_item_keep_active_checkbutton");
	 GtkWidget *checkbutton = widget_from_builder("delete_item_keep_active_checkbutton");
	 if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbutton))) {
	    // don't destroy it.
	 } else {
	    gint upositionx, upositiony;
            std::cout << "GTK-FIXME gdk_window_get_root_origin A " << std::endl;
	    // gdk_window_get_root_origin (delete_item_widget->window, &upositionx, &upositiony);
	    // delete_item_widget_x_position = upositionx;
	    // delete_item_widget_y_position = upositiony;
	    gtk_widget_set_visible(delete_item_widget, FALSE);
	    delete_item_widget = 0;
	    normal_cursor();
	 }
      }

      if (graphics_info_t::go_to_atom_window)
	 update_go_to_atom_window_on_changed_mol(imol);

      // faster is passing a blank asc, but to do that needs to check that
      // updating other geometry graphs will work (not crash) with residues/mol
      // unset.
      //
      // atom_selection_container_t asc = molecules[imol].atom_sel;
      atom_selection_container_t asc;
      update_validation(imol);
   }
   graphics_draw();
}

void
graphics_info_t::delete_sidechain_range(int imol,
					const coot::residue_spec_t &res_1,
					const coot::residue_spec_t &res_2) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].delete_sidechain_range(res_1, res_2);
      if (delete_item_widget) {
	 // GtkWidget *checkbutton = lookup_widget(graphics_info_t::delete_item_widget, "delete_item_keep_active_checkbutton");
	 GtkWidget *checkbutton = widget_from_builder("delete_item_keep_active_checkbutton");
	 if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbutton))) {
	    // don't destroy it.
	 } else {
	    // gtk_widget_destroy(delete_item_widget);
	    gtk_widget_set_visible(delete_item_widget, FALSE);
	    delete_item_widget = 0;
	    normal_cursor();
	 }
      }

      if (graphics_info_t::go_to_atom_window)
	 update_go_to_atom_window_on_changed_mol(imol);

      // faster is passing a blank asc, but to do that needs to check that
      // updating other geometry graphs will work (not crash) with residues/mol
      // unset.
      // It seems that I have done that now.
      //
      // atom_selection_container_t asc = molecules[imol].atom_sel;
      // atom_selection_container_t asc;
      update_validation(imol);
   }
   graphics_draw();

}

void
graphics_info_t::delete_active_residue() {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > aa = active_atom_spec();
   if (aa.first) {
      int imol = aa.second.first;
      coot::residue_spec_t res_spec(aa.second.second);
      molecules[imol].delete_residue(res_spec);
      delete_residue_from_geometry_graphs(imol, res_spec);
   }
   graphics_draw();
}

void
graphics_info_t::delete_active_residue_alt_conf_atoms() {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > aa = active_atom_spec();
   if (aa.first) {
      int imol = aa.second.first;
      const coot::atom_spec_t &atom_spec = aa.second.second;
      coot::residue_spec_t res_spec(aa.second.second);
      molecules[imol].delete_residue_with_full_spec(atom_spec.model_number, atom_spec.chain_id, atom_spec.res_no,
                                                    atom_spec.ins_code, atom_spec.alt_conf);
      if (atom_spec.alt_conf.empty())
         delete_residue_from_geometry_graphs(imol, res_spec);
   }
   graphics_draw();
}




// // static
// void
// graphics_info_t::move_molecule_here_item_select(GtkWidget *item,
// 						GtkPositionType pos) {
//    std::cout << "----------------- move_molecule_here_item_select! " << pos << std::endl;
//    graphics_info_t::move_molecule_here_molecule_number = pos;
// }

// static
void graphics_info_t::move_molecule_here_combobox_changed(GtkWidget *combobox, gpointer data) {

   GtkTreeIter iter;
   gboolean state = gtk_combo_box_get_active_iter(GTK_COMBO_BOX(combobox), &iter);
   if (state) {
      GtkTreeModel *model = gtk_combo_box_get_model(GTK_COMBO_BOX(combobox));
      GValue label_as_value = { 0, };
      gtk_tree_model_get_value(model, &iter, 0, &label_as_value);
      int imol = g_value_get_int(&label_as_value);
      std::cout << "move_molecule_here_combobox_changed() imol: " << imol << std::endl;
      graphics_info_t::move_molecule_here_molecule_number = imol;
   } else {
      std::cout << "bad state" << std::endl;
   }
}



void
graphics_info_t::do_probe_dots_on_rotamers_and_chis() {

   do_interactive_probe();
}


void
graphics_info_t::do_interactive_probe() const {

   // we need to test for GUILE_GTK use, because interactive-guile is
   // defined in a file that depends on guile-gtk (and of course, if
   // we do not have guile-gtk, then that file is not loaded).

#if defined USE_GUILE && defined USE_GUILE_GTK && !defined WINDOWS_MINGW
   if (moving_atoms_asc->n_selected_atoms > 0) {
      if (moving_atoms_asc->mol) {
	 moving_atoms_asc->mol->WritePDBASCII("molprobity-tmp-moving-file.pdb");
	 std::string c = "(";
	 c += "interactive-probe ";
	 c += float_to_string(probe_dots_on_chis_molprobity_centre.x());
	 c += " ";
	 c += float_to_string(probe_dots_on_chis_molprobity_centre.y());
	 c += " ";
	 c += float_to_string(probe_dots_on_chis_molprobity_centre.z());
	 c += " ";
	 c += float_to_string(probe_dots_on_chis_molprobity_radius);
	 c += " \"";
	 c += moving_atoms_asc->atom_selection[0]->GetChainID();
	 c += "\" ";
	 c += int_to_string(moving_atoms_asc->atom_selection[0]->GetSeqNum());
	 c += ")";
	 std::cout << "interactive probe debug: " << c << std::endl;
	 scm_c_eval_string(c.c_str());
      }
   }

#else
#ifdef USE_PYTHON

   if (moving_atoms_asc->n_selected_atoms > 0) {
      if (moving_atoms_asc->mol) {
	 moving_atoms_asc->mol->WritePDBASCII("molprobity-tmp-moving-file.pdb");
	 std::string c = "";
	 c += "interactive_probe(";
	 c += float_to_string(probe_dots_on_chis_molprobity_centre.x());
	 c += ", ";
	 c += float_to_string(probe_dots_on_chis_molprobity_centre.y());
	 c += ", ";
	 c += float_to_string(probe_dots_on_chis_molprobity_centre.z());
	 c += ", ";
	 c += float_to_string(probe_dots_on_chis_molprobity_radius);
	 c += ", \"";
	 c += moving_atoms_asc->atom_selection[0]->GetChainID();
	 c += "\", ";
	 c += int_to_string(moving_atoms_asc->atom_selection[0]->GetSeqNum());
	 c += ")";
	 PyRun_SimpleString((char *) c.c_str());
      }
   }

#endif // USE_PYTHON
#endif // USE_GUILE

}

void
graphics_info_t::check_and_warn_inverted_chirals_and_cis_peptides() const {

   if (moving_atoms_asc) {
      if (moving_atoms_asc_type == coot::NEW_COORDS_REPLACE ||
	  moving_atoms_asc_type == coot::NEW_COORDS_REPLACE_CHANGE_ALTCONF) { // needed?
	 if (moving_atoms_asc->mol) {

	    std::string message_string = "Unset";

	    // ================= chirals ================================

	    std::pair<std::vector<std::string> , std::vector <coot::atom_spec_t> >
	       bv = coot::inverted_chiral_volumes(imol_moving_atoms, moving_atoms_asc->mol, geom_p,
						  cif_dictionary_read_number);
	    if (bv.second.size() > 0) {
	       if (bv.second.size() == 1) {
		  int i = 0;
		  message_string = "There is one residue with an\n";
		  message_string += "incorrect chiral volume\n";
		  message_string += bv.second[i].chain_id;
		  message_string += " ";
		  message_string += coot::util::int_to_string(bv.second[i].res_no);
		  message_string += bv.second[i].ins_code;
		  message_string += " ";
		  message_string += bv.second[i].atom_name;
		  message_string += " ";
		  message_string += bv.second[i].alt_conf;
		  message_string += "\n";
	       } else {
		  message_string = "There are ";
		  message_string += coot::util::int_to_string(bv.second.size());
		  message_string += " residues with \n";
		  message_string += "incorrect chiral volumes\n";
		  for (unsigned int i=0; i<bv.second.size(); i++) {
		     message_string += bv.second[i].chain_id;
		     message_string += " ";
		     message_string += coot::util::int_to_string(bv.second[i].res_no);
		     message_string += bv.second[i].ins_code;
		     message_string += " ";
		     message_string += bv.second[i].atom_name;
		     message_string += " ";
		     message_string += bv.second[i].alt_conf;
		     message_string += "\n";
		  }
	       }
	    }

	    // ================== cis peptides ==========================

	    std::vector<coot::util::cis_peptide_info_t> cis_pep_info_vec =
	       coot::util::cis_peptides_info_from_coords(moving_atoms_asc->mol);

	    int n_cis = cis_pep_info_vec.size();

	    if (false)
	       std::cout << "here with n_cis " << n_cis << " and g.moving_atoms_n_cis_peptides"
			 << graphics_info_t::moving_atoms_n_cis_peptides << std::endl;

	    if (n_cis > graphics_info_t::moving_atoms_n_cis_peptides) {
	       if (message_string == "Unset")
		  message_string.clear();
	       if (n_cis == 1) {
		  message_string += "\nWARNING: A cis-peptide ";
		  message_string += cis_pep_info_vec[0].string();
		  message_string += " has been introduced\n";
	       } else {
		  if ((n_cis - graphics_info_t::moving_atoms_n_cis_peptides) > 1) {
		     message_string += "\nWARNING: Extra cis-peptides have been introduced\n";
		     message_string += "\nWARNING: We now have these cis-peptides:\n";
		     for (unsigned int i=0; i<cis_pep_info_vec.size(); i++) {
			message_string += cis_pep_info_vec[i].string();
			message_string += "\n";
		     }
		  } else {
		     message_string += "\nWARNING: We now have these cis-peptides:\n";
		     for (unsigned int i=0; i<cis_pep_info_vec.size(); i++) {
			message_string += cis_pep_info_vec[i].string();
			message_string += "\n";
		     }
		  }
	       }
	    }

	    // ===================== show it? ==============================

	    if (show_chiral_volume_errors_dialog_flag) {
#if 0 // 20230923-PE accept_reject_dialog is not a thing now. I should delete this widget
	       if (accept_reject_dialog) {
		  if (message_string != "Unset") {
		     if (false)
			std::cout << "debug:: here in check_and_warn_inverted_chirals_and_cis_peptides() A calling "
				  << "update_accept_reject_dialog_with_results() with message string \""
				  << message_string << "\"" << std::endl;
		     update_accept_reject_dialog_with_results(accept_reject_dialog,
							      coot::CHIRAL_CENTRES,
							      coot::refinement_results_t(message_string));
		  } else {
		     coot::refinement_results_t rr("");
                     
		     if (false)
			std::cout << "debug:: here in check_and_warn_inverted_chirals_and_cis_peptides() B calling "
				  << "update_accept_reject_dialog_with_results() with rr.info \""
				  << rr.info_text << "\"" << std::endl;
		     update_accept_reject_dialog_with_results(accept_reject_dialog,
							      coot::CHIRAL_CENTRES,
							      rr);
		  }
	       }
	       if (message_string != "Unset") {
		  std::cout << message_string << std::endl;
	       }
#endif
	    }
	 }
      }
   }
}



// calling function must have the restraints lock.
void
graphics_info_t::tabulate_geometric_distortions(coot::restraints_container_t &rr) const {

   coot::geometry_distortion_info_container_t gdic = rr.geometric_distortions();

   std::ofstream f("coot-refinement-debug.tab");

   if (f) {

      std::vector<std::pair<double, std::string> > rest_info;
      for (unsigned int ii=0; ii<gdic.geometry_distortion.size(); ii++) {
	 const coot::geometry_distortion_info_t &gd = gdic.geometry_distortion[ii];
	 const coot::simple_restraint &rest = gd.restraint;

	 if (rest.restraint_type == coot::BOND_RESTRAINT) {
	    std::string s = "bond  " + coot::util::float_to_string(gd.distortion_score);
	    for (unsigned int iat=0; iat<gd.atom_indices.size(); iat++)
	       s += " " + rr.get_atom_spec(gd.atom_indices[iat]).format();
	    s += " indices: ";
	    for (unsigned int iat=0; iat<gd.atom_indices.size(); iat++)
	       s += " " + coot::util::int_to_string(gd.atom_indices[iat]);
	    s += " target: ";
	    s += coot::util::float_to_string(rest.target_value);
	    s += " sigma: ";
	    s += coot::util::float_to_string(rest.sigma);
	    rest_info.push_back(std::pair<double, std::string> (gd.distortion_score, s));
	 }
	 if (rest.restraint_type == coot::ANGLE_RESTRAINT) {
	    std::string s = "angle " + coot::util::float_to_string(gd.distortion_score);
	    for (unsigned int iat=0; iat<gd.atom_indices.size(); iat++)
	       s += " " + rr.get_atom_spec(gd.atom_indices[iat]).format();
	    s += " indices: ";
	    for (unsigned int iat=0; iat<gd.atom_indices.size(); iat++)
	       s += " " + coot::util::int_to_string(gd.atom_indices[iat]);
	    s += " target: ";
	    s += coot::util::float_to_string(rest.target_value);
	    s += " ";
	    s += coot::util::float_to_string(rest.sigma);
	    rest_info.push_back(std::pair<double, std::string> (gd.distortion_score, s));
	 }
	 if (rest.restraint_type == coot::TORSION_RESTRAINT) {
	    std::string s = "torsion ";
	    s += coot::util::float_to_string_using_dec_pl(gd.distortion_score, 4);
	    s += " ";
	    s += " " + rr.get_atom_spec(rest.atom_index_1).format();
	    s += " ";
	    s += " " + rr.get_atom_spec(rest.atom_index_2).format();
	    s += " ";
	    s += " " + rr.get_atom_spec(rest.atom_index_3).format();
	    s += " ";
	    s += " " + rr.get_atom_spec(rest.atom_index_4).format();
	    s += " idx: ";
	    s += coot::util::int_to_string(ii);
	    s += " target: ";
	    s += coot::util::float_to_string(rest.target_value);
	    rest_info.push_back(std::pair<double, std::string> (gd.distortion_score, s));
	 }
	 if (rest.restraint_type == coot::TRANS_PEPTIDE_RESTRAINT) {
	    std::string s = "trans " + coot::util::float_to_string(gd.distortion_score);
	    s += " " + rr.get_atom_spec(rest.atom_index_2).format();
	    s += " " + rr.get_atom_spec(rest.atom_index_3).format();
	    s += " ";
	    s += coot::util::float_to_string(rest.target_value);
	    s += " ";
	    rest_info.push_back(std::pair<double, std::string> (gd.distortion_score, s));
	 }
	 if (rest.restraint_type == coot::PLANE_RESTRAINT) {
	    std::string s = "plane " + coot::util::float_to_string(gd.distortion_score);
	    for (unsigned int iat=0; iat<gd.atom_indices.size(); iat++)
	       s += " " + rr.get_atom_spec(gd.atom_indices[iat]).format();
	    rest_info.push_back(std::pair<double, std::string> (gd.distortion_score, s));
	 }
	 if (rest.restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT) {
	    std::string s = "nbc   " + coot::util::float_to_string(gd.distortion_score);
	    for (unsigned int iat=0; iat<gd.atom_indices.size(); iat++)
	       s += " " + rr.get_atom_spec(gd.atom_indices[iat]).format();
	    s += " indices: ";
	    for (unsigned int iat=0; iat<gd.atom_indices.size(); iat++)
	       s += " " + coot::util::int_to_string(gd.atom_indices[iat]);
	    s += " target: ";
	    s += coot::util::float_to_string(rest.target_value);
	    s += " ";
	    s += coot::util::float_to_string(rest.sigma);
	    rest_info.push_back(std::pair<double, std::string> (gd.distortion_score, s));
	 }
	 if (rest.restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) {
	    std::string s = "chiral " + coot::util::float_to_string(gd.distortion_score);
	    for (unsigned int iat=0; iat<gd.atom_indices.size(); iat++)
	       s +=  " " + rr.get_atom_spec(gd.atom_indices[iat]).format();
	    rest_info.push_back(std::pair<double, std::string> (gd.distortion_score, s));
	 }
	 if (rest.restraint_type == coot::RAMACHANDRAN_RESTRAINT) {
	    std::string s = "rama " + coot::util::float_to_string(gd.distortion_score);
	    for (unsigned int iat=0; iat<gd.atom_indices.size(); iat++)
	       s +=  " " + rr.get_atom_spec(gd.atom_indices[iat]).format();
	    rest_info.push_back(std::pair<double, std::string> (gd.distortion_score, s));
	 }
	 if (rest.restraint_type == coot::START_POS_RESTRAINT) {
	    std::string s = "start-pos " + coot::util::float_to_string(gd.distortion_score);
	    for (unsigned int iat=0; iat<gd.atom_indices.size(); iat++)
	       s += " " + rr.get_atom_spec(gd.atom_indices[iat]).format();
	    rest_info.push_back(std::pair<double, std::string> (gd.distortion_score, s));
	 }
	 if (rest.restraint_type == coot::PARALLEL_PLANES_RESTRAINT) {
	    std::string s = "parallel-plane " + coot::util::float_to_string(gd.distortion_score);
	    for (unsigned int iat=0; iat<gd.atom_indices.size(); iat++)
	       s += " " + rr.get_atom_spec(gd.atom_indices[iat]).format();
	    rest_info.push_back(std::pair<double, std::string> (gd.distortion_score, s));
	 }
	 if (rest.restraint_type == coot::GEMAN_MCCLURE_DISTANCE_RESTRAINT) {
	    std::string s = "geman-mcclure " + coot::util::float_to_string(gd.distortion_score);
	    for (unsigned int iat=0; iat<gd.atom_indices.size(); iat++)
	       s += " " + rr.get_atom_spec(gd.atom_indices[iat]).format();
	    s += " indices: ";
	    for (unsigned int iat=0; iat<gd.atom_indices.size(); iat++)
	       s += " " + coot::util::int_to_string(gd.atom_indices[iat]);
	    s += " target: ";
	    s += coot::util::float_to_string(rest.target_value);
	    rest_info.push_back(std::pair<double, std::string> (gd.distortion_score, s));
	 }
	 if (rest.restraint_type == coot::TARGET_POS_RESTRAINT) {
	    std::string s = "pull-atom " + coot::util::float_to_string(gd.distortion_score);
	    s += " " + rr.get_atom_spec(gd.atom_indices[0]).format();
	    rest_info.push_back(std::pair<double, std::string> (gd.distortion_score, s));
	 }
      }
      std::sort(   rest_info.begin(), rest_info.end());
      std::reverse(rest_info.begin(), rest_info.end());
      for (unsigned int i=0; i<rest_info.size(); i++) {
	 f << rest_info[i].second + "\n";
      }
      f.close();
   }
}


void
graphics_info_t::auto_fit_rotamer_ng(int imol, const coot::residue_spec_t &res_spec, const std::string &alt_conf) {

   int imol_map = Imol_Refinement_Map();
   if (is_valid_map_molecule(imol_map)) {
      int res_no = res_spec.res_no;
      std::string chain_id = res_spec.chain_id;
      std::string ins_code = res_spec.ins_code;
      mmdb::Residue *residue_p = get_residue(imol, res_spec);
      float f = molecules[imol].auto_fit_best_rotamer(rotamer_search_mode,
                                                      res_no, alt_conf, ins_code, chain_id,
                                                      imol_map, rotamer_fit_clash_flag,
                                                      rotamer_lowest_probability, *Geom_p());
      if (rotamer_auto_fit_do_post_refine_flag) {
         // Run refine zone with autoaccept, autorange on the "clicked" atom:
         // refine_auto_range(imol, chain_id.c_str(), res_no, alt_conf.c_str());
         // usee the refine() function
         std::cout << "ERROR:: auto_fit_rotamer_ng Missing refine() function"
                   << std::endl;
      }
      if (graphics_info_t::reset_b_factor_moved_atoms_flag) {
         // reset_b_factor_residue_range(imol, chain_id.c_str(), res_no, res_no);
         std::cout << "ERROR:: auto_fit_rotamer_ng Missing reset B-factor residue range function"
                   << std::endl;
      }
      update_geometry_graphs(&residue_p, 1, imol, imol_map); // yikes!
      std::cout << "Fitting score for best rotamer: " << f << std::endl;
      graphics_draw();
      run_post_manipulation_hook(imol, 0);
   } else {
      show_select_map_frame();
   }
}
