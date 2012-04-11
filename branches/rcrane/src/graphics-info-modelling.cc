/* src/graphics-info.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 by The University of York
 * Copyright 2007, 2008, 2009 by the University of Oxford
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


#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif

#include <algorithm>

#include <iostream>
#include <stdexcept>

#include <gtk/gtk.h>  // must come after mmdb_manager on MacOS X Darwin
#include <GL/glut.h>  // for some reason...  // Eh?

#include <gdk/gdkkeysyms.h> // for keyboarding (in this case nudge_active_residue) added 20091101

#include <sys/types.h> // for stating
#include <sys/stat.h>

#if !defined _MSC_VER
#include <unistd.h>
#endif

#include "guile-fixups.h"

#include "coot-sysdep.h"


#include <mmdb/mmdb_manager.h>
#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb-crystal.h"

#include "Cartesian.h"
#include "Bond_lines.h"
#ifdef USE_DUNBRACK_ROTAMERS
#include "dunbrack.hh"
#else 
#include "richardson-rotamer.hh"
#endif 

#include "clipper/core/map_utils.h" // Map_stats
#include "graphical_skel.h"

#include "interface.h"

#include "molecule-class-info.h"
#include "coot-coord-extras.hh"


#include "globjects.h"
#include "torsion-general.hh"
#include "ligand.hh"
#include "ideal-rna.hh"
#include "graphics-info.h"
#include "rotate-translate-modes.hh"

// Including python needs to come after graphics-info.h, because
// something in Python.h (2.4 - chihiro) is redefining FF1 (in
// ssm_superpose.h) to be 0x00004000 (Grrr).
//
// 20100813: Python.h needs to come before to stop"_POSIX_C_SOURCE" redefined problems 
//
// #ifdef USE_PYTHON
// #include "Python.h"
// #endif // USE_PYTHON


#include "coot-utils.hh"


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


//    std::cout << "DEBUG:: In copy_mol_and_refine() refine range: "
// 	     << "chain  :" << chain_id_1 << ": "
// 	     << resno_1 << " :" << inscode_1 << ": "
// 	     << resno_2 << " :" << inscode_2 << ": "
// 	     << "coords mol: " << imol_for_atoms << " map mol: " << imol_for_map
// 	     << std::endl;
      
#ifdef HAVE_GSL

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

   CMMDBManager *mol = molecules[imol].atom_sel.mol; // short-hand usage

   // We want to check for flanking atoms if the dictionary "group"
   // entry is not non-polymer.  So let's do a quick residue selection
   // of the first residue and find its residue type, look it up and
   // get the group.  If it is "non-polymer", then we can tinker with
   // the have_flanking_residue_at_* flags.

   int SelHnd_first = mol->NewSelection();
   int n_residue_first;
   PCResidue *residue_first = NULL;
   mol->Select(SelHnd_first, STYPE_RESIDUE, 0,
	       chain_id_1.c_str(),
	       resno_1, inscode_1.c_str(),
	       resno_1, inscode_1.c_str(),
	       "*",  // residue name
	       "*",  // Residue must contain this atom name?
	       "*",  // Residue must contain this Element?
	       "*",  // altLocs
	       SKEY_NEW); // selection key
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
	 geom_p->get_monomer_restraints(residue_type_first);
      if (p.first) {
	 group = p.second.residue_info.group;
      }
      cif_dictionary_read_number++;
   }
   mol->DeleteSelection(SelHnd_first);

   if (group != "non-polymer") { // i.e. it is (or can be) a polymer
      int SelHnd_ends = mol->NewSelection();
      int n_atoms_ends;
      PPCAtom atoms_end = 0;
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


   // Consider as the altconf the altconf of one of the residues (we
   // must test that the altlocs of the selected atoms to be either
   // the same as each other (A = A) or one of them is "".  We need to
   // know the mmdb syntax for "either".
   // 
   // 
   // 
   int iselection_resno_start = resno_1;
   int iselection_resno_end   = resno_2;
   if (have_flanking_residue_at_start) iselection_resno_start--;
   if (have_flanking_residue_at_end)   iselection_resno_end++;
   //
   int selHnd = mol->NewSelection();
   int nSelResidues; 
   PCResidue *SelResidues = NULL;
   mol->Select(selHnd, STYPE_RESIDUE, 0,
	       chain_id_1.c_str(),
	       iselection_resno_start, "*",
	       iselection_resno_end, "*",
	       "*",  // residue name
	       "*",  // Residue must contain this atom name?
	       "*",  // Residue must contain this Element?
	       "*",  // altLocs
	       SKEY_NEW // selection key
	       );
   molecules[imol].atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

   // Return 0 (first) if any of the residues don't have a dictionary
   // entry and a list of the residue type that don't have restraints.
   std::pair<int, std::vector<std::string> > icheck = 
      check_dictionary_for_residue_restraints(SelResidues, nSelResidues);


   if (0) {  // debugging.
      std::cout << "Selecting from chain id " << chain_id_1 << std::endl;
      std::cout << "selecting from residue " << iselection_resno_start
		<< " to " << iselection_resno_end << " selects "
		<< nSelResidues << " residues" << std::endl;
      std::cout << "=============== icheck: " << icheck.first << std::endl;
   }

   if (icheck.first == 0) { 

      std::cout << "INFO:: check_dictionary_for_residues - problem..." << std::endl;
      std::string problem_residues = "Refinement setup failure.\nFailed to find restraints for:\n";
      for (unsigned int icheck_res=0; icheck_res<icheck.second.size(); icheck_res++) { 
	 std::cout << "WARNING:: Failed to find restraints for :" 
		   << icheck.second[icheck_res] << ":" << std::endl;
	 problem_residues+= " ";
	 problem_residues+= icheck.second[icheck_res];
      }
      info_dialog(problem_residues);
      return rr; // fail
   } else {

      // 20100201
      
      bool check_hydrogens_too_flag = 0;
      // convert to CResidues vector
      std::vector<CResidue *> residues;
      for (int ires=0; ires<nSelResidues; ires++)
	 residues.push_back(SelResidues[ires]);
      std::pair<bool, std::vector<std::pair<std::string, std::vector<std::string> > > >
	 icheck_atoms = Geom_p()->atoms_match_dictionary(residues, check_hydrogens_too_flag);

      if (! icheck_atoms.first) {
	 std::cout << "Fail atom check" << std::endl;
	 info_dialog_refinement_non_matching_atoms(icheck_atoms.second);
	 return rr; // fail
      } else { 
	 
	 return copy_mol_and_refine_inner(imol_for_atoms,
					  resno_1, resno_2,
					  nSelResidues, SelResidues,
					  chain_id_1, altconf,
					  have_flanking_residue_at_start,
					  have_flanking_residue_at_end,
					  imol_for_map);
      }
   }
}

// static
void
graphics_info_t::info_dialog_missing_refinement_residues(const std::vector<std::string> &res_names) {

   std::string problem_residues = "Refinement setup failure.\nFailed to find restraints for:\n";
   for (unsigned int icheck_res=0; icheck_res<res_names.size(); icheck_res++) { 
      problem_residues+= " ";
      problem_residues+= res_names[icheck_res];
   }
   info_dialog(problem_residues);
}



coot::refinement_results_t
graphics_info_t::copy_mol_and_refine_inner(int imol_for_atoms,
					   int resno_1,
					   int resno_2,
					   int nSelResidues,
					   PCResidue *SelResidues,
					   const std::string &chain_id_1,
					   const std::string &altconf,
					   short int have_flanking_residue_at_start,
					   short int have_flanking_residue_at_end,
					   int imol_for_map
					   ) {

   coot::refinement_results_t rr(0, GSL_CONTINUE, "");
   short int have_disulfide_residues = 0; // of course not in linear mode.

   if (nSelResidues > 0) {

      std::vector<coot::atom_spec_t> fixed_atom_specs = molecules[imol_for_atoms].get_fixed_atoms();
      
      const char *chn = chain_id_1.c_str();

      if (nSelResidues > refine_regularize_max_residues) { 

	 std::cout << "WARNING:: Hit heuristic fencepost! Too many residues "
		   << "to refine\n "
		   << "          FYI: " << nSelResidues
		   << " > " << refine_regularize_max_residues
		   << " (which is your current maximum).\n";
	 std::cout << "Use (set-refine-max-residues "
		   << 2*refine_regularize_max_residues
		   << ") to increase limit\n";

      } else { 

	 // notice that we have to make 2 atom selections, one, which includes
	 // flanking (and disulphide) residues that is used for the restraints
	 // (restraints_container_t constructor) and one that is the moving atoms
	 // (which does not have flanking atoms).
	 // 
	 // The restraints_container_t moves the atom of the mol that is passes to
	 // it.  This must be the same mol as the moving atoms mol so that the
	 // changed atom positions can be seen.  However (as I said) the moving
	 // atom mol should not have flanking residues shown.  So we make an asc
	 // that has the same mol as that passed to the restraints but a different
	 // atom selection (it is the atom selection that is used in the bond
	 // generation).
	 //
	 short int in_alt_conf_split_flag = 0;
	 if (altconf != "")
	    in_alt_conf_split_flag = 1;
	    
	 CMMDBManager *residues_mol = 
	    create_mmdbmanager_from_res_selection(SelResidues, nSelResidues, 
						  have_flanking_residue_at_start,
						  have_flanking_residue_at_end,
						  altconf,
						  chain_id_1,
						  // 0, // 0 because we are not in alt conf split
						  in_alt_conf_split_flag, 
						  imol_for_atoms);

	 coot::restraints_container_t restraints(resno_1,
						 resno_2,
						 have_flanking_residue_at_start,
						 have_flanking_residue_at_end,
						 have_disulfide_residues,
						 altconf,
						 chn,
						 residues_mol,
						 fixed_atom_specs);

	 // this is where regularize and refine differ:
	 if (imol_for_map != -1)
	    restraints.add_map(molecules[imol_for_map].xmap_list[0],
			       geometry_vs_map_weight);

	 atom_selection_container_t local_moving_atoms_asc =
	    make_moving_atoms_asc(residues_mol, resno_1, resno_2);

	 // coot::restraint_usage_Flags flags = coot::BONDS;
	 // coot::restraint_usage_Flags flags = coot::BONDS_AND_ANGLES;
	 // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_AND_PLANES;
	 // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_TORSIONS_AND_PLANES; 
	 coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_AND_NON_BONDED;
	 // flags = coot::BONDS_ANGLES_PLANES_AND_NON_BONDED; 20071124
	 flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;

	 short int do_residue_internal_torsions = 0;

	 if (do_torsion_restraints) { 
	    do_residue_internal_torsions = 1;
	    // flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS; // fail
	    flags = coot::BONDS_ANGLES_AND_TORSIONS; // OK
	    flags = coot::BONDS_ANGLES_TORSIONS_AND_PLANES;  // OK
	    // flags = coot::BONDS_ANGLES_TORSIONS_PLANES_AND_NON_BONDED; // fail
	    flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS; // fail
	 } 

	 if (do_rama_restraints) 
	    flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_RAMA;

	 // coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;

	 // 20080108 Recall that we do secondary structure restraints
	 // with pseudo bonds now.  We don't do it by torsion
	 // refinement of the phi and psi.
	 //
	 // However, ramachandran goodness will use phi and psi
	 //
	 if (molecules[imol_for_atoms].extra_restraints.has_restraints())
	    restraints.add_extra_restraints(molecules[imol_for_atoms].extra_restraints);
	 int nrestraints = 
	    restraints.make_restraints(*geom_p, flags,
				       do_residue_internal_torsions,
				       rama_plot_restraint_weight,
				       do_rama_restraints,
				       pseudo_bonds_type);

	 if (do_numerical_gradients)
	    restraints.set_do_numerical_gradients();

	 rr = update_refinement_atoms(nrestraints, restraints, rr, local_moving_atoms_asc,
				      1, imol_for_atoms, chain_id_1);

	 // local_moving_atoms_asc.clear_up(); // crash.
	 // 
	 // Hmm... local_moving_atoms_asc gets transfered to (class
	 // variable) moving_atoms_asc in update_refinement_atoms().
	 // Do we delete old moving_atoms_asc mol and selection there
	 // before they are replaced?
	 //
	 // No.
	 //
	 // Should we?
	 //
	 // Yes.
	 //
	 // Problem is that we are not sure that every time the old
	 // moving_atoms_asc is replaced in that manner
	 // (make_moving_atoms_asc()) that moving_atoms_asc is
	 // expired/unreferenced.
	 
	 
      }
      
   } else {
      std::cout << "No Atoms!!!!  This should never happen: " << std::endl;
      std::cout << "  in create_regularized_graphical_object" << std::endl;
   } 
   return rr;
#else 

   std::cout << "Cannot refine without compilation with GSL" << std::endl;
   return coot::refinement_results_t(0, 0, "");

#endif // HAVE_GSL
}


coot::refinement_results_t
graphics_info_t::update_refinement_atoms(int n_restraints,
					 coot::restraints_container_t &restraints,
					 coot::refinement_results_t rr_in,
					 atom_selection_container_t local_moving_atoms_asc,
					 bool need_residue_order_check,
					 int imol,
					 std::string chain_id_1) {

   coot::refinement_results_t rr = rr_in;
   
   if (n_restraints > 0) {
      moving_atoms_asc_type = coot::NEW_COORDS_REPLACE;
      last_restraints = restraints;

      regularize_object_bonds_box.clear_up();
      make_moving_atoms_graphics_object(local_moving_atoms_asc); // sets moving_atoms_asc
      int n_cis = coot::util::count_cis_peptides(moving_atoms_asc->mol);
      graphics_info_t::moving_atoms_n_cis_peptides = n_cis;
      // std::cout << "DEBUG:: start of ref have: " << n_cis << " cis peptides"
      // << std::endl;
      short int continue_flag = 1;
      int step_count = 0; 
      print_initial_chi_squareds_flag = 1; // unset by drag_refine_idle_function
      while ((step_count < 10000) && continue_flag) {
	 int retval = drag_refine_idle_function(NULL);
	 step_count += dragged_refinement_steps_per_frame;
	 if (retval == GSL_SUCCESS) { 
	    continue_flag = 0;
	    rr = graphics_info_t::saved_dragged_refinement_results;
	 }
	 if (retval == GSL_ENOPROG) {
	    continue_flag = 0;
	    rr = graphics_info_t::saved_dragged_refinement_results;
	 }
      }

      // if we reach here with continue_flag == 1, then we
      // were refining (regularizing more like) and ran out
      // of steps before convergence.  We still want to give
      // the use a dialog though.
      //
      if (continue_flag == 1) {
	 rr = graphics_info_t::saved_dragged_refinement_results;
	 rr.info = "Time's up...";
      }
      
   } else { 
      if (use_graphics_interface_flag) {

	 // Residues might be out of order (causing problems in atom
	 // selection).  But the sphere selection doesn't need this
	 // check.

	 GtkWidget *widget = create_no_restraints_info_dialog();
	 if (! need_residue_order_check) {
	    gtk_widget_show(widget);
	 } else { 
	    bool residues_ordered_flag =
	       molecules[imol].progressive_residues_in_chain_check_by_chain(chain_id_1.c_str());
	    
	    GtkWidget *l = lookup_widget(widget, "no_restraints_extra_label");
	    if (l) {
	       if (!residues_ordered_flag) {
		  gtk_widget_show(l);
	       } else {
		  gtk_widget_hide(l); // by default it is show?
	       }
	    }
	    gtk_widget_show(widget);
	 }
      }
   } 
   return rr;
} 



coot::refinement_results_t 
graphics_info_t::refine_residues_vec(int imol, 
				     const std::vector<CResidue *> &residues,
				     const char *alt_conf, 
				     CMMDBManager *mol) {

   bool use_map_flag = 1;
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
graphics_info_t::regularize_residues_vec(int imol, 
					 const std::vector<CResidue *> &residues,
					 const char *alt_conf, 
					 CMMDBManager *mol) {

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

// simple CResidue * interface to refinement.  20081216
//
// Needs use_map flag, I guess
// 
coot::refinement_results_t
graphics_info_t::generate_molecule_and_refine(int imol,
					      const std::vector<CResidue *> &residues,
					      const char *alt_conf,
					      CMMDBManager *mol,
					      bool use_map_flag) { 

   coot::refinement_results_t rr(0, GSL_CONTINUE, "");
   
#ifdef HAVE_GSL

   if (is_valid_map_molecule(Imol_Refinement_Map()) || (! use_map_flag)) {
      float weight = geometry_vs_map_weight;
      coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
      short int do_residue_internal_torsions = 0;
      if (do_torsion_restraints) { 
	 do_residue_internal_torsions = 1;
	 flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS;
      } 
      
      if (do_rama_restraints) 
	 flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_RAMA;
      
      std::vector<coot::atom_spec_t> fixed_atom_specs = molecules[imol].get_fixed_atoms();

      // OK, so the passed residues are the residues in the graphics_info_t::molecules[imol]
      // molecule.  We need to do 2 things:
      //
      // convert the CResidue *s of the passed residues to the CResidue *s of residues mol
      //
      // and
      //
      // in create_mmdbmanager_from_res_vector() make sure that that contains the flanking atoms.
      //
      // The flanking atoms are fixed the passed residues are not fixed.
      // Keep a clear head.
      
      std::vector<std::string> residue_types = coot::util::residue_types_in_residue_vec(residues);
      // use try_dynamic_add()
      bool have_restraints = geom_p->have_dictionary_for_residue_types(residue_types);

      if (have_restraints) { 
      
	 std::string residues_alt_conf = alt_conf;
	 imol_moving_atoms = imol;
	 std::pair<CMMDBManager *, std::vector<CResidue *> > residues_mol_and_res_vec =
	    create_mmdbmanager_from_res_vector(residues, imol, mol, residues_alt_conf);

	 // We only want to act on these new residues and molecule, if
	 // there is something there.
	 // 
	 if (residues_mol_and_res_vec.first) {

	    // Now we want to do an atom name check.  This stops exploding residues.
	    //
	    bool check_hydrogens_too_flag = 0;
	    std::pair<bool, std::vector<std::pair<std::string, std::vector<std::string> > > >
	       icheck_atoms = Geom_p()->atoms_match_dictionary(residues, check_hydrogens_too_flag); 
	    
	    if (! icheck_atoms.first) {
	       // Oops. Just give us a dialog and don't start the refinement
	       info_dialog_refinement_non_matching_atoms(icheck_atoms.second);
	       
	    } else { 
	    
	       atom_selection_container_t local_moving_atoms_asc =
		  make_moving_atoms_asc(residues_mol_and_res_vec.first, residues);
	       std::vector<std::pair<bool,CResidue *> > local_residues;  // not fixed.
	       for (unsigned int i=0; i<residues_mol_and_res_vec.second.size(); i++)
		  local_residues.push_back(std::pair<bool, CResidue *>(0, residues_mol_and_res_vec.second[i]));
	       coot::restraints_container_t restraints(local_residues, *Geom_p(),
						       residues_mol_and_res_vec.first,
						       fixed_atom_specs);

	       if (use_map_flag) { 
		  clipper::Xmap<float> &xmap = molecules[Imol_Refinement_Map()].xmap_list[0];
		  restraints.add_map(xmap, weight);
	       }
	 
	       if (molecules[imol].extra_restraints.has_restraints())
		  restraints.add_extra_restraints(molecules[imol].extra_restraints);
	       int n_restraints = restraints.make_restraints(*Geom_p(), flags,
							     do_residue_internal_torsions,
							     rama_plot_restraint_weight,
							     do_rama_restraints,
							     pseudo_bonds_type);
	       
	 
	       std::string dummy_chain = ""; // not used
	       rr = update_refinement_atoms(n_restraints, restraints, rr, local_moving_atoms_asc,
					    0, imol, dummy_chain);
	    }
	 }
      } else {
	 // we didn't have restraints for everything.
	 // 
	 std::pair<int, std::vector<std::string> > icheck = check_dictionary_for_residue_restraints(residues);
	 if (icheck.first == 0) { 
	    info_dialog_missing_refinement_residues(icheck.second);
	 }

      } 
   }

   return rr;

#else
   
   std::cout << "Cannot refine without compilation with GSL" << std::endl;
   return coot::refinement_results_t(0, 0, "");

#endif    
}

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
graphics_info_t::make_moving_atoms_asc(CMMDBManager *residues_mol,
				       int resno_1,
				       int resno_2) const {

   atom_selection_container_t local_moving_atoms_asc;
   local_moving_atoms_asc.mol = (MyCMMDBManager *) residues_mol;
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

   return local_moving_atoms_asc;
}


atom_selection_container_t
graphics_info_t::make_moving_atoms_asc(CMMDBManager *residues_mol,
				       const std::vector<CResidue *> &residues) const {

   atom_selection_container_t local_moving_atoms_asc;
   local_moving_atoms_asc.UDDOldAtomIndexHandle = -1;  // true?
   local_moving_atoms_asc.UDDAtomIndexHandle = -1;

   int SelHnd = residues_mol->NewSelection();

   for (unsigned int ir=0; ir<residues.size(); ir++) {
      const char *chain_id = residues[ir]->GetChainID();
      const char *inscode = residues[ir]->GetInsCode();
      int resno = residues[ir]->GetSeqNum();
      residues_mol->Select(SelHnd, STYPE_ATOM,
			   0, chain_id,
			   resno, // starting resno, an int
			   inscode, // any insertion code
			   resno, // ending resno
			   inscode, // ending insertion code
			   "*", // any residue name
			   "*", // atom name
			   "*", // elements
			   "*",  // alt loc.	
			   SKEY_OR);
   }

   local_moving_atoms_asc.mol = (MyCMMDBManager *) residues_mol;
   local_moving_atoms_asc.SelectionHandle = SelHnd;
   residues_mol->GetSelIndex(local_moving_atoms_asc.SelectionHandle,
			     local_moving_atoms_asc.atom_selection,
			     local_moving_atoms_asc.n_selected_atoms);
   // std::cout << "returning a atom selection for all moving atoms "
   // << local_moving_atoms_asc.n_selected_atoms << " atoms "
   // << std::endl;
   return local_moving_atoms_asc;
} 


// Return 0 (first) if any of the residues don't have a dictionary
// entry and a list of the residue type that don't have restraints.
// 
std::pair<int, std::vector<std::string> >
graphics_info_t::check_dictionary_for_residue_restraints(PCResidue *SelResidues, int nSelResidues) {

   int status;
   bool status_OK = 1; // pass, by default
   std::vector<std::string> res_name_vec;

   for (int ires=0; ires<nSelResidues; ires++) {
      std::string resn(SelResidues[ires]->GetResName());
      std::string resname = adjust_refinement_residue_name(resn);
      int status = geom_p->have_dictionary_for_residue_type(resname, cif_dictionary_read_number);
      if (! status) { 
	 status_OK = 0;
	 res_name_vec.push_back(resname);
      } 

      if (0)
	 std::cout << "DEBUG:: have_dictionary_for_residues() on residue "
		   << ires << " of " << nSelResidues << ", "
		   << resname << " returns "
		   << status << std::endl;
      cif_dictionary_read_number++;
   }
   return std::pair<int, std::vector<std::string> > (status_OK, res_name_vec);
}

std::pair<int, std::vector<std::string> >
graphics_info_t::check_dictionary_for_residue_restraints(const std::vector<CResidue *> &residues) {

   std::vector<std::string> res_name_vec;
   std::pair<int, std::vector<std::string> > r(0, res_name_vec);
   for (unsigned int i=0; i<residues.size(); i++) {
      std::string resname = adjust_refinement_residue_name(residues[i]->GetResName());
      int status = geom_p->have_dictionary_for_residue_type(resname, cif_dictionary_read_number);
      if (! status) {
	 r.first = 0;
	 r.second.push_back(resname);
      }
      cif_dictionary_read_number++; // not sure why this is needed.
   }
   return r;
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
// 
CMMDBManager *
graphics_info_t::create_mmdbmanager_from_res_selection(PCResidue *SelResidues, 
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

   CMMDBManager *residues_mol = new CMMDBManager;
   CModel *model = new CModel;
   CChain *chain = new CChain;
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

   CResidue *r;
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
      r = coot::deep_copy_this_residue(SelResidues[ires], altconf, whole_res_flag, 
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
   residues_mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
   residues_mol->FinishStructEdit();

   return residues_mol;
}

// called by simple_refine_residues (a refinement from a vector of CResidues).
//
// The returned mol should have flanking residues too.
// 
// return a NULL in the first of the pair if the past residue vector is of size 0.
// 
std::pair<CMMDBManager *, std::vector<CResidue *> >
graphics_info_t::create_mmdbmanager_from_res_vector(const std::vector<CResidue *> &residues,
						    int imol, 
						    CMMDBManager *mol_in,
						    std::string alt_conf) {
   
   CMMDBManager *new_mol = 0;
   std::vector<CResidue *> rv;
   int n_flanker = 0; // a info/debugging counter

   if (residues.size() > 0) { 

      new CMMDBManager;
      CModel *model_p = new CModel;
      CChain *chain_p = new CChain;
   
      float dist_crit = 3.0;

      // First add new versions of the passed residues:
      // 
      std::string chain_id_1 = residues[0]->GetChainID();
      short int whole_res_flag = 0;
      int atom_index_udd = molecules[imol].atom_sel.UDDAtomIndexHandle;

      // FIXME, we need to check that we are adding r to the right
      // chain, now residues can be in many different chains.
      //
      
      for (unsigned int ires=0; ires<residues.size(); ires++) {
	 CResidue *r;

	 std::string ref_res_chain_id = residues[ires]->GetChainID();

	 CChain *chain_p = NULL;
	 int n_new_mol_chains = model_p->GetNumberOfChains();
	 for (int ich=0; ich<n_new_mol_chains; ich++) {
	    if (ref_res_chain_id == model_p->GetChain(ich)->GetChainID()) {
	       chain_p = model_p->GetChain(ich);
	       break;
	    }
	 }

	 // Add a new one then.
	 if (! chain_p) {
	    chain_p = new CChain;
	    chain_p->SetChainID(ref_res_chain_id.c_str());
	    model_p->AddChain(chain_p);
	 }

	 // found in mmdb-extras
	 r = coot::deep_copy_this_residue(residues[ires], alt_conf, whole_res_flag, 
					  atom_index_udd);
	 if (r) { 
	    chain_p->AddResidue(r);
	    r->seqNum = residues[ires]->GetSeqNum();
	    r->SetResName(residues[ires]->GetResName());
	    // 	 std::cout << " adding moving residue " << " " << coot::residue_spec_t(r)
	    // 		   << std::endl;
	    rv.push_back(r);
	 } 
      }

      new_mol = new CMMDBManager;
      new_mol->AddModel(model_p);
      new_mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
      new_mol->FinishStructEdit();

      
      // Now the flanking residues:
      //
      std::vector<CResidue *> flankers_in_reference_mol;
      
      for (unsigned int ires=0; ires<residues.size(); ires++) {
	 CResidue *res_ref = residues[ires];
	 std::vector<CResidue *> neighbours =
	    coot::residues_near_residue(res_ref, mol_in, dist_crit);
	 // now add the elements of neighbours if they are not already
	 // in flankers_in_reference_mol (and not in residues either of
	 // course)
	 // 
	 for (unsigned int in=0; in<neighbours.size(); in++) { 
	    bool found = 0;
	    for (unsigned int iflank=0; iflank<flankers_in_reference_mol.size(); iflank++) {
	       if (neighbours[in] == flankers_in_reference_mol[iflank]) {
		  found = 1;
		  break;
	       }
	    }

	    // in residues?
	    if (! found) {
	       for (unsigned int ii=0; ii<residues.size(); ii++) {
		  if (neighbours[in] == residues[ii]) {
		     found = 1;
		     break;
		  } 
	       } 
	    } 
	    
	    if (! found) {
	       flankers_in_reference_mol.push_back(neighbours[in]);
	    }
	 } 
      }

      // So we have a vector of residues that were flankers in the
      // reference molecule, we need to add copies of those to
      // new_mol (making sure that they go into the correct chain).
      for (unsigned int ires=0; ires<flankers_in_reference_mol.size(); ires++) {
	 CResidue *r;

	 std::string ref_res_chain_id = flankers_in_reference_mol[ires]->GetChainID();

	 CChain *chain_p = NULL;
	 int n_new_mol_chains = model_p->GetNumberOfChains();
	 for (int ich=0; ich<n_new_mol_chains; ich++) {
	    if (ref_res_chain_id == model_p->GetChain(ich)->GetChainID()) {
	       chain_p = model_p->GetChain(ich);
	       break;
	    }
	 }

	 if (! chain_p) {
	    // Add a new one then.
	    chain_p = new CChain;
	    chain_p->SetChainID(ref_res_chain_id.c_str());
	    model_p->AddChain(chain_p);
	 }

	 r = coot::deep_copy_this_residue(flankers_in_reference_mol[ires],
					  alt_conf, whole_res_flag, 
					  atom_index_udd);
	 if (r) { 
	    chain_p->AddResidue(r);
	    r->seqNum = flankers_in_reference_mol[ires]->GetSeqNum();
	    r->SetResName(flankers_in_reference_mol[ires]->GetResName());
	    // 	 std::cout << " adding flanking residue " << " " << coot::residue_spec_t(r)
	    // 		   << std::endl;
	    n_flanker++;
	 }
      }
   }
//    std::cout << "DEBUG:: in create_mmdbmanager_from_res_vector: " << rv.size()
// 	     << " free residues and " << n_flanker << " flankers" << std::endl;
   return std::pair <CMMDBManager *, std::vector<CResidue *> > (new_mol, rv);
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
      short int ifound = 0;
      for (int igeom=0; igeom<geom_p->size(); igeom++) {
	 if (pdb_residue_types[ipdb] == (*geom_p)[igeom].comp_id()) { 
	    ifound = 1;
	    break;
	 }
      }
      if (ifound == 0) { 
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

   PPCAtom SelAtom = molecules[imol].atom_sel.atom_selection; 

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
      cout << "Picked atoms are not in the same chain.  Failure" << endl; 
   } else { 
      flash_selection(imol, resno_1, inscode_1, resno_2, inscode_2, altconf, chain_id_1);
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
   
   CAtom *this_atom =  molecules[imol].atom_sel.atom_selection[atom_index];
   CResidue *this_res = this_atom->residue;
   CChain *this_chain = this_res->chain;
   int resno = this_res->GetSeqNum();
   char *inscode = this_res->GetInsCode();
   
   CResidue *prev_res = this_chain->GetResidue(resno-refine_auto_range_step, inscode);
   CResidue *next_res = this_chain->GetResidue(resno+refine_auto_range_step, inscode);

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
graphics_info_t::refinement_results_to_scm(coot::refinement_results_t &rr) {

   SCM r = SCM_BOOL_F;

   if (rr.found_restraints_flag) {
      SCM lights_scm = SCM_EOL;
      SCM progress_scm = SCM_MAKINUM(rr.progress);
      SCM info_scm = scm_from_locale_string(rr.info.c_str());
      for (int il=rr.lights.size()-1; il>=0; il--) {
	 SCM light_scm = SCM_EOL;
	 SCM value_scm = scm_double2num(rr.lights[il].value);
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
graphics_info_t::refinement_results_to_py(coot::refinement_results_t &rr) {
   PyObject *r = Py_False;

   if (rr.found_restraints_flag) {
      PyObject *lights_py = Py_False;
      PyObject *progress_py = PyInt_FromLong(rr.progress);
      PyObject *info_py = PyString_FromString(rr.info.c_str());
      if (rr.lights.size())
	lights_py = PyList_New(rr.lights.size());
      for (int il=0; il<rr.lights.size(); il++) {
	PyObject *light_py = PyList_New(3);
	PyObject *value_py = PyFloat_FromDouble(rr.lights[il].value);
	PyObject *label_py = PyString_FromString(rr.lights[il].label.c_str());
	PyObject *name_py  = PyString_FromString(rr.lights[il].name.c_str());
	
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



void
graphics_info_t::flash_selection(int imol,
				 int resno_1, 
				 std::string ins_code_1,
				 int resno_2, 
				 std::string ins_code_2,
				 std::string altconf,
				 std::string chain_id_1) { 

   // First make an atom selection of the residues selected to regularize.
   // 
   int selHnd = ((CMMDBManager *)molecules[imol].atom_sel.mol)->NewSelection();
   int nSelAtoms;
   PPCAtom SelAtom;
   const char *chn  = chain_id_1.c_str();
   const char *ins1 = ins_code_1.c_str();
   const char *ins2 = ins_code_2.c_str();

   ((CMMDBManager *)molecules[imol].atom_sel.mol)->SelectAtoms(selHnd, 0, 
							       chn,
							       resno_1, ins1,
							       resno_2, ins2,
							       "*",      // RNames
							       "*","*",  // ANames, Elements
							       "*" );    // Alternate locations.


   ((CMMDBManager *)molecules[imol].atom_sel.mol)->GetSelIndex(selHnd,
							       SelAtom,
							       nSelAtoms);
//    cout << nSelAtoms << " atoms selected to regularize from residue "
// 	<< resno_1 << " to " << resno_2 << " chain " << chn << endl;

   if (nSelAtoms) { 
      // now we can make an atom_selection_container_t with our new
      // atom selection that we will use to find bonds.

      atom_selection_container_t asc; 
      asc.mol = molecules[imol].atom_sel.mol; 
      asc.atom_selection = SelAtom; 
      asc.n_selected_atoms = nSelAtoms; 

      int fld = 0;
      Bond_lines_container bonds(asc, fld); // don't flash disulfides

      graphical_bonds_container empty_box; 
      graphical_bonds_container regular_box = bonds.make_graphical_bonds();
      
      int flash_length = residue_selection_flash_frames_number;

      if (glarea) { 
	 for (int iflash=0; iflash<2; iflash++) { 
	    regularize_object_bonds_box = regular_box; 
	    for (int i=0; i<flash_length; i++)
	       graphics_draw();
	    regularize_object_bonds_box = empty_box; 
	    for (int i=0; i<flash_length; i++)
	       graphics_draw();
	 }
      }
      regularize_object_bonds_box = empty_box; 

   } // atoms selected
   molecules[imol].atom_sel.mol->DeleteSelection(selHnd);
   graphics_draw();
}

// static
void
graphics_info_t::flash_position(const clipper::Coord_orth &pos) {

   if (glarea) {
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

//    cout << "refine (fit to map): molecule " << imol
// 	<< " atom index " << i_atom_no_1
// 	<< " to " << i_atom_no_2 << endl; 

   int resno_1, resno_2;

   int imol_map = Imol_Refinement_Map();
   if (imol_map == -1) { // magic number check,
      // if not -1, then it has been set by user

      show_select_map_dialog();

   } else { 

      PPCAtom SelAtom = molecules[imol].atom_sel.atom_selection; 

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
graphics_info_t::check_for_single_hetatom(CResidue *res_p) const {

   bool r = 0;

   int n_atoms = res_p->GetNumberOfAtoms();
   if (n_atoms == 1) {
      PPCAtom residue_atoms;
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
// chain_id_2 are not the same chain (CChain *), then we don't call
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

//    std::cout << "DEBUG:: ================ refine_residue_range: "
// 	     << imol << " " << chain_id_1
//  	     << " " <<  resno_1 << ":" << ins_code_1 << ":"
//  	     << " " <<  resno_2 << ":" << ins_code_2 << ":"
//  	     << " " << ":" << altconf << ": " << is_water_like_flag << std::endl;

   coot::refinement_results_t rr;
   
   int imol_map = Imol_Refinement_Map();
   if (imol_map == -1) { // magic number check,
      // if not -1, then it has been set by user
      show_select_map_dialog();
   } else { 

      // if ( chain_id_1 != chain_id_2 ) {
      // Used to be pointer comparison, let that be done in the calling function.
      if (chain_id_1 != chain_id_2) {

	 // for now we will bug out.  In futre, we will want to be
	 // able to refine glycosylation.
	 //
	 std::cout << "Picked atoms are not in the same chain.  Failure" << std::endl;
	 std::cout << "FYI: chain ids are: \"" << chain_id_1
		   << "\" and \"" << chain_id_2 << "\"" << std::endl;
      } else {
	 if (molecules[imol_map].has_map()) {  // it may have been
					       // closed after it was
					       // selected.
	    short int simple_water = 0;
	    if (resno_1 == resno_2) {
	       if (is_water_like_flag) {
		  simple_water = 1;
// 		  std::string s = "That water molecule doesn't have restraints.\n";
// 		  s += "Using Rigid Body Fit Zone";
// 		  GtkWidget *w = wrapped_nothing_bad_dialog(s);
// 		  gtk_widget_show(w);
		  
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
	       flash_selection(imol, resno_1, ins_code_1, resno_2, ins_code_2, altconf, chain_id_1);
	       long t0 = glutGet(GLUT_ELAPSED_TIME);
	       rr = copy_mol_and_refine(imol, imol_map, resno_1, ins_code_1, resno_2, ins_code_2,
					altconf, chain_id_1);
	       short int istat = rr.found_restraints_flag;
	       long t1 = glutGet(GLUT_ELAPSED_TIME);
	       std::cout << "Refinement elapsed time: " << float(t1-t0)/1000.0 << std::endl;
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
	    show_select_map_dialog();
	 } 
      }
   } // same chains test
   return rr;
}


// Question to self: Are you sure that imol_rigid_body_refine (the
// coordinates molecule) is set when we get here?
// 
// Also: are residue_range_atom_index_1 and residue_range_atom_index_2 set?
// They should be.
// 
void
graphics_info_t::execute_rigid_body_refine(short int auto_range_flag) { /* atom picking has happened.
						Actually do it */

   CAtom *atom1;
   CAtom *atom2;

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
   
   std::string chain(chain_id_1);
   std::string altconf = atom1->altLoc;

//    std::cout << "-----------------------------------------------------" << std::endl;
//    std::cout << "-----------------------------------------------------" << std::endl;
//    std::cout << " Rigid Body Refinement "
// 	     << " imol: " << imol_rigid_body_refine << " residue "
// 	     << ires1 << " to " << ires2 << " chain " << chain << std::endl;

   int imol_ref_map = Imol_Refinement_Map();  // -1 is a magic number

   if (Imol_Refinement_Map() == -1 ) { // magic number
      //
      std::cout << "Please set a map against which the refimentment should occur"
		<< std::endl;
      show_select_map_dialog();  // protected
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
		  std::vector<int> from_ref_delete_atom_indices;

		  // a vector of atoms that should be deleted from the
		  // moving mol, because they don't match the altconf
		  std::vector<int> from_mov_delete_atom_indices;

		  try { // (yes, belt and braces in this case)
		     if (! mol[ifrag][ires].is_undefined()) { 
			// std::cout << "adding residue mol[" << ifrag<< "][" << ires << "]: "
			// << mol[ifrag][ires] << std::endl;
			range_mol[ir].addresidue(mol[ifrag][ires], 1);
		     }
		  }
		  catch (std::runtime_error rte) {
		     std::cout << "ERROR:: execute_rigid_body_refine() " << rte.what() << std::endl;
		  } 
		  

		  for (unsigned int iat=0; iat<mol[ifrag][ires].atoms.size(); iat++) {
		     if (mol[ifrag][ires][iat].altLoc == altconf) {
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

   bool success = 0; // fail initially
   bool debug = 0;
   
   std::vector<coot::minimol::atom *> range_atoms = range_mol.select_atoms_serial();
   //       std::cout << "There are " << range_atoms.size() << " atoms from initial ligand "
   // 		<< std::endl;


   // debugging
   if (debug) {
      range_mol.write_file("rigid-body-range-mol.pdb", 44);
      mol_without_moving_zone.write_file("rigid-body-without-moving-zone.pdb", 44);
   } 
   
   coot::ligand lig;
   lig.import_map_from(molecules[imol_ref_map].xmap_list[0], 
		       molecules[imol_ref_map].map_sigma());
   
   lig.install_ligand(range_mol);
   lig.find_centre_by_ligand(0); // don't test ligand size
   lig.set_map_atom_mask_radius(0.5);
   lig.mask_map(mol_without_moving_zone, mask_water_flag);
   if (debug)
      lig.output_map("rigid-body.map");
   lig.set_dont_write_solutions();
   lig.set_dont_test_rotations();
   lig.set_acceptable_fit_fraction(graphics_info_t::rigid_body_fit_acceptable_fit_fraction);
   lig.fit_ligands_to_clusters(1);
   coot::minimol::molecule moved_mol = lig.get_solution(0);
   
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
      // 	 rigid_body_asc.mol = (MyCMMDBManager *) moved_mol.pcmmdbmanager();

      // 	 int SelHnd = rigid_body_asc.mol->NewSelection();
      // 	 rigid_body_asc.mol->SelectAtoms(SelHnd, 0, "*",
      // 					 ANY_RES, // starting resno, an int
      // 					 "*", // any insertion code
      // 					 ANY_RES, // ending resno
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
      rigid_body_asc = make_asc(moved_mol.pcmmdbmanager());

      moving_atoms_asc_type = coot::NEW_COORDS_REPLACE;
      imol_moving_atoms = imol_rigid_body_refine;
      make_moving_atoms_graphics_object(rigid_body_asc);
      // 	 std::cout << "DEBUG:: execute_rigid_body_refine " 
      // 		   << " make_moving_atoms_graphics_object UDDOldAtomIndexHandle " 
      // 		   << moving_atoms_asc->UDDOldAtomIndexHandle << std::endl;
      graphics_draw();
      if (! refinement_immediate_replacement_flag) { 
	 coot::refinement_results_t dummy;
	 if (use_graphics_interface_flag) { 
	    do_accept_reject_dialog("Rigid Body Fit", dummy); // constructed ref res
	 }
      }
      // 
   } else {
      if (use_graphics_interface_flag) { 
	 GtkWidget *w = create_rigid_body_refinement_failed_dialog();
	 gtk_widget_show(w);
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
	 PPCAtom selatoms;
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
	    if (selatoms[0]->GetUDData(molecules[imol].atom_sel.UDDAtomIndexHandle, ind_1) == UDDATA_Ok) {
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
	    if (selatoms[0]->GetUDData(molecules[imol].atom_sel.UDDAtomIndexHandle, ind_2) == UDDATA_Ok) {
	       residue_range_atom_index_2 = ind_2;
	    }
	 }
	 molecules[imol].atom_sel.mol->DeleteSelection(SelHnd);

	 //   if (ind_1 >= 0 && ind_2 >= 0)
	 //     execute_rigid_body_refine(0); // not autorange
      }
   }
}


#include "residue_by_phi_psi.hh"

// The passed residue type is either N, C or (now [20031222]) M.
// 
void 
graphics_info_t::execute_add_terminal_residue(int imol, 
					      const std::string &terminus_type,
					      const CResidue *res_p,
					      const std::string &chain_id, 
					      const std::string &res_type_in,
					      short int immediate_addition_flag) {

   std::string res_type = res_type_in; // const
   int imol_map = Imol_Refinement_Map();
   if (imol_map == -1) { 
      if (1) { 
	 show_select_map_dialog();
      } else { 
	 // just shove it on without a map
	 if (molecules[imol].has_model()) {
	    // float phi = graphics_info_t::terminal_residue_addition_direct_phi;
	    // float psi = graphics_info_t::terminal_residue_addition_direct_psi;
	    CMMDBManager *orig_mol = graphics_info_t::molecules[imol].atom_sel.mol;
	    //	    CResidue *res_new = add_terminal_residue_directly(terminus_type, res_p,
	    // chain_id, res_type, phi, psi);
	    CResidue *res_new = 0;
	    CMMDBManager *new_mol = coot::util::create_mmdbmanager_from_residue(orig_mol,
										res_new);
	    if (new_mol) { 
	       atom_selection_container_t extra_residue_asc = make_asc(new_mol);
	       graphics_info_t::molecules[imol].add_coords(extra_residue_asc);
	    }
	 }
      }
   } else { 


      if (terminus_type == "not-terminal-residue") {
	 std::cout << "that residue was not at a terminus" << std::endl;
      } else {

	 imol_moving_atoms = imol;

	 std::string residue_type_string = res_type;
	 CResidue *unconst_res_p = (CResidue *) res_p;     // bleugh.
	 int residue_number = unconst_res_p->GetSeqNum(); // bleugh.
	 if (residue_type_string == "auto") {
	    int resno_added = -1; // was unset
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

	 float bf = default_new_atoms_b_factor;
	 coot::residue_by_phi_psi addres(molecules[imol].atom_sel.mol,
					 terminus_type, res_p, chain_id, res_type, bf);

	 // std::cout << "DEBUG:: term_type: " << terminus_type << std::endl;

	 // This was for debugging, so that we get *some* solutions at least.
	 // 	
	 addres.set_acceptable_fit_fraction(0.0); //  the default is 0.5, I think

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
	 if (terminus_type == "MC" || terminus_type == "MN" ||
	     terminus_type == "singleton")
	    masked_map_val = 0.0; 
	 addres.set_masked_map_value(masked_map_val);   
	 addres.import_map_from(molecules[imol_map].xmap_list[0], 
				molecules[imol_map].map_sigma());

	 // This masked map will be the one that is used for rigid
	 // body refinement, unlike normal ligand class usage which
	 // uses the xmap_pristine.
	 // 

	 // old style mask by all the protein atoms - slow? (especially on sgis?)
	 // 	 short int mask_waters_flag = 1; 
	 //      addres.mask_map(molecules[imol].atom_sel.mol, mask_waters_flag);
	 //
	 PPCAtom atom_sel = NULL;
	 int n_selected_atoms = 0;
	 realtype radius = 8.0;  // more than enough for 2 residue mainchains.
	 int SelHndSphere = molecules[imol].atom_sel.mol->NewSelection();
	 CAtom *terminal_at = NULL;
	 std::string atom_name = "Unassigned";
	 if (terminus_type == "MC" || terminus_type == "C" ||
	     terminus_type == "singleton")
	    atom_name = " C  ";
	 if (terminus_type == "MN" || terminus_type == "N")
	    atom_name = " N  ";
	 if (atom_name != "Unassigned") { 
	    PPCAtom residue_atoms;
	    int nResidueAtoms;
	    CResidue *res_tmp_p = (CResidue *) res_p;
	    res_tmp_p->GetAtomTable(residue_atoms, nResidueAtoms);
	    for (int i=0; i<nResidueAtoms; i++)
	       if (atom_name == residue_atoms[i]->name) {
		  terminal_at = residue_atoms[i];
		  break;
	       }

	    if (terminal_at) { 
	       molecules[imol].atom_sel.mol->SelectSphere(SelHndSphere, STYPE_ATOM,
							  terminal_at->x,
							  terminal_at->y,
							  terminal_at->z,
							  radius, SKEY_NEW);
	       molecules[imol].atom_sel.mol->GetSelIndex(SelHndSphere, atom_sel, n_selected_atoms);
	       int invert_flag = 0;
	       addres.mask_map(molecules[imol].atom_sel.mol, SelHndSphere, invert_flag);
	       molecules[imol].atom_sel.mol->DeleteSelection(SelHndSphere);
	    }
	 } else {
	    std::cout << "WARNING:: terminal atom not assigned - no masking!" << std::endl;
	 }

 	 // addres.output_map("terminal-residue.map");
	 // This bit can be deleted later:
	 // 
	 if (terminal_residue_do_rigid_body_refine) 
	    std::cout << "fitting terminal residue with rigid body refinement using " 
		      << add_terminal_residue_n_phi_psi_trials << " trials " 
		      << std::endl;
	 else
	    std::cout << "fitting terminal residue with " 
		      << add_terminal_residue_n_phi_psi_trials << " random trials" 
		      << std::endl;

  	 coot::minimol::molecule mmol = 
	    addres.best_fit_phi_psi(add_terminal_residue_n_phi_psi_trials, 
				    terminal_residue_do_rigid_body_refine,
				    add_terminal_residue_add_other_residue_flag);

	 std::vector<coot::minimol::atom *> mmatoms = mmol.select_atoms_serial();
	 // std::cout << "---- ----- mmol has " << mmatoms.size() 
	 // << " atoms" << std::endl;
	 mmol.check();

	 if (mmol.is_empty()) {
	    
	    // this should not happen:
	    std::cout <<  "WARNING: empty molecule: "
		      << "failed to find a fit for terminal residue"
		      << std::endl;

	 } else { 

	    // check that we are adding some atoms:
	    // 
	    std::vector<coot::minimol::atom *> mmatoms = mmol.select_atoms_serial();
	    if (mmatoms.size() == 0) { 
	       std::cout << "WARNING: failed to find a fit for terminal residue"
			 << std::endl;
	       if (use_graphics_interface_flag) { 
		  GtkWidget *w = create_add_terminal_residue_finds_none_dialog();
		  gtk_widget_show(w);
	       }

	    } else {

	       atom_selection_container_t terminal_res_asc;
	       float bf = default_new_atoms_b_factor;

	       // if this is begin added to a shelx molecule, then we
	       // need to set the occs to 11.0
	       //
	       if (graphics_info_t::molecules[imol].is_from_shelx_ins()) {
		  bf = 11.0;
	       } 
	       terminal_res_asc.mol = (MyCMMDBManager *) mmol.pcmmdbmanager();

	       int SelHnd = terminal_res_asc.mol->NewSelection();
	       terminal_res_asc.mol->SelectAtoms(SelHnd, 0, "*",
						 ANY_RES, // starting resno, an int
						 "*", // any insertion code
						 ANY_RES, // ending resno
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
// 	       std::cout << "-------------- terminal_res_asc --------" << std::endl;
// 	       debug_atom_selection_container(terminal_res_asc);

	       // terminal_res_asc.mol->WritePDBASCII("terminal_res_asc.pdb");
	       
	       // atom_selection_container_t tmp_asc = add_cb_to_terminal_res(terminal_res_asc);

 	       atom_selection_container_t tmp_asc =
 		  add_side_chain_to_terminal_res(terminal_res_asc, res_type);


// 	       std::cout << "-------------- tmp_asc --------" << std::endl;
// 	       debug_atom_selection_container(tmp_asc);

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

	       coot::residue_spec_t rs(unconst_res_p);
	       graphics_info_t::molecules[imol].remove_ter_atoms(rs);

	       if (! immediate_addition_flag) { 
		  make_moving_atoms_graphics_object(tmp_asc);
		  moving_atoms_asc_type = coot::NEW_COORDS_INSERT;
		  graphics_draw();
		  coot::refinement_results_t dummy;
		  if (use_graphics_interface_flag) { 
		     do_accept_reject_dialog("Terminal Residue", dummy);
		  } 
	       } else {

		  if (molecules[imol].is_from_shelx_ins()) { 
		     for (unsigned int i=0; i<tmp_asc.n_selected_atoms; i++) {
			tmp_asc.atom_selection[i]->occupancy = 11.0;
		     }
		  } 
		  molecules[imol_moving_atoms].insert_coords(tmp_asc);
		  graphics_draw();
	       }
	    }
	 }
      }
   }
}


void
graphics_info_t::execute_simple_nucleotide_addition(int imol, const std::string &term_type, 
						    CResidue *res_p, const std::string &chain_id) {

   // debug maybe?
   std::cout << "INFO:: Adding a nucleotide, with term_type: \"" << term_type << "\""
	     << std::endl;

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
   
      if (coot::util::nucleotide_is_DNA(res_p)) { 
	 RNA_or_DNA_str = "DNA";
	 form_str = "B";
      }

      coot::ideal_rna ir(RNA_or_DNA_str, form_str, single_stranded_flag,
			 seq, graphics_info_t::standard_residues_asc.mol);
      ir.use_v3_names();
      CMMDBManager *mol = ir.make_molecule();

      int match_resno;
      int interesting_resno;
      if (term_type == "C" || term_type == "MC") {
	 match_resno = 1;
	 interesting_resno = 2;
      } else {
	 interesting_resno = 1;
	 match_resno = 2;
      }

      CResidue *moving_residue_p = NULL;
      CResidue *interesting_residue_p = NULL;
      int imod = 1;
      // now set moving_residue_p and interesting_residue_p:
      CModel *model_p = mol->GetModel(imod);
      CChain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 PCResidue residue_p;
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
	       coot::util::transform_mol(mol, rtop_pair.second);
	       // byte gz = GZM_NONE;
	       // mol->WritePDBASCII("overlapped.pdb", gz);
	       CMMDBManager *residue_mol =
		  coot::util::create_mmdbmanager_from_residue(mol, interesting_residue_p);

	       atom_selection_container_t asc = make_asc(residue_mol);
	       // set the chain id of the chain that contains interesting_residue_p:
	       CModel *model_p = residue_mol->GetModel(imod);
	       CChain *chain_p;
	       // run over chains of the existing mol
	       int nchains = model_p->GetNumberOfChains();
	       for (int ichain=0; ichain<nchains; ichain++) {
		  chain_p = model_p->GetChain(ichain);
		  int nres = chain_p->GetNumberOfResidues();
		  PCResidue residue_p;
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

   // now we are called by chain and molecule pick (as well as the old
   // zone pick).

   // We use the rot_trans_object_type to distinguish.

   const char *chain_id = "*";
   int resno_1 = ANY_RES;
   int resno_2 = ANY_RES;
   std::string insertion_code_selection = "*"; // reset on start and stop residue in range being the same
   bool good_settings = 0; // fail initially
   CAtom *atom1 = molecules[imol_rot_trans_object].atom_sel.atom_selection[rot_trans_atom_index_1];
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

   if (rot_trans_object_type == ROT_TRANS_TYPE_ZONE) {
      CAtom *atom2 = molecules[imol_rot_trans_object].atom_sel.atom_selection[rot_trans_atom_index_2];
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

      GtkWidget *widget = create_rotate_translate_obj_dialog();
      GtkWindow *main_window = GTK_WINDOW(lookup_widget(graphics_info_t::glarea,
							"window1"));
      gtk_window_set_transient_for(GTK_WINDOW(widget), main_window);

      do_rot_trans_adjustments(widget);

      // set its position if it was shown before
      if (rotate_translate_x_position > -100) {
	 gtk_widget_set_uposition(widget,
				  rotate_translate_x_position,
				  rotate_translate_y_position);
      }
      gtk_widget_show(widget);

      atom_selection_container_t rt_asc;
      // No! It cannot point to the same CAtoms.
      // rt_asc.mol = molecules[imol_rot_trans_object].atom_sel.mol; 
      // MyCMMDBManager *mol = new MyCMMDBManager;
      // mol->Copy(molecules[imol_rot_trans_object].atom_sel.mol, MMDBFCM_All);
      // how about we instead use:
      // CMMDBManager *mol = create_mmdbmanager_from_res_selection();
      //
      PCResidue *sel_residues = NULL;
      int n_sel_residues;
      int selHnd = molecules[imol_rot_trans_object].atom_sel.mol->NewSelection();
      molecules[imol_rot_trans_object].atom_sel.mol->Select(selHnd, STYPE_RESIDUE, 0,
							    chain_id,
							    resno_1, insertion_code_selection.c_str(),
							    resno_2, insertion_code_selection.c_str(),
							    "*",  // residue name
							    "*",  // Residue must contain this atom name?
							    "*",  // Residue must contain this Element?
							    "*",  // altLocs
							    SKEY_NEW // selection key
							    );
      molecules[imol_rot_trans_object].atom_sel.mol->GetSelIndex(selHnd, sel_residues, n_sel_residues);

      
      short int alt_conf_split_flag = 0;
      std::string altloc_string(altLoc);
      if (altloc_string != "")
	 alt_conf_split_flag = 1;

      // create a complete new clean copy of chains/residues/atoms
      std::pair<CMMDBManager *, int> mp(0, 0);


      if (rot_trans_object_type == ROT_TRANS_TYPE_ZONE) 
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
      make_moving_atoms_graphics_object(rt_asc); // shallow copy rt_asc to moving_atoms_asc

      // set the rotation centre atom index:
      //   rot_trans_atom_index_rotation_origin_atom = 
      //       find_atom_index_in_moving_atoms(chain_id, 
      // 				      atom2->GetSeqNum(),
      // 				      atom2->name);  // uses moving_atoms_asc


      rot_trans_rotation_origin_atom = find_atom_in_moving_atoms(origin_atom_spec);

      if (0) { 
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
	       
	       CAtom *atom_1 = 0; 
	       CAtom *atom_2 = 0; 
	       CAtom *atom_3 = 0; 
	       CAtom *atom_4 = 0;
	       int im = torsion_general_atom_index_1_mol_no;

	       if (torsion_general_atom_index_1 < molecules[im].atom_sel.n_selected_atoms) { 
		  if (torsion_general_atom_index_2 < molecules[im].atom_sel.n_selected_atoms) { 
		     if (torsion_general_atom_index_3 < molecules[im].atom_sel.n_selected_atoms) { 
			if (torsion_general_atom_index_4 < molecules[im].atom_sel.n_selected_atoms) {

			   atom_1 = molecules[im].atom_sel.atom_selection[torsion_general_atom_index_1];
			   atom_2 = molecules[im].atom_sel.atom_selection[torsion_general_atom_index_2];
			   atom_3 = molecules[im].atom_sel.atom_selection[torsion_general_atom_index_3];
			   atom_4 = molecules[im].atom_sel.atom_selection[torsion_general_atom_index_4];

			   CResidue *r1 = atom_1->GetResidue();
			   CResidue *r2 = atom_2->GetResidue();
			   CResidue *r3 = atom_3->GetResidue();
			   CResidue *r4 = atom_4->GetResidue();

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
				    make_moving_atoms_graphics_object(residue_asc);

				    std::vector<coot::atom_spec_t> as;
				    as.push_back(atom_1);
				    as.push_back(atom_2);
				    as.push_back(atom_3);
				    as.push_back(atom_4);
				    torsion_general_atom_specs = as;
				    graphics_draw();
				    torsion_general_reverse_flag = 0;
				    CResidue *res_local = get_first_res_of_moving_atoms();
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

CResidue *
graphics_info_t::get_first_res_of_moving_atoms() {

   CResidue *r = 0;
   CModel *model_p = moving_atoms_asc->mol->GetModel(1);
   if (model_p) {
      CChain *chain_p = model_p->GetChain(0);
      if (chain_p) {
	 CResidue *residue_p = chain_p->GetResidue(0);
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
      GtkWidget *hscale = lookup_widget(dialog, hscale_lab[i].c_str());
      GtkAdjustment *adj = GTK_ADJUSTMENT(gtk_adjustment_new(0.0, -180.0, 360.0, 0.1, 1.0, 0));
      gtk_range_set_adjustment(GTK_RANGE(hscale), GTK_ADJUSTMENT(adj));
      gtk_signal_connect(GTK_OBJECT(adj), 
			 "value_changed",
			 GTK_SIGNAL_FUNC(graphics_info_t::rot_trans_adjustment_changed), 
			 GINT_TO_POINTER(i));
   }
}


coot::ScreenVectors::ScreenVectors() {
   
   coot::Cartesian centre = unproject_xyz(0, 0, 0.5);
   coot::Cartesian front  = unproject_xyz(0, 0, 0.0);
   coot::Cartesian right  = unproject_xyz(1, 0, 0.5);
   coot::Cartesian top    = unproject_xyz(0, 1, 0.5);

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

   graphics_info_t g;  // because rotate_round_vector is not static - it should be.  
                       // FIXME at some stage.
   double v = adj->value;
   
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
      CAtom *rot_centre = rot_trans_rotation_origin_atom;
      clipper::Coord_orth rotation_centre(0,0,0); // updated.

      // But! maybe we have a different rotation centre
      if (rot_trans_zone_rotates_about_zone_centre) {
	 if (moving_atoms_asc->n_selected_atoms  > 0) {
	    graphics_info_t g;
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
	 clipper::Coord_orth new_pos = 
	    g.rotate_round_vector(screen_vector, co, rotation_centre, x_diff * 0.018);
	 moving_atoms_asc->atom_selection[i]->x = new_pos.x();
	 moving_atoms_asc->atom_selection[i]->y = new_pos.y();
	 moving_atoms_asc->atom_selection[i]->z = new_pos.z();
      }
   } else { 

      for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
	 moving_atoms_asc->atom_selection[i]->x += x_add;
	 moving_atoms_asc->atom_selection[i]->y += y_add;
	 moving_atoms_asc->atom_selection[i]->z += z_add;
      }
   }
   int do_disulphide_flag = 0;

   if (molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS ||
       molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS_PLUS_LIGANDS ||
       molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR ||
       molecules[imol_moving_atoms].Bonds_box_type() == coot::COLOUR_BY_RAINBOW_BONDS) {
      
      Bond_lines_container bonds;
      bonds.do_Ca_plus_ligands_bonds(*moving_atoms_asc, 1.0, 4.7);
      regularize_object_bonds_box.clear_up();
      regularize_object_bonds_box = bonds.make_graphical_bonds();
   } else {
      Bond_lines_container bonds(*moving_atoms_asc, do_disulphide_flag);
      regularize_object_bonds_box.clear_up();
      regularize_object_bonds_box = bonds.make_graphical_bonds();
   }
   graphics_draw();
}


// --- nudge active residue
// static
void
graphics_info_t::nudge_active_residue(guint direction) {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = graphics_info_t::active_atom_spec();
   if (active_atom.first) {
      clipper::Coord_orth shift(0,0,0);
      clipper::Mat33<double> mat(1,0,0,0,1,0,0,0,1);
      double shift_scale_factor = 0.01 * zoom; // needs to be 0.04 for funny mode?
      coot::ScreenVectors screen_vectors;

      if (direction == GDK_Left) { 
	 // std::cout << "Left nudge residue" << std::endl;
	 shift = clipper::Coord_orth(-shift_scale_factor * screen_vectors.screen_x.x(),
				     -shift_scale_factor * screen_vectors.screen_x.y(),
				     -shift_scale_factor * screen_vectors.screen_x.z());
      } 
      if (direction == GDK_Right) { 
	 // std::cout << "Right nudge residue" << std::endl;
	 shift = clipper::Coord_orth(shift_scale_factor * screen_vectors.screen_x.x(),
				     shift_scale_factor * screen_vectors.screen_x.y(),
				     shift_scale_factor * screen_vectors.screen_x.z());
      } 
      if (direction == GDK_Up) { 
	 // std::cout << "Up nudge residue" << std::endl;
	 shift = clipper::Coord_orth(-shift_scale_factor * screen_vectors.screen_y.x(),
				     -shift_scale_factor * screen_vectors.screen_y.y(),
				     -shift_scale_factor * screen_vectors.screen_y.z());
      } 
      if (direction == GDK_Down) { 
	 // std::cout << "Down nudge residue" << std::endl;
	 shift = clipper::Coord_orth(shift_scale_factor * screen_vectors.screen_y.x(),
				     shift_scale_factor * screen_vectors.screen_y.y(),
				     shift_scale_factor * screen_vectors.screen_y.z());
      } 

      // all constructed.  Apply it
      clipper::RTop_orth rtop(mat, shift);
      int imol = active_atom.second.first;
      graphics_info_t::molecules[imol].transform_zone_by(active_atom.second.second.chain,
							 active_atom.second.second.resno,
							 active_atom.second.second.resno,
							 active_atom.second.second.insertion_code,
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
 


void
graphics_info_t::execute_db_main() { 

   int imol = db_main_imol;
   CAtom *at1 = molecules[imol].atom_sel.atom_selection[db_main_atom_index_1];
   CAtom *at2 = molecules[imol].atom_sel.atom_selection[db_main_atom_index_2];
   std::string chain_id = at1->GetChainID();
   int iresno_start = at1->GetSeqNum();
   int iresno_end   = at2->GetSeqNum();

   std::string direction_string("forwards"); // forwards
   execute_db_main(imol, chain_id, iresno_start, iresno_end, direction_string);
}

// 
void
graphics_info_t::execute_db_main(int imol,
				 std::string chain_id,
				 int iresno_start,
				 int iresno_end,
				 std::string direction_string) { 

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

      // mt is a minimol of the Baton Atoms:
      coot::minimol::molecule mt(molecules[imol].atom_sel.mol);
      coot::minimol::molecule target_ca_coords;

      if (direction_string != "backwards") { 
	 for (unsigned int i=0; i<mt.fragments.size(); i++)
	    if (mt.fragments[i].fragment_id == chain_id)
	       target_ca_coords.fragments.push_back(mt.fragments[i]);
      } else { // backwards code.
	 for (unsigned int i=0; i<mt.fragments.size(); i++) {
	    if (mt[i].fragment_id == chain_id) {
	       // put in the residues of mt.fragments[i] backwards:
	       
	       // The seqnum of the residues is ignored, the only
	       // important thing is the ires.
	       
	       int ifrag = target_ca_coords.fragment_for_chain(chain_id);
	       if (mt[i].max_residue_number() > 1) { 
		  for (int ires=mt[i].max_residue_number(); ires>=mt[i].min_res_no(); ires--) {
		     target_ca_coords[ifrag].residues.push_back(mt[ifrag][ires]);
		  }
		  break;
	       }
	    }
	 }
      }

      // now target_ca_coords has only one chain, the chain of the zone.
      // Note that match_target_fragment selects CAs from target_ca_coords
      // so we don't need to filter them out here.


      // write out target_ca_coords:
      // target_ca_coords.write_file("target_ca_coords.pdb", 20);

      main_chain.match_target_fragment(target_ca_coords,
				       iresno_start,
				       iresno_end,
				       ilength);

      main_chain.merge_fragments();
      coot::minimol::molecule mol;
      mol.fragments.push_back(main_chain.mainchain_fragment());
      mol.write_file("db-mainchain.pdb", 20.0);

      // std::cout << "DEBUG:: mol.is_empty() returns " << mol.is_empty() << std::endl;
      std::vector<coot::minimol::atom *> serial_atoms = mol.select_atoms_serial();
      // std::cout << "DEBUG:: serial_atoms.size() returns " << serial_atoms.size() << std::endl;
      
      if (serial_atoms.size() > 0) {
	 std::pair<std::vector<float>, std::string> cell_spgr = 
	    molecules[imol].get_cell_and_symm();
	 float bf = default_new_atoms_b_factor;
	 atom_selection_container_t asc = make_asc(mol.pcmmdbmanager());
	 set_mmdb_cell_and_symm(asc, cell_spgr); // tinker with asc. 
	                                         // Consider asc as an object.
	 int imol_new = create_molecule();
	 molecules[imol_new].install_model(imol_new, asc, "mainchain", 1);
	 graphics_draw();
      } else {
	 std::string s("Sorry, failed to convert that residue range.\nToo short, perhaps?");
	 GtkWidget *w = wrapped_nothing_bad_dialog(s);
	 gtk_widget_show(w);
      }
      main_chain.clear_results();
   }
}

// --------------------------------------------------------------------------------
//                 Rotamer stuff
// --------------------------------------------------------------------------------

void
graphics_info_t::do_rotamers(int atom_index, int imol) {


   if (use_graphics_interface_flag) { 
      // display the buttons for the rotamer options and display
      // the most likely in the graphics as a
      // moving_atoms_asc.

      rotamer_residue_atom_index = atom_index;  // save for button
      // callbacks, so that we
      // can get the residue.
      rotamer_residue_imol = imol;
      std::string altconf = molecules[imol].atom_sel.atom_selection[atom_index]->altLoc;

//       std::cout << "DEBUG:: in do_rotamers() atom_index is " << atom_index
// 		<< " and alconf is :" <<  altconf << ":" << std::endl;
   
      GtkWidget *window = create_rotamer_selection_dialog();
      set_transient_and_position(COOT_ROTAMER_SELECTION_DIALOG, window);
      rotamer_dialog = window;

      // Test if this was an alt confed atom.
      // If it was, then we should set up the hscale.
      // It it was not, then we should destroy the hscale
      if (altconf.length() > 0) {
	 GtkWidget *hscale = lookup_widget(window, "new_alt_conf_occ_hscale");
	 float v = add_alt_conf_new_atoms_occupancy;
	 // The max value is 3rd arg - 6th arg (here 2 and 1 is the same as 1 and 0)
	 GtkAdjustment *adj = GTK_ADJUSTMENT(gtk_adjustment_new(v, 0.0, 2.0, 0.01, 0.1, 1.0));
	 gtk_range_set_adjustment(GTK_RANGE(hscale), GTK_ADJUSTMENT(adj));
	 gtk_signal_connect(GTK_OBJECT(adj), 
			    "value_changed",
			    GTK_SIGNAL_FUNC(graphics_info_t::new_alt_conf_occ_adjustment_changed), 
			    NULL);
      
      } else {
	 GtkWidget *frame = lookup_widget(window, "new_alt_conf_occ_frame");
	 gtk_widget_destroy(frame);
      }

   
      /* Events for widget must be set before X Window is created */
      gtk_widget_set_events(GTK_WIDGET(window),
			    GDK_KEY_PRESS_MASK);
      /* Capture keypress events */
      //    rotamer_key_press_event is not defined (yet)
      //    gtk_signal_connect(GTK_OBJECT(window), "key_press_event",
      // 		      GTK_SIGNAL_FUNC(rotamer_key_press_event), NULL);
      /* set focus to glarea widget - we need this to get key presses. */
      GTK_WIDGET_SET_FLAGS(window, GTK_CAN_FOCUS);
      gtk_widget_grab_focus(GTK_WIDGET(glarea)); // but set focus to the graphics.
   
      fill_rotamer_selection_buttons(window, atom_index, imol);

      // act as if the button for the first rotamer was pressed
      short int stat = generate_moving_atoms_from_rotamer(0);

      if (stat)
	 gtk_widget_show(window);
   }
}


// static
void graphics_info_t::new_alt_conf_occ_adjustment_changed(GtkAdjustment *adj,
							  gpointer user_data) {

   graphics_info_t g;
   g.add_alt_conf_new_atoms_occupancy = adj->value;

   // Change the occupancies of the intermediate atoms:
   //
   if (moving_atoms_asc) {
      for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
	 // this if test is a kludge!
	 // Don't change the alt conf for fully occupied atoms.
	 if (moving_atoms_asc->atom_selection[i]->occupancy < 0.99) 
	    moving_atoms_asc->atom_selection[i]->occupancy = adj->value;
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
      int imod = 1;
      CModel *model_p = moving_atoms_asc->mol->GetModel(imod);
      CChain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 PCResidue residue_p;
	 CAtom *at;
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
   }
   Bond_lines_container bonds(*moving_atoms_asc, geom_p, 0, 1, 0);
   regularize_object_bonds_box.clear_up();
   regularize_object_bonds_box = bonds.make_graphical_bonds();
   graphics_draw();
}


// static
void
graphics_info_t::mark_atom_as_fixed(int imol, const coot::atom_spec_t &atom_spec, bool state) {
   if (!moving_atoms_asc) {
      std::cout << "WARNING:: No intermediate atoms - fail" << std::endl;
   } else {
      if ((imol >=0) && (imol < n_molecules())) {
	 if (graphics_info_t::molecules[imol].has_model()) {
	    graphics_info_t::molecules[imol].mark_atom_as_fixed(atom_spec, state);
	 }
      }
   }      
}

void
graphics_info_t::fill_rotamer_selection_buttons(GtkWidget *window, int atom_index, int imol) const {

   // for each rotamer do this:

   GSList *gr_group = NULL;
   GtkWidget *rotamer_selection_radio_button;
   GtkWidget *rotamer_selection_dialog = window;
   GtkWidget *rotamer_selection_button_vbox =
      lookup_widget(window, "rotamer_selection_button_vbox");
   graphics_info_t g;
   std::string alt_conf = g.molecules[imol].atom_sel.atom_selection[atom_index]->altLoc;
   CResidue *residue = g.molecules[imol].atom_sel.atom_selection[atom_index]->residue;
      
#ifdef USE_DUNBRACK_ROTAMERS			
      coot::dunbrack d(residue, g.molecules[imol].atom_sel.mol, g.rotamer_lowest_probability, 0);
#else
      std::cout << " in fill_rotamer_selection_buttons altconf is :" << alt_conf << ":"
		<< std::endl;
      coot::richardson_rotamer d(residue, alt_conf,
				 g.molecules[imol].atom_sel.mol, g.rotamer_lowest_probability, 0);
#endif // USE_DUNBRACK_ROTAMERS

   std::vector<float> probabilities = d.probabilities();
   // std::cout << "There are " << probabilities.size() << " probabilities"
   // << std::endl;

   // Attach the number of residues to the dialog so that we can get
   // that data item when we make a synthetic key press due to
   // keyboard (arrow?) key press:
   gtk_object_set_user_data(GTK_OBJECT(window), GINT_TO_POINTER(probabilities.size()));

   GtkWidget *frame;
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
   
      rotamer_selection_radio_button =
	 gtk_radio_button_new_with_label (gr_group, button_label.c_str());
      gr_group = gtk_radio_button_group (GTK_RADIO_BUTTON (rotamer_selection_radio_button));
      gtk_widget_ref (rotamer_selection_radio_button);
      gtk_object_set_data_full (GTK_OBJECT (rotamer_selection_dialog),
				button_name.c_str(), rotamer_selection_radio_button,
				(GtkDestroyNotify) gtk_widget_unref);
      
      int *iuser_data = new int;
      *iuser_data = i;
      gtk_signal_connect (GTK_OBJECT (rotamer_selection_radio_button), "toggled",
			  GTK_SIGNAL_FUNC (on_rotamer_selection_button_toggled),
			  iuser_data);
       
       gtk_widget_show (rotamer_selection_radio_button);
       frame = gtk_frame_new(NULL);
       gtk_container_add(GTK_CONTAINER(frame), rotamer_selection_radio_button);
       gtk_box_pack_start (GTK_BOX (rotamer_selection_button_vbox),
			   frame, FALSE, FALSE, 0);
       gtk_container_set_border_width (GTK_CONTAINER (frame), 2);
       gtk_widget_show(frame);
   }
} 


void
graphics_info_t::on_rotamer_selection_button_toggled (GtkButton       *button,
						      gpointer         user_data) {

   int *i_tmp = (int *) user_data;
   int i = *i_tmp;
   
   graphics_info_t g;
   g.generate_moving_atoms_from_rotamer(i);
   
} 

// Return 1 for valid (i.e. non-GLY, non-ALA) residue, 0 otherwise
// (including residue type not found).
// 
short int 
graphics_info_t::generate_moving_atoms_from_rotamer(int irot) {

   int imol = rotamer_residue_imol;
   int atom_index = rotamer_residue_atom_index; 

   CAtom    *at_rot   = molecules[imol].atom_sel.atom_selection[atom_index];
   CResidue *residue  = molecules[imol].atom_sel.atom_selection[atom_index]->residue;
   int atom_index_udd = molecules[imol].atom_sel.UDDAtomIndexHandle;
   std::string altconf = at_rot->altLoc;

   if (std::string(residue->name) == "GLY" ||
       std::string(residue->name) == "ALA") {
      std::cout << "INFO:: This residue ("<< residue->name
		<< ") doesn't have rotamers\n";
      return 0;
   }

   // We need to filter out atoms that are not (either the same
   // altconf as atom_index or "")
   // 
   CResidue *tres = coot::deep_copy_this_residue(residue, 
						 std::string(at_rot->altLoc),
						 0, atom_index_udd);
   if (!tres) {
      return 0;
   } else { 
      PPCAtom residue_atoms;
      int nResidueAtoms;
      std::string mol_atom_altloc;
      std::string atom_altloc = molecules[imol].atom_sel.atom_selection[atom_index]->altLoc;
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
	 Geom_p()->get_monomer_restraints(monomer_type);

#ifdef USE_DUNBRACK_ROTAMERS			
      coot::dunbrack d(tres, molecules[imol].atom_sel.mol,
		       rotamer_lowest_probability, 0);
#else
      coot::richardson_rotamer d(tres, altconf, molecules[imol].atom_sel.mol,
				 rotamer_lowest_probability, 0);
#endif // USE_DUNBRACK_ROTAMERS

      if (p.first) { 
	 // std::cout << "generate_moving_atoms_from_rotamer " << irot << std::endl;
	 // The magic happens here:
	 CResidue *moving_res = d.GetResidue(p.second, irot);

	 //
	 if (moving_res == NULL) {
	    std::cout << "Failure to find rotamer for residue type: "
		      << residue->name << std::endl;
	    return 0;
	 } else { 

	    MyCMMDBManager *mol = new MyCMMDBManager;
	    CModel *model_p = new CModel;
	    CChain *chain_p = new CChain;
	    CResidue *res_p = new CResidue;
	    res_p->SetResID(residue->GetResName(),
			    residue->GetSeqNum(),
			    residue->GetInsCode());
   
	    PPCAtom residue_atoms_2;
	    int nResidueAtoms_2;
	    ((CResidue *)moving_res)->GetAtomTable(residue_atoms_2, nResidueAtoms_2);
	    CAtom *atom_p;
	    int i_add;
	    for(int iat=0; iat<nResidueAtoms_2; iat++) {
	       atom_p = new CAtom;
	       atom_p->Copy(residue_atoms_2[iat]);
	       i_add = res_p->AddAtom(atom_p);
	    }
	    chain_p->AddResidue(res_p);
	    chain_p->SetChainID(residue->GetChainID());
	    model_p->AddChain(chain_p);
	    mol->AddModel(model_p);
	    mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
	    mol->FinishStructEdit();

	    imol_moving_atoms = imol;
	    *moving_atoms_asc = make_asc(mol);
	    //    std::cout << "there are " << moving_atoms_asc->n_selected_atoms
	    // 	     << " selected atoms in the moving_atoms_asc" << std::endl;

	    moving_atoms_asc_type = coot::NEW_COORDS_REPLACE_CHANGE_ALTCONF;
	    make_moving_atoms_graphics_object(*moving_atoms_asc);
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
graphics_info_t::get_rotamer_probability(CResidue *res,
					 const std::string &altconf,
					 CMMDBManager *mol,
					 float lowest_probability,
					 short int add_extra_PHE_and_TYR_rotamers_flag) {

   coot::rotamer_probability_info_t r(coot::rotamer_probability_info_t::MISSING_ATOMS,0,"");
#ifdef USE_DUNBRACK_ROTAMERS			
   coot::dunbrack d(res, mol, rotamer_lowest_probability, 1);
#else
   if (!rot_prob_tables.is_well_formatted()) {
      rot_prob_tables.fill_tables();
   }
   if (rot_prob_tables.is_well_formatted()) {
      try {
	 std::vector<coot::rotamer_probability_info_t> v = rot_prob_tables.probability_this_rotamer(res);
	 if (v.size() > 0) {
	    r = v[0];
	 } 
      }
      catch (std::runtime_error e) {
	 std::cout << "get_rotamer_probability: caught: " << e.what() << std::endl;
      } 
   } else {
      coot::richardson_rotamer d(res, altconf, mol, rotamer_lowest_probability, 1);
      r  = d.probability_of_this_rotamer();
   } 
#endif // USE_DUNBRACK_ROTAMERS

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

   int imol = create_pointer_atom_molecule_maybe();
   if (imol >= 0 && imol < n_molecules()) {
      if (molecules[imol].open_molecule_p()) { 
	 molecules[imol].add_typed_pointer_atom(RotationCentre(), type); // update bonds

	 update_environment_distances_by_rotation_centre_maybe(imol);
	 
	 graphics_draw();
      } else {
	 std::cout << "ERROR:: invalid (closed) molecule in place_typed_atom_at_pointer: "
		   << imol << std::endl;
      }
   } else {
      std::cout << "WARNING:: invalid molecule number in place_typed_atom_at_pointer"
		<< imol << std::endl;
   }
}

// Tinker with the atom positions of residue
// Return 1 on success.
// We need to pass the asc for the mol because we need it for seekcontacts()
// Of course the asc that is passed is the moving atoms asc.
//
// nth_chi is 1-based (i.e. rotating about CA-CB, nth_chi is 1).
// 
short int 
graphics_info_t::update_residue_by_chi_change(CResidue *residue,
					      atom_selection_container_t &asc,
					      int nth_chi, double diff) {
   short int istat = 0;
   double angle = diff/60.0;
   bool reverse = edit_chi_angles_reverse_fragment; 

   std::string monomer_type = residue->GetResName();
   // this can throw an exception
   std::pair<short int, coot::dictionary_residue_restraints_t> p =
      Geom_p()->get_monomer_restraints(monomer_type);
   
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
	 catch (std::runtime_error rte) {
	    // std::cout << rte.what() << std::endl;
	    int base_atom_index = 0;

	    // tmp hack for testing.
	    // coot::contact_info contact = coot::getcontacts(*moving_atoms_asc);

	    // c.f. get_contact_indices_from_restraints().  urgh. Same
	    // functionality: "written twice".
	    // 
	    coot::contact_info contact = coot::getcontacts(*moving_atoms_asc, monomer_type, Geom_p());
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
	    catch (std::runtime_error rte) {
	       std::cout << "Update chi - contact fall-back fails - " << rte.what() << std::endl;
	    }
	 }
      }
      catch (std::runtime_error rte) {
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

// this can throw an exception.
std::pair<std::string, std::string>
graphics_info_t::get_chi_atom_names(CResidue *residue,
				    const coot::dictionary_residue_restraints_t &rest,
				    int nth_chi) const {

   std::pair<std::string, std::string> r(" CA ", " CB ");
   int torsion_index = nth_chi -1;
   std::vector <coot::dict_torsion_restraint_t> torsion_restraints =
      rest.get_non_const_torsions(find_hydrogen_torsions_flag);
   
   if ((torsion_index >=0) && (torsion_index < torsion_restraints.size())) {
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






// Called by mouse motion callback (in_edit_chi_mode_flag)
// 
void
graphics_info_t::rotate_chi(double x, double y) {

   // real values start at 1:
   int chi = edit_chi_current_chi;

   mouse_current_x = x;
   mouse_current_y = y;
   double diff;

   diff  = mouse_current_x - GetMouseBeginX();
   diff += mouse_current_y - GetMouseBeginY();

   diff *= 15;

   // std::cout << "graphics_info_t::rotate_chi " << chi << " by "
   // << diff << std::endl;

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
	 CModel *model_p = moving_atoms_asc->mol->GetModel(1);
	 if (model_p) {
	    CChain *chain_p = model_p->GetChain(0);
	    if (chain_p) {
	       CResidue *residue_p = chain_p->GetResidue(0);
	       if (residue_p) {
		  istat = update_residue_by_chi_change(residue_p, *moving_atoms_asc, chi, diff);
	       }
	    }
	 }
      }
   }
   
   if (istat == 0) {
      // std::cout << "regenerating object" << std::endl;
      regularize_object_bonds_box.clear_up();
      make_moving_atoms_graphics_object(*moving_atoms_asc); // make new bonds
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
      CResidue *residue_p = get_first_res_of_moving_atoms();
      if (residue_p) {

	 std::string altconf = chi_angle_alt_conf;
	 try {
	    // use class variable (previous saved)
	    int base_atom_index = 0;
	    coot::atom_tree_t tree(torsion_general_contact_indices, base_atom_index, residue_p, altconf);
	    tree.rotate_about(specs_local[1].atom_name, specs_local[2].atom_name,
			      diff, torsion_general_reverse_flag);
	    regularize_object_bonds_box.clear_up();
	    make_moving_atoms_graphics_object(*moving_atoms_asc);
	    graphics_draw();
	 }
	 catch (std::runtime_error rte) {
	    std::cout << "INFO:: tree by contacts failed " << rte.what() << std::endl;
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

   std::vector<CResidue *> residues;
   for (unsigned int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
      CResidue *r = moving_atoms_asc->atom_selection[i]->residue;
      if (std::find(residues.begin(), residues.end(), r) == residues.end())
	 residues.push_back(r);
   }

   std::vector<std::pair<CAtom *, CAtom *> > link_bond_atom_pairs = 
      coot::torsionable_link_bonds(residues, moving_atoms_asc->mol, Geom_p());
   coot::contact_info contacts(*moving_atoms_asc, geom_p, link_bond_atom_pairs);
   std::vector<std::vector<int> > contact_indices =
      contacts.get_contact_indices_with_reverse_contacts();
   try { 
      coot::atom_tree_t tree(contact_indices, 0,
			     moving_atoms_asc->mol,
			     moving_atoms_asc->SelectionHandle);

      int index_1 = multi_residue_torsion_rotating_atom_index_pair.first;
      int index_2 = multi_residue_torsion_rotating_atom_index_pair.second;
      tree.rotate_about(index_1, index_2, diff, multi_residue_torsion_reverse_fragment_mode);
      make_moving_atoms_graphics_object(*moving_atoms_asc);
      graphics_draw();
   }
   catch (std::runtime_error rte) {
      std::cout << "WARNING:: " << rte.what() << std::endl;
   }
}



// 	 //  debug:
// 	 std::cout << "DEBUG:: residue_mol: ----------------- " << std::endl;
// 	 CModel *model_p = residues_mol->GetModel(1);
// 	 CChain *chain_p;
// 	 int nchains = model_p->GetNumberOfChains();
// 	 std::cout << "DEBUG:: residue_mol: nchains " << nchains << std::endl;
// 	 for (int ichain=0; ichain<nchains; ichain++) { 
// 	    chain_p = model_p->GetChain(ichain);
// 	    int nres = chain_p->GetNumberOfResidues();
// 	    for (int ires=0; ires<nres; ires++) { 
// 	       PCResidue residue_p = chain_p->GetResidue(ires);
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
   
   CResidue *r = molecules[imol].get_residue(chain_id, resno, ins_code);
   if (!r) {
      std::cout << "WARNING:: Residue " << " chain-id :" << chain_id << ":  resno: " << resno
		<< " inscode :" << ins_code << ": not found" << std::endl;
   } else {
      PPCAtom residue_atoms;
      int n_residue_atoms;
      int at_index = -1;
      r->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++) {
	 std::string atom_name(residue_atoms[i]->name);
	 std::string atom_alt_conf(residue_atoms[i]->altLoc);
	 std::cout << "   " << i << " " << atom_name << " :" << atom_alt_conf << ":" << std::endl;
	 if (atom_alt_conf == altconf) {
	    CAtom *at = residue_atoms[i];
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
				      const coot::residue_spec_t &res1,
				      const coot::residue_spec_t &res2) {

   molecules[imol].delete_zone(res1, res2);
   if (delete_item_widget) {
      GtkWidget *checkbutton = lookup_widget(graphics_info_t::delete_item_widget,
					     "delete_item_keep_active_checkbutton");
      if (GTK_TOGGLE_BUTTON(checkbutton)->active) {
	 // don't destroy it.
      } else {
	 gint upositionx, upositiony;
	 gdk_window_get_root_origin (delete_item_widget->window, &upositionx, &upositiony);
	 delete_item_widget_x_position = upositionx;
	 delete_item_widget_y_position = upositiony;
	 gtk_widget_destroy(delete_item_widget);
	 delete_item_widget = 0;
	 normal_cursor();
      }
   }

   if ((imol >=0) && (imol < n_molecules())) {
      graphics_info_t::molecules[imol].delete_zone(res1, res2);
      if (graphics_info_t::go_to_atom_window) {
	 update_go_to_atom_window_on_changed_mol(imol);
      }
   }
   graphics_draw();
}


// static
void
graphics_info_t::move_molecule_here_item_select(GtkWidget *item,
						GtkPositionType pos) {

   graphics_info_t::move_molecule_here_molecule_number = pos;

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
   
#if defined USE_GUILE_GTK && !defined WINDOWS_MINGW
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

#ifdef HAVE_GSL

   if (moving_atoms_asc) { 
      if (moving_atoms_asc_type == coot::NEW_COORDS_REPLACE ||
	  moving_atoms_asc_type == coot::NEW_COORDS_REPLACE_CHANGE_ALTCONF) { // needed?
	 if (moving_atoms_asc->mol) {

	    std::string s = "Unset";

	    // Chirals:
	    std::pair<std::vector<std::string> , std::vector <coot::atom_spec_t> >
	       bv = coot::inverted_chiral_volumes(moving_atoms_asc->mol,
						  geom_p,
						  cif_dictionary_read_number);
	    if (bv.second.size() > 0) {
	       if (bv.second.size() == 1) {
		  int i = 0;
		  s = "There is one residue with an\n";
		  s += "incorrect chiral volume\n";
		  s += bv.second[i].chain;
		  s += " ";
		  s += coot::util::int_to_string(bv.second[i].resno);
		  s += bv.second[i].insertion_code;
		  s += " ";
		  s += bv.second[i].atom_name;
		  s += " ";
		  s += bv.second[i].alt_conf;
		  s += "\n";
	       } else {
		  s = "There are ";
		  s += coot::util::int_to_string(bv.second.size());
		  s += " residues with \n";
		  s += "incorrect chiral volumes\n";
		  for (unsigned int i=0; i<bv.second.size(); i++) { 
		     s += bv.second[i].chain;
		     s += " ";
		     s += coot::util::int_to_string(bv.second[i].resno);
		     s += bv.second[i].insertion_code;
		     s += " ";
		     s += bv.second[i].atom_name;
		     s += " ";
		     s += bv.second[i].alt_conf;
		     s += "\n";
		  }
	       }
	    }

	    // Cis peptides:

	    int n_cis = coot::util::count_cis_peptides(moving_atoms_asc->mol);
	    // std::cout << "DEBUG:: End ref: have " << n_cis << " CIS peptides "
	    // << std::endl;

	    if (n_cis > graphics_info_t::moving_atoms_n_cis_peptides) {
	       if (n_cis == 1) {
		  s += "\nWARNING: A CIS peptide has been introduced\n";
	       } else {
		  if ((n_cis - graphics_info_t::moving_atoms_n_cis_peptides) > 1) {
		     s += "\nWARNING: An extra CIS peptide has been introduced\n";
		  } else { 
		     s += "\nWARNING: Extra CIS peptides have been introduced\n";
		  }
	       }
	    }
	    
	    if (show_chiral_volume_errors_dialog_flag) {
	       if (accept_reject_dialog) { 
		  // info_dialog(s);
		  if (s != "Unset") {
		     update_accept_reject_dialog_with_results(accept_reject_dialog,
							      coot::CHIRAL_CENTRES,
							      s);
		  } else { 
		     update_accept_reject_dialog_with_results(accept_reject_dialog,
							      coot::CHIRAL_CENTRES,
							      coot::refinement_results_t(" "));
		  }
	       }
	       if (s != "Unset") {
		  std::cout << s << std::endl;
	       }
	    } 
	 }
      }
   }
#endif // HAVE_GSL   
}



