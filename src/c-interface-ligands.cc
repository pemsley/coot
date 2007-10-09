/* src/c-interface-ligands.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */
#include <stdlib.h>
#include <iostream>

#if defined _MSC_VER
#include <windows.h>
#endif
 
#include "globjects.h" //includes gtk/gtk.h

#include "callbacks.h"
#include "interface.h" // now that we are moving callback
		       // functionality to the file, we need this
		       // header since some of the callbacks call
		       // fuctions built by glade.
#include <vector>
#include <string>

#include "mmdb_manager.h"
#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb-crystal.h"

#include "graphics-info.h"
#include "c-interface.h"
#include "cc-interface.hh"
#include "wligand.hh"


/*  ----------------------------------------------------------------------- */
/*                  ligand overlay                                          */
/*  ----------------------------------------------------------------------- */
/*! \brief "Template"-based matching.  Overlap the first residue in
  imol_ligand onto the residue specified by the reference parameters.
  Use graph matching, not atom names.  */

#ifdef USE_GUILE
SCM
overlap_ligands(int imol_ligand, int imol_ref, const char *chain_id_ref,
		int resno_ref) {

   SCM scm_status = SCM_BOOL_F;
   coot::graph_match_info_t rtop_info =
      overlap_ligands_internal(imol_ligand, imol_ref, chain_id_ref, resno_ref, 1);

   if (rtop_info.success) { 
      SCM match_info = scm_cons(scm_int2num(rtop_info.n_match), SCM_EOL);
      match_info = scm_cons(scm_double2num(rtop_info.dist_score), match_info);
      SCM s = scm_cons(match_info, SCM_EOL);
      scm_status = scm_cons(rtop_to_scm(rtop_info.rtop), s);
   }
   return scm_status;
}
#endif


#ifdef USE_GUILE
SCM
analyse_ligand_differences(int imol_ligand, int imol_ref, const char *chain_id_ref,
			   int resno_ref) {

   SCM scm_status = SCM_BOOL_F;
   coot::graph_match_info_t rtop_info =
      overlap_ligands_internal(imol_ligand, imol_ref, chain_id_ref, resno_ref, 0);

   std::cout << "analyse_ligand_differences: success       " << rtop_info.success << std::endl;
   std::cout << "analyse_ligand_differences: n_match       " << rtop_info.n_match << std::endl;
   std::cout << "analyse_ligand_differences: dist_score    " << rtop_info.dist_score << std::endl;
   std::cout << "analyse_ligand_differences: atoms matched " << rtop_info.matching_atom_names.size() << std::endl;
   std::cout << "analyse_ligand_differences: rtop: \n" << rtop_info.rtop.format() << std::endl;
   
   if (rtop_info.success) {
      SCM match_info = scm_cons(scm_int2num(rtop_info.n_match), SCM_EOL);
      match_info = scm_cons(scm_double2num(rtop_info.dist_score), match_info);
      SCM s = scm_cons(match_info, SCM_EOL);
      scm_status = scm_cons(rtop_to_scm(rtop_info.rtop), s);
   }
   return scm_status;
}
#endif   


coot::graph_match_info_t
overlap_ligands_internal(int imol_ligand, int imol_ref, const char *chain_id_ref,
			 int resno_ref, bool apply_rtop_flag) {

   coot::graph_match_info_t graph_info;
   
   int istat = 0;

   CResidue *residue_moving = 0;
   CResidue *residue_reference = 0;

   // running best ligands:
   CResidue *best_residue_moving = NULL;
   double best_score = 99999999.9; // low score good.
   clipper::RTop_orth best_rtop;
   

   if (! is_valid_model_molecule(imol_ligand))
      return graph_info;

   if (! is_valid_model_molecule(imol_ref))
      return graph_info;

   CMMDBManager *mol_moving = graphics_info_t::molecules[imol_ligand].atom_sel.mol;
   CMMDBManager *mol_ref    = graphics_info_t::molecules[imol_ref].atom_sel.mol;
   
   for (int imod=1; imod<=mol_moving->GetNumberOfModels(); imod++) { 
      CModel *model_p = mol_moving->GetModel(imod);
      CChain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 PCResidue residue_p;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    if (residue_p) { 
	       int n_atoms = residue_p->GetNumberOfAtoms();
	       if (n_atoms > 0) {
		  residue_moving = residue_p;
		  break;
	       }
	    }
	 }
	 if (residue_moving)
	    break;
      }

      if (! residue_moving) {
	 std::cout << "Oops.  Failed to find moving residue" << std::endl;
      } else { 
	 int imodel_ref = 1;
	 CModel *model_ref_p = mol_ref->GetModel(imodel_ref);
	 CChain *chain_p;
	 int nchains = model_ref_p->GetNumberOfChains();
	 for (int ichain=0; ichain<nchains; ichain++) {
	    chain_p = model_ref_p->GetChain(ichain);
	    if (std::string(chain_p->GetChainID()) == std::string(chain_id_ref)) { 
	       int nres = chain_p->GetNumberOfResidues();
	       PCResidue residue_p;
	       for (int ires=0; ires<nres; ires++) { 
		  residue_p = chain_p->GetResidue(ires);
		  if (residue_p) {
		     int seqnum = residue_p->GetSeqNum();
		     if (seqnum == resno_ref) {
			residue_reference = residue_p;
			break;
		     }
		  }
	       }
	       if (residue_reference)
		  break;
	    }
	 }

	 if (!residue_reference) {
	    std::cout << "Oops.  Failed to find reference residue" << std::endl;
	 } else { 
	    coot::graph_match_info_t rtop_info =
	       coot::graph_match(residue_moving, residue_reference, apply_rtop_flag);
	    if (rtop_info.success) {
	       if (rtop_info.dist_score < best_score) { // low score good.
		  best_score = rtop_info.dist_score;
		  best_residue_moving = residue_moving;
		  best_rtop = rtop_info.rtop;
		  graph_info = rtop_info;
	       }
	    } else {
	       // std::cout << "Oops.  Match failed somehow" << std::endl;
	    }
	 }
      }
   }
   if (apply_rtop_flag) { 
      if (best_residue_moving) {
	 // move just the best ligand:
	 graphics_info_t::molecules[imol_ligand].transform_by(best_rtop, best_residue_moving);
	 // delete everything except the best ligand:
	 graphics_info_t::molecules[imol_ligand].delete_all_except_res(best_residue_moving);
	 graphics_draw();
      }
   }
   return graph_info;
}



/*  ----------------------------------------------------------------------- */
/*                  ligand fitting stuff                                    */
/*  ----------------------------------------------------------------------- */

// We need to look up the vboxes and add items, like we did in code in
// gtk-manual.c
// 
// This is called by callback.c's find ligand button clicked callback.
// 
// Return the state.  If there is not ligands, protein or map, then
// return 0.  This is checked by the calling function (in callbacks.c)
// and the dialog is not displayed if 0 is returned.
// 
int fill_ligands_dialog(GtkWidget *find_ligand_dialog) {

   int ifound_map, ifound_coords, ifound_ligand; 
   short int diff_maps_only_flag = 0;
   ifound_map = fill_ligands_dialog_map_bits(find_ligand_dialog, diff_maps_only_flag);
   if (ifound_map == 0) {
      std::cout << "WARNING:: you must have a map to search for ligands!"
		<< std::endl;
      std::string s("WARNING:: you must have a map to\n search for ligands!");
      GtkWidget *w = wrapped_nothing_bad_dialog(s);
      gtk_widget_show(w);
   } 
   ifound_coords = fill_ligands_dialog_protein_bits(find_ligand_dialog);
   if (ifound_coords == 0) {
      std::cout << "Error: you must have a protein to mask the map!"
		<< std::endl;
      std::string s("WARNING:: you must have a protein\n to mask the map");
      GtkWidget *w = wrapped_nothing_bad_dialog(s);
      gtk_widget_show(w);
   } 
   ifound_ligand = fill_ligands_dialog_ligands_bits(find_ligand_dialog);
   if (ifound_ligand == 0) {
      std::cout << "Error: you must have at least one  ligand to search for!"
		<< std::endl;
      std::string s("WARNING:: you must have at least one\n          ligand to search for!\n");
      s += "         Ligands have less than ";
      s += graphics_info_t::int_to_string(graphics_info_t::find_ligand_ligand_atom_limit);
      s += " atoms\n";
      GtkWidget *w = wrapped_nothing_bad_dialog(s);
      gtk_widget_show(w);
   }

   // The mask waters toggle buttons:

   GtkWidget *togglebutton;
   togglebutton = lookup_widget(find_ligand_dialog,
				"find_ligand_mask_waters_yes_radiobutton");
   if (graphics_info_t::find_ligand_mask_waters_flag)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(togglebutton), TRUE);
   else 
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(togglebutton), FALSE);

      
   // 040211: fill new sigma level entry
   fill_ligands_sigma_level_entry(find_ligand_dialog);

   // expert options
   fill_ligands_expert_options(find_ligand_dialog);
   // shall we see the expert option frame?
   if (graphics_info_t::ligand_expert_flag == 0) {
      GtkWidget *frame = lookup_widget(find_ligand_dialog, "ligand_expert_frame");
      gtk_widget_hide(frame);
   }
      
   return ifound_ligand * ifound_map * ifound_coords;

   // 050924 New "Expert Options" entries:
   

}

// 110204: fill new sigma level entry
// 
void fill_ligands_sigma_level_entry(GtkWidget *dialog) {
   
   GtkWidget *entry = lookup_widget(dialog, "find_ligand_sigma_level_entry");
   gtk_entry_set_text(GTK_ENTRY(entry), graphics_info_t::float_to_string(graphics_info_t::ligand_cluster_sigma_level).c_str());
   
}

void fill_ligands_expert_options(GtkWidget *find_ligand_dialog) {

   GtkWidget *entry = lookup_widget(find_ligand_dialog, "ligand_n_samples_entry");
   gtk_entry_set_text(GTK_ENTRY(entry),
		      graphics_info_t::int_to_string(graphics_info_t::ligand_wiggly_ligand_n_samples).c_str());

   entry = lookup_widget(find_ligand_dialog, "ligand_n_top_ligands_entry");
   gtk_entry_set_text(GTK_ENTRY(entry),
		      graphics_info_t::int_to_string(graphics_info_t::find_ligand_n_top_ligands).c_str());

}


int fill_ligands_dialog_map_bits(GtkWidget *find_ligand_dialog, 
				 short int diff_maps_only_flag) {

   return fill_ligands_dialog_map_bits_by_dialog_name(find_ligand_dialog,
						      "find_ligand_map",
						      diff_maps_only_flag);

}

// dialog_name is typically find_ligand_map
// 
int fill_ligands_dialog_map_bits_by_dialog_name(GtkWidget *find_ligand_dialog,
						const char *dialog_name, 
						short int diff_maps_only_flag) {

   int ifound = 0;
   graphics_info_t g;
   // Add map elements:
   GSList *find_ligand_map_group = NULL;
   //
   std::string vbox_name = dialog_name;
   vbox_name += "_vbox";
   
   GtkWidget *find_ligand_map_vbox =
      lookup_widget(find_ligand_dialog,vbox_name.c_str());
   if (find_ligand_map_vbox == NULL) {
      std::cout << "disaster! find_ligand map vbox not found " << std::endl;
   } else { 
      for (int imol=0; imol<g.n_molecules; imol++) {
	 if (g.molecules[imol].has_map()) {

	    if ((!diff_maps_only_flag) ||
		g.molecules[imol].is_difference_map_p()) {

	       ifound++; // there was a map
	       std::string map_str(dialog_name);
	       map_str += "_radiobutton_";
	       map_str += g.int_to_string(imol);
	       std::string map_button_label = g.int_to_string(imol);
	       map_button_label += " ";
	       map_button_label += g.molecules[imol].name_;

	       GtkWidget *find_ligand_map_radiobutton_imol =
		  gtk_radio_button_new_with_label (find_ligand_map_group,
						   map_button_label.c_str());
	       find_ligand_map_group =
		  gtk_radio_button_group (GTK_RADIO_BUTTON (find_ligand_map_radiobutton_imol));
	       gtk_widget_ref (find_ligand_map_radiobutton_imol);
	       gtk_object_set_data_full (GTK_OBJECT (find_ligand_dialog),
					 map_str.c_str(),
					 find_ligand_map_radiobutton_imol,
					 (GtkDestroyNotify) gtk_widget_unref);
	       gtk_widget_show (find_ligand_map_radiobutton_imol);
	       gtk_box_pack_start (GTK_BOX (find_ligand_map_vbox),
				   find_ligand_map_radiobutton_imol, FALSE, FALSE, 0);
	    }
	 }
      }
   }
   return ifound; 
}

int fill_ligands_dialog_protein_bits(GtkWidget *find_ligand_dialog) {
   
   return fill_ligands_dialog_protein_bits_by_dialog_name(find_ligand_dialog,
							  "find_ligand_protein");
   
}


int fill_ligands_dialog_protein_bits_by_dialog_name(GtkWidget *find_ligand_dialog,
						    const char *dialog_name) {
   
   
   int ifound = 0;
   graphics_info_t g;
   // Add protein elements:
   GSList *find_ligand_protein_group = NULL;
   //
   std::string vbox_name(dialog_name);
   vbox_name += "_vbox";
   
   GtkWidget *find_ligand_protein_vbox =
      lookup_widget(find_ligand_dialog,vbox_name.c_str());
   if (find_ligand_protein_vbox == NULL) {
      std::cout << "disaster! find_ligand protein vbox not found " << std::endl;
   } else { 
      for (int imol=0; imol<g.n_molecules; imol++) {
	 if (g.molecules[imol].atom_sel.n_selected_atoms) {
	    //
	    ifound = 1; // there was a protein
	    std::string protein_str(dialog_name);
	    protein_str += "_radiobutton_";
	    protein_str += g.int_to_string(imol);
	    std::string protein_button_label = g.int_to_string(imol);
	    protein_button_label += " ";
	    protein_button_label += g.molecules[imol].name_;
	    GtkWidget *find_ligand_protein_radiobutton_imol =
	       gtk_radio_button_new_with_label (find_ligand_protein_group,
						protein_button_label.c_str());
	    find_ligand_protein_group =
	       gtk_radio_button_group (GTK_RADIO_BUTTON
				       (find_ligand_protein_radiobutton_imol));
	    gtk_widget_ref (find_ligand_protein_radiobutton_imol);
	    gtk_object_set_data_full (GTK_OBJECT (find_ligand_dialog),
				      protein_str.c_str(),
				      find_ligand_protein_radiobutton_imol,
				      (GtkDestroyNotify) gtk_widget_unref);
	    gtk_widget_show (find_ligand_protein_radiobutton_imol);
	    gtk_box_pack_start (GTK_BOX (find_ligand_protein_vbox),
				find_ligand_protein_radiobutton_imol, FALSE, FALSE, 0);
	 }
      }
   }
   return ifound; 
}

int fill_vbox_with_coords_options_by_dialog_name(GtkWidget *find_ligand_dialog,
						 const char *dialog_name,
						 short int have_ncs_flag) {
   
   
   int ifound = 0;
   graphics_info_t g;
   // Add protein elements:
   GSList *find_ligand_protein_group = NULL;
   //
   std::string vbox_name(dialog_name);
   vbox_name += "_vbox";
   
   GtkWidget *find_ligand_protein_vbox =
      lookup_widget(find_ligand_dialog,vbox_name.c_str());
   if (find_ligand_protein_vbox == NULL) {
      std::cout << "disaster! fill_vbox_with_coords_options_by_dialog_name coords"
		<< " vbox not found " << std::endl;
   } else { 
      for (int imol=0; imol<g.n_molecules; imol++) {
	 if (g.molecules[imol].has_model()) {

	    if (!have_ncs_flag || g.molecules[imol].has_ncs_p()) { 
	       //
	       ifound = 1; // there was a protein
	       std::string protein_str(dialog_name);
	       protein_str += "_radiobutton_";
	       protein_str += g.int_to_string(imol);
	       std::string protein_button_label = g.int_to_string(imol);
	       protein_button_label += " ";
	       protein_button_label += g.molecules[imol].name_;
	       GtkWidget *find_ligand_protein_radiobutton_imol =
		  gtk_radio_button_new_with_label (find_ligand_protein_group,
						   protein_button_label.c_str());
	       find_ligand_protein_group =
		  gtk_radio_button_group (GTK_RADIO_BUTTON
					  (find_ligand_protein_radiobutton_imol));
	       gtk_widget_ref (find_ligand_protein_radiobutton_imol);
	       gtk_object_set_data_full (GTK_OBJECT (find_ligand_dialog),
					 protein_str.c_str(),
					 find_ligand_protein_radiobutton_imol,
					 (GtkDestroyNotify) gtk_widget_unref);
	       gtk_widget_show (find_ligand_protein_radiobutton_imol);
	       gtk_box_pack_start (GTK_BOX (find_ligand_protein_vbox),
				   find_ligand_protein_radiobutton_imol, FALSE, FALSE, 0);
	    }
	 }
      }
   }
   return ifound; 
}


int fill_ligands_dialog_ligands_bits(GtkWidget *find_ligand_dialog) {

   int ifound = 0;
   graphics_info_t g;
   // Add ligand elements
   //
   GtkWidget *hbox;

   GtkWidget *find_ligand_ligands_vbox =
      lookup_widget(find_ligand_dialog, "find_ligand_ligands_vbox");
   if (find_ligand_ligands_vbox == NULL) {
      std::cout << "disaster! find_ligand protein vbox not found " << std::endl;
   } else { 
      for (int imol=0; imol<g.n_molecules; imol++) {
	 if ((g.molecules[imol].atom_sel.n_selected_atoms < graphics_info_t::find_ligand_ligand_atom_limit) &&
	     g.molecules[imol].has_model()) {
	    ifound = 1; // there was a ligand

	    // create an hbox:
	    hbox = gtk_hbox_new (FALSE, 0);
	    
	    std::string ligands_str("find_ligand_ligand_checkbutton_");
	    ligands_str += g.int_to_string(imol);
	    std::string ligands_button_label = g.int_to_string(imol);
	    ligands_button_label += " ";
	    ligands_button_label += g.molecules[imol].name_;

	    // ligand on/off check button
	    // 
	    GtkWidget *find_ligand_ligands_checkbutton_imol =
	       gtk_check_button_new_with_label (ligands_button_label.c_str());
	    gtk_widget_ref (find_ligand_ligands_checkbutton_imol);
	    gtk_object_set_data_full (GTK_OBJECT (find_ligand_dialog),
				      ligands_str.c_str(),
				      find_ligand_ligands_checkbutton_imol,
				      (GtkDestroyNotify) gtk_widget_unref);

	    gtk_widget_show (find_ligand_ligands_checkbutton_imol);
	    gtk_box_pack_start (GTK_BOX (hbox),
				find_ligand_ligands_checkbutton_imol, FALSE, FALSE, 0);

	    // wligand on/off check button
	    //
	    std::string wligands_str("find_ligand_wligand_checkbutton_");
	    wligands_str += g.int_to_string(imol);
	    GtkWidget *find_ligand_wligands_checkbutton_imol =
	       gtk_check_button_new_with_label ("Flexible?");
	    gtk_widget_ref (find_ligand_wligands_checkbutton_imol);
	    gtk_object_set_data_full (GTK_OBJECT (find_ligand_dialog),
				      wligands_str.c_str(),
				      find_ligand_wligands_checkbutton_imol,
				      (GtkDestroyNotify) gtk_widget_unref);

	    gtk_widget_show (find_ligand_wligands_checkbutton_imol);
	    gtk_box_pack_start (GTK_BOX (hbox),
				find_ligand_wligands_checkbutton_imol, FALSE, FALSE, 0);

	    // pack the hbox into the vbox
	    gtk_box_pack_start (GTK_BOX (find_ligand_ligands_vbox),
				hbox, FALSE, FALSE, 0);
	    gtk_widget_show(hbox);
	 }
      }
   }
   return ifound; 
}


/* get which map to search, protein mask and ligands from button and
   then do it*/
// Recall that the map and the protein are radio buttons (i.e. we want
// exactly one of each).  Hoewever, the ligands are check buttons, we
// can have as many of those on that we like.
//
// We should add a check to the install_ligand that checks the number
// of atoms in the ligand and gives a warning... hmmm... that invovles
// another popup-sigh...  Oh well, plan for it, even if you don't
// implement it now - i.e. invoke a post-check function that creates
// the coot::ligand object.
//
void execute_get_mols_ligand_search(GtkWidget *button) {

   graphics_info_t g;
   GtkWidget *ligand_button;
   std::vector<int> chief_ligand_many_atoms; // caches imols with lots
					     // of atoms.
   std::vector<std::pair<int, bool> > wiggly_ligand_info; 

   // extract the sigma level and stick it in
   // graphics_info_t::ligand_cluster_sigma_level
   // 
   set_ligand_cluster_sigma_level_from_widget(button);
   set_ligand_expert_options_from_widget(button); // ligand_n_top_sites
						  // and wiggly ligand
						  // n (conformer)
						  // samples
   
   // flags to say if we have found things:
   short int found_active_button_for_map = 0;
   short int found_active_button_for_protein = 0;
   short int found_active_button_for_ligands = 0;
   // This is where we store the mols:
   int find_ligand_map_mol = -1; // get set?
   int find_ligand_protein_mol;
   // std::vector<int> find_ligand_ligand_mols;  old.
   
   // Find the first active map radiobutton
   found_active_button_for_map = 0;
   for (int imol=0; imol<g.n_molecules; imol++) {
      if (g.molecules[imol].xmap_is_filled[0]) { 
	 std::string map_str = "find_ligand_map_radiobutton_";
	 map_str += g.int_to_string(imol);
	 ligand_button = lookup_widget(button, map_str.c_str());
	 if (ligand_button) { 
	    if (GTK_TOGGLE_BUTTON(ligand_button)->active) {
	       find_ligand_map_mol = imol;
	       found_active_button_for_map = 1;
	       break;
	    }
	 } else {
	    std::cout << map_str << " widget not found in "
		      << "execute_get_mols_ligand_search" << std::endl;
	 }
      }
   }

   // Find the first active protein radiobutton
   found_active_button_for_protein = 0;
   for (int imol=0; imol<g.n_molecules; imol++) {
      if (g.molecules[imol].atom_sel.n_selected_atoms > 0) { 
	 std::string protein_str = "find_ligand_protein_radiobutton_";
	 protein_str += g.int_to_string(imol);
	 ligand_button = lookup_widget(button, protein_str.c_str());
	 if (ligand_button) { 
	    if (GTK_TOGGLE_BUTTON(ligand_button)->active) {
	       find_ligand_protein_mol = imol;
	       found_active_button_for_protein = 1;
	       break;
	    }
	 } else {
	    std::cout << protein_str << " widget not found in "
		      << "execute_get_mols_ligand_search" << std::endl;
	 }
      }
   }

   // Now, do we mask waters for the protein mask?
   
   GtkWidget *togglebutton;
   togglebutton = lookup_widget(button,
				"find_ligand_mask_waters_yes_radiobutton");
   if (GTK_TOGGLE_BUTTON(togglebutton)->active)
      graphics_info_t::find_ligand_mask_waters_flag = 1;
   else 
      graphics_info_t::find_ligand_mask_waters_flag = 0;

   // For each imol in molecules, if we have selected coordinates,
   // construct a string begining "find_ligand_ligands_imol_" then the
   // molecule number
   
   found_active_button_for_ligands = 0;
   for (int imol=0; imol<g.n_molecules; imol++) {
      if (g.molecules[imol].has_model() &&
	  g.molecules[imol].atom_sel.n_selected_atoms < graphics_info_t::find_ligand_ligand_atom_limit) {
	 std::string ligand_str = "find_ligand_ligand_checkbutton_";
	 ligand_str += g.int_to_string(imol);
	 ligand_button = lookup_widget(button, ligand_str.c_str());

	 std::string wiggly_str = "find_ligand_wligand_checkbutton_";
	 wiggly_str += g.int_to_string(imol);
	 GtkWidget *wiggly_button = lookup_widget(button, wiggly_str.c_str());
		  
	 if (ligand_button && wiggly_button) { 
	    if (GTK_TOGGLE_BUTTON(ligand_button)->active) {
	       
	       bool wiggly_state = 0;
	       if (GTK_TOGGLE_BUTTON(wiggly_button)->active)
		  wiggly_state = 1;
	       wiggly_ligand_info.push_back(std::pair<int, bool> (imol, wiggly_state));
// 	       std::cout << "DEBUG:: wiggly info: " << imol <<  " " << wiggly_state
// 			 << " pushed back" << std::endl;
	       found_active_button_for_ligands = 1;
	       int n_atoms = g.molecules[imol].atom_sel.n_selected_atoms;
	       if (n_atoms > 100) {
		  std::cout << "WARNING:: molecule " << imol
			    << " has unexpectedly many atoms ("
			    << n_atoms << ")" << std::endl;
		  chief_ligand_many_atoms.push_back(imol);
	       } 
// 	    } else {
// 	       std::cout << "DEBUG:: button " << ligand_str
// 			 << " was not active" << std::endl;
	    }
	 } else {
	    std::cout << ligand_str << " widget not found in "
		      << "execute_get_mols_ligand_search" << std::endl;
	 }
      }
   }

   if ( found_active_button_for_map &&
	found_active_button_for_protein &&
	found_active_button_for_ligands) {

      // So store the clicked info in a graphics_info_t static:
      // because we need this info when we click the OK button of the
      // create_find_ligand_many_atoms_dialog() widget.  We don't want
      // to mess with set_user_data for many data.
      // 
      graphics_info_t g;
      g.set_find_ligands_mols(find_ligand_map_mol,
			      find_ligand_protein_mol,
			      wiggly_ligand_info);

      if (chief_ligand_many_atoms.size() == 0 ) {
	 execute_ligand_search();
      } else {
	 // we need to delete this widget when OK and cancel of the
	 // many atoms widget is pressed, so let's attach it as user
	 // data.
	 GtkWidget *widget = lookup_widget(button, "find_ligand_dialog");
	 do_find_ligand_many_atoms_in_ligands(widget);
      }
   } else {
	 std::cout << "Something wrong in the selection of map/molecules"
		   << std::endl;
   }

}

// q_ligands is questionable ligands (i.e. very large)
//   
void do_find_ligand_many_atoms_in_ligands(GtkWidget *find_ligand_dialog) {

   GtkWidget *widget = create_find_ligand_many_atoms_dialog();
   gtk_object_set_user_data(GTK_OBJECT(widget), (char *) find_ligand_dialog);
   gtk_widget_show(widget); 
}

void set_ligand_expert_options_from_widget(GtkWidget *button) {

   GtkWidget *entry = lookup_widget(button, "ligand_n_samples_entry");
   const gchar *text = gtk_entry_get_text(GTK_ENTRY(entry));
   if (text) {
      int isample = atoi(text);
      if ((isample > 0) && (isample < 1000000))
	 graphics_info_t::ligand_wiggly_ligand_n_samples = isample;
   }
   entry = lookup_widget(button, "ligand_n_top_ligands_entry");
   text = gtk_entry_get_text(GTK_ENTRY(entry));
   if (text) {
      int itop = atoi(text);
      if ((itop > 0) && (itop < 1000000))
	 graphics_info_t::find_ligand_n_top_ligands = itop;
   }
} 


/*  extract the sigma level and stick it in */
/*  graphics_info_t::ligand_cluster_sigma_level */
void set_ligand_cluster_sigma_level_from_widget(GtkWidget *button) { 
   
   GtkWidget *entry = lookup_widget(button, "find_ligand_sigma_level_entry");
   short int setit = 0;

   if (entry) {
      const gchar *text = gtk_entry_get_text(GTK_ENTRY(entry));
      if (text) { 
	 float f = atof(text);
	 if (f > 0.0 && f < 1000.0) { 
	    graphics_info_t::ligand_cluster_sigma_level = f;
	    setit = 1;
	 } 
      }
   }
   if (setit == 0) 
      std::cout << "INFO:: ignoring bogus attempt to set "
		<< "the ligand search sigma level" 
		<< std::endl;
} 


void ligand_expert() {  /* sets the flag to have expert option ligand entries in the GUI */

   graphics_info_t::ligand_expert_flag = 1;
}

/*! \brief set the protein molecule for ligand searching */
void set_ligand_search_protein_molecule(int imol) {
   graphics_info_t::set_ligand_protein_mol(imol);
   
}

/*! \brief set the map molecule for ligand searching */
void set_ligand_search_map_molecule(int imol_map) {
   graphics_info_t::set_ligand_map_mol(imol_map);
}

/*! \brief add a ligand molecule to the list of ligands to search for
  in ligand searching */
void add_ligand_search_wiggly_ligand_molecule(int imol_ligand) {
   graphics_info_t::find_ligand_add_flexible_ligand(imol_ligand);
}

/*! \brief add a ligand molecule to the list of ligands to search for
  in ligand searching */
void add_ligand_search_ligand_molecule(int imol_ligand) {
   if (is_valid_model_molecule(imol_ligand))
      graphics_info_t::find_ligand_add_rigid_ligand(imol_ligand);

   graphics_info_t g;
//    std::cout << "DEBUG:: graphics_info_t::find_ligand_wiggly_ligands()["
// 	     << imol_ligand << "] is ("
// 	     << g.find_ligand_ligand_mols()[imol_ligand].first  << ", " 
// 	     << g.find_ligand_ligand_mols()[imol_ligand].second << ")" << std::endl;

}

void add_ligand_clear_ligands() {

   graphics_info_t g;
   g.find_ligand_clear_ligand_mols();
} 




// execute_find_ligands_real, you might say
// 
#ifdef USE_GUILE
SCM execute_ligand_search() {

   std::vector<int> solutions = execute_ligand_search_internal();
   return generic_int_vector_to_list_internal(solutions);
}
#else
// Fixme Bernhard
void execute_ligand_search() {

   std::vector<int> solutions = execute_ligand_search_internal();
}
#endif 

std::vector<int>
execute_ligand_search_internal() {
   
   std::cout << "Executing ligand search..." << std::endl;
   std::vector<int> solutions;

   graphics_info_t g;

   if (! is_valid_model_molecule(g.find_ligand_protein_mol())) {
      std::cout << "Protein molecule for ligand search not set" << std::endl;
      std::cout << "Aborting ligand search" << std::endl;
      return solutions; 
   }
   if (! is_valid_map_molecule(g.find_ligand_map_mol())) {
      std::cout << "Map molecule for ligand search not set" << std::endl;
      std::cout << "Aborting ligand search" << std::endl;
      return solutions; 
   }
   if (g.find_ligand_ligand_mols().size() == 0) {
      std::cout << "No defined ligand molecules" << std::endl;
      std::cout << "Aborting ligand search" << std::endl;
      return solutions; 
   } 
   
   CMMDBManager *protein_mol = 
      g.molecules[g.find_ligand_protein_mol()].atom_sel.mol;

   coot::wligand wlig;
   if (g.ligand_verbose_reporting_flag)
      wlig.set_verbose_reporting();
   wlig.import_map_from(g.molecules[g.find_ligand_map_mol()].xmap_list[0]);
   std::vector<std::pair<int, bool> > ligands = g.find_ligand_ligand_mols();
   for(unsigned int i=0; i<ligands.size(); i++) {

      std::cout << "ligand number " << i << " is molecule number "
		<< g.find_ligand_ligand_mols()[i].first << "  " 
		<< " with wiggly flag: "
		<< g.find_ligand_ligand_mols()[i].second << std::endl;

      if (ligands[i].second) {
	 // argh (i).
	 coot::minimol::molecule mmol(g.molecules[ligands[i].first].atom_sel.mol);

	 for(unsigned int ifrag=0; ifrag<mmol.fragments.size(); ifrag++) {
	    for (int ires=mmol[ifrag].min_res_no(); ires<=mmol[ifrag].max_residue_number();
		 ires++) {
// 	       if (mmol[ifrag][ires].n_atoms() > 0) {
// 		  std::cout << "DEBUG:: in execute_ligand_search:  mmol["
// 			    << ifrag << "][" << ires << "].name :"
// 			    <<  mmol[ifrag][ires].name << ":" << std::endl;
// 	       }
	    }
	 }

	 std::pair<short int, std::string> istat_pair =
	    wlig.install_simple_wiggly_ligands(g.Geom_p(), mmol,
					       g.ligand_wiggly_ligand_n_samples);
	 if (! istat_pair.first) {
	    // if any fail, it all fails.
	    std::cout << "Error in flexible ligand definition.\n";
	    GtkWidget *w = wrapped_nothing_bad_dialog(istat_pair.second);
	    gtk_widget_show(w);
	    return solutions;
	 }
      } else { 
	 // argh (ii).
	 wlig.install_ligand(g.molecules[ligands[i].first].atom_sel.mol);
      }
   }

   int imol = g.n_molecules;
   
   short int mask_waters_flag; // treat waters like other atoms?
   mask_waters_flag = g.find_ligand_mask_waters_flag;
   if (graphics_info_t::map_mask_atom_radius > 0) {
      // only do this if it was set by the user.
      wlig.set_map_atom_mask_radius(graphics_info_t::map_mask_atom_radius);
   } else { 
      wlig.set_map_atom_mask_radius(2.0);  // Angstroms
   } 

   std::string name("ligand masked map");
   // std::cout << "DEBUG:: calling mask_map\n";
   wlig.mask_map(protein_mol, mask_waters_flag); // mask by protein
   // std::cout << "DEBUG:: done mask_map\n";
   g.expand_molecule_space_maybe(); 
   g.molecules[imol].new_map(wlig.masked_map(), wlig.masked_map_name());
   wlig.set_acceptable_fit_fraction(g.ligand_acceptable_fit_fraction);
   wlig.find_clusters(g.ligand_cluster_sigma_level);  // trashes the xmap
   wlig.fit_ligands_to_clusters(g.find_ligand_n_top_ligands); // 10 clusters
   wlig.make_pseudo_atoms(); // put anisotropic atoms at the ligand sites

   g.scroll_wheel_map = imol;  // change the current scrollable map to
			       // the masked map.

   g.n_molecules++; 

   // now add in the solution ligands:
   int n_final_ligands = wlig.n_final_ligands(); 
   int n_new_ligand = 0;
   coot::minimol::molecule m;
   for (int ilig=0; ilig<n_final_ligands; ilig++) { 
      m = wlig.get_solution(ilig);
      if (! m.is_empty()) {
	 float bf = graphics_info_t::default_new_atoms_b_factor;
	 atom_selection_container_t asc = make_asc(m.pcmmdbmanager(bf));
	 int g_mol = g.n_molecules;
	 std::string label = "Fitted ligand #";
	 label += g.int_to_string(ilig);
	 g.expand_molecule_space_maybe(); 
	 g.molecules[g_mol].install_model(asc, label, 1);
	 solutions.push_back(g_mol);
	 g.n_molecules++;
	 n_new_ligand++;
	 if (g.go_to_atom_window){
	    g.update_go_to_atom_window_on_new_mol();
	    g.update_go_to_atom_window_on_changed_mol(g_mol);
	 }
      }
   }

   // We need some python code here to match post-ligand-fit-gui
#ifdef USE_GUILE
      safe_scheme_command("(post-ligand-fit-gui)");
#endif // USE_GUILE

   if (n_new_ligand) { 
      GtkWidget *w = create_new_ligands_info_dialog();
      GtkWidget *label = lookup_widget(w, "new_ligands_info_dialog_label");
      std::string label_str("  Found ");
      label_str += graphics_info_t::int_to_string(n_new_ligand);
      if (n_new_ligand == 1) 
	 label_str += " acceptable ligand  ";
      else 
	 label_str += " acceptable ligands  ";
      gtk_label_set_text(GTK_LABEL(label), label_str.c_str());
      gtk_widget_show(w);
   } else { 
      GtkWidget *w = create_no_new_ligands_info_dialog();
      gtk_widget_show(w);
   }

   graphics_draw();
   return solutions;
}

void set_ligand_cluster_sigma_level(float f) { /* default 2.2 */
   graphics_info_t::ligand_cluster_sigma_level = f;
}

void set_find_ligand_mask_waters(int i) {

   graphics_info_t::find_ligand_mask_waters_flag = i;

}

/*! \brief set the atom radius for map masking */
void set_map_mask_atom_radius(float rad) {
   graphics_info_t::map_mask_atom_radius = rad; 

}



void set_ligand_flexible_ligand_n_samples(int i) {
   graphics_info_t::ligand_wiggly_ligand_n_samples = i;
} 


void set_ligand_acceptable_fit_fraction(float f) {
   if (f >= 0.0 && f<= 1.0) {
      graphics_info_t::ligand_acceptable_fit_fraction = f;
      // std::cout << "ligand acceptable fit fraction set to " << f << "\n";
   } else {
      // std::cout << "ligand acceptable fit fraction " << f << " ignored\n";
   }
}

void set_find_ligand_n_top_ligands(int n) { /* fit the top n ligands,
					      not all of them, default
					      10. */
   graphics_info_t g;
   g.find_ligand_n_top_ligands = n;

}


/*  ----------------------------------------------------------------------- */
/*                  Mask                                                    */
/*  ----------------------------------------------------------------------- */
/* The idea is to generate a new map that has been masked by some
   coordinates. */
/*!     (mask-map-by-molecule map-no mol-no invert?)  creates and
        displays a masked map, cuts down density where the coordinates
        are not.  If invert? is #t, cut the density down where the
        atoms are.  */
int mask_map_by_molecule(int map_mol_no, int coord_mol_no, short int invert_flag) {

   int imol_new_map = -1;

   // create a new map
   //
   // c.f. a new map that is created "Masked by Protein"
   // we use molecule_class_info_t's new_map()
   //

   // Where should the bulk of this function be?  Let's put it here for now.
   coot::ligand lig;
   graphics_info_t g;
   if (map_mol_no >= g.n_molecules) {
      std::cout << "No such molecule (no map) at molecule number " << map_mol_no << std::endl;
   } else {

      if (coord_mol_no >= g.n_molecules) {
	 std::cout << "No such molecule (no coords) at molecule number " << map_mol_no << std::endl;
      } else { 
      
	 if (!g.molecules[map_mol_no].has_map()) {
	    std::cout << "No map in molecule number " << map_mol_no << std::endl;
	 } else {
	    
	    if (! g.molecules[coord_mol_no].has_model()) {
	       std::cout << "No model in molecule number " << map_mol_no << std::endl;
	    } else {
	       short int mask_waters_flag; // treat the waters like protein atoms?
	       mask_waters_flag = graphics_info_t::find_ligand_mask_waters_flag;
	       lig.import_map_from(g.molecules[map_mol_no].xmap_list[0]);
	       int selectionhandle = g.molecules[coord_mol_no].atom_sel.mol->NewSelection();

	       if (graphics_info_t::map_mask_atom_radius > 0) {
		  lig.set_map_atom_mask_radius(graphics_info_t::map_mask_atom_radius);
	       }

	       // make a selection:
	       std::string rnames = "*";
	       if (!mask_waters_flag)
		  rnames = "!HOH,WAT"; // treat waters differently to regular atoms.
	       g.molecules[coord_mol_no].atom_sel.mol->SelectAtoms(selectionhandle, 0, "*",
								   ANY_RES, "*",
								   ANY_RES, "*",
								   (char *) rnames.c_str(),
								   "*", "*", "*");
	       
	       lig.mask_map(g.molecules[coord_mol_no].atom_sel.mol, selectionhandle, invert_flag);
	       g.molecules[coord_mol_no].atom_sel.mol->DeleteSelection(selectionhandle);
	       std::cout << "INFO:: Creating masked  map in molecule number "
			 << g.n_molecules << std::endl;
	       g.molecules[g.n_molecules].new_map(lig.masked_map(), "Generic Masked Map");
	       imol_new_map = g.n_molecules; 
	       g.n_molecules++;
	       graphics_draw();
	    }
	 }
      }
   }
   return imol_new_map; 
}

int
mask_map_by_atom_selection(int map_mol_no, int coords_mol_no, const char *mmdb_atom_selection, short int invert_flag) {

   int imol_new_map = -1;
   graphics_info_t g;
   if (is_valid_map_molecule(map_mol_no)) {
      if (is_valid_model_molecule(coords_mol_no)) {
	 std::cout << "making lig..." << std::endl;
	 coot::ligand lig;
	 short int mask_waters_flag; // treat the waters like protein atoms?
	 lig.import_map_from(g.molecules[map_mol_no].xmap_list[0]);

	 if (graphics_info_t::map_mask_atom_radius > 0) {
	    lig.set_map_atom_mask_radius(graphics_info_t::map_mask_atom_radius);
	 }
	 int selectionhandle = g.molecules[coords_mol_no].atom_sel.mol->NewSelection();
	 g.molecules[coords_mol_no].atom_sel.mol->Select(selectionhandle, STYPE_ATOM,
							 (char *) mmdb_atom_selection,
							 SKEY_NEW);
	 lig.mask_map(g.molecules[coords_mol_no].atom_sel.mol, selectionhandle, invert_flag);
	 g.molecules[g.n_molecules].new_map(lig.masked_map(), "Generic Masked Map");
	 imol_new_map = g.n_molecules; 
	 g.n_molecules++;
	 graphics_draw();
      } else {
	 std::cout << "No model molecule in " << coords_mol_no << std::endl;
      }
   } else {
      std::cout << "No map molecule in " << map_mol_no << std::endl;
   } 
   return imol_new_map; 
}




void do_find_ligands_dialog() {
   GtkWidget *dialog;
   int istate;
   dialog = create_find_ligand_dialog();
   istate = fill_ligands_dialog(dialog); /* return OK, we have map(s), ligand(s), masking(s) */
   if (istate == 0) {
      gtk_widget_destroy(dialog);
      std::string s("Problem finding maps, coords or ligands!");
      graphics_info_t g;
      g.statusbar_text(s);
      std::cout << s << std::endl;
   }
   else 
     gtk_widget_show(dialog); 

}

/*  ----------------------------------------------------------------------- */
/*                  monomer lib                                             */
/*  ----------------------------------------------------------------------- */
std::vector<std::pair<std::string, std::string> > monomer_lib_3_letter_codes_matching(const std::string &search_string, short int allow_minimal_descriptions_flag) {

   graphics_info_t g;
   std::vector<std::pair<std::string, std::string > > v = g.Geom_p()->matching_names(search_string, allow_minimal_descriptions_flag);

   // search the list of monomers for text search_text:

   return v;
} 


// we allocate new memory here, without ever giving it back.  The
// memory should be freed when the dialog is destroyed.
// 
int
handle_make_monomer_search(const char *text, GtkWidget *viewport) {

   int stat = 0;
   std::string t(text);

   // std::cout << "DEBUG:: handle_make_monomer_search " << text << std::endl;

   GtkWidget *vbox_current = lookup_widget(viewport, "monomer_search_results_vbox");
   GtkWidget *checkbutton =
      lookup_widget(viewport, "monomer_search_minimal_descriptions_checkbutton");
   short int allow_minimal_descriptions_flag = 0;
   GtkWidget *dialog = lookup_widget(viewport, "monomer_search_dialog");

   if (GTK_TOGGLE_BUTTON(checkbutton)->active)
      allow_minimal_descriptions_flag = 1;
   
   std::vector<std::pair<std::string, std::string> > v =
      monomer_lib_3_letter_codes_matching(t, allow_minimal_descriptions_flag);

   // std::cout << "DEBUG:: " << v.size() << " solutions matching" << std::endl;

   // here clear the current contents of the monomer vbox:
   // delete the user_data assocated with the buttons too.
    GList *children = gtk_container_children(GTK_CONTAINER(vbox_current));
    int nchild = 0; 
    while (children) {
       // std::cout << "child " << nchild << "  " << (GtkWidget *) children->data << std::endl;
       gtk_widget_destroy((GtkWidget *) children->data); 
       nchild++;
       children = g_list_remove_link(children, children);
       
    }

   GtkWidget *vbox = vbox_current;

   // std::cout << "DEBUG:: monomers v.size() " << v.size() << std::endl;
   // add new buttons
   for (unsigned int i=0; i<v.size(); i++) {
      // std::cout << i << " " << v[i].first << std::endl;
      std::string l = v[i].first;
      l += " : ";
      l += v[i].second;
      // std::cout << "Giving the button the label :" << l << ":" << std::endl;
      GtkWidget *button = gtk_button_new_with_label(l.c_str());
      // std::cout << "Adding button: " << button << std::endl;
      std::string button_name = "monomer_button_";
      string *s = new string(v[i].first); // the 3-letter-code/comp_id (for user data).
      button_name += v[i].first;
      gtk_widget_ref (button);
      gtk_object_set_data_full (GTK_OBJECT (dialog), 
				button_name.c_str(), button,
				(GtkDestroyNotify) gtk_widget_unref);
      gtk_box_pack_start (GTK_BOX (vbox), button, FALSE, FALSE, 0);
      gtk_container_set_border_width (GTK_CONTAINER (button), 2);

      gtk_signal_connect(GTK_OBJECT(button), "clicked",
			 GTK_SIGNAL_FUNC (on_monomer_lib_search_results_button_press),
			 s);
      
      gtk_widget_show(button);
   }

   int new_box_size = v.size() * 28 + 120; // plus extra for static parts.
   if (new_box_size > 520)
      new_box_size = 520;
   
   // we need to set widget size to new_box_size.  On the dialog?

   gtk_widget_set_usize(dialog, dialog->allocation.width, new_box_size);
      
   // a box of 14 is 400 pixels.  400 is about max size, I'd say 
   gtk_signal_emit_by_name(GTK_OBJECT(vbox), "check_resize");
   gtk_widget_show (vbox);
   return stat;

} 

void
on_monomer_lib_search_results_button_press (GtkButton *button,
					    gpointer user_data) {

   std::string *s = (std::string *) user_data;
   get_monomer(s->c_str());

} 

// Sigh.  This is wrong. We don't want to replace the coordinates - we
// do that in molecule_class_info_t::replace_coords()
//
// But this may be useful in future
// (not not yet compile tested)
// 
// atom_selection_container_t rigid_body_asc;

//       std::vector<CAtom *> atoms;
//       CAtom *at;

//       for (int n=0; n<g.molecules[imol].atom_sel.n_selected_atoms; n++) {
// 	 at = g.molecules[imol].atom_sel.atom_selection[n];
// 	 std::string mmdb_chain(at->residue-GetChainID());
// 	 for(int ifrag=0; ifrag<mol.fragments.size(); ifrag++) {
// 	    for (int ires=1; ires<=mol[ifrag].n_residues(); ires++) {
// 	       for (int iatom=0; iatom<mol[ifrag][ires].atoms.size(); iatom++) {
// 		  if (mmdb_chain == mol[ifrag].fragment_id) {
// 		     if (ires == at->residue->seqNum) {
// 			std::string mmdb_name(at->name);
// 			if (mmdb_name == mol[ifrag][ires][iatom].name) {
// 			   atoms.push_back(at);
// 			}
// 		     }
// 		  }
// 	       }
// 	    }
// 	 }
//       }

//       PPCAtom atom_selection = new CAtom *[atoms.size()];
//       for (int i=0; i<atoms.size(); i++) 
// 	 atom_selection[i] = atoms[i];

//       rigid_body_asc.atom_selection = atom_selection;
//       rigid_body_asc.n_selected_atoms = atoms.size();
//       std::cout << "INFO:: make_atom_selection found " << atoms.size()
// 		<< " atoms" << std::endl;

//       rigid_body_asc.mol = g.molecules[g.imol_rigid_body_refine].atom_sel.mol;


// c.f. molecule_class_info_t::atom_index
//
// This is not a member class of molecule_class_info_t because we dont
// want minimol dependency in molecule_class_info_t and because we
// want the number of atoms *and* the pointer to be returned... messy
// messy - lets do that here.
//
// We only use the n_selected_atoms and atom_selection fields of this
// atom_selection_container_t.
// 
// atom_selection_container_t 
// make_atom_selection(int imol, const coot::minimol::molecule &mol) {

// }

//       std::vector<coot::minimol::atom *> atoms = range_mol.select_atoms_serial();
//       for (int i=0; i<atoms.size(); i++) {
// 	 std::cout << "range mol atom: " << atoms[i]->pos.format() << std::endl;
//       }
