/* src/c-interface.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 The University of York
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#include <stdlib.h>
#include <iostream>

#ifdef USE_GUILE
#include <guile/gh.h>
#endif // USE_GUILE

#ifdef USE_PYTHON
#include "Python.h"
#endif // USE_PYTHON

#define HAVE_CIF  // will become unnessary at some stage.

#include <sys/types.h> // for stating
#include <sys/stat.h>
#if !defined _MSC_VER
#include <unistd.h>
#else
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

#include "Cartesian.h"
#include "Bond_lines.h"

// #include "xmap-interface.h"
#include "graphics-info.h"

// #include "atom-utils.h" // asc_to_graphics
// #include "db-main.h" not yet

#include "coot-coord-utils.hh"
#include "coot-fasta.hh"

#include "BuildCas.h"
#include "helix-placement.hh"

#include "trackball.h" // adding exportable rotate interface

#include "coot-utils.hh"  // for is_member_p
#include "coot-map-heavy.hh"  // for fffear
#include "c-interface.h"
#include "cc-interface.hh"

#include "ideal-rna.hh"
#include "cmtz-interface.hh" // for valid columns mtz_column_types_info_t

/*  ------------------------------------------------------------------------ */
/*                   Maps - (somewhere else?):                               */
/*  ------------------------------------------------------------------------ */
/*! \brief Calculate SFs from an MTZ file and generate a map. 
 @return the new molecule number. */
int map_from_mtz_by_calc_phases(const char *mtz_file_name, 
				const char *f_col, 
				const char *sigf_col,
				int imol_coords) {

   int ir = -1; // return value
   graphics_info_t g;
   int imol_map = g.n_molecules;
   if (is_valid_model_molecule(imol_coords)) { 
      std::string m(mtz_file_name);
      std::string f(f_col);
      std::string s(sigf_col);
      atom_selection_container_t a = g.molecules[imol_coords].atom_sel;
      short int t = molecule_map_type::TYPE_2FO_FC;
      int istat = g.molecules[imol_map].make_map_from_mtz_by_calc_phases(m,f,s,a,t);
      if (istat != -1) {
	 g.n_molecules++;
	 graphics_draw();
	 ir = imol_map;
      } else {
	 ir = -1; // error
      }
   } 
   return ir;
} 


/*! \brief fire up a GUI, which asks us which model molecule we want
  to calc phases from.  On "OK" button there, we call
  map_from_mtz_by_refmac_calc_phases() */
void calc_phases_generic(const char *mtz_file_name) {

   graphics_info_t g;
   coot::mtz_column_types_info_t r = coot::get_f_phi_columns(mtz_file_name);
   if (r.f_cols.size() == 0) {
      std::cout << "No Fobs found in " << mtz_file_name << std::endl;
      std::string s =  "No Fobs found in ";
      s += mtz_file_name;
      g.statusbar_text(s);
   } else { 
      if (r.sigf_cols.size() == 0) {
	 std::cout << "No SigFobs found in " << mtz_file_name << std::endl;
	 std::string s =  "No SigFobs found in ";
	 s += mtz_file_name;
	 g.statusbar_text(s);
      } else {
	 // normal path:
	 std::string f_obs_col = r.f_cols[0].column_label;
	 std::string sigfobs_col = r.sigf_cols[0].column_label;
	 std::vector<std::string> v;
	 v.push_back("refmac-for-phases-and-make-map");
	 v.push_back(coot::util::single_quote(mtz_file_name));
	 v.push_back(coot::util::single_quote(f_obs_col));
	 v.push_back(coot::util::single_quote(sigfobs_col));
	 std::string c = languagize_command(v);
	 std::cout << "command: " << c << std::endl;
#ifdef USE_GUILE
	 safe_scheme_command(c);
#else
#ifdef USE_PYTHON
	 safe_python_command(c);
#endif
#endif
      }
   }
}

/*! \brief Calculate SFs (using refmac optionally) from an MTZ file
  and generate a map. Get F and SIGF automatically (first of their
  type) from the mtz file.

@return the new molecule number, -1 on a problem. */
int map_from_mtz_by_refmac_calc_phases(const char *mtz_file_name, 
				       const char *f_col, 
				       const char *sigf_col, 
				       int imol_coords) {

   int istat = -1;


   return istat;
} 


/*  ------------------------------------------------------------------------ */
/*                   model/fit/refine functions:                             */
/*  ------------------------------------------------------------------------ */
void set_model_fit_refine_rotate_translate_zone_label(const char *txt) {
   graphics_info_t::model_fit_refine_rotate_translate_zone_string = txt;
}

void set_model_fit_refine_place_atom_at_pointer_label(const char *txt) {
   graphics_info_t::model_fit_refine_place_atom_at_pointer_string = txt;
}

int copy_molecule(int imol) {
   int iret = -1;
   if (is_valid_model_molecule(imol)) {
      int new_mol_number = graphics_info_t::n_molecules;
      CMMDBManager *m = graphics_info_t::molecules[imol].atom_sel.mol;
      CMMDBManager *n = new CMMDBManager;
      n->Copy(m, MMDBFCM_All);
      atom_selection_container_t asc = make_asc(n);
      std::string label = "Copy of ";
      label += graphics_info_t::molecules[imol].name_;
      graphics_info_t g;
      g.expand_molecule_space_maybe(); 
      graphics_info_t::molecules[new_mol_number].install_model(asc, label, 1);
      iret = graphics_info_t::n_molecules;
      graphics_info_t::n_molecules++;
   }
   if (is_valid_map_molecule(imol)) {
      int new_mol_number = graphics_info_t::n_molecules;
      std::string label = "Copy of ";
      label += graphics_info_t::molecules[imol].name_;
      graphics_info_t g;
      g.expand_molecule_space_maybe(); 
      graphics_info_t::molecules[new_mol_number].new_map(graphics_info_t::molecules[imol].xmap_list[0], label);
      if (graphics_info_t::molecules[imol].is_difference_map_p()) {
	 graphics_info_t::molecules[new_mol_number].set_map_is_difference_map();
      }
      iret = graphics_info_t::n_molecules;
      graphics_info_t::n_molecules++;
   }
   if (iret != -1) 
      graphics_draw();
   return iret;
}

/*! \brief replace the parts of molecule number imol that are
  duplicated in molecule number imol_frag */
int replace_fragment(int imol_target, int imol_fragment,
		     const char *mmdb_atom_selection_str) {

   int istate = 0;
   if (is_valid_model_molecule(imol_target)) {
      if (is_valid_model_molecule(imol_fragment)) {
	 PPCAtom atom_selection = NULL;
	 int n_atoms;
	 CMMDBManager *mol = graphics_info_t::molecules[imol_fragment].atom_sel.mol;
	 int SelHnd = mol->NewSelection();
	 mol->Select(SelHnd, STYPE_ATOM, (char *) mmdb_atom_selection_str, SKEY_OR);
	 CMMDBManager *mol_new =
	    coot::util::create_mmdbmanager_from_atom_selection(mol, SelHnd);
	 atom_selection_container_t asc = make_asc(mol_new);
	 istate = graphics_info_t::molecules[imol_target].replace_fragment(asc);
	 mol->DeleteSelection(SelHnd);
	 graphics_draw();
      }
   }
   return istate;
} 

/*  ------------------------------------------------------------------------ */
/*                         backup/undo functions:                            */
/*  ------------------------------------------------------------------------ */

void turn_off_backup(int imol) {
   if (imol < graphics_n_molecules())
      graphics_info_t::molecules[imol].turn_off_backup();
} 

void turn_on_backup(int imol) {
   if (imol < graphics_n_molecules())
      graphics_info_t::molecules[imol].turn_on_backup();
} 

void apply_undo() {		/* "Undo" button callback */
   graphics_info_t g;
   g.apply_undo();
}

void apply_redo() { 
   graphics_info_t g;
   g.apply_redo();
}

void set_undo_molecule(int imol) {

   // Potentially a problem here.  What if we set the undo molecule to
   // a molecule that has accidentally just deleted the only ligand
   // residue.
   //
   // Too bad.
   //
//    if (is_valid_model_molecule(imol)) {
//       graphics_info_t g;
//       g.set_undo_molecule_number(imol);
//    }
   
   // 20060522 so how about I check that the index is within limits?
   //          and then ask if if the mol is valid (rather than the
   //          number of atoms selected):

   if ((imol >= 0) && (imol < graphics_info_t::n_molecules)) {
      graphics_info_t g;
      if (g.molecules[imol].atom_sel.mol) {
	 std::cout << "INFO:: undo molecule number set to: " << imol << std::endl;
	 g.set_undo_molecule_number(imol);
      }
   }
}

/*! \brief show the Undo Molecule chooser - i.e. choose the molecule
  to which the "Undo" button applies. */
void show_set_undo_molecule_chooser() {

#ifdef USE_GUILE

   std::string s("(molecule-chooser-gui ");
   s += "\"Choose Molecule for Undo operations\" ";
   s += "set-undo-molecule)";
   // std::cout << s << std::endl;

   safe_scheme_command(s);

#endif // USE_GUILE   

} 



void set_unpathed_backup_file_names(int state) {
   graphics_info_t::unpathed_backup_file_names_flag = state;
}

int  unpathed_backup_file_names_state() {
   return graphics_info_t::unpathed_backup_file_names_flag;
}



/*  ------------------------------------------------------------------------ */
/*                    terminal residue functions:                            */
/*  ------------------------------------------------------------------------ */

void do_add_terminal_residue(short int state) { 

   graphics_info_t g;
   g.in_terminal_residue_define = state;
   if (state) {
      int imol_map = g.Imol_Refinement_Map();
      if (imol_map >= 0) { 
	 std::cout << "click on an atom of a terminal residue" << std::endl;
	 g.pick_cursor_maybe();
	 g.pick_pending_flag = 1;
      } else {
	 g.show_select_map_dialog();
	 g.in_terminal_residue_define = 0;
	 g.model_fit_refine_unactive_togglebutton("model_refine_dialog_fit_terminal_residue_togglebutton");
	 g.normal_cursor();
      }
   } else {
      g.normal_cursor();
   } 
}

void 
set_add_terminal_residue_n_phi_psi_trials(int n) { 
   graphics_info_t g;
   g.add_terminal_residue_n_phi_psi_trials = n;
}

void
set_add_terminal_residue_add_other_residue_flag(int i) {
   graphics_info_t::add_terminal_residue_add_other_residue_flag = i;
}

void set_terminal_residue_do_rigid_body_refine(short int v) { 

   graphics_info_t g;
   g.terminal_residue_do_rigid_body_refine = v;

}

void set_add_terminal_residue_do_post_refine(short int istat) {
   graphics_info_t::add_terminal_residue_do_post_refine = istat;
}



int add_terminal_residue_immediate_addition_state() {
   int i = graphics_info_t::add_terminal_residue_immediate_addition_flag;
   return i;
} 

void set_add_terminal_residue_immediate_addition(int i) {
   graphics_info_t::add_terminal_residue_immediate_addition_flag = i;
}

// return 0 on failure, 1 on success
//
int add_terminal_residue(int imol, 
			 char *chain_id, 
			 int residue_number,
			 char *residue_type,
			 int immediate_add) {

   int istate = 0;
   if (imol >= 0 ) {
      if (imol <= graphics_info_t::n_molecules) {
	 if (graphics_info_t::molecules[imol].has_model()) {

	    // We don't do this as a member function of
	    // molecule_class_info_t because we are using
	    // graphics_info_t::execute_add_terminal_residue, which
	    // does the molecule manipulations outside of the
	    // molecule_class_info_t class.
	    //

	    graphics_info_t g;
	    int atom_indx = atom_index(imol, chain_id, residue_number, " CA ");
	    if (atom_indx >= 0) {
	       std::string term_type = g.molecules[imol].get_term_type(atom_indx);
	       std::string inscode = "";
	       CResidue *res_p =
		  g.molecules[imol].residue_from_external(residue_number, inscode,
							  std::string(chain_id));
	       g.execute_add_terminal_residue(imol, term_type, res_p, chain_id,
					      residue_type, immediate_add);
	       istate = 1;
	    } else {
	       std::cout << "WARNING:: in add_terminal_residue: "
			 << " Can't find atom index for CA in residue "
			 << residue_number << " " << chain_id << std::endl;
	    }
	 }
      }
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("add-terminal-residue");
   command_strings.push_back(graphics_info_t::int_to_string(imol));
   command_strings.push_back(single_quote(chain_id));
   command_strings.push_back(graphics_info_t::int_to_string(residue_number));
   command_strings.push_back(graphics_info_t::int_to_string(immediate_add));
   add_to_history(command_strings);
   return istate;
}

void set_add_terminal_residue_default_residue_type(const char *type) {

   if (type) 
      graphics_info_t::add_terminal_residue_type = type;
}



/*  ----------------------------------------------------------------------- */
/*                  rotate/translate buttons                                */
/*  ----------------------------------------------------------------------- */

void do_rot_trans_setup(short int state) { 
   graphics_info_t g;
   g.in_rot_trans_object_define = state;
   if (state){ 
      g.pick_cursor_maybe();
      std::cout << "click on 2 atoms to define a zone" << std::endl;
      g.pick_pending_flag = 1;
   } else {
      g.pick_pending_flag = 0;
      g.normal_cursor();
   }
}

void do_rot_trans_adjustments(GtkWidget *dialog) { 
   graphics_info_t g;
   g.do_rot_trans_adjustments(dialog);
}


void rot_trans_reset_previous() { 
   graphics_info_t g;
   // rot_trans adjustments:
   for (int i=0; i<6; i++) 
      g.previous_rot_trans_adjustment[i] = -10000;
}

void set_rotate_translate_zone_rotates_about_zone_centre(int istate) {
   graphics_info_t::rot_trans_zone_rotates_about_zone_centre = istate;
} 



/*  ----------------------------------------------------------------------- */
/*                  delete residue                                          */
/*  ----------------------------------------------------------------------- */
void delete_residue(int imol, const char *chain_id, int resno, const char *inscode) {

   graphics_info_t g;
   short int istat = g.molecules[imol].delete_residue(chain_id, resno, std::string(inscode));
   if (istat) { 
      // now if the go to atom widget was being displayed, we need to
      // redraw the residue list and atom list (if the molecule of the
      // residue and atom list is the molecule that has just been
      // deleted)

      g.update_go_to_atom_window_on_changed_mol(imol);

      graphics_draw();
   } else { 
      std::cout << "failed to delete residue " << chain_id 
		<< " " << resno << "\n";
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("delete-residue");
   command_strings.push_back(g.int_to_string(imol));
   command_strings.push_back(single_quote(chain_id));
   command_strings.push_back(g.int_to_string(resno));
   command_strings.push_back(single_quote(std::string(inscode)));
   add_to_history(command_strings);
}

void delete_residue_hydrogens(int imol,
			      const char *chain_id,
			      int resno,
			      const char *inscode,
			      const char *altloc) {

   graphics_info_t g;
   if (is_valid_model_molecule(imol)) { 
      short int istat = g.molecules[imol].delete_residue_hydrogens(chain_id, resno, inscode, altloc);
      if (istat) { 
	 // now if the go to atom widget was being displayed, we need to
	 // redraw the residue list and atom list (if the molecule of the
	 // residue and atom list is the molecule that has just been
	 // deleted)

	 g.update_go_to_atom_window_on_changed_mol(imol);
	 graphics_draw();

      } else { 
	 std::cout << "failed to delete residue hydrogens " << chain_id 
		   << " " << resno << "\n";
      }
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("delete-residue-hydrogens");
   command_strings.push_back(g.int_to_string(imol));
   command_strings.push_back(single_quote(chain_id));
   command_strings.push_back(g.int_to_string(resno));
   command_strings.push_back(single_quote(inscode));
   command_strings.push_back(single_quote(altloc));
   add_to_history(command_strings);
} 


void delete_residue_with_altconf(int imol,
				 const char *chain_id,
				 int resno,
				 const char *inscode,
				 const char *altloc) {
   std::string altconf(altloc);
   graphics_info_t g;
   short int istat =
      g.molecules[imol].delete_residue_with_altconf(chain_id, resno, inscode, altconf);
   
   if (istat) { 
      // now if the go to atom widget was being displayed, we need to
      // redraw the residue list and atom list (if the molecule of the
      // residue and atom list is the molecule that has just been
      // deleted)
      // 
      g.update_go_to_atom_window_on_changed_mol(imol);

      graphics_draw();
   } else { 
      std::cout << "failed to delete residue atoms " << chain_id 
		<< " " << resno << " :" << altconf << ":\n";
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("delete-residue-with-altconf");
   command_strings.push_back(g.int_to_string(imol));
   command_strings.push_back(single_quote(chain_id));
   command_strings.push_back(g.int_to_string(resno));
   command_strings.push_back(single_quote(inscode));
   command_strings.push_back(single_quote(altloc));
   add_to_history(command_strings);
}

void delete_residue_sidechain(int imol, const char *chain_id, int resno, const char *ins_code) {

   std::string inscode(ins_code);
   graphics_info_t g;

   if (is_valid_model_molecule(imol)) { 
      short int istat =
	 g.molecules[imol].delete_residue_sidechain(std::string(chain_id), resno,
						    inscode);
      
      if (istat) {
      g.update_go_to_atom_window_on_changed_mol(imol);
      graphics_draw();
      }
      
      if (graphics_info_t::delete_item_widget != NULL) {
	 GtkWidget *checkbutton = lookup_widget(graphics_info_t::delete_item_widget,
						"delete_item_keep_active_checkbutton");
	 if (GTK_TOGGLE_BUTTON(checkbutton)->active) {
	    // dont destroy it
	 } else {
	    gint upositionx, upositiony;
	    gdk_window_get_root_origin (graphics_info_t::delete_item_widget->window,
					&upositionx, &upositiony);
	    graphics_info_t::delete_item_widget_x_position = upositionx;
	    graphics_info_t::delete_item_widget_y_position = upositiony;
	    gtk_widget_destroy(graphics_info_t::delete_item_widget);
	    graphics_info_t::delete_item_widget = NULL;
	 }
      }
   }

   std::string cmd = "delete-residue-sidechain";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(chain_id);
   args.push_back(resno);
   args.push_back(ins_code);
   add_to_history_typed(cmd, args);
}


void set_add_alt_conf_new_atoms_occupancy(float f) {  /* default 0.5 */

   graphics_info_t g;
   g.add_alt_conf_new_atoms_occupancy = f;
}


int set_atom_attribute(int imol, const char *chain_id, int resno, const char *ins_code, const char *atom_name, const char*alt_conf, const char *attribute_name, float val) {
   int istat = 0;
   if (is_valid_model_molecule(imol)) {
      istat = graphics_info_t::molecules[imol].set_atom_attribute(chain_id, resno, ins_code, atom_name, alt_conf, attribute_name, val);
   }
   graphics_draw();
   return istat;
} 

int set_atom_string_attribute(int imol, const char *chain_id, int resno, const char *ins_code, const char *atom_name, const char*alt_conf, const char *attribute_name, const char *attribute_value) {
   int istat = 0; 
   if (is_valid_model_molecule(imol)) {
      istat = graphics_info_t::molecules[imol].set_atom_string_attribute(chain_id, resno, ins_code, atom_name, alt_conf, attribute_name, attribute_value);
      graphics_draw();
   }
   return istat;
}




// imol has changed.
// Now fix up the Go_To_Atom window to match:
// 
void update_go_to_atom_window_on_changed_mol(int imol) {

   // now if the go to atom widget was being displayed, we need to
   // redraw the residue list and atom list (if the molecule of the
   // residue and atom list is the molecule that has just been
   // deleted)
   graphics_info_t g;
   g.update_go_to_atom_window_on_changed_mol(imol);
}

// a new molecule has has been read in.
// 
// Now fix up the Go_To_Atom window to match by changing the option menu
// 
void update_go_to_atom_window_on_new_mol() {

   graphics_info_t g;
   g.update_go_to_atom_window_on_new_mol();
   add_to_history_simple("update-go-to-atom-window-on-new-mol");
}


void update_go_to_atom_window_on_other_molecule_chosen(int imol) {
   graphics_info_t g;
   g.update_go_to_atom_window_on_other_molecule_chosen(imol);
   add_to_history_simple("update-go-to-atom-window-on-other-molecule-chosen");

} 

void delete_atom(int imol, const char *chain_id, int resno, const char *at_name, 
		 const char *altLoc) {

   graphics_info_t g;
   short int istat = g.molecules[imol].delete_atom(chain_id, resno, at_name, altLoc);
   if (istat) { 
      // now if the go to atom widget was being displayed, we need to
      // redraw the residue list and atom list (if the molecule of the
      // residue and atom list is the molecule that has just been
      // deleted)
      //

      g.update_go_to_atom_window_on_changed_mol(imol);
      if (g.go_to_atom_window) {
	 int go_to_atom_imol = g.go_to_atom_molecule();
	 if (go_to_atom_imol == imol) { 

	    // The go to atom molecule matched this molecule, so we
	    // need to regenerate the residue and atom lists.
	    GtkWidget *gtktree = lookup_widget(g.go_to_atom_window,
					       "go_to_atom_residue_tree");
#if (GTK_MAJOR_VERSION == 1) || defined (GTK_ENABLE_BROKEN)
	    g.fill_go_to_atom_residue_list_gtk1(gtktree);
#else 	    
	    g.fill_go_to_atom_residue_tree_gtk2(gtktree);
#endif	    
	 } else { 
	    std::cout << "DEBUG:: no molecule number match\n";
	 } 
      }
      graphics_draw();
   } else { 
      std::cout << "failed to delete atom " << chain_id 
		<< " " << resno << " " <<  at_name << "\n";
   }

   std::string cmd = "delete-atom";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(chain_id);
   args.push_back(resno);
   args.push_back(at_name);
   args.push_back(altLoc);
   add_to_history_typed(cmd, args);

} 

void set_delete_atom_mode() {
   graphics_info_t g;
   g.delete_item_atom = 1;
   g.delete_item_residue_zone = 0;
   g.delete_item_residue_hydrogens = 0;
   g.delete_item_residue = 0;
   g.delete_item_sidechain = 0; 
}

void set_delete_residue_mode() {
   graphics_info_t g;
   g.delete_item_atom = 0;
   g.delete_item_residue_zone = 0;
   g.delete_item_residue_hydrogens = 0;
   g.delete_item_water = 0;
   g.delete_item_residue = 1;
   g.delete_item_sidechain = 0; 
}

void set_delete_residue_hydrogens_mode() {

   graphics_info_t g;
   g.delete_item_residue = 0;
   g.delete_item_residue_zone = 0;
   g.delete_item_atom = 0;
   g.delete_item_water = 0;
   g.delete_item_residue_hydrogens = 1;
   g.delete_item_sidechain = 0; 

}

void set_delete_residue_zone_mode() {

   graphics_info_t g;
   g.delete_item_residue = 0;
   g.delete_item_residue_zone = 1;
   g.delete_item_atom = 0;
   g.delete_item_water = 0;
   g.delete_item_residue_hydrogens = 0;
   g.delete_item_sidechain = 0; 

} 

void set_delete_water_mode() {

   graphics_info_t g;
   g.delete_item_residue = 0;
   g.delete_item_residue_zone = 0;
   g.delete_item_water = 1;
   g.delete_item_atom = 0;
   g.delete_item_residue_hydrogens = 0;
   g.delete_item_sidechain = 0; 

} 

void set_delete_sidechain_mode() {

   graphics_info_t g;
   g.delete_item_residue = 0;
   g.delete_item_residue_zone = 0;
   g.delete_item_water = 0;
   g.delete_item_atom = 0;
   g.delete_item_residue_hydrogens = 0;
   g.delete_item_sidechain = 1; 

}

// (predicate) a boolean (or else it's residue)
//
// Used by on_model_refine_dialog_delete_button_clicked callback to
// determine if the Atom checkbutton should be active when the new
// dialog is displayed.
// 
short int delete_item_mode_is_atom_p() {
   short int v=0;
   if (graphics_info_t::delete_item_atom == 1)
      v = 1;
   if (graphics_info_t::delete_item_residue == 1)
      v = 0;
   if (graphics_info_t::delete_item_water == 1)
      v = 0;
   return v;
}

short int delete_item_mode_is_residue_p() {
   short int v=0;
   if (graphics_info_t::delete_item_residue == 1)
      v = 1;
   return v;
}

short int delete_item_mode_is_water_p() {
   short int v=0;
   if (graphics_info_t::delete_item_water == 1)
      v = 1;
   return v;
}

short int delete_item_mode_is_sidechain_p() {
   short int v=0;
   if (graphics_info_t::delete_item_sidechain == 1)
      v = 1;
   return v;
}

void delete_atom_by_atom_index(int imol, int index, short int do_delete_dialog) {
   graphics_info_t g;

   std::string atom_name = g.molecules[imol].atom_sel.atom_selection[index]->name;
   std::string chain_id  = g.molecules[imol].atom_sel.atom_selection[index]->GetChainID();
   int resno             = g.molecules[imol].atom_sel.atom_selection[index]->GetSeqNum();
   std::string altconf   = g.molecules[imol].atom_sel.atom_selection[index]->altLoc;

   // g.molecules[imol].delete_atom(chain_id, resno, atom_name);
   delete_atom(imol, chain_id.c_str(), resno, atom_name.c_str(), altconf.c_str());
   if (graphics_info_t::delete_item_widget != NULL) {
      if (do_delete_dialog) { // via ctrl

	 // another check is needed, is the check button active?
	 // 
	 // If not we can go ahead and delete the dialog
	 //
	 GtkWidget *checkbutton = lookup_widget(graphics_info_t::delete_item_widget,
						"delete_item_keep_active_checkbutton");
	 if (GTK_TOGGLE_BUTTON(checkbutton)->active) {
	    // don't kill the widget
	    pick_cursor_maybe(); // it was set to normal_cursor() in
                                 // graphics-info-define's delete_item().
	 } else {
	 
	    gint upositionx, upositiony;
	    gdk_window_get_root_origin (graphics_info_t::delete_item_widget->window,
					&upositionx, &upositiony);
	    graphics_info_t::delete_item_widget_x_position = upositionx;
	    graphics_info_t::delete_item_widget_y_position = upositiony;
	    gtk_widget_destroy(graphics_info_t::delete_item_widget);
	    graphics_info_t::delete_item_widget = NULL;
	    graphics_draw();
	 }
      }
   }
   std::string cmd = "delete-atom-by-atom-index";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(index);
   args.push_back(do_delete_dialog);
   add_to_history_typed(cmd, args);

}

void delete_residue_by_atom_index(int imol, int index, short int do_delete_dialog_by_ctrl) {

   graphics_info_t g;
   std::string chain_id  = g.molecules[imol].atom_sel.atom_selection[index]->GetChainID();
   int resno             = g.molecules[imol].atom_sel.atom_selection[index]->GetSeqNum();
   std::string altloc    = g.molecules[imol].atom_sel.atom_selection[index]->altLoc;
   std::string inscode   = g.molecules[imol].atom_sel.atom_selection[index]->GetInsCode();

   if (altloc == "") 
      delete_residue(imol, chain_id.c_str(), resno, inscode.c_str());
   else
      delete_residue_with_altconf(imol, chain_id.c_str(), resno, inscode.c_str(), altloc.c_str());
   
   if (graphics_info_t::delete_item_widget != NULL) {
      short int do_delete_dialog = do_delete_dialog_by_ctrl;
      GtkWidget *checkbutton = lookup_widget(graphics_info_t::delete_item_widget,
					     "delete_item_keep_active_checkbutton");
      if (GTK_TOGGLE_BUTTON(checkbutton)->active) { 
	 do_delete_dialog = 0;
	 pick_cursor_maybe(); // it was set to normal_cursor() in
			      // graphics-info-define's delete_item().
      }
      if (do_delete_dialog) { 
	 gint upositionx, upositiony;
	 gdk_window_get_root_origin (graphics_info_t::delete_item_widget->window,
				     &upositionx, &upositiony);
	 graphics_info_t::delete_item_widget_x_position = upositionx;
	 graphics_info_t::delete_item_widget_y_position = upositiony;
	 gtk_widget_destroy(graphics_info_t::delete_item_widget);
	 graphics_info_t::delete_item_widget = NULL;
      }
   }
   graphics_draw();
}

void delete_residue_hydrogens_by_atom_index(int imol, int index, short int do_delete_dialog) {

   graphics_info_t g;
   std::string chain_id  = g.molecules[imol].atom_sel.atom_selection[index]->GetChainID();
   int resno             = g.molecules[imol].atom_sel.atom_selection[index]->GetSeqNum();
   std::string altloc    = g.molecules[imol].atom_sel.atom_selection[index]->altLoc;
   std::string inscode   = g.molecules[imol].atom_sel.atom_selection[index]->GetInsCode();


   delete_residue_hydrogens(imol, chain_id.c_str(), resno, inscode.c_str(), altloc.c_str());
   
   if (graphics_info_t::delete_item_widget != NULL) {
      if (do_delete_dialog) { 
	 GtkWidget *checkbutton = lookup_widget(graphics_info_t::delete_item_widget,
						"delete_item_keep_active_checkbutton");
	 if (GTK_TOGGLE_BUTTON(checkbutton)->active) {
	    // dont destroy it
	 } else {
	    // save the position of the window and kill it off.
	    gint upositionx, upositiony;
	    gdk_window_get_root_origin (graphics_info_t::delete_item_widget->window,
					&upositionx, &upositiony);
	    graphics_info_t::delete_item_widget_x_position = upositionx;
	    graphics_info_t::delete_item_widget_y_position = upositiony;
	    gtk_widget_destroy(graphics_info_t::delete_item_widget);
	    graphics_info_t::delete_item_widget = NULL;
	 }
      }
   }
   graphics_draw();
}

// Deletes all altconfs, the whole residue goes.
// 
void delete_residue_range(int imol, const char *chain_id, int resno_start, int resno_end) {

   // altconf is ignored currently.
   // 
   coot::residue_spec_t res1(chain_id, resno_start);
   coot::residue_spec_t res2(chain_id, resno_end);

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      g.delete_residue_range(imol, res1, res2);
      if (graphics_info_t::go_to_atom_window) {
	 update_go_to_atom_window_on_changed_mol(imol);
      }
   }
   graphics_draw();
}



void store_delete_item_widget(GtkWidget *widget) {
   graphics_info_t::delete_item_widget = widget;
}

void clear_pending_delete_item() { 

   graphics_info_t g;
   g.delete_item_atom = 0;
   g.delete_item_residue = 0;
   g.delete_item_residue_zone = 0;
   g.delete_item_residue_hydrogens = 0;

}

// Do we want this function really.  Consider binning it.
//
void set_keep_delete_item_active_state(int istate) {

   // graphics_info_t g;
} 


/* We need to set the pending delete flag and that can't be done in
   callback, so this wrapper does it */
GtkWidget *wrapped_create_delete_item_dialog() {

   GtkWidget *widget = create_delete_item_dialog();
   GtkWidget *atom_toggle_button;

   if (delete_item_mode_is_atom_p()) { 
      atom_toggle_button = lookup_widget(GTK_WIDGET(widget),
					 "delete_item_atom_radiobutton");
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(atom_toggle_button), TRUE);
      std::cout << "Click on the atom that you wish to delete\n";
   } else {
      if (delete_item_mode_is_water_p()) {
	 GtkWidget *water_toggle_button = lookup_widget(widget,
							"delete_item_water_radiobutton");
	 gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(water_toggle_button), TRUE);
      } else { 
	 set_delete_residue_mode(); // The widget default radio button
	 std::cout << "Click on an atom in the residue that you wish to delete\n";
      }
   }
   graphics_info_t::pick_pending_flag = 1;
   pick_cursor_maybe();
   set_transient_and_position(COOT_DELETE_WINDOW, widget);
   store_delete_item_widget(widget);
   return widget; 
}

// -----------------------------------------------------
//  move molecule here widget
// -----------------------------------------------------
GtkWidget *wrapped_create_move_molecule_here_dialog() {

   GtkWidget *w = create_move_molecule_here_dialog();
   GtkWidget *option_menu = lookup_widget(w, "move_molecule_here_optionmenu"); 
   int imol = first_coords_imol();
   graphics_info_t::move_molecule_here_molecule_number = imol;
   GtkSignalFunc callback_func = GTK_SIGNAL_FUNC(graphics_info_t::move_molecule_here_item_select);
   
   graphics_info_t g;
   g.fill_option_menu_with_coordinates_options(option_menu, callback_func, imol);

   return w;
}


void move_molecule_here_by_widget(GtkWidget *w) {

   int imol = graphics_info_t::move_molecule_here_molecule_number;
   // std::cout << "move_molecule_here imol: " << imol << std::endl;
   if (is_valid_model_molecule(imol)) {
      
      // (move-molecule-here imol)
      coot::Cartesian cen =
	 centre_of_molecule(graphics_info_t::molecules[imol].atom_sel);

      graphics_info_t g;
      coot::Cartesian rc = g.RotationCentre();

      float x = rc.x() - cen.x();
      float y = rc.y() - cen.y();
      float z = rc.z() - cen.z();

      translate_molecule_by(imol, x, y, z);
      
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("move-molecule-here");
   command_strings.push_back(clipper::String(imol));
   add_to_history(command_strings);
}


// ---------------------------------------------------------------------
//                 rotamer
// ---------------------------------------------------------------------
// 

void set_write_peaksearched_waters() {
   graphics_info_t g;
   g.ligand_water_write_peaksearched_atoms = 1;
} 


// Called by the Model/Fit/Refine Rotamers button callback.
void setup_rotamers(short int state) {
   graphics_info_t g;
   g.in_rotamer_define = state;
   if (state) { 
      g.pick_cursor_maybe();
      g.pick_pending_flag = 1;
      std::cout << "Click on an atom in a residue for which you wish to see rotamers"
		<< std::endl;
   } else {
      g.normal_cursor();
   }
} 


void do_rotamers(int atom_index, int imol) {

   std::cout << "     Rotamer library:" << std::endl;
   std::cout << "     R. L. Dunbrack, Jr. and F. E. Cohen." << std::endl;
   std::cout << "     Bayesian statistical analysis of ";
   std::cout << "protein sidechain rotamer preferences" << std::endl;
   std::cout << "     Protein Science, 6, 1661-1681 (1997)." << std::endl;
   std::cout << "" << std::endl;

   graphics_info_t g;
   g.do_rotamers(atom_index, imol); 

}

void
set_rotamer_lowest_probability(float f) {
   graphics_info_t g;
   g.rotamer_lowest_probability = f;
}

void
set_rotamer_check_clashes(int i) {
   graphics_info_t::rotamer_fit_clash_flag = i;
}

// Return -999 on imol indexing error
// -99.9 in class function error
// 
float
auto_fit_best_rotamer(int resno,
		      const char *altloc,
		      const char *insertion_code, 
		      const char *chain_id, int imol_coords, int imol_map,
		      int clash_flag, float lowest_probability) {

   float f = -999.9;

   if (imol_coords < graphics_n_molecules()) {
      if (imol_map < graphics_n_molecules()) {
	 // if (graphics_info_t::molecules[imol_map].has_map()) {
	 std::string ins(insertion_code);
	 std::string chain(chain_id);
	 if (imol_map < 0 ) {
	    std::cout << "INFO:: fitting rotamers by clash score only " << std::endl;
	    f = graphics_info_t::molecules[imol_coords].auto_fit_best_rotamer(resno, altloc, ins,
									      chain, imol_map,
									      1,
									      lowest_probability);
	 } else {
	    if (graphics_info_t::molecules[imol_map].has_map()) {
	       f = graphics_info_t::molecules[imol_coords].auto_fit_best_rotamer(resno, altloc, ins,
										 chain, imol_map,
										 clash_flag,
										 lowest_probability);

	       // get the residue so that it can update the geometry graph
	       CResidue *residue_p =
		  graphics_info_t::molecules[imol_coords].get_residue(resno, ins, chain);
	       if (residue_p) {
		  graphics_info_t g;
		  g.update_geometry_graphs(&residue_p, 1, imol_coords, imol_map);
	       }
	       std::cout << "Fitting score for best rotamer: " << f << std::endl;
	    }
	 }
	 graphics_draw();
      }
   }
   return f;
}


void
set_auto_fit_best_rotamer_clash_flag(int i) { /* 0 off, 1 on */
   graphics_info_t::rotamer_fit_clash_flag = i;
} 

void
setup_auto_fit_rotamer(short int state) {
   graphics_info_t::in_auto_fit_define = state;
   if (state) { 
      graphics_info_t::pick_cursor_maybe();
      graphics_info_t::pick_pending_flag = 1;
      std::cout << "Click on an atom in the residue that you wish to fit\n";
   } else {
      graphics_info_t::normal_cursor();
   }
}


// FIXME  (autofit rotamer seems to return a score OK).
float
rotamer_score(int resno, const char *insertion_code,
	      const char *chain_id, int imol_coords, int imol_map,
	      int rotamer_number) {
   float f = 0;
   return f;
}

void
set_graphics_rotamer_dialog(GtkWidget *w) {
   std::cout << "DEBUG:: setting graphics rotamer_dialog to " << w << std::endl;
   graphics_info_t::rotamer_dialog = w;
}




// ---------------------------------------------------------------------
//                 mutation 
// ---------------------------------------------------------------------
// 
void
setup_mutate(short int state) {

   graphics_info_t g;
   g.in_mutate_define = state;
   if (state) { 
      g.pick_cursor_maybe();
      g.pick_pending_flag = 1;
      std::cout << "Click on an atom in a residue which you wish to mutate"
		<< std::endl;
   } else {
      g.normal_cursor();
   } 
}

void
setup_mutate_auto_fit(short int state) { 

   graphics_info_t g;

   if (state) { 
      int imol_map = g.Imol_Refinement_Map(); 
      if (imol_map >= 0) { 
	 std::cout << "Click on an atom in a residue which you wish to mutate"
		   << std::endl;
	 g.in_mutate_auto_fit_define = state;
	 g.pick_cursor_maybe();
	 g.pick_pending_flag = 1;
      } else { 
	 // map chooser dialog
	 g.show_select_map_dialog();
	 g.in_mutate_auto_fit_define = 0;
	 normal_cursor();
	 g.model_fit_refine_unactive_togglebutton("model_refine_dialog_mutate_auto_fit_togglebutton");
      }
   } else {
      g.in_mutate_auto_fit_define = state;
      g.normal_cursor();
   }
}


/* 1 for yes, 0 for no. */
void set_mutate_auto_fit_do_post_refine(short int istate) {

   graphics_info_t::mutate_auto_fit_do_post_refine_flag = istate;
} 

/*! \brief what is the value of the previous flag? */
int mutate_auto_fit_do_post_refine_state() {
   return graphics_info_t::mutate_auto_fit_do_post_refine_flag;
} 



/*! \brief set a flag saying that the chosen residue should only be
  added as a stub (mainchain + CB) */
void set_residue_type_chooser_stub_state(short int istat) {

   graphics_info_t::residue_type_chooser_stub_flag = istat;
}


void
do_mutation(const char *type, short int stub_flag) {
   graphics_info_t g;
   // use g.mutate_residue_atom_index and g.mutate_residue_imol
   g.do_mutation(type, stub_flag);

}

void
place_atom_at_pointer() {
   graphics_info_t g;
   if (g.pointer_atom_is_dummy)
      g.place_dummy_atom_at_pointer();
   else {
      // put up a widget which has a OK callback button which does a 
      // g.place_typed_atom_at_pointer();
      GtkWidget *window = create_pointer_atom_type_dialog();
      
      GtkWidget *optionmenu = lookup_widget(window, "pointer_atom_molecule_optionmenu");
//       GtkSignalFunc callback_func =
// 	 GTK_SIGNAL_FUNC(graphics_info_t::pointer_atom_molecule_menu_item_activate);

      fill_place_atom_molecule_option_menu(optionmenu);

      gtk_widget_show(window);
   }
}

// This is a copy - more or less - of
// fill_option_menu_with_coordinates_options, except we also add at
// the top "New Molecule" if a molecule by the name of "Pointer Atoms"
// is not found.
// 
// Note that we can't use fill_option_menu_with_coordinates_options
// and add to it because gtk_menu_set_active fails/is ignored.
// 
// fill_pointer_atom_molecule_option_menu
// 
void fill_place_atom_molecule_option_menu(GtkWidget *optionmenu) { 

   GtkSignalFunc callback_func = 
      GTK_SIGNAL_FUNC(graphics_info_t::pointer_atom_molecule_menu_item_activate);
   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(optionmenu));
   gtk_widget_destroy(menu);
   menu = gtk_menu_new();

   GtkWidget *menuitem;
   int pointer_atoms_mol = -1;

   for (int imol=0; imol<graphics_n_molecules(); imol++) {
      if (graphics_info_t::molecules[imol].has_model()) { 
	 if (graphics_info_t::molecules[imol].name_ == "Pointer Atoms") { 
	    pointer_atoms_mol = imol;
	 } 
      }
   }

   int menu_index = 0;

   if (pointer_atoms_mol == -1) { 
      // There were no pointer atoms so let's create "New Molecule" at
      // the top of the list.
      // 
      GtkWidget *menu_item = gtk_menu_item_new_with_label("New Molecule");
      int imol_new = -10;
      gtk_signal_connect(GTK_OBJECT(menu_item), "activate",
			 callback_func,
			 GINT_TO_POINTER(imol_new));
      gtk_menu_append(GTK_MENU(menu), menu_item); 
      gtk_widget_show(menu_item); 
      gtk_menu_set_active(GTK_MENU(menu), 0);
      menu_index = 0;
   }

   for (int imol=0; imol<graphics_n_molecules(); imol++) {
      if (graphics_info_t::molecules[imol].has_model()) {
	 std::string ss = graphics_info_t::int_to_string(imol);
	 ss += " " ;
	 int ilen = graphics_info_t::molecules[imol].name_.length();
	 int left_size = ilen-graphics_info_t::go_to_atom_menu_label_n_chars_max;
	 if (left_size <= 0) {
	    // no chop
	    left_size = 0;
	 } else {
	    // chop
	    ss += "...";
	 }
	 ss += graphics_info_t::molecules[imol].name_.substr(left_size, ilen);
	 menuitem =  gtk_menu_item_new_with_label (ss.c_str());
	 menu_index++;
	 gtk_signal_connect (GTK_OBJECT (menuitem), "activate",
			     callback_func,
			     GINT_TO_POINTER(imol));
	 gtk_menu_append(GTK_MENU(menu), menuitem); 
	 gtk_widget_show(menuitem);

	 // set any previously saved active position:
	 if (graphics_info_t::user_pointer_atom_molecule == imol) {
	    std::cout << "setting active menu item to "
		      << menu_index << std::endl;
	    gtk_menu_set_active(GTK_MENU(menu),menu_index);
	 }
      }
   }
   
   /* Link the new menu to the optionmenu widget */
   gtk_option_menu_set_menu(GTK_OPTION_MENU(optionmenu),
			    menu);
} 

void
place_typed_atom_at_pointer(const char *type) {
   graphics_info_t g;
   g.place_typed_atom_at_pointer(std::string(type));
}

void set_pointer_atom_is_dummy(int i) { 
   graphics_info_t::pointer_atom_is_dummy = i;
} 

      

void display_where_is_pointer() {
   graphics_info_t g;
   g.display_where_is_pointer();
}

// draw the baton?
void set_draw_baton(short int i) {
   graphics_info_t g;
   g.draw_baton_flag = i;
   if (i == 1)
      g.start_baton_here();
   graphics_draw();
}

// Mouse movement moves the baton not the view?
void set_baton_mode(short int i) {
   graphics_info_t::baton_mode = i;
}

void accept_baton_position() { 	/* put an atom at the tip and move baton */

   graphics_info_t g;
   g.accept_baton_position();

}

void baton_try_another() {
   graphics_info_t g;
   g.baton_try_another();
}

void shorten_baton() {
   graphics_info_t g;
   g.shorten_baton();
}

void lengthen_baton() {
   graphics_info_t g;
   g.lengthen_baton();
}

void baton_build_delete_last_residue() {

   graphics_info_t g;
   g.baton_build_delete_last_residue();

} 

void set_baton_build_params(int istart_resno, 
			   const char *chain_id, 
			   const char *backwards) { 

   graphics_info_t g;
   g.set_baton_build_params(istart_resno, chain_id, backwards); 
}

// User data has been placed in the window - we use it to get the
// molecule number.
void baton_mode_calculate_skeleton(GtkWidget *window) {

   int imol = -1;

   int *i;

   std::cout << "getting intermediate data in baton_mode_calculate_skeleton "
	     << std::endl;
   i = (int *) gtk_object_get_user_data(GTK_OBJECT(window));

   std::cout << "got intermediate int: " << i << " " << *i << std::endl;

   imol = *i;

   std::cout << "calculating map for molecule " << imol << std::endl;
   if (imol < graphics_info_t::n_molecules && imol >= 0) { 
      skeletonize_map(0, imol);
   }

}


/* Reverse the direction of a the fragment of the clicked on
   atom/residue.  A fragment is a consequitive range of residues -
   where there is a gap in the numbering, that marks breaks between
   fragments in a chain.  There also needs to be a distance break - if
   the CA of the next/previous residue is more than 5A away, that also
   marks a break. Thow away all atoms in fragment other than CAs*/
void reverse_direction_of_fragment(int imol, const char *chain_id, int resno) {

   if (is_valid_model_molecule(imol)) {
      // return 1 if we did it.
      int istatus = graphics_info_t::molecules[imol].reverse_direction_of_fragment(std::string(chain_id), resno);
      if (istatus)
	 graphics_draw();
   }
}





// -----------------------------------------------------------------------------
//                               Automutation stuff 
// -----------------------------------------------------------------------------
// 
short int progressive_residues_in_chain_check(const char *chain_id, int imol) {
   
   graphics_info_t g;
   if (imol < graphics_n_molecules()) {
      return g.molecules[imol].progressive_residues_in_chain_check_by_chain(chain_id);
   } else {
      std::cout << "no such molecule number in progressive_residues_in_chain_check\n";
      return 0;
   }
} 

// return success on residue type match
// success: 1, failure: 0.
int
mutate_internal(int ires_serial, const char *chain_id, int imol, std::string &target_res_type) {

   graphics_info_t g;
   int istate = 0;
   if (imol < graphics_n_molecules()) {
      istate = g.molecules[imol].mutate_single_multipart(ires_serial, chain_id, target_res_type);
      if (istate == 0) {
	 std::cout << "ERROR: got bad state in mutate_internal" << std::endl;
      }
      graphics_draw();
   }
   return istate;
}

// causes a make_backup()
int
mutate(int ires, const char *chain_id, int imol, const char *target_res_type) {

   std::string target_type(target_res_type);

   std::string inscode("");
   if (imol < graphics_n_molecules()) {
      graphics_info_t::molecules[imol].mutate(ires, inscode, std::string(chain_id), std::string(target_res_type));
      graphics_draw();
   }
   return 0;
} 


// Return success on residue type match
// success: 1, failure: 0.
//
// Does not cause a make_backup().
//
int
mutate_single_residue_by_serial_number(int ires, const char *chain_id, int imol,
			   char target_res_type) {

   std::string target_as_str = coot::util::single_letter_to_3_letter_code(target_res_type);
   std::cout << "INFO:: mutate target_res_type :" << target_as_str << ":" << std::endl;
      
   return mutate_internal(ires, chain_id, imol, target_as_str);

}

// Previously, I was using mutate_single_residue_by_seqno to be a
// wrapper for mutate_single_residue_by_serial_number.
//
// But that fails when the residue is perfectly reasonably except that
// the serial number is -1 (I don't know wny this happens but it does
// in terminal residue addition).  So I will need new functionally
// that does the residue at root by seqnum not serial_number.
// 
int mutate_single_residue_by_seqno(int ires, const char *inscode,
				   const char *chain_id, 
				   int imol, char target_res_type) { 

   int status = -1; 
   std::string target_as_str = coot::util::single_letter_to_3_letter_code(target_res_type);
   
   if (imol < graphics_n_molecules()) {
      if (imol >= 0) { 
	 status = graphics_info_t::molecules[imol].mutate(ires,
							  std::string(inscode),
							  std::string(chain_id),
							  target_as_str);
      }
   }
   return status;
}


// return -1 on error:
// 
int chain_n_residues(const char *chain_id, int imol) {

   graphics_info_t g;
   if (imol < graphics_n_molecules()) {
      return g.molecules[imol].chain_n_residues(chain_id);
   } else { 
      return -1;
   }
}

// Return NULL (#f) on failure.
// 
char *resname_from_serial_number(int imol, const char *chain_id, int serial_num) {

   char *r = NULL;
   if (is_valid_model_molecule(imol)) {
      CMMDBManager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      int nchains = mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) {
	 CChain *chain_p = mol->GetChain(1,ichain);
	 std::string mol_chain_id(chain_p->GetChainID());
	 if (mol_chain_id == std::string(chain_id)) {
	    int nres = chain_p->GetNumberOfResidues();
	    if (serial_num < nres) {
	       int ch_n_res;
	       PCResidue *residues;
	       chain_p->GetResidueTable(residues, ch_n_res);
	       CResidue *this_res = residues[serial_num];
	       r = this_res->GetResName();
	    }
	 }
      }
   }
   return r;

}

// Return < -9999 on failure
int  seqnum_from_serial_number(int imol, const char *chain_id, int serial_num) {

   int iseqnum = -10000;
   if (is_valid_model_molecule(imol)) {
      CMMDBManager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      int nchains = mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) {
	 CChain *chain_p = mol->GetChain(1,ichain);
	 std::string mol_chain_id(chain_p->GetChainID());
	 if (mol_chain_id == std::string(chain_id)) {
	    int nres = chain_p->GetNumberOfResidues();
	    if (serial_num < nres) {
	       int ch_n_res;
	       PCResidue *residues;
	       chain_p->GetResidueTable(residues, ch_n_res);
	       CResidue *this_res = residues[serial_num];
	       iseqnum = this_res->GetSeqNum();
	    }
	 }
      }
   }
   return iseqnum;
}

char *insertion_code_from_serial_number(int imol, const char *chain_id, int serial_num) {

   char *r = NULL;
   if (is_valid_model_molecule(imol)) {
      CMMDBManager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      int nchains = mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) {
	 CChain *chain_p = mol->GetChain(1,ichain);
	 std::string mol_chain_id(chain_p->GetChainID());
	 if (mol_chain_id == std::string(chain_id)) {
	    int nres = chain_p->GetNumberOfResidues();
	    if (serial_num < nres) {
	       int ch_n_res;
	       PCResidue *residues;
	       chain_p->GetResidueTable(residues, ch_n_res);
	       CResidue *this_res = residues[serial_num];
	       r = this_res->GetInsCode();
	    }
	 }
      }
   }
   return r;
}


char *chain_id(int imol, int ichain) {

   char *r = NULL;
   if (is_valid_model_molecule(imol)) {
      CMMDBManager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      CChain *chain_p = mol->GetChain(1,ichain);
      r = chain_p->GetChainID();
   }
   return r;
}


int n_chains(int imol) {

   int nchains = -1;
   if (is_valid_model_molecule(imol)) {
      CMMDBManager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      nchains = mol->GetNumberOfChains(1);
   }
   return nchains;
}

// return -1 on error (e.g. chain_id not found, or molecule is not a
// model), 0 for no, 1 for is.
int is_solvent_chain_p(int imol, const char *chain_id) {

   int r = -1;
   if (is_valid_model_molecule(imol)) {
      CMMDBManager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      int nchains = mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) {
	 CChain *chain_p = mol->GetChain(1,ichain);
	 std::string mol_chain_id(chain_p->GetChainID());
	 if (mol_chain_id == std::string(chain_id)) {
	    r = chain_p->isSolventChain();
	 }
      }
   }
   return r;
}


/*  ----------------------------------------------------------------------- */
/*                         Renumber residue range                           */
/*  ----------------------------------------------------------------------- */

int renumber_residue_range(int imol, const char *chain_id,
			   int start_res, int last_res, int offset) {

   int i=0;
   if (imol >= 0) {
      if (imol <= graphics_info_t::n_molecules) {
	 if (graphics_info_t::molecules[imol].has_model()) {
	    i = graphics_info_t::molecules[imol].renumber_residue_range(chain_id,
									start_res,
									last_res,
									offset);
	    if (i) {
	       graphics_info_t g;
	       graphics_draw();
	       g.update_go_to_atom_window_on_changed_mol(imol);
	    }
	 }
      }
   }
   return i;
}

GtkWidget *wrapped_create_renumber_residue_range_dialog() {

   GtkWidget *w = create_renumber_residue_range_dialog();
   int imol = first_coords_imol();
   graphics_info_t::renumber_residue_range_molecule = imol;
   if (graphics_info_t::renumber_residue_range_molecule >= 0) { 
      graphics_info_t::fill_renumber_residue_range_dialog(w);
      graphics_info_t g;
      g.fill_renumber_residue_range_internal(w, imol);
   }
   

   return w;
}

void renumber_residues_from_widget(GtkWidget *window) {

   int imol = graphics_info_t::renumber_residue_range_molecule;

   GtkWidget *e1 = lookup_widget(window, "renumber_residue_range_resno_1_entry");
   GtkWidget *e2 = lookup_widget(window, "renumber_residue_range_resno_2_entry");
   GtkWidget *offent = lookup_widget(window, "renumber_residue_range_offset_entry");
   

   std::pair<short int, int> r1  = int_from_entry(e1);
   std::pair<short int, int> r2  = int_from_entry(e2);
   std::pair<short int, int> off = int_from_entry(offent);

   if (r1.first && r2.first && off.first) {
      int start_res = r1.second;
      int last_res =  r2.second;
      int offset   = off.second;

      if (imol >= 0) {
	 if (imol < graphics_info_t::n_molecules) {
	    if (graphics_info_t::molecules[imol].has_model()) {
	       std::string chain = graphics_info_t::renumber_residue_range_chain;
	       
	       renumber_residue_range(imol, chain.c_str(),
				      start_res, last_res, offset);
	       
	    }
	 }
      }
   } else {
      std::cout << "Sorry. Couldn't read residue or offset from entry widget\n";
   } 
}

int change_residue_number(int imol, const char *chain_id, int current_resno, const char *current_inscode, int new_resno, const char *new_inscode) {

   int idone = -1;
   if (is_valid_model_molecule(imol)) {
      std::string chain_id_str(chain_id);
      std::string current_inscode_str(current_inscode);
      std::string new_inscode_str(new_inscode);
      graphics_info_t::molecules[imol].change_residue_number(chain_id_str, current_resno, current_inscode_str, new_resno, new_inscode_str);
      graphics_draw();
      idone = 1;
   } 
   return idone;

} 


/*  ----------------------------------------------------------------------- */
/*                         backup                                           */
/*  ----------------------------------------------------------------------- */

void make_backup(int imol) {

   if (imol<graphics_info_t::n_molecules) {
      if (imol >= 0) {
	 if (graphics_info_t::molecules[imol].has_model()) {
	    graphics_info_t::molecules[imol].make_backup_from_outside();
	 } else {
	    std::cout << "No model for this molecule" << std::endl;
	 } 
      } else {
	 std::cout << "No model :" << imol << std::endl;
      }
   } else {
      std::cout << "No model :" << imol << std::endl;
   }
}

int backup_state(int imol) {

   int istate = -1;

   if (imol<graphics_info_t::n_molecules) {
      if (imol >= 0) {
	 if (graphics_info_t::molecules[imol].has_model()) {
	    istate = graphics_info_t::molecules[imol].backups_state();
	 } else {
	    std::cout << "No model for this molecule" << std::endl;
	 } 
      } else {
	 std::cout << "No model :" << imol << std::endl;
      }
   } else {
      std::cout << "No model :" << imol << std::endl;
   }
   return istate;
} 

void set_have_unsaved_changes(int imol) {

   if (imol<graphics_info_t::n_molecules) {
      if (imol >= 0) {
	 if (graphics_info_t::molecules[imol].has_model()) {
	    graphics_info_t::molecules[imol].set_have_unsaved_changes_from_outside();
	 }
      }
   }
}


/*  ------------------------------------------------------------------------ */
/*                         Write PDB file:                                   */
/*  ------------------------------------------------------------------------ */

// return 1 on error, 0 on success
int
write_pdb_file(int imol, const char *file_name) {

   graphics_info_t g;
   int istat = 0;
   if (is_valid_model_molecule(imol)) {
      istat = g.molecules[imol].write_pdb_file(std::string(file_name));
   }
   return istat;
}

/*! \brief write molecule number imol's residue range as a PDB to file
  file_name */
/*  return 0 on success, 1 on error. */
int
write_residue_range_to_pdb_file(int imol, const char *chainid, 
				int resno_start, int resno_end,
				const char *filename) {

   int istat = 0;
   if (is_valid_model_molecule(imol)) {
      std::string chain(chainid);
      if (resno_end < resno_start) {
	 int tmp = resno_end;
	 resno_end = resno_start;
	 resno_start = tmp;
      } 
      CMMDBManager *mol =
	 graphics_info_t::molecules[imol].get_residue_range_as_mol(chain, resno_start, resno_end);
      if (mol) {
	 istat = mol->WritePDBASCII((char *)filename);
	 delete mol; // give back the memory.
      }
   }
   return istat;
} 


/*  ------------------------------------------------------------------------ */
/*                         refmac stuff                                      */
/*  ------------------------------------------------------------------------ */

void execute_refmac(GtkWidget *window) {  /* lookup stuff here. */

   // The passed window, is the refmac dialog, where one selects the
   // coords molecule and the map molecule.

   std::cout << "DEUBG here 1" << std::endl;

   GtkWidget *option_menu = lookup_widget(window,
					  "run_refmac_coords_optionmenu");

   GtkWidget *menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));

   GtkWidget *active_item = gtk_menu_get_active(GTK_MENU(menu));

   int imol_coords = graphics_info_t::refmac_molecule; 
   if (imol_coords < 0) {
      std::cout << "No coordinates selected for refmac\n";
   } else { 

//       std::cout << " Running refmac coords molecule number "
// 		<< imol_coords << std::endl;

      option_menu = lookup_widget(window, "run_refmac_map_optionmenu");
      menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(option_menu));
      active_item = gtk_menu_get_active(GTK_MENU(menu));

      // active_item is set if there was at least one map with refmac params:
      // if none, it is null.
      
      if (active_item == 0) {
	 add_status_bar_text("No map has associated Refmac Parameters - no REFMAC!");
      } else { 
	 int imol_window = GPOINTER_TO_INT(gtk_object_get_user_data(GTK_OBJECT(active_item)));
	 if (imol_window < 0) {
	    std::cout << "No map data selected for refmac\n";
	 } else { 
	 
	    int imol_map_refmac = imol_window;
	    if (!is_valid_map_molecule(imol_map_refmac)) {
	       std::string s = "Invalid molecule number: ";
	       s += graphics_info_t::int_to_string(imol_map_refmac);
	       std::cout << s << std::endl;
	       graphics_info_t g;
	       g.statusbar_text(s);
	    } else {
	       // normal path
	       if (graphics_info_t::molecules[imol_map_refmac].Have_sensible_refmac_params()) { 

		  std::cout << " Running refmac refmac params molecule number "
			    << imol_map_refmac << std::endl;

		  graphics_info_t g;

		  std::string refmac_dir("coot-refmac");
		  short int have_ccp4i_project = 0;
		  if (graphics_info_t::refmac_ccp4i_project_dir != "") { 
		     refmac_dir = graphics_info_t::refmac_ccp4i_project_dir;
		     have_ccp4i_project = 1;
		  }
		  int istat = make_directory_maybe(refmac_dir.c_str());
		  if (istat != 0) { // fails
		     std::cout << "WARNING failed to make directory for refmac -"
			       << " run refmac fails\n" << std::endl;
		  } else {

		     // now lookup the active state of the difference map and
		     // the phase combine buttons:
		     //
		     int diff_map_flag;
		     int phase_combine_flag;
		     GtkWidget *checkbutton =
			lookup_widget(window, "run_refmac_phase_combine_checkbutton");
		     if (GTK_TOGGLE_BUTTON(checkbutton)->active) {
			phase_combine_flag = 1;
		     } else {
			phase_combine_flag = 0;
		     }
	       
		     checkbutton =  lookup_widget(window,"run_refmac_diff_map_checkbutton");
		     if (GTK_TOGGLE_BUTTON(checkbutton)->active) {
			diff_map_flag = 1;
		     } else {
			diff_map_flag = 0;
		     }


		     // g.molecules[imol_coords].increment_refmac_count();
      
		     std::string pdb_in_filename  = refmac_dir;
		     std::string pdb_out_filename = refmac_dir;
		     std::string mtz_out_filename = refmac_dir;
		     std::cout << "DEBUG:: pdb_in_filename is now 1 " <<  pdb_in_filename << std::endl;
		     if (! have_ccp4i_project) { 
			pdb_in_filename  += "/";
			pdb_out_filename += "/";
			mtz_out_filename += "/";
		     }
		     std::cout << "DEBUG:: pdb_in_filename is now 2 " <<  pdb_in_filename << std::endl;
		     pdb_in_filename += g.molecules[imol_coords].Refmac_in_name();
		     std::cout << "DEBUG:: pdb_in_filename is now 3 " <<  pdb_in_filename << std::endl;

		     // cleverness happens in Refmac_out_name:
		     pdb_out_filename += g.molecules[imol_coords].Refmac_out_name();
		     mtz_out_filename += g.molecules[imol_coords].Refmac_mtz_out_name();

		     std::string mtz_in_filename = g.molecules[imol_map_refmac].Refmac_mtz_filename();
		     std::string refmac_count_string =
			g.int_to_string(g.molecules[imol_coords].Refmac_count());

		     std::cout << "DEBUG:: mtz_out_filename: " <<
			mtz_out_filename << std::endl;
		     std::cout << "DEBUG:: pdb_out_filename: " <<
			pdb_out_filename << std::endl;

		     std::string phib_string;
		     std::string fom_string;

		     if (g.molecules[imol_map_refmac].Fourier_weight_label() != "") {
			phib_string = g.molecules[imol_map_refmac].Fourier_phi_label();
			fom_string  = g.molecules[imol_map_refmac].Fourier_weight_label();
		     } else {
			std::cout << "WARNING:: Can't do phase combination if we don't use FOMs ";
			std::cout << "to make the map" << std::endl;
			std::cout << "WARNING:: Turning off phase combination." << std::endl;
			phase_combine_flag = 0; 
		     }
		     // 	    std::cout << "DEBUG:: fom_string " << fom_string << " "
		     // 		      << g.molecules[imol_map_refmac].Fourier_weight_label()
		     // 		      << std::endl;

		     std::string cif_lib_filename = ""; // default, none
		     if (graphics_info_t::cif_dictionary_filename_vec->size() > 0) {
			cif_lib_filename = (*graphics_info_t::cif_dictionary_filename_vec)[0];
		     }

		     // 	    std::cout << "DEBUG:: attempting to write pdb input file "
		     // 		      << pdb_in_filename << std::endl;
		     int ierr = g.molecules[imol_coords].write_pdb_file(pdb_in_filename);
		     if (!ierr) { 
			std::cout << "refmac ccp4i project dir " 
				  << graphics_info_t::refmac_ccp4i_project_dir 
				  << std::endl;
			execute_refmac_real(pdb_in_filename, pdb_out_filename,
					    mtz_in_filename, mtz_out_filename,
					    cif_lib_filename,
					    g.molecules[imol_map_refmac].Refmac_fobs_col(),
					    g.molecules[imol_map_refmac].Refmac_sigfobs_col(),
					    g.molecules[imol_map_refmac].Refmac_r_free_col(),
					    g.molecules[imol_map_refmac].Refmac_r_free_sensible(),
					    refmac_count_string,
					    g.swap_pre_post_refmac_map_colours_flag,
					    imol_map_refmac,
					    diff_map_flag,
					    phase_combine_flag, phib_string, fom_string,
					    graphics_info_t::refmac_ccp4i_project_dir);
		     } else {
			std::cout << "WARNING:: fatal error in writing pdb input file"
				  << pdb_in_filename << " for refmac.  Can't run refmac"
				  << std::endl;
		     } 
		  }
	       }
	    }
	 }
      }
   }
}

//       int slen = mtz_in_filename.length(); c
//       if (slen > 4) {
// 	 mtz_out_filename = mtz_in_filename.substr(0,slen - 4) + "-refmac-";
// 	 mtz_out_filename += g.int_to_string(g.molecules[imol_coords].Refmac_count());
// 	 mtz_out_filename += ".mtz";
//       } else {
// 	 mtz_out_filename = "post-refmac";
// 	 mtz_out_filename += g.int_to_string(g.molecules[imol_coords].Refmac_count());
// 	 mtz_out_filename += ".mtz";
//       } 

// If ccp4i_project_dir is "", then carry on and put the log file in
// this directory.  If not, put it in the appropriate project dir. The
// pdb_in etc filename are manipulated in the calling routine.
//
// if swap_map_colours_post_refmac_flag is not 1 then imol_refmac_map is ignored.
// 
void
execute_refmac_real(std::string pdb_in_filename,
		    std::string pdb_out_filename,
		    std::string mtz_in_filename,
		    std::string mtz_out_filename,
		    std::string cif_lib_filename,
		    std::string fobs_col_name,
		    std::string sigfobs_col_name,
		    std::string r_free_col_name,
		    short int have_sensible_free_r_flag,
		    std::string refmac_count_str,
		    int swap_map_colours_post_refmac_flag,
		    int imol_refmac_map,
		    int diff_map_flag,
		    int phase_combine_flag,
		    std::string phib_string,
		    std::string fom_string, 
		    std::string ccp4i_project_dir) {


   std::cout << "DEUBG here 2 in execute_refmac_real " << std::endl;

   std::vector<std::string> cmds;

   cmds.push_back(std::string("run-refmac-by-filename"));
   cmds.push_back(single_quote(pdb_in_filename));
   cmds.push_back(single_quote(pdb_out_filename));
   cmds.push_back(single_quote(mtz_in_filename));
   cmds.push_back(single_quote(mtz_out_filename));
   cmds.push_back(single_quote(cif_lib_filename));
   cmds.push_back(refmac_count_str);
   cmds.push_back(graphics_info_t::int_to_string(swap_map_colours_post_refmac_flag));
   cmds.push_back(graphics_info_t::int_to_string(imol_refmac_map));
   cmds.push_back(graphics_info_t::int_to_string(diff_map_flag));
   cmds.push_back(graphics_info_t::int_to_string(phase_combine_flag));

   std::string phase_combine_cmd;
   if (phase_combine_flag) {
      phase_combine_cmd += "(cons \"";
      phase_combine_cmd += phib_string;
      phase_combine_cmd += "\" \"";
      phase_combine_cmd += fom_string;
      phase_combine_cmd += "\")";
      phase_combine_cmd += " ";
   } else {
      phase_combine_cmd += "'dummy ";
   }
   cmds.push_back(phase_combine_cmd);

   cmds.push_back(graphics_info_t::int_to_string(-1)); // don't use NCYCLES
   cmds.push_back(single_quote(ccp4i_project_dir));
   cmds.push_back(single_quote(fobs_col_name));
   cmds.push_back(single_quote(sigfobs_col_name));

   if (have_sensible_free_r_flag) { 
      cmds.push_back(single_quote(r_free_col_name));
   }

   graphics_info_t g;
   short int ilang = coot::STATE_SCM;
   std::string cmd;

#ifdef USE_PYTHON
#ifndef USE_GUILE
   ilang = coot::STATE_PYTHON;
#endif
#endif
   if (ilang == coot::STATE_PYTHON) { 
      cmd = g.state_command(cmds, ilang);
      safe_python_command(cmd);
   } else {
      cmd = g.state_command(cmds, ilang);
      safe_scheme_command(cmd);
   } 

#ifdef USE_PYTHON

   std::cout << "put construction of python command code here\n"; 

#endif // USE_PYTHON   
} 

// Not needed? because we look at the active menu item at OK button-press time?
//
// Well, that was indeeed the way that we used to do it, now (WDW)
// that we rationalize the coordinates molecules option menu filling
// we have to set the the active molecule in this callback.
// 
void
refmac_molecule_button_select(GtkWidget *item, GtkPositionType pos) {

   graphics_info_t::refmac_molecule = pos;

}

int set_refmac_molecule(int imol) {
   graphics_info_t::refmac_molecule = imol;
   return imol;
}


void fill_option_menu_with_refmac_options(GtkWidget *optionmenu) {

   graphics_info_t g;
   g.fill_option_menu_with_refmac_options(optionmenu);

} 

void set_refmac_counter(int imol, int refmac_count) {

   graphics_info_t g;
   if (imol< g.n_molecules) {
      g.molecules[imol].set_refmac_counter(refmac_count);
      std::cout << "INFO:: refmac counter of molecule number " << imol
		<< " incremented to " << refmac_count << std::endl;
   } else {
      std::cout << "WARNING:: refmac counter of molecule number " << imol
		<< " not incremented to " << refmac_count << std::endl;
   } 

} 


const char *refmac_name(int imol) {

   graphics_info_t g;
   return g.molecules[imol].Refmac_in_name().c_str();

} 


short int 
add_OXT_to_residue(int imol, int reso, const char *insertion_code, const char *chain_id) {

   short int istat = -1; 
   if (imol < graphics_n_molecules()) { 
      istat = graphics_info_t::molecules[imol].add_OXT_to_residue(reso, std::string(insertion_code),
								  std::string(chain_id));
      graphics_info_t::molecules[imol].update_symmetry();
      graphics_draw();
   }
   return istat;
}

void apply_add_OXT_from_widget(GtkWidget *w) {

   int imol = graphics_info_t::add_OXT_molecule;
   int resno = -9999;
   std::string chain_id = graphics_info_t::add_OXT_chain;

   GtkWidget *terminal_checkbutton = lookup_widget(w, "add_OXT_c_terminus_radiobutton");
   GtkWidget *residue_number_entry = lookup_widget(w, "add_OXT_residue_entry");

   if (GTK_TOGGLE_BUTTON(terminal_checkbutton)->active) {
      std::cout << "DEBUG:: auto determine C terminus..." << std::endl;
      // we need to determine the last residue in this chain:
      if (imol >= 0) {
	 if (imol < graphics_info_t::n_molecules) {
	    if (graphics_info_t::molecules[imol].has_model()) {
	       std::pair<short int, int> p =
		  graphics_info_t::molecules[imol].last_residue_in_chain(chain_id);
	       if (p.first) {
		  resno = p.second;
	       } 
	    }
	 }
      }
   } else {
      // we get the resno from the widget
      std::pair<short int, int> p = int_from_entry(residue_number_entry);
      if (p.first) {
	 resno = p.second;
      }
   }

   if (resno > -9999) { 
      if (imol >= 0) {
	 if (imol < graphics_info_t::n_molecules) {
	    if (graphics_info_t::molecules[imol].has_model()) { 
	       std::cout << "DEBUG:: adding OXT to " << imol << " "
			 << chain_id << " " << resno << std::endl;
	       
	       add_OXT_to_residue(imol, resno, "", chain_id.c_str());
	    }
	 }
      }
   } else {
      std::cout << "WARNING:: Could not determine last residue - not adding OXT\n";
   } 
}


GtkWidget *wrapped_create_add_OXT_dialog() {

   GtkWidget *w = create_add_OXT_dialog();

   GtkWidget *option_menu = lookup_widget(w, "add_OXT_molecule_optionmenu");

   GtkSignalFunc callback_func = GTK_SIGNAL_FUNC(graphics_info_t::add_OXT_molecule_item_select);

   graphics_info_t g;
   int imol = first_coords_imol();
   graphics_info_t::add_OXT_molecule = imol;
   g.fill_option_menu_with_coordinates_options(option_menu, callback_func, imol);

   g.fill_add_OXT_dialog_internal(w, imol); // function needs object (not static)

   return w;
}



/*  ------------------------------------------------------------------------ */
/*                         db_mainchain:                                     */
/*  ------------------------------------------------------------------------ */

// The button callback:
void 
do_db_main(short int state) {
   graphics_info_t g;
   g.in_db_main_define = state;
   if (state) { 
      g.pick_cursor_maybe();
      g.pick_pending_flag = 1;
      std::cout << "click on 2 atoms to define a the range of residues to"
		<< " convert to mainchain" << std::endl;
   } else {
      g.pick_pending_flag = 0;
      g.normal_cursor();
   }
}


void
db_mainchain(int imol,
	     const char *chain_id,
	     int iresno_start,
	     int iresno_end,
	     const char *direction_string) {
   
   if (imol < graphics_n_molecules()) {
      graphics_info_t g;
      g.execute_db_main(imol, std::string(chain_id), iresno_start, iresno_end,
			std::string(direction_string));
   } else {
      std::cout << "WARNING molecule index error" << std::endl;
   } 
}

/*  ----------------------------------------------------------------------- */
/*                  alternate conformation                                  */
/*  ----------------------------------------------------------------------- */

// so shall we define a nomenclature for the split type:
// 0: split at Ca
// 1: split whole residue
// 2: split residue range

/* c-interface-build function */
short int alt_conf_split_type_number() { 
   return graphics_info_t::alt_conf_split_type_number();
}

void setup_alt_conf_with_dialog(GtkWidget *dialog) {

   GtkWidget *widget_ca = lookup_widget(dialog, 
					"add_alt_conf_ca_radiobutton");
   GtkWidget *widget_whole = lookup_widget(dialog, 
					   "add_alt_conf_whole_single_residue_radiobutton");
   GtkWidget *widget_range = lookup_widget(dialog, 
					   "add_alt_conf_residue_range_radiobutton");

   if (graphics_info_t::alt_conf_split_type_number() == 0)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widget_ca), TRUE);
   if (graphics_info_t::alt_conf_split_type_number() == 1)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widget_whole), TRUE);
   if (graphics_info_t::alt_conf_split_type_number() == 2)
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(widget_range), TRUE);

   if (graphics_info_t::alt_conf_split_type_number() < 2) { 
      std::cout << "Click on the residue you want to split" << std::endl;
   } else { 
      std::cout << "Click on a residue range you want to split" << std::endl;
   }

   graphics_info_t::in_add_alt_conf_define = 1;
   graphics_info_t::pick_cursor_maybe();
   graphics_info_t::pick_pending_flag = 1;
   graphics_info_t::add_alt_conf_dialog = dialog;
} 


void set_add_alt_conf_split_type_number(short int i) { 
   graphics_info_t::alt_conf_split_type = i;
} 

void unset_add_alt_conf_dialog()  { /* set the static dialog holder in
				     graphics info to NULL */
   graphics_info_t::add_alt_conf_dialog = NULL;
}

void unset_add_alt_conf_define() {  /* turn off pending atom pick */

   graphics_info_t g;
   graphics_info_t::in_add_alt_conf_define = 0;
   g.normal_cursor();
}


void altconf() { 

   GtkWidget *widget = create_add_alt_conf_dialog();
   setup_alt_conf_with_dialog(widget);
   gtk_widget_show(widget);
}

void set_show_alt_conf_intermediate_atoms(int i) {
   graphics_info_t::show_alt_conf_intermediate_atoms_flag = i;
}

int show_alt_conf_intermediate_atoms_state() {
   return graphics_info_t::show_alt_conf_intermediate_atoms_flag;
}

void zero_occupancy_residue_range(int imol, const char *chain_id, int ires1, int ires2) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_occupancy_residue_range(std::string(chain_id), ires1, ires2, 0.0);
   } else {
      std::cout << "WARNING:: invalid model molecule number in zero_occupancy_residue_range "
		<< imol << std::endl;
   }
   graphics_draw();
}

void fill_occupancy_residue_range(int imol, const char *chain_id, int ires1, int ires2) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_occupancy_residue_range(std::string(chain_id), ires1, ires2, 1.0);
   } else {
      std::cout << "WARNING:: invalid model molecule number in fill_occupancy_residue_range "
		<< imol << std::endl;
   }
   graphics_draw();
}


/*  ----------------------------------------------------------------------- */
/*                  Edit Chi Angles                                         */
/*  ----------------------------------------------------------------------- */
void setup_edit_chi_angles(short int state) {

   graphics_info_t g;

   if (state) { 
      g.in_edit_chi_angles_define = 1;
      std::cout << "Click on an atom in the residue that you want to edit" << std::endl;
      g.pick_cursor_maybe();
      g.statusbar_text("Click on a atom. The clicked atom affects the torsion's wagging dog/tail...");
      g.pick_pending_flag = 1;
   } else {
      g.in_edit_chi_angles_define = 0;
   }
}

int set_show_chi_angle_bond(int imode) {

   graphics_info_t::draw_chi_angle_flash_bond_flag = imode;
   graphics_draw();
   return 0; // should be a void function, I imagine.
} 



// Set a flag: Should torsions that move hydrogens be
// considered/displayed in button box?
// 
void set_find_hydrogen_torsion(short int state) {
   graphics_info_t g;
   g.find_hydrogen_torsions = state;
}

void set_graphics_edit_current_chi(int ichi) { /* button callback */

   graphics_info_t::edit_chi_current_chi = 0;
   graphics_info_t::in_edit_chi_mode_flag = 0; // off
}

void unset_moving_atom_move_chis() { 
   graphics_info_t::moving_atoms_move_chis_flag = 0; // keyboard 1,2,3
						     // etc cant put
						     // graphics/mouse
						     // into rotate
						     // chi mode.
} 


// altconf is ignored here
int edit_chi_angles(int imol, const char *chain_id, int resno, 
		    const char *ins_code, const char *altconf) {

   int status = 0;

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      // type void
      int atom_index = atom_index_first_atom_in_residue(imol, chain_id, resno, ins_code);
      if (atom_index > -1) {
	 g.execute_edit_chi_angles(atom_index, imol);
	 status = 1;
      }
   }
   return status;
}



/*  ------------------------------------------------------------------------ */
/*                         recover session:                                  */
/*  ------------------------------------------------------------------------ */
/* section Recover Session Function */
/* After a crash (shock horror!) we provide this convenient interface
   to restore the session.  It runs through all the molecules with
   models and looks at the coot backup directory looking for related
   backup files that are more recent that the read file. */
void recover_session() { 

   int i_rec = 0;
   for (int imol=0; imol<graphics_info_t::n_molecules; imol++) { 
      if (graphics_info_t::molecules[imol].has_model()) { 
	 coot::backup_file_info info = 
	    graphics_info_t::molecules[imol].recent_backup_file_info();
	 if (info.status) { 

	    coot::backup_file_info *info_copy = new coot::backup_file_info;
	    *info_copy = info;
	    info_copy->imol = imol;
	    
	    GtkWidget *widget = create_recover_coordinates_dialog();
	    gtk_object_set_user_data(GTK_OBJECT(widget), info_copy);
	    
	    GtkWidget *label1, *label2;
	    label1 = lookup_widget(widget, "recover_coordinates_read_coords_label");
	    label2 = lookup_widget(widget, "recover_coordinates_backup_coordinates_label");

	    gtk_label_set_text(GTK_LABEL(label1), info.name.c_str());
	    gtk_label_set_text(GTK_LABEL(label2), info.backup_file_name.c_str());

	    gtk_widget_show(widget);
	    i_rec++;
	 }
      }
   }
   if (i_rec == 0) {
      GtkWidget *w = create_nothing_to_recover_dialog();
      gtk_widget_show(w);
   }
}

// widget needed for lookup of user data:
// 
void execute_recover_session(GtkWidget *widget) { 

   coot::backup_file_info *info = (coot::backup_file_info *) gtk_object_get_user_data(GTK_OBJECT(widget));

   if (info) { 
      
      graphics_info_t g;
      if (info->imol >= 0 && info->imol < g.n_molecules) { 
	 g.molecules[info->imol].execute_restore_from_recent_backup(info->backup_file_name);
	 graphics_draw();
      }
   } else { 
      std::cout << "ERROR:: couldn't find user data in execute_recover_session\n";
   } 
} 


void translate_molecule_by(int imol, float x, float y, float z) { 

   if (imol < graphics_info_t::n_molecules && imol >= 0) { 
      if (graphics_info_t::molecules[imol].has_model()) { 
	 graphics_info_t::molecules[imol].translate_by(x, y, z);
      }
   }
   graphics_draw();
} 



void assign_fasta_sequence(int imol, const char *chain_id_in, const char *seq) { 

   // format "> name \n <sequence>"
   if (imol < graphics_info_t::n_molecules && imol >= 0) {
      const std::string chain_id = chain_id_in;
      graphics_info_t::molecules[imol].assign_fasta_sequence(chain_id, std::string(seq));
   }
} 

/*  ----------------------------------------------------------------------- */
/*                  trim                                                    */
/*  ----------------------------------------------------------------------- */
void 
trim_molecule_by_map(int imol_coords, int imol_map, 
		     float map_level, int delete_or_zero_occ_flag) {

   graphics_info_t g;
   if (imol_coords < g.n_molecules) { 
      if (imol_map < g.n_molecules) {
	 if (g.molecules[imol_map].has_map()) { 
	    int iv = g.molecules[imol_coords].trim_by_map(g.molecules[imol_map].xmap_list[0], 
							  map_level,
							  delete_or_zero_occ_flag);
	    if (iv) 
	       graphics_draw();
	 } else { 
	    std::cout << "molecule " << imol_map << " has no map" << std::endl;
	 } 
      } else { 
	 std::cout << "No such molecule for map as " << imol_map << std::endl;
      } 
   } else { 
      std::cout << "No such molecule for model as " << imol_coords << std::endl;
   } 
}



/*  ----------------------------------------------------------------------- */
/*               Simplex Refinement                                         */
/*  ----------------------------------------------------------------------- */
//
// Perhaps this should be a just a call to a graphics_info_t function?
// 
void
fit_residue_range_to_map_by_simplex(int res1, int res2, char *altloc,
				    char *chain_id, int imol, int imol_for_map) {


   // The molecule_class_info_t updates its bonds.
   // 
   if (imol < graphics_info_t::n_molecules) {
      if (graphics_info_t::molecules[imol].has_model()) {
	 if (imol_for_map < graphics_info_t::n_molecules) {
	    if (graphics_info_t::molecules[imol_for_map].has_map()) { 
	       graphics_info_t::molecules[imol].fit_residue_range_to_map_by_simplex(res1, res2, altloc, chain_id, imol_for_map);
	    } else {
	       std::cout << "No map for molecule " << imol_for_map << std::endl;
	    }
	 } else {
	    std::cout << "No molecule " << imol_for_map << std::endl;
	 }
      } else {
	 std::cout << "No coordinates for molecule " << imol << std::endl;
      } 
   } else {
      std::cout << "No molecule " << imol << std::endl;
   }

   graphics_draw();

} 

// Return a score of the fit to the map.
// 
float
score_residue_range_fit_to_map(int res1, int res2, char *altloc,
			       char *chain_id, int imol, int imol_for_map) {

   float f = 0.0;

   if (imol < graphics_info_t::n_molecules) {
      if (graphics_info_t::molecules[imol].has_model()) {
	 if (imol_for_map < graphics_info_t::n_molecules) {
	    if (graphics_info_t::molecules[imol_for_map].has_map()) { 
	       f = graphics_info_t::molecules[imol].score_residue_range_fit_to_map(res1, res2, altloc, chain_id, imol_for_map);
	    } else {
	       std::cout << "No map for molecule " << imol_for_map << std::endl;
	    }
	 } else {
	    std::cout << "No molecule " << imol_for_map << std::endl;
	 }
      } else {
	 std::cout << "No coordinates for molecule " << imol << std::endl;
      } 
   } else {
      std::cout << "No molecule " << imol << std::endl;
   }
   return f;
}
									 
/*  ----------------------------------------------------------------------- */
/*                         Planar Peptide Restraints                        */
/*  ----------------------------------------------------------------------- */

void add_planar_peptide_restraints() {

   graphics_info_t g;
   g.Geom_p()->add_planar_peptide_restraint();
} 

void remove_planar_peptide_restraints() {

   graphics_info_t g;
   g.Geom_p()->remove_planar_peptide_restraint();
}

void add_omega_torsion_restriants() {

   graphics_info_t g;
   g.Geom_p()->add_omega_peptide_restraints();
} 

/*! \brief remove omega restraints on CIS and TRANS linked residues. */
void remove_omega_torsion_restriants() {

   graphics_info_t g;
   g.Geom_p()->remove_omega_peptide_restraints();
}



/*  ----------------------------------------------------------------------- */
/*                         Merge Molecules                                  */
/*  ----------------------------------------------------------------------- */

GtkWidget *wrapped_create_merge_molecules_dialog() {

   GtkWidget *w = create_merge_molecules_dialog();
   // fill the dialog here
   GtkWidget *molecule_option_menu = lookup_widget(w, "merge_molecules_optionmenu");
   GtkWidget *molecules_vbox       = lookup_widget(w, "merge_molecules_vbox");

   GtkSignalFunc callback_func = GTK_SIGNAL_FUNC(merge_molecules_menu_item_activate);
   GtkSignalFunc checkbox_callback_func = GTK_SIGNAL_FUNC(on_merge_molecules_check_button_toggled);


   fill_vbox_with_coordinates_options(molecules_vbox, checkbox_callback_func);

   int imol_master = graphics_info_t::merge_molecules_master_molecule;
   if (imol_master == -1) { 
      for (int i=0; i<graphics_info_t::n_molecules; i++) {
	 if (graphics_info_t::molecules[i].has_model()) {
	    graphics_info_t::merge_molecules_master_molecule = i;
	    imol_master = i;
	    break;
	 }
      }
   }

   graphics_info_t g;
   g.fill_option_menu_with_coordinates_options(molecule_option_menu,
					       callback_func, imol_master);
   return w;
}

void merge_molecules_menu_item_activate(GtkWidget *item, 
					GtkPositionType pos) {

   graphics_info_t::merge_molecules_master_molecule = pos;
}

void fill_vbox_with_coordinates_options(GtkWidget *dialog,
					GtkSignalFunc checkbox_callback_func) {

   GtkWidget *checkbutton;
   std::string button_label;
   GtkWidget *molecules_vbox = lookup_widget(dialog, "merge_molecules_vbox");

   // Unset any preconcieved notion of merging molecules:
   // 
   graphics_info_t::merge_molecules_merging_molecules->resize(0);

   for (int i=0; i<graphics_info_t::n_molecules; i++) {
      if (graphics_info_t::molecules[i].has_model()) {
	 button_label = graphics_info_t::int_to_string(i);
	 button_label += " ";
	 button_label += graphics_info_t::molecules[i].name_for_display_manager();
	 std::string button_name = "merge_molecules_checkbutton_";
	 button_name += graphics_info_t::int_to_string(i);

	 checkbutton = gtk_check_button_new_with_label(button_label.c_str());
  	 gtk_widget_ref (checkbutton);
  	 gtk_object_set_data_full (GTK_OBJECT (dialog),
  				   button_name.c_str(), checkbutton,
  				   (GtkDestroyNotify) gtk_widget_unref);
	 // The callback (if active) adds this molecule to the merging molecules list.
	 // If not active, it tries to remove it from the list.
	 //
	 // Why am I doing it like this instead of just looking at the
	 // state of the checkbutton when the OK button is pressed?
	 // Because (for the life of me) I can't seem to correctly
	 // lookup the checkbuttons from the button (or dialog for
	 // that matter).
	 // 
	 //  We look at the state when the
	 // "Merge" button is pressed - we don't need a callback to do
	 // that.
	 // 
  	 gtk_signal_connect (GTK_OBJECT (checkbutton), "toggled",
  			     GTK_SIGNAL_FUNC (checkbox_callback_func),
  			     GINT_TO_POINTER(i));
	 gtk_widget_show (checkbutton);
	 gtk_box_pack_start (GTK_BOX (molecules_vbox), checkbutton, FALSE, FALSE, 0);
	 gtk_container_set_border_width (GTK_CONTAINER (checkbutton), 2);
      }
   }
}

// The callback (if active) adds this molecule to the merging molecules list.
// If not active, it tries to remove it from the list.
// 
void on_merge_molecules_check_button_toggled (GtkToggleButton *togglebutton,
					      gpointer         user_data) {

   int imol = GPOINTER_TO_INT(user_data);
   if (togglebutton->active) {
      std::cout << "INFO:: adding molecule " << imol << " to merging list\n";
      graphics_info_t::merge_molecules_merging_molecules->push_back(imol);
   } else {
      std::cout << "INFO:: removing molecule " << imol << " from merging list\n";
      if (coot::is_member_p(*graphics_info_t::merge_molecules_merging_molecules, imol)) {
	 // passing a pointer
	 coot::remove_member(graphics_info_t::merge_molecules_merging_molecules, imol);
      }
   }
}


// Display the gui
void do_merge_molecules_gui() {

   GtkWidget *w = wrapped_create_merge_molecules_dialog();
   gtk_widget_show(w);
} 

// The action on Merge button press:
// 
void do_merge_molecules(GtkWidget *dialog) {

   std::vector<int> add_molecules = *graphics_info_t::merge_molecules_merging_molecules;
   if (add_molecules.size() > 0) { 
      int istat = merge_molecules(add_molecules, graphics_info_t::merge_molecules_master_molecule);
      if (istat)
	 graphics_draw();
   }
}


int
merge_molecules(const std::vector<int> &add_molecules, int imol) {

   int imerged = 0;
   std::vector<atom_selection_container_t> add_molecules_at_sels;
   if (imol < graphics_info_t::n_molecules) {
      if (imol >= 0) {
	 if (graphics_info_t::molecules[imol].has_model()) {

	    for (unsigned int i=0; i<add_molecules.size(); i++) {
	       if (add_molecules[i]<graphics_info_t::n_molecules) {
		  if (i >= 0) {
		     if (graphics_info_t::molecules[add_molecules[i]].has_model()) {
			if (add_molecules[i] != imol) { 
			   add_molecules_at_sels.push_back(graphics_info_t::molecules[add_molecules[i]].atom_sel);
			}
		     }
		  }
	       }
	    }
	    if (add_molecules_at_sels.size() > 0) { 
	       imerged = graphics_info_t::molecules[imol].merge_molecules(add_molecules_at_sels);
	    }
	 }
      }
   }
   return imerged;
}



/*  ----------------------------------------------------------------------- */
/*                         Mutate Sequence GUI                              */
/*  ----------------------------------------------------------------------- */

GtkWidget *wrapped_create_mutate_sequence_dialog() {

   GtkWidget *w = create_mutate_sequence_dialog();

   set_transient_and_position(COOT_MUTATE_RESIDUE_RANGE_WINDOW, w);

   GtkWidget *molecule_option_menu = lookup_widget(w, "mutate_molecule_optionmenu");
   GtkWidget *chain_option_menu    = lookup_widget(w, "mutate_molecule_chain_optionmenu");
//    GtkWidget *entry1 = lookup_widget(w, "mutate_molecule_resno_1_entry");
//    GtkWidget *entry2 = lookup_widget(w, "mutate_molecule_resno_2_entry");
//    GtkWidget *textwindow = lookup_widget(w, "mutate_molecule_sequence_text");
   GtkSignalFunc callback_func = GTK_SIGNAL_FUNC(mutate_sequence_molecule_menu_item_activate);


   // Get the default molecule and fill chain optionmenu with the molecules chains:
   int imol = -1; 
   for (int i=0; i<graphics_info_t::n_molecules; i++) {
      if (graphics_info_t::molecules[i].has_model()) {
	 imol = i;
	 break;
      }
   }
   if (imol >= 0) {
      graphics_info_t::mutate_sequence_imol = imol;
      GtkSignalFunc callback =
	 GTK_SIGNAL_FUNC(mutate_sequence_chain_option_menu_item_activate);
      std::string set_chain = graphics_info_t::fill_chain_option_menu(chain_option_menu, imol,
								      callback);
      graphics_info_t::mutate_sequence_chain_from_optionmenu = set_chain;
   } else {
      graphics_info_t::mutate_sequence_imol = -1; // flag for can't mutate
   }
   graphics_info_t g;
   // std::cout << "DEBUG:: filling option menu with default molecule " << imol << std::endl;
   g.fill_option_menu_with_coordinates_options(molecule_option_menu, callback_func, imol);
   return w;
}

void mutate_sequence_molecule_menu_item_activate(GtkWidget *item, 
						 GtkPositionType pos) {

   // change the chain id option menu here...
   std::cout << "DEBUG:: mutate_sequence_molecule_menu_item_activate got pos:"
	     << pos << std::endl;

   graphics_info_t::mutate_sequence_imol = pos;

   GtkWidget *chain_option_menu =
      lookup_widget(item, "mutate_molecule_chain_optionmenu");

   GtkSignalFunc callback_func =
      GTK_SIGNAL_FUNC(mutate_sequence_chain_option_menu_item_activate);
   
   std::string set_chain = graphics_info_t::fill_chain_option_menu(chain_option_menu,
								   pos, callback_func);

   graphics_info_t::mutate_sequence_chain_from_optionmenu = set_chain;
}


void mutate_sequence_chain_option_menu_item_activate (GtkWidget *item,
						      GtkPositionType pos) { 

   char *data = NULL;
   data = (char *)pos;
   std::cout << "INFO:: mutate_sequence_chain_option_menu_item_activate "
	     << " got data: " << data << std::endl;
   // this can fail when more than one sequence mutate is used at the same time:
   if (data) 
      graphics_info_t::mutate_sequence_chain_from_optionmenu = data;
}


// Don't use this function.  Use the one in graphics_info_t which you
// pass the callback function and get back a chain id.
// nvoid fill_chain_option_menu(GtkWidget *chain_option_menu, int imol) {

//   GtkSignalFunc callback_func =
//      GTK_SIGNAL_FUNC(mutate_sequence_chain_option_menu_item_activate);

   // fill_chain_option_menu_with_callback(chain_option_menu, imol, callback_func);
// }

// the generic form of the above
// void fill_chain_option_menu_with_callback(GtkWidget *chain_option_menu, int imol,
//  					  GtkSignalFunc callback_func) {

   // junk this function and use the one that returns a string.
// }




// The "Mutate" button action:
// 
void do_mutate_sequence(GtkWidget *dialog) {

#ifdef USE_PYTHON
#ifdef USE_GUILE
   short int state_lang = coot::STATE_SCM;
#else    
   short int state_lang = coot::STATE_PYTHON;
#endif
#else // python not used
#ifdef USE_GUILE
   short int state_lang = coot::STATE_SCM;
#else    
   short int state_lang = 0;
#endif
#endif   
   
   
   // decode the dialog here

   GtkWidget *entry1 = lookup_widget(dialog, "mutate_molecule_resno_1_entry");
   GtkWidget *entry2 = lookup_widget(dialog, "mutate_molecule_resno_2_entry");

   int t;
   int res1 = -9999, res2 = -99999;
   graphics_info_t g;
   
   const gchar *entry_text = gtk_entry_get_text(GTK_ENTRY(entry1));
   t = atoi(entry_text);
   if ((t > -999) && (t < 9999))
      res1 = t;
   entry_text = gtk_entry_get_text(GTK_ENTRY(entry2));
   t = atoi(entry_text);
   if ((t > -999) && (t < 9999))
      res2 = t;

   if (res2 < res1) {
      t = res1;
      res1 = res2;
      res2 = t;
   }


   // set the imol and chain_id:
   // 
   int imol = graphics_info_t::mutate_sequence_imol;
   std::string chain_id = graphics_info_t::mutate_sequence_chain_from_optionmenu;

   // Auto fit?
   GtkWidget *checkbutton = lookup_widget(dialog, "mutate_sequence_do_autofit_checkbutton"); 
   short int autofit_flag = 0;

   if (GTK_TOGGLE_BUTTON(checkbutton)->active)
      autofit_flag = 1;
      

   if (imol>= 0) {
      if (imol < graphics_info_t::n_molecules) {

	 // get the sequence:
	 GtkWidget *text = lookup_widget(dialog, "mutate_molecule_sequence_text");
	 char *txt = NULL;

	 gint start_pos = 0;
	 gint end_pos = -1;

	 txt = gtk_editable_get_chars(GTK_EDITABLE(text), start_pos, end_pos);

	 if (txt) {
	    std::string sequence(txt);
	    sequence = coot::util::plain_text_to_sequence(sequence);
	    std::cout << "we got the sequence: " << sequence << std::endl;

	    if (int(sequence.length()) == (res2 - res1 + 1)) {
	       std::vector<std::string> cmd_strings;
	       if (autofit_flag)
		  cmd_strings.push_back("mutate-and-autofit-residue-range");
	       else 
		  cmd_strings.push_back("mutate-residue-range");
	       cmd_strings.push_back(graphics_info_t::int_to_string(imol));
	       cmd_strings.push_back(single_quote(chain_id));
	       cmd_strings.push_back(graphics_info_t::int_to_string(res1));
	       cmd_strings.push_back(graphics_info_t::int_to_string(res2));
	       cmd_strings.push_back(single_quote(sequence));
	       std::string cmd = g.state_command(cmd_strings, state_lang);
	       if (state_lang == coot::STATE_SCM) {
		  safe_scheme_command(cmd);
	       }
	       update_go_to_atom_window_on_changed_mol(imol);
	    } else {
	       std::cout << "WARNING:: can't mutate.  Sequence of length: "
			 << sequence.length() << " but residue range size: "
			 << res2 - res1 + 1 << "\n";
	    } 
	 } else {
	    std::cout << "WARNING:: can't mutate.  No sequence\n";
	 } 
      } else {
	 std::cout << "WARNING:: Bad molecule number: " << imol << std::endl;
	 std::cout << "          Can't mutate." << std::endl;
      }
   } else {
      std::cout << "WARNING:: unassigned molecule number: " << imol << std::endl;
      std::cout << "          Can't mutate." << std::endl;
   }
}

GtkWidget *wrapped_fit_loop_dialog() {

   GtkWidget *w = wrapped_create_mutate_sequence_dialog();

   GtkWidget *label              = lookup_widget(w, "function_for_molecule_label");
   GtkWidget *method_frame       = lookup_widget(w, "loop_fit_method_frame");
   GtkWidget *mutate_ok_button   = lookup_widget(w, "mutate_sequence_ok_button");
   GtkWidget *fit_loop_ok_button = lookup_widget(w, "fit_loop_ok_button");
   GtkWidget *checkbutton        = lookup_widget(w, "mutate_sequence_do_autofit_checkbutton");
   
   gtk_label_set_text(GTK_LABEL(label), "\nFit loop in Molecule:\n");
   gtk_widget_hide(mutate_ok_button);
   gtk_widget_hide(checkbutton);
   gtk_widget_show(fit_loop_ok_button);

   gtk_widget_show(method_frame);

   return w;
}

// And the function called by the Fit Loop (OK) button.
// 
void fit_loop_from_widget(GtkWidget *dialog) {

#ifdef USE_PYTHON
#ifdef USE_GUILE
   short int state_lang = coot::STATE_SCM;
#else    
   short int state_lang = coot::STATE_PYTHON;
#endif
#else // python not used
#ifdef USE_GUILE
   short int state_lang = coot::STATE_SCM;
#else    
   short int state_lang = 0;
#endif
#endif

   // decode the dialog here

   GtkWidget *entry1 = lookup_widget(dialog, "mutate_molecule_resno_1_entry");
   GtkWidget *entry2 = lookup_widget(dialog, "mutate_molecule_resno_2_entry");

   int t;
   int res1 = -9999, res2 = -99999;
   graphics_info_t g;
   
   const gchar *entry_text = gtk_entry_get_text(GTK_ENTRY(entry1));
   t = atoi(entry_text);
   if ((t > -999) && (t < 9999))
      res1 = t;
   entry_text = gtk_entry_get_text(GTK_ENTRY(entry2));
   t = atoi(entry_text);
   if ((t > -999) && (t < 9999))
      res2 = t;

   if (res2 < res1) {
      t = res1;
      res1 = res2;
      res2 = t;
   }


   // set the imol and chain_id:
   // 
   int imol = graphics_info_t::mutate_sequence_imol;
   std::string chain_id = graphics_info_t::mutate_sequence_chain_from_optionmenu;

   // Auto fit?
   GtkWidget *checkbutton = lookup_widget(dialog, "mutate_sequence_do_autofit_checkbutton"); 
   short int autofit_flag = 0;

   if (GTK_TOGGLE_BUTTON(checkbutton)->active)
      autofit_flag = 1;
      

   if (imol>= 0) {
      if (imol < graphics_info_t::n_molecules) {

	 // get the sequence:
	 GtkWidget *text = lookup_widget(dialog, "mutate_molecule_sequence_text");
	 char *txt = NULL;

	 gint start_pos = 0;
	 gint end_pos = -1;

	 txt = gtk_editable_get_chars(GTK_EDITABLE(text), start_pos, end_pos);

	 if (txt) {
	    std::string sequence(txt);
	    sequence = coot::util::plain_text_to_sequence(sequence);
	    int text_widget_sequence_length = sequence.length();
	    std::cout << "INFO:: mutating to the sequence :" << sequence
		      << ":" << std::endl;

	    if (int(sequence.length()) == (res2 - res1 + 1)) {
	    } else {
	       // so set sequence to poly-ala and give us a message:
	       sequence = "";
	       for (int i=0; i<(res2 - res1 + 1); i++)
		  sequence += "A";

	       std::cout << "WARNING:: Sequence of length: "
			 << text_widget_sequence_length << " but residue range size: "
			 << res2 - res1 + 1 << ".  Using Poly-Ala\n";
	       std::string s("WARNING:: Mis-matched sequence length\nUsing Poly Ala");
	       GtkWidget *w = wrapped_nothing_bad_dialog(s);
	       gtk_widget_show(w);
	    }
	    std::vector<std::string> cmd_strings;
	    cmd_strings.push_back("fit-gap");
	    cmd_strings.push_back(graphics_info_t::int_to_string(imol));
	    cmd_strings.push_back(single_quote(chain_id));
	    cmd_strings.push_back(graphics_info_t::int_to_string(res1));
	    cmd_strings.push_back(graphics_info_t::int_to_string(res2));
	    cmd_strings.push_back(single_quote(sequence));
	    std::string cmd = g.state_command(cmd_strings, state_lang);
	    if (state_lang == coot::STATE_SCM) {
	       safe_scheme_command(cmd);
	    }
	 }
      }
   }
}



/*  ----------------------------------------------------------------------- */
/*                         Align and Mutate GUI                             */
/*  ----------------------------------------------------------------------- */
GtkWidget *wrapped_create_align_and_mutate_dialog() {

   graphics_info_t g;
   GtkWidget *w = create_align_and_mutate_dialog();
   // fill w
   GtkWidget *mol_optionmenu   = lookup_widget(w, "align_and_mutate_molecule_optionmenu");
   GtkWidget *chain_optionmenu = lookup_widget(w, "align_and_mutate_chain_optionmenu");
   GtkSignalFunc callback = GTK_SIGNAL_FUNC(align_and_mutate_molecule_menu_item_activate);
   GtkSignalFunc chain_callback = GTK_SIGNAL_FUNC(align_and_mutate_chain_option_menu_item_activate);

   int imol = graphics_info_t::align_and_mutate_imol;
   if (imol == -1 || (! g.molecules[imol].has_model())) { 
      for (int i=0; i<g.n_molecules; i++) {
	 if (g.molecules[i].has_model()) {
	    imol = i;
	    break;
	 }
      }
   }

   if (imol >= 0) {
      g.fill_option_menu_with_coordinates_options(mol_optionmenu, callback, imol);
      std::string set_chain = graphics_info_t::fill_chain_option_menu(chain_optionmenu, imol,
								      chain_callback);
      graphics_info_t::align_and_mutate_chain_from_optionmenu = set_chain;
   }
   
   return w;
}


int do_align_mutate_sequence(GtkWidget *w) {

   //
   int handled_state = 0;  // initially unhandled (return value).
   int imol = graphics_info_t::align_and_mutate_imol;
   std::string chain_id = graphics_info_t::align_and_mutate_chain_from_optionmenu;
   GtkWidget *autofit_checkbutton = lookup_widget(w, "align_and_mutate_autofit_checkbutton");

   short int do_auto_fit = 0;
   if (GTK_TOGGLE_BUTTON(autofit_checkbutton)->active)
      do_auto_fit = 1;

   graphics_info_t g;
   int imol_refinement_map = g.Imol_Refinement_Map();

   short int early_stop = 0;
   if (do_auto_fit == 1)
      if (imol_refinement_map == -1)
	 early_stop = 1;

   if (early_stop) {
      std::string s = "WARNING:: autofit requested, but \n   refinement map not set!";
      std::cout << s << "\n";
      GtkWidget *warn = wrapped_nothing_bad_dialog(s);
      gtk_widget_show(warn);
   } else { 

      handled_state = 1;
      if (imol >= 0) {
	 GtkWidget *text_w = lookup_widget(w, "align_and_mutate_sequence_text");
	 char *txt = NULL;
	 gint start_pos = 0;
	 gint end_pos = -1;
      
	 txt = gtk_editable_get_chars(GTK_EDITABLE(text_w), start_pos, end_pos);
      
	 if (txt) {
	    std::string sequence(txt);

	    if (is_valid_model_molecule(imol)) {
	       graphics_info_t g;
	       g.mutate_chain(imol, chain_id, sequence);
	       if (do_auto_fit) {
		  std::vector<std::string> s;
		  s.push_back("fit-chain");
		  s.push_back(coot::util::int_to_string(imol));
		  s.push_back(single_quote(chain_id));
#ifdef USE_GUILE
		  safe_scheme_command(languagize_command(s));
#endif // USE_GUILE		  
	       }
	       graphics_draw();
	    }
	 }
      } else {
	 std::cout << "WARNING:: inapproproate molecule number " << imol << std::endl;
      }
   }
   return handled_state;
}

void align_and_mutate(int imol, const char *chain_id, const char *fasta_maybe) {

   if (is_valid_model_molecule(imol)) {
      if (chain_id) { 
	 graphics_info_t g;
	 g.mutate_chain(imol, std::string(chain_id), std::string(fasta_maybe));
	 graphics_draw();
	 g.update_go_to_atom_window_on_changed_mol(imol);
      } else {
	 std::cout << "WARNING:: bad (NULL) chain_id - no alignment" << std::endl;
      }
   } else {
      std::cout << "WARNING:: inapproproate molecule number " << imol << std::endl;
   }
}



void align_and_mutate_molecule_menu_item_activate(GtkWidget *item, 
						  GtkPositionType pos) {

   GtkWidget *chain_optionmenu = lookup_widget(item, "align_and_mutate_chain_optionmenu");
   GtkSignalFunc chain_callback = GTK_SIGNAL_FUNC(align_and_mutate_chain_option_menu_item_activate);
   graphics_info_t::align_and_mutate_imol = pos;
   int imol = pos;
   std::string set_chain = graphics_info_t::fill_chain_option_menu(chain_optionmenu, imol,
								   chain_callback);
   
}

void align_and_mutate_chain_option_menu_item_activate (GtkWidget *item,
						       GtkPositionType pos) {

   char *data = NULL;
   data = (char *)pos;
   if (data)
      graphics_info_t::align_and_mutate_chain_from_optionmenu = data;
   std::cout << "align_and_mutate_chain_from_optionmenu is now "
	     << graphics_info_t::align_and_mutate_chain_from_optionmenu
	     << std::endl;
}


// --------------------------------------------------------------
//                 symmetry
// --------------------------------------------------------------

/* for shelx FA pdb files, there is no space group.  So allow the user
   to set it.  This can be initted with a HM symbol or a symm list for
   clipper */
void set_space_group(int imol, const char *spg) {

   // we should test that this is a clipper string here...
   // does it contain X Y and Z chars?
   // 
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_mmdb_symm(spg);
   }
}


// 
void setup_save_symmetry_coords() {

   graphics_info_t::in_save_symmetry_define = 1;
   //    GtkWidget *w = wrapped_nothing_bad_dialog(std::string("Now click on a symmetry atom"));
   std::string s = "Now click on a symmetry atom";
   graphics_info_t g;
   g.statusbar_text(s);
   pick_cursor_maybe();

}

void save_symmetry_coords_from_fileselection(GtkWidget *fileselection) {

   coot::Symm_Atom_Pick_Info_t *symm_info =
      (coot::Symm_Atom_Pick_Info_t *) gtk_object_get_user_data(GTK_OBJECT(fileselection));

   const gchar *filename = gtk_file_selection_get_filename(GTK_FILE_SELECTION(fileselection));

   if (symm_info) {
      // std::cout << "Preshift to origin:  " << symm_info->pre_shift_to_origin << std::endl;
      save_symmetry_coords(filename,
			   symm_info->imol,
			   symm_info->symm_trans.isym(),
			   symm_info->symm_trans.x(), 
			   symm_info->symm_trans.y(), 
			   symm_info->symm_trans.z(), 
			   symm_info->pre_shift_to_origin.us,			   
			   symm_info->pre_shift_to_origin.vs,
			   symm_info->pre_shift_to_origin.ws);
   } else {
      std::cout << "ERROR:: failed to get user data from save symmetry coords fileselection"
		<< std::endl;
      std::cout << "ERROR:: saving of symmetry coordinates failed" << std::endl;
   }
}

void save_symmetry_coords(const char *filename,
			  int imol, 
			  int symop_no, 
			  int shift_a, 
			  int shift_b, 
			  int shift_c,
			  int pre_shift_to_origin_na,
			  int pre_shift_to_origin_nb,
			  int pre_shift_to_origin_nc) {

   // Copy the coordinates molecule manager
   // Transform them
   // write them out

   if (imol >= 0) { 
      if (imol < graphics_info_t::n_molecules) { 
	 if (graphics_info_t::molecules[imol].has_model()) { 
	    CMMDBManager *mol2 = new CMMDBManager;
	    mol2->Copy(graphics_info_t::molecules[imol].atom_sel.mol, MMDBFCM_All);
	    
	    atom_selection_container_t asc = make_asc(mol2);
	    mat44 mat;
	    mat44 mat_origin_shift;

	    mol2->GetTMatrix(mat_origin_shift, 0,
			     -pre_shift_to_origin_na,
			     -pre_shift_to_origin_nb,
			     -pre_shift_to_origin_nc);
			     
	    mol2->GetTMatrix(mat, symop_no, shift_a, shift_b, shift_c);

	    clipper::RTop_orth to_origin_rtop(clipper::Mat33<double>(1,0,0,0,1,0,0,0,1),
					      clipper::Coord_orth(mat_origin_shift[0][3],
								  mat_origin_shift[1][3],
								  mat_origin_shift[2][3]));
	    
	    for (int i=0; i<asc.n_selected_atoms; i++) {
	       clipper::Coord_orth co;
	       clipper::Coord_orth trans_pos; 

	       clipper::Mat33<double> clipper_mat(mat[0][0], mat[0][1], mat[0][2],
						  mat[1][0], mat[1][1], mat[1][2],
						  mat[2][0], mat[2][1], mat[2][2]);
	       clipper::Coord_orth  cco(mat[0][3], mat[1][3], mat[2][3]);
	       clipper::RTop_orth rtop(clipper_mat, cco);
	       co = clipper::Coord_orth(asc.atom_selection[i]->x, 
					asc.atom_selection[i]->y, 
					asc.atom_selection[i]->z);
	       trans_pos = co.transform(to_origin_rtop);
	       trans_pos = trans_pos.transform(rtop);
	       asc.atom_selection[i]->x = trans_pos.x();
	       asc.atom_selection[i]->y = trans_pos.y();
	       asc.atom_selection[i]->z = trans_pos.z();
	    } 
	    asc.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
	    asc.mol->FinishStructEdit();
	    
	    int ierr = mol2->WritePDBASCII((char *)filename);
	    if (ierr) {
	       std::cout << "WARNING:: WritePDBASCII to " << filename << " failed." << std::endl;
	       std::string s = "WARNING:: WritePDBASCII to file ";
	       s += filename;
	       s += " failed.";
	       graphics_info_t g;
	       g.statusbar_text(s);
	    } else {
	       std::cout << "INFO:: Wrote symmetry atoms to " << filename << "." << std::endl;
	       std::string s = "INFO:: Wrote symmetry atoms to file ";
	       s += filename;
	       s += ".";
	       graphics_info_t g;
	       g.statusbar_text(s);
	    }
	 }
      }
   }
}




/*  ----------------------------------------------------------------------- */
/*                  sequence (assignment)                                   */
/*  ----------------------------------------------------------------------- */
/* section Sequence (Assignment) */

void assign_sequence(int imol_coords, int imol_map, const char *chain_id) {

   if (is_valid_model_molecule(imol_coords))
      if (is_valid_map_molecule(imol_map))
	 graphics_info_t::molecules[imol_coords].assign_sequence(graphics_info_t::molecules[imol_map].xmap_list[0], std::string(chain_id));

}


/*  ----------------------------------------------------------------------- */
/*                  base mutation                                           */
/*  ----------------------------------------------------------------------- */

void do_base_mutation(const char *type) {

   graphics_info_t g;
   int imol = g.mutate_residue_imol;

   if (is_valid_model_molecule(imol)) {

      // This is dangerous (in a pathological case). Really we should
      // save a residue spec in graphics-info-defines.cc not generate it here.
      // 
      int idx = g.mutate_residue_atom_index;
      CAtom *at = graphics_info_t::molecules[imol].atom_sel.atom_selection[idx];
      CResidue *r = at->residue;
      if (r) {
	 coot::residue_spec_t res_spec_t(r);
	 int istat = graphics_info_t::molecules[imol].mutate_base(res_spec_t, std::string(type));
	 if (istat)
	    graphics_draw();
      }
   }
}

/*  ----------------------------------------------------------------------- */
/*                  Change chain ID                                         */
/*  ----------------------------------------------------------------------- */

GtkWidget *wrapped_create_change_chain_id_dialog() {

   GtkWidget *w = create_change_chain_id_dialog();
   GtkWidget *mol_option_menu =  lookup_widget(w, "change_chain_id_molecule_optionmenu");
   GtkWidget *chain_option_menu =  lookup_widget(w, "change_chain_id_chain_optionmenu");
   GtkWidget *residue_range_no_radiobutton =  lookup_widget(w, "change_chain_residue_range_no_radiobutton");

   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(residue_range_no_radiobutton), TRUE);

   GtkSignalFunc callback_func = GTK_SIGNAL_FUNC(change_chain_ids_mol_option_menu_item_activate);

   int imol = first_coords_imol();
   if (imol >= 0) {
      graphics_info_t::change_chain_id_molecule = imol;
      GtkSignalFunc chain_callback_func =
	 GTK_SIGNAL_FUNC(change_chain_ids_chain_menu_item_activate);
      std::string set_chain = graphics_info_t::fill_chain_option_menu(chain_option_menu,
								       imol,
								      chain_callback_func);
      graphics_info_t::change_chain_id_from_chain = set_chain;
   }
   graphics_info_t g; 
   g.fill_option_menu_with_coordinates_options(mol_option_menu, callback_func, imol);
   return w;
}

void
change_chain_ids_mol_option_menu_item_activate(GtkWidget *item,
					       GtkPositionType pos) {
   graphics_info_t::change_chain_id_molecule = pos;
   int imol = pos;
   GtkWidget *chain_option_menu =  lookup_widget(item, "change_chain_id_chain_optionmenu");
   GtkSignalFunc chain_callback_func =
      GTK_SIGNAL_FUNC(change_chain_ids_chain_menu_item_activate);
   std::string set_chain = graphics_info_t::fill_chain_option_menu(chain_option_menu,
								   imol,
								   chain_callback_func);
   graphics_info_t::change_chain_id_from_chain = set_chain;
}

void
change_chain_ids_chain_menu_item_activate(GtkWidget *item,
					  GtkPositionType pos) {
   char *data = NULL;
   data = (char *)pos;
   // this can fail when more than one sequence mutate is used at the same time:
   if (data) 
      graphics_info_t::change_chain_id_from_chain = data;
}


void
change_chain_id_by_widget(GtkWidget *w) {

   GtkWidget *residue_range_yes_radiobutton =  lookup_widget(w, "change_chain_residue_range_yes_radiobutton");

   GtkWidget *residue_range_from_entry  =  lookup_widget(w, "change_chain_residues_from_entry");
   GtkWidget *residue_range_to_entry  =  lookup_widget(w, "change_chains_residues_to_entry");
   GtkWidget *change_chains_new_chain_entry  =  lookup_widget(w, "change_chains_new_chain_id");

   int imol = graphics_info_t::change_chain_id_molecule;
   short int use_res_range_flag = 0;
   int from_resno = -9999;
   int to_resno = -9999;

   if (GTK_TOGGLE_BUTTON(residue_range_yes_radiobutton)->active) { 
      use_res_range_flag = 1;
      std::pair<short int, int> p1 = int_from_entry(residue_range_from_entry);
      std::pair<short int, int> p2 = int_from_entry(residue_range_to_entry);
      if (p1.first)
	 from_resno = p1.second;
      if (p2.first)
	 to_resno = p2.second;
   }

   const gchar *txt = gtk_entry_get_text(GTK_ENTRY(change_chains_new_chain_entry));

   if (txt) { 
   
      if (is_valid_model_molecule(imol)) {
	 std::string to_chain_id(txt);
	 std::string from_chain_id = graphics_info_t::change_chain_id_from_chain;
	 std::pair<int, std::string> r = 
	    graphics_info_t::molecules[imol].change_chain_id(from_chain_id,
							     to_chain_id,
							     use_res_range_flag,
							     from_resno,
							     to_resno);
	 if (r.first == 1) { // it went OK
	    update_go_to_atom_window_on_changed_mol(imol);
	    graphics_draw();
	 } else {
	    GtkWidget *ws = wrapped_nothing_bad_dialog(r.second);
	    gtk_widget_show(ws);
	 }
      }
   } else {
      std::cout << "ERROR: Couldn't get txt in change_chain_id_by_widget\n";
   }
}


/*  ----------------------------------------------------------------------- */
/*               Nomenclature Errors                                        */
/*  ----------------------------------------------------------------------- */
/* return the number of resides altered. */
int fix_nomenclature_errors(int imol) {

   int ifixed = 0;
   if (is_valid_model_molecule(imol)) {
      std::vector<CResidue *> vr =
	 graphics_info_t::molecules[imol].fix_nomenclature_errors();
      ifixed = vr.size();
      graphics_info_t g;
      g.update_geometry_graphs(graphics_info_t::molecules[imol].atom_sel,
			       imol);
      graphics_draw();
   }
   // update geometry graphs (not least rotamer graph).
   // but we have no intermediate atoms...
   
   return ifixed; 

} 

/*  ----------------------------------------------------------------------- */
/*                  cis <-> trans conversion                                */
/*  ----------------------------------------------------------------------- */
void do_cis_trans_conversion_setup(int istate) {

   if (istate == 1) { 
      graphics_info_t::in_cis_trans_convert_define = 1;
      pick_cursor_maybe(); // depends on ctrl key for rotate
   } else {
      graphics_info_t::in_cis_trans_convert_define = 0;
      normal_cursor(); // depends on ctrl key for rotate
   }
} 

/*  ----------------------------------------------------------------------- */
/*                  180 degree flip                                         */
/*  ----------------------------------------------------------------------- */
/* rotate 180 degrees round the last chi angle */
void do_180_degree_side_chain_flip(int imol, const char* chain_id, int resno, 
			const char *inscode, const char *altconf) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      int istatus =
	 g.molecules[imol].do_180_degree_side_chain_flip(std::string(chain_id),
							 resno,
							 std::string(inscode),
							 std::string(altconf),
							 g.Geom_p());
      std::string s;
      if (istatus == 0) {
	 s = "Problem flipping chi angle on residue ";
	 s += chain_id;
	 s += graphics_info_t::int_to_string(resno);
	 s += ". Not done.";
      } else {
	 s = "Chi angle on residue ";
	 s += chain_id;
	 s += graphics_info_t::int_to_string(resno);
	 s += " successfully flipped.";
	 graphics_draw();
      }
      g.statusbar_text(s);
   }
}

// graphics click stuff
void setup_180_degree_flip(short int state) {

   graphics_info_t g;
   graphics_info_t::in_180_degree_flip_define = state;
   if (state) {
      g.in_180_degree_flip_define = 1;
      std::cout << "Click on a residue that you want to flip" << std::endl;
      g.pick_cursor_maybe();
      g.statusbar_text("Click on an atom in the residue that you want to flip");
      g.pick_pending_flag = 1;
   } else {
      g.normal_cursor();
      g.pick_pending_flag = 0;
   }
}

/*  ----------------------------------------------------------------------- */
/*                  reverse direction                                       */
/*  ----------------------------------------------------------------------- */
void
setup_reverse_direction(short int istate) {

   graphics_info_t::in_reverse_direction_define = istate;
   graphics_info_t g;
   if (istate == 1) {
      g.pick_cursor_maybe();
      g.statusbar_text("Click on an atom in the fragment that you want to reverse");
      g.pick_pending_flag = 1;
   } else {
      g.normal_cursor();
   }

}



/*  ----------------------------------------------------------------------- */
/*                  De-chainsaw                                             */
/*  ----------------------------------------------------------------------- */
// Fill amino acid residues
void fill_partial_residues(int imol) {

   if (is_valid_model_molecule(imol)) { 
      graphics_info_t g;
      int imol_map = g.Imol_Refinement_Map();
      std::cout << "DEBUG:: calling mci fill_partial_residues " << std::endl;
      coot::util::missing_atom_info m_i_info =
	 g.molecules[imol].fill_partial_residues(g.Geom_p(), imol_map);
      std::cout << "DEBUG:: done    mci fill_partial_residues " << std::endl;
      graphics_draw();

      if (imol_map > -1) { 
	 int refinement_replacement_state = refinement_immediate_replacement_state();
	 set_refinement_immediate_replacement(1);
      	 for (unsigned int i=0; i<m_i_info.residues_with_missing_atoms.size(); i++) {
      	    int resno =  m_i_info.residues_with_missing_atoms[i]->GetSeqNum();
      	    std::string chain_id = m_i_info.residues_with_missing_atoms[i]->GetChainID();
      	    std::string residue_type = m_i_info.residues_with_missing_atoms[i]->GetResName();
      	    std::string inscode = m_i_info.residues_with_missing_atoms[i]->GetInsCode();
      	    std::string altconf("");
	    short int is_water = 0;
	    // hmmm backups are being done....
      	    g.refine_residue_range(imol, chain_id, chain_id, resno, resno, altconf, is_water);
	    accept_regularizement();
      	 }
	 set_refinement_immediate_replacement(refinement_replacement_state);
      } else {
	 g.show_select_map_dialog();
      } 
   }
}

void
copy_chain(int imol, const char *from_chain, const char *to_chain) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].copy_chain(std::string(from_chain),
						  std::string(to_chain));
      graphics_draw();
   }
} 


void copy_from_ncs_master_to_others(int imol, const char *chain_id) {

   if (is_valid_model_molecule(imol)) {
      std::string c(chain_id);
      int ncopied =
	 graphics_info_t::molecules[imol].copy_from_ncs_master_to_others(c);
      graphics_draw();
   }
}



void copy_residue_range_from_ncs_master_to_others(int imol, 
			const char *master_chain_id,
			int res_range_start,
			int res_range_end) {

   if (is_valid_model_molecule(imol)) {
      std::string c(master_chain_id);
      int ncopied =
	 graphics_info_t::molecules[imol].copy_residue_range_from_ncs_master_to_others(c,
										       res_range_start,
										       res_range_end);
      graphics_draw();
   }
}


/*  ----------------------------------------------------------------------- */
/*                  Helices                                                 */
/*  ----------------------------------------------------------------------- */
/* Add a helix somewhere close to this point in the map, try to fit
   the orientation. Add to a molecule called "Helices", create it if needed.*/
int place_helix_here() {

   graphics_info_t g;
   clipper::Coord_orth pt(g.RotationCentre_x(),
			  g.RotationCentre_y(),
			  g.RotationCentre_z());
   int imol_map = g.Imol_Refinement_Map();
   if (imol_map != -1) {
      // perhaps we should add to the constructor the user-setable
      // g.helix_placement_cb_density_descrimination_ratio
      // defaults to 1.0.
      int imol = -1;
      coot::helix_placement p(graphics_info_t::molecules[imol_map].xmap_list[0]);
      coot::helix_placement_info_t n = p.place_alpha_helix_near_kc_version(pt, 20);
       if (! n.success) {
 	 n = p.place_alpha_helix_near_kc_version(pt, 9);
       }

       if (n.success) {
	  float bf = graphics_info_t::default_new_atoms_b_factor;
	 atom_selection_container_t asc = make_asc(n.mol[0].pcmmdbmanager(bf));
	 g.expand_molecule_space_maybe();
	 imol = g.n_molecules;
	 graphics_info_t::molecules[imol].install_model(asc, "Helix", 1);
	 g.n_molecules++;

	 if (n.mol.size() > 1) { 
	    atom_selection_container_t asc2 = make_asc(n.mol[1].pcmmdbmanager(bf));
	    g.expand_molecule_space_maybe();
	    imol = g.n_molecules;
	    graphics_info_t::molecules[imol].install_model(asc2, "Reverse Helix", 1);
	    g.n_molecules++;
	 }
	 
	 if (g.go_to_atom_window) {
	    g.set_go_to_atom_molecule(imol);
	    g.update_go_to_atom_window_on_new_mol();
	 } else {
	    g.set_go_to_atom_molecule(imol);
	 }
	 g.statusbar_text("Helix added");
      } else {
	 std::cout << "Helix addition failure: message: " << n.failure_message << "\n";
	 g.statusbar_text(n.failure_message);
      }
      std::vector<std::string> command_strings;
      command_strings.push_back("set-rotation-centre");
      command_strings.push_back(coot::util::float_to_string(g.RotationCentre_x()));
      command_strings.push_back(coot::util::float_to_string(g.RotationCentre_y()));
      command_strings.push_back(coot::util::float_to_string(g.RotationCentre_z()));
      add_to_history(command_strings);
      
      command_strings.resize(0);
      command_strings.push_back("place-helix-here");
      add_to_history(command_strings);
      graphics_draw();
      return imol;
   } else {
      std::cout << " You need to set the map to fit against\n";
      g.statusbar_text("You need to set the map to fit against");
      g.show_select_map_dialog();
      return -1;
   }
} 


/*  ----------------------------------------------------------------------- */
/*                  Place a strand                                          */
/*  ----------------------------------------------------------------------- */
int place_strand_here(int n_residues, int n_sample_strands) {

   int imol = -1; // failure status 
   graphics_info_t g;
   clipper::Coord_orth pt(g.RotationCentre_x(),
			  g.RotationCentre_y(),
			  g.RotationCentre_z());
   int imol_map = g.Imol_Refinement_Map();
   if (imol_map != -1) {

      coot::helix_placement p(graphics_info_t::molecules[imol_map].xmap_list[0]);
      coot::helix_placement_info_t si = p.place_strand(pt, n_residues, n_sample_strands);
      if (si.success) {
	 // nice to refine the fragment here, but the interface
	 // doesn't work that way, so put the refinement after the
	 // molecule has been accepted.
	 float bf = graphics_info_t::default_new_atoms_b_factor;
	 atom_selection_container_t asc = make_asc(si.mol[0].pcmmdbmanager(bf));
	 g.expand_molecule_space_maybe();
	 imol = g.n_molecules;
	 graphics_info_t::molecules[imol].install_model(asc, "Strand", 1);
	 g.statusbar_text("Strand added");
	 g.n_molecules++;

	 // Now refine.
	 coot::minimol::zone_info_t zi = si.mol[0].zone_info();
	 if (zi.is_simple_zone) {
	    graphics_info_t g;
	    int save_rirf = g.refinement_immediate_replacement_flag;
	    coot::pseudo_restraint_bond_type save_pseudos = g.pseudo_bonds_type;
	    g.pseudo_bonds_type = coot::STRAND_PSEUDO_BONDS;
	    g.refinement_immediate_replacement_flag = 1;
	    g.refine_residue_range(imol, zi.chain_id, zi.chain_id, zi.resno_1, zi.resno_2,
				   "", 0);
	    g.refinement_immediate_replacement_flag = save_rirf;
	    accept_regularizement();
	    g.pseudo_bonds_type = save_pseudos;
	 } 
      } else {
	 std::cout << "Strand addition failure: message: " << si.failure_message << "\n";
	 g.statusbar_text(si.failure_message);
      }
      if (g.go_to_atom_window) {
	 g.set_go_to_atom_molecule(imol);
	 g.update_go_to_atom_window_on_new_mol();
      }
      
      std::vector<std::string> command_strings;
      command_strings.push_back("set-rotation-centre");
      command_strings.push_back(coot::util::float_to_string(g.RotationCentre_x()));
      command_strings.push_back(coot::util::float_to_string(g.RotationCentre_y()));
      command_strings.push_back(coot::util::float_to_string(g.RotationCentre_z()));
      add_to_history(command_strings);
      
      command_strings.resize(0);
      command_strings.push_back("place-strand-here");
      command_strings.push_back(coot::util::int_to_string(n_residues));
      command_strings.push_back(coot::util::int_to_string(n_sample_strands));
      add_to_history(command_strings);
      graphics_draw();
      return imol;
   } else {
      std::cout << " You need to set the map to fit against\n";
      g.statusbar_text("You need to set the map to fit against");
      g.show_select_map_dialog();
      return -1;
   }
   

} 



/*  ----------------------------------------------------------------------- */
/*                  FFfearing                                               */
/*  ----------------------------------------------------------------------- */


// return the new model number
// 
int fffear_search(int imol_model, int imol_map) {

   float angular_resolution = graphics_info_t::fffear_angular_resolution;
   int imol_new = -1;
   if (is_valid_model_molecule(imol_model))
      if (is_valid_map_molecule(imol_map)) { 
	 coot::util::fffear_search f(graphics_info_t::molecules[imol_model].atom_sel.mol,
				     graphics_info_t::molecules[imol_model].atom_sel.SelectionHandle,
				     graphics_info_t::molecules[imol_map].xmap_list[0],
				     angular_resolution);

	 imol_new = graphics_info_t::n_molecules;
	 std::string name("FFFear search results");
	 graphics_info_t::molecules[imol_new].new_map(f.get_results_map(), name);
	 graphics_info_t::n_molecules++;

	 std::vector<std::pair<float, clipper::RTop_orth> > p = f.scored_orientations();
	 if (p.size() > 0) {
	    // install new molecule(s) that has been rotated.
	 }

      }

   return imol_new;

}


void set_fffear_angular_resolution(float f) {

   graphics_info_t::fffear_angular_resolution = f; 
}

float fffear_angular_resolution() {

   return graphics_info_t::fffear_angular_resolution; 
}

/*  \brief create a molecule of idea nucleotides 

use the given sequence (single letter code)

RNA_or_DNA is either "RNA" or "DNA"

form is either "A" or "B"

@return the new molecule number or -1 if a problem */
int ideal_nucleic_acid(const char *RNA_or_DNA, const char *form,
		       short int single_stranded_flag,
		       const char *sequence) {

   int istat = -1; 
   short int do_rna_flag = -1;
   short int form_flag = -1;

   float here_x = graphics_info_t::RotationCentre_x();
   float here_y = graphics_info_t::RotationCentre_y();
   float here_z = graphics_info_t::RotationCentre_z();

   std::string RNA_or_DNA_str(RNA_or_DNA);
   std::string form_str(form);

   if (RNA_or_DNA_str == "RNA")
      do_rna_flag = 1;
   if (RNA_or_DNA_str == "DNA")
      do_rna_flag = 0;

   if (form_str == "A")
      form_flag = 1;
   
   if (form_str == "B")
      form_flag = 1;

   if (! (form_flag > 0)) {
      std::cout << "Problem in nucleic acid form, use only either \"A\" or \"B\"."
		<< std::endl;
   } else {
      if (! (do_rna_flag >= 0)) {
	 std::cout << "Problem in nucleic acid type, use only either \"RNA\" or \"DNA\"."
		   << "You said: \"" << RNA_or_DNA << "\"" << std::endl;
      } else {
	 // proceed, input is good

	 std::string down_sequence(coot::util::downcase(sequence));
	 if (graphics_info_t::standard_residues_asc.read_success) {
	    coot::ideal_rna ir(RNA_or_DNA_str, form_str, single_stranded_flag,
			       down_sequence,
			       graphics_info_t::standard_residues_asc.mol);
	    CMMDBManager *mol = ir.make_molecule();

	    if (mol) { 
	       int imol = graphics_info_t::n_molecules;
	       istat = imol;
	       std::string label = "Ideal " + form_str;
	       label += " form ";
	       label += RNA_or_DNA_str;
	       atom_selection_container_t asc = make_asc(mol);
	       graphics_info_t::molecules[imol].install_model(asc, label, 1);
	       graphics_info_t::molecules[imol].translate_by(here_x, here_y, here_z);
	       graphics_info_t::n_molecules++;
	       graphics_draw();
	       if (graphics_info_t::go_to_atom_window) {
		  graphics_info_t g;
		  g.update_go_to_atom_window_on_new_mol();
		  g.update_go_to_atom_window_on_changed_mol(imol);
	       }
	    }
	 } else {
	    std::string s("WARNING:: Can't proceed with Idea RNA - no standard residues!");
	    std::cout << s << std::endl;
	    graphics_info_t g;
	    g.statusbar_text(s);
	 } 
      }
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("ideal-nucleic-acid");
   command_strings.push_back(single_quote(RNA_or_DNA_str));
   command_strings.push_back(single_quote(form_str));
   command_strings.push_back(coot::util::int_to_string(single_stranded_flag));
   command_strings.push_back(single_quote(sequence));
   add_to_history(command_strings);

   return istat;
}


GtkWidget *wrapped_nucleotide_builder_dialog() {

   GtkWidget *w = create_nucleotide_builder_dialog(); 
   return w;
} 

void ideal_nucleic_acid_by_widget(GtkWidget *builder_dialog) {

   std::string type = "RNA";
   std::string form = "A";
   short int single_stranded_flag = 0;
   GtkWidget *entry = lookup_widget(builder_dialog, "nucleotide_sequence");
   GtkWidget *type_optionmenu = lookup_widget(builder_dialog,
					      "nucleotide_builder_type_optionmenu");
   GtkWidget *form_optionmenu = lookup_widget(builder_dialog,
					      "nucleotide_builder_form_optionmenu");
   GtkWidget *strand_optionmenu = lookup_widget(builder_dialog,
						"nucleotide_builder_strand_optionmenu");


   GtkWidget *menu;
   GtkWidget *active_item;
   int active_index;

   menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(type_optionmenu));
   active_item = gtk_menu_get_active(GTK_MENU(menu));
   active_index = g_list_index(GTK_MENU_SHELL(menu)->children, active_item);
   std::cout << "DEBUG:: active_index for type: " << active_index << std::endl;
   if (active_index == 1)
      type = "DNA";

   menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(form_optionmenu));
   active_item = gtk_menu_get_active(GTK_MENU(menu));
   active_index = g_list_index(GTK_MENU_SHELL(menu)->children, active_item);
   std::cout << "DEBUG:: active_index for form: " << active_index << std::endl;
   if (active_index == 1)
      form = "B";

   menu = gtk_option_menu_get_menu(GTK_OPTION_MENU(strand_optionmenu));
   active_item = gtk_menu_get_active(GTK_MENU(menu));
   active_index = g_list_index(GTK_MENU_SHELL(menu)->children, active_item);
   std::cout << "DEBUG:: active_index for strand: " << active_index << std::endl;
   if (active_index == 1)
      single_stranded_flag = 1;

   
   const char *txt = gtk_entry_get_text(GTK_ENTRY(entry));
   if (txt) {
      ideal_nucleic_acid(type.c_str(), form.c_str(), single_stranded_flag, txt);
   }
}

/*  ----------------------------------------------------------------------- */
/*                  rigid body refinement                                   */
/*  ----------------------------------------------------------------------- */

void do_rigid_body_refine(short int state){

   graphics_info_t g;
   
   g.set_in_rigid_body_refine(state);
   if (state) { 
      g.pick_cursor_maybe();
      g.pick_pending_flag = 1;
      std::cout << "click on 2 atoms to define a range of residue "
		<< "to rigid body refine: " << std::endl;
   } else {
      g.normal_cursor();
   }
}

void execute_rigid_body_refine(short int auto_range_flag){
   graphics_info_t g;
   g.execute_rigid_body_refine(auto_range_flag);
}

void rigid_body_refine_zone(int resno_start, int resno_end, 
			    const char *chain_id, int imol) {
   graphics_info_t g;

   // need to set graphics_info's residue_range_atom_index_1,
   // residue_range_atom_index_2, imol_rigid_body_refine

   if (imol < g.n_molecules) {
      if (g.molecules[imol].has_model()) { 
	 g.imol_rigid_body_refine = imol;

	 g.set_residue_range_refine_atoms(resno_start, resno_end,
					  std::string(chain_id), imol);
	 g.execute_rigid_body_refine(0);
      }
   }
}

void fill_option_menu_with_refine_options(GtkWidget *option_menu) { 

   graphics_info_t g;

   g.fill_option_menu_with_map_options(option_menu, 
				       GTK_SIGNAL_FUNC(graphics_info_t::refinement_map_select));
}

void
set_rigid_body_fit_acceptable_fit_fraction(float f) {
   if (f >= 0.0 && f<= 1.0) { 
      graphics_info_t::rigid_body_fit_acceptable_fit_fraction = f;
   } else {
      std::cout << "ignoring set_rigid_body_fit_acceptable_fit_fraction"
		<< " of " << f << std::endl;
   } 
} 

/*  ----------------------------------------------------------------------- */
//                                 refinement
/*  ----------------------------------------------------------------------- */

void set_matrix(float f) {

   graphics_info_t::geometry_vs_map_weight = f;
}

float matrix_state() {
   return graphics_info_t::geometry_vs_map_weight;
}

void set_refine_auto_range_step(int i) { 
   graphics_info_t::refine_auto_range_step = i;
} 

void set_refine_max_residues(int n) { 
   graphics_info_t::refine_regularize_max_residues = n;
}

// (refine-zone-atom-index-define 0 688 688)
// 
void refine_zone_atom_index_define(int imol, int ind1, int ind2) {
   
   graphics_info_t g;

   if (g.n_molecules > imol) {
      if (g.molecules[imol].has_model()) {
	 if (g.molecules[imol].atom_sel.n_selected_atoms > ind1 &&
	     g.molecules[imol].atom_sel.n_selected_atoms > ind2) {
	    g.refine(imol, 0, ind1, ind2);
	 } else {
	    std::cout << "WARNING: atom index error in "
		      << "refine_zone_atom_index_define\n";
	 }
      } else {
	 std::cout << "WARNING: no model for molecule " << imol << " in "
		   << "refine_zone_atom_index_define\n";
      }
   } else {
      std::cout << "WARNING: no molecule " << imol << " in "
		<< "refine_zone_atom_index_define\n";
   }
}

void refine_zone(int imol, const char *chain_id,
		 int resno1,
		 int resno2,
		 const char *altconf) {

   if (imol >= 0) {
      if (is_valid_model_molecule(imol)) {
	 graphics_info_t g;
// 	 int index1 = atom_index(imol, chain_id, resno1, " CA ");
// 	 int index2 = atom_index(imol, chain_id, resno2, " CA ");
	 // the "" is the insertion code (not passed to this function (yet)
	 int index1 = graphics_info_t::molecules[imol].atom_index_first_atom_in_residue(chain_id, resno1, ""); 
	 int index2 = graphics_info_t::molecules[imol].atom_index_first_atom_in_residue(chain_id, resno2, "");
	 short int auto_range = 0;
	 if (index1 >= 0) {
	    if (index2 >= 0) { 
	       g.refine(imol, auto_range, index1, index2);
	    } else {
	       std::cout << "WARNING:: refine_zone: Can't get index for resno2: "
			 << resno2 << std::endl;
	    } 
	 } else {
	    std::cout << "WARNING:: refine_zone: Can't get index for resno1: "
		      << resno1 << std::endl;
	 }
      } else {
	 std::cout << "Not a valid model molecule" << std::endl;
      }
   }
}


void refine_auto_range(int imol, const char *chain_id, int resno1, const char *altconf) {

   if (imol >= 0) {
      if (imol < graphics_info_t::n_molecules) {
	 graphics_info_t g;
	 int index1 = atom_index(imol, chain_id, resno1, " CA ");
	 short int auto_range = 1;
	 if (index1 >= 0) { 
	    g.refine(imol, auto_range, index1, index1);
	 } else {
	    std::cout << "WARNING:: refine_auto_range: Can't get index for resno1: "
		      << resno1 << std::endl;
	 }
      }
   }
}

// This does not control if the atoms are accepted immediately, just
// whether the Accept Refinemnt gui is shown.
// 
void set_refinement_immediate_replacement(int istate) {
   graphics_info_t::refinement_immediate_replacement_flag = istate;
}

int  refinement_immediate_replacement_state() {
   return graphics_info_t::refinement_immediate_replacement_flag; 
} 


int imol_refinement_map() {

   graphics_info_t g;
   return g.Imol_Refinement_Map();

}

int set_imol_refinement_map(int imol) {

   int r = -1; 
   if (is_valid_map_molecule(imol)) { 
      graphics_info_t g;
      r = g.set_imol_refinement_map(imol);
   }
   return r;
}

/*  ----------------------------------------------------------------------- */
/*                  regularize/refine                                       */
/*  ----------------------------------------------------------------------- */

void do_regularize(short int state) { 

   //
   graphics_info_t g; 

   g.set_in_range_define_for_regularize(state);  // TRUE or FALSE
   if (state) { 
      g.untoggle_model_fit_refine_buttons_except("model_refine_dialog_regularize_zone_togglebutton");
      // and kill the delete dialog if it is there
      if (g.delete_item_widget) { 
	 gtk_widget_destroy(g.delete_item_widget);
	 g.delete_item_widget = NULL;
	 // hopefully superfluous:
	 g.delete_item_atom = 0;
	 g.delete_item_residue = 0;
	 g.delete_item_residue_hydrogens = 0;
      }
      std::cout << "click on 2 atoms (in the same molecule)" << std::endl; 
      g.pick_cursor_maybe();
      g.pick_pending_flag = 1;
   } else { 
      g.normal_cursor();
   }
}

void do_refine(short int state) { 

   //
   graphics_info_t g; 

   g.set_in_range_define_for_refine(state);  // TRUE or Not...

   // std::cout << "DEBUG:: in do_refine" << std::endl;
   
   if (state) { 
      g.untoggle_model_fit_refine_buttons_except("model_refine_dialog_refine_togglebutton");
      // and kill the delete dialog if it is there
      if (g.delete_item_widget) { 
	 gtk_widget_destroy(g.delete_item_widget);
	 g.delete_item_widget = NULL;
	 // hopefully superfluous:
	 g.delete_item_atom = 0;
	 g.delete_item_residue = 0;
	 g.delete_item_residue_hydrogens = 0;
      }
      
      int imol_map = g.Imol_Refinement_Map();
      // std::cout << "DEBUG:: in do_refine, imol_map: " << imol_map << std::endl;
      if (imol_map >= 0) {
	 if (g.molecules[imol_map].has_map()) { 
	    std::cout << "click on 2 atoms (in the same molecule)" << std::endl; 
	    g.pick_cursor_maybe();
	    g.pick_pending_flag = 1;
	    std::string s = "Pick 2 atoms or Autozone (pick 1 atom the press the A key)...";
	    g.statusbar_text(s);
	 } else {
	    g.show_select_map_dialog();
	    g.in_range_define_for_refine = 0;
	    g.model_fit_refine_unactive_togglebutton("model_refine_dialog_refine_togglebutton");
	 }
      } else {
	 // map chooser dialog
	 g.show_select_map_dialog();
	 g.in_range_define_for_refine = 0;
	 g.model_fit_refine_unactive_togglebutton("model_refine_dialog_refine_togglebutton");
      }
   } else { 
      g.normal_cursor();
      g.in_range_define_for_refine = 0;
      // g.pick_pending_flag = 0;
   }

}

void set_residue_selection_flash_frames_number(int i) {

   graphics_info_t::residue_selection_flash_frames_number = i;
}

void set_refinement_refine_per_frame(int i) {

   graphics_info_t::dragged_refinement_refine_per_frame_flag = i; 
}

int  refinement_refine_per_frame_state() {

   return graphics_info_t::dragged_refinement_refine_per_frame_flag;
}

void set_secondary_structure_restraints_type(int itype) {

#ifdef HAVE_GSL   

   if (itype == 0)
      graphics_info_t::pseudo_bonds_type = coot::NO_PSEUDO_BONDS;
   if (itype == 1)
      graphics_info_t::pseudo_bonds_type = coot::HELIX_PSEUDO_BONDS;
   if (itype == 2)
      graphics_info_t::pseudo_bonds_type = coot::STRAND_PSEUDO_BONDS;
#endif // HAVE_GSL   
} 

/*! \brief return the secondary structure restraints type */
int secondary_structure_restraints_type() {

   // cast a pseudo_restraint_bond_type to an int
#ifdef HAVE_GSL   
   return graphics_info_t::pseudo_bonds_type;
#else
   return 0;
#endif // HAVE_GSL   
} 


void accept_regularizement() {

   graphics_info_t g;
   g.accept_moving_atoms();	// does a g.clear_up_moving_atoms();
   g.clear_moving_atoms_object();
}


/* \brief Experimental interface for Ribosome People. 

Ribosome People have many chains in their pdb file, they prefer segids
to chainids (chainids are only 1 character).  But coot uses the
concept of chain ids and not seg-ids.  mmdb allow us to use more than
one char in the chainid, so after we read in a pdb, let's replace the
chain ids with the segids. Will that help? */
int exchange_chain_ids_for_seg_ids(int imol) {

   int istat = 0;

   if (is_valid_model_molecule(imol)) {
      istat = graphics_info_t::molecules[imol].exchange_chain_ids_for_seg_ids();
      graphics_draw();
      graphics_info_t g;
      g.update_go_to_atom_window_on_changed_mol(imol);
   }
   return istat;
}

/*  ----------------------------------------------------------------------- */
/*             New Molecule by Various Selection                            */
/*  ----------------------------------------------------------------------- */
/*! \brief create a new molecule that consists of only the residue of
  type residue_type in molecule number imol

@return the noew molecule number, -1 means an error. */
int new_molecule_by_residue_type_selection(int imol_orig, const char *residue_type) {

   int imol = -1;

   if (is_valid_model_molecule(imol_orig)) {

      PCAtom *atom_selection = 0;
      int nSelAtoms = 0;
      imol = graphics_info_t::n_molecules;
      CMMDBManager *mol_orig = graphics_info_t::molecules[imol_orig].atom_sel.mol;
      int SelectionHandle = mol_orig->NewSelection();
      mol_orig->SelectAtoms(SelectionHandle, 0, "*",
			    ANY_RES, "*",
			    ANY_RES, "*",
			    (char *) residue_type,
			    "*", "*", "*");
      CMMDBManager *mol =
	 coot::util::create_mmdbmanager_from_atom_selection(mol_orig, SelectionHandle);

      if (mol) {
	 std::string name = "residue type ";
	 name += residue_type;
	 name += " from ";
	 name += graphics_info_t::molecules[imol_orig].name_for_display_manager();
	 atom_selection_container_t asc = make_asc(mol);
	 graphics_info_t::molecules[imol].install_model(asc, name, 1);
	 graphics_info_t::n_molecules++;
      } else {
	 std::cout << "in new_molecule_by_residue_type_selection "
		   << "Something bad happened - null molecule" << std::endl;
      } 
      mol_orig->DeleteSelection(SelectionHandle);
      graphics_draw();
   } else {
      std::cout << "Molecule number " << imol_orig << " is not a valid "
		<< "model molecule" << std::endl;
   } 
   return imol;
}

int new_molecule_by_atom_selection(int imol_orig, const char* atom_selection_str) {

   int imol = -1;
   if (is_valid_model_molecule(imol_orig)) {
      PCAtom *atom_selection = 0;
      int nSelAtoms = 0;
      imol = graphics_info_t::n_molecules;
      CMMDBManager *mol_orig = graphics_info_t::molecules[imol_orig].atom_sel.mol;
      int SelectionHandle = mol_orig->NewSelection();
      mol_orig->Select(SelectionHandle, STYPE_ATOM,
		       (char *) atom_selection_str, // sigh... Why do I have to do this?
		                                    // mmdb_selmngr.h says this arg should
		                                    // be const pstr CID.  Hmmm...
		       SKEY_OR);
      CMMDBManager *mol =
	 coot::util::create_mmdbmanager_from_atom_selection(mol_orig, SelectionHandle);

      { // debug code 
	 int imod = 1;
	 
	 CModel *model_p = mol->GetModel(imod);
	 CChain *chain_p;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<nchains; ichain++) {
	    chain_p = model_p->GetChain(ichain);
	    std::cout << "In new_molecule_by_atom_selection has contains chain id "
		      << chain_p->GetChainID() << std::endl;
	 }
      }

      if (mol) {
	 std::string name = "atom selection from ";
	 name += graphics_info_t::molecules[imol_orig].name_for_display_manager();
	 atom_selection_container_t asc = make_asc(mol);
	 graphics_info_t::molecules[imol].install_model(asc, name, 1);
	 graphics_info_t::n_molecules++;
      } else {
	 std::cout << "in new_molecule_by_atom_selection "
		   << "Something bad happened - null molecule" << std::endl;
      } 
      mol_orig->DeleteSelection(SelectionHandle);
      graphics_draw();
   } else {
      std::cout << "Molecule number " << imol_orig << " is not a valid "
		<< "model molecule" << std::endl;
   }
   return imol;
} 


// ---------------------------------------------------------------------
// b-factor
// ---------------------------------------------------------------------


void set_default_temperature_factor_for_new_atoms(float new_b) {

   graphics_info_t::default_new_atoms_b_factor = new_b;
} 

float default_new_atoms_b_factor() {
   return graphics_info_t::default_new_atoms_b_factor;
} 
