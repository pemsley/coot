/* src/c-interface-build.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007, 2008 The University of York
 * Author: Paul Emsley
 * Copyright 2007 by Paul Emsley
 * Copyright 2007 by Bernhard Lohkamp
 * Copyright 2008 by Kevin Cowtan
 * Copyright 2007, 2008, 2009, 2010, 2011 The University of Oxford
 * Copyright 2016 by Medical Research Council
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

#include "geometry/residue-and-atom-specs.hh"
#include "gtk/gtk.h"
#ifdef USE_PYTHON
#include <Python.h>  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#include "python-3-interface.hh"
#endif

#include "compat/coot-sysdep.h"

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#define HAVE_CIF  // will become unnessary at some stage.

#include <sys/types.h> // for stating
#include <sys/stat.h>
#include <string.h> // strncmp
#if !defined _MSC_VER
#include <unistd.h>
#else
#define S_ISDIR(m) (((m) & S_IFMT) == S_IFDIR)
#define S_ISREG(m) (((m) & S_IFMT) == S_IFREG)
#include <windows.h>
#endif

#include <mmdb2/mmdb_manager.h>

#include "coords/mmdb-extras.hh"
#include "coords/mmdb.hh"
#include "coords/Cartesian.hh"
#include "coords/Bond_lines.hh"

// just delete this header?
#include "globjects.h" //includes gtk/gtk.h

#include "coords/mmdb-crystal.hh"


#include "graphics-info.h"

#include "coot-utils/coot-coord-utils.hh"
#include "utils/coot-fasta.hh"

#include "skeleton/BuildCas.h"
#include "ligand/helix-placement.hh"
#include "ligand/fast-ss-search.hh"

#include "utils/coot-utils.hh"  // for is_member_p
#include "coot-utils/coot-map-heavy.hh"  // for fffear

#include "guile-fixups.h"


#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "cc-interface.hh"
#include "cc-interface-scripting.hh"

#include "ligand/ligand.hh" // for rigid body fit by atom selection.

#include "cmtz-interface.hh" // for valid columns mtz_column_types_info_t
#include "c-interface-mmdb.hh"
#include "c-interface-scm.hh"
#include "c-interface-python.hh"

#ifdef USE_DUNBRACK_ROTAMERS
#include "ligand/dunbrack.hh"
#else 
#include "ligand/richardson-rotamer.hh"
#endif

#include "ligand/backrub-rotamer.hh"
#include "rotamer-search-modes.hh"

#include "protein_db/protein_db_utils.h"
#include "protein_db-interface.hh"

#include "cootilus/cootilus-build.h"

#include "c-interface-refine.hh"



void delete_residue_sidechain(int imol, const char *chain_id, int resno, const char *ins_code,
			      short int do_delete_dialog) {

   std::string inscode(ins_code);
   graphics_info_t g;

   if (is_valid_model_molecule(imol)) { 
      mmdb::Residue *residue_p =
	 graphics_info_t::molecules[imol].get_residue(chain_id, resno, ins_code);
      if (residue_p) {
	 coot::residue_spec_t spec(residue_p);
	 g.delete_residue_from_geometry_graphs(imol, spec);
      }
      short int istat =
	 g.molecules[imol].delete_residue_sidechain(std::string(chain_id), resno,
						    inscode);
      
      if (istat) {
	 g.update_go_to_atom_window_on_changed_mol(imol);
	 graphics_draw();
      }
   }

   std::string cmd = "delete-residue-sidechain";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(resno);
   args.push_back(ins_code);
   add_to_history_typed(cmd, args);
}

// ---------------------------------------------------------------------
//                 rotamer
// ---------------------------------------------------------------------
// 


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
   std::string cmd = "setup-rotamers";
   std::vector<coot::command_arg_t> args;
   args.push_back(state);
   add_to_history_typed(cmd, args);
} 


void do_rotamers(int atom_index, int imol) {

//    std::cout << "     Rotamer library:" << std::endl;
//    std::cout << "     R. L. Dunbrack, Jr. and F. E. Cohen." << std::endl;
//    std::cout << "     Bayesian statistical analysis of ";
//    std::cout << "protein sidechain rotamer preferences" << std::endl;
//    std::cout << "     Protein Science, 6, 1661-1681 (1997)." << std::endl;
//    std::cout << "" << std::endl;

   graphics_info_t g;
   g.do_rotamers(atom_index, imol); 
   std::string cmd = "do-rotamers";
   std::vector<coot::command_arg_t> args;
   args.push_back(atom_index);
   args.push_back(imol);
   add_to_history_typed(cmd, args);
}

// same as do_rotamers, except, a better name and we give residue
// specs, so that we can use the active residue.
void show_rotamers_dialog(int imol, const char *chain_id, int resno, const char *ins_code, const char *altconf) {

   int atom_index = -1;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      atom_index = g.molecules[imol].atom_index_first_atom_in_residue(chain_id, resno, ins_code, altconf);
      if (atom_index != -1) {
	 g.do_rotamers(atom_index, imol); 
      } else {
	 std::cout << "No atom index found in molecule " << imol << std::endl;
      }
   }
} 


void
set_rotamer_lowest_probability(float f) {
#ifdef USE_DUNBRACK_ROTAMERS
   graphics_info_t g;
   g.rotamer_lowest_probability = f;
   std::string cmd = "set-rotamer-lowest-probability";
   std::vector<coot::command_arg_t> args;
   args.push_back(f);
   add_to_history_typed(cmd, args);
#endif    
}

void
set_rotamer_check_clashes(int i) {
   graphics_info_t::rotamer_fit_clash_flag = i;
   std::string cmd = "set-rotamer-check-clashes";
   std::vector<coot::command_arg_t> args;
   args.push_back(i);
   add_to_history_typed(cmd, args);
}

// Return -999 on imol indexing error
// -99.9 in class function error
// 
float
auto_fit_best_rotamer(int imol_coords,
                      const char *chain_id,
                      int resno,
                      const char *insertion_code,
                      const char *altloc,
                      int imol_map, int clash_flag, float lowest_probability) {

   float f = -999.9;

   if (is_valid_model_molecule(imol_coords)) {

      std::string ins(insertion_code);
      std::string chain(chain_id);
      graphics_info_t g;
      int mode = graphics_info_t::rotamer_search_mode;
      if (! is_valid_map_molecule(imol_map)) {
         std::cout << "INFO:: fitting rotamers by clash score only " << std::endl;
         imol_map = -1;
         f = graphics_info_t::molecules[imol_coords].auto_fit_best_rotamer(mode,
                                                                           resno, altloc, ins,
                                                                           chain, imol_map,
                                                                           1,
                                                                           lowest_probability,
                                                                           *g.Geom_p());
      } else {
         f = g.molecules[imol_coords].auto_fit_best_rotamer(mode,
                                                            resno, altloc, ins,
                                                            chain, imol_map,
                                                            clash_flag,
                                                            lowest_probability,
                                                            *g.Geom_p());

         // first do real space refine if requested
         if (g.rotamer_auto_fit_do_post_refine_flag) {
            // Run refine zone with autoaccept, autorange on
            // the "clicked" atom:
            // BL says:: dont think we do autoaccept!?
            //
            // 20230115-PE Yes. This needs an auto-accept wrapper.
            //
            refine_auto_range(imol_coords, chain_id, resno, altloc);
         }

         // get the residue so that it can update the geometry graph
         mmdb::Residue *residue_p = g.molecules[imol_coords].get_residue(chain, resno, ins);
         if (residue_p)
            g.update_geometry_graphs(&residue_p, 1, imol_coords, imol_map); // in auto_fit_best_rotamer
         std::cout << "Fitting score for best rotamer: " << f << std::endl;
      }
      g.run_post_manipulation_hook(imol_coords, 0);
      graphics_draw();
   }
   std::string cmd = "auto-fit-best-rotamer";
   std::vector<coot::command_arg_t> args;
   args.push_back(resno);
   args.push_back(coot::util::single_quote(altloc));
   args.push_back(coot::util::single_quote(insertion_code));
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(imol_coords);
   args.push_back(imol_map);
   args.push_back(clash_flag);
   args.push_back(lowest_probability);
   add_to_history_typed(cmd, args);
   return f;
}


/*! auto-fit the rotamer for the active residue */
float auto_fit_rotamer_active_residue() {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      const auto &atom_spec = pp.second.second;
      const std::string &chain_id = atom_spec.chain_id;
      int res_no = atom_spec.res_no;
      const std::string &ins_code = atom_spec.ins_code;
      const std::string &alt_loc = atom_spec.alt_conf;
      int imol_map = imol_refinement_map();
      bool clash_flag = 1; // true
      float lowest_probability = 2.0;
      return auto_fit_best_rotamer(imol, chain_id.c_str(), res_no, ins_code.c_str(), alt_loc.c_str(), imol_map, clash_flag, lowest_probability);
   } else {
      return -1;
   }
}




void
set_auto_fit_best_rotamer_clash_flag(int i) { /* 0 off, 1 on */
   graphics_info_t::rotamer_fit_clash_flag = i;
   std::string cmd = "set-auto-fit-best-rotamer-clash-flag";
   std::vector<coot::command_arg_t> args;
   args.push_back(i);
   add_to_history_typed(cmd, args);
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
   std::string cmd = "setup-auto-fit-rotamer";
   std::vector<coot::command_arg_t> args;
   args.push_back(state);
   add_to_history_typed(cmd, args);
}


// FIXME  (autofit rotamer seems to return a score OK).
float
rotamer_score(int imol, const char *chain_id, int res_no, const char *insertion_code,
	      const char *alt_conf) {

   float f = 0;

   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      mmdb::Residue *residue_p =
	 graphics_info_t::molecules[imol].get_residue(chain_id, res_no, insertion_code);
      if (residue_p) {
	 float lp = graphics_info_t::rotamer_lowest_probability;
	 graphics_info_t g;
	 coot::rotamer_probability_info_t d_score = 
	    g.get_rotamer_probability(residue_p, alt_conf, mol, lp, 1);
	 if (d_score.state == coot::rotamer_probability_info_t::OK) 
	    f = d_score.probability;
      }
   }
   
   
   std::string cmd = "rotamer-score";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(res_no);
   args.push_back(coot::util::single_quote(insertion_code));
   add_to_history_typed(cmd, args);
   return f;
}


/*! \brief return the number of rotamers for this residue */
int n_rotamers(int imol, const char *chain_id, int resno, const char *ins_code) {

   int r = -1; 
   if (is_valid_model_molecule(imol)) { 
      mmdb::Residue *res = graphics_info_t::molecules[imol].get_residue(chain_id, resno, ins_code);
      if (res) {
	 graphics_info_t g;
#ifdef USE_DUNBRACK_ROTAMERS
	 coot::dunbrack d(res, g.molecules[imol].atom_sel.mol, g.rotamer_lowest_probability, 0);
#else			
	 std::string alt_conf = "";
	 coot::richardson_rotamer d(res, alt_conf, g.molecules[imol].atom_sel.mol,
				    g.rotamer_lowest_probability, 0);
#endif // USE_DUNBRACK_ROTAMERS
	 
	 std::vector<float> probabilities = d.probabilities();
	 r = probabilities.size();
      }
   }
   return r;
} 

/*! \brief set the residue specified to the rotamer number specifed. */
int set_residue_to_rotamer_number(int imol, const char *chain_id, int resno, const char *ins_code,
				  const char *alt_conf, int rotamer_number) {

   int i_done = 0;
   if (is_valid_model_molecule(imol)) {
      int n = rotamer_number;
      coot::residue_spec_t res_spec(chain_id, resno, ins_code);
      graphics_info_t g;
      i_done = g.molecules[imol].set_residue_to_rotamer_number(res_spec, alt_conf, n, *g.Geom_p());
      graphics_draw();
   }
   return i_done; 
}

int set_residue_to_rotamer_name(int imol, const char *chain_id, int resno, const char *ins_code,
				const char *alt_conf,
				const char *rotamer_name) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t res_spec(chain_id, resno, ins_code);
      graphics_info_t g;
      status = g.molecules[imol].set_residue_to_rotamer_name(res_spec, alt_conf, rotamer_name, *g.Geom_p());
      graphics_draw();
   } 
   return status;
} 


#ifdef USE_GUILE
SCM get_rotamer_name_scm(int imol, const char *chain_id, int resno, const char *ins_code) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t res_spec(chain_id, resno, ins_code);
      mmdb::Residue *res = graphics_info_t::molecules[imol].get_residue(res_spec);
      if (res) {
#ifdef USE_DUNBRACK_ROTAMERS
#else
	 // we are not passed an alt conf.  We should be, shouldn't we?
	 std::string alt_conf = "";
	 mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
	 coot::richardson_rotamer d(res, alt_conf, mol, 0.0, 1);
	 coot::rotamer_probability_info_t prob = d.probability_of_this_rotamer();
	 std::cout << "INFO:: " << coot::residue_spec_t(res) << " " << prob << std::endl;
	 r = scm_from_locale_string(prob.rotamer_name.c_str());
#endif      
      }
   }
   return r;
} 
#endif 

#ifdef USE_PYTHON
PyObject *get_rotamer_name_py(int imol, const char *chain_id, int resno, const char *ins_code) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t res_spec(chain_id, resno, ins_code);
      mmdb::Residue *res = graphics_info_t::molecules[imol].get_residue(res_spec);
      if (res) {
#ifdef USE_DUNBRACK_ROTAMERS
#else
	 // we are not passed an alt conf.  We should be, shouldn't we?
	 std::string alt_conf = "";
	 coot::richardson_rotamer d(res, alt_conf, graphics_info_t::molecules[imol].atom_sel.mol,
				    0.0, 1);
	 coot::rotamer_probability_info_t prob = d.probability_of_this_rotamer();
	 r = myPyString_FromString(prob.rotamer_name.c_str());
#endif      
      }
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
} 
#endif 




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
   std::string cmd = "setup-mutate";
   std::vector<coot::command_arg_t> args;
   args.push_back(state);
   add_to_history_typed(cmd, args);
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
	 g.show_select_map_frame();
	 g.in_mutate_auto_fit_define = 0;
	 normal_cursor();
	 g.model_fit_refine_unactive_togglebutton("model_refine_dialog_mutate_auto_fit_togglebutton");
      }
   } else {
      g.in_mutate_auto_fit_define = state;
      g.normal_cursor();
   }
   std::string cmd = "setup-mutate-auto-fit";
   std::vector<coot::command_arg_t> args;
   args.push_back(state);
   add_to_history_typed(cmd, args);
}


/* 1 for yes, 0 for no. */
void set_mutate_auto_fit_do_post_refine(short int istate) {

   graphics_info_t::mutate_auto_fit_do_post_refine_flag = istate;
   std::string cmd = "set-mutate-auto-fit-do-post-refine";
   std::vector<coot::command_arg_t> args;
   args.push_back(istate);
   add_to_history_typed(cmd, args);
} 

/*! \brief what is the value of the previous flag? */
int mutate_auto_fit_do_post_refine_state() {
   add_to_history_simple("mutate-auto-fit-do-post-refine-state");
   return graphics_info_t::mutate_auto_fit_do_post_refine_flag;
} 

/* 1 for yes, 0 for no. */
void set_rotamer_auto_fit_do_post_refine(short int istate) {

   graphics_info_t::rotamer_auto_fit_do_post_refine_flag = istate;
   std::string cmd = "set-rotamer-auto-fit-do-post-refine";
   std::vector<coot::command_arg_t> args;
   args.push_back(istate);
   add_to_history_typed(cmd, args);
} 

/*! \brief what is the value of the previous flag? */
int rotamer_auto_fit_do_post_refine_state() {
   add_to_history_simple("rotamer-auto-fit-do-post-refine-state");
   return graphics_info_t::rotamer_auto_fit_do_post_refine_flag;
} 

/*! \brief set a flag saying that the chosen residue should only be
  added as a stub (mainchain + CB) */
void set_residue_type_chooser_stub_state(short int istat) {
   graphics_info_t::residue_type_chooser_stub_flag = istat;
   std::string cmd = "set-residue-type-chooser-stub-state";
   std::vector<coot::command_arg_t> args;
   args.push_back(istat);
   add_to_history_typed(cmd, args);
}

void handle_residue_type_chooser_entry_chose_type(const char *entry_text,
                                                  short int stub_mode) {

   // mutate_and_auto_fit

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      coot::atom_spec_t atom_spec = pp.second.second;
      if (entry_text) {
         std::string et(entry_text);
         if (! et.empty()) {
            char c = et[0];
            c = std::toupper(c);
            std::string res_type = coot::util::single_letter_to_3_letter_code(c);
            if (coot::util::is_standard_amino_acid_name(res_type)) {
               graphics_info_t g;
               // set the imol and the atom index - bleugh (it's because we are
               // using a function called from graphics-info-defines.cc - which
               // uses an atom click. Sigh.
               mmdb::Atom *at = g.molecules[imol].get_atom(atom_spec);
               if (at) {
                  int atom_index = -1;
                  int UDDAtomIndexHandle = g.molecules[imol].atom_sel.UDDAtomIndexHandle;
                  int status = at->GetUDData(UDDAtomIndexHandle, atom_index);
                  if (status == mmdb::UDDATA_Ok) {
                     g.mutate_auto_fit_residue_atom_index = atom_index;
                     g.mutate_auto_fit_residue_imol = imol;
                     coot::residue_spec_t res_spec = coot::residue_spec_t(coot::atom_spec_t(at));
                     g.do_mutation(imol, res_spec, res_type, stub_mode);
                  } else {
                     std::cout << "UDData not OK " << std::endl;
                  }
               } else {
                  std::cout << "No atom in molecule" << std::endl;
               }
            } else {
               std::cout << "Non-standard residue " << c << " " << res_type << std::endl;
            }
         } else {
            std::cout << "empty entry text" << std::endl;
         }
      } else {
         std::cout << "No entry text" << std::endl;
      }
   } else {
      std::cout << "No active atom" << std::endl;
   }
}



void
do_mutation(const char *type, short int stub_button_state_flag) {
   graphics_info_t g;
   // use g.mutate_residue_atom_index and g.mutate_residue_imol
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      coot::residue_spec_t res_spec = coot::residue_spec_t(pp.second.second);
      g.do_mutation(imol, res_spec, type, stub_button_state_flag);
   }
   std::string cmd = "do-mutation";
   std::vector<coot::command_arg_t> args;
   args.push_back(coot::util::single_quote(type));
   args.push_back(stub_button_state_flag);
   add_to_history_typed(cmd, args);
}

void mutate_active_residue() {

   graphics_info_t g;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      bool show_stub_flag = true;
      GtkWidget *widget = g.wrapped_create_residue_type_chooser_window(show_stub_flag);
      gtk_widget_set_visible(widget, TRUE);
      g.in_mutate_auto_fit_define = 0;
      g.residue_type_chooser_auto_fit_flag = 1;
      g.pick_pending_flag = 0;
   }

}

void
mutate_active_residue_to_single_letter_code(const std::string &slc) {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      coot::atom_spec_t atom_spec = pp.second.second;
      std::string u_slc = coot::util::upcase(slc);
      std::string res_type = coot::util::single_letter_to_3_letter_code(u_slc);
      std::cout << "debug:: single_letter_to_3_letter_code(): scl: " << u_slc << " res_type: " << res_type << std::endl;
      if (! res_type.empty()) {
         graphics_info_t g;
         int ires = atom_spec.res_no;
         std::string ins_code = atom_spec.ins_code;
         std::string chain_id = atom_spec.chain_id;
         g.molecules[imol].mutate(ires, ins_code, std::string(chain_id), res_type);
         g.update_validation(imol);
         graphics_draw();
      }
   }
}

// \brief show keyboard mutate dialog
void show_keyboard_mutate_frame() {

   GtkWidget *w = widget_from_builder("keyboard_mutate_frame");
   GtkWidget *e = widget_from_builder("keyboard_mutate_entry");
   gtk_widget_grab_focus(e);
   gtk_widget_set_visible(w, TRUE);

}


//! mutate by overlap - 20241211-PE welcome into Coot at long last.
int mutate_by_overlap(int imol, const std::string &chain_id, int res_no, const std::string &new_type) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      status = g.molecules[imol].mutate_by_overlap(chain_id, res_no, new_type);
   }
   return status;
}


// return success on residue type match
// success: 1, failure: 0.
int
mutate_internal(int ires_serial, const char *chain_id, int imol, const std::string &target_res_type) {

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
mutate(int imol, const char *chain_id, int ires, const char *inscode,  const char *target_res_type) {

   int istate = 0;
   std::string target_type(target_res_type);

   if (is_valid_model_molecule(imol)) {
      istate = graphics_info_t::molecules[imol].mutate(ires, inscode, std::string(chain_id), std::string(target_res_type));
      graphics_info_t g;
      g.update_validation(imol);
      graphics_draw();
   }
   std::string cmd = "mutate";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(ires);
   args.push_back(coot::util::single_quote(inscode));
   args.push_back(coot::util::single_quote(target_res_type));
   add_to_history_typed(cmd, args);

   return istate;
}

// return success status.
int mutate_base(int imol, const char *chain_id, int res_no, const char *ins_code, const char *res_type) {

   int istate = 0;
   graphics_info_t g;
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t r(chain_id, res_no, ins_code);
      bool use_old_style_naming = false;
      istate = g.molecules[imol].mutate_base(r, res_type, use_old_style_naming);
      graphics_draw();
   }
   std::string cmd = "mutate-base";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(res_no);
   args.push_back(coot::util::single_quote(ins_code));
   args.push_back(coot::util::single_quote(res_type));
   add_to_history_typed(cmd, args);

   return istate;
}



// Return success on residue type match
// success: 1, failure: 0.
//
// Does not cause a make_backup().
//
int
mutate_single_residue_by_serial_number(int ires_serial, const char *chain_id, int imol,
			   char target_res_type) {

   std::string target_as_str = coot::util::single_letter_to_3_letter_code(target_res_type);
   std::cout << "INFO:: mutate target_res_type :" << target_as_str << ":" << std::endl;

   return mutate_internal(ires_serial, chain_id, imol, target_as_str);

}

// Previously, I was using mutate_single_residue_by_seqno to be a
// wrapper for mutate_single_residue_by_serial_number.
//
// But that fails when the residue is perfectly reasonably except that
// the serial number is -1 (I don't know wny this happens but it does
// in terminal residue addition).  So I will need new functionally
// that does the residue at root by seqnum not serial_number.
//
int mutate_single_residue_by_seqno(int imol, const char *chain_id, int ires, const char *inscode,
				   char target_res_type) {

   int status = -1;
   std::string target_as_str = coot::util::single_letter_to_3_letter_code(target_res_type);

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      status = g.molecules[imol].mutate(ires, std::string(inscode),
					std::string(chain_id), target_as_str);
   }
   return status;
}

int mutate_residue_range(int imol, const std::string &chain_id, int res_no_start, int res_no_end,
                         const std::string &target_sequence) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {

      int delta_residue_range = res_no_end - res_no_start  + 1;
      int target_sequence_length = target_sequence.length();
      if (target_sequence_length == delta_residue_range) {
         for (int ii=0; ii<target_sequence_length; ii++) {
            int res_no = res_no_start + ii;
            std::string ins_code;
            char target_res_type = target_sequence[ii];
            mutate_single_residue_by_seqno(imol, chain_id.c_str(), res_no, ins_code.c_str(), target_res_type);
            status = 1;
         }
      } else {
         // log this
         std::cout << "WARNING:: in mutate_residue_range() mismatch sequence and range "
                   << delta_residue_range << " " << target_sequence_length << std::endl;
      }
   }
   return status;
}

int mutate_and_autofit_residue_range(int imol, const char *chain_id, int start_res_no, int stop_res_no,
                                     const char *target_sequence) {

  int status = 0;
  if (is_valid_model_molecule(imol)) {
     int imol_map = imol_refinement_map();
     int delta_residue_range = stop_res_no - start_res_no  + 1;
     int target_sequence_length = std::string(target_sequence).length();
     if (target_sequence_length == delta_residue_range) {
        int backup_mode = backup_state(imol);
        if (backup_mode)
           turn_off_backup(imol);
        mutate_residue_range(imol, chain_id, start_res_no, stop_res_no, std::string(target_sequence));
        int clash_flag = true;
        float lowest_probability = 0.3; // changed from 0.5
        for (int ii=0; ii<target_sequence_length; ii++) {
           int res_no = start_res_no + ii;
           std::string ins_code;
           std::string alt_conf;
           auto_fit_best_rotamer(imol, chain_id, res_no, ins_code.c_str(), alt_conf.c_str(),
                                 imol_map, clash_flag, lowest_probability);
        }
        if (backup_mode)
           turn_on_backup(imol);
     }
  }
  return status;
}



/* push the residues along a bit

e.g. if nudge_by is 1, then the sidechain of residue 20 is moved up
onto what is currently residue 21.  The mainchain numbering and atoms is not changed. */
int nudge_residue_sequence(int imol, const char *chain_id, int res_no_range_start,
			   int res_no_range_end,
			   int nudge_by,
			   short int nudge_residue_numbers_also) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = graphics_info_t::molecules[imol].nudge_residue_sequence(chain_id,
								       res_no_range_start,
								       res_no_range_end,
								       nudge_by,
								       nudge_residue_numbers_also);
      if (status)
	 graphics_draw();
   }
   return status;
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
      g.add_status_bar_text("Click on a atom. The clicked atom affects the torsion's wagging dog/tail...");
      g.pick_pending_flag = 1;
   } else {
      g.in_edit_chi_angles_define = 0;
   }
   std::string cmd = "setup-edit-chi-angles";
   std::vector<coot::command_arg_t> args;
   args.push_back(state);
   add_to_history_typed(cmd, args);
}

void rotate_chi(float am) {

   graphics_info_t g;
   if (g.in_edit_chi_mode_flag || g.in_edit_torsion_general_flag) {
      g.rotate_chi(am, am);
   }
} 


void setup_torsion_general(short int state) {

   graphics_info_t g;
   if (state) {
      g.in_torsion_general_define = 1;
      g.pick_cursor_maybe();
      g.add_status_bar_text("Click on a atom. The order of the clicked atoms affects the torsion's wagging dog/tail...");
      g.pick_pending_flag = 1;
   } else {
      g.in_torsion_general_define = 0;
   }
   std::string cmd = "setup-torsion-general";
   std::vector<coot::command_arg_t> args;
   args.push_back(state);
   add_to_history_typed(cmd, args);
}

void toggle_torsion_general_reverse()  { /* a bool really */
   if (graphics_info_t::torsion_general_reverse_flag)
      graphics_info_t::torsion_general_reverse_flag = 0;
   else 
      graphics_info_t::torsion_general_reverse_flag = 1;
}

void setup_residue_partial_alt_locs(short int state) {

   graphics_info_t g;
   g.in_residue_partial_alt_locs_define = state;
   g.pick_cursor_maybe();
   g.add_status_bar_text("Click on an atom to identify the residue for alt confs");

   std::string cmd = "setup-residue-partial-alt-locs";
   std::vector<coot::command_arg_t> args;
   args.push_back(state);
   add_to_history_typed(cmd, args);
}


int set_show_chi_angle_bond(int imode) {

   graphics_info_t::draw_chi_angle_flash_bond_flag = imode;
   graphics_draw();
   std::string cmd = "set-show-chi-angle-bond";
   std::vector<coot::command_arg_t> args;
   args.push_back(imode);
   add_to_history_typed(cmd, args);
   return 0; // should be a void function, I imagine.
} 


/* a callback from the callbacks.c, setting the state of
   graphics_info_t::edit_chi_angles_reverse_fragment flag */
void set_edit_chi_angles_reverse_fragment_state(short int istate) {

   graphics_info_t::edit_chi_angles_reverse_fragment = istate;

} 


// Set a flag: Should torsions that move hydrogens be
// considered/displayed in button box?
// 
void set_find_hydrogen_torsions(short int state) {
   graphics_info_t g;
   g.find_hydrogen_torsions_flag = static_cast<bool>(state);
   std::string cmd = "set-find-hydrogen-torsion";
   std::vector<coot::command_arg_t> args;
   args.push_back(state);
   add_to_history_typed(cmd, args);
}

void set_graphics_edit_current_chi(int ichi) { /* button callback */

   graphics_info_t g;
   graphics_info_t::edit_chi_current_chi = ichi;
   if (ichi == 0) {
      graphics_info_t::in_edit_chi_mode_flag = 0; // off
   } else {
      graphics_info_t::in_edit_chi_mode_flag = 1; // on
      g.setup_flash_bond_using_moving_atom_internal(ichi-1);
   }

}

void unset_moving_atom_move_chis() { 
   graphics_info_t::moving_atoms_move_chis_flag = 0; // keyboard 1,2,3
						     // etc cant put
						     // graphics/mouse
						     // into rotate
						     // chi mode.
} 

// BL says:: maybe this should rather go via a set state function!?
void set_moving_atom_move_chis() { 
   graphics_info_t::moving_atoms_move_chis_flag = 1; // keyboard 1,2,3
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
	 // 20090815
	 // Do I want this?
	 // g.chi_angle_alt_conf = molecules[imol].atom_sel.atom_selection[atom_index]->altLoc;
	 // or this?
	 g.chi_angle_alt_conf = altconf;
	 g.execute_edit_chi_angles(atom_index, imol);
	 status = 1;
      }
   }
   return status;
}




/*  ----------------------------------------------------------------------- */
/*                  sequence (assignment)                                   */
/*  ----------------------------------------------------------------------- */
/* section Sequence (Assignment) */

void assign_sequence(int imol_coords, int imol_map, const char *chain_id) {

   if (is_valid_model_molecule(imol_coords))
      if (is_valid_map_molecule(imol_map))
	 graphics_info_t::molecules[imol_coords].assign_sequence(graphics_info_t::molecules[imol_map].xmap, std::string(chain_id));

  std::string cmd = "assign-sequence";
  std::vector<coot::command_arg_t> args;
  args.push_back(imol_coords);
  args.push_back(imol_map);
  args.push_back(single_quote(chain_id));
  add_to_history_typed(cmd, args);
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
      if (idx > 0 && idx < g.molecules[imol].atom_sel.n_selected_atoms) {
	 mmdb::Atom *at = graphics_info_t::molecules[imol].atom_sel.atom_selection[idx];
	 mmdb::Residue *r = at->residue;
	 if (r) {
	    std::string cbn = "";
	    if (coot::util::nucleotide_is_DNA(r)) {
	       cbn = coot::util::canonical_base_name(type, coot::DNA);
	    } else {
	       cbn = coot::util::canonical_base_name(type, coot::RNA);
	    }

	    if (cbn != "") {
	       bool old = g.convert_to_v2_atom_names_flag;
	       coot::residue_spec_t res_spec(r);
	       int istat = graphics_info_t::molecules[imol].mutate_base(res_spec, cbn, old);
	       if (istat)
		  graphics_draw();
	       // Is this the right function?
	       update_go_to_atom_window_on_changed_mol(imol);
	    } else {
	       std::string s = "No canonical base name found";
	       std::cout << "WARNING:: " << s << std::endl;
	       add_status_bar_text(s.c_str());
	    } 
	 }
      } else {
	 std::cout << "ERROR:: out of range atom index " << idx << std::endl;
      }
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
      g.add_status_bar_text(s);
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
      g.add_status_bar_text("Click on an atom in the residue that you want to flip");
      g.pick_pending_flag = 1;
   } else {
      g.normal_cursor();
      g.pick_pending_flag = 0;
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
      coot::util::missing_atom_info m_i_info =
	 g.molecules[imol].fill_partial_residues(g.Geom_p(), imol_map);

      // We do refinement here because we can't do refinement (easily) from inside molecule_class_info_t.
      // Not today, anyway.

      if (imol_map > -1) {

	 int backup_mode = backup_state(imol);
	 if (backup_mode)
	    g.molecules[imol].make_backup_from_outside();
	 turn_off_backup(imol);

	 int refinement_replacement_state = refinement_immediate_replacement_state();
	 set_refinement_immediate_replacement(1);

         std::string alt_conf("");
         std::vector<mmdb::Residue *> residues;

      	 for (unsigned int i=0; i<m_i_info.residues_with_missing_atoms.size(); i++) {
      	    int resno =  m_i_info.residues_with_missing_atoms[i]->GetSeqNum();
      	    std::string chain_id = m_i_info.residues_with_missing_atoms[i]->GetChainID();
      	    std::string residue_type = m_i_info.residues_with_missing_atoms[i]->GetResName();
      	    std::string inscode = m_i_info.residues_with_missing_atoms[i]->GetInsCode();
            residues.push_back(m_i_info.residues_with_missing_atoms[i]);
      	 }
         mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
         coot::refinement_results_t rr = g.refine_residues_vec(imol, residues, alt_conf.c_str(), mol);
         c_accept_moving_atoms();
	 set_refinement_immediate_replacement(refinement_replacement_state);

	 if (backup_mode)
	    turn_on_backup(imol);

      } else {
	 g.show_select_map_frame();
      }
      graphics_draw();
   }
}

void fill_partial_residue(int imol, const char *chain_id, int resno, const char* inscode) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      int imol_map = g.Imol_Refinement_Map();
      if (imol_map > -1) {
	 coot::residue_spec_t rs(chain_id, resno, inscode);
	 g.molecules[imol].fill_partial_residue(rs, g.Geom_p(), imol_map);
	 // post process...
	 int refinement_replacement_state = refinement_immediate_replacement_state();
	 set_refinement_immediate_replacement(1);
	 std::string altconf("");
	 short int is_water = 0;
	 // hmmm backups are being done....
	 g.refine_residue_range(imol, chain_id, chain_id, resno, inscode, resno, inscode,
				altconf, is_water);
	 accept_regularizement();
	 set_refinement_immediate_replacement(refinement_replacement_state);

      } else {
	 g.show_select_map_frame();
      }
   }
}

// Fill amino acid residues, do backrub rotamer search for residues, but don't do refinement
//
void simple_fill_partial_residues(int imol) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      int imol_map = g.Imol_Refinement_Map();
      if (imol_map > -1) {
	 coot::util::missing_atom_info m_i_info =
	    g.molecules[imol].fill_partial_residues(g.Geom_p(), imol_map);
	 graphics_draw();
      } else {
	 add_status_bar_text("Refinement map not set");
      }
   }
}


#include "c-interface-sequence.hh"
#include "ligand/side-chain-densities.hh"

std::string sequence_from_map(int imol, const std::string &chain_id, int resno_start, int resno_end, int imol_map) {

   std::string guessed_sequence;
   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {
         mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
         const clipper::Xmap<float> &xmap = graphics_info_t::molecules[imol_map].xmap;
         coot::side_chain_densities scd;
         scd.fill_residue_blocks(mol, chain_id, resno_start, resno_end, xmap);
         guessed_sequence =
            scd.guess_the_sequence(mol, chain_id, resno_start, resno_end, xmap);
         // std::cout << "guessed sequence " << guessed_sequence << std::endl;
      }
   }
   return guessed_sequence;
}

void apply_fasta_multi_to_fragment(int imol, const std::string &chain_id, int resno_start, int resno_end, int imol_map,
                                   const coot::fasta_multi &fam) {

   // change the residue number if you can

   if (is_valid_model_molecule(imol)) {
      if (is_valid_map_molecule(imol_map)) {
         mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
         const clipper::Xmap<float> &xmap = graphics_info_t::molecules[imol_map].xmap;
         coot::side_chain_densities scd;
         unsigned int n_sequences = fam.size();
         if (n_sequences > 0) {
            for (unsigned int idx=0; idx<n_sequences; idx++) {
               std::string sequence = fam[idx].sequence;
               const std::string &name = fam[idx].name;
               std::pair<std::string, std::vector<mmdb::Residue *> > a_run_of_residues =
                  scd.setup_test_sequence(mol, chain_id, resno_start, resno_end, xmap);
               // scd.test_sequence(mol, chain_id, resno_start, resno_end, xmap, name, sequence);
               if (a_run_of_residues.first.empty()) {
                  bool print_slider_results = true;
                  scd.test_sequence(a_run_of_residues.second, xmap, name, sequence, print_slider_results);
               } else {
                  std::cout << "ERROR:: when generating a run of residues " << std::endl;
                  std::cout << a_run_of_residues.first << std::endl;
                  add_status_bar_text(a_run_of_residues.first.c_str());
               }
            }

            bool only_probable = false; // just give me your best shot - I know that it might be wrong.
            bool print_sequencing_solutions_flag = true;
            coot::side_chain_densities::results_t new_sequence_result = scd.get_result(only_probable, print_sequencing_solutions_flag);
            std::string new_sequence = new_sequence_result.sequence;
            std::cout << "debug:: apply_fasta_multi_to_fragment(): get_result(): new_sequence " << new_sequence << std::endl;
            int offset = new_sequence_result.offset;
            if (! new_sequence.empty()) {
               int sl = new_sequence.length();
               int residue_count = resno_end - resno_start + 1;
               std::cout << "DEBUG:: new_sequence length " << sl << " residue_count " << residue_count
                         << " offset " << offset << std::endl;
               if (sl == residue_count) {
                  molecule_class_info_t &m = graphics_info_t::molecules[imol];
                  m.make_backup_from_outside();
                  bool backup_state = m.backups_state();
                  m.turn_off_backup();
                  int ires_serial_first = m.residue_serial_number(chain_id, resno_start, "");
                  if (ires_serial_first != -1) {
                     for (int ires=resno_start; ires<=resno_end; ires++) {
                        int idx_offset = ires-resno_start;
                        // need to convert ires to a serial number
                        int ires_serial = ires_serial_first + idx_offset;
                        char letter = new_sequence[idx_offset];
                        std::string new_residue_type = coot::util::single_letter_to_3_letter_code(letter);
                        m.mutate_single_multipart(ires_serial, chain_id, new_residue_type);
                     }
                     m.fill_partial_residues(graphics_info_t::Geom_p(), imol_map);
                     m.backrub_rotamer_residue_range(chain_id, resno_start, resno_end, *graphics_info_t::Geom_p());
                     // now, can I change the residue numbers?
                     // residue-1 in the target (sequence) number scheme has offset 0
                     int resno_offset = offset + 1 - resno_start;
                     m.renumber_residue_range(chain_id, resno_start, resno_end, resno_offset);
                  } else {
                     std::cout << "WARNING:: apply_sequence_to_fragment() failed to find serial number of residue "
                               << chain_id << " with res-no " << resno_start << std::endl;
                  }
                  if (backup_state)
                     m.turn_on_backup();
               } else {
                  std::cout << "WARNING:: residue count (based on given resno start and end) does not match "
                            << "new sequence length" << std::endl;
               }
            }
            graphics_draw();
         }
      } else {
         std::cout << "WARNING:: not a valid map " << imol_map << std::endl;
      }

   } else {
      std::cout << "WARNING:: not a valid model molecule " << imol << std::endl;
   }
}


void apply_sequence_to_fragment(int imol, const std::string &chain_id, int resno_start, int resno_end, int imol_map,
                                const std::string &multi_sequence_file_name) {

   coot::fasta_multi fam(multi_sequence_file_name);
   apply_fasta_multi_to_fragment(imol, chain_id, resno_start, resno_end, imol_map, fam);
}


void assign_sequence_to_active_fragment() {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      coot::atom_spec_t atom_spec = pp.second.second;
      int imol_map = imol_refinement_map();
      if (is_valid_model_molecule(imol)) {
         if (is_valid_map_molecule(imol_map)) {
            const molecule_class_info_t &m = graphics_info_t::molecules[imol];
            mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
            coot::residue_spec_t residue_spec(atom_spec);
            mmdb::Residue *residue_p = graphics_info_t::molecules[imol].get_residue(residue_spec);
            if (residue_p) {
               float close_dist_max = 1.7;
               // Yes, simple_residue_tree() is very slow: ~2 seconds for 1000 residues
               // std::cout << "getting the residue tree ... " << std::endl;
               std::vector<mmdb::Residue *> v = coot::simple_residue_tree(residue_p, mol, close_dist_max);
               // std::cout << "done getting the residue tree ... " << std::endl;
               if (! v.empty()) {
                  std::string chain_id = atom_spec.chain_id;
                  int resno_low  =  10000000;
                  int resno_high = -10000000;
                  for (unsigned int i=0; i<v.size(); i++) {
                     int resno = v[i]->GetSeqNum();
                     // check that v[i] is in the same chain as the active atom
                     if (v[i]->GetChain() == residue_p->GetChain()) {
                        if (resno < resno_low)  resno_low  = resno;
                        if (resno > resno_high) resno_high = resno;
                     }
                  }
                  coot::fasta_multi fam;
                  std::vector<std::pair<std::string, std::string> > sequence_info = m.sequence_info();
                  for (auto seq : sequence_info) {
                     coot::fasta f(seq.first, seq.second, coot::fasta::SIMPLE_STRING);
                     fam.add(f);
                  }
                  std::cout << "debug:: calling apply_fasta_multi_to_fragment() " << chain_id
                            << " " << resno_low << " " << resno_high << " with imol_map: " << imol_map << std::endl;
                  apply_fasta_multi_to_fragment(imol, chain_id, resno_low, resno_high, imol_map, fam);
               } else {
                  std::cout << "empty v from simple_residue_tree() " << std::endl;
               }
            } else {
               std::cout << "residue not found in molecules " << residue_spec << std::endl;
            }
         } else {
            std::cout << "Not a valid map molecule " << imol_map << std::endl;
         }
      } else {
         std::cout << "Not a valid model molecule " << imol << std::endl;
      }
   } else {
      std::cout << "No active atom" << std::endl;
   }
}


//! @return a python list of residue specs for the residues in the given chain
PyObject *get_residues_in_chain_py(int imol, const std::string &chain_id) {

   PyObject *l = PyList_New(0);
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      if (mol) {
         int imod = 1;
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               if (chain_id == chain_p->GetChainID()) {
                  int n_res = chain_p->GetNumberOfResidues();
                  std::vector<coot::residue_spec_t> res_specs;
                  for (int ires=0; ires<n_res; ires++) {
                     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                     if (residue_p) {
                        coot::residue_spec_t rs(residue_p);
                        res_specs.push_back(rs);
                     }
                  }
                  if (! res_specs.empty()) {
                     // delete l here
                     l = PyList_New(res_specs.size());
                     for (unsigned int i=0; i<res_specs.size(); i++) {
                        const auto &rs = res_specs[i];
                        PyObject *o = residue_spec_to_py(rs);
                        PyList_SetItem(l, i, o);
                     }
                  }
               }
            }
         }
      }
   }
   return l;
}
