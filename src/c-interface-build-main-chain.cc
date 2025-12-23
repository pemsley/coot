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

#ifdef USE_PYTHON
#include <Python.h>  // before system includes to stop "POSIX_C_SOURCE" redefined problems
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
#include <windows.h>
#endif


#include <mmdb2/mmdb_manager.h>

#include "coords/mmdb-extras.hh"
#include "coords/mmdb.hh"
#include "coords/Cartesian.hh"
#include "coords/Bond_lines.hh"
#include "coords/mmdb-crystal.hh"

// 20220723-PE maybe just delete this header altogether.
#include "globjects.h" //includes gtk/gtk.h


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
	 g.show_select_map_frame();
	 g.in_terminal_residue_define = 0;
	 g.model_fit_refine_unactive_togglebutton("model_refine_dialog_fit_terminal_residue_togglebutton");
	 g.normal_cursor();
      }
   } else {
      g.normal_cursor();
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("do-add-terminal-residue");
   command_strings.push_back(graphics_info_t::int_to_string(state));
   add_to_history(command_strings);
}

void
set_add_terminal_residue_n_phi_psi_trials(int n) {
   graphics_info_t g;
   g.add_terminal_residue_n_phi_psi_trials = n;
   std::vector<std::string> command_strings;
   command_strings.push_back("set-add-terminal-residue-n-phi-psi-trials");
   command_strings.push_back(graphics_info_t::int_to_string(n));
   add_to_history(command_strings);
}

void
set_add_terminal_residue_add_other_residue_flag(int i) {
   graphics_info_t::add_terminal_residue_add_other_residue_flag = i;
   std::vector<std::string> command_strings;
   command_strings.push_back("set-add-terminal-residue-add-other-residue-flag");
   command_strings.push_back(graphics_info_t::int_to_string(i));
   add_to_history(command_strings);
}

void set_add_terminal_residue_do_rigid_body_refine(short int v) { 

   graphics_info_t g;
   g.add_terminal_residue_do_rigid_body_refine = v;
   std::vector<std::string> command_strings;
   command_strings.push_back("set-terminal-residue-do-rigid-body-refine");
   command_strings.push_back(graphics_info_t::int_to_string(v));
   add_to_history(command_strings);

}

// deprecate this at some stage (its 201709025) now that we have the canonically
// named function above.
//
void set_terminal_residue_do_rigid_body_refine(short int v) { 
   set_add_terminal_residue_do_rigid_body_refine(v);
}

void set_add_terminal_residue_debug_trials(short int debug_state) {

   graphics_info_t g;
   g.add_terminal_residue_debug_trials = debug_state;
   std::vector<std::string> command_strings;
   command_strings.push_back("set-terminal-residue-debug-trials");
   command_strings.push_back(graphics_info_t::int_to_string(debug_state));
   add_to_history(command_strings);

}


void set_add_terminal_residue_do_post_refine(short int istat) {
   graphics_info_t::add_terminal_residue_do_post_refine = istat;
   std::vector<std::string> command_strings;
   command_strings.push_back("set-add-terminal-residue-do-post-refine");
   command_strings.push_back(graphics_info_t::int_to_string(istat));
   add_to_history(command_strings);
}

int add_terminal_residue_do_post_refine_state() {
  int i = graphics_info_t::add_terminal_residue_do_post_refine;
  return i;
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
			 const char *chain_id, 
			 int residue_number,
			 const char *residue_type,
			 int immediate_add) {

   int istate = 0;
   graphics_info_t g;
   std::string residue_type_string = residue_type;
   
   int imol_map = g.Imol_Refinement_Map();
   if (imol_map == -1) {
      std::cout << "WARNING:: Refinement/Fitting map is not set." << std::endl;
      std::cout << "          addition of terminal residue terminated." << std::endl;
   } else { 
      if (is_valid_model_molecule(imol)) { 
	 // We don't do this as a member function of
	 // molecule_class_info_t because we are using
	 // graphics_info_t::execute_add_terminal_residue, which
	 // does the molecule manipulations outside of the
	 // molecule_class_info_t class.
	 //

	 int atom_indx = atom_index(imol, chain_id, residue_number, " CA ");
	 if (atom_indx >= 0) {
	    std::string term_type = g.molecules[imol].get_term_type(atom_indx);
	    std::string inscode = "";
	    mmdb::Residue *res_p =
	       g.molecules[imol].get_residue(chain_id, residue_number, inscode);

	    if (res_p)
	       istate = g.execute_add_terminal_residue(imol, term_type, res_p, chain_id,
						       residue_type_string, immediate_add);
	    
	 } else {
	    std::cout << "WARNING:: in add_terminal_residue: "
		      << " Can't find atom index for CA in residue "
		      << residue_number << " " << chain_id << std::endl;
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

/*! \brief Add a terminal residue using given phi and psi angles
 */
int add_terminal_residue_using_phi_psi(int imol, const char *chain_id, int res_no, 
				       const char *residue_type, float phi, float psi) {

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      status = graphics_info_t::molecules[imol].add_terminal_residue_using_phi_psi(chain_id, res_no, residue_type, phi, psi);
      if (status)
	 graphics_draw();
   } 
   return status;

} 


void set_add_terminal_residue_default_residue_type(const char *type) {

   if (type) 
      graphics_info_t::add_terminal_residue_type = type;
   std::string cmd = "set-add-terminal-residue-default-residue-type";
   std::vector<coot::command_arg_t> args;
   args.push_back(single_quote(type));
   add_to_history_typed(cmd, args);
}



// draw the baton?
int try_set_draw_baton(short int i) {

   graphics_info_t g;

   g.try_set_draw_baton(i);
   std::string cmd = "set-draw-baton";
   std::vector<coot::command_arg_t> args;
   args.push_back(i);
   add_to_history_typed(cmd, args);

   return g.draw_baton_flag;
}

// Mouse movement moves the baton not the view?
void set_baton_mode(short int i) {
   graphics_info_t::baton_mode = i;
   std::string cmd = "set-baton-mode";
   std::vector<coot::command_arg_t> args;
   args.push_back(i);
   add_to_history_typed(cmd, args);
}

void accept_baton_position() { 	/* put an atom at the tip and move baton */

   graphics_info_t g;
   g.accept_baton_position();
   add_to_history_simple("accept-baton-position");
}

void baton_tip_try_another() {
   graphics_info_t g;
   g.baton_tip_try_another();
   add_to_history_simple("baton-try-another");
}

void baton_tip_previous() {
   graphics_info_t g;
   g.baton_tip_previous();
   add_to_history_simple("baton-try-another");
}

void shorten_baton() {
   graphics_info_t g;
   g.shorten_baton();
   add_to_history_simple("shorten-baton");
}

void lengthen_baton() {
   graphics_info_t g;
   g.lengthen_baton();
   add_to_history_simple("lengthen-baton");
}

void baton_build_delete_last_residue() {

   graphics_info_t g;
   g.baton_build_delete_last_residue();
   add_to_history_simple("baton-build-delete-last-residue");
} 

void set_baton_build_params(int istart_resno, 
			   const char *chain_id, 
			   const char *backwards) { 

   graphics_info_t g;
   g.set_baton_build_params(istart_resno, chain_id, backwards); 
   std::string cmd = "set-baton-build-params";
   std::vector<coot::command_arg_t> args;
   args.push_back(istart_resno);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(coot::util::single_quote(backwards));
   add_to_history_typed(cmd, args);
}



/* Reverse the direction of a the fragment of the clicked on
   atom/residue.  A fragment is a consequitive range of residues -
   where there is a gap in the numbering, that marks breaks between
   fragments in a chain.  There also needs to be a distance break - if
   the CA of the next/previous residue is more than 5A away, that also
   marks a break. Thow away all atoms in fragment other than CAs */
void reverse_direction_of_fragment(int imol, const char *chain_id, int resno) {

   if (is_valid_model_molecule(imol)) {
      // return 1 if we did it.
      int istatus = graphics_info_t::molecules[imol].reverse_direction_of_fragment(std::string(chain_id), resno);
      if (istatus)
	 graphics_draw();
   }
   std::string cmd = "reverse-direction-of-fragment";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(resno);
   add_to_history_typed(cmd, args);
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
   std::string cmd = "do-db-main";
   std::vector<coot::command_arg_t> args;
   args.push_back(state);
   add_to_history_typed(cmd, args);
}


int
db_mainchain(int imol,
	     const char *chain_id,
	     int iresno_start,
	     int iresno_end,
	     const char *direction_string) {

   int imol_new = -1;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      imol_new = g.execute_db_main(imol, std::string(chain_id), iresno_start, iresno_end,
				   std::string(direction_string));
   } else {
      std::cout << "WARNING:: molecule index error" << std::endl;
   } 
   std::string cmd = "db-mainchain";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(iresno_start);
   args.push_back(iresno_end);
   args.push_back(coot::util::single_quote(direction_string));
   add_to_history_typed(cmd, args);

   return imol_new;
}

/*! \brief CA-Zone to Mainchain for a fragment based on the given residue.

Both directions are built.
 */
int db_mainchains_fragment(int imol, const char *chain_id, int res_no) {

   int imol_new = -1;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      coot::residue_spec_t spec(chain_id, res_no, "");
      std::pair<int, int> n = g.execute_db_main_fragment(imol, spec); // builds both
      imol_new = n.first; // because we return an int not a SCM/PyObject *.
   }
   return imol_new;
}


/*  ----------------------------------------------------------------------- */
/*                  pepflip                                                 */
/*  ----------------------------------------------------------------------- */
// use the values that are in graphics_info
void do_pepflip(short int state) {

   graphics_info_t g;

   g.set_in_pepflip_define(state);
   if (state) { 
      g.pick_cursor_maybe();
      g.pick_pending_flag = 1;
      std::cout << "click on a atom in the peptide you wish to flip: "
		<< std::endl;
   } else {
      g.normal_cursor();
   } 
      
} 

void pepflip(int imol, const char *chain_id, int resno,
	     const char *inscode, const char *alt_conf) {
   /* the residue with CO,
      for scripting interface. */

   if (is_valid_model_molecule(imol)) { 
      graphics_info_t g;
      g.molecules[imol].pepflip_residue(chain_id, resno, inscode, alt_conf);
      g.update_validation(imol);
      graphics_draw();
   } 
} 

int pepflip_intermediate_atoms() {

   graphics_info_t g;
   return g.pepflip_intermediate_atoms();
} 

int pepflip_intermediate_atoms_other_peptide() {

   graphics_info_t g;
   return g.pepflip_intermediate_atoms_other_peptide();
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


/*! \brief make the planar peptide restraints tight

Useful when refining models with cryo-EM maps */
void make_tight_planar_peptide_restraints() {
   graphics_info_t g;
   g.Geom_p()->make_tight_planar_peptide_restraint();
}


/* return 1 if planar peptide restraints are on, 0 if off */
int planar_peptide_restraints_state() {

   graphics_info_t g;
   bool r = g.Geom_p()->planar_peptide_restraint_state();
   int rr = r;
   return rr;
}

/*  ----------------------------------------------------------------------- */
/*                         Trans Peptide Restraints                         */
/*  ----------------------------------------------------------------------- */


/*! \brief add a restraint on peptides to keep trans peptides trans 

i.e. omega in trans-peptides is restraints to 180 degrees.
 */
void set_use_trans_peptide_restraints(short int on_off_state) {

   graphics_info_t g;
   g.do_trans_peptide_restraints = on_off_state;
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
/*               Use Cowtan's protein_db to discover loops                  */
/*  ----------------------------------------------------------------------- */

#ifdef USE_GUILE
SCM
protein_db_loops_scm(int imol_coords, SCM residue_specs_scm, int imol_map, int nfrags,
		     bool preserve_residue_names) {

   SCM r = SCM_BOOL_F;
   std::vector<coot::residue_spec_t> specs = scm_to_residue_specs(residue_specs_scm);
   if (!specs.size()) {
      std::cout << "WARNING:: Ooops - no specs in " 
		<< scm_to_locale_string(display_scm(residue_specs_scm))
		<< std::endl;
   } else {
      std::pair<std::pair<int, int>, std::vector<int> > p = 
	 protein_db_loops(imol_coords, specs, imol_map, nfrags, preserve_residue_names);
      SCM mol_list_scm = SCM_EOL;
      // add backwards (for scheme)

      std::cout << "debug:: protein_db_loops() here are the results from protein_db_loops()" << std::endl;
      for (std::size_t ii=0; ii<p.second.size(); ii++) {
	 std::cout << "     " << p.second[ii] << std::endl;
      }

      for (int i=(p.second.size()-1); i>=0; i--)
	 mol_list_scm = scm_cons(scm_from_int(p.second[i]), mol_list_scm);
      SCM first_pair_scm = scm_list_2(scm_from_int(p.first.first), scm_from_int(p.first.second));
      r = scm_list_2(first_pair_scm, mol_list_scm);
   } 
   return r;
} 
#endif 

#ifdef USE_PYTHON
PyObject *
protein_db_loops_py(int imol_coords, PyObject *residue_specs_py, int imol_map, int nfrags,
		    bool preserve_residue_names) {

   PyObject *r = Py_False;
   std::vector<coot::residue_spec_t> specs = py_to_residue_specs(residue_specs_py);
   if (specs.empty()) {
      std::cout << "WARNING:: protein_db_loops_py(): Ooops - no specs in "
		<< PyUnicode_AsUTF8String(display_python(residue_specs_py))
		<< std::endl;
   } else {
      std::pair<std::pair<int, int>, std::vector<int> > p = protein_db_loops(imol_coords, specs, imol_map, nfrags, preserve_residue_names);
      PyObject *mol_list_py = PyList_New(p.second.size());
      for (unsigned int i=0; i< p.second.size(); i++)
         PyList_SetItem(mol_list_py, i, PyLong_FromLong(p.second[i]));
      r = PyList_New(2);
      PyObject *loop_molecules = PyList_New(2);
      PyList_SetItem(loop_molecules, 0, PyLong_FromLong(p.first.first));
      PyList_SetItem(loop_molecules, 1, PyLong_FromLong(p.first.second));

      PyList_SetItem(r, 0, loop_molecules);
      PyList_SetItem(r, 1, mol_list_py);
   } 
   if (PyBool_Check(r)) {
      Py_INCREF(r);
   }
   return r;
} 
#endif //PYTHON 

// return in the first pair, the imol of the new molecule generated
// from an atom selection of the imol_coords for the residue selection
// of the loop and the molecule number of the consolidated solutions
// (displayed in purple).  and the second of the outer pair, there is
// vector of molecule indices for each of the candidate loops.
// 
// return -1 in the first of the pair on failure
// 
std::pair<std::pair<int, int> , std::vector<int> > 
protein_db_loops(int imol_coords, const std::vector<coot::residue_spec_t> &residue_specs, int imol_map,
		 int nfrags, bool preserve_residue_names) {
   
   int imol_consolodated = -1;
   int imol_loop_orig = -1; // set later hopefully
   std::vector<int> vec_chain_mols;

   if (! is_valid_model_molecule(imol_coords)) {
      std::cout << "WARNING:: molecule number " << imol_coords 
		<< " is not a valid molecule " << std::endl;
   } else { 
      if (! is_valid_map_molecule(imol_map)) {
	 std::cout << "WARNING:: molecule number " << imol_coords << " is not a valid map " 
		   << std::endl;
      } else { 
	 
	 if (residue_specs.size()) { 
	    std::string chain_id = residue_specs[0].chain_id;
	    // what is the first resno?
	    std::vector<coot::residue_spec_t> rs = residue_specs;
	    std::sort(rs.begin(), rs.end());
	    int first_res_no = rs[0].res_no;

	    const clipper::Xmap<float> &xmap = graphics_info_t::molecules[imol_map].xmap;

	    std::vector<ProteinDB::Chain> chains =
	       graphics_info_t::molecules[imol_coords].protein_db_loops(residue_specs, nfrags, xmap);

	    graphics_info_t g;
	    if (! chains.empty()) {

	       // a molecule for each chain
	       for(unsigned int ich=0; ich<chains.size(); ich++) { 
		  mmdb::Manager *mol = make_mol(chains[ich], chain_id, first_res_no, preserve_residue_names);
		  coot::util::delete_anomalous_atoms(mol); //CBs in GLY etc
		  int imol = graphics_info_t::create_molecule();
		  std::string name = "Loop candidate #"; 
		  name += coot::util::int_to_string(ich);
		  std::cout << "INFO:: installing molecule number " << imol
			    << " with name " << name << std::endl;
		  g.molecules[imol].install_model(imol, mol, g.Geom_p(), name, 1);
		  vec_chain_mols.push_back(imol);
		  set_mol_displayed(imol, 0);
		  set_mol_active(imol, 0);
	       }

	       // The consolodated molecule
	       imol_consolodated = graphics_info_t::create_molecule();
	       mmdb::Manager *mol = make_mol(chains, chain_id, first_res_no, preserve_residue_names);
	       std::string name = "All Loop candidates "; 
	       graphics_info_t::molecules[imol_consolodated].install_model(imol_consolodated, 
									   mol, g.Geom_p(), name, 1);
	       graphics_info_t::molecules[imol_consolodated].set_bond_thickness(2);
	       graphics_info_t::molecules[imol_consolodated].bonds_colour_map_rotation = 260; // purple

	       // now create a new molecule that is the loop with a
	       // copy of the original coordinates
	       // 
	       std::string ass = protein_db_loop_specs_to_atom_selection_string(residue_specs);

	       if (true) {
		  imol_loop_orig = new_molecule_by_atom_selection(imol_coords, ass.c_str());
		  set_mol_active(imol_loop_orig, 0);
		  set_mol_displayed(imol_loop_orig, 0);
	       }
	       graphics_draw();
	    }
	 }
      }
   }
   std::pair<int, int> p(imol_loop_orig, imol_consolodated);
   return std::pair<std::pair<int, int>, std::vector<int> > (p, vec_chain_mols); 
} 

// so that we can create a "original loop" molecule from the atom
// specs picked (i.e. the atom selection string should extend over the
// range from the smallest residue number to the largest (in the same
// chain)).
std::string
protein_db_loop_specs_to_atom_selection_string(const std::vector<coot::residue_spec_t> &specs) { 

   std::string r = "////"; // fail

   // check that we have the same chain id in all the specs, if yes,
   // then continue.
   std::map<std::string, int> chain_ids;
   for (unsigned int i=0; i<specs.size(); i++)
      chain_ids[specs[i].chain_id]++;
   if (chain_ids.size() == 1) { 
      std::map<std::string, int>::const_iterator it = chain_ids.begin();
      std::string chain_id = it->first;
      int lowest_resno = 9999;
      int highest_resno = -999;
      for (unsigned int i=0; i<specs.size(); i++) { 
	 if (specs[i].res_no < lowest_resno)
	    lowest_resno = specs[i].res_no;
	 if (specs[i].res_no > highest_resno)
	    highest_resno = specs[i].res_no;
      }
      r = "//";
      r += chain_id;
      r += "/";
      r += coot::util::int_to_string(lowest_resno);
      r += "-";
      r += coot::util::int_to_string(highest_resno);
   }
   return r;
}


/* ------------------------------------------------------------------------- */
/*                      LINKs                                                */
/* ------------------------------------------------------------------------- */
void
make_link(int imol, const coot::atom_spec_t &spec_1, const coot::atom_spec_t &spec_2,
	  const std::string &link_name, float length) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      g.molecules[imol].make_link(spec_1, spec_2, link_name, length, *g.Geom_p());
      graphics_draw();
   }
}

#ifdef USE_GUILE
void make_link_scm(int imol, SCM spec_1, SCM spec_2,
		   const std::string &link_name, float length) {

   // the link name and length are currently dummy values - not used.

   coot::atom_spec_t s1 = atom_spec_from_scm_expression(spec_1);
   coot::atom_spec_t s2 = atom_spec_from_scm_expression(spec_2);
   if (s1.string_user_data != "OK")
      std::cout << "WARNING:: problem with atom spec "
		<< scm_to_locale_string(display_scm(spec_1)) << std::endl;
   else
      if (s2.string_user_data != "OK")
	 std::cout << "WARNING:: make_link_scm(): problem with atom spec "
		   << scm_to_locale_string(display_scm(spec_2)) << std::endl;
      else
	 make_link(imol, s1, s2, link_name, length);
}
#endif

#ifdef USE_PYTHON
void make_link_py(int imol, PyObject *spec_1, PyObject *spec_2,
                  const std::string &link_name, float length) {

   coot::atom_spec_t s1 = atom_spec_from_python_expression(spec_1);
   coot::atom_spec_t s2 = atom_spec_from_python_expression(spec_2);
   if (s1.string_user_data != "OK")
     std::cout << "WARNING:: make_link_py(): A problem with atom spec "
               << PyUnicode_AsUTF8String(display_python(spec_1)) << " " << s1 << std::endl;
   else
     if (s2.string_user_data != "OK")
       std::cout << "WARNING:: make_link_py() B problem with atom spec "
                 << PyUnicode_AsUTF8String(display_python(spec_2)) << " " << s2 << std::endl;
     else
       make_link(imol, s1, s2, link_name, length);
}
#endif

int add_nucleotide(int imol, const char *chain_id, int res_no) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      g.execute_simple_nucleotide_addition(imol, chain_id, res_no);
      graphics_draw();
      return 1;
   }
   return 0;

}

#include "coot-utils/pepflip-using-difference-map.hh"

#ifdef USE_GUILE
SCM pepflip_using_difference_map_scm(int imol_coords, int imol_difference_map, float n_sigma) {

   SCM r = SCM_EOL;

   if (is_valid_model_molecule(imol_coords)) {
      if (is_valid_map_molecule(imol_difference_map)) {
         graphics_info_t g;
         if (g.molecules[imol_difference_map].is_difference_map_p()) {
            const clipper::Xmap<float> &diff_xmap = g.molecules[imol_difference_map].xmap;
            mmdb::Manager *mol = g.molecules[imol_coords].atom_sel.mol;
            coot::pepflip_using_difference_map pf(mol, diff_xmap);
            std::vector<coot::residue_spec_t> flips = pf.get_suggested_flips(n_sigma);
            for (std::size_t i=0; i<flips.size(); i++) {
               SCM flip_scm = residue_spec_to_scm(flips[i]);
               r = scm_cons(flip_scm, r);
	    }
	 }
      }
   }
   r = scm_reverse(r);

   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *pepflip_using_difference_map_py(int imol_coords, int imol_difference_map, float n_sigma) {

   PyObject *o = PyList_New(0);

   if (is_valid_model_molecule(imol_coords)) {
      if (is_valid_map_molecule(imol_difference_map)) {
         graphics_info_t g;
         if (g.molecules[imol_difference_map].is_difference_map_p()) {
            const clipper::Xmap<float> &diff_xmap = g.molecules[imol_difference_map].xmap;
            mmdb::Manager *mol = g.molecules[imol_coords].atom_sel.mol;
            coot::pepflip_using_difference_map pf(mol, diff_xmap);
            std::vector<coot::residue_spec_t> flips = pf.get_suggested_flips(n_sigma);
            if (flips.size() > 0) {
               o = PyList_New(flips.size());
               for (std::size_t i=0; i<flips.size(); i++) {
                  PyList_SetItem(o, i, residue_spec_to_py(flips[i]));
               }
            }
	 }
      }
   }
   return o;
}
#endif
