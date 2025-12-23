/* src/c-interface-build.cc
 *
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007, 2008 The University of York
 * Author: Paul Emsley
 * Copyright 2007 by Paul Emsley
 * Copyright 2007 by Bernhard Lohkamp
 * Copyright 2008 by Kevin Cowtan
 * Copyright 2007, 2008, 2009, 2010, 2011 The University of Oxford
 * Copyright 2013, 2014, 2015, 2016 by Medical Research Council
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

#include <cstddef>
#include "geometry/residue-and-atom-specs.hh"
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

#define HAVE_CIF  // will become unnecessary at some stage.

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

#include "utils/coot-fasta.hh"
#include "utils/coot-utils.hh"  // for is_member_p
#include "coords/mmdb-extras.hh"
#include "coords/mmdb.hh"
#include "coords/mmdb-crystal.hh"
#include "coords/Cartesian.hh"
#include "coords/Bond_lines.hh"
#include "coot-utils/coot-map-heavy.hh"  // for fffear
#include "coot-utils/coot-coord-utils.hh"

// 20220723-PE perhaps delete (the use of) this include file completely?
#include "globjects.h" //includes gtk/gtk.h


#include "graphics-info.h"

#include "widget-headers.hh"

#include "skeleton/BuildCas.h"
#include "ligand/helix-placement.hh"
#include "ligand/fast-ss-search.hh"

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

#include "utils/logging.hh"
extern logging logger;


/*  ------------------------------------------------------------------------ */
/*                   model/fit/refine functions:                             */
/*  ------------------------------------------------------------------------ */
 void set_model_fit_refine_rotate_translate_zone_label(const char *txt) {
   graphics_info_t::model_fit_refine_rotate_translate_zone_string = txt;
   // if we have the dialog open we shall change the label
   std::vector<std::string> command_strings;
   command_strings.push_back("set-model-fit-refine-rotate-translate-zone-label");
   command_strings.push_back(txt);
   add_to_history(command_strings);
}

void set_model_fit_refine_place_atom_at_pointer_label(const char *txt) {
   graphics_info_t::model_fit_refine_place_atom_at_pointer_string = txt;
   std::vector<std::string> command_strings;
   command_strings.push_back("set-model-fit-refine-atom-at-pointer-label");
   command_strings.push_back(txt);
   add_to_history(command_strings);
}

int copy_molecule(int imol) {
   int iret = -1;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      iret = g.copy_model_molecule(imol);
      if (is_valid_model_molecule(iret))
         g.molecules[iret].set_have_unsaved_changes_from_outside();
   }
   if (is_valid_map_molecule(imol)) {
      int new_mol_number = graphics_info_t::create_molecule();
      std::string label = "Copy_of_";
      label += graphics_info_t::molecules[imol].name_; // use get_name()
      bool is_em_flag = graphics_info_t::molecules[imol].is_EM_map();
      graphics_info_t::molecules[new_mol_number].install_new_map(graphics_info_t::molecules[imol].xmap, label, is_em_flag);
      if (graphics_info_t::molecules[imol].is_difference_map_p()) {
         graphics_info_t::molecules[new_mol_number].set_map_is_difference_map(true);
      }
      iret = new_mol_number;
   }
   if (iret != -1)
      graphics_draw();
   std::vector<std::string> command_strings;
   command_strings.push_back("copy-molecule");
   command_strings.push_back(graphics_info_t::int_to_string(imol));
   add_to_history(command_strings);
   return iret;
}

int
add_ligand_delete_residue_copy_molecule(int imol_ligand_new,
					const char *chain_id_ligand_new,
					int res_no_ligand_new,
					int imol_current,
					const char *chain_id_ligand_current,
					int res_no_ligand_current) {

   int r = -1;
   bool created_flag = 0; // we only want to do this once

//    std::cout << "debug:: searching for residue :"
// 	     << chain_id_ligand_current << ": " << res_no_ligand_current
// 	     << " in molecule " << imol_current
// 	     << " replacing it with atom of :" << chain_id_ligand_new << ": "
// 	     << res_no_ligand_new << " of molecule " << imol_ligand_new
// 	     << std::endl;

   if (! is_valid_model_molecule(imol_ligand_new)) {
      std::cout << "WARNING:: ligand molecule " << imol_ligand_new << " is not a valid molecule"
		<< std::endl;
   } else {
      if (! is_valid_model_molecule(imol_current)) {
	 std::cout << "WARNING:: (surrounding) molecule " << imol_current
		   << " is not a valid molecule" << std::endl;
      } else {
	 graphics_info_t g;
	 mmdb::Residue *res_ligand_new =
	    g.molecules[imol_ligand_new].get_residue(chain_id_ligand_new,
						     res_no_ligand_new, "");
	 mmdb::Residue *res_ligand_current =
	    g.molecules[imol_current].get_residue(chain_id_ligand_current,
						  res_no_ligand_current, "");
	 if (!res_ligand_current || !res_ligand_new) {

	    // so which was it then?
	    if (! res_ligand_current)
	       std::cout << "WARNING:: Oops, reference residue (being replaced) not found"
			 << std::endl;
	    if (! res_ligand_new)
	       std::cout << "WARNING:: Oops, new residue (replacing other) not found"
			 << std::endl;

	 } else {
	    mmdb::Manager *n = new mmdb::Manager;
	    n->Copy(g.molecules[imol_current].atom_sel.mol, mmdb::MMDBFCM_All);

	    // now find the residue in imol_ligand_current.
	    // and replace its atoms.
	    int imodel = 1;
	    mmdb::Model *model_p = n->GetModel(imodel);
	    mmdb::Chain *chain_p;
	    int nchains = model_p->GetNumberOfChains();
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain_p = model_p->GetChain(ichain);
	       if (! strncmp(chain_id_ligand_current, chain_p->GetChainID(), 4)) {
		  int nres = chain_p->GetNumberOfResidues();
		  mmdb::Residue *residue_p;
		  mmdb::Atom *at;
		  for (int ires=0; ires<nres; ires++) {
		     residue_p = chain_p->GetResidue(ires);
		     if (residue_p->GetSeqNum() == res_no_ligand_current) {

			// delete the current atoms (backwards so that
			// we don't have reindexing problems)
			//
			int n_atoms = residue_p->GetNumberOfAtoms();
			for (int iat=n_atoms-1; iat>=0; iat--) {
			   residue_p->DeleteAtom(iat);
			}

			n_atoms = res_ligand_new->GetNumberOfAtoms();
			for (int iat=0; iat<n_atoms; iat++) {
			   mmdb::Atom *at_copy = new mmdb::Atom;
			   at_copy->Copy(res_ligand_new->GetAtom(iat));
			   residue_p->AddAtom(at_copy);
			}
			residue_p->SetResName(res_ligand_new->GetResName());
			n->FinishStructEdit();

			r = graphics_info_t::create_molecule();
			atom_selection_container_t asc = make_asc(n);
			std::string label = "Copy_of_";
			label += coot::util::int_to_string(imol_current);
			label += "_with_";
			label += chain_id_ligand_current;
			label += coot::util::int_to_string(res_no_ligand_current);
			label += "_replaced";
			g.molecules[r].install_model(r, asc, g.Geom_p(), label, 1);
			created_flag = 1;
			break;
		     }

		     if (created_flag)
			break;
		  }
	       }
	       if (created_flag)
		  break;
	    }
	 }
      }
   }

   if (created_flag) {
      graphics_draw();
   }
   std::cout << "add_ligand_delete_residue_copy_molecule() returns " << r << std::endl;
   return r;
}


/*! \brief replace the parts of molecule number imol that are
  duplicated in molecule number imol_frag */
int replace_fragment(int imol_target, int imol_fragment,
		     const char *mmdb_atom_selection_str) {

   int istate = 0;
   if (is_valid_model_molecule(imol_target)) {
      if (is_valid_model_molecule(imol_fragment)) {
	 mmdb::Manager *mol = graphics_info_t::molecules[imol_fragment].atom_sel.mol;

         std::vector<std::string> parts = coot::util::split_string(mmdb_atom_selection_str, "||");
	 int SelHnd = mol->NewSelection();
         for (const auto &part : parts)
            mol->Select(SelHnd, mmdb::STYPE_ATOM, part.c_str(), mmdb::SKEY_OR);
	 mmdb::Manager *mol_new = coot::util::create_mmdbmanager_from_atom_selection(mol, SelHnd);
	 atom_selection_container_t asc = make_asc(mol_new);
	 istate = graphics_info_t::molecules[imol_target].replace_fragment(asc);
	 mol->DeleteSelection(SelHnd);
	 graphics_draw();
      }
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("replace-fragement");
   command_strings.push_back(graphics_info_t::int_to_string(imol_target));
   command_strings.push_back(graphics_info_t::int_to_string(imol_fragment));
   command_strings.push_back(single_quote(mmdb_atom_selection_str));
   add_to_history(command_strings);
   return istate;
}

#ifdef USE_GUILE
int replace_residues_from_mol_scm(int imol_target,
				 int imol_ref,
				 SCM residue_specs_list_ref_scm) {

   int status = 0;
   if (is_valid_model_molecule(imol_target)) {
      if (is_valid_model_molecule(imol_ref)) {
	 mmdb::Manager *mol = graphics_info_t::molecules[imol_ref].atom_sel.mol;
	 std::vector<coot::residue_spec_t> specs = scm_to_residue_specs(residue_specs_list_ref_scm);
	 if (specs.size()) {
	    mmdb::Manager *mol_new = coot::util::create_mmdbmanager_from_residue_specs(specs, mol);
	    if (mol_new) {
	       atom_selection_container_t asc = make_asc(mol_new);
	       status = graphics_info_t::molecules[imol_target].replace_fragment(asc);
	       graphics_draw();
	    }
	 }
      }
   }

   if (false) {
      // we can't yet add the residue_spec_list_scm
      std::vector<std::string> command_strings;
      command_strings.push_back("replace-fragement-from-mol-scm");
      command_strings.push_back(graphics_info_t::int_to_string(imol_target));
      add_to_history(command_strings);
   }
   return status;
}
#endif // USE_GUILE

/*! \brief replace the given residues from the reference molecule to the target molecule
*/
#ifdef USE_PYTHON
int replace_residues_from_mol_py(int imol_target,
				 int imol_ref,
				 PyObject *residue_specs_list_ref_py) {

   int status = 0;
   if (is_valid_model_molecule(imol_target)) {
      if (is_valid_model_molecule(imol_ref)) {
	 mmdb::Manager *mol = graphics_info_t::molecules[imol_ref].atom_sel.mol;
	 std::vector<coot::residue_spec_t> specs = py_to_residue_specs(residue_specs_list_ref_py);
	 if (specs.size()) {
	    mmdb::Manager *mol_new = coot::util::create_mmdbmanager_from_residue_specs(specs, mol);
	    atom_selection_container_t asc = make_asc(mol_new);
	    status = graphics_info_t::molecules[imol_target].replace_fragment(asc);
	    graphics_draw();
	 }
      }
   }

   if (false) {
      // we can't yet add the residue_spec_list_py
      std::vector<std::string> command_strings;
      command_strings.push_back("replace-fragement-from-mol-py");
      command_strings.push_back(graphics_info_t::int_to_string(imol_target));
      add_to_history(command_strings);
   }
   return status;
}
#endif /* USE_PYTHON */


/*! \brief copy the given residue range from the reference chain to the target chain

resno_range_start and resno_range_end are inclusive. */
int copy_residue_range(int imol_target,    const char *chain_id_target,
		       int imol_reference, const char *chain_id_reference,
		       int resno_range_start, int resno_range_end) {

   int status = 0;
   if (! (is_valid_model_molecule(imol_target))) {
      std::cout << "WARNING:: not a valid model molecule "
		<< imol_target << std::endl;
   } else {
      if (! (is_valid_model_molecule(imol_reference))) {
	 std::cout << "WARNING:: not a valid model molecule "
		   << imol_reference << std::endl;
      } else {
	 mmdb::Chain *chain_p = graphics_info_t::molecules[imol_reference].get_chain(chain_id_reference);
	 if (! chain_p) {
	    std::cout << "WARNING:: not chain " << chain_id_reference << " in molecule "
		      << imol_reference << std::endl;
	 } else {
	    mmdb::Chain *chain_pt = graphics_info_t::molecules[imol_target].get_chain(chain_id_target);
	    if (! chain_pt) {
	       std::cout << "WARNING:: not chain " << chain_id_target << " in molecule "
			 << imol_target << std::endl;
	    } else {
	       clipper::RTop_orth rtop = clipper::RTop_orth::identity();
	       status = graphics_info_t::molecules[imol_target].copy_residue_range(chain_p, chain_pt,
								  resno_range_start, resno_range_end,
								  rtop);
	       graphics_draw();
	    }
	 }
      }
   }
   return status;
}


void set_refinement_move_atoms_with_zero_occupancy(int state) {
   // convert a int to a bool.
   graphics_info_t::refinement_move_atoms_with_zero_occupancy_flag = state;
}

int refinement_move_atoms_with_zero_occupancy_state() {
   // convert a bool to an int.
   return graphics_info_t::refinement_move_atoms_with_zero_occupancy_flag;
}


/*  ------------------------------------------------------------------------ */
/*                         backup/undo functions:                            */
/*  ------------------------------------------------------------------------ */

void turn_off_backup(int imol) {

   if (is_valid_model_molecule(imol))
      graphics_info_t::molecules[imol].turn_off_backup();
   std::vector<std::string> command_strings;
   command_strings.push_back("turn-off-backup");
   command_strings.push_back(graphics_info_t::int_to_string(imol));
   add_to_history(command_strings);
}

void turn_on_backup(int imol) {
   if (is_valid_model_molecule(imol))
      graphics_info_t::molecules[imol].turn_on_backup();
   std::vector<std::string> command_strings;
   command_strings.push_back("turn-on-backup");
   command_strings.push_back(graphics_info_t::int_to_string(imol));
   add_to_history(command_strings);
}

int apply_undo() {		/* "Undo" button callback */
   graphics_info_t g;
   int r = g.apply_undo();
   add_to_history_simple("apply-undo");
   return r;
}

int  apply_redo() {
   graphics_info_t g;
   int r = g.apply_redo();
   add_to_history_simple("apply-redo");
   return r;
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

   if ((imol >= 0) && (imol < graphics_info_t::n_molecules())) {
      graphics_info_t g;
      if (g.molecules[imol].atom_sel.mol) {
	 std::cout << "INFO:: undo molecule number set to: " << imol << std::endl;
	 g.set_undo_molecule_number(imol);
      }
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("set-undo-molecule");
   command_strings.push_back(graphics_info_t::int_to_string(imol));
   add_to_history(command_strings);
}

void set_unpathed_backup_file_names(int state) {
   graphics_info_t::unpathed_backup_file_names_flag = state;
   std::vector<std::string> command_strings;
   command_strings.push_back("set-unpathed-backup-file-names");
   command_strings.push_back(graphics_info_t::int_to_string(state));
   add_to_history(command_strings);
}

int  unpathed_backup_file_names_state() {
   add_to_history_simple("unpathed-backup-file-names-state");
   return graphics_info_t::unpathed_backup_file_names_flag;
}

void set_backup_compress_files(int state) {

  graphics_info_t::backup_compress_files_flag = state;
  std::vector<std::string> command_strings;
  command_strings.push_back("set-backup-compress-files");
  command_strings.push_back(graphics_info_t::int_to_string(state));
  add_to_history(command_strings);
}

int backup_compress_files_state() {

  int state = graphics_info_t::backup_compress_files_flag;
  return state;
}

void set_decoloned_backup_file_names(int state) {
   graphics_info_t::decoloned_backup_file_names_flag = state;
   std::vector<std::string> command_strings;
   command_strings.push_back("set-decoloned-backup-file-names");
   command_strings.push_back(graphics_info_t::int_to_string(state));
   add_to_history(command_strings);
}

int decoloned_backup_file_names_state() {
   add_to_history_simple("decoloned-backup-file-names-state");
   return graphics_info_t::decoloned_backup_file_names_flag;
}


/*! \brief Make a backup for a model molecule
 *
 * @param imol the model molecule index
 * @description a description that goes along with this back point
 */
int make_backup_checkpoint(int imol, const char *description) {

   int backup_index = -1;
   if (is_valid_model_molecule(imol)) {
      std::string ss(description);
      backup_index = graphics_info_t::molecules[imol].make_backup_checkpoint(ss);
   }
   return backup_index;
}

/*! \brief Restore molecule from backup
 * 
 * restore model @p imol to checkpoint backup @p backup_index
 *
 * @param imol the model molecule index
 * @param backup_index the backup index to restore to
 */
int restore_to_backup_checkpoint(int imol, int backup_index) {

   backup_index = -1;
   if (is_valid_model_molecule(imol)) {
      backup_index = graphics_info_t::molecules[imol].restore_to_backup_checkpoint(backup_index);
   }
   return backup_index;
}

#ifdef USE_PYTHON
/*! \brief Compare current model to backup
 * 
 * @param imol the model molecule index
 * @param backup_index the backup index to restore to
 * @return a Python dict, with 2 items, a "status" which is either "ok" 
 *         or "fail" or "bad-index" and a list of residue specs for residues
 *         that have at least one atom in a different place (which might be empty).
 */
PyObject *compare_current_model_to_backup(int imol, int backup_index) {

   PyObject *d = PyDict_New();
   if (is_valid_model_molecule(imol)) {
      // How do I return "bad backup_index?" Use a pair.
      std::pair<bool, std::vector<coot::residue_spec_t> > mvp = graphics_info_t::molecules[imol].compare_current_model_to_backup(backup_index);
      if (mvp.first) {
         std::vector<coot::residue_spec_t> mv = mvp.second;
         PyObject *l = PyList_New(mv.size());
         for (unsigned int i=0; i<mv.size(); i++) {
             const auto &rs(mv[i]);
             PyObject *s = residue_spec_to_py(rs);
             PyList_SetItem(l, i, s);
         }
         PyDict_SetItemString(d, "moved-residues-list", l);
         PyDict_SetItemString(d, "status", myPyString_FromString("ok"));
      } else {
         PyDict_SetItemString(d, "status", myPyString_FromString("bad-index"));
      }
   } else {
      PyDict_SetItemString(d, "status", myPyString_FromString("fail"));
   }
   return d;
}
#endif

#ifdef USE_PYTHON
/*! \brief Get backup info
 * 
 * @param imol the model molecule index
 * @param backup_index the backup index to restore to
 * @return a Python list of the given description (str)
 *         and a timestamp (str).
 */
PyObject *get_backup_info(int imol, int backup_index) {

   PyObject *r = PyList_New(0);
   if (is_valid_model_molecule(imol)) {
      auto backup_info = graphics_info_t::molecules[imol].get_backup_info(backup_index);
      r = PyList_New(2);
      PyObject *d  = myPyString_FromString(backup_info.description.c_str());
      PyObject *dt = myPyString_FromString(backup_info.get_timespec_string().c_str());
      PyList_SetItem(r, 0, d);
      PyList_SetItem(r, 1, dt);
   }
   return r;
}
#endif

void print_backup_history_info(int imol) {
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].print_backup_history_info();
   }
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
   std::string cmd = "do-rot-trans-setup";
   std::vector<coot::command_arg_t> args;
   args.push_back(state);
   add_to_history_typed(cmd, args);
}


void rot_trans_reset_previous() {
   graphics_info_t g;
   // rot_trans adjustments:
   for (int i=0; i<6; i++)
      g.previous_rot_trans_adjustment[i] = -10000;
   add_to_history_simple("rot-trans-reset-previous");
}

void set_rotate_translate_zone_rotates_about_zone_centre(int istate) {
   graphics_info_t::rot_trans_zone_rotates_about_zone_centre = istate;
   std::string cmd = "set-rotate-translate-zone-rotates-about-zone-centre";
   std::vector<coot::command_arg_t> args;
   args.push_back(istate);
   add_to_history_typed(cmd, args);
}

void set_rot_trans_object_type(short int rt_type) { /* zone, chain, mol */

   graphics_info_t::rot_trans_object_type = rt_type;
}

int
get_rot_trans_object_type() {

   return graphics_info_t::rot_trans_object_type;
}

/*  ----------------------------------------------------------------------- */
/*                  spin search                                             */
/*  ----------------------------------------------------------------------- */
void spin_search_by_atom_vectors(int imol_map, int imol, const std::string &chain_id,
				 int resno, const std::string &ins_code,
				 const std::pair<std::string, std::string> &direction_atoms,
				 const std::vector<std::string> &moving_atoms_list) {


   if (is_valid_map_molecule(imol_map)) {

      if (is_valid_model_molecule(imol)) {

	 graphics_info_t::molecules[imol].spin_search(graphics_info_t::molecules[imol_map].xmap,
						      chain_id, resno, ins_code,
						      direction_atoms, moving_atoms_list);
	 graphics_draw();

      } else {
	 std::cout << "Molecule number " << imol << " is not a valid model" << std::endl;
      }
   } else {
      std::cout << "Molecule number " << imol_map << " is not a valid map" << std::endl;
   }
}

#ifdef USE_GUILE
/*! \brief for the given residue, spin the atoms in moving_atom_list
  around the bond defined by direction_atoms_list looking for the best
  fit to density of imom_map map of the first atom in
  moving_atom_list.  Works (only) with atoms in altconf "" */
void spin_search(int imol_map, int imol, const char *chain_id, int resno,
		 const char *ins_code, SCM direction_atoms_list, SCM moving_atoms_list) {

   std::vector<std::string> s = generic_list_to_string_vector_internal(direction_atoms_list);

   if (s.size() == 2) {
      std::pair<std::string, std::string> p(s[0], s[1]);

      spin_search_by_atom_vectors(imol_map, imol, chain_id, resno, ins_code, p,
				  generic_list_to_string_vector_internal(moving_atoms_list));
   } else {
      std::cout << "bad direction atom pair" << std::endl;
   }
}
#endif
#ifdef USE_PYTHON
void spin_search_py(int imol_map, int imol, const char *chain_id, int resno,
                 const char *ins_code, PyObject *direction_atoms_list, PyObject *moving_atoms_list) {

   std::vector<std::string> s = generic_list_to_string_vector_internal_py(direction_atoms_list);

   if (s.size() == 2) {
      std::pair<std::string, std::string> p(s[0], s[1]);

      spin_search_by_atom_vectors(imol_map, imol, chain_id, resno, ins_code, p,
                                  generic_list_to_string_vector_internal_py(moving_atoms_list));
   } else {
      std::cout << "bad direction atom pair" << std::endl;
   }
}
#endif // PYTHON

/*  ----------------------------------------------------------------------- */
/*                  Spin around N and CB for N-termal addition              */
/*  ----------------------------------------------------------------------- */

#ifdef USE_PYTHON
void spin_N_py(int imol, PyObject *residue_spec_py, float angle) {

   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t residue_spec = residue_spec_from_py(residue_spec_py);
      graphics_info_t::molecules[imol].spin_N(residue_spec, angle);
      graphics_draw();
   }
}
#endif // USE_PYTHON

#ifdef USE_GUILE
void spin_N_scm(int imol, SCM residue_spec_scm, float angle) {

   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t residue_spec = residue_spec_from_scm(residue_spec_scm);
      graphics_info_t::molecules[imol].spin_N(residue_spec, angle);
      graphics_draw();
   }
}
#endif // USE_GUILE

#ifdef USE_PYTHON
//! \brief Spin search the density based on possible positions of CG of a side-chain
PyObject *CG_spin_search_py(int imol_model, int imol_map) {

   PyObject *r = Py_False;

   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {
	 graphics_info_t g;
	 const clipper::Xmap<float> &xmap = g.molecules[imol_map].xmap;
	 std::vector<std::pair<coot::residue_spec_t, float> > rv =
	    g.molecules[imol_model].em_ringer(xmap);
	 r = PyList_New(rv.size());
	 for (std::size_t i=0; i<rv.size(); i++) {
	    const coot::residue_spec_t &spec = rv[i].first;
	    double delta_angle = rv[i].second;
	    PyObject *item_py = PyList_New(2);
	    PyList_SetItem(item_py, 0, PyFloat_FromDouble(delta_angle));
	    PyList_SetItem(item_py, 1, residue_spec_to_py(spec));
	    PyList_SetItem(r, i, item_py);
	 }
      }
   }

   if (PyBool_Check(r))
     Py_INCREF(r);

   return r;
}
#endif // USE_PYTHON

#ifdef USE_GUILE
//! \brief Spin search the density based on possible positions of CG of a side-chain
SCM CG_spin_search_scm(int imol_model, int imol_map) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {
	 graphics_info_t g;
	 const clipper::Xmap<float> &xmap = g.molecules[imol_map].xmap;
	 std::vector<std::pair<coot::residue_spec_t, float> > rv =
	    g.molecules[imol_model].em_ringer(xmap);
	 r = SCM_EOL;
	 for (std::size_t i=0; i<rv.size(); i++) {
	    const coot::residue_spec_t &spec = rv[i].first;
	    double delta_angle = rv[i].second;
	    SCM res_spec_scm = residue_spec_to_scm(spec);
	    SCM item_scm = scm_list_2(res_spec_scm, scm_from_double(delta_angle));
	    r = scm_cons(item_scm, r);
	 }
	 r = scm_reverse(r);
      }
   }
   return r;
}
#endif // USE_GUILE




/*  ----------------------------------------------------------------------- */
/*                  delete residue                                          */
/*  ----------------------------------------------------------------------- */
void delete_residue(int imol, const char *chain_id, int resno, const char *inscode) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      int model_number_ANY = mmdb::MinInt4;
      std::string ic(inscode);
      short int istat = g.molecules[imol].delete_residue(model_number_ANY, chain_id, resno, ic);

      g.update_validation(imol);

      if (istat) {
	 // now if the go to atom widget was being displayed, we need to
	 // redraw the residue list and atom list (if the molecule of the
	 // residue and atom list is the molecule that has just been
	 // deleted)

	 g.update_go_to_atom_window_on_changed_mol(imol);

	 if (! is_valid_model_molecule(imol)) {

	    g.delete_molecule_from_display_manager(imol, false);
         }

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
   } else {
      add_status_bar_text("Oops bad molecule from whcih to delete a residue");
   }
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


void
delete_residue_with_full_spec(int imol,
			      int imodel,
			      const char *chain_id,
			      int resno,
			      const char *inscode,
			      const char *altloc) {
   graphics_info_t g;
   if (is_valid_model_molecule(imol)) {
      std::string altconf(altloc);
      short int istat =
	 g.molecules[imol].delete_residue_with_full_spec(imodel, chain_id, resno, inscode, altconf);

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
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("delete-residue-with-full_spec");
   command_strings.push_back(g.int_to_string(imol));
   command_strings.push_back(g.int_to_string(imodel));
   command_strings.push_back(single_quote(chain_id));
   command_strings.push_back(g.int_to_string(resno));
   command_strings.push_back(single_quote(inscode));
   command_strings.push_back(single_quote(altloc));
   add_to_history(command_strings);
}


#ifdef USE_GUILE
/*! \brief delete residues in the residue spec list */
void delete_residues_scm(int imol, SCM residue_specs_scm) {
   if (is_valid_model_molecule(imol)) {
      std::vector<coot::residue_spec_t> specs = scm_to_residue_specs(residue_specs_scm);
      graphics_info_t::molecules[imol].delete_residues(specs);
      graphics_draw();
   }
}
#endif

#ifdef USE_PYTHON
/*! \brief delete residues in the residue spec list */
void delete_residues_py(int imol, PyObject *residue_specs_py) {
   if (is_valid_model_molecule(imol)) {
      std::vector<coot::residue_spec_t> specs = py_to_residue_specs(residue_specs_py);
      graphics_info_t::molecules[imol].delete_residues(specs);
      graphics_draw();
   }
}
#endif




/*! \brief delete all hydrogens in molecule */
int delete_hydrogens(int imol) {

   int n_deleted = 0;
   if (is_valid_model_molecule(imol)) {
      n_deleted = graphics_info_t::molecules[imol].delete_hydrogens();
      if (n_deleted)
	 graphics_draw();
   }
   return n_deleted;
}

/*! \brief delete all hydrogens in molecule */
int delete_hydrogen_atoms(int imol) {
   return delete_hydrogens(imol);
}

int delete_waters(int imol) {

   int n_deleted = 0;
   if (is_valid_model_molecule(imol)) {
      n_deleted = graphics_info_t::molecules[imol].delete_waters();
      if (n_deleted)
	 graphics_draw();
   }
   return n_deleted;
}

void delete_chain(int imol, const std::string &chain_id_in) {

   std::string chain_id(chain_id_in);
   graphics_info_t g;
   if (is_valid_model_molecule(imol)) {
      g.delete_chain_from_geometry_graphs(imol, chain_id);
      short int istat = g.molecules[imol].delete_chain(std::string(chain_id));
      if (istat) {
	 g.update_go_to_atom_window_on_changed_mol(imol);
	 graphics_draw();
      }

#if 0 // 20220609-PE we don't need to do this sort of thing any more
      if (delete_item_widget_is_being_shown()) {
	 if (delete_item_widget_keep_active_on()) {
	    // dont destroy it
	 } else {
	    //store_delete_item_widget_position(); // and destroy it.
	 }
      }
#endif
      if (! is_valid_model_molecule(imol))
	 g.delete_molecule_from_display_manager(imol, false);

      std::string cmd = "delete-chain";
      std::vector<coot::command_arg_t> args;
      args.push_back(imol);
      args.push_back(coot::util::single_quote(chain_id));
      add_to_history_typed(cmd, args);
      graphics_draw();
   }

}

/*! \brief delete the chain  */
void delete_sidechains_for_chain(int imol, const std::string &chain_id_in) {

   std::string chain_id(chain_id_in);
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;

      g.molecules[imol].delete_sidechains_for_chain(std::string(chain_id));
      std::string cmd = "delete-sidechains-for-chain";
      std::vector<coot::command_arg_t> args;
      args.push_back(imol);
      args.push_back(coot::util::single_quote(chain_id));
      add_to_history_typed(cmd, args);
      graphics_draw();
   }
}





void set_add_alt_conf_new_atoms_occupancy(float f) {  /* default 0.5 */

   graphics_info_t g;
   g.add_alt_conf_new_atoms_occupancy = f;
   std::string cmd = "set-add-alt-conf-new-atoms-occupancy";
   std::vector<coot::command_arg_t> args;
   args.push_back(f);
   add_to_history_typed(cmd, args);
}


float get_add_alt_conf_new_atoms_occupancy() {

   graphics_info_t g;
   return g.add_alt_conf_new_atoms_occupancy;
}


void set_numerical_gradients(int istate) {

   graphics_info_t::do_numerical_gradients = istate;
}

void set_debug_refinement(int state) {
   graphics_info_t::do_debug_refinement = state;
}

#ifdef USE_GUILE
SCM get_residue_alt_confs_scm(int imol, const char *chain_id, int res_no, const char *ins_code) {

   SCM r = SCM_EOL;
   std::cout << "get_residue_alt_confs_scm(): Needs to be implemented" << std::endl;
   return r;

}
#endif

#ifdef USE_PYTHON
/*! \brief Return either None (on failure) or a list of alt-conf strings (might be [""]) */
PyObject *get_residue_alt_confs_py(int imol, const char *chain_id, int res_no, const char *ins_code) {

   PyObject *r = Py_False;

   if (is_valid_model_molecule(imol)) {
      mmdb::Residue *residue_p = graphics_info_t::molecules[imol].get_residue(chain_id, res_no, ins_code);
      if (residue_p) {
         std::vector<std::string> ac = graphics_info_t::molecules[imol].get_residue_alt_confs(residue_p);
         if (! ac.empty()) {
            r = PyList_New(ac.size());
            for (unsigned int i=0; i<ac.size(); i++) {
               PyObject *s = myPyString_FromString(ac[i].c_str());
               PyList_SetItem(r, i, s);
            }
         }
      }
   }

   if (PyBool_Check(r))
      Py_INCREF(r);

   return r;
}
#endif


/*! \brief swap atom alt-confs */
int swap_residue_alt_confs(int imol, const char *chain_id, int res_no, const char *ins_code) {

   int state = 0;

   if (is_valid_model_molecule(imol)) {
      state = graphics_info_t::molecules[imol].swap_residue_alt_confs(chain_id, res_no, ins_code);
      graphics_info_t::graphics_draw();
   }

   std::string cmd = "swap-residue-alt-confs";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(res_no);
   args.push_back(coot::util::single_quote(ins_code));
   add_to_history_typed(cmd, args);
   return state;
}


int swap_atom_alt_conf(int imol, const char *chain_id, int res_no, const char *ins_code, const char *atom_name, const char*alt_conf) {

   int istat = 0;
   if (is_valid_model_molecule(imol)) {
      istat = graphics_info_t::molecules[imol].swap_atom_alt_conf(chain_id, res_no, ins_code, atom_name, alt_conf);
   }
   graphics_draw();
   std::string cmd = "swap-atom-alt-conf";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(res_no);
   args.push_back(coot::util::single_quote(ins_code));
   args.push_back(coot::util::single_quote(atom_name));
   args.push_back(coot::util::single_quote(alt_conf));
   add_to_history_typed(cmd, args);
   return istat;
}


int set_atom_attribute(int imol, const char *chain_id, int resno, const char *ins_code, const char *atom_name, const char*alt_conf, const char *attribute_name, float val) {
   int istat = 0;
   if (is_valid_model_molecule(imol)) {
      istat = graphics_info_t::molecules[imol].set_atom_attribute(chain_id, resno, ins_code, atom_name, alt_conf, attribute_name, val);
   }
   graphics_draw();
   std::string cmd = "set-atom-attribute";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(resno);
   args.push_back(coot::util::single_quote(ins_code));
   args.push_back(coot::util::single_quote(atom_name));
   args.push_back(coot::util::single_quote(alt_conf));
   args.push_back(coot::util::single_quote(attribute_name));
   args.push_back(val);
   add_to_history_typed(cmd, args);
   return istat;
}

int set_atom_string_attribute(int imol, const char *chain_id, int resno, const char *ins_code, const char *atom_name, const char*alt_conf, const char *attribute_name, const char *attribute_value) {
   int istat = 0;
   if (is_valid_model_molecule(imol)) {
      istat = graphics_info_t::molecules[imol].set_atom_string_attribute(chain_id, resno, ins_code, atom_name, alt_conf, attribute_name, attribute_value);
      graphics_draw();
   }
   std::string cmd = "set-atom-string-attribute";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(resno);
   args.push_back(coot::util::single_quote(ins_code));
   args.push_back(coot::util::single_quote(atom_name));
   args.push_back(coot::util::single_quote(alt_conf));
   args.push_back(coot::util::single_quote(attribute_name));
   args.push_back(coot::util::single_quote(attribute_value));
   add_to_history_typed(cmd, args);
   return istat;
}

#ifdef USE_GUILE
int set_atom_attributes(SCM attribute_expression_list) {

   int r= 0;
   SCM list_length_scm = scm_length(attribute_expression_list);
   int list_length = scm_to_int(list_length_scm);
   int n = graphics_info_t::n_molecules();
   std::vector<std::vector<coot::atom_attribute_setting_t> > v(n);

   if (list_length > 0) {
      for (int iattr=0; iattr<list_length; iattr++) {
	 SCM iattr_scm = scm_from_int(iattr);
	 SCM attribute_expression = scm_list_ref(attribute_expression_list, iattr_scm);
	 if (scm_is_true(scm_list_p(attribute_expression))) {
	    SCM attr_expression_length_scm = scm_length(attribute_expression);
	    int attr_expression_length = scm_to_int(attr_expression_length_scm);
	    if (attr_expression_length != 8) {
	       std::cout << "Incomplete attribute expression: "
			 << scm_to_locale_string(display_scm(attribute_expression))
			 << std::endl;
	    } else {
	       SCM imol_scm            = scm_list_ref(attribute_expression, scm_from_int(0));
	       SCM chain_id_scm        = scm_list_ref(attribute_expression, scm_from_int(1));
	       SCM resno_scm           = scm_list_ref(attribute_expression, scm_from_int(2));
	       SCM ins_code_scm        = scm_list_ref(attribute_expression, scm_from_int(3));
	       SCM atom_name_scm       = scm_list_ref(attribute_expression, scm_from_int(4));
	       SCM alt_conf_scm        = scm_list_ref(attribute_expression, scm_from_int(5));
	       SCM attribute_name_scm  = scm_list_ref(attribute_expression, scm_from_int(6));
	       SCM attribute_value_scm = scm_list_ref(attribute_expression, scm_from_int(7));
	       int imol = scm_to_int(imol_scm);
	       if (is_valid_model_molecule(imol)) {
		  std::string chain_id = scm_to_locale_string(chain_id_scm);
		  int resno = scm_to_int(resno_scm);

		  std::string inscode        = "-*-unset-*-:";
		  std::string atom_name      = "-*-unset-*-:";
		  std::string alt_conf       = "-*-unset-*-:";
		  std::string attribute_name = "-*-unset-*-:";

		  if (scm_is_true(scm_string_p(ins_code_scm)))
		      inscode        = scm_to_locale_string(ins_code_scm);
		  if (scm_is_true(scm_string_p(atom_name_scm)))
		     atom_name      = scm_to_locale_string(atom_name_scm);
		  if (scm_is_true(scm_string_p(alt_conf_scm)))
		     alt_conf       = scm_to_locale_string(alt_conf_scm);
		  if (scm_is_true(scm_string_p(attribute_name_scm)))
		     attribute_name = scm_to_locale_string(attribute_name_scm);

		  if ((inscode        == "-*-unset-*-:") ||
		      (atom_name      == "-*-unset-*-:") ||
		      (alt_conf       == "-*-unset-*-:") ||
		      (attribute_name == "-*-unset-*-:")) {

		     std::cout << "WARNING:: bad attribute expression: "
			       << scm_to_locale_string(display_scm(attribute_expression))
			       << std::endl;

		  } else {

		     coot::atom_attribute_setting_help_t att_val;
		     if (scm_is_true(scm_string_p(attribute_value_scm))) {
			// std::cout << "a string value :" << att_val.s << ":" << std::endl;
			att_val = coot::atom_attribute_setting_help_t(scm_to_locale_string(attribute_value_scm));
		     } else {
			att_val = coot::atom_attribute_setting_help_t(float(scm_to_double(attribute_value_scm)));
			// std::cout << "a float value :" << att_val.val << ":" << std::endl;
		     }
		     v[imol].push_back(coot::atom_attribute_setting_t(chain_id, resno, inscode, atom_name, alt_conf, attribute_name, att_val));
		     //		     std::cout << "DEBUG:: Added attribute: "
		     //                        << scm_to_locale_string(display_scm(attribute_expression))
		     //        << std::endl;
		  }
	       }
	    }
	 }
      }
   }

   for (int i=0; i<n; i++) {
      if (v[i].size() > 0) {
	 // std::cout << "DEBUG:: setting atom attributes for molecule " << i << " " << v[i].size()
	 //           << " attributes to set " << std::endl;
	 graphics_info_t::molecules[i].set_atom_attributes(v[i]);
      }
   }
   if (v.size() > 0)
      graphics_draw();
   return r;
}
#endif // USE_GUILE

#ifdef USE_PYTHON
int set_atom_attributes_py(PyObject *attribute_expression_list) {

   int r= 0;
   int list_length = PyObject_Length(attribute_expression_list);
   int n = graphics_info_t::n_molecules();
   std::vector<std::vector<coot::atom_attribute_setting_t> > v(n);
   PyObject *attribute_expression;
   PyObject *imol_py;
   PyObject *chain_id_py;
   PyObject *resno_py;
   PyObject *ins_code_py;
   PyObject *atom_name_py;
   PyObject *alt_conf_py;
   PyObject *attribute_name_py;
   PyObject *attribute_value_py;

   if (list_length > 0) {
      for (int iattr=0; iattr<list_length; iattr++) {
	 attribute_expression = PyList_GetItem(attribute_expression_list, iattr);
	 if (PyList_Check(attribute_expression)) {
	    int attr_expression_length = PyObject_Length(attribute_expression);
	    if (attr_expression_length != 8) {
               // char *ps = PyUnicode_AsUTF8String(display_python(attribute_expression));
	       char *ps = 0; // FIXME Python3
	       if (ps) {
		  std::string ae(ps);
		  std::cout << "Incomplete attribute expression: " << ae << std::endl;
	       }
	    } else {
	       imol_py            = PyList_GetItem(attribute_expression, 0);
	       chain_id_py        = PyList_GetItem(attribute_expression, 1);
	       resno_py           = PyList_GetItem(attribute_expression, 2);
	       ins_code_py        = PyList_GetItem(attribute_expression, 3);
	       atom_name_py       = PyList_GetItem(attribute_expression, 4);
	       alt_conf_py        = PyList_GetItem(attribute_expression, 5);
	       attribute_name_py  = PyList_GetItem(attribute_expression, 6);
	       attribute_value_py = PyList_GetItem(attribute_expression, 7);
	       int imol = PyLong_AsLong(imol_py);
	       if (is_valid_model_molecule(imol)) {

		  if (! PyUnicode_Check(chain_id_py)) {
		     std::cout << "WARNING:: bad chain " << chain_id_py << std::endl;
		  } else {
                     // std::string chain_id = PyUnicode_AsUTF8String(chain_id_py);
		     std::string chain_id = PyBytes_AS_STRING(PyUnicode_AsEncodedString(chain_id_py, "UTF-8", "strict"));
		     int resno = PyLong_AsLong(resno_py);

		     std::string inscode        = "-*-unset-*-:";
		     std::string atom_name      = "-*-unset-*-:";
		     std::string alt_conf       = "-*-unset-*-:";
		     std::string attribute_name = "-*-unset-*-:";

		     if (PyUnicode_Check(ins_code_py))
			inscode        = PyBytes_AS_STRING(PyUnicode_AsUTF8String(ins_code_py));
		     if (PyUnicode_Check(atom_name_py))
			atom_name      = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_name_py));
		     if (PyUnicode_Check(alt_conf_py))
			alt_conf       = PyBytes_AS_STRING(PyUnicode_AsUTF8String(alt_conf_py));
		     if (PyUnicode_Check(attribute_name_py))
			attribute_name = PyBytes_AS_STRING(PyUnicode_AsUTF8String(attribute_name_py));

		     if ((inscode        == "-*-unset-*-:") ||
			 (atom_name      == "-*-unset-*-:") ||
			 (alt_conf       == "-*-unset-*-:") ||
			 (attribute_name == "-*-unset-*-:")) {

                        std::string ss = myPyString_AsString(display_python(attribute_expression));
			std::cout << "WARNING:: bad attribute expression: "
				  << PyUnicode_AsUTF8String(attribute_expression)
				  << std::endl;

		     } else {

			coot::atom_attribute_setting_help_t att_val;
			if (PyUnicode_Check(attribute_value_py)) {
			   att_val = coot::atom_attribute_setting_help_t(myPyString_AsString(attribute_value_py));
			} else {
			   att_val = coot::atom_attribute_setting_help_t(float(PyFloat_AsDouble(attribute_value_py)));
			   // std::cout << "debug:: a float value :" << att_val.val << ":" << std::endl;
			}
                        coot::atom_attribute_setting_t as(chain_id, resno, inscode, atom_name, alt_conf,
                                                          attribute_name, att_val);
			v[imol].push_back(as);

                        if (false)
                           std::cout << "DEBUG:: Added attribute: "
                                     << myPyString_AsString(display_python(attribute_expression));
		     }
		  }
	       }
	    }
	 }
      }
   }

   for (int i=0; i<n; i++) {
      if (v[i].size() > 0){
	 graphics_info_t::molecules[i].set_atom_attributes(v[i]);
      }
   }
   if (v.size() > 0)
      graphics_draw();
   return r;
}
#endif // USE_PYTHON



void set_residue_name(int imol, const char *chain_id, int res_no, const char *ins_code, const char *new_residue_name) {

   if (chain_id && ins_code && new_residue_name) {
      if (is_valid_model_molecule(imol)) {
	 graphics_info_t::molecules[imol].set_residue_name(chain_id, res_no, ins_code, new_residue_name);
	 graphics_draw();
      }
      std::string cmd = "set-residue-name";
      std::vector<coot::command_arg_t> args;
      args.push_back(imol);
      args.push_back(coot::util::single_quote(chain_id));
      args.push_back(res_no);
      args.push_back(coot::util::single_quote(ins_code));
      args.push_back(coot::util::single_quote(new_residue_name));
      add_to_history_typed(cmd, args);
   }
}

#ifdef USE_GUILE
SCM all_residues_with_serial_numbers_scm(int imol) {

   SCM r = SCM_BOOL_F;

   if (is_valid_model_molecule(imol)) {
      std::vector<coot::residue_spec_t> specs = graphics_info_t::molecules[imol].all_residues();
      r = SCM_EOL;
      for (std::size_t i=0; i<specs.size(); i++) {
	 SCM spec_scm = residue_spec_to_scm(specs[i]);
	 int iserial = specs[i].int_user_data;
	 spec_scm = scm_cons(scm_from_int(iserial), spec_scm);
	 r = scm_cons(spec_scm, r);
      }
      r = scm_reverse(r);
   }
   return r;
}
#endif


#ifdef USE_PYTHON
PyObject *all_residues_with_serial_numbers_py(int imol) {

   PyObject *r = Py_False;

   if (is_valid_model_molecule(imol)) {
      std::vector<coot::residue_spec_t> specs = graphics_info_t::molecules[imol].all_residues();
      r = PyList_New(specs.size());
      for (std::size_t i=0; i<specs.size(); i++) {
	 PyObject *spec_py = residue_spec_to_py(specs[i]);
	 int iserial = specs[i].int_user_data;
	 PyList_Insert(spec_py, 0, PyLong_FromLong(iserial));
	 PyList_SetItem(r, i, spec_py);
      }
   }
   if (PyBool_Check(r))
     Py_INCREF(r);
   return r;
}
#endif

void
regularize_residues(int imol, const std::vector<coot::residue_spec_t> &residue_specs) {

   std::string alt_conf;
   if (is_valid_model_molecule(imol)) {
      if (residue_specs.size() > 0) {
	 std::vector<mmdb::Residue *> residues;
	 for (unsigned int i=0; i<residue_specs.size(); i++) {
	    coot::residue_spec_t rs = residue_specs[i];
	    mmdb::Residue *r = graphics_info_t::molecules[imol].get_residue(rs);
	    if (r) {
	       residues.push_back(r);
	    }
	 }

	 if (residues.size() > 0) {
	    graphics_info_t g;
	    // normal
	    mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
	    g.regularize_residues_vec(imol, residues, alt_conf.c_str(), mol);
	 }
      } else {
	 std::cout << "No residue specs found" << std::endl;
      }
   }
}


/*! \brief If there is a refinement on-going already, we don't want to start a new one

The is the means to ask if that is the case. This needs a scheme wrapper to provide refinement-already-ongoing?
  */
short int refinement_already_ongoing_p() {

   short int state = 0;
   if (graphics_info_t::moving_atoms_displayed_p())
      state = 1;
   return state;
}


#ifdef USE_GUILE
SCM refine_residues_scm(int imol, SCM r) {
   return refine_residues_with_alt_conf_scm(imol, r, "");
}
#endif // USE_GUILE

#ifdef USE_GUILE
SCM regularize_residues_scm(int imol, SCM r) { /* presumes the alt_conf is "". */
   return regularize_residues_with_alt_conf_scm(imol, r, "");
}
#endif // USE_GUILE

coot::refinement_results_t
refine_residues_with_alt_conf(int imol, const std::vector<coot::residue_spec_t> &residue_specs,
			      const std::string &alt_conf) {

   coot::refinement_results_t rr;
   if (graphics_info_t::moving_atoms_displayed_p()) {
      add_status_bar_text("No refinement - a modelling/refinement operation is already underway");
   } else {
      if (is_valid_model_molecule(imol)) {
	 if (residue_specs.size() > 0) {
	    std::vector<mmdb::Residue *> residues;
	    for (unsigned int i=0; i<residue_specs.size(); i++) {
	       coot::residue_spec_t rs = residue_specs[i];
	       mmdb::Residue *r = graphics_info_t::molecules[imol].get_residue(rs);
	       if (r) {
		  residues.push_back(r);
	       }
	    }

	    if (residues.size() > 0) {
	       graphics_info_t g;
	       int imol_map = g.Imol_Refinement_Map();
	       if (! is_valid_map_molecule(imol_map)) {
		  add_status_bar_text("Refinement map not set");
	       } else {
		  // normal
		  mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
		  rr = g.refine_residues_vec(imol, residues, alt_conf.c_str(), mol);
	       }
               g.conditionally_wait_for_refinement_to_finish();
	    }
	 } else {
	    std::cout << "No residue specs found" << std::endl;
	 }
      }
   }
   return rr;
}





#ifdef USE_GUILE
SCM refine_residues_with_alt_conf_scm(int imol, SCM r, const char *alt_conf) { /* to be renamed later. */

   SCM rv = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      std::vector<coot::residue_spec_t> residue_specs = scm_to_residue_specs(r);
      g.residue_type_selection_was_user_picked_residue_range = false;
      coot::refinement_results_t rr =
	 refine_residues_with_alt_conf(imol, residue_specs, alt_conf);
      rv = g.refinement_results_to_scm(rr);
   }
   return rv;
}
#endif // USE_GUILE

// these functions need to call each other the other way round
//
#ifdef USE_GUILE
SCM refine_residues_with_modes_with_alt_conf_scm(int imol, SCM residues_spec_list_scm,
						 const char *alt_conf,
						 SCM mode_1,
						 SCM mode_2,
						 SCM mode_3) {
   return refine_residues_with_alt_conf_scm(imol, residues_spec_list_scm, alt_conf);
}
#endif // USE_GUILE

#ifdef USE_GUILE
SCM regularize_residues_with_alt_conf_scm(int imol, SCM res_spec_scm, const char *alt_conf) {

   SCM rv = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      std::vector<coot::residue_spec_t> residue_specs = scm_to_residue_specs(res_spec_scm);

      if (residue_specs.size() > 0) {
	 std::vector<mmdb::Residue *> residues;
	 for (unsigned int i=0; i<residue_specs.size(); i++) {
	    coot::residue_spec_t rs = residue_specs[i];
	    mmdb::Residue *r = graphics_info_t::molecules[imol].get_residue(rs);
	    if (r) {
	       residues.push_back(r);
	    }
	 }

	 if (residues.size() > 0) {
	    graphics_info_t g;
	    mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
	    g.residue_type_selection_was_user_picked_residue_range = false;
	    coot::refinement_results_t rr =
	       g.regularize_residues_vec(imol, residues, alt_conf, mol);
            g.conditionally_wait_for_refinement_to_finish();
	    rv = g.refinement_results_to_scm(rr);
	 }
      }
   }
   return rv;
}
#endif // USE_GUILE

#ifdef USE_GUILE
std::vector<coot::residue_spec_t> scm_to_residue_specs(SCM r) {
   std::vector<coot::residue_spec_t> residue_specs;
   SCM r_length_scm = scm_length(r);
   int r_length = scm_to_int(r_length_scm);
   for (int i=0; i<r_length; i++) {
      SCM res_spec_scm = scm_list_ref(r, scm_from_int(i));
      std::pair<bool, coot::residue_spec_t> res_spec =
	 make_residue_spec(res_spec_scm);
      if (res_spec.first) {
	 residue_specs.push_back(res_spec.second);
      }
   }
   return residue_specs;
}
#endif // USE_GUILE



#ifdef USE_PYTHON
PyObject *refine_residues_py(int imol, PyObject *r) {
   return refine_residues_with_alt_conf_py(imol, r, "");
}
#endif // USE_PYTHON



#ifdef USE_PYTHON
PyObject *regularize_residues_py(int imol, PyObject *r) { /* presumes the alt_conf is "". */
   return regularize_residues_with_alt_conf_py(imol, r, "");
}
#endif // USE_PYTHON


#ifdef USE_PYTHON
PyObject *refine_residues_with_alt_conf_py(int imol, PyObject *r,
					   const char *alt_conf) {
   PyObject *m1 = Py_False;
   PyObject *m2 = Py_False;
   PyObject *m3 = Py_False;
   Py_INCREF(m1);
   Py_INCREF(m2);
   Py_INCREF(m3);
   return refine_residues_with_modes_with_alt_conf_py(imol, r, alt_conf, m1, m2, m3);
}
#endif // USE_PYTHON

#ifdef USE_PYTHON
PyObject *refine_residues_with_modes_with_alt_conf_py(int imol, PyObject *res_specs_py,
						      const char *alt_conf,
						      PyObject *mode_1,
						      PyObject *mode_2,
						      PyObject *mode_3) {

   PyObject *rv = Py_False;
   if (is_valid_model_molecule(imol)) {
      std::vector<coot::residue_spec_t> residue_specs = py_to_residue_specs(res_specs_py);

      if (residue_specs.size() > 0) {
	 std::vector<mmdb::Residue *> residues;
	 for (unsigned int i=0; i<residue_specs.size(); i++) {
	    coot::residue_spec_t rs = residue_specs[i];
	    mmdb::Residue *r = graphics_info_t::molecules[imol].get_residue(rs);
	    if (r) {
	       residues.push_back(r);
	    }
	 }

	 if (residues.size() > 0) {
	    graphics_info_t g;
	    int imol_map = g.Imol_Refinement_Map();
	    if (! is_valid_map_molecule(imol_map)) {
	       add_status_bar_text("Refinement map not set");
	    } else {
	       // normal
	       mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;

	       bool soft_mode_hard_mode = false;
	       if (PyUnicode_Check(mode_1)) {
                  std::string s = PyBytes_AS_STRING(PyUnicode_AsUTF8String(mode_1));
		  if (s == "soft-mode/hard-mode")
		     soft_mode_hard_mode = true;
	       }

	     if (soft_mode_hard_mode) {
// 		double w_orig = g.geometry_vs_map_weight;
// 		g.geometry_vs_map_weight /= 50;
// 		coot::refinement_results_t rr =
// 		   g.refine_residues_vec(imol, residues, alt_conf, mol);
// 		g.geometry_vs_map_weight = w_orig;
// 		g.last_restraints.set_map_weight(w_orig);
// 		// continue with normal/hard mode
// 		g.drag_refine_refine_intermediate_atoms();
// 		// rv = g.refinement_results_to_py(rr);
	     } else {
		// normal
		g.residue_type_selection_was_user_picked_residue_range = false;
		coot::refinement_results_t rr =
		   g.refine_residues_vec(imol, residues, alt_conf, mol);
                g.conditionally_wait_for_refinement_to_finish();
		rv = g.refinement_results_to_py(rr);
	     }
          }
        }
      } else {
        std::cout << "No residue specs found" << std::endl;
      }
   }

   if (PyBool_Check(rv)) {
     Py_INCREF(rv);
   }

   return rv;
}
#endif // USE_PYTHON

#ifdef USE_PYTHON
PyObject *regularize_residues_with_alt_conf_py(int imol, PyObject *res_specs_py, const char *alt_conf) {

   PyObject *rv = Py_False;
   if (is_valid_model_molecule(imol)) {
      std::vector<coot::residue_spec_t> residue_specs = py_to_residue_specs(res_specs_py);

      if (residue_specs.size() > 0) {
	 std::vector<mmdb::Residue *> residues;
	 for (unsigned int i=0; i<residue_specs.size(); i++) {
	    coot::residue_spec_t rs = residue_specs[i];
	    mmdb::Residue *r = graphics_info_t::molecules[imol].get_residue(rs);
	    if (r) {
	       residues.push_back(r);
	    }
	 }

	 if (residues.size() > 0) {
	    graphics_info_t g;
        mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
        coot::refinement_results_t rr =
		  g.regularize_residues_vec(imol, residues, alt_conf, mol);
	       rv = g.refinement_results_to_py(rr);
	    }
	 }
      } else {
	 std::cout << "No residue specs found" << std::endl;
      }

   if (PyBool_Check(rv)) {
     Py_INCREF(rv);
   }
   return rv;
}
#endif // USE_PYTHON

/* Used by on_accept_reject_refinement_reject_button_clicked() */
void stop_refinement_internal() {

   graphics_info_t g;
   g.stop_refinement_internal();

}

void set_refinement_use_soft_mode_nbc_restraints(short int flag) {

   graphics_info_t g;
   g.set_use_harmonic_approximations_for_nbcs(flag);
   
}




#ifdef USE_PYTHON
std::vector<coot::residue_spec_t> py_to_residue_specs(PyObject *r) {

   std::vector<coot::residue_spec_t> residue_specs;
   int r_length = PyObject_Length(r);
   for (int i=0; i<r_length; i++) {
      PyObject *res_spec_py = PyList_GetItem(r, i);
      std::pair<bool, coot::residue_spec_t> res_spec = make_residue_spec_py(res_spec_py);
      if (res_spec.first) {
	 residue_specs.push_back(res_spec.second);
      }
   }
   return residue_specs;
}
#endif // USE_PYTHON

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
   std::string cmd = "update-go-to-atom-window-on-changed-mol";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
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

void delete_atom(int imol, const char *chain_id, int resno, const char *ins_code,
		 const char *at_name, const char *altLoc) {


   if (is_valid_model_molecule(imol)) {
      // after the atom is deleted from the molecule then the calling
      // char * pointers are out of date!  We can't use them.  So, save
      // them as strings before writing the history.
      graphics_info_t g;

      if (! chain_id) {
	 std::cout << "ERROR:: in delete_atom() trapped null chain_id\n";
	 return;
      }
      if (! ins_code) {
	 std::cout << "ERROR:: in delete_atom() trapped null ins_code\n";
	 return;
      }
      if (! at_name) {
	 std::cout << "ERROR:: in delete_atom() trapped null at_name\n";
	 return;
      }
      if (! altLoc) {
	 std::cout << "ERROR:: in delete_atom() trapped null altLoc\n";
	 return;
      }

      //
      std::string chain_id_string = chain_id;
      std::string ins_code_string = ins_code;
      std::string atom_name_string = at_name;
      std::string altloc_string = altLoc;

      mmdb::Residue *residue_p =
	 graphics_info_t::molecules[imol].get_residue(chain_id, resno, ins_code);
      if (residue_p) {
	 if (residue_p->GetNumberOfAtoms() > 1) {
	    coot::residue_spec_t spec(residue_p);
	    g.delete_residue_from_geometry_graphs(imol, spec);
	 } else {
	    // only one atom left, so better delete whole residue.
	    delete_residue(imol, chain_id, resno, ins_code);
	    return;
	 }
      }

      short int istat = g.molecules[imol].delete_atom(chain_id, resno, ins_code, at_name, altLoc);
      if (istat) {
	 // now if the go to atom widget was being displayed, we need to
	 // redraw the residue list and atom list (if the molecule of the
	 // residue and atom list is the molecule that has just been
	 // deleted)
	 //

	 g.update_go_to_atom_window_on_changed_mol(imol);
	 update_go_to_atom_residue_list(imol);
	 graphics_draw();
      } else {
	 std::cout << "failed to delete atom  chain_id: :" << chain_id
		   << ": " << resno << " incode :" << ins_code
		   << ": atom-name :" <<  at_name << ": altloc :" <<  altLoc << ":" << "\n";
      }

      std::string cmd = "delete-atom";
      std::vector<coot::command_arg_t> args;
      args.push_back(imol);
      args.push_back(coot::util::single_quote(chain_id_string));
      args.push_back(resno);
      args.push_back(coot::util::single_quote(ins_code_string));
      args.push_back(coot::util::single_quote(atom_name_string));
      args.push_back(coot::util::single_quote(altloc_string));
      add_to_history_typed(cmd, args);

   } else {
      std::cout << "ERROR:: Model number " << imol << " is not a valid molecule" << std::endl;
   }

}

void set_delete_atom_mode() {

   graphics_info_t g;
   g.delete_item_atom = 1;
   g.delete_item_residue_zone = 0;
   g.delete_item_residue_hydrogens = 0;
   g.delete_item_residue = 0;
   g.delete_item_sidechain = 0;
   g.delete_item_sidechain_range = 0;
   g.delete_item_chain = 0;
   pick_cursor_maybe();
   add_to_history_simple("set-delete-atom-mode");
}

void set_delete_residue_mode() {

   graphics_info_t g;
   g.delete_item_atom = 0;
   g.delete_item_residue_zone = 0;
   g.delete_item_residue_hydrogens = 0;
   g.delete_item_water = 0;
   g.delete_item_residue = 1;
   g.delete_item_sidechain = 0;
   g.delete_item_sidechain_range = 0;
   g.delete_item_chain = 0;
   pick_cursor_maybe();
   add_to_history_simple("set-delete-residue-mode");
}

void set_delete_residue_hydrogens_mode() {

   graphics_info_t g;
   g.delete_item_residue = 0;
   g.delete_item_residue_zone = 0;
   g.delete_item_atom = 0;
   g.delete_item_water = 0;
   g.delete_item_residue_hydrogens = 1;
   g.delete_item_sidechain = 0;
   g.delete_item_sidechain_range = 0;
   g.delete_item_chain = 0;
   pick_cursor_maybe();
   add_to_history_simple("set-delete-residue-hydrogens-mode");

}

void set_delete_residue_zone_mode() {

   graphics_info_t g;
   g.delete_item_residue = 0;
   g.delete_item_residue_zone = 1;
   g.delete_item_atom = 0;
   g.delete_item_water = 0;
   g.delete_item_residue_hydrogens = 0;
   g.delete_item_sidechain = 0;
   g.delete_item_sidechain_range = 0;
   g.delete_item_chain = 0;
   add_to_history_simple("set-delete-residue-zone-mode");
}

void set_delete_water_mode() {

   graphics_info_t g;
   g.delete_item_residue = 0;
   g.delete_item_residue_zone = 0;
   g.delete_item_water = 1;
   g.delete_item_atom = 0;
   g.delete_item_residue_hydrogens = 0;
   g.delete_item_sidechain = 0;
   g.delete_item_sidechain_range = 0;
   g.delete_item_chain = 0;
   pick_cursor_maybe();
   add_to_history_simple("set-delete-residue-water-mode");

}

void set_delete_sidechain_mode() {

   graphics_info_t g;
   std::cout << "set_delete_sidechain_mode " << std::endl;
   g.delete_item_residue = 0;
   g.delete_item_residue_zone = 0;
   g.delete_item_water = 0;
   g.delete_item_atom = 0;
   g.delete_item_residue_hydrogens = 0;
   g.delete_item_chain = 0;
   g.delete_item_sidechain = 1;
   g.delete_item_sidechain_range = 0;
   pick_cursor_maybe();
   add_to_history_simple("set-delete-sidechain-mode");

}

void set_delete_sidechain_range_mode() {

   graphics_info_t g;
   std::cout << "set_delete_sidechain_range_mode " << std::endl;
   g.delete_item_residue = 0;
   g.delete_item_residue_zone = 0;
   g.delete_item_water = 0;
   g.delete_item_atom = 0;
   g.delete_item_residue_hydrogens = 0;
   g.delete_item_chain = 0;
   g.delete_item_sidechain = 0;
   g.delete_item_sidechain_range = 1;
   pick_cursor_maybe();
   add_to_history_simple("set-delete-sidechain-range-mode");
}


void set_delete_chain_mode() {

   std::cout << "in set_delete_chain_mode()! " << std::endl;
   graphics_info_t g;
   g.delete_item_residue = 0;
   g.delete_item_residue_zone = 0;
   g.delete_item_water = 0;
   g.delete_item_atom = 0;
   g.delete_item_residue_hydrogens = 0;
   g.delete_item_sidechain = 0;
   g.delete_item_sidechain_range = 0;
   g.delete_item_chain = 1;
   pick_cursor_maybe();
   add_to_history_simple("set-delete-sidechain-mode");

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

short int delete_item_mode_is_sidechain_range_p() {
   short int v=0;
   if (graphics_info_t::delete_item_sidechain_range == 1)
      v = 1;
   return v;
}

short int delete_item_mode_is_chain_p() {
   short int v=0;
   if (graphics_info_t::delete_item_chain == 1)
      v = 1;
   return v;
}

// 20090911 Dangerous?  This used to be called by
// check_if_in_delete_item_define water mode.  But it caused a crash.
// It seemed to be accessing bad memory for the atom name.  That was
// because we were doing a delete_residue_hydrogens() between getting
// the atom index from the pick and calling this function.  Therefore
// the atom index/selection was out of date when we got here.  That
// has been fixed now, but it is a concern for other function that
// might want to do something similar.
//
// Delete this function when doing a cleanup.
//
void delete_atom_by_atom_index(int imol, int index, short int do_delete_dialog) {
   graphics_info_t g;

   if (index < g.molecules[imol].atom_sel.n_selected_atoms) {
      const char *atom_name = g.molecules[imol].atom_sel.atom_selection[index]->name;
      const char *chain_id  = g.molecules[imol].atom_sel.atom_selection[index]->GetChainID();
      const char *altconf   = g.molecules[imol].atom_sel.atom_selection[index]->altLoc;
      const char *ins_code  = g.molecules[imol].atom_sel.atom_selection[index]->GetInsCode();
      int resno             = g.molecules[imol].atom_sel.atom_selection[index]->GetSeqNum();

      mmdb::Residue *residue_p =
	 graphics_info_t::molecules[imol].get_residue(chain_id, resno, ins_code);
      if (residue_p) {
	 coot::residue_spec_t spec(residue_p);
	 g.delete_residue_from_geometry_graphs(imol, spec);
      }

      std::cout << "calling delete_atom() with args chain_id :"
		<< chain_id << ": resno " << resno << " inscode :"
		<< ins_code << ": atom-name " << atom_name << ": altconf :"
		<< altconf << ":" << std::endl;
      delete_atom(imol, chain_id, resno, ins_code, atom_name, altconf);
      // delete_object_handle_delete_dialog(do_delete_dialog);
   }

   // no need for this, the called delete_atom() does it.
//    std::string cmd = "delete-atom-by-atom-index";
//    std::vector<coot::command_arg_t> args;
//    args.push_back(imol);
//    args.push_back(index);
//    args.push_back(do_delete_dialog);
//    add_to_history_typed(cmd, args);

}

void delete_residue_by_atom_index(int imol, int index, short int do_delete_dialog_by_ctrl) {

   graphics_info_t g;
   int imodel            = g.molecules[imol].atom_sel.atom_selection[index]->GetModelNum();
   std::string chain_id  = g.molecules[imol].atom_sel.atom_selection[index]->GetChainID();
   int resno             = g.molecules[imol].atom_sel.atom_selection[index]->GetSeqNum();
   std::string altloc    = g.molecules[imol].atom_sel.atom_selection[index]->altLoc;
   std::string inscode   = g.molecules[imol].atom_sel.atom_selection[index]->GetInsCode();

   // I don't think that there is any need to call get_residue() here,
   // we can simply construct spec from chain_id, resno and inscode.
   // There are other places where we do this too (to delete a residue
   // from the geometry graphs).
   mmdb::Residue *residue_p = g.molecules[imol].get_residue(chain_id, resno, inscode);
   if (residue_p) {
      coot::residue_spec_t spec(residue_p);
      g.delete_residue_from_geometry_graphs(imol, spec);
   }

   if ((altloc == "") && (g.molecules[imol].atom_sel.mol->GetNumberOfModels() == 1))
      delete_residue(imol, chain_id.c_str(), resno, inscode.c_str());
   else
      delete_residue_with_full_spec(imol, imodel, chain_id.c_str(), resno,
				    inscode.c_str(), altloc.c_str());

   short int do_delete_dialog = do_delete_dialog_by_ctrl;
   // delete_object_handle_delete_dialog(do_delete_dialog);

   graphics_draw();
   std::string cmd = "delete-residue-by-atom-index";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(index);
   args.push_back(do_delete_dialog_by_ctrl);
   add_to_history_typed(cmd, args);
}

void delete_residue_hydrogens_by_atom_index(int imol, int index, short int do_delete_dialog) {

   graphics_info_t g;
   std::string chain_id  = g.molecules[imol].atom_sel.atom_selection[index]->GetChainID();
   int resno             = g.molecules[imol].atom_sel.atom_selection[index]->GetSeqNum();
   std::string altloc    = g.molecules[imol].atom_sel.atom_selection[index]->altLoc;
   std::string inscode   = g.molecules[imol].atom_sel.atom_selection[index]->GetInsCode();


   delete_residue_hydrogens(imol, chain_id.c_str(), resno, inscode.c_str(), altloc.c_str());

   // delete_object_handle_delete_dialog(do_delete_dialog);
   graphics_draw();
   std::string cmd = "delete-residue-hydrogens-by-atom-index";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(index);
   args.push_back(do_delete_dialog);
   add_to_history_typed(cmd, args);
}

// Deletes all altconfs, the whole residue goes.
//
void delete_residue_range(int imol, const char *chain_id, int resno_start, int resno_end) {

   // Note to self: do you want this or the graphics_info_t version?

   // altconf is ignored currently.
   //

   if (resno_start > resno_end) {
      std::swap(resno_start, resno_end);
   }
   coot::residue_spec_t res1(chain_id, resno_start);
   coot::residue_spec_t res2(chain_id, resno_end);

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      g.delete_residue_range(imol, res1, res2);

      // cheap! I should find the residues with insertion codes in this range too.
      // How to do that? Hmm... Needs a class function. This will do for now
      //
      std::vector<coot::residue_spec_t> res_specs;
      for (int i=resno_start; i<=resno_end; i++) {
	 coot::residue_spec_t r(chain_id, i, "");
	 res_specs.push_back(r);
      }

      g.delete_residues_from_geometry_graphs(imol, res_specs);
      if (graphics_info_t::go_to_atom_window) {
	 update_go_to_atom_window_on_changed_mol(imol);
      }
      if (! is_valid_model_molecule(imol))
	 g.delete_molecule_from_display_manager(imol, false);
   }
   graphics_draw();
   std::string cmd = "delete-residue-range";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(resno_start);
   args.push_back(resno_end);
   add_to_history_typed(cmd, args);
}




void clear_pending_delete_item() {

   graphics_info_t g;
   g.delete_item_atom = 0;
   g.delete_item_residue = 0;
   g.delete_item_residue_zone = 0;
   g.delete_item_residue_hydrogens = 0;
   add_to_history_simple("clear-pending-delete-item");
}




int move_molecule_to_screen_centre_internal(int imol) {

   int imoved_stat = 0;
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
      imoved_stat = 1;
      // finally show and activate
      set_mol_displayed(imol, 1);
      set_mol_active(imol, 1);

      if (0) {
	 std::cout << "-------------------- move_molecule_to_screen_centre_internal() "
		   << imol << std::endl;
	 std::cout << "           calling g.setup_graphics_ligand_view_aa() "
		   << std::endl;
      }
      g.setup_graphics_ligand_view_using_active_atom(imol); // only in imol
   }
   return imoved_stat;
}


void set_write_peaksearched_waters() {
   graphics_info_t g;
   g.ligand_water_write_peaksearched_atoms = 1;
   add_to_history_simple("set-write-peaksearched-waters");
}


void
place_atom_at_pointer() {
   graphics_info_t g;
   if (g.pointer_atom_is_dummy)
      g.place_dummy_atom_at_pointer();
   else {
      place_atom_at_pointer_by_window();
   }
   add_to_history_simple("place-atom-at-pointer");
}


int pointer_atom_molecule() {
   return graphics_info_t::user_pointer_atom_molecule;
}

int create_pointer_atom_molecule_maybe() {
   graphics_info_t g;
   return g.create_pointer_atom_molecule_maybe();
}

void
set_pointer_atom_molecule(int imol) {
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::user_pointer_atom_molecule = imol;
   }
}

void
place_typed_atom_at_pointer(const char *type) {
   graphics_info_t g;
   g.place_typed_atom_at_pointer(std::string(type));
   std::string cmd = "place-typed-atom-at-pointer";
   std::vector<coot::command_arg_t> args;
   args.push_back(single_quote(type));
   add_to_history_typed(cmd, args);
}

void add_an_atom(const std::string &element) {
   // same as above? but modern?
   graphics_info_t g;
   g.place_typed_atom_at_pointer(element);
   std::string cmd = "add-an-atom";
   std::vector<coot::command_arg_t> args;
   args.push_back(single_quote(element));
   add_to_history_typed(cmd, args);
}

#ifdef USE_PYTHON
void nudge_the_temperature_factors_py(int imol, PyObject *residue_spec_py, float amount) {
   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t residue_spec = residue_spec_from_py(residue_spec_py);
      graphics_info_t::molecules[imol].change_b_factors_of_residue_by(residue_spec, amount);
   }
}
#endif


void set_pointer_atom_is_dummy(int i) {
   graphics_info_t::pointer_atom_is_dummy = i;
   std::string cmd = "set-pointer-atom-is-dummy";
   std::vector<coot::command_arg_t> args;
   args.push_back(i);
   add_to_history_typed(cmd, args);
}



void display_where_is_pointer() {
   graphics_info_t g;
   g.display_where_is_pointer();
   add_to_history_simple("display-where-is-pointer");
}


// -----------------------------------------------------------------------------
//                               Automutation stuff
// -----------------------------------------------------------------------------
//
short int progressive_residues_in_chain_check(const char *chain_id, int imol) {

   std::string cmd = "progressive-residues-in-chain-check";
   std::vector<coot::command_arg_t> args;
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(imol);
   add_to_history_typed(cmd, args);

   graphics_info_t g;
   if (imol < graphics_n_molecules()) {
      return g.molecules[imol].progressive_residues_in_chain_check_by_chain(chain_id);
   } else {
      std::cout << "no such molecule number in progressive_residues_in_chain_check\n";
      return 0;
   }
}



// return -1 on error:
//
int chain_n_residues(const char *chain_id, int imol) {

   graphics_info_t g;
   std::string cmd = "chain-n-residues";
   std::vector<coot::command_arg_t> args;
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(imol);
   add_to_history_typed(cmd, args);
   
   if (is_valid_model_molecule(imol)) {
      return g.molecules[imol].chain_n_residues(chain_id);
   } else {
      return -1;
   }

}

// Return "" on failure.
//
std::string resname_from_serial_number(int imol, const char *chain_id, int serial_num) {

   std::string r;
   if (is_valid_model_molecule(imol)) {
      if (serial_num >= 0) {
	 unsigned int uisn = serial_num;
	 r = graphics_info_t::molecules[imol].res_name_from_serial_number(chain_id, uisn);
      }
   }
   std::string cmd = "resname-from-serial-number";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(serial_num);
   add_to_history_typed(cmd, args);
   return r;
}

// Return < -9999 on failure
int  seqnum_from_serial_number(int imol, const char *chain_id, int serial_num) {

   int UNSET_SERIAL_NUMBER = -10000;
   int iseqnum = UNSET_SERIAL_NUMBER;

   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      int nchains = mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) {
	 mmdb::Chain *chain_p = mol->GetChain(1,ichain);
	 std::string mol_chain_id(chain_p->GetChainID());
	 if (mol_chain_id == std::string(chain_id)) {
	    int nres = chain_p->GetNumberOfResidues();
	    if ((serial_num < nres) && (serial_num >= 0)) {
	       int ch_n_res;
	       mmdb::PResidue *residues;
	       chain_p->GetResidueTable(residues, ch_n_res);
	       mmdb::Residue *this_res = residues[serial_num];
	       iseqnum = this_res->GetSeqNum();
	    } else {
	       std::cout << "WARNING:: seqnum_from_serial_number: requested residue with serial_num "
			 << serial_num << " but only " << nres << " residues in chain "
			 << mol_chain_id << std::endl;
	    }
	 }
      }
      if (iseqnum == UNSET_SERIAL_NUMBER) {
	 std::cout << "WARNING: seqnum_from_serial_number: returning UNSET serial number "
		   << std::endl;
      }
   } else {
      std::cout << "WARNING molecule number " << imol << " is not a valid model molecule "
		<< std::endl;
   }
   std::string cmd = "seqnum-from-serial-number";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(serial_num);
   add_to_history_typed(cmd, args);
   return iseqnum;
}

//! \brief return the serial number of the specified residue
//!
//! @return -1 on failure to find the residue
//
int serial_number_from_residue_specs(int imol, const std::string &chain_id, int res_no, const std::string &ins_code) {

   int serial_number = -1;

   if (is_valid_model_molecule(imol)) {
      serial_number = graphics_info_t::molecules[imol].residue_serial_number(chain_id, res_no, ins_code);
   }

   return serial_number;

}


char *insertion_code_from_serial_number(int imol, const char *chain_id, int serial_num) {

   char *r = NULL;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      int nchains = mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) {
	 mmdb::Chain *chain_p = mol->GetChain(1,ichain);
	 std::string mol_chain_id(chain_p->GetChainID());
	 if (mol_chain_id == std::string(chain_id)) {
	    int nres = chain_p->GetNumberOfResidues();
	    if (serial_num < nres) {
	       int ch_n_res;
	       mmdb::PResidue *residues;
	       chain_p->GetResidueTable(residues, ch_n_res);
	       mmdb::Residue *this_res = residues[serial_num];
	       r = this_res->GetInsCode();
	    }
	 }
      }
   }
   std::string cmd = "insertion-code-from-serial-number";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(serial_num);
   add_to_history_typed(cmd, args);
   return r;
}

#ifdef USE_GUILE
SCM
chain_id_scm(int imol, int ichain) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      mmdb::Chain *chain_p = mol->GetChain(1,ichain);
      if (chain_p)
	 r = scm_from_locale_string(chain_p->GetChainID());
   }
   std::string cmd = "chain_id";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(ichain);
   add_to_history_typed(cmd, args);
   return r;
}
#endif // USE_GUILE


#ifdef USE_PYTHON
PyObject *
chain_id_py(int imol, int ichain) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      mmdb::Chain *chain_p = mol->GetChain(1,ichain);
      if (chain_p)
	 r = myPyString_FromString(chain_p->GetChainID());
   }
   std::string cmd = "chain_id";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(ichain);
   add_to_history_typed(cmd, args);

   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif // USE_GUILE

int n_chains(int imol) {

   int nchains = -1;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      nchains = mol->GetNumberOfChains(1);
   }
   std::string cmd = "n-chains";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
   return nchains;
}

#ifdef USE_PYTHON
/*! \brief get the chain ids of molecule number imol

  @param imol is the molecule index
  @return a list of the the chain ids or None on failure
*/
PyObject *get_chain_ids_py(int imol) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      std::vector<std::string> chain_ids;
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            chain_ids.push_back(chain_p->GetChainID());
         }
      }
      if (! chain_ids.empty()) {
         r = PyList_New(chain_ids.size());
         for (unsigned int i=0; i<chain_ids.size(); i++) {
            PyObject *o = myPyString_FromString(chain_ids[i].c_str());
            PyList_SetItem(r, i, o);
         }
      }
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif



/*! \brief return the number of models in molecule number imol

useful for NMR or other such multi-model molecules.

return the number of models or -1 if there was a problem with the
given molecule.
*/
int n_models(int imol) {

   int r = -1; // fail;
   if (is_valid_model_molecule(imol)) {
      r = graphics_info_t::molecules[imol].n_models();
   }
   std::string cmd = "n-models";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
   return r;
}

/*!\brief return the number of residues in the molecule,

return -1 if this is a map or closed.
 */
int n_residues(int imol) {

   int r = -1;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      r = g.molecules[imol].n_residues();
   }
   return r;
}

int n_atoms(int imol) {

   int r = -1;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      r = g.molecules[imol].n_atoms();
   }
   return r;
}



// return -1 on error (e.g. chain_id not found, or molecule is not a
// model), 0 for no, 1 for is.
int is_solvent_chain_p(int imol, const char *chain_id) {

   int r = -1;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      int nchains = mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) {
	 mmdb::Chain *chain_p = mol->GetChain(1,ichain);
	 std::string mol_chain_id(chain_p->GetChainID());
	 if (mol_chain_id == std::string(chain_id)) {
	    r = chain_p->isSolventChain();
	 }
      }
   }
   std::string cmd = "is-solvent-chain-p";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   add_to_history_typed(cmd, args);
   return r;
}

// return -1 on error (e.g. chain_id not found, or molecule is not a
// model), 0 for no, 1 for is.
int is_protein_chain_p(int imol, const char *chain_id) {

   int r = -1;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      int nchains = mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) {
	 mmdb::Chain *chain_p = mol->GetChain(1,ichain);
	 std::string mol_chain_id(chain_p->GetChainID());
	 if (mol_chain_id == std::string(chain_id)) {
	    r = chain_p->isAminoacidChain();
	 }
      }
   }
   std::string cmd = "is-protein-chain-p";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   add_to_history_typed(cmd, args);

   return r;
}

// return -1 on error (e.g. chain_id not found, or molecule is not a
// model), 0 for no, 1 for is.
int is_nucleotide_chain_p(int imol, const char *chain_id) {

   int r = 0;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      int nchains = mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) {
	 mmdb::Chain *chain_p = mol->GetChain(1,ichain);
	 std::string mol_chain_id(chain_p->GetChainID());

	 if (mol_chain_id == std::string(chain_id)) {
	    r = chain_p->isNucleotideChain();
	    break;
	 }
      }
   }


   std::string cmd = "is-nucleotide-chain-p";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   add_to_history_typed(cmd, args);
   return r;
}

//
// /*! \brief sort the chain ids of the imol-th molecule in lexographical order */
void sort_chains(int imol) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].sort_chains();
      if (graphics_info_t::use_graphics_interface_flag) {
	 graphics_info_t g;
	 if (g.go_to_atom_window) {
	    g.update_go_to_atom_window_on_changed_mol(imol);
	 }
      }
   }
}

// /*! \brief sort the residues of the imol-th molecule */
void sort_residues(int imol) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].sort_residues();
      if (graphics_info_t::use_graphics_interface_flag) {
        graphics_info_t g;
        if (g.go_to_atom_window) {
           g.update_go_to_atom_window_on_changed_mol(imol);
        }
      }
   }
}


/*! \brief simply print secondardy structure info to the
  terminal/console.  In future, this could/should return the info.  */
void print_header_secondary_structure_info(int imol) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].print_secondary_structure_info();
   }
}

#include "coot-utils/secondary-structure-headers.hh"

//
void write_header_secondary_structure_info(int imol, const char *file_name) {

   if (is_valid_model_molecule(imol)) {
      mmdb::io::File f;
      bool Text = true;
      f.assign(file_name, Text);
      if (f.rewrite()) {
	 mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
	 int n_models = mol->GetNumberOfModels();
	 int imod = 1;
	 mmdb::Model *model_p = mol->GetModel(imod);
	 int ss_status = model_p->CalcSecStructure(1);

	 // now do something clever with the residues of model_p to
	 // build mmdb sheet and strand objects and attach them to the
	 // model (should be part of mmdb).

	 coot::secondary_structure_header_records ssr(mol, false);

	 if (ss_status == mmdb::SSERC_Ok) {
	    std::cout << "INFO:: SSE status was OK\n";
	    model_p->PDBASCIIDumpPS(f); // dump CHelix and CStrand records.
	 }
      }
   }
}

void add_header_secondary_structure_info(int imol) {

   if (is_valid_model_molecule(imol))
      graphics_info_t::molecules[imol].add_secondary_structure_header_records();
}





/*  ----------------------------------------------------------------------- */
/*                         Renumber residue range                           */
/*  ----------------------------------------------------------------------- */

int renumber_residue_range(int imol, const char *chain_id,
			   int start_res, int last_res, int offset) {

   int i=0;
   if (imol >= 0) {
      if (imol <= graphics_info_t::n_molecules()) {
	 if (graphics_info_t::molecules[imol].has_model()) {
	    i = graphics_info_t::molecules[imol].renumber_residue_range(chain_id,
									start_res,
									last_res,
									offset);

	    if (i) {
	       graphics_info_t g;
	       graphics_draw();
	       g.update_go_to_atom_window_on_changed_mol(imol);
	       g.update_validation(imol);
	    }
	 }
      }
   }
   std::string cmd = "renumber-residue-range";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(start_res);
   args.push_back(last_res);
   args.push_back(offset);
   add_to_history_typed(cmd, args);
   return i;
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
      graphics_info_t g;
      g.update_go_to_atom_window_on_changed_mol(imol);
      g.update_validation(imol);
   }
   std::string cmd = "change-residue-number";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(current_resno);
   args.push_back(coot::util::single_quote(current_inscode));
   args.push_back(new_resno);
   args.push_back(coot::util::single_quote(new_inscode));
   add_to_history_typed(cmd, args);
   return idone;
}


/*  ----------------------------------------------------------------------- */
/*                         backup                                           */
/*  ----------------------------------------------------------------------- */

void make_backup(int imol) {

   if (is_valid_model_molecule(imol)) {
      if (graphics_info_t::molecules[imol].has_model()) {
	 graphics_info_t::molecules[imol].make_backup_from_outside();
      } else {
	 std::cout << "No model for this molecule" << std::endl;
      }
   } else {
      std::cout << "No model :" << imol << std::endl;
   }
   std::string cmd = "make-backup";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
}

int backup_state(int imol) {

   int istate = -1;

   if (is_valid_model_molecule(imol)) {
      if (graphics_info_t::molecules[imol].has_model()) {
	 istate = graphics_info_t::molecules[imol].backups_state();
      } else {
	 std::cout << "No model for this molecule" << std::endl;
      }
      } else {
      std::cout << "No model :" << imol << std::endl;
   }
   std::string cmd = "backup-state";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
   return istate;
}

void set_have_unsaved_changes(int imol) {

   if (is_valid_model_molecule(imol)) {
      if (graphics_info_t::molecules[imol].has_model()) {
	 graphics_info_t::molecules[imol].set_have_unsaved_changes_from_outside();
      }
   }
   std::string cmd = "set-have-unsaved-changes";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
}

int have_unsaved_changes_p(int imol) {

   int r = -1;
   if (is_valid_model_molecule(imol)) {
      if (graphics_info_t::molecules[imol].has_model()) {
	 r = graphics_info_t::molecules[imol].Have_unsaved_changes_p();
      }
   }
   return r;

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
   std::string cmd = "write-pdb-file";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(file_name));
   add_to_history_typed(cmd, args);
   return istat;
}

// return 1 on error, 0 on success
int
write_cif_file(int imol, const char *file_name) {

   graphics_info_t g;
   int istat = 0;
   if (is_valid_model_molecule(imol)) {
      istat = g.molecules[imol].write_cif_file(std::string(file_name));
   }
   std::string cmd = "write-cif-file";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(file_name));
   add_to_history_typed(cmd, args);
   return istat;
}

/*! \brief write molecule number imol's residue range as a PDB to file
  file_name */
/*  return 0 on success, -1 on error. */
int
write_residue_range_to_pdb_file(int imol, const char *chain_id,
				int resno_start, int resno_end,
				const char *filename) {

   int istat = 1;
   if (is_valid_model_molecule(imol)) {
      std::string chain(chain_id);
      if (resno_end < resno_start) {
	 int tmp = resno_end;
	 resno_end = resno_start;
	 resno_start = tmp;
      }
      mmdb::Manager *mol =
	 graphics_info_t::molecules[imol].get_residue_range_as_mol(chain, resno_start, resno_end);
      if (mol) {
        mmdb_manager_delete_conect(mol);
	 istat = mol->WritePDBASCII(filename);
	 delete mol; // give back the memory.
      }
   }
   std::string cmd = "write-residue-range-to-pdb-file";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(resno_start);
   args.push_back(resno_end);
   args.push_back(coot::util::single_quote(filename));
   add_to_history_typed(cmd, args);
   return istat;
}

int write_chain_to_pdb_file(int imol, const char *chain_id, const char *filename) {

   int istat = 0;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
      int SelHnd = mol->NewSelection(); // d
      mol->SelectAtoms(SelHnd, 1,
		       chain_id,
		       mmdb::ANY_RES, "*",
		       mmdb::ANY_RES, "*",
		       "*", // any residue name
		       "*", // atom name
		       "*", // elements
		       "*",  // alt loc.
		       mmdb::SKEY_NEW);
      mmdb::Manager *new_mol = coot::util::create_mmdbmanager_from_atom_selection(mol, SelHnd);
      if (new_mol) {
	 istat = new_mol->WritePDBASCII(filename);
	 delete new_mol;
      }
      mol->DeleteSelection(SelHnd);
   }
   std::string cmd = "write-chain-to-pdb-file";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(coot::util::single_quote(filename));
   add_to_history_typed(cmd, args);
   return istat;
}


/*! \brief save all modified coordinates molecules to the default
  names and save the state too. */
int quick_save() {

   // std::cout << "Quick save..." << std::endl;
   graphics_info_t g;
   g.quick_save();
   return 0;
}


/*! \brief return the state of the write_conect_records_flag.
  */
int get_write_conect_record_state() {
  graphics_info_t g;
  return g.write_conect_records_flag;
}
/*! \} */

/*! \brief set the flag to write (or not) conect records to the PDB file.
  */
void set_write_conect_record_state(int state) {
  graphics_info_t g;
  if (state == 1) {
    g.write_conect_records_flag = 1;
  } else {
    g.write_conect_records_flag = 0;
  }
}

/*! \} */




short int
add_OXT_to_residue(int imol, const char *chain_id, int resno, const char *insertion_code) {

   short int istat = -1;
   if (is_valid_model_molecule(imol)) {
      if (insertion_code) {
	 if (chain_id) {
	    graphics_info_t g;
	    istat = graphics_info_t::molecules[imol].add_OXT_to_residue(resno, std::string(insertion_code),
									std::string(chain_id),
									g.Geom_p());
	    g.molecules[imol].update_symmetry();
	    graphics_draw();
	 }
      }
   } else {
      std::cout << "WARNING:: in add_OXT_to_residue() imol " << imol << " is not valid" << std::endl;
   }
   std::string cmd = "add-OXT-to-residue";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(resno);
   args.push_back(coot::util::single_quote(insertion_code));
   args.push_back(coot::util::single_quote(chain_id));
   add_to_history_typed(cmd, args);
   std::cout << "debug:: add_OXT_to_residue() returns istat " << istat << std::endl;
   return istat;
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
   add_to_history_simple("alt-conf-split-type-number");
   return graphics_info_t::alt_conf_split_type_number();
}


void set_add_alt_conf_split_type_number(short int i) {
   graphics_info_t::alt_conf_split_type = i;
   std::string cmd = "set-add-alt-conf-split-type-number";
   std::vector<coot::command_arg_t> args;
   args.push_back(i);
   add_to_history_typed(cmd, args);

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


void set_show_alt_conf_intermediate_atoms(int i) {
   graphics_info_t::show_alt_conf_intermediate_atoms_flag = i;
   std::string cmd = "set-show-alt-conf-intermediate-atoms";
   std::vector<coot::command_arg_t> args;
   args.push_back(i);
   add_to_history_typed(cmd, args);
}

int show_alt_conf_intermediate_atoms_state() {
   add_to_history_simple("show-alt-conf-intermediate-atoms-state");
   return graphics_info_t::show_alt_conf_intermediate_atoms_flag;
}

void zero_occupancy_residue_range(int imol, const char *chain_id, int ires1, int ires2) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_occupancy_residue_range(std::string(chain_id), ires1, ires2, 0.0);
   } else {
      std::cout << "WARNING:: invalid model molecule number in zero_occupancy_residue_range "
		<< imol << std::endl;
   }
   std::string cmd = "zero-occupancy-residue-range";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(ires1);
   args.push_back(ires2);
   add_to_history_typed(cmd, args);
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
   std::string cmd = "fill-occupancy-residue-range";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(ires1);
   args.push_back(ires2);
   add_to_history_typed(cmd, args);
}


void set_occupancy_residue_range(int imol, const char *chain_id, int ires1, int ires2,
				 float occ) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_occupancy_residue_range(std::string(chain_id), ires1, ires2, occ);
   } else {
      std::cout << "WARNING:: invalid model molecule number in set_occupancy_residue_range "
		<< imol << std::endl;
   }
   graphics_draw();
   std::string cmd = "set-occupancy-residue-range";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(ires1);
   args.push_back(ires2);
   args.push_back(occ);
   add_to_history_typed(cmd, args);
}


void
set_b_factor_residue_range(int imol, const char *chain_id, int ires1, int ires2, float bval) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_b_factor_residue_range(std::string(chain_id), ires1, ires2, bval);
   } else {
      std::cout << "WARNING:: invalid model molecule number in set_b_factor_residue_range "
		<< imol << std::endl;
   }
   graphics_draw();
   std::string cmd = "set-b-factor-residue-range";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(ires1);
   args.push_back(ires2);
   args.push_back(bval);
   add_to_history_typed(cmd, args);

}

void
reset_b_factor_residue_range(int imol, const char *chain_id, int ires1, int ires2) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_b_factor_residue_range(std::string(chain_id), ires1, ires2,
								  graphics_info_t::default_new_atoms_b_factor);
   } else {
      std::cout << "WARNING:: invalid model molecule number in reset_b_factor_residue_range "
		<< imol << std::endl;
   }
   graphics_draw();
   std::string cmd = "reset-b-factor-residue-range";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(ires1);
   args.push_back(ires2);
   add_to_history_typed(cmd, args);

}

#ifdef USE_PYTHON
void set_b_factor_residues_py(int imol, PyObject *residue_specs_b_value_tuple_list_py) {

   if (is_valid_model_molecule(imol)) {
      if (PyList_Check(residue_specs_b_value_tuple_list_py)) {
	 unsigned int l = PyObject_Length(residue_specs_b_value_tuple_list_py);
	 if (l > 0) {
	    std::vector<std::pair<coot::residue_spec_t, double> > rbs;
	    for (unsigned int i=0; i<l; i++) {
	       PyObject *tuple_py = PyList_GetItem(residue_specs_b_value_tuple_list_py, i);
	       if (PyTuple_Check(tuple_py)) {
		  unsigned int l2 = PyObject_Length(tuple_py);
		  if (l2 == 2) {
		     PyObject *spec_py = PyTuple_GetItem(tuple_py, 0);
		     PyObject *bfac_py = PyTuple_GetItem(tuple_py, 1);
		     if (PyFloat_Check(bfac_py) || PyLong_Check(bfac_py)) {
			coot::residue_spec_t spec = residue_spec_from_py(spec_py);
			double b = PyFloat_AsDouble(bfac_py);
			std::pair<coot::residue_spec_t, double> p(spec, b);
			rbs.push_back(p);
		     }
		  }
	       }
	    }
	    graphics_info_t::molecules[imol].set_b_factor_residues(rbs);
	 }
      }
   }
}
#endif // USE_PYTHON

#ifdef USE_GUILE
void set_b_factor_residues_scm(int imol, SCM residue_specs_b_value_tuple_list_scm) {

   if (is_valid_model_molecule(imol)) {
      if (scm_is_true(scm_list_p(residue_specs_b_value_tuple_list_scm))) {
	 SCM l_scm = scm_length(residue_specs_b_value_tuple_list_scm);
	 unsigned int l = scm_to_int(l_scm);
	 if (l > 0) {
	    std::vector<std::pair<coot::residue_spec_t, double> > rbs;
	    for (unsigned int i=0; i<l; i++) {
	       SCM item_scm = scm_list_ref(residue_specs_b_value_tuple_list_scm,
					   scm_from_int(l));
	       if (scm_is_true(scm_list_p(item_scm))) {
		  SCM l2_scm = scm_length(item_scm);
		  unsigned int l2 = scm_to_int(l2_scm);
		  if (l2 == 2) {
		     SCM spec_scm = scm_list_ref(item_scm, scm_from_int(0));
		     SCM    b_scm = scm_list_ref(item_scm, scm_from_int(1));
		     coot::residue_spec_t spec = residue_spec_from_scm(spec_scm);
		     double b = scm_to_double(b_scm);
		     std::pair<coot::residue_spec_t, double> p(spec, b);
		     rbs.push_back(p);
		  }
	       }
	    }
	    graphics_info_t::molecules[imol].set_b_factor_residues(rbs);
	 }
      }
   }
}
#endif // USE_GUILE




void translate_molecule_by(int imol, float x, float y, float z) {

   if (is_valid_model_molecule(imol)) {
      if (graphics_info_t::molecules[imol].has_model()) {
	 graphics_info_t::molecules[imol].translate_by(x, y, z);
      }
   }
   graphics_draw();
}

/*! \brief transform molecule number imol by the given rotation
  matrix, then translate by (x,y,z) in Angstroms  */
void transform_molecule_by(int imol,
			   float m11, float m12, float m13,
			   float m21, float m22, float m23,
			   float m31, float m32, float m33,
			   float x, float y, float z) {

   if (is_valid_model_molecule(imol)) {
      clipper::Mat33<double> clipper_mat(m11, m12, m13,
					 m21, m22, m23,
					 m31, m32, m33);
      clipper::Coord_orth cco(x,y,z);
      clipper::RTop_orth rtop(clipper_mat, cco);
      graphics_info_t::molecules[imol].transform_by(rtop);
   }
   graphics_draw();

}

void transform_zone(int imol, const char *chain_id, int resno_start, int resno_end, const char *ins_code,
		    float m11, float m12, float m13,
		    float m21, float m22, float m23,
		    float m31, float m32, float m33,
		    float x, float y, float z) {

   if (is_valid_model_molecule(imol)) {

      clipper::Mat33<double> clipper_mat(m11, m12, m13,
					 m21, m22, m23,
					 m31, m32, m33);
      clipper::Coord_orth cco(x,y,z);
      clipper::RTop_orth rtop(clipper_mat, cco);
      bool do_backup = 1;
      graphics_info_t::molecules[imol].transform_zone_by(chain_id, resno_start, resno_end, ins_code, rtop,
							 do_backup);

      std::string cmd = "transform-zone";
      std::vector<coot::command_arg_t> args;
      args.push_back(imol);
      args.push_back(chain_id);
      args.push_back(resno_start);
      args.push_back(resno_end);
      args.push_back(ins_code);
      args.push_back(m11);
      args.push_back(m12);
      args.push_back(m13);
      args.push_back(m21);
      args.push_back(m22);
      args.push_back(m23);
      args.push_back(m31);
      args.push_back(m32);
      args.push_back(m33);
      args.push_back(x);
      args.push_back(y);
      args.push_back(z);
      add_to_history_typed(cmd, args);
   }
}



// Sequenc utils

void assign_fasta_sequence(int imol, const char *chain_id_in, const char *seq) {

   // format "> name \n <sequence>"
   if (is_valid_model_molecule(imol)) {
      const std::string chain_id = chain_id_in;
      graphics_info_t::molecules[imol].assign_fasta_sequence(chain_id, std::string(seq));
   }
}

void assign_pir_sequence(int imol, const char *chain_id_in, const char *seq) {

   if (is_valid_model_molecule(imol)) {
      const std::string chain_id = chain_id_in;
      graphics_info_t::molecules[imol].assign_pir_sequence(chain_id, std::string(seq));
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("assign-pir-sequence");
   command_strings.push_back(coot::util::int_to_string(imol));
   command_strings.push_back(single_quote(chain_id_in));
   command_strings.push_back(single_quote(seq));
   add_to_history(command_strings);
}

/*! \brief Associate the sequence to the molecule - to be used later for sequence assignment (.c.f pir file)   */
void associate_sequence_from_file(int imol, const char *file_name) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].associate_sequence_from_file(std::string(file_name));
   }
   std::vector<std::string> command_strings;
   command_strings.push_back("associate-sequence-from-file");
   command_strings.push_back(coot::util::int_to_string(imol));
   command_strings.push_back(single_quote(file_name));
   add_to_history(command_strings);
}


void assign_sequence_from_file(int imol, const char *file) {
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].assign_sequence_from_file(std::string(file));
   } else {
      std::cout << "WARNING:: assign_sequence_from_file() molecule number " << imol
                << " is not a valid molecule" << std::endl;
   }
   std::string cmd = "assign-sequence-from-file";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(single_quote(file));
   add_to_history_typed(cmd, args);
}

void assign_sequence_from_string(int imol, const char *chain_id_in, const char *seq) {
   if (is_valid_model_molecule(imol)) {
    const std::string chain_id = chain_id_in;
    graphics_info_t::molecules[imol].assign_sequence_to_NCS_related_chains_from_string(chain_id, std::string(seq));
  }
  std::string cmd = "assign-sequence-from-string";
  std::vector<coot::command_arg_t> args;
  args.push_back(imol);
  args.push_back(single_quote(chain_id_in));
  args.push_back(single_quote(seq));
  add_to_history_typed(cmd, args);
}

void delete_all_sequences_from_molecule(int imol) {

   if (is_valid_model_molecule(imol)) {
      if ((graphics_info_t::molecules[imol].sequence_info()).size() > 0) {
	 std::cout <<"BL DEBUG:: we have sequence info"<<std::endl;
	 graphics_info_t::molecules[imol].delete_all_sequences_from_molecule();
      } else {
	 std::cout <<"BL DEBUG:: no sequence info"<<std::endl;
      }
   }
}

void delete_sequence_by_chain_id(int imol, const char *chain_id_in) {
   if (is_valid_model_molecule(imol)) {
      if ((graphics_info_t::molecules[imol].sequence_info()).size() > 0) {
	 std::cout <<"BL DEBUG:: we have sequence info"<<std::endl;
	 const std::string chain_id = chain_id_in;
	 graphics_info_t::molecules[imol].delete_sequence_by_chain_id(chain_id);
      } else {
	 std::cout <<"BL DEBUG:: no sequence info"<<std::endl;
      }
   }
}


/*  ----------------------------------------------------------------------- */
/*                  trim                                                    */
/*  ----------------------------------------------------------------------- */
void
trim_molecule_by_map(int imol_coords, int imol_map,
		     float map_level, int delete_or_zero_occ_flag) {

   graphics_info_t g;
   if (is_valid_model_molecule(imol_coords)) {
      if (is_valid_map_molecule(imol_map)) {
	 if (g.molecules[imol_map].has_xmap()) {
	    int iv = g.molecules[imol_coords].trim_by_map(g.molecules[imol_map].xmap,
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


/*! \brief trim the molecule by the value in the B-factor column */
void trim_molecule_by_b_factor(int imol, float limit, short int keep_higher) {

   if (is_valid_model_molecule(imol)) {
      bool keep_higher_flag = keep_higher;
      graphics_info_t::molecules[imol].trim_molecule_by_b_factor(limit, keep_higher_flag);
   } else {
      std::cout << "WARNING:: " << imol << " is not a valid model molecule" << std::endl;
   }
   graphics_draw();

}

/*! \brief convert the value in the B-factor column (typically pLDDT for AlphaFold models) to a temperature factor */
void pLDDT_to_b_factor(int imol) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].pLDDT_to_b_factor();
   } else {
      std::cout << "WARNING:: " << imol << " is not a valid model molecule" << std::endl;
   }
   graphics_draw();
}


/*  ----------------------------------------------------------------------- */
/*               Simplex Refinement                                         */
/*  ----------------------------------------------------------------------- */
//
// Perhaps this should be a just a call to a graphics_info_t function?
//
void
fit_residue_range_to_map_by_simplex(int res1, int res2, const char *altloc,
				    const char *chain_id, int imol, int imol_for_map) {


   // The molecule_class_info_t updates its bonds.
   //
   if (is_valid_model_molecule(imol)) {
      if (graphics_info_t::molecules[imol].has_model()) {
	 if (is_valid_map_molecule(imol_for_map)) {
	    if (graphics_info_t::molecules[imol_for_map].has_xmap()) {
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
score_residue_range_fit_to_map(int res1, int res2, const char *altloc,
			       const char *chain_id, int imol, int imol_for_map) {

   float f = 0.0;

   if (is_valid_model_molecule(imol)) {
      if (graphics_info_t::molecules[imol].has_model()) {
	 if (is_valid_map_molecule(imol_for_map)) {
	    if (graphics_info_t::molecules[imol_for_map].has_xmap()) {
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

#ifdef USE_GUILE
// e.g. atom_spec: '("A" 81 "" " CA " "")
//      position   '(2.3 3.4 5.6)
SCM drag_intermediate_atom_scm(SCM atom_spec, SCM position) {
   SCM retval = SCM_BOOL_F;
   std::pair<bool, coot::atom_spec_t> p = make_atom_spec(atom_spec);
   if (p.first) {
      SCM pos_length_scm = scm_length(position);
      int pos_length = scm_to_int(pos_length_scm);
      if (pos_length == 3) {
	 SCM x_scm = scm_list_ref(position, scm_from_int(0));
	 SCM y_scm = scm_list_ref(position, scm_from_int(1));
	 SCM z_scm = scm_list_ref(position, scm_from_int(2));
	 double x = scm_to_double(x_scm);
	 double y = scm_to_double(y_scm);
	 double z = scm_to_double(z_scm);
	 clipper::Coord_orth pt(x,y,z);
	 graphics_info_t::drag_intermediate_atom(p.second, pt);
      }
   } else {
      std::cout << "WARNING:: bad atom spec in drag_intermediate_atom_scm() " << std::endl;
   }
   return retval;
}
#endif

#ifdef USE_GUILE
SCM mark_atom_as_fixed_scm(int imol, SCM atom_spec, int state) {
   SCM retval = SCM_BOOL_F;
   std::pair<bool, coot::atom_spec_t> p = make_atom_spec(atom_spec);
   if (p.first) {
      graphics_info_t::mark_atom_as_fixed(imol, p.second, state);
      graphics_draw();
   }
   return retval;
}
#endif


#ifdef USE_PYTHON
PyObject *drag_intermediate_atom_py(PyObject *atom_spec, PyObject *position) {
// e.g. atom_spec: ["A", 81, "", " CA ", ""]
//      position   [2.3, 3.4, 5.6]
   PyObject *retval = Py_False;
   PyObject *x_py;
   PyObject *y_py;
   PyObject *z_py;
   std::pair<bool, coot::atom_spec_t> p = make_atom_spec_py(atom_spec);
   if (p.first) {
      int pos_length = PyObject_Length(position);
      if (pos_length == 3) {
	 x_py = PyList_GetItem(position, 0);
	 y_py = PyList_GetItem(position, 1);
	 z_py = PyList_GetItem(position, 2);
	 double x = PyFloat_AsDouble(x_py);
	 double y = PyFloat_AsDouble(y_py);
	 double z = PyFloat_AsDouble(z_py);
	 clipper::Coord_orth pt(x,y,z);
	 graphics_info_t::drag_intermediate_atom(p.second, pt);
	 retval = Py_True; // Shall we return True if atom is dragged?
      }
   }

   Py_INCREF(retval);
   return retval;
}
#endif // USE_PYTHON

#ifdef USE_PYTHON
//! \brief add a target position for an intermediate atom and refine
//
// A function requested by Hamish.
// This applies to intermediate atoms (add_extra_target_position_restraint)
// does not. This activates refinement after the restraint is added (add_extra_target_position_restraint
// does not).
//
// We need a vector (of atom specs) version of this so that we don't keep stopping and starting
// the refinement as we add new pull restraints (maybe there will be 50 of them or so)
//
PyObject *add_target_position_restraint_for_intermediate_atom_py(PyObject *atom_spec, PyObject *position) {

// e.g. atom_spec: ["A", 81, "", " CA ", ""]
//      position   [2.3, 3.4, 5.6]
   PyObject *retval = Py_False;
   std::pair<bool, coot::atom_spec_t> p = make_atom_spec_py(atom_spec);
   if (p.first) {
      int pos_length = PyObject_Length(position);
      if (pos_length == 3) {
	 PyObject *x_py = PyList_GetItem(position, 0);
	 PyObject *y_py = PyList_GetItem(position, 1);
	 PyObject *z_py = PyList_GetItem(position, 2);
	 double x = PyFloat_AsDouble(x_py);
	 double y = PyFloat_AsDouble(y_py);
	 double z = PyFloat_AsDouble(z_py);
	 clipper::Coord_orth pt(x,y,z);
	 graphics_info_t g;
	 g.add_target_position_restraint_for_intermediate_atom(p.second, pt); // refines after added

	 retval = Py_True;
      }
   }

   Py_INCREF(retval);
   return retval;
}
#endif

// and the multiple-atom version of that (so that they can be applied at the same time)
#ifdef USE_PYTHON
PyObject *add_target_position_restraints_for_intermediate_atoms_py(PyObject *atom_spec_position_list) {

   PyObject *ret_val = Py_False; // not changed by function at the moment

   if (PyList_Check(atom_spec_position_list)) {
      graphics_info_t g;
      if (false) // debug
	 std::cout << "add_target_position_restraints_for_intermediate_atoms_py processing "
		   << PyBytes_AS_STRING(PyUnicode_AsUTF8String(display_python(atom_spec_position_list))) << std::endl;
      std::vector<std::pair<coot::atom_spec_t, clipper::Coord_orth> > atom_spec_position_vec;
      unsigned int len = PyObject_Length(atom_spec_position_list);
      for (std::size_t i=0; i<len; i++) {
	 PyObject *list_item = PyList_GetItem(atom_spec_position_list, i);
	 PyObject *atom_spec_py = PyList_GetItem(list_item, 0);
	 PyObject *position_py  = PyList_GetItem(list_item, 1);
	 std::pair<bool, coot::atom_spec_t> p = make_atom_spec_py(atom_spec_py);
	 if (p.first) {
	    int pos_length = PyObject_Length(position_py);
	    if (PyList_Check(position_py)) {
	       if (pos_length == 3) {
		  PyObject *x_py = PyList_GetItem(position_py, 0);
		  PyObject *y_py = PyList_GetItem(position_py, 1);
		  PyObject *z_py = PyList_GetItem(position_py, 2);
		  double x = PyFloat_AsDouble(x_py);
		  double y = PyFloat_AsDouble(y_py);
		  double z = PyFloat_AsDouble(z_py);
		  clipper::Coord_orth pt(x,y,z);
		  std::pair<coot::atom_spec_t, clipper::Coord_orth> pp(p.second, pt);
		  atom_spec_position_vec.push_back(pp);
	       }
	    } else {
	       PyObject *ds = display_python(position_py);
	       if (ds)
		  std::cout << "WARNING:: position is not a list "
			    << PyUnicode_AsUTF8String(ds) << std::endl;
	       else
		  std::cout << "WARNING:: position is not a list - null from display_python() with input"
			    << position_py << std::endl;
	    }
	 }
      }
      g.add_target_position_restraints_for_intermediate_atoms(atom_spec_position_vec); // refines after added

   } else {
      std::cout << "WARNING:: add_target_position_restraints_for_intermediate_atoms_py() Not a list" << std::endl;
   }
   Py_INCREF(ret_val);
   return ret_val;
}
#endif // USE_PYTHON



#ifdef USE_PYTHON
PyObject *mark_atom_as_fixed_py(int imol, PyObject *atom_spec, int state) {
   PyObject *retval = Py_False;
   std::pair<bool, coot::atom_spec_t> p = make_atom_spec_py(atom_spec);
   if (p.first) {
      graphics_info_t::mark_atom_as_fixed(imol, p.second, state);
      graphics_draw();
      retval = Py_True; // Shall we return True if atom got marked?
   }
   Py_INCREF(retval);
   return retval;
}
#endif // USE_PYTHON



/*! \brief clear all fixed atoms */
void clear_all_fixed_atoms(int imol) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].clear_all_fixed_atoms();
      graphics_draw();
   }
   graphics_info_t g;
   g.setup_draw_for_anchored_atom_markers(); // update the instancing buffer
}


void clear_fixed_atoms_all() {

   for (int i=0; i<graphics_info_t::n_molecules(); i++) {
      if (is_valid_model_molecule(i)) {
	 clear_all_fixed_atoms(i);
      }
   }
   graphics_draw();
}


// ipick is on/off, is_unpick is when we are picking a fixed atom to
// be unfixed.
//
void setup_fixed_atom_pick(short int ipick, short int is_unpick) {

   graphics_info_t g;
   if (ipick == 0) {
      graphics_info_t::in_fixed_atom_define = coot::FIXED_ATOM_NO_PICK;
   } else {
      g.pick_cursor_maybe();
      if (is_unpick) {
	 graphics_info_t::in_fixed_atom_define = coot::FIXED_ATOM_UNFIX;
      } else {
	 graphics_info_t::in_fixed_atom_define = coot::FIXED_ATOM_FIX;
      }
   }
}



#ifdef USE_GUILE
// return e.g (list 1 "C" "D")
//
SCM merge_molecules(SCM add_molecules, int imol) {
   SCM r = SCM_BOOL_F;

   std::vector<int> vam;
   SCM l_length_scm = scm_length(add_molecules);

   int l_length = scm_to_int(l_length_scm);
   for (int i=0; i<l_length; i++) {
      SCM le = scm_list_ref(add_molecules, scm_from_int(i));
      int ii = scm_to_int(le);
      vam.push_back(ii);
   }

   std::pair<int, std::vector<merge_molecule_results_info_t> > v =
      merge_molecules_by_vector(vam, imol);

   SCM v_scm = SCM_EOL;
   for (std::size_t i=0; i<v.second.size(); i++) {
      if (v.second[i].is_chain) {
	 SCM item_scm = scm_from_locale_string(v.second[i].chain_id.c_str());
	 v_scm = scm_cons(item_scm, v_scm);
      } else {
	 SCM item_scm = residue_spec_to_scm(v.second[i].spec);
	 v_scm = scm_cons(item_scm, v_scm);
      }
   }
   v_scm = scm_reverse(v_scm);

   r = SCM_EOL;
   r = scm_cons(v_scm, r);
   r = scm_cons(scm_from_int(v.first), r);

   return r;
}
#endif

#ifdef USE_PYTHON
// some python version of the merge_molecules()
// return e.g [1,"C","D"]
//
PyObject *merge_molecules_py(PyObject *add_molecules, int imol) {

   PyObject *r = Py_False;
   PyObject *le;

   std::vector<int> vam;

   int l_length = PyObject_Length(add_molecules);
   for (int i=0; i<l_length; i++) {
      le = PyList_GetItem(add_molecules, i);
//      int ii = (int)le;
      int ii = PyLong_AsLong(le);
      vam.push_back(ii);
   }

   std::pair<int, std::vector<merge_molecule_results_info_t> > v =
      merge_molecules_by_vector(vam, imol);

   r = PyList_New(v.second.size() + 1);
   PyList_SetItem(r, 0, PyLong_FromLong(v.first));

   // 20180529-PE return a residue spec on merging if we can, else return a
   // chain id as before.
   //
   for (unsigned int i=0; i<v.second.size(); i++) {
      if (v.second[i].is_chain) {
	 PyObject *o = myPyString_FromString(v.second[i].chain_id.c_str());
	 PyList_SetItem(r, i+1, o);
      } else {
	 PyObject *o = residue_spec_to_py(v.second[i].spec);
	 PyList_SetItem(r, i+1, o);
      }
   }

   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif

std::pair<int, std::vector<merge_molecule_results_info_t> >
merge_molecules_by_vector(const std::vector<int> &add_molecules, int imol) {

   std::pair<int, std::vector<merge_molecule_results_info_t> >  merged_info;

   std::vector<atom_selection_container_t> add_molecules_at_sels;
   if (is_valid_model_molecule(imol)) {
      for (unsigned int i=0; i<add_molecules.size(); i++) {
	 if (is_valid_model_molecule(add_molecules[i])) {
	    if (add_molecules[i] != imol) {
	       add_molecules_at_sels.push_back(graphics_info_t::molecules[add_molecules[i]].atom_sel);
	       set_mol_displayed(add_molecules[i], 0);
	       set_mol_active(add_molecules[i], 0);
	    }
	 }
      }
   }
   if (add_molecules_at_sels.size() > 0) {
      merged_info = graphics_info_t::molecules[imol].merge_molecules(add_molecules_at_sels);
   }

   if (graphics_info_t::use_graphics_interface_flag) {
      graphics_info_t g;
      g.update_go_to_atom_window_on_changed_mol(imol);
   }
   return merged_info;
}


#ifdef USE_GUILE
void set_merge_molecules_ligand_spec_scm(SCM ligand_spec_scm) {

   coot::residue_spec_t spec = residue_spec_from_scm(ligand_spec_scm);
   graphics_info_t g;
   g.set_merge_molecules_ligand_spec(spec);
}
#endif // USE_GUILE

#ifdef USE_PYTHON
void set_merge_molecules_ligand_spec_py(PyObject *ligand_spec_py) {

   coot::residue_spec_t spec = residue_spec_from_py(ligand_spec_py);
   graphics_info_t g;
   g.set_merge_molecules_ligand_spec(spec);
}
#endif // USE_PYTHON



/*  ----------------------------------------------------------------------- */
/*                         construct a molecule and update                  */
/*  ----------------------------------------------------------------------- */

#ifdef USE_GUILE
int clear_and_update_molecule(int molecule_number, SCM molecule_expression) {

   int state = 0;
   if (is_valid_model_molecule(molecule_number)) {

      mmdb::Manager *mol = mmdb_manager_from_scheme_expression(molecule_expression);
      if (mol) {
	 state = 1;
	 graphics_info_t::molecules[molecule_number].replace_molecule(mol);
	 graphics_draw();
	 graphics_info_t g;
	 g.update_validation(molecule_number);
      }
   } else {
      std::cout << "WARNING:: " << molecule_number << " is not a valid model molecule"
		<< std::endl;
   }
   return state;
}
#endif // USE_GUILE

#ifdef USE_PYTHON
int clear_and_update_molecule_py(int molecule_number, PyObject *molecule_expression) {

   int state = 0;
   if (is_valid_model_molecule(molecule_number)) {

      std::deque<mmdb::Model *> model_list = mmdb_models_from_python_expression(molecule_expression);
      if (!model_list.empty()) {
         state = 1;
         graphics_info_t::molecules[molecule_number].replace_models(model_list);
         graphics_info_t g;
	 g.update_validation(molecule_number);
         graphics_draw();
      }
   }
   return state;
}
#endif // USE_GUILE

#ifdef USE_GUILE
// Return a molecule number, -1 on error.
int add_molecule(SCM molecule_expression, const char *name) {

   int imol = -1;
   mmdb::Manager *mol =
      mmdb_manager_from_scheme_expression(molecule_expression);
   if (mol) {
      imol = graphics_info_t::create_molecule();
      atom_selection_container_t asc = make_asc(mol);
      graphics_info_t g;
      g.molecules[imol].install_model(imol, asc, g.Geom_p(), name, 1);
      graphics_draw();
   } else {
      std::cout << "WARNING:: bad format, no molecule created"
		<< std::endl;
   }
   return imol;
}
#endif // USE_GUILE

#ifdef USE_PYTHON
// Return a molecule number, -1 on error.
int add_molecule_py(PyObject *molecule_expression, const char *name) {

   int imol = -1;
   mmdb::Manager *mol =
      mmdb_manager_from_python_expression(molecule_expression);
   if (mol) {
      imol = graphics_info_t::create_molecule();
      atom_selection_container_t asc = make_asc(mol);
      graphics_info_t g;
      g.molecules[imol].install_model(imol, asc, g.Geom_p(), name, 1);
      graphics_draw();
   } else {
      std::cout << "WARNING:: bad format, no molecule created"
                << std::endl;
   }
   return imol;
}
#endif // USE_PYTHON





void change_chain_id(int imol, const char *from_chain_id, const char *to_chain_id,
		     short int use_res_range_flag, int from_resno, int to_resno) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      std::pair<int, std::string> r =
	 graphics_info_t::molecules[imol].change_chain_id(from_chain_id,
							  to_chain_id,
							  use_res_range_flag,
							  from_resno,
							  to_resno);
      graphics_draw();
      g.update_go_to_atom_window_on_changed_mol(imol);
      g.update_validation(imol);
   }
}

#ifdef USE_GUILE
SCM change_chain_id_with_result_scm(int imol, const char *from_chain_id, const char *to_chain_id,
				    short int use_res_range_flag, int from_resno, int to_resno){

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      std::pair<int, std::string> p =
	 g.molecules[imol].change_chain_id(from_chain_id,
					   to_chain_id,
					   use_res_range_flag,
					   from_resno,
					   to_resno);
      graphics_draw();
      g.update_go_to_atom_window_on_changed_mol(imol);
      g.update_validation_graphs(imol);
      r = SCM_EOL;
      r = scm_cons(scm_from_locale_string(p.second.c_str()), r);
      r = scm_cons(scm_from_int(p.first), r);
   }
   return r;
}
#endif // USE_GUILE


#ifdef USE_PYTHON
PyObject *change_chain_id_with_result_py(int imol, const char *from_chain_id, const char *to_chain_id,
					 short int use_res_range_flag, int from_resno, int to_resno){

   PyObject *v = Py_False;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      std::pair<int, std::string> r =
	 g.molecules[imol].change_chain_id(from_chain_id,
					   to_chain_id,
					   use_res_range_flag,
					   from_resno,
					   to_resno);

      graphics_draw();
      g.update_go_to_atom_window_on_changed_mol(imol);
      g.update_validation(imol);
      v = PyList_New(2);
      PyList_SetItem(v, 0, PyLong_FromLong(r.first));
      PyList_SetItem(v, 1, myPyString_FromString(r.second.c_str()));
   }
   return v;
}
#endif // USE_PYTHON




/*! \brief set way nomenclature errors should be handled on reading
  coordinates.

  mode should be "auto-correct", "ignore", "prompt".  The
  default is "prompt" */
void set_nomenclature_errors_on_read(const char *mode) {

   std::string m(mode);
   if (m == "auto-correct")
      graphics_info_t::nomenclature_errors_mode = coot::AUTO_CORRECT;
   if (m == "ignore")
      graphics_info_t::nomenclature_errors_mode = coot::IGNORE;
   if (m == "prompt")
      graphics_info_t::nomenclature_errors_mode = coot::PROMPT;

}


/*  ----------------------------------------------------------------------- */
/*               Nomenclature Errors                                        */
/*  ----------------------------------------------------------------------- */
/* return the number of resides altered. */
int fix_nomenclature_errors(int imol) {

   int ifixed = 0;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      std::vector<mmdb::Residue *> vr =
	 graphics_info_t::molecules[imol].fix_nomenclature_errors(g.Geom_p());
      ifixed = vr.size();
      g.update_validation(imol);
      graphics_draw();
   }
   // update geometry graphs (not least rotamer graph).
   // but we have no intermediate atoms...

   return ifixed;

}

// the residue type and the spec.
//
std::vector<std::pair<std::string, coot::residue_spec_t> >
list_nomenclature_errors(int imol) {

   std::vector<std::pair<std::string, coot::residue_spec_t> > r;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      r = g.molecules[imol].list_nomenclature_errors(g.Geom_p());
   }
   return r;
}

#ifdef USE_GUILE
SCM list_nomenclature_errors_scm(int imol) {

   std::vector<std::pair<std::string, coot::residue_spec_t> > v = list_nomenclature_errors(imol);
   SCM r = SCM_EOL;
   if (v.size()) {
      for(int i=v.size()-1; i>=0; i--) {
	 r = scm_cons(residue_spec_to_scm(v[i].second), r);
      }
   }
   return r;

}
#endif // USE_GUILE


#ifdef USE_PYTHON
PyObject *list_nomenclature_errors_py(int imol) {

   PyObject *r = PyList_New(0);
   std::vector<std::pair<std::string, coot::residue_spec_t> > v = list_nomenclature_errors(imol);
   if (v.size()) {
      r = PyList_New(v.size());
      for (unsigned int i=0; i<v.size(); i++) {
	 PyList_SetItem(r, i, residue_spec_to_py(v[i].second));
      }
   }
   return r;
}
#endif // USE_PYTHON


#ifdef USE_GUILE

SCM missing_atom_info_scm(int imol) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      r = SCM_EOL;
      graphics_info_t g;
      bool missing_hydrogens_flag = 0;
      coot::util::missing_atom_info m_i_info =
	 g.molecules[imol].missing_atoms(missing_hydrogens_flag, g.Geom_p());
      for (unsigned int i=0; i<m_i_info.residues_with_missing_atoms.size(); i++) {
	 mmdb::Residue *residue_p = m_i_info.residues_with_missing_atoms[i];
	 int resno                = residue_p->GetSeqNum();
	 std::string chain_id     = residue_p->GetChainID();
	 std::string residue_type = residue_p->GetResName();
	 std::string inscode      = residue_p->GetInsCode();
	 std::string altconf("");
	 SCM l = SCM_EOL;
	 l = scm_cons(scm_from_locale_string(inscode.c_str()), l);
	 l = scm_cons(scm_from_int(resno), l);
	 l = scm_cons(scm_from_locale_string(chain_id.c_str()), l);
	 r = scm_cons(l, r);

         std::map<mmdb::Residue *, std::vector<std::string> >::const_iterator it;
         it = m_i_info.residue_missing_atom_names_map.find(residue_p);
         if (it != m_i_info.residue_missing_atom_names_map.end()) {
            const std::vector<std::string> &missing_atom_names = it->second;
            if (! missing_atom_names.empty()) {
               std::cout << "INFO:: residue " << coot::residue_spec_t(residue_p) << " has missing atoms ";
               for (unsigned int iat=0; iat<missing_atom_names.size(); iat++)
                  std::cout << single_quote(missing_atom_names[iat]) << " ";
               std::cout << std::endl;
            }
         }
      }
      r = scm_reverse(r);
   }
   return r;
}
#endif // USE_GUILE

#ifdef USE_PYTHON

PyObject *missing_atom_info_py(int imol) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      r = PyList_New(0);
      graphics_info_t g;
      short int missing_hydrogens_flag = 0;
      coot::util::missing_atom_info m_i_info =
	 g.molecules[imol].missing_atoms(missing_hydrogens_flag, g.Geom_p());
      for (unsigned int i=0; i<m_i_info.residues_with_missing_atoms.size(); i++) {
	 int resno =  m_i_info.residues_with_missing_atoms[i]->GetSeqNum();
	 std::string chain_id = m_i_info.residues_with_missing_atoms[i]->GetChainID();
	 std::string residue_type = m_i_info.residues_with_missing_atoms[i]->GetResName();
	 std::string inscode = m_i_info.residues_with_missing_atoms[i]->GetInsCode();
	 std::string altconf("");
	 PyObject *l = PyList_New(0);
	 PyList_Append(l, myPyString_FromString(chain_id.c_str()));
	 PyList_Append(l, PyLong_FromLong(resno));
	 PyList_Append(l, myPyString_FromString(inscode.c_str()));
	 PyList_Append(r, l);
	 Py_XDECREF(l);
      }
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif // USE_PYTHON

void
copy_chain(int imol, const char *from_chain, const char *to_chain) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].copy_chain(std::string(from_chain),
						  std::string(to_chain));
      graphics_draw();
   }
}


void set_place_helix_here_fudge_factor(float ff) {
   graphics_info_t::place_helix_here_fudge_factor = ff;
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
   if (is_valid_map_molecule(imol_map)) {
      // perhaps we should add to the constructor the user-setable
      // g.helix_placement_cb_density_descrimination_ratio
      // defaults to 1.0.
      int imol = -1;
      const clipper::Xmap<float> &xmap = graphics_info_t::molecules[imol_map].xmap;
      coot::helix_placement p(xmap);

      // test
      std::vector<float> dv_test(100);
      for (unsigned int i=0; i<100; i++) {
	 dv_test[i] = 100-i;
      }

      std::vector<float> dv = coot::util::density_map_points_in_sphere(pt, 20, xmap);
      float f_iqr = coot::util::interquartile_range(dv);
      float iqr_2 = 0.5*f_iqr;

      float min_density_limit = iqr_2;
      min_density_limit = iqr_2 * g.place_helix_here_fudge_factor; // 20171117
      float high_density_turning_point = f_iqr * 4 * g.place_helix_here_fudge_factor;
      std::cout << "DEBUG:: choosing map min_density limit: " << min_density_limit << std::endl;
      std::cout << "DEBUG:: choosing map high_density_turning_point limit: " << high_density_turning_point
		<< std::endl;

      float map_rmsd = g.molecules[imol_map].map_sigma();
      float bf = g.default_new_atoms_b_factor;

      // 20171117
      // min-density_limit min_density_limit_for_trim should be higher for cryo-EM maps
      //
      coot::helix_placement_info_t n =
	 p.place_alpha_helix_near_kc_version(pt, 20, min_density_limit, high_density_turning_point, bf, map_rmsd);
       if (! n.success) {
	  n = p.place_alpha_helix_near_kc_version(pt, 9, min_density_limit, high_density_turning_point, bf, map_rmsd);
       }

       if (n.success) {
	  atom_selection_container_t asc = make_asc(n.mol[0].pcmmdbmanager());

	  imol = g.create_molecule();
	  std::string mol_name = "Helix-";
	  mol_name+= coot::util::int_to_string(imol);
	  graphics_info_t::molecules[imol].install_model(imol, asc, g.Geom_p(), mol_name, 1);

	  if (n.mol.size() > 1) {
	     atom_selection_container_t asc2 = make_asc(n.mol[1].pcmmdbmanager());
	     imol = g.create_molecule();
	     mol_name = "Reverse-Helix-";
	     mol_name+= coot::util::int_to_string(imol);
	     graphics_info_t::molecules[imol].install_model(imol, asc2, g.Geom_p(), mol_name, 1);
	  }

	  if (g.go_to_atom_window) {
	     g.set_go_to_atom_molecule(imol);
	     g.update_go_to_atom_window_on_new_mol();
	  } else {
	     g.set_go_to_atom_molecule(imol);
	  }

	  g.add_status_bar_text("Helix added");
       } else {
	  std::cout << "Helix addition failure: message: " << n.failure_message << "\n";
	  g.add_status_bar_text(n.failure_message);
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
      g.add_status_bar_text("You need to set the map to fit against");
      g.show_select_map_frame();
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

      coot::helix_placement p(graphics_info_t::molecules[imol_map].xmap);
      float s = graphics_info_t::molecules[imol_map].map_sigma();
      float ff = graphics_info_t::place_helix_here_fudge_factor;
      if (graphics_info_t::molecules[imol_map].is_EM_map())
	 ff = 3.0;
      float s_with_ff = s * ff;
      coot::helix_placement_info_t si =
	 p.place_strand(pt, n_residues, n_sample_strands, s_with_ff);
      if (si.success) {
	 // nice to refine the fragment here, but the interface
	 // doesn't work that way, so put the refinement after the
	 // molecule has been accepted.
	 atom_selection_container_t asc = make_asc(si.mol[0].pcmmdbmanager());
	 imol = g.create_molecule();
	 std::string mol_name = "Strand-";
	 mol_name+= coot::util::int_to_string(imol);
	 graphics_info_t::molecules[imol].install_model(imol, asc, g.Geom_p(), mol_name, 1);
	 g.add_status_bar_text("Strand added");

	 // Now refine.
	 coot::minimol::zone_info_t zi = si.mol[0].zone_info();
	 if (zi.is_simple_zone) {
	    int save_rirf = g.refinement_immediate_replacement_flag;

	    coot::pseudo_restraint_bond_type save_pseudos = g.pseudo_bonds_type;
	    g.pseudo_bonds_type = coot::STRAND_PSEUDO_BONDS;
	    g.refinement_immediate_replacement_flag = 1;
	    g.refine_residue_range(imol, zi.chain_id, zi.chain_id, zi.resno_1, "",
				   zi.resno_2, "", "", 0);
	    accept_regularizement();
	    g.pseudo_bonds_type = save_pseudos;

	    g.refinement_immediate_replacement_flag = save_rirf;
	 }
      } else {
	 std::cout << "Strand addition failure: message: " << si.failure_message << "\n";
	 g.add_status_bar_text(si.failure_message);
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
      g.add_status_bar_text("You need to set the map to fit against");
      g.show_select_map_frame();
      return -1;
   }
}


/*  ----------------------------------------------------------------------- */
/*                  Autofind Helices                                        */
/*  ----------------------------------------------------------------------- */
/* Find secondary structure in the current map.
   Add to a molecule called "Helices", create it if
   needed. */
int find_helices() {
  return find_secondary_structure_local( 1, 7, FIND_SECSTRUC_NORMAL,
					 0, 0, FIND_SECSTRUC_NORMAL, 0.0 );
}


/*  ----------------------------------------------------------------------- */
/*                  Autofind Strands                                        */
/*  ----------------------------------------------------------------------- */
/* Find secondary structure in the current map.
   Add to a molecule called "Strands", create it if
   needed. */
int find_strands() {
  return find_secondary_structure_local( 0, 0, FIND_SECSTRUC_STRICT,
					 1, 5, FIND_SECSTRUC_STRICT, 0.0 );
}


/*  ----------------------------------------------------------------------- */
/*                  Autofind Secondary Structure                            */
/*  ----------------------------------------------------------------------- */
/* Find secondary structure in the current map.
   Add to a molecule called "SecStruc", create it if
   needed. */
int find_secondary_structure(
    short int use_helix,  int helix_length,  int helix_target,
    short int use_strand, int strand_length, int strand_target )
{
  return find_secondary_structure_local(
      use_helix,  helix_length,  helix_target,
      use_strand, strand_length, strand_target,
      0.0 );
}

/*  ----------------------------------------------------------------------- */
/*                  Autofind Secondary Structure                            */
/*  ----------------------------------------------------------------------- */
/* Find secondary structure in the current map.
   Add to a molecule called "SecStruc", create it if
   needed. */
int find_secondary_structure_local(
    short int use_helix,  int helix_length,  int helix_target,
    short int use_strand, int strand_length, int strand_target,
    float radius )
{
   graphics_info_t g;
   clipper::Coord_orth pt(g.RotationCentre_x(),
			  g.RotationCentre_y(),
			  g.RotationCentre_z());
   int imol_map = g.Imol_Refinement_Map();
   if (imol_map != -1) {
      using coot::SSfind;
      coot::SSfind::SSTYPE ta[] = { SSfind::ALPHA3, SSfind::ALPHA3S,
				    SSfind::ALPHA2, SSfind::ALPHA4 };
      coot::SSfind::SSTYPE tb[] = { SSfind::BETA3,  SSfind::BETA3S,
				    SSfind::BETA2,  SSfind::BETA4  };
      int imol = -1;
      std::vector<SSfind::Target> tgtvec;
      if ( use_helix )
	tgtvec.push_back( SSfind::Target(ta[helix_target %4],helix_length ) );
      if ( use_strand )
	tgtvec.push_back( SSfind::Target(tb[strand_target%4],strand_length) );
      coot::fast_secondary_structure_search ssfind;
      ssfind( graphics_info_t::molecules[imol_map].xmap, pt, 0.0, tgtvec );
      if (ssfind.success) {
	 atom_selection_container_t asc = make_asc(ssfind.mol.pcmmdbmanager());
	 imol = g.create_molecule();
	 graphics_info_t::molecules[imol].install_model(imol,asc,g.Geom_p(),"SecStruc",1);
	 g.molecules[imol].ca_representation(true);
	 if (g.go_to_atom_window) {
	    g.set_go_to_atom_molecule(imol);
	    g.update_go_to_atom_window_on_new_mol();
	 } else {
	    g.set_go_to_atom_molecule(imol);
	 }
	 g.add_status_bar_text("Secondary structure added");
      } else {
	 std::cout << "No secondary structure found\n";
	 g.add_status_bar_text("No secondary structure found" );
      }
      std::vector<std::string> command_strings;
      command_strings.resize(0);
      command_strings.push_back("find-secondary-structure");
      add_to_history(command_strings);
      graphics_draw();
      return imol;
   } else {
      std::cout << " You need to set the map to fit against\n";
      g.add_status_bar_text("You need to set the map to fit against");
      g.show_select_map_frame();
      return -1;
   }
}


/*  ----------------------------------------------------------------------- */
/*                  Autofind Nucleic acid chains                            */
/*  ----------------------------------------------------------------------- */
/* Find secondary structure local to current view in the current map.
   Add to a molecule called "NuclAcid", create it if needed. */
int find_nucleic_acids_local(float radius) {

   // check for data file
   std::string nafile;
   const char *cp = getenv("COOT_PREFIX");
   if (cp) nafile = std::string(cp) + "/share/coot/nautilus_lib.pdb";
   else    nafile = std::string(PKGDATADIR) + "/nautilus_lib.pdb";
   if (!coot::file_exists(nafile)) {
      std::cout << "Ooops! Can't find nautilus data! - fail" << std::endl;
      return -1;
   }

   // check for map
   graphics_info_t g;
   clipper::Coord_orth pt(g.RotationCentre_x(),
			  g.RotationCentre_y(),
			  g.RotationCentre_z());
   int imol_map = g.Imol_Refinement_Map();
   if (imol_map == -1) {
      std::cout << " You need to set the map to fit against\n";
      g.add_status_bar_text("You need to set the map to fit against");
      g.show_select_map_frame();
      return -1;
   }

   // FIND OR CREATE A MODEL FOR THE NAs
   int imol = -1;
   mmdb::Manager *mol;
   for ( int i = 0 ; i < graphics_n_molecules(); i++)
     if ( graphics_info_t::molecules[i].has_model() &&
	  graphics_info_t::molecules[i].name_ == "NuclAcid" ) {
       imol = i;
       mol = graphics_info_t::molecules[i].atom_sel.mol;
       break;
     }
   if ( imol < 0 ) {
     imol = g.create_molecule();
     mol = new mmdb::Manager;
     graphics_info_t::molecules[imol].install_model( imol, mol, g.Geom_p(), "NuclAcid", 1 );
   }

   // build the model
   Coot_nucleic_acid_build nafind( nafile );
   bool success = nafind.build( mol, graphics_info_t::molecules[imol_map].xmap, pt, radius );
   graphics_info_t::molecules[imol].update_molecule_after_additions();

   if (success) {
   	 if (g.go_to_atom_window) {
   	    g.set_go_to_atom_molecule(imol);
   	    g.update_go_to_atom_window_on_new_mol();
   	 } else {
   	    g.set_go_to_atom_molecule(imol);
   	 }
   	 std::cout << "Nucleic acids found" << std::endl;;
   	 g.add_status_bar_text("Nucleic acids added");
   } else {
   	 std::cout << "No nucleic acids found\n";
   	 g.add_status_bar_text("No nucleic acids found" );
   }
   std::vector<std::string> command_strings;
   command_strings.resize(0);
   command_strings.push_back("find-nucleic-acids-local");
   add_to_history(command_strings);
   graphics_draw();
   return imol;
}


/*  ----------------------------------------------------------------------- */
/*                  FFfearing                                               */
/*  ----------------------------------------------------------------------- */


// return the new model number
//
int fffear_search(int imol_model, int imol_map) {

   float angular_resolution = graphics_info_t::fffear_angular_resolution;
   int imol_new = -1;
   if (!is_valid_model_molecule(imol_model)) {
      std::cout << "WARNING:: this is not a valid model: " << imol_model << std::endl;
   } else {
      if (!is_valid_map_molecule(imol_map)) {
	 std::cout << "WARNING:: this is not a valid map: " << imol_map << std::endl;
      } else {
	 coot::util::fffear_search f(graphics_info_t::molecules[imol_model].atom_sel.mol,
				     graphics_info_t::molecules[imol_model].atom_sel.SelectionHandle,
				     graphics_info_t::molecules[imol_map].xmap,
				     angular_resolution);

	 imol_new = graphics_info_t::create_molecule();
	 std::string name("FFFear search results");
	 bool is_em_flag = graphics_info_t::molecules[imol_map].is_EM_map();
	 graphics_info_t::molecules[imol_new].install_new_map(f.get_results_map(), name, is_em_flag);

	 std::vector<std::pair<float, clipper::RTop_orth> > p = f.scored_orientations();
	 if (p.size() > 0) {
	    // install new molecule(s) that has been rotated.
	 }
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

void rigid_body_refine_zone(int imol, const char *chain_id, int resno_start, int resno_end) {

   graphics_info_t g;
   std::string altconf = ""; // should be passed?

   // need to set graphics_info's residue_range_atom_index_1,
   // residue_range_atom_index_2, imol_rigid_body_refine

   if (imol < g.n_molecules()) {
      if (g.molecules[imol].has_model()) {
	 g.imol_rigid_body_refine = imol;

	 g.set_residue_range_refine_atoms(std::string(chain_id),
					  resno_start, resno_end,
					  altconf,
					  imol);
	 g.execute_rigid_body_refine(0);
      }
   }
}


void
rigid_body_refine_by_atom_selection(int imol, const char *atom_selection_string) {

   graphics_info_t g;
   int imol_ref_map = g.Imol_Refinement_Map();
   if (is_valid_map_molecule(imol_ref_map)) {
      if (is_valid_model_molecule(imol)) {
	 bool mask_waters_flag = 0;

	 // so the bulk of this function is to generate
	 // mol_without_moving_zone and range_mol from
	 // atom_selection_string.
	 //
	 // Let's make a (ligand) utility function for this.
	 //
	 bool fill_mask = true;
	 mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;

	 if (false)
	    std::cout << "debug in rigid_body_refine_by_atom_selection() start: here "
                      << "UDDAtomIndexHandle is "
		      << g.molecules[imol].atom_sel.UDDAtomIndexHandle << std::endl;

	 std::string atom_selection_str(atom_selection_string);
	 // first is the atoms of the mask (not in the selection
	 // second is the atoms of the selection
	 std::pair<coot::minimol::molecule, coot::minimol::molecule> p =
	    coot::make_mols_from_atom_selection_string(mol, atom_selection_str, fill_mask);

	 g.imol_rigid_body_refine = imol;
	 g.rigid_body_fit(p.first,   // without selected region.
			  p.second,  // selected region.
			  imol_ref_map,
			  mask_waters_flag);
      } else {
	 std::cout << "WARNING:: model molecule " << imol << " is not valid " << std::endl;
      }
   } else {
      std::cout << "WARNING:: refinement map not defined. " << std::endl;
   }
}


#ifdef USE_GUILE
/* ! \brief rigid body refine using residue ranges.  residue_ranges is
   a list of residue ranges.  A residue range is (list chain-id
   resno-start resno-end). */
SCM
rigid_body_refine_by_residue_ranges_scm(int imol, SCM residue_ranges) {

   SCM ret_val = SCM_BOOL_F;
   std::vector<coot::high_res_residue_range_t> res_ranges;
   if (scm_is_true(scm_list_p(residue_ranges))) {
      SCM rr_length_scm = scm_length(residue_ranges);
      int rr_length = scm_to_int(rr_length_scm);
      if (rr_length > 0) {
	 for (int irange=0; irange<rr_length; irange++) {
	    SCM range_scm = scm_list_ref(residue_ranges, scm_from_int(irange));
	    if (scm_is_true(scm_list_p(range_scm))) {
	       SCM range_length_scm = scm_length(range_scm);
	       int range_length = scm_to_int(range_length_scm);
	       if (range_length == 3) {
		  SCM chain_id_scm    = scm_list_ref(range_scm, scm_from_int(0));
		  SCM resno_start_scm = scm_list_ref(range_scm, scm_from_int(1));
		  SCM resno_end_scm   = scm_list_ref(range_scm, scm_from_int(2));
		  if (scm_is_string(chain_id_scm)) {
		     std::string chain_id = scm_to_locale_string(chain_id_scm);
		     if (scm_is_true(scm_number_p(resno_start_scm))) {
			int resno_start = scm_to_int(resno_start_scm);
			if (scm_is_true(scm_number_p(resno_end_scm))) {
			   int resno_end = scm_to_int(resno_end_scm);
			   // recall that mmdb does crazy things with
			   // the residue selection if the second
			   // residue is before the first residue in
			   // the chain.  So, check and swap if needed.
			   if (resno_end < resno_start)
			      std::swap<int>(resno_start, resno_end);
			   coot::high_res_residue_range_t rr(chain_id, resno_start, resno_end);
			   res_ranges.push_back(rr);
			}
		     }
		  }
	       }
	    }
	 }
	 int status = rigid_body_fit_with_residue_ranges(imol, res_ranges); // test for res_ranges
	                                                                    // length in here.
	 if (status == 1) // good
	    ret_val = SCM_BOOL_T;
      } else {
	 std::cout << "incomprehensible input to rigid_body_refine_by_residue_ranges_scm"
		   << " null list" << std::endl;
      }
   } else {
      std::cout << "incomprehensible input to rigid_body_refine_by_residue_ranges_scm"
		<< " not a list" << std::endl;
   }
   return ret_val;
}
#endif /* USE_GUILE */

#ifdef USE_PYTHON
/* ! \brief rigid body refine using residue ranges.  residue_ranges is
   a list of residue ranges.  A residue range is [chain-id,
   resno-start, resno-end]. */
PyObject *
rigid_body_refine_by_residue_ranges_py(int imol, PyObject *residue_ranges) {

   PyObject *ret_val = Py_False;
   std::vector<coot::high_res_residue_range_t> res_ranges;
   if (PyList_Check(residue_ranges)) {
      int rr_length = PyObject_Length(residue_ranges);
      if (rr_length > 0) {
	 for (int irange=0; irange<rr_length; irange++) {
	   PyObject *range_py = PyList_GetItem(residue_ranges, irange);
	   if (PyList_Check(range_py)) {
	     int range_length = PyObject_Length(range_py);
	     if (range_length == 3) {
	       PyObject *chain_id_py    = PyList_GetItem(range_py, 0);
	       PyObject *resno_start_py = PyList_GetItem(range_py, 1);
	       PyObject *resno_end_py   = PyList_GetItem(range_py, 2);
	       if (PyUnicode_Check(chain_id_py)) {
		 std::string chain_id = PyBytes_AS_STRING(PyUnicode_AsUTF8String(chain_id_py));
		 if (PyLong_Check(resno_start_py)) {
		   int resno_start = PyLong_AsLong(resno_start_py);
		   if (PyLong_Check(resno_end_py)) {
		     int resno_end = PyLong_AsLong(resno_end_py);
		     // recall that mmdb does crazy things with
		     // the residue selection if the second
		     // residue is before the first residue in
		     // the chain.  So, check and swap if needed.
		     if (resno_end < resno_start)
		       std::swap<int>(resno_start, resno_end);
		     coot::high_res_residue_range_t rr(chain_id, resno_start, resno_end);
		     res_ranges.push_back(rr);
		   }
		 }
	       }
	     }
	   }
	 }
	 int status = rigid_body_fit_with_residue_ranges(imol, res_ranges); // test for res_ranges
	                                                                    // length in here.
	 if (status == 1) // good
	   ret_val = Py_True;
      } else {
	std::cout << "incomprehensible input to rigid_body_refine_by_residue_ranges_scm"
		  << " null list" << std::endl;
      }
   } else {
     std::cout << "incomprehensible input to rigid_body_refine_by_residue_ranges_scm"
	       << " not a list" << std::endl;
   }
   if (PyBool_Check(ret_val)) {
     Py_INCREF(ret_val);
   }
   return ret_val;
}
#endif /* USE_PYTHON */


/*  ----------------------------------------------------------------------- */
/*                  rigid body fitting (multiple residue ranges)            */
/*  ----------------------------------------------------------------------- */
int rigid_body_fit_with_residue_ranges(int imol,
					const std::vector<coot::high_res_residue_range_t> &residue_ranges) {

   int success = 0;
   graphics_info_t g;
   int imol_ref_map = g.Imol_Refinement_Map();
   if (is_valid_map_molecule(imol_ref_map)) {
      if (is_valid_model_molecule(imol)) {

	 if (residue_ranges.size()) {
	    bool mask_waters_flag = 0;

	    mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
	    int SelHnd = mol->NewSelection();

	    for (unsigned int ir=0; ir<residue_ranges.size(); ir++) {
	       mol->SelectAtoms(SelHnd, 0,
				residue_ranges[ir].chain_id.c_str(),
				residue_ranges[ir].start_resno, "*",
				residue_ranges[ir].end_resno, "*",
				"*","*","*","*", mmdb::SKEY_OR);
	    }

	    mmdb::Manager *mol_from_selected =
	       coot::util::create_mmdbmanager_from_atom_selection(mol, SelHnd);
	    // atom selection in mol gets inverted by this function:
	    mmdb::Manager *mol_from_non_selected =
	       coot::util::create_mmdbmanager_from_atom_selection(mol, SelHnd, 1);

	    coot::minimol::molecule range_mol  = coot::minimol::molecule(mol_from_selected);
	    coot::minimol::molecule masked_mol = coot::minimol::molecule(mol_from_non_selected);
	    delete mol_from_selected;
	    delete mol_from_non_selected;
	    mol->DeleteSelection(SelHnd);
	    g.imol_rigid_body_refine = imol;
	    success = g.rigid_body_fit(masked_mol,   // without selected region.
				       range_mol,  // selected region.
				       imol_ref_map,
				       mask_waters_flag);
	 }
      }
   }
   return success;
}



void set_secondary_structure_restraints_type(int itype) {

   // Remember that Rama restraints are not secondary structure restraints.

   if (itype == 0)
      graphics_info_t::pseudo_bonds_type = coot::NO_PSEUDO_BONDS;
   if (itype == 1)
      graphics_info_t::pseudo_bonds_type = coot::HELIX_PSEUDO_BONDS;
   if (itype == 2)
      graphics_info_t::pseudo_bonds_type = coot::STRAND_PSEUDO_BONDS;

}

/*! \brief return the secondary structure restraints type */
int secondary_structure_restraints_type() {

   // cast a pseudo_restraint_bond_type to an int
   return graphics_info_t::pseudo_bonds_type;
}


void accept_regularizement() {

   c_accept_moving_atoms();
}

void c_accept_moving_atoms() {

   graphics_info_t g;
   while (g.continue_threaded_refinement_loop)
      std::this_thread::sleep_for(std::chrono::milliseconds(200));
   g.clear_hud_buttons();
   g.accept_moving_atoms();
   g.clear_moving_atoms_object();

}

#ifdef USE_GUILE
SCM accept_moving_atoms_scm() {

   graphics_info_t g;
   while (g.continue_threaded_refinement_loop) {
      std::this_thread::sleep_for(std::chrono::milliseconds(200));
   }
   coot::refinement_results_t rr = g.accept_moving_atoms(); // does a g.clear_up_moving_atoms();
   rr.show();
   g.clear_moving_atoms_object();
   return g.refinement_results_to_scm(rr);
}
#endif

#ifdef USE_PYTHON
PyObject *accept_moving_atoms_py() {

   graphics_info_t g;
   while (g.continue_threaded_refinement_loop) {
      // std::cout << "wait ..." << std::endl;
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
   }

   coot::refinement_results_t rr = g.accept_moving_atoms(); // does a g.clear_up_moving_atoms();
   rr.show();
   if (g.use_graphics_interface_flag) {
      g.clear_hud_buttons();
      g.clear_moving_atoms_object();
   }
   PyObject *o = g.refinement_results_to_py(rr);
   return o;
}
#endif


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

      imol = graphics_info_t::create_molecule();
      mmdb::Manager *mol_orig = graphics_info_t::molecules[imol_orig].atom_sel.mol;
      int SelectionHandle = mol_orig->NewSelection();
      mol_orig->SelectAtoms(SelectionHandle, 0, "*",
			    mmdb::ANY_RES, "*",
			    mmdb::ANY_RES, "*",
			    (char *) residue_type,
			    "*", "*", "*");
      mmdb::Manager *mol =
	 coot::util::create_mmdbmanager_from_atom_selection(mol_orig, SelectionHandle);

      if (mol) {
	 std::string name = "residue type ";
	 name += residue_type;
	 name += " from ";
	 name += graphics_info_t::molecules[imol_orig].name_for_display_manager();
	 atom_selection_container_t asc = make_asc(mol);
	 if (asc.n_selected_atoms > 0) {
	    bool shelx_flag = 0;
	    if (graphics_info_t::molecules[imol_orig].is_from_shelx_ins())
	       shelx_flag = 1;
	    graphics_info_t g;
	    g.molecules[imol].install_model(imol, asc, g.Geom_p(), name, 1, shelx_flag);
	 } else {
            std::cout << "in new_molecule_by_residue_type_selection "
                      << "Something bad happened - No residues selected"
                      << std::endl;
            std::string s = "Oops, failed to select residue type. ";
            s += "No residues selected\n";
            s += "Residue ";
            s += "\"";
            s += residue_type;
            s += "\" does not exist in molecule ";
	    s += coot::util::int_to_string(imol_orig);
	    s += "!?";
            info_dialog(s.c_str());
            imol = -1;
            graphics_info_t::erase_last_molecule();
	 }
      } else {
	 std::cout << "in new_molecule_by_residue_type_selection "
		   << "Something bad happened - null molecule" << std::endl;
	 graphics_info_t::erase_last_molecule();
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

   auto recentre_on_new_fragment = [] (int imol) {
      graphics_info_t g;
      coot::view_info_t this_view(g.view_quaternion, g.RotationCentre(), g.zoom, "");
      float new_zoom = 100.0;
      coot::Cartesian new_rotation_centre = g.molecules[imol].centre_of_molecule();
      coot::view_info_t  new_view(g.view_quaternion, new_rotation_centre, new_zoom, "");
      int nsteps = int(1000.0/g.views_play_speed);
      coot::view_info_t::interpolate(this_view, new_view, nsteps);
   };

   int imol = -1;
   if (is_valid_model_molecule(imol_orig)) {
      imol = graphics_info_t::create_molecule();
      mmdb::Manager *mol_orig = graphics_info_t::molecules[imol_orig].atom_sel.mol;
      int SelectionHandle = mol_orig->NewSelection();

      std::vector<std::string> parts = coot::util::split_string(std::string(atom_selection_str), "||");
      for (const auto &part : parts)
         mol_orig->Select(SelectionHandle, mmdb::STYPE_ATOM, part.c_str(), mmdb::SKEY_OR);

      mmdb::Manager *mol = coot::util::create_mmdbmanager_from_atom_selection(mol_orig, SelectionHandle);

      if (mol) {
	 std::string name = "Atom selection ";
         name += atom_selection_str;
         name += " from ";
	 name += graphics_info_t::molecules[imol_orig].name_for_display_manager();
	 atom_selection_container_t asc = make_asc(mol);
	 if (asc.n_selected_atoms > 0){
	    bool shelx_flag = 0;
	    if (graphics_info_t::molecules[imol_orig].is_from_shelx_ins())
	       shelx_flag = 1;
	    graphics_info_t g;
	    g.molecules[imol].install_model(imol, asc, g.Geom_p(), name, 1, shelx_flag);
	    g.molecules[imol].set_have_unsaved_changes_from_outside();
	    update_go_to_atom_window_on_new_mol();
            recentre_on_new_fragment(imol);
	 } else {
	    std::cout << "in new_molecule_by_atom_selection "
		      << "Something bad happened - No atoms selected"
		      << std::endl;
	    std::string s = "WARNING:: Oops! failed to create fragment.  ";
	    s += "No atoms selected\n";
	    s += "Incorrect atom specifier? ";
	    s += "\"";
	    s += atom_selection_str;
	    s += "\"";
	    info_dialog(s.c_str());
	    imol = -1;
	    graphics_info_t::erase_last_molecule();
	 }
      } else {
	 // mol will (currently) never be null,
	 // create_mmdbmanager_from_atom_selection() always returns a
	 // good mmdb::Manager pointer.
	 std::cout << "in new_molecule_by_atom_selection "
		   << "Something bad happened - null molecule" << std::endl;
	 std::string s = "WARNING:: Oops! failed to create fragment.  ";
	 s += "Incorrect atom specifier?\n";
	 s += "\"";
	 s += atom_selection_str;
	 s += "\"";
	 info_dialog(s.c_str());
	 imol = -1;
	 graphics_info_t::erase_last_molecule();
      }
      mol_orig->DeleteSelection(SelectionHandle);
      graphics_draw();
   } else {
      std::cout << "Molecule number " << imol_orig << " is not a valid "
		<< "model molecule" << std::endl;
   }
   return imol;
}

int new_molecule_by_sphere_selection(int imol_orig, float x, float y, float z, float r,
				     short int allow_partial_residues_flag) {

   int imol = -1;
   if (is_valid_model_molecule(imol_orig)) {
      imol = graphics_info_t::create_molecule();
      mmdb::Manager *mol_orig = graphics_info_t::molecules[imol_orig].atom_sel.mol;
      int SelectionHandle = mol_orig->NewSelection();

      mmdb::Manager *mol = NULL;
      if (allow_partial_residues_flag) {
	 mol_orig->SelectSphere(SelectionHandle, mmdb::STYPE_ATOM,
				x, y, z, r, mmdb::SKEY_OR);
	 mol = coot::util::create_mmdbmanager_from_atom_selection(mol_orig,
								  SelectionHandle);
      } else {
	 graphics_info_t g;
	 mol_orig->SelectSphere(SelectionHandle, mmdb::STYPE_RESIDUE,
				x, y, z, r, mmdb::SKEY_OR);
	 std::string alt_conf = "";

	 // convert for mmdb residue list to std::vector or residues.
	 std::vector<mmdb::Residue *> residues;
	 mmdb::PPResidue SelResidues = 0;
	 int nSelResidues;
	 mol_orig->GetSelIndex(SelectionHandle, SelResidues, nSelResidues);
	 for (int i=0; i<nSelResidues; i++)
	    residues.push_back(SelResidues[i]);

	 std::pair<mmdb::Manager *, std::vector<mmdb::Residue *> > mp =
	    g.create_mmdbmanager_from_res_vector(residues,
						 imol_orig, // for uddatom index.
						 mol_orig, alt_conf);
	 mol = mp.first;
      }

      if (mol) {
	 std::string name = "sphere selection from ";
	 name += graphics_info_t::molecules[imol_orig].name_for_display_manager();
	 atom_selection_container_t asc = make_asc(mol);
	 if (asc.n_selected_atoms > 0){
	    bool shelx_flag = 0;
	    if (graphics_info_t::molecules[imol_orig].is_from_shelx_ins())
	       shelx_flag = 1;
	    graphics_info_t g;
	    g.molecules[imol].install_model(imol, asc, g.Geom_p(), name, 1, shelx_flag);
	    g.molecules[imol].set_have_unsaved_changes_from_outside();
	 } else {
	    graphics_info_t::erase_last_molecule();
	    std::cout << "in new_molecule_by_atom_selection "
		      << "Something bad happened - No atoms selected"
		      << std::endl;
	    std::string s = "Oops! failed to create fragment.  ";
	    s += "No atoms selected\n";
	    s += "Incorrect position or radius? ";
	    s += "Radius ";
	    s += coot::util::float_to_string(r);
	    s += " at (";
	    s += coot::util::float_to_string(x);
	    s += ", ";
	    s += coot::util::float_to_string(y);
	    s += ", ";
	    s += coot::util::float_to_string(z);
	    s += ")";
	    info_dialog(s.c_str());
	    imol = -1;
	 }
      } else {
	 // mol will (currently) never be null,
	 // create_mmdbmanager_from_atom_selection() always returns a
	 // good mmdb::Manager pointer.
	 graphics_info_t::erase_last_molecule();
	 std::cout << "in new_molecule_by_atom_selection "
		   << "Something bad happened - null molecule" << std::endl;
	 std::string s = "Oops! failed to create fragment.  ";
	 s += "No atoms selected\n";
	 s += "Incorrect position or radius? ";
	 s += "Radius ";
	 s += coot::util::float_to_string(r);
	 s += " at (";
	 s += coot::util::float_to_string(x);
	 s += ", ";
	 s += coot::util::float_to_string(y);
	 s += ", ";
	 s += coot::util::float_to_string(z);
	 s += ")";
	 info_dialog(s.c_str());
	 imol = -1;
      }
      mol_orig->DeleteSelection(SelectionHandle);
      graphics_draw();
   } else {
      std::cout << "Molecule number " << imol_orig << " is not a valid "
		<< "model molecule" << std::endl;
   }
   return imol;
}

#ifdef USE_PYTHON
/*! \brief create a new molecule that consists of only the atoms
  of the specified list of residues
@return the new molecule number, -1 means an error. */
int new_molecule_by_residue_specs_py(int imol, PyObject *residue_spec_list_py) {

   int imol_new = -1;
   if (is_valid_model_molecule(imol)) {
      std::vector<coot::residue_spec_t> specs = py_to_residue_specs(residue_spec_list_py);
      if (specs.size()) {
	 graphics_info_t g;
	 mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
	 mmdb::Manager *mol_new = coot::util::create_mmdbmanager_from_residue_specs(specs, mol);
	 if (mol_new) {
	    imol_new = graphics_info_t::create_molecule();
	    atom_selection_container_t asc = make_asc(mol_new);
	    std::string label = "residues-selected-from-mol-";
	    label += coot::util::int_to_string(imol);
	    g.molecules[imol_new].install_model(imol_new, asc, g.Geom_p(), label, 1);
	    graphics_draw();
	 }
      }
   }
   return imol_new;

}
#endif /* USE_PYTHON */

#ifdef USE_GUILE
/*! \brief create a new molecule that consists of only the atoms
  of the specified list of residues
@return the new molecule number, -1 means an error. */
int new_molecule_by_residue_specs_scm(int imol, SCM residue_spec_list_scm) {

   int imol_new = -1;
   if (is_valid_model_molecule(imol)) {
      if (scm_is_true(scm_list_p(residue_spec_list_scm))) {
	 SCM len_scm = scm_length(residue_spec_list_scm);
	 int len = scm_to_int(len_scm);
	 if (len > 0) {
	    std::vector<coot::residue_spec_t> residue_specs;
	    for (int i=0; i<len; i++) {
	       SCM spec_scm = scm_list_ref(residue_spec_list_scm, scm_from_int(i));
	       coot::residue_spec_t spec = residue_spec_from_scm(spec_scm);
	       if (! spec.empty()) {
		  residue_specs.push_back(spec);
	       }
	    }
	    if (! residue_specs.empty()) {
	       mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
	       mmdb::Manager *mol_new = coot::util::create_mmdbmanager_from_residue_specs(residue_specs,mol);
	       if (mol_new) {
		  graphics_info_t g;
		  imol_new = graphics_info_t::create_molecule();
		  atom_selection_container_t asc = make_asc(mol_new);
		  std::string label = "residues-selected-from-mol-";
		  label += coot::util::int_to_string(imol);
		  g.molecules[imol_new].install_model(imol_new, asc, g.Geom_p(), label, 1);
		  graphics_draw();
	       }
	    }
	 }
      }
   }
   return imol_new;
}
#endif /* USE_GUILE */


// ---------------------------------------------------------------------
// b-factor
// ---------------------------------------------------------------------


void set_default_temperature_factor_for_new_atoms(float new_b) {

   graphics_info_t::default_new_atoms_b_factor = new_b;
}

float default_new_atoms_b_factor() {
   return graphics_info_t::default_new_atoms_b_factor;
}

void set_reset_b_factor_moved_atoms(int state) {

    graphics_info_t::reset_b_factor_moved_atoms_flag = state;
}

int get_reset_b_factor_moved_atoms_state() {

    return graphics_info_t::reset_b_factor_moved_atoms_flag;
}

#ifdef USE_GUILE
void set_temperature_factors_for_atoms_in_residue_scm(int imol, SCM residue_spec_scm, float b_factor) {

   if (is_valid_model_molecule(imol)) {
      std::pair<bool, coot::residue_spec_t> res_spec = make_residue_spec(residue_spec_scm);
      if (res_spec.first) {
	 std::vector<coot::residue_spec_t> res_spec_vec;
	 res_spec_vec.push_back(res_spec.second);
	 graphics_info_t::molecules[imol].set_b_factor_residue(res_spec.second, b_factor);
      }
   }
}
#endif // USE_GUILE

/*  ----------------------------------------------------------------------- */
/*                  SHELX stuff                                             */
/*  ----------------------------------------------------------------------- */

/* section SHELXL Functions */
// return
int read_shelx_ins_file(const std::string &filename, short int recentre_flag) {

   int istat = -1;
   graphics_info_t g;
   if (true) {
      int imol = graphics_info_t::create_molecule();

      // ugly method to recente the molecule on read
      // (save, set and reset state).
      //
      short int reset_centre_flag = g.recentre_on_read_pdb;
      g.recentre_on_read_pdb = recentre_flag;

      istat = g.molecules[imol].read_shelx_ins_file(std::string(filename));
      if (istat != 1) {
	 graphics_info_t::erase_last_molecule();
	 std::cout << "WARNING:: " << istat << " on read_shelx_ins_file "
		   << filename << std::endl;
      } else {
	 // std::cout << "DEBUG:: Molecule " << imol << " read successfully\n";
         logger.log(log_t::DEBUG, "Molecule", imol, "read successfully");
	 istat = imol; // for return status
	 if (g.go_to_atom_window) {

	    // See comments in
	    // handle_read_draw_molecule_with_recentre() about this.
	    // g.set_go_to_atom_molecule(imol); // No. 20090620

	    g.update_go_to_atom_window_on_new_mol();
	 }
	 graphics_draw();
	 std::vector<std::string> command_strings;
	 command_strings.push_back("read-shelx-ins-file");
	 command_strings.push_back(single_quote(coot::util::intelligent_debackslash(filename)));
	 add_to_history(command_strings);
      }
      g.recentre_on_read_pdb = reset_centre_flag;
   } else {
      std::cout << "ERROR:: null filename in read_shelx_ins_file" << std::endl;
   }
   return istat;
}

int write_shelx_ins_file(int imol, const char *filename) {

   int istat = 0;
   if (filename) {

      if (is_valid_model_molecule(imol)) {
	 std::pair<int, std::string> stat =
	    graphics_info_t::molecules[imol].write_shelx_ins_file(std::string(filename));
	 istat = stat.first;
	 graphics_info_t g;
	 g.add_status_bar_text(stat.second);
	 std::cout << stat.second << std::endl;
	 if (istat != 1) {
	    // wrapped_nothing_bad_dialog(stat.second);
            info_dialog(stat.second.c_str());
	 }
      } else {
	 std::cout << "WARNING:: invalid molecule (" << imol
		   << ") for write_shelx_ins_file" << std::endl;
      }
   }
   return istat;
}

void add_shelx_string_to_molecule(int imol, const char *str) {

   if (is_valid_model_molecule(imol))
      graphics_info_t::molecules[imol].add_shelx_string_to_molecule(str);
}



#ifdef USE_GUILE
SCM chain_id_for_shelxl_residue_number(int imol, int resno) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      std::pair<bool, std::string> ch =
	 graphics_info_t::molecules[imol].chain_id_for_shelxl_residue_number(resno);
      if (ch.first)
	 r = scm_from_locale_string(ch.second.c_str());
   }
   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *chain_id_for_shelxl_residue_number_py(int imol, int resno) {

   PyObject *r;
   r = Py_False;
   if (is_valid_model_molecule(imol)) {
      std::pair<bool, std::string> ch =
	 graphics_info_t::molecules[imol].chain_id_for_shelxl_residue_number(resno);
      if (ch.first)
	 r = myPyString_FromString(ch.second.c_str());
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif // USE_PYTHON








/*  ----------------------------------------------------------------------- */
/*                  SMILES                                                  */
/*  ----------------------------------------------------------------------- */
// 20230605-PE Old hideous scripting thing. Delete.
void do_smiles_gui() {

#if defined USE_GUILE
#if defined USE_GUILE_GTK
   safe_scheme_command("(smiles-gui)");
#else
#ifdef USE_PYGTK
   safe_python_command("smiles_gui()");
#endif // USE_PYGTK
#endif // USE_GUILE_GTK
#else
#ifdef USE_PYGTK
   safe_python_command("smiles_gui()");
#endif // USE_PYGTK
#endif // USE_GUILE

}

void set_residue_density_fit_scale_factor(float f) {

   graphics_info_t::residue_density_fit_scale_factor = f;
}

float residue_density_fit_scale_factor() {
   return graphics_info_t::residue_density_fit_scale_factor;
}


// dictionary
int handle_cif_dictionary(const std::string &filename) {

   short int new_molecule_flag = 0; // no
   return handle_cif_dictionary_for_molecule(filename, coot::protein_geometry::IMOL_ENC_ANY,
					     new_molecule_flag);
}

// imol_enc can be the model molecule number or
// IMOL_ENC_ANY for all
// IMOL_ENC_AUTO for auto
// IMOL_ENC_UNSET for unset - not useful.
//
// This function now handles the optional generation of a new molecule
// based on the value of new_molecule_from_dictionary_cif_checkbutton_state
//
int handle_cif_dictionary_for_molecule(const std::string &filename, int imol_enc,
				       short int new_molecule_from_dictionary_cif_checkbutton_state) {

   graphics_info_t g;
   short int show_dialog_flag = 0;
   if (graphics_info_t::use_graphics_interface_flag)
      show_dialog_flag = 1;
   // add_cif_dictionary() returns the comp_id index
   coot::read_refmac_mon_lib_info_t rmit =
      g.add_cif_dictionary(coot::util::intelligent_debackslash(filename),
			   imol_enc, show_dialog_flag);


   /* if we choose a specific molecule for this dictionary, then it doesn't make
       sense to create a new molecule to which this dictionary does not refer!
   */
   bool do_new_molecule = true;

   const char *s1 = "Molecule Select type Auto disables Generate a Molecule for non-auto-load residue type";
   const char *s2 = "Molecule Select type for a specific molecule disables Generate a Molecule";

   if (rmit.success) {
      if (imol_enc >= 0 || imol_enc == coot::protein_geometry::IMOL_ENC_AUTO) {
	 if (imol_enc == coot::protein_geometry::IMOL_ENC_AUTO) {
	    if (g.Geom_p()->is_non_auto_load_ligand(rmit.comp_id)) {
	       std::cout << "INFO:: " << s1 << std::endl;
	       add_status_bar_text(s1);
	       do_new_molecule = false;
	    }
	 } else {
	    std::cout << "INFO:: " << s2 << std::endl;
	    add_status_bar_text(s2);
	    do_new_molecule = false;
	 }

	 if (false)
	    std::cout << "handle_cif_dictionary_for_molecule here with rmit "
		      << rmit.success << " " << rmit.comp_id << " "
		      << do_new_molecule << std::endl;
      }

      if (do_new_molecule)
	 if (new_molecule_from_dictionary_cif_checkbutton_state)
	    get_monomer_for_molecule_by_index(rmit.monomer_idx, imol_enc);
   }

   graphics_draw();
   return rmit.monomer_idx;
}

int read_cif_dictionary(const std::string &filename) {

   return handle_cif_dictionary(filename);

}

/*! \brief some programs produce PDB files with ATOMs where there
  should be HETATMs.  This is a function to assign HETATMs as per the
  PDB definition. */
int assign_hetatms(int imol) {

   int r = 0;
   if (is_valid_model_molecule(imol)) {
      r = graphics_info_t::molecules[imol].assign_hetatms();
   }
   return r;
}

/*! \brief if this is not a standard group, then turn the atoms to HETATMs.

Return 1 on atoms changes, 0 on not. Return -1 if residue not found.
*/
int hetify_residue(int imol, const char * chain_id, int resno, const char *ins_code) {

   int r = -1;
   if (is_valid_model_molecule(imol)) {
      r = graphics_info_t::molecules[imol].hetify_residue_atoms(chain_id, resno, ins_code);
      graphics_draw();
   }
   return r;
}

/*! \brief residue has HETATMs?
return 1 if all atoms of the specified residue are HETATMs, else,
return 0.  If residue not found, return -1. */
int residue_has_hetatms(int imol, const char * chain_id, int resno, const char *ins_code) {

   int r = -1;
   if (is_valid_model_molecule(imol)) {
      r = graphics_info_t::molecules[imol].residue_has_hetatms(chain_id, resno, ins_code);
   }
   return r;

}


//! Add hydrogens to imol from the given pdb file
void add_hydrogens_from_file(int imol, std::string pdb_with_Hs_file_name) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      g.molecules[imol].add_hydrogens_from_file(pdb_with_Hs_file_name);
      graphics_draw();
   }
}

//! \brief add hydrogen atoms to the specified residue
void add_hydrogen_atoms_to_residue(int imol, std::string chain_id, int res_no, std::string ins_code) {

   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t rs(chain_id, res_no, ins_code);
      graphics_info_t::molecules[imol].add_hydrogen_atoms_to_residue(rs);
      graphics_draw();
   }
}

#ifdef USE_PYTHON
//! \brief add hydrogen atoms to the specified residue
void add_hydrogen_atoms_to_residue_py(int imol, PyObject *residue_spec_py) {

   if (is_valid_model_molecule(imol)) {
      coot::residue_spec_t rs = residue_spec_from_py(residue_spec_py);
      graphics_info_t::molecules[imol].add_hydrogen_atoms_to_residue(rs);
      graphics_draw();
   }
}
#endif




/*  ----------------------------------------------------------------------- */
/*                  CNS data stuff                                          */
/*  ----------------------------------------------------------------------- */
int handle_cns_data_file_with_cell(const char *filename, int imol, float a, float b, float c, float alpha, float beta, float gamma, const char *spg_info) {

   clipper::Spacegroup sg;
   clipper::Cell cell;
   clipper::Cell_descr cell_d(a, b, c,
			      clipper::Util::d2rad(alpha),
			      clipper::Util::d2rad(beta),
			      clipper::Util::d2rad(alpha));
   clipper::Spgr_descr sg_d(spg_info);
   cell.init(cell_d);
   sg.init(sg_d);
   int imol_new = graphics_info_t::create_molecule();
   int istat = graphics_info_t::molecules[imol_new].make_map_from_cns_data(sg, cell, filename);
   if (istat != -1) {
      graphics_draw();
   }
   return istat;
}


int handle_cns_data_file(const char *filename, int imol_coords) {

   int istat = -1; // returned int
   // first, does the file exist?
   struct stat s;
   int status = stat(filename, &s);
   // stat check the link targets not the link itself, lstat stats the
   // link itself.
   //
   if (status != 0 || !S_ISREG (s.st_mode)) {
      std::cout << "Error reading " << filename << std::endl;
      return -1; // which is status in an error
   } else {
      if (S_ISDIR(s.st_mode)) {
	 std::cout << filename << " is a directory." << std::endl;
      } else {
	 if (is_valid_model_molecule(imol_coords)) {
	    int imol = graphics_info_t::create_molecule();
	    std::pair<bool, clipper::Spacegroup> sg =
	       graphics_info_t::molecules[imol_coords].space_group();
	    std::pair<bool,clipper::Cell> cell =  graphics_info_t::molecules[imol_coords].cell();
	    if (sg.first && cell.first) {
	       istat = graphics_info_t::molecules[imol].make_map_from_cns_data(sg.second,
									       cell.second,
									       filename);
	       if (istat != -1) {
		  graphics_draw();
	       } else {
		  graphics_info_t::erase_last_molecule();
	       }
	    } else {
	       graphics_info_t::erase_last_molecule();
	    }
	 }
      }
   }
   return istat;
}


void
set_moving_atoms(double phi, double psi) {

   graphics_info_t g;
   g.set_edit_phi_psi_to(phi, psi);
}

void
accept_phi_psi_moving_atoms() {

   graphics_info_t g;
   g.accept_moving_atoms();
   clear_moving_atoms_object();

}

void
setup_edit_phi_psi(short int state) {

   graphics_info_t g;
   g.in_edit_phi_psi_define = state;
   if (state) {
      g.pick_cursor_maybe();
      g.pick_pending_flag = 1;

      std::cout << "click on an atom in the residue for phi/psi editting"
		<< std::endl;
   } else {
      g.normal_cursor();
   }
}




// not the write place for this function.  c-interface-map.cc would be better.
int laplacian (int imol) {

   int iret = -1;
   if (is_valid_map_molecule(imol)) {
      clipper::Xmap<float> xmap = coot::util::laplacian_transform(graphics_info_t::molecules[imol].xmap);
      int new_molecule_number = graphics_info_t::create_molecule();
      bool is_em_flag = graphics_info_t::molecules[imol].is_EM_map();
      std::string label = "Laplacian of ";
      label += graphics_info_t::molecules[imol].name_;
      graphics_info_t::molecules[new_molecule_number].install_new_map(xmap, label, is_em_flag);
      iret = new_molecule_number;
   }
   return iret;
}





void
show_partial_charge_info(int imol, const char *chain_id, int resno, const char *ins_code) {

   if (is_valid_model_molecule(imol)) {
      mmdb::Residue *residue =
	 graphics_info_t::molecules[imol].get_residue(chain_id, resno, ins_code);
      if (residue) {
	 std::string resname = residue->GetResName();
	 int read_number = graphics_info_t::cif_dictionary_read_number;
	 graphics_info_t g;
	 if (g.Geom_p()->have_dictionary_for_residue_type(resname, imol, read_number)) {

	 }
	 graphics_info_t::cif_dictionary_read_number++;
      }
   }
}

// -----------------------------------------
//
// -----------------------------------------
/*! \brief add an alternative conformer to a residue.  Add it in
  conformation rotamer number rotamer_number.

 Return #f on fail, the altconf string on success */
#ifdef USE_GUILE
SCM add_alt_conf_scm(int imol, const char*chain_id, int res_no, const char *ins_code,
		     const char *alt_conf, int rotamer_number) {

   SCM r=SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      // returns a bool,string pair (done-correct-flag, new atoms alt conf string)
      std::pair<bool, std::string> p =
	 g.split_residue(imol, std::string(chain_id), res_no,
			 std::string(ins_code), std::string(alt_conf));
      std::cout << "debug:: split_residue() returned " << p.first << " \"" << p.second << "\"" << std::endl;
      if (p.first) {
	 r = scm_from_locale_string(p.second.c_str());
      }
   }
   return r;
}
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *add_alt_conf_py(int imol, const char*chain_id, int res_no, const char *ins_code,
			  const char *alt_conf, int rotamer_number) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      // returns a bool,string pair (done-correct-flag, new atoms alt conf string)
      std::pair<bool, std::string> p =
	 g.split_residue(imol, std::string(chain_id), res_no,
			 std::string(ins_code), std::string(alt_conf));
      if (p.first) {
	 r = myPyString_FromString(p.second.c_str());
      }
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif // USE_PYTHON


/*  ----------------------------------------------------------------------- */
/*                  Backrub                                                 */
/*  ----------------------------------------------------------------------- */

/*! \brief set the mode of rotamer search, options are (ROTAMERSEARCHAUTOMATIC),
  (ROTAMERSEARCHLOWRES) (aka. "backrub rotamers),
  (ROTAMERSEARCHHIGHRES) (with rigid body fitting) */
void set_rotamer_search_mode(int mode) {

   if ((mode == ROTAMERSEARCHAUTOMATIC) ||
       (mode == ROTAMERSEARCHLOWRES) ||
       (mode == ROTAMERSEARCHHIGHRES)) {
      graphics_info_t::rotamer_search_mode = mode;
   } else {
      std::string m = "Rotamer Mode ";
      m += coot::util::int_to_string(mode);
      m += " not found";
      add_status_bar_text(m.c_str());
      std::cout << m << std::endl;
   }
}

int rotamer_search_mode_state() {

   return graphics_info_t::rotamer_search_mode;
}


/*! \name Backrubbing function */
/*! \{ */
/* \brief do a back-rub rotamer search (with autoaccept) */
int backrub_rotamer(int imol, const char *chain_id, int res_no,
		    const char *ins_code, const char *alt_conf) {

  int status = 0;

  if (is_valid_model_molecule(imol)) {
     graphics_info_t g;
     int imol_map = g.Imol_Refinement_Map();
     if (is_valid_map_molecule(imol_map)) {

	std::pair<bool,float> brs  =
	   g.molecules[imol].backrub_rotamer(chain_id, res_no, ins_code, alt_conf,
					     *g.Geom_p());
	status = brs.first;
	graphics_draw();

     } else {
	std::cout << "   WARNING:: " << imol_map << " is not a valid map molecule"
		  << std::endl;
     }
  } else {
     std::cout << "   WARNING:: " << imol << " is not a valid model molecule"
	       << std::endl;
  }
  return status;
}
/*! \} */


int backrub_rotamer_intermediate_atoms() {

   graphics_info_t g;
   return g.backrub_rotamer_intermediate_atoms();
}


/* add a linked residue based purely on dictionary templete.
   For addition of NAG to ASNs typically.

   This doesn't work with residues with alt confs.

   return success status (0 = fail).
*/
int add_linked_residue(int imol, const char *chain_id, int resno, const char *ins_code,
		       const char *new_residue_comp_id, const char *link_type, int n_trials) {

   // Are you sure that this is the function that you want to edit?

   int status = 0;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      if (g.Geom_p()->have_dictionary_for_residue_type_no_dynamic_add(new_residue_comp_id, imol)) {
      } else {
	 g.Geom_p()->try_dynamic_add(new_residue_comp_id, g.cif_dictionary_read_number);
      }
      g.cif_dictionary_read_number++;
      coot::residue_spec_t res_spec(chain_id, resno, ins_code);
// 	 g.molecules[imol].add_linked_residue(res_spec, new_residue_comp_id,
// 					      link_type, g.Geom_p());
      float new_b = g.default_new_atoms_b_factor;
      // 20140429
      coot::residue_spec_t new_res_spec =
	 g.molecules[imol].add_linked_residue_by_atom_torsions(res_spec, new_residue_comp_id,
							       link_type, g.Geom_p(), new_b);

      if (! new_res_spec.unset_p()) {
	 if (is_valid_map_molecule(imol_refinement_map())) {
	    const clipper::Xmap<float> &xmap = g.molecules[imol_refinement_map()].xmap;
	    std::vector<coot::residue_spec_t> residue_specs;
	    residue_specs.push_back(res_spec);
	    residue_specs.push_back(new_res_spec);
	    g.molecules[imol].multi_residue_torsion_fit(residue_specs, xmap, n_trials, g.Geom_p());
	 }
      }
      graphics_draw();
   }
   return status;
}

/* add a linked residue based purely on dictionary template.
   For addition of NAG to ASNs typically.

   This doesn't work with residues with alt confs.

   return status is #f for fail and the spec of the added residue on success.
*/
// mode is either 1: add  2: add and fit  3: add, fit and refine
#ifdef USE_GUILE
SCM add_linked_residue_scm(int imol, const char *chain_id, int resno, const char *ins_code,
			   const char *new_residue_comp_id, const char *link_type, int mode) {

   int n_trials = 5000;
   SCM r = SCM_BOOL_F;
   // bool do_fit_and_refine = graphics_info_t::linked_residue_fit_and_refine_state;

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;

      if (g.Geom_p()->have_dictionary_for_residue_type_no_dynamic_add(new_residue_comp_id, imol)) {
      } else {
	 std::cout << "INFO:: dictionary does not already have " << new_residue_comp_id
		   << " dynamic add it now" << std::endl;
	 int status = g.Geom_p()->try_dynamic_add(new_residue_comp_id, g.cif_dictionary_read_number);
	 if (status == 0) { // fail
	    std::cout << "WARNING:: failed to add dictionary for " << new_residue_comp_id << std::endl;
	 }
      }
      g.cif_dictionary_read_number++;
      coot::residue_spec_t res_spec(chain_id, resno, ins_code);

      float new_b = g.default_new_atoms_b_factor;
      // 20140429
      coot::residue_spec_t new_res_spec =
	 g.molecules[imol].add_linked_residue_by_atom_torsions(res_spec, new_residue_comp_id,
							       link_type, g.Geom_p(), new_b);

      // this is a new residue and so we don't want to use extra restraints
      // from an old residue (different type and different position) that might
      // have the same spec.
      //
      graphics_info_t::molecules[imol].delete_extra_restraints_for_residue(new_res_spec);

      if (mode > 1) {
	 if (! new_res_spec.unset_p()) {
	    r = residue_spec_to_scm(new_res_spec);
	    if (is_valid_map_molecule(imol_refinement_map())) {
	       const clipper::Xmap<float> &xmap = g.molecules[imol_refinement_map()].xmap;
	       std::vector<coot::residue_spec_t> residue_specs;
	       residue_specs.push_back(res_spec);
	       residue_specs.push_back(new_res_spec);

	       int n_rounds_of_fit_and_refine = 1;

	       for (int ii=0; ii<n_rounds_of_fit_and_refine; ii++) {
		  g.molecules[imol].multi_residue_torsion_fit(residue_specs, xmap, n_trials, g.Geom_p());

		  if (mode > 2) {
		     // refine and re-torsion-fit
		     int replace_mode = graphics_info_t::refinement_immediate_replacement_flag;
		     std::string alt_conf;
		     graphics_info_t::refinement_immediate_replacement_flag = 1;
		     refine_residues_with_alt_conf(imol, residue_specs, alt_conf);
		     accept_regularizement();
		     remove_initial_position_restraints(imol, residue_specs);
		     graphics_info_t::refinement_immediate_replacement_flag = replace_mode;
		  }
	       }
	    }
	 }
      }
      graphics_draw();
   }
   return r;
}
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *add_linked_residue_py(int imol, const char *chain_id, int resno, const char *ins_code,
				const char *new_residue_comp_id, const char *link_type, int mode) {

   int n_trials = 6000;
   PyObject *r = Py_False;
   bool do_fit_and_refine = graphics_info_t::linked_residue_fit_and_refine_state;

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      if (g.Geom_p()->have_dictionary_for_residue_type_no_dynamic_add(new_residue_comp_id, imol)) {
      } else {
	 std::cout << "INFO:: dictionary does not already have " << new_residue_comp_id
		   << " dynamic add it now" << std::endl;
	 g.Geom_p()->try_dynamic_add(new_residue_comp_id, g.cif_dictionary_read_number);
      }
      g.cif_dictionary_read_number++;
      coot::residue_spec_t res_spec(chain_id, resno, ins_code);
      float new_b = g.default_new_atoms_b_factor;

      // 20140429
      coot::residue_spec_t new_res_spec =
	 g.molecules[imol].add_linked_residue_by_atom_torsions(res_spec, new_residue_comp_id,
							       link_type, g.Geom_p(), new_b);

      // this is a new residue and so we don't want to use extra restraints
      // from an old residue (different type and different position) that might
      // have the same spec.
      //
      graphics_info_t::molecules[imol].delete_extra_restraints_for_residue(new_res_spec);

      if (do_fit_and_refine) {
         if (! new_res_spec.unset_p()) {
            r = residue_spec_to_py(new_res_spec);
            if (is_valid_map_molecule(imol_refinement_map())) {
               const clipper::Xmap<float> &xmap =
                  g.molecules[imol_refinement_map()].xmap;
               std::vector<coot::residue_spec_t> residue_specs;
               residue_specs.push_back(res_spec);
               residue_specs.push_back(new_res_spec);

               // 2 rounds of fit then refine
               for (int ii=0; ii<2; ii++) {
                  g.molecules[imol].multi_residue_torsion_fit(residue_specs, xmap, n_trials, g.Geom_p());

                  // refine and re-torsion-fit
                  int rep_mode = graphics_info_t::refinement_immediate_replacement_flag;
                  std::string alt_conf;
                  graphics_info_t::refinement_immediate_replacement_flag = 1;
                  refine_residues_with_alt_conf(imol, residue_specs, alt_conf);
                  accept_regularizement();
                  remove_initial_position_restraints(imol, residue_specs);
                  graphics_info_t::refinement_immediate_replacement_flag = rep_mode;
               }
            }
         }
      }
      graphics_draw();
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif

void set_add_linked_residue_do_fit_and_refine(int state) {

   graphics_info_t::linked_residue_fit_and_refine_state = state;
}

#ifdef USE_PYTHON
// return the number of atoms added
int
add_residue_with_atoms_py(int imol, PyObject *residue_spec_py, const std::string &res_name, PyObject *list_of_atoms_py) {

   int n_added = 0;
   if (is_valid_model_molecule(imol)) {
     coot::residue_spec_t res_spec = residue_spec_from_py(residue_spec_py);
     std::vector<coot::minimol::atom> list_of_atoms;
     graphics_info_t g;
     if (PyList_Check(list_of_atoms_py)) {
        Py_ssize_t p_len = PyList_Size(list_of_atoms_py);
        for (unsigned int i=0; i<p_len; i++) {
           PyObject *atom_py = PyList_GetItem(list_of_atoms_py, i);
           if (PyList_Check(atom_py)) {
              Py_ssize_t a_len = PyList_Size(atom_py);
              if (a_len == 3 || a_len == 4) {
                 PyObject *name_list_py    = PyList_GetItem(atom_py, 0);
                 PyObject *ele_occ_list_py = PyList_GetItem(atom_py, 1);
                 PyObject *pos_list_py     = PyList_GetItem(atom_py, 2);
                 if (PyList_Check(name_list_py)) {
                    if (PyList_Check(ele_occ_list_py)) {
                       if (PyList_Check(pos_list_py)) {
                          Py_ssize_t name_list_len    = PyList_Size(name_list_py);
                          Py_ssize_t ele_occ_list_len = PyList_Size(ele_occ_list_py);
                          Py_ssize_t pos_list_len     = PyList_Size(pos_list_py);
                          if (name_list_len == 2) {
                             if (ele_occ_list_len == 4) {
                                if (pos_list_len == 3) {
                                   PyObject *name_py     = PyList_GetItem(name_list_py, 0);
                                   PyObject *alt_conf_py = PyList_GetItem(name_list_py, 1);
                                   PyObject *occ_py      = PyList_GetItem(ele_occ_list_py, 0);
                                   PyObject *b_py        = PyList_GetItem(ele_occ_list_py, 1);
                                   PyObject *ele_py      = PyList_GetItem(ele_occ_list_py, 2);
                                   PyObject *seg_id_py   = PyList_GetItem(ele_occ_list_py, 3);
                                   PyObject *pos_x_py    = PyList_GetItem(pos_list_py, 0);
                                   PyObject *pos_y_py    = PyList_GetItem(pos_list_py, 1);
                                   PyObject *pos_z_py    = PyList_GetItem(pos_list_py, 2);
                                   std::string name     = PyBytes_AS_STRING(PyUnicode_AsUTF8String(name_py));
                                   std::string alt_conf = PyBytes_AS_STRING(PyUnicode_AsUTF8String(alt_conf_py));
                                   std::string ele      = PyBytes_AS_STRING(PyUnicode_AsUTF8String(ele_py));
                                   std::string seg_id   = PyBytes_AS_STRING(PyUnicode_AsUTF8String(seg_id_py));
                                   float x = PyFloat_AsDouble(pos_x_py);
                                   float y = PyFloat_AsDouble(pos_y_py);
                                   float z = PyFloat_AsDouble(pos_z_py);
                                   float o = PyFloat_AsDouble(occ_py);
                                   float b = PyFloat_AsDouble(b_py);
                                   clipper::Coord_orth pos(x,y,z);
                                   coot::minimol::atom at(name, ele, pos, alt_conf, o, b);
                                   list_of_atoms.push_back(at);
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
     std::cout << "extracted " << list_of_atoms.size() << " atoms from Python expression" << std::endl;
     n_added = g.molecules[imol].add_residue_with_atoms(res_spec, res_name, list_of_atoms);
   }
   return n_added;
}
#endif // USE_PYTHON


/*! \brief add or remove auto H-bond restraints */
void set_auto_h_bond_restraints(int state) {

   graphics_info_t g;
   g.make_auto_h_bond_restraints_flag = state;

}

void set_refine_hydrogen_bonds(int state) {
   set_auto_h_bond_restraints(state);
}

void set_refine_use_noughties_physics(short int state) {
   graphics_info_t::noughties_physics = state;
}

int get_refine_use_noughties_physics_state() {
   return graphics_info_t::noughties_physics;
}




void get_mol_edit_lock(std::atomic<bool> &mol_edit_lock) {
   // std::cout << "debug:: test_function_scm() trying to get the lock with mol_edit_lock " << mol_edit_lock << std::endl;
   bool unlocked = false;
   while (! mol_edit_lock.compare_exchange_weak(unlocked, true)) {
      // std::cout << "test_function_scm() failed to get the mol_edit_lock" << std::endl;
      std::this_thread::sleep_for(std::chrono::microseconds(100));
      unlocked = false;
   }
   // std::cout << "debug:: test_function_scm() got the lock" << std::endl;
}

void release_mol_edit_lock(std::atomic<bool> &mol_edit_lock) {
   mol_edit_lock = false;
   // std::cout << "debug:: test_function_scm() released the lock" << std::endl;
};

#include "ligand/libres-tracer.hh"

void res_tracer(int imol_map, const std::string &pir_file_name) {

   if (! is_valid_map_molecule(imol_map)) {
      std::cout << "not a valid map: " << imol_map << std::endl;
      return;
   }

   // std::string hklin_file_name = "coot-download/1gwd_map.mtz";
   // std::string f_col_label   = "FWT";
   // std::string phi_col_label = "PHWT";
   // std::cout << "Read mtz file " << hklin_file_name << " " << f_col_label << " " << phi_col_label << std::endl;
   // bool use_weights = false;
   // bool is_diff_map = false;

   // std::string pir_file_name = "1gwd.pir";

   coot::fasta_multi fam;
   fam.read(pir_file_name);
   double variation = 0.4; // speed
   unsigned int n_top_spin_pairs = 1000; // Use for tracing at most this many spin score pairs (which have been sorted).
   // This and variation affect the run-time (and results?)
   // n_top_spin_pairs = 1000; // was 1000

   unsigned int n_top_fragments = 2000; // was 4000 // The top 1000 fragments at least are all the same trace for no-side-chain lyso test
   float flood_atom_mask_radius = 1.0; // was 0.6 for emdb
   unsigned int n_phi_psi_trials = 100000; // was 5000
   float weight = 20.0f; // calculate this (using rmsd)
   bool with_ncs = false;
   float rmsd_cuffoff = 2.3;

   mmdb::Manager *working_mol = new mmdb::Manager;

   int imol_new = graphics_info_t::create_molecule();
   atom_selection_container_t asc = make_asc(working_mol);
   std::string label = "Building Molecule";
   const std::vector<coot::ghost_molecule_display_t> ghosts;
   bool shelx_flag = false;
   graphics_info_t g;
   g.molecules[imol_new].install_model_with_ghosts(imol_new, asc, g.Geom_p(), label, 1, ghosts,
                                                   shelx_flag, false, false);
   update_go_to_atom_window_on_new_mol();

   const clipper::Xmap<float> &xmap = g.molecules[imol_map].xmap;
   // coot::util::map_fill_from_mtz(&xmap, hklin_file_name, f_col_label, phi_col_label, "", use_weights, is_diff_map);
   float xmap_rmsd = g.molecules[imol_map].map_sigma();

   if (true) { // 20221216-PE what's going wrong with xmap?
      clipper::Cell c = xmap.cell();
      std::cout << "debug:: in res_tracer() xmap cell " << c.format() << std::endl;
   }

   int imol_new_map = g.create_molecule();
   label = "Map";
   bool is_em_map_flag = false;
   g.molecules[imol_new_map].install_new_map(xmap, label, is_em_map_flag);
   g.graphics_draw();

   watch_res_tracer_data_t *watch_data_p = new watch_res_tracer_data_t(working_mol, imol_new);
   std::cout << "post-constructor with mol_edit_lock: " << watch_data_p->mol_edit_lock << std::endl;

   // pass geom to this too.
   std::thread t(res_tracer_proc, xmap, xmap_rmsd, fam, variation, n_top_spin_pairs, n_top_fragments, rmsd_cuffoff, flood_atom_mask_radius,
                 weight, n_phi_psi_trials, with_ncs, watch_data_p);

   auto watching_timeout_func = [] (gpointer data) {
      watch_res_tracer_data_t *watch_data_p = static_cast<watch_res_tracer_data_t *>(data);
      if (false)
         std::cout << "debug:: watching_timeout_func runs... finished: " << watch_data_p->finished
                   << " lock: " << watch_data_p->mol_edit_lock
                   << " update_flag: " << watch_data_p->update_flag << std::endl;
      if (watch_data_p->update_flag) {
         watch_data_p->update_flag = false;
         graphics_info_t g;
         get_mol_edit_lock(watch_data_p->mol_edit_lock);
         atom_selection_container_t asc_new = make_asc(watch_data_p->working_mol);
         g.molecules[watch_data_p->imol_new].atom_sel = asc_new;
         g.molecules[watch_data_p->imol_new].make_bonds_type_checked();
         release_mol_edit_lock(watch_data_p->mol_edit_lock);
         if (watch_data_p->update_count == 1) {
            auto rc = g.molecules[watch_data_p->imol_new].centre_of_molecule();
            g.setRotationCentreSimple(rc);
            update_maps();
         }
         g.graphics_draw();
      }
      if (watch_data_p->finished) {
         std::cout << "Final update of working_mol..." << std::endl;
         get_mol_edit_lock(watch_data_p->mol_edit_lock);
         atom_selection_container_t asc_new = make_asc(watch_data_p->working_mol);
         graphics_info_t g;
         g.molecules[watch_data_p->imol_new].atom_sel = asc_new;
         g.molecules[watch_data_p->imol_new].make_bonds_type_checked();
         release_mol_edit_lock(watch_data_p->mol_edit_lock);
         g.graphics_draw();
      }
      int return_status = TRUE;
      if (watch_data_p->finished)
         return_status = FALSE; // don't continue
      return return_status;
   };

   g_timeout_add(500, watching_timeout_func, watch_data_p);

   t.detach();
}


// not sure where this function should live...

void to_generic_object_add_mesh(int object_number, PyObject *mesh_py) {
   // a mesh has 2 elements
   // 0: list of vncs
   // 1: list of triangles
   // easily translated to a Mesh
   if (PyList_Check(mesh_py)) {
      long l = PyObject_Length(mesh_py);
      if (l == 2) {
         PyObject *vnc_list = PyList_GetItem(mesh_py, 0);
         PyObject *triangles_list = PyList_GetItem(mesh_py, 1);
         long lv = PyObject_Length(vnc_list);
         long lt = PyObject_Length(triangles_list);

         std::vector<s_generic_vertex> vertices;
         std::vector<g_triangle> triangles;

         for (long i=0; i<lv; i++) {
            PyObject *vnc = PyList_GetItem(vnc_list, i);
            if (PyList_Check(vnc)) {
               long l_vnc = PyObject_Length(vnc);
               if (l_vnc == 3) {
                  PyObject *vertex = PyList_GetItem(vnc, 0);
                  PyObject *normal = PyList_GetItem(vnc, 1);
                  PyObject *colour = PyList_GetItem(vnc, 2);
                  long l_vert = PyObject_Length(vertex);
                  long l_norm = PyObject_Length(normal);
                  long l_col  = PyObject_Length(colour);
                  if (l_vert == 3) {
                     if (l_norm == 3) {
                        if (l_col == 4) {
                           double x = PyFloat_AsDouble(PyList_GetItem(vertex, 0));
                           double y = PyFloat_AsDouble(PyList_GetItem(vertex, 1));
                           double z = PyFloat_AsDouble(PyList_GetItem(vertex, 2));

                           double n1 = PyFloat_AsDouble(PyList_GetItem(normal, 0));
                           double n2 = PyFloat_AsDouble(PyList_GetItem(normal, 1));
                           double n3 = PyFloat_AsDouble(PyList_GetItem(normal, 2));

                           double c0 = PyFloat_AsDouble(PyList_GetItem(colour, 0));
                           double c1 = PyFloat_AsDouble(PyList_GetItem(colour, 1));
                           double c2 = PyFloat_AsDouble(PyList_GetItem(colour, 2));
                           double c3 = PyFloat_AsDouble(PyList_GetItem(colour, 3));

                           glm::vec3 pos(x,y,z);
                           glm::vec3 n(n1, n2, n3);
                           glm::vec4 c(c0, c1, c2, c3);
                           s_generic_vertex v(pos, n, c);
                           vertices.push_back(v);
                        }
                     }
                  }
               }
            }
         }
         for (long i=0; i<lt; i++) {
            PyObject *tri = PyList_GetItem(triangles_list, i);
            if (PyList_Check(tri)) {
               long l_tri = PyObject_Length(tri);
               if (l_tri == 3) {
                  int t0 = PyLong_AsLong(PyList_GetItem(tri, 0));
                  int t1 = PyLong_AsLong(PyList_GetItem(tri, 1));
                  int t2 = PyLong_AsLong(PyList_GetItem(tri, 2));
                  int vertices_size = vertices.size();
                  if (t0 < vertices_size) {
                     if (t1 < vertices_size) {
                        if (t2 < vertices_size) {
                           g_triangle t(t0, t1, t2);
                           triangles.push_back(t);
                        }
                     }
                  }
               }
            }
         }

         std::cout << "Debug:: to_generic_object_add_mesh() found "
                   << vertices.size() << " vertices and " << triangles.size() << " triangles\n";
         if (! vertices.empty()) {
            if (! triangles.empty()) {
               Mesh m(vertices, triangles);
               m.set_material_specularity(1,64);
               m.setup_buffers();
               meshed_generic_display_object o(m);
               graphics_info_t g;
               g.attach_buffers();
               g.generic_display_objects.push_back(o);
            }
         }
      }
   }
}


// move this function
void to_generic_object_attach_translation_gizmo(int object_number) {

   if (object_number >= 0) {
      graphics_info_t g;
      int ss = g.generic_display_objects.size(); // type conversion
      if (object_number < ss) {
         g.translation_gizmo.attached_to_generic_display_object_number = object_number;
         // now move the gizmo to the middle of the object
         const meshed_generic_display_object &gdo = g.generic_display_objects[object_number];
         std::optional<glm::vec3> centre = gdo.mesh.get_centre_of_mesh();
         std::optional<float> r = gdo.mesh.get_radius_of_gyration();
         if (centre) {
            if (r) {
               std::cout << "debug:: got radius of gyration r " << r.value() << std::endl;
               const glm::vec3 p(centre.value());
               coot::Cartesian pc(p.x, p.y, p.z);
               g.translation_gizmo.set_scale_absolute(r.value());
               g.translation_gizmo.set_position(pc);
               // this has been done before now
               // g.attach_buffers();
               // g.setup_draw_for_translation_gizmo();
            }
         }
         // should we draw it?
         g.translation_gizmo.attached_to_molecule_number = translation_gizmo_t::UNATTACHED;
         if (g.generic_display_objects[object_number].mesh.get_draw_this_mesh()) {
            g.translation_gizmo_mesh.set_draw_mesh_state(true);
         }
      }
      g.graphics_draw();
   }
}



// move this function
void generic_object_mesh_calculate_normals(int object_number) {

   graphics_info_t g;
   unsigned int object_number_u(object_number);
   if (object_number >= 0) {
      if (object_number_u < g.generic_display_objects.size()) {
         g.generic_display_objects[object_number].mesh.calculate_normals();
      }
   }

}
