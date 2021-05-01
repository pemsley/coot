/* src/c-interface-mutate.cc
 * 
 * Copyright 2007 by Paul Emsley
 * Copyright 2007 by The University of Oxford
 * Copyright 2007 by Bernhard Lohkamp
 * Copyright 2007 by The University of York
 * Copyright 2013 by Medical Research Council
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

#include "compat/coot-sysdep.h"


#include <stdlib.h>
#include <iostream>

// ---------------------------------------------------------------
// Note to self:
//      Be very careful if you reorder include files here that the
//      new order compiles on Mac and other autobulding system.
//      There have been problems with preprocessor strangeness
//      20071010
// ---------------------------------------------------------------

#include <gtk/gtk.h>

#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include "cootaneer/cootaneer-sequence.h"
#include "graphics-info.h"

// Including python needs to come after graphics-info.h, because
// something in Python.h (2.4 - chihiro) is redefining FF1 (in
// ssm_superpose.h) to be 0x00004000 (Grrr).
// BL says:: and (2.3 - dewinter), i.e. is a Mac - Python issue
// since the follwing two include python graphics-info.h is moved up
#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "cc-interface.hh"  // includes coot-coord-utils.hh
#include "utils/coot-utils.hh"
#include "guile-fixups.h"


#ifdef USE_GUILE
coot::atom_spec_t
atom_spec_from_scm_expression(SCM expr) {

   coot::atom_spec_t atom_spec;
   atom_spec.string_user_data = "Bad Spec";

   SCM len_expr = scm_length(expr);
   int len_view = scm_to_int(len_expr);
   if (len_view == 5) {

      SCM idx_scm = SCM_MAKINUM(0);
      SCM chain_id_scm = scm_list_ref(expr, idx_scm);
      std::string chain_id = scm_to_locale_string(chain_id_scm);

      idx_scm = SCM_MAKINUM(1);
      SCM resno_scm = scm_list_ref(expr, idx_scm);
      int resno = scm_to_int(resno_scm);

      idx_scm = SCM_MAKINUM(2);
      SCM ins_code_scm = scm_list_ref(expr, idx_scm);
      std::string ins_code = scm_to_locale_string(ins_code_scm);
      
      idx_scm = SCM_MAKINUM(3);
      SCM atom_name_scm = scm_list_ref(expr, idx_scm);
      std::string atom_name = scm_to_locale_string(atom_name_scm);
      
      idx_scm = SCM_MAKINUM(4);
      SCM alt_conf_scm = scm_list_ref(expr, idx_scm);
      std::string alt_conf = scm_to_locale_string(alt_conf_scm);

      //        std::cout << "decoding spec :" << chain_id << ": " << resno << " :"
      // << ins_code // << ": :" << atom_name << ": :" << alt_conf << ":" << std::endl;
       
      atom_spec = coot::atom_spec_t(chain_id, resno, ins_code, atom_name, alt_conf);
      // if all OK:
      atom_spec.string_user_data = "OK";
   }
   
   return atom_spec;
}
#endif // USE_GUILE

#ifdef USE_PYTHON
coot::atom_spec_t
atom_spec_from_python_expression(PyObject *expr) {

   coot::atom_spec_t atom_spec;
   atom_spec.string_user_data = "Bad Spec";

   if (PyList_Check(expr)) {
      int len = PyList_Size(expr);

      // atom expressions of size 6 are prefixed by (i.e. is the 0th element) the
      // model molecule number

      if (len == 5 || len == 6) {

	 int offset = 0;
	 if (len == 6) offset = 1;
      
	 PyObject *chain_id_python = PyList_GetItem(expr, 0+offset);
	 std::string chain_id = PyString_AsString(chain_id_python);

	 PyObject *resno_python = PyList_GetItem(expr, 1+offset);
	 int resno = PyInt_AsLong(resno_python);

	 PyObject *ins_code_python = PyList_GetItem(expr, 2+offset);
	 std::string ins_code = PyString_AsString(ins_code_python);
      
	 PyObject *atom_name_python = PyList_GetItem(expr, 3+offset);
	 std::string atom_name = PyString_AsString(atom_name_python);
      
	 PyObject *alt_conf_python = PyList_GetItem(expr, 4+offset);
	 std::string alt_conf = PyString_AsString(alt_conf_python);

	 if (false)
	    std::cout << "decoding spec chain-id :" << chain_id << " resno :" << resno << " :" << ins_code
		      << ": :" << atom_name << ": :" << alt_conf << ":" << std::endl;

	 atom_spec = coot::atom_spec_t(chain_id, resno, ins_code, atom_name, alt_conf);

	 // currently atom_spec.int_user_data is -1 - use it to store the molecule
	 // number if we can.
	 //
	 if (len == 6) {
	    PyObject *o = PyList_GetItem(expr, 0);
	    if (PyInt_Check(o)) {
	       long imol = PyInt_AsLong(o);
	       atom_spec.int_user_data = imol;
	    }
	 }

	    // if all OK:
	 atom_spec.string_user_data = "OK";
      }
   }
   
   return atom_spec;
}
#endif // USE_PYTHON

/*  ----------------------------------------------------------------------- */
/*                  Cootaneer                                               */
/*  ----------------------------------------------------------------------- */

#ifdef USE_GUILE
int cootaneer(int imol_map, int imol_model, SCM atom_in_fragment_atom_spec) {

   int istat = 0;
   coot::atom_spec_t atom_spec =
      atom_spec_from_scm_expression(atom_in_fragment_atom_spec);
   if (atom_spec.string_user_data == "Bad Spec") {
      std::cout << "Bad SCM expression for atom spec" << std::endl;
      return -1; 
   } else { 
      istat = cootaneer_internal(imol_map, imol_model, atom_spec);
      graphics_draw();
   }
   return istat;
}
#endif // USE_GUILE

#ifdef USE_PYTHON
int cootaneer_py(int imol_map, int imol_model, PyObject *atom_in_fragment_atom_spec) {

   int istat = 0;
   coot::atom_spec_t atom_spec =
      atom_spec_from_python_expression(atom_in_fragment_atom_spec);
   if (atom_spec.string_user_data == "Bad Spec") {
      std::cout << "Bad Python expression for atom spec" << std::endl;
      return -1; 
   } else {
      istat = cootaneer_internal(imol_map, imol_model, atom_spec);
      graphics_draw();
   }
   return istat;
}
#endif // USE_PYTHON

#include "utils/coot-utils.hh"

int cootaneer_internal(int imol_map, int imol_model, const coot::atom_spec_t &atom_spec) {
   int istat = 0;
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {

	 std::string llkdfile = coot::package_data_dir();
	 llkdfile += "/cootaneer-llk-2.40.dat";

	 // was that over-ridden?
	 const char *cp = getenv("COOT_PREFIX");
	 if (cp) {
	    llkdfile = cp;
	    llkdfile += "/share/coot/cootaneer-llk-2.40.dat";
	 }

	 if (!coot::file_exists(llkdfile)) {

	    std::cout << "Ooops! Can't find cootaneer likelihood data! - failure"
		      << std::endl;

	 } else {

	    std::string chain_id = atom_spec.chain_id;

	    mmdb::Manager *mol = graphics_info_t::molecules[imol_model].atom_sel.mol;
	    std::pair<mmdb::Manager *, std::vector<coot::residue_spec_t> > mmdb_info =
	       coot::util::get_fragment_from_atom_spec(atom_spec, mol);

	    if (!mmdb_info.first) {
	       std::cout << "Bad - no fragment from atom spec" << std::endl;
	    } else {

	       // mmdb needs a spacegroup, it seems, or crash in sequence_chain().
	       // OK.  Let's bolt one in.
	       //
	       mmdb::realtype a[6];
	       mmdb::realtype vol;
	       int orthcode;
	       mol->GetCell(a[0], a[1], a[2], a[3], a[4], a[5], vol, orthcode);
	       char *sg = mol->GetSpaceGroup();
	       mmdb_info.first->SetCell(a[0], a[1], a[2], a[3], a[4], a[5]);
	       if (sg)
		  mmdb_info.first->SetSpaceGroup(sg);

      
	       // create sequencer
	       Coot_sequence sequencer( llkdfile );
	       std::vector<std::pair<std::string, std::string> > seq =
		  graphics_info_t::molecules[imol_model].sequence_info();

	       if (seq.size() > 0) { 

		  // and apply
		  sequencer.sequence_chain(graphics_info_t::molecules[imol_map].xmap,
					   seq, *mmdb_info.first, chain_id);
		  std::string bestseq = sequencer.best_sequence();
		  std::string fullseq = sequencer.full_sequence();
		  double conf = sequencer.confidence();
		  int chnnum = sequencer.chain_number();
		  int chnoff = sequencer.chain_offset();
   
		  // write results
		  std::cout << "\nSequence: " << bestseq << "\nConfidence: " << conf << "\n";
		  if ( chnnum >= 0 ) { 
		     std::cout << "\nFrom    : " << fullseq << "\nChain id: "
			       << chnnum << "\tOffset: " << chnoff+1 << "\n";

		     if (conf > 0.9) {
		       // BL and KC says:: we do not want to change the chain_id of the sequenced 
		       // fragment here!!
		       //			std::string chain_id = seq[chnnum].first;
			std::vector<coot::residue_spec_t> mmdb_residues = mmdb_info.second;
			graphics_info_t g;
			istat = g.molecules[imol_model].apply_sequence(imol_map,
								       mmdb_info.first,
								       mmdb_residues, bestseq,
								       chain_id, chnoff+1,
								       *g.Geom_p());
		     }
		  }
	       } else {
		  std::string s = "Oops - no sequence information has been given to molecule\n";
		  s += "number ";
		  s += coot::util::int_to_string(imol_model);
		  info_dialog(s.c_str());
	       }
	       delete mmdb_info.first;
	    }
	 }
	 
      } else {
	 std::cout << "Not a valid map molecule " << imol_model << std::endl;
      }
   } else {
      std::cout << "Not a valid model molecule " << imol_model << std::endl;
   }

   return istat;
}

#ifdef USE_GUILE
SCM find_terminal_residue_type(int imol, const char *chain_id, int resno) {
   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      std::pair<bool, std::string> p = 
	 graphics_info_t::molecules[imol].find_terminal_residue_type(chain_id, resno,
								     graphics_info_t::alignment_wgap,
								     graphics_info_t::alignment_wspace);
      if (p.first) {
	 r = scm_makfrom0str(p.second.c_str());
      }
   }
   return r;
}
#endif // GUILE

#ifdef USE_PYTHON
PyObject *find_terminal_residue_type_py(int imol, const char *chain_id, int resno) {
  PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      std::pair<bool, std::string> p = 
	 graphics_info_t::molecules[imol].find_terminal_residue_type(chain_id, resno,
								     graphics_info_t::alignment_wgap,
								     graphics_info_t::alignment_wspace);
      if (p.first) {
	 r = PyString_FromString(p.second.c_str());
      }
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif // PYTHON


void align_and_mutate(int imol, const char *chain_id, const char *fasta_maybe, short int renumber_residues_flag) {

   if (is_valid_model_molecule(imol)) {
      if (chain_id) { 
	 graphics_info_t g;
         bool auto_fit = false; // allow this to be a passed variable at some stage.
	 g.mutate_chain(imol, std::string(chain_id), std::string(fasta_maybe), auto_fit, renumber_residues_flag);
	 graphics_draw();
	 g.update_go_to_atom_window_on_changed_mol(imol);
      } else {
	 std::cout << "WARNING:: bad (NULL) chain_id - no alignment" << std::endl;
      }
   } else {
      std::cout << "WARNING:: inapproproate molecule number " << imol << std::endl;
   }
}

/*! \brief set the penalty for affine gap and space when aligning */
void set_alignment_gap_and_space_penalty(float wgap, float wspace) {

   graphics_info_t::alignment_wgap = wgap;
   graphics_info_t::alignment_wspace = wspace;
}


#ifdef USE_GUILE
SCM alignment_results_scm(int imol, const char *chain_id, const char *seq) {

   SCM r = SCM_BOOL_F;
   return r;
}
#endif /* USE_GUILE */

/*! \brief align sequence to closest chain (compare across all chains in all molecules).  

Typically match_fraction is 0.95 or so.

Return 1 if we were successful, 0 if not. */
std::pair<int, std::string>
align_to_closest_chain(std::string target_seq_in, float match_fraction_crit) {

   int status = 0;
   float match_fragment_best = 1.0;
   std::string chain_id_best;
   int imol_best = -1;
   std::string target(target_seq_in);
   
   for (int imol=0; imol<graphics_n_molecules(); imol++) {
      if (is_valid_model_molecule(imol)) {
	 if (target.length() > 0) {
	    std::pair<bool, std::pair<std::string, coot::chain_mutation_info_container_t> > r = 
	       graphics_info_t::molecules[imol].try_align_on_all_chains(target, match_fraction_crit,
									graphics_info_t::alignment_wgap,
									graphics_info_t::alignment_wspace);
	    if (r.first) {
	       float sum_changes =
		  r.second.second.single_insertions.size() +
		  r.second.second.mutations.size() +
		  r.second.second.deletions.size();
	       float match_frag = sum_changes/float(target.length());
	       if (match_frag < match_fragment_best) {
		  status = 1;
		  imol_best = imol;
		  chain_id_best = r.second.first;
	       } 
	    }
	 } 
      } 
   }

   if (status) {
      assign_sequence_from_string(imol_best, chain_id_best.c_str(), target_seq_in.c_str());
      std::cout << "INFO:: sequence assigned to chain \"" << chain_id_best
		<< "\" of molecule " << imol_best << std::endl;
   } 
   return std::pair<int, std::string> (imol_best, chain_id_best);
}

#ifdef USE_PYTHON
PyObject *align_to_closest_chain_py(std::string target_seq, float match_fraction) {

   PyObject *r = Py_False;
   std::pair<int, std::string> result = align_to_closest_chain(target_seq, match_fraction);
   if (is_valid_model_molecule(result.first)) {
      r = PyList_New(2);
      PyList_SetItem(r, 0, PyInt_FromLong(result.first));
      PyList_SetItem(r, 1, PyString_FromString(result.second.c_str()));
   }

   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
} 
#endif /* USE_PYTHON */

#ifdef USE_GUILE
SCM align_to_closest_chain_scm(std::string target_seq, float match_fraction) {

   SCM r = SCM_BOOL_F;
   std::pair<int, std::string> result = align_to_closest_chain(target_seq, match_fraction);
   if (is_valid_model_molecule(result.first)) {
      SCM a = SCM_MAKINUM(result.first);
      SCM b = scm_makfrom0str(result.second.c_str());
      r = scm_list_2(a,b);
   }
   return r;
}
#endif /* USE_GUILE */




#ifdef USE_GUILE
SCM  nearest_residue_by_sequence_scm(int imol, const char* chain_id, int resno, const char *ins_code) { 

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      coot::residue_spec_t spec(chain_id, resno, ins_code);
      mmdb::Residue *residue_p = coot::nearest_residue_by_sequence(mol, spec);
      if (residue_p) {
	 r = residue_spec_to_scm(coot::residue_spec_t(residue_p));
      }
   }
   return r;
}
#endif /* USE_GUILE */



#ifdef USE_PYTHON
PyObject *alignment_results_py(int imol, const char *chain_id, const char *seq) {

   PyObject *r = Py_False;

   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
} 
#endif /* USE_PYTHON */

#ifdef USE_PYTHON
PyObject *nearest_residue_by_sequence_py(int imol, const char* chain_id, int resno, const char *ins_code) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      coot::residue_spec_t spec(chain_id, resno, ins_code);
      mmdb::Residue *residue_p = coot::nearest_residue_by_sequence(mol, spec);
      if (residue_p) {
         r = residue_spec_to_py(coot::residue_spec_t(residue_p));
      }
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif /* USE_PYTHON */



void resolve_clashing_sidechains_by_deletion(int imol) {

   if (is_valid_model_molecule(imol)) {
      coot::protein_geometry *geom_p = graphics_info_t::Geom_p();
      graphics_info_t::molecules[imol].resolve_clashing_sidechains_by_deletion(geom_p);
      graphics_draw();
   }
}

void resolve_clashing_sidechains_by_rebuilding(int imol) {

   if (is_valid_model_molecule(imol)) {
      coot::protein_geometry *geom_p = graphics_info_t::Geom_p();
      graphics_info_t g;
      int imol_refinement_map = g.Imol_Refinement_Map();
      graphics_info_t::molecules[imol].resolve_clashing_sidechains_by_rebuilding(geom_p, imol_refinement_map);
      graphics_draw();
   }
}
