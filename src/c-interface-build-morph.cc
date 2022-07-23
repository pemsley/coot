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
#define S_ISDIR(m) (((m) & S_IFMT) == S_IFDIR)
#define S_ISREG(m) (((m) & S_IFMT) == S_IFREG)
#include <windows.h>
#endif


#include <mmdb2/mmdb_manager.h>
#include "coords/mmdb-extras.h"
#include "coords/mmdb.h"

#ifdef EMSCRIPTEN_THING
// just delete it?
#include "globjects.h" //includes gtk/gtk.h
#endif

#include "coords/mmdb-crystal.h"

#include "coords/Cartesian.h"
#include "coords/Bond_lines.h"

#include "graphics-info.h"

#include "coot-utils/coot-coord-utils.hh"
#include "utils/coot-fasta.hh"

#include "skeleton/BuildCas.h"
#include "ligand/helix-placement.hh"
#include "ligand/fast-ss-search.hh"

#include "trackball.h" // adding exportable rotate interface

#include "utils/coot-utils.hh"  // for is_member_p
#include "coot-utils/coot-map-heavy.hh"  // for fffear

#include "guile-fixups.h"


#include "c-interface.h"
#ifdef EMSCRIPTEN_THING
#include "c-interface-gtk-widgets.h"
#endif
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


int morph_fit_all(int imol, float transformation_averaging_radius) {

   int success = 0;
   graphics_info_t g;
   int imol_ref_map = g.Imol_Refinement_Map();
   if (is_valid_map_molecule(imol_ref_map)) {
      if (is_valid_model_molecule(imol)) {
	 success = g.molecules[imol].morph_fit_all(g.molecules[imol_ref_map].xmap,
						   transformation_averaging_radius);
	 graphics_draw();
      }
   }
   return success;
}

int morph_fit_chain(int imol, std::string chain_id, float transformation_averaging_radius) {

   int success = 0;
   graphics_info_t g;
   int imol_ref_map = g.Imol_Refinement_Map();
   if (is_valid_map_molecule(imol_ref_map)) {
      if (is_valid_model_molecule(imol)) {
	 success = g.molecules[imol].morph_fit_chain(chain_id,
						     g.molecules[imol_ref_map].xmap,
						     transformation_averaging_radius);
	 graphics_draw();
      }
   }
   return success;
}


#ifdef USE_GUILE
int morph_fit_residues_scm(int imol, SCM residue_specs_scm, float transformation_averaging_radius) {

   std::vector<coot::residue_spec_t> residue_specs = scm_to_residue_specs(residue_specs_scm);
   return morph_fit_residues(imol, residue_specs, transformation_averaging_radius);
}
#endif // USE_GUILE

#ifdef USE_PYTHON
int morph_fit_residues_py( int imol, PyObject *residue_specs_py, float transformation_averaging_radius) {

   std::vector<coot::residue_spec_t> residue_specs = py_to_residue_specs(residue_specs_py);
   return morph_fit_residues(imol, residue_specs, transformation_averaging_radius);
}
#endif // USE_PYTHON

int morph_fit_residues(int imol, const std::vector<coot::residue_spec_t> &residue_specs,
		       float transformation_averaging_radius) {

   graphics_info_t g;
   int success = 0;
   int imol_ref_map = g.Imol_Refinement_Map();
   if (is_valid_map_molecule(imol_ref_map)) {
      if (is_valid_model_molecule(imol)) {
	 const clipper::Xmap<float> &xmap = g.molecules[imol_ref_map].xmap;
	 success = g.molecules[imol].morph_fit_residues(residue_specs, xmap,
							transformation_averaging_radius);
	 graphics_draw();
      } else {
	 std::cout << "WARNING:: no valid map. Stopping now" << std::endl;
      }
   } else {
      std::cout << "WARNING:: " << imol << " is not a valid map molecule " << std::endl;
   }
   return success;
}

int morph_fit_by_secondary_structure_elements(int imol, const std::string &chain_id) {

   graphics_info_t g;
   int success = 0;
   int imol_ref_map = g.Imol_Refinement_Map();
   if (is_valid_map_molecule(imol_ref_map)) {
      if (is_valid_model_molecule(imol)) {
	 const clipper::Xmap<float> &xmap = g.molecules[imol_ref_map].xmap;
	 float map_rmsd                   = g.molecules[imol_ref_map].map_sigma();
	 g.molecules[imol].add_secondary_structure_header_records();
	 success = g.molecules[imol].morph_fit_by_secondary_structure_elements(chain_id, xmap, map_rmsd);
	 graphics_draw();
      } else {
	 std::cout << "WARNING:: no valid map. Stopping now" << std::endl;
      }
   } else {
      std::cout << "WARNING:: " << imol << " is not a valid map molecule " << std::endl;
   }
   return success;
}


// --------------------------------------------------------------------------
//                   jiggle fit
// --------------------------------------------------------------------------
//
float fit_to_map_by_random_jiggle(int imol, const char *chain_id, int resno, const char *ins_code,
				  int n_trials,
				  float jiggle_scale_factor) {
   float val = -999.0;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      int imol_map = g.Imol_Refinement_Map();
      if (is_valid_map_molecule(imol_map)) {
	 coot::residue_spec_t rs(chain_id, resno, ins_code);
	 float map_sigma = g.molecules[imol_map].map_sigma();
	 val = g.molecules[imol].fit_to_map_by_random_jiggle(rs,
							     g.molecules[imol_map].xmap,
							     map_sigma,
							     n_trials, jiggle_scale_factor);
	 graphics_draw();
      } else {
	 std::cout << "WARNING:: Refinement map not set" << std::endl;
	 add_status_bar_text("Refinement map not set.");
      }
   }
   return val;
}

float fit_molecule_to_map_by_random_jiggle(int imol, int n_trials, float jiggle_scale_factor) {

   float val = -999.0;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      int imol_map = g.Imol_Refinement_Map();
      if (! is_valid_map_molecule(imol_map)) {
         info_dialog("WARNING:: Refinement map is not set");
      } else {
         float map_sigma = g.molecules[imol_map].map_sigma();

         mmdb::Atom **atom_selection = 0; // g.molecules[imol].atom_sel.atom_selection;
         int n_atoms = 0; // g.molecules[imol].atom_sel.n_selected_atoms;

         // select only main-chain or nucleotide atoms
         mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
         int SelHnd = mol->NewSelection(); // d
         mol->SelectAtoms(SelHnd, 0, "*",
                          mmdb::ANY_RES, "*",
                          mmdb::ANY_RES, "*", "*",
                          "CA,C,N,O,CB,P,C1',N1,C2,N3,C4,N4,O2,C5,C6,O4,N9,C8,N7,N6","*","*",mmdb::SKEY_NEW);
         mol->GetSelIndex(SelHnd, atom_selection, n_atoms);

         // fill the chains - we want to apply the tranformation to the chains, not the atom selection
         std::vector<mmdb::Chain *> chains;
         mmdb::Model *model_p = mol->GetModel(1);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++)
               chains.push_back(model_p->GetChain(ichain));
         }

         bool use_biased_density_scoring = false; // not for all-molecule
         val = g.molecules[imol].fit_to_map_by_random_jiggle(atom_selection, n_atoms,
                                                             g.molecules[imol_map].xmap,
                                                             map_sigma,
                                                             n_trials, jiggle_scale_factor,
                                                             use_biased_density_scoring, chains);
         mol->DeleteSelection(SelHnd);
         graphics_draw();
      }
   }
   return val;
}

float fit_chain_to_map_by_random_jiggle(int imol, const char *chain_id, int n_trials, float jiggle_scale_factor) {

   float val = -999.0;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      int imol_map = g.Imol_Refinement_Map();
      mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
      if (! is_valid_map_molecule(imol_map)) {
         info_dialog("WARNING:: Refinement map is not set");
      } else {
         float map_sigma = g.molecules[imol_map].map_sigma();
         const clipper::Xmap<float> &xmap = g.molecules[imol_map].xmap;
         g.molecules[imol].fit_chain_to_map_by_random_jiggle(chain_id, xmap, map_sigma, n_trials, jiggle_scale_factor);
      }
   } else {
      add_status_bar_text("Jiggle Fit: No molecule");
   }

   graphics_draw();

   return val;
}


float fit_molecule_to_map_by_random_jiggle_and_blur(int imol, int n_trials, float jiggle_scale_factor, float map_blur_factor) {
   float r = -100;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      int imol_map = g.Imol_Refinement_Map();
      mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
      if (! is_valid_map_molecule(imol_map)) {
         info_dialog("WARNING:: Refinement map is not set");
      } else {
         // happy path
         const clipper::Xmap<float> &xmap = g.molecules[imol_map].xmap;
         clipper::Xmap<float> xmap_blur = coot::util::sharpen_blur_map(xmap, map_blur_factor);
         float map_sigma = g.molecules[imol_map].map_sigma(); // is that the right map?

         mmdb::Atom **atom_selection = 0;
         int n_atoms = 0;

         // select only main-chain or nucleotide atoms
         int SelHnd = mol->NewSelection(); // d
         mol->SelectAtoms(SelHnd, 0, "*",
                          mmdb::ANY_RES, "*",
                          mmdb::ANY_RES, "*", "*",
                          "CA,C,N,O,CB,P,C1',N1,C2,N3,C4,N4,O2,C5,C6,O4,N9,C8,N7,N6","*","*",mmdb::SKEY_NEW);
         mol->GetSelIndex(SelHnd, atom_selection, n_atoms);

         // fill the chains - we want to apply the tranformation to the chains, not the atom selection
         std::vector<mmdb::Chain *> chains;
         mmdb::Model *model_p = mol->GetModel(1);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++)
            chains.push_back(model_p->GetChain(ichain));
         }

         bool use_biased_density_scoring = false; // not for all-molecule
         r = g.molecules[imol].fit_to_map_by_random_jiggle(atom_selection, n_atoms,
                                                           xmap_blur, map_sigma,
                                                           n_trials, jiggle_scale_factor,
                                                           use_biased_density_scoring, chains);
         // now fit to original map, just a few trials will do, becuase the rigid body fit
         // is what we really want
         r = g.molecules[imol].fit_to_map_by_random_jiggle(atom_selection, n_atoms,
                                                           xmap, map_sigma,
                                                           12, 0.4,
                                                           use_biased_density_scoring, chains);

         mol->DeleteSelection(SelHnd);
         graphics_draw();
      }
   }
   return r;
}

/*!  \brief jiggle fit the chain to the current refinment map
 *
 * Use a map that is blurred by the give factor for fitting.
 * @return < -100 if not possible, else return the new best fit for this chain.  */
float fit_chain_to_map_by_random_jiggle_and_blur(int imol, const char *chain_id, int n_trials, float jiggle_scale_factor, float map_blur_factor) {
   float r = -100;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      int imol_map = g.Imol_Refinement_Map();
      mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
      if (! is_valid_map_molecule(imol_map)) {
         info_dialog("WARNING:: Refinement map is not set");
      } else {
         // happy path
         const clipper::Xmap<float> &xmap = g.molecules[imol_map].xmap;
         clipper::Xmap<float> xmap_blur = coot::util::sharpen_blur_map(xmap, map_blur_factor);
         float map_sigma = g.molecules[imol_map].map_sigma(); // is that the right map?
         g.molecules[imol].fit_chain_to_map_by_random_jiggle(chain_id, xmap_blur, map_sigma, n_trials, jiggle_scale_factor);
      }
   }
   graphics_draw();
   return r;
}
