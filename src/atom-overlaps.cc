/* src/atom-overlaps.cc
 * 
 * Copyright 2015 by Medical Research Council
 * Author: Paul Emsley
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
#include <Python.h>
#endif // USE_PYTHON

// #include <stdio.h>
// #include <string.h>

// #include <string>
// #include <vector>
#include <iostream>
// #include <algorithm>
// #include <map>

// #include "coot-utils/coot-coord-utils.hh"
// #include "utils/coot-utils.hh"

#include <gtk/gtk.h>
#include "graphics-info.h"
#include "c-interface.h"   // is_valid_model_molecule()
#include "cc-interface.hh" // residue_spec_from_scm()
#include "c-interface-ligands-swig.hh" // where these functions are declared.

#include "coot-utils/atom-overlaps.hh"

#include "guile-fixups.h"

#ifdef USE_GUILE
// internal bumps scoring, sphere overlap
SCM ligand_atom_overlaps_scm(int imol, SCM ligand_spec, double neighb_radius) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      coot::residue_spec_t rs = residue_spec_from_scm(ligand_spec);
      mmdb::Residue *residue_p = g.molecules[imol].get_residue(rs);
      if (residue_p) {
	 mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
	 std::vector<mmdb::Residue *> neighb_residues =
	    coot::residues_near_residue(residue_p, mol, neighb_radius);
	 
	 if (neighb_residues.size()) {
	    coot::atom_overlaps_container_t ol(residue_p, neighb_residues, mol, g.Geom_p());
	    ol.make_overlaps();
	    r = SCM_EOL;
	    for (unsigned int i=0; i<ol.overlaps.size(); i++) {
	       SCM spec_1_scm = atom_spec_to_scm(coot::atom_spec_t(ol.overlaps[i].atom_1));
	       SCM spec_2_scm = atom_spec_to_scm(coot::atom_spec_t(ol.overlaps[i].atom_2));
	       SCM sb = SCM_BOOL_F;
	       if (ol.overlaps[i].is_h_bond) sb = SCM_BOOL_T;
	       SCM l = scm_list_4(spec_1_scm,
				  spec_2_scm,
				  scm_double2num(ol.overlaps[i].overlap_volume),
				  sb);
	       r = scm_cons(l, r);
	    }
	    r = scm_reverse(r);
	 }
      }
   }; 
   return r;
} 
#endif

#ifdef USE_PYTHON
// internal bumps scoring, sphere overlap
PyObject *ligand_atom_overlaps_py(int imol, PyObject *ligand_spec, double neighb_radius) {

   PyObject *r = Py_False;

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      coot::residue_spec_t rs = residue_spec_from_py(ligand_spec);
      mmdb::Residue *r = g.molecules[imol].get_residue(rs);
      if (r) {
	 mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
	 std::vector<mmdb::Residue *> neighb_residues =
	    coot::residues_near_residue(r, mol, neighb_radius);
	 
	 if (neighb_residues.size()) {
	    coot::atom_overlaps_container_t ol(r, neighb_residues, mol, g.Geom_p());
	    ol.make_overlaps();
	 }
      } 
   };
   
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
} 
#endif
