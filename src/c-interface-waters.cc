/* src/c-interface-waters.cc
 * 
 * Copyright 2004, 2005 by The University of York
 * Author: Paul Emsley
 * Copyright 2008 by The University of Oxford
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
 

#include <stdlib.h>
#include <iostream>

#if defined _MSC_VER
#include <windows.h>
#define snprintf _snprintf
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
#include "guile-fixups.h"

// Including python needs to come after graphics-info.h, because
// something in Python.h (2.4 - chihiro) is redefining FF1 (in
// ssm_superpose.h) to be 0x00004000 (Grrr).
//
#ifdef USE_PYTHON
#include "Python.h"
#endif // USE_PYTHON

#include "c-interface.h"
#include "cc-interface.hh"

#include "wligand.hh"

char *get_text_for_find_waters_sigma_cut_off() {

   graphics_info_t g;
   char *text = (char *) malloc (100);
   snprintf(text, 99, "%5.1f", g.find_waters_sigma_cut_off);
   return text;

}

void set_value_for_find_waters_sigma_cut_off(float f) {
   graphics_info_t::find_waters_sigma_cut_off = f;
}

/* default 0.07, I think. */
void set_water_check_spherical_variance_limit(float f) { 
   graphics_info_t g;
   g.ligand_water_variance_limit = f;
}

void set_ligand_water_to_protein_distance_limits(float f1, float f2) {
   graphics_info_t g;
   g.ligand_water_to_protein_distance_lim_min = f1;
   g.ligand_water_to_protein_distance_lim_max = f2;
}

void set_ligand_water_n_cycles(int i) {

   graphics_info_t::ligand_water_n_cycles = i;
}

/* 0 off, 1 on */
void set_ligand_verbose_reporting(int i) {
   graphics_info_t::ligand_verbose_reporting_flag = i;
}


//    short int new_waters_mol_flag = 1; // 1 mean a new molecule,
// 				         // 0 means the masking molecule.
// 
void
execute_find_waters_real(int imol_for_map,
			 int imol_for_protein,
			 short int new_waters_mol_flag, 
			 float sigma_cut_off) {

   find_waters(imol_for_map, imol_for_protein, new_waters_mol_flag, sigma_cut_off, 1);
}



/* Let's give access to the sigma level (default 4) */
float check_waters_by_difference_map_sigma_level_state() { 
   return graphics_info_t::check_waters_by_difference_map_sigma_level;
} 

void set_check_waters_by_difference_map_sigma_level(float f) { 
   graphics_info_t::check_waters_by_difference_map_sigma_level = f;
}


void renumber_waters(int imol) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].renumber_waters();
      graphics_draw();
      if (graphics_info_t::go_to_atom_window) {
	 update_go_to_atom_window_on_changed_mol(imol);
      }
   }
} 

/* Put the blob under the cursor to the screen centre.  Check only
   positive blobs.  Useful function if bound to a key. */
int blob_under_pointer_to_screen_centre() {

   int r = 0;
   int imol_map = imol_refinement_map();
   graphics_info_t g;
   if (imol_map != -1) {
      // OK we have a map to search.
      coot::Cartesian front = unproject(0.0);
      coot::Cartesian back  = unproject(1.0);
      clipper::Coord_orth p1(front.x(), front.y(), front.z());
      clipper::Coord_orth p2( back.x(),  back.y(),  back.z());
      try { 
	 clipper::Coord_orth blob =
	    graphics_info_t::molecules[imol_map].find_peak_along_line_favour_front(p1, p2);
	 coot::Cartesian cc(blob.x(), blob.y(), blob.z());
	 g.setRotationCentre(cc);
	 for(int ii=0; ii<graphics_info_t::n_molecules(); ii++) {
	    graphics_info_t::molecules[ii].update_map();
	    graphics_info_t::molecules[ii].update_symmetry();
	 }
	 g.make_pointer_distance_objects();
	 graphics_draw();
      }
      catch (std::runtime_error mess) {
	 std::cout << mess.what() << std::endl;
      }
   } else {
      std::string s = "Refinement map not selected - no action";
      std::cout << s << std::endl;
      // add_status_bar_text(s.c_str());
      info_dialog(s.c_str());
   }
   return r;
}


#ifdef USE_GUILE
/*! return the chain id of the water chain from a shelx molecule.  Raw interface
  Return #f if no chain or bad imol*/
SCM water_chain_from_shelx_ins_scm(int imol) {
   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      CChain *water_chain =
	 graphics_info_t::molecules[imol].water_chain_from_shelx_ins();
      if (water_chain) {
	 r = scm_makfrom0str(water_chain->GetChainID());
      } 
   }
   return r;
}
/*! return the chain id of the water chain. Raw interface */
SCM water_chain_scm(int imol) {
   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      CChain *water_chain = graphics_info_t::molecules[imol].water_chain();
      if (water_chain) {
	 r = scm_makfrom0str(water_chain->GetChainID());
      } 
   }
   return r;
}
#endif 

#ifdef USE_PYTHON
/*! return the chain id of the water chain from a shelx molecule.  Raw interface.
Return False if no chain or bad imol*/
PyObject *water_chain_from_shelx_ins_py(int imol) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      CChain *water_chain =
	 graphics_info_t::molecules[imol].water_chain_from_shelx_ins();
      if (water_chain) {
	 r = PyString_FromString(water_chain->GetChainID());
      } 
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
/*! return the chain id of the water chain. Raw interface */
PyObject *water_chain_py(int imol) {
   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      CChain *water_chain = graphics_info_t::molecules[imol].water_chain();
      if (water_chain) {
	 r = PyString_FromString(water_chain->GetChainID());
      } 
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif 

/*! \brief move waters of molecule number imol so that they are around the protein.

@return the number of moved waters. */
int move_waters_to_around_protein(int imol) {

   int r = 0;
   if (is_valid_model_molecule(imol)) {
      r = graphics_info_t::molecules[imol].move_waters_to_around_protein();
      std::cout << "INFO:: moved " << r << " water molecules" << std::endl;
      graphics_draw();
   } 
   return r;
}

/*! \brief return the maximum minimum distance of any water atom to
  any protein atom - used in validation of
  move_waters_to_around_protein() funtion.*/
float max_water_distance(int imol) {

   float f = -1;

   if (is_valid_model_molecule(imol)) {
      f = graphics_info_t::molecules[imol].max_water_distance();
   } 
   return f;
} 



// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//                       this is water validation (really)
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//
#ifdef USE_GUILE
/*! \brief return an improper list first is list of metals, second is
  list of waters that are coordinated with at least
  coordination_number of other atoms at distances less than or equal
  to dist_max */
// (list (list atom-spec-central (list atom-spec-contactor-1 atom-spec-contactor-2 ..)))
// 
SCM highly_coordinated_waters_scm(int imol, int coordination_number, float dist_max) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      CMMDBManager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      bool mol_has_symmetry = coot::mol_has_symmetry(mol);
      if (mol_has_symmetry) {
	 CMMDBManager *new_mol = coot::util::copy_molecule(mol);
	 coot::util::translate_close_to_origin(new_mol);
	 mol = new_mol; // do water coordination check with a molecule
			// that has been translated as close as
			// possible to the origin.  mol needs to be
			// deleted at the end of this function.
      }

      // this one happens for free.
      //
      coot::util::water_coordination_t wc_metals(mol, 4.0); // 4.0 magic number from Nayal and Di Cera 1996
      std::vector<std::pair<coot::util::contact_atoms_info_t, coot::util::contact_atoms_info_t::ele_index_t> >
	 metals = wc_metals.metals();

      std::cout << "Found " << metals.size() << " Metals (Na+, Li+, K+, Mg2+, Ca2+) "
		<< "amongst the waters" << std::endl;
      for (unsigned int i=0; i<metals.size(); i++) { 
	 std::cout << metals[i].first.central_atom() << " .... type: " << metals[i].second << std::endl;
      }
	 
      coot::util::water_coordination_t wc;
      SCM metal_results = SCM_EOL;
      for (unsigned int i=0; i<metals.size(); i++) {
	 std::string ele = coot::util::contact_atoms_info_t::ele_to_string(metals[i].second);
	 SCM metal_str_scm = scm_from_locale_string(ele.c_str());
	 SCM metal_results_ele = SCM_EOL;
	 SCM spec = atom_spec_to_scm(coot::atom_spec_t(metals[i].first.central_atom()));
	 metal_results_ele = scm_cons(metal_str_scm, metal_results_ele);
	 metal_results_ele = scm_cons(spec, metal_results_ele);
	 metal_results = scm_cons(metal_results_ele, metal_results);
      }
      metal_results = scm_reverse(metal_results);
      

      if (dist_max < 4.0) {
	 wc = wc_metals;
      } else {
	 wc = coot::util::water_coordination_t(mol, dist_max);
      }
      std::vector<coot::util::contact_atoms_info_t> water_contacts = 
	 wc.get_highly_coordinated_waters(coordination_number, dist_max);
      r = SCM_EOL; // a list (at least) because we didn't fail.
      for (unsigned int j=0; j<water_contacts.size(); j++) {
	 SCM atom_and_neighbours = SCM_EOL; // combines atom_spec_central and it neighbours
	 SCM atom_spec_central_scm = atom_spec_to_scm(water_contacts[j].central_atom());
	 SCM neighbours = SCM_EOL;
	 for (unsigned int k=0; k<water_contacts[j].size(); k++) {
	    coot::util::contact_atoms_info_t::contact_atom_t at = water_contacts[j][k];
	    if (at.dist < dist_max) {
	       SCM contactor_scm = atom_spec_to_scm(at.at);
	       neighbours = scm_cons(contactor_scm, neighbours);
	    }
	 }
	 atom_and_neighbours = scm_cons(scm_reverse(neighbours), atom_and_neighbours);
	 atom_and_neighbours = scm_cons(atom_spec_central_scm, atom_and_neighbours);
	 r = scm_cons(atom_and_neighbours, r);
      }
      r = scm_reverse(r);
      if (mol_has_symmetry)
	 delete mol; // it was a copy.
      return scm_cons(metal_results, r);
   }
   return SCM_BOOL_F;

}
#endif


#ifdef USE_PYTHON
/*! \brief return a list of waters that are coordinated with at least
  coordination_number of other atoms at distances less than or equal
  to dist_max */
PyObject *highly_coordinated_waters_py(int imol, int coordination_number, float dist_max) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
     CMMDBManager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
     coot::util::water_coordination_t wc(mol, dist_max);
     std::vector<coot::util::contact_atoms_info_t> water_contacts = 
       wc.get_highly_coordinated_waters(coordination_number, dist_max);
     r = PyList_New(0); // a list (at least) because we didn't fail.
     for (unsigned int j=0; j<water_contacts.size(); j++) {
       PyObject *atom_and_neighbours = PyList_New(2);
       PyObject *atom_spec_central_py = atom_spec_to_py(water_contacts[j].central_atom());
       PyObject *neighbours = PyList_New(0);
       for (unsigned int k=0; k<water_contacts[j].size(); k++) {
         coot::util::contact_atoms_info_t::contact_atom_t at = water_contacts[j][k];
         if (at.dist < dist_max) {
           PyObject *contactor_py = atom_spec_to_py(at.at);
           PyList_Append(neighbours, contactor_py);
         }
       }
       PyList_SetItem(atom_and_neighbours, 0, atom_spec_central_py);
       PyList_SetItem(atom_and_neighbours, 1, neighbours);
       PyList_Append(r, atom_and_neighbours);
     }
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
} 
#endif
