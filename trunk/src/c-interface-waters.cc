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

#ifdef USE_GUILE
#include <guile/gh.h>
#endif // USE_GUILE

// Including python needs to come after graphics-info.h, because
// something in Python.h (2.4 - chihiro) is redefining FF1 (in
// ssm_superpose.h) to be 0x00004000 (Grrr).
//
#ifdef USE_PYTHON
#include "Python.h"
#endif // USE_PYTHON

#include "c-interface.h"

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

