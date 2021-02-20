/* src/c-interface-ncs.cc
 * 
 * Copyright 2004, 2005 The University of York
 * Copyright 2007, 2008 The University of Oxford
 * Copyright 2015, 2016 by Medical Research Council
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
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#include "python-3-interface.hh"
#endif

#include "compat/coot-sysdep.h"


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
		       // functions are built by glade.

#include <vector>
#include <string>
#include <algorithm>

#include <mmdb2/mmdb_manager.h>
#include "coords/mmdb-extras.h"
#include "coords/mmdb.h"
#include "coords/mmdb-crystal.h"

#include "graphics-info.h"

#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "cc-interface.hh"

#include "guile-fixups.h"

// Including python needs to come after graphics-info.h, because
// something in Python.h (2.4 - chihiro) is redefining FF1 (in
// ssm_superpose.h) to be 0x00004000 (Grrr).
//
// 20100813: Python.h needs to come before to stop"_POSIX_C_SOURCE" redefined problems 
//
// #ifdef USE_PYTHON
// #include "Python.h"
// #endif // USE_PYTHON

int add_strict_ncs_matrix(int imol,
			  const char *this_chain_id,
			  const char *target_chain_id,
			  float m11, float m12, float m13, 
			  float m21, float m22, float m23, 
			  float m31, float m32, float m33, 
			  float t1,  float t2,  float t3)
{

   int istat = 0;
   if (is_valid_model_molecule(imol)) {
      clipper::Mat33<double> m(m11, m12, m13,
			       m21, m22, m23,
			       m31, m32, m33);
      clipper::Coord_orth c(t1, t2, t3);
      clipper::RTop_orth rtop(m, c);

      coot::coot_mat44 cm44;
      // The orientation is a guess: FIXME or testme.
      cm44.m[0].v4[0] = m11;
      cm44.m[0].v4[1] = m12;
      cm44.m[0].v4[2] = m13;
      cm44.m[1].v4[0] = m21;
      cm44.m[1].v4[1] = m22;
      cm44.m[1].v4[2] = m23;
      cm44.m[2].v4[0] = m31;
      cm44.m[2].v4[1] = m32;
      cm44.m[2].v4[2] = m33;
      // translation
      cm44.m[0].v4[3] = t1;
      cm44.m[1].v4[3] = t2;
      cm44.m[2].v4[3] = t3;
      // sensibles
      cm44.m[3].v4[0] = 0.0;
      cm44.m[3].v4[1] = 0.0;
      cm44.m[3].v4[2] = 0.0;
      cm44.m[3].v4[3] = 1.0;
      
      istat = 1;
      std::string tch = target_chain_id;
      std::string chain_id = this_chain_id;

      graphics_info_t::molecules[imol].add_strict_ncs_matrix(chain_id, tch, cm44);

      graphics_draw();
   }
   return istat;
}

int add_strict_ncs_from_mtrix_from_self_file(int imol) {

   int istat = 0; 
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      g.molecules[imol].add_strict_ncs_from_mtrix_from_self_file();
   }
   return istat;
}



GtkWidget *wrapped_create_ncs_maps_dialog() {

   GtkWidget *dialog = create_ncs_maps_dialog();
   short int diff_maps_only_flag = 0;
   int ifound;

   // Maps: 
   ifound = fill_ligands_dialog_map_bits_by_dialog_name(dialog,
							"ncs_maps_maps",
							diff_maps_only_flag);
   if (ifound == 0) {
      std::cout << "Error: you must have a difference map to analyse!" << std::endl;
      GtkWidget *none_frame = lookup_widget(dialog, "no_maps_frame");
      gtk_widget_show(none_frame);
   }

   // Models:
   short int have_ncs_flag = 1;
   ifound = fill_vbox_with_coords_options_by_dialog_name(dialog, "ncs_maps_models",
							 have_ncs_flag);
   if (ifound == 0) {
      std::cout << "You must have molecules with NCS to use this function\n";
      GtkWidget *none_frame = lookup_widget(dialog, "no_models_frame");
      gtk_widget_show(none_frame);
   }


   return dialog;
}

void add_ncs_matrix(int imol,
		    const char *chain_id, 
		    const char *target_chain_id, 
		    float m11, float m12, float m13, 
		    float m21, float m22, float m23, 
		    float m31, float m32, float m33, 
		    float t1,  float t2,  float t3) {

   if (is_valid_model_molecule(imol)) {
      clipper::Mat33<double> m(m11, m12, m13,
			       m21, m22, m23,
			       m31, m32, m33);
      clipper::Coord_orth c(t1, t2, t3);
      clipper::RTop_orth rtop(m, c);

//       std::cout  << "DEBUG:: add_ncs_matrix for imol " << imol << " chain_id " << chain_id
// 		 << " target_chain_id " << target_chain_id << std::endl;

      graphics_info_t::molecules[imol].add_ncs_ghost(chain_id, target_chain_id, rtop);
      graphics_draw();
   }
}

void clear_ncs_ghost_matrices(int imol) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].clear_ncs_ghost_matrices();
   }
} 


int show_strict_ncs_state(int imol) {

   if (is_valid_model_molecule(imol)) {
      return graphics_info_t::molecules[imol].show_strict_ncs_flag;
   }
   return 0;
}

void set_show_strict_ncs(int imol, int state) {

   short int istate = 0;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].show_strict_ncs_flag = state;
      istate = 1;
      graphics_draw();
   }
}


/* At what level of homology should we say that we can't see homology
   for NCS calculation? (default 0.8) */
void set_ncs_homology_level(float flev) {

   graphics_info_t::ncs_homology_level = flev;

}

/*  ----------------------------------------------------------------------- */
/*                  NCS                                                     */
/*  ----------------------------------------------------------------------- */
/* section NCS */
void set_draw_ncs_ghosts(int imol, int istate) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_show_ghosts(istate);
      graphics_draw();
   }
}

/*! \brief return the drawing state of NCS ghosts for molecule number imol   */
int draw_ncs_ghosts_state(int imol) {

   int r = -1;
   if (is_valid_model_molecule(imol)) {
      r = graphics_info_t::molecules[imol].draw_ncs_ghosts_p();
   }
   return r;
} 


void set_ncs_ghost_bond_thickness(int imol, float f) {

   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_ghost_bond_thickness(f);
      graphics_draw();
   }
}


void ncs_update_ghosts(int imol) {

   if (is_valid_model_molecule(imol)) {
      int incs = graphics_info_t::molecules[imol].update_ncs_ghosts();
      if (incs)
	 graphics_draw();
   }
}

GtkWidget *wrapped_create_ncs_control_dialog() {

   GtkWidget *w = create_ncs_control_dialog();

   for (int imol=0; imol<graphics_info_t::n_molecules(); imol++)
      if (is_valid_model_molecule(imol))
	 graphics_info_t::molecules[imol].fill_ncs_control_frame(w);
   return w; 
}


void ncs_control_change_ncs_master_to_chain(int imol, int ichain) {

   std::cout << "DEBUG:: ncs_control_change_ncs_master_to_chain imol: " << imol
	     << " and ichain: " << ichain << std::endl;
   if (is_valid_model_molecule(imol)) {
      std::vector<std::string> chain_ids =
	 coot::util::chains_in_molecule(graphics_info_t::molecules[imol].atom_sel.mol);
      if (ichain < int(chain_ids.size()))
	 graphics_info_t::molecules[imol].set_ncs_master_chain(chain_ids[ichain], graphics_info_t::ncs_homology_level);
      graphics_draw();
   } 
}

void ncs_control_change_ncs_master_to_chain_id(int imol, const char *chain_id) {

   std::cout << "DEBUG ncs_control_change_ncs_master_to_chain_id imol: " << imol
	     << " and chain_id: " << chain_id << std::endl;
   if (is_valid_model_molecule(imol)) {
     std::vector<std::string> chain_ids =
       coot::util::chains_in_molecule(graphics_info_t::molecules[imol].atom_sel.mol);
     std::vector<std::string>::iterator match = find(chain_ids.begin(), chain_ids.end(), chain_id);
     if (match != chain_ids.end())
	graphics_info_t::molecules[imol].set_ncs_master_chain(chain_id, graphics_info_t::ncs_homology_level);
     graphics_draw();
   }
}

void ncs_control_change_ncs_master_to_chain_update_widget(GtkWidget *w, int imol, int ichain) {

   std::cout << "DEBUG ncs_control_change_ncs_master_to_chain_update_widget imol: " << imol
	     << " and ichain: " << ichain << std::endl;
   if (is_valid_model_molecule(imol)) {
      ncs_control_change_ncs_master_to_chain(imol, ichain); 

      // Now we want to update the widget.  We need to change the sensitivity of
      // all the Chain check boxes in the dispaly ncs chain vbox.
      // 
      // We need to change to desensitve the chain that matches ichain.
      //
      // This is a -widget-work function
      graphics_info_t::molecules[imol].ncs_control_change_ncs_master_to_chain_update_widget(w, ichain);
   } 
}


void ncs_control_display_chain(int imol, int ichain, int state) {

   std::cout << "%%%% ncs_control_display_chain " << std::endl;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t::molecules[imol].set_display_ncs_ghost_chain(ichain, state);
      graphics_draw();
   }
}

// The NCS "Yes" was pressed, do we need to make the rtops for this
// molecule?  We do if there ncs_ghosts.size() > 0.
//
// If ncs_ghosts.size() == 0, then this widget should not be active.
// 
// 
void
make_ncs_ghosts_maybe(int imol) {

   // c.f. fill_bond_parameters_internals()
   
   // int imol = graphics_info_t::bond_parameters_molecule; old way.  broken.
   
   if (is_valid_model_molecule(imol)) {  // it should be!
      if (graphics_info_t::molecules[imol].has_ncs_p()) {
	 if (graphics_info_t::molecules[imol].ncs_ghosts_have_rtops_p() == 0) {
	    graphics_info_t::molecules[imol].fill_ghost_info(1, graphics_info_t::ncs_homology_level);
	 }
      }
   } 
}

// 
int make_dynamically_transformed_ncs_maps(int imol_model, int imol_map, int overwrite_maps_of_same_name_flag) {

   int nmaps=0;

   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {
	 float homology_lev = graphics_info_t::ncs_homology_level;

	 graphics_info_t g;
	 if (g.molecules[imol_model].has_ncs_p() > 0) {
	    if (g.molecules[imol_model].ncs_ghosts_have_rtops_p() == 0) {
	       // fill the rtops and set the flag
	       g.molecules[imol_model].fill_ghost_info(1, homology_lev); 
	    }
	 }

	 std::vector<coot::ghost_molecule_display_t> local_ncs_ghosts =
	    g.molecules[imol_model].NCS_ghosts();
	 for (unsigned int ighost=0; ighost<local_ncs_ghosts.size(); ighost++) {
	    int imol = -1;// gets reset
	    std::string ncs_map_name = "Map ";
	    ncs_map_name += coot::util::int_to_string(imol_map);
	    ncs_map_name += " ";
	    ncs_map_name += local_ncs_ghosts[ighost].name;
	    
	    if (overwrite_maps_of_same_name_flag) { 
	       int imap_ghost_existing = g.lookup_molecule_name(ncs_map_name);

	       if (is_valid_map_molecule(imap_ghost_existing)) {
		  imol = imap_ghost_existing;
	       } else {
		  imol = graphics_info_t::create_molecule();
	       }
	    } else {
	       imol = graphics_info_t::create_molecule();
	    } 
	    g.molecules[imol].install_ghost_map(g.molecules[imol_map].xmap,
						ncs_map_name, // was local_ncs_ghosts[ighost].name?
						local_ncs_ghosts[ighost],
						g.molecules[imol_map].is_difference_map_p(),
						g.swap_difference_map_colours,
						g.molecules[imol_map].get_contour_level());
	    nmaps++;
	 }

	 if (graphics_info_t::ncs_maps_do_average_flag) {
        std:: string imol_map_name = coot::util::int_to_string(imol_map);
	    std::vector<std::pair<clipper::Xmap<float>, std::string> > xmaps  = 
          g.molecules[imol_model].ncs_averaged_maps(g.molecules[imol_map].xmap, homology_lev, imol_map_name);
	    std::cout << "INFO:: made " << xmaps.size() << " averaged map(s)" << std::endl;
	    for (unsigned int i=0; i<xmaps.size(); i++) {
	       std::string name;
	       name += xmaps[i].second;
	       int imol = -1;// gets reset
	       int imap_ghost_existing = g.lookup_molecule_name(name);
	       if (overwrite_maps_of_same_name_flag) { 
		  if (is_valid_map_molecule(imap_ghost_existing)) {
		     imol = imap_ghost_existing;
		  } else {
		     imol = graphics_info_t::create_molecule();
		  }
	       } else {
		  imol = graphics_info_t::create_molecule();
	       }

	       bool is_em_flag = graphics_info_t::molecules[imol_map].is_EM_map();
	       g.molecules[imol].install_new_map(xmaps[i].first, name, is_em_flag);
	       if (g.molecules[imol_map].is_difference_map_p())
		  g.molecules[imol].set_map_is_difference_map();
	       nmaps++;
	    } 
	 }
	 
      } else {
	 std::cout << "WARNING:: molecule number " << imol_map
		   << " is not a valid map molecule\n";
      }
   } else {
      std::cout << "WARNING:: molecule number " << imol_model
		<< " is not a valid model molecule\n";
   }

   if (nmaps > 0)
      graphics_draw();
   
   return nmaps;
}

// Also this should go elsewhere
// 
int make_dynamically_transformed_ncs_maps_by_widget(GtkWidget *dialog) {

   // Decode the options and call the above function:

   int imol_map = -1;
   int imol_coords = -1;
 
   GtkWidget *map_button;
   short int found_active_button_for_map = 0;
   for (int imol=0; imol<graphics_info_t::n_molecules(); imol++) {
      if (graphics_info_t::molecules[imol].has_xmap()) {
	 std::string map_str = "ncs_maps_maps_radiobutton_";
	 map_str += graphics_info_t::int_to_string(imol);
	 map_button = lookup_widget(dialog, map_str.c_str());
	 if (map_button) {
	    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(map_button))) { 
	       imol_map = imol;
	       found_active_button_for_map = 1;
	       break;
	    }
	 } else {
	    std::cout << "WARNING:: (error) " << map_str << " map button not found in "
		      << "make_dynamically_transformed_ncs_maps_by_widget" << std::endl;
	 }
      }
   }

   GtkWidget *coords_button;
   short int found_active_button_for_coords = 0;
   for (int imol=0; imol<graphics_info_t::n_molecules(); imol++) {
      if (graphics_info_t::molecules[imol].has_model()) {
	 if (graphics_info_t::molecules[imol].has_ncs_p()) {
	    std::string coords_str = "ncs_maps_models_radiobutton_";
	    coords_str += graphics_info_t::int_to_string(imol);
	    coords_button = lookup_widget(dialog, coords_str.c_str());
	    if (coords_button) {
	       if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(coords_button))) {
		  imol_coords = imol;
		  found_active_button_for_coords = 1;
	       }
	    } else {
	       std::cout << "WARNING:: (error) " << coords_str
			 << " coords button not found in "
			 << "make_dynamically_transformed_ncs_maps_by_widget"
			 << std::endl;
	    }
	 }
      }
   }

   // Use a button to get this value at some stage.
   int overwrite_maps_of_same_name_flag = 0;

   if (!found_active_button_for_coords) {
      std::cout << "You need to define a set of coordinates for NCS maping\n";
   } else {
      if (!found_active_button_for_map) {
	 std::cout << "You need to define a map for NCS maping\n";
      } else {
	 make_dynamically_transformed_ncs_maps(imol_coords, imol_map,
					       overwrite_maps_of_same_name_flag);
      }
   }

   return 0;
}


// Should be in c-interface-em.cc, perhaps?
int scale_cell(int imol_map, float fac_u, float fac_v, float fac_w) {

   int retval = 0;
   if (is_valid_map_molecule(imol_map)) {
      retval = graphics_info_t::molecules[imol_map].scale_cell(fac_u, fac_v, fac_w);
      graphics_draw();
   }
   return retval; 
}


#ifdef USE_GUILE
SCM ncs_chain_differences_scm(int imol, const char *master_chain_id) {

   float mc_weight = 1.0;
   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      coot::ncs_differences_t diffs = 
	 graphics_info_t::molecules[imol].ncs_chain_differences(master_chain_id,
								mc_weight);
      if (diffs.size() == 0) {
	 std::cout << "no diffs" << std::endl;
      } else {
	 r = SCM_EOL;
	 for (int idiff=int(diffs.size()-1); idiff>=0; idiff--) {
	    SCM l_residue_data = SCM_EOL;
	    coot::ncs_chain_difference_t cd = diffs.diffs[idiff];
	    if (cd.residue_info.size() > 0) {
	       std::cout << "NCS target chain has " << cd.residue_info.size()
			 << " peers." << std::endl;
	       //	       for (int iresinf=0; iresinf<cd.residue_info.size(); iresinf++) {
	       for (int iresinf=(int(cd.residue_info.size())-1); iresinf>=0; iresinf--) {
		  if (0)
		     std::cout << "resinfo: "
			       << cd.residue_info[iresinf].resno << " "
			       << cd.residue_info[iresinf].inscode << " "
			       << cd.residue_info[iresinf].serial_number << " to "
			       << cd.residue_info[iresinf].target_resno << " "
			       << cd.residue_info[iresinf].target_inscode << " "
			       << cd.residue_info[iresinf].target_serial_number << " diff: "
			       << cd.residue_info[iresinf].mean_diff
			       << std::endl;
		  coot::residue_spec_t this_res(cd.peer_chain_id,
						cd.residue_info[iresinf].resno,
						cd.residue_info[iresinf].inscode);
		  coot::residue_spec_t target_res(diffs.target_chain_id,
						  cd.residue_info[iresinf].target_resno,
						  cd.residue_info[iresinf].target_inscode);
		  SCM res_l = SCM_EOL;
		  res_l = scm_cons(scm_double2num(cd.residue_info[iresinf].mean_diff), res_l);
//		  res_l = scm_cons(scm_cdr(scm_residue(target_res)), res_l);
//		  res_l = scm_cons(scm_cdr(scm_residue(this_res)), res_l);
		  l_residue_data = scm_cons(res_l, l_residue_data);
	       }
	       r = scm_cons(l_residue_data, SCM_EOL);
	       r = scm_cons(scm_makfrom0str(diffs.target_chain_id.c_str()), r);
	       r = scm_cons(scm_makfrom0str(cd.peer_chain_id.c_str()), r);
	    }
	 }
      }
   }
   return r;
}
#endif	/* USE_GUILE */
#ifdef USE_PYTHON
PyObject *ncs_chain_differences_py(int imol, const char *master_chain_id) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      coot::ncs_differences_t diffs = 
	 graphics_info_t::molecules[imol].ncs_chain_differences(master_chain_id, 1.0);
      if (diffs.size() == 0) {
	 std::cout << "no diffs" << std::endl;
      } else {
	 r = PyList_New(0);
	 for (unsigned int idiff=0; idiff<diffs.size(); idiff++) {
	    PyObject *l_residue_data = PyList_New(0);
	    coot::ncs_chain_difference_t cd = diffs.diffs[idiff];
	    if (cd.residue_info.size() > 0) {
	       std::cout << "NCS target chain has " << cd.residue_info.size()
			 << " peers." << std::endl;
	       //	       for (int iresinf=0; iresinf<cd.residue_info.size(); iresinf++) {
	       for (unsigned int iresinf=0; iresinf<cd.residue_info.size(); iresinf++) {
		  if (0)
		     std::cout << "resinfo: "
			       << cd.residue_info[iresinf].resno << " "
			       << cd.residue_info[iresinf].inscode << " "
			       << cd.residue_info[iresinf].serial_number << " to "
			       << cd.residue_info[iresinf].target_resno << " "
			       << cd.residue_info[iresinf].target_inscode << " "
			       << cd.residue_info[iresinf].target_serial_number << " diff: "
			       << cd.residue_info[iresinf].mean_diff
			       << std::endl;
		  coot::residue_spec_t this_res(cd.peer_chain_id,
						cd.residue_info[iresinf].resno,
						cd.residue_info[iresinf].inscode);
		  coot::residue_spec_t target_res(diffs.target_chain_id,
						  cd.residue_info[iresinf].target_resno,
						  cd.residue_info[iresinf].target_inscode);
		  // according to Paul's documentation we should have resno and inscode
		  // for both residues here too
		  PyObject *thisr = PyList_GetSlice(residue_spec_to_py(this_res), 2, 4);
		  
		  PyObject *masta = PyList_GetSlice(residue_spec_to_py(target_res), 2, 4);

		  // res_l list seems only to have one element in Paul's scm code
		  // currently?! Correct?!
		  PyObject *res_l = PyList_New(3);
		  PyList_SetItem(res_l, 0, thisr);
		  PyList_SetItem(res_l, 1, masta);
		  PyList_SetItem(res_l, 2, PyFloat_FromDouble(cd.residue_info[iresinf].mean_diff));
		  PyList_Append(l_residue_data, res_l);
		  Py_XDECREF(res_l);
	       }
	       PyList_Append(r, myPyString_FromString(cd.peer_chain_id.c_str()));
	       PyList_Append(r, myPyString_FromString(diffs.target_chain_id.c_str()));
	       PyList_Append(r, l_residue_data);
	       Py_XDECREF(l_residue_data);
	    }
	 }
      }
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif	/* USE_PYTHON */


#ifdef USE_GUILE
SCM ncs_chain_ids_scm(int imol) {
   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      if (graphics_info_t::molecules[imol].has_ncs_p()) {
	 std::vector<std::vector<std::string> > ncs_ghost_chains =
	    graphics_info_t::molecules[imol].ncs_ghost_chains();
// 	 std::cout << "There are " << ncs_ghost_chains.size() << " ncs ghost chains"
// 		<< std::endl;
	 if (ncs_ghost_chains.size() > 0) {
	    r = SCM_EOL;
	    for (int i=(int(ncs_ghost_chains.size())-1); i>=0; i--) {
	       SCM string_list_scm =
		  generic_string_vector_to_list_internal(ncs_ghost_chains[i]);
	       r = scm_cons(string_list_scm, r);
	    }
	 }
      }
   }
   return r;
}
#endif	/* USE_GUILE */

#ifdef USE_PYTHON
PyObject *ncs_chain_ids_py(int imol) {
   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      if (graphics_info_t::molecules[imol].has_ncs_p()) {
	 std::vector<std::vector<std::string> > ncs_ghost_chains =
	    graphics_info_t::molecules[imol].ncs_ghost_chains();
// 	 std::cout << "There are " << ncs_ghost_chains.size() << " ncs ghost chains"
// 		<< std::endl;
	 if (ncs_ghost_chains.size() > 0) {
	    r = PyList_New(ncs_ghost_chains.size());
	    for (unsigned int i=0; i<ncs_ghost_chains.size(); i++) {
	       PyObject *string_list_py =
		  generic_string_vector_to_list_internal_py(ncs_ghost_chains[i]);
	       PyList_SetItem(r, i, string_list_py);
	    }
	 }
      }
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif	/* USE_PYTHON */


#ifdef USE_GUILE
/* return #f on bad imol or a list of ghosts on good imol.  Can
   include NCS rtops if they are available, else the rtops are #f */
SCM ncs_ghosts_scm(int imol) {
   SCM r = SCM_BOOL_F;

   if (!is_valid_model_molecule(imol)) {
      std::cout << "WARNING:: molecule number " << imol << " is not valid"
		<< std::endl;
   } else {
      r = SCM_EOL;
      std::vector<coot::ghost_molecule_display_t> ncs_ghosts =
	 graphics_info_t::molecules[imol].NCS_ghosts();
      for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
	 SCM ghost_scm = SCM_EOL;
	 SCM display_it_flag_scm = SCM_BOOL_F;
	 if (ncs_ghosts[ighost].display_it_flag)
	    display_it_flag_scm = SCM_BOOL_T;
	 SCM rtop_scm = SCM_BOOL_F;
	 if (graphics_info_t::molecules[imol].ncs_ghosts_have_rtops_p())
	    rtop_scm = rtop_to_scm(ncs_ghosts[ighost].rtop);
	 SCM target_chain_id_scm = scm_makfrom0str(ncs_ghosts[ighost].target_chain_id.c_str());
	 SCM chain_id_scm = scm_makfrom0str(ncs_ghosts[ighost].chain_id.c_str());
	 SCM name_scm = scm_makfrom0str(ncs_ghosts[ighost].name.c_str());
	 
	 ghost_scm = scm_cons(display_it_flag_scm, ghost_scm);
	 ghost_scm = scm_cons(rtop_scm,            ghost_scm);
	 ghost_scm = scm_cons(target_chain_id_scm, ghost_scm);
	 ghost_scm = scm_cons(chain_id_scm,        ghost_scm);
	 ghost_scm = scm_cons(name_scm,            ghost_scm);
	 r = scm_cons(ghost_scm, r);
      }
      r = scm_reverse(r);
   }
   return r;
}
#endif	/* USE_GUILE */

#ifdef USE_PYTHON
/* return False on bad imol or a list of ghosts on good imol.  Can
   include NCS rtops if they are available, else the rtops are False */
PyObject *ncs_ghosts_py(int imol) {

   PyObject *r = Py_False;
   if (!is_valid_model_molecule(imol)) {
      std::cout << "WARNING:: molecule number " << imol << " is not valid"
		<< std::endl;
   } else {
     r = PyList_New(0);
     std::vector<coot::ghost_molecule_display_t> ncs_ghosts =
       graphics_info_t::molecules[imol].NCS_ghosts();
     for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
       PyObject *ghost_py = PyList_New(5);
       PyObject *display_it_flag_py = Py_False;
       if (ncs_ghosts[ighost].display_it_flag)
	 display_it_flag_py = Py_True;
       PyObject *rtop_py = Py_False;
       if (graphics_info_t::molecules[imol].ncs_ghosts_have_rtops_p()) {
	 rtop_py = rtop_to_python(ncs_ghosts[ighost].rtop);
       }
       PyObject *target_chain_id_py = myPyString_FromString(ncs_ghosts[ighost].target_chain_id.c_str());
       PyObject *chain_id_py = myPyString_FromString(ncs_ghosts[ighost].chain_id.c_str());
       PyObject *name_py = myPyString_FromString(ncs_ghosts[ighost].name.c_str());

       PyList_SetItem(ghost_py, 0, name_py);
       PyList_SetItem(ghost_py, 1, chain_id_py);
       PyList_SetItem(ghost_py, 2, target_chain_id_py);
       PyList_SetItem(ghost_py, 3, rtop_py);
       PyList_SetItem(ghost_py, 4, display_it_flag_py);
       PyList_Append(r, ghost_py);
       Py_XDECREF(ghost_py);
     }
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif	/* USE_PYTHON */


// This should be  in c-interface-ncs-gui.cc
void validation_graph_ncs_diffs_mol_selector_activate (GtkMenuItem     *menuitem,
						      gpointer         user_data) {
   
   int imol = GPOINTER_TO_INT(user_data);
   graphics_info_t g;
   // g.ncs_diffs_from_mol(imol);
   std::cout << "fixme in validation_graph_ncs_diffs_mol_selector_activate() " << std::endl;

}

void
set_ncs_matrix_type(int flag) {

   graphics_info_t g;
   if (flag == coot::NCS_SSM) {
#ifdef HAVE_SSMLIB
      g.ncs_matrix_flag = coot::NCS_SSM;
#else
      std::cout<< "WARNING:: not compiled with SSM, so will use LSQ for NCS matrix determination" <<std::endl;
      g.ncs_matrix_flag = coot::NCS_LSQ;
#endif // HAVE_SSMLIB
   } else {
      if (flag == coot::NCS_LSQ2) {
	g.ncs_matrix_flag = coot::NCS_LSQ2;
      } else {
        g.ncs_matrix_flag = coot::NCS_LSQ;
      }
   }
}

int
get_ncs_matrix_state() {

   graphics_info_t g;
   return g.ncs_matrix_flag;
}

/*! \brief Given that we are in chain current_chain, apply the NCS
  operator that maps current_chain on to next_ncs_chain, so that the
  relative view is preserved.  For NCS skipping. */
void apply_ncs_to_view_orientation(int imol, const char *current_chain, const char *next_ncs_chain) {

   if (is_valid_model_molecule(imol)) {

#if 0

      // we don't use quat - but we do use glm::quat glm_quat - so this
      // function could be rewritten one day.

      short int forward_flag = 1; // emulate previous behaviour.  Not
				  // sure that this is what is needed.
      coot::util::quaternion q(graphics_info_t::quat[0],
			       graphics_info_t::quat[1],
			       graphics_info_t::quat[2],
			       graphics_info_t::quat[3]);
      clipper::Mat33<double> current_view_mat = q.matrix();
      clipper::Coord_orth current_centre(graphics_info_t::RotationCentre_x(), 
					 graphics_info_t::RotationCentre_y(),
					 graphics_info_t::RotationCentre_z());
      std::pair<bool, clipper::RTop_orth> new_ori = 
	 graphics_info_t::molecules[imol].apply_ncs_to_view_orientation(current_view_mat,
									current_centre,
									current_chain,
									next_ncs_chain,
									forward_flag);

      std::cout << "DEBUG::   NCS view in:  \n" << current_view_mat.format() << std::endl;

      std::cout << "DEBUG::   NCS view out: " << new_ori.first << std::endl;
      std::cout << "DEBUG::   NCS view out: \n" << new_ori.second.format() << "\n";
      
      if (new_ori.first) {
	 coot::util::quaternion vq(new_ori.second.rot());
	 graphics_info_t::quat[0] = vq.q0;
	 graphics_info_t::quat[1] = vq.q1;
	 graphics_info_t::quat[2] = vq.q2;
	 graphics_info_t::quat[3] = vq.q3;
      }
      graphics_draw();
#endif
   } 
}


/*! \brief Given that we are in chain current_chain, apply the NCS
  operator that maps current_chain on to next_ncs_chain, so that the
  relative view is preserved.  For NCS skipping. */
void apply_ncs_to_view_orientation_and_screen_centre(int imol,
						     const char *current_chain,
						     const char *next_ncs_chain,
						     short int forward_flag) {

   if (is_valid_model_molecule(imol)) {
      
#if 0

      // we don't use quat - but we do use glm::quat glm_quat - so this
      // function could be rewritten one day.

      coot::util::quaternion q(graphics_info_t::quat[0],
			       graphics_info_t::quat[1],
			       graphics_info_t::quat[2],
			       graphics_info_t::quat[3]);
      clipper::Coord_orth current_centre(graphics_info_t::RotationCentre_x(), 
					 graphics_info_t::RotationCentre_y(),
					 graphics_info_t::RotationCentre_z());
      clipper::Mat33<double> current_view_mat = q.matrix();
      std::pair<bool, clipper::RTop_orth> new_ori = 
	 graphics_info_t::molecules[imol].apply_ncs_to_view_orientation(current_view_mat,
									current_centre,
									current_chain,
									next_ncs_chain,
									forward_flag);

//       std::cout << "   NCS view in:  \n" << current_view_mat.format() << std::endl;
//       std::cout << "   NCS view out: " << new_ori.first << std::endl;
//       std::cout << "   NCS view out: \n" << new_ori.second.format() << "\n";
      
      if (new_ori.first) {
	 coot::util::quaternion vq(new_ori.second.rot());
	 graphics_info_t::quat[0] = vq.q0;
	 graphics_info_t::quat[1] = vq.q1;
	 graphics_info_t::quat[2] = vq.q2;
	 graphics_info_t::quat[3] = vq.q3;

	 clipper::Coord_orth new_centre(new_ori.second.trn());
	 graphics_info_t g;
	 g.setRotationCentre(coot::Cartesian(new_centre.x(),
					     new_centre.y(),
					     new_centre.z()));
	 g.update_things_on_move();
	 if (graphics_info_t::environment_show_distances) {
	    std::pair<int, int> r =  g.get_closest_atom();
	    // std::cout << "got closest atom " << r.first << " " << r.second << std::endl;
	    if (r.first >= 0) {
	       g.mol_no_for_environment_distances = r.second;
	       g.update_environment_distances_maybe(r.first, r.second);
	    }
	 }

      }
      graphics_draw();
#endif
   } 
}

#include "cc-interface-ncs.hh"

std::vector<int>
make_ncs_maps(int imol_model, int imol_map) {

   double border = 2.5; // A: the additional distance between the centre of the atom
			// and the most-extented atom (so that we capture the
			// density of the most-extented atom.
   std::vector<int> new_map_molecule_numbers;
   
   int status = -1;

   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {
	 
	 // get_model_ncs_info(): returns ghost info.

	 if (graphics_info_t::molecules[imol_model].ncs_ghosts_have_rtops_p() == 0) {
	    graphics_info_t::molecules[imol_model].fill_ghost_info(1, graphics_info_t::ncs_homology_level);
	 }

	 std::vector<coot::ghost_molecule_display_t> ncs_ghosts =
	    graphics_info_t::molecules[imol_model].NCS_ghosts();

	 if (! ncs_ghosts.empty()) {

	    // get the middle of and target chain_id

	    for (unsigned int ighst=0; ighst<ncs_ghosts.size(); ighst++) {
	       	       
	       std::pair<clipper::Coord_orth, double> ccr =
		  graphics_info_t::molecules[imol_model].chain_centre_and_radius(ncs_ghosts[ighst].target_chain_id);

	       const clipper::Coord_orth chain_centre = ccr.first;
	       double ncs_chain_radius = ccr.second + border; // if we are are the
							      // centre of a given
							      // chain_id, how big a
							      // radius do we need
							      // to encompass all
							      // atoms of that chain?
	       clipper::RTop_orth rt = ncs_ghosts[ighst].rtop;

	       std::cout << "rtop for ghost is \n" << rt.format() << std::endl;
	       const clipper::Cell &map_cell =
		  graphics_info_t::molecules[imol_map].xmap.cell();
	       const clipper::Spacegroup &space_group =
		  graphics_info_t::molecules[imol_map].xmap.spacegroup();
	       std::string space_group_name = space_group.symbol_xhm();

	       int imol_new = transform_map_raw(imol_map,
						rt.rot()(0,0), rt.rot()(0,1), rt.rot()(0,2),
						rt.rot()(1,0), rt.rot()(1,1), rt.rot()(1,2),
						rt.rot()(2,0), rt.rot()(2,1), rt.rot()(2,2),
						rt.trn()[0], 
						rt.trn()[1], 
						rt.trn()[2],
						chain_centre.x(),
						chain_centre.y(),
						chain_centre.z(),
						ncs_chain_radius,
						space_group_name.c_str(),
						map_cell.descr().a(),
						map_cell.descr().b(),
						map_cell.descr().c(),
						clipper::Util::rad2d(map_cell.descr().alpha()),
						clipper::Util::rad2d(map_cell.descr().beta()),
						clipper::Util::rad2d(map_cell.descr().gamma()));

	       new_map_molecule_numbers.push_back(imol_new);
	       
	    }
	 } 
      }
   } 

   return new_map_molecule_numbers;
}



void copy_from_ncs_master_to_others(int imol, const char *chain_id) {

   if (is_valid_model_molecule(imol)) {
      std::string c(chain_id);
      graphics_info_t::molecules[imol].copy_from_ncs_master_to_others(c);
      graphics_draw();
   }
}



void
copy_residue_range_from_ncs_master_to_others(int imol, 
					     const char *master_chain_id,
					     int res_range_start,
					     int res_range_end) {
   
   if (!is_valid_model_molecule(imol)) {
      std::cout << " molecule " << imol << " is not a valid model molecule"
		<< std::endl;
   } else {


      GtkWidget *ncs_dialog = 0; // we should store/look this up somehow
                                 // (using graphics_info_t::ncs_control_dialog, say)

      int ich = -1;
      std::vector<std::string> chain_ids =
	 coot::util::chains_in_molecule(graphics_info_t::molecules[imol].atom_sel.mol);
      for (unsigned int ii=0; ii<chain_ids.size(); ii++) { 
	 if (chain_ids[ii] == master_chain_id) {
	    ich = ii;
	    break;
	 } 
      }

      if (ich != -1) {
	 ncs_control_change_ncs_master_to_chain_update_widget(ncs_dialog, imol, ich);
	 std::string mc(master_chain_id);
	 graphics_info_t::molecules[imol].copy_residue_range_from_ncs_master_to_others(mc,
										       res_range_start,
										       res_range_end);
      }
      graphics_draw();
   }
}

#ifdef USE_GUILE
SCM ncs_master_chains_scm(int imol) {
   SCM r = SCM_BOOL_F;

   if (is_valid_model_molecule(imol)) {
      std::vector<std::string> v = graphics_info_t::molecules[imol].ncs_master_chains();
      if (v.size()) {
	 r = generic_string_vector_to_list_internal(v);
      } 
   } 
   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *ncs_master_chains_py(int imol) {

   PyObject *r = Py_False;

   if (is_valid_model_molecule(imol)) {
      std::vector<std::string> v = graphics_info_t::molecules[imol].ncs_master_chains();
      if (v.size()) {
	 r = generic_string_vector_to_list_internal_py(v);
      }
   } 

   if (PyBool_Check(r)) {
      Py_INCREF(r);
   }
   return r;

}
#endif


#ifdef USE_GUILE
void copy_residue_range_from_ncs_master_to_chains_scm(int imol, const char *master_chain_id, 
						      int residue_range_start, int residue_range_end, 
						      SCM chain_id_list_in) {

   if (is_valid_model_molecule(imol)) {
      std::string c(master_chain_id);
      std::vector<std::string> chain_id_list = generic_list_to_string_vector_internal(chain_id_list_in);
      graphics_info_t::molecules[imol].copy_residue_range_from_ncs_master_to_chains(c,
										    residue_range_start,
										    residue_range_end,
										    chain_id_list);
      graphics_draw();
   }
}

void copy_from_ncs_master_to_chains_scm(int imol, const char *master_chain_id, 
                                        SCM chain_id_list_in) {

   if (is_valid_model_molecule(imol)) {
      std::string c(master_chain_id);
      std::vector<std::string> chain_id_list = generic_list_to_string_vector_internal(chain_id_list_in);
      graphics_info_t::molecules[imol].copy_from_ncs_master_to_chains(c,
										    chain_id_list);
      graphics_draw();
   }
} 
#endif 
#ifdef USE_PYTHON
void copy_residue_range_from_ncs_master_to_chains_py(int imol, const char *master_chain_id, 
						     int residue_range_start, int residue_range_end,
						     PyObject *chain_id_list_in) {
   
   if (is_valid_model_molecule(imol)) {
      std::string c(master_chain_id);
      std::vector<std::string> chain_id_list = generic_list_to_string_vector_internal_py(chain_id_list_in);
      graphics_info_t::molecules[imol].copy_residue_range_from_ncs_master_to_chains(c,
										    residue_range_start,
										    residue_range_end,
										    chain_id_list);
      graphics_draw();
   }
}

void copy_from_ncs_master_to_chains_py(int imol, const char *master_chain_id, 
                                       PyObject *chain_id_list_in) {
   
   if (is_valid_model_molecule(imol)) {
      std::string c(master_chain_id);
      std::vector<std::string> chain_id_list = generic_list_to_string_vector_internal_py(chain_id_list_in);
      graphics_info_t::molecules[imol].copy_from_ncs_master_to_chains(c,
										    chain_id_list);
      graphics_draw();
   }
} 
#endif 
