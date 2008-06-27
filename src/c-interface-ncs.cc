/* src/c-interface-ncs.cc
 * 
 * Copyright 2004, 2005 The University of York
 * Copyright 2007, 2008 The University of Oxford
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

#include "mmdb_manager.h"
#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb-crystal.h"

#include "graphics-info.h"
#include "atom-utils.h" // asc_to_graphics

#include "c-interface.h"
#include "cc-interface.hh"

#ifdef USE_GUILE
#include <guile/gh.h>
#if (SCM_MAJOR_VERSION > 1) || (SCM_MINOR_VERSION > 7)
// no fix up needed 
#else    
#define scm_to_locale_string SCM_STRING_CHARS
#define scm_to_int     gh_scm2int
#define scm_to_double  gh_scm2double
#define scm_is_true    gh_scm2bool
#define scm_cdr        SCM_CDR
#endif // SCM version
#endif // USE_GUILE

// Including python needs to come after graphics-info.h, because
// something in Python.h (2.4 - chihiro) is redefining FF1 (in
// ssm_superpose.h) to be 0x00004000 (Grrr).
//
#ifdef USE_PYTHON
#include "Python.h"
#endif // USE_PYTHON

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

      graphics_info_t::molecules[imol].add_strict_ncs_matrix(chain_id,
							     tch,
							     cm44);

      graphics_draw();
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

   std::cout << "DEBUG ncs_control_change_ncs_master_to_chain imol: " << imol
	     << " and ichain: " << ichain << std::endl;
   if (is_valid_model_molecule(imol)) {
      std::vector<std::string> chain_ids =
	 coot::util::chains_in_molecule(graphics_info_t::molecules[imol].atom_sel.mol);
      if (ichain < int(chain_ids.size()))
	 graphics_info_t::molecules[imol].set_ncs_master_chain(chain_ids[ichain]);
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
	    std::cout << "%%%%%%%%% calling fill_ghost_info from c-interfac.cc make_ncs_ghosts_maybe"
		      << std::endl;
	    graphics_info_t::molecules[imol].fill_ghost_info(1, graphics_info_t::ncs_homology_level);
	    // gtk_widget_destroy(w);
	 }
      }
   } 
}

// 
int make_dynamically_transformed_ncs_maps(int imol_model, int imol_map) {

   int nmaps=0;

   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {
	 float homology_lev = graphics_info_t::ncs_homology_level;

	 // nmaps = graphics_info_t::molecules[imol_model].make_dynamically_transformed_maps(imol_map, graphics_info_t::ncs_maps_do_average_flag, graphics_info_t::ncs_homology_level);

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
	    int imol = graphics_info_t::create_molecule();
	    g.molecules[imol].install_ghost_map(g.molecules[imol_map].xmap_list[0],
						local_ncs_ghosts[ighost].name,
						local_ncs_ghosts[ighost],
						g.molecules[imol_map].is_difference_map_p(),
						g.swap_difference_map_colours,
						g.molecules[imol_map].map_sigma());
	    nmaps++;
	 }

	 if (graphics_info_t::ncs_maps_do_average_flag) {
	    std::vector<std::pair<clipper::Xmap<float>, std::string> > xmaps  = 
	       g.molecules[imol_model].ncs_averaged_maps(g.molecules[imol_map].xmap_list[0], homology_lev);
	    std::cout << "INFO:: made " << xmaps.size() << " averaged map(s)" << std::endl;
	    for (unsigned int i=0; i<xmaps.size(); i++) { 
	       std::string name;
	       name += xmaps[i].second;
	       int imol = g.create_molecule();
	       g.molecules[imol].new_map(xmaps[i].first, name);
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
      if (graphics_info_t::molecules[imol].has_map()) {
	 std::string map_str = "ncs_maps_maps_radiobutton_";
	 map_str += graphics_info_t::int_to_string(imol);
	 map_button = lookup_widget(dialog, map_str.c_str());
	 if (map_button) {
	    if (GTK_TOGGLE_BUTTON(map_button)->active) { 
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
	       if (GTK_TOGGLE_BUTTON(coords_button)->active) {
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

   if (!found_active_button_for_coords) {
      std::cout << "You need to define a set of coordinates for NCS maping\n";
   } else {
      if (!found_active_button_for_map) {
	 std::cout << "You need to define a map for NCS maping\n";
      } else {
	 make_dynamically_transformed_ncs_maps(imol_coords, imol_map);
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
		  PyObject *thisr = PyList_GetSlice(py_residue(this_res), 2, 4);
		  
		  PyObject *masta = PyList_GetSlice(py_residue(target_res), 2, 4);

		  // res_l list seems only to have one element in Paul's scm code
		  // currently?! Correct?!
		  PyObject *res_l = PyList_New(0);
		  PyList_Append(res_l, thisr);
		  PyList_Append(res_l, masta);
		  PyList_Append(res_l, PyFloat_FromDouble(cd.residue_info[iresinf].mean_diff));
		  PyList_Append(l_residue_data, res_l);
	       }
	       PyList_Append(r, PyString_FromString(cd.peer_chain_id.c_str()));
	       PyList_Append(r, PyString_FromString(diffs.target_chain_id.c_str()));
	       PyList_Append(r, l_residue_data); 
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
       PyObject *ghost_py = PyList_New(0);
       PyObject *display_it_flag_py = Py_False;
       if (ncs_ghosts[ighost].display_it_flag)
	 display_it_flag_py = Py_True;
       PyObject *rtop_py = Py_False;
       if (graphics_info_t::molecules[imol].ncs_ghosts_have_rtops_p()) {
	 rtop_py = rtop_to_python(ncs_ghosts[ighost].rtop);
       }
       PyObject *target_chain_id_py = PyString_FromString(ncs_ghosts[ighost].target_chain_id.c_str());
       PyObject *chain_id_py = PyString_FromString(ncs_ghosts[ighost].chain_id.c_str());
       PyObject *name_py = PyString_FromString(ncs_ghosts[ighost].name.c_str());

       PyList_Append(ghost_py, name_py);
       PyList_Append(ghost_py, chain_id_py);
       PyList_Append(ghost_py, target_chain_id_py);
       PyList_Append(ghost_py, rtop_py);
       PyList_Append(ghost_py, display_it_flag_py);
       PyList_Append(r, ghost_py);
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
#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
   graphics_info_t g;
   g.ncs_diffs_from_mol(imol);
#else    
   printf("not compiled with HAVE_GTK_CANVAS/GNOME_CANVAS - remake\n"); 
#endif /* HAVE_GTK_CANVAS */

}
