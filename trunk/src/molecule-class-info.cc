
/* src/molecule-class-info.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2007 by Paul Emsley
 * Copyright 2007, 2008, 2009 by The University of Oxford
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

#ifndef _MSC_VER
#include <unistd.h>
#else
#include <windows.h>
#include <direct.h>
#endif

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>

#include "mmdb_manager.h"
#include "mmdb_tables.h"
#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb-crystal.h"
#include "gtk-manual.hh"

#include "clipper/ccp4/ccp4_mtz_io.h"
#include "clipper/ccp4/ccp4_map_io.h"
#include "clipper/core/xmap.h"
#include "clipper/cns/cns_map_io.h"
#include "clipper/core/hkl_compute.h"
#include "clipper/core/map_utils.h" // Map_stats
#include "clipper/core/resol_basisfn.h"
#include "clipper/core/resol_targetfn.h"
#include "clipper/mmdb/clipper_mmdb.h"
#include "clipper/clipper-phs.h"
#include "clipper/contrib/sfcalc_obs.h"
#include "clipper/contrib/sfscale.h"
#include "clipper/contrib/sfweight.h"

#include "coot-sysdep.h"

// For stat, mkdir:
#include <sys/types.h>
#include <sys/stat.h>

#include "Bond_lines.h"

#include "gl-matrix.h"
#include "graphics-info.h"

#include "Bond_lines_ext.h"  
#include "globjects.h" // for set_bond_colour(), r_50

#include "coot-coord-utils.hh"
#include "coot-utils.hh"

#include <GL/glut.h> // needed (only?) for wirecube

#ifndef CLIPPER_MAP_INTERP
#include "clipper/core/map_interp.h"
#endif 

#include "ligand.hh"

// for debugging
#include "c-interface.h"





// Return the molecule number of the molecule that we just filled.
// Return -1 if there was a failure.
// 
int
molecule_class_info_t::handle_read_draw_molecule(int imol_no_in,
						 std::string filename,
						 std::string cwd,
						 short int reset_rotation_centre,
						 short int is_undo_or_redo,
						 float bond_width_in) {

   //
   graphics_info_t g;
   imol_no = imol_no_in;

   if (! is_undo_or_redo) 
      bond_width = bond_width_in;

   // std::cout << "DEBUG:: ---- imol_no is now " << imol_no << std::endl;

   // need to check that filename exists and is a file.
   //
   struct stat s;
   int status = stat(filename.c_str(), &s);

   // stat check the link targets not the link itself, lstat stats the
   // link itself.
   // 
   if (status != 0 || !S_ISREG (s.st_mode)) {
      std::cout << "Error reading " << filename << std::endl;
      if (S_ISDIR(s.st_mode)) {
	 std::cout << filename << " is a directory." << endl;
      }
      return -1; // which is status in an error
   }

				    
   // Read in pdb, [shelx files use the read_shelx_ins_file method]
   //

   atom_sel = get_atom_selection(filename);
   
   if (atom_sel.read_success == 1) {


      // LINK info:
      int n_models = atom_sel.mol->GetNumberOfModels();
      std::cout << "INFO:: Found " << n_models << " models\n";
      for (int imod=1; imod<=n_models; imod++) {
	 CModel *model_p = atom_sel.mol->GetModel(imod);
	 if (model_p) { 
	    int n_links = model_p->GetNumberOfLinks();
	    std::cout << "   Model "  << imod << " had " << n_links
		      << " links\n";} 
      }
      

      //
      // and move mol_class_info to indexed molecule[n_molecules];
      // note indexing difficulties/considerations.
      
      // save it in the static molecule_class_info_t
      //
      
      // CMMDBCryst *cryst_p =  (atom_sel.mol)->get_cell_p();
      
      mat44 my_matt;

      // 
      // 
      int err = atom_sel.mol->GetTMatrix(my_matt, 0, 0, 0, 0);
      if (err != SYMOP_Ok) {
	 cout << "!! Warning:: No symmetry available for this molecule"
	      << endl;
      } else { 
	 cout << "Symmetry available for this molecule" << endl;
      }
      
      // initialize some things.
      //
//       std::cout << "initialize_coordinate_things_on_read_molecule_internal for mol no "
// 		<< imol_no << std::endl;
      initialize_coordinate_things_on_read_molecule_internal(filename, is_undo_or_redo);

      set_have_unit_cell_flag_maybe();

      if (molecule_is_all_c_alphas()) {
	 ca_representation();
      } else {

	 if (! is_undo_or_redo) { 
	    // 	 std::cout << "DEBUG:: filling ghost info in
	    // 	 handle_read_draw_molecule" << std::endl;
	    short int do_rtops_flag = 0;
	    // 0.7 is not used (I think) if do_rtops_flag is 0.
	    // hack to fix Mac bug/strangeness

	    // It only makes sense to do NCS chain searching on
	    // crystallographic models.  Which generally only have one
	    // model per manager.  This may change in future...
	    int nmodels = atom_sel.mol->GetNumberOfModels();
	    if (nmodels == 1) { 
	       fill_ghost_info(do_rtops_flag, 0.7); // returns nghosts
	       // std::cout << "INFO:: found " << nghosts << " ghosts\n";
	    }
	 } else {
	    update_mols_in_additional_representations(); //uses atom_sel.mol
	 } 
	 // Generate bonds and save them in the graphical_bonds_container
	 // which has static data members.
	 //
	 if (bonds_box_type == coot::UNSET_TYPE)
	    bonds_box_type = coot::NORMAL_BONDS;
	 make_bonds_type_checked();
      }

      drawit = 1;
      if (g.show_symmetry == 1) {
	 if (show_symmetry) {  // internal
	    update_symmetry();
	 }
      }

      // Now, we have no map assocaited with this molecule, 
      // so we set the draw_vects to zero.
      // 
      // However, in future, we will have maps associated with
      // coordinates, so we should put a test here first before
      // we:
      //

      // debug();
      
      //
      if (g.recentre_on_read_pdb || imol_no_in == 0)  // n_molecules is updated
						      // in c-interface.cc
	 if (reset_rotation_centre) 
	    g.setRotationCentre(::centre_of_molecule(atom_sel)); 

      // update the maps so that they appear around the new centre. 
      // 
      if (reset_rotation_centre) 
	 for (int ii=0; ii<g.n_molecules(); ii++) {
	    g.molecules[ii].update_map(); 
	 }

      // save state strings
      save_state_command_strings_.push_back("handle-read-draw-molecule");
      std::string f1 = coot::util::intelligent_debackslash(filename);
      std::string f2 = coot::util::relativise_file_name(f1, cwd);
      save_state_command_strings_.push_back(single_quote(f2));

      return 1;

   } else {
      std::cout << "There was a coordinates read error\n";
      return -1;
   }
}

// Cleaner interface to molecule's attributes:
std::pair<bool, clipper::Spacegroup>
molecule_class_info_t::space_group() const {

   clipper::Spacegroup sg;
   std::pair<bool, clipper::Spacegroup> p(0, sg);

   if (has_model()) {
      // I want just the symmetry
      
      try { // for now
	 std::pair<clipper::Cell, clipper::Spacegroup> cell_sg =
	    coot::util::get_cell_symm(atom_sel.mol);
	 if (!cell_sg.second.is_null()) { 
	    p.first = 1;
	    p.second = cell_sg.second;
	 }
      }
      catch (std::runtime_error rte) {
	 std::cout << "ERROR:: " << rte.what() << std::endl;
      }
   }
   return p;
}

std::pair<bool, clipper::Cell>
molecule_class_info_t::cell() const {

   clipper::Cell cell;
   std::pair<bool, clipper::Cell> p(0, cell);
   if (has_map()) {
      p = std::pair<bool, clipper::Cell> (1, xmap_list[0].cell());
   } 

   if (has_model()) { 
      realtype a[6];
      realtype vol;
      int orthcode;
      atom_sel.mol->GetCell(a[0], a[1], a[2], a[3], a[4], a[5], vol, orthcode);
      clipper::Cell_descr cdr(a[0], a[1], a[2],
			      clipper::Util::d2rad(a[3]),
			      clipper::Util::d2rad(a[4]),
			      clipper::Util::d2rad(a[5]));
      p.first = 1;
      p.second = clipper::Cell(cdr);
   } 
   return p;
}


coot::Cartesian
molecule_class_info_t::centre_of_molecule() const {

   double xs=0, ys=0, zs=0;
   if (atom_sel.n_selected_atoms > 0) {
      for (int i=0; i<atom_sel.n_selected_atoms; i++) {
	 xs += atom_sel.atom_selection[i]->x;
	 ys += atom_sel.atom_selection[i]->y;
	 zs += atom_sel.atom_selection[i]->z;
      }
      xs /= double(atom_sel.n_selected_atoms);
      ys /= double(atom_sel.n_selected_atoms);
      zs /= double(atom_sel.n_selected_atoms);
   }
   return coot::Cartesian(xs, ys, zs);
}

float
molecule_class_info_t::size_of_molecule() const {
   
   double d2_sum = 0;
   float r = 0;
   coot::Cartesian centre = centre_of_molecule();
   if (atom_sel.n_selected_atoms > 0) {
      for (int i=0; i<atom_sel.n_selected_atoms; i++) {
	 double d =
	    (atom_sel.atom_selection[i]->x - centre.x()) *
	    (atom_sel.atom_selection[i]->x - centre.x()) + 
	    (atom_sel.atom_selection[i]->y - centre.y()) *
	    (atom_sel.atom_selection[i]->y - centre.y()) + 
	    (atom_sel.atom_selection[i]->z - centre.z()) *
	    (atom_sel.atom_selection[i]->z - centre.z());
	 // std::cout << i << " adding in d = " << d << std::endl;
	 d2_sum += d;
      }
      double msd = d2_sum/double(atom_sel.n_selected_atoms);
      r = sqrt(msd);
//       std::cout << " for imol = " << imol_no << " " 
// 	       << r << " = sqrt(" << msd << ") " << " = sqrt("
// 	       << d2_sum << "/" << atom_sel.n_selected_atoms << ")" << std::endl;
   }
   return r;
}


std::string
molecule_class_info_t::show_spacegroup() const { 

   std::string s("No spacegroup");
   
   if (has_model()) {
      char *st = atom_sel.mol->GetSpaceGroup();
      if (st)
	 s = st;
   }
   
   if (has_map()) 
      s = xmap_list[0].spacegroup().symbol_hm();

   return s;
}


coot::at_dist_info_t
molecule_class_info_t::closest_atom(const coot::Cartesian &pt) const {

   coot::at_dist_info_t at_info = closest_atom(pt, 1);
   if (at_info.atom)
      return at_info;
   else
      return closest_atom(pt, 0);

}

coot::at_dist_info_t
molecule_class_info_t::closest_atom(const coot::Cartesian &pt, bool ca_check_flag) const {
   return closest_atom(pt, ca_check_flag, "", 0); // don't use chain-id
}

coot::at_dist_info_t
molecule_class_info_t::closest_atom(const coot::Cartesian &pt, bool ca_check_flag,
				    const std::string &chain_id,
				    bool use_this_chain_id) const {

   coot::at_dist_info_t at_info(0,0,0);
   CAtom *at_best = 0;
   float dist_best = 99999999999.9;

   for (int iat=0; iat<atom_sel.n_selected_atoms; iat++) {
      CAtom *at = atom_sel.atom_selection[iat];
      std::string chain_id_from_at(at->GetChainID());
      if ((chain_id_from_at == chain_id) || !use_this_chain_id) { 
	 float d2 = (at->x - pt.x()) * (at->x - pt.x());
	 d2 += (at->y - pt.y()) * (at->y - pt.y());
	 d2 += (at->z - pt.z()) * (at->z - pt.z());
	 if (d2 < dist_best) {
	    dist_best = d2;
	    at_best = at;
	    // Now, does this at belong to a residue that has a CA?  If
	    // it does, reset at_best to be the CA of the residue, but
	    // keep dist_best as it was, of course.
	    if (ca_check_flag == 1) {
	       CResidue *res = at->residue;
	       int natoms;
	       PPCAtom residue_atoms;
	       res->GetAtomTable(residue_atoms, natoms);
	       for (int iatom=0; iatom<natoms; iatom++) {
		  if (! residue_atoms[iatom]->isTer()) { 
		     if (! strcmp(residue_atoms[iatom]->name, " CA ")) {
			if (! strcmp(residue_atoms[iatom]->altLoc, at->altLoc)) {
			   at_best = residue_atoms[iatom];
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   if (at_best) { 
      at_info.dist = sqrt(dist_best);
      at_info.atom = at_best;
      at_info.imol = imol_no;
   }
   return at_info;
}




std::string
molecule_class_info_t::single_quote(const std::string &s) const {
   std::string r("\"");
   r += s;
   r += "\"";
   return r;
}

void
molecule_class_info_t::install_model(int imol_no_in,
				     atom_selection_container_t asc,
				     const std::string &name, 
				     short int display_in_display_control_widget_status,
				     bool is_from_shelx_ins) {

   imol_no = imol_no_in;
   graphics_info_t g;
   bond_width = g.default_bond_width; // bleugh, perhaps this should
				      // be a passed parameter?
   is_from_shelx_ins_flag = is_from_shelx_ins;

   atom_sel = asc;

   CMMDBCryst *cryst_p =  (atom_sel.mol)->get_cell_p();
   mat44 my_matt;
   
   int err = cryst_p->GetTMatrix(my_matt, 0, 0, 0, 0);
   if (err != 0) {
      std::cout << "!! Warning:: No symmetry available for this molecule"
		<< std::endl;
   } else { 
      std::cout << "Symmetry available for this molecule" << std::endl;
   }
   set_have_unit_cell_flag_maybe();
   
   makebonds();
   if (g.show_symmetry == 1)
      if (show_symmetry) 
	 update_symmetry();

   have_unsaved_changes_flag = 1; 

   short int is_undo_or_redo = 0;

   if (display_in_display_control_widget_status == 0) {
      // treat as undo then (e.g. "terminal residue") made by add_cb_to_terminal_res().
      is_undo_or_redo = 1;
   } else {
      pickable_atom_selection = 1;
   }
   initialize_coordinate_things_on_read_molecule_internal(name, is_undo_or_redo);
}


void
molecule_class_info_t::label_atoms(int brief_atom_labels_flag) {

   if (drawit) {
      
      if (has_model()) { 

	 // keep a list of atoms that have labels (either in graphics_info
	 // or mol_class_info) and loop over them calling label_atom(i)
	 // which labels the i'th atom of the atom selection in
	 // mol_class_info.
	 //
	 //
	 //
	 int n_atoms_to_label = labelled_atom_index_list.size();

	 // also remove labels from atom indexes list of over the end.
	 for (int ii=0; ii<n_atoms_to_label ; ii++)
	    label_atom(labelled_atom_index_list[ii], brief_atom_labels_flag);
      
	 n_atoms_to_label = labelled_symm_atom_index_list.size();

	 for (int ii=0; ii<n_atoms_to_label ; ii++) {

	       // symm_trans_t st = g.labelled_symm_atom_symm_trans(ii);

	       // cout << "symm label for atom: " <<  ii << " ppc index: "
	       //     << g.labelled_symm_atom(ii) << " " << st << endl;

	       //label_symm_atom(g.labelled_symm_atom(ii), st);
	       test_label_symm_atom(ii);
	 }
      }
   }
}

void
molecule_class_info_t::trim_atom_label_table() {

   int new_max_atom_index = atom_sel.n_selected_atoms;

   // Man, this is a clumsy way of removing specific ints from a vector<int>.
   // It's needed because as we do an erase, the labelled_atom_index_list changes
   // so that we were double erasing an element:
   // erasing *it: 141 limit: 135
   // erasing *it: 141 limit: 135
   //    -> crash.
   // So do one at a time.
   //
   // Perhaps it would be better to use a queue or construct a new
   // vector<int> and reassign labelled_atom_index_list at the end.
   //
   // Anyway, this seems to work as it is, so I'll leave it for now.
   // 
   bool have_over_end_labels = 1;
   while (have_over_end_labels) {
      std::vector<int>::iterator it;
      have_over_end_labels = 0;
      for (it = labelled_atom_index_list.begin();
	   it != labelled_atom_index_list.end();
	   it++) {
	 if (*it >= new_max_atom_index) {
	    labelled_atom_index_list.erase(it);
	    have_over_end_labels = 1;
	    break;
	 }
      }
   }

   // and now for symmetry index
   //
   have_over_end_labels = 1;
   while (have_over_end_labels) {
      std::vector<int>::iterator it;
      have_over_end_labels = 0;
      for (it = labelled_symm_atom_index_list.begin();
	   it != labelled_symm_atom_index_list.end();
	   it++) {
	 if (*it >= new_max_atom_index) {
	    labelled_symm_atom_index_list.erase(it);
	    have_over_end_labels = 1;
	    break;
	 }
      }
   }
}


void
molecule_class_info_t::anisotropic_atoms() {
   
   int c; // atom colour

   if (has_model()) { 
      graphics_info_t g;
      if (drawit) {
	 if (g.show_aniso_atoms_flag == 1 ) {
	    glPushMatrix();

	    float rx = g.X();
	    float ry = g.Y();
	    float rz = g.Z();

	    float x1, y1, z1;
	    float x_diff, y_diff, z_diff;
	    float d2, mc_r2 = g.show_aniso_atoms_radius*g.show_aniso_atoms_radius;
	    float rad_50, r;

	    for (int i=0; i<atom_sel.n_selected_atoms; i++) {
      
	       // put a wiresphere at the atom positions
	 
	       if (atom_sel.atom_selection[i]->u11 > 0) {
	 
		  glLineWidth(1.0);
		  glPushMatrix();
	 
		  x1 = atom_sel.atom_selection[i]->x;
		  y1 = atom_sel.atom_selection[i]->y;
		  z1 = atom_sel.atom_selection[i]->z;

		  x_diff = x1 - rx;
		  y_diff = y1 - ry;
		  z_diff = z1 - rz;

		  d2 = x_diff*x_diff + y_diff*y_diff + z_diff*z_diff;

		  // are we either inside the distance or there is no distance set?
		  //
		  if ( (d2 <= mc_r2) || (g.show_aniso_atoms_radius_flag == 0) ) { 
	    
		     c = atom_colour(atom_sel.atom_selection[i]->element);
		     set_bond_colour_by_mol_no(c);
	       
		     GL_matrix mat(atom_sel.atom_selection[i]->u11, 
				   atom_sel.atom_selection[i]->u12,
				   atom_sel.atom_selection[i]->u13,
				   atom_sel.atom_selection[i]->u12, 
				   atom_sel.atom_selection[i]->u22,
				   atom_sel.atom_selection[i]->u23,
				   atom_sel.atom_selection[i]->u13, 
				   atom_sel.atom_selection[i]->u23,
				   atom_sel.atom_selection[i]->u33);

		     glTranslatef(x1, y1, z1);
		     // glMultMatrixf(mat.get());
		     // std::cout << "Atom Us " << std::endl;
		     // mat.print_matrix();
		     // std::cout << "Choleskied: " << std::endl;
		     // mat.cholesky().print_matrix();
		     std::pair<bool,GL_matrix> chol_pair = mat.cholesky();
		     if (chol_pair.first) { 
			glMultMatrixf(chol_pair.second.get());
			rad_50 = r_50(atom_sel.atom_selection[i]->element);
			r = rad_50_and_prob_to_radius(rad_50,
						      g.show_aniso_atoms_probability);
			// note: g.show_aniso_atoms_probability is in the range
			// 0.0 -> 100.0
			glutWireSphere(r, 10, 10);
		     } else {
			std::cout << "Bad Anistropic Us for " << atom_sel.atom_selection[i]
				  << std::endl;
		     }
		  } 
		  glPopMatrix();
	       }
	    }
	    glPopMatrix();
	 }
      }
   }
}

// fix the name to something involving rotation perhaps?
// 
// not const because bond_colour_internal is set.
void 
molecule_class_info_t::set_bond_colour_by_mol_no(int i) {

   if (bonds_rotate_colour_map_flag == 0) {
      set_bond_colour(i);
   } else {
      std::vector<float> rgb(3);
      // float rotation_size = float(imol_no + 1) * bonds_colour_map_rotation/360.0;
      float rotation_size = bonds_colour_map_rotation/360.0;

//       std::cout << " ::::::::::: in set_bond_colour bonds_colour_map_rotation is "
// 		<< bonds_colour_map_rotation << " for imol " << imol_no << std::endl;

      // rotation_size typically then: 2*32/360 = 0.178
      
      while (rotation_size > 1.0) { // no more black bonds?
	 rotation_size -= 1.0;
      } 
      switch (i) {
      case YELLOW_BOND: 
	  rgb[0] = 0.8; rgb[1] =  0.8; rgb[2] =  0.3;
	 break;
      case BLUE_BOND: 
	  rgb[0] = 0.5; rgb[1] =  0.5; rgb[2] =  1.0;
	 break;
      case RED_BOND: 
	  rgb[0] = 1.0; rgb[1] =  0.3; rgb[2] =  0.3;
	 break;
      case GREEN_BOND:
	 rgb[0] = 0.1; rgb[1] =  0.99; rgb[2] =  0.1;
	 break;
      case GREY_BOND: 
	 rgb[0] = 0.7; rgb[1] =  0.7; rgb[2] =  0.7;
	 break;
      case HYDROGEN_GREY_BOND: 
	 rgb[0] = 0.6; rgb[1] =  0.6; rgb[2] =  0.6;
	 break;
// replaced in mmdb-extras.h
//       case white:   
// 	 rgb[0] = 0.99; rgb[1] =  0.99; rgb[2] = 0.99;
// 	 break;
      case MAGENTA_BOND:
	 rgb[0] = 0.99; rgb[1] =  0.2; rgb[2] = 0.99;
	 break;
      case ORANGE_BOND:
	 rgb[0] = 0.89; rgb[1] =  0.89; rgb[2] = 0.1;
	 break;
      case CYAN_BOND:
	 rgb[0] = 0.1; rgb[1] =  0.89; rgb[2] = 0.89;
	 break;
      
      default:
	 rgb[0] = 0.8; rgb[1] =  0.2; rgb[2] =  0.2;
	 rgb = rotate_rgb(rgb, float(i*26.0/360.0));

      }

      // "correct" for the +1 added in the calculation of the rotation
      // size.
      // 21. is the default colour map rotation
      rgb = rotate_rgb(rgb, float(1.0 - 21.0/360.0));

      if (graphics_info_t::rotate_colour_map_on_read_pdb_c_only_flag) {
	 if (i == YELLOW_BOND) { 
	    std::vector<float> rgb_new = rotate_rgb(rgb, rotation_size);
	    bond_colour_internal = rgb_new;
	    glColor3f(rgb_new[0],rgb_new[1], rgb_new[2]);
	 } else {
	    bond_colour_internal = rgb;
	    glColor3f(rgb[0],rgb[1], rgb[2]);
	 }
      } else {
//  	 std::cout << "DEBUG: rotating coordinates colour map by "
//  		   << rotation_size * 360.0 << " degrees " << std::endl;
	 std::vector<float> rgb_new = rotate_rgb(rgb, rotation_size);
	 bond_colour_internal = rgb_new;
	 if (graphics_info_t::use_graphics_interface_flag) { 
	    glColor3f(rgb_new[0], rgb_new[1], rgb_new[2]);
	 }
      }
   }
}


// aka rainbow - or maybe b factor, occupancy
void
molecule_class_info_t::set_bond_colour_by_colour_wheel_position(int i, int bond_type) {

   float max_colour = 30;
   std::vector<float> rgb(3);
   rgb[0] = 0.2; rgb[1] =  0.2; rgb[2] =  0.8; // blue
   
   // 30 is the size of rainbow colours, 0 -> 1.0 is the range of rainbow colours
   
   float rotation_size = 1.0 - float(i) * 0.7/max_colour;
   rgb = rotate_rgb(rgb, rotation_size);
   bond_colour_internal = rgb;
   glColor3f(rgb[0], rgb[1], rgb[2]);
}



// We find a box (symm: 2 trans: 0 0 0), but we don't find any atoms
// in it, but we still get though to bonds,
// bonds.make_graphical_symmetry_bonds
// 
// But there are no bonds generated.
// 
void
molecule_class_info_t::update_symmetry() {

   graphics_info_t g;

   // a bit of a hack...
   int shift_search_size = graphics_info_t::symmetry_shift_search_size;
      
   if ((graphics_info_t::show_symmetry == 1) &&
       (show_symmetry == 1)) {

      // don't do stuff until we have read in a molecule.
      //
      if (drawit == 1) {
      
	 molecule_extents_t extents(atom_sel, g.symmetry_search_radius);
	 graphics_info_t g;
	 coot::Cartesian point = g.RotationCentre();

      
	 // cout << "extents " << extents << endl;
	 // cout << "point:  " << point << endl;
	 std::vector<std::pair<symm_trans_t, Cell_Translation> > symm_trans_boxes =
	    extents.which_boxes(point, atom_sel, shift_search_size);

//   	 std::cout << "DEBUG:: symm_trans_boxes.size() is " 
//   		   << symm_trans_boxes.size() << std::endl;
//   	 std::cout << "Here are the symms we should check:" << std::endl;
//    	 for(int ii=0; ii<symm_trans_boxes.size(); ii++)
//   	    std::cout << ii << " " << symm_trans_boxes[ii].first << " "
//  		      << symm_trans_boxes[ii].second << std::endl;
	 

	 if (symm_trans_boxes.size() > 0) {

	    // when bonds goes out of scope (i.e. immediate after
	    // this) then the class data member vector "bonds" of the
	    // Bond_lines_container gets given back.
	    //
	    // It is with the "new"ly allocated graphical_symmetry_bonds
	    // that we need to concern ourselves.
	    //
	    Bond_lines_container bonds;

	    //
	    // delete the old symmetry_bonds_box
	    //
	    // symmetry_bonds_box.clear_up();
	    clear_up_all_symmetry();
 	    symmetry_bonds_box.clear();

// 	    for (unsigned int ibox=0; ibox<symm_trans_boxes.size(); ibox++)
// 	       std::cout << "box " << ibox << "/" << symm_trans_boxes.size()
// 			 << " " << symm_trans_boxes[ibox] << "\n";

	    symmetry_bonds_box = 
	       bonds.addSymmetry_vector_symms(atom_sel,
					      point,
					      graphics_info_t::symmetry_search_radius,
					      symm_trans_boxes,
					      symmetry_as_calphas,
					      symmetry_whole_chain_flag,
					      draw_hydrogens_flag);
	    
	 } else {
	    Bond_lines_container bonds(NO_SYMMETRY_BONDS);
	 }

	 if (show_strict_ncs_flag == 1) {
	    if (strict_ncs_matrices.size() > 0) {
	       update_strict_ncs_symmetry(point, extents);
	    }
	 }

      } else { 
	 // cout << "update_symmetry: no molecule yet" << endl;
      }
   }
}

//
void
molecule_class_info_t::draw_coord_unit_cell(const coot::colour_holder &cell_colour) {

   // Don't display if we have closed this molecule
   // (perhaps use (atom_sel.mol==NULL) instead?) (no).

   // (same test as has_model()):
   if (atom_sel.n_selected_atoms > 0) { 

      if (show_unit_cell_flag == 1) {

	 if (drawit) { 
	 
	    if (have_unit_cell == 1) {

	       glLineWidth(2.0);
	       glColor3f(cell_colour.red, cell_colour.green, cell_colour.blue);
	 
	       float corners[8][3] = {
		  {0,0,0}, //0 
		  {0,0,1}, //1 
		  {0,1,0}, //2 
		  {0,1,1}, //3 
		  {1,0,0}, //4 
		  {1,0,1}, //5 
		  {1,1,0}, //6 
		  {1,1,1}};//7 

	       realtype x_orth, y_orth, z_orth;
	       // rsc = real_space_corners
	       float rsc[8][3];

	       for (int ii=0; ii<8; ii++) {
	    
		  atom_sel.mol->Frac2Orth(corners[ii][0], corners[ii][1], corners[ii][2],
					  x_orth, y_orth, z_orth);
	    
		  rsc[ii][0] = x_orth;
		  rsc[ii][1] = y_orth;
		  rsc[ii][2] = z_orth;
	       }
	 
	       draw_unit_cell_internal(rsc);
	    }
	 }
      }
   }
}

//
void
molecule_class_info_t::draw_map_unit_cell(const coot::colour_holder &cell_colour) {

   if (has_map()) { 
      if (show_unit_cell_flag == 1) {
      
	 if ( max_xmaps > 0 ) {
	    if ( xmap_is_filled[0] ) {

	       if (drawit_for_map) { 
	    
		  // rsc = real_space_corners
		  float rsc[8][3];
	    
		  glLineWidth(2.0);
		  glColor3f(cell_colour.red, cell_colour.green, cell_colour.blue);
	    
		  float corners[8][3] = {
		     {0,0,0}, //0 
		     {0,0,1}, //1 
		     {0,1,0}, //2 
		     {0,1,1}, //3 
		     {1,0,0}, //4 
		     {1,0,1}, //5 
		     {1,1,0}, //6 
		     {1,1,1}};//7 
	    
		  for (int ii=0; ii<8; ii++) {
	       
		     clipper::Coord_frac c_f(corners[ii][0],corners[ii][1],corners[ii][2]);
	       
		     clipper::Coord_orth c_o = c_f.coord_orth( xmap_list[0].cell());
	       
		     rsc[ii][0] = c_o.x();
		     rsc[ii][1] = c_o.y();
		     rsc[ii][2] = c_o.z();
		  }
		  draw_unit_cell_internal(rsc); 
	       }
	    }
	 }
      }
   }
}

//
void
molecule_class_info_t::draw_unit_cell_internal(float rsc[8][3]) {

   std::vector<coot::CartesianPair> p;

   // bottom left connections
   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[0][0], rsc[0][1], rsc[0][2]),
				   coot::Cartesian(rsc[1][0], rsc[1][1], rsc[1][2])));

   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[0][0], rsc[0][1], rsc[0][2]),
				   coot::Cartesian(rsc[2][0], rsc[2][1], rsc[2][2])));
	 
   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[0][0], rsc[0][1], rsc[0][2]),
				   coot::Cartesian(rsc[4][0], rsc[4][1], rsc[4][2])));

   // top right front connections
   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[6][0], rsc[6][1], rsc[6][2]),
				   coot::Cartesian(rsc[4][0], rsc[4][1], rsc[4][2])));
	 
   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[6][0], rsc[6][1], rsc[6][2]),
				   coot::Cartesian(rsc[2][0], rsc[2][1], rsc[2][2])));
	 
   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[6][0], rsc[6][1], rsc[6][2]),
				   coot::Cartesian(rsc[7][0], rsc[7][1], rsc[7][2])));
	 
   // from 5
   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[5][0], rsc[5][1], rsc[5][2]), 
				   coot::Cartesian(rsc[7][0], rsc[7][1], rsc[7][2])));

   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[5][0], rsc[5][1], rsc[5][2]), 
				   coot::Cartesian(rsc[4][0], rsc[4][1], rsc[4][2])));

   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[5][0], rsc[5][1], rsc[5][2]), 
				   coot::Cartesian(rsc[1][0], rsc[1][1], rsc[1][2])));

   // from 3
   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[3][0], rsc[3][1], rsc[3][2]), 
				   coot::Cartesian(rsc[1][0], rsc[1][1], rsc[1][2])));

   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[3][0], rsc[3][1], rsc[3][2]),
				   coot::Cartesian(rsc[7][0], rsc[7][1], rsc[7][2])));

   p.push_back(coot::CartesianPair(coot::Cartesian(rsc[3][0], rsc[3][1], rsc[3][2]), 
				   coot::Cartesian(rsc[2][0], rsc[2][1], rsc[2][2])));

   float x1, y1, z1;
   float x2, y2, z2;
   glBegin(GL_LINES);
   for (unsigned i=0; i<p.size(); i++) {
//       glVertex3f(p[i].getStart().x(),  p[i].getStart().y(),  p[i].getStart().z());
//       glVertex3f(p[i].getFinish().x(), p[i].getFinish().y(), p[i].getFinish().z());
      coot::Cartesian diff = p[i].getFinish() - p[i].getStart();
      for (float j=0.0; j<0.999; j+=0.1) {
	 x1 = p[i].getStart().x() + (j)*diff.x();
	 y1 = p[i].getStart().y() + (j)*diff.y();
	 z1 = p[i].getStart().z() + (j)*diff.z();
	 x2 = p[i].getStart().x() + (j+0.1)*diff.x();
	 y2 = p[i].getStart().y() + (j+0.1)*diff.y();
	 z2 = p[i].getStart().z() + (j+0.1)*diff.z();
	 
	 glVertex3f(x1, y1, z1);
	 glVertex3f(x2, y2, z2);
      }
   }
   glEnd();

	   
	 // add a label
// 	 glColor3f(1.0, 1.0, 1.0);
// 	 glColor3f(1.0, 0.2, 1.0);
	 glRasterPos3f(-1.6, -1.6,-1.6);      
	 printString("0");
	 glRasterPos3f(rsc[1][0]-1, rsc[1][1], rsc[1][2]+1);
	 printString("C");
	 glRasterPos3f(rsc[2][0]+1, rsc[2][1], rsc[2][2]+1);      
	 printString("B");
	 glRasterPos3f(rsc[4][0]+1, rsc[4][1]+1, rsc[4][2]-1);      
	 printString("A");


}

// --------------------------------------------------------------------
//   Conversion functions
// --------------------------------------------------------------------
//
void
molecule_class_info_t::initialize_coordinate_things_on_read_molecule(std::string molecule_name) {
   
   // presume not an undo/redo by default:
   initialize_coordinate_things_on_read_molecule_internal(molecule_name, 0);
}

// If we are a redo/undo, then we don't want to update (add a) mol in
// display control widget
// 
// Or a non graphics_info_t::molecules[] usage of this class.
// 
void
molecule_class_info_t::initialize_coordinate_things_on_read_molecule_internal(std::string molecule_name,
									      short int is_undo_or_redo) {

   // we use xmap_is_filled[0] to see if this molecule is a map
   //
   // FIXME.  Delete these lines, max_xmaps should be used instead.
   //
   xmap_is_filled = new int[10];
   xmap_is_filled[0] = 0;

   //
   name_ = molecule_name;

   // 
   drawit = 1; // by default, display it, we change change this later, if we want.

   //
   if (! is_undo_or_redo) { 
      bonds_colour_map_rotation = (imol_no + 1) * graphics_info_t::rotate_colour_map_on_read_pdb;
      while (bonds_colour_map_rotation > 360.0)
	 bonds_colour_map_rotation -= 360.0;
      bonds_rotate_colour_map_flag = graphics_info_t::rotate_colour_map_on_read_pdb_flag;
//       std::cout << "::::::: in initialization setting bonds_colour_map_rotation "
// 		<< bonds_colour_map_rotation << " for imol no " << imol_no << std::endl;
   }

   if (! is_undo_or_redo) { 
      // std::cout << "DEBUG:: not an undo/redo!\n";
      new_coords_mol_in_display_control_widget(); // uses drawit
   }
}

void
molecule_class_info_t::set_symm_bond_colour_mol(int icol) {

   switch (icol) {
      case GREEN_BOND:
	 glColor3f (combine_colour(0.1,0),
		    combine_colour(0.8,1),
		    combine_colour(0.1,2));
	 break;
      case BLUE_BOND: 
	 glColor3f (combine_colour(0.2,0),
		    combine_colour(0.2,1),
		    combine_colour(0.8,2));
	 break;
      case RED_BOND: 
	 glColor3f (combine_colour(0.8,0),
		    combine_colour(0.1,1),
		    combine_colour(0.1,2));
	 break;
      case YELLOW_BOND: 
	 glColor3f (combine_colour(0.7,0),
		    combine_colour(0.7,1),
		    combine_colour(0.0,2));
	 break;
      
      default:
	 glColor3f (combine_colour(0.7, 0),
		    combine_colour(0.8, 1),
		    combine_colour(0.8, 2));
   }
}

void
molecule_class_info_t::set_symm_bond_colour_mol_and_symop(int icol, int isymop) {

//    std::cout << "in set_symm_bond_colour_mol_and_symop " << imol_no << " " << icol << " "
// 	     << isymop << " symmetry_rotate_colour_map_flag: "
// 	     << symmetry_rotate_colour_map_flag << "\n";
   
   if (symmetry_rotate_colour_map_flag) { 
      if (symmetry_colour_by_symop_flag) { 
	 set_symm_bond_colour_mol_rotate_colour_map(icol, isymop);
      } else { 
	 set_symm_bond_colour_mol_rotate_colour_map(icol, 0);
      }
   } else {
      set_symm_bond_colour_mol(icol);
   } 

}

void
molecule_class_info_t::set_symm_bond_colour_mol_rotate_colour_map(int icol, int isymop) {

   float step = graphics_info_t::rotate_colour_map_on_read_pdb/360.0;
   float rotation_size = float(icol + isymop) * step;
   std::vector<float> orig_colours(3);
   std::vector<float> rgb_new(3);
   std::vector<float> t_colours(3);

   switch (icol) {
   case GREEN_BOND:
      t_colours[0] = combine_colour(0.1, 0);
      t_colours[1] = combine_colour(0.8, 1);
      t_colours[2] = combine_colour(0.1, 2);
      rgb_new = rotate_rgb(t_colours, rotation_size);
      glColor3f (rgb_new[0], rgb_new[1], rgb_new[2]);
      break;
   case BLUE_BOND: 
      t_colours[0] = combine_colour(0.2, 0);
      t_colours[1] = combine_colour(0.2, 1);
      t_colours[2] = combine_colour(0.8, 2);
      rgb_new = rotate_rgb(t_colours, rotation_size);
      glColor3f (rgb_new[0], rgb_new[1], rgb_new[2]);
      break;
   case RED_BOND: 
      t_colours[0] = combine_colour(0.8, 0);
      t_colours[1] = combine_colour(0.1, 1);
      t_colours[2] = combine_colour(0.1, 2);
      rgb_new = rotate_rgb(t_colours, rotation_size);
      glColor3f (rgb_new[0], rgb_new[1], rgb_new[2]);
      break;
   case YELLOW_BOND: 
      t_colours[0] = combine_colour(0.7, 0);
      t_colours[1] = combine_colour(0.7, 1);
      t_colours[2] = combine_colour(0.0, 2);
      rgb_new = rotate_rgb(t_colours, rotation_size);
      glColor3f (rgb_new[0], rgb_new[1], rgb_new[2]);
      break;
      
   default:
      t_colours[0] = combine_colour(0.6, 0);
      t_colours[1] = combine_colour(0.7, 1);
      t_colours[2] = combine_colour(0.7, 2);
      rgb_new = rotate_rgb(t_colours, rotation_size);
      glColor3f (rgb_new[0], rgb_new[1], rgb_new[2]);
   }
}



float
molecule_class_info_t::combine_colour(float v, int col_part_index) {

   // col_part_index is 0,1,2 for red gree blue components of the colour
   double w = graphics_info_t::symmetry_colour_merge_weight;
   return w*graphics_info_t::symmetry_colour[col_part_index] + v*(1.0-w);
}

// amount is not in degrees, it is in fractions of a circle, e.g. 10/360.
// 
void
molecule_class_info_t::rotate_rgb_in_place(float *rgb, const float &amount) const {

   float *hsv;
   hsv = new float[3];
   convert_rgb_to_hsv_in_place(rgb, hsv);
   hsv[0] += amount;
   if (hsv[0] > 1.0) hsv[0] -= 1.0;
   convert_hsv_to_rgb_in_place(hsv, rgb);
   delete [] hsv;
}

// This allocated memory for xmap_is_diff_map, xmap_is_filled and
// contour_level, but *does not filll them!*
// So they need to be filled after calling this function.
void 
molecule_class_info_t::initialize_map_things_on_read_molecule(std::string molecule_name,
							      int is_diff_map,
							      short int swap_difference_map_colours) {

   // unset coordinates, this is not a set of coordinates:
   atom_sel.n_selected_atoms = 0;
   atom_sel.mol = 0;  // tested (in set_undo_molecule()) to see if this
                      // was a coordinates molecule.  So maps have to
                      // set this to NULL.

   // Map initialization:
   n_draw_vectors = 0;
   n_diff_map_draw_vectors = 0;
   xmap_is_filled = new int[10];

   // give it some memory:
   xmap_is_diff_map = new int[10];
   xmap_is_diff_map[0] = is_diff_map;  // and set the first (only) one. 

   contour_level = new float[10];

   draw_vectors = NULL;
   diff_map_draw_vectors = NULL;

   show_unit_cell_flag = 0;
   have_unit_cell      = 0; // hmmm - CHECKME.

   map_colour = new double*[10];

   map_colour[0] = new double[4];
   if (is_diff_map == 1) { 
      if (! swap_difference_map_colours) { 
	 map_colour[0][0] = 0.2; 
	 map_colour[0][1] = 0.6; 
	 map_colour[0][2] = 0.2;
      } else { 
	 map_colour[0][0] = 0.6; 
	 map_colour[0][1] = 0.2; 
	 map_colour[0][2] = 0.2; 
      } 
   } else {
      std::vector<float> orig_colours(3);
      orig_colours[0] =  0.2;
      orig_colours[1] =  0.5;
      orig_colours[2] =  0.7;
      float rotation_size = float(imol_no) * graphics_info_t::rotate_colour_map_for_map/360.0;
      // std::cout << "rotating map colour by " << rotation_size * 360.0 << std::endl;
      std::vector<float> rgb_new = rotate_rgb(orig_colours, rotation_size);
      map_colour[0][0] = rgb_new[0];
      map_colour[0][1] = rgb_new[1];
      map_colour[0][2] = rgb_new[2];
   } 
      
   // negative contour level
   // 
   map_colour[1] = new double[4];
   if (! swap_difference_map_colours) { 
      map_colour[1][0] = 0.6; 
      map_colour[1][1] = 0.2; 
      map_colour[1][2] = 0.2; 
   } else { 
      map_colour[1][0] = 0.2; 
      map_colour[1][1] = 0.6; 
      map_colour[1][2] = 0.2;
   } 
   name_ = molecule_name;

   drawit_for_map = 1; // display the map initially, by default

   // We can't call this untill xmap_is_filled[0] has been assigned,
   // and here we only make room for it.
   // 
   // update_map_in_display_control_widget();
   
}


void
molecule_class_info_t::update_mol_in_display_control_widget() const { 

   graphics_info_t g;

   // we don't want to add a display control hbox if we are simply
   // doing an undo: This is now deal with by the calling function.
   // 
//    std::cout << "update_mol_in_display_control_widget() now" << std::endl;
//    std::cout << "update_mol_in_display_control_widget() passed derefrerence imol_no_ptr: "
// 	     << *imol_no_ptr << std::endl;
   std::string dmn = name_for_display_manager();
   if (g.display_control_window()) 
      update_name_in_display_control_molecule_combo_box(g.display_control_window(), 
							dmn.c_str(), 
							imol_no);
}

void
molecule_class_info_t::new_coords_mol_in_display_control_widget() const { 

   graphics_info_t g;

   // we don't want to add a display control hbox if we are simply
   // doing an undo: This is now handled by the calling function.
   //
   bool show_add_reps_flag = 0;
   if (add_reps.size() > 0)
      show_add_reps_flag = 1;
   
   std::string dmn = name_for_display_manager();
   if (g.display_control_window()) {
      display_control_molecule_combo_box(g.display_control_window(), 
					 dmn.c_str(), 
					 imol_no, show_add_reps_flag);
      if (add_reps.size() > 0) {
	 GtkWidget *vbox = display_control_add_reps_container(g.display_control_window(), imol_no);
	 for (unsigned int iar=0; iar<add_reps.size(); iar++) {
	    std::string name = coot::util::int_to_string(iar);
	    name += " ";
	    name += add_reps[iar].info_string();
	    display_control_add_reps(vbox, imol_no, iar, add_reps[iar].show_it,
				     add_reps[iar].bonds_box_type, name);
	 } 
      } 
   }
   
}

std::string
molecule_class_info_t::name_for_display_manager() const { 

   std::string s("");
   if (graphics_info_t::show_paths_in_display_manager_flag) { 
      s = name_;
   } else {
      if (has_model()) {
#ifdef WINDOWS_MINGW
	 std::string::size_type islash = coot::util::intelligent_debackslash(name_).find_last_of("/");
#else
	 std::string::size_type islash = name_.find_last_of("/");
#endif // MINGW
	 if (islash == std::string::npos) { 
	    s = name_;
	 } else {
	    s = name_.substr(islash+1, name_.length());
	 }
      } else {
	 // This is a map, so we want to strip of xxx/ from each of
	 // the (space separated) strings.
	 // e.g.:
	 // thing/other.mtz -> other.mtz
	 // but
	 // Averaged -> Averaged

	 std::vector<std::string> v = coot::util::split_string(name_, " ");
	 for (unsigned int i=0; i<v.size(); i++) {
	    if (i > 0) 
	       s += " ";
	    std::pair<std::string, std::string> p = coot::util::split_string_on_last_slash(v[i]);
	    if (p.second == "")
	       s += v[i];
	    else 
	       s += p.second;
	 }
      } 
   }
   return s;
}

std::string
molecule_class_info_t::dotted_chopped_name() const {

   std::string ss = coot::util::int_to_string(imol_no);
   ss += " " ;
   int ilen = name_.length();
   int left_size = ilen-graphics_info_t::go_to_atom_menu_label_n_chars_max;
   if (left_size <= 0) {
      // no chop
      left_size = 0;
   } else {
      // chop
      ss += "...";
   } 
   ss += name_.substr(left_size, ilen);
   return ss;
}



void
molecule_class_info_t::add_to_labelled_atom_list(int atom_index) {

   // note initialization n_labelled_atoms is 0;
   // 
   if (is_in_labelled_list(atom_index) == 1) {
      unlabel_atom(atom_index);
   } else { 
      labelled_atom_index_list.push_back(atom_index);
   }
}

// int
// molecule_class_info_t::labelled_atom(int i) {

//    return labelled_atom_index_list[i];

// } 


// Olde Pointere Stuff.  Delete when you feel it should go.
// int
// molecule_class_info_t::max_labelled_atom() {

//    return n_labelled_atoms;
// }


// or as we would say in lisp: rember
void
molecule_class_info_t::unlabel_atom(int atom_index) {

   // 
   // Remove atom_index from the list of atoms to be labelled.
   //
   std::vector<int>::iterator it;
   for (it = labelled_atom_index_list.begin(); it != labelled_atom_index_list.end(); it++) {
      if ( *it == atom_index) {
	 labelled_atom_index_list.erase(it);
	 break;
      }
   }
}

void
molecule_class_info_t::unlabel_last_atom() {
   // remove the last atom from the list (if
   // there *are* atoms in the list, else do
   // nothing).
   unsigned int las = labelled_atom_index_list.size();
   if (las > 0) {
      int atom_index = labelled_atom_index_list[las-1];
      // std::cout << "unlabelling atom index" << atom_index << std::endl;
      unlabel_atom(atom_index);
   }
}

// or as we would say in lisp: member?
bool
molecule_class_info_t::is_in_labelled_list(int i) {

   // is the i'th atom in the list of atoms to be labelled?

   for (unsigned int ii=0; ii<labelled_atom_index_list.size(); ii++) {
      if (labelled_atom_index_list[ii] == i) {
	 return 1;
      }
   }
   return 0;
}

// ------------------------- residue exists? -------------------------

int
molecule_class_info_t::does_residue_exist_p(const std::string &chain_id,
					    int resno,
					    const std::string &inscode) const {
   int state = 0;
   if (atom_sel.n_selected_atoms > 0) {
      int n_models = atom_sel.mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) { 
	 
	 CModel *model_p = atom_sel.mol->GetModel(imod);
	 CChain *chain_p;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 if (nchains <= 0) { 
	    std::cout << "ERROR:: bad nchains in molecule " << nchains
		      << std::endl;
	 } else { 
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain_p = model_p->GetChain(ichain);
	       if (chain_p == NULL) {  
		  // This should not be necessary. It seem to be a
		  // result of mmdb corruption elsewhere - possibly
		  // DeleteChain in update_molecule_to().
		  std::cout << "NULL chain in ... " << std::endl;
	       } else {
		  PCResidue residue_p;
		  if (chain_id == chain_p->GetChainID()) { 
		     int nres = chain_p->GetNumberOfResidues();
		     for (int ires=0; ires<nres; ires++) { 
			residue_p = chain_p->GetResidue(ires);
			if (resno == residue_p->seqNum) {
			   if (inscode == residue_p->GetInsCode()) {
			      state = 1;
			      break;
			   }
			}
		     }
		  }
	       }
	    }
	 }
	 if (state)
	    break;
      }
   }
   return state;
}



// ------------------------- symmmetry atom labels -------------------------

void
molecule_class_info_t::add_atom_to_labelled_symm_atom_list(int atom_index,
							   const symm_trans_t &symm_trans,
							   const Cell_Translation &pre_shift_cell_trans) {

   if ( is_in_labelled_symm_list(atom_index) == 1 ) {
      unlabel_symm_atom(atom_index);
   } else { 
      labelled_symm_atom_index_list.push_back(atom_index);
      std::pair<symm_trans_t, Cell_Translation> p(symm_trans, pre_shift_cell_trans);
      labelled_symm_atom_symm_trans_.push_back(p);
   }
}

// no need for this now we are using vectors.
// int
// molecule_class_info_t::labelled_symm_atom(int i) {
//    return labelled_symm_atom_index_list[i];
// }



std::pair<symm_trans_t, Cell_Translation> 
molecule_class_info_t::labelled_symm_atom_symm_trans(int i) {

   return labelled_symm_atom_symm_trans_[i];
}


// old syle pointer using function
// int
// molecule_class_info_t::max_labelled_symm_atom() {

//    return n_labelled_symm_atoms;
// }

void
molecule_class_info_t::unlabel_symm_atom(int atom_index) {
   
   std::vector<int>::iterator it;
   for (it = labelled_symm_atom_index_list.begin();
	it != labelled_symm_atom_index_list.end();
	it++) {
      if ( *it == atom_index) { 
	 labelled_symm_atom_index_list.erase(it);
	 break;
      } 
   }
}

// shall we pass the symm_trans too?  Ideally we should, I think.
//
bool
molecule_class_info_t::is_in_labelled_symm_list(int i) {

   // is the i'th atom in the list of atoms to be labelled?

   for (unsigned int ii=0; ii<labelled_symm_atom_index_list.size(); ii++) {
      if (labelled_symm_atom_index_list[ii] == i) {
	 return 1;
      }
   }
   return 0;
}


int molecule_class_info_t::add_atom_label(char *chain_id, int iresno, char *atom_id) { 

   // int i = atom_index(chain_id, iresno, atom_id);
   int i = atom_spec_to_atom_index(std::string(chain_id),
				   iresno,
				   std::string(atom_id));
   if (i > 0) 
      add_to_labelled_atom_list(i);
   else
      std::cout << atom_id << "/" << iresno << "/" << chain_id
		<< " is not found in this molecule: (" <<  imol_no << ") "
		<< name_ << std::endl; 

   return i; 
}


int molecule_class_info_t::remove_atom_label(char *chain_id, int iresno, char *atom_id) {

   int i = atom_index(chain_id, iresno, atom_id);
   if (i > 0) 
      unlabel_atom(i);
   return i; 
}

 
void
molecule_class_info_t::draw_molecule(short int do_zero_occ_spots) {


   //
   //
   //  We also need a c++ object to store molecular information
   //

   if (has_model()) { 
      if (drawit == 1) {
	 if (!cootsurface) { 
	    if (do_zero_occ_spots)
	       zero_occupancy_spots();
	    display_bonds();
	    draw_fixed_atom_positions();

	    // ghosts
// 	    std::cout << "debug ghosts: " << show_ghosts_flag << " " << ncs_ghosts.size()
// 		      << std::endl;
	    if (show_ghosts_flag) {
	       if (ncs_ghosts.size() > 0) {
		  for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
		     display_ghost_bonds(ighost);
		  }
	       }
	    }
	 }
      }
   }
}


void
molecule_class_info_t::zero_occupancy_spots() const {

   if (bonds_box.n_zero_occ_spot > 0) { 

      glColor3f(0.8, 0.7, 0.7);
      glPointSize(6.5);
      glBegin(GL_POINTS); 
      for (int i=0; i<bonds_box.n_zero_occ_spot; i++) { 
	 glVertex3f(bonds_box.zero_occ_spot[i].x(),
		    bonds_box.zero_occ_spot[i].y(),
		    bonds_box.zero_occ_spot[i].z());
      }
      glEnd();
   }
}

void
molecule_class_info_t::draw_fixed_atom_positions() const {

   if (fixed_atom_positions.size() > 0) {

      glColor3f(0.6, 0.95, 0.6);
      glPointSize(10.5);
      glBegin(GL_POINTS); 
      for (unsigned int i=0; i<fixed_atom_positions.size(); i++) {
	 glVertex3f(fixed_atom_positions[i].x(),
		    fixed_atom_positions[i].y(),
		    fixed_atom_positions[i].z());
      }
      glEnd();
   }
}

void
molecule_class_info_t::display_ghost_bonds(int ighost) {


   // debug
//     int ilines = 0;
//     for (int i=0; i<ncs_ghosts[ighost].bonds_box.num_colours; i++)
//        ilines += ncs_ghosts[ighost].bonds_box.bonds_[i].num_lines;
//     std::cout << " ghost " << ighost << " has "
// 	      << ncs_ghosts[ighost].bonds_box.num_colours
// 	      << " colours and " << ilines << " lines\n";

//    std::cout << "debug ighost: " << ighost << " ncs_ghosts.size(): " << ncs_ghosts.size()
// 	     << std::endl;
   if (ighost<int(ncs_ghosts.size())) {
//       std::cout << "debug ncs_ghosts[" << ighost << "].display_it_flag "
// 		<< ncs_ghosts[ighost].display_it_flag << std::endl;
      if (ncs_ghosts[ighost].display_it_flag) {
	 glLineWidth(ghost_bond_width);
	 int c;
	 for (int i=0; i<ncs_ghosts[ighost].bonds_box.num_colours; i++) {
	    c = atom_colour(atom_sel.atom_selection[i]->element);
	    if (ncs_ghosts[ighost].bonds_box.bonds_[i].num_lines > 0) 
	       set_bond_colour_by_mol_no(ighost);
	    glBegin(GL_LINES);
	    for (int j=0; j< ncs_ghosts[ighost].bonds_box.bonds_[i].num_lines; j++) {
	       glVertex3f(ncs_ghosts[ighost].bonds_box.bonds_[i].pair_list[j].getStart().get_x(),
			  ncs_ghosts[ighost].bonds_box.bonds_[i].pair_list[j].getStart().get_y(),
			  ncs_ghosts[ighost].bonds_box.bonds_[i].pair_list[j].getStart().get_z());
	       glVertex3f(ncs_ghosts[ighost].bonds_box.bonds_[i].pair_list[j].getFinish().get_x(),
			  ncs_ghosts[ighost].bonds_box.bonds_[i].pair_list[j].getFinish().get_y(),
			  ncs_ghosts[ighost].bonds_box.bonds_[i].pair_list[j].getFinish().get_z());
	    }
	    glEnd();
	 }
      }
   }
}


// This used to use an int_grid.  Goodbye.  Clipper skeletonization is better.
// 
// void
// molecule_class_info_t::update_skeleton() {

// }


void
molecule_class_info_t::display_bonds() {

   display_bonds(bonds_box, bond_width);
   for (unsigned int i=0; i<add_reps.size(); i++) {
      if (add_reps[i].show_it) {
	 display_bonds(add_reps[i].bonds_box, add_reps[i].bond_width);
      }
   }
   display_symmetry_bonds();
}



void
molecule_class_info_t::display_bonds(const graphical_bonds_container &bonds_box,
				     float p_bond_width) {

   //

   coot::CartesianPair pair;
   Lines_list ll;
   
   glLineWidth(p_bond_width);
   //

   for (int i=0; i<bonds_box.num_colours; i++) {

      ll = bonds_box.bonds_[i];
      //cout << "j range: for i = " << i << " is "
      //	   << bonds_box.bonds_[i].num_lines << endl;

//       std::cout << "DEBUG:: bonds_index: " << i << " has " << bonds_box.bonds_[i].num_lines
// 		<< " lines, thin_flag:  " << bonds_box.bonds_[i].thin_lines_flag  << std::endl;

      if (bonds_box.bonds_[i].thin_lines_flag)
 	 glLineWidth(p_bond_width/2.0);
      else 
 	 glLineWidth(p_bond_width);

      if ( bonds_box.bonds_[i].num_lines > 512000) {
	 std::cout << "Fencepost heuristic failure bonds_box.bonds_[i].num_lines "
	      << bonds_box.bonds_[i].num_lines << std::endl;
      }
      // std::cout << "debug:: bonds_box_type " << bonds_box_type << std::endl;
      if (bonds_box_type != coot::COLOUR_BY_RAINBOW_BONDS) {
	 // if test suggested by Ezra Peisach.
	 if (bonds_box.bonds_[i].num_lines > 0)
	    set_bond_colour_by_mol_no(i); // outside inner loop
      } else {
	 set_bond_colour_by_colour_wheel_position(i, coot::COLOUR_BY_RAINBOW);
      }
      int linesdrawn = 0;

      if (1) { 
	 glBegin(GL_LINES); 
	 for (int j=0; j< bonds_box.bonds_[i].num_lines; j++) {
	    
	    // 	 if ( j > 200000) {
	    // 	    cout << "Heuristics fencepost failure j " << j << endl;
	    // 	    exit(1);
	    // 	 }
	    
	    glVertex3f(ll.pair_list[j].getStart().get_x(),
		       ll.pair_list[j].getStart().get_y(),
		       ll.pair_list[j].getStart().get_z());
	    glVertex3f(ll.pair_list[j].getFinish().get_x(),
		       ll.pair_list[j].getFinish().get_y(),
		       ll.pair_list[j].getFinish().get_z());
	    if ( (++linesdrawn & 1023) == 0) {
	       glEnd();
	       glBegin(GL_LINES);
	       linesdrawn = 0;
	    }
	 }
	 glEnd();
      }

      // quads.  This is poor currently.  Needs to be based on view
      // rotation matrix and zoom, so that we can draw the quad
      // perpendicular to the view direction (and of the right
      // thickness).
      if (0) {
	 glBegin(GL_QUADS);
	 for (int j=0; j< bonds_box.bonds_[i].num_lines; j++) {
	    glVertex3f(ll.pair_list[j].getStart().get_x()+0.05,
		       ll.pair_list[j].getStart().get_y(),
		       ll.pair_list[j].getStart().get_z());
	    glVertex3f(ll.pair_list[j].getStart().get_x()-0.05,
		       ll.pair_list[j].getStart().get_y(),
		       ll.pair_list[j].getStart().get_z());

	    glVertex3f(ll.pair_list[j].getFinish().get_x(),
		       ll.pair_list[j].getFinish().get_y()+0.05,
		       ll.pair_list[j].getFinish().get_z());
	    glVertex3f(ll.pair_list[j].getFinish().get_x(),
		       ll.pair_list[j].getFinish().get_y()-0.05,
		       ll.pair_list[j].getFinish().get_z());
	 }
	 glEnd();
      } 

      
   }

}

void
molecule_class_info_t::display_symmetry_bonds() {

   // We may come here after having done additional_representations -
   // which would change the line width.
   // 
   glLineWidth(bond_width);
   
   Lines_list ll;
   if ((show_symmetry == 1) && (graphics_info_t::show_symmetry == 1)) {
      int isymop;

      for (unsigned int isym=0; isym<symmetry_bonds_box.size(); isym++) { 
	 // isymop = isym;
	 isymop = symmetry_bonds_box[isym].second.first.isym();

	 if (symmetry_bonds_box[isym].first.symmetry_has_been_created == 1) {
	 
	    for (int icol=0; icol<symmetry_bonds_box[isym].first.num_colours; icol++) {
	 
	       set_symm_bond_colour_mol_and_symop(icol, isymop);
	       int linesdrawn = 0;
	    
	       ll = symmetry_bonds_box[isym].first.symmetry_bonds_[icol];
	 
	       glBegin(GL_LINES); 
	       for (int j=0; j< symmetry_bonds_box[isym].first.symmetry_bonds_[icol].num_lines; j++) {

		  // pair = ll.pair_list[j];
	    
		  glVertex3f(ll.pair_list[j].getStart().get_x(),
			     ll.pair_list[j].getStart().get_y(),
			     ll.pair_list[j].getStart().get_z());
		  glVertex3f(ll.pair_list[j].getFinish().get_x(),
			     ll.pair_list[j].getFinish().get_y(),
			     ll.pair_list[j].getFinish().get_z());
		  if ( (++linesdrawn & 1023) == 0) {
		     glEnd();
		     glBegin(GL_LINES);
		     linesdrawn = 0;
		  }
	       }
	       glEnd();
	    }
	 }
      }

      if (show_strict_ncs_flag == 1) {
	 // isn -> i_strict_ncs
	 for (unsigned int isn=0; isn<strict_ncs_bonds_box.size(); isn++) {

// 	    std::cout << "here 1 "
// 		      << isn << " " 
// 		      << strict_ncs_bonds_box[isn].first.symmetry_has_been_created
// 		      << " \n" ;

	    if (strict_ncs_bonds_box[isn].first.symmetry_has_been_created == 1) {
	 
	       for (int icol=0; icol<strict_ncs_bonds_box[isn].first.num_colours; icol++) {
	 
// 		  std::cout << "here 3 - isn: " << isn << " "
// 			    << "icol: " << icol << " num lines: "
// 			    << strict_ncs_bonds_box[isn].first.symmetry_bonds_[icol].num_lines
// 			    << "\n" ;
	       
		  set_symm_bond_colour_mol_and_symop(icol, isn);
		  int linesdrawn = 0;
	    
		  ll = strict_ncs_bonds_box[isn].first.symmetry_bonds_[icol];
	 
		  glBegin(GL_LINES); 
		  for (int j=0; j< strict_ncs_bonds_box[isn].first.symmetry_bonds_[icol].num_lines; j++) {

		     // pair = ll.pair_list[j];
	    
		     glVertex3f(ll.pair_list[j].getStart().get_x(),
				ll.pair_list[j].getStart().get_y(),
				ll.pair_list[j].getStart().get_z());
		     glVertex3f(ll.pair_list[j].getFinish().get_x(),
				ll.pair_list[j].getFinish().get_y(),
				ll.pair_list[j].getFinish().get_z());
		     if ( (++linesdrawn & 1023) == 0) {
			glEnd();
			glBegin(GL_LINES);
			linesdrawn = 0;
		     }
		  }
		  glEnd();
	       }
	    }
	 }
      }
   }
}

// publically accessible
std::pair<coot::dipole, int>
molecule_class_info_t::add_dipole(const std::vector<coot::residue_spec_t> &res_specs,
				  const coot::protein_geometry &geom) {

   int id = -1;
   coot::dipole d;
   std::vector<std::pair<coot::dictionary_residue_restraints_t, CResidue *> > pairs;
   
   for (unsigned int ires=0; ires<res_specs.size(); ires++) { 
      CResidue *residue_p = get_residue(res_specs[ires]);
      
      if (residue_p) {
	 try {
	    std::string res_type = residue_p->GetResName();
	    std::pair<short int, coot::dictionary_residue_restraints_t> rp = 
	       geom.get_monomer_restraints(res_type);
	    if (rp.first) {
	       std::pair<coot::dictionary_residue_restraints_t, CResidue *> p(rp.second, residue_p);
	       pairs.push_back(p);
	    } else {
	       std::cout << "INFO:: no monomer restraints found for "
			 << coot::residue_spec_t(residue_p) << " type: " << res_type << std::endl;
	    } 
	 }
	 catch (std::runtime_error mess) {
	    std::cout << mess.what() << std::endl;
	 }
      } else {
	 std::cout << "   add_dipole: trapped null residue" << std::endl;
      } 
   }

   if (pairs.size() > 0) {
      try { 
	 coot::dipole dl(pairs);
	 dipoles.push_back(dl);
	 id = dipoles.size() -1;
	 d = dl;
      }
      catch (std::runtime_error mess) {
	    std::cout << mess.what() << std::endl;
      }
   }
   
   return std::pair<coot::dipole, int> (d,id); 
}

void
molecule_class_info_t::delete_dipole(int dipole_number) {

   if (dipole_number < int(dipoles.size())) {
      std::vector<coot::dipole>::iterator it;
      int n=0;
      for (it=dipoles.begin(); it!=dipoles.end(); it++) {
	 if (n == dipole_number) {
	    dipoles.erase(it);
	    break;
	 }
	 n++;
      }
   }
}


void
molecule_class_info_t::draw_dipoles() const {

   if (! drawit)
      return;
   
   if (dipoles.size() > 0) { 
      glPushMatrix();
      glLineWidth(2.0);
      std::vector<clipper::Coord_orth> arrow_points;
      arrow_points.push_back(clipper::Coord_orth(0,0,0));
      arrow_points.push_back(clipper::Coord_orth(0.13,0,0));
      arrow_points.push_back(clipper::Coord_orth(0.13,0,0.66));
      arrow_points.push_back(clipper::Coord_orth(0.3,0,0.66));
      arrow_points.push_back(clipper::Coord_orth(0,0,1));
      // shift so that the centre of the arrow is at the origin
      for (unsigned int i=0; i<arrow_points.size(); i++)
	 arrow_points[i] -= clipper::Coord_orth(0,0,0.4); // not 0.5,
                                                          // esthetics

      // Guide-line object
      if (0) { 
	 glColor3f(0.8,0.6,0.4);
	 for (unsigned int i=0; i<dipoles.size(); i++) {
	    clipper::Coord_orth pt = dipoles[i].position();
	    clipper::Coord_orth d = dipoles[i].get_dipole();
	    double sc = 3.0;
	    glBegin(GL_LINES);
	    glVertex3d(pt.x(), pt.y(), pt.z());
	    glVertex3d(pt.x() + d.x() * sc,
		       pt.y() + d.y() * sc,
		       pt.z() + d.z() * sc);
	    glVertex3d(pt.x(), pt.y(), pt.z());
	    glVertex3d(pt.x() - d.x() * sc,
		       pt.y() - d.y() * sc,
		       pt.z() - d.z() * sc);
	    glEnd();
	 }
      }

      glColor3f(0.9,0.6,0.8);
      for (unsigned int i=0; i<dipoles.size(); i++) {
	 clipper::Coord_orth pt = dipoles[i].position();
	 clipper::Coord_orth  d = dipoles[i].get_dipole();
	 clipper::Coord_orth d_unit = dipoles[i].get_unit_dipole();

	 // make an arbitrary vector not parallel to d_unit.
	 // 
	 clipper::Coord_orth arb(0,0.1,0.9);
	 if (d_unit.y() < d_unit.z())
	    arb = clipper::Coord_orth(0.0, 0.9, 0.1);
	 if (d_unit.x() < d_unit.y())
	    arb = clipper::Coord_orth(0.9, 0.0, 0.1);

	 clipper::Coord_orth p1(clipper::Coord_orth::cross(arb, d_unit).unit());
	 clipper::Coord_orth p2(clipper::Coord_orth::cross( p1, d_unit).unit());
	 clipper::Coord_orth p3 = d_unit;
	 
	 GL_matrix m(p1.x(), p1.y(), p1.z(),
		     p2.x(), p2.y(), p2.z(),
		     p3.x(), p3.y(), p3.z());

	 glPushMatrix();
	 glTranslated(pt.x(), pt.y(), pt.z());
	 glMultMatrixf(m.get());
	 glBegin(GL_LINES);

	 // scale the dipole
	 //
	 double ds = sqrt(d.lengthsq());
	 std::vector<clipper::Coord_orth> local_arrow_points = arrow_points;
	 for (unsigned int i=0; i<arrow_points.size(); i++)
	    local_arrow_points[i] = 0.181818181818181818181818181818181818 * ds * arrow_points[i];
	       
	 for (unsigned int i=0; i<local_arrow_points.size()-1; i++) {
	    glVertex3d(local_arrow_points[i].x(),
		       local_arrow_points[i].y(),
		       local_arrow_points[i].z());
	    glVertex3d(local_arrow_points[i+1].x(),
		       local_arrow_points[i+1].y(),
		       local_arrow_points[i+1].z());
	    glVertex3d(-local_arrow_points[i].x(),
		        local_arrow_points[i].y(),
		        local_arrow_points[i].z());
	    glVertex3d(-local_arrow_points[i+1].x(),
		        local_arrow_points[i+1].y(),
		        local_arrow_points[i+1].z());
	 }
	 glEnd();
	 glPopMatrix();
      }
      glPopMatrix();
   }
}

std::string
coot::atom_selection_info_t::name () const {

   std::string s = "An Add Rep atom sel string";
   return s;
}


// return the atom selection and the number of atoms
int
coot::atom_selection_info_t::select_atoms(CMMDBManager *mol) const {

   int SelHnd = -1;
   const char *alt_conf_str = "*";
   if (alt_conf_is_set)
      alt_conf_str = altconf.c_str();
   if (type == BY_ATTRIBUTES) {
      SelHnd = mol->NewSelection();
      mol->SelectAtoms(SelHnd, 0, chain_id.c_str(),
		       resno_start, // starting resno, an int
		       ins_code.c_str(), // any insertion code
		       resno_start, // ending resno
		       ins_code.c_str(), // ending insertion code
		       "*", // any residue name
		       "*", // atom name
		       "*", // elements
		       alt_conf_str  // alt loc.
		       );
   }
   if (type == BY_STRING) {
      SelHnd = mol->NewSelection();
      mol->Select(SelHnd, STYPE_ATOM, atom_selection_str.c_str(), SKEY_NEW);
   } 
   return SelHnd;
}


std::string
coot::atom_selection_info_t::mmdb_string() const {

   std::string s = atom_selection_str;
   if (type == BY_ATTRIBUTES) {
      s = "//";
      s += chain_id;
      s += "/";
      s += coot::util::int_to_string(resno_start);
      s += "-";
      s += coot::util::int_to_string(resno_end);
   } 
   return s;
}


void
coot::additional_representations_t::fill_bonds_box() {

   if (representation_type != coot::BALL_AND_STICK) { 
      atom_selection_container_t atom_sel;

      atom_sel.mol =  (MyCMMDBManager *) mol;
      atom_sel.SelectionHandle = mol->NewSelection();
   
      if (atom_sel_info.type == coot::atom_selection_info_t::BY_ATTRIBUTES) {
      
	 mol->SelectAtoms(atom_sel.SelectionHandle,
			  0, (char *) atom_sel_info.chain_id.c_str(),
			  atom_sel_info.resno_start, (char *) atom_sel_info.ins_code.c_str(),
			  atom_sel_info.resno_end,   (char *) atom_sel_info.ins_code.c_str(),
			  "*", "*", "*", "*");
      }
      if (atom_sel_info.type == coot::atom_selection_info_t::BY_STRING) {
	 mol->Select(atom_sel.SelectionHandle, STYPE_ATOM, 
		     (char *) atom_sel_info.atom_selection_str.c_str(),
		     SKEY_NEW);
      
      }
      mol->GetSelIndex(atom_sel.SelectionHandle,
		       atom_sel.atom_selection,
		       atom_sel.n_selected_atoms);

      if (bonds_box_type == coot::NORMAL_BONDS) { 
	 Bond_lines_container bonds(atom_sel, 1, draw_hydrogens_flag);
	 bonds_box.clear_up();
	 bonds_box = bonds.make_graphical_bonds();
      }
      mol->DeleteSelection(atom_sel.SelectionHandle);
   }
}

std::string
coot::additional_representations_t::info_string() const {

   std::string s("Fat Bonds: ");

   if (representation_type == coot::BALL_AND_STICK) {
     s = "Ball and Stick: "; 
   }

   if (representation_type == coot::STICKS) {
     s = "Sticks: "; 
   }

   if (atom_sel_info.type == coot::atom_selection_info_t::BY_STRING)
      s += atom_sel_info.atom_selection_str;
   if (atom_sel_info.type == coot::atom_selection_info_t::BY_ATTRIBUTES) { 
      s += atom_sel_info.chain_id;
      s += " ";
      s += coot::util::int_to_string(atom_sel_info.resno_start);
      if (atom_sel_info.resno_end != atom_sel_info.resno_start) {
	 s += "-";
	 s += coot::util::int_to_string(atom_sel_info.resno_end);
      }
      s += atom_sel_info.ins_code;
   }
   return s;
} 


int
molecule_class_info_t::add_additional_representation(int representation_type,
						     const int &bonds_box_type_in, 
						     float bonds_width,
						     bool draw_hydrogens_flag,
						     const coot::atom_selection_info_t &info,
						     GtkWidget *display_control_window,
						     const gl_context_info_t &glci) {

   coot::additional_representations_t rep(atom_sel.mol,
					  representation_type,
					  bonds_box_type_in,
					  bonds_width, draw_hydrogens_flag, info);
   add_reps.push_back(rep);
   int n_rep = add_reps.size() -1;
   std::string name = rep.info_string();
   GtkWidget *vbox = display_control_add_reps_container(display_control_window, imol_no);
   display_control_add_reps(vbox, imol_no, n_rep, rep.show_it, rep.bonds_box_type, name);
   if (representation_type == coot::BALL_AND_STICK) {
      int display_list_handle_index = make_ball_and_stick(info.mmdb_string(), 0.12, 0.28, 1, glci);
      if ((display_list_handle_index >= 0) && (display_list_handle_index < display_list_tags.size())) {
	 add_reps[n_rep].add_display_list_handle(display_list_handle_index);
      } 
   }
     
   return n_rep;
}

int
molecule_class_info_t::adjust_additional_representation(int represenation_number, 
							const int &bonds_box_type_in, 
							float bonds_width,
							bool draw_hydrogens_flag,
							const coot::atom_selection_info_t &info, 
							bool show_it_flag_in) {
   return -1;
}


void
molecule_class_info_t::clear_additional_representation(int representation_number) {

   if (add_reps.size() > representation_number) {
      if (representation_number >= 0) {
	 add_reps[representation_number].clear(); 
      } 
   } 
} 

void
molecule_class_info_t::set_show_additional_representation(int representation_number,
							  bool on_off_flag) {

   if (add_reps.size() > representation_number) {
      if (representation_number >= 0) {
	 add_reps[representation_number].show_it = on_off_flag;
	 if (add_reps[representation_number].representation_type == coot::BALL_AND_STICK ||
	     add_reps[representation_number].representation_type == coot::STICKS) { 
	   int dl_index = add_reps[representation_number].display_list_handle;
	   // std::cout << "Ball and stick add rep toggled to " << on_off_flag << std::endl;
	   display_list_tags[dl_index].display_it = on_off_flag;
	 }
      }
   }
}

void
molecule_class_info_t::set_show_all_additional_representations(bool on_off_flag) {
   int n_reps = add_reps.size();
   for (unsigned int i=0; i<n_reps; i++)
      set_show_additional_representation(i, on_off_flag);
} 



// Return a pair.
// 
// If first string of length 0 on error to construct dataname(s).
std::pair<std::string, std::string>
molecule_class_info_t::make_import_datanames(const std::string &f_col_in,
					     const std::string &phi_col_in,
					     const std::string &weight_col_in,
					     int use_weights) const {

   // If use_weights return 2 strings, else set something useful only for pair.first

   std::string f_col = f_col_in;
   std::string phi_col = phi_col_in;
   std::string weight_col = weight_col_in;

#ifdef WINDOWS_MINGW
   std::string::size_type islash_f   =      coot::util::intelligent_debackslash(f_col).find_last_of("/");
   std::string::size_type islash_phi =    coot::util::intelligent_debackslash(phi_col).find_last_of("/");
#else
   std::string::size_type islash_f   =      f_col.find_last_of("/");
   std::string::size_type islash_phi =    phi_col.find_last_of("/");
#endif // MINGW
   short int label_error = 0; 

   if (islash_f != std::string::npos) {
      // f_col is of form e.g. xxx/yyy/FWT
      if (f_col.length() > islash_f)
	 f_col = f_col.substr(islash_f+1);
      else
	 label_error = 1;
   }

   if (islash_phi != std::string::npos) {
      // phi_col is of form e.g. xxx/yyy/PHWT
      if (phi_col.length() > islash_phi)
	 phi_col = phi_col.substr(islash_phi+1);
      else
	 label_error = 1;
   }

   if (use_weights) { 
#ifdef WINDOWS_MINGW
      std::string::size_type islash_fom = coot::util::intelligent_debackslash(weight_col).find_last_of("/");
#else
      std::string::size_type islash_fom = weight_col.find_last_of("/");
#endif
      if (islash_fom != std::string::npos) {
	 // weight_col is of form e.g. xxx/yyy/WT
	 if (weight_col.length() > islash_fom)
	    weight_col = weight_col.substr(islash_fom+1);
	 else
	    label_error = 1;
      }
   }

  
   std::pair<std::string, std::string> p("", "");

   if (!label_error) {
      std::string no_xtal_dataset_prefix= "/*/*/";
      if (use_weights) { 
	 p.first  = no_xtal_dataset_prefix + "[" +   f_col + " " +      f_col + "]";
	 p.second = no_xtal_dataset_prefix + "[" + phi_col + " " + weight_col + "]";
      } else {
	 p.first  = no_xtal_dataset_prefix + "[" +   f_col + " " + phi_col + "]";
      }
   }
   return p;
}



void
molecule_class_info_t::filter_by_resolution(clipper::HKL_data< clipper::datatypes::F_phi<float> > *fphidata,
					    const float &reso_low,
					    const float &reso_high) const {

   float inv_low  = 1.0/(reso_low*reso_low);
   float inv_high = 1.0/(reso_high*reso_high);
   int n_data = 0;
   int n_reset = 0;

      
   for (clipper::HKL_info::HKL_reference_index hri = fphidata->first(); !hri.last(); hri.next()) {
//        std::cout << "high: " << inv_high << " low: " << inv_low
//  		<< " data: " << hri.invresolsq() << std::endl;
      n_data++;

      if ( hri.invresolsq() > inv_low &&
	   hri.invresolsq() < inv_high) {
      } else {
	 (*fphidata)[hri].f() = 0.0;
	 n_reset++;
      } 
   }
   std::cout << "Chopped " << n_reset << " data out of " << n_data << std::endl;
}


void
molecule_class_info_t::test_label_symm_atom(int i) {
   //
 
   // same test as has_model():
   if (has_model()) {
      
      if (i < labelled_symm_atom_index_list.size()) { 

	 int iatom_index = labelled_symm_atom_index_list[i];

	 if (iatom_index < atom_sel.n_selected_atoms) { 
	    std::pair <symm_trans_t, Cell_Translation> st = labelled_symm_atom_symm_trans_[i];
	 
	    std::string label =
	       make_symm_atom_label_string(atom_sel.atom_selection[iatom_index], st.first);

	    GLfloat blueish[3] = { 0.7, 0.7, 1.0 };
	 
	    glColor3fv(blueish);
	 
	    coot::Cartesian symm_point = translate_atom_with_pre_shift(atom_sel, iatom_index, st);

	    glRasterPos3f(symm_point.get_x(),
			  symm_point.get_y()+0.02,
			  symm_point.get_z()+0.02);
	 
	    printString(label);
	 }
      }
   }
}

void
molecule_class_info_t::label_symm_atom(int i, symm_trans_t symm_trans) {
   //

   // same test as has_model():
   if (atom_sel.n_selected_atoms > 0 ) { 

      if (atom_sel.n_selected_atoms > 0) {

	 std::string label =
	    make_symm_atom_label_string(atom_sel.atom_selection[i], symm_trans);

	 GLfloat blueish[3] = { 0.7, 0.8, 1.0 };

	 glColor3fv(blueish);

	 std::cout << "label_symm_atom :" << symm_trans << std::endl;
	 coot::Cartesian symm_point =
	    translate_atom(atom_sel, i, symm_trans);
            
	 glRasterPos3f(symm_point.get_x(),
		       symm_point.get_y()+0.02,
		       symm_point.get_z()+0.02);

	 std::cout << "adding symm label: " << label << std::endl;
	 printString(label);
	 std::cout << "done  symm label: in  label_symm_atom(..): " << label << std::endl;

      }
   }   
}

// Put a label at the ith atom of mol_class_info::atom_selection. 
//
void
molecule_class_info_t::label_atom(int i, int brief_atom_labels_flag) {

   if (has_model()) { 

      if (i < atom_sel.n_selected_atoms) { 

	 PCAtom atom = (atom_sel.atom_selection)[i];

	 if (atom) { 

	    std::string label = make_atom_label_string(atom, brief_atom_labels_flag);
      
	    // GLfloat white[3] = { 1.0, 1.0, 1.0 };
	    GLfloat pink[3] =  { graphics_info_t::font_colour.red,
				 graphics_info_t::font_colour.green,
				 graphics_info_t::font_colour.blue };
      
	    // glClear(GL_COLOR_BUFFER_BIT);
	    glColor3fv(pink);
	    // glShadeModel (GL_FLAT);
      
	    glRasterPos3f((atom)->x, (atom)->y+0.02, (atom)->z +0.02);
	    printString(label);
      
	 }
      } else { 
	 std::cout << "INFO:: trying to label atom out of range: " 
		   << i << " " << atom_sel.n_selected_atoms 
		   << " Removing label\n";
	 unlabel_atom(i);
      }
   }
}


void
molecule_class_info_t::set_have_unit_cell_flag_maybe() {
   
   CMMDBCryst *cryst_p = atom_sel.mol->get_cell_p();

   mat44 my_matt;

   int err = cryst_p->GetTMatrix(my_matt, 0, 0, 0, 0);

   if (err != 0) {
      have_unit_cell = 0;
      std::cout << "No Symmetry for this model" << std::endl;
   } else { 
      have_unit_cell = 1;
   }
}

// This function is not used in anger, right?
//
// (a debugging function).
// 
void
check_static_vecs_extents() {

   //
   int imol = 0;
   
   graphics_info_t g;

   cout << "checking extents of the " << g.molecules[imol].n_draw_vectors
	<< " in the graphics static" << endl;

   coot::Cartesian first, second;
   float max_x = -9999, min_x = 9999;
   float max_y = -9999, min_y = 9999;
   float max_z = -9999, min_z = 9999;
   
   for (int i=0; i<g.molecules[imol].n_draw_vectors; i++) {
      first  = g.molecules[imol].draw_vectors[i].getStart();
      second = g.molecules[imol].draw_vectors[i].getFinish();

      if (first.get_x() < min_x) min_x = first.get_x();
      if (first.get_y() < min_y) min_y = first.get_y();
      if (first.get_z() < min_z) min_z = first.get_z();

      if (second.get_x() < min_x) min_x = second.get_x();
      if (second.get_y() < min_y) min_y = second.get_y();
      if (second.get_z() < min_z) min_z = second.get_z();

      if (first.get_x() > max_x) max_x = first.get_x();
      if (first.get_y() > max_y) max_y = first.get_y();
      if (first.get_z() > max_z) max_z = first.get_z();

      if (second.get_x() > max_x) max_x = second.get_x();
      if (second.get_y() > max_y) max_y = second.get_y();
      if (second.get_z() > max_z) max_z = second.get_z();
      
   }
   cout << min_x << " " << max_x << endl
	<< min_y << " " << max_y << endl
	<< min_z << " " << max_z << endl;
}


void
molecule_class_info_t::makebonds(float min_dist, float max_dist) {

   Bond_lines_container bonds(atom_sel, min_dist, max_dist);
   bonds_box.clear_up();
   bonds_box = bonds.make_graphical_bonds();
   bonds_box_type = coot::NORMAL_BONDS;
}

void
molecule_class_info_t::makebonds(float max_dist) {
   
   Bond_lines_container bonds(atom_sel, max_dist);

   // bonds.check(); 

   bonds_box.clear_up();
   bonds_box = bonds.make_graphical_bonds();

   // cout << "makebonds bonds_box.num_colours "
   // << bonds_box.num_colours << endl;
}

void
molecule_class_info_t::makebonds() {

   int do_disulphide_flag = 1;
   Bond_lines_container bonds(atom_sel, do_disulphide_flag, draw_hydrogens_flag);
   bonds_box.clear_up();
   bonds_box = bonds.make_graphical_bonds();
   bonds_box_type = coot::NORMAL_BONDS;

}

void
molecule_class_info_t::make_ca_bonds(float min_dist, float max_dist) {

   Bond_lines_container bonds;
   bonds.do_Ca_bonds(atom_sel, min_dist, max_dist);
   bonds_box = bonds.make_graphical_bonds();
   bonds_box_type = coot::CA_BONDS;
   // std::cout << "ca: bonds_box_type is now " << bonds_box_type << std::endl;

}

void
molecule_class_info_t::make_ca_bonds() { 
   make_ca_bonds(2.4, 4.7); 
}

void
molecule_class_info_t::make_ca_plus_ligands_bonds() { 

   Bond_lines_container bonds;
   bonds.do_Ca_plus_ligands_bonds(atom_sel, 2.4, 4.7);
   bonds_box = bonds.make_graphical_bonds();
   bonds_box_type = coot::CA_BONDS_PLUS_LIGANDS;
   
   // std::cout << "ca: bonds_box_type is now " << bonds_box_type << std::endl;
}

void
molecule_class_info_t::make_colour_by_chain_bonds(short int change_c_only_flag) {
   // 
   Bond_lines_container bonds;
   bonds.do_colour_by_chain_bonds(atom_sel, draw_hydrogens_flag, change_c_only_flag);
   bonds_box = bonds.make_graphical_bonds();
   bonds_box_type = coot::COLOUR_BY_CHAIN_BONDS;

   if (graphics_info_t::glarea) 
      graphics_info_t::graphics_draw();
} 

void
molecule_class_info_t::make_colour_by_molecule_bonds() { 

   // 
   Bond_lines_container bonds;
   bonds.do_colour_by_molecule_bonds(atom_sel, draw_hydrogens_flag);
   bonds_box = bonds.make_graphical_bonds();
   bonds_box_type = coot::COLOUR_BY_MOLECULE_BONDS;

   if (graphics_info_t::glarea) 
      graphics_info_t::graphics_draw();

} 



void 
molecule_class_info_t::make_bonds_type_checked() { 

   if (bonds_box_type == coot::NORMAL_BONDS)
      makebonds();
   if (bonds_box_type == coot::CA_BONDS)
      make_ca_bonds();
   if (bonds_box_type == coot::COLOUR_BY_CHAIN_BONDS)
      // Baah, we have to use the static in graphics_info_t here as it
      // is not a per-molecule property.
      make_colour_by_chain_bonds(graphics_info_t::rotate_colour_map_on_read_pdb_c_only_flag);
   if (bonds_box_type == coot::COLOUR_BY_MOLECULE_BONDS)
      make_colour_by_molecule_bonds();
   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS)
      make_ca_plus_ligands_bonds();
   if (bonds_box_type == coot::BONDS_NO_WATERS)
      bonds_no_waters_representation();
   if (bonds_box_type == coot::BONDS_SEC_STRUCT_COLOUR)
      bonds_sec_struct_representation();
   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR)
      ca_plus_ligands_sec_struct_representation();
   if (bonds_box_type == coot::COLOUR_BY_RAINBOW_BONDS)
      ca_plus_ligands_rainbow_representation();
   if (bonds_box_type == coot::COLOUR_BY_OCCUPANCY_BONDS)
      occupancy_representation();
   if (bonds_box_type == coot::COLOUR_BY_B_FACTOR_BONDS)
      b_factor_representation();

   // bleugh. But if we don't do this here, where *do* we do it?
   // Should the glci be passed to make_bonds_type_checked()?  Urgh.
   // That is called from many places....
   // 
   gl_context_info_t glci(graphics_info_t::glarea, graphics_info_t::glarea_2);
   
   update_additional_representations(glci);
   update_fixed_atom_positions();
   update_ghosts();
}

void
molecule_class_info_t::update_additional_representations(const gl_context_info_t &gl_info) {

   for (unsigned int i=0; i<add_reps.size(); i++) {
      // make_ball_and_stick is not available from inside an add_rep,
      // so we do it outside.
      if (add_reps[i].representation_type != coot::BALL_AND_STICK) {
	 add_reps[i].update_self();
      } else {

	 // let's remove the old/current dlo from the tags list before we add a new one.
	 
	 int old_handle = add_reps[i].display_list_handle;
	 remove_display_list_object_with_handle(old_handle);
	 int handle = make_ball_and_stick(add_reps[i].atom_sel_info.mmdb_string(), 0.11, 0.24, 1, gl_info);

	 std::cout << " update a ball and stick rep " << i << " "
		   << add_reps[i].show_it << std::endl;
	 if ((handle >= 0) && (handle < display_list_tags.size())) 
	    add_reps[i].update_self_display_list_entity(handle);
	 display_list_tags[handle].display_it = add_reps[i].show_it;
      }
   }
}


void
molecule_class_info_t::remove_display_list_object_with_handle(int handle_index) {

//    for (unsigned int i=0; i<display_list_tags.size(); i++) { 
//       std::cout << "   display_list_tags: index " << i << " tag " << display_list_tags[i].tag;
//       if (display_list_tags[i].is_closed)
// 	 std::cout << " closed";
//       std::cout << std::endl;
//    }
   
//    std::cout << "closing display_list_tags number..." << handle_index
// 	     << " which has GL tag " << display_list_tags[handle_index].tag << std::endl;
   display_list_tags[handle_index].close_yourself();
} 

void
molecule_class_info_t::update_mols_in_additional_representations() {
   
   for (unsigned int i=0; i<add_reps.size(); i++) {
      add_reps[i].change_mol(atom_sel.mol);
   }
}


void
molecule_class_info_t::update_fixed_atom_positions() {

   fixed_atom_positions.clear();
   bool found_match = 0; 
   for(unsigned int i=0; i<fixed_atom_specs.size(); i++) {
      int ifast_index = fixed_atom_specs[i].int_user_data;
      if (ifast_index != -1) {
	 if (ifast_index < atom_sel.n_selected_atoms) {
	    CAtom *at = atom_sel.atom_selection[ifast_index];
	    if (fixed_atom_specs[i].matches_spec(at)) {
	       found_match = 1;
	       coot::Cartesian pos(at->x, at->y, at->z);
	       fixed_atom_positions.push_back(pos);
	    } 
	 } 
      }
      if (! found_match) {
	 // use a slower method to find atom
	 int idx = full_atom_spec_to_atom_index(fixed_atom_specs[i]);
	 if (idx != -1) {
	    CAtom *at = atom_sel.atom_selection[idx];
	    if (fixed_atom_specs[i].matches_spec(at)) {
	       coot::Cartesian pos(at->x, at->y, at->z);
	       fixed_atom_positions.push_back(pos);
	    }
	 }
      }
   }
}

std::vector<coot::atom_spec_t>
molecule_class_info_t::get_fixed_atoms() const {
   return fixed_atom_specs;
}


// export the molecule in atom_selection_container_t atom_sel;
// 
int
molecule_class_info_t::export_coordinates(std::string filename) const { 

   //
   int err = atom_sel.mol->WritePDBASCII((char *)filename.c_str()); 
   
   if (err) { 
      std::cout << "WARNING:: export coords: There was an error in writing "
		<< filename << std::endl; 
      std::cout << GetErrorDescription(err) << std::endl;
      graphics_info_t g;
      std::string s = "ERROR:: writing coordinates file ";
      s += filename;
      g.statusbar_text(s);
   } else {
      std::string s = "INFO:: coordinates file ";
      s += filename;
      s += " saved successfully";
      graphics_info_t g;
      g.statusbar_text(s);
   } 
   return err;
}

// Perhaps this should be a util function?
CMMDBManager
*molecule_class_info_t::get_residue_range_as_mol(const std::string &chain_id,
						 int resno_start,
						 int resno_end) const {

   int imod = 1;

   CMMDBManager *mol_new = new CMMDBManager;
   CModel *model_new = new CModel;
   CChain *chain_new = new CChain;
   
   realtype cell[6];
   realtype vol;
   int orthcode;
   char *spacegroup_str = atom_sel.mol->GetSpaceGroup();
   atom_sel.mol->GetCell(cell[0], cell[1], cell[2],
			 cell[3], cell[4], cell[5],
			 vol, orthcode); 
   mol_new->SetCell(cell[0], cell[1], cell[2],
		    cell[3], cell[4], cell[5], orthcode);
   mol_new->SetSpaceGroup(spacegroup_str);

   CModel *model_p = atom_sel.mol->GetModel(imod);
   CChain *chain_p;
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      if (std::string(chain_p->GetChainID()) == chain_id) { 
	 int nres = chain_p->GetNumberOfResidues();
	 PCResidue residue_p;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    if (residue_p->GetSeqNum() >= resno_start) {
	       if (residue_p->GetSeqNum() <= resno_end) {
		  CResidue *res_new =
		     coot::util::deep_copy_this_residue(residue_p, "", 1);
		  chain_new->AddResidue(res_new);
	       }
	    }
	 }
      }
   }

   chain_new->SetChainID(chain_id.c_str());
   model_new->AddChain(chain_new);
   mol_new->AddModel(model_new);
   mol_new->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
   mol_new->FinishStructEdit();
   return mol_new;
} 


std::string
molecule_class_info_t::make_symm_atom_label_string(PCAtom atom, symm_trans_t symm_trans) const {

   std::string s = make_atom_label_string(atom, 0);
   s += " ";
   s += symm_trans.str(graphics_info_t::symmetry_atom_labels_expanded_flag);
   return s;
}

// This took one minute to write.
//
// And then another 2 to add the alt location code (and test it).
// 
std::string
molecule_class_info_t::make_atom_label_string(PCAtom atom, int brief_atom_labels_flag) const {

   char *chain_id  = atom->GetChainID();
   char *res_name  = atom->GetResName();
   int   res_no    = atom->GetSeqNum();
   char *atom_name = atom->name;
   char *ins_code  = atom->GetInsCode();

   // format: atom_name/res_no res_name/chain_id
   // new format: atom_name,alt_conf/res_no res_name/chain_id if altconf != ""
   graphics_info_t g;

   std::string s(atom_name);
   std::string alt_loc(atom->altLoc);
   if (alt_loc != "") {
      int slen = s.length();
      if (slen > 0) {
	 if (s[slen-1] == ' ') {
	    s = s.substr(0,slen-1) + ",";
	 } else {
	    s += ",";
	 }
      } else {
	 s += ",";
      }
      s += alt_loc;
   }

   if (brief_atom_labels_flag) {
      s += g.int_to_string(res_no);
      if (strlen(ins_code) > 0) {
	 s += ins_code;
	 s += " ";
      }
      s += chain_id;
   } else { 
      s += "/";
      s += g.int_to_string(res_no);
      s += ins_code;
      s += " ";
      s += res_name;
      s += "/";
      s += chain_id;
   }

   return s;
}

// Don't use this function.
// 
// For labelling from guile and replacing coordinates and others.
// 
// Return -1 on failure to find match
int 
molecule_class_info_t::atom_spec_to_atom_index(std::string chain, int resno, 
					       std::string atom_name) const { 

   int iatom_index = -1; 
   int selHnd = atom_sel.mol->NewSelection();

   atom_sel.mol->SelectAtoms(selHnd, 0, (char *) chain.c_str(), 
			    resno, "*", // start, insertion code
			    resno, "*", // end, insertion code
			    "*", // residue name
			    (char *) atom_name.c_str(),
			    "*", // elements
			    "*"); // alt locs

   int nSelAtoms;
   PPCAtom local_SelAtom; 
   atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);

   if (nSelAtoms == 0) {
      std::cout << "Sorry (atom_spec_to_atom_index): Could not find " << atom_name << "/"
		<< resno << "/" << chain << " in this molecule: ("
		<<  imol_no << ") " << name_ << std::endl; 

      // debug:
      selHnd = atom_sel.mol->NewSelection();
      
      atom_sel.mol->SelectAtoms(selHnd, 0, "*",
				ANY_RES, "*", // start, insertion code
				ANY_RES, "*", // end, insertion code
				"*", // residue name
				(char *) atom_name.c_str(),
				"*", // elements
				"*"); // alt locs

      atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);

      std::cout << "There were " << nSelAtoms << " atoms with resno "
		<< resno << std::endl;

      
   } else {
      // compare pointers
      for (int i=0; i<atom_sel.n_selected_atoms; i++) {
	 if (atom_sel.atom_selection[i] == local_SelAtom[0]) {
	    iatom_index = i;
	    break; 
	 }
      }
   } 
   return iatom_index; 
}			    

int
molecule_class_info_t::atom_to_atom_index(CAtom *at) const {

   int iatom_index_udd = -1;
   int ic;
   if (at->GetUDData(atom_sel.UDDAtomIndexHandle, ic) == UDDATA_Ok) {
      iatom_index_udd = ic;
   }
   if (iatom_index_udd == -1)
      iatom_index_udd=full_atom_spec_to_atom_index(at);
   return iatom_index_udd;
}

int
molecule_class_info_t::full_atom_spec_to_atom_index(const coot::atom_spec_t &atom_spec) const {

   return full_atom_spec_to_atom_index(atom_spec.chain,
				       atom_spec.resno,
				       atom_spec.insertion_code,
				       atom_spec.atom_name,
				       atom_spec.alt_conf);

} 


// return -1 on no atom found.
int
molecule_class_info_t::full_atom_spec_to_atom_index(const std::string &chain,
						    int resno,
						    const std::string &insertion_code,
						    const std::string &atom_name,
						    const std::string &alt_conf) const {

   int iatom_index = -1; 

   // some protection for null molecule.
   if (! atom_sel.mol) {
      std::cout << "ERROR:: null molecule " << imol_no << " " << atom_sel.mol
		<< " (in full_atom_spec_to_atom_index)" << std::endl;
      return -1;
   }

   int selHnd = atom_sel.mol->NewSelection();
   int idx = 0;

   atom_sel.mol->SelectAtoms(selHnd, 0, (char *) chain.c_str(), 
			    resno, insertion_code.c_str(), // start, insertion code
			    resno, insertion_code.c_str(), // end, insertion code
			    "*", // residue name
			    atom_name.c_str(),
			    "*", // elements
			    alt_conf.c_str()); // alt locs

   int nSelAtoms;
   PPCAtom local_SelAtom; 
   atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);

//    std::cout << "DEBUG:: full_atom_spec_to_atom_index for :" << chain << ": "
// 	     << resno << " :" << insertion_code << ": :" 
//  	     << atom_name << ": :" << alt_conf << ": finds " << nSelAtoms <<  " atoms\n";

   if (nSelAtoms == 0) { 

      std::cout << "Sorry (full_atom_spec_to_atom_index) Could not find "
		<< atom_name << "," << "\"" << alt_conf  << "\"" << "/"
		<< resno << insertion_code << "/" << chain << " in this molecule: ("
		<<  imol_no << ") " << name_ << std::endl; 
      // debug:
      int selHnd2 = atom_sel.mol->NewSelection();
      
      atom_sel.mol->SelectAtoms(selHnd2, 0, 
				chain.c_str(),
				resno, "*", // start, insertion code
				resno, "*", // end, insertion code
				"*", // residue name
				"*", // atom name
				"*", // elements
				"*"); // alt locs

      atom_sel.mol->GetSelIndex(selHnd2, local_SelAtom, nSelAtoms);

      std::cout << "There were " << nSelAtoms << " atoms in that residue:\n";
      for (int i=0; i<nSelAtoms; i++) { 
	 std::cout << "      " << local_SelAtom[i] << "\n";
      }

      atom_sel.mol->DeleteSelection(selHnd2);

   } else { 

      if (nSelAtoms != 1) {
	 // the wildcard atom selection case "*HO2"
	 short int found = 0;
	 for (int i=0; i<nSelAtoms; i++) { 
	    if (std::string(local_SelAtom[i]->GetChainID()) == chain) { 
	       if (local_SelAtom[i]->residue->seqNum == resno) { 
		  if (std::string(local_SelAtom[i]->GetInsCode()) == insertion_code) { 
		     if (std::string(local_SelAtom[i]->name) == atom_name) { 
			if (std::string(local_SelAtom[i]->altLoc) == alt_conf) { 
			   found = 0;
			   idx = i;
			   break;
			} 
		     }
		  }
	       }
	    }
	 } 
      }
      
      int iatom_index_udd = -1;
      int ic;
      if (local_SelAtom[idx]->GetUDData(atom_sel.UDDAtomIndexHandle, ic) == UDDATA_Ok) {
	 iatom_index_udd = ic;
      }
      iatom_index = iatom_index_udd;
   }
   atom_sel.mol->DeleteSelection(selHnd); // Oh dear, this should have
					  // been in place for years
					  // (shouldn't it?) 20071121
   return iatom_index; 
}
//       // compare pointers
//       for (int i=0; i<atom_sel.n_selected_atoms; i++) {
// 	 if (atom_sel.atom_selection[i] == local_SelAtom[0]) {
// 	    iatom_index = i;
// 	    break; 
// 	 }
//       }
	 // 	 if (iatom_index != iatom_index_udd) { 
	 // 	    std::cout << "ERROR: atom indexes (UDD test) dont match "
	 // 		      << iatom_index << " " << iatom_index_udd << std::endl;
	 // 	 }


// Does atom at from moving atoms match atom_sel.atom_selection[this_mol_index_maybe]?
// or has atom_sel changed in the mean time?
bool
molecule_class_info_t::moving_atom_matches(CAtom *at, int this_mol_index_maybe) const {

   bool matches = 0;
   if (atom_sel.n_selected_atoms > 0) { 
      if (this_mol_index_maybe >= atom_sel.n_selected_atoms) {
	 return 0;
      } else {
	 std::string atom_name_mov = at->name;
	 std::string ins_code_mov  = at->GetInsCode();
	 std::string alt_conf_mov  = at->altLoc;
	 std::string chain_id_mov  = at->GetChainID();
	 int resno_mov = at->GetSeqNum();

	 std::string atom_name_ref = atom_sel.atom_selection[this_mol_index_maybe]->name;
	 std::string ins_code_ref  = atom_sel.atom_selection[this_mol_index_maybe]->GetInsCode();
	 std::string alt_conf_ref  = atom_sel.atom_selection[this_mol_index_maybe]->altLoc;
	 std::string chain_id_ref  = atom_sel.atom_selection[this_mol_index_maybe]->GetChainID();
	 int resno_ref = atom_sel.atom_selection[this_mol_index_maybe]->GetSeqNum();

	 if (atom_name_ref == atom_name_mov) {
	    if (ins_code_ref == ins_code_mov) {
	       if (resno_ref == resno_mov) {
		  if (alt_conf_ref == alt_conf_mov) {
		     matches = 1;
		  }
	       }
	    }
	 }
      }
   }
   return matches;
}



// Another attempt at that:
// 
// find "by hand" the atom with the given characteristics in
// the atom selection.
//
// return -1 if atom not found.
// Note we have to search for " CA " etc
// 
int molecule_class_info_t::atom_index(const char *chain_id, int iresno, const char *atom_id) {

   int n = atom_sel.n_selected_atoms;
   for (int i=0; i<n; i++) {
      if ( ( ! strcmp(atom_id,atom_sel.atom_selection[i]->name) ) &&
	   (atom_sel.atom_selection[i]->residue->seqNum == iresno)  &&
	   ( ! strcmp(chain_id,atom_sel.atom_selection[i]->residue->GetChainID()) )
	   ) {
	 return i;
      }
   }

   return -1; 
}
int 
molecule_class_info_t::atom_index_first_atom_in_residue(const std::string &chain_id,
		                  		        int iresno,
				                        const std::string &ins_code) const {

   bool tacf = 0; // test_alt_conf_flag
   return atom_index_first_atom_in_residue_internal(chain_id, iresno, ins_code, "", tacf); 
}

int 
molecule_class_info_t::atom_index_first_atom_in_residue(const std::string &chain_id,
		                  		        int iresno,
				                        const std::string &ins_code,
							const std::string &alt_conf) const {
   
   bool tacf = 1; // test_alt_conf_flag
   return atom_index_first_atom_in_residue_internal(chain_id, iresno, ins_code, alt_conf, tacf); 

}

int 
molecule_class_info_t::atom_index_first_atom_in_residue_internal(const std::string &chain_id,
								 int iresno,
								 const std::string &ins_code,
								 const std::string &alt_conf,
								 bool test_alt_conf_flag) const {

   int index = -1; // failure
   int selHnd = atom_sel.mol->NewSelection();
   int nSelResidues;
   PPCResidue SelResidues;
   atom_sel.mol->Select(selHnd, STYPE_RESIDUE, 1,
			(char *) chain_id.c_str(), 
			iresno, (char *) ins_code.c_str(),
			iresno, (char *) ins_code.c_str(),
			"*",  // residue name
			"*",  // Residue must contain this atom name?
			"*",  // Residue must contain this Element?
			"*",  // altLocs
			SKEY_NEW // selection key
			);
   atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
   if (nSelResidues > 0) {
      int ic = -1;
      for (int ires=0; ires<nSelResidues; ires++) {
	 int natoms;
	 PPCAtom residue_atoms;
	 SelResidues[ires]->GetAtomTable(residue_atoms, natoms);
	 for (int iatom=0; iatom<natoms; iatom++) {
	    if (test_alt_conf_flag == 0
		|| alt_conf == residue_atoms[iatom]->altLoc) { 
	       if (residue_atoms[iatom]->GetUDData(atom_sel.UDDAtomIndexHandle, ic) == UDDATA_Ok) {
		  index = ic;
		  break;
	       }
	    }
	 }
	 if (index > -1)
	    break;
      }
   }
   atom_sel.mol->DeleteSelection(selHnd);
   // std::cout << "DEBUG:: atom_index_first_atom_in_residue returns " << index << std::endl;
   return index;

}

// Put the regularization results back into the molecule:
//
//// Recall that regularized_asc contains an atom_selection_container_t
// with the new coordinates in.  the mol contains all the molecule and
// the atom_selection contains just the moving parts.
//
void
molecule_class_info_t::replace_coords(const atom_selection_container_t &asc,
				      bool change_altconf_occs_flag,
				      bool replace_coords_with_zero_occ_flag) {

   int n_atom = 0;
   int tmp_index;

   make_backup();

   // debug::
   if (0) { 
      std::cout << "DEBUG:: --------------- replace_coords replacing "
		<< asc.n_selected_atoms << " atoms " << std::endl;
      for (int i=0; i<asc.n_selected_atoms; i++) {
	 CAtom *atom = asc.atom_selection[i];
	 bool is_ter_state = atom->isTer();
	 std::cout << "DEBUG:: in replace_coords, intermediate atom: chain-id :"
		   << atom->residue->GetChainID() <<  ": "
		   << atom->residue->seqNum << " inscode :" 
		   << atom->GetInsCode() << ": name :" 
		   << atom->name << ": altloc :"
		   << atom->altLoc << ": occupancy: "
		   << atom->occupancy << " :"
		   << " ter state: " << is_ter_state << std::endl;
      }
   }

   // For each atom in the new set of atoms:
   // 
   for (int i=0; i<asc.n_selected_atoms; i++) {
      int idx = -1;
      CAtom *atom = asc.atom_selection[i];
      if (! atom->isTer()) { 
//       idx = atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
// 				    atom->residue->seqNum,
// 				    std::string(atom->name));
	 if (asc.UDDOldAtomIndexHandle >= 0) { // OK for fast atom indexing
	    if (atom->GetUDData(asc.UDDOldAtomIndexHandle, tmp_index) == UDDATA_Ok) {
	       if (tmp_index >= 0) { 
		  if (moving_atom_matches(atom, tmp_index)) { 
		     // std::cout << "DEBUG:: successfully found old atom index" << std::endl;
		     idx = tmp_index;
		  } else {
		     // std::cout << "DEBUG:: atom index mismatch (this molecule was changed)"
		     // << std::endl;
		     idx = full_atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
							atom->residue->seqNum,
							std::string(atom->GetInsCode()),
							std::string(atom->name),
							std::string(atom->altLoc));
		     // std::cout << "full_atom_spec_to_atom_index gives index: " << idx << std::endl;
		  }
	       } else {
		  // This shouldn't happen.
		  std::cout << "Good Handle, bad index found for old atom: specing" << std::endl;
		  idx = full_atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
						     atom->residue->seqNum,
						     std::string(atom->GetInsCode()),
						     std::string(atom->name),
						     std::string(atom->altLoc));
	       }
	    } else { 
	       std::cout << "ERROR:: non-bad handle (" << asc.UDDOldAtomIndexHandle 
			 <<  "), bad GetUDData for this atom " << std::endl;
	    } 
	 } else {
	    //  	 std::cout << "DEBUG:: asc.UDDOldAtomIndexHandle is " 
	    //   		   << asc.UDDOldAtomIndexHandle << " using full atom spec to atom index..."
	    // 		   << std::endl;
	    
	    idx = full_atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
					       atom->residue->seqNum,
					       std::string(atom->GetInsCode()),
					       std::string(atom->name),
					       std::string(atom->altLoc));
	    if (idx == -1) {
	       std::cout << "DEBUG:: idx: " << idx << "\n";
	       std::cout << "ERROR:: failed to find spec for chain-id :"
			 << std::string(atom->residue->GetChainID()) <<  ": "
			 << atom->residue->seqNum << " inscode :" 
			 << std::string(atom->GetInsCode()) << ": name :" 
			 << std::string(atom->name) << ": altloc :"
			 << std::string(atom->altLoc) << ":" << std::endl;
	    }
	 }

	 if (change_altconf_occs_flag) { 
	    if (idx >= 0) {
	       n_atom++;
	       CAtom *mol_atom = atom_sel.atom_selection[idx];
	       float atom_occ = atom->occupancy;
	       // if this is a shelx molecule, then we don't change
	       // occupancies this way.  We do it by changing the FVAR
	       if (is_from_shelx_ins_flag) { 
		  atom_occ = mol_atom->occupancy;

		  // OK, one more go.  We have an occupancy of 31 or -31
		  // say.  Now, the alt conf atoms has been immmediately
		  // added with the old occupancy for the actual FVAR number
		  // - this happens before we get to twiddle the occupancy
		  // slider.  So here we have to find out the index of the
		  // replaced atom and set it's fvar to whatever the slider
		  // value had been set to.

		  int fvar_number = coot::ShelxIns::shelx_occ_to_fvar(atom_occ);
		  if (fvar_number > 1) { 
		     // 	       std::cout << "DEBUG:: replace_coords: setting fvar number "
		     // 			 <<  fvar_number << " (generated from occ " << atom_occ << ") to "
		     // 			 << graphics_info_t::add_alt_conf_new_atoms_occupancy << std::endl;
		     shelxins.set_fvar(fvar_number, graphics_info_t::add_alt_conf_new_atoms_occupancy);
		  }

		  if (movable_atom(mol_atom, replace_coords_with_zero_occ_flag))
		     mol_atom->SetCoordinates(atom->x,
					      atom->y,
					      atom->z,
					      atom_occ,
					      mol_atom->tempFactor);
	       } else { 
		  if (movable_atom(mol_atom, replace_coords_with_zero_occ_flag))
		     mol_atom->SetCoordinates(atom->x,
					      atom->y,
					      atom->z,
					      atom_occ,
					      mol_atom->tempFactor);
	       }

	       // similarly we adjust occupancy if this is not a shelx molecule
	       if (! is_from_shelx_ins_flag) { 
		  adjust_occupancy_other_residue_atoms(mol_atom, mol_atom->residue, 0);
	       } 
	       // std::cout << atom << " coords replace " << idx << " " << mol_atom << std::endl;
	    } else {
	       std::cout << "ERROR:: bad atom index in replace_coords replacing atom: "
			 << atom << std::endl;
	    }
	 } else {

	    // don't change alt confs.

	    // 	    std::cout << "DEBUG:: no change of alt conf occs for atom index "
	    // 		      << idx << std::endl;

	    if (idx != -1 ) {  // enable this text when fixed.
	       CAtom *mol_atom = atom_sel.atom_selection[idx];
	       if (movable_atom(mol_atom, replace_coords_with_zero_occ_flag)) { 
		  mol_atom->SetCoordinates(atom->x,
					   atom->y,
					   atom->z,
					   mol_atom->occupancy,
					   mol_atom->tempFactor);
		  n_atom++;
	       }
	    }
	 }
      }
   }
   std::cout << "INFO:: replace_coords: " << n_atom << " atoms updated." << std::endl;
   have_unsaved_changes_flag = 1; 

   make_bonds_type_checked();
}


// helper function for above function
bool
molecule_class_info_t::movable_atom(CAtom *mol_atom, bool replace_coords_with_zero_occ_flag) const {

   bool m = 1;

   if ((mol_atom->occupancy < 0.0001) &&
       (mol_atom->occupancy > -0.0001))
      if (replace_coords_with_zero_occ_flag == 0)
	 m = 0; // zero occupancy and "dont move zero occ atoms is set"
   return m;
}


// This relies on the mol of the asc being different to the mol of the
// atom_sel.
//
// If it is the same, then error and do nothing.
//
// Typically, this is called by fit terminal residue, which has its
// mol created from a pcmmdbmanager() from the molecule of the
// residue, so this is fine in this case.
// 
void
molecule_class_info_t::insert_coords(const atom_selection_container_t &asc) {

   // for each residue in the asc, do a InsResidue into its chain:

   if (! (atom_sel.n_selected_atoms > 0) ) {
      std::cout << "ERROR: Can't insert_coords this asc  - no atoms in molecule!\n"; 
   } else {

      // pointer comparison
      if (asc.mol == atom_sel.mol) {
	 std::cout << "ERROR:: matching asc.mol and atom_sel.mol in insert_coords\n";
	 std::cout << "ERROR:: new algorithm required\n";
      } else {
         make_backup(); // checks backup_this_molecule
	 insert_coords_internal(asc);
      }
   }
} 

void 
molecule_class_info_t::insert_coords_internal(const atom_selection_container_t &asc) {

   // run over each chain, residue of the asc (if terminal residue
   // fit only one chain, one residue, of course).

   short int inserted = 0; // not inserted yet
   CChain *asc_chain;
   int imod = 1; 
   CModel *asc_model_p = asc.mol->GetModel(imod);
   int asc_n_chains = asc_model_p->GetNumberOfChains();
   for (int i_asc_chain=0; i_asc_chain<asc_n_chains; i_asc_chain++) {
      asc_chain = asc.mol->GetChain(1,i_asc_chain);
      int nres_asc = asc_chain->GetNumberOfResidues();
//       std::cout << "DEBUG:: There are " << nres_asc << " residues in "
// 		<< "asc_chain (chain id: " << asc_chain->GetChainID()
// 		<< ")." << std::endl;

      int udd_atom_index = asc.UDDAtomIndexHandle; 

      for (int ires_asc=0; ires_asc<nres_asc; ires_asc++) {
	 CResidue *asc_residue = asc_chain->GetResidue(ires_asc);

	 // Now find the corresponding chain in our atom_sel.mol:

	 CChain *chain;
	 int imodel = 1;
	 int n_chains = atom_sel.mol->GetNumberOfChains(imodel);
	 for (int i_chain=0; i_chain<n_chains; i_chain++) {

	    chain = atom_sel.mol->GetChain(1,i_chain);

	    // test chains
	    std::string asc_chain_str(asc_chain->GetChainID());
	    std::string mol_chain_str(    chain->GetChainID());
// 	    std::cout << "comparig chain ids :" << asc_chain_str << ": :"
// 		      << mol_chain_str << ":" << std::endl;
	    if (asc_chain_str == mol_chain_str) {

	       // insert that residue!
	       PCResidue res = 
		  coot::deep_copy_this_residue(asc_residue, "", 1, udd_atom_index);
// 	       std::cout  << "DEBUG:: inserting residue in chain "
// 			  << mol_chain << " residue number "
// 			  << asc_residue->GetSeqNum()
// 			  << " i_chain = " << i_chain
// 			  << " ires_asc = " << ires_asc
// 			  << std::endl;

	       std::pair<int, CResidue *> serial_number =
		  find_serial_number_for_insert(asc_residue->GetSeqNum(),
						mol_chain_str);

// 	       std::cout << "DEBUG:: returned serial_number: " << serial_number
// 			 << std::endl;

	       if (res) { 

		  if (serial_number.first != -1) {
		     // insert at this position (other residues are
		     // shifted up).
		     chain->InsResidue(res, serial_number.first);
		     coot::copy_segid(serial_number.second, res);
		     inserted = 1;
		  } else { 

		     std::cout << "DEBUG:: insert_coords_internal() add residue\n";
		     CResidue *last_residue = last_residue_in_chain(chain);
		     if (last_residue) { 

			if (0) {   // debug::
			   int nat = last_residue->GetNumberOfAtoms();
			   for (int iat=0; iat<nat; iat++) {
			      std::cout << iat << " of " << nat << " "
					<< last_residue->GetAtom(iat) << std::endl;
			   }
			}
			
			chain->AddResidue(res);
			coot::copy_segid(last_residue, res);
			inserted = 1;
		     }
		  }
	       }
	    }
	    //if (inserted) break;
	 }
      

	 if (! inserted) {
	    // OK, there was no chain in the current mol that matches
	    // the chain of the asc.
	    // Let's copy the asc chain and add it to atom_sel.mol
	    CChain *new_chain = new CChain;
	    int imodel = 1;
	    CModel *this_model = atom_sel.mol->GetModel(imodel);
	    this_model->AddChain(new_chain);
	    new_chain->SetChainID(asc_chain->GetChainID());

	    std::cout << "DEBUG:: Creating a new chain " << asc_chain->GetChainID()
		      << std::endl;

	    CResidue *res = 
	       coot::deep_copy_this_residue(asc_residue, "", 1, udd_atom_index);
	    if (res) { 
	       new_chain->AddResidue(res);
	       atom_sel.mol->FinishStructEdit(); // so that we don't keep adding a
	                                         // new Chain to atom_sel.mol
	    }

	 } 
	 //if (inserted) break;
      }
      //if (inserted) break;
   }
   atom_sel.mol->FinishStructEdit();
   update_molecule_after_additions();
}



void 
molecule_class_info_t::insert_coords_change_altconf(const atom_selection_container_t &asc) {

   // There are 2 things we want to do here.
   //
   // For matching atoms, if there is only one atom that matches the
   // spec (appart from the altconf), we move the original atoms
   // altconf to "A".  If they are not "" we leave don't change the
   // altconf.
   // (see table in // molecule_class_info_t::make_new_alt_conf()
   // [molecule-class-info-other.cc]
   // 

   // The second thing is the change the occ of the existing atoms:
   // 1 atom:  new_occ_existing_atom = 1 - occ_new_atom;
   // general: new_occ_existing_atom = current_occ_existing_atom - occ_new_atom/n_similar_atoms
   // where n_similar_atoms is the number of alt confs we have for that atom.
//    std::cout << "DEBUG:: ----------------- in insert_coords_change_altconf ------ " << std::endl;
//    std::cout << "DEBUG:: IN insert_coords_change_altconf" << std::endl;
   make_backup();
   
   // OK if we were from a shelx ins file, then we have to create a
   // new FVAR for this new alt conf thingy.
   int shelx_occ_fvar_number = -1; 
   if (is_from_shelx_ins_flag) {
      // OK, what was the occupancy?
      if (asc.n_selected_atoms > 0) {
	 float occ = asc.atom_selection[0]->occupancy;
// 	 std::cout << "DEBUG:: IN insert_coords_change_altconf adding fvar "
// 		   << occ << std::endl;
	 shelx_occ_fvar_number = 10 * shelxins.add_fvar(occ); // FVAR 1 is not written
	 // to SHELX file so if shelx ins file has 1 FVAR value, then we've just
	 // created shelx FVAR 3.
	 shelx_occ_fvar_number += 1;  // so thats (e.g.) 1 x the 20th FVAR
      } 
   }
   
   char *chain_id;
   char *atom_name;
   int  resno;
   float occ; 
   CAtom *at;
   for(int i=0; i<asc.n_selected_atoms; i++) { 
      at = asc.atom_selection[i];
      chain_id  = at->GetChainID();
      atom_name = at->GetAtomName();
      resno     = at->GetSeqNum();
      occ       = at->occupancy;
      char *inscode = at->GetInsCode();

      // Now find that atom the corresponding atom (with altconf "" in
      // the original atoms).  We skip over atoms that don't have
      // altconf "").

      int selHnd = atom_sel.mol->NewSelection();
      atom_sel.mol->SelectAtoms(selHnd, 0, chain_id,
				resno, inscode,
				resno, inscode,
				"*", // residue name
				atom_name,
				"*", 
				"*"); // alt-loc
      int nSelAtoms;
      PPCAtom local_SelAtom;
      atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);

      if (nSelAtoms == 0) { 
	 // debugging
// 	       std::cout << "add alt conf: skipping " << atom_name << "/"
// 		<< resno << "/" << chain_id << " in this molecule: ("
// 		<<  imol_no << ") " << name_ << std::endl; 
      } else {

	 // Let's deal with every atom in this residue that that has
	 // the same atom name.

	 // first check for atoms with alt conf "" and move them over if needed
	 for (int iat=0; iat<nSelAtoms; iat++) {
	    std::string current_alt_conf = local_SelAtom[iat]->altLoc;
	    if (current_alt_conf == "") { 
	       std::string new_alt_conf("A");
	       // force it down the atom's throat :)
	       strncpy(local_SelAtom[0]->altLoc, new_alt_conf.c_str(), 2);
	    }
	 }

	 if (shelx_occ_fvar_number == -1) {
	    // i.e. was not from a shelx ins file (normal case):
	    // 
	    // now stuff occupancies in...
	    for (int iat=0; iat<nSelAtoms; iat++) {
	       // local_SelAtom[0]->occupancy = 1.0 - occ; // complemetary (1+1 atom case)
	       local_SelAtom[iat]->occupancy -= occ/float(nSelAtoms); // complemetary (general case)

	       // 20091227 But don't add atoms with negative
	       // occupancy, e.g. the residue before the split was at
	       // zero occupancy. 
	       if (local_SelAtom[iat]->occupancy < 0.0)
		  local_SelAtom[iat]->occupancy = 0.0;
	    }

	 } else {

	    // This is the SHELX case:
	    if (nSelAtoms > 1) {
	       // so we added an alt conf to a residue that already
	       // has an alt conf.  This involves messing with SUMP.
	       // Let's not handle that for now.
	       std::cout << "WARNING:: SHELX occupancy handler under-resourced on handling "
			 << at << std::endl;
	    }  else {
	       local_SelAtom[0]->occupancy = -shelx_occ_fvar_number;
	    }
	 } 
      }
      atom_sel.mol->DeleteSelection(selHnd);
   } 
   insert_coords_atoms_into_residue_internal(asc, shelx_occ_fvar_number);

}

// In this instance, we don't want to install a whole residue, we want
// to install atoms in this residue (alt conf B) into a a atom_sel mol
// residue that contains (say) "" and "A".
//
// -1 is passed as shelx_occ_fvar_number if this atom was not from a
// SHELX ins file.  If shelx_occ_fvar_number > 1, then use this as the
// new atom's occupancy.
// 
void 
molecule_class_info_t::insert_coords_atoms_into_residue_internal(const atom_selection_container_t &asc,
								 int shelx_occ_fvar_number) {

   char *chain_id;
   char *atom_name;
   int  resno;
   CAtom *at;
   int afix_handle_this_mol = -1;
   int afix_handle_intermediate_mol = -1;

   afix_handle_this_mol    = atom_sel.mol->GetUDDHandle(UDR_ATOM, "shelx afix");
   afix_handle_intermediate_mol = asc.mol->GetUDDHandle(UDR_ATOM, "shelx afix");

//    std::cout << "DEBUG in insert_coords_atoms_into_residue_internal afix handles:"
// 	     << afix_handle_this_mol << " " << afix_handle_intermediate_mol << std::endl;
   
   for(int i=0; i<asc.n_selected_atoms; i++) { 
      at = asc.atom_selection[i];
      chain_id  = at->GetChainID();
      atom_name = at->GetAtomName();
      resno     = at->GetSeqNum();

      // Now find the corresponding residue in atom_sel.mol;

      int selHnd = atom_sel.mol->NewSelection();
      int nSelResidues;
      PPCResidue SelResidues;
      atom_sel.mol->Select(selHnd, STYPE_RESIDUE, 1,
			   chain_id, 
			   resno, "*",
			   resno, "*",
			   "*",  // residue name
			   "*",  // Residue must contain this atom name?
			   "*",  // Residue must contain this Element?
			   "*",  // altLocs
			   SKEY_NEW // selection key
			   );
      atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
      
      if (nSelResidues != 1) { 
	 std::cout << "ERROR:: something broken in residue selection in ";
	 std::cout << "insert_coords_atoms_into_residue_internal: got " << nSelResidues
		   << " residues." << std::endl;
      } else { 
	 CAtom *t = new CAtom;
	 t->Copy(at);
	 SelResidues[0]->AddAtom(t);
	 // if these coords were from a shelx ins file, then we are
	 // passed the free varible number (FVAR) for this atom's
	 // occupancy.
	 if (shelx_occ_fvar_number > 1)
	    t->occupancy = shelx_occ_fvar_number;

	 // If these coords were from a shelx ins file, then we need
	 // to copy the AFIX numbers too:
	 int afix_number; // set by getUDD
	 if (afix_handle_intermediate_mol > -1) {
	    int ierr = at->GetUDData(afix_handle_intermediate_mol, afix_number);
	    if (ierr == UDDATA_Ok) {
	       if (afix_handle_this_mol > -1) {
		  t->PutUDData(afix_handle_this_mol, afix_number);
	       } else {
		  std::cout << "ERROR:: bad afix handle for this molecule in "
			    << "insert_coords_atoms_into_residue_internal"
			    << afix_handle_this_mol << " " << at << std::endl;
	       }

 	    } else {
	       if (is_from_shelx_ins_flag) 
		  std::cout << "ERROR:: attempt to get UDD afix number from "
			    << "intermediate molecule failed " << at << std::endl;
 	    } 
	 } else {
	    std::cout << "ERROR:: bad afix handle for intermediate molecule in "
		      << "insert_coords_atoms_into_residue_internal"
		      << afix_handle_intermediate_mol << " " << at << std::endl;
	 } 
      } 
      atom_sel.mol->DeleteSelection(selHnd);
   } 
   atom_sel.mol->FinishStructEdit();
   update_molecule_after_additions();
} 

// We need to find the serial number of the residue after the residue
// we want to insert (i.e. the new residue will be inserted just
// before the residue whose serial number we return).
//
// return -1 on error.
std::pair<int, CResidue *>
molecule_class_info_t::find_serial_number_for_insert(int seqnum_new,
						     const std::string &chain_id) const {

   int iserial_no = -1;
   int current_diff = 999999;
   CChain *chain;
   int n_chains = atom_sel.mol->GetNumberOfChains(1);
   CResidue *res = NULL;
   
   for (int i_chain=0; i_chain<n_chains; i_chain++) {
      
      chain = atom_sel.mol->GetChain(1,i_chain);
      
      // test chains
      std::string mol_chain(    chain->GetChainID());
      if (chain_id == mol_chain) {

	 // find the 
	 int nres = chain->GetNumberOfResidues();
	 
	 for (int ires=0; ires<nres; ires++) { // ires is a serial number
	    CResidue *residue = chain->GetResidue(ires);

	    // we are looking for the smallest negative diff:
	    // 
	    int diff = residue->GetSeqNum() - seqnum_new;

	    if ( (diff > 0) && (diff < current_diff) ) {
	       iserial_no = ires;
	       res = residue;
	       current_diff = diff;
	    }
	 }
      }
   }
   return std::pair<int, CResidue *> (iserial_no, res);
} 


// Put the regularization results back into the molecule:
//
// Recall that regularized_asc contains an atom_selection_container_t
// with the new coordinates in.  the mol contains all the molecule and
// the atom_selection contains just the moving parts.
//
void
molecule_class_info_t::add_coords(const atom_selection_container_t &asc) {

   CAtom *atom;
   CAtom *mol_atom; // an atom already existing in mol
   int n_atom = 0;
   CChain *chain;

   // std::cout << "DEBUG:: ----------------- in add_coords ----------- " << std::endl;

   make_backup();
   
   for (int i=0; i<asc.n_selected_atoms; i++) {
      int idone = 0;
      atom = asc.atom_selection[i];
      // chain = atom->GetChain();
      
      // run over chains of the existing mol
      int nchains = atom_sel.mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) { 
	 
	 chain = atom_sel.mol->GetChain(1,ichain);
	 std::string atom_chain_id(atom->GetChainID());
	 std::string mol_chain_id(chain->GetChainID());
	 
	 if (atom_chain_id == mol_chain_id) {

	    // int iseqno_at = atom->GetSeqNum();
	    int nres = chain->GetNumberOfResidues();
	    for (int ires=0; ires<nres; ires++) { 
	       PCResidue res = chain->GetResidue(ires);
	       if (res) { // fixes bug 030813, but best way?
		          // 030814, ah, we discover FinishStructEdit().
		  if (res->GetSeqNum() == atom->GetSeqNum()) { 
		     int natom = res->GetNumberOfAtoms();
		     for (int iat=0; iat<natom; iat++) {
			
			mol_atom = res->GetAtom(atom->GetAtomName());
			if (mol_atom) { // we got a match
			   // replace the coordinates then
			   
			   // This should not happen very often
			   std::cout << "add_coords: replacing " << mol_atom 
				     << " with new atom " << atom << std::endl;
			   mol_atom->SetCoordinates(atom->x,
						    atom->y,
						    atom->z,
						    mol_atom->occupancy,
						    mol_atom->tempFactor);
			   idone = 1;
			   break;
			   
			} else { 
			   
			   std::cout << "adding atom to existing residue " 
				     << atom << " (already has " 
				     << res->GetNumberOfAtoms() << " atoms)" 
				     << std::endl;
			   CAtom *new_atom = new CAtom;
			   new_atom->Copy(atom);
			   res->AddAtom(new_atom);
			   new_atom->occupancy = 1.0;
			   new_atom->tempFactor = 10.0;
			   // chain id:
			   std::cout << "setting chainid of this new atom from "
				     << new_atom->GetChainID() << " to : "
				     << atom->GetChainID() << std::endl;
			   new_atom->residue->chain->SetChainID(atom->GetChainID());
			   idone = 1;
			   n_atom++;
			   break;
			}
		     }
		  }
	       }
	       if (idone == 1) break;
	    } // residue loop
	 }
      }

      if (idone == 0) { 

	 std::cout << "adding whole residue triggered by atom " 
		   << atom << std::endl;
	 std::cout << "     with element " << atom->element << std::endl;

	 // in this bit of code, atom is an atom from the asc and
	 // atom_p is a new atom that we are adding to a new residue
	 // (that we are adding to an existing chain (that we do a
	 // lookup to find)).
	 // 
	 CResidue *res_p = new CResidue;
	 CAtom *atom_p = new CAtom;
	 // CChain *chain_p = atom_sel.mol->GetChain(1,0);
	 PCChain chain_p = atom_sel.mol->GetChain(1,atom->GetChainID());
	 chain_p->AddResidue(res_p);
	 atom_p->SetAtomName(atom->name);
	 atom_p->SetCoordinates(atom->x, atom->y, atom->z,
				atom->occupancy, atom->tempFactor);
	 atom_p->SetElementName(atom->element);
	 res_p->AddAtom(atom_p);
	 res_p->seqNum = atom->GetSeqNum();
	 res_p->SetResID(atom->residue->name,
			 atom->GetSeqNum(),
			 atom->GetInsCode());
	 
	 // add to end:
	 atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
	 atom_sel.mol->FinishStructEdit();
      }
   }

   atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
   atom_sel.mol->FinishStructEdit();
   std::cout << "INFO:: " << n_atom << " atoms added to molecule." << std::endl;

   // now regenerate the atom_selection
   //

   // Uncomment when we have fixed the bug, this seems not to be it.
   // 
   // clear out the old
   // if (atom_sel.atom_selection != NULL)
   // delete [] atom_sel.atom_selection;
   
   // and in with the new:
   int selHnd = atom_sel.mol->NewSelection(); 
   atom_sel.mol->SelectAtoms(selHnd, 0,"*",ANY_RES,"*",ANY_RES,
			    "*","*",  // EndInsertionCode, RNames
			    "*","*",  // ANames, Elements
			    "*" );    // Alternate locations.

   int old_n_atoms = atom_sel.n_selected_atoms;
   atom_sel.mol->GetSelIndex(selHnd, 
			     atom_sel.atom_selection, 
			     atom_sel.n_selected_atoms);

   std::cout << "INFO:: old n_atoms: " << old_n_atoms << " new: " 
	     << atom_sel.n_selected_atoms << std::endl;

   debug_selection();
   have_unsaved_changes_flag = 1; 

   make_bonds_type_checked();
   // std::cout << "DEBUG:: ---------------- done add_coords ----------- " << std::endl;
} 



void
molecule_class_info_t::close_yourself() {

   // Deletion causing problems on application closure

   short int was_map = 0;
   short int was_coords = 0;

   name_ = ""; // not "Baton Atoms" or anything now.
   
   if (atom_sel.n_selected_atoms > 0)
      was_coords = 1;

   if (xmap_is_filled[0])
      was_map = 1;

   // delete from display manager combo box
   // 
   graphics_info_t g;
   GtkWidget *display_control_window = g.display_control_window();
   // 
   if (display_control_window) { // is being displayed
      std::string display_frame_name = "display_mol_frame_";
      if (was_map)
	 display_frame_name = "display_map_frame_";
      display_frame_name += g.int_to_string(imol_no);
      // std::cout << "DEBUG:: looking up " << display_frame_name << std::endl;
      GtkWidget *display_frame = lookup_widget(display_control_window,
					       display_frame_name.c_str());
      if (display_frame) {

	 gtk_widget_destroy(display_frame);

      }
   } else {
      // std::cout << "close: display_control_window is not active" << std::endl;
   }

   if (was_coords) { 
      atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);
      delete atom_sel.mol;
      // atom_sel.mol = 0; done later
   }

   if (was_map) {
      fc_skeleton_draw_on = 0; // turn off the skeleton
      delete [] xmap_list;
      xmap_list = NULL;
      xmap_is_filled[0] = 0;
      max_xmaps = 0;
   }

   bonds_box.clear_up();
   // symmetry_bonds_box?  (It is a vector of pairs)
   drawit = 0;
   drawit_for_map = 0;

   // Do these whatever the molecule type:
   atom_sel.n_selected_atoms = 0;
   atom_sel.atom_selection = NULL;
   atom_sel.mol = NULL;

   //
   // gl widget redraw is done in close_molecule
}


// Return the atom index of the "next" atom
// -1 on failure.
int
molecule_class_info_t::intelligent_next_atom(const std::string &chain_id,
					     int resno,
					     const std::string &atom_name,
					     const std::string &ins_code) {

   // This is really a problem of "what is the next residue?", the
   // actual atom is a superficial problem that is handled by
   // intelligent_this_residue_atom().  We simply have to find this
   // residue in the chain, and return the residue after that.
   //
   // If there is no next residue, use the residue at the beginning of
   // the next chain.
   //
   // If there is no next chain, use the residue at the start of the
   // first chain.
   //
   // If this residue can't be found, then go through chain and look
   // for the first residue that has higher residue number than resno.

   int i_atom_index = -1; // failure initially.
   if (atom_sel.n_selected_atoms <= 0 || atom_sel.mol == NULL) { 
      std::cout << "ERROR:: trying to move to (next) atom of a closed molecule!\n";
   } else {

      CResidue *first_residue = NULL;
      CResidue *next_residue = NULL;
      bool found_this_residue = 0; 
      int imod = 1;
      CModel *model_p = atom_sel.mol->GetModel(imod);
      CChain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 if (chain_id == chain_p->GetChainID() || (found_this_residue && !next_residue)) {
	    int nres = chain_p->GetNumberOfResidues();
	    PCResidue residue_p;
	    for (int ires=0; ires<nres; ires++) { 
	       residue_p = chain_p->GetResidue(ires);
	       if (! first_residue)
		  first_residue = residue_p;
	       if (found_this_residue) { 
		  next_residue = residue_p;
		  break;
	       }
	       if (residue_p->GetSeqNum() == resno) {
		  if (ins_code == residue_p->GetInsCode()) {
		     found_this_residue = 1; // setup for next loop
		  }
	       }
	    }
	 }
	 if (next_residue)
	    break;
      }

      // Now the case where this residue was not in the molecule
      // (e.g. we just deleted "this" water)
      //
      if (!next_residue) { 
	 for (int ichain=0; ichain<nchains; ichain++) {
	    chain_p = model_p->GetChain(ichain);
	    if (chain_id == chain_p->GetChainID()) { 
	       int nres = chain_p->GetNumberOfResidues();
	       PCResidue residue_p;
	       for (int ires=0; ires<nres; ires++) { 
		  residue_p = chain_p->GetResidue(ires);
		  if (residue_p->GetSeqNum() > resno) {
		     next_residue = residue_p;
		     break;
		  }
	       }
	    }
	    if (next_residue)
	       break;
	 }
      }
	       
      // The case where this residue was the last in the molecule:
      // 
      if (found_this_residue && !next_residue) {
	 int nchains = model_p->GetNumberOfChains();
	 if (nchains > 0) { 
	    chain_p = model_p->GetChain(0);
	    int nres = chain_p->GetNumberOfResidues();
	    if (nres > 0) 
	       next_residue = chain_p->GetResidue(0);
	 }
      }

      if (next_residue)
	 i_atom_index = intelligent_this_residue_atom(next_residue);
      
   }
//    std::cout << "DEBUG:: Intelligent:: returning index " << i_atom_index << std::endl;
//    CAtom *at = atom_sel.atom_selection[i_atom_index];
//    std::cout << "        Returing index of " << at << std::endl;
   return i_atom_index;
} 


// Return the atom index of the "next" atom
// -1 on failure.
int
molecule_class_info_t::intelligent_previous_atom(const std::string &chain_id,
						 int resno,
						 const std::string &atom_name,
						 const std::string &ins_code) {

   // This is quite similar to intelligent_next_atom() (see comments
   // there).  However, this is a bit more complex, because we keep a
   // hold on the previous residue and only if "this" residue is found
   // do we make the previous residue previous_residue.
   //
   // We don't handle the "last_residue" in molecule case, the analog
   // of first_residue in the intelligent_next_atom() function.
   // (Maybe we should?)
   
   int i_atom_index = -1;
   if (atom_sel.n_selected_atoms <= 0 || atom_sel.mol == NULL) { 
      std::cout << "ERROR:: trying to move to (next) atom of a closed molecule!\n";
   } else {

      CResidue *prev_residue_candidate = NULL;
      CResidue *prev_residue = NULL;
      bool found_this_residue = 0; 
      int imod = 1;
      CModel *model_p = atom_sel.mol->GetModel(imod);
      CChain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 if (chain_id == chain_p->GetChainID() || (found_this_residue && !prev_residue)) {
	    int nres = chain_p->GetNumberOfResidues();
	    PCResidue residue_p;
	    for (int ires=0; ires<nres; ires++) { 
	       residue_p = chain_p->GetResidue(ires);
	       if (residue_p->GetSeqNum() == resno) {
		  if (ins_code == residue_p->GetInsCode()) {
		     found_this_residue = 1;
		     if (prev_residue_candidate) { 
			prev_residue = prev_residue_candidate;
			break;
		     }
		  }
	       }
	       if (prev_residue)
		  break;
	       prev_residue_candidate = residue_p;
	    }
	 }
	 if (prev_residue)
	    break;
      }

      // Handle the case where we are on first atom of water chain and
      // want to go back to protein (in that case
      // prev_residue_candidate would not have been set because it
      // never passes the chain id test)
      // 
      if (! prev_residue) {
	 CChain *prev_chain = NULL;
	 CChain *prev_chain_candidate = NULL;
	 int nchains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<nchains; ichain++) {
	    chain_p = model_p->GetChain(ichain);
	    if (chain_id == chain_p->GetChainID()) {
	       if (prev_chain_candidate) { 
		  prev_chain = prev_chain_candidate;
		  break;
	       }
	    }
	    prev_chain_candidate = chain_p; // for next loop
	 }
	 if (prev_chain) {
	    // the last residue in prev_chain then:
	    int nres = prev_chain->GetNumberOfResidues();
	    if (nres > 0) 
	       prev_residue = prev_chain->GetResidue(nres-1);
	 } 
      }
      
      if (prev_residue)
	 i_atom_index = intelligent_this_residue_atom(prev_residue);
   }
   return i_atom_index;
}


// If there is a CA in this residue then return the index of that
// atom, if not, then return the index of the first atom in the
// residue.
// 
// Return -1 on no atoms in residue.
//
int
molecule_class_info_t::intelligent_this_residue_atom(CResidue *res_p) const {

   PPCAtom residue_atoms;
   int nResidueAtoms;
   int ir = -1; 
   
   res_p->GetAtomTable(residue_atoms, nResidueAtoms);
   for (int i=0; i<nResidueAtoms; i++) {
      std::string atom_name(residue_atoms[i]->name);
      if (atom_name == " CA ") {
	 ir = atom_to_atom_index(residue_atoms[i]);
	 if (ir == -1)
	    ir = full_atom_spec_to_atom_index(residue_atoms[i]->GetChainID(),
					      residue_atoms[i]->GetSeqNum(),
					      residue_atoms[i]->GetInsCode(),
					      residue_atoms[i]->name,
					      residue_atoms[i]->altLoc);
      }
   }

   if (ir == -1) { 
      if (nResidueAtoms > 0) {
	 ir = atom_to_atom_index(residue_atoms[0]);
	 if (ir == -1)
	    ir =  full_atom_spec_to_atom_index(residue_atoms[0]->GetChainID(),
					       residue_atoms[0]->GetSeqNum(),
					       residue_atoms[0]->GetInsCode(),
					       residue_atoms[0]->name,
					       residue_atoms[0]->altLoc);
      }
   }
   return ir;
}

// If there is a CA in this residue then return that atom (pointer)
// atom, if not, then return the index of the first atom in the
// residue.
// 
// Return NULL on no atoms in residue.
//
CAtom *
molecule_class_info_t::intelligent_this_residue_mmdb_atom(CResidue *res_p) const {
   PPCAtom residue_atoms;
   int nResidueAtoms;
   
   res_p->GetAtomTable(residue_atoms, nResidueAtoms);
   for (int i=0; i<nResidueAtoms; i++) {
      std::string atom_name(residue_atoms[i]->name);
      if (atom_name == " CA ") {
	 return residue_atoms[i];
      }
   }

   if (nResidueAtoms > 0) {
      return residue_atoms[0];
   }

   // failure
   return NULL;

}

// Return pointer to atom " CA ", or the first atom in the residue, or
// null (no residue or atoms error):
// 
CAtom *
molecule_class_info_t::atom_intelligent(const std::string &chain_id, int resno,
					const std::string &ins_code) const { 

   CAtom *at = NULL; 

   if (atom_sel.n_selected_atoms > 0) { 
      int selHnd = atom_sel.mol->NewSelection();
      PPCResidue SelResidue;
      int nSelResidues;

      char *ins_code_search = (char *) ins_code.c_str(); // bleugh (as usual)
      
      atom_sel.mol->Select (selHnd, STYPE_RESIDUE, 0,
			    (char *)chain_id.c_str(), 
			    resno, ins_code_search,
			    resno, ins_code_search,
			    "*",  // residue name
			    "*",  // Residue must contain this atom name?
			    "*",  // Residue must contain this Element?
			    "*",  // altLocs
			    SKEY_NEW // selection key
			    );

      atom_sel.mol->GetSelIndex(selHnd, SelResidue, nSelResidues);
      
      if (nSelResidues == 0) {
	 std::cout << "INFO:: No selected residues" << std::endl;
      } else {
	 
	 PPCAtom residue_atoms;
	 int nResidueAtoms;
	 SelResidue[0]->GetAtomTable(residue_atoms, nResidueAtoms);
	 if (nResidueAtoms == 0) {
	    std::cout << "INFO:: No atoms in residue" << std::endl;
	 } else {
	    short int found_it = 0;
	    for (int i=0; i<nResidueAtoms; i++) { 
	       if (std::string(residue_atoms[i]->name) == std::string(" CA ")) { 
		  at = residue_atoms[i];
		  found_it = 1;
		  break;
	       }
	    } 
	    if (! found_it) 
	       at = residue_atoms[0];
	 }
      }
      atom_sel.mol->DeleteSelection(selHnd); // Safe to put it here?
   }
   return at;
} 


// ----------------------------------------------------------------------
//               Pointer Atoms
// ----------------------------------------------------------------------
void
molecule_class_info_t::add_pointer_atom(coot::Cartesian pos) {


   if (atom_sel.mol) { 
      CChain *chain_p = water_chain();
      
      if (! chain_p) {
	 // we have to make one then
	 chain_p = new CChain;
	 std::pair<short int, std::string> p = unused_chain_id();
	 if (p.first)
	    chain_p->SetChainID(p.second.c_str());
	 CModel *model_p = atom_sel.mol->GetModel(1);
	 model_p->AddChain(chain_p);
      }
   
      make_backup();
      std::string mol_chain_id(chain_p->GetChainID());
      // int ires_prev = chain_p->GetNumberOfResidues();
      int ires_prev = coot::util::max_resno_in_chain(chain_p).second;

      CResidue *res_p = new CResidue;
      CAtom *atom_p = new CAtom;
      chain_p->AddResidue(res_p);
      atom_p->SetAtomName(" O  ");
      float bf = graphics_info_t::default_new_atoms_b_factor;
      atom_p->SetCoordinates(pos.x(), pos.y(), pos.z(), 1.0, bf);

      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
      res_p->seqNum = ires_prev + 1;
      res_p->SetResName("HOH");
      coot::hetify_residue_atoms(res_p);

      atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      std::cout << atom_p << " added to molecule" << std::endl;

      have_unsaved_changes_flag = 1; 
      make_bonds_type_checked();
   }
}

// This is a bit messy, I'm afraid - we test single atom twice. If you use this for 
// a multiatom other than SO4 and P04, you will need to add it to the type test
// 
void 
molecule_class_info_t::add_typed_pointer_atom(coot::Cartesian pos, const std::string &type) { 

   short int single_atom = 1; // true

   std::cout << "INFO:: adding atom of type " << type << " at " << pos << std::endl;
   make_backup();

   // we get a chain pointer or NULL, if there is not yet a chain only
   // of the given type:
   CChain *single_type = coot::util::chain_only_of_type(atom_sel.mol, type);

   // We do different things (e.g adding the chain) if this is a new
   // chain or a pre-existing one, let's set a flag.
   short int pre_existing_chain_flag;
   CChain *chain_p;
   if (single_type) {
      chain_p = single_type;
      pre_existing_chain_flag = 1;
   } else {
      chain_p = new CChain;
      pre_existing_chain_flag = 0;
   }
   
   std::pair<short int, std::string> mol_chain_id = unused_chain_id();
   CResidue *res_p = new CResidue;

   // type test
   if (type == "PO4") single_atom = 0;
   if (type == "SO4") single_atom = 0;

   if (single_atom) { 
      CAtom *atom_p = new CAtom;
      float occ;
      if (is_from_shelx_ins_flag)
	 occ = 11.0;
      else
	 occ = 1.0;
      atom_p->SetCoordinates(pos.x(), pos.y(), pos.z(), occ,
			     graphics_info_t::default_new_atoms_b_factor);
      atom_p->Het = 1; // it's a HETATM.

      if (type == "Water") { 

	 // special rule for water: we add a water to a water chain if
	 // possible

	 atom_p->SetAtomName(" O  ");
	 atom_p->SetElementName(" O");
	 res_p->SetResName("HOH");

	 CChain *w = water_chain();
	 int wresno = 1;
	 
	 if (w) {
	    // remove a TER atom if it exists on the last residue
	    // prior to insertion of a new residue.
	    remove_TER_on_last_residue(w);
	    
	    // Now add atom to chain w.
	    std::pair<short int, int> wresno_pair = next_residue_in_chain(w);
	    if (wresno_pair.first) { 
	       wresno = wresno_pair.second;
	    } else { 
	       wresno = 1;
	    } 
	    res_p->seqNum = wresno;
	    res_p->AddAtom(atom_p);
	    w->AddResidue(res_p);
	    std::cout << atom_p << " added to molecule" << std::endl;
	    atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
	    atom_sel.mol->FinishStructEdit();
	    atom_sel = make_asc(atom_sel.mol);
	    have_unsaved_changes_flag = 1;
	    make_bonds_type_checked();
	    
	 } else { 
	    // There was no water chain
	    res_p->AddAtom(atom_p);
	    std::cout << atom_p << " added to molecule (and new chain)" << std::endl;
	    if (!pre_existing_chain_flag) { 
	       chain_p->SetChainID(mol_chain_id.second.c_str());
	       atom_sel.mol->GetModel(1)->AddChain(chain_p);
	    }
	    res_p->seqNum = 1; // start of a new chain.
	    chain_p->AddResidue(res_p);
	    atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
	    atom_sel.mol->FinishStructEdit();
	    atom_sel = make_asc(atom_sel.mol);
	    have_unsaved_changes_flag = 1;
	    make_bonds_type_checked();
	 } 
      } else { 
	
  	 // Not water
	 std::string element = "";

	 if (mol_chain_id.first || pre_existing_chain_flag) { 

	    if (type == "Br") { 
	       atom_p->SetAtomName("BR  ");
	       atom_p->SetElementName("BR");
	       res_p->SetResName(" BR");
	       element = "BR";
	    } else { 
	       if (type == "Ca") { 
		  atom_p->SetAtomName("CA  ");
		  atom_p->SetElementName("CA");
		  res_p->SetResName(" CA");
		  element = "CA";
	       } else { 
		  if (type == "Na") { 
		     atom_p->SetAtomName("NA  ");
		     atom_p->SetElementName("NA");
		     res_p->SetResName(" NA");
		     element = "NA";
		  } else { 
		     if (type == "Cl") { 
			atom_p->SetAtomName("CL  ");
			atom_p->SetElementName("CL");
			res_p->SetResName(" CL");
			element = "CL";
		     } else { 
			if (type == "Mg") { 
			   atom_p->SetAtomName("MG  ");
			   atom_p->SetElementName("MG");
			   res_p->SetResName(" MG");
			   element = "MG";
			} else { 

			   // User Typed atom:

			   // make up (guess) the residue type and element
			   std::string at_name = type;
			   std::string ele = type;
			   std::string resname = type;
			   if (type.length() > 4)
			      at_name = at_name.substr(0,4);
			   if (type.length() > 3)
			      resname = at_name.substr(0,3);
			   if (type.length() > 2)
			      ele = at_name.substr(0,2);

			   element = ele;
			   atom_p->SetAtomName(at_name.c_str());
			   atom_p->SetElementName(ele.c_str());
			   res_p->SetResName(resname.c_str());
			} 
		     }
		  }
	       }
	    }

	    res_p->AddAtom(atom_p);
	    std::cout << atom_p << " added to molecule" << std::endl;
	    if (! pre_existing_chain_flag) { 
	       chain_p->SetChainID(mol_chain_id.second.c_str());
	       atom_sel.mol->GetModel(1)->AddChain(chain_p);
	    }
	    std::pair<short int, int> ires_prev_pair = coot::util::max_resno_in_chain(chain_p);
	    int previous_max = 0;
	    if (ires_prev_pair.first) { // was not an empty chain
	       previous_max =  ires_prev_pair.second;
	       res_p->seqNum = previous_max + 1;
	    } else {

	       // was an empty chain.  Handle the shelx case:

	       if (! is_from_shelx_ins_flag) { 
		  res_p->seqNum = 1 ; // start of a new chain.
	       } else {
		  // in a shelx molecule, we can't make the residue
		  // number 1 because there are no chains.  We need to
		  // make the residue number bigger than the biggest
		  // residue number so far.
		  std::pair<short int, int> ires_prev_pair =
		     coot::util::max_resno_in_molecule(atom_sel.mol);
		  if (ires_prev_pair.first) {
		     res_p->seqNum = ires_prev_pair.second + 1;
		  } else {
		     res_p->seqNum = 1;
		  }
	       }
	       
	    }

	    // Add this element to the sfac (redundancy check in the addition function
	    if (is_from_shelx_ins_flag) {
	       shelxins.add_sfac(element);
	    } 
	    chain_p->AddResidue(res_p);
	    atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
	    atom_sel.mol->FinishStructEdit();
	    atom_sel = make_asc(atom_sel.mol);
	    have_unsaved_changes_flag = 1;
	    make_bonds_type_checked();
	 } else { 
	    std::cout << "WARNING:: Can't find new chain for new atom\n";
	 } 
      } // type was water, or not

   } else { 
      // multi atom:

      if (mol_chain_id.first || pre_existing_chain_flag) { 
	 add_pointer_multiatom(res_p, pos, type);
	 coot::hetify_residue_atoms(res_p);
	 if (! pre_existing_chain_flag) { 
	    chain_p->SetChainID(mol_chain_id.second.c_str());
	    atom_sel.mol->GetModel(1)->AddChain(chain_p);
	 }
	 std::pair<short int, int> ires_prev_pair = coot::util::max_resno_in_chain(chain_p);
	 int previous_max = 0;
	 if (ires_prev_pair.first) { // was not an empty chain
	    previous_max =  ires_prev_pair.second;
	 } 
	 res_p->seqNum = previous_max + 1;
	 
	 chain_p->AddResidue(res_p);
	 atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
	 atom_sel.mol->FinishStructEdit();
	 atom_sel = make_asc(atom_sel.mol);
	 have_unsaved_changes_flag = 1;
	 make_bonds_type_checked();
      } else { 
	 std::cout << "WARNING:: Can't find new chain for new atom\n";
      } 
   }
   // or we could just use update_molecule_after_additions() there.
}

// return status [1 means "usable"] and a chain id [status = 0 when
// there are 2*26 chains...]
// 
std::pair<short int, std::string>
molecule_class_info_t::unused_chain_id() const { 
   
   std::string r("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
   std::pair<short int, std::string> s(0,""); 
   CChain *chain_p;
   if (atom_sel.n_selected_atoms > 0) { 
      CModel *model_p = atom_sel.mol->GetModel(1);
      int nchains = model_p->GetNumberOfChains();
      int idx;
   
      for (int ich=0; ich<nchains; ich++) {
	 chain_p = model_p->GetChain(ich);
	 idx = r.find(std::string(chain_p->GetChainID()));
	 // 	 while (idx != std::string::npos) { 
	    r = r.substr(0, idx) + r.substr(idx+1);
	    idx = r.find(std::string(chain_p->GetChainID()));
	    s.first = 1;
	    // 	 } // Take out the while, as per Ezra's suggestion.
      }
      std::string tstring = r.substr(0,1);
      s.second = tstring;
   } else {
      s.first = 1;
      s.second = "A";
   } 
   return s;
} 

void 
molecule_class_info_t::add_pointer_multiatom(CResidue *res_p, 
					     const coot::Cartesian &pos, const std::string &type) {
   
   coot::Cartesian p;
   float bf = graphics_info_t::default_new_atoms_b_factor;
   res_p->SetResName(type.c_str());
   if (type == "SO4") { 
      CAtom *atom_p;

      atom_p = new CAtom;
      p = pos + coot::Cartesian(0.000, 0.000, 0.088);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" S  ");
      atom_p->SetElementName(" S");
      res_p->AddAtom(atom_p);
      
      atom_p = new CAtom;
      p = pos + coot::Cartesian( 1.227, 0.000, -0.813);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O1 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
      
      atom_p = new CAtom;
      p = pos + coot::Cartesian(-1.227, 0.000, -0.813);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O2 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
      
      atom_p = new CAtom;
      p = pos + coot::Cartesian(  0.000, -1.263, 0.740);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O3 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
      
      atom_p = new CAtom;
      p = pos + coot::Cartesian( 0.000, 1.263, 0.740);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O4 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
   } 

   if (type == "PO4") { 
      CAtom *atom_p;

      atom_p = new CAtom;
      p = pos + coot::Cartesian(0.000, 0.021, 0.036);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" P  ");
      atom_p->SetElementName(" P");
      res_p->AddAtom(atom_p);
      
      atom_p = new CAtom;
      p = pos + coot::Cartesian( 1.315,   0.599,  -0.691 );
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O1 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
      
      atom_p = new CAtom;
      p = pos + coot::Cartesian( -1.315,   0.599,  -0.691);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O2 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
      
      atom_p = new CAtom;
      p = pos + coot::Cartesian( 0.000,  -1.587,  -0.055);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O3 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
      
      atom_p = new CAtom;
      p = pos + coot::Cartesian( 0.000,   0.434,   1.457);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O4 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
   } 

} 


// ----------------------------------------------------------------------
//                 Save yourself
// ----------------------------------------------------------------------
//
// return 0 on success.
int
molecule_class_info_t::save_coordinates(const std::string filename) {

   int ierr = 0;
   std::string ext = coot::util::file_name_extension(filename);
   if (coot::util::extension_is_for_shelx_coords(ext)) {
      write_shelx_ins_file(filename);
   } else {
      byte bz = GZM_NONE;
      ierr = write_atom_selection_file(atom_sel, filename, bz);
   }

   if (ierr) {
      std::cout << "WARNING!! Coordinates write to " << filename
		<< " failed!" << std::endl;
      std::string ws = "WARNING:: export coords: There was an error ";
      ws += "in writing ";
      ws += filename;
      GtkWidget *w = graphics_info_t::wrapped_nothing_bad_dialog(ws);
      gtk_widget_show(w);
   } else {
      have_unsaved_changes_flag = 0;

      // Now we have updated the molecule name, how shall we restore
      // this from the state file?
      std::vector<std::string> strings;
      strings.push_back("handle-read-draw-molecule");
      strings.push_back(single_quote(coot::util::intelligent_debackslash(filename)));
      save_state_command_strings_ = strings;

      name_ = filename;  // hmm... // update go to atom widget now? FIXME.
      std::string::size_type icoot = filename.rfind("-coot-");
      if (icoot != std::string::npos) { 
	 coot_save_index++;
      }
      update_mol_in_display_control_widget();  // FIXME. 
   }
   return ierr;
} 



// Return 1 on yes, unsaved changes present,
//        0 on no
int
molecule_class_info_t::Have_unsaved_changes_p() const {
   if (has_model())
      return have_unsaved_changes_flag;
   else
      return 0;
}


// ----------------------------------------------------------------------
//               Baton Atoms
// ----------------------------------------------------------------------


// Recall that the chain is set by the creation of the empty molecule in 
// graphics_info_t::baton_build_atoms_molecule();
// direction_flag is +1 for forward building, -1 for backwards direction.
//
// In the case of direction_flag being negative, I think that residues
// (atoms) will go into the chain with their seqNums in decreasing
// order.  The may need to be sorted later, I think (maybe not).
// 
CAtom *
molecule_class_info_t::add_baton_atom(coot::Cartesian pos, 
				      int istart_resno,
				      short int iresno_active,
				      short int direction_flag) {

   int nchains = atom_sel.mol->GetNumberOfChains(1);

   if (nchains != 1) {
      std::cout << "failed to add baton atom" << std::endl;
      return NULL;
   }
      
   make_backup();
   CChain *chain_p = atom_sel.mol->GetChain(1,0);
   
   std::string mol_chain_id(chain_p->GetChainID());
   int n_res = chain_p->GetNumberOfResidues();


   // if this is the first atom to be added in the chain, we get the
   // seqnum from the passed istart_resno otherwise it is the seqnum
   // of the previous residue, plus or minus one.
   // 
   int this_res_seqnum; 
   if (n_res == 0) {
      this_res_seqnum = istart_resno;
   } else { 

      if (iresno_active == 0) { 
	 int ires_prev = chain_p->GetResidue(n_res-1)->seqNum; // seqnum of the last
                                                               // residue in chain.
	 this_res_seqnum = ires_prev + 1*direction_flag;
      } else { 
	 this_res_seqnum = istart_resno;
      } 
   }

   CResidue *res_p = new CResidue;
   CAtom *atom_p = new CAtom;
   chain_p->AddResidue(res_p);
   atom_p->SetAtomName(" CA ");
   atom_p->SetCoordinates(pos.x(), pos.y(), pos.z(), 1.0,
			  graphics_info_t::default_new_atoms_b_factor);

   atom_p->SetElementName(" C");
   res_p->AddAtom(atom_p);
   res_p->seqNum = this_res_seqnum;
   res_p->SetResName("ALA");

   atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
   atom_sel.mol->FinishStructEdit();
   atom_sel = make_asc(atom_sel.mol);
   std::cout << atom_p << " added to molecule" << std::endl;

   have_unsaved_changes_flag = 1; 
   make_ca_bonds(2.4, 4.7);

   return atom_p;
}

// Return a vector of upto 3 positions of the most latestly added
// atoms with the most lastest atom addition (that is the passed atom)
// in the back() slot of the vector.
// 
std::vector<clipper::Coord_orth>
molecule_class_info_t::previous_baton_atom(const CAtom* latest_atom_addition,
					   short int direction) const {
   std::vector<clipper::Coord_orth> positions;
   int direction_sign = +1; 

   if (direction == 1) { // building forward, look in negative
			 // direction for previously build atoms.
      direction_sign = +1;
   } else { 
      direction_sign = -1; // building backward, look in positive
			   // direction for previously build atoms.
   }
   int ires_last_atom = ((CAtom *) latest_atom_addition)->GetSeqNum();

   char *chain = ((CAtom *) latest_atom_addition)->GetChainID();
   // does the CA for the (ires_last_atom-2) exist?
   int selHnd = atom_sel.mol->NewSelection();
	 
   atom_sel.mol->SelectAtoms(selHnd, 0, chain,
			     ires_last_atom-2*direction_sign, "*", 
			     ires_last_atom-2*direction_sign, "*",
			     "*", " CA ", "*", "*"); 
      
   int nSelAtoms;
   PPCAtom local_SelAtom; 
   atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);
      
   if (nSelAtoms == 0) {
      std::cout << "residue with sequence number " << ires_last_atom - 2*direction_sign
		<< " not found for ires_last_atom = " << ires_last_atom 
		<< " with direction_sign = " << direction_sign << "\n";
   } else {
      positions.push_back(clipper::Coord_orth(local_SelAtom[0]->x,
					      local_SelAtom[0]->y,
					      local_SelAtom[0]->z));
   }
   atom_sel.mol->DeleteSelection(selHnd);

   // rinse, lather, repeat...

   // does the CA for the (ires_last_atom-1) exist?
   selHnd = atom_sel.mol->NewSelection();
      
   atom_sel.mol->SelectAtoms(selHnd, 0, chain,
			     ires_last_atom-direction_sign, "*", 
			     ires_last_atom-direction_sign, "*",
			     "*", " CA ", "*", "*"); 

   atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);
      
   if (nSelAtoms == 0) {
      std::cout << "residue with sequence number " << ires_last_atom - direction_sign
		<< " not found\n";
   } else {
      positions.push_back(clipper::Coord_orth(local_SelAtom[0]->x,
					      local_SelAtom[0]->y,
					      local_SelAtom[0]->z));
   }
   atom_sel.mol->DeleteSelection(selHnd);

   // And finally this one is guaranteed to exist:
   // 
   positions.push_back(clipper::Coord_orth(latest_atom_addition->x,
					   latest_atom_addition->y,
					   latest_atom_addition->z));

   return positions;
   
} 

#include "CalphaBuild.hh"

std::vector<coot::scored_skel_coord>
molecule_class_info_t::next_ca_by_skel(const std::vector<clipper::Coord_orth> &previous_ca_positions,
				       const clipper::Coord_grid &coord_grid_start,
				       short int use_coord_grid_start_flag,
				       float ca_ca_bond_length,
				       float map_cut_off,
				       int max_skeleton_search_depth) const {

   std::vector<coot::scored_skel_coord> t; 
   coot::CalphaBuild buildca(max_skeleton_search_depth);

    std::cout << "DEBUG:: ------ "
 	     << "in molecule_class_info_t::next_ca_by_skel skeleton_treenodemap_is_filled is "
	      << skeleton_treenodemap_is_filled << " for molecule " << imol_no << std::endl;

    
   if (skeleton_treenodemap_is_filled) { 
      t = buildca.next_ca_by_skel(previous_ca_positions,
				  coord_grid_start,
				  use_coord_grid_start_flag,
				  ca_ca_bond_length,
				  xskel_cowtan, xmap_list[0],
				  map_cut_off,
				  skeleton_treenodemap);
   } else {
      std::cout << "treenodemap is not filled" << std::endl;
   }
   return t;
} 

#include <time.h>

// ----------------------------------------------------------------------
//               Dummy Atoms (not bonded)
// ----------------------------------------------------------------------
void
molecule_class_info_t::add_dummy_atom(coot::Cartesian pos) {

   int nchains = atom_sel.mol->GetNumberOfChains(1);

   if (nchains != 1) {
      std::cout << "failed to add dummy atom" << std::endl;
      return;
   }

   make_backup();

   CChain *chain_p = atom_sel.mol->GetChain(1,0);
   
   std::string mol_chain_id(chain_p->GetChainID());
   int ires_prev = chain_p->GetNumberOfResidues();

   CResidue *res_p = new CResidue;
   CAtom *atom_p = new CAtom;
   chain_p->AddResidue(res_p);
   atom_p->SetAtomName(" DUM");
   atom_p->SetCoordinates(pos.x(), pos.y(), pos.z(), 1.0,
			  graphics_info_t::default_new_atoms_b_factor);

   atom_p->SetElementName(" O");
   res_p->AddAtom(atom_p);
   res_p->seqNum = ires_prev + 1;
   res_p->SetResName("DUM");

   atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
   atom_sel.mol->FinishStructEdit();
   atom_sel = make_asc(atom_sel.mol);
   // std::cout << atom_p << " added to molecule" << std::endl;

   have_unsaved_changes_flag = 1; 
   makebonds(0.0, 0.0);

}



// backups:

// Backup filename: return a stub.
// 
std::string
molecule_class_info_t::save_molecule_filename(const std::string &dir) { 

   std::string time_string = save_time_string;
   graphics_info_t g;
   
   if ((history_index == 0) ||
       history_index != max_history_index) { 

      time_string = dir;

      // unix dependent logic here:  Don't know how to do this on other systems...
      // We want a filename proceeded by a directory name:
      // i.e. we end up with something like
      // "coot-backup/a.pdb_Tues_Aug_19_20:16:00_2003_modification_0.mmdbbin"

      time_string += "/";

      std::string clean_name = name_;
      if (g.unpathed_backup_file_names_flag) {
	 clean_name = name_for_display_manager();
      }
      // convert "/" to "_"
      int slen = clean_name.length();
      for (int i=0; i<slen; i++)
#if defined(__WIN32__) || defined(__CYGWIN__) || defined(WINDOWS_MINGW)
// BL says: we change /, \ and : to _ in windows
         if (clean_name[i] == '/' || clean_name[i] == '\\' 
                             || clean_name[i] == ':')
            clean_name[i] = '_';
#else
	 if (clean_name[i] == '/')
	    clean_name[i] = '_';
#endif // win32 things

      time_string += clean_name;
      time_string += "_";
      
      // add in the time component:

#if defined(__CYGWIN__) || defined(_MSC_VER)

      // but not if we are in windows:

#else      
      time_t t;
      time(&t);
      char *chars_time = ctime(&t);
#ifdef WINDOWS_MINGW
// BL says: why not? We can fix this. I show you how it's done in MINGW:
// dunno if it works in other win32 systems. Havent checked
// we just convert the : to _
      for (int i=0; i<24; i++) {
         if (chars_time[i] == ':') {
             chars_time[i] = '_';
         }
      }
#endif // MINGW
      time_string += chars_time;
#endif // other WIN32

      // strip off the trailing newline:
      slen = time_string.length();
      if (slen > 2) 
	 time_string = time_string.substr(0,slen-1);
      
      // convert spaces to underscores
      // 
      for (unsigned int i=0; i<time_string.length(); i++)
	 if (time_string[i] == ' ')
	    time_string[i] = '_';
	    
#if defined(__WIN32__) || defined(__CYGWIN__) || defined(WINDOWS_MINGW) || defined(_MSC_VER)

      // convert : to underscores in windows
      // 
#ifndef WINDOWS_MINGW
      // BL say: nonsense since we would transform the directory C: here.
      // we have done it before already
      for (int i=0; i<time_string.length(); i++)
	 if (time_string[i] == ':')
	    time_string[i] = '_';
#endif // MINGW
#endif // other win32
      
      time_string += "_modification_";

      save_time_string = time_string; // why do we do this?  Ah, because we want the
                                      // time to calculated at the start:
                                      // and use that as a stub.

      time_string += g.int_to_string(history_index);
      //time_string += ".mmdbbin";
      if (! is_from_shelx_ins_flag) 
	 time_string += ".pdb";
      else 
	 time_string += ".res";

#if defined(__WIN32__) || defined(__CYGWIN__) || defined(WINDOWS_MINGW) || defined(_MSC_VER)
      // we can do now too (I hope for all of them?!?)
      if (! is_from_shelx_ins_flag)
	 time_string += ".gz"; // 'cos we can do compression.  Groovy baby!
#else
      if (! is_from_shelx_ins_flag)
	 time_string += ".gz"; // 'cos we can do compression.  Groovy baby!
#endif
      
   } else {
      // (this is not the first save molecule that we have done)
      
      // add to the stub that we have previously generated.
      //
      time_string += g.int_to_string(history_index);
      if (! is_from_shelx_ins_flag) 
	 time_string += ".pdb";
      else 
	 time_string += ".res";
#if defined(__WIN32__) || defined(__CYGWIN__) || defined(WINDOWS_MINGW) || defined(_MSC_VER)
      // same here
      if (! is_from_shelx_ins_flag) 
	 time_string += ".gz";
#else
      if (! is_from_shelx_ins_flag) 
	 time_string += ".gz";
#endif
   } 
   return time_string;
}

// Return like mkdir: mkdir returns zero on success, or -1 if an  error  occurred
//
// if it already exists as a dir, return 0 of course.
// 
int
molecule_class_info_t::make_maybe_backup_dir(const std::string &backup_dir) const {

   return coot::util::create_directory(backup_dir);
}

// Ignore return value.
//
// If successful, increase history_index and if not in a backup
// increase max_history_index too.
// 
int
molecule_class_info_t::make_backup() { // changes history details

   if (backup_this_molecule) { 
      std::string backup_dir("coot-backup");

      //shall we use the environment variable instead?
      char *env_var = getenv("COOT_BACKUP_DIR");
      if (env_var) { 
	 struct stat buf;
#ifdef WINDOWS_MINGW
         // we better debackslash the directory
	 std::string tmp_dir = env_var;
         tmp_dir = coot::util::intelligent_debackslash(tmp_dir);
	 int err = stat(tmp_dir.c_str(), &buf);
#else
	 int err = stat(env_var, &buf);
#endif // MINGW
	 if (!err) {
	    if (! S_ISDIR(buf.st_mode)) {
	       env_var = NULL;
	    }
	 } else {
	    env_var = NULL;
	 }
      }
      if (env_var)
	 backup_dir = env_var;

      if (atom_sel.mol) {
	 int dirstat = make_maybe_backup_dir(backup_dir);

	 // all is hunkey-dorey.  Directory exists.
	 if (dirstat == 0) { 
	    
	    std::string backup_file_name = save_molecule_filename(backup_dir);
 	    std::cout << "INFO:: backup file " << backup_file_name << std::endl;

#if defined(__WIN32__) || defined(__CYGWIN__) || defined(WINDOWS_MINGW) || defined(_MSC_VER)
            // and again, although not used any more!?
	    byte gz = GZM_ENFORCE;
#else
	    byte gz = GZM_ENFORCE;
#endif
	    // Writing out a modified binary mmdb like this results in the
	    // file being unreadable (crash in mmdb read).
	    // 
	    // int istat = atom_sel.mol->WriteMMDBF((char *)backup_file_name.c_str(), gz);
	    int istat;
	    if (! is_from_shelx_ins_flag) {
	       // istat = atom_sel.mol->WritePDBASCII((char *)backup_file_name.c_str(), gz);
	       istat = write_atom_selection_file(atom_sel, backup_file_name, gz);
	       // WriteMMDBF returns 0 on success, else mmdb:Error_CantOpenFile (15)
	       if (istat) { 
		  std::cout<< "WARNING:: WritePDBASCII failed! Return status "
			   << istat << std::endl;
	       }
	    } else { 
	       std::pair<int, std::string> p = write_shelx_ins_file(backup_file_name);
	       istat = p.first;
	    }
		  
	    save_history_file_name(backup_file_name);
	    if (history_index == max_history_index)
	       max_history_index++;
	    history_index++;
	 } else {
	    std::cout << "BACKUP:: directory "<< backup_dir << " failure" << std::endl;
	 } 
      } else {
	 std::cout << "BACKUP:: Ooops - no atoms to backup for this empty molecule"
		   << std::endl;
      }
   } else {
      // Occasionally useful but mostly tedious...
      // std::cout << "INFO:: backups turned off on this molecule"
      // << std::endl;
   }
   return 0;
}


void
molecule_class_info_t::save_history_file_name(const std::string &file) {

   // First, history_index is zero and the vec is zero,
   // normal service, then another backup: history_index is 1 and vec is 1.
   // 
   if (history_index == int(history_filename_vec.size())) {
      history_filename_vec.push_back(file);
   } else {
      // we have gone back in history.
      // 
      if (history_index < int(history_filename_vec.size())) {
	 history_filename_vec[history_index] = file;
      }
   }
} 

// restore from (previous) backup
void
molecule_class_info_t::restore_from_backup(int history_offset, const std::string &cwd) {

   int hist_vec_index = history_index + history_offset;
   if (int(history_filename_vec.size()) > hist_vec_index) {
      std::cout << "restoring from backup " << history_filename_vec.size()
		<< " " << history_index << std::endl;
      std::string save_name = name_;
      if (hist_vec_index < history_filename_vec.size() &&
	  hist_vec_index >= 0) { 
	 std::string filename = history_filename_vec[hist_vec_index];
	 //      history_index = hist_index;
	 short int reset_rotation_centre = 0;
	 // handle_read_draw_molecule uses graphics_info_t::n_molecules
	 // to determine its molecule number.  We don't want it to
	 // change.
	 int save_imol = imol_no;
	 // similarly, it messes with the save_state_command_strings_, we
	 // don't want that either:
	 std::vector<std::string> save_save_state = save_state_command_strings_;
	 short int is_undo_or_redo = 1;
	 handle_read_draw_molecule(imol_no, filename, cwd, reset_rotation_centre, is_undo_or_redo,
				   bond_width);
	 save_state_command_strings_ = save_save_state;
	 imol_no = save_imol; 
	 name_ = save_name;
      }
   } else {
      std::cout << "not restoring from backup because "
		<< history_filename_vec.size()
		<< " " << history_index << std::endl;
   }
}


// I need to write an essay on how that backup system works.
//
// Insight: if we are at hist_index = max_hist_index (i.e. not in a
// backup) then when an undo is requested, we should make a backup.
// This makes the indexing a bit tricky,
//
// So imagine this situation:
//
//             hist_index max_hist_index                filenames filled
// pepflip         1         1                              [0]
// rotate          2         2                              [0,1]
// undo            3         3 [first step is a backup]     [0,1,2]
//
// filenames get routinely pushed back onto history_filename_vec
// (i.e. not in an undo situation)
//
// So imagine we make one mod, then the history_filename_vec size() is 1.
// on undo:
//    we restore from backup using history_filename_vec index 0.
//    we have added to history_filename_vec [now has size 2] in this proceedure
//    history_index was 1 on starting the undo
//    at end of undo it is 0.
//
// how do we redo that?
//    (obviously) we restore from backup using history_filename_vec index 1.
//    we have not added to history_filename_vec [now has size 2]
//    history_index was 0 on starting the redo.
//    at end of redo, history_index is 1
//
// So having done 2 mods:
//
// on undo:
//    restore from backup using history_filename_vec index 1
//    we have added to history_filename_vec [now has size 3] in this proceedure
//    history_index was 2 on starting undo
//    at end of undo it is 1.
//
// on redo:
//    restore from backup using history_filename_vec index 2
//    we have not added to history_filename_vec [now has size 3]
//    history_index was 1 on starting the redo
//    at end of redo, history_index was 2


//
// [It would be cool to have the Redo button greyed out when there are
// no redos availabile (set its state when either it or undo is
// pressed)]
// 
// initially it should be greyed out (insensitive).

// restore from (next) backup
void
molecule_class_info_t::apply_undo(const std::string &cwd) {

//    std::cout << std::endl << "DEBUG:: in apply undo start hist_index: "
// 	     << history_index
// 	     << " max_history_index: " << max_history_index << std::endl;

   if (history_index > 0) {
      int offset = -1;
      if (history_index == max_history_index) { 
	 make_backup(); // increments history_index
	 offset--;
      }
      restore_from_backup(offset, cwd);
      history_index += offset;

      // So that we don't get asked to save the molecule on exist when
      // we have reverted all our modifications:
      // 
      if (history_index == 0) { 
	 have_unsaved_changes_flag = 0;
      }
   }

   std::cout << "DEBUG:: apply_undo: (end) history_index: " <<
      history_index << " max_history_index: " << max_history_index << std::endl;

}

void
molecule_class_info_t::apply_redo(const std::string &cwd) {

   if (history_index < max_history_index) {
      std::cout << "DEBUG:: molecule applying redo " << history_index << std::endl;

      // When there are 3 backups made and we are viewing molecule 2,
      // we don't want to restore from history_filename_vec[3]:
      // 
      if (int(history_filename_vec.size()) > (history_index + 1)) { 
	 restore_from_backup(+1, cwd); 
	 history_index++; 
	 have_unsaved_changes_flag = 1;
      } else {
	 std::cout << "Not redoing history file vec: " << history_filename_vec.size()
		   << " " << history_index << std::endl;
      } 
   } else {
      std::cout << "Not redoing history: " << max_history_index
		<< " " << history_index << std::endl;
   }
} 



// For model view (go to atom)
//
std::vector<coot::model_view_residue_button_info_t>
molecule_class_info_t::model_view_residue_button_labels() const {

   std::vector<coot::model_view_residue_button_info_t> v;

   if (atom_sel.n_selected_atoms > 0) { 

      int nchains = atom_sel.mol->GetNumberOfChains(1);

      if (nchains < 1) {
	 std::cout << "failed to find chains for atom in "
		   << " model_view_residue_button_info_t" << std::endl;
      } else {

	 graphics_info_t g;
	 CModel *model_p = atom_sel.mol->GetModel(1);

	 CChain *chain;
	 // run over chains of the existing mol
	 int nchains = model_p->GetNumberOfChains();
	 if (nchains <= 0) { 
	    std::cout << "bad nchains in model_view_residue_button_info_t: "
		      << nchains << std::endl;
	 } else { 
	    for (int ichain=0; ichain<nchains; ichain++) {
	       chain = model_p->GetChain(ichain);
	       if (chain == NULL) {  
		  // This should not be necessary. It seem to be a
		  // result of mmdb corruption elsewhere - possibly
		  // DeleteChain in update_molecule_to().
		  std::cout << "NULL chain in model_view_residue_button_info_t: "
			    << std::endl;
	       } else { 
		  int nres = chain->GetNumberOfResidues();
		  for (int ires=0; ires<nres; ires++) { 
		     PCResidue residue_p = chain->GetResidue(ires);
		     std::string button_label = 
			g.int_to_string(residue_p->GetSeqNum());
		     button_label += " ";
		     button_label += residue_p->GetChainID();
		     button_label += " ";
		     button_label += residue_p->name;

		     v.push_back(coot::model_view_residue_button_info_t(button_label,
									residue_p));
		  }
	       }
	    }
	 }
      }
   }
   return v;
}

// return vector of atom list (aka button) info for this residue
// 
std::vector<coot::model_view_atom_button_info_t>
molecule_class_info_t::model_view_atom_button_labels(char *chain_id, int seqno) const {

   graphics_info_t g;
   std::vector<coot::model_view_atom_button_info_t> v;

   // protection against the molecule having been deleted after the
   // gtklist widget was created:
   if (atom_sel.n_selected_atoms > 0) { 

      CChain *chain;
      
      // first we have to find the residue res_p (from which we wil get the atoms)
      //
      int nchains = atom_sel.mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain = atom_sel.mol->GetChain(1, ichain);
	 if (chain == NULL) { 
	    // This should not be necessary. It seem to be a result of
	    // mmdb corruption elsewhere - possibly DeleteChain in
	    // update_molecule_to().
	    std::cout << "ERROR getting chain in model_view_atom_button_info_t\n";
	 } else { 
	    std::string residue_chain_id(chain_id); // passed from residue list
	    std::string mol_chain_id(chain->GetChainID());
	    if (residue_chain_id == mol_chain_id) {
	       int nres = chain->GetNumberOfResidues();
	       for (int ires=0; ires<nres; ires++) { 
		  PCResidue res_p = chain->GetResidue(ires);
		  if (res_p->GetSeqNum() == seqno) {
      
		     PPCAtom residue_atoms;
		     int nResidueAtoms;

		     res_p->GetAtomTable(residue_atoms, nResidueAtoms);
		     for (int i=0; i<nResidueAtoms; i++) {
			if (! residue_atoms[i]->isTer()) { 
			   std::string button_label = residue_atoms[i]->name;
			   std::string altConf = residue_atoms[i]->altLoc;
			   if (altConf != "") { 
			      button_label += ",";
			      button_label += altConf;
			   }
			   button_label += " occ=";
			   button_label += g.float_to_string(residue_atoms[i]->occupancy);
			   button_label += " bf=";
			   button_label += g.float_to_string(residue_atoms[i]->tempFactor);
			   
			   v.push_back(coot::model_view_atom_button_info_t(button_label, residue_atoms[i]));
			}
		     }
		  }
	       }
	    } 
	 } 
      }
   }
   return v;
}


std::vector<coot::model_view_atom_tree_chain_t>
molecule_class_info_t::model_view_residue_tree_labels() const {

   std::vector<coot::model_view_atom_tree_chain_t> v;

   if (atom_sel.n_selected_atoms > 0) {

      CChain *chain_p;
      int im = 1;
      int nchains = atom_sel.mol->GetNumberOfChains(im);
      for (int ichain=0; ichain<nchains; ichain++) {

	 chain_p = atom_sel.mol->GetChain(im, ichain);
	 std::string chain_label("Chain ");
	 chain_label += chain_p->GetChainID();
	 v.push_back(coot::model_view_atom_tree_chain_t(chain_label));
	 
	 if (! chain_p) {
	    std::cout << "ERROR getting chain in model_view_residue_tree_labels\n";
	 } else {
	    int nres = chain_p->GetNumberOfResidues();
	    PCResidue residue_p;
	    for (int ires=0; ires<nres; ires++) {
	       residue_p = chain_p->GetResidue(ires);
	       std::string label = residue_p->GetChainID();
	       label += " ";
	       label += coot::util::int_to_string(residue_p->GetSeqNum());
	       label += residue_p->GetInsCode();
	       label += " ";
	       label += residue_p->name;
	       coot::model_view_atom_tree_item_info_t res(label, residue_p);
	       v.back().add_residue(res);
	    }
	 }
      }
   }
   return v;
}




// Return 0 on failure.
short int 
molecule_class_info_t::move_std_residue(CResidue *moving_residue,
					const CResidue *reference_residue) const {

   std::pair<clipper::RTop_orth, short int> pair = 
      coot::util::get_ori_to_this_res((CResidue *)reference_residue); 

   short int istat = 0; // success

   if (!reference_residue) { 
      std::cout << "This should not happen!" << std::endl;
      std::cout << "null reference residue in move_std_residue" << std::endl;
   } else { 

      if (pair.second == 1) { // successful attempt to get the matrix
	 PPCAtom residue_atoms = NULL;
	 int nResidueAtoms;
	 moving_residue->GetAtomTable(residue_atoms, nResidueAtoms);
	 if (nResidueAtoms == 0) {
	    std::cout << " something broken in atom residue selection in ";
	    std::cout << "mutate, got 0 atoms" << std::endl;
	    istat = 0;
	 } else {
	    istat = 1;
// 	    std::cout << "DEBUG:: move_std_residue: " << nResidueAtoms
// 		      << " atoms in residue " 
// 		      << moving_residue << " " << moving_residue->seqNum << " " 
// 		      << moving_residue->GetChainID() << std::endl;
	    for(int iat=0; iat<nResidueAtoms; iat++) {
	       if (residue_atoms[iat]) { 
// 		  std::cout  << "residue atom " << iat << " coords: " 
// 			     << residue_atoms[iat]->x << " "
// 			     << residue_atoms[iat]->y << " "
// 			     << residue_atoms[iat]->z << std::endl;
		  clipper::Coord_orth co(residue_atoms[iat]->x,
					 residue_atoms[iat]->y,
					 residue_atoms[iat]->z);
		  clipper::Coord_orth rotted = co.transform(pair.first); // an rtop
		  residue_atoms[iat]->x = rotted.x();
		  residue_atoms[iat]->y = rotted.y();
		  residue_atoms[iat]->z = rotted.z();
	       } else { 
		  istat = 0;
		  std::cout << "ERROR:: bad residue atom in move_std_residue: iat: "
			    << iat << std::endl;
	       }
	    }
	 }
      } else { 
	 istat = 0; // failure
	 std::cout << "DISASTER - failed to generate RTop for move_std_residue\n";
	 if (reference_residue) { 
	    // 	 molecule-class-info.cc:4184: passing `const CResidue' as `this' 
	    // argument of `int CResidue::GetSeqNum ()' discards qualifiers
	    CResidue *tmp = (CResidue *) reference_residue;
	    std::cout << "mainchain atoms missing from residue " 
		      << tmp->GetSeqNum() 
		      << tmp->GetChainID() << std::endl;
	 } else { 
	    std::cout << "This should not happen!" << std::endl;
	    std::cout << "null residue in move_std_residue" << std::endl;
	 }
      }
   }
   return istat;
}


void
molecule_class_info_t::make_backup_from_outside() {  // when we have a multi mutate, we
				    // want the wrapper to make a
				    // backup when we start and set
				    // changes when when finish.
				    // Rather crap that this needs to
				    // be done externally, I think.

   make_backup();
}


void
molecule_class_info_t::set_have_unsaved_changes_from_outside() {

   have_unsaved_changes_flag = 1;

}


//
// Get a deep copy:
// return NULL on failure
// 
CResidue *
molecule_class_info_t::get_standard_residue_instance(const std::string &residue_type) {

   graphics_info_t g;
   CResidue *std_residue = NULL;
   
   if (g.standard_residues_asc.read_success) { 
//      std::cout << "DEBUG:: There are " << g.standard_residues_asc.n_selected_atoms
// 	       << " atoms in standard_residues_asc" << std::endl;
     int selHnd = g.standard_residues_asc.mol->NewSelection();
     g.standard_residues_asc.mol->Select ( selHnd,STYPE_RESIDUE, 1, // .. TYPE, iModel
					   "*", // Chain(s) it's "A" in this case.
					   ANY_RES,"*",  // starting res
					   ANY_RES,"*",  // ending res
					   (char *) residue_type.c_str(),  // residue name
					   "*",  // Residue must contain this atom name?
					   "*",  // Residue must contain this Element?
					   "*",  // altLocs
					   SKEY_NEW // selection key
					   );
     // get the standard orientation residue for this residue type
     PPCResidue SelResidue;
     int nSelResidues;

     g.standard_residues_asc.mol->GetSelIndex(selHnd, SelResidue, nSelResidues);
   
     if (nSelResidues != 1) {
       std::cout << "This should never happen - ";
       std::cout << "badness in get_standard_residue_instance, we selected " << nSelResidues
		 << " residues looking for residues of type :" << residue_type << ":\n";
     } else {
       std_residue = coot::deep_copy_this_residue(SelResidue[0], "", 1, 
						  g.standard_residues_asc.UDDAtomIndexHandle);
     }
     g.standard_residues_asc.mol->DeleteSelection(selHnd);
   }
   return std_residue;
}


// return the number of residues in chain with chain_id, return -1 on error
// 
int
molecule_class_info_t::chain_n_residues(const char *chain_id) const {

   int r = -1;
   
   if (atom_sel.n_selected_atoms > 0) {
      int nchains = atom_sel.mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) {
	 const CChain *chain_p = (const CChain *)atom_sel.mol->GetChain(1,ichain);
	 std::string mol_chain_id(((CChain*)chain_p)->GetChainID());
	 if (mol_chain_id == std::string(chain_id)) {
	    r = ((CChain *)chain_p)->GetNumberOfResidues();
	 }
      }
   }
   return r; 
} 


void
molecule_class_info_t::store_refmac_params(const std::string &mtz_filename,
					   const std::string &fobs_col,
					   const std::string &sigfobs_col,
					   const std::string &r_free_col, 
					   int r_free_flag) { 

   have_sensible_refmac_params = 1; // true
   refmac_mtz_filename = mtz_filename; 
   refmac_fobs_col = fobs_col;
   refmac_sigfobs_col = sigfobs_col;
   refmac_r_free_col = r_free_col;
   refmac_r_free_flag_sensible = r_free_flag;

   std::cout << "INFO:: Stored refmac parameters: " 
	     << refmac_fobs_col << " "
	     << refmac_sigfobs_col;
   if (r_free_flag)
      std::cout << " " << refmac_r_free_col << " is sensible." << std::endl;
   else
      std::cout << " the r-free-flag is not sensible" << std::endl;
} 

void
molecule_class_info_t::store_refmac_mtz_filename(const std::string &mtz_filename) { 

   refmac_mtz_filename = mtz_filename; 
}

void
molecule_class_info_t::store_refmac_phase_params(const std::string &phi,
						 const std::string &fom,
						 const std::string &hla,
						 const std::string &hlb,
						 const std::string &hlc,
						 const std::string &hld) {
  
  std::cout << "BL DEBUG:: in store refmac pahse params" <<std::endl;
  have_refmac_phase_params = 1; // true
  refmac_phi_col = phi;
  refmac_fom_col = fom;
  refmac_hla_col = hla; 
  refmac_hlb_col = hlb; 
  refmac_hlc_col = hlc; 
  refmac_hld_col = hld; 
}

void
molecule_class_info_t::store_refmac_file_mtz_filename(const std::string &mtz_filename) { 

   refmac_file_mtz_filename = mtz_filename; 
}


// return 0 on success
// 
int
molecule_class_info_t::write_pdb_file(const std::string &filename) {

   int err = 1; // fail
   if (atom_sel.n_selected_atoms > 0) { 
      std::string ext = coot::util::file_name_extension(filename);
      if (coot::util::extension_is_for_shelx_coords(ext)) {
	 write_shelx_ins_file(filename);
      } else {
	 byte bz = GZM_NONE;
	 err = write_atom_selection_file(atom_sel, filename, bz);
      }
   }
   return err; 
} 


// Add this molecule (typically of waters to this molecule by trying
// to put them into an already-existing solvent chain).  If a solvent
// chain does not already exist, put create a new chain id for the
// water_mol atoms.
//
// All the atoms of water_mol need to be in a
// chain that has a different chain id to all the chains in this
// molecule.  Else fail (return status 0).
// 
int
molecule_class_info_t::insert_waters_into_molecule(const coot::minimol::molecule &water_mol) {

   int istat = 0;  // set to failure initially

   // So run over the the chains of the existing molecule looking for
   // a solvent chain.  If there isn't one we simply use
   // append_to_molecule()
   //
   int nchains = atom_sel.mol->GetNumberOfChains(1);
   CChain *chain_p = NULL;
   CChain *solvent_chain_p = NULL; 
   short int i_have_solvent_chain_flag = 0;
   for (int ichain=0; ichain<nchains; ichain++) { 
      
      chain_p = atom_sel.mol->GetChain(1,ichain);
      if (chain_p->isSolventChain()) {
	 solvent_chain_p = chain_p;
	 std::string mol_chain_id(chain_p->GetChainID());
	 i_have_solvent_chain_flag = 1;
      }
   }


   // For every atom in water_mol, create a new atom and a new residue
   // for it. Add the residue to our model's solvent chain and the
   // atom the the residue (of course).
   //
   if (i_have_solvent_chain_flag == 0) {
      
      // We didn't manage to find a solvent chain.
      // We need to create a new chain.
      chain_p = new CChain;
      atom_sel.mol->GetModel(1)->AddChain(chain_p);
      std::pair<short int, std::string> u = unused_chain_id();
      if (u.first)
	 chain_p->SetChainID(u.second.c_str());
      else 
	 chain_p->SetChainID("Z");
   } else {
      chain_p = solvent_chain_p; // put it back, (kludgey, should use
				 // solvent_chain_p from here, not chain_p).
      // OK, we also need to remove any TER cards that are in that chain_p
      remove_TER_on_last_residue(solvent_chain_p);
   } 

//    std::cout << "Debug:: choose chain " << chain_p->GetChainID()
// 	     << " with have_solvent flag: " << i_have_solvent_chain_flag
// 	     << std::endl;
//    std::cout << "Debug:: isSolvent for each residue of chain: " << std::endl;
//    for (int tmp_r=0; tmp_r<chain_p->GetNumberOfResidues(); tmp_r++) {
//       CResidue *rtmp = chain_p->GetResidue(tmp_r);
//       short int flag = isSolvent(rtmp->name);
// 	 std::cout << rtmp->name << " is solvent? " << flag << std::endl;
//    }
   
   std::pair<short int, int> p = coot::util::max_resno_in_chain(chain_p);
   float bf = graphics_info_t::default_new_atoms_b_factor; // 20.0 by default
   int max_resno;
   if (p.first) { 
      max_resno = p.second;
   } else {
      max_resno = 0;
   }
   if (p.first || (i_have_solvent_chain_flag == 0)) {
      make_backup();
      std::cout << "INFO:: Adding to solvent chain: " << chain_p->GetChainID()
		<< std::endl;
      int prev_max_resno = max_resno;
      CResidue *new_residue_p = NULL;
      CAtom    *new_atom_p = NULL;
      int water_count = 0;
      float occ = 1.0;
      if (is_from_shelx_ins_flag)
	 occ = 11.0;
      for (unsigned int ifrag=0; ifrag<water_mol.fragments.size(); ifrag++) {
	 for (int ires=water_mol[ifrag].min_res_no();
	      ires<=water_mol[ifrag].max_residue_number();
	      ires++) {
	    for (unsigned int iatom=0; iatom<water_mol[ifrag][ires].atoms.size(); iatom++) {
	       new_residue_p = new CResidue;
	       new_residue_p->SetResName("HOH");
	       new_residue_p->seqNum = prev_max_resno + 1 + water_count; 
	       water_count++; 
	       new_atom_p = new CAtom;
	       new_atom_p->SetCoordinates(water_mol[ifrag][ires][iatom].pos.x(),
					  water_mol[ifrag][ires][iatom].pos.y(),
					  water_mol[ifrag][ires][iatom].pos.z(), occ, bf);
	       new_atom_p->SetAtomName(water_mol[ifrag][ires][iatom].name.c_str());
	       new_atom_p->Het = 1; // waters are now HETATMs
	       strncpy(new_atom_p->element, water_mol[ifrag][ires][iatom].element.c_str(), 3);
	       strncpy(new_atom_p->altLoc, water_mol[ifrag][ires][iatom].altLoc.c_str(), 2);

	       // residue number, atom name, occ, coords, b factor

	       // add the atom to the residue and the residue to the chain
	       new_residue_p->AddAtom(new_atom_p);
	       chain_p->AddResidue(new_residue_p);
	    }
	 }
      }
      atom_sel.mol->FinishStructEdit();
      update_molecule_after_additions(); // sets unsaved changes flag
   }
   
   return istat;
}

// Add this molecule (typically of waters to this
// molecule... somehow).  All the atoms of water_mol need to be in a
// chain that has a different chain id to all the chains in this
// molecule.  Else fail (return status 0).
// 
int
molecule_class_info_t::append_to_molecule(const coot::minimol::molecule &water_mol) {

   int istat = 0; // fail status initially.
   int n_atom = 0;  // 0 new atoms added initially.
   
   if (atom_sel.n_selected_atoms > 0) {

      make_backup();

      // run over the chains in water_mol (there is only one for waters)
      //
      for (unsigned int ifrag=0; ifrag<water_mol.fragments.size(); ifrag++) {

// 	 std::cout << "DEBUG:: append_to_molecule: fragment id id for frag " << ifrag
// 		   << " is " << water_mol[ifrag].fragment_id << std::endl;

	 short int imatch = 0;
	 
	 // Run over chains of the existing mol, to see if there
	 // already exists a chain with the same chain id as the
	 // waters we want to add.  Only if imatch is 0 does this
	 // function do anything.
	 // 
	 int nchains = atom_sel.mol->GetNumberOfChains(1);
	 CChain *chain;
	 for (int ichain=0; ichain<nchains; ichain++) { 
	 
	    chain = atom_sel.mol->GetChain(1,ichain);
	    std::string mol_chain_id(chain->GetChainID());
	 
	    if (water_mol.fragments[ifrag].fragment_id == mol_chain_id) {
	       //
	       imatch = 1;
	       istat = 1;
	       std::cout << "INFO:: Can't add waters from additional molecule "
			 << "chain id = " << mol_chain_id << std::endl
			 << "INFO:: That chain id already exists in this molecule"
			 << std::endl;
	       break;
	    }
	 }

	 CModel *model_p = atom_sel.mol->GetModel(1);
	 if (imatch == 0) {
	    // There was not already a chain in this molecule of that name.

	    CChain *new_chain_p;
	    CAtom *new_atom_p;
	    CResidue *new_residue_p;

	    new_chain_p = new CChain;
	    std::cout << "DEBUG INFO:: chain id of new chain :"
		      << water_mol[ifrag].fragment_id << ":" << std::endl;
	    new_chain_p->SetChainID(water_mol[ifrag].fragment_id.c_str());
	    model_p->AddChain(new_chain_p);

	    for (int ires=water_mol[ifrag].min_res_no();
		 ires<=water_mol[ifrag].max_residue_number();
		 ires++) {

	       if (water_mol[ifrag][ires].atoms.size() > 0) {
		  new_residue_p = new CResidue;
		  new_residue_p->seqNum = ires;
		  strcpy(new_residue_p->name, water_mol[ifrag][ires].name.c_str());
		  new_chain_p->AddResidue(new_residue_p);
		  for (unsigned int iatom=0; iatom<water_mol[ifrag][ires].atoms.size(); iatom++) {
		     
		     new_atom_p = new CAtom;
		     new_atom_p->SetAtomName(water_mol[ifrag][ires][iatom].name.c_str());
		     new_atom_p->SetElementName(water_mol[ifrag][ires][iatom].element.c_str());
		     new_atom_p->SetCoordinates(water_mol[ifrag][ires][iatom].pos.x(),
						water_mol[ifrag][ires][iatom].pos.y(),
						water_mol[ifrag][ires][iatom].pos.z(),
						1.0, graphics_info_t::default_new_atoms_b_factor);
		     new_residue_p->AddAtom(new_atom_p);
		     n_atom++; 
		  }
	       }
	    }
	 }
      }

      std::cout << "INFO:: " << n_atom << " atoms added to molecule." << std::endl;
      if (n_atom > 0) { 
	 atom_sel.mol->FinishStructEdit();
	 update_molecule_after_additions(); // sets unsaved changes flag
      }
   }

   return istat;
} 


void
molecule_class_info_t::update_molecule_after_additions() {

   atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);

   atom_sel = make_asc(atom_sel.mol); // does the udd stuff too.

//    std::cout << "old_n_atoms: " << old_n_atoms << " new: " 
// 	     << atom_sel.n_selected_atoms << std::endl;

   have_unsaved_changes_flag = 1;
   make_bonds_type_checked();
} 

std::string
molecule_class_info_t::Refmac_in_name() const {

   return Refmac_name_stub() + "-pre.pdb";

}

std::string
molecule_class_info_t::Refmac_out_name() const {

   return Refmac_name_stub() + ".pdb";

}

std::string
molecule_class_info_t::Refmac_mtz_out_name() const {

   return Refmac_name_stub() + ".mtz";

}

// combine and strip molecule and refmac count to come up with a pdb filename
// for refmac
std::string
molecule_class_info_t::Refmac_name_stub() const {

   // shall we try to take into account ccp4i refmac naming?
   // OK:
   // Here is an example: 
   // demo.pdb         -> demo_refmac1.pdb 
   // demo_refmac1.pdb -> demo_refmac2.pdb 

   std::string refmac_name = "pre-refmac.pdb"; // default

   // First strip off the path of name_:
   std::string stripped_name; 
   // /a/b.mtz -> b.mtz
#ifdef WINDOWS_MINGW
   std::string::size_type islash = coot::util::intelligent_debackslash(name_).find_last_of("/");
#else
   std::string::size_type islash = name_.find_last_of("/");
#endif // MINGW
   if (islash == string::npos) {
      // std::cout << "DEBUG:: slash not found in " << name_ << std::endl;
      stripped_name = name_;
   } else {
      // std::cout << "DEBUG:: slash found at " << islash << std::endl;
      // stripped_name = name_.substr(islash+1, name_.length());
      stripped_name = name_.substr(islash+1);
   } 
   // std::cout << "DEBUG:: stripped_name: " << stripped_name << std::endl;      

   
   std::string::size_type irefmac = stripped_name.rfind("-refmac");
   std::string::size_type irefmac_ccp4i = stripped_name.rfind("_refmac");
   
   if (irefmac == string::npos) { // not found

      // so was it a ccp4i refmac pdb file?

      if ( ! (irefmac_ccp4i == string::npos) ) { 
	 // it *was* a ccp4i pdb file:
	 //
	 refmac_name = stripped_name.substr(0,irefmac_ccp4i) + "_refmac";
	 refmac_name += graphics_info_t::int_to_string(refmac_count);
      }
      // std::cout << "DEBUG:: irefmac not found in " << stripped_name  << std::endl;
      // lets strip off ".pdb", ".pdb.gz"
      std::string::size_type ipdb = stripped_name.rfind(".pdb");

      if (ipdb == string::npos) { // not a pdb

	 // std::cout << "DEBUG:: ipdb not found" << std::endl;
	 // just tack "refmac-2.pdb" on to the name then
	 refmac_name = stripped_name + "_refmac"; 
	 refmac_name += graphics_info_t::int_to_string(refmac_count); 

      } else {
	 // is a pdb:

	 // std::cout << "DEBUG:: ipdb *was* found" << std::endl;
	 refmac_name = stripped_name.substr(0,ipdb) + "_refmac"; 
	 refmac_name += graphics_info_t::int_to_string(refmac_count); 
      }
   } else {

      // refmac *was* found as part of the name
      // std::cout << "DEBUG:: irefmac *was* found in " << stripped_name << std::endl;
      refmac_name = stripped_name.substr(0,irefmac) + "_refmac";
      refmac_name += graphics_info_t::int_to_string(refmac_count);
   }

   // std::cout << "DEBUG:: returning refmac_name: " << refmac_name << std::endl;
   return refmac_name;

}


std::string
molecule_class_info_t::name_sans_extension(short int include_path_flag) const {

   std::string outstring = name_;
   
   std::string::size_type ipdb = name_.rfind(".pdb");
   if (ipdb != std::string::npos)
      outstring = name_.substr(0, ipdb);

   std::string::size_type islash = outstring.rfind("/");
   if (islash != std::string::npos)
      outstring = outstring.substr(islash+1);
      
   return outstring;
}

void
molecule_class_info_t::update_molecule_to(std::vector<coot::scored_skel_coord> &pos_position) {

   std::cout << "DEBUG:: molecule_class_info_t update_molecule_to() with " << pos_position.size()
	     << " skeleton positions" << std::endl;

   if (has_model()) {
      CModel *model_p = atom_sel.mol->GetModel(1); 
      
      if (! model_p) { 
	 std::cout << "ERROR:: Disaster in finding model_p in update_molecule_to" 
		   << std::endl;
      } else {
	 CChain *chain_p;
	 int n_chains = atom_sel.mol->GetNumberOfChains(1);
	 for (int i_chain=0; i_chain<n_chains; i_chain++) {
	    model_p->DeleteChain(i_chain);
	 }

	 // Now add new chain to model:
	 chain_p = new CChain;
	 if (chain_p) {
	    model_p->AddChain(chain_p);
	    add_multiple_dummies(chain_p, pos_position);
	 } else { 
	    std::cout << "ERROR:: creating chain in mol::update_molecule_to" << std::endl;
	 }
      }
   } else {
      std::cout << "WARNING:: strange! This is not a valid model molecule. " << std::endl;
   }
}

// function callable from graphics_info_t::create_molecule_and_display()
// 
void
molecule_class_info_t::add_multiple_dummies(const std::vector<coot::scored_skel_coord> &pos_position) {

   if (has_model()) {
      CModel *model_p = atom_sel.mol->GetModel(1);
      int n_chains = atom_sel.mol->GetNumberOfChains(1);
      if (n_chains > 0) {
	 CChain *chain_p = model_p->GetChain(0);
	 add_multiple_dummies(chain_p, pos_position);
      }
   }
}


// we presume that the chain exists.  This exists so that we dont do a
// backup every time we add a dummy atom (as is done using
// add_dummy_atom().
// 
void
molecule_class_info_t::add_multiple_dummies(CChain *chain_p,
					    const std::vector<coot::scored_skel_coord> &pos_position) {


   if (pos_position.size() > 0) {
      make_backup(); // maybe
   }
   
   for (unsigned int i=0; i<pos_position.size(); i++) {
      CResidue *res_p = new CResidue;
      CAtom *atom_p = new CAtom;
      chain_p->AddResidue(res_p);
      atom_p->SetAtomName(" DUM");
      atom_p->SetCoordinates(pos_position[i].position.x(),
			     pos_position[i].position.y(),
			     pos_position[i].position.z(), 1.0,
			     graphics_info_t::default_new_atoms_b_factor);

      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
      res_p->seqNum = i + 1;
      res_p->SetResName("DUM");

      // std::cout << atom_p << " added to molecule" << std::endl;
   }

   // std::cout << "DEBUG:: add_multiple_dummies finishing.. "
   // << pos_position.size() << std::endl;
   // if (pos_position.size() > 0) {
   
   // Actually, we want to run this code when there are no new guide
   // points too.  This sets atom_sel.SelectionHandle properly, which
   // is needed in close_yourself, where a DeleteSelection() is done
   // to give back the memory.
      atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      have_unsaved_changes_flag = 1; 
      makebonds(0.0, 0.0);
      // }
} 

void
molecule_class_info_t::add_multiple_dummies(const std::vector<coot::Cartesian> &pos_position) {

   if (atom_sel.mol) {
      CModel *model_p = atom_sel.mol->GetModel(1);
      int n_chains = atom_sel.mol->GetNumberOfChains(1);
      if (n_chains > 0) {
	 CChain *chain_p = model_p->GetChain(0);
	 if (pos_position.size() > 0) {
	    make_backup(); // maybe
	    
	    for (unsigned int i=0; i< pos_position.size(); i++) { 
	       CResidue *res_p = new CResidue;
	       CAtom *atom_p = new CAtom;
	       chain_p->AddResidue(res_p);
	       atom_p->SetAtomName(" DUM");
	       atom_p->SetCoordinates(pos_position[i].x(),
				      pos_position[i].y(),
				      pos_position[i].z(), 1.0,
				      graphics_info_t::default_new_atoms_b_factor);
	    
	       atom_p->SetElementName(" O");
	       res_p->AddAtom(atom_p);
	       res_p->seqNum = i + 1;
	       res_p->SetResName("DUM");
	    }
	    atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
	    atom_sel.mol->FinishStructEdit();
	    atom_sel = make_asc(atom_sel.mol);
	    have_unsaved_changes_flag = 1; 
	    makebonds(0.0, 0.0);
	 }
      }
   }
}


// return an empty vector on failure, a vector of size 6 on success:
// 
std::pair<std::vector<float>, std::string> 
molecule_class_info_t::get_cell_and_symm() const { 

   std::pair<std::vector<float>, std::string> cell_spgr;
   
   mat44 my_matt;
   if (atom_sel.mol) { 
      int err = atom_sel.mol->GetTMatrix(my_matt, 0, 0, 0, 0);
      if (err != 0) {
	 std::cout << "!! Warning:: No symmetry available for this template molecule"
		   << std::endl;
      } else {
	 realtype a[6];
	 realtype vol;
	 int orthcode;
	 atom_sel.mol->GetCell(a[0], a[1], a[2], a[3], a[4], a[5], vol, orthcode);
	 for (int i=0; i<6; i++) cell_spgr.first.push_back(a[i]);
	 cell_spgr.second = std::string(atom_sel.mol->GetSpaceGroup());
      }
   }
   return cell_spgr;
} 

void 
molecule_class_info_t::set_mmdb_cell_and_symm(std::pair<std::vector<float>, std::string> cell_spgr) {

   if (cell_spgr.first.size() == 6) { 
      std::vector<float> a = cell_spgr.first; // short name
      atom_sel.mol->SetCell(a[0], a[1], a[2], a[3], a[4], a[5]);
      atom_sel.mol->SetSpaceGroup((char *)cell_spgr.second.c_str());
      std::cout << "successfully set cell and symmetry" << std::endl;
   } else { 
      std::cout << "WARNING:: failure to set cell on this molecule" << std::endl;
   } 
}

bool 
molecule_class_info_t::set_mmdb_symm(const std::string &spg) {

   atom_sel.mol->SetSpaceGroup(spg.c_str());
   std::string new_sg;
   const char *new_sg_chars = atom_sel.mol->GetSpaceGroup();
   if (new_sg_chars)
      new_sg = new_sg_chars;
   return (new_sg == spg);
} 

// Return atom_index of -1 when no nearest atom.
// 
std::pair<float, int>
molecule_class_info_t::nearest_atom(const coot::Cartesian &pos) const { 

   float min_dist = 999999999; 
   float d;
   int atom_index = -1;

   for(int i=0; i<atom_sel.n_selected_atoms; i++) { 
      coot::Cartesian a(atom_sel.atom_selection[i]->x, atom_sel.atom_selection[i]->y, atom_sel.atom_selection[i]->z);
      d =  fabs((pos-a).length());
      if (d < min_dist) { 
	 min_dist = d;
	 atom_index = i;
      }
   }

   std::pair<float, int> r;
   r.first = min_dist;
   r.second = atom_index;
   return r;
}

// return an empty
// vector
// if closed or is a coords
// mol, 3 elements and 6
// elements for a
// difference map.
std::vector<float>
molecule_class_info_t::map_colours() const {

   std::vector<float> v;
   if (has_map()) {
      if (is_difference_map_p()) {
	 v.resize(6,0.3);
	 v[0] = map_colour[0][0];
	 v[1] = map_colour[0][1];
	 v[2] = map_colour[0][2];
	 v[3] = map_colour[1][0];
	 v[4] = map_colour[1][1];
	 v[5] = map_colour[1][2];
      } else {
	 v.resize(3,0.3);
	 v[0] = map_colour[0][0];
	 v[1] = map_colour[0][1];
	 v[2] = map_colour[0][2];
      }
   }
   return v;
}

// perhaps there is a better place for this?
// 
std::vector<std::string> 
molecule_class_info_t::set_map_colour_strings() const { 

   // return something like
   // (list "set_last_map_colour" "0.2" "0.3" "0.4")


   std::vector<std::string> r;

   r.push_back("set-last-map-colour");
   r.push_back(graphics_info_t::float_to_string(map_colour[0][0]));
   r.push_back(graphics_info_t::float_to_string(map_colour[0][1]));
   r.push_back(graphics_info_t::float_to_string(map_colour[0][2]));

   return r;   
} 


// and symm labels.
void
molecule_class_info_t::remove_atom_labels() { 

   labelled_atom_index_list.clear();
   labelled_symm_atom_index_list.clear();
}

coot::minimol::molecule
molecule_class_info_t::eigen_flip_residue(const std::string &chain_id, int resno) {


   coot::minimol::molecule m;
   
   CResidue *res = get_residue(chain_id, resno, "");
   if (!res) {
      std::cout << "DEBUG:: residue not found " << chain_id << " " << resno
		<< " in molecule number " << MoleculeNumber()
		<< std::endl;
   } else { 
      // make_backup();

      coot::ligand lig;
      coot::minimol::residue r(res);
      coot::minimol::fragment f(res->GetChainID());
      f.residues.push_back(res);
      coot::minimol::molecule ligand;
      ligand.fragments.push_back(f);
      
      ligand_flip_number++;
      if (ligand_flip_number == 4)
	 ligand_flip_number = 0;

      lig.install_ligand(ligand);
      m = lig.flip_ligand(ligand_flip_number);
      
      make_bonds_type_checked();
      have_unsaved_changes_flag = 1;

      replace_coords(make_asc(m.pcmmdbmanager()), 0, 1);
   }
   return m;
} 



// So that we can move around all the atoms of a ligand (typically)
// 
void 
molecule_class_info_t::translate_by(float x, float y, float z) { 

   if (has_model()) {
      make_backup();
      for (int i=0; i<atom_sel.n_selected_atoms; i++) {
	 atom_sel.atom_selection[i]->x += x;
	 atom_sel.atom_selection[i]->y += y;
	 atom_sel.atom_selection[i]->z += z;
      }
      make_bonds_type_checked();
      have_unsaved_changes_flag = 1;
   }
} 

// Sets coot_save_index maybe (using set_coot_save_index()).
// 
std::string 
molecule_class_info_t::stripped_save_name_suggestion() { 

   std::string s;
   
   std::string stripped_name1;
#ifdef WINDOWS_MINGW
   std::string::size_type islash = coot::util::intelligent_debackslash(name_).find_last_of("/");
#else
   std::string::size_type islash = name_.find_last_of("/");
#endif // MINGW
   if (islash == std::string::npos) { 
      stripped_name1 = name_;
   } else {
      stripped_name1 = name_.substr(islash+1, name_.length());
   }
   // so we have got rid of the pathname.
   // now lets get rid of the extension
   // 
   std::string::size_type ibrk   = stripped_name1.rfind(".brk");
   std::string::size_type ibrkgz = stripped_name1.rfind(".brk.gz");
   std::string::size_type ipdb   = stripped_name1.rfind(".pdb");
   std::string::size_type ires   = stripped_name1.rfind(".res");
   std::string::size_type ipdbgz = stripped_name1.rfind(".pdb.gz");
   std::string::size_type icoot  = stripped_name1.rfind("-coot-");

   std::string stripped_name2;
   if (icoot == std::string::npos) { 
      if (ibrk == std::string::npos) { 
	 if (ibrkgz == std::string::npos) { 
	    if (ipdb == std::string::npos) { 
	       if (ires == std::string::npos) { 
		  if (ipdbgz == std::string::npos) { 
		     stripped_name2 = stripped_name1;
		  } else {
		     stripped_name2 = stripped_name1.substr(0, ipdbgz);
		  }
	       } else {
		  stripped_name2 = stripped_name1.substr(0, ires);
	       }
	    } else {
	       stripped_name2 = stripped_name1.substr(0, ipdb);
	    }
	 } else {
	    stripped_name2 = stripped_name1.substr(0, ibrkgz);
	 }
      } else { 
	 stripped_name2 = stripped_name1.substr(0, ibrk);
      }
   } else {
      set_coot_save_index(stripped_name1.substr(icoot));
      stripped_name2 = stripped_name1.substr(0,icoot);
   }

   stripped_name2 += "-coot-";
   stripped_name2 += graphics_info_t::int_to_string(coot_save_index);
   // As per George Sheldrick's suggestion, if this was from shelx,
   // suggest a .ins extension, not .pdb
   if (!is_from_shelx_ins_flag) {
      stripped_name2 += ".pdb";
   } else { 
      stripped_name2 += ".ins";
   }

//    std::cout << "DEBUG:: stripped_save_name_suggestion: " 
// 	     << stripped_name2 << std::endl;

   return stripped_name2;
} 

int 
molecule_class_info_t::set_coot_save_index(const std::string &filename) { 

   // std::cout << "extracting from :" << filename << std::endl;

   // filename is something like: "-coot-12.pdb".
   // 
   // We want to find 12 and set coot_save_index to 12 + 1, which gets
   // used to suggest the next saved filename.
   // 

   std::string twelve_pdb = filename.substr(6);
   // std::cout << "twelve_pdb:"<< twelve_pdb << std::endl;
   
   std::string::size_type ipdb   = twelve_pdb.rfind(".pdb");
   if (ipdb != std::string::npos) { 
      // .pdb was found
      std::string twelve = twelve_pdb.substr(0,ipdb);
      int i = atoi(twelve.c_str());
      // std::cout << "found i: " << i << std::endl;
      if (i >= 0 && i<100000)
	 coot_save_index = i+1;
   } 
   return coot_save_index; 
} 


void
molecule_class_info_t::transform_by(mat44 mat) { 

#ifdef HAVE_GSL
   if (has_model()) { 
      clipper::Coord_orth co;
      clipper::Coord_orth trans_pos; 
      make_backup();
      clipper::Mat33<double> clipper_mat(mat[0][0], mat[0][1], mat[0][2],
					 mat[1][0], mat[1][1], mat[1][2],
					 mat[2][0], mat[2][1], mat[2][2]);
      clipper::Coord_orth cco(mat[0][3], mat[1][3], mat[2][3]);
      clipper::RTop_orth rtop(clipper_mat, cco);
      std::cout << "INFO:: coordinates transformed by orthonal matrix: \n"
		<< rtop.format() << std::endl;
      clipper::Rotation rtn( clipper_mat );
      clipper::Polar_ccp4 polar = rtn.polar_ccp4();
      clipper::Euler_ccp4 euler = rtn.euler_ccp4();
      std::cout << "  Rotation - polar (omega,phi,kappa)  " << clipper::Util::rad2d(polar.omega()) << " " << clipper::Util::rad2d(polar.phi()) << " " << clipper::Util::rad2d(polar.kappa()) << std::endl;
      std::cout << "  Rotation - euler (alpha,beta,gamma) " << clipper::Util::rad2d(euler.alpha()) << " " << clipper::Util::rad2d(euler.beta()) << " " << clipper::Util::rad2d(euler.gamma()) << std::endl;
      std::cout << "  Translation - Angstroms             " << cco.x() << " " << cco.y() << " " << cco.z() << " " << std::endl;
      for (int i=0; i<atom_sel.n_selected_atoms; i++) { 
	 // atom_sel.atom_selection[i]->Transform(mat); // doesn't compile!
	 // Argh.  sigh.  Use clipper. c.f. graphics_info_t::fill_hybrid_atoms()
	 co = clipper::Coord_orth(atom_sel.atom_selection[i]->x, 
				  atom_sel.atom_selection[i]->y, 
				  atom_sel.atom_selection[i]->z);
	 trans_pos = co.transform(rtop);
	 atom_sel.atom_selection[i]->x = trans_pos.x();
	 atom_sel.atom_selection[i]->y = trans_pos.y();
	 atom_sel.atom_selection[i]->z = trans_pos.z();
      } 
      atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked();
   }
#endif // HAVE_GSL
}


void
molecule_class_info_t::transform_by(const clipper::RTop_orth &rtop) {

   make_backup();
   std::cout << "INFO:: coordinates transformed by orthonal matrix: \n"
	     << rtop.format() << std::endl;
   if (have_unit_cell) {
      clipper::Cell cell(clipper::Cell_descr(atom_sel.mol->get_cell().a,
					     atom_sel.mol->get_cell().b,
					     atom_sel.mol->get_cell().c,
					     clipper::Util::d2rad(atom_sel.mol->get_cell().alpha),
					     clipper::Util::d2rad(atom_sel.mol->get_cell().beta),
					     clipper::Util::d2rad(atom_sel.mol->get_cell().gamma))); 
      std::cout << "INFO:: fractional coordinates matrix:" << std::endl;
      std::cout << rtop.rtop_frac(cell).format() << std::endl;
   } else {
      std::cout << "No unit cell for this molecule, hence no fractional matrix." << std::endl;
   }
   clipper::Coord_orth co;
   clipper::Coord_orth trans_pos; 
   if (has_model()) { 
      for (int i=0; i<atom_sel.n_selected_atoms; i++) {
	 co = clipper::Coord_orth(atom_sel.atom_selection[i]->x, 
				  atom_sel.atom_selection[i]->y, 
				  atom_sel.atom_selection[i]->z);
	 trans_pos = co.transform(rtop);
	 atom_sel.atom_selection[i]->x = trans_pos.x();
	 atom_sel.atom_selection[i]->y = trans_pos.y();
	 atom_sel.atom_selection[i]->z = trans_pos.z();
      }
      atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked();
   }
}

void
molecule_class_info_t::transform_by(const clipper::RTop_orth &rtop, CResidue *residue_moving) {

   make_backup();
   std::cout << "INFO:: coordinates transformed_by: \n"
	     << rtop.format() << std::endl;
   if (has_model()) {
      transform_by_internal(rtop, residue_moving);
      atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked();
   }
}



void
molecule_class_info_t::transform_by_internal(const clipper::RTop_orth &rtop, CResidue *residue_moving) {

   if (has_model()) {
      clipper::Coord_orth co;
      clipper::Coord_orth trans_pos; 
      PPCAtom residue_atoms;
      int n_residue_atoms;
      residue_moving->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iatom=0; iatom<n_residue_atoms; iatom++) { 
	 clipper::Coord_orth p(residue_atoms[iatom]->x,
			       residue_atoms[iatom]->y,
			       residue_atoms[iatom]->z);
	 clipper::Coord_orth p2 = p.transform(rtop);
	 residue_atoms[iatom]->x = p2.x();
	 residue_atoms[iatom]->y = p2.y();
	 residue_atoms[iatom]->z = p2.z();
      }
   }
}


void
molecule_class_info_t::transform_zone_by(const std::string &chain_id, int resno_start, int resno_end,
					 const std::string &ins_code,
					 const clipper::RTop_orth &rtop,
					 bool make_backup_flag) {

   if (make_backup_flag)
      make_backup();
   
   bool transformed_something = 0;
   if (resno_end < resno_start)
      std::swap(resno_end, resno_start);

   int imod = 1;
   CModel *model_p = atom_sel.mol->GetModel(imod);
   CChain *chain_p;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      if (chain_id == chain_p->GetChainID()) { 
	 int nres = chain_p->GetNumberOfResidues();
	 CResidue *residue_p;
	 CAtom *at;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    int this_resno = residue_p->GetSeqNum();
	    std::string this_ins_code = residue_p->GetInsCode();
	    if ((this_resno >= resno_start) &&
		(this_resno <= resno_end)) { 
	       if (this_ins_code == ins_code) {
		  transform_by_internal(rtop, residue_p);
		  transformed_something = 1;
	       }
	    }
	 }
      }
   }

   if (transformed_something) {
      atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked();
   } 
} 



void
molecule_class_info_t::set_contour_level(float f) {
   if (has_map()) { 
      contour_level[0] = f;
      update_map();
   }
}

void
molecule_class_info_t::set_contour_level_by_sigma(float f) {
   if (has_map()) { 
      contour_level[0] = f * map_sigma_;
      update_map();
   }
}

std::vector <std::string>
molecule_class_info_t::get_map_contour_strings() const {

   std::vector <std::string> s; 
   s.push_back("set-last-map-contour-level");
   char cs[100];
   snprintf(cs, 99, "%e", contour_level[0]);
   s.push_back(cs);

   return s;
} 

std::vector <std::string>
molecule_class_info_t::get_map_contour_sigma_step_strings() const {

   std::vector <std::string> s; 
   s.push_back("set-last-map-sigma-step");
   s.push_back(graphics_info_t::float_to_string(contour_sigma_step));
   
   return s;
}

short int
molecule_class_info_t::contoured_by_sigma_p() const { 
   return contour_by_sigma_flag;
} 



CChain *
molecule_class_info_t::water_chain() const { 

   CChain *water_chain = 0;

   if (has_model()) { 

      CModel *model_p = atom_sel.mol->GetModel(1);
      CResidue *residue_p;
      CChain *chain_p;
      
      if (model_p) {

	 if (is_from_shelx_ins_flag) {
	    water_chain = water_chain_from_shelx_ins();
	 } else { 
	    int nchains = model_p->GetNumberOfChains();
	    for (int ich=0; ich<nchains; ich++) {
	       chain_p = model_p->GetChain(ich);
	       int nres = chain_p->GetNumberOfResidues();
	       short int all_water_flag = 1; 
	       for (int ires=0; ires<nres; ires++) { 
		  residue_p = chain_p->GetResidue(ires);
		  std::string resname(residue_p->name);
		  if (! ( (resname == "WAT") || (resname == "HOH"))) {
		     all_water_flag = 0;
		     break;
		  }
	       }
	       if (all_water_flag) { 
		  water_chain = chain_p;
		  break;
	       }
	    }
	 }
      } 
   }
   return water_chain;
}


// there is only one chain from a shelxl ins file.
CChain *
molecule_class_info_t::water_chain_from_shelx_ins() const {

   CChain *water_chain = 0;
   CModel *model_p = atom_sel.mol->GetModel(1);

   if (has_model()) { 
      int nchains = model_p->GetNumberOfChains();
      for (int ich=0; ich<nchains; ich++) {
	 water_chain = model_p->GetChain(ich);
      }
   }
   return water_chain;
}


// return state, max_resno + 1, or 0, 1 of no residues in chain.
// 
std::pair<short int, int> 
molecule_class_info_t::next_residue_in_chain(CChain *w) const { 

   std::pair<short int, int> p(0,1);
   int max_res_no = -9999;

   if (w) { 
      int nres = w->GetNumberOfResidues();
      CResidue *residue_p;
      if (nres > 0) { 
	 for (int ires=nres-1; ires>=0; ires--) { 
	    residue_p = w->GetResidue(ires);
	    if (residue_p->seqNum > max_res_no) {
	       max_res_no = residue_p->seqNum;
	       p = std::pair<short int, int>(1, max_res_no+1);
	    }
	 }
      }
   }
   return p;
}



// add a factor to scale the colours in b factor representation:.
// It goes into the atom_sel.mol
void
molecule_class_info_t::set_b_factor_bonds_scale_factor(float f) {

   std::cout << "Here Adding b-factor scale " << f << std::endl;
   if (atom_sel.mol) {
      // bleugh, casting.
      int udd_handle =
	 atom_sel.mol->RegisterUDReal(UDR_HIERARCHY,
				      (char *) coot::b_factor_bonds_scale_handle_name.c_str());
      if (udd_handle > 0) {
// 	 std::cout << "Adding b-factor scale " << f << " with handle "
// 		   << udd_handle << std::endl;
	 atom_sel.mol->PutUDData(udd_handle, f);

	 // test getting the uddata:
	 int udd_b_factor_handle =
	    atom_sel.mol->GetUDDHandle(UDR_HIERARCHY, (char *) coot::b_factor_bonds_scale_handle_name.c_str());
// 	 std::cout << "debug:: test Got b factor udd handle: "
// 		   << udd_b_factor_handle << std::endl;
	 if (udd_b_factor_handle > 0) {
	    realtype scale;
	    if (atom_sel.mol->GetUDData(udd_b_factor_handle, scale) == UDDATA_Ok) {
// 	       std::cout << " test got b factor scale: " << scale << std::endl;
	    } else {
 	       std::cout << "ERROR:: bad get b factor scale " << std::endl;
	    }
	 }
      }
   }
   make_bonds_type_checked();
}

std::pair<bool, std::string>
molecule_class_info_t::chain_id_for_shelxl_residue_number(int shelxl_resno) const {

   int imod = 1;
   bool found_it = 0;
   std::string chain_id_unshelxed = "not-found";
      
   CModel *model_p = atom_sel.mol->GetModel(imod);
   CChain *chain_p;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      PCResidue residue_p;
      for (int ires=0; ires<nres; ires++) { 
	 residue_p = chain_p->GetResidue(ires);
	 int resno = residue_p->GetSeqNum();
	 if (resno == shelxl_resno) {
	    chain_id_unshelxed = chain_p->GetChainID();
	    found_it = 1;
	 }

	 if (found_it)
	    break;
      }
      if (found_it)
	 break;
   }

   return std::pair<bool, std::string> (found_it, chain_id_unshelxed);
} 


void
molecule_class_info_t::debug() const {

   int imod = 1;
      
   CModel *model_p = atom_sel.mol->GetModel(imod);
   CChain *chain_p;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      PCResidue residue_p;
      for (int ires=0; ires<nres; ires++) { 
	 residue_p = chain_p->GetResidue(ires);
	 if (residue_p) {
	    std::cout << "   " << chain_p->GetChainID() << " " << residue_p->GetSeqNum()
		      << " " << residue_p->index << std::endl;
	 }
      }
   }
}

void
molecule_class_info_t::clear_all_fixed_atoms() {
   fixed_atom_specs.clear();
   fixed_atom_positions.clear();
}

void
molecule_class_info_t::mark_atom_as_fixed(const coot::atom_spec_t &atom_spec, bool state) {

   if (has_model()) {
      int imod = 1;
      bool found = 0;
      CModel *model_p = atom_sel.mol->GetModel(imod);
      CChain *chain_p;
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 std::string chain_id_model = chain_p->GetChainID();
	 if (atom_spec.chain == chain_id_model) { 
	    int nres = chain_p->GetNumberOfResidues();
	    PCResidue residue_p;
	    CAtom *at;
	    for (int ires=0; ires<nres; ires++) { 
	       residue_p = chain_p->GetResidue(ires);
	       int resno_model = residue_p->GetSeqNum();
	       if (resno_model == atom_spec.resno) { 
		  int n_atoms = residue_p->GetNumberOfAtoms();
		  
		  for (int iat=0; iat<n_atoms; iat++) {
		     at = residue_p->GetAtom(iat);
		     if (atom_spec.matches_spec(at)) {
			if (state) {
			   // try to get the atom index of this atom
			   // (at) and make it part of the spec.  So
			   // that we can use it again in
			   // update_fixed_atom_positions().  If the
			   // atom index is correct, we don't need to
			   // search the molecule for the spec.
			   coot::atom_spec_t atom_spec_local = atom_spec;
			   int idx = get_atom_index(at);
			   atom_spec_local.int_user_data = idx;
			   fixed_atom_specs.push_back(atom_spec_local);
			   std::cout << "INFO:: " << atom_spec << " marked as fixed"
				     << std::endl;
			   found = 1; 
			} else {
			   //  try to remove at from marked list
			   if (fixed_atom_specs.size() > 0) {
			      std::vector<coot::atom_spec_t>::iterator it;
			      for (it=fixed_atom_specs.begin();
				   it != fixed_atom_specs.end();
				   it++) {
				 if (atom_spec == *it) {
				    std::cout << "INFO:: removed " << atom_spec
					      << " from fixed atom." << std::endl;
				    fixed_atom_specs.erase(it);
				    found = 1;
				    break;
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
      if (found) {
	 update_fixed_atom_positions();
      }
   }
}


int
molecule_class_info_t::move_waters_to_around_protein() {

   make_backup();
   int r = coot::util::move_waters_around_protein(atom_sel.mol);
   have_unsaved_changes_flag = 1; 
   make_bonds_type_checked();
   return r;
} 



// Return the maximum minimum distance of waters to protein atoms.
// return something negative when we can't do above (no protein
// atoms or no water atoms).
float
molecule_class_info_t::max_water_distance() {

   // Do not account for alt confs of the water positions. That
   // sophistication can ome later if we use this function for
   // something other than testing a greg-test.

   float f = -1.0;

   std::vector<clipper::Coord_orth> protein_positions;
   std::vector<clipper::Coord_orth> water_positions;

   for (int iat=0; iat<atom_sel.n_selected_atoms; iat++) {
      clipper::Coord_orth pt(atom_sel.atom_selection[iat]->x,
			     atom_sel.atom_selection[iat]->y,
			     atom_sel.atom_selection[iat]->z);
      std::string res_name(atom_sel.atom_selection[iat]->GetResName());
      if (res_name == "HOH" || res_name == "WAT") {
	 water_positions.push_back(pt);
      } else {
	 protein_positions.push_back(pt);
      }
   }

   // now protein_positions and water_positions are filled.
   if (protein_positions.size() > 0) { 
      if (water_positions.size() > 0) {
	 double max_water_dist_2 = -1.0;
	 for (unsigned int iw=0; iw<water_positions.size(); iw++) { 
	    double best_dist_2 = 999999999.9;
	    for (unsigned int ip=0; ip<protein_positions.size(); ip++) {
	       double d2 = (water_positions[iw]-protein_positions[ip]).lengthsq();
	       if (d2 < best_dist_2) {
		  best_dist_2 = d2;
	       }
	    }
	    if (best_dist_2 > max_water_dist_2) {
	       max_water_dist_2 = best_dist_2;
	    } 
	 }
	 if (max_water_dist_2 > 0.0)
	    f = sqrt(max_water_dist_2);
      }
   }
   return f;
} 


// ---- utility function --- (so that we know to delete hydrogens
// from HETATM molecule before merging with this one)
//
bool
molecule_class_info_t::molecule_has_hydrogens() const {

   bool r = 0;
   for (int i=0; i<atom_sel.n_selected_atoms; i++) {
      std::string ele(atom_sel.atom_selection[i]->element);
      if (ele == " H") {
	 r = 1;
	 break;
      } 
      if (ele == " D") {
	 r = 1;
	 break;
      }
   }
   return r;
}

// -------- simply print it (at the moment) --------------
void
molecule_class_info_t::print_secondary_structure_info() {

   int n_models = atom_sel.mol->GetNumberOfModels();
   for (int imod=1; imod<=n_models; imod++) { 
      CModel *model_p = atom_sel.mol->GetModel(imod);
      coot::util::print_secondary_structure_info(model_p);
   }
}

