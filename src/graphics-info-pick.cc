/* src/graphics-info-pick.cc
 * 
 * Copyright 2004, 2005 by The University of York
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
#endif

#include "compat/coot-sysdep.h"

#include <string>
#include <vector> // for mmdb-crystal

#include <math.h>

#include <string.h> // strncmp

#include <mmdb2/mmdb_manager.h>
#include "coords/mmdb-extras.h"
#include "coords/mmdb.h"
#include "coords/mmdb-crystal.h" //need for Bond_lines now
#include "coords/Cartesian.h"
#include "coords/Bond_lines.h"

#include "graphics-info.h"

#if __APPLE__
#   include <OpenGL/glu.h>
#else
#   include <GL/glu.h>
#endif

#include "molecule-class-info.h"

#include "coords/cos-sin.h"

#include "globjects.h"

#include "pick.h"

#include "cc-interface.hh"  // for atom_info_as_text_for_statusbar, maybe not here then

// #define DEBUG_SYMM_PICK 1

#ifdef DEBUG_SYMM_PICK
#include <ostream>
void write_symm_search_point(std::ofstream&s , const coot::Cartesian &cart) {
   s << "(" << cart.x() << " " << cart.y() << " " << cart.z() << ")\n";
}
#endif 


// 
coot::Symm_Atom_Pick_Info_t
graphics_info_t::symmetry_atom_pick() const { 

   coot::Cartesian front = unproject(0.0);
   coot::Cartesian back  = unproject(1.0);
   return symmetry_atom_pick(front, back);
}

// 
coot::Symm_Atom_Pick_Info_t
graphics_info_t::symmetry_atom_pick(const coot::Cartesian &front, const coot::Cartesian &back) const {

   
   coot::Cartesian screen_centre = RotationCentre();

   // Cartesian centre_unproj = unproject(0.5); // not needed, use midpoint of front and back.
   coot::Cartesian mid_point = front.mid_point(back);
   float dist_front_to_back = (front - back).amplitude();

#ifdef DEBUG_SYMM_PICK
   std::ofstream sps("symm-points.txt");
#endif 
   coot::Symm_Atom_Pick_Info_t p_i;
   p_i.success = GL_FALSE; // no atom found initially.
   float dist, min_dist = 0.4;

   for (int imol=0; imol<n_molecules(); imol++) {

      if (graphics_info_t::molecules[imol].atom_selection_is_pickable()) {
	 if (graphics_info_t::molecules[imol].show_symmetry) { 
      
	    if (graphics_info_t::molecules[imol].has_model()) {
	       atom_selection_container_t atom_sel = graphics_info_t::molecules[imol].atom_sel;
	       molecule_extents_t mol_extents(atom_sel, symmetry_search_radius);
	       std::vector<std::pair<symm_trans_t, Cell_Translation> > boxes =
		  mol_extents.which_boxes(screen_centre, atom_sel);

	       if (boxes.size() > 0) { // used to be 1 (which was a bug, I'm pretty sure).

		  char *spacegroup_str = atom_sel.mol->GetSpaceGroup();
		  if (!spacegroup_str) {
		     std::cout << "ERROR:: null spacegroup_str in symmetry pick\n";
		  } else {
		     clipper::Spacegroup spg;
		     clipper::Cell cell;
		     short int spacegroup_ok = 0;
		     try {
			std::pair<clipper::Cell,clipper::Spacegroup> xtal =
			   coot::util::get_cell_symm( atom_sel.mol );
			cell = xtal.first;
			spg  = xtal.second;
			spacegroup_ok = 1;
		     } catch (const std::runtime_error &except) {
			cout << "!! get_cell_symm() fails in symmetry_atom_pick"
			     << endl;
		     }
		     if (spacegroup_ok == 1) { 
			// 		     std::cout << "DEBUG:: Initing clipper::Spacegroup: "
			// 			       << spg.symbol_hm() << std::endl;

			// we want to generate all the symmetry atoms for each box
			// and find the atom with the closest approach.

			std::vector<coot::clip_hybrid_atom> hybrid_atom(atom_sel.n_selected_atoms);

#ifdef DEBUG_SYMM_PICK			   
			std::cout << "symm_atom_pick: there are " << boxes.size() << " boxes" << std::endl;
			std::cout << "Here are the boxes: " << std::endl;
			for (int ii=0; ii< boxes.size(); ii++) {
			   std::cout << ii << " " << boxes[ii].first << " " << boxes[ii].second  << std::endl;
			} 
#endif // DEBUG_SYMM_PICK			   

			for (unsigned int ii=0; ii< boxes.size(); ii++) {

			   // What are the symmetry atoms for this box?

			   // Notice that we do the unusual step of creating the
			   // vector here and modifying it via a pointer - this is
			   // for speed.
			   // 
			   fill_hybrid_atoms(&hybrid_atom, atom_sel, spg, cell, boxes[ii]);
			   int n = hybrid_atom.size();
			   // std::cout << "DEBUG:: n hybrid_atoms" << n << std::endl;
			   for(int i=0; i<n; i++) { 
			      dist = hybrid_atom[i].pos.distance_to_line(front, back);

			      // 		  if (boxes[ii].isym() == 3 && boxes[ii].x() == 0 && boxes[ii].y() == 0 && boxes[ii].z() == 0) {
			      // 			std::cout << "selected box: " << hybrid_atom[i].atom << " "
			      // 				  << hybrid_atom[i].pos << " "
			      // 				  << "scrcent " << screen_centre << std::endl;
			      // 		  }

#ifdef DEBUG_SYMM_PICK			   
			      write_symm_search_point(sps, hybrid_atom[i].pos);
#endif			   
		  
			      if (dist < min_dist) {

				 float dist_to_rotation_centre = (hybrid_atom[i].pos - screen_centre).amplitude();
				 // std::cout << "dist_to rotation_centre " << dist_to_rotation_centre << " " 
				 //  			       << hybrid_atom[i].atom << " "
				 // 			       << hybrid_atom[i].pos << " "
				 // 			       << "scrcent " << screen_centre << " "
				 //  			       << boxes[ii].isym() << " " 
				 //  			       << boxes[ii].x() << " " << boxes[ii].y() << " " 
				 //  			       << boxes[ii].z() << " "  << std::endl;

				 if (dist_to_rotation_centre < symmetry_search_radius ||
				     molecules[imol].symmetry_whole_chain_flag || molecules[imol].symmetry_as_calphas) {
				    // 			std::cout << "updating min_dist to " << dist << " " 
				    // 				  << hybrid_atom[i].atom << " "
				    // 				  << boxes[ii].isym() << " " 
				    // 				  << boxes[ii].x() << " " << boxes[ii].y() << " " 
				    // 				  << boxes[ii].z() << " "  << std::endl;

				    // how do we know we are not picking
				    // something behind the back clipping
				    // plane?  e.g. we are zoomed in or have
				    // narrow clipping planes, we don't want to
				    // pick symmetry atoms that we can't see.

				    float front_picked_atom_line_length = (front - hybrid_atom[i].pos).amplitude();
				    float  back_picked_atom_line_length = (back  - hybrid_atom[i].pos).amplitude();
			   
				    // std::cout << "comparing " << front_picked_atom_line_length << " to "
				    // << dist_front_to_back << " to " << back_picked_atom_line_length
				    // << std::endl;

				    // This test make things better, but still
				    // not right, it seem, so lets artificially
				    // slim down the limits - so that we don't
				    // get so many false positives.
			   
				    if ( (front_picked_atom_line_length < dist_front_to_back) &&
					 ( back_picked_atom_line_length < dist_front_to_back)
					 ) { 
				       min_dist = dist;
				       p_i.success = GL_TRUE;
				       p_i.hybrid_atom = hybrid_atom[i];
				       // and these for labelling:
				       p_i.atom_index = i;
				       p_i.imol = imol; 
				       p_i.symm_trans = boxes[ii].first;
				       p_i.pre_shift_to_origin = boxes[ii].second;
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
      }
   }

   if (p_i.success == GL_TRUE) {
      
      // make sure that we are looking at the molecule that had
      // the nearest contact:
      // 
      atom_selection_container_t SelAtom = graphics_info_t::molecules[p_i.imol].atom_sel; 
      
      std::string alt_conf_bit("");
      mmdb::Atom *at = molecules[p_i.imol].atom_sel.atom_selection[p_i.atom_index];
      if (strncmp(at->altLoc, "", 1))
	 alt_conf_bit=std::string(",") + std::string(at->altLoc);
      
      std::cout << "(" << p_i.imol << ") " 
		<< SelAtom.atom_selection[p_i.atom_index]->name 
		<< alt_conf_bit << "/"
		<< SelAtom.atom_selection[p_i.atom_index]->GetModelNum()
		<< "/chainid=\""
		<< SelAtom.atom_selection[p_i.atom_index]->GetChainID()  << "\"/"
		<< SelAtom.atom_selection[p_i.atom_index]->GetSeqNum()   
		<< SelAtom.atom_selection[p_i.atom_index]->GetInsCode()  << "/"
		<< SelAtom.atom_selection[p_i.atom_index]->GetResName()
		<< ", occ: " 
		<< SelAtom.atom_selection[p_i.atom_index]->occupancy 
		<< " with B-factor: "
		<< SelAtom.atom_selection[p_i.atom_index]->tempFactor
		<< " element: "
		<< SelAtom.atom_selection[p_i.atom_index]->element
		<< " at " << translate_atom(SelAtom, p_i.atom_index, p_i.symm_trans) << std::endl;
   }

   if (p_i.success) {
      std::pair<symm_trans_t, Cell_Translation> symm_trans(p_i.symm_trans, p_i.pre_shift_to_origin);      
      coot::Cartesian tpos = translate_atom_with_pre_shift(molecules[p_i.imol].atom_sel,
							   p_i.atom_index,
							   symm_trans);

      std::string ai;
      int expanded_flag = 1;
      ai = atom_info_as_text_for_statusbar(p_i.atom_index, p_i.imol, symm_trans);
      add_status_bar_text(ai);
      //gtk_statusbar_push(GTK_STATUSBAR(graphics_info_t::statusbar),
      //		 graphics_info_t::statusbar_context_id,
      //		 ai.c_str());
   }
   
//    std::cout << "Picked: status: " << p_i.success << " index: "<< p_i.atom_index
//  	     << " " << p_i.symm_trans << " "
//  	     << p_i.pre_shift_to_origin << std::endl;
   
   return p_i;
}

// This presumes that the mmdb symop order is the same as the clipper
// symop order [no longer true 12Nov2003]
// 
// private 
//
// So, to recap what happened here: It turns out, that unlike I
// expected, the clipper::Spacegroup::symops() don't come in the same
// order as those of mmdb.  So we have to construct the symop (really
// an RTop_orth here) from the mmdb matrix and transform the
// coordinates with this rtop.  A bit ugly, eh? 
// 
// But on the plus side, it works, no mind-boggling memory allocation
// problems and it's quite fast compared to the old version.
//
void 
graphics_info_t::fill_hybrid_atoms(std::vector<coot::clip_hybrid_atom> *hybrid_atoms, 
				   const atom_selection_container_t &asc,
				   const clipper::Spacegroup &spg, 
				   const clipper::Cell &cell,
				   const std::pair<symm_trans_t, Cell_Translation> &symm_trans) const { 

   clipper::Coord_orth trans_pos; 
   clipper::Coord_orth co; 

   // if (symm_trans.isym() == 3 && symm_trans.x() == 0 && symm_trans.y() == 0 && symm_trans.z() == 0) {
   // this does not return what I expect!  The [3] of P 21 21 21 is
   // 1/2+X,1/2-Y,-Z ,but this returns: -x+1/2, -y, z+1/2, the
   // [1] element of the [zero indexed] symops list.
   // 
   // std::cout << "selected symm_trans: " << spg.symop(symm_trans.isym()).format() << std::endl;

   // so let's create a symop from symm_trans.isym() by creating an rtop_frac
   // 
   // Either that or make which_box return clipper symop.  Yes, I like that.
   // 
   mmdb::mat44 my_matt;
   mmdb::mat44 mol_to_origin_mat;

   asc.mol->GetTMatrix(mol_to_origin_mat, 0, 
		       -symm_trans.second.us,
		       -symm_trans.second.vs,
		       -symm_trans.second.ws);
      
   asc.mol->GetTMatrix(my_matt,
		       symm_trans.first.isym(),
		       symm_trans.first.x(),
		       symm_trans.first.y(),
		       symm_trans.first.z());
      
   //       std::cout << "my_matt is: " << std::endl
   // 		<< my_matt[0][0] << " "  << my_matt[0][1] << " "
   // 		<< my_matt[0][2] << " "  << my_matt[0][3] << " "  << std::endl
   // 		<< my_matt[1][0] << " "  << my_matt[1][1] << " "
   // 		<< my_matt[1][2] << " "  << my_matt[1][3] << " "  << std::endl
   // 		<< my_matt[2][0] << " "  << my_matt[2][1] << " "
   // 		<< my_matt[2][2] << " "  << my_matt[2][3] << " "  << std::endl
   // 		<< my_matt[3][0] << " "  << my_matt[3][1] << " "
   // 		<< my_matt[3][2] << " "  << my_matt[3][3] << " "  << std::endl; 
      
   clipper::Mat33<double> clipper_mol_to_ori_mat(mol_to_origin_mat[0][0], mol_to_origin_mat[0][1], mol_to_origin_mat[0][2],
						 mol_to_origin_mat[1][0], mol_to_origin_mat[1][1], mol_to_origin_mat[1][2],
						 mol_to_origin_mat[2][0], mol_to_origin_mat[2][1], mol_to_origin_mat[2][2]);
   clipper::Mat33<double> clipper_mat(my_matt[0][0], my_matt[0][1], my_matt[0][2],
				      my_matt[1][0], my_matt[1][1], my_matt[1][2],
				      my_matt[2][0], my_matt[2][1], my_matt[2][2]);
   clipper::Coord_orth  cco(my_matt[0][3], my_matt[1][3], my_matt[2][3]);
   clipper::Coord_orth  cco_mol_to_ori(mol_to_origin_mat[0][3], mol_to_origin_mat[1][3], mol_to_origin_mat[2][3]);
   clipper::RTop_orth rtop_mol_to_ori(clipper_mol_to_ori_mat, cco_mol_to_ori);
   clipper::RTop_orth rtop(clipper_mat, cco);
      
   for(int i=0; i<asc.n_selected_atoms; i++) {
      co = clipper::Coord_orth(asc.atom_selection[i]->x, 
			       asc.atom_selection[i]->y, 
			       asc.atom_selection[i]->z);
      trans_pos = co.transform(rtop_mol_to_ori).transform(rtop);

//       std::cout << i << "  " << cf.format() << "  " << trans_pos_cf.format() 
// 		<< "  " << spg.symop(symm_trans.isym()).format() << std::endl;
      (*hybrid_atoms)[i] = coot::clip_hybrid_atom(asc.atom_selection[i], 
						  coot::Cartesian(trans_pos.x(),
								  trans_pos.y(),
								  trans_pos.z()));
   }
}


// pickable moving atoms molecule
//
pick_info
graphics_info_t::moving_atoms_atom_pick(short int pick_mode) const {

   pick_info p_i;
   p_i.min_dist = 0; // keep compiler happy
   p_i.atom_index = -1; // ditto
   p_i.imol = -1; // ditto
   p_i.success = GL_FALSE;
   coot::Cartesian front = unproject(0.0);
   coot::Cartesian back  = unproject(1.0);

   float m_front = 0.8;   // pickable distance
   float m_back =  0.04;

   std::pair<float, float> min_dist(999,999);
   float close_score_best = 999.9;

   float d_front_to_back = (back-front).amplitude();

   // This is the signal that moving_atoms_asc is clear
   if (moving_atoms_asc->n_selected_atoms > 0) {

      for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
	 mmdb::Atom *at = moving_atoms_asc->atom_selection[i];
	 coot::Cartesian atom_pos(at->x, at->y, at->z);
	 //
	 bool at_is_hydrogen = coot::is_hydrogen_p(at);

	 if ((! at_is_hydrogen && pick_mode == PICK_ATOM_NON_HYDROGEN) ||
	     (pick_mode == PICK_ATOM_ALL_ATOM)) {

	    if (atom_pos.within_box(front,back)) {
	       float dist = atom_pos.distance_to_line(front, back);
	       float d_atom_to_front = (atom_pos-front).amplitude();

	       float frac = d_atom_to_front/d_front_to_back;

	       // we are are looking for a low score, so the
	       // front has a lower score than the back
	       float m = m_back + frac * (m_front-m_back);

	       float m_limit = m_back + (1.0-frac) * (m_front-m_back);

	       float close_score = dist * m;

	       if (false) // debugging algorithm
		  std::cout << i << " " << coot::atom_spec_t(at) << " "
			    << close_score_best << " "
			    << close_score << " dist " << dist
			    << " d-atom-front-to-back " << d_atom_to_front
			    << " m " << m << " m_limit " << m_limit
			    << " frac " << frac
			    << std::endl;

	       if (dist < m_limit) {
		  if (close_score < close_score_best) {

		     close_score_best = close_score;
		     p_i.success = GL_TRUE;
		     p_i.atom_index = i;
		  }
	       }
	    }
	 }
      }
   }
   return p_i;
}

// examines the imol_moving_atoms molecule for correspondence
// 
// static
bool
graphics_info_t::fixed_atom_for_refinement_p(mmdb::Atom *at) {

   bool r = 0; 
   if (is_valid_model_molecule(imol_moving_atoms)) {
      std::vector<coot::atom_spec_t> fixed = molecules[imol_moving_atoms].get_fixed_atoms();
      for (unsigned int ifixed=0; ifixed<fixed.size(); ifixed++) {
	 if (fixed[ifixed].matches_spec(at)) {
// 	    std::cout << " fixed_atom_for_refinement_p found a matcher "
// 		      << fixed[ifixed] << std::endl;
	    r = 1;
	    break;
	 } 
      }
   }
   return r;
} 


// Setup moving atom-drag if we are.
// (we are acting on button1 down)
// note to self, maybe you wanted this (check_if_moving_atom_pull()), not pick_intermediate_atom()?
//
void
graphics_info_t::check_if_moving_atom_pull(bool was_a_double_click) {

   short int atom_pick_mode = PICK_ATOM_ALL_ATOM;
   if (! moving_atoms_have_hydrogens_displayed)
      atom_pick_mode = PICK_ATOM_NON_HYDROGEN;
   pick_info pi = moving_atoms_atom_pick(atom_pick_mode);
   if (pi.success == GL_TRUE) {

      // Flash picked atom.
      mmdb::Atom *at = moving_atoms_asc->atom_selection[pi.atom_index];
      clipper::Coord_orth co(at->x, at->y, at->z);
      graphics_info_t::flash_position(co);

      // quite possibly we will have success if moving_atoms_asc is
      // not null...

      // setup for move (pull) of that atom then (and later we will
      // move the other atoms in some sort of stretch system)
      //
      // Recall that we came here from a button click event, no a
      // motion event - this is a setup function - not a "doit"
      // function
      //
      // This is a bit tedious, because we have to search in the
      // imol_moving_atoms molecule.
      //

      // 20180219. Now we have many.

      moving_atoms_currently_dragged_atom_index = pi.atom_index;

      if (! was_a_double_click) {
	 //
	 moving_atoms_dragged_atom_indices.insert(pi.atom_index);

	 if (false)
	    std::cout << "moving_atoms_currently_dragged_atom_index " << moving_atoms_currently_dragged_atom_index
		      << std::endl;

	 in_moving_atoms_drag_atom_mode_flag = 1;

         // This is no longer relevant I think
         //
	 // Find the fixed_points_sheared_drag from the imol_moving_atoms.
	 // and their flag, have_fixed_points_sheared_drag_flag
	 // set_fixed_points_for_sheared_drag(); // superceeded
	 //
	 // if (drag_refine_idle_function_token != -1) {
	 //    // 	 std::cout << "DEBUG:: removing token " << drag_refine_idle_function_token
	 //    // 		   << std::endl;
	 //    gtk_idle_remove(drag_refine_idle_function_token);
	 //    drag_refine_idle_function_token = -1; // magic "not in use" value
	 // }

      } else {
	 // a double click on a dragged intermediate atom means "remove this pull restraint"

	 std::set<int>::const_iterator it = moving_atoms_dragged_atom_indices.find(pi.atom_index);
	 if (it != moving_atoms_dragged_atom_indices.end()) {
	    std::cout << "DEBUG:: erasing moving atoms dragged atom with index: " << *it << std::endl;
	    moving_atoms_dragged_atom_indices.erase(it);
	    if (*it < moving_atoms_asc->n_selected_atoms) {
	       mmdb::Atom *at = moving_atoms_asc->atom_selection[*it];
	       if (at) {
		  coot::atom_spec_t spec(at);
		  atom_pull_off(spec);
		  clear_atom_pull_restraint(spec, true); // refine again
	       }
	    } else {
	       std::cout << "ERROR:: dragged atom out of range " << *it << " "
			 << moving_atoms_asc->n_selected_atoms << std::endl;
	    }
	 }
      }

   } else {
      in_moving_atoms_drag_atom_mode_flag = 0;
   }
}


// static
void
graphics_info_t::atom_pull_off(const coot::atom_spec_t &spec) {
    for (std::size_t i=0; i<atom_pulls.size(); i++) {
       if (atom_pulls[i].spec == spec)
          atom_pulls[i].off();
    }
}

// static
void
graphics_info_t::atom_pulls_off(const std::vector<coot::atom_spec_t> &specs) {
    for (std::size_t j=0; j<specs.size(); j++)
       for (std::size_t i=0; i<atom_pulls.size(); i++)
          if (atom_pulls[i].spec == specs[j])
             atom_pulls[i].off();
}

// static
void
graphics_info_t::all_atom_pulls_off() {
     for (std::size_t i=0; i<atom_pulls.size(); i++)
       atom_pulls[i].off();
     atom_pulls.clear();
}

void
graphics_info_t::set_fixed_points_for_sheared_drag() {

}


// 
void
graphics_info_t::move_moving_atoms_by_shear(int screenx, int screeny,
					    short int squared_flag) {
   
   // First we have to find the "fixed" points connected to them most
   // distant from the moving_atoms_dragged_atom.
   //
   // So now we ask, how much should we move the moving atom?  We have
   // the current mouse position: x, y and the previous mouse
   // position: GetMouseBeginX(), GetMouseBeginY(). Let's use
   // unproject to find the shift in world coordinates...
   //
   coot::Cartesian old_mouse_real_world = unproject_xyz(int(GetMouseBeginX()),
							int(GetMouseBeginY()),
							0.5);
   coot::Cartesian current_mouse_real_world = unproject_xyz(screenx, screeny, 0.5);
   
   coot::Cartesian diff = current_mouse_real_world - old_mouse_real_world;
   
   // now tinker with the moving atoms coordinates...

   // 20180219 - I can't be bothered to fix this non-used code

   /*

   if (moving_atoms_dragged_atom_index >= 0) {
      if (moving_atoms_dragged_atom_index < moving_atoms_asc->n_selected_atoms) {
	 mmdb::Atom *at = moving_atoms_asc->atom_selection[moving_atoms_dragged_atom_index];
	 if (at) {
	    
	    move_moving_atoms_by_shear_internal(diff, squared_flag);
	    
	    // and regenerate the bonds of the moving atoms:
	    // 
	    int do_disulphide_flag = 0;
	    Bond_lines_container bonds(*moving_atoms_asc, do_disulphide_flag);
	    regularize_object_bonds_box.clear_up();
	    regularize_object_bonds_box = bonds.make_graphical_bonds();
	    graphics_draw();
	 } else {
	    std::cout << "ERROR: null atom in move_moving_atoms_by_shear\n";
	 }
      } else {
	 std::cout << "ERROR: out of index (over) in move_moving_atoms_by_shear\n";
      }
   } else {
      std::cout << "ERROR: out of index (under) in move_moving_atoms_by_shear\n";
   }
   */
}


// For rotate/translate moving atoms dragged movement
// 
void
graphics_info_t::move_moving_atoms_by_simple_translation(int screenx, int screeny) {


   coot::Cartesian old_mouse_real_world = unproject_xyz(int(GetMouseBeginX()),
						  int(GetMouseBeginY()),
						  0.5);
   coot::Cartesian current_mouse_real_world = unproject_xyz(screenx, screeny, 0.5);

   coot::Cartesian diff = current_mouse_real_world - old_mouse_real_world;
   mmdb::Atom *at;
   for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
      at = moving_atoms_asc->atom_selection[i];
      at->x += diff.x();
      at->y += diff.y();
      at->z += diff.z();
   }
   int do_disulphide_flag = 0;

   if (molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS ||
       molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS_PLUS_LIGANDS ||
       molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS_PLUS_LIGANDS_AND_SIDECHAINS ||
       molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR ||
       molecules[imol_moving_atoms].Bonds_box_type() == coot::COLOUR_BY_RAINBOW_BONDS) {
      
      Bond_lines_container bonds;
      bonds.do_Ca_plus_ligands_bonds(*moving_atoms_asc, imol_moving_atoms, Geom_p(), 1.0, 4.7,
                                     draw_missing_loops_flag, false);
      regularize_object_bonds_box.clear_up();
      regularize_object_bonds_box = bonds.make_graphical_bonds();
   } else {
      Bond_lines_container bonds(*moving_atoms_asc, do_disulphide_flag);
      regularize_object_bonds_box.clear_up();
      regularize_object_bonds_box = bonds.make_graphical_bonds();
   }
   graphics_draw();
}

// uses moving_atoms_dragged_atom_index
void
graphics_info_t::move_single_atom_of_moving_atoms(int screenx, int screeny) {

   
   coot::Cartesian old_mouse_real_world = unproject_xyz(int(GetMouseBeginX()),
							int(GetMouseBeginY()), 0.5);
   coot::Cartesian current_mouse_real_world = unproject_xyz(screenx, screeny, 0.5);

   coot::Cartesian diff = current_mouse_real_world - old_mouse_real_world;
   mmdb::Atom *at = moving_atoms_asc->atom_selection[moving_atoms_currently_dragged_atom_index];
   at->x += diff.x();
   at->y += diff.y();
   at->z += diff.z();

   int do_disulphide_flag = 0;
   Bond_lines_container bonds(*moving_atoms_asc, do_disulphide_flag);
   regularize_object_bonds_box.clear_up();
   regularize_object_bonds_box = bonds.make_graphical_bonds();
   graphics_draw();

}

void
graphics_info_t::move_atom_pull_target_position(int screen_x, int screen_y) {

   coot::Cartesian front = unproject_xyz(screen_x, screen_y, 0);
   coot::Cartesian back  = unproject_xyz(screen_x, screen_y, 1);

   mmdb::Atom *at = moving_atoms_asc->atom_selection[moving_atoms_currently_dragged_atom_index];
   coot::Cartesian c_at(at->x, at->y, at->z);
   coot::Cartesian fb = back - front;
   coot::Cartesian  h = c_at - front;

   float d_fb = fb.length();
   float d_h  =  h.length();

   double cos_theta = coot::dot_product(fb, h)/(d_fb*d_h);
   double d = d_h * cos_theta;
   float z_depth = d/d_fb;

   if (z_depth > 1.0) z_depth = 1.0;
   if (z_depth < 0.0) z_depth = 0.0;

   coot::Cartesian current_mouse_real_world = unproject_xyz(screen_x, screen_y, z_depth);
   clipper::Coord_orth c_pos(current_mouse_real_world.x(),
			     current_mouse_real_world.y(),
			     current_mouse_real_world.z());

   atom_pull_info_t atom_pull_local = atom_pull_info_t(coot::atom_spec_t(at), c_pos);

   if (false)
      std::cout << "graphics_info_t::move_atom_pull_target_position() atom_pull.status: "
		<< atom_pull_local.get_status() << " "
		<< " atom pull:: idx " << moving_atoms_currently_dragged_atom_index << " "
		<< coot::co(at).format() << " to " << current_mouse_real_world
		<< std::endl;

   continue_threaded_refinement_loop = false;
   while (restraints_lock) {
        std::this_thread::sleep_for(std::chrono::milliseconds(2)); // not sure about the delay
   }

   // add to or replace in atom_pulls (for representation)
   add_or_replace_current(atom_pull_local);
   //
   last_restraints->add_atom_pull_restraint(atom_pull_local.spec, c_pos);
   thread_for_refinement_loop_threaded();
   graphics_draw();  // needed
}

void graphics_info_t::add_or_replace_current(const atom_pull_info_t &atom_pull_in) {

   bool done = false;

   std::vector<atom_pull_info_t>::iterator it;
   for(it=atom_pulls.begin(); it!=atom_pulls.end(); it++) {
      if (it->spec == atom_pull_in.spec) {
	 it->pos = atom_pull_in.pos;
	 it->on(); // do do do be do... turn it o-o-o-on
	 done = true;
	 break;
      }
   }

   if (! done) {
      // std::cout << "Adding to atom_pulls: " << atom_pull_in.spec << " " << atom_pull_in.pos.format() << std::endl;
      atom_pulls.push_back(atom_pull_in);
   }

}


// we need a vector version of this.
void
graphics_info_t::add_target_position_restraint_for_intermediate_atom(const coot::atom_spec_t &spec,
								     const clipper::Coord_orth &target_pos) {


   // get the restraints lock before adding these

   get_restraints_lock(__FUNCTION__);

   atom_pull_info_t atom_pull_local = atom_pull_info_t(spec, target_pos);
   add_or_replace_current(atom_pull_local);
   if (last_restraints) {
      last_restraints->add_atom_pull_restraint(spec, target_pos);
   }
   release_restraints_lock(__FUNCTION__);

   thread_for_refinement_loop_threaded();

}

// vector version of above
void
graphics_info_t::add_target_position_restraints_for_intermediate_atoms(const std::vector<std::pair<coot::atom_spec_t, clipper::Coord_orth> > &atom_spec_position_vec) {

   // get the restraints lock before adding this

   if (last_restraints) {

      get_restraints_lock(__FUNCTION__);
      for (std::size_t i=0; i<atom_spec_position_vec.size(); i++) {
	 coot::atom_spec_t spec  = atom_spec_position_vec[i].first;
	 clipper::Coord_orth pos = atom_spec_position_vec[i].second;
	 atom_pull_info_t atom_pull_local = atom_pull_info_t(spec, pos);
	 add_or_replace_current(atom_pull_local);
	 last_restraints->add_atom_pull_restraint(spec, pos);
      }

      // unlock restraints and start refinement again
      release_restraints_lock(__FUNCTION__);
      thread_for_refinement_loop_threaded();

   } else {
      std::cout << "WARNING:: in add_target_position_restraints_for_intermediate_atoms() no restraints"
		<< std::endl;
   }

}




// diff_std is the difference in position of the moving atoms, the
// other atoms of the moving_atoms_asc are moved relative to it by
// different amounts...
//
// This function currently seems to be called with linear_movement_scaling_flag as false.
// 
void
graphics_info_t::move_moving_atoms_by_shear_internal(const coot::Cartesian &diff_std, 
                                                     short int linear_movement_scaling_flag) {

   coot::Cartesian diff = diff_std;
   mmdb::Atom *mat = moving_atoms_asc->atom_selection[moving_atoms_currently_dragged_atom_index];
   coot::Cartesian moving_atom(mat->x, mat->y, mat->z);
   float d_to_moving_at_max = -9999999.9;
   int d_array_size = moving_atoms_asc->n_selected_atoms;
   float *d_to_moving_at = new float[d_array_size];
   float frac;
   mmdb::Atom *at;
   
   for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {

      at = moving_atoms_asc->atom_selection[i];
      coot::Cartesian this_atom(at->x, at->y, at->z);
      d_to_moving_at[i] = (this_atom - moving_atom).length();
      if (d_to_moving_at[i] > d_to_moving_at_max) {
	 d_to_moving_at_max = d_to_moving_at[i];
      }
   }

   double dr;
   for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {

      at = moving_atoms_asc->atom_selection[i];
      coot::Cartesian atom_pt(at->x, at->y, at->z);
      if (linear_movement_scaling_flag == 0) {
	 // 	 frac = (1.0 - d_to_moving_at[i]/d_to_moving_at_max); old
	 dr = d_to_moving_at[i]/d_to_moving_at_max;
	 frac = (1-pow(dr, refinement_drag_elasticity)); // 0.1 by default
      } else {
	 dr = d_to_moving_at[i]/d_to_moving_at_max;
	 // frac = (1.0 - dr)*(1.0 - dr);
	 frac = (1.0 - dr);
      }
	    
      if (! fixed_atom_for_refinement_p(at)) { 
	       
	 at->x += frac*diff.x();
	 at->y += frac*diff.y();
	 at->z += frac*diff.z();
      }
   }
   delete [] d_to_moving_at;
}

void
graphics_info_t::do_post_drag_refinement_maybe() {

#ifdef HAVE_GSL

   // std::cout << "Here in do_post_drag_refinement_maybe() with last_restraints_size() "
   // << last_restraints_size() << std::endl;

   if (last_restraints_size() > 0) {
       thread_for_refinement_loop_threaded();
    } else {
       std::cout << "DEBUG:: not doing refinement - no restraints." << std::endl;
    }

#endif // HAVE_GSL   
}


// static
void
graphics_info_t::statusbar_ctrl_key_info() { // Ctrl to rotate or pick?

   graphics_info_t g;
   if (graphics_info_t::control_key_for_rotate_flag) {
      g.add_status_bar_text("Use Ctrl Left-mouse to rotate the view.");
   } else {
      g.add_status_bar_text("Use Ctrl Left-mouse to pick an atom...");
   }
}

clipper::Coord_orth
graphics_info_t::moving_atoms_centre() const {

   clipper::Coord_orth moving_middle(0,0,0);

   // Let's find the middle of the moving atoms and set
   // rotation_centre to that:
   int n = moving_atoms_asc->n_selected_atoms;
   if (n > 0) {
      float sum_x = 0.0; float sum_y = 0.0; float sum_z = 0.0;
      for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
	 sum_x += moving_atoms_asc->atom_selection[i]->x;
	 sum_y += moving_atoms_asc->atom_selection[i]->y;
	 sum_z += moving_atoms_asc->atom_selection[i]->z;
      }
      moving_middle = clipper::Coord_orth(sum_x/float(n), sum_y/float(n), sum_z/float(n));
   }
   return moving_middle;
} 

								       
// Presumes that rotation centre can be got from mmdb::Atom *rot_trans_rotation_origin_atom;
// 
void
graphics_info_t::rotate_intermediate_atoms_round_screen_z(double angle) {

   if (rot_trans_rotation_origin_atom) { 
      if (moving_atoms_asc->mol) {
	 if (moving_atoms_asc->n_selected_atoms > 0) {
	    coot::Cartesian front  = unproject_xyz(0, 0, 0.0);
	    coot::Cartesian centre = unproject_xyz(0, 0, 0.5);
	    coot::Cartesian screen_z = (front - centre);

	    clipper::Coord_orth screen_vector =  clipper::Coord_orth(screen_z.x(), 
								     screen_z.y(), 
								     screen_z.z());
	    mmdb::Atom *rot_centre = rot_trans_rotation_origin_atom;
	    clipper::Coord_orth rotation_centre(rot_centre->x, 
						rot_centre->y, 
						rot_centre->z);
	    // But! maybe we have a different rotation centre
	    if (rot_trans_zone_rotates_about_zone_centre) {
	       if (moving_atoms_asc->n_selected_atoms  > 0) {
		  rotation_centre = moving_atoms_centre();
	       }
	    }
	 
	    for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
	       clipper::Coord_orth co(moving_atoms_asc->atom_selection[i]->x,
				      moving_atoms_asc->atom_selection[i]->y,
				      moving_atoms_asc->atom_selection[i]->z);
	       clipper::Coord_orth new_pos = 
		  coot::util::rotate_around_vector(screen_vector, co, rotation_centre, angle);
	       moving_atoms_asc->atom_selection[i]->x = new_pos.x();
	       moving_atoms_asc->atom_selection[i]->y = new_pos.y();
	       moving_atoms_asc->atom_selection[i]->z = new_pos.z();
	    }
	    int do_disulphide_flag = 0;

	    if (molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS ||
		molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS_PLUS_LIGANDS ||
		molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS_PLUS_LIGANDS_AND_SIDECHAINS ||
		molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR ||
		molecules[imol_moving_atoms].Bonds_box_type() == coot::COLOUR_BY_RAINBOW_BONDS) {
	       
	       Bond_lines_container bonds;
	       bonds.do_Ca_plus_ligands_bonds(*moving_atoms_asc, imol_moving_atoms, Geom_p(), 1.0, 4.7,
                                              draw_missing_loops_flag, false);
	       regularize_object_bonds_box.clear_up();
	       regularize_object_bonds_box = bonds.make_graphical_bonds();
	    } else {
	       Bond_lines_container bonds(*moving_atoms_asc, do_disulphide_flag);
	       regularize_object_bonds_box.clear_up();
	       regularize_object_bonds_box = bonds.make_graphical_bonds();
	    }
	    graphics_draw();
	 }
      }
   }
} 

// Presumes that rotation centre can be got from mmdb::Atom *rot_trans_rotation_origin_atom;
// 
void
graphics_info_t::rotate_intermediate_atoms_round_screen_x(double angle) {

   if (rot_trans_rotation_origin_atom) { 
      if (moving_atoms_asc->mol) {
	 if (moving_atoms_asc->n_selected_atoms > 0) {
	    coot::Cartesian front  = unproject_xyz(0, 0, 0.0);
	    coot::Cartesian centre = unproject_xyz(0, 0, 0.5);
	    coot::Cartesian right  = unproject_xyz(1, 0, 0.5);
	    coot::Cartesian screen_z = (front - centre);
	    coot::Cartesian screen_x = (right - centre);

	    clipper::Coord_orth screen_vector =  clipper::Coord_orth(screen_x.x(), 
								     screen_x.y(), 
								     screen_x.z());
	    mmdb::Atom *rot_centre = rot_trans_rotation_origin_atom;
	    clipper::Coord_orth rotation_centre(rot_centre->x, 
						rot_centre->y, 
						rot_centre->z);
	 
	    // But! maybe we have a different rotation centre
	    if (rot_trans_zone_rotates_about_zone_centre)
	       if (moving_atoms_asc->n_selected_atoms  > 0)
		  rotation_centre = moving_atoms_centre();
	    
	    for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
	       clipper::Coord_orth co(moving_atoms_asc->atom_selection[i]->x,
				      moving_atoms_asc->atom_selection[i]->y,
				      moving_atoms_asc->atom_selection[i]->z);
	       clipper::Coord_orth new_pos = 
		  coot::util::rotate_around_vector(screen_vector, co, rotation_centre, angle);
	       moving_atoms_asc->atom_selection[i]->x = new_pos.x();
	       moving_atoms_asc->atom_selection[i]->y = new_pos.y();
	       moving_atoms_asc->atom_selection[i]->z = new_pos.z();
	    }
	    int do_disulphide_flag = 0;
	    if (molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS ||
		molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS_PLUS_LIGANDS ||
		molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR ||
		molecules[imol_moving_atoms].Bonds_box_type() == coot::COLOUR_BY_RAINBOW_BONDS) {
	       
	       Bond_lines_container bonds;
	       bonds.do_Ca_plus_ligands_bonds(*moving_atoms_asc, imol_moving_atoms, Geom_p(), 1.0, 4.7,
                                              draw_missing_loops_flag, false);
	       regularize_object_bonds_box.clear_up();
	       regularize_object_bonds_box = bonds.make_graphical_bonds();
	    } else {
	       Bond_lines_container bonds(*moving_atoms_asc, do_disulphide_flag);
	       regularize_object_bonds_box.clear_up();
	       regularize_object_bonds_box = bonds.make_graphical_bonds();
	    }
	    graphics_draw();
	 }
      }
   }
} 



int graphics_info_t::move_reference_chain_to_symm_chain_position() {

   int r = 0;
   if (use_graphics_interface_flag) { 
      int iw = graphics_info_t::glarea->allocation.width;
      int ih = graphics_info_t::glarea->allocation.height;
      coot::Cartesian front = unproject_xyz(iw/2, ih/2, 0);
      coot::Cartesian back  = unproject_xyz(iw/2, ih/2, 1);
      coot::Symm_Atom_Pick_Info_t naii = symmetry_atom_pick(front, back);
      if (naii.success == GL_TRUE) {
	 if (is_valid_model_molecule(naii.imol)) {
	    graphics_info_t::molecules[naii.imol].move_reference_chain_to_symm_chain_position(naii);
	    graphics_draw();
	 } else {
	    std::cout << "not valid mol" << std::endl;
	 } 
      } else {
	 std::cout << "bad pick " << std::endl;
	 std::string s = "Symm Atom not found at centre.  Are you centred on a symm atom?";
	 add_status_bar_text(s);
	 gdk_beep();
      }
   }
   return r;
}

