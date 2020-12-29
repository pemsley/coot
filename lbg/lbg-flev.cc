/* lbg/lbg.cc
 * 
 * Copyright 2010, 2011, 2012 by The University of Oxford
 * Copyright 2012, 2013, 2014, 2015, 2016 by Medical Research Council
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

#ifdef HAVE_GOOCANVAS

#ifdef USE_PYTHON
#include <Python.h> // this is here get round header warnings
#endif

#include <sys/types.h>  // for stating
#include <sys/stat.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdexcept>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include "compat/coot-sysdep.h"

#include <cairo.h>
#if CAIRO_HAS_PDF_SURFACE
#include <cairo-pdf.h>
#endif
#include <cairo-svg.h>
#include "lbg.hh"
#include "lbg-drag-and-drop.hh"
#include "qed-interface.hh" // interface to silicos-it biscu-it python function

#include <gdk/gdkkeysyms.h> // for keyboarding.

void
lbg_info_t::read_draw_residues(const std::string &file_name) {

   bool show_dynamics = 0;

   // read residues need to happen after the ligand has been placed on the canvas
   // i.e. mol.input_coords_to_canvas_coords() is called in read_residues();
   // 
   std::vector<lbg_info_t::residue_circle_t> file_residue_circles = read_residues(file_name);

   double max_dist_water_to_ligand_atom  = 3.2; // don't draw waters that are far from ligand
   double max_dist_water_to_protein_atom = 3.2; // don't draw waters that are not somehow 
                                                // attached to the protein.

   residue_circles = filter_residue_waters(file_residue_circles,
					   max_dist_water_to_ligand_atom,
					   max_dist_water_to_protein_atom);

   std::vector<int> primary_indices = get_primary_indices();

   initial_residues_circles_layout(); // twiddle residue_circles
   
   // cycles of optimisation
   // 
   std::vector<residue_circle_t> current_circles = residue_circles;
   std::vector<int> empty_vector; // for add_rep_handles, can't do it from a file.
   if (show_dynamics)
      draw_residue_circles(current_circles, empty_vector);

   // For debugging the minimizer (vs original positions).
   // 
   current_circles = offset_residues_from_orig_positions();
   
   for (int iround=0; iround<120; iround++) {
      std::cout << ":::::::::::::::::::: minimization round " << iround
		<< " :::::::::::::::::::::::::::::::::::::::::::\n";
      std::pair<int, std::vector<lbg_info_t::residue_circle_t> > new_c =
	 optimise_residue_circle_positions(residue_circles, current_circles, primary_indices);
      current_circles = new_c.second;
      if (show_dynamics)
	 draw_residue_circles(current_circles, empty_vector);
      if (new_c.first == GSL_ENOPROG)
	 break;
      if (new_c.first == GSL_SUCCESS) { 
	 break;
      }
   }

   // save the class member data
   residue_circles = current_circles;

   // has the current solution problems due to residues too close to the ligand?
   std::pair<bool, std::vector<int> > problem_status = solution_has_problems_p();
   // std::cout << "::::::::: problem status: " << problem_status.first << std::endl;
   for (unsigned int ip=0; ip<problem_status.second.size(); ip++) {
      std::cout << ":::::::::::: "
		<< residue_circles[problem_status.second[ip]].residue_label << " "
		<< residue_circles[problem_status.second[ip]].residue_type << " "
		<< residue_circles[problem_status.second[ip]].pos
		<< std::endl;
   }
   
   if (problem_status.second.size()) {
      
      // fiddle with residue_circles and reoptimise.
      // 
      reposition_problematics_and_reoptimise(problem_status.second,
					     primary_indices);
   }
   
   draw_all_flev_annotations(); // contour debugging.
}

std::vector<lbg_info_t::residue_circle_t>
lbg_info_t::filter_residue_waters(const std::vector<lbg_info_t::residue_circle_t> &r_in,
				  double max_dist_water_to_ligand_atom,
				  double max_dist_water_to_protein_atom) const {

   std::vector<lbg_info_t::residue_circle_t> v;

   for (unsigned int i=0; i<r_in.size(); i++) { 
      if (r_in[i].residue_type != "HOH") {
	 v.push_back(r_in[i]);
      } else { 
	 if (r_in[i].water_dist_to_protein > max_dist_water_to_protein_atom) {
	    // pass
	 } else {
	    bool all_too_long_distance = 1;
	    std::vector<lbg_info_t::bond_to_ligand_t> sasisfactory_bonds;
	    for (unsigned int j=0; j<r_in[i].bonds_to_ligand.size(); j++) {
	       if (r_in[i].bonds_to_ligand[j].bond_length < max_dist_water_to_ligand_atom) { 
		  all_too_long_distance = 0;
		  sasisfactory_bonds.push_back(r_in[i].bonds_to_ligand[j]);
	       }
	    }
	    if (all_too_long_distance) {
	       // pass
	    } else {
	       // replace input waters with sasisfactory_bonds.
	       // 
	       lbg_info_t::residue_circle_t n = r_in[i];
	       n.bonds_to_ligand = sasisfactory_bonds;
	       v.push_back(n);
	    } 
	 } 
      }
   }

   return v;
} 


// fiddle with the position of some residues in residue_circles
//
void
lbg_info_t::reposition_problematics_and_reoptimise(const std::vector<int> &problematics,
						   const std::vector<int> &primary_indices) {

   std::pair<lig_build::pos_t, lig_build::pos_t> l_e_pair =
      mol.ligand_extents();
   lbg_info_t::ligand_grid grid(l_e_pair.first, l_e_pair.second);
   grid.fill(mol);
   for (unsigned int ip=0; ip<problematics.size(); ip++) {

      // a trapped residue is now treated as a primary (bonding)
      // residue with a distance to the "bad" solution of 3.5A (this
      // distance does not matter) - for initial positioning.  It is
      // not added to the primaries, so this bond length is not
      // optimised.
      
      std::vector<std::pair<lig_build::pos_t, double> > attachment_points;
      lig_build::pos_t attach_pos = residue_circles[problematics[ip]].pos;
      attach_pos+= lig_build::pos_t (5* double(rand())/double(RAND_MAX),
				     5* double(rand())/double(RAND_MAX));
      
      std::pair<lig_build::pos_t, double> p(attach_pos, 3.5);

      attachment_points.push_back(p);
      initial_primary_residue_circles_layout(grid, problematics[ip],
					     attachment_points);
   }

   
   // set the initial positions for optimisation, now that we have
   // fiddled with residue_circles
   std::vector<residue_circle_t> current_circles = residue_circles;

   for (int iround=0; iround<120; iround++) {
      std::pair<int, std::vector<lbg_info_t::residue_circle_t> > new_c =
	 optimise_residue_circle_positions(residue_circles, current_circles, primary_indices);
      current_circles = new_c.second;
      if (new_c.first == GSL_ENOPROG)
	 break;
      if (new_c.first == GSL_SUCCESS)
	 break;
   }
   // now transfer results from the updating/tmp current_circles to the class data item:
   residue_circles = current_circles;
}

void
lbg_info_t::recentre_considering_residue_centres() {

   std::pair<bool,lig_build::pos_t> tl = get_residue_circles_top_left();

   // std::cout << "----------------------- top left correction " << tl.first << " "
   // << tl.second << std::endl;
  
    if (tl.first) { 
       top_left_correction = lig_build::pos_t(-tl.second.x+40, -tl.second.y+30);
       clear_and_redraw(top_left_correction);
    }
} 


// for debugging the minimization vs original positions
std::vector<lbg_info_t::residue_circle_t>
lbg_info_t::offset_residues_from_orig_positions() {

   std::vector<residue_circle_t> current_circles = residue_circles;
   for (unsigned int i=0; i<current_circles.size(); i++) {
      double delta_x = 5 * double(rand())/double(RAND_MAX);
      double delta_y = 5 * double(rand())/double(RAND_MAX);
      current_circles[i].pos.x += delta_x;
      current_circles[i].pos.y += delta_y;
   }
   return current_circles;
}

void
lbg_info_t::draw_all_flev_annotations() {

   draw_all_flev_residue_attribs();
   draw_all_flev_ligand_annotations();
} 

void
lbg_info_t::draw_all_flev_residue_attribs() {

   if (draw_flev_annotations_flag) {
      draw_residue_circles(residue_circles, additional_representation_handles);
      draw_bonds_to_ligand();
      draw_stacking_interactions(residue_circles);
   }
}

// top left and bottom right corners.
//
std::pair<lig_build::pos_t, lig_build::pos_t>
lbg_info_t::flev_residues_extents() const {

   std::pair<lig_build::pos_t, lig_build::pos_t> p; // defaults with (-1, -1) for coordinates
   if (draw_flev_annotations_flag) {
      lig_build::pos_t ligand_centre = mol.get_ligand_centre();
      GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));
      if (residue_circles.size()) {
	 p.first  = lig_build::pos_t( 9999, 9999);
	 p.second = lig_build::pos_t(-9999, -9999);
	 for (unsigned int i=0; i<residue_circles.size(); i++) {
	    const lig_build::pos_t &pos = residue_circles[i].pos;
	    if (pos.x < p.first.x) p.first.x = pos.x;
	    if (pos.y < p.first.y) p.first.y = pos.y;
	    if (pos.x > p.second.x) p.second.x = pos.x;
	    if (pos.y > p.second.y) p.second.y = pos.y;
	 }
      }
   }
   return p;
}

void
lbg_info_t::draw_all_flev_ligand_annotations() {

   if (draw_flev_annotations_flag) {
      draw_substitution_contour();
      draw_solvent_accessibility_of_atoms();
   } 
} 

std::vector<int>
lbg_info_t::get_primary_indices() const {

   std::vector<int> primary_indices;  // this primary_indices needs to
                                      // get passed to the
                                      // primary_indices used in
                                      // residue cirlce optimization.
   
   for(unsigned int ic=0; ic<residue_circles.size(); ic++) {
      if (residue_circles[ic].is_a_primary_residue()) {
         primary_indices.push_back(ic);
      }
   }
   return primary_indices;
} 


// twiddle residue_circles, taking in to account the residues that
// bond to the ligand (including stacking residues).
// 
void
lbg_info_t::initial_residues_circles_layout() {

   // when we move a primary, we want to know it's index in
   // residue_circles, because that's what we really want to move.
   // 
   // std::vector<std::pair<int, lbg_info_t::residue_circle_t> > primaries;

   // now a class data member because it is used in the layout penalty
   // function (we want to have nice bond lengths for residues bonded
   // to the ligand).
   // 
   std::vector<int> primary_indices;  // this primary_indices needs to
				      // get passed to the
				      // primary_indices used in
				      // residue circle optimization.
   
   for(unsigned int ic=0; ic<residue_circles.size(); ic++) {
      if (residue_circles[ic].is_a_primary_residue()) {
	 primary_indices.push_back(ic);
      }
   }

   // primaries get placed first and are checked for non-crossing
   // ligand interaction bonds (they have penalty scored added if they
   // cross).
   //
   try { 
      std::pair<lig_build::pos_t, lig_build::pos_t> l_e_pair =
	 mol.ligand_extents();
      lbg_info_t::ligand_grid grid(l_e_pair.first, l_e_pair.second);
      grid.fill(mol);

      for (unsigned int iprimary=0; iprimary<primary_indices.size(); iprimary++) {
	 int idx = primary_indices[iprimary];
	 std::vector<std::pair<lig_build::pos_t, double> > attachment_points =
	    residue_circles[idx].get_attachment_points(mol);
	 initial_primary_residue_circles_layout(grid, idx, attachment_points);
      }
      // show_mol_ring_centres();

      position_non_primaries(grid, primary_indices); // untrap residues as needed.
   }
   catch (const std::runtime_error &rte) {
      std::cout << rte.what() << std::endl;
   }
}

// fiddle with the position of the residue_circles[primary_index].
//
void
lbg_info_t::initial_primary_residue_circles_layout(const lbg_info_t::ligand_grid &grid,
						   int primary_index,
						   const std::vector<std::pair<lig_build::pos_t, double> > &attachment_points) {

   if (0) 
      std::cout << "DEBUG:: starting initial_primary_residue_circles_layout() primary_index " << primary_index
		<< " " << residue_circles[primary_index].residue_label << " "
		<< residue_circles[primary_index].residue_type
		<< " has position " << residue_circles[primary_index].pos
		<< std::endl;
   
   if (0)
      std::cout << " =========== adding quadratic for residue " 
		<< residue_circles[primary_index].residue_label
		<< " ============================"
		<< std::endl;
   lbg_info_t::ligand_grid primary_grid = grid;
	 
   // attachment points are points on the ligand, in ligand
   // coordinates to which this primary residue circle is
   // attached (often only one attachment, but can be 2
   // sometimes).
   primary_grid.add_quadratic(attachment_points);

   if (0)
      if (primary_index == 0) 
	 show_grid(grid);
   
   lig_build::pos_t best_pos = primary_grid.find_minimum_position();

   // OK, consider the case where there are 2 residues bonding to the
   // same atom in the ligand.  They wil both be given exactly the
   // same best_pos and they will subsequently refine together,
   // resulting in one residue sitting on top of another - bad news.
   //
   // So let's shift the residue a bit in the direction that it came
   // from so that the residues don't refine together.
   //
   lig_build::pos_t a_to_b_uv = (residue_circles[primary_index].pos - best_pos).unit_vector();
   
   residue_circles[primary_index].pos = best_pos + a_to_b_uv * 4;

   if (0)
      std::cout << "DEBUG::  ending initial_primary_residue_circles_layout() primary_index "
		<< primary_index
		<< " " << residue_circles[primary_index].residue_label << " "
		<< residue_circles[primary_index].residue_type
		<< " has position " << residue_circles[primary_index].pos
		<< std::endl;
   
}

// untrap residues as needed.
void
lbg_info_t::position_non_primaries(const lbg_info_t::ligand_grid &grid,
				   const std::vector<int> &primary_indices) {

   std::vector<lbg_info_t::grid_index_t> already_positioned;
   
   for (unsigned int irc=0; irc<residue_circles.size(); irc++) {
      // if not a primary...
      if (std::find(primary_indices.begin(), primary_indices.end(), irc) ==
	  primary_indices.end()) {
	 
	 // if this point has rejection underneath it, find the
	 // position in grid that doesn't have any rejection.
	 //

	 std::pair<lbg_info_t::grid_index_t, lig_build::pos_t> pos =
	    grid.find_nearest_zero(residue_circles[irc].pos, already_positioned);
	 residue_circles[irc].pos = pos.second;
	 if (pos.first.is_valid_p())
	    already_positioned.push_back(pos.first);
			 
      }
   }
}


// attachment points are points on the ligand, in ligand coordinates
// to which this primary residue circle is attached (often only one
// attachment, but can be 2 sometimes).
//
void
lbg_info_t::ligand_grid::add_quadratic(const std::vector<std::pair<lig_build::pos_t, double> > &attachment_points) {

   if (attachment_points.size()) {
      double scale_by_n_attach = 1.0/double(attachment_points.size());
      
      for (unsigned int iattach=0; iattach<attachment_points.size(); iattach++) {
	 for (int ix=0; ix<x_size(); ix++) {
	    for (int iy=0; iy<y_size(); iy++) {
	       lig_build::pos_t pos = to_canvas_pos(ix, iy);
	       double d2 = (pos-attachment_points[iattach].first).lengthsq();
	       double val = 0.00002 * d2 * scale_by_n_attach;
	       grid_[ix][iy] += val;
	    }
	 }
      }
   }
}

// can throw an exception
lig_build::pos_t
lbg_info_t::ligand_grid::find_minimum_position() const {

   double best_pos_score = 1000000;
   lig_build::pos_t best_pos;
   for (int ix=0; ix<x_size(); ix++) {
      for (int iy=0; iy<y_size(); iy++) {
	 if (grid_[ix][iy] < best_pos_score) {
	    best_pos_score = grid_[ix][iy];
	    best_pos = to_canvas_pos(ix,iy);
	 }
      }
   }
   if (best_pos_score > (1000000-1))
      throw std::runtime_error("failed to get minimum position from ligand grid");
   return best_pos;
}

// actually, not exactly zero but something small.
//
// Don't return a grid-point/position that matches anything in
// already_positioned.
// 
std::pair<lbg_info_t::grid_index_t, lig_build::pos_t>
lbg_info_t::ligand_grid::find_nearest_zero(const lig_build::pos_t &pos,
					   const std::vector<lbg_info_t::grid_index_t> &already_positioned) const {

   lig_build::pos_t p;
   lbg_info_t::grid_index_t rgi; // initially invalid
   double shortest_dist = 43e23;
   double crit = 0.05; // less than this is effectively zero.

   try { 
      lbg_info_t::grid_index_t gi=grid_pos_nearest(pos);
      if (grid_[gi.i()][gi.j()] < crit) {
	 p = pos; // fine, no change
      } else {
	 // search for someplace else
	 for (int ix=0; ix<x_size(); ix++) {
	    for (int iy=0; iy<y_size(); iy++) {
	       if (grid_[ix][iy] < crit) {
		  lig_build::pos_t gp = to_canvas_pos(ix, iy);
		  double d = (gp - pos).lengthsq();
		  if (d < shortest_dist) {
		     lbg_info_t::grid_index_t candidate_grid_index(ix, iy);
		     // This is OK if there is no other previous
		     // solution at the same position (if there is,
		     // keep trying, of course).
		     if (std::find(already_positioned.begin(), already_positioned.end(),
				   candidate_grid_index) == already_positioned.end()) { 
			shortest_dist = d;
			p = gp;
			rgi = candidate_grid_index;
		     }
		  } 
	       }
	    }
	 }
      }
   }
   catch (const std::runtime_error &rte) {
      // the pos was off the grid.  It won't be trapped inside the
      // ligand, so just return what we were given.
      p = pos;
   } 
   return std::pair<lbg_info_t::grid_index_t, lig_build::pos_t> (rgi, p);
} 


lig_build::pos_t
lbg_info_t::ligand_grid::to_canvas_pos(const double &ii, const double &jj) const {

   lig_build::pos_t p(ii*scale_fac, jj*scale_fac);
   p += top_left;
   return p;
}

// can throw a std::runtime_error.
// 
lbg_info_t::grid_index_t
lbg_info_t::ligand_grid::grid_pos_nearest(const lig_build::pos_t &pos) const {

   lig_build::pos_t p = pos - top_left;
   int idx_x = int(p.x/scale_fac+0.5);
   int idx_y = int(p.y/scale_fac+0.5);

   if ((idx_x < 0) || (idx_x >= x_size()) || (idx_y < 0) || (idx_y >= y_size()))
       throw std::runtime_error("out of grid index");
   
   return lbg_info_t::grid_index_t(idx_x, idx_y);
}

 


std::pair<int, int>
lbg_info_t::ligand_grid::canvas_pos_to_grid_pos(const lig_build::pos_t &pos) const {

   lig_build::pos_t p = pos - top_left;
   int ix = (int)(p.x/scale_fac);
   int iy = (int)(p.y/scale_fac);
   return std::pair<int, int> (ix, iy);
}

// for marching squares, ii and jj are the indices of the bottom left-hand side.
int 
lbg_info_t::ligand_grid::square_type(int ii, int jj, float contour_level) const {

   int square_type = lbg_info_t::ligand_grid::MS_NO_SQUARE;
   if ((ii+1) >= x_size_) {
      return lbg_info_t::ligand_grid::MS_NO_SQUARE;
   } else { 
      if ((jj+1) >= y_size_) {
	 return lbg_info_t::ligand_grid::MS_NO_SQUARE;
      } else {
	 float v00 = get(ii, jj);
	 float v01 = get(ii, jj+1);
	 float v10 = get(ii+1, jj);
	 float v11 = get(ii+1, jj+1);
	 if (v00 > contour_level) { 
	    if (v01 > contour_level) { 
	       if (v10 > contour_level) { 
		  if (v11 > contour_level) {
		     return lbg_info_t::ligand_grid::MS_NO_CROSSING;
		  }
	       }
	    }
	 }
	 if (v00 < contour_level) { 
	    if (v01 < contour_level) { 
	       if (v10 < contour_level) { 
		  if (v11 < contour_level) {
		     return lbg_info_t::ligand_grid::MS_NO_CROSSING;
		  }
	       }
	    }
	 }

	 // OK, so it is not either of the trivial cases (no
	 // crossing), there are 14 other variants.
	 // 
	 if (v00 < contour_level) { 
	    if (v01 < contour_level) { 
	       if (v10 < contour_level) { 
		  if (v11 < contour_level) {
		     return lbg_info_t::ligand_grid::MS_NO_CROSSING;
		  } else {
		     return lbg_info_t::ligand_grid::MS_UP_1_1;
		  }
	       } else {
		  if (v11 < contour_level) {
		     return lbg_info_t::ligand_grid::MS_UP_1_0;
		  } else {
		     return lbg_info_t::ligand_grid::MS_UP_1_0_and_1_1;
		  }
	       }
	    } else {

	       // 0,1 is up
	       
	       if (v10 < contour_level) { 
		  if (v11 < contour_level) {
		     return lbg_info_t::ligand_grid::MS_UP_0_1;
		  } else {
		     return lbg_info_t::ligand_grid::MS_UP_0_1_and_1_1;
		  }
	       } else {
		  if (v11 < contour_level) {
		     return lbg_info_t::ligand_grid::MS_UP_0_1_and_1_0;      // hideous valley
		  } else {
		     return lbg_info_t::ligand_grid::MS_UP_0_1_and_1_0_and_1_1; // (only 0,0 down)
		  }
	       }
	    }
	 } else {

	    // 0,0 is up
	    
	    if (v01 < contour_level) { 
	       if (v10 < contour_level) { 
		  if (v11 < contour_level) {
		     return lbg_info_t::ligand_grid::MS_UP_0_0;
		  } else {
		     return lbg_info_t::ligand_grid::MS_UP_0_0_and_1_1; // another hideous valley
		  }
	       } else {
		  // 1,0 is up
		  if (v11 < contour_level) {
		     return lbg_info_t::ligand_grid::MS_UP_0_0_and_1_0;
		  } else {
		     return lbg_info_t::ligand_grid::MS_UP_0_0_and_1_0_and_1_1; // 0,1 is down
		  }
	       }
	    } else {

	       // 0,1 is up
	       
	       if (v10 < contour_level) { 
		  if (v11 < contour_level) {
		     return lbg_info_t::ligand_grid::MS_UP_0_0_and_0_1;
		  } else {
		     return lbg_info_t::ligand_grid::MS_UP_0_0_and_0_1_and_1_1; // 1,0 is down
		  }
	       } else {
		  // if we get here, this test must pass.
		  if (v11 < contour_level) {
		     return lbg_info_t::ligand_grid::MS_UP_0_0_and_0_1_and_1_0; // only 1,1 is down
		  }
	       }
	    }
	 } 
      }
   }
   return square_type;
}

lbg_info_t::contour_fragment::contour_fragment(int ms_type,
					       const float &contour_level, 
					       const lbg_info_t::grid_index_t &grid_index_prev,
					       const lbg_info_t::grid_index_t &grid_index,
					       const lbg_info_t::ligand_grid &grid) {

   int ii_next = lbg_info_t::grid_index_t::INVALID_INDEX;
   int jj_next = lbg_info_t::grid_index_t::INVALID_INDEX;

   float v00 = grid.get(grid_index.i(),   grid_index.j());
   float v01 = grid.get(grid_index.i(),   grid_index.j()+1);
   float v10 = grid.get(grid_index.i()+1, grid_index.j());
   float v11 = grid.get(grid_index.i()+1, grid_index.j()+1);

   float frac_x1 = -1; 
   float frac_y1 = -1;
   float frac_x2 = -1;  // for hideous valley
   float frac_y2 = -1;

   lbg_info_t::contour_fragment::coordinates c1(0.0, X_AXIS_LOW);
   lbg_info_t::contour_fragment::coordinates c2(Y_AXIS_LOW, 0.0);
   cp_t p(c1,c2);
   
   switch (ms_type) {

   case lbg_info_t::ligand_grid::MS_UP_0_0:
   case lbg_info_t::ligand_grid::MS_UP_0_1_and_1_0_and_1_1:

      frac_x1 = (v00-contour_level)/(v00-v10);
      frac_y1 = (v00-contour_level)/(v00-v01);
      c1 = lbg_info_t::contour_fragment::coordinates(frac_x1, X_AXIS_LOW);
      c2 = lbg_info_t::contour_fragment::coordinates(Y_AXIS_LOW, frac_y1);
      p = cp_t(c1,c2);
      coords.push_back(p);
      break;

   case lbg_info_t::ligand_grid::MS_UP_0_1:
   case lbg_info_t::ligand_grid::MS_UP_0_0_and_1_0_and_1_1:
      
      // these look fine
      // std::cout << " ----- case MS_UP_0,1 " << std::endl;
      frac_y1 = (v00 - contour_level)/(v00-v01);
      frac_x2 = (v01 - contour_level)/(v01-v11);
      c1 = lbg_info_t::contour_fragment::coordinates(Y_AXIS_LOW, frac_y1);
      c2 = lbg_info_t::contour_fragment::coordinates(frac_x2, X_AXIS_HIGH);
      p = cp_t(c1,c2);
      coords.push_back(p);
      break;

      
   case lbg_info_t::ligand_grid::MS_UP_1_0:
   case lbg_info_t::ligand_grid::MS_UP_0_0_and_0_1_and_1_1:
      
      // std::cout << " ----- case MS_UP_1,0 " << std::endl;
      frac_x1 = (contour_level - v00)/(v10-v00);
      frac_y1 = (contour_level - v10)/(v11-v10);
      c1 = lbg_info_t::contour_fragment::coordinates(frac_x1, X_AXIS_LOW);
      c2 = lbg_info_t::contour_fragment::coordinates(Y_AXIS_HIGH, frac_y1);
      p = cp_t(c1,c2);
      coords.push_back(p);
      break;

      
      
   case lbg_info_t::ligand_grid::MS_UP_1_1:
   case lbg_info_t::ligand_grid::MS_UP_0_0_and_0_1_and_1_0:

      // std::cout << " ----- case MS_UP_1,1 " << std::endl;
      frac_x1 = (v01-contour_level)/(v01-v11);
      frac_y1 = (v10-contour_level)/(v10-v11);
      c1 = lbg_info_t::contour_fragment::coordinates(frac_x1, X_AXIS_HIGH);
      c2 = lbg_info_t::contour_fragment::coordinates(Y_AXIS_HIGH, frac_y1);
      p = cp_t(c1,c2);
      coords.push_back(p);
      break;

   case lbg_info_t::ligand_grid::MS_UP_0_0_and_0_1:
   case lbg_info_t::ligand_grid::MS_UP_1_0_and_1_1:
      
      // std::cout << " ----- case MS_UP_0,0 and 0,1 " << std::endl;
      frac_x1 = (v00-contour_level)/(v00-v10);
      frac_x2 = (v01-contour_level)/(v01-v11);
      c1 = lbg_info_t::contour_fragment::coordinates(frac_x1, X_AXIS_LOW);
      c2 = lbg_info_t::contour_fragment::coordinates(frac_x2, X_AXIS_HIGH);
      p = cp_t(c1,c2);
      coords.push_back(p);
      break;

   case lbg_info_t::ligand_grid::MS_UP_0_0_and_1_0:
   case lbg_info_t::ligand_grid::MS_UP_0_1_and_1_1:
      
      // std::cout << " ----- case MS_UP_0,0 and 1,0 " << std::endl;
      frac_y1 = (v00-contour_level)/(v00-v01);
      frac_y2 = (v10-contour_level)/(v10-v11);
      c1 = lbg_info_t::contour_fragment::coordinates(Y_AXIS_LOW,  frac_y1);
      c2 = lbg_info_t::contour_fragment::coordinates(Y_AXIS_HIGH, frac_y2);
      p = cp_t(c1,c2);
      coords.push_back(p);
      break;
      

   default:
      std::cout << "ERROR:: unhandled square type: " << ms_type << std::endl;

   } 

}

std::vector<std::pair<lig_build::pos_t, double> >
lbg_info_t::residue_circle_t::get_attachment_points(const widgeted_molecule_t &mol) const {

   // debug
//    std::cout << "........ in get_attachment_points() mol has " << mol.atoms.size()
// 	     << " atoms " << std::endl;
//    for (unsigned int i=0; i<mol.atoms.size(); i++) { 
//       std::cout << "    " << i << "   " << mol.atoms[i] << " "
// 		<< mol.atoms[i].get_atom_name() <<  std::endl;
//    }
   
   
   std::vector<std::pair<lig_build::pos_t, double> > v;

   for (unsigned int i=0; i<bonds_to_ligand.size(); i++) {
      if (bonds_to_ligand[i].is_set()) {
	 try {
	    lig_build::pos_t pos =
	       mol.get_atom_canvas_position(bonds_to_ligand[i].ligand_atom_name);
	    if (bonds_to_ligand[i].is_set()) { 
	       std::pair<lig_build::pos_t, double> p(pos, bonds_to_ligand[i].bond_length);
	       v.push_back(p);
	    }
	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << "WARNING:: " << rte.what() << std::endl;
	 }
      }
   }

   // with a ring system on the ligand, that is.
   // 
   if (has_ring_stacking_interaction()) {
      try {
	 lig_build::pos_t pos = mol.get_ring_centre(ligand_ring_atom_names);
	 double stacking_dist = 4.5; // A
	 std::pair<lig_build::pos_t, double> p(pos, stacking_dist);
	 v.push_back(p);
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "WARNING:: " << rte.what() << std::endl;
      }
   }

   if (get_stacking_type() == lbg_info_t::residue_circle_t::CATION_PI_STACKING) {
      try {
	 std::string at_name = get_ligand_cation_atom_name();
	 double stacking_dist = 4.2; // a bit shorter, because we
				     // don't have to go from the
				     // middle of a ring system, the
				     // target point (an atom) is more
				     // accessible.  Still 3.5 was too
				     // short (->crowded) when
				     // showing 3 cation-pi interactions
				     // on a alkylated N.
	 
	 lig_build::pos_t pos = mol.get_atom_canvas_position(at_name);
	 std::pair<lig_build::pos_t, double> p(pos, stacking_dist);
	 v.push_back(p);
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "WARNING:: " << rte.what() << std::endl;
      }
   }

   return v;
}


// return 1 on solution having problems, 0 for no problems, also
// return a list of the residue circles with problems.
// 
std::pair<bool, std::vector<int> >
lbg_info_t::solution_has_problems_p() const {

   std::vector<int> v;
   bool status = 0;

   // scale factor of 1.2 is too small, 2.0 seems good for H102 in 2aei 
   double crit_dist_2 = 2.0 * SINGLE_BOND_CANVAS_LENGTH * SINGLE_BOND_CANVAS_LENGTH;
   double crit_dist_wide_2 = 4 * crit_dist_2;

   // a problem can be 2 atoms < crit_dist_2 or
   //                  5 atoms < crit_dist_wide_2
   // 
   // we need the wider check to find atoms that are enclosed by a
   // system of 10 or so atoms.
   // 

   for (unsigned int i=0; i<residue_circles.size(); i++) {
      int n_close = 0;
      int n_close_wide = 0;
      for (unsigned int j=0; j<mol.atoms.size(); j++) {
	 double d2 = (residue_circles[i].pos-mol.atoms[j].atom_position).lengthsq();
//  	 std::cout << "comparing " << d2 << " " << crit_dist_2 << " "
//  		   << residue_circles[i].residue_label
//  		   << std::endl;
	 if (d2 < crit_dist_2) {
	    n_close++;
	    if (n_close > 1) {
	       v.push_back(i);
	       status = 1;
	       break;
	    }
	 }

	 if (d2 < crit_dist_wide_2) {
	    n_close_wide++;
	    if (n_close_wide > 5) {
	       v.push_back(i);
	       status = 1;
	       break;
	    } 
	 } 
      }
   }
   return std::pair<bool, std::vector<int> > (status, v);
} 
 



void
lbg_info_t::show_grid(const lbg_info_t::ligand_grid &grid) {

   int n_objs = 0;
   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));
   for (int ix=0; ix<grid.x_size(); ix+=1) {
      for (int iy=0; iy<grid.y_size(); iy+=1) {
	 lig_build::pos_t pos = grid.to_canvas_pos(ix, iy);
	 double intensity = grid.get(ix, iy);
	 int int_as_int = int(255 * intensity);
	 std::string colour = grid_intensity_to_colour(int_as_int);
 	 GooCanvasItem *rect =
 	    goo_canvas_rect_new(root, pos.x, pos.y, 3.5, 3.5,
 				"line-width", 1.0, // in show_grid()
 				"fill_color", colour.c_str(),
				"stroke-color", colour.c_str(),
 				NULL);
	 n_objs++;
      }
   }

   if (1) {  // debugging
      grid.show_contour(root, 0.1);
      grid.show_contour(root, 0.2);
      grid.show_contour(root, 0.3);
      grid.show_contour(root, 0.4);
      grid.show_contour(root, 0.5);
      grid.show_contour(root, 0.6);
      grid.show_contour(root, 0.7);
      grid.show_contour(root, 0.8);
      grid.show_contour(root, 0.9);
      grid.show_contour(root, 1.0);
   }
   
}

void
lbg_info_t::ligand_grid::avoid_ring_centres(std::vector<std::vector<std::string> > ring_atoms_list,
					    const widgeted_molecule_t &mol) {

   // For the substitution contour we blob in a circle of radius
   // 1/(2*sin(180/n_ring_atoms)) bond lengths.
   
   for (unsigned int iring=0; iring<ring_atoms_list.size(); iring++) {
      try { 
	 lig_build::pos_t centre = mol.get_ring_centre(ring_atoms_list[iring]);
	 int n_atoms = ring_atoms_list[iring].size();
	 // just in case we get here with n_atoms = 1 for some reason...
	 if (n_atoms < 3) n_atoms = 3;
	 
	 double radius = 1/(2*sin(M_PI/double(n_atoms))) * 1.5; // in "A" or close
	 // std::cout << "avoid_ring_centres() adding ring centre at " << centre
	 // << " n_atoms: " << n_atoms << " radius " << radius << std::endl;
	 add_for_accessibility(radius, centre);
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "Opps - failed to find ring centre for ring atom name "
		   << iring << std::endl;
      }
   }
}

void
lbg_info_t::ligand_grid::show_contour(GooCanvasItem *root, float contour_level) const {

   std::vector<widgeted_atom_ring_centre_info_t> dummy_unlimited_atoms;
   std::vector<std::vector<std::string> > dummy_ring_atom_names;   
   show_contour(root, contour_level, dummy_unlimited_atoms, dummy_ring_atom_names);
}

void
lbg_info_t::ligand_grid::show_contour(GooCanvasItem *root, float contour_level,
				      const std::vector<widgeted_atom_ring_centre_info_t> &unlimited_atoms,
				      const std::vector<std::vector<std::string> > &ring_atoms_list) const {


   GooCanvasItem *group = goo_canvas_group_new (root, "stroke-color", "#880000",
						NULL);
   bool debug = 0;
   int ii=0;
   int jj=0;

   // fill the ring centre vector, if the unlimited atom have ring centres
   std::vector<std::pair<bool, lig_build::pos_t> > ring_centres(unlimited_atoms.size());
   for (unsigned int i=0; i<unlimited_atoms.size(); i++) { 
      try {
	 lig_build::pos_t p;
      }
      catch (const std::runtime_error &rte) {
	 std::cout << rte.what() << std::endl;
      } 
   }

   std::vector<std::pair<lig_build::pos_t, lig_build::pos_t> > line_fragments;

   lbg_info_t::grid_index_t grid_index_prev(0,0);

   for (int ix=0; ix<x_size(); ix+=1) {
      for (int iy=0; iy<y_size(); iy+=1) {
	 int ms_type = square_type(ix, iy, contour_level);

	 lbg_info_t::grid_index_t grid_index(ix,iy);

	 if ((ms_type != MS_NO_CROSSING) && (ms_type != MS_NO_SQUARE)) { 
	    lbg_info_t::contour_fragment cf(ms_type, contour_level,
					    grid_index_prev,
					    grid_index,
					    *this); // sign of bad architecture
	    if (cf.coords.size() == 1) {

	       if (debug)
		  std::cout << "plot contour ("
			    << cf.get_coords(ix, iy, 0).first << " "
			    << cf.get_coords(ix, iy, 0).second << ") to ("
			    << cf.get_coords(ix, iy, 1).first << " "
			    << cf.get_coords(ix, iy, 1).second << ")" << std::endl;
			
	       std::pair<double, double> xy_1 = cf.get_coords(ix, iy, 0);
	       std::pair<double, double> xy_2 = cf.get_coords(ix, iy, 1);
	       lig_build::pos_t pos_1 = to_canvas_pos(xy_1.first, xy_1.second);
	       lig_build::pos_t pos_2 = to_canvas_pos(xy_2.first, xy_2.second);

	       lig_build::pos_t p1 = to_canvas_pos(cf.get_coords(ix, iy, 0).first,
						   cf.get_coords(ix, iy, 0).second);
	       lig_build::pos_t p2 = to_canvas_pos(cf.get_coords(ix, iy, 1).first,
						   cf.get_coords(ix, iy, 1).second);
	       std::pair<lig_build::pos_t, lig_build::pos_t> fragment_pair(p1, p2);

	       // Now filter out this fragment pair if it is too close
	       // to an unlimited_atom_positions
	       bool plot_it = 1;
	       double dist_crit = 4.0 * LIGAND_TO_CANVAS_SCALE_FACTOR;

	       
	       for (unsigned int i=0; i<unlimited_atoms.size(); i++) { 
// 		  lig_build::pos_t p = to_canvas_pos(unlimited_atom_positions[i].x,
// 						     unlimited_atom_positions[i].y);
		  
		  lig_build::pos_t p = unlimited_atoms[i].atom.atom_position;

		  // if this atom has a ring centre, use the ring
		  // centre to atom vector to unplot vectors only in a
		  // particular direction.
		  // 
		  if ((p - p1).lengthsq() < (dist_crit * dist_crit)) {
		     if (1) { // for debugging
			if (unlimited_atoms[i].has_ring_centre_flag) {
			   // std::cout << " atom " << i << " has ring_centre ";
			   lig_build::pos_t d_1 =
			      unlimited_atoms[i].ring_centre - unlimited_atoms[i].atom.atom_position;
			   lig_build::pos_t d_2 = unlimited_atoms[i].atom.atom_position - p1;
			   double cos_theta =
			      lig_build::pos_t::dot(d_1, d_2)/(d_1.length()*d_2.length());
			   // std::cout << " cos_theta " << cos_theta << " for unlimited atom " << i << std::endl;
			   if (cos_theta > 0.3) { // only cut in the "forwards" direction

// 			      std::cout << " cutting by ring-centred unlimited atom " << i << " "
// 				 << unlimited_atoms[i].atom.get_atom_name()
// 					<< std::endl;
			      
			      plot_it = 0;
			      break;
			   }
			   // std::cout << std::endl;
			
			} else {
			   plot_it = 0;
//  			   std::cout << " cutting by unlimited atom " << i << " "
//  				     << unlimited_atoms[i].atom.get_atom_name()
//  				     << std::endl;
			   break;
			}
		     }
		  }
		  
	       } // end unlimited atoms loop

	       if (plot_it)
		  line_fragments.push_back(fragment_pair);
	    } 
	 }
      }
   }

   std::vector<std::vector<lig_build::pos_t> > contour_lines = make_contour_lines(line_fragments);

   plot_contour_lines(contour_lines, root);

   // check the orientation of the canvas
   if (0) { 
      lig_build::pos_t grid_ori = to_canvas_pos(0.0, 0.0);
      goo_canvas_rect_new (group,
			   grid_ori.x, grid_ori.y, 5.0, 5.0,
			   "line-width", 1.0, // in show_contour()
			   "stroke-color", "green",
			   "fill_color", "blue",
			   NULL);
   }
} 

std::vector<std::vector<lig_build::pos_t> >
lbg_info_t::ligand_grid::make_contour_lines(const std::vector<std::pair<lig_build::pos_t, lig_build::pos_t> > &line_fragments) const {

   // Look for neighboring cells so that a continuous line fragment is
   // generated. That is fine for closed contours, but this handled
   // badly contour line fragments that go over the edge, they will be
   // split into several line fragment - depending on which fragment
   // is picked first.
   //
   // Perhaps the last items in the line_frag_queue should be (any)
   // points that are on the edge of the grid - because they are the
   // best starting place for dealing with contour line fragments that
   // go over the edge.  If you do this, then the input arg
   // line_fragments should carry with it info about this line
   // fragment being on an edge (generated in show_contour())..
   //
   
   std::vector<std::vector<lig_build::pos_t> > v; // returned item.

   std::deque<std::pair<lig_build::pos_t, lig_build::pos_t> > line_frag_queue;
   for (unsigned int i=0; i<line_fragments.size(); i++)
      line_frag_queue.push_back(line_fragments[i]);

   while (!line_frag_queue.empty()) {

      // start a new line fragment
      // 
      std::vector<lig_build::pos_t> working_points;
      std::pair<lig_build::pos_t, lig_build::pos_t> start_frag = line_frag_queue.back();
      working_points.push_back(start_frag.first);
      working_points.push_back(start_frag.second);
      line_frag_queue.pop_back();
      bool line_fragment_terminated = 0;

      while (! line_fragment_terminated) { 
	 // Is there a fragment in line_frag_queue that has a start that
	 // is working_points.back()?
	 //
	 bool found = 0;
	 std::deque<std::pair<lig_build::pos_t, lig_build::pos_t> >::iterator it;
	 for (it=line_frag_queue.begin(); it!=line_frag_queue.end(); it++) { 
	    if (it->first.near_point(working_points.back(), 0.1)) {
	       working_points.push_back(it->second);
	       line_frag_queue.erase(it);
	       found = 1;
	       break;
	    }
	 }
	 if (! found)
	    line_fragment_terminated = 1;
      }

      v.push_back(working_points);
   }
   
   return v;
} 

void
lbg_info_t::ligand_grid::plot_contour_lines(const std::vector<std::vector<lig_build::pos_t> > &contour_lines, GooCanvasItem *root) const {

   GooCanvasItem *group = goo_canvas_group_new (root, "stroke-color",
						// "#111111", not so dark
						"#777777",
						NULL);
   GooCanvasLineDash *dash = goo_canvas_line_dash_new (2, 1.5, 2.5);

   for (unsigned int i=0; i<contour_lines.size(); i++) { 
      for (int j=0; j<int(contour_lines[i].size()-1); j++) {
	 
	 goo_canvas_polyline_new_line(group,
				      contour_lines[i][j].x,
				      contour_lines[i][j].y,
				      contour_lines[i][j+1].x,
				      contour_lines[i][j+1].y,
				      "line_width", 1.0,
				      "line-dash", dash,
				      NULL);
	 
      }
   }
   goo_canvas_line_dash_unref(dash);
   
} 




void
lbg_info_t::show_mol_ring_centres() {

   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));
   std::vector<lig_build::pos_t> c = mol.get_ring_centres();
   for (unsigned int i=0; i<c.size(); i++) {
      GooCanvasItem *rect_item = goo_canvas_rect_new (root,
						      c[i].x, c[i].y, 4.0, 4.0,
						      "line-width", 1.0, // in show_mol_ring_centres()
						      "stroke-color", "blue",
						      "fill_color", "blue",
						      NULL);
   }
}

void
lbg_info_t::show_unlimited_atoms(const std::vector<widgeted_atom_ring_centre_info_t> &ua) {
   
   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));
   for (unsigned int i=0; i<ua.size(); i++) {
      GooCanvasItem *rect_item =
	 goo_canvas_rect_new(root,
			     ua[i].atom.atom_position.x -6.0,
			     ua[i].atom.atom_position.y -6.0,
			     12.0, 12.0,
			     "line-width", 1.0, // in show_unlimited_atoms()
			     "stroke-color", "lightblue",
			     "fill_color", "lightblue",
			     NULL);
   }
}

void
lbg_info_t::show_ring_centres(std::vector<std::vector<std::string> > ring_atoms_list,
			      const widgeted_molecule_t &mol) {

   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));
   for (unsigned int iring=0; iring<ring_atoms_list.size(); iring++) {
      try { 
	 lig_build::pos_t ring_centre = mol.get_ring_centre(ring_atoms_list[iring]);
	 GooCanvasItem *rect_item =
	    goo_canvas_rect_new(root,
				ring_centre.x -6.0,
				ring_centre.y -6.0,
				12.0, 12.0,
				"line-width", 1.0, // in show_ring_centres()
				"stroke-color", "purple",
				"fill_color", "purple",
				NULL);

      }
      catch (const std::runtime_error &rte) {
	 std::cout << "Opps - failed to find ring centre for ring atom name "
		   << iring << std::endl;
      }
   }
}


// this can cache ring centres in mol if they are not there already.
void
lbg_info_t::show_ring_centres(widgeted_molecule_t &mol) {

   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));
   std::vector<lig_build::pos_t> ring_centre_list = mol.get_ring_centres();

   std::cout << "found " << ring_centre_list.size() << " ring centres" << std::endl;
   for (unsigned int iring=0; iring<ring_centre_list.size(); iring++) {
      std::cout << "   " << ring_centre_list[iring] << std::endl;
      GooCanvasItem *item = goo_canvas_ellipse_new(root,
						   ring_centre_list[iring].x, ring_centre_list[iring].y,
						   12.0, 12.0,
						   "line-width", 1.0,
						   "stroke-color", "purple",
						   "fill_color", "purple",
						   NULL);
   }

}



// convert val [0->255] to a hex colour string in red
// 
std::string
lbg_info_t::grid_intensity_to_colour(int val) const {

   std::string colour = "#000000";
   int step = 8;
      
   // higest value -> ff0000
   // lowest value -> ffffff

   int inv_val = 255 - val;
   if (inv_val < 0)
      inv_val = 0;
   
   int val_high = (inv_val & 240)/16;
   int val_low  = inv_val &  15;

   std::string s1 = sixteen_to_hex_let(val_high);
   std::string s2 = sixteen_to_hex_let(val_low);

   colour = "#ff";
   colour += s1;
   colour += s2;
   colour += s1;
   colour += s2;

   // std::cout << "val " << val << " returns colour " << colour << "\n";
   return colour;
}

std::string
lbg_info_t::sixteen_to_hex_let(int v) const {

   std::string l = "0";
   if (v == 15)
      l = "f";
   if (v == 14)
      l = "e";
   if (v == 13)
      l = "d";
   if (v == 12)
      l = "c";
   if (v == 11)
      l = "b";
   if (v == 10)
      l = "a";
   if (v == 9)
      l = "9";
   if (v == 8)
      l = "8";
   if (v == 7)
      l = "7";
   if (v == 6)
      l = "6";
   if (v == 5)
      l = "5";
   if (v == 4)
      l = "4";
   if (v == 3)
      l = "3";
   if (v == 2)
      l = "2";
   if (v == 1)
      l = "1";

   return l;

}

lbg_info_t::ligand_grid::ligand_grid(const lig_build::pos_t &low_x_and_y,
				     const lig_build::pos_t &high_x_and_y) {

   double extra_extents = 100;
   top_left     = low_x_and_y  - lig_build::pos_t(extra_extents, extra_extents);
   bottom_right = high_x_and_y + lig_build::pos_t(extra_extents, extra_extents);
   scale_fac = 5; // seems good
   double delta_x = bottom_right.x - top_left.x;
   double delta_y = bottom_right.y - top_left.y;
   if (0)
      std::cout << " in making grid, got delta_x and delta_y "
		<< delta_x << " " << delta_y << std::endl;
   x_size_ = int(delta_x/scale_fac+1);
   y_size_ = int(delta_y/scale_fac+1);

   std::vector<double> tmp_y (y_size_, 0.0);
   grid_.resize(x_size_);
   for (int i=0; i<x_size_; i++)
      grid_[i] = tmp_y;
   if (0) 
      std::cout << "--- ligand_grid constructor, grid has extents "
		<< x_size_ << " " << y_size_ << " real " << grid_.size()
		<< " " << grid_[0].size()
		<< std::endl;
   
}

// arg is not const reference because get_ring_centres() caches the return value.
void
lbg_info_t::ligand_grid::fill(widgeted_molecule_t mol) {

   double exp_scale = 0.0011;
   double rk = 3000.0;

   // int grid_extent = 15; // 10, 12 is not enough
   int grid_extent = 50 ; // untraps 2wot residues?
   
   for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
      for (int ipos_x= -grid_extent; ipos_x<=grid_extent; ipos_x++) {
	 for (int ipos_y= -grid_extent; ipos_y<=grid_extent; ipos_y++) {
	    std::pair<int, int> p = canvas_pos_to_grid_pos(mol.atoms[iat].atom_position);
	    int ix_grid = ipos_x + p.first;
	    int iy_grid = ipos_y + p.second;
	    if ((ix_grid >= 0) && (ix_grid < x_size())) {
	       if ((iy_grid >= 0) && (iy_grid < y_size())) {
		  double d2 = (to_canvas_pos(ix_grid, iy_grid) - mol.atoms[iat].atom_position).lengthsq();
		  double val =  rk * exp(-0.5*exp_scale*d2);
		  grid_[ix_grid][iy_grid] += val;
	       } else {
// 		  std::cout << "ERROR:: out of range in y: " << ix_grid << "," << iy_grid << " "
// 			    << "and grid size: " << x_size() << "," << y_size() << std::endl;
	       }
	    } else {
// 	       std::cout << "ERROR:: out of range in x: " << ix_grid << "," << iy_grid << " "
// 			 << "and grid size: " << x_size() << "," << y_size() << std::endl;
	    }
	 }
      }
   }

   std::vector<lig_build::pos_t> mol_ring_centres = mol.get_ring_centres();

   // std::cout << "DEBUG:: found " << mol_ring_centres.size() << " ring centres " << std::endl;
   
   for (unsigned int ir=0; ir<mol_ring_centres.size(); ir++) { 
      for (int ipos_x= -10; ipos_x<=10; ipos_x++) {
	 for (int ipos_y= -10; ipos_y<=10; ipos_y++) {
	    std::pair<int, int> p = canvas_pos_to_grid_pos(mol_ring_centres[ir]);
	    int ix_grid = ipos_x + p.first;
	    int iy_grid = ipos_y + p.second;
	    if ((ix_grid >= 0) && (ix_grid < x_size())) {
	       if ((iy_grid >= 0) && (iy_grid < y_size())) {
		  double d2 = (to_canvas_pos(ix_grid, iy_grid) - mol_ring_centres[ir]).lengthsq();
		  double val = rk * exp(-0.5* exp_scale * d2);
		  grid_[ix_grid][iy_grid] += val;
	       }
	    }
	 }
      }
   }
   normalize(); // scaled peak value to 1.
}


void
lbg_info_t::ligand_grid::add_for_accessibility(double bash_dist, const lig_build::pos_t &atom_pos) {

   bool debug = 0;
   int grid_extent = 45;

   double inv_scale_factor = 1.0/double(LIGAND_TO_CANVAS_SCALE_FACTOR);
   
   for (int ipos_x= -grid_extent; ipos_x<=grid_extent; ipos_x++) {
      for (int ipos_y= -grid_extent; ipos_y<=grid_extent; ipos_y++) {
	 std::pair<int, int> p = canvas_pos_to_grid_pos(atom_pos);
	 int ix_grid = ipos_x + p.first;
	 int iy_grid = ipos_y + p.second;
	 if ((ix_grid >= 0) && (ix_grid < x_size())) {
	    if ((iy_grid >= 0) && (iy_grid < y_size())) {
	       
	       double d2 = (to_canvas_pos(ix_grid, iy_grid) - atom_pos).lengthsq();
	       d2 *= (inv_scale_factor * inv_scale_factor);
	       double val = substitution_value(d2, bash_dist);
	       if (debug)
		  if (val > 0)
		     std::cout << "adding " << val << " to grid " << ix_grid << " " << iy_grid
			       << " from " << sqrt(d2) << " vs " << bash_dist << std::endl;
	       grid_[ix_grid][iy_grid] += val;
	    }
	 }
      }
   }
}

void
lbg_info_t::ligand_grid::add_for_accessibility_no_bash_dist_atom(double scale,
								 const lig_build::pos_t &atom_pos) { 

   bool debug = 1;
   int grid_extent = 40;

   double inv_scale_factor = 1.0/double(LIGAND_TO_CANVAS_SCALE_FACTOR);
   
   for (int ipos_x= -grid_extent; ipos_x<=grid_extent; ipos_x++) {
      for (int ipos_y= -grid_extent; ipos_y<=grid_extent; ipos_y++) {
	 std::pair<int, int> p = canvas_pos_to_grid_pos(atom_pos);
	 int ix_grid = ipos_x + p.first;
	 int iy_grid = ipos_y + p.second;
	 if ((ix_grid >= 0) && (ix_grid < x_size())) {
	    if ((iy_grid >= 0) && (iy_grid < y_size())) {
	       
	       double d2 = (to_canvas_pos(ix_grid, iy_grid) - atom_pos).lengthsq();
	       d2 *= (inv_scale_factor * inv_scale_factor);
	       // a triangle function, 1 at the atom centre, 0 at 1.5A and beyond
	       //
	       double dist_at_zero = 2.6;
	       
	       double val = 0.0;
	       if (d2< dist_at_zero * dist_at_zero)
		  val = - 1.0/dist_at_zero * sqrt(d2) + 1.0;
	       grid_[ix_grid][iy_grid] += val;
	    }
	 }
      }
   }
}



// scale peak value to 1.0
void
lbg_info_t::ligand_grid::normalize() {

   double max_int = 0.0;

   // std::cout << "normalizing grid " << x_size() << " by " << y_size() << std::endl;
   for (int ix=0; ix<x_size(); ix++) {
      for (int iy=0; iy<y_size(); iy++) {
	 double intensity = grid_[ix][iy];
	 if (intensity > max_int)
	    max_int = intensity;
      }
   }
   if (max_int > 0.0) {
      double sc_fac = 1.0/max_int;
      for (int ix=0; ix<x_size(); ix++) {
	 for (int iy=0; iy<y_size(); iy++) {
	    grid_[ix][iy] *= sc_fac;
	 }
      }
   }
}



// minimise layout energy
std::pair<int, std::vector<lbg_info_t::residue_circle_t> >
lbg_info_t::optimise_residue_circle_positions(const std::vector<lbg_info_t::residue_circle_t> &r,
					      const std::vector<lbg_info_t::residue_circle_t> &c,
					      const std::vector<int> &primary_indices) const { 

   if (r.size() > 0) {
      if (c.size() == r.size()) { 
	 optimise_residue_circles orc(r, c, mol, primary_indices);
	 int status = orc.get_gsl_min_status();
	 return orc.solution();
      }
   }
   std::vector<lbg_info_t::residue_circle_t> dv; // dummy 
   return std::pair<int, std::vector<lbg_info_t::residue_circle_t> > (0, dv);
}


std::vector<lbg_info_t::residue_circle_t>
lbg_info_t::read_residues(const std::string &file_name) const {
      
   std::vector<residue_circle_t> v;
   std::ifstream f(file_name.c_str());
   if (!f) {
      std::cout << "Failed to open " << file_name << std::endl;
   } else {

      std::cout << "opened residue_circle file: " << file_name << std::endl;
      
      std::vector<std::string> lines;
      std::string line;
      while (std::getline(f, line)) { 
	 lines.push_back(line);
      }

      for (unsigned int i=0; i<lines.size(); i++) {
	 std::vector<std::string> words = coot::util::split_string_no_blanks(lines[i], " ");
	 
	 // debug input
	 if (0) { 
	    std::cout << i << " " << lines[i] << "\n";
	    for (unsigned int j=0; j<words.size(); j++) { 
	       std::cout << "  " << j << " " << words[j] << " ";
	    }
	    std::cout << "\n";
	 }

	 if (words.size() > 5) {
	    if (words[0] == "RES") {
	       try { 
		  double pos_x = lig_build::string_to_float(words[1]);
		  double pos_y = lig_build::string_to_float(words[2]);
		  double pos_z = lig_build::string_to_float(words[3]);
		  clipper::Coord_orth cp(pos_x, pos_y, pos_z);
		  std::string res_type = words[4];
		  std::string label = words[5];
		  coot::residue_spec_t dum_spec;
		  clipper::Coord_orth dum_click_pos(0,0,0);
		  residue_circle_t rc(cp, dum_click_pos, dum_spec, res_type, label);
		  lig_build::pos_t pos = mol.input_coords_to_canvas_coords(cp);
		  if (res_type == "HOH") {
		     if (words.size() > 7) {
			try {
			   double dist_to_protein = lig_build::string_to_float(words[7]);
			   rc.set_water_dist_to_protein(dist_to_protein);
			}
			catch (const std::runtime_error &rte) {
			}
		     } 
		  } 
		  rc.set_canvas_pos(pos);
		  v.push_back(residue_circle_t(rc)); // why is there a constructor here?
	       }
	       catch (const std::runtime_error &rte) {
		  std::cout << "failed to parse :" << lines[i] << ":" << std::endl;
	       }
	    }
	 }
	 if (words.size() == 6) {
	    if (words[0] == "BOND") { // written with space
	       try {
		  double bond_l = lig_build::string_to_float(words[3]);
		  std::string atom_name = lines[i].substr(5,4);
		  int bond_type = lig_build::string_to_int(words[5]);
		  lbg_info_t::bond_to_ligand_t btl(atom_name, bond_l);
		  btl.bond_type = bond_type;
// 		  std::cout << "adding bond " << v.size() << " to :"
// 			    << atom_name << ": " << bond_l << std::endl;
		  if (v.size())
		     v.back().add_bond_to_ligand(btl);
	       }
	       catch (const std::runtime_error &rte) {
		  std::cout << "failed to parse :" << lines[i] << ":" << std::endl;
	       }
	    }
	 }

	 if (words.size() == 3) {
	    if (words[0] == "SOLVENT_ACCESSIBILITY_DIFFERENCE") {
	       try {
		  double se_holo = lig_build::string_to_float(words[1]);
		  double se_apo  = lig_build::string_to_float(words[2]);
		  if (v.size()) { 
		     v.back().set_solvent_exposure_diff(se_holo, se_apo);
		  }
	       }
	       catch (const std::runtime_error &rte) {
		  std::cout << "failed to parse :" << lines[i] << ":" << std::endl;
	       }
	    }
	 }

	 if (words.size() > 2) {
	    if (words[0] == "STACKING") {
	       std::vector<std::string> ligand_ring_atom_names;
	       std::string type = words[1];
	       std::string::size_type l = lines[i].length();
	       int n_atoms_max = 20;
	       for (int iat=0; iat<n_atoms_max; iat++) {
		  std::string::size_type name_index = 19+iat*6;
		  if (l > (name_index+4)) {
		     std::string ss = lines[i].substr(name_index, 4);
		     if (0)
			std::cout << "   Found ligand ring stacking atom :"
				  << ss << ":" << std::endl;
		     ligand_ring_atom_names.push_back(ss);
		  }
	       }
	       if (type == "pi-pi")
		  if (v.size())
		     v.back().set_stacking(type, ligand_ring_atom_names, "");
	       if (type == "pi-cation")
		  if (v.size())
		     v.back().set_stacking(type, ligand_ring_atom_names, "");
	       if (type == "cation-pi") {
		  std::string ligand_cation_atom_name = lines[i].substr(19,4);
		  // std::cout << "DEBUG:: on reading residue info file found ligand cation "
		  // << "name :" << ligand_cation_atom_name << ":" << std::endl;
		  v.back().set_stacking(type, ligand_ring_atom_names, ligand_cation_atom_name);
	       }
	    }
	 }
      }
   }
   std::cout << "found " << v.size() << " residue centres" << std::endl;
   return v;
} 

// "must take exactly one argument" problem
// 
// std::ostream &
// lbg_info_t::operator<<(std::ostream &s, residue_circle_t rc) {

//    s << "res-circ{" << rc.pos << " " << rc.label << " with bond_to_ligand length "
//      << bond_to_ligand.bond_length << "}";

// }


// if you don't have add_rep_handles, then pass a vector or size 0.
// 
void
lbg_info_t::draw_residue_circles(const std::vector<residue_circle_t> &l_residue_circles,
				 const std::vector<int> &add_rep_handles) {

   double max_dist_water_to_ligand_atom  = 3.3; // don't draw waters that are far from ligand
   double max_dist_water_to_protein_atom = 3.3; // don't draw waters that are not somehow 
                                                // attached to the protein.

   if (draw_flev_annotations_flag) { 
      bool draw_solvent_exposures = 1;
      try { 
	 lig_build::pos_t ligand_centre = mol.get_ligand_centre();
	 GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));

	 if (draw_solvent_exposures)
	    for (unsigned int i=0; i<l_residue_circles.size(); i++)
	       draw_solvent_exposure_circle(l_residue_circles[i], ligand_centre, root);

	 for (unsigned int i=0; i<l_residue_circles.size(); i++) {
	    lig_build::pos_t pos = l_residue_circles[i].pos;
	    int add_rep_handle = -1; // default, no handle
	    if (add_rep_handles.size() == l_residue_circles.size())
	       add_rep_handle = add_rep_handles[i];

	    draw_residue_circle_top_layer(l_residue_circles[i], ligand_centre, add_rep_handle);
	 }
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "WARNING:: draw_residue_circles: " << rte.what() << std::endl;
      }
   }
}

// static
/* This handles button presses on the canvas */
// extern "C" G_MODULE_EXPORT 
static bool
on_residue_circle_clicked(GooCanvasItem  *item,
			  GooCanvasItem  *target_item,
			  GdkEventButton *event,
			  gpointer        user_data) { 

   gboolean r = true;

   std::cout << "on_residue_circle_clicked()!!!!! " << std::endl;

   if (target_item) { 
      coot::residue_spec_t *spec_p = (coot::residue_spec_t *) g_object_get_data (G_OBJECT (target_item), "spec");
   
      if (spec_p) { 
      
	 std::cout << "on_residue_circle_clicked(): clicked on " << *spec_p << std::endl;
      
      } else {
	 std::cout << "on_residue_circle_clicked(): null spec" << std::endl;
      }
   } else {
      std::cout << "on_residue_circle_clicked(): NULL target item" << std::endl;
   } 

   return r;
}


// if you don't have an add_rep_handle, then pass -1 (something negative)
// 
void
lbg_info_t::draw_residue_circle_top_layer(const residue_circle_t &residue_circle,
					  const lig_build::pos_t &ligand_centre,
					  int add_rep_handle) {


   if (0)
      std::cout << "   adding cirles " << residue_circle.residue_type
		<< " at init pos " << residue_circle.pos << " and canvas_drag_offset "
		<< canvas_drag_offset << std::endl;

   lig_build::pos_t circle_pos = residue_circle.pos;
      
   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));

   GooCanvasItem *group = goo_canvas_group_new (root, "stroke-color", "#111111",
						NULL);

   // don't draw waters that don't have bonds to the ligand (or the
   // bonds to the ligand atoms are too long, or the water is too far
   // from any protein atom).
   //
   if (residue_circle.residue_type == "HOH") {
      if (residue_circle.bonds_to_ligand.size() == 0) {
	 return;
      }
   }

   // Capitalise the residue type (takes less space than upper case).
   std::string rt = residue_circle.residue_type.substr(0,1);
   rt += coot::util::downcase(residue_circle.residue_type.substr(1));
   
   // correct that if we are looking at dna: DA, DT, DC, DG
   if (residue_circle.residue_type == "DA" ||
       residue_circle.residue_type == "DT" || 
       residue_circle.residue_type == "DC" || 
       residue_circle.residue_type == "DG") {
      rt = "d";
      rt += residue_circle.residue_type.substr(1);
   } 

   // fill colour and stroke colour of the residue circle
   std::pair<std::string, std::string> col = get_residue_circle_colour(residue_circle.residue_type);
   double line_width = 1.0;
   if (col.second != "#111111") // needs checking, FIXME
      line_width = 3.0;
   GooCanvasItem *circle = NULL;
   GooCanvasItem *text_1 = NULL;
   GooCanvasItem *text_2 = NULL;
   

   if (col.first != "") {
      circle = goo_canvas_ellipse_new(group,
				      circle_pos.x, circle_pos.y,
				      standard_residue_circle_radius,
				      standard_residue_circle_radius,
				      "line_width", line_width,
				      "fill-color",   col.first.c_str(),
				      "stroke-color", col.second.c_str(),
				      NULL);
   } else {
      circle = goo_canvas_ellipse_new(group,
				      circle_pos.x, circle_pos.y,
				      standard_residue_circle_radius,
				      standard_residue_circle_radius,
				      "line_width", line_width,
				      "stroke-color", col.second.c_str(),
				      NULL);
   }
   
   text_1 = goo_canvas_text_new(group, rt.c_str(),
				circle_pos.x, circle_pos.y-6, -1,
				GOO_CANVAS_ANCHOR_CENTER,
				"font", "Sans 9",
				"fill_color", "#111111",
				NULL);

   text_2 = goo_canvas_text_new(group, residue_circle.residue_label.c_str(),
				circle_pos.x, circle_pos.y+6.5, -1,
				GOO_CANVAS_ANCHOR_CENTER,
				"font", "Sans 7",
				"fill_color", "#111111",
				NULL);

   if (circle && text_1 && text_2) {
      // so that if a residue is clicked on (or moused over or
      // whatever, then the target_item has some interesting attached
      // data for use in the callback)
      coot::residue_spec_t *sp_p_1 = new coot::residue_spec_t(residue_circle.spec);
      coot::residue_spec_t *sp_p_2 = new coot::residue_spec_t(residue_circle.spec);
      coot::residue_spec_t *sp_p_3 = new coot::residue_spec_t(residue_circle.spec);
      sp_p_1->int_user_data = imol;
      sp_p_2->int_user_data = imol;
      sp_p_3->int_user_data = imol;
      g_object_set_data_full(G_OBJECT(circle), "spec", sp_p_1, g_free);
      g_object_set_data_full(G_OBJECT(text_1), "spec", sp_p_2, g_free);
      g_object_set_data_full(G_OBJECT(text_2), "spec", sp_p_3, g_free);
      int *add_rep_handle_p_1 = new int(add_rep_handle);
      int *add_rep_handle_p_2 = new int(add_rep_handle);
      int *add_rep_handle_p_3 = new int(add_rep_handle);
      g_object_set_data_full(G_OBJECT(circle), "add_rep_handle", add_rep_handle_p_1, g_free);
      g_object_set_data_full(G_OBJECT(text_1), "add_rep_handle", add_rep_handle_p_2, g_free);
      g_object_set_data_full(G_OBJECT(text_2), "add_rep_handle", add_rep_handle_p_3, g_free);
   }
}


// Return the fill colour and the stroke colour.
// 
std::pair<std::string, std::string>
lbg_info_t::get_residue_circle_colour(const std::string &residue_type) const {

   std::string fill_colour = "";
   std::string stroke_colour = "#111111";

   std::string grease = "#ccffbb";
   std::string purple = "#eeccee";
   std::string red    = "#cc0000";
   std::string blue   = "#0000cc";
   std::string metalic_grey = "#d9d9d9";

   if (residue_type == "ALA")
      fill_colour = grease;
   if (residue_type == "TRP")
      fill_colour = grease;
   if (residue_type == "PHE")
      fill_colour = grease;
   if (residue_type == "LEU")
      fill_colour = grease;
   if (residue_type == "PRO")
      fill_colour = grease;
   if (residue_type == "ILE")
      fill_colour = grease;
   if (residue_type == "VAL")
      fill_colour = grease;
   if (residue_type == "MET")
      fill_colour = grease;
   if (residue_type == "MSE")
      fill_colour = grease;

   if (residue_type == "GLY")
      fill_colour = purple;
   if (residue_type == "ASP")
      fill_colour = purple;
   if (residue_type == "ASN")
      fill_colour = purple;
   if (residue_type == "CYS")
      fill_colour = purple;
   if (residue_type == "GLN")
      fill_colour = purple;
   if (residue_type == "GLU")
      fill_colour = purple;
   if (residue_type == "HIS")
      fill_colour = purple;
   if (residue_type == "LYS")
      fill_colour = purple;
   if (residue_type == "LYS")
      fill_colour = purple;
   if (residue_type == "ARG")
      fill_colour = purple;
   if (residue_type == "SER")
      fill_colour = purple;
   if (residue_type == "THR")
      fill_colour = purple;
   if (residue_type == "TYR")
      fill_colour = purple;
   if (residue_type == "HOH")
      fill_colour = "white";

   if (residue_type == "ASP")
      stroke_colour = red;
   if (residue_type == "GLU")
      stroke_colour = red;
   if (residue_type == "LYS")
      stroke_colour = blue;
   if (residue_type == "ARG")
      stroke_colour = blue;
   if (residue_type == "HIS")
      stroke_colour = blue;

   // Bases
   if (residue_type == "U")  fill_colour = purple;
   if (residue_type == "T")  fill_colour = purple;
   if (residue_type == "C")  fill_colour = purple;
   if (residue_type == "A")  fill_colour = purple;
   if (residue_type == "G")  fill_colour = purple;
   if (residue_type == "DT")  fill_colour = purple;
   if (residue_type == "DC")  fill_colour = purple;
   if (residue_type == "DA")  fill_colour = purple;
   if (residue_type == "DG")  fill_colour = purple;
   

   // metals
   if (residue_type == "ZN") 
      fill_colour = metalic_grey;
   if (residue_type == "MG")
      fill_colour = metalic_grey;
   if (residue_type == "NA")
      fill_colour = metalic_grey;
   if (residue_type == "CA")
      fill_colour = metalic_grey;
   if (residue_type == "K")
      fill_colour = metalic_grey;
	 
   return std::pair<std::string, std::string> (fill_colour, stroke_colour);
}




// solvent exposure difference of the residue due to ligand binding
void 
lbg_info_t::draw_solvent_exposure_circle(const residue_circle_t &residue_circle,
					 const lig_build::pos_t &ligand_centre,
					 GooCanvasItem *group) {

   if (residue_circle.residue_type != "HOH") { 
      if (residue_circle.se_diff_set()) {
	 std::pair<double, double> se_pair = residue_circle.solvent_exposures();
	 double radius_extra = (se_pair.second - se_pair.first) * 19;  // was 18, was 14, was 22.
	 if (radius_extra > 0) {
	    lig_build::pos_t to_lig_centre_uv = (ligand_centre - residue_circle.pos).unit_vector();
	    lig_build::pos_t se_circle_centre = residue_circle.pos - to_lig_centre_uv * radius_extra;

	    std::string fill_colour = get_residue_solvent_exposure_fill_colour(radius_extra);
	    double r = standard_residue_circle_radius + radius_extra;

	    GooCanvasItem *circle = goo_canvas_ellipse_new(group,
							   se_circle_centre.x, se_circle_centre.y,
							   r, r,
							   "line_width", 0.0,
							   "fill-color", fill_colour.c_str(),
							   NULL);
	 }
      }
   } 
}

std::string
lbg_info_t::get_residue_solvent_exposure_fill_colour(double r) const {

   std::string colour = "#8080ff";
   double step = 0.7;
   if (r > step)
      colour = "#e0e0ff";
   if (r > step * 2)
      colour = "#d8d8ff";
   if (r > step * 3)
      colour = "#d0d0ff";
   if (r > step * 4)
      colour = "#c0c8ff";
   if (r > step * 5)
      colour = "#b0c0ff";
   if (r > step * 6)
      colour = "#a0b8ff";
   if (r > step * 7)
      colour = "#90b0ff";
   if (r > step * 8)
      colour = "#80a8ff";
   if (r > step * 9)
      colour = "#70a0ff";

   return colour;
} 


std::vector<solvent_accessible_atom_t>
lbg_info_t::read_solvent_accessibilities(const std::string &file_name) const {

   // return this
   std::vector<solvent_accessible_atom_t> solvent_accessible_atoms;
   
   std::ifstream f(file_name.c_str());
   if (!f) {
      std::cout << "Failed to open " << file_name << std::endl;
   } else {

      std::cout << "reading solvent accessibilites file: " << file_name << std::endl;
      
      std::vector<std::string> lines;
      std::string line;
      while (std::getline(f, line)) { 
	 lines.push_back(line);
      }

      for (unsigned int i=0; i<lines.size(); i++) {
	 std::vector<std::string> words = coot::util::split_string_no_blanks(lines[i], " ");
	 if (words.size() > 5) {
	    
	    if (words[0] == "ATOM:") {
	       std::string atom_name = lines[i].substr(5,4);
	       try {
		  double pos_x = lig_build::string_to_float(words[2]);
		  double pos_y = lig_build::string_to_float(words[3]);
		  double pos_z = lig_build::string_to_float(words[4]);
		  double sa    = lig_build::string_to_float(words[5]);
		  clipper::Coord_orth pt(pos_x, pos_y, pos_z);
		  
		  if (0) 
		     std::cout << "got atom name :" << atom_name << ": and pos "
			       << pt.format() << " and accessibility: " << sa
			       << std::endl;
		  solvent_accessible_atom_t saa(atom_name, pt, sa);
		  solvent_accessible_atoms.push_back(saa);
	       }
	       catch (const std::runtime_error &rte) {
		  std::cout << "failed to parse :" << lines[i] << ":" << std::endl;
	       }
	    }
	 }

	 if (words[0] == "BASH:") {
	    if (words.size() == 2) {
	       try {
		  if (words[1] == "unlimited") {
		     if (solvent_accessible_atoms.size())
			solvent_accessible_atoms.back().add_unlimited();
		  } else {
		     if (solvent_accessible_atoms.size()) { 
			double bash_dist = lig_build::string_to_float(words[1]);
			solvent_accessible_atoms.back().add_bash_dist(bash_dist);
		     }
		  }
	       }
	       catch (const std::runtime_error &rte) {
		  std::cout << "failed to parse :" << lines[i] << ":" << std::endl;
	       }
	    } 
	 }
      }
   }
   return solvent_accessible_atoms;
}

void
lbg_info_t::draw_solvent_accessibility_of_atoms() {

   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));
   
   for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
      lig_build::pos_t pos = mol.atoms[iat].atom_position;
      double sa = mol.atoms[iat].get_solvent_accessibility();
      if (sa  > 0) 
	 draw_solvent_accessibility_of_atom(pos, sa, root);
   }
}

void
lbg_info_t::draw_solvent_accessibility_of_atom(const lig_build::pos_t &pos, double sa,
					       GooCanvasItem *root) {

   if (draw_flev_annotations_flag) { 
      int n_circles = int(sa*40) + 1;    // needs fiddling?
      if (n_circles> 10) n_circles = 10; // needs fiddling?
      
      GooCanvasItem *group = goo_canvas_group_new (root, NULL);
      
      for (int i=0; i<n_circles; i++) { 
	 double rad =  LIGAND_TO_CANVAS_SCALE_FACTOR/23.0 * 3.0 * double(i+1); // needs fiddling?
	 GooCanvasItem *circle = goo_canvas_ellipse_new(group,
							pos.x, pos.y,
							rad, rad,
							"line_width", 0.0,
							"fill-color-rgba", 0x5555cc30,
							NULL);

	 goo_canvas_item_lower(group, NULL); // to the bottom
      }
   }
}

void
lbg_info_t::draw_substitution_contour() {

   bool debug = 0;

   if (draw_flev_annotations_flag) { 
      if (mol.atoms.size() > 0) {


	 // first of all, do we have any bash distances for the atoms of this molecule?
	 bool have_bash_distances = 0;
	 for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
	    if (mol.atoms[iat].bash_distances.size()) {
	       have_bash_distances = 1;
	       break;
	    }
	 }

	 // If we don't have bash distances, then don't grid and contour anything.  If
	 // we do, we do....
	 // 
	 if (have_bash_distances) {
      
	    GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));
	    try { 

	       std::pair<lig_build::pos_t, lig_build::pos_t> l_e_pair =
		  mol.ligand_extents();
	       ligand_grid grid(l_e_pair.first, l_e_pair.second);
   
	       if (0) { // debug
		  for (unsigned int i=0; i<mol.atoms.size(); i++) { 
		     std::cout << "in draw_substitution_contour() atom " << i << " "
			       << mol.atoms[i].get_atom_name()
			       << " has "
			       << mol.atoms[i].bash_distances.size() << " bash distances"
			       << std::endl;
		     for (unsigned int j=0; j<mol.atoms[i].bash_distances.size(); j++) {
			std::cout << "  " << mol.atoms[i].bash_distances[j];
		     }
		     if (mol.atoms[i].bash_distances.size())
			std::cout << std::endl;
		  }
	       }

	          

	       std::vector<widgeted_atom_ring_centre_info_t> unlimited_atoms;
	       for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
		  int n_bash_distances = 0;
		  double sum_bash = 0.0;
		  bool unlimited = 0;
		  int n_unlimited = 0;
	       
		  if (mol.atoms[iat].bash_distances.size()) { 
		     for (unsigned int j=0; j<mol.atoms[iat].bash_distances.size(); j++) {
			if (! mol.atoms[iat].bash_distances[j].unlimited()) {
			   sum_bash += mol.atoms[iat].bash_distances[j].dist;
			   n_bash_distances++;
			} else {
			   unlimited = 1;
			   n_unlimited++;
			} 
		     }

		     if (! unlimited) { 
			if (n_bash_distances > 0) {
			   double bash_av = sum_bash/double(n_bash_distances);
			   if (debug)
			      std::cout << "   none unlimited, using bash_av " << bash_av
					<< " for atom " << mol.atoms[iat].get_atom_name()
					<< std::endl;
			   grid.add_for_accessibility(bash_av, mol.atoms[iat].atom_position);
			}
		     } else {

			// Now, were more than half of the bash distances
			// unlimited?  If yes, then this is unlimited.

			if (double(n_unlimited)/double(mol.atoms[iat].bash_distances.size()) > 0.5) { 
	    
			   // just shove some value in to make the grid values vaguely
			   // correct - the unlimited atoms properly assert themselves
			   // in the drawing of the contour (that is, the selection of
			   // the contour fragments).
			   //
			   grid.add_for_accessibility(1.2, mol.atoms[iat].atom_position);
			   // 	       std::cout << "adding unlimited_atom_position "
			   // 			 << iat << " "
			   // 			 << mol.atoms[iat].atom_position

			   // not elegant because no constructor for
			   // widgeted_atom_ring_centre_info_t (because no simple
			   // constructor for lig_build::atom_t).
			   // 
			   widgeted_atom_ring_centre_info_t ua(mol.atoms[iat]);
			   unlimited_atoms.push_back(ua);
			} else {

			   // treat as a limited:
			   double bash_av =
			      (sum_bash + 4.0 * n_unlimited)/double(n_bash_distances+n_unlimited);
			   if (debug)
			      std::cout << "   few unlimited, as limited using bash_av "
					<< bash_av << " for atom "
					<< mol.atoms[iat].get_atom_name()
					<< std::endl;
			
			   grid.add_for_accessibility(bash_av, mol.atoms[iat].atom_position);
			} 
		     }
		  
		  } else {

		     // we don't get here currently, now that there is an
		     // outer test for having bash distances. Current way
		     // is OK, I think.
		  
		     // there were no bash distancs - what do we do?  Leaving out
		     // atoms means gaps over the ligand - bleugh.  Shove some
		     // value in?  1.0?  if they are not hydrogens, of course.
		     // 
		     if (mol.atoms[iat].element != "H") // checked.
			grid.add_for_accessibility_no_bash_dist_atom(1.0, mol.atoms[iat].atom_position);
		  } 
	       }

	       // Put some values around the ring centres too.
	       // 
	       grid.avoid_ring_centres(ring_atoms_list, mol);

	       // for debugging
	       // show_grid(grid);

	       // std::vector<widgeted_atom_ring_centre_info_t> dummy_unlimited_atoms;
	       grid.show_contour(root, 0.5, unlimited_atoms, ring_atoms_list);
	       // debug
	       // show_unlimited_atoms(unlimited_atoms);
	       // show_ring_centres(ring_atoms_list, mol);

	       if (debug) { 
		  std::cout << "Here are the "<< unlimited_atoms.size()
			    << " unlimited atoms: " << std::endl;
		  for (unsigned int iat=0; iat<unlimited_atoms.size(); iat++)
		     std::cout << "   " << unlimited_atoms[iat] << std::endl;
	       }

	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << rte.what() << std::endl;
	    }
	 }
      }
   }
}


// r_squared the squared distance between the atom and the grid point.
//
double
lbg_info_t::ligand_grid::substitution_value(double r_squared, double bash_dist) const {

   double D = bash_dist;
   double r = sqrt(r_squared);

   if (bash_dist < 1) {
      // this is not in C&L, so this is an enhancement of their
      // algorithm.  Without this there is non-smoothness as the
      // contour follows the interface between 0 and 1.

      // So we add a nice smooth slope (similar to that of
      // conventional bash distances (below).  However we add a slope
      // between +/- 0.2 of the bash distance.
      //
      double small = 0.2;
      if (r > (D + small)) { 
	 return 0;
      } else {
	 if (r < (D - small)) {
	    return 1;
	 } else {
	    double m = (r-(D-small))/(2*small);
	    return (0.5 * (1 + cos(M_PI * m)));
	 } 
      }

   } else { 

      if (r<1)
	 return 1; 
      if (r<(D-1)) {
	 return 1;
      } else {
	 if (r > D) { 
	    return 0;
	 } else {
	    // double v = 0.5*(1 + cos(0.5 *M_PI * (D-r))); // C&L - eh? typo/bug
	    double v = 0.5 * (1.0 + cos(M_PI * (D-1-r)));
	    return v;
	 }
      }
   }
}

void
lbg_info_t::draw_stacking_interactions(const std::vector<residue_circle_t> &rc) {

   for (unsigned int ires=0; ires<rc.size(); ires++) {
      int st = rc[ires].get_stacking_type();
      clipper::Coord_orth click_pos = rc[ires].residue_centre_real;
	 
      if (rc[ires].has_ring_stacking_interaction()) {
	 
	 std::vector<std::string> ligand_ring_atom_names =
	    rc[ires].get_ligand_ring_atom_names();
	 if ((st == lbg_info_t::residue_circle_t::PI_PI_STACKING) ||
	     (st == lbg_info_t::residue_circle_t::PI_CATION_STACKING)) {
	    try {
	       lig_build::pos_t lc = mol.get_ring_centre(ligand_ring_atom_names);
	       draw_annotated_stacking_line(lc, rc[ires].pos, st, click_pos);
	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << rte.what() << std::endl;
	    }
	 }
      }
      if (st == lbg_info_t::residue_circle_t::CATION_PI_STACKING) {
	 std::string at_name = rc[ires].get_ligand_cation_atom_name();
	 lig_build::pos_t atom_pos = mol.get_atom_canvas_position(at_name);
	 draw_annotated_stacking_line(atom_pos, rc[ires].pos, st, click_pos);
      }
   }
}

// click_pos is where we recentre in 3D graphics when the annotation
// (line) is clicked.
void
lbg_info_t::draw_annotated_stacking_line(const lig_build::pos_t &ligand_ring_centre,
					 const lig_build::pos_t &residue_pos,
					 int stacking_type,
					 const clipper::Coord_orth &click_pos) {	

   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));

   double ligand_target_shortening_factor = 3;
   if (stacking_type == lbg_info_t::residue_circle_t::CATION_PI_STACKING)
      ligand_target_shortening_factor = 6;
   lig_build::pos_t a_to_b_uv = (ligand_ring_centre - residue_pos).unit_vector();
   lig_build::pos_t A = residue_pos + a_to_b_uv * 20;
   lig_build::pos_t B = ligand_ring_centre;
   // just short of the middle of the ring
   lig_build::pos_t C = B - a_to_b_uv * ligand_target_shortening_factor; 
   lig_build::pos_t mid_pt = (A + B) * 0.5;
   lig_build::pos_t close_mid_pt_1 = mid_pt - a_to_b_uv * 17;
   lig_build::pos_t close_mid_pt_2 = mid_pt + a_to_b_uv * 17;

   gboolean start_arrow = 0;
   gboolean end_arrow = 1;
   std::string stroke_colour = "#008000";
   GooCanvasItem *group = goo_canvas_group_new (root,
						"stroke-color", stroke_colour.c_str(),
						NULL);

   std::vector<lig_build::pos_t> hex_and_ring_centre(2);
   hex_and_ring_centre[0] = mid_pt - a_to_b_uv * 7;
   hex_and_ring_centre[1] = mid_pt + a_to_b_uv * 7;

   // for cations, we draw a 6-member ring and a "+" text.
   //
   // for pi-pi, we draw 2 6-member rings
   // 
   for(int ir=0; ir<2; ir++) {

      bool do_ring = 0;
      bool do_anion = 0;
      
      if (stacking_type == lbg_info_t::residue_circle_t::PI_PI_STACKING) {
	 do_ring = 1; 
      } 
      if (stacking_type == lbg_info_t::residue_circle_t::PI_CATION_STACKING) {
	 if (ir == 0) {
	    do_ring = 0;
	    do_anion = 1;
	 } else {
	    do_ring = 1;
	    do_anion = 0;
	 }
      }
      if (stacking_type == lbg_info_t::residue_circle_t::CATION_PI_STACKING) {
	 if (ir == 0) {
	    do_ring = 1;
	    do_anion = 0;
	 } else {
	    do_ring = 0;
	    do_anion = 1;
	 }
      }
      
      if (do_ring) { 
	 double angle_step = 60;
	 double r = 8; // radius
	 for (int ipt=1; ipt<=6; ipt++) {
	    int ipt_0 = ipt - 1;
	    double theta_deg_1 = 30 + angle_step * ipt;
	    double theta_1 = theta_deg_1 * DEG_TO_RAD;
	    lig_build::pos_t pt_1 =
	       hex_and_ring_centre[ir] + lig_build::pos_t(sin(theta_1), cos(theta_1)) * r;
      
	    double theta_deg_0 = 30 + angle_step * ipt_0;
	    double theta_0 = theta_deg_0 * DEG_TO_RAD;
	    lig_build::pos_t pt_0 =
	       hex_and_ring_centre[ir] + lig_build::pos_t(sin(theta_0), cos(theta_0)) * r;
	    GooCanvasItem *line = goo_canvas_polyline_new_line(group,
							       pt_1.x, pt_1.y,
							       pt_0.x, pt_0.y,
							       "line_width", 1.8,
							       NULL);
	    clipper::Coord_orth *pos_l = new clipper::Coord_orth(click_pos);
	    g_object_set_data_full(G_OBJECT(line), "position", pos_l, g_free);
	 }
	 // Now the circle in the annotation aromatic ring:
	 GooCanvasItem *ring = goo_canvas_ellipse_new(group,
						      hex_and_ring_centre[ir].x,
						      hex_and_ring_centre[ir].y,
						      4.0, 4.0,
						      "line_width", 1.0,
						      NULL);
	 clipper::Coord_orth *pos_r = new clipper::Coord_orth(click_pos);
	 g_object_set_data_full(G_OBJECT(ring), "position", pos_r, g_free);
      }
      
      if (do_anion) {
	 // the "+" symbol for the anion
	 // 
	 GooCanvasItem *text_1 = goo_canvas_text_new(group,
						     "+",
						     hex_and_ring_centre[ir].x,
						     hex_and_ring_centre[ir].y,
						     -1,
						     GOO_CANVAS_ANCHOR_CENTER,
						     "font", "Sans 12",
						     "fill_color", stroke_colour.c_str(),
						     NULL);
	 clipper::Coord_orth *pos_p = new clipper::Coord_orth(click_pos);
	 g_object_set_data_full(G_OBJECT(text_1), "position", pos_p, g_free);
      }
   }


   GooCanvasLineDash *dash = goo_canvas_line_dash_new (2, 2.5, 2.5);
   GooCanvasItem *item_1 =
      goo_canvas_polyline_new_line(group,
				   A.x, A.y,
				   close_mid_pt_1.x, close_mid_pt_1.y,
				   "line-width", 2.5, // in draw_annotated_stacking_line()
				   "line-dash", dash,
				   "stroke-color", stroke_colour.c_str(),
				   NULL);

   GooCanvasItem *item_2 =
      goo_canvas_polyline_new_line(group,
				   close_mid_pt_2.x, close_mid_pt_2.y,
				   C.x, C.y,
				   "line-width", 2.5, // in draw_annotated_stacking_line()
				   "line-dash", dash,
				   // "end_arrow",   end_arrow,
				   "stroke-color", stroke_colour.c_str(),
				   NULL);

   clipper::Coord_orth *pos_1 = new clipper::Coord_orth(click_pos);
   clipper::Coord_orth *pos_2 = new clipper::Coord_orth(click_pos);
   g_object_set_data_full(G_OBJECT(item_1), "position", pos_1, g_free);
   g_object_set_data_full(G_OBJECT(item_2), "position", pos_2, g_free);
   
   // Now the circle blob at the centre of the aromatic ligand ring:
   if (stacking_type != lbg_info_t::residue_circle_t::CATION_PI_STACKING) { 
      GooCanvasItem *item_o = goo_canvas_ellipse_new(group,
						     B.x, B.y,
						     3.0, 3.0,
						     "line_width", 1.0,
						     "fill_color", stroke_colour.c_str(),
						     NULL);
      clipper::Coord_orth *pos_o = new clipper::Coord_orth(click_pos);
      g_object_set_data_full(G_OBJECT(item_o), "position", pos_o, g_free);
   }

   clipper::Coord_orth *pos_p = new clipper::Coord_orth(click_pos);
   g_object_set_data_full(G_OBJECT(group), "position", pos_p, g_free);


}

void
lbg_info_t::draw_bonds_to_ligand() {

   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));

   GooCanvasLineDash *dash_dash     = goo_canvas_line_dash_new (2, 2.1, 2.1);
   GooCanvasLineDash *dash_unbroken = goo_canvas_line_dash_new (1, 3);

   GooCanvasLineDash *dash = dash_dash; // re-assigned later (hopefully)
   for (unsigned int ic=0; ic<residue_circles.size(); ic++) { 
      if (residue_circles[ic].bonds_to_ligand.size()) {
	 for (unsigned int ib=0; ib<residue_circles[ic].bonds_to_ligand.size(); ib++) { 

	    try {

	       lig_build::pos_t pos = residue_circles[ic].pos;
	       std::string at_name = residue_circles[ic].bonds_to_ligand[ib].ligand_atom_name;
	       lig_build::pos_t lig_atom_pos = mol.get_atom_canvas_position(at_name);
	       lig_build::pos_t rc_to_lig_atom = lig_atom_pos - pos;
	       lig_build::pos_t rc_to_lig_atom_uv = rc_to_lig_atom.unit_vector();
	       lig_build::pos_t B = lig_atom_pos;

	       if (residue_circles[ic].bonds_to_ligand[ib].bond_type != coot::fle_ligand_bond_t::BOND_COVALENT)
		  B -= rc_to_lig_atom_uv * 8; // put the end point (say) of a hydrogen bond
		                              // some pixels away from the atom centre - for better
              		                      // aesthetics.
		  
	       lig_build::pos_t A = pos + rc_to_lig_atom_uv * 20;

	       // some colours
	       std::string blue = "blue";
	       std::string green = "darkgreen";
	       std::string lime = "#888820";
	       std::string stroke_colour = "#111111"; // unset

	       // arrows (acceptor/donor) and stroke colour (depending
	       // on mainchain or sidechain interaction)
	       // 
	       gboolean start_arrow = 0;
	       gboolean   end_arrow = 0;
	       if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::H_BOND_DONOR_SIDECHAIN) {
		  end_arrow = 1;
		  stroke_colour = green;
		  dash = dash_dash;
	       }
	       if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::H_BOND_DONOR_MAINCHAIN) {
		  end_arrow = 1;
		  stroke_colour = blue;
		  dash = dash_dash;
	       }
	       if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::H_BOND_ACCEPTOR_SIDECHAIN) {
		  start_arrow = 1;
		  stroke_colour = green;
		  dash = dash_dash;
	       }
	       if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::H_BOND_ACCEPTOR_MAINCHAIN) {
		  start_arrow = 1;
		  stroke_colour = blue;
		  dash = dash_dash;
	       }
	       if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::METAL_CONTACT_BOND) {
		  stroke_colour = "#990099";
		  dash = dash_dash;
	       }
	       if (residue_circles[ic].bonds_to_ligand[ib].bond_type == bond_to_ligand_t::BOND_COVALENT) {
		  dash = goo_canvas_line_dash_new (2, 2.7, 0.1);
		  stroke_colour = "#bb00bb";
	       }
	       
	       if (residue_circles[ic].residue_type == "HOH") { 
		  stroke_colour = lime;
		  start_arrow = 0;
		  end_arrow = 0;
		  dash = dash_dash;
	       }
	       GooCanvasItem *item = goo_canvas_polyline_new_line(root,
								  A.x, A.y,
								  B.x, B.y,
								  "line-width", 2.5, // in draw_bonds_to_ligand()
								  "line-dash", dash,
 								  "start_arrow", start_arrow,
 								  "end_arrow",   end_arrow,
								  "stroke-color", stroke_colour.c_str(),
								  NULL);
	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << "WARNING:: " << rte.what() << std::endl;
	    }
	 }
	 
      } else {
	 if (0) 
	    std::cout << "... no bond to ligand from residue circle "
		      << residue_circles[ic].residue_label << " "
		      << residue_circles[ic].residue_type << std::endl;
      }
   }
}

double
lbg_info_t::bottom_of_flev_items() {

   double biggest_y = 0;
   for (unsigned int ic=0; ic<residue_circles.size(); ic++) {
      // std::cout << "   residue cirlce " << ic << " " << residue_circles[ic].pos << std::endl;
      if (residue_circles[ic].pos.y > biggest_y)
	 biggest_y = residue_circles[ic].pos.y;
   }
   // also consider atom position here - for the (very?) rare case
   // where the ligand atoms are below residue circles?
   
   return biggest_y;
} 


// This is for reading SMILES, to use elsewhere, the lines
// need to have newlines added.
std::string
lbg_info_t::file_to_string(const std::string &file_name) const {

   std::string s;
   std::string line;
   std::ifstream f(file_name.c_str());
   if (!f) {
      std::cout << "Failed to open " << file_name << std::endl;
   } else {
      while (std::getline(f, line)) { 
	 s += line;
      }
   }
   return s;
}

void
lbg_info_t::show_key() {

   GooCanvasItem *root = goo_canvas_get_root_item (GOO_CANVAS(canvas));

   // class variable.
   if (! key_group)
      key_group = goo_canvas_group_new (root, NULL);

   GooCanvasLineDash *dash = goo_canvas_line_dash_new (2, 2.5, 2.5);

   double dy  =  30.0;
   double dxc = 250.0; // column offset
   lig_build::pos_t co(dxc, 0);
   lig_build::pos_t ro(0, dy);
   double bf = bottom_of_flev_items();
   lig_build::pos_t tl(100, bf+100); // top left of key

   // Lines
   lig_build::pos_t A       = tl;
   lig_build::pos_t B       = tl + lig_build::pos_t(40,  0);
   lig_build::pos_t AB_txt  = tl + lig_build::pos_t(50,  0);
   lig_build::pos_t Ad      = A +  ro;
   lig_build::pos_t Bd      = B +  ro;
   lig_build::pos_t ABd_txt = AB_txt + ro;
   lig_build::pos_t C       = tl + lig_build::pos_t(0,  2*dy);
   lig_build::pos_t D       = tl + lig_build::pos_t(40, 2*dy);
   lig_build::pos_t CD_txt  = tl + lig_build::pos_t(50, 2*dy);
   lig_build::pos_t Cd      = C +  lig_build::pos_t(0,  dy);
   lig_build::pos_t Dd      = D +  lig_build::pos_t(0,  dy);
   lig_build::pos_t CDd_txt = CD_txt + lig_build::pos_t(0,  dy);

   lig_build::pos_t WatA    = tl + ro * 4;
   lig_build::pos_t WatB    = WatA + lig_build::pos_t(40,  0);
   lig_build::pos_t Wat_txt = WatA + lig_build::pos_t(50,  0);
   
   lig_build::pos_t MetalA    = tl + ro * 5;
   lig_build::pos_t MetalB    = MetalA + lig_build::pos_t(40,  0);
   lig_build::pos_t Metal_txt = MetalA + lig_build::pos_t(50,  0);
   // Circles
   lig_build::pos_t E       = tl + co;
   lig_build::pos_t F       = E + ro;
   lig_build::pos_t G       = E + ro * 2;
   lig_build::pos_t H       = E + ro * 3;
   lig_build::pos_t I       = E + ro * 4;
   lig_build::pos_t J       = E + ro * 5;
   lig_build::pos_t Et      = E + lig_build::pos_t(25, 0);
   lig_build::pos_t Ft      = F + lig_build::pos_t(25, 0);
   lig_build::pos_t Gt      = G + lig_build::pos_t(25, 0);
   lig_build::pos_t Ht      = H + lig_build::pos_t(25, 0);
   lig_build::pos_t It      = I + lig_build::pos_t(25, 0);
   lig_build::pos_t Jt      = J + lig_build::pos_t(25, 0);
   // Rest
   lig_build::pos_t K       = tl + co + lig_build::pos_t(150, 0);
   lig_build::pos_t L       = K + ro * 1.5;
   lig_build::pos_t M       = K + ro * 3;
   lig_build::pos_t Kt      = K + lig_build::pos_t(25, 0);
   lig_build::pos_t Lt      = L + lig_build::pos_t(25, 0);
   lig_build::pos_t Mt      = M + lig_build::pos_t(25, 0);
   
   std::string stroke_colour = "#008000";
   stroke_colour = "#0000cc"; // blue
   gboolean end_arrow = 1;
   GooCanvasItem *item_1 =
      goo_canvas_polyline_new_line(key_group,
				   A.x, A.y,
				   B.x, B.y,
				   "line-width", 2.5, // in show_key()
				   "line-dash", dash,
				   "stroke-color", stroke_colour.c_str(),
				   "end_arrow",   end_arrow,
				   NULL);
   GooCanvasItem *text_1 = goo_canvas_text_new(key_group,
					       "Main-chain acceptor",
					       AB_txt.x, AB_txt.y,
					       -1,
					       GOO_CANVAS_ANCHOR_WEST,
					       "font", "Sans 10",
					       "fill_color", "#111111",
					       NULL);

   GooCanvasItem *item_2 =
      goo_canvas_polyline_new_line(key_group,
				   Ad.x, Ad.y,
				   Bd.x, Bd.y,
				   "line-width", 2.5, // in show_key()
				   "line-dash", dash,
				   "stroke-color", stroke_colour.c_str(),
				   "start_arrow",   end_arrow,
				   NULL);
   
   GooCanvasItem *text_2 = goo_canvas_text_new(key_group,
					       "Main-chain donor",
					       ABd_txt.x, ABd_txt.y,
					       -1,
					       GOO_CANVAS_ANCHOR_WEST,
					       "font", "Sans 10",
					       "fill_color", "#111111",
					       NULL);
   
   stroke_colour = "#008000"; // green
   GooCanvasItem *item_3 =
      goo_canvas_polyline_new_line(key_group,
				   C.x, C.y,
				   D.x, D.y,
				   "line-width", 2.5, // in show_key()
				   "line-dash", dash,
				   "stroke-color", stroke_colour.c_str(),
				   "end_arrow",   end_arrow,
				   NULL);

   GooCanvasItem *text_3 = goo_canvas_text_new(key_group,
					       "Side-chain acceptor",
					       CD_txt.x, CD_txt.y,
					       -1,
					       GOO_CANVAS_ANCHOR_WEST,
					       "font", "Sans 10",
					       "fill_color", "#111111",
					       NULL);
   
   GooCanvasItem *item_4 =
      goo_canvas_polyline_new_line(key_group,
				   Cd.x, Cd.y,
				   Dd.x, Dd.y,
				   "line-width", 2.5, // in show_key()
				   "line-dash", dash,
				   "stroke-color", stroke_colour.c_str(),
				   "start_arrow",   end_arrow,
				   NULL);

   GooCanvasItem *text_4 = goo_canvas_text_new(key_group,
					       "Side-chain donor",
					       CDd_txt.x, CDd_txt.y,
					       -1,
					       GOO_CANVAS_ANCHOR_WEST,
					       "font", "Sans 10",
					       "fill_color", "#111111",
					       NULL);

   std::string lime = "#888820";
   GooCanvasItem *item_5 =
      goo_canvas_polyline_new_line(key_group,
				   WatA.x, WatA.y,
				   WatB.x, WatB.y,
				   "line-width", 2.5, // in show_key()
				   "line-dash", dash,
				   "stroke-color", lime.c_str(),
				   NULL);

   GooCanvasItem *text_5 = goo_canvas_text_new(key_group,
					       "H-bond to Water",
					       Wat_txt.x, Wat_txt.y,
					       -1,
					       GOO_CANVAS_ANCHOR_WEST,
					       "font", "Sans 10",
					       "fill_color", "#111111",
					       NULL);
   GooCanvasItem *item_6 =
      goo_canvas_polyline_new_line(key_group,
				   MetalA.x, MetalA.y,
				   MetalB.x, MetalB.y,
				   "line-width", 2.5, // in show_key()
				   "line-dash", dash,
				   "stroke-color",  "#990099",
				   NULL);

   GooCanvasItem *text_6 = goo_canvas_text_new(key_group,
					       "Metal Bond",
					       Metal_txt.x, Metal_txt.y,
					       -1,
					       GOO_CANVAS_ANCHOR_WEST,
					       "font", "Sans 10",
					       "fill_color", "#111111",
					       NULL);

   
   // circles
   double key_residue_radius = 13;
   double line_width = 1.0;
   std::string grease = "#ccffbb";
   std::string purple = "#eeccee";
   std::string red    = "#cc0000";
   std::string blue   = "#0000cc";
   std::string metalic_grey = "#d9d9d9";

   GooCanvasItem *circle_grease = 
      goo_canvas_ellipse_new(key_group,
			     E.x, E.y, key_residue_radius, key_residue_radius,
			     "line_width", line_width,
			     "fill-color", grease.c_str(),
			     NULL);
   GooCanvasItem *circle_polar = 
      goo_canvas_ellipse_new(key_group,
			     F.x, F.y, key_residue_radius, key_residue_radius,
			     "line_width", line_width,
			     "fill-color", purple.c_str(),
			     NULL);
   GooCanvasItem *circle_acidic = 
      goo_canvas_ellipse_new(key_group,
			     G.x, G.y, key_residue_radius, key_residue_radius,
			     "stroke-color", red.c_str(),
			     "fill-color", purple.c_str(),
			     "line_width", 2.0,
			     NULL);
   GooCanvasItem *circle_basic = 
      goo_canvas_ellipse_new(key_group,
			     H.x, H.y, key_residue_radius, key_residue_radius,
			     "fill-color", purple.c_str(),
			     "stroke-color", blue.c_str(),
			     "line_width", 2.0,
			     NULL);
   GooCanvasItem *circle_water = 
      goo_canvas_ellipse_new(key_group,
			     I.x, I.y, key_residue_radius, key_residue_radius,
			     "line_width", line_width,
			     "fill-color", "white",
			     NULL);
   GooCanvasItem *circle_metal = 
      goo_canvas_ellipse_new(key_group,
			     J.x, J.y, key_residue_radius, key_residue_radius,
			     "line_width", line_width,
			     "fill-color", metalic_grey.c_str(),
			     NULL);

   GooCanvasItem *text_grease = goo_canvas_text_new(key_group,
						    "Grease",
						    Et.x, Et.y,
						    -1,
						    GOO_CANVAS_ANCHOR_WEST,
						    "font", "Sans 10",
						    "fill_color", "#111111",
						    NULL);
   GooCanvasItem *text_polar = goo_canvas_text_new(key_group,
						   "Polar",
						   Ft.x, Ft.y,
						   -1,
						   GOO_CANVAS_ANCHOR_WEST,
						   "font", "Sans 10",
						   "fill_color", "#111111",
						   NULL);
   GooCanvasItem *text_acidic = goo_canvas_text_new(key_group,
						   "Acidic",
						   Gt.x, Gt.y,
						   -1,
						   GOO_CANVAS_ANCHOR_WEST,
						   "font", "Sans 10",
						   "fill_color", "#111111",
						   NULL);
   GooCanvasItem *text_basic = goo_canvas_text_new(key_group,
						   "Basic",
						   Ht.x, Ht.y,
						   -1,
						   GOO_CANVAS_ANCHOR_WEST,
						   "font", "Sans 10",
						   "fill_color", "#111111",
						   NULL);
   GooCanvasItem *text_water = goo_canvas_text_new(key_group,
						   "Water",
						   It.x, It.y,
						   -1,
						   GOO_CANVAS_ANCHOR_WEST,
						   "font", "Sans 10",
						   "fill_color", "#111111",
						   NULL);
   
   GooCanvasItem *text_metal = goo_canvas_text_new(key_group,
						   "Metal",
						   Jt.x, Jt.y,
						   -1,
						   GOO_CANVAS_ANCHOR_WEST,
						   "font", "Sans 10",
						   "fill_color", "#111111",
						   NULL);

   // Solvent accessibility of atoms
   bool save_state = draw_flev_annotations_flag;
   draw_flev_annotations_flag = true;
   lig_build::pos_t K_dx = K + lig_build::pos_t(-4, 0); // K + delta x to make objects line up...
   draw_solvent_accessibility_of_atom(K_dx, 0.15, key_group);
   draw_flev_annotations_flag = save_state;
   GooCanvasItem *text_solvent_atom = goo_canvas_text_new(key_group,
							  "Solvent exposure",
							  Kt.x, Kt.y,
							  -1,
							  GOO_CANVAS_ANCHOR_WEST,
							  "font", "Sans 10",
							  "fill_color", "#111111",
							  NULL);
   
   // solvent exposue ligand difference on environment/receptor residue
   GooCanvasItem *circle_residue_exposure_1 = 
      goo_canvas_ellipse_new(key_group,
			     L.x-5.0, L.y-5.0, key_residue_radius * 1.3 , key_residue_radius * 1.3,
			     "line_width", 0.0,
			     "fill-color", "#bbbbff",
			     NULL);
   GooCanvasItem *circle_residue_exposure_2 = 
      goo_canvas_ellipse_new(key_group,
			     L.x, L.y, key_residue_radius * 0.9, key_residue_radius * 0.9,
			     "line_width", line_width,
			     "fill-color", "white",
			     NULL);
   GooCanvasItem *text_solvent_residue = goo_canvas_text_new(key_group,
							     "Residue protection",
							     Lt.x, Lt.y,
							     -1,
							     GOO_CANVAS_ANCHOR_WEST,
							     "font", "Sans 10",
							     "fill_color", "#111111",
							     NULL);

   std::vector<std::vector<lig_build::pos_t> > contour_lines;
   std::vector<lig_build::pos_t> contour_set;

   for (unsigned int i=0; i<11; i++) {
      double delta_x = 15 * sin(i * 2 * M_PI * 0.1);
      double delta_y = 12 * cos(i * 2 * M_PI * 0.1);
      lig_build::pos_t p = M + lig_build::pos_t(delta_x-3, delta_y);
      contour_set.push_back(p);
   }
   contour_lines.push_back(contour_set);
   lig_build::pos_t dummy_pos;
   ligand_grid grid(dummy_pos, dummy_pos);
   grid.plot_contour_lines(contour_lines, key_group);

   GooCanvasItem *text_substitution_contour = goo_canvas_text_new(key_group,
								  "Substitution Contour",
								  Mt.x, Mt.y,
								  -1,
								  GOO_CANVAS_ANCHOR_WEST,
								  "font", "Sans 10",
								  "fill_color", "#111111",
								  NULL);
}


void
lbg_info_t::hide_key() {

   if (key_group) { 
      gint n_children = goo_canvas_item_get_n_children (key_group);
      for (int i=0; i<n_children; i++)
	 goo_canvas_item_remove_child(key_group, 0);
      key_group = NULL;
   }
}

bool
lbg_info_t::annotate(const std::vector<std::pair<coot::atom_spec_t, float> > &s_a_v,
		     const std::vector<coot::fle_residues_helper_t> &centres,
		     const std::vector<int> &additional_representation_handles_in,
		     const std::vector<coot::fle_ligand_bond_t> &bonds_to_ligand,
		     const std::vector<coot::solvent_exposure_difference_helper_t> &sed,
		     const coot::flev_attached_hydrogens_t &ah,
		     const coot::pi_stacking_container_t &pi_stack_info,
		     const coot::dictionary_residue_restraints_t &restraints) {


   if (0) {
      for (unsigned int ib=0; ib<bonds_to_ligand.size(); ib++) { 
	 std::cout << "  =============== lbg::annotate() bond to ligand " << ib << " "
		   << bonds_to_ligand[ib].ligand_atom_spec.atom_name
		   << " by residue " << bonds_to_ligand[ib].res_spec.chain_id << " "
		   << bonds_to_ligand[ib].res_spec.res_no << " type: "
		   << bonds_to_ligand[ib].bond_type
		   << std::endl;
      }
   }

   if (0) { 
      std::cout << "--------------------------------------------------------------" << std::endl;
      std::cout << "======== lbg_info_t::annotate() here are bash distances for atoms:" << std::endl;
      std::map<std::string, std::vector<coot::bash_distance_t> >::const_iterator it;
      for (it=ah.atom_bashes.begin(); it!=ah.atom_bashes.end(); it++) {
	 std::cout << it->first << " " << it->second.size() << " bashes " << std::endl;
	 for (unsigned int i=0; i<it->second.size(); i++) { 
	    std::cout << "   " << it->second[i] << std::endl;
	 }
      }
      std::cout << "--------------------------------------------------------------" << std::endl;
   }

   bool r = false;
   std::vector<solvent_accessible_atom_t> solvent_accessible_atoms =
      convert(s_a_v, ah);
   widgeted_molecule_t new_mol = mol;
   new_mol.map_solvent_accessibilities_to_atoms(solvent_accessible_atoms);

   // fill class data item std::vector<residue_circle_t> residue_circles;
   // 
   residue_circles.clear();
   for (unsigned int i=0; i<centres.size(); i++) {

      if (0)
	 std::cout << "debugg:: in lbg_info_t::annotate() handling circle " << i << " of "
		   << centres.size() << std::endl;
      
      std::string label = centres[i].spec.chain_id;
      label += coot::util::int_to_string(centres[i].spec.res_no);
      label += centres[i].spec.ins_code;
      
      clipper::Coord_orth cp(centres[i].transformed_relative_centre.x(),
			     centres[i].transformed_relative_centre.y(),
			     centres[i].transformed_relative_centre.z());

      residue_circle_t circle(cp,
			      centres[i].interaction_position,
			      centres[i].spec,
			      centres[i].residue_name,
			      label);
      
      lig_build::pos_t pos = mol.input_coords_to_canvas_coords(cp);
      circle.set_canvas_pos(pos);
      if (centres[i].residue_name == "HOH") {
	 for (unsigned int ib=0; ib<bonds_to_ligand.size(); ib++) { 
	    if (bonds_to_ligand[ib].res_spec == centres[i].spec) {
	       circle.set_water_dist_to_protein(bonds_to_ligand[ib].water_protein_length);
	    }
	 }
      }

      // now add any H-bonds to the ligand:
      // 
      for (unsigned int ib=0; ib<bonds_to_ligand.size(); ib++) { 
	 if (bonds_to_ligand[ib].res_spec == centres[i].spec) {
	    double bond_l = bonds_to_ligand[ib].bond_length;
	    std::string ligand_atom_name = bonds_to_ligand[ib].ligand_atom_spec.atom_name;
	    lbg_info_t::bond_to_ligand_t btl(ligand_atom_name, bond_l);
	    btl.bond_type = bonds_to_ligand[ib].bond_type;
	    circle.add_bond_to_ligand(btl);
	 }
      }

      // Now add any pi-stacking interactions to/from the ligand:
      //
      for (unsigned int istack=0; istack<pi_stack_info.stackings.size(); istack++) {
	 coot::residue_spec_t spec(pi_stack_info.stackings[istack].res);
	 if (spec == centres[i].spec) {
	    if (pi_stack_info.stackings[istack].type == coot::pi_stacking_instance_t::PI_PI_STACKING) {
	       std::vector<std::string> lra = pi_stack_info.stackings[istack].ligand_ring_atom_names;
	       circle.set_stacking("pi-pi", lra, "");
	    }
	    if (pi_stack_info.stackings[istack].type == coot::pi_stacking_instance_t::PI_CATION_STACKING) {
	       std::vector<std::string> lra = pi_stack_info.stackings[istack].ligand_ring_atom_names;
	       circle.set_stacking("pi-cation", lra, "");
	    }
	    // cation in ligand, PHE (say)
	    if (pi_stack_info.stackings[istack].type == coot::pi_stacking_instance_t::CATION_PI_STACKING) {
	       std::vector<std::string> lra_null;
	       circle.set_stacking("cation-pi", lra_null, pi_stack_info.stackings[istack].ligand_cationic_atom_name);
	    }
	 }
      }

      // solvent exposure difference annotation.
      //
      // std::vector<coot::solvent_exposure_difference_helper_t> &sed,
      //
      for (unsigned int ised=0; ised<sed.size(); ised++) { 
	 if (sed[ised].res_spec == centres[i].spec) {
	    circle.set_solvent_exposure_diff(sed[ised].exposure_fraction_holo,
					     sed[ised].exposure_fraction_apo);
	 }
      }
      residue_circles.push_back(circle);
   }

   // additional_representation_handles transfer
   additional_representation_handles = additional_representation_handles_in;


   // ligand ring centres (some of which may be aromatic and are in pi_stack_info too).
   //
   // However, for the substitution contour and the initial layout
   // ring avoidance, we blob in a circle of radius 1/(2*sin(180/n_ring_atoms)) bond lengths.

   // sets the object variable
   ring_atoms_list = restraints.get_ligand_ring_list();
   
   // just checking that it was passed correctly - 
   //
   if (0) {
      std::cout << "------------------- ring list ------------" << std::endl;
      for (unsigned int i=0; i<ring_atoms_list.size(); i++) {
	 std::cout << "ring list " << i << "   ";
	 for (unsigned int j=0; j<ring_atoms_list[i].size(); j++) { 
	 std::cout << ring_atoms_list[i][j] << "  ";
	 }
	 std::cout << std::endl;
      }
   }
   
   import_from_widgeted_molecule(new_mol);
   render();
   refine_residue_circle_positions();

   // has the current solution problems due to residues too close to the ligand?
   std::pair<bool, std::vector<int> > problem_status = solution_has_problems_p();
   // std::cout << "::::::::: problem status: " << problem_status.first << std::endl;
   for (unsigned int ip=0; ip<problem_status.second.size(); ip++) {
      std::cout << ":::::::::::: "
		<< residue_circles[problem_status.second[ip]].residue_label << " "
		<< residue_circles[problem_status.second[ip]].residue_type << " "
		<< residue_circles[problem_status.second[ip]].pos
		<< std::endl;
   }
   
   if (problem_status.second.size()) {
      // fiddle with residue_circles and reoptimise.
      std::vector<int> primary_indices = get_primary_indices();
      reposition_problematics_and_reoptimise(problem_status.second, primary_indices);
   }

   recentre_considering_residue_centres();
   return r;
}


std::vector<solvent_accessible_atom_t>
lbg_info_t::convert(const std::vector<std::pair<coot::atom_spec_t, float> > &s_a_v,
		    const coot::flev_attached_hydrogens_t &ah) const {

   std::vector<solvent_accessible_atom_t> r(s_a_v.size());

   for (unsigned int i=0; i<s_a_v.size(); i++) {
      std::string name = s_a_v[i].first.atom_name;
      r[i].atom_name = name;
      r[i].solvent_accessibility = s_a_v[i].second;

      std::map<std::string, std::vector<coot::bash_distance_t> >::const_iterator it =
	 ah.atom_bashes.find(name);
      if (it != ah.atom_bashes.end()) 
	 r[i].bash_distances = it->second;
   }

   return r;
} 

#endif // HAVE_GOOCANVAS
