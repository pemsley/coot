/* lbg/residue-circles.cc
 * 
 * Author: Paul Emsley
 * Copyright 2010 by The University of Oxford
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
#endif

#include "flev.hh"
#include "optimize-residue-circles.hh"
#include "flev-scale-factor.hh"

bool pli::optimise_residue_circles::score_vs_ligand_atoms       = 1;
bool pli::optimise_residue_circles::score_vs_ring_centres       = 1;
bool pli::optimise_residue_circles::score_vs_other_residues     = 1;
bool pli::optimise_residue_circles::score_vs_original_positions = 1;
bool pli::optimise_residue_circles::score_vs_other_residues_for_angles = 0;
bool pli::optimise_residue_circles::score_vs_ligand_atom_bond_length = 1;


pli::optimise_residue_circles::optimise_residue_circles(const std::vector<residue_circle_t> &r,
							       const std::vector<residue_circle_t> &c,
							       const svg_molecule_t &mol_in,
							       const std::vector<int> &primary_indices_in) {

   bool show_dynamics = 0;

   mol = mol_in;
   starting_circles = r;
   current_circles = c;
   primary_indices = primary_indices_in;

   setup_angles(); // for many residues bound to the same ligand atom.
		   // Uses current_circles and mol.

   if (r.size() == 0)
      return;

   gsl_multimin_function_fdf my_func;

   int n_var = starting_circles.size()*2;
   my_func.n = n_var;  /* number of function components */
   my_func.f   = &optimise_residue_circles::f;
   my_func.df  = &optimise_residue_circles::df;
   my_func.fdf = &optimise_residue_circles::fdf;
   my_func.params = static_cast<void *> (this);

   // set the starting points
   gsl_vector *x = gsl_vector_alloc(n_var);

   for (unsigned int i=0; i<starting_circles.size(); i++) {
      gsl_vector_set(x, 2*i,   current_circles[i].pos.x);
      gsl_vector_set(x, 2*i+1, current_circles[i].pos.y);
   }

   gsl_multimin_fdfminimizer *s;
   const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_pr;

   s = gsl_multimin_fdfminimizer_alloc (T, n_var);
   gsl_multimin_fdfminimizer_set (s, &my_func, x, 1, 1e-4);
   size_t iter = 0;
   size_t n_steps = 500; // increase? (was 400) trying 500
   if (show_dynamics)
      n_steps = 60;

   do {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate(s);

      if (status)
	 break;

      status = gsl_multimin_test_gradient (s->gradient, 2e-3);

      if (status == GSL_SUCCESS)
	 std::cout << "INFO:: Success:: Minimum found at iter " << iter << "\n";

   } while (status == GSL_CONTINUE && iter < n_steps);

   if (false)
      std::cout << "============= terminating at iter " << iter << " with status "
		<< status;
   if (status == GSL_ENOPROG)
      std::cout << " NO_PROGRESS";

   for (unsigned int i=0; i<current_circles.size(); i++) {
      current_circles[i].pos.x = gsl_vector_get(s->x, 2*i  );
      current_circles[i].pos.y = gsl_vector_get(s->x, 2*i+1);
   }

   gsl_multimin_fdfminimizer_free (s);
   gsl_vector_free (x);

}

std::vector<int>
flev_t::get_primary_indices() const {

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

void
flev_t::initial_residues_circles_layout() {

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
      ligand_grid grid(l_e_pair.first, l_e_pair.second);
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


void
flev_t::refine_residue_circle_positions() { // changes the positions of residue_circles

   std::cout << "--------------- refine_residue_circle_positions() --- start --- " << std::endl;

   std::vector<int> primary_indices = get_primary_indices();
   initial_residues_circles_layout(); // twiddle residue_circles
   std::vector<residue_circle_t> current_circles = residue_circles;

   for (int iround=0; iround<30; iround++) {
      std::cout << "flev_t::refine_residue_circle_positions(): iround      " << iround << std::endl;
      std::pair<int, std::vector<residue_circle_t> > new_c =
	 optimise_residue_circle_positions(residue_circles, current_circles, primary_indices);
      current_circles = new_c.second;
      if (new_c.first == GSL_ENOPROG)
	 break;
      if (new_c.first == GSL_SUCCESS) {
	 break;
      }
   }
   residue_circles = current_circles;
}

std::pair<bool,lig_build::pos_t>
flev_t::get_residue_circles_top_left() const {

   bool status = false;
   lig_build::pos_t p(999999,999999);

   if (residue_circles.size()) {
      status = true;
      for (unsigned int i=0; i<residue_circles.size(); i++) {
	 if (residue_circles[i].pos.x < p.x) p.x = residue_circles[i].pos.x;
	 if (residue_circles[i].pos.y < p.y) p.y = residue_circles[i].pos.y;
      }
   }
   return std::pair<bool, lig_build::pos_t> (status, p);
}




// Return the status and updated positions.
//
std::pair<int, std::vector<residue_circle_t> >
pli::optimise_residue_circles::solution() const {
   return std::pair<int, std::vector<residue_circle_t> > (status, current_circles);
}


void
pli::optimise_residue_circles::setup_angles() {

   for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
      std::vector<int> residue_indexes;
      for (unsigned int ires=0; ires<current_circles.size(); ires++) {
	 for (unsigned int ibond=0;
	      ibond<current_circles[ires].bonds_to_ligand.size();
	      ibond++) {
	    if (current_circles[ires].bonds_to_ligand[ibond].ligand_atom_name == mol.atoms[iat].get_atom_name()) {
	       residue_indexes.push_back(ires);
	    }
	 }
      }
      if (residue_indexes.size() > 1) {
	 // a bit complicated.
	 angle a1(iat, residue_indexes[0], residue_indexes[1]);
	 angles.push_back(a1);
	 if (residue_indexes.size() > 2) {
	    angle a2(iat, residue_indexes[1], residue_indexes[2]);
	    angles.push_back(a2);
	    angle a3(iat, residue_indexes[0], residue_indexes[2]);
	    angles.push_back(a3);
	 }
      }
   }

}


// static
double
pli::optimise_residue_circles::f(const gsl_vector *v, void *params) {

   double score = 0.0;

   pli::optimise_residue_circles *orc = static_cast<optimise_residue_circles *>(params);

   double rk = 3000.0;
   double exp_scale = 0.0021; // 0.002 is too much kickage

   for (unsigned int i=0; i<orc->current_circles.size(); i++) {

      if (score_vs_ligand_atoms) {
	 // first score against the fixed (ligand) atoms and ring centres
	 // (Do I need to test for not being a hydrogen?)
	 ///
	 for (unsigned int iat=0; iat<orc->mol.atoms.size(); iat++) {
	    double d_pt_1 = gsl_vector_get(v, 2*i  ) - orc->mol.atoms[iat].atom_position.x;
	    double d_pt_2 = gsl_vector_get(v, 2*i+1) - orc->mol.atoms[iat].atom_position.y;
	    double d2 = d_pt_1 * d_pt_1 + d_pt_2 * d_pt_2;
	    score += rk * exp(-0.5*exp_scale*d2);
	 }
      }

      if (score_vs_ring_centres) {
	 // score against the ligand ring centres
	 //

	 std::vector<lig_build::pos_t> mol_ring_centres = orc->mol.get_ring_centres();

	 for (unsigned int irc=0; irc<mol_ring_centres.size(); irc++) {
	    double d_pt_1 = gsl_vector_get(v, 2*i  ) - mol_ring_centres[irc].x;
	    double d_pt_2 = gsl_vector_get(v, 2*i+1) - mol_ring_centres[irc].y;
	    double d2 = d_pt_1 * d_pt_1 + d_pt_2 * d_pt_2;
	    score += rk * exp(-0.5*exp_scale*d2);
	 }
      }


      if (score_vs_other_residues) {
	 // score against the other residue centres
	 //
	 double kk = 300.0;
	 for (unsigned int ic=0; ic<orc->current_circles.size(); ic++) {
	    if (ic != i) {
	       double d_pt_1 = gsl_vector_get(v, 2*i  ) - gsl_vector_get(v, 2*ic  );
	       double d_pt_2 = gsl_vector_get(v, 2*i+1) - gsl_vector_get(v, 2*ic+1);
	       double d2 = d_pt_1 * d_pt_1 + d_pt_2 * d_pt_2;
	       score += kk * exp(-0.5*exp_scale*d2);
	    }
	 }
      }

      if (score_vs_original_positions) {
	 // score against the original position (d^2 function)
	 //
	 double k = 0.002;
	 double d_1 = gsl_vector_get(v, 2*i  ) - orc->starting_circles[i].pos.x;
	 double d_2 = gsl_vector_get(v, 2*i+1) - orc->starting_circles[i].pos.y;
	 score += k * (d_1*d_1 + d_2*d_2);
      }
   }

   if (score_vs_other_residues_for_angles) {
      double k_angle_scale = 1.0;
      for (unsigned int iang=0; iang<orc->angles.size(); iang++) {
	 // we want to get cos_theta and shove it in the penalty
	 // score function. cos_theta comes from the dot product.

	 lig_build::pos_t &at_pos = orc->mol.atoms[orc->angles[iang].i_atom_index].atom_position;
	 int idx_1 = orc->angles[iang].ires_1_index;
	 int idx_2 = orc->angles[iang].ires_2_index;
	 lig_build::pos_t r1_pos(gsl_vector_get(v, 2*idx_1), gsl_vector_get(v, 2*idx_1+1));
	 lig_build::pos_t r2_pos(gsl_vector_get(v, 2*idx_2), gsl_vector_get(v, 2*idx_2+1));
	 lig_build::pos_t A = r1_pos - at_pos;
	 lig_build::pos_t B = r2_pos - at_pos;

	 double a_dot_b = lig_build::pos_t::dot(A,B);
	 double cos_theta = a_dot_b/(A.length()*B.length());
	 double one_minus_ct = 1.0-cos_theta;
	 double angle_penalty = k_angle_scale * exp(-2.5*one_minus_ct*one_minus_ct);
      }
   }

   // score vs attachment length, only for primary residues. What
   // is the index list of primary residues?
   //
   // If we do this, we can scrap score_vs_original_positions
   // perhaps.  Or reduce its weight.
   //
   // A simple quadratic function.
   //

   if (score_vs_ligand_atom_bond_length) {
      double kk_bond_length_scale = 10;
      for (unsigned int iprimary=0; iprimary<orc->primary_indices.size(); iprimary++) {
	 int idx = orc->primary_indices[iprimary];
	 std::vector<std::pair<lig_build::pos_t, double> > attachment_points =
	    orc->current_circles[idx].get_attachment_points(orc->mol);
	 for (unsigned int iattach=0; iattach<attachment_points.size(); iattach++) {
	    //
	    double target_length = 1.5 * LIGAND_TO_CANVAS_SCALE_FACTOR *
	       attachment_points[iattach].second;
	    lig_build::pos_t current_pos(gsl_vector_get(v, 2*idx),
					 gsl_vector_get(v, 2*idx+1));
	    lig_build::pos_t bond_vector = (attachment_points[iattach].first - current_pos);
	    double dist_to_attachment_point = bond_vector.length();
	    double d = dist_to_attachment_point - target_length;
	    double bond_length_penalty = kk_bond_length_scale * d * d;
            std::cout << "bond_vector: " << iprimary << " " << iattach << " "
                      << bond_vector.x << " " << bond_vector.y << " d  " << d << std::endl;
	    score += bond_length_penalty;
	 }
      }
   }
   std::cout << "optimise_residue_circles() returning score " << score << std::endl;
   return score;
}

