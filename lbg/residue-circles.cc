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


#include <Python.h>

#include "lbg.hh"
bool lbg_info_t::optimise_residue_circles::score_vs_ligand_atoms       = 1;
bool lbg_info_t::optimise_residue_circles::score_vs_ring_centres       = 1;
bool lbg_info_t::optimise_residue_circles::score_vs_other_residues     = 1;
bool lbg_info_t::optimise_residue_circles::score_vs_original_positions = 1;
bool lbg_info_t::optimise_residue_circles::score_vs_other_residues_for_angles = 0;
bool lbg_info_t::optimise_residue_circles::score_vs_ligand_atom_bond_length = 1;


lbg_info_t::optimise_residue_circles::optimise_residue_circles(const std::vector<lbg_info_t::residue_circle_t> &r,
							       const std::vector<lbg_info_t::residue_circle_t> &c,
							       const widgeted_molecule_t &mol_in,
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
   my_func.f   = &lbg_info_t::optimise_residue_circles::f;
   my_func.df  = &lbg_info_t::optimise_residue_circles::df;
   my_func.fdf = &lbg_info_t::optimise_residue_circles::fdf;
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

   if (0) 
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


void
lbg_info_t::refine_residue_circle_positions() { // changes the positions of residue_circles

   std::vector<int> primary_indices = get_primary_indices();
   initial_residues_circles_layout(); // twiddle residue_circles
   std::vector<lbg_info_t::residue_circle_t> current_circles = residue_circles;
   
   for (int iround=0; iround<30; iround++) {
      std::pair<int, std::vector<lbg_info_t::residue_circle_t> > new_c =
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
lbg_info_t::get_residue_circles_top_left() const {

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
std::pair<int, std::vector<lbg_info_t::residue_circle_t> > 
lbg_info_t::optimise_residue_circles::solution() const {
   return std::pair<int, std::vector<lbg_info_t::residue_circle_t> > (status, current_circles);
}


void
lbg_info_t::optimise_residue_circles::setup_angles() {

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
lbg_info_t::optimise_residue_circles::f(const gsl_vector *v, void *params) {

   double score = 0.0;

   double atom_scale = 1.0;

   lbg_info_t::optimise_residue_circles *orc = static_cast<optimise_residue_circles *>(params);

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
	    score += bond_length_penalty;
	 }
      }
   }
   return score;
}


void
lbg_info_t::optimise_residue_circles::df(const gsl_vector *v, void *params, gsl_vector *df) {

   double atom_scale = 1.0;

   lbg_info_t::optimise_residue_circles *orc = static_cast<optimise_residue_circles *> (params);

   // initialise the derivatives vector
   for (unsigned int i=0; i<orc->current_circles.size(); i++) { 
      gsl_vector_set(df, 2*i,   0.0);
      gsl_vector_set(df, 2*i+1, 0.0);
   }

   for (unsigned int i=0; i<orc->current_circles.size(); i++) { 

      double rk = 3000.0;
      double exp_scale = 0.0021;
      
      if (score_vs_ligand_atoms) { 
	 // first score against the fixed (ligand) atoms...
	 // 
	 for (unsigned int iat=0; iat<orc->mol.atoms.size(); iat++) {
	    double d_pt_1 = gsl_vector_get(v, 2*i  ) - orc->mol.atoms[iat].atom_position.x;
	    double d_pt_2 = gsl_vector_get(v, 2*i+1) - orc->mol.atoms[iat].atom_position.y;
	    double d2 = d_pt_1 * d_pt_1 + d_pt_2 * d_pt_2;
	    
	    double df_part_x =  rk * d_pt_1 * -exp_scale  * exp(-0.5*exp_scale*d2);
	    double df_part_y =  rk * d_pt_2 * -exp_scale  * exp(-0.5*exp_scale*d2);
	    
	    gsl_vector_set(df, 2*i,   gsl_vector_get(df, 2*i  ) + df_part_x);
	    gsl_vector_set(df, 2*i+1, gsl_vector_get(df, 2*i+1) + df_part_y);
	 }
      }

      if (score_vs_ring_centres) { 

	 //  score against the ligand ring centres
	 //
	 std::vector<lig_build::pos_t> mol_ring_centres = orc->mol.get_ring_centres();

	 for (unsigned int irc=0; irc<mol_ring_centres.size(); irc++) { 
	    double d_pt_1 = gsl_vector_get(v, 2*i  ) - mol_ring_centres[irc].x;
	    double d_pt_2 = gsl_vector_get(v, 2*i+1) - mol_ring_centres[irc].y;
	    double d2 = d_pt_1 * d_pt_1 + d_pt_2 * d_pt_2;
	    
	    double df_part_x =  rk * d_pt_1 * -exp_scale  * exp(-0.5*exp_scale*d2);
	    double df_part_y =  rk * d_pt_2 * -exp_scale  * exp(-0.5*exp_scale*d2);
	    
	    gsl_vector_set(df, 2*i,   gsl_vector_get(df, 2*i  ) + df_part_x);
	    gsl_vector_set(df, 2*i+1, gsl_vector_get(df, 2*i+1) + df_part_y);
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

	       double df_part_x =  2 * kk * d_pt_1 * -exp_scale  * exp(-0.5*exp_scale*d2);
	       double df_part_y =  2 * kk * d_pt_2 * -exp_scale  * exp(-0.5*exp_scale*d2);
	       gsl_vector_set(df, 2*i,   gsl_vector_get(df, 2*i  ) + df_part_x);
	       gsl_vector_set(df, 2*i+1, gsl_vector_get(df, 2*i+1) + df_part_y);	    
	    }
	 }
      }

      if (score_vs_original_positions) { 
	  // score against the original position (1/d^2 function)
	  //
	  double k = 0.002;
	  double df_part_x = 2.0 * (gsl_vector_get(v, 2*i  ) - orc->starting_circles[i].pos.x);
	  double df_part_y = 2.0 * (gsl_vector_get(v, 2*i+1) - orc->starting_circles[i].pos.y);
	  // std::cout << "  gradient " << i << "  " << df_part_x << "  " << df_part_y << std::endl;

	  gsl_vector_set(df, 2*i,   gsl_vector_get(df, 2*i  ) + k * df_part_x);
	  gsl_vector_set(df, 2*i+1, gsl_vector_get(df, 2*i+1) + k * df_part_y);	    
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
	 double A_length = A.length();
	 double B_length = B.length();
	 double A_length_recip = 1.0/A_length;
	 double B_length_recip = 1.0/B_length;
	    
	 double a_dot_b = lig_build::pos_t::dot(A,B);
	 double cos_theta = a_dot_b/(A.length()*B.length());

	 double gamma = cos_theta;
	 double one_minus_ct = 1.0-cos_theta;
	 double angle_penalty = k_angle_scale * exp(-2.5*one_minus_ct*one_minus_ct);

	 // note: S1 contains k_angle_scale
	 double S1 = -2.5*angle_penalty*(-2*(1-gamma));
	 double S2_x1 = (r2_pos.x - at_pos.x) * A_length_recip * B_length_recip + (r1_pos.x - at_pos.x) * (-1.0 * A_length_recip * A_length_recip * A_length_recip * B_length_recip * a_dot_b);
	 double S2_x2 = (r1_pos.x - at_pos.x) * B_length_recip * A_length_recip + (r1_pos.x - at_pos.x) * (-1.0 * B_length_recip * B_length_recip * B_length_recip * A_length_recip * a_dot_b);
	 double S2_y1 = (r2_pos.y - at_pos.y) * A_length_recip * B_length_recip + (r1_pos.y - at_pos.y) * (-1.0 * A_length_recip * A_length_recip * A_length_recip * B_length_recip * a_dot_b);
	 double S2_y2 = (r1_pos.y - at_pos.y) * B_length_recip * A_length_recip + (r1_pos.y - at_pos.y) * (-1.0 * B_length_recip * B_length_recip * B_length_recip * A_length_recip * a_dot_b);

	 gsl_vector_set(df, 2*idx_1,   gsl_vector_get(df, 2*idx_1)   * S1 * S2_x1);
	 gsl_vector_set(df, 2*idx_2,   gsl_vector_get(df, 2*idx_2)   * S1 * S2_x2);
	 gsl_vector_set(df, 2*idx_1+1, gsl_vector_get(df, 2*idx_1+1) * S1 * S2_y1);
	 gsl_vector_set(df, 2*idx_2+1, gsl_vector_get(df, 2*idx_2+1) * S1 * S2_y2);
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
	    // 1.5 is a fudge factor to make sure that the residue
	    // circle is aesthetically distanced from the ligand
	    // atom.
	    // 
	    double target_length = 1.5 * LIGAND_TO_CANVAS_SCALE_FACTOR *
	       attachment_points[iattach].second;
	    lig_build::pos_t current_pos(gsl_vector_get(v, 2*idx),
					 gsl_vector_get(v, 2*idx+1));
	    lig_build::pos_t bond_vector = (attachment_points[iattach].first - current_pos);
	    double dist_to_attachment_point = bond_vector.length();
	    double frac_bond_length_dev = (dist_to_attachment_point - target_length)/dist_to_attachment_point;

	    // debuging
	    double debug_old_v1 =   gsl_vector_get(df, 2*idx);
	    double debug_old_v2 =   gsl_vector_get(df, 2*idx+1);
	       
	    gsl_vector_set(df, 2*idx,   gsl_vector_get(df, 2*idx)
			   + 2.0 * kk_bond_length_scale * frac_bond_length_dev
			   * (current_pos.x - attachment_points[iattach].first.x));
	    gsl_vector_set(df, 2*idx+1, gsl_vector_get(df, 2*idx+1)
			   + 2.0 * kk_bond_length_scale * frac_bond_length_dev
			   * (current_pos.y - attachment_points[iattach].first.y));

	    if (0) { 
	       std::cout << "some numbers: " << current_pos << " " << bond_vector << " "
			 << dist_to_attachment_point << " " << target_length << " "
			 << frac_bond_length_dev << " "
			 << (current_pos.x - attachment_points[iattach].first.x) << " "
			 << (current_pos.y - attachment_points[iattach].first.y) << "\n";
	       double debug_new_v1 =   gsl_vector_get(df, 2*idx);
	       double debug_new_v2 =   gsl_vector_get(df, 2*idx+1);
	       std::cout << "reset df for " << 2*idx    << " from " << debug_old_v1
			 << " to " << debug_new_v1 << "\n";
	       std::cout << "reset df for " << 2*idx+1  << " from " << debug_old_v2
			 << " to " << debug_new_v2 << "\n";
	    }
	 }
      }
   }

   // orc->numerical_gradients((gsl_vector *)v, df, params);
   
}

void
lbg_info_t::optimise_residue_circles::fdf(const gsl_vector *x, void *params, double *f_in, gsl_vector *df_in) {

   *f_in = f(x, params);
   df(x, params, df_in);
}

void
lbg_info_t::optimise_residue_circles::numerical_gradients(gsl_vector *x, gsl_vector *df, void *params) const {

   
   double micro_step = 0.0001;  // the difference between the gradients
			        // seems not to depend on the
			        // micro_step size (0.0001 vs 0.001)
   
//    std::cout << "analytical_gradients" << std::endl; 
//    for (unsigned int i=0; i<df->size; i++) {
//       double tmp = gsl_vector_get(df, i); 
//       std::cout << i << "  " << tmp << std::endl; 
//    }

   for (unsigned int i=0; i<x->size; i++) { 
      
      double tmp = gsl_vector_get(x, i);
      gsl_vector_set(x, i, tmp+micro_step); 
      double v1 = f(x, params);
      gsl_vector_set(x, i, tmp-micro_step); 
      double v2 = f(x, params);
      gsl_vector_set(x, i, tmp);
      double v_av = 0.5 * (v1 - v2);
      std::cout << "gradient_comparison " << i << " " << gsl_vector_get(df, i) << "    " << v_av/micro_step << std::endl; 
   }
}
