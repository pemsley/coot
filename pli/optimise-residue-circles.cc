
#include "flev-scale-factor.hh"
#include "optimise-residue-circles.hh"

// I don't know why these are static.
bool pli::optimise_residue_circles::score_vs_ligand_atoms              = true;
bool pli::optimise_residue_circles::score_vs_ring_centres              = false;
bool pli::optimise_residue_circles::score_vs_other_residues            = true;
bool pli::optimise_residue_circles::score_vs_original_positions        = true;
bool pli::optimise_residue_circles::score_vs_other_residues_for_angles = false; //  problem. Fix later
bool pli::optimise_residue_circles::score_vs_ligand_atom_bond_length   = true;

pli::optimise_residue_circles::optimise_residue_circles(const std::vector<residue_circle_t> &r,
                                                        const std::vector<residue_circle_t> &c,
                                                        const svg_molecule_t &mol_in,
                                                        const std::vector<int> &primary_indices_in) {

   bool show_dynamics = false;

   score_vs_ligand_atoms_rk = 3000.0;
   score_vs_ligand_atoms_exp_scale = 0.0008; // quite sensitive value - 0.0002 is big kick
   score_vs_original_positions_kk = 1.1;
   score_vs_ligand_atom_bond_length_kk = 0.001;

   score_vs_other_residues_kk = 0.5;
   score_vs_other_residues_exp_scale = 0.001;

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
      std::cout << "starting_circle " << i << " is at " << current_circles[i].pos.x << " " << current_circles[i].pos.y
                << std::endl;
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
      std::cout << " NO_PROGRESS\n";

   for (unsigned int i=0; i<current_circles.size(); i++) {
      current_circles[i].pos.x = gsl_vector_get(s->x, 2*i  );
      current_circles[i].pos.y = gsl_vector_get(s->x, 2*i+1);
   }

   gsl_multimin_fdfminimizer_free (s);
   gsl_vector_free (x);

}

// static
double
pli::optimise_residue_circles::f(const gsl_vector *v, void *params) {

   double score = 0.0;

   pli::optimise_residue_circles *orc = static_cast<optimise_residue_circles *>(params);

   for (unsigned int i=0; i<orc->current_circles.size(); i++) {

      if (orc->score_vs_ligand_atoms) {

         double rk = orc->score_vs_ligand_atoms_rk;
         double exp_scale = orc->score_vs_ligand_atoms_exp_scale; // 0.002 is too much kickage

	 // first score against the fixed (ligand) atoms and ring centres
	 // (Do I need to test for not being a hydrogen?)
	 ///
	 for (unsigned int iat=0; iat<orc->mol.atoms.size(); iat++) {
	    const double &d_pt_1 = gsl_vector_get(v, 2*i  ) - orc->mol.atoms[iat].atom_position.x;
	    const double &d_pt_2 = gsl_vector_get(v, 2*i+1) - orc->mol.atoms[iat].atom_position.y;
	    double d2 = d_pt_1 * d_pt_1 + d_pt_2 * d_pt_2;
	    score += rk * exp(-0.5*exp_scale*d2);
            std::cout << "circles-vs-ligand-atoms circle " << i << " iat: " << iat << " at "
                      << orc->mol.atoms[iat].atom_position.x << " "
                      << orc->mol.atoms[iat].atom_position.y << " circle: "
                      << d_pt_1 << " " << d_pt_2 << " d " << std::sqrt(d2) << " score " << score << std::endl;
	 }
      }

      if (orc->score_vs_ring_centres) {
	 // score against the ligand ring centres
	 //

	 std::vector<lig_build::pos_t> mol_ring_centres = orc->mol.get_ring_centres();

         double rk = orc->score_vs_ligand_atoms_rk; // uses same params as for atoms
         double exp_scale = orc->score_vs_ligand_atoms_exp_scale; // 0.002 is too much kickage

	 for (unsigned int irc=0; irc<mol_ring_centres.size(); irc++) {
	    double d_pt_1 = gsl_vector_get(v, 2*i  ) - mol_ring_centres[irc].x;
	    double d_pt_2 = gsl_vector_get(v, 2*i+1) - mol_ring_centres[irc].y;
	    double d2 = d_pt_1 * d_pt_1 + d_pt_2 * d_pt_2;
	    score += rk * exp(-0.5*exp_scale*d2);
	 }
      }


      if (orc->score_vs_other_residues) {
	 // score against the other residue centres
	 //
	 double kk        = orc->score_vs_other_residues_kk;
         double exp_scale = orc->score_vs_other_residues_exp_scale;
	 for (unsigned int ic=0; ic<orc->current_circles.size(); ic++) {
	    if (ic != i) {
	       double d_pt_1 = gsl_vector_get(v, 2*i  ) - gsl_vector_get(v, 2*ic  );
	       double d_pt_2 = gsl_vector_get(v, 2*i+1) - gsl_vector_get(v, 2*ic+1);
	       double d2 = d_pt_1 * d_pt_1 + d_pt_2 * d_pt_2;
	       score += kk * exp(-0.5*exp_scale*d2);
	    }
	 }
      }

      if (orc->score_vs_original_positions) {
	 // score against the original position (d^2 function)
	 //
	 double k = orc->score_vs_original_positions_kk;
	 double d_1 = gsl_vector_get(v, 2*i  ) - orc->starting_circles[i].pos.x;
	 double d_2 = gsl_vector_get(v, 2*i+1) - orc->starting_circles[i].pos.y;
	 score += k * (d_1*d_1 + d_2*d_2);
      }
   }

   if (orc->score_vs_other_residues_for_angles) {
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

   if (orc->score_vs_ligand_atom_bond_length) {
      double kk = orc->score_vs_ligand_atom_bond_length_kk;
      for (unsigned int iprimary=0; iprimary<orc->primary_indices.size(); iprimary++) {
	 int idx = orc->primary_indices[iprimary];
	 std::vector<std::pair<lig_build::pos_t, double> > attachment_points =
	    orc->current_circles[idx].get_attachment_points(orc->mol);
	 for (unsigned int iattach=0; iattach<attachment_points.size(); iattach++) {
	    //
	    double target_length = 1.5 * LIGAND_TO_CANVAS_SCALE_FACTOR *
	       attachment_points[iattach].second * 3.0;
	    lig_build::pos_t current_pos(gsl_vector_get(v, 2*idx),
					 gsl_vector_get(v, 2*idx+1));
	    lig_build::pos_t bond_vector = (attachment_points[iattach].first - current_pos);
	    double dist_to_attachment_point = bond_vector.length();
	    double d = dist_to_attachment_point - target_length;
	    double bond_length_penalty = kk * d * d;
            std::cout << "bond_vector: attachment " << iprimary << " " << iattach << " "
                      << bond_vector.x << " " << bond_vector.y << " d  " << d << std::endl;
	    score += bond_length_penalty;
	 }
      }
   }
   std::cout << "optimise_residue_circles() returning score " << score << std::endl;
   return score;
}

// should this be a static member function?
void
pli::optimise_residue_circles::df(const gsl_vector *v, void *params, gsl_vector *df) {

   optimise_residue_circles *orc = static_cast<optimise_residue_circles *> (params);

   // initialise the derivatives vector
   for (unsigned int i=0; i<orc->current_circles.size(); i++) {
      gsl_vector_set(df, 2*i,   0.0);
      gsl_vector_set(df, 2*i+1, 0.0);
   }

   for (unsigned int i=0; i<orc->current_circles.size(); i++) {

      double rk        = orc->score_vs_ligand_atoms_rk;
      double exp_scale = orc->score_vs_ligand_atoms_exp_scale;

      if (orc->score_vs_ligand_atoms) {
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

      if (orc->score_vs_ring_centres) {

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

      if (orc->score_vs_other_residues) {
         // score against the other residue centres
         //
	 double kk        = orc->score_vs_other_residues_kk;
         double exp_scale = orc->score_vs_other_residues_exp_scale;
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

      if (orc->score_vs_original_positions) {
          // score against the original position (1/d^2 function)
          //
          double k = orc->score_vs_original_positions_kk;
          double df_part_x = 2.0 * (gsl_vector_get(v, 2*i  ) - orc->starting_circles[i].pos.x);
          double df_part_y = 2.0 * (gsl_vector_get(v, 2*i+1) - orc->starting_circles[i].pos.y);
          // std::cout << "  gradient " << i << "  " << df_part_x << "  " << df_part_y << std::endl;

          gsl_vector_set(df, 2*i,   gsl_vector_get(df, 2*i  ) + k * df_part_x);
          gsl_vector_set(df, 2*i+1, gsl_vector_get(df, 2*i+1) + k * df_part_y);
      }
   }


   if (orc->score_vs_other_residues_for_angles) {
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

   if (orc->score_vs_ligand_atom_bond_length) {
      double kk = orc->score_vs_ligand_atom_bond_length_kk;
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
               attachment_points[iattach].second * 3.0;
            lig_build::pos_t current_pos(gsl_vector_get(v, 2*idx),
                                         gsl_vector_get(v, 2*idx+1));
            lig_build::pos_t bond_vector = (attachment_points[iattach].first - current_pos);
            double dist_to_attachment_point = bond_vector.length();
            double frac_bond_length_dev = (dist_to_attachment_point - target_length)/dist_to_attachment_point;

            // debuging
            double debug_old_v1 =   gsl_vector_get(df, 2*idx);
            double debug_old_v2 =   gsl_vector_get(df, 2*idx+1);

            gsl_vector_set(df, 2*idx,   gsl_vector_get(df, 2*idx)
                           + 2.0 * kk * frac_bond_length_dev
                           * (current_pos.x - attachment_points[iattach].first.x));
            gsl_vector_set(df, 2*idx+1, gsl_vector_get(df, 2*idx+1)
                           + 2.0 * kk * frac_bond_length_dev
                           * (current_pos.y - attachment_points[iattach].first.y));

            if (false) {
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
pli::optimise_residue_circles::fdf(const gsl_vector *x, void *params, double *f_in, gsl_vector *df_in) {

   *f_in = f(x, params);
   df(x, params, df_in);
}

void
pli::optimise_residue_circles::numerical_gradients(gsl_vector *x, gsl_vector *df, void *params) const {

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


std::vector<std::pair<lig_build::pos_t, double> >
residue_circle_t::get_attachment_points(const svg_molecule_t &mol) const {

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

   if (get_stacking_type() == residue_circle_t::CATION_PI_STACKING) {
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

