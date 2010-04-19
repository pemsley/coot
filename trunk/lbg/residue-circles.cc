
#include "lbg.hh"

bool lbg_info_t::optimise_residue_circles::score_vs_ligand_atoms = 1;
bool lbg_info_t::optimise_residue_circles::score_vs_ring_centres = 1;
bool lbg_info_t::optimise_residue_circles::score_vs_other_residues = 1;
bool lbg_info_t::optimise_residue_circles::score_vs_original_positions = 1;

lbg_info_t::optimise_residue_circles::optimise_residue_circles(const std::vector<lbg_info_t::residue_circle_t> &r,
							       const std::vector<lbg_info_t::residue_circle_t> &c,
							       const widgeted_molecule_t &mol_in) {
   mol = mol_in;
   starting_circles = r;
   current_circles = c;

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

   bool score_vs_ligand_atoms = 1;
   bool score_vs_ring_centres = 1;
   bool score_vs_other_residues = 1;
   bool score_vs_original_positions = 1;
   
   do {
	 iter++;
	 status = gsl_multimin_fdfminimizer_iterate (s);
	 
	 if (status)
	    break;
	 
	 status = gsl_multimin_test_gradient (s->gradient, 1e-3);
	 
	 if (status == GSL_SUCCESS)
	    printf ("Minimum found at:\n");
	 
// 	 printf ("iteration: %d %.5f %.5f %10.5f\n", int(iter),
// 		 gsl_vector_get (s->x, 0), 
// 		 gsl_vector_get (s->x, 1), 
// 		 s->f);

	 
   } while (status == GSL_CONTINUE && iter < 60);

   std::cout << "============= terminating at iter " << iter << " with status "
	     << status;
   if (status == GSL_ENOPROG)
      std::cout << " NO_PROGRESS";
   std::cout << " ============ " << std::endl;

   for (unsigned int i=0; i<current_circles.size(); i++) {
      if (i==0) { 
	 std::cout << "replacing " << current_circles[i].pos.x << " with " << gsl_vector_get(s->x, 2*i  ) << "    ";
	 std::cout << "replacing " << current_circles[i].pos.y << " with " << gsl_vector_get(s->x, 2*i+1) << std::endl;
      }
      current_circles[i].pos.x = gsl_vector_get(s->x, 2*i  );
      current_circles[i].pos.y = gsl_vector_get(s->x, 2*i+1);
   }
   
   gsl_multimin_fdfminimizer_free (s);
   gsl_vector_free (x);
   
}

// Return the status and updated positions.
// 
std::pair<int, std::vector<lbg_info_t::residue_circle_t> > 
lbg_info_t::optimise_residue_circles::solution() const {
   return std::pair<int, std::vector<lbg_info_t::residue_circle_t> > (status, current_circles);
}

double
lbg_info_t::optimise_residue_circles::f(const gsl_vector *v, void *params) {

   double score = 0.0;

   double atom_scale = 1.0;

   lbg_info_t::optimise_residue_circles *orc = static_cast<optimise_residue_circles *>(params);

   double rk = 3000.0;
   double exp_scale = 0.0025; // 0.002 is too much kickage
   
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
	 
	 for (unsigned int ib=0; ib<orc->mol.bonds.size(); ib++) { 
	    if (orc->mol.bonds[ib].have_centre_pos()) {
	       double d_pt_1 = gsl_vector_get(v, 2*i  ) - orc->mol.bonds[ib].centre_pos().x;
	       double d_pt_2 = gsl_vector_get(v, 2*i+1) - orc->mol.bonds[ib].centre_pos().y;
	       double d2 = d_pt_1 * d_pt_1 + d_pt_2 * d_pt_2;
	       score += rk * exp(-0.5*exp_scale*d2);
	    }
	 }
      }


      if (score_vs_other_residues) { 
	 // score against the other residue centres
	 //
	 double kk = 30.0;
	 for (unsigned int ic=0; ic<orc->current_circles.size(); ic++) {
	    if (ic != i) {
	       double d_pt_1 = gsl_vector_get(v, 2*i  ) - gsl_vector_get(v, 2*ic  );
	       double d_pt_2 = gsl_vector_get(v, 2*i+1) - gsl_vector_get(v, 2*ic+1);
	       double d2 = d_pt_1 * d_pt_1 + d_pt_2 * d_pt_2;
	       // if (d2 < 0.001) d2 = 0.001;
	       score += kk * 1.0/sqrt(d2);
	    }
	 }
      }

      if (score_vs_original_positions) { 
	 // score against the original position (1/d^2 function)
	 // 
	 double k = 0.0002;
	 double d_1 = gsl_vector_get(v, 2*i  ) - orc->starting_circles[i].pos.x;
	 double d_2 = gsl_vector_get(v, 2*i+1) - orc->starting_circles[i].pos.y;
	 // std::cout << "score addition " << i << "  " << d_1*d_1 << " + " << d_2*d_2 << std::endl;
	 score += k * (d_1*d_1 + d_2*d_2);
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
      double exp_scale = 0.0025;
      
      if (score_vs_ligand_atoms) { 
	 // first score against the fixed (ligand) atoms...
	 // 
	 for (unsigned int iat=0; iat<orc->mol.atoms.size(); iat++) {
	    double d_pt_1 = gsl_vector_get(v, 2*i  ) - orc->mol.atoms[iat].atom_position.x;
	    double d_pt_2 = gsl_vector_get(v, 2*i+1) - orc->mol.atoms[iat].atom_position.y;
	    double d2 = d_pt_1 * d_pt_1 + d_pt_2 * d_pt_2;
	    
// 	    double d = sqrt(d2);
// 	    double df_part_x = (-1.0/d2)*(d_pt_1/d);
// 	    double df_part_y = (-1.0/d2)*(d_pt_2/d);

	    double df_part_x =  rk * d_pt_1 * -exp_scale  * exp(-0.5*exp_scale*d2);
	    double df_part_y =  rk * d_pt_2 * -exp_scale  * exp(-0.5*exp_scale*d2);
	    
	    gsl_vector_set(df, 2*i,   gsl_vector_get(df, 2*i  ) + df_part_x);
	    gsl_vector_set(df, 2*i+1, gsl_vector_get(df, 2*i+1) + df_part_y);
	 }
      }

      if (score_vs_ring_centres) { 

	 //  score against the ligand ring centres
	 //
	 for (unsigned int ib=0; ib<orc->mol.bonds.size(); ib++) { 
	    if (orc->mol.bonds[ib].have_centre_pos()) {
	       double d_pt_1 = gsl_vector_get(v, 2*i  ) - orc->mol.bonds[ib].centre_pos().x;
	       double d_pt_2 = gsl_vector_get(v, 2*i+1) - orc->mol.bonds[ib].centre_pos().y;
	       double d2 = d_pt_1 * d_pt_1 + d_pt_2 * d_pt_2;
	       // double d = sqrt(d2);
	       // double df_part_x = (-1.0/d2)*(d_pt_1/d);
	       // double df_part_y = (-1.0/d2)*(d_pt_2/d);

	       double df_part_x =  rk * d_pt_1 * -exp_scale  * exp(-0.5*exp_scale*d2);
	       double df_part_y =  rk * d_pt_2 * -exp_scale  * exp(-0.5*exp_scale*d2);
	       
	       gsl_vector_set(df, 2*i,   gsl_vector_get(df, 2*i  ) + df_part_x);
	       gsl_vector_set(df, 2*i+1, gsl_vector_get(df, 2*i+1) + df_part_y);
	    }
	 }
      }

      if (score_vs_other_residues) { 
	 // score against the other residue centres
	 //
	 double kk = 10.0;
	 for (unsigned int ic=0; ic<orc->current_circles.size(); ic++) {
	    if (ic != i) {
	       double d_pt_1 = gsl_vector_get(v, 2*i  ) - gsl_vector_get(v, 2*ic  );
	       double d_pt_2 = gsl_vector_get(v, 2*i+1) - gsl_vector_get(v, 2*ic+1);
	       double d2 = d_pt_1 * d_pt_1 + d_pt_2 * d_pt_2;
	       // if (d2 < 0.001) d2 = 0.001;
	       double d = sqrt(d2);
	       double df_part_x = -2*d_pt_1/(d*d*d);
	       double df_part_y = -2*d_pt_2/(d*d*d);
	       gsl_vector_set(df, 2*i,   gsl_vector_get(df, 2*i  ) + kk * df_part_x);
	       gsl_vector_set(df, 2*i+1, gsl_vector_get(df, 2*i+1) + kk * df_part_y);	    
	    }
	 }
      }

      if (score_vs_original_positions) { 
	  // score against the original position (1/d^2 function)
	  //
	  double k = 0.0002;
	  double df_part_x = 2.0 * (gsl_vector_get(v, 2*i  ) - orc->starting_circles[i].pos.x);
	  double df_part_y = 2.0 * (gsl_vector_get(v, 2*i+1) - orc->starting_circles[i].pos.y);
	  // std::cout << "  gradient " << i << "  " << df_part_x << "  " << df_part_y << std::endl;

	  gsl_vector_set(df, 2*i,   gsl_vector_get(df, 2*i  ) + k * df_part_x);
	  gsl_vector_set(df, 2*i+1, gsl_vector_get(df, 2*i+1) + k * df_part_y);	    
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
      std::cout << "gradient_comparison " << gsl_vector_get(df, i) << "    " << v_av/micro_step << std::endl; 
   }
}
