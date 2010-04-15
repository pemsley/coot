
#include "lbg.hh"


lbg_info_t::optimise_residue_circles::optimise_residue_circles(const std::vector<lbg_info_t::residue_circle_t> &r,
							       const widgeted_molecule_t &mol_in) {
   mol = mol_in;
   input_circles = r;

   gsl_multimin_function_fdf my_func;

   int n_var = input_circles.size()*2;
   my_func.n = n_var;  /* number of function components */
   my_func.f   = &lbg_info_t::optimise_residue_circles::f;
   my_func.df  = &lbg_info_t::optimise_residue_circles::df;
   my_func.fdf = &lbg_info_t::optimise_residue_circles::fdf;
   my_func.params = static_cast<void *> (this);
   
   // set the starting points
   gsl_vector *x = gsl_vector_alloc(n_var);
   
   for (unsigned int i=0; i<input_circles.size(); i++) {
      gsl_vector_set(x, 2*i,   input_circles[i].pos.x);
      gsl_vector_set(x, 2*i+1, input_circles[i].pos.y);
   }

   gsl_multimin_fdfminimizer *s;
   const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_fr;
   
   s = gsl_multimin_fdfminimizer_alloc (T, n_var);
   gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-4);
   size_t iter = 0;
   
   do {
	 iter++;
	 status = gsl_multimin_fdfminimizer_iterate (s);
	 
	 if (status)
	    break;
	 
	 status = gsl_multimin_test_gradient (s->gradient, 1e-3);
	 
	 if (status == GSL_SUCCESS)
	    printf ("Minimum found at:\n");
	 
	 printf ("iteration: %d %.5f %.5f %10.5f\n", int(iter),
		 gsl_vector_get (s->x, 0), 
		 gsl_vector_get (s->x, 1), 
		 s->f);
	 
   } while (status == GSL_CONTINUE && iter < 14);

   std::cout << "============= terminating at iter " << iter << " with status "
	     << status;
   if (status == GSL_ENOPROG)
      std::cout << " NO_PROGRESS";
   std::cout << " ============ " << std::endl;

   // input becomes output
   for (unsigned int i=0; i<input_circles.size(); i++) {
      std::cout << "replacing " << input_circles[i].pos.x << " with " << gsl_vector_get(s->x, 2*i  ) << "    ";
      std::cout << "replacing " << input_circles[i].pos.y << " with " << gsl_vector_get(s->x, 2*i+1) << std::endl;
      input_circles[i].pos.x = gsl_vector_get(s->x, 2*i  );
      input_circles[i].pos.y = gsl_vector_get(s->x, 2*i+1);
   }
   
   gsl_multimin_fdfminimizer_free (s);
   gsl_vector_free (x);
   
}


std::vector<lbg_info_t::residue_circle_t>
lbg_info_t::optimise_residue_circles::solution() const {
   return input_circles;
}

double
lbg_info_t::optimise_residue_circles::f(const gsl_vector *v, void *params) {

   double score = 0.0;

   double atom_scale = 1.0;

   lbg_info_t::optimise_residue_circles *orc = static_cast<optimise_residue_circles *>(params);

   for (unsigned int i=0; i<orc->input_circles.size(); i++) { 

      // first score against the fixed (ligand) atoms and ring centres
      // (Do I need to test for not being a hydrogen?)
      /// 
      for (unsigned int iat=0; iat<orc->mol.atoms.size(); iat++) {
	 double d_pt_1 = gsl_vector_get(v, 2*i  ) - orc->mol.atoms[iat].atom_position.x;
	 double d_pt_2 = gsl_vector_get(v, 2*i+1) - orc->mol.atoms[iat].atom_position.y;
	 double d2 = d_pt_1 * d_pt_1 + d_pt_2 * d_pt_2;
	 if (d2 < 1) d2 = 1;
	 if (0) 
	    std::cout << "mol atom " << iat << " adding score " << 1.0/sqrt(d2)
		      << " d2 is " << d2 << " from " << gsl_vector_get(v, 2*i  )
		      << " - " << orc->mol.atoms[iat].atom_position.x 
		      << "\n";
	 score += 1.0/sqrt(d2);
      }

      // score against the ligand ring centres
      //
      for (unsigned int ib=0; ib<orc->mol.bonds.size(); ib++) { 
	 if (orc->mol.bonds[ib].have_centre_pos()) {
	    double d_pt_1 = gsl_vector_get(v, 2*i  ) - orc->mol.bonds[ib].centre_pos().x;
	    double d_pt_2 = gsl_vector_get(v, 2*i+1) - orc->mol.bonds[ib].centre_pos().y;
	    double d2 = d_pt_1 * d_pt_1 + d_pt_2 * d_pt_2;
	    if (d2 < 1) d2 = 1;
	    score += 1.0/sqrt(d2);
	 }
      }


      // score against the other atoms
      //
      for (unsigned int ic=0; ic<orc->input_circles.size(); ic++) {
	 if (ic != i) {
	    double d_pt_1 = gsl_vector_get(v, 2*i  ) - gsl_vector_get(v, 2*ic  );
	    double d_pt_2 = gsl_vector_get(v, 2*i+1) - gsl_vector_get(v, 2*ic+1);
	    double d2 = d_pt_1 * d_pt_1 + d_pt_2 * d_pt_2;
	    if (d2 < 1) d2 = 1;
	    score += 1.0/sqrt(d2);
	 }
      }

      // score against the original position (1/d^2 function)
      // 
      double k = 0.0003;
      double d_1 = gsl_vector_get(v, 2*i  ) - orc->input_circles[i].pos.x;
      double d_2 = gsl_vector_get(v, 2*i+1) - orc->input_circles[i].pos.y;
      score += k * (d_1*d_1 + d_2*d_2);

   }

   return score;
}


void
lbg_info_t::optimise_residue_circles::df(const gsl_vector *v, void *params, gsl_vector *df) {

   double atom_scale = 1.0;

   lbg_info_t::optimise_residue_circles *orc = static_cast<optimise_residue_circles *> (params);

   // initialise the derivatives vector
   for (unsigned int i=0; i<orc->input_circles.size(); i++) { 
      gsl_vector_set(df, 2*i,   0.0);
      gsl_vector_set(df, 2*i+1, 0.0);
   }

   for (unsigned int i=0; i<orc->input_circles.size(); i++) { 

      // first score against the fixed (ligand) atoms...
      // 
      for (unsigned int iat=0; iat<orc->mol.atoms.size(); iat++) {
	 double d_pt_1 = gsl_vector_get(v, 2*i  ) - orc->mol.atoms[iat].atom_position.x;
	 double d_pt_2 = gsl_vector_get(v, 2*i+1) - orc->mol.atoms[iat].atom_position.y;
	 double d2 = d_pt_1 * d_pt_1 + d_pt_2 * d_pt_2;
	 double d = sqrt(d2);
	 double df_part_x = (-1.0/d2)*(d_pt_1/d);
	 double df_part_y = (-1.0/d2)*(d_pt_2/d);
	 gsl_vector_set(df, 2*i,   gsl_vector_get(df, 2*i  ) + df_part_x);
	 gsl_vector_set(df, 2*i+1, gsl_vector_get(df, 2*i+1) + df_part_y);
      }

       //  score against the ligand ring centres
       //
       for (unsigned int ib=0; ib<orc->mol.bonds.size(); ib++) { 
 	 if (orc->mol.bonds[ib].have_centre_pos()) {
 	    double d_pt_1 = gsl_vector_get(v, 2*i  ) - orc->mol.bonds[ib].centre_pos().x;
 	    double d_pt_2 = gsl_vector_get(v, 2*i+1) - orc->mol.bonds[ib].centre_pos().y;
 	    double d2 = d_pt_1 * d_pt_1 + d_pt_2 * d_pt_2;
 	    double d = sqrt(d2);
 	    double df_part_x = (-1.0/d2)*(d_pt_1/d);
 	    double df_part_y = (-1.0/d2)*(d_pt_2/d);
	    gsl_vector_set(df, 2*i,   gsl_vector_get(df, 2*i  ) + df_part_x);
 	    gsl_vector_set(df, 2*i+1, gsl_vector_get(df, 2*i+1) + df_part_y);
	 }
       }

       // score against the other atoms
       //
       for (unsigned int ic=0; ic<orc->input_circles.size(); ic++) {
	  if (ic != i) {
	     double d_pt_1 = gsl_vector_get(v, 2*i  ) - gsl_vector_get(v, 2*ic  );
	     double d_pt_2 = gsl_vector_get(v, 2*i+1) - gsl_vector_get(v, 2*ic+1);
	     double d2 = d_pt_1 * d_pt_1 + d_pt_2 * d_pt_2;
	     double d = sqrt(d2);
	     double df_part_x = (-1.0/d2)*(d_pt_1/d);
	     double df_part_y = (-1.0/d2)*(d_pt_2/d);
	     gsl_vector_set(df, 2*i,   gsl_vector_get(df, 2*i  ) + df_part_x);
	     gsl_vector_set(df, 2*i+1, gsl_vector_get(df, 2*i+1) + df_part_y);	    
	  }
       }

       // score against the original position (1/d^2 function)
       //
       double k = 0.0003;
       double df_part_x = -2.0 * (gsl_vector_get(v, 2*i  ) - orc->input_circles[i].pos.x);
       double df_part_y = -2.0 * (gsl_vector_get(v, 2*i+1) - orc->input_circles[i].pos.y);
       gsl_vector_set(df, 2*i,   gsl_vector_get(df, 2*i  ) - k * df_part_x);
       gsl_vector_set(df, 2*i+1, gsl_vector_get(df, 2*i+1) - k * df_part_y);	    
   }
}

void
lbg_info_t::optimise_residue_circles::fdf(const gsl_vector *x, void *params, double *f_in, gsl_vector *df_in) {

   *f_in = f(x, params);
   df(x, params, df_in);
}
