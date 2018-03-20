
#ifdef HAVE_GSL
#ifdef HAVE_CXX_THREAD

#include <thread>
#include <chrono>

#include <iostream>
#include <iomanip>
#include <algorithm> // for sort
#include <stdexcept>

#include "simple-restraint.hh"
#include "process-dfs-in-range.hh" // put this in simple-restraint.hh FIXME

// parallel version of my_df()
void
coot::split_the_gradients_with_threads(const gsl_vector *v,
				       restraints_container_t *restraints_p,
				       gsl_vector *df) {

   // noisy but useful?
   // std::cout << "----------------- split_the_gradients_with_threads() " << std::endl;

#ifdef HAVE_CXX_THREAD

   // --------------------------------------- restraints --------------------------------------

   // auto tp_1a = std::chrono::high_resolution_clock::now();
   unsigned int n_t = restraints_p->n_threads;
   unsigned int restraints_size = restraints_p->size();
   std::vector<std::vector<std::size_t> > restraints_indices(n_t);
   std::vector<std::vector<double> > results(n_t);

   // each restraint_index vector will contain about
   // restraints_size/n_threads restraint indices
   //
   int r_reserve_size = std::lround(static_cast<float>(restraints_size)/static_cast<float>(n_t)) + 2;

   // First fill the restraints indices vectors
   //
   for (std::size_t n=0; n<restraints_indices.size(); n++)
      restraints_indices[n].reserve(r_reserve_size);

   unsigned int i_thread = 0; // insert to vector for this thread
   for (unsigned int ir=0; ir<restraints_size; ir++) {
      restraints_indices[i_thread].push_back(ir);
      ++i_thread;
      if (i_thread==n_t) i_thread=0;
   }

   if (false) { // debug thread-based restraints splitting
      for (std::size_t ii=0; ii<restraints_indices.size(); ii++) {
	 std::cout << "::: thread " << ii << " has restraints ";
	 for (std::size_t jj=0; jj<restraints_indices[ii].size(); jj++)
	    std::cout << " " << restraints_indices[ii][jj];
	 std::cout << std::endl;
      }
   }

   // Now make space for the df results, vectors of size 3*n_atoms
   //
   unsigned int n_variables = restraints_p->n_variables();
   for (std::size_t ii=0; ii<n_t; ii++)
      results[ii] = std::vector<double>(n_variables, 0);

   // use restraints in the range of restraints_indices[i}
   // to fill results
   //
   std::vector<std::thread> threads;
   for (std::size_t ii=0; ii<restraints_indices.size(); ii++) {
      threads.push_back(std::thread(process_dfs_in_range,
				    std::cref(restraints_indices[ii]),
				    restraints_p, v,
				    std::ref(results[ii])));
   }
   for (std::size_t ii=0; ii<restraints_indices.size(); ii++)
      threads.at(ii).join();

   // consolidate
   for (std::size_t i_thread=0; i_thread<n_t; i_thread++) {
      for (unsigned int i=0; i<n_variables; i++) {
	 if (results[i_thread][i] != 1110.0) {
	    // std::cout << "to df with idx " << i << " adding " << results[i_thread][i] << std::endl;
	    *gsl_vector_ptr(df, i) += results[i_thread][i];
	 }
      }
   }

   // auto tp_1b = std::chrono::high_resolution_clock::now();

   // --------------------------------------- map --------------------------------------

   if (restraints_p->include_map_terms()) {

      threads.clear();
      // like above, but this time we split the set of atoms into sets for each thread
      std::vector<std::vector<std::size_t> > atom_indices(n_t);
      i_thread = 0;
      unsigned int n = restraints_p->get_n_atoms();
      for (unsigned int ir=0; ir<n; ir++) {
	 atom_indices[i_thread].push_back(ir);
	 ++i_thread;
	 if (i_thread==n_t) i_thread=0;
      }

      auto tp_2 = std::chrono::high_resolution_clock::now();
      // we can manipulate df directly because (unlike restraints) each atom can
      // only touch 3 indices in the df vector - and do so uniquely.
      //
      for (std::size_t ii=0; ii<atom_indices.size(); ii++)
	 threads.push_back(std::thread(process_electron_density_dfs_for_atoms,
				       atom_indices[ii], restraints_p, v, df));
      for (std::size_t ii=0; ii<atom_indices.size(); ii++)
	 threads.at(ii).join();

      /*
      auto tp_3 = std::chrono::high_resolution_clock::now();
      auto d32 = chrono::duration_cast<chrono::microseconds>(tp_3 - tp_2).count();
      auto d_1ab = chrono::duration_cast<chrono::microseconds>(tp_1b - tp_1a).count();
      std::cout << "info:: df geometry: " << d_1ab << " microseconds\n";
      std::cout << "info:: df electron_density: " << d32 << " microseconds\n";
      */
   }

#endif // HAVE_CXX_THREAD

}


// fill results
void
coot::process_dfs_in_range(const std::vector<std::size_t> &restraints_indices,
			   coot::restraints_container_t *restraints_p,
			   const gsl_vector *v,
			   std::vector<double> &results) {

   for (std::size_t i=0; i<restraints_indices.size(); i++) {

      const simple_restraint &rest = restraints_p->at(restraints_indices[i]);

      if (false)
	 std::cout << "process_dfs_in_range() i " << i << " restraint index " << restraints_indices[i]
		   << " " << rest << std::endl;

      if (restraints_p->restraints_usage_flag & coot::BONDS_MASK)
	 if (rest.restraint_type == coot::BOND_RESTRAINT)
	    process_dfs_bond(rest, v, results);
	 
      if (restraints_p->restraints_usage_flag & coot::ANGLES_MASK)
	 if (rest.restraint_type == coot::ANGLE_RESTRAINT)
	    process_dfs_angle(rest, v, results);

      // Torsions are not yet turned on in the constructor
      if (restraints_p->restraints_usage_flag & coot::TORSIONS_MASK)
	 if (rest.restraint_type == coot::TORSION_RESTRAINT)
	    process_dfs_torsion(rest, v, results);

      if (restraints_p->restraints_usage_flag & coot::CHIRAL_VOLUME_MASK)
	 if (rest.restraint_type == coot::CHIRAL_VOLUME_RESTRAINT)
	    process_dfs_chiral_volume(rest, v, results);

      if (restraints_p->restraints_usage_flag & coot::PLANES_MASK)
	 if (rest.restraint_type == coot::PLANE_RESTRAINT)
	    process_dfs_plane(rest, v, results);

      if (restraints_p->restraints_usage_flag & coot::GEMAN_MCCLURE_DISTANCE_MASK)
	 if (rest.restraint_type == coot::GEMAN_MCCLURE_DISTANCE_RESTRAINT)
	    process_dfs_geman_mcclure_distance(rest, restraints_p->geman_mcclure_alpha, v, results);

      if (restraints_p->restraints_usage_flag & coot::NON_BONDED_MASK) {
	 if (rest.restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT) {
	    if (rest.nbc_function == simple_restraint::LENNARD_JONES) {
	       process_dfs_non_bonded_lennard_jones(rest, restraints_p->lennard_jones_epsilon, v, results);
	    } else {
	       process_dfs_non_bonded(rest, v, results);
	    }
	 }
      }

      if (restraints_p->restraints_usage_flag & coot::TRANS_PEPTIDE_MASK)
	 if (rest.restraint_type == coot::TRANS_PEPTIDE_RESTRAINT)
	    process_dfs_trans_peptide(rest, v, results);

      if (restraints_p->restraints_usage_flag & coot::RAMA_PLOT_MASK)
	 if (rest.restraint_type == coot::RAMACHANDRAN_RESTRAINT)
	    process_dfs_rama(rest, restraints_p, v, results);

      if (rest.restraint_type == coot::TARGET_POS_RESTRANT)
	 process_dfs_target_position(rest, restraints_p->log_cosh_target_distance_scale_factor, v, results);


   }
}

void
coot::process_dfs_bond(const coot::simple_restraint &restraint,
		       const gsl_vector *v,
		       std::vector<double> &results) { // fill results

   const double &target_val = restraint.target_value;

   int idx_1 = 3*restraint.atom_index_1;
   int idx_2 = 3*restraint.atom_index_2;
   clipper::Coord_orth a1(gsl_vector_get(v,idx_1),
			  gsl_vector_get(v,idx_1+1),
			  gsl_vector_get(v,idx_1+2));
   clipper::Coord_orth a2(gsl_vector_get(v,idx_2),
			  gsl_vector_get(v,idx_2+1),
			  gsl_vector_get(v,idx_2+2));

   double b_i_sqrd = (a1-a2).lengthsq();
   b_i_sqrd = b_i_sqrd > 0.01 ? b_i_sqrd : 0.01;  // Garib's stabilization

   double weight = 1/(restraint.sigma * restraint.sigma);

   double constant_part = 2.0*weight * (1 - target_val * f_inv_fsqrt(b_i_sqrd));

   double x_k_contrib = constant_part*(a1.x()-a2.x());
   double y_k_contrib = constant_part*(a1.y()-a2.y());
   double z_k_contrib = constant_part*(a1.z()-a2.z());

   double x_l_contrib = constant_part*(a2.x()-a1.x());
   double y_l_contrib = constant_part*(a2.y()-a1.y());
   double z_l_contrib = constant_part*(a2.z()-a1.z());

   if (! restraint.fixed_atom_flags[0]) {
      results[idx_1  ] += x_k_contrib;
      results[idx_1+1] += y_k_contrib;
      results[idx_1+2] += z_k_contrib;
   }

   if (! restraint.fixed_atom_flags[1]) {
      results[idx_2  ] += x_l_contrib;
      results[idx_2+1] += y_l_contrib;
      results[idx_2+2] += z_l_contrib;
   }

}

void
coot::process_dfs_angle(const coot::simple_restraint &restraint,
		       const gsl_vector *v,
		       std::vector<double> &results) { // fill results

   int idx; 

   double a;
   double b;
   double l_over_a_sqd;
   double l_over_b_sqd;
   double l_ab;

   double a_dot_b;
   double cos_theta;
   double theta;
   double prem;

   double x_k_contrib;
   double y_k_contrib;
   double z_k_contrib;

   double x_m_contrib;
   double y_m_contrib;
   double z_m_contrib;

   double term1x;
   double term1y;
   double term1z;

   double term2x;
   double term2y;
   double term2z;

   double x_l_mid_contrib;
   double y_l_mid_contrib;
   double z_l_mid_contrib;

   double weight;
   double ds_dth;
   double w_ds_dth;

   double target_value = restraint.target_value*DEGTORAD;

   idx = 3*(restraint.atom_index_1); 
   clipper::Coord_orth k(gsl_vector_get(v,idx), 
			 gsl_vector_get(v,idx+1), 
			 gsl_vector_get(v,idx+2));
   idx = 3*(restraint.atom_index_2); 
   clipper::Coord_orth l(gsl_vector_get(v,idx), 
			 gsl_vector_get(v,idx+1), 
			 gsl_vector_get(v,idx+2));
   idx = 3*(restraint.atom_index_3); 
   clipper::Coord_orth m(gsl_vector_get(v,idx), 
			 gsl_vector_get(v,idx+1), 
			 gsl_vector_get(v,idx+2));

   clipper::Coord_orth a_vec = (k - l); 
   clipper::Coord_orth b_vec = (m - l);  
   a = sqrt(a_vec.lengthsq());
   b = sqrt(b_vec.lengthsq()); 

   // Garib's stabilization
   if (a < 0.01) { 
      a = 0.01;
      a_vec = clipper::Coord_orth(0.01, 0.01, 0.01);
   }
   if (b < 0.01) { 
      b = 0.01;
      b_vec = clipper::Coord_orth(0.01, 0.01, 0.01);
   }
	    
   l_over_a_sqd = 1/(a*a);
   l_over_b_sqd = 1/(b*b);
   l_ab = 1/(a*b); 

   // for the end atoms: 
   // \frac{\partial \theta}{\partial x_k} =
   //    -\frac{1}{sin\theta} [(x_l-x_k)cos\theta + \frac{x_m-x_l}{ab}]
	 
   a_dot_b = clipper::Coord_orth::dot(a_vec,b_vec);
   cos_theta = a_dot_b/(a*b);
   // we need to stabilize cos_theta
   if (cos_theta < -1.0) cos_theta = -1.0;
   if (cos_theta >  1.0) cos_theta =  1.0;
   theta = acos(cos_theta); 

   // we need to stabilize $\theta$ too.
   a_dot_b = clipper::Coord_orth::dot(a_vec, b_vec); 
   cos_theta = a_dot_b/(a*b);
   if (cos_theta < -1) cos_theta = -1.0;
   if (cos_theta >  1) cos_theta =  1.0;
   theta = acos(cos_theta);

   // theta = theta > 0.001 ? theta : 0.001;
   if (theta < 0.001) theta = 0.001; // it was never -ve.

   prem = -1/sin(theta); 

   // The end atoms:
   x_k_contrib = prem*(cos_theta*(l.x()-k.x())*l_over_a_sqd + l_ab*(m.x()-l.x()));
   y_k_contrib = prem*(cos_theta*(l.y()-k.y())*l_over_a_sqd + l_ab*(m.y()-l.y()));
   z_k_contrib = prem*(cos_theta*(l.z()-k.z())*l_over_a_sqd + l_ab*(m.z()-l.z()));
	    
   x_m_contrib = prem*(cos_theta*(l.x()-m.x())*l_over_b_sqd + l_ab*(k.x()-l.x()));
   y_m_contrib = prem*(cos_theta*(l.y()-m.y())*l_over_b_sqd + l_ab*(k.y()-l.y()));
   z_m_contrib = prem*(cos_theta*(l.z()-m.z())*l_over_b_sqd + l_ab*(k.z()-l.z()));

   // For the middle atom, we have more cross terms in 
   // the derivatives of ab and a_dot_b.
   // 
   // I will split it up so that it is easier to read: 
   // 
   term1x = (-cos_theta*(l.x()-k.x())*l_over_a_sqd) -cos_theta*(l.x()-m.x())*l_over_b_sqd;
   term1y = (-cos_theta*(l.y()-k.y())*l_over_a_sqd) -cos_theta*(l.y()-m.y())*l_over_b_sqd;
   term1z = (-cos_theta*(l.z()-k.z())*l_over_a_sqd) -cos_theta*(l.z()-m.z())*l_over_b_sqd;

   term2x = (-(k.x()-l.x())-(m.x()-l.x()))*l_ab;
   term2y = (-(k.y()-l.y())-(m.y()-l.y()))*l_ab; 
   term2z = (-(k.z()-l.z())-(m.z()-l.z()))*l_ab; 

   x_l_mid_contrib = prem*(term1x + term2x); 
   y_l_mid_contrib = prem*(term1y + term2y); 
   z_l_mid_contrib = prem*(term1z + term2z);

   // and finally the term that is common to all, $\frac{\partial S}{\partial \theta}
   // dS/d(th).
   //
   weight = 1/(restraint.sigma * restraint.sigma);
   ds_dth = 2*(theta - target_value)*RADTODEG*RADTODEG;
   w_ds_dth = weight * ds_dth; 

   if (!restraint.fixed_atom_flags[0]) { 
      idx = 3*(restraint.atom_index_1);
      // gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + x_k_contrib*w_ds_dth); 
      // gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + y_k_contrib*w_ds_dth); 
      // gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + z_k_contrib*w_ds_dth);
      results[idx  ] += x_k_contrib*w_ds_dth;
      results[idx+1] += y_k_contrib*w_ds_dth;
      results[idx+2] += z_k_contrib*w_ds_dth;
   }
   if (!restraint.fixed_atom_flags[2]) { 
      idx = 3*(restraint.atom_index_3);
      // gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + x_m_contrib*w_ds_dth); 
      // gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + y_m_contrib*w_ds_dth); 
      // gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + z_m_contrib*w_ds_dth); 
      results[idx  ] += x_m_contrib*w_ds_dth;
      results[idx+1] += y_m_contrib*w_ds_dth;
      results[idx+2] += z_m_contrib*w_ds_dth;
   }

   // and mid atom
   if (!restraint.fixed_atom_flags[1]) { 
      idx = 3*(restraint.atom_index_2);
      // gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + x_l_mid_contrib*w_ds_dth); 
      // gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + y_l_mid_contrib*w_ds_dth); 
      // gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + z_l_mid_contrib*w_ds_dth); 
      results[idx  ] += x_l_mid_contrib*w_ds_dth;
      results[idx+1] += y_l_mid_contrib*w_ds_dth;
      results[idx+2] += z_l_mid_contrib*w_ds_dth;
   }

}

void
coot::process_dfs_torsion(const coot::simple_restraint &this_restraint,
			  const gsl_vector *v,
			  std::vector<double> &results) { // fill results

   int n_torsion_restr = 0; 
   int idx; 

   idx = 3*(this_restraint.atom_index_1); 
   clipper::Coord_orth P1(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(this_restraint.atom_index_2); 
   clipper::Coord_orth P2(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(this_restraint.atom_index_3); 
   clipper::Coord_orth P3(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(this_restraint.atom_index_4); 
   clipper::Coord_orth P4(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));

   try {
      coot::distortion_torsion_gradients_t dtg =
	 fill_distortion_torsion_gradients(P1, P2, P3, P4);

      if (dtg.zero_gradients) {

	 std::cout << "debug:: in process_dfs_torsion zero_gradients " << std::endl;

      } else {

	 //
	 // use period

	 double diff = 99999.9;
	 double tdiff;
	 double trial_target;
	 int per = this_restraint.periodicity;

	 if (dtg.theta < 0.0) dtg.theta += 360.0; 

	 for(int iper=0; iper<per; iper++) { 
	    trial_target = this_restraint.target_value + double(iper)*360.0/double(per); 
	    if (trial_target >= 360.0) trial_target -= 360.0; 
	    tdiff = dtg.theta - trial_target;
	    if (tdiff < -180) tdiff += 360;
	    if (tdiff >  180) tdiff -= 360;
	    // std::cout << "   iper: " << iper << "   " << dtg.theta << "   " << trial_target << "   " << tdiff << "   " << diff << std::endl;
	    if (fabs(tdiff) < fabs(diff)) { 
	       diff = tdiff;
	    }
	 }
	 if (diff < -180.0) { 
	    diff += 360.; 
	 } else { 
	    if (diff > 180.0) { 
	       diff -= 360.0; 
	    }
	 }
		  
	 if (true)
	    std::cout << "in df_torsion: dtg.theta is " << dtg.theta 
		      <<  " and target is " << this_restraint.target_value 
		      << " and diff is " << diff
		      << " and periodicity: " << this_restraint.periodicity << std::endl;

	 double tt = tan(clipper::Util::d2rad(dtg.theta));
	 double torsion_scale = (1.0/(1+tt*tt)) *
	    clipper::Util::rad2d(1.0);

	 double weight = 1/(this_restraint.sigma * this_restraint.sigma);

	 // 	       std::cout << "torsion weight: " << weight << std::endl;
	 // 	       std::cout << "torsion_scale : " << torsion_scale << std::endl; 
	 // 	       std::cout << "diff          : " << torsion_scale << std::endl; 	       

	 double xP1_contrib = 2.0*diff*dtg.dD_dxP1*torsion_scale * weight;
	 double xP2_contrib = 2.0*diff*dtg.dD_dxP2*torsion_scale * weight;
	 double xP3_contrib = 2.0*diff*dtg.dD_dxP3*torsion_scale * weight;
	 double xP4_contrib = 2.0*diff*dtg.dD_dxP4*torsion_scale * weight;

	 double yP1_contrib = 2.0*diff*dtg.dD_dyP1*torsion_scale * weight;
	 double yP2_contrib = 2.0*diff*dtg.dD_dyP2*torsion_scale * weight;
	 double yP3_contrib = 2.0*diff*dtg.dD_dyP3*torsion_scale * weight;
	 double yP4_contrib = 2.0*diff*dtg.dD_dyP4*torsion_scale * weight;

	 double zP1_contrib = 2.0*diff*dtg.dD_dzP1*torsion_scale * weight;
	 double zP2_contrib = 2.0*diff*dtg.dD_dzP2*torsion_scale * weight;
	 double zP3_contrib = 2.0*diff*dtg.dD_dzP3*torsion_scale * weight;
	 double zP4_contrib = 2.0*diff*dtg.dD_dzP4*torsion_scale * weight;
	    
	 if (! this_restraint.fixed_atom_flags[0]) { 
	    idx = 3*(this_restraint.atom_index_1);
	    // *gsl_vector_ptr(df, idx  ) += xP1_contrib;
	    // *gsl_vector_ptr(df, idx+1) += yP1_contrib;
	    // *gsl_vector_ptr(df, idx+2) += zP1_contrib;

	    results[idx  ] += xP1_contrib;
	    results[idx+1] += yP1_contrib;
	    results[idx+2] += zP1_contrib;
	 }

	 if (! this_restraint.fixed_atom_flags[1]) { 
	    idx = 3*(this_restraint.atom_index_2);
	    // *gsl_vector_ptr(df, idx  ) += xP2_contrib;
	    // *gsl_vector_ptr(df, idx+1) += yP2_contrib;
	    // *gsl_vector_ptr(df, idx+2) += zP2_contrib;

	    results[idx  ] += xP2_contrib;
	    results[idx+1] += yP2_contrib;
	    results[idx+2] += zP2_contrib;
	 }

	 if (! this_restraint.fixed_atom_flags[2]) { 
	    idx = 3*(this_restraint.atom_index_3);
	    // *gsl_vector_ptr(df, idx  ) += xP3_contrib;
	    // *gsl_vector_ptr(df, idx+1) += yP3_contrib;
	    // *gsl_vector_ptr(df, idx+2) += zP3_contrib;

	    results[idx  ] += xP3_contrib;
	    results[idx+1] += yP3_contrib;
	    results[idx+2] += zP3_contrib;
	 }

	 if (! this_restraint.fixed_atom_flags[3]) { 
	    idx = 3*(this_restraint.atom_index_4);
	    // *gsl_vector_ptr(df, idx  ) += xP4_contrib;
	    // *gsl_vector_ptr(df, idx+1) += yP4_contrib;
	    // *gsl_vector_ptr(df, idx+2) += zP4_contrib;

	    results[idx  ] += xP4_contrib;
	    results[idx+1] += yP4_contrib;
	    results[idx+2] += zP4_contrib;
	 }
      }
   }
   catch (const std::runtime_error &rte) {
      std::cout << "Caught runtime_error" << rte.what() << std::endl;
   }
} 

void
coot::process_dfs_chiral_volume(const coot::simple_restraint &restraint,
				const gsl_vector *v,
				std::vector<double> &results) { // fill results


   double cv;
   double distortion;
   
   int idx = 3*restraint.atom_index_centre;
   clipper::Coord_orth centre(gsl_vector_get(v, idx),
			      gsl_vector_get(v, idx+1),
			      gsl_vector_get(v, idx+2));

   idx = 3*( restraint.atom_index_1);
   clipper::Coord_orth a1(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));
   idx = 3*( restraint.atom_index_2);
   clipper::Coord_orth a2(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));
   idx = 3*( restraint.atom_index_3);
   clipper::Coord_orth a3(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));

   clipper::Coord_orth a = a1 - centre;
   clipper::Coord_orth b = a2 - centre;
   clipper::Coord_orth c = a3 - centre;

   cv = clipper::Coord_orth::dot(a, clipper::Coord_orth::cross(b,c));

   distortion = cv - restraint.target_chiral_volume;
	    
   // 	    std::cout << "---- xxx ---- DEBUG:: chiral volume deriv: " 
   // 		      << cv << " chiral distortion " 
   // 		      << distortion << "\n";
   // distortion /= (restraint.sigma * restraint.sigma);

   double P0_x_contrib =
      - (b.y()*c.z() - b.z()*c.y())
      - (a.z()*c.y() - a.y()*c.z())
      - (a.y()*b.z() - a.z()*b.y());
		  
   double P0_y_contrib = 
      - (b.z()*c.x() - b.x()*c.z())
      - (a.x()*c.z() - a.z()*c.x())
      - (a.z()*b.x() - a.x()*b.z());

   double P0_z_contrib = 
      - (b.x()*c.y() - b.y()*c.x())
      - (a.y()*c.x() - a.x()*c.y())
      - (a.x()*b.y() - a.y()*b.x());

   double P1_x_contrib = b.y()*c.z() - b.z()*c.y();
   double P1_y_contrib = b.z()*c.x() - b.x()*c.z();
   double P1_z_contrib = b.x()*c.y() - b.y()*c.x();

   double P2_x_contrib = a.z()*c.y() - a.y()*c.z();
   double P2_y_contrib = a.x()*c.z() - a.z()*c.x();
   double P2_z_contrib = a.y()*c.x() - a.x()*c.y();

   double P3_x_contrib = a.y()*b.z() - a.z()*b.y();
   double P3_y_contrib = a.z()*b.x() - a.x()*b.z();
   double P3_z_contrib = a.x()*b.y() - a.y()*b.x();

   double s = 2*distortion/(restraint.sigma * restraint.sigma);

   if (!restraint.fixed_atom_flags[0]) { 
      idx = 3*( restraint.atom_index_centre);
      // gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + s * P0_x_contrib);
      // gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + s * P0_y_contrib);
      // gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + s * P0_z_contrib);
      results[idx  ] += s * P0_x_contrib;
      results[idx+1] += s * P0_y_contrib;
      results[idx+2] += s * P0_z_contrib;
   }
	       
   if (!restraint.fixed_atom_flags[1]) { 
      idx = 3*( restraint.atom_index_1);
      // gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + s * P1_x_contrib); 
      // gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + s * P1_y_contrib); 
      // gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + s * P1_z_contrib);
      results[idx  ] += s * P1_x_contrib;
      results[idx+1] += s * P1_y_contrib;
      results[idx+2] += s * P1_z_contrib;
   }

   if (!restraint.fixed_atom_flags[2]) { 
      idx = 3*( restraint.atom_index_2);
      // gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + s * P2_x_contrib); 
      // gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + s * P2_y_contrib); 
      // gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + s * P2_z_contrib);
      results[idx  ] += s * P2_x_contrib;
      results[idx+1] += s * P2_y_contrib;
      results[idx+2] += s * P2_z_contrib;
   }

   if (!restraint.fixed_atom_flags[3]) { 
      idx = 3*( restraint.atom_index_3);
      // gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + s * P3_x_contrib); 
      // gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + s * P3_y_contrib); 
      // gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + s * P3_z_contrib);
      results[idx  ] += s * P3_x_contrib;
      results[idx+1] += s * P3_y_contrib;
      results[idx+2] += s * P3_z_contrib;
   }
}

void
coot::process_dfs_plane(const coot::simple_restraint &plane_restraint,
		       const gsl_vector *v,
		       std::vector<double> &results) { // fill results

   int idx; 

   double devi_len;
   coot::plane_distortion_info_t plane_info =
      distortion_score_plane_internal(plane_restraint, v);
   int n_plane_atoms = plane_restraint.plane_atom_index.size();

   for (int j=0; j<n_plane_atoms; j++) {
      if (! plane_restraint.fixed_atom_flags[j] ) {
	 idx = 3*plane_restraint.plane_atom_index[j].first;
	 devi_len =
	    plane_info.abcd[0]*gsl_vector_get(v,idx  ) + 
	    plane_info.abcd[1]*gsl_vector_get(v,idx+1) +
	    plane_info.abcd[2]*gsl_vector_get(v,idx+2) -
	    plane_info.abcd[3];

	 double weight = 1/(plane_restraint.plane_atom_index[j].second *
			    plane_restraint.plane_atom_index[j].second);

	 clipper::Grad_orth<double> d(2.0 * weight * devi_len * plane_info.abcd[0],
				      2.0 * weight * devi_len * plane_info.abcd[1],
				      2.0 * weight * devi_len * plane_info.abcd[2]);

	 // *gsl_vector_ptr(df, idx  ) += d.dx();
	 // *gsl_vector_ptr(df, idx+1) += d.dy();
	 // *gsl_vector_ptr(df, idx+2) += d.dz();

	 results[idx  ] += d.dx();
	 results[idx+1] += d.dy();
	 results[idx+2] += d.dz();
      }
   }
}

void
coot::process_dfs_non_bonded(const coot::simple_restraint &this_restraint,
			     const gsl_vector *v,
			     std::vector<double> &results) { // fill results

   const double &target_val = this_restraint.target_value;

   // no need to calculate anything if both these atoms are non-moving
   //
   if (this_restraint.fixed_atom_flags[0] && this_restraint.fixed_atom_flags[1])
      return;

   int idx_1 = 3*this_restraint.atom_index_1;
   int idx_2 = 3*this_restraint.atom_index_2;

   // no need to calculate anything if both these atoms are non-moving
   //
   if (this_restraint.fixed_atom_flags[0] && this_restraint.fixed_atom_flags[1])
      return;

   // check for both-ways nbcs (seems OK)
   //
   // std::cout << "nbc: idx_1 " << idx_1 << " idx_2 " << idx_2 << std::endl;

   clipper::Coord_orth a1(gsl_vector_get(v,idx_1),
			  gsl_vector_get(v,idx_1+1),
			  gsl_vector_get(v,idx_1+2));

   clipper::Coord_orth a2(gsl_vector_get(v,idx_2),
			  gsl_vector_get(v,idx_2+1),
			  gsl_vector_get(v,idx_2+2));

   double b_i_sqrd = (a1-a2).lengthsq();

   if (b_i_sqrd < this_restraint.target_value * this_restraint.target_value) {

      double weight = 1.0/(this_restraint.sigma * this_restraint.sigma);

      double b_i = sqrt(b_i_sqrd);
      // b_i_sqrd = b_i_sqrd > 0.01 ? b_i_sqrd : 0.01;  // Garib's stabilization
      double constant_part = 2.0*weight * (1 - target_val * f_inv_fsqrt(b_i_sqrd));

      if (! this_restraint.fixed_atom_flags[0]) {

	 double x_k_contrib = constant_part*(a1.x()-a2.x());
	 double y_k_contrib = constant_part*(a1.y()-a2.y());
	 double z_k_contrib = constant_part*(a1.z()-a2.z());

 	 // *gsl_vector_ptr(df, idx_1  ) += x_k_contrib;
 	 // *gsl_vector_ptr(df, idx_1+1) += y_k_contrib;
 	 // *gsl_vector_ptr(df, idx_1+2) += z_k_contrib;

	 results[idx_1  ] += x_k_contrib;
	 results[idx_1+1] += y_k_contrib;
	 results[idx_1+2] += z_k_contrib;
      }

      if (! this_restraint.fixed_atom_flags[1]) {

 	 // *gsl_vector_ptr(df, idx_2  ) += x_l_contrib;
 	 // *gsl_vector_ptr(df, idx_2+1) += y_l_contrib;
 	 // *gsl_vector_ptr(df, idx_2+2) += z_l_contrib;

	 double x_l_contrib = constant_part*(a2.x()-a1.x());
	 double y_l_contrib = constant_part*(a2.y()-a1.y());
	 double z_l_contrib = constant_part*(a2.z()-a1.z());

	 results[idx_2  ] += x_l_contrib;
	 results[idx_2+1] += y_l_contrib;
	 results[idx_2+2] += z_l_contrib;
      }
   }
}

void
coot::process_dfs_non_bonded_lennard_jones(const coot::simple_restraint &this_restraint,
					   const double &lj_epsilon,
					   const gsl_vector *v,
					   std::vector<double> &results) { // fill results

   // theres a problem here somewhere.  Need to check this - change the NBC function
   // in make_non_bonded_contact_restraints()
   
   int idx_1 = 3*this_restraint.atom_index_1;
   int idx_2 = 3*this_restraint.atom_index_2;

   clipper::Coord_orth a1(gsl_vector_get(v,idx_1),
			  gsl_vector_get(v,idx_1+1),
			  gsl_vector_get(v,idx_1+2));

   clipper::Coord_orth a2(gsl_vector_get(v,idx_2),
			  gsl_vector_get(v,idx_2+1),
			  gsl_vector_get(v,idx_2+2));

   double lj_sigma = this_restraint.target_value;
   double max_dist = lj_sigma * 2.5; // 2.5 is conventional limit, i.e. ~3.5 * 2.5
   max_dist = 999.9; // does this match the one in the gradients? And the one in distortion score?
   
   double b_i_sqrd = (a1-a2).lengthsq();
   if (b_i_sqrd < 0.81) b_i_sqrd = 0.81; // stabilize (as per distortion score lj)

   if (b_i_sqrd < (max_dist * max_dist)) {

      // double lj_r_min = pow(2.0, 1.0/6.0) * lj_sigma; // precalculate this pow value - done
      double lj_r_min = 1.122462048309373 * lj_sigma;
      double lj_r = std::sqrt(b_i_sqrd);
      double alpha = lj_r_min/lj_r;
      double dalpha_dr = -lj_r_min/b_i_sqrd;

      double alpha_sqrd = lj_r_min*lj_r_min/b_i_sqrd;
      double alpha_up_5  = alpha_sqrd * alpha_sqrd * alpha;
      double alpha_up_6  = alpha_sqrd * alpha_sqrd * alpha_sqrd;
      double alpha_up_11 = alpha_up_6 * alpha_up_5;
      // double dVlj_dalpha = 12.0 * lj_epsilon * (std::pow(alpha, 11) - std::pow(alpha, 5));
      double dVlj_dalpha = 12.0 * lj_epsilon * (alpha_up_11 - alpha_up_5);

      double dVlj_dr = dVlj_dalpha * dalpha_dr;

      double constant_part = dVlj_dr/lj_r; // this is why we need the square root

      double delta_x = a1.x() - a2.x();
      double delta_y = a1.y() - a2.y();
      double delta_z = a1.z() - a2.z();

      if (! this_restraint.fixed_atom_flags[0]) {
 	 // *gsl_vector_ptr(df, idx_1  ) += x_k_contrib;
 	 // *gsl_vector_ptr(df, idx_1+1) += y_k_contrib;
 	 // *gsl_vector_ptr(df, idx_1+2) += z_k_contrib;

	 double x_k_contrib = constant_part*(a1.x()-a2.x());
	 double y_k_contrib = constant_part*(a1.y()-a2.y());
	 double z_k_contrib = constant_part*(a1.z()-a2.z());

	 results[idx_1  ] += x_k_contrib;
	 results[idx_1+1] += y_k_contrib;
	 results[idx_1+2] += z_k_contrib;

      }
      if (! this_restraint.fixed_atom_flags[1]) {
 	 // *gsl_vector_ptr(df, idx_2  ) += x_l_contrib;
 	 // *gsl_vector_ptr(df, idx_2+1) += y_l_contrib;
 	 // *gsl_vector_ptr(df, idx_2+2) += z_l_contrib;

	 double x_l_contrib = constant_part*(a2.x()-a1.x());
	 double y_l_contrib = constant_part*(a2.y()-a1.y());
	 double z_l_contrib = constant_part*(a2.z()-a1.z());

	 results[idx_2  ] += x_l_contrib;
	 results[idx_2+1] += y_l_contrib;
	 results[idx_2+2] += z_l_contrib;
      }
   }
}


void
coot::process_dfs_target_position(const coot::simple_restraint &restraint,
				  const double &log_cosh_target_distance_scale_factor,
				  const gsl_vector *v,
				  std::vector<double> &results) {

   double sigma = 0.04;
   int idx = 3*(restraint.atom_index_1);

   double scale = log_cosh_target_distance_scale_factor;
   double top_out_dist = 4.0;  // Angstroms, needs tweaking?
   double k = 1.0 / top_out_dist;

   clipper::Coord_orth current_pos(gsl_vector_get(v,idx),
				   gsl_vector_get(v,idx+1),
				   gsl_vector_get(v,idx+2));
   double dist_x = gsl_vector_get(v, idx)   - restraint.atom_pull_target_pos[0];
   double dist_y = gsl_vector_get(v, idx+1) - restraint.atom_pull_target_pos[1];
   double dist_z = gsl_vector_get(v, idx+2) - restraint.atom_pull_target_pos[2];

   double dist = clipper::Coord_orth::length(current_pos, restraint.atom_pull_target_pos);
   double z = dist/top_out_dist;
   double e_2z = exp(2.0 * z);
   double tanh_z = (e_2z - 1) / (e_2z + 1);
   double constant_part = scale * k * tanh_z / dist;

   results[idx  ] += constant_part * dist_x;
   results[idx+1] += constant_part * dist_y;
   results[idx+2] += constant_part * dist_z;
   
}


void
coot::process_electron_density_dfs_for_atoms(const std::vector<std::size_t> &atom_indices,
					     const restraints_container_t *restraints_p,
					     const gsl_vector *v, gsl_vector *df) {

   for (std::size_t i=0; i<atom_indices.size(); i++) {
      const std::size_t &atom_idx = atom_indices[i];
      if (restraints_p->use_map_gradient_for_atom[atom_idx]) {
	 int idx = 3 * atom_idx;
	 clipper::Coord_orth ao(gsl_vector_get(v,idx),
				gsl_vector_get(v,idx+1),
				gsl_vector_get(v,idx+2));
	 clipper::Grad_orth<double> grad_orth = restraints_p->electron_density_gradient_at_point(ao);
	 double zs = restraints_p->Map_weight() * restraints_p->atom_z_occ_weight[atom_idx];
	 *gsl_vector_ptr(df, idx  ) -= zs * grad_orth.dx();
	 *gsl_vector_ptr(df, idx+1) -= zs * grad_orth.dy();
	 *gsl_vector_ptr(df, idx+2) -= zs * grad_orth.dz();
      }
   }
}


void
coot::process_dfs_trans_peptide(const coot::simple_restraint &restraint,
					 const gsl_vector *v,
					 std::vector<double> &results) {

   int idx;
   idx = 3*(restraint.atom_index_1); 
   clipper::Coord_orth P1(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(restraint.atom_index_2); 
   clipper::Coord_orth P2(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(restraint.atom_index_3); 
   clipper::Coord_orth P3(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(restraint.atom_index_4); 
   clipper::Coord_orth P4(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));

   try {

      // if the bond angles are (near) linear, then the distortion gradients
      // are zero.
      //
      distortion_torsion_gradients_t dtg =
	 fill_distortion_torsion_gradients(P1, P2, P3, P4);

      if (dtg.zero_gradients) {

	 std::cout << "debug:: in process_dfs_trans_peptide zero_gradients " << std::endl;

      } else {

	 double diff = dtg.theta - restraint.target_value;
	 // because trans restraints - 180, (see distortion score notes)
	 if (diff >  180) diff -= 360;
	 if (diff < -180) diff += 360;

	 if (false)
	    std::cout << "in process_dfs_trans_peptide: dtg.theta is " << dtg.theta 
		      <<  " and target is " << restraint.target_value 
		      << " and diff is " << diff 
		      << " and periodicity: " << restraint.periodicity << std::endl;

	 double tt = tan(clipper::Util::d2rad(dtg.theta));
	 double trans_peptide_scale = (1.0/(1+tt*tt)) *
	    clipper::Util::rad2d(1.0);

	 double weight = 1/(restraint.sigma * restraint.sigma);

	 // std::cout << "trans_peptide weight: " << weight << std::endl;
	 // std::cout << "trans_peptide_scale : " << trans_peptide_scale << std::endl; 
	 // std::cout << "diff          : " << trans_peptide_scale << std::endl; 	       

	 double xP1_contrib = 2.0*diff*dtg.dD_dxP1*trans_peptide_scale * weight;
	 double xP2_contrib = 2.0*diff*dtg.dD_dxP2*trans_peptide_scale * weight;
	 double xP3_contrib = 2.0*diff*dtg.dD_dxP3*trans_peptide_scale * weight;
	 double xP4_contrib = 2.0*diff*dtg.dD_dxP4*trans_peptide_scale * weight;

	 double yP1_contrib = 2.0*diff*dtg.dD_dyP1*trans_peptide_scale * weight;
	 double yP2_contrib = 2.0*diff*dtg.dD_dyP2*trans_peptide_scale * weight;
	 double yP3_contrib = 2.0*diff*dtg.dD_dyP3*trans_peptide_scale * weight;
	 double yP4_contrib = 2.0*diff*dtg.dD_dyP4*trans_peptide_scale * weight;

	 double zP1_contrib = 2.0*diff*dtg.dD_dzP1*trans_peptide_scale * weight;
	 double zP2_contrib = 2.0*diff*dtg.dD_dzP2*trans_peptide_scale * weight;
	 double zP3_contrib = 2.0*diff*dtg.dD_dzP3*trans_peptide_scale * weight;
	 double zP4_contrib = 2.0*diff*dtg.dD_dzP4*trans_peptide_scale * weight;

	 if (! restraint.fixed_atom_flags[0]) {
	    idx = 3*(restraint.atom_index_1);
	    // *gsl_vector_ptr(df, idx  ) += xP1_contrib;
	    // *gsl_vector_ptr(df, idx+1) += yP1_contrib;
	    // *gsl_vector_ptr(df, idx+2) += zP1_contrib;

	    results[idx  ] += xP1_contrib;
	    results[idx+1] += yP1_contrib;
	    results[idx+2] += zP1_contrib;
	 }

	 if (! restraint.fixed_atom_flags[1]) {
	    idx = 3*(restraint.atom_index_2);
	    // *gsl_vector_ptr(df, idx  ) += xP2_contrib;
	    // *gsl_vector_ptr(df, idx+1) += yP2_contrib;
	    // *gsl_vector_ptr(df, idx+2) += zP2_contrib;

	    results[idx  ] += xP2_contrib;
	    results[idx+1] += yP2_contrib;
	    results[idx+2] += zP2_contrib;
	 }

	 if (! restraint.fixed_atom_flags[2]) { 
	    idx = 3*(restraint.atom_index_3);
	    // *gsl_vector_ptr(df, idx  ) += xP3_contrib;
	    // *gsl_vector_ptr(df, idx+1) += yP3_contrib;
	    // *gsl_vector_ptr(df, idx+2) += zP3_contrib;

	    results[idx  ] += xP3_contrib;
	    results[idx+1] += yP3_contrib;
	    results[idx+2] += zP3_contrib;
	 }

	 if (! restraint.fixed_atom_flags[3]) { 
	    idx = 3*(restraint.atom_index_4);
	    // *gsl_vector_ptr(df, idx  ) += xP4_contrib;
	    // *gsl_vector_ptr(df, idx+1) += yP4_contrib;
	    // *gsl_vector_ptr(df, idx+2) += zP4_contrib;

	    results[idx  ] += xP4_contrib;
	    results[idx+1] += yP4_contrib;
	    results[idx+2] += zP4_contrib;
	 }
      }
   }
   catch (const std::runtime_error &rte) {
      std::cout << "Caught runtime_error" << rte.what() << std::endl;
   }
}

void
coot::process_dfs_geman_mcclure_distance(const coot::simple_restraint &this_restraint,
					 const double &alpha,
					 const gsl_vector *v,
					 std::vector<double> &results) {


   int idx_1 = 3*this_restraint.atom_index_1;
   int idx_2 = 3*this_restraint.atom_index_2;
   clipper::Coord_orth a1(gsl_vector_get(v,idx_1),
			  gsl_vector_get(v,idx_1+1),
			  gsl_vector_get(v,idx_1+2));
   clipper::Coord_orth a2(gsl_vector_get(v,idx_2),
			  gsl_vector_get(v,idx_2+1),
			  gsl_vector_get(v,idx_2+2));

   double b_i_sqrd = (a1-a2).lengthsq();
   b_i_sqrd = b_i_sqrd > 0.01 ? b_i_sqrd : 0.01;  // Garib's stabilization

   double b_i = sqrt(b_i_sqrd);
   const double &target_val = this_restraint.target_value;
   double weight = 1.0/(this_restraint.sigma * this_restraint.sigma);

   // Let z = (boi - bi)/sigma
   //    S_i = z^2/(1 + alpha * z^2)
   //
   double bit = b_i - this_restraint.target_value;
   double z = bit/this_restraint.sigma;

   double beta  = 1 + alpha * z * z;
   double d_Si_d_zi = 2.0 * z  / (beta * beta);
   double d_zi_d_bi = 1.0/this_restraint.sigma;
   double d_b_d_x_m = 1.0/b_i;

   double constant_part_gm = d_Si_d_zi * d_zi_d_bi * d_b_d_x_m;
   double constant_part = constant_part_gm;

   double constant_part_lsq = 2.0*weight * (1 - target_val * f_inv_fsqrt(b_i_sqrd));

   constant_part = constant_part_lsq / (beta * beta);

   // constant_part = constant_part_lsq; // force least squares

   // The final part is dependent on the coordinates:

   if (! this_restraint.fixed_atom_flags[0]) {
      // *gsl_vector_ptr(df, idx_1  ) += x_k_contrib;
      // *gsl_vector_ptr(df, idx_1+1) += y_k_contrib;
      // *gsl_vector_ptr(df, idx_1+2) += z_k_contrib;

      double x_k_contrib = constant_part*(a1.x()-a2.x());
      double y_k_contrib = constant_part*(a1.y()-a2.y());
      double z_k_contrib = constant_part*(a1.z()-a2.z());
      results[idx_1  ] += x_k_contrib;
      results[idx_1+1] += y_k_contrib;
      results[idx_1+2] += z_k_contrib;
   }

   if (! this_restraint.fixed_atom_flags[1]) {
      // *gsl_vector_ptr(df, idx_2  ) += x_l_contrib;
      // *gsl_vector_ptr(df, idx_2+1) += y_l_contrib;
      // *gsl_vector_ptr(df, idx_2+2) += z_l_contrib;

      double x_l_contrib = constant_part*(a2.x()-a1.x());
      double y_l_contrib = constant_part*(a2.y()-a1.y());
      double z_l_contrib = constant_part*(a2.z()-a1.z());

      results[idx_2  ] += x_l_contrib;
      results[idx_2+1] += y_l_contrib;
      results[idx_2+2] += z_l_contrib;
   }
}


void
coot::process_dfs_rama(const coot::simple_restraint &rama_restraint,
		       const coot::restraints_container_t *restraints,
		       const gsl_vector *v,
		       std::vector<double> &results) {

   try {

      int idx;

      idx = 3*(rama_restraint.atom_index_1);
      clipper::Coord_orth P1(gsl_vector_get(v,idx), 
			     gsl_vector_get(v,idx+1), 
			     gsl_vector_get(v,idx+2));
      idx = 3*(rama_restraint.atom_index_2); 
      clipper::Coord_orth P2(gsl_vector_get(v,idx), 
			     gsl_vector_get(v,idx+1), 
			     gsl_vector_get(v,idx+2));
      idx = 3*(rama_restraint.atom_index_3); 
      clipper::Coord_orth P3(gsl_vector_get(v,idx), 
			     gsl_vector_get(v,idx+1), 
			     gsl_vector_get(v,idx+2));
      idx = 3*(rama_restraint.atom_index_4); 
      clipper::Coord_orth P4(gsl_vector_get(v,idx), 
			     gsl_vector_get(v,idx+1), 
			     gsl_vector_get(v,idx+2));
      idx = 3*(rama_restraint.atom_index_5); 
      clipper::Coord_orth P5(gsl_vector_get(v,idx), 
			     gsl_vector_get(v,idx+1), 
			     gsl_vector_get(v,idx+2));

      clipper::Coord_orth a = P2 - P1; 
      clipper::Coord_orth b = P3 - P2; 
      clipper::Coord_orth c = P4 - P3;
      clipper::Coord_orth d = P5 - P4;

      // New assignements:
      // TRANS    psi    (2nd N) (2nd CA) (2nd C ) (3nd N)
      // TRANS    phi    (1st C) (2nd N ) (2nd CA) (2nd C) 
      // 
      // So Rama_atoms in this order:
      //   0       1        2      3         4
      //  P1      P2       P3     P4        P5
      // (1st C) (2nd N) (2nd CA) (2nd C) (3rd N)

      // ---------- phi ------------------
      // b*b * [ a.(bxc)/b ]
      double E = clipper::Coord_orth::dot(a,clipper::Coord_orth::cross(b,c)) *
	 sqrt( b.lengthsq() );

      // b*b * [ -a.c+(a.b)(b.c)/(b*b) ] = -a.c*b*b + (a.b)(b.c)
      double G = - clipper::Coord_orth::dot(a,c)*b.lengthsq()
	 + clipper::Coord_orth::dot(a,b)*clipper::Coord_orth::dot(b,c);

      double phi = clipper::Util::rad2d(atan2(E,G));
      if (phi < 180.0)
	 phi += 360.0;
      if (phi > 180.0)
	 phi -= 360.0;

      // ---------- psi ------------------
      // b*b * [ a.(bxc)/b ]
      double H = clipper::Coord_orth::dot(b, clipper::Coord_orth::cross(c,d)) *
	 sqrt( c.lengthsq() );

      // b*b * [ -a.c+(a.b)(b.c)/(b*b) ] = -a.c*b*b + (a.b)(b.c)
      double I = - clipper::Coord_orth::dot(b,d)*c.lengthsq()
	 + clipper::Coord_orth::dot(b,c)*clipper::Coord_orth::dot(c,d);

      double psi = clipper::Util::rad2d(atan2(H,I));
      if (psi < 180.0)
	 psi += 360.0;
      if (psi > 180.0)
	 psi -= 360.0;


      if ( clipper::Util::isnan(phi) ) {
	 std::cout << "WARNING: observed torsion phi is a NAN!" << std::endl;
	 // throw an exception
      } 
      if ( clipper::Util::isnan(psi) ) {
	 std::cout << "WARNING: observed torsion psi is a NAN!" << std::endl;
	 // throw an exception
      }
	    
      double phir = clipper::Util::d2rad(phi);
      double psir = clipper::Util::d2rad(psi);
      double R = restraints->rama_prob(phir, psir);
	    
      // std::cout << "df rama distortion for " << phi << " " << psi << " is "
      // << R << std::endl;

      // this can throw an exception
      coot::distortion_torsion_gradients_t dtg_phi =
	 fill_distortion_torsion_gradients(P1, P2, P3, P4);

      // this can throw an exception
      coot::distortion_torsion_gradients_t dtg_psi =
	 fill_distortion_torsion_gradients(P2, P3, P4, P5);

      // Faster to use these, not calculate them above?
      //

      double tan_phir = tan(phir);
      double tan_psir = tan(psir);

      double multiplier_phi = 1.0;
      double multiplier_psi = 1.0;

      if (restraints->rama_type == restraints_container_t::RAMA_TYPE_ZO) {
	 std::pair<float,float> zo_rama_pair = restraints->zo_rama_grad(rama_restraint.rama_plot_residue_type, phir, psir);
	 if (false)
	    std::cout << "debug:: in my_df_rama() rama_plot_residue_type is "
		      << rama_restraint.rama_plot_residue_type << " gradients "
		      << zo_rama_pair.first << " " << zo_rama_pair.second
		      << std::endl;

	 multiplier_phi = -restraints->get_rama_plot_weight()/(1.0 + tan_phir*tan_phir) * zo_rama_pair.first;
	 multiplier_psi = -restraints->get_rama_plot_weight()/(1.0 + tan_psir*tan_psir) * zo_rama_pair.second;
      } else {
	 LogRamachandran::Lgrad lgrd = restraints->rama_grad(phir, psir);
	 multiplier_phi = 10.0/(1.0 + tan_phir*tan_phir) * lgrd.DlogpDphi;
	 multiplier_psi = 10.0/(1.0 + tan_psir*tan_psir) * lgrd.DlogpDpsi;
      }

      double xP1_contrib = multiplier_phi*dtg_phi.dD_dxP1;
      double yP1_contrib = multiplier_phi*dtg_phi.dD_dyP1;
      double zP1_contrib = multiplier_phi*dtg_phi.dD_dzP1;

      double xP2_contrib = multiplier_phi*dtg_phi.dD_dxP2;
      double yP2_contrib = multiplier_phi*dtg_phi.dD_dyP2;
      double zP2_contrib = multiplier_phi*dtg_phi.dD_dzP2;

      double xP3_contrib = multiplier_phi*dtg_phi.dD_dxP3;
      double yP3_contrib = multiplier_phi*dtg_phi.dD_dyP3;
      double zP3_contrib = multiplier_phi*dtg_phi.dD_dzP3;

      double xP4_contrib = multiplier_phi*dtg_phi.dD_dxP4;
      double yP4_contrib = multiplier_phi*dtg_phi.dD_dyP4;
      double zP4_contrib = multiplier_phi*dtg_phi.dD_dzP4;

      // The class variable gives a misleading name here for the
      // follwing blocks. P2 is in postion 1 for dtg_phi, P3 is
      // in position 2, P4 is called in the 3rd position (and
      // P5 in 4th).

      xP2_contrib += multiplier_psi * dtg_psi.dD_dxP1;
      yP2_contrib += multiplier_psi * dtg_psi.dD_dyP1;
      zP2_contrib += multiplier_psi * dtg_psi.dD_dzP1;

      xP3_contrib += multiplier_psi * dtg_psi.dD_dxP2;
      yP3_contrib += multiplier_psi * dtg_psi.dD_dyP2;
      zP3_contrib += multiplier_psi * dtg_psi.dD_dzP2;

      xP4_contrib += multiplier_psi * dtg_psi.dD_dxP3;
      yP4_contrib += multiplier_psi * dtg_psi.dD_dyP3;
      zP4_contrib += multiplier_psi * dtg_psi.dD_dzP3;

      if (0) { 
	 xP2_contrib = 0.0;
	 yP2_contrib = 0.0;
	 zP2_contrib = 0.0;

	 xP3_contrib = 0.0;
	 yP3_contrib = 0.0;
	 zP3_contrib = 0.0;

	 xP4_contrib = 0.0;
	 yP4_contrib = 0.0;
	 zP4_contrib = 0.0;
      }

      double xP5_contrib = multiplier_psi*dtg_psi.dD_dxP4;
      double yP5_contrib = multiplier_psi*dtg_psi.dD_dyP4;
      double zP5_contrib = multiplier_psi*dtg_psi.dD_dzP4;

      if (! rama_restraint.fixed_atom_flags[0]) { 
	 idx = 3*(rama_restraint.atom_index_1);
	 // gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + xP1_contrib);
	 // gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + yP1_contrib);
	 // gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + zP1_contrib);

	 results[idx  ] += xP1_contrib;
	 results[idx+1] += yP1_contrib;
	 results[idx+2] += zP1_contrib;
      }

      if (! rama_restraint.fixed_atom_flags[1]) { 
	 idx = 3*(rama_restraint.atom_index_2);
	 // gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + xP2_contrib);
	 // gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + yP2_contrib);
	 // gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + zP2_contrib);

	 results[idx  ] += xP2_contrib;
	 results[idx+1] += yP2_contrib;
	 results[idx+2] += zP2_contrib;
      }

      if (! rama_restraint.fixed_atom_flags[2]) { 
	 idx = 3*(rama_restraint.atom_index_3);
	 // gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + xP3_contrib);
	 // gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + yP3_contrib);
	 // gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + zP3_contrib);

	 results[idx  ] += xP3_contrib;
	 results[idx+1] += yP3_contrib;
	 results[idx+2] += zP3_contrib;
      }

      if (! rama_restraint.fixed_atom_flags[3]) { 
	 idx = 3*(rama_restraint.atom_index_4);
	 // gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + xP4_contrib);
	 // gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + yP4_contrib);
	 // gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + zP4_contrib);

	 results[idx  ] += xP4_contrib;
	 results[idx+1] += yP4_contrib;
	 results[idx+2] += zP4_contrib;
      }

      if (! rama_restraint.fixed_atom_flags[4]) { 
	 idx = 3*(rama_restraint.atom_index_5);
	 // gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + xP5_contrib);
	 // gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + yP5_contrib);
	 // gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + zP5_contrib);

	 results[idx  ] += xP5_contrib;
	 results[idx+1] += yP5_contrib;
	 results[idx+2] += zP5_contrib;
      }
   }
   catch (const std::runtime_error &rte) {
      std::cout << "ERROR:: my_df_rama() caught " << rte.what() << std::endl;
   } 

}

#endif // threads
#endif // HAVE_GSL
