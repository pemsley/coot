
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

// house-keeping functions for the holders of gradient info for threads

#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

void
coot::restraints_container_t::make_df_restraints_indices() {

   // std::cout << "---------------------------------------------------------------" << std::endl;
   // std::cout << "            make_df_restraints_indices() " << size() << std::endl;
   // std::cout << "---------------------------------------------------------------" << std::endl;

   // does restraints index vectors and df_by_thread_results

   // this can be called more than once. It is called, for example, after
   // extra_restraints are added (geman-mcclure).

   // no longer use n_t : now use n_r_s - the number of restraints-indices sets.

   // I think that the restaints sets were not finishing evenly, so other restraints
   // thread were waiting for one other to finish (typically). Let's have more
   // sets and they all get shoved on to the thread pool queue - hopefully
   // that will address some timing issues.
   // This (or something similar) should probably for the evaluation of
   // distortion too.
   //
   unsigned int n_r_s = n_threads; // needs optimizing
   unsigned int restraints_size = size();

   // Perhaps restraints_container_t was constructed without setting the thread pool.
   // So we get here with n_t == 0.
   // But the refinement still needs somewhere to put the results. Make single vectors
   //
   if (n_r_s == 0) n_r_s = 1;


   // these are now class variables
   // std::vector<std::vector<std::size_t> > restraints_indices(n_t);
   // std::vector<std::vector<double> > df_by_thread_results(n_t);
   restraints_indices.clear();
   restraints_indices.resize(n_r_s);
   if (df_by_thread_results.size() > 0) {
      if (false)
	 std::cout << "currently df_by_thread_results has size " << df_by_thread_results.size()
		   << " - so clearing" << std::endl;
      df_by_thread_results.clear();
   }
   std::vector<mmdb::Link> links;
   df_by_thread_results.resize(n_r_s); // this may not be a good idea, needs optimization.

   // each restraint_index vector will contain about
   // restraints_size/n_threads restraint indices
   //
   int r_reserve_size = std::lround(static_cast<float>(restraints_size)/static_cast<float>(n_r_s)) + 2;

   // First fill the restraints indices vectors
   //
   for (std::size_t n=0; n<restraints_indices.size(); n++) {
      // std::cout << "DEBUG:: reserve size " << r_reserve_size << " for restraints_index set "
      //           << n << std::endl;
      restraints_indices[n].reserve(r_reserve_size);
   }

   unsigned int i_thread = 0; // insert to vector for this thread
   for (unsigned int ir=0; ir<restraints_size; ir++) {
      // std::cout << "pushing back ir " << ir << " to i_thread " << i_thread << " of " << restraints_size << std::endl;
      restraints_indices[i_thread].push_back(ir);
      ++i_thread;
      if (i_thread==n_r_s) i_thread=0;
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
   unsigned int n_var = n_variables();
   for (std::size_t ii=0; ii<n_r_s; ii++)
      df_by_thread_results[ii] = std::vector<double>(n_var, 0);


   // like above, but this time we split the set of _atoms_ into sets for each thread
   //
   // but unlike restraints_vec, this does not dynamically change size.

   df_by_thread_atom_indices.clear();
   df_by_thread_atom_indices.resize(n_r_s);
   i_thread = 0; // not really threads - now index for sets of restraints-indices
   unsigned int n = get_n_atoms();
   for (unsigned int ir=0; ir<n; ir++) {
      // std::cout << "adding atom ir " << ir << " to thread indices vec " << i_thread << " of " << n_t << std::endl;
      df_by_thread_atom_indices[i_thread].push_back(ir);
      ++i_thread;
      if (i_thread==n_r_s) i_thread=0;
   }

   // add this to the class when it works again
   // threaded_distortion_container.clear();
   // threaded_distortion_container.resize(size());
}
#endif // HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

void
coot::restraints_container_t::clear_df_by_thread_results() {

   for (std::size_t i=0; i<df_by_thread_results.size(); i++) {
      std::vector<double> &v = df_by_thread_results[i];
      for (std::size_t j=0; j<v.size(); j++) {
	 v[j] = 0.0;
      }
   }
}
#endif // HAVE_BOOST_BASED_THREAD_POOL_LIBRARY


// parallel version of my_df()
void
coot::split_the_gradients_with_threads(const gsl_vector *v,
				       restraints_container_t *restraints_p,
				       gsl_vector *df) {

#ifdef HAVE_CXX_THREAD
#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
   // --------------------------------------- restraints --------------------------------------

   // Does restraints_p->n_threads depends on having the thread pool?
   //
   // I will try to make it not so.

   //x auto tp_0 = std::chrono::high_resolution_clock::now();

   if (! restraints_p->thread_pool_p)
      return; // no derivatives if we don't have a thread pool. Modern only. BOOST only.

   restraints_p->clear_df_by_thread_results();

   // use restraints in the range of restraints_indices[i}
   // to fill results
   //

   //x auto tp_1 = std::chrono::high_resolution_clock::now();

   std::atomic<unsigned int> done_count_for_threads(0);
   for (std::size_t ii=0; ii<restraints_p->restraints_indices.size(); ii++) {

      restraints_p->thread_pool_p->push(process_dfs_in_range,
					std::cref(restraints_p->restraints_indices[ii]),
					restraints_p, v,
					std::ref(restraints_p->df_by_thread_results[ii]),
					std::ref(done_count_for_threads)
					);
   }
   //x auto tp_2 = std::chrono::high_resolution_clock::now();

   // wait for the threads in the thread pool
   while (done_count_for_threads != restraints_p->restraints_indices.size()) {
      std::this_thread::sleep_for(std::chrono::microseconds(1));
   }

   //x auto tp_3 = std::chrono::high_resolution_clock::now();

   unsigned int n_r_s = restraints_p->restraints_indices.size();
   unsigned int n_variables = restraints_p->n_variables();

   /*
      using threads slows things down by ~250us (baah!)

      unsigned int n_var_split = n_variables/2;
      done_count_for_threads = 0;
      restraints_p->thread_pool_p->push(consolidate_derivatives, n_r_s,           0, n_var_split, restraints_p->df_by_thread_results, df, std::ref(done_count_for_threads));
      restraints_p->thread_pool_p->push(consolidate_derivatives, n_r_s, n_var_split, n_variables, restraints_p->df_by_thread_results, df, std::ref(done_count_for_threads));
      while (done_count_for_threads != 2)
         std::this_thread::sleep_for(std::chrono::microseconds(1));
   */

   // consolidate - ~300us, GM restraints don't slow things down!? How can that be? Cache misses?
   //
   for (std::size_t i_r_s=0; i_r_s<n_r_s; i_r_s++) {
      const std::vector<double> &results_block = restraints_p->df_by_thread_results[i_r_s];
      for (unsigned int i=0; i<n_variables; i++) {
	 if (results_block[i] != 0.0) { // this does speed things up a bit
	    *gsl_vector_ptr(df, i) += results_block[i];
	 }
      }
   }

   //x auto tp_4 = std::chrono::high_resolution_clock::now();

   // --------------------------------------- map --------------------------------------

   if (restraints_p->include_map_terms()) {

      //x auto tp_5 = std::chrono::high_resolution_clock::now();

      // we can manipulate df directly because (unlike restraints) each atom can
      // only touch 3 indices in the df vector - and do so uniquely.
      //

      // reset the done count - previously we used this for the restraints, now
      // we will use it for the electron density score of the atoms
      //
      done_count_for_threads = 0;

      for (std::size_t ii=0; ii<restraints_p->df_by_thread_atom_indices.size(); ii++) {
	 restraints_p->thread_pool_p->push(process_electron_density_dfs_for_atoms,
					   restraints_p->df_by_thread_atom_indices[ii],
					   restraints_p, v, df,
					   std::ref(done_count_for_threads));
      }

      //x auto tp_6 = std::chrono::high_resolution_clock::now();

      // wait for the threads in the thread pool (~20us for threads to complete)
      while (done_count_for_threads != restraints_p->df_by_thread_atom_indices.size()) {
	 std::this_thread::sleep_for(std::chrono::microseconds(1));
      }
      //x auto tp_7 = std::chrono::high_resolution_clock::now();

      //x auto d43 = chrono::duration_cast<chrono::microseconds>(tp_4 - tp_3).count();
      //x std::cout << "timings:: distortion consolidation d43 " << std::setw(5) << d43 << " " << std::endl;

      /*
      auto d10 = chrono::duration_cast<chrono::microseconds>(tp_1 - tp_0).count();
      auto d21 = chrono::duration_cast<chrono::microseconds>(tp_2 - tp_1).count();
      auto d32 = chrono::duration_cast<chrono::microseconds>(tp_3 - tp_2).count();
      auto d43 = chrono::duration_cast<chrono::microseconds>(tp_4 - tp_3).count();
      auto d54 = chrono::duration_cast<chrono::microseconds>(tp_5 - tp_4).count();
      auto d65 = chrono::duration_cast<chrono::microseconds>(tp_6 - tp_5).count();
      auto d76 = chrono::duration_cast<chrono::microseconds>(tp_7 - tp_6).count();
      if (true)
	 std::cout << "timings:: distortion "
		   << "d10 " << std::setw(5) << d10 << " "
		   << "d21 " << std::setw(5) << d21 << " "
		   << "d32 " << std::setw(5) << d32 << " "
		   << "d43 " << std::setw(5) << d43 << " "
		   << "d54 " << std::setw(5) << d54 << " "
		   << "d65 " << std::setw(5) << d65 << " "
		   << "d76 " << std::setw(5) << d76 << " "
		   << "\n";
      */
   }

#endif // HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
#endif // HAVE_CXX_THREAD

}

void
coot::consolidate_derivatives(unsigned int thread_index,
                              unsigned int n_restraints_sets,
                              unsigned int variable_idx_start,
                              unsigned int variable_idx_end,  // stop before this end, e.g. 0, 10
                              const std::vector<std::vector<double> > &df_sets_from, gsl_vector *df,
                              std::atomic<unsigned int> &done_count_for_threads) {

   for (unsigned int i=variable_idx_start; i<variable_idx_end; i++) {
      for (std::size_t i_r_s=0; i_r_s<n_restraints_sets; i_r_s++) {
	 if (df_sets_from[i_r_s][i] != 0.0) { // this test does speed things up (a bit)
	    *gsl_vector_ptr(df, i) += df_sets_from[i_r_s][i];
	 }
      }
   }
   done_count_for_threads++;
}



// fill results
void
coot::process_dfs_in_range(int thread_idx,
			   const std::vector<std::size_t> &restraints_indices,
			   coot::restraints_container_t *restraints_p,
			   const gsl_vector *v,
			   std::vector<double> &results,  // fill results
			   std::atomic<unsigned int> &done_count_for_threads
			   ) {

   unsigned int n_restraints = restraints_p->size(); // signed change

   for (std::size_t i=0; i<restraints_indices.size(); i++) {

      // restraints_vec can change size due to pull atom restraints
      if (restraints_indices[i] >= n_restraints)
	 continue;

      const simple_restraint &rest = restraints_p->at(restraints_indices[i]);

      if (false)
	 std::cout << "process_dfs_in_range() i " << i << " restraint index " << restraints_indices[i]
		   << " " << rest << std::endl;

      if (restraints_p->restraints_usage_flag & coot::GEMAN_MCCLURE_DISTANCE_MASK) {
	 if (rest.restraint_type == coot::GEMAN_MCCLURE_DISTANCE_RESTRAINT) {
	    process_dfs_geman_mcclure_distance(rest, restraints_p->geman_mcclure_alpha, v, results);
	    continue;
	 }
      }

      if (restraints_p->restraints_usage_flag & coot::NON_BONDED_MASK) {
	 if (rest.restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT) {
	    if (! rest.is_H_non_bonded_contact || restraints_p->apply_H_non_bonded_contacts_state()) {
	       if (rest.nbc_function == simple_restraint::LENNARD_JONES) {
		  process_dfs_non_bonded_lennard_jones(rest, restraints_p->lennard_jones_epsilon, v, results);
	       } else {
		  process_dfs_non_bonded(rest, v, results);
	       }
	    }
	    continue;
	 }
      }

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

      if (restraints_p->restraints_usage_flag & coot::TRANS_PEPTIDE_MASK)
	 if (rest.restraint_type == coot::TRANS_PEPTIDE_RESTRAINT)
	    process_dfs_trans_peptide(rest, v, results);

      if (restraints_p->restraints_usage_flag & coot::RAMA_PLOT_MASK)
	 if (rest.restraint_type == coot::RAMACHANDRAN_RESTRAINT)
	    process_dfs_rama(rest, restraints_p, v, results);

      if (restraints_p->restraints_usage_flag & coot::PARALLEL_PLANES_MASK)
	 if (rest.restraint_type == coot::PARALLEL_PLANES_RESTRAINT)
	    process_dfs_parallel_planes(rest, v, results);

      if (rest.restraint_type == coot::TARGET_POS_RESTRAINT)
	 process_dfs_target_position(rest, restraints_p->log_cosh_target_distance_scale_factor, v, results);

   }

   done_count_for_threads++; // atomic operation
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

   double weight = 1.0/(restraint.sigma * restraint.sigma);

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

   double a = sqrt(a_vec.lengthsq());
   double b = sqrt(b_vec.lengthsq());

   // Garib's stabilization
   if (a < 0.01) {
      a = 0.01;
      a_vec = clipper::Coord_orth(0.01, 0.01, 0.01);
   }
   if (b < 0.01) {
      b = 0.01;
      b_vec = clipper::Coord_orth(0.01, 0.01, 0.01);
   }
	    
   double l_over_a_sqd = 1.0/(a*a);
   double l_over_b_sqd = 1.0/(b*b);
   double l_ab         = 1.0/(a*b);

   // for the end atoms: 
   // \frac{\partial \theta}{\partial x_k} =
   //    -\frac{1}{sin\theta} [(x_l-x_k)cos\theta + \frac{x_m-x_l}{ab}]
	 
   double a_dot_b = clipper::Coord_orth::dot(a_vec,b_vec);
   double cos_theta = a_dot_b/(a*b);
   // we need to stabilize cos_theta
   if (cos_theta < -1.0) cos_theta = -1.0;
   if (cos_theta >  1.0) cos_theta =  1.0;
   double theta = acos(cos_theta);

   // we need to stabilize $\theta$ too.
   // theta = theta > 0.001 ? theta : 0.001;
   if (theta < 0.001) theta = 0.001; // it was never -ve.

   double prem = -1.0/sin(theta);

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
   double weight = 1.0/(restraint.sigma * restraint.sigma);
   double ds_dth = 2.0*(theta - target_value)*RADTODEG*RADTODEG;
   double w_ds_dth = weight * ds_dth;

   if (!restraint.fixed_atom_flags[0]) { 
      idx = 3*(restraint.atom_index_1);
      // gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + x_k_contrib*w_ds_dth); 
      // gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + y_k_contrib*w_ds_dth); 
      // gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + z_k_contrib*w_ds_dth);

      if (false)
	 std::cout << "debug angle gradient: " << idx << " "
		   << " theta " << theta << " target_value " << target_value << " "
		   << std::setw(12) << x_k_contrib << " "
		   << std::setw(12) << y_k_contrib << " "
		   << std::setw(12) << z_k_contrib << " "
		   << std::setw(12) << " w_ds_dth " << w_ds_dth << std::endl;

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
		  
	 if (false)
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
   // this calculates plane_info.distortion_score, but we don't need it here.
   // should I skip that?
   //
   coot::plane_distortion_info_t plane_info =
      distortion_score_plane_internal(plane_restraint, v, false);
   int n_plane_atoms = plane_restraint.plane_atom_index.size();

   for (int j=0; j<n_plane_atoms; j++) {
      if (! plane_restraint.fixed_atom_flags[j] ) {
	 idx = 3*plane_restraint.plane_atom_index[j].first;
	 devi_len =
	    plane_info.abcd[0]*gsl_vector_get(v,idx  ) + 
	    plane_info.abcd[1]*gsl_vector_get(v,idx+1) +
	    plane_info.abcd[2]*gsl_vector_get(v,idx+2) -
	    plane_info.abcd[3];

	 double weight = 1.0/(plane_restraint.plane_atom_index[j].second *
			      plane_restraint.plane_atom_index[j].second);

	 clipper::Grad_orth<double> d(2.0 * weight * devi_len * plane_info.abcd[0],
				      2.0 * weight * devi_len * plane_info.abcd[1],
				      2.0 * weight * devi_len * plane_info.abcd[2]);
	 results[idx  ] += d.dx();
	 results[idx+1] += d.dy();
	 results[idx+2] += d.dz();
      }
   }
}

void
coot::process_dfs_parallel_planes(const coot::simple_restraint &ppr,
				  const gsl_vector *v,
				  std::vector<double> &results) { // fill results

   unsigned int first_atoms_size = ppr.plane_atom_index.size();

   // this calculates plane_info.distortion_score, but we don't need it here.
   // should I skip that? No. We don't need to optimize parallel plane restraints
   // at the moment.
   //
   plane_distortion_info_t plane_info =
      distortion_score_2_planes(ppr.plane_atom_index, ppr.atom_index_other_plane, ppr.sigma, v);

   unsigned int n_plane_atoms = ppr.plane_atom_index.size();

   double weight = 0.25 /(ppr.sigma * ppr.sigma); // hack down the weight, c.f. distortion_score_2_planes

   // hack the weight - needs a better fix than this
   // weight *= 0.1;

   for (unsigned int j=0; j<n_plane_atoms; j++) {
      if (! ppr.fixed_atom_flags[j] ) {
	 unsigned int idx = 3*ppr.plane_atom_index[j].first;
	 double devi_len =
	    plane_info.abcd[0]*(gsl_vector_get(v,idx  ) - plane_info.centre_1.x()) +
	    plane_info.abcd[1]*(gsl_vector_get(v,idx+1) - plane_info.centre_1.y()) +
	    plane_info.abcd[2]*(gsl_vector_get(v,idx+2) - plane_info.centre_1.z()) -
	    plane_info.abcd[3];

	 clipper::Grad_orth<double> d(2.0 * weight * devi_len * plane_info.abcd[0],
				      2.0 * weight * devi_len * plane_info.abcd[1],
				      2.0 * weight * devi_len * plane_info.abcd[2]);

	 results[idx  ] += d.dx();
	 results[idx+1] += d.dy();
	 results[idx+2] += d.dz();
      }
   }

   // second plane
   n_plane_atoms = ppr.atom_index_other_plane.size();
   for (unsigned int j=0; j<n_plane_atoms; j++) {
      if (! ppr.fixed_atom_flags_other_plane[j] ) {
	 unsigned int idx = 3*ppr.atom_index_other_plane[j].first;
	 double devi_len =
	    plane_info.abcd[0]*(gsl_vector_get(v,idx  ) - plane_info.centre_2.x()) +
	    plane_info.abcd[1]*(gsl_vector_get(v,idx+1) - plane_info.centre_2.y()) +
	    plane_info.abcd[2]*(gsl_vector_get(v,idx+2) - plane_info.centre_2.z()) -
	    plane_info.abcd[3];

	 clipper::Grad_orth<double> d(2.0 * weight * devi_len * plane_info.abcd[0],
				      2.0 * weight * devi_len * plane_info.abcd[1],
				      2.0 * weight * devi_len * plane_info.abcd[2]);

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

   double sigma = 0.03;
   int idx = 3*(restraint.atom_index_1);

   bool harmonic_restraint = true;

   if (harmonic_restraint) {
      
      double constant_part = 2.0 / (sigma * sigma);

      double dist_x = gsl_vector_get(v, idx)   - restraint.atom_pull_target_pos[0];
      double dist_y = gsl_vector_get(v, idx+1) - restraint.atom_pull_target_pos[1];
      double dist_z = gsl_vector_get(v, idx+2) - restraint.atom_pull_target_pos[2];

      // *gsl_vector_ptr(df, idx  ) += constant_part * dist_x;
      // *gsl_vector_ptr(df, idx+1) += constant_part * dist_y;
      // *gsl_vector_ptr(df, idx+2) += constant_part * dist_z;

      results[idx  ] += constant_part * dist_x;
      results[idx+1] += constant_part * dist_y;
      results[idx+2] += constant_part * dist_z;

   } else {

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
}


void
coot::process_electron_density_dfs_for_atoms(int thread_idx,
					     const std::vector<std::size_t> &atom_indices,
					     const restraints_container_t *restraints_p,
					     const gsl_vector *v, gsl_vector *df,
					     std::atomic<unsigned int> &done_count_for_threads) {

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
   done_count_for_threads++;
}


void
coot::process_dfs_trans_peptide(const coot::simple_restraint &restraint,
					 const gsl_vector *v,
				std::vector<double> &results) {

   int idx;

   // checked:    P1 is CA_1
   //             P2 is C_1
   //             P3 is N_2
   //             P4 is CA_2

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

   // mid-point is a misnomer here
   // I mean "point at the "closest-approach" fraction -
   // which is close to the real mid-point, but not quite.

   // i.e. if closest_approach_fraction_CA_CA is 0.9, that means
   // it is a lot closer to CA_2 than CA_1
   //
   double closest_approach_fraction_CA_CA = 0.5;
   double closest_approach_fraction_C_N   = 0.5;
   double best_closest_approach = 0.055;

   const double &p_CA_CA = closest_approach_fraction_CA_CA; // shorthand aliases
   const double &p_C_N   = closest_approach_fraction_C_N;

   double q_CA_CA = 1.0 - closest_approach_fraction_CA_CA;
   double q_C_N   = 1.0 - closest_approach_fraction_C_N;

   clipper::Coord_orth mid_pt_1 = q_CA_CA * P1 + closest_approach_fraction_CA_CA * P4;
   clipper::Coord_orth mid_pt_2 = q_C_N   * P2 + closest_approach_fraction_C_N * P3;

   double dist_sqrd = (mid_pt_2-mid_pt_1).lengthsq();

   double trans_pep_dist_scale_factor = 4000.0; // needs tweaking
   double weight = trans_pep_dist_scale_factor;

   // d is the distance from the "mid-points" to the expected distance
   // between "midpoints" for an ideal trans-peptide
   double b = sqrt(dist_sqrd);
   double delta = b - best_closest_approach;

   double dS_ddelta = weight * 2.0 * delta;
   double db_da = 0.5 / b;

   // std::cout << "b: " << b << " delta: " << delta << std::endl;

   double constant_part = dS_ddelta * db_da;

   double xP1_contrib = constant_part * 2.0 * q_CA_CA * ( mid_pt_1.x() - mid_pt_2.x());
   double yP1_contrib = constant_part * 2.0 * q_CA_CA * ( mid_pt_1.y() - mid_pt_2.y());
   double zP1_contrib = constant_part * 2.0 * q_CA_CA * ( mid_pt_1.z() - mid_pt_2.z());

   double xP2_contrib = constant_part * 2.0 * q_C_N * ( mid_pt_2.x() - mid_pt_1.x());
   double yP2_contrib = constant_part * 2.0 * q_C_N * ( mid_pt_2.y() - mid_pt_1.y());
   double zP2_contrib = constant_part * 2.0 * q_C_N * ( mid_pt_2.z() - mid_pt_1.z());

   double xP3_contrib = constant_part * 2.0 * p_C_N * ( mid_pt_2.x() - mid_pt_1.x());
   double yP3_contrib = constant_part * 2.0 * p_C_N * ( mid_pt_2.y() - mid_pt_1.y());
   double zP3_contrib = constant_part * 2.0 * p_C_N * ( mid_pt_2.z() - mid_pt_1.z());

   double xP4_contrib = constant_part * 2.0 * p_CA_CA * ( mid_pt_1.x() - mid_pt_2.x());
   double yP4_contrib = constant_part * 2.0 * p_CA_CA * ( mid_pt_1.y() - mid_pt_2.y());
   double zP4_contrib = constant_part * 2.0 * p_CA_CA * ( mid_pt_1.z() - mid_pt_2.z());

   if (false)
      std::cout << "fixed_flags: "
		<< restraint.fixed_atom_flags[0] << " "
		<< restraint.fixed_atom_flags[1] << " "
		<< restraint.fixed_atom_flags[2] << " "
		<< restraint.fixed_atom_flags[3] << " "
		<< std::endl;

   if (! restraint.fixed_atom_flags[0]) {
      idx = 3*(restraint.atom_index_1);
      results[idx  ] += xP1_contrib;
      results[idx+1] += yP1_contrib;
      results[idx+2] += zP1_contrib;
   }

   if (! restraint.fixed_atom_flags[1]) {
      idx = 3*(restraint.atom_index_2);
      results[idx  ] += xP2_contrib;
      results[idx+1] += yP2_contrib;
      results[idx+2] += zP2_contrib;
   }

   if (! restraint.fixed_atom_flags[2]) {
      idx = 3*(restraint.atom_index_3);
      results[idx  ] += xP3_contrib;
      results[idx+1] += yP3_contrib;
      results[idx+2] += zP3_contrib;
   }

   if (! restraint.fixed_atom_flags[3]) {
      idx = 3*(restraint.atom_index_4);
      results[idx  ] += xP4_contrib;
      results[idx+1] += yP4_contrib;
      results[idx+2] += zP4_contrib;
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

   double beta  = 1.0 + alpha * z * z;
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
