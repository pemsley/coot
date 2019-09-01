/* ideal/simple-restraint.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 by The University of York
 * Copyright 2008, 2009, 2010  by The University of Oxford
 * Copyright 2016 by Medical Research Council
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

// #define ANALYSE_REFINEMENT_TIMING

// we don't want to compile anything if we don't have gsl
#ifdef HAVE_GSL

#include <string.h> // for strcmp

#ifdef ANALYSE_REFINEMENT_TIMING
#include <sys/time.h> // for gettimeofday()
#endif // ANALYSE_REFINEMENT_TIMING


#ifdef HAVE_CXX_THREAD
#include <thread>
#include <chrono>
#endif

#include <fstream>
#include <algorithm> // for sort
#include <stdexcept>

#include "simple-restraint.hh"

void 
coot::numerical_gradients(gsl_vector *v, 
			  void *params, 
			  gsl_vector *df) {

   // Not that (at present) user-only fixed atoms re removed from numerical gradient calculations,
   // Other fixed atom (in flanking residues) continue to have numerical gradients calculated for them.
   // This is confusing and undesirable.
   // What are the flanking residues?

//    clipper::Coord_orth a(0,0,2); 
//    cout << "length test: " << a.lengthsq() << endl; 

   // double S = coot::distortion_score(v, params); 
   double tmp; 
   double new_S_plus; 
   double new_S_minu; 
   double val;
   double micro_step = 0.0001;  // the difference between the gradients
			        // seems not to depend on the
			        // micro_step size (0.0001 vs 0.001)
                                // ... Hmm... is does for rama restraints?

   if (false) {
      std::cout << "in numerical_gradients: df->size is " << df->size << std::endl;
      std::cout << "in numerical_gradients: v ->size is " <<  v->size << std::endl;
   }

   coot::restraints_container_t *restraints = (coot::restraints_container_t *)params;

   std::vector<double> analytical_derivs(v->size);
   std::vector<double>  numerical_derivs(v->size);

   for (unsigned int i=0; i<df->size; i++)
      analytical_derivs[i] = gsl_vector_get(df, i);

   if (false) {
      std::cout << "debug:: in numerical_gradients() here are the " << restraints->fixed_atom_indices.size()
		<< " fixed_atom indices: \n";
      for (std::size_t ii=0; ii<restraints->fixed_atom_indices.size(); ii++)
	 std::cout << " " << restraints->fixed_atom_indices[ii];
      std::cout << "\n";
   }

   for (unsigned int i=0; i<v->size; i++) { 

      int iat = i/3; // if iat is a fixed atom, then we shouldn't generate val
      bool make_val = true;
      std::vector<int>::const_iterator it;
      if (std::find(restraints->fixed_atom_indices.begin(),
		    restraints->fixed_atom_indices.end(), iat) != restraints->fixed_atom_indices.end())
	 make_val = false;

      if (make_val) {
	 tmp = gsl_vector_get(v, i);
	 gsl_vector_set(v, i, tmp+micro_step);
	 new_S_plus = coot::distortion_score(v, params);
	 gsl_vector_set(v, i, tmp-micro_step);
	 new_S_minu = coot::distortion_score(v, params);
	 // new_S_minu = 2*tmp - new_S_plus;

	 // now put v[i] back to tmp
	 gsl_vector_set(v, i, tmp);

	 val = (new_S_plus - new_S_minu)/(2*micro_step);
      } else {
	 val = 0; // does the constructor of numerical_derivs make this unnecessary?
      }
      numerical_derivs[i] = val;

      // overwrite the analytical gradients with numerical ones:
      // gsl_vector_set(df, i, val);
   }

   for (unsigned int i=0; i<v->size; i++) {
      std::cout << i << " analytical: " << analytical_derivs[i]
		<< " numerical: " << numerical_derivs[i] << "\n";
   }
   
} // Note to self: try 5 atoms and doctor the .rst file if necessary.
  // Comment out the bond gradients.
  // Try getting a perfect structure (from refmac idealise) and
  // distorting an atom by adding 0.5 to one of its coordinates.  We
  // should be able to trace the gradients that way.


void coot::my_df(const gsl_vector *v, 
		 void *params, 
		 gsl_vector *df) {

#ifdef ANALYSE_REFINEMENT_TIMING
#endif // ANALYSE_REFINEMENT_TIMING

   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params;
   int n_var = restraints->n_variables();

   // first initialize the derivative vector:
   for (int i=0; i<n_var; i++) {
      gsl_vector_set(df,i,0);
   }

   // std::cout << "debug:: in my_df() usage_flags " << restraints->restraints_usage_flag
   // << std::endl;
   my_df_bonds     (v, params, df); 
   my_df_angles    (v, params, df);
   my_df_torsions  (v, params, df);
   my_df_trans_peptides(v, params, df);
   my_df_rama      (v, params, df);
   my_df_planes    (v, params, df);
   my_df_non_bonded(v, params, df);
   my_df_chiral_vol(v, params, df);
   my_df_start_pos (v, params, df);
   my_df_parallel_planes(v, params, df);
   my_df_geman_mcclure_distances(v, params, df);
   
   if (restraints->include_map_terms()) {
      // std::cout << "Using map terms " << std::endl;
      coot::my_df_electron_density(v, params, df);
   } 

   if (restraints->do_numerical_gradients_status())
      coot::numerical_gradients((gsl_vector *)v, params, df); 

#ifdef ANALYSE_REFINEMENT_TIMING
#endif // ANALYSE_REFINEMENT_TIMING

}
   
/* The gradients of f, df = (df/dx(k), df/dy(k) .. df/dx(l) .. ). */
void coot::my_df_bonds(const gsl_vector *v,
		       void *params,
		       gsl_vector *df) {

   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params; 

   // the length of gsl_vector should be equal to n_var: 
   // 
   // int n_var = restraints->n_variables();
   //    float derivative_value; 
   int idx; 
   int n_bond_restr = 0; // debugging counter

   //     for (int i=0; i<n_var; i++) { 
   //       gsl_vector_set(df, i, derivative_value); 
   //     } 

    // Now run over the bonds
    // and add the contribution from this bond/restraint to 
    // dS/dx_k dS/dy_k dS/dz_k dS/dx_l dS/dy_l dS/dz_l for each bond
    // 

   if (restraints->restraints_usage_flag & coot::BONDS_MASK) {

      double target_val;
      double b_i_sqrd;
      double weight;
      double x_k_contrib;
      double y_k_contrib;
      double z_k_contrib;

      double x_l_contrib;
      double y_l_contrib;
      double z_l_contrib;

      for (unsigned int i=restraints->restraints_limits_bonds.first; i<=restraints->restraints_limits_bonds.second; i++) {

	 if ( (*restraints)[i].restraint_type == coot::BOND_RESTRAINT) { 

// 	    std::cout << "DEBUG bond restraint fixed flags: "
// 		      << (*restraints)[i].fixed_atom_flags[0] << " "
// 		      << (*restraints)[i].fixed_atom_flags[1] << " "
// 		      << restraints->get_atom((*restraints)[i].atom_index_1)->GetSeqNum() << " "
// 		      << restraints->get_atom((*restraints)[i].atom_index_1)->name
// 		      << " to " 
// 		      << restraints->get_atom((*restraints)[i].atom_index_2)->GetSeqNum() << " "
// 		      << restraints->get_atom((*restraints)[i].atom_index_2)->name
// 		      << std::endl;
	    
// 	    int n_fixed=0;
// 	    if  ((*restraints)[i].fixed_atom_flags[0])
// 	       n_fixed++;
// 	    if  ((*restraints)[i].fixed_atom_flags[1])
// 	       n_fixed++;
	    
	    n_bond_restr++; 

	    target_val = (*restraints)[i].target_value;

	    // what is the index of x_k?
	    idx = 3*((*restraints)[i].atom_index_1); 
	    clipper::Coord_orth a1(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));
	    idx = 3*((*restraints)[i].atom_index_2); 
	    clipper::Coord_orth a2(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));

	    // what is b_i?
	    b_i_sqrd = (a1-a2).lengthsq();
	    b_i_sqrd = b_i_sqrd > 0.01 ? b_i_sqrd : 0.01;  // Garib's stabilization

	    // b_i = clipper::Coord_orth::length(a1,a2); 
	    // b_i = b_i > 0.1 ? b_i : 0.1;  // Garib's stabilization

	    weight = 1/((*restraints)[i].sigma * (*restraints)[i].sigma);
	    
	    // weight = 1.0;
	    // weight = pow(0.021, -2.0);
	    // std::cout << "df weight is " << weight << std::endl;

	    // double constant_part = 2.0*weight*(b_i - target_val)/b_i;
	    double constant_part = 2.0*weight * (1.0 - target_val * f_inv_fsqrt(b_i_sqrd));
	    
	    x_k_contrib = constant_part*(a1.x()-a2.x());
	    y_k_contrib = constant_part*(a1.y()-a2.y());
	    z_k_contrib = constant_part*(a1.z()-a2.z());
	    
	    x_l_contrib = constant_part*(a2.x()-a1.x());
	    y_l_contrib = constant_part*(a2.y()-a1.y());
	    z_l_contrib = constant_part*(a2.z()-a1.z());

	    if (!restraints->at(i).fixed_atom_flags[0]) {
	       idx = 3*restraints->at(i).atom_index_1;

	       *gsl_vector_ptr(df, idx  ) += x_k_contrib;
	       *gsl_vector_ptr(df, idx+1) += y_k_contrib;
	       *gsl_vector_ptr(df, idx+2) += z_k_contrib;

	    } else {
	       // debug
	       if (0) {
		  idx = 3*((*restraints)[i].atom_index_1 - 0);  
		  std::cout << "BOND Fixed atom[0] "
			    << restraints->get_atom((*restraints)[i].atom_index_1)->GetSeqNum() << " " 
			    << restraints->get_atom((*restraints)[i].atom_index_1)->name << " " 
			    << ", Not adding " << x_k_contrib << " "
			    << y_k_contrib << " "
			    << z_k_contrib << " to " << gsl_vector_get(df, idx) << " "
			    << gsl_vector_get(df, idx+1) << " "
			    << gsl_vector_get(df, idx+2) << std::endl;
	       }
	    }

	    if (!(*restraints)[i].fixed_atom_flags[1]) { 
	       idx = 3*((*restraints)[i].atom_index_2 - 0); 
	       // std::cout << "bond second non-fixed  idx is " << idx << std::endl; 
	       // cout << "second idx is " << idx << endl;

	       *gsl_vector_ptr(df, idx  ) += x_l_contrib;
	       *gsl_vector_ptr(df, idx+1) += y_l_contrib;
	       *gsl_vector_ptr(df, idx+2) += z_l_contrib;

	    } else {
	       // debug
	       if (false) {
		  idx = 3*((*restraints)[i].atom_index_2 - 0);  
		  std::cout << "BOND Fixed atom[1] "
			    << restraints->get_atom((*restraints)[i].atom_index_2)->GetSeqNum() << " " 
			    << restraints->get_atom((*restraints)[i].atom_index_2)->name << " " 
			    << ", Not adding " << x_k_contrib << " "
			    << y_k_contrib << " "
			    << z_k_contrib << " to "
			    << gsl_vector_get(df, idx) << " "
			    << gsl_vector_get(df, idx+1) << " "
			    << gsl_vector_get(df, idx+2) << std::endl;
	       }
	    } 
	 }
      }
   }
}

// the calling function makes sure that this_restraint is a non-bonded contact restraint
//
void
coot::my_df_non_bonded_single(const gsl_vector *v,
			      gsl_vector *df,
			      const simple_restraint &this_restraint) {

   // pass const restraints_container_t &restraints to debug atoms

   const double &target_val = this_restraint.target_value;

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

   double b_i_sqrd;
   double weight;
   double x_k_contrib;
   double y_k_contrib;
   double z_k_contrib;

   double x_l_contrib;
   double y_l_contrib;
   double z_l_contrib;

   b_i_sqrd = (a1-a2).lengthsq();

   if (b_i_sqrd < this_restraint.target_value * this_restraint.target_value) {

      weight = 1.0/(this_restraint.sigma * this_restraint.sigma);

      double b_i = sqrt(b_i_sqrd);
      // b_i_sqrd = b_i_sqrd > 0.01 ? b_i_sqrd : 0.01;  // Garib's stabilization
      double constant_part = 2.0*weight * (1 - target_val * f_inv_fsqrt(b_i_sqrd));

      x_k_contrib = constant_part*(a1.x()-a2.x());
      y_k_contrib = constant_part*(a1.y()-a2.y());
      z_k_contrib = constant_part*(a1.z()-a2.z());

      x_l_contrib = constant_part*(a2.x()-a1.x());
      y_l_contrib = constant_part*(a2.y()-a1.y());
      z_l_contrib = constant_part*(a2.z()-a1.z());

      if (! this_restraint.fixed_atom_flags[0]) {
	 idx_1 = 3*this_restraint.atom_index_1;

// 	 std::cout << "df: " << this_restraint.atom_index_1 << " " << this_restraint.atom_index_2 << " "
// 		   << x_k_contrib << " " << y_k_contrib << " " << z_k_contrib << " "
// 		   << restraints.get_atom_spec(this_restraint.atom_index_1) << " " 
// 		   << restraints.get_atom_spec(this_restraint.atom_index_2) << std::endl;

 	 *gsl_vector_ptr(df, idx_1  ) += x_k_contrib;
 	 *gsl_vector_ptr(df, idx_1+1) += y_k_contrib;
 	 *gsl_vector_ptr(df, idx_1+2) += z_k_contrib;
      }

      if (! this_restraint.fixed_atom_flags[1]) {
	 idx_2 = 3*this_restraint.atom_index_2;

// 	 std::cout << "df: " << this_restraint.atom_index_2 << " " << this_restraint.atom_index_1 << " "
// 		   << x_l_contrib << " " << y_l_contrib << " " << z_l_contrib << std::endl;

 	 *gsl_vector_ptr(df, idx_2  ) += x_l_contrib;
 	 *gsl_vector_ptr(df, idx_2+1) += y_l_contrib;
 	 *gsl_vector_ptr(df, idx_2+2) += z_l_contrib;
      }
   }
}

#ifdef HAVE_CXX_THREAD
void
coot::my_df_non_bonded_thread_dispatcher(int thread_idx,
					 const gsl_vector *v,
					 gsl_vector *df,
					 restraints_container_t *restraints_p,
					 int idx_start,
					 int idx_end,
					 std::atomic<unsigned int> &done_count_for_threads) {

   // std::cout << "my_df_non_bonded_thread_dispatcher() start " << idx_start << " " << idx_end
   // << std::endl;
   for (int i=idx_start; i<idx_end; i++) {
      // std::cout << "dispatching i " << i << " in range " << idx_start << " to " << idx_end
      // << std::endl;
      const simple_restraint &this_restraint = (*restraints_p)[i];
      if (this_restraint.restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT)
	 my_df_non_bonded_single(v, df, this_restraint);
   }
   done_count_for_threads++;
}
#endif // HAVE_CXX_THREAD

void
coot::my_df_non_bonded(const  gsl_vector *v, 
			void *params, 
			gsl_vector *df) {

   // first extract the object from params 
   //
   restraints_container_t *restraints_p = static_cast<restraints_container_t *>(params);

   // the length of gsl_vector should be equal to n_var: 
   // 
   // int n_var = restraints->n_variables();
   // float derivative_value; 
   int idx; 
   int n_non_bonded_restr = 0; // debugging counter

   if (restraints_p->restraints_usage_flag & coot::NON_BONDED_MASK) { 

      unsigned int restraints_size = restraints_p->size();

#ifdef HAVE_CXX_THREAD

      std::atomic<unsigned int> done_count_for_threads(0); // updated by my_df_non_bonded_thread_dispatcher

      if ((restraints_p->thread_pool_p) && (restraints_p->n_threads > 0)) {

	 unsigned int n_per_thread = restraints_size/restraints_p->n_threads;

	 for (unsigned int i_thread=0; i_thread<restraints_p->n_threads; i_thread++) {
	    int idx_start = i_thread * n_per_thread;
	    int idx_end   = idx_start + n_per_thread;
	    // for the last thread, set the end atom index
	    if (i_thread == (restraints_p->n_threads - 1))
	       idx_end = restraints_size; // for loop uses iat_start and tests for < iat_end

	    restraints_p->thread_pool_p->push(my_df_non_bonded_thread_dispatcher,
					      v, df, restraints_p, idx_start, idx_end,
					      std::ref(done_count_for_threads));

	 }
	 // restraints->thread_pool_p->stop(true); // wait
	 while (done_count_for_threads != restraints_p->n_threads) {
	    std::this_thread::sleep_for(std::chrono::microseconds(1));
	 }
	 
      } else {
	 for (unsigned int i=0; i<restraints_size; i++) {
	    const simple_restraint &this_restraint = (*restraints_p)[i];
	    if (this_restraint.restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT) {
	       if (this_restraint.fixed_atom_flags[0]==false || this_restraint.fixed_atom_flags[1]==false)
		  my_df_non_bonded_single(v, df, this_restraint);
	    }
	 }
      }

#else
      for (unsigned int i=0; i<restraints_size; i++) {
	 const simple_restraint &this_restraint = (*restraints_p)[i];
	 if (this_restraint.restraint_type == coot::NON_BONDED_CONTACT_RESTRAINT) {
	    // no need to calculate anything if both these atoms are non-moving
	    //
	    if (this_restraint.fixed_atom_flags[0]==false || this_restraint.fixed_atom_flags[1]==false)
	       my_df_non_bonded_single(v, df, this_restraint);
	 }
      }
#endif

   }
}

void
coot::my_df_geman_mcclure_distances(const  gsl_vector *v, 
				    void *params, 
				    gsl_vector *df) {

   restraints_container_t *restraints = (restraints_container_t *)params;

   // the length of gsl_vector should be equal to n_var: 
   // 
   // int n_var = restraints->n_variables();
   // float derivative_value; 
   int idx;

   if (restraints->restraints_usage_flag & GEMAN_MCCLURE_DISTANCE_MASK) { 

      double target_val;
      double b_i_sqrd;
      double x_k_contrib;
      double y_k_contrib;
      double z_k_contrib;
      
      double x_l_contrib;
      double y_l_contrib;
      double z_l_contrib;

      if (false)
	 std::cout << "in my_df_geman_mcclure_distances restraints limits "
		   << restraints->restraints_limits_geman_mclure.first  << " "
		   << restraints->restraints_limits_geman_mclure.second << std::endl;

      for (unsigned int i=restraints->restraints_limits_geman_mclure.first; i<=restraints->restraints_limits_geman_mclure.second; i++) {

	 const simple_restraint &rest = (*restraints)[i];
      
	 if (rest.restraint_type == GEMAN_MCCLURE_DISTANCE_RESTRAINT) { 

	    idx = 3*rest.atom_index_1;
	    clipper::Coord_orth a1(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));
	    idx = 3*rest.atom_index_2;
	    clipper::Coord_orth a2(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));

	    // what is b_i?
	    // b_i = clipper::Coord_orth::length(a1,a2);
	    b_i_sqrd = (a1-a2).lengthsq();
	    b_i_sqrd = b_i_sqrd > 0.01 ? b_i_sqrd : 0.01;  // Garib's stabilization

	    if (true) {

	       double b_i = sqrt(b_i_sqrd);
	       double weight = 1.0/(rest.sigma * rest.sigma);

	       // double constant_part = 2.0*weight*(b_i - target_val)/b_i;
	       // double constant_part = 2.0*weight * (1 - target_val * f_inv_fsqrt(b_i_sqrd));

	       // Let z = (boi - bi)/sigma
	       //    S_i = z^2/(1 + alpha * z^2)
	       //
	       double bit = b_i - rest.target_value;
	       double z = bit/rest.sigma;

	       const double &alpha = restraints->geman_mcclure_alpha;
	       double beta  = 1 + alpha * z * z;
	       double d_Si_d_zi = 2.0 * z  / (beta * beta);
	       double d_zi_d_bi = 1.0/rest.sigma;
	       double d_b_d_x_m = 1.0/b_i;

	       double constant_part_gm = d_Si_d_zi * d_zi_d_bi * d_b_d_x_m;
	       double constant_part = constant_part_gm;

	       {
		  const double &target_val = rest.target_value;
		  double constant_part_lsq = 2.0*weight * (1 - target_val * f_inv_fsqrt(b_i_sqrd));

		  if (false)
		     std::cout << "debug compare idx " << i << " add lsq "
			       << constant_part_lsq << " GM: " << constant_part_gm
			       << " from beta " << beta
			       << " d_Si_d_zi " << d_Si_d_zi
			       << " d_zi_d_bi " << d_zi_d_bi
			       << " d_b_d_x_m " << d_b_d_x_m
			       << " b_i " << b_i
			       << " target " << rest.target_value
			       << " sigma " << rest.sigma
			       << " z " << z
			       << std::endl;

		  constant_part = constant_part_lsq / (beta * beta);

		  // constant_part = constant_part_lsq; // force least squares
	       }


	       // The final part is dependent on the coordinates:
 	       x_k_contrib = constant_part*(a1.x()-a2.x());
	       y_k_contrib = constant_part*(a1.y()-a2.y());
 	       z_k_contrib = constant_part*(a1.z()-a2.z());
 	       x_l_contrib = constant_part*(a2.x()-a1.x());
 	       y_l_contrib = constant_part*(a2.y()-a1.y());
 	       z_l_contrib = constant_part*(a2.z()-a1.z());

	       if (! rest.fixed_atom_flags[0]) {
#ifdef HAVE_CXX_THREAD
		  // use atomic lock to access derivs of atom atom_idx_1
		  unsigned int unlocked = 0;
		  while (! restraints->gsl_vector_atom_pos_deriv_locks.get()[rest.atom_index_1].compare_exchange_weak(unlocked, 1)) {
		     std::cout << "oops locked! [0] " << rest.atom_index_1 << std::endl;
		     std::this_thread::sleep_for(std::chrono::nanoseconds(10));
		     unlocked = 0;
		  }
#endif
		  idx = 3*rest.atom_index_1;
		  *gsl_vector_ptr(df, idx  ) += x_k_contrib;
		  *gsl_vector_ptr(df, idx+1) += y_k_contrib;
		  *gsl_vector_ptr(df, idx+2) += z_k_contrib;
		  // std::cout << "unlock [0] " << rest.atom_index_1 << std::endl;
#ifdef HAVE_CXX_THREAD
		  restraints->gsl_vector_atom_pos_deriv_locks.get()[rest.atom_index_1] = 0; // unlock
#endif
	       }

	       if (! rest.fixed_atom_flags[1]) { 
#ifdef HAVE_CXX_THREAD
		  // use atomic lock to access derivs of atom atom_idx_2
		  unsigned int unlocked = 0;
		  while (! restraints->gsl_vector_atom_pos_deriv_locks.get()[rest.atom_index_2].compare_exchange_weak(unlocked, 1)) {
		     std::cout << "oops locked! [1] " << rest.atom_index_2 << std::endl;
		     std::this_thread::sleep_for(std::chrono::nanoseconds(10));
		     unlocked = 0;
		  }
#endif
		  idx = 3*rest.atom_index_2;
		  *gsl_vector_ptr(df, idx  ) += x_l_contrib;
		  *gsl_vector_ptr(df, idx+1) += y_l_contrib;
		  *gsl_vector_ptr(df, idx+2) += z_l_contrib;
		  // std::cout << "unlock [1] " << rest.atom_index_1 << std::endl;
#ifdef HAVE_CXX_THREAD
		  restraints->gsl_vector_atom_pos_deriv_locks.get()[rest.atom_index_2] = 0; // unlock
#endif
	       }
	    }
	 }
      }
   }
}
   


// Add in the angle gradients
//
void coot::my_df_angles(const gsl_vector *v, 
			void *params, 
			gsl_vector *df) {

   int n_angle_restr = 0; 
   int idx; 

   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params;

   if (restraints->restraints_usage_flag & coot::ANGLES_MASK) { // 2: angles

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

      for (unsigned int i=restraints->restraints_limits_angles.first; i<=restraints->restraints_limits_angles.second; i++) {

	 if ( (*restraints)[i].restraint_type == coot::ANGLE_RESTRAINT) {

	    n_angle_restr++;

	    double target_value = (*restraints)[i].target_value*DEGTORAD;

	    idx = 3*((*restraints)[i].atom_index_1); 
	    clipper::Coord_orth k(gsl_vector_get(v,idx), 
			gsl_vector_get(v,idx+1), 
			gsl_vector_get(v,idx+2));
	    idx = 3*((*restraints)[i].atom_index_2); 
	    clipper::Coord_orth l(gsl_vector_get(v,idx), 
			gsl_vector_get(v,idx+1), 
			gsl_vector_get(v,idx+2));
	    idx = 3*((*restraints)[i].atom_index_3); 
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
	    weight = 1/((*restraints)[i].sigma * (*restraints)[i].sigma);
	    ds_dth = 2*(theta - target_value)*RADTODEG*RADTODEG;
	    w_ds_dth = weight * ds_dth; 

	    if (!(*restraints)[i].fixed_atom_flags[0]) { 
	       idx = 3*((*restraints)[i].atom_index_1);
	       gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + x_k_contrib*w_ds_dth); 
	       gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + y_k_contrib*w_ds_dth); 
	       gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + z_k_contrib*w_ds_dth); 
	    }
	    if (!(*restraints)[i].fixed_atom_flags[2]) { 
	       idx = 3*((*restraints)[i].atom_index_3);
	       gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + x_m_contrib*w_ds_dth); 
	       gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + y_m_contrib*w_ds_dth); 
	       gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + z_m_contrib*w_ds_dth); 
	    }

	    // and mid atom
	    if (!(*restraints)[i].fixed_atom_flags[1]) { 
	       idx = 3*((*restraints)[i].atom_index_2);
	       gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + x_l_mid_contrib*w_ds_dth); 
	       gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + y_l_mid_contrib*w_ds_dth); 
	       gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + z_l_mid_contrib*w_ds_dth); 
	    }
	 }
      }
   }
   // cout << "added " << n_angle_restr << " angle restraint gradients" << endl; 
} 

// Add in the torsion gradients
//
void coot::my_df_torsions(const gsl_vector *v, 
			  void *params, 
			  gsl_vector *df) {

   my_df_torsions_internal(v, params, df, 0);
}


// this can throw a std::runtime_error if there is a problem calculating the torsion.
// 
coot::distortion_torsion_gradients_t
coot::fill_distortion_torsion_gradients(const clipper::Coord_orth &P1,
					const clipper::Coord_orth &P2,
					const clipper::Coord_orth &P3,
					const clipper::Coord_orth &P4) {

   coot::distortion_torsion_gradients_t dtg; 
   clipper::Coord_orth a = P2 - P1;
   clipper::Coord_orth b = P3 - P2; 
   clipper::Coord_orth c = P4 - P3;

   double b_lengthsq = b.lengthsq();
   double b_length = sqrt(b_lengthsq); 
   if (b_length < 0.01) { 
      b_length = 0.01; // Garib's stabilization
      b_lengthsq = b_length * b_length; 
   }

   if (true) {
      if (b_length < 0.5)
	 std::cout << "ERROR:: fill_distortion_torsion_gradients() problem with b_length "
		   << b_length << std::endl;
   }

   double H = -clipper::Coord_orth::dot(a,c);
   double J =  clipper::Coord_orth::dot(a,b); 
   double K =  clipper::Coord_orth::dot(b,c); 
   double L = 1/b_lengthsq;
   double one_over_b = 1/b_length;

   // a.(bxc)/b
   double E = one_over_b*clipper::Coord_orth::dot(a,clipper::Coord_orth::cross(b,c)); 
   // -a.c+(a.b)(b.c)/(b*b)
   double G = H+J*K*L;
   double F = 1/G;
   if (G == 0.0) F = 999999999.9;

   dtg.theta = clipper::Util::rad2d(atan2(E,G));
   if ( clipper::Util::isnan(dtg.theta) ) {
      std::cout << "oops: bad torsion: " << E << "/" << G << std::endl;
      std::string mess = "WARNING: fill_distortion_torsion_gradients() observed torsion theta is a NAN!";
      throw std::runtime_error(mess);
   }

   double al = sqrt(clipper::Coord_orth::dot(a,a));
   double cl = sqrt(clipper::Coord_orth::dot(c,c));
   double cos_a1 = clipper::Coord_orth::dot(a,b)/(al*b_length);
   double cos_a2 = clipper::Coord_orth::dot(b,c)/(b_length*cl);

   if (false)
      std::cout << "F " << F << " G " << G << " E " << E << " theta " << dtg.theta
		<< " cos(a1) " << cos_a1 << " cos(a2) " << cos_a2
		<< std::endl;

   // instabilty when the P2-P3-P4 or P1-P2-p3 angle is linear. Give up with the derivatives
   // similar escape in the distortion score
   
   if (cos_a1 > 0.9 || cos_a2> 0.9) {

      dtg.zero_gradients = true;

      // x
      dtg.dD_dxP1 = 0;
      dtg.dD_dxP2 = 0;
      dtg.dD_dxP3 = 0;
      dtg.dD_dxP4 = 0;

      // y
      dtg.dD_dyP1 = 0;
      dtg.dD_dyP2 = 0;
      dtg.dD_dyP3 = 0;
      dtg.dD_dyP4 = 0;

      // z
      dtg.dD_dzP1 = 0;
      dtg.dD_dzP2 = 0;
      dtg.dD_dzP3 = 0;
      dtg.dD_dzP4 = 0;

      return dtg;
   }


   // 	    double clipper_theta = 
   // 	       clipper::Util::rad2d(clipper::Coord_orth::torsion(P1, P2, P3, P4));

   dtg.zero_gradients = false;

   // x
   double dH_dxP1 =  c.x(); 
   double dH_dxP2 = -c.x(); 
   double dH_dxP3 =  a.x(); 
   double dH_dxP4 = -a.x(); 

   double dK_dxP1 = 0; 
   double dK_dxP2 = -c.x(); 
   double dK_dxP3 =  c.x() - b.x(); 
   double dK_dxP4 =  b.x(); 

   double dJ_dxP1 = -b.x();
   double dJ_dxP2 =  b.x() - a.x();
   double dJ_dxP3 =  a.x();
   double dJ_dxP4 =  0;

   double dL_dxP1 = 0;
   double dL_dxP2 =  2*(P3.x()-P2.x())*L*L; // check sign
   double dL_dxP3 = -2*(P3.x()-P2.x())*L*L;
   double dL_dxP4 = 0;

   // y 
   double dH_dyP1 =  c.y(); 
   double dH_dyP2 = -c.y(); 
   double dH_dyP3 =  a.y(); 
   double dH_dyP4 = -a.y(); 

   double dK_dyP1 = 0; 
   double dK_dyP2 = -c.y(); 
   double dK_dyP3 =  c.y() - b.y(); 
   double dK_dyP4 =  b.y(); 

   double dJ_dyP1 = -b.y();
   double dJ_dyP2 =  b.y() - a.y();
   double dJ_dyP3 =  a.y();
   double dJ_dyP4 =  0;

   double dL_dyP1 = 0;
   double dL_dyP2 =  2*(P3.y()-P2.y())*L*L; // check sign
   double dL_dyP3 = -2*(P3.y()-P2.y())*L*L;
   double dL_dyP4 = 0;

   // z 
   double dH_dzP1 =  c.z(); 
   double dH_dzP2 = -c.z(); 
   double dH_dzP3 =  a.z(); 
   double dH_dzP4 = -a.z(); 

   double dK_dzP1 = 0; 
   double dK_dzP2 = -c.z(); 
   double dK_dzP3 =  c.z() - b.z(); 
   double dK_dzP4 =  b.z(); 

   double dJ_dzP1 = -b.z();
   double dJ_dzP2 =  b.z() - a.z();
   double dJ_dzP3 =  a.z();
   double dJ_dzP4 =  0;

   double dL_dzP1 = 0;
   double dL_dzP2 =  2*(P3.z()-P2.z())*L*L;
   double dL_dzP3 = -2*(P3.z()-P2.z())*L*L;
   double dL_dzP4 = 0;

   // M
   double dM_dxP1 = -(b.y()*c.z() - b.z()*c.y());
   double dM_dxP2 =  (b.y()*c.z() - b.z()*c.y()) + (a.y()*c.z() - a.z()*c.y());
   double dM_dxP3 =  (b.y()*a.z() - b.z()*a.y()) - (a.y()*c.z() - a.z()*c.y());
   double dM_dxP4 = -(b.y()*a.z() - b.z()*a.y());

   double dM_dyP1 = -(b.z()*c.x() - b.x()*c.z());
   double dM_dyP2 =  (b.z()*c.x() - b.x()*c.z()) + (a.z()*c.x() - a.x()*c.z());
   double dM_dyP3 = -(a.z()*c.x() - a.x()*c.z()) + (b.z()*a.x() - b.x()*a.z());
   double dM_dyP4 = -(b.z()*a.x() - b.x()*a.z());

   double dM_dzP1 = -(b.x()*c.y() - b.y()*c.x());
   double dM_dzP2 =  (b.x()*c.y() - b.y()*c.x()) + (a.x()*c.y() - a.y()*c.x());
   double dM_dzP3 = -(a.x()*c.y() - a.y()*c.x()) + (a.y()*b.x() - a.x()*b.y());
   double dM_dzP4 = -(a.y()*b.x() - a.x()*b.y());

   // E
   double dE_dxP1 = dM_dxP1*one_over_b;
   double dE_dyP1 = dM_dyP1*one_over_b;
   double dE_dzP1 = dM_dzP1*one_over_b; 

   // M = Eb
   double dE_dxP2 = dM_dxP2*one_over_b + E*(P3.x() - P2.x())*L;
   double dE_dyP2 = dM_dyP2*one_over_b + E*(P3.y() - P2.y())*L;
   double dE_dzP2 = dM_dzP2*one_over_b + E*(P3.z() - P2.z())*L;
	    
   double dE_dxP3 = dM_dxP3*one_over_b - E*(P3.x() - P2.x())*L;
   double dE_dyP3 = dM_dyP3*one_over_b - E*(P3.y() - P2.y())*L;
   double dE_dzP3 = dM_dzP3*one_over_b - E*(P3.z() - P2.z())*L;
	    
   double dE_dxP4 = dM_dxP4*one_over_b;
   double dE_dyP4 = dM_dyP4*one_over_b;
   double dE_dzP4 = dM_dzP4*one_over_b;

   double EFF = E*F*F;
   double JL = J*L;
   double KL = K*L;
   double JK = J*K;

   // x
   dtg.dD_dxP1 = F*dE_dxP1 - EFF*(dH_dxP1 + JL*dK_dxP1 + KL*dJ_dxP1 + JK*dL_dxP1);
   dtg.dD_dxP2 = F*dE_dxP2 - EFF*(dH_dxP2 + JL*dK_dxP2 + KL*dJ_dxP2 + JK*dL_dxP2);
   dtg.dD_dxP3 = F*dE_dxP3 - EFF*(dH_dxP3 + JL*dK_dxP3 + KL*dJ_dxP3 + JK*dL_dxP3);
   dtg.dD_dxP4 = F*dE_dxP4 - EFF*(dH_dxP4 + JL*dK_dxP4 + KL*dJ_dxP4 + JK*dL_dxP4);

   // y
   dtg.dD_dyP1 = F*dE_dyP1 - EFF*(dH_dyP1 + JL*dK_dyP1 + KL*dJ_dyP1 + JK*dL_dyP1);
   dtg.dD_dyP2 = F*dE_dyP2 - EFF*(dH_dyP2 + JL*dK_dyP2 + KL*dJ_dyP2 + JK*dL_dyP2);
   dtg.dD_dyP3 = F*dE_dyP3 - EFF*(dH_dyP3 + JL*dK_dyP3 + KL*dJ_dyP3 + JK*dL_dyP3);
   dtg.dD_dyP4 = F*dE_dyP4 - EFF*(dH_dyP4 + JL*dK_dyP4 + KL*dJ_dyP4 + JK*dL_dyP4);

   // z
   dtg.dD_dzP1 = F*dE_dzP1 - EFF*(dH_dzP1 + JL*dK_dzP1 + KL*dJ_dzP1 + JK*dL_dzP1);
   dtg.dD_dzP2 = F*dE_dzP2 - EFF*(dH_dzP2 + JL*dK_dzP2 + KL*dJ_dzP2 + JK*dL_dzP2);
   dtg.dD_dzP3 = F*dE_dzP3 - EFF*(dH_dzP3 + JL*dK_dzP3 + KL*dJ_dzP3 + JK*dL_dzP3);
   dtg.dD_dzP4 = F*dE_dzP4 - EFF*(dH_dzP4 + JL*dK_dzP4 + KL*dJ_dzP4 + JK*dL_dzP4);

   return dtg;
} 

// Add in the torsion gradients
//
void coot::my_df_torsions_internal(const gsl_vector *v, 
				   void *params, 
				   gsl_vector *df,
				   bool do_rama_torsions) {
   
   int n_torsion_restr = 0; 
   int idx; 

   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params;

   if (restraints->restraints_usage_flag & coot::TORSIONS_MASK) { 

      for (unsigned int i=restraints->restraints_limits_torsions.first; i<=restraints->restraints_limits_torsions.second; i++) {
      
	 if ( (*restraints)[i].restraint_type == coot::TORSION_RESTRAINT) {

	    n_torsion_restr++;

	    idx = 3*((*restraints)[i].atom_index_1); 
	    clipper::Coord_orth P1(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));
	    idx = 3*((*restraints)[i].atom_index_2); 
	    clipper::Coord_orth P2(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));
	    idx = 3*((*restraints)[i].atom_index_3); 
	    clipper::Coord_orth P3(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));
	    idx = 3*((*restraints)[i].atom_index_4); 
	    clipper::Coord_orth P4(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));

	    try { 
	       coot::distortion_torsion_gradients_t dtg =
		  fill_distortion_torsion_gradients(P1, P2, P3, P4);

	       if (! do_rama_torsions) { 
		  // 
		  // use period 

		  double diff = 99999.9; 
		  double tdiff; 
		  double trial_target; 
		  int per = (*restraints)[i].periodicity;

		  if (dtg.theta < 0.0) dtg.theta += 360.0; 

		  for(int iper=0; iper<per; iper++) { 
		     trial_target = (*restraints)[i].target_value + double(iper)*360.0/double(per); 
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
			       <<  " and target is " << (*restraints)[i].target_value 
			       << " and diff is " << diff
			       << " and periodicity: " << (*restraints)[i].periodicity << std::endl;

		  double tt = tan(clipper::Util::d2rad(dtg.theta));
		  double torsion_scale = (1.0/(1+tt*tt)) *
		     clipper::Util::rad2d(1.0);

		  double weight = 1/((*restraints)[i].sigma * (*restraints)[i].sigma);

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
	    
		  if (! (*restraints)[i].fixed_atom_flags[0]) { 
		     idx = 3*((*restraints)[i].atom_index_1);
		     *gsl_vector_ptr(df, idx  ) += xP1_contrib;
		     *gsl_vector_ptr(df, idx+1) += yP1_contrib;
		     *gsl_vector_ptr(df, idx+2) += zP1_contrib;
		  }

		  if (! (*restraints)[i].fixed_atom_flags[1]) { 
		     idx = 3*((*restraints)[i].atom_index_2);
		     *gsl_vector_ptr(df, idx  ) += xP2_contrib;
		     *gsl_vector_ptr(df, idx+1) += yP2_contrib;
		     *gsl_vector_ptr(df, idx+2) += zP2_contrib;
		  }

		  if (! (*restraints)[i].fixed_atom_flags[2]) { 
		     idx = 3*((*restraints)[i].atom_index_3);
		     *gsl_vector_ptr(df, idx  ) += xP3_contrib;
		     *gsl_vector_ptr(df, idx+1) += yP3_contrib;
		     *gsl_vector_ptr(df, idx+2) += zP3_contrib;
		  }

		  if (! (*restraints)[i].fixed_atom_flags[3]) { 
		     idx = 3*((*restraints)[i].atom_index_4);
		     *gsl_vector_ptr(df, idx  ) += xP4_contrib;
		     *gsl_vector_ptr(df, idx+1) += yP4_contrib;
		     *gsl_vector_ptr(df, idx+2) += zP4_contrib;
		  }
	       }
	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << "Caught runtime_error" << rte.what() << std::endl;
	    } 
	 } 
      }
   }
}

// Add in the torsion gradients
//
void coot::my_df_rama(const gsl_vector *v, 
		      void *params, 
		      gsl_vector *df) {

   // First calculate the torsions:
   // theta = arctan(E/G); 
   // where E = a.(bxc) and G = -a.c + (a.b)(b.c)

   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params;

   try { 

      if (restraints->restraints_usage_flag & coot::RAMA_PLOT_MASK) { 
     
	 for (int i=0; i<restraints->size(); i++) {
      
	    if ( (*restraints)[i].restraint_type == coot::RAMACHANDRAN_RESTRAINT) {

	       int idx;
	       coot::simple_restraint rama_restraint = (*restraints)[i];

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

	       if (! (*restraints)[i].fixed_atom_flags[0]) { 
		  idx = 3*((*restraints)[i].atom_index_1);
		  gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + xP1_contrib);
		  gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + yP1_contrib);
		  gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + zP1_contrib);
	       }

	       if (! (*restraints)[i].fixed_atom_flags[1]) { 
		  idx = 3*((*restraints)[i].atom_index_2);
		  gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + xP2_contrib);
		  gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + yP2_contrib);
		  gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + zP2_contrib);
	       }

	       if (! (*restraints)[i].fixed_atom_flags[2]) { 
		  idx = 3*((*restraints)[i].atom_index_3);
		  gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + xP3_contrib);
		  gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + yP3_contrib);
		  gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + zP3_contrib);
	       }

	       if (! (*restraints)[i].fixed_atom_flags[3]) { 
		  idx = 3*((*restraints)[i].atom_index_4);
		  gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + xP4_contrib);
		  gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + yP4_contrib);
		  gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + zP4_contrib);
	       }

	       if (! (*restraints)[i].fixed_atom_flags[4]) { 
		  idx = 3*((*restraints)[i].atom_index_5);
		  gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + xP5_contrib);
		  gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + yP5_contrib);
		  gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + zP5_contrib);
	       }
	    }
	 }
      }
   }
   catch (const std::runtime_error &rte) {
      std::cout << "ERROR:: my_df_rama() caught " << rte.what() << std::endl;
   } 
}

//  the chiral volumes
void 
coot::my_df_chiral_vol(const gsl_vector *v, void *params, gsl_vector *df) { 

   int n_chiral_vol_restr = 0;
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params;
   int idx;
   double cv;
   double distortion;
   
   if (restraints->restraints_usage_flag & coot::CHIRAL_VOLUME_MASK) {
      
      for (unsigned int i=restraints->restraints_limits_chirals.first; i<=restraints->restraints_limits_chirals.second; i++) {
	 
	 if ( (*restraints)[i].restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) {

	    n_chiral_vol_restr++;

	    idx = 3*(*restraints)[i].atom_index_centre;
	    clipper::Coord_orth centre(gsl_vector_get(v, idx),
				       gsl_vector_get(v, idx+1),
				       gsl_vector_get(v, idx+2));

	    idx = 3*( (*restraints)[i].atom_index_1);
	    clipper::Coord_orth a1(gsl_vector_get(v, idx),
				   gsl_vector_get(v, idx+1),
				   gsl_vector_get(v, idx+2));
	    idx = 3*( (*restraints)[i].atom_index_2);
	    clipper::Coord_orth a2(gsl_vector_get(v, idx),
				   gsl_vector_get(v, idx+1),
				   gsl_vector_get(v, idx+2));
	    idx = 3*( (*restraints)[i].atom_index_3);
	    clipper::Coord_orth a3(gsl_vector_get(v, idx),
				   gsl_vector_get(v, idx+1),
				   gsl_vector_get(v, idx+2));

	    clipper::Coord_orth a = a1 - centre;
	    clipper::Coord_orth b = a2 - centre;
	    clipper::Coord_orth c = a3 - centre;

	    cv = clipper::Coord_orth::dot(a, clipper::Coord_orth::cross(b,c));

	    distortion = cv - (*restraints)[i].target_chiral_volume;
	    
// 	    std::cout << "---- xxx ---- DEBUG:: chiral volume deriv: " 
// 		      << cv << " chiral distortion " 
// 		      << distortion << "\n";
	    // distortion /= ((*restraints)[i].sigma * (*restraints)[i].sigma);

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

	    double s = 2*distortion/((*restraints)[i].sigma * (*restraints)[i].sigma);

	    if (!(*restraints)[i].fixed_atom_flags[0]) { 
	       idx = 3*( (*restraints)[i].atom_index_centre);
	       gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + s * P0_x_contrib); 
	       gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + s * P0_y_contrib); 
	       gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + s * P0_z_contrib);
	    }
	       
	    if (!(*restraints)[i].fixed_atom_flags[1]) { 
	       idx = 3*( (*restraints)[i].atom_index_1);
	       gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + s * P1_x_contrib); 
	       gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + s * P1_y_contrib); 
	       gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + s * P1_z_contrib);
	    }

	    if (!(*restraints)[i].fixed_atom_flags[2]) { 
	       idx = 3*( (*restraints)[i].atom_index_2);
	       gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + s * P2_x_contrib); 
	       gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + s * P2_y_contrib); 
	       gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + s * P2_z_contrib);
	    }

	    if (!(*restraints)[i].fixed_atom_flags[3]) { 
	       idx = 3*( (*restraints)[i].atom_index_3);
	       gsl_vector_set(df, idx,   gsl_vector_get(df, idx)   + s * P3_x_contrib); 
	       gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + s * P3_y_contrib); 
	       gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + s * P3_z_contrib);
	    }
	 }
      }
   }
} 

// manipulate the gsl_vector_get *df
// 
void
coot::my_df_planes(const gsl_vector *v, 
		   void *params, 
		   gsl_vector *df) {

   
   // first extract the object from params 
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params;

   int idx; 

   if (restraints->restraints_usage_flag & coot::PLANES_MASK) {
      
      int n_plane_atoms;
      double devi_len;
      double weight;

      for (unsigned int i=restraints->restraints_limits_planes.first; i<=restraints->restraints_limits_planes.second; i++) {
       
	 if ( (*restraints)[i].restraint_type == coot::PLANE_RESTRAINT) {

	    const simple_restraint &plane_restraint = (*restraints)[i];

	    if (false) { // debug
	       std::cout << "restraint index " << i << " plane with fixed atom indices:\n";
	       for (std::size_t jj=0; jj<plane_restraint.fixed_atom_flags.size(); jj++) {
		  std::cout << " " << plane_restraint.fixed_atom_flags[jj];
		  mmdb::Atom *at = restraints->get_atom(plane_restraint.plane_atom_index[jj].first);
		  std::cout << "    " << atom_spec_t(at) << std::endl;
	       }
	       std::cout << "\n";
	    }

	    coot::plane_distortion_info_t plane_info =
	       distortion_score_plane_internal(plane_restraint, v);
	    n_plane_atoms = plane_restraint.plane_atom_index.size();
	    // weight = 1/((*restraints)[i].sigma * (*restraints)[i].sigma);
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

		  if (false) { // debug
		     if (plane_restraint.plane_atom_index.size() >= 4) {
			std::cout << "   gradients plane_restraint " << plane_restraint << " "
				  << d.format() << std::endl;
			int idxp = plane_restraint.plane_atom_index[j].first;
			mmdb::Atom *at = restraints->get_atom(idxp);
			std::cout << "   plane atom " << j << " " << atom_spec_t(at) << std::endl;
		     }
		  }


		  *gsl_vector_ptr(df, idx  ) += d.dx();
		  *gsl_vector_ptr(df, idx+1) += d.dy();
		  *gsl_vector_ptr(df, idx+2) += d.dz();
	       }
	    }
	 }
      }
   }
}


// manipulate the gsl_vector_get *df
// 
void
coot::my_df_parallel_planes(const gsl_vector *v, 
			    void *params, 
			    gsl_vector *df) {
   int idx;
   double devi_len;

   // first extract the object from params 
   //
   coot::restraints_container_t *restraints = (coot::restraints_container_t *)params;
   
   if (restraints->restraints_usage_flag & coot::PARALLEL_PLANES_MASK) {
      for (unsigned int i=restraints->restraints_limits_parallel_planes.first; i<=restraints->restraints_limits_parallel_planes.second; i++) {
       
	 if ( (*restraints)[i].restraint_type == coot::PARALLEL_PLANES_RESTRAINT) {
	    const simple_restraint &ppr = (*restraints)[i];

	    unsigned int first_atoms_size = ppr.plane_atom_index.size();
	    
	    plane_distortion_info_t plane_info =
	       distortion_score_2_planes(ppr.plane_atom_index, ppr.atom_index_other_plane, ppr.sigma, v);

	    // first plane
	    unsigned int n_plane_atoms = ppr.plane_atom_index.size();
	    double weight = 1/(ppr.sigma * ppr.sigma);
	    // hack the weight - needs a better fix than this
	    weight *= 8.0;
	    for (unsigned int j=0; j<n_plane_atoms; j++) {
	       if (! ppr.fixed_atom_flags[j] ) { 
		  idx = 3*ppr.plane_atom_index[j].first;
		  devi_len =
		     plane_info.abcd[0]*(gsl_vector_get(v,idx  ) - plane_info.centre_1.x()) + 
		     plane_info.abcd[1]*(gsl_vector_get(v,idx+1) - plane_info.centre_1.y()) +
		     plane_info.abcd[2]*(gsl_vector_get(v,idx+2) - plane_info.centre_1.z()) -
		     plane_info.abcd[3];

		  clipper::Grad_orth<double> d(2.0 * weight * devi_len * plane_info.abcd[0],
					       2.0 * weight * devi_len * plane_info.abcd[1],
					       2.0 * weight * devi_len * plane_info.abcd[2]);

		  // std::cout << "pp first plane devi len " << devi_len << " and gradients "
		  // << d.format() << std::endl;

		  gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + d.dx());
		  gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + d.dy());
		  gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + d.dz());
	       }
	    }

	    // second plane
	    n_plane_atoms = ppr.atom_index_other_plane.size();
	    for (unsigned int j=0; j<n_plane_atoms; j++) {
	       if (! ppr.fixed_atom_flags_other_plane[j] ) {
		  idx = 3*ppr.atom_index_other_plane[j].first;
		  devi_len =
		     plane_info.abcd[0]*(gsl_vector_get(v,idx  ) - plane_info.centre_2.x()) + 
		     plane_info.abcd[1]*(gsl_vector_get(v,idx+1) - plane_info.centre_2.y()) +
		     plane_info.abcd[2]*(gsl_vector_get(v,idx+2) - plane_info.centre_2.z()) -
		     plane_info.abcd[3];

		  clipper::Grad_orth<double> d(2.0 * weight * devi_len * plane_info.abcd[0],
					       2.0 * weight * devi_len * plane_info.abcd[1],
					       2.0 * weight * devi_len * plane_info.abcd[2]);

		  gsl_vector_set(df, idx,   gsl_vector_get(df, idx  ) + d.dx());
		  gsl_vector_set(df, idx+1, gsl_vector_get(df, idx+1) + d.dy());
		  gsl_vector_set(df, idx+2, gsl_vector_get(df, idx+2) + d.dz());
	       }
	    }
	 }
      }
   }
}



double
coot::restraints_container_t::electron_density_score_at_point(const clipper::Coord_orth &ao) const {
      double dv; 
      
      clipper::Coord_frac af = ao.coord_frac(xmap.cell()); 
      clipper::Coord_map  am = af.coord_map(xmap.grid_sampling()); 
      // clipper::Interp_linear::interp(map, am, dv); 
      clipper::Interp_cubic::interp(xmap, am, dv); 
      
      return dv;  
}

clipper::Grad_orth<double>
coot::restraints_container_t::electron_density_gradient_at_point(const clipper::Coord_orth &ao) const {
   
   clipper::Grad_map<double> grad;
   double dv;
   
   clipper::Coord_frac af = ao.coord_frac(xmap.cell()); 
   clipper::Coord_map  am = af.coord_map(xmap.grid_sampling()); 
   clipper::Interp_cubic::interp_grad(xmap, am, dv, grad);
   clipper::Grad_frac<double> grad_frac = grad.grad_frac(xmap.grid_sampling());
   return grad_frac.grad_orth(xmap.cell());
} 

#endif // HAVE_GSL
