/* ideal/simple-restraint-density.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 by The University of York
 * Copyright 2008, 2009, 2010  by The University of Oxford
 * Copyright 2013, 2014, 2015, 2016 by Medical Research Council
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

// #define ANALYSE_REFINEMENT_TIMING

#include <string.h> // for strcmp

#ifdef ANALYSE_REFINEMENT_TIMING
#ifdef _MSC_VER
#include <time.h>
#else
#include <sys/time.h>
#endif
#endif // ANALYSE_REFINEMENT_TIMING


#include <fstream>
#include <algorithm> // for sort
#include <stdexcept>
#include <iomanip>

#ifdef HAVE_CXX_THREAD
#include <thread>
#include <chrono>
#endif // HAVE_CXX_THREAD

#include "utils/split-indices.hh"
#include "geometry/mol-utils.hh"
#include "geometry/main-chain.hh"
#include "simple-restraint.hh"

//
#include "coot-utils/coot-coord-extras.hh"  // is_nucleotide_by_dict

// #include "mmdb.h" // for printing of mmdb::Atom pointers as info not raw
                     // pointers.  Removed. Too much (linking issues in)
                     // Makefile pain.

#include "compat/coot-sysdep.h"

#include "utils/logging.hh"
extern logging logger;

// Electron-density scoring and gradients for restraints_container_t

// Ah, but (c.f. distortion) we want to return a low value for a good
// fit and a high one for a bad.
double
coot::electron_density_score(const gsl_vector *v, void *params) {

   coot::restraints_container_t *restraints_p = static_cast<restraints_container_t *> (params);
   if (restraints_p->include_map_terms() == 1) {
      return electron_density_score_from_restraints(v, restraints_p);
   } else {
      return 0.0;
   }
}


// Ah, but (c.f. distortion) we want to return a low value for a good
// fit and a high one for a bad.
double
coot::electron_density_score_from_restraints_simple(const gsl_vector *v,
					     coot::restraints_container_t *restraints_p) {

#ifdef HAVE_CXX_THREAD
   auto tp_1 = std::chrono::high_resolution_clock::now();
#endif
   // We weight and sum to get the score and negate.
   //
   double score = 0.0;

   if (restraints_p->include_map_terms() == 1) {

      unsigned int n_atoms = restraints_p->get_n_atoms();
      for (unsigned int iat=0; iat<n_atoms; iat++) {
	 bool use_it = false;
	 if (restraints_p->use_map_gradient_for_atom[iat]) {

	    int idx = 3 * iat;
	    clipper::Coord_orth ao(gsl_vector_get(v,idx),
				   gsl_vector_get(v,idx+1),
				   gsl_vector_get(v,idx+2));

	    score += restraints_p->Map_weight() *
	       restraints_p->atom_z_occ_weight[iat] *
	       restraints_p->electron_density_score_at_point(ao);
	 }
      }
   }
#ifdef HAVE_CXX_THREAD
   auto tp_2 = std::chrono::high_resolution_clock::now();
   auto d21 = std::chrono::duration_cast<std::chrono::microseconds>(tp_2 - tp_1).count();
   // std::cout << "info:: f electron_density: " << d21 << " microseconds\n";
#endif // HAVE_CXX_THREAD

   return -score;
}


double
coot::electron_density_score_from_restraints(const gsl_vector *v,
					     coot::restraints_container_t *restraints_p) {

   double score = 0.0;
#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

   std::vector<std::pair<unsigned int, unsigned int> > &ranges = restraints_p->m_atom_index_ranges;

   // they are no longer "threads" - they are restraints sets - to prevent
   // "thread starvation"
   std::atomic<unsigned int> done_count_for_restraints_sets(0);

   if (restraints_p->thread_pool_p) {

      double results[1024]; // we will always have less than 1024 threads

      unsigned int n_ranges = ranges.size(); // clang scan-build fix.
      for(unsigned int i=0; i<n_ranges; i++) {
         results[i] = 0.0;
	 restraints_p->thread_pool_p->push(electron_density_score_from_restraints_using_atom_index_range,
					   v, std::cref(ranges[i]), restraints_p, &(results[i]),
					   std::ref(done_count_for_restraints_sets));
      }
      while (done_count_for_restraints_sets < ranges.size())
	 std::this_thread::sleep_for(std::chrono::nanoseconds(300));

      // consolidate
      for(unsigned int i=0; i<n_ranges; i++)
	 score += results[i];
   } else {
      // null thread pool. restraints_container_t was created without a call to
      // set the thread_pool. Happens in crankshafting currently.
      electron_density_score_from_restraints_using_atom_index_range(0, v, ranges[0], restraints_p, &score,
								    done_count_for_restraints_sets);
   }

#else
   std::cout << __FUNCTION__ << " no thread pool" << std::endl;
#endif // HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
   return score;
}

#ifdef HAVE_CXX_THREAD
// atom_index_range works "as expected"
// so given atom_index_range of 0,10 we start at the first value (0) and check that the
// current value is less than the atom_index_range.second (10):
// ie. density values for atom indices 0 to 9 inclusive are added.
//
void
coot::electron_density_score_from_restraints_using_atom_index_range(int thread_index,
                                             const gsl_vector *v,
					     const std::pair<unsigned int, unsigned int> &atom_index_range,
					     coot::restraints_container_t *restraints_p,
					     double *result,
					     std::atomic<unsigned int> &done_count_for_threads) {

   // auto tp_1 = std::chrono::high_resolution_clock::now();

   // We weight and sum to get the score and negate.
   //
   double score = 0;

   if (restraints_p->include_map_terms() == 1) {

      if (false)
	 std::cout << "debug:: " << __FUNCTION__ << "() range: "
		   << atom_index_range.first << " " << atom_index_range.second << std::endl;

      for (unsigned int iat=atom_index_range.first; iat<atom_index_range.second; iat++) {

	 // we get a crash around here. Protect from wrong index into v vector
	 if (iat < restraints_p->get_n_atoms()) {
	    if (restraints_p->use_map_gradient_for_atom[iat]) {

	       int idx = 3 * iat;
	       clipper::Coord_orth ao(gsl_vector_get(v,idx),
				      gsl_vector_get(v,idx+1),
				      gsl_vector_get(v,idx+2));

               //               std::cout << "ao:" << ao.format() << std::endl; // prograam terminated before
                                                                                // sphere-refine had finished.

	       score += restraints_p->Map_weight() *
		  restraints_p->atom_z_occ_weight[iat] *
		  restraints_p->electron_density_score_at_point(ao);
	    }
	 } else {
	    std::cout << "ERROR:: electron_density_score_from_restraints_using_atom_index_range "
		      << " caught bad atom index " << iat << " " << restraints_p->get_n_atoms()
		      << std::endl;
	 }
      }
   }

   // auto tp_2 = std::chrono::high_resolution_clock::now();
   // auto d21 = std::chrono::duration_cast<std::chrono::microseconds>(tp_2 - tp_1).count();
   // std::cout << "info:: f electron_density: " << d21 << " microseconds\n";
   // return -score;

   *result = -score;
   done_count_for_threads++; // atomic
}
#endif




// Note that the gradient for the electron density is opposite to that
// of the gradient for the geometry (consider a short bond on the edge
// of a peak - in that case the geometry gradient will be negative as
// the bond is lengthened and the electron density gradient will be
// positive).
//
// So we want to change that positive gradient for a low score when
// the atoms coinside with the density - hence the contributions that
// we add are negated.
//
void coot::my_df_electron_density(const gsl_vector *v,
				  void *params,
				  gsl_vector *df) {

   // first extract the object from params
   //
   coot::restraints_container_t *restraints_p = static_cast<restraints_container_t *> (params);
   if (restraints_p->include_map_terms() == 1) {
      my_df_electron_density_single(v, restraints_p, df, 0, v->size/3);
   }
}



#ifdef HAVE_CXX_THREAD


// restraints are modified by atomic done_count_for_threads changing.
//
// multi-threaded inner
//
void coot::my_df_electron_density_threaded_single(int thread_idx, const gsl_vector *v,
						  coot::restraints_container_t *restraints,
						  gsl_vector *df,
						  int atom_idx_start, int atom_idx_end,
						  std::atomic<unsigned int> &done_count_for_threads) {

   for (int iat=atom_idx_start; iat<atom_idx_end; ++iat) {
      if (restraints->use_map_gradient_for_atom[iat]) {

	 int idx = 3 * iat;
	 clipper::Coord_orth ao(gsl_vector_get(v,idx),
				gsl_vector_get(v,idx+1),
				gsl_vector_get(v,idx+2));

	 clipper::Grad_orth<double> grad_orth = restraints->electron_density_gradient_at_point(ao);
	 float zs = restraints->Map_weight() * restraints->atom_z_occ_weight[iat];

	 if (false) {
	    std::cout << "electron density df: adding "
		      <<  - zs * grad_orth.dx() << " "
		      <<  - zs * grad_orth.dy() << " "
		      <<  - zs * grad_orth.dz() << " to "
		      <<  gsl_vector_get(df, idx  ) << " "
		      <<  gsl_vector_get(df, idx+1) << " "
		      <<  gsl_vector_get(df, idx+2) << "\n";
	 }

	 // 	    gsl_vector_set(df, i,   gsl_vector_get(df, i  ) - zs * grad_orth.dx());
	 // 	    gsl_vector_set(df, i+1, gsl_vector_get(df, i+1) - zs * grad_orth.dy());
	 // 	    gsl_vector_set(df, i+2, gsl_vector_get(df, i+2) - zs * grad_orth.dz());


         // we no longer do this sort of locking

	 // use atomic lock to access derivs of atom idx
	 // unsigned int unlocked = 0;
	 // while (! restraints->gsl_vector_atom_pos_deriv_locks.get()[idx].compare_exchange_weak(unlocked, 1)) {
	 //    std::this_thread::sleep_for(std::chrono::nanoseconds(10));
	 //    unlocked = 0;
	 // }

	 *gsl_vector_ptr(df, idx  ) -= zs * grad_orth.dx();
	 *gsl_vector_ptr(df, idx+1) -= zs * grad_orth.dy();
	 *gsl_vector_ptr(df, idx+2) -= zs * grad_orth.dz();
	 // restraints->gsl_vector_atom_pos_deriv_locks.get()[idx] = 0; // unlock
      }
   }

   // int sleepy_time = 10000 + int(2000*float(coot::util::random())/float(RAND_MAX));
   // std::this_thread::sleep_for(std::chrono::microseconds(sleepy_time));
   ++done_count_for_threads; // atomic
}
#endif


//
void coot::my_df_electron_density_single(const gsl_vector *v,
					 coot::restraints_container_t *restraints,
					 gsl_vector *df,
					 int atom_idx_start, int atom_idx_end) {

   //    std::cout << "debug:: in my_df_electron_density_single() " << atom_idx_start << " " << atom_idx_end
   // << std::endl;

   for (int iat=atom_idx_start; iat<atom_idx_end; ++iat) {
      if (restraints->use_map_gradient_for_atom[iat]) {

	 int idx = 3 * iat;
	 clipper::Coord_orth ao(gsl_vector_get(v,idx),
				gsl_vector_get(v,idx+1),
				gsl_vector_get(v,idx+2));

	 clipper::Grad_orth<double> grad_orth = restraints->electron_density_gradient_at_point(ao);
	 float zs = restraints->Map_weight() * restraints->atom_z_occ_weight[iat];

	 if (false) {
	    std::cout << "electron density df: adding "
		      <<  - zs * grad_orth.dx() << " "
		      <<  - zs * grad_orth.dy() << " "
		      <<  - zs * grad_orth.dz() << " to "
		      <<  gsl_vector_get(df, idx  ) << " "
		      <<  gsl_vector_get(df, idx+1) << " "
		      <<  gsl_vector_get(df, idx+2) << "\n";
	 }

	 // 	    gsl_vector_set(df, i,   gsl_vector_get(df, i  ) - zs * grad_orth.dx());
	 // 	    gsl_vector_set(df, i+1, gsl_vector_get(df, i+1) - zs * grad_orth.dy());
	 // 	    gsl_vector_set(df, i+2, gsl_vector_get(df, i+2) - zs * grad_orth.dz());

	 *gsl_vector_ptr(df, idx  ) -= zs * grad_orth.dx();
	 *gsl_vector_ptr(df, idx+1) -= zs * grad_orth.dy();
	 *gsl_vector_ptr(df, idx+2) -= zs * grad_orth.dz();
      }
   }
}

void coot::my_df_electron_density_old (gsl_vector *v,
				       void *params,
				       gsl_vector *df) {

   // first extract the object from params
   //
   coot::restraints_container_t *restraints =
      (coot::restraints_container_t *)params;

   if (restraints->include_map_terms() == 1) {

      double new_S_minu, new_S_plus, tmp, val;

      std::cout << "density_gradients" << std::endl;
      for (unsigned int i=0; i<v->size; i++) {

	 tmp = gsl_vector_get(v, i);
	 gsl_vector_set(v, i, tmp+0.01);
	 new_S_plus = coot::electron_density_score(v, params);
	 gsl_vector_set(v, i, tmp-0.01);
	 new_S_minu = coot::electron_density_score(v, params);
	 // new_S_minu = 2*tmp - new_S_plus;

	 // restore the initial value:
	 gsl_vector_set(v, i, tmp);

	 val = (new_S_plus - new_S_minu)/(2*0.01);
	 std::cout << "density gradient: " << i << " " << val << std::endl;

	 // add this density term to the gradient
	 gsl_vector_set(df, i, gsl_vector_get(df, i) + val);
      }
   }
}

// Compute both f and df together.
void coot::my_fdf(const gsl_vector *x, void *params,
		  double *f, gsl_vector *df) {

   // 20170423 these can be done in parallel? ... check the timings at least.
   *f = coot::distortion_score(x, params);
    coot::my_df(x, params, df);
}


