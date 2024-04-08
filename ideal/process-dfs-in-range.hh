/*
 * ideal/process-df-in-range.hh
 * 
 * Copyright 2018 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 */


#ifndef COOT_PROCESS_DFS_IN_RANGE_HH
#define COOT_PROCESS_DFS_IN_RANGE_HH

#include <vector>
#include <atomic>

#include <gsl/gsl_multimin.h>
#include "simple-restraint.hh"

namespace coot {

   void process_dfs_in_range(int thread_index,
			     const std::vector<std::size_t> &restraints_indices,
			     coot::restraints_container_t *restraints_p,
			     const gsl_vector *v,
			     std::vector<double> &results,  // fill results
			     std::atomic<unsigned int> &done_count_for_threads);

   void
   process_dfs_bond(const coot::simple_restraint &rest,
		    const gsl_vector *v,
		    std::vector<double> &result); // fill results

   void
   process_dfs_angle(const coot::simple_restraint &restraint,
		     const gsl_vector *v,
		     std::vector<double> &results); // fill results

   void
   process_dfs_torsion(const coot::simple_restraint &restraint,
		       const gsl_vector *v,
		       std::vector<double> &results);

   void
   process_dfs_chiral_volume(const simple_restraint &restraint,
			     const gsl_vector *v,
			     std::vector<double> &results); // fill results

   void
   process_dfs_plane(const simple_restraint &restraint,
		     const gsl_vector *v,
		     std::vector<double> &results);
   void
   process_dfs_improper_dihedral(const simple_restraint &restraint,
		     const gsl_vector *v,
		     std::vector<double> &results);

   void
   process_dfs_parallel_planes(const simple_restraint &restraint,
			       const gsl_vector *v,
			       std::vector<double> &results);

   void
   process_dfs_non_bonded(const simple_restraint &restraint,
			  const gsl_vector *v,
			  std::vector<double> &results); // fill results

   void
   process_dfs_non_bonded_lennard_jones(const simple_restraint &restraint,
					const double &lj_epsilon,
					const gsl_vector *v,
					std::vector<double> &results);

   void
   process_dfs_geman_mcclure_distance(const simple_restraint &restraint,
				      const double &alpha,
				      const gsl_vector *v,
				      std::vector<double> &results);

   void process_dfs_trans_peptide(const coot::simple_restraint &restraint,
				  const gsl_vector *v,
				  std::vector<double> &results);

   void process_dfs_target_position(const simple_restraint &restraint,
				    const double &sf,
				    const gsl_vector *v,
				    std::vector<double> &results);

   void process_dfs_start_position(const simple_restraint &restraint,
				   const gsl_vector *v,
				   std::vector<double> &results);

   void process_dfs_rama(const simple_restraint &restraint,
			 const restraints_container_t *restraints_p, // for plots and multipliers
			 const gsl_vector *v,
			 std::vector<double> &results);

   void process_electron_density_dfs_for_atoms(int thread_idx,
					       const std::vector<std::size_t> &atom_indices,
					       const restraints_container_t *restraints_p,
					       const gsl_vector *v, gsl_vector *df,
					       std::atomic<unsigned int> &done_count_for_threads);


}

#endif // COOT_PROCESS_DFS_IN_RANGE_HH
