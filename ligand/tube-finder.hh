/*
 * ligand/tube-finder.hh
 *
 * Copyright 2021 by Medical Research Council
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
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */
#ifndef TUBE_FINDER_HH
#define TUBE_FINDER_HH

#include <clipper/core/xmap.h>
#include "gsl/gsl_multimin.h"

namespace coot {

   class tube_finder_t {

      class simplex_param_t {
      public:
         clipper::Coord_orth test_point_centre_orig;
         std::vector<clipper::Coord_orth> original_positions;
         const clipper::Xmap<float> *xmap;
      };
      
      // return the fitted rtop_orth
      std::pair<bool, clipper::Coord_orth> 
      fit_stick_in_map(const clipper::Coord_orth &test_point,
                       const clipper::Mat33<double> &orientation,
                       const clipper::Xmap<float> &xmap,
                       float density_level_crit) const;
      std::pair<bool, clipper::Coord_orth>
      fit_to_map_by_simplex_rigid(const clipper::Coord_orth &test_point,
                                  const clipper::Mat33<double> &orientation,
                                  const clipper::Xmap<float> &xmap,
                                  float density_level_crit) const;
      static double my_f_simplex_rigid_internal(const gsl_vector *v, void *params);
      clipper::Coord_orth apply_shifts_rigid_internal(gsl_vector *x, const simplex_param_t &par) const;
      static clipper::RTop_orth construct_matrix(const gsl_vector *v);
      // return negative on failure
      static float sphere_variance(const clipper::Coord_orth &centre_point,
                                   const std::vector<clipper::Coord_orth> &sphere_points,
                                   float radius,
                                   const clipper::Xmap<float> &xmap);

   public:
      explicit tube_finder_t(const clipper::Xmap<float> &xmap);
      std::vector<clipper::Coord_orth> get_positions() const;
      std::vector<clipper::Coord_orth> positions; // may be separate in the strands and helices

   };
}

#endif // TUBE_FINDER_HH
