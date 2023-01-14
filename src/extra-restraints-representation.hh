/* src/extra-restraints-representation.hh
 * 
 * Copyright 2015 by Medical Research Council
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

#ifndef EXTRA_RESTRAINTS_REPRESENTATION_HH
#define EXTRA_RESTRAINTS_REPRESENTATION_HH

#include <vector>
#include <clipper/core/coords.h>
#include "coot-utils/coot-coord-utils.hh" // for lsq

namespace coot {
   // ------------ extra restraints (e.g. user-defined) ------------
   //
   class extra_restraints_representation_t {
   public:

      class extra_bond_restraints_respresentation_t {
      public:
	 clipper::Coord_orth first;
	 clipper::Coord_orth second;
	 double target_dist;
	 double esd;
	 extra_bond_restraints_respresentation_t(const clipper::Coord_orth &f,
						 const clipper::Coord_orth &s,
						 double d,
						 double e) : first(f), second(s) {
	    target_dist = d;
	    esd = e;
	 }
         double length_delta() const {
            double dd = clipper::Coord_orth(first - second).lengthsq();
            return std::sqrt(dd) - target_dist;
         }
         double distortion_score_GM(const double &sigma, const double &alpha) const {
            double bl = clipper::Coord_orth::length(first, second);
            double bit = bl - target_dist;
            double z = bit/sigma;
            double distortion = z*z/(1+alpha*z*z);
            return distortion;
         }
      };

      class extra_parallel_planes_restraints_representation_t {
      public:
	 clipper::Coord_orth ring_centre;
	 clipper::Coord_orth plane_projection_point;
	 clipper::Coord_orth normal; // for the ring plane
	 double ring_radius;
	 double pp_radius; // projection point
	 extra_parallel_planes_restraints_representation_t(const clipper::Coord_orth &rc,
							   const clipper::Coord_orth &ppp,
							   const clipper::Coord_orth &norm,
							   double r1, double r2) : ring_centre(rc), plane_projection_point(ppp), normal(norm) {
	    ring_radius = r1;
	    pp_radius = r2;
	 }
      };

      extra_restraints_representation_t() {
	 prosmart_restraint_display_limit_high = 0;
	 prosmart_restraint_display_limit_low  = 0;
      }
      
      std::vector<extra_bond_restraints_respresentation_t> bonds;
      double prosmart_restraint_display_limit_high; // show only those below this
      double prosmart_restraint_display_limit_low;  // and above this (n-sigma)
      std::vector<extra_parallel_planes_restraints_representation_t> parallel_planes;
      
      void clear() {
	 bonds.clear();
	 parallel_planes.clear();
      }
      void add_bond(const clipper::Coord_orth &pt1, const clipper::Coord_orth &pt2) {
	 extra_bond_restraints_respresentation_t br(pt1, pt2, -1, -1);
	 bonds.push_back(br);

      }
      void add_bond(const clipper::Coord_orth &pt1, const clipper::Coord_orth &pt2,
		    const double &d, const double &e) {
	 extra_bond_restraints_respresentation_t br(pt1, pt2, d, e);
	 bonds.push_back(br);
      }

      // maybe extra parameters are needed here (e.g. for colouring later, perhaps).
      void add_parallel_plane(const lsq_plane_info_t &pi_1,
			      const lsq_plane_info_t &pi_2);
   };
}


#endif // EXTRA_RESTRAINTS_REPRESENTATION_HH
