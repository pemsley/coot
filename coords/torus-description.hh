/*
 * coords/torus-description.hh
 * 
 * Copyright 2020 by Medical Research Council
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
 */



#ifndef COORDS_TORUS_DESCRIPTION_HH
#define COORDS_TORUS_DESCRIPTION_HH

namespace coot {

   // For OpenGL solid model, it is much better to draw a torus than a
   // set of short sticks (particularly for the ring representing
   // aromaticity).  So now (20100831 Bond_lines_container contains a
   // number of torus descriptions).
   //
   class torus_description_t {
   public:
      double inner_radius;
      double outer_radius;
      int n_sides;
      int n_rings;
      clipper::Coord_orth centre;
      clipper::Coord_orth normal;
      torus_description_t(const clipper::Coord_orth &pt,
			  const clipper::Coord_orth &normal_in,
			  double ir1, double ir2, int n1, int n2) {
	 inner_radius = ir1;
	 outer_radius = ir2;
	 n_sides = n1;
	 n_rings = n2;
	 centre = pt;
	 normal = normal_in;
      }
   };

}

#endif // COORDS_TORUS_DESCRIPTION_HH
