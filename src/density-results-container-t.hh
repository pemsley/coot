/*
 * src/density-results-container-t.hh
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
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */


#ifndef DENSITY_RESULTS_CONTAINER_T
#define DENSITY_RESULTS_CONTAINER_T

class density_results_t {
   public:

   // angle in degrees
   density_results_t(const clipper::Coord_orth &p, double angle_in, float density_value_in) :
         position(p), angle(angle_in), density_value(density_value_in)  {}
   clipper::Coord_orth position;
   float angle; // degrees
   float density_value;
};

class density_results_container_t {

   public:
     std::vector<density_results_t> scored_points;

};
#endif // DENSITY_RESULTS_CONTAINER_T
