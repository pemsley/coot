/*
 * coot-utils/fib-sphere.cc
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
 *
 */

#include <math.h>
#include "fib-sphere.hh"

std::vector<clipper::Coord_orth>
coot::fibonacci_sphere(unsigned int n_samples) {

   // I watched the Coding Adventure on Boids
   double phi = M_PI * (3.0 - sqrt(5.0));

   std::vector<clipper::Coord_orth> points(n_samples);
   for (unsigned int i=0; i<n_samples; i++) {
      double d(n_samples - 1.0);
      double y = 1.0 - (static_cast<double>(i)/d) * 2.0;
      double radius = sqrt(1.0 - y * y);
      double theta = phi * static_cast<double>(i);
      double x = cos(theta) * radius;
      double z = sin(theta) * radius;
      points[i] = clipper::Coord_orth(x,y,z);
   }
   return points;
}

