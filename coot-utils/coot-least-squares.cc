/* coot-utils/coot-least-squares.cc
 * 
 * Copyright 2010 by The University of Oxford
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

#include <iostream>

#include "coot-least-squares.hh"

coot::least_squares_fit::least_squares_fit(const std::vector<std::pair<double, double> > &data) {

   a_ = 0;
   b_ = 0;
   
   double sum_x = 0;
   double sum_y = 0;
   
   double sum_x_sqrd = 0;
   double sum_y_sqrd = 0;

   double sum_xy = 0.0;

   if (data.size() > 2) {
      double n = data.size();
      for (unsigned int i=0; i<data.size(); i++) {
	 const double &x = data[i].first;
	 const double &y = data[i].second;
	 sum_x += x;
	 sum_y += y;
	 sum_x_sqrd += x * x;
	 sum_y_sqrd += y * y;
	 sum_xy += x * y;
      }

      double x_mean = sum_x/n;
      double y_mean = sum_y/n;
      // double cov_xy = sum_xy/n;

      double ss_xx = sum_x_sqrd - n * x_mean * x_mean;
      double ss_xy = sum_xy     - n * x_mean * y_mean;

      if (false) { 
	 std::cout << "n: " << n << std::endl;
	 std::cout << "mean x " << x_mean << std::endl;
	 std::cout << "mean y " << y_mean << std::endl;
	 std::cout << "ss_xx: " << ss_xx << std::endl;
	 std::cout << "ss_xy: " << ss_xy << std::endl;
      }

      b_ = ss_xy / ss_xx;                 // m
      a_ = y_mean - b_ * x_mean;          // c
   }
} 
