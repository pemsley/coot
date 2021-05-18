/* coot-utils/coot-coord-rama.hh
 * 
 * Copyright 2014 by Medical Research Council
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

#ifndef COOT_DENSITY_STATS_HH
#define COOT_DENSITY_STATS_HH

#include <vector>
#include <cmath>

namespace coot {

   enum map_stats_t {
      SIMPLE,
      WITH_KOLMOGOROV_SMIRNOV_DIFFERENCE_MAP_TEST};

   namespace util {

      class density_stats_info_t {
	 
      public:
	 // these are doubles because we can potentially do lots (millions) of tiny adds.
	 double n;
	 double sum_sq; // sum of squared elements
	 double sum;
	 double sum_weight;

	 density_stats_info_t() {
	    n = 0.0;
	    sum = 0.0;
	    sum_sq = 0.0;
	    sum_weight = 0.0;
	 }
	 void add(const double &v) {
	    n += 1.0;
	    sum += v;
	    sum_sq += v*v;
	    sum_weight += 1.0;
	 } 
	 void add(const double &v, const double &weight) {
	    n += weight;
	    sum += weight*v;
	    sum_sq += weight*v*v;
	    sum_weight += 1.0;
	 }
	 std::pair<double, double> mean_and_variance() const {
	    double mean = -1;
	    double var  = -1;
	    if (n > 0) {
	       mean = sum/sum_weight;
	       var = sum_sq/sum_weight - mean*mean;
	    }
	    return std::pair<double, double> (mean, var);
	 }
      };


      class density_correlation_stats_info_t {
      public:
	 double n;
	 double sum_xy;
	 double sum_sqrd_x;
	 double sum_sqrd_y;
	 double sum_x;
	 double sum_y;
	 // for doing KS tests (against normal distribution) , we want
	 // all the density samples
	 std::vector<double> density_values;
	 density_correlation_stats_info_t() {
	    n = 0;
	    sum_xy = 0;
	    sum_sqrd_x = 0;
	    sum_sqrd_y = 0;
	    sum_x = 0;
	    sum_y = 0;
	 }
	 density_correlation_stats_info_t(double n_in,
					  double sum_xy_in,
					  double sum_sqrd_x_in,
					  double sum_sqrd_y_in,
					  double sum_x_in,
					  double sum_y_in) {
	    n = n_in;
	    sum_xy = sum_xy_in;
	    sum_sqrd_x = sum_sqrd_x_in;
	    sum_sqrd_y = sum_sqrd_y_in;
	    sum_x = sum_x_in;
	    sum_y = sum_y_in;
	 }
	 double var_x() const {
	    double mean_x = sum_x/n;
	    return (sum_sqrd_x/n - mean_x * mean_x);
	 } 
	 double var_y() const {
	    double mean_y = sum_y/n;
	    return (sum_sqrd_y/n - mean_y * mean_y);
	 }
         void add(const double &x, const double &y) {
            n += 1;
            sum_x  += x;
            sum_y  += y;
            sum_xy += x * y;
            sum_sqrd_x += x * x;
            sum_sqrd_y += y * y;
         }
	 double correlation() const {
	    double top = n * sum_xy - sum_x * sum_y;
	    double b_1 = n * sum_sqrd_x - sum_x * sum_x;
	    double b_2 = n * sum_sqrd_y - sum_y * sum_y;
	    if (b_1 < 0) b_1 = 0;
	    if (b_2 < 0) b_2 = 0;
	    double c = top/(std::sqrt(b_1) * std::sqrt(b_2));
	    return c;
	 }
      }; 

   }
}

#endif // COOT_DENSITY_STATS_HH
