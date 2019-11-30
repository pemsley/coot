/* analysis/stats.cc
 * 
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

#ifndef INCLUDE_STATS_HH
#define INCLUDE_STATS_HH

#include <vector>
#include <utility>
#include <algorithm>

#include <math.h>

#include <gsl/gsl_sf_erf.h>

namespace coot {

   namespace stats {

      // 1-d data
      class single {
	 
	 // double cached_kurtosis;
	 // bool have_cached_kurtosis;

      public:
	 std::vector<double> v;

	 single() {
	    // have_cached_kurtosis = false;
	 }
	 single(const std::vector<double> &v_in) {
	    v = v_in;
	    // have_cached_kurtosis = false;
	 }
	 unsigned int size() const { return v.size(); }
	 void add(const double &a) {
	    v.push_back(a);
	    // have_cached_kurtosis = false;
	 }
	 
	 double mean() const {
	    double m = 0;
	    double sum = 0;
	    if (v.size() ) { 
	       for (unsigned int i=0; i<v.size(); i++)
		  sum += v[i];
	       m = sum/double(v.size());
	    }
	    return m;
	 }

	 double variance() const {
	    double var = 0;
	    double sum = 0;
	    double sum_sq = 0;
	    if (v.size() ) { 
	       for (unsigned int i=0; i<v.size(); i++) { 
		  sum += v[i];
		  sum_sq += v[i] * v[i];
	       }
	       double m = sum/double(v.size());
	       var = sum_sq/double(v.size()) - m*m;
	    }
	    if (var < 0) var = 0; // numerical stability
	    return var;
	 }

	 double skew() const {
	    double skew = 0;
	    double sum = 0;
	    double sum_cubed = 0;
	    double m = mean();
	    double var = variance();
	    double sigma = sqrt(var);
	    double s3 = sigma * sigma * sigma;
	    if (v.size() > 0) {
	       for (unsigned int i=0; i<v.size(); i++) {
		  double delta = v[i] - m;
		  sum_cubed += delta * delta * delta;
	       }
	       skew = (sum_cubed/double(v.size()))/s3;
	    }
	    return skew;
	 }

	 double kurtosis() const {

	    // recall kurtosis, $k$ of $N$ observations:
	    // k = \frac{\Sigma(x_i - \mu)^4} {N \sigma^4} - 3    
	    // (x_i - \mu)^4 = x_i^4 + 4x_i^3\mu + 6x_i^2\mu^2 + 4x_i\mu^3 + \mu^4

	    // Can't enable this! A compiler bug (maybe) is apparent when sorting
	    // (bonds_vec_k_sorter()).
	    // g++ (Ubuntu 4.4.3-4ubuntu5.1) 4.4.3
	    // 
	    // if (have_cached_kurtosis)
	    // return cached_kurtosis;
	    
	    double k = -999;
	    if (v.size() ) {

	       double sum_to_the_4 = 0;
	       double m = mean();
	       double var = variance();

	       if (var > 0) {
		  for (unsigned int i=0; i<v.size(); i++) { 
		     double t = v[i] - m;
		     sum_to_the_4 += t * t * t * t;
		  }
		  k = sum_to_the_4/(double(v.size()) * var * var);
		  // cached_kurtosis = k;
		  // have_cached_kurtosis = true;
	       }
	    }
	    return k;
	 }

	 std::pair<double, double> median_and_iqr() const {

	    std::vector<double> vv = v;
	    std::sort(vv.begin(), vv.end());
	    int n = vv.size();

	    int idx_q1 = static_cast<int>(0.25 * static_cast<double>(n));
	    int idx_q2 = static_cast<int>(0.50 * static_cast<double>(n));
	    int idx_q3 = static_cast<int>(0.75 * static_cast<double>(n));
	    double iqr = vv[idx_q3] - vv[idx_q1];
	    double m = vv[idx_q2];
	    if (n%2 == 0) {
	       int idx_q2a = idx_q2 + 1;
	       if (idx_q2a < n)
		  m = (m + vv[idx_q2a]) * 0.5;
	    }
	    return std::pair<double, double> (m,iqr);
	 }
      };

      class pnorm {

	 // return a cumulative probability for this number of standard deviations from the mean
	 // e.g. return for 0, return 0.5 and -1 return 0.1586

	 void init() { }
      public:
	 pnorm() { init(); }
	 double erf(const double &z) const; // public for testing.
	 double get(const double &x) const {
	    // return 0.5 * (1 + erf(x/sqrt(2.0)));
	    return 0.5 * (1 + gsl_sf_erf(x/sqrt(2.0)));
	 } 
      };

      // 20150807-PE
      //
      double get_kolmogorov_smirnov_vs_normal(const std::vector<double> &v1,
					      const double &reference_mean,
					      const double &reference_sd); 
   }
}

#endif // INCLUDE_STATS_HH
