/* coot-utils/xmap-stats.cc
 * 
 * Copyright 2002, The University of York
 * Copyright 2014, 2015 by Medical Research Council
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

#include <iomanip>
#include "xmap-stats.hh"



template <class T>
mean_and_variance<T>
map_density_distribution(const clipper::Xmap<T> &map,
			 unsigned int n_bins,     // default 40
			 bool write_output_flag,  // default false 
			 bool ignore_pseudo_zeros // default false
			 ) {

   mean_and_variance<T> mv;  // returned object

   double density_min = 1e10, density_max = -1e10, sum = 0.0, sum_sq = 0;
   double n_point = 0.0;
   T rho;

   clipper::Xmap_base::Map_reference_index ix;
   for (ix=map.first(); !ix.last(); ix.next()) {

      n_point += 1.0;
      rho = map[ix];
      if (! clipper::Util::is_nan(rho)) {
	 double v = double (rho);
	 if (v < density_min) density_min = v;
	 if (v > density_max) density_max = v;
	 sum += v;
	 sum_sq += v*v;
      }
   }

   float mean = sum/n_point;

   float var = float( (n_point*sum_sq - sum*sum) / (n_point*n_point) );
   float range = float(density_max - density_min);
   float inv_range = (range>0.0) ? (1.0/range) : (1.0);

   mv.mean = mean;
   mv.variance = var; 
   mv.range = range;
   mv.max_density = density_max;
   mv.min_density = density_min;

   std::vector<int> bin(n_bins+1, 0);
   mv.bin_width = range/double(n_bins);

   for (ix=map.first(); !ix.last(); ix.next()) {
      rho = map[ix];
      if (! clipper::Util::is_nan(rho)) { 
	 int bin_no = int (n_bins*(rho - density_min)*inv_range);
	 // std::cout << "bin_no: " << bin_no << " from " << nbins << "*(" << map[ix] << "- " << min << ") * "
	 // << inv_range << std::endl;
	 bin[bin_no]++;
	 // std::cout << "bin[" << bin_no << "] is now " << bin[bin_no] << std::endl;
      }
   }

   unsigned int ibin_max_counts = 0;
   for (unsigned int ii=0; ii <= n_bins; ii++) {
      if (bin[ii] > bin[ibin_max_counts])
         ibin_max_counts = ii;
   }
   mv.histogram_max = bin[ibin_max_counts];
   if (write_output_flag) {
      for (unsigned int ii=0; ii <= n_bins; ii++) {
	 std::cout.width(10);
	 std::cout << std::right << ((float) (ii) + 0.5)*range/float(n_bins) + density_min << "    "
		   << std::right << bin[ii] << std::endl;
      }
   }
   unsigned int ibin_second_highest_counts = 0;
   for (unsigned int i=0; i<bin.size(); i++) {
      if (i != ibin_max_counts)
         if (bin[i] > bin[ibin_second_highest_counts])
            ibin_second_highest_counts = i;
   }
   mv.histogram_max = bin[ibin_second_highest_counts] * 1.25;
   mv.bins = bin;

   std::cout << "Pre-filter Map statistics: mean: " << mean << " st.d: " << sqrt(var) << std::endl;
   std::cout << "Pre-filter Map statistics:  min: " << density_min << " max: " << density_max << std::endl;

   if (ignore_pseudo_zeros) {

      // Print to screen and fill bins in mv.
      
      int nbins_filter = 10000;
      long n = 0;
      bin = std::vector<int>(nbins_filter+1, 0);

      for (ix=map.first(); !ix.last(); ix.next()) {
	 rho = map[ix];
	 if (! clipper::Util::is_nan(rho)) {
	    int bin_no = static_cast<int>(nbins_filter*(rho - density_min)*inv_range);
	    bin[bin_no]++;
	    n++;
	 }
      }
      ibin_max_counts = 0;
      for (unsigned int i=0; i<bin.size(); i++) {
	 if (bin[i] > bin[ibin_max_counts])
	    ibin_max_counts = i;
      }

      int n_remainder = n - bin[ibin_max_counts];

      std::cout << "INFO:: n grid points:             " << n << std::endl;
      std::cout << "INFO:: mean before filtering:     " << mv.mean     << std::endl;
      std::cout << "INFO:: variance before filtering: " << mv.variance << std::endl;

      double di = static_cast<double>(ibin_max_counts);
      double average_density_for_ibin_max_counts      = (di + 0.5) * range/double(nbins_filter) + density_min;
      double lower_bound_density_for_ibin_max_counts  =  di        * range/double(nbins_filter) + density_min;
      double uppoer_bound_density_for_ibin_max_counts = (di+1.0)   * range/double(nbins_filter) + density_min;

      std::cout << "INFO:: filter by ignoring " << bin[ibin_max_counts]
		<< " of " << n << " counts ( = " << std::setprecision(4)
		<< 100*float(bin[ibin_max_counts])/float(n) << "%)"
		<< " with values around "
		<< std::setprecision(4)
		<< average_density_for_ibin_max_counts
                << " bounds " << lower_bound_density_for_ibin_max_counts << " "
                << uppoer_bound_density_for_ibin_max_counts
		<< " from bin-number " << ibin_max_counts << " of "
		<< nbins_filter<< std::endl;
      sum = 0;
      sum_sq = 0;
      for (unsigned int i=0; i<bin.size(); i++) {
	 if (i != ibin_max_counts) { 
	    double v = (static_cast<double>(i) + 0.5) * range/static_cast<double>(nbins_filter) + density_min;
            int counts = bin[i];
            double weight = static_cast<double>(counts);
	    sum += v * weight;
	    sum_sq += v * v * weight;
	 }
      }

      mv.mean = sum/static_cast<float>(n_remainder);
      mv.variance = sum_sq/double(n_remainder) - mv.mean * mv.mean;
      if (mv.variance < 0.0) mv.variance = 0.0;

      std::cout << "Post-filter Map statistics: mean: " << mv.mean << " st.d: " << sqrt(mv.variance) << std::endl;
      std::cout << "Post-filter Map statistics: min: " << density_min << " max: " << density_max << std::endl;

   }
   return mv; 
}

// instantiate that 
template mean_and_variance<float> map_density_distribution(const clipper::Xmap<float> &map,
                                                           unsigned int n_bins,
                                                           bool write_output_flag,
                                                           bool ignore_pseude_zeros);

// template mean_and_variance<int> map_density_distribution(const clipper::Xmap<int> &map, short int write_output_flag); 


// ostream& operator<<(ostream& s, mean_and_variance<T> mv) { 

//    s << "mean: " << mv.mean << ", sigma: " << sqrt(mv.variance); 
//    return s; 
// }

