/* coot-utils/xmap-stats.cc
 * 
 * Copyright 2002, The University of York
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
map_density_distribution(const clipper::Xmap<T> &map, bool write_output_flag, bool ignore_pseude_zeros) { 


   mean_and_variance<T> mv;  // returned object

   double min = 1e10, max = -1e10, sum = 0.0, sum_sq = 0;
   double v; // optimised?
   double n_point = 0.0;
   
   clipper::Xmap_base::Map_reference_index ix;
   for (ix=map.first(); !ix.last(); ix.next()) {

      n_point += 1.0; 
      v = double (map[ix]); 

      if (v < min) min = v; 
      if (v > max) max = v;

      sum += v;
      sum_sq += v*v; 

   }

   float mean = float( sum/n_point );

   float var = float( (n_point*sum_sq - sum*sum) / (n_point*n_point) );
   float range = float( max - min );
   float inv_range = (range>0.0) ? (1.0/range) : (1.0);


   mv.mean = mean;
   mv.variance = var; 
   mv.range = range;
   mv.max_density = max;
   mv.min_density = min;


   if (write_output_flag) { 
      std::cout << "Map statistics: mean: " << mean << " st.d: " << sqrt(var)
	   << std::endl;
      std::cout <<  "Map statistics: min: " << min << ", max: " << max << std::endl;

      int nbins = 40;
      int bin_no; 
      std::vector<int> bin(nbins+1, 0);

      for (ix=map.first(); !ix.last(); ix.next()) {
	 bin_no = int (nbins*(map[ix] - min)*inv_range);
	 bin[bin_no]++; 
      }

      for (int ii=0; ii<= nbins; ii++) {
	 std::cout.width(10);
	 std::cout << std::right << ((float) (ii) + 0.5)*range/float(nbins) + min << "    "
		   << std::right << bin[ii] << std::endl;
      }
   }

   if (ignore_pseude_zeros) {
      int nbins = 10000;
      int bin_no; 
      long n = 0;
      std::vector<int> bin(nbins+1, 0);

      for (ix=map.first(); !ix.last(); ix.next()) {
	 bin_no = int (nbins*(map[ix] - min)*inv_range);
	 bin[bin_no]++; 
	 n++;
      }
      double sum = 0;
      double sum_sq = 0;
      unsigned int ibin_max_counts = 0;
      for (unsigned int i=0; i<bin.size(); i++) {
	 if (bin[i] > bin[ibin_max_counts])
	    ibin_max_counts = i;
      }
      std::cout << "INFO:: ignoring " << bin[ibin_max_counts]
		<< " of " << n << " counts ( = " << std::setprecision(2)
		<< 100*float(bin[ibin_max_counts])/float(n) << "%)"
		<< " with values around "
		<< std::setprecision(4)
		<< (static_cast<double>(ibin_max_counts) + 0.5) * range /double(nbins) + min
		<< " from bin " << ibin_max_counts << " of "
		<< nbins<< std::endl;
      for (unsigned int i=0; i<bin.size(); i++) {
	 if (i != ibin_max_counts) { 
	    double v = (static_cast<double>(i) + 0.5) * range /double(nbins) + min;
	    sum += v * bin[i];
	    sum_sq += v * v * bin[i];
	 }
      }

      if (0) { 
	 std::cout << "=== n is " << n << std::endl;
	 std::cout << "=== previous mean "     << mv.mean     << std::endl;
	 std::cout << "=== previous variance " << mv.variance << std::endl;
      } 
      mv.mean = sum/double(n);
      mv.variance = sum_sq/double(n) - mv.mean * mv.mean;
      if (mv.variance < 0)
	 mv.variance = 0;
   } 
   return mv; 
}

// instantiate that 
template mean_and_variance<float> map_density_distribution(const clipper::Xmap<float> &map,
							   bool write_output_flag,
							   bool ignore_pseude_zeros);

// template mean_and_variance<int> map_density_distribution(const clipper::Xmap<int> &map, short int write_output_flag); 


// ostream& operator<<(ostream& s, mean_and_variance<T> mv) { 

//    s << "mean: " << mv.mean << ", sigma: " << sqrt(mv.variance); 
//    return s; 
// }

