/*
 * coot-utils/xmap-stats.hh
 *
 * Copyright 2008 by University of York
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


#ifndef XMAP_STATS_H
#define XMAP_STATS_H

#include "clipper/core/xmap.h"

template<class T>
class mean_and_variance { 
   
public:
   T mean; 
   T variance; 
   T range; 
   T max_density;
   T min_density;
   T bin_width;
   int histogram_max; // highest counts number (frequency), avoiding the "white line"
   std::vector<int> bins;

   std::size_t size() { return bins.size(); }
   
   // Fix this in a spare moment: 
   //
   // friend ostream& operator<<(ostream&, mean_and_variance<T>);
};

  
template<class T>
mean_and_variance<T>
map_density_distribution(const clipper::Xmap<T> &map,
			 unsigned int hist_n_bins=40,
			 bool write_output_flag=false,
			 bool ignore_pseude_zeros=false);

#endif // XMAP_STATS_H
