

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
