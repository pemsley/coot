

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
   
   // Fix this in a spare moment: 
   //
   // friend ostream& operator<<(ostream&, mean_and_variance<T>);
};

  
template<class T>
mean_and_variance<T>
map_density_distribution(const clipper::Xmap<T> &map,
			 bool write_output_flag=false,
			 bool ignore_pseude_zeros=false);

#endif // XMAP_STATS_H
