

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

  // Fix this in a spare moment:  To complex (or I am too stupid
  // to do it while rushing now). 
  // 
  // friend ostream& operator<<(ostream&, mean_and_variance<T>);


};

  
template<class T>
mean_and_variance<T>
map_density_distribution(const clipper::Xmap<T> &map, short int write_output_flag); 

#endif // XMAP_STATS_H
