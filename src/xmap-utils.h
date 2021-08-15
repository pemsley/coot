#ifndef XMAP_UTILS_H
#define XMAP_UTILS_H

#include <string>
#include <clipper/core/xmap.h>

// I suppose that this could be templated to (i.e. an int is also useful)
// 
float nearest_step(float val, float step); 

int 
nint(float val);


class Xmap_plus_bits {
 public:
   clipper::Xmap<float> xmap;
   std::string name;
   float min, max, mean, std_dev; 
}; 


#endif // XMAP_UTILS_H
