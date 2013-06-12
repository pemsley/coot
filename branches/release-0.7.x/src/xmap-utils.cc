
#include <iostream>

#include "clipper/core/xmap.h"

#include "xmap-utils.h"



float nearest_step(float val, float step) { 

   return step*nint(val/step); 

} 


int 
nint(float val) { 

  if (val > 0) {
    return int (val + 0.5); 
  } else { 
    return int (val - 0.5); 
  }

} 


