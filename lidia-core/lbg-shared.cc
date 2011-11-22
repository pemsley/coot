
#include <iostream>"

#include "lbg-shared.hh"

std::ostream&
coot::operator<< (std::ostream& s, const coot::bash_distance_t &bd) {

   if (bd.unlimited()) { 
      s << "unlimited"; // C&L.
   } else {
      s << bd.dist;
   } 
   return s;
} 
