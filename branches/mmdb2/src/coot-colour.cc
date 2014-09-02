
#include <vector>
#include <iostream>

#include "coot-colour.hh"

namespace coot { 
   std::ostream& operator<<(std::ostream &s, colour_t col) {
      s << col.col[0] << " " << col.col[1] << " " << col.col[2];
      return s;
   }
}
