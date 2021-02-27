
#ifndef RADIX_HH
#define RADIX_HH

#include <set>

namespace coot {

   unsigned int suggest_radix(unsigned int nx_in, bool multiple_4_flag=false);

   std::set<unsigned int> unique_factors(unsigned int nx);

}


#endif // RADIX_HH

