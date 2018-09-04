
#include <set>

#include "radix.hh"

std::set<unsigned int> coot::unique_factors(unsigned int nx) {

   std::set<unsigned int> s;
   unsigned int r = nx;
   std::set<unsigned int> good_factors;
   good_factors.insert(2);
   good_factors.insert(3);
   good_factors.insert(5); // maybe 7 also
   bool done = false;
   while (! done) {
      done = true;
      // for (std::size_t ii=0; ii<good_factors.size(); ii++) {
      std::set<unsigned int>::const_iterator it;
      for (it=good_factors.begin(); it!=good_factors.end(); it++) {
	 if (r%*it == 0) {
	    s.insert(*it);
	    r = r / *it;
	    done = false;
	 }
      }
   }
   if (r != 1)
      s.insert(r);
   return s;
}

unsigned int
coot::suggest_radix(unsigned int nx_in, bool multiple_4_flag) {

   unsigned int nx = nx_in;
   std::set<unsigned int> good_factors;
   good_factors.insert(2);
   good_factors.insert(3);
   good_factors.insert(5); // maybe 7 also
   bool is_good = false;
   while (! is_good) {
      is_good = true;
      std::set<unsigned int> factors = unique_factors(nx);
      std::set<unsigned int>::const_iterator it;
      for (it=factors.begin(); it!=factors.end(); it++) {
	 if (good_factors.find(*it) == good_factors.end()) {
	    is_good = false;
	    break;
	 }
      }
      if (multiple_4_flag)
	 if (nx%4 != 0)
	    is_good = false;
      if (! is_good)
	 nx++;
   }
   return nx;
}
