
#include <iostream>
#include <vector>
#include <algorithm>
#include "align-utils.hh"
#ifdef HAVE_CXX11
#include <algorithm>
#endif

// strings with - in them
std::string
coot::alignment_matches(const std::string &aligned,
			const std::string &target) {

   unsigned int la = aligned.size();
   unsigned int lt = target.size();

   if (lt<la) la = lt;

   std::string matches(la, ' ');

   for (unsigned int i=0; i<lt; i++) {
      if (aligned[i] == target[i])
	 matches[i] = '|';
   }

#ifdef HAVE_CXX11
   std::vector<char> aromatics  = {'F', 'W', 'Y'}; // not H interestingly
   std::vector<char> aliphatics = {'A', 'G', 'L', 'V', 'I'};
   std::vector<char> hydroxyls  = {'T', 'Y', 'S', 'C', 'M'}; // and sulfhydryl
   std::vector<char> acidics    = {'D', 'L', 'H'};
   std::vector<char> basics     = {'R', 'E', 'N', 'Q'};
   // std::vector<char> cyclics    = {'P'}; won't make non-P matches

   std::vector<std::vector<char> > types(5);
   types[0] = aromatics;
   types[1] = aliphatics;
   types[2] = hydroxyls;
   types[3] = acidics;
   types[4] = basics;

   for (unsigned int i=0; i<lt; i++) {
      char ca = aligned[i];
      char ta = target[i];
      if (ca != ta) { // we would found it above if this were true
	 if (ca != '-') {
	    if (ta != '-') {
	       for (std::size_t j=0; j<types.size(); j++) {
		  const std::vector<char> &codes = types[j];
		  // are both ca and ta in codes?
		  std::vector<char>::const_iterator it_1;
		  std::vector<char>::const_iterator it_2;
		  it_1 = std::find(codes.begin(), codes.end(), ca);
		  it_2 = std::find(codes.begin(), codes.end(), ta);
		  if (it_1 != codes.end()) {
		     if (it_2 != codes.end()) {
			matches[i] = ':';
			break;
		     }
		  }
	       }
	    }
	 }
      }
   }

#endif

   return matches;
}
