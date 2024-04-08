/*
 * utils/align-utils.cc
 *
 * Copyright 2009 by University of Oxford
 * Author: Paul Emsley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


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
