/* lidia-core/primes.hh
 * 
 * Copyright 2016 by Medical Research Council
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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

#ifndef PRIMES_HH
#define PRIMES_HH

#include <vector>

namespace cod {


   class primes {
      std::vector<unsigned int> prime_numbers;
   public:

      // generate primes up to and including pr_max_in
      primes(unsigned int pr_max_in) {
	 unsigned int pr_max = pr_max_in + 1;
	 std::vector<bool> pr(pr_max+1, true); // max index pr_max_in
	 pr[0] = false; pr[1] = false;
	 for (unsigned int i=2; i<pr_max; i++) {
	    for (unsigned int z=i*2; z<pr_max; z+=i) {
	       pr[z] = false;
	    }
	 }

	 unsigned int n_primes = 0;
	 for (unsigned int i=0; i<pr_max; i++) {
	    if (pr[i])
	       n_primes++;
	 }

	 prime_numbers.reserve(n_primes);
	 
	 for (unsigned int i=0; i<pr_max; i++)
	    if (pr[i]) prime_numbers.push_back(i);
	 
      }
      std::vector<unsigned int> get_primes() const { return prime_numbers; }
   };
}

#endif // PRIMES_HH
