
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
