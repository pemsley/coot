
#include <iostream>
#include "split-indices.hh"

// split and fill indices
void
coot::split_indices(std::vector<std::vector<unsigned int> > *indices,
		    unsigned int n_items,
		    unsigned int n_threads) {

   indices->resize(n_threads);
   unsigned int n_per_thread = n_items/n_threads + 1;
   if ((n_items % n_threads) == 0)
      n_per_thread = n_items/n_threads;
   for (std::size_t i=0; i<n_threads; i++)
      indices->at(i).reserve(n_per_thread);

   unsigned int idx = 0;
   for (unsigned int i=0; i<n_items; i++) {
      indices->at(idx).push_back(i);
      idx++;
      if (idx == n_threads)
	 idx = 0;
   }
}


std::vector<std::pair<unsigned int, unsigned int> >
coot::atom_index_ranges(unsigned int n_atoms, unsigned int n_threads) {

   std::vector<std::pair<unsigned int, unsigned int> > v;

   if (n_threads == 0) n_threads = 1;

   unsigned int npt = n_atoms/n_threads;

   if ((n_atoms % n_threads) == 0) {

      // divides exactly

      for (std::size_t i=0; i<n_threads; i++) {
	 unsigned int v1 = i*npt;
	 unsigned int v2 = (i+1)*npt;
	 std::pair<unsigned int, unsigned int> p(v1, v2);
	 v.push_back(p);
      }

   } else {

      if (n_threads == 1) {

	 unsigned int v1 = 0;
	 unsigned int v2 = n_atoms;
	 std::pair<unsigned int, unsigned int> p(v1, v2);
	 v.push_back(p);

      } else {

	 if (n_threads == 2) {

	    // we know that it doesn't divide exactly

	    unsigned int v1 = 0;
	    unsigned int v2 = n_atoms/2;
	    std::pair<unsigned int, unsigned int> p0(v1, v2);
	    v.push_back(p0);
	    v1 = n_atoms/2;
	    v2 = n_atoms;
	    std::pair<unsigned int, unsigned int> p1(v1, v2);
	    v.push_back(p1);

	 } else {

	    // the final thread picks up the scraps (if any) (consider 25 atoms, 6 threads)
            // 20210812-PE  also consider 10 atoms 8 threads
            // Also  what about 4 atoms and 10 threds?
            // I need a test function for this.

	    npt = n_atoms/(n_threads-1);

	    if (false)
	       std::cout << "debug:: atom_index_ranges() n_atoms " << n_atoms
			 << " n_threads " << n_threads << " n_per_thread " << npt << std::endl;

	    for (std::size_t i=0; i<(n_threads-1); i++) {
	       unsigned int v1 = i*npt;
	       unsigned int v2 = (i+1)*npt;
	       std::pair<unsigned int, unsigned int> p(v1, v2);
	       v.push_back(p);
	    }
	    // scraps for the last thread
	    unsigned int v1 = (n_threads-1)*npt;
	    unsigned int v2 = n_atoms;
	    if (v2 > v1) {
	       std::pair<unsigned int, unsigned int> p(v1, v2);
	       v.push_back(p);
	    }
	 }
      }
   }

   return v;

}
