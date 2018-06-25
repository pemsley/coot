
#include "split-indices.hh"

// split and fill indices
void
coot::split_indices(std::vector<std::vector<unsigned int> > *indices,
		    unsigned int n_items,
		    unsigned int n_threads) {

   indices->resize(n_threads);
   unsigned int n_per_thread = n_items/n_threads + 1;
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
