
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


std::vector<std::pair<unsigned int, unsigned int> >
coot::atom_index_ranges(unsigned int n_atoms, unsigned int n_threads) {

   std::vector<std::pair<unsigned int, unsigned int> > v;

   unsigned int npt = n_atoms/n_threads + 1;
   for (std::size_t i=0; i<n_threads; i++) {
      unsigned int v1 = i*npt;
      unsigned int v2 = (i+1)*npt;
      if (i == (n_threads - 1))
	 v2 = n_atoms;
      std::pair<unsigned int, unsigned int> p(v1, v2);
      v.push_back(p);
   }

   return v;

}
