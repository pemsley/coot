
#include <vector>

namespace coot {

   // split and fill indices
   void split_indices(std::vector<std::vector<unsigned int> > *indices,
		      unsigned int n_items,
		      unsigned int n_threads);

   std::vector<std::pair<unsigned int, unsigned int> >
   atom_index_ranges(unsigned int n_atoms, unsigned int n_threads);

}
