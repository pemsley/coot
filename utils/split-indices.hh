
#include <vector>

namespace coot {

   // split and fill indices
   void split_indices(std::vector<std::vector<unsigned int> > *indices,
		      unsigned int n_items,
		      unsigned int n_threads);

}
