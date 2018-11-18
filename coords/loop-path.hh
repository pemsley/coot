
#ifndef LOOP_PATH_HH
#define LOOP_PATH_HH
#include <mmdb2/mmdb_manager.h>

#include "Cartesian.h"

namespace coot {

   std::vector<coot::CartesianPair>
   loop_path(mmdb::Atom *start_back_2,
	     mmdb::Atom *start,
	     mmdb::Atom *end,
	     mmdb::Atom *end_plus_2,
	     unsigned int n_line_segments);

}

#endif // LOOP_PATH_HH
