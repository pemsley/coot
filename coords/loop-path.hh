
#ifndef LOOP_PATH_HH
#define LOOP_PATH_HH
#include <mmdb2/mmdb_manager.h>

#include "Cartesian.h"

namespace coot {

   // first is: these points need have bad CA-CA distance spots added
   // second is the vector of line segments
   std::pair<bool, std::vector<coot::CartesianPair> >
   loop_path(mmdb::Atom *start_back_2,
	     mmdb::Atom *start,
	     mmdb::Atom *end,
	     mmdb::Atom *end_plus_2,
	     unsigned int n_line_segments);

   bool is_sane_inter_residue_distance(double dist_between_residues, int res_no_delta, bool is_NA);

}

#endif // LOOP_PATH_HH
