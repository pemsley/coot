
// # put this in mol-utils.cc when that arrives in this branch

#ifndef MOL_UTILS_2_HH
#define MOL_UTILS_2_HH

#include <utility>
#include <map>
#include <string>

#include "mmdb2/mmdb_manager.h"

namespace coot {

   std::map<std::string, std::pair<int, int> > get_residue_number_limits(mmdb::Manager *mol);

}


#endif // MOL_UTILS_2_HH
