#ifndef POLAR_ATOMS_HH
#define POLAR_ATOMS_HH

#include "compat/coot-sysdep.h"
#include <mmdb2/mmdb_manager.h>

namespace coot {
   void buried_unsatisfied_polar_atoms(mmdb::Manager *mol);
}

#endif
