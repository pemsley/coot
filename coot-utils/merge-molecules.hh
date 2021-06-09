#ifndef MERGE_MOLECULES_HH
#define MERGE_MOLECULES_HH

#include <vector>
#include <mmdb2/mmdb_manager.h>

namespace coot {
   // only copy the first chain of the mol_otheres
   void merge_molecules(mmdb::Manager *mol, std::vector<mmdb::Manager *> mol_others);
}



#endif // MERGE_MOLECULES_HH
