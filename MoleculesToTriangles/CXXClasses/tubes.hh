#include <mmdb2/mmdb_manager.h>
#include "coot-utils/simple-mesh.hh"

coot::simple_mesh_t
make_tubes_representation(mmdb::Manager *mol,
                          const std::string &atom_selection_str,
                          const std::string &colour_scheme,
                          int secondaryStructureUsageFlag);
