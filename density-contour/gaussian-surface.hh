#ifndef GAUSSIAN_SURFACE_HH
#define GAUSSIAN_SURFACE_HH

#include <mmdb2/mmdb_manager.h>

#include "coot-utils/simple-mesh.hh"

namespace coot {

   class gaussian_surface_t {
      simple_mesh_t mesh;
      void using_an_nxmap(mmdb::Manager *mol);
      void using_an_xmap(mmdb::Manager *mol, const std::string &chain_id);
      void using_calc_density(mmdb::Manager *mol);
   public:
      explicit gaussian_surface_t(mmdb::Manager *mol, const std::string &chain_id);
      simple_mesh_t get_surface() const;
   };

}

#endif // GAUSSIAN_SURFACE_HH
