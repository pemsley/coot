#ifndef GAUSSIAN_SURFACE_HH
#define GAUSSIAN_SURFACE_HH

#include <mmdb2/mmdb_manager.h>

#include "coot-utils/simple-mesh.hh"

namespace coot {

   class gaussian_surface_t {
      simple_mesh_t mesh;
      void using_an_nxmap(mmdb::Manager *mol);
      void using_an_xmap(mmdb::Manager *mol, const std::string &chain_id,
                         float sigma, float contour_level, float box_radius, float grid_scale);
      void using_calc_density(mmdb::Manager *mol);
   public:
      // explicit gaussian_surface_t(mmdb::Manager *mol, const std::string &chain_id);
      explicit gaussian_surface_t(mmdb::Manager *mol, const std::string &chain_id,
                                  float sigma=4.4, float contour_level=4.0, float box_radius=5.0, float grid_scale=0.7);
      simple_mesh_t get_surface() const;
   };

}

#endif // GAUSSIAN_SURFACE_HH
