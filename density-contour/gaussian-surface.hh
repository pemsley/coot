#ifndef GAUSSIAN_SURFACE_HH
#define GAUSSIAN_SURFACE_HH

#include <mmdb2/mmdb_manager.h>
#include <clipper/core/xmap.h>

#include "coot-utils/simple-mesh.hh"

namespace coot {

   class gaussian_surface_t {
      simple_mesh_t mesh;
      void using_an_nxmap(mmdb::Manager *mol);
      void using_an_xmap(mmdb::Manager *mol, const std::string &chain_id,
                         float sigma, float contour_level, float box_radius, float grid_scale,
                         float b_factor);
      void using_calc_density(mmdb::Manager *mol);
      void normals_from_function_gradient(const clipper::Xmap<float> &xmap,
                                          const glm::vec3 &cb); // changes mesh normals
   public:
      // explicit gaussian_surface_t(mmdb::Manager *mol, const std::string &chain_id);
      explicit gaussian_surface_t(mmdb::Manager *mol, const std::string &chain_id,
                                  float sigma=4.4, float contour_level=4.0, float box_radius=5.0,
                                  float grid_scale=0.7, float fft_b_factor=100.0);
      simple_mesh_t get_surface() const;
   };

}

#endif // GAUSSIAN_SURFACE_HH
