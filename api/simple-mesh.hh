#ifndef SIMPLE_MESH_HH
#define SIMPLE_MESH_HH

#include <vector>
#include "vertex.hh"
#include "g_triangle.hh"

namespace coot {

   class simple_mesh_t {
   public:
      std::vector<api::vnc_vertex> vertices;
      std::vector<g_triangle> triangles;
      simple_mesh_t() {}
      simple_mesh_t(const std::vector<api::vnc_vertex> &vertices_in,
                    const std::vector<g_triangle> &triangles_in) : vertices(vertices_in), triangles(triangles_in) {}
      void translate(const glm::vec3 &t);
   };
}

#endif // SIMPLE_MESH_HH
