
#ifndef COOT_CYLINDER_HH
#define COOT_CYLINDER_HH

#include <vector>
#include "generic-vertex.hh"
#include "g_triangle.hh"

class cylinder {
public:
   // std::vector<vertex_with_rotation_translation> vertices;
   std::vector<s_generic_vertex> vertices;
   std::vector<g_triangle> triangle_indices_vec;
   cylinder() {} // needs to be filled
   cylinder(const std::pair<glm::vec3, glm::vec3> &cart_pair,
            float base_radius, float top_radius, float height,
            unsigned int n_slices=8, unsigned int n_stacks=2);
};

#endif
