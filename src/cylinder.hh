
#ifndef COOT_CYLINDER_HH
#define COOT_CYLINDER_HH

#include <vector>
#include "generic-vertex.hh"
#include "g_triangle.hh"

class cylinder {
   void add_flat_cap(float z);
   float height;
   float base_radius;
   float top_radius;
   unsigned int n_slices;
   glm::mat4 ori;
   glm::vec3 start;
   float unstubby_rounded_cap_factor;
public:
   std::vector<s_generic_vertex> vertices;
   std::vector<g_triangle> triangle_indices_vec;
   cylinder() { unstubby_rounded_cap_factor = 1.0; } // needs to be filled
   cylinder(const std::pair<glm::vec3, glm::vec3> &cart_pair,
            float base_radius, float top_radius, float height,
            unsigned int n_slices=8, unsigned int n_stacks=2);
   void add_flat_end_cap();
   void add_flat_start_cap();
   void add_octahemisphere_end_cap();
   void add_octahemisphere_start_cap();
   void set_unstubby_rounded_cap_factor(float f) { unstubby_rounded_cap_factor = f; }
};

#endif
