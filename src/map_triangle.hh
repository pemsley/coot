#ifndef MAP_TRIANGLE_HH
#define MAP_TRIANGLE_HH

#include "glm/glm.hpp"
#include "g_triangle.hh"

class map_triangle_t : public g_triangle {
public:
   glm::vec3 mid_point;
   float back_front_projection_distance;
   map_triangle_t(const unsigned int &a0,
                  const unsigned int &a1,
                  const unsigned int &a2,
                  const glm::vec3 &mp) : g_triangle(a0, a1, a2), mid_point(mp), back_front_projection_distance(0.0f) {}
   map_triangle_t(const g_triangle &gt, const glm::vec3 &mp) : g_triangle(gt), mid_point(mp), back_front_projection_distance(0.0f) {}
};


#endif // MAP_TRIANGLE_HH
