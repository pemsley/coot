
#ifndef GENERIC_VERTEX
#define GENERIC_VERTEX

#define GLM_ENABLE_EXPERIMENTAL // needed?
#include <glm/glm.hpp>

class generic_vertex {
public:
   glm::mat3 model_rotation_matrix; // orientation
   glm::vec3 model_translation; // the coordinates of the first atom of the bond
   glm::vec3 pos;
   glm::vec3 normal; // normalized when set
   glm::vec4 colour;
   generic_vertex(const glm::vec3 &p, const glm::vec3 &n, const glm::vec4 &c) : pos(p), normal(n), colour(c) {}
   generic_vertex() {}
};

// class graphical_triangle {
// public:
//    graphical_triangle(const unsigned int &a0,
//                       const unsigned int &a1,
//                       const unsigned int &a2) {
//       point_id[0] = a0;
//       point_id[1] = a1;
//       point_id[2] = a2;
//    }
//    unsigned int point_id[3];
// };

// simple version of above 
class s_generic_vertex {
public:
   glm::vec3 pos;
   glm::vec3 normal; // normalized on input
   glm::vec4 color;
   s_generic_vertex(const glm::vec3 pos_in,
                    const glm::vec3 norm_in,
                    const glm::vec4 col_in) : pos(pos_in), normal(norm_in), color(col_in)  {}
   s_generic_vertex() {}
};


#endif
