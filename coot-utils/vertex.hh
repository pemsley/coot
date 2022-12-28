#ifndef VERTEX_HH
#define VERTEX_HH

#include <glm/glm.hpp>

// for standard objects at the origin - typically used in instancing
namespace coot {

   namespace api {

      class vn_vertex {
      public:
         glm::vec3 pos;
         glm::vec3 normal; // normalized on input
         vn_vertex(const glm::vec3 &pos_in,
                   const glm::vec3 &norm_in) :
            pos(pos_in), normal(norm_in) {}
         vn_vertex() {}
      };

      class vnc_vertex {
      public:
         glm::vec3 pos;
         glm::vec3 normal; // normalized on input
         glm::vec4 color;  // or colour?
         vnc_vertex(const glm::vec3 &pos_in,
                    const glm::vec3 &norm_in,
                    const glm::vec4 &col_in) : pos(pos_in), normal(norm_in), color(col_in) {}
         explicit vnc_vertex(const vn_vertex &vn) :
            pos(vn.pos), normal(vn.normal), color(glm::vec4(0.5, 0.5, 0.5, 1.0)) {}
         vnc_vertex() {}
      };
   }
}


#endif // VERTEX_HH
