
#include "simple-mesh.hh"


void
coot::simple_mesh_t::translate(const glm::vec3 &t) {

   for (auto &item : vertices) {
      item.pos += t;
   }

}
