
#include "simple-mesh.hh"


void
coot::simple_mesh_t::translate(const glm::vec3 &t) {

   for (auto &item : vertices) {
      item.pos += t;
   }

}


void
coot::simple_mesh_t::add_submesh(const simple_mesh_t &submesh) {

   unsigned int idx_base = vertices.size();
   unsigned int idx_base_tri = triangles.size();
   vertices.insert(vertices.end(), submesh.vertices.begin(), submesh.vertices.end());
   triangles.insert(triangles.end(), submesh.triangles.begin(), submesh.triangles.end());
   for (unsigned int i=idx_base_tri; i<triangles.size(); i++)
      triangles[i].rebase(idx_base);
   
}
