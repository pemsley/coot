#ifndef MOLECULAR_TRIANGLES_MESH_HH
#define MOLECULAR_TRIANGLES_MESH_HH

#include <vector>
#include <string>
#include "generic-vertex.hh"
#include "g_triangle.hh"

class molecular_triangles_mesh_t {
public:
   molecular_triangles_mesh_t() {}
   molecular_triangles_mesh_t(const std::vector<s_generic_vertex> &vertices_in,
                              const std::vector<g_triangle> &triangles_in,
                              const std::string &name_in) :
      vertices(vertices_in), triangles(triangles_in), name(name_in) { type_index = 0; }
   std::vector<s_generic_vertex> vertices;
   std::vector<g_triangle> triangles;
   std::string name;
   unsigned int type_index;
   void add_to_mesh(const std::vector<s_generic_vertex> &gv,
                    const std::vector<g_triangle> &tris);
   void add_to_mesh(const molecular_triangles_mesh_t &new_mesh);
   bool empty() const { return vertices.empty(); }

};


#endif // MOLECULAR_TRIANGLES_MESH_HH
