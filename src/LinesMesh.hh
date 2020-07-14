
#ifndef LINES_MESH_HH
#define LINES_MESH_HH

#include <vector>
#include <clipper/core/coords.h>

#include "generic-vertex.hh"
#include "Shader.hh"

// Currently only ambient-lit - e.g. unit cell

class LinesMesh {
   GLuint vao;
   GLuint buffer_id;
   GLuint index_buffer_id;
public:
   LinesMesh() {}
   // e.g. a box will have 8 vertices and 12 * 2 indices
   LinesMesh(const std::vector<s_generic_vertex> &vertices_in,
             const std::vector<unsigned int> &indices_in) : vertices(vertices_in), indices(indices_in) {}
   LinesMesh(const clipper::Cell &cell);
   std::vector<s_generic_vertex> vertices;
   std::vector<unsigned int> indices;
   void setup(Shader *shader_p);
   void draw(Shader *shader_p, const glm::mat4 &mvp);
   bool empty() const { return (vertices.size() == 0); }
};


#endif // LINES_MESH_HH
