
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
   void init();
   void make_vertices_for_pulse(float radius);
   glm::vec3 central_position;
   std::string name;

public:
   LinesMesh() { init(); }
   // e.g. a box will have 8 vertices and 12 * 2 indices
   LinesMesh(const std::vector<s_generic_vertex> &vertices_in,
             const std::vector<unsigned int> &indices_in) : vertices(vertices_in), indices(indices_in) {
      init();
   }
   LinesMesh(const clipper::Cell &cell);
   std::vector<s_generic_vertex> vertices;
   std::vector<unsigned int> indices;
   void set_name(const std::string &n) { name = n; }
   void setup(Shader *shader_p);
   void setup_pulse(const glm::vec3 &position, Shader *shader_p);
   void update_buffers_for_pulse(float delta_time); // delta time in ms.
   void draw(Shader *shader_p, const glm::mat4 &mvp, const glm::mat4 &view_rotation, bool use_view_rotation=false);
   void clear();
   bool empty() const { return (vertices.size() == 0); }
};


#endif // LINES_MESH_HH
