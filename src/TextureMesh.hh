
#ifndef TEXTURE_MESH_HH
#define TEXTURE_MESH_HH

#include <vector>
#include <string>
#include <epoxy/gl.h>
#include <glm/glm.hpp>
#include "g_triangle.hh"

#include "obj_loader.h"

#include "Shader.hh"

class TextureMeshVertex {
public:
   glm::vec3 position;
   glm::vec3 normal;
   glm::vec4 color;
   glm::vec2 texCoord;
   TextureMeshVertex(const glm::vec3 &p, const glm::vec3 &n, const glm::vec4 &col, const glm::vec2 &tc) :
      position(p), normal(n), color(col), texCoord(tc) { }
};

class TextureMesh {
   GLuint vao;
   GLuint buffer_id;
   GLuint index_buffer_id;
   std::vector<TextureMeshVertex> vertices;
   std::vector<g_triangle> triangles;
   std::string name;

public:
   TextureMesh() { draw_this_mesh = true; vao = 99999999; index_buffer_id = 99999999; }
   TextureMesh(const std::string &n) { name = n;
      draw_this_mesh = true; vao = 99999999; index_buffer_id = 99999999; }
   bool draw_this_mesh;
   void import(const IndexedModel &ind_model, float scale);
   void setup_camera_facing_quad(Shader *shader_p);
   void setup_buffers();
   void set_colour(const glm::vec4 &col_in);
   void draw(Shader *shader,
             const glm::mat4 &mvp,
             const glm::mat4 &view_rotation_matrix,
             const std::map<unsigned int, lights_info_t> &lights,
             const glm::vec3 &eye_position, // eye position in view space (not molecule space)
             const glm::vec4 &background_colour,
             bool do_depth_fog);
   void draw_atom_label(const std::string &atom_label,
                        const glm::vec3 &atom_label_position,
                        const glm::vec4 &text_colour, // set using subbufferdata
                        Shader *shader,
                        const glm::mat4 &mvp,
                        const glm::mat4 &view_rotation_matrix,
                        const std::map<unsigned int, lights_info_t> &lights,
                        const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                        const glm::vec4 &background_colour,
                        bool do_depth_fog);

};

#endif // TEXTURE_MESH_HH
