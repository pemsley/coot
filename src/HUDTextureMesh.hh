
#ifndef HUD_TEXTURE_MESH_HH
#define HUD_TEXTURE_MESH_HH

#include <string>
#include <vector>
#include <epoxy/gl.h>
#include <glm/glm.hpp>
#include "g_triangle.hh"
#include "Shader.hh"

class HUDTextureMesh_attribs_t {
public:
   glm::vec2 position;
   glm::vec2 texture_coords;
   HUDTextureMesh_attribs_t(const glm::vec2 &p, const glm::vec2 &tc) : position(p), texture_coords(tc) {}
};

class HUDTextureMesh {
   glm::vec2 position; // uniforms
   glm::vec2 scales;
   GLuint vao;
   GLuint buffer_id;
   GLuint index_buffer_id;
   bool first_time;
   std::vector<HUDTextureMesh_attribs_t> vertices;
   std::vector<g_triangle> triangles;
   std::string name;
   void init();
   void setup_buffers();
   bool draw_this_mesh;

public:
   HUDTextureMesh(const std::string &n) : name(n) { init(); }
   void setup_quad(); // camera-facing, of course
   void set_position_and_scales(const glm::vec2 &pos, const glm::vec2 &scales);
   void setup_texture_coords_for_nbcs_only();
   void setup_texture_coords_for_nbcs_and_rama();
   void draw(Shader *shader_p);
   void close() { draw_this_mesh = false; }
};


#endif // HUD_TEXTURE_MESH_HH
