
#ifndef HUD_TEXTURE_MESH_HH
#define HUD_TEXTURE_MESH_HH

#include <string>
#include <vector>
#include <epoxy/gl.h>
#include <glm/glm.hpp>
#include "g_triangle.hh"
#include "Shader.hh"

#include "ft-character.hh"
// #include "TextureMesh.hh"

class HUDTextureMesh_attribs_t {
public:
   glm::vec2 position;
   glm::vec2 texture_coords;
   HUDTextureMesh_attribs_t(const glm::vec2 &p, const glm::vec2 &tc) : position(p), texture_coords(tc) {}
};

class HUDTextureMesh {
   enum { VAO_NOT_SET = 99999999 };
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
   HUDTextureMesh() { init(); }
   explicit HUDTextureMesh(const std::string &n) : name(n) { init(); }
   void setup_quad(); // camera-facing, of course
   // for the tooltip background, the position is dynamic (depending on the mouse position)
   // but the scale is fixed - we shouldn't be setting the scale of the tooltop background
   // when the mouse moves.
   void set_position(const glm::vec2 &pos);
   void set_scales(const glm::vec2 &scales);
   void set_position_and_scales(const glm::vec2 &pos, const glm::vec2 &scales);
   void setup_texture_coords_for_nbcs_only();
   void setup_texture_coords_for_nbcs_and_rama();
   void draw(Shader *shader_p);
   void draw_label(const std::string &label, bool highlight_label_flag, Shader *shader_p,
                   const std::map<GLchar, FT_character> &ft_characters);
   void draw_label(const std::string &label, glm::vec4 &text_colour, Shader *shader_p,
                   const std::map<GLchar, FT_character> &ft_characters);
   void close() { draw_this_mesh = false; }
};


#endif // HUD_TEXTURE_MESH_HH
