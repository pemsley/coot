
#ifndef HUD_MESH_HH
#define HUD_MESH_HH

#include <string>
#include <vector>
#include <epoxy/gl.h>
#include <glm/glm.hpp>
#include "g_triangle.hh"
#include "Shader.hh"

class HUD_bar_attribs_t {
public:
   glm::vec4 colour;
   glm::vec2 position_offset;
   // if this_is_a_non_moving_atoms_residue then we want to make 2 bars, an 80% thickness bar
   // and a small "grey" bar below. same x coords.
   float scale_x;
   float scale_y;
   HUD_bar_attribs_t(const glm::vec4 &c, const glm::vec2 &p, const float &l) :
      colour(c), position_offset(p), scale_x(l) {
      scale_y = 1.0;
   }
};

class HUDMesh {
   void setup_buffers();
   void init();
   bool first_time;
   unsigned int max_n_instances;
   unsigned int n_instances;
   unsigned int inst_hud_bar_attribs_buffer_id;
public:
   GLuint vao;
   GLuint buffer_id;
   GLuint index_buffer_id;
   bool this_mesh_is_closed;
   bool use_blending;
   std::vector<glm::vec2> vertices;
   std::vector<g_triangle> triangles;
   std::string name;
   HUDMesh() { init(); }
   HUDMesh(const std::string &n) : name(n) { init(); }
   void setup_camera_facing_quad_for_bar();
   void setup_instancing_buffer(unsigned int n_boxes);
   void update_instancing_buffer_data(const std::vector<HUD_bar_attribs_t> &new_bars);
   void draw(Shader *shader);
   void close() { this_mesh_is_closed = true; }

};

#endif // HUD_MESH_HH
