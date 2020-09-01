
#ifndef HUD_TEXTURE_MESH_HH
#define HUD_TEXTURE_MESH_HH

#include "g_triangle.hh"

#include "lights-info.hh"
#include "Shader.hh"
#include "Material.hh"

class HUDTextureMesh {
   GLuint vao;
   GLuint buffer_id;
   GLuint index_buffer_id;
   std::vector<glm::vec2> vertices;
   std::vector<g_triangle> triangles;
   std::string name;
   void init();
   void setup_buffers();
   bool draw_this_mesh;
   Material material;

public:
   HUDTextureMesh(const std::string &n) : name(n) { init(); }
   void setup_quad();
   void set_material;
   void draw(Shader *shader_p);
   void close() { draw_this_mesh = false; }
};


#endif // HUD_TEXTURE_MESH_HH
