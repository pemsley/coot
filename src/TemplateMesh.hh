
#ifndef HUD_TEXTURE_MESH_HH
#define HUD_TEXTURE_MESH_HH

#include "g_triangle.hh"

#include "lights-info.hh"
#include "Shader.hh"
#include "Material.hh"

class XYMesh {
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
   bool this_mesh_is_closed;

public:
   XYMesh(const std::string &n) : name(n) { init(); }
   std::string get_name() const { return name; }
   bool get_draw_status() const { return draw_this_mesh; }
   void setup_quad();
   void set_material(Material mat) { material = mat; }
   void draw(Shader *shader_p);
   void close() {
      draw_this_mesh = false;
      this_mesh_is_closed = true; // and delete the buffers if not first time,
                                  //  so don't inline this function
   }
};


#endif // HUD_TEXTURE_MESH_HH
