
#ifndef HUD_TEXTURE_MESH_HH
#define HUD_TEXTURE_MESH_HH

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

public:
   HUDTextureMesh(const std::string &n) : name(n) { init(); }
   void setup_quad();
   void draw(Shader *shader_p);
   void close() { draw_this_mesh = false; }
};


#endif // HUD_TEXTURE_MESH_HH
