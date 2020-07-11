
#include "generic-vertex.hh"
#include "gl-lights-info.hh"
#include "Shader.hh"
#include "Material.hh"
#include "g_triangle.hh"

class graphical_molecule {
   void setup_instancing_buffers(); // or buffers, when we add rotation
   Material material;
   void setup_instanced_balls(Shader *shader_p, const Material &material_in);
   void setup_instanced_dodecs(Shader *shader_p, const Material &material_in);

public:
   graphical_molecule() { is_instanced = false; normals_setup = false; }
   graphical_molecule(const std::vector<s_generic_vertex> &v,
                      const std::vector<g_triangle> &t)
      : vertices(v), triangle_vertex_indices(t) {
      is_instanced = false;
      normals_setup = false;
      setup_buffers();
   }
   GLuint vao;
   GLuint buffer_id;
   GLuint inst_model_translation_buffer_id;
   GLuint inst_colour_buffer_id;
   GLuint index_buffer_id;
   bool is_instanced;
   std::vector<s_generic_vertex> vertices;
   std::vector<g_triangle> triangle_vertex_indices;
   void setup_simple_triangles(Shader *shaderp, const Material &material_in);
   // setup_simple_triangles() calls fill_with_simple_triangles_vertices()
   void fill_with_simple_triangles_vertices();
   void fill_with_direction_triangles();
   void setup_instanced_debugging_objects(Shader *shader_p, const Material &material_in);

   bool normals_setup;
   GLuint vao_normals;
   GLuint normals_buffer_id;

   // can be considered as "draw_self()"
   void draw(Shader *shader,
             const glm::mat4 &mvp,
             const glm::mat4 &world_rotation_matrix,
             const glm::mat4 &world_rotation_translation_matrix,
             const std::map<unsigned int, gl_lights_info_t> &lights,
             const glm::vec3 &eye_position);
   void draw_normals(const glm::mat4 &mvp); // debugging
   Shader shader_for_draw_normals;

   void setup_rama_balls(Shader *shader_p, const Material &material_in); // call fill_rama_balls() and setup_buffers()
   void fill_rama_balls(); // make up some balls
   void add_one_ball(float scale, const glm::vec3 &centre);
   void add_one_origin_ball();
   void add_one_origin_dodec();
   void fill_one_dodec();
   void setup_buffers();
};


