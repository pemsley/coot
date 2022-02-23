
#ifndef TEXTURE_MESH_HH
#define TEXTURE_MESH_HH

#include <vector>
#include <string>
#include <epoxy/gl.h>
#include <glm/glm.hpp>
#include "g_triangle.hh"

#include "obj_loader.h"

#include "Shader.hh"
#include "Texture.hh" // now TextureMesh contains a vector of Textures - I am not sure this is a good
                      // arrangement (it is a result of load_from_glTF() which also fill/creates textures)

class TextureMeshVertex {
public:
   glm::vec3 position;
   glm::vec3 normal;
   glm::vec3 tangent;   // calculated later, using triangles
   glm::vec3 bitangent; // (same)
   glm::vec4 color;
   glm::vec2 texCoord;
   TextureMeshVertex(const glm::vec3 &p, const glm::vec3 &n, const glm::vec4 &col, const glm::vec2 &tc) :
      position(p), normal(n), color(col), texCoord(tc) {}
};

class TextureInfoType {
public:
   Texture texture;
   std::string name;
   // std::string sampler_name; // e.g. "base_texture", sampler is now in Texture::type
   GLuint unit;
   TextureInfoType(const Texture &t, const std::string &n) :
      texture(t), name(n) {}
};

class TextureMesh {
   enum { VAO_NOT_SET = 99999999 };
   GLuint vao;
   GLuint buffer_id;
   GLuint index_buffer_id;
   std::vector<TextureMeshVertex> vertices;
   std::vector<g_triangle> triangles;
   std::string name;
   std::string file_name;

   int n_instances; // number of instances to be _drawn_
   int n_instances_allocated; // that we made space for in glBufferData()
   bool is_instanced;
   // note: the instanced data for happy-face-residue-markers is just the position
   // (at least until this gets combined with sad-face-residue-markers/textures).
   // A combined (next to each other?) texture? But do we need that microoptimization?
   // The opacity will be a uniform.
   unsigned int draw_count; // so that I can animate the happy faces depending on the draw_count
   unsigned int inst_positions_id;

public:
   TextureMesh() : vao(VAO_NOT_SET), index_buffer_id(VAO_NOT_SET), draw_this_mesh(true) {
      n_instances_allocated = 0;
      n_instances = 0;
      is_instanced = false;
      inst_positions_id = -1;
      draw_count = 0;
   }
   explicit TextureMesh(const std::string &n):
      vao(VAO_NOT_SET), index_buffer_id(VAO_NOT_SET), name(n), draw_this_mesh(true) {
      n_instances_allocated = 0;
      n_instances = 0;
      is_instanced = false;
      inst_positions_id = -1;
      draw_count = 0;
   }
   bool draw_this_mesh;
   std::vector<TextureInfoType> textures;
   void import(const IndexedModel &ind_model, float scale);
   void import(const std::vector<TextureMeshVertex> &vertices, const std::vector<g_triangle> &triangles_in);
   bool have_instances() const { return is_instanced; }
   void setup_tbn(unsigned int n_vertices); // tangent bitangent normal, pass the n_vertices for validation of indices.
   void setup_camera_facing_quad(Shader *shader_p, float scale_x, float scale_y);
   void setup_buffers();
   void set_colour(const glm::vec4 &col_in);
   void setup_instancing_buffers(unsigned int n_happy_faces_max); // setup the buffer, don't add data
   std::string get_name() const { return name; }
   // this is for an ephemeral instanced texturemesh
   void update_instancing_buffer_data_for_happy_faces(const std::vector<glm::vec3> &positions, // in 3D space (of the CAs)
                                                      unsigned int draw_count_in,
                                                      unsigned int draw_count_max,
                                                      const glm::vec3 &screen_y_uv);
   void update_instancing_buffer_data(const std::vector<glm::vec3> &positions);
   void add_texture(const TextureInfoType &ti) { textures.push_back(ti); }
   void translate(const glm::vec3 &t); // include call to setup_buffers();
   void apply_scale(const float &sf); // include call to setup_buffers();
   void apply_transformation(const glm::mat4 &m);  // transform the positions in the vertices
   void draw(Shader *shader,
             const glm::mat4 &mvp,
             const glm::mat4 &view_rotation_matrix,
             const std::map<unsigned int, lights_info_t> &lights,
             const glm::vec3 &eye_position, // eye position in view space (not molecule space)
             const glm::vec4 &background_colour,
             bool do_depth_fog);
   void draw_with_shadows(Shader *shader,
                          const glm::mat4 &mvp,
                          const glm::mat4 &view_rotation_matrix,
                          const std::map<unsigned int, lights_info_t> &lights,
                          const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                          const glm::vec4 &background_colour,
                          bool do_depth_fog,
                          const glm::mat4 &light_space_mvp,
                          unsigned int depthMap,
                          float shadow_strength,
                          unsigned int shadow_softness,
                          bool show_just_shadows);
   void draw_for_ssao(Shader *shader_for_tmeshes_p,
                      const glm::mat4 &model,
                      const glm::mat4 &view,
                      const glm::mat4 &projection);
   void draw_atom_label(const std::string &atom_label,
                        const glm::vec3 &atom_label_position,
                        const glm::vec4 &text_colour, // set using subbufferdata
                        Shader *shader,
                        const glm::mat4 &mvp,
                        const glm::mat4 &view_rotation_matrix,
                        const glm::vec4 &background_colour,
                        bool do_depth_fog,
                        bool is_perspective_projection);
   // draw an ephemeral instanced opacity-varying texturemesh.
   // Other draw_instances() functions may be needed in future, if so change the name of this one.
   void draw_instances(Shader *shader_p, const glm::mat4 &mvp, const glm::mat4 &view_rotation,
                       unsigned int draw_count, unsigned int draw_count_max);

   bool load_from_glTF(const std::string &file_name, bool include_call_to_setup_buffers=true);
};

#endif // TEXTURE_MESH_HH
