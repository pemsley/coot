
#ifndef INSTANCED_MARKUP_MESH_HH
#define INSTANCED_MARKUP_MESH_HH

#include "generic-vertex.hh"
#include "g_triangle.hh"

#include "lights-info.hh"
#include "Shader.hh"
#include "Material.hh"
#include "Particle.hh"

// to be used (in the first case for rama balls) and later hopefully
// for contact dots

// the vertices are glm::vec3
// and the instanced variables are colour size position

class Instanced_Markup_Mesh_Vertex_attrib_t {
public:
   glm::vec3 position;
   float displacement; // so the balls can be "rough"
};

// Are you using the rama-balls.shader?

class Instanced_Markup_Mesh_attrib_t {
public:
   glm::vec4 colour;
   glm::vec3 position;
   float size;
   float specular_strength;
   float shininess;
   Instanced_Markup_Mesh_attrib_t() {}
   Instanced_Markup_Mesh_attrib_t(const glm::vec4 &c, glm::vec3 &p, float s) : colour(c), position(p), size(s) {
      specular_strength = 0.5;
      shininess = 64.0;
   }
};



class Instanced_Markup_Mesh {
   GLuint vao;
   GLuint buffer_id;
   GLuint index_buffer_id;
   GLuint inst_attribs_buffer_id;
   unsigned int n_instances;
   unsigned int max_n_instances;
   std::vector<Instanced_Markup_Mesh_Vertex_attrib_t> vertices;
   std::vector<g_triangle> triangles;
   std::string name;
   void init();
   void setup_buffers();
   bool draw_this_mesh;
   bool first_time;
   bool this_mesh_is_closed;

public:
   Instanced_Markup_Mesh(const std::string &n) : name(n) { init(); }
   bool is_closed() const {return this_mesh_is_closed; }  // once closed, it's gone.
   std::string get_name() const { return name; }
   void setup_octasphere(unsigned int num_subdivisions);
   void setup_instancing_buffers(unsigned int max_nun_instances); // allocate space, don't fill it with data yet
   void update_instancing_buffers(const std::vector<Instanced_Markup_Mesh_attrib_t> &balls);
   bool get_draw_status() const { return draw_this_mesh; }
   void set_draw_status(bool status) { draw_this_mesh = status; }
   void draw(Shader *shader_p,
             const glm::mat4 &mvp,
             const glm::mat4 &view_rotation_matrix,
             const std::map<unsigned int, lights_info_t> &lights,
             const glm::vec3 &eye_position, // eye position in view space (not molecule space)
             const glm::vec4 &background_colour,
             bool do_depth_fog);
   void close();
   void clear();
};


#endif // INSTANCED_MARKUP_MESH_HH
