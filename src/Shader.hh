
#ifndef SHADER_HH
#define SHADER_HH

#include "lights-info.hh"

#include <map>
class Shader {
public:
   enum class Entity_t { NONE = -1, MODEL, MAP, INFRASTRUCTURE, VALIDATION, HUD_TEXT, TEXT_3D,
                         SCREEN, /* i.e. draw the image texture into the screen frame buffer using a quad */
                         MOLECULAR_TRIANGLES,
                         GENERIC_DISPLAY_OBJECT,
                         INSTANCED_DISPLAY_OBJECT};
private:
   enum class ShaderType { NONE = -1, VERTEX = 0, FRAGMENT = 1, GEOMETRY = 2 };
   std::string default_directory; // if not in current directory try to find the shader here
   void parse(const std::string &file_name);
   unsigned int compile_shader(const std::string &source, ShaderType type) const;
   void set_uniform_locations();
   void set_attribute_locations();
   unsigned int create() const;
   std::map<std::string, unsigned int> uniform_location_map;
   GLint glGetUniformLocation_internal(const std::string &key);
   std::string VertexSource;
   std::string FragmentSource;
   Entity_t entity_type;
public:
   Shader();
   Shader(const std::string &file_name, Entity_t e);
   Shader(const std::string &vs_file_name,  const std::string &fs_file_name);
   void init(const std::string &file_name, Entity_t e);
   void set_default_directory(const std::string &dir);
   unsigned int program_id;
   std::string name; // for warning/error messages
   void Use();
   unsigned int get_program_id() const { return program_id; } // or use above function
   // consider a map of these - the variable name is the same as the uniform name
   // in the shader.
   virtual void set_more_uniforms_for_molecular_triangles();
   unsigned int view_rotation_uniform_location;
   unsigned int mvp_uniform_location;
   unsigned int background_colour_uniform_location;
   unsigned int line_colour_uniform_location;
   unsigned int eye_position_uniform_location;
   unsigned int hud_projection_uniform_location;
   unsigned int atom_label_projection_uniform_location;
   unsigned int atom_label_textColour_uniform_location;

   // general purpose, not class member uniform locations
   void set_int_for_uniform(const std::string &uniform_name, int value);
   void set_bool_for_uniform(const std::string &uniform_name, bool value);
   void set_float_for_uniform(const std::string &uniform_name, float v);
   void set_vec4_for_uniform(const std::string &uniform_name, float f0, float f1, float f2, float f3);
   void set_vec4_for_uniform(const std::string &uniform_name, const glm::vec4 &v);
   void set_vec3_for_uniform(const std::string &uniform_name, const glm::vec3 &v);
   void setup_light(unsigned int light_index,
                    const lights_info_t &light,
                    const glm::mat4 &view_rotation_matrix); // mouse trackball rotation
   void close();
};

#endif // SHADER_HH
