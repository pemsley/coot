
#ifndef SHADER_HH
#define SHADER_HH

#include <map>
class Shader {
public:
   enum class Entity_t { NONE = -1, MODEL, MAP, INFRASTRUCTURE, VALIDATION, HUD_TEXT, TEXT_3D,
                         SCREEN, /* i.e. draw the image texture into the screen frame buffer using a quad*/
                         MOLECULAR_TRIANGLES, GENERIC_GRAPHICS_OBJECT};
private:
   enum class ShaderType { NONE = -1, VERTEX = 0, FRAGMENT = 1, GEOMETRY = 2 };
   void parse(const std::string &file_name);
   unsigned int compile_shader(const std::string &source, ShaderType type) const;
   void set_uniform_locations();
   void set_attribute_locations();
   unsigned int create() const;
   std::map<std::string, unsigned int> uniform_location_map;
   unsigned int glGetUniformLocation_internal(const std::string &key);
   std::string VertexSource;
   std::string FragmentSource;
   Entity_t entity_type;
public:
   Shader();
   Shader(const std::string &file_name, Entity_t e);
   Shader(const std::string &vs_file_name,  const std::string &fs_file_name);
   void init(const std::string &file_name, Entity_t e);
   unsigned int program_id;
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
   unsigned int zoom_uniform_location;
   unsigned int map_opacity_uniform_location;
   unsigned int light_0_is_on_uniform_location;
   unsigned int light_1_is_on_uniform_location;
   unsigned int light_0_position_uniform_location;
   unsigned int light_1_position_uniform_location;
   unsigned int light_0_diffuse_colour_uniform_location;
   unsigned int light_1_diffuse_colour_uniform_location;
   // general purpose, not class member uniform locations
   void set_int_for_uniform(const std::string &uniform_name, int value) const;
   void close();
};

#endif // SHADER_HH
