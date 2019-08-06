
#ifndef SHADER_HH
#define SHADER_HH

#include <map>
class Shader {
public:
   enum class Entity_t { NONE = -1, MODEL, MAP, INFRASTRUCTURE, VALIDATION, HUD_TEXT, TEXT_3D};
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
   void init(const std::string &file_name, Entity_t e);
   unsigned int program_id;
   unsigned int get_program_id() const { return program_id; }
   unsigned int view_rotation_uniform_location;
   unsigned int mvp_uniform_location;
   unsigned int background_colour_uniform_location;
   unsigned int eye_position_uniform_location;
};

#endif // SHADER_HH
