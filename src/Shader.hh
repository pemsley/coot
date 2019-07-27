

class Shader {
      enum class ShaderType { NONE = -1, VERTEX = 0, FRAGMENT = 1, GEOMETRY = 2 };
      void parse(const std::string &file_name);
      unsigned int compile_shader(const std::string &source, ShaderType type) const;
      void set_uniform_locations();
      void set_attribute_locations();
      unsigned int create() const;
      std::string VertexSource;
      std::string FragmentSource;
public:
   Shader();
   Shader(const std::string &file_name);
   void init(const std::string &file_name);
   unsigned int program_id;
   unsigned int get_program_id() const { return program_id; }
   unsigned int view_rotation_uniform_location;
   unsigned int mvp_uniform_location;
};
