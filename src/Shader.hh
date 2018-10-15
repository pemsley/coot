

#ifndef SHADER_HH
#define SHADER_HH

#include <string>

struct shader_program_source {
   std::string VertexSource;
   std::string FragmentSource;
   std::string GeometrySource;
};

class Shader {
   shader_program_source parse_shader(const std::string &file_name) const;
   unsigned int compile_shader(const std::string &source, unsigned int type);
   std::string file_to_string(const std::string &file_name) const;
   unsigned int CreateShader(const std::string &vertex_shader,
			     const std::string &fragment_shader);
   unsigned int CreateShader(const std::string &geometry_shader,
			     const std::string &vertex_shader,
			     const std::string &fragment_shader);
   GLuint m_programID;

   enum { TRANSFORM_U, NUM_UNIFORMS };
   GLuint m_uniforms[NUM_UNIFORMS];

public:
   Shader() {}
   Shader(const std::string &file_name);
   void Bind();
   GLuint get_programID() const { return m_programID; }

};

#endif // SHADER_HH
