
#include <iostream>
#include <fstream>

#include <GL/glew.h>

#include "Shader.hh"

Shader::Shader(const std::string &file_name) {

   shader_program_source sps = parse_shader(file_name);

   if (sps.GeometrySource.empty())
      m_programID = CreateShader(sps.VertexSource, sps.FragmentSource);
   else
      m_programID = CreateShader(sps.GeometrySource, sps.VertexSource, sps.FragmentSource);

   glBindAttribLocation(m_programID, 0, "position");
   glBindAttribLocation(m_programID, 1, "texCoord");
   glBindAttribLocation(m_programID, 2, "normal");

   m_uniforms[TRANSFORM_U] = glGetUniformLocation(m_programID, "transform");
}

void
Shader::Bind() {

      glUseProgram(m_programID);
}


shader_program_source
Shader::parse_shader(const std::string &file_name) const {

   enum class ShaderType { NONE = -1, VERTEX = 0, FRAGMENT = 1, GEOMETRY = 2 };

   ShaderType type = ShaderType::NONE;
   shader_program_source ss;
   std::ifstream f(file_name.c_str());
   if (f) {
      std::string line;
      while(std::getline(f, line)) {
	 if (line.find("#shader") != std::string::npos) {
	    if (line.find("vertex") != std::string::npos)
	       type = ShaderType::VERTEX;
	    if (line.find("fragment") != std::string::npos)
	       type = ShaderType::FRAGMENT;
	    if (line.find("geometry") != std::string::npos)
	       type = ShaderType::GEOMETRY;
	 } else {
	    if (type == ShaderType::VERTEX)
	       ss.VertexSource += line + "\n";
	    if (type == ShaderType::FRAGMENT)
	       ss.FragmentSource += line + "\n";
	    if (type == ShaderType::GEOMETRY)
	       ss.GeometrySource += line + "\n";
	 }
      }
   } else {
      std::cout << "Failed to open " << file_name  << std::endl;
   }

   return ss;
}

unsigned int
Shader::compile_shader(const std::string &source, unsigned int type) {

   std::string type_s = "vertex";
   if (type == GL_FRAGMENT_SHADER)
      type_s = "fragment";
   if (type == GL_GEOMETRY_SHADER)
      type_s = "geometry";
   unsigned int id = glCreateShader(type);
   const char *s = source.c_str();
   int l = source.size() + 1;
   glShaderSource(id,  1,  &s, &l);
   glCompileShader(id);

   int result;
   glGetShaderiv(id, GL_COMPILE_STATUS, &result);
   if (result == GL_FALSE) {
      int length;
      glGetShaderiv(id, GL_INFO_LOG_LENGTH, &length);
      char message[length+1];
      glGetShaderInfoLog(id, length, &length, message);
      std::cout << "Failed to compile " << type_s << " shader: " << message << std::endl;
   } else {
      std::cout << "glCompileShader() result was good for " << type_s << " shader " << std::endl;
   } 

   return id;
}

std::string
Shader::file_to_string(const std::string &file_name) const {

   std::ifstream f(file_name.c_str());
   if (f) {
      std::string s((std::istreambuf_iterator<char>(f)),
		    std::istreambuf_iterator<char>());
      return s;
   } else {
      return std::string("");
   }
}

unsigned int
Shader::CreateShader(const std::string &geometry_shader, const std::string &vertex_shader, const std::string &fragment_shader) {

   unsigned int program  = glCreateProgram();
   unsigned int vs = compile_shader(  vertex_shader,   GL_VERTEX_SHADER);
   unsigned int fs = compile_shader(fragment_shader, GL_FRAGMENT_SHADER);
   unsigned int gs = compile_shader(geometry_shader, GL_GEOMETRY_SHADER);

   glAttachShader(program, gs);
   glAttachShader(program, vs);
   glAttachShader(program, fs);
   glLinkProgram(program);
   glValidateProgram(program);

   glDeleteShader(gs);
   glDeleteShader(vs);
   glDeleteShader(fs);

   return program;
}

unsigned int
Shader::CreateShader(const std::string &vertex_shader, const std::string &fragment_shader) {

   unsigned int program  = glCreateProgram();
   std::cout << "CreateShader returned program " << program << std::endl;
   unsigned int vs = compile_shader(vertex_shader, GL_VERTEX_SHADER);
   unsigned int fs = compile_shader(fragment_shader, GL_FRAGMENT_SHADER);

   glAttachShader(program, vs);
   glAttachShader(program, fs);
   glLinkProgram(program);
   glValidateProgram(program);

   glDeleteShader(vs);
   glDeleteShader(fs);

   return program;

}

