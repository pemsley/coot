
#include <string>
#include <iostream>
#include <fstream>
#include <epoxy/gl.h>

#include "Shader.hh"

Shader::Shader() {
   program_id = 0; // unset
}

Shader::Shader(const std::string &file_name) {
   init(file_name);
}

void Shader::init(const std::string &file_name) {
   parse(file_name);
   if (! VertexSource.empty()) {
      if (! FragmentSource.empty()) {
         program_id = create();
         std::cout << "debug() Shader::init() program_id " << program_id << std::endl;
         if (true) {
            set_uniform_locations();
            glBindAttribLocation(program_id, 0, "position"); // use set_attribute_locations()
            glBindAttribLocation(program_id, 1, "normal");   // for consistency
            glBindAttribLocation(program_id, 2, "colour");
            glBindAttribLocation(program_id, 3, "translate_position");
            glBindAttribLocation(program_id, 4, "model_matrix");
         }
      } else {
         std::cout << "Empty Fragment Shader source\n";
      }
   } else {
      std::cout << "Empty Vertex Shader source\n";
   }
}

void Shader::set_uniform_locations() {
   GLuint err;
   mvp_uniform_location           = glGetUniformLocation(program_id, "mvp");
   err = glGetError(); if (err) std::cout << "error:: set_uniform_locations() error 1: " << err << std::endl;
   view_rotation_uniform_location = glGetUniformLocation(program_id, "view_rotation");
   err = glGetError(); if (err) std::cout << "error:: set_uniform_locations() error 2: " << err << std::endl;
   std::cout << "debug:: set_uniform_locations() " << mvp_uniform_location << " " << view_rotation_uniform_location
             << std::endl;
}

void Shader::parse(const std::string &file_name) {
   std::ifstream f(file_name.c_str());
   if (f) {
      std::string line;
      ShaderType type(ShaderType::NONE);
      while(std::getline(f, line)) {
         if (line.find("#shader") != std::string::npos) {
            if (line.find("vertex") != std::string::npos)
            type = ShaderType::VERTEX;
            if (line.find("fragment") != std::string::npos)
            type = ShaderType::FRAGMENT;
         } else {
            if (type == ShaderType::VERTEX)
            VertexSource += line + "\n";
            if (type == ShaderType::FRAGMENT)
            FragmentSource += line + "\n";
         }
      }
   } else {
      std::cout << "Failed to open " << file_name  << std::endl;
   }
}

unsigned int Shader::compile_shader(const std::string &source, ShaderType type) const {

   unsigned int i_type = 0;
   std::string type_s = "-";
   if (type == ShaderType::FRAGMENT) {
      i_type = GL_FRAGMENT_SHADER;
      type_s = "Fragment";
   }
   if (type == ShaderType::VERTEX) {
      i_type = GL_VERTEX_SHADER;
      type_s = "Vertex";
   }
   unsigned int id = glCreateShader(i_type);
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
      std::cout << "error:: Failed to compile " << type_s << " shader: " << message << std::endl;
   } else {
      std::cout << "glCompileShader() result was good for " << type_s << " shader " << std::endl;
   }
   return id;
}

unsigned int Shader::create() const {
   unsigned int program = glCreateProgram();
   unsigned int vs = compile_shader(  VertexSource, ShaderType::VERTEX);
   unsigned int fs = compile_shader(FragmentSource, ShaderType::FRAGMENT);

   glAttachShader(program, vs);
   glAttachShader(program, fs);
   glLinkProgram(program);
   glValidateProgram(program);

   glDeleteShader(vs);
   glDeleteShader(fs);

   return program;
}
