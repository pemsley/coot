
#include <string>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <epoxy/gl.h>

#include "Shader.hh"

Shader::Shader() {
   program_id = 0; // unset
   zoom_uniform_location = -999; // for debugging
   map_opacity_uniform_location = -999;
}

Shader::Shader(const std::string &file_name, Shader::Entity_t e) {
   init(file_name, e);
}

Shader::Shader(const std::string &vs_file_name, const std::string &fs_file_name) {

   entity_type = Entity_t::HUD_TEXT; // hackety-hack
   program_id = glCreateProgram();

   parse(vs_file_name);
   if (! VertexSource.empty()) {
      unsigned int vs = compile_shader(VertexSource, ShaderType::VERTEX);
      parse(fs_file_name);
      if (! FragmentSource.empty()) {
         unsigned int fs = compile_shader(FragmentSource, ShaderType::FRAGMENT);

         glAttachShader(program_id, vs);
         glAttachShader(program_id, fs);
         glLinkProgram(program_id);
         glValidateProgram(program_id);
      } else {
        std::cout << "Oops - empty Fragment shader" << std::endl;
      }
   }
}

void Shader::init(const std::string &file_name_in, Shader::Entity_t e) {

   // clear then go
   VertexSource.clear();
   FragmentSource.clear();

   // is the shader in the current directory? If not, try to find
   // it in default_directory

   // I don't want this file to depend on utils, so include this in-line
   std::string file_name(file_name_in);
   bool file_exists = false;
   struct stat s;
   int fstat = stat(file_name.c_str(), &s);
   if ( fstat != -1 )
      file_exists = true;
   if (! file_exists)
      file_name = default_directory + "/" + file_name;
   fstat = stat(file_name.c_str(), &s);
   if ( fstat != -1 )
      file_exists = true;

   std::cout << "::: Shader compile " << file_name << std::endl;
   if (! file_exists) {
      std::cout << "WARNING:: Missing file " << file_name << std::endl;
      return;
   }

   entity_type = e;
   parse(file_name);
   if (! VertexSource.empty()) {
      if (! FragmentSource.empty()) {
         program_id = create();
         std::cout << "debug:: Shader::init() " << file_name << " program_id " << program_id << std::endl;
         if (true) {
            Use();
            set_uniform_locations();
            set_attribute_locations();
         }
      } else {
         std::cout << "Empty Fragment Shader source\n";
      }
   } else {
      std::cout << "Empty Vertex Shader source\n";
   }
}

void
Shader::set_default_directory(const std::string &s) {

   if (!s.empty())
      default_directory = s;
}

void
Shader::close() {
   glDeleteProgram (program_id);
}

void
Shader::set_int_for_uniform(const std::string &uniform_name, int value) const {
   GLuint err = glGetError();
   if (err) std::cout << "set_int_for_uniform() start err " << err << std::endl;
   GLuint loc = glGetUniformLocation(program_id, uniform_name.c_str());
   err = glGetError(); if (err) std::cout << "set_int_for_uniform() A err " << err << std::endl;
   glUniform1i(loc,value);
   err = glGetError(); if (err) std::cout << "set_int_for_uniform() B err " << err << std::endl;
}


void
Shader::Use() {
   GLuint err = glGetError();
   if (err) std::cout << "Shader::Use() pre glUseProgram() err " << err << std::endl;
   glUseProgram(program_id);
   err = glGetError();
   if (err) std::cout << "Shader::Use() err " << err << " for program_id " << program_id << std::endl;
}

void
Shader::set_attribute_locations() {

   if (entity_type == Entity_t::MODEL) {
      glBindAttribLocation(program_id, 0, "model_rotation_matrix_0");
      glBindAttribLocation(program_id, 1, "model_rotation_matrix_1");
      glBindAttribLocation(program_id, 2, "model_rotation_matrix_2");
      glBindAttribLocation(program_id, 3, "model_translation");
      glBindAttribLocation(program_id, 4, "position");
      glBindAttribLocation(program_id, 5, "normal");
      glBindAttribLocation(program_id, 6, "colour");
   }
   if (entity_type == Entity_t::MAP) {
      glBindAttribLocation(program_id, 0, "position");
      glBindAttribLocation(program_id, 1, "normal");
      glBindAttribLocation(program_id, 2, "colour");
   }
   if (entity_type == Entity_t::HUD_TEXT) {
      glBindAttribLocation(program_id, 0, "vertex"); // 2 x 2 pos, texture
   }
   if (entity_type == Entity_t::TEXT_3D) { // atom label
      glBindAttribLocation(program_id, 0, "vertex"); // 2 x 2 pos, texture
   }
}

unsigned int
Shader::glGetUniformLocation_internal(const std::string &key) {

   // don't ask the hardware about the location of the uniform if we
   // have asked before.

   // sort of strange thing - because I store the locations - Hmm.

   std::map<std::string, GLuint>::const_iterator it = uniform_location_map.find(key);
   if (it != uniform_location_map.end()) {
      return it->second;
   } else {
      GLuint l = glGetUniformLocation(program_id, key.c_str());
      uniform_location_map[key] = l;
      return l;
   }
}

void Shader::set_uniform_locations() {
   GLuint err;

   if (entity_type == Entity_t::MODEL ||
       entity_type == Entity_t::MAP ||
       entity_type == Entity_t::MOLECULAR_TRIANGLES ||
       entity_type == Entity_t::GENERIC_GRAPHICS_OBJECT) {
      err = glGetError(); if (err) std::cout << "error:: set_uniform_locations() error 0: " << err << std::endl;
      mvp_uniform_location           = glGetUniformLocation_internal("mvp");
      err = glGetError(); if (err) std::cout << "error:: set_uniform_locations() error 1: " << err << std::endl;
      view_rotation_uniform_location = glGetUniformLocation_internal("view_rotation");
      err = glGetError(); if (err) std::cout << "error:: set_uniform_locations() error 2: " << err << std::endl;
      background_colour_uniform_location = glGetUniformLocation_internal("background_colour");
      err = glGetError(); if (err) std::cout << "error:: set_uniform_locations() error 3: " << err << std::endl;
      eye_position_uniform_location = glGetUniformLocation_internal("eye_position");
      err = glGetError(); if (err) std::cout << "error:: set_uniform_locations() error 4: " << err << std::endl;

      is_perspective_projection_uniform_location = glGetUniformLocation_internal("is_perspective_projection");
      light_0_is_on_uniform_location = glGetUniformLocation_internal("light_0_is_on");
      light_1_is_on_uniform_location = glGetUniformLocation_internal("light_1_is_on");
      light_0_position_uniform_location = glGetUniformLocation_internal("light_0_position");
      light_1_position_uniform_location = glGetUniformLocation_internal("light_1_position");
      light_0_diffuse_colour_uniform_location = glGetUniformLocation_internal("light_0_diffuse_colour");
      light_1_diffuse_colour_uniform_location = glGetUniformLocation_internal("light_1_diffuse_colour");
      err = glGetError(); if (err) std::cout << "error:: set_uniform_locations() error 5: " << err << std::endl;

      // the compiler can "throw these away" 4294967295 if they are not used in the fragment shader (it optimizes)
      if (false)
         std::cout << "debug:: set_uniform_locations() mvp: "
                   << mvp_uniform_location << " view-rot: "
                   << view_rotation_uniform_location << " bg: "
                   << background_colour_uniform_location << " eye_pos: "
                   << eye_position_uniform_location << std::endl;
   }
   if (entity_type == Entity_t::MOLECULAR_TRIANGLES) {
      set_more_uniforms_for_molecular_triangles();
   }
   if (entity_type == Entity_t::MAP) {
      map_opacity_uniform_location = glGetUniformLocation_internal("map_opacity");
      err = glGetError(); if (err) std::cout << "error:: set_uniform_locations() error 1b: " << err << std::endl;
   }
   if (entity_type == Entity_t::INFRASTRUCTURE) {
      mvp_uniform_location           = glGetUniformLocation_internal("mvp");
      err = glGetError(); if (err) std::cout << "error:: set_uniform_locations() error 1c: " << err << std::endl;
      view_rotation_uniform_location = glGetUniformLocation_internal("view_rotation");
      err = glGetError(); if (err) std::cout << "error:: set_uniform_locations() error 2c: " << err << std::endl;
      line_colour_uniform_location = glGetUniformLocation_internal("line_colour");
      err = glGetError(); if (err) std::cout << "error:: set_uniform_locations() error 3c: " << err << std::endl;
      background_colour_uniform_location = glGetUniformLocation_internal("background_colour");
      err = glGetError(); if (err) std::cout << "error:: set_uniform_locations() error 4c: " << err << std::endl;
   }
   if (entity_type == Entity_t::HUD_TEXT) {
      hud_projection_uniform_location           = glGetUniformLocation_internal("projection");
      err = glGetError(); if (err) std::cout << "error:: set_uniform_locations() error 5d: " << err << std::endl;
   }
   if (entity_type == Entity_t::TEXT_3D) {
      atom_label_projection_uniform_location           = glGetUniformLocation_internal("projection"); // may change
      err = glGetError(); if (err) std::cout << "error:: set_uniform_locations() error 6a: " << err << std::endl;
      atom_label_textColour_uniform_location           = glGetUniformLocation_internal("textColour");
      err = glGetError(); if (err) std::cout << "error:: set_uniform_locations() error 6b: " << err << std::endl;
   }
   if (entity_type == Entity_t::SCREEN) {
      zoom_uniform_location = glGetUniformLocation_internal("zoom");
      err = glGetError(); if (err) std::cout << "error:: set_uniform_locations() error 7a: " << err << std::endl;
      is_perspective_projection_uniform_location = glGetUniformLocation_internal("is_perspective_projection");
      err = glGetError(); if (err) std::cout << "error:: set_uniform_locations() error 7b: " << err << std::endl;
   }
}

void Shader::set_more_uniforms_for_molecular_triangles() {

   // put more uniforms here
}

void Shader::parse(const std::string &file_name) {
   std::ifstream f(file_name.c_str());
   if (f) {
      VertexSource.clear();
      FragmentSource.clear();
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
      std::cout << "WARNING:: Shade::parse(): Failed to open " << file_name  << std::endl;
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
      std::cout << "   glCompileShader() result was good for " << type_s << " shader " << std::endl;
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
   GLuint err = glGetError();
   if (err)
      std::cout << "Shader::create() err " << err << std::endl;
   else
      std::cout << "   Shader::create() link was good " << std::endl;

   glDeleteShader(vs);
   glDeleteShader(fs);

   return program;
}
