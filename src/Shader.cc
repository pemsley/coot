
#include <string>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <epoxy/gl.h>
#include <glm/gtc/type_ptr.hpp>
#define GLM_ENABLE_EXPERIMENTAL 
#include <glm/gtx/string_cast.hpp>

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

void Shader::init(const std::string &file_name, Shader::Entity_t e) {

   // clear then go
   VertexSource.clear();
   FragmentSource.clear();
   std::string::size_type pos = file_name.find_first_of(".shader");
   name = file_name;
   std::cout << "::: Shader compile " << file_name << std::endl;
   bool file_exists = true; // fixme
   if (! file_exists) {
      std::cout << "WARNING:: Missing file " << file_name << std::endl;
      return;
   }

   entity_type = e;
   parse(file_name);
   if (! VertexSource.empty()) {
      if (! FragmentSource.empty()) {
         program_id = create();
         Use();
         set_uniform_locations();
         set_attribute_locations();
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
Shader::set_int_for_uniform(const std::string &uniform_name, int value) {
   GLuint err = glGetError();
   if (err) std::cout << "set_int_for_uniform() start err " << err << std::endl;
   GLint loc = glGetUniformLocation_internal(uniform_name.c_str());
   err = glGetError(); if (err) std::cout << "set_int_for_uniform() A err " << err << std::endl;
   glUniform1i(loc,value);
   err = glGetError(); if (err) std::cout << "set_int_for_uniform() B err " << err << std::endl;
}

void
Shader::set_bool_for_uniform(const std::string &uniform_name, bool value) {

   GLuint err = glGetError();
   if (err) std::cout << "set_bool_for_uniform() " << uniform_name << " start err " << err << std::endl;
   GLint loc = glGetUniformLocation_internal(uniform_name.c_str());
   // std::cout << "set_bool_for_uniform() got loc " << loc << std::endl;
   err = glGetError();
   if (err) std::cout << "ERROR:: " << name << " set_bool_for_uniform() " << uniform_name << " A err "
                      << err << std::endl;
   glUniform1i(loc, value);
   err = glGetError();
   if (err) std::cout << "ERROR:: " << name << " set_bool_for_uniform() " << uniform_name << " B err "
                      << err << std::endl;
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

GLint
Shader::glGetUniformLocation_internal(const std::string &key) {

   // don't ask the hardware about the location of the uniform if we
   // have asked before.

   // sort of strange thing - because I store the locations - Hmm.

   std::map<std::string, GLuint>::const_iterator it = uniform_location_map.find(key);
   if (it != uniform_location_map.end()) {
      return it->second;
   } else {
      GLint l = glGetUniformLocation(program_id, key.c_str());
      if (l == -1)
         if (false)
            std::cout << "INFO/WARNING:: " << name << " can't get a uniform location for " << key << std::endl;
      uniform_location_map[key] = l;
      return l;
   }
}

void Shader::set_uniform_locations() {
   GLuint err;

   if (entity_type == Entity_t::MODEL ||
       entity_type == Entity_t::MAP ||
       entity_type == Entity_t::MOLECULAR_TRIANGLES ||
       entity_type == Entity_t::GENERIC_DISPLAY_OBJECT) {
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


void
Shader::set_float_for_uniform(const std::string &u_name, float f) {

   GLuint idx = glGetUniformLocation_internal(u_name);
   GLenum err = glGetError(); if (err) std::cout << "error:: set_float_for_uniform() error 1a: "
                                                 << err << std::endl;
   glUniform1f(idx, f);
   err = glGetError(); if (err) std::cout << "error:: set_float_for_uniform() error 1b: "
                                          << err << std::endl;
}

void
Shader::set_vec4_for_uniform(const std::string &u_name, float f0, float f1, float f2, float f3) {

   GLuint idx = glGetUniformLocation_internal(u_name);
   float v[4];
   v[0] = f0; v[1] = f1; v[2] = f2; v[3] = f3;
   glUniform4fv(idx, 1, v);
}

void
Shader::set_vec4_for_uniform(const std::string &u_name, const glm::vec4 &v) {

   GLuint idx = glGetUniformLocation_internal(u_name);
   glUniform4fv(idx, 1, glm::value_ptr(v));
}

void
Shader::set_vec3_for_uniform(const std::string &u_name, const glm::vec3 &v) {

   GLuint idx = glGetUniformLocation_internal(u_name);
   glUniform3fv(idx, 1, glm::value_ptr(v));
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

void
Shader::setup_light(unsigned int light_index, const gl_lights_info_t &light,
                    const glm::mat4 &wrtm) {

   std::string s = "light_sources[" + std::to_string(light_index) + std::string("]");
   std::string a;

   a = s + ".is_on";
   set_bool_for_uniform(a, light.is_on);
   a = s + ".ambient";
   set_vec4_for_uniform(a, light.ambient);
   a = s + ".diffuse";
   set_vec4_for_uniform(a, light.diffuse);
   a = s + ".specular";
   set_vec4_for_uniform(a, light.specular);

   // the lights are in view coordinates (1,1,2) say, they need to move as the
   // world is rotated by the mouse (analogous to the eye position, which also
   // moves in the world coordinates when the view is rotated by the mouse

   glm::mat4 iwrtm = glm::inverse(wrtm);
   glm::vec4 p4 = light.position * iwrtm;
   glm::vec4 p3 = p4 / p4.w;

   p3 = light.position;
   p3.z = 1;
   p3.x = 1;

   std::cout << "sending light position "  << glm::to_string(p3) << std::endl;

   // std::cout << "lights: using wrm "  << glm::to_string(wrm) << std::endl;

   a = s + ".position";
   set_vec4_for_uniform(a, p3);

}
