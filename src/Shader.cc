
#include <string>
#include <iostream>
#include <fstream>
#include <epoxy/gl.h>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>  // to_string()

#include "Shader.hh"

Shader::Shader() {
   program_id = 0; // unset
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
   name = file_name;
   std::cout << "::: Shader compile " << file_name << std::endl;

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
   if (err) std::cout << "error:: Shader::set_int_for_uniform() " << name << " start err " << err << std::endl;
   GLint loc = glGetUniformLocation_internal(uniform_name.c_str());
   err = glGetError(); if (err) std::cout << "error:: Shader::set_int_for_uniform() " << name
                                          << " A err " << err << std::endl;
   glUniform1i(loc,value);
   err = glGetError(); if (err) std::cout << "error:: Shader::set_int_for_uniform() " << name
                                          << " B glUniform1i for uniform " << uniform_name
                                          << " loc: " << loc << " value: " << value
                                          << " err " << err << std::endl;
}

void
Shader::set_bool_for_uniform(const std::string &uniform_name, bool value) {

   GLuint err = glGetError();
   if (err)
      std::cout << "error:: Shader::set_bool_for_uniform() " << name << " "
                << uniform_name << " start err " << err << std::endl;

   GLint loc = glGetUniformLocation_internal(uniform_name.c_str());
   // std::cout << "set_bool_for_uniform() got loc " << loc << std::endl;
   err = glGetError();
   if (err)
      std::cout << "ERROR:: " << name << " Shader::set_bool_for_uniform() " << uniform_name << " A err "
                << err << std::endl;
   glUniform1i(loc, value);
   err = glGetError();
   if (err)
      std::cout << "ERROR:: " << name << " Shader::set_bool_for_uniform() " << uniform_name << " B err "
                << err << std::endl;
}



void
Shader::Use() {
   GLuint err = glGetError();
   if (err) std::cout << "Shader::Use() " << name << " pre glUseProgram() err " << err << std::endl;
   glUseProgram(program_id);
   err = glGetError();
   if (err) std::cout << "Shader::Use() " << name << " err " << err
                      << " for program_id " << program_id << std::endl;
}

void
Shader::set_attribute_locations() {

   // Is this needed? instanced-balls was working without it.

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
   if (entity_type == Entity_t::GENERIC_DISPLAY_OBJECT) {
      glBindAttribLocation(program_id, 0, "position");
      glBindAttribLocation(program_id, 1, "normal");
      glBindAttribLocation(program_id, 2, "colour");
      glBindAttribLocation(program_id, 3, "model_translation");
   }
   if (entity_type == Entity_t::INSTANCED_DISPLAY_OBJECT) {
      glBindAttribLocation(program_id, 0, "position");
      glBindAttribLocation(program_id, 1, "normal");
      glBindAttribLocation(program_id, 2, "colour");
      glBindAttribLocation(program_id, 3, "model_rotation_translation_scale_0");
      glBindAttribLocation(program_id, 4, "model_rotation_translation_scale_1");
      glBindAttribLocation(program_id, 5, "model_rotation_translation_scale_2");
      glBindAttribLocation(program_id, 6, "model_rotation_translation_scale_3");
   }
   if (entity_type == Entity_t::HUD_TEXT) {
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
      // std::cout << "creating a new location for key " << key << " " << l << std::endl;
      return l;
   }
}

void Shader::set_uniform_locations() {
   GLuint err;

   if (entity_type == Entity_t::MODEL ||
       entity_type == Entity_t::MAP ||
       entity_type == Entity_t::MOLECULAR_TRIANGLES ||
       entity_type == Entity_t::INSTANCED_DISPLAY_OBJECT ||
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
      atom_label_projection_uniform_location = glGetUniformLocation_internal("projection");
      err = glGetError(); if (err) std::cout << "error:: set_uniform_locations() error 6a: " << err << std::endl;
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

void
Shader::set_vec2_for_uniform(const std::string &u_name, const glm::vec2 &v) {

   GLuint idx = glGetUniformLocation_internal(u_name);
   glUniform2fv(idx, 1, glm::value_ptr(v));
}

void Shader::set_more_uniforms_for_molecular_triangles() {

   // put more uniforms here
}

#include <sys/stat.h>

void Shader::parse(const std::string &file_name_in) {

   std::string file_name = file_name_in;
   bool file_exists = true;
   struct stat buffer;
   if (stat (name.c_str(), &buffer) != 0) file_exists = false;

   if (! file_exists)
      if (! default_directory.empty())
         file_name = default_directory + "/" + file_name;

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
            if (type == ShaderType::VERTEX) {
               VertexSource += line;
               VertexSource += "\n";
            }
            if (type == ShaderType::FRAGMENT) {
               FragmentSource += line;
               FragmentSource += "\n";
            }
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
Shader::setup_light(unsigned int light_index, const lights_info_t &light,
                    const glm::mat4 &vrm) {

   GLenum err = glGetError();
   if (err) std::cout << "error setup_light() " << name << " -- start -- " << err << std::endl;

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

   glm::mat4 ivrm = glm::inverse(vrm);
   glm::vec4 p4   = glm::vec4(light.direction,1.0) * vrm;
   glm::vec4 p4_i = glm::vec4(light.direction,1.0) * ivrm;
   glm::vec4 p4_wc(glm::vec3(p4_i / p4_i.w), 1.0);

   err = glGetError();
   if (err) std::cout << "error setup_light() " << light_index << " "
                      << name << " A " << err << std::endl;

   if (false)
      std::cout << "sending light direction_in_molecule_coordinates_space orig: "
                << glm::to_string(light.direction) << " now: "
                << glm::to_string(glm::vec3(p4)) << std::endl;

   a = s + ".direction";
   set_vec3_for_uniform(a, light.direction);

   err = glGetError();
   if (err) std::cout << "error setup_light() " << name << " B " << err << std::endl;

   a = s + ".direction_in_molecule_coordinates_space";
   set_vec3_for_uniform(a, glm::vec3(p4));
   err = glGetError();
   if (err) std::cout << "error setup_light() " << light_index << " "
                      << name << " -- end -- " << err << std::endl;

}
