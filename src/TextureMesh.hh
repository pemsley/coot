/*
 * src/TextureMesh.hh
 *
 * Copyright 2020 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#ifndef TEXTURE_MESH_HH
#define TEXTURE_MESH_HH

#include <chrono>
#include <vector>
#include <string>
#include <epoxy/gl.h>
#include <glm/glm.hpp>
#include "coot-utils/g_triangle.hh"

#include "obj_loader.h"

#include "Shader.hh"
#include "Texture.hh" // now TextureMesh contains a vector of Textures - I am not sure this is a good
#include "stereo-eye.hh"
                      // arrangement (it is a result of load_from_glTF() which also fill/creates textures)

class TextureMeshVertex {
public:
   glm::vec3 position;
   glm::vec3 normal;
   glm::vec3 tangent;   // calculated later, using triangles
   glm::vec3 bitangent; // (same)
   glm::vec4 color;
   glm::vec2 texCoord;
   TextureMeshVertex(const glm::vec3 &p, const glm::vec3 &n, const glm::vec4 &col, const glm::vec2 &tc) :
      position(p), normal(n), color(col), texCoord(tc) {}
};

class TextureInfoType {
public:
   Texture texture;
   std::string name;
   // std::string sampler_name; // e.g. "base_texture", sampler is now in Texture::type
   GLuint unit;
   // TextureInfoType(const Texture &t, const std::string &n) : texture(t), name(n) {}
   TextureInfoType(const Texture &t, const std::string &n) {
      GLenum err = glGetError();
      if (err) std::cout << "GL ERROR:: TextureInfoType() A " << (err) << "\n";
      texture = t;
      err = glGetError();
      if (err) std::cout << "GL ERROR:: TextureInfoType() B " << (err) << "\n";
      name = n;
      err = glGetError();
      if (err) std::cout << "GL ERROR:: TextureInfoType() C " << (err) << "\n";
   }
};

class TextureMesh {
   enum { VAO_NOT_SET = 99999999 };
   GLuint vao;
   GLuint buffer_id;
   GLuint index_buffer_id;
   std::vector<TextureMeshVertex> vertices;
   std::vector<g_triangle> triangles;
   std::string name;
   std::string file_name;

   int n_instances; // number of instances to be _drawn_
   int n_instances_allocated; // that we made space for in glBufferData()
   bool is_instanced;
   // note: the instanced data for happy-face-residue-markers is just the position
   // (at least until this gets combined with sad-face-residue-markers/textures).
   // A combined (next to each other?) texture? But do we need that microoptimization?
   // The opacity will be a uniform.
   unsigned int draw_count; // so that I can animate the happy faces depending on the draw_count
   unsigned int inst_positions_id;
   static std::string _(int err);
   std::chrono::time_point<std::chrono::system_clock>  time_constructed;
   bool do_animation;
   float animation_A; // amplitude
   float animation_k; // wave number (scaled by position along the spine)
   float animation_w; // frequency

public:
   TextureMesh() : vao(VAO_NOT_SET), index_buffer_id(VAO_NOT_SET), draw_this_mesh(true) {
      n_instances_allocated = 0;
      n_instances = 0;
      is_instanced = false;
      inst_positions_id = -1;
      draw_count = 0;
      std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
      time_constructed = now;
      do_animation = false;
      set_animation_paramaters(0.05, 1.0, 1.0);
   }
   explicit TextureMesh(const std::string &n):
      vao(VAO_NOT_SET), index_buffer_id(VAO_NOT_SET), name(n), draw_this_mesh(true) {
      n_instances_allocated = 0;
      n_instances = 0;
      is_instanced = false;
      inst_positions_id = -1;
      draw_count = 0;
      std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
      time_constructed = now;
      do_animation = false;
      set_animation_paramaters(0.05, 1.0, 1.0);
   }
   bool draw_this_mesh;
   std::vector<TextureInfoType> textures;
   void import(const IndexedModel &ind_model, float scale);
   void import(const std::vector<TextureMeshVertex> &vertices, const std::vector<g_triangle> &triangles_in);
   bool have_instances() const { return is_instanced; }
   int get_n_instances() const { return n_instances; }
   void setup_tbn(unsigned int n_vertices); // tangent bitangent normal, pass the n_vertices for validation of indices.
   void setup_camera_facing_quad(float scale_x, float scale_y, float offset_x, float offset_y);
   // for the Z section x-offset_y and y_offset should be zero
   // but for the other sections they need to be offset so the layout looks pretty.
   void setup_tomo_quad(float x_scale, float y_scale, float x_offset, float y_offset,float z_pos, bool texture_x_y_swap_flag);
   void setup_buffers();
   void set_colour(const glm::vec4 &col_in);
   void setup_instancing_buffers(unsigned int n_happy_faces_max); // setup the buffer, don't add data
   std::string get_name() const { return name; }
   // this is for an ephemeral instanced texturemesh
   void update_instancing_buffer_data_for_happy_faces(const std::vector<glm::vec3> &positions, // in 3D space (of the CAs)
                                                      unsigned int draw_count_in,
                                                      unsigned int draw_count_max,
                                                      const glm::vec3 &screen_y_uv);
   void update_instancing_buffer_data(const std::vector<glm::vec3> &positions);
   void add_texture(const TextureInfoType &ti) { textures.push_back(ti); }
   void translate(const glm::vec3 &t); // include call to setup_buffers();
   void apply_scale(const float &sf); // include call to setup_buffers();
   void apply_transformation(const glm::mat4 &m);  // transform the positions in the vertices
   void draw(Shader *shader,
             const glm::mat4 &mvp,
             const glm::mat4 &view_rotation_matrix,
             const std::map<unsigned int, lights_info_t> &lights,
             const glm::vec3 &eye_position, // eye position in view space (not molecule space)
             const glm::vec4 &background_colour,
             bool do_depth_fog);
   void draw_with_shadows(Shader *shader,
                          const glm::mat4 &mvp,
                          const glm::mat4 &view_rotation_matrix,
                          const std::map<unsigned int, lights_info_t> &lights,
                          const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                          const glm::vec4 &background_colour,
                          bool do_depth_fog,
                          const glm::mat4 &light_space_mvp,
                          unsigned int depthMap,
                          float shadow_strength,
                          unsigned int shadow_softness,
                          bool show_just_shadows);
   void draw_for_ssao(Shader *shader_for_tmeshes_p,
                      const glm::mat4 &model,
                      const glm::mat4 &view,
                      const glm::mat4 &projection);
   void draw_atom_label(const std::string &atom_label,
                        const glm::vec3 &atom_label_position,
                        const glm::vec4 &text_colour, // set using subbufferdata
                        Shader *shader,
                        stereo_eye_t eye,
                        const glm::mat4 &mvp,
                        const glm::mat4 &view_rotation_matrix,
                        const glm::vec4 &background_colour,
                        bool do_depth_fog,
                        bool is_perspective_projection);
   // draw an ephemeral instanced opacity-varying texturemesh.
   // Other draw_instances() functions may be needed in future, if so change the name of this one.
   void draw_instances(Shader *shader_p, const glm::mat4 &mvp, const glm::mat4 &view_rotation,
                       const glm::vec4 &bg_col, bool is_perspective_projection);

   void draw_instances_for_ssao(Shader *shader_p,
                                const glm::mat4 &model,
                                const glm::mat4 &view,
                                const glm::mat4 &projection);

   void draw_fading_instances(Shader *shader_p, const glm::mat4 &mvp, const glm::mat4 &view_rotation,
                              unsigned int draw_count, unsigned int draw_count_max);

   bool load_from_glTF(const std::string &file_name, bool include_call_to_setup_buffers=true);

   // Holy Zarquon swimming fish
   void set_animation_paramaters(float amplitide_overall, float wave_number, float freq);

   //! return the time in milliseconds since this Model was constructed
   float duration() const {
      std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
      float time = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_constructed).count();
      return time;
   }

   void set_do_animation(bool state) { do_animation = state; }
};

#endif // TEXTURE_MESH_HH
