/*
 * src/Model.hh
 *
 * Copyright 2021 by Medical Research Council
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

#ifndef MODEL_HH
#define MODEL_HH

#include <string>

#include <gtk/gtk.h>
#include <epoxy/gl.h>

#include "stereo-eye.hh"
#include "Shader.hh"
#include "Texture.hh"
#include "Mesh.hh"
#include "TextureMesh.hh"
#include "molecular-mesh-generator.hh"

#ifdef USE_ASSIMP
#include <assimp/scene.h>
#endif

class Model {
   // 20211023-PE in future, tmeshes will need to take a set/vector/map of Texture (or reference/indices)
   // for normal maps/specularity maps/occlusion maps and stuff.
   //
   std::string directory;
   std::string append_dir_file(const std::string &directory, const std::string &tfn);
   void assimp_import(const std::string &file_name_in);
   std::chrono::time_point<std::chrono::system_clock> time_constructed;

#ifdef USE_ASSIMP
   void import(const aiScene* scene);
   void assimp_import(const std::string &file_name_in);
   void processNode(aiNode *node, const aiScene *scene);
   TextureMesh processMesh(aiMesh *mesh, const aiScene *scene);
   std::vector<Texture> loadMaterialTextures(aiMaterial *mat, aiTextureType type,
                                             const std::string &typeName);
#endif

   void draw_tmeshes_with_shadows(Shader *shader_p,
				  const glm::mat4 &mvp,
				  const glm::mat4 &view_rotation_matrix,
				  const std::map<unsigned int, lights_info_t> &lights,
				  const glm::vec3 &eye_position, // eye position in view space (not molecule space)
				  const glm::vec4 &background_colour,
				  bool do_depth_fog,
				  const glm::mat4 &light_view_mvp,
				  unsigned int depthMap,
				  float shadow_strength,
                                  unsigned int shadow_softness,
				  bool show_Just_shadows);
   void draw_tmesh_with_shadows(unsigned int tmesh_index, Shader *shader_p, // for plain meshes (e.g. molecular triangles)
                                const glm::mat4 &mvp,
                                const glm::mat4 &view_rotation_matrix,
                                const std::map<unsigned int, lights_info_t> &lights,
                                const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                                const glm::vec4 &background_colour,
                                bool do_depth_fog,
                                const glm::mat4 &light_view_mvp,
                                unsigned int depthMap,
                                float shadow_strength,
                                unsigned int shadow_softness,
                                bool show_Just_shadows);
   void draw_mesh_with_shadows(unsigned int mesh_index, Shader *shader_p, // for plain meshes (e.g. molecular triangles)
                               const glm::mat4 &mvp,
                               const glm::mat4 &view_rotation_matrix,
                               const std::map<unsigned int, lights_info_t> &lights,
                               const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                               float opacity,
                               const glm::vec4 &background_colour,
                               bool do_depth_fog,
                               const glm::mat4 &light_view_mvp,
                               unsigned int depthMap,
                               float shadow_strength,
                               unsigned int shadow_softness,
                               bool show_Just_shadows);

public:
   Model() {
      draw_this_model = true;
      time_constructed = std::chrono::system_clock::now();
      do_animation = false;
   }
   explicit Model(const std::string &obj_file_name);
   explicit Model(const std::vector<molecular_triangles_mesh_t> &mtm,
                  Material material, GtkWidget *gl_area);
   explicit Model(const Mesh &m) {
      draw_this_model = true;
      meshes.push_back(m);
      time_constructed = std::chrono::system_clock::now();
      do_animation = false;
   }

   // std::vector<std::pair<TextureMesh, std::string> > tmeshes;
   std::vector<TextureMesh> tmeshes;
   std::vector<Mesh> meshes;
   // should models have names?

   bool do_animation; // so that we don't have to query the texture meshes
                      // i.e. if the model animates, the meshes animate.

   void add_mesh(const Mesh &m) { meshes.push_back(m); }
   void add_tmesh(const TextureMesh &tm) { tmeshes.push_back(tm); }

   // ---- drawing ----
   bool draw_this_model; // overall control
   void draw_mesh(unsigned int mesh_index,
		  Shader *shader_p, // for plain meshes (e.g. molecular triangles)
                  stereo_eye_t eye,
                  const glm::mat4 &mvp,
                  const glm::mat4 &view_rotation_matrix,
                  const std::map<unsigned int, lights_info_t> &lights,
                  const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                  const glm::vec3 &rotation_centre,
                  float opacity,
                  const glm::vec4 &background_colour,
                  bool do_depth_fog,
		  bool draw_just_shadows);
   void draw_tmesh(unsigned int tmesh_index, Shader *shader_p, // for texture meshes
                  const glm::mat4 &mvp,
                  const glm::mat4 &view_rotation_matrix,
                  const std::map<unsigned int, lights_info_t> &lights,
                  const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                  const glm::vec4 &background_colour,
                  bool do_depth_fog);
   void draw_with_shadows(Shader *tmeshes_shader_p,
			  Shader *meshes_shader_p,
			  const glm::mat4 &mvp,
			  const glm::mat4 &view_rotation_matrix,
			  const std::map<unsigned int, lights_info_t> &lights,
			  const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                          float opacity,
			  const glm::vec4 &background_colour,
			  bool do_depth_fog,
			  const glm::mat4 &light_view_mvp,
			  unsigned int depthMap,
			  float shadow_strength,
                          unsigned int shadow_softness, // 1, 2 or 3
			  bool show_just_shadows);
   void draw_for_ssao(Shader *shader_for_tmeshes_p, Shader *shader_for_meshes_p,
                      const glm::mat4 &model,
                      const glm::mat4 &view,
                      const glm::mat4 &projection);
   void draw_meshes(Shader *shader_p, // for plain meshes (e.g. molecular triangles)
                    stereo_eye_t eye,
                    const glm::mat4 &mvp,
                    const glm::mat4 &view_rotation_matrix,
                    const std::map<unsigned int, lights_info_t> &lights,
                    const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                    const glm::vec3 &rotation_centre,
                    float opacity,
                    const glm::vec4 &background_colour,
                    bool do_depth_fog);
   void draw_tmeshes(Shader *shader_tmeshes_p,
                     const glm::mat4 &mvp,
                     const glm::mat4 &view_rotation_matrix,
                     const std::map<unsigned int, lights_info_t> &lights,
                     const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                     const glm::vec4 &background_colour,
                     bool do_depth_fog);

   void translate(const glm::vec3 &t);
   void scale(const float &sf);
   // debugging function
   bool export_as_obj(const std::string &file_name) const;
   //! return the time in milliseconds since this Model was constructed
   float duration() const {
      std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
      float time = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_constructed).count();
      return time;
   }

   // 20250527-PE Holy Zarquon swimming fish
   void set_animation_parameters(float amplitude_overall, float wave_number, float freq);
   void set_do_animation(bool state);
};

#endif // MODEL_HH
