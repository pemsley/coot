/*
 * src/Model.cc
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
//
#include <iostream>
#include "stereo-eye.hh"
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>  // to_string()

#ifdef USE_ASSIMP
#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing flags
#endif

#include "Model.hh"

#ifdef USE_ASSIMP
Model::Model(const std::string &obj_file_name) {

   std::cout << "file " << obj_file_name << std::endl;
   assimp_import(obj_file_name);
   draw_this_model = true;
   time_constructed = std::chrono::system_clock::now();
   do_animation = false;

}
#endif

Model::Model(const std::vector<molecular_triangles_mesh_t> &mtm,
             Material material, GtkWidget *gl_area) {

   // make a mesh for each of the elements of the molecular_triangles_mesh
   draw_this_model = true;

   for (unsigned int i=0; i<mtm.size(); i++) {
      const molecular_triangles_mesh_t &m = mtm[i];
      gtk_gl_area_attach_buffers(GTK_GL_AREA(gl_area));
      Mesh mesh(m);
      mesh.setup(material);
      meshes.push_back(mesh);
   }
   time_constructed = std::chrono::system_clock::now();
   do_animation = false;
}


void
Model::translate(const glm::vec3 &t) {

   for (auto &tm : tmeshes) {
      tm.translate(t);
   }
   for (auto &m : meshes) {
      m.translate_by(t); // nice and conistent. Hmm.
   }
}

void
Model::scale(const float &sf) {

   for (auto &tm : tmeshes) {
      tm.apply_scale(sf);
   }
   for (auto &m : meshes) {
      m.apply_scale(sf);
   }
}


#ifdef USE_ASSIMP
void
Model::assimp_import(const std::string &file_name_in) {

   std::string file_name = file_name_in;
   // file_name = "cube.obj";
   // file_name = "cessna.obj";

   Assimp::Importer importer;
   // And have it read the given file with some example postprocessing
   // Usually - if speed is not the most important aspect for you - you'll
   // probably to request more postprocessing than we do in this example.
   const aiScene* scene = importer.ReadFile(file_name,
                                            aiProcess_CalcTangentSpace       |
                                            aiProcess_Triangulate            |
                                            aiProcess_GenNormals             |
                                            aiProcess_JoinIdenticalVertices  |
                                            aiProcess_SortByPType);


   // If the import failed, report it
   if( !scene) {
      std::cout << "Error in read of " <<  file_name << " " << importer.GetErrorString()
                << std::endl;
      scene = 0;
   } else {
      //happy path
      std::cout << "------------ scene read OK from " << file_name << std::endl;
      // directory = file_name_in.substr(0, file_name_in.find_last_of('/'));
      std::string::size_type pos = file_name.find_last_of('/');
      if (pos != std::string::npos) {
         directory = file_name.substr(0, file_name.find_last_of('/'));
      }
      processNode(scene->mRootNode, scene);
   }
}
#endif

#ifdef USE_ASSIMP
void
Model::processNode(aiNode *node, const aiScene *scene) {

   std::cout << "processNode() " << node->mName.C_Str() << std::endl;

    for(unsigned int i = 0; i < node->mNumMeshes; i++) {
        aiMesh *mesh = scene->mMeshes[node->mMeshes[i]];
        tmeshes.push_back(processMesh(mesh, scene));
    }

    // recuse if needed
    for(unsigned int i = 0; i < node->mNumChildren; i++) {
       processNode(node->mChildren[i], scene);
    }
}
#endif

#ifdef USE_ASSIMP
TextureMesh
Model::processMesh(aiMesh *assimp_mesh, const aiScene *scene) {

   TextureMesh tmesh;
   bool get_textures = false;

   if (scene == 0) {
      std::cout << "null scene" << std::endl;
      return tmesh;
   }

   std::vector<TextureMeshVertex> vertices;
   std::vector<g_triangle> triangles;
   std::vector<Texture> textures;

   if (assimp_mesh->mNormals == nullptr) {
      std::cout << "assimp load provided no normals " << std::endl;
   }

   // --- vertices ---
   for (unsigned int i = 0; i < assimp_mesh->mNumVertices; i++) {
      const auto &v = assimp_mesh->mVertices[i];
      glm::vec3 normal(0,0,1);
      glm::vec3 pos(v.x, v.y, v.z);
      glm::vec4 col(0.8, 0.9, 0.99, 1.0);
      glm::vec2 tc(0,0);
      if (assimp_mesh->mNormals) {
         const auto &n = assimp_mesh->mNormals[i];
         normal = glm::vec3(n.x, n.y, n.z);
      }
      if (assimp_mesh->mTextureCoords[0]) {
         tc.x = assimp_mesh->mTextureCoords[0][i].x;
         tc.y = assimp_mesh->mTextureCoords[0][i].y;
      }
      TextureMeshVertex vertex(pos, normal, col, tc);
      vertices.push_back(vertex);
   }

   // --- indices ---
   for (unsigned int i = 0; i < assimp_mesh->mNumFaces; i++) {
      const aiFace &face = assimp_mesh->mFaces[i];
      g_triangle t(face.mIndices[0], face.mIndices[1], face.mIndices[2]);
      triangles.push_back(t);
   }

   if (get_textures) {
      // -- textures ---
      aiMaterial *material = scene->mMaterials[assimp_mesh->mMaterialIndex];
      if (material) {
         std::vector<Texture> diffuse_maps =
            loadMaterialTextures(material, aiTextureType_DIFFUSE, "texture_diffuse");
         std::vector<Texture> specular_maps =
            loadMaterialTextures(material, aiTextureType_SPECULAR, "texture_specular");
         std::vector<Texture> normal_maps =
            loadMaterialTextures(material, aiTextureType_NORMALS, "texture_normal");
         std::vector<Texture> ao_maps =
            loadMaterialTextures(material, aiTextureType_AMBIENT_OCCLUSION, "texture_ao");
         textures.insert(textures.end(),  diffuse_maps.begin(),  diffuse_maps.end());
         textures.insert(textures.end(), specular_maps.begin(), specular_maps.end());
         textures.insert(textures.end(),   normal_maps.begin(),   normal_maps.end());
         textures.insert(textures.end(),       ao_maps.begin(),       ao_maps.end());
      }
   }

   tmesh.import(vertices, triangles);
   tmesh.setup_buffers();

   return tmesh;
}
#endif

std::string
Model::append_dir_file(const std::string &directory, const std::string &fn) {

   std::string full_path = fn;
   if (! directory.empty())
      full_path = directory + "/" + fn;
   return full_path;
}

#include <sys/stat.h>

#ifdef USE_ASSIMP
std::vector<Texture>
Model::loadMaterialTextures(aiMaterial *mat,
                            aiTextureType type,
                            const std::string &typeName) {

   std::vector<Texture> textures;

   for (unsigned int i=0; i<mat->GetTextureCount(type); i++) {
      aiString texture_file_name;
      mat->GetTexture(type, i, &texture_file_name);
      std::string tfn = texture_file_name.C_Str();
      std::string full_path = append_dir_file(directory, tfn);
      // std::cout << "here with tfn " << tfn << " directory \"" << directory << "\"" << std::endl;

      bool file_exists = true;
      struct stat buffer;
      if (stat (full_path.c_str(), &buffer) != 0)
         file_exists = false;
      if (! file_exists) {
         std::cout << "failed to find texture image file " << full_path << std::endl;
      } else {
         // std::cout << "getting texture! " << tfn << " " << directory << std::endl;
         Texture texture;
         texture.init(tfn, directory);
         texture.set_type(typeName);
         texture.set_file_name(tfn);
         textures.push_back(texture);
      }
   }
   return textures;
}

#endif // USE_ASSIMP


void
Model::draw_meshes(Shader *shader_p,  // e.g. molecular_triangles_shader
                   stereo_eye_t eye,
                   const glm::mat4 &mvp,
                   const glm::mat4 &view_rotation_matrix,
                   const std::map<unsigned int, lights_info_t> &lights,
                   const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                   const glm::vec3 &rotation_centre,
                   float opacity,
                   const glm::vec4 &background_colour,
                   bool do_depth_fog) {

   if (! draw_this_model) return;

   // meshes should be draw opaque, right?
   glDisable(GL_BLEND);

   for (unsigned int i=0; i<meshes.size(); i++) {
      if (false) {
         const std::string &mesh_name = meshes[i].name;
         std::cout << "draw_meshes() drawing mesh " << i << " " << mesh_name << std::endl;
      }
      bool draw_just_shadows = false; // pass this if needed.
      bool wireframe_mode = false;
      meshes[i].draw(shader_p, eye, mvp, view_rotation_matrix, lights, eye_position, rotation_centre,
                     opacity, background_colour, wireframe_mode, do_depth_fog, draw_just_shadows);
      if (false)
         meshes[i].draw_normals(mvp, 0.1);
   }
}

void
Model::draw_mesh(unsigned int mesh_index,
                 Shader *shader_p,
                 stereo_eye_t eye,
                 const glm::mat4 &mvp,
                 const glm::mat4 &view_rotation_matrix,
                 const std::map<unsigned int, lights_info_t> &lights,
                 const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                 const glm::vec3 &rotation_centre,
                 float opacity,
                 const glm::vec4 &background_colour,
                 bool do_depth_fog,
		 bool draw_just_shadows) {

   glDisable(GL_BLEND);
   // std::cout << "Model draw_mesh() " << mesh_index << " \"" << meshes[mesh_index].name << "\"" << std::endl;
   bool wireframe_mode = false;
   meshes[mesh_index].draw(shader_p, eye, mvp, view_rotation_matrix, lights, eye_position, rotation_centre,
                           opacity, background_colour, wireframe_mode, do_depth_fog, draw_just_shadows);
}

void
Model::draw_tmesh(unsigned int mesh_index,
                 Shader *shader_p,
                 const glm::mat4 &mvp,
                 const glm::mat4 &view_rotation_matrix,
                 const std::map<unsigned int, lights_info_t> &lights,
                 const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                 const glm::vec4 &background_colour,
                 bool do_depth_fog) {

   tmeshes[mesh_index].draw(shader_p, mvp, view_rotation_matrix, lights, eye_position, background_colour, do_depth_fog);
}



void
Model::draw_with_shadows(Shader *shader_for_tmeshes_with_shadows_p,
			 Shader *shader_for_meshes_with_shadows_p,
			 const glm::mat4 &mvp,
			 const glm::mat4 &view_rotation,
			 const std::map<unsigned int, lights_info_t> &lights,
			 const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                         float opacity,
			 const glm::vec4 &bg_col,
			 bool do_depth_fog,
			 const glm::mat4 &light_view_mvp,
			 unsigned int depthMap_texture,
			 float shadow_strength,
                         unsigned int shadow_softness,  // 1, 2 or 3.
			 bool show_just_shadows) {

   // std::cout << "debug:: mvp in Model::draw_with_shadows() is " << glm::to_string(mvp) << std::endl;

   if (shader_for_tmeshes_with_shadows_p) {
      for (unsigned int j=0; j<tmeshes.size(); j++) {
	 draw_tmesh_with_shadows(j, shader_for_tmeshes_with_shadows_p,  mvp, view_rotation, lights, eye_position, // opacity
				 bg_col, do_depth_fog, light_view_mvp,
                                 depthMap_texture, shadow_strength, shadow_softness, show_just_shadows);
      }
   }

   // now the coloured vertices mesh (for molecule things, not textured things)

   // currently I won't draw the shadow, I think that it might be better to combine
   // textures later, rather than have each of the Mesh shaders/draw calls sample the shadow map
   //
   if (shader_for_meshes_with_shadows_p) {
      if (! meshes.empty()) {
	 for (unsigned int j=0; j<meshes.size(); j++) {
	    draw_mesh_with_shadows(j, shader_for_meshes_with_shadows_p,
				   mvp, view_rotation, lights, eye_position, opacity, bg_col, do_depth_fog,
				   light_view_mvp, depthMap_texture, shadow_strength, shadow_softness,
				   show_just_shadows);
	 }
      }
   }
}



void
Model::draw_tmeshes_with_shadows(Shader *shader_p,
				 const glm::mat4 &mvp,
				 const glm::mat4 &view_rotation_matrix,
				 const std::map<unsigned int, lights_info_t> &lights,
				 const glm::vec3 &eye_position, // eye position in view space (not molecule space)
				 const glm::vec4 &bg_col_v4,
				 bool do_depth_fog,
				 const glm::mat4 &light_view_mvp,
				 unsigned int shadow_depthMap_texture,
				 float shadow_strength,
				 unsigned int shadow_softness,
				 bool show_just_shadows) {

   if (! tmeshes.empty()) {
      for (unsigned int j=0; j<tmeshes.size(); j++) {
         float opacity = 1.0f;
	 draw_mesh_with_shadows(j, shader_p,  mvp, view_rotation_matrix, lights, eye_position, opacity,
				bg_col_v4, do_depth_fog, light_view_mvp,
                                shadow_depthMap_texture, shadow_strength, shadow_softness, show_just_shadows);
      }
   }
}


void
Model::draw_tmesh_with_shadows(unsigned int mesh_index,
                               Shader *shader_p,
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
                               bool show_just_shadows) {

   tmeshes[mesh_index].draw_with_shadows(shader_p, mvp, view_rotation_matrix, lights, eye_position,
                                         background_colour, do_depth_fog, light_view_mvp,
                                         depthMap, shadow_strength, shadow_softness,
                                         show_just_shadows);
}

void
Model::draw_mesh_with_shadows(unsigned int mesh_index, Shader *shader_p,
                              const glm::mat4 &mvp,
                              const glm::mat4 &view_rotation_matrix,
                              const std::map<unsigned int, lights_info_t> &lights,
                              const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                              float opacity,
                              const glm::vec4 &background_colour,
                              bool do_depth_fog,
                              const glm::mat4 &light_view_mvp,
                              unsigned int shadow_depthMap,
                              float shadow_strength,
                              unsigned int shadow_softness,
                              bool show_just_shadows) {
   

   // meshes should be draw opaque, right?
   glDisable(GL_BLEND);

   meshes[mesh_index].draw_with_shadows(shader_p, mvp, view_rotation_matrix, lights, eye_position, opacity,
					background_colour, do_depth_fog, light_view_mvp,
					shadow_depthMap, shadow_strength, shadow_softness, show_just_shadows);
}

void
Model::draw_tmeshes(Shader *shader_p,
                    const glm::mat4 &mvp,
                    const glm::mat4 &view_rotation_matrix,
                    const std::map<unsigned int, lights_info_t> &lights,
                    const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                    const glm::vec4 &background_colour,
                    bool do_depth_fog) {

   if (! draw_this_model) return;

   for (unsigned int i=0; i<tmeshes.size(); i++)
      tmeshes[i].draw(shader_p, mvp, view_rotation_matrix, lights,
                      eye_position, background_colour, do_depth_fog);
}

void
Model::draw_for_ssao(Shader *shader_for_tmeshes_p, Shader *shader_for_meshes_p,
                     const glm::mat4 &model, // matrix of course
                     const glm::mat4 &view,
                     const glm::mat4 &projection)  {

   if (shader_for_tmeshes_p)
      for (unsigned int i=0; i<tmeshes.size(); i++) {
         if (shader_for_tmeshes_p->unset_p()) {
            std::cout << "ERROR:: in draw_for_ssao() Ooopps! skipping draw_for_ssao() because shader_for_tmeshes_p is not setup "
                      << std::endl;
         } else {
            tmeshes[i].draw_for_ssao(shader_for_tmeshes_p, model, view, projection);
         }
      }

   if (shader_for_meshes_p) {
      if (shader_for_meshes_p->unset_p()) {
         std::cout << "ERROR:: in draw_for_ssao() Ooopps! skippping draw_for_ssao() because shader_for_meshes_p is not setup "
                   << std::endl;
      } else {
         for (unsigned int i=0; i<meshes.size(); i++) {
            meshes[i].draw_for_ssao(shader_for_meshes_p, model, view, projection);
         }
      }
   }

}


#include <fstream>

bool
Model::export_as_obj(const std::string &file_name) const {

   unsigned int status = true;

   std::cout << "Model::export_as_obj() exporting " << meshes.size() << " meshes " << std::endl;

   auto file_name_sans_extension = [] (const std::string &f) {
                                      std::string r = f;
                                      std::string::size_type iext = f.find_last_of(".");
                                      if (iext != std::string::npos)
                                         r = f.substr(0, iext);
                                      return r;
                                   };

   std::ofstream f(file_name.c_str());
   std::string stub = file_name_sans_extension(file_name);
   std::string mtl_file_name = stub + ".mtl";
   if (f) {
      f << "# Some model\n";
      f << "# " << "\n";
      f << "" << "\n";
      f << "mtllib ";
      f << mtl_file_name;
      f << "\n";
      f << "g exported_obj\n";

      unsigned int vertex_index_offset = 0;
      for (unsigned int i=0; i<meshes.size(); i++) {
         const Mesh &mesh(meshes[i]);
         // check the mesh.type here
         std::string use_mtl = "usemtl ";
         // std::cout << "Model::export_as_obj() mesh.type " << mesh.type << std::endl;
         if (mesh.type == 1)
            use_mtl += "red";
         if (mesh.type == 2)
            use_mtl += "green";
         if (mesh.type == 4)
            use_mtl += "blue";
         f << use_mtl;
         f << "\n";

         meshes[i].export_as_obj(f, vertex_index_offset);
         vertex_index_offset += meshes[i].vertices.size();
      }
      f.close();

      std::ofstream f_mtl(mtl_file_name.c_str());
      if (f_mtl) {
         f_mtl << "newmtl red\n";
         f_mtl << "  Ka 0.800 0.300 0.300\n";
         f_mtl << "  Kd 0.800 0.300 0.300\n";
         f_mtl << "  Ks 1.0 1.0 1.0\n";
         f_mtl << "  Ns 5.0\n";
         f_mtl << "  illum 2\n";
         f_mtl << "newmtl green\n";
         f_mtl << "  Ka 0.300 0.800 0.300\n";
         f_mtl << "  Kd 0.300 0.800 0.300\n";
         f_mtl << "  Ks 1.0 1.0 1.0\n";
         f_mtl << "  Ns 5.0\n";
         f_mtl << "  illum 2\n";
         f_mtl << "newmtl blue\n";
         f_mtl << "  Ka 0.300 0.300 0.800\n";
         f_mtl << "  Kd 0.300 0.300 0.800\n";
         f_mtl << "  Ks 1.0 1.0 1.0\n";
         f_mtl << "  Ns 5.0\n";
         f_mtl << "  illum 2\n";
      }
   }

   return status;
}

// 20250527-PE Holy Zarquon swimming fish
void
Model::set_animation_parameters(float amplitide_overall, float wave_number, float freq) {

   for (unsigned int i=0; i<tmeshes.size(); i++) {
      tmeshes[i].set_animation_paramaters(amplitide_overall, wave_number, freq);
   }
}

void
Model::set_do_animation(bool state) {

   do_animation = state;
   for (unsigned int i=0; i<tmeshes.size(); i++) {
      tmeshes[i].set_do_animation(state);
   }

}
