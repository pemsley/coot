
#ifndef MODEL_HH
#define MODEL_HH

#include <string>

#include <gtk/gtk.h>
#include <epoxy/gl.h>

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
   Model() { draw_this_model = true; }
   explicit Model(const std::string &obj_file_name);
   explicit Model(const std::vector<molecular_triangles_mesh_t> &mtm,
                  Material material, GtkWidget *gl_area);
   explicit Model(const Mesh &m) { draw_this_model = true; meshes.push_back(m); }

   // std::vector<std::pair<TextureMesh, std::string> > tmeshes;
   std::vector<TextureMesh> tmeshes;
   std::vector<Mesh> meshes;
   // should models have names?

   void add_mesh(const Mesh &m) { meshes.push_back(m); }
   void add_tmesh(const TextureMesh &tm) { tmeshes.push_back(tm); }

   // ---- drawing ----
   bool draw_this_model; // overall control
   void draw_mesh(unsigned int mesh_index,
		  Shader *shader_p, // for plain meshes (e.g. molecular triangles)
                  const glm::mat4 &mvp,
                  const glm::mat4 &view_rotation_matrix,
                  const std::map<unsigned int, lights_info_t> &lights,
                  const glm::vec3 &eye_position, // eye position in view space (not molecule space)
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
                    const glm::mat4 &mvp,
                    const glm::mat4 &view_rotation_matrix,
                    const std::map<unsigned int, lights_info_t> &lights,
                    const glm::vec3 &eye_position, // eye position in view space (not molecule space)
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

};

#endif // MODEL_HH
