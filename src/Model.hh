
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
   std::vector<TextureMesh> tmeshes;
   std::vector<Mesh> meshes;
   std::string directory;
   std::string append_dir_file(const std::string &directory, const std::string &tfn);
   void assimp_import(const std::string &file_name_in);

#ifdef USE_ASSIMP
   void import(const aiScene* scene);
   void processNode(aiNode *node, const aiScene *scene);
   TextureMesh processMesh(aiMesh *mesh, const aiScene *scene);
   std::vector<Texture> loadMaterialTextures(aiMaterial *mat, aiTextureType type,
                                             const std::string &typeName);
#endif

public:
   Model() { draw_this_model = true; }
   explicit Model(const std::string &obj_file_name);
   explicit Model(const std::vector<molecular_triangles_mesh_t> &mtm,
                  Shader *shader_p, Material material,
                  GtkWidget *gl_area); // remove this arg when bug is fixed.
   bool draw_this_model; // overall control
   void draw_meshes(Shader *shader_p, // for plain meshes (e.g. molecular triangles)
                    const glm::mat4 &mvp,
                    const glm::mat4 &view_rotation_matrix,
                    const std::map<unsigned int, lights_info_t> &lights,
                    const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                    const glm::vec4 &background_colour,
                    bool do_depth_fog);
   void draw_tmeshes(Shader *shader_tmeshes_p,
                     const glm::mat4 &mvp,
                     const glm::mat4 &view_rotation_matrix,
                     const std::map<unsigned int, lights_info_t> &lights,
                     const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                     const glm::vec4 &background_colour,
                     bool do_depth_fog);
   // debugging function
   void add_mesh(const Mesh &m) { meshes.push_back(m); }
   bool export_as_obj(const std::string &file_name) const;

};

#endif // MODEL_HH
