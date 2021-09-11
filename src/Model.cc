//
#include <iostream>

#ifdef USE_ASSIMP
#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing flags
#endif

#include "Model.hh"

Model::Model(const std::string &obj_file_name) {

#if USE_ASSIMP
   std::cout << "file " << obj_file_name << std::endl;
   assimp_import(obj_file_name);
   draw_this_model = true;
#else
   std::cout << "Model not compiled with ASSIMP" << std::endl;
#endif

}

Model::Model(const std::vector<molecular_triangles_mesh_t> &mtm,
             Shader *shader_p, Material material,
             GtkWidget *gl_area) { // remove arg

   // make a mesh for each of the elements of the molecular_triangles_mesh
   draw_this_model = true;

   for (unsigned int i=0; i<mtm.size(); i++) {
      const molecular_triangles_mesh_t &m = mtm[i];
      gtk_gl_area_attach_buffers(GTK_GL_AREA(gl_area)); // no use.
      Mesh mesh(m);
      // mesh.setup(shader_p, material); 20210910-PE
      mesh.setup(material);
      meshes.push_back(mesh);
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

#include <sys/stat.h>

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

std::string
Model::append_dir_file(const std::string &directory, const std::string &fn) {

   std::string full_path = fn;
   if (! directory.empty())
      full_path = directory + "/" + fn;
   return full_path;
}

void
Model::draw_meshes(Shader *shader_p,  // e.g. molecular_triangles_shader
                   const glm::mat4 &mvp,
                   const glm::mat4 &view_rotation_matrix,
                   const std::map<unsigned int, lights_info_t> &lights,
                   const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                   const glm::vec4 &background_colour,
                   bool do_depth_fog) {

   if (! draw_this_model) return;

   for (unsigned int i=0; i<meshes.size(); i++) {
      // std::string mesh_name = meshes[i].name;
      meshes[i].draw(shader_p, mvp, view_rotation_matrix, lights,
                     eye_position, background_colour, do_depth_fog);
   }

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
