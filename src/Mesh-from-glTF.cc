
// Define these only in *one* .cc file.
#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
// #define TINYGLTF_NOEXCEPTION // optional. disable exception handling.
#include "tiny_gltf.h"

#include <iostream>
#include <iomanip>

#include "Mesh.hh"

void
Mesh::load_from_glTF(const std::string &file_name) {

   auto dump_primitive = [] (const tinygltf::Primitive &primitive) {
                            const std::vector<std::map<std::string, int> > targets = primitive.targets;
                            for (unsigned int i=0; i<targets.size(); i++) {
                               const auto &m = targets[i];
                               std::map<std::string, int>::const_iterator it;
                               for (it=m.begin(); it!=m.end(); ++it) {
                                  const std::string &key = it->first;
                                  int idx = it->second;
                                  std::cout << "vector: " << i << " key: " << key << " " << idx << std::endl;
                               }
                            }
                         };

   auto dump_bufferview = [] (tinygltf::Model &model,
                              const tinygltf::BufferView &buffer_view) {

                             const tinygltf::Buffer &buffer = model.buffers[buffer_view.buffer];
                             const std::string &name = buffer.name;
                             std::vector<unsigned char> data = buffer.data;
                             std::cout << "buffer \"" << name << "\" has " << data.size() << " data items" << std::endl;

                             std::map<std::string, int>::const_iterator it;
                             // argh.
                          };

   auto proc_mesh = [] (tinygltf::Model &model, const tinygltf::Mesh &mesh) {
                       for (unsigned int i=0; i<mesh.primitives.size(); i++) {
                          const tinygltf::Primitive &primitive = mesh.primitives[i];
                          if (primitive.indices < 0) return;
                          std::map<std::string, int>::const_iterator it;
                          for (it=primitive.attributes.begin(); it!=primitive.attributes.end(); ++it) {
                             const std::string &key = it->first;
                             int index = it->second;
                             if (index < 0) return;
                             const tinygltf::Accessor &accessor = model.accessors[index];
                             int size = -1;
                             if (accessor.type == TINYGLTF_TYPE_SCALAR) size = 1;
                             if (accessor.type == TINYGLTF_TYPE_VEC2)   size = 2;
                             if (accessor.type == TINYGLTF_TYPE_VEC3)   size = 3;
                             if (accessor.type == TINYGLTF_TYPE_VEC4)   size = 4;
                             if (size != -1) {
                                if ((key == "POSITION") ||
                                    (key == "NORMAL")   ||
                                    (key == "TEXCOORD_0")) {
                                   int byte_stride = accessor.ByteStride(model.bufferViews[accessor.bufferView]);
                                   if (byte_stride != -1) {
                                      size_t buffer_offset = accessor.byteOffset;
                                      // ... the data size have been set/are accessible
                                      // glVertexAttribPointer here
                                      std::cout << "mesh primitive " << i << " key: " << std::setw(10) << key
                                                << " index:" << index << " type " << accessor.type
                                                << " byte_stride " << byte_stride
                                                << " buffer offset " << buffer_offset << std::endl;
                                   }
                                }
                             }
                          }

                          const tinygltf::Accessor &accessor_index = model.accessors[primitive.indices];
                          int mode = -1;
                          if (primitive.mode == TINYGLTF_MODE_TRIANGLES)      mode = GL_TRIANGLES;
                          if (primitive.mode == TINYGLTF_MODE_TRIANGLE_STRIP) mode = GL_TRIANGLE_STRIP;
                          if (primitive.mode == TINYGLTF_MODE_TRIANGLE_FAN)   mode = GL_TRIANGLE_FAN;
                          if (primitive.mode == TINYGLTF_MODE_POINTS)         mode = GL_POINTS;
                          if (primitive.mode == TINYGLTF_MODE_LINE)           mode = GL_LINES;
                          if (primitive.mode == TINYGLTF_MODE_LINE_LOOP)      mode = GL_LINE_LOOP;
                          if (mode != -1) {
                             // glDrawElements(), ie, make triangles
                          }
                       }
                   };
      

   auto proc_node = [proc_mesh] (tinygltf::Model &model, const tinygltf::Node &node) {
                       if (node.mesh > -1)  {
                          proc_mesh(model, model.meshes[node.mesh]);
                       }
                   };

   auto proc_model = [proc_node]  (tinygltf::Model &model) {

                        for (unsigned int iscene=0; iscene<model.scenes.size(); iscene++) {
                           const tinygltf::Scene &scene = model.scenes[iscene];
                           for (size_t i = 0; i < scene.nodes.size(); i++) {
                              proc_node(model, model.nodes[scene.nodes[i]]);
                           }
                        }
                     };

   auto setup_mesh_state = [] (tinygltf::Model &model) {

                              int bufferviews_size = model.bufferViews.size();
                              for (int i=0; i<bufferviews_size; i++) {
                                 const tinygltf::BufferView &bufferView = model.bufferViews[i];

                                 if (bufferView.target == 0) continue;
                                 int sparse_accessor = -1;
                                 for (unsigned int a_i=0; a_i<model.accessors.size(); a_i++) {
                                    const auto &accessor = model.accessors[a_i];
                                    if (accessor.bufferView == i) {
                                       // std::cout << "Index " << i << " is used by accessor " << a_i << std::endl;
                                       if (accessor.sparse.isSparse) {
                                          std::cout << "WARNING:: this bufferView has at least one sparse accessor" 
                                                    << std::endl;
                                          sparse_accessor = a_i;
                                          break;
                                       }

                                       const tinygltf::Buffer &buffer = model.buffers[bufferView.buffer];
                                       const float *positions = reinterpret_cast<const float *>(&buffer.data[bufferView.byteOffset + accessor.byteOffset]);
                                       std::cout << "accessor.count " << accessor.count << std::endl;
                                       int size = -1;
                                       if (accessor.type == TINYGLTF_TYPE_SCALAR) size = 1;
                                       if (accessor.type == TINYGLTF_TYPE_VEC2)   size = 2;
                                       if (accessor.type == TINYGLTF_TYPE_VEC3)   size = 3;
                                       if (accessor.type == TINYGLTF_TYPE_VEC4)   size = 4;
                                       if (accessor.type == TINYGLTF_TYPE_VEC3) {
                                          for (size_t ii=0; ii<accessor.count; ii++) {
                                             float x = positions[ii*3 + 0];
                                             float y = positions[ii*3 + 1];
                                             float z = positions[ii*3 + 2];
                                             std::cout << "VEC3 bufferview " << i << " accessor " << a_i << " ii " << ii
                                                       << " accessor-type " << accessor.type << " size " << size 
                                                       << " x " << x << " y " << y << " z " << z << std::endl;
                                          }
                                       }
                                       if (accessor.type == TINYGLTF_TYPE_VEC2) {
                                          for (size_t ii=0; ii<accessor.count; ii++) {
                                             float x = positions[ii*3 + 0];
                                             float y = positions[ii*3 + 1];
                                             std::cout << "VEC2 bufferview " << i << " accessor " << a_i << " ii " << ii
                                                       << " accessor-type " << accessor.type << " size " << size 
                                                       << " x " << x << " y " << y << std::endl;
                                          }
                                       }
                                       if (accessor.type == TINYGLTF_TYPE_SCALAR) {
                                          const int *indices = reinterpret_cast<const int *>(&buffer.data[bufferView.byteOffset + accessor.byteOffset]);
                                          for (size_t ii=0; ii<accessor.count; ii++) {
                                             float idx = indices[ii];
                                             std::cout << "SCALAR bufferview " << i << " accessor " << a_i << " ii " << ii
                                                       << " accessor-type " << accessor.type << " size " << size << " idx " << idx << std::endl;
                                          }
                                       }
                                    }
                                 }
                              }
                           };

   tinygltf::Model model;
   tinygltf::TinyGLTF loader;
   std::string err;
   std::string warn;

   // use the extension to check which function to use
   //
   // bool ret = loader.LoadASCIIFromFile(&model, &err, &warn, file_name);
   bool ret = loader.LoadBinaryFromFile(&model, &err, &warn, file_name); // for binary glTF(.glb)

   if (!warn.empty()) {
      std::cout << "WARNING:: " << warn << std::endl;
   }

   if (!err.empty()) {
      std::cout << "ERROR:: " << err << std::endl;
   }

   if (!ret) {
      std::cout << "WARNING:: failed to parse glTF from " << file_name << std::endl;
   }

   if (err.empty() && ret) {
      // OK, so what's in the tinygltf model?
      if (! model.meshes.empty()) {

         for (unsigned int i=0; i<model.meshes.size(); i++) {
            std::cout  << "Model mesh name     : " << model.meshes[i].name << std::endl;
            for (unsigned int k=0; k<model.meshes[i].primitives.size(); k++) {
               dump_primitive(model.meshes[i].primitives[k]);
            }
         }
      }

      if (! model.bufferViews.empty()) {
         for (unsigned int i=0; i<model.bufferViews.size(); i++) {
            const tinygltf::BufferView &buffer_view = model.bufferViews[i];
            dump_bufferview(model, buffer_view);
         }
      }
   }

   std::cout << "proc_model" << std::endl;
   proc_model(model);
   std::cout << "setup_mesh_state" << std::endl;
   setup_mesh_state(model);

}
