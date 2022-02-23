
// Define these only in *one* .cc file.

#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION


// #define TINYGLTF_NOEXCEPTION // optional. disable exception handling.
#include "tiny_gltf.h"

#include <iostream>
#include <iomanip>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>  // to_string()

#include "utils/coot-utils.hh"
#include "Mesh.hh"


// this should have its own file.
std::ostream& operator <<(std::ostream &s, const g_triangle &t) {
   s << "[" << t.point_id[0] << " " << t.point_id[1] << " " << t.point_id[2] << "]";
   return s;
}

// include_call_to_setup_buffers is by default true
bool
Mesh::load_from_glTF(const std::string &file_name_in, bool include_call_to_setup_buffers) {

   // is this what's in a node?
   class extracted_buffer_info_t {
   public:
      int buffer_view_index;
      std::vector<glm::vec3> positions;
      std::vector<glm::vec3> normals;
      std::vector<glm::vec4> vertex_colours;
      std::vector<glm::vec2> texture_coords;
      std::vector<g_triangle> triangles;
      glm::vec4 base_colour;
   };


   auto indent = [] (unsigned int il) {
                   std::string s;
                   for (unsigned int i=0; i<il; i++) s += "   ";
                   return s;
                   };

   auto proc_indices = [indent] (tinygltf::Model &model, const tinygltf::Accessor &indices_accessor) {

                          std::vector<g_triangle> triangles;
                          std::cout << indent(3) << "proc_indices()" << std::endl;

                          std::vector<unsigned int> mesh_indices;
                          const tinygltf::BufferView &buffer_view = model.bufferViews[indices_accessor.bufferView];
                          const tinygltf::Buffer &buffer = model.buffers[buffer_view.buffer];
                          const uint8_t *base = &buffer.data.at(buffer_view.byteOffset + indices_accessor.byteOffset);

                          if (indices_accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT) {
                             const uint32_t *p = reinterpret_cast<const uint32_t *>(base);
                             for (size_t i=0; i<indices_accessor.count; i++)
                                mesh_indices.push_back(p[i]);
                          }

                          if (indices_accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT) {
                             const uint16_t *p = reinterpret_cast<const uint16_t *>(base);
                             for (size_t i=0; i<indices_accessor.count; i++)
                                mesh_indices.push_back(p[i]);
                          }

                          if (indices_accessor.componentType == TINYGLTF_COMPONENT_TYPE_BYTE) {
                             const uint8_t *p = reinterpret_cast<const uint8_t *>(base);
                             for (size_t i=0; i<indices_accessor.count; i++)
                                mesh_indices.push_back(p[i]);
                          }

                          unsigned int n_triangles = mesh_indices.size()/3;
                          if (n_triangles * 3 == mesh_indices.size()) {
                             for (unsigned int ii=0; ii<mesh_indices.size(); ii+=3) {
                                g_triangle t(mesh_indices[ii], mesh_indices[ii+1], mesh_indices[ii+2]);
                                triangles.push_back(t);
                             }
                          }
                          return triangles;
                       };

   auto proc_material = [] (tinygltf::Model &model, int material_index) {

                           std::cout << "debug:: proc_material(): material_index " << material_index
                                     << " materials.size() " << model.materials.size() << std::endl;

                           const tinygltf::Material &material = model.materials[material_index];

                           std::cout << "debug:: proc_material(): material.pbrMetallicRoughness() basecolour size"
                                     << material.pbrMetallicRoughness.baseColorFactor.size() << std::endl;

                           if (false)
                              std::cout << "Material " << material_index << " name: \"" << material.name << "\""
                                        << " alphaMode" << material.alphaMode
                                        << " pbr-metallicroughness colour "
                                        << material.pbrMetallicRoughness.baseColorFactor[0] << " "
                                        << material.pbrMetallicRoughness.baseColorFactor[1] << " "
                                        << material.pbrMetallicRoughness.baseColorFactor[2] << " "
                                        << material.pbrMetallicRoughness.baseColorFactor[3] << " "
                                        << " metalicFactor " << material.pbrMetallicRoughness.metallicFactor
                                        << " roughnessFactor " << material.pbrMetallicRoughness.roughnessFactor
                                        << std::endl;
                           glm::vec4 colour(material.pbrMetallicRoughness.baseColorFactor[0],
                                            material.pbrMetallicRoughness.baseColorFactor[1],
                                            material.pbrMetallicRoughness.baseColorFactor[2],
                                            material.pbrMetallicRoughness.baseColorFactor[3]);
                           return colour;
                   };

   auto proc_primitive = [proc_indices, proc_material, indent] (tinygltf::Model &model, const tinygltf::Primitive &primitive) {

                            std::cout << indent(2) << "proc_primitive()" << std::endl;

                            extracted_buffer_info_t ebi; // returned object

                            // ------------------- primitive material -------------------

                            std::cout << indent(3) << "primitive material " << primitive.material << std::endl;
                            if (primitive.material >= 0) {
                               glm::vec4 colour = proc_material(model, primitive.material);
                               ebi.base_colour = colour;
                            } else {
			      std::cout << indent(3) << "info;: skipping proc_material()" << std::endl;
                            }

                            // ------------------- primitive attributes -------------------

                            std::map<std::string, int>::const_iterator it;
                            for (it=primitive.attributes.begin(); it!=primitive.attributes.end(); ++it) {
                               const std::string &key = it->first;
                               int index = it->second;
                               std::cout << indent(3) << "primvitive attribute " << key << " " << index << std::endl;

                               const tinygltf::Accessor &accessor = model.accessors[index];
                               const tinygltf::BufferView &buffer_view = model.bufferViews[accessor.bufferView];
                               const tinygltf::Buffer &buffer = model.buffers[buffer_view.buffer];
                               if (buffer_view.target == TINYGLTF_TARGET_ARRAY_BUFFER)
                                  std::cout << indent(3) << "an array buffer" << std::endl;
                               if (buffer_view.target == TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER)
                                  std::cout << indent(3) << "an index buffer" << std::endl;

                               if (false)
                                  std::cout << indent(2) << "accessor componentType " << accessor.componentType
                                            << " accessor type " << accessor.type << std::endl;

                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_BYTE)
                                  std::cout << "      component type byte " << std::endl;
                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE)
                                  std::cout << "      component type unsigned byte " << std::endl;
                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_SHORT)
                                  std::cout << "      component type short" << std::endl;
                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT)
                                  std::cout << "      component type unsigned short" << std::endl;
                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_INT)
                                  std::cout << "      component type int" << std::endl;
                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT)
                                  std::cout << "      component type unsigned int" << std::endl;
                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT)
                                  std::cout << "      component type float" << std::endl;
                               if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_DOUBLE)
                                  std::cout << "      component type double" << std::endl;

                               if (buffer_view.target == TINYGLTF_TARGET_ARRAY_BUFFER) {
                                  if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT) {
                                     if (key == "POSITION") {
                                        if (accessor.type == TINYGLTF_TYPE_VEC3) {
                                           // std::cout << "      ---- float array for positions " << std::endl;
                                           const float *positions = reinterpret_cast<const float *>(&buffer.data[buffer_view.byteOffset + accessor.byteOffset]);
                                           for (size_t ii=0; ii<accessor.count; ii++) {
                                              float x = positions[ii*3 + 0];
                                              float y = positions[ii*3 + 1];
                                              float z = positions[ii*3 + 2];
                                              ebi.positions.push_back(glm::vec3(x,y,z));
                                           }
                                        }
                                     }
                                  }
                               }

                               if (buffer_view.target == TINYGLTF_TARGET_ARRAY_BUFFER) {
                                  if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT) {
                                     if (key == "NORMAL") {
                                        if (accessor.type == TINYGLTF_TYPE_VEC3) {
                                           // std::cout << "      ---- float array for normals " << std::endl;
                                           const float *normals = reinterpret_cast<const float *>(&buffer.data[buffer_view.byteOffset + accessor.byteOffset]);
                                           for (size_t ii=0; ii<accessor.count; ii++) {
                                              float x = normals[ii*3 + 0];
                                              float y = normals[ii*3 + 1];
                                              float z = normals[ii*3 + 2];
                                              ebi.normals.push_back(glm::vec3(x,y,z));
                                           }
                                        }
                                     }
                                  }
                               }

			       std::cout << "Here A " << key << std::endl;
			       if (buffer_view.target == TINYGLTF_TARGET_ARRAY_BUFFER) {
				  std::cout << "Here B " << key << " " << accessor.componentType << std::endl;
				  if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT) {
				     std::cout << "Here C1 " << key << std::endl;
				     if (key =="COLOR_0") {
					std::cout << "Here D1 " << std::endl;
                                        if (accessor.type == TINYGLTF_TYPE_VEC4) {
                                           const float *cols = reinterpret_cast<const float *>(&buffer.data[buffer_view.byteOffset + accessor.byteOffset]);
                                           for (size_t ii=0; ii<accessor.count; ii++) {
                                              float r = cols[ii*4 + 0];
                                              float g = cols[ii*4 + 1];
                                              float b = cols[ii*4 + 2];
                                              float a = cols[ii*4 + 3];
					      // std::cout << "Found col " << r << " " << g << " " << b << std::endl;
					      glm::vec4 col(r,g,b,a);
					      ebi.vertex_colours.push_back(col);
					   }
					}
				     }
				  }
				  if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT) {
				     std::cout << "Here C2 " << key <<  " type unsigned short with type "
					       << accessor.type << std::endl;
				     if (accessor.type == TINYGLTF_TYPE_VEC4) {
					if (key == "COLOR_0") {
                                           const unsigned short *cols = reinterpret_cast<const unsigned short *>(&buffer.data[buffer_view.byteOffset + accessor.byteOffset]);
					   const size_t s = sizeof(unsigned short);
                                           for (size_t ii=0; ii<accessor.count; ii++) {
                                              int r = cols[ii*4 + 0];
                                              int g = cols[ii*4 + 1];
                                              int b = cols[ii*4 + 2];
                                              int a = cols[ii*4 + 3];
					      int max = 256 * 256 -1;
					      float rr = static_cast<float>(r)/static_cast<float>(max);
					      float gg = static_cast<float>(g)/static_cast<float>(max);
					      float bb = static_cast<float>(b)/static_cast<float>(max);
					      float aa = static_cast<float>(a)/static_cast<float>(max);
					      //std::cout << "Found col from 'unsigned short' "
                                              // << rr << " " << gg << " " << bb << " " << aa << std::endl;
					      glm::vec4 col(rr,gg,bb,aa);
					      ebi.vertex_colours.push_back(col);
					   }
					}
				     }
				  }
				  if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE) {
				     std::cout << "Here C2 " << key <<  " type unsigned byte with type "
					       << accessor.type << std::endl;
				     if (accessor.type == TINYGLTF_TYPE_VEC4) {
					if (key == "COLOR_0") {
                                           const unsigned char *cols = reinterpret_cast<const unsigned char *>(&buffer.data[buffer_view.byteOffset + accessor.byteOffset]);
					   const size_t s = sizeof(unsigned char);
                                           for (size_t ii=0; ii<accessor.count; ii++) {
                                              int r = cols[ii*4 + 0];
                                              int g = cols[ii*4 + 1];
                                              int b = cols[ii*4 + 2];
                                              int a = cols[ii*4 + 3];
					      int max = 256 -1;
					      float rr = static_cast<float>(r)/static_cast<float>(max);
					      float gg = static_cast<float>(g)/static_cast<float>(max);
					      float bb = static_cast<float>(b)/static_cast<float>(max);
					      float aa = static_cast<float>(a)/static_cast<float>(max);
					      // std::cout << "Found col from 'unsigned byte' "
                                              // << rr << " " << gg << " " << bb << " " << aa << std::endl;
					      glm::vec4 col(rr,gg,bb,aa);
					      ebi.vertex_colours.push_back(col);
					   }
					}
				     }
				  }
			       }

                               if (buffer_view.target == TINYGLTF_TARGET_ARRAY_BUFFER) {
                                  if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT) {
                                     if (key == "TEXCOORD_0") {
                                        if (accessor.type == TINYGLTF_TYPE_VEC2) {
                                           // std::cout << "      ---- float array for texture coordinates " << std::endl;
                                           const float *tcs = reinterpret_cast<const float *>(&buffer.data[buffer_view.byteOffset + accessor.byteOffset]);
                                           for (size_t ii=0; ii<accessor.count; ii++) {
                                              float x = tcs[ii*2 + 0];
                                              float y = tcs[ii*2 + 1];
                                              ebi.texture_coords.push_back(glm::vec2(x,y));
                                           }
                                        }
                                     }
                                  }
                               }
                            }

                            // ------------------- primitive indices -------------------

                            const tinygltf::Accessor &indices_accessor = model.accessors[primitive.indices];
                            std::vector<g_triangle> triangles = proc_indices(model, indices_accessor);
                            ebi.triangles = triangles; // std::move here?

                            std::cout << indent(2) << "debug:: proc_primitive() returns ebi with " << ebi.positions.size() << " positions "
                                      << ebi.normals.size() << " normals " << ebi.texture_coords.size() << " texture-coords and "
                                      << triangles.size() << " triangles" << std::endl;

                            return ebi;
                         };

   auto proc_mesh = [proc_primitive, indent] (tinygltf::Model &model, const tinygltf::Mesh &mesh) {

                       std::cout << indent(1) << "proc_mesh()" << std::endl;
                       std::vector<extracted_buffer_info_t> ebi_vec;
                       for (unsigned int i=0; i<mesh.primitives.size(); i++) {
                          std::cout << indent(1) << "proc_mesh() primitive number " << i << std::endl;
                          const tinygltf::Primitive &primitive = mesh.primitives[i];
                          if (primitive.indices >= 0) {
                             auto new_mesh_primitive = proc_primitive(model, primitive);
                             ebi_vec.push_back(new_mesh_primitive);
                          }
                       }
                       return ebi_vec;
                   };


   auto proc_node = [proc_mesh] (tinygltf::Model &model, size_t i_node_index, const tinygltf::Node &node) {

                       std::vector<extracted_buffer_info_t> r;

                       std::cout << "debug:: proc_node() Node " << i_node_index << " info:" << std::endl;
                       std::cout << "   Node name: " << node.name << std::endl;
                       std::cout << "   Node mesh index: " << node.mesh << std::endl;
                       std::cout << "   Node rotation vec elements count: " << node.rotation.size() << std::endl;
                       std::cout << "   Node scale vec elements count:    " << node.scale.size() << std::endl;
                       std::cout << "   Node scale vec translation count: " << node.translation.size() << std::endl;
                       if (node.mesh > -1) { // not a light or a camera
                          r = proc_mesh(model, model.meshes[node.mesh]);
                       } else {
                       }
                       return r;
                   };

   auto proc_model_v2 = [proc_node] (tinygltf::Model &model) {

                           std::vector<extracted_buffer_info_t> r;

                           unsigned int n_accessors = model.accessors.size();
                           std::cout << "debug:: model has " << n_accessors << " accessors" << std::endl;

                           // --- Materials ---

                           std::cout << "INFO:: this model contains " << model.materials.size() << " materials" << std::endl;
                           for (unsigned int imat=0; imat<model.materials.size(); imat++) {
                              const tinygltf::Material &material = model.materials[imat];
                              std::cout << "Material " << imat << " name: " << material.name << " alphaMode" << material.alphaMode
                                        << " pbr-metallicroughness colour "
                                        << material.pbrMetallicRoughness.baseColorFactor[0] << " "
                                        << material.pbrMetallicRoughness.baseColorFactor[1] << " "
                                        << material.pbrMetallicRoughness.baseColorFactor[2] << " "
                                        << material.pbrMetallicRoughness.baseColorFactor[3] << " "
                                        << " metalicFactor " << material.pbrMetallicRoughness.metallicFactor
                                        << " roughnessFactor " << material.pbrMetallicRoughness.roughnessFactor
                                        << std::endl;
                           }

                           // --- Vertices and Indices ---
#if 0
                           // scenes

                           for (unsigned int iscene=0; iscene<model.scenes.size(); iscene++) {
                              const tinygltf::Scene &scene = model.scenes[iscene];
                              std::cout << ":::: This scene contains " << scene.nodes.size() << " nodes "<< std::endl;
                              for (size_t i = 0; i < scene.nodes.size(); i++) {
                                 auto nodes = proc_node(model, i, model.nodes[scene.nodes[i]]);
                                 r.insert(r.end(), nodes.begin(), nodes.end());
                              }
                           }
#endif

                           // model

                           std::cout << "This model contains " << model.nodes.size() << " nodes" << std::endl;

                           for (size_t i = 0; i < model.nodes.size(); i++) {
                              auto nodes = proc_node(model, i, model.nodes[i]);
                              r.insert(r.end(), nodes.begin(), nodes.end());
                           }

                           return r;
                        };

   bool status = true;
   tinygltf::Model model;
   tinygltf::TinyGLTF loader;
   std::string err;
   std::string warn;

   std::string file_name = file_name_in;
   if (coot::file_exists(file_name)) {
      // do nothing
      std::cout << "debug:: file " << file_name << " exists and will be  read" << std::endl;
   } else {
      std::string dir = coot::package_data_dir();
      std::string dir_2 = coot::util::append_dir_dir(dir, "glTF");
      file_name = coot::util::append_dir_file(dir_2, file_name);
   }

   bool use_binary = false;
   std::string ext = coot::util::file_name_extension(file_name_in);
   if (ext == ".glb") use_binary = true;

   // use the extension to check which function to use
   //
   // bool ret = loader.LoadASCIIFromFile(&model, &err, &warn, file_name);
   bool read_success = false;

   if (use_binary) {
      read_success = loader.LoadBinaryFromFile(&model, &err, &warn, file_name); // for binary glTF(.glb)
   } else {
      std::cout << ":::::::::::::::::: read ASCII :::::::::::::::::: " << std::endl;
      read_success = loader.LoadASCIIFromFile(&model, &err, &warn, file_name);
      std::cout << ":::::::::::::::::: read ASCII success " << read_success << "  :::::::::::::::::: " << std::endl;
   }

   if (!warn.empty()) {
      std::cout << "WARNING:: load_from_glTF():" << warn << std::endl;
   }

   if (!err.empty()) {
      std::cout << "ERROR:: load_from_glTF(): " << err << std::endl;
   }

   if (read_success == false) {
      std::cout << "WARNING:: failed to parse glTF from " << file_name << std::endl;
   }

   if (err.empty() && read_success) { // OK, good, let's go!

      std::cout << "load_from_glTF() A ::::::::::::: success path, calling proc_model_v2()" << std::endl;
      std::vector<extracted_buffer_info_t> r = proc_model_v2(model);
      std::cout << "load_from_glTF(): B ::::::::::: found " << r.size() << " mesh primitives" << std::endl;
      for (unsigned int i=0; i<r.size(); i++) {
         const extracted_buffer_info_t &ebi = r[i];
         if (ebi.normals.size() > 0) {
	    std::cout << "debug ebi.normals size  " << ebi.normals.size() << " ebi positions size " << ebi.positions.size() << " "
		      << "ebi vertex colours size " << ebi.vertex_colours.size() << std::endl;
            if (ebi.normals.size() == ebi.positions.size()) {
	       if (ebi.normals.size() == ebi.vertex_colours.size()) {
		  unsigned int max_vertex_index = 0; // set this
		  if (max_vertex_index < ebi.normals.size()) {
		     unsigned int idx_vert_base = vertices.size();
		     unsigned int idx_tri_base = triangles.size();
		     for (unsigned int j=0; j<ebi.normals.size(); j++) {
			s_generic_vertex g(ebi.positions[j], ebi.normals[j], ebi.vertex_colours[j]);
			// std::cout << "s_generic_vertex colour " << glm::to_string(g.color) << std::endl;
			vertices.push_back(g);
		     }
		     triangles.insert(triangles.end(), ebi.triangles.begin(), ebi.triangles.end());
		     if (idx_vert_base != 0) {
			for (unsigned int ii=idx_tri_base; ii<triangles.size(); ii++) {
			   triangles[ii].rebase(idx_vert_base);
			}
		     }
		  }
	       } else {
		  unsigned int max_vertex_index = 0; // set this
		  if (max_vertex_index < ebi.normals.size()) {
		     unsigned int idx_vert_base = vertices.size();
		     unsigned int idx_tri_base = triangles.size();
		     for (unsigned int j=0; j<ebi.normals.size(); j++) {
			s_generic_vertex g(ebi.positions[j], ebi.normals[j], ebi.base_colour);
			vertices.push_back(g);
		     }
		     triangles.insert(triangles.end(), ebi.triangles.begin(), ebi.triangles.end());
		     if (idx_vert_base != 0) {
			for (unsigned int ii=idx_tri_base; ii<triangles.size(); ii++) {
			   triangles[ii].rebase(idx_vert_base);
			}
		     }
		  }
	       }
            }
         }
      }

      if (include_call_to_setup_buffers) {
         std::cout << "pre-setup_buffer() " << vertices.size() << " vertices "
		   << triangles.size() << " triangles "  << std::endl;
         setup_buffers();
      }
   } else {
      std::cout << "::::::::::::: non-success path" << std::endl;
      status = false; // boo
   }

   std::cout << "load_from_glTF() returns status " << status << std::endl;
   return status;
}


// pass the name (that should be visible in blender?)
// pass one of "SHINY_PLASTIC, CLAY or METALLIC" as a hint for how to consstruct the glTF material.
void
Mesh::export_to_glTF(const std::string &file_name, bool use_binary_format) const {

   // make a tinygltf::Model and use that for export

   // perhaps this can go in the header - it's the same as in the reader
   class extracted_buffer_info_t {
   public:
      int buffer_view_index;
      std::vector<glm::vec3> positions;
      std::vector<glm::vec3> normals;
      std::vector<glm::vec4> vertex_colours;
      std::vector<glm::vec2> texture_coords;
      std::vector<g_triangle> triangles;
      glm::vec4 base_colour;
   };


   auto fill_buffer_info = [] (extracted_buffer_info_t &ebi,  // reference
                               const std::vector<s_generic_vertex> &vertices,
                               const std::vector<g_triangle> &triangles) {

                              ebi.positions.resize(vertices.size());
                              ebi.normals.resize(vertices.size());
                              ebi.vertex_colours.resize(vertices.size());
                              for (unsigned int i=0; i<vertices.size(); i++)
                                 ebi.positions[i] = vertices[i].pos;
                              for (unsigned int i=0; i<vertices.size(); i++)
                                 ebi.normals[i] = vertices[i].normal;
                              for (unsigned int i=0; i<vertices.size(); i++)
                                 ebi.vertex_colours[i] = vertices[i].color;
                              ebi.triangles = triangles;
                           };

   extracted_buffer_info_t ebi;
   fill_buffer_info(ebi, vertices, triangles); // fills ebi

   tinygltf::Model model;
   model.asset.generator = "Coot 1.0-pre 20220205";
   model.defaultScene = 0;

   // --- Materials ---

   struct tinygltf::Material mat = tinygltf::Material();

   // and/or set the data items:

   mat.name = name;
   mat.emissiveFactor = std::vector<double>(3,0.0);
   mat.alphaMode = "OPAQUE";
   mat.alphaCutoff = 0.5;
   mat.doubleSided = true;
   mat.pbrMetallicRoughness = tinygltf::PbrMetallicRoughness();
   mat.normalTexture        = tinygltf::NormalTextureInfo();
   mat.occlusionTexture     = tinygltf::OcclusionTextureInfo();
   mat.emissiveTexture      = tinygltf::TextureInfo();
   mat.pbrMetallicRoughness.metallicFactor = 0.0;  // shiny plastic
   mat.pbrMetallicRoughness.roughnessFactor = 0.2;

   model.materials.push_back(mat); // index 0
   int material_index = 0;

   // --- Vertices and Indices ---

   tinygltf::Scene scene;
   scene.name = "A Molecular Scene from Coot";
   tinygltf::Node node;

   int mesh_index = 0;
   node.name = name;
   node.mesh = mesh_index;

   tinygltf::Mesh mesh;
   mesh.name = name;

   int indices_index = 0; // index of the access for the indices view

   tinygltf::Buffer start_buffer;
   model.buffers.push_back(start_buffer); // position in vector is the buffer_index
   tinygltf::Buffer &buffer = model.buffers.back();
   buffer.name = name;
   buffer.uri = "";
   int buffer_index = 0;

   std::cout << "debug:: buffer.data.size() is now " << buffer.data.size() << std::endl;

   {  // ---- positions ----

      tinygltf::Accessor accessor_for_positions;
      int buffer_view_index = 0;
      accessor_for_positions.bufferView = buffer_view_index;
      accessor_for_positions.name = name;
      accessor_for_positions.byteOffset = 0;
      accessor_for_positions.normalized = false;
      accessor_for_positions.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
      size_t n_items = ebi.positions.size();
      size_t n_bytes = n_items * sizeof(glm::vec3);
      accessor_for_positions.count = n_items; // not bytes, probably
      accessor_for_positions.type = TINYGLTF_TYPE_VEC3;

      tinygltf::BufferView buffer_view;
      buffer_view.name = name;
      buffer_view.buffer = buffer_index;
      buffer_view.byteOffset = 0; // this is the first addition to buffer data
      buffer_view.byteLength = n_bytes;
      // buffer_view.byteStride;

      // let's make a unsigned char vector for these postions and append them to
      // buffer

      std::vector<glm::vec3> *vec3_data = &ebi.positions;
      std::vector<unsigned char> *uchar_data = reinterpret_cast<std::vector<unsigned char> *> (vec3_data);
      buffer.data.insert(buffer.data.end(), uchar_data->begin(), uchar_data->end());

      buffer_view.target = TINYGLTF_TARGET_ARRAY_BUFFER;

      model.bufferViews.push_back(buffer_view);

      // min and max values for postions
      double big = 9e9;
      std::vector<double> min_values = { big,  big,  big};
      std::vector<double> max_values = {-big, -big, -big};

      for (const auto &v : ebi.positions) {
	 for(unsigned int i=0; i<3; i++) {
	    if (v[i] > max_values[i]) max_values[i] = v[i];
	    if (v[i] < min_values[i]) min_values[i] = v[i];
	 }
      }
      accessor_for_positions.minValues = min_values;
      accessor_for_positions.maxValues = max_values;
      model.accessors.push_back(accessor_for_positions);

   }
   std::cout << "debug:: buffer.data.size() is now " << buffer.data.size() << std::endl;

   {  // ---- normals ----

      tinygltf::Accessor accessor_for_normals;
      int buffer_view_index = 1;
      accessor_for_normals.bufferView = buffer_view_index;
      accessor_for_normals.name = name;
      accessor_for_normals.byteOffset = 0;
      accessor_for_normals.normalized = false;
      accessor_for_normals.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
      size_t n_items = ebi.positions.size();
      size_t n_bytes = n_items * sizeof(glm::vec3);
      accessor_for_normals.count = n_items;
      accessor_for_normals.type = TINYGLTF_TYPE_VEC3;

      // now add to buffer.data (std::vector<unsigned char>)

      tinygltf::BufferView buffer_view;
      buffer_view.name = name;
      buffer_view.buffer = buffer_index;
      buffer_view.byteOffset = buffer.data.size();
      buffer_view.byteLength = n_bytes;
      buffer_view.byteStride = 0; // tightly packed.
      buffer_view.target = TINYGLTF_TARGET_ARRAY_BUFFER;

      model.bufferViews.push_back(buffer_view);
      std::vector<glm::vec3> *vec3_data = &ebi.normals;
      std::vector<unsigned char> *uchar_data = reinterpret_cast<std::vector<unsigned char> *> (vec3_data);
      buffer.data.insert(buffer.data.end(), uchar_data->begin(), uchar_data->end());
      model.accessors.push_back(accessor_for_normals);
   }
   std::cout << "debug:: buffer.data.size() is now " << buffer.data.size() << std::endl;

   {  // ---- colours ----

      tinygltf::Accessor accessor_for_colours;
      int buffer_view_index = 2;
      accessor_for_colours.bufferView = buffer_view_index;
      accessor_for_colours.name = name;
      accessor_for_colours.byteOffset = 0;
      accessor_for_colours.normalized = false;
      accessor_for_colours.componentType = TINYGLTF_COMPONENT_TYPE_FLOAT;
      size_t n_items = ebi.positions.size();
      size_t n_bytes = n_items * sizeof(glm::vec4);
      accessor_for_colours.count = n_items;
      accessor_for_colours.type = TINYGLTF_TYPE_VEC4;

      // now add to buffer.data (std::vector<unsigned char>)

      tinygltf::BufferView buffer_view;
      buffer_view.name = name;
      buffer_view.buffer = buffer_index;
      buffer_view.byteOffset = buffer.data.size();
      buffer_view.byteLength = n_bytes;
      buffer_view.byteStride = 0; // tightly packed.
      buffer_view.target = TINYGLTF_TARGET_ARRAY_BUFFER;

      model.bufferViews.push_back(buffer_view);
      std::vector<glm::vec4> *vec4_data = &ebi.vertex_colours;
      std::vector<unsigned char> *uchar_data = reinterpret_cast<std::vector<unsigned char> *> (vec4_data);
      buffer.data.insert(buffer.data.end(), uchar_data->begin(), uchar_data->end());
      model.accessors.push_back(accessor_for_colours);
   }
   std::cout << "debug:: buffer.data.size() is now " << buffer.data.size() << std::endl;

   int accessor_index_for_indices = -1;
   { // ---------- indices/triangles --------------

      tinygltf::Accessor accessor_for_indices;
      int buffer_view_index = 3; // index of the buffer
      accessor_for_indices.bufferView = buffer_view_index;
      accessor_for_indices.name = name;
      accessor_for_indices.byteOffset = 0;
      accessor_for_indices.normalized = false;  // not that this means much.
      accessor_for_indices.componentType = TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT;
      size_t n_bytes = ebi.triangles.size() * 3 * sizeof(unsigned int);
      accessor_for_indices.count = ebi.triangles.size() * 3;
      accessor_for_indices.type = TINYGLTF_TYPE_SCALAR;

      tinygltf::BufferView buffer_view;
      buffer_view.name = name;
      buffer_view.buffer = buffer_index;
      buffer_view.byteOffset = buffer.data.size();
      buffer_view.byteLength = n_bytes;
      buffer_view.byteStride = 0; // tightly packed
      buffer_view.target = TINYGLTF_TARGET_ELEMENT_ARRAY_BUFFER;
      std::vector<g_triangle> *index_data = &ebi.triangles;
      std::vector<unsigned char> *uchar_data = reinterpret_cast<std::vector<unsigned char> *> (index_data);
      buffer.data.insert(buffer.data.end(), uchar_data->begin(), uchar_data->end());

      std::cout << "DEBUG:: buffer_view for triangles triangles.size " << triangles.size() << std::endl;
      std::cout << "DEBUG:: buffer_view for triangles: n_bytes " << n_bytes << std::endl;
      std::cout << "DEBUG:: buffer_view for triangles: accessor_for_indices.count "
		<< accessor_for_indices.count << std::endl;

      model.accessors.push_back(accessor_for_indices);
      accessor_index_for_indices = model.accessors.size() - 1;
      model.bufferViews.push_back(buffer_view);
   }

   std::cout << "debug:: buffer.data.size() is now " << buffer.data.size() << std::endl;

   tinygltf::Primitive prim;
   int position_idx = 0;
   int normal_idx = 1;
   int colour_idx = 2;
   prim.attributes["POSITION"] = position_idx;
   prim.attributes["NORMAL"]   =   normal_idx;
   prim.attributes["COLOR_0"]  =   colour_idx;
   prim.material = material_index;
   prim.indices = indices_index;
   prim.mode = TINYGLTF_MODE_TRIANGLES; // this will need some thought for a map as lines.
   prim.indices = accessor_index_for_indices;

   mesh.primitives.push_back(prim);

   model.meshes.push_back(mesh);
   model.nodes.push_back(node); // index 0
   int node_index = 0;
   scene.nodes.push_back(node_index);

   model.scenes.push_back(scene);

   // Now write the beasty!

   bool embedImages  = true;
   bool embedBuffers = true;
   bool prettyPrint  = true;
   tinygltf::TinyGLTF tgltf;
   std::cout << "--------- calling WriteGltfSceneToFile() with uri \"" << model.buffers[0].uri << "\"" << std::endl;
   std::cout << "--------- calling WriteGltfSceneToFile() with use_binary_format " << use_binary_format << std::endl;
   tgltf.WriteGltfSceneToFile(&model, file_name, embedImages, embedBuffers, prettyPrint, use_binary_format);

}


