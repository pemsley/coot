
#include "coot-utils/tiny_gltf.h"

#include "molecules-container.hh"
#include "coot-utils/oct.hh"
#include "coot-utils/vertex.hh"

//! @params `n_divisions` is a number divisble by 2, at least 4 (typically 16)
//! @return a unit-vector end-cap octohemisphere mesh
coot::simple_mesh_t
molecules_container_t::get_octahemisphere(unsigned int n_slices) const {

   coot::simple_mesh_t m;
   unsigned int num_subdivisions = 3;
   if (n_slices ==  4) num_subdivisions = 1;
   if (n_slices ==  8) num_subdivisions = 2;
   if (n_slices == 16) num_subdivisions = 3;
   if (n_slices == 32) num_subdivisions = 4;

   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > h = tessellate_hemisphere_patch(num_subdivisions);
   const auto &hemisphere_vertices = h.first;
   std::vector<coot::api::vnc_vertex> vertices(hemisphere_vertices.size());
   const std::vector<g_triangle> &triangles = h.second;
   glm::vec4 c(0.5f, 0.5f, 0.5f, 1.0f);
   for (unsigned int i=0; i<hemisphere_vertices.size(); i++) {
      const auto &p = hemisphere_vertices[i];
      const auto &n = hemisphere_vertices[i];
      vertices[i] = coot::api::vnc_vertex(p, n, c);
   }
   m = coot::simple_mesh_t(vertices, triangles);
   return m;
}


coot::simple_mesh_t
molecules_container_t::make_mesh_from_gltf_file(const std::string &file_name_in) {

   coot::simple_mesh_t sm;

   auto &triangles = sm.triangles;
   auto &vertices = sm.vertices;

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

      std::vector<extracted_buffer_info_t> r = proc_model_v2(model);
      for (unsigned int i=0; i<r.size(); i++) {
         const extracted_buffer_info_t &ebi = r[i];
         if (ebi.normals.size() > 0) {
            if (false)
               std::cout << "debug ebi normals size  " << ebi.normals.size()
                         << " ebi positions size " << ebi.positions.size() << " "
                         << "ebi vertex colours size " << ebi.vertex_colours.size() << std::endl;
            if (ebi.normals.size() == ebi.positions.size()) {
	       if (ebi.normals.size() == ebi.vertex_colours.size()) {
		  unsigned int max_vertex_index = 0; // set this
		  if (max_vertex_index < ebi.normals.size()) {
		     unsigned int idx_vert_base = vertices.size();
		     unsigned int idx_tri_base = triangles.size();
		     for (unsigned int j=0; j<ebi.normals.size(); j++) {
                        coot::api::vnc_vertex g(ebi.positions[j], ebi.normals[j], ebi.vertex_colours[j]);
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
                        coot::api::vnc_vertex g(ebi.positions[j], ebi.normals[j], ebi.base_colour);
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

   } else {
      std::cout << "::::::::::::: non-success path" << std::endl;
      status = false; // boo
   }

   std::cout << "load_from_glTF() status " << status << std::endl;

   return sm;
}
