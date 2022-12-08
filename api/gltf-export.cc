
#include "simple-mesh.hh"
#include "tiny_gltf.h"

// pass the name (that should be visible in blender?)
// pass one of "SHINY_PLASTIC, CLAY or METALLIC" as a hint for how to consstruct the glTF material.
void
coot::simple_mesh_t::export_to_gltf(const std::string &file_name, bool use_binary_format) const {

   std::string name = "Test Export";

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
                               const std::vector<api::vnc_vertex> &vertices,
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


