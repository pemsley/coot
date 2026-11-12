/*
 * coot-utils/textured-mesh.cc
 *
 * Copyright 2026 by Medical Research Council
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

#include <cstring>
#include <iostream>

#include "tiny_gltf.h"
#include "utils/coot-utils.hh"
#include "textured-mesh.hh"

namespace {

   // The elements of an accessor as floats, n_components per element, or an empty vector
   // if the accessor can't be read. Integer components are accepted only if they are
   // normalized (as glTF requires of integer colours and texture coordinates).
   std::vector<float>
   accessor_floats(const tinygltf::Model &model, int accessor_index, int *n_components_p) {

      std::vector<float> v;
      *n_components_p = 0;
      if (accessor_index < 0 || accessor_index >= static_cast<int>(model.accessors.size())) return v;
      const tinygltf::Accessor &accessor = model.accessors[accessor_index];
      if (accessor.bufferView < 0) return v;
      const tinygltf::BufferView &buffer_view = model.bufferViews[accessor.bufferView];
      const tinygltf::Buffer &buffer = model.buffers[buffer_view.buffer];
      int n_comp = tinygltf::GetNumComponentsInType(accessor.type);
      int comp_size = tinygltf::GetComponentSizeInBytes(accessor.componentType);
      int stride = accessor.ByteStride(buffer_view);
      if (n_comp <= 0 || comp_size <= 0 || stride <= 0) return v;
      if (accessor.componentType != TINYGLTF_COMPONENT_TYPE_FLOAT && ! accessor.normalized) return v;
      size_t start = buffer_view.byteOffset + accessor.byteOffset;
      if (accessor.count > 0)
         if (start + (accessor.count - 1) * stride + n_comp * comp_size > buffer.data.size()) return v;

      v.resize(accessor.count * n_comp);
      for (size_t i=0; i<accessor.count; i++) {
         const unsigned char *p = buffer.data.data() + start + i * stride;
         for (int c=0; c<n_comp; c++) {
            float f = 0.0f;
            if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_FLOAT) {
               std::memcpy(&f, p + 4 * c, 4);
            } else if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE) {
               f = static_cast<float>(p[c]) / 255.0f;
            } else if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT) {
               uint16_t s;
               std::memcpy(&s, p + 2 * c, 2);
               f = static_cast<float>(s) / 65535.0f;
            } else {
               return std::vector<float>();
            }
            v[i * n_comp + c] = f;
         }
      }
      *n_components_p = n_comp;
      return v;
   }

   std::vector<unsigned int>
   accessor_indices(const tinygltf::Model &model, int accessor_index) {

      std::vector<unsigned int> v;
      if (accessor_index < 0 || accessor_index >= static_cast<int>(model.accessors.size())) return v;
      const tinygltf::Accessor &accessor = model.accessors[accessor_index];
      if (accessor.bufferView < 0) return v;
      const tinygltf::BufferView &buffer_view = model.bufferViews[accessor.bufferView];
      const tinygltf::Buffer &buffer = model.buffers[buffer_view.buffer];
      int comp_size = tinygltf::GetComponentSizeInBytes(accessor.componentType);
      if (comp_size <= 0) return v;
      size_t start = buffer_view.byteOffset + accessor.byteOffset;
      if (start + accessor.count * comp_size > buffer.data.size()) return v;

      v.resize(accessor.count);
      const unsigned char *p = buffer.data.data() + start;
      for (size_t i=0; i<accessor.count; i++) {
         if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_BYTE) {
            v[i] = p[i];
         } else if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT) {
            uint16_t s;
            std::memcpy(&s, p + 2 * i, 2);
            v[i] = s;
         } else if (accessor.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT) {
            uint32_t u;
            std::memcpy(&u, p + 4 * i, 4);
            v[i] = u;
         } else {
            return std::vector<unsigned int>();
         }
      }
      return v;
   }

   // tinygltf has already decoded the image, expanded to 4 channels by default.
   // Return an empty texture if it is in a form we don't handle.
   coot::texture_rgba_t
   image_to_texture(const tinygltf::Image &image) {

      coot::texture_rgba_t t;
      if (image.component != 4) return t;
      size_t n_pixels = static_cast<size_t>(image.width) * static_cast<size_t>(image.height);
      if (image.bits == 8) {
         if (image.image.size() != n_pixels * 4) return t;
         t.data = image.image;
      } else if (image.bits == 16) {
         if (image.image.size() != n_pixels * 8) return t;
         t.data.resize(n_pixels * 4);
         for (size_t i=0; i<n_pixels * 4; i++) {
            uint16_t s;
            std::memcpy(&s, image.image.data() + 2 * i, 2);
            t.data[i] = static_cast<unsigned char>(s >> 8);
         }
      } else {
         return t;
      }
      t.width  = image.width;
      t.height = image.height;
      return t;
   }
}


bool
coot::textured_mesh_t::read_gltf(const std::string &file_name_in) {

   std::string file_name = file_name_in;
   if (! coot::file_exists(file_name)) {
      std::string dir = coot::util::append_dir_dir(coot::package_data_dir(), "glTF");
      file_name = coot::util::append_dir_file(dir, file_name_in);
   }

   tinygltf::Model model;
   tinygltf::TinyGLTF loader;
   std::string err;
   std::string warn;
   bool read_success = false;
   if (coot::util::file_name_extension(file_name) == ".glb")
      read_success = loader.LoadBinaryFromFile(&model, &err, &warn, file_name);
   else
      read_success = loader.LoadASCIIFromFile(&model, &err, &warn, file_name);

   if (! warn.empty())
      std::cout << "WARNING:: textured_mesh_t::read_gltf(): " << warn << std::endl;
   if (! err.empty())
      std::cout << "ERROR:: textured_mesh_t::read_gltf(): " << err << std::endl;
   if (! read_success) {
      std::cout << "WARNING:: textured_mesh_t::read_gltf(): failed to read " << file_name << std::endl;
      return false;
   }

   int texture_index_used = -1;

   for (const auto &gltf_mesh : model.meshes) {
      for (const auto &primitive : gltf_mesh.primitives) {

         if (primitive.mode != TINYGLTF_MODE_TRIANGLES && primitive.mode != -1) continue;

         auto it_pos  = primitive.attributes.find("POSITION");
         auto it_norm = primitive.attributes.find("NORMAL");
         if (it_pos == primitive.attributes.end()) continue;
         if (it_norm == primitive.attributes.end()) {
            std::cout << "WARNING:: textured_mesh_t::read_gltf(): skipping a primitive without normals in "
                      << file_name << std::endl;
            continue;
         }

         int n_comp_pos = 0;
         int n_comp_norm = 0;
         std::vector<float> positions = accessor_floats(model, it_pos->second,  &n_comp_pos);
         std::vector<float> normals   = accessor_floats(model, it_norm->second, &n_comp_norm);
         if (n_comp_pos != 3 || n_comp_norm != 3) continue;
         if (positions.size() != normals.size()) continue;
         unsigned int n_vertices = positions.size() / 3;

         // the material: the base colour (which tints the texture) and the texture

         glm::vec4 base_colour(1.0f, 1.0f, 1.0f, 1.0f);
         std::string texcoord_name = "TEXCOORD_0";
         if (primitive.material >= 0 && primitive.material < static_cast<int>(model.materials.size())) {
            const tinygltf::PbrMetallicRoughness &pbr = model.materials[primitive.material].pbrMetallicRoughness;
            if (pbr.baseColorFactor.size() == 4)
               base_colour = glm::vec4(pbr.baseColorFactor[0], pbr.baseColorFactor[1],
                                       pbr.baseColorFactor[2], pbr.baseColorFactor[3]);
            int texture_index = pbr.baseColorTexture.index;
            if (texture_index >= 0 && texture_index < static_cast<int>(model.textures.size())) {
               texcoord_name = "TEXCOORD_" + std::to_string(pbr.baseColorTexture.texCoord);
               if (texture_index_used == -1) {
                  int image_index = model.textures[texture_index].source;
                  if (image_index >= 0 && image_index < static_cast<int>(model.images.size())) {
                     texture = image_to_texture(model.images[image_index]);
                     if (texture.empty())
                        std::cout << "WARNING:: textured_mesh_t::read_gltf(): unhandled image format in "
                                  << file_name << std::endl;
                     else
                        texture_index_used = texture_index;
                  }
               } else {
                  if (texture_index != texture_index_used)
                     std::cout << "WARNING:: textured_mesh_t::read_gltf(): only one texture is used; ignoring texture "
                               << texture_index << " in " << file_name << std::endl;
               }
            }
         }

         // optional per-vertex attributes

         std::vector<float> colours;
         int n_comp_col = 0;
         auto it_col = primitive.attributes.find("COLOR_0");
         if (it_col != primitive.attributes.end()) {
            colours = accessor_floats(model, it_col->second, &n_comp_col);
            if (colours.size() != n_vertices * n_comp_col) colours.clear();
         }

         std::vector<float> texcoords;
         int n_comp_uv = 0;
         auto it_uv = primitive.attributes.find(texcoord_name);
         if (it_uv != primitive.attributes.end()) {
            texcoords = accessor_floats(model, it_uv->second, &n_comp_uv);
            if (n_comp_uv != 2 || texcoords.size() != n_vertices * 2) texcoords.clear();
         }

         std::vector<unsigned int> indices = accessor_indices(model, primitive.indices);
         if (primitive.indices < 0) {
            indices.resize(n_vertices); // non-indexed: consecutive triples
            for (unsigned int i=0; i<n_vertices; i++) indices[i] = i;
         }
         if (indices.size() % 3 != 0) continue;

         unsigned int idx_base = mesh.vertices.size();
         for (unsigned int i=0; i<n_vertices; i++) {
            glm::vec3 pos(positions[3*i], positions[3*i+1], positions[3*i+2]);
            glm::vec3 normal(normals[3*i], normals[3*i+1], normals[3*i+2]);
            glm::vec4 col = base_colour;
            if (! colours.empty()) {
               col.r *= colours[n_comp_col*i];
               col.g *= colours[n_comp_col*i+1];
               col.b *= colours[n_comp_col*i+2];
               if (n_comp_col == 4) col.a *= colours[n_comp_col*i+3];
            }
            mesh.vertices.push_back(api::vnc_vertex(pos, normal, col));
            if (texcoords.empty())
               uvs.push_back(glm::vec2(0.0f, 0.0f));
            else
               uvs.push_back(glm::vec2(texcoords[2*i], texcoords[2*i+1]));
         }
         for (unsigned int i=0; i<indices.size(); i+=3) {
            if (indices[i] >= n_vertices || indices[i+1] >= n_vertices || indices[i+2] >= n_vertices) continue;
            mesh.triangles.push_back(g_triangle(idx_base + indices[i], idx_base + indices[i+1], idx_base + indices[i+2]));
         }
      }
   }

   if (mesh.name.empty())
      mesh.name = coot::util::file_name_non_directory(file_name_in);
   return ! mesh.vertices.empty();
}
