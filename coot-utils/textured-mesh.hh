/*
 * coot-utils/textured-mesh.hh
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
#ifndef TEXTURED_MESH_HH
#define TEXTURED_MESH_HH

#include <vector>
#include <string>
#include <glm/glm.hpp>
#include "simple-mesh.hh"

namespace coot {

   //! An 8-bit RGBA image, used as the base-colour texture of a `textured_mesh_t`
   class texture_rgba_t {
   public:
      //! width in pixels
      int width;
      //! height in pixels
      int height;
      //! width * height * 4 bytes, R, G, B, A, rows from the top of the image down -
      //! so row 0 is at texture coordinate v = 0 (the glTF convention)
      std::vector<unsigned char> data;
      texture_rgba_t() : width(0), height(0) {}
      //! @return true if there is no image
      bool empty() const { return data.empty(); }
   };

   //! A `simple_mesh_t` with texture coordinates and a base-colour texture
   //!
   //! The vertex colours of `mesh` tint the texture (the texture colour is multiplied by them),
   //! as COLOR_0 and baseColorFactor do in glTF. So they are white for an untinted texture.
   class textured_mesh_t {
   public:
      //! the geometry: positions, normals, (tint) colours and triangles. The name is `mesh.name`.
      simple_mesh_t mesh;
      //! texture coordinates, one per vertex of `mesh`
      std::vector<glm::vec2> uvs;
      //! the base-colour texture (empty if the mesh has none)
      texture_rgba_t texture;

      textured_mesh_t() {}

      //! @return true if there is a texture and every vertex has texture coordinates
      bool has_texture() const { return ! texture.empty() && uvs.size() == mesh.vertices.size(); }

      //! Read the triangle primitives of the meshes in a glTF (.gltf or .glb) file.
      //!
      //! If `file_name` does not exist, it is looked for in the glTF directory of Coot's package data.
      //!
      //! Node transformations are not applied. Only one texture is used: that of the first primitive
      //! that has a base-colour texture. Primitives without normals are skipped. Vertices without
      //! texture coordinates are given (0,0).
      //!
      //! @param file_name the glTF file name
      //! @return true on success
      bool read_gltf(const std::string &file_name);
   };
}

#endif // TEXTURED_MESH_HH
