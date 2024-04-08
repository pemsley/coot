/*
 * api/blender-mesh.hh
 * 
 * Copyright 2020 by Medical Research Council
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
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifndef COOT_API_BLENDER_MESH_HH
#define COOT_API_BLENDER_MESH_HH

#include "instancing.hh"

namespace coot {

   class blender_triangle_t {
   public:
      blender_triangle_t() : colour_index(-1) {};
      blender_triangle_t(const g_triangle &t, int col_idx) : triangle(t), colour_index(col_idx) {}
      g_triangle triangle;
      int colour_index;

      void rebase(const unsigned int &idx_base) {
         triangle.rebase(idx_base);
      }
   };

   class blender_mesh_t {
      static size_t make_colour_hash(const glm::vec4 &col);
   public:
      blender_mesh_t() {}
      explicit blender_mesh_t(const instanced_mesh_t &im);
      explicit blender_mesh_t(const simple_mesh_t &sm);
      std::map<int, glm::vec4> colour_table; // include alpha
      std::vector<glm::vec3> vertices;
      std::vector<glm::vec3> normals;
      std::vector<blender_triangle_t> triangles;
   };

   // the array/lists that we need to send to Blender Python are:
   // PyList_SetItem(r_py, 0, vertices_py);                  // get_vertices_for_blender()
   // PyList_SetItem(r_py, 1, tris_py);                      // get_triangles_for_blender()
   // PyList_SetItem(r_py, 2, face_colours_dict_py);         // get_colour_table_for_blender()

   // where
   // vertices_py:          a list of (x,y,z)
   // tris_py:              a 4-member (int) list of vertex indices for a triangle and a colour_index (int)
   // face_colours_dict_py: a list of index and 3 members (idx, (r,g,b))
}


#endif // BLENDER_MESH_HH
