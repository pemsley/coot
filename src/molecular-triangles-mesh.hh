/*
 * src/molecular-triangles-mesh.hh
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
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */
#ifndef MOLECULAR_TRIANGLES_MESH_HH
#define MOLECULAR_TRIANGLES_MESH_HH

#include <vector>
#include <string>
#include "generic-vertex.hh"
#include "coot-utils/g_triangle.hh"

class molecular_triangles_mesh_t {
public:
   molecular_triangles_mesh_t() {}
   molecular_triangles_mesh_t(const std::vector<s_generic_vertex> &vertices_in,
                              const std::vector<g_triangle> &triangles_in,
                              const std::string &name_in) :
      vertices(vertices_in), triangles(triangles_in), name(name_in) { type_index = 0; }
   std::vector<s_generic_vertex> vertices;
   std::vector<g_triangle> triangles;
   std::string name;
   unsigned int type_index;
   void add_to_mesh(const std::vector<s_generic_vertex> &gv,
                    const std::vector<g_triangle> &tris);
   void add_to_mesh(const molecular_triangles_mesh_t &new_mesh);
   bool empty() const { return vertices.empty(); }

};


#endif // MOLECULAR_TRIANGLES_MESH_HH
