/*
 * MoleculesToTriangles/CXXClasses/m2t-mesh.hh
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
#ifndef M2T_MESH_HH
#define M2T_MESH_HH

#include <vector>
#include <algorithm>
#include <glm/glm.hpp>

// A minimal vertex/triangle/mesh trio, native to MoleculesToTriangles.
//
// MoleculesToTriangles is low in the link order, so classes here must not
// depend on coot-utils (coot::api::vnc_vertex, g_triangle, coot::simple_mesh_t).
// Callers higher up (e.g. api/) that are allowed to depend on coot-utils
// should convert one of these into a coot::simple_mesh_t.

namespace coot {
   namespace m2t {

      class mesh_vertex_t {
      public:
         glm::vec3 pos;
         glm::vec3 normal;
         glm::vec4 color;
         mesh_vertex_t() {}
         mesh_vertex_t(const glm::vec3 &pos_in, const glm::vec3 &normal_in, const glm::vec4 &color_in) :
            pos(pos_in), normal(normal_in), color(color_in) {}
      };

      class mesh_triangle_t {
      public:
         unsigned int point_id[3];
         mesh_triangle_t() {}
         mesh_triangle_t(unsigned int a0, unsigned int a1, unsigned int a2) {
            point_id[0] = a0;
            point_id[1] = a1;
            point_id[2] = a2;
         }
         void rebase(unsigned int idx_base) {
            point_id[0] += idx_base;
            point_id[1] += idx_base;
            point_id[2] += idx_base;
         }
         void reverse_winding() {
            std::swap(point_id[0], point_id[1]);
         }
      };

      class simple_mesh_t {
      public:
         std::vector<mesh_vertex_t> vertices;
         std::vector<mesh_triangle_t> triangles;
         simple_mesh_t() {}
         simple_mesh_t(const std::vector<mesh_vertex_t> &vertices_in,
                       const std::vector<mesh_triangle_t> &triangles_in) :
            vertices(vertices_in), triangles(triangles_in) {}
         void add_submesh(const simple_mesh_t &submesh) {
            unsigned int idx_base = vertices.size();
            unsigned int idx_base_tri = triangles.size();
            vertices.insert(vertices.end(), submesh.vertices.begin(), submesh.vertices.end());
            triangles.insert(triangles.end(), submesh.triangles.begin(), submesh.triangles.end());
            for (unsigned int i=idx_base_tri; i<triangles.size(); i++)
               triangles[i].rebase(idx_base);
         }
      };

   }
}

#endif // M2T_MESH_HH
