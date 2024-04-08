/*
 * coot-utils/oct.hh
 *
 * Copyright 2023 by Medical Research Council
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

#ifndef OCT_HH
#define OCT_HH

#include <map>
#include <set>
#include <vector>

#include "vertex.hh"
#include "coot-utils/g_triangle.hh"

// should be in its own file?
class ortep_t {
public:
   std::vector<glm::vec3> vertices;
   std::vector<glm::vec3> normals;
   std::vector<g_triangle> triangles;

   std::vector<glm::vec3> vertices_for_lines;
   std::vector<std::pair<unsigned int, unsigned int> > indices_for_lines;

   void make_lines();
   void transform(const glm::mat4 &m);
};

std::vector<g_triangle>
make_octasphere_triangles(unsigned int i_row,
                          unsigned int geodesic_verts_size,
                          unsigned int verts_size);

std::pair<std::vector<glm::vec3>, std::vector<g_triangle> >
tessellate_hemisphere_patch(unsigned int num_subdivisions);

std::pair<std::vector<glm::vec3>, std::vector<g_triangle> >
tessellate_octasphere(unsigned int num_subdivisions, bool remove_redundant_vertices_flag=false);

std::pair<std::vector<coot::api::vnc_vertex>, std::vector<g_triangle> >
make_octasphere(unsigned int num_subdivisions, const glm::vec3 &position,
                float radius, const glm::vec4 &colour_in, bool remove_redundant_vertices=false);

// put this inside ortep_t?
ortep_t
tessellate_sphere_sans_octant();

std::pair<std::vector<coot::api::vnc_vertex>, std::vector<g_triangle> >
make_octasphere_dish(unsigned int num_subdivisions, const glm::vec3 &position,
                     float radius, float radiusAlongNormal,
                     const glm::vec3 &dish_normal,
                     const glm::vec4 &colour_in);

std::map<unsigned int, std::set<unsigned int> >
find_same_vertices(const std::vector<glm::vec3> &verts);



#endif // OCT_HH

