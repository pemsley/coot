/*
 * coot-utils/ortep.hh
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

#include <vector>
#include <glm/glm.hpp>
#include "g_triangle.hh"

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

// put this inside ortep_t?
ortep_t
tessellate_sphere_sans_octant();

