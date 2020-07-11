
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software (oct.cc) and associated
// documentation files (the "Software"), to deal in the Software
// without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// # This script can generate spheres, rounded cubes, and capsules.
// # For more information, see https://prideout.net/blog/octasphere/
// # Copyright 2019 Philip Rideout
// # Translated into C++/glm by Paul Emsley

#include <iostream>
#include <cmath>
#include <vector>
#include <glm/ext.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/string_cast.hpp>

#include "g_triangle.hh"
#include "oct.hh"

std::vector<glm::vec3>
compute_geodesic(const glm::vec3 &point_a, const glm::vec3 &point_b,
                 unsigned int num_segments) {

   float angle_between_endpoints = std::acos(glm::dot(point_a, point_b));
   glm::vec3 rotation_axis = glm::cross(point_a, point_b);
   std::vector<glm::vec3> point_list;
   point_list.push_back(point_a);
   if (num_segments == 0)
      return point_list;
   float dtheta = angle_between_endpoints / static_cast<float>(num_segments);
   for (unsigned int point_index=1; point_index<num_segments; point_index++) {
      float theta = static_cast<float>(point_index) * dtheta;
      glm::vec3 rotated_point = glm::rotate(point_a, theta, rotation_axis);
      point_list.push_back(rotated_point);
   }
   point_list.push_back(point_b);
   return point_list;
}

std::pair<std::vector<glm::vec3>, std::vector<g_triangle> >
tessellate_octasphere_patch(unsigned int num_subdivisions) {
   unsigned int n = std::pow(2, num_subdivisions) + 1;
   unsigned int num_verts = n * (n + 1); // 2
   std::vector<glm::vec3> verts;
   std::vector<g_triangle> triangles;
   for (unsigned int i=0; i<n; i++) {
      float theta = M_PI * 0.5 * static_cast<float>(i) / static_cast<float>(n - 1);
      glm::vec3 point_a(0, std::sin(theta), std::cos(theta));
      glm::vec3 point_b(std::cos(theta), std::sin(theta), 0);
      unsigned int num_segments = n - 1 - i;
      std::vector<glm::vec3> geodesic_verts = compute_geodesic(point_a, point_b, num_segments);
      std::vector<g_triangle> geo_triangles = make_octasphere_triangles(i, geodesic_verts.size(), verts.size());
      triangles.insert(triangles.end(), geo_triangles.begin(), geo_triangles.end());
      verts.insert(verts.end(), geodesic_verts.begin(), geodesic_verts.end());
   }

   return std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > (verts,triangles);
}
