/*
 * src/pumpkin.cc
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

#include <iostream>
#include <fstream>
#include <map>
#define GLM_ENABLE_EXPERIMENTAL
// #include <glm/ext.hpp>
#include "pumpkin.hh"

std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> > pumpkin() {

   std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> > p;
   std::vector<position_normal_vertex> &vertices   = p.first;
   std::vector<g_triangle> &triangles = p.second;

   unsigned int n_slices_raw = 10 * 3; // every 3rd one we add an extra one (starting at 0)
   unsigned int n_slices = n_slices_raw + n_slices_raw/3;

   unsigned int n_for_theta = 109; // for all the way round, which we don't do,
                                   // because we use the inner loop to rotate 360 degrees
   float pi = 3.1415926f;

   // std::ofstream ff("cardioid.table"); // debugging
   unsigned int n_per_longitude_line = 0;
   
   // vertices
   float r_base = 1.0;
   for (unsigned int i=0; i<n_for_theta; i++) {
      float f = static_cast<float>(i)/static_cast<float>(n_for_theta);
      float t = 0.5 * pi * (-1.0 + 4.0 * f); // -pi/2 to  3pi/2
      if (t > 0.7 * 0.5 * pi) continue; // stop looping over the top
      n_per_longitude_line++;
      float alpha = 1.0;
      float R = alpha * (1.0 - sinf(t));
      float x = R * cosf(t);
      float z = R * sinf(t);
      // std::cout << " i " << i << " vertices.size() " << vertices.size() << std::endl;

      for (unsigned int j=0; j<n_slices_raw; j++) {

         bool is_inner_longitude = false;
         // every 3rd sice, decrease the radius
         float r = r_base;
         if (j%3 == 0) {
            is_inner_longitude = true;
            r = 0.88 * r_base;
         }

         // now rotate x and y around the Z axis and make vertices in a ring
         //
         float phi = 2.0 * pi * static_cast<float>(j)/static_cast<float>(n_slices_raw);
         if (is_inner_longitude) {
            float delta_phi = 0.03; // radians of course
            for (int jj=-1; jj<2; jj+=2) {
               float phi_this = phi + delta_phi * static_cast<float>(jj);
               float xx = r * x * cosf(phi_this);
               float yy = r * x * sinf(phi_this);
               float zz = z;
               glm::vec3 pos(xx, yy, zz);
               // ff << pos.x << " " << pos.y << " " << pos.z << "\n";
               position_normal_vertex v;
               v.pos = pos;
               v.normal = glm::vec3(0,0,1);
               vertices.push_back(v);
            }
            
         } else {
            float xx = r * x * cosf(phi);
            float yy = r * x * sinf(phi);
            float zz = z;
            glm::vec3 pos(xx, yy, zz);
            // ff << pos.x << " " << pos.y << " " << pos.z << "\n";
            position_normal_vertex v;
            v.pos = pos;
            v.normal = glm::vec3(0,0,1);
            vertices.push_back(v);
         }
      }
   }

   if (false) {
      std::cout << "---- n_per_longitude_line " << n_per_longitude_line << std::endl;
      std::cout << "------------- n_slices " << n_slices << std::endl;
      std::cout << " final vertices.size() " << vertices.size() << std::endl;
   }

   // triangles
   for (unsigned int i=0; i<n_for_theta; i++) {
      if (i == n_per_longitude_line) break;
      for (unsigned int j=0; j<n_slices; j++) {
         unsigned int idx_base = i * n_slices + j;
         unsigned int idx_a = idx_base;
         unsigned int idx_b = idx_base + 1;
         unsigned int idx_c = idx_base + n_slices;
         unsigned int idx_d = idx_base + n_slices + 1;
         if (j == (n_slices - 1)) {
            idx_b = i * n_slices;
            idx_d = i * n_slices + n_slices;
         }
         if (idx_a < vertices.size() && idx_b < vertices.size() &&
             idx_c < vertices.size() && idx_d < vertices.size()) {
            triangles.push_back(g_triangle(idx_a, idx_b, idx_d));
            triangles.push_back(g_triangle(idx_d, idx_c, idx_a));
         } else {
            if (false)
               std::cout << "pumpkin(): bad index i " << i << " j " << j << " "
                         << idx_a << " " << idx_b << " " << idx_c << " " << idx_d
                         << std::endl;
         }
      }
   }

   // normals of the vertices, average them.
   std::map<unsigned int, std::vector<glm::vec3> > normals_map;
   // collect the face normals, per vertex
   for (unsigned int i=0; i<triangles.size(); i++) {
      const g_triangle &t = triangles[i];
      glm::vec3 sum(0,0,0);
      for (unsigned int j=0; j<3; j++)
         sum += vertices[triangles[i].point_id[j]].pos;
      glm::vec3 n = glm::normalize(glm::cross(vertices[t[1]].pos-vertices[t[0]].pos,
                                              vertices[t[2]].pos-vertices[t[0]].pos));
      for (unsigned int j=0; j<3; j++)
         normals_map[t[j]].push_back(n);
   }
   // average the normals, for each vertex in the normals map
   for (unsigned int i=0; i<vertices.size(); i++) {
      glm::vec3 sum(0,0,0);
      std::map<unsigned int, std::vector<glm::vec3> >::const_iterator it =
         normals_map.find(i);
      if (it != normals_map.end()) {
         const std::vector<glm::vec3> &vv = it->second;
         for (unsigned int j=0; j<vv.size(); j++) {
            const glm::vec3 &v = vv[j];
            sum += v;
         }
      }
      glm::vec3 n = glm::normalize(sum);
      vertices[i].normal = n;
   }
   return p;
}

#include <glm/gtx/rotate_vector.hpp>

std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> > pumpkin_stalk() {

   float length = 1.4;
   float twist_1_d = 90.0; // degrees, around z axis
   float twist_2_d = 20.0; // degrees, around y axis
   float delta_phi_for_bevel_d = 1.0;
   float r_start = 0.2;
   float r_end   = 0.11;
   const float pi = 3.1415926f;
   
   std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> > vi;
   std::vector<position_normal_vertex>  &vertices  = vi.first;
   std::vector<g_triangle> &triangles = vi.second;

   unsigned int n_layers = 5;
   unsigned int n_slices = 28;
   float height_per_slice = length/static_cast<float>(n_layers);
   glm::vec3 y_axis(0,1,0);

   // basically a cone, with bevels at the ridge and vallies.
   for (unsigned int i=0; i<n_layers; i++) {
      // fat end at the bottom
      float f_i = static_cast<float>(i)/static_cast<float>(n_layers);
      float r_layer = r_start + f_i * (r_end - r_start);
      for (unsigned int j=0; j<n_slices; j++) {

         float f_j = static_cast<float>(j)/static_cast<float>(n_slices);
         float theta = 2.0 * pi * f_j;
         float r_this = r_layer;
         if (j%2 == 0)
            r_this = r_layer * 0.9;

         // with bevels
         float delta_theta = 2.0 * pi * delta_phi_for_bevel_d / 180.0;
         for (int jj=-1; jj<2; jj+=2) {
            float theta_this = theta + delta_theta * static_cast<float>(jj);
            float twist_this = f_i * twist_1_d * (pi/180.0);
            float theta_total = theta_this + twist_this;
            float x = r_this * sinf(theta_total);
            float y = r_this * cosf(theta_total);
            float z = static_cast<float>(i) * height_per_slice;
            glm::vec3 p(x,y,z);
            glm::vec3 n(0,0,1);  // reassigned late

            // now curl up the stalk a bit
            float angle = f_i * twist_2_d * (pi/180.0);
            glm::vec3 a = glm::rotate(p, angle, y_axis);
            p = a;
            vertices.push_back(position_normal_vertex(p,n));
         }
      }
   }

   // triangles
   for (unsigned int i=0; i<n_layers; i++) {
      for (unsigned int j=0; j<2*n_slices; j++) {
         unsigned int idx_base = i * 2 * n_slices + j;
         unsigned int idx_a = idx_base;
         unsigned int idx_b = idx_base + 1;
         unsigned int idx_c = idx_base + 2 * n_slices;
         unsigned int idx_d = idx_base + 2 * n_slices + 1;
         if (j == (2 * n_slices - 1)) {
            idx_b = i * 2 * n_slices + 1;
            idx_d = i * 2 * n_slices + 2 * n_slices + 1;
         }
         if (idx_a < vertices.size() && idx_b < vertices.size() &&
             idx_c < vertices.size() && idx_d < vertices.size()) {
            g_triangle t1(idx_a, idx_d, idx_b);
            g_triangle t2(idx_d, idx_a, idx_c);
            triangles.push_back(t1);
            triangles.push_back(t2);
         } else {
            if (false)
               std::cout << "pumpkin_stalk(): bad index "
                         << idx_a << " " << idx_b << " " << idx_c << " " << idx_d
                         << std::endl;
         }
      }
   }

   // stalk top triangles
   unsigned int n_vertices = vertices.size();
   unsigned int n_vertices_per_slice = n_slices * 2;
   unsigned int idx_stalk_top_start = n_vertices - n_vertices_per_slice;
   for (unsigned int i=idx_stalk_top_start; i<(n_vertices-2); i++) {
      unsigned int i1 = idx_stalk_top_start;
      unsigned int i2 = i+1;
      unsigned int i3 = i+2;
      g_triangle t(i1, i3, i2);
      triangles.push_back(t);
   }

   // normals of the vertices, average them.
   std::map<unsigned int, std::vector<glm::vec3> > normals_map;
   // collect the face normals, per vertex
   for (unsigned int i=0; i<triangles.size(); i++) {
      const g_triangle &t = triangles[i];
      glm::vec3 sum(0,0,0);
      for (unsigned int j=0; j<3; j++)
         sum += vertices[triangles[i].point_id[j]].pos;
      glm::vec3 n = glm::normalize(glm::cross(vertices[t[1]].pos-vertices[t[0]].pos,
                                              vertices[t[2]].pos-vertices[t[0]].pos));
      for (unsigned int j=0; j<3; j++)
         normals_map[t[j]].push_back(n);
   }
   // average the normals, for each vertex in the normals map
   for (unsigned int i=0; i<vertices.size(); i++) {
      glm::vec3 sum(0,0,0);
      std::map<unsigned int, std::vector<glm::vec3> >::const_iterator it =
         normals_map.find(i);
      if (it != normals_map.end()) {
         const std::vector<glm::vec3> &vv = it->second;
         for (unsigned int j=0; j<vv.size(); j++) {
            const glm::vec3 &v = vv[j];
            sum += v;
         }
      }
      glm::vec3 n = glm::normalize(sum);
      vertices[i].normal = n;
      // std::cout << "debug: normal for vertex " << i << " is " << glm::to_string(n) << std::endl;
   }

   return vi;
}



void
export_pumpkin_as_obj(const std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> > &vi) {

   const std::vector<position_normal_vertex>  &vertices  = vi.first;
   const std::vector<g_triangle> &triangles = vi.second;
   std::string name = "pumpkin";
   std::string file_name = "pumpkin.obj";

   glm::vec3 n(0,0,1);
   glm::vec4 c(0.7, 0.5,0.2, 1.0);

   std::cout << "export_pumpkin() vertices.size() " << vertices.size() << std::endl;
   std::cout << "export_pumpkin() triangles.size() " << triangles.size() << std::endl;

   std::vector<s_generic_vertex> g_vertices(vertices.size());
   for (unsigned int i=0; i<vertices.size(); i++) {
      g_vertices[i].pos    = vertices[i].pos;
      g_vertices[i].normal = vertices[i].normal;
      g_vertices[i].color  = c;
   }

   std::ofstream f(file_name.c_str());
   if (f) {
      f << "# " << name << "\n";
      f << "# " << "\n";
      f << "" << "\n";
      f << "g exported_obj\n";
      for (unsigned int i=0; i<vertices.size(); i++) {
         const s_generic_vertex &vert = g_vertices[i];
         f << "v " << vert.pos.x << " " << vert.pos.y << " " << vert.pos.z;
         f << " " << vert.color.r << " " << vert.color.g << " " << vert.color.b;
         f << "\n";
      }
      for (unsigned int i=0; i<vertices.size(); i++) {
         const s_generic_vertex &vert = g_vertices[i];
         f << "vn " << vert.normal.x << " " << vert.normal.y << " " << vert.normal.z << "\n";
      }
      for (unsigned int i=0; i<triangles.size(); i++) {
         const g_triangle &tri = triangles[i];
         f << "f "
           << tri.point_id[0]+1 << "//" << tri.point_id[0]+1 << " "
           << tri.point_id[1]+1 << "//" << tri.point_id[1]+1 << " "
           << tri.point_id[2]+1 << "//" << tri.point_id[2]+1 << "\n";
      }
   }
}


#ifdef STANDALONE_PUMPKIN

int
main(int argc, char **argv) {

   std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> > vi1 = pumpkin();
   std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> > vi2 = pumpkin_stalk();

   unsigned int idx_base = vi1.first.size();
   unsigned int tri_base = vi1.second.size();

    vi1.first.insert(vi1.first.end(), vi2.first.begin(), vi2.first.end());
    vi1.second.insert(vi1.second.end(), vi2.second.begin(), vi2.second.end());
    for (unsigned int i=tri_base; i<vi1.second.size(); i++)
       vi1.second[i].rebase(idx_base);

   export_pumkin(vi1);
   return 0;

}

#endif
