
#define _USE_MATH_DEFINES
#include <cmath>
const double pi = M_PI;

#include <iostream>
#include "cylinder.hh"
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/rotate_vector.hpp>

#include "oct.hh"

// Make triangles for a cylinder along the z axis
//
// top_radius is currently ignored, cylinders only
// This should have a less generic name.
//
cylinder::cylinder(const std::pair<glm::vec3, glm::vec3> &pos_pair,
                   float base_radius_in, float top_radius_in, float height_in,
                   unsigned int n_slices_in, unsigned int n_stacks) {

      // n_stacks*n_slices = 12
      // cylinder_vertex vertices[12];
      // n_triangles = n_stacks * n_slices * 2;
      // tri_indices cylinder_indices[24];

   n_slices = n_slices_in;
   height = height_in;
   base_radius = base_radius_in;
   top_radius = top_radius_in;

   start = pos_pair.first; // save start
   const glm::vec3 &finish = pos_pair.second;
   glm::vec3 b = finish - start;
   glm::vec3 normalized = glm::normalize(b);
   ori = glm::orientation(normalized, glm::vec3(0.0, 0.0, 1.0));
   glm::mat4 tori = glm::transpose(ori);

   // std::cout << "ori " << glm::to_string(ori) << "\n";

   triangle_indices_vec.resize(n_stacks * n_slices * 2);
   vertices.resize((n_stacks + 1) * n_slices);
   unsigned int idx = 0;

   float one_over_n_slices = 1.0/static_cast<float>(n_slices);
   float one_over_n_stacks = 1.0/static_cast<float>(n_stacks-1);
   float height_step = height/static_cast<float>(n_stacks);

   for (unsigned int i_stack=0; i_stack<=n_stacks; i_stack++) {
      for (unsigned int i_slice=0; i_slice<n_slices; i_slice++) {
         float z_this = i_stack * height_step;
         if (i_stack == n_stacks) z_this = height;
         float theta_this = 2.0 * pi * static_cast<float>(i_slice) * one_over_n_slices;

         float x = cosf(theta_this);
         float y = sinf(theta_this);

         int idx = i_stack*n_slices + i_slice;

         float delta_radius = top_radius - base_radius;
         int int_stack = n_stacks - 1 - i_stack;
         float interpolated_radius = base_radius + delta_radius * one_over_n_stacks * static_cast<float>(int_stack);

         // glm::vec4 p_1(x*interpolated_radius, y*interpolated_radius, z_this, 1.0f);
         glm::vec4 p_1(x*top_radius, y*top_radius, z_this, 1.0f);
         glm::vec4 p_n(x, y, 0.0f, 1.0f);

         s_generic_vertex &v = vertices[idx];
         v.pos = glm::vec3(ori * p_1);
         v.pos += start;
         v.normal = ori * p_n;
         idx++;
      }
   }

   // indices
   for (unsigned int i_stack=0; i_stack<n_stacks; i_stack++) {
      unsigned int i_slice_last = n_slices-1;
      for (unsigned int i_slice=0; i_slice<n_slices; i_slice++) {
         unsigned int idx_0 = i_stack*n_slices + i_slice;
         unsigned int idx_1 = i_stack*n_slices + i_slice + 1;
         unsigned int idx_2 = (i_stack+1)*n_slices + i_slice;
         unsigned int idx_3 = (i_stack+1)*n_slices + i_slice + 1;
         if (i_slice == i_slice_last) {
            idx_1 = i_stack*n_slices;
            idx_3 = (i_stack+1)*n_slices;
         }

         // std::cout << "t1: " << idx_0 << " " << idx_1 << " " << idx_2 << std::endl;
         // std::cout << "t2: " << idx_1 << " " << idx_3 << " " << idx_2 << std::endl;
         g_triangle ti_1(idx_0, idx_1, idx_2);
         g_triangle ti_2(idx_1, idx_3, idx_2); // both clockwise
         int idx = i_stack*n_slices + i_slice;
         // std::cout << "   indices idx " << idx << " " << triangle_indices_vec.size() << std::endl;
         triangle_indices_vec[2*idx  ] = ti_1;
         triangle_indices_vec[2*idx+1] = ti_2;
      }
   }
   // std::cout << "Finished cylinder constructor" << std::endl;
}




void
cylinder::add_flat_start_cap() {
   add_flat_cap(0.0f);
}


void
cylinder::add_flat_end_cap() {

   add_flat_cap(height);
}

void
cylinder::add_flat_cap(float z) {

   glm::vec3 n(0,0,1);
   if (z == 0.0f) n = -n;
   glm::vec4 n4(n, 1.0f);

   unsigned int idx_base = vertices.size();
   unsigned int idx_base_tri = triangle_indices_vec.size();

   s_generic_vertex vert;
   vert.pos    = glm::vec3(ori * glm::vec4(0,0,z,1.0f)) + start;
   vert.normal = glm::vec3(ori * n4);
   vertices.push_back(vert);

   float one_over_n_slices = 1.0/static_cast<float>(n_slices);
   float radius = base_radius;

   for (unsigned int i=0; i<n_slices; i++) {
      float theta_this = 2.0 * pi * static_cast<float>(i) * one_over_n_slices;
      float x = cosf(theta_this);
      float y = sinf(theta_this);
      glm::vec4 p_1(x*radius, y*radius, z, 1.0f);
      s_generic_vertex v;
      v.pos = glm::vec3(ori * p_1);
      v.pos += start;
      v.normal = glm::vec3(ori * n4);
      vertices.push_back(v);
   }

   for (unsigned int i=0; i<n_slices; i++) {
      unsigned int i_next = idx_base + i + 1 + 1;
      if (i == (n_slices-1)) i_next = idx_base + 1;
      g_triangle triangle(idx_base, idx_base + i + 1, i_next);
      triangle_indices_vec.push_back(triangle);
   }


}


void
cylinder::add_octahemisphere_end_cap() {

   float radius = base_radius;
   unsigned int num_subdivisions = 2;
   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > hemi =
      tessellate_hemisphere_patch(num_subdivisions);

   std::vector<glm::vec3> &vv = hemi.first;

   std::vector<s_generic_vertex> nv(vv.size());
   for (unsigned int i=0; i<vv.size(); i++) {
      nv[i].normal = glm::vec3(ori * glm::vec4(vv[i], 1.0f));
      vv[i] *= radius;
      vv[i].z += height;
      glm::vec4 p_1(ori * glm::vec4(vv[i], 1.0f));
      nv[i].pos = glm::vec3(p_1);
      nv[i].pos += start;
   }

   unsigned int idx_base = vertices.size();
   unsigned int idx_base_tri = triangle_indices_vec.size();
   vertices.insert(vertices.end(), nv.begin(), nv.end());
   triangle_indices_vec.insert(triangle_indices_vec.end(), hemi.second.begin(), hemi.second.end());
   for (unsigned int i=idx_base_tri; i<triangle_indices_vec.size(); i++)
      triangle_indices_vec[i].rebase(idx_base);

}

void
cylinder::add_octahemisphere_start_cap() {

   unsigned int num_subdivisions = 2;
   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > hemi =
      tessellate_hemisphere_patch(num_subdivisions);

   std::vector<glm::vec3> &vv = hemi.first;

   std::vector<s_generic_vertex> nv(vv.size());
   for (unsigned int i=0; i<vv.size(); i++) {
      vv[i].z = -vv[i].z;
      glm::vec4 p_1(ori * glm::vec4(vv[i], 1.0f));
      nv[i].pos = glm::vec3(p_1);
      nv[i].pos *= base_radius;
      nv[i].pos += start;
      nv[i].normal = glm::vec3(p_1);
   }

   unsigned int idx_base = vertices.size();
   unsigned int idx_base_tri = triangle_indices_vec.size();
   vertices.insert(vertices.end(), nv.begin(), nv.end());
   triangle_indices_vec.insert(triangle_indices_vec.end(), hemi.second.begin(), hemi.second.end());
   for (unsigned int i=idx_base_tri; i<triangle_indices_vec.size(); i++)
      triangle_indices_vec[i].rebase(idx_base);

}
