
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
const double pi = M_PI;

#define GLM_ENABLE_EXPERIMENTAL
#include "g_triangle.hh"
#include "cylinder-with-rotation-translation.hh"
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/string_cast.hpp>  // to_string()


// Make triangles for a cylinder along the z axis
//
// top_radius is currently ignored, cylinders only
// This should have a less generic name.
//
cylinder_with_rotation_translation::cylinder_with_rotation_translation(const std::pair<glm::vec3,
                                                                       glm::vec3> &pos_pair,
                                                                       float base_radius_in,
                                                                       float top_radius_in,
                                                                       float height_in,
                                                                       unsigned int n_slices_in,
                                                                       unsigned int n_stacks) {

   // n_stacks*n_slices = 12
   // cylinder_vertex vertices[12];
   // n_triangles = n_stacks * n_slices * 2;
   // tri_indices cylinder_indices[24];

   height = height_in;
   n_slices = n_slices_in;
   base_radius = base_radius_in;
   top_radius = top_radius_in;
   const glm::vec3 &start  = pos_pair.first;
   const glm::vec3 &finish = pos_pair.second;
   glm::vec3 b = finish - start;
   b = glm::normalize(b);
   glm::vec3 normalized_bond_orientation(b.x, b.y, b.z);
   glm::mat4 ori = glm::orientation(normalized_bond_orientation, glm::vec3(0.0, 0.0, 1.0));

   model_rotation_matrix = glm::transpose(ori);
   model_translation = start;

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

         // guess and fiddle - not proud
         float delta_radius = top_radius - base_radius;
         int int_stack = n_stacks - 1 - i_stack;
         float interpolated_radius = base_radius + delta_radius * one_over_n_stacks * static_cast<float>(int_stack);

         if (false)
            if (base_radius != top_radius)
               std::cout << "debug radii: " << i_stack << " " << base_radius << " "
                         << top_radius << " " << interpolated_radius << std::endl;
         glm::vec3 sp(x*interpolated_radius, y*interpolated_radius, z_this);
         glm::vec3 t = start;

         vertex_with_rotation_translation &v = vertices[idx];

         v.pos = sp;
         v.normal = glm::vec3(x,y,0.0f);
         v.model_rotation_matrix = glm::transpose(ori);
         v.model_translation = t;
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
cylinder_with_rotation_translation::add_flat_start_cap() {
   add_flat_cap(0.0f);
}


void
cylinder_with_rotation_translation::add_flat_end_cap() {
   add_flat_cap(height);
}

void
cylinder_with_rotation_translation::add_flat_cap(float z) {

   glm::vec3 n(0,0,1);
   if (z == 0.0f) n = -n;

   unsigned int idx_base = vertices.size();
   unsigned int idx_base_tri = triangle_indices_vec.size();

   vertex_with_rotation_translation vertex;
   vertex.pos    = glm::vec3(0,0,z);
   vertex.normal = n;
   vertex.model_rotation_matrix = model_rotation_matrix;
   vertex.model_translation = model_translation;
   vertices.push_back(vertex);

   float one_over_n_slices = 1.0/static_cast<float>(n_slices);
   float radius = base_radius;

   for (unsigned int i=0; i<n_slices; i++) {
      float theta_this = 2.0 * pi * static_cast<float>(i) * one_over_n_slices;
      float x = cosf(theta_this);
      float y = sinf(theta_this);
      glm::vec3 sp(x*radius, y*radius, z);
      vertex_with_rotation_translation v;
      v.pos    = sp;
      v.normal = n;
      v.model_rotation_matrix = model_rotation_matrix;
      v.model_translation = model_translation;
      vertices.push_back(v);
   }

   for (unsigned int i=0; i<n_slices; i++) {
      unsigned int i_next = idx_base + i + 1 + 1;
      if (i == (n_slices-1)) i_next = idx_base + 1;
      g_triangle triangle(idx_base, idx_base + i + 1, i_next);
      triangle_indices_vec.push_back(triangle);
   }


}
