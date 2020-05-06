
#define _USE_MATH_DEFINES
#include <cmath>
const double pi = M_PI;

#include "cylinder.hh"
#include <glm/gtx/rotate_vector.hpp>

// Make triangles for a cylinder along the z axis
//
// top_radius is currently ignored, cylinders only
// This should have a less generic name.
//
cylinder::cylinder(const coot::CartesianPair &pospair,
                   float base_radius, float top_radius, float height,
                   unsigned int n_slices, unsigned int n_stacks) {

      // n_stacks*n_slices = 12
      // cylinder_vertex vertices[12];
      // n_triangles = n_stacks * n_slices * 2;
      // tri_indices cylinder_indices[24];

      const coot::Cartesian &start  = pospair.getStart();
      const coot::Cartesian &finish = pospair.getFinish();
      coot::Cartesian b = finish - start;
      b.unit_vector_yourself();
      glm::vec3 normal(b.x(), b.y(), b.z());
      glm::mat4 ori = glm::orientation(normal, glm::vec3(0.0, 0.0, 1.0));

      // std::cout << "ori " << glm::to_string(ori) << "\n";

      triangle_indices_vec.resize(n_stacks * n_slices * 2);
      vertices.resize((n_stacks + 1) * n_slices);
      unsigned int idx = 0;

      float one_over_n_slices = 1.0/static_cast<float>(n_slices);
      float one_over_n_stacks = 1.0/static_cast<float>(n_stacks-1);
      float height_step = height/static_cast<float>(n_stacks);
      for (unsigned int i_stack=0; i_stack<=n_stacks; i_stack++) {
         // std::cout << "i_stack " << i_stack << std::endl;
         for (unsigned int i_slice=0; i_slice<n_slices; i_slice++) {
            // std::cout << "   i_slice " << i_slice << std::endl;
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

            if (base_radius != top_radius)
               std::cout << "debug radii: " << i_stack << " " << base_radius << " "
                         << top_radius << " " << interpolated_radius << std::endl;
            glm::vec3 sp(x*interpolated_radius, y*interpolated_radius, z_this);
            glm::vec3 t(start.x(), start.y(), start.z());

            generic_vertex &v = vertices[idx];

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
            tri_indices ti_1(idx_0, idx_1, idx_2);
            tri_indices ti_2(idx_1, idx_3, idx_2); // both clockwise
            int idx = i_stack*n_slices + i_slice;
            // std::cout << "   indices idx " << idx << " " << triangle_indices_vec.size() << std::endl;
            triangle_indices_vec[2*idx  ] = ti_1;
            triangle_indices_vec[2*idx+1] = ti_2;
         }
      }
      // std::cout << "Finished cylinder constructor" << std::endl;
   }


