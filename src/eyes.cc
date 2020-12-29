
#include <iostream>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include "eyes.hh"

// for instanced objects (those animated/moving) the position in molecular space
// and the colour is dictacted by the instancing matrices (and colours)
// We don't need to create them here.
std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> >
make_spherical_surface_circular_patch(float radius_ball,
                                      float radius_patch,
                                      float h_scale,
                                      float v_scale,
                                      unsigned int n_slices) {

   std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> > p;
   std::vector<position_normal_vertex> &vertices = p.first;
   std::vector<g_triangle> &triangles = p.second;

   unsigned int n_rings = 3;

   auto checked_add = [] (const g_triangle &g, std::vector<g_triangle> &triangles,
                          unsigned int n_vertices) {
                         bool done = false;
                         if (g[0] < n_vertices)
                            if (g[1] < n_vertices)
                               if (g[2] < n_vertices) {
                                  triangles.push_back(g);
                                  done = true;
                               }
                         if (!done)
                            std::cout << "checked_add: reject "
                                      << g[0] << " "
                                      << g[1] << " "
                                      << g[2] << " "
                                      << std::endl;
                      };

   // vertices

   glm::vec3 o(0, 0, radius_ball);
   glm::vec3 n(0, 0, 1);
   vertices.push_back(position_normal_vertex(o,n)); // 0th index

   float sf = 1.0/static_cast<float>(n_slices);
   for (unsigned int iring=1; iring<=2; iring++) {
      float ring_scale = radius_patch * static_cast<float>(iring);
      for (unsigned int islice=0; islice<n_slices; islice++) {
         float theta = 2.0f * M_PI * sf * static_cast<float>(islice);

         float x = ring_scale * h_scale * sinf(theta);
         float y = ring_scale * v_scale * cosf(theta);
         float z = radius_ball;
         glm::vec3 vec_uv = glm::normalize(glm::vec3(x,y,z));
         glm::vec3 pos = radius_ball * vec_uv;
         position_normal_vertex v(pos, vec_uv);
         vertices.push_back(v);
      }
   }
   unsigned int n_vertices = vertices.size();

   // indices

   for (unsigned int iring=1; iring<=2; iring++) {
      unsigned int ring_offset_inner = (iring-2) * n_slices + 1;
      unsigned int ring_offset_outer = (iring-1) * (n_slices + 1);
      for (unsigned int islice=0; islice<n_slices; islice++) {
         if (iring == 1) {
            // the first/inside ring connects to the centre
            unsigned int i_this = islice + 1;
            unsigned int i_next = i_this + 1;
            if (i_this == n_slices) i_next = 1;
            g_triangle t(0, i_this, i_next);
            checked_add(t, p.second, n_vertices);
         } else {
            // outer rings connect to the ring inside it, quads
            unsigned int i_this_inner = ring_offset_inner + islice;
            unsigned int i_next_inner = i_this_inner + 1;
            unsigned int i_this_outer = ring_offset_outer + islice;
            unsigned int i_next_outer = i_this_outer + 1;

            if (islice == (n_slices - 1)) {
               i_next_inner = ring_offset_inner;
               i_next_outer = ring_offset_outer;
            }

            g_triangle g1(i_this_inner, i_next_inner, i_next_outer);
            g_triangle g2(i_next_outer, i_this_outer, i_this_inner);

            checked_add(g1, p.second, n_vertices);
            checked_add(g2, p.second, n_vertices);

         }
      }
   }

   return p;

}

// a curved bar on the surface of a sphere
//
std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> >
make_curved_bar_patch(float atom_radius, float phi_start, float phi_end, float m_height,
                      unsigned int n_slices, float curve) {

   // m_height is the mouth "width" :-) (up/down, not side to side)

   // phi is azimuthal angle
   // theta is inclination (polar) angle

   std::pair<std::vector<position_normal_vertex>, std::vector<g_triangle> > p;
   std::vector<position_normal_vertex> &vertices = p.first;
   std::vector<g_triangle> &triangles = p.second;

   float phi_step = (phi_end - phi_start)/static_cast<float>(n_slices);
   for (float phi=phi_start; phi<phi_end; phi+=phi_step) {
      float phi_this = phi;
      float phi_next = phi + phi_step;
      float theta_this = -0.71f;
      float theta_next = theta_this - m_height;
      unsigned int idx_base = vertices.size();

      glm::vec3 p1(sinf(theta_this) * cos(phi_this),
                   sinf(theta_this) * sin(phi_this),
                   cosf(theta_this));
      glm::vec3 p2(sinf(theta_this) * cos(phi_next),
                   sinf(theta_this) * sin(phi_next),
                   cosf(theta_this));
      glm::vec3 p3(sinf(theta_next) * cos(phi_this),
                   sinf(theta_next) * sin(phi_this),
                   cosf(theta_next));
      glm::vec3 p4(sinf(theta_next) * cos(phi_next),
                   sinf(theta_next) * sin(phi_next),
                   cosf(theta_next));

      p1 = glm::rotate(p1, glm::radians(90.0f), glm::vec3(0,0,1));
      p2 = glm::rotate(p2, glm::radians(90.0f), glm::vec3(0,0,1));
      p3 = glm::rotate(p3, glm::radians(90.0f), glm::vec3(0,0,1));
      p4 = glm::rotate(p4, glm::radians(90.0f), glm::vec3(0,0,1));

      position_normal_vertex v1(atom_radius * p1, p1);
      position_normal_vertex v2(atom_radius * p2, p2);
      position_normal_vertex v3(atom_radius * p3, p3);
      position_normal_vertex v4(atom_radius * p4, p4);

      vertices.push_back(v1);
      vertices.push_back(v2);
      vertices.push_back(v3);
      vertices.push_back(v4);

      if (false) {
         std::cout << "phi " << phi << " p1 " << glm::to_string(p1) << std::endl;
         std::cout << "phi " << phi << " p2 " << glm::to_string(p2) << std::endl;
         std::cout << "phi " << phi << " p3 " << glm::to_string(p3) << std::endl;
         std::cout << "phi " << phi << " p4 " << glm::to_string(p4) << std::endl;
      }

      // maybe triangles should be done differently
      triangles.push_back(g_triangle(idx_base, idx_base + 2, idx_base + 1));
      triangles.push_back(g_triangle(idx_base + 1, idx_base + 2, idx_base + 3));
   }

   return p;
}

// peritonial dialysis
//
