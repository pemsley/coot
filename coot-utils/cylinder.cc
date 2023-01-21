
#define _USE_MATH_DEFINES
#include <cmath>
const double pi = M_PI;

#include <iostream>
#include "cylinder.hh"
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/string_cast.hpp>  // to_string()

#include "oct.hh"

// Make triangles for a cylinder along the z axis
//
// top_radius is currently ignored, cylinders only
// This should have a less generic name.
//
cylinder::cylinder(const std::pair<glm::vec3, glm::vec3> &pos_pair,
                   float base_radius_in, float top_radius_in, float height_in,
                   unsigned int n_slices_in, unsigned int n_stacks) {

   basic_colour = glm::vec4(0.5, 0.5, 0.5, 1.0);
   init(pos_pair, base_radius_in, top_radius_in, height_in, basic_colour, n_slices_in, n_stacks);

}


cylinder::cylinder(const std::pair<glm::vec3, glm::vec3> &cart_pair,
                   float base_radius, float top_radius, float height,
                   const glm::vec4 &basic_colour_in,
                   unsigned int n_slices, unsigned int n_stacks) {

   basic_colour = basic_colour_in;
   init(cart_pair, base_radius, top_radius, height, basic_colour_in, n_slices, n_stacks);
}


void
cylinder::init(const std::pair<glm::vec3, glm::vec3> &pos_pair,
               float base_radius_in, float top_radius_in, float height_in,
               const glm::vec4 &basic_colour,
               unsigned int n_slices_in, unsigned int n_stacks) {

   if (n_stacks < 2) n_stacks = 2;

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
   float bond_length = glm::distance(start, finish);
   glm::vec3 normalized = glm::normalize(b);
   // glm::orientation fails for normalized = glm::vec3(0.0, 0.0, -1.0);
   bool is_pathological = false;
   if (fabsf(normalized[0]) < 0.00001)
      if (fabsf(normalized[1]) < 0.00001)
         if (fabsf(normalized[2]) > 0.98)
            is_pathological = true;
   if (is_pathological)
      ori = glm::mat4(1,0,0,0, 0,-1,0,0, 0,0,-1,0, 0,0,0,1);
   else
      ori = glm::orientation(normalized, glm::vec3(0.0, 0.0, 1.0));

   // glm::mat4 tori = glm::transpose(ori);
   // std::cout << "ori " << glm::to_string(ori) << "\n";

   triangles.resize(n_stacks * n_slices * 2);
   vertices.resize((n_stacks + 1) * n_slices);

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
         int stack_int = n_stacks - i_stack;
         float stack_scale = 0.5 * static_cast<float>(stack_int);
         float interpolated_radius = base_radius + delta_radius * one_over_n_stacks * stack_scale;

         if (false)
            std::cout << i_stack << " stack_int " << stack_int << " stack_scale " << stack_scale
                      << " interpolated_radius " << interpolated_radius << std::endl;

         glm::vec4 p_1(x*interpolated_radius, y*interpolated_radius, z_this, 1.0f);
         // glm::vec4 p_1(x*top_radius, y*top_radius, z_this, 1.0f);
         glm::vec4 p_n(x, y, 0.0f, 1.0f);

         coot::api::vnc_vertex &v = vertices[idx];
         v.color = basic_colour;
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
         triangles[2*idx  ] = ti_1;
         triangles[2*idx+1] = ti_2;
      }
   }
   // std::cout << "Finished cylinder constructor" << std::endl;
}



void
cylinder::add_flat_start_cap() {
   add_flat_cap(0);
}


void
cylinder::add_flat_end_cap() {

   add_flat_cap(1);
}

void
cylinder::add_flat_cap(int end_type) {

   glm::vec3 n(0,0,-1);
   if (end_type == 1) n = -n;
   glm::vec4 n4(n, 1.0f);

   float z = 0.0f;
   if (end_type == 1) z = height;

   unsigned int idx_base = vertices.size();

   coot::api::vnc_vertex vert;
   vert.pos    = glm::vec3(ori * glm::vec4(0,0,z,1.0f)) + start;
   vert.normal = glm::vec3(ori * n4);
   vert.color  = basic_colour;
   vertices.push_back(vert);
   
   float one_over_n_slices = 1.0/static_cast<float>(n_slices);
   float radius = base_radius;

   for (unsigned int i=0; i<n_slices; i++) {
      float theta_this = 2.0 * pi * static_cast<float>(i) * one_over_n_slices;
      float x = cosf(theta_this);
      float y = sinf(theta_this);
      glm::vec4 p_1(x * radius, y *radius, z, 1.0);
      coot::api::vnc_vertex v;
      v.pos    = glm::vec3(ori * p_1);
      v.pos+= start;
      v.normal = glm::vec3(ori * n4);
      v.color  = basic_colour;
      vertices.push_back(v);
      // std::cout << "add flat cap vertex " << end_type << " " << i << " " << glm::to_string(v.pos) << std::endl;
   }

   for (unsigned int i=0; i<n_slices; i++) {
      unsigned int i_next = idx_base + i + 1 + 1;
      if (i == (n_slices-1)) i_next = idx_base + 1;
      // hmmm.. .which is the correct winding? Both seem to work OK. // 20211006-PE but then again I am not using backface_culling
      // g_triangle triangle(idx_base, i_next, idx_base + i + 1);
      g_triangle triangle(idx_base, idx_base + i + 1, i_next);
      triangles.push_back(triangle);
   }

}


void
cylinder::add_octahemisphere_end_cap() {

   float radius = base_radius;
   unsigned int num_subdivisions = 2;
   if (n_slices ==  8) num_subdivisions = 1;
   if (n_slices == 32) num_subdivisions = 3;
   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > hemi = tessellate_hemisphere_patch(num_subdivisions);

   std::vector<glm::vec3> &vv = hemi.first;

   std::vector<coot::api::vnc_vertex> nv(vv.size());
   for (unsigned int i=0; i<vv.size(); i++) {
      nv[i].normal = glm::vec3(ori * glm::vec4(vv[i], 1.0f));
      vv[i] *= radius;
      vv[i].z *= unstubby_rounded_cap_factor; // 20210911-PE so that when I make them short and fat scaling z down
                                              // more than x and y, the end cap can still be rounded (like a pill)
      vv[i].z += height;
      glm::vec4 p_1(ori * glm::vec4(vv[i], 1.0f));
      nv[i].pos = glm::vec3(p_1);
      nv[i].pos += start;

      nv[i].color = basic_colour;
   }

   unsigned int idx_base = vertices.size();
   unsigned int idx_base_tri = triangles.size();
   vertices.insert(vertices.end(), nv.begin(), nv.end());
   triangles.insert(triangles.end(), hemi.second.begin(), hemi.second.end());
   for (unsigned int i=idx_base_tri; i<triangles.size(); i++)
      triangles[i].rebase(idx_base);

}

void
cylinder::add_octahemisphere_start_cap() {

   unsigned int num_subdivisions = 2;
   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > hemi =
      tessellate_hemisphere_patch(num_subdivisions);

   std::vector<glm::vec3> &vv = hemi.first;

   std::vector<coot::api::vnc_vertex> nv(vv.size());
   for (unsigned int i=0; i<vv.size(); i++) {
      vv[i].z = -vv[i].z;  // 20210911-PE dangerous? (changes the winding)
      vv[i].z *= unstubby_rounded_cap_factor;
      glm::vec4 p_1(ori * glm::vec4(vv[i], 1.0f));
      nv[i].pos = glm::vec3(p_1);
      nv[i].pos *= base_radius;
      nv[i].pos += start;
      nv[i].normal = glm::vec3(p_1);
      nv[i].color = basic_colour;
   }

   unsigned int idx_base = vertices.size();
   unsigned int idx_base_tri = triangles.size();
   vertices.insert(vertices.end(), nv.begin(), nv.end());
   auto rw_triangles = hemi.second;
   for (auto &t : rw_triangles)
      t.reverse_winding();
   triangles.insert(triangles.end(), rw_triangles.begin(), rw_triangles.end());
   for (unsigned int i=idx_base_tri; i<triangles.size(); i++)
      triangles[i].rebase(idx_base);

}

void
cylinder::add_vertices_and_triangles(const std::pair<std::vector<coot::api::vnc_vertex>, std::vector<g_triangle> > &vt) {

   unsigned int idx_base = vertices.size();
   unsigned int idx_base_tri = triangles.size();
   vertices.insert(vertices.end(), vt.first.begin(), vt.first.end());
   triangles.insert(triangles.end(), vt.second.begin(), vt.second.end());
   unsigned int n_tri = triangles.size();
   for (unsigned int i=idx_base_tri; i<n_tri; i++)
      triangles[i].rebase(idx_base);
}

#include <glm/gtx/string_cast.hpp>  // to_string()

void
cylinder::add_sad_face()  {

   // I want a circle wrapped around a cylinder

   auto make_curved_circular_mesh = [] (const glm::vec4 &light_colour, const glm::vec4 &dark_colour,
                                        float x_base) {

                              std::vector<coot::api::vnc_vertex> points;
                              std::vector<g_triangle> triangles;

                              float eye_middle_z = 0.9;
                              glm::vec3 circle_middle(x_base, 0.0, eye_middle_z);
                              glm::vec3 n_middle(1,0,0);
                              coot::api::vnc_vertex cp(circle_middle, n_middle, light_colour);
                              points.push_back(cp);
                              unsigned int n_slices = 24;
                              float scale = 0.03;
                              float z_scale = 1.1/0.7; // compensates for stubby cylinder
                              for (unsigned int i=0; i<n_slices; i++) {
                                 double theta = 2.0 * M_PI * static_cast<double>(i)/static_cast<double>(n_slices);
                                 float x = x_base * (0.9 + 0.1 * fabsf(sinf(theta)));
                                 float y = scale * cosf(theta);
                                 float z = scale * sinf(theta) * z_scale + eye_middle_z;
                                 glm::vec3 pt_outer_ring(x, y, z);
                                 glm::vec3 delta = pt_outer_ring - circle_middle;
                                 glm::vec3 pt_inner_ring = circle_middle + glm::vec3(0.6, 0.6, 0.6) * delta;
                                 glm::vec3 n = glm::normalize(pt_outer_ring - circle_middle * glm::vec3(0.9, 1.0, 1.0)); // not googly-eyed
                                 coot::api::vnc_vertex v_o(pt_outer_ring, n, light_colour);
                                 coot::api::vnc_vertex v_i(pt_inner_ring, n, light_colour); // they cant be the same normal
                                 points.push_back(v_o);
                                 points.push_back(v_i);
                              }

                              // outer ring triangles (whites)
                              //
                              for (unsigned int i=0; i<n_slices; i++) {
                                 if (i+1 == n_slices) {
                                    g_triangle t1(2*i+2, 2*i+1, 1);
                                    g_triangle t2(2*i+2, 1, 2);
                                    triangles.push_back(t1);
                                    triangles.push_back(t2);
                                 } else {
                                    g_triangle t1(2*i+2, 2*i+1, 2*i+3);
                                    g_triangle t2(2*i+2, 2*i+3, 2*i+4);
                                    triangles.push_back(t1);
                                    triangles.push_back(t2);
                                 }
                              }

                              // inner ring vertices (pupils)

                              unsigned int idx_dark_mid_point = points.size();
                              glm::vec3 circle_dark_middle(x_base * 1.06, 0.0, eye_middle_z);
                              glm::vec3 n_dark_middle(1,0,0);
                              coot::api::vnc_vertex cdm(circle_dark_middle, n_middle, dark_colour);
                              points.push_back(cdm);

                              for (unsigned int i=0; i<n_slices; i++) {
                                 double theta = 2.0 * M_PI * static_cast<double>(i)/static_cast<double>(n_slices);
                                 float x = x_base * (0.9 + 0.1 * fabsf(sinf(theta)));
                                 float y = scale * cosf(theta);
                                 float z = scale * sinf(theta) * z_scale + eye_middle_z;
                                 glm::vec3 pt_outer_ring(x, y, z);
                                 glm::vec3 delta = pt_outer_ring - circle_middle;
                                 glm::vec3 pt_inner_ring = circle_middle + glm::vec3(0.6, 0.6, 0.6) * delta;
                                 glm::vec3 n = glm::normalize(pt_outer_ring - circle_middle * glm::vec3(0.9, 1.0, 1.0));
                                 coot::api::vnc_vertex v_i(pt_inner_ring, n, dark_colour);
                                 points.push_back(v_i);
                              }

                              // innner ring triangles (pupils)

                              for (unsigned int i=0; i<n_slices; i++) {
                                 if (i+1 == n_slices) {
                                    g_triangle t(idx_dark_mid_point, idx_dark_mid_point + i + 1, idx_dark_mid_point + 1);
                                    triangles.push_back(t);
                                 } else {
                                    g_triangle t(idx_dark_mid_point, idx_dark_mid_point + i + 1, idx_dark_mid_point + i + 2);
                                    triangles.push_back(t);
                                 }
                              }

                              return std::make_pair(points, triangles);
                           };

   auto rotate_points_about_z_axis = [] (std::pair<std::vector<coot::api::vnc_vertex>, std::vector<g_triangle> > &vertices_and_triangles,
                                         float angle) {  // in radians

                                        glm::vec3 z(0,0,1);
                                        for (unsigned int i=0; i<vertices_and_triangles.first.size(); i++) {
                                           auto &vertex = vertices_and_triangles.first[i];
                                           auto &pos = vertex.pos;
                                           pos = glm::rotate(pos, angle, z);
                                        }
                                     };

   auto make_mouth = [] () {
                        std::vector<coot::api::vnc_vertex> points;
                        std::vector<g_triangle> triangles;

                        float x_base = 0.065;
                        int n_per_side = 8;
                        // float z_base = 0.55; // for frown
                        float z_base = 0.85;
                        glm::vec3 n(1,0,0);
                        glm::vec4 col(0.6, 0.2, 0.1, 1.0);
                        for (int i = -n_per_side; i<=n_per_side; i++) {
                           float theta = 0.9 * static_cast<float>(i)/static_cast<float>(n_per_side);
                           float x = x_base * cosf(0.75 * theta);
                           float y_i = 0.05 * sinf(theta);
                           float y_o = y_i;
                           float z_i = z_base - 0.20 * cosf(0.5 * theta); // was + 0.2... for frown
                           float z_o = z_i + 0.03;
                           glm::vec3 v_i(x, y_i, z_i);
                           glm::vec3 v_o(x, y_o, z_o);
                           glm::vec3 mno(0.04, 0.0, 0.5 * (z_i + z_o));
                           glm::vec3 delta_i = v_i - mno;
                           glm::vec3 delta_o = v_o - mno;
                           glm::vec3 n_i = glm::normalize(delta_i);
                           glm::vec3 n_o = glm::normalize(delta_o);
                           points.push_back(coot::api::vnc_vertex(v_i, n_i, col));
                           points.push_back(coot::api::vnc_vertex(v_o, n_o, col));
                        }
                        for (int i = 0; i<2*n_per_side; i++) {
                           g_triangle t1(i*2,   i*2+1, i*2+2);
                           g_triangle t2(i*2+1, i*2+3, i*2+2);
                           triangles.push_back(t1);
                           triangles.push_back(t2);
                        }
                        return std::make_pair(points, triangles);
                   };

   glm::vec4 light_colour(0.9, 0.9, 0.9, 1.0);
   glm::vec4 dark_colour(0.1, 0.1, 0.1, 1.0);
   float x_base = 0.065;
   auto vertices_and_triangles_1 = make_curved_circular_mesh(light_colour, dark_colour, x_base);
   auto vertices_and_triangles_2 = vertices_and_triangles_1;
   rotate_points_about_z_axis(vertices_and_triangles_1,  30.0 * M_PI/180.0);
   rotate_points_about_z_axis(vertices_and_triangles_2, -30.0 * M_PI/180.0);
   add_vertices_and_triangles(vertices_and_triangles_1);
   add_vertices_and_triangles(vertices_and_triangles_2);

   float inv_fac = 1.0/255.0;
   glm::vec4 no(232 * inv_fac, 190 * inv_fac , 172 * inv_fac, 1.0f);
   no = glm::vec4(0.9, 0.4, 0.6, 1.0);
   std::vector<coot::api::vnc_vertex> no_vertices;
   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > hemi = tessellate_hemisphere_patch(2);
   glm::vec3 y(0,1,0);
   for (auto &vert : hemi.first) {
      glm::vec3 spos = vert * glm::vec3(0.03f, 0.03f, 0.03f);
      glm::vec3 rpos = glm::rotate(spos, static_cast<float>(0.5 * M_PI), y);
      glm::vec3 tpos = rpos + glm::vec3(0.06, 0, 0.8);
      glm::vec3 norm = glm::rotate(glm::normalize(vert), static_cast<float>(0.5 * M_PI), y);
      coot::api::vnc_vertex g(tpos, norm, no);
      no_vertices.push_back(g);
   }
   add_vertices_and_triangles(std::make_pair(no_vertices, hemi.second));

   auto vertices_and_triangles_m = make_mouth();
   add_vertices_and_triangles(vertices_and_triangles_m);

}

void
cylinder::z_translate(float f) {

   for (auto &v : vertices)
      v.pos.z += f;
}

void
cylinder::crenulations() {

   glm::vec4 col(0.2, 0.2, 0.2, 1.0);

   // vertices.clear();
   // triangles.clear();
   unsigned int idx_base = vertices.size();

   // top and bottom rings
   unsigned int n_stacks = 2;
   unsigned int n_slices = 80; // or so

   // add_flat_start_cap();

   // each vertex has 3 occurrances, because (although the positions are the same), the normals are different
   // (one is tangential the other is from connecting points with different radii), the third is to
   // fill the gaps on the z = (0,1) planes
   //
   std::vector<coot::api::vnc_vertex> vertices_local(n_slices * n_stacks * 3);

   for (unsigned int i_stack=0; i_stack<n_stacks; i_stack++) {
      glm::vec3 n(0,0,-1);
      if (i_stack == 1) n = glm::vec3(0,0,1);
      for (unsigned int i_slice=0; i_slice<n_slices; i_slice++) {
         unsigned int idx_0 = i_stack * n_slices * 3 + 3 * i_slice;
         unsigned int idx_1 = i_stack * n_slices * 3 + 3 * i_slice + 1;
         unsigned int idx_2 = i_stack * n_slices * 3 + 3 * i_slice + 2;

         float f_this = static_cast<float>(i_slice  )/static_cast<float>(n_slices);
         float f_prev = static_cast<float>(i_slice-1)/static_cast<float>(n_slices);
         float f_next = static_cast<float>(i_slice+1)/static_cast<float>(n_slices);
         float x_this = cosf(2.0 * M_PI * f_this);
         float y_this = sinf(2.0 * M_PI * f_this);
         float x_prev = cosf(2.0 * M_PI * f_prev);
         float y_prev = sinf(2.0 * M_PI * f_prev);
         float x_next = cosf(2.0 * M_PI * f_next);
         float y_next = sinf(2.0 * M_PI * f_next);
         float z = static_cast<float>(i_stack);
         int rem = i_slice % 4;
         float out_fac = 1.1f;
         if (rem == 0) {
            x_prev *= out_fac;
            y_prev *= out_fac;
         }
         if (rem == 1) {
            x_next *= out_fac;
            y_next *= out_fac;
         }
         if (rem == 2) {
            x_this *= out_fac;
            y_this *= out_fac;
            x_next *= out_fac;
            y_next *= out_fac;
         }
         if (rem == 3) {
            x_this *= out_fac;
            y_this *= out_fac;
         }

         // positions
         glm::vec3 pos(x_this, y_this, z);
         vertices_local[idx_0].pos = pos;
         vertices_local[idx_1].pos = pos;
         vertices_local[idx_2].pos = pos;

         // normals
         glm::vec3 delta_to_prev(x_this-x_prev, y_this-y_prev, 0.0f); // order will need fiddling
         glm::vec3 delta_to_next(x_next-x_this, y_next-y_this, 0.0f);

         float angle = - M_PI/2.0;
         glm::vec3 n1 = glm::rotate(delta_to_prev, angle, glm::vec3(0,0,1));
         glm::vec3 n2 = glm::rotate(delta_to_next, angle, glm::vec3(0,0,1));

         vertices_local[idx_0].normal = glm::normalize(n1);
         vertices_local[idx_1].normal = glm::normalize(n2);
         vertices_local[idx_2].normal = n; // for top and bottom face completion

         if (false) {
            std::cout << "vertex_local idx " << idx_0 << " position " << glm::to_string(pos) << "     normal " << glm::to_string(vertices_local[idx_0].normal) << std::endl;
            std::cout << "vertex_local idx " << idx_1 << " position " << glm::to_string(pos) << "     normal " << glm::to_string(vertices_local[idx_1].normal) << std::endl;
            std::cout << "vertex_local idx " << idx_2 << " position " << glm::to_string(pos) << "     normal " << glm::to_string(vertices_local[idx_2].normal) << std::endl;
         }

         vertices_local[idx_0].color = col;
         vertices_local[idx_1].color = col;
         vertices_local[idx_2].color = col;
      }
   }

   float z_height = 0.08;
   for (unsigned int i=0; i<vertices_local.size(); i++) {
      auto &v = vertices_local[i];
      v.pos.z *= z_height;
      v.pos.z -= z_height;
      v.pos.x *= base_radius;
      v.pos.y *= base_radius;
      // std::cout << " scale-down " << i << " " << glm::to_string(v.pos) << std::endl;
   }

   coot::api::vnc_vertex offset_origin(glm::vec3(0.0f, 0.0f, -z_height), glm::vec3(0,0,-1), col);
   vertices_local.push_back(offset_origin);

   vertices.insert(vertices.end(), vertices_local.begin(), vertices_local.end());
   unsigned int idx_offset_origin = vertices.size() - 1;

   // -------- triangles -----------

   for (unsigned int i_slice=0; i_slice<n_slices; i_slice++) {
      unsigned int idx_this_z0 = idx_base +                3 * i_slice;
      unsigned int idx_this_z1 = idx_base + 3 * n_slices + 3 * i_slice;

      unsigned int idx_1 = idx_this_z0 + 1;
      unsigned int idx_2 = idx_this_z0 + 3;
      unsigned int idx_3 = idx_this_z1 + 1;
      unsigned int idx_4 = idx_this_z1 + 3;

      if ((i_slice+1) == n_slices) { // last slice
         idx_2 = idx_base;
         idx_4 = idx_base + 3 * n_slices;
      }
      if (false)
         std::cout << "adding triangle " << idx_1 << " " << idx_2 << " " << idx_3 << " "
                   << glm::to_string(vertices[idx_1].pos) << " "
                   << glm::to_string(vertices[idx_2].pos) << " "
                   << glm::to_string(vertices[idx_3].pos) << " "
                   << std::endl;
      triangles.push_back(g_triangle(idx_1, idx_2, idx_3));
      triangles.push_back(g_triangle(idx_2, idx_3, idx_4));
   }

   // fill in the horizontal gaps at the top and bottom (face completion)
   //
   for (unsigned int i_stack=0; i_stack<n_stacks; i_stack++) {
      for (unsigned int i_slice=1; i_slice<n_slices; i_slice += 4) {

         unsigned int idx_1 = idx_base + i_stack * 3 * n_slices + 3 * i_slice       + 2;
         unsigned int idx_2 = idx_base + i_stack * 3 * n_slices + 3 * (i_slice + 1) + 2;
         unsigned int idx_3 = idx_base + i_stack * 3 * n_slices + 3 * (i_slice + 2) + 2;
         unsigned int idx_4 = idx_base + i_stack * 3 * n_slices + 3 * (i_slice + 3) + 2;

         if (i_slice >= (n_slices-3)) {
            idx_4 = idx_base + i_stack * 3 * n_slices + 3 * (0) + 2;
         }

         if ((i_slice+1) == n_slices) { // last slice
            idx_2 = idx_base + i_slice * 3 * n_slices + 2;
            idx_4 = idx_base + i_slice * 3 * n_slices;
         }

         // std::cout << "adding triangle " << idx_1 << " " << idx_2 << " " << idx_3 << " "
         // << glm::to_string(vertices[idx_1].pos) << " "
         // << glm::to_string(vertices[idx_2].pos) << " "
         // << glm::to_string(vertices[idx_3].pos) << " "
         // << std::endl;

         triangles.push_back(g_triangle(idx_1, idx_2, idx_3));
         triangles.push_back(g_triangle(idx_1, idx_3, idx_4));
      }
   }

   // the flat start cap
   for (unsigned int i_slice=0; i_slice<n_slices; i_slice++) {
      unsigned int idx_this = idx_base + (i_slice    ) * 3 + 2;
      unsigned int idx_next = idx_base + (i_slice + 1) * 3 + 2;
      if ((i_slice+1) == n_slices) { // last slice
         idx_next = idx_base + (0) * 3 + 2;
      }
      triangles.push_back(g_triangle(idx_offset_origin, idx_this, idx_next));
   }
}
