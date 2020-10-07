
#include <mmdb2/mmdb_manager.h>
#include "Python.h"

// Perhaps this (needed to compile meshed-generic-display-object.hh) should go in Mesh.hh.
#include <epoxy/gl.h>
#ifndef __APPLE__
#include <epoxy/glx.h>
#endif

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>  // to_string()
#include <glm/gtx/rotate_vector.hpp>

#include "meshed-generic-display-object.hh"
#include "oct.hh"
#include "graphics-info.h"
#include "cylinder.hh"

glm::vec3 coord_orth_to_glm(const clipper::Coord_orth &co) {
   return glm::vec3(co.x(), co.y(), co.z());
}

void
meshed_generic_display_object::add_line(const coot::colour_holder &colour,
                                        const std::string &colour_name, int line_width,
                                        const std::pair<clipper::Coord_orth, clipper::Coord_orth> &coords) {
}


void
meshed_generic_display_object::add_point(const coot::colour_holder &colour_in,
                                         const std::string &colour_name,
                                         const int &size_in,
                                         const clipper::Coord_orth &coords_in) {

   unsigned int num_subdivisions = 3;
   float radius = 0.03 * size_in; // changing the scaling is fun
   glm::vec4 col(colour_in.red, colour_in.green, colour_in.blue, 1.0);
   glm::vec3 position = coord_orth_to_glm(coords_in);
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
      oct = wrapped_make_octasphere(num_subdivisions, position, radius, col);
   mesh.import(oct);

}

void
meshed_generic_display_object::add_arrow(const arrow_t &arrow) {

   unsigned int n_slices = 8;
   float arrow_radius = 0.1;
   std::pair<glm::vec3, glm::vec3 > start_end(coord_orth_to_glm(arrow.start_point),
                                              coord_orth_to_glm(arrow.end_point));
   add_cylinder(start_end, arrow.col, arrow.radius, n_slices, true, true, FLAT_CAP, FLAT_CAP);

}


std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
meshed_generic_display_object::wrapped_make_octasphere(unsigned int num_subdivisions,
                                                       const glm::vec3 &position,
                                                       float radius,
                                                       const glm::vec4 &col) {

   std::map<unsigned int, std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > >::iterator it;

   unsigned int map_index = num_subdivisions + static_cast<int> (radius * 100.0f);
   it = origin_octasphere_map.find(map_index);
   if (it != origin_octasphere_map.end()) {
      std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > oct = it->second;
      for (unsigned int i=0; i<oct.first.size(); i++) {
         s_generic_vertex &v(oct.first[i]);
         v.pos += position;
         v.color = col;
      }
      return oct;
   } else {
      glm::vec3 origin(0,0,0);
      std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > oct =
         make_octasphere(num_subdivisions, origin, radius, col);
      origin_octasphere_map[map_index] = oct;
      for (unsigned int i=0; i<oct.first.size(); i++) {
         s_generic_vertex &v(oct.first[i]);
         v.pos += position;
         v.color = col;
      }
      return oct;
   }
}

void
meshed_generic_display_object::init(const graphical_bonds_container &bonds_box,
                                    bool background_is_black_flag) {

   mesh.clear();
   for (int i=0; i<bonds_box.num_colours; i++) {
      unsigned int n_slices = 16;
      graphical_bonds_lines_list<graphics_line_t> &ll = bonds_box.bonds_[i];
      bool do_thinning = ll.thin_lines_flag;
      float dark_bg_cor = 0.0;
      if (! background_is_black_flag)
         dark_bg_cor = 0.29;
      float fidx = static_cast<float>(i);
      coot::colour_holder col (0.8-dark_bg_cor, 0.8-0.4*fidx-dark_bg_cor, 0.4+0.5*fidx-dark_bg_cor);
      unsigned int n_segments = 8;
      for (int j=0; j< ll.num_lines; j++) {
         glm::vec3 s = cartesian_to_glm(ll.pair_list[j].positions.getStart());
         glm::vec3 e = cartesian_to_glm(ll.pair_list[j].positions.getFinish());
         glm::vec3 delta = e - s;
         glm::vec3 delta_frag = (1.0f / static_cast<float>(n_segments)) * delta;
         for(unsigned int iseg=0; iseg<n_segments; iseg++) {
            glm::vec3 pos_1 = s + static_cast<float>(iseg) * delta_frag;
            glm::vec3 pos_2 = pos_1 + 0.5f * delta_frag;
            std::pair<glm::vec3, glm::vec3> start_end(pos_1, pos_2);
            add_cylinder(start_end, col, 0.06, n_slices, true, true, FLAT_CAP, FLAT_CAP);
         }
      }
   }
}


void
meshed_generic_display_object::add_cylinder(const std::pair<glm::vec3, glm::vec3 > &start_end,
                                            const coot::colour_holder &col, float line_radius,
                                            unsigned int n_slices,
                                            bool cap_start, bool cap_end,
                                            cap_type start_cap_type, cap_type end_cap_type) {

   float h = glm::distance(start_end.first, start_end.second);
   cylinder c(start_end, line_radius, line_radius, h, n_slices, 2);
   glm::vec4 colour(col.red, col.green, col.blue, 1.0f);
   if (false)
      std::cout << "add_cylinder: " << glm::to_string(start_end.first) << " "
                << glm::to_string(start_end.second) << " "
                << c.vertices.size() << " " << c.triangle_indices_vec.size()
                << " with height " << h << std::endl;

   if (cap_start) {
      if (start_cap_type == FLAT_CAP)
         c.add_flat_start_cap();
      if (start_cap_type == ROUNDED_CAP)
         c.add_octahemisphere_start_cap();

   }
   if (cap_end) {
      if (end_cap_type == FLAT_CAP)
         c.add_flat_end_cap();
      if (end_cap_type == ROUNDED_CAP)
         c.add_octahemisphere_end_cap();
   }

   for (unsigned int i=0; i<c.vertices.size(); i++)
      c.vertices[i].color = colour;
   mesh.import(c.vertices, c.triangle_indices_vec);

}

void
meshed_generic_display_object::add_cone(const std::pair<glm::vec3, glm::vec3> &start_end,
                                        const coot::colour_holder &col, float base_radius, float top_radius,
                                        unsigned int n_slices,
                                        bool cap_start, bool cap_end,
                                        cap_type start_cap_type, cap_type end_cap_type) {

   float h = glm::distance(start_end.first, start_end.second);
   cylinder c(start_end, base_radius, top_radius, h, n_slices, 2);
   glm::vec4 colour(col.red, col.green, col.blue, 1.0f);

   // are your start and end points where you think they are?
   if (false)
      std::cout << "add_cone: " << glm::to_string(start_end.first) << " "
                << glm::to_string(start_end.second)
                << " base_radius " << base_radius << " top_radius " << top_radius << " "
                << c.vertices.size() << " " << c.triangle_indices_vec.size()
                << " with height " << h << std::endl;

   if (cap_start) {
      if (start_cap_type == FLAT_CAP)
         c.add_flat_start_cap();
      if (start_cap_type == ROUNDED_CAP)
         c.add_octahemisphere_start_cap();

   }
   if (cap_end) {
      if (end_cap_type == FLAT_CAP)
         c.add_flat_end_cap();
      if (end_cap_type == ROUNDED_CAP)
         c.add_octahemisphere_end_cap();
   }

   for (unsigned int i=0; i<c.vertices.size(); i++)
      c.vertices[i].color = colour;

   mesh.import(c.vertices, c.triangle_indices_vec);

}



void
meshed_generic_display_object::add_dodecahedron(const coot::colour_holder &colour_in,
                                                const std::string &colour_name,
                                                double radius, const clipper::Coord_orth &pos) {
}

void
meshed_generic_display_object::add_pentakis_dodecahedron(const coot::colour_holder &colour_in,
                                                         const std::string &colour_name,
                                                         double stellation_factor,
                                                         double radius,
                                                         const clipper::Coord_orth &pos) {

}

glm::vec3 rotate_around_vector(const glm::vec3 &direction,
                               const glm::vec3 &position,
                               const glm::vec3 &origin_shift,
                               float angle) {

   glm::vec3 p1 = position - origin_shift;
   glm::vec3 p2 = glm::rotate(p1, angle, direction);
   glm::vec3 p3 = p2 + origin_shift;
   return p3;
}

glm::vec4
colour_holder_to_glm(const coot::colour_holder &ch) {
   return glm::vec4(ch.red, ch.green, ch.blue, 1.0f);
}


void
meshed_generic_display_object::add_arc(const arc_t &arc) {

   // arc is now an elaboration on the torus code.

   // phi goes around the ring viewed from the top (down the z axis)
   // theta goes around eaach segment

   // torus.radius_1 is R
   // torus.radius_2 is r

   // x = (R + r * cos(theta)) * cos(phi);
   // y = (R + r * cos(theta)) * sin(phi);
   // z = r * sin(theta)
   //
   // the ring around the middle of the doughnut:
   // x_ring = R * cos(phi)
   // y_ring = R * sin(phi)
   // z_ring = 0;

   const unsigned int n_phi_steps = 20;
   const unsigned int n_theta_steps = 12;
   std::vector<s_generic_vertex> vertices((n_theta_steps + 1) * n_phi_steps);
   std::vector<g_triangle> triangles;
   const float R = arc.radius;
   const float r = arc.radius_inner;
   const float pi = 3.1415926535;

   float angle_delta_deg = arc.delta_angle;
   float angle_delta = angle_delta_deg * pi / 180.0f;
   glm::vec3 start_point = coord_orth_to_glm(arc.start_point);
   glm::vec4 col = colour_holder_to_glm(arc.col);
   glm::mat4 rot_mat = glm::orientation(coord_orth_to_glm(arc.normal), glm::vec3(0.0, 0.0, 1.0));
   std::cout << "rot_mat: " << glm::to_string(rot_mat)<< std::endl;

   const clipper::Mat33<double> &m = arc.orientation_matrix;
   glm::mat3 ori_mat(m(0,0), m(0,1), m(0,2),
                     m(1,0), m(1,1), m(1,2),
                     m(2,0), m(2,1), m(2,2));

   for (unsigned int ip=0; ip<=n_phi_steps; ip++) {
      float phi_raw = angle_delta * static_cast<float>(ip)/static_cast<float>(n_phi_steps);
      float phi = phi_raw;
      for (unsigned int it=0; it<n_theta_steps; it++) {
         float theta = 2.0f * pi * static_cast<float>(it)/static_cast<float>(n_theta_steps);
         s_generic_vertex v;
         v.pos.x = (R + r * cosf(theta)) * cosf(phi);
         v.pos.y = (R + r * cosf(theta)) * sinf(phi);
         v.pos.z = r * sinf(theta);
         v.normal.x = cosf(theta) * cosf(phi);
         v.normal.y = cosf(theta) * sinf(phi);
         v.normal.z = sinf(theta);
         // now move orient the fragment and normal
         // v.pos    = glm::vec3(rot_mat * glm::vec4(v.pos,    1.0f));
         // v.normal = glm::vec3(rot_mat * glm::vec4(v.normal, 1.0f));;
         v.pos    = ori_mat * v.pos;
         v.normal = ori_mat * v.normal;
         v.pos += start_point;
         v.color = col;
         unsigned int vertex_idx = ip * n_theta_steps + it;
         vertices[vertex_idx] = v;
      }
   }

   // carefully, carefully :-)
   for (unsigned int ip=0; ip<n_phi_steps; ip++) {
      unsigned int ip_this = ip;
      unsigned int ip_next = ip + 1;
      unsigned int idx_base_phi_00 = ip_this * n_theta_steps;
      unsigned int idx_base_phi_10 = ip_next * n_theta_steps;
      for (unsigned int it=0; it<n_theta_steps; it++) {
         unsigned int it_this = it;
         unsigned int it_next = it + 1;
         if (it_next == n_theta_steps) it_next = 0;

         unsigned int idx_00 = idx_base_phi_00 + it_this;
         unsigned int idx_01 = idx_base_phi_00 + it_next;
         unsigned int idx_10 = idx_base_phi_10 + it_this;
         unsigned int idx_11 = idx_base_phi_10 + it_next;

         g_triangle t1(idx_00, idx_10, idx_11);
         g_triangle t2(idx_00, idx_11, idx_01);
         triangles.push_back(t1);
         triangles.push_back(t2);
      }
   }

   mesh.import(vertices, triangles);

}

void meshed_generic_display_object::add_torus(const meshed_generic_display_object::torus_t &torus) {

   // phi goes around the ring viewed from the top (down the z axis)
   // theta goes around eaach segment

   // torus.radius_1 is R
   // torus.radius_2 is r

   // x = (R + r * cos(theta)) * cos(phi);
   // y = (R + r * cos(theta)) * sin(phi);
   // z = r * sin(theta)
   //
   // the ring around the middle of the doughnut:
   // x_ring = R * cos(phi)
   // y_ring = R * sin(phi)
   // z_ring = 0;

   const unsigned int n_theta_steps = 60;
   const unsigned int n_phi_steps = 60;
   std::vector<s_generic_vertex> vertices(n_theta_steps * n_phi_steps);
   std::vector<g_triangle> triangles;
   const float R = torus.radius_1;
   const float r = torus.radius_2;
   const float pi = 3.1415926535;
   glm::vec4 col = colour_holder_to_glm(torus.col);

   glm::vec4 centre(torus.position.x(), torus.position.y(), torus.position.z(), 1.0);
   glm::vec3 ring_normal(torus.normal.x(), torus.normal.y(), torus.normal.z());
   glm::mat4 ori = glm::orientation(ring_normal, glm::vec3(0,0,1));

   for (unsigned int ip=0; ip<n_phi_steps; ip++) {
      float phi = 2.0f * pi * static_cast<float>(ip)/static_cast<float>(n_phi_steps);
      for (unsigned int it=0; it<n_theta_steps; it++) {
         float theta = 2.0f * pi * static_cast<float>(it)/static_cast<float>(n_theta_steps);
         s_generic_vertex v;
         glm::vec4 pos;
         pos.x = (R + r * cosf(theta)) * cosf(phi);
         pos.y = (R + r * cosf(theta)) * sinf(phi);
         pos.z = r * sinf(theta);
         pos.w = 1.0; // or 0?
         v.pos = glm::vec3(pos * ori) + glm::vec3(centre);
         glm::vec4 normal;
         normal.x = cosf(theta) * cosf(phi);
         normal.y = cosf(theta) * sinf(phi);
         normal.z = sinf(theta);
         normal.w = 1.0;
         v.normal = glm::vec3(normal * ori);
         v.color = col;
         vertices[ip * n_theta_steps + it] = v;
      }
   }

   // carefullly, carefully :-)
   for (unsigned int ip=0; ip<n_phi_steps; ip++) {
      unsigned int ip_this = ip;
      unsigned int ip_next = ip + 1;
      if (ip_next == n_phi_steps) ip_next = 0;
      unsigned int idx_base_phi_00 = ip_this * n_theta_steps;
      unsigned int idx_base_phi_10 = ip_next * n_theta_steps;
      for (unsigned int it=0; it<n_theta_steps; it++) {
         unsigned int it_this = it;
         unsigned int it_next = it + 1;
         if (it_next == n_theta_steps) it_next = 0;

         unsigned int idx_00 = idx_base_phi_00 + it_this;
         unsigned int idx_01 = idx_base_phi_00 + it_next;
         unsigned int idx_10 = idx_base_phi_10 + it_this;
         unsigned int idx_11 = idx_base_phi_10 + it_next;

         g_triangle t1(idx_00, idx_10, idx_11);
         g_triangle t2(idx_00, idx_11, idx_01);
         triangles.push_back(t1);
         triangles.push_back(t2);
      }
   }

   mesh.import(vertices, triangles);

}


coot::colour_holder
colour_values_from_colour_name(const std::string &c) {

   coot::colour_holder colour;
   colour.red = 0.4;
   colour.green = 0.4;
   colour.blue = 0.4;

   if (c.length() == 7) {
      if (c[0] == '#') {
         return coot::colour_holder(c); // hex colour string
      }
   }

   if (c == "blue") {
      colour.red = 0.1; colour.green = 0.1;
      colour.blue = 0.8;
   } else {
      if (c == "sky") {
         colour.red = 0.53 * 0.6;
         colour.green = 0.81 * 0.6;
         colour.blue = 0.92 * 0.6;
      } else {
         if (c == "green") {
            colour.red   = 0.05;
            colour.green = 0.8;
            colour.blue  = 0.05;
         } else {
            if (c == "greentint") {  // old Hydrogen-bond colour
               colour.red = 0.3;
               colour.green = 0.35;
               colour.blue = 0.3;
            } else {
               if (c == "darkpurple") { // new Hydrogen-bond colour
                  colour.red   = 0.38;
                  colour.green = 0.05;
                  colour.blue  = 0.4;
               } else {
                  if (c == "sea") {
                     colour.red = 0.1;
                     colour.green = 0.6;
                     colour.blue = 0.6;
                  } else {
                     if (c == "yellow") {
                        colour.red = 0.8;
                        colour.green = 0.8;
                        colour.blue = 0.0;
                     } else {
                        if (c == "yellowtint") {
                           colour.red = 0.65;
                           colour.green = 0.65;
                           colour.blue = 0.4;
                        } else {
                           if (c == "orange") {
                              colour.red = 0.9;
                              colour.green = 0.6;
                              colour.blue = 0.1;
                           } else {
                              if (c == "red") {
                                 colour.red = 0.9;
                                 colour.green = 0.1;
                                 colour.blue = 0.1;
                              } else {
                                 if (c == "hotpink") {
                                    colour.red = 0.9;
                                    colour.green = 0.2;
                                    colour.blue = 0.6;
                                 } else {
                                    if (c == "pink") {
                                       colour.red = 0.9;
                                       colour.green = 0.3;
                                       colour.blue = 0.3;
                                    } else {
                                       if (c == "cyan") {
                                          colour.red = 0.1;
                                          colour.green = 0.7;
                                          colour.blue = 0.7;
                                       } else {
                                          if (c == "aquamarine") {
                                             colour.red = 0.1;
                                             colour.green = 0.8;
                                             colour.blue = 0.6;
                                          } else {
                                             if (c == "forestgreen") {
                                                colour.red   = 0.6;
                                                colour.green = 0.8;
                                                colour.blue  = 0.1;
                                             } else {
                                                if (c == "yellowgreen") {
                                                   colour.red   = 0.6;
                                                   colour.green = 0.8;
                                                   colour.blue  = 0.2;
                                                } else {
                                                   if (c == "goldenrod") {
                                                      colour.red   = 0.85;
                                                      colour.green = 0.65;
                                                      colour.blue  = 0.12;
                                                   } else {
                                                      if (c == "orangered") {
                                                         colour.red   = 0.9;
                                                         colour.green = 0.27;
                                                         colour.blue  = 0.0;
                                                      } else {
                                                         if (c == "magenta") {
                                                            colour.red   = 0.7;
                                                            colour.green = 0.2;
                                                            colour.blue  = 0.7;
                                                         } else {
                                                            if (c == "cornflower") {
                                                               colour.red   = 0.38;
                                                               colour.green = 0.58;
                                                               colour.blue  = 0.93;
                                                            } else {
                                                               if (c == "royalblue") {
                                                                  colour.red   = 0.25;
                                                                  colour.green = 0.41;
                                                                  colour.blue  = 0.88;
                                                               }
                                                            }
                                                         }
                                                      }
                                                   }
                                                }
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

//    std::cout << "debug:: in colour_values_from_colour_name from colour " << c
// 	     << " we assign colour values "
// 	     << colour[0] << " "
// 	     << colour[1] << " "
// 	     << colour[2] << "\n";
   return colour;
}

