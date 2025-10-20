/*
 * src/meshed-generic-display-object.cc
 *
 * Copyright 2020 by Medical Research Council
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

#include <mmdb2/mmdb_manager.h>
#include "Python.h"

// Perhaps this (needed to compile meshed-generic-display-object.hh) should go in Mesh.hh.
#include <epoxy/gl.h>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>  // to_string()
#include <glm/gtx/rotate_vector.hpp>

#include "graphics-info.h"
#include "meshed-generic-display-object.hh"
#include "coot-utils/oct.hh"
#include "coot-utils/cylinder.hh"

glm::vec3 coord_orth_to_glm(const clipper::Coord_orth &co) {
   return glm::vec3(co.x(), co.y(), co.z());
}


glm::vec4
colour_holder_to_glm(const coot::colour_holder &ch) {
   return glm::vec4(ch.red, ch.green, ch.blue, ch.alpha);
}

void
meshed_generic_display_object::add_line(const coot::colour_holder &colour,
                                        const std::string &colour_name, int line_width,
                                        const std::pair<clipper::Coord_orth, clipper::Coord_orth> &coords) {

   auto vnc_vertex_to_generic_vertex = [] (const coot::api::vnc_vertex &v) {
      return s_generic_vertex(v.pos, v.normal, v.color);
   };

   auto vnc_vertex_vector_to_generic_vertex_vector = [vnc_vertex_to_generic_vertex] (const std::vector<coot::api::vnc_vertex> &vv) {
      std::vector<s_generic_vertex> vo(vv.size());
      for (unsigned int i=0; i<vv.size(); i++)
         vo[i] = vnc_vertex_to_generic_vertex(vv[i]);
      return vo;
   };

   glm::vec3 start = coord_orth_to_glm(coords.first);
   glm::vec3 end   = coord_orth_to_glm(coords.second);
   auto cart_pair = std::make_pair(start, end);
   float bl = glm::distance(start, end);
   auto col = colour_holder_to_glm(colour);
   cylinder c(cart_pair, line_width, line_width, bl, col);
   c.add_flat_start_cap();
   c.add_flat_end_cap();
   std::vector<s_generic_vertex> converted_vertices = vnc_vertex_vector_to_generic_vertex_vector(c.vertices);
   mesh.import(converted_vertices, c.triangles);
}

void
meshed_generic_display_object::add_lines(std::vector<line_info_t> &liv) {

   auto vnc_vertex_to_generic_vertex = [] (const coot::api::vnc_vertex &v) {
      return s_generic_vertex(v.pos, v.normal, v.color);
   };

   auto vnc_vertex_vector_to_generic_vertex_vector = [vnc_vertex_to_generic_vertex] (const std::vector<coot::api::vnc_vertex> &vv) {
      std::vector<s_generic_vertex> vo(vv.size());
      for (unsigned int i=0; i<vv.size(); i++)
         vo[i] = vnc_vertex_to_generic_vertex(vv[i]);
      return vo;
   };

   for (unsigned int i=0; i<liv.size(); i++) {
      const auto &li = liv[i];
      glm::vec4 col = colour_holder_to_glm(li.colour);
      glm::vec3 start = coord_orth_to_glm(li.position_start);
      glm::vec3 end   = coord_orth_to_glm(li.position_end);
      float bl = glm::distance(start, end);
      auto cart_pair = std::make_pair(start, end);
      cylinder c(cart_pair, li.radius, li.radius, bl, col);
      c.add_flat_start_cap();
      c.add_flat_end_cap();
      std::vector<s_generic_vertex> converted_vertices = vnc_vertex_vector_to_generic_vertex_vector(c.vertices);
      mesh.import(converted_vertices, c.triangles);
   }

}


void
meshed_generic_display_object::add_sphere(const meshed_generic_display_object::sphere_t &sphere) {

   unsigned int num_subdivisions = 2;
   float radius = sphere.radius;
   glm::vec4 col = sphere.col;
   glm::vec3 position = coord_orth_to_glm(sphere.centre);
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
      oct = wrapped_make_octasphere(num_subdivisions, position, radius, col);
   mesh.import(oct);
}

void
meshed_generic_display_object::add_point(const coot::colour_holder &colour_in,
                                         const std::string &colour_name,
                                         const int &size_in,
                                         const clipper::Coord_orth &coords_in,
                                         unsigned int num_subdivisions) {

   // unsigned int num_subdivisions = 1;
   float radius = 0.03 * size_in; // changing the scaling is fun
   glm::vec4 col(colour_in.red, colour_in.green, colour_in.blue, colour_in.alpha);
   object_info_t oi; // 20240414-PE I need to add an oi for other object types too.
   oi.position = coords_in;
   oi.colour = colour_in;
   info.push_back(oi);
   glm::vec3 position_glm = coord_orth_to_glm(coords_in);
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > oct =
       wrapped_make_octasphere(num_subdivisions, position_glm, radius, col);
   if (false)
      std::cout << "debug:: mgdo::add_point() with colour " << glm::to_string(col)
                << " adding " << oct.first.size() << " " << oct.second.size()
                << " vertices and triangles " << std::endl;
   mesh.import(oct);

}

void
meshed_generic_display_object::add_points(std::vector<point_info_t> &piv, unsigned int num_subdivisions) {

   for (unsigned int i=0; i<piv.size(); i++) {
      const auto &pi = piv[i];
      glm::vec3 position_glm = coord_orth_to_glm(pi.position);
      float radius = 0.03 * static_cast<float>(pi.width);
      glm::vec4 col = colour_holder_to_glm(pi.colour);
      // std::cout << "in add_points() with colour " << glm::to_string(col) << std::endl;
      std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > oct =
         wrapped_make_octasphere(num_subdivisions, position_glm, radius, col);
      mesh.import(oct);
   }

}

void
meshed_generic_display_object::add_arrow(const arrow_t &arrow) {

   // unsigned int n_slices = 8;
   // std::pair<glm::vec3, glm::vec3 > start_end(coord_orth_to_glm(arrow.start_point),
   //                                            coord_orth_to_glm(arrow.end_point));
   // add_cylinder(start_end, arrow.col, arrow.radius, n_slices, true, true, FLAT_CAP, FLAT_CAP);

   coot::colour_holder col = arrow.col;

   float h = glm::distance(coord_orth_to_glm(arrow.start_point), coord_orth_to_glm(arrow.end_point));
   glm::vec4 base_colour(arrow.col.red, arrow.col.green, arrow.col.blue, 1.0f);
   unsigned int n_slices = 30;
   std::pair<glm::vec3, glm::vec3> start_end(coord_orth_to_glm(arrow.start_point), coord_orth_to_glm(arrow.end_point));
   cylinder c(start_end, arrow.radius, arrow.radius, h, base_colour, n_slices, 2);
   c.add_flat_start_cap();
   add_cylinder(start_end, col, arrow.radius, n_slices, true, false, FLAT_CAP, FLAT_CAP);

   clipper::Coord_orth delta_uv((arrow.end_point - arrow.start_point).unit());
   clipper::Coord_orth cone_start = arrow.end_point + 1.3 * delta_uv;
   clipper::Coord_orth cone_end   = arrow.end_point;
   float base_radius = arrow.radius * 3.0;
   float top_radius = 0.0;
   std::pair<glm::vec3, glm::vec3> cone_start_end(coord_orth_to_glm(cone_start), coord_orth_to_glm(cone_end));
   add_cone(cone_start_end, col, base_radius, top_radius, n_slices, false, true, FLAT_CAP, FLAT_CAP);

}


std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
meshed_generic_display_object::wrapped_make_octasphere(unsigned int num_subdivisions,
                                                       const glm::vec3 &position,
                                                       float radius,
                                                       const glm::vec4 &col) {

   auto vnc_vertex_to_generic_vertex = [] (const coot::api::vnc_vertex &v) {
      return s_generic_vertex(v.pos, v.normal, v.color);
   };

   auto vnc_vertex_vector_to_generic_vertex_vector = [vnc_vertex_to_generic_vertex] (const std::vector<coot::api::vnc_vertex> &vv) {
      std::vector<s_generic_vertex> vo(vv.size());
      for (unsigned int i=0; i<vv.size(); i++)
         vo[i] = vnc_vertex_to_generic_vertex(vv[i]);
      return vo;
   };

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
      std::pair<std::vector<coot::api::vnc_vertex>, std::vector<g_triangle> > oct = make_octasphere(num_subdivisions, origin, radius, col);
      std::vector<s_generic_vertex> converted_vertices = vnc_vertex_vector_to_generic_vertex_vector(oct.first);
      std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > oct_c(converted_vertices, oct.second);

      origin_octasphere_map[map_index] = oct_c;
      for (unsigned int i=0; i<oct.first.size(); i++) {
         s_generic_vertex &v(oct_c.first[i]);
         v.pos += position;
         v.color = col;
      }
      return oct_c;
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
            float line_radius = 0.06;
            if (do_thinning) line_radius *= 0.5;
            add_cylinder(start_end, col, line_radius, n_slices, true, true, FLAT_CAP, FLAT_CAP);
         }
      }
   }
}


void
meshed_generic_display_object::add_cylinder(const std::pair<glm::vec3, glm::vec3 > &start_end,
                                            const coot::colour_holder &col, float line_radius,
                                            unsigned int n_slices,
                                            bool cap_start, bool cap_end,
                                            cap_type start_cap_type, cap_type end_cap_type,
                                            bool do_faces,
                                            float unstubby_cap_factor) {

   auto vnc_vertex_to_generic_vertex = [] (const coot::api::vnc_vertex &v) {
      return s_generic_vertex(v.pos, v.normal, v.color);
   };

   auto vnc_vertex_vector_to_generic_vertex_vector = [vnc_vertex_to_generic_vertex] (const std::vector<coot::api::vnc_vertex> &vv) {
      std::vector<s_generic_vertex> vo(vv.size());
      for (unsigned int i=0; i<vv.size(); i++)
         vo[i] = vnc_vertex_to_generic_vertex(vv[i]);
      return vo;
   };

   if (false) {
      std::cout << "add_cylinder(): cap_start: " << cap_start  << std::endl;
      std::cout << "add_cylinder(): cap_end: "   << cap_end  << std::endl;
      std::cout << "add_cylinder(): do_faces: "  << do_faces  << std::endl;
      std::cout << "add_cylinder(): unstubby_cap_factor: " << unstubby_cap_factor  << std::endl;
   }

   float h = glm::distance(start_end.first, start_end.second);
   glm::vec4 base_colour(col.red, col.green, col.blue, 1.0f);
   if (do_faces) cap_start = false;
   cylinder c(start_end, line_radius, line_radius, h, base_colour, n_slices, 2); // not colour of base
   if (do_faces) c.crenulations();
   c.set_unstubby_rounded_cap_factor(unstubby_cap_factor);
   if (false)
      std::cout << "add_cylinder: " << glm::to_string(start_end.first) << " "
                << glm::to_string(start_end.second) << " "
                << c.vertices.size() << " " << c.triangles.size()
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
   if (do_faces) {
      c.add_sad_face();
   }

   // for (unsigned int i=0; i<c.vertices.size(); i++)
   // c.vertices[i].color = colour;

   std::vector<s_generic_vertex> converted_vertices = vnc_vertex_vector_to_generic_vertex_vector(c.vertices);
   mesh.import(converted_vertices, c.triangles);

   // caches eyelashes!

   if (do_faces) {
      // std::string file_name("grey-eyelashes-many-lashes.glb");
      // eyelash_r.load_from_glTF(file_name, false); // tries local directory first
      Mesh eyelash_r = graphics_info_t::get_mesh_for_eyelashes();
      eyelash_r.apply_scale(0.026);
      eyelash_r.translate_by(glm::vec3(0.07, 0, 0.93));
      Mesh eyelash_l = eyelash_r;
      glm::mat4 rm(1.0f);
      glm::mat4 mirror_y(1.0f); mirror_y[1][1] = -1.0f; // does this change the normals also?
      eyelash_l.apply_transformation(mirror_y);
      eyelash_l.invert_normals();
      glm::mat4 m_r = glm::rotate(rm,  0.5f, glm::vec3(0,0,1));
      glm::mat4 m_l = glm::rotate(rm, -0.5f, glm::vec3(0,0,1));
      eyelash_r.apply_transformation(m_r);
      eyelash_l.apply_transformation(m_l);
      mesh.import(eyelash_r.vertices, eyelash_r.triangles);
      mesh.import(eyelash_l.vertices, eyelash_l.triangles);
   }

}

void
meshed_generic_display_object::add_cone(const std::pair<glm::vec3, glm::vec3> &start_end,
                                        const coot::colour_holder &col, float base_radius, float top_radius,
                                        unsigned int n_slices,
                                        bool cap_start, bool cap_end,
                                        cap_type start_cap_type, cap_type end_cap_type) {

   auto vnc_vertex_to_generic_vertex = [] (const coot::api::vnc_vertex &v) {
      return s_generic_vertex(v.pos, v.normal, v.color);
   };

   auto vnc_vertex_vector_to_generic_vertex_vector = [vnc_vertex_to_generic_vertex] (const std::vector<coot::api::vnc_vertex> &vv) {
      std::vector<s_generic_vertex> vo(vv.size());
      for (unsigned int i=0; i<vv.size(); i++)
         vo[i] = vnc_vertex_to_generic_vertex(vv[i]);
      return vo;
   };

   float h = glm::distance(start_end.first, start_end.second);
   cylinder c(start_end, base_radius, top_radius, h, n_slices, 2);
   glm::vec4 colour(col.red, col.green, col.blue, 1.0f);

   // are your start and end points where you think they are?
   if (false)
      std::cout << "add_cone: " << glm::to_string(start_end.first) << " "
                << glm::to_string(start_end.second)
                << " base_radius " << base_radius << " top_radius " << top_radius << " "
                << c.vertices.size() << " " << c.triangles.size()
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

   std::vector<s_generic_vertex> converted_vertices = vnc_vertex_vector_to_generic_vertex_vector(c.vertices);
   mesh.import(converted_vertices, c.triangles);

}



void
meshed_generic_display_object::add_dodecahedron(const coot::colour_holder &colour_in,
                                                const std::string &colour_name,
                                                double radius, const clipper::Coord_orth &pos) {
}

#include "make-a-dodec.hh"

void
meshed_generic_display_object::add_pentakis_dodecahedron(const coot::colour_holder &colour_in,
                                                         const std::string &colour_name,
                                                         double stellation_factor,
                                                         double radius,
                                                         const clipper::Coord_orth &pos_c) {

   // these lambdas could/should be in the class

   auto clipper_to_glm = [] (const clipper::Coord_orth &c) {
      return glm::vec3(c.x(), c.y(), c.z());
   };

   auto vn_vertex_to_generic_vertex = [] (const vn_vertex &v,
					  const glm::vec3 &pos,
					  float radius,
					  const glm::vec4 &color) {
      return s_generic_vertex(v.pos * radius + pos, v.normal, color);
   };

   auto vn_vertex_vector_to_generic_vertex_vector = [vn_vertex_to_generic_vertex]
      (const std::vector<vn_vertex> &vv,
       const glm::vec3 &pos,
       float radius,
       const glm::vec4 &color) {

      std::vector<s_generic_vertex> vo(vv.size());
      for (unsigned int i=0; i<vv.size(); i++)
         vo[i] = vn_vertex_to_generic_vertex(vv[i], pos, radius, color);
      return vo;
   };

   int spikey_mode = DODEC_SPIKEY_MODE;
   std::pair<std::vector<vn_vertex>, std::vector<g_triangle> > dodec =
      make_pentakis_dodec(spikey_mode);
   const std::vector<vn_vertex>   &vertices = dodec.first;
   const std::vector<g_triangle> &triangles = dodec.second;

   glm::vec4 col = colour_holder_to_glm(colour_in);
   glm::vec3 pos = clipper_to_glm(pos_c);

   std::vector<s_generic_vertex> converted_vertices =
      vn_vertex_vector_to_generic_vertex_vector(vertices, pos, radius, col);
   mesh.import(converted_vertices, triangles);

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

   const unsigned int n_phi_steps = 60;
   const unsigned int n_theta_steps = 60;
   std::vector<s_generic_vertex> vertices((n_theta_steps + 1) * n_phi_steps); // 20211007-PE this doesn't look right!
                                                                              // c.f. add_dashed_angle_markup()
   std::vector<g_triangle> triangles;
   const float R = arc.radius;
   const float r = arc.radius_inner;
   const float pi = 3.1415926535;

   float angle_delta_deg = arc.delta_angle;
   float angle_delta = angle_delta_deg * pi / 180.0f;
   glm::vec3 start_point = coord_orth_to_glm(arc.start_point); // this is the central atom position
   glm::vec4 col = colour_holder_to_glm(arc.col);
   glm::mat4 rot_mat = glm::orientation(coord_orth_to_glm(arc.normal), glm::vec3(0.0, 0.0, 1.0));
   // std::cout << "add_arc: rot_mat: " << glm::to_string(rot_mat)<< std::endl;

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
         v.pos += start_point; // central atom position
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

   // std::cout << "torus.normal: " << torus.normal.format() << std::endl;

   glm::vec4 centre(torus.position.x(), torus.position.y(), torus.position.z(), 1.0);
   glm::vec3 ring_normal(glm::normalize(glm::vec3(torus.normal.x(), torus.normal.y(), torus.normal.z())));
   glm::mat4 ori = glm::orientation(ring_normal, glm::vec3(0.0, 0.0, 1.0));

   // std::cout << "ring_normal " << glm::to_string(ring_normal) << std::endl;
   // std::cout << "ori " << glm::to_string(ori) << std::endl;

   for (unsigned int ip=0; ip<n_phi_steps; ip++) {
      float phi = 2.0f * pi * static_cast<float>(ip)/static_cast<float>(n_phi_steps);
      for (unsigned int it=0; it<n_theta_steps; it++) {
         float theta = 2.0f * pi * static_cast<float>(it)/static_cast<float>(n_theta_steps);
         s_generic_vertex v;
         glm::vec4 pos;
         pos.x = (R + r * cosf(theta)) * cosf(phi);
         pos.y = (R + r * cosf(theta)) * sinf(phi);
         pos.z = torus.height_scale * r * sinf(theta);
         pos.w = 1.0; // or 0?
         glm::mat4 tori = glm::transpose(ori);
         v.pos = glm::vec3(pos * tori) + glm::vec3(centre);
         glm::vec4 normal;
         normal.x = cosf(theta) * cosf(phi);
         normal.y = cosf(theta) * sinf(phi);
         normal.z = sinf(theta);
         normal.w = 1.0;
         v.normal = glm::vec3(normal * tori);
         v.color = col;
         vertices[ip * n_theta_steps + it] = v;
      }
   }

   // carefully, carefully :-)
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

void
meshed_generic_display_object::add_dashed_line(const coot::colour_holder &colour, const std::string &colour_name,
                                               const std::pair<clipper::Coord_orth, clipper::Coord_orth> &coords,
                                               const Material &material,
                                               float line_width_scale, unsigned int n_segments) // default args 1.0, 5
{

   coot::simple_distance_object_t sdo(-1, coords.first, -1, coords.second);
   glm::vec4 glm_colour(colour.red, colour.green, colour.blue, 1.0);
   mesh.add_dashed_line(sdo, material, glm_colour);
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
                  colour.red   = 0.48;
                  colour.green = 0.05;
                  colour.blue  = 0.5;
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

void
meshed_generic_display_object::translate(const coot::Cartesian &t) {

   glm::vec3 tt(t.x(), t.y(), t.z());
   mesh.translate_by(tt);
}


// remove from info vector and remove 182 triangles from the mesh (that's a bit of a hack)
void
meshed_generic_display_object::remove_last_object() {

   if (! info.empty()) {
      info.pop_back();
   }

   mesh.remove_last_subobject(74, 128);
}
