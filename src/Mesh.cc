/*
 * src/Mesh.cc
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

#include <vector>
#include <iostream>
#include <map>
#include <set>
#include <chrono>
#include "stereo-eye.hh"

#ifdef USE_BACKWARD
#include <utils/backward.hpp>
#endif


// Having this up here...
#include <mmdb2/mmdb_manager.h>
// Fixes this:
// In file included from /usr/include/X11/Xlib.h:44,
//                 from /usr/include/epoxy/glx.h:36,
//                 from ../../coot/src/Mesh.cc:9:
///home/paule/autobuild/Linux-penelope-pre-release-gtk3-python/include/mmdb2/mmdb_io_file.h:131:22: error: expected unqualified-id before numeric constant
//  131 |         inline bool  Success   () { return IOSuccess; }
//      |                      ^~~~~~~


#include <epoxy/gl.h>

// #define THIS_IS_HMT

#include "Mesh.hh"
#include "Shader.hh"
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>  // to_string()
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/vector_angle.hpp>
#include <glm/gtc/type_ptr.hpp>  // for value_ptr() 20240326-PE

#include "coot-utils/oct.hh"
#include "coot-utils/cylinder.hh"

#ifdef THIS_IS_HMT

// Get rid of this at some stage - needs thinking about meshed objects
#include "generic-display-object.hh"
#else
#include "old-generic-display-object.hh"
#endif

#include "utils/logging.hh"
extern logging logger;


// this can be outside of Mesh
std::string stringify_error_code(GLenum err) {

   std::string r = std::to_string(err);
   if (err == GL_INVALID_ENUM)  r = "GL_INVALID_ENUM";
   if (err == GL_INVALID_VALUE) r = "GL_INVALID_VALUE";
   if (err == GL_INVALID_OPERATION) r = "GL_INVALID_OPERATION";
   return r;
}

void
Mesh::init() {

   vao = VAO_NOT_SET;
   clear(false);
   first_time = true;
   is_instanced = false;
   is_instanced_colours = false;
   is_instanced_with_rts_matrix = false;
   use_blending = false;
   draw_this_mesh = true;
   hydrogen_bond_cylinders_angle = 0.0;
   normals_are_setup = false;
   this_mesh_is_closed = false;
   n_instances = 0;
   n_instances_allocated = 0;
   particle_draw_count = 0;
   gl_lines_mode = false;
   is_headless = false;
   mesh_is_semi_transparent = false;
   buffer_id = 0; // not valid
   index_buffer_id = 0; // not valid
   debug_mode = false;
   inst_colour_buffer_id = -1;
   inst_model_translation_buffer_id = -1;

   std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
   time_constructed = now;
}

Mesh::Mesh(const std::pair<std::vector<s_generic_vertex>, const std::vector<g_triangle> > &indexed_vertices) {

   init();
   vertices  = indexed_vertices.first;
   triangles = indexed_vertices.second;
}

Mesh::Mesh(const std::vector<s_generic_vertex> &vertices_in, const std::vector<g_triangle> &triangles_in) {

   init();
   vertices  = vertices_in;
   triangles = triangles_in;
}


// a molecular_triangles_mesh_t is a poor man's Mesh. Why does it exist?
Mesh::Mesh(const molecular_triangles_mesh_t &mtm) {

   init();
   vertices  = mtm.vertices;
   triangles = mtm.triangles;
   name = mtm.name;
}


Mesh::Mesh(const std::string &name_in, const coot::simple_mesh_t &mesh) {

   name = name_in;
   vertices.resize(mesh.vertices.size());
   for (unsigned int i = 0; i < mesh.vertices.size(); i++) {
      const auto &vv = mesh.vertices[i];
      s_generic_vertex v(vv.pos, vv.normal, vv.color);
      vertices[i] = v;
   }
   triangles = mesh.triangles;
}

Mesh::~Mesh() {

   if (false)
      std::cout << "DEBUG:: Mesh destructor " << name << " with vao " << vao << std::endl;

   bool gl_buffers_flag = false;
   clear(gl_buffers_flag);
}



void
Mesh::set_is_headless() {
   is_headless = true;
}


// static 
std::pair<float, float> Mesh::get_stereo_x_scale_and_offset(stereo_eye_t eye) {

   // for side by side stereo, of course
   float stereo_x_scale  = 1.0;
   float stereo_x_offset = 0.0;
   if (eye == stereo_eye_t::LEFT_EYE) {
      stereo_x_scale = 2.0f;
   }
   if (eye == stereo_eye_t::RIGHT_EYE) {
      stereo_x_offset = -0.5f;
   }
   return std::pair<float, float>(stereo_x_scale, stereo_x_offset);
}

void
Mesh::set_draw_this_mesh(bool state) {

   if (state)
      if (! vertices.empty())
         if (! triangles.empty())
            draw_this_mesh = true;

   if (! state)
      draw_this_mesh = false;
}

void
Mesh::close() {

   clear();
   draw_this_mesh = false;
   this_mesh_is_closed = true;
}

void
Mesh::import(const std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > &indexed_vertices,
             bool fill_lines_vertex_indices) { // default false

   is_instanced = false;
   is_instanced_colours = false;
   is_instanced_with_rts_matrix = false;
   use_blending = false;
   if (fill_lines_vertex_indices) gl_lines_mode =  true; // flag used for drawing.

   unsigned int idx_base = vertices.size();
   unsigned int idx_tri_base = triangles.size();
   vertices.insert(vertices.end(), indexed_vertices.first.begin(), indexed_vertices.first.end());
   triangles.insert(triangles.end(),
                    indexed_vertices.second.begin(),
                    indexed_vertices.second.end());
   for (unsigned int i=idx_tri_base; i<triangles.size(); i++)
      triangles[i].rebase(idx_base);

   if (fill_lines_vertex_indices) {
      // every triangle has 6 lines
      const auto &triangles = indexed_vertices.second;
      lines_vertex_indices.resize(triangles.size() * 6);
      // replace, not add - is that a good idea?
      for (unsigned int i=0; i<triangles.size(); i++) {
         lines_vertex_indices[i*6  ] = triangles[i].point_id[0];
         lines_vertex_indices[i*6+1] = triangles[i].point_id[1];
         lines_vertex_indices[i*6+2] = triangles[i].point_id[1];
         lines_vertex_indices[i*6+3] = triangles[i].point_id[2];
         lines_vertex_indices[i*6+4] = triangles[i].point_id[2];
         lines_vertex_indices[i*6+5] = triangles[i].point_id[0];
      }
   }
}

void
Mesh::import(const std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > &indexed_vertices,
             const std::vector<std::pair<int, map_triangle_t> > &map_triangle_centres_in,
             bool fill_lines_vertex_indices) { // default false

   import(indexed_vertices, fill_lines_vertex_indices);
   map_triangle_centres = map_triangle_centres_in;
}


void
Mesh::import(const std::vector<s_generic_vertex> &gv, const std::vector<g_triangle> &indexed_vertices) {

   is_instanced = false;
   is_instanced_colours = false;
   is_instanced_with_rts_matrix = false;

   unsigned int idx_base = vertices.size();
   unsigned int idx_tri_base = triangles.size();

   // std::cout << "on Mesh::import idx_base " << idx_base << " idx_tri_base " << idx_tri_base << std::endl;
   vertices.insert(vertices.end(), gv.begin(), gv.end());
   triangles.insert(triangles.end(),
                    indexed_vertices.begin(),
                    indexed_vertices.end());
   if (idx_base > 0)
      for (unsigned int i=idx_tri_base; i<triangles.size(); i++)
         triangles[i].rebase(idx_base);
}

void
Mesh::import(const std::vector<position_normal_vertex> &verts, const std::vector<g_triangle> &indexed_vertices,
             const glm::vec4 &colour) {

   is_instanced = false;
   is_instanced_with_rts_matrix = false;

   unsigned int idx_base = vertices.size();
   unsigned int idx_tri_base = triangles.size();

   std::vector<s_generic_vertex> gv(verts.size());
   for (unsigned int ii=0; ii<verts.size(); ii++) {
      gv[ii].pos    = verts[ii].pos;
      gv[ii].normal = verts[ii].normal;
      gv[ii].color  = colour;
   }

   vertices.insert(vertices.end(), gv.begin(), gv.end());
   triangles.insert(triangles.end(),
                    indexed_vertices.begin(),
                    indexed_vertices.end());
   if (idx_base != 0)
      for (unsigned int i=idx_tri_base; i<triangles.size(); i++)
         triangles[i].rebase(idx_base);

   // now the caller function should call setup(shader, material) so that setup_buffers() is run

}

void
Mesh::import_lines(const std::vector<s_generic_vertex> &v,
                   const std::vector<unsigned int> &line_indices) {

   vertices = v;

   std::cout << "::::::::::::::::::: import_lines vertices.size " << vertices.size() << std::endl;
   std::cout << "::::::::::::::::::: import_lines lines_indices .size " << line_indices.size() << std::endl;
   lines_vertex_indices = line_indices;
   gl_lines_mode = true;
}


void
Mesh::apply_scale(float s) {

   glm::vec3 sf(s, s, s);
   for (unsigned int ii=0; ii<vertices.size(); ii++)
      vertices[ii].pos *= s;
   setup_buffers(); // transfer the coordinates
}

void
Mesh::apply_transformation(const glm::mat4 &m) {

   for (unsigned int ii=0; ii<vertices.size(); ii++) {
      glm::vec4 p_1(vertices[ii].pos, 1.0);
      glm::vec4 p_2 = p_1 * m;
      glm::vec3 p_3(p_2);
      vertices[ii].pos = p_3;
   }
   setup_buffers(); // transfer the coordinates
}

void
Mesh::translate_by(const glm::vec3 &t) {

   for (unsigned int ii=0; ii<vertices.size(); ii++)
      vertices[ii].pos += t;

   setup_buffers();
}

// this does brightening too. Pass a number between 0 and 1 - and that will change the ammount that each of the
// colour components becomes full. So degree 0 is no change and degree 1 is full red, green, blue (white, obviously).
void
Mesh::pastelize(float degree) {

   for (unsigned int ii=0; ii<vertices.size(); ii++) {
      for (unsigned int i=0; i<3; i++) {
         const float &cc = vertices[ii].color[i];
         float r = 1.0f - cc;
         vertices[ii].color[i] += r * degree;
      }
   }
   setup_buffers();
}

void
Mesh::brighten(float scale_factor) {

   for (unsigned int ii=0; ii<vertices.size(); ii++) {
      for (unsigned int i=0; i<3; i++) {
         const float &cc = vertices[ii].color[i];
         vertices[ii].color[i] *= scale_factor;
      }
   }
   setup_buffers();
}




void
Mesh::debug() const {

   std::cout << "Mesh::debug() " << name << " n-vertices " << vertices.size()
             << " n-triangles " << triangles.size() << std::endl;
}


// 20210910-PE remove the Shader argument, surely not needed.
void
Mesh::setup(const Material &material_in) {

   // 20250224-PE I don't like the idea of mixing up material
   // and buffers. Rename this to set_material()0
   // and don't setup_buffers().

   // Generic objects need a lot of reworking
   // because I want to "setup" only after all objects have been added - not
   // after adding every object (say 1000 spheres/dots).
   // But I don't want to draw if the object is not setup. Hmm.
   // Also adding thousands of balls is slow!
   // HOLE dots is a good example of where I should be using instancing

   // if (setup_has_been_done)
   // return;

   material = material_in;
   // shader_p->Use();
   setup_buffers();
}


// maybe should not be const
void
Mesh::setup_simple_triangles(Shader *shader_p, const Material &material_in) {

   material = material_in;
   shader_p->Use();
   fill_with_simple_triangles_vertices();
   fill_with_direction_triangles();
   setup_buffers();

}

void
Mesh::setup_rama_balls(Shader *shader_p, const Material &material_in) {

   material = material_in;
   shader_p->Use();
   fill_rama_balls();
   setup_buffers();
}


void
Mesh::fill_one_dodec() {

   // These are dodec balls for the moment, not flat shaded real dodecs.
   // In a bit I will transform them to pentakis dodecs.

   dodec d;
   std::vector<clipper::Coord_orth> v = d.coords();
   vertices.resize(v.size()); // should be 20

   for (unsigned int i=0; i<v.size(); i++) {
      glm::vec3 a = glm::vec3(v[i].x(), v[i].y(), v[i].z());
      float f = static_cast<float>(i)/static_cast<float>(v.size());
      vertices[i].pos    = 0.2f * a;
      vertices[i].normal = a;
      vertices[i].color  = glm::vec4(0.3 + 0.4 * f, 0.8 - 0.8 * f, 0.1 + 0.9 * f, 1.0);
   }

   for (unsigned int i=0; i<12; i++) {
      const std::vector<unsigned int> &face = d.face(i);
      g_triangle gt_0(face[0], face[1], face[2]);
      g_triangle gt_1(face[0], face[2], face[3]);
      g_triangle gt_2(face[0], face[3], face[4]);
      triangles.push_back(gt_0);
      triangles.push_back(gt_1);
      triangles.push_back(gt_2);
   }

}

#ifdef THIS_IS_HMT
#else
#include "old-generic-display-object.hh"
#endif


void
Mesh::add_one_ball(float scale, const glm::vec3 &centre) { // i.e. a smooth-shaded pentakis dodec

   pentakis_dodec pkdd(1.05);
   // coordninate system finagalling - baah - convert dodecs to glm.
   clipper::Coord_orth centre_c(centre.x, centre.y, centre.z);
   coot::old_generic_display_object_t::pentakis_dodec_t penta_dodec(pkdd, 0.01, centre_c);
   std::vector<clipper::Coord_orth> v = penta_dodec.pkdd.d.coords();

   unsigned int vertex_index_start_base   = vertices.size();

   const std::vector<clipper::Coord_orth> &pv = penta_dodec.pkdd.pyrimid_vertices;
   for (unsigned int i=0; i<12; i++) {

      const std::vector<unsigned int> &face = penta_dodec.pkdd.d.face(i);

      // first the base point (tip of the triangles/pyrimid)
      //
      clipper::Coord_orth pvu(pv[i].unit());
      s_generic_vertex gv;
      gv.pos    = scale * (glm::vec3(pv[i].x(), pv[i].y(), pv[i].z())) + centre;
      gv.normal = glm::vec3(pvu.x(), pvu.y(), pvu.z());
      vertices.push_back(gv);

      for (unsigned int j=0; j<5; j++) {
         const clipper::Coord_orth &pt = v[face[j]];
         clipper::Coord_orth ptu(pt.unit());
         gv.pos    = scale * (glm::vec3(pt.x(),  pt.y(),  pt.z())) + centre;
         gv.normal = glm::vec3(ptu.x(), ptu.y(), ptu.z());
         vertices.push_back(gv);
      }
      unsigned int idx_base = vertex_index_start_base + i*6;
      g_triangle gt_0(idx_base, idx_base + 1, idx_base + 2);
      g_triangle gt_1(idx_base, idx_base + 2, idx_base + 3);
      g_triangle gt_2(idx_base, idx_base + 3, idx_base + 4);
      g_triangle gt_3(idx_base, idx_base + 4, idx_base + 5);
      g_triangle gt_4(idx_base, idx_base + 5, idx_base + 1);
      triangles.push_back(gt_0);
      triangles.push_back(gt_1);
      triangles.push_back(gt_2);
      triangles.push_back(gt_3);
      triangles.push_back(gt_4);
   }
}

void
Mesh::add_one_origin_ball() { // i.e. a smooth-shaded pentakis dodec

   pentakis_dodec pkdd(1.05);
   // coordninate system finagalling - baah - convert dodecs to glm.
   clipper::Coord_orth centre_c(0,0,0);
   coot::old_generic_display_object_t::pentakis_dodec_t penta_dodec(pkdd, 0.01, centre_c);
   std::vector<clipper::Coord_orth> v = penta_dodec.pkdd.d.coords();

   unsigned int vertex_index_start_base   = vertices.size();

   const std::vector<clipper::Coord_orth> &pv = penta_dodec.pkdd.pyrimid_vertices;
   for (unsigned int i=0; i<12; i++) {

      const std::vector<unsigned int> &face = penta_dodec.pkdd.d.face(i);
      float scale = 0.5;

      // first the base point (tip of the triangles/pyrimid)
      //
      clipper::Coord_orth pvu(pv[i].unit());
      s_generic_vertex gv;
      gv.pos    = scale * 0.1f * (glm::vec3(pv[i].x(), pv[i].y(), pv[i].z()));
      gv.normal = glm::vec3(pvu.x(), pvu.y(), pvu.z());
      vertices.push_back(gv);

      for (unsigned int j=0; j<5; j++) {
         const clipper::Coord_orth &pt = v[face[j]];
         clipper::Coord_orth ptu(pt.unit());
         gv.pos    = scale * 0.1f * (glm::vec3(pt.x(),  pt.y(),  pt.z()));
         gv.normal = glm::vec3(ptu.x(), ptu.y(), ptu.z());
         vertices.push_back(gv);
      }
      unsigned int idx_base = vertex_index_start_base + i*6;
      g_triangle gt_0(idx_base, idx_base + 1, idx_base + 2);
      g_triangle gt_1(idx_base, idx_base + 2, idx_base + 3);
      g_triangle gt_2(idx_base, idx_base + 3, idx_base + 4);
      g_triangle gt_3(idx_base, idx_base + 4, idx_base + 5);
      g_triangle gt_4(idx_base, idx_base + 5, idx_base + 1);
      triangles.push_back(gt_0);
      triangles.push_back(gt_1);
      triangles.push_back(gt_2);
      triangles.push_back(gt_3);
      triangles.push_back(gt_4);
   }
}


void
Mesh::add_one_origin_dodec() { // i.e. a smooth-shaded pentakis dodec

   if (false) {
      glm::vec3 p0(0,0,0);
      glm::vec3 p1(1,0,0);
      glm::vec3 p2(1,1,0);
      glm::vec3 n(0,0,0);
      glm::vec4 c(1,0,0,1);
      s_generic_vertex gv0(p0,n,c);
      s_generic_vertex gv1(p1,n,c);
      s_generic_vertex gv2(p2,n,c);
      vertices.push_back(gv0);
      vertices.push_back(gv1);
      vertices.push_back(gv2);
      triangles.push_back(g_triangle(0,1,2));
   }

   if (true) {
      float scale = 0.05;
      dodec d;
      float radius = 0.2;
      clipper::Coord_orth position(0,0,0);
      coot::old_generic_display_object_t::dodec_t dod(d, radius, position);
      // we can put a colour into dod here
   
      std::vector<clipper::Coord_orth> v = dod.d.coords();

      unsigned int vertex_index_start_base = vertices.size();
      for (unsigned int i=0; i<12; i++) {
         const std::vector<unsigned int> &face = dod.d.face(i);
         clipper::Coord_orth sum_vertex(0,0,0);
         for (unsigned int j=0; j<5; j++)
            sum_vertex += v[face[j]];
         clipper::Coord_orth face_normal(sum_vertex.unit());
         for (unsigned int j=0; j<5; j++) {
            glm::vec3 p(v[face[j]].x(), v[face[j]].y(), v[face[j]].z());
            glm::vec3 n(face_normal.x(), face_normal.y(), face_normal.z());
            glm::vec4 c(1,0,0,1);
            s_generic_vertex gv(scale * p,n,c);
            vertices.push_back(gv);
         }
         unsigned int idx_base = vertex_index_start_base + i*5;
         g_triangle gt_0(idx_base, idx_base + 1, idx_base + 2);
         g_triangle gt_1(idx_base, idx_base + 2, idx_base + 3);
         g_triangle gt_2(idx_base, idx_base + 3, idx_base + 4);
         triangles.push_back(gt_0);
         triangles.push_back(gt_1);
         triangles.push_back(gt_2);
      }
   }
}

void
Mesh::add_one_origin_cylinder(unsigned int n_slices, unsigned int n_stacks) {

   auto vnc_vertex_to_generic_vertex = [] (const coot::api::vnc_vertex &v) {
      return s_generic_vertex(v.pos, v.normal, v.color);
   };

   auto vnc_vertex_vector_to_generic_vertex_vector = [vnc_vertex_to_generic_vertex] (const std::vector<coot::api::vnc_vertex> &vv) {
      std::vector<s_generic_vertex> vo(vv.size());
      for (unsigned int i=0; i<vv.size(); i++)
         vo[i] = vnc_vertex_to_generic_vertex(vv[i]);
      return vo;
   };

   // short fat, radius 1, height 1.

   cylinder c(std::pair<glm::vec3, glm::vec3>(glm::vec3(0,0,0), glm::vec3(0,0,1)),
              1.0, 1.0, 1.0, n_slices, n_stacks);

   // vertices = std::move(c.vertices);
   // triangle_vertex_indices = std::move(c.triangle_indices_vec);

   unsigned int idx_base = vertices.size();
   unsigned int idx_tri_base = triangles.size();
   std::vector<s_generic_vertex> converted_vertices = vnc_vertex_vector_to_generic_vertex_vector(c.vertices);
   vertices.insert(vertices.end(), converted_vertices.begin(), converted_vertices.end());
   triangles.insert(triangles.end(), c.triangles.begin(), c.triangles.end());
   for (unsigned int i=idx_tri_base; i<triangles.size(); i++)
      triangles[i].rebase(idx_base);

}

void
Mesh::add_one_origin_octahemisphere(unsigned int num_subdivisions) {

   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > oct =
      tessellate_hemisphere_patch(num_subdivisions);

   float turn_fraction = 0.125; // default num_subdivisions = 1
   float angle = 2.0f * M_PI * turn_fraction;
   glm::vec3 z_axis(0,0,1);
   glm::vec4 atom_colour(0.2f, 0.6f, 0.3f, 1.0f);

   vertices.resize(oct.first.size());
   for (unsigned int i=0; i<oct.first.size(); i++) {
      vertices[i].pos    = glm::rotate(oct.first[i], angle, z_axis);
      vertices[i].normal = glm::rotate(oct.first[i], angle, z_axis);
      vertices[i].color  = atom_colour;
   }
   triangles = oct.second;

}

void
Mesh::add_one_origin_octasphere(unsigned int num_subdivisions) {

   // what is this?

   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > oct =
      tessellate_hemisphere_patch(num_subdivisions);

   float turn_fraction = 0.125; // default num_subdivisions = 1
   float angle = 2.0f * M_PI * turn_fraction;
   glm::vec3 z_axis(0,0,1);
   glm::vec4 atom_colour(0.2f, 0.6f, 0.3f, 1.0f);

   vertices.resize(oct.first.size());
   for (unsigned int i=0; i<oct.first.size(); i++) {
      vertices[i].pos    = glm::rotate(oct.first[i], angle, z_axis);
      vertices[i].normal = glm::rotate(oct.first[i], angle, z_axis);
      vertices[i].color  = atom_colour;
   }
   triangles = oct.second;

}

void
Mesh::fill_rama_balls() {

   // glm::vec3 centre(0,0,0.0);
   // add_one_ball(centre);

   for (unsigned int i=0; i<6; i++) {
      float f = static_cast<float>(i)/10.0;
      glm::vec3 c(0.44, 0.05, f - 0.2);
      float scale = 0.04;
      add_one_ball(scale, c);
   }
}

void
Mesh::setup_debugging_instancing_buffers() {

   // put these in their own file - make inst_matrices be member data

   is_instanced = true;
   is_instanced_colours = true;
   std::vector<glm::vec3> inst_trans_matrices;
   inst_trans_matrices.push_back(glm::vec3(0.25, 0.25, -0.2));
   inst_trans_matrices.push_back(glm::vec3(0.25, 0.25, -0.1));
   inst_trans_matrices.push_back(glm::vec3(0.25, 0.25,  0.0));
   inst_trans_matrices.push_back(glm::vec3(0.25, 0.25,  0.1));
   inst_trans_matrices.push_back(glm::vec3(0.25, 0.25,  0.2));
   inst_trans_matrices.push_back(glm::vec3(0.25, 0.25,  0.3));

   std::vector<glm::vec4> inst_col_matrices;
   inst_col_matrices.push_back(glm::vec4(0.8, 0.0, 0.0, 1.0));
   inst_col_matrices.push_back(glm::vec4(0.6, 0.3, 0.0, 1.0));
   inst_col_matrices.push_back(glm::vec4(0.5, 0.5, 0.0, 1.0));
   inst_col_matrices.push_back(glm::vec4(0.3, 0.6, 0.0, 1.0));
   inst_col_matrices.push_back(glm::vec4(0.0, 0.9, 0.0, 1.0));
   inst_col_matrices.push_back(glm::vec4(0.0, 0.7, 1.0, 1.0));

   if (false) { // debugging. one red ball at the origin.
      inst_trans_matrices.clear();
      inst_col_matrices.clear();
      inst_col_matrices.push_back(glm::vec4(0.8, 0.0, 0.0, 1.0));
      inst_trans_matrices.push_back(glm::vec3(0,0,0));
   }

   n_instances = inst_trans_matrices.size();
   n_instances_allocated = n_instances;

   // has the vao been genvertexarrayed before this is called?
   glBindVertexArray(vao);

   glGenBuffers(1, &inst_colour_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, inst_colour_buffer_id);
   glBufferData(GL_ARRAY_BUFFER, n_instances * sizeof(glm::vec4), &(inst_col_matrices[0]), GL_STATIC_DRAW);
   glEnableVertexAttribArray(2);
   glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), 0);
   glVertexAttribDivisor(2, 1);

   glGenBuffers(1, &inst_model_translation_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, inst_model_translation_buffer_id);
   glBufferData(GL_ARRAY_BUFFER, n_instances * sizeof (glm::vec3), &(inst_trans_matrices[0]), GL_STATIC_DRAW);
   glEnableVertexAttribArray(3);
   glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), 0);
   glVertexAttribDivisor(3, 1);

   glBindVertexArray(0);
}

void
Mesh::delete_gl_buffers() {

   // 20230828-PE Because Mesh is copied around, this can be called many times if it is
   // called from the destructor. But deleting the GL buffers more than
   // once causes as crash.
   // So let's not delete it at all (Hmm)


   if (false)
      // std::cout << "INFO:: Mesh::delete_gl_buffers() called for mesh \"" << name << "\""
      //           << " with vao " << vao << std::endl;
      logger.log(log_t::INFO, "Mesh::delete_gl_buffers() called for mesh", name, "with vao", vao);

   if (vao == VAO_NOT_SET) {
      std::cout << "ERROR:: Mesh::delete_gl_buffers() called without the VAO set for mesh \""
                << name << "\"" << std::endl;
   } else {
      glBindVertexArray(vao);
      if (buffer_id != 0) { // 0 is not valid
         // useful but noisy
         if (false)
            std::cout << "delete_gl_buffers() deleting buffer_id " << buffer_id << std::endl;
         glDeleteBuffers(1, &buffer_id);
         buffer_id = 0;
      }
      glDeleteBuffers(1, &index_buffer_id);

      if (is_instanced) {
         glDeleteBuffers(1, &inst_model_translation_buffer_id);
         glDeleteBuffers(1, &inst_colour_buffer_id);
         if (is_instanced_with_rts_matrix)
            glDeleteBuffers(1, &inst_rts_buffer_id);
      }
      glDeleteVertexArrays(1, &vao);
      vao = VAO_NOT_SET;
   }
}

void
Mesh::setup_buffers() {

   if (is_headless) return;

   // 20250524-PE this is often not an error with EM maps, say.
#if 0
   if (vertices.empty())
      std::cout << "WARNING:: Mesh::setup_buffers() zero vertices -  probably an error"
                << std::endl;
   if (triangles.empty())
      std::cout << "WARNING:: Mesh::setup_buffers() zero triangles - probably an error"
                << std::endl;
#endif

   if (vertices.empty()) {

#if 0
      // from where was this called then?
#if USE_BACKWARD
               backward::StackTrace st;
               backward::Printer p;
               st.load_here(32);
               p.print(st);
#endif
#endif
   }

   if (false) {
      std::cout << "DEBUG:: in setup_buffers() " << name << " with first_time: " << first_time << std::endl;
      std::cout << "DEBUG:: in setup_buffers() n_vertices: "  << vertices.size()  << std::endl;
      std::cout << "DEBUG:: in setup_buffers() n_triangles: " << triangles.size() << std::endl;
   }

   if (vertices.empty()) return;
   if (triangles.empty() && lines_vertex_indices.empty()) return;

   GLenum err = glGetError();
   if (err) {
      logger.log(log_t::GL_ERROR, logging::function_name_t("Mesh::setup_buffers"),
                 name, stringify_error_code(err));
      err = glGetError();
      if (err != 0)
         std::cout << "GL ERROR:: Mesh::setup_buffers() \"" << name << "\" --- start --- stack-clear "
                   << stringify_error_code(err) << std::endl;
      err = glGetError();
      if (err != 0)
         std::cout << "GL ERROR:: Mesh::setup_buffers() \"" << name << "\" --- start --- stack-clear "
                   << stringify_error_code(err) << std::endl;
   }

   bool setup_buffers_for_gl_lines = false;
   if (! lines_vertex_indices.empty())
      setup_buffers_for_gl_lines = true;

   if (first_time) {
      glGenVertexArrays(1, &vao);
      // std::cout << "DEBUG:: Mesh::setup_buffers() \"" << name << "\" first time: generated VAO " << vao << std::endl;
      // note: don't return before we set first_time = false at the end
   } else {
      // std::cout << "DEBUG:: Mesh::setup_buffers() ######### not first time \"" << name << "\" using VAO " << vao << std::endl;
   }

   // 20220304-PE wondering why binding of this VAO fails? It's because you forgot to setup_buffers()
   // before making new geometry. (Hopefully I will never read this again)
   //
   glBindVertexArray(vao);
   err = glGetError();
   if (err) {
      // 20220803-PE did you forget to attach_buffers() beforehand again?
      std::cout << "GL ERROR:: Mesh::setup_buffers() on binding vao " << vao << " error " << _(err) << std::endl;
      logger.log(log_t::GL_ERROR, logging::function_name_t("Mesh::setup_buffers()"),
		 {"on binding vao", vao, stringify_error_code(err)});
   }

   unsigned int n_vertices = vertices.size();

   // std::cout << "setup_buffers() setup_buffers_for_gl_lines " << setup_buffers_for_gl_lines << std::endl;

   if (first_time) {
      glGenBuffers(1, &buffer_id);
      glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
      glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(s_generic_vertex), &(vertices[0]), GL_STATIC_DRAW);
   } else {
      glDeleteBuffers(1, &buffer_id);
      glGenBuffers(1, &buffer_id);


      if (false)
         std::cout << "DEBUG setup_buffers: buffer_id=" << buffer_id << " about to bind" << std::endl;

      // Check if we have a valid GL context
      GLenum err = glGetError();
      if (err)
         std::cout << "ERROR:: setup_buffers(): GL error before bind: " << err << std::endl;

      GLint current_context = 0;
      glGetIntegerv(GL_CURRENT_PROGRAM, &current_context);
      if (false)
         std::cout << "Current GL context check (program): " << current_context << std::endl;

      glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
      glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(s_generic_vertex), &(vertices[0]), GL_STATIC_DRAW);
   }

   // position
   glEnableVertexAttribArray(0);
   glVertexAttribPointer (0, 3, GL_FLOAT, GL_FALSE, sizeof(s_generic_vertex), 0);

   // normal
   glEnableVertexAttribArray (1);
   glVertexAttribPointer (1, 3, GL_FLOAT, GL_FALSE, sizeof(s_generic_vertex),
                          reinterpret_cast<void *>(sizeof(glm::vec3)));

   // they are s_generic_vertex:
   // layout(location = 0) in vec3 position;
   // layout(location = 1) in vec3 normal;
   // layout(location = 2) in vec4 colour;


   if (false)
      std::cout << "debug:: in Mesh::setup_buffers() is_instanced_colours for mesh with name \"" << name << "\""
                << " is_instanced_colours: " << is_instanced_colours
		<< " (not that that should matter any more)" << std::endl;

   glEnableVertexAttribArray(2);
   glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(s_generic_vertex),
                         reinterpret_cast<void *>(2 * sizeof(glm::vec3)));

   unsigned int n_triangles = triangles.size();
   unsigned int n_bytes_for_triangles = n_triangles * sizeof(g_triangle);
   unsigned int n_bytes_for_lines = lines_vertex_indices.size() * sizeof(unsigned int);

   if (first_time) {
      glGenBuffers(1, &index_buffer_id);
      err = glGetError();
      if (err) {
	 std::cout << "GL ERROR:: Mesh::setup_buffers()\n";
	 logger.log(log_t::GL_ERROR, logging::function_name_t("Mesh::setup_buffers()"),
		    {"on glGenBuffers()", stringify_error_code(err)});
      }
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
      err = glGetError();
      if (err) {
	 std::cout << "GL ERROR:: Mesh::setup_buffers()\n";
	 logger.log(log_t::GL_ERROR, logging::function_name_t("Mesh::setup_buffers()"),
		    {"on glBindBuffer()", stringify_error_code(err)});
      }
   } else {
      glDeleteBuffers(1, &index_buffer_id);
      glGenBuffers(1, &index_buffer_id);
      err = glGetError();
      if (err) {
	 std::cout << "GL ERROR:: Mesh::setup_buffers()\n";
	 logger.log(log_t::GL_ERROR, logging::function_name_t("Mesh::setup_buffers()"),
		    {"on delete and gen", stringify_error_code(err)});
      }
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
      err = glGetError();
      if (err) {
	 std::cout << "GL ERROR:: Mesh::setup_buffers()\n";
	 logger.log(log_t::GL_ERROR, logging::function_name_t("Mesh::setup_buffers()"),
		    {"on glBindBuffer() - not first time:", stringify_error_code(err)});
      }
   }

   if (false)
      std::cout << "debug:: glBufferData for index buffer_id " << index_buffer_id
                << " n_triangles: " << n_triangles
                << " allocating with size: " << n_bytes_for_triangles << " bytes_for_triangles" << std::endl;

   if (setup_buffers_for_gl_lines) {
      // unsigned int n_bytes_for_gl_lines = n_bytes_for_triangles * 2; // 20220424-PE was this until now
      unsigned int n_bytes_for_gl_lines = n_bytes_for_lines;
      // std::cout << "setup_buffers() allocating " << n_bytes_for_gl_lines << " bytes for gl-lines " << std::endl;
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_bytes_for_gl_lines, &lines_vertex_indices[0], GL_STATIC_DRAW);
      err = glGetError();
      if (err)
	 logger.log(log_t::GL_ERROR, logging::function_name_t("Mesh::setup_buffers"),
		    {"on gl_lines glBufferData() ", stringify_error_code(err), "for mesh", name});
   } else {
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_bytes_for_triangles, &triangles[0], GL_STATIC_DRAW);
      err = glGetError();
      if (err) {
	 std::cout << "GL ERROR:: Mesh::setup_buffers()\n";
	 logger.log(log_t::GL_ERROR, logging::function_name_t("Mesh::setup_buffers"),
		    {"on glBufferData() ", stringify_error_code(err), "for mesh", name});
      }
   }

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);

   glBindBuffer(GL_ARRAY_BUFFER, 0);
   glUseProgram(0);
   glBindVertexArray(0);

   first_time = false;

   if (false) {
      std::cout << "********** ending setup_buffers() " << name << " first_time: " << first_time << std::endl;
      std::cout << "********** ending setup_buffers() " << name << " vao: " << vao << std::endl;
   }

}

std::optional<glm::vec3>
Mesh::get_centre_of_mesh() const {

   if (vertices.empty()) {
      return std::nullopt;
   } else {
      size_t n_vertices = vertices.size();
      glm::vec3 sum(0,0,0);
      for (const auto &vertex : vertices) {
         sum += vertex.pos;
      }
      glm::vec3 av(sum.x/static_cast<double>(n_vertices),
                   sum.y/static_cast<double>(n_vertices),
                   sum.z/static_cast<double>(n_vertices));
      return av;
   }
}

std::optional<float>
Mesh::get_radius_of_gyration() const {

   if (vertices.empty()) {
      return std::nullopt;
   } else {
      std::optional<glm::vec3> centre = get_centre_of_mesh();
      if (centre) {
         size_t n_vertices = vertices.size();
         double sum_dd = 0.0;
         for (const auto &vertex : vertices) {
            double dd = glm::distance2(vertex.pos, centre.value());
            sum_dd += dd;
         }
         double rr = sum_dd / static_cast<double>(n_vertices);
         double r = std::sqrt(rr);
         return r;
      } else {
         return std::nullopt;
      }
   }
}



void
Mesh::setup_instanced_debugging_objects(Shader *shader_p, const Material &material_in) {

   material = material_in; // not currently used in the shader.
   shader_p->Use();
   unsigned int n_vert = vertices.size();
   add_one_origin_dodec(); // flat-faced for debugging
   // move the dodec
   for (unsigned int i=n_vert; i<vertices.size(); i++) {
      s_generic_vertex &v = vertices[i];
      v.pos += glm::vec3(-0.05,0.15,0);
   }
   add_one_origin_ball();
   setup_buffers(); // sets the vao
   setup_debugging_instancing_buffers();

}

void
Mesh::setup_instanced_dodecs(Shader *shader_p, const Material &material_in) {

   is_instanced = true;
   is_instanced_colours =  true; //  I think.
   material = material_in; // not currently used in the shader.
   shader_p->Use();
   add_one_origin_dodec(); // flat-faced for debugging
   setup_buffers(); // sets the vao

   setup_debugging_instancing_buffers(); // a few translated objects

}

// rts rotation, translation & scale
void
Mesh::import_and_setup_instanced_cylinders(Shader *shader_p,
                                           const Material &material_in,
                                           const std::vector<glm::mat4> &mats,
                                           const std::vector<glm::vec4> &colours) {

   GLenum err = glGetError();
   if (err) std::cout << "error import_and_setup_instanced_cylinders() -- start -- "
                      << err << std::endl;
   is_instanced = true;
   is_instanced_colours = true;
   is_instanced_with_rts_matrix = true;

   shader_p->Use();
   material = material_in; // not currently used in the shader.
   add_one_origin_cylinder(16);
   setup_buffers();

   n_instances = mats.size();
   n_instances_allocated = n_instances;
   if (false)
      std::cout << "::::::::::::: debug:: import_and_setup_instanced_cylinders() calls "
                << "setup_matrix_and_colour_instancing_buffers_standard" << std::endl;
   setup_matrix_and_colour_instancing_buffers_standard(mats, colours);
   err = glGetError(); if (err) std::cout << "error import_and_setup_instanced_cylinders() -- end -- "
                                          << err << std::endl;

}

void
Mesh::setup_rtsc_instancing(Shader *shader_p,
                            const std::vector<glm::mat4> &mats,
                            const std::vector<glm::vec4> &colours,
                            unsigned int n_instances_in,
                            const Material &material_in) {

   is_instanced = true;
   is_instanced_colours = true;
   is_instanced_with_rts_matrix = true;

   if (shader_p)
      shader_p->Use();
   material = material_in;

   setup_buffers();
   n_instances = n_instances_in;
   n_instances_allocated = n_instances;

   if (false)
      std::cout << "::::::::::::: debug:: setup_rtsc_instancing() calls "
                << "setup_matrix_and_colour_instancing_buffers_standard()"
                << std::endl;
   setup_matrix_and_colour_instancing_buffers_standard(mats, colours);
   GLenum err = glGetError();
   if (err) std::cout << "   error setup_instanced_cylinders() -- end -- "
                      << err << std::endl;

   if (false) {
      std::cout << "setup_rtsc_instancing(): " << vertices.size() << " vertices" << std::endl;
      std::cout << "setup_rtsc_instancing(): " << triangles.size() << " triangles" << std::endl;
   }
}

void
Mesh::setup_instanced_octahemispheres(Shader *shader_p,
                                      const Material &material_in,
                                      const std::vector<glm::mat4> &mats,
                                      const std::vector<glm::vec4> &colours) {

   GLenum err = glGetError(); if (err) std::cout << "   error setup_instanced_octahemispheres() "
                                                 << " -- start -- " << err << std::endl;
   is_instanced = true;
   is_instanced_colours = true;
   is_instanced_with_rts_matrix = true;

   material = material_in; // not currently used in the shader.
   shader_p->Use();
   add_one_origin_octahemisphere(2);
   setup_buffers();
   n_instances = mats.size();
   n_instances_allocated = n_instances;

   std::cout << "::::::::::::: debug:: setup_instanced_octahemispheres() calls"
             << " setup_matrix_and_colour_instancing_buffers_standard()" << std::endl;

   setup_matrix_and_colour_instancing_buffers_standard(mats, colours); // maybe pass a flag MATS_AND_COLOURS
                                            // because we might have
                                            // other instanced geometry that doesn't change colour.
                                            // How about HOLE balls?

   err = glGetError(); if (err) std::cout << "   error setup_instanced_octahemispheres() -- end -- "
                                          << err << std::endl;

}


// instancing buffer for particles. Make *space* for n_particles, but set
// n_instances = 0.
void
Mesh::setup_vertex_and_instancing_buffers_for_particles(unsigned int n_instances_in) {

   bool debug = false;

   // we want to allocate space for n_particles instances, but until the particles
   // are created, and the particle bufffer data updated, we don't want to draw
   // them, so n_instances is 0.
   //
   n_instances = 0;
   n_instances_allocated = n_instances_in;
   particle_draw_count = 0;

   // glm::vec3 n(0,0,1);
   // glm::vec4 c(0.8, 0.4, 0.8, 0.8);

   setup_camera_facing_polygon(5, 0.3, true, 0.3); // calls setup_buffers() for the vertices and sets the VAO

   glBindVertexArray(vao);
   GLenum err = glGetError();
   if (err) {
      std::cout << "GL error ####"
		<< " setup_vertex_and_instancing_buffers_for_particles() B "
		<< err << std::endl;
      logger.log(log_t::GL_ERROR, "setup_vertex_and_instancing_buffers_for_particles() B");
   }


   // a Particle has position, velocity and colour. We need position and colour

   // I shouldn't need to make 2 buffers (ie. 2 calls to glBufferData) here!
   // Look at how the Mesh for ribbons does it.

   // instanced position

   unsigned int n_bytes = n_instances_allocated * sizeof(Particle);

   glGenBuffers(1, &inst_model_translation_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, inst_model_translation_buffer_id);
   glBufferData(GL_ARRAY_BUFFER, n_bytes, nullptr, GL_DYNAMIC_DRAW);
   glEnableVertexAttribArray(3);
   glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(Particle), 0);
   glVertexAttribDivisor(3, 1);
   err = glGetError();
   if (err) {
      std::cout << "GL error #####"
		<< " setup_instancing_buffers_for_particles() B "
		<< err << std::endl;
      logger.log(log_t::GL_ERROR, "setup_instancing_buffers_for_particles() B");
   }

   // instanced colours - setup another buffer - extravagent.
   glGenBuffers(1, &inst_colour_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, inst_colour_buffer_id);
   glBufferData(GL_ARRAY_BUFFER, n_instances_allocated * sizeof(Particle), nullptr, GL_DYNAMIC_DRAW);
   glEnableVertexAttribArray(4);
   glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, sizeof(Particle),
                         reinterpret_cast<void *>(2 * sizeof(glm::vec3)));
   glVertexAttribDivisor(4, 1);

   // index the quad/hex/polygon

   glGenBuffers(1, &index_buffer_id);
   err = glGetError(); if (err) std::cout << "GL error setup_instancing_buffers_for_particles()\n";
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
   err = glGetError(); if (err) std::cout << "GL error setup_instancing_buffers_for_particles()\n";
   unsigned int n_triangles = triangles.size();
   n_bytes = n_triangles * 3 * sizeof(unsigned int);
   if (debug)
      std::cout << "debug:: setup_vertex_and_instancing_buffers_for_particles() "
                << "vao " << vao
                << " glBufferData for index buffer_id " << index_buffer_id
                << " n_triangles: " << n_triangles
                << " allocating with size: " << n_bytes << " bytes for triangles"
                << std::endl;
   glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_bytes, &triangles[0], GL_DYNAMIC_DRAW);
   err = glGetError(); if (err) std::cout << "GL error setup_instancing_buffers_for_particles()\n";

   err = glGetError();
   if (err) std::cout << "GL error #####"
                      << " setup_vertex_and_instancing_buffers_for_particles() --- end --- "
                      << err << std::endl;
   glBindVertexArray(0);

}


void
Mesh::setup_matrix_and_colour_instancing_buffers(const std::vector<glm::mat4> &mats,
                                                 const std::vector<glm::vec4> &colours) {

   // Do you want this function or setup_matrix_and_colour_instancing_buffers_old()?
   // (maybe _old is the wrong suffix, should be _standard)

   std::cout << "--------------------------------------------------------------------------- "
             << "setup_matrix_and_colour_instancing_buffers(): mats size " << mats.size()
             << " colours size " << colours.size() << std::endl;

   GLenum err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() -- start -- "
                      << err << std::endl;

   n_instances = mats.size();
   n_instances_allocated = n_instances;

   if (vao == VAO_NOT_SET)
      std::cout << "GL ERROR:: in setup_matrix_and_colour_instancing_buffers() You didn't correctly setup this Mesh "
                << name << " " << std::endl;

   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "ERROR:: setup_matrix_and_colour_instancing_buffers() B binding-vao, with vao "
                      << vao << " err: " << err << std::endl;


   // ------------- orientation/position/scale matrices ----------------------------

   glGenBuffers(1, &inst_rts_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, inst_rts_buffer_id);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() C0 "
                      << err << std::endl;
   glBufferData(GL_ARRAY_BUFFER, n_instances * sizeof(glm::mat4), &(mats[0]), GL_DYNAMIC_DRAW);
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() C1 " << err << std::endl;

   glEnableVertexAttribArray(3);
   glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), 0);
   glVertexAttribDivisor(3, 1);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() C3 " << err << std::endl;
   glEnableVertexAttribArray(4);
   glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), reinterpret_cast<void *>(sizeof(glm::vec4)));
   glVertexAttribDivisor(4, 1);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() C4 " << err << std::endl;
   glEnableVertexAttribArray(5);
   glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), reinterpret_cast<void *>(2 * sizeof(glm::vec4)));
   glVertexAttribDivisor(5, 1);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() C5 " << err << std::endl;
   glEnableVertexAttribArray(6);
   glVertexAttribPointer(6, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), reinterpret_cast<void *>(3 * sizeof(glm::vec4)));
   glVertexAttribDivisor(6, 1);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() C6 " << err << std::endl;

   // ------------- colours ----------------------------

   glGenBuffers(1, &inst_colour_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, inst_colour_buffer_id);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() B0 "
                      << err << std::endl;
   glBufferData(GL_ARRAY_BUFFER, n_instances * sizeof(glm::vec4), &(colours[0]), GL_DYNAMIC_DRAW);
   glEnableVertexAttribArray(7);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() B1 "
                      << err << std::endl;
   glVertexAttribPointer(7, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), 0);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() B2 "
                      << err << std::endl;
   glVertexAttribDivisor(7, 1);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers() B3 "
                      << err << std::endl;

}


void
Mesh::setup_matrix_and_colour_instancing_buffers_standard(const std::vector<glm::mat4> &mats,
                                                          const std::vector<glm::vec4> &colours) {

   // this function doesn't make sense in it's current form.
   // Instead, I need to *make space* for n_mats and n_colours.
   //
   // I don't need matrices with values here - because they will be updated/replaced by
   // the tick function. So this function should be changed to pass the size of
   // matrix and colour vectors, not the actual values.
   //
   // send only the size (they should both be the same size) (they can be updated differently)

   if (false)
      std::cout << "Mesh::setup_matrix_and_colour_instancing_buffers_standard(): mats size " << mats.size()
                << " colours size " << colours.size() << " and is_instanced_colours " << is_instanced_colours
                << std::endl;

   GLenum err = glGetError();
   if (err) std::cout << "Error setup_matrix_and_colour_instancing_buffers_standard() -- start -- "
                      << err << std::endl;

   n_instances = mats.size();
   n_instances_allocated = n_instances;

   if (false)
      std::cout << "in setup_matrix_and_colour_instancing_buffers_standard() "
                << "n_instances " << n_instances << std::endl;

   const std::vector<glm::mat4> &inst_rts_matrices = mats;
   const std::vector<glm::vec4> &inst_col_matrices = colours;

   // these vectors should be the same size. Add a check.

   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers_standard() A "
                      << err << std::endl;

   if (vao == VAO_NOT_SET)
      std::cout << "GL ERROR:: in setup_matrix_and_colour_instancing_buffers_standard() You didn't correctly setup this Mesh "
                << name << " " << std::endl;

   glBindVertexArray(vao);

   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::setup_matrix_and_colour_instancing_buffers_standard() B binding-vao "
                      << err << " with vao " << vao << std::endl;

   // -------- rotation/translation/scale matrices -----------

   // 20210827-PE Here I should clean out the old buffer data before overriding it.
   // glDeleteBuffers() if this is not the first time that this function has been called.

   glGenBuffers(1, &inst_rts_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, inst_rts_buffer_id);
   if (false)
      std::cout << "setup_matrix_and_colour_instancing_buffers_standard() allocating matrix buffer data "
                << n_instances * 4 * sizeof(glm::mat4) << std::endl;
   glBufferData(GL_ARRAY_BUFFER, n_instances * 4 * sizeof (glm::vec4), &(inst_rts_matrices[0]), GL_DYNAMIC_DRAW); // dynamic

   err = glGetError(); if (err) std::cout << "GL ERROR:: setup_instancing_buffers() C1 " << err << std::endl;
   glEnableVertexAttribArray(3);
   glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), 0);
   glVertexAttribDivisor(3, 1);
   err = glGetError(); if (err) std::cout << "   error setup_instancing_buffers() C2 " << err << std::endl;
   glEnableVertexAttribArray(4);
   glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), reinterpret_cast<void *>(sizeof(glm::vec4)));
   glVertexAttribDivisor(4, 1);
   err = glGetError(); if (err) std::cout << "   error setup_instancing_buffers() C2 " << err << std::endl;
   glEnableVertexAttribArray(5);
   glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), reinterpret_cast<void *>(2 * sizeof(glm::vec4)));
   glVertexAttribDivisor(5, 1);
   err = glGetError(); if (err) std::cout << "   error setup_instancing_buffers() C3 " << err << std::endl;
   glEnableVertexAttribArray(6);
   glVertexAttribPointer(6, 4, GL_FLOAT, GL_FALSE, sizeof(glm::mat4), reinterpret_cast<void *>(3 * sizeof(glm::vec4)));
   glVertexAttribDivisor(6, 1);

   err = glGetError(); if (err) std::cout << "   error setup_instancing_buffers() E " << err << std::endl;

   // rama balls we want to have instanced colours but hydrogen bond rotating cylinders we do not.

   // -------- colours -----------

   glGenBuffers(1, &inst_colour_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, inst_colour_buffer_id);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: setup_matrix_and_colour_instancing_buffers_standard() B0 "
                      << err << std::endl;
   if (false)
      std::cout << "setup_matrix_and_colour_instancing_buffers_old() allocating colour buffer data "
                << n_instances * sizeof(glm::vec4) << std::endl;
   glBufferData(GL_ARRAY_BUFFER, n_instances * sizeof(glm::vec4), &(inst_col_matrices[0]), GL_DYNAMIC_DRAW); // dynamic
   glEnableVertexAttribArray(7);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers_standard() B1 "
                      << err << std::endl;
   glVertexAttribPointer(7, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), 0);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers_standard() B2 "
                      << err << std::endl;
   glVertexAttribDivisor(7, 1);
   err = glGetError();
   if (err) std::cout << "error setup_matrix_and_colour_instancing_buffers_standard() B3 "
                      << err << std::endl;

   
}


// void
// Mesh::setup_instanced_balls(Shader *shader_p, const Material &material_in) {

// }




void
Mesh::fill_with_simple_triangles_vertices() {

   float scale = 0.4;

   std::vector<s_generic_vertex> v(6);
   float z_sep = 0.4;

   v[0].pos = scale * glm::vec3(  0.0f,  0.5f,  -z_sep);
   v[1].pos = scale * glm::vec3(  0.5f, -0.36f, -z_sep);
   v[2].pos = scale * glm::vec3( -0.5f, -0.36f, -z_sep);
   v[3].pos = scale * glm::vec3(  0.0f,  0.5f,   z_sep);
   v[4].pos = scale * glm::vec3(  0.5f, -0.36f,  z_sep);
   v[5].pos = scale * glm::vec3( -0.5f, -0.36f,  z_sep);

   v[0].normal = glm::vec3( 0.2f, 0.2f,  0.9f);
   v[1].normal = glm::vec3( 0.2f, 0.9f,  0.2f);
   v[2].normal = glm::vec3( 0.9f, 0.3f,  0.1f);
   v[3].normal = glm::vec3( 0.0f, 0.9f, -0.1f);
   v[4].normal = glm::vec3( 0.9f, 0.3f, -0.2f);
   v[5].normal = glm::vec3( 0.2f, 0.6f, -0.9f);

   v[0].color = glm::vec4(0.0f, 0.0f, 0.0f, 1.f);
   v[1].color = glm::vec4(0.2f, 0.3f, 1.0f, 1.f);
   v[2].color = glm::vec4(0.5f, 0.9f, 0.2f, 1.f);
   v[3].color = glm::vec4(0.2f, 0.2f, 0.9f, 1.f);
   v[4].color = glm::vec4(0.1f, 0.9f, 0.2f, 1.f);
   v[5].color = glm::vec4(0.9f, 0.3f, 0.2f, 1.f);

   unsigned int idx_base = vertices.size();
   vertices.insert(vertices.end(), v.begin(), v.end());

   g_triangle gt_0(idx_base,   idx_base+1, idx_base+2);
   g_triangle gt_1(idx_base+3, idx_base+4, idx_base+5);

   triangles.push_back(gt_0);
   triangles.push_back(gt_1);

}

void
Mesh::fill_with_direction_triangles() {

   std::vector<s_generic_vertex> v;
   unsigned int v_size = 3 * 3; // bypass weird flycheck bug
   v.resize(v_size);
   float scale = 0.25;

   v[0].pos = scale * glm::vec3( -0.2f, 0.0f, 0.0f);
   v[1].pos = scale * glm::vec3(  0.2f, 0.0f, 0.0f);
   v[2].pos = scale * glm::vec3(  0.0f, 0.0f, 0.5f);

   v[3].pos = scale * glm::vec3( -0.2f,  0.0f, 0.0f);
   v[4].pos = scale * glm::vec3(  0.2f,  0.0f, 0.0f);
   v[5].pos = scale * glm::vec3(  0.0f,  0.5f, 0.0f);

   v[6].pos = scale * glm::vec3(  0.0f, -0.2f, 0.01f);
   v[7].pos = scale * glm::vec3(  0.0f,  0.2f, 0.01f);
   v[8].pos = scale * glm::vec3(  0.5f,  0.0f, 0.0f);

   v[0].normal = glm::vec3( 0.2f, 0.2f,  0.9f);
   v[1].normal = glm::vec3( 0.2f, 0.9f,  0.2f);
   v[2].normal = glm::vec3( 0.9f, 0.1f,  0.1f);
   v[3].normal = glm::vec3( 0.0f, 0.9f, -0.1f);
   v[4].normal = glm::vec3( 0.9f, 0.3f, -0.2f);
   v[5].normal = glm::vec3( 0.1f, 0.9f, -0.1f);
   v[6].normal = glm::vec3( 0.0f, 0.9f, -0.1f);
   v[7].normal = glm::vec3( 0.9f, 0.3f, -0.2f);
   v[8].normal = glm::vec3( 0.1f, 0.1f, -0.9f);

   v[0].color = glm::vec4(0.8f, 0.0f, 0.0f, 1.f);
   v[1].color = glm::vec4(0.8f, 0.3f, 1.0f, 1.f);
   v[2].color = glm::vec4(0.8f, 0.1f, 0.1f, 1.f);
   v[3].color = glm::vec4(0.2f, 0.8f, 0.9f, 1.f);
   v[4].color = glm::vec4(0.1f, 0.9f, 0.2f, 1.f);
   v[5].color = glm::vec4(0.1f, 0.8f, 0.1f, 1.f);
   v[6].color = glm::vec4(0.4f, 0.2f, 0.3f, 1.f);
   v[7].color = glm::vec4(0.1f, 0.4f, 0.3f, 1.f);
   v[8].color = glm::vec4(0.1f, 0.1f, 0.9f, 1.f);

   unsigned int base = vertices.size();
   vertices.insert(vertices.end(), v.begin(), v.end());
 
   g_triangle gt_0(base,base+1,base+2);
   g_triangle gt_1(base+3,base+4,base+5);
   g_triangle gt_2(base+6,base+7,base+8);
   triangles.push_back(gt_0);
   triangles.push_back(gt_1);
   triangles.push_back(gt_2);

}


void
Mesh::draw_normals(const glm::mat4 &mvp, float normal_scaling) {

   GLenum err = glGetError(); if (err) std::cout << "   error draw_normals() -- start -- "
                                                 << err << std::endl;

   if (! normals_are_setup) {
      glGenVertexArrays(1, &normals_vao);
      std::cout << "####### draw_normals() new normals_vao " << normals_vao << std::endl;
   }
   glBindVertexArray(normals_vao);

   auto vec_length = [] (const glm::vec3 &v) {
                        float s = v.x * v.x + v.y * v.y + v.z * v.z;
                        return sqrtf(s);
                     };

   if (shader_for_draw_normals.name.empty())
      shader_for_draw_normals.init("draw-normals.shader", Shader::Entity_t::GENERIC_DISPLAY_OBJECT);

   shader_for_draw_normals.Use();
   glUniformMatrix4fv(shader_for_draw_normals.mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);

   err = glGetError(); if (err) std::cout << "   error draw_normals() post mvp uniform "
                                          << err << std::endl;

   std::vector<glm::vec3> tmp_normals; // these are positions, 2 vertices per normal
   std::vector<glm::vec4> tmp_colours; // 2 colours per normal, one for each of the tmp_normals

   if (! normals_are_setup) {
      for (unsigned int i=0; i<vertices.size(); i++) {
         float f = vec_length(vertices[i].normal);
         if (false)
            std::cout << i
                      << " vertex: " << glm::to_string(vertices[i].pos)
                      << " normal: " << glm::to_string(vertices[i].normal)
                      << " length " << f << std::endl;
         tmp_normals.push_back(vertices[i].pos);
         tmp_normals.push_back(vertices[i].pos + normal_scaling * vertices[i].normal);
         tmp_colours.push_back(glm::vec4(0.4, 0.4, 0.7, 1.0));
         tmp_colours.push_back(glm::vec4(0.4, 0.4, 0.7, 1.0));
      }

      for (unsigned int i=0; i<triangles.size(); i++) {
         const glm::vec3 &p0 = vertices[triangles[i].point_id[0]].pos;
         const glm::vec3 &p1 = vertices[triangles[i].point_id[1]].pos;
         const glm::vec3 &p2 = vertices[triangles[i].point_id[2]].pos;
         glm::vec3 delta_0(0.0001, 0.0001, 0.0001);
         glm::vec3 delta_1(0.0001, 0.0001, 0.0001);
         glm::vec3 delta_2(0.0001, 0.0001, 0.0001);
         if (p0.x < 0) delta_0.x *= -1.0;
         if (p0.y < 0) delta_0.y *= -1.0;
         if (p0.z < 0) delta_0.z *= -1.0;
         if (p1.x < 0) delta_1.x *= -1.0;
         if (p1.y < 0) delta_1.y *= -1.0;
         if (p1.z < 0) delta_1.z *= -1.0;
         if (p2.x < 0) delta_2.x *= -1.0;
         if (p2.y < 0) delta_2.y *= -1.0;
         if (p2.z < 0) delta_2.z *= -1.0;
         tmp_normals.push_back(p0+delta_0);
         tmp_normals.push_back(p1+delta_1);
         tmp_normals.push_back(p0+delta_0);
         tmp_normals.push_back(p2+delta_2);
         tmp_normals.push_back(p1+delta_1);
         tmp_normals.push_back(p2+delta_2);
         tmp_colours.push_back(glm::vec4(0.5, 0.5, 0.5, 1.0));
         tmp_colours.push_back(glm::vec4(0.5, 0.5, 0.5, 1.0));
         tmp_colours.push_back(glm::vec4(0.5, 0.5, 0.5, 1.0));
         tmp_colours.push_back(glm::vec4(0.5, 0.5, 0.5, 1.0));
         tmp_colours.push_back(glm::vec4(0.5, 0.5, 0.5, 1.0));
         tmp_colours.push_back(glm::vec4(0.5, 0.5, 0.5, 1.0));
      }
   }

   unsigned int n_normals = triangles.size() * 6 + 2 * vertices.size();

   if (! normals_are_setup)
      glGenBuffers(1, &normals_buffer_id);

   glBindBuffer(GL_ARRAY_BUFFER, normals_buffer_id);
   if (! normals_are_setup)
      glBufferData(GL_ARRAY_BUFFER, n_normals * sizeof(glm::vec3), &(tmp_normals[0]), GL_STATIC_DRAW);
   glEnableVertexAttribArray(0);
   glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), 0);

   if (! normals_are_setup)
      glGenBuffers(1, &normals_colour_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, normals_colour_buffer_id);
   if (! normals_are_setup)
      glBufferData(GL_ARRAY_BUFFER, n_normals * sizeof(glm::vec4), &(tmp_colours[0]), GL_STATIC_DRAW);
   glEnableVertexAttribArray(1);
   glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(glm::vec4), 0);
   unsigned int n_vertices = 2 * n_normals;
   glDrawArrays(GL_LINES, 0, n_vertices); // number of vertices, not number of lines
   err = glGetError(); if (err) std::cout << "   error draw_normals() post gldrawarrays "
                                          << err << std::endl;

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);

   if (! normals_are_setup)
      normals_are_setup = true;

}


// This function is used (only) by molecules_as_meshes (currently disabled - 20210824). It has a different
// layout to the draw() (which can/does draw other instances - e.g. hydrogen bond cylinders).
// i.e. instanced models are not the same as instanced objects - don't mix them up.
// (I am not sure why this function is neeeded now - perhaps I didn't want to mess up the working
/// draw() function as I was developing it.)
//
// molecules_as_meshes seems faster than current molecules, but I can't get the (specular) lighting right.
//
// There are several things we should not pass to the shader when pass_type is PASS_TYPE_FOR_SHADOWS
// so we need a flag to know to turn them off.
void
Mesh::draw_instanced(int pass_type,
                     Shader *shader_p,
                     stereo_eye_t eye,
                     const glm::mat4 &mvp,
                     const glm::mat4 &view_rotation_matrix,
                     const std::map<unsigned int, lights_info_t> &lights,
                     const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                     const glm::vec4 &background_colour,
                     bool do_depth_fog,
                     bool transferred_colour_is_instanced,
                     bool do_pulse,                    // optional args
                     bool do_rotate_z,
                     float pulsing_amplitude,
                     float pulsing_frequency,
                     float pulsing_phase_distribution,
                     float z_rotation_angle) {

   // debug_mode = true;
   if (debug_mode)
      std::cout << "DEBUG:: Mesh::draw_instanced() Mesh " << name << " -- start -- with shader " << shader_p->name
                << "pass_type: " << pass_type << " and do_pulse " << do_pulse
		<< " and draw_this_mesh " << draw_this_mesh << std::endl;

   if (! draw_this_mesh) return;

   unsigned int n_triangles = triangles.size();
   GLuint n_verts = 3 * n_triangles;

   if (n_triangles == 0) return;

   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw_instanced() " << name << " " << shader_p->name
                      << " -- start -- " << err << std::endl;
   shader_p->Use();
   const std::string &shader_name = shader_p->name;

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: " << shader_p->name << " draw_instanced() post mvp uniform "
                      << err << std::endl;

   // the view_rotation_matrix is not used to calculate glPosition in instanced-objects.shader
   //
   std::pair<float, float> stereo_x_scale_and_offset = get_stereo_x_scale_and_offset(eye);
   const float &stereo_x_scale  = stereo_x_scale_and_offset.first;
   const float &stereo_x_offset = stereo_x_scale_and_offset.second;

   glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE, &view_rotation_matrix[0][0]);
   err = glGetError();
   if (err)
      std::cout << "GL ERROR:: Mesh::draw_instanced() " << name << " " << shader_p->name
                << " draw_instanced() post view rotation uniform " << err << std::endl;

   std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
   float time = std::chrono::duration_cast<std::chrono::milliseconds>(now - time_constructed).count();

   // std::cout << "sending time " << time << std::endl;
   shader_p->set_float_for_uniform("time", time);

   std::map<unsigned int, lights_info_t>::const_iterator it;
   unsigned int light_idx = 0;
   it = lights.find(light_idx);
   if (it != lights.end()) {
      shader_p->setup_light(light_idx, it->second, view_rotation_matrix);
   }
   light_idx = 1;
   it = lights.find(light_idx);
   if (it != lights.end()) {
      shader_p->setup_light(light_idx, it->second, view_rotation_matrix);
   }

   shader_p->set_float_for_uniform("stereo_x_scale",  stereo_x_scale);
   shader_p->set_float_for_uniform("stereo_x_offset", stereo_x_offset);
   shader_p->set_vec4_for_uniform("background_colour", background_colour);
   shader_p->set_bool_for_uniform("do_depth_fog", do_depth_fog);

   shader_p->set_bool_for_uniform("transferred_colour_is_instanced", transferred_colour_is_instanced);

   err = glGetError();
   if (err) std::cout << "GL ERROR:: draw_instanced() pre-setting material " << err << std::endl;
   shader_p->set_vec4_for_uniform( "material.ambient",   material.ambient);
   shader_p->set_vec4_for_uniform( "material.diffuse",   material.diffuse);
   shader_p->set_vec4_for_uniform( "material.specular",  material.specular);
   shader_p->set_float_for_uniform("material.shininess", material.shininess);
   shader_p->set_float_for_uniform("material.specular_strength", material.specular_strength);
   err = glGetError();
   if (err) std::cout << "GL ERROR draw_instanced(): " << shader_name << " post-material "
                      << " with GL err " << err << std::endl;
   shader_p->set_vec3_for_uniform("eye_position", eye_position);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw_instanced() \"" << name << "\" \"" << shader_name << "\" post-set eye position "
                      << " with GL err " << err << std::endl;
   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw_instanced() " << shader_name << " pre-glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;

   if (vao == VAO_NOT_SET)
      std::cout << "GL ERROR:: You forgot to setup this Mesh: \"" << name << "\" " << shader_p->name << std::endl;

   // std::cout << "Mesh::draw_instanced() using vao " << vao << std::endl;
   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw_instanced() " << shader_name << " " << name
                      << " glBindVertexArray() vao " << vao << " with GL err " << err << std::endl;

   shader_p->set_bool_for_uniform("do_pulse", do_pulse);
   shader_p->set_bool_for_uniform("do_rotate_z", do_rotate_z);
   shader_p->set_float_for_uniform("pulsing_amplitude", pulsing_amplitude);
   shader_p->set_float_for_uniform("pulsing_frequency", pulsing_frequency);
   shader_p->set_float_for_uniform("pulsing_phase_distribution", pulsing_phase_distribution);
   shader_p->set_float_for_uniform("z_rotation_angle", z_rotation_angle);

   // fix these                                       
   // glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   // err = glGetError(); if (err) std::cout << "error draw_instanced() glBindBuffer() v " << err << std::endl;
   // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
   // err = glGetError(); if (err) std::cout << "error draw_instanced() glBindBuffer() i " << err << std::endl;

   glEnableVertexAttribArray(0);  // vec3 position
   glEnableVertexAttribArray(1);  // vec3 normal
   glEnableVertexAttribArray(2);  // vec4 colour
   glEnableVertexAttribArray(3);  // instanced mat-row-0 from setup_matrix_and_colour_instancing_buffers_standard()
   glEnableVertexAttribArray(4);  // instanced mat-row-1
   glEnableVertexAttribArray(5);  // instanced mat-row-2
   glEnableVertexAttribArray(6);  // instanced mat-row-3
   glEnableVertexAttribArray(7);  // instanced colour

   // glBindBuffer(GL_ARRAY_BUFFER, inst_rts_buffer_id); // needed?
   // err = glGetError(); if (err) std::cout << "error draw_instanced() glBindBuffer() inst_rts_buffer_id" << std::endl;

   if (debug_mode)
      std::cout << "Mesh::draw_instanced() Mesh \"" << name << "\" drawing n_verts " << n_verts << " n_instances " << n_instances
                << " with shader " << shader_p->name << " and vao " << vao << std::endl;

   err = glGetError();
   if (err) std::cout << "GL_ERROR pre glDrawElementsInstanced()" << std::endl;

   if (false)
      std::cout << "draw_instanced() pre-glDrawElementsInstanced()"
                << " shader: " << shader_p->name << " vao: " << vao
                << " n_triangle_verts: " << n_verts << " n_instances: " << n_instances
                << " with GL err " << err << std::endl;

   glDrawElementsInstanced(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr, n_instances);
   err = glGetError();
   if (err) std::cout << "error draw_instanced() glDrawElementsInstanced()"
                      << " shader: " << shader_p->name << " vao: " << vao
                      << " n_triangle_verts: " << n_verts << " n_instances: " << n_instances
                      << " with GL err " << err << std::endl;
   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);
   glDisableVertexAttribArray(4);
   glDisableVertexAttribArray(5);
   glDisableVertexAttribArray(6);
   glDisableVertexAttribArray(7);
   glUseProgram(0);

}

// // keep these because they are s_generic_vertex, made in setup_buffers()
// layout(location = 0) in vec3 position;
// layout(location = 1) in vec3 normal;
// layout(location = 2) in vec4 colour;

// // c.f. extra-distance-restraints-markup.hh
// layout(location = 3) in float width;  // these are all instanced.
// layout(location = 4) in float length;
// layout(location = 5) in vec3 position;
// layout(location = 6) in mat3 orientation;
// layout(location = 7) in vec4 colour_instanced;
void
Mesh::draw_extra_distance_restraint_instances(Shader *shader_p,
                                              const glm::mat4 &mvp,
                                              const glm::mat4 &view_rotation_matrix,
                                              const std::map<unsigned int, lights_info_t> &lights,
                                              const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                                              const glm::vec4 &background_colour,
                                              bool do_depth_fog) {

   if (! draw_this_mesh) return;
   unsigned int n_triangles = triangles.size();
   GLuint n_verts = 3 * n_triangles;

   if (n_triangles == 0) return;

   GLenum err = glGetError();
   if (err) std::cout << "error Mesh::draw_instanced() " << name << " " << shader_p->name
                      << " -- start -- " << err << std::endl;
   shader_p->Use();
   const std::string &shader_name = shader_p->name;

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   err = glGetError();
   if (err) std::cout << "error:: " << shader_p->name << " draw_extra_distance_restraint_instances() post mvp uniform "
                      << err << std::endl;

   glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE, &view_rotation_matrix[0][0]);
   err = glGetError();
   if (err) std::cout << "error:: Mesh::draw_instanced() " << name << " " << shader_p->name
                      << " draw_instanced() post view rotation uniform " << err << std::endl;

   std::map<unsigned int, lights_info_t>::const_iterator it;
   unsigned int light_idx = 0;
   it = lights.find(light_idx);
   if (it != lights.end())
      shader_p->setup_light(light_idx, it->second, view_rotation_matrix);
   light_idx = 1;
   it = lights.find(light_idx);
   if (it != lights.end())
      shader_p->setup_light(light_idx, it->second, view_rotation_matrix);

   shader_p->set_vec4_for_uniform("background_colour", background_colour);
   shader_p->set_bool_for_uniform("do_depth_fog", do_depth_fog);

   err = glGetError();
   if (err) std::cout << "GL ERROR:: draw_instanced() pre-setting material " << err << std::endl;
   shader_p->set_vec4_for_uniform( "material.ambient",   material.ambient);
   shader_p->set_vec4_for_uniform( "material.diffuse",   material.diffuse);
   shader_p->set_vec4_for_uniform( "material.specular",  material.specular);
   shader_p->set_float_for_uniform("material.shininess", material.shininess);
   shader_p->set_float_for_uniform("material.specular_strength", material.specular_strength);
   err = glGetError();
   if (err) std::cout << "GL ERROR draw_instanced(): " << shader_name << " post-material "
                      << " with GL err " << err << std::endl;
   shader_p->set_vec3_for_uniform("eye_position", eye_position);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw_instanced() \"" << name << "\" \"" << shader_name << "\" post-set eye position "
                      << " with GL err " << err << std::endl;
   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw_instanced() " << shader_name << " pre-glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;

   if (vao == VAO_NOT_SET)
      std::cout << "GL ERROR:: You forgot to setup this Mesh " << name << " " << shader_p->name << std::endl;

   // std::cout << "Mesh::draw_instanced() using vao " << vao << std::endl;
   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw_instanced() " << shader_name << " " << name
                      << " glBindVertexArray() vao " << vao << " with GL err " << err << std::endl;

   glEnableVertexAttribArray(0);  // vec3 position
   glEnableVertexAttribArray(1);  // vec3 normal
   glEnableVertexAttribArray(2);  // vec4 colour // not used in the shader - but assigned in setup_buffers().
   glEnableVertexAttribArray(3);  // width
   glEnableVertexAttribArray(4);  // length
   glEnableVertexAttribArray(5);  // position
   glEnableVertexAttribArray(6);  // mat3 orientation 0
   glEnableVertexAttribArray(7);  // mat3 orientation 1
   glEnableVertexAttribArray(8);  // mat3 orientation 2
   glEnableVertexAttribArray(9);  // colour_instanced

   if (false)
      std::cout << "draw_extra_distance_restraint_instances() n_verts: " << n_verts
                << " n_instances " << n_instances << " with shader " << shader_p->name << std::endl;

   // hack/test
   // n_instances = 200;

   glDrawElementsInstanced(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr, n_instances);
   err = glGetError();
   if (err) std::cout << "error draw_instanced() glDrawElementsInstanced()"
                      << " shader: " << shader_p->name << " vao: " << vao
                      << " n_triangle_verts: " << n_verts << " n_instances: " << n_instances
                      << " with GL err " << err << std::endl;

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);
   glDisableVertexAttribArray(4);
   glDisableVertexAttribArray(5);
   glDisableVertexAttribArray(6);
   glDisableVertexAttribArray(7);
   glDisableVertexAttribArray(8);
   glDisableVertexAttribArray(9);

   glUseProgram(0);

}



void
Mesh::draw_particles(Shader *shader_p, const glm::mat4 &mvp, const glm::mat4 &view_rotation) {

   debug_mode = false; // this is not the place for this.

   if (debug_mode)
      std::cout << "in draw_particles() with n_instances " << n_instances << " and n_triangles: "
                << triangles.size() << std::endl;

   // this can happen when all the particles have life 0 - and have been removed.
   if (n_instances == 0) return;
   if (triangles.empty()) return;

   particle_draw_count += 1;
   shader_p->Use();
   glBindVertexArray(vao);
   GLenum err = glGetError();
   if (err) std::cout << "error draw_particles() " << shader_p->name
                      << " glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;

   glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   err = glGetError(); if (err) std::cout << "GL ERROR:: Mesh::draw_particles() glBindBuffer() v "
                                          << err << std::endl;
   glEnableVertexAttribArray(0); // vertex positions
   glEnableVertexAttribArray(1); // vertex normal
   glEnableVertexAttribArray(2); // vertex colours

   glBindBuffer(GL_ARRAY_BUFFER, inst_colour_buffer_id);
   glEnableVertexAttribArray(3); // instanced model/particle colours (time varying)

   glBindBuffer(GL_ARRAY_BUFFER, inst_model_translation_buffer_id);
   glEnableVertexAttribArray(4); // instanced model/particle translations (time varying)

   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw_particles() " << shader_p->name
                      << " draw() ___PRE___ mvp uniform " << err << std::endl;

   if (debug_mode)
      std::cout << "debug:: Mesh::draw_particles()    sending mvp " << glm::to_string(mvp) << std::endl;
   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, glm::value_ptr(mvp));
   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw_particles() " << shader_p->name
                      << " draw_particles() post mvp uniform " << err << std::endl;
   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw_particles() " << shader_p->name
                      << " draw_particles() post mvp uniform 2 " << err << std::endl;

   if (debug_mode)
      std::cout << "debug sending view_rotation " << glm::to_string(view_rotation) << std::endl;

   glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE, glm::value_ptr(view_rotation));
   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw_particles() " << shader_p->name
                      << " draw_particles() post view_rotation uniform " << err << std::endl;
   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw_particles() " << shader_p->name
                      << " draw_particles() post view_rotation uniform 2 " << err << std::endl;
   //
   float rotation_angle = 0.05f * static_cast<float>(particle_draw_count);
   // std::cout << "Mesh::draw_particles() sending rotation_angle " << rotation_angle << std::endl;
   shader_p->set_float_for_uniform("rotation_angle", rotation_angle);

   glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   unsigned int n_verts = 3 * triangles.size();
   if (debug_mode)
      std::cout << "debug:: Mesh::draw_particles() " << name << " with shader \"" << shader_p->name << "\""
                << " drawing n_instances " << n_instances << std::endl;
   glDrawElementsInstanced(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr, n_instances);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw_particles() " << shader_p->name
                      << " glDrawElementsInstanced() vao " << vao
                      << " with GL err " << err << std::endl;

   glDisable(GL_BLEND);
}

// set the glLineWidth before this draw call (if drawing with lines)
void
Mesh::draw(Shader *shader_p,
           stereo_eye_t eye,
           const glm::mat4 &mvp,
           const glm::mat4 &mouse_based_rotation_matrix,
           const std::map<unsigned int, lights_info_t> &lights,
           const glm::vec3 &eye_position,
           const glm::vec3 &rotation_centre,
           float opacity, // map_opacity
           const glm::vec4 &background_colour,
           bool draw_as_lines_flag, // or surface mesh
           bool do_depth_fog,
           bool show_just_shadows) {

   if (false)
      std::cout << "debug:: Mesh::draw() \"" << name << "\" shader: " << shader_p->name
                << " draw_this_mesh: " << draw_this_mesh
                << " n-vertices:" << vertices.size()
                << " n-tris:" << triangles.size()
                << " is_instanced " <<  is_instanced
                << " is_instanced_colours " <<  is_instanced_colours
                << " is_instanced_with_rts_matrix " <<  is_instanced_with_rts_matrix
                << " draw_as_lines_flag " << draw_as_lines_flag
                << std::endl;

   if (! draw_this_mesh) return;

   const std::string &shader_name = shader_p->name;

   unsigned int n_triangles = triangles.size();
   GLuint n_verts = 3 * n_triangles;

   // At the moment, I don't think that it's an error to come here if there are no triangles (yet)
   // Just quietly do nothing and return.
   if (n_triangles == 0 && lines_vertex_indices.empty()) return;

   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw() " << name << " shader: " << shader_p->name
                      << " -- start -- " << err << std::endl;
   shader_p->Use();
   err = glGetError();
   if (err)
      std::cout << "GL ERROR:: Mesh::draw() post Use() " << name
                << " shader: " << shader_p->name
                << " error: " << stringify_error_code(err) << std::endl;

   GLint mvp_uniform_location = glGetUniformLocation(shader_p->program_id, "mvp");
   glUniformMatrix4fv(mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   err = glGetError();
   if (err)
      std::cout << "GL ERROR:: Mesh::draw() " << name
                << " shader: " << shader_p->name << " post mvp uniform"
                << " location " << shader_p->mvp_uniform_location
                << " error: " << stringify_error_code(err) << std::endl;

   glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE,
                      &mouse_based_rotation_matrix[0][0]);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw() " << name << " " << shader_p->name
                      << " draw() post view rotation uniform " << err << std::endl;

   std::pair<float, float> stereo_x_scale_and_offset = get_stereo_x_scale_and_offset(eye);
   const float &stereo_x_scale  = stereo_x_scale_and_offset.first;
   const float &stereo_x_offset = stereo_x_scale_and_offset.second;

   std::map<unsigned int, lights_info_t>::const_iterator it;
   unsigned int light_idx = 0;
   it = lights.find(light_idx);
   if (it != lights.end()) {
      shader_p->setup_light(light_idx, it->second, mouse_based_rotation_matrix);
   }
   light_idx = 1;
   it = lights.find(light_idx);
   if (it != lights.end()) {
      shader_p->setup_light(light_idx, it->second, mouse_based_rotation_matrix);
   }
   shader_p->setup_eye_position(eye_position, rotation_centre, mouse_based_rotation_matrix);

   // add material properties, use class built-ins this time.

   // not using the built-in - hmm.
   shader_p->set_vec4_for_uniform("background_colour", background_colour);

   shader_p->set_bool_for_uniform("do_depth_fog", do_depth_fog);

   err = glGetError(); if (err) std::cout << "   error draw() pre-setting material "
                                          << err << std::endl;
   shader_p->set_bool_for_uniform( "do_specular",        material.do_specularity); // change from object in C++ to global unform in GLSL
   shader_p->set_vec4_for_uniform( "material.ambient",   material.ambient);
   shader_p->set_vec4_for_uniform( "material.diffuse",   material.diffuse);
   shader_p->set_vec4_for_uniform( "material.specular",  material.specular);
   shader_p->set_float_for_uniform("material.shininess", material.shininess);
   shader_p->set_float_for_uniform("material.specular_strength", material.specular_strength);

   if (false)
      std::cout << "DEBUG:: Mesh::draw() name: " << name << " shader-name: " << shader_p->name
                << " eye " << static_cast<int>(eye) 
                << " stereo_x_scale " << stereo_x_scale << " stereo_x_offset " << stereo_x_offset << std::endl;

   shader_p->set_float_for_uniform("stereo_x_scale",  stereo_x_scale);
   shader_p->set_float_for_uniform("stereo_x_offset", stereo_x_offset);

   if (false) {
      std::cout << "debug:: Mesh::draw(): " << name << " " << shader_p->name << " do_specular "
                << material.do_specularity << std::endl;
      std::cout << "debug:: Mesh::draw(): " << name << " " << shader_p->name << " material.ambient "
                << glm::to_string(material.ambient) << std::endl;
      std::cout << "debug:: Mesh::draw(): " << name << " " << shader_p->name << " material.ambient "
                << glm::to_string(material.ambient) << std::endl;
      std::cout << "debug:: Mesh::draw(): " << name << " " << shader_p->name << " material.diffuse "
                << glm::to_string(material.diffuse) << std::endl;
      std::cout << "debug:: Mesh::draw(): " << name << " " << shader_p->name << " material.specular "
                << glm::to_string(material.specular) << std::endl;
      std::cout << "debug:: Mesh::draw(): " << name << " " << shader_p->name << " material.shininess "
                << material.shininess << std::endl;
      std::cout << "debug:: Mesh::draw(): " << name << " " << shader_p->name << " sent material.shininess "
                << material.shininess << std::endl;
      std::cout << "debug:: Mesh::draw(): " << name << " " << shader_p->name << " sent material.specular_strength "
                << material.specular_strength << std::endl;
   }

   shader_p->set_float_for_uniform("opacity", opacity);
   if (false)
      std::cout << "sending opacity " << opacity << std::endl;

   err = glGetError();
   if (err) std::cout << "GL ERROR:: draw() " << shader_name << " set float for uniform opacity "
                      << " with GL err " << stringify_error_code(err) << std::endl;

   // this is not useful - eye_position_in_molecule_coordinates_space is what is needed
   // for correct specular reflections.
   glm::vec4 eye_position_v4(eye_position, 1.0);
   if (false)
      std::cout << "sending eye_position " << glm::to_string(eye_position) << std::endl;
   shader_p->set_vec3_for_uniform("eye_position", eye_position);

   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw() \"" << name << "\" \"" << shader_name << "\" post-set eye position "
                      << " with GL err " << err << std::endl;

   shader_p->set_bool_for_uniform("show_shadows", show_just_shadows);

   // bind the vertices and their indices

   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw() " << shader_name << " pre-glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;

   if (vao == VAO_NOT_SET)
      std::cout << "GL ERROR:: You forgot to setup this Mesh: \"" << name << "\" "
                << shader_p->name << std::endl; // Or was not set because no vertices or no triangles.. so crash then

   // if (vao == VAO_NOT_SET) return; // maybe?

   // std::cout << "Mesh::draw() using vao " << vao << std::endl;
   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "error:: Mesh::draw() " << shader_name << " " << name
                      << " glBindVertexArray() vao " << vao << " with GL err " << err << std::endl;

   if (true) { // I doubt that the block is needed. 20210823-PE
      glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
      err = glGetError(); if (err) std::cout << "GL ERROR:: Mesh::draw() glBindBuffer() v "
                                             << err << std::endl;
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
      err = glGetError(); if (err) std::cout << "GL ERROR:: Mesh::draw() glBindBuffer() i "
                                             << err << std::endl;
   }

   glEnableVertexAttribArray(0);
   glEnableVertexAttribArray(1);
   glEnableVertexAttribArray(2);
   if (is_instanced_colours) {
      glBindBuffer(GL_ARRAY_BUFFER, inst_colour_buffer_id);
      err = glGetError();
      if (err)
         std::cout << "GL ERROR:: draw() glBindBuffer() inst col " << err << std::endl;
   }

   if (is_instanced)
      glEnableVertexAttribArray(3);

   err = glGetError();
   if (err) std::cout << "GL ERROR:: draw() glBindBuffer() Mesh::draw() post-vertex arrays "
                      << "shader " << shader_p->name << " error " << err << std::endl;

   if (is_instanced_with_rts_matrix) {
      glBindBuffer(GL_ARRAY_BUFFER, inst_rts_buffer_id);
      glEnableVertexAttribArray(4);
      glEnableVertexAttribArray(5);
      glEnableVertexAttribArray(6);
   }

   err = glGetError();
   if (err) std::cout << "   error draw() " << name << " pre-draw " << err << std::endl;

   // std::cout << "Here with use_blending " << use_blending << std::endl;
   if (use_blending) {
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   }

   if (is_instanced) {

      // If you are here, did you remember to use gtk_gl_area_attach_buffers(GTK_GL_AREA(di.gl_area));
      // before making a new VAO?

      if (false)
         std::cout << "debug:: Mesh::draw() instanced: " << name << " " << shader_p->name
                   << " drawing " << n_verts
                   << " triangle vertices"  << " in " << n_instances << " instances" << std::endl;
      glDrawElementsInstanced(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr, n_instances);
      err = glGetError();
      if (err) std::cout << "GL ERROR:: draw() glDrawElementsInstanced()"
                         << " shader: " << shader_p->name
                         << " vao: " << vao
                         << " n_triangle_verts: " << n_verts
                         << " n_instances: " << n_instances
                         << " with GL err " << err << std::endl;
   } else {

      // If you are here, did you remember to use gtk_gl_area_attach_buffers(GTK_GL_AREA(di.gl_area));
      // before making a new VAO?
      // 20220212-PE Sigh... No... no I didn't (map-as-mesh change)

      if (false)
         std::cout << "debug:: Mesh::draw() " << name << " shader " << shader_p->name
                   << " vao " << vao
                   << " drawing " << n_verts << " vertices from triangles with gl_lines_mode " << gl_lines_mode
                   << std::endl;

      if (draw_as_lines_flag) {

         unsigned int n_verts_for_gl_lines = lines_vertex_indices.size();

         if (false)
            std::cout << "debug:: Mesh::draw() " << name << " shader " << shader_p->name
                      << " vao " << vao
                      << " drawing " << n_verts_for_gl_lines << " indices "
                      << std::endl;

         glDrawElements(GL_LINES, n_verts_for_gl_lines, GL_UNSIGNED_INT, nullptr);
         err = glGetError();
         if (err) std::cout << "GL ERROR:: draw() glDrawElements()"
                            << " of Mesh \"" << name << "\""
                            << " shader: " << shader_p->name
                            << " vao " << vao
                            << " n_verts_for_gl_lines " << n_verts_for_gl_lines
                            << " with GL err " << err << std::endl;
      } else {

         glDrawElements(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr);
         err = glGetError();
         if (err) std::cout << "GL ERROR:: draw() glDrawElements()"
                            << " of Mesh \"" << name << "\""
                            << " shader: " << shader_p->name
                            << " vao " << vao
                            << " n_triangle_verts " << n_verts
                            << " with GL err " << err << std::endl;
      }
   }

   if (use_blending) {
      glDisable(GL_BLEND);
   }
   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   if (is_instanced) glDisableVertexAttribArray(3);
   if (is_instanced) glDisableVertexAttribArray(4);
   if (is_instanced) glDisableVertexAttribArray(5);
   if (is_instanced) glDisableVertexAttribArray(6);
   glUseProgram(0);

}


void
Mesh::draw_with_shadows(Shader *shader_p,
                        const glm::mat4 &mvp,
                        const glm::mat4 &view_rotation_matrix,
                        const std::map<unsigned int, lights_info_t> &lights,
                        const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                        float opacity,
                        const glm::vec4 &background_colour,
                        bool do_depth_fog,
                        const glm::mat4 &light_view_mvp,
                        unsigned int shadow_depthMap,
                        float shadow_strength,
                        unsigned int shadow_softness,
                        bool show_just_shadows) {

   // this is largely a copy of draw() - then shadow params added.

   if (false)
      std::cout << "debug:: Mesh::draw_with_shadows() \"" << name << "\" shader \"" << shader_p->name
                << " shadow-depth-map " << shadow_depthMap << "\" --- start ---"
                << std::endl;

   if (false)
      std::cout << "debug:: Mesh::draw_with_shadows() \"" << name << "\" shader: " << shader_p->name
                << " draw_this_mesh: " << draw_this_mesh << " n-tris:" << triangles.size()
                << " is_instanced " <<  is_instanced
                << " is_instanced_colours " <<  is_instanced_colours
                << " is_instanced_with_rts_matrix " <<  is_instanced_with_rts_matrix
                << std::endl;

   if (! draw_this_mesh) return;

   unsigned int n_triangles = triangles.size();
   GLuint n_verts = 3 * n_triangles;

   // At the moment, I don't think that it's an error to come here if there are no triangles (yet)
   // Just quietly do nothing and return.
   if (n_triangles == 0) return;

   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw() " << name << " " << shader_p->name
                      << " -- start -- " << err << std::endl;
   shader_p->Use();
   const std::string &shader_name = shader_p->name;

   // std::cout << "debug:: Mesh::draw_with_shadows() sending mvp " << glm::to_string(mvp) << std::endl;
   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, glm::value_ptr(mvp));
   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw_with_shadows() shader: " << shader_p->name << " post mvp uniform "
                      << err << std::endl;

   glUniformMatrix4fv(shader_p->view_rotation_uniform_location, 1, GL_FALSE, &view_rotation_matrix[0][0]);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw_with_shadows() " << name << " " << shader_p->name
                      << " draw() post view rotation uniform err " << _(err) << std::endl;

   shader_p->set_mat4_for_uniform("light_space_mvp", light_view_mvp);
   err = glGetError();
   if (err) std::cout << "GL ERROR: TextureMesh::draw_with_shadows(): "
                      << shader_p->name << " post light-space-mvp err " << _(err) << std::endl;

   std::map<unsigned int, lights_info_t>::const_iterator it;
   unsigned int light_idx = 0;
   it = lights.find(light_idx);
   if (it != lights.end())
      shader_p->setup_light(light_idx, it->second, view_rotation_matrix);
   light_idx = 1;
   it = lights.find(light_idx);
   if (it != lights.end())
      shader_p->setup_light(light_idx, it->second, view_rotation_matrix);

   glActiveTexture(GL_TEXTURE0);
   err = glGetError(); if (err) std::cout << "GL ERROR:: Mesh::draw_with_shadows() A4 err " << _(err) << std::endl;
   glBindTexture(GL_TEXTURE_2D, shadow_depthMap);
   shader_p->set_int_for_uniform("shadow_map", 0); // there is only one sampler2D
   err = glGetError(); if (err) std::cout << "GL ERROR:: Mesh::draw_with_shadows() A5 err " << _(err) << std::endl;

   shader_p->set_float_for_uniform("opacity", opacity);

   // add material properties, use class built-ins this time.

   // not using the built-in - hmm.
   shader_p->set_vec4_for_uniform("background_colour", background_colour);

   shader_p->set_bool_for_uniform("do_depth_fog", do_depth_fog);

   bool transferred_colour_is_instanced = true; // pass this?
   shader_p->set_bool_for_uniform("transferred_colour_is_instanced", transferred_colour_is_instanced);

   err = glGetError(); if (err) std::cout << "   error draw() pre-setting material "
                                          << err << std::endl;
   shader_p->set_bool_for_uniform( "do_specular",        material.do_specularity); // change from object in C++ to global unform in GLSL
   shader_p->set_vec4_for_uniform( "material.ambient",   material.ambient);
   shader_p->set_vec4_for_uniform( "material.diffuse",   material.diffuse);
   shader_p->set_vec4_for_uniform( "material.specular",  material.specular);
   shader_p->set_float_for_uniform("material.shininess", material.shininess);
   shader_p->set_float_for_uniform("material.specular_strength", material.specular_strength);

   if (false) {
      std::cout << "debug:: Mesh::draw_with_shadows(): " << name << " " << shader_p->name << " do_specular "
                << material.do_specularity << std::endl;
      std::cout << "debug:: Mesh::draw_with_shadows(): " << name << " " << shader_p->name << " material.ambient "
                << glm::to_string(material.ambient) << std::endl;
      std::cout << "debug:: Mesh::draw_with_shadows(): " << name << " " << shader_p->name << " material.diffuse "
                << glm::to_string(material.diffuse) << std::endl;
      std::cout << "debug:: Mesh::draw_with_shadows(): " << name << " " << shader_p->name << " material.specular "
                << glm::to_string(material.specular) << std::endl;
      std::cout << "debug:: Mesh::draw_with_shadows(): " << name << " " << shader_p->name << " sent material.shininess "
                << material.shininess << std::endl;
      std::cout << "debug:: Mesh::draw_with_shadows(): " << name << " " << shader_p->name << " sent material.specular_strength "
                << material.specular_strength << std::endl;
   }

   err = glGetError();
   if (err) std::cout << "GL ERROR:: draw_with_shadows() " << shader_name << " pre-set eye position "
                      << " with GL err " << err << std::endl;

   if (false)
      std::cout << "sending eye_position " << glm::to_string(eye_position) << std::endl;
   shader_p->set_vec3_for_uniform("eye_position", eye_position);

   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw_with_shadows() \"" << name << "\" \"" << shader_name << "\" post-set eye position"
                      << " with GL err " << err << std::endl;

   shader_p->set_bool_for_uniform("show_shadows", show_just_shadows);
   shader_p->set_float_for_uniform("shadow_strength", shadow_strength);
   shader_p->set_int_for_uniform("shadow_softness", shadow_softness);

   // bind the vertices and their indices

   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw_with_shadows() " << shader_name << " pre-glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;

   if (vao == VAO_NOT_SET)
      std::cout << "ERROR:: You forgot to setup this Mesh \"" << name << "\" "
                << shader_p->name << std::endl; // Or was not set because no vertices or no triangles.. so crash then

   // if (vao == VAO_NOT_SET) return; // maybe?

   // std::cout << "Mesh::draw() using vao " << vao << std::endl;
   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "error:: Mesh::draw_with_shadows() " << shader_name << " " << name
                      << " glBindVertexArray() vao " << vao << " with GL err " << err << std::endl;

   if (true) {
      glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
      err = glGetError(); if (err) std::cout << "GL ERROR:: Mesh::draw_with_shadows() glBindBuffer() v "
                                             << err << std::endl;
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
      err = glGetError(); if (err) std::cout << "GL ERROR:: Mesh::draw_with_shadows() glBindBuffer() i "
                                             << err << std::endl;
   }

   glEnableVertexAttribArray(0);
   glEnableVertexAttribArray(1);
   glEnableVertexAttribArray(2);
   if (is_instanced_colours) {
      glBindBuffer(GL_ARRAY_BUFFER, inst_colour_buffer_id);
      err = glGetError();
      if (err)
         std::cout << "GL ERROR:: draw_with_shadows() glBindBuffer() inst col " << err << std::endl;
   }

   if (is_instanced)
      glEnableVertexAttribArray(3);

   err = glGetError();
   if (err) std::cout << "GL ERROR:: draw_with_shadows() glBindBuffer() Mesh::draw() post-vertex arrays "
                      << "shader " << shader_p->name << " error " << err << std::endl;

   if (is_instanced_with_rts_matrix) {
      glEnableVertexAttribArray(4);
      glEnableVertexAttribArray(5);
      glEnableVertexAttribArray(6);
      glEnableVertexAttribArray(7);  // instanced colour
   }

   err = glGetError();
   if (err) std::cout << "   error draw_with_shadows() " << name << " pre-draw " << err << std::endl;

   if (use_blending) {
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   }

   // std::cout << "is_instanced " << is_instanced << " " << name << std::endl;

   if (is_instanced) {

      // If you are here, did you remember to use gtk_gl_area_attach_buffers(GTK_GL_AREA(di.gl_area));
      // before making a new VAO?

      if (false)
         std::cout << "debug:: Mesh::draw_with_shadows() instanced: " << name << " " << shader_p->name
                   << " drawing " << n_verts
                   << " triangle vertices"  << " in " << n_instances << " instances" << std::endl;

      glDrawElementsInstanced(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr, n_instances);
      err = glGetError();
      if (err) std::cout << "GL ERROR:: draw() glDrawElementsInstanced()"
                         << " shader: " << shader_p->name
                         << " vao: " << vao
                         << " n_triangle_verts: " << n_verts
                         << " n_instances: " << n_instances
                         << " with GL err " << err << std::endl;
   } else {

      // If you are here, did you remember to use gtk_gl_area_attach_buffers(GTK_GL_AREA(di.gl_area));
      // before making a new VAO?

      if (false)
         std::cout << "debug:: Mesh::draw_with_shadows() \"" << name << "\" shader \"" << shader_p->name << "\""
                   << " vao " << vao
                   << " drawing " << n_verts << " triangle vertices"  << std::endl;

      if (gl_lines_mode) {

         unsigned int n_verts_for_gl_lines = lines_vertex_indices.size();
         glDrawElements(GL_LINES, n_verts_for_gl_lines, GL_UNSIGNED_INT, nullptr);
         err = glGetError();
         if (err) std::cout << "GL ERROR:: draw_with_shadows() glDrawElements()"
                            << " of Mesh \"" << name << "\""
                            << " shader: " << shader_p->name
                            << " vao " << vao
                            << " n_verts_for_gl_lines " << n_verts_for_gl_lines
                            << " with GL err " << err << std::endl;
      } else {

         glDrawElements(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr);
         err = glGetError();
         if (err) std::cout << "GL ERROR:: draw_with_shadows() glDrawElements()"
                            << " of Mesh \"" << name << "\""
                            << " shader: " << shader_p->name
                            << " vao " << vao
                            << " n_triangle_verts " << n_verts
                            << " with GL err " << err << std::endl;
      }
   }

   if (use_blending) {
      glDisable(GL_BLEND);
   }
   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   if (is_instanced) glDisableVertexAttribArray(3);
   if (is_instanced) glDisableVertexAttribArray(4);
   if (is_instanced) glDisableVertexAttribArray(5);
   if (is_instanced) glDisableVertexAttribArray(6);
   if (is_instanced) glDisableVertexAttribArray(7);
   glUseProgram(0);

   // std::cout << "Mesh::draw_with_shadows() \"" << name << "\" \"" << shader_p->name << "\" --- done --- " << std::endl;

}




void
Mesh::draw_for_ssao(Shader *shader_p,
                    const glm::mat4 &model,
                    const glm::mat4 &view,
                    const glm::mat4 &projection) {

   if (false)
      std::cout << "debug:: start Mesh::draw_for_ssao() this mesh: " << name << " with shader " << shader_p->name
                << std::endl;

   if (! draw_this_mesh) return;

   if (! shader_p) return; // if we don't want this mesh to be drawn a null shader is passed

   unsigned int n_triangles = triangles.size();
   GLuint n_verts = 3 * n_triangles;
   if (n_triangles == 0) return;

   if (false)
      std::cout << "debug:: in Mesh::draw_for_ssao() " << name << " " << shader_p->name
                << " n_verts " << n_verts << " n_triangles " << n_triangles << std::endl;

   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw_for_ssao() " << shader_p->name << " -- start -- "
                      << err << std::endl;

   shader_p->Use();
   const std::string &shader_name = shader_p->name;

   // vertex Shader:
   //
   // layout(location = 0) in vec3 position;
   // layout(location = 1) in vec3 normal;
   // layout(location = 2) in vec4 colour;

   // uniform mat4 model; // include the view rotation and item translation
   // uniform mat4 view; // the lookat matrix
   // uniform mat4 projection; // projection matrix

   shader_p->set_mat4_for_uniform("model",      model);
   shader_p->set_mat4_for_uniform("view",       view);
   shader_p->set_mat4_for_uniform("projection", projection);

   if (false) {
      glm::mat4 mvp = projection * view * model;
      std::cout << "debug:: in Mesh::draw_for_ssao() model "  << glm::to_string(model) << std::endl;
      std::cout << "debug:: in Mesh::draw_for_ssao() view  "  << glm::to_string(view) << std::endl;
      std::cout << "debug:: in Mesh::draw_for_ssao() proj  "  << glm::to_string(projection) << std::endl;
      std::cout << "debug:: in Mesh::draw_for_ssao() (calc) " << glm::to_string(mvp) << std::endl;
   }

   err = glGetError();
   if (err)
     std::cout << "GL ERROR:: Mesh::draw_for_ssao() " << shader_name << " post uniforms" << std::endl;

   if (vao == VAO_NOT_SET)
      std::cout << "Mesh::draw_for_ssao() You forgot to setup this mesh "
                << "(or setup with empty vertices or triangles) "
                << "\"" << name << "\" \"" << shader_p->name << "\"" << std::endl;

   glBindVertexArray(vao);
   err = glGetError();
   if (err)
     std::cout << "   error draw_for_ssao() \"" << shader_name << "\" \"" << name << "\""
               << " glBindVertexArray() vao " << vao << " with GL err " << err << std::endl;

   glEnableVertexAttribArray(0); // position
   glEnableVertexAttribArray(1); // normal
   glEnableVertexAttribArray(2); // colour

   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw_for_ssao() " << name << " pre-draw " << err << std::endl;

   // 20220804-PE did you forget to attach_buffers() before making this mesh?
   glDrawElements(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw_for_ssao() glDrawElements() of Mesh "
                      << "\"" << name << "\""
                      << " shader: " << shader_p->name
                      << " vao " << vao
                      << " n_triangle_verts " << n_verts
                      << " with GL err " << err << std::endl;

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);

   glUseProgram(0);

}



void
Mesh::draw_instances_for_ssao(Shader *shader_p,
                              const glm::mat4 &model,
                              const glm::mat4 &view,
                              const glm::mat4 &projection) {

   if (! draw_this_mesh) return;
   if (n_instances == 0) return;
   if (triangles.empty()) return;

   if (false)
      std::cout << "debug:: Mesh::draw_instances_for_ssao() " << name << " \"" << shader_p->name << "\""
                << " " << n_instances << std::endl;

   shader_p->Use();
   glBindVertexArray(vao);
   GLenum err = glGetError();
   if (err) std::cout << "error draw_instances() " << shader_p->name
                      << " glBindVertexArray() vao " << vao
                      << " with GL err " << err << std::endl;

   glEnableVertexAttribArray(0); // vertex positions
   glEnableVertexAttribArray(1); // vertex normal
   glEnableVertexAttribArray(2); // vertex colour
   glEnableVertexAttribArray(3); // 4xvec4 for model rotation,scale,translation
   glEnableVertexAttribArray(4);
   glEnableVertexAttribArray(5);
   glEnableVertexAttribArray(6);
   glEnableVertexAttribArray(7); // instanced colour - not used of course

   shader_p->set_mat4_for_uniform("model",      model);
   shader_p->set_mat4_for_uniform("view",       view);
   shader_p->set_mat4_for_uniform("projection", projection);

   GLuint n_verts = 3 * triangles.size();
   glDrawElementsInstanced(GL_TRIANGLES, n_verts, GL_UNSIGNED_INT, nullptr, n_instances);

   glDisableVertexAttribArray(0);
   glDisableVertexAttribArray(1);
   glDisableVertexAttribArray(2);
   glDisableVertexAttribArray(3);
   glDisableVertexAttribArray(4);
   glDisableVertexAttribArray(5);
   glDisableVertexAttribArray(6);
   glDisableVertexAttribArray(7);

}



// draw symmetry with lines
void
Mesh::draw_simple_bond_lines(Shader *shader_p,
                             const glm::mat4 &mvp,
                             const glm::vec4 &background_colour,
                             float line_width,
                             bool do_depth_fog) {

   if (vao == VAO_NOT_SET)
      std::cout << "ERROR:: You forgot to setup this Mesh " << name << " "
                << shader_p->name << std::endl;

   shader_p->Use();
   GLenum err = glGetError();
   if (err) std::cout << "error:: Mesh::draw_symmetry() " << shader_p->name << " " << name
                      << " use shader with GL err " << err << std::endl;

   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::draw_simple_bond_lines() " << shader_p->name << " " << name
                      << " glBindVertexArray() vao " << vao << " with GL err " << err << std::endl;

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   err = glGetError();
   if (err) std::cout << "error:: " << shader_p->name << "GL ERROR:: Mesh::draw_simple_bond_lines() post mvp uniform "
                      << err << std::endl;

   glEnableVertexAttribArray(0);
   glEnableVertexAttribArray(1);

   shader_p->set_vec4_for_uniform("background_colour", background_colour);
   shader_p->set_bool_for_uniform("do_depth_fog", do_depth_fog);

#ifdef __APPLE__
#else
   glLineWidth(line_width); // may not be respected
#endif

   unsigned int n_verts = n_simple_bond_lines * 2;
   unsigned int first = 0;

   glDrawArrays(GL_LINES, first, n_verts); // first and count
   err = glGetError();
   if (err) std::cout << "error:: Mesh::draw_simple_bond_lines() " << shader_p->name << " " << name
                      << " post glDrawArrays() " << vao << " with GL err " << err << std::endl;

   if (false) {
      // why does this do damage? Anyway, the "switch" is done by binding the VAO
      glDisableVertexAttribArray(0);
      glDisableVertexAttribArray(1);
   }
   glBindVertexArray(0);
}

// draw symmetry with lines
void
Mesh::draw_symmetry(Shader *shader_p,
                    const glm::mat4 &mvp,
                    const glm::mat4 &mouse_based_model_rotation_matrix,
                    const std::map<unsigned int, lights_info_t> &lights,
                    const glm::vec3 &eye_position,
                    const glm::vec4 &background_colour,
                    bool do_depth_fog) {

   float line_width = 5.0; // pass this

   if (vao == VAO_NOT_SET)
      std::cout << "ERROR:: You forgot to setup this Mesh " << name << " "
                << shader_p->name << std::endl;

   shader_p->Use();
   GLenum err = glGetError();
   if (err) std::cout << "error:: Mesh::draw_symmetry() " << shader_p->name << " " << name
                      << " use shader with GL err " << err << std::endl;

   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "error:: Mesh::draw_symmetry() " << shader_p->name << " " << name
                      << " glBindVertexArray() vao " << vao << " with GL err " << err << std::endl;

   glUniformMatrix4fv(shader_p->mvp_uniform_location, 1, GL_FALSE, &mvp[0][0]);
   err = glGetError();
   if (err) std::cout << "error:: " << shader_p->name << " Mesh::draw_symmetry() post mvp uniform "
                      << err << std::endl;

   // is this needed?
   // glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   // err = glGetError(); if (err) std::cout << "   error draw() glBindBuffer() " << err << std::endl;

   glEnableVertexAttribArray(0);
   glEnableVertexAttribArray(1);

   shader_p->set_vec4_for_uniform("background_colour", background_colour);
   shader_p->set_bool_for_uniform("do_depth_fog", do_depth_fog);
#ifdef __APPLE__
#else
   glLineWidth(line_width);
#endif
   unsigned int n_verts = n_symmetry_atom_lines_vertices;
   unsigned int first = 0;

   glDrawArrays(GL_LINES, first, n_verts); // first and count
   err = glGetError();
   if (err) std::cout << "error:: Mesh::draw_symmetry() " << shader_p->name << " " << name
                      << " post glDrawArrays() " << vao << " with GL err " << err << std::endl;

   if (false) {
      // why does this do damage? Anyway, the "switch" is done by binding the VAO
      glDisableVertexAttribArray(0);
      glDisableVertexAttribArray(1);
   }
   glBindVertexArray(0);
}

void
Mesh::update_vertices() {

   unsigned int vs = vertices.size();
   if (vs > 0) {
      glBindVertexArray(vao); // needed?
      glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
      glBufferSubData(GL_ARRAY_BUFFER, 0, vs * sizeof(s_generic_vertex), &(vertices[0]));
   }
}

void
Mesh::update_instancing_buffer_data(const std::vector<glm::mat4> &mats,
                                    const std::vector<glm::vec4> &colours) {

   // glBufferSubData(	GLenum        target,
   //                   GLintptr      offset,
   //                   GLsizeiptr    size,
   //                   const GLvoid *data);

   unsigned int n_mats =    mats.size();
   unsigned int n_cols = colours.size();

   //is_instanced = true; 20210910-PE don't do this here. Do it when these buffers are allocated.

   if (vao == VAO_NOT_SET)
      std::cout << "You forgot to setup this Mesh " << name << std::endl;
   glBindVertexArray(vao);

   if (false) {
      std::cout << "-------- update_instancing_buffer_data() mats " << mats.size() << std::endl;
      std::cout << "-------- update_instancing_buffer_data() cols " << colours.size() << std::endl;
   }

   if (n_mats > 0) {
      glBindBuffer(GL_ARRAY_BUFFER, inst_rts_buffer_id);
      glBufferSubData(GL_ARRAY_BUFFER, 0, n_mats * 4 * sizeof(glm::vec4), &(mats[0]));
   }
   if (n_cols > 0) {
      glBindBuffer(GL_ARRAY_BUFFER, inst_colour_buffer_id);
      glBufferSubData(GL_ARRAY_BUFFER, 0, n_cols * sizeof(glm::vec4), &(colours[0]));
   }
}

#if 0
void
Mesh::update_instancing_buffer_data(const std::vector<glm::mat4> &mats) {

   // glBufferSubData(	GLenum        target,
   //                   GLintptr      offset,
   //                   GLsizeiptr    size,
   //                   const GLvoid *data);

   if (vao == VAO_NOT_SET)
      std::cout << "You forgot to setup this Mesh " << name << std::endl;
   glBindVertexArray(vao); // needed? 20210910-PE I imagine so!
   GLenum err = glGetError();
   if (err)
      std::cout << "GL error Mesh::update_instancing_buffer_data() --start-- " << "binding vao " << vao
                << " error " << _(err) << std::endl;

   int n_mats =    mats.size();
   if (n_mats > n_instances_allocated) {
      std::vector<glm::vec4> dummy;
      std::cout << "::::::::::::: debug:: update_instancing_buffer_data(mats) calls setup_matrix_and_colour_instancing_buffers()"
                << std::endl;
      setup_matrix_and_colour_instancing_buffers_standard(mats, dummy);
      std::cout << "::::::::::::: debug:: update_instancing_buffer_data(mats) returned from setup_matrix_and_colour_instancing_buffers()"
                << std::endl;
   }

   if (n_mats > 0) {
      for (unsigned int i=0; i<mats.size(); i++) {
         std::cout << "debug:: update_instancing_buffer_data()" << inst_rts_buffer_id << " "
                   << glm::to_string(mats[i]) << std::endl;
      }
      glBindBuffer(GL_ARRAY_BUFFER, inst_rts_buffer_id);
      glBufferSubData(GL_ARRAY_BUFFER, 0, n_mats * 4 * sizeof(glm::vec4), &(mats[0]));
   }
}
#endif

void
Mesh::update_instancing_buffer_data_standard(const std::vector<glm::mat4> &mats) {

   // glBufferSubData(	GLenum        target,
   //                   GLintptr      offset,
   //                   GLsizeiptr    size,
   //                   const GLvoid *data);

   GLenum err = glGetError();
   if (err) std::cout << "GL Error Mesh::update_instancing_buffer_data_standard() --start-- error: " << err << std::endl;

   if (vao == VAO_NOT_SET)
      std::cout << "You forgot to setup this Mesh " << name << std::endl;

   glBindVertexArray(vao);
   err = glGetError();
   if (err)
      std::cout << "GL error Mesh::update_instancing_buffer_data_standard() A1 "
                << "binding vao " << vao << " error " << _(err) << std::endl;
   if (err == GL_INVALID_OPERATION)
      std::cout << "Because vao was not the name of a vertex array object previously returned from a call to glGenVertexArrays (or zero)"
                << std::endl;

   int n_mats = mats.size();
   if (n_mats > n_instances_allocated) {
      std::vector<glm::vec4> dummy;
      std::cout << "::::::::::::: debug:: update_instancing_buffer_data_standard(mats) calls setup_matrix_and_colour_instancing_buffers_standard()"
                << std::endl;
      setup_matrix_and_colour_instancing_buffers_standard(mats, dummy);
      std::cout << "::::::::::::: debug:: update_instancing_buffer_data(mats) returned from setup_matrix_and_colour_instancing_buffers()"
                << std::endl;
      // if (n_mats > 1000) n_mats = 1000;
   }

   if (n_mats > 0) {
      glBindBuffer(GL_ARRAY_BUFFER, inst_rts_buffer_id);
      glBufferSubData(GL_ARRAY_BUFFER, 0, n_mats * 4 * sizeof(glm::vec4), &(mats[0]));
   }
}

void
Mesh::setup_instancing_buffer_data_for_extra_distance_restraints(unsigned int n_matrices) {

   GLenum err = glGetError();
   if (err) std::cout << "Error setup_matrix_and_colour_instancing_buffers_standard() -- start -- "
                      << err << std::endl;

   n_instances = n_matrices;
   n_instances_allocated = n_instances;

   // these vectors should be the same size. Add a check.

   err = glGetError();
   if (err) std::cout << "error setup_instancing_buffer_data_for_extra_distance_restraints() A "
                      << err << std::endl;

   if (vao == VAO_NOT_SET)
      std::cout << "ERROR:: in setup_instancing_buffer_data_for_extra_distance_restraints() You didn't correctly setup this Mesh "
                << name << " " << std::endl;

   glBindVertexArray(vao);

   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::setup_instancing_buffer_data_for_extra_distance_restraints() B binding-vao "
                      << err << " with vao " << vao << std::endl;

   if (! first_time) {
      // std::cout << "in setup_instancing_buffer_data_for_extra_distance_restraints() deleting inst_rts_buffer_id " << std::endl;
      glDeleteBuffers(1, &inst_rts_buffer_id); // setup_buffers() sets first_time to false but doesn't set inst_rts_buffer_id.
   }
   glGenBuffers(1, &inst_rts_buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, inst_rts_buffer_id);
   unsigned int size_of_edrmidt = sizeof(extra_distance_restraint_markup_instancing_data_t);
   if (false)
      std::cout << "Mesh::setup_instancing_buffer_data_for_extra_distance_restraints() allocating matrix buffer data "
                << n_instances * size_of_edrmidt << std::endl;
   glBufferData(GL_ARRAY_BUFFER, n_instances * size_of_edrmidt, nullptr, GL_DYNAMIC_DRAW); // dynamic

   // 3 float width
   // 4 float length
   // 5 vec3 position
   // 6 mat3 ori 0
   // 7 mat3 ori 1
   // 8 mat3 ori 2
   // 9 vec4 colour inst
   err = glGetError(); if (err) std::cout << "   ERROR setup_instancing_buffer_data_for_extra_distance_restraints() C0 " << err << std::endl;
   glEnableVertexAttribArray(3);
   glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, size_of_edrmidt, 0);
   glVertexAttribDivisor(3, 1);
   err = glGetError(); if (err) std::cout << "   error setup_instancing_buffer_data_for_extra_distance_restraints() C3 " << err << std::endl;
   glEnableVertexAttribArray(4);
   glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, size_of_edrmidt, reinterpret_cast<void *>(sizeof(float)));
   glVertexAttribDivisor(4, 1);
   err = glGetError(); if (err) std::cout << "   error setup_instancing_buffer_data_for_extra_distance_restraints() C4 " << err << std::endl;

   glEnableVertexAttribArray(5);
   glVertexAttribPointer(5, 3, GL_FLOAT, GL_FALSE, size_of_edrmidt, reinterpret_cast<void *>(2 * sizeof(float)));
   glVertexAttribDivisor(5, 1);
   err = glGetError(); if (err) std::cout << "   error setup_instancing_buffer_data_for_extra_distance_restraints() C5 " << err << std::endl;

   // orientation 3 x vec3
   glEnableVertexAttribArray(6);
   glVertexAttribPointer(6, 3, GL_FLOAT, GL_FALSE, size_of_edrmidt, reinterpret_cast<void *>(2 * sizeof(float) + sizeof(glm::vec3)));
   glVertexAttribDivisor(6, 1);
   err = glGetError(); if (err) std::cout << "   error setup_instancing_buffer_data_for_extra_distance_restraints() C6 " << err << std::endl;
   glEnableVertexAttribArray(7);
   glVertexAttribPointer(7, 3, GL_FLOAT, GL_FALSE, size_of_edrmidt, reinterpret_cast<void *>(2 * sizeof(float) + 2 * sizeof(glm::vec3)));
   glVertexAttribDivisor(7, 1);
   err = glGetError(); if (err) std::cout << "   error setup_instancing_buffer_data_for_extra_distance_restraints() C7 " << err << std::endl;
   glEnableVertexAttribArray(8);
   glVertexAttribPointer(8, 3, GL_FLOAT, GL_FALSE, size_of_edrmidt, reinterpret_cast<void *>(2 * sizeof(float) + 3 * sizeof(glm::vec3)));
   glVertexAttribDivisor(8, 1);
   err = glGetError(); if (err) std::cout << "   error setup_instancing_buffer_data_for_extra_distance_restraints() C8 " << err << std::endl;

   // this is the colour
   glEnableVertexAttribArray(9);
   glVertexAttribPointer(9, 4, GL_FLOAT, GL_FALSE, size_of_edrmidt, reinterpret_cast<void *>(2 * sizeof(float) + 4 * sizeof(glm::vec3)));
   glVertexAttribDivisor(9, 1);
   err = glGetError(); if (err) std::cout << "   error setup_instancing_buffer_data_for_extra_distance_restraints() C9 " << err << std::endl;

}

//! GM restraints are not just orienation an colour there is a width and length too
//! (hmm - I mean, maybe that can be captured in the mat4).
void
Mesh::update_instancing_buffer_data_for_extra_distance_restraints(const std::vector<extra_distance_restraint_markup_instancing_data_t> &edrmid) {

   GLenum err = glGetError();
   if (err) std::cout << "GL Error Mesh::update_instancing_buffer_data_standard() --start-- error: " << err << std::endl;

   if (vao == VAO_NOT_SET)
      std::cout << "You forgot to setup this Mesh " << name << std::endl;

   glBindVertexArray(vao);
   err = glGetError();
   if (err)
      std::cout << "GL error Mesh::update_instancing_buffer_data_standard() A1 "
                << "binding vao " << vao << " error " << _(err) << std::endl;
   if (err == GL_INVALID_OPERATION)
      std::cout << "Because vao was not the name of a vertex array object previously returned from a call to glGenVertexArrays (or zero)"
                << std::endl;

   unsigned int size_of_edrmidt = sizeof(extra_distance_restraint_markup_instancing_data_t);

   int n_mats = edrmid.size();

   if (n_mats > n_instances_allocated) n_mats = n_instances_allocated;

   if (n_mats > 0) {
      glBindBuffer(GL_ARRAY_BUFFER, inst_rts_buffer_id);
      glBufferSubData(GL_ARRAY_BUFFER, 0, n_mats * size_of_edrmidt, &(edrmid[0]));
      n_instances = n_mats;
   }

}

// static
std::string
Mesh::_(int err) {

   std::string s = std::to_string(err);
   if (err == GL_INVALID_ENUM)      s = "GL_INVALID_ENUM";
   if (err == GL_INVALID_OPERATION) s = "GL_INVALID_OPERATION";
   if (err == GL_INVALID_VALUE)     s = "GL_INVALID_VALUE";
   return s;
}

void
Mesh::update_instancing_buffer_data_for_particles(const particle_container_t &particles) {

   if (false) {
      std::cout << "debug:: update_instancing_buffer_data_for_particles()" << std::endl;
      for (unsigned int i=0; i<particles.size(); i++) {
         std::cout << "    " << i << " "
                   << glm::to_string(particles.particles[i].position) << " "
                   << glm::to_string(particles.particles[i].velocity) << " "
                   << glm::to_string(particles.particles[i].colour)   << std::endl;
      }
   }

   is_instanced = true;
   is_instanced_colours = true;

   if (false) // debugging
      std::cout << "DEBUG:: Mesh::update_instancing_buffer_data_for_particles() " << name
                << " --- start --- " << std::endl;

   GLenum err = glGetError();
   if (err) {
      std::cout << "GL ERROR:: Mesh::update_instancing_buffer_data_for_particles() A0 "
		<< "--- start --- " << _(err) << std::endl;
      logger.log(log_t::GL_ERROR, "Mesh::update_instancing_buffer_data_for_particles() A0 ",
		 "--- start --- ", _(err));
   }

   if (vao == VAO_NOT_SET) {
      std::cout << "GL ERROR:: You forgot to setup this Mesh " << name << std::endl;
      logger.log(log_t::GL_ERROR, "You forgot to setup this Mesh ", name);
   }

   glBindVertexArray(vao);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::update_instancing_buffer_data_for_particles() A1 "
                      << "binding vao " << vao << " " << _(err) << std::endl;

   n_instances = particles.size();
   if (n_instances > n_instances_allocated) {
      std::cout << "OOPPS! Too many particles! " << n_instances << " " << n_instances_allocated << std::endl;
      n_instances = n_instances_allocated;
   }

   if (n_instances > 0) {

      if (false)
         std::cout << "DEBUG:: update_instancing_buffer_data_for_particles() transfering " << n_instances
                   << " particle/instances " << std::endl;

      glBindBuffer(GL_ARRAY_BUFFER, inst_model_translation_buffer_id);
      err = glGetError();
      if (err) std::cout << "GL ERROR:: Mesh::update_instancing_buffer_data_for_particles() A3 "
                         << " vao " << vao
                         << " inst_model_translation_buffer_id " << inst_model_translation_buffer_id
                         << "\n";
      glBufferSubData(GL_ARRAY_BUFFER, 0, n_instances * sizeof(Particle), &(particles.particles[0]));
      err = glGetError();
      if (err) std::cout << "GL ERROR:: Mesh::update_instancing_buffer_data_for_particles() B " << _(err) << "\n";
      glBindBuffer(GL_ARRAY_BUFFER, inst_colour_buffer_id);
      err = glGetError();
      if (err) std::cout << "GL ERROR:: Mesh::update_instancing_buffer_data_for_particles() C\n";
      glBufferSubData(GL_ARRAY_BUFFER, 0, n_instances * sizeof(Particle), &(particles.particles[0]));
      err = glGetError();
      if (err) std::cout << "GL ERROR:: Mesh::update_instancing_buffer_data_for_particles() D " << _(err) << "\n";
   }

   if (false)
      std::cout << "DEBUG:: Mesh::update_instancing_buffer_data_for_particles() " << name
                << " --- returning --- " << std::endl;

}

void
Mesh::smooth_triangles() {

   // What other vertices have the same position as this one?
   // std::map<unsigned int, std::set<unsigned int> >

   // I need to test if aa.pos == bb.pos

   // If there are lots of triangles, that will be slow
   // so instead test aa_hash == bb_hash
   // aa_hash = make_has(aa.pos)
   // which needs a hash map
   // and a glm::vec3 -> hash for position

   // but for now, do it the simple-minded way

   int n = vertices.size();
   std::map<unsigned int, std::set<unsigned int> > m;
   for (unsigned int i=0; i<triangles.size(); i++) {
      for (unsigned int j=i; j<triangles.size(); j++) {
         if (i != j) {
            // raw compare glm::vec3
            if (vertices[triangles[i].point_id[0]].pos == vertices[triangles[j].point_id[0]].pos) {
               m[triangles[i].point_id[0]].insert(triangles[j].point_id[0]);
               m[triangles[j].point_id[0]].insert(triangles[i].point_id[0]);
            }
            if (vertices[triangles[i].point_id[0]].pos == vertices[triangles[j].point_id[1]].pos) {
               m[triangles[i].point_id[0]].insert(triangles[j].point_id[1]);
               m[triangles[j].point_id[1]].insert(triangles[i].point_id[0]);
            }
            if (vertices[triangles[i].point_id[0]].pos == vertices[triangles[j].point_id[2]].pos) {
               m[triangles[i].point_id[0]].insert(triangles[j].point_id[2]);
               m[triangles[j].point_id[2]].insert(triangles[i].point_id[0]);
            }
            if (vertices[triangles[i].point_id[1]].pos == vertices[triangles[j].point_id[0]].pos) {
               m[triangles[i].point_id[1]].insert(triangles[j].point_id[0]);
               m[triangles[j].point_id[0]].insert(triangles[i].point_id[1]);
            }
            if (vertices[triangles[i].point_id[1]].pos == vertices[triangles[j].point_id[1]].pos) {
               m[triangles[i].point_id[1]].insert(triangles[j].point_id[1]);
               m[triangles[j].point_id[1]].insert(triangles[i].point_id[1]);
            }
            if (vertices[triangles[i].point_id[1]].pos == vertices[triangles[j].point_id[2]].pos) {
               m[triangles[i].point_id[1]].insert(triangles[j].point_id[2]);
               m[triangles[j].point_id[2]].insert(triangles[i].point_id[1]);
            }
            if (vertices[triangles[i].point_id[2]].pos == vertices[triangles[j].point_id[0]].pos) {
               m[triangles[i].point_id[2]].insert(triangles[j].point_id[0]);
               m[triangles[j].point_id[0]].insert(triangles[i].point_id[2]);
            }
            if (vertices[triangles[i].point_id[2]].pos == vertices[triangles[j].point_id[1]].pos) {
               m[triangles[i].point_id[2]].insert(triangles[j].point_id[1]);
               m[triangles[j].point_id[1]].insert(triangles[i].point_id[2]);
            }
            if (vertices[triangles[i].point_id[2]].pos == vertices[triangles[j].point_id[2]].pos) {
               m[triangles[i].point_id[2]].insert(triangles[j].point_id[2]);
               m[triangles[j].point_id[2]].insert(triangles[i].point_id[2]);
            }

         }
      }
   }

   float cos_min = cosf(45.0 * 2.0 * M_PI);
   std::map<unsigned int, std::set<unsigned int> >::const_iterator it;
   for (it=m.begin(); it!=m.end(); ++it) {
      glm::vec3 sum = vertices[it->first].pos;
      glm::vec3 base_normal = vertices[it->first].normal;
      unsigned int count = 1;
      std::set<unsigned int>::const_iterator it_s;
      for (it_s=it->second.begin(); it_s!=it->second.end(); ++it_s) {
         float dp = glm::dot(base_normal, vertices[*it_s].normal);
         if (dp > cos_min) {
            sum += vertices[*it_s].pos;
            count += 1;
         }
      }
      float sc = 1.0f / static_cast<float>(count);
      glm::vec3 av = sc * sum;
      vertices[it->first].pos = av;
   }

}

#include <glm/gtx/rotate_vector.hpp> // for orientation

void
Mesh::setup_camera_facing_outline() {

   if (true) {

      unsigned int idx_base = vertices.size();
      add_one_origin_cylinder();
      float half = 0.5;
      float angle = 0.5f * M_PI;
      for (unsigned int i=idx_base; i<vertices.size(); i++) {
         vertices[i].pos.x *= 0.03f;
         vertices[i].pos.y *= 0.03f;
         vertices[i].pos    = glm::rotate(vertices[i].pos,    angle, glm::vec3(1,0,0));
         vertices[i].normal = glm::rotate(vertices[i].normal, angle, glm::vec3(1,0,0));
         vertices[i].pos.x += half;
         vertices[i].pos.y += half;
         vertices[i].color = glm::vec4(0.3f, 0.4f, 0.5f, 1.0f);
      }

      idx_base = vertices.size();
      add_one_origin_cylinder();
      for (unsigned int i=idx_base; i<vertices.size(); i++) {
         vertices[i].pos.x *= 0.03f;
         vertices[i].pos.y *= 0.03f;
         vertices[i].pos    = glm::rotate(vertices[i].pos,    angle, glm::vec3(1,0,0));
         vertices[i].normal = glm::rotate(vertices[i].normal, angle, glm::vec3(1,0,0));
         vertices[i].pos.x -= half;
         vertices[i].pos.y += half;
         vertices[i].color = glm::vec4(0.3f, 0.4f, 0.5f, 1.0f);
      }

      idx_base = vertices.size();
      add_one_origin_cylinder();
      for (unsigned int i=idx_base; i<vertices.size(); i++) {
         vertices[i].pos.x *= 0.03f;
         vertices[i].pos.y *= 0.03f;
         vertices[i].pos    = glm::rotate(vertices[i].pos,    angle, glm::vec3(0,1,0));
         vertices[i].normal = glm::rotate(vertices[i].normal, angle, glm::vec3(0,1,0));
         vertices[i].pos.x += -half;
         vertices[i].pos.y += -half;
         vertices[i].color = glm::vec4(0.3f, 0.4f, 0.5f, 1.0f);
      }

      idx_base = vertices.size();
      add_one_origin_cylinder();
      for (unsigned int i=idx_base; i<vertices.size(); i++) {
         vertices[i].pos.x *= 0.03f;
         vertices[i].pos.y *= 0.03f;
         vertices[i].pos    = glm::rotate(vertices[i].pos,    angle, glm::vec3(0,1,0));
         vertices[i].normal = glm::rotate(vertices[i].normal, angle, glm::vec3(0,1,0));
         vertices[i].pos.x += -half;
         vertices[i].pos.y +=  half;
         vertices[i].color = glm::vec4(0.3f, 0.4f, 0.5f, 1.0f);
      }
   }
   setup_buffers();
}

void
Mesh::setup_camera_facing_quad() {

   float scale_x = 0.4; // pass?
   float scale_y = 0.2;

   glm::vec3 n(0,0,1);
   glm::vec4 col(1.0, 1.0, 1.0, 1.0);

   vertices.clear();
   triangles.clear();

   vertices.push_back(s_generic_vertex(glm::vec3( scale_x,  scale_y, 0.0f), n, col));
   vertices.push_back(s_generic_vertex(glm::vec3(-scale_x,  scale_y, 0.0f), n, col));
   vertices.push_back(s_generic_vertex(glm::vec3(-scale_x, -scale_y, 0.0f), n, col));
   vertices.push_back(s_generic_vertex(glm::vec3( scale_x, -scale_y, 0.0f), n, col));

   triangles.push_back(g_triangle(0,1,2));
   triangles.push_back(g_triangle(2,3,0));

   setup_buffers();

}


void
Mesh::setup_camera_facing_hex() {

   // for hexagonal coloured particles

   float s = 0.1;
   glm::vec3 n(0,0,1);
   glm::vec4 cc(0.4, 0.4, 0.4, 1.0);
   glm::vec4  c(0.4, 0.4, 0.4, 0.1);
   float ot = 0.5;
   float tt = 0.7;
   s_generic_vertex g0(s * glm::vec3(  0,   0, 0), n, cc);
   s_generic_vertex g1(s * glm::vec3( tt,  ot, 0), n, c);
   s_generic_vertex g2(s * glm::vec3( tt, -ot, 0), n, c);
   s_generic_vertex g3(s * glm::vec3(  0,  -1, 0), n, c);
   s_generic_vertex g4(s * glm::vec3(-tt, -ot, 0), n, c);
   s_generic_vertex g5(s * glm::vec3(-tt,  ot, 0), n, c);
   s_generic_vertex g6(s * glm::vec3(  0,   1, 0), n, c);
   unsigned int idx_base = vertices.size();
   unsigned int idx_tri_base = triangles.size();
   vertices.push_back(g0);
   vertices.push_back(g1);
   vertices.push_back(g2);
   vertices.push_back(g3);
   vertices.push_back(g4);
   vertices.push_back(g5);
   vertices.push_back(g6);
   triangles.push_back(g_triangle(0,1,2));
   triangles.push_back(g_triangle(0,2,3));
   triangles.push_back(g_triangle(0,3,4));
   triangles.push_back(g_triangle(0,4,5));
   triangles.push_back(g_triangle(0,5,6));
   triangles.push_back(g_triangle(0,6,1));
   if (idx_tri_base != 0)
      for (unsigned int i=idx_tri_base; i<triangles.size(); i++)
         triangles[i].rebase(idx_base);

   // setup_buffers();
 
}

void
Mesh::setup_camera_facing_polygon(unsigned int n_sides, float scale, bool do_stellation, float stellation_factor) {

   float turn_per_step = 2.0f * M_PI / static_cast<float>(n_sides);
   glm::vec3 n(0,0,1);
   glm::vec4 ccol(1.0, 1.0, 1.0, 1.00);
   glm::vec4  col(0.4, 0.4, 0.4, 0.951);

   unsigned int idx_base = vertices.size();
   unsigned int idx_tri_base = triangles.size();
   vertices.push_back(s_generic_vertex(glm::vec3(0.0f, 0.0f, 0.0f), n, ccol));

   if (do_stellation) {
      for (unsigned int i=0; i<n_sides; i++) {
         float a1 = static_cast<float>(i) * turn_per_step;
         float a2 = (static_cast<float>(i) + 0.5f) * turn_per_step;
         float s1 = sinf(a1);
         float c1 = cosf(a1);
         float s2 = sinf(a2);
         float c2 = cosf(a2);
         glm::vec3 v1 =                     scale * glm::vec3(s1, c1, 0.0f);
         glm::vec3 v2 = stellation_factor * scale * glm::vec3(s2, c2, 0.0f);
         vertices.push_back(s_generic_vertex(v1, n, col));
         vertices.push_back(s_generic_vertex(v2, n, col));
      }
      for (unsigned int idx=0; idx<n_sides; idx++) {
         unsigned int idx_this                = 2 * idx + 1;
         unsigned int idx_for_in_vertex       = 2 * idx + 2;
         unsigned int idx_for_next_out_vertex = 2 * idx + 3;
         if (idx_for_next_out_vertex == 2 * n_sides + 1) idx_for_next_out_vertex = 1;
         triangles.push_back(g_triangle(0, idx_this, idx_for_in_vertex));
         triangles.push_back(g_triangle(0, idx_for_in_vertex, idx_for_next_out_vertex));
      }

   } else {
      for (unsigned int i=0; i<n_sides; i++) {
         float a = static_cast<float>(i) * turn_per_step;
         float s = sinf(a);
         float c = cosf(a);
         vertices.push_back(s_generic_vertex(scale * glm::vec3(s, c, 0.0f), n, col));
      }
      for (unsigned int idx=1; idx<=n_sides; idx++) {
         unsigned int idx_next = idx + 1;
         if (idx == n_sides) idx_next = 1;
         triangles.push_back(g_triangle(0, idx, idx_next));
      }
      if (idx_tri_base != 0)
         for (unsigned int i=idx_tri_base; i<triangles.size(); i++)
            triangles[i].rebase(idx_base);
   }

   setup_buffers(); // this was not here - why not?

}


#include "coot-utils/cylinder-with-rotation-translation.hh"

void
Mesh::setup_hydrogen_bond_cyclinders(Shader *shader_p, const Material &material_in) {

   // call this
   // gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0]));
   // before calling this function

   // Make some space (allocate the buffers) for the rotatation/translation instances also.
   // (the colours are not instanced)

   is_instanced_colours = false; // keep the original colours
   is_instanced = true;
   is_instanced_with_rts_matrix = true;

   material = material_in;
   shader_p->Use();

   unsigned int n_slices = 20;
   unsigned int n_stacks = 80;
   glm::vec3 start_pos(0, 0, 0);
   glm::vec3 end_pos(1, 0, 0); // along the x axis - hmm!
   std::pair<glm::vec3, glm::vec3> pp(start_pos, end_pos);
   float height = glm::distance(start_pos, end_pos);
   float radius = 0.03;
   radius = 0.1;
   cylinder_with_rotation_translation c(pp, radius, radius, height, n_slices, n_stacks);
   c.add_spiral(); // take 2 colours

   // now convert the vertices
   std::vector<s_generic_vertex> new_vertices(c.vertices.size());
   for (unsigned int ii=0; ii<c.vertices.size(); ii++) {
      const coot::api::vertex_with_rotation_translation &v = c.vertices[ii];
      s_generic_vertex gv(v.pos, v.normal, v.colour);
      // gv.color = glm::vec4(0,1,0,1);
      new_vertices[ii] = gv;
   }
   unsigned int idx_base = vertices.size();
   unsigned int idx_tri_base = triangles.size();
   vertices.insert(vertices.end(), new_vertices.begin(), new_vertices.end());
   triangles.insert(triangles.end(), c.triangle_indices_vec.begin(), c.triangle_indices_vec.end());
   for (unsigned int ii=idx_tri_base; ii<triangles.size(); ii++)
      triangles[ii].rebase(idx_base);

   // std::cout << ":::::::::::::::: setup_hydrogen_bond_cyclinders() calls setup_buffers()" << std::endl;
   setup_buffers();
   // std::cout << ":::::::::::::::: setup_hydrogen_bond_cyclinders() returns from setup_buffers()" << std::endl;

   std::vector<glm::mat4> mats(1000, glm::mat4(1.0f));
   std::vector<glm::vec4> colours; //dummy
   n_instances = mats.size();
   setup_matrix_and_colour_instancing_buffers_standard(mats, colours);
}

void
Mesh::setup_extra_distance_restraint_cylinder(const Material &material_in) { // make a cylinder for instancing

   GLenum err = glGetError();
   if (err) {
      std::cout << "GL ERROR:: Mesh::setup_extra_distance_restraint_cylinder() \""
                << name << "\" --- start --- "
                << stringify_error_code(err) << std::endl;
      err = glGetError();
      if (err != 0)
         std::cout << "GL ERROR:: Mesh::setup_extra_distance_restraint_cylinder() \""
                   << name << "\" --- start --- stack-clear "
                   << stringify_error_code(err) << std::endl;
   }

   auto vnc_vertex_to_generic_vertex = [] (const coot::api::vnc_vertex &v) {
      return s_generic_vertex(v.pos, v.normal, v.color);
   };

   auto vnc_vertex_vector_to_generic_vertex_vector = [vnc_vertex_to_generic_vertex] (const std::vector<coot::api::vnc_vertex> &vv) {
      std::vector<s_generic_vertex> vo(vv.size());
      for (unsigned int i=0; i<vv.size(); i++)
         vo[i] = vnc_vertex_to_generic_vertex(vv[i]);
      return vo;
   };

   material = material_in;

   is_instanced = true;
   is_instanced_with_rts_matrix = false;

   unsigned int n_slices = 8; // 32; is too many
   unsigned int n_stacks = 2;

   glm::vec3 start_pos(0, 0, 0);
   glm::vec3 end_pos(0, 0, 1);
   std::pair<glm::vec3, glm::vec3> pp(start_pos, end_pos);
   float height = 1.0;
   float radius = 1.0;
   cylinder c(pp, radius, radius, height, n_slices, n_stacks);

   std::vector<s_generic_vertex> new_vertices(c.vertices.size());
   for (unsigned int ii=0; ii<c.vertices.size(); ii++) {
      const coot::api::vnc_vertex &v = c.vertices[ii];
      new_vertices[ii] = vnc_vertex_to_generic_vertex(v);
   }
   unsigned int idx_base = vertices.size();
   unsigned int idx_tri_base = triangles.size();
   vertices.insert(vertices.end(), new_vertices.begin(), new_vertices.end());
   triangles.insert(triangles.end(), c.triangles.begin(), c.triangles.end());
   for (unsigned int ii=idx_tri_base; ii<triangles.size(); ii++)
      triangles[ii].rebase(idx_base);

   setup_buffers();
}



void
Mesh::test_cyclinders(Shader *shader_p, const Material &material_in) {

   // call this
   // gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0]));
   // before calling this function

   is_instanced_colours = false; // keep the original colours
   is_instanced = true;
   is_instanced_with_rts_matrix = true;

   material = material_in;
   shader_p->Use();

   unsigned int n_slices = 20;
   unsigned int n_stacks = 80;
   glm::vec3 start_pos(0, 0, 0);
   glm::vec3 end_pos(1, 0, 0);
   std::pair<glm::vec3, glm::vec3> pp(start_pos, end_pos);
   float height = glm::distance(start_pos, end_pos);
   float radius = 0.05;
   cylinder_with_rotation_translation c(pp, radius, radius, height, n_slices, n_stacks);
   if (true)
      c.add_spiral(); // take 2 colours

   // now convert the vertices
   std::vector<s_generic_vertex> new_vertices(c.vertices.size());
   for (unsigned int ii=0; ii<c.vertices.size(); ii++) {
      const coot::api::vertex_with_rotation_translation &v = c.vertices[ii];
      s_generic_vertex gv(v.pos, v.normal, v.colour);
      // gv.color = glm::vec4(0,1,0,1);
      new_vertices[ii] = gv;
   }
   unsigned int idx_base = vertices.size();
   unsigned int idx_tri_base = triangles.size();
   vertices.insert(vertices.end(), new_vertices.begin(), new_vertices.end());
   triangles.insert(triangles.end(), c.triangle_indices_vec.begin(), c.triangle_indices_vec.end());
   for (unsigned int ii=idx_tri_base; ii<triangles.size(); ii++)
      triangles[ii].rebase(idx_base);


   setup_buffers();

   // some test cylinders
   std::vector<glm::mat4> mats;
   std::vector<glm::vec4> colours; // empty

   for (auto i=0; i<10; i++) {
      float x = static_cast<float>(i);
      glm::mat4 m(1.0);
      glm::vec3 t(2 * x, 0, 0);
      m += glm::translate(t);
      mats.push_back(m);

      // see get_bond_matrix() in update_mats_and_colours() in the generator
   }

   glm::vec3 p1(42.08, 9.67, 14.42);
   glm::vec3 p2(40.59, 5.68, 13.24);

   glm::vec3 p3(44.88, 12.95, 8.76);
   glm::vec3 p4(46.13, 10.59, 9.97);

   float theta = 1.0;

   mats.push_back(make_hydrogen_bond_cylinder_orientation(p1, p2, theta));
   mats.push_back(make_hydrogen_bond_cylinder_orientation(p3, p4, theta));

   n_instances = mats.size();

   std::cout << "::::::::::::: debug:: test_cyclinders() calls setup_matrix_and_colour_instancing_buffers"
             << std::endl;

   setup_matrix_and_colour_instancing_buffers(mats, colours);


}

// static
glm::mat4
Mesh::make_hydrogen_bond_cylinder_orientation(const glm::vec3 &p1, const glm::vec3 &p2, float theta) {

   glm::mat4 u(1.0f);
   glm::vec3 delta(p2-p1);
   float h = glm::distance(p1, p2);
   glm::mat4 sc = glm::scale(u, glm::vec3(1.0f, 1.0f, h));
   glm::mat4 rot = glm::rotate(u, theta, glm::vec3(0.0f, 0.0f, 1.0f));
   glm::vec3 normalized_bond_orientation(glm::normalize(delta));
   glm::mat4 ori = glm::orientation(normalized_bond_orientation,
                                    glm::vec3(0.0f, 0.0f, 1.0f));
   glm::mat4 t = glm::translate(u, p1);
   glm::mat4 m = t * ori * rot * sc;
   return m;
}

void
Mesh::add_dashed_line(const coot::simple_distance_object_t &sdo, const Material &material,
                      const glm::vec4 &colour) {

   auto coord_orth_to_glm = [] (const clipper::Coord_orth &co) {
                               return glm::vec3(co.x(), co.y(), co.z());
                            };

   auto vnc_vertex_to_generic_vertex = [] (const coot::api::vnc_vertex &v) {
      return s_generic_vertex(v.pos, v.normal, v.color);
   };

   auto vnc_vertex_vector_to_generic_vertex_vector = [vnc_vertex_to_generic_vertex] (const std::vector<coot::api::vnc_vertex> &vv) {
      std::vector<s_generic_vertex> vo(vv.size());
      for (unsigned int i=0; i<vv.size(); i++)
         vo[i] = vnc_vertex_to_generic_vertex(vv[i]);
      return vo;
   };

   double l = sdo.length();
   unsigned int n_seg = static_cast<unsigned int>(l) * 3;
   if (n_seg < 3) n_seg = 3;

   double segment_length = l/(2.0 * static_cast<double>(n_seg));
   clipper::Coord_orth uv = sdo.uv();

   for (unsigned int iseg=0; iseg<n_seg; iseg++) {
      float frac_path_pos_1 = (0.5f + static_cast<float>(2 * iseg)) * segment_length / l;
      float frac_path_pos_2 = (1.5f + static_cast<float>(2 * iseg)) * segment_length / l;
      clipper::Coord_orth pos_1(sdo.start_pos + frac_path_pos_1 * uv * l);
      clipper::Coord_orth pos_2(sdo.end_pos   + frac_path_pos_2 * uv * l);

      auto pp = std::make_pair(coord_orth_to_glm(pos_1), coord_orth_to_glm(pos_2));
      cylinder c(pp, 0.04, 0.04, segment_length, colour);
      c.add_flat_start_cap();
      c.add_flat_end_cap();
      import(vnc_vertex_vector_to_generic_vertex_vector(c.vertices), c.triangles);
   }
   setup(material);

}

void
Mesh::add_dashed_angle_markup(const glm::vec3 &pos_1, const glm::vec3 &pos_2, const glm::vec3 &pos_3,
                              const glm::vec4 &col, const Material &material) {

   // currently not dashed and no end caps.

   auto calculate_normal = [] (const glm::vec3 &d_1, const glm::vec3 &d_3) {
                              auto uv_1 = glm::normalize(d_1);
                              auto uv_3 = glm::normalize(d_3);
                              auto uv_n = glm::cross(uv_1, uv_3);
                              return uv_n;
                           };

   const unsigned int n_phi_steps = 20;
   const unsigned int n_theta_steps = 12;
   // std::vector<s_generic_vertex> local_vertices((n_theta_steps + 1) * n_phi_steps);
   std::vector<s_generic_vertex> local_vertices(n_theta_steps * (n_phi_steps + 1));

   // std::cout << "debug:: local_vertices resize to " << local_vertices.size() << std::endl;
   std::vector<g_triangle> local_triangles;
   float arc_radius = 0.7;
   float arc_radius_inner = 0.05;
   const float R = arc_radius;
   const float r = arc_radius_inner;
   const float pi = 3.1415926535;

   glm::vec3 d_1 = pos_1 - pos_2; // direction to pos_1
   glm::vec3 d_3 = pos_3 - pos_2; // direction to pos_3
   glm::vec3 z(0,0,1);
   glm::vec3 normal = calculate_normal(d_1, d_3);

   glm::vec3 centre_atom_position = pos_2;
   glm::mat4 ori_mat_1 = glm::orientation(normal, z);

   glm::mat3 rot_mat(ori_mat_1); // nope

   // try again for the rotation matrix
   //
   { // it's just magic (see arc_info_type constructor)
      glm::vec3 v1 = glm::normalize(d_3);
      glm::vec3 v2 = glm::normalize(d_1);
      glm::vec3 v3 = glm::cross(v1, v2);
      glm::vec3 v4 = glm::cross(v3, v1);
      rot_mat = glm::mat3(v1, v4, v3);
   }

   float theta_n = glm::angle(glm::normalize(d_1), glm::normalize(d_3));
   for (unsigned int ip=0; ip<=n_phi_steps; ip++) {
      float phi = theta_n * static_cast<float>(ip)/static_cast<float>(n_phi_steps);
      for (unsigned int it=0; it<n_theta_steps; it++) {
         float theta = 2.0 * pi * static_cast<float>(it)/static_cast<float>(n_theta_steps);
         s_generic_vertex v;
         float x = (R + r * cosf(theta)) * cosf(phi);
         float y = (R + r * cosf(theta)) * sinf(phi);
         float z = r * sinf(theta);
         glm::vec3 basic_pos(x,y,z);
         float normal_x = cosf(theta) * cosf(phi);
         float normal_y = cosf(theta) * sinf(phi);
         float normal_z = sinf(theta);
         glm::vec3 basic_normal(normal_x, normal_y, normal_z);
         v.pos = rot_mat * basic_pos + centre_atom_position;
         v.normal = rot_mat * basic_normal;
         v.color = col;
         unsigned int vertex_idx = ip * n_theta_steps + it;
         local_vertices[vertex_idx] = v;
      }
   }

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
         local_triangles.push_back(t1);
         local_triangles.push_back(t2);
      }
   }
   import(local_vertices, local_triangles);
   setup(material);


}


#include <fstream>

bool
Mesh::export_as_obj(const std::string &file_name) const {
   return export_as_obj_internal(file_name);
}


bool
Mesh::export_as_obj_internal(const std::string &file_name) const {

   bool status = true;

   std::cout << "debug:: export_as_obj_internal: n vertices:  " <<  vertices.size() << std::endl;
   std::cout << "debug:: export_as_obj_internal: n triangles: " << triangles.size() << std::endl;

   std::ofstream f(file_name.c_str());
   if (f) {
      f << "# " << name << "\n";
      f << "# " << "\n";
      f << "" << "\n";
      f << "g exported_obj\n";
      for (unsigned int i=0; i<vertices.size(); i++) {
         const s_generic_vertex &vert = vertices[i];
         f << "v " << vert.pos.x << " " << vert.pos.y << " " << vert.pos.z;
         f << " " << vert.color.r << " " << vert.color.g << " " << vert.color.b;
         f << "\n";
      }
      for (unsigned int i=0; i<vertices.size(); i++) {
         const s_generic_vertex &vert = vertices[i];
         f << "vn " << -vert.normal.x << " " << -vert.normal.y << " " << -vert.normal.z << "\n";
      }
      for (unsigned int i=0; i<triangles.size(); i++) {
         const g_triangle &tri = triangles[i];
         f << "f "
           << tri.point_id[0]+1 << "//" << tri.point_id[0]+1 << " "
           << tri.point_id[1]+1 << "//" << tri.point_id[1]+1 << " "
           << tri.point_id[2]+1 << "//" << tri.point_id[2]+1 << "\n";
      }
   } else {
      status = false;
   }
   return status;
}

bool
Mesh::export_as_obj(std::ofstream &f, unsigned int vertex_index_offset) const {

   bool status = true;
   if (f) {
      for (unsigned int i=0; i<vertices.size(); i++) {
         const s_generic_vertex &vert = vertices[i];
         f << "v " << vert.pos.x << " " << vert.pos.y << " " << vert.pos.z;
         // f << " " << vert.color.r << " " << vert.color.g << " " << vert.color.b;
         f << "\n";
      }
      for (unsigned int i=0; i<vertices.size(); i++) {
         const s_generic_vertex &vert = vertices[i];
         f << "vn " << vert.normal.x << " " << vert.normal.y << " " << vert.normal.z << "\n";
      }
      for (unsigned int i=0; i<triangles.size(); i++) {
         const g_triangle &tri = triangles[i];
         f << "f "
           << tri.point_id[0]+1+vertex_index_offset << "//" << tri.point_id[0]+1+vertex_index_offset << " "
           << tri.point_id[1]+1+vertex_index_offset << "//" << tri.point_id[1]+1+vertex_index_offset << " "
           << tri.point_id[2]+1+vertex_index_offset << "//" << tri.point_id[2]+1+vertex_index_offset << "\n";
      }
   }
   return status;
}



// We import from assimp from Model and export from Mesh - at the moment
//
#ifdef USE_ASSIMP
#include <assimp/Exporter.hpp>
#endif // USE_ASSIMP

bool
Mesh::export_as_obj_via_assimp(const std::string &file_name) const {

   unsigned int status = false;

   std::cout << "exporting to " << file_name << std::endl;

#ifdef USE_ASSIMP

   if (! vertices.empty()) {
      // this will do a copy. Is that what I want?
      aiScene scene = generate_scene();
      // now export scener to file_name;
      Assimp::Exporter ae;

      if (false) { // 20 formats - "obj" is one of them
         size_t n_formats = ae.GetExportFormatCount();
         for (size_t i=0; i<n_formats; i++) {
            const aiExportFormatDesc *fd = ae.GetExportFormatDescription(i);
            std::cout << i << " " << fd->id << " " << fd->description << std::endl;
         }
      }

      std::string format_id = "obj";

      std::cout << "----- calling ae.Export() " << std::endl;
      std::cout << "----- calling ae.Export() scene: " << &scene << std::endl;
      std::cout << "----- calling ae.Export() mMaterials " << scene.mMaterials << std::endl;
      std::cout << "----- calling ae.Export() mMaterials[0] " << scene.mMaterials[0] << std::endl;
      // Oh, scene gets copies in Export()!
      aiReturn air = ae.Export(&scene, format_id.c_str(), file_name.c_str(), 0);
      std::cout << "export status " << air << std::endl;

   }
#endif
   return status;
}

#if USE_ASSIMP
// Use make_shared and a shared or unique? pointer as the return value
aiScene
Mesh::generate_scene() const {

   aiScene scene;
   scene.mRootNode = new aiNode();

   std::cout << "debug:: scene.mNumMaterials " << scene.mNumMaterials << std::endl;

   // Materials? Maybe we have to.
   // else when Exporter::Export() does the scene copy the scene.mMaterials gets
   // set to null and GetMaterialName() fails

   scene.mNumMaterials = 1;
   scene.mMaterials = new aiMaterial *[1];
   scene.mMaterials[0] = new aiMaterial();

   aiMesh mesh;
   scene.mMeshes = new aiMesh *[1];
   scene.mNumMeshes = 1;
   scene.mMeshes[0] = new aiMesh();
   std::cout << "new aimesh mMaterialIndex " << scene.mMeshes[0]->mMaterialIndex << std::endl;
   scene.mMeshes[0]->mMaterialIndex = 0; // Hmm? Dangerous?
   scene.mRootNode->mMeshes = new unsigned int[1];
   scene.mRootNode->mMeshes[0] = 0;
   scene.mRootNode->mNumMeshes = 1;
   aiMesh *pMesh = scene.mMeshes[0];

   //  --- vertices ---

   pMesh->mVertices = new aiVector3D[vertices.size()];
   pMesh->mNumVertices = vertices.size();
   pMesh->mTextureCoords[0] = new aiVector3D[vertices.size()];       // needed?
   pMesh->mNumUVComponents[0] = vertices.size();                     // needed?
   for (unsigned int i=0; i<vertices.size(); i++) {
      const s_generic_vertex &vert = vertices[i];
      aiVector3D aiv(vert.pos.x, vert.pos.y, vert.pos.z);
      // std::cout << "vertex " << i << " " << glm::to_string(vert.pos)  << std::endl;
      pMesh->mVertices[i] = aiv;
      pMesh->mTextureCoords[0][i] = aiVector3D(0,0,0); // for now
   }

   //  --- normals ---


   //  --- triangles ---

   pMesh->mFaces = new aiFace[triangles.size()];
   pMesh->mNumFaces = triangles.size();
   for (unsigned int i=0; i<triangles.size(); i++) {
      aiFace &face = pMesh->mFaces[i];
      face.mIndices = new unsigned int[3];
      face.mIndices[0] = triangles[i].point_id[0];
      face.mIndices[1] = triangles[i].point_id[1];
      face.mIndices[2] = triangles[i].point_id[2];
   }

   return scene;
}
#endif


void
Mesh::debug_to_file() const {

   std::string n = "debug-mesh-" + name;
   std::ofstream f(n);
   std::cout << "# number of vertices: " << vertices.size() << std::endl;
   std::cout << "# number of triangles: " << triangles.size() << std::endl;

   for (unsigned int i=0; i<vertices.size(); i++) {
      const auto &vertex = vertices[i];
      f << "vertex " << i << " : v: " << glm::to_string(vertex.pos) << " n: "
        << glm::to_string(vertex.normal) << " c: " << glm::to_string(vertex.color) << "\n";
   }
   f.close();

}

#include <glm/gtx/norm.hpp>

void
Mesh::sort_map_triangles(const glm::vec3 &eye_position) {

   bool do_the_sort = false;
   if (glm::length2(eye_position - previous_eye_position) > 0.0001) do_the_sort = true;

   if (! do_the_sort) return;

   if (false) { // debug
      for (unsigned int i=0; i<map_triangle_centres.size(); i++) {
         auto delta(map_triangle_centres[i].second.mid_point - eye_position);
         std::cout << "triangle " << i << " " << glm::to_string(map_triangle_centres[i].second.mid_point) << " "
                   << sqrt(glm::length2(delta)) << std::endl;
      }
   }

   if (false)
      std::cout << "sorting map triangles " << map_triangle_centres.size() << " triangles with eye position"
                << glm::to_string(eye_position) << std::endl;
   auto tp_0 = std::chrono::high_resolution_clock::now();

   for (unsigned int i=0; i<map_triangle_centres.size(); i++) {
      auto delta(map_triangle_centres[i].second.mid_point - eye_position);
      float dd = glm::length2(delta);
      map_triangle_centres[i].second.back_front_projection_distance = dd;
   }

   // this sign needs checking (I did, I think that it's right now).
   auto map_triangle_sorter = [] (const std::pair<int, map_triangle_t> &t1,
                                  const std::pair<int, map_triangle_t> &t2) {
                                 return (t1.second.back_front_projection_distance < t2.second.back_front_projection_distance);
                              };

   std::sort(map_triangle_centres.begin(), map_triangle_centres.end(), map_triangle_sorter);

   unsigned int n_triangle_centres = map_triangle_centres.size();

   int *indices_for_triangles = new int[3 * n_triangle_centres]; // d (hot-path mallocing! Woo!)
   for (unsigned int i=0; i<map_triangle_centres.size(); i++) {
      indices_for_triangles[3*i  ] = map_triangle_centres[i].second.point_id[0];
      indices_for_triangles[3*i+1] = map_triangle_centres[i].second.point_id[1];
      indices_for_triangles[3*i+2] = map_triangle_centres[i].second.point_id[2];
   }

   // if (xmap_is_diff_map)
   // return;

   GLenum err = glGetError();

   glBindVertexArray(vao);
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer_id);
   err = glGetError(); if (err) std::cout << "GL error: sorting triangles: " << err << std::endl;

   unsigned int n_bytes = 3 * n_triangle_centres * sizeof(unsigned int);
   glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, n_bytes, &indices_for_triangles[0]);
   err = glGetError(); if (err) std::cout << "GL error: sorting triangles: " << err << std::endl;
   glBindVertexArray(0);

   delete [] indices_for_triangles;

   // for next time
   previous_eye_position = eye_position;

   auto tp_1 = std::chrono::high_resolution_clock::now();
   auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
   // std::cout << "INFO:: sorting triangles " << d10 << " milliseconds" << std::endl;

}


// export suitable for Blender - just the vertices.
std::vector<glm::vec3>
Mesh::just_vertices() const {

   std::vector<glm::vec3> v(vertices.size());
   for (unsigned int i=0; i<vertices.size(); i++) {
      v[i] = vertices[i].pos;
   }
   return v;
}




void
Mesh::invert_normals() { // flip normals

   for (auto &vertex : vertices)
      vertex.normal = -vertex.normal;

}

void
Mesh::calculate_normals() {

   std::map<unsigned int, std::vector<glm::vec3> > normal_map;
   for (unsigned int i=0; i<triangles.size(); i++) {
      const auto &triangle = triangles[i];
      const glm::vec3 &v0 = vertices[triangle.point_id[0]].pos;
      const glm::vec3 &v1 = vertices[triangle.point_id[1]].pos;
      const glm::vec3 &v2 = vertices[triangle.point_id[2]].pos;
      glm::vec3 d1 = v1 - v0;
      glm::vec3 d2 = v2 - v0;
      glm::vec3 c = glm::cross(d1, d2);
      glm::vec3 n = -glm::normalize(c);
      bool b = glm::any(glm::isnan(n));
      if (!b) {
         normal_map[triangle.point_id[0]].push_back(n);
         normal_map[triangle.point_id[1]].push_back(n);
         normal_map[triangle.point_id[2]].push_back(n);
      }
   }
   std::map<unsigned int, std::vector<glm::vec3> >::const_iterator it;
   for (it=normal_map.begin(); it!=normal_map.end(); ++it) {
      const unsigned int &idx = it->first;
      const std::vector<glm::vec3> &vecs = it->second;
      glm::vec3 sum(0,0,0);
      for (unsigned int j=0; j<vecs.size(); j++)
         sum += vecs[j];
      float fact = 1.0f/static_cast<float>(vecs.size());
      glm::vec3 av_norm = sum *= fact;
      std::cout << glm::to_string(av_norm) << "\n";
      vertices[idx].normal = av_norm;
   }
   setup_buffers();

}


void
Mesh::remove_last_subobject(unsigned int n_vertices, unsigned int n_triangles) {

   if (triangles.size() >= n_triangles)
      triangles.erase(triangles.end() - n_triangles, triangles.end());
}
