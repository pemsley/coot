/*
 * src/model-molecule-meshes.cc
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

// #include <glm/ext.hpp>
#define GLM_ENABLE_EXPERIMENTAL // # for to_string()
#include <glm/gtx/string_cast.hpp>  // to_string()
#include "model-molecule-meshes.hh"
#include "api/make-instanced-graphical-bonds.hh" // make_instanced_graphical_bonds_spherical_atoms() etc
#include "generic-vertex.hh"

#include "graphics-info.h" // for get_residue - used for Rotamer markup.

bool
model_molecule_meshes_t::empty() const {

   bool state = false;
   if (instanced_meshes.empty())
      if (simple_mesh.empty())
         state = true;

   return state;
}

// transfer this to the meshes
void
model_molecule_meshes_t::set_debug_mode(bool state) {

   // std::vector<Mesh> instanced_meshes;
   // Mesh simple_mesh;

   simple_mesh.debug_mode = state; // make this a function?
   for (unsigned int i=0; i<instanced_meshes.size(); i++)
      instanced_meshes[i].debug_mode = state;

}


void
model_molecule_meshes_t::draw_simple_bond_lines(Shader *shader_p,
                                                const glm::mat4 &mvp,
                                                const glm::vec4 &background_colour,
                                                float line_width,
                                                bool do_depth_fog) {

   simple_mesh.draw_simple_bond_lines(shader_p, mvp, background_colour, line_width, do_depth_fog);
}


// wrap draw_simple_for_ssao() and draw_instances_for_ssao()
void
model_molecule_meshes_t::draw_for_ssao(Shader *shader_for_meshes_p,
                                       Shader *shader_for_instanced_meshes_p,
                                       const glm::mat4 &model,
                                       const glm::mat4 &view,
                                       const glm::mat4 &projection) { // draw into the gbuffer framebuffer.

   simple_mesh.draw_for_ssao(shader_for_meshes_p, model, view, projection);
   for (unsigned int i=0; i<instanced_meshes.size(); i++) {
      instanced_meshes[i].draw_instances_for_ssao(shader_for_instanced_meshes_p, model, view, projection);
   }
}


void
model_molecule_meshes_t::convert_and_fill_meshes(const coot::instanced_mesh_t &im) {

   for (unsigned int i_g=0; i_g<im.geom.size(); i_g++) {
      const coot::instanced_geometry_t &ig = im.geom[i_g];
      // mesh
      std::vector<s_generic_vertex> src_vertices(ig.vertices.size());
      glm::vec4 col(0.5f, 0.5f, 0.5f, 1.0f);
      for (unsigned int iv=0; iv<ig.vertices.size(); iv++)
         src_vertices[iv] = s_generic_vertex(ig.vertices[iv].pos,
                                             ig.vertices[iv].normal,
                                             col);
      Mesh m(src_vertices, ig.triangles);
      std::string mesh_name("make_graphical_bonds instancing-mesh-" + std::to_string(i_g) + " " + ig.name);
      m.set_name(mesh_name);
      m.setup_buffers();

      if (! ig.instancing_data_A.empty()) {
         // instancing data
         glm::mat4 unit_matrix(1.0f);
         std::vector<glm::mat4> matrices(ig.instancing_data_A.size(), unit_matrix);
         std::vector<glm::vec4> colours(ig.instancing_data_A.size());
         for (unsigned int i_A=0; i_A<ig.instancing_data_A.size(); i_A++) {
            const coot::instancing_data_type_A_t &Atd = ig.instancing_data_A[i_A];
            // 20230828-PE Atom sizes (sphere radius) can vary. Waters and metals are bigger.
            // size for a sphere should be a vector of size 3 of the same number.
	    float sphere_radius = Atd.size[0];
            colours[i_A] = Atd.colour;
	    glm::mat4 mm = glm::scale(unit_matrix, Atd.size);
            matrices[i_A] = glm::translate(mm, Atd.position/sphere_radius);
         }
         Material material;
         m.setup_instancing_buffer_data(material, matrices, colours);
         instanced_meshes.push_back(m);
      }

      if (! ig.instancing_data_B.empty()) {
         // instancing data
         glm::mat4 unit_matrix(1.0f);
         std::vector<glm::mat4> matrices(ig.instancing_data_B.size(), unit_matrix);
         std::vector<glm::vec4> colours(ig.instancing_data_B.size());
         for (unsigned int i_B=0; i_B<ig.instancing_data_B.size(); i_B++) {
            const coot::instancing_data_type_B_t &Btd = ig.instancing_data_B[i_B];
            colours[i_B] = Btd.colour;
	    glm::mat4 sc = glm::scale(unit_matrix, Btd.size);
            glm::mat4 ori = Btd.orientation; // make ref
            glm::mat4 t = glm::translate(unit_matrix, Btd.position);
            glm::mat4 mm = t * ori * sc;
            matrices[i_B] = mm;
         }
         Material material;
         m.setup_instancing_buffer_data(material, matrices, colours);
         instanced_meshes.push_back(m);
      }

   }

   // the "markup" mesh is "simple" i.e not instanced

   if (! im.markup.vertices.empty()) {
      Material material;
      material.turn_specularity_on(true);
      material.shininess = 256;
      material.specular_strength = 0.8;
      std::vector<s_generic_vertex> src_vertices(im.markup.vertices.size());
      const auto &triangles = im.markup.triangles;
      for (unsigned int i=0; i<im.markup.vertices.size(); i++) {
         const auto &v = im.markup.vertices[i];
         src_vertices[i] = s_generic_vertex(v.pos, v.normal, v.color);
      }
      simple_mesh.name = "markup";
      simple_mesh.set_material(material);
      simple_mesh.import(src_vertices, triangles);
      simple_mesh.setup_buffers();
   }
}


void
model_molecule_meshes_t::make_graphical_bonds(int imol, const graphical_bonds_container &bonds_box,
                                              float atom_radius, float bond_radius,
                                              bool show_atoms_as_aniso_flag,
                                              float aniso_probability,
                                              bool show_aniso_atoms_as_ortep_flag,
                                              int num_subdivisions, int n_slices, int n_stacks,
                                              const std::vector<glm::vec4> &colour_table) {

   //    std::cout << "DEBUG:: in model_molecule_meshes_t::make_graphical_bonds()) with show_aniso_atoms_as_ortep_flag "
   // << show_aniso_atoms_as_ortep_flag << std::endl;

   // who calls this function?

   // first clear the buffers of what we ready (might) have
   instanced_meshes.clear();
   simple_mesh.clear();
   im.clear();

   // api functions
   coot::api_bond_colour_t dummy_bonds_box_type(coot::api_bond_colour_t::NORMAL_BONDS);

   // std::cout << "calling make_instanced_graphical_bonds_spherical_atoms() with atom_radius " << atom_radius << std::endl;
   // Atom radii are limited to 2.0

   // these are in api/coot-molecule-bonds-instanced
   make_instanced_graphical_bonds_spherical_atoms(im, bonds_box, dummy_bonds_box_type, atom_radius, bond_radius,
                                                  show_atoms_as_aniso_flag, aniso_probability,
                                                  show_aniso_atoms_as_ortep_flag,
                                                  num_subdivisions, colour_table);
   make_instanced_graphical_bonds_bonds(im, bonds_box, bond_radius, n_slices, n_stacks, colour_table);
   make_graphical_bonds_cis_peptides(im.markup, bonds_box);
   add_rotamer_dodecs(imol, bonds_box);
   add_ramachandran_spheres(imol, bonds_box);


   // ===================================== now convert instancing.hh meshes to src style "Mesh"es =======================

   convert_and_fill_meshes(im); // fill class data meshes

}


void
model_molecule_meshes_t::make_symmetry_bonds(int imol,
                                             const std::vector<std::pair<graphical_bonds_container, std::pair<symm_trans_t, Cell_Translation> > > &symmetry_bonds_box,
                                             float atom_radius, float bond_radius,
                                             int num_subdivisions, int n_slices, int n_stacks,
                                             const std::vector<glm::vec4> &colour_table) {

   auto debug_instancing_meshes = [] (const coot::instanced_mesh_t &im) {
      std::cout << "im has " << im.geom.size() << " geoms" << std::endl;
      for (unsigned int ig=0; ig<im.geom.size(); ig++) {
         unsigned int n_A = im.geom[ig].instancing_data_A.size();
         unsigned int n_B = im.geom[ig].instancing_data_B.size();
         std::cout << "   geom " << ig << " has " << im.geom[ig].vertices.size() << " vertices and "
                   << im.geom[ig].triangles.size() << " triangles "
                   << n_A << " type A instances and " << n_B << " type B instances" << std::endl;
      }
   };

   // each symmetry_bonds item the vector needs to be tested with:
   // if (gbc.symmetry_has_been_created == 1)
   // before looping through the colours

   instanced_meshes.clear();
   im.clear();
   for (unsigned int i=0; i<symmetry_bonds_box.size(); i++) {
      const std::pair<graphical_bonds_container, std::pair<symm_trans_t, Cell_Translation> > &s_box = symmetry_bonds_box[i];
      const graphical_bonds_container &gbc = s_box.first;
      if (gbc.symmetry_has_been_created == 1) {
         make_instanced_graphical_symmetry_bonds_bonds(im, s_box, bond_radius, n_slices, n_stacks, colour_table);
      } else {
         // std::cout << "DEBUG:: skip this box" << std::endl;
      }
   }

   // debug_instancing_meshes(im);
   convert_and_fill_meshes(im); // fill class data meshes
}


mmdb::Residue *
model_molecule_meshes_t::get_residue(int imol, const coot::residue_spec_t &spec) const {

   mmdb::Residue *residue_p = graphics_info_t::get_residue(imol, spec);
   return residue_p;
}

std::pair<bool, coot::Cartesian>
model_molecule_meshes_t::get_HA_unit_vector(mmdb::Residue *r) const {
   bool status = false;
   coot::Cartesian dir;
   mmdb::Atom *CA = r->GetAtom(" CA ");
   mmdb::Atom *C  = r->GetAtom(" C  ");
   mmdb::Atom *N  = r->GetAtom(" N  ");
   mmdb::Atom *CB = r->GetAtom(" CB ");

   if (CA && C && N && CB) {
      coot::Cartesian ca_pos(CA->x, CA->y, CA->z);
      coot::Cartesian  c_pos( C->x,  C->y,  C->z);
      coot::Cartesian  n_pos( N->x,  N->y,  N->z);
      coot::Cartesian cb_pos(CB->x, CB->y, CB->z);
      coot::Cartesian dir_1 = ca_pos - c_pos;
      coot::Cartesian dir_2 = ca_pos - n_pos;
      coot::Cartesian dir_3 = ca_pos - cb_pos;
      coot::Cartesian r = dir_1 + dir_2 + dir_3;
      dir = r.unit();
      status = true;
   } else {
      if (CA && C && N) {
         coot::Cartesian ca_pos(CA->x, CA->y, CA->z);
         coot::Cartesian  c_pos( C->x,  C->y,  C->z);
         coot::Cartesian  n_pos( N->x,  N->y,  N->z);
         coot::Cartesian dir_1 = ca_pos - c_pos;
         coot::Cartesian dir_2 = ca_pos - n_pos;
         coot::Cartesian r = dir_1 + dir_2;
         dir = r.unit();
         status = true;
      }
   }
   return std::make_pair(status, dir);
}

// this function should live somewhere more basic, I think. perhaps coot-utils (shapes.hh?)
std::pair<std::vector<coot::api::vn_vertex>, std::vector<g_triangle> >
model_molecule_meshes_t::get_dodec_vertices_and_triangles() const {

   auto clipper_to_glm = [] (const clipper::Coord_orth &c) {
                              return glm::vec3(c.x(), c.y(), c.z());
                           };

   dodec d;
   std::vector<clipper::Coord_orth> coords = d.coords();
   std::vector<glm::vec3> dodec_postions(coords.size());
   for (unsigned int i=0; i<coords.size(); i++)
      dodec_postions[i] = clipper_to_glm(coords[i]);

   std::vector<coot::api::vn_vertex> dodec_vertices;
   std::vector<g_triangle> dodec_triangles;
   dodec_triangles.reserve(36);

   for (unsigned int iface=0; iface<12; iface++) {

      std::vector<coot::api::vn_vertex> face_verts;
      std::vector<g_triangle> face_triangles;
      face_triangles.reserve(3);

      std::vector<unsigned int> indices_for_face = d.face(iface);
      glm::vec3 ns(0,0,0);
      for (unsigned int j=0; j<5; j++)
         ns += dodec_postions[indices_for_face[j]];
      glm::vec3 normal = glm::normalize(ns);

      for (unsigned int j=0; j<5; j++) {
         glm::vec3 &pos = dodec_postions[indices_for_face[j]];
         coot::api::vn_vertex v(0.5f * pos, normal);
         face_verts.push_back(v);
      }

      face_triangles.push_back(g_triangle(0,1,2));
      face_triangles.push_back(g_triangle(0,2,3));
      face_triangles.push_back(g_triangle(0,3,4));

      unsigned int idx_base = dodec_vertices.size();
      unsigned int idx_tri_base = dodec_triangles.size();
      dodec_vertices.insert(dodec_vertices.end(), face_verts.begin(), face_verts.end());
      dodec_triangles.insert(dodec_triangles.end(), face_triangles.begin(), face_triangles.end());
      for (unsigned int jj=idx_tri_base; jj<dodec_triangles.size(); jj++)
         dodec_triangles[jj].rebase(idx_base);
   }

   return std::make_pair(dodec_vertices, dodec_triangles);
}

void
model_molecule_meshes_t::add_rotamer_dodecs(int imol, const graphical_bonds_container &bonds_box) {

   auto cartesian_to_glm = [] (const coot::Cartesian &c) {
                              return glm::vec3(c.x(), c.y(), c.z());
                           };
   auto clipper_to_glm = [] (const clipper::Coord_orth &c) {
                              return glm::vec3(c.x(), c.y(), c.z());
                           };
   auto colour_holder_to_glm = [] (const coot::colour_holder &ch) {
                                  return glm::vec4(ch.red, ch.green, ch.blue, 1.0f);
                               };
   // std::cout << "DEBUG:: in make_graphical_bonds(): n_rotamer_markups: " << bonds_box.n_rotamer_markups << std::endl;
   if (bonds_box.n_rotamer_markups > 0) {

      std::pair<std::vector<coot::api::vn_vertex>, std::vector<g_triangle> > vtd = get_dodec_vertices_and_triangles();
      const std::vector<coot::api::vn_vertex> &dodec_vertices = vtd.first;
      const std::vector<g_triangle> &dodec_triangles = vtd.second;

      coot::instanced_geometry_t ig(dodec_vertices, dodec_triangles);
      ig.name = "Rotamer dodecs";

      double rama_ball_pos_offset_scale = 1.5;
      float size = 1.0;
      glm::vec3 size_3(size, size, size);
      for (int i=0; i<bonds_box.n_rotamer_markups; i++) {
         const rotamer_markup_container_t &rm = bonds_box.rotamer_markups[i];
         if (rm.rpi.state == coot::rotamer_probability_info_t::OK) { // should be coot::rotamer_probability_info_t::OK
            const coot::residue_spec_t &residue_spec = rm.spec;
            // std::cout << "Rotamer-markup " << i << " " << residue_spec << std::endl;
            mmdb::Residue *residue_p = get_residue(imol, residue_spec);
            coot::Cartesian offset(0,0,rama_ball_pos_offset_scale);
            if (residue_p) {
               std::pair<bool, coot::Cartesian> hav = get_HA_unit_vector(residue_p);
               if (hav.first) offset = hav.second * rama_ball_pos_offset_scale;
            } else {
               std::cout << "in add_rotamer_dodecs() sadge null residue_p " << imol << " " << residue_spec << std::endl;
            }
            glm::vec3 atom_pos = clipper_to_glm(rm.pos) + cartesian_to_glm(offset);
            auto rm_col = rm.col;
            rm_col.scale_intensity(0.75); // was 0.6 in Mesh-from-graphical-bonds.cc
            auto this_dodec_colour = colour_holder_to_glm(rm_col);
            coot::instancing_data_type_A_t idA(atom_pos, this_dodec_colour, size_3);
            ig.instancing_data_A.push_back(idA);
         }
      }
      im.add(ig);
   }
}

#include "coot-utils/oct.hh"                                        

void
model_molecule_meshes_t::add_ramachandran_spheres(int imol, const graphical_bonds_container &gbc) {

   // 20230828-PE I am "doing this" yet again because this version uses instancing.

   auto colour_holder_to_glm = [] (const coot::colour_holder &ch) {
                                  return glm::vec4(ch.red, ch.green, ch.blue, 1.0f);
                               };

   // This should match the HUD rotamer bar colour.
   //
   auto prob_raw_to_colour_rotation = [] (float prob) {
                                         if (prob > 0.5) prob = 0.5; // 0.4 and 2.5 f(for q) might be better (not tested)
                                         // good probabilities have q = 0
                                         // bad probabilities have q 0.66
                                         double q = (1.0 - 2.0 * prob);
                                         q = pow(q, 20.0);
                                         return q;
                                      };

   if (gbc.n_ramachandran_goodness_spots > 0) {
      unsigned int num_subdivisions = 2;
      std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > octaball = tessellate_octasphere(num_subdivisions);
      const auto &octaball_vertices_raw  = octaball.first;
      const auto &octaball_triangles_raw = octaball.second;
      std::vector<coot::api::vn_vertex> octaball_vertices(octaball_vertices_raw.size()); // cooked
      for (unsigned int i=0; i<octaball_vertices_raw.size(); i++) {
         const auto &v = octaball_vertices_raw[i];
         octaball_vertices[i] = coot::api::vn_vertex(v,v);
      }
      coot::instanced_geometry_t ig(octaball_vertices, octaball_triangles_raw);
      float size = 0.5f;
      glm::vec3 size_3(size, size, size);
      for (int i=0; i<gbc.n_ramachandran_goodness_spots; i++) {
         const coot::Cartesian &position = gbc.ramachandran_goodness_spots_ptr[i].first;
         const float &prob_raw = gbc.ramachandran_goodness_spots_ptr[i].second;
         double q = prob_raw_to_colour_rotation(prob_raw);
         coot::colour_holder col = coot::colour_holder(q, 0.0, 1.0, false, std::string(""));
         col.scale_intensity(0.6); // calm down the otherwise super-bright Rama ball colours
         glm::vec4 col_glm = colour_holder_to_glm(col);
         glm::vec3 ball_position = cartesian_to_glm(position);
         coot::instancing_data_type_A_t idA(ball_position, col_glm, size_3);
         ig.instancing_data_A.push_back(idA);
      }
      im.add(ig);
   }
}


// the simple-lines option for the main molecule
void
model_molecule_meshes_t::make_bond_lines(const graphical_bonds_container &bonds_box, const std::vector<glm::vec4> &colour_table) {

   std::cout << "model_molecule_meshes_t::make_bond_lines() missing function" << std::endl;
}


// wrapper for both the instanced and simple meshes
void
model_molecule_meshes_t::draw(Shader *shader_mesh,
                              Shader *shader_for_instanced_mesh,
                              stereo_eye_t eye,
                              const glm::mat4 &mvp,
                              const glm::mat4 &view_rotation_matrix,
                              const std::map<unsigned int, lights_info_t> &lights,
                              const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                              float opacity,
                              const glm::vec4 &background_colour,
                              bool gl_lines_mode, // i.e. as chickenwire
                              bool do_depth_fog,
                              bool show_just_shadows) {

   bool transferred_colour_is_instanced = true;
   glm::vec3 dummy_rc(0,0,0);
   draw_instances(shader_for_instanced_mesh, mvp, view_rotation_matrix, lights, eye_position, background_colour, do_depth_fog, transferred_colour_is_instanced);
   draw_simple(shader_mesh, mvp, view_rotation_matrix, lights, eye_position, dummy_rc, opacity, background_colour, gl_lines_mode, do_depth_fog, show_just_shadows);

}

void
model_molecule_meshes_t::draw_instances(Shader *shader_for_instanced_meshes_p,
                                        const glm::mat4 &mvp,
                                        const glm::mat4 &view_rotation_matrix,
                                        const std::map<unsigned int, lights_info_t> &lights,
                                        const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                                        const glm::vec4 &background_colour,
                                        bool do_depth_fog,
                                        bool transferred_colour_is_instanced,
                                        bool do_pulse,
                                        bool do_rotate_z,
                                        float pulsing_amplitude,
                                        float pulsing_frequency,
                                        float pulsing_phase_distribution,
                                        float z_rotation_angle) {

   if (false)
      std::cout << "model_molecule_meshes_t::draw_instances(): n-instanced meshes: " << instanced_meshes.size() << std::endl;
   for (unsigned int i=0; i<instanced_meshes.size(); i++) {
      auto &mesh = instanced_meshes[i];
      if (false)
         std::cout << "   calling mesh.draw_instanced() \"" << mesh.name << "\" with shader "
                   << "\"" << shader_for_instanced_meshes_p->name << "\"" << std::endl;
      int pass_type = graphics_info_t::PASS_TYPE_WITH_SHADOWS; // is this right?
      mesh.draw_instanced(pass_type, shader_for_instanced_meshes_p, mvp, view_rotation_matrix, lights, eye_position, background_colour,
                          do_depth_fog, transferred_colour_is_instanced,
                          do_pulse, do_rotate_z, pulsing_amplitude, pulsing_frequency,
                          pulsing_phase_distribution, z_rotation_angle);
   }

}

// the (previous) Mesh draw function (simple here means simple mesh, not lines)
void
model_molecule_meshes_t::draw_simple(Shader *shader,
                                     const glm::mat4 &mvp,
                                     const glm::mat4 &view_rotation_matrix,
                                     const std::map<unsigned int, lights_info_t> &lights,
                                     const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                                     const glm::vec3 &rotation_centre,
                                     float opacity,
                                     const glm::vec4 &background_colour,
                                     bool gl_lines_mode, // i.e. as chickenwire
                                     bool do_depth_fog,
                                     bool show_just_shadows) {

   // std::cout << "model_molecule_meshes_t::draw_simple()" << std::endl;
   simple_mesh.draw(shader, mvp, view_rotation_matrix, lights, eye_position, rotation_centre,
                    opacity, background_colour, gl_lines_mode, do_depth_fog, show_just_shadows);
}

// instanced models
void
model_molecule_meshes_t::draw_molecule_with_shadows(Shader *shader,
                                                    const glm::mat4 &mvp,
                                                    const glm::mat4 &model_rotation_matrix,
                                                    const std::map<unsigned int, lights_info_t> &lights,
                                                    const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                                                    float opacity,
                                                    const glm::vec4 &background_colour,
                                                    bool do_depth_fog,
                                                    const glm::mat4 &light_view_mvp,
                                                    unsigned int shadow_depthMap,
                                                    float shadow_strength,
                                                    unsigned int shadow_softness, // 1, 2 or 3.
                                                    bool show_just_shadows) {

   std::vector<Mesh>::iterator it;
   for (it=instanced_meshes.begin(); it!=instanced_meshes.end(); ++it) {
      auto &mesh(*it);
      // std::cout << "in draw_molecule_as_meshes_with_shadows() drawing instanced_mesh \"" << mesh.name << "\"" << std::endl;
      mesh.draw_with_shadows(shader, mvp, model_rotation_matrix, lights,
                             eye_position, opacity, background_colour, do_depth_fog, light_view_mvp,
                             shadow_depthMap, shadow_strength, shadow_softness, show_just_shadows);
   }

}
