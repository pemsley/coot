

// this is the non-instanced version

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>  // to_string()

#include "Mesh.hh"
#include "bond-colour-mode.hh"
#include "coot-utils/oct.hh"
#include "coot-utils/cylinder.hh"

// a wrapper for the following functions
void
Mesh::make_graphical_bonds(const graphical_bonds_container &gbc,
                           int bonds_box_type,
                           unsigned int representation_type, // BALL_AND_STICK or BALLS_NOT_BONDS
                           int udd_handle_bonded_type,
                           bool draw_cis_peptides,
                           float atom_radius,
                           float bond_radius,
                           unsigned int num_subdivisions,
                           unsigned int n_slices,
                           unsigned int n_stacks,
                           const std::vector<glm::vec4> &colour_table) {

   // need to add:
   // cis peptides,
   // missing residue loops
   // and rama balls if intermediate atoms.

   if (colour_table.empty())
      std::cout << "ERROR:: :::::::::::::::::::::: empty colour_table() in Mesh::make_graphical_bonds() " << std::endl;

   clear();

   unsigned int n_bonds = 0;
   for (int icol_bond=0; icol_bond<gbc.num_colours; icol_bond++) {
      graphical_bonds_lines_list<graphics_line_t> &ll = gbc.bonds_[icol_bond];
      n_bonds += ll.num_lines;
   }
   unsigned int allocation_for_vertices  = 68 * n_bonds;
   unsigned int allocation_for_triangles = 80 * n_bonds;

   vertices.reserve(allocation_for_vertices);
   triangles.reserve(allocation_for_triangles);

   make_graphical_bonds_spherical_atoms(gbc, bonds_box_type, udd_handle_bonded_type, atom_radius, bond_radius, num_subdivisions, colour_table);
   make_graphical_bonds_hemispherical_atoms(gbc, bonds_box_type, udd_handle_bonded_type, atom_radius, bond_radius, num_subdivisions, colour_table);
   if (representation_type == BALL_AND_STICK)
      make_graphical_bonds_bonds(gbc, bond_radius, n_slices, n_stacks, colour_table);

   if (draw_cis_peptides)
      make_graphical_bonds_cis_peptides(gbc);

   // We shouldn't make Rama balls here - they should be made in an InstancedMesh - rama balls have their own
   // (at the moment unused) layout and shader.
   //
   glm::vec3 screen_up_dir(0,1,0); // for now
   make_graphical_bonds_rama_balls(gbc, screen_up_dir);

   make_graphical_bonds_rotamer_dodecs(gbc, screen_up_dir);

   setup_buffers();
}

void
Mesh::make_graphical_bonds_spherical_atoms(const graphical_bonds_container &gbc,
                                           int bonds_box_type,
                                           int udd_handle_bonded_type,
                                           float atom_radius,
                                           float bond_radius,
                                           unsigned int num_subdivisions,
                                           const std::vector<glm::vec4> &colour_table) {

   // udd_handle_bonded_type can be NO_BOND, BONDED_WITH_STANDARD_ATOM_BOND, BONDED_WITH_BOND_TO_HYDROGEN
   // BONDED_WITH_HETATM_BOND.

   auto cartesian_to_glm = [] (const coot::Cartesian &co) {
                            return glm::vec3(co.x(), co.y(), co.z());
                         };
   GLenum err = glGetError();
   if (err) std::cout << "error make_graphical_bonds_spherical_atoms() --start-- error "
                      << err << std::endl;

   bool atoms_have_bigger_radius_than_bonds = false;
   if (atom_radius > bond_radius)
      atoms_have_bigger_radius_than_bonds = true;

   // ----------------------- setup the vertices and triangles ----------------------

   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > octasphere_geom =
      tessellate_octasphere(num_subdivisions);

   is_instanced = false;

   glm::mat4 unit(1.0);
   for (int icol=0; icol<gbc.n_consolidated_atom_centres; icol++) {
      glm::vec4 col = colour_table[icol];
      for (unsigned int i=0; i<gbc.consolidated_atom_centres[icol].num_points; i++) {
         const graphical_bonds_atom_info_t &at_info = gbc.consolidated_atom_centres[icol].points[i];
         bool do_it = atoms_have_bigger_radius_than_bonds;
         mmdb::Atom *at = at_info.atom_p;

         if (! do_it) {
            if (at) {
               int state = -1;
               at->GetUDData(udd_handle_bonded_type, state);
               if (state == graphical_bonds_container::NO_BOND) {
                  do_it = true;
               }
            }
         }

         if (do_it) {
            unsigned int idx_base = vertices.size();
            unsigned int idx_tri_base = triangles.size();
            float scale = 1.0;
            if (at_info.is_hydrogen_atom) {
               if (atoms_have_bigger_radius_than_bonds) {
                  scale *= 0.66;
               } else {
                  scale *= 0.5; // bonds go to half-width, so should atoms.
               }
            }
            glm::vec3 t = cartesian_to_glm(at_info.position);  // (at->x, at->y, at->z);
            float sar = scale * atom_radius * at_info.radius_scale;
            if (at_info.is_water) {
               if (atoms_have_bigger_radius_than_bonds) {
                  float f = 1.33; // with a radius_scale of 2.0 waters are too chonky
                  sar = scale * atom_radius * f;
               }
            }
            glm::vec3 sc(sar, sar, sar);
            glm::mat4 mm = glm::scale(unit, sc);
            mm = glm::translate(mm, t);

            std::vector<s_generic_vertex> local_vertices(octasphere_geom.first.size());

            for (unsigned int ii=0; ii<local_vertices.size(); ii++) {
               auto &vert = local_vertices[ii];
               glm::vec3 p = octasphere_geom.first[ii] * sc + t;
               vert = s_generic_vertex(p, octasphere_geom.first[ii], col);
            }
            vertices.insert(vertices.end(), local_vertices.begin(), local_vertices.end());
            triangles.insert(triangles.end(), octasphere_geom.second.begin(), octasphere_geom.second.end());
            for (unsigned int k=idx_tri_base; k<triangles.size(); k++)
               triangles[k].rebase(idx_base);
         }
      }
   }

   err = glGetError();
   if (err) std::cout << "GL ERROR:: make_graphical_bonds_spherical_atoms() error --end-- " << err << std::endl;


}

void
Mesh::make_graphical_bonds_hemispherical_atoms(const graphical_bonds_container &gbc,
                                               int bonds_box_type,
                                               int udd_handle_bonded_type,
                                               float atom_radius,
                                               float bond_radius,
                                               unsigned int num_subdivisions,
                                               const std::vector<glm::vec4> &colour_table) {


   // BONDED_WITH_HETATM_BOND.

   // do these need to be passed to get_glm_colour_for_bonds()?

   bool atoms_have_bigger_radius_than_bonds = false;
   if (atom_radius > bond_radius)
      atoms_have_bigger_radius_than_bonds = true;

   // like above, different z axis because we want the hemisphere to extend outside the cylinder - and we don't need to
   // scale to bond length
   auto get_octahemi_matrix = [] (const glm::vec3 &pos_1, const glm::vec3 &pos_2, float radius) {
                                 glm::vec3 delta = pos_2 - pos_1;
                                 glm::mat4 u(1.0f);
                                 // orient
                                 glm::vec3 normalized_bond_orientation(glm::normalize(delta));
                                 glm::mat4 ori = glm::orientation(normalized_bond_orientation, glm::vec3(0.0, 0.0, -1.0));

                                 // 20220204-PE no translation because I want to use this matrix for the vertex normals
                                 // translate
                                 // glm::mat4 t = glm::translate(u, pos_1);
                                 // glm::mat4 m = t * ori * sc;
                                 return ori;
                              };

   auto cartesian_to_glm = [] (const coot::Cartesian &co) {
                            return glm::vec3(co.x(), co.y(), co.z());
                         };

   // ----------------------- setup the vertices and triangles ----------------------

   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > octasphere_geom =
      tessellate_hemisphere_patch(num_subdivisions);

   is_instanced = false;

   GLenum err = glGetError();
   if (err) std::cout << "error make_graphical_bonds_hemispherical_atoms() post-setup_buffers() error "
                      << err << std::endl;

   // ----------------------- setup the instances ----------------------

   // what is the other atom of the bond to this atom (if it has one). Just find the first one
   // if there are many
   //
   // this is hideous. And will be slow on a big molecule. Replace it.
   // When the bonds are calculated, the at_info should store
   // a pointer to a bonded atom (mmdb::Atom *bonded_atom_other_atom)
   //
   std::map<int, int> bonded_atom_other_atom;
   std::map<int, mmdb::Atom *> index_to_atom;
   for (int icol_bond=0; icol_bond<gbc.num_colours; icol_bond++) {
      graphical_bonds_lines_list<graphics_line_t> &ll = gbc.bonds_[icol_bond];
      for (int j=0; j<ll.num_lines; j++) {
         int idx_1 = ll.pair_list[j].atom_index_1;
         int idx_2 = ll.pair_list[j].atom_index_2;
         std::map<int, int>::const_iterator it = bonded_atom_other_atom.find(idx_1);
         if (it == bonded_atom_other_atom.end()) {
            bonded_atom_other_atom[idx_1] = idx_2;
         } else {
            // std::cout << "found other atom bond " << idx_1 << " " << it->first << " " << it->second << std::endl;
         }
         it = bonded_atom_other_atom.find(idx_2);
         if (it == bonded_atom_other_atom.end()) {
            bonded_atom_other_atom[idx_2] = idx_1;
         } else {
            //std::cout << "found other atom bond " << idx_2 << " " << it->first << " " << it->second << std::endl;
         }
      }
   }
   for (int icol=0; icol<gbc.n_consolidated_atom_centres; icol++) {
      for (unsigned int i=0; i<gbc.consolidated_atom_centres[icol].num_points; i++) {
         const graphical_bonds_atom_info_t &at_info = gbc.consolidated_atom_centres[icol].points[i];
         index_to_atom[at_info.atom_index] = at_info.atom_p;
      }
   }

   if (false) {
      for (const auto &item : index_to_atom) {
         std::cout << "   index_to_atom: " << item.first << " " << item.second << " " << coot::atom_spec_t(item.second)
                   << std::endl;
      }
   }

   glm::mat4 unit(1.0);
   for (int icol=0; icol<gbc.n_consolidated_atom_centres; icol++) {
      glm::vec4 col = colour_table[icol];
      for (unsigned int i=0; i<gbc.consolidated_atom_centres[icol].num_points; i++) {
         const graphical_bonds_atom_info_t &at_info = gbc.consolidated_atom_centres[icol].points[i];
         bool do_it = true;
         if (atoms_have_bigger_radius_than_bonds) do_it = false; // not as hemispheres, that is
         mmdb::Atom *at = at_info.atom_p;

         if (do_it) {
            if (at) {
               int state = -1;
               at->GetUDData(udd_handle_bonded_type, state);
               if (state == graphical_bonds_container::NO_BOND) {
                  do_it = false;
               }
            }
         }

         // Oh dear! I need to know where the atom at the other end of the bond is!
         // Which bond? Any bond. The first one. We need to do some pre-processing to know that.
         // That's done now.

         // There is an opportunity for missing atoms here. Things that are not NO_BOND but have not correctly
         // set the atom_index

         if (do_it) {

            unsigned int idx_base = vertices.size();
            unsigned int idx_tri_base = triangles.size();

            float scale = 1.0;
            if (at_info.is_hydrogen_atom) scale *= 0.5;
            glm::vec3 t = cartesian_to_glm(at_info.position); // (at->x, at->y, at->z);
            float sar = scale * atom_radius;
            glm::vec3 sc(sar, sar, sar);

            int atom_index = at_info.atom_index;
            std::map<int, int>::const_iterator it = bonded_atom_other_atom.find(atom_index);
            if (it != bonded_atom_other_atom.end()) {
               // check the index before using it?
               int other_atom_index = it->second;
               if (false)
                  std::cout << "here with other_atom_index " << other_atom_index << " with index_to_atom size "
                            << index_to_atom.size() << std::endl;
               mmdb::Atom *other_at = index_to_atom[other_atom_index];
               if (other_at) {
                  // std::cout << "   other_at " << other_at << " " << coot::atom_spec_t(other_at) << std::endl;
                  glm::vec3 other_atom_pos(other_at->x, other_at->y, other_at->z);
                  glm::mat4 mm = get_octahemi_matrix(t, other_atom_pos, bond_radius); // a rotation matrix

                  std::vector<s_generic_vertex> local_vertices(octasphere_geom.first.size());

                  for (unsigned int ii=0; ii<local_vertices.size(); ii++) {
                     glm::vec3 p0(octasphere_geom.first[ii]);
                     glm::vec4 p1(p0, 1.0);
                     glm::vec4 p2 = mm * p1;
                     glm::vec3 p3(p2);
                     glm::vec3 p4 = p3 * sc;
                     glm::vec3 p5 = p4 + t;
                     auto &vert = local_vertices[ii];
                     vert = s_generic_vertex(p5, p3, col);
                     if (false)
                        std::cout << "make_graphical_bonds_hemispherical_atoms() vertex " << ii << " of " << local_vertices.size() << " "
                                  << coot::atom_spec_t(at) << " vertex position " << glm::to_string(p5) << std::endl;
                  }

                  vertices.insert(vertices.end(), local_vertices.begin(), local_vertices.end());
                  triangles.insert(triangles.end(), octasphere_geom.second.begin(), octasphere_geom.second.end());
                  for (unsigned int k=idx_tri_base; k<triangles.size(); k++)
                     triangles[k].rebase(idx_base);
               } else {
                  std::cout << "WARNING:: other_at is null for atom_index " << atom_index << " " << coot::atom_spec_t(at) << std::endl;
               }
            } else {
               // 20220220-PE Is there an error that causes us to get here? I think so, but I have other things to do so
               // inhibit this message for the moment.
               if (false)
                  std::cout << "DEBUG:: failed to find other_atom for atom index " << atom_index << " " << coot::atom_spec_t(at)
                            << " in bonded_atom_other_atom map of size " << bonded_atom_other_atom.size()  << std::endl;
            }
         }
      }
   }

   // std::cout << " End of hemispheres" << std::endl;
   // for (unsigned int i = 0; i < vertices.size(); i++) {
   //    std::cout << i << " " << glm::to_string(vertices[i].pos) << std::endl;
   // }
}

void
Mesh::make_graphical_bonds_bonds(const graphical_bonds_container &gbc,
                                 float bond_radius,
                                 unsigned int n_slices,
                                 unsigned int n_stacks,
                                 const std::vector<glm::vec4> &colour_table) {

   // unsigned int n_slices = 8;
   // unsigned int n_stacks = 2; // try 1 later.

   // do we need 4 version of this for the end caps? None, left only, right only and both.

   auto get_bond_matrix = [] (const glm::vec3 &pos_1, const glm::vec3 &pos_2, float radius) {
                             glm::vec3 delta = pos_2 - pos_1;
                             float l_delta = glm::distance(pos_2, pos_1);
                             glm::mat4 u(1.0f);
                             glm::mat4 sc = glm::scale(u, glm::vec3(radius, radius, l_delta));
                             // orient
                             glm::vec3 normalized_bond_orientation(glm::normalize(delta));
                             glm::mat4 ori = glm::orientation(normalized_bond_orientation, glm::vec3(0.0, 0.0, 1.0)); // nice
                             // translate
                             glm::mat4 t = glm::translate(u, pos_1);
                             glm::mat4 m = t * ori * sc;
                             return m;
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

   std::pair<glm::vec3, glm::vec3> pp(glm::vec3(0,0,0), glm::vec3(0,0,1));
   is_instanced = false;

   // ----------------------- add the vertices and triangles ----------------------

   // std::cout << "debug:: in make_graphical_bonds_bonds() there are " << gbc.num_colours << " colours " << std::endl;
   for (int icol=0; icol<gbc.num_colours; icol++) {
      glm::vec4 col = colour_table[icol];
      graphical_bonds_lines_list<graphics_line_t> &ll = gbc.bonds_[icol];
      // std::cout << "debug:: in make_graphical_bonds_bonds() colour icol " << icol  << " has " << ll.num_lines << " bond lines " << std::endl;
      for (int j=0; j<ll.num_lines; j++) {
         const coot::Cartesian &start  = ll.pair_list[j].positions.getStart();
         const coot::Cartesian &finish = ll.pair_list[j].positions.getFinish();
         bool thin = ll.thin_lines_flag;
         float bond_radius_this = bond_radius;
         if (thin)
            bond_radius_this *= 0.5;
         if (ll.pair_list[j].cylinder_class == graphics_line_t::KEK_DOUBLE_BOND_INNER_BOND)
            bond_radius_this *= 0.8;
         float bl = ll.pair_list[j].positions.amplitude();
         glm::vec3 pos_1(start.x(),   start.y(),  start.z());
         glm::vec3 pos_2(finish.x(), finish.y(), finish.z());
         // glm::mat4 mm = get_bond_matrix(pos_1, pos_2, bond_radius);
         if (false)
            std::cout << "making bond between " << glm::to_string(pos_1) << " " << glm::to_string(pos_2) << " width " << bond_radius_this
                      << std::endl;
         cylinder cc(std::make_pair(pos_1, pos_2), bond_radius_this, bond_radius_this, bl, col, n_slices, n_stacks);
         cc.set_unstubby_rounded_cap_factor(1.0);

         // we can't add the octasphere caps yet because they appear in the middle of double bonds
         // there should be a means to say "no end caps at all"  - or cap with a hemisphere or
         // cap flat (the end of a double bond).

         if (ll.pair_list[j].has_begin_cap) {
            cc.add_flat_start_cap();
         } else {
            // cc.add_octahemisphere_start_cap();
         }
         if (ll.pair_list[j].has_end_cap) {
            cc.add_flat_end_cap();
         } else {
            // cc.add_octahemisphere_end_cap();
         }
         unsigned int idx_base = vertices.size();
         unsigned int idx_tri_base = triangles.size();
         std::vector<s_generic_vertex> converted_vertices = vnc_vertex_vector_to_generic_vertex_vector(cc.vertices);
         vertices.insert(vertices.end(), converted_vertices.begin(), converted_vertices.end());
         triangles.insert(triangles.end(), cc.triangles.begin(), cc.triangles.end());
         for (unsigned int k=idx_tri_base; k<triangles.size(); k++)
            triangles[k].rebase(idx_base);
      }
   }

   // std::cout << " End of bonds" << std::endl;
   // for (unsigned int i = 0; i < vertices.size(); i++) {
   //    std::cout << i << " " << glm::to_string(vertices[i].pos) << std::endl;
   // }
}


void
Mesh::make_graphical_bonds_rama_balls(const graphical_bonds_container &gbc,
                                      const glm::vec3 &screen_up_dir) {

   auto cartesian_to_glm = [] (const coot::Cartesian &c) {
                              return glm::vec3(c.x(), c.y(), c.z());
                           };

   auto prob_raw_to_colour_rotation = [] (float prob) {
                                         if (prob > 0.5) prob = 0.5; // 0.4 and 2.5 f(for q) might be better (not tested)
                                         // good probabilities have q = 0
                                         // bad probabilities have q 0.66
                                         double q = (1.0 - 2.0 * prob);
                                         q = pow(q, 20.0);
                                         return q;
                                      };

   if (true) { // maybe I want an "is-intermediate-atoms" test here (and rama restraints are not used)?

      unsigned int num_subdivisions = 2;
      glm::vec3 origin(0,0,0);
      float rama_ball_pos_offset_scale = 0.6;
      float rama_ball_radius = 0.5;

      std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > octaball = tessellate_octasphere(num_subdivisions);

      for (int i=0; i<gbc.n_ramachandran_goodness_spots; i++) {
         const coot::Cartesian &position = gbc.ramachandran_goodness_spots_ptr[i].first;
         const float &prob_raw    = gbc.ramachandran_goodness_spots_ptr[i].second;
         double q = prob_raw_to_colour_rotation(prob_raw);
         coot::colour_holder col = coot::colour_holder(q, 0.0, 1.0, false, std::string(""));
         col.scale_intensity(0.6); // calm down the otherwise super-bright Rama ball colours
         glm::vec3 atom_position = cartesian_to_glm(position);
         glm::vec3 ball_position = atom_position + rama_ball_pos_offset_scale * screen_up_dir;
         unsigned int idx_base = vertices.size();
         unsigned int idx_tri_base = triangles.size();
         for (unsigned int ibv=0; ibv<octaball.first.size(); ibv++) {
            glm::vec4 col_v4(col.red, col.green, col.blue, 1.0f);
            const glm::vec3 &vertex_position = octaball.first[ibv];
            s_generic_vertex vertex(ball_position + rama_ball_radius * vertex_position, vertex_position, col_v4);
            vertices.push_back(vertex);
         }
         std::vector<g_triangle> octaball_triangles = octaball.second;
         triangles.insert(triangles.end(), octaball_triangles.begin(), octaball_triangles.end());

         for (unsigned int k=idx_tri_base; k<triangles.size(); k++)
            triangles[k].rebase(idx_base);
      }

   }
}

#include "utils/dodec.hh"

void
Mesh::make_graphical_bonds_rotamer_dodecs(const graphical_bonds_container &gbc,
                                          const glm::vec3 &screen_up_dir) {

   auto cartesian_to_glm = [] (const coot::Cartesian &c) {
                              return glm::vec3(c.x(), c.y(), c.z());
                           };
   auto clipper_to_glm = [] (const clipper::Coord_orth &c) {
                              return glm::vec3(c.x(), c.y(), c.z());
                           };

   auto colour_holder_to_glm = [] (const coot::colour_holder &ch) {
                                  return glm::vec4(ch.red, ch.green, ch.blue, 1.0f);
                               };

   if (false)
      std::cout << "in make_graphical_bonds_rotamer_dodecs() we have " << gbc.n_rotamer_markups
                << " rotamer markups" << std::endl;

   glm::vec4 col(0.6, 0.2, 0.8, 1.0);
   if (gbc.n_rotamer_markups > 0) {

      dodec d;
      std::vector<clipper::Coord_orth> coords = d.coords();
      std::vector<glm::vec3> dodec_postions(coords.size());
      for (unsigned int i=0; i<coords.size(); i++)
         dodec_postions[i] = clipper_to_glm(coords[i]);

      std::vector<s_generic_vertex> dodec_vertices;
      std::vector<g_triangle> dodec_triangles;
      dodec_triangles.reserve(36);

      for (unsigned int iface=0; iface<12; iface++) {

         std::vector<s_generic_vertex> face_verts;
         std::vector<g_triangle> face_triangles;
         face_triangles.reserve(3);

         std::vector<unsigned int> indices_for_face = d.face(iface);
         glm::vec3 ns(0,0,0);
         for (unsigned int j=0; j<5; j++)
            ns += dodec_postions[indices_for_face[j]];
         glm::vec3 normal = glm::normalize(ns);

         for (unsigned int j=0; j<5; j++) {
            glm::vec3 &pos = dodec_postions[indices_for_face[j]];
            s_generic_vertex v(0.5f * pos, normal, col);
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

#if 0
      // Just one - as a test.
      //
      unsigned int idx_base = vertices.size();
      unsigned int idx_tri_base = triangles.size();
      vertices.insert(vertices.end(), dodec_vertices.begin(), dodec_vertices.end());
      triangles.insert(triangles.end(), dodec_triangles.begin(), dodec_triangles.end());
      for (unsigned int jj=idx_tri_base; jj<triangles.size(); jj++)
         triangles[jj].rebase(idx_base);
#endif

      // now there is a dodec at the origin, dodec_vertices and dodec_triangle

      // let's make copies of that and move them around the protein

      for (int i=0; i<gbc.n_rotamer_markups; i++) {
         const rotamer_markup_container_t &rm = gbc.rotamer_markups[i];
         glm::vec3 atom_pos = cartesian_to_glm(rm.pos);

         std::vector<s_generic_vertex> this_dodec_vertices = dodec_vertices; // at the origin to start
         // now move it.
         for (unsigned int j=0; j<dodec_vertices.size(); j++) {
            auto rm_col = rm.col;
            rm_col.scale_intensity(0.6);
            this_dodec_vertices[j].pos  += atom_pos;
            this_dodec_vertices[j].pos  += 1.5f * screen_up_dir;
            this_dodec_vertices[j].color = colour_holder_to_glm(rm_col);
         }
         unsigned int idx_base = vertices.size();
         unsigned int idx_tri_base = triangles.size();
         vertices.insert(vertices.end(), this_dodec_vertices.begin(), this_dodec_vertices.end());
         triangles.insert(triangles.end(), dodec_triangles.begin(), dodec_triangles.end());
         for (unsigned int jj=idx_tri_base; jj<triangles.size(); jj++)
            triangles[jj].rebase(idx_base);
      }
   }
}

#include "molecular-mesh-generator.hh"

void
Mesh::make_graphical_bonds_cis_peptides(const graphical_bonds_container &gbc) {

   auto cartesian_to_glm = [] (const coot::Cartesian &c) {
                              return glm::vec3(c.x(), c.y(), c.z());
                           };
   for (int i=0; i<gbc.n_cis_peptide_markups; i++) {
      const graphical_bonds_cis_peptide_markup &m = gbc.cis_peptide_markups[i];
      int single_model_view_current_model_number = 0; // 20220220-PE hack - this number comes from
                                                      // a graphics_info_t static, so I suppose
                                                      // that it should be passed here.
      if ((single_model_view_current_model_number == 0) ||
          (single_model_view_current_model_number == m.model_number)) {

         std::vector<glm::vec3> glm_quad(4);
         glm_quad[0] = cartesian_to_glm(m.pt_ca_1);
         glm_quad[1] = cartesian_to_glm(m.pt_c_1);
         glm_quad[2] = cartesian_to_glm(m.pt_n_2);
         glm_quad[3] = cartesian_to_glm(m.pt_ca_2);
         molecular_mesh_generator_t mmg;
         coot::util::cis_peptide_quad_info_t::type_t type = coot::util::cis_peptide_quad_info_t::CIS;
         if (m.is_pre_pro_cis_peptide) type = coot::util::cis_peptide_quad_info_t::PRE_PRO_CIS;
         if (m.is_twisted)             type = coot::util::cis_peptide_quad_info_t::TWISTED_TRANS;
         std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > cpg =
            mmg.make_cis_peptide_geom(glm_quad, type);

         unsigned int idx_base = vertices.size();
         unsigned int idx_tri_base = triangles.size();

         vertices.insert(vertices.end(), cpg.first.begin(), cpg.first.end());
         triangles.insert(triangles.end(), cpg.second.begin(), cpg.second.end());
         for (unsigned int jj=idx_tri_base; jj<triangles.size(); jj++)
            triangles[jj].rebase(idx_base);
      }
   }
}




// the simple-lines option for the main molecule
void
Mesh::make_bond_lines(const graphical_bonds_container &bonds_box, const std::vector<glm::vec4> &colour_table) {

   // this is "on the fly" - i.e. does a built-in "setup_buffers" and throws away the vertices when done

   if (first_time)
      glGenVertexArrays(1, &vao);

   glBindVertexArray(vao);

   n_simple_bond_lines = 0;
   for (int icol=0; icol<bonds_box.num_colours; icol++) {
      glm::vec4 col = colour_table[icol];
      graphical_bonds_lines_list<graphics_line_t> &ll = bonds_box.bonds_[icol];
      n_simple_bond_lines += ll.num_lines;
   }

   std::vector<simple_atoms_line_vertex> vertices_for_atoms_as_lines;
   vertices_for_atoms_as_lines.reserve(n_simple_bond_lines);

   for (int icol=0; icol<bonds_box.num_colours; icol++) {
      glm::vec4 col = colour_table[icol];
      graphical_bonds_lines_list<graphics_line_t> &ll = bonds_box.bonds_[icol];
      for (int j=0; j<ll.num_lines; j++) {
         const coot::Cartesian &start  = ll.pair_list[j].positions.getStart();
         const coot::Cartesian &finish = ll.pair_list[j].positions.getFinish();
         glm::vec3 pos_1(start.x(),   start.y(),  start.z());
         glm::vec3 pos_2(finish.x(), finish.y(), finish.z());
         vertices_for_atoms_as_lines.push_back(simple_atoms_line_vertex(pos_1, col));
         vertices_for_atoms_as_lines.push_back(simple_atoms_line_vertex(pos_2, col));
      }
   }

   if (first_time) {
   } else {
      glDeleteBuffers(1, &buffer_id);
   }

   unsigned int buffer_size = n_simple_bond_lines * 2 * sizeof(simple_atoms_line_vertex);
   glGenBuffers(1, &buffer_id);
   glBindBuffer(GL_ARRAY_BUFFER, buffer_id);
   glBufferData(GL_ARRAY_BUFFER, buffer_size, &(vertices_for_atoms_as_lines[0]), GL_STATIC_DRAW);

   glEnableVertexAttribArray(0); // position
   glEnableVertexAttribArray(1); // colour

   glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(simple_atoms_line_vertex), 0);
   glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(simple_atoms_line_vertex), reinterpret_cast<void *>(sizeof(glm::vec3)));

   // done
   glBindVertexArray(0); // we don't need to disable the vertexattribarray if we unbind the vertex array
   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: Mesh::make_bond_lines() check-before-return error " << err << std::endl;

}

// the draw function draw_simple_bond_lines() is in Mesh.cc
