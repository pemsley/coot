
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>  // to_string()
#include <glm/gtx/rotate_vector.hpp>
#include <glm/ext.hpp>

#include "coot_molecule.hh"
#include "coot_molecule_bonds.hh"
#include "coot-utils/oct.hh"
#include "coot-utils/cylinder.hh"

void
make_instanced_graphical_bonds_spherical_atoms(coot::instanced_mesh_t &m, // fill this
                                     const graphical_bonds_container &gbc,
                                     int bonds_box_type,
                                     int udd_handle_bonded_type,
                                     float atom_radius,
                                     float bond_radius,
                                     unsigned int num_subdivisions,
                                     const std::vector<glm::vec4> &colour_table) {

   // 20230114-PE
   // copied and edited from from src/Mesh-from-graphical-bonds-instanced.cc

   coot::instanced_geometry_t ig;

   bool atoms_have_bigger_radius_than_bonds = false;
   if (atom_radius > bond_radius) atoms_have_bigger_radius_than_bonds = true;

   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > octaphere_geom =
      tessellate_octasphere(num_subdivisions);

   // ----------------------- setup the vertices and triangles ---------------------

   std::vector<coot::api::vn_vertex> local_vertices(octaphere_geom.first.size());
   for (unsigned int i=0; i<octaphere_geom.first.size(); i++) {
      const glm::vec3 &v(octaphere_geom.first[i]);
      local_vertices[i] = coot::api::vn_vertex(v, v);
   }
   ig.vertices = local_vertices;
   ig.triangles = octaphere_geom.second;

   // ----------------------- setup the instances ----------------------

   int cts = colour_table.size(); // change type   
   for (int icol=0; icol<gbc.n_consolidated_atom_centres; icol++) {
      // from src/Mesh-from-graphical-bonds-instanced.cc: glm::vec4 col = get_glm_colour_for_bonds_func(icol, bonds_box_type);
      glm::vec4 col(0.4, 0.4, 0.4, 1.0);
      if (icol<cts)
         col = colour_table[icol];
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
            float scale = 1.0;
            if (at_info.is_hydrogen_atom) scale *= 0.5;
            glm::vec3 t(at->x, at->y, at->z);
            float sar = scale * atom_radius;
            glm::vec3 sc(sar, sar, sar);
            coot::instancing_data_type_A_t idA(t, col, sc);
            ig.instancing_data_A.push_back(idA);
         }
      }
   }

   m.geom.push_back(ig);

}


void
make_instanced_graphical_bonds_hemispherical_atoms(coot::instanced_mesh_t &m, // fill this
                                     const graphical_bonds_container &gbc,
                                     int bonds_box_type,
                                     int udd_handle_bonded_type,
                                     float atom_radius,
                                     float bond_radius,
                                     unsigned int num_subdivisions,
                                     const std::vector<glm::vec4> &colour_table) {

   // copied and edited from Mesh::make_graphical_bonds_hemispherical_atoms_instanced_version

   coot::instanced_geometry_t ig;

   bool atoms_have_bigger_radius_than_bonds = false;
   if (atom_radius > bond_radius)
      atoms_have_bigger_radius_than_bonds = true;

   // like above, different z axis because we want the hemisphere to extend outside the cylinder - and we don't need to
   // scale to bond length
   auto get_octahemi_matrix = [] (const glm::vec3 &pos_1, const glm::vec3 &pos_2, float radius) {
                                 glm::vec3 delta = pos_2 - pos_1;
                                 glm::mat4 u(1.0f);
                                 glm::mat4 sc = glm::scale(u, glm::vec3(radius, radius, radius));
                                 // orient
                                 glm::vec3 normalized_bond_orientation(glm::normalize(delta));
                                 glm::mat4 ori = glm::orientation(normalized_bond_orientation, glm::vec3(0.0, 0.0, -1.0));
                                 // translate
                                 glm::mat4 t = glm::translate(u, pos_1);
                                 glm::mat4 m = t * ori * sc;
                                 return m;
                              };

   // ----------------------- setup the vertices and triangles ---------------------
   
   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > octaphere_geom =
      tessellate_hemisphere_patch(num_subdivisions);

   std::vector<coot::api::vn_vertex> local_vertices(octaphere_geom.first.size());
   for (unsigned int i=0; i<octaphere_geom.first.size(); i++) {
      const glm::vec3 &v(octaphere_geom.first[i]);
      local_vertices[i] = coot::api::vn_vertex(v, v);
   }
   ig.vertices = local_vertices;
   ig.triangles = octaphere_geom.second;

   // ----------------------- setup the instances ----------------------

   std::vector<glm::mat4> instanced_matrices;
   std::vector<glm::vec4> instanced_colours;

   glm::mat4 unit(1.0);
   int cts = colour_table.size(); // change type   
   for (int icol=0; icol<gbc.n_consolidated_atom_centres; icol++) {
      // glm::vec4 col = get_glm_colour_for_bonds_func(icol, bonds_box_type);
      glm::vec4 col(0.4, 0.4, 0.4, 1.0);
      if (icol<cts) // it will be of course!
         col = colour_table[icol];
      for (unsigned int i=0; i<gbc.consolidated_atom_centres[icol].num_points; i++) {
         const graphical_bonds_atom_info_t &at_info = gbc.consolidated_atom_centres[icol].points[i];
         mmdb::Atom *at = at_info.atom_p;
         glm::vec3 t(at->x, at->y, at->z);
         glm::mat4 ori(1.0); // 20230114-PE needs fixing.
         float scale = 1.0;
         if (at_info.is_hydrogen_atom) scale *= 0.5;
         float sar = scale * atom_radius;
         glm::vec3 sc(sar, sar, sar);
         coot::instancing_data_type_B_t id(t, col, sc, ori);
         ig.instancing_data_B.push_back(id);
      }
   }

   m.geom.push_back(ig);
}

void
make_instanced_graphical_bonds_bonds(coot::instanced_mesh_t &m,
                           const graphical_bonds_container &gbc,
                           float bond_radius,
                           unsigned int n_slices,
                           unsigned int n_stacks,
                           const std::vector<glm::vec4> &colour_table) {

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

   auto convert_vertices = [] (const std::vector<coot::api::vnc_vertex> &v_in) {
      std::vector<coot::api::vn_vertex> v_out;
      return v_out;
   };

   // 20230114-PE
   // copied and edited from src/Mesh::Mesh-from-graphical-bonds-instanced.cc
   // make_graphical_bonds_bonds_instanced_version()

   coot::instanced_geometry_t ig_00;
   coot::instanced_geometry_t ig_01;
   coot::instanced_geometry_t ig_10;
   coot::instanced_geometry_t ig_11;

   // ----------------------- setup the vertices and triangles ----------------------

   std::pair<glm::vec3, glm::vec3> pp(glm::vec3(0,0,0), glm::vec3(0,0,1));
   cylinder c_00(pp, 1.0, 1.0, 1.0, n_slices, n_stacks);
   cylinder c_01(pp, 1.0, 1.0, 1.0, n_slices, n_stacks);
   cylinder c_10(pp, 1.0, 1.0, 1.0, n_slices, n_stacks); // possibly none of these actually
   cylinder c_11(pp, 1.0, 1.0, 1.0, n_slices, n_stacks);
   c_01.add_flat_end_cap();   // this and below might be the wrong way around
   c_10.add_flat_start_cap();
   c_11.add_flat_end_cap();
   c_11.add_flat_end_cap();

   ig_00.vertices = convert_vertices(c_00.vertices);
   ig_01.vertices = convert_vertices(c_01.vertices);
   ig_10.vertices = convert_vertices(c_10.vertices);
   ig_11.vertices = convert_vertices(c_11.vertices);

   ig_00.triangles = c_00.triangles;
   ig_01.triangles = c_01.triangles;
   ig_10.triangles = c_10.triangles;
   ig_11.triangles = c_11.triangles;

   int cts = colour_table.size(); // change type   
   for (int icol=0; icol<gbc.num_colours; icol++) {

      glm::vec4 col(0.4, 0.4, 0.4, 1.0);
      if (icol<cts) // it will be of course!
         col = colour_table[icol];
      graphical_bonds_lines_list<graphics_line_t> &ll = gbc.bonds_[icol];
      for (int j=0; j<ll.num_lines; j++) {
         const coot::Cartesian &start  = ll.pair_list[j].positions.getStart();
         const coot::Cartesian &finish = ll.pair_list[j].positions.getFinish();
         float bl = ll.pair_list[j].positions.amplitude();
         glm::vec3 pos_1(start.x(),   start.y(),  start.z());
         glm::vec3 pos_2(finish.x(), finish.y(), finish.z());
         glm::mat4 ori = get_bond_matrix(pos_1, pos_2, bond_radius);
         float scale = 1.0;
         if (ll.thin_lines_flag) scale *= 0.5;
         float sar = scale * bond_radius;
         glm::vec3 sc(sar, sar, bl);
         coot::instancing_data_type_B_t id(pos_1, col, sc, ori); // perhaps use pos_1?
         int cappiness = 0;
         if (ll.pair_list[j].has_begin_cap) cappiness += 1;
         if (ll.pair_list[j].has_end_cap)   cappiness += 2;

         if (cappiness == 0) ig_00.instancing_data_B.push_back(id);
         if (cappiness == 1) ig_01.instancing_data_B.push_back(id);
         if (cappiness == 2) ig_10.instancing_data_B.push_back(id);
         if (cappiness == 3) ig_11.instancing_data_B.push_back(id);

      }
   }

   if (false) {
      // there  are no "_10" bonds cylinders
      std::cout << "debug:: bonds: 00 " << ig_00.instancing_data_B.size() << std::endl;
      std::cout << "debug:: bonds: 01 " << ig_01.instancing_data_B.size() << std::endl;
      std::cout << "debug:: bonds: 10 " << ig_10.instancing_data_B.size() << std::endl;
      std::cout << "debug:: bonds: 11 " << ig_11.instancing_data_B.size() << std::endl;
   }

   if (! ig_00.instancing_data_B.empty()) m.geom.push_back(ig_00);
   if (! ig_01.instancing_data_B.empty()) m.geom.push_back(ig_01);
   if (! ig_10.instancing_data_B.empty()) m.geom.push_back(ig_10);
   if (! ig_11.instancing_data_B.empty()) m.geom.push_back(ig_11);

}



coot::instanced_mesh_t
coot::molecule_t::get_bonds_mesh_instanced(const std::string &mode, coot::protein_geometry *geom,
                                           bool against_a_dark_background, float bonds_width, float atom_radius_to_bond_width_ratio,
                                           int smoothness_factor,
                                           bool draw_hydrogen_atoms_flag,
                                           bool draw_missing_residue_loops_flag) {

   coot::instanced_mesh_t m;

   // 20230113-PE there is a hunk here that is the same as get_bonds_mesh()

   //
   float bond_radius = bonds_width;
   float atom_radius = bond_radius;
   if (atom_radius_to_bond_width_ratio > 1.0)
      atom_radius = bond_radius * atom_radius_to_bond_width_ratio;

   unsigned int num_subdivisions = 1;
   unsigned int n_slices = 8;
   unsigned int n_stacks = 2;
   // int representation_type = BALL_AND_STICK;

   if (smoothness_factor == 2) {
      num_subdivisions = 2;
      n_slices = 16; // was 18
   }

   if (smoothness_factor == 3) {
      num_subdivisions = 3;
      n_slices = 32;
   }

   bonds_box_type = coot::COLOUR_BY_CHAIN_BONDS;

   if (mode == "CA+LIGANDS") {
      // something
   }

   std::set<int> no_bonds_to_these_atoms; // empty
   // we don't make rotamer dodecs in this function
   makebonds(geom, nullptr, no_bonds_to_these_atoms, draw_hydrogen_atoms_flag, draw_missing_residue_loops_flag); // this makes the bonds_box.

   // get the udd_handle_bonded_type after making the bonds (because the handle is made by making the bond)
   int udd_handle_bonded_type = atom_sel.mol->GetUDDHandle(mmdb::UDR_ATOM, "found bond");
   if (udd_handle_bonded_type == mmdb::UDDATA_WrongUDRType) {
      std::cout << "ERROR:: in get_bonds_mesh() wrong udd data type " << udd_handle_bonded_type << std::endl;
      return m;
   } else {
      // std::cout << "debug:: OK, udd_handle_bonded_type is " << udd_handle_bonded_type
      // << " not " << mmdb::UDDATA_WrongUDRType << std::endl;
   }

   std::vector<glm::vec4> colour_table = make_colour_table(against_a_dark_background);
   if (colour_table.empty()) {
      std::cout << "ERROR:: you need to make the bonds before getting the bonds mesh" << std::endl;
   }

   const graphical_bonds_container &gbc = bonds_box; // alias because it's named like that in Mesh-from-graphical-bonds

   unsigned int n_bonds = 0;
   for (int icol_bond=0; icol_bond<gbc.num_colours; icol_bond++) {
      graphical_bonds_lines_list<graphics_line_t> &ll = gbc.bonds_[icol_bond];
      n_bonds += ll.num_lines;
   }

   make_instanced_graphical_bonds_spherical_atoms(m, gbc, bonds_box_type, udd_handle_bonded_type,
                                        atom_radius, bond_radius, num_subdivisions, colour_table);
   make_instanced_graphical_bonds_hemispherical_atoms(m, gbc, bonds_box_type, udd_handle_bonded_type, atom_radius,
                                            bond_radius, num_subdivisions, colour_table);

   make_instanced_graphical_bonds_bonds(m, gbc, bond_radius, n_slices, n_stacks, colour_table);

   make_graphical_bonds_cis_peptides(m.markup, gbc);

   return m;

}
