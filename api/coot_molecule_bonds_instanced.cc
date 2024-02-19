
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>  // to_string()
#include <glm/gtx/rotate_vector.hpp>
#include <glm/ext.hpp>

#include "coot_molecule.hh"
#include "coot_molecule_bonds.hh"
#include "coot-utils/oct.hh"
#include "coot-utils/cylinder.hh"

// Atom radii are limited to 2.0
void
make_instanced_graphical_bonds_spherical_atoms(coot::instanced_mesh_t &m, // add to this
                                               const graphical_bonds_container &gbc,
                                               coot::api_bond_colour_t bonds_box_type, // remove these one day
                                               float base_atom_radius,
                                               float base_bond_radius,
                                               unsigned int num_subdivisions,
                                               const std::vector<glm::vec4> &colour_table) {


   // 20230114-PE
   // copied and edited from from src/Mesh-from-graphical-bonds-instanced.cc

   coot::instanced_geometry_t ig("spherical-atoms");

   bool atoms_have_bigger_radius_than_bonds = false;
   if (base_atom_radius > base_bond_radius) atoms_have_bigger_radius_than_bonds = true;

   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > octosphere_geom =
      tessellate_octasphere(num_subdivisions);

   // ----------------------- setup the vertices and triangles ---------------------

   std::vector<coot::api::vn_vertex> local_vertices(octosphere_geom.first.size());
   for (unsigned int i=0; i<octosphere_geom.first.size(); i++) {
      const glm::vec3 &v(octosphere_geom.first[i]);
      local_vertices[i] = coot::api::vn_vertex(v, v);
   }
   ig.vertices = local_vertices;
   ig.triangles = octosphere_geom.second;

   // ----------------------- setup the instances ----------------------

   int cts = colour_table.size(); // change type   
   for (int icol=0; icol<gbc.n_consolidated_atom_centres; icol++) {
      // from src/Mesh-from-graphical-bonds-instanced.cc: glm::vec4 col = get_glm_colour_for_bonds_func(icol, bonds_box_type);
      glm::vec4 col(0.4, 0.4, 0.4, 1.0);
      if (icol<cts)
         col = colour_table[icol];
      for (unsigned int i=0; i<gbc.consolidated_atom_centres[icol].num_points; i++) {
         const graphical_bonds_atom_info_t &at_info = gbc.consolidated_atom_centres[icol].points[i];
         mmdb::Atom *at = at_info.atom_p;
         // bool do_it = atoms_have_bigger_radius_than_bonds;

         // if (! do_it) {
         //    int state = -1;
         //    at->GetUDData(udd_handle_bonded_type, state);
         //    if (state == graphical_bonds_container::NO_BOND) {
         //       do_it = true;
         //    }
         // }

         bool do_it = true;  // everything is spherical for the moment.

         if (do_it) {

            // radius_scale is 2 for waters and 4 for metals
            // 20240218-PE base_atom_radius is typically 0.12 but can be 1.67 for "Goodsell" model.
            // 4 * 1.67 is 6.68 and that is too big. So let's just add a limit to the size of sar
            float scale = at_info.radius_scale;
            float sar = scale * base_atom_radius;
            if (sar > 2.2) sar = 2.2; // atom radius limit
            // 20231113-PE should I check for waters for this limit?
            if (at_info.is_water)
               if (sar > 0.65) sar = 0.65f;
            glm::vec3 sc(sar, sar, sar);
            glm::vec3 t(at->x, at->y, at->z);
            coot::instancing_data_type_A_t idA(t, col, sc);
            ig.instancing_data_A.push_back(idA);
            // std::cout << "at: " << coot::atom_spec_t(at) << " scale: " << scale << " sar " << sar << " " << std::endl;
         }
      }
   }

   m.add(ig);

}


void
make_instanced_graphical_bonds_hemispherical_atoms(coot::instanced_mesh_t &m, // add to this
                                                   const graphical_bonds_container &gbc,
                                                   coot::api_bond_colour_t bonds_box_type,
                                                   float atom_radius,
                                                   float bond_radius,
                                                   unsigned int num_subdivisions,
                                                   const std::vector<glm::vec4> &colour_table) {

   return; // 20230224-PE every atom is spherical for the moment.

   // copied and edited from Mesh::make_graphical_bonds_hemispherical_atoms_instanced_version

   coot::instanced_geometry_t ig("hemispherical-atoms");

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
         // 20231113-PE get rid
         // if (at_info.is_water) scale *= 3.33;
         float sar = scale * atom_radius;
         glm::vec3 sc(sar, sar, sar);
         coot::instancing_data_type_B_t id(t, col, sc, ori);
         ig.instancing_data_B.push_back(id);
      }
   }

   m.add(ig);
}

void make_graphical_bonds_spherical_atoms_with_vdw_radii_instanced(coot::instanced_mesh_t &m, const graphical_bonds_container &gbc,
                                                                   unsigned int num_subdivisions,
                                                                   const std::vector<glm::vec4> &colour_table,
                                                                   const coot::protein_geometry &geom) {

   coot::instanced_geometry_t ig("vdW Balls");

   // ----------------------- setup the vertices and triangles ---------------------

   std::pair<std::vector<glm::vec3>, std::vector<g_triangle> > octosphere_geom =
      tessellate_octasphere(num_subdivisions);

   std::vector<coot::api::vn_vertex> local_vertices(octosphere_geom.first.size());
   for (unsigned int i=0; i<octosphere_geom.first.size(); i++) {
      const glm::vec3 &v(octosphere_geom.first[i]);
      local_vertices[i] = coot::api::vn_vertex(v, v);
   }
   ig.vertices = local_vertices;
   ig.triangles = octosphere_geom.second;

   // ----------------------- setup the instances ----------------------

   std::map<std::string, float> ele_to_radius_map;
   glm::mat4 unit(1.0);
   for (int icol=0; icol<gbc.n_consolidated_atom_centres; icol++) {
      glm::vec4 col = colour_table[icol];
      for (unsigned int i=0; i<gbc.consolidated_atom_centres[icol].num_points; i++) {
         const graphical_bonds_atom_info_t &at_info = gbc.consolidated_atom_centres[icol].points[i];
         mmdb::Atom *at = at_info.atom_p;
         std::string ele(at->element);
         std::map<std::string, float>::const_iterator it = ele_to_radius_map.find(ele);
         float atom_radius = 1.0;
         if (it != ele_to_radius_map.end()) {
            atom_radius = it->second;
         } else {
            std::string atom_name(at->GetAtomName());
            std::string residue_name(at->GetResName());
            atom_radius = geom.get_vdw_radius(atom_name, residue_name, coot::protein_geometry::IMOL_ENC_ANY, false);
            ele_to_radius_map[ele] = atom_radius;
         }

         glm::vec3 t(at->x, at->y, at->z);
         glm::vec3 sc(atom_radius, atom_radius, atom_radius);
         coot::instancing_data_type_A_t id(t, col, sc);
         ig.instancing_data_A.push_back(id);
      }
   }
   m.add(ig);
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
                             // 20230116-PE no scaling and translation of the orientation matrix.
                             // we need a rotation matrix because that is used for the normals.
                             // glm::mat4 t = glm::translate(u, pos_1);
                             // glm::mat4 m = t * ori * sc;
                             return ori;
                          };

   auto convert_vertices = [] (const std::vector<coot::api::vnc_vertex> &v_in) {
      std::vector<coot::api::vn_vertex> v_out(v_in.size());
      for (unsigned int i=0; i<v_in.size(); i++) {
         const auto &v = v_in[i];
         v_out[i] = coot::api::vn_vertex(v.pos, v.normal);
      }
      return v_out;
   };

   // 20230114-PE
   // copied and edited from src/Mesh::Mesh-from-graphical-bonds-instanced.cc
   // make_graphical_bonds_bonds_instanced_version()

   coot::instanced_geometry_t ig_00("cylinder-00");
   coot::instanced_geometry_t ig_01("cylinder-01");
   coot::instanced_geometry_t ig_10("cylinder-10");
   coot::instanced_geometry_t ig_11("cylinder-11");

   // ----------------------- setup the vertices and triangles ----------------------

   std::pair<glm::vec3, glm::vec3> pp(glm::vec3(0,0,0), glm::vec3(0,0,1));
   cylinder c_00(pp, 1.0, 1.0, 1.0, n_slices, n_stacks);
   cylinder c_01(pp, 1.0, 1.0, 1.0, n_slices, n_stacks);
   cylinder c_10(pp, 1.0, 1.0, 1.0, n_slices, n_stacks); // possibly none of these actually
   cylinder c_11(pp, 1.0, 1.0, 1.0, n_slices, n_stacks);
   c_01.add_flat_start_cap();
   c_10.add_flat_end_cap();   // 20230122-PE these orientations have now been checked.
   c_11.add_flat_start_cap();
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

      // std::cout << "_bonds_bonds icol " << icol << " has n_lines " << ll.num_lines << std::endl;

      for (int j=0; j<ll.num_lines; j++) {
         const coot::Cartesian &start  = ll.pair_list[j].positions.getStart();
         const coot::Cartesian &finish = ll.pair_list[j].positions.getFinish();
         float bl = ll.pair_list[j].positions.amplitude();
         glm::vec3 pos_1(start.x(),   start.y(),  start.z());
         glm::vec3 pos_2(finish.x(), finish.y(), finish.z());
         glm::mat4 ori = get_bond_matrix(pos_2, pos_1, bond_radius);
         float scale = 1.0;
         if (ll.thin_lines_flag) scale *= 0.5;
         if (ll.pair_list[j].cylinder_class == graphics_line_t::KEK_DOUBLE_BOND_INNER_BOND)
            scale *= 0.7;
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

   if (! ig_00.instancing_data_B.empty()) m.add(ig_00);
   if (! ig_01.instancing_data_B.empty()) m.add(ig_01);
   if (! ig_10.instancing_data_B.empty()) m.add(ig_10);
   if (! ig_11.instancing_data_B.empty()) m.add(ig_11);

}



coot::instanced_mesh_t
coot::molecule_t::get_bonds_mesh_instanced(const std::string &mode, coot::protein_geometry *geom,
                                           bool against_a_dark_background, float bonds_width,
                                           float atom_radius_to_bond_width_ratio,
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

   if (smoothness_factor == 4) {
      num_subdivisions = 4;
      n_slices = 64;
   }


   bonds_box_type = coot::api_bond_colour_t::COLOUR_BY_CHAIN_BONDS;

   std::set<int> no_bonds_to_these_atoms = no_bonds_to_these_atom_indices; // weird that this is not passed.

   if (mode == "CA+LIGANDS") {

      bool do_bonds_to_hydrogens = false;
      Bond_lines_container bonds(geom, "dummy-CA-mode", no_bonds_to_these_atom_indices, do_bonds_to_hydrogens);
      float min_dist = 2.4;
      float max_dist = 4.7;
      bonds.do_Ca_plus_ligands_bonds(atom_sel, imol_no, geom, min_dist, max_dist, draw_missing_residue_loops_flag, draw_hydrogen_atoms_flag);
      bonds_box = bonds.make_graphical_bonds_no_thinning();
      std::vector<glm::vec4> colour_table = make_colour_table(against_a_dark_background);
      make_instanced_graphical_bonds_bonds(m, bonds_box, bond_radius, n_slices, n_stacks, colour_table);
   }

   if (mode == "COLOUR-BY-CHAIN-AND-DICTIONARY") {

      // makebonds() is a trivial wrapper of make_colour_by_chain_bonds(), so just remove the makebonds() call.

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

      // print_colour_table("from get_bonds_mesh_instanced()");
      if (colour_table.empty()) {
         std::cout << "ERROR:: you need to make the bonds before getting the bonds mesh" << std::endl;
      }

      const graphical_bonds_container &gbc = bonds_box; // alias because it's named like that in Mesh-from-graphical-bonds

      make_instanced_graphical_bonds_spherical_atoms(m, gbc, bonds_box_type, atom_radius, bond_radius,
                                                     num_subdivisions, colour_table);
      make_instanced_graphical_bonds_hemispherical_atoms(m, gbc, bonds_box_type, atom_radius,
                                                         bond_radius, num_subdivisions, colour_table);

      make_instanced_graphical_bonds_bonds(m, gbc, bond_radius, n_slices, n_stacks, colour_table);

      make_graphical_bonds_cis_peptides(m.markup, gbc);
   }

   if (mode == "VDW_BALLS" || mode == "VDW-BALLS") {

      // we don't make rotamer dodecs in this function
      // makebonds() is a trivial wrapper of make_colour_by_chain_bonds(), so just remove the makebonds() call.
      makebonds(geom, nullptr, no_bonds_to_these_atoms, draw_hydrogen_atoms_flag, draw_missing_residue_loops_flag); // this makes the bonds_box.
      std::vector<glm::vec4> colour_table = make_colour_table(against_a_dark_background);
      make_graphical_bonds_spherical_atoms_with_vdw_radii_instanced(m, bonds_box, num_subdivisions, colour_table, *geom);

   }

   if (mode == "USER-DEFINED-COLOURS") {

      std::cout << "---------------  in get_bonds_mesh_instanced() C mode is " << mode << " bonds_box_type is " << int(bonds_box_type) << std::endl;
      bool force_rebond = true;
      bool do_rotamer_markup = false; // pass this
      // don't do rotamer dodecs
      int udd_handle_bonded_type = atom_sel.mol->GetUDDHandle(mmdb::UDR_ATOM, "found bond");
      // as above
      if (udd_handle_bonded_type == mmdb::UDDATA_WrongUDRType) {
         std::cout << "ERROR:: in get_bonds_mesh() wrong udd data type " << udd_handle_bonded_type << std::endl;
         return m;
      }
      std::cout << "---------------  in get_bonds_mesh_instanced() D mode is " << mode << " bonds_box_type is " << int(bonds_box_type) << std::endl;
      make_colour_by_chain_bonds(geom, no_bonds_to_these_atoms, true, false, draw_hydrogen_atoms_flag,
                                 draw_missing_residue_loops_flag, do_rotamer_markup, nullptr, force_rebond);
      bonds_box_type = coot::api_bond_colour_t::COLOUR_BY_USER_DEFINED_COLOURS____BONDS;
      std::cout << "---------------  in get_bonds_mesh_instanced() E mode is " << mode << " bonds_box_type is " << int(bonds_box_type) << std::endl;
      std::vector<glm::vec4> colour_table = make_colour_table(against_a_dark_background);
      make_instanced_graphical_bonds_spherical_atoms(m, bonds_box, bonds_box_type, atom_radius, bond_radius,
                                                     num_subdivisions, colour_table);
      make_instanced_graphical_bonds_hemispherical_atoms(m, bonds_box, bonds_box_type,
                                                         atom_radius, bond_radius,
                                                         num_subdivisions, colour_table);
      make_instanced_graphical_bonds_bonds(m, bonds_box, bond_radius, n_slices, n_stacks, colour_table);
      make_graphical_bonds_cis_peptides(m.markup, bonds_box);
   }

   return m;

}


coot::instanced_mesh_t
coot::molecule_t::get_bonds_mesh_for_selection_instanced(const std::string &mode, const std::string &multi_cids,
                                                         coot::protein_geometry *geom, bool against_a_dark_background,
                                                         float bond_radius, float atom_radius_to_bond_width_ratio,
                                                         int  num_subdivisions,
                                                         bool draw_hydrogen_atoms_flag,
                                                         bool draw_missing_residue_loops) {

   auto count_atoms_in_mol = [] (mmdb::Manager *mol) {
      unsigned int n = 0;
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  int n_atoms = residue_p->GetNumberOfAtoms();
                  for (int iat=0; iat<n_atoms; iat++) {
                     mmdb::Atom *at = residue_p->GetAtom(iat);
                     if (! at->isTer()) {
                        n++;
                     }
                  }
               }
            }
         }
      }
      return n;
   };

   coot::instanced_mesh_t m;

   std::vector<std::string> v = coot::util::split_string(multi_cids, "||");

   int udd_atom_selection = atom_sel.mol->NewSelection(); // d
   if (! v.empty()) {
      for (const auto &cid : v) {
         atom_sel.mol->Select(udd_atom_selection, mmdb::STYPE_ATOM, cid.c_str(), mmdb::SKEY_OR);
      }
   }

   // this does atom index transfer
   mmdb::Manager *new_mol = util::create_mmdbmanager_from_atom_selection(atom_sel.mol, udd_atom_selection, false);
   int transfer_atom_index_handle = new_mol->GetUDDHandle(mmdb::UDR_ATOM, "transfer atom index");
   atom_selection_container_t atom_sel_ligand = make_asc(new_mol); // cleared up at end of function
   atom_sel_ligand.UDDAtomIndexHandle = transfer_atom_index_handle;
   atom_sel.mol->DeleteSelection(udd_atom_selection);

   if (false) {
      unsigned int n_atoms = count_atoms_in_mol(new_mol);
      std::cout << "debug:: in get_bonds_mesh_for_selection_instanced() there are " << n_atoms
                << " atoms in the atom selection: " << multi_cids << std::endl;
   }

   apply_user_defined_atom_colour_selections(indexed_user_defined_colour_selection_cids,
                                             indexed_user_defined_colour_selection_cids_apply_to_non_carbon_atoms_also,
                                             atom_sel_ligand.mol);

   float atom_radius = bond_radius * atom_radius_to_bond_width_ratio;

   unsigned int n_slices = 8;
   unsigned int n_stacks = 2;

   if (num_subdivisions == 2) {
      n_slices = 16;
   }

   if (num_subdivisions == 3) {
      n_slices = 32;
   }

   const std::set<int> &no_bonds_to_these_atoms = no_bonds_to_these_atom_indices;
   int udd_handle_bonded_type = atom_sel_ligand.mol->GetUDDHandle(mmdb::UDR_ATOM, "found bond");
   if (udd_handle_bonded_type == mmdb::UDDATA_WrongUDRType) {
      std::cout << "ERROR:: in get_bonds_mesh() wrong udd data type " << udd_handle_bonded_type << std::endl;
      return m;
   }

   bool change_c_only_flag = true;
   bool goodsell_mode = false;
   bool do_rota_markup = false;
   bonds_box_type = coot::api_bond_colour_t::COLOUR_BY_CHAIN_BONDS; // used in colour table?

   if (mode == "COLOUR-BY-CHAIN-AND-DICTIONARY") {

      Bond_lines_container bonds(geom, no_bonds_to_these_atoms, draw_hydrogen_atoms_flag);
      bonds.do_colour_by_chain_bonds(atom_sel_ligand, false, imol_no, draw_hydrogen_atoms_flag,
                                     draw_missing_residue_loops, change_c_only_flag, goodsell_mode,
                                     do_rota_markup);

      // 20240108-PE We need to fill the bonds box before making the colour table.
      //             Why was this not here before?
      bonds_box.clear_up();
      bonds_box = bonds.make_graphical_bonds_no_thinning();

      // std::cout << "------------------------------------- calling make_colour_table() " << std::endl;
      std::vector<glm::vec4> colour_table = make_colour_table(against_a_dark_background);
      // std::cout << "------------------------------------- done make_colour_table() " << std::endl;

      // print_colour_table("from get_bonds_mesh_for_selection_instanced()");

      auto gbc = bonds.make_graphical_bonds();

      make_instanced_graphical_bonds_spherical_atoms(m, gbc, bonds_box_type, atom_radius, bond_radius,
                                                     num_subdivisions, colour_table);
      make_instanced_graphical_bonds_hemispherical_atoms(m, gbc, bonds_box_type, atom_radius, bond_radius,
                                                         num_subdivisions, colour_table);

      make_instanced_graphical_bonds_bonds(m, gbc, bond_radius, n_slices, n_stacks, colour_table);

      make_graphical_bonds_cis_peptides(m.markup, gbc);

      atom_sel_ligand.clear_up();
   }

   if (mode == "CA+LIGANDS") {

      bool do_bonds_to_hydrogens = false;
      Bond_lines_container bonds(geom, "dummy-CA-mode", no_bonds_to_these_atoms, do_bonds_to_hydrogens);
      float min_dist = 2.4;
      float max_dist = 4.7;
      bool draw_missing_residue_loops_flag = true;
      bonds.do_Ca_plus_ligands_bonds(atom_sel_ligand, imol_no, geom, min_dist, max_dist, draw_hydrogen_atoms_flag,
                                     draw_missing_residue_loops_flag);
      bonds_box.clear_up();
      bonds_box = bonds.make_graphical_bonds_no_thinning();
      std::vector<glm::vec4> colour_table = make_colour_table(against_a_dark_background);
      make_instanced_graphical_bonds_bonds(m, bonds_box, bond_radius, n_slices, n_stacks, colour_table);
      bonds_box.clear_up();
   }

   if (mode == "VDW-BALLS") {

      Bond_lines_container bonds(geom, no_bonds_to_these_atoms, draw_hydrogen_atoms_flag);
      bonds.do_colour_by_chain_bonds(atom_sel_ligand, false, imol_no, draw_hydrogen_atoms_flag,
                                     draw_missing_residue_loops, change_c_only_flag, goodsell_mode, do_rota_markup);
      bonds_box.clear_up();
      bonds_box = bonds.make_graphical_bonds();

      std::vector<glm::vec4> colour_table = make_colour_table(against_a_dark_background); // uses bonds_box
      make_graphical_bonds_spherical_atoms_with_vdw_radii_instanced(m, bonds_box, num_subdivisions, colour_table, *geom);
      bonds_box.clear_up();
   }

   return m;
}


coot::instanced_mesh_t
coot::molecule_t::get_extra_restraints_mesh(int mode) const {

   auto  distortion_score_harmonic = [] (const double &bl, const double &target_dist,
                                   const double &sigma) {
      double bit = bl - target_dist;
      double z = bit/sigma;
      double distortion = z*z;
      return distortion;
   };

   auto  distortion_score_GM = [] (const double &bl, const double &target_dist,
                                   const double &sigma, const double &alpha) {
      double bit = bl - target_dist;
      double z = bit/sigma;
      double distortion = z*z/(1+alpha*z*z);
      return distortion;
   };

   auto convert_vertices = [] (const std::vector<coot::api::vnc_vertex> &v_in) {
      std::vector<coot::api::vn_vertex> v_out(v_in.size());
      for (unsigned int i=0; i<v_in.size(); i++) {
         const auto &v = v_in[i];
         v_out[i] = coot::api::vn_vertex(v.pos, v.normal);
      }
      return v_out;
   };

   auto clipper_to_glm = [] (const clipper::Coord_orth &co) {
                            return glm::vec3(co.x(), co.y(), co.z());
                         };

   // type is 1 for harmonic and 2 for GM.
   auto make_instancing_data = [this, distortion_score_GM, distortion_score_harmonic, clipper_to_glm] (int atom_index_udd_handle, mmdb::Atom *at_1, mmdb::Atom *at_2,
                                                 double target_bond_dist, double target_sigma, double geman_mcclure_alpha,
                                                 int type) {

                  glm::vec3 z0(0,0,0);
                  glm::vec3 z1(0,0,1);
                  glm::vec4 col_base(0.5, 0.5, 0.5, 1.0);
                  clipper::Coord_orth p_1 = coot::co(at_1);
                  clipper::Coord_orth p_2 = coot::co(at_2);
                  clipper::Coord_orth delta = p_2 - p_1;
                  clipper::Coord_orth delta_uv = clipper::Coord_orth(delta.unit());
                  double bl = std::sqrt(delta.lengthsq());
                  glm::vec3 delta_uv_glm = clipper_to_glm(delta_uv);
                  glm::mat4 ori = glm::orientation(delta_uv_glm, z1);
                  glm::vec3 p = clipper_to_glm(p_2);
                  double delta_length = target_bond_dist - bl;
                  if (delta_length >  1.0) delta_length =  1.0;
                  if (delta_length < -1.0) delta_length = -1.0;

                  double penalty = 0;
                  if (type == 1) penalty = distortion_score_harmonic(bl, target_bond_dist, target_sigma);
                  if (type == 2) penalty = distortion_score_GM(bl, target_bond_dist, target_sigma, geman_mcclure_alpha);
                  double width = 0.23 * penalty;
                  if (width < 0.01) width = 0.01;
                  if (width > 0.08) width = 0.08;

		  width *= 2.5;

                  glm::vec3 s(width, width, bl); // a function of delta_length?
                  glm::vec4 col = col_base + delta_length * glm::vec4(-0.8f, 0.8f, -0.8, 0.0f);
                  col = 0.8 * col; // calm down
		  col.a = 1.0;
                  return instancing_data_type_B_t(p, col, s, ori);
   };


   double geman_mcclure_alpha = 1.0;

   instanced_mesh_t im;
   if (mode > -999) { // check the mode here

      coot::instanced_geometry_t igeom;
      unsigned int n_slices = 6;
      unsigned int n_stacks = 6;
      glm::vec3 z0(0,0,0);
      glm::vec3 z1(0,0,1);
      std::pair<glm::vec3, glm::vec3> pp(z0, z1);
      cylinder c_00(pp, 1.0, 1.0, 1.0, n_slices, n_stacks);
      c_00.add_flat_end_cap();
      c_00.add_flat_start_cap();
      igeom.vertices = convert_vertices(c_00.vertices);
      igeom.triangles = c_00.triangles;
      int atom_index_udd_handle = atom_sel.UDDAtomIndexHandle;

      if (! extra_restraints.bond_restraints.empty()) {
         for (unsigned int i=0; i<extra_restraints.bond_restraints.size(); i++) {
            const auto &r = extra_restraints.bond_restraints[i];
            mmdb::Atom *at_1 = get_atom(r.atom_1);
            mmdb::Atom *at_2 = get_atom(r.atom_2);
            if (at_1) {
               if (at_2) {
                  int idx_1;
                  int idx_2;
                  at_1->GetUDData(atom_index_udd_handle, idx_1);
                  at_2->GetUDData(atom_index_udd_handle, idx_2);
                  if (no_bonds_to_these_atom_indices.find(idx_1) != no_bonds_to_these_atom_indices.end()) continue;
                  if (no_bonds_to_these_atom_indices.find(idx_2) != no_bonds_to_these_atom_indices.end()) continue;
                  auto idB = make_instancing_data(atom_index_udd_handle, at_1, at_2, r.bond_dist, r.esd, geman_mcclure_alpha, 1);
                  igeom.instancing_data_B.push_back(idB);
               } else {
                  std::cout << "WARNING:: no atom found " << r.atom_2 << std::endl;
               }
            } else {
               std::cout << "WARNING:: no atom found " << r.atom_1 << std::endl;
            }
         }
         im.add(igeom);
      }

      if (! extra_restraints.geman_mcclure_restraints.empty()) {
         for (unsigned int i=0; i<extra_restraints.geman_mcclure_restraints.size(); i++) {
            const auto &r = extra_restraints.geman_mcclure_restraints[i];
            mmdb::Atom *at_1 = get_atom(r.atom_1);
            mmdb::Atom *at_2 = get_atom(r.atom_2);
            if (at_1) {
               if (at_2) {
                  int idx_1;
                  int idx_2;
                  at_1->GetUDData(atom_index_udd_handle, idx_1);
                  at_2->GetUDData(atom_index_udd_handle, idx_2);
                  if (no_bonds_to_these_atom_indices.find(idx_1) != no_bonds_to_these_atom_indices.end()) continue;
                  if (no_bonds_to_these_atom_indices.find(idx_2) != no_bonds_to_these_atom_indices.end()) continue;
                  auto idB = make_instancing_data(atom_index_udd_handle, at_1, at_2, r.bond_dist, r.esd, geman_mcclure_alpha, 2);
                  igeom.instancing_data_B.push_back(idB);
               } else {
                  std::cout << "WARNING:: no atom found " << r.atom_2 << std::endl;
               }
            } else {
               std::cout << "WARNING:: no atom found " << r.atom_1 << std::endl;
            }
         }
         im.add(igeom);
      }
   }
   return im;
}
