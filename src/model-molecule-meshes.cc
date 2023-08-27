
#include "model-molecule-meshes.hh"
#include "api/make-instanced-graphical-bonds.hh" // make_instanced_graphical_bonds_spherical_atoms() etc
#include "generic-vertex.hh"
void
model_molecule_meshes_t::draw_simple_bond_lines(Shader *shader,
                                                const glm::mat4 &glm,
                                                const glm::vec4 &background_colour,
                                                float line_width,
                                                bool do_depth_fog) {

}


// wrap draw_simple_for_ssao() and draw_instances_for_ssao()
void
model_molecule_meshes_t::draw_for_ssao(Shader *shader_for_meshes_p,
                                       const glm::mat4 &model,
                                       const glm::mat4 &view,
                                       const glm::mat4 &projection) { // draw into the gbuffer framebuffer.

}


void
model_molecule_meshes_t::make_graphical_bonds(const graphical_bonds_container &bonds_box,
                                              coot::api_bond_colour_t bonds_box_type,
                                              const std::string &model_representation_mode,
                                              int udd_handle_bonded_type,
                                              bool draw_cis_peptide_markups,
                                              float atom_radius, float bond_radius,
                                              int num_subdivisions, int n_slices, int n_stacks,
                                              const std::vector<glm::vec4> &colour_table,
                                              coot::protein_geometry *prtein_geometry_p) {

   std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ model_molecule_meshes_t::make_graphical_bonds() @@@@@@@@@@@@@@@@@@@@@@@@@@@"
             << std::endl;

   if (true) { // test the style
      make_instanced_graphical_bonds_spherical_atoms(im, bonds_box, bonds_box_type, udd_handle_bonded_type, atom_radius, bond_radius,
                                                     num_subdivisions, colour_table);
      make_instanced_graphical_bonds_bonds(im, bonds_box, bond_radius, n_slices, n_stacks, colour_table);
      make_graphical_bonds_cis_peptides(im.markup, bonds_box);
   }

   // now convert im meshes to src style "Mesh"es.
   std::cout << "in make_graphical_bonds() found im.geom.size()" << im.geom.size() << std::endl;
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
      m.setup_buffers();

      if (! ig.instancing_data_A.empty()) {
         // instancing data
         glm::mat4 unit_matrix(1.0f);
         std::vector<glm::mat4> matrices(ig.instancing_data_A.size(), unit_matrix);
         std::vector<glm::vec4> colours(ig.instancing_data_A.size());
         for (unsigned int i_A=0; i_A<ig.instancing_data_A.size(); i_A++) {
	    float sar = atom_radius;
            const coot::instancing_data_type_A_t &Atd = ig.instancing_data_A[i_A];
            colours[i_A] = Atd.colour;
	    glm::mat4 mm = glm::scale(unit_matrix, Atd.size);
            matrices[i_A] = glm::translate(mm, Atd.position/sar);
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

   std::cout << "2222222222222222222222 im.markup.vertices size " << im.markup.vertices.size() << std::endl;
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


// the simple-lines option for the main molecule
void
model_molecule_meshes_t::make_bond_lines(const graphical_bonds_container &bonds_box, const std::vector<glm::vec4> &colour_table) {

}


// wrapper for both the instanced and simple meshes
void
model_molecule_meshes_t::draw(Shader *shader_mesh,
                              Shader *shader_for_instanced_mesh,
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
   draw_instances(shader_for_instanced_mesh, mvp, view_rotation_matrix, lights, eye_position, background_colour, do_depth_fog, transferred_colour_is_instanced);
   draw_simple(shader_mesh, mvp, view_rotation_matrix, lights, eye_position, opacity, background_colour, gl_lines_mode, do_depth_fog, show_just_shadows);

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

   std::cout << "model_molecule_meshes_t::draw_instances() " << instanced_meshes.size() << std::endl;
   for (unsigned int i=0; i<instanced_meshes.size(); i++) {
      auto &mesh = instanced_meshes[i];
      std::cout << "   calling mesh.draw_instanced() with shader " << shader_for_instanced_meshes_p->name << std::endl;
      mesh.draw_instanced(shader_for_instanced_meshes_p, mvp, view_rotation_matrix, lights, eye_position, background_colour,
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
                                     float opacity,
                                     const glm::vec4 &background_colour,
                                     bool gl_lines_mode, // i.e. as chickenwire
                                     bool do_depth_fog,
                                     bool show_just_shadows) {

   std::cout << "model_molecule_meshes_t::draw_simple()" << std::endl;
   simple_mesh.draw(shader, mvp, view_rotation_matrix, lights, eye_position, opacity, background_colour,
                    gl_lines_mode, do_depth_fog, show_just_shadows);
}
