/*
 * src/model-molecule-meshes.hh
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
#ifndef MODEL_MOLECULE_MESHES_HH
#define MODEL_MOLECULE_MESHES_HH

#include "coords/graphical-bonds-container.hh"
#include "api/instancing.hh"
#include "Shader.hh"
#include "Material.hh"
#include "Mesh.hh"
#include "stereo-eye.hh"

// This class draws the meshes in instanced_mesh_t
// (and that is a vector of instanced meshes and a simple_mesh_t)

class model_molecule_meshes_t {

   void convert_and_fill_meshes(const coot::instanced_mesh_t &im);

public:
   model_molecule_meshes_t() : simple_mesh(Mesh("model_molecule_meshes_t constructor")) {}
   coot::instanced_mesh_t im;
   std::vector<Mesh> instanced_meshes;
   Mesh simple_mesh;
   Material material;
   std::string name;
   void set_debug_mode(bool state); // transfer this to the meshes

   void import(const coot::instanced_mesh_t &api_mol_mesh);
   void set_name(const std::string &n) { name = n; }

   bool empty() const;

   // wrapper for both the instanced and simple meshes
   void draw(Shader *shader_mesh,
             Shader *shader_instanced_mesh,
             stereo_eye_t eye,
             const glm::mat4 &mvp,
             const glm::mat4 &view_rotation_matrix,
             const std::map<unsigned int, lights_info_t> &lights,
             const glm::vec3 &eye_position, // eye position in view space (not molecule space)
             float opacity,
             const glm::vec4 &background_colour,
             bool gl_lines_mode, // i.e. as chickenwire
             bool do_depth_fog,
             bool show_just_shadows);

   void draw_instances(Shader *shader_for_instanced_meshes_p,
                       stereo_eye_t eye,
                       const glm::mat4 &mvp,
                       const glm::mat4 &view_rotation_matrix,
                       const std::map<unsigned int, lights_info_t> &lights,
                       const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                       const glm::vec4 &background_colour,
                       bool do_depth_fog,
                       bool transferred_colour_is_instanced,
                       bool do_pulse = false,
                       bool do_rotate_z = false,
                       float pulsing_amplitude = 0.0f,
                       float pulsing_frequency = 0.f,
                       float pulsing_phase_distribution = 0.0f,
                       float z_rotation_angle = 0.0f);

   // the (previous) Mesh draw function
   void draw_simple(Shader *shader,
                    stereo_eye_t eye,
                    const glm::mat4 &mvp,
                    const glm::mat4 &view_rotation_matrix,
                    const std::map<unsigned int, lights_info_t> &lights,
                    const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                    const glm::vec3 &rotation_centre, // maybe not needed (was for fresnel)
                    float opacity,
                    const glm::vec4 &background_colour,
                    bool gl_lines_mode, // i.e. as chickenwire
                    bool do_depth_fog,
                    bool show_just_shadows);

   // this class takes over the rendering of these meshes.
   //
   // Mesh molecule_as_mesh; // non-instancing
   //
   // Mesh molecule_as_mesh_atoms_1; // for instancing (sphere)
   // Mesh molecule_as_mesh_atoms_2; // for instancing (hemisphere - not currently used)
   // Mesh molecule_as_mesh_bonds_c00; // for instancing, no-end-caps
   // Mesh molecule_as_mesh_bonds_c10; // for instancing, start end cap
   // Mesh molecule_as_mesh_bonds_c01; // for instancing, end end-cap
   // Mesh molecule_as_mesh_bonds_c11; // for instancing, both end-caps
   // Mesh molecule_as_mesh_bonds_round_cap_start // instancing

   void set_material(const Material &material_in) { material = material_in; }
   void set_material_specularity(float specular_strength, float shininess) {
      material.specular_strength = specular_strength;
      material.shininess = shininess; }
   void set_material_diffuse(const glm::vec4 &diffuse) { material.diffuse = diffuse; }
   void set_material_ambient(const glm::vec4 &ambient) { material.ambient = ambient; }

   // --------------------- drawing function ------------------------------------

   void draw_simple_bond_lines(Shader *shader,
                               const glm::mat4 &glm,
                               const glm::vec4 &background_colour,
                               float line_width,
                               bool do_depth_fog);

   // wrap draw_simple_for_ssao() and draw_instances_for_ssao()
   void draw_for_ssao(Shader *shader_for_meshes_p,
                      Shader *shader_for_instanced_meshes_p,
                      const glm::mat4 &model,
                      const glm::mat4 &view,
                      const glm::mat4 &projection);  // draw into the gbuffer framebuffer.

   // basic means "non-instanced" here
   void draw_simple_for_ssao(Shader *shader_for_meshes_p,
                             const glm::mat4 &model,
                             const glm::mat4 &view,
                             const glm::mat4 &projection);  // draw into the gbuffer framebuffer.

   void draw_instances_for_ssao(Shader *shader_p,
                                const glm::mat4 &model,
                                const glm::mat4 &view,
                                const glm::mat4 &projection);

   // instanced models
   void draw_molecule_with_shadows(Shader *shader,
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
                                   bool show_just_shadows);


   // --------------------------- mesh-making functions -----------------------------------

   // the simple-lines option for the main molecule
   void make_bond_lines(const graphical_bonds_container &bonds_box, const std::vector<glm::vec4> &colour_table);

   // 20230828-PE it seem sthat udd_handle_bonded_type is not used at the moment
   void make_graphical_bonds(int imol,
                             const graphical_bonds_container &bonds_box,
                             float atom_radius, float bond_radius,
                             bool show_atoms_as_aniso_flag,
                             float aniso_probability,
                             bool show_aniso_atoms_as_ortep_flag,
                             int num_subdivisions, int n_slices, int n_stacks,
                             const std::vector<glm::vec4> &colour_table);

   void make_symmetry_bonds(int imol,
                            const std::vector<std::pair<graphical_bonds_container, std::pair<symm_trans_t, Cell_Translation> > > &symmetry_bonds_box,
                            float atom_radius, float bond_radius,
                            int num_subdivisions, int n_slices, int n_stacks,
                            const std::vector<glm::vec4> &colour_table);

   void add_rotamer_dodecs(int imol, const graphical_bonds_container &bonds_box);
   void add_ramachandran_spheres(int imol, const graphical_bonds_container &gbc);

   std::pair<bool, coot::Cartesian> get_HA_unit_vector(mmdb::Residue *r) const;

   mmdb::Residue *get_residue(int imol, const coot::residue_spec_t &spec) const;

   std::pair<std::vector<coot::api::vn_vertex>, std::vector<g_triangle> > get_dodec_vertices_and_triangles() const;

};

#endif // MODEL_MOLECULE_MESHES_HH
