#ifndef MODEL_MOLECULE_MESHES_HH
#define MODEL_MOLECULE_MESHES_HH

#include "coords/graphical-bonds-container.hh"
#include "api/instancing.hh"
#include "Shader.hh"
#include "Material.hh"
#include "Mesh.hh"
#include "api/bond-colour.hh"

// This class draws the meshes in instanced_mesh_t
// (and that is a vector of instanced meshes and a simple_mesh_t)

class model_molecule_meshes_t {

public:
   model_molecule_meshes_t() : simple_mesh(Mesh("model_molecule_meshes_t constructor")) {}
   coot::instanced_mesh_t im;
   std::vector<Mesh> instanced_meshes;
   Mesh simple_mesh;
   Material material;
   std::string name;

   void import(const coot::instanced_mesh_t &api_mol_mesh);
   void set_name(const std::string &n) { name = n; }

   // wrapper for both the instanced and simple meshes
   void draw(Shader *shader_mesh,
             Shader *shader_instanced_mesh,
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
                    const glm::mat4 &mvp,
                    const glm::mat4 &view_rotation_matrix,
                    const std::map<unsigned int, lights_info_t> &lights,
                    const glm::vec3 &eye_position, // eye position in view space (not molecule space)
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
   void draw_simple_bond_lines(Shader *shader,
                               const glm::mat4 &glm,
                               const glm::vec4 &background_colour,
                               float line_width,
                               bool do_depth_fog);

   // wrap draw_simple_for_ssao() and draw_instances_for_ssao()
   void draw_for_ssao(Shader *shader_for_meshes_p,
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

   // the simple-lines option for the main molecule
   void make_bond_lines(const graphical_bonds_container &bonds_box, const std::vector<glm::vec4> &colour_table);

   // 20230828-PE it seem sthat udd_handle_bonded_type is not used at the moment
   void make_graphical_bonds(const graphical_bonds_container &bonds_box,
                             bool draw_cis_peptide_markups,
                             float atom_radius, float bond_radius,
                             int num_subdivisions, int n_slices, int n_stacks,
                             const std::vector<glm::vec4> &colour_table);

};

#endif // MODEL_MOLECULE_MESHES_HH
