//
#ifndef MESH_HH
#define MESH_HH

#include <chrono>

#include "generic-vertex.hh"
#include "g_triangle.hh"

#include "lights-info.hh"
#include "Shader.hh"
#include "Material.hh"
#include "Particle.hh"
#include "molecular-triangles-mesh.hh"
#include "simple-distance-object.hh"

#ifdef USE_ASSIMP
#include <assimp/scene.h>
#endif

#ifdef THIS_IS_HMT
#else
#include "coords/graphical-bonds-container.hh"
#include "coords/mmdb-crystal.h"
#endif

class Mesh {
   enum { VAO_NOT_SET = 99999999 };
   void setup_debugging_instancing_buffers(); // or buffers, when we add rotation
   Material material;
   void setup_instanced_balls( Shader *shader_p, const Material &material_in);
   void setup_instanced_dodecs(Shader *shader_p, const Material &material_in);
   void setup_buffers();
   // rts rotation, translation & scale
   void setup_matrix_and_colour_instancing_buffers(const std::vector<glm::mat4> &mats, const std::vector<glm::vec4> &colours);
   // this is when we learning/testing/playing - I don't want to delete it (yet) just in case I want
   // to go back and test out previous example code.
   void setup_matrix_and_colour_instancing_buffers_standard(const std::vector<glm::mat4> &mats, const std::vector<glm::vec4> &colours);
   int n_instances; // instances to be drawn
   int n_instances_allocated; // that we made space for in glBufferData()
   bool normals_are_setup;
   GLuint normals_vao;
   GLuint normals_buffer_id;
   GLuint normals_colour_buffer_id;
   bool first_time;
   float hydrogen_bond_cylinders_angle;
   unsigned int particle_draw_count;
   void init();
   bool draw_this_mesh;
#if USE_ASSIMP
   aiScene generate_scene() const;
#endif

public:
   GLuint vao;
   GLuint buffer_id;
   GLuint index_buffer_id;
   GLuint inst_rts_buffer_id;
   GLuint inst_model_translation_buffer_id;
   GLuint inst_colour_buffer_id;
   GLuint inst_scales_buffer_id;
   unsigned int n_symmetry_atom_lines_vertices;
   bool this_mesh_is_closed;
   bool is_instanced;
   bool is_instanced_colours;
   bool is_instanced_with_rts_matrix;
   bool use_blending;
   std::vector<s_generic_vertex> vertices;
   std::vector<g_triangle> triangles;
   Shader shader_for_draw_normals;
   std::string name;
   unsigned int type; // from molecular triangles object type
   std::chrono::time_point<std::chrono::system_clock>  time_constructed;

   Mesh() { init(); }
   // import from somewhere else
   explicit Mesh(const std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > &indexed_vertices);
   explicit Mesh(const std::string &name_in) : name(name_in) { init(); }
   explicit Mesh(const molecular_triangles_mesh_t &mtm);
   // If this mesh will become part of another mesh, then we don't wan to setup buffers for this one
   void load_from_glTF(const std::string &file_name, bool include_call_to_setup_buffers=true);
   void export_to_glTF(const std::string &file_name, bool use_binary_format); // 20210927-PE I'd rather not include tiny_gltf.h in Mesh.hh

   void debug() const;
   void debug_to_file() const;
   void clear() {
      // delete some gl buffers here
      is_instanced = false;
      is_instanced_colours = false;
      is_instanced_with_rts_matrix = false;
      vertices.clear();
      triangles.clear();
      use_blending = false;
      normals_are_setup = false;
   }
   void close();
   void set_draw_mesh_state(bool state) { draw_this_mesh = state; }
   void set_name(const std::string &n) { name = n; }
   void import(const std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > &indexed_vertices); // adds to the mesh
   void import(const std::vector<s_generic_vertex> &gv, const std::vector<g_triangle> &indexed_vertices);   // adds to the mesh
   void import(const std::vector<position_normal_vertex> &verts,
               const std::vector<g_triangle> &indexed_vertices,
               const glm::vec4 &colour);
   // void setup(Shader *shader_p, const Material &material_in);  I don't need the shader, do I?
   void setup(const Material &material_in);
   // can be considered as "draw_self()"
   void draw(Shader *shader,
             const glm::mat4 &mvp,
             const glm::mat4 &view_rotation_matrix,
             const std::map<unsigned int, lights_info_t> &lights,
             const glm::vec3 &eye_position, // eye position in view space (not molecule space)
             const glm::vec4 &background_colour,
             bool do_depth_fog);
   void draw_particles(Shader *shader_p, const glm::mat4 &mvp, const glm::mat4 &view_rotation);
   void draw_normals(const glm::mat4 &mvp, float normal_scaling); // debugging
   void draw_symmetry(Shader *shader_p,
                      const glm::mat4 &mvp,
                      const glm::mat4 &mouse_based_rotation_matrix,
                      const std::map<unsigned int, lights_info_t> &lights,
                      const glm::vec3 &eye_position,
                      const glm::vec4 &background_colour,
                      bool do_depth_fog);


   // the current draw method is a mess with changing between colours and instancing adn idex of the attribute array
   // the new system has this layout:
   // 0 vertex position
   // 1 vertex normal
   // 2 vertex colour
   // 3 instance colour
   // 4 instance rot-trans-1  // these include the scale
   // 5 instance rot-trans-2
   // 6 instance rot-trans-3
   // 7 instance rot-trans-4
   void draw_instanced(Shader *shader,
                       const glm::mat4 &mvp,
                       const glm::mat4 &view_rotation_matrix,
                       const std::map<unsigned int, lights_info_t> &lights,
                       const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                       const glm::vec4 &background_colour,
                       bool do_depth_fog);

   // testing/example functions
   void setup_rama_balls(Shader *shader_p, const Material &material_in); // call fill_rama_balls() and setup_buffers()
   void setup_simple_triangles(Shader *shader_p, const Material &material_in);
   // setup_simple_triangles() calls fill_with_simple_triangles_vertices()
   void fill_with_simple_triangles_vertices();
   void fill_with_direction_triangles();
   void setup_instanced_debugging_objects(Shader *shader_p, const Material &material_in);

   // this is an import function, name it so.
   void import_and_setup_instanced_cylinders(Shader *shader_p,
                                             const Material &material_in,
                                             const std::vector<glm::mat4> &mats,
                                             const std::vector<glm::vec4> &colours);
   // this is an import function, name it so.
   void setup_instanced_octahemispheres(Shader *shader_p,
                                        const Material &material_in,
                                        const std::vector<glm::mat4> &mats,
                                        const std::vector<glm::vec4> &colours);

   void test_cyclinders(Shader *shader_p, const Material &material_in);
   void setup_camera_facing_outline();
   void setup_camera_facing_quad();
   void setup_camera_facing_hex();
   void setup_camera_facing_polygon(unsigned int n_sides = 8, float scale=1.0);
   void setup_hydrogen_bond_cyclinders(Shader *shader_p, const Material &material_in);

   void setup_rtsc_instancing(Shader *shader_p,
                              const std::vector<glm::mat4> &mats,
                              const std::vector<glm::vec4> &colours,
                              unsigned int n_instances_in,
                              const Material &material_in);

   // These are separated because they have different vertices
   // They use a sphere, hemisphere or cylinder at the origin a the vertices and each bond and atom
   // has a translation and rotation (so that the geometry of the bonds matches the geom of the
   // hemisphere (octaballs)).
   void make_graphical_bonds_spherical_atoms(Shader *shader_p,
                                             const Material &material,
                                             const graphical_bonds_container &gbc, int udd_handle_bonded_type,
                                             float atom_radius,
                                             float bond_radius,
                                             unsigned int num_subdivisions,
                                             glm::vec4 (*get_glm_colour_for_bonds) (int, int));
   void make_graphical_bonds_hemispherical_atoms(Shader *shader_p,
                                                 const Material &material,
                                                 const graphical_bonds_container &gbc, int udd_handle_bonded_type,
                                                 float atom_radius,
                                                 float bond_radius,
                                                 unsigned int num_subdivisions,
                                                 glm::vec4 (*get_glm_colour_for_bonds) (int, int));
   void make_graphical_bonds_bonds(Shader *shader_p,
                                   const Material &material,
                                   const graphical_bonds_container &gbc,
                                   float bond_radius,
                                   unsigned int n_slices,
                                   unsigned int n_stacks,
                                   glm::vec4 (*get_glm_colour_for_bonds) (int, int));
   void make_symmetry_atoms_bond_lines(Shader *shader_p,
                                  const std::vector<std::pair<graphical_bonds_container, std::pair<symm_trans_t, Cell_Translation> > > &symmetry_bonds_boxes);

   // update vertices
   // the vertices have been updated (externally) (say by position and colour) so they need to
   // be re-pushed to the graphics card
   //
   void update_vertices(); // push to graphics.

   // when the position, orientation and colour change:
   // transfer to the graphics_card with glBufferSubData
   void update_instancing_buffer_data(const std::vector<glm::mat4> &mats,
                                      const std::vector<glm::vec4> &colours);
   // when the positions change: 20210826-PE Hmm. Thisis not clear
   void update_instancing_buffer_data_standard(const std::vector<glm::mat4> &mats);

   // void setup_instancing_buffers(const particle_container_t &particles);
   void setup_vertex_and_instancing_buffers_for_particles(unsigned int n_particles); // setup the buffer, don't add data
   void update_instancing_buffer_data_for_particles(const particle_container_t &particles);

   void set_draw_this_mesh(bool state); // only set true if there are vertices and triangles
   bool get_draw_this_mesh() const { return draw_this_mesh; };

   // make space
   void setup_instancing_buffer_data(Shader *shader_p,
                                     const Material &mat,
                                     const std::vector<glm::mat4> &instanced_matrices,
                                     const std::vector<glm::vec4> &instanced_colours);

   void fill_rama_balls(); // make up some balls
   void add_one_ball(float scale, const glm::vec3 &centre);
   void add_one_origin_ball();
   void add_one_origin_cylinder(unsigned int nslices=8, unsigned int n_stacks=2);
   void add_one_origin_octasphere(unsigned int num_subdivisions=1); // default matches 8 sided cylinder
   void add_one_origin_octahemisphere(unsigned int num_subdivisions=1); // default matches 8 sided cylinder
   void add_one_origin_dodec();
   void fill_one_dodec();
   void flatten_triangles(); // needs implementation (will generate new vertices).
   void smooth_triangles();  // needs implementation.
   bool is_closed() const { return this_mesh_is_closed; }
   bool have_instances() const { return (n_instances > 0); }
   void translate_by(const glm::vec3 &t);
   bool export_as_obj_via_assimp(const std::string &file_name) const;
   bool export_as_obj_internal(const std::string &file_name) const;
   bool export_as_obj(const std::string &file_name) const;
   bool export_as_obj(std::ofstream &f, unsigned int vertex_index_offset) const;

   // make the matrix (called several times). After which, the calling function calls update_instancing_buffer_data()
   static glm::mat4 make_hydrogen_bond_cylinder_orientation(const glm::vec3 &p1, const glm::vec3 &p2, float theta);

   void add_dashed_line(const coot::simple_distance_object_t &l, const Material &material, const glm::vec4 &colour);
   void add_dashed_angle_markup(const glm::vec3 &pos_1, const glm::vec3 &pos_2, const glm::vec3 &pos_3,
                                const glm::vec4 &colour_in, const Material &mat);

   void apply_scale(float scale_factor);  // scale the positions in the vertices
   void apply_transformation(const glm::mat4 &m);  // transform the positions in the vertices

};

#endif
