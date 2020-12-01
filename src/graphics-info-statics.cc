
#ifdef USE_PYTHON
#include "Python.h"
#endif
#include "graphics-info.h"

bool graphics_info_t::residue_type_selection_was_user_picked_residue_range = false;

bool graphics_info_t::make_auto_h_bond_restraints_flag = false;

bool graphics_info_t::do_rotamer_restraints = false;

bool graphics_info_t::do_debug_refinement = false;

std::atomic<bool> graphics_info_t::on_going_updating_map_lock(false);

std::map<keyboard_key_t, key_bindings_t> graphics_info_t::key_bindings_map;

std::string graphics_info_t::mtz_file_for_refmac;

bool graphics_info_t::convert_dictionary_planes_to_improper_dihedrals_flag = false;

GtkWidget *graphics_info_t::main_window = NULL;

// now the clipping planes are scale, not offsets
float graphics_info_t::clipping_front = 1.0;
float graphics_info_t::clipping_back  = 1.0;

#ifdef USE_MOLECULES_TO_TRIANGLES
std::shared_ptr<Renderer>   graphics_info_t::mol_tri_renderer    = 0;
std::shared_ptr<SceneSetup> graphics_info_t::mol_tri_scene_setup = 0;
#endif // USE_MOLECULES_TO_TRIANGLES

// --------------------------------------------------------------------------------------------

float *graphics_info_t::mvp = new float[16];
int    graphics_info_t::mvp_location = -1;
int    graphics_info_t::view_rotation_location = -1;
glm::quat graphics_info_t::glm_quat = glm::quat(1,0,0,0);
GLuint graphics_info_t::programID_for_central_cube = 0;
GLuint graphics_info_t::central_cube_vertexarray_id = 0;
GLuint graphics_info_t::central_cube_array_buffer_id = 0;
GLuint graphics_info_t::central_cube_index_buffer_id = 0;
GLuint graphics_info_t::hud_text_vertexarray_id = 0;
GLuint graphics_info_t::hud_text_array_buffer_id = 0;
GLuint graphics_info_t::screen_quad_vertex_array_id = 0;
GLuint graphics_info_t::blur_quad_vertex_array_id = 0;
GLuint graphics_info_t::textureColorbuffer_screen = 0;
GLuint graphics_info_t::textureColorbuffer_blur = 0;
framebuffer graphics_info_t::screen_framebuffer;
framebuffer graphics_info_t::blur_framebuffer;
bool graphics_info_t::perspective_projection_flag = false;


// --------------------------------------------------------------------------------------------

// GLuint graphics_info_t::programID_for_maps = 0; in a shader now  - as
//programID_for_central_cube should be
Shader graphics_info_t::shader_for_maps;
Shader graphics_info_t::shader_for_map_caps;
Shader graphics_info_t::shader_for_models;
Shader graphics_info_t::shader_for_model_as_meshes;
Shader graphics_info_t::shader_for_moleculestotriangles;
Shader graphics_info_t::shader_for_central_cube;
Shader graphics_info_t::shader_for_origin_cube;
Shader graphics_info_t::shader_for_hud_text;
Shader graphics_info_t::shader_for_hud_geometry_bars;
Shader graphics_info_t::shader_for_hud_geometry_labels;
Shader graphics_info_t::shader_for_atom_labels;
Shader graphics_info_t::shader_for_rama_balls;
Shader graphics_info_t::shader_for_screen;
Shader graphics_info_t::shader_for_blur;
Shader graphics_info_t::shader_for_lines;
Shader graphics_info_t::shader_for_lines_pulse;
Shader graphics_info_t::shader_for_particles;
Shader graphics_info_t::shader_for_instanced_objects; // used for boids - also HOLE
Shader graphics_info_t::shader_for_hud_geometry_tooltip_text;
Shader graphics_info_t::shader_for_happy_face_residue_markers;
meshed_generic_display_object graphics_info_t::mesh_for_environment_distances;
std::chrono::time_point<std::chrono::system_clock> graphics_info_t::previous_frame_time = std::chrono::high_resolution_clock::now();
std::chrono::time_point<std::chrono::system_clock> graphics_info_t::previous_frame_time_for_per_second_counter = std::chrono::high_resolution_clock::now();
long graphics_info_t::frame_counter = 0;
long graphics_info_t::frame_counter_at_last_display = 0;
std::queue<std::chrono::time_point<std::chrono::system_clock> > graphics_info_t::frame_draw_queue;
std::set<mmdb::Residue *> graphics_info_t::moving_atoms_visited_residues;
mmdb::Atom *graphics_info_t::active_atom_for_hud_geometry_bar = 0;

unsigned int graphics_info_t::draw_count_for_happy_face_residue_markers = 0;

glm::vec3 graphics_info_t::eye_position = glm::vec3(0,0,95);
float graphics_info_t::screen_z_near_perspective =  76.0; // was 83
float graphics_info_t::screen_z_far_perspective  = 125.0;

std::map<unsigned int, lights_info_t> graphics_info_t::lights;

molecule_class_info_t graphics_info_t::moving_atoms_molecule;

bool graphics_info_t::vera_font_loaded = false;

unsigned int graphics_info_t::n_atom_pulls = 0;
GLuint graphics_info_t::m_VertexArray_for_pull_restraints_ID = 0;
GLuint graphics_info_t::m_VertexBuffer_for_pull_restraints_ID = 0;
GLuint graphics_info_t::m_IndexBuffer_for_atom_pull_restraints_ID = 0;
unsigned int graphics_info_t::n_triangles_for_atom_pull_restraints = 0;
unsigned int graphics_info_t::n_vertices_for_atom_pull_restraints = 0;

coot::view_info_t graphics_info_t::reorienting_residue_start_view;
coot::view_info_t graphics_info_t::reorienting_residue_end_view;

bool graphics_info_t::smooth_scroll_on_going = false;

bool graphics_info_t::shader_do_ambient_occlusion_flag = true;
bool graphics_info_t::shader_do_depth_blur_flag = true;
bool graphics_info_t::shader_do_depth_fog_flag = true;
bool graphics_info_t::shader_do_outline_flag = false;
bool graphics_info_t::draw_normals_flag = false;

// static
void
graphics_info_t::make_gl_context_current(bool gl_context_current_request_index) {

   // what does this do now?

#if 0
   if (glareas.empty()) return;
   if (display_mode_use_secondary_p()) {
      if (gl_context_current_request_index == GL_CONTEXT_SECONDARY) {
         if (glareas.size() > 1) {
            GtkWidget *glarea = glareas[1];
            if (glarea) {
               make_current_gl_context(glarea);
            }
         }
      }
      if (gl_context_current_request_index == GL_CONTEXT_MAIN) {
         GtkWidget *glarea = glareas[0];
	 if (glarea) {
            make_current_gl_context(glarea);
	 }
      }
   } else {
      if (gl_context_current_request_index == GL_CONTEXT_MAIN) {
         GtkWidget *glarea = glareas[0];
	 if (glarea) {
            make_current_gl_context(glarea);
	 }
      }
   }
#endif
}

float graphics_info_t::contact_dots_density = 1.0;

bool graphics_info_t::draw_missing_loops_flag = true;

bool graphics_info_t::sequence_view_is_docked_flag = true;

unsigned int graphics_info_t::framebuffer_scale = 1; // on supersampling by default.


bool graphics_info_t::do_tick_particles = false;
bool graphics_info_t::do_tick_spin = false;
bool graphics_info_t::do_tick_boids = false;
bool graphics_info_t::do_tick_hydrogen_bonds_mesh = false;
bool graphics_info_t::do_tick_happy_face_residue_markers = false;
int graphics_info_t::n_particles = 600;
Mesh graphics_info_t::mesh_for_particles = Mesh("mesh for particles");
particle_container_t graphics_info_t::particles;
bool graphics_info_t::setup_draw_for_particles_semaphore = false;
glm::vec3 graphics_info_t::identification_pulse_centre;
bool graphics_info_t::particles_have_been_shown_already_for_this_round_flag = false;
int graphics_info_t::wait_for_hooray_refinement_tick_id = -1; // delete this tick function on refinement
                                                              // shutdown

std::vector<std::pair<glm::vec3, glm::vec3> > graphics_info_t::hydrogen_bonds_atom_position_pairs;

std::chrono::time_point<std::chrono::high_resolution_clock> graphics_info_t::tick_hydrogen_bond_mesh_t_previous = std::chrono::high_resolution_clock::now();


fun::boids_container_t graphics_info_t::boids;
Mesh graphics_info_t::mesh_for_boids;
LinesMesh graphics_info_t::lines_mesh_for_boids_box;

Mesh graphics_info_t::mesh_for_hydrogen_bonds;;

LinesMesh graphics_info_t::lines_mesh_for_identification_pulse;
LinesMesh graphics_info_t::lines_mesh_for_delete_item_pulse;
std::vector<glm::vec3> graphics_info_t::delete_item_pulse_centres;

std::vector<atom_label_info_t> graphics_info_t::labels;
TextureMesh graphics_info_t::tmesh_for_labels = TextureMesh("tmesh-for-labels");
HUDMesh graphics_info_t::mesh_for_hud_geometry = HUDMesh("hud-geometry");

TextureMesh graphics_info_t::tmesh_for_happy_face_residues_markers =
   TextureMesh("tmesh-for-happy-faces");
Texture graphics_info_t::texture_for_happy_face_residue_marker;
std::vector<glm::vec3> graphics_info_t::happy_face_residue_marker_starting_positions;

HUDTextureMesh graphics_info_t::tmesh_for_hud_geometry_tooltip_label =
   HUDTextureMesh("tmesh-for-hud-geometry-tooltip-labels");

float graphics_info_t::pull_restraint_neighbour_displacement_max_radius = 1.0;

coot::command_history_t graphics_info_t::command_history;

std::vector<Instanced_Markup_Mesh> graphics_info_t::instanced_meshes;

Texture graphics_info_t::texture_for_hud_geometry_labels;
Texture graphics_info_t::texture_for_hud_tooltip_background;
bool graphics_info_t::draw_hud_tooltip_flag = false;

HUDTextureMesh graphics_info_t::mesh_for_hud_geometry_labels = HUDTextureMesh("tmesh-for-hud-geometry-labels");
HUDTextureMesh graphics_info_t::mesh_for_hud_tooltip_background = HUDTextureMesh("tmesh-for-hud-tooltip-background");

Instanced_Markup_Mesh graphics_info_t::rama_balls_mesh = Instanced_Markup_Mesh("rama-balls");
bool graphics_info_t::draw_stick_mode_atoms_default = true;

std::string graphics_info_t::label_for_hud_geometry_tooltip;
