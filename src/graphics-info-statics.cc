
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
meshed_generic_display_object graphics_info_t::mesh_for_environment_distances;
std::chrono::time_point<std::chrono::system_clock> graphics_info_t::previous_frame_time = std::chrono::high_resolution_clock::now();
std::chrono::time_point<std::chrono::system_clock> graphics_info_t::previous_frame_time_for_per_second_counter = std::chrono::high_resolution_clock::now();
long graphics_info_t::frame_counter = 0;
long graphics_info_t::frame_counter_at_last_display = 0;
std::queue<std::chrono::time_point<std::chrono::system_clock> > graphics_info_t::frame_draw_queue;

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
int graphics_info_t::n_particles = 300;
Mesh graphics_info_t::mesh_for_particles;
particle_container_t graphics_info_t::particles;
glm::vec3 graphics_info_t::identification_pulse_centre;

fun::boids_container_t graphics_info_t::boids;
Mesh graphics_info_t::mesh_for_boids;
LinesMesh graphics_info_t::lines_mesh_for_boids_box;

LinesMesh graphics_info_t::lines_mesh_for_identification_pulse;
LinesMesh graphics_info_t::lines_mesh_for_delete_item_pulse;
std::vector<glm::vec3> graphics_info_t::delete_item_pulse_centres;

std::vector<atom_label_info_t> graphics_info_t::labels;
TextureMesh graphics_info_t::tmesh_for_labels = TextureMesh("tmesh-for-labels");
HUDMesh graphics_info_t::mesh_for_hud_geometry = HUDMesh("hud-geometry");

float graphics_info_t::pull_restraint_neighbour_displacement_max_radius = 1.0;

coot::command_history_t graphics_info_t::command_history;

std::vector<Instanced_Markup_Mesh> graphics_info_t::instanced_meshes;

Texture graphics_info_t::texture_for_hud_geometry_labels;

HUDTextureMesh graphics_info_t::mesh_for_hud_geometry_labels = HUDTextureMesh("tmesh-for-hud-geometry-labels");

Instanced_Markup_Mesh graphics_info_t::rama_balls_mesh = Instanced_Markup_Mesh("rama-balls");
bool graphics_info_t::draw_stick_mode_atoms_default = true;

