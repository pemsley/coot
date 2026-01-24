/*
 * src/graphics-info-statics.cc
 *
 * Copyright 2019 by Medical Research Council
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

#ifdef USE_PYTHON
#include "Python.h"
#endif

#include "utils/xdg-base.hh"
#include "graphics-info.h"
#include "rotate-translate-modes.hh"
#include "rotamer-search-modes.hh"

#if !defined WINDOWS_MINGW

int graphics_info_t::scale_down_graphics = 1;
int graphics_info_t::scale_up_graphics = 1;

#ifdef USE_GUILE
#ifdef USE_GUILE_GTK
bool graphics_info_t::prefer_python = 0; // prefer python scripts when
					 // scripting (if we have a
					 // choice). Default: no.
#else
#ifdef USE_PYTHON
bool graphics_info_t::prefer_python = 1; // prefer python (gui) when no
                                         // guile gui. Fixes (place-strand-here-gui)
                                         // problem when there is no guile-gtk
#else
// guile but not guile-gui, and not python
bool graphics_info_t::prefer_python = 0;
#endif // USE_PYTHON test
#endif // USE_GUILE_GTK test

#ifdef USE_GUILE
bool graphics_info_t::scm_boot_guile_booted = false; // false until my_wrap_scm_boot_guile() has been run
#endif

#else // USE_GUILE test (no guile path)
#ifdef USE_PYTHON
bool graphics_info_t::prefer_python = 1; // Python, not guile
#else
bool graphics_info_t::prefer_python = 0; // no GUILE or PYTHON
#endif // python test
#endif // USE_GUILE test
#else // Windows test (windows path)
bool graphics_info_t::prefer_python = 1; // Default: yes in Windows
#endif // windows test

bool graphics_info_t::graphics_is_gl_es = false;

// bool graphics_info_t::using_trackpad = false;
bool graphics_info_t::use_primary_mouse_for_view_rotation_flag = false;

bool graphics_info_t::use_gemmi = false;
short int graphics_info_t::python_at_prompt_flag = 0;

GtkApplication *graphics_info_t::application = 0;

int graphics_info_t::show_paths_in_display_manager_flag = 0;
std::vector<std::string> graphics_info_t::command_line_scripts;
coot::command_line_commands_t graphics_info_t::command_line_commands;
std::vector<std::string> graphics_info_t::command_line_accession_codes;

std::vector<coot::lsq_range_match_info_t> *graphics_info_t::lsq_matchers;
std::vector<coot::generic_text_object_t> graphics_info_t::generic_texts;
std::vector<coot::view_info_t> graphics_info_t::views;
bool graphics_info_t::do_expose_swap_buffers_flag = 1;

#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
ctpl::thread_pool graphics_info_t::static_thread_pool(coot::get_max_number_of_threads());
#endif // HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

clipper::Xmap<float> *graphics_info_t::dummy_xmap = new clipper::Xmap<float>;
std::vector<std::pair<std::string, clipper::Xmap<float> > > graphics_info_t::map_partition_results;
int graphics_info_t::map_partition_results_state = 0; // inactive
std::string graphics_info_t::map_partition_results_state_string; // "Done A Chain" etc.

// logging
unsigned int graphics_info_t::logging_line_index = 0;


//WII
#ifdef WII_INTERFACE_WIIUSE
wiimote** graphics_info_t::wiimotes = NULL;
#endif

float graphics_info_t::view_rotation_per_pixel_scale_factor = 1.0;

// Views
float graphics_info_t::views_play_speed = 10.0;

// movies
std::string graphics_info_t::movie_file_prefix = "movie_";
int graphics_info_t::movie_frame_number = 0;
int graphics_info_t::make_movie_flag = 0;


// LSQ
short int graphics_info_t::in_lsq_plane_deviation = 0;
short int graphics_info_t::in_lsq_plane_define    = 0;
GtkWidget *graphics_info_t::lsq_plane_dialog = 0;
std::vector<clipper::Coord_orth> *graphics_info_t::lsq_plane_atom_positions;
std::string graphics_info_t::lsq_match_chain_id_ref;
std::string graphics_info_t::lsq_match_chain_id_mov;
int graphics_info_t::lsq_ref_imol = -1;
int graphics_info_t::lsq_mov_imol = -1;
// LSQ dialog values
lsq_dialog_values_t graphics_info_t::lsq_dialog_values;



// side by side stereo?
bool graphics_info_t::in_side_by_side_stereo_mode = false;
bool graphics_info_t::in_wall_eyed_side_by_side_stereo_mode = false;

// display list for maps?
short int graphics_info_t::display_lists_for_maps_flag = 0;

//
int graphics_info_t::save_imol = -1;

bool graphics_info_t::cryo_EM_refinement_flag = false;
double graphics_info_t::geman_mcclure_alpha = 1; // soft, (20180230-PE was 2, too soft, I think)
double graphics_info_t::lennard_jones_epsilon = 2.0; // 20181008-PE less soft than 0.5
double graphics_info_t::log_cosh_target_distance_scale_factor = 2000.0;

// accept/reject
GtkWidget *graphics_info_t::accept_reject_dialog = 0;

// refinement control dialog (so that we don't get multiple copies)
GtkWidget *graphics_info_t::refine_params_dialog = 0;

// flag to display the accep/reject dialog in the toolbar
int graphics_info_t::accept_reject_dialog_docked_flag = coot::DIALOG;

// flag to show/hide/sensitise the docked accept/reject dialog in the toolbar
int graphics_info_t::accept_reject_dialog_docked_show_flag = coot::DIALOG_DOCKED_SHOW;

// the refinement toolbar show/hide
short int graphics_info_t::model_toolbar_show_hide_state = 1;

// the refinement toolbar position
short int graphics_info_t::model_toolbar_position_state = coot::model_toolbar::RIGHT;

// the refinement toolbar style
short int graphics_info_t::model_toolbar_style_state = 1;

// the main toolbar show/hide
short int graphics_info_t::main_toolbar_show_hide_state = 1;

// the main toolbar position
// // not using (yet)
//short int graphics_info_t::main_toolbar_position_state = coot::main_toolbar::TOP;

// the main toolbar style
short int graphics_info_t::main_toolbar_style_state = 2;

// refmac option menu
int graphics_info_t::refmac_molecule = -1;

//
short int graphics_info_t::active_map_drag_flag = 1; // true
long int graphics_info_t::Frames = 0;  // These 2 are to measure graphics speed.
long int graphics_info_t::T0 = 0;
bool    graphics_info_t::show_fps_flag = false;
int    graphics_info_t::control_is_pressed = 0 ; // false
short int graphics_info_t::control_key_for_rotate_flag = 1;
short int graphics_info_t::pick_pending_flag = 0;
int    graphics_info_t::a_is_pressed = 0;
int    graphics_info_t::shift_is_pressed = 0 ;       // false
int    graphics_info_t::y_is_pressed = 0 ;       // false
int    graphics_info_t::z_is_pressed = 0 ;       // false
std::pair<double, double> graphics_info_t::mouse_begin         = std::pair<double, double> (0,0);
std::pair<double, double> graphics_info_t::mouse_clicked_begin = std::pair<double, double> (0,0);
float  graphics_info_t::rotation_centre_x = 0.0;
float  graphics_info_t::rotation_centre_y = 0.0;
float  graphics_info_t::rotation_centre_z = 0.0;
// float  graphics_info_t::old_rotation_centre_x = 0.0;
// float  graphics_info_t::old_rotation_centre_y = 0.0;
// float  graphics_info_t::old_rotation_centre_z = 0.0;
coot::Cartesian graphics_info_t::old_rotation_centre(0,0,0);
clipper::Coord_orth graphics_info_t::hole_start = clipper::Coord_orth(0,0,0);
clipper::Coord_orth graphics_info_t::hole_end   = clipper::Coord_orth(0,0,0);

float  graphics_info_t::zoom                = 100;
int    graphics_info_t::smooth_scroll       =   1; // flag: default is ..
int    graphics_info_t::smooth_scroll_n_steps =  20;
float  graphics_info_t::smooth_scroll_limit =  20.0; // A
float  graphics_info_t::smooth_scroll_zoom_limit = 30.0; // A
int    graphics_info_t::smooth_scroll_do_zoom = 0;  // initially no, too ugly ATM.
short int graphics_info_t::smooth_scroll_on = 0;
int    graphics_info_t::smooth_scroll_current_step = 0;
coot::Cartesian graphics_info_t::smooth_scroll_delta;
int    graphics_info_t::mouse_just_cliked     = 0;
float  graphics_info_t::user_defined_rotation_centre_crosshairs_size_scale_factor = 0.05;
glm::vec4 graphics_info_t::rotation_centre_cross_hairs_colour = glm::vec4(0.8, 0.8, 0.8, 1.0);
short int graphics_info_t::quanta_like_zoom_flag = 0;
int    graphics_info_t::go_to_ligand_animate_view_n_steps = 50;

// graphics display size:
int graphics_info_t::graphics_x_size = 500; // 20220602-PE these are olden values
int graphics_info_t::graphics_y_size = 500; // how do they tie into new graphics?

int graphics_info_t::graphics_x_position = 0; // gets overwritten at start
int graphics_info_t::graphics_y_position = 0;
int graphics_info_t::model_fit_refine_dialog_stays_on_top_flag = 1;

bool graphics_info_t::use_graphics_interface_flag = false;

// display control size and position and the vboxes for the maps and
// the molecules in the display control
int graphics_info_t::display_manager_x_size = -1;
int graphics_info_t::display_manager_y_size = -1;
int graphics_info_t::display_manager_x_position = -1;
int graphics_info_t::display_manager_y_position = -1;

int graphics_info_t::display_manager_molecules_vbox_x_size = -1;
int graphics_info_t::display_manager_molecules_vbox_y_size = -1;
int graphics_info_t::display_manager_maps_vbox_x_size = -1;
int graphics_info_t::display_manager_maps_vbox_y_size = -1;

int graphics_info_t::display_manager_paned_position = -1;


int graphics_info_t::go_to_atom_window_x_position = -100;
int graphics_info_t::go_to_atom_window_y_position = -100;

int graphics_info_t::delete_item_widget_x_position = -100;
int graphics_info_t::delete_item_widget_y_position = -100;

int graphics_info_t::accept_reject_dialog_x_position = -100;
int graphics_info_t::accept_reject_dialog_y_position = -100;

int graphics_info_t::edit_chi_angles_dialog_x_position = -1;
int graphics_info_t::edit_chi_angles_dialog_y_position = -1;

short int graphics_info_t::draw_chi_angle_flash_bond_flag = 0;
std::pair<clipper::Coord_orth, clipper::Coord_orth> graphics_info_t::flash_bond =
   std::pair<clipper::Coord_orth, clipper::Coord_orth> (clipper::Coord_orth(0,0,0),
							clipper::Coord_orth(0,0,0));
bool graphics_info_t::flash_intermediate_atom_pick_flag = 0;
clipper::Coord_orth graphics_info_t::intermediate_flash_point;


int graphics_info_t::default_bond_width = 5;
bool graphics_info_t::use_variable_bond_width = false;


int graphics_info_t::rotamer_selection_dialog_x_position = -100;
int graphics_info_t::rotamer_selection_dialog_y_position = -100;

int graphics_info_t::model_fit_refine_x_position = -100;  // initially unset
int graphics_info_t::model_fit_refine_y_position = -100;

// save rotate/translate widget position:
int graphics_info_t::rotate_translate_x_position = -100;  // initially unset
int graphics_info_t::rotate_translate_y_position = -100;

// save Ramachandran widget position:
int graphics_info_t::ramachandran_plot_x_position = -100;  // initially unset
int graphics_info_t::ramachandran_plot_y_position = -100;

coot::Cartesian graphics_info_t::distance_pos_1 = coot::Cartesian(0.0, 0.0, 0.0);
coot::Cartesian graphics_info_t::angle_tor_pos_1 = coot::Cartesian(0.0, 0.0, 0.0);
coot::Cartesian graphics_info_t::angle_tor_pos_2 = coot::Cartesian(0.0, 0.0, 0.0);
coot::Cartesian graphics_info_t::angle_tor_pos_3 = coot::Cartesian(0.0, 0.0, 0.0);
coot::Cartesian graphics_info_t::angle_tor_pos_4 = coot::Cartesian(0.0, 0.0, 0.0);

// Shall we have file name filtering (i.e. before fileselection is
// displayed) on by default?
int graphics_info_t::filter_fileselection_filenames_flag = 1; // yes
int graphics_info_t::file_chooser_dialog_x_size = -1; // unset
int graphics_info_t::file_chooser_dialog_y_size = -1;


// things for quaternion-based view rotation:
double graphics_info_t::mouse_current_x = 0.0;
double graphics_info_t::mouse_current_y = 0.0;
// float* graphics_info_t::quat = new float[4]; // gone. Use glm_quat
float graphics_info_t::trackball_size = 12.8; // better for me for now. was 0.8

// 20220606-PE new style mouse positioning - because of event controllers
double graphics_info_t::drag_begin_x = 0.0;
double graphics_info_t::drag_begin_y = 0.0;
double graphics_info_t::mouse_x = 0.0;
double graphics_info_t::mouse_y = 0.0;
double graphics_info_t::mouse_speed = 2.0;
std::pair<double, double> graphics_info_t::mouse_previous_position = std::make_pair(0.0, 0.0);

// residue reorientation on "space"
bool graphics_info_t::reorienting_next_residue_mode = false;


// things for baton quaternion rotation: Must use a c++ class at some
// stage:
// baton
// float* graphics_info_t::baton_quat = new float[4]; // gone - needs glm::quat version
float  graphics_info_t::baton_length = 3.8; // A;
int    graphics_info_t::baton_next_ca_options_index = 0;
int    graphics_info_t::user_set_baton_atom_molecule = -1; // initially unset
short int graphics_info_t::baton_tmp_atoms_to_new_molecule = 0;
short int graphics_info_t::draw_baton_flag = 0; // off initially
short int graphics_info_t::baton_mode = 0;      // rotation mode
coot::Cartesian graphics_info_t::baton_root = coot::Cartesian(0.0, 0.0, 0.0);
coot::Cartesian graphics_info_t::baton_tip =  coot::Cartesian(0.0, 0.0, 3.8);
std::vector<coot::scored_skel_coord> graphics_info_t::baton_next_ca_options;
std::vector<clipper::Coord_orth> graphics_info_t::baton_previous_ca_positions;
short int graphics_info_t::baton_build_direction_flag = 1; // forwards by default
int graphics_info_t::baton_build_start_resno = 1;
short int graphics_info_t::baton_build_params_active = 0; // not active initially.
std::string graphics_info_t::baton_build_chain_id = std::string("");

// place helix here
float graphics_info_t::place_helix_here_fudge_factor = 1.0; // (it's multiplicative)

// double*  graphics_info_t::symm_colour_merge_weight = new double[10];
// double **graphics_info_t::symm_colour = new double*[10];

double graphics_info_t::symmetry_colour_merge_weight = 0.5; // 0.0 -> 1.0

glm::vec4 graphics_info_t::symmetry_colour = glm::vec4(0.4, 0.4, 0.4, 1.0);

double*  graphics_info_t::skeleton_colour = new double[4];
int      graphics_info_t::map_for_skeletonize = -1;
int      graphics_info_t::max_skeleton_search_depth = 10;

short int graphics_info_t::swap_pre_post_refmac_map_colours_flag = 0;

float     graphics_info_t::symmetry_search_radius = 13.0;
int       graphics_info_t::symmetry_shift_search_size = 1; // which_boxes search hack

// short int graphics_info_t::symmetry_as_calphas = 0; // moved to per molecule basis
// short int graphics_info_t::symmetry_rotate_colour_map_flag = 0; // moved to per molecule basis
float     graphics_info_t::symmetry_operator_rotate_colour_map = 60.0; //degrees
// int       graphics_info_t::symmetry_colour_by_symop_flag = 1; // moved to per molecule basis
// int       graphics_info_t::symmetry_whole_chain_flag = 0;  // moved to per molecule basis

// esoteric depth cue?
int       graphics_info_t::esoteric_depth_cue_flag = 1; // on by default.

// save coords fileselection dir
int graphics_info_t::save_coordinates_in_original_dir_flag = 0;

// save CONECT records, by default we dont
int graphics_info_t::write_conect_records_flag = 0;

// by default convert nucleic acid names to match the (currently v2)
// dictionary.
//
bool graphics_info_t::convert_to_v2_atom_names_flag = 0; // changed 20110505

//
short int graphics_info_t::print_initial_chi_squareds_flag = 0;

short int graphics_info_t::show_symmetry = 0;

float    graphics_info_t::box_radius_xray =  20.0;
float    graphics_info_t::box_radius_em   =  80.0;


int      graphics_info_t::debug_atom_picking = 0;

int   graphics_info_t::imol_map_sharpening = -1;
float graphics_info_t::map_sharpening_scale_limit = 200.0;

//
int graphics_info_t::imol_remarks_browswer = -1;

// dragged moving atom:
int       graphics_info_t::moving_atoms_currently_dragged_atom_index = -1;
short int graphics_info_t::in_moving_atoms_drag_atom_mode_flag = 0;
std::set<int> empty_int_set;
std::set<int> graphics_info_t::moving_atoms_dragged_atom_indices = empty_int_set;

// validate moving atoms
int       graphics_info_t::moving_atoms_n_cis_peptides = -1;  // unset

// for picking intermediate atoms
bool      graphics_info_t::moving_atoms_have_hydrogens_displayed = false;


std::string graphics_info_t::model_fit_refine_place_atom_at_pointer_string = "";
std::string graphics_info_t::model_fit_refine_rotate_translate_zone_string = "";

// Change of plan, so that we are more compatible with Stuart.
//
// We shall make maps and skeletons be objects of molecules.
// Hmmm... skeletons are objects of maps.
//
// And molecules are instances of a molecule_class_info_t
//
// We need to store how many molecules we have and where to find them.

int graphics_info_t::n_molecules_max = 60;
// int graphics_info_t::n_molecules = 0; // gets incremented on pdb reading

// generic display objects, gets set in init.
// std::vector<coot::old_generic_display_object_t> *graphics_info_t::generic_objects_p = NULL;
std::vector<meshed_generic_display_object> graphics_info_t::generic_display_objects;
bool graphics_info_t::display_generic_objects_as_solid_flag = 0;
bool graphics_info_t::display_environment_graphics_object_as_solid_flag = 0;
GtkWidget *graphics_info_t::generic_objects_dialog = NULL;


coot::console_display_commands_t graphics_info_t::console_display_commands;

// 20 molecules is enough for anyone, surely?
// FIXME - we should be using dynamic allocation, and allocating
// an new array of molecule_class_info_ts and moving over the pointers.
//
// Fixed ( see graphics_info_t::initialize_graphics_molecules(); )
//
// molecule_class_info_t* graphics_info_t::molecules = NULL; yesterday's array

// The molecule for undoing
int graphics_info_t::undo_molecule = -1;

// backup filenames
bool graphics_info_t::unpathed_backup_file_names_flag = 0;
#ifdef WINDOWS_MINGW
bool graphics_info_t::decoloned_backup_file_names_flag = 1;
#else
bool graphics_info_t::decoloned_backup_file_names_flag = 0;
#endif

// backup compress files (default: compress)
int graphics_info_t::backup_compress_files_flag = 1;

// Auto read
int graphics_info_t::auto_read_do_difference_map_too_flag = 1;

// nomenclature errors
coot::nomenclature_error_handle_type graphics_info_t::nomenclature_errors_mode = coot::PROMPT;




// Tip of the Day?
short int graphics_info_t::do_tip_of_the_day_flag = 1;

// Browser URL
// Also gets get in group? scheme code
std::string graphics_info_t::browser_open_command = "firefox -remote";

// should be: each map, each contour level, each colour
// (triple star)
//

std::vector<GtkWidget *> graphics_info_t::glareas;
int graphics_info_t::hud_start_graphics_window_x_width = -1;  // unset
int graphics_info_t::hud_start_graphics_window_x_height = -1;

GtkWidget *graphics_info_t::statusbar = NULL;
guint      graphics_info_t::statusbar_context_id = 0;
std::string graphics_info_t::main_window_title;

//
int       graphics_info_t::atom_label_font_size = 2; // medium
void     *graphics_info_t::atom_label_font = 0; //GLUT_BITMAP_HELVETICA_12;
int       graphics_info_t::label_atom_on_recentre_flag = 1;
int       graphics_info_t::symmetry_atom_labels_expanded_flag = 0;
coot::colour_holder graphics_info_t::font_colour = coot::colour_holder(1.0, 0.8, 0.8);

bool      graphics_info_t::stroke_characters = false;

short int graphics_info_t::brief_atom_labels_flag = 0;
short int graphics_info_t::seg_ids_in_atom_labels_flag = 0;

// scroll wheel
int       graphics_info_t::scroll_wheel_map = -1; // (initial magic value)
                                                  // updated on new read-map.

//
GLuint theMapContours = 0;

float graphics_info_t::iso_level_increment = 0.05;
float graphics_info_t::diff_map_iso_level_increment =  0.005;
short int graphics_info_t::swap_difference_map_colours = 0; // default: not in Jan-Dohnalek-mode.

// No idle functions to start (but setting them to zero doesn't set that - the
// idle functions are added by gtk_idle_add()).
int   graphics_info_t::idle_function_spin_rock_token = -1; // magic "unset" value
std::chrono::time_point<std::chrono::high_resolution_clock> graphics_info_t::time_holder_for_rocking = std::chrono::high_resolution_clock::now();

double graphics_info_t::idle_function_rock_amplitude_scale_factor = 1.0;
double graphics_info_t::idle_function_rock_freq_scale_factor = 1.0;
double graphics_info_t::idle_function_rock_angle_previous = 0;

std::vector<std::chrono::time_point<std::chrono::high_resolution_clock> > graphics_info_t::leftquote_press_times;

#ifdef USE_PYTHON
// Hamish python
std::string graphics_info_t::python_draw_function_string;
#endif

// new style (20110505 ligand interactions)
//
long  graphics_info_t::time_holder_for_ligand_interactions = 0;
int   graphics_info_t::idle_function_ligand_interactions_token = 0;
double graphics_info_t::ligand_interaction_pulse_previous = 0;


int   graphics_info_t::drag_refine_idle_function_token = -1; // magic unused value
coot::refinement_results_t graphics_info_t::saved_dragged_refinement_results(0, -2, "");
float graphics_info_t::idle_function_rotate_angle = 1.0; // degrees
bool  graphics_info_t::refinement_move_atoms_with_zero_occupancy_flag = 1; // yes

bool graphics_info_t::post_intermediate_atoms_moved_ready;
#ifdef USE_PYTHON
PyObject *graphics_info_t::post_intermediate_atoms_moved_hook = 0;
#endif


float graphics_info_t::map_sampling_rate = 2.5;

// Initialise the static atom_sel.
//
//We use this to store the atom selections of the molecules (not that
//an atom_selection_container_t contains the atom selection and the
//molecular manager (or at least, pointers to them).
//
// Only one atom_selection_container_t at the moment.  Later perhaps
// we will make the atom_selection_container_t* for multiple
// molecules.
//
//
//atom_selection_container_t aaa;
//atom_selection_container_t molecule_class_info_t::atom_sel = aaa;

//

short int graphics_info_t::show_aniso_atoms_flag = 0; // initially don't show.
short int graphics_info_t::show_aniso_atoms_radius_flag = 0;
float     graphics_info_t::show_aniso_atoms_radius = 12.0;
float     graphics_info_t::show_aniso_atoms_probability = 0.5; // 20250602-PE 0.0 to 1.10 now

// initialise the molecule (scene) rotation axis statics.
//
float molecule_rot_t::x_axis_angle = 0.0;
float molecule_rot_t::y_axis_angle = 0.0;

// 0: never run it
// 1: ask to run it
// 2: always run it
short int graphics_info_t::run_state_file_status = 1;
bool      graphics_info_t::state_file_was_run_flag = false;
// did we start with --no-startup-scripts?
bool      graphics_info_t::run_startup_scripts_flag = true;

GtkWidget *graphics_info_t::preferences_widget = NULL;
int        graphics_info_t::mark_cis_peptides_as_bad_flag = 1;

std::vector<std::string> graphics_info_t::preferences_general_tabs;
std::vector<std::string> graphics_info_t::preferences_bond_tabs;
std::vector<std::string> graphics_info_t::preferences_geometry_tabs;
std::vector<std::string> graphics_info_t::preferences_colour_tabs;
std::vector<std::string> graphics_info_t::preferences_map_tabs;
std::vector<std::string> graphics_info_t::preferences_other_tabs;
std::vector<coot::preferences_icon_info_t> *graphics_info_t::model_toolbar_icons;
std::vector<coot::preferences_icon_info_t> *graphics_info_t::main_toolbar_icons;

std::vector<coot::preference_info_t> graphics_info_t::preferences_internal;
std::vector<coot::preference_info_t> graphics_info_t::preferences_internal_default;

// Torsions with hydrogens?
bool graphics_info_t::find_hydrogen_torsions_flag = 0; // no

// Which ccp4i project shall we put at the top of a fileselection optionmenu?
//
int graphics_info_t::ccp4_projects_index_last = 0;

// phs reading
std::string graphics_info_t::phs_filename = "";

// Pointer Distances
float graphics_info_t::pointer_min_dist = 0.1;
float graphics_info_t::pointer_max_dist = 3.6;
bool graphics_info_t::show_pointer_distances_flag = false;

// restraints that have less penalty/energy than
// this are not worth drawing.
float graphics_info_t::extra_distance_restraint_penalty_cutoff = 0.01; // draw them all to start with
bool graphics_info_t::show_extra_distance_restraints_flag = true;
std::vector<extra_distance_restraint_markup_instancing_data_t> graphics_info_t::extra_distance_restraints_markup_data;
Mesh graphics_info_t::mesh_for_extra_distance_restraints = Mesh("mesh_for_extra_distance_restraints"); // draw this with instancing


// Go to Atom widget:
//
std::string graphics_info_t::go_to_atom_chain_     =  "A";
std::string graphics_info_t::go_to_atom_atom_name_ = "CA";
std::string graphics_info_t::go_to_atom_atom_altLoc_ = "";
std::string graphics_info_t::go_to_atom_inscode_ = "";
int         graphics_info_t::go_to_atom_residue_   = -9999; // magic
							    // number. unset
							    // initially.
int         graphics_info_t::go_to_atom_molecule_  = 0;
int         graphics_info_t::go_to_ligand_n_atoms_limit = 6;

int         graphics_info_t::go_to_atom_mol_menu_active_position = -1; // unset
                                                                       // initially.
GtkWidget  *graphics_info_t::go_to_atom_window = NULL;
int         graphics_info_t::go_to_atom_menu_label_n_chars_max = 40;

GtkWidget *graphics_info_t::model_fit_refine_dialog = NULL;
short int  graphics_info_t::model_fit_refine_dialog_was_sucked = 0;
GtkWidget *graphics_info_t::residue_info_dialog = NULL;
GtkWidget *graphics_info_t::difference_map_peaks_dialog = NULL;
GtkWidget *graphics_info_t::checked_waters_baddies_dialog = NULL;

int graphics_info_t::in_base_paring_define = 0;

GtkWidget *graphics_info_t::other_modelling_tools_dialog = 0;

coot::residue_spec_t graphics_info_t::current_residue = coot::residue_spec_t();

// Skeleton Widgets:
float graphics_info_t::skeleton_level = 0.2;
float graphics_info_t::skeleton_box_radius = 40.0;

// Autobuild control
short int graphics_info_t::autobuild_flag = 0;

// Fileselection sorting:
short int graphics_info_t::sticky_sort_by_date = 0; // initally not.

// Maps and Molecule Display Control window:
GtkWidget *graphics_info_t::display_control_window_ = NULL;

int graphics_info_t::draw_axes_flag = 1; // on by default now.

bool graphics_info_t::regenerate_bonds_needs_make_bonds_type_checked_flag = true;

//
short int graphics_info_t::draw_crosshairs_flag = 0;

// For defining a range (to, say, regularize)
//
int       graphics_info_t::refine_regularize_max_residues = 20;
int       graphics_info_t::refine_auto_range_step = 1; // +/- 1 about the clicked residue
coot::atom_spec_t graphics_info_t::in_range_first_picked_atom  = coot::atom_spec_t(nullptr);
coot::atom_spec_t graphics_info_t::in_range_second_picked_atom = coot::atom_spec_t(nullptr);
short int graphics_info_t::in_range_define = 0; // regularization
short int graphics_info_t::in_range_define_for_refine = 0;// refine (i.e. with map)
short int graphics_info_t::in_distance_define = 0;
short int graphics_info_t::in_angle_define = 0;
short int graphics_info_t::in_torsion_define = 0;
short int graphics_info_t::fix_chiral_volume_before_refinement_flag = 1;
int       graphics_info_t::check_chiral_volume_molecule = 0;
int       graphics_info_t::add_reps_molecule_option_menu_item_select_molecule = 0; // option menu
int       graphics_info_t::add_reps_molecule_combobox_molecule = 0;

short int graphics_info_t::refinement_immediate_replacement_flag = 0;
short int graphics_info_t::show_chiral_volume_errors_dialog_flag = 1; // on by default


int graphics_info_t::residue_range_mol_no  = 0;
int graphics_info_t::residue_range_atom_index_1 = 0;
int graphics_info_t::residue_range_atom_index_2 = 0;
int graphics_info_t::geometry_atom_index_1 = 0;
int graphics_info_t::geometry_atom_index_2 = 0;
int graphics_info_t::geometry_atom_index_3 = 0;
int graphics_info_t::geometry_atom_index_4 = 0;
int graphics_info_t::geometry_atom_index_1_mol_no = -1; // must be set before use.
int graphics_info_t::geometry_atom_index_2_mol_no = -1;
int graphics_info_t::geometry_atom_index_3_mol_no = -1;
int graphics_info_t::geometry_atom_index_4_mol_no = -1;

// torsion general
int graphics_info_t::torsion_general_atom_index_1 = -1;
int graphics_info_t::torsion_general_atom_index_2 = -1;
int graphics_info_t::torsion_general_atom_index_3 = -1;
int graphics_info_t::torsion_general_atom_index_4 = -1;
int graphics_info_t::torsion_general_atom_index_1_mol_no = -1;
int graphics_info_t::torsion_general_atom_index_2_mol_no = -1;
int graphics_info_t::torsion_general_atom_index_3_mol_no = -1;
int graphics_info_t::torsion_general_atom_index_4_mol_no = -1;
short int graphics_info_t::in_edit_torsion_general_flag = 0;
short int graphics_info_t::in_residue_partial_alt_locs_define = 0;
int graphics_info_t::imol_residue_partial_alt_locs = 0;
coot::residue_spec_t graphics_info_t::residue_partial_alt_locs_spec;
double graphics_info_t::residue_partial_alt_locs_rotate_fragment_angle = 40; // degrees

std::vector<coot::atom_spec_t> graphics_info_t::torsion_general_atom_specs;
bool graphics_info_t::torsion_general_reverse_flag = 0;
Tree graphics_info_t::torsion_general_tree;
std::vector<std::vector<int> > graphics_info_t::torsion_general_contact_indices;



//
short int graphics_info_t::in_residue_info_define = 0;
int graphics_info_t::residue_selection_flash_frames_number = 2; // was 3
short int graphics_info_t::in_torsion_general_define = 0;

short int graphics_info_t::in_save_symmetry_define = 0;

short int graphics_info_t::in_rot_trans_object_define = 0;
short int graphics_info_t::rot_trans_object_type = ROT_TRANS_TYPE_ZONE;
short int graphics_info_t::rot_trans_zone_rotates_about_zone_centre = 0;
int graphics_info_t::rot_trans_atom_index_1 = -1;
int graphics_info_t::rot_trans_atom_index_2 = -1;
int graphics_info_t::imol_rot_trans_object = -1;
// int graphics_info_t::rot_trans_atom_index_rotation_origin_atom = -1; // old
mmdb::Atom *graphics_info_t::rot_trans_rotation_origin_atom = NULL;

float *graphics_info_t::previous_rot_trans_adjustment = new float[6];

short int graphics_info_t::in_db_main_define = 0;
int graphics_info_t::db_main_imol = -1;
int graphics_info_t::db_main_atom_index_1 = -1;
int graphics_info_t::db_main_atom_index_2 = -1;
coot::db_main graphics_info_t::main_chain;

coot::fixed_atom_pick_state_t graphics_info_t::in_fixed_atom_define = coot::FIXED_ATOM_NO_PICK;
// GtkWidget *graphics_info_t::fixed_atom_dialog = 0; old style. Now we look it up each time

std::vector<coot::simple_distance_object_t> graphics_info_t::measure_distance_object_vec;
std::vector<std::pair<clipper::Coord_orth, clipper::Coord_orth> > graphics_info_t::pointer_distances_object_vec;
std::vector<coot::coord_orth_triple> graphics_info_t::measure_angle_object_vec;
Mesh graphics_info_t::mesh_for_measure_distance_object_vec = Mesh("mesh-for-measure-distance-object-vec");
Mesh graphics_info_t::mesh_for_measure_angle_object_vec    = Mesh("mesh-for-measure-angle-object-vec");
std::vector<atom_label_info_t> graphics_info_t::labels_for_measure_distances_and_angles;
Mesh graphics_info_t::mesh_for_eyelashes("eyelashes");

// Mesh graphics_info_t::mesh_for_pointer_distances = Mesh("mesh-for-pointer-distances");
meshed_generic_display_object graphics_info_t::mesh_for_pointer_distances("mesh-in-generic-display-object for pointer distances");
std::vector<atom_label_info_t> graphics_info_t::labels_for_pointer_distances;

int graphics_info_t::show_origin_marker_flag = 1;

translation_gizmo_t graphics_info_t::translation_gizmo; // axes with cones at the centre of an object.
                                                        // This thing should be clickable
Mesh graphics_info_t::translation_gizmo_mesh = Mesh("translation gizmo");
bool graphics_info_t::translation_gizmo_is_being_dragged = false;
translation_gizmo_t::pick_info_t graphics_info_t::translation_gizmo_axis_dragged = translation_gizmo_t::pick_info_t::NONE;

//
float graphics_info_t::geometry_vs_map_weight = 60.0;

int graphics_info_t::refine_params_dialog_geman_mcclure_alpha_combobox_position    = 3;
int graphics_info_t::refine_params_dialog_lennard_jones_epsilon_combobox_position  = 3;
int graphics_info_t::refine_params_dialog_rama_restraints_weight_combobox_position = 3;
int graphics_info_t::refine_params_dialog_torsions_weight_combox_position = 2;
bool graphics_info_t::refine_params_dialog_extra_control_frame_is_visible = false;

int   graphics_info_t::rama_n_diffs = 50;

atom_selection_container_t *graphics_info_t::moving_atoms_asc = NULL;
short int graphics_info_t::moving_atoms_asc_type = coot::NEW_COORDS_UNSET; // unset
int graphics_info_t::imol_moving_atoms = 0;
coot::extra_restraints_representation_t graphics_info_t::moving_atoms_extra_restraints_representation;

bool graphics_info_t::noughties_physics = false;

bool graphics_info_t::draw_it_for_moving_atoms_restraints_graphics_object = false;
bool graphics_info_t::draw_it_for_moving_atoms_restraints_graphics_object_user_control = false;
int graphics_info_t::imol_refinement_map = -1; // magic initial value "None set"
                                               // checked in graphics_info_t::refine()

graphical_bonds_container graphics_info_t::regularize_object_bonds_box;
graphical_bonds_container graphics_info_t::environment_object_bonds_box;
graphical_bonds_container graphics_info_t::symmetry_environment_object_bonds_box;
int graphics_info_t::default_bonds_box_type = coot::NORMAL_BONDS; // Phil wants to change this.
int graphics_info_t::mol_no_for_environment_distances = -1;

int   graphics_info_t::bond_parameters_molecule = -1; // unset


short int graphics_info_t::do_torsion_restraints = 0;
short int graphics_info_t::do_peptide_omega_torsion_restraints = 0;
bool      graphics_info_t::do_rama_restraints = 0; // No.
bool      graphics_info_t::do_numerical_gradients = 0; // No.
int       graphics_info_t::restraints_rama_type = coot::RAMA_TYPE_LOGRAMA;
// float     graphics_info_t::rama_restraints_weight = 10.0; // clipper-rama weight, gets reset on set_refine_ramachandran_restraints_type()
float     graphics_info_t::rama_plot_restraints_weight = 1.0;  // what is this? This is what is actually used in
                                                               // in on_refine_params_rama_restraints_weight_combobox_changed_gtkbuilder_callback()
                                                               // and in make_restraints - 20220329-PE delete the above one then.

// for Kevin Keating
bool      graphics_info_t::use_only_extra_torsion_restraints_for_torsions_flag = 0;

//
bool      graphics_info_t::do_trans_peptide_restraints = true;

//
bool graphics_info_t::do_intermediate_atoms_rama_markup = true;
bool graphics_info_t::do_intermediate_atoms_rota_markup = true;

//
short int graphics_info_t::guile_gui_loaded_flag = FALSE;
short int graphics_info_t::python_gui_loaded_flag = FALSE;

//
bool  graphics_info_t::find_ligand_do_real_space_refine_ = true; // default on
int   graphics_info_t::find_ligand_map_mol_ = -1;
int   graphics_info_t::find_ligand_protein_mol_ = -1;
bool  graphics_info_t::find_ligand_here_cluster_flag = 0;
int   graphics_info_t::find_ligand_n_top_ligands = 10;
bool  graphics_info_t::find_ligand_multiple_solutions_per_cluster_flag = false;
float graphics_info_t::find_ligand_score_by_correl_frac_limit = 0.7;
float graphics_info_t::find_ligand_score_correl_frac_interesting_limit = 0.9;

double graphics_info_t::map_to_model_correlation_atom_radius = 1.5;

short int graphics_info_t::find_ligand_mask_waters_flag = 0;
float graphics_info_t::map_mask_atom_radius = -99; // unset
// std::vector<int> *mol_tmp;
std::vector<std::pair<int, bool> > *graphics_info_t::find_ligand_ligand_mols_;
float graphics_info_t::find_waters_sigma_cut_off = 1.4;
float graphics_info_t::ligand_acceptable_fit_fraction = 0.75;
float graphics_info_t::ligand_cluster_sigma_level = 1.0; // sigma
int   graphics_info_t::ligand_wiggly_ligand_n_samples = 50;
int   graphics_info_t::ligand_wiggly_ligand_count = 0; // dummy
int   graphics_info_t::ligand_verbose_reporting_flag = 0;
// std::vector<short int> *graphics_info_t::find_ligand_wiggly_ligands_; bye!
short int graphics_info_t::ligand_expert_flag = 0;

float graphics_info_t::ligand_water_to_protein_distance_lim_max = 3.2;
float graphics_info_t::ligand_water_to_protein_distance_lim_min = 2.4;
float graphics_info_t::ligand_water_variance_limit = 0.12;
int   graphics_info_t::ligand_water_n_cycles = 3;
int   graphics_info_t::find_ligand_ligand_atom_limit = 400;
// EJD wants it, but off by default:
short int graphics_info_t::ligand_water_write_peaksearched_atoms = 0;
std::vector<clipper::Coord_orth> *graphics_info_t::ligand_big_blobs = NULL;

bool graphics_info_t::graphics_ligand_view_flag = false;
int  graphics_info_t::graphics_ligand_view_imol = -1;


short int graphics_info_t::do_probe_dots_on_rotamers_and_chis_flag = 0;
short int graphics_info_t::do_probe_dots_post_refine_flag = 0;
coot::Cartesian graphics_info_t::probe_dots_on_chis_molprobity_centre = coot::Cartesian(0.0, 0.0, 0.0);
float graphics_info_t::probe_dots_on_chis_molprobity_radius = 6.0;
bool graphics_info_t::do_coot_probe_dots_during_refine_flag = false;

float grey_level = 0.24;
float norm_255 = 1.0/255.0;

// this background is too "light" when we zoom in - the depth-cueing looks bad.
// glm::vec3 graphics_info_t::background_colour = glm::vec3(51.0 * norm_255, 57.0 * norm_255, 59.0 * norm_255);
glm::vec3 graphics_info_t::background_colour = glm::vec3(17.0 * norm_255,
                                                         17.0 * norm_255,
                                                         17.0 * norm_255);
//
short int graphics_info_t::delete_item_atom = 0;
short int graphics_info_t::delete_item_residue = 0;
short int graphics_info_t::delete_item_residue_zone = 0;
short int graphics_info_t::delete_item_residue_hydrogens = 0;
short int graphics_info_t::delete_item_water = 0;
short int graphics_info_t::delete_item_sidechain = 0;
short int graphics_info_t::delete_item_sidechain_range = 0;
short int graphics_info_t::delete_item_chain = 0;
GtkWidget *graphics_info_t::delete_item_widget = NULL;
int       graphics_info_t::keep_delete_item_active_flag = 0;
coot::residue_spec_t graphics_info_t::delete_item_residue_zone_1;
coot::residue_spec_t graphics_info_t::delete_item_sidechain_range_1;
int graphics_info_t::delete_item_residue_zone_1_imol = -1;
int graphics_info_t::delete_item_sidechain_range_1_imol = -1;


GtkWidget *graphics_info_t::symmetry_controller_dialog = 0;

short int graphics_info_t::do_scroll_by_wheel_mouse_flag = 1;


// dummy atom typed (dialog) or dummy (forced/auto)
short int graphics_info_t::pointer_atom_is_dummy = 0;
int       graphics_info_t::user_pointer_atom_molecule = -1; // unset.

// read a pdb, shall we recentre?
short int graphics_info_t::recentre_on_read_pdb = 1;

// Buttons - now are configurable
//
GdkModifierType graphics_info_t::button_1_mask_ = GDK_BUTTON1_MASK;
GdkModifierType graphics_info_t::button_2_mask_ = GDK_BUTTON2_MASK;
GdkModifierType graphics_info_t::button_3_mask_ = GDK_BUTTON3_MASK;

std::map<GLchar, FT_character> graphics_info_t::ft_characters;

// shall we show the density leven on screen?
short int graphics_info_t::display_density_level_on_screen = 1;
short int graphics_info_t::display_density_level_this_image = 1;
std::string graphics_info_t::display_density_level_screen_string =
   "Welcome to Coot";

// This kills the compiler:  Move the allocation to init.
// GtkWidget **graphics_info_t::dynarama_is_displayed = new GtkWidget *[graphics_info_t::n_molecules_max];
float       graphics_info_t::residue_density_fit_scale_factor = 1.0;

// cif dictionary
std::vector<std::string> *graphics_info_t::cif_dictionary_filename_vec = NULL;
int  graphics_info_t::cif_dictionary_read_number = 1;
// bool graphics_info_t::cif_dictionary_file_selector_create_molecule_flag = true; // Too annoying
bool graphics_info_t::cif_dictionary_file_selector_create_molecule_flag = false;


// map radius slider
float graphics_info_t::map_radius_slider_max = 50.0;

// Rotate colour map
short int graphics_info_t::rotate_colour_map_on_read_pdb_flag = 1; // do it.
short int graphics_info_t::rotate_colour_map_on_read_pdb_c_only_flag = 1; // rotate Cs only by default
float     graphics_info_t::rotate_colour_map_on_read_pdb = 21.0;  // degrees
float     graphics_info_t::rotate_colour_map_for_map = 31.0; // 20240907-PE 14.0;  // degrees

float graphics_info_t::goodsell_chain_colour_wheel_rotation_step = 0.221;


// cell colour
coot::colour_holder graphics_info_t::cell_colour =
   coot::colour_holder(0.8, 0.8, 0.2);


// regulariziation
//
coot::protein_geometry *graphics_info_t::geom_p = NULL;


// rotamer probabilities
int graphics_info_t::rotamer_search_mode = ROTAMERSEARCHAUTOMATIC;
coot::rotamer_probability_tables graphics_info_t::rot_prob_tables;
float graphics_info_t::rotamer_distortion_scale = 0.3;

// PHENIX support
std::string graphics_info_t::external_refinement_program_button_label = "*-*";


// pepflip
int graphics_info_t::atom_index_pepflip = 0;
int graphics_info_t::imol_pepflip = 0;
int graphics_info_t::iresno_pepflip = 0;
short int graphics_info_t::in_pepflip_define = 0;

// rigid body refinement
short int graphics_info_t::in_rigid_body_define = 0;
int       graphics_info_t::imol_rigid_body_refine = 0;

// terminal residue define
short int graphics_info_t::in_terminal_residue_define = 0;
short int graphics_info_t::add_terminal_residue_immediate_addition_flag = 1;
short int graphics_info_t::add_terminal_residue_do_post_refine = 0;
float graphics_info_t::terminal_residue_addition_direct_phi = -135.0;
float graphics_info_t::terminal_residue_addition_direct_psi =  135.0;


// CIS <-> TRANS conversion
int graphics_info_t::in_cis_trans_convert_define = 0;

// add OXT atom
int         graphics_info_t::add_OXT_molecule = -1;
std::string graphics_info_t::add_OXT_chain;

float graphics_info_t::default_new_atoms_b_factor = 30.0;

int graphics_info_t::reset_b_factor_moved_atoms_flag = 0;

// show environment
float graphics_info_t::environment_min_distance = 0.0;
float graphics_info_t::environment_max_distance = 3.2;
bool graphics_info_t::environment_show_distances = 0;
bool graphics_info_t::environment_distances_show_bumps = 1;
bool graphics_info_t::environment_distances_show_h_bonds = 1;

short int graphics_info_t::environment_distance_label_atom = 0;

// dynamic distances to intermediate atoms:
short int graphics_info_t::in_dynamic_distance_define = 0;
coot::intermediate_atom_distance_t graphics_info_t::running_dynamic_distance;
std::vector<coot::intermediate_atom_distance_t> graphics_info_t::dynamic_distances;

//
bool graphics_info_t::disable_state_script_writing = 0;

// Dynamic map resampling and sizing
short int graphics_info_t::dynamic_map_resampling = 0;
short int graphics_info_t::dynamic_map_size_display = 0;
int       graphics_info_t::graphics_sample_step = 1;
int       graphics_info_t::dynamic_map_zoom_offset = 0;

// history
short int graphics_info_t::python_history = 1; // on
short int graphics_info_t::guile_history  = 1; // on
coot::history_list_t graphics_info_t::history_list;

// build one residue, n trials:
int graphics_info_t::add_terminal_residue_n_phi_psi_trials = 5000;
int graphics_info_t::add_terminal_residue_add_other_residue_flag = 0; // no.
std::string graphics_info_t::add_terminal_residue_type = "auto"; // was "ALA" before 20080601
short int graphics_info_t::add_terminal_residue_do_rigid_body_refine = 0; // off by default
bool      graphics_info_t::add_terminal_residue_debug_trials = false;
float graphics_info_t::rigid_body_fit_acceptable_fit_fraction = 0.75;

// rotamer
#ifdef USE_DUNBRACK_ROTAMERS
float graphics_info_t::rotamer_lowest_probability = 2.0; // percent
#else
float graphics_info_t::rotamer_lowest_probability = 0.0; // compatibility.  Limit
                                                         // not used, practically.
#endif
short int graphics_info_t::in_rotamer_define = 0;
int graphics_info_t::rotamer_residue_atom_index = -1; // unset initially.
int graphics_info_t::rotamer_residue_imol = -1;       // unset initially.
int graphics_info_t::rotamer_fit_clash_flag = 1;      //  check clashes initially.
short int graphics_info_t::in_auto_fit_define = 0;    // not in auto fit initially.
coot::atom_spec_t graphics_info_t::rotamer_residue_atom_spec;

// mutation
short int graphics_info_t::in_mutate_define = 0; // not, initially.
int graphics_info_t::mutate_residue_atom_index = -1;
int graphics_info_t::mutate_residue_imol = -1;
atom_selection_container_t asc;
atom_selection_container_t graphics_info_t::standard_residues_asc = asc;
std::string graphics_info_t::mutate_sequence_chain_from_combobox;
int         graphics_info_t::mutate_sequence_imol;
int         graphics_info_t::align_and_mutate_imol;
// std::string graphics_info_t::align_and_mutate_chain_from_optionmenu;
std::string graphics_info_t::align_and_mutate_chain_from_combobox;
int         graphics_info_t::nsv_canvas_pixel_limit = 22500;

mmdb::realtype    graphics_info_t::alignment_wgap   = -3.0; // was -0.5 (Bob) // was -3.0;
mmdb::realtype    graphics_info_t::alignment_wspace = -0.4;

//
short int graphics_info_t::in_reverse_direction_define = 0;

// user defined
short int graphics_info_t::in_user_defined_define = 0;
#ifdef USE_GUILE
SCM graphics_info_t::user_defined_click_scm_func = SCM_BOOL_F;
#endif // GUILE
#ifdef USE_PYTHON
PyObject *graphics_info_t::user_defined_click_py_func = NULL;
#endif // PYTHON
std::vector<coot::atom_spec_t> graphics_info_t::user_defined_atom_pick_specs;


// Miguel's axis orientation matrix ---------------
GL_matrix graphics_info_t::axes_orientation = GL_matrix();
short int graphics_info_t::use_axes_orientation_matrix_flag = 0;

// ---- NCS -----
float graphics_info_t::ncs_min_hit_ratio = 0.9;
short int graphics_info_t::ncs_maps_do_average_flag = 1;
short int graphics_info_t::ncs_matrix_flag = coot::NCS_LSQ;
// mutate auto fit
short int graphics_info_t::in_mutate_auto_fit_define = 0;
int graphics_info_t::mutate_auto_fit_residue_atom_index = -1;
int graphics_info_t::mutate_auto_fit_residue_imol = -1;
short int graphics_info_t::residue_type_chooser_auto_fit_flag = 0;
short int graphics_info_t::residue_type_chooser_stub_flag = 0;
short int graphics_info_t::mutate_auto_fit_do_post_refine_flag = 0;
short int graphics_info_t::rotamer_auto_fit_do_post_refine_flag = 0;

short int graphics_info_t::in_add_alt_conf_define = 0;
GtkWidget *graphics_info_t::add_alt_conf_dialog = NULL;
short int graphics_info_t::alt_conf_split_type = 1; // usually it is after Ca?
int graphics_info_t::add_alt_conf_atom_index = -1;
int graphics_info_t::add_alt_conf_imol = -1;
float graphics_info_t::add_alt_conf_new_atoms_occupancy = 0.5;
short int graphics_info_t::show_alt_conf_intermediate_atoms_flag = 0;
float graphics_info_t::ncs_homology_level = 0.7;

// dynarama
//
// edit phi/psi
short int graphics_info_t::in_edit_phi_psi_define = 0;
int graphics_info_t::edit_phi_psi_atom_index = -1;
int graphics_info_t::edit_phi_psi_imol = -1;
short int graphics_info_t::in_backbone_torsion_define = 0;
#ifdef DO_RAMA_PLOT
coot::rama_plot  *graphics_info_t::edit_phi_psi_plot = NULL;
int graphics_info_t::rama_psi_axis_mode = coot::rama_plot::PSI_CLASSIC;
#endif
float graphics_info_t::rama_level_prefered = 0.02;
float graphics_info_t::rama_level_allowed = 0.002;
float graphics_info_t::rama_plot_background_block_size = 2; // divisible into 360 preferably.
coot::ramachandran_points_container_t graphics_info_t::rama_points = coot::ramachandran_points_container_t();

ramachandrans_container_t graphics_info_t::ramachandrans_container = ramachandrans_container_t();

// edit chi
short int graphics_info_t::in_edit_chi_angles_define = 0;
int       graphics_info_t::edit_chi_current_chi = 0; // real values start at 1.
short int graphics_info_t::in_edit_chi_mode_flag = 0;
short int graphics_info_t::in_edit_chi_mode_view_rotate_mode = 0;
bool      graphics_info_t::edit_chi_angles_reverse_fragment = 0;
short int graphics_info_t::moving_atoms_move_chis_flag = 0;
coot::atom_spec_t graphics_info_t::chi_angles_clicked_atom_spec;
std::string graphics_info_t::chi_angle_alt_conf = "";
// multi-residue torsion
bool graphics_info_t::in_multi_residue_torsion_define = false;
bool graphics_info_t::in_multi_residue_torsion_mode = false;
bool graphics_info_t::multi_residue_torsion_reverse_fragment_mode = false;
int  graphics_info_t::multi_residue_torsion_picked_residues_imol = -1;
std::pair<int, int> graphics_info_t::multi_residue_torsion_rotating_atom_index_pair = std::pair<int,int>(-1,-1);
std::vector<coot::residue_spec_t> graphics_info_t::multi_residue_torsion_picked_residue_specs;


// 180 degree flip
short int graphics_info_t::in_180_degree_flip_define = 0;

int graphics_info_t::ramachandran_plot_differences_imol1 = -1;
int graphics_info_t::ramachandran_plot_differences_imol2 = -1;
std::string graphics_info_t::ramachandran_plot_differences_imol1_chain = "";
std::string graphics_info_t::ramachandran_plot_differences_imol2_chain = "";


// state language variable:
short int graphics_info_t::state_language = 1; // scheme

// directory saving:
std::string graphics_info_t::directory_for_fileselection = "";
std::string graphics_info_t::directory_for_saving_for_fileselection = "";
std::string graphics_info_t::directory_for_filechooser = "";
std::string graphics_info_t::directory_for_saving_for_filechooser = "";


// Residue info
short int graphics_info_t::residue_info_pending_edit_b_factor = 0;
short int graphics_info_t::residue_info_pending_edit_occ = 0;
int       graphics_info_t::residue_info_n_atoms = -1;
std::vector<coot::select_atom_info> graphics_info_t::residue_info_edits; // used to be static (written in the very olden days)

//  Backbone torsions:
clipper::Coord_orth graphics_info_t::backbone_torsion_end_ca_1;
clipper::Coord_orth graphics_info_t::backbone_torsion_end_ca_2;
int graphics_info_t::backbone_torsion_peptide_button_start_pos_x;
int graphics_info_t::backbone_torsion_peptide_button_start_pos_y;
int graphics_info_t::backbone_torsion_carbonyl_button_start_pos_x;
int graphics_info_t::backbone_torsion_carbonyl_button_start_pos_y;

// show citation notice?
// short int graphics_info_t::show_citation_notice = 1; // on by default :)
short int graphics_info_t::show_citation_notice = 0; // on by default :)

// we have dragged shear fixed points?
short int graphics_info_t::have_fixed_points_sheared_drag_flag = 0;
// smaller is smoother and less jerky - especially for big molecules
int       graphics_info_t::dragged_refinement_steps_per_frame = 20;
short int graphics_info_t::dragged_refinement_refine_per_frame_flag = 0;
double    graphics_info_t::refinement_drag_elasticity = 0.25;

// save the restraints:
//
coot::restraints_container_t *graphics_info_t::last_restraints = 0;
// 20220504-PE so that I can check for cleared/removed non-bonded contact baddies
std::map<int, std::vector<int> > graphics_info_t::previous_round_nbc_baddies_atom_index_map;

// clipper::Xmap<float> blank_dummy_xmap;
// ref version: coot::restraints_container_t(blank_dummy_xmap);
//
//
bool graphics_info_t::draw_zero_occ_spots_flag = true; // on by default

bool graphics_info_t::draw_cis_peptide_markups = true; // on by default

int   graphics_info_t::check_waters_molecule = -1; // unset initially.
float graphics_info_t::check_waters_b_factor_limit  = 80.0;
float graphics_info_t::check_waters_map_sigma_limit = 1.0;
float graphics_info_t::check_waters_min_dist_limit = 2.3;
float graphics_info_t::check_waters_max_dist_limit = 3.5;
float graphics_info_t::check_waters_by_difference_map_sigma_level = 3.5;
int   graphics_info_t::check_waters_by_difference_map_map_number = -1;

// default map sigma level:
float graphics_info_t::default_sigma_level_for_map = 1.5;
float graphics_info_t::default_sigma_level_for_fofc_map = 3.0;


// geometry widget:
GtkWidget *graphics_info_t::geometry_dialog = NULL;

bool graphics_info_t::add_ccp4i_projects_to_optionmenu_flag = true;

// run refmac widget:
std::string graphics_info_t::refmac_ccp4i_project_dir = std::string("");
std::string graphics_info_t::libcheck_ccp4i_project_dir = std::string("");

std::vector<int> *graphics_info_t::preset_number_refmac_cycles;
coot::refmac::refmac_refinement_method_type graphics_info_t::refmac_refinement_method = coot::refmac::RESTRAINED;
coot::refmac::refmac_phase_input_type graphics_info_t::refmac_phase_input   = coot::refmac::NO_PHASES;
coot::refmac::refmac_use_tls_type     graphics_info_t::refmac_use_tls_flag  = coot::refmac::TLS_ON;
coot::refmac::refmac_use_twin_type    graphics_info_t::refmac_use_twin_flag = coot::refmac::TWIN_OFF;
coot::refmac::refmac_use_sad_type     graphics_info_t::refmac_use_sad_flag  = coot::refmac::SAD_OFF;
coot::refmac::refmac_use_ncs_type     graphics_info_t::refmac_use_ncs_flag  = coot::refmac::NCS_ON;
coot::refmac::refmac_use_intensities_type graphics_info_t::refmac_use_intensities_flag  = coot::refmac::AMPLITUDES;
coot::refmac::refmac_used_mtz_file_type graphics_info_t::refmac_used_mtz_file_flag = coot::refmac::MTZ;
const gchar *graphics_info_t::saved_refmac_file_filename = NULL;
int graphics_info_t::refmac_ncycles = 5;
GtkWidget *graphics_info_t::refmac_dialog_mtz_file_label = NULL;
std::vector<coot::refmac::sad_atom_info_t> graphics_info_t::refmac_sad_atoms;
short int graphics_info_t::have_sensible_refmac_params = 0;
std::string graphics_info_t::refmac_mtz_file_filename = "";
std::string graphics_info_t::refmac_fobs_col = "";
std::string graphics_info_t::refmac_sigfobs_col = "";
std::string graphics_info_t::refmac_r_free_col = "";
int graphics_info_t::refmac_r_free_flag_sensible = 0;

// scrollin' scrollin' scrollin'... Shall we stop? When shall we stop?
short int graphics_info_t::stop_scroll_diff_map_flag = 1; // stop on
short int graphics_info_t::stop_scroll_iso_map_flag  = 1; // ditto.
float     graphics_info_t::stop_scroll_diff_map_level = 0.0;
float     graphics_info_t::stop_scroll_iso_map_level = 0.0;

// globing
//
std::vector<std::string> *graphics_info_t::coordinates_glob_extensions;
std::vector<std::string> *graphics_info_t::data_glob_extensions;
std::vector<std::string> *graphics_info_t::map_glob_extensions;
std::vector<std::string> *graphics_info_t::dictionary_glob_extensions;

// superposition
int graphics_info_t::superpose_imol1 = -1;
int graphics_info_t::superpose_imol2 = -1;
std::string graphics_info_t::superpose_imol1_chain = "";
std::string graphics_info_t::superpose_imol2_chain = "";

// unbonded star size [doesn't work yet]
float graphics_info_t::unbonded_atom_star_size = 0.5;

// Raster3D
float graphics_info_t::raster3d_bond_thickness    = 0.18;
float graphics_info_t::raster3d_atom_radius    = 0.25;
float graphics_info_t::raster3d_density_thickness = 0.015;
bool  graphics_info_t::raster3d_enable_shadows = 1;
int   graphics_info_t::renderer_show_atoms_flag = 1;
float graphics_info_t::raster3d_bone_thickness    = 0.05;
int   graphics_info_t::raster3d_water_sphere_flag = 0;
std::string graphics_info_t::raster3d_font_size = "4";

// map (density) line thickness:
int graphics_info_t::map_line_width = 1;

// bonding stuff
int   graphics_info_t::bond_thickness_intermediate_value = -1;
float graphics_info_t::bond_thickness_intermediate_atoms = 5; // (no so) thick white atom bonds

// merge molecules
int graphics_info_t::merge_molecules_master_molecule = -1;
std::vector<int> *graphics_info_t::merge_molecules_merging_molecules;
coot::residue_spec_t graphics_info_t::merge_molecules_ligand_spec;

// change chain ids:
int graphics_info_t::change_chain_id_molecule = -1;
std::string graphics_info_t::change_chain_id_from_chain = "";

// renumber residues
int         graphics_info_t::renumber_residue_range_molecule = -1;
std::string graphics_info_t::renumber_residue_range_chain;

// antialiasing:
short int graphics_info_t::do_anti_aliasing_flag = 0;

// lighting
short int graphics_info_t::do_lighting_flag = 0;
bool      graphics_info_t::do_flat_shading_for_solid_density_surface = 1;


// stereo?
int graphics_info_t::display_mode = coot::MONO_MODE;
float graphics_info_t::hardware_stereo_angle_factor = 1.0;

// remote controlled coot
int graphics_info_t::try_port_listener = 0;
int graphics_info_t::remote_control_port_number;
std::string graphics_info_t::remote_control_hostname;
int graphics_info_t::coot_socket_listener_idle_function_token = -1; //  default (off)
// Did we get a good socket when we tried to open it?  If so, set
// something non-zero here (which is done as a scheme command).
int graphics_info_t::listener_socket_have_good_socket_state = 0;

// I don't think that we need the mutex stuff when using waiting strings
// - so python version doesn't have them (at the moment).
std::string graphics_info_t::socket_string_waiting = "";
std::string graphics_info_t::socket_python_string_waiting = "";
volatile bool graphics_info_t::have_socket_string_waiting_flag = false;
volatile bool graphics_info_t::have_socket_python_string_waiting_flag = false;
volatile bool graphics_info_t::socket_string_waiting_mutex_lock = false;


// validation
std::vector<clipper::Coord_orth> *graphics_info_t::diff_map_peaks = new std::vector<clipper::Coord_orth>;
int   graphics_info_t::max_diff_map_peaks = 0;
float graphics_info_t::difference_map_peaks_sigma_level = 5.6; // 20220419-PE 5.0 give too many peaks
float graphics_info_t::difference_map_peaks_max_closeness = 2.0; // A

// save state file name
#ifdef USE_GUILE
std::string graphics_info_t::save_state_file_name = "0-coot.state.scm";
#else
std::string graphics_info_t::save_state_file_name = "0-coot.state.py";
#endif

// auto-read mtz columns
// std::string graphics_info_t::auto_read_MTZ_FWT_col = "FWT";
// std::string graphics_info_t::auto_read_MTZ_PHWT_col = "PHWT";
// std::string graphics_info_t::auto_read_MTZ_DELFWT_col = "DELFWT";
// std::string graphics_info_t::auto_read_MTZ_PHDELWT_col = "PHDELWT";

std::vector<coot::mtz_column_trials_info_t> graphics_info_t::user_defined_auto_mtz_pairs;

short int graphics_info_t::probe_available = -1; // don't know yet

// fffearing
float graphics_info_t::fffear_angular_resolution = 15.0; // degrees

// move molecule here
int graphics_info_t::move_molecule_here_molecule_number = -1;

// pseudo bond for sec str restraints
coot::pseudo_restraint_bond_type graphics_info_t::pseudo_bonds_type = coot::NO_PSEUDO_BONDS;


// MYSQL database
#ifdef USE_MYSQL_DATABASE
MYSQL *graphics_info_t::mysql = 0;
int    graphics_info_t::query_number = 1;
std::string graphics_info_t::sessionid = "";
std::pair<std::string, std::string> graphics_info_t::db_userid_username("no-userid","no-user-name");
std::string graphics_info_t::mysql_host   = "localhost";
std::string graphics_info_t::mysql_user   = "cootuser";
std::string graphics_info_t::mysql_passwd = "password";
#endif // USE_MYSQL_DATABASE

//
int graphics_info_t::ncs_next_chain_skip_key = GDK_KEY_o;
int graphics_info_t::ncs_prev_chain_skip_key = GDK_KEY_O;
int graphics_info_t::update_go_to_atom_from_current_residue_key = GDK_KEY_p;

//
#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
// 20220528-PE FIXME
#else
GdkCursorType graphics_info_t::pick_cursor_index = GDK_CROSSHAIR;
#endif

// update self?
bool graphics_info_t::update_self = 0;  // Set by command line arg --update-self

float graphics_info_t::electrostatic_surface_charge_range = 0.5;

// network
volatile bool graphics_info_t::curl_handlers_lock = 0; // not locked.


// Defaults for the file chooser
int graphics_info_t::gtk2_file_chooser_selector_flag = coot::CHOOSER_STYLE;
int graphics_info_t::gtk2_chooser_overwrite_flag = coot::CHOOSER_OVERWRITE_PROTECT;


// Graphics ligand view
// graphics_ligand_molecule graphics_info_t::graphics_ligand_mol;
int                      graphics_info_t::show_graphics_ligand_view_flag = 1; // user control

// FLEV
// use-defined flev params
float graphics_info_t::fle_water_dist_max = 3.25;   // 3.25
float graphics_info_t::fle_h_bond_dist_max = 3.9;   // 3.9

// Mogul
float graphics_info_t::mogul_max_badness = 5.0;   // The z value colour at which the bond is fully red.
                                                  // Some might say that 5.0 is to liberal (it allow
                                                  // too much badness, that is).

// Glyco fit and refine
bool graphics_info_t::linked_residue_fit_and_refine_state = true;

//
bool graphics_info_t::allow_duplseqnum = true; // 20181214-PE - I presume that this is safe now?

std::map<std::string, std::string> graphics_info_t::extensions_registry;

//
std::map<std::string, std::pair<std::string, std::string> > graphics_info_t::user_name_passwd_map;

std::vector<std::pair<clipper::Coord_orth, std::string> > graphics_info_t::user_defined_interesting_positions;
unsigned int graphics_info_t::user_defined_interesting_positions_idx = 0;

// atom pull
// atom_pull_info_t graphics_info_t:: atom_pull = atom_pull_info_t();
std::vector<atom_pull_info_t> graphics_info_t:: atom_pulls;
bool graphics_info_t::auto_clear_atom_pull_restraint_flag = true;

bool graphics_info_t::continue_update_refinement_atoms_flag = false;

std::pair<bool, float> graphics_info_t::model_display_radius = std::pair<bool, float> (false, 15);

// need to configure for this!
// #define GDKGLEXT_HAVE_MODE_SAMPLES_SHIFT true

// Chemical Feature Clusters, cfc
GtkWidget *graphics_info_t::cfc_dialog = NULL;

cfc_gui_t graphics_info_t::cfc_gui;

bool graphics_info_t::coot_is_a_python_module = true;

bool graphics_info_t::residue_type_selection_was_user_picked_residue_range = false;

bool graphics_info_t::make_auto_h_bond_restraints_flag = false;

bool graphics_info_t::do_rotamer_restraints = false;

bool graphics_info_t::do_debug_refinement = false;

std::atomic<bool> graphics_info_t::on_going_updating_map_lock(false);

std::map<keyboard_key_t, key_bindings_t> graphics_info_t::key_bindings_map;
std::vector<keyboard_key_t> graphics_info_t::keyboard_key_history;

std::string graphics_info_t::mtz_file_for_refmac;

bool graphics_info_t::convert_dictionary_planes_to_improper_dihedrals_flag = false;

GtkBuilder *graphics_info_t::gtkbuilder  = NULL;
GtkBuilder *graphics_info_t::preferences_gtkbuilder  = NULL;
GtkWidget  *graphics_info_t::main_window = NULL;

// now the clipping planes are scale, not offsets
float graphics_info_t::clipping_front = 1.0;
float graphics_info_t::clipping_back  = 1.0;

#ifdef USE_MOLECULES_TO_TRIANGLES
std::shared_ptr<Renderer>   graphics_info_t::mol_tri_renderer    = 0;
std::shared_ptr<SceneSetup> graphics_info_t::mol_tri_scene_setup = 0;
#endif // USE_MOLECULES_TO_TRIANGLES

graphics_ligand_mesh_molecule_t graphics_info_t::graphics_ligand_mesh_molecule;

// --------------------------------------------------------------------------------------------

float *graphics_info_t::mvp = new float[16];
int    graphics_info_t::mvp_location = -1;
int    graphics_info_t::view_rotation_location = -1;
glm::quat graphics_info_t::view_quaternion = glm::quat(1,0,0,0);
GLuint graphics_info_t::programID_for_central_cube = 0;
GLuint graphics_info_t::central_cube_vertexarray_id = 0;
GLuint graphics_info_t::central_cube_array_buffer_id = 0;
GLuint graphics_info_t::central_cube_index_buffer_id = 0;
GLuint graphics_info_t::hud_text_vertexarray_id = 0;
GLuint graphics_info_t::hud_text_array_buffer_id = 0;
GLuint graphics_info_t::screen_quad_vertex_array_id = 0;
GLuint graphics_info_t::blur_x_quad_vertex_array_id = 0;
GLuint graphics_info_t::blur_y_quad_vertex_array_id = 0;
GLuint graphics_info_t::combine_textures_using_depth_quad_vertex_array_id = 0;
GLuint graphics_info_t::blur_quad_vertex_array_id = 0;
GLuint graphics_info_t::textureColorbuffer_screen = 0;
GLuint graphics_info_t::textureColorbuffer_blur = 0;

GLuint graphics_info_t::rotation_centre_crosshairs_vertexarray_id = 0;
GLuint graphics_info_t::rotation_centre_crosshairs_vertex_buffer_id = 0;
GLuint graphics_info_t::rotation_centre_crosshairs_index_buffer_id = 0;

bool graphics_info_t::use_framebuffers = true;
framebuffer graphics_info_t::screen_framebuffer;
framebuffer graphics_info_t::blur_x_framebuffer;
framebuffer graphics_info_t::blur_y_framebuffer;
framebuffer graphics_info_t::combine_textures_using_depth_framebuffer;
framebuffer graphics_info_t::blur_framebuffer; // 2020
unsigned int graphics_info_t::framebuffer_scale = 1; // on supersampling by default.

bool graphics_info_t::perspective_projection_flag = false;
float graphics_info_t::perspective_fov = 26.0; // was 30.0


// --------------------------------------------------------------------------------------------

// GLuint graphics_info_t::programID_for_maps = 0; in a shader now  - as
//programID_for_central_cube should be
Shader graphics_info_t::shader_for_maps;
Shader graphics_info_t::shader_for_map_caps;
Shader graphics_info_t::shader_for_models;
Shader graphics_info_t::shader_for_model_as_meshes;
Shader graphics_info_t::shader_for_moleculestotriangles;
Shader graphics_info_t::shader_for_symmetry_atoms_bond_lines;
Shader graphics_info_t::shader_for_central_cube;
Shader graphics_info_t::shader_for_origin_cube;
Shader graphics_info_t::shader_for_hud_text;
Shader graphics_info_t::shader_for_hud_geometry_bars;
Shader graphics_info_t::shader_for_hud_geometry_labels;
Shader graphics_info_t::shader_for_hud_buttons;
Shader graphics_info_t::shader_for_hud_image_texture;
Shader graphics_info_t::shader_for_atom_labels;
Shader graphics_info_t::shader_for_rama_balls;
Shader graphics_info_t::shader_for_x_blur;
Shader graphics_info_t::shader_for_y_blur;
Shader graphics_info_t::shader_for_dof_blur_by_texture_combination;
Shader graphics_info_t::shader_for_effects;
Shader graphics_info_t::shader_for_hud_lines;
Shader graphics_info_t::shader_for_lines;
// Shader graphics_info_t::shader_for_anchored_atom_markers; // use happy face markers instead
Shader graphics_info_t::shader_for_lines_pulse;
Shader graphics_info_t::shader_for_particles;
Shader graphics_info_t::shader_for_ligand_view;
Shader graphics_info_t::shader_for_instanced_objects; // used for boids - also HOLE
Shader graphics_info_t::shader_for_extra_distance_restraints;
Shader graphics_info_t::shader_for_hud_geometry_tooltip_text;
Shader graphics_info_t::shader_for_happy_face_residue_markers;
Shader graphics_info_t::shader_for_happy_face_residue_markers_for_ssao;
Shader graphics_info_t::shader_for_rama_plot_axes_and_ticks;
Shader graphics_info_t::shader_for_rama_plot_phi_phis_markers;
bool   graphics_info_t::shaders_have_been_compiled = false;
meshed_generic_display_object graphics_info_t::mesh_for_environment_distances;
std::chrono::time_point<std::chrono::high_resolution_clock> graphics_info_t::previous_frame_time = std::chrono::high_resolution_clock::now();
std::chrono::time_point<std::chrono::high_resolution_clock> graphics_info_t::previous_frame_time_for_per_second_counter = std::chrono::high_resolution_clock::now();
long graphics_info_t::frame_counter = 0;
long graphics_info_t::frame_counter_at_last_display = 0;
float graphics_info_t::fps = 0.0;
float graphics_info_t::fps_std_dev = -1.0;
float graphics_info_t::fps_times_scale_factor = 0.002;
std::list<std::chrono::time_point<std::chrono::high_resolution_clock> > graphics_info_t::frame_time_history_list;

std::set<mmdb::Residue *> graphics_info_t::moving_atoms_visited_residues;
mmdb::Atom *graphics_info_t::active_atom_for_hud_geometry_bar = 0;

unsigned int graphics_info_t::draw_count_for_happy_face_residue_markers = 0;

glm::vec3 graphics_info_t::eye_position = glm::vec3(0,0,95);
float graphics_info_t::screen_z_near_perspective =  76.0; // was 83
float graphics_info_t::screen_z_far_perspective  = 125.0;

float graphics_info_t::goodselliness = 0.3; // the pastelization factor

std::map<unsigned int, lights_info_t> graphics_info_t::lights;

std::vector<molecule_class_info_t> graphics_info_t::molecules;
molecule_class_info_t graphics_info_t::moving_atoms_molecule;
std::atomic<bool> molecule_class_info_t::draw_vector_sets_lock(false);

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
coot::Cartesian graphics_info_t::smooth_scroll_target_point = coot::Cartesian(0,0,0);
coot::Cartesian graphics_info_t::smooth_scroll_start_point = coot::Cartesian(0,0,0);

bool graphics_info_t::shader_do_ambient_occlusion_flag = true;
bool graphics_info_t::shader_do_depth_blur_flag = true; // 2020 version
bool graphics_info_t::shader_do_depth_fog_flag = true;
bool graphics_info_t::shader_do_outline_flag = false;
bool graphics_info_t::shader_do_depth_of_field_blur_flag = false;
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


bool graphics_info_t::draw_missing_loops_flag = true;

bool graphics_info_t::sequence_view_is_docked_flag = true;


int graphics_info_t::tick_function_id = -1; // unset
bool graphics_info_t::do_tick_particles = false;
bool graphics_info_t::do_tick_spin = false;
bool graphics_info_t::do_tick_rock = false;
bool graphics_info_t::do_tick_boids = false;
bool graphics_info_t::do_tick_constant_draw = false; // was true: 20220606-PE hack because keyboard controller not working yet
bool graphics_info_t::do_tick_hydrogen_bonds_mesh = false;
bool graphics_info_t::do_tick_happy_face_residue_markers = false;
bool graphics_info_t::do_tick_outline_for_active_residue = false;
bool graphics_info_t::do_tick_gone_diegos = false;
bool graphics_info_t::do_tick_gone_diff_map_peaks = false;
int graphics_info_t::n_particles = 220;
Mesh graphics_info_t::mesh_for_particles = Mesh("mesh for particles");
particle_container_t graphics_info_t::particles;
bool graphics_info_t::setup_draw_for_particles_semaphore = false;
glm::vec3 graphics_info_t::identification_pulse_centre;
bool graphics_info_t::particles_have_been_shown_already_for_this_round_flag = false;
int graphics_info_t::wait_for_hooray_refinement_tick_id = -1; // delete this tick function on refinement
                                                              // shutdown

std::vector<std::pair<glm::vec3, glm::vec3> > graphics_info_t::hydrogen_bonds_atom_position_pairs;

std::chrono::time_point<std::chrono::high_resolution_clock> graphics_info_t::tick_hydrogen_bond_mesh_t_previous = std::chrono::high_resolution_clock::now();

Mesh graphics_info_t::mesh_for_outline_of_active_residue = Mesh("mesh for active residue");
Shader graphics_info_t::shader_for_outline_of_active_residue;
unsigned int graphics_info_t::outline_for_active_residue_frame_count = 0;

fun::boids_container_t graphics_info_t::boids;
Mesh graphics_info_t::mesh_for_boids = Mesh("mesh for boids");
LinesMesh graphics_info_t::lines_mesh_for_boids_box;

Mesh graphics_info_t::mesh_for_hydrogen_bonds = Mesh("mesh for hydrogen bonds");

LinesMesh graphics_info_t::lines_mesh_for_identification_pulse;
LinesMesh graphics_info_t::lines_mesh_for_generic_pulse;
std::vector<glm::vec3> graphics_info_t::generic_pulse_centres;

LinesMesh graphics_info_t::lines_mesh_for_hud_lines;
LinesMesh graphics_info_t::lines_mesh_for_pull_restraint_neighbour_displacement_max_radius_ring;

std::vector<atom_label_info_t> graphics_info_t::labels;
TextureMesh graphics_info_t::tmesh_for_labels = TextureMesh("tmesh-for-labels");
HUDMesh graphics_info_t::mesh_for_hud_geometry = HUDMesh("hud-geometry");
// these 3 for testing images (20210831-PE may be usefule for mip mappsing testing later)
Texture graphics_info_t::texture_for_camera_facing_quad;
TextureMesh graphics_info_t::tmesh_for_camera_facing_quad;
Shader graphics_info_t::camera_facing_quad_shader;
HUDTextureMesh graphics_info_t::tmesh_for_hud_image_testing = HUDTextureMesh("tmesh_for_hud_image_testing");

bool graphics_info_t::draw_bad_nbc_atom_pair_markers_flag = false;
TextureMesh graphics_info_t::tmesh_for_happy_face_residues_markers = TextureMesh("tmesh-for-happy-faces");
Texture graphics_info_t::texture_for_happy_face_residue_marker;
std::vector<glm::vec3> graphics_info_t::happy_face_residue_marker_starting_positions;

TextureMesh graphics_info_t::tmesh_for_bad_nbc_atom_pair_markers = TextureMesh("tmesh-for-angry-diego");
Texture graphics_info_t::texture_for_bad_nbc_atom_pair_markers;
std::vector<glm::vec3> graphics_info_t::bad_nbc_atom_pair_marker_positions;

Mesh graphics_info_t::bad_nbc_atom_pair_dashed_line = Mesh("bad nbc atom_pair dashed line instanced mesh");

TextureMesh graphics_info_t::tmesh_for_unhappy_atom_markers = TextureMesh("tmesh-unhappy-atom-outliers");
Texture graphics_info_t::texture_for_unhappy_atom_markers;

TextureMesh graphics_info_t::tmesh_for_chiral_volume_outlier_markers = TextureMesh("tmesh-chiral-volume-outliers");
Texture graphics_info_t::texture_for_chiral_volume_outlier_markers;

TextureMesh graphics_info_t::tmesh_for_anchored_atom_markers = TextureMesh("tmesh-for-anchored-atoms");
Texture graphics_info_t::texture_for_anchored_atom_markers;
std::vector<glm::vec3> graphics_info_t::anchored_atom_marker_texture_positions;

HUDTextureMesh graphics_info_t::tmesh_for_hud_geometry_tooltip_label = HUDTextureMesh("tmesh-for-hud-geometry-tooltip-labels");

HUDTextureMesh graphics_info_t::tmesh_for_hud_refinement_dialog_arrow = HUDTextureMesh("tmesh-for-hud-refinement-dialog-arrow");
Texture graphics_info_t::texture_for_hud_refinement_dialog_arrow;
Texture graphics_info_t::texture_for_hud_refinement_dialog_arrow_highlighted;
bool graphics_info_t::hud_refinement_dialog_arrow_is_moused_over = false;

Texture        graphics_info_t::texture_for_background_image;
HUDTextureMesh graphics_info_t::tmesh_for_background_image = HUDTextureMesh("tmesh-for-background-image");
Shader         graphics_info_t::shader_for_background_image;
bool           graphics_info_t::draw_background_image_flag = false; // uses "background-image.png"

float graphics_info_t::pull_restraint_neighbour_displacement_max_radius = 0.0; // we don't see it initially.

xdg_t xdg;
coot::command_history_t graphics_info_t::command_history = coot::command_history_t(xdg);

std::vector<Instanced_Markup_Mesh> graphics_info_t::instanced_meshes;

std::vector<TextureMesh> graphics_info_t::texture_meshes;
Shader graphics_info_t::shader_for_texture_meshes;

std::map<std::string, Texture> graphics_info_t::texture_for_hud_geometry_labels_map;
Texture graphics_info_t::texture_for_hud_tooltip_background;
bool graphics_info_t::draw_hud_tooltip_flag = false;

HUDTextureMesh graphics_info_t::mesh_for_hud_geometry_labels = HUDTextureMesh("tmesh-for-hud-geometry-labels");
HUDTextureMesh graphics_info_t::mesh_for_hud_tooltip_background = HUDTextureMesh("tmesh-for-hud-tooltip-background");

Instanced_Markup_Mesh graphics_info_t::rama_balls_mesh = Instanced_Markup_Mesh("rama-balls");
bool graphics_info_t::draw_stick_mode_atoms_default = true;

std::vector<HUD_button_info_t> graphics_info_t::hud_button_info;
HUDMesh graphics_info_t::mesh_for_hud_buttons = HUDMesh("mesh-for-hud-buttons");

std::string graphics_info_t::label_for_hud_geometry_tooltip;

bool graphics_info_t::auto_recontour_map_flag = true;

bool graphics_info_t::mol_displayed_toggle_do_redraw = true; // normally true

double graphics_info_t::torsion_restraints_weight = 1.0;

bool graphics_info_t::use_harmonic_approximation_for_NBCs = false;

bool graphics_info_t::draw_hud_colour_bar_flag = false;
std::vector<std::pair<unsigned int, coot::colour_holder> > graphics_info_t::user_defined_colours; // initially empty

unsigned int graphics_info_t::bond_smoothness_factor = 1; // changes num_subdivisions and n_slices

float graphics_info_t::contact_dots_density = 0.6; // 20220308-PE was 1.0, 20230613-PE was 0.4 (too low for ligand contact dots)
float graphics_info_t::contact_dot_sphere_subdivisions = 1;
bool graphics_info_t::all_atom_contact_dots_ignore_water_flag = false;
bool graphics_info_t::all_atom_contact_dots_do_vdw_surface = false;

bool graphics_info_t::refinement_has_finished_moving_atoms_representation_update_needed_flag = false;

gl_rama_plot_t graphics_info_t::gl_rama_plot;
bool graphics_info_t::draw_gl_ramachandran_plot_flag = false;
bool graphics_info_t::draw_gl_ramachandran_plot_user_control_flag = true;


float graphics_info_t::focus_blur_z_depth = 0.15;
float graphics_info_t::focus_blur_strength = 1.0;

// 20220129-PE crows

std::vector<Model> graphics_info_t::models;

unsigned int graphics_info_t::noiseTexture = 0;
unsigned int graphics_info_t::ssaoColorBuffer = 0;
unsigned int graphics_info_t::ssaoColorBufferBlur = 0;
std::vector<glm::vec3> graphics_info_t::ssaoKernel;
unsigned int graphics_info_t::ssao_blur_size = 1;

Shader graphics_info_t::shader_for_meshes;
Shader graphics_info_t::shader_for_tmeshes;
Shader graphics_info_t::shader_for_meshes_shadow_map;
Shader graphics_info_t::shader_for_instanced_meshes_shadow_map;
Shader graphics_info_t::shader_for_texture_meshes_shadow_map;
Shader graphics_info_t::shader_for_shadow_map_image_texture_mesh;
float graphics_info_t::shadow_box_size = 120.0;

float graphics_info_t::SSAO_bias = 0.02;
float graphics_info_t::SSAO_radius = 30.0;
float graphics_info_t::ssao_strength = 0.4;
bool graphics_info_t::use_ssao  = true;  // in the effects filter, adds (or not) the SSAO effects
bool graphics_info_t::show_just_ssao = false; // in the effects filter, shows *just* the SSAO effects
unsigned int graphics_info_t::effects_shader_output_type(EFFECTS_SHADER_STANDARD);

float graphics_info_t::effects_brightness = 1.0f;
float graphics_info_t::effects_gamma = 1.0f;

Shader graphics_info_t::shader_for_tmeshes_for_ssao;
Shader graphics_info_t::shader_for_meshes_for_ssao;
Shader graphics_info_t::shader_for_instanced_meshes_for_ssao;
Shader graphics_info_t::shader_for_rotation_centre_cross_hairs_for_ssao;

unsigned int graphics_info_t::rboDepth = 0;
unsigned int graphics_info_t::ssaoFBO = 0;
unsigned int graphics_info_t::ssaoBlurFBO = 0;

framebuffer graphics_info_t::framebuffer_for_ssao_gbuffer;
framebuffer graphics_info_t::framebuffer_for_ssao;
framebuffer graphics_info_t::framebuffer_for_ssao_blur;

unsigned int graphics_info_t::n_ssao_kernel_samples = 64;
Shader graphics_info_t::shaderGeometryPass;
Shader graphics_info_t::shaderSSAO;
Shader graphics_info_t::shaderSSAOBlur;

framebuffer graphics_info_t::framebuffer_for_effects;

unsigned int graphics_info_t::shadow_depthMap_framebuffer = 0;
unsigned int graphics_info_t::shadow_depthMap_texture = 0; // the texture
float graphics_info_t::shadow_strength = 0.0; // 0 to 1 (strong shadows don't look good)
unsigned int graphics_info_t::shadow_softness = 2; // 1, 2 or 3
unsigned int graphics_info_t::shadow_texture_multiplier = 2;
unsigned int graphics_info_t::shadow_texture_width  = graphics_info_t::shadow_texture_multiplier * 1024;
unsigned int graphics_info_t::shadow_texture_height = graphics_info_t::shadow_texture_multiplier * 1024;
bool graphics_info_t::show_just_shadows = false; // show *just* the shadows in the texture-mesh-with-shadows shader

// unsigned short int graphics_info_t::displayed_image_type = graphics_info_t::SHOW_AO_SCENE;
unsigned short int graphics_info_t::displayed_image_type = graphics_info_t::SHOW_BASIC_SCENE;

GLuint graphics_info_t::screen_AO_quad_vertex_array_id = 0;
GLuint graphics_info_t::screen_AO_quad_VBO = 0;

Shader graphics_info_t::shader_for_tmeshes_with_shadows;
Shader graphics_info_t::shader_for_meshes_with_shadows;
Shader graphics_info_t::shader_for_instanced_meshes_with_shadows;
HUDTextureMesh graphics_info_t::tmesh_for_shadow_map = HUDTextureMesh("tmesh-for-shadow-map");

bool graphics_info_t::stereo_style_2010 = false;

bool graphics_info_t::ignore_pseudo_zeros_for_map_stats = true;

HUDTextureMesh graphics_info_t::tmesh_for_hud_colour_bar = HUDTextureMesh("tmesh for HUD colour bar");
Texture graphics_info_t::texture_for_hud_colour_bar;
GtkListStore* graphics_info_t::validation_graph_model_list = nullptr;
int graphics_info_t::active_validation_graph_model_idx = -1;
graphics_info_t::validation_graph_map_t graphics_info_t::validation_graph_widgets = graphics_info_t::validation_graph_map_t();
graphics_info_t::validation_data_map_t graphics_info_t::validation_graph_data = graphics_info_t::validation_data_map_t();
std::string graphics_info_t::active_validation_graph_chain_id = std::string();

GtkListStore* graphics_info_t::ramachandran_plot_model_list = nullptr;
std::vector<graphics_info_t::widgeted_rama_plot_t> graphics_info_t::rama_plot_boxes;

// 20230430-PE updating maps
int graphics_info_t::updating_maps_imol_map      = -1;
int graphics_info_t::updating_maps_imol_diff_map = -1;
std::vector<api::rail_points_t> graphics_info_t::rail_point_history;
coot::util::sfcalc_genmap_stats_t graphics_info_t::latest_sfcalc_stats;

bool graphics_info_t::curmudgeon_mode = false;
bool graphics_info_t::use_sounds = true;
guint graphics_info_t::updating_maps_timeout_function_idx = UPDATING_MAPS_TIMEOUT_FUNCTION_IDX_UNSET;

std::vector<meshed_particle_container_t> graphics_info_t::meshed_particles_for_gone_diegos;
meshed_particle_container_t graphics_info_t::meshed_particles_for_gone_diff_map_peaks(Mesh("gone diff map peaks"), particle_container_t());

float graphics_info_t::gaussian_surface_sigma = 4.4;
float graphics_info_t::gaussian_surface_contour_level = 4.0;
float graphics_info_t::gaussian_surface_box_radius = 5.0;
float graphics_info_t::gaussian_surface_grid_scale = 0.7;
float graphics_info_t::gaussian_surface_fft_b_factor = 100.0;
short int graphics_info_t::gaussian_surface_chain_colour_mode = 1; // 1 for "by chain" , 2 for "by NCS"

std::vector<coot::positron_metadata_t> graphics_info_t::positron_metadata;

bool graphics_info_t::tomo_picker_flag = false;
graphics_info_t::tomo_view_info_t graphics_info_t::tomo_view_info;

coot::inchikey_store_t graphics_info_t::inchikey_store;

std::pair<bool, std::string> graphics_info_t::servalcat_fofc    = std::pair<bool, std::string> (false, "");
std::pair<bool, std::string> graphics_info_t::servalcat_refine  = std::pair<bool, std::string> (false, "");

std::string graphics_info_t::current_alt_conf = "";
