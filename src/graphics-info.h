/* src/graphics-info.h
 * -*-c++-*-
 *
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2007 by Paul Emsley
 * Copyright 2007, 2008 by The University of Oxford
 * Copyright 2016 by Medical Research Council
 *
 * Author: Paul Emsley
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
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

// -*-c++-*-
#ifndef GRAPHICS_INFO_H
#define GRAPHICS_INFO_H

#include "compat/coot-sysdep.h"
#include "geometry/residue-and-atom-specs.hh"
#include "validation-graphs/validation-information.hh"
#include "validation-graphs/validation-graphs.hh"
// need gtk things
#include <gtk/gtk.h>
#include <epoxy/gl.h>
#include <map>
#include <memory>

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif // HAVE_VECTOR

#define GLM_ENABLE_EXPERIMENTAL // # for norm things
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>

#include <utils/ctpl.h>
#include <coot-utils/inchikey-store.hh>

#ifdef USE_MOLECULES_TO_TRIANGLES
// #include <MoleculesToTriangles/CXXClasses/RendererGLSL.hpp>
#endif // USE_MOLECULES_TO_TRIANGLES

#include "ft-character.hh"

#include "clipper/core/xmap.h"

#include "coords/Cartesian.hh"
#include "ccp4mg-utils/mgtree.h"

#include "pick.hh" // 20220723-PE no picking for WebAssembly build

#include "compat/coot-sysdep.h"
#include "command-arg.hh"

#include "ligand/rotamer.hh"

#include "geometry/protein-geometry.hh"

#include "stereo-eye.hh"
#include "molecule-class-info.h"

#ifdef HAVE_SSMLIB
#include <ssm/ssm_align.h>
#endif

#include "clip-hybrid-atom.hh"
#include "atom-label-info.hh"

#include "db-main/db-main.hh"
#include "build/CalphaBuild.hh"
#include "ideal/simple-restraint.hh"
#include "coot-utils/positron.hh"
#include "api/cell.hh"

// #ifdef DO_GEOMETRY_GRAPHS
// #include "test-validation"
// #endif

#include "history_list.hh"

// 20220406-PE Temporarily remove GOOCANVAS-dependent functions
#define DO_SEQUENCE_VIEW
#define DO_RAMA_PLOT
#define DO_GEOMETRY_GRAPHS
#if (GTK_MAJOR_VERSION == 4) || (GTK_MINOR_VERSION == 94)
#undef DO_SEQUENCE_VIEW
#undef DO_RAMA_PLOT
#undef DO_GEOMETRY_GRAPHS
#endif

#undef DO_SEQUENCE_VIEW
#undef DO_RAMA_PLOT
#undef DO_GEOMETRY_GRAPHS

#ifdef DO_SEQUENCE_VIEW
#include "sequence-view.hh"
#endif
#ifdef DO_RAMA_PLOT
#include "rama_plot.hh"
#endif



#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/coot-coord-extras.hh"
#include "coot-utils/atom-overlaps.hh"

#include "positioned-widgets.h"

#include "coot-database.hh"

#include "mtz-column-auto-read.hh"

#include "atom-pull.hh"

#include "key-bindings.hh"

#ifdef USE_LIBCURL
#ifndef HAVE_CURL_H
#define HAVE_CURL_H
// defined in new python!?
#ifdef socklen_t
#undef socklen_t
#endif
#include <curl/curl.h>
#endif // HAVE_CURL_H
#endif

#include "Texture.hh"
#include "TextureMesh.hh"
#include "HUDMesh.hh"
#include "HUDTextureMesh.hh"
#include "Instanced-Markup-Mesh.hh"
#include "Model.hh"
#include "boids.hh"
#include "meshed-particle-container.hh"
// #include "graphics-ligand-view.hh"
#include "graphics-ligand-mesh-molecule.hh"
#include "restraints-editor.hh"
#include "framebuffer.hh"
#include "lights-info.hh"
#include "atom-label-info.hh"
#include "translation-gizmo.hh"

#ifdef USE_GUILE
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvolatile"
#include <libguile.h>
#pragma GCC diagnostic pop
#endif

enum { N_ATOMS_MEANS_BIG_MOLECULE = 400 };

#include "old-generic-display-object.hh"

#include "meshed-generic-display-object.hh"

#include "simple-distance-object.hh"

#include "gl-rama-plot.hh"
#include "glarea_tick_function.hh"

#include "utils/logging.hh" // 20221008-PE - aftter years of thinking about it, let's do it.

#include "rail-points.hh"

#include "extra-distance-restraint-markup.hh"

#include "labelled-button-info.hh"

#include "cfc-gui.hh"

#ifdef USE_BACKWARD
#include <utils/backward.hpp>
#endif

// /usr/include/libintl.h:51:14: error: expected unqualified-id before ‘const’
//    51 | extern char *dcgettext (const char *__domainname,
//       |              ^~~~~~~~~
// /usr/include/libintl.h:51:14: error: expected ‘)’ before ‘const’
// ../../coot/src/support.h:39:42: note: to match this ‘(’
//    39 | #  define dcgettext(Domain,Message,Type) (Message)
//       |                                          ^
// /usr/include/libintl.h:82:14: error: expected unqualified-id before ‘const’
//   82 | extern char *textdomain (const char *__domainname) __THROW;
//       |              ^~~~~~~~~~

#include "graphics-info-enums.hh" // still needs splitting

#include "view.hh"

#include "lsq-dialog-values.hh"
#include "select-atom-info.hh"
#include "gl-bits.hh"

const guint UPDATING_MAPS_TIMEOUT_FUNCTION_IDX_UNSET = 99999999;

class graphics_info_t {

   static int n_molecules_max;

   static bool in_side_by_side_stereo_mode; // user interfae should toggle this and it gets used
                                            // in render_scene()

   static std::pair<double, double> mouse_begin;
   static std::pair<double, double> mouse_clicked_begin;

   static float rotation_centre_x;
   static float rotation_centre_y;
   static float rotation_centre_z;

   static coot::Cartesian old_rotation_centre;
   static void set_old_rotation_centre(const coot::Cartesian &rc) {
     old_rotation_centre = rc;
   }
   static coot::Cartesian get_old_rotation_centre() {
     return old_rotation_centre;
   }

   // delete these when working
   // static float old_rotation_centre_x;
   // static float old_rotation_centre_y;
   // static float old_rotation_centre_z;

   static long int T0;
   static long int Frames;
   //
   static short int active_map_drag_flag;


   // return the value "did the scroll happen/was it set up" - or not?
   bool smooth_scroll_maybe(float x, float y, float z,
			    bool do_zoom_and_move_flag,
			    float target_zoom);
   bool smooth_scroll_maybe_sinusoidal_acceleration(float x, float y, float z,
                                                    short int do_zoom_and_move_flag,
                                                    float target_zoom);
   void smooth_scroll_maybe_stepped_acceleration(float x, float y, float z,
			    short int do_zoom_and_move_flag,
			    float target_zoom);


   static float trackball_size;

   // Go To Atom privates:
   static std::string go_to_atom_chain_;
   static std::string go_to_atom_atom_name_;
   static int         go_to_atom_residue_;
   static int         go_to_atom_molecule_;
   static std::string go_to_atom_atom_altLoc_;
   static std::string go_to_atom_inscode_;

   // distance static
   static coot::Cartesian distance_pos_1;

   // angle/torsion static
   static coot::Cartesian angle_tor_pos_1;
   static coot::Cartesian angle_tor_pos_2;
   static coot::Cartesian angle_tor_pos_3;
   static coot::Cartesian angle_tor_pos_4;

   static GtkBuilder *gtkbuilder; // use this for widget lookups
   static GtkBuilder *preferences_gtkbuilder; // use this for widget lookups in the preferences dialog
   // public static GtkApplication *application; // below. used for finding the menu bar
   //
   static GtkWidget *display_control_window_;
   static GtkWidget *get_display_control_window();

   //
   // static int mol_no_for_environment_distances is public
   static graphical_bonds_container regularize_object_bonds_box;
   static graphical_bonds_container environment_object_bonds_box;
   static graphical_bonds_container symmetry_environment_object_bonds_box;

   //
   static GdkModifierType button_1_mask_;
   static GdkModifierType button_2_mask_;
   static GdkModifierType button_3_mask_;

   static bool find_ligand_do_real_space_refine_;
   static int find_ligand_protein_mol_;
   static int find_ligand_map_mol_;
   static std::vector<std::pair<int, bool> > *find_ligand_ligand_mols_; // contain a molecule number
                                                                        // and flag for is_wiggly?
   // 20230417-PE I don't want this to be a pointer these days. Change it
   static coot::protein_geometry* geom_p;

   static coot::rotamer_probability_tables rot_prob_tables;

   // now used for residue_score (from c-interface.h)
   // coot::rotamer_probability_info_t get_rotamer_probability(...)

   static std::atomic<bool> moving_atoms_lock;
   static std::atomic<unsigned int> moving_atoms_bonds_lock; // regularize_object_bonds_box is being updated

   // make this private again when bugs are fixed.
   static atom_selection_container_t *moving_atoms_asc;
   // static molecule_class_info_t moving_atoms_molecule; // used as a container for glsl variables. public access
   mmdb::Residue *get_first_res_of_moving_atoms();
   static int imol_moving_atoms;
   static int imol_refinement_map;
   static int moving_atoms_n_cis_peptides;
   static bool moving_atoms_have_hydrogens_displayed;

   //
   static int undo_molecule; // -1 initially

   // No, we don't want MMDBManager in the include files.
   // (or indeed mmdb-extras, etc...)
   //
/*    //  */
/*    void create_regularized_graphical_object(const std::string chain_id_1, */
/* 						     int resno_1, */
/* 						     int resno_2, */
/* 						     MMDBManager *results);  */


   // db-main
   //
   static coot::db_main main_chain;

   // flash picked intermediate atom - try to stop Erik going insane.
   static bool flash_intermediate_atom_pick_flag;
   static clipper::Coord_orth intermediate_flash_point;
   // function that uses this, picked_intermediate_atom_graphics_object() is public

   void environment_graphics_object_internal(const graphical_bonds_container &env_bonds_box) const;
   void environment_graphics_object_internal_lines(const graphical_bonds_container &env_bonds_box) const;
   void environment_graphics_object_internal_tubes(const graphical_bonds_container &env_bonds_box) const;
   void environment_graphics_object_internal_tube(const coot::CartesianPair &pair,
						  int ipart, int n_parts) const;
   void graphics_object_internal_single_tube(const coot::Cartesian &base_point,
					     const coot::Cartesian &end_point,
					     const double &radius,
					     const coot::tube_end_t &end_type) const;
   void graphics_object_internal_arrow(const coot::Cartesian &base_point,
				       const coot::Cartesian &end_point,
				       float fraction_head_size,
				       const double &radius) const;

   void graphics_object_internal_torus(const coot::Cartesian &base_point,
				       const coot::Cartesian &end_point,
				       const double &radius_1,
				       const double &radius_2,
				       int n_ring_atoms) const;

   void graphics_object_internal_arc(float start_angle,
				     float end_angle,
				     const coot::Cartesian &start_point,
				     const coot::Cartesian &start_dir,
				     const coot::Cartesian &normal,
				     float radius, float radius_inner);

   void graphics_object_internal_dodec(const coot::old_generic_display_object_t::dodec_t &dodec);

   void graphics_object_internal_pentakis_dodec(const coot::old_generic_display_object_t::pentakis_dodec_t &penta_dodec);

   void read_standard_residues();   // for mutation, we have
				    // pre-prepared a pdb file with
				    // residues in Buccaneer "Standard
				    // Orientation", idealized and in
				    // the most likely rotamer.

   // state
   static short int state_language; // a bit-tested variable, 1 = scheme,
                                    // 2 = python, 3 = both.

   std::string state_command(const std::string &name_space, const std::string &str,                 short int state_lang) const;
   std::string state_command(const std::string &name_space, const std::string &str, int i,          short int state_lang) const;
   std::string state_command(const std::string &name_space, const std::string &str, int i1, int i2, short int state_lang) const;
   std::string state_command(const std::string &name_space, const std::string &str, float f,        short int state_lang) const;
   std::string state_command(const std::string &str, float f,        short int state_lang, short unsigned int v) const;
   std::string state_command(const std::string &name_space, const std::string &str, float f1, float f2, float f3, short int state_lang) const;
   std::string state_command(const std::string &name_space, const std::string &str, const std::string &str2, short int state_lang);

   // baton stuff
   static coot::Cartesian baton_root;
   static coot::Cartesian baton_tip;
   static float baton_length;
   static std::vector<coot::scored_skel_coord> baton_next_ca_options;
   // baton_previous_ca_positions->back() is the point closest to the new
   // baton tip (is the baton root)
   static std::vector<clipper::Coord_orth> baton_previous_ca_positions; // up to 3.
   coot::Cartesian non_skeleton_tip_pos() const;
   void baton_next_directions(int imol_for_skel, mmdb::Atom *atom, const coot::Cartesian& pos,
			      const clipper::Coord_grid &cg_start,
			      short int use_cg_start);
   coot::Cartesian baton_tip_by_ca_option(int index) const;
   static int baton_next_ca_options_index;
   static int user_set_baton_atom_molecule; // -1 if not set (default).
   // imol_for_baton_atoms?  see baton_build_atoms_molecule();

   // starting point and direction parameters:
   static int baton_build_start_resno;
   static std::string baton_build_chain_id;
   static short int baton_build_direction_flag; // +1 for forwards, -1 for backwards
   static short int baton_build_params_active; // usually they are
					       // ignored because we
					       // get the residue
					       // number, etc from the
					       // previous atoms.


   int imol_for_skeleton() const;

   void create_molecule_and_display(std::vector<coot::scored_skel_coord> &pos_position,
				    const std::string &molname);
   // as above, except we update molecule with name molname to the
   // pos_positions (and delete everything else).
   //
   void update_molecule_to(std::vector<coot::scored_skel_coord> &pos_position,
			   const std::string &molname);
   int create_empty_molecule(const std::string &molname);

   // for multi-threading
   static void update_maps_for_mols(const std::vector<int> &mol_idxs);

   // symm_atom_pick (public) uses this (private) function:

   void
   fill_hybrid_atoms(std::vector<coot::clip_hybrid_atom> *hybrid_atoms,
		     const atom_selection_container_t &asc,
		     const clipper::Spacegroup &spg,
		     const clipper::Cell &cell,
		     const std::pair<symm_trans_t, Cell_Translation> &symm_trans) const;

   // directory saving for fileselection
   //
   static std::string directory_for_fileselection;
   static std::string directory_for_saving_for_fileselection;
   static std::string directory_for_filechooser;
   static std::string directory_for_saving_for_filechooser;

   // distance object vector, and angle,
   // 20211006-PE both the vectors and the Meshes are needed because I will need to
   // rebuild the mesh if a distance is deleted. (Angle same should that happen one day)
   static std::vector<coot::simple_distance_object_t> measure_distance_object_vec;
   static std::vector<coot::coord_orth_triple> measure_angle_object_vec;
   static Mesh mesh_for_measure_distance_object_vec;
   static Mesh mesh_for_measure_angle_object_vec;

   // 20180217 moving_atoms_dragged_atom_index -> moving_atoms_dragged_atom_indices
   //          Now we can have many dragged atoms
   //
   static int moving_atoms_currently_dragged_atom_index;
   static std::set<int> moving_atoms_dragged_atom_indices;
   static void remove_moving_atoms_dragged_atom_index(int idx);
   static void    add_moving_atoms_dragged_atom_index(int idx);
   // make unset_moving_atoms_currently_dragged_atom_index() public


   static coot::restraints_container_t *last_restraints;
   // the mode flag is public:
   // 20220504-PE so that I can check for cleared/removed non-bonded contact baddies
   static std::map<int, std::vector<int> > previous_round_nbc_baddies_atom_index_map;

   // return the state of having found restraints.
   bool make_last_restraints(const std::vector<std::pair<bool,mmdb::Residue *> > &local_resiudes,
			     const std::vector<mmdb::Link> &links,
			     const coot::protein_geometry &geom,
			     mmdb::Manager *mol_for_residue_selection,
			     const std::vector<coot::atom_spec_t> &fixed_atom_specs,
			     coot::restraint_usage_Flags flags,
			     bool use_map_flag,
			     const clipper::Xmap<float> *xmap_p);

   // which uses the following...
#ifdef USE_GUILE
   void run_post_manipulation_hook_scm(int imol, int mode);
#endif
#ifdef USE_PYTHON
   void run_post_manipulation_hook_py(int imol, int mode);
#endif

   // void run_post_read_model_hook(int imol); // now public as it is called from handle_read_draw_molecule()

   void run_post_set_rotation_centre_hook();
   // which uses the following...
#ifdef USE_GUILE
   void run_post_set_rotation_centre_hook_scm();
#endif
#ifdef USE_PYTHON
   void run_post_set_rotation_centre_hook_py();
#endif

   // edit ramachandran store:
   static coot::ramachandran_points_container_t rama_points;
   std::pair<std::pair<double, double>, std::pair<double, double> >
   phi_psi_pairs_from_moving_atoms();

   // end points of backbone torsion moving
   static clipper::Coord_orth backbone_torsion_end_ca_1;
   static clipper::Coord_orth backbone_torsion_end_ca_2;
   static int backbone_torsion_peptide_button_start_pos_x;
   static int backbone_torsion_peptide_button_start_pos_y;
   static int backbone_torsion_carbonyl_button_start_pos_x;
   static int backbone_torsion_carbonyl_button_start_pos_y;

   // We use this ramachandran_points_container to pass ramachan plots
   // to the bond (and markup atom) generator (Bond-lines).
   //
   static ramachandrans_container_t ramachandrans_container;

   clipper::Coord_orth moving_atoms_centre() const;
   void set_edit_backbone_adjustments(GtkWidget *widget);
   static void edit_backbone_peptide_changed_func (GtkAdjustment *adj, GtkWidget *window); // callback
   static void edit_backbone_carbonyl_changed_func(GtkAdjustment *adj, GtkWidget *window); // callback

   void check_and_warn_inverted_chirals_and_cis_peptides() const;

#ifdef DO_RAMA_PLOT
   void handle_rama_plot_update(coot::rama_plot *plot);
#endif

   // Geometry Graphs:
#ifdef DO_GEOMETRY_GRAPHS
   coot::geometry_graphs *geometry_graph_dialog_to_object(GtkWidget *w) const;
#endif

   std::vector<coot::geometry_distortion_info_container_t>
     geometric_distortions_from_mol(int imol, const atom_selection_container_t &asc, bool with_nbcs);
   void print_geometry_distortion(const std::vector<coot::geometry_distortion_info_container_t> &v) const;

   void check_if_in_range_defines();
   void check_if_in_rotamer_define_gtk4(const pick_info &naii);

   static std::vector<std::string> model_fit_refine_toggle_button_name_list();
   static std::vector<std::string> model_fit_refine_button_name_list();
   static std::vector<std::string> other_modelling_tools_toggle_button_name_list();
   static std::vector<std::string> other_modelling_tools_button_name_list();

   // refinement_results_t is in ideal/simple-restraints.hh

   coot::refinement_results_t
     copy_mol_and_refine_inner(int imol_for_atoms,
			       int resno_1,
			       int resno_2,
			       int nSelResidues,
			       mmdb::PResidue *SelResidues,
			       const std::string &chain_id_1,
			       const std::string &altconf,
			       short int have_flanking_residue_at_start,
			       short int have_flanking_residue_at_end,
			       int imol_for_map);


   // rename me
   coot::refinement_results_t
     update_refinement_atoms(int n_restraints,
			     coot::restraints_container_t *restraints, // actually last_restraints
			     coot::refinement_results_t rr,
			     atom_selection_container_t local_mov_ats,
			     bool need_residue_order_check,
			     int imol,
			     std::string chain_id);

   static void refinement_loop_threaded();
   // several function move the atoms and need some refinement afterward (pepflip,JED-flip)
   // they all want to continue the refinement of last_restraints - but don't need to
   // start (and detach) a new thread to do so if a refinement thread is already running.
   // Hence thread_for_refinement_loop_threaded() which replaces use of
   // add_drag_refine_idle_function()
   // drag_refine_refine_intermediate_atoms()
   //
   static void thread_for_refinement_loop_threaded();
   static void refinement_of_last_restraints_needs_reset();
   static bool refinement_of_last_restraints_needs_reset_flag;
   void update_restraints_with_atom_pull_restraints(); // make static also?
   coot::restraint_usage_Flags set_refinement_flags() const; // make static?
   void debug_refinement();

   static void get_restraints_lock(const std::string &calling_function_name);
   static void release_restraints_lock(const std::string &calling_function_name);
   static std::string restraints_locking_function_name; //  static because it is set by above

   // similar for moving atoms:
   static void get_moving_atoms_lock(const std::string &calling_function_name);
   static void release_moving_atoms_lock(const std::string &calling_function_name);
   static std::string moving_atoms_locking_function_name; //  static because it is set by above


   // 201803004:
   // refinement now uses references to Xmaps.
   // A dummy_map is created and a reference to that is created. Then
   // the reference is reset to a real xmap in a molecule (imol_for_map).
   // But, for a reason I don't understand, the refinement crashes when I do that.
   // When the initial dummy_xmap doesn't go out of scope, then the refinement is OK.
   // So this static dummy map is the map that doesn't go out of scope.
   // We only need one of it, so it goes here, rather than get created every
   // time we do a refinement. It may need to be public in future.
   static clipper::Xmap<float> *dummy_xmap;

   std::string adjust_refinement_residue_name(const std::string &resname) const;
   static void info_dialog_missing_refinement_residues(const std::vector<std::string> &res_names);
   // replaced by public
   // static void show_missing_refinement_residues_dialog(const std::vector<std::string> &res_names);

   void info_dialog_alignment(coot::chain_mutation_info_container_t mutation_info) const;
   void info_dialog_refinement_non_matching_atoms(std::vector<std::pair<mmdb::Residue *, std::vector<std::string> > > nma);

   // bottom left flat ligand view:
   //
   static bool graphics_ligand_view_flag;
   static int  graphics_ligand_view_imol; // the molecle of the ligand. We don't want to show the HUG ligand
                                          // if the molecule of which it is part is not displayed.

   // ----------------------------------------------------------------
   //             public:
   // ----------------------------------------------------------------

public:

   enum { USE_PYTHON_STATE_COMMANDS = 2, USE_SCM_STATE_COMMANDS = 1 };

   //
   void init();
   void setup_key_bindings();

   static bool graphics_is_gl_es; // can we make a fallback coot where this is true?

   static bool use_gemmi;
   void set_use_gemmi(bool state) { use_gemmi = state; }

   static bool coot_is_a_python_module; //turned off in main()
   static bool prefer_python;

   static bool do_expose_swap_buffers_flag;

   static std::vector<coot::positron_metadata_t> positron_metadata;

#ifdef USE_GUILE
   static bool scm_boot_guile_booted; // false until my_wrap_scm_boot_guile() has been run
#endif

   // static std::queue<std::chrono::time_point<std::chrono::high_resolution_clock> > frame_draw_queue;
   static std::chrono::time_point<std::chrono::high_resolution_clock> previous_frame_time;
   static std::chrono::time_point<std::chrono::high_resolution_clock> previous_frame_time_for_per_second_counter;

   static void graphics_grab_focus();
   static void graphics_draw();

   // sometimes (when we have 100s of molecules, we don't want to redraw when a molecule
   // is displayed or undisplayed)
   static bool mol_displayed_toggle_do_redraw; // normally true

   static void toggle_display_of_last_model();

   static bool is_valid_model_molecule(int imol);
   static bool is_valid_map_molecule(int imol);
   static bool is_difference_map(int imol);
   static bool display_mode_use_secondary_p();

   // New-style Application!
   static GtkApplication *application; // used for finding the menu bar

   // New-style gtkbuilder!
   static void set_gtkbuilder(GtkBuilder *builder) { gtkbuilder = builder; }
   static void set_preferences_gtkbuilder(GtkBuilder *builder) { preferences_gtkbuilder = builder; }
   static GtkWidget *get_widget_from_builder(const std::string &w_name); // use gtkbuilder to do new-style lookup_widget();
   static GtkWidget *get_widget_from_preferences_builder(const std::string &w_name); // use gtkbuilder to do new-style lookup_widget();
   static bool gui_from_gtkbuilder() { return (gtkbuilder == NULL) ? false : true; }
   static GObject *get_gobject_from_builder(const std::string &w_name); // use gtkbuilder but return a gobject (for menus)

   enum {GL_CONTEXT_MAIN = 0, GL_CONTEXT_SECONDARY = 1};

   static void make_gl_context_current(bool gl_context_current_request_index);

   // ------------- main window -----------------------
   static GtkWidget *main_window;
   static void set_main_window(GtkWidget *w) { main_window = w; }
   static GtkWidget *get_main_window() { return main_window; }
   // ------------- glareas -----------------------
   std::chrono::time_point<std::chrono::system_clock> tp_now;

   static std::vector<GtkWidget *> glareas;
   static GtkAllocation get_glarea_allocation();
   static gl_context_info_t get_gl_context_info();

   // we need to store these because when we want to correct the size
   // and position of HUD objects (refinemetn arrow, HUD refinement buttons)
   // we need to know them.
   static int hud_start_graphics_window_x_width;
   static int hud_start_graphics_window_x_height;
   static std::string main_window_title;
   static short int model_fit_refine_dialog_was_sucked;

   // ------------- statusbar -----------------------
   static GtkWidget *statusbar;
   static guint statusbar_context_id;
   static void add_status_bar_text(const std::string &text);

   static void statusbar_ctrl_key_info(); // Ctrl to rotate or pick?

   // -------------------------------------------------


   static bool convert_to_v2_atom_names_flag; // shall we convert nucletotides to
				       // match the dictionary names? Often
				       // we want to do this (give current
				       // Coot architecture).  Sometimes
				       // not, though.

   // To be used to (typically) get the menu item text label from chain
   // option menus (rather than the ugly/broken casting of
   // GtkPositionType data.
   static std::string menu_item_label(GtkWidget *menu_item);

   // accept/reject window, now controlled by keyboarding in main window.
   static GtkWidget *accept_reject_dialog;
   static GtkWidget *refine_params_dialog;

   void save_accept_reject_dialog_window_position(GtkWidget *acc_rej_dialog);

   // flag to display the accept/reject dialog in the toolbar
   static int accept_reject_dialog_docked_flag;

   // flag to show/hide/sensitise docked accept/reject dialog
   static int accept_reject_dialog_docked_show_flag;

   // flag for the refinement toolbar show/hide
   static short int model_toolbar_show_hide_state;

   // flag for the refinement toolbar position
   static short int model_toolbar_position_state;

   // flag for the refinement toolbar style
   static short int model_toolbar_style_state;

   // flag for the main toolbar show/hide
   static short int main_toolbar_show_hide_state;

   // flag for the main toolbar position
   static short int main_toolbar_position_state;

   // flag for the main toolbar style
   static short int main_toolbar_style_state;

   static short int do_lighting_flag;
   static bool do_flat_shading_for_solid_density_surface;

   static bool sequence_view_is_docked_flag;

   static bool validation_graphs_is_docked; // is in main window

   static short int do_anti_aliasing_flag; // BL feature
   void set_do_anti_aliasing(int state);
   void draw_anti_aliasing();
   static int display_mode; // e.g. HARDWARE_STEREO_MODE, DTI_SIDE_BY_SIDE_STEREO
   static float stereo_angle;  // 6 degrees default
   static bool in_wall_eyed_side_by_side_stereo_mode;
   // static stereo_eye_t which_eye;
   static glm::vec3 eye_position; // useful in projection (testing)
   static bool stereo_style_2010;

   // return a vector of the current valid map molecules
   std::vector<int> valid_map_molecules() const;

   static int create_molecule(); // return the new molecule number
   static void erase_last_molecule();
   static int get_latest_model_molecule();
   static int get_biggest_model_molecule();

   static bool use_graphics_interface_flag; // 20220409-PE now defaults is false!

   // Display size
   static int graphics_x_size;
   static int graphics_y_size;
   static int graphics_x_position;
   static int graphics_y_position;

   static int model_fit_refine_x_position;
   static int model_fit_refine_y_position;

   static int display_manager_x_size;
   static int display_manager_y_size;
   static int display_manager_x_position;
   static int display_manager_y_position;

   static int display_manager_molecules_vbox_x_size;
   static int display_manager_molecules_vbox_y_size;
   static int display_manager_maps_vbox_x_size;
   static int display_manager_maps_vbox_y_size;

   static int display_manager_paned_position;

   static int go_to_atom_window_x_position;
   static int go_to_atom_window_y_position;

   static int rotate_translate_x_position;
   static int rotate_translate_y_position;

   static int delete_item_widget_x_position;
   static int delete_item_widget_y_position;

   static int accept_reject_dialog_x_position;
   static int accept_reject_dialog_y_position;

   static int model_fit_refine_dialog_stays_on_top_flag;
   static std::string model_fit_refine_place_atom_at_pointer_string;
   static std::string model_fit_refine_rotate_translate_zone_string;

   static int edit_chi_angles_dialog_x_position;
   static int edit_chi_angles_dialog_y_position;

   static int rotamer_selection_dialog_x_position;
   static int rotamer_selection_dialog_y_position;

   static int ramachandran_plot_x_position;
   static int ramachandran_plot_y_position;

   static bool stroke_characters;
   static void printString(const std::string &s,
			   const double &x, const double &y, const double &z);
   static void printString_for_axes(const std::string &s,
				    const double &x, const double &y, const double &z);
   static void printString_for_density_level(const std::string &s,
					     const double &x, const double &y, const double &z);
   static void printString_internal(const std::string &s,
				    const double &x, const double &y, const double &z,
				    bool do_unproject, bool mono_font, double scale_factor);

   std::string get_directory_for_fileselection() const; // 20250725-PE can this be removed now?
   std::string get_directory_for_filechooser() const;

   static int map_line_width;

   static double mouse_current_x, mouse_current_y;
   // static float* quat;  // rotation quaternion, allocatated [4]
   // static float* baton_quat;

   static int mouse_just_cliked; // delete this?
   static float zoom;
   static short int quanta_like_zoom_flag;
   static void mouse_zoom(double delta_x, double delta_y);
   static void mouse_zoom_by_scale_factor(double sf);
   static void mouse_zoom_by_scale_factor_inner(double sf);
   static void scroll_zoom(int direction);
   static void handle_delete_item_curor_change(GtkWidget *widget);

   static float box_radius_xray;
   static float box_radius_em;

   static float iso_level_increment;
   static float diff_map_iso_level_increment;
   static short int swap_difference_map_colours; // for Jan Dohnalek

   static float map_sampling_rate; // Shannon sampling rate multiplier (1.5 default)

   //
   static int control_is_pressed;
   static int shift_is_pressed;
   static int y_is_pressed;
   static int z_is_pressed;
   static int a_is_pressed;
   static short int control_key_for_rotate_flag;
   static short int pick_pending_flag;
   // 20051004 now use use this interface rather direct access to pick_pending_flag
   //          because a overal pick_pending_flag doesn't work when we consider turning
   //          off various different model/fit/refine toggle-buttons.  So
   //          static_graphics_pick_pending() checks each of the pending picks which
   //          would make the graphics not rotate (ctrl key issues).
   short int static_graphics_pick_pending() const;

   static coot::residue_spec_t current_residue; // to be updated on set_go_to_atom and
		   			              // middle-mouse recentring.  check for
					              // unset_p() when used.
   static int current_residue_imol;


   // -------------- gtk4 event controller -----------------------------

   void on_glarea_drag_begin_primary(GtkGestureDrag *gesture, double x, double y, GtkWidget *gl_area);
   void on_glarea_drag_update_primary(GtkGestureDrag *gesture, double delta_x, double delta_y, GtkWidget *gl_area);
   void on_glarea_drag_end_primary(GtkGestureDrag *gesture, double x, double y, GtkWidget *gl_area);

   void on_glarea_drag_begin_secondary(GtkGestureDrag *gesture, double x, double y, GtkWidget *gl_area);
   void on_glarea_drag_update_secondary(GtkGestureDrag *gesture, double delta_x, double delta_y, GtkWidget *gl_area);
   void on_glarea_drag_end_secondary(GtkGestureDrag *gesture, double x, double y, GtkWidget *gl_area);

   void on_glarea_drag_begin_middle(GtkGestureDrag *gesture, double x, double y, GtkWidget *gl_area);
   void on_glarea_drag_update_middle(GtkGestureDrag *gesture, double delta_x, double delta_y, GtkWidget *gl_area);
   void on_glarea_drag_end_middle(GtkGestureDrag *gesture, double x, double y, GtkWidget *gl_area);

   void on_glarea_click(GtkGestureClick* self,
                        gint n_press,
                        gdouble x,
                        gdouble y,
                        gpointer user_data);

   void on_glarea_scrolled(GtkEventControllerScroll *controller,
                           double                    dx,
                           double                    dy,
                           gpointer                  user_data);

   void on_glarea_motion(GtkEventControllerMotion* controller,
                         gdouble x,
                         gdouble y,
                         gpointer user_data);

   // step is usually 1 or -1, but can be 2 or 3 (with fast scrolling)
   void change_model_molecule_representation_mode(int step);

   void do_drag_pan_gtk3(GtkWidget *widget, double drag_delta_x, double drag_delta_y);
   // here drag_delta_y and drag_delta_y are the differences from
   // drag_begin_x and drag_begin_y.
   void do_drag_pan_gtk4(GtkWidget *widget, double drag_delta_x, double drag_delta_y);

   gboolean on_glarea_key_controller_key_pressed(GtkEventControllerKey *controller,
                                                 guint                  keyval,
                                                 guint                  keycode,
                                                 guint                  modifiers);
   void on_glarea_key_controller_key_released(GtkEventControllerKey *controller,
                                              guint                  keyval,
                                              guint                  keycode,
                                              guint                  modifiers);

   //
   static coot::colour_holder cell_colour;

   //
   // Frames per seconds control
   // (all we need to do is get and set).
   //

   // symm colour is a part of the molecule now
   static double symmetry_colour_merge_weight;
   // this is merged with the colour from standard atom colours
   static glm::vec4 symmetry_colour;
   static GdkRGBA symmetry_colour_to_rgba();
   static void rgba_to_symmetry_colour(GdkRGBA rgba);

   // Rotate colour map?
   static short int rotate_colour_map_on_read_pdb_flag;
   static float rotate_colour_map_on_read_pdb; // e.g. 5.0 (degrees)
   static short int rotate_colour_map_on_read_pdb_c_only_flag;

   static float rotate_colour_map_for_map; // e.g. 31.0 (degrees)

   //
   static float symmetry_search_radius;
   // static short int symmetry_as_calphas; // moved to per molecule basis
   // static short int symmetry_rotate_colour_map_flag; // do we want symmetry of other
						     // molecules to have a different
						     // colour [MOL]?
                                                     // moved to per molecule basis

   static int symmetry_shift_search_size; // the shift size for which_boxes. A hack.


   static float symmetry_operator_rotate_colour_map;
   // static int   symmetry_colour_by_symop_flag; // moved to per molecule basis
   // static int   symmetry_whole_chain_flag; moved to molecule_class_info_t
   static int   symmetry_atom_labels_expanded_flag;

   //
   static short int show_symmetry;
   static void update_symmetry(); // of models

   // Clipping Planes:
   static float clipping_front;
   static float clipping_back;
   void adjust_clipping(float d);
   void set_clipping_front(float v);
   void set_clipping_back(float v);
   float get_clipping_plane_front() const;
   float get_clipping_plane_back() const;

   void increase_clipping_front();
   void increase_clipping_back();
   void decrease_clipping_front();
   void decrease_clipping_back();

   // This is for the display object
   static short int display_lists_for_maps_flag;

   // expose this so that it can be seen in draw();
   static short int smooth_scroll_on; // flag used to show wirecube centre
                                      // (not anything else so rename it?)
   std::vector<int> displayed_map_imols() const;

   // expose so that they can be used in c-interface.cc
   static int smooth_scroll;
   static int smooth_scroll_n_steps;
   static float smooth_scroll_limit;
   static float smooth_scroll_zoom_limit; // above this value we zoom, if zoom is on.
   static int   smooth_scroll_do_zoom;
   static bool smooth_scroll_on_going;
   static int go_to_ligand_animate_view_n_steps; // 50 default

   static
   gboolean smooth_scroll_animation_func(GtkWidget *widget,
                                         GdkFrameClock *frame_clock,
                                         gpointer data);
   static
   gboolean smooth_sinusoidal_scroll_animation_func(GtkWidget *widget,
                                                    GdkFrameClock *frame_clock,
                                                    gpointer data);

   // in that function, we need to know the current step
   static int smooth_scroll_current_step;
   // and the position delta (position at the end of the animation - the postion at the start of the animation)
   static coot::Cartesian smooth_scroll_delta;

   // for reorienting-next-residue, we need to store the start and end view
   static coot::view_info_t reorienting_residue_start_view;
   static coot::view_info_t reorienting_residue_end_view;

   // atom label font size
   static int atom_label_font_size; // range of 1->3. small, medium, large.
   static void *atom_label_font;
   static int label_atom_on_recentre_flag;
   static coot::colour_holder font_colour;
   void try_label_unlabel_active_atom();// l key binding
   void remove_all_atom_labels();

   static int n_molecules() { return molecules.size();}
   static int n_map_molecules();
   static std::vector<molecule_class_info_t> molecules;

   // To which map is the mouse scroll wheel attached?
   //
   static int scroll_wheel_map;
   void update_scroll_wheel_map_on_molecule_close();
   int intelligent_get_scroll_wheel_map() const;

   void contour_level_scroll_scrollable_map(int direction);

   // a static utility function (excised for portability)
   static std::string add_dir_file(const std::string &dirname, const std::string &filename);

   /*! \brief is the given file name suitable to be read as coordinates? */
   short int file_type_coords(const std::string &file_name);

   // state_command (now public, it's called from c-interface-build (mutate sequence)
   //
   std::string state_command(const std::vector<std::string> &str, short int state_lang) const;
   //
   // 20220828-PE Let's pass the module information
   std::string state_command(const std::string &module, const std::string &func_name,
                             const std::vector<coot::command_arg_t> &args, short int state_lang) const;

   // esoteric depth cue on/off  (on by default)
   static int esoteric_depth_cue_flag;

   // external c interface get/set functions.
   //
   static void ShowFPS();
   static int  GetFPSFlag() {return show_fps_flag; }
   static void SetShowFPS(int t);

   void SetActiveMapDrag(int t);
   short int GetActiveMapDrag() const;

   // return -1 on error
   //
   int lookup_molecule_name(const std::string &molecule_name) const;

   void SetMouseBegin(double x, double y);
   void SetMouseClicked(double x, double y);

   double GetMouseBeginX() const;
   double GetMouseBeginY() const;
   double GetMouseClickedX() const {return mouse_clicked_begin.first;}
   double GetMouseClickedY() const {return mouse_clicked_begin.second;}

   // We are given atom atom index, we will use this to look up
   // the atom using molecule_class_info and find its coordinates.
   // Those coordinates will get used in draw() to centre on that
   // atom.
   //
   void setRotationCentre(int atom_index, int imol);

   // if dir is true, we are going forward
   void reorienting_next_residue(bool dir);
   static bool reorienting_next_residue_mode;

   // return the "I did a jump" status
   bool setRotationCentre(coot::Cartesian centre, bool force_jump=false);
   void setRotationCentreAndZoom(coot::Cartesian centre,
				 float target_zoom);
   void setRotationCentreSimple(const coot::Cartesian &c);

   void setRotationCentre(const coot::clip_hybrid_atom &hybrid_atom);

   // old style: soon to be redundent
   void setRotationCentre(const symm_atom_info_t &symm_atom_info);

   static void set_rotation_centre(const clipper::Coord_orth &pt);
   void run_post_manipulation_hook(int imol, int mode);

   void run_post_read_model_hook(int imol); // conditional compilation handled internally

   void update_things_on_move();
   void update_things_on_move_and_redraw();

#ifdef DO_RAMA_PLOT
   void update_ramachandran_plot_point_maybe(int imol, mmdb::Atom *atom);
   void update_ramachandran_plot_point_maybe(int imol, const coot::residue_spec_t &res_spec);
   void update_ramachandran_plot_point_maybe(int imol, atom_selection_container_t moving_atoms);
   void update_ramachandran_plot_background_from_res_spec(coot::rama_plot *plot, int imol,
                                                          const coot::residue_spec_t &res_spec);
#endif

   float X() { return rotation_centre_x; };
   float Y() { return rotation_centre_y; };
   float Z() { return rotation_centre_z; };

   static coot::Cartesian RotationCentre() { return coot::Cartesian(rotation_centre_x,
                                                                    rotation_centre_y,
                                                                    rotation_centre_z);}
   // modern API
   static coot::Cartesian get_rotation_centre_cart() { return RotationCentre(); }

   // we need static, so that we don't need to instance a
   // graphics_info_t for every frame draw.
   static float RotationCentre_x() { return rotation_centre_x; }
   static float RotationCentre_y() { return rotation_centre_y; }
   static float RotationCentre_z() { return rotation_centre_z; }

   static clipper::Coord_orth hole_start;
   static clipper::Coord_orth hole_end;
   void set_hole_start() { hole_start = get_rotation_centre_co(); }
   void set_hole_end()   { hole_end   = get_rotation_centre_co(); }

   static coot::Cartesian smooth_scroll_start_point;
   static coot::Cartesian smooth_scroll_target_point;

   // possibly for multi-threading, public access.
   void update_maps();

   // pointer: aka rotation centre:
   //
   void display_where_is_pointer() const {
      std::cout << "Pointer at" << RotationCentre() << std::endl;
   }

   std::pair<double, double> get_pointer_position_frac() const;

   // x_diff and y_diff are the scale factors to the x and y
   // drag vectors.
   //
   void add_to_RotationCentre(coot::CartesianPair x_y, double x_diff, double y_diff) {

      // x_drag and y_drag are the model space vector due
      // to screen x and y differences respectively.
      //
      coot::Cartesian x_drag = x_y.getStart();
      coot::Cartesian y_drag = x_y.getFinish();

      rotation_centre_x += x_drag.get_x()*x_diff + y_drag.get_x()*y_diff;
      rotation_centre_y += x_drag.get_y()*x_diff + y_drag.get_y()*y_diff;
      rotation_centre_z += x_drag.get_z()*x_diff + y_drag.get_z()*y_diff;

   };

   // simple
   static void add_vector_to_rotation_centre(const coot::Cartesian &vec);
   // old method - not sure that it's useful or even works now...
   void add_vector_to_RotationCentre(const coot::Cartesian &vec); // do the updates and redraw too.

   // map colour by scripting
   void set_last_map_colour(double f1, double f2, double f3) const;

   // sigma contour stepping?
   // Calling this turns it on.
   void set_last_map_contour_level(float f);
   void set_last_map_contour_level_by_sigma(float f);
   static bool auto_recontour_map_flag;

   //
   static float user_defined_rotation_centre_crosshairs_size_scale_factor;
   static glm::vec4 rotation_centre_cross_hairs_colour;
   void set_rotation_centre_cross_hairs_colour(const glm::vec4 &c);

   static void Increment_Frames();

   void set_font_size(int size);

   void update_map_colour_menu();

   static short int do_scroll_by_wheel_mouse_flag;
   //
   void set_scrollable_map(int imol);

   // Transfered from molecule_class_info, because they are a parameter
   // of the *graphis* at the moment, not each molecule (there is no
   // molecule control in Anisotropic Atoms at the momemnt).
   //
   static short int show_aniso_atoms_flag;
   static float     show_aniso_atoms_radius;
   static short int show_aniso_atoms_radius_flag; // shall the atoms be
                                           // limited to a certain distance?
   static float     show_aniso_atoms_probability;

   //
   void set_vt_surface(int v); // virtual trackball
   float get_trackball_size() const {return trackball_size; };
   int vt_surface_status() const;

   // skeleton colour
   static double* skeleton_colour;

   static gint idle_contour_function(gpointer data);

   // idle function token (holder)
   //
   static int idle_function_spin_rock_token;
   // static long time_holder_for_rocking;
   static std::chrono::time_point<std::chrono::high_resolution_clock> time_holder_for_rocking;
   // drag refine idle function token:
   static int drag_refine_idle_function_token;
   static float idle_function_rotate_angle; // degrees
   static double idle_function_rock_amplitude_scale_factor;
   static double idle_function_rock_freq_scale_factor;
   static double idle_function_rock_angle_previous;
   static gint drag_refine_idle_function(GtkWidget *widget);
   static void add_drag_refine_idle_function();
   static void remove_drag_refine_idle_function();

   static std::vector<std::chrono::time_point<std::chrono::high_resolution_clock> > leftquote_press_times;
   int get_n_pressed_for_leftquote_tap(std::chrono::time_point<std::chrono::high_resolution_clock> tp);

   static int drag_refine_refine_intermediate_atoms();

   static double refinement_drag_elasticity;
   static coot::refinement_results_t saved_dragged_refinement_results;
   static bool post_intermediate_atoms_moved_ready;
#ifdef USE_PYTHON
   static PyObject *post_intermediate_atoms_moved_hook;
   void register_post_intermediate_atoms_moved_hook(PyObject *function_name);
#endif
   void run_post_intermediate_atoms_moved_hook_maybe(); // set a python variable when the intermediate
                                                        // atoms move

   static bool showing_intermediate_atoms_from_refinement();

#ifdef USE_GUILE
   SCM refinement_results_to_scm(const coot::refinement_results_t &rr) const;
#endif
#ifdef USE_PYTHON
   PyObject *refinement_results_to_py(const coot::refinement_results_t &rr) const;
#endif
   static bool cryo_EM_refinement_flag;

   // ligand interactions (pulsing cylindrical bonds, or whatever)
   static long time_holder_for_ligand_interactions;
   static double ligand_interaction_pulse_previous; // a timing holder
   static int   idle_function_ligand_interactions_token;

   // now that the refinement goes via the idle function callback
   // (which gets called several times) we need to external control of
   // when to print the chi_squareds... and we turn them on when we
   // set the idle function.  Then once the function has run, we turn
   // them off.  Internally (in the restraints minimize()), on
   // GSL_SUCCESS and GSL_ENOPROG we print the chi squareds.
   static short int print_initial_chi_squareds_flag;

   // phs filename (because we actually do things when we OK on the
   // *coordinates* filename).  So give store and retrieve functions.
   //
   static std::string phs_filename; // make me private.
   std::string get_phs_filename() const;
   void set_phs_filename( std::string filename);
   // static int phs_cell_from_molecule;

   void clear_up_moving_atoms_wrapper(); // wraps the below functions - and is the Esc key function

   void clear_up_glsl_buffers_for_moving_atoms();
   // get rid of the actual molecule (as opposed to
   // clear_moving_atoms_object which removed the bonds).
   void clear_up_moving_atoms();

   // if the imol for moving atoms is imol, delete the moving atoms (called from close_molecule)
   void clear_up_moving_atoms_maybe(int imol);

   void delete_pointers_to_map_in_other_molecules(int imol_map);

   // 0: never run it
   // 1: ask to run it
   // 2: run it without asking
   static short int run_state_file_status;
   static bool state_file_was_run_flag;
   static bool run_startup_scripts_flag;


   // Go To Atom
   //

   const char *go_to_atom_chain();
   const char *go_to_atom_atom_name();
   const char *go_to_atom_ins_code();
   const char *go_to_atom_alt_conf();
   int go_to_atom_residue();
   int go_to_atom_molecule();
   static int go_to_atom_mol_menu_active_position;
   static int go_to_atom_menu_label_n_chars_max; // the last 30 (or so) chars
   std::string make_mmdb_atom_string_from_go_to_atom();
   static int         go_to_ligand_n_atoms_limit; // ligands must have at least this
						  // number of atoms for the "go to ligand"
                                                  // button and function to see it.
   static std::vector<std::string> go_to_ligand_non_interesting_comp_ids;

   void set_go_to_atom_chain_residue_atom_name(const std::string &chain_id,
					      int resno, const std::string &ins_code,
					      const std::string &atom_name, const std::string &altLoc);
   // 20220723-PE these look like GUI callbacks
   void set_go_to_atom_chain_residue_atom_name(const char *t1,
					       int it2, const char *t3);
   void set_go_to_atom_chain_residue_atom_name(const char *chain_id,
					       int resno, const char *atom_name, const char *altLoc);
   // 20230611-PE this is indeed a gui callback
   static void set_go_to_atom(int imol, const coot::atom_spec_t &spec);

   // 20230520-PE why isn't this here?
   // void set_go_to_atom(const coot::atom_spec_t &atom_spec);
   // I think it should be. Next time you come here, make it so and use it in
   // on_generic_atom_spec_toggle_button_toggled().

   void set_go_to_residue_intelligent(const std::string &chain_id, int resno,
				      const std::string &ins_code);

   // 20211015-PE why doesn't this function exist already?
   void go_to_residue(int imol, const coot::residue_spec_t &rs);
   static std::pair<std::string, std::string> split_atom_name(const std::string &atom_name);
   static std::pair<std::string, std::string> split_resno_inscode(const std::string &atom_name);

   mmdb::Atom *get_atom(int imol, const coot::atom_spec_t &spec) const;
   static mmdb::Residue *get_residue(int imol, const coot::residue_spec_t &spec);

   void set_go_to_atom_molecule(int pos);

   int try_centre_from_new_go_to_atom();

   // 20220723-PE There functions have (will have) internals that are blanked out by
   // EMSCRIPTEN conditional compilation
   void update_go_to_atom_window_on_changed_mol(int imol);
   void update_go_to_atom_window_on_new_mol();
   void update_go_to_atom_window_on_other_molecule_chosen(int imol);
   int update_go_to_atom_molecule_on_go_to_atom_molecule_deleted(); // return new gotoatom mol


   void update_widget_go_to_atom_values(GtkWidget *window, mmdb::Atom *atom);
   void make_synthetic_select_on_residue_list(GtkWidget *residue_list, mmdb::Atom *atom_p) const;

   void make_synthetic_select_on_residue_tree(GtkWidget *residue_list, mmdb::Atom *atom_p) const;
   void make_synthetic_select_on_residue_tree_gtk1(GtkWidget *residue_list, mmdb::Atom *atom_p) const;

   //int go_to_atom_molecule_optionmenu_active_molecule(GtkWidget *widget); // DELETE-ME

   static void fill_go_to_atom_window_residue_and_atom_lists_gtk4();
   static void fill_go_to_atom_window_gtk3(GtkWidget *widget);
   static void fill_go_to_atom_window_gtk2(GtkWidget *go_to_atom_window,
					   GtkWidget *residue_tree_scrolled_window,
					   GtkWidget *atom_list_scrolled_window);
   static void fill_go_to_atom_atom_list_gtk2(GtkWidget *atom_tree, int imol,
					      char *chain_id, int seqno, char *ins_code);

   static void go_to_atom_residue_tree_destroy(gpointer data);
   static void go_to_atom_list_destroy(gpointer data);


   static void clear_atom_list(GtkWidget *atom_gtklist);
   static void fill_go_to_atom_residue_list_gtk1(GtkWidget *gtklist);

   // goodbye gtk2
   // static void fill_go_to_atom_residue_tree_and_atom_list_gtk2(int imol,
   // GtkWidget *gtktree,
   // GtkWidget *atom_list);

   void fill_go_to_atom_option_menu(GtkWidget *option_menu);


   // goodbye my friends.

    void fill_option_menu_with_coordinates_options(GtkWidget *option_menu,
 						  GCallback callback_func);
   void fill_option_menu_with_coordinates_options(GtkWidget *option_menu,
						  GCallback signal_func,
						  int imol_active_position);
   void fill_option_menu_with_coordinates_options_internal(GtkWidget *option_menu,
							   GCallback callback_func,
							   short int set_last_active_flag);
   void fill_option_menu_with_coordinates_options_internal_2(GtkWidget *option_menu,
							     GCallback callback_func,
							     short int set_last_active_flag,
							     int imol_active);
   void fill_option_menu_with_coordinates_options_internal_3(GtkWidget *option_menu,
							     GCallback callback_func,
							     std::vector<int> fill_with_these_molecules,
							     short int set_last_active_flag,
							     int imol_active);
   void fill_option_menu_with_coordinates_options_internal_with_active_mol(GtkWidget *option_menu,
									   GCallback callback_func,
									   int imol_active);
   void fill_option_menu_with_coordinates_options_possibly_small(GtkWidget *option_menu,
								 GCallback callback_func,
								 int imol,
								 bool fill_with_small_molecule_only_flag);

   // use this.
   void new_fill_combobox_with_coordinates_options(GtkWidget *combobox_molecule,
                                                   GCallback callback_func,
                                                   int imol_active);

   void fill_combobox_with_coordinates_options(GtkWidget *combobox,
					       GCallback callback_func,
					       int imol_active);

   void fill_combobox_with_coordinates_options_with_set_last(GtkWidget *combobox,
							     GCallback callback_func,
							     bool set_last_active_flag);

   int combobox_get_imol(GtkComboBox *combobox) const;

   static void go_to_atom_mol_combobox_changed(GtkWidget *combobox, gpointer data);

   GtkWidget *dialog_box_of_buttons_internal(const std::string &window_title,
                                             const std::vector<std::tuple<std::string, GCallback, gpointer> > &buttons,
                                             const std::string &close_button_label);

#if 0

   // where is this used? If anywhere...
   static void on_go_to_atom_residue_list_selection_changed (GtkList *gtklist,
							     gpointer user_data);
#endif

   static void on_go_to_atom_residue_tree_selection_changed (GtkTreeView *gtklist,
							     gpointer user_data);
   static void
     residue_tree_residue_row_activated(GtkTreeView        *treeview,
					GtkTreePath        *path,
					GtkTreeViewColumn  *col,
					gpointer            userdata);
   static void
     atom_tree_atom_row_activated(GtkTreeView        *treeview,
				  GtkTreePath        *path,
				  GtkTreeViewColumn  *col,
				  gpointer            userdata);
   static gboolean
     residue_tree_selection_func(GtkTreeSelection *selection,
				 GtkTreeModel *model,
				 GtkTreePath *path,
				 gboolean path_currently_selected,
				 gpointer data);

   static gboolean
     atom_tree_selection_func(GtkTreeSelection *selection,
			      GtkTreeModel *model,
			      GtkTreePath *path,
			      gboolean path_currently_selected,
			      gpointer data);

// BL says:: put my gtk2 stuff in here too:
   static int gtk2_file_chooser_selector_flag;
   static int gtk2_chooser_overwrite_flag;

   void apply_go_to_atom_from_widget(GtkWidget *widget);
   static void pointer_atom_molecule_combobox_changed(GtkWidget *combobox, gpointer data);

   // return success status
   int intelligent_next_atom_centring(GtkWidget *widget);
   int intelligent_previous_atom_centring(GtkWidget *widget);
   int intelligent_near_atom_centring(GtkWidget *widget, const std::string &direction);

   // this can be used for symmetry atom pick:
   std::pair<coot::Cartesian, coot::Cartesian> get_front_and_back_for_pick() const;

   // int find_atom_index_in_moving_atoms(char *chain_id, int resno, char *atom_name) const;
   mmdb::Atom *find_atom_in_moving_atoms(const coot::atom_spec_t &at) const;

   pick_info find_atom_index_from_goto_info(int imol);
   pick_info atom_pick_gtk3(bool intermediate_atoms_only_flag) const;
   pick_info pick_moving_atoms(const coot::Cartesian &front, const coot::Cartesian &back) const;
   mmdb::Atom *get_moving_atom(const pick_info &pi) const; // return 0 on lookup failure

   coot::Symm_Atom_Pick_Info_t symmetry_atom_pick() const;
   coot::Symm_Atom_Pick_Info_t symmetry_atom_pick(const coot::Cartesian &front, const coot::Cartesian &back) const;
   coot::Symm_Atom_Pick_Info_t symmetry_atom_close_to_screen_centre() const;

   bool tomo_pick(double x, double y, gint n_press, bool shift_is_pressed);

   // map skeletonization level (and (different widget) boxsize).
   //
   static float skeleton_level;
   static float skeleton_box_radius;

   // autobuild
   static short int autobuild_flag;

   // file selection should be sorted by date?
   static short int sticky_sort_by_date;

   // file filter should be on?
   static short int sticky_file_filter;

   //
   void stop_refinement_internal();

   void show_refine_params_dialog(); // not used for map selection now.
   void show_select_map_frame(); // it's an overlay now
   void show_select_map_frame_gtkbuilder();
   void show_select_map_dialog_old_style();

   // Map and molecule display.  We need this so that we can look up
   // the names of the boxes so that we can add extra entries to them
   // when we create new maps and molecules
   //
   void save_display_control_widget_in_graphics(GtkWidget *widget);
   GtkWidget *display_control_window();

   static void activate_scroll_radio_button_in_display_manager(int imol);
   static float find_waters_sigma_cut_off;


   // for defining a range (to, say, regularize), we click on atoms
   // and manipulate atom indexes
   //
   // in_range_define has special values 0, 1, 2.  0 means not in a
   // range definition, 1 mean that we have just clicked the menu
   // item, 2 means that we are waiting for the last atom pick
   // (e.g. we have clicked on the first atom and are waiting for the
   // user to select the second).
   //
   static coot::atom_spec_t in_range_first_picked_atom;
   static coot::atom_spec_t in_range_second_picked_atom;
   static short int in_range_define; // initially 0
   static short int in_range_define_for_refine; // initially 0
   static int refine_regularize_max_residues;
   static int residue_range_mol_no;
   static int residue_range_atom_index_1;
   static int residue_range_atom_index_2;
   static int geometry_atom_index_1;
   static int geometry_atom_index_2;
   static int geometry_atom_index_3;
   static int geometry_atom_index_4;
   static int geometry_atom_index_1_mol_no;
   static int geometry_atom_index_2_mol_no;
   static int geometry_atom_index_3_mol_no;
   static int geometry_atom_index_4_mol_no;
   static short int fix_chiral_volume_before_refinement_flag;
   static short int show_chiral_volume_errors_dialog_flag;

   // torsion general
   static int torsion_general_atom_index_1;
   static int torsion_general_atom_index_2;
   static int torsion_general_atom_index_3;
   static int torsion_general_atom_index_4;
   static int torsion_general_atom_index_1_mol_no;
   static int torsion_general_atom_index_2_mol_no;
   static int torsion_general_atom_index_3_mol_no;
   static int torsion_general_atom_index_4_mol_no;
   static std::vector<coot::atom_spec_t> torsion_general_atom_specs;
   static bool torsion_general_reverse_flag;
   static Tree torsion_general_tree;
   static std::vector<std::vector<int> > torsion_general_contact_indices;

   //
   static int imol_pepflip;
   static int iresno_pepflip;
   static int atom_index_pepflip;
   static short int in_pepflip_define;
   //
   static int imol_rigid_body_refine;
   static short int in_rigid_body_define;
   // uses imol_refinement_map;
   static short int in_terminal_residue_define;
   static std::string add_terminal_residue_type;

   // CIS <-> TRANS conversion
   static int in_cis_trans_convert_define;

   // rotate/translate object mode
   static short int in_rot_trans_object_define;
   static short int rot_trans_object_type;
   static int rot_trans_atom_index_1;
   static int rot_trans_atom_index_2;
   static int imol_rot_trans_object;
   static short int rot_trans_zone_rotates_about_zone_centre;

   // additional representation
   static int add_reps_molecule_option_menu_item_select_molecule; // option menu
   static int add_reps_molecule_combobox_molecule;

   static coot::fixed_atom_pick_state_t in_fixed_atom_define;
   // static GtkWidget *fixed_atom_dialog; 20211202-PE extract this from glade each time.

   static short int in_torsion_general_define;
   // static int rot_trans_atom_index_rotation_origin_atom; old naive way.
   static mmdb::Atom *rot_trans_rotation_origin_atom; // "Eugene's way"

   static int imol_residue_partial_alt_locs;
   static short int in_residue_partial_alt_locs_define;
   static coot::residue_spec_t residue_partial_alt_locs_spec;
   void residue_partial_alt_locs_split_residue(int i_bond, bool wag_the_dog);
   static double residue_partial_alt_locs_rotate_fragment_angle;

   static short int in_user_defined_define;

   // save symmetry?
   static short int in_save_symmetry_define;

   // Where should we open up the save coords fileselection?
   static int save_coordinates_in_original_dir_flag;

   // Was private, but need to be used by auto_fit_best_rotamer() scripting function.
   void update_geometry_graphs(mmdb::PResidue *SelResidues, int nSelResidues, int imol_coords, int imol_map);
   void delete_residue_from_geometry_graphs(int imol, coot::residue_spec_t res_spec);
   void delete_residues_from_geometry_graphs(int imol, const std::vector<coot::residue_spec_t> &res_specs);
   void delete_chain_from_geometry_graphs(int imol, const std::string &chain_id);


   void execute_rotate_translate_ready(); // manual movement
   void unsetup_rotate_translate_buttons(GtkWidget *window); /* delete the user data */
   void do_rot_trans_adjustments(GtkWidget *dialog);
   static void rot_trans_adjustment_changed(GtkAdjustment *adj, gpointer user_data);
   static float *previous_rot_trans_adjustment;

   // rottrans_buttons class calls back this function on button pressed mouse motion
   //
   // void rot_trans_obj(int xdiff, const std::string &button_label); old

   void set_in_range_define_for_regularize(short int state) { in_range_define = state; } // true
   void set_in_range_define_for_refine(short int state) { in_range_define_for_refine = state; } // true
   void set_in_pepflip_define(short int state) { in_pepflip_define = state; }
   void set_in_rigid_body_refine(short int state) { in_rigid_body_define = state; }
   static float rigid_body_fit_acceptable_fit_fraction;
   static int refine_auto_range_step;  // +/- (1) about the clicked residue

   // called by scripting interface to rigid_body_refine_zone.  This
   // is not a great ame for the function (or for the residues that it
   // sets) since they do obviously refer to rigid body refinement.
   //
   void set_residue_range_refine_atoms(const std::string &chain_id,
					int resno_start, int resno_end,
				       const std::string &alt_conf,
				       int imol);
   void execute_rigid_body_refine(short int auto_range_flag);
   // called by above:
   // return the sucess status (0 for fail).
   // Replacing atom positions in imol_rigid_body_refine, so make sure
   // that you set that correctly before calling this function.
   bool rigid_body_fit(const coot::minimol::molecule &mol_without_moving_zone,
		       const coot::minimol::molecule &range_mol,
		       int imol_ref_map,
		       bool mask_water_flag);

   static short int in_residue_info_define; // initially 0
   static float geometry_vs_map_weight; // actually it's the other way around, isn't it? rename this.
   static float rama_plot_restraints_weight;
   static int  restraints_rama_type;
   static int rama_n_diffs;
   static double torsion_restraints_weight;
   static int refine_params_dialog_geman_mcclure_alpha_combobox_position;
   static int refine_params_dialog_lennard_jones_epsilon_combobox_position;
   static int refine_params_dialog_rama_restraints_weight_combobox_position;
   static int refine_params_dialog_torsions_weight_combox_position;
   static bool refine_params_dialog_extra_control_frame_is_visible;

   // similarly for distance and angles:
   //
   static short int in_distance_define;
   static short int in_angle_define;
   static short int in_torsion_define;

   void set_in_distance_define() { in_distance_define = 1; }
   void set_in_angle_define()    { in_angle_define = 1; }

   static short int in_db_main_define;
   static int db_main_imol;
   static int db_main_atom_index_1;
   static int db_main_atom_index_2;
   int load_db_main();

   static short int in_reverse_direction_define;

   // wrappers for regularization and refinement:
   //
   static short int do_torsion_restraints; // all, including side chains
   static short int do_peptide_omega_torsion_restraints;
   static bool do_rama_restraints;
   static bool do_trans_peptide_restraints;
   static bool do_numerical_gradients; // for debugging
   static bool do_rotamer_restraints;
   static bool do_debug_refinement;

   std::vector<std::pair<mmdb::Residue *, std::vector<coot::dict_torsion_restraint_t> > > make_rotamer_torsions(const std::vector<std::pair<bool, mmdb::Residue *> > &local_residues) const;

   coot::refinement_results_t regularize(int imol, short int auto_range_flag, int i_atom_start, int i_atom_end);
   coot::refinement_results_t refine    (int imol, short int auto_range_flag, int i_atom_start, int i_atom_end);
   // a more modern interface to refine:
   coot::refinement_results_t refine_residue_range(int imol,
						   const std::string &chain_id1,
						   const std::string &chain_id2,
						   int resno_1,
						   const std::string &ins_code_1,
						   int resno_2,
						   const std::string &ins_code_2,
						   const std::string &altconf,
						   short int is_water_flag);
   coot::refinement_results_t get_refinement_results() const;

   static void flash_position(const clipper::Coord_orth &pos);

   void repeat_refine_zone(); // no interesting return value because it uses refine()
                              // as in check_if_in_refine_define().
   std::pair<int, int> auto_range_residues(int atom_index, int imol) const;


   // Idealize the geometry without considering the map.
   //
   // return 1 if restraints were found, 0 if not.
   //
   coot::refinement_results_t
     copy_mol_and_regularize(int imol,
			     int resno_1,
			     std::string inscode_1,
			     int resno_2,
			     std::string inscode_2,
			     std::string altconf, // use this altconf or "" atoms.
			     std::string chain_id_1);

   // Regularize *and* fit to density.
   //
   // return 1 if restraints were found, 0 if not.
   //
   coot::refinement_results_t
   copy_mol_and_refine(int imol_for_atoms,
		       int imol_for_map,
		       int resno_1,
		       std::string inscode_1,
		       int resno_2,
		       std::string inscode_2,
		       std::string altconf, // use this altconf or "" atoms.
		       std::string chain_id_1);

   // no refinement, simple copy, can return -1 for a problem.
   int copy_model_molecule(int imol);

   bool check_for_no_restraints_object(std::string &resname_1, std::string &resname_2) const;
   bool check_for_single_hetatom(mmdb::Residue *res_p) const;


   // void (int imol,
   // const std::vector<mmdb::Residue *> &residues);


   // return 0 if any of the residues in selection don't have (at least) bond
   // restraints.  Try to auto-load the dictionary cifs and try again.
   // The vector is a list of residues for which no restraints could be found.
   std::pair<int, std::vector<std::string> >
     check_dictionary_for_residue_restraints(int imol, mmdb::PResidue *SelResidues, int nSelResidues);
   std::pair<int, std::vector<std::string> >
     check_dictionary_for_residue_restraints(int imol, const std::vector<mmdb::Residue *> &residues);

   static void show_missing_refinement_residues_dialog(const std::vector<std::string> &res_names,
						       bool run_get_monomer_post_fetch_flag);

   // called by copy_mol_and_refine and copy_mol_and_regularize
   //
   mmdb::Manager *create_mmdbmanager_from_res_selection(mmdb::PResidue *SelResidues,
						       int nSelResidues,
						       int have_flanking_residue_at_start,
						       int have_flanking_residue_at_end,
						       const std::string &altconf,
						       const std::string &chain_id_1,
						       short int residue_from_alt_conf_split_flag,
						       int imol); // imol is for uddatom index

   // called by simple_refine_residues (a refinement from a vector of mmdb::Residues).
   //
   // The returned mol should have flanking residues too.
   //
   // return also a vector of residues that correspond to the residues
   // that were input - the non-fixed residues.
   //
   std::pair<mmdb::Manager *, std::vector<mmdb::Residue *> >
   create_mmdbmanager_from_res_vector(const std::vector<mmdb::Residue *> &residues,
				      int imol, // for uddatom index.
				      mmdb::Manager *mol,
				      std::string alt_conf);
   // which uses
   int find_serial_number_for_insert(int seqnum_new,
				     const std::string &ins_code,
				     mmdb::Chain *chain_p) const;

   // simple mmdb::Residue * interface to refinement.  20081216
   coot::refinement_results_t
     generate_molecule_and_refine(int imol,  // needed for UDD Atom handle transfer
				  const std::vector<mmdb::Residue *> &residues,
				  const std::string &alt_conf,
				  mmdb::Manager *mol,
				  bool use_map_flag);
   coot::refinement_results_t
     generate_molecule_from_molecule_and_refine(int imol, mmdb::Manager *mol, bool use_map_flag);

   coot::refinement_results_t
     refine_residues_vec(int imol,
			 const std::vector<mmdb::Residue *> &residues,
			 const std::string &alt_conf,
			 mmdb::Manager *mol);
   coot::refinement_results_t
     regularize_residues_vec(int imol,
			     const std::vector<mmdb::Residue *> &residues,
			     const std::string &alt_conf,
			     mmdb::Manager *mol);
   coot::refinement_results_t triple_refine_auto_accept();

   coot::refinement_results_t refine_molecule(int imol, mmdb::Manager *mol);
   coot::refinement_results_t refine_chain(int imol, const std::string &chain_id, mmdb::Manager *mol);

   static bool use_harmonic_approximation_for_NBCs;
   void set_use_harmonic_approximations_for_nbcs(bool flag);

   // on reading a pdb file, we get a list of residues, use these to
   // load monomers from the dictionary
   int load_needed_monomers(const std::vector<std::string> &pdb_residue_types);


   // geometry graphs
   void update_geometry_graphs(const atom_selection_container_t &asc, int imol_moving_atoms);
   void update_geometry_graphs(int imol_moving_atoms); // convenience function - includes sequence view too!

   void update_validation(int imol);
   // which wraps:
   void update_validation_graphs(int imol);
   // 20211201-PE currently update_geometry
   void update_ramachandran_plot(int imol);


   // Display the graphical object of the regularization

   static void draw_moving_atoms_graphics_object(bool against_a_dark_background);
   static void draw_ramachandran_goodness_spots();
   static void draw_rotamer_probability_object();
   static void draw_moving_atoms_peptide_markup();
   static void draw_moving_atoms_atoms(bool against_a_dark_background);
   static void draw_moving_atoms_restraints_graphics_object();
   std::vector<coot::old_generic_display_object_t::dodec_t> get_rotamer_dodecs();

   static int mol_no_for_environment_distances;
   static bool display_environment_graphics_object_as_solid_flag;
   static void draw_environment_graphics_object(stereo_eye_t eye);
   // void symmetry_environment_graphics_object() const;

   // for flashing the picked intermediate atom.
   static void picked_intermediate_atom_graphics_object();

   void update_environment_distances_maybe(int index, int imol);
   void update_environment_distances_by_rotation_centre_maybe(int imol);

   // cif dictionary read number.  Update on reading (or attempting to
   // read) a cif dictionary file.  Public because reading cif file in
   // molecule-class-info-other.cc (bad chiral volumes) needs it.
   static int cif_dictionary_read_number;

   static int residue_selection_flash_frames_number;

   // scripting
   static short int guile_gui_loaded_flag;
   static short int python_gui_loaded_flag;
   static std::vector<std::string> command_line_scripts;
   static coot::command_line_commands_t command_line_commands;
   static std::vector<std::string> command_line_accession_codes;

   // background colour
   static glm::vec3 background_colour;
   static glm::vec4 get_background_colour();
   static bool background_is_black_p();

   // dynarama: a list of dynarama canvases, each of which has
   // attached a pointer to a rama_plot class object.
   //
   static GtkWidget **dynarama_is_displayed;
   void set_dynarama_is_displayed(GtkWidget *dyna_canvas, int imol);
   void destroy_edit_backbone_rama_plot(); // only one of these.

   // sequence view
   static GtkWidget **sequence_view_is_displayed;
   void set_sequence_view_is_displayed(GtkWidget *seq_view_canvas, int imol);
   GtkWidget * get_sequence_view_is_displayed(int imol) const;
   static int nsv_canvas_pixel_limit;
   void sequence_view_highlight_residue_maybe(mmdb::Atom *next_atom, GtkWidget *svc);

   // Geometry Graphs:

/* Old style array of molecules code */
/*    static GtkWidget **geometry_graph; */
/*    static GtkWidget **b_factor_variance_graph; */
/*    static GtkWidget **residue_density_fit_graph; */
/*    static GtkWidget **omega_distortion_graph; */
/*    static GtkWidget **rotamer_graph; */
/*    static GtkWidget **ncs_diffs_graph; */

#ifdef DO_GEOMETRY_GRAPHS
   std::vector<coot::geometry_graph_block_info_generic>
   density_fit_from_mol(const atom_selection_container_t &asc, int imol_moving_atoms,
			int imol_for_map);
   std::vector<coot::geometry_graph_block_info_generic>
   density_fit_from_residues(mmdb::PResidue *SelResidues, int nSelResidues,
			     int imol_moving_atoms,
			     int imol_for_map) const;
#endif

   coot::omega_distortion_info_container_t
     omega_distortions_from_mol(const atom_selection_container_t &asc, const std::string &chain_id);

#ifdef DO_GEOMETRY_GRAPHS
   std::vector<coot::geometry_graph_block_info_generic>
     rotamers_from_mol(const atom_selection_container_t &asc, int imol_moving_atoms);
   std::vector<coot::geometry_graph_block_info_generic>
     rotamers_from_residue_selection(mmdb::PResidue *SelResidues,
				   int nSelResidues, int imol);

   std::vector<coot::geometry_graph_block_info_generic> ncs_diffs_from_mol(int imol);
   std::vector<coot::geometry_graph_block_info_generic> ncs_diffs(int imol, const coot::ncs_chain_difference_t &d);

#endif

   // now used for rotamer_score (from c-interface.h), so it is now not GTK2-only 20090817
   coot::rotamer_probability_info_t get_rotamer_probability(mmdb::Residue *res,
							    const std::string &alt_conf,
							    mmdb::Manager *mol,
							    float lowest_probability,
							    short int add_extra_PHE_and_TYR_rotamers_flag);

   static float residue_density_fit_scale_factor; // 1.0 usually, adjustable for CNS/mapman map users.

   // shall we recentre when a new molecule is added (default no)
   static short int recentre_on_read_pdb;

   // map radius slider maximum
   static float map_radius_slider_max;

   // mouse buttons
   GdkModifierType gdk_button1_mask();
   GdkModifierType gdk_button2_mask();
   GdkModifierType gdk_button3_mask();

   void quanta_buttons();

   static int draw_axes_flag;
   //
   static short int display_density_level_on_screen;
   static short int display_density_level_this_image;
   static std::string display_density_level_screen_string;
   void set_density_level_string(int imol, float dlevel);

   static coot::Cartesian to_cartesian(const clipper::Coord_orth &co);
   static clipper::Coord_orth to_coord_orth(const coot::Cartesian &c);

   static std::string int_to_string(int i);
   static std::string float_to_string(float f);
   static std::string float_to_string_using_dec_pl(float f, unsigned short n_dec_pl);
   static std::string backslash_filename(const std::string &s); // needed for windows?

   void set_find_ligands_mols(int map, int protein,
			      const std::vector<std::pair<int, bool> > &ligand_wiggly_info);
   int find_ligand_map_mol() const;
   int find_ligand_protein_mol() const;
   static void set_ligand_protein_mol(int imol);
   static void set_ligand_map_mol(int imol);
   static void find_ligand_add_rigid_ligand(int imol);
   static void find_ligand_add_flexible_ligand(int imol);
   void set_find_ligand_do_real_space_refine_state(bool state);
   bool find_ligand_do_real_space_refine_state();
   std::vector<std::pair<int, bool> > find_ligand_ligand_mols() const;
   void find_ligand_clear_ligand_mols();

   static std::vector<clipper::Coord_orth> *ligand_big_blobs;
   static int find_ligand_n_top_ligands;
   static short int find_ligand_mask_waters_flag;
   static int find_ligand_ligand_atom_limit;
   static short int ligand_expert_flag;
   static bool find_ligand_here_cluster_flag;
   static float map_mask_atom_radius;
   static bool find_ligand_multiple_solutions_per_cluster_flag; // default is false
   static float find_ligand_score_by_correl_frac_limit; // 0.7
   static float find_ligand_score_correl_frac_interesting_limit; // 0.9;

   static void rebond_molecule_corresponding_to_moving_atoms();

   // Geometry issues:

   static coot::protein_geometry *Geom_p();

   //
   std::vector <coot::dict_torsion_restraint_t> get_monomer_torsions_from_geometry(const std::string &monomer_type) const;

   // make regularize_object_bonds_box be the empty box (and redisplay)
   //
   void clear_moving_atoms_object();
   // copy the contents of moving_atoms_asc into the molecule being refined.
   //
   coot::refinement_results_t accept_moving_atoms();

   void update_moving_atoms_from_molecule_atoms(const coot::minimol::molecule &mm);

   void set_refinement_map(int imol);

   // public access to the clear the in range defines
   void clear_pending_picks();


   // For updating the goto atom widget, let's store the window
   // pointer here, so that on recentering (e.g. on "space"), we check
   // to see if this variable is non-null, if it is, we use it to
   // update the vaules in the entries.
   //
   // Also, this is needed for deleting molecules (menuitems) from the
   // menu when molecules are closed.
   //
   // FIXME Needs to be implemented.  close window button callback and
   // destroy window need to set this to NULL.
   //

   // 20220723-PE has (will have) internals blanked by EMSCRIPTEN
   void model_fit_refine_unactive_togglebutton(const std::string &button_name) const;

   static GtkWidget *go_to_atom_window;


   // For updating the model/fit/refine widget (turning on/off
   // checkbuttons), let's store the Model/Fit/Refine window and
   // create only one of them.
   //
   static GtkWidget *model_fit_refine_dialog;
   static GtkWidget *other_modelling_tools_dialog;

   // And similarly for the rotamers, we now want to be able to change
   // them via keyboard arrow keys, so we need to register the rotamer
   // dialog here, the number of rotamers will be attached to the
   // dialog and the dialog will get sent :next rotamer:, :previous
   // rotamer: signals.

   // static GtkWidget *rotamer_dialog; // 20240304-PE goodbye, these days, look it up.

   // And also for the difference map peaks dialog, which people want
   // to scroll through using . and ,

   static GtkWidget *difference_map_peaks_dialog;
   static float difference_map_peaks_max_closeness;

   //
   void other_modelling_tools_unactive_togglebutton(const std::string &button_name) const;


   // And also for the residue info button, we don't want to mess
   // around with having lots of these, life is too short to sort out
   // that mess.
   static GtkWidget *residue_info_dialog;


   // geometry:
   // (text message in the console, currently)

   // not const because we add to distance_object_vec
   // return the distance
   float add_measure_distance(const coot::Cartesian &p1,
                              const coot::Cartesian &p2);
   static std::vector<atom_label_info_t> labels_for_measure_distances_and_angles;
   void add_measure_angle() const; // uses class variables
   double get_geometry_torsion() const;
   void display_geometry_torsion() const;
   // return the distance
   // double display_geometry_distance_symm(int imol1, const coot::Cartesian &p1,
   //                                       int imol2, const coot::Cartesian &p2);
   //
   void pepflip(); // 20220814-PE old style flip interface
   void pepflip(int imol, const coot::atom_spec_t &spec);

   // return true if moving_atoms_asc was not null (more or less if
   // the pepflip was made)
   //
   bool pepflip_intermediate_atoms();
   bool pepflip_intermediate_atoms_other_peptide();
   bool pepflip_intermediate_atoms(mmdb::Atom *at_close);

   // return true if moving_atoms_asc was not null (more or less if
   // the rotamer fit was made)
   //
   bool backrub_rotamer_intermediate_atoms();

   // we need this in c-interface.cc for the rigid body refinement
   // which refines against a map
   //
   // return -1 on no map found.
   //
   // This maybe sets imol_refinement_map (if it was -1 or if the map
   // it was previously is now closed).
   int Imol_Refinement_Map() const;
   //
   int set_imol_refinement_map(int imol);
   float get_estimated_map_weight(int imol_map);

   enum moving_atoms_do_rama_markup { MOVING_ATOMS_DO_RAMA_MARKUP_FALSE, MOVING_ATOMS_DO_RAMA_MARKUP_TRUE, MOVING_ATOMS_DO_RAMA_MARKUP_USE_INTERNAL_SETTING};
   enum moving_atoms_do_rota_markup { MOVING_ATOMS_DO_ROTA_MARKUP_FALSE, MOVING_ATOMS_DO_ROTA_MARKUP_TRUE, MOVING_ATOMS_DO_ROTA_MARKUP_USE_INTERNAL_SETTING};

   void make_moving_atoms_graphics_object(int imol, const atom_selection_container_t &asc,
                                          unsigned int do_rama_markup=MOVING_ATOMS_DO_RAMA_MARKUP_USE_INTERNAL_SETTING,
                                          unsigned int do_rota_markup=MOVING_ATOMS_DO_ROTA_MARKUP_USE_INTERNAL_SETTING);
   static short int moving_atoms_asc_type;
   void make_moving_atoms_restraints_graphics_object();
   static coot::extra_restraints_representation_t moving_atoms_extra_restraints_representation;
   static bool draw_it_for_moving_atoms_restraints_graphics_object;
   static bool draw_it_for_moving_atoms_restraints_graphics_object_user_control;
   static bool draw_missing_loops_flag;

   //
   static float environment_min_distance;
   static float environment_max_distance;
   static bool environment_show_distances;
   static bool environment_distances_show_bumps;
   static bool environment_distances_show_h_bonds;
   static short int environment_distance_label_atom;

   // private?
   static std::set<mmdb::Residue *> moving_atoms_visited_residues; // visited during an rsr. Reset on "make_restraints"
   static mmdb::Atom *active_atom_for_hud_geometry_bar;
   void update_environment_graphics_object(int atom_index, int imol);
   void update_symmetry_environment_graphics_object(int atom_index, int imol);
   void add_distance_labels_for_environment_distances();
   static std::vector<atom_label_info_t> labels;  // environment distances, maybe other things too.
   static TextureMesh tmesh_for_labels;
   static HUDMesh mesh_for_hud_geometry;
   static std::string label_for_hud_geometry_tooltip;
   static std::map<std::string, Texture> texture_for_hud_geometry_labels_map;     // image to texture for
   static HUDTextureMesh mesh_for_hud_geometry_labels; // labels for the bars
   static HUDTextureMesh mesh_for_hud_tooltip_background;
   static Texture texture_for_hud_tooltip_background;
   static HUDTextureMesh tmesh_for_hud_geometry_tooltip_label;
   static HUDTextureMesh tmesh_for_hud_image_testing;
   static Shader shader_for_hud_geometry_tooltip_text; // shader for the above tmesh (not like atom labels
                                                       // HUD labels are in 2D, don't need mvp, eye position
                                                       // etc.).
   static HUDTextureMesh tmesh_for_hud_refinement_dialog_arrow;
   static Texture texture_for_hud_refinement_dialog_arrow;
   static Texture texture_for_hud_refinement_dialog_arrow_highlighted;
   static bool hud_refinement_dialog_arrow_is_moused_over; // maybe can be outside EMSCRIPTEN

   static void show_refinement_and_regularization_parameters_frame(); // 20211102-PE it's here because the dialog can be
                                                                // shown by pressing the green arrow during refinement.

   static float get_x_base_for_hud_geometry_bars();

   static Texture texture_for_camera_facing_quad; // debugging
   static TextureMesh tmesh_for_camera_facing_quad;
   static Shader camera_facing_quad_shader;  // uses camera-facing-quad-shader-for-testing.shader

   static Texture texture_for_background_image;
   static HUDTextureMesh tmesh_for_background_image; // a HUDTextureMesh because we don't need tanget, bitangent
   static Shader shader_for_background_image;
   static bool draw_background_image_flag;

   void add_label(const std::string &l, const glm::vec3 &p, const glm::vec4 &c);

   //
   static short int dynamic_map_resampling;
   static short int dynamic_map_size_display;
   static int graphics_sample_step;
   static int dynamic_map_zoom_offset;
   static float goodsell_chain_colour_wheel_rotation_step;

   // uses cif_dictionary_filename_vec.
   //imol_enc can be the model molecule number or
   // IMOL_ENC_ANY = -999999, IMOL_ENC_AUTO = -999998, IMOL_ENC_UNSET = -999997.
   //
   // @return the index of the monomer in the geometry store.
   //
   coot::read_refmac_mon_lib_info_t
   add_cif_dictionary(std::string cif_dictionary_filename,
			  int imol_enc,
			  short int show_no_bonds_dialog_maybe_flag);
   void import_all_refmac_cifs();
   static std::vector<std::string> *cif_dictionary_filename_vec;

   // make private?
   static short int python_history;
   static short int guile_history;
   void add_history_command(const std::vector<std::string> &command_strings);
   static coot::history_list_t history_list;
   static coot::command_history_t command_history;

   // this does not quote strings - it just copies out the arguments
   // "bare".  If the arguments are strings they should be quoted
   // before calling this.
   static std::string pythonize_command_strings(const std::vector<std::string> &command_strings);
   // likewise, no string quoting.
   static std::string schemize_command_strings(const std::vector<std::string> &command_strings);

   static coot::console_display_commands_t console_display_commands;

   // build one residue by phi psi
   static int add_terminal_residue_n_phi_psi_trials;
   static int add_terminal_residue_add_other_residue_flag;
   static short int add_terminal_residue_do_post_refine;
   static float terminal_residue_addition_direct_phi;
   static float terminal_residue_addition_direct_psi;
   static bool add_terminal_residue_debug_trials;
   // we allow terminal fitting without rigid body refinement
   static short int add_terminal_residue_do_rigid_body_refine;

   // return success status: 1 for success
   int execute_add_terminal_residue(int imol,
				     const std::string &terminus,
				     mmdb::Residue *residue_p,
				     const std::string &chain_id,
				     const std::string &res_type,
				     bool immediate_addition_flag);
   void add_terminal_residue_using_active_atom(); // wraps above
   void execute_simple_nucleotide_addition(int imol, const std::string &term_type,
                                           mmdb::Residue *res_p, const std::string &chain_id);
   void execute_simple_nucleotide_addition(int imol, const std::string &chain_id, int res_no);

   static short int add_terminal_residue_immediate_addition_flag;
   static short int refinement_immediate_replacement_flag;  // don't dialog me please
   // called by above (private)
   atom_selection_container_t add_side_chain_to_terminal_res(atom_selection_container_t asc,
							     const std::string &res_type,
							     const std::string &terminus_type,
							     bool add_other_residue_flag);


   // public (from globjects);
   //
   void execute_db_main();
   // direction string (flag): "forwards" or "backwards" (actually, if
   // not "backwards", "forwards" is presumed).
   // return the new molecule number (can be -1)
   int execute_db_main(int imol, std::string chain_id,
		       int iresno_start, int iresno_end, std::string direction_string);
   std::pair<int, int> execute_db_main_fragment(int imol, coot::residue_spec_t spec); // build both directions.

   static float ligand_acceptable_fit_fraction;
   static float ligand_cluster_sigma_level; // was 2.2 default
   static int   ligand_wiggly_ligand_n_samples;
   static int   ligand_wiggly_ligand_count;
   static int   ligand_verbose_reporting_flag;

   // Eleanor wants control over water parameters
   //
   static float ligand_water_to_protein_distance_lim_max;
   static float ligand_water_to_protein_distance_lim_min;
   static float ligand_water_variance_limit;
   static int   ligand_water_n_cycles;
   static short int ligand_water_write_peaksearched_atoms;

   // Delete item mode: (this is a terrible way of doing it). Use one variable!
   static short int delete_item_atom;
   static short int delete_item_residue;
   static short int delete_item_residue_zone;
   static short int delete_item_residue_hydrogens;
   static short int delete_item_water;
   static short int delete_item_sidechain;
   static short int delete_item_sidechain_range;
   static short int delete_item_chain;
   // must save the widget so that it can be deleted when the item is selected.

   static GtkWidget *delete_item_widget; // 20220723-PE do I need this now?
   static int keep_delete_item_active_flag;
   // really, we should save pick data with atom or residue specs, so
   // let's start with delete zones' first click:
   static int delete_item_residue_zone_1_imol;
   static int delete_item_sidechain_range_1_imol;
   static coot::residue_spec_t delete_item_residue_zone_1;
   static coot::residue_spec_t delete_item_sidechain_range_1;
   void delete_residue_range(int imol, const coot::residue_spec_t &res1,
			     const coot::residue_spec_t &res2);
   void delete_sidechain_range(int imol,
			       const coot::residue_spec_t &res_1,
			       const coot::residue_spec_t &res_2);
   void delete_active_residue();
   void delete_active_residue_alt_conf_atoms();
   // c-info functions really, but we cant have mmdb_manager there, so the are moved here.

   static void output_residue_info_as_text(int atom_index, int imol);

   static void output_residue_info_dialog(int imol, int atom_index);

   static void output_residue_info_dialog(int imol, const coot::residue_spec_t &rs);
   //

   static void fill_output_residue_info_widget(GtkWidget *widget, int imol,
					       std::string residue_name,
					       mmdb::PPAtom atoms, int n_atoms);
   static void fill_output_residue_info_widget_atom(GtkWidget *dialog, GtkWidget *widget,
						    int imol, mmdb::PAtom atom, int iat);
   // and the keypress callbacks for the above
#if (GTK_MAJOR_VERSION >= 4)
#else
   static gboolean on_residue_info_occ_entry_key_release_event (GtkWidget       *widget,
								GdkEventKey     *event,
								gpointer         user_data);
   static gboolean on_residue_info_master_atom_occ_changed (GtkWidget       *widget,
							    GdkEventKey     *event,
							    gpointer         user_data);

   static gboolean on_residue_info_master_atom_b_factor_changed (GtkWidget       *widget,
								 GdkEventKey     *event,
								 gpointer         user_data);
#endif

   // Return the molecule number of the selected map (I mean, top of
   // the list, in the option menu)

   static void refinement_map_combobox_changed(GtkWidget *c, gpointer data);
   int fill_combobox_with_map_mtz_options(GtkWidget *combobox, GCallback signal_func, int imol);

   void fill_combobox_with_molecule_options(GtkWidget *combobox,
					    GCallback signal_func,
					    int imol_active_position,
					    const std::vector<int> &molecules_index_vec);

   // return the imol of the active item molecule
   int fill_combobox_with_map_options(GtkWidget *combobox,
				      GCallback signal_func,
				      int imol_active_position);

   void fill_combobox_with_difference_map_options(GtkWidget *combobox,
						  GCallback signal_func,
						  int imol_active_position);

   GtkWidget *wrapped_create_skeleton_dialog(bool show_ca_mode_needs_skel_label);
   void skeletonize_map_by_optionmenu(GtkWidget *optionmenu);
   void skeletonize_map_by_combobox(GtkWidget *combobox);

   static void on_skeleton_ok_button_dynamic_clicked (GtkButton       *button,
						      gpointer         user_data);
   int try_set_draw_baton(short int i);


   void fill_combobox_with_skeleton_options(GtkWidget *combobox);

   void set_on_off_skeleton_radio_buttons(GtkWidget *skeleton_frame);
   void set_on_off_single_map_skeleton_radio_buttons(GtkWidget *skeleton_frame,
						     int i);
   void set_contour_sigma_button_and_entry(GtkWidget *window, int imol);

   void fill_combobox_with_refmac_methods_options(GtkWidget *combobox);
   void fill_combobox_with_refmac_phase_input_options(GtkWidget *combobox);
   void fill_combobox_with_refmac_labels_options(GtkWidget *combobox);
   void fill_combobox_with_refmac_mtz_file_options(GtkWidget *combobox);
   void fill_combobox_with_refmac_file_labels_options(GtkWidget *combobox);
   void fill_combobox_with_refmac_ncycles_options(GtkWidget *combobox);

   void update_refmac_column_labels_frame(GtkWidget *optionmenu,
					  GtkWidget *fobs_menu, GtkWidget *fiobs_menu, GtkWidget *fpm_menu,
					  GtkWidget *r_free_menu,
					  GtkWidget *phases_menu, GtkWidget *fom_menu, GtkWidget *hl_menu);
   static void refinement_map_select(GtkWidget *item, GtkPositionType pos);
   static void refinement_map_select_add_columns(GtkWidget *item, GtkPositionType pos);
   static void select_refinement_map_combobox_changed(GtkWidget *combobox, gpointer data);
   // static void   skeleton_map_select(GtkWidget *item, GtkPositionType pos);
   static void skeleton_map_combobox_changed(GtkWidget *combobox, gpointer data);

   static int map_for_skeletonize; // used by skeletonize_map;
   static void   skeletonize_map(int imol, short int prune_flag);
   static void unskeletonize_map(int imol);
   static void set_initial_map_for_skeletonize();
   static bool add_ccp4i_projects_to_optionmenu_flag;
   static std::string refmac_ccp4i_project_dir;
   static std::string libcheck_ccp4i_project_dir;

   // refmac stuff
   static std::vector<int> *preset_number_refmac_cycles;
   static coot::refmac::refmac_refinement_method_type refmac_refinement_method;
   static coot::refmac::refmac_phase_input_type refmac_phase_input;
   static coot::refmac::refmac_use_tls_type     refmac_use_tls_flag;
   static coot::refmac::refmac_use_twin_type    refmac_use_twin_flag;
   static coot::refmac::refmac_use_sad_type     refmac_use_sad_flag;
   static coot::refmac::refmac_use_ncs_type     refmac_use_ncs_flag;
   static coot::refmac::refmac_use_intensities_type refmac_use_intensities_flag;
   static coot::refmac::refmac_used_mtz_file_type refmac_used_mtz_file_flag;
   static const char *saved_refmac_file_filename;
   static int refmac_ncycles;
   static void set_refmac_refinement_method(int method);
   static void refmac_change_refinement_method(GtkWidget *item, GtkPositionType pos);
   static void refmac_refinement_method_combobox_changed(GtkWidget *combobox, gpointer data);
   static void set_refmac_phase_input(int phase_flag);
   static void refmac_change_phase_input(GtkWidget *item, GtkPositionType pos);
   static void refmac_refinement_phase_info_combobox_changed(GtkWidget *combobox, gpointer data);
   static void set_refmac_use_tls(int state);
   static void set_refmac_use_twin(int state);
   static void set_refmac_use_sad(int state);
   static void set_refmac_use_intensities(int state);
   static void set_refmac_used_mtz_file(int state);
   static void set_refmac_n_cycles(int no_cycles);
   static void set_refmac_use_ncs(int state);
   static std::vector<coot::refmac::sad_atom_info_t> refmac_sad_atoms;
   static void refmac_change_ncycles(GtkWidget *item, GtkPositionType pos);
   static GtkWidget *refmac_dialog_mtz_file_label;
   void add_refmac_sad_atom(const char *atom_name, float fp, float fpp, float lambda);
   void add_refmac_ncycle_no(int &cycle);
   static short int have_sensible_refmac_params;
   static std::string refmac_mtz_file_filename;
   static std::string refmac_fobs_col;
   static std::string refmac_sigfobs_col;
   static std::string refmac_r_free_col;
   static int refmac_r_free_flag_sensible;

   std::string Refmac_mtz_file_filename() const { return refmac_mtz_file_filename; }
   std::string Refmac_fobs_col() const { return refmac_fobs_col; }
   std::string Refmac_sigfobs_col() const { return refmac_sigfobs_col; }
   std::string Refmac_r_free_col() const { return refmac_r_free_col; }
   short int Refmac_r_free_sensible() const { return refmac_r_free_flag_sensible; }
   void store_refmac_params(const std::string &mtz_filename,
                            const std::string &fobs_col,
                            const std::string &sigfobs_col,
                            const std::string &r_free_col,
                            int r_free_flag_sensible);
   void store_refmac_phase_params(const std::string &phi,
                                  const std::string &fom,
                                  const std::string &hla,
                                  const std::string &hlb,
                                  const std::string &hlc,
                                  const std::string &hld);
   void store_refmac_file_mtz_filename(const std::string &mtz_filename);

   static int max_skeleton_search_depth;
   // Rotamer stuff
   static int rotamer_search_mode;
   static float rotamer_lowest_probability;  // Expressed as a percentage (2.0 default).
                                      // This is the P_r1_r2_r3_r4 probability, no bayesian
                                      // stuff involved, currently.
   static float rotamer_distortion_scale; // for the validation graphs
   static int rotamer_fit_clash_flag;
   static short int in_rotamer_define;
   static int rotamer_residue_atom_index;  // not used because we use the atom spec now
   static coot::atom_spec_t rotamer_residue_atom_spec;
   static int rotamer_residue_imol;

   // 20220812-PE old interface
   void do_rotamers(int atom_index, int imol) ; // display the rotamer option and display
                         			// the most likely in the graphics as a
			                        // moving_atoms_asc
   // 20220812-PE new interface
   void do_rotamers(int imol, mmdb::Atom *active_atom);

   // 20220812-PE void fill_rotamer_selection_buttons(GtkWidget *window, int atom_index, int imol) const;
   void fill_rotamer_selection_buttons(GtkWidget *window, mmdb::Atom *atom, int imol) const;

   short int generate_moving_atoms_from_rotamer(int imol, coot::atom_spec_t &spec, int irot);
   static void on_rotamer_selection_button_toggled(GtkCheckButton       *button,
                                                   gpointer         user_data);
   void set_rotamer_fit_clash_flag(int i) { rotamer_fit_clash_flag = i; }
   // autofit rotamer:
   static short int in_auto_fit_define;

   void auto_fit_rotamer_ng(int imol, const coot::residue_spec_t &res_spec, const std::string &alt_conf);

   // Mutation stuff
   static short int in_mutate_define;
   static int mutate_residue_atom_index;
   static int mutate_residue_imol;
   void do_mutation(int imol, const coot::residue_spec_t &res_spec,
                    const std::string &residue_type, short int do_stub_flag);
   static atom_selection_container_t standard_residues_asc;
   GtkWidget *wrapped_create_residue_type_chooser_window(bool show_stub_option_flag) const; /* set the stub checkbox */

   // mutate then auto fit:
   static short int in_mutate_auto_fit_define;
   static int mutate_auto_fit_residue_atom_index;
   static int mutate_auto_fit_residue_imol;
   void do_mutation_auto_fit(int imol, const coot::residue_spec_t &res_spec, const std::string &residue_type, short int do_stub_flag);
   static short int residue_type_chooser_auto_fit_flag;
   static short int residue_type_chooser_stub_flag;
   static short int mutate_auto_fit_do_post_refine_flag;
   static short int rotamer_auto_fit_do_post_refine_flag;

   // mutate sequence
   static std::string mutate_sequence_chain_from_combobox;
   static int         mutate_sequence_imol;

   // align and mutate
   static int         align_and_mutate_imol;
   // static std::string align_and_mutate_chain_from_optionmenu;
   static std::string align_and_mutate_chain_from_combobox;
   static mmdb::realtype alignment_wgap;
   static mmdb::realtype alignment_wspace;

   void mutate_chain(int imol, const std::string &chain_id,
		     const std::string &seq, bool do_auto_fit_flag,
		     bool renumber_residues_flag);

   // save molecule [option menu usage]
   static int save_imol;

   // --- Pointer Distances ---
   static float pointer_min_dist;
   static float pointer_max_dist;
   static bool show_pointer_distances_flag;
   void clear_pointer_distances();
   static std::vector<std::pair<clipper::Coord_orth, clipper::Coord_orth> > pointer_distances_object_vec;
   // static Mesh mesh_for_pointer_distances; // here for future-Paul
   static void draw_pointer_distances_objects(stereo_eye_t eye); // draw them
   void make_pointer_distance_objects(); // (re)generate them
   static std::vector<atom_label_info_t> labels_for_pointer_distances;

   // --- Extra Distance Restraints ---

   static void draw_extra_distance_restraints(stereo_eye_t eye, int pass_type);
   void make_extra_distance_restraints_objects(); // (re)generate them
   static float extra_distance_restraint_penalty_cutoff; // restraints that have less penalty/energy than
                                                         // this are not worth drawing.
   void set_extra_distance_restraint_penalty_cutoff(float c) {
      extra_distance_restraint_penalty_cutoff = c;
      make_extra_distance_restraints_objects();
   }
   static bool show_extra_distance_restraints_flag;
   void set_show_extra_distance_restriants(bool s) { show_extra_distance_restraints_flag = s; }
   static std::vector<extra_distance_restraint_markup_instancing_data_t> extra_distance_restraints_markup_data;
   static Mesh mesh_for_extra_distance_restraints; // draw this with instancing

   // --- Dynamic distances to intermediate atoms: ---

   static short int in_dynamic_distance_define;
   static coot::intermediate_atom_distance_t running_dynamic_distance;
   static std::vector<coot::intermediate_atom_distance_t> dynamic_distances;

   // Pointer atoms
   static short int pointer_atom_is_dummy; // force dummy atom, no atom type choice.
   void place_dummy_atom_at_pointer();
   void place_typed_atom_at_pointer(const std::string &type);
   int create_pointer_atom_molecule_maybe() const; // create a pointer atoms molecule if it
                                                   // does not exist, otherwise return the
                                                   // molecule number of the already existing
                                                   // pointer atoms molecule.
   static int user_pointer_atom_molecule; // set to -1;

   // save state, return success status of writing
   int save_state();
   int save_state_file(const std::string &filename);
   int save_state_file(const std::string &filename, short int il);
   int save_history() const;
   std::vector<std::string> save_state_data_and_models(short int lang_flag) const;
   std::vector<std::string> save_state_data_and_models(const std::string &filename,
						       short int lang_flag) const;
   void quick_save();

   // and the generalization of that: - display the label overlay for 2 seconds
   static void ephemeral_overlay_label_from_id(const std::string &overlay_label_id);

   // and the generalization of that! Just pass the text of the ephemeral overlay label
   static void ephemeral_overlay_label(const std::string &overlay_label_text);

   static std::string save_state_file_name;

   // show citation?
   static short int show_citation_notice;

   static void info_dialog(const std::string &s, bool use_markup=false);

   // static GtkWidget *info_dialog(const std::string &s, bool use_markup=false);
   void fill_unsaved_changes_dialog(GtkWidget *dialog) const;
   // makes an info_dialog and writes text
   void info_dialog_and_text(const std::string &s, bool use_markup=false);

   // Return success status.
   //
   short int write_state(const std::vector<std::string> &commands,
			 const std::string &filename) const;
   // 20090914 (crazy) attempt to get Coot to build on 32-bit MacOSX.
   //
   short int write_state_c_mode(const std::vector<std::string> &commands,
				const std::string &filename) const;
   short int write_state_fstream_mode(const std::vector<std::string> &commands,
				      const std::string &filename) const;

   //
   int check_for_unsaved_changes() const; // in state

   // File selection directory saving:
   // (uses private directory_for_fileselection)
   //
   // 20220723-PE these should be deleted (but not while I'm doing EMSCRIPTEN things)
   void set_directory_for_fileselection(GtkWidget *fileselection) const;
   void save_directory_from_fileselection(const GtkWidget *fileselection);
   void save_directory_for_saving_from_fileselection(const GtkWidget *fileselection);
   void set_file_for_save_fileselection(GtkWidget *fileselection) const;

   // for file_chooser

   void set_directory_for_filechooser(GtkWidget *fileselection) const;
   void save_directory_from_filechooser(const GtkWidget *fileselection);
   void save_directory_for_saving_from_filechooser(const GtkWidget *fileselection);
   void set_file_for_save_filechooser(GtkWidget *fileselection) const;

   // saving temporary files (undo)
   //
   std::string save_molecule_dir(int imol) const;
   int apply_undo();
   int apply_redo();
   int apply_undo_or_redo(bool is_undo);
   void activate_redo_button();
   // used by option menu item callback which sets the molecule for undoing
   void set_undo_molecule_number(int i) { undo_molecule = i; }
   // undo_molecule_select uses set_undo_molecule_number()

   // who calls this?
   //   static void undo_molecule_select(GtkWidget *item, GtkPositionType pos);

   static void undo_molecule_combobox_changed(GtkWidget *c, gpointer data);
   void fill_combobox_with_undo_options(GtkWidget *option_menu); // not const


   int Undo_molecule(coot::undo_type) const; // return -2 on ambiguity, -1 on unset
			      // and a molecule number >=0 for no
			      // ambiguity (or undo_molecule has been
			      // set already).

   // save CONECT records?
   //
   static int write_conect_records_flag;

   bool check_if_moving_atom_pull(bool was_a_double_click); // and setup moving atom-drag if we are.

   bool check_if_hud_bar_clicked(double x, double y); // if true, set draw_hud_tooltip_flag
   std::pair<bool, mmdb::Atom *> check_if_moused_over_hud_bar(double mouse_x, double mouse_y);
   // most of the above function is comment, so put it here:
   std::pair<bool, mmdb::Atom *> check_if_hud_bar_moused_over_or_act_on_hud_bar_clicked(double mouse_x, double mouse_y, bool act_on_hit);

   bool check_if_hud_button_clicked(double x, double y);
   bool check_if_hud_button_moused_over(double x, double y, bool button_1_is_down);
   bool check_if_hud_button_moused_over_or_act_on_hit(double x, double y, bool act_on_hit, bool button_1_is_down);

   bool check_if_hud_rama_plot_clicked(double mouse_x, double mouse_y);

   void unset_moving_atoms_currently_dragged_atom_index() {
     moving_atoms_currently_dragged_atom_index = -1;
   }

   // baton stuff:
   //
   static short int draw_baton_flag;
   static short int baton_mode; // if set, rotation moves the baton, not the view
   static short int baton_tmp_atoms_to_new_molecule;
   static void draw_baton_object();
   // return a boolean, shall we really draw the baton or not (for
   // example, we don't want to do that if there is no skeletonized
   // map
   bool start_baton_here(); // modify baton root and tip.
   void accept_baton_position(); /* put an atom at the tip and move baton */
   int baton_build_atoms_molecule() const; // -1 on no such molecule
					   // with name "Baton Atoms"
   void rotate_baton(const double &x, const double &y);
   void toggle_baton_mode();  // on a "b" button press

   void baton_tip_try_another();
   void baton_tip_previous();
   void shorten_baton();
   void lengthen_baton();
   void baton_build_delete_last_residue();
   void set_max_skeleton_search_depth(int v) { max_skeleton_search_depth = v; }
   void set_baton_build_params(int istart_resno,
			      const char *chain_id,
			      const char *backwards);


   // save parameters for running refmac: not here.  Now in molecule class
//    void save_refmac_params(std::string fobs,
// 			   std::string sig_fobs, std::string r_free,
// 			   short int r_free_flag);
//    static std::string refmac_fobs_col;
//    static std::string refmac_sigfobs_col;
//    static std::string refmac_r_free_col;
//    short int refmac_r_free_flag; // is sensible r-free column?

   // bond thickness
   void set_bond_thickness(int imol, float thick);
   static int bond_thickness_intermediate_value; // not intermediate atoms
   // static void bond_parameters_molecule_menu_item_select(GtkWidget *item, GtkPositionType pos);
   static void bond_parameters_molecule_combobox_changed(GtkWidget *combobox, gpointer data);
   static void bond_parameters_bond_width_combobox_changed(GtkWidget *combobox, gpointer data);
   static void fill_bond_parameters_internals(GtkWidget *w, int imol);
   static void bond_width_item_select(GtkWidget *item, GtkPositionType pos);
   // colour map rotation adjustment change:
   static void bond_parameters_colour_rotation_adjustment_changed(GtkAdjustment *adj,
								  GtkWidget *window);
   void fill_bond_colours_dialog_internal(GtkWidget *w);
   // new style: each molecule has its own bond rotation value: Here
   // is where the individual colour rotation steps are changed:
   static void bonds_colour_rotation_adjustment_changed(GtkAdjustment *adj,
							GtkWidget *window);
   static int bond_parameters_molecule;
   static float bond_thickness_intermediate_atoms; // white atoms
   void set_bond_thickness_intermediate_atoms(float f);


   // residue info
   //
   // Why not use just the atom index, you might ask.  Because, from
   // scripting we will pass the atom attributes, not the atom index so
   // we don't have to pick an atom in the c-interface and then expand
   // it to the residue in graphics-info.
   //
   static short int   residue_info_pending_edit_b_factor;
   static short int   residue_info_pending_edit_occ;
   static int         residue_info_n_atoms; // so that we can release the memory and
					    // propogate the changes to the other widgets
   static std::vector<coot::select_atom_info> residue_info_edits;

   void reset_residue_info_edits() { // set residue_info_edits to zero elements
      residue_info_edits.clear();
   }

   void residue_info_release_memory(GtkWidget *dialog);

   // 20230218-PE webassembly merge: restoring these as they are in graphics-info-gui.cc
   //             why were they removed?
   static void  residue_info_add_b_factor_edit(coot::select_atom_info sai, float val);
   static void  residue_info_add_occ_edit(     coot::select_atom_info sai, float val);

   void apply_residue_info_changes(); // and clear the vector
   static void residue_info_edit_b_factor_apply_to_other_entries_maybe(GtkWidget *widget);
   static void residue_info_edit_occ_apply_to_other_entries_maybe(GtkWidget *widget);

   void add_picked_atom_info_to_status_bar(int imol, int atom_index);

   // crosshairs
   static short int draw_crosshairs_flag;
   void crosshairs_text() const;

   // environment_show_distances
   // Return imol = -1 if no (close) atoms found.
   //
   // index, imol
   std::pair<int, int> get_closest_atom() const;

   //

   static GdkCursor pick_cursor_index; // user setable

   static void pick_cursor_maybe();
   static void pick_cursor_real();
   static void normal_cursor();
   static void watch_cursor();
   static void fleur_cursor();


   // Alternate Conformation
   static short int alt_conf_split_type;
   static short int alt_conf_split_type_number();
   static short int in_add_alt_conf_define;
   static int add_alt_conf_atom_index;
   static int add_alt_conf_imol;
   std::pair<bool,std::string> split_residue(int imol, int atom_index);
   std::pair<bool,std::string> split_residue(int imol, const std::string &chain_id, int resno, const std::string &ins_code, const std::string &altconf);
   void split_residue_range(int imol, int index_1, int index2);
   // we add altconf because we want to triple-split a residue
   // sometimes (which we do by splitting one of the already existing
   // altconfs - we need to know which one).
   static float add_alt_conf_new_atoms_occupancy;
   static short int show_alt_conf_intermediate_atoms_flag;
   static GtkWidget *add_alt_conf_dialog;
   static void new_alt_conf_occ_adjustment_changed(GtkAdjustment *adj, gpointer user_data);

   // Backbone torsion
   //
   static short int in_backbone_torsion_define;
   void execute_setup_backbone_torsion_edit(int imol, int atom_index);

   // Edit Phi/Psi
   static short int in_edit_phi_psi_define; // set by button callback: setup_edit_phi_psi()
   static int edit_phi_psi_atom_index;
   static int edit_phi_psi_imol; // needed, yes (for the ramachandran widget)
   void edit_phi_psi();
   void execute_edit_phi_psi(int atom_index, int imol); // c.f. execute_rotate_translate_ready()
   void set_edit_phi_psi_to(double phi, double psi); // a callback
						     // from the
						     // ramachandran
						     // widget
   void rama_plot_for_single_phi_psi(int imol, int atom_index);
   void rama_plot_for_2_phi_psis(int imol, int atom_index);

   // Edit Chi
   static short int in_edit_chi_angles_define;
   static short int in_edit_chi_mode_flag; // c.f. in baton_mode
   static short int in_edit_chi_mode_view_rotate_mode;
   static bool      edit_chi_angles_reverse_fragment;
   static coot::atom_spec_t chi_angles_clicked_atom_spec;
   void execute_edit_chi_angles(int atom_index, int imol);
   enum edit_chi_edit_type { UNSET, EDIT_CHI, RESIDUE_PARTIAL_ALT_LOCS};
   int wrapped_create_edit_chi_angles_dialog(const std::string &res_type,
					     edit_chi_edit_type mode);
   // used by above:
   // (imol should be encoded into vbox - it isn't yet) // FIXME
   int fill_chi_angles_vbox(GtkWidget *vbox, std::string res_type, edit_chi_edit_type mode);
   static void clear_out_container(GtkWidget *vbox);
   static std::string chi_angle_alt_conf;

   // multi-residue torsion
   static bool in_multi_residue_torsion_mode;   // for rotating atoms (not view)
   static bool in_multi_residue_torsion_define; // for picking atoms
   static bool multi_residue_torsion_reverse_fragment_mode;
   static int multi_residue_torsion_picked_residues_imol;
   static std::pair<int, int> multi_residue_torsion_rotating_atom_index_pair;
   static std::vector<coot::residue_spec_t> multi_residue_torsion_picked_residue_specs;

   static int edit_chi_current_chi;

   // real values start at 1:
   static void on_change_current_chi_button_clicked(GtkButton *button,
						    gpointer user_data);
   static void on_change_current_chi_button_entered(GtkButton *button,
						    gpointer user_data);

   static short int moving_atoms_move_chis_flag;
   void setup_flash_bond_using_moving_atom_internal(int ibond);

   void setup_flash_bond(int imol, coot::residue_spec_t residue_spec, int i_bond);

   //! angle in degrees.
   void rotate_chi(double x, double y); // a 'callback' from the
					     // chi editting
   					     // function

   // Torsion General
   static short int in_edit_torsion_general_flag;
   void rotate_chi_torsion_general(double x, double y);

   void execute_torsion_general();

   // which bond are we rotating about in the ligand?

   void add_flash_bond(const std::pair<clipper::Coord_orth, clipper::Coord_orth> &p);
   static std::pair<clipper::Coord_orth, clipper::Coord_orth> flash_bond;
   static short int draw_chi_angle_flash_bond_flag;
   static void draw_chi_angles_flash_bond();

   // Tinker with the atom positions of residue
   // (used by rotate_chi)
   // We need to pass the asc for the mol because we need it for seekcontacts()
   //
   short int update_residue_by_chi_change(int imol,
					  mmdb::Residue *residue,
					  atom_selection_container_t &asc,
					  int chi, double diff);
   // temporary storage, during the change-over
   short int update_residue_by_chi_change_old(mmdb::Residue *residue,
					  atom_selection_container_t &asc,
					  int chi, double diff);
   // this can throw an std::runtime_error exception.
   std::pair<std::string, std::string> get_chi_atom_names(mmdb::Residue *residue,
							  const coot::dictionary_residue_restraints_t &rest,
							  int nth_chi) const;


   std::vector<std::vector<int> > get_contact_indices_from_restraints(mmdb::Residue *residue,
								      const atom_selection_container_t &asc,
								      short int is_regular_residue_flag) const;

   // Do 180 degree sidechain flip stuff
   //
   static short int in_180_degree_flip_define;

   // More moving backbone stuff:
   //
   void set_moving_atoms(atom_selection_container_t asc, int imol, int new_coords_type);
   //
   // changes moving atoms:
   //
   void change_peptide_carbonyl_by(double angle); // in degrees.
   void change_peptide_peptide_by(double angle); // in degress

   // button call backs for moving backbone (this is the alternative
   // method to using rottrans buttons method (pointer to info class
   // attached to rot/trans window))
   // These change private data:
   void set_backbone_torsion_peptide_button_start_pos(int ix, int iy);
   void change_peptide_peptide_by_current_button_pos(int ix, int iy);
   void set_backbone_torsion_carbonyl_button_start_pos(int ix, int iy);
   void change_peptide_carbonyl_by_current_button_pos(int ix, int iy);

   // sequence view:
   //
   // return NULL on failure (this molecule does not have a sequence
   // view)
   //

#ifdef DO_SEQUENCE_VIEW
   coot::sequence_view *get_sequence_view(int imol);
#endif
#ifdef DO_RAMA_PLOT
   static coot::rama_plot *edit_phi_psi_plot;
#endif

   // distances and angles displayed on screen
   // uses distance_objects vector
   static bool display_generic_objects_as_solid_flag;
   static void draw_geometry_objects();
   static void draw_dynamic_distances();
   static void draw_generic_objects(unsigned int pass_type, stereo_eye_t eye);
   static void draw_generic_objects_simple();
   static void draw_generic_objects_solid();
   static void draw_generic_text();
   static void draw_particles(stereo_eye_t eye);
   static void draw_molecules_atom_labels(stereo_eye_t eye);
   static void draw_boids(stereo_eye_t eye);
   static void draw_happy_face_residue_markers(stereo_eye_t eye);
   static void draw_anchored_atom_markers(stereo_eye_t eye);
   static void draw_unhappy_atom_markers(stereo_eye_t eye, unsigned int pass_type);
   static void draw_hydrogen_bonds_mesh(stereo_eye_t eye); // like boids
   void setup_draw_for_particles();
   void clear_measure_distances();
   void clear_last_measure_distance();
   static GtkWidget *geometry_dialog;
   void unset_geometry_dialog_distance_togglebutton();
   void unset_geometry_dialog_angle_togglebutton();
   void unset_geometry_dialog_torsion_togglebutton();
   void unset_geometry_dialog_dynamic_distance_togglebutton();
   static bool setup_draw_for_particles_semaphore;
   static bool particles_have_been_shown_already_for_this_round_flag;

   static std::vector<glm::vec3> happy_face_residue_marker_starting_positions; // based on residues
                                                                               // and filled by
                                                                               // setup_draw_for_happy_face_residue_markers()
                                                                               // which is called once per refinement.
   static bool do_tick_happy_face_residue_markers;
   static unsigned int draw_count_for_happy_face_residue_markers;
   const unsigned int draw_count_max_for_happy_face_residue_markers = 160;
   void setup_draw_for_happy_face_residue_markers_init(); // run this once to setup instancing buffer
   void setup_draw_for_happy_face_residue_markers(); // run this every time we want to see faces,
                                                     // it sets the start position of the faces.
   static Texture texture_for_happy_face_residue_marker;
   // likewise
   static std::vector<glm::vec3> anchored_atom_marker_texture_positions; // based on residues
                                                                         // and filled by
                                                                         // setup_draw_for_anchored_atom_markers()
                                                                          // which is called once per refinement.
   // move these inside EMSCRIPTEN?
   void setup_draw_for_anchored_atom_markers_init(); // run this once to setup instancing buffer
   void setup_draw_for_anchored_atom_markers();     // run this every time we want to anchored atoms
                                                    // it sets the start position of the textures
   static Texture texture_for_anchored_atom_markers;
   static TextureMesh tmesh_for_anchored_atom_markers;

   void setup_draw_for_unhappy_atom_markers();
   static Texture texture_for_unhappy_atom_markers;
   static TextureMesh tmesh_for_unhappy_atom_markers;
   static void add_unhappy_atom_marker(int imol, const coot::atom_spec_t &atom_spec);
   static void remove_all_unhappy_atom_markers();

   static std::vector<meshed_particle_container_t> meshed_particles_for_gone_diegos;
   static void setup_draw_for_particles_for_new_gone_diegos(const std::vector<glm::vec3> &positions);
   // static void setup_draw_for_particles_for_gone_diegos(); // unused atm

   static meshed_particle_container_t meshed_particles_for_gone_diff_map_peaks;
   static void setup_draw_for_particles_for_gone_diff_map_peaks(const std::vector<std::pair<glm::vec3, float> > &positions);

   static bool draw_bad_nbc_atom_pair_markers_flag; // user can turn them off
   static void setup_draw_for_bad_nbc_atom_pair_markers();
   static void draw_bad_nbc_atom_pair_markers(stereo_eye_t eye, unsigned int pass_type);
   static void update_bad_nbc_atom_pair_marker_positions();
   static Texture texture_for_bad_nbc_atom_pair_markers;
   static TextureMesh tmesh_for_bad_nbc_atom_pair_markers;
   static std::vector<glm::vec3> bad_nbc_atom_pair_marker_positions;
   const unsigned int draw_count_max_for_bad_nbc_atom_pair_markers = 100; // needed?

   static void setup_draw_for_bad_nbc_atom_pair_dashed_line();
   static void update_bad_nbc_atom_pair_dashed_lines();
   static Mesh bad_nbc_atom_pair_dashed_line; // instanced mesh
   static void draw_bad_nbc_atom_pair_dashed_lines(stereo_eye_t eye, unsigned int pass_type);

   void setup_draw_for_chiral_volume_outlier_markers();
   static void draw_chiral_volume_outlier_markers(stereo_eye_t eye, unsigned int pass_type);
   static void update_chiral_volume_outlier_marker_positions();
   static Texture texture_for_chiral_volume_outlier_markers;
   static TextureMesh tmesh_for_chiral_volume_outlier_markers;

   static void update_hydrogen_bond_positions(); // if the intermediate atoms had hydrogen bond restraints, we can have dynamic
                                                 // hydrogen bonds. c.f. update_bad_nbc_atom_pair_marker_positions()

   static bool find_hydrogen_torsions_flag;

   static Texture texture_for_hud_colour_bar;
   static HUDTextureMesh tmesh_for_hud_colour_bar;

   // pickable moving atoms molecule
   // (we want to be able to avoid picking hydrogen atoms if the
   // are not displayed)
   pick_info moving_atoms_atom_pick(short int pick_mode) const;
   static short int in_moving_atoms_drag_atom_mode_flag;
   // when shift is pressed move (more or less) just the "local"
   // moving atoms atoms, we do this by making the shift proportional
   // to the square of the distance ratio.
   void move_moving_atoms_by_shear(int screenx, int screeny, short int squared_flag);
   void move_moving_atoms_by_shear_internal(const coot::Cartesian &diff_std,
					    short int squared_flag);
   void move_moving_atoms_by_simple_translation(int screenx, int screeny); // for rot/trans
   void move_single_atom_of_moving_atoms(int screenx, int screeny);
   // if control_is_pressed is true, then we only want to move the dragged atom
   void move_atom_pull_target_position(double screenx, double screeny, bool control_is_pressed);
   void move_dragged_anchored_atom(double screen_x, double screen_y);
   void add_target_position_restraint_for_intermediate_atom(const coot::atom_spec_t &spec,
							    const clipper::Coord_orth &target_pos);
   void add_target_position_restraints_for_intermediate_atoms(const std::vector<std::pair<coot::atom_spec_t, clipper::Coord_orth> > &atom_spec_position_vec); // refines after added
   short int rotate_intermediate_atoms_maybe(unsigned int widget_height, unsigned int widget_width);
                                                 // do it if have intermediate atoms
                                                 // and ctrl is pressed.
   // axis: 0 for Z, 1 for X.
   void rotate_intermediate_atoms_round_screen_z(double angle);
   void rotate_intermediate_atoms_round_screen_x(double angle);

   static void drag_intermediate_atom(const coot::atom_spec_t &atom_spec, const clipper::Coord_orth &pt);
   static void mark_atom_as_fixed(int imol, const coot::atom_spec_t &atom_spec, bool state);
   // this used to be called mark_atom_as_fixed() - I don't know why one would want to use it.
   static void while_moving_atoms_active_mark_atom_as_fixed(int imol, const coot::atom_spec_t &atom_spec, bool state);
   // static std::vector<mmdb::Atom *> fixed_intermediate_atoms;
   static bool fixed_atom_for_refinement_p(mmdb::Atom *); // examines the imol_moving_atoms molecule
                                                          // for correspondence

   // mol is new (not from molecules[imol]) molecule for the moving atoms.
   //
   atom_selection_container_t make_moving_atoms_asc(mmdb::Manager *mol,
						    int resno_1,
						    int resno_2) const;
   atom_selection_container_t make_moving_atoms_asc(mmdb::Manager *mol,
						    const std::vector<mmdb::Residue *> &residues) const;
   static bool moving_atoms_displayed_p() {
      if (moving_atoms_asc)
         if (moving_atoms_asc->mol)
            return true;
      return false;
   }
   // so that we know that fixed_points_sheared_drag_1 and
   // fixed_points_sheared_drag_2 are sensible:
   //
   static short int have_fixed_points_sheared_drag_flag;
   void set_fixed_points_for_sheared_drag();
   void do_post_drag_refinement_maybe();

   int last_restraints_size() const;

   static int dragged_refinement_steps_per_frame;
   static short int dragged_refinement_refine_per_frame_flag;
   static bool refinement_move_atoms_with_zero_occupancy_flag;

   //
   static bool draw_zero_occ_spots_flag;
   static bool draw_cis_peptide_markups;

   // static
   static std::string ccp4_defs_file_name();
   static int ccp4_projects_index_last;
   void set_directory_for_fileselection_string(std::string filename);
   void set_directory_for_saving_for_fileselection_string(std::string filename);

   void set_directory_for_filechooser_string(std::string filename);
   void set_directory_for_saving_for_filechooser_string(std::string filename);

   // 20220319-PE what do these do these days?
   static int file_chooser_dialog_x_size;
   static int file_chooser_dialog_y_size;

   // Origin marker for Johan
   static int show_origin_marker_flag;


   // Check waters interface for Florence:
   static int   check_waters_molecule;
   static float check_waters_b_factor_limit;
   static float check_waters_map_sigma_limit;
   static float check_waters_min_dist_limit;
   static float check_waters_max_dist_limit;
   static float check_waters_by_difference_map_sigma_level;
   static int   check_waters_by_difference_map_map_number;
   // save the dialog so that . and , can be used on it
   static GtkWidget *checked_waters_baddies_dialog;

   // zoom widget:
   //
   static void zoom_adj_changed(GtkAdjustment *adj, GtkWidget *window);
   static void set_zoom_adjustment(GtkWidget *widget);

   //
   static int show_paths_in_display_manager_flag;

   // scrollin' scrollin' scrollin'... Shall we stop? When shall we stop?
   static short int stop_scroll_diff_map_flag;
   static short int stop_scroll_iso_map_flag;
   static float stop_scroll_iso_map_level;
   static float stop_scroll_diff_map_level;

   // globbing [we want a "Filter" button so we need a list of
   // extensions that are coordinates].
   //
   // Same for "data".
   //
   static std::vector<std::string> *coordinates_glob_extensions;
   static std::vector<std::string> *data_glob_extensions;
   static std::vector<std::string> *map_glob_extensions;
   static std::vector<std::string> *dictionary_glob_extensions;
   void add_coordinates_glob_extension(const std::string &extension);
   void add_data_glob_extension(const std::string &extension);
   void add_map_glob_extension(const std::string &extension);
   void add_dictionary_glob_extension(const std::string &extension);
   void remove_coordinates_glob_extension(const std::string &extension);
   void remove_data_glob_extension(const std::string &extension);
   void remove_map_glob_extension(const std::string &extension);
   void remove_dictionary_glob_extension(const std::string &extension);
   static int filter_fileselection_filenames_flag;

   // superposition
   static int superpose_imol1;
   static int superpose_imol2;
   static std::string superpose_imol1_chain;
   static std::string superpose_imol2_chain;

   /*
   static void superpose_optionmenu_activate_mol1(GtkWidget *item, GtkPositionType pos);
   static void superpose_optionmenu_activate_mol2(GtkWidget *item, GtkPositionType pos);
   static void superpose_moving_chain_option_menu_item_activate (GtkWidget *item,
								 GtkPositionType pos);
   static void superpose_reference_chain_option_menu_item_activate (GtkWidget *item,
								    GtkPositionType pos);
   static void fill_superpose_option_menu_with_chain_options(GtkWidget *chain_optionmenu,
							     int is_reference_structure_flag);
   */

   static void superpose_reference_chain_combobox_changed(GtkWidget *combobox, gpointer data);
   static void superpose_moving_chain_combobox_changed(GtkWidget *combobox, gpointer data);

   static void superpose_combobox_changed_mol1(GtkWidget *c, gpointer data);
   static void superpose_combobox_changed_mol2(GtkWidget *c, gpointer data);
   static void fill_superpose_combobox_with_chain_options(int imol_active,
							  bool is_reference_structure_flag);
   static int         ramachandran_plot_differences_imol1;
   static int         ramachandran_plot_differences_imol2;
   static std::string ramachandran_plot_differences_imol1_chain;
   static std::string ramachandran_plot_differences_imol2_chain;
   static float rama_level_prefered;
   static float rama_level_allowed;
   static float rama_plot_background_block_size; // divisible into 360 preferably.
#ifdef DO_RAMA_PLOT
   static int rama_psi_axis_mode;
#endif

   /*
     Return the index of the superposed molecule - which could either be a
     new molecule (if move_imol2_flag was 1) or the imol2 or -1 (signifying
     failure to do the SMM superposition).
   */
   int superpose_with_atom_selection(atom_selection_container_t asc_ref,
				      atom_selection_container_t asc_mov,
				      int imol_mov,
				      std::string moving_mol_name,
				      std::string referennce_mol_name,
				      bool move_copy_of_imol2_flag);

#ifdef HAVE_SSMLIB
   void print_ssm_sequence_alignment(ssm::Align *SSMAlign,
				     atom_selection_container_t asc_ref,
				     atom_selection_container_t asc_mov,
				     mmdb::PAtom *atom_selection1,
				     mmdb::PAtom *atom_selection2,
				     int n_selected_atoms_1, int n_selected_atoms_2,
				     bool move_copy_of_imol2_flag);

   void make_and_print_horizontal_ssm_sequence_alignment(ssm::Align *SSMAlign,
							 atom_selection_container_t asc_ref,
							 atom_selection_container_t asc_mov,
							 mmdb::PAtom *atom_selection1,
							 mmdb::PAtom *atom_selection2,
							 int n_selected_atoms_1, int n_selected_atoms_2) const;

   void map_secondary_structure_headers(ssm::Align *SSMAlign,
					atom_selection_container_t asc_ref,
					atom_selection_container_t asc_mov,
					mmdb::PAtom *atom_selection1,
					mmdb::PAtom *atom_selection2,
					int n_selected_atoms_1, int n_selected_atoms_2) const;
   //
   void print_horizontal_ssm_sequence_alignment(std::pair<std::string, std::string> aligned_sequences) const;

   std::pair<std::string, std::string>
      get_horizontal_ssm_sequence_alignment(ssm::Align *SSMAlign,
					   atom_selection_container_t asc_ref,
					   atom_selection_container_t asc_mov,
					   mmdb::PAtom *atom_selection1, mmdb::PAtom *atom_selection2,
					   int n_selected_atoms_1, int n_selected_atoms_2) const;

#endif  // HAVE_SSMLIB


   // widget stuff:

   static std::pair<short int, float> float_from_entry(GtkWidget *entry);
   static std::pair<short int, int>   int_from_entry(GtkWidget *entry);

   // for widgets that are created by graphics_info_t functions (we
   // can't call c-interface.cc versions of these functions from here);
   static void      store_window_position(int window_type, GtkWidget *widget);
   static void set_transient_and_position(int widget_type, GtkWidget *window);
   static void set_transient_for_main_window(GtkWidget *dialog);

   // contour level saving
   void set_last_map_sigma_step(float level);

   // ----- Refmac state params: ----
   //
   // Normally a new map gets its own colour depending on the molecule
   // number, but often when you run refmac you want the new map to
   // have the colour of the old map (blue say) and the old map gets
   // coloured with the new colour scheme (horrible red, often).
   //
   static short int swap_pre_post_refmac_map_colours_flag;

   // ------ (More) water checking: -------
   void check_waters_by_difference_map(int imol_waters, int imol_diff_map,
				       int interactive_flag);

   // and a gui wrapper for that:
   GtkWidget *wrapped_create_checked_waters_by_variance_dialog(const std::vector <coot::atom_spec_t> &v, int imol);

   // and its buttons callbacks:
   // static void on_generic_atom_spec_button_clicked (GtkButton *button,
   //                                                  gpointer user_data);
   //
   // 20230520-PE the above should be toggle buttons (shouldn't they?) - they have the "toggled" signal,
   // not the "clicked" signal.
   //
   static void on_generic_atom_spec_toggle_button_toggled(GtkToggleButton *toggle_button, gpointer user_data);

   // ----- chiral volumes: ----

   void check_chiral_volumes(int imol);
   GtkWidget *wrapped_check_chiral_volumes_dialog(const std::vector <coot::atom_spec_t> &v,
						  int imol);
   static void on_inverted_chiral_volume_button_clicked(GtkButton *button,
							gpointer user_data);
   static int check_chiral_volume_molecule;
   // Tell us which residue types for chiral volumes restraints were missing:
   GtkWidget *wrapped_create_chiral_restraints_problem_dialog(const std::vector<std::string> &sv) const;
   static void check_chiral_volume_molecule_combobox_changed(GtkWidget *w, gpointer data);


   // unbonded star size - actually too messy to fix properly - so not used.
   static float unbonded_atom_star_size;

   // -------- povray/raster3d --------
   short int raster3d(std::string filename);
   short int povray(std::string filename);
   static float raster3d_bond_thickness;
   static float raster3d_atom_radius;
   static float raster3d_density_thickness;
   static int renderer_show_atoms_flag;
   static float raster3d_bone_thickness;
   static bool  raster3d_enable_shadows;
   static int raster3d_water_sphere_flag;
   static std::string raster3d_font_size;

   short int renderman(std::string filename);


   // ---- simple torsion ------
   // return true if all atoms found
   bool set_angle_tors(int imol,
		       const coot::atom_spec_t &as1,
		       const coot::atom_spec_t &as2,
		       const coot::atom_spec_t &as3,
		       const coot::atom_spec_t &as4);

   // ------- auto read mtz file: --------
   static int auto_read_do_difference_map_too_flag;

   // ------- refmac molecules option menu  -----
   static int refmac_molecule;

   // ------- new style refmac usage ---
   // above needs to be cleaned or fixed.
   static std::string mtz_file_for_refmac;

   // ------ new style combobox usage -------

   // the top one of this is probably what you want.
   std::string get_active_label_in_comboboxtext(GtkComboBoxText *combobox);
   std::string get_active_label_in_combobox(GtkComboBox *combobox) const;

   // ------ add OXT -------
   void fill_add_OXT_dialog_internal(GtkWidget *w);
   static int add_OXT_molecule;
   void fill_add_OXT_dialog_internal(GtkWidget *widget, int imol);
   static void add_OXT_molecule_combobox_changed(GtkWidget *widget, gpointer data);
   // return the default chain string (top of the list).
   // (return "no-chain" if it was not assigned (nothing in the list)).

   static std::string fill_option_menu_with_chain_options(GtkWidget *option_menu,
							  int imol,
							  GCallback signal_func);

   static std::string fill_combobox_with_chain_options(GtkWidget *combobox,
						       int imol,
						       GCallback f);
   static std::string fill_combobox_with_chain_options(GtkWidget *combobox,
						       int imol,
						       GCallback f,
						       const std::string &active_chain_id);

   // as above, except if one of the chain options is active_chain_id,
   // then set the active menu item to that.
   static std::string fill_option_menu_with_chain_options(GtkWidget *option_menu,
							  int imol,
							  GCallback signal_func,
							  const std::string &active_chain_id);
   // static void add_OXT_chain_menu_item_activate (GtkWidget *item,
   //GtkPositionType pos);
   static void add_OXT_chain_combobox_changed(GtkWidget *combobox, gpointer data);

   //
   static GtkWidget *wrapped_nothing_bad_dialog(const std::string &label, bool use_markup=false);

   static std::string add_OXT_chain;

   // ----- merge molecules ------
   static int merge_molecules_master_molecule;
   static std::vector<int> *merge_molecules_merging_molecules;
   static coot::residue_spec_t merge_molecules_ligand_spec; // JED feature
   void set_merge_molecules_ligand_spec(const coot::residue_spec_t &spec_in);

   // ------ change chain ids:
   static int change_chain_id_molecule;
   static std::string change_chain_id_from_chain;


   // ----- renumber residue range -------
   static void fill_renumber_residue_range_dialog(GtkWidget *w);
   static void renumber_residue_range_molecule_combobox_changed(GtkWidget *combobox,
								gpointer data);
   void fill_renumber_residue_range_internal(GtkWidget *w, int imol);
   static void renumber_residue_range_chain_menu_item_select(GtkWidget *item, GtkPositionType pos);
   static void renumber_residue_range_chain_combobox_changed(GtkWidget *combobox, gpointer data);
   static int renumber_residue_range_molecule;
   static std::string renumber_residue_range_chain;

   // -------- public toggle button control: ---------
   void untoggle_model_fit_refine_buttons_except(const std::string &button_name);
   static void set_model_fit_refine_button_names(GtkWidget *w);
   static void set_other_modelling_tools_button_names(GtkWidget *w);

   static std::vector<keyboard_key_t> keyboard_key_history;
   static void add_key_to_history(const keyboard_key_t &kk) { keyboard_key_history.push_back(kk); };
   static void check_keyboard_history_for_easter_egg_codes();

   // -------- keyboard rotamer control: ---------
   static void rotamer_dialog_next_rotamer();
   static void rotamer_dialog_previous_rotamer();
   static void rotamer_dialog_neighbour_rotamer(int istep); // could be private

   // -------- keyboard difference map peak control: ------------
   static void difference_map_peaks_next_peak();
   static void difference_map_peaks_previous_peak();
   static void difference_map_peaks_neighbour_peak(int istep); // could be private

   //  -------- keyboard baddie waters control: ------
   static void checked_waters_next_baddie(int dir);

   // ------ NCS ------
   static float ncs_min_hit_ratio;
   static short int ncs_maps_do_average_flag;
   /* At what level of homology should we say that we can't see homology
      for NCS calculation? (default 0.8) */
   static float ncs_homology_level;
   static short int ncs_matrix_flag;

   // 20230621-PE not yet implemented, although see l28 lambda (which doesn't work)
   static void ncs_jump_forward();  // A -> B -> C -> A
   static void ncs_jump_backward(); // B -> A

   // ------------- validation: -------------------

   // pretty graphs on a canvas: (not const because geom_p may be added to)
   void geometric_distortion(int imol);
   void b_factor_graphs(int imol);
   void calc_b_factor_graphs(int imol); //// NEW ADDITION BY RICHARD
   void omega_graphs(int imol);
   coot::rotamer_graphs_info_t rotamer_graphs(int imol); // give results back to scripting layer
   void density_fit_graphs(int imol);
   static void diff_map_peaks_dialog_update_button_clicked_func(GtkButton *button, gpointer user_data); // called by below
   static void fill_difference_map_peaks_button_box();

   static void show_diff_map_peaks_vbox(int imol_map, int imol_coords,
                                                          const std::vector<std::pair<clipper::Coord_orth, float> > &centres,
                                                          float n_sigma,
                                                          bool do_positive_level_flag,
                                                          bool do_negative_level_flag,
                                                          bool around_model_only_flag);
   static void hide_vertical_validation_frame_if_appropriate();
   // the buttons callback for above:
   static void on_diff_map_peak_button_selection_toggled (GtkToggleButton *button,
							  gpointer         user_data);

   static std::vector<clipper::Coord_orth> *diff_map_peaks;
   static int max_diff_map_peaks;
   void clear_diff_map_peaks();
   static float difference_map_peaks_sigma_level;

   // ---------------- backup filenames ----------------------
   static bool unpathed_backup_file_names_flag;
   static bool decoloned_backup_file_names_flag;
   static int backup_compress_files_flag;

   // --------- Miguel's axis orientation matrix ---------------
   static GL_matrix axes_orientation;
   static short int use_axes_orientation_matrix_flag;

   // --------- preferences ---------------
   static GtkWidget *preferences_widget;
   static int mark_cis_peptides_as_bad_flag;

   static std::vector<std::string> preferences_general_tabs;
   static std::vector<std::string> preferences_bond_tabs;
   static std::vector<std::string> preferences_geometry_tabs;
   static std::vector<std::string> preferences_colour_tabs;
   static std::vector<std::string> preferences_map_tabs;
   static std::vector<std::string> preferences_other_tabs;

   static std::vector<coot::preferences_icon_info_t> *model_toolbar_icons;
   static std::vector<coot::preferences_icon_info_t> *main_toolbar_icons;

   short int save_preference_file(const std::string &filename, short int il);
   static std::vector<coot::preference_info_t> preferences_internal;
   static std::vector<coot::preference_info_t> preferences_internal_default;
   void make_preferences_internal();
   void preferences_internal_change_value(int preference_type, int ivalue);
   void preferences_internal_change_value(int preference_type, float fvalue);
   void preferences_internal_change_value(int preference_type,
					  float fvalue1, float fvalue2, float fvalue3);
   void preferences_internal_change_value(int preference_type, int ivalue1, int ivalue);

   static void preferences_model_toolbar_icon_toggled(GtkCellRendererToggle *button,
					      gchar *path,
		    			      gpointer data);
   static void update_toolbar_icons(GtkTreeModel *model, int toolbar_index);
   static void update_main_toolbar_icons(GtkTreeModel *model);
   static void update_model_toolbar_icons(GtkTreeModel *model);
   void fill_preferences_model_toolbar_icons(GtkWidget *preferences,
					     GtkWidget *scrolled_window);
   static void preferences_main_toolbar_icon_toggled(GtkCellRendererToggle *button,
					      gchar *path,
		    			      gpointer data);
   static void preferences_toolbar_icon_toggled(GtkCellRendererToggle *button,
					      gchar *path,
		    			      gpointer data,
                                              int toolbar_index);
   void fill_preferences_main_toolbar_icons(GtkWidget *preferences,
					     GtkWidget *scrolled_window);
   void fill_preferences_toolbar_icons(GtkWidget *preferences,
				       GtkWidget *scrolled_window,
				       int toolbar_index);

   void show_hide_toolbar_icon_pos(int pos, int show_hide_flag, int toolbar_index);
   std::vector<int> get_model_toolbar_icons_list();
   std::vector<int> get_main_toolbar_icons_list();
   void add_to_preferences(const std::string &file_name, const std::string &contents) const;
   std::string get_preferences_directory() const;

   // --- remote controlled coot: ----
   static int try_port_listener;
   static int remote_control_port_number;
   static std::string remote_control_hostname;
   static int coot_socket_listener_idle_function_token; // -1 default (off)
   static std::string socket_string_waiting;
   static std::string socket_python_string_waiting;
   static volatile bool have_socket_string_waiting_flag;
   static volatile bool have_socket_python_string_waiting_flag;
   static volatile bool socket_string_waiting_mutex_lock;
#ifdef USE_GUILE
   static SCM safe_scheme_command(const std::string &scheme_command);
   static SCM process_socket_string_waiting();
#endif
   static gboolean process_socket_string_waiting_bool(gpointer user_data);
   static gboolean process_socket_python_string_waiting_bool(gpointer user_data);

   // --------- Tip of the Day ---------------
   static short int do_tip_of_the_day_flag;

   // --------- LSQing ---------------
   //
   // if cell and space group are both not null then change the cell
   // and space group of the moving molecule to that of the reference
   // (presumably).
   std::pair<int, clipper::RTop_orth> apply_lsq(int imol_ref, int imol_moving,
						const std::vector<coot::lsq_range_match_info_t> &matches);
   // sometimes we just want the matrix but not to actually move the coordinates.
   std::pair<int, clipper::RTop_orth> lsq_get_and_apply_matrix_maybe(int imol_ref, int imol_moving,
								     const std::vector<coot::lsq_range_match_info_t> &matches,
								     bool apply_matrix);
   static std::vector<coot::lsq_range_match_info_t> *lsq_matchers;
   // the simple widget for LSQing: (perhaps these should be vectors of strings
   // in the general case (more complex widget)?)
   static std::string lsq_match_chain_id_ref;
   static std::string lsq_match_chain_id_mov;
   static int lsq_ref_imol;
   static int lsq_mov_imol;
   static lsq_dialog_values_t lsq_dialog_values;

   // -------- some people don't want restype and slashes in atom label:
   static short int brief_atom_labels_flag;
   //          and some people want seg-ids in their atom labels (Francesca)
   static short int seg_ids_in_atom_labels_flag;

   GtkWidget *wrapped_create_display_control_window();

   void delete_molecule_from_display_manager(int imol, bool was_map_flag);

   // -------- undo move: suggested by Frank von Delft -----
   void undo_last_move();

   // ---- auto read MTZ file ---

   // kill these off when the new system works.
/*    static std::string auto_read_MTZ_FWT_col; */
/*    static std::string auto_read_MTZ_PHWT_col; */
/*    static std::string auto_read_MTZ_DELFWT_col; */
/*    static std::string auto_read_MTZ_PHDELWT_col; */

   static std::vector<coot::mtz_column_trials_info_t> user_defined_auto_mtz_pairs;

   // ---- cis trans conversion ---
   void cis_trans_conversion(mmdb::Atom *at, int imol, short int is_N_flag);

   // return true if the isomerisation was made
   //
   bool cis_trans_conversion_intermediate_atoms();

   // symmetry control dialog:
   GtkWidget *wrapped_create_symmetry_controller_dialog() const;
   static GtkWidget *symmetry_controller_dialog;  // returned by above

   // ----- LSQ Plane -----
   static std::vector<clipper::Coord_orth> *lsq_plane_atom_positions;
   int add_lsq_plane_atom(int imol, int atom_index);
   int measure_lsq_plane_deviant_atom(int imol, int atom_index);
   int remove_last_lsq_plane_atom();
   void render_lsq_plane_atoms(); // put a blob at atoms in lsq_plane_atom_positions
   static short int in_lsq_plane_deviation;
   static short int in_lsq_plane_define;

   GtkWidget *wrapped_create_lsq_plane_dialog();
   static GtkWidget *lsq_plane_dialog;


   // -------- Fffearing -------------
   static float fffear_angular_resolution;

   // -------- Base Pairing (Watson Crick) -------------
   static int in_base_paring_define;

   // -------- Meshes control (i.e. the Meshes of molecule_class_info)
   void set_show_molecular_representation(int imol, unsigned int mesh_idx, bool on_off);

   void update_molecular_representation_widgets();
   static void molecular_representation_meshes_checkbutton_toggled(GtkCheckButton *button, gpointer *user_data);
   static void undisplay_all_molecule_meshes(int imol); // for imol, of course

   int add_molecular_representation(int imol,
                                    const std::string &atom_selection,
				    const std::string &colour_scheme,
				    const std::string &style,
                                    int secondary_structure_usage_flag);
   int add_ribbon_representation_with_user_defined_colours(int imol, const std::string &name);
   void remove_molecular_representation(int imol, int idx);

   // -------- Texture Meshes (for importing glTF models) -------------
   static std::vector<TextureMesh> texture_meshes;
   static void draw_texture_meshes(stereo_eye_t eye);

   static Mesh mesh_for_eyelashes;
   static Mesh &get_mesh_for_eyelashes();

   // these are for non-molecule based generic display objects using instancing
   static std::vector<Instanced_Markup_Mesh> instanced_meshes;

   static meshed_generic_display_object mesh_for_environment_distances;
   static meshed_generic_display_object mesh_for_pointer_distances;

   static GtkWidget *generic_objects_dialog;
   static std::vector<meshed_generic_display_object> generic_display_objects;
   static bool is_valid_generic_display_object_number(int obj_no);
   static void from_generic_object_remove_last_item(int object_number);

   static void set_display_generic_object_simple(int object_number, short int istate) {
      // std::cout << "in set_display_generic_object_simple " << object_number << " " << istate << std::endl;
      graphics_info_t g;
      if (object_number >= 0) {
         if (object_number < int(g.generic_display_objects.size())) {
            g.generic_display_objects[object_number].mesh.set_draw_this_mesh(istate);
         } else {
            std::cout << "ERROR:: BAD object_number in set_display_generic_object_simple: "
                      << object_number << std::endl;
         }
      } else {
         std::cout << "ERROR:: BAD object_number in set_display_generic_object_simple: "
                   << object_number << std::endl;
      }
   }

   static void
   on_generic_objects_dialog_object_check_button_toggled(GtkButton       *button,
                                                         gpointer         user_data) ;
   void generic_objects_dialog_grid_add_object_internal(const meshed_generic_display_object &gdo,
                                                        GtkWidget *dialog,
                                                        GtkWidget *grid,
                                                        int io);
   int new_generic_object_number(const std::string &name);
   int new_generic_object_number_for_molecule(const std::string &name, int imol);
   static int generic_object_index(const std::string &name);

   static void myglLineWidth(int n_pixels);

   // ---- active atom:
   static std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom_spec();
   static std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom_spec(int imol);
   static std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom_spec_internal(int imol);
   static std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom_spec_simple();
   // direct for immediate usage, imol and atom,
   // return (-1, null) on not found
   std::pair<int, mmdb::Atom *> get_active_atom() const;

   // this can return -1 if there is no active atom molecule.
   int copy_active_atom_molecule();

   // ---- open url
   // Hmm.. do we need a vector here?
   static std::string browser_open_command;

   // -- variable bond width (lines get thinner as we zoom out)
   static bool use_variable_bond_width;

   // -- default bond width
   static int default_bond_width;
   static int default_bonds_box_type; // Phil want to configure this.
   static bool draw_stick_mode_atoms_default; // true,

   // ---- default sigma level:
   static float default_sigma_level_for_map;
   static float default_sigma_level_for_fofc_map;

   // probe dots on intermediate atoms (we need to have hydrogens)
   static short int do_probe_dots_on_rotamers_and_chis_flag;
   static short int do_probe_dots_post_refine_flag;
   static bool do_coot_probe_dots_during_refine_flag;
   void do_probe_dots_on_rotamers_and_chis();
   void do_probe_dots_post_refine();
   void do_interactive_probe() const; // molprobity probe
   // not const because it manipulates generic graphics objects
   void do_interactive_coot_probe(); // coot probe

   // can be private?
   void setup_for_probe_dots_on_chis_molprobity(int imol);
   static coot::Cartesian probe_dots_on_chis_molprobity_centre;
   static float probe_dots_on_chis_molprobity_radius;

   // a text string and a handle (so that it can be removed)
   static std::vector<coot::generic_text_object_t> generic_texts;
   // not that draw_generic_texts() is part of draw_molecule_atom_labels().

   // -- move molecule here
   static int move_molecule_here_molecule_number;
   // static void move_molecule_here_item_select(GtkWidget *item,
   // GtkPositionType pos);
   static void move_molecule_here_combobox_changed(GtkWidget *combobox, gpointer data);


   // -- make the key user changable.  A template for other bound
   // functions:
   static int ncs_next_chain_skip_key;
   static int ncs_prev_chain_skip_key;
   static int update_go_to_atom_from_current_residue_key;

   // --- keyboarding the Go To Atom (or residue, really)
   static void apply_go_to_residue_keyboading_string(const std::string &text);
   void apply_go_to_residue_keyboading_string_inner(int imol, mmdb::Atom *new_centre_atom);
   // go to the middle residue of the first occurance of the sequence triplet if you can
   // seq_trip is of course something like "ACE"
   // return the "found the triplet and moved there" status: 0 for fail.
   //
   int apply_go_to_residue_from_sequence_triplet(int imol, const std::string &seq_trip);

   // -- PHENIX support
   static std::string external_refinement_program_button_label;

   // ---- pseudo bond for sec str restraints
   static coot::pseudo_restraint_bond_type pseudo_bonds_type;


   // ---- a debugging thing for user to give me feedback
   static int debug_atom_picking;

   // ---- default new atoms b factor
   static float default_new_atoms_b_factor;

   // reset the b factors after atoms were moved, default no
   static int reset_b_factor_moved_atoms_flag;

   // ---- disable state script writing on exiting (no-guano)
   static bool disable_state_script_writing;

   // ---- man this is tricky.  Store the listener socket state, so
   // that we can call coot-listener-idle-function-proc if there is a
   // good socket.  The actual socket is a scheme variable.
   //
   static int listener_socket_have_good_socket_state;

   void post_recentre_update_and_redraw();

   // --- map sharpening
   static int imol_map_sharpening;
   static float map_sharpening_scale_limit;

   // -- remarks browswer
   static int imol_remarks_browswer;

   // --- user defined picks
   static std::vector<coot::atom_spec_t> user_defined_atom_pick_specs;
   static bool residue_type_selection_was_user_picked_residue_range;

   // --- electrostatic charnge range scale (typically 0.5)
   static float electrostatic_surface_charge_range;

   // --- nudge active residue
   static void nudge_active_residue(guint direction);
   static void nudge_active_residue_by_rotate(guint direction);

   // --- curl handlers
#ifdef USE_LIBCURL
   static std::vector<coot::simple_curl_handler_t> curl_handlers;
   // we have inner function for the add and remove handles so that we can use a lock.
   static volatile bool curl_handlers_lock;
   bool add_curl_handle_and_file_name(std::pair<CURL *, std::string> p);
   bool add_curl_handle_and_file_name_inner(std::pair<CURL *, std::string> p);
   bool remove_curl_handle_with_file_name(std::string file_name);
   bool remove_curl_handle_with_file_name_inner(std::string file_name);
   // return NULL on no such filename being transfered.
   CURL *get_curl_handle_for_file_name(const std::string &filename) const;
   CURL *get_curl_handle_for_file_name_inner(const std::string &filename) const;
   static bool curl_handler_stop_it_flag_set(CURL *c);
   static bool curl_handler_stop_it_flag_set_inner(CURL *c);
   static void set_stop_curl_download_flag(const std::string &file_name);
   static void set_stop_curl_download_flag_inner(const std::string &file_name);
#endif

#ifdef USE_GUILE
   static SCM user_defined_click_scm_func;
   SCM atom_spec_to_scm(const coot::atom_spec_t &spec) const;
#endif
#ifdef USE_PYTHON
   static PyObject *user_defined_click_py_func;
   PyObject *atom_spec_to_py(const coot::atom_spec_t &spec) const;
#endif
   void run_user_defined_click_func();

   // --- unapply symmetry to current view, (we are looking at
   // symmetry and we want to get back to the main molecule,
   // preserving the orientation, if possible.
   //
   // pre_translation is the translation that needed to be applied to
   // the molecule so that it was close to the origin (from there we
   // do the symmetry expansion, and it seems that st generates the
   // symmetry-related molecule that we are looking at now).
   int unapply_symmetry_to_view(int imol, const std::vector<std::pair<clipper::RTop_orth, clipper::Coord_orth> > &symm_mat_and_pre_shift);

   //
   int move_reference_chain_to_symm_chain_position();


#ifdef USE_MYSQL_DATABASE
   // MYSQL database
   static MYSQL *mysql;
   static int query_number;
   static std::string mysql_host;
   static std::string mysql_user;
   static std::string mysql_passwd;
string   static std::string sessionid;
   static std::pair<std::string, std::string> db_userid_username;
#endif

   static std::vector<coot::view_info_t> views;
   static float views_play_speed;

   static std::string movie_file_prefix;
   static int movie_frame_number;
   static int make_movie_flag;
   static void dump_a_movie_image(); // should be private
   static int screendump_image(const std::string &file_name);


   // ------------------------- restraints editor ----------------------------
   //
   static std::vector<coot::restraints_editor> restraints_editors;
   coot::restraints_editor get_restraints_editor(GtkWidget *w);
   void clear_restraints_editor_by_dialog(GtkWidget *dialog);

   // Kevin Keating (for example) wants to be able set torsion
   // restraints but not have those "fight" the built-in torsion
   // restraints. i.e. the torsion restraints should be the
   // user_defined ones only.  This is off (0) by default, but the use
   // can turn it on - and then the user-defined torsion restraints
   // will take effect (and not the built-in ones).  Maybe we'd want
   // to do this sort of things with bonds and angles too - but I
   // don't see it yet.  (That may require a rework).
   //
   static bool use_only_extra_torsion_restraints_for_torsions_flag;


   // We want --python to give us a python prompt with --no-graphics.
   // To do that, c_inner_main looks at the "python at the prompt"
   // flag.  This flag gets set in main(), from the command line
   // argument processing.
   static short int python_at_prompt_flag;

   // all molecule rotamer score, (depends on private rotamer probability tables)
   coot::rotamer_score_t all_molecule_rotamer_score(int imol) const;

   // update self?
   static bool update_self;

   // what shall we do with nomenclature errors on reading pdb files?
   //
   static coot::nomenclature_error_handle_type nomenclature_errors_mode;

   void multi_torsion_residues(int imol, const std::vector<coot::residue_spec_t> &v);
   static void on_multi_residue_torsion_button_clicked(GtkButton *button, gpointer user_data);
   void rotate_multi_residue_torsion(double x, double y);

   // Old style OpenGL:
   // bottom left ligand view
   // void setup_graphics_ligand_view_aa();
   // void setup_graphics_ligand_view_aa(int imol); // only allow imol to be potential active residue.
   // void setup_graphics_ligand_view(int imol, mmdb::Residue *residue, const std::string &alt_conf);
   // which stores in:
   // static graphics_ligand_molecule graphics_ligand_mol;
   //
   // void close_graphics_ligand_view_for_mol(int imol_in) {
   // if (graphics_ligand_mol.imol == imol_in)
   // graphics_ligand_view_flag = false;
   // }
   // static void graphics_ligand_view();  // actually draw it

   // replaced by
   //

   // 20220723-PE There functions have (will have) internals that are blanked out by
   // EMSCRIPTEN conditional compilation
   void setup_graphics_ligand_view_using_active_atom(int imol);
   void setup_graphics_ligand_view_using_active_atom();

   static graphics_ligand_mesh_molecule_t graphics_ligand_mesh_molecule;
   void setup_draw_for_ligand_view();
   void close_graphics_ligand_view_for_mol(int imol_in) {
      if (graphics_ligand_mesh_molecule.imol == imol_in)
         graphics_ligand_view_flag = false;
   }
   void setup_graphics_ligand_view(int imol, mmdb::Residue *residue, const std::string &alt_conf);

   static int show_graphics_ligand_view_flag; // user control, default 1 (on).

   // don't redraw everything, just those that have a residue with name res_name
   //
   void redraw_molecules_with_residue(const std::string &res_name);

   // use-defined flev params
   static float fle_water_dist_max;   // 3.25
   static float fle_h_bond_dist_max;  // 3.9

   // e.g. user_name_passwd_map["proxy"] -> "fred", "bill3"
   static std::map<std::string, std::pair<std::string, std::string> > user_name_passwd_map;

   // Mogul (default is 5.0)
   static float mogul_max_badness;

   // place helix here fudge factor (for EM maps?)
   static float place_helix_here_fudge_factor;

   coot::geometry_distortion_info_container_t geometric_distortions(int imol, mmdb::Residue *residue_p,
								    bool with_nbcs);

   // tabulate_geometric_distortions runs geometric_distortions() on restraints.
   void tabulate_geometric_distortions(coot::restraints_container_t &restraints) const;

   static bool linked_residue_fit_and_refine_state;

   static bool allow_duplseqnum;

   static short int probe_available; // need no, yes, don't-know-yet

   void perpendicular_ligand_view(int imol, const coot::residue_spec_t &residue_spec);

   static double map_to_model_correlation_atom_radius;

   static std::vector<std::pair<clipper::Coord_orth, std::string> >
     user_defined_interesting_positions;
   static unsigned int user_defined_interesting_positions_idx;

   void register_user_defined_interesting_positions(const std::vector<std::pair<clipper::Coord_orth, std::string> > &udip);

   // atom pull restraint
   // static atom_pull_info_t atom_pull; 20180218 just one
   static std::vector<atom_pull_info_t> atom_pulls;
   static unsigned int n_atom_pulls; // look inside atom_pulls and count the number of pulls with atoms
   static void all_atom_pulls_off();
   static void atom_pull_off(const coot::atom_spec_t &spec);
   static void atom_pulls_off(const std::vector<coot::atom_spec_t> &specs);
   void add_or_replace_current(const atom_pull_info_t &atom_pull_in);
   static void draw_atom_pull_restraints(stereo_eye_t eye);
   // we don't want to refine_again if the accept/reject dialog "Accept" button was clicked
   // (not least because now the refined atoms have gone out of scope)
   void clear_atom_pull_restraint(const coot::atom_spec_t &spec, bool refine_again_flag);
   void clear_atom_pull_restraints(const std::vector<coot::atom_spec_t> &specs,
				   bool refine_again_flag) {
      for (std::size_t i=0; i<specs.size(); i++)
	 clear_atom_pull_restraint(specs[i], false);
      if (refine_again_flag)
	 if (last_restraints)
	    drag_refine_refine_intermediate_atoms();
   }

   void clear_all_atom_pull_restraints(bool refine_again_flag);
   static bool auto_clear_atom_pull_restraint_flag;
   static LinesMesh lines_mesh_for_pull_restraint_neighbour_displacement_max_radius_ring;
   static void draw_intermediate_atoms_pull_restraint_neighbour_displacement_max_radius_ring();
   static void setup_lines_mesh_for_proportional_editing(); // called from gl widget realize function

   static bool continue_update_refinement_atoms_flag;

   // for CFC, no graphics_draw()
   void display_all_model_molecules();
   void undisplay_all_model_molecules_except(int imol);
   void undisplay_all_model_molecules_except(const std::vector<int> &keep_these);
   static GtkWidget *cfc_dialog;

   // new CFC - can scrap the above when this works:
   static cfc_gui_t cfc_gui;

   static bool do_intermediate_atoms_rama_markup; // true
   static bool do_intermediate_atoms_rota_markup; // false

   static Instanced_Markup_Mesh rama_balls_mesh;
   void setup_rama_balls();
   void update_rama_balls(std::vector<Instanced_Markup_Mesh_attrib_t> *balls_p);

   static bool all_atom_contact_dots_ignore_water_flag; // false by default
   static bool all_atom_contact_dots_do_vdw_surface; // false by default
   static float contact_dots_density; // 1 by default
   static void setup_cylinder_clashes(const coot::atom_overlaps_dots_container_t &c,
                                      int imol, float tube_radius, bool extra_annotation=false);

   static unsigned int bond_smoothness_factor; // default 1, changes num_subdivisions and n_slices

   bool coot_all_atom_contact_dots_are_begin_displayed_for(int imol) const;
   void coot_all_atom_contact_dots_instanced(mmdb::Manager *mol, int imol); // creates/updates
   // meshes in molecules.
   static float contact_dot_sphere_subdivisions;
   static bool get_exta_annotation_state();

   static void fill_rotamer_probability_tables();
   static std::pair<bool, float> model_display_radius;

   void set_model_display_radius(bool on_off, float radius_in);
   // - maybe this is not the best place for this function?
   static bool is_within_display_radius(const coot::CartesianPair &p);
   static bool is_within_display_radius(const coot::Cartesian &p);

   static bool make_auto_h_bond_restraints_flag;

   static bool cif_dictionary_file_selector_create_molecule_flag;

   static double geman_mcclure_alpha;

   static void set_geman_mcclure_alpha(float alpha); // reruns refinement if we have restraints

   static bool update_maps_on_recentre_flag;

   static double lennard_jones_epsilon;

   static void set_lennard_jones_epsilon(float epsilon); // reruns refinement if we have restraints

   static double log_cosh_target_distance_scale_factor;

   static std::pair<bool,float> coords_centre_radius;  // should the display radius limit be applied? And
                                                       // if so, what is it? (say 20A)
                                                       // used in draw_bonds().

   static translation_gizmo_t translation_gizmo;
   static Mesh translation_gizmo_mesh;
   static void setup_draw_for_translation_gizmo();
   static void draw_translation_gizmo(stereo_eye_t eye);
   static bool translation_gizmo_is_being_dragged;
   static translation_gizmo_t::pick_info_t translation_gizmo_axis_dragged;
   static translation_gizmo_t::pick_info_t translation_gizmo_picked();
   static void translate_things_on_translation_gizmo_dragged(); // including the gizmo

   // extensions registry
   // a name (a script file name) and a version number/identifier as a string
   //
   static std::map<std::string, std::string> extensions_registry;
   // return empty string on extension-not-found
   void register_extension(const std::string &extension,
			   const std::string &version);
   std::string get_version_for_extension(const std::string &extension_name) const;

   void eigen_flip_active_residue(); // is this right?
   // Use this in future. Can return NULL.
   static mmdb::Atom *get_moving_atoms_active_atom(const coot::Cartesian &rc, float within_radius_limit);

   static int jed_flip_intermediate_atoms(bool invert_atom_selection);
   static int crankshaft_peptide_rotation_optimization_intermediate_atoms();
   static int side_chain_flip_180_intermediate_atoms();
   // this one is used by clicking on the fip-180 modelling tooglebutton and then picking on an active atom
   void side_chain_flip_180_moving_atoms_residue(const coot::residue_spec_t &spec,
                                                 const std::string &alt_conf);

   // ---------------------------------------------

   static float view_rotation_per_pixel_scale_factor;

   // ---------------------------------------------

   /* model-view-projection matrix */
   static float *mvp; // needed?
   static int mvp_location;            // GLSL
   static int view_rotation_location;  // GLSL
   // static glm::quat glm_quat;
   void set_view_quaternion(float i, float j, float k, float l);
   static glm::vec3 get_rotation_centre();
   static clipper::Coord_orth get_rotation_centre_co();
   static void add_to_rotation_centre(const glm::vec3 &offset);

   // static bool using_trackpad;
   static bool use_primary_mouse_for_view_rotation_flag;
   static double mouse_x;
   static double mouse_y;
   static double drag_begin_x; // gtk pixels
   static double drag_begin_y;
   static double mouse_speed; // default 1.0, but adjusted to be bigger for wide display
   static std::pair<double, double> mouse_previous_position;
   static void set_mouse_previous_position(double x, double y);
   static double get_mouse_previous_position_x();
   static double get_mouse_previous_position_y();
   static glm::quat view_quaternion;
   // static void update_view_quaternion(int area_width, int area_height);
   // let's copy the one from crows:
   static void update_view_quaternion(int glarea_width, int glarea_height,
                                      double delta_x_drag, double delta_y_drag);

   // static benny::Camera camera;
   // static Transform transform;

   static float trackball_project_to_sphere(float r, float x, float y);
   static glm::quat trackball_to_quaternion(float p1x, float p1y, float p2x, float p2y, float trackball_size);

   // Shader things
   static bool draw_the_other_things;
   static GLuint other_triangles_vertexarray_id;
   static GLuint other_triangles_array_buffer_id;
   static GLuint other_triangles_index_buffer_id;
   static GLuint programID_for_central_cube;
   static GLuint central_cube_vertexarray_id;
   static GLuint central_cube_array_buffer_id;
   static GLuint central_cube_index_buffer_id;
   static GLuint hud_text_vertexarray_id;
   static GLuint rotation_centre_crosshairs_vertexarray_id;
   static GLuint rotation_centre_crosshairs_vertex_buffer_id;
   static GLuint rotation_centre_crosshairs_index_buffer_id;
   // STATIC GLuint framebuffer_id;  // now we get access from the framebuffer class
   static GLuint blur_x_quad_vertex_array_id;
   static GLuint blur_y_quad_vertex_array_id;
   static GLuint combine_textures_using_depth_quad_vertex_array_id;
   static GLuint screen_quad_vertex_array_id;
   static GLuint blur_quad_vertex_array_id; // old
   static GLuint textureColorbuffer_screen;
   static GLuint textureColorbuffer_blur;
   static GLuint hud_text_array_buffer_id;
   static Shader shader_for_outline_of_active_residue;
   static Shader shader_for_maps;
   static Shader shader_for_map_caps;
   static Shader shader_for_models;
   static Shader shader_for_model_as_meshes; // _as_Model in due course
   static Shader shader_for_moleculestotriangles;
   static Shader shader_for_moleculestotriangles_with_shadows;
   static Shader shader_for_origin_cube;
   static Shader shader_for_central_cube;
   static Shader shader_for_rotation_centre_cross_hairs_for_ssao; // central_cube by a modern name
   static Shader shader_for_hud_text;
   static Shader shader_for_hud_buttons;
   static Shader shader_for_hud_image_texture;
   static Shader shader_for_atom_labels;
   static Shader shader_for_x_blur;
   static Shader shader_for_y_blur;
   static Shader shader_for_dof_blur_by_texture_combination;
   static Shader shader_for_symmetry_atoms_bond_lines;
   static Shader shader_for_lines;
   static Shader shader_for_rama_balls;
   static Shader shader_for_particles;
   static Shader shader_for_instanced_objects;
   static Shader shader_for_extra_distance_restraints;
   static Shader shader_for_hud_geometry_bars;
   static Shader shader_for_hud_geometry_labels; // for labels image
   static Shader shader_for_lines_pulse; // "you are here" pulse
   static Shader shader_for_ligand_view;
   static Shader shader_for_happy_face_residue_markers;
   static Shader shader_for_happy_face_residue_markers_for_ssao;
   // static Shader shader_for_anchored_atom_markers; // this now uses the happy face residue marker shader
   static Shader shader_for_rama_plot_axes_and_ticks;
   static Shader shader_for_rama_plot_phi_phis_markers;
   static Shader shader_for_hud_lines; // actally in 3D because it uses LinesMesh class
   static Shader shader_for_texture_meshes;
   static long frame_counter;
   static float fps; // for on-screen FPS (fps is not calculated every frame)
   static float fps_std_dev; // for on-screen FPS IQR (fps is not calculated every frame)
   static long frame_counter_at_last_display;
   static bool perspective_projection_flag;
   static float perspective_fov;
   static float screen_z_near_perspective;
   static float screen_z_far_perspective;
   static float goodselliness;
   static bool shader_do_ambient_occlusion_flag;
   static bool shader_do_depth_blur_flag;
   static bool shader_do_depth_fog_flag;
   static bool shader_do_outline_flag;
   static bool shader_do_depth_of_field_blur_flag;
   static bool draw_normals_flag;
   static bool use_framebuffers;
   static framebuffer screen_framebuffer;
   static framebuffer blur_x_framebuffer;
   static framebuffer blur_y_framebuffer;
   static framebuffer blur_framebuffer; // from 2020
   static framebuffer combine_textures_using_depth_framebuffer;
   static unsigned int framebuffer_scale;

   void set_perspective_fov(float angle) { perspective_fov = angle; } // in degress (typically 30 or so)

   // ---------------------------------------------
   static bool shaders_have_been_compiled;
   bool init_shaders(); // return status (true = OK)
   bool init_shader(const std::string &shader_file_name);
   void init_framebuffers(unsigned int width, unsigned int height);// 20220129-PE a crows thing

   static int scale_up_graphics;    // default 1
   static int scale_down_graphics;  // default 1

   // draw-2 functions
   void init_screen_quads();
   void init_blur_quads();
   void init_central_cube();
   void init_buffers();
   void init_hud_text();
   static void draw_hud_geometry_bars();
   static void draw_hud_geometry_tooltip();
   static bool draw_hud_tooltip_flag;
   static glm::mat4 get_molecule_mvp(stereo_eye_t eye, bool debug_matrices=false);
   static glm::mat4 get_model_view_matrix();
   glm::mat4 get_mvp_for_shadow_map(const glm::vec3 &light_position) const;
   static glm::mat4 get_light_space_mvp(int light_index);
   static glm::vec3 get_world_space_eye_position();
   static glm::vec4 unproject(float z);
   static glm::vec4 unproject(float x, float y, float z);
   static glm::vec3 unproject_to_world_coordinates(const glm::vec3 &projected_coords);
   static glm::vec3 get_screen_y_uv();
   static glm::vec3 get_screen_x_uv();

   static glm::mat4 get_model_rotation(); // the quaterion from the mouse now rotates the model (not the view)
   static void setup_map_uniforms(Shader *shader_p, // in the draw loop
                                  const glm::mat4 &mvp,
                                  const glm::mat4 &view_rotation,
                                  float density_surface_opacity);
   static gboolean render(bool render_to_screendump_framebuffer_flag=false,
                          const std::string &output_file_name="coot-screendump.tga");
   static gboolean render_scene(); // like crows
   static gboolean render_scene_for_eye_internal(stereo_eye_t e); // like crows
   enum { PASS_TYPE_STANDARD, PASS_TYPE_GEN_SHADOW_MAP, PASS_TYPE_SSAO, PASS_TYPE_WITH_SHADOWS};
   static void render_scene_with_x_blur();
   static void render_scene_with_y_blur();
   static void render_scene_with_texture_combination_for_depth_blur();
   static void draw_map_molecules(stereo_eye_t eye, bool draw_transparent_maps);
   static void draw_map_molecules_with_shadows();
   static void draw_model_molecules(stereo_eye_t eye);
   static void draw_model_molecules_symmetry_with_shadows(stereo_eye_t eye);
   static void draw_intermediate_atoms(stereo_eye_t eye, unsigned int pass_type);
   static void draw_intermediate_atoms_rama_balls(stereo_eye_t eye, unsigned int pass_type);
   static void draw_molecule_atom_labels(molecule_class_info_t &m,
                                         stereo_eye_t eye,
                                         const glm::mat4 &mvp,
                                         const glm::mat4 &view_rotation);
   static void draw_hud_refinement_dialog_arrow_tab();
   static void draw_hud_colour_bar();
   static void draw_molecular_triangles();
   static void draw_molecules(stereo_eye_t eye);
   static void draw_meshes();
   static void draw_meshed_generic_display_object_meshes(stereo_eye_t eye, unsigned int pass_type);
   static void draw_molecules_other_meshes(stereo_eye_t eye, unsigned int pass_type);
   static void draw_instanced_meshes(stereo_eye_t eye);
   static void draw_unit_cells(stereo_eye_t eye);
   static void draw_cube(stereo_eye_t eye, GtkGLArea *glarea, unsigned int cube_type);
   static void draw_central_cube(stereo_eye_t eye, GtkGLArea *glarea);
   static void draw_origin_cube(stereo_eye_t eye, GtkGLArea *glarea);
   static void draw_rotation_centre_crosshairs(stereo_eye_t eye, GtkGLArea *glarea, unsigned int pass_type);
   static void draw_outlined_active_residue(stereo_eye_t eye);
   static void draw_hud_ligand_view();
   static void draw_hud_buttons();
   static void draw_hud_fps();
   static void draw_measure_distance_and_angles(stereo_eye_t eye);
   static void draw_ncs_ghosts();
   static std::list<std::chrono::time_point<std::chrono::high_resolution_clock> > frame_time_history_list;
   void set_do_ambient_occlusion(bool s) { shader_do_ambient_occlusion_flag = s; } // caller redraws

   static gl_rama_plot_t gl_rama_plot;
   static void draw_hud_ramachandran_plot(); // OpenGL rama plot
   void clear_gl_rama_plot();  // why not static also?
   static bool draw_gl_ramachandran_plot_flag;    // used for RSR vs translation or rotamers
   static bool draw_gl_ramachandran_plot_user_control_flag; // used for user not wanting to see it during RSR

   void reset_frame_buffers(int width, int height);
   void setup_lights();
   void translate_in_screen_x(float step_size);
   void translate_in_screen_z(float step_size);
   void move_forwards();
   void move_backwards();
   void step_screen_left();
   void step_screen_right();

   int blob_under_pointer_to_screen_centre();

   void display_next_map(); // one at a time, all, none.

   // make this private when the glarea render function is moved into graphics_info_t
   //
   static molecule_class_info_t moving_atoms_molecule; // used as a container for glsl variables. public access

   static bool regenerate_bonds_needs_make_bonds_type_checked_flag;
   void set_regenerate_bonds_needs_make_bonds_type_checked(bool state);
   bool get_regenerate_bonds_needs_make_bonds_type_checked_state();

   void sfcalc_genmap(int imol_model,
                      int imol_map_with_data_attached,
                      int imol_updating_difference_map);

   // this has a different API!
   coot::util::sfcalc_genmap_stats_t
   sfcalc_genmaps_using_bulk_solvent(int imol_model,
                                     int imol_map_with_data_attached,
                                     clipper::Xmap<float> *xmap_2fofc_p,
                                     clipper::Xmap<float> *xmap_fofc_p);
   static coot::util::sfcalc_genmap_stats_t latest_sfcalc_stats;

   static bool refinement_has_finished_moving_atoms_representation_update_needed_flag;

   static bool ignore_pseudo_zeros_for_map_stats;

   static std::atomic<bool> restraints_lock;
   static bool continue_threaded_refinement_loop; // so that the ESC key can stop the refinement
   static int  threaded_refinement_loop_counter;
   static int  threaded_refinement_loop_counter_bonds_gen;
   static bool threaded_refinement_needs_to_clear_up; // because ESC was pressed
   static bool threaded_refinement_needs_to_accept_moving_atoms; // because Return was pressed
   static int  threaded_refinement_redraw_timeout_fn_id; // -1 initially

   static int regenerate_intermediate_atoms_bonds_timeout_function();
   // static int regenerate_intermediate_atoms_bonds_timeout_function_and_draw(); old timeout style
   static gint regenerate_intermediate_atoms_bonds_timeout_function_and_draw(gpointer data);
   // we need to wait for the refinement to finish when we are in
   // immediate accept mode or no-gui.  In scripted (e.g. sphere-refine)
   // we should not wait
   void conditionally_wait_for_refinement_to_finish();

   static bool convert_dictionary_planes_to_improper_dihedrals_flag;
   void set_convert_dictionary_planes_to_improper_dihedrals(bool state);

   // for updating (difference) maps - we don't want to set 2 of these (or more) off
   //  at the same time.
   static std::atomic<bool> on_going_updating_map_lock;

   /*! \brief shiftfield B-factor refinement */
   void shiftfield_b_factor_refinement(int imol);

   /*! \brief shiftfield xyz refinement */
   void shiftfield_xyz_factor_refinement(int imol);

   // if pull_restraint_neighbour_displacement_max_radius < 1.5 (say) then
   // turn off proportional editing.
   static float pull_restraint_neighbour_displacement_max_radius;
   void pull_restraint_neighbour_displacement_change_max_radius(bool up_or_down); // change above
   static void draw_pull_restraint_neighbour_displacement_max_radius_circle();

   static void poke_the_refinement();

   // -------------------------------- Rail Points ------------------------------------------
   //! \name Rail Points!

   // 20230430-PE did I define the somewhere else also?
   static int updating_maps_imol_map;
   static int updating_maps_imol_diff_map;
   static std::vector<api::rail_points_t> rail_point_history; // map and model (model currently not used)

   //! calling this adds to the rail_points history. Make this pairs when we add model scoring.
   //! @returns the new rail points (since last modification)
   int calculate_new_rail_points(const updating_model_molecule_parameters_t &ummp);

   //! the total rail points
   //! @returns the sum of all rail points accumulated since the maps were connected.
   int rail_points_total() const;

   void updating_maps_update_the_coot_points_overlay();

   static guint updating_maps_timeout_function_idx;

   // --------------------------------------------------------------------------------------

   // by default, user-defined colours are on a colour wheel, but we can overwride that
   // by setting actual user defined colours for give colour indices
   //

   static bool draw_hud_colour_bar_flag;
   static std::vector<std::pair<unsigned int, coot::colour_holder> > user_defined_colours;
   // this function sets up the colour bar too and enables its drawing. It will need extra args for
   // the tick marks.
   static void set_user_defined_colours(const std::vector<std::pair<unsigned int, coot::colour_holder> > &user_defined_colours_in);
   static bool have_user_defined_colours();
   // run glColor3f())
   static void set_bond_colour_from_user_defined_colours(int icol);
   static void print_user_defined_colour_table();

#ifdef USE_PYTHON
   PyObject *pyobject_from_graphical_bonds_container(int imol,
						     const graphical_bonds_container &bonds_box) const;
   PyObject *get_intermediate_atoms_bonds_representation();
   PyObject *get_intermediate_atoms_distortions_py();
   PyObject *restraint_to_py(const coot::simple_restraint &restraint) const;
   PyObject *geometry_distortion_to_py(const coot::geometry_distortion_info_t &gd) const;
#endif

#ifdef USE_GUILE
   SCM geometry_distortion_to_scm(const coot::geometry_distortion_info_t &gd) const;
   SCM restraint_to_scm(const coot::simple_restraint &restraint) const;
#endif // USE_GUILE

#ifdef USE_PYTHON
   // Python function, called per frame draw - for Hamish
   static std::string python_draw_function_string;
   void set_python_draw_function(const std::string &f);
#endif // USE_PYTHON

   static ctpl::thread_pool static_thread_pool;

#ifdef USE_MOLECULES_TO_TRIANGLES
   static std::shared_ptr<Renderer> mol_tri_renderer;
   static std::shared_ptr<SceneSetup> mol_tri_scene_setup;
#endif // USE_MOLECULES_TO_TRIANGLES

   // text
   static std::map<GLchar, FT_character> ft_characters;
   void load_freetype_font_textures();
   static bool vera_font_loaded;

   // key-bindings
   static std::map<keyboard_key_t, key_bindings_t> key_bindings_map;
   static void add_key_binding(keyboard_key_t k, key_bindings_t kb) {
      key_bindings_map[k] = kb;
   }

   static void print_key_bindings();

   // GL IDs go here

   // intermediate atom pull restraints (note to self: intermediate atoms have
   // their own molecule_class_info_t now).
   // Generic Display objects should have their own class (that contains these IDs)
   // How about Rama balls and rota dodecs?
   // Where do contact dots fit in?
   // Overlap volumes?

   void setup_atom_pull_restraints_glsl();
   static GLuint m_VertexArray_for_pull_restraints_ID;
   static GLuint m_VertexBuffer_for_pull_restraints_ID;
   static GLuint m_IndexBuffer_for_atom_pull_restraints_ID;
   static unsigned int n_triangles_for_atom_pull_restraints; // index triangles, obviously
   static unsigned int n_vertices_for_atom_pull_restraints;

   // Lights
   //
   static std::map<unsigned int, lights_info_t> lights;

   void setup_hud_geometry_bars();
   static float hud_geometry_distortion_to_bar_size_nbc(float distortion);
   static float hud_geometry_distortion_to_bar_size_atom_pull(float distortion);
   static float hud_geometry_distortion_to_bar_size_rama(float distortion);
   static float hud_geometry_distortion_to_rotation_amount_rama(float distortion);

   void setup_hud_buttons();
   static HUDMesh mesh_for_hud_buttons;
   static std::vector<HUD_button_info_t> hud_button_info;


   void show_test_buttons(); // testing function
   // when the HUD buttons are shown, we want to show the atom pull buttons too
   // (they are in the GTK toolbar, not in OpenGL)

   void show_atom_pull_toolbar_buttons();
   void hide_atom_pull_toolbar_buttons();

   void show_accept_reject_hud_buttons();
   void reset_hud_buttons_size_and_position();
   void clear_hud_buttons(); // called by clear_up_moving_atoms_wrapper();

   static Mesh mesh_for_outline_of_active_residue;
   void update_mesh_for_outline_of_active_residue(int imol, const coot::atom_spec_t &spec, int n_press);
   static unsigned int outline_for_active_residue_frame_count;

   // Mesh mesh_for_particles("mesh-for-particles");
   // int n_particles = 100;
   static Mesh mesh_for_particles;
   static int n_particles;
   static particle_container_t particles;
   static std::vector<glm::vec3> get_particle_centre_positions();

   static TextureMesh tmesh_for_happy_face_residues_markers; // in 3d, like atom labels
   std::vector<glm::vec3> get_happy_face_residue_marker_positions();

   // these are "setup" by the function that starts them
   static LinesMesh lines_mesh_for_identification_pulse;
   static LinesMesh lines_mesh_for_generic_pulse; // loop through generic_pulse_centres
   static glm::vec3 identification_pulse_centre;
   static void draw_at_screen_centre_pulse(stereo_eye_t eye); // green sonar ping
   // these are variations of the (typically) multi-centre identification pulse.
   static void draw_generic_pulses(stereo_eye_t eye);
   static void draw_invalid_residue_pulse(stereo_eye_t eye);
   static void draw_delete_item_pulse(stereo_eye_t eye);
   static std::vector<glm::vec3> generic_pulse_centres;
   std::vector<glm::vec3> residue_to_positions(mmdb::Residue *residue_p) const;
   std::vector<glm::vec3> residue_to_side_chain_positions(mmdb::Residue *residue_p) const;
   void setup_delete_item_pulse(mmdb::Residue *residue_p);
   void setup_delete_residues_pulse(const std::vector<mmdb::Residue *> &residues);
   void setup_invalid_residue_pulse(mmdb::Residue *residue_p);
   static void pulse_marked_positions(const std::vector<glm::vec3> &positions, // generalization of above
                                      bool broken_lines_mode, unsigned int n_rings, float radius_overall,
                                      unsigned int n_ticks, const glm::vec4 &col, float resize_factor = 1.005f);

   static gboolean screen_centre_pulse_function(GtkWidget *widget,
                                                GdkFrameClock *frame_clock,
                                                gpointer data);
   static gboolean generic_pulse_function(GtkWidget *widget,
                                          GdkFrameClock *frame_clock,
                                          gpointer data);
   // this should wrap generic_pulse_function
   // return the continue-status
   static gboolean invalid_residue_pulse_function(GtkWidget *widget,
                                                  GdkFrameClock *frame_clock,
                                                  gpointer data);
   static gboolean wait_for_hooray_refinement_tick_func(GtkWidget *widget,
                                                        GdkFrameClock *frame_clock,
                                                        gpointer data);
   static int wait_for_hooray_refinement_tick_id;

   static Mesh mesh_for_boids; // with instancing
   static fun::boids_container_t boids;
   static LinesMesh lines_mesh_for_boids_box;
   void setup_draw_for_boids();

   static LinesMesh lines_mesh_for_hud_lines; // Used 3d but actually is just HUD lines, but there is no
                                              // HUDLinesMesh class - so I will use this (for now)

   // Let's base dynamic hydrogen bonds on how boids worked.
   static Mesh mesh_for_hydrogen_bonds; // with instancing, because dynamic
   void setup_draw_for_hydrogen_bonds();

   static float focus_blur_z_depth;
   static float focus_blur_strength;

   static std::chrono::time_point<std::chrono::high_resolution_clock> tick_hydrogen_bond_mesh_t_previous;

   static void add_a_tick();

   static gboolean glarea_tick_func(GtkWidget *widget,
                                    G_GNUC_UNUSED GdkFrameClock *frame_clock,
                                    G_GNUC_UNUSED gpointer data);

   static gboolean tick_function_is_active();

   static std::vector<std::pair<glm::vec3, glm::vec3> > hydrogen_bonds_atom_position_pairs;
   // this can be made more sophisticated later
   static void update_hydrogen_bond_mesh(const std::string &label);

   static int tick_function_id; // store the return value from gtk_widget_add_tick_callback()
   static bool do_tick_particles;
   static bool do_tick_spin;
   static bool do_tick_rock;
   static bool do_tick_boids;
   static bool do_tick_hydrogen_bonds_mesh;
   static bool do_tick_outline_for_active_residue;
   static bool do_tick_constant_draw;
   static bool do_tick_gone_diegos;
   static bool do_tick_gone_diff_map_peaks;

   static void fullscreen();
   static void unfullscreen();

   // 20220129-PE integrating crows

   static void render_scene_sans_depth_blur(stereo_eye_t eye,
                                     Shader *shader_for_tmeshes_p, Shader *shader_for_meshes_p,
                                     Shader *shader_for_tmeshes_with_shadows_p,
                                     Shader *shader_for_meshes_with_shadows_p,
                                     int width, int height);

   static void render_scene_with_depth_blur(stereo_eye_t eye,
                                            Shader *shader_for_tmeshes_p, Shader *shader_for_meshes_p,
                                            Shader *shader_for_tmeshes_with_shadows_p,
                                            Shader *shader_for_meshes_with_shadows_p,
                                            int width, int height);

   static Shader shader_for_meshes; // 20220208-PE temporary shader for molecules-as-meshes while I merge the crow code. No shadows.
   static Shader shader_for_tmeshes;
   static Shader shader_for_tmeshes_with_shadows;
   static Shader shader_for_meshes_with_shadows;
   static Shader shader_for_tmeshes_for_ssao; // render to white, don't use texture (like ssao_geometry)
   static Shader shader_for_meshes_for_ssao;  // render to white
   static Shader shader_for_instanced_meshes_for_ssao;
   static Shader shader_for_texture_meshes_shadow_map;
   static Shader shader_for_meshes_shadow_map;
   static Shader shader_for_instanced_meshes_shadow_map;
   static Shader shader_for_shadow_map_image_texture_mesh;
   static Shader shader_for_instanced_meshes_with_shadows;
   static Shader shader_for_effects; // colour balance or gamma ramp, say.
   static float shadow_box_size; // 20220402-PE needs to be big enough to cover the molecule. How big is that? I don't know how to
                                 // calculate it now, so let the user decide.
   // static bool show_shadow_map;
   enum { SHOW_SHADOW_MAP, SHOW_BASIC_SCENE, SHOW_SHADOWED_SCENE, SHOW_DEPTH_BLUR_SCENE, SHOW_AO_SCENE };
   static unsigned short int displayed_image_type;
   static HUDTextureMesh tmesh_for_shadow_map;
   static bool show_just_shadows; // show *just* the shadows in the texture-mesh-with-shadows shader

   void draw_models(stereo_eye_t eye,
                    Shader *shader_for_tmeshes_p,
                    Shader *shader_for_meshes_p,
                    Shader *shader_for_tmeshes_with_shadows_p,
                    Shader *shader_for_meshes_with_shadows_p,
                    int graphics_x_size,
                    int graphics_y_size,
                    bool draw_shadows=false,
                    float shadow_strength = 0.4,
                    bool show_just_shadows = false);

   void draw_models_with_shadows(stereo_eye_t eye,
				 Shader *shader_for_tmeshes_with_shadows_p,
				 Shader *shader_for_meshes_with_shadows_p,
				 int graphics_x_size,
				 int graphics_y_size,
				 bool draw_shadows=false,
				 float shadow_strength = 0.4,
				 bool show_just_shadows = false);

   void draw_Models_for_shadow_map(unsigned int light_index); // Models with a capital M (Models, not molecular models)
                                                              // Maybe use a better name?
   void draw_molecules_for_shadow_map(unsigned int light_index);

   static void draw_models_for_ssao();
   static void draw_molecules_for_ssao(stereo_eye_t eye);
   static void draw_molecules_with_shadows(stereo_eye_t eye); // use the above created shadow map to colour the pixels


   // DOF blur

   // static bool use_depth_blur_state;

   // void render_scene_with_texture_combination_for_depth_blur();
   // void render_scene_with_x_blur();
   // void render_scene_with_y_blur();

   static GtkAdjustment *focus_blur_z_depth_adjustment;
   // static float focus_blur_z_depth;
   // static float focus_blur_strength;
   // static GLuint blur_x_quad_vertex_array_id;
   //    static GLuint blur_y_quad_vertex_array_id;
   // static GLuint combine_textures_using_depth_quad_vertex_array_id;
   // static Shader shader_for_x_blur;
   // static Shader shader_for_y_blur;
   // static Shader shader_for_dof_blur_by_texture_combination;
   // static framebuffer combine_textures_using_depth_framebuffer;
   // static framebuffer blur_x_framebuffer;

   static GLuint screen_AO_quad_vertex_array_id;
   static GLuint screen_AO_quad_VBO;

   static unsigned int shadow_texture_multiplier;
   static unsigned int shadow_texture_width;  //  = 4 * 1024; // too big?      // derived from the above (n x 1024)
   static unsigned int shadow_texture_height; //  = 4 * 1024;
   void set_shadow_texture_resolution_multiplier(unsigned int m);
   static float shadow_strength;
   static unsigned int shadow_softness; // 1, 2 or 3
   static unsigned int shadow_depthMap_framebuffer; // change this to a real framebuffer
                                             // when things are working
   static unsigned int shadow_depthMap_texture; // for the framebuffer texture

   static unsigned int rboDepth;
   static unsigned int n_ssao_kernel_samples;
   static void generate_ssao_kernel_samples();
   // static framebuffer blur_y_framebuffer;
   static framebuffer framebuffer_for_ssao_gbuffer;
   static framebuffer framebuffer_for_ssao;
   static framebuffer framebuffer_for_ssao_blur;
   static framebuffer framebuffer_for_effects;
   enum { EFFECTS_SHADER_STANDARD, EFFECTS_SHADER_INPUT_TEXTURE, EFFECTS_SHADER_SSAO, EFFECTS_SHADER_DEPTH};
   static unsigned int effects_shader_output_type;

   // AO new try (dec-2021)

   static void draw_hud_elements();
   static void render_3d_scene(GtkGLArea *gl_area, stereo_eye_t eye);
   static void render_3d_scene_for_ssao(); // c.f. above, this one doesn't pass the gl_area.
                                           // I don't know which is best.
   static void render_3d_scene_with_shadows(stereo_eye_t eye); // change the shader from meshes.shader to meshes-with-shadows.shader
   void init_joey_ssao_stuff(int w, int h);
   void read_some_test_models();

   static std::vector<Model> models; // from crows 20220129-PE
   static void add_model(const Model &model);
   void scale_model(unsigned int idx, float scale_factor);

   // static unsigned int gBufferFBO;
   static unsigned int ssaoFBO;
   static unsigned int ssaoBlurFBO;
   static Shader shaderGeometryPass;
   static Shader shaderSSAO;
   static Shader shaderSSAOBlur;

   // now in the framebuffer
   // static unsigned int gPosition;
   // static unsigned int gNormal;
   // static unsigned int gAlbedo;
   static unsigned int noiseTexture;
   static unsigned int ssaoColorBuffer;
   static unsigned int ssaoColorBufferBlur;
   static unsigned int ssao_blur_size; // in the blur shader, how big should the box be
                                       // 0,1,2.

   // now uniforms
   static float SSAO_radius;
   static float SSAO_bias;
   static bool use_ssao;
   static float ssao_strength;
   static bool show_just_ssao; // show *just* the SSAO colour in the effects shader

   static float effects_brightness;
   static float effects_gamma;

   static void renderQuad();
   static std::vector<glm::vec3> ssaoKernel;
   // Camera camera(glm::vec3(0.0f, 0.0f, 5.0f));
   // static Camera camera;

   static bool show_fps_flag;
   static float fps_times_scale_factor;
   static glm::mat4 get_projection_matrix(bool do_orthographic_projection,
                                          int graphics_x_size, int graphics_y_size);

   static glm::mat4 get_view_matrix(); // the lookAt() matrix
   static glm::mat4 get_model_matrix();

   void resize_framebuffers_textures_renderbuffers(int width, int height);

   void read_test_gltf_models();

   //! return the model index
   int load_gltf_model(const std::string &gltf_file_name);

   //! \brief set the model animation parameters
   void set_model_animation_parameters(unsigned int model_index, float amplitude, float wave_numer, float freq) {
      if (model_index < models.size()) {
         auto &model = models[model_index];
         model.set_animation_parameters(amplitude, wave_numer, freq);
      }
   }

   //! \brief enable/disable the model animation (on or off)
   void set_model_animation_state(unsigned int model_index, bool state);

   static std::string stringify_error_code(GLenum err);

   static void attach_buffers(const char *s = __builtin_FUNCTION()) ;

   /// List of label strings (col 0) and model indices (int) (col 1)
   static GtkListStore* validation_graph_model_list;
   /// Should be called when a model gets added or deleted.
   /// Updates the GtkListStore for the validation graph model combobox
   static void refresh_validation_graph_model_list();
   // 20230420-PE and just like the above, let's have a model that we can use for the Ramachandran plot molecule chooser
   static GtkListStore* ramachandran_plot_model_list;
   // like the above, this should be called when a model gets added or deleted.
   static void refresh_ramachandran_plot_model_list();

   /// -1 if none
   static int active_validation_graph_model_idx;

   private:
   static std::string active_validation_graph_chain_id;
   typedef std::map<coot::validation_graph_type,GtkWidget*> validation_graph_map_t;
   static validation_graph_map_t validation_graph_widgets;
   typedef std::map<coot::validation_graph_type,std::shared_ptr<coot::validation_information_t>> validation_data_map_t;
   static validation_data_map_t validation_graph_data;
   public:
   static void update_active_validation_graph_model(int new_model_idx);
   static void change_validation_graph_chain(const std::string& chain_id);

   static void create_validation_graph(int imol, coot::validation_graph_type type);
   static void destroy_validation_graph(coot::validation_graph_type type);

   // 20230417-PE functions to fill the validation information for the new valiadtionn graphs
   coot::validation_information_t get_validation_data_for_geometry_analysis(int imol);

   // "Coot: " will be prepended to the dialog title before use
   void fill_generic_validation_box_of_buttons(const std::string &dialog_title,
                                               const std::vector<labelled_button_info_t> &v);

   void fill_atoms_with_zero_occupancy_box_of_buttons(const std::vector<labelled_button_info_t> &lbv);



   // 20230419-PE ----- a holder for the OpenGL-based Ramachandran Plots
   // each rama plot is held in a GtkBox - and that box has a "Close" button
   // which destroys the plot - and removes it from this vector
   //
   class widgeted_rama_plot_t {
   public:
      int imol;
      gl_rama_plot_t rama;
      std::string residue_selection;
      GtkWidget *gtk_gl_area;
      GtkWidget *close_button;
      GtkWidget *box;
      widgeted_rama_plot_t(int imol, const std::string &residue_selection,
                           const gl_rama_plot_t &rama, GtkWidget *gtk_gl_area, GtkWidget *button, GtkWidget *box) :
         imol(imol), rama(rama), residue_selection(residue_selection),
         gtk_gl_area(gtk_gl_area), close_button(button), box(box) {}
      bool matches_gl_area(GtkWidget *w) const { return (w == gtk_gl_area); }
   };
   static std::vector<widgeted_rama_plot_t> rama_plot_boxes;
   static void draw_rama_plots(); // draw the rama plots in the above vector
   static void remove_plot_from_rama_plots(GtkWidget *rama_plot);
   static void rama_plot_boxes_handle_close_molecule(int imol);
   static void rama_plot_boxes_handle_molecule_update(int imol);
   static void rama_plot_boxes_handle_molecule_update(GtkWidget *box);
   static void rama_plot_boxes_handle_molecule_update(GtkWidget *box, const std::string &entry_string); // used in the residue selection
                                                                                                        // entry change callback

   static void add_shortcuts_to_window(GtkWidget *shortcuts_window);

   static bool noughties_physics; // false by default

   static float gaussian_surface_sigma;
   static float gaussian_surface_contour_level;
   static float gaussian_surface_box_radius;
   static float gaussian_surface_grid_scale;
   static float gaussian_surface_fft_b_factor;
   static short int gaussian_surface_chain_colour_mode;

   static std::pair<bool, std::string> servalcat_fofc;
   static std::pair<bool, std::string> servalcat_refine; // output "pdb" file name

   static std::pair<bool, std::string> acedrg_link;
   static bool acedrg_running; // not link acedrg - this is acedrg from CCD

   static bool curmudgeon_mode; // default false, particles and faces
   static bool use_sounds; // default true

   static std::vector<std::pair<std::string, clipper::Xmap<float> > > map_partition_results;
   static int map_partition_results_state;
   static std::string map_partition_results_state_string; // "Done A Chain" etc.

   static unsigned int logging_line_index;

   static coot::inchikey_store_t inchikey_store;
   static void read_inchikeys();

   // 20251219-PE Do I need this?
   static std::string current_alt_conf; // nobody wants this "per molecule" right?

   // add a pumpkin as a graphics object and draw it.
   void pumpkin();

   // tomogram stuff
   void set_tomo_section_view_section(int imol, int section_index);
   static bool tomo_picker_flag;
   class tomo_view_info_t {
      public:
      tomo_view_info_t() : imol(-1), section_index(-1), axis(-1) {}
      tomo_view_info_t(int imol, coot::Cell c, int section, int axis_in) :
         imol(imol), cell(c), section_index(section), axis(axis_in) {}
      int imol;
      coot::Cell cell;
      int section_index;
      int axis;
      float get_P_z() { return cell.c * static_cast<float>(section_index)/static_cast<float>(molecules[imol].xmap.grid_sampling().nw()); }
   };
   static tomo_view_info_t tomo_view_info;

};


class molecule_rot_t {

 public:
   static float x_axis_angle;
   static float y_axis_angle;
};


void do_accept_reject_dialog(std::string fit_type, const coot::refinement_results_t &ref_results);
// old real GTK dialog interface:
void do_accept_reject_dialog_with_a_dialog(std::string fit_type, const coot::refinement_results_t &ref_results);
// new OpenGL buttons interface
void do_accept_reject_hud_buttons(std::string fit_type, const coot::refinement_results_t &ref_results);

void add_accept_reject_lights(GtkWidget *window, const coot::refinement_results_t &ref_results);
// return a pointer to a "new" object
GdkRGBA colour_by_distortion(float dist);
GdkRGBA colour_by_rama_plot_distortion(float plot_value, int rama_plot_type);
void set_colour_accept_reject_event_box(GtkWidget *eventbox, GdkRGBA *col);
GtkWidget *wrapped_create_accept_reject_refinement_dialog();
void update_accept_reject_dialog_with_results(GtkWidget *accept_reject_dialog,
					      coot::accept_reject_text_type text_type,
					      const coot::refinement_results_t &rr);
GtkWidget *wrapped_create_multi_residue_torsion_dialog(const std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > &pairs);

// Some currently useless Perspective View definition
//
#define VIEW_ASPECT 1.3


#endif // GRAPHICS_INFO_H
