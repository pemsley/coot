/* src/callbacks.h
 * 
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include <gtk/gtk.h>

typedef const char entry_char_type;

void
on_window1_destroy                     (GtkObject       *object,
                                        gpointer         user_data);

void
on_file1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_open_coordinates1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_open_dataset1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_exit1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_clipping1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_coords_fileselection1_destroy       (GtkObject       *object,
                                        gpointer         user_data);

void
on_dataset_fileselection1_destroy      (GtkObject       *object,
                                        gpointer         user_data);

void
on_ok_button_coordinates_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_ok_button_dataset_clicked           (GtkButton       *button,
                                        gpointer         user_data);

void
on_cancel_button1_clicked              (GtkButton       *button,
                                        gpointer         user_data);

void
on_cancel_button1_clicked              (GtkButton       *button,
                                        gpointer         user_data);

void
on_cancel_button2_clicked              (GtkButton       *button,
                                        gpointer         user_data);

void
on_dataset_fileselection1_destroy      (GtkObject       *object,
                                        gpointer         user_data);

void
on_cancel_coords_button1_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_cancel_dataset_button1_clicked      (GtkButton       *button,
                                        gpointer         user_data);

void
on_column_label_ok_button_clicked      (GtkButton       *button,
                                        gpointer         user_data);

void
on_button2_clicked                     (GtkButton       *button,
                                        gpointer         user_data);

void
on_column_label_cancel_button_clicked  (GtkButton       *button,
                                        gpointer         user_data);

void
on_about1_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_clipping_button_clicked             (GtkButton       *button,
                                        gpointer         user_data);

void
on_density_ok_button_clicked           (GtkButton       *button,
                                        gpointer         user_data);

void
on_density_ok_cancel_clicked           (GtkButton       *button,
                                        gpointer         user_data);

void
on_density_size1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_density_cancel_clicked              (GtkButton       *button,
                                        gpointer         user_data);

gboolean
on_hscale1_key_release_event           (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data);

gboolean
on_hscale1_key_press_event             (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data);

gboolean
on_hscale1_button_press_event          (GtkWidget       *widget,
                                        GdkEventButton  *event,
                                        gpointer         user_data);

gboolean
on_hscale1_button_release_event        (GtkWidget       *widget,
                                        GdkEventButton  *event,
                                        gpointer         user_data);

void
on_fps1_activate                       (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_fps_window_ok_button_clicked        (GtkButton       *button,
                                        gpointer         user_data);

void
on_active_map_ok_button_clicked        (GtkButton       *button,
                                        gpointer         user_data);

void
on_dragged_map1_activate               (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_map_colour1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_item1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_show_symmetry_ok_button_pressed     (GtkButton       *button,
                                        gpointer         user_data);

void
on_show_symmetry1_activate             (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_button2_clicked                     (GtkButton       *button,
                                        gpointer         user_data);

void
on_show_symmetry_ok_button_clicked     (GtkButton       *button,
                                        gpointer         user_data);

void
on_about_ok_button_clicked             (GtkButton       *button,
                                        gpointer         user_data);

void
on_symmetry_colour_patch_button_clicked (GtkButton       *button,
                                        gpointer         user_data);

/* void  */
/* symmetry_colour_adjustment_changed (GtkAdjustment *adj,  */
/* 				    GtkWidget *window); */

void
on_anisotropic_atoms1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_show_aniso_ok_button_clicked        (GtkButton       *button,
                                        gpointer         user_data);

void aniso_probability_adjustment_changed(GtkAdjustment *adj, 
					  GtkWidget *window);

void 
on_smooth_scrolling_window_ok_button_clicked (GtkButton       *button,
					      gpointer         user_data);

void
on_recentring1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_font_size1_activate                 (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_button1_clicked                     (GtkButton       *button,
                                        gpointer         user_data);

void
on_font_size_ok_button_clicked         (GtkButton       *button,
                                        gpointer         user_data);

void
on_greer_skeleton1_activate            (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_cowtan_foadi_skeleton1_activate     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_these_are_placeholders1_activate    (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_mapcolourmap1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_map1_activate                       (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_map2_activate                       (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_map1_activate                       (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_map2_activate                       (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_attach_scroll_wheel_to_which_map_1_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_greer_on_activate                   (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_greer_off_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_foadi_on_activate                   (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_foadi_off_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_open_map1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_ok_button_map_name_clicked          (GtkButton       *button,
                                        gpointer         user_data);

void
on_cancel_button_map_name_clicked      (GtkButton       *button,
                                        gpointer         user_data);

void
on_vt_flat1_activate                   (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_vt_spherical_surface1_activate      (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_skeleton_colour1_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_fps1_activate                       (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_phs_info_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data);

void
on_phs_info_cancel_button_clicked      (GtkButton       *button,
                                        gpointer         user_data);

void
on_ok_phs_coord_button_clicked         (GtkButton       *button,
                                        gpointer         user_data);

void
on_cancel_phs_coord_button_clicked     (GtkButton       *button,
                                        gpointer         user_data);

void
on_map_and_mol_control1_activate       (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_display_only_active1_activate       (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_go_to_atom_ok_button_clicked        (GtkButton       *button,
                                        gpointer         user_data);

void
on_go_to_atom_cancel_button_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_go_to_atom1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_go_to_atom_next_residue_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_go_to_atom_previous_residue_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_skeleton_box_radius1_activate       (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_skeletonization_level1_activate     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_skeleton_box_size_ok_button_clicked (GtkButton       *button,
                                        gpointer         user_data);

void
on_skel_box_radius_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data);

void
on_skel_box_radius_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_skeletonization_level_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_skeletonization_level_apply_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_skeletonization_level_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_on1_activate                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_off1_activate                       (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_autobuild_ca_on_activate            (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_autobuild_ca_off_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_display_control_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data);


void
on_display_control_window_glade_destroy
                                        (GtkObject       *object,
                                        gpointer         user_data);

void
on_rotation_centre_size_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_pink_pointer_size1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_phs_cell_choice_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data);

void
on_phs_cell_choice_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_debug_testing_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_test_thing1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_regularize1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_scripting_window_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_get_pdb_using_code1_activate        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_get_pdb_and_sf_using_code1_activate (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

gboolean
on_accession_code_entry_key_press_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data);

void
on_dynarama_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data);

void
on_ramachandran_plot1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_background_black1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_background_white1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_dynarama_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data);

void
on_dynarama_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data);

void
on_dynarama_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data);

void
on_go_to_atom_apply_button_clicked     (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligands1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_find_ligand_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligand_cancel_button_clicked   (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligand_many_atoms_continue_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligand_many_atoms_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligands1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_find_ligand_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligand_cancel_button_clicked   (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligand_many_atoms_continue_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligand_many_atoms_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
_unimplemented_                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_find_ligands1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_find_ligand_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligand_cancel_button_clicked   (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligand_many_atoms_continue_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligand_many_atoms_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_test_thing1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_find_ligands1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_find_ligand_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligand_cancel_button_clicked   (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligand_many_atoms_continue_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligand_many_atoms_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_prune_and_colour1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_find_ligands1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_find_ligand_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligand_cancel_button_clicked   (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligand_many_atoms_continue_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligand_many_atoms_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligands1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_find_ligand_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligand_cancel_button_clicked   (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligand_many_atoms_continue_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligand_many_atoms_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_info1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_model_refine_activate               (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

/* void */
/* on_model_refine_regularize_zone_button_clicked */
/*                                         (GtkButton       *button, */
/*                                         gpointer         user_data); */

void
on_model_refine_dialog_dismiss_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_save_coordinates1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_save_coords_dialog_save_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_save_coord_ok_button_clicked        (GtkButton       *button,
                                        gpointer         user_data);

void
on_save_coords_cancel_button_clicked   (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_refine_togglebutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_regularize_togglebutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_fixed_atoms_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_rigid_body_togglebutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_refine_params_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_refine_params_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_goto_atom_window_destroy            (GtkObject       *object,
                                        gpointer         user_data);

void
on_distance1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_angle1_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_model_refine_dialog_pepflip_togglebutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_accept_reject_refinement_accept_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_accept_reject_refinement_reject_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_accept_reject_refinement_docked_accept_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_accept_reject_refinement_docked_reject_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_accept_reject_dialog_frame_docked_hide
                                        (GtkWidget       *widget,
                                        gpointer         user_data);

void
on_accept_reject_refinement_atom_pull_autoclear_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_accept_reject_atom_pull_clear_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_find_waters_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_waters_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_waters_cancel_button_clicked   (GtkButton       *button,
                                        gpointer         user_data);

void
on_fast_sss_dialog_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_fast_sss_dialog_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data);


void
on_environment_distances1_activate     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_enviroment_distance_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_enviroment_distance_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

gboolean
on_environment_distance_min_entry_key_press_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data);

gboolean
on_environment_distance_max_entry_key_press_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data);


void
on_refine_params_use_torsions_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_refine_params_use_initial_pos_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_map_dynamic_map_sampling_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_map_dynamic_map_size_display_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_refine_params_use_peptide_torsions_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_import_cif_dictionary1_activate     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_cif_fileselection_ok_button_clicked (GtkButton       *button,
                                        gpointer         user_data);

void
on_cif_fileselection_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_cif_dictionary_fileselection_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_cif_dictionary_fileselection_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_fit_terminal_residue_togglebutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_type_chooser_ALA_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_type_chooser_ARG_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_type_chooser_ASN_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_type_chooser_ASP_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_type_chooser_CYS_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_type_chooser_GLN_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_type_chooser_GLU_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_type_chooser_GLY_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_type_chooser_HIS_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_type_chooser_ILE_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_type_chooser_LEU_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_type_chooser_LYS_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_type_chooser_MET_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_type_chooser_PHE_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_type_chooser_PRO_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_type_chooser_SER_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_type_chooser_THR_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_type_chooser_TRP_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_type_chooser_TYR_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_type_chooser_VAL_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_rot_trans_togglebutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_rot_trans_by_residue_range_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_model_refine_dialog_rot_trans_by_chain_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_model_refine_dialog_rot_trans_by_molecule_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_rotate_translate_obj_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_rotate_translate_obj_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_run_script1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_run_script_fileselection_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_run_script_fileselection_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_db_main_togglebutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_delete_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_delete_item_residue_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_delete_item_atom_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_close_molecule1_activate            (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_close_molecule_close_button_clicked (GtkButton       *button,
                                        gpointer         user_data);

void
on_close_molecule_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_delete_item_cancel_button_clicked   (GtkButton       *button,
                                        gpointer         user_data);

void
on_reset_view_button_clicked           (GtkButton       *button,
                                        gpointer         user_data);

void
on_reset_view_toolbutton_clicked       (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_hints1_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_residue_info_ok_button_clicked      (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_info_cancel_button_clicked  (GtkButton       *button,
                                        gpointer         user_data);

void
on_hints_dialog_ok_button_clicked      (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_find_ligands_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_rotamer_selection_ok_button_clicked (GtkButton       *button,
                                        gpointer         user_data);

void
on_rotamer_selection_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_rotamer_selection_button_rot_0_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_rotamer_togglebutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_mutate_togglebutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data);

gboolean
on_go_to_atom_chain_entry_key_press_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data);

gboolean
on_go_to_atom_residue_entry_key_press_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data);

gboolean
on_go_to_atom_atom_name_entry_key_press_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data);

void
on_model_refine_dialog_pointer_atom_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_unsaved_changes_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_unsaved_changes_continue_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_environment_distance_label_atom_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_baton_accept_button_clicked         (GtkButton       *button,
                                        gpointer         user_data);

void
on_baton_try_again_button_clicked      (GtkButton       *button,
                                        gpointer         user_data);

void
on_baton_tip_previous_button_clicked      (GtkButton       *button,
                                        gpointer         user_data);

void
on_baton_shorten_button_clicked        (GtkButton       *button,
                                        gpointer         user_data);

void
on_baton_dialog_ok_button_clicked      (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_baton_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_baton_lengthen_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_use_weights_checkbutton_toggled     (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_environment_distance_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_environment_distance_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_show_symmetry_as_calphas_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

/* void */
/* on_go_to_atom_residue_list_select_child */
/*                                         (GtkList         *list, */
/*                                         GtkWidget       *widget, */
/*                                         gpointer         user_data); */

/* void */
/* on_go_to_atom_residue_list_selection_changed */
/*                                         (GtkList         *list, */
/*                                         gpointer         user_data); */

/* void */
/* on_go_to_atom_residue_list_unselect_child */
/*                                         (GtkList         *list, */
/*                                         GtkWidget       *widget, */
/*                                         gpointer         user_data); */

void
on_coordinates_recentring1_activate    (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_read_pdb_recentre_yes_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_read_pdb_recentre_no_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_read_pdb_recentre_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_pointer_atom_type_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_pointer_atom_type_ok_button_clicked (GtkButton       *button,
                                        gpointer         user_data);

void
on_refmac_column_labels_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_run_refmac_run_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_run_refmac_cancel_button_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_refmac_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_single_map_properties_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_single_map_properties_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_single_map_properties_colour_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_single_map_properties_colour_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_run_refmac_phase_input_optionmenu_changed
                                        (GtkOptionMenu   *optionmenu,
                                        gpointer         user_data);

void
on_run_refmac_tls_checkbutton_toggled  (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_run_refmac_twin_checkbutton_toggled (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_run_refmac_sad_checkbutton_toggled  (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_run_refmac_map_mtz_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_run_refmac_mtz_file_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_run_refmac_mtz_filechooserdialog_response
                                        (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data);

void
on_run_refmac_mtz_filechooserdialog_destroy
                                        (GtkObject       *object,
                                        gpointer         user_data);

void
on_run_refmac_mtz_filechooser_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_run_refmac_file_help_button_clicked (GtkButton       *button,
                                        gpointer         user_data);

void
on_run_refmac_file_help_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_run_refmac_sad_help_button_clicked  (GtkButton       *button,
                                        gpointer         user_data);

void
on_run_refmac_sad_help_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_run_refmac_nolabels_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_run_refmac_map_optionmenu_changed   (GtkOptionMenu   *optionmenu,
                                        gpointer         user_data);

void
on_run_refmac_diff_map_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_run_refmac_ncs_checkbutton_toggled  (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_run_refmac_phase_combine_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_baton_undo_button_clicked           (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_undo_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_undo_molecule_chooser_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_clear_pending_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_skeleton_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data);

void
on_skeleton_cancel_button_clicked      (GtkButton       *button,
                                        gpointer         user_data);

void
on_skeleton1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_virtual_trackball_menu_pops         (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_model_refine_dialog_auto_fit_rotamer_togglebutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_residue_info2_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_import_all_dictionary_cifs1_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_model_refine_dialog_destroy         (GtkObject       *object,
                                        gpointer         user_data);

void
on_residue_info_apply_all_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_residue_parameters1_activate        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_residue_info_dialog_destroy         (GtkObject       *object,
                                        gpointer         user_data);

void
on_crosshairs1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_crosshairs_on_radiobutton_toggled   (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_crosshairs_off_radiobutton_toggled  (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_display_crosshairs_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_mutate_auto_fit_togglebutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_add_alt_conf_cb_radiobutton_toggled (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_add_alt_conf_residue_range_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_add_alt_conf_cancel_button_clicked  (GtkButton       *button,
                                        gpointer         user_data);


void
on_validation_dialog_next_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_validation_dialog_prev_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_validation_dialog_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_edit_phi_psi_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_dynarama_cancel_button_clicked      (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_redo_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_edit_chi_angles_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

gboolean
on_window1_configure_event             (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data);

gboolean
on_window1_configure_event             (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data);

void
on_run_state_file_ok_button_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_run_state_file_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_edit_backbone_torsions_dialog_destroy
                                        (GtkObject       *object,
                                        gpointer         user_data);

void
on_edit_backbone_torsion_rotate_peptide_button_pressed
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_edit_backbone_torsion_rotate_peptide_button_released
                                        (GtkButton       *button,
                                        gpointer         user_data);

gboolean
on_edit_backbone_torsion_rotate_peptide_button_motion_notify_event
                                        (GtkWidget       *widget,
                                        GdkEventMotion  *event,
                                        gpointer         user_data);

void
on_edit_backbone_torsion_rotate_peptide_carbonyl_button_pressed
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_edit_backbone_torsion_rotate_peptide_carbonyl_button_released
                                        (GtkButton       *button,
                                        gpointer         user_data);

gboolean
on_edit_backbone_torsion_rotate_peptide_carbonyl_button_motion_notify_event
                                        (GtkWidget       *widget,
                                        GdkEventMotion  *event,
                                        gpointer         user_data);

void
on_model_refine_dialog_edit_backbone_torsions_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_edit_backbone_torsions_button_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_button6_clicked                     (GtkButton       *button,
                                        gpointer         user_data);

void
on_edit_backbone_torsion_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_edit_backbone_torsion_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_sequence_view1_activate             (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_clear_simple_distances2_activate    (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_workflow_radiobutton_coords_only_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_workflow_radiobutton_coords_and_dataset_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_workflow_radiobutton_coords_and_map_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_workflow_radiobutton_map_from_dataset_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_workstate_radiobutton_clean_slate_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_workflow_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data);

void
on_workflow_cancel_button_clicked      (GtkButton       *button,
                                        gpointer         user_data);

void
on_select_map_for_fitting_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_map_select_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_clear_atom_labels1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_add_alt_conf_ca_radiobutton_toggled (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_add_alt_conf_whole_single_residue_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_model_refine_dialog_add_alt_conf_togglebutton_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_add_alt_conf_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_add_alt_conf_dialog_destroy         (GtkObject       *object,
                                        gpointer         user_data);

void
on_run_refmac_help_button_clicked      (GtkButton       *button,
                                        gpointer         user_data);

void
on_run_refmac_help_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_run_refmac_nolabels_help_button_clicked      
					(GtkButton       *button,
                                        gpointer         user_data);

void
on_run_refmac_nolabels_help_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_no_restraints_info_dialog_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_no_cif_dictionary_bonds_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_pointer_atom_type_other_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_ligand_big_blob_dismiss_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_edit_chi_angles_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_edit_chi_angles_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_check_waters_ok_button_clicked      (GtkButton       *button,
                                        gpointer         user_data);

void
on_check_waters_cancel_button_clicked  (GtkButton       *button,
                                        gpointer         user_data);

void
on_edit_chi_angles_normal_rotation_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_geometry_distance_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_geometry_clear_last_distance_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_geometry_clear_all_distances_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_geometry_clear_atom_labels_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_geometry_dialog_close_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_geometry_angle_togglebutton_toggled (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_distances_and_angles1_activate      (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_geometry_dialog_destroy             (GtkObject       *object,
                                        gpointer         user_data);

void
on_new_ligands_info_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_no_new_ligands_info_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_zoom1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_zoom_dialog_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_edit_chi_angles_dialog_destroy      (GtkObject       *object,
                                        gpointer         user_data);

void
on_check_waters_low_occ_dist_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_check_waters_zero_occ_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);
void
on_check_waters_by_difference_map_active_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
					 gpointer         user_data);

void
on_get_monomer1_activate               (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_libcheck_monomer_ok_button_clicked  (GtkButton       *button,
                                        gpointer         user_data);

void
on_libcheck_monomer_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

gboolean
on_libcheck_monomer_entry_key_press_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data);

void
on_recover_coordinates_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_recover_coordinates_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_recover_session1_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_centre_atom_label1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_centre_atom_label_ok_button_clicked (GtkButton       *button,
                                        gpointer         user_data);

void
on_edit_chi_angles_help_button_clicked (GtkButton       *button,
                                        gpointer         user_data);

void
on_help_chi_angles_dismiss_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_rotate_translate_obj_dialog_destroy (GtkObject       *object,
                                        gpointer         user_data);

void
on_no_symmetry_warning_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_nothing_to_recover_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_superpose_dialog_superpose_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_superpose_dialog_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_superpose_nonsense_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_ssm_superposition1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_clipping1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_add_terminal_residue_finds_none_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_single_map_sigma_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_no_bad_chiral_volumes_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_check_chiral_volumes_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_check_chiral_volumes_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_bad_chiral_atoms1_activate     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_chiral_volume_baddies_dialog_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_dynarama_window_destroy             (GtkObject       *object,
                                        gpointer         user_data);

void
on_delete_item_residue_hydrogens_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_rigid_body_refinement_failed_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_baton_mode_calculate_skeleton_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_baton_mode_calculate_skeleton_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_column_label_expert_mode_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_columns_label_use_resolution_limits_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_column_labels_use_resolution_limits_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_auto_open_mtz_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_merge_molecules_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data);

void
on_merge_molecules_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_mutate_sequence_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data);

void
on_mutate_sequence_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_merge_molecules1_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_mutate_molecule1_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_draw_hydrogens_yes_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_draw_hydrogens_no_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_bond_parameters_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data);

void
on_renumber_residue_range_radiobutton_1_toggled
                                        (GtkToggleButton *togglebutton,
					 gpointer         user_data);

void
on_renumber_residue_range_radiobutton_3_toggled
                                        (GtkToggleButton *togglebutton,
					 gpointer         user_data);
void
on_renumber_residue_range_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_renumber_residue_range_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_add_OXT_c_terminus_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_add_OXT_residue_radiobutton_toggled (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_add_OXT_ok_button_clicked           (GtkButton       *button,
                                        gpointer         user_data);

void
on_add_OXT_cancel_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_add_OXT_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_renumber_residues1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_bond_parameters1_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_bond_parameters_apply_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_bond_parameters_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data);

void
on_bond_parameters_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_background_colour1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_renumber_residues1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_ligand_no_blobs_OK_button_clicked   (GtkButton       *button,
                                        gpointer         user_data);

void
on_new_delete_molecules_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_new_delete_molecules_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_ramachandran_plot2_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_incorrect_chiral_volumes1_activate  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_unmodelled_blobs1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_find_blobs_ok_button_clicked        (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_blobs_cancel_button_clicked    (GtkButton       *button,
                                        gpointer         user_data);


void
on_chiral_restraints_problem_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_check_waters_diff_map_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_check_waters_diff_map_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_check_waters_by_difference_map_variance1_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_interesting_waters_by_difference_map_check_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_nothing_bad_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_skeletonize_map_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_skeletonize_map_dialog_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_antialiasing1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_antialias_dialog_yes_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_antialias_dialog_no_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_antialiasing_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_geometry_graphs_ok_button_clicked   (GtkButton       *button,
                                        gpointer         user_data);

void
on_geometry_graphs_dialog_destroy      (GtkObject       *object,
                                        gpointer         user_data);

void
on_save_symmetry_coords_fileselection_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_save_symmetry_coords_fileselection_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_save_symmetry_coordinates1_activate (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_experimental1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_geometry_analysis1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_peptide_omega_analysis1_activate    (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_temp_fact_variance_analysis1_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);
////B 
void
on_temp_fact_analysis1_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);
////E 
void
on_rotamer_analysis1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_density_fit_analysis1_activate      (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_stereo1_activate                    (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_stereo_dialog_mon_button_toggled    (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_stereo_dialog_hardware_stereo_button_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_stereo_dialog_ok_button_clicked     (GtkButton       *button,
                                        gpointer         user_data);

void
on_stereo_dialog_mono_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_stereo_dialog_hardware_stereo_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_stereo_dialog_zalman_stereo_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_residue_type_chooser_MSE_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_preferences1_activate               (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_preferences_general_radiotoolbutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data);

void
on_preferences_bond_radiotoolbutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data);

void
on_preferences_map_radiotoolbutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data);

void
on_preferences_geometry_radiotoolbutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data);

void
on_preferences_colour_radiotoolbutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data);

void
on_preferences_other_radiotoolbutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data);

void
on_preferences_geometry_cis_peptide_bad_yes_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_geometry_cis_peptide_bad_no_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_bond_colours_hscale_value_changed
                                        (GtkRange        *range,
                                        gpointer         user_data);

void
on_preferences_bond_colours_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_bg_colour_black_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_bg_colour_white_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_bg_colour_own_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_bg_colour_colorbutton_color_set
                                        (GtkColorButton  *colorbutton,
                                        gpointer         user_data);

void
on_preferences_bg_colour_colorbutton_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_preferences_map_radius_entry_activate
                                        (GtkEntry        *entry,
                                        gpointer         user_data);

void
on_preferences_map_radius_entry_changed
                                        (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_preferences_map_increment_size_entry_activate
                                        (GtkEntry        *entry,
                                        gpointer         user_data);

void
on_preferences_map_increment_size_entry_changed
                                        (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_preferences_map_diff_increment_entry_activate
                                        (GtkEntry        *entry,
                                        gpointer         user_data);

void
on_preferences_map_diff_increment_entry_changed
                                        (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_preferences_map_sampling_entry_activate
                                        (GtkEntry        *entry,
                                        gpointer         user_data);

void
on_preferences_map_sampling_entry_changed
                                        (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_preferences_map_dynamic_sampling_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_map_dynamic_size_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_diff_map_colours_coot_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_diff_map_colours_o_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_map_colours_hscale_value_changed
                                        (GtkRange        *range,
                                        gpointer         user_data);

void
on_preferences_smooth_scroll_on_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_smooth_scroll_off_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_smooth_scroll_steps_entry_activate
                                        (GtkEntry        *entry,
                                        gpointer         user_data);

void
on_preferences_smooth_scroll_steps_entry_changed
                                        (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_preferences_smooth_scroll_limit_entry_activate
                                        (GtkEntry        *entry,
                                        gpointer         user_data);

void
on_preferences_smooth_scroll_limit_entry_changed
                                        (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_preferences_map_drag_on_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_map_drag_off_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_antialias_on_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_antialias_off_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_hid_spherical_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_hid_flat_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_filechooser_off_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_filechooser_on_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_file_overwrite_yes_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_file_overwrite_no_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_file_filter_on_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_file_filter_off_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_file_sort_by_date_on_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_file_sort_by_date_off_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_dialog_accept_docked_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_dialog_accept_detouched_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_dialog_accept_docked_show_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_dialog_accept_docked_hide_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_dialog_accept_on_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_dialog_accept_off_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_recentre_pdb_on_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_recentre_pdb_off_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_console_info_on_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_console_info_off_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_tips_on_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_tips_off_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_refinement_speed_molasses_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_refinement_speed_crock_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_refinement_speed_default_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_refinement_speed_own_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_refinement_speed_entry_activate
                                        (GtkEntry        *entry,
                                        gpointer         user_data);

void
on_preferences_refinement_speed_entry_changed
                                        (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_preferences_spin_speed_entry_activate
                                        (GtkEntry        *entry,
                                        gpointer         user_data);

void
on_preferences_spin_speed_entry_changed
                                        (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_preferences_font_size_small_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_font_size_medium_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_font_size_large_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_font_size_others_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

#if (GTK_MAJOR_VERSION > 1)
void
on_preferences_font_size_combobox_changed
                                        (GtkComboBox     *combobox,
                                        gpointer         user_data);

void
on_preferences_font_colour_default_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_font_colour_own_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_font_colorbutton_color_set
                                        (GtkColorButton  *colorbutton,
                                        gpointer         user_data);
#endif

void
on_preferences_font_colorbutton_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_preferences_pink_pointer_entry_activate
                                        (GtkEntry        *entry,
                                        gpointer         user_data);

void
on_preferences_pink_pointer_entry_changed
                                        (GtkEditable     *editable,
                                        gpointer         user_data);

#if (GTK_MAJOR_VERSION > 1)
void
on_preferences_bond_width_combobox_changed
                                        (GtkComboBox     *combobox,
                                        gpointer         user_data);
#endif

void
on_preferences_model_toolbar_show_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_model_toolbar_right_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_model_toolbar_left_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_model_toolbar_top_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_model_toolbar_bottom_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_model_toolbar_hide_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_model_toolbar_main_icons_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_model_toolbar_all_icons_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_model_toolbar_user_defined1_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_preferences_model_toolbar_style_icons_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_model_toolbar_style_both_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_model_toolbar_style_text_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_model_toolbar_show_icon_all_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_preferences_model_toolbar_show_icon_selection_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_preferences_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_preferences_reset_button_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_preferences_main_toolbar_show_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_main_toolbar_hide_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_main_toolbar_top_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_main_toolbar_bottom_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_main_toolbar_right_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_main_toolbar_left_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_main_toolbar_style_icons_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_main_toolbar_style_both_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_main_toolbar_style_text_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_preferences_destroy                 (GtkObject       *object,
                                        gpointer         user_data);


void
on_diff_map_peaks_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_generate_diff_map_peaks_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_generate_diff_map_peaks_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_difference_map_peaks1_activate      (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_superpose_reference_chain_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_superpose_moving_chain_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_draw_ncs_ghosts_yes_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_draw_ncs_ghosts_no_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_ncs_maps_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data);

void
on_ncs_maps_cancel_button_clicked      (GtkButton       *button,
                                        gpointer         user_data);

void
on_ncs_maps1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_rotamer_selection_dialog_destroy    (GtkObject       *object,
                                        gpointer         user_data);

void
on_pointer_distances_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_pointer_distances_ok_button_clicked (GtkButton       *button,
                                        gpointer         user_data);

void
on_pointer_distances1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_align_and_mutate_ok_button_clicked  (GtkButton       *button,
                                        gpointer         user_data);

void
on_align_and_mutate_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_ramachandran_plot_differences_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_ramachandran_plot_differences_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_ramachandran_differences_plot1_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_ramachandran_plot_differences_first_chain_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_ramachandran_plot_differences_second_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_check_waters1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_checked_waters_baddies_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_checked_waters_baddies_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_align___mutate1_activate            (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_about_other_button_clicked          (GtkButton       *button,
                                        gpointer         user_data);

void
on_align___mutate1_activate            (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_delete_item_keep_active_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_delete_item_dialog_destroy          (GtkObject       *object,
                                        gpointer         user_data);

void
on_save_state1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_geometry_torsion_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_main_window_statusbar_text_popped   (GtkStatusbar    *statusbar,
                                        guint            context_id,
                                        gchar           *text,
                                        gpointer         user_data);

void
on_main_window_statusbar_text_pushed   (GtkStatusbar    *statusbar,
                                        guint            context_id,
                                        gchar           *text,
                                        gpointer         user_data);

void
on_fit_loop1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_fit_loop_by_database_search1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_fit_loop_by_rama_search1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_fit_loop1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_fit_loop_ok_button_clicked          (GtkButton       *button,
                                        gpointer         user_data);

void
on_reset_view1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_base_chooser_A_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_base_chooser_C_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_base_chooser_G_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_base_chooser_T_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_base_chooser_U_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_base_chooser_cancel_button_clicked  (GtkButton       *button,
                                        gpointer         user_data);

void
on_nucleic_acid_base_chooser_dialog_destroy
                                        (GtkObject       *object,
                                        gpointer         user_data);

void
on_change_chain_ids2_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_change_chains_rechain_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_change_chain_cancel_button_clicked  (GtkButton       *button,
                                        gpointer         user_data);

void
on_delete_item_residue_range_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_on_line_docs_url1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_on_line_documentation_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_save_state_ok_button1_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_save_state_cancel_button1_clicked   (GtkButton       *button,
                                        gpointer         user_data);

void
on_change_chain_residue_range_no_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_change_chain_residue_range_yes_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_mutate_sequence_do_autofit_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_mutate_sequence_use_ramachandran_restraints_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_check_waters_b_factor_entry_active_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_check_waters_min_dist_entry_active_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_check_waters_max_dist_entry_active_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_check_waters_map_sigma_entry_active_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_residue_info_occ_apply_all_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_residue_info_occ_apply_to_altconf_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
					 gpointer         user_data);


void
on_residue_info_b_factor_apply_all_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_delete_item_water_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_modelling_cis_trans_button_clicked  (GtkButton       *button,
                                        gpointer         user_data);

void
on_other_modelling_tools_close_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_other_modelling_tools1_activate     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_cis_trans_conversion_toggle_button_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_other_model_tools_dialog_destroy    (GtkObject       *object,
                                        gpointer         user_data);

void
on_undo_last_navigation1_activate      (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_get_pdb_and_map_using_eds1_activate (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_fetch_pdb_and_map_using_pdbredo1_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_ligand_big_blob_dialog_destroy      (GtkObject       *object,
                                        gpointer         user_data);

void
on_do_180_degree_sidechain_flip_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_model_refine_dialog_do_180_degree_sidechain_flip_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_simple1_activate                    (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_povray1_activate                    (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_raster3d1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_screendump_image_ok_button_clicked  (GtkButton       *button,
                                        gpointer         user_data);

void
on_screendump_image_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_reverse_fragment_direction_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_place_helix_here_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_place_helix_here_button_clicked     (GtkButton       *button,
                                        gpointer         user_data);

void
on_other_tools_place_strand_here_button_clicked     (GtkButton       *button,
                                        gpointer         user_data);


void
on_diff_map_peaks_dialog_destroy       (GtkObject       *object,
                                        gpointer         user_data);

void
on_symmetry_controller_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_molecule_0_checkbutton_toggled      (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_display_sphere_radiobutton_molecule_0_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_display_all_radiobutton_molecule_0_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_display_CA_radiobutton_molecule_0_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_colour_symm_std_molecule_0_toggled  (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_colour_symm_by_symop_molecule_0_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_colour_symm_by_molecule_molecule_0_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_show_symmetry_molecule_control_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_show_symmetry_apply_button_clicked  (GtkButton       *button,
                                        gpointer         user_data);

void
on_ncs_controller_molecule_n_display_ncs_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_ncs_controller_molecule_n_display_chain_ich_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_ncs_controller_ncs_master_chain_ich_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_ncs_control_ok_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_lsq_plane_add_atom_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_lsq_plane_deviant_atom_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_lsq_plane_ok_button_clicked         (GtkButton       *button,
                                        gpointer         user_data);

void
on_lsq_plane_dialog_destroy            (GtkObject       *object,
                                        gpointer         user_data);

void
on_plane_distances1_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_lsq_plane_delete_last_atom_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_ncs_ghost_control1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_coords_colour_control_dialog_destroy
                                        (GtkObject       *object,
                                        gpointer         user_data);

void
on_coord_colour_control_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_bond_colours1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_accept_reject_refinement_dialog_destroy
                                        (GtkObject       *object,
                                        gpointer         user_data);

void
on_coot_doc_url_monolithic_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_coot_doc_url_sectioned_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

gboolean
on_coot_online_doc_search_entry_key_press_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data);

void
on_coot_online_doc_search_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_smiles1_activate                    (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_bond_parameters_rotate_colour_map_c_only_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_generic_display_objects1_activate   (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

gboolean
on_entry1_key_press_event              (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data);

void
on_map_radius_apply_button_clicked     (GtkButton       *button,
                                        gpointer         user_data);

void
on_refine_params_use_peptide_torsions_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_refine_params_use_helix_peptide_torsions_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_refine_params_use_beta_strand_peptide_torsions_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_refine_params_use_ramachandran_goodness_torsions_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_refine_params_use_planar_peptides_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_refine_params_use_trans_peptide_restraints_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_refine_params_use_peptide_omegas_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_probe_clashes1_activate             (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_validate1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_spin_view_on_off1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_rock_view_on_off1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_other_tools_RNA_button_clicked      (GtkButton       *button,
                                        gpointer         user_data);

void
on_other_tools_base_pair_toggle_button_toggled      (GtkToggleButton       *button,
                                        gpointer         user_data);

void
on_ideal_rna_ok_button_clicked         (GtkButton       *button,
                                        gpointer         user_data);

void
on_ideal_rna_cancel_button_clicked     (GtkButton       *button,
                                        gpointer         user_data);

void
on_unit_cell_yes_radiobutton_toggled   (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_unit_cell_no_radiobutton_toggled    (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_clear_atom_labels1_activate         (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_move_molecule_here1_activate        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_move_molecule_here_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_move_molecule_here_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_monomer_library_search_dialog_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_search_monomer_library1_activate    (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_monomer_library_search_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_monomer_search_entry_changed        (GtkEditable     *editable,
                                        gpointer         user_data);

gboolean
on_monomer_search_entry_key_press_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data);

void
on_least_squares_match_type_all_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_least_squares_match_type_main_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_least_squares_ok_button_clicked     (GtkButton       *button,
                                        gpointer         user_data);

void
on_least_squares_cancel_button_clicked (GtkButton       *button,
                                        gpointer         user_data);


void
on_least_squares_close_button_clicked  (GtkButton       *button,
                                        gpointer         user_data);

void
on_lsq_superpose1_activate             (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_goto_atom_window_state_changed      (GtkWidget       *widget,
                                        GtkStateType     state,
                                        gpointer         user_data);

gboolean
on_goto_atom_window_configure_event    (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data);

gboolean
on_model_refine_dialog_configure_event (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data);

gboolean
on_display_control_window_glade_configure_event
                                        (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data);

gboolean
on_delete_item_dialog_configure_event  (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data);

gboolean
on_rotate_translate_obj_dialog_configure_event
                                        (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data);

gboolean
on_accept_reject_refinement_dialog_configure_event
                                        (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data);

gboolean
on_dynarama_window_configure_event     (GtkWidget       *widget,
                                        GdkEventConfigure *event,
                                        gpointer         user_data);

void
on_coords_fileselection1_destroy       (GtkObject       *object,
                                        gpointer         user_data);

void
on_dataset_fileselection1_destroy      (GtkObject       *object,
                                        gpointer         user_data);

void
on_map_name_fileselection1_destroy     (GtkObject       *object,
                                        gpointer         user_data);

void
on_phs_coordinates_fileselection_destroy
                                        (GtkObject       *object,
                                        gpointer         user_data);

void
on_save_coords_fileselection1_destroy  (GtkObject       *object,
                                        gpointer         user_data);

void
on_save_symmetry_coords_fileselection_destroy
                                        (GtkObject       *object,
                                        gpointer         user_data);

void
on_save_state_fileselection_destroy    (GtkObject       *object,
                                        gpointer         user_data);

void
on_screendump_fileselection_destroy    (GtkObject       *object,
                                        gpointer         user_data);

gboolean
on_window1_delete_event                (GtkWidget       *widget,
                                        GdkEvent        *event,
                                        gpointer         user_data);

void
on_delete_item_sidechain_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);
void
on_delete_item_sidechain_range_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_delete_item_chain_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);
void
on_residue_type_chooser_stub_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_set_undo_molecule_button_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_gln_and_asn_b_factor_outiers1_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_gln_and_asn_b_factor_outliers1_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_stereo_dialog_side_by_side_stereo_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_sec_str_rest_no_rest_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_sec_str_rest_helix_rest_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_sec_str_rest_strand_rest_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_superpose_nonsense_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_refine_params_dialog_destroy        (GtkObject       *object,
                                        gpointer         user_data);

void
on_stereo_dialog_side_by_side_stereo_crosseyed_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_stereo_dialog_side_by_side_stereo_walleyed_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_update_go_to_atom_from_current_position_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_show_aniso_close_button_clicked     (GtkButton       *button,
                                        gpointer         user_data);

void
on_display_manager_button_clicked      (GtkButton       *button,
                                        gpointer         user_data);

#if (GTK_MAJOR_VERSION == 2)

void
on_display_manager_toolbutton_clicked  (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_symmetry_colorbutton_color_set      (GtkColorButton  *colorbutton,
                                        gpointer         user_data);

void
on_coords_filechooserdialog1_response  (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data);

void
on_coords_filechooserdialog1_destroy
                                        (GtkObject       *object,
                                        gpointer         user_data);

void
on_coords_filechooserdialog1_recentre_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_dataset_filechooserdialog1_response (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data);

void
on_dataset_filechooserdialog1_destroy  (GtkObject       *object,
                                        gpointer         user_data);

void
on_map_name_filechooserdialog1_response
                                        (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data);

void
on_map_name_filechooserdialog1_destroy (GtkObject       *object,
                                        gpointer         user_data);

void
on_map_filechooser_is_difference_map_button_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_phs_coordinates_filechooserdialog1_response
                                        (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data);

void
on_phs_coordinates_filechooserdialog1_destroy
                                        (GtkObject       *object,
                                        gpointer         user_data);

void
on_save_coords_filechooserdialog1_response
					(GtkDialog * dialog, 
					gint response_id, 
					gpointer user_data);

void
on_save_coords_filechooserdialog1_destroy
					(GtkObject * object, 
					gpointer user_data);

void
on_cif_dictionary_filechooserdialog1_response
					(GtkDialog * dialog, 
					gint response_id, 
					gpointer user_data);

void
on_cif_dictionary_filechooserdialog1_destroy
					(GtkObject * object, 
					gpointer user_data);

void
on_run_script_filechooserdialog1_response
					(GtkDialog * dialog, 
					gint response_id, 	
					gpointer user_data);

void
on_run_script_filechooserdialog1_destroy
					(GtkObject * object, 
					gpointer user_data);

#if ((GTK_MAJOR_VERSION == 2) && (GTK_MINOR_VERSION > 9))
GtkFileChooserConfirmation
on_save_coords_filechooserdialog1_confirm_overwrite
					(GtkFileChooser * filechooser, 
					gpointer user_data);

GtkFileChooserConfirmation
on_save_symmetry_coords_filechooserdialog1_confirm_overwrite
					(GtkFileChooser * filechooser, 
					gpointer user_data);

GtkFileChooserConfirmation
on_save_state_filechooserdialog1_confirm_overwrite
					(GtkFileChooser * filechooser, 
					gpointer user_data);

GtkFileChooserConfirmation
on_screendump_filechooserdialog1_confirm_overwrite
					(GtkFileChooser * filechooser, 
					gpointer user_data);
#endif	/* GTK_MINOR_VERSION */

void
on_save_symmetry_coords_filechooserdialog1_response
					(GtkDialog * dialog, 
					gint response_id, 
					gpointer user_data);

void
on_save_symmetry_coords_filechooserdialog1_destroy
					(GtkObject * object, 
					gpointer user_data);


void
on_save_state_filechooserdialog1_response
					(GtkDialog * dialog, 
					gint response_id, 
					gpointer user_data);

void
on_save_state_filechooserdialog1_destroy
					(GtkObject * object, 
					gpointer user_data);


void
on_screendump_filechooserdialog1_response
					(GtkDialog * dialog, 
					gint response_id, 
					gpointer user_data);

void
on_screendump_filechooserdialog1_destroy
					(GtkObject * object, 
					gpointer user_data);

#endif 


void
on_display_control_all_maps_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_display_control_all_models_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_single_map_properties_contour_level_apply_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_checked_waters_baddies_dialog_destroy
                                        (GtkObject       *object,
                                        gpointer         user_data);

void
on_model_toolbar_style_changed         (GtkToolbar      *toolbar,
                                        GtkToolbarStyle  style,
                                        gpointer         user_data);

gboolean
on_model_toolbar_button_press_event    (GtkWidget       *widget,
                                        GdkEventButton  *event,
                                        gpointer         user_data);

void
on_model_toolbar_refine_control_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_toolbar_select_map_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

#if (GTK_MAJOR_VERSION > 1)
void
on_model_toolbar_refine_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data);

void
on_model_toolbar_regularize_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data);


void
on_model_toolbar_fixed_atoms_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);


void
on_model_toolbar_rigid_body_fit_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data);

#ifdef GTK_TYPE_MENU_TOOL_BUTTON
void
on_model_toolbar_rot_trans_toolbutton_show_menu
                                        (GtkMenuToolButton *menutoolbutton,
                                        gpointer         user_data);
#endif

#ifdef GTK_TYPE_MENU_TOOL_BUTTON
void
on_model_toolbar_rot_trans_toolbutton_clicked
                                        (GtkMenuToolButton *menutoolbutton,
                                        gpointer         user_data);
#endif

void
on_model_toolbar_rot_trans_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data);

void
on_model_toolbar_auto_fit_rotamer_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data);

void
on_model_toolbar_rotamers_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data);

void
on_model_toolbar_edit_chi_angles_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data);

void
on_model_toolbar_torsion_general_toggletoolbutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data);

void
on_model_toolbar_flip_peptide_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data);

void
on_model_toolbar_sidechain_180_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data);

void
on_model_toolbar_edit_backbone_torsions_toggletoolbutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data);

void
on_model_toolbar_mutate_and_autofit_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data);

void
on_model_toolbar_simple_mutate_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data);

void
on_model_toolbar_add_terminal_residue_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data);

void
on_model_toolbar_add_alt_conf_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_model_toolbar_refmac_button_clicked (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_model_toolbar_display_manager_togglebutton_toggled
                                        (GtkToggleToolButton *toggletoolbutton,
                                        gpointer         user_data);

#endif /* GTK_MAJOR_VERSION */

void
on_model_toolbar_find_water_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_toolbar_add_atom_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_toolbar_clear_pending_picks_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_toolbar_delete_button_clicked (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_toolbar_undo_button_clicked   (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_toolbar_redo_button_clicked   (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_toolbar_icons1_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_model_toolbar_text1_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_model_toolbar_icons_and_text1_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_model_toolbar_main_icons_activate   (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_model_toolbar_all_icons_activate    (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_model_toolbar_setting1_activate     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_toolbar_display_manager_maps_all_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_toolbar_display_manager_molecules_all_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_scripting_python1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_scripting_scheme1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_scheme_window_close_button_clicked  (GtkButton       *button,
                                        gpointer         user_data);

void
on_aboutdialog_close                   (GtkDialog       *dialog,
                                        gpointer         user_data);
void
on_aboutdialog_response                (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data);
void
on_check_waters_check1_activate        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_check_waters_delete1_activate       (GtkMenuItem     *menuitem,
                                        gpointer         user_data);



void
on_model_refine_dialog_torsion_general_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);
void
on_accept_reject_reverse_button_clicked
                                        (GtkButton       *button,
					 gpointer         user_data);

void
on_accept_reject_reverse_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_geometry_dynamic_distance_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);
void
on_accession_code_ok_button_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_ncs_differences1_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_fix_atom_togglebutton_toggled       (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_unfix_atom_togglebutton_toggled     (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_clear_fixed_atoms_button_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_fixed_atom_close_button_clicked     (GtkButton       *button,
                                        gpointer         user_data);

void
on_fixed_atom_dialog_destroy           (GtkObject       *object,
                                        gpointer         user_data);

void
on_add_rep_add_rep_button_clicked      (GtkButton       *button,
                                        gpointer         user_data);

void
on_add_rep_cancel_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);


void
on_add_rep_add_rep_button_clicked      (GtkButton       *button,
                                        gpointer         user_data);

void
on_add_rep_cancel_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_additional_representation1_activate (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_display_additional_representations_close_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_all1_activate                       (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_all2_activate                       (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_residue_editor_select_monomer_type_ok_button_clicked (GtkButton       *button,
						    gpointer         user_data);

void
on_residue_editor_select_monomer_type_cancel_button_clicked (GtkButton       *button,
						    gpointer         user_data);
void
on_restraints1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data);


void
on_restraint_editor_add_restraint_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_restraints_editor_close_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_restraints_editor_apply_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);
void
on_restraints_editor_save_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_restraint_editor_delete_restraint_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

#if ((GTK_MAJOR_VERSION == 2) && (GTK_MINOR_VERSION > 9))
GtkFileChooserConfirmation
on_save_restraint_chooserdialog_confirm_overwrite
                                        (GtkFileChooser  *filechooser,
                                        gpointer         user_data);
#endif

void
on_save_restraint_chooserdialog_response
                                        (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data);

void
on_save_restraint_chooserdialog_close  (GtkDialog       *dialog,
                                        gpointer         user_data);

void
on_run_refmac_nolabels_help_dialog_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_fix_atoms_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_model_refine_dialog_fast_sss_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_other_tools_build_na_button_clicked (GtkButton       *button,
                                        gpointer         user_data);

void
on_build_na_dialog_radius_entry_activate
                                        (GtkEntry        *entry,
                                        gpointer         user_data);

void
on_build_na_dialog_cancelbutton_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_build_na_dialog_okbutton_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_fast_sss_dialog_citation_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

#if (GTK_MAJOR_VERSION > 1)
void
on_coot_references_coot_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_coot_references_wincoot_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_coot_references_refmac_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_coot_references_ssm_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_coot_references_mmdb_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_coot_references_clipper_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_coot_references_buccaneer_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_coot_references_molprobity_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_coot_references_calpha_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_coot_references_xligand_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_coot_references_eds_toolbutton_clicked (GtkToolButton   *toolbutton,
					   gpointer         user_data);

void
on_coot_references_others_toolbutton_clicked
                                        (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_toolbar_multi_refine_continue_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_toolbar_multi_refine_stop_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_toolbar_multi_refine_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void on_coords_toolbutton_clicked (GtkToolButton   *toolbutton,
				   gpointer         user_data);

void
on_go_to_atom_toolbutton_clicked       (GtkToolButton   *toolbutton,
                                        gpointer         user_data);


void
on_go_to_ligand_toolbutton_clicked     (GtkToolButton   *toolbutton,
                                        gpointer         user_data);


#endif /* GTK_MAJOR_VERSION */

void
on_coot_references_closebutton_clicked (GtkButton       *button,
                                        gpointer         user_data);


void
on_refine_params_use_ramachandran_goodness_torsions_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_edit_chi_angles_add_hydrogen_torsions_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_find_ligands_search_all_radiobutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_find_ligands_search_here_radiobutton_toggled
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_symmetry_controller_dialog_destroy  (GtkObject       *object,
                                        gpointer         user_data);

void
on_map_sharpening1_activate            (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_map_sharpening_ok_button_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_map_sharpening_optimize_button_clicked    (GtkButton       *button,
                                        gpointer         user_data);

void
on_map_sharpening_reset_button_clicked (GtkButton       *button,
                                        gpointer         user_data);

void
on_map_sharpening_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);


void
on_baton_build_params_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_baton_build_params_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);
void
on_baton_build_set_params_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);
void
on_edit_chi_angles_reverse_fragment_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_move_molecule_here_big_molecules_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_environment_distances_h_bonds_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_environment_distances_bumps_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_check_water_by_difference_map_optionmenu_changed
                                        (GtkOptionMenu   *optionmenu,
                                        gpointer         user_data);

void
on_pisa_interfces_close_button_clicked (GtkButton       *button,
                                        gpointer         user_data);

void
on_displayed_map_style_as_lines_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_displayed_map_style_as_cut_glass_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_displayed_map_style_as_transparent_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_map_opacity_hscale_value_changed    (GtkRange        *range,
                                        gpointer         user_data);


void
on_refine_params_weight_matrix_entry_changed
                                        (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_remarks_browser1_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data);
void
on_remarks_browser_molecule_chooser_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_remarks_browser_molecule_chooser_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_fix_nomenclature_errors_ok_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_fix_nomenclature_errors_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);
void
on_ligand_builder1_activate            (GtkMenuItem     *menuitem,
                                        gpointer         user_data);
void
on_multi_residue_torsion_cancel_button_clicked
                                        (GtkButton       *button,
					 gpointer         user_data);

void
on_multi_residue_torsion_OK_button_clicked
                                        (GtkButton       *button,
					 gpointer         user_data);

void
on_multi_residue_torsion_reverse_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
					 gpointer         user_data);


void
on_multi_residue_torsion_pick_apply_button_clicked
                                        (GtkButton       *button,
					 gpointer         user_data);


void
on_multi_residue_torsion_pick_cancel_button_activate
                                        (GtkButton       *button,
					 gpointer         user_data); /* needed? No. */

void
on_multi_residue_torsion_pick_cancel_button_clicked
                                        (GtkButton       *button,
					 gpointer         user_data);


void
on_multi_residue_torsion_start_button_clicked
                                        (GtkButton       *button,
					 gpointer         user_data);



void
on_keyboard_go_to_residue_entry_changed   (GtkEditable     *editable,
                                        gpointer         user_data);


gboolean
on_keyboard_go_to_residue_entry_key_press_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
					 gpointer         user_data);
void
on_mogul_geometry_dialog_close_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void 
on_ligand_check_okbutton_clicked(GtkButton       *button,
				 gpointer         user_data);

void
on_generic_objects_dialog_closebutton_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_generic_objects_dialog_close        (GtkDialog       *dialog,
                                        gpointer         user_data);

void
on_generic_objects_dialog_destroy      (GtkObject       *object,
                                        gpointer         user_data);

void
on_generic_objects_display_all_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
					 gpointer         user_data);
void
on_generic_objects_close_all_button_clicked
                                        (GtkButton       *button,
					 gpointer         user_data);
void
on_export_map1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_export_map_fragment1_activate       (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_export_map_dialog_ok_button_clicked (GtkButton       *button,
                                        gpointer         user_data);

void
on_export_map_dialog_cancel_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);
/* void */
/* on_export_map_filechooserdialog_cancel_button_clicked */
/*                                         (GtkButton       *button, */
/*                                         gpointer         user_data); */

/* void */
/* on_export_map_filechooserdialog_save_button_clicked */
/*                                         (GtkButton       *button, */
/*                                         gpointer         user_data); */
void
on_export_map_filechooserdialog_response
                                        (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data);

void
on_cfc_dialog_response                 (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data);


void
on_dynarama_outliers_only_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
					 gpointer         user_data);

void
on_find_ligand_real_space_refine_solutions_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
					 gpointer         user_data);

void
on_edit_copy_molecule1_activate        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_edit_copy_fragment1_activate        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_edit_replace_residue1_activate      (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_edit_replace_fragment1_activate     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_edit_renumber_residues1_activate    (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_edit_change_chain_ids1_activate     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_edit_merge_molecules1_activate      (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_weight_maxtrix_estimate_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_mutate_molecule_resno_1_entry_changed
                                        (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_mutate_molecule_resno_2_entry_changed
                                        (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_mutate_molecule_sequence_text_insert_at_cursor
                                        (GtkTextView     *textview,
                                        gchar           *string,
                                        gpointer         user_data);

gboolean
on_mutate_molecule_sequence_text_key_release_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data);

gboolean
on_mutate_molecule_sequence_text_button_release_event
                                        (GtkWidget       *widget,
                                        GdkEventButton  *event,
                                        gpointer         user_data);

void
on_display_control_last_model_only_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_display_control_align_labels_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_delete_item_sidechain_range_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_symmetry_colorbutton_color_set      (GtkColorButton  *colorbutton,
                                        gpointer         user_data);

void
on_symmetry_colorbutton_color_set      (GtkColorButton  *colorbutton,
                                        gpointer         user_data);

void
on_curlew_install_button_clicked       (GtkButton       *button,
                                        gpointer         user_data);

void
on_curlew_dialog_closebutton_clicked   (GtkButton       *button,
                                        gpointer         user_data);

void
on_curlew_dialog_close                 (GtkDialog       *dialog,
                                        gpointer         user_data);

void
on_curlew_dialog_response              (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data);

void
on_modelling_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);


void
on_edit_settings_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_align_and_mutate1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_calculate_all_molecule_activate     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_calculate_dock_sequence_activate    (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_calculate_map_tools_activate        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_calculate_modules_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_calculate_ncs_tools_activate        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_calculate_pisa_activate             (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_draw_representation_tools_activate  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_calculate_views_activate            (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_calculate_load_tutorial_model_and_data1_activate
                                        (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_refine_params_geman_mcclure_alpha_combobox_changed
                                        (GtkComboBox     *combobox,
                                        gpointer         user_data);

void
on_refine_params_lennard_jones_epsilon_combobox_changed
                                        (GtkComboBox     *combobox,
                                        gpointer         user_data);

void
on_refine_params_rama_restraints_weight_combobox_changed
                                        (GtkComboBox     *combobox,
                                        gpointer         user_data);

void
on_refine_params_more_control_togglebutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_accept_reject_flip_this_peptide_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_accept_reject_flip_next_peptide_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_accept_reject_crankshaft_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_accept_reject_backrub_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_accept_reject_backrub_rotamer_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_symmetry_always_on_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_curlew1_activate              (GtkMenuItem     *menuitem,
                                  gpointer         user_data);

void
on_show_symmetry_no_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_show_symmetry_yes_radiobutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

gboolean
on_symmetry_radius_entry_key_release_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data);

void
on_hscale_symmetry_colour_value_changed
                                        (GtkRange        *range,
                                        gpointer         user_data);

void
on_show_symmetry_expanded_labels_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_map_radius_em_entry_changed         (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_map_radius_em_button_clicked        (GtkButton       *button,
                                        gpointer         user_data);

gboolean
on_map_radius_em_entry_key_press_event (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data);

void
on_simple_refmac_dialog_response       (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data);

void
on_simple_refmac_dialog_close          (GtkDialog       *dialog,
                                        gpointer         user_data);

void
on_simple_refmac_mtz_file_button_clicked
                                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_simple_refmac_filechooserdialog_response
                                        (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data);

void
on_label_neighbours1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

gboolean
on_residue_type_chooser_entry_key_press_event
                                        (GtkWidget       *widget,
                                        GdkEventKey     *event,
                                        gpointer         user_data);

void
on_map_properties_dialog_specularity_state_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_map_properties_dialog_fresnel_state_checkbutton_toggled
                                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_map_properties_dialog_specularity_strength_entry_activate
                                        (GtkEntry        *entry,
                                        gpointer         user_data);

void
on_map_properties_dialog_specularity_shininess_entry_activate
                                        (GtkEntry        *entry,
                                        gpointer         user_data);

void
on_map_properties_dialog_fresnel_bias_entry_activate
                                        (GtkEntry        *entry,
                                        gpointer         user_data);

void
on_map_properties_dialog_fresnel_scale_entry_activate
                                        (GtkEntry        *entry,
                                        gpointer         user_data);

void
on_map_properties_dialog_fresnel_power_entry_activate
                                        (GtkEntry        *entry,
                                        gpointer         user_data);
