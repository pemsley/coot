#############
# some re-definitions from coot python functions
############
# sorted by header files:

import coot

# c-interface.h:
chain_id               = coot.chain_id_py
remarks                = coot.remarks_py
residue_centre         = coot.residue_centre_py
coot_sys_build_type    = coot.coot_sys_build_type_py
accept_moving_atoms    = coot.accept_moving_atoms_py
molecule_name_stub     = coot.molecule_name_stub_py
run_clear_backups      = coot.run_clear_backups_py
het_group_residues     = coot.het_group_residues_py
replace_residues_from_mol = coot.replace_residues_from_mol_py
test_internal          = coot.test_internal_py
test_internal_single   = coot.test_internal_single_py
select_atom_under_pointer = coot.select_atom_under_pointer_py
get_map_colour         = coot.get_map_colour_py
average_map            = coot.average_map_py
refmac_parameters      = coot.refmac_parameters_py
map_parameters         = coot.map_parameters_py
cell                   = coot.cell_py
map_cell               = coot.cell_py
get_refmac_sad_atom_info = coot.get_refmac_sad_atom_info_py
origin_pre_shift       = coot.origin_pre_shift_py
save_coords_name_suggestion = coot.save_coords_name_suggestion_py
save_state_file        = coot.save_state_file_py
save_state_file_name   = coot.save_state_file_name_py
run_state_file         = coot.run_state_file_py
go_to_ligand           = coot.go_to_ligand_py
centre_of_mass_string  = coot.centre_of_mass_string_py
set_atom_attributes    = coot.set_atom_attributes_py
test_function          = coot.test_function_py
glyco_tree             = coot.glyco_tree_py
glyco_tree_residues    = coot.glyco_tree_residues_py
glyco_tree_residue_id  = coot.glyco_tree_residue_id_py
glyco_tree_compare_trees = coot.glyco_tree_compare_trees_py
glyco_tree_matched_residue_pairs = coot.glyco_tree_matched_residue_pairs_py
display_maps           = coot.display_maps_py
space_group            = coot.space_group_py
symmetry_operators     = coot.symmetry_operators_py
symmetry_operators_to_xHM = coot.symmetry_operators_to_xHM_py
merge_molecules        = coot.merge_molecules_py
set_merge_molecules_ligand_spec = coot.set_merge_molecules_ligand_spec_py
alignment_results      = coot.alignment_results_py
nearest_residue_by_sequence = coot.nearest_residue_by_sequence_py
change_chain_id_with_result        = coot.change_chain_id_with_result_py
probe_available_p      = coot.probe_available_p_py
# to be renamed later
matching_compound_names_from_dictionary = coot.matching_compound_names_from_dictionary_py
# keep for compatibility reasons
comp_id_to_name           = coot.comp_id_to_name_py
refine_residues        = coot.refine_residues_py
refine_residues_with_modes_with_alt_conf = coot.refine_residues_with_modes_with_alt_conf_py
refine_residues_with_alt_conf        = coot.refine_residues_with_alt_conf_py
regularize_residues    = coot.regularize_residues_py
regularize_residues_with_alt_conf  = coot.regularize_residues_with_alt_conf_py
refine_zone_with_score = coot.refine_zone_with_score_py
regularize_zone_with_score = coot.regularize_zone_with_score_py

chiral_volume_errors = coot.chiral_volume_errors_py
add_extra_bond_restraints = coot.add_extra_bond_restraints_py
delete_extra_restraints_for_residue_spec = coot.delete_extra_restraints_for_residue_spec_py
generate_local_self_restraints_by_residues = coot.generate_local_self_restraints_by_residues_py
delete_extra_restraint = coot.delete_extra_restraint_py
list_extra_restraints  = coot.list_extra_restraints_py
add_atom_geometry_distance = coot.add_atom_geometry_distance_py
get_pointer_position_frac = coot.get_pointer_position_frac_py
non_standard_residue_names         = coot.non_standard_residue_names_py
chain_id_for_shelxl_residue_number = coot.chain_id_for_shelxl_residue_number_py
gln_asn_b_factor_outliers          = coot.gln_asn_b_factor_outliers_py
map_peaks              = coot.map_peaks_py
map_peaks_near_point   = coot.map_peaks_near_point_py
map_peaks_near_point_from_list = coot.map_peaks_near_point_from_list_py
map_peaks_around_molecule = coot.map_peaks_around_molecule_py
screen_vectors         = coot.screen_vectors_py
get_torsion            = coot.get_torsion_py
set_torsion            = coot.set_torsion_py
multi_residue_torsion  = coot.multi_residue_torsion_py
view_name              = coot.view_name_py
view_description       = coot.view_description_py
go_to_view             = coot.go_to_view_py
movie_file_name_prefix = coot.movie_file_name_prefix_py
execute_ligand_search  = coot.execute_ligand_search_py
overlap_ligands        = coot.overlap_ligands_py
analyse_ligand_differences = coot.analyse_ligand_differences_py
# compare_ligand_atom_types = coot.compare_ligand_atom_types_py # comment out hack to load this smoothly.
additional_representation_info = coot.additional_representation_info_py
rigid_body_refine_by_residue_ranges = coot.rigid_body_refine_by_residue_ranges_py
find_terminal_residue_type = coot.find_terminal_residue_type_py
delete_residues        = coot.delete_residues_py
cis_peptides           = coot.cis_peptides_py
twisted_trans_peptides = coot.twisted_trans_peptides_py
get_rotamer_name       = coot.get_rotamer_name_py
missing_atom_info      = coot.missing_atom_info_py
rotamer_graphs         = coot.rotamer_graphs_py
add_alt_conf           = coot.add_alt_conf_py
highly_coordinated_waters = coot.highly_coordinated_waters_py
metal_coordination     = coot.metal_coordination_py
add_lsq_atom_pair      = coot.add_lsq_atom_pair_py
apply_lsq_matches      = coot.apply_lsq_matches_py
get_lsq_matrix         = coot.get_lsq_matrix_py
make_image_raster3d    = coot.make_image_raster3d_py
make_image_povray      = coot.make_image_povray_py
raster_screen_shot     = coot.raster_screen_shot_py
ncs_master_chains      = coot.ncs_master_chains_py
copy_residue_range_from_ncs_master_to_chains = coot.copy_residue_range_from_ncs_master_to_chains_py
copy_from_ncs_master_to_chains = coot.copy_from_ncs_master_to_chains_py
ncs_chain_differences  = coot.ncs_chain_differences_py
ncs_chain_ids          = coot.ncs_chain_ids_py
ncs_ghosts             = coot.ncs_ghosts_py
new_molecule_by_residue_specs = coot.new_molecule_by_residue_specs_py
pucker_info            = coot.pucker_info_py
sequence_info          = coot.sequence_info_py
alignment_mismatches   = coot.alignment_mismatches_py
do_clipped_surface     = coot.do_clipped_surface_py
user_mods              = coot.user_mods_py
drag_intermediate_atom = coot.drag_intermediate_atom_py
mark_atom_as_fixed     = coot.mark_atom_as_fixed_py
mark_multiple_atoms_as_fixed = coot.mark_multiple_atoms_as_fixed_py
ccp4i_projects         = coot.ccp4i_projects_py
add_dipole             = coot.add_dipole_py
add_dipole_for_residues = coot.add_dipole_for_residues_py
get_pkgdatadir         = coot.get_pkgdatadir_py
pkgdatadir             = coot.get_pkgdatadir_py
handle_pisa_interfaces = coot.handle_pisa_interfaces_py
add_pisa_interface_bond = coot.add_pisa_interface_bond_py
matching_compound_names_from_sbase = coot.matching_compound_names_from_sbase_py
add_linked_residue     = coot.add_linked_residue_py
all_molecule_rotamer_score = coot.all_molecule_rotamer_score_py
all_molecule_ramachandran_score = coot.all_molecule_ramachandran_score_py
all_molecule_ramachandran_region = coot.all_molecule_ramachandran_region_py
user_defined_click     = coot.user_defined_click_py

# graphics_info.h:
#run_post_manipulation_hook = coot.run_post_manipulation_hook_py

# c-interface-python.hh:
make_atom_spec         = coot.make_atom_spec_py

# cc-interface.hh:
goto_next_atom_maybe   = coot.goto_next_atom_maybe_py
goto_prev_atom_maybe   = coot.goto_prev_atom_maybe_py
set_go_to_atom_from_res_spec = coot.set_go_to_atom_from_res_spec_py
set_go_to_atom_from_atom_spec = coot.set_go_to_atom_from_atom_spec_py
active_atom_spec       = coot.active_atom_spec_py
get_symmetry           = coot.get_symmetry_py
map_colour_components  = coot.map_colour_components_py
multi_residue_torsion_fit  = coot.multi_residue_torsion_fit_py
dictionaries_read      = coot.dictionaries_read_py
cif_file_for_comp_id   = coot.cif_file_for_comp_id_py
dictionary_entries     = coot.dictionary_entries_py
SMILES_for_comp_id     = coot.SMILES_for_comp_id_py
monomer_restraints     = coot.monomer_restraints_py
monomer_restraints_for_molecule = coot.monomer_restraints_for_molecule_py
set_monomer_restraints = coot.set_monomer_restraints_py
list_nomenclature_errors = coot.list_nomenclature_errors_py
residue_spec_make_triple = coot.residue_spec_make_triple_py
get_residue_specs_in_mol = coot.get_residue_specs_in_mol_py
get_residue_by_type    = coot.get_residue_by_type_py
atom_info_string       = coot.atom_info_string_py
residue_info           = coot.residue_info_py
residue_name           = coot.residue_name_py
all_residues_with_serial_numbers = coot.all_residues_with_serial_numbers_py
residue_centre_from_spec = coot.residue_centre_from_spec_py
chain_fragments        = coot.chain_fragments_py
set_b_factor_residues  = coot.set_b_factor_residues_py
clear_and_update_molecule = coot.clear_and_update_molecule_py
add_molecule           = coot.add_molecule_py
active_residue         = coot.active_residue_py
closest_atom_simple    = coot.closest_atom_simple_py
closest_atom           = coot.closest_atom_py
closest_atom_raw       = coot.closest_atom_raw_py
residues_near_residue  = coot.residues_near_residue_py
residues_near_residues = coot.residues_near_residues_py
residues_near_position = coot.residues_near_position_py
get_environment_distances_representation = coot.get_environment_distances_representation_py
all_residues_with_serial_numbers = coot.all_residues_with_serial_numbers_py
refine_zone_with_full_residue_spec = coot.refine_zone_with_full_residue_spec_py
morph_fit_residues     = coot.morph_fit_residues_py
find_blobs             = coot.find_blobs_py
water_chain_from_shelx_ins = coot.water_chain_from_shelx_ins_py
water_chain            = coot.water_chain_py
make_variance_map      = coot.make_variance_map_py
spin_search            = coot.spin_search_py
spin_N                 = coot.spin_N_py
ligand_search_make_conformers = coot.ligand_search_make_conformers_py
cootaneer              = coot.cootaneer_py
generic_string_vector_to_list_internal = coot.generic_string_vector_to_list_internal_py
generic_list_to_string_vector_internal = coot.generic_list_to_string_vector_internal_py
generic_int_vector_to_list_internal = coot.generic_int_vector_to_list_internal_py
inverse_rtop           = coot.inverse_rtop_py
# curl and hence coot_get_url_as_string are conditionally compiled.
try:
    coot_get_url_as_string = coot.coot_get_url_as_string_py
    curl_progress_info     = coot.curl_progress_info_py
except:
    pass # or print a message
score_rotamers         = coot.score_rotamers_py
protein_db_loops       = coot.protein_db_loops_py
make_link              = coot.make_link_py
link_info              = coot.link_info_py
map_to_model_correlation = coot.map_to_model_correlation_py
map_to_model_correlation_stats = coot.map_to_model_correlation_stats_py
map_to_model_correlation_per_residue = coot.map_to_model_correlation_per_residue_py
qq_plot_map_and_model  = coot.qq_plot_map_and_model_py
density_score_residue  = coot.density_score_residue_py
map_mean               = coot.map_mean_py
map_sigma              = coot.map_sigma_py
map_statistics         = coot.map_statistics_py
register_interesting_positions_list = coot.register_interesting_positions_list_py
molecule_atom_overlaps = coot.molecule_atom_overlaps_py
align_to_closest_chain = coot.align_to_closest_chain_py
key_sym_code           = coot.key_sym_code_py

# c-interface-ligands-swig.hh:
# new_molecule_sans_biggest_ligand = coot.new_molecule_sans_biggest_ligand_py # comment out hack to load this smoothly.
# gui_ligand_metrics     = coot.gui_ligand_metrics_py # comment out hack to load this smoothly.
residues_torsions_match = coot.residues_torsions_match_py
get_intermediate_atoms_distortions = coot.get_intermediate_atoms_distortions_py
ligand_atom_overlaps   = coot.ligand_atom_overlaps_py
residues_torsions_match = coot.residues_torsions_match_py
kolmogorov_smirnov     = coot.kolmogorov_smirnov_py
kolmogorov_smirnov_vs_normal = coot.kolmogorov_smirnov_vs_normal_py
kullback_liebler       = coot.kullback_liebler_py
linked_residues        = coot.linked_residues_py
add_target_position_restraint_for_intermediate_atom = coot.add_target_position_restraint_for_intermediate_atom_py
add_target_position_restraints_for_intermediate_atoms = coot.add_target_position_restraints_for_intermediate_atoms_py

coot_contact_dots_for_ligand = coot.coot_contact_dots_for_ligand_py
switch_HIS_protonation = coot.switch_HIS_protonation_py

# probe-clash-score.hh
probe_clash_score      = coot.probe_clash_score_py
# probe_clash_score_as   = coot.probe_clash_score_as_py # comment out hack to load this smoothly.

# c-interface-generic-objects.h
generic_object_name    = coot.generic_object_name_py

# and some acronyms
de_chainsaw                    = coot.fill_partial_residues

# fix typo of set_find_hydrogen_torsions (for backwards compatibility?!)
set_find_hydrogen_torsion = coot.set_find_hydrogen_torsions

# redefine in case someone has old function
toggle_idle_ligand_interactions = coot.toggle_flev_idle_ligand_interactions
