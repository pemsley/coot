#############
# some re-definitions from coot python functions
############
# sorted by header files:

# c-interface.h:
chain_id               = chain_id_py
remarks                = remarks_py
residue_centre         = residue_centre_py
coot_sys_build_type    = coot_sys_build_type_py
accept_moving_atoms    = accept_moving_atoms_py
molecule_name_stub     = molecule_name_stub_py
run_clear_backups      = run_clear_backups_py
het_group_residues     = het_group_residues_py
replace_residues_from_mol = replace_residues_from_mol_py
test_internal          = test_internal_py
test_internal_single   = test_internal_single_py
select_atom_under_pointer = select_atom_under_pointer_py
get_map_colour         = get_map_colour_py
average_map            = average_map_py
refmac_parameters      = refmac_parameters_py
map_parameters         = map_parameters_py
cell                   = cell_py
map_cell               = cell_py
get_refmac_sad_atom_info = get_refmac_sad_atom_info_py
origin_pre_shift       = origin_pre_shift_py
save_coords_name_suggestion = save_coords_name_suggestion_py
save_state_file        = save_state_file_py
save_state_file_name   = save_state_file_name_py
run_state_file         = run_state_file_py
go_to_ligand           = go_to_ligand_py
centre_of_mass_string  = centre_of_mass_string_py
set_atom_attributes    = set_atom_attributes_py
test_function          = test_function_py
glyco_tree             = glyco_tree_py
glyco_tree_residues    = glyco_tree_residues_py
glyco_tree_residue_id  = glyco_tree_residue_id_py
glyco_tree_compare_trees = glyco_tree_compare_trees_py
glyco_tree_matched_residue_pairs = glyco_tree_matched_residue_pairs_py
display_maps           = display_maps_py
space_group            = space_group_py
symmetry_operators     = symmetry_operators_py
symmetry_operators_to_xHM = symmetry_operators_to_xHM_py
merge_molecules        = merge_molecules_py
set_merge_molecules_ligand_spec = set_merge_molecules_ligand_spec_py
alignment_results      = alignment_results_py
nearest_residue_by_sequence = nearest_residue_by_sequence_py
change_chain_id_with_result        = change_chain_id_with_result_py
probe_available_p      = probe_available_p_py
# to be renamed later
matching_compound_names_from_dictionary = matching_compound_names_from_dictionary_py
# keep for compatibility reasons
comp_id2name           = comp_id_to_name_py
comp_id_to_name           = comp_id_to_name_py
refine_residues        = refine_residues_py
refine_residues_with_modes_with_alt_conf = refine_residues_with_modes_with_alt_conf_py
refine_residues_with_alt_conf        = refine_residues_with_alt_conf_py
regularize_residues    = regularize_residues_py
regularize_residues_with_alt_conf  = regularize_residues_with_alt_conf_py
refine_zone_with_score = refine_zone_with_score_py
regularize_zone_with_score = regularize_zone_with_score_py

chiral_volume_errors = chiral_volume_errors_py
add_extra_bond_restraints = add_extra_bond_restraints_py
delete_extra_restraints_for_residue_spec = delete_extra_restraints_for_residue_spec_py
generate_local_self_restraints_by_residues = generate_local_self_restraints_by_residues_py
delete_extra_restraint = delete_extra_restraint_py
list_extra_restraints  = list_extra_restraints_py
add_atom_geometry_distance = add_atom_geometry_distance_py
get_pointer_position_frac = get_pointer_position_frac_py
non_standard_residue_names         = non_standard_residue_names_py
chain_id_for_shelxl_residue_number = chain_id_for_shelxl_residue_number_py
gln_asn_b_factor_outliers          = gln_asn_b_factor_outliers_py
map_peaks              = map_peaks_py
map_peaks_near_point   = map_peaks_near_point_py
map_peaks_near_point_from_list = map_peaks_near_point_from_list_py
map_peaks_around_molecule = map_peaks_around_molecule_py
screen_vectors         = screen_vectors_py
get_torsion            = get_torsion_py
set_torsion            = set_torsion_py
multi_residue_torsion  = multi_residue_torsion_py
view_name              = view_name_py
view_description       = view_description_py
go_to_view             = go_to_view_py
movie_file_name_prefix = movie_file_name_prefix_py
execute_ligand_search  = execute_ligand_search_py
overlap_ligands        = overlap_ligands_py
analyse_ligand_differences = analyse_ligand_differences_py
# compare_ligand_atom_types = compare_ligand_atom_types_py # comment out hack to load this smoothly.
additional_representation_info = additional_representation_info_py
pepflip_using_difference_map = pepflip_using_difference_map_py
rigid_body_refine_by_residue_ranges = rigid_body_refine_by_residue_ranges_py
find_terminal_residue_type = find_terminal_residue_type_py
delete_residues        = delete_residues_py
cis_peptides           = cis_peptides_py
twisted_trans_peptides = twisted_trans_peptides_py
get_rotamer_name       = get_rotamer_name_py
missing_atom_info      = missing_atom_info_py
rotamer_graphs         = rotamer_graphs_py
add_alt_conf           = add_alt_conf_py
highly_coordinated_waters = highly_coordinated_waters_py
metal_coordination     = metal_coordination_py
add_lsq_atom_pair      = add_lsq_atom_pair_py
apply_lsq_matches      = apply_lsq_matches_py
get_lsq_matrix         = get_lsq_matrix_py
make_image_raster3d    = make_image_raster3d_py
make_image_povray      = make_image_povray_py
raster_screen_shot     = raster_screen_shot_py
ncs_master_chains      = ncs_master_chains_py
copy_residue_range_from_ncs_master_to_chains = copy_residue_range_from_ncs_master_to_chains_py
copy_from_ncs_master_to_chains = copy_from_ncs_master_to_chains_py
ncs_chain_differences  = ncs_chain_differences_py
ncs_chain_ids          = ncs_chain_ids_py
ncs_ghosts             = ncs_ghosts_py
new_molecule_by_residue_specs = new_molecule_by_residue_specs_py
pucker_info            = pucker_info_py
sequence_info          = sequence_info_py
alignment_mismatches   = alignment_mismatches_py
do_clipped_surface     = do_clipped_surface_py
user_mods              = user_mods_py
drag_intermediate_atom = drag_intermediate_atom_py
mark_atom_as_fixed     = mark_atom_as_fixed_py
mark_multiple_atoms_as_fixed = mark_multiple_atoms_as_fixed_py
ccp4i_projects         = ccp4i_projects_py
add_dipole             = add_dipole_py
add_dipole_for_residues = add_dipole_for_residues_py
get_pkgdatadir         = get_pkgdatadir_py
pkgdatadir             = get_pkgdatadir_py
handle_pisa_interfaces = handle_pisa_interfaces_py
add_pisa_interface_bond = add_pisa_interface_bond_py
matching_compound_names_from_sbase = matching_compound_names_from_sbase_py
add_linked_residue     = add_linked_residue_py
all_molecule_rotamer_score = all_molecule_rotamer_score_py
all_molecule_ramachandran_score = all_molecule_ramachandran_score_py
all_molecule_ramachandran_region = all_molecule_ramachandran_region_py
add_cablam_markup      = add_cablam_markup_py
user_defined_click     = user_defined_click_py

# graphics_info.h:
#run_post_manipulation_hook = run_post_manipulation_hook_py

# c-interface-python.hh:
make_atom_spec         = make_atom_spec_py

# cc-interface.hh:
goto_next_atom_maybe   = goto_next_atom_maybe_py
goto_prev_atom_maybe   = goto_prev_atom_maybe_py
set_go_to_atom_from_res_spec = set_go_to_atom_from_res_spec_py
set_go_to_atom_from_atom_spec = set_go_to_atom_from_atom_spec_py
active_atom_spec       = active_atom_spec_py
get_symmetry           = get_symmetry_py
map_colour_components  = map_colour_components_py
multi_residue_torsion_fit  = multi_residue_torsion_fit_py
dictionaries_read      = dictionaries_read_py
cif_file_for_comp_id   = cif_file_for_comp_id_py
dictionary_entries     = dictionary_entries_py
SMILES_for_comp_id     = SMILES_for_comp_id_py
monomer_restraints     = monomer_restraints_py
monomer_restraints_for_molecule = monomer_restraints_for_molecule_py
set_monomer_restraints = set_monomer_restraints_py
list_nomenclature_errors = list_nomenclature_errors_py
residue_spec_make_triple = residue_spec_make_triple_py
get_residue_specs_in_mol = get_residue_specs_in_mol_py
get_residue_by_type    = get_residue_by_type_py
atom_info_string       = atom_info_string_py
residue_info           = residue_info_py
residue_name           = residue_name_py
all_residues_with_serial_numbers = all_residues_with_serial_numbers_py
residue_centre_from_spec = residue_centre_from_spec_py
chain_fragments        = chain_fragments_py
set_b_factor_residues  = set_b_factor_residues_py
clear_and_update_molecule = clear_and_update_molecule_py
add_molecule           = add_molecule_py
active_residue         = active_residue_py
closest_atom_simple    = closest_atom_simple_py
closest_atom           = closest_atom_py
closest_atom_raw       = closest_atom_raw_py
residues_near_residue  = residues_near_residue_py
residues_near_residues = residues_near_residues_py
residues_near_position = residues_near_position_py
get_environment_distances_representation = get_environment_distances_representation_py
all_residues_with_serial_numbers = all_residues_with_serial_numbers_py
refine_zone_with_full_residue_spec = refine_zone_with_full_residue_spec_py
morph_fit_residues     = morph_fit_residues_py
find_blobs             = find_blobs_py
water_chain_from_shelx_ins = water_chain_from_shelx_ins_py
water_chain            = water_chain_py
make_variance_map      = make_variance_map_py
spin_search            = spin_search_py
spin_N                 = spin_N_py
ligand_search_make_conformers = ligand_search_make_conformers_py
cootaneer              = cootaneer_py
generic_string_vector_to_list_internal = generic_string_vector_to_list_internal_py
generic_list_to_string_vector_internal = generic_list_to_string_vector_internal_py
generic_int_vector_to_list_internal = generic_int_vector_to_list_internal_py
inverse_rtop           = inverse_rtop_py
# curl and hence coot_get_url_as_string are conditionally compiled.
try:
    coot_get_url_as_string = coot_get_url_as_string_py
    curl_progress_info     = curl_progress_info_py
except:
    pass # or print a message
score_rotamers         = score_rotamers_py
protein_db_loops       = protein_db_loops_py
make_link              = make_link_py
link_info              = link_info_py
map_to_model_correlation = map_to_model_correlation_py
map_to_model_correlation_stats = map_to_model_correlation_stats_py
map_to_model_correlation_per_residue = map_to_model_correlation_per_residue_py
qq_plot_map_and_model  = qq_plot_map_and_model_py
density_score_residue  = density_score_residue_py
map_mean               = map_mean_py
map_sigma              = map_sigma_py
map_statistics         = map_statistics_py
register_interesting_positions_list = register_interesting_positions_list_py
molecule_atom_overlaps = molecule_atom_overlaps_py
align_to_closest_chain = align_to_closest_chain_py
key_sym_code           = key_sym_code_py

# c-interface-ligands-swig.hh:
if (use_gui_qm != 2):
  new_molecule_sans_biggest_ligand = new_molecule_sans_biggest_ligand_py # comment out hack to load this smoothly.
  gui_ligand_metrics     = gui_ligand_metrics_py # comment out hack to load this smoothly.
residues_distortions   = residues_distortions_py
get_intermediate_atoms_distortions = get_intermediate_atoms_distortions_py
ligand_atom_overlaps   = ligand_atom_overlaps_py
residues_torsions_match = residues_torsions_match_py
kolmogorov_smirnov     = kolmogorov_smirnov_py
kolmogorov_smirnov_vs_normal = kolmogorov_smirnov_vs_normal_py
kullback_liebler       = kullback_liebler_py
coot_contact_dots_for_ligand =coot_contact_dots_for_ligand_py
switch_HIS_protonation = switch_HIS_protonation_py
linked_residues        = linked_residues_py
get_ligand_distortion_summary_info = get_ligand_distortion_summary_info_py

add_target_position_restraint_for_intermediate_atom = add_target_position_restraint_for_intermediate_atom_py
add_target_position_restraints_for_intermediate_atoms = add_target_position_restraints_for_intermediate_atoms_py

coot_contact_dots_for_ligand = coot_contact_dots_for_ligand_py
switch_HIS_protonation = switch_HIS_protonation_py

# probe-clash-score.hh
probe_clash_score      = probe_clash_score_py
# probe_clash_score_as   = probe_clash_score_as_py # comment out hack to load this smoothly.

# c-interface-generic-objects.h
generic_object_name    = generic_object_name_py

# and some acronyms
de_chainsaw                    = fill_partial_residues

# fix typo of set_find_hydrogen_torsions (for backwards compatibility?!)
set_find_hydrogen_torsion = set_find_hydrogen_torsions

# redefine in case someone has old function
toggle_idle_ligand_interactions = toggle_flev_idle_ligand_interactions
