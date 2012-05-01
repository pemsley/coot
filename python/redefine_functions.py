#############
# some re-definitions from coot python functions
############
# sorted by header files:

# c-interface.h:
chain_id               = chain_id_py
remarks                = remarks_py
residue_centre         = residue_centre_py
coot_sys_build_type    = coot_sys_build_type_py
run_clear_backups      = run_clear_backups_py
test_internal          = test_internal_py
test_internal_single   = test_internal_single_py
get_map_colour         = get_map_colour_py
average_map            = average_map_py
refmac_parameters      = refmac_parameters_py
map_sigma              = map_sigma_py
map_parameters         = map_parameters_py
cell                   = cell_py
map_cell               = cell_py
get_refmac_sad_atom_info = get_refmac_sad_atom_info_py
origin_pre_shift       = origin_pre_shift_py
save_state_file_name   = save_state_file_name_py
run_state_file         = run_state_file_py
wrapped_create_run_state_file_dialog =  wrapped_create_run_state_file_dialog_py
centre_of_mass_string  = centre_of_mass_string_py
set_atom_attributes    = set_atom_attributes_py
symmetry_operators     = symmetry_operators_py
symmetry_operators_to_xHM = symmetry_operators_to_xHM_py
merge_molecules        = merge_molecules_py
alignment_results      = alignment_results_py
nearest_residue_by_sequence = nearest_residue_by_sequence_py
refine_residues        = refine_residues_py
copy_residue_range_from_ncs_master_to_chains = copy_residue_range_from_ncs_master_to_chains_py
# to be renamed later
matching_compound_names_from_dictionary = matching_compound_names_from_dictionary_py
comp_id2name           = comp_id_to_name_py
refine_residues_with_alt_conf        = refine_residues_with_alt_conf_py
refine_zone_with_score = refine_zone_with_score_py
regularize_residues    = regularize_residues_py
regularize_residues_with_alt_conf  = regularize_residues_with_alt_conf_py
refine_zone_with_score = refine_zone_with_score_py
regularize_zone_with_score = regularize_zone_with_score_py
change_chain_id_with_result        = change_chain_id_with_result_py
non_standard_residue_names         = non_standard_residue_names_py
chain_id_for_shelxl_residue_number = chain_id_for_shelxl_residue_number_py
gln_asn_b_factor_outliers          = gln_asn_b_factor_outliers_py
map_peaks              = map_peaks_py
map_peaks_near_point   = map_peaks_near_point_py
view_name              = view_name_py
view_description       = view_description_py
go_to_view             = go_to_view_py
movie_file_name_prefix = movie_file_name_prefix_py
execute_ligand_search  = execute_ligand_search_py
overlap_ligands        = overlap_ligands_py
analyse_ligand_differences = analyse_ligand_differences_py
additional_representation_info = additional_representation_info_py
rigid_body_refine_by_residue_ranges = rigid_body_refine_by_residue_ranges_py
find_terminal_residue_type = find_terminal_residue_type_py
cis_peptides           = cis_peptides_py
get_rotamer_name       = get_rotamer_name_py
missing_atom_info      = missing_atom_info_py
rotamer_graphs         = rotamer_graphs_py
add_alt_conf           = add_alt_conf_py
highly_coordinated_waters = highly_coordinated_waters_py
add_lsq_atom_pair      = add_lsq_atom_pair_py
apply_lsq_matches      = apply_lsq_matches_py
get_lsq_matrix         = get_lsq_matrix_py
make_image_raster3d    = make_image_raster3d_py
make_image_povray      = make_image_povray_py
raster_screen_shot     = raster_screen_shot_py
ncs_master_chains      = ncs_master_chains_py
ncs_chain_differences  = ncs_chain_differences_py
ncs_chain_ids          = ncs_chain_ids_py
ncs_ghosts             = ncs_ghosts_py
pucker_info            = pucker_info_py
sequence_info          = sequence_info_py
alignment_mismatches   = alignment_mismatches_py
do_clipped_surface     = do_clipped_surface_py
generic_object_name    = generic_object_name_py
user_mods              = user_mods_py
probe_available_p      = probe_available_p_py
drag_intermediate_atom = drag_intermediate_atom_py
mark_atom_as_fixed     = mark_atom_as_fixed_py
mark_multiple_atoms_as_fixed = mark_multiple_atoms_as_fixed_py
ccp4i_projects         = ccp4i_projects_py
add_dipole             = add_dipole_py
add_dipole_for_residues = add_dipole_for_residues_py
get_pkgdatadir         = get_pkgdatadir_py
handle_pisa_interfaces = handle_pisa_interfaces_py
pkgdatadir             = get_pkgdatadir_py
matching_compound_names_from_sbase = matching_compound_names_from_sbase_py
add_linked_residue     = add_linked_residue_py
all_molecule_rotamer_score = all_molecule_rotamer_score_py
all_molecule_ramachandran_score = all_molecule_ramachandran_score_py
user_defined_click     = user_defined_click_py

get_torsion            = get_torsion_py
set_torsion            = set_torsion_py
multi_residue_torsion  = multi_residue_torsion_py
multi_residue_torsion_fit  = multi_residue_torsion_fit_py
test_function          = test_function_py
display_maps           = display_maps_py
space_group            = space_group_py

# graphics_info.h:
#run_post_manipulation_hook = run_post_manipulation_hook_py

# c-interface-python.hh:
make_atom_spec         = make_atom_spec_py

# cc-interface.hh:
goto_next_atom_maybe   = goto_next_atom_maybe_py
goto_prev_atom_maybe   = goto_prev_atom_maybe_py
set_go_to_atom_from_res_spec = set_go_to_atom_from_res_spec_py
get_symmetry           = get_symmetry_py
map_colour_components  = map_colour_components_py
dictionaries_read      = dictionaries_read_py
cif_file_for_comp_id   = cif_file_for_comp_id_py
monomer_restraints     = monomer_restraints_py
set_monomer_restraints = set_monomer_restraints_py
list_nomenclature_errors = list_nomenclature_errors_py
atom_info_string       = atom_info_string_py
residue_info           = residue_info_py
residue_name           = residue_name_py
clear_and_update_molecule = clear_and_update_molecule_py
add_molecule           = add_molecule_py
active_residue         = active_residue_py
closest_atom           = closest_atom_py
residues_near_residue  = residues_near_residue_py
residues_near_position = residues_near_position_py
refine_zone_with_full_residue_spec = refine_zone_with_full_residue_spec_py
list_extra_restraints  = list_extra_restraints_py
delete_extra_restraint = delete_extra_restraint_py
water_chain_from_shelx_ins = water_chain_from_shelx_ins_py
water_chain            = water_chain_py
make_variance_map      = make_variance_map_py
spin_search            = spin_search_py
cootaneer              = cootaneer_py
generic_string_vector_to_list_internal = generic_string_vector_to_list_internal_py
generic_list_to_string_vector_internal = generic_list_to_string_vector_internal_py
generic_int_vector_to_list_internal = generic_int_vector_to_list_internal_py
inverse_rtop           = inverse_rtop_py
protein_db_loops       = protein_db_loops_py
make_link              = make_link_py
key_sym_code           = key_sym_code_py
screen_vectors         = screen_vectors_py

# curl and hence coot_get_url_as_string are conditionally compiled.
try:
    coot_get_url_as_string = coot_get_url_as_string_py
    curl_progress_info     = curl_progress_info_py
except:
    pass # or print a message


# and some acronyms
de_chainsaw                    = fill_partial_residues

# fix typo of set_find_hydrogen_torsions (for backwards compatibility?!)
set_find_hydrogen_torsion = set_find_hydrogen_torsions
