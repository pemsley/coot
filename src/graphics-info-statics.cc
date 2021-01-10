
#include "Python.h"
#include "graphics-info.h"

bool graphics_info_t::residue_type_selection_was_user_picked_residue_range = false;

bool graphics_info_t::make_auto_h_bond_restraints_flag = false;

bool graphics_info_t::do_rotamer_restraints = false;

bool graphics_info_t::do_debug_refinement = false;

std::atomic<bool> graphics_info_t::on_going_updating_map_lock(false);

std::string graphics_info_t::mtz_file_for_refmac;

bool graphics_info_t::convert_dictionary_planes_to_improper_dihedrals_flag = false;

bool graphics_info_t::draw_missing_loops_flag = true;

bool graphics_info_t::sequence_view_is_docked_flag = true;

float graphics_info_t::pull_restraint_neighbour_displacement_max_radius = 1.0;

bool graphics_info_t::draw_stick_mode_atoms_default = true;

bool graphics_info_t::auto_recontour_map_flag = true;
