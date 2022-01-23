
add_key_binding("Refine Active Residue", "r", lambda: manual_refine_residues(0))
add_key_binding("Refine Active Residue AA", "x", lambda: refine_active_residue())
add_key_binding("Triple Refine", "t", lambda: manual_refine_residues(1))
add_key_binding("Autofit Rotamer", "j", lambda: auto_fit_rotamer_active_residue())
add_key_binding("Pepflip", "q", lambda: pepflip_active_residue())
add_key_binding("Go To Blob", "g", lambda: blob_under_pointer_to_screen_centre())
add_key_binding("Eigen-flip Ligand", "e", lambda: flip_active_ligand())
add_key_binding("Add Water", "w", lambda: place_typed_atom_at_pointer("Water"))

def add_water_in_blob():
    blob_under_pointer_to_screen_centre()
    place_typed_atom_at_pointer("Water")
    refine_active_residue()
add_key_binding("Add Water + (centre+refine)", "W", lambda: add_water_in_blob())

def key_binding_func_1():
    active_atom = active_residue()
    if (not active_atom):
        print "No active atom"
    else:
        imol      = active_atom[0]
        chain_id  = active_atom[1]
        res_no    = active_atom[2]
        ins_code  = active_atom[3]
        atom_name = active_atom[4]
        alt_conf  = active_atom[5]
        add_terminal_residue(imol, chain_id, res_no, "auto", 1)
add_key_binding("Add terminal residue", "y", lambda: key_binding_func_1())

def key_binding_func_2():
    active_atom = active_residue()
    if (not active_atom):
        print "No active atom"
    else:
        imol      = active_atom[0]
        chain_id  = active_atom[1]
        res_no    = active_atom[2]
        ins_code  = active_atom[3]
        atom_name = active_atom[4]
        alt_conf  = active_atom[5]
        fill_partial_residue(imol, chain_id, res_no, ins_code)
add_key_binding("Fill Partial", "k", lambda: key_binding_func_2())

def delete_residue_sidechain_key():
    with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                               aa_ins_code, aa_atom_name, aa_alt_conf]:
        delete_residue_sidechain(aa_imol, aa_chain_id, aa_res_no,
                                 aa_ins_code, 0)
add_key_binding("Delete Residue Sidechain", "K",
                lambda: delete_residue_sidechain_key())

def rotamer_dialog_for_ar():
    with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                               aa_ins_code, aa_atom_name, aa_alt_conf]:
        show_rotamers_dialog(aa_imol, aa_chain_id, aa_res_no,
                             aa_ins_code, aa_alt_conf)
add_key_binding("Rotamers dialog for Active Residue", "Q",
                lambda: rotamer_dialog_for_ar())

def key_binding_func_5():
    active_atom = active_residue()
    if (not active_atom):
        add_status_bar_text("No active residue")
    else:
        imol      = active_atom[0]
        chain_id  = active_atom[1]
        res_no    = active_atom[2]
        ins_code  = active_atom[3]
        atom_name = active_atom[4]
        alt_conf  = active_atom[5]
        name = get_rotamer_name(imol, chain_id, res_no, ins_code)
        if (not name):
            add_status_bar_text("No Name found")
        else:
            if (name == ""):
                add_status_bar_text("No name for this")
            else:
                add_status_bar_text("Rotamer name: " + name)
add_key_binding("Rotamer name in Status Bar", "~", lambda: key_binding_func_5())

refine_residue_sphere_radius_key = 4.5  # Angstroms
add_key_binding("Refine residue in a sphere", "R",
                lambda: sphere_refine(refine_residue_sphere_radius_key))

def key_binding_func_21():
    if not valid_map_molecule_qm(imol_refinement_map()):
        info_dialog("Must set the refinement map")
    else:
        # not using active atom
        active_atom = active_residue()
        if (not active_atom):
            add_status_bar_text("No active residue")
        else:
            imol      = active_atom[0]
            chain_id  = active_atom[1]
            res_no    = active_atom[2]
            ins_code  = active_atom[3]
            atom_name = active_atom[4]
            alt_conf  = active_atom[5]

            rc_spec = [chain_id, res_no, ins_code]
            ls = residues_near_residue(imol, rc_spec, 1.9)
            with_auto_accept([refine_residues, imol, [rc_spec] + ls])
add_key_binding("Neighbours Refine", "h", lambda: key_binding_func_21())

add_key_binding("Regularize Residues in sphere", "B",
                lambda: sphere_regularize(refine_residue_sphere_radius_key))

def edit_chi_angles_key_func():
    with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                               aa_ins_code, aa_atom_name, aa_alt_conf]:
        edit_chi_angles(aa_imol, aa_chain_id, aa_res_no,
                        aa_ins_code, aa_alt_conf)
        set_moving_atom_move_chis()
        set_graphics_edit_current_chi(1)
add_key_binding("Edit Chi Angles", "X", lambda: edit_chi_angles_key_func())

# BL says:: I like to keep edit chi, so use "Y" for something potentially
# more useful (according to Paule)
add_key_binding("Just One or Next Map", "Y", lambda: display_cycle_through_maps())

def jiggle_fit_residue_key():
    with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                               aa_ins_code, aa_atom_name, aa_alt_conf]:
        fit_to_map_by_random_jiggle(aa_imol, aa_chain_id, aa_res_no,
                                    aa_ins_code, 100, 1.0)
add_key_binding("Jiggle Fit Residue", "J", lambda: jiggle_fit_residue_key())

def step_scrollable_map_number():
    maps = map_molecule_list()
    current = scroll_wheel_map()
    if maps:
        l = maps.index(current)
        if (l == len(maps) - 1):
            new = maps[0]
        else:
            new = maps[l+1]
        set_scroll_wheel_map(new)
add_key_binding("Step scrollable map number", "M",
                lambda: step_scrollable_map_number())

def delete_residue_hydrogens_key():
    with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                               aa_ins_code, aa_atom_name, aa_alt_conf]:
        delete_residue_hydrogens(aa_imol, aa_chain_id, aa_res_no,
                                 aa_ins_code, aa_alt_conf)
add_key_binding("Delete Residue Hydrogens", "P",
                lambda: delete_residue_hydrogens_key())

def key_binding_func_9():
    with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                               aa_ins_code, aa_atom_name, aa_alt_conf]:
        additional_representation_by_attributes(aa_imol, aa_chain_id,
                                                aa_res_no, aa_res_no, aa_ins_code,
                                                2, 2, 0.15, 1)
add_key_binding("ball-and-stickify residue", "$", lambda: key_binding_func_9())

add_key_binding("Undo Symmetry View", "V", lambda: undo_symmetry_view())

add_key_binding("Accept Baton Position", "A", lambda: accept_baton_position())

add_key_binding("Cootilus here", "N", lambda: find_nucleic_acids_local(6.0))

def jed_flip_key_func(dir):
    status = jed_flip_intermediate_atoms()
    if status == 0:
       active_atom = active_residue()
       if (not active_atom):
           print "No active atom"
       else:
           imol      = active_atom[0]
           chain_id  = active_atom[1]
           res_no    = active_atom[2]
           ins_code  = active_atom[3]
           atom_name = active_atom[4]
           alt_conf  = active_atom[5]
           jed_flip(imol, chain_id, res_no, ins_code, atom_name, alt_conf, dir)

add_key_binding("JED-Flip", "F", lambda: jed_flip_key_func(0))

add_key_binding("Reverse JED-Flip", "G", lambda: jed_flip_key_func(1))

# Paul's not sure about this one. I likey!
# add_key_binding("Delete this water", "D", lambda: delete_atom(*active_residue()))

