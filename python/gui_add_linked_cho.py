
def add_pyranose_pseudo_ring_plane_restraints(comp_id):

    def filter_out(plane_name_sub_string, plane_restraints):
        # BL says:: still something missing!?
        return plane_restraints
        
    restraints = monomer_restraints(comp_id)

    if (not isinstance(restraints, dict)):
        print "failed to get %s restraints" %comp_id
    else:
        plane_restraints = restraints["_chem_comp_plane_atom"]
        new_plane_restraints = ["_chem_comp_plane_atom",
                                [["pseudo-ring-1", [" C1 ", " C2 ", " C4 ", " C5 "], 0.01],
                                 ["pseudo-ring-2", [" C2 ", " C3 ", " C5 ", " O5 "], 0.01],
                                 ["pseudo-ring-3", [" C3 ", " C4 ", " O5 ", " C1 "], 0.01]]] + \
                                filter_out("pseudo-ring-", plane_restraints) #should be list already
        restraints["_chem_comp_plane_atom"] = new_plane_restraints
        print "BL DEBUG:: new plane restraints", new_plane_restraints

        set_monomer_restraints(comp_id, restraints)
        
def add_synthetic_pyranose_planes():
    for comp_id in ["NAG", "BMA", "MAN", "GAL", "GLC", "FUC", "XYP"]:
        add_pyranose_pseudo_ring_plane_restraints(comp_id)

def use_unimodal_pyranose_ring_torsions():
    for tlc in ["NAG", "BMA", "MAN", "GAL", "GLC", "FUC", "XYP"]:
        use_unimodal_ring_torsion_restraints(tlc)

# to be filled later
def add_cho_restraints_for_residue(imol, new_res_spec):
    return False


def multi_add_linked_residue(imol, res_spec, residues_to_add):

    "---------------- multi-add-linked-residue", imol, res_spec
    set_go_to_atom_molecule(imol)
    current_matrix = matrix_state() # BL says:: maybe should save
                                    # and return to original weight?!
    set_matrix(8)

    current_refinement_rate = dragged_refinement_steps_per_frame()

    current_residue_spec = res_spec
    for residue_to_add in residues_to_add:

        set_go_to_atom_from_res_spec(current_residue_spec)
        set_dragged_refinement_steps_per_frame(300)

        if (not isinstance(current_residue_spec, list)):
            print "OOps not a proper res_spec %s with residue_to_add: %s" \
                  %(current_residue_spec, residue_to_add)
            return False
        else:
            if not residue_to_add:
                return "done" # odd return?! why?
            else:
                if (not isinstance(residue_to_add, list) and \
                    not len(residue_to_add) == 2):
                    print "Oops - not a residue_link string pair when adding new_res_pair"
                else:
                    new_res = residue_to_add[0]
                    new_link = residue_to_add[1]

                    new_res_spec = add_linked_residue(imol,
                                                      res_spec2chain_id(current_residue_spec),
                                                      res_spec2res_no(current_residue_spec),
                                                      res_spec2ins_code(current_residue_spec),
                                                      new_res, new_link, 2)


                    # residues_near_residue takes a 3-part spec and makes 3-part specs
                    if new_res_spec:
                        ls = residues_near_residue(imol, current_residue_spec, 1.9)

                        add_cho_restraints_for_residue(imol, new_res_spec)
                        
                        with AutoAccept():
                            refine_residues(imol, [current_residue_spec] + ls)
                        # make the new res the one to add to (remove starting bool)
                        current_residue_spec = new_res_spec[1:]  

    # restore refinement mode and matrix weight
    set_dragged_refinement_steps_per_frame(current_refinement_rate)
    set_matrix(current_matrix)


# also "precursor"
# high mannose was used for human
# high mannose is now used for human too
# call this "High Mannose" in the GUI
# 
def oligomannose_tree():
    ret = [["NAG-ASN", "NAG"],
           ["BETA1-4", "NAG"],
           ["BETA1-4", "BMA"],
           [
               ["ALPHA1-6", "MAN"],
               [
                   ["ALPHA1-6", "MAN"],
                   ["ALPHA1-3", "MAN"]
                   ],
               [
                   ["ALPHA1-3", "MAN"],
                   ["ALPHA1-2", "MAN"]
                   ]
               ],
           [
               ["ALPHA1-3", "MAN"],
               ["ALPHA1-2", "MAN"],
               ["ALPHA1-2", "MAN"],
               ["ALPHA1-3", "GLC"],
               ["ALPHA1-3", "GLC"],
               ["ALPHA1-2", "GLC"]
               ]
           ]
    return ret

# Hybrid is for any system
#
# Plant Hybrid als allows an alpha1-3 FUC
#

# hybrid mammal
#
def hybrid_mammal_tree():
    ret = [["NAG-ASN", "NAG"],
           [
               [["BETA1-4", "NAG"],
                ["BETA1-4", "BMA"],
                [
                    [
                        ["ALPHA1-6", "MAN"],
                        [
                            ["ALPHA1-6", "MAN"],
                            ["ALPHA1-3", "MAN"]
                        ]
                    ],
                    [
                        ["ALPHA1-3", "MAN"],
                        ["BETA1-2", "NAG"],
                        ["BETA1-4", "GAL"],
                        [
                            ["ALPHA2-3, SIA"],
                            ["ALPHA2-6, SIA"]
                        ]
                    ],
                    ["BETA1-4", "NAG"]
                ]
               ],
               ["ALPHA1-6", "FUC"]
           ]
    ]
    return ret

# hybrid plant
#
def hybrid_plant_derived_tree():
    ret = [["NAG-ASN", "NAG"],
           [
               [["BETA1-4", "NAG"],
                ["BETA1-4", "BMA"],
                [
                    [
                        ["ALPHA1-6", "MAN"],
                        [
                            ["ALPHA1-6", "MAN"],
                            ["ALPHA1-3", "MAN"]
                        ]
                    ],
                    [
                        ["ALPHA1-3", "MAN"],
                        ["BETA1-2", "NAG"],
                        ["BETA1-4", "GAL"],
                        [
                            ["ALPHA2-3, SIA"],
                            ["ALPHA2-6, SIA"]
                        ]
                    ],
                    ["XYP-BMA", "XYP"],
                    ["BETA1-4", "NAG"]
                ]
               ],
               ["ALPHA1-6", "FUC"]
               ["ALPHA1-3", "FUC"]
           ]
    ]
    return ret

# complex mammal
# biantennary mammal
def complex_mammal_tree():
    ret = [["NAG-ASN", "NAG"],
           [
               [["BETA1-4", "NAG"],
                ["BETA1-4", "BMA"],
                [
                    ["ALPHA1-6", "MAN"],
                    ["BETA1-2", "NAG"],
                    ["BETA1-4", "GAL"],
                    [["ALPHA2-3", "SIA"],
                     ["ALPHA2-6", "SIA"]
                    ]
                    
                ],
                [
                    ["ALPHA1-3", "MAN"],
                    ["BETA1-2", "NAG"],
                    ["BETA1-4", "GAL"],
                    [["ALPHA2-3", "SIA"],
                     ["ALPHA2-6", "SIA"]
                    ]
                ],
                [
                    ["BETA1-4", "NAG"]
                ]
               ],
               ["ALPHA1-6", "FUC"]
           ]
    ]
    return ret

# complex plant
# plant biantennary
def complex_plant_tree():
    ret = [["NAG-ASN", "NAG"],
           [
               [["BETA1-4", "NAG"],
                ["BETA1-4", "BMA"],
                [
                    ["ALPHA1-6", "MAN"],
                    ["BETA1-2", "NAG"],
                    ["BETA1-4", "GAL"],
                    [["ALPHA2-3", "SIA"],
                     ["ALPHA2-6", "SIA"]
                    ]
                    
                ],
                [
                    ["ALPHA1-3", "MAN"],
                    ["BETA1-2", "NAG"],
                    ["BETA1-4", "GAL"],
                    [["ALPHA2-3", "SIA"],
                     ["ALPHA2-6", "SIA"]
                    ]
                ],
                [
                    ["BETA1-4", "NAG"],
                    ["BETA1-2", "XYP"]  # change link!?
                ]
               ],
               ["ALPHA1-6", "FUC"],
               ["ALPHA1-3", "FUC"]
           ]
    ]
    return ret

# old
def paucimannose_tree():
    ret = [["NAG-ASN", "NAG"],
           [  ["ALPHA1-3", "FUC"] ],
           [  ["BETA1-4", "NAG"],
              ["BETA1-4", "BMA"],
              [ ["XYP-BMA",  "XYP"] ],
              [ ["ALPHA1-3", "MAN"] ],
              [ ["ALPHA1-6", "MAN"] ]
           ]
           ]
    return ret

# now can be user defined
global add_linked_residue_tree_correlation_cut_off
add_linked_residue_tree_correlation_cut_off = 0.6

def add_linked_residue_tree(imol, parent, tree):

    global add_linked_residue_tree_correlation_cut_off
    
    def centre_view_on_residue_centre(res_spec):
        res_centre = residue_centre(imol,
                                    residue_spec2chain_id(res_spec),
                                    residue_spec2res_no(res_spec),
                                    residue_spec2ins_code(res_spec))
        if (isinstance(res_centre, list)):
            set_rotation_centre(*res_centre)
    
    def well_fitting_qm(res_spec):
        with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                   aa_ins_code, aa_atom_name, aa_alt_conf]:
            neighbs = residues_near_residue(aa_imol, res_spec, 4)
            c = map_to_model_correlation(imol, [res_spec], neighbs, 0,
                                         imol_refinement_map())
            print "######## new residue %s density correlation: %s" %(res_spec, c)
            if (not isNumber(c)):
                return False
            else:
                return c > add_linked_residue_tree_correlation_cut_off
    
    def func(parent, res_pair):
        if (not isinstance(parent, list)):
            print "WARNING:: OOps not a proper res_spec %s with residue_to_add: %s" %(parent, res_pair)
            return False
        else:
            # OK. go go go
            new_link = res_pair[0]
            new_res_type = res_pair[1]

            centre_view_on_residue_centre(parent)

            # tree_residues is empty (excludes single ASN) on startup,
            # so should add this
            tree_residues = glyco_tree_residues(imol, parent)
            if is_just_an_ASN_qm(imol, tree_residues):
                tree_residues = [parent]
            imol_save = new_molecule_by_residue_specs(imol,
                                                      tree_residues)
            # old
#            active_atom = active_residue_py()
#            active_residue = active_atom[:4]
#            save_residue_specs = glyco_tree_residues_py(imol, active_residue)
#
#            imol_glyco_pre = new_molecule_by_residue_specs_py(imol, save_residue_specs)
#            set_mol_displayed(imol_glyco_pre, 0)
#            set_mol_active(imol_glyco_pre, 0)

            new_res_spec = add_linked_residue(imol,
                                              res_spec2chain_id(parent),
                                              res_spec2res_no(parent),
                                              res_spec2ins_code(parent),
                                              new_res_type, new_link, 2)
            set_mol_displayed(imol_save, 0)
            set_mol_active(imol_save, 0)
            ls = residues_near_residue(imol, parent, 1.9)
            local_ls = [parent] + ls
            add_cho_restraints_for_residue(imol, new_res_spec)
            rotate_y_scene(100, 0.5)
            with AutoAccept():
                refine_residues(imol, local_ls)
            if (not isinstance(new_res_spec, list)):
                # badness
                return False
            else:
                # okay!?
                preped_new_res_spec = new_res_spec[1:]  # strip off leading result
                if well_fitting_qm(preped_new_res_spec):
                    # close_molecule(imol_glyco_pre)
                    return preped_new_res_spec
                else:
                    # ------------ bad fit -----------------
                    # delete residue and restore others
                    print "----------- That was not well fitting. Deleting:", preped_new_res_spec
                    delete_extra_restraints_for_residue_spec(imol,
                                                             preped_new_res_spec)
                    delete_residue_by_spec(imol, preped_new_res_spec)
                    # restore glyco-tree residues from imol_save
                    write_pdb_file(imol_save, "test_glyco.pdb")
                    replace_fragment(imol, imol_save, "//")
                    return False
#                    with AutoAccept():
#                        # Note: may not get rid of screwed up refinement from
#                        # adding a carb too much...
#                        refine_residues(imol, local_ls)
#                        return False
                    
    def process_tree(parent, tree, proc_func):
        for branch in tree:
            if ((len(branch) == 2) and
                isinstance(branch[0], str) and
                isinstance(branch[0], str)):
                # have the link, res pair!
                new_res = proc_func(parent, branch)
                if new_res:
                    parent = new_res
                else:
                    break
            else:
                if isinstance(branch, list):
                    process_tree(parent, branch, proc_func)

    def is_just_an_ASN_qm(imol, glyco_tree):
        if not isinstance(glyco_tree, list):
            return False
        else:
            if not (len(glyco_tree) == 0):
                return False
            else:
                with UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no,
                                               aa_ins_code, aa_atom_name, aa_alt_conf, aa_res_spec]:
                    res_spec = aa_res_spec
                    rn = residue_name(imol,
                                      residue_spec2chain_id(res_spec),
                                      residue_spec2res_no(res_spec),
                                      residue_spec2ins_code(res_spec))
                    if not isinstance(rn, str):
                        return False
                    else:
                        return rn == "ASN"

    # main line of add_linked_residue_tree
    #
    add_synthetic_pyranose_planes()
    use_unimodal_pyranose_ring_torsions()
    set_refine_with_torsion_restraints(1)
    current_weight = matrix_state()
    set_matrix(8)
    set_residue_selection_flash_frames_number(1)
    set_go_to_atom_molecule(imol)
    previous_m = default_new_atoms_b_factor()
    m = median_temperature_factor(imol)
    try:
        new_m = m * 1.55
        if new_m > 0:
            set_default_temperature_factor_for_new_atoms(new_m)
    except:
        print "BL INFO:: not changing default B for new carb atoms"

    # prevent tree building if there is already a partial tree here
    # (only proceed with the one ASN)
    #
    with UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
                                   aa_atom_name, aa_alt_conf, aa_res_spec]:
        start_tree = glyco_tree_residues(aa_imol, aa_res_spec)
        print "::::::::::::::::: start-tree:", start_tree

        if not is_just_an_ASN_qm(aa_imol, start_tree):
            info_dialog("Must start on Single ASN")
        else:
            # ok, continue
            start_pos_view = add_view_here("Glyo Tree Start Pos")
            process_tree(parent, tree, func)
            go_to_view_number(start_pos_view, 0)
            with AutoAccept():
                refine_residues(aa_imol, glyco_tree_residues(aa_imol, aa_res_spec))
            # validate build
            g = glyco_validate()
            g.auto_delete_residues()
    # reset
    set_default_temperature_factor_for_new_atoms(previous_m)
    set_matrix(current_weight)


def add_linked_residue_with_extra_restraints_to_active_residue(new_res_type,
                                                               link_type):
    current_weight = matrix_state()
    set_matrix(8)
    set_refine_with_torsion_restraints(1)
    set_add_linked_residue_do_fit_and_refine(0)
    with UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
                                   aa_atom_name, aa_alt_conf, aa_res_spec]:
        new_res_spec = add_linked_residue(aa_imol, aa_chain_id, aa_res_no,
                                          aa_ins_code, new_res_type, link_type, 2)
        add_cho_restraints_for_residue(aa_imol, new_res_spec)
        # refine that
        with AutoAccept():
            residues = [aa_res_spec] + residues_near_residue(aa_imol, aa_res_spec, 1.9)
            refine_residues(aa_imol, residues)

def delete_all_cho():
    delete_cho_ls = []
    with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
                               aa_atom_name, aa_alt_conf]:
        if valid_model_molecule_qm(aa_imol):
            with NoBackups(aa_imol):
                for chain_id in chain_ids(aa_imol):
                    for res_serial in range(chain_n_residues(chain_id, aa_imol)):
                        res_no = seqnum_from_serial_number(aa_imol,
                                                           chain_id, res_serial)
                        rn = residue_name(aa_imol, chain_id, res_no, "")
                        if (isinstance(rn, str)):
                            if rn in ["NAG", "MAN", "BMA", "FUC", "XYP",
                                      "SIA", "GLC"]:
                                residue_spec = [chain_id, res_no, ""]
                                delete_cho_ls.append(residue_spec)
#                now we have delete_residues, we don't need to delete them one by one
#                for cho_res_spec in delete_cho_ls:
#                    delete_residue(aa_imol,
#                                   residue_spec2chain_id(cho_res_spec),
#                                   residue_spec2res_no(cho_res_spec), "")
                delete_residues(aa_imol, delete_cho_ls)

def interactive_add_cho_dialog():

    add_synthetic_pyranose_planes()
    use_unimodal_pyranose_ring_torsions()
    # button list with [label, function]
    buttons = [
        ["Update for Current Residue", lambda func: printf("dummy")],
        ["Add a NAG-ASN NAG",
         lambda func:
         add_linked_residue_with_extra_restraints_to_active_residue("NAG", "NAG-ASN")],
        ["Add a BETA1-4 NAG",
         lambda func:
         add_linked_residue_with_extra_restraints_to_active_residue("NAG", "BETA1-4")],
        ["Add a BETA1-4 BMA",
         lambda func:
         add_linked_residue_with_extra_restraints_to_active_residue("BMA", "BETA1-4")],
        ["Add an ALPHA1-2 MAN",
         lambda func:
         add_linked_residue_with_extra_restraints_to_active_residue("MAN", "ALPHA1-2")],
        ["Add an ALPHA1-3 MAN",
         lambda func:
         add_linked_residue_with_extra_restraints_to_active_residue("MAN", "ALPHA1-3")],
        ["Add an ALPHA2-3 MAN",
         lambda func:
         add_linked_residue_with_extra_restraints_to_active_residue("MAN", "ALPHA2-3")],
        ["Add an ALPHA2-3 GAL",
         lambda func:
         add_linked_residue_with_extra_restraints_to_active_residue("GAL", "ALPHA2-3")],
        ["Add an ALPHA1-6 MAN",
         lambda func:
         add_linked_residue_with_extra_restraints_to_active_residue("MAN", "ALPHA1-6")],
        ["Add a BETA1-2 NAG",
         lambda func:
         add_linked_residue_with_extra_restraints_to_active_residue("NAG", "BETA1-2")],
        ["Add a BETA1-4 GAL",
         lambda func:
         add_linked_residue_with_extra_restraints_to_active_residue("GAL", "BETA1-4")],
        ["Add an ALPHA1-2 FUC",
         lambda func:
         add_linked_residue_with_extra_restraints_to_active_residue("FUC", "ALPHA1-2")],
        ["Add an ALPHA1-3 FUC",
         lambda func:
         add_linked_residue_with_extra_restraints_to_active_residue("FUC", "ALPHA1-3")],
        ["Add an ALPHA1-6 FUC",
         lambda func:
         add_linked_residue_with_extra_restraints_to_active_residue("FUC", "ALPHA1-6")],
        ["Add an BETA1-6 FUL",
         lambda func:
         add_linked_residue_with_extra_restraints_to_active_residue("FUL", "BETA1-6")],
        ["Add an XYP-BMA XYP",
         lambda func:
         add_linked_residue_with_extra_restraints_to_active_residue("XYP", "XYP-BMA")]
    ]

    vbox = dialog_box_of_buttons("Add N-linked Glycan",
                                 [360, 520], buttons, "Close")
    gui_add_linked_cho_dialog_vbox_set_rotation_centre_hook(vbox)
    # set the callback on the first button
    children = vbox.get_children()
    if isinstance(children, list):
        if len(children) > 0:
            first_button = children[0]
            first_button.connect("clicked", lambda func:
                                 gui_add_linked_cho_dialog_vbox_set_rotation_centre_hook(vbox))

    # add a widget to allow the user to choose the tree type
    hbox_1 = gtk.HBox(False, 2)
    hbox_2 = gtk.HBox(False, 2)
    butt_1 = gtk.RadioButton(None, "Oligomannose")
    butt_2 = gtk.RadioButton(butt_1, "Hybrid")
    butt_3 = gtk.RadioButton(butt_1, "Complex")
    butt_4 = gtk.RadioButton(butt_1, "Expert Mode")

    hbox_1.pack_start(butt_1, False, False, 2)
    hbox_1.pack_start(butt_2, False, False, 2)
    hbox_2.pack_start(butt_3, False, False, 2)
    hbox_2.pack_start(butt_4, False, False, 2)

    butt_1.show()
    butt_2.show()
    butt_3.show()
    butt_4.show()
    hbox_1.show()
    hbox_2.show()
    vbox.pack_start(hbox_1, False, False, 2)
    vbox.pack_start(hbox_2, False, False, 2)
    hbox_1.set_homogeneous(True)
    hbox_2.set_homogeneous(True)
    vbox.reorder_child(hbox_1, 0)
    vbox.reorder_child(hbox_2, 1)

    # "global" var post-set-rotation-centre-hook
    # BL Note:: maybe should be a global!?
    post_set_rotation_centre_script = gui_add_linked_cho_dialog_vbox_set_rotation_centre_hook(vbox)

#
def glyco_tree_dialog_set_button_active_state(button, glyco_id, tree_type):

    def get_sensitive_button_list(glyco_id, tree_type):

        if not isinstance(glyco_id, list):
            return []
        else:
            level_number = glyco_id[0]
            residue_type = glyco_id[1]
            link_type = glyco_id[2]
            parent_residue_type = glyco_id[3]
            residue_spec = glyco_id[4]
            active_button_label_ls = []
            
            if tree_type == 'expert-mode':
                active_button_label_ls = "expert-mode" # ???

            if tree_type == 'oligomannose':

                if level_number == 0:
                    if residue_type == "ASN":
                        active_button_label_ls = ["Add a NAG-ASN NAG"]
                        
                if level_number == 1:
                    if residue_type == "NAG":
                        active_button_label_ls = ["Add a BETA1-4 NAG"]
                        
                if level_number == 2:
                    if residue_type == "NAG":
                        active_button_label_ls = ["Add a BETA1-4 BMA"]

                if level_number == 3:
                    if residue_type == "BMA":
                        active_button_label_ls = ["Add an ALPHA1-3 MAN",
                                                  "Add an ALPHA1-6 MAN"]

                if level_number == 4:
                    if residue_type == "MAN":
                        if link_type == "ALPHA1-3":
                            active_button_label_ls = ["Add an ALPHA1-2 MAN"]
                        if link_type == "ALPHA1-6":
                            active_button_label_ls = ["Add an ALPHA1-3 MAN",
                                                      "Add an ALPHA1-6 MAN"]
                        
                if level_number == 5:
                    if residue_type == "MAN":
                        # active_button_label_ls is the same, so no if required?
                        # or wrong tree branches. FIXME.
                        if link_type == "ALPHA1-2":
                            active_button_label_ls = ["Add an ALPHA1-2 MAN"]
                        if link_type == "ALPHA1-6":
                            active_button_label_ls = ["Add an ALPHA1-2 MAN"]
                        if link_type == "ALPHA1-3":
                            active_button_label_ls = ["Add an ALPHA1-2 MAN"]
                            
                if level_number == 6:
                    if residue_type == "MAN":
                        if link_type == "ALPHA1-2":
                            active_button_label_ls = ["Add an ALPHA1-3 GLC"]

                # inconsistencies between links here. FIXME
                if level_number == 7:
                    if residue_type == "GLC":
                        if link_type == "ALPHA1-2":
                            active_button_label_ls = ["Add an ALPHA1-3 GLC"]

                if level_number == 8:
                    if residue_type == "GLC":
                        if link_type == "ALPHA1-2":
                            active_button_label_ls = ["Add an ALPHA1-2 GLC"]

            # hybrid
            if tree_type == 'hybrid':

                if level_number == 0:
                    if residue_type == "ASN":
                        active_button_label_ls = ["Add a NAG-ASN NAG"]
                        
                if level_number == 1:
                    if residue_type == "NAG":
                        active_button_label_ls = ["Add a BETA1-4 NAG",
                                                  "Add an ALPHA1-3 FUC"]
                        
                if level_number == 2:
                    if residue_type == "NAG":
                        active_button_label_ls = ["Add a BETA1-4 BMA"]

                if level_number == 3:
                    if residue_type == "BMA":
                        active_button_label_ls = ["Add an ALPHA1-3 MAN",
                                                  "Add an ALPHA1-6 MAN"]

                if level_number == 4:
                    if residue_type == "MAN":
                        active_button_label_ls = ["Add an ALPHA1-3 MAN",
                                                  "Add an ALPHA1-6 MAN"]
                
            # complex
            if tree_type == 'complex':

                if level_number == 0:
                    if residue_type == "ASN":
                        active_button_label_ls = ["Add a NAG-ASN NAG"]
                        
                if level_number == 1:
                    if residue_type == "NAG":
                        active_button_label_ls = ["Add a BETA1-4 NAG",
                                                  "Add an ALPHA1-6 FUC"]
                        
                if level_number == 2:
                    if residue_type == "NAG":
                        active_button_label_ls = ["Add a BETA1-4 BMA"]

                if level_number == 3:
                    if residue_type == "BMA":
                        active_button_label_ls = ["Add an ALPHA1-3 MAN",
                                                  "Add an ALPHA1-6 MAN",
                                                  "Add an XYP-BMA XYP"]
                        
            return active_button_label_ls

    # main line
    #
    l = button.get_label()
    active_button_label_ls = get_sensitive_button_list(glyco_id, tree_type)
    if (active_button_label_ls == "expert-mode"):
        button.set_sensitive(True)
        if not (l == "Update for Current Residue"):
            button.set_sensitive(l in active_button_label_ls)

# return a value!?
#
def gui_add_linked_cho_dialog_vbox_set_rotation_centre_hook(vbox):

    def get_tree_type():
        tree_type = "oligomannose"
        children = vbox.get_children()
        for child in children:
            if type(child) == gtk.Box:
                for box_child in child:
                    if type(box_child) == gtk.RadioButton:
                        if box_child.get_active():
                            l = box_child.get_label()
                            if l == "Oligomannose":
                                tree_type = 'oligomannose'
                            if l == "Hybrid":
                                tree_type = 'hybrid'
                            if l == "Expert Mode":
                                tree_type = 'expert-mode'
                            if l == "Complex":
                                tree_type = 'complex'
        return tree_type

    with UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
                                   aa_atom_name, aa_alt_conf, aa_res_spec]:
        glyco_id = glyco_tree_residue_id(aa_imol, aa_res_spec)
        # Paule says:
        # if it was an ASP create a level-0 glyco-id for that (glyco-tree-residue-id doesn't
        # do that (not sure why)).
        if not glyco_id:
            rn = residue_name(aa_imol, aa_chain_id, aa_res_no, aa_ins_code)
            if isinstance(rn, str):
                if rn == "ASN":
                    glyco_id = [0, "ASN", "", "", aa_res_spec]
        if isinstance(glyco_id, list):
            tree_type = get_tree_type()
            children = vbox.get_children()
            for child in children:
                if (type(child) == gtk.Button):
                    glyco_tree_dialog_set_button_active_state(child, glyco_id,
                                                              tree_type)
            return True
        return False

# return the new molecule number
#
def new_molecule_from_this_glyco_tree():
    with UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
                                   aa_atom_name, aa_alt_conf, aa_res_spec]:
        tree_residues = glyco_tree_residues(aa_imol, aa_res_spec)
        imol = new_molecule_by_residue_specs(aa_imol, tree_residues)
        return imol

                                
def delete_glyco_tree():
    
    active_atom = active_residue_py()
    try:
        imol = active_atom[0]
        active_residue = active_atom[:4]
        print "active_residue", active_residue
        glyco_tree_residues = glyco_tree_residues_py(imol, active_residue)
        print active_residue
        print "glyco_tree_residues", glyco_tree_residues
        for res in glyco_tree_residues:
            rn = residue_name_by_spec(imol, res)
            if rn != "ASN":
                delete_residue_by_spec(imol, res)
    except KeyError as e:
        print e
    except TypeError as e:
        print e

# this has to be at the top level for reasons that I don't understand
# (python function scoping) - I'll just put it down to python being shit as usual.
#
def glyco_validation_dialog_set_go_to_residue(imol, residue_spec):
    rc = residue_centre(imol,
                        res_spec_utils.residue_spec_to_chain_id(residue_spec),
                        res_spec_utils.residue_spec_to_res_no(residue_spec),
                        '')
    set_rotation_centre(*rc)

def load_privateer_dictionary():
    if os.path.exists("privateer-lib.cif"):
        read_cif_dictionary("privateer-lib.cif")
        set_refine_with_torsion_restraints(1)
    
class glyco_validate:

    def run_privateer(self, imol, glyco_tree_residues, hklin_fn, fp_col, sigfp_col, pdbin_fn, privateer_log):
        if is_valid_model_molecule(imol):
            args = ['-mtzin', hklin_fn, '-colin-fo', fp_col+','+sigfp_col,
                    '-pdbin', pdbin_fn]
            # using'-mode', 'ccp4i2' makes a file with no residue nambers (afaics)
            popen_command("privateer", args, [], privateer_log, False)

    def make_privateer_validation_info(self, imol, fp_col, sigfp_col, glyco_tree_residues):

        imol_map = imol_refinement_map()

        if len(glyco_tree_residues) > 0:
            d = get_directory("coot-ccp4")
            spid = str(os.getpid())
            fn_pdb = "coot-privateer-" + spid + ".pdb"
            fn_log = "coot-privateer-" + spid + ".log"
            privateer_pdb = os.path.join(d, fn_pdb)
            privateer_log = os.path.join(d, fn_log)
            hklin_fn = mtz_file_name(imol_map);
            write_pdb_file(imol, privateer_pdb)
            self.run_privateer(imol, glyco_tree_residues, hklin_fn, fp_col, sigfp_col, privateer_pdb, privateer_log)
            pvi = self.parse_privateer_log(privateer_log, imol, glyco_tree_residues)
            return pvi
        else:
            return []

    # return a list of residues with privateer validation info
    #
    def parse_privateer_log(self, log_file_name, imol, glyco_tree_residues):

        print 'parse_privateer_log', log_file_name, imol, glyco_tree_residues

        pvi = []
        f = open(log_file_name)
        lines = f.readlines()
        for line in lines:
            l = line.rstrip()
            if len(l) > 10:
                if l[:4] == 'coot':
                    words = l.split()
                    if len(words) > 12:
                        for r in glyco_tree_residues:
                            rn = residue_name_by_spec(imol, r)
                            try:
                                res_id = rn + "-" + res_spec_utils.residue_spec_to_chain_id(r) + \
                                         '-' + str(res_spec_utils.residue_spec_to_res_no(r))
                                # print "res_id", res_id
                                if words[1] == res_id:
                                    # print words[12] , yes or check
                                    new_item = (r, words)
                                    pvi.append(new_item)
                            except TypeError as e:
                                print e
        print "parsed", log_file_name
        return pvi

    def make_validation_dialog(self, imol, privateer_validation_info):

        print "make_validation_dialog", privateer_validation_info
        buttons = []
        for vi in privateer_validation_info:
            residue_spec = vi[0]
            words = vi[1]
            state = words[12]
            if state == "yes":
                state = " OK  " # label spacing
            button_text = words[1] + "   Q=" + words[3] + " RSCC=" + words[6] + \
                          " Cnf=" + words[8] + "   " + state
            # string->function? Hideous
            func = "glyco_validation_dialog_set_go_to_residue(" + str(imol) + "," + \
                   str(residue_spec) + ")"
            button = [button_text, func]
            buttons.append(button)
        if len(buttons) > 0:
            button = ["Load Privateer Dictionary", "load_privateer_dictionary()"]
            buttons.append(button)
            dialog_box_of_buttons("Privateer Validation", (400, 220), buttons, " Close ")

    def validation_dialog(self):
        
        active_atom = active_residue_py()
        try:
            imol = active_atom[0]
            active_residue = active_atom[:4]
            glyco_tree_residues = glyco_tree_residues_py(imol, active_residue)

            print 'imol', imol
            print 'active_residue', active_residue
            print 'glyco_tree_residues', glyco_tree_residues

            fp_col='FP'
            sigfp_col='SIGFP'    
            pvi = self.make_privateer_validation_info(imol, fp_col, sigfp_col, glyco_tree_residues)
            # print "validation_dialog()", pvi
            if len(pvi) > 0:
                # now make a gui
                self.make_validation_dialog(imol, pvi)

        except KeyError as e:
            print e

        # no active atom
        except TypeError as e:
            print e

    def auto_delete_residues_internal(self, imol, glyco_tree_residues):
        fp_col='FP'
        sigfp_col='SIGFP'
        pvi = self.make_privateer_validation_info(imol, fp_col, sigfp_col, glyco_tree_residues)
        for res_info in pvi:
            print "test", res_info
            res_status = res_info[1][12]
            if res_status == 'check' or res_status == 'no':
                delete_residue_by_spec(imol, res_info[0])

    def auto_delete_residues(self):
        try:
            active_atom = active_residue_py()
            imol = active_atom[0]
            active_residue = active_atom[:4]
            glyco_tree_residues = glyco_tree_residues_py(imol, active_residue)
            self.auto_delete_residues_internal(imol, glyco_tree_residues)
        except TypeError as e:
            print e
    
            
# graphics...

def add_module_carbohydrate():
    if (have_coot_python):
        if coot_python.main_menubar():
            menu = coot_menubar_menu("Glyco")

            add_simple_coot_menu_menuitem(
                menu, "N-linked Glycan Addition Dialog...",
                lambda func:
                interactive_add_cho_dialog())
            
            def add_multi_carbo_link_func(link_list):
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    multi_add_linked_residue(aa_imol,
                                             [aa_chain_id, aa_res_no, aa_ins_code],
                                             link_list)
                
            add_simple_coot_menu_menuitem(
                menu, "N-link add NAG, NAG, BMA",
                lambda func: add_multi_carbo_link_func([["NAG", "NAG-ASN"],
                                                        ["NAG", "BETA1-4"],
                                                        ["BMA", "BETA1-4"]]))
            
            # add_simple_coot_menu_menuitem(
            #     menu, "Add a ASN-NAG NAG",
            #     lambda func:
            #     add_linked_residue_with_extra_restraints_to_active_residue("NAG", "NAG-ASN"))

            # add_simple_coot_menu_menuitem(
            #     menu, "Add a BETA1-4 NAG",
            #     lambda func:
            #     add_linked_residue_with_extra_restraints_to_active_residue("NAG", "BETA1-4"))

            # add_simple_coot_menu_menuitem(
            #     menu, "Add a BETA1-4 BMA",
            #     lambda func:
            #     add_linked_residue_with_extra_restraints_to_active_residue("BMA", "BETA1-4"))
            
            # add_simple_coot_menu_menuitem(
            #     menu, "Add an ALPHA1-2 MAN",
            #     lambda func:
            #     add_linked_residue_with_extra_restraints_to_active_residue("MAN", "ALPHA1-2"))
            
            # add_simple_coot_menu_menuitem(
            #     menu, "Add an ALPHA1-3 MAN",
            #     lambda func:
            #     add_linked_residue_with_extra_restraints_to_active_residue("MAN", "ALPHA1-3"))

            # we should do this only if we are sitting on an SIA.
            # Attaching a SIA to a MAN (i.e. reverse order) would be a
            # good test too...
            # add_simple_coot_menu_menuitem(
            #     menu, "Add an ALPHA2-3 MAN",
            #     lambda func:
            #     add_linked_residue_with_extra_restraints_to_active_residue("MAN", "ALPHA2-3"))

            # # same consideration as above
            # add_simple_coot_menu_menuitem(
            #     menu, "Add an ALPHA2-3 GAL",
            #     lambda func:
            #     add_linked_residue_with_extra_restraints_to_active_residue("GAL", "ALPHA2-3"))

            # add_simple_coot_menu_menuitem(
            #     menu, "Add an ALPHA1-6 MAN",
            #     lambda func:
            #     add_linked_residue_with_extra_restraints_to_active_residue("MAN", "ALPHA1-6"))

            # add_simple_coot_menu_menuitem(
            #     menu, "Add an ALPHA1-3 FUC",
            #     lambda func:
            #     add_linked_residue_with_extra_restraints_to_active_residue("FUC", "ALPHA1-3"))

            # add_simple_coot_menu_menuitem(
            #     menu, "Add an ALPHA1-6 FUC",
            #     lambda func:
            #     add_linked_residue_with_extra_restraints_to_active_residue("FUC", "ALPHA1-6"))

            # add_simple_coot_menu_menuitem(
            #     menu, "Add an XYP-BMA XYP",
            #     lambda func:
            #     add_linked_residue_with_extra_restraints_to_active_residue("XYP", "XYP-BMA"))


            # the mode in the function call now takes take of this
            # add_simple_coot_menu_menuitem(
            #     menu, "Auto Fit & Refine On for Link Addition",
            #     lambda func: set_add_linked_residue_do_fit_and_refine(1))

            # add_simple_coot_menu_menuitem(
            #     menu, "Auto Fit & Refine Off for Link Addition",
            #     lambda func: set_add_linked_residue_do_fit_and_refine(0))

            def add_oligo_tree_func(oligo_tree):
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    make_backup(aa_imol)
                    # switch backup off?!
                    add_linked_residue_tree(aa_imol,
                                            [aa_chain_id, aa_res_no, aa_ins_code],
                                            oligo_tree)
                
            add_simple_coot_menu_menuitem(
                menu, "Add Oligomannose",
                lambda func: add_oligo_tree_func(oligomannose_tree()))

            add_simple_coot_menu_menuitem(
                menu, "Add Hybrid (Mammal)",
                lambda func: add_oligo_tree_func(hybrid_mammal_tree()))

            add_simple_coot_menu_menuitem(
                menu, "Add Hybrid (Plant)",
                lambda func: add_oligo_tree_func(hybrid_plant_derived_tree()))

            add_simple_coot_menu_menuitem(
                menu, "Add Complex",
                lambda func: add_oligo_tree_func(complex_mammal_tree()))

            add_simple_coot_menu_menuitem(
                menu, "Add Complex (Plant)",
                lambda func: add_oligo_tree_func(complex_plant_tree()))

            add_simple_coot_menu_menuitem(
                menu, "Delete All Carbohydrate",
                lambda func: delete_all_cho())

            def torsion_fit_this_func(refine = False):
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    centre_residue = [aa_chain_id,aa_res_no, aa_ins_code]
                    multi_residue_torsion_fit(aa_imol,
                                              [centre_residue],
                                              30000)
                    if refine:
                        with AutoAccept():
                            refine_residues(aa_imol, [centre_residue])

            def torsion_fit_this_and_neighbours_func(refine = False):
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    centre_residue = [aa_chain_id,aa_res_no, aa_ins_code]
                    residues = residues_near_residue(aa_imol, centre_residue, 1.9)
                    residues.append(centre_residue)
                    multi_residue_torsion_fit(aa_imol, residues, 30000)
                    if refine:
                        with AutoAccept():
                            refine_residues(aa_imol, [centre_residue])

            add_simple_coot_menu_menuitem(
                menu, "Torsion Fit this residue",
                lambda func: torsion_fit_this_func())

            # add_simple_coot_menu_menuitem(
            #     menu, "Torsion Fit This Residue and Neighbours",
            #     lambda func: torsion_fit_this_and_neighbours_func())

            add_simple_coot_menu_menuitem(
                menu, "Torsion Fit & Refine this residue",
                lambda func: torsion_fit_this_func(True))

            add_simple_coot_menu_menuitem(
                menu, "Add synthetic pyranose plane restraints",
                lambda func: add_synthetic_pyranose_planes())

            add_simple_coot_menu_menuitem(
                menu, "Use Unimodal ring torsion restraints",
                lambda func: use_unimodal_pyranose_ring_torsions())

            add_simple_coot_menu_menuitem(
                menu, "Extract this Tree",
                lambda func:
                new_molecule_from_this_glyco_tree())
