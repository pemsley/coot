
def add_pyranose_pseudo_ring_plane_restraints(comp_id):

    import re
    
    def filter_out(plane_name_sub_string, plane_restraints):
        return filter(lambda s: not plane_name_sub_string in s[0],
                      plane_restraints)
        
    restraints = monomer_restraints(comp_id)

    if (not isinstance(restraints, dict)):
        print "failed to get %s restraints" %comp_id
    else:
        plane_restraints = restraints["_chem_comp_plane_atom"]
        new_plane_restraints = [["pseudo-ring-1", [" C1 ", " C2 ", " C4 ", " C5 "], 0.01],
                                ["pseudo-ring-2", [" C2 ", " C3 ", " C5 ", " O5 "], 0.01],
                                ["pseudo-ring-3", [" C3 ", " C4 ", " O5 ", " C1 "], 0.01]] + \
                                filter_out("pseudo-ring-", plane_restraints) # should be list already 
        restraints["_chem_comp_plane_atom"] = new_plane_restraints

        set_monomer_restraints(comp_id, restraints)
        
def add_synthetic_pyranose_planes():
    for comp_id in ["NAG", "BMA", "MAN", "GAL", "GLC", "FUC", "XYP"]:
        add_pyranose_pseudo_ring_plane_restraints(comp_id)

def use_unimodal_pyranose_ring_torsions():
    for tlc in ["NAG", "BMA", "MAN", "GAL", "GLC", "FUC", "XYP"]:
        use_unimodal_ring_torsion_restraints(tlc)


def multi_add_linked_residue(imol, res_spec, residues_to_add):

    "---------------- multi-add-linked-residue", imol, res_spec
    set_go_to_atom_molecule(imol)
    wm = matrix_state() 
    set_matrix(wm/4.)

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
                                                      res_spec_to_chain_id(current_residue_spec),
                                                      res_spec_to_res_no(current_residue_spec),
                                                      res_spec_to_ins_code(current_residue_spec),
                                                      new_res, new_link, 2) # add and fit mode


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
    set_matrix(wm)
    
# return the new molecule number
#
def new_molecule_from_this_glyco_tree():
    with UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
                                   aa_atom_name, aa_alt_conf, aa_res_spec]:
        tree_residues = glyco_tree_residues(aa_imol, aa_res_spec)
        imol = new_molecule_by_residue_specs(aa_imol, tree_residues)
        return imol


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
# Plant Hybrid also allows an alpha1-3 FUC
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
               ["ALPHA1-6", "FUC"],
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
                    [
                        [["BETA1-2", "NAG"],
                         ["BETA1-4", "GAL"],
                         [["ALPHA2-3", "SIA"],
                          ["ALPHA2-6", "SIA"]
                        ]],
                        [["BETA1-6", "NAG"],
                         ["BETA1-4", "GAL"]]
                    ]
                    
                ],
                [
                    ["ALPHA1-3", "MAN"],
                    [
                        [["BETA1-2", "NAG"],
                         ["BETA1-4", "GAL"],
                         [["ALPHA2-3", "SIA"],
                          ["ALPHA2-6", "SIA"]
                         ]],
                        [["BETA1-4", "NAG"],
                         ["BETA1-4", "GAL"]]
                    ]
                ],
                [
                    ["BETA1-4", "NAG"]
                ],
                [
                    ["XYP-BMA", "XYP"]
                ]
               ],
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
add_linked_residue_tree_correlation_cut_off = 0.50

def add_linked_residue_tree(imol, parent, tree):

    global add_linked_residue_tree_correlation_cut_off
    
    def centre_view_on_residue_centre(res_spec):
        res_centre = residue_centre(imol,
                                    residue_spec_to_chain_id(res_spec),
                                    residue_spec_to_res_no(res_spec),
                                    residue_spec_to_ins_code(res_spec))
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
                                              res_spec_to_chain_id(parent),
                                              res_spec_to_res_no(parent),
                                              res_spec_to_ins_code(parent),
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
            if not (len(glyco_tree) == 1):
                return False
            else:
                with UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no,
                                               aa_ins_code, aa_atom_name, aa_alt_conf, aa_res_spec]:
                    res_spec = aa_res_spec
                    rn = residue_name(imol,
                                      residue_spec_to_chain_id(res_spec),
                                      residue_spec_to_res_no(res_spec),
                                      residue_spec_to_ins_code(res_spec))
                    if not isinstance(rn, str):
                        return False
                    else:
                        return rn == "ASN"

    # main line of add_linked_residue_tree
    #
    add_synthetic_pyranose_planes()
    use_unimodal_pyranose_ring_torsions()
    set_refine_with_torsion_restraints(1)
    wm = matrix_state()
    set_matrix(wm / 4.)
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
        # print "::::::::::::::::: start-tree:", start_tree

        # why do we need both test here? 
        # 5n11 needs the is-just-an-ASN? test
        # 5n09/5wzy needs the null test. 
        # Hmmm.
        if not ((start_tree == []) or (is_just_an_ASN_qm(aa_imol, start_tree))):
            info_dialog("Must start on Single ASN")
            print "start_tree:", start_tree
        else:
            # ok, continue
            start_pos_view = add_view_here("Glyco Tree Start Pos")
            process_tree(parent, tree, func)
            go_to_view_number(start_pos_view, 0)
            with AutoAccept():
                refine_residues(aa_imol, glyco_tree_residues(aa_imol, aa_res_spec))

            # add a test here that the tree here (centre of screen) matches a known tree.
            #
            # and that each is 4C1 (or 1C4 for FUC?) (XYP?)
            #
            # validate build - no more
            # g = glyco_validate()
            # BL says:: lets not auto delete since we may remove sugars in
            # the middle of the tree (OK? doesnt take RSCC into account)

            # 20180330-PE OK.
            #
            # g.auto_delete_residues()
            # g.validation_dialog()

    # reset
    set_default_temperature_factor_for_new_atoms(previous_m)
    set_matrix(current_weight)


def add_linked_residue_with_extra_restraints_to_active_residue(new_res_type,
                                                               link_type):
    wm = matrix_state()
    set_matrix(wm / 8.)
    set_refine_with_torsion_restraints(1)
    set_add_linked_residue_do_fit_and_refine(0)
    with UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
                                   aa_atom_name, aa_alt_conf, aa_res_spec]:
        new_res_spec = add_linked_residue(aa_imol, aa_chain_id, aa_res_no,
                                          aa_ins_code, new_res_type, link_type, 2)
        if isinstance(new_res_spec, list):
            add_cho_restraints_for_residue(aa_imol, new_res_spec)
            # refine that
            with AutoAccept():
                residues = [aa_res_spec] + residues_near_residue(aa_imol, aa_res_spec, 1.9)
                refine_residues(aa_imol, residues)
    set_matrix(wm)

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
                                      "SIA", "GLC", "GAL"]:
                                residue_spec = [chain_id, res_no, ""]
                                delete_cho_ls.append(residue_spec)
#                now we have delete_residues, we don't need to delete them one by one
#                for cho_res_spec in delete_cho_ls:
#                    delete_residue(aa_imol,
#                                   residue_spec_to_chain_id(cho_res_spec),
#                                   residue_spec_to_res_no(cho_res_spec), "")
                delete_residues(aa_imol, delete_cho_ls)


# Extra functionsnto in scheme....
                                
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

