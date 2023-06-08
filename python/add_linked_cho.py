import coot
import coot_utils

def add_pyranose_pseudo_ring_plane_restraints(comp_id):

    import re

    def filter_out(plane_name_sub_string, plane_restraints):
        return [s for s in plane_restraints if not plane_name_sub_string in s[0]]

    restraints = coot.monomer_restraints_py(comp_id)

    if (not isinstance(restraints, dict)):
        print("failed to get %s restraints" %comp_id)
    else:
        plane_restraints = restraints["_chem_comp_plane_atom"]
        new_plane_restraints = [["pseudo-ring-1", [" C1 ", " C2 ", " C4 ", " C5 "], 0.01],
                                ["pseudo-ring-2", [" C2 ", " C3 ", " C5 ", " O5 "], 0.01],
                                ["pseudo-ring-3", [" C3 ", " C4 ", " O5 ", " C1 "], 0.01]] + \
                                filter_out("pseudo-ring-", plane_restraints) # should be list already
        restraints["_chem_comp_plane_atom"] = new_plane_restraints

        coot.set_monomer_restraints_py(comp_id, restraints)

def add_synthetic_pyranose_planes():
    for comp_id in ["NAG", "BMA", "MAN", "GAL", "GLC", "FUC", "XYP", "BGC"]:
        add_pyranose_pseudo_ring_plane_restraints(comp_id)

def use_unimodal_pyranose_ring_torsions():
    for tlc in ["NAG", "BMA", "MAN", "GAL", "GLC", "FUC", "XYP", "BGC"]:
        coot.use_unimodal_ring_torsion_restraints(tlc)


def multi_add_linked_residue(imol, res_spec, residues_to_add):

    "---------------- multi-add-linked-residue", imol, res_spec
    coot.set_go_to_atom_molecule(imol)
    wm = coot.matrix_state()
    coot.set_matrix(wm/4.)

    current_refinement_rate = coot.dragged_refinement_steps_per_frame()

    current_residue_spec = res_spec

    for residue_to_add in residues_to_add:

        coot_utils.set_go_to_atom_from_res_spec(current_residue_spec)
        coot.set_dragged_refinement_steps_per_frame(300)

        if (not isinstance(current_residue_spec, list)):
            print("OOps not a proper res_spec %s with residue_to_add: %s" %(current_residue_spec, residue_to_add))
            return False
        else:
            if not residue_to_add:
                return "done" # odd return?! why?
            else:
                if (not isinstance(residue_to_add, list) and \
                    not len(residue_to_add) == 2):
                    print("Oops - not a residue_link string pair when adding new_res_pair")
                else:
                    new_res = residue_to_add[0]
                    new_link = residue_to_add[1]

                    new_res_spec = coot.add_linked_residue(imol,
                                                      coot_utils.res_spec_to_chain_id(current_residue_spec),
                                                      coot_utils.res_spec_to_res_no(current_residue_spec),
                                                      coot_utils.res_spec_to_ins_code(current_residue_spec),
                                                      new_res, new_link, 2) # add and fit mode


                    # residues_near_residue takes a 3-part spec and makes 3-part specs
                    if new_res_spec:
                        ls = coot.residues_near_residue_py(imol, current_residue_spec, 1.9)

                        cho_restraints_from_models.add_cho_restraints_for_residue(imol, new_res_spec)

                        with AutoAccept():
                            coot.refine_residues_py(imol, [current_residue_spec] + ls)
                        # make the new res the one to add to (remove starting bool)
                        current_residue_spec = new_res_spec[1:]

    # restore refinement mode and matrix weight
    coot.set_dragged_refinement_steps_per_frame(current_refinement_rate)
    coot.set_matrix(wm)

# return the new molecule number
#
def new_molecule_from_this_glyco_tree():
    with coot_utils.UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
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

def add_linked_residue_add_cho_function(imol, parent, res_pair):

    global add_linked_residue_tree_correlation_cut_off

    def well_fitting_qm(res_spec):
        with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                   aa_ins_code, aa_atom_name, aa_alt_conf]:
            neighbs = coot.residues_near_residue_py(aa_imol, res_spec, 4)
            c = coot.map_to_model_correlation(imol, [res_spec], neighbs, 0,
                                         coot.imol_refinement_map())
            print("######## new residue %s density correlation: %s" %(res_spec, c))
            if (not coot_utils.isNumber(c)):
                return False
            else:
                if c > add_linked_residue_tree_correlation_cut_off:
                    symm_clash = coot.clashes_with_symmetry(imol,
                                                       res_spec_utils.residue_spec_to_chain_id(res_spec),
                                                       res_spec_utils.residue_spec_to_res_no(res_spec),
                                                       coot_utils.residue_spec_to_ins_code(res_spec),
                                                       2.0)
                    if (symm_clash == 1):
                        return False
                    else:
                        return True
                else:
                    return False

    def centre_view_on_residue_centre(res_spec):
        res_centre = residue_centre(imol,
                                    res_spec_utils.residue_spec_to_chain_id(res_spec),
                                    res_spec_utils.residue_spec_to_res_no(res_spec),
                                    coot_utils.residue_spec_to_ins_code(res_spec))
        if (isinstance(res_centre, list)):
            coot.set_rotation_centre(*res_centre)

    # main line
    #
    if (not isinstance(parent, list)):
        print("WARNING:: OOps not a proper res_spec %s with residue_to_add: %s" %(parent, res_pair))
        return False
    else:
        if not (len(res_pair) == 2):
            print("Oops - not a residue-link string pair when adding res-pair", res_pair)
            return False
        else:
            # OK! Go!
            new_link = res_pair[0]
            new_res_type = res_pair[1]

            centre_view_on_residue_centre(parent)

            tree_residues = glyco_tree_residues(imol, parent)
            imol_save = new_molecule_by_residue_specs(imol,
                                                      tree_residues)

            new_res_spec = coot.add_linked_residue(imol,
                                              coot_utils.res_spec_to_chain_id(parent),
                                              coot_utils.res_spec_to_res_no(parent),
                                              coot_utils.res_spec_to_ins_code(parent),
                                              new_res_type, new_link, 2)
            # 2 = add and link mode
            coot.set_mol_displayed(imol_save, 0)
            coot.set_mol_active(imol_save, 0)
            ls = coot.residues_near_residue_py(imol, parent, 1.9)
            local_ls = [parent] + ls
            cho_restraints_from_models.add_cho_restraints_for_residue(imol, new_res_spec)
            coot.rotate_y_scene(100, 0.5)
            with AutoAccept():
                coot.refine_residues_py(imol, local_ls)
            if (not isinstance(new_res_spec, list)):
                # badness
                return False
            else:
                # okay!?
                preped_new_res_spec = new_res_spec[1:]  # strip off leading result
                if well_fitting_qm(preped_new_res_spec):
                    return preped_new_res_spec
                else:
                    # ------------ bad fit -----------------
                    # delete residue and restore others
                    print("----------- That was not well fitting. Deleting:", preped_new_res_spec)
                    delete_extra_restraints_for_residue_spec(imol,
                                                             preped_new_res_spec)
                    coot_utils.delete_residue_by_spec(*preped_new_res_spec)
                    # restore glyco-tree residues from imol_save
                    coot.replace_fragment(imol, imol_save, "//")
                    return False

def add_linked_residue_tree(imol, parent, tree):

    def func_test():
        count = 10
        pass # needed?

    def process_tree(parent, tree, proc_func):

        # helper function to test for a link, resname pair
        #
        def link_res_pair(pair):
            # test if pair i.e. two elements
            if len(pair) == 2:
                # test if both are strings
                res = [isinstance(x, str) for x in pair]
                # could add tests for len pair[1] == 3 and "-" in pair[0]
                return all(res)
            return False

        if not tree:
            return []
        if isinstance(tree[0], list) and not link_res_pair(tree[0]):
            part_1 = process_tree(parent, tree[0], proc_func)
            part_2 = process_tree(parent, tree[1:], proc_func)
            return [part_1, part_2]
        else:
            new_res = proc_func(imol, parent, tree[0])
            return [new_res, process_tree(new_res, tree[1:], proc_func)]

    def is_just_an_ASN_qm(imol, glyco_tree):

        print("glyco-tree:", glyco_tree)
        if not isinstance(glyco_tree, list):
            return False
        else:
            if not (len(glyco_tree) == 1):
                return False
            else:
                with coot_utils.UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no,
                                               aa_ins_code, aa_atom_name, aa_alt_conf, aa_res_spec]:
                    res_spec = aa_res_spec
                    rn = coot.residue_name(imol,
                                      res_spec_utils.residue_spec_to_chain_id(res_spec),
                                      res_spec_utils.residue_spec_to_res_no(res_spec),
                                      coot_utils.residue_spec_to_ins_code(res_spec))
                    if not isinstance(rn, str):
                        return False
                    else:
                        return rn == "ASN"

    # main line of add_linked_residue_tree
    #
    add_synthetic_pyranose_planes()
    use_unimodal_pyranose_ring_torsions()
    coot.set_refine_with_torsion_restraints(1)
    wm = coot.matrix_state()
    coot.set_matrix(wm / 4.)
    coot.set_residue_selection_flash_frames_number(1)
    coot.set_go_to_atom_molecule(imol)
    coot_utils.set_go_to_atom_from_res_spec(parent)
    previous_m = coot.default_new_atoms_b_factor()
    m = coot.median_temperature_factor(imol)
    try:
        new_m = m * 1.55
        if new_m > 0:
            coot.set_default_temperature_factor_for_new_atoms(new_m)
    except:
        print("BL INFO:: not changing default B for new carb atoms")

    # prevent tree building if there is already a partial tree here
    # (only proceed with the one ASN)
    #
    with coot_utils.UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
                                   aa_atom_name, aa_alt_conf, aa_res_spec]:
        start_tree = glyco_tree_residues(aa_imol, aa_res_spec)
        # print "::::::::::::::::: start-tree:", start_tree

        # why do we need both test here? 
        # 5n11 needs the is-just-an-ASN? test
        # 5n09/5wzy needs the null test. 
        # Hmmm.
        if not ((start_tree == []) or (is_just_an_ASN_qm(aa_imol, start_tree))):
            coot.info_dialog("Must start on Single ASN")
            print("start_tree:", start_tree)
        else:
            # OK, continue
            start_pos_view = coot_utils.add_view_here("Glyco Tree Start Pos")
            process_tree(parent, tree, add_linked_residue_add_cho_function)
            coot.go_to_view_number(start_pos_view, 0)
            with AutoAccept():
                coot.refine_residues_py(aa_imol, glyco_tree_residues(aa_imol, aa_res_spec))

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
    coot.set_default_temperature_factor_for_new_atoms(previous_m)
    coot.set_matrix(wm)


def add_linked_residue_with_extra_restraints_to_active_residue(new_res_type,
                                                               link_type):
    wm = coot.matrix_state()
    coot.set_matrix(wm / 8.)
    coot.set_refine_with_torsion_restraints(1)
    coot.set_add_linked_residue_do_fit_and_refine(0)
    with coot_utils.UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
                                   aa_atom_name, aa_alt_conf, aa_res_spec]:
        new_res_spec = coot.add_linked_residue(aa_imol, aa_chain_id, aa_res_no,
                                          aa_ins_code, new_res_type, link_type, 2)
        if isinstance(new_res_spec, list):
            cho_restraints_from_models.add_cho_restraints_for_residue(aa_imol, new_res_spec)
            # refine that
            with AutoAccept():
                residues = [aa_res_spec] + coot.residues_near_residue_py(aa_imol, aa_res_spec, 1.9)
                coot.refine_residues_py(aa_imol, residues)
    coot.set_matrix(wm)

def delete_all_cho():
    delete_cho_ls = []
    with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
                               aa_atom_name, aa_alt_conf]:
        if coot_utils.valid_model_molecule_qm(aa_imol):
            with coot_utils.NoBackups(aa_imol):
                for chain_id in coot_utils.chain_ids(aa_imol):
                    for res_serial in range(coot.chain_n_residues(chain_id, aa_imol)):
                        res_no = coot.seqnum_from_serial_number(aa_imol,
                                                           chain_id, res_serial)
                        rn = coot.residue_name(aa_imol, chain_id, res_no, "")
                        if (isinstance(rn, str)):
                            if rn in ["NAG", "MAN", "BMA", "FUC", "XYP",
                                      "SIA", "GLC", "GAL"]:
                                residue_spec = [chain_id, res_no, ""]
                                delete_cho_ls.append(residue_spec)
#                now we have delete_residues, we don't need to delete them one by one
#                for cho_res_spec in delete_cho_ls:
#                    coot.delete_residue(aa_imol,
#                                   res_spec_utils.residue_spec_to_chain_id(cho_res_spec),
#                                   res_spec_utils.residue_spec_to_res_no(cho_res_spec), "")
                coot.delete_residues_py(aa_imol, delete_cho_ls)


