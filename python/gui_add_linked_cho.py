
def add_pyranose_pseudo_ring_plane_restraints(comp_id):

    def filter_out(plane_name_sub_string, plane_restraints):
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
                                filter_out("pseudo-ring-3", plane_restraints) #should be list already
        restraints["_chem_comp_plane_atom"] = new_plane_restraints
        print "BL DEBUG:: new plane restraints", new_plane_restraints

        set_monomer_restraints(comp_id, restraints)
        
def add_synthetic_carbohydrate_planes():
    add_pyranose_pseudo_ring_plane_restraints("NAG")
    add_pyranose_pseudo_ring_plane_restraints("BMA")
    add_pyranose_pseudo_ring_plane_restraints("MAN")
    add_pyranose_pseudo_ring_plane_restraints("GAL")
    add_pyranose_pseudo_ring_plane_restraints("GLC")
    add_pyranose_pseudo_ring_plane_restraints("XYP")
    add_pyranose_pseudo_ring_plane_restraints("FUC")


def multi_add_linked_residue(imol, res_spec, residues_to_add):

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
                                                      new_res, new_link, 3)


                    # residues_near_residue takes a 3-part spec and makes 3-part specs
                    ls = residues_near_residue(imol, current_residue_spec, 1.9)
                    with AutoAccept():
                        refine_residues(imol, [current_residue_spec] + ls)
                    # make the new res the one to add to (remove starting bool)
                    current_residue_spec = new_res_spec[1:]  

    # restore refinement mode and matrix weight
    set_dragged_refinement_steps_per_frame(current_refinement_rate)
    set_matrix(current_matrix)


def oligomannose_tree():
    ret = [["NAG-ASN", "NAG"],
           ["BETA1-4", "NAG"],
           ["BETA1-4", "BMA"],
           [
               ["ALPHA1-6", "MAN"],
               [
                   ["ALPHA1-6", "MAN"],
                   ["ALPHA1-2", "MAN"]
                   ],
               [
                   ["ALPHA1-3", "MAN"],
                   ["ALPHA1-2", "MAN"]
                   ]
               ],
           [
               ["ALPHA1-3", "MAN"],
               ["ALPHA1-2", "MAN"],
               ["ALPHA1-2", "MAN"]
               ]
           ]
    return ret

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

            active_atom = active_residue_py()
            active_residue = active_atom[:4]
            save_residue_specs = glyco_tree_residues_py(imol, active_residue)

            imol_glyco_pre = new_molecule_by_residue_specs_py(imol, save_residue_specs)
            set_mol_displayed(imol_glyco_pre, 0)
            set_mol_active(imol_glyco_pre, 0)
            new_res_spec = add_linked_residue(imol,
                                              res_spec2chain_id(parent),
                                              res_spec2res_no(parent),
                                              res_spec2ins_code(parent),
                                              new_res_type, new_link, 3)
            ls = residues_near_residue(imol, parent, 1.9)
            local_ls = [parent] + ls
            with AutoAccept():
                refine_residues(imol, local_ls)
            # make this optional
            rotate_y_scene(200, 0.5)
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
                    print "----------- That was not well fitting. Deleting:", preped_new_res_spec
                    delete_residue_by_spec(imol, preped_new_res_spec)
                    replace_residues_from_mol_py(imol, imol_glyco_pre, save_residue_specs)
                    with AutoAccept():
                        # Note: may not get rid of screwed up refinement from
                        # adding a carb too much...
                        refine_residues(imol, local_ls)
                        return False
                    
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

    # main line of add_linked_residue_tree
    #
    add_synthetic_carbohydrate_planes()
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
    process_tree(parent, tree, func)
    g = glyco_validate()
    g.auto_delete_residues()
    # reset
    set_default_temperature_factor_for_new_atoms(previous_m)


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

            def add_carbo_link_func(sugar, link):
                current_weight = matrix_state()
                set_matrix(8)
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    add_linked_residue(aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
                                       sugar, link, 3)
                    set_matrix(current_weight)
                    
            add_simple_coot_menu_menuitem(
                menu, "Add a ASN-NAG NAG",
                lambda func: add_carbo_link_func("NAG", "NAG-ASN"))

            add_simple_coot_menu_menuitem(
                menu, "Add a BETA1-4 NAG",
                lambda func: add_carbo_link_func("NAG", "BETA1-4"))

            add_simple_coot_menu_menuitem(
                menu, "Add a BETA1-4 BMA",
                lambda func: add_carbo_link_func("BMA", "BETA1-4"))
            
            add_simple_coot_menu_menuitem(
                menu, "Add an ALPHA1-2 MAN",
                lambda func: add_carbo_link_func("MAN", "ALPHA1-2"))
            
            add_simple_coot_menu_menuitem(
                menu, "Add an ALPHA1-3 MAN",
                lambda func: add_carbo_link_func("MAN", "ALPHA1-3"))

            # we should do this only if we are sitting on an SIA.
            # Attaching a SIA to a MAN (i.e. reverse order) would be a
            # good test too...
            add_simple_coot_menu_menuitem(
                menu, "Add an ALPHA2-3 MAN",
                lambda func: add_carbo_link_func("MAN", "ALPHA2-3"))

            # same consideration as above
            add_simple_coot_menu_menuitem(
                menu, "Add an ALPHA2-3 GAL",
                lambda func: add_carbo_link_func("GAL", "ALPHA2-3"))

            add_simple_coot_menu_menuitem(
                menu, "Add an ALPHA1-6 MAN",
                lambda func: add_carbo_link_func("MAN", "ALPHA1-6"))

            add_simple_coot_menu_menuitem(
                menu, "Add an ALPHA1-3 FUC",
                lambda func: add_carbo_link_func("FUC", "ALPHA1-3"))

            add_simple_coot_menu_menuitem(
                menu, "Add an ALPHA1-6 FUC",
                lambda func: add_carbo_link_func("FUC", "ALPHA1-6"))

            add_simple_coot_menu_menuitem(
                menu, "Add an XYP-BMA XYP",
                lambda func: add_carbo_link_func("XYP", "XYP-BMA"))

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

            add_simple_coot_menu_menuitem(
                menu, "Auto Fit & Refine On for Link Addition",
                lambda func: set_add_linked_residue_do_fit_and_refine(1))

            add_simple_coot_menu_menuitem(
                menu, "Auto Fit & Refine Off for Link Addition",
                lambda func: set_add_linked_residue_do_fit_and_refine(0))

            def add_oligo_mannose_func():
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    make_backup(aa_imol)
                    # switch backup off?!
                    add_linked_residue_tree(aa_imol,
                                            [aa_chain_id, aa_res_no, aa_ins_code],
                                            oligomannose_tree())
                
            def add_paucimannose_func():
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    make_backup(aa_imol)
                    # switch backup off?!
                    add_linked_residue_tree(aa_imol,
                                            [aa_chain_id, aa_res_no, aa_ins_code],
                                            paucimannose_tree())

            add_simple_coot_menu_menuitem(
                menu, "Add Oligomannose",
                lambda func: add_oligo_mannose_func())

            add_simple_coot_menu_menuitem(
                menu, "Add Paucimannose",
                lambda func: add_paucimannose_func())

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

            add_simple_coot_menu_menuitem(
                menu, "Torsion Fit This Residue and Neighbours",
                lambda func: torsion_fit_this_and_neighbours_func())


            add_simple_coot_menu_menuitem(
                menu, "Torsion Fit & Refine this residue",
                lambda func: torsion_fit_this_func(True))

            add_simple_coot_menu_menuitem(
                menu, "Add synthetic pyranose plane restraints",
                lambda func: add_synthetic_carbohydrate_planes())

