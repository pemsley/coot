
import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk
import coot
import coot_gui
import coot_gui_api
import coot_utils
import add_linked_cho

def interactive_add_cho_dialog():

    def dummy_func():
        print("dummy_func()")
        pass

    def refine_tree_func():
        with coot_utils.UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
                                       aa_atom_name, aa_alt_conf, aa_res_spec]:
            refine_residues(aa_imol, glyco_tree_residues(aa_imol, aa_res_spec))

    add_linked_cho.add_synthetic_pyranose_planes()
    add_linked_cho.use_unimodal_pyranose_ring_torsions()
    # button list with [label, function]
    buttons = [
        ["Update for Current Residue", lambda func: dummy_func()], # dummy_func replaced later
        ["Refine Tree", lambda func: refine_tree_func()],
        ["Add a NAG-ASN NAG",
         lambda func:
         add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("NAG", "NAG-ASN")],
        ["Add a BETA1-4 NAG",
         lambda func:
         add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("NAG", "BETA1-4")],
        ["Add a BETA1-4 BMA",
         lambda func:
         add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("BMA", "BETA1-4")],
        ["Add an ALPHA1-2 MAN",
         lambda func:
         add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("MAN", "ALPHA1-2")],
        ["Add an ALPHA1-3 MAN",
         lambda func:
         add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("MAN", "ALPHA1-3")],
        ["Add an ALPHA2-3 MAN",
         lambda func:
         add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("MAN", "ALPHA2-3")],
        ["Add an ALPHA2-3 GAL",
         lambda func:
         add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("GAL", "ALPHA2-3")],
        ["Add an ALPHA1-6 MAN",
         lambda func:
         add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("MAN", "ALPHA1-6")],
        ["Add a BETA1-2 NAG",
         lambda func:
         add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("NAG", "BETA1-2")],
        ["Add a BETA1-4 GAL",
         lambda func:
         add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("GAL", "BETA1-4")],
        ["Add an ALPHA1-2 FUC",
         lambda func:
         add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("FUC", "ALPHA1-2")],
        ["Add an ALPHA1-3 FUC",
         lambda func:
         add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("FUC", "ALPHA1-3")],
        ["Add an ALPHA1-6 FUC",
         lambda func:
         add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("FUC", "ALPHA1-6")],
        ["Add an BETA1-6 FUL",
         lambda func:
         add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("FUL", "BETA1-6")],
        ["Add an XYP-BMA XYP",
         lambda func:
         add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("XYP", "XYP-BMA")],
        ["Add an ALPHA2-3 SIA",
         lambda func:
         add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("SIA", "ALPHA2-3")],
        ["Add an ALPHA2-6 SIA",
         lambda func:
         add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("SIA", "ALPHA2-6")]
    ]

    vbox = coot_gui.dialog_box_of_buttons("Add N-linked Glycan",
                                 [420, 600], buttons, "Close")[0]
    gui_add_linked_cho_dialog_vbox_set_rotation_centre_hook(vbox)
    # set the callback on the first button
    children = vbox.get_children()
    if isinstance(children, list):
        if len(children) > 0:
            first_button = children[0]
            first_button.connect("clicked", lambda func:
                                 gui_add_linked_cho_dialog_vbox_set_rotation_centre_hook(vbox))

    # add a widget to allow the user to choose the tree type
    table = Gtk.Table(3, 2, False)
    butt_1 = Gtk.RadioButton(None, "High Mannose")
    butt_2 = Gtk.RadioButton(butt_1, "Hybrid (Mammal)")
#    butt_3 = Gtk.RadioButton(butt_1, "Hybrid (Plant)")
    butt_4 = Gtk.RadioButton(butt_1, "Complex (Mammal)")
    butt_5 = Gtk.RadioButton(butt_1, "Complex (Plant)")
    butt_6 = Gtk.RadioButton(butt_1, "Expert User Mode")

    butt_1.show()
    butt_2.show()
 #   butt_3.show()
    butt_4.show()
    butt_5.show()
    butt_6.show()

    # this is now how we do it these days. FIXME-PE
    #     # add buttons for nice(?) layout/order
    #     table.attach(butt_1, 0, 1, 0, 1, Gtk.EXPAND|Gtk.FILL, Gtk.EXPAND|Gtk.FILL, 0, 0) # high mannose
    #     table.attach(butt_4, 1, 2, 0, 1, Gtk.EXPAND|Gtk.FILL, Gtk.EXPAND|Gtk.FILL, 0, 0) # complex mammal
    #     table.attach(butt_6, 0, 1, 1, 2, Gtk.EXPAND|Gtk.FILL, Gtk.EXPAND|Gtk.FILL, 0, 0) # Expert
    # #    table.attach(butt_3, 1, 2, 1, 2, Gtk.EXPAND|Gtk.FILL, Gtk.EXPAND|Gtk.FILL, 0, 0)
    #     table.attach(butt_5, 1, 2, 1, 2, Gtk.EXPAND|Gtk.FILL, Gtk.EXPAND|Gtk.FILL, 0, 0) # complex plant
    #     table.attach(butt_2, 2, 3, 0, 1, Gtk.EXPAND|Gtk.FILL, Gtk.EXPAND|Gtk.FILL, 0, 0) # hybrid mammal
    
    vbox.pack_start(table, True, True, 2)
    table.show()
    vbox.reorder_child(table, 0)

    for butt in [butt_1, butt_2, butt_4, butt_5, butt_6]:
        butt.connect("toggled", lambda func:
                     gui_add_linked_cho_dialog_vbox_set_rotation_centre_hook(vbox))

    # "global" var post-set-rotation-centre-hook
    # BL Note:: maybe should be a global!?
    global post_set_rotation_centre_script
    def post_set_rotation_centre_script():
        gui_add_linked_cho_dialog_vbox_set_rotation_centre_hook(vbox)

#
def glyco_tree_dialog_set_button_active_state(button, glyco_id, tree_type):

    def glyco_id2level_number(glyco_id):
        return glyco_id[0]

    def glyco_id2prime_arm_sym(glyco_id):
        return glyco_id[1]

    def glyco_id2residue_type(glyco_id):
        return glyco_id[2]

    def glyco_id2link_type(glyco_id):
        return glyco_id[3]

    def glyco_id2paren_residue_type(glyco_id):
        return glyco_id[4]

    def glyco_id2residue_spec(glyco_id):
        return glyco_id[5]
    
    def get_sensitive_button_list(glyco_id, tree_type):

        if not isinstance(glyco_id, list):
            return []
        else:
            level_number = glyco_id2level_number(glyco_id)
            prim_arm_sym = glyco_id2prime_arm_sym(glyco_id)
            residue_type = glyco_id2residue_type(glyco_id)
            link_type = glyco_id2link_type(glyco_id)
            parent_residue_type = glyco_id2paren_residue_type(glyco_id)
            residue_spec = glyco_id2residue_spec(glyco_id)
            active_button_label_ls = []
            
            if tree_type == 'expert-user-mode':
                active_button_label_ls = "expert-user-mode" # ???

            # -----------------------------------------------------------------------------------
            # Note the trees tested here match those from (get-tree-type) which examines the button
            # label of the active radio button in the dialog (these are not the auto-build trees)
            # ------------------------------------------------------------------------------------
            # BL says:: test and control the trees

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

            # hybrid mammal
            if tree_type == 'hybrid-mammal':

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

                if level_number == 2:
                    if residue_type == "FUC":
                        active_button_label_ls = ["Add a BETA1-4 GAL"]

                if level_number == 3:
                    if residue_type == "BMA":
                        active_button_label_ls = ["Add an ALPHA1-3 MAN",
                                                  "Add an ALPHA1-6 MAN"]

                if level_number == 3:
                    if residue_type == "GAL":
                        active_button_label_ls = ["Add an ALPHA1-2 FUC"]

                if level_number == 4:
                    if residue_type == "MAN":
                        active_button_label_ls = ["Add an ALPHA1-3 MAN",
                                                  "Add an ALPHA1-6 MAN",
                                                  "Add a BETA1-2 NAG"]
                
                if level_number == 5:
                    if residue_type == "NAG":
                        active_button_label_ls = ["Add a BETA1-4 GAL"]
                
                if level_number == 5:
                    if residue_type == "GAL":
                        active_button_label_ls = ["Add an ALPHA2-3 SIA",
                                                  "Add an ALPHA2-6 SIA"]
                
            # hybrid plant
            if tree_type == 'hybrid-plant':

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

                if level_number == 2:
                    if residue_type == "FUC":
                        active_button_label_ls = ["Add a BETA1-4 GAL"]

                if level_number == 3:
                    if residue_type == "BMA":
                        active_button_label_ls = ["Add an ALPHA1-3 MAN",
                                                  "Add an ALPHA1-6 MAN",
                                                  "Add an XYP-BMA XYP"]

                if level_number == 3:
                    if residue_type == "GAL":
                        active_button_label_ls = ["Add an ALPHA1-2 FUC"]

                # note that level-number is not enough for complete disambiguation,
                # we need to know we are 4 or 4' (Vliegenthart et al 1983 nomenclature).
                if level_number == 4:
                    if residue_type == "MAN":
                        active_button_label_ls = ["Add an ALPHA1-3 MAN",
                                                  "Add an ALPHA1-6 MAN",
                                                  "Add a BETA1-2 NAG"]
                
                if level_number == 5:
                    if residue_type == "NAG":
                        active_button_label_ls = ["Add a BETA1-4 GAL"]
                
                if level_number == 5:
                    if residue_type == "GAL":
                        active_button_label_ls = ["Add an ALPHA2-3 SIA",
                                                  "Add an ALPHA2-6 SIA"]
                
            # complex mammal
            #
            if tree_type == 'complex-mammal':

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
                                                  "Add a BETA1-4 NAG"]

                #  note that level-number is not enough for complete disambiguation, 
                # we need to know we are 4 or 4' (Vl...? nomenclature).
                #
                if level_number == 4:
                    if residue_type == "MAN":
                        active_button_label_ls = ["Add an ALPHA1-3 MAN",
                                                  "Add an ALPHA1-6 MAN",
                                                  "Add a BETA1-2 NAG"]

                if level_number == 5:
                    if residue_type == "NAG":
                        active_button_label_ls = ["Add a BETA1-4 GAL"]

                if level_number == 6:
                    if residue_type == "GAL":
                        active_button_label_ls = ["Add an ALPHA2-3 SIA",
                                                  "Add an ALPHA2-6 SIA"]

                        
            # complex plant
            #
            if tree_type == 'complex-plant':

                if level_number == 0:
                    if residue_type == "ASN":
                        active_button_label_ls = ["Add a NAG-ASN NAG"]
                        
                if level_number == 1:
                    if residue_type == "NAG":
                        active_button_label_ls = ["Add a BETA1-4 NAG",
                                                  "Add an ALPHA1-3 FUC",
                                                  "Add an ALPHA1-6 FUC"]
                        
                if level_number == 2:
                    if residue_type == "NAG":
                        active_button_label_ls = ["Add a BETA1-4 BMA"]

                if level_number == 3:
                    if residue_type == "BMA":
                        active_button_label_ls = ["Add an ALPHA1-3 MAN",
                                                  "Add an ALPHA1-6 MAN",
                                                  "Add an XYP-BMA XYP",
                                                  "Add a BETA1-4 NAG"]

                if level_number == 4:
                    if residue_type == "MAN":
                        active_button_label_ls = ["Add an ALPHA1-3 MAN",
                                                  "Add an ALPHA1-6 MAN",
                                                  "Add a BETA1-2 NAG"]

                if level_number == 5:
                    if residue_type == "NAG":
                        active_button_label_ls = ["Add a BETA1-4 GAL"]

                if level_number == 6:
                    if residue_type == "GAL":
                        active_button_label_ls = ["Add an ALPHA2-3 SIA",
                                                  "Add an ALPHA2-6 SIA"]
                        
            return active_button_label_ls

    # main line
    #
    l = button.get_label()
    active_button_label_ls = get_sensitive_button_list(glyco_id, tree_type)
    if (active_button_label_ls == "expert-user-mode"):
        button.set_sensitive(True)
    else:
        if not (l == "Update for Current Residue") and not (l == "Refine Tree"):
            button.set_sensitive(l in active_button_label_ls)

# vbox is the vbox of the dialog box of buttons. One of the children of the vbox
# is the table that contains the buttons
#
def gui_add_linked_cho_dialog_vbox_set_rotation_centre_hook(vbox):

    def get_tree_type():
        tree_type = "oligomannose"
        children = vbox.get_children()
        for child in children:
            if type(child) == Gtk.Table:
                for table_child in child:
                    if type(table_child) == Gtk.RadioButton:
                        if table_child.get_active():
                            l = table_child.get_label()
                            # this is a bit ugly because we are testing that these strings
                            # match button labels (set in interactive-add-cho-dialog)
                            if l == "High Mannose":
                                tree_type = 'oligomannose'
                            if l == "Hybrid (Mammal)":
                                tree_type = 'hybrid-mammal'
                            if l == "Hybrid (Plant)":
                                tree_type = 'hybrid-plant'
                            if l == "Expert User Mode":
                                tree_type = 'expert-user-mode'
                            if l == "Complex (Mammal)":
                                tree_type = 'complex-mammal'
                            if l == "Complex (Plant)":
                                tree_type = 'complex-plant'
        return tree_type

    with coot_utils.UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
                                   aa_atom_name, aa_alt_conf, aa_res_spec]:
        glyco_id = coot.glyco_tree_residue_id_py(aa_imol, aa_res_spec)
        # Paule says:
        # if it was an ASP create a level-0 glyco-id for that (glyco-tree-residue-id doesn't
        # do that (not sure why)).
        if not glyco_id:
            rn = coot.residue_name(aa_imol, aa_chain_id, aa_res_no, aa_ins_code)
            if isinstance(rn, str):
                if rn == "ASN":
                    glyco_id = [0, "unset", "ASN", "", "", aa_res_spec]
        if isinstance(glyco_id, list):
            tree_type = get_tree_type()
            children = vbox.get_children()
            for child in children:
                if (type(child) == Gtk.Button):
                    glyco_tree_dialog_set_button_active_state(child, glyco_id,
                                                              tree_type)
            return True
        return False


# There is some graphics to the validation here, so keep it in the graphics part
# wont be missed in the non-grpahical part for now, I think (BL that is).
    
# this has to be at the top level for reasons that I don't understand
# (python function scoping) - I'll just put it down to python being shit as usual.
#
def glyco_validation_dialog_set_go_to_residue(imol, residue_spec):
    rc = residue_centre(imol,
                        res_spec_utils.residue_spec_to_chain_id(residue_spec),
                        res_spec_utils.residue_spec_to_res_no(residue_spec),
                        '')
    coot.set_rotation_centre(*rc)

def load_privateer_dictionary():
    if os.path.exists("privateer-lib.cif"):
        coot.read_cif_dictionary("privateer-lib.cif")
        coot.set_refine_with_torsion_restraints(1)
    
class glyco_validate:

    def run_privateer(self, imol, glyco_tree_residues, hklin_fn, fp_col, sigfp_col, pdbin_fn, privateer_log):
        if coot.is_valid_model_molecule(imol):
            args = ['-mtzin', hklin_fn, '-colin-fo', fp_col+','+sigfp_col,
                    '-pdbin', pdbin_fn]
            # using'-mode', 'ccp4i2' makes a file with no residue nambers (afaics)
            coot_utils.popen_command("privateer", args, [], privateer_log, False)

    def make_privateer_validation_info(self, imol, fp_col, sigfp_col, glyco_tree_residues):

        imol_map = coot.imol_refinement_map()

        if len(glyco_tree_residues) > 0:
            d = coot_utils.get_directory("coot-ccp4")
            spid = str(os.getpid())
            fn_pdb = "coot-privateer-" + spid + ".pdb"
            fn_log = "coot-privateer-" + spid + ".log"
            privateer_pdb = os.path.join(d, fn_pdb)
            privateer_log = os.path.join(d, fn_log)
            hklin_fn = coot.mtz_file_name(imol_map);
            coot.write_pdb_file(imol, privateer_pdb)
            self.run_privateer(imol, glyco_tree_residues, hklin_fn, fp_col, sigfp_col, privateer_pdb, privateer_log)
            pvi = self.parse_privateer_log(privateer_log, imol, glyco_tree_residues)
            return pvi
        else:
            return []

    # return a list of residues with privateer validation info
    #
    def parse_privateer_log(self, log_file_name, imol, glyco_tree_residues):

        print('parse_privateer_log', log_file_name, imol, glyco_tree_residues)

        pvi = []
        f = open(log_file_name)
        lines = f.readlines()
        f.close()
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
                                         '-' + str(coot_utils.residue_spec_to_res_no(r))
                                # print "res_id", res_id
                                if words[1] == res_id:
                                    # print words[12] , yes or check
                                    new_item = (r, words)
                                    pvi.append(new_item)
                            except TypeError as e:
                                print(e)
        print("parsed", log_file_name)
        return pvi

    def make_validation_dialog(self, imol, privateer_validation_info):

        print("make_validation_dialog", privateer_validation_info)
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
            coot_gui.dialog_box_of_buttons("Privateer Validation", (400, 220), buttons, " Close ")

    def validation_dialog(self):
        
        active_atom = coot.active_residue_py()
        try:
            imol = active_atom[0]
            active_residue = active_atom[:4]
            glyco_tree_residues = coot.glyco_tree_residues_py(imol, active_residue)

            print('imol', imol)
            print('active_residue', active_residue)
            print('glyco_tree_residues', glyco_tree_residues)

            fp_col='FP'
            sigfp_col='SIGFP'    
            pvi = self.make_privateer_validation_info(imol, fp_col, sigfp_col, glyco_tree_residues)
            # print "validation_dialog()", pvi
            if len(pvi) > 0:
                # now make a gui
                self.make_validation_dialog(imol, pvi)

        except KeyError as e:
            print(e)

        # no active atom
        except TypeError as e:
            print(e)

    def auto_delete_residues_internal(self, imol, glyco_tree_residues):
        fp_col='FP'
        sigfp_col='SIGFP'
        pvi = self.make_privateer_validation_info(imol, fp_col, sigfp_col, glyco_tree_residues)
        for res_info in pvi:
            print("test", res_info)
            res_status = res_info[1][12]
            if res_status == 'check' or res_status == 'no':
                coot_utils.delete_residue_by_spec(imol, res_info[0])

    def auto_delete_residues(self):
        try:
            active_atom = coot.active_residue_py()
            imol = active_atom[0]
            active_residue = active_atom[:4]
            glyco_tree_residues = coot.glyco_tree_residues_py(imol, active_residue)
            self.auto_delete_residues_internal(imol, glyco_tree_residues)
        except TypeError as e:
            print(e)
    
            
# graphics...

def add_module_carbohydrate_gui():
    menu = coot_gui.attach_module_menu_button("Glyco")

    coot_gui.add_simple_action_to_menu(
        menu, "N-linked Glycan Addition...","interactive_add_cho_dialog",
        lambda _simple_action, _two:
        interactive_add_cho_dialog())
    
    def add_multi_carbo_link_func(link_list):
        with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                    aa_ins_code, aa_atom_name, aa_alt_conf]:
            add_linked_cho.multi_add_linked_residue(aa_imol,
                                        [aa_chain_id, aa_res_no, aa_ins_code],
                                        link_list)


    def set_default_cho_b_factor_func():
        with coot_utils.UsingActiveAtom(True) as [aa_imol, aa_chain_id, aa_res_no,
                                        aa_ins_code, aa_atom_name,
                                        aa_alt_conf, aa_res_spec]:
            residues = coot.residues_near_residue_py(aa_imol, aa_res_spec, 10)
            imol_region = coot.new_molecule_by_residue_specs(aa_imol, residues)
            # BL says:: why not do a new mol by sphere selection?!
            m = coot.median_temperature_factor(imol_region)
            coot.close_molecule(imol_region)
            if coot_utils.isNumber(m):
                new_m = m* 1.55
                coot.set_default_temperature_factor_for_new_atoms(new_m)
                s = "New Temperature Factor set to " + str(new_m)
                coot.info_dialog(s)
        
    coot_gui.add_simple_action_to_menu(
        menu, "Set Default N-linked CHO Atoms B-factor","set_default_cho_b_factor",
        lambda _simple_action, _two: set_default_cho_b_factor_func()
    )

    
    coot_gui.add_simple_action_to_menu(
        menu, "N-link add NAG, NAG, BMA","multi_carbo_link_NAG_NAG_BMA",
        lambda _simple_action, _two: add_multi_carbo_link_func([["NAG", "NAG-ASN"],
                                                ["NAG", "BETA1-4"],
                                                ["BMA", "BETA1-4"]]))
    
    # coot_gui.add_simple_coot_menu_menuitem(
    #     menu, "Add a ASN-NAG NAG",
    #     lambda _simple_action, _two:
    #     add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("NAG", "NAG-ASN"))

    # coot_gui.add_simple_coot_menu_menuitem(
    #     menu, "Add a BETA1-4 NAG",
    #     lambda _simple_action, _two:
    #     add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("NAG", "BETA1-4"))

    # coot_gui.add_simple_coot_menu_menuitem(
    #     menu, "Add a BETA1-4 BMA",
    #     lambda _simple_action, _two:
    #     add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("BMA", "BETA1-4"))
    
    # coot_gui.add_simple_coot_menu_menuitem(
    #     menu, "Add an ALPHA1-2 MAN",
    #     lambda _simple_action, _two:
    #     add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("MAN", "ALPHA1-2"))
    
    # coot_gui.add_simple_coot_menu_menuitem(
    #     menu, "Add an ALPHA1-3 MAN",
    #     lambda _simple_action, _two:
    #     add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("MAN", "ALPHA1-3"))

    # we should do this only if we are sitting on an SIA.
    # Attaching a SIA to a MAN (i.e. reverse order) would be a
    # good test too...
    # coot_gui.add_simple_coot_menu_menuitem(
    #     menu, "Add an ALPHA2-3 MAN",
    #     lambda _simple_action, _two:
    #     add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("MAN", "ALPHA2-3"))

    # # same consideration as above
    # coot_gui.add_simple_coot_menu_menuitem(
    #     menu, "Add an ALPHA2-3 GAL",
    #     lambda _simple_action, _two:
    #     add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("GAL", "ALPHA2-3"))

    # coot_gui.add_simple_coot_menu_menuitem(
    #     menu, "Add an ALPHA1-6 MAN",
    #     lambda _simple_action, _two:
    #     add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("MAN", "ALPHA1-6"))

    # coot_gui.add_simple_coot_menu_menuitem(
    #     menu, "Add an ALPHA1-3 FUC",
    #     lambda _simple_action, _two:
    #     add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("FUC", "ALPHA1-3"))

    # coot_gui.add_simple_coot_menu_menuitem(
    #     menu, "Add an ALPHA1-6 FUC",
    #     lambda _simple_action, _two:
    #     add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("FUC", "ALPHA1-6"))

    # coot_gui.add_simple_coot_menu_menuitem(
    #     menu, "Add an XYP-BMA XYP",
    #     lambda _simple_action, _two:
    #     add_linked_cho.add_linked_residue_with_extra_restraints_to_active_residue("XYP", "XYP-BMA"))


    # the mode in the function call now takes take of this
    # coot_gui.add_simple_coot_menu_menuitem(
    #     menu, "Auto Fit & Refine On for Link Addition",
    #     lambda _simple_action, _two: coot.set_add_linked_residue_do_fit_and_refine(1))

    # coot_gui.add_simple_coot_menu_menuitem(
    #     menu, "Auto Fit & Refine Off for Link Addition",
    #     lambda _simple_action, _two: coot.set_add_linked_residue_do_fit_and_refine(0))

    def add_oligo_tree_func(oligo_tree):
        with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                    aa_ins_code, aa_atom_name, aa_alt_conf]:
            coot.make_backup(aa_imol)
            # switch backup off?!
            add_linked_cho.add_linked_residue_tree(aa_imol,
                                    [aa_chain_id, aa_res_no, aa_ins_code],
                                    oligo_tree)
        
    coot_gui.add_simple_action_to_menu(
        menu, "Add High Mannose","add_high_mannose",
        lambda _simple_action, _two: add_oligo_tree_func(add_linked_cho.oligomannose_tree()))

    coot_gui.add_simple_action_to_menu(
        menu, "Add Hybrid (Mammal)","add_hybrid_mammal",
        lambda _simple_action, _two: add_oligo_tree_func(add_linked_cho.hybrid_mammal_tree()))
    
#            in practice, no one will be doing this. 
#            coot_gui.add_simple_action_to_menu(
#                menu, "Add Hybrid (Plant)",
#                lambda _simple_action, _two: add_oligo_tree_func(add_linked_cho.hybrid_plant_derived_tree()))

    coot_gui.add_simple_action_to_menu(
        menu, "Add Complex (Mammal)","add_complex_mammal",
        lambda _simple_action, _two: add_oligo_tree_func(add_linked_cho.complex_mammal_tree()))

    coot_gui.add_simple_action_to_menu(
        menu, "Add Complex (Plant)","add_complex_plant",
        lambda _simple_action, _two: add_oligo_tree_func(add_linked_cho.complex_plant_tree()))

    coot_gui.add_simple_action_to_menu(
        menu, "Delete All Carbohydrate","delete_all_cho",
        lambda _simple_action, _two: add_linked_cho.delete_all_cho())

    def torsion_fit_this_func(refine = False):
        with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                    aa_ins_code, aa_atom_name, aa_alt_conf]:
            centre_residue = [aa_chain_id,aa_res_no, aa_ins_code]
            coot.multi_residue_torsion_fit(aa_imol,
                                        [centre_residue],
                                        30000)
            if refine:
                with AutoAccept():
                    refine_residues(aa_imol, [centre_residue])

    def torsion_fit_this_and_neighbours_func(refine = False):
        with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                    aa_ins_code, aa_atom_name, aa_alt_conf]:
            centre_residue = [aa_chain_id,aa_res_no, aa_ins_code]
            residues = coot.residues_near_residue_py(aa_imol, centre_residue, 1.9)
            residues.append(centre_residue)
            coot.multi_residue_torsion_fit(aa_imol, residues, 30000)
            if refine:
                with AutoAccept():
                    refine_residues(aa_imol, [centre_residue])

    coot_gui.add_simple_action_to_menu(
        menu, "Torsion Fit this residue","torsion_fit_this",
        lambda _simple_action, _two: torsion_fit_this_func())

    # coot_gui.add_simple_action_to_menu(
    #     menu, "Torsion Fit This Residue and Neighbours",
    #     lambda _simple_action, _two: torsion_fit_this_and_neighbours_func())

    coot_gui.add_simple_action_to_menu(
        menu, "Torsion Fit & Refine this residue","torsion_fit_refine_this",
        lambda _simple_action, _two: torsion_fit_this_func(True))

    coot_gui.add_simple_action_to_menu(
        menu, "Add synthetic pyranose plane restraints","add_synthetic_pyranose_planes",
        lambda _simple_action, _two: add_linked_cho.add_synthetic_pyranose_planes())

    coot_gui.add_simple_action_to_menu(
        menu, "Use Unimodal ring torsion restraints","use_unimodal_pyranose_ring_torsions",
        lambda _simple_action, _two: add_linked_cho.use_unimodal_pyranose_ring_torsions())

    # This should probably become a checkbox??
    # note: this is duplicated in Restraints

    coot_gui.add_simple_action_to_menu(
        menu, "Display Extra Restraints","cho__display_extra_restraints",
        lambda _simple_action, _two: coot_utils.using_active_atom(set_show_extra_restraints, "aa_imol", 1))

    coot_gui.add_simple_action_to_menu(
        menu, "Undisplay Extra Restraints","cho__undisplay_extra_restraints",
        lambda _simple_action, _two: coot_utils.using_active_atom(set_show_extra_restraints, "aa_imol", 0))

    coot_gui.add_simple_action_to_menu(
        menu, "Extract this Tree","new_molecule_from_this_glyco_tree",
        lambda _simple_action, _two:
        add_linked_cho.new_molecule_from_this_glyco_tree())
