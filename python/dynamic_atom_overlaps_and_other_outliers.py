#
import coot
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
import coot_gui
import coot_gui_api
import coot_utils
import res_spec_utils
import redefine_functions as rdf

def make_quick_test_validation_info(imol):

    def find_rama_baddies():
        rs = coot.all_molecule_ramachandran_score_py(imol)
        # print "rs:", rs
        scored_residues = rs[5]
        interesting = list(filter(lambda item: item[2] < 0.02, scored_residues))

        # remove phi-psi and res-names from the return value
        munged = list(map(lambda item: [item[1], item[2]], interesting))
        return ["Ramachandran Improbables", munged]

    # a list of atom specs
    def find_chiral_volume_baddies():
        r = coot.chiral_volume_errors_py(imol)
        if not isinstance(r, list):
            return []
        else:
            return r

    def find_c_beta_baddies():
        try:
            l = coot.c_beta_deviations_py(imol)
            # print("debug:: find_c_beta_baddies l", l)
            return l
        except:
            return []

    def rotamer_score_residues(imol):

        residues = coot_utils.all_residues_sans_water(imol)
        ret = map(lambda residue_spec: [residue_spec,
                                        coot.rotamer_score(imol,
                                                      res_spec_utils.residue_spec_to_chain_id(residue_spec),
                                                      res_spec_utils.residue_spec_to_res_no(residue_spec),
                                                      res_spec_utils.residue_spec_to_ins_code(residue_spec),
                                                      "")],
                  residues)
        return ret

    def filter_rotamer_baddies(baddies):

        het_groups_in_mol = coot.het_group_residues_py(imol)

        ret = []

        for baddie in baddies:
            spec, score = baddie
            res_name = coot.residue_name(imol,
                                    res_spec_utils.residue_spec_to_chain_id(spec),
                                    res_spec_utils.residue_spec_to_res_no(spec),
                                    res_spec_utils.residue_spec_to_ins_code(spec))

            # print "filter-rotamers testing baddie:", baddie
            if res_name in ["ALA", "GLY", "UNK", "HOH",
                            "DA", "DG", "DT", "DC",
                            "A", "G", "U", "C"]:
                pass
            else:
                # if spec is a het-group then no rotamers for that (return False)
                is_het = any(map(lambda item: coot_utils.residue_specs_match_qm(item, spec),
                                 het_groups_in_mol))
                if is_het:
                    pass
                else:
                    if res_name == "LEU":
                        # there is something strange with LEU rotamers, D 427?
                        if score < 1.1:
                            ret += [baddie]
                    else:
                        if score < 0.021:
                            ret += [baddie]
        return ret

    def molecule_atom_overlap_baddies():
        return coot.molecule_atom_overlaps_py(imol)

    def filter_molecule_atom_overlap_baddies(mao_baddies):
        baddie_limit = 2.0  # more than this is marked as a baddie, was 2.2. Is 2.0 good?
        return filter(lambda mao_item: mao_item['overlap-volume'] > baddie_limit, mao_baddies)

    def non_pro_cis_peptide_baddies():  # do the filter here - just for consistency
        cis_peps = coot.cis_peptides_py(imol)
        return filter(lambda peptide: "PRO" != coot_utils.residue_spec_to_residue_name(imol, peptide[1]), cis_peps)

    def sort_buttons_inner(buttons_with_spec):
        if buttons_with_spec:
			# sorted takes care of mixed lists...
            ret = sorted(buttons_with_spec)
        else:
            ret = []
        return ret

    def sort_buttons(buttons):
        def remove_from_list(ls):
            del ls[0][0]
            return ls
        buttons_1 = map(lambda button: button if len(button[0]) == 3
                        else remove_from_list(button), buttons)
        buttons_2 = sort_buttons_inner(buttons_1)
        return buttons_2

    def twisted_trans_peptide_baddies():
        return coot.twisted_trans_peptides_py(imol)

    def set_go_to_imol_and_spec(imol, atom_spec):
        coot.set_go_to_atom_molecule(imol)
        # if we were actually given a residue spec, we should use that!
        if len(atom_spec) == 3:
            rdf.set_go_to_atom_from_res_spec(atom_spec) # it's a res spec
        else:
            rdf.set_go_to_atom_from_atom_spec(atom_spec)

    # now do something....
    frb = find_rama_baddies()
    fcbb = find_c_beta_baddies()
    maob = molecule_atom_overlap_baddies()
    filtered_mao_baddies = filter_molecule_atom_overlap_baddies(maob)

    # rama
    baddies = list(filter(lambda baddie: baddie[1] < 0.002, list(frb[1])))
    def get_rama_prob(val):
        return val[1]
    sorted_filtered_rama_baddies = sorted(baddies, key=get_rama_prob)

    c_beta_baddies = filter(lambda baddie: baddie[1][''] > 0.25, fcbb)
    def get_c_beta_score(val):
        distortion = val[1]['']
        return distortion
    sorted_filtered_c_beta_baddies = sorted(c_beta_baddies, key=get_c_beta_score, reverse=True)

    # Rotamers
    #
    rotamer_baddies = rotamer_score_residues(imol)
    filtered_rotamer_baddies = filter_rotamer_baddies(rotamer_baddies)
    def get_rotamer_score(val):
        if val[1] == 0.0:
            # score things with score 0.0 (meaning missing sidechain)
            # as if they are better than low probability outliers
            # give high score in sorting
            return 10
        else:
            return val[1]
    sorted_filtered_rotamer_baddies = sorted(filtered_rotamer_baddies, key=get_rotamer_score)

    rama_buttons = []
    for baddie in sorted_filtered_rama_baddies:
        spec, rama_prob = baddie
        baddie_score = 0.5
        score_string = '{:6.2f} %'.format(100 * rama_prob)
        button_label = "Ramachandran Outlier " + \
                       res_spec_utils.residue_spec_to_chain_id(spec) + \
                       " " + \
                       str(res_spec_utils.residue_spec_to_res_no(spec)) + \
                       res_spec_utils.residue_spec_to_ins_code(spec) + \
                       " " + \
                       coot_utils.residue_spec_to_residue_name(imol, spec) + \
                       " " + \
                       score_string
        def callback_generator(imol, atom_spec):
            return lambda thing : set_go_to_imol_and_spec(imol, atom_spec)
        rama_buttons.append([spec, baddie_score, button_label, callback_generator(imol, spec)])

    c_beta_buttons = []
    for baddie in sorted_filtered_c_beta_baddies:
        spec = baddie[0]
        # score = baddie[1][0][1]  # only the first score
        score = baddie[1]['']  # only the first score
        score_string = '{:6.2f}'.format(score)
        baddie_score = 0.5
        # button_label = "C-beta deviant " + \
        #                coot_utils.residue_spec_to_string(spec) + \
        #                " " + \
        #                coot_utils.residue_spec_to_residue_name(imol, spec) + \
        #                " " + \
        #                score_string + u'\u212B'.encode('utf-8')
        lab = "C-beta deviant "
        lab += coot_utils.residue_spec_to_string(spec)
        lab += " "
        lab += coot_utils.residue_spec_to_residue_name(imol, spec)
        lab += " "
        lab += score_string
        button_label = lab
        def callback_generator(imol, atom_spec_1):
            return lambda thing : set_go_to_imol_and_spec(imol, atom_spec_1)
        c_beta_buttons.append([spec, baddie_score, button_label, callback_generator(imol, spec)])

    non_pro_cis_peptide_buttons = []
    for baddie in non_pro_cis_peptide_baddies():
        spec_1 = baddie[0]
        spec_2 = baddie[1]
        omega = baddie[2]
        baddie_score = 0.5
        button_label = "Non-PRO cis-peptide " + \
                       coot_utils.residue_spec_to_string(spec_1) + \
                       " - " + \
                       coot_utils.residue_spec_to_string(spec_2)
        def callback_generator(imol, atom_spec_1):
            return lambda thing : set_go_to_imol_and_spec(imol, atom_spec_1)
        non_pro_cis_peptide_buttons.append([spec_1, baddie_score, button_label, callback_generator(imol, spec_1)])

    twisted_trans_peptide_buttons = []
    for baddie in twisted_trans_peptide_baddies():
        spec_1 = baddie[0]
        spec_2 = baddie[1]
        omega = baddie[2]
        baddie_score = 1.0
        button_label = "Twisted trans-peptide " + \
                       coot_utils.residue_spec_to_string(spec_1) + \
                       " - " + \
                       coot_utils.residue_spec_to_string(spec_2)
        def callback_generator(imol, atom_spec_1):
            return lambda thing : set_go_to_imol_and_spec(imol, atom_spec_1)
        twisted_trans_peptide_buttons.append([spec_1, baddie_score, button_label, callback_generator(imol, spec_1)])

    rota_buttons = []
    for baddie in sorted_filtered_rotamer_baddies:
        spec, score = baddie

        # Paul is not sure that he likes a score of
        # 0.0 meaning "Missing sidechain"
        # we have lost some information on the way.
        # 20200511-PE Yeah, like the fact that the residue was RNA!
        #
        score_string = '{:6.2f} %'.format(score)
        ms_string = "Missing Sidechain" if score == 0.0 else "Rotamer Outlier"
        baddie_score = 0.5
        rot_name = coot.get_rotamer_name_py(imol,
                                    res_spec_utils.residue_spec_to_chain_id(spec),
                                    res_spec_utils.residue_spec_to_res_no(spec),
                                    res_spec_utils.residue_spec_to_ins_code(spec))
        button_label = ms_string + " " + \
                       coot_utils.residue_spec_to_string(spec) + \
                       " " + \
                       coot_utils.residue_spec_to_residue_name(imol, spec) + \
                        " "
        button_label += rot_name if isinstance(rot_name, str) else " "
        button_label += "" if score == 0.0 else score_string
        def callback_generator(imol, atom_spec):
            return lambda thing : set_go_to_imol_and_spec(imol, atom_spec)
        rota_buttons.append([spec, baddie_score, button_label, callback_generator(imol, spec)])

    chiral_volume_buttons = []
    for baddie_atom_spec_6 in find_chiral_volume_baddies():
        # strip off leading, incorrect imol
        baddie_atom_spec = baddie_atom_spec_6[1:]
        baddie_score = 1.0
        residue_spec = coot_utils.atom_spec_to_residue_spec(baddie_atom_spec)
        button_label = "Chiral Volume Error " + coot_utils.atom_spec_to_string(baddie_atom_spec)
        def callback_generator(imol, atom_spec):
            return lambda thing : set_go_to_imol_and_spec(imol, atom_spec)
        chiral_volume_buttons.append([residue_spec, baddie_score, button_label, callback_generator(imol, baddie_atom_spec)])

    atom_overlap_buttons = []
    for baddie in filtered_mao_baddies:
        atom_spec_1 = baddie['atom-1-spec']
        atom_spec_2 = baddie['atom-2-spec']
        overlap = baddie['overlap-volume']
        baddie_score = 0.5
        residue_spec = coot_utils.atom_spec_to_residue_spec(atom_spec_1)
        button_label = "Atom Overlap " + \
                       coot_utils.atom_spec_to_string(atom_spec_1) + \
                       " on " + \
                       coot_utils.atom_spec_to_string(atom_spec_2) + \
                       " OV: " + \
                       '{:5.2f}'.format(overlap)
        def callback_generator(imol, atom_spec_1):
            return lambda thing : set_go_to_imol_and_spec(imol, atom_spec_1)
        atom_overlap_buttons.append([residue_spec, baddie_score, button_label, callback_generator(imol, atom_spec_1)])


    # This gives a list in "baddie-type" order.
    # If we want a list in Chain/Residue order,
    #    each baddie now is associated with (prefixed by)
    #    a residue spec - and use those to sort residues.

    # these buttons have 3 fields - spec label func

    buttons = c_beta_buttons + atom_overlap_buttons + chiral_volume_buttons + rama_buttons + rota_buttons + \
        non_pro_cis_peptide_buttons + twisted_trans_peptide_buttons

    sorted_buttons = sort_buttons(buttons)

    return sorted_buttons

# find bad things in the structure - rama, C-beta, rotamer, atom clashes baddies
#
global dialog_vbox, window, missing_sidechains_checkbutton
dialog_vbox = False
window = False
missing_sidechains_checkbutton = False

def quick_test_validation_outliers_dialog(imol):

    global dialog_vbox, window, missing_sidechains_checkbutton

    def make_window_title(n):
        return "Coot Interesting/Outliers/Problems: " + str(n)

    def regenerate_button_fn(*args):

        def cb_func(button, callback):
            callback(button)

        if dialog_vbox:
            buttons = make_quick_test_validation_buttons(imol)
            old_buttons = dialog_vbox.get_children()
            for button_info in buttons:
                button_spec = button_info[0]
                f           = button_info[1]
                label       = button_info[2]
                cb_fun      = button_info[3]
                button = Gtk.Button(label)
                button.connect("clicked", cb_func, cb_fun)
                dialog_vbox.pack_start(button, False, False, 2)
                button.show()
            if window:
                window.set_title(make_window_title(len(buttons)))
            for butt in old_buttons:
                butt.destroy()
        else:
            quick_test_validation_outliers_dialog(imol)

    def ok_to_have_missing_sidechain_buttons_qm():
        if missing_sidechains_checkbutton:
            return missing_sidechains_checkbutton.get_active()
        else:
            return True

    def missing_sidechains_checkbutton_toggled(widget):
        state = widget.get_active()
        if not state:
            destroy_buttons_with_label("Missing Sidechain", dialog_vbox)
        else:
            # update the dialog
            regenerate_button_local.emit("clicked")

    def destroy_buttons_with_label(label_fragment_txt, dialog_vbox):
        current_buttons = dialog_vbox.get_children()
        for button in current_buttons:
            label = button.get_label()
            if label_fragment_txt in label:
                button.destroy()

    # main line

    # return the update button so that the caller can emit a "clicked" event
    buttons_and_other_info = make_quick_test_validation_buttons(imol)
    buttons = [but[2:] for but in buttons_and_other_info]

    dialog_vbox, window = coot_gui.dialog_box_of_buttons(make_window_title(len(buttons)),
                                                         [360, 400], buttons, " Close ")

    window_bits = window.get_children()
    vbox_outer = window_bits[0]
    control_button_vbox_1 = Gtk.HBox(False, 2)
    missing_sidechains_checkbutton_local = Gtk.CheckButton("Missing Sidechains")
    regenerate_button_local = Gtk.Button("Update")

    missing_sidechains_checkbutton = missing_sidechains_checkbutton_local

    missing_sidechains_checkbutton.set_active(True)
    control_button_vbox_1.pack_start(missing_sidechains_checkbutton, False, False, 2)

    vbox_outer.pack_start(regenerate_button_local, False, False, 6)
    vbox_outer.pack_start(control_button_vbox_1, False, False, 2)
    regenerate_button_local.connect("clicked", regenerate_button_fn)
    missing_sidechains_checkbutton.connect("toggled", missing_sidechains_checkbutton_toggled)

    control_button_vbox_1.show()
    missing_sidechains_checkbutton.show()
    regenerate_button_local.show()
    return regenerate_button_local

def make_quick_test_validation_buttons(imol):
    i = make_quick_test_validation_info(imol)
    return i

if True:

    def get_existing_submenu(menu, submenu_label):
       label_dict = {}
       for menu_child in menu.get_children():
         for c in menu_child.get_children():
           try:
             t = c.get_text()
             if t == submenu_label:
               return menu_child
           except KeyError as e:
             pass

    if coot_gui_api.main_menubar():

        def make_quick_test_validation_dialog_func():
            with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code, aa_atom_name, aa_alt_conf]:
                update_button = quick_test_validation_outliers_dialog(aa_imol)
                global post_manipulation_script
                def emit_update_button_clicked_script(*args):
                    update_button.emit("clicked")
                post_manipulation_script = emit_update_button_clicked_script

        menu = coot_gui.coot_menubar_menu("Validate")
        coot_gui.add_simple_coot_menu_menuitem(menu, "Overlaps, Peptides, CBeta, Rama & Rota Outliers",
                                               lambda func: make_quick_test_validation_dialog_func())



