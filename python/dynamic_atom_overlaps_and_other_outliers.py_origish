#
def make_quick_test_validation_info(imol):

    def find_rama_baddies():
        rs = all_molecule_ramachandran_score(imol)
        # print "rs:", rs
        scored_residues = rs[5]
        interesting = filter(lambda item: item[2] < 0.02, scored_residues)

        # remove phi-psi and res-names from the return value
        munged = map(lambda item: [item[1], item[2]], interesting)
        return ["Ramachandran Improbables", munged]

    # a list of atom specs
    def find_chiral_volume_baddies():
        r = chiral_volume_errors(imol)
        if not isinstance(r, list):
            return []
        else:
            return r

    def find_c_beta_baddies():
        try:
            c_beta_deviations
            return c_beta_deviations(imol)
        except:
            return []

    def rotamer_score_residues(imol):

        residues = all_residues_sans_water(imol)
        ret = map(lambda residue_spec: [residue_spec,
                                        rotamer_score(imol,
                                                      residue_spec_to_chain_id(residue_spec),
                                                      residue_spec_to_res_no(residue_spec),
                                                      residue_spec_to_ins_code(residue_spec),
                                                      "")],
                  residues)
        return ret

    def filter_rotamer_baddies(baddies):

        het_groups_in_mol = het_group_residues(imol)

        ret = []

        for baddie in baddies:
            spec, score = baddie
            res_name = residue_name(imol,
                                    residue_spec_to_chain_id(spec),
                                    residue_spec_to_res_no(spec),
                                    residue_spec_to_ins_code(spec))

            # print "filter-rotamers testing baddie:", baddie
            if res_name in ["ALA", "GLY", "UNK", "HOH"]:
                pass
            else:
                # if spec is a het-group then no rotamers for that (return False)
                is_het = any(map(lambda item: residue_specs_match_qm(item, spec),
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
        return molecule_atom_overlaps(imol)

    def filter_molecule_atom_overlap_baddies(mao_baddies):
        baddie_limit = 2.0  # more than this is marked as a baddie, was 2.2. Is 2.0 good?
        return filter(lambda mao_item: mao_item['overlap-volume'] > baddie_limit, mao_baddies)

    def non_pro_cis_peptide_baddies():  # do the filter here - just for consistency
        cis_peps = cis_peptides(imol)
        return filter(lambda peptide: "PRO" !=
                      residue_spec_to_residue_name(imol,
                                                   peptide[1]),
                      cis_peps)

    def sort_buttons_inner(buttons_with_spec):
        def sort_fn(b_1, b_2):
            spec_1 = b_1[0]
            spec_2 = b_2[0]
            return residue_spec_less_than(b_1, b_2)

        if buttons_with_spec:
            ret = sorted(buttons_with_spec, key=lambda x: (x[0][0], x[0][1]))
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
        return twisted_trans_peptides(imol)

    def destroy_buttons_with_label(label_fragment_txt, dialog_vbox):
        current_buttons = dialog_vbox.get_children()
        for button in current_buttons:
            label = button.get_label()
            if label_fragment_txt in label:
                button.destroy()

    # now do something....
    frb = find_rama_baddies()
    fcbb = find_c_beta_baddies()
    maob = molecule_atom_overlap_baddies()
    filtered_mao_baddies = filter_molecule_atom_overlap_baddies(maob)

    # rama
    baddies = filter(lambda baddie: baddie[1] < 0.002, frb[1])
    def get_rama_prob(val):
        return val[1]
    baddies.sort(key=get_rama_prob)
    sorted_filtered_rama_baddies = baddies

    # c-beta
    c_beta_baddies = filter(lambda baddie: baddie[1][0][1] > 0.25, fcbb)
    def get_c_beta_score(val):
        return val[1][0][1]
    c_beta_baddies.sort(key=get_c_beta_score, reverse=True)
    sorted_filtered_c_beta_baddies = c_beta_baddies

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
    sorted_filtered_rotamer_baddies = filtered_rotamer_baddies.sort(key=get_rotamer_score)

    rama_buttons = []
    for baddie in sorted_filtered_rama_baddies:
        spec, rama_prob = baddie
        baddie_score = 0.5
        score_string = '{:6.2f} %'.format(100 * rama_prob)
        button_label = "Ramachandran Outlier " + \
                       residue_spec_to_chain_id(spec) + \
                       " " + \
                       str(residue_spec_to_res_no(spec)) + \
                       residue_spec_to_ins_code(spec) + \
                       " " + \
                       residue_spec_to_residue_name(imol, spec) + \
                       " " + \
                       score_string
        rama_buttons.append([spec, baddie_score, button_label,
                             [[set_go_to_atom_molecule, imol],
                              [set_go_to_atom_from_res_spec, spec]]])

    c_beta_buttons = []
    for baddie in sorted_filtered_c_beta_baddies:
        spec = baddie[0]
        score = baddie[1][0][1]  # only the first score
        score_string = '{:6.2f}'.format(score)
        baddie_score = 0.5
        button_label = "C-beta deviant " + \
                       residue_spec_to_string(spec) + \
                       " " + \
                       residue_spec_to_residue_name(imol, spec) + \
                       " " + \
                       score_string + u'\u212B'.encode('utf-8')
        c_beta_buttons.append([spec, baddie_score, button_label,
                               [[set_go_to_atom_molecule, imol],
                                [set_go_to_atom_from_res_spec, spec]]])

    non_pro_cis_peptide_buttons = []
    for baddie in non_pro_cis_peptide_baddies():
        spec_1 = baddie[0]
        spec_2 = baddie[1]
        omega = baddie[2]
        baddie_score = 0.5
        button_label = "Non-PRO cis-peptide " + \
                       residue_spec_to_string(spec_1) + \
                       " - " + \
                       residue_spec_to_string(spec_2)
        non_pro_cis_peptide_buttons.append([spec_1, baddie_score, button_label,
                                            [[set_go_to_atom_molecule, imol],
                                             [set_go_to_atom_from_res_spec, spec_1]]])

    twisted_trans_peptide_buttons = []
    for baddie in twisted_trans_peptide_baddies():
        spec_1 = baddie[0]
        spec_2 = baddie[1]
        omega = baddie[2]
        baddie_score = 1.0
        button_label = "Twisted trans-peptide " + \
                       residue_spec_to_string(spec_1) + \
                       " - " + \
                       residue_spec_to_string(spec_2)
        twisted_trans_peptide_buttons.append([spec_1, baddie_score, button_label,
                                              [[set_go_to_atom_molecule, imol],
                                               [set_go_to_atom_from_res_spec, spec_1]]])

    rota_buttons = []
    for baddie in sorted_filtered_rama_baddies:
        spec, score = baddie

        # Paul is not sure that he likes a score of
        # 0.0 meaning "Missing sidechain"
        # we have lost some information on the way
        #
        score_string = '{:6.2f} %'.format(score)
        ms_string = "Missing Sidechain" if score == 0.0 else "Rotamer Outlier"
        baddie_score = 0.5
        rot_name = get_rotamer_name(imol,
                                    residue_spec_to_chain_id(spec),
                                    residue_spec_to_res_no(spec),
                                    residue_spec_to_ins_code(spec))
        button_label = ms_string + " " + \
                       residue_spec_to_string(spec) + \
                       " " + \
                       residue_spec_to_residue_name(imol, spec) + \
                        " "
        button_label += rot_name if isinstance(rot_name, str) else " "
        button_label += "" if score == 0.0 else score_string
        rota_buttons.append([spec, baddie_score, button_label,
                             [[set_go_to_atom_molecule, imol],
                              [set_go_to_atom_from_res_spec, spec]]])

    chiral_volume_buttons = []
    for baddie_atom_spec_6 in find_chiral_volume_baddies():
        # strip off leading, incorrect imol
        baddie_atom_spec = baddie_atom_spec_6[1:]
        baddie_score = 1.0
        residue_spec = atom_spec_to_residue_spec(baddie_atom_spec)
        button_label = "Chiral Volume Error " + \
                       atom_spec_to_string(baddie_atom_spec)
        chiral_volume_buttons.append([residue_spec, baddie_score, button_label,
                                      [[set_go_to_atom_molecule, imol],
                                       [set_go_to_atom_from_atom_spec, baddie_atom_spec]]])

    atom_overlap_buttons = []
    for baddie in filtered_mao_baddies:
        atom_spec_1 = baddie['atom-1-spec']
        atom_spec_2 = baddie['atom-2-spec']
        overlap = baddie['overlap-volume']
        baddie_score = 0.5
        residue_spec = atom_spec_to_residue_spec(atom_spec_1)
        button_label = "Atom Overlap " + \
                       atom_spec_to_string(atom_spec_1) + \
                       " on " + \
                       atom_spec_to_string(atom_spec_2) + \
                       " OV: " + \
                       '{:5.2f}'.format(overlap)
        atom_overlap_buttons.append([residue_spec, baddie_score, button_label,
                                     [[set_go_to_atom_molecule, imol],
                                      [set_go_to_atom_from_atom_spec, atom_spec_1]]])


    # This gives a list in "baddie-type" order.
    # If we want a list in Chain/Residue order,
    #    each baddie now is associated with (prefixed by)
    #    a residue spec - and use those to sort residues.

    # these buttons have 3 fields - spec label func

    buttons =  chiral_volume_buttons + \
               rama_buttons + \
               rota_buttons + \
               non_pro_cis_peptide_buttons + \
               twisted_trans_peptide_buttons + \
               c_beta_buttons + \
               atom_overlap_buttons

    sorted_buttons = sort_buttons(buttons)
    sorted_buttons.reverse()
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
            for item in callback:
                item[0](*item[1:])

        if dialog_vbox:
            buttons = make_quick_test_validation_buttons(imol)
            old_buttons = dialog_vbox.get_children()
            for button_spec in buttons:
                button = gtk.Button(button_spec[0])
                button.connect("clicked", cb_func, button_spec[1])
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
            # i.e. no buttons with "Missing Sidechain"
            # Mmmh, really?!?
            destroy_buttons_with_label("Missing Sidechain", dialog_vbox)

    # main line

    #  return the update button so that the caller can emit a "clicked" event
    buttons = make_quick_test_validation_buttons(imol)

    dialog_vbox, window = dialog_box_of_buttons(make_window_title(len(buttons)),
                                                [360, 400], buttons, " Close ")

    window_bits = window.get_children()
    vbox_outer = window_bits[0]
    control_button_vbox_1 = gtk.HBox(False, 2)
    missing_sidechains_checkbutton_local = gtk.CheckButton("Missing Sidechains")
    regenerate_button_local = gtk.Button("Update")

    missing_sidechains_checkbutton = missing_sidechains_checkbutton_local

    vbox_outer.pack_start(regenerate_button_local, False, False, 6)
    vbox_outer.pack_start(control_button_vbox_1, False, False, 2)
    regenerate_button_local.connect("clicked", regenerate_button_fn)
    missing_sidechains_checkbutton.connect("toggled",
                                           missing_sidechains_checkbutton_toggled)

    control_button_vbox_1.show()
    missing_sidechains_checkbutton.show()
    regenerate_button_local.show()
    return regenerate_button_local

def make_quick_test_validation_buttons(imol):
    i = make_quick_test_validation_info(imol)
    # remove the spec (and the baddie score for now)
    return map(lambda x: x[2:], i)

if (have_coot_python):
    if coot_python.main_menubar():
        menu = coot_menubar_menu("Validate")

        def make_quick_test_validation_dialog_func():
            with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                   aa_ins_code, aa_atom_name, aa_alt_conf]:
                update_button = quick_test_validation_outliers_dialog(aa_imol)
                global post_manipulation_script
                def post_manipulation_script(*args):
                    update_button.emit("clicked")

        add_simple_coot_menu_menuitem(
        menu, "Overlaps, Peptides, CBeta, Rama & Rota Outliers",
        lambda func: make_quick_test_validation_dialog_func()
        )


