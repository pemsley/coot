
import types
import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk

import coot
import coot_utils
import coot_gui
import res_spec_utils
import redefine_functions as rdf

# consider for coot_utils
def residue_spec_to_string(spec):
    # print("DEBUG:: in residue_spec_to_string() with spec", spec)
    if len(spec) == 4:
        return (spec[1] + str(spec[2]) + spec[3])
    else:
        return (spec[0] + str(spec[1]) + spec[2])

# find bad things in the structure - rama, C-beta, rotamer, atom clashes baddies
#
def validation_outliers_dialog(imol, imol_map):

    dialog_vbox = False
    window = False
    missing_sidechains_checkbutton = False
    cg_torsion_diff_checkbutton = False
    poor_density_checkbutton =  False

    def find_rama_baddies():
        rs = coot.all_molecule_ramachandran_score_py(imol)
        scored_residues = rs[5]
        interesting = filter(lambda item: item[2] < 0.02, scored_residues)

        # remove phi-psi and res-names from the return value
        munged = map(lambda item: [item[1], item[2]], interesting)
        return ["Ramachandran Improbables", munged]

    # a list of atom specs
    def find_chiral_volume_baddies():
        r = coot.chiral_volume_errors_py(imol)
        if not isinstance(r, list):
            return []
        else:
            return r

    def make_window_title(n):
        return "Coot Interesting/Outliers/Problems: " + str(n)

    def find_c_beta_baddies():
        try:
            c_beta_deviations
        except:
            return []
        return c_beta_deviations(imol)

    def find_em_ringer_baddies():
        try:
            CG_spin_research
        except:
            return []
        if not isinstance(scored_residues, list):
            return []
        else:
            scored_residues = CG_spin_research(imol, imol_map)
            interesting = filter(lambda item: item[1] > 30 or item[1] < -30,
                                 scored_residues)
            # return a list of residue specs - not what is expected?
            return interesting

    def rotamer_score_residues(imol):

        residues = coot_utils.all_residues_sans_water(imol)
        ret = map(lambda residue_spec: [residue_spec,
                                        coot.rotamer_score(imol,
                                                      res_spec_utils.residue_spec_to_chain_id(residue_spec),
                                                      res_spec_utils.residue_spec_to_res_no(residue_spec),
                                                      coot_utils.residue_spec_to_ins_code(residue_spec),
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
                                    coot_utils.residue_spec_to_ins_code(spec))

            if res_name in ["ALA", "GLY", "UNK", "HOH"]:
                pass
            elif (score == 0.0 and not ok_to_have_missing_sidechain_buttons_qm()):
                pass
            else:
                # if spec is a het-group then no rotamers for that (return False)
                is_het = any(map(lambda item: residue_specs_match_qm(item, spec),
                                 het_groups_in_mol))
                if is_het:
                    pass
                else:
                    if res_name == "LEU":
                        if score < 1.1:
                            ret += [baddie]
                    else:
                        if score < 0.021:
                            ret += [baddie]
        return ret

    def molecule_atom_overlap_baddies():
        l = coot.molecule_atom_overlaps_py(imol)
        if isinstance(l, list):
            return l
        else:
            info_dialog(l)
            return []

    def filter_molecule_atom_overlap_baddies(mao_baddies):
        baddie_limit = 2.0  # more than this is marked as a baddie, was 2.2. Is 2.0 good?
        return filter(lambda mao_item: mao_item['overlap-volume'] > baddie_limit, mao_baddies)

    def non_pro_cis_peptide_baddies():  # do the filter here - just for consistency
        cis_peps = coot.cis_peptides_py(imol)
        return filter(lambda peptide: "PRO" != coot_utils.residue_spec_to_residue_name(imol, peptide[1]), cis_peps)

    def twisted_trans_peptide_baddies():
        return coot.twisted_trans_peptides_py(imol)

    def destroy_buttons_with_label(label_fragment_txt, dialog_vbox):
        current_buttons = dialog_vbox.get_children()
        for button in current_buttons:
            label = button.get_label()
            if label_fragment_txt in label:
                button.destroy()

    def regenerate_button_fn(*args):

        def cb_func(button, callback):
            callback("some-arg")

        if dialog_vbox:
            old_buttons = dialog_vbox.get_children()
            for butt in old_buttons:
                butt.destroy()
            buttons = make_buttons()
            for button_parts in buttons:
                # print("button parts 0", button_parts[0])
                # print("button parts 1", button_parts[1])
                button = Gtk.Button(button_parts[0])
                button.connect("clicked", cb_func, button_parts[1])
                dialog_vbox.pack_start(button, False, False, 2)
                button.show()
            if window:
                window.set_title(make_window_title(len(buttons)))
            # validation_outliers_dialog(imol, imol_map)

    def ok_to_do_density_correlations_qm():
        if poor_density_checkbutton:
            return poor_density_checkbutton.get_active()
        else:
            return True

    def ok_to_have_missing_sidechain_buttons_qm():
        if missing_sidechains_checkbutton:
            return missing_sidechains_checkbutton.get_active()
        else:
            return True

    def ok_to_do_CG_torsion_diffs_qm():
        try:
            CG_spin_search
        except:
            return False
        return cg_torsion_diff_checkbutton.get_active()

    def go_to_atom_in_molecule(imol, atom_spec):
        coot.set_go_to_atom_molecule(imol)
        rdf.set_go_to_atom_from_atom_spec(atom_spec)

    def go_to_residue_in_molecule(imol, res_spec):
        coot.set_go_to_atom_molecule(imol)
        rdf.set_go_to_atom_from_res_spec(res_spec)

    def generator(res_spec):
        func = lambda res_spec_c=res_spec : go_to_residue_in_molecule(imol, res_spec_c)
        def action(arg):
            func(res_spec)
        return action

    def generator_for_atom_spec(atom_spec):
        func = lambda atom_spec_c=atom_spec : go_to_atom_in_molecule(imol, atom_spec_c)
        def action(arg):
            func(atom_spec)
        return action

    def make_buttons():

        frb = find_rama_baddies()
        fcbb = find_c_beta_baddies()
        filtered_mao_baddies = filter_molecule_atom_overlap_baddies(molecule_atom_overlap_baddies())
        residue_correlations = [] if not ok_to_do_density_correlations_qm() \
                               else coot.map_to_model_correlation_per_residue_py(imol,
                                                                                 coot_utils.all_residues_sans_water(imol),
                                                                                 0,
                                                                                 imol_map)

        def get_rama_prob(val):
            return val[1]

        # rama
        baddies = filter(lambda baddie: baddie[1] < 0.002, frb[1])
        # baddies.sort(key=get_rama_prob)
        sorted_filtered_rama_baddies = baddies

        # c-beta
        c_beta_baddies = filter(lambda baddie: baddie[1][0][1] > 0.25, fcbb)
        def get_c_beta_score(val):
            return val[1][0][1]
        # c_beta_baddies.sort(key=get_c_beta_score, reverse=True)
        sorted_filtered_c_beta_baddies = c_beta_baddies

        # Density correlations
        #
        density_baddies = filter(lambda baddie: baddie[1] < 0.8, residue_correlations)

        # CG Torsion
        #
        cg_torsion_baddies = find_em_ringer_baddies()

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
            score_string = '{:6.2f} %'.format(100 * rama_prob)
            button_label = "Ramachandran Outlier " + \
                           res_spec_utils.residue_spec_to_chain_id(spec) + \
                           " " + \
                           str(coot_utils.residue_spec_to_res_no(spec)) + \
                           coot_utils.residue_spec_to_ins_code(spec) + \
                           " " + \
                           coot_utils.residue_spec_to_residue_name(imol, spec) + \
                           " " + \
                           score_string
            rama_buttons.append([button_label, generator(spec)])

        c_beta_buttons = []
        for baddie in sorted_filtered_c_beta_baddies:
            spec = baddie[0]
            score = baddie[1][0][1]  # only the first score
            score_string = '{:6.2f}'.format(score)
            button_label = "C-beta deviant " + \
                           residue_spec_to_string(spec) + \
                           " " + \
                           coot_utils.residue_spec_to_residue_name(imol, spec) + \
                           " " + \
                           score_string + u'\u212B'.encode('utf-8')
            c_beta_buttons.append([button_label, generator(residue_spec)])

        non_pro_cis_peptide_buttons = []
        for baddie in non_pro_cis_peptide_baddies():
            spec_1 = baddie[0]
            spec_2 = baddie[1]
            omega = baddie[2]
            button_label = "Non-PRO cis-peptide " + \
                           residue_spec_to_string(spec_1) + \
                           " - " + \
                           residue_spec_to_string(spec_2)
            non_pro_cis_peptide_buttons.append([button_label, generator(spec_1)])

        twisted_trans_peptide_buttons = []
        for baddie in twisted_trans_peptide_baddies():
            spec_1 = baddie[0]
            spec_2 = baddie[1]
            omega = baddie[2]
            button_label = "Twisted trans-peptide " + \
                           residue_spec_to_string(spec_1) + \
                           " - " + \
                           residue_spec_to_string(spec_2)
            twisted_trans_peptide_buttons.append([button_label, generator(spec_1)])

        rota_buttons = []
        for baddie in filtered_rotamer_baddies:
            spec, score = baddie

            # Paul is not sure that he likes a score of
            # 0.0 meaning "Missing sidechain"
            # we have lost some information on the way.
            # 20200511-PE Yeah, like the fact that the residue was RNA!
            #
            score_string = '{:6.2f} %'.format(score)
            ms_string = "Missing Sidechain" if score == 0.0 else "Rotamer Outlier"
            rot_name = coot.get_rotamer_name_py(imol,
                                        res_spec_utils.residue_spec_to_chain_id(spec),
                                        res_spec_utils.residue_spec_to_res_no(spec),
                                        coot_utils.residue_spec_to_ins_code(spec))
            button_label = ms_string + " " + \
                           residue_spec_to_string(spec) + \
                           " " + \
                           coot_utils.residue_spec_to_residue_name(imol, spec) + \
                            " "
            button_label += rot_name if isinstance(rot_name, str) else " "
            button_label += "" if score == 0.0 else score_string
            rota_buttons.append([button_label, generator(spec)])

        density_baddies_buttons = []
        for baddie in density_baddies:
            spec, score = baddie
            button_label = "Poor Density Fit " + \
                           residue_spec_to_string(spec) + \
                           " " + \
                           '{:5.2f}'.format(score)
            density_baddies_buttons.append([button_label, generator(spec)])

        cg_torsion_buttons = []
        for baddie in cg_torsion_baddies:
            spec, score = baddie
            button_label = "CG Torsion Diff " + \
                           residue_spec_to_string(spec) + \
                           " " + \
                           '{:5.2f}'.format(score)
            cg_torsion_buttons.append([button_label, generator(spec)])

        chiral_volume_buttons = []
        for baddie_atom_spec in find_chiral_volume_baddies():
            atom_spec = baddie_atom_spec
            if len(atom_spec) == 6:
                atom_spec = atom_spec[1:]
            button_label = "Chiral Volume Error " + coot_utils.atom_spec_to_string(atom_spec)
            func = generator_for_atom_spec(atom_spec)
            chiral_volume_buttons.append([button_label, func])

        atom_overlap_buttons = []
        for baddie in filtered_mao_baddies:
            atom_spec_1 = baddie['atom-1-spec']
            atom_spec_2 = baddie['atom-2-spec']
            overlap = baddie['overlap-volume']
            button_label = "Atom Overlap " + \
                           coot_utils.atom_spec_to_string(atom_spec_1) + \
                           " on " + \
                           coot_utils.atom_spec_to_string(atom_spec_2) + \
                           " OV: " + \
                           '{:5.2f}'.format(overlap)
            atom_overlap_buttons.append([button_label, generator_for_atom_spec(atom_spec_1)])

        buttons =  chiral_volume_buttons + \
                   rama_buttons + \
                   rota_buttons + \
                   non_pro_cis_peptide_buttons + \
                   twisted_trans_peptide_buttons + \
                   density_baddies_buttons + \
                   c_beta_buttons + \
                   cg_torsion_buttons + \
                   atom_overlap_buttons
        return buttons

    # main line

    buttons = make_buttons()

    dialog_vbox, window = coot_gui.dialog_box_of_buttons(make_window_title(len(buttons)),
                                                [350, 400], buttons, " Close ")

    window_bits = window.get_children()
    vbox_outer = window_bits[0]
    control_button_vbox_1 = Gtk.HBox(False, 2)
    missing_sidechains_checkbutton_local = Gtk.CheckButton("Missing Sidechains")
    poor_density_checkbutton_local = Gtk.CheckButton("Poor Density Fit")
    cg_torsion_diff_checkbutton_local = Gtk.CheckButton("CG Torsion Diff.")
    regenerate_button_local = Gtk.Button("Update")

    missing_sidechains_checkbutton = missing_sidechains_checkbutton_local
    poor_density_checkbutton = poor_density_checkbutton_local
    cg_torsion_diff_checkbutton = cg_torsion_diff_checkbutton_local

    missing_sidechains_checkbutton.set_active(True)
    poor_density_checkbutton.set_active(True)
    cg_torsion_diff_checkbutton.set_active(True)
    control_button_vbox_1.pack_start(missing_sidechains_checkbutton, False, False, 2)
    control_button_vbox_1.pack_start(poor_density_checkbutton, False, False, 2)

    vbox_outer.pack_start(regenerate_button_local, False, False, 6)
    vbox_outer.pack_start(control_button_vbox_1, False, False, 2)

    regenerate_button_local.connect("clicked", regenerate_button_fn)

    def missing_sidechains_checkbutton_toggled(widget):
        state = widget.get_active()
        if not state:
            # i.e. no buttons with "Missing Sidechain"
            # Mmmh, really?!?
            destroy_buttons_with_label("Missing Sidechain", dialog_vbox)

    missing_sidechains_checkbutton.connect("toggled",
                                           missing_sidechains_checkbutton_toggled)

    def poor_density_checkbutton_toggled(widget):
        state = widget.get_active()
        if not state:
            destroy_buttons_with_label("Poor Density", dialog_vbox)

    poor_density_checkbutton.connect("toggled",
                                     poor_density_checkbutton_toggled)

    # 20190102-PE depends on the version of coot that we are using
    #
    have_cg_spin = False
    try:
        CG_spin_search
        have_cg_spin = True
    except:
        pass
    if have_cg_spin:
        control_button_vbox_1.pack_start(cg_torsion_diff_checkbutton, False, False, 2)
        cg_torsion_diff_checkbutton.show()

        def cg_torsion_diff_checkbutton_toggled(widget):
            state = widget.get_active()
            if not state:
                destroy_buttons_with_label("CG TOrsion", dialog_vbox)
        cg_torsion_diff_checkbutton.connect("toggled",
                                            cg_torsion_diff_checkbutton_toggled)

    control_button_vbox_1.show()
    missing_sidechains_checkbutton.show()
    poor_density_checkbutton.show()
    regenerate_button_local.show()

