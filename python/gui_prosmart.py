# gui_prosmart.py
# Copyright 2007 by Paul Emsley
# Copyright 2007 by The University of Oxford
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import os
import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, GObject
from gi.repository import Gio
from gi.repository import GLib
import coot
import coot_gui_api
import coot_utils
import coot_gui
from coot_gui import add_simple_action_to_menu, attach_module_menu_button

# target is my molecule, ref is the homologous (high-res) model
#
# extra arg: include_side_chains=False
#
def run_prosmart(imol_target, imol_ref, include_side_chains=False):
    """
    target is my molecule, ref is the homologous (high-res) model

    extra arg: include_side_chains=False
    """

    dir_stub = "coot-ccp4"
    coot.make_directory_maybe(dir_stub)
    r_t_name = coot.molecule_name_stub_py(imol_target, 0).replace(" ", "_");
    r_r_name = coot.molecule_name_stub_py(imol_ref,    0).replace(" ", "_");
    target_pdb_file_name =    os.path.join(dir_stub, r_t_name + "-prosmart.pdb")
    reference_pdb_file_name = os.path.join(dir_stub, r_r_name + "-prosmart-ref.pdb")
    prosmart_out = os.path.join("ProSMART_Output", r_t_name + "-prosmart.txt")

    coot.write_pdb_file(imol_target, target_pdb_file_name)
    coot.write_pdb_file(imol_ref, reference_pdb_file_name)
    prosmart_exe = coot_utils.find_exe("prosmart")
    if prosmart_exe:
        l = ["-p1", target_pdb_file_name,
             "-p2", reference_pdb_file_name,
             "-restrain_seqid", "30"]
        if include_side_chains:
            l += ["-side"]
        coot_utils.popen_command(prosmart_exe,
                                 l,
                                 [],
                                 os.path.join(dir_stub, "prosmart.log"),
                                 False)
        if (not os.path.isfile(prosmart_out)):
            print("file not found", prosmart_out)
        else:
            print("Reading ProSMART restraints from", prosmart_out)
            coot.add_refmac_extra_restraints(imol_target, prosmart_out)
    else:
        coot.info_dialog("No prosmart")



def add_module_restraints():
    menu = attach_module_menu_button("Restraints")

    def generate_all_molecule_self_restraints(val):
        with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                    aa_ins_code, aa_atom_name, aa_alt_conf]:
            generate_self_restraints(aa_imol, val)

    def generate_self_restraint_func(sig):
        with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                    aa_ins_code, aa_atom_name, aa_alt_conf]:
            coot.generate_local_self_restraints(aa_imol, aa_chain_id, sig)

    # Generate self restraints for residues around sphere
    def generate_self_restraint_in_sphere_func():
        with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                    aa_ins_code, aa_atom_name, aa_alt_conf]:
            centred_residue = [aa_chain_id, aa_res_no, aa_ins_code]
            radius = 10
            local_dist_max = 4.2
            other_residues = residues_near_residue(aa_imol,
                                                    centred_residue,
                                                    radius)
            residue_specs = (centred_residue + other_residues) if isinstance(other_residues, list) else centred_residue
            coot.generate_local_self_restraints_by_residues_py(aa_imol,
                                                            residue_specs,
                                                            local_dist_max)

    def display_extra_restraints_func(state):
        with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                    aa_ins_code, aa_atom_name, aa_alt_conf]:
            coot.set_show_extra_restraints(aa_imol, state)

    add_simple_action_to_menu(
        menu, "Generate Self Restraints 3.7 for Chain","generate_self_restraint_37",
        lambda _simple_action, _arg: generate_self_restraint_func(3.7))

    add_simple_action_to_menu(
        menu, "Generate Chain Self Restraints 4.3 for Chain","generate_self_restraint_43",
        lambda _simple_action, _arg: generate_self_restraint_func(4.3))

    add_simple_action_to_menu(
        menu, "Generate Chain Self Restraints 6 for Chain","generate_self_restraint_60",
        lambda _simple_action, _arg: generate_self_restraint_func(6))

    add_simple_action_to_menu(
        menu, "Generate Local Self Restraints 6","generate_self_restraint_in_sphere",
        lambda _simple_action, _arg: generate_self_restraint_in_sphere_func())

    add_simple_action_to_menu(
        menu, "Generate All-molecule Self Restraints 4.3","generate_all_molecule_self_restraints_43",
        lambda _simple_action, _arg: generate_all_molecule_self_restraints(4.3))

    add_simple_action_to_menu(
        menu, "Generate All-molecule Self Restraints 5.0","generate_all_molecule_self_restraints_50",
        lambda _simple_action, _arg: generate_all_molecule_self_restraints(5.0))

    add_simple_action_to_menu(
        menu, "Generate All-molecule Self Restraints 6.0","generate_all_molecule_self_restraints_60",
        lambda _simple_action, _arg: generate_all_molecule_self_restraints(6.0))

    # todo: Move it under R/RC and make it a single toggle

    add_simple_action_to_menu(
        menu, "Undisplay Extra Restraints","undisplay_extra_restraints",
        lambda _simple_action, _arg: display_extra_restraints_func(0))

    add_simple_action_to_menu(
        menu, "Display Extra Restraints","display_extra_restraints",
        lambda _simple_action, _arg: display_extra_restraints_func(1))


def add_module_prosmart():
    
    menu = attach_module_menu_button("ProSMART")

    def combobox_to_molecule_number(combobox):
        imol = -1
        tree_iter = combobox.get_active_iter()
        if tree_iter is not None:
            model = combobox.get_model()
            it = model[tree_iter]
            imol = it[0]
        return imol

    def launch_prosmart_gui():
        def go_button_cb(*args):
            imol_tar = combobox_to_molecule_number(combobox_tar)
            imol_ref = combobox_to_molecule_number(combobox_ref)
            do_side_chains = check_button.get_active()
            run_prosmart(imol_tar, imol_ref, do_side_chains)
            window.destroy()
        window = Gtk.Window()
        hbox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
        vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        h_sep = Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL)
        chooser_hint_text_1 = "Target molecule"
        chooser_hint_text_2 = "Reference (high-res) molecule"
        go_button = Gtk.Button(label="ProSMART")
        cancel_button = Gtk.Button(label="Cancel")
        check_button = Gtk.CheckButton(label="Include Side-chains")

        combobox_tar = coot_gui.generic_molecule_chooser(vbox, chooser_hint_text_1)
        combobox_ref = coot_gui.generic_molecule_chooser(vbox, chooser_hint_text_2)

        vbox.append(check_button)
        vbox.append(h_sep)
        vbox.append(hbox)
        hbox.append(go_button)
        hbox.append(cancel_button)
        window.set_child(vbox)

        cancel_button.connect("clicked", lambda w: w.destroy())

        go_button.connect("clicked", go_button_cb)
        window.show()

    def generate_self_restraint_func(sig):
        with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                    aa_ins_code, aa_atom_name, aa_alt_conf]:
            generate_local_self_restraints(aa_imol, aa_chain_id, sig)

    def prosmart_cut_to_func(sig_low, sig_high):
        with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                    aa_ins_code, aa_atom_name, aa_alt_conf]:
            coot.set_extra_restraints_prosmart_sigma_limits(aa_imol,
                                                        sig_low, sig_high)
    add_simple_action_to_menu(
        menu, "Show Only Deviant Distances Beyond 6","prosmart_cut_to_6",
        lambda _simple_action, _arg: prosmart_cut_to_func(-6, 6))

    add_simple_action_to_menu(
        menu, "Show Only Deviant Distances Beyond 4","prosmart_cut_to_4",
        lambda _simple_action, _arg: prosmart_cut_to_func(-4, 4))

    add_simple_action_to_menu(
        menu, "Show Only Deviant Distances Beyond 2.0","prosmart_cut_to_2",
        lambda _simple_action, _arg: prosmart_cut_to_func(-2, 2))

    add_simple_action_to_menu(
        menu, "Show Only Deviant Distances Beyond 1.0","prosmart_cut_to_1",
        lambda _simple_action, _arg: prosmart_cut_to_func(-1, 1))

    add_simple_action_to_menu(
        menu, "Undisplay All Extra Distance Restraints","prosmart_cut_to_0",
        lambda _simple_action, _arg: prosmart_cut_to_func(0, 0))

    def restraint_to_ca_func(state):
        with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                    aa_ins_code, aa_atom_name, aa_alt_conf]:
            coot.set_extra_restraints_representation_for_bonds_go_to_CA(aa_imol, state)

    # # clutter
    # coot_gui.add_simple_coot_menu_menuitem(
    #     menu, "Restraint Representation To CA",
    #     lambda func: restraint_to_ca_func(1))

    # coot_gui.add_simple_coot_menu_menuitem(
    #     menu, "Restraint Representation To Home Atom",
    #     lambda func: restraint_to_ca_func(0))

    add_simple_action_to_menu(
        menu, "Run ProSMART...","run_prosmart_gui", lambda _simple_action, _arg: launch_prosmart_gui())

    ## extra
    
    # def delete_all_extra_restraints_func():
    #     with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
    #                                aa_ins_code, aa_atom_name, aa_alt_conf]:
    #         coot.delete_all_extra_restraints(aa_imol)

    # coot_gui.add_simple_coot_menu_menuitem(
    #     menu, "Delete All Extra Restraints",
    #     lambda func: delete_all_extra_restraints_func())

