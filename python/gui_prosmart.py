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
from gi.repository import Gtk
import gi
gi.require_version('Gtk', '3.0')
import coot
import coot_gui_api
import coot_utils
import coot_gui

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



def add_module_prosmart_old():
    
    if True:
        if coot_gui_api.main_menubar():
            menu = coot_gui.coot_menubar_menu("ProSMART")

            def generate_all_molecule_self_restraints(val):
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    generate_self_restraints(aa_imol, val)

            def generate_self_restraint_func(sig):
                with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    coot.generate_local_self_restraints(aa_imol, aa_chain_id, sig)

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

            add_simple_coot_menu_menuitem(
                menu, "Generate Chain Self Restraints 3.7 for Chain",
                lambda func: generate_self_restraint_func(3.7))

            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Generate Chain Self Restraints 4.3 for Chain",
                lambda func: generate_self_restraint_func(4.3))

            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Generate Chain Self Restraints 6 for Chain",
                lambda func: generate_self_restraint_func(6))

            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Generate Local Self Restraints 6",
                lambda func: generate_self_restraint_in_sphere_func())

            add_simple_coot_menu_menuitem(
                menu, "Generate All-molecule Self Restraints 4.3",
                lambda func: generate_all_molecule_self_restraints(4.3))

            add_simple_coot_menu_menuitem(
                menu, "Generate All-molecule Self Restraints 5.0",
                lambda func: generate_all_molecule_self_restraints(5.0))

            def display_extra_restraints_func(state):
                with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    coot.set_show_extra_restraints(aa_imol, state)

            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Undisplay Extra Restraints",
                lambda func: display_extra_restraints_func(0))

            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Display Extra Restraints",
                lambda func: display_extra_restraints_func(1))


def add_module_prosmart():
    
    if True:
        if coot_gui_api.main_menubar():
            menu = coot_gui.coot_menubar_menu("ProSMART")

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
                hbox = Gtk.HBox(False, 0)
                vbox = Gtk.VBox(False, 0)
                h_sep = Gtk.HSeparator()
                chooser_hint_text_1 = " Target molecule "
                chooser_hint_text_2 = " Reference (high-res) molecule "
                go_button = Gtk.Button(" ProSMART ")
                cancel_button = Gtk.Button("  Cancel  ")
                check_button = Gtk.CheckButton("Include Side-chains")

                combobox_tar = coot_gui.generic_molecule_chooser(vbox, chooser_hint_text_1)
                combobox_ref = coot_gui.generic_molecule_chooser(vbox, chooser_hint_text_2)

                vbox.pack_start(check_button, False, False, 2)
                vbox.pack_start(h_sep, False, False, 2)
                vbox.pack_start(hbox, False, False, 2)
                hbox.pack_start(go_button, False, False, 6)
                hbox.pack_start(cancel_button, False, False, 6)
                window.add(vbox)

                cancel_button.connect("clicked", lambda w: window.destroy())

                go_button.connect("clicked", go_button_cb)
                window.show_all()

            def generate_self_restraint_func(sig):
                with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    generate_local_self_restraints(aa_imol, aa_chain_id, sig)

            def prosmart_cut_to_func(sig_low, sig_high):
                with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                           aa_ins_code, aa_atom_name, aa_alt_conf]:
                    coot.set_extra_restraints_prosmart_sigma_limits(aa_imol,
                                                               sig_low, sig_high)
            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Show Only Deviant Distances Beyond 6",
                lambda func: prosmart_cut_to_func(-6, 6))

            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Show Only Deviant Distances Beyond 4",
                lambda func: prosmart_cut_to_func(-4, 4))

            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Show Only Deviant Distances Beyond 2.0",
                lambda func: prosmart_cut_to_func(-2, 2))

            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Show Only Deviant Distances Beyond 1.0",
                lambda func: prosmart_cut_to_func(-1, 1))

            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Undisplay All Extra Distance Restraints",
                lambda func: prosmart_cut_to_func(0, 0))

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

            coot_gui.add_simple_coot_menu_menuitem(
                menu, "Run ProSMART...", lambda func: launch_prosmart_gui())

            ## extra
            
            # def delete_all_extra_restraints_func():
            #     with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
            #                                aa_ins_code, aa_atom_name, aa_alt_conf]:
            #         coot.delete_all_extra_restraints(aa_imol)

            # coot_gui.add_simple_coot_menu_menuitem(
            #     menu, "Delete All Extra Restraints",
            #     lambda func: delete_all_extra_restraints_func())

