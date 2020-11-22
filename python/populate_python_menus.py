# populate_python_menus.py
# Copyright 2007, 2008 by Bernhard Lohkamp
# Copyright 2006, 2007, 2008 by The University of York
# Copyright 2015 by Medical Research Council
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
# along with this program; if not, write to the Free Software
# Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA

# import pygtk, gtk, pango

import coot
import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk
import coot_gui_api
import time
import numbers
import coot_utils
import coot_gui
import fitting

# rename this file populate_python_menus.py

if coot_gui_api.main_menubar():

    def coot_menubar_menu_simple(s):
        menu = Gtk.Menu()
        menuitem = Gtk.MenuItem(s)
        menuitem.set_submenu(menu)
        main_menubar = coot_gui_api.main_menubar()
        main_menubar.append(menuitem)
        menuitem.show()
        return menu

    def coot_menubar_menu(menu_label_string):

        def get_menu_bar_label_dict():
            label_dict = {}
            for menu_child in coot_gui_api.main_menubar().get_children():
                label_dict[menu_child.get_children()[0].get_text()] = menu_child
            return label_dict

        try:
            menu_label_dict = get_menu_bar_label_dict()
            return menu_label_dict[menu_label_string].get_submenu()
        except KeyError as e:
            return coot_menubar_menu_simple(menu_label_string)

    def get_existing_submenu(menu, submenu_label):
        label_dict = {}
        for menu_child in menu.get_children():
            for c in menu_child.get_children():
                try:
                    t = c.get_text()
                    # print("########### get_existing_submenu get_text on c:", t)
                    if t == submenu_label:
                        return menu_child
                except KeyError as e:
                    pass
        return None


    # this should be added to the "Refine" module. It is not used here
    # (it doesn't belong in "Calculate")
    def add_refinement_options(submenu_refine):

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_refine,
            "Set Refinement Options...",
            lambda func: coot_gui.refinement_options_gui())

        submenu = Gtk.Menu()
        menuitem2 = Gtk.MenuItem("Peptide Restraints...")

        menuitem2.set_submenu(submenu)
        submenu_refine.append(menuitem2)
        menuitem2.show()

        def add_restr_func1():
            print('Planar Peptide Restraints added')
            add_planar_peptide_restraints()

        coot_gui.add_simple_coot_menu_menuitem(
            submenu,
            "Add Planar Peptide Restraints",
            lambda func: add_restr_func1())


        def add_restr_func2():
            print('Planar Peptide Restraints removed')
            remove_planar_peptide_restraints()

            coot_gui.add_simple_coot_menu_menuitem(
                submenu, "Remove Planar Peptide Restraints",
                lambda func: add_restr_func2())



    calculate_menu = coot_menubar_menu("Calculate")
    draw_menu      = coot_menubar_menu("Draw")
    edit_menu      = coot_menubar_menu("Edit")

    submenu_all_molecule   = Gtk.Menu()
    submenu_maps           = Gtk.Menu()
    submenu_modelling      = Gtk.Menu()
    submenu_ncs            = Gtk.Menu()
    submenu_views          = Gtk.Menu()
    submenu_representation = Gtk.Menu()
    submenu_pisa           = Gtk.Menu()
    submenu_modules        = Gtk.Menu()
    submenu_settings       = Gtk.Menu()

    calculate_all_molecule_menu_item = get_existing_submenu(calculate_menu, "All Molecule...")
    calculate_all_molecule_menu_item.set_submenu(submenu_all_molecule)

    calculate_maps_menu_item = get_existing_submenu(calculate_menu, "Map Tools...")
    calculate_maps_menu_item.set_submenu(submenu_maps)

    calculate_model_menu_item = get_existing_submenu(calculate_menu, "Modelling...")
    calculate_model_menu_item.set_submenu(submenu_modelling)

    calculate_ncs_menu_item = get_existing_submenu(calculate_menu, "NCS Tools...")
    calculate_ncs_menu_item.set_submenu(submenu_ncs)

    calculate_pisa_menu_item = get_existing_submenu(calculate_menu, "PISA...")
    calculate_pisa_menu_item.set_submenu(submenu_pisa)

    calculate_modules_menu_item = get_existing_submenu(calculate_menu, "Modules...")
    calculate_modules_menu_item.set_submenu(submenu_modules)

    draw_representation_menu_item = get_existing_submenu(draw_menu, "Representation Tools...")
    draw_representation_menu_item.set_submenu(submenu_representation)

    edit_settings_menu_item = get_existing_submenu(edit_menu, "Settings...")
    edit_settings_menu_item.set_submenu(submenu_settings)

    draw_views_menu_item = Gtk.MenuItem("Views")
    draw_views_menu_item.set_submenu(submenu_views)
    draw_menu.append(draw_views_menu_item)
    draw_views_menu_item.show()


    # ---------------------------------------------------------------------
    #     Modules
    # ---------------------------------------------------------------------

    def add_calculate_modules_menu(submenu_modules):

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modules, "CCP4...",
            lambda func: coot_gui.add_module_ccp4())

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modules, "SHELX...",
            lambda func: shelx_extensions.add_module_shelx())

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modules, "User-defined Restraints...",
            lambda func: add_module_user_defined_restraints())

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modules, "ProSMART",
            lambda func: gui_prosmart.add_module_prosmart())

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modules, "Carbohydrate",
            lambda func: gui_add_linked_cho.add_module_carbohydrate_gui())

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modules, "Cryo-EM",
            lambda func: coot_gui.add_module_cryo_em())

    # ---------------------------------------------------------------------
    #     Settings
    # ---------------------------------------------------------------------

    def add_edit_settings_menu(submenu_settings):

        submenu = Gtk.Menu()
        menuitem2 = Gtk.MenuItem("Rotate Translate Zone Mode...")

        menuitem2.set_submenu(submenu)
        submenu_settings.append(menuitem2)
        menuitem2.show()

        coot_gui.add_simple_coot_menu_menuitem(
            submenu, "Rotate About Fragment Centre",
            lambda func: set_rotate_translate_zone_rotates_about_zone_centre(1))


        coot_gui.add_simple_coot_menu_menuitem(
            submenu, "Rotate About Second Clicked Atom",
            lambda func: set_rotate_translate_zone_rotates_about_zone_centre(0))


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_settings,
            "Set Density Fit Graph Weight...",
            lambda func: coot_gui.generic_single_entry("set weight (smaller means apparently better fit)",
                                                       str("%.2f" %residue_density_fit_scale_factor()), "Set it",
                                                       lambda text: set_den_gra_func(text)))


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_settings,
            "Set Spin Speed",
            lambda func: coot_gui.generic_single_entry("Set Spin Speed (smaller is slower)",
                                                       str(idle_function_rotate_angle()), "Set it",
                                                       lambda text: set_idle_function_rotate_angle(float(text))))


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_settings, "Nudge Centre...",
            lambda func: coot_gui.nudge_screen_centre_gui())


        def all_mol_symm_func():
            for imol in coot_utils.molecule_number_list():
                if coot_utils.valid_model_molecule_qm(imol):
                    set_symmetry_whole_chain(imol, 1)

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_settings, "All Molecules use \"Near Chains\" Symmetry",
            lambda func: [valid_model_molecule_qm(imol) and
                          set_symmetry_whole_chain(imol, 1) for imol in coot_utils.molecule_number_list()])


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_settings, "Save Graphics Size and Positions",
            lambda func: graphics_window_size_and_position_to_preferences())


        def save_dialog_func():
            post_model_fit_refine_dialog()
            post_go_to_atom_window()

            def delete_event(*args):
                window.destroy()
                return False

            def go_func(*args):
                save_dialog_positions_to_preferences_file()
                window.destroy()

            window = Gtk.Window(gtk.WINDOW_TOPLEVEL)
            label_text = "   When happy, press \"Save\" to save   \n" + "   dialog positions"
            label = Gtk.Label(label_text)
            h_sep = Gtk.HSeparator()
            cancel_button = Gtk.Button("  Cancel  ")
            go_button = Gtk.Button("  Save  ")
            vbox = Gtk.VBox(False, 4)
            hbox = Gtk.HBox(False, 4)

            hbox.pack_start(go_button,     False, False, 6)
            hbox.pack_start(cancel_button, False, False, 6)
            vbox.pack_start(label, False, False, 6)
            vbox.pack_start(h_sep, False, False, 6)
            vbox.pack_start(hbox,  False, False, 6)
            window.add(vbox)

            go_button.connect("clicked", go_func)
            cancel_button.connect("clicked", lambda w: window.destroy())

            window.show_all()

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_settings, "Save Dialog Positions...",
            lambda func: save_dialog_func())


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_settings, "Key Bindings...", lambda func: coot_gui.key_bindings_gui())

        def install_and_show_key_bindings():
            coot_utils.file_to_preferences("template_key_bindings.py") # copy and evaluate
            coot_gui.key_bindings_gui()

            coot_gui.add_simple_coot_menu_menuitem(
                submenu_settings, "Python: Install Template Keybindings",
                lambda func: install_and_show_key_bindings())

        def quick_save_func(txt):
            try:
                n = int(txt)
                gobject.timeout_add(1000*n, quick_save)
            except:
                print("BL INFO:: could not add timer for auto save!")

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_settings, "Enable Quick-Save checkpointing...",
            lambda func: coot_gui.generic_single_entry("Checkpoint interval (seconds)",
                                                       "30",
                                                       " Start Auto-saving ",
                                                       lambda txt:
                                                       quick_save_func(txt)))

        def fit_protein_func1(imol):
            if coot.imol_refinement_map() == -1:
                add_status_bar_text("oops. Must set a map to fit")
            else:
                global continue_multi_refine
                continue_multi_refine = True
                fitting.interruptible_fit_protein(imol, fitting.fit_protein_fit_function)

        def fit_protein_func2(imol):
            if coot.imol_refinement_map() == -1:
                add_status_bar_text("oops. Must set a map to fit")
            else:
                global continue_multi_refine
                continue_multi_refine = True
                fitting.interruptible_fit_protein(imol, fitting.fit_protein_stepped_refine_function)

        def fit_protein_func3(imol):
            if coot.imol_refinement_map() == -1:
                add_status_bar_text("oops. Must set a map to fit")
            else:
                global continue_multi_refine
                continue_multi_refine = True
                fitting.interruptible_fit_protein(imol, fitting.fit_protein_rama_fit_function)

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_all_molecule,
            "[Post MR] Fill Partial Residues...",
            lambda func: coot_gui.molecule_chooser_gui("Find and Fill residues with missing atoms",
                                                       lambda imol: coot.fill_partial_residues(imol)))

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_all_molecule,
            "Fit Protein...",
            lambda func: coot_gui.molecule_chooser_gui("Fit Protein using Rotamer Search",
                                                       lambda imol: fit_protein_func1(imol)))

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_all_molecule,
            "Stepped Refine...",
            lambda func: coot_gui.molecule_chooser_gui("Fit Protein using Real-Space Refinement",
                                                       lambda imol: fit_protein_func2(imol)))

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_all_molecule,
            "Refine/Improve Ramachandran Plot...",
            lambda func: coot_gui.molecule_chooser_gui("Refine Protein with Ramachanran Plot Optimization: ",
                                                       lambda imol: fit_protein_func3(imol)))

        # --- add_edit_settings_menu() ends here

    def add_calculate_maps_menu(submenu_maps):

        #---------------------------------------------------------------------
        # ----    Map Tools...
        #---------------------------------------------------------------------

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_maps,
            "Mask Map by Atom Selection...",
            lambda func: coot_gui.molecule_chooser_gui("Define the molecule that has atoms to mask the map",
                                                       lambda imol: coot_gui.generic_multiple_entries_with_check_button(
                                                           [[" Map molecule number: ", mask_map_func()],
                                                            [" Atom selection: ", "//A/1"],
                                                            ["Radius around atoms: ", mask_map_radius_func()]],
                                                           [" Invert Masking? ", lambda active_state: mask_map_func1(active_state)],
                                                           "  Mask Map  ", lambda text_list, invert: mask_map_func2(imol, text_list, invert))))

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_maps,
            "Copy Map...",
            lambda func: coot_gui.map_molecule_chooser_gui("Map to Copy...",
                                                           lambda imol: copy_molecule(imol)))


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_maps,
            "Make a Smoother Copy...",
            lambda func: coot_gui.map_molecule_chooser_gui("Map Molecule to Smoothenize...",
                                                           lambda imol: smooth_map(imol, 1.25)))


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_maps,
            "Make a Very Smooth Copy...",
            lambda func: coot_gui.map_molecule_chooser_gui("Map Molecule to Smoothenize...",
                                                           lambda imol: smooth_map(imol, 2.0)))


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_maps,
            "Make a Difference Map...",
            lambda func: coot_gui.make_difference_map_gui())

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_maps,
            "Transform map by LSQ model fit...",
            lambda func: coot_gui.transform_map_using_lsq_matrix_gui())


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_maps,
            "Average Maps...",
            lambda func: coot_gui.average_map_gui())


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_maps,
            "Map Density Histogram...",
            lambda func: coot_gui.map_molecule_chooser_gui("Choose the map",
                                                           lambda imol: map_histogram(imol)))

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_maps,
            "Brighten Maps",
            lambda func: coot_utils.brighten_maps())

        def set_diff_map_func(imol):
            print("setting map number %s to be a difference map" %imol)
            set_map_is_difference_map(imol)

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_maps,
            "Set map is a difference map...",
            # do we need this wrapper?
            lambda func: coot_gui.map_molecule_chooser_gui("Which map should be considered a difference map?",
                                                           lambda imol: set_diff_map_func(imol)))

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_maps,
            "Another (contour) level...",
            lambda func: another_level())


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_maps,
            "Multi-chicken...",
            lambda func: coot_gui.map_molecule_chooser_gui("Choose a molecule for multiple contouring",
                                                           lambda imol: (set_map_displayed(imol, 0), coot_utils.multi_chicken(imol))))




    def add_calculate_modelling_menu(submenu_maps):

        #---------------------------------------------------------------------
        #     Molecule functions/Modelling
        #
        #---------------------------------------------------------------------

        def add_hydrogens_with_coot_reduce():
            with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                       aa_ins_code, aa_atom_name, aa_alt_conf]:
                coot_reduce(aa_imol)

        def add_hydrogens_refmac_func():
            with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                       aa_ins_code, aa_atom_name, aa_alt_conf]:
                coot_utils.add_hydrogens_using_refmac(aa_imol)

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Add Hydrogens",
            lambda func: add_hydrogens_with_coot_reduce())

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Add Hydrogens using Refmac",
            lambda func: add_hydrogens_refmac_func())

        coot_gui.add_simple_coot_menu_menuitem(
        submenu_modelling,
            "Add Other Solvent Molecules...",
            lambda func: coot_gui.solvent_ligands_gui())

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Arrange Waters Around Protein...",
            lambda func: coot_gui.molecule_chooser_gui(
                "Arrange waters in molecule: ",
                lambda imol: move_waters_to_around_protein(imol)))

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Assign (force) HETATMs for this Residue",
            lambda func: coot_utils.using_active_atom(hetify_residue,
                                                      "aa_imol", "aa_chain_id", "aa_res_no", "aa_ins_code"))


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Assign HETATM to molecule...",
            lambda func: coot_gui.molecule_chooser_gui("Assign HETATMs as per PDB definition",
                                                       lambda imol: assign_hetatms(imol)))

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Atoms with Zero Occupancies...",
            lambda func: coot_gui.molecule_chooser_gui(
                "Which molecule to check for Atoms with zero occupancies?",
                lambda imol: coot_gui.zero_occ_atoms_gui(imol)))

        # --- D --------

        def delete_sidechains_for_active_chain():
            with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code, aa_atom_name, aa_alt_conf]:
                coot.delete_sidechains_for_chain(aa_imol, aa_chain_id)

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Delete Side-chain Atoms for Active Chain",
            lambda x: delete_sidechains_for_active_chain())

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Duplicate range (pick atoms)",
            lambda func: coot_gui.duplicate_range_by_atom_pick())


        # --- F --------

        def get_smiles_pdbe_func():
            with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code, aa_atom_name, aa_alt_conf]:
                comp_id = residue_name(aa_imol, aa_chain_id,
                                       aa_res_no, aa_ins_code)
                coot_utils.get_SMILES_for_comp_id_from_pdbe(comp_id)

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Fetch PDBe description for this ligand",
            lambda func: get_smiles_pdbe_func())

        def get_pdbe_ligand_func(comp_id):
            status = coot_utils.get_SMILES_for_comp_id_from_pdbe(comp_id)
            get_monomer(comp_id)

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Fetch PDBe Ligand Description",
            lambda func: coot_gui.generic_single_entry("Fetch PDBe Ligand Desciption for comp_id:",
                                                       "", " Fetch ", lambda comp_id: get_pdbe_ligand_func(comp_id)))

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Fix Nomenclature Errors...",
            lambda func: coot_gui.molecule_chooser_gui("Fix Nomenclature Error in molecule:",
                                                       lambda imol: fix_nomenclature_errors(imol)))


        # --- I --------

        def chiral_centre_inverter_func():
            with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                       aa_ins_code, aa_atom_name, aa_alt_conf]:
                invert_chiral_centre(aa_imol, aa_chain_id, aa_res_no,
                                     aa_ins_code, aa_atom_name)

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Invert This Chiral Centre",
            lambda func:
            chiral_centre_inverter_func())


        # --- J --------

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "JLigand launch",
            lambda func: jligand_gui.launch_jligand_function())


        # --- M ---

        def make_link_ext_func(*args):
            m_spec_1 = args[0]
            m_spec_2 = args[1]
            imol_1 = coot_utils.atom_spec_to_imol(m_spec_1)
            imol_2 = coot_utils.atom_spec_to_imol(m_spec_2)
            spec_1 = m_spec_1[2:]
            spec_2 = m_spec_2[2:]
            if not (imol_1 == imol_2):
                print("Mismatch molecules")
            else:
                make_link(imol_1, spec_1, spec_2, "dummy", 0.1)

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Make Link (click 2 atoms)...",
            lambda func:
            user_defined_click(2, make_link_ext_func))


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Merge Water Chains...",
            lambda func: coot_gui.molecule_chooser_gui("Merge Water Chains in molecule:",
                                                       lambda imol: coot_utils.merge_solvent_chains(imol)))


        def mon_dict_func(text):
            idealized = 0
            new_model = get_monomer_from_dictionary(text, idealized)
            if not coot_utils.valid_model_molecule_qm(new_model):
                get_monomer(text)

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Monomer from Dictionary",
            lambda func:
            coot_gui.generic_single_entry("Pull coordinates from CIF dictionary for 3-letter-code:", "",
                                          " Get Coords ",
                                          lambda text: mon_dict_func(text)))


        # -- N --

        def new_mol_sphere_func1(imol, text):
            try:
                radius = float(text)
                args = [imol] + coot_utils.rotation_centre() + [radius, 0]
                new_molecule_by_sphere_selection(*args)
            except:
                print("WARNING:: no valid radius", text)

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "New Molecule by Sphere...",
            lambda func: coot_gui.generic_chooser_and_entry(
                "Choose a molecule from which to select a sphere of atoms:",
                "Radius:", "10.0",
                lambda imol, text: new_mol_sphere_func1(imol, text)))

        def new_mol_sym_func1(imol, text):
            from types import ListType
            pre_shift = origin_pre_shift(imol)
            if (type(pre_shift) is not ListType):
                print("bad pre-shift aborting")
            else:
                new_molecule_by_symop(imol, text,
                                      pre_shift[0],
                                      pre_shift[1],
                                      pre_shift[2])

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "New Molecule from Symmetry Op...",
            lambda func: coot_gui.generic_chooser_and_entry(
                "Molecule from which to generate a symmetry copy",
                "SymOp", "X,Y,Z",
                lambda imol, text: new_mol_sym_func1(imol, text)))


        # --- P ---

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Phosphorylate this residue",
            lambda func: coot_utils.phosphorylate_active_residue())


        # ---- R ---------

        # --- Ren ---

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Rename Residue...",
            lambda func: coot_gui.rename_residue_gui())


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Renumber Waters...",
            lambda func: coot_gui.molecule_chooser_gui(
                "Renumber waters of which molecule?",
                lambda imol: renumber_waters(imol)))

        # --- Reo ---

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Reorder Chains...",
            lambda func: coot_gui.molecule_chooser_gui("Sort Chain IDs in molecule:",
                                                       lambda imol: sort_chains(imol))) # an internal function

        # --- Res ---

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Residue Type Selection...",
            lambda func: coot_gui.generic_chooser_and_entry("Choose a molecule to select residues from: ","Residue Type:","",
                                                            lambda imol, text: (new_molecule_by_residue_type_selection(imol, text),
                                                                                update_go_to_atom_window_on_new_mol())))


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Residues with Alt Confs...",
            lambda func: coot_gui.molecule_chooser_gui(
                "Which molecule to check for Alt Confs?",
                lambda imol: coot_gui.alt_confs_gui(imol)))


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Residues with Cis Peptide Bonds...",
            lambda func: coot_gui.molecule_chooser_gui("Choose a molecule for checking for Cis Peptides",
                                                       lambda imol: coot_gui.cis_peptides_gui(imol)))


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Residues with Missing Atoms...",
            lambda func: coot_gui.molecule_chooser_gui(
                "Which molecule to check for Missing Atoms?",
                lambda imol: coot_gui.missing_atoms_gui(imol)))

        # --- Rig ---

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Rigid Body Fit Residue Ranges...",
            lambda func:
            coot_gui.residue_range_gui(lambda imol, ls: rigid_body_refine_by_residue_ranges(imol, ls),
                                       "Rigid Body Refine",
                                       "  Fit  "))

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Rigid Body Fit Molecule...",
            lambda func: coot_gui.molecule_chooser_gui("Rigid Body Fit Molecule",
                                                       lambda imol: rigid_body_refine_by_atom_selection(imol, "//")))

        # ---- S --------

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Superpose ligands",
            lambda func: coot_gui.superpose_ligand_gui())


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "Symm Shift Reference Chain Here",
            lambda func: move_reference_chain_to_symm_chain_position())

        # ---- W ---------

        def whats_this():
            central_residue = active_residue()
            res_name = residue_name(*central_residue[0:4])
            mol_no = central_residue[0]
            n = comp_id2name(res_name)
            s = "(mol. no: " + str(mol_no) + ")  " + \
                res_name  + ":  " + \
                n if isinstance(n, str) else " <no-name-found>"
            add_status_bar_text(s)

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_modelling,
            "What's this?",
            lambda func: whats_this())


    def add_calculate_ncs_menu(submenu_ncs):

    #---------------------------------------------------------------------
    #     NCS functions
    #
    #---------------------------------------------------------------------

        def copy_ncs_range_func(imol, chain_id, text1, text2):
            try:
                r1 = int(text1)
                r2 = int(text2)
                if (chain_id == ncs.ncs_master_chain_id(imol)):
                    copy_residue_range_from_ncs_master_to_others(imol, chain_id, r1, r2)
                else:
                    # different given master to current master
                    # ask what to do.
                    txt = "Current master chain is %s, but you asked to copy from %s.\n" \
                        %(ncs_master_chain_id(imol), chain_id)
                    txt += "Apply this master change?\n\n"
                    txt += "N.B. if no, then nothing is copied."
                    r = coot_gui.yes_no_dialog(txt, "Change Master")
                    if r:
                        ncs_control_change_ncs_master_to_chain_id(imol, chain_id)
                        copy_residue_range_from_ncs_master_to_others(imol, chain_id, r1, r2)
                        ## could change master back?!
                    else:
                        info_dialog("Master chain was not changed and copy not applied.")
            except:
                print("WARNING:: no valid number input")

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_ncs,
            "Copy NCS Residue Range...",
            lambda func: coot_gui.generic_chooser_and_entry("Apply NCS Range from Master",
                                                            "Master Chain ID",
                                                            coot_utils.get_first_ncs_master_chain(),  # returns "" on fail
                                                            lambda imol, chain_id: coot_gui.generic_double_entry("Start of Residue Number Range",
                                                                                                                 "End of Residue Number Range",
                                                                                                                 "", "", False, False,
                                                                                                                 "Apply NCS Residue Range",
                                                                                                                 lambda text1, text2: copy_ncs_range_func(imol, chain_id, text1, text2))))


        def copy_ncs_chain_func(imol, chain_id):
            ncs_chains = ncs_chain_ids(imol)
            if (ncs_chains):
                # maybe this could be a function to avoid repetition.
                if (chain_id == ncs.ncs_master_chain_id(imol)):
                    # hunkey dorey
                    copy_from_ncs_master_to_others(imol, chain_id)
                else:
                    # different given master to current master
                    # ask what to do.
                    txt = "Current master chain is %s, but you asked to copy from %s.\n" \
                        %(ncs_master_chain_id(imol), chain_id)
                    txt += "Apply this master change?\n\n"
                    txt += "N.B. if no, then nothing is copied."
                    r = coot_gui.yes_no_dialog(txt, "Change Master")
                    if r:
                        ncs_control_change_ncs_master_to_chain_id(imol, chain_id)
                        copy_from_ncs_master_to_others(imol, chain_id)
                        ## could change master back?!
                    else:
                        info_dialog("Master chain was not changed and copy not applied.")
            else:
                s = "You need to define NCS operators for molecule " + str(imol)
                info_dialog(s)

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_ncs,
            "Copy NCS Chain...",
            lambda func: coot_gui.generic_chooser_and_entry("Apply NCS edits from NCS Master Chain to Other Chains",
                                                            "Master Chain ID",
                                                            coot_utils.get_first_ncs_master_chain(),  # can return  "".
                                                            lambda imol, chain_id: copy_ncs_chain_func(imol, chain_id)))


        def ncs_ghost_res_range_func(imol):

            def handle_go_function(resno_1_text, resno_2_text):
                cont = False
                try:
                    resno_1 = int(resno_1_text)
                    resno_2 = int(resno_2_text)
                    cont = True
                except:
                    print("WARNING:: input residue numbers have to be integers")
                if cont:
                    ghost_ncs_chain_ids = ncs_chain_ids(imol)
                    if (type(ghost_ncs_chain_ids) is ListType):
                        # because we can have hetero-NCS,
                        # but we ignore NCS other that
                        # that of the first type.
                        ghost_chain_list = ghost_ncs_chain_ids[0]
                        ncs.manual_ncs_ghosts(imol, resno_1, resno_2, ghost_chain_list)

            from types import ListType
            label_1 = "Start ResNumber"
            label_2 = "End ResNumber"
            entry_1_default_text = "10"
            entry_2_default_text = "20"
            go_button_label = "Regenerate Ghosts"

            coot_gui.generic_double_entry(label_1, label_2, entry_1_default_text,
                                          entry_2_default_text, False, False,
                                          go_button_label, handle_go_function)


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_ncs,
            "NCS Ghosts by Residue Range...",
            lambda func: coot_gui.molecule_chooser_gui("Make local NCS ghosts for molecule:",
                                                       lambda imol: ncs_ghost_res_range_func(imol)))

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_ncs,
            "Update NCS Ghosts using Local Match",
            lambda func: ncs.update_ncs_ghosts_by_local_sphere())


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_ncs,
            "NCS Jumping...",
            lambda func: coot_gui.ncs_jumping_gui())


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_ncs,
            "NCS ligands...",
            lambda func: coot_gui.ncs_ligand_gui())

        submenu = Gtk.Menu()
        menuitem2 = Gtk.MenuItem("NCS matrix type...")

        menuitem2.set_submenu(submenu)
        submenu_ncs.append(menuitem2)
        menuitem2.show()

        coot_gui.add_simple_coot_menu_menuitem(
            submenu,
            "Accurate (SSM)",
            lambda func: set_ncs_matrix_type(0))


        coot_gui.add_simple_coot_menu_menuitem(
            submenu,
            "Fast (LSQ)",
            lambda func: set_ncs_matrix_type(1))


        coot_gui.add_simple_coot_menu_menuitem(
            submenu,
            "Extra Fast (LSQ, only every 2nd CA)",
            lambda func: set_ncs_matrix_type(2))




    def add_draw_representations_menu(submenu_views):

        # ---------------------------------------------------------------------
        #     Views/Representations
        # ---------------------------------------------------------------------
        #

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_representation,
            "Undo Symmetry View",
            lambda func: undo_symmetry_view())


        def make_ball_n_stick_func(imol, text):
            bns_handle = make_ball_and_stick(imol, text, 0.18, 0.3, 1)
            print("handle: ", bns_handle)


        global default_ball_and_stick_selection    # maybe should be at the top of the file

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_representation,
            "Ball & Stick...",
            lambda func: coot_gui.generic_chooser_and_entry("Ball & Stick",
                                                            "Atom Selection:",
                                                            default_ball_and_stick_selection,
                                                            lambda imol, text: make_ball_n_stick_func(imol, text)))

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_representation,
            "Add Balls to Simple Sticks",
            lambda func: [set_draw_stick_mode_atoms(imol, 1) for imol in coot_utils.molecule_number_list()])


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_representation,
            "Simple Sticks (No Balls)",
            lambda func: [set_draw_stick_mode_atoms(imol, 0) for imol in coot_utils.molecule_number_list()])


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_representation,
            "Clear Ball & Stick...",
            lambda func: coot_gui.molecule_chooser_gui(
                "Choose a molecule from which to clear Ball&Stick objects",
                lambda imol: clear_ball_and_stick(imol)))


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_representation,
            "Electrostatic Surface...",
            lambda func: coot_gui.molecule_chooser_gui(
                "Choose a molecule to represent as a surface..." + \
                "\n" + \
                "Can be SLOW",
                lambda imol: do_surface(imol, 1)))  # shall we switch on the light too?!


        def surface_func1(clipped = 0):
            active_atom = active_residue()
            aa_imol      = active_atom[0]
            aa_chain_id  = active_atom[1]
            aa_res_no    = active_atom[2]
            aa_ins_code  = active_atom[3]
            aa_atom_name = active_atom[4]
            aa_alt_conf  = active_atom[5]
            central_residue = active_residue()
            residues = residues_near_residue(aa_imol, central_residue[1:4], 6.0)
            # no waters in surface, thanks.
            # but what if they have different names?! (only HOH so far)
            filtered_residues = []
            for res in residues:
                if (residue_name(aa_imol, *res) != "HOH"):
                    filtered_residues.append(res)
            imol_copy = copy_molecule(aa_imol)
            # delete the interesting residue from the copy (so that
            # it is not surfaced).
            delete_residue(imol_copy, aa_chain_id, aa_res_no, aa_ins_code)
            if clipped:
                do_clipped_surface(imol_copy, filtered_residues)
            else:
                do_surface(imol_copy, 1)

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_representation,
            "Clipped Surface Here (This Residue)",
            lambda func:
            surface_func1(1))

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_representation,
            "Full Surface Around Here (This Residue)",
            lambda func:
            surface_func1())


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_representation,
            "Un-Surface...",
            lambda func: coot_gui.molecule_chooser_gui(
                "Choose a molecule to represent conventionally...",
                lambda imol: do_surface(imol, 0)))

        def hilight_site_func():
            active_atom = active_residue()
            if (active_atom):
                imol = active_atom[0]
                centre_residue_spec = [active_atom[1],
                                       active_atom[2],
                                       active_atom[3]]
                coot_utils.hilight_binding_site(imol, centre_residue_spec, 230,4)

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_representation,
            "Highlight Interesting Site (here)...",
            lambda func: hilight_site_func())


        def make_dot_surf_func(imol,text):
            # I think a single colour is better than colour by atom
            set_dots_colour(imol, 0.5, 0.5, 0.5)
            dots_handle = dots(imol, text, text, 2, 1)
            print("dots handle: ", dots_handle)

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_representation,
            "Dotted Surface...",
            lambda func: coot_gui.generic_chooser_and_entry("Surface for molecule",
                                                            "Atom Selection:", "//A/1-2",
                                                            lambda imol, text: make_dot_surf_func(imol, text)))


        def clear_dot_surf_func(imol,text):
            try:
                n = int(text)
                clear_dots(imol,n)
            except:
              print("BL WARNING:: dots handle number shall be an integer!!")

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_representation,
            "Clear Surface Dots...",
            lambda func: coot_gui.generic_chooser_and_entry("Molecule with Dotted Surface",
                                                            "Dots Handle Number:", "0",
                                                            lambda imol, text: clear_dot_surf_func(imol, text)))


        def limit_model_disp_func(text):
            try:
                f = float(text)
                if f < 0.1:
                    set_model_display_radius(0, 10)
                else:
                    set_model_display_radius(1, f)
            except:
                set_model_display_radius(0, 10)

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_representation,
            "Limit Model Display Radius...",
            lambda func: coot_gui.generic_single_entry("Display Radius Limit (0 for \'no limit\') ",
                                                       #  "15.0" ;; maybe this should be the map radius
                                                       # BL says:: I think it should be the current one
                                                       str(get_map_radius()),
                                                       "Set: ",
                                                       lambda text: limit_model_disp_func(text)))


        coot_gui.add_simple_coot_menu_menuitem(
            submenu_representation,
            "HOLE...",
            lambda func: test_hole.hole_ify())

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_representation,
            "Label All CAs...",
            lambda func: coot_gui.molecule_chooser_gui("Choose a molecule to label",
                                                       lambda imol: coot_utils.label_all_CAs(imol)))



        #---------------------------------------------------------------------
        #     3D annotations
        #---------------------------------------------------------------------

        submenu = Gtk.Menu()
        menuitem2 = Gtk.MenuItem("3D Annotations...")

        menuitem2.set_submenu(submenu)
        submenu_representation.append(menuitem2)
        menuitem2.show()

        coot_gui.add_simple_coot_menu_menuitem(
            submenu,
            "Annotate position...",
            lambda func: coot_gui.generic_single_entry("Annotation: ", "",
                                                       "Make Annotation",
                                                       lambda txt: coot_utils.add_annotation_here(txt)))


        coot_gui.add_simple_coot_menu_menuitem(
            submenu,
            "Save Annotations...",
            lambda func: coot_gui.generic_single_entry("Save Annotations",
                                                       "coot_annotations.py",
                                                       " Save ",
                                                       lambda file_name: coot_utils.save_annotations(file_name)))


        coot_gui.add_simple_coot_menu_menuitem(
            submenu,
            "Load Annotations...",
            lambda func: coot_gui.generic_single_entry("Load Annotations",
                                                       "coot_annotations.py",
                                                       " Load ",
                                                       lambda file_name: coot_utils.load_annotations(file_name)))


        coot_gui.add_simple_coot_menu_menuitem(
            submenu,
            "Remove annotation here",
            lambda func: coot_utils.remove_annotation_here())


        coot_gui.add_simple_coot_menu_menuitem(
            submenu,
            "Remove annotation near click",
            lambda func: coot_utils.remove_annotation_at_click())



    def add_draw_views_menu(submenu_views):

           coot_gui.add_simple_coot_menu_menuitem(
               submenu_views,
               "Add View...",
               lambda func: coot_gui.view_saver_gui())

           # BL says:: maybe check if number at some point
           coot_gui.add_simple_coot_menu_menuitem(
               submenu_views,
               "Add a Spin View...",
               lambda func: coot_gui.generic_double_entry("Number of Steps", 
                                                          "Number of Degrees (total)", "3600", "360", 
                                                          False, False,  #check button text and callback
                                                          "  Add Spin  ",
                                                          lambda text_1, text_2: add_spin_view("Spin", int(text_1), float(text_2))))

           coot_gui.add_simple_coot_menu_menuitem(
               submenu_views,
               "Views Panel...",
               lambda func: coot_gui.views_panel_gui())

           coot_gui.add_simple_coot_menu_menuitem(
               submenu_views,
               "Play Views",
               lambda func: list(map(eval,["go_to_first_view(1)",
                                           "time.sleep(1)", "play_views()"])))

           # BL says:: maybe check if number at some point
           coot_gui.add_simple_coot_menu_menuitem(
               submenu_views, "Set Views Play Speed...",
               lambda func: coot_gui.generic_single_entry("Set Views Play Speed",
                                                          str(views_play_speed()), "  Set it  ",
                                                          lambda text: set_views_play_speed(float(text))))


           coot_gui.add_simple_coot_menu_menuitem(
               submenu_views, "Save Views...",
               lambda func: coot_gui.generic_single_entry("Save Views",
                                                          "coot-views.py", " Save ",
                                                          lambda txt: save_views(txt)))


    def add_calculate_pisa_menu(submenu_pisa):
        # ---------------------------------------------------------------------
        #     PISA Interface and Assemblies
        # ---------------------------------------------------------------------

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_pisa, "PISA Assemblies...",
            lambda func:
            coot_gui.molecule_chooser_gui("Choose molecule for PISA assembly analysis",
                                          lambda imol:
                                          parse_pisa_xml.pisa_assemblies(imol)))

        coot_gui.add_simple_coot_menu_menuitem(
            submenu_pisa, "PISA Interfaces...",
            lambda func:
            coot_gui.molecule_chooser_gui("Choose molecule for PISA interface analysis",
                                          lambda imol:
                                          parse_pisa_xml.pisa_interfaces(imol)))



    #---------------------------------------------------------------------
    #  Run those functions to populate the menus
    #---------------------------------------------------------------------

    add_calculate_modules_menu(submenu_modules)
    add_calculate_maps_menu(submenu_maps)
    add_calculate_modelling_menu(submenu_modelling)
    add_calculate_ncs_menu(submenu_ncs)
    add_edit_settings_menu(submenu_settings)
    add_draw_representations_menu(submenu_representation)
    add_draw_views_menu(submenu_views)
    add_calculate_pisa_menu(submenu_pisa)
