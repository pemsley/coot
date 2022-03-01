# extensions.py
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
import gui_add_linked_cho
import gui_prosmart
import shelx_extensions


def add_coot_menu_separator(menu):
    sep = Gtk.MenuItem()
    menu.add(sep)
    #   sep.props.sensitive = False
    sep.show()
   
# if have_coot_python: # how do I get to have_coot_python now?
#    if coot_python.main_menubar():

if True:

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
             # print("########### extensions get_existing_submenu get_text on c:", t)
             if t == submenu_label:
               return menu_child
           except KeyError as e:
             pass
           
       return False
       
     # --------------------------------------------------
     #           Validation
     # --------------------------------------------------

     menu = coot_gui.coot_menubar_menu("Validate")
     if menu:
         import generic_objects
         import dynamic_atom_overlaps_and_other_outliers as dao
         coot_gui.add_simple_coot_menu_menuitem(menu, "Highly coordinated waters...",
                                                lambda func: coot_gui.water_coordination_gui())
         coot_gui.add_simple_coot_menu_menuitem(menu, "Atom Overlaps (Coot)",
                                                lambda func: coot_utils.using_active_atom(coot.coot_all_atom_contact_dots, "aa_imol"))
         coot_gui.add_simple_coot_menu_menuitem(menu, "All-Atom Contact Dots (Molprobity)",
                                                lambda func: coot_utils.using_active_atom(generic_objects.probe, "aa_imol"))
         #coot_gui.add_simple_coot_menu_menuitem(menu,"Atom Overlaps Dialog",
         #                                      lambda func: coot_utils.using_active_atom(dao.make_quick_test_validation_buttons, "aa_imol"))
         coot_gui.add_simple_coot_menu_menuitem(menu, "Highly coordinated waters...",
                                                lambda func: coot_gui.water_coordination_gui())

       # this does not exist. Let's do in in C++ with glade
       #coot_gui.add_simple_coot_menu_menuitem(menu, "Pepflips from Difference Map...",
       #                               lambda func: coot_gui.pepflips_by_difference_map_gui())


         def validation_outliers_func():
             with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                    aa_ins_code, aa_atom_name, aa_alt_conf]:
                 imol_map = coot.imol_refinement_map()
                 if not coot_utils.valid_map_molecule_qm(imol_map):
                     coot.info_dialog_and_text("Refinement Map is currently not set")
                 else:
                     find_baddies.validation_outliers_dialog(aa_imol, imol_map)

         coot_gui.add_simple_coot_menu_menuitem(menu, "Validation Outliers",
                                                lambda func: validation_outliers_func())

         coot_gui.add_simple_coot_menu_menuitem(menu, "List Ramachandran outliers...",
                                                lambda func: rama_outlier_gui())


     # --------------------------------------------------
     #           user_define_restraints plugin
     # --------------------------------------------------

     def add_module_user_defined_restraints():
         import user_define_restraints
         menu = coot_gui.coot_menubar_menu("Restraints")
     
     # ---------------------------------------------
     #           extensions
     # ---------------------------------------------

     # menu = coot_gui.coot_menubar_menu("Extensions")

     # calculate_menu = menu
     # draw_menu = menu
     # edit_menu = menu

     calculate_menu = coot_gui.coot_menubar_menu("Calculate")
     draw_menu      = coot_gui.coot_menubar_menu("Draw")
     edit_menu      = coot_gui.coot_menubar_menu("Edit")
     edit_settings_menu = get_existing_submenu(edit_menu, "Settings...")
     calculate_all_molecule_menu = get_existing_submenu(calculate_menu, "All Molecule...")
     calculate_map_tools_menu    = get_existing_submenu(calculate_menu, "Map Tools...")
     calculate_ncs_tools_menu    = get_existing_submenu(calculate_menu, "NCS Model Tools...")
     calculate_ncs_maps_menu     = get_existing_submenu(calculate_menu, "NCS Maps...")
     calculate_modelling_menu    = get_existing_submenu(calculate_menu, "Modelling...")
     calculate_modules_menu      = get_existing_submenu(calculate_menu, "Modules...")
     draw_representation_menu    = get_existing_submenu(draw_menu, "Representation Tools...")

     # make submenus:
     submenu_all_molecule = Gtk.Menu()

     # print("DEBUG::::::::::::::::::::: edit_settings_menu:", edit_settings_menu)
     # print("DEBUG::::::::::::::::::::: calculate_all_molecule_menu:", calculate_menu)
     # print("DEBUG::::::::::::::::::::: submenu_all_molecule:", submenu_all_molecule)

     # menuitem_2 = Gtk.MenuItem("All Molecule...")
     submenu_maps = Gtk.Menu()
     menuitem_3 = Gtk.MenuItem("Maps...")
     submenu_models = Gtk.Menu() # submenu_modelling really
     # menuitem_4 = Gtk.MenuItem("Modelling...")
     submenu_refine = Gtk.Menu()
     menuitem_5 = Gtk.MenuItem("Refine...")
     submenu_representation = Gtk.Menu()
     menuitem_6 = Gtk.MenuItem("Representations")
     submenu_settings = Gtk.Menu()
     menuitem_7 = Gtk.MenuItem("Settings...")
     submenu_pisa = Gtk.Menu()
     menuitem_pisa = Gtk.MenuItem("PISA...")
     submenu_pdbe = Gtk.Menu()
     menuitem_pdbe = Gtk.MenuItem("PDBe...")
     submenu_modules = Gtk.Menu()
     # menuitem_modules = Gtk.MenuItem("Modules...")
     submenu_ncs = Gtk.Menu()
     # menuitem_ncs = Gtk.MenuItem("NCS Tools...") # find existing...

     # menuitem_2.set_submenu(submenu_all_molecule) replace by the following
     calculate_all_molecule_menu.set_submenu(submenu_all_molecule)
     # calculate_menu.append(menuitem_2)
     # menuitem_2.show()

     # menuitem_3.set_submenu(submenu_maps)
     # calculate_menu.append(menuitem_3)
     # menuitem_3.show()

     calculate_map_tools_menu.set_submenu(submenu_maps)
     calculate_ncs_tools_menu.set_submenu(submenu_ncs)

     # menuitem_4.set_submenu(submenu_models)
     #calculate_menu.append(menuitem_4)
     # menuitem_4.show()

     calculate_modelling_menu.set_submenu(submenu_models)

     # menuitem_ncs.set_submenu(submenu_ncs) 20210928-PE
     # calculate_menu.append(menuitem_ncs)
     # menuitem_ncs.show()

     calculate_ncs_tools_menu.set_submenu(submenu_ncs)

     # menuitem_6.set_submenu(submenu_representation)
     # draw_menu.append(menuitem_6)
     # menuitem_6.show()

     draw_representation_menu.set_submenu(submenu_representation)

     menuitem_pisa.set_submenu(submenu_pisa)
     draw_menu.append(menuitem_pisa)
     menuitem_pisa.show()

     # menuitem_7.set_submenu(submenu_settings)
     edit_settings_menu.set_submenu(submenu_settings)
     menu.append(menuitem_7)
     menuitem_7.show()

     calculate_modules_menu.set_submenu(submenu_modules)

     # menuitem_modules.set_submenu(submenu_modules)
     # calculate_menu.append(menuitem_modules)
     # menuitem_modules.show()

     # menuitem_7.set_submenu(submenu_settings) # already set
     # edit_menu.append(menuitem_7)
     # menuitem_7.show()

     # where does the Refine submenu go? In Edit -> Settings
     edit_settings_submenu = menuitem_7
     # menuitem_5.set_submenu(submenu_refine)
     # edit_settings_submenu.append(menuitem_5)
     # menuitem_5.show()

     # give edit_settings_menu a submenu
     # submenu_settings = Gtk.Menu()
     # edit_settings_menu.set_submenu(submenu_settings)
     # submenu_settings.show()
     # edit_settings_submenu = edit_settings_menu.get_submenu()
     # print("###### debug edit_settings_submenu", edit_settings_submenu)


     
     # menuitem_pdbe.set_submenu(submenu_pdbe)
     # menu.append(menuitem_pdbe)
     # menuitem_pdbe.show()
     

     #---------------------------------------------------------------------
     #     Post MR
     #
     #---------------------------------------------------------------------

     coot_gui.add_simple_coot_menu_menuitem(
       submenu_all_molecule,
       "[Post MR] Fill Partial Residues...",
       lambda func: coot_gui.molecule_chooser_gui("Find and Fill residues with missing atoms",
		lambda imol: coot.fill_partial_residues(imol)))


     # old style - not interruptable
     #     coot_gui.add_simple_coot_menu_menuitem(
     #       submenu_all_molecule,
     #       "[Post MR] Fit Protein...",
     #       lambda func: coot_gui.molecule_chooser_gui("Fit Protein using Rotamer Search",
     #		lambda imol: (coot.imol_refinement_map() == -1 and
     #                              coot.add_status_bar_text("oops. Must set a map to fit") or
     #                              fitting.fit_protein(imol))))
     #
     #
     #     coot_gui.add_simple_coot_menu_menuitem(
     #       submenu_all_molecule,
     #       "[Post MR] Stepped Refine...",
     #       lambda func: coot_gui.molecule_chooser_gui("Stepped Refine: ",
     #		lambda imol: (coot.imol_refinement_map() == -1 and
     #                              coot.add_status_bar_text("oops. Must set a map to fit") or
     #                              fitting.stepped_refine_protein(imol))))
     #
     #     coot_gui.add_simple_coot_menu_menuitem(
     #       submenu_all_molecule,
     #       "Refine/Improve Ramachandran Plot...",
     #       lambda func: coot_gui.molecule_chooser_gui("Refine Protein with Ramachanran Plot Optimization: ",
     #                lambda imol: (coot.imol_refinement_map() == -1 and
     #                              coot.add_status_bar_text("oops. Must set a map to fit") or
     #                              fitting.stepped_refine_protein_for_rama(imol))))
     #                                         

     # BL says:: we cannot do this with lambda functions in python
     # statements are not allowed! Needed for globals!!
     def fit_protein_func1(imol):
       if (coot.imol_refinement_map() == -1):
         coot.add_status_bar_text("oops. Must set a map to fit")
       else:
         global continue_multi_refine
         continue_multi_refine = True
         interruptible_fitting.fit_protein(imol, fitting.fit_protein_fit_function)

     coot_gui.add_simple_coot_menu_menuitem(
       submenu_all_molecule,
       "Fit Protein...",
       lambda func: coot_gui.molecule_chooser_gui("Fit Protein using Rotamer Search",
		lambda imol: fitting.fit_protein_func1(imol)))


     def fit_protein_func2(imol):
       if (coot.imol_refinement_map() == -1):
         coot.add_status_bar_text("oops. Must set a map to fit")
       else:
         global continue_multi_refine
         continue_multi_refine = True
         interruptible_fitting.fit_protein(imol, fitting.fit_protein_stepped_refine_function)

     coot_gui.add_simple_coot_menu_menuitem(
       submenu_all_molecule,
       "Stepped Refine...",
       lambda func: coot_gui.molecule_chooser_gui("Fit Protein using Real-Space Refinement",
		lambda imol: fitting.fit_protein_func2(imol)))


     def fit_protein_func3(imol):
       if (coot.imol_refinement_map() == -1):
         coot.add_status_bar_text("oops. Must set a map to fit")
       else:
         global continue_multi_refine
         continue_multi_refine = True
         interruptible_fitting.fit_protein(imol, fitting.fit_protein_rama_fit_function)
         
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_all_molecule,
       "Refine/Improve Ramachandran Plot...",
       lambda func: coot_gui.molecule_chooser_gui("Refine Protein with Ramachanran Plot Optimization: ",
                lambda imol: fitting.fit_protein_func3(imol)))
                                         

     #---------------------------------------------------------------------
     #     Map functions
     #
     #---------------------------------------------------------------------

     def mask_map_func():
         f = ""
         molecule_list = coot_utils.molecule_number_list()
         if not molecule_list == []:
             for i in molecule_list:
                 if coot.is_valid_map_molecule(molecule_list[i]):
                     print("%s is a valid map molecule" %molecule_list[i])
                     f = str(molecule_list[i])
                     break
         else:
             print("BL WARNING:: dunno what to do!? No map found")
             f = False
         return f

     def mask_map_func1(active_state):
         print("changed active_state to ", active_state)

     def mask_map_func2(imol, texts_list, invert_mask_qm):
       # map imol
       text_1 = texts_list[0]
       # atom selection
       text_2 = texts_list[1]
       # radius
       text_3 = texts_list[2]

       continue_qm = False
       try:
         n = int(text_1)
         continue_qm = True
       except:
         print("BL WARNING:: input %s for Map molecule number is not an integer.\nBailing out" %n)

       if (continue_qm):
         if (invert_mask_qm):
           invert = 1
         else:
           invert = 0
         print("debug:: invert-mask? is", invert)
         if (text_3 != "default"):
           try:
             new_radius = float(text_3)
             coot.set_map_mask_atom_radius(new_radius)
           except:
             print("BL WARNING:: could not set map mask radius to %s. It's not a number" %text_3)
         coot.mask_map_by_atom_selection(n, imol, text_2, invert)
                
     def mask_map_radius_func():
       rad = coot.map_mask_atom_radius()
       if (rad > 0):
         return str(rad)
       else:
         return "default"
       
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
                                             lambda imol: coot.copy_molecule(imol)))

     
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_maps,
       "Make a Smoother Copy...", 
       lambda func: coot_gui.map_molecule_chooser_gui("Map Molecule to Smoothenize...", 
                                             lambda imol: coot.smooth_map(imol, 1.25)))

     
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_maps,
       "Make a Very Smooth Copy...", 
       lambda func: coot_gui.map_molecule_chooser_gui("Map Molecule to Smoothenize...", 
                                             lambda imol: coot.smooth_map(imol, 2.0)))

     
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

#      coot_gui.add_simple_coot_menu_menuitem(
#        submenu_maps,
#        "Export map...",
#        lambda func: coot_gui.generic_chooser_and_file_selector("Export Map: ",
#                 coot_utils.valid_map_molecule_qm, "File_name: ", "",
#                 lambda imol, text: (coot.export_map(imol, text) == 1 and
#                                     coot.add_status_bar_text("Map " + str(imol) + " exported to " + text))))
                                    

#      def export_map_func(imol, radius_string, file_name):
#        radius = -1.
#        try:
#          radius = float(radius_string)
#        except:
#          print "BL WARNING:: radius %s was no number!" %(radius_string)
#        if (radius >= 0):
#          coot.export_map_fragment(*([imol] + coot_utils.rotation_centre() + \
#                                [radius, file_name]))
         
#      coot_gui.add_simple_coot_menu_menuitem(
#        submenu_maps,
#        "Export Local Map Fragment...", # (for Pymol, say)
#        lambda func: coot_gui.generic_chooser_entry_and_file_selector(
#                "Export Map: ",
#                coot_utils.valid_map_molecule_qm,
#                "Radius (A): ", "10",
#                "File-name: ",
#                lambda imol, radius_string, file_name:
#                    export_map_func(imol, radius_string, file_name)
#                )
#        )


     coot_gui.add_simple_coot_menu_menuitem(
       submenu_maps,
       "Map Density Histogram...",
       lambda func: coot_gui.map_molecule_chooser_gui("Choose the map",
		lambda imol: coot.map_histogram(imol)))


     coot_gui.add_simple_coot_menu_menuitem(
       submenu_maps,
       "Brighten Maps",
       lambda func: coot_utils.brighten_maps())


     def set_diff_map_func(imol):
       print("setting map number %s to be a difference map" %imol)
       coot.set_map_is_difference_map(imol)
        
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_maps,
       "Set map is a difference map...",
       lambda func: coot_gui.map_molecule_chooser_gui("Which map should be considered a difference map?",
		lambda imol: set_diff_map_func(imol)))


     coot_gui.add_simple_coot_menu_menuitem(
       submenu_maps,
       "Another (contour) level...",
       lambda func: coot.another_level())


     coot_gui.add_simple_coot_menu_menuitem(
       submenu_maps,
       "Multi-chicken...",
       lambda func: coot_gui.map_molecule_chooser_gui("Choose a molecule for multiple contouring",
		lambda imol: (coot.set_map_displayed(imol, 0), coot_utils.multi_chicken(imol))))

     
     #---------------------------------------------------------------------
     #     Molecule functions/Modelling
     #
     #---------------------------------------------------------------------

     def add_hydrogens_with_coot_reduce():
       with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                  aa_ins_code, aa_atom_name, aa_alt_conf]:
         coot.coot_reduce(aa_imol)

     def add_hydrogens_refmac_func():
       with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                  aa_ins_code, aa_atom_name, aa_alt_conf]:
         coot_utils.add_hydrogens_using_refmac(aa_imol)
     
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Add Hydrogen Atoms",
       lambda func: add_hydrogens_with_coot_reduce())


     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Add Hydrogens using Refmac",
       lambda func: add_hydrogens_refmac_func())


     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Add Other Solvent Molecules...",
       lambda func: coot_gui.solvent_ligands_gui())

     
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Arrange Waters Around Protein...",
       lambda func: coot_gui.molecule_chooser_gui(
          "Arrange waters in molecule: ",
          lambda imol: coot.move_waters_to_around_protein(imol)))


     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Assign (force) HETATMs for this Residue",
       lambda func: coot_utils.using_active_atom(hetify_residue,
                                      "aa_imol", "aa_chain_id", "aa_res_no", "aa_ins_code"))

     
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Assign HETATM to molecule...", 
       lambda func: coot_gui.molecule_chooser_gui("Assign HETATMs as per PDB definition", 
		lambda imol: coot.assign_hetatms(imol)))

     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Atoms with Zero Occupancies...",
       lambda func: coot_gui.molecule_chooser_gui(
         "Which molecule to check for Atoms with zero occupancies?",
         lambda imol: coot_gui.zero_occ_atoms_gui(imol)))

     # in main menu now
     # coot_gui.add_simple_coot_menu_menuitem(
     #   submenu_models,
     #   "Copy Coordinates Molecule...", 
     #   lambda func: coot_gui.molecule_chooser_gui("Molecule to Copy...", 
     #    	lambda imol: coot.copy_molecule(imol)))


     # moved to main menu now
     # should stay open if helper function returns False
     # def atom_selection_from_fragmemt_func(imol, text, button_state):
     #   print "BL DEBUG:: imol, text, button_state", imol, text, button_state
     #   jmol = coot.new_molecule_by_atom_selection(imol, text)
     #   if button_state:
     #     coot_utils.move_molecule_to_screen_centre(jmol)
     #   return coot_utils.valid_model_molecule_qm(jmol)
     # coot_gui.add_simple_coot_menu_menuitem(
     #   submenu_models,
     #   "Copy Fragment...", 
     #   lambda func: coot_gui.generic_chooser_and_entry_and_check_button("Create a new Molecule\n \
     #                              From which molecule shall we copy the fragment?", 
     #                                                           "Atom selection for fragment", "//A/1-10", "Move molecule here?", 
     #    	                                               lambda imol, text, button_state: atom_selection_from_fragmemt_func(imol, text, button_state),
     #                                          False))

     # --- D ---

     coot_gui.add_simple_coot_menu_menuitem(submenu_models, "Delete Hydrogen Atoms",
       lambda func: coot_utils.using_active_atom(coot.delete_hydrogens, "aa_imol"))

     coot_gui.add_simple_coot_menu_menuitem(submenu_models, "Delete Side-chains for Active Chain",
                                            lambda func: coot_utils.using_active_atom(delete_sidechains_for_chain, "aa_imol", "aa_chain_id"))

     # now in main menu
##     coot_gui.add_simple_coot_menu_menuitem(
##       submenu_models,
##       "DB Loop...",
##       lambda func: coot_gui.click_protein_db_loop_gui())

     
     # errr... move this...(??) . OK, but comment it out for now.
     # submenu = Gtk.Menu()
     # menuitem2 = Gtk.MenuItem("Dock Sequence...")
     # menuitem2.set_submenu(submenu)
     # menu.append(menuitem2)
     # menuitem2.show()
     
     # coot_gui.add_simple_coot_menu_menuitem(
     #  submenu,
     #  "Dock Sequence...", 
     #  lambda func: coot_gui.cootaneer_gui_bl())


     # def associate_seq_func(imol, chain_id, pir_file):
     #   import os, re
     #   chain_count = 0
     #   reg_chain = re.compile("chain", re.IGNORECASE)
     #   print("assoc seq:", imol, chain_id, pir_file)
     #   if (os.path.isfile(pir_file)):
     #     fin = open(pir_file, 'r')
     #     seq_text = fin.read()
     #     fin.close()
     #     coot.assign_pir_sequence(imol, chain_id, seq_text)
     #   else:
     #     print("BL WARNING:: could not find", pir_file)

     # old
     #add_simple_coot_menu_menuitem(
     #  submenu,
     #  "Associate Sequence...",
     #  lambda func: coot_gui.generic_chooser_entry_and_file_selector(
     #    "Associate Sequence to Model: ",
     #    coot_utils.valid_model_molecule_qm,
     #    "Chain ID",
     #    "",
     #    "Select PIR file",
     #    lambda imol, chain_id, seq_file_name: associate_seq_func(imol, chain_id, seq_file_name)))

     # FIXME
     # need to do this differently. Kevin's assign_sequence_from_file may
     # not be ideal (seems to only read the chain from format and wont
     # assign more than one chain...
     #
     #add_simple_coot_menu_menuitem(
     #  submenu,
     #  "Associate Sequence to Molecule...",
     #  lambda func:    coot_gui.generic_chooser_and_file_selector(
     #                  "Associate Sequence with Molecule: ",
     #                  coot_utils.valid_model_molecule_qm,
     #                  "Select Sequence File",
     #                  "",
     #                  lambda imol, sequence_file_name:
     #                  coot.assign_sequence_from_file(imol,
     #                                            sequence_file_name)
     #  ))
     
     # coot_gui.add_simple_coot_menu_menuitem(
     #   submenu,
     #   "Associate Sequence to Chain...",
     #   lambda func: coot_gui.associate_sequence_with_chain_gui()) # no alignment on OK press


     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Duplicate range (pick atoms)",
       lambda func: coot_gui.duplicate_range_by_atom_pick())
     
     # ---- F ---------

     # doublication to entry in main gtk code!
     # submenu = gtk.Menu()
     # menuitem2 = gtk.MenuItem("Find Secondary Structure...")

     # menuitem2.set_submenu(submenu)
     # submenu_models.append(menuitem2)
     # menuitem2.show()
     
     # coot_gui.add_simple_coot_menu_menuitem(
     #   submenu,
     #  "Find Helices", 
     #   lambda func: coot.find_helices())

     # coot_gui.add_simple_coot_menu_menuitem(
     #  submenu,
     #  "Find Strands", 
     #  lambda func: coot.find_strands())

     def get_smiles_pdbe_func():
       with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code, aa_atom_name, aa_alt_conf]:
         comp_id = coot.residue_name(aa_imol, aa_chain_id,
                                aa_res_no, aa_ins_code)
         # print "BL INFO:: here with residue name", comp_id
         coot_utils.get_SMILES_for_comp_id_from_pdbe(comp_id)
                         
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Fetch PDBe description for this ligand",
       lambda func: get_smiles_pdbe_func())


     def get_pdbe_ligand_func(comp_id):
         status = coot_utils.get_SMILES_for_comp_id_from_pdbe(comp_id)
         coot.get_monomer(comp_id)

     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Fetch PDBe Ligand Description",
       lambda func: coot_gui.generic_single_entry("Fetch PDBe Ligand Desciption for comp_id:",
                                         "", " Fetch ", lambda comp_id: get_pdbe_ligand_func(comp_id)))


     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Fix Nomenclature Errors...",
       lambda func: coot_gui.molecule_chooser_gui("Fix Nomenclature Error in molecule:",
                                         lambda imol: coot.fix_nomenclature_errors(imol)))


     # --- I --------

     def chiral_centre_inverter_func():
         with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                    aa_ins_code, aa_atom_name, aa_alt_conf]:
             coot.invert_chiral_centre(aa_imol, aa_chain_id, aa_res_no,
                                  aa_ins_code, aa_atom_name)
             
     coot_gui.add_simple_coot_menu_menuitem(
         submenu_models,
         "Invert This Chiral Centre",
         lambda func:
            chiral_centre_inverter_func())
     

     # --- J --------

     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "JLigand launch",
       lambda func:
          jligand_gui.launch_jligand_function())
     

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
         coot.make_link(imol_1, spec_1, spec_2, "dummy", 0.1)
       
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Make Link (click 2 atoms)...",
       lambda func:
       user_defined_click(2, make_link_ext_func))

       
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Merge Water Chains...",
       lambda func: coot_gui.molecule_chooser_gui("Merge Water Chains in molecule:",
                                         lambda imol: coot_utils.merge_solvent_chains(imol)))


     def mon_dict_func(text):
       idealized = 0
       new_model = coot.get_monomer_from_dictionary(text, idealized)
       if not coot_utils.valid_model_molecule_qm(new_model):
         coot.get_monomer(text)
         
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Monomer from Dictionary",
       lambda func:
         coot_gui.generic_single_entry("Pull coordinates from CIF dictionary for 3-letter-code:", "",
                              " Get Coords ",
                              lambda text: mon_dict_func(text)))
                              

     def morph_fit_chain_func(radius):
       with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                  aa_ins_code, aa_atom_name, aa_alt_conf]:
         coot.morph_fit_chain(aa_imol, aa_chain_id, radius)
         
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Morph Fit Chain (Averaging Radius 7)",
       lambda func: morph_fit_chain_func(7)
       )
     
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Morph Fit Chain (Averaging Radius 11)",
       lambda func: morph_fit_chain_func(11)
       )


     def morph_fit_ss_func(radius):
       with coot_utils.UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                  aa_ins_code, aa_atom_name, aa_alt_conf]:
         coot.morph_fit_by_secondary_structure_elements(aa_imol, aa_chain_id)
         
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Morph Fit Chain based on Secondary Structure",
       lambda func: morph_fit_chain_func(7)
       )
     

     # -- N --

     def new_mol_sphere_func1(imol, text):
       try:
         radius = float(text)
       except:
         print("WARNING:: no valid radius", text)
       args = [imol] + coot_utils.rotation_centre() + [radius, 0]
       coot.new_molecule_by_sphere_selection(*args)

     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
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
         coot.new_molecule_by_symop(imol, text,
                               pre_shift[0],
                               pre_shift[1],
                               pre_shift[2])
       
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "New Molecule from Symmetry Op...",
       lambda func: coot_gui.generic_chooser_and_entry(
          "Molecule from which to generate a symmetry copy",
          "SymOp", "X,Y,Z",
          lambda imol, text: new_mol_sym_func1(imol, text)))


     # --- P ---
     
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Phosphorylate this residue",
       lambda func: coot_utils.phosphorylate_active_residue())

     # ---- R ---------

     # --- Ren --- 

     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Rename Residue...",
       lambda func: coot_gui.rename_residue_gui())


     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Renumber Waters...",
       lambda func: coot_gui.molecule_chooser_gui(
          "Renumber waters of which molecule?",
          lambda imol: coot.renumber_waters(imol)))

     # --- Reo ---

     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Reorder Chains...",
       lambda func: coot_gui.molecule_chooser_gui("Sort Chain IDs in molecule:",
                                         lambda imol: coot.sort_chains(imol))) # an internal function

     # --- Rep ---

     # BL says:: may work, not sure about function entirely
     # this is in main menu now
     # coot_gui.add_simple_coot_menu_menuitem(
     #   submenu_models,
     #   "Replace Fragment...",
     #   lambda func: coot_gui.molecule_chooser_gui("Define the molecule that needs updating",
     #    	lambda imol_base: coot_gui.generic_chooser_and_entry(
     #    			"Molecule that contains the new fragment:",
     #    			"Atom Selection","//",
     #    			lambda imol_fragment, atom_selection_str:
     #    			replace_fragment(imol_base, imol_fragment, atom_selection_str))))

     # in main menu now
     # coot_gui.add_simple_coot_menu_menuitem(
     #   submenu_models,
     #   "Replace Residue...",
     #   lambda func: coot_gui.generic_single_entry("Replace this residue with residue of type:",
     #                                     "ALA", "Mutate",
     #                                     lambda text: coot_utils.using_active_atom(mutate_by_overlap,
     #                                                                    "aa_imol", "aa_chain_id", "aa_res_no",
     #                                                                    text)))

     # --- Res ---
     
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Residue Type Selection...",
	lambda func: coot_gui.generic_chooser_and_entry("Choose a molecule to select residues from: ","Residue Type:","",
                                               lambda imol, text: (coot.new_molecule_by_residue_type_selection(imol, text),
                                                                   coot.update_go_to_atom_window_on_new_mol())))


     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Residues with Alt Confs...",
       lambda func: coot_gui.molecule_chooser_gui(
         "Which molecule to check for Alt Confs?",
         lambda imol: coot_gui.alt_confs_gui(imol)))


     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Residues with Cis Peptide Bonds...",
       lambda func: coot_gui.molecule_chooser_gui("Choose a molecule for checking for Cis Peptides",
                                         lambda imol: coot_gui.cis_peptides_gui(imol)))


     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Residues with Missing Atoms...",
       lambda func: coot_gui.molecule_chooser_gui(
         "Which molecule to check for Missing Atoms?",
         lambda imol: coot_gui.missing_atoms_gui(imol)))

     # --- Rig ---

     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Rigid Body Fit Residue Ranges...",
       lambda func:
       coot_gui.residue_range_gui(lambda imol, ls: rigid_body_refine_by_residue_ranges(imol, ls),
                         "Rigid Body Refine",
                         "  Fit  "))

     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Rigid Body Fit Molecule...",
       lambda func: coot_gui.molecule_chooser_gui("Rigid Body Fit Molecule",
                lambda imol: coot.rigid_body_refine_by_atom_selection(imol, "//")))                                         
       
     # ---- S --------

     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Superpose ligands",
       lambda func: coot_gui.superpose_ligand_gui())
     

     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "Symm Shift Reference Chain Here",
       lambda func: coot.move_reference_chain_to_symm_chain_position())

     
     # ---- U ---------
     
     # submenu = gtk.Menu()
     # menuitem2 = gtk.MenuItem("Rotamer Search...")
     
     # menuitem2.set_submenu(submenu)
     # submenu_models.append(menuitem2)
     # menuitem2.show()

     # coot_gui.add_simple_coot_menu_menuitem(
     #  submenu,
     #  "Use \"Backrub\" Rotamers",
     # lambda func: coot.set_rotamer_search_mode(ROTAMERSEARCHLOWRES)
     #  )
     
     #coot_gui.add_simple_coot_menu_menuitem(
     #  submenu,
     #  "DONT use \"Backrub\" Rotamers",
     #  lambda func: coot.set_rotamer_search_mode(ROTAMERSEARCHHIGHRES)
     #  )


     # coot_gui.add_simple_coot_menu_menuitem(
     #  submenu_models,
     #  "Use SEGIDs...",
     #  lambda func: coot_gui.molecule_chooser_gui("Exchange the Chain IDs, replace with SEG IDs",
     #  lambda imol: coot.exchange_chain_ids_for_seg_ids(imol)))

     
     # ---- W ---------

     def whats_this():
       central_residue = active_residue()
       res_name = coot.residue_name(*central_residue[0:4])
       mol_no = central_residue[0]
       n = comp_id2name(res_name)
       s = "(mol. no: " + str(mol_no) + ")  " + \
           res_name  + ":  " + \
           n if isinstance(n, str) else " <no-name-found>"
       coot.add_status_bar_text(s)
       
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_models,
       "What's this?",
       lambda func: whats_this()
       )
     

     #---------------------------------------------------------------------
     #     NCS functions
     #
     #---------------------------------------------------------------------

     def copy_ncs_range_func(imol, chain_id, text1, text2):
       try:
         r1 = int(text1)
         r2 = int(text2)
         if (chain_id == ncs.ncs_master_chain_id(imol)):
           # hunkey dorey
           coot.copy_residue_range_from_ncs_master_to_others(imol, chain_id, r1, r2)
         else:
           # different given master to current master
           # ask what to do.
           txt = "Current master chain is %s, but you asked to copy from %s.\n" \
                 %(ncs.ncs_master_chain_id(imol), chain_id)
           txt += "Apply this master change?\n\n"
           txt += "N.B. if no, then nothing is copied."
           r = coot_gui.yes_no_dialog(txt, "Change Master")
           if r:
             coot.ncs_control_change_ncs_master_to_chain_id(imol, chain_id)
             coot.copy_residue_range_from_ncs_master_to_others(imol, chain_id, r1, r2)
             ## could change master back?!
           else:
             coot.info_dialog("Master chain was not changed and copy not applied.")
       except:
         print("BL WARNING:: no valid number input")

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
           coot.copy_from_ncs_master_to_others(imol, chain_id)
         else:
           # different given master to current master
           # ask what to do.
           txt = "Current master chain is %s, but you asked to copy from %s.\n" \
                 %(ncs.ncs_master_chain_id(imol), chain_id)
           txt += "Apply this master change?\n\n"
           txt += "N.B. if no, then nothing is copied."
           r = coot_gui.yes_no_dialog(txt, "Change Master")
           if r:
             coot.ncs_control_change_ncs_master_to_chain_id(imol, chain_id)
             coot.copy_from_ncs_master_to_others(imol, chain_id)
             ## could change master back?!
           else:
             coot.info_dialog("Master chain was not changed and copy not applied.")
       else:
         s = "You need to define NCS operators for molecule " + str(imol)
         coot.info_dialog(s)
           
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_ncs,
       "Copy NCS Chain...",
       lambda func: coot_gui.generic_chooser_and_entry("Apply NCS edits from NCS Master Chain to Other Chains",
                                              "Master Chain ID",
                                              coot_utils.get_first_ncs_master_chain(),  # can return  "".
                                              lambda imol, chain_id: copy_ncs_chain_func(imol, chain_id)))


     def ncs_ghost_res_range_func(imol):
       from types import ListType
       label_1 = "Start ResNumber"
       label_2 = "End ResNumber"
       entry_1_default_text = "10"
       entry_2_default_text = "20"
       go_button_label = "Regenerate Ghosts"
       def handle_go_function(resno_1_text, resno_2_text):
         cont = False
         try:
           resno_1 = int(resno_1_text)
           resno_2 = int(resno_2_text)
           cont = True
         except:
           print("BL WARNING:: input residue numbers have to be integers")
         if (cont):
           ghost_ncs_chain_ids = ncs_chain_ids(imol)
           if (type(ghost_ncs_chain_ids) is ListType):
             # because we can have hetero-NCS,
             # but we ignore NCS other that
             # that of the first type.
             ghost_chain_list = ghost_ncs_chain_ids[0]
             ncs.manual_ncs_ghosts(imol, resno_1, resno_2, ghost_chain_list)
         
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
       lambda func: coot.set_ncs_matrix_type(0))


     coot_gui.add_simple_coot_menu_menuitem(
       submenu,
       "Fast (LSQ)",
       lambda func: coot.set_ncs_matrix_type(1))


     coot_gui.add_simple_coot_menu_menuitem(
       submenu,
       "Extra Fast (LSQ, only every 2nd CA)",
       lambda func: coot.set_ncs_matrix_type(2))


     # ---------------------------------------------------------------------
     #     Refinement
     # ---------------------------------------------------------------------
     #

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
        coot.add_planar_peptide_restraints()

     coot_gui.add_simple_coot_menu_menuitem(
       submenu,
       "Add Planar Peptide Restraints",
       lambda func: add_restr_func1())


     def add_restr_func2():
        print('Planar Peptide Restraints removed')
        coot.remove_planar_peptide_restraints()

     coot_gui.add_simple_coot_menu_menuitem(
       submenu, "Remove Planar Peptide Restraints",
       lambda func: add_restr_func2())


     def shelx_ref_func():
       window = Gtk.Window(gtk.WINDOW_TOPLEVEL)
       vbox = Gtk.VBox(False, 0)
       hbox = Gtk.HBox(False, 0)
       go_button = Gtk.Button("  Refine  ")
       cancel_button = Gtk.Button("  Cancel  ")
       entry_hint_text = "HKL data filename \n(leave blank for default)"
       chooser_hint_text = " Choose molecule for SHELX refinement  "
       h_sep = Gtk.HSeparator()

       window.add(vbox)
       option_menu_mol_list_pair = coot_gui.generic_molecule_chooser(vbox, chooser_hint_text)
       entry = coot_gui.file_selector_entry(vbox, entry_hint_text)

       def shelx_delete_event(*args):
         window.destroy()
         return False

       def shelx_go_funcn_event(*args):
         import operator
         txt = entry.get_text()
         imol = coot_gui.get_option_menu_active_molecule(*option_menu_mol_list_pair)
         if (isinstance(imol, numbers.Number)):
           if (len(txt) == 0):
             shelx.shelxl_refine(imol)
           else:
             shelx.shelxl_refine(imol, txt)
         window.destroy()
         return False

       go_button.connect("clicked", shelx_go_funcn_event)
       cancel_button.connect("clicked", shelx_delete_event)

       vbox.pack_start(h_sep, False, False, 2)
       vbox.pack_start(hbox, False, False, 2)
       hbox.pack_start(go_button, True, False, 0)
       hbox.pack_start(cancel_button, True, False, 0)
       window.show_all()
       
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_refine,
       "SHELXL Refine...", 
       lambda func: shelx_ref_func())


     coot_gui.add_simple_coot_menu_menuitem(
       submenu_refine,
       "Read REFMAC logfile...",
       lambda func: coot_gui.generic_chooser_and_file_selector("Read Refmac log file",
                       coot_utils.valid_model_molecule_qm, "Logfile name: ", "",
                       lambda imol, text: refmac.read_refmac_log(imol, text)))
                       
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_refine,
       "Occupancy refinement input for REFMAC...",
       lambda func: coot_gui.generic_chooser_and_file_selector("Extra restraints file",
                       coot_utils.valid_model_molecule_qm, "Restraints file name: ",
                       "refmac_extra_params.txt",
                       lambda imol, text: refmac.restraints_for_occupancy_refinement(imol, text)))
       

     coot_gui.add_simple_coot_menu_menuitem(
       submenu_refine,
       "Auto-weight refinement",
       lambda func: coot_utils.auto_weight_for_refinement())
     

     coot_gui.add_simple_coot_menu_menuitem(
       submenu_refine,
       "Set Undo Molecule...",
       lambda func: coot_gui.molecule_chooser_gui("Set the Molecule for 'Undo' Operations",
                    lambda imol: coot.set_undo_molecule(imol)))


     # BL says: has no checking for text = number yet
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_refine,
       "B factor bonds scale factor...",
       lambda func: coot_gui.generic_chooser_and_entry("Choose a molecule to which the B-factor colour scale is applied:",
                "B-factor scale:", "1.0", 
                lambda imol, text: coot.set_b_factor_bonds_scale_factor(imol,float(text))))


     def set_mat_func(text):
        import operator
        t = float(text)
        if isinstance(t, numbers.Number):
                s = "Matrix set to " + text
                coot.set_matrix(t)
                coot.add_status_bar_text(s)
        else:
                coot.add_status_bar_text("Failed to read a number")

     coot_gui.add_simple_coot_menu_menuitem(
       submenu_refine,
       "Set Matrix (Refinement Weight)...",
       lambda func: coot_gui.generic_single_entry("set matrix: (smaller means better geometry)", 
                str(coot.matrix_state()), "Set it", 
                lambda text: set_mat_func(text)))


     def set_den_gra_func(text):
        import operator
        t = float(text)
        if isinstance(t, numbers.Number):
                s = "Density Fit scale factor set to " + text
                coot.set_residue_density_fit_scale_factor(t)
                coot.add_status_bar_text(s)
        else:
                coot.add_status_bar_text("Failed to read a number")



#      # ---------------------------------------------------------------------
#      #     Recent structures from the PDBe
#      # ---------------------------------------------------------------------
#      #
#      coot_gui.add_simple_coot_menu_menuitem(
#        submenu_pdbe, "PDBe recent structures...",
#        lambda func: get_recent_pdbe.pdbe_latest_releases_gui())

#      # we do test for refmac at startup not runtime (for simplicity)
#      if coot_utils.command_in_path_qm("refmac5"):
#        mess = " Get it "
#      else:
#        mess = "\n  WARNING::refmac5 not in the path - SF calculation will fail  \n\n"
       
#      coot_gui.add_simple_coot_menu_menuitem(
#        submenu_pdbe, "Get from PDBe...",
#        lambda func: coot_gui.generic_single_entry("Get PDBe accession code",
#                                          "", " Get it ",
#                                          lambda text:
#                                          get_recent_pdbe.pdbe_get_pdb_and_sfs_cif("include-sfs", text.rstrip().lstrip())))

#      # ---------------------------------------------------------------------
#      #     Tutorial data
#      # ---------------------------------------------------------------------
#      #
#      def load_tutorial_data_func():
#        data_dir = False
#        prefix_dir = os.getenv("COOT_PREFIX")
#        if not prefix_dir:
#          pkg_data_dir = pkgdatadir()
#        else:
#          pkg_data_dir = os.path.join(prefix_dir, "share", "coot")
#        if os.path.isdir(pkg_data_dir):  
#          data_dir = os.path.join(pkg_data_dir, "data")
#        if data_dir:
#          pdb_file_name = os.path.join(data_dir, "tutorial-modern.pdb")
#          mtz_file_name = os.path.join(data_dir, "rnasa-1.8-all_refmac1.mtz")

#          if os.path.isfile(pdb_file_name):
#            coot.read_pdb(pdb_file_name)
#          if os.path.isfile(mtz_file_name):
#            coot.make_and_draw_map(mtz_file_name, "FWT", "PHWT", "", 0, 0)
#            coot.make_and_draw_map(mtz_file_name, "DELFWT", "PHDELWT", "", 0, 1)
       
#      coot_gui.add_simple_coot_menu_menuitem(
#        menu,
#        "Load tutorial model and data",
#        lambda func: load_tutorial_data_func()
#        )

     
     # ---------------------------------------------------------------------
     #     Views/Representations
     # ---------------------------------------------------------------------
     #

     def make_dot_surf_func(imol,text):
        # I think a single colour is better than colour by atom
        coot.set_dots_colour(imol, 0.5, 0.5, 0.5)
        density = 1.0
        dots_handle = coot.dots(imol, text, text, density, 1)
        print("INFO::dots handle: ", dots_handle)

     coot_gui.add_simple_coot_menu_menuitem(
       submenu_representation,
       "Dotted Surface...",
       lambda func: coot_gui.generic_chooser_and_entry("Surface for molecule", 
                "Atom Selection:", "//A/1-2", 
                lambda imol, text: make_dot_surf_func(imol, text)))


     def clear_dot_surf_func(imol,text):
         try:
             n = int(text)
             coot.clear_dots(imol,n)
         except:
             print("WARNING:: dots handle number should be an integer")

     coot_gui.add_simple_coot_menu_menuitem(
       submenu_representation,
       "Clear Surface Dots...",
       lambda func: coot_gui.generic_chooser_and_entry("Molecule with Dotted Surface", 
                "Dots Handle Number:", "0", 
                lambda imol, text: clear_dot_surf_func(imol, text)))

     import coot_hole
     coot_gui.add_simple_coot_menu_menuitem(
         submenu_representation,
         "HOLE...",
         lambda func: coot_hole.hole_ify())

     # Views submenu
     submenu = Gtk.Menu()
     menuitem2 = Gtk.MenuItem("Views")
 
     menuitem2.set_submenu(submenu)
     draw_menu.append(menuitem2)
     menuitem2.show()

     coot_gui.add_simple_coot_menu_menuitem(
       submenu,
       "Add View...",
       lambda func: coot_gui.view_saver_gui())

     # BL says:: maybe check if number at some point
     coot_gui.add_simple_coot_menu_menuitem(
       submenu,
       "Add a Spin View...",
       lambda func: coot_gui.generic_double_entry("Number of Steps", 
                                                  "Number of Degrees (total)", "3600", "360", 
                                                  False, False,                 #check button text and callback
                                                  "  Add Spin-xxx  ",
                                                  lambda text_1, text_2, active_state: coot.add_spin_view("Spin", int(text_1), float(text_2))))

     coot_gui.add_simple_coot_menu_menuitem(
       submenu,
       "Views Panel...",
       lambda func: coot_gui.views_panel_gui())
 
     coot_gui.add_simple_coot_menu_menuitem(
       submenu,
       "Play Views",
       lambda func: list(map(eval,["go_to_first_view(1)",
                              "time.sleep(1)", "play_views()"])))
 
     # BL says:: maybe check if number at some point
     coot_gui.add_simple_coot_menu_menuitem(
       submenu, "Set Views Play Speed...",
       lambda func: coot_gui.generic_single_entry("Set Views Play Speed",
                        str(coot.views_play_speed()), "  Set it  ",
                        lambda text: coot.set_views_play_speed(float(text))))


     coot_gui.add_simple_coot_menu_menuitem(
       submenu, "Save Views...",
       lambda func: coot_gui.generic_single_entry("Save Views",
                                         "coot-views.py", " Save ",
                                         lambda txt: coot.save_views(txt)))

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


     # ---------------------------------------------------------------------
     #     Modules
     # ---------------------------------------------------------------------

     import refine

     coot_gui.add_simple_coot_menu_menuitem(submenu_modules, "CCP4", lambda func: coot_gui.add_module_ccp4())

     coot_gui.add_simple_coot_menu_menuitem(submenu_modules, "Carbohydrate", lambda func: gui_add_linked_cho.add_module_carbohydrate_gui())

     coot_gui.add_simple_coot_menu_menuitem(submenu_modules, "Cryo-EM", lambda func: coot_gui.add_module_cryo_em())

     coot_gui.add_simple_coot_menu_menuitem(submenu_modules, "ProSMART", lambda func: gui_prosmart.add_module_prosmart())

     coot_gui.add_simple_coot_menu_menuitem(submenu_modules, "SHELX", lambda func: shelx_extensions.add_module_shelx())

     coot_gui.add_simple_coot_menu_menuitem(submenu_modules, "Refine", lambda func: refine.add_module_refine())

     coot_gui.add_simple_coot_menu_menuitem(submenu_modules, "Restraints", lambda func: add_module_user_defined_restraints())

     # ---------------------------------------------------------------------
     #     Settings
     # ---------------------------------------------------------------------

     submenu = Gtk.Menu()
     menuitem2 = Gtk.MenuItem("Rotate Translate Zone Mode...")
 
     menuitem2.set_submenu(submenu)
     submenu_settings.append(menuitem2)
     menuitem2.show()

     coot_gui.add_simple_coot_menu_menuitem(
       submenu, "Rotate About Fragment Centre",
       lambda func: coot.set_rotate_translate_zone_rotates_about_zone_centre(1))


     coot_gui.add_simple_coot_menu_menuitem(
       submenu, "Rotate About Second Clicked Atom",
       lambda func: coot.set_rotate_translate_zone_rotates_about_zone_centre(0))


     coot_gui.add_simple_coot_menu_menuitem(
       submenu_settings,
       "Set Density Fit Graph Weight...",
       lambda func: coot_gui.generic_single_entry("set weight (smaller means apparently better fit)",
                str("%.2f" %residue_density_fit_scale_factor()), "Set it",
                lambda text: set_den_gra_func(text)))

     # BL says:: maybe check if number at some point
     coot_gui.add_simple_coot_menu_menuitem(
       submenu_settings,
       "Set Spin Speed",
       lambda func: coot_gui.generic_single_entry("Set Spin Speed (smaller is slower)",
                        str(coot.idle_function_rotate_angle()), "Set it",
                        lambda text: coot.set_idle_function_rotate_angle(float(text))))


     coot_gui.add_simple_coot_menu_menuitem(
       submenu_settings, "Nudge Centre...",
       lambda func: coot_gui.nudge_screen_centre_gui())


     def all_mol_symm_func():
        for imol in coot_utils.molecule_number_list():
                if coot_utils.valid_model_molecule_qm(imol):
                        coot.set_symmetry_whole_chain(imol, 1)

     coot_gui.add_simple_coot_menu_menuitem(
       submenu_settings, "All Molecules use \"Near Chains\" Symmetry", 
       lambda func: [valid_model_molecule_qm(imol) and
                         coot.set_symmetry_whole_chain(imol, 1) for imol in coot_utils.molecule_number_list()])


     coot_gui.add_simple_coot_menu_menuitem(
       submenu_refine, "Question Accept Refinement", 
       lambda func: coot.set_refinement_immediate_replacement(0))


     coot_gui.add_simple_coot_menu_menuitem(
       submenu_settings, "Save Graphics Size and Positions",
       lambda func: coot.graphics_window_size_and_position_to_preferences())


     def save_dialog_func():
       coot.post_model_fit_refine_dialog()
       coot.post_go_to_atom_window()

       def delete_event(*args):
         window.destroy()
         return False
       
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

       def go_func(*args):
         save_dialog_positions_to_preferences_file()
         window.destroy()
         
       go_button.connect("clicked", go_func)
       cancel_button.connect("clicked", lambda w: window.destroy())

       window.show_all()

     coot_gui.add_simple_coot_menu_menuitem(
       submenu_settings, "Save Dialog Positions...",
       lambda func: save_dialog_func())


     coot_gui.add_simple_coot_menu_menuitem(
       submenu_settings, "Key Bindings...",
       lambda func: coot_gui.key_bindings_gui())

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
       lambda func:
          coot_gui.generic_single_entry("Checkpoint interval (seconds)",
                               "30",
                               " Start Auto-saving ",
                               lambda txt:
                                  quick_save_func(txt)))

     # Doesnt seem to be working right currently, so comment out?! Not any more?!
     # add to validate menu
     menu = coot_gui.coot_menubar_menu("Validate")

     coot_gui.add_simple_coot_menu_menuitem(menu, "Pukka Puckers...?",
                                   lambda func: coot_gui.molecule_chooser_gui(
       "Choose a molecule for ribose pucker analysis",
       lambda imol: coot_utils.pukka_puckers_qm(imol)))

     coot_gui.add_simple_coot_menu_menuitem(
         menu,
         "Alignment vs PIR...",
         #lambda func: coot_gui.molecule_chooser_gui("Alignment vs PIR info for molecule:",
         #   lambda imol: coot_gui.wrapper_alignment_mismatches_gui(imol)))
         lambda func: coot_gui.associate_pir_wih_molecule_gui(True))

  else:
    print("BL WARNING:: could not find the main_menubar! Sorry, no extensions menu!")

