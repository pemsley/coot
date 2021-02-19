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

import pygtk, gtk, pango
import time


def add_coot_menu_separator(menu):
  sep = gtk.MenuItem()
  menu.add(sep)
  #   sep.props.sensitive = False
  sep.show()
   
if (have_coot_python):
  if coot_python.main_menubar():

     # --------------------------------------------------
     #           coot news dialog and updates dialog
     # --------------------------------------------------
#   comment these out for now     
#      menu = coot_menubar_menu("About")
#      if (menu):
#        add_simple_coot_menu_menuitem(menu, "Coot News...",
#                                      lambda func: whats_new_dialog())

#        os_type = os.name
#        if not os_type == 'mac':
#          add_simple_coot_menu_menuitem(menu, "Check for Updates...",
#                                        lambda func: (printf("checking for updates..."),
#                                                      check_for_updates_gui()))


     # --------------------------------------------------
     #           coordinated water validation dialog
     # --------------------------------------------------

     menu = coot_menubar_menu("Validate")
     if menu:
       add_simple_coot_menu_menuitem(menu, "Highly coordinated waters...",
                                     lambda func: water_coordination_gui())

       def validation_outliers_func():
         with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                    aa_ins_code, aa_atom_name, aa_alt_conf]:
           imol_map = imol_refinement_map()
           if not valid_map_molecule_qm(imol_map):
             info_dialog_and_text("Refinement Map is currently not set")
           else:
             validation_outliers_dialog(aa_imol, imol_map)
       add_simple_coot_menu_menuitem(menu, "Validation Outliers",
                                     lambda func: validation_outliers_func())

       add_simple_coot_menu_menuitem(menu, "List Ramachandran outliers...",
                                     lambda func: rama_outlier_gui())


     # --------------------------------------------------
     #           user_define_restraints plugin
     # --------------------------------------------------

     def add_module_user_defined_restraints():
       menu = coot_menubar_menu("Restraints")
       load_from_search_load_path("user_define_restraints.py")
     

     
     # ---------------------------------------------
     #           extensions
     # ---------------------------------------------

     menu = coot_menubar_menu("E_xtensions")

     # make submenus:
     submenu_all_molecule = gtk.Menu()
     menuitem_2 = gtk.MenuItem("All Molecule...")
     submenu_maps = gtk.Menu()
     menuitem_3 = gtk.MenuItem("Maps...")
     submenu_models = gtk.Menu()
     menuitem_4 = gtk.MenuItem("Modelling...")
     submenu_refine = gtk.Menu()
     menuitem_5 = gtk.MenuItem("Refine...")
     submenu_representation = gtk.Menu()
     menuitem_6 = gtk.MenuItem("Representations")
     submenu_settings = gtk.Menu()
     menuitem_7 = gtk.MenuItem("Settings...")
     submenu_pisa = gtk.Menu()
     menuitem_pisa = gtk.MenuItem("PISA...")
     submenu_pdbe = gtk.Menu()
     menuitem_pdbe = gtk.MenuItem("PDBe...")
     submenu_modules = gtk.Menu()
     menuitem_modules = gtk.MenuItem("Modules...")
     submenu_ncs = gtk.Menu()
     menuitem_ncs = gtk.MenuItem("NCS...")

     menuitem_2.set_submenu(submenu_all_molecule)
     menu.append(menuitem_2)
     menuitem_2.show()
     
     menuitem_3.set_submenu(submenu_maps)
     menu.append(menuitem_3)
     menuitem_3.show()
     
     menuitem_4.set_submenu(submenu_models)
     menu.append(menuitem_4)
     menuitem_4.show()
     
     menuitem_ncs.set_submenu(submenu_ncs)
     menu.append(menuitem_ncs)
     menuitem_ncs.show()
     
     menuitem_5.set_submenu(submenu_refine)
     menu.append(menuitem_5)
     menuitem_5.show()
     
     menuitem_6.set_submenu(submenu_representation)
     menu.append(menuitem_6)
     menuitem_6.show()
     
     menuitem_pisa.set_submenu(submenu_pisa)
     menu.append(menuitem_pisa)
     menuitem_pisa.show()

     menuitem_7.set_submenu(submenu_settings)
     menu.append(menuitem_7)
     menuitem_7.show()

     menuitem_modules.set_submenu(submenu_modules)
     menu.append(menuitem_modules)
     menuitem_modules.show()

     menuitem_pdbe.set_submenu(submenu_pdbe)
     menu.append(menuitem_pdbe)
     menuitem_pdbe.show()
     
     

     #---------------------------------------------------------------------
     #     Post MR
     #
     #---------------------------------------------------------------------

     add_simple_coot_menu_menuitem(
       submenu_all_molecule,
       "[Post MR] Fill Partial Residues...",
       lambda func: molecule_chooser_gui("Find and Fill residues with missing atoms",
		lambda imol: fill_partial_residues(imol)))


     # old style - not interruptable
     #     add_simple_coot_menu_menuitem(
     #       submenu_all_molecule,
     #       "[Post MR] Fit Protein...",
     #       lambda func: molecule_chooser_gui("Fit Protein using Rotamer Search",
     #		lambda imol: (imol_refinement_map() == -1 and
     #                              add_status_bar_text("oops. Must set a map to fit") or
     #                              fit_protein(imol))))
     #
     #
     #     add_simple_coot_menu_menuitem(
     #       submenu_all_molecule,
     #       "[Post MR] Stepped Refine...",
     #       lambda func: molecule_chooser_gui("Stepped Refine: ",
     #		lambda imol: (imol_refinement_map() == -1 and
     #                              add_status_bar_text("oops. Must set a map to fit") or
     #                              stepped_refine_protein(imol))))
     #
     #     add_simple_coot_menu_menuitem(
     #       submenu_all_molecule,
     #       "Refine/Improve Ramachandran Plot...",
     #       lambda func: molecule_chooser_gui("Refine Protein with Ramachanran Plot Optimization: ",
     #                lambda imol: (imol_refinement_map() == -1 and
     #                              add_status_bar_text("oops. Must set a map to fit") or
     #                              stepped_refine_protein_for_rama(imol))))
     #                                         

     # BL says:: we cannot do this with lambda functions in python
     # statements are not allowed! Needed for globals!!
     def fit_protein_func1(imol):
       if (imol_refinement_map() == -1):
         add_status_bar_text("oops. Must set a map to fit")
       else:
         global continue_multi_refine
         continue_multi_refine = True
         interruptible_fit_protein(imol, fit_protein_fit_function)

     add_simple_coot_menu_menuitem(
       submenu_all_molecule,
       "Fit Protein...",
       lambda func: molecule_chooser_gui("Fit Protein using Rotamer Search",
		lambda imol: fit_protein_func1(imol)))


     def fit_protein_func2(imol):
       if (imol_refinement_map() == -1):
         add_status_bar_text("oops. Must set a map to fit")
       else:
         global continue_multi_refine
         continue_multi_refine = True
         interruptible_fit_protein(imol, fit_protein_stepped_refine_function)

     add_simple_coot_menu_menuitem(
       submenu_all_molecule,
       "Stepped Refine...",
       lambda func: molecule_chooser_gui("Fit Protein using Real-Space Refinement",
		lambda imol: fit_protein_func2(imol)))


     def fit_protein_func3(imol):
       if (imol_refinement_map() == -1):
         add_status_bar_text("oops. Must set a map to fit")
       else:
         global continue_multi_refine
         continue_multi_refine = True
         interruptible_fit_protein(imol, fit_protein_rama_fit_function)
         
     add_simple_coot_menu_menuitem(
       submenu_all_molecule,
       "Refine/Improve Ramachandran Plot...",
       lambda func: molecule_chooser_gui("Refine Protein with Ramachanran Plot Optimization: ",
                lambda imol: fit_protein_func3(imol)))
                                         

     #---------------------------------------------------------------------
     #     Map functions
     #
     #---------------------------------------------------------------------

     def mask_map_func():
	f = ""
	molecule_list = molecule_number_list()
	if not molecule_list == []:
          for i in molecule_list:
            if is_valid_map_molecule(molecule_list[i]):
              print "%s is a valid map molecule" %molecule_list[i]
              f = str(molecule_list[i])
              break
	else:
		print "BL WARNING:: dunno what to do!? No map found"
                f = False
	return f
      
     def mask_map_func1(active_state):
	print "changed active_state to ", active_state
        
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
         print "BL WARNING:: input %s for Map molecule number is not an integer.\nBailing out" %n

       if (continue_qm):
         if (invert_mask_qm):
           invert = 1
         else:
           invert = 0
         print "debug:: invert-mask? is", invert
         if (text_3 != "default"):
           try:
             new_radius = float(text_3)
             set_map_mask_atom_radius(new_radius)
           except:
             print "BL WARNING:: could not set map mask radius to %s. It's not a number" %text_3
         mask_map_by_atom_selection(n, imol, text_2, invert)
                
     def mask_map_radius_func():
       rad = map_mask_atom_radius()
       if (rad > 0):
         return str(rad)
       else:
         return "default"
       
     add_simple_coot_menu_menuitem(
       submenu_maps,
       "Mask Map by Atom Selection...",
       lambda func: molecule_chooser_gui("Define the molecule that has atoms to mask the map",
             lambda imol: generic_multiple_entries_with_check_button(
                         [[" Map molecule number: ", mask_map_func()],
                          [" Atom selection: ", "//A/1"],
                          ["Radius around atoms: ", mask_map_radius_func()]],
                         [" Invert Masking? ", lambda active_state: mask_map_func1(active_state)],
                         "  Mask Map  ", lambda text_list, invert: mask_map_func2(imol, text_list, invert))))


     add_simple_coot_menu_menuitem(
       submenu_maps,
       "Copy Map...",
       lambda func: map_molecule_chooser_gui("Map to Copy...", 
		lambda imol: copy_molecule(imol)))

     
     add_simple_coot_menu_menuitem(
       submenu_maps,
       "Make a Smoother Copy...", 
       lambda func: map_molecule_chooser_gui("Map Molecule to Smoothenize...", 
		lambda imol: smooth_map(imol, 1.25)))

     
     add_simple_coot_menu_menuitem(
       submenu_maps,
       "Make a Very Smooth Copy...", 
       lambda func: map_molecule_chooser_gui("Map Molecule to Smoothenize...", 
		lambda imol: smooth_map(imol, 2.0)))

     
     add_simple_coot_menu_menuitem(
       submenu_maps,
       "Make a Difference Map...",
       lambda func: make_difference_map_gui())


     add_simple_coot_menu_menuitem(
       submenu_maps,
       "Transform map by LSQ model fit...",
       lambda func: transform_map_using_lsq_matrix_gui())


     add_simple_coot_menu_menuitem(
       submenu_maps,
       "Average Maps...",
       lambda func: average_map_gui())

#      add_simple_coot_menu_menuitem(
#        submenu_maps,
#        "Export map...",
#        lambda func: generic_chooser_and_file_selector("Export Map: ",
#                 valid_map_molecule_qm, "File_name: ", "",
#                 lambda imol, text: (export_map(imol, text) == 1 and
#                                     add_status_bar_text("Map " + str(imol) + " exported to " + text))))
                                    

#      def export_map_func(imol, radius_string, file_name):
#        radius = -1.
#        try:
#          radius = float(radius_string)
#        except:
#          print "BL WARNING:: radius %s was no number!" %(radius_string)
#        if (radius >= 0):
#          export_map_fragment(*([imol] + rotation_centre() + \
#                                [radius, file_name]))
         
#      add_simple_coot_menu_menuitem(
#        submenu_maps,
#        "Export Local Map Fragment...", # (for Pymol, say)
#        lambda func: generic_chooser_entry_and_file_selector(
#                "Export Map: ",
#                valid_map_molecule_qm,
#                "Radius (A): ", "10",
#                "File-name: ",
#                lambda imol, radius_string, file_name:
#                    export_map_func(imol, radius_string, file_name)
#                )
#        )


     add_simple_coot_menu_menuitem(
       submenu_maps,
       "Map Density Histogram...",
       lambda func: map_molecule_chooser_gui("Choose the map",
		lambda imol: map_histogram(imol)))


     add_simple_coot_menu_menuitem(
       submenu_maps,
       "Brighten Maps",
       lambda func: brighten_maps())


     def set_diff_map_func(imol):
       print "setting map number %s to be a difference map" %imol
       set_map_is_difference_map(imol, 1)
        
     add_simple_coot_menu_menuitem(
       submenu_maps,
       "Set map is a difference map...",
       lambda func: map_molecule_chooser_gui("Which map should be considered a difference map?",
		lambda imol: set_diff_map_func(imol)))


     add_simple_coot_menu_menuitem(
       submenu_maps,
       "Another (contour) level...",
       lambda func: another_level())


     add_simple_coot_menu_menuitem(
       submenu_maps,
       "Multi-chicken...",
       lambda func: map_molecule_chooser_gui("Choose a molecule for multiple contouring",
		lambda imol: (set_map_displayed(imol, 0), multi_chicken(imol))))

     
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
         add_hydrogens_using_refmac(aa_imol)
     
     add_simple_coot_menu_menuitem(
       submenu_models,
       "Add Hydrogens",
       lambda func: add_hydrogens_with_coot_reduce())


     add_simple_coot_menu_menuitem(
       submenu_models,
       "Add Hydrogens using Refmac",
       lambda func: add_hydrogens_refmac_func())


     add_simple_coot_menu_menuitem(
       submenu_models,
       "Add Other Solvent Molecules...",
       lambda func: solvent_ligands_gui())

     
     add_simple_coot_menu_menuitem(
       submenu_models,
       "Arrange Waters Around Protein...",
       lambda func: molecule_chooser_gui(
          "Arrange waters in molecule: ",
          lambda imol: move_waters_to_around_protein(imol)))


     add_simple_coot_menu_menuitem(
       submenu_models,
       "Assign (force) HETATMs for this Residue",
       lambda func: using_active_atom(hetify_residue,
                                      "aa_imol", "aa_chain_id", "aa_res_no", "aa_ins_code"))

     
     add_simple_coot_menu_menuitem(
       submenu_models,
       "Assign HETATM to molecule...", 
       lambda func: molecule_chooser_gui("Assign HETATMs as per PDB definition", 
		lambda imol: assign_hetatms(imol)))

     # in main menu now
     # add_simple_coot_menu_menuitem(
     #   submenu_models,
     #   "Copy Coordinates Molecule...", 
     #   lambda func: molecule_chooser_gui("Molecule to Copy...", 
     #    	lambda imol: copy_molecule(imol)))


     # moved to main menu now
     # should stay open if helper function returns False
     # def atom_selection_from_fragmemt_func(imol, text, button_state):
     #   print "BL DEBUG:: imol, text, button_state", imol, text, button_state
     #   jmol = new_molecule_by_atom_selection(imol, text)
     #   if button_state:
     #     move_molecule_to_screen_centre(jmol)
     #   return valid_model_molecule_qm(jmol)
     # add_simple_coot_menu_menuitem(
     #   submenu_models,
     #   "Copy Fragment...", 
     #   lambda func: generic_chooser_and_entry_and_check_button("Create a new Molecule\n \
     #                              From which molecule shall we copy the fragment?", 
     #                                                           "Atom selection for fragment", "//A/1-10", "Move molecule here?", 
     #    	                                               lambda imol, text, button_state: atom_selection_from_fragmemt_func(imol, text, button_state),
     #                                          False))

     # --- D ---

     # now in main menu
##     add_simple_coot_menu_menuitem(
##       submenu_models,
##       "DB Loop...",
##       lambda func: click_protein_db_loop_gui())

     
     # errr... move this...(??)
     submenu = gtk.Menu()
     menuitem2 = gtk.MenuItem("Assign Sequence...")

     menuitem2.set_submenu(submenu)
     menu.append(menuitem2)
     menuitem2.show()
     
     add_simple_coot_menu_menuitem(
       submenu,
       "1: Associate Sequence to Chain...",
       lambda func: associate_sequence_with_chain_gui()) # no alignment on OK press

     add_simple_coot_menu_menuitem(
       submenu,
       "2: Assign Sequence...", 
       lambda func: cootaneer_gui_bl())

     def associate_seq_func(imol, chain_id, pir_file):
       import os, re
       chain_count = 0
       reg_chain = re.compile("chain", re.IGNORECASE)
       print "assoc seq:", imol, chain_id, pir_file
       if (os.path.isfile(pir_file)):
         fin = open(pir_file, 'r')
         seq_text = fin.read()
         fin.close()
         assign_pir_sequence(imol, chain_id, seq_text)
       else:
         print "BL WARNING:: could not find", pir_file

     # old
     #add_simple_coot_menu_menuitem(
     #  submenu,
     #  "Associate Sequence...",
     #  lambda func: generic_chooser_entry_and_file_selector(
     #    "Associate Sequence to Model: ",
     #    valid_model_molecule_qm,
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
     #  lambda func:    generic_chooser_and_file_selector(
     #                  "Associate Sequence with Molecule: ",
     #                  valid_model_molecule_qm,
     #                  "Select Sequence File",
     #                  "",
     #                  lambda imol, sequence_file_name:
     #                  assign_sequence_from_file(imol,
     #                                            sequence_file_name)
     #  ))
     

     add_simple_coot_menu_menuitem(
       submenu_models,
       "Duplicate range (pick atoms)",
       lambda func: duplicate_range_by_atom_pick())
     
     # ---- F ---------

     # doublication to entry in main gtk code!
     submenu = gtk.Menu()
     menuitem2 = gtk.MenuItem("Find Secondary Structure...")

     menuitem2.set_submenu(submenu)
     submenu_models.append(menuitem2)
     menuitem2.show()
     
     add_simple_coot_menu_menuitem(
       submenu,
       "Find Helices", 
       lambda func: find_helices())

     add_simple_coot_menu_menuitem(
       submenu,
       "Find Strands", 
       lambda func: find_strands())


     def get_smiles_pdbe_func():
       with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code, aa_atom_name, aa_alt_conf]:
         comp_id = residue_name(aa_imol, aa_chain_id,
                                aa_res_no, aa_ins_code)
         # print "BL INFO:: here with residue name", comp_id
         get_SMILES_for_comp_id_from_pdbe(comp_id)
                         
     add_simple_coot_menu_menuitem(
       submenu_models,
       "Fetch PDBe description for this ligand",
       lambda func: get_smiles_pdbe_func())


     def get_pdbe_ligand_func(comp_id):
         status = get_SMILES_for_comp_id_from_pdbe(comp_id)
         get_monomer(comp_id)

     add_simple_coot_menu_menuitem(
       submenu_models,
       "Fetch PDBe Ligand Description",
       lambda func: generic_single_entry("Fetch PDBe Ligand Desciption for comp_id:",
                                         "", " Fetch ", lambda comp_id: get_pdbe_ligand_func(comp_id)))


     add_simple_coot_menu_menuitem(
       submenu_models,
       "Fix Nomenclature Errors...",
       lambda func: molecule_chooser_gui("Fix Nomenclature Error in molecule:",
                                         lambda imol: fix_nomenclature_errors(imol)))


     # --- I --------

     def chiral_centre_inverter_func():
         with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                    aa_ins_code, aa_atom_name, aa_alt_conf]:
             invert_chiral_centre(aa_imol, aa_chain_id, aa_res_no,
                                  aa_ins_code, aa_atom_name)
             
     add_simple_coot_menu_menuitem(
         submenu_models,
         "Invert This Chiral Centre",
         lambda func:
            chiral_centre_inverter_func())
     

     # --- J --------

     add_simple_coot_menu_menuitem(
       submenu_models,
       "JLigand launch",
       lambda func:
          launch_jligand_function())
     

     # --- M ---

     def make_link_ext_func(*args):
       m_spec_1 = args[0]
       m_spec_2 = args[1]
       imol_1 = atom_spec_to_imol(m_spec_1)
       imol_2 = atom_spec_to_imol(m_spec_2)
       spec_1 = m_spec_1[2:]
       spec_2 = m_spec_2[2:]
       if not (imol_1 == imol_2):
         print "Mismatch molecules"
       else:
         make_link(imol_1, spec_1, spec_2, "dummy", 0.1)
       
     add_simple_coot_menu_menuitem(
       submenu_models,
       "Make Link (click 2 atoms)...",
       lambda func:
       user_defined_click(2, make_link_ext_func))

       
     add_simple_coot_menu_menuitem(
       submenu_models,
       "Merge Water Chains...",
       lambda func: molecule_chooser_gui("Merge Water Chains in molecule:",
                                         lambda imol: merge_solvent_chains(imol)))


     def mon_dict_func(text):
       idealized = 0
       new_model = get_monomer_from_dictionary(text, idealized)
       if not valid_model_molecule_qm(new_model):
         get_monomer(text)
         
     add_simple_coot_menu_menuitem(
       submenu_models,
       "Monomer from Dictionary",
       lambda func:
         generic_single_entry("Pull coordinates from CIF dictionary for 3-letter-code:", "",
                              " Get Coords ",
                              lambda text: mon_dict_func(text)))
                              

     def morph_fit_chain_func(radius):
       with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                  aa_ins_code, aa_atom_name, aa_alt_conf]:
         morph_fit_chain(aa_imol, aa_chain_id, radius)
         
     add_simple_coot_menu_menuitem(
       submenu_models,
       "Morph Fit Chain (Averaging Radius 7)",
       lambda func: morph_fit_chain_func(7)
       )
     
     add_simple_coot_menu_menuitem(
       submenu_models,
       "Morph Fit Chain (Averaging Radius 11)",
       lambda func: morph_fit_chain_func(11)
       )


     def morph_fit_ss_func(radius):
       with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no,
                                  aa_ins_code, aa_atom_name, aa_alt_conf]:
         morph_fit_by_secondary_structure_elements(aa_imol, aa_chain_id)
         
     add_simple_coot_menu_menuitem(
       submenu_models,
       "Morph Fit Chain based on Secondary Structure",
       lambda func: morph_fit_chain_func(7)
       )
     
     
     # -- N --

     def new_mol_sphere_func1(imol, text):
       try:
         radius = float(text)
       except:
         print "WARNING:: no valid radius", text
       args = [imol] + rotation_centre() + [radius, 0]
       new_molecule_by_sphere_selection(*args)

     add_simple_coot_menu_menuitem(
       submenu_models,
       "New Molecule by Sphere...",
       lambda func: generic_chooser_and_entry(
          "Choose a molecule from which to select a sphere of atoms:",
          "Radius:", "10.0",
          lambda imol, text: new_mol_sphere_func1(imol, text)))


     def new_mol_sym_func1(imol, text):
       from types import ListType
       pre_shift = origin_pre_shift(imol)
       if (type(pre_shift) is not ListType):
         print "bad pre-shift aborting"
       else:
         new_molecule_by_symop(imol, text,
                               pre_shift[0],
                               pre_shift[1],
                               pre_shift[2])
       
     add_simple_coot_menu_menuitem(
       submenu_models,
       "New Molecule from Symmetry Op...",
       lambda func: generic_chooser_and_entry(
          "Molecule from which to generate a symmetry copy",
          "SymOp", "X,Y,Z",
          lambda imol, text: new_mol_sym_func1(imol, text)))


     # --- P ---
     
     add_simple_coot_menu_menuitem(
       submenu_models,
       "Phosphorylate this residue",
       lambda func: phosphorylate_active_residue())


     # # FIXME:: for now only when prodrg in path
     # if (command_in_path_qm("cprodrg")):
     #   add_simple_coot_menu_menuitem(
     #   submenu_models,
     #   "Prodrg-ify this residue (generate restraints)",
     #   lambda func: using_active_atom(prodrg_ify,
     #                                  "aa_imol", "aa_chain_id", "aa_res_no", "aa_ins_code"))
     

     # ---- R ---------

     # --- Ren --- 

     add_simple_coot_menu_menuitem(
       submenu_models,
       "Rename Residue...",
       lambda func: rename_residue_gui())


     add_simple_coot_menu_menuitem(
       submenu_models,
       "Renumber Waters...",
       lambda func: molecule_chooser_gui(
          "Renumber waters of which molecule?",
          lambda imol: renumber_waters(imol)))

     # --- Reo ---

     add_simple_coot_menu_menuitem(
       submenu_models,
       "Reorder Chains...",
       lambda func: molecule_chooser_gui("Sort Chain IDs in molecule:",
                                         lambda imol: sort_chains(imol))) # an internal function

     # --- Rep ---

     # BL says:: may work, not sure about function entirely
     # this is in main menu now
     # add_simple_coot_menu_menuitem(
     #   submenu_models,
     #   "Replace Fragment...",
     #   lambda func: molecule_chooser_gui("Define the molecule that needs updating",
     #    	lambda imol_base: generic_chooser_and_entry(
     #    			"Molecule that contains the new fragment:",
     #    			"Atom Selection","//",
     #    			lambda imol_fragment, atom_selection_str:
     #    			replace_fragment(imol_base, imol_fragment, atom_selection_str))))

     # in main menu now
     # add_simple_coot_menu_menuitem(
     #   submenu_models,
     #   "Replace Residue...",
     #   lambda func: generic_single_entry("Replace this residue with residue of type:",
     #                                     "ALA", "Mutate",
     #                                     lambda text: using_active_atom(mutate_by_overlap,
     #                                                                    "aa_imol", "aa_chain_id", "aa_res_no",
     #                                                                    text)))

     # --- Res ---
     
     add_simple_coot_menu_menuitem(
       submenu_models,
       "Residue Type Selection...",
	lambda func: generic_chooser_and_entry("Choose a molecule to select residues from: ","Residue Type:","",
                                               lambda imol, text: (new_molecule_by_residue_type_selection(imol, text),
                                                                   update_go_to_atom_window_on_new_mol())))


     add_simple_coot_menu_menuitem(
       submenu_models,
       "Residues with Alt Confs...",
       lambda func: molecule_chooser_gui(
         "Which molecule to check for Alt Confs?",
         lambda imol: alt_confs_gui(imol)))


     add_simple_coot_menu_menuitem(
       submenu_models,
       "Residues with Cis Peptide Bonds...",
       lambda func: molecule_chooser_gui("Choose a molecule for checking for Cis Peptides",
                                         lambda imol: cis_peptides_gui(imol)))


     add_simple_coot_menu_menuitem(
       submenu_models,
       "Residues with Missing Atoms...",
       lambda func: molecule_chooser_gui(
         "Which molecule to check for Missing Atoms?",
         lambda imol: missing_atoms_gui(imol)))

     # --- Rig ---

     add_simple_coot_menu_menuitem(
       submenu_models,
       "Rigid Body Fit Residue Ranges...",
       lambda func:
       residue_range_gui(lambda imol, ls: rigid_body_refine_by_residue_ranges(imol, ls),
                         "Rigid Body Refine",
                         "  Fit  "))

     add_simple_coot_menu_menuitem(
       submenu_models,
       "Rigid Body Fit Molecule...",
       lambda func: molecule_chooser_gui("Rigid Body Fit Molecule",
                lambda imol: rigid_body_refine_by_atom_selection(imol, "//")))                                         
       
     # ---- S --------

     add_simple_coot_menu_menuitem(
       submenu_models,
       "Superpose ligands",
       lambda func: superpose_ligand_gui())
     

     add_simple_coot_menu_menuitem(
       submenu_models,
       "Symm Shift Reference Chain Here",
       lambda func: move_reference_chain_to_symm_chain_position())

     
     # ---- U ---------
     
     submenu = gtk.Menu()
     menuitem2 = gtk.MenuItem("Rotamer Search...")
     
     menuitem2.set_submenu(submenu)
     submenu_models.append(menuitem2)
     menuitem2.show()

     add_simple_coot_menu_menuitem(
       submenu,
       "Use \"Backrub\" Rotamers",
       lambda func: set_rotamer_search_mode(ROTAMERSEARCHLOWRES)
       )
     
     add_simple_coot_menu_menuitem(
       submenu,
       "DONT use \"Backrub\" Rotamers",
       lambda func: set_rotamer_search_mode(ROTAMERSEARCHHIGHRES)
       )


     add_simple_coot_menu_menuitem(
       submenu_models,
       "Use SEGIDs...",
       lambda func: molecule_chooser_gui("Exchange the Chain IDs, replace with SEG IDs",
		lambda imol: exchange_chain_ids_for_seg_ids(imol)))

     
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
       
     add_simple_coot_menu_menuitem(
       submenu_models,
       "What's this?",
       lambda func: whats_this()
       )
     
     
     # an python extra, this is
     add_simple_coot_menu_menuitem(
       submenu_models,
       "Atoms with Zero Occupancies...",
       lambda func: molecule_chooser_gui(
         "Which molecule to check for Atoms with zero occupancies?",
         lambda imol: zero_occ_atoms_gui(imol)))

     
     #---------------------------------------------------------------------
     #     NCS functions
     #
     #---------------------------------------------------------------------

     def copy_ncs_range_func(imol, chain_id, text1, text2):
       try:
         r1 = int(text1)
         r2 = int(text2)
         if (chain_id == ncs_master_chain_id(imol)):
           # hunkey dorey
           copy_residue_range_from_ncs_master_to_others(imol, chain_id, r1, r2)
         else:
           # different given master to current master
           # ask what to do.
           txt = "Current master chain is %s, but you asked to copy from %s.\n" \
                 %(ncs_master_chain_id(imol), chain_id)
           txt += "Apply this master change?\n\n"
           txt += "N.B. if no, then nothing is copied."
           r = yes_no_dialog(txt, "Change Master")
           if r:
             ncs_control_change_ncs_master_to_chain_id(imol, chain_id)
             copy_residue_range_from_ncs_master_to_others(imol, chain_id, r1, r2)
             ## could change master back?!
           else:
             info_dialog("Master chain was not changed and copy not applied.")
       except:
         print "BL WARNING:: no valid number input"

     add_simple_coot_menu_menuitem(
       submenu_ncs,
       "Copy NCS Residue Range...",
       lambda func: generic_chooser_and_entry("Apply NCS Range from Master",
                                              "Master Chain ID",
                                              get_first_ncs_master_chain(),  # returns "" on fail
                lambda imol, chain_id: generic_double_entry("Start of Residue Number Range",
                                       "End of Residue Number Range",
                                       "", "", False, False,
                                       "Apply NCS Residue Range",
                                       lambda text1, text2: copy_ncs_range_func(imol, chain_id, text1, text2))))


     def copy_ncs_chain_func(imol, chain_id):
       ncs_chains = ncs_chain_ids(imol)
       if (ncs_chains):
         # maybe this could be a function to avoid repetition.
         if (chain_id == ncs_master_chain_id(imol)):
           # hunkey dorey
           copy_from_ncs_master_to_others(imol, chain_id)
         else:
           # different given master to current master
           # ask what to do.
           txt = "Current master chain is %s, but you asked to copy from %s.\n" \
                 %(ncs_master_chain_id(imol), chain_id)
           txt += "Apply this master change?\n\n"
           txt += "N.B. if no, then nothing is copied."
           r = yes_no_dialog(txt, "Change Master")
           if r:
             ncs_control_change_ncs_master_to_chain_id(imol, chain_id)
             copy_from_ncs_master_to_others(imol, chain_id)
             ## could change master back?!
           else:
             info_dialog("Master chain was not changed and copy not applied.")
       else:
         s = "You need to define NCS operators for molecule " + str(imol)
         info_dialog(s)
           
     add_simple_coot_menu_menuitem(
       submenu_ncs,
       "Copy NCS Chain...",
       lambda func: generic_chooser_and_entry("Apply NCS edits from NCS Master Chain to Other Chains",
                                              "Master Chain ID",
                                              get_first_ncs_master_chain(),  # can return  "".
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
           print "BL WARNING:: input residue numbers have to be integers"
         if (cont):
           ghost_ncs_chain_ids = ncs_chain_ids(imol)
           if (type(ghost_ncs_chain_ids) is ListType):
             # because we can have hetero-NCS,
             # but we ignore NCS other that
             # that of the first type.
             ghost_chain_list = ghost_ncs_chain_ids[0]
             manual_ncs_ghosts(imol, resno_1, resno_2, ghost_chain_list)
         
       generic_double_entry(label_1, label_2, entry_1_default_text,
                            entry_2_default_text, False, False,
                            go_button_label, handle_go_function)
       

     add_simple_coot_menu_menuitem(
       submenu_ncs,
       "NCS Ghosts by Residue Range...",
       lambda func: molecule_chooser_gui("Make local NCS ghosts for molecule:",
                                         lambda imol: ncs_ghost_res_range_func(imol)))

     add_simple_coot_menu_menuitem(
       submenu_ncs,
       "Update NCS Ghosts using Local Match",
       lambda func: update_ncs_ghosts_by_local_sphere())


     add_simple_coot_menu_menuitem(
       submenu_ncs,
       "NCS Jumping...",
       lambda func: ncs_jumping_gui())


     add_simple_coot_menu_menuitem(
       submenu_ncs,
       "NCS ligands...",
       lambda func: ncs_ligand_gui())

     submenu = gtk.Menu()
     menuitem2 = gtk.MenuItem("NCS matrix type...")
     
     menuitem2.set_submenu(submenu)
     submenu_ncs.append(menuitem2)
     #tooltip = gtk.Tooltips()
     if gtk.pygtk_version >= (2,12):
       menuitem2.set_tooltip_text("use to change the way the NCS matrix is calculated")
     else:
       coot_tooltips.set_tip(menuitem2, "use to change the way the NCS matrix is calculated")
     menuitem2.show()

     add_simple_coot_menu_menuitem(
       submenu,
       "Accurate (SSM)",
       lambda func: set_ncs_matrix_type(0))


     add_simple_coot_menu_menuitem(
       submenu,
       "Fast (LSQ)",
       lambda func: set_ncs_matrix_type(1))


     add_simple_coot_menu_menuitem(
       submenu,
       "Extra Fast (LSQ, only every 2nd CA)",
       lambda func: set_ncs_matrix_type(2))


     # ---------------------------------------------------------------------
     #     Refinement
     # ---------------------------------------------------------------------
     #

     add_simple_coot_menu_menuitem(
       submenu_refine,
       "Set Refinement Options...",
       lambda func: refinement_options_gui())
 
     submenu = gtk.Menu()
     menuitem2 = gtk.MenuItem("Peptide Restraints...")

     menuitem2.set_submenu(submenu)
     submenu_refine.append(menuitem2)
     menuitem2.show()

     def add_restr_func1():
	print 'Planar Peptide Restraints added'
	add_planar_peptide_restraints()

     add_simple_coot_menu_menuitem(
       submenu,
       "Add Planar Peptide Restraints",
       lambda func: add_restr_func1())


     def add_restr_func2():
	print 'Planar Peptide Restraints removed'
	remove_planar_peptide_restraints()

     add_simple_coot_menu_menuitem(
       submenu, "Remove Planar Peptide Restraints",
       lambda func: add_restr_func2())


     def shelx_ref_func():
       window = gtk.Window(gtk.WINDOW_TOPLEVEL)
       vbox = gtk.VBox(False, 0)
       hbox = gtk.HBox(False, 0)
       go_button = gtk.Button("  Refine  ")
       cancel_button = gtk.Button("  Cancel  ")
       entry_hint_text = "HKL data filename \n(leave blank for default)"
       chooser_hint_text = " Choose molecule for SHELX refinement  "
       h_sep = gtk.HSeparator()

       window.add(vbox)
       option_menu_mol_list_pair = generic_molecule_chooser(vbox, chooser_hint_text)
       entry = file_selector_entry(vbox, entry_hint_text)

       def shelx_delete_event(*args):
         window.destroy()
         return False

       def shelx_go_funcn_event(*args):
         import operator
         txt = entry.get_text()
         imol = get_option_menu_active_molecule(*option_menu_mol_list_pair)
         if (operator.isNumberType(imol)):
           if (len(txt) == 0):
             shelxl_refine(imol)
           else:
             shelxl_refine(imol, txt)
         window.destroy()
         return False

       go_button.connect("clicked", shelx_go_funcn_event)
       cancel_button.connect("clicked", shelx_delete_event)

       vbox.pack_start(h_sep, False, False, 2)
       vbox.pack_start(hbox, False, False, 2)
       hbox.pack_start(go_button, True, False, 0)
       hbox.pack_start(cancel_button, True, False, 0)
       window.show_all()
       
     add_simple_coot_menu_menuitem(
       submenu_refine,
       "SHELXL Refine...", 
       lambda func: shelx_ref_func())


     add_simple_coot_menu_menuitem(
       submenu_refine,
       "Read REFMAC logfile...",
       lambda func: generic_chooser_and_file_selector("Read Refmac log file",
                       valid_model_molecule_qm, "Logfile name: ", "",
                       lambda imol, text: read_refmac_log(imol, text)))
                       
     add_simple_coot_menu_menuitem(
       submenu_refine,
       "Occupancy refinement input for REFMAC...",
       lambda func: generic_chooser_and_file_selector("Extra restraints file",
                       valid_model_molecule_qm, "Restraints file name: ",
                       "refmac_extra_params.txt",
                       lambda imol, text: restraints_for_occupancy_refinement(imol, text)))
       

     # An example with a submenu:
     #
     submenu = gtk.Menu()
     menuitem2 = gtk.MenuItem("Refinement Speed...")
     menuitem2.set_submenu(submenu)
     submenu_refine.append(menuitem2)
     menuitem2.show()

     add_simple_coot_menu_menuitem(
       submenu,
       "Molasses Refinement mode", 
       lambda func: (printf("Molasses..."),
                     set_dragged_refinement_steps_per_frame(4)))


     add_simple_coot_menu_menuitem(
       submenu,
       "Smooth Refinement mode",
       lambda func: (set_dragged_refinement_steps_per_frame(42)))


     add_simple_coot_menu_menuitem(
       submenu,
       "Crocodile Refinement mode", 
       lambda func: (printf("Crock..."),
                     set_dragged_refinement_steps_per_frame(220)))


     add_simple_coot_menu_menuitem(
       submenu,
       "Normal Refinement mode (1 Emsley)", 
       lambda func: (printf("Default Speed (1 Emsley)..."),
                     set_dragged_refinement_steps_per_frame(140)))


     add_simple_coot_menu_menuitem(
       submenu_refine,
       "Auto-weight refinement",
       lambda func: auto_weight_for_refinement())
     

     add_simple_coot_menu_menuitem(
       submenu_refine,
       "Set Undo Molecule...",
       lambda func: molecule_chooser_gui("Set the Molecule for 'Undo' Operations",
                    lambda imol: set_undo_molecule(imol)))


     # BL says: has no checking for text = number yet
     add_simple_coot_menu_menuitem(
       submenu_refine,
       "B factor bonds scale factor...",
       lambda func: generic_chooser_and_entry("Choose a molecule to which the B-factor colour scale is applied:",
		"B-factor scale:", "1.0", 
		lambda imol, text: set_b_factor_bonds_scale_factor(imol,float(text))))


     def set_mat_func(text):
	import operator
	t = float(text)
	if operator.isNumberType(t):
		s = "Matrix set to " + text
		set_matrix(t)
		add_status_bar_text(s)
	else:
		add_status_bar_text("Failed to read a number")

     add_simple_coot_menu_menuitem(
       submenu_refine,
       "Set Matrix (Refinement Weight)...",
       lambda func: generic_single_entry("set matrix: (smaller means better geometry)", 
		str(matrix_state()), "Set it", 
		lambda text: set_mat_func(text)))


     def set_den_gra_func(text):
	import operator
	t = float(text)
	if operator.isNumberType(t):
		s = "Density Fit scale factor set to " + text
		set_residue_density_fit_scale_factor(t)
		add_status_bar_text(s)
	else:
		add_status_bar_text("Failed to read a number")



#      # ---------------------------------------------------------------------
#      #     Recent structures from the PDBe
#      # ---------------------------------------------------------------------
#      #
#      add_simple_coot_menu_menuitem(
#        submenu_pdbe, "PDBe recent structures...",
#        lambda func: pdbe_latest_releases_gui())

#      # we do test for refmac at startup not runtime (for simplicity)
#      if command_in_path_qm("refmac5"):
#        mess = " Get it "
#      else:
#        mess = "\n  WARNING::refmac5 not in the path - SF calculation will fail  \n\n"
       
#      add_simple_coot_menu_menuitem(
#        submenu_pdbe, "Get from PDBe...",
#        lambda func: generic_single_entry("Get PDBe accession code",
#                                          "", " Get it ",
#                                          lambda text:
#                                          pdbe_get_pdb_and_sfs_cif("include-sfs", text.rstrip().lstrip())))

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
#            read_pdb(pdb_file_name)
#          if os.path.isfile(mtz_file_name):
#            make_and_draw_map(mtz_file_name, "FWT", "PHWT", "", 0, 0)
#            make_and_draw_map(mtz_file_name, "DELFWT", "PHDELWT", "", 0, 1)
       
#      add_simple_coot_menu_menuitem(
#        menu,
#        "Load tutorial model and data",
#        lambda func: load_tutorial_data_func()
#        )

     
     # ---------------------------------------------------------------------
     #     Views/Representations
     # ---------------------------------------------------------------------
     #

     add_simple_coot_menu_menuitem(
       submenu_representation,
       "Undo Symmetry View",
       lambda func: undo_symmetry_view())


     def make_ball_n_stick_func(imol, text):
       bns_handle = make_ball_and_stick(imol, text, 0.18, 0.3, 1)
       print "handle: ", bns_handle


     global default_ball_and_stick_selection    # maybe should be at the top of the file
     add_simple_coot_menu_menuitem(
       submenu_representation,
       "Ball & Stick...",
       lambda func: generic_chooser_and_entry("Ball & Stick",
                                              "Atom Selection:",
                                              default_ball_and_stick_selection,
                                              lambda imol, text: make_ball_n_stick_func(imol, text)))

     add_simple_coot_menu_menuitem(
       submenu_representation,
       "Add Balls to Simple Sticks",
       lambda func: map(lambda imol: set_draw_stick_mode_atoms(imol, 1), molecule_number_list()))


     add_simple_coot_menu_menuitem(
       submenu_representation,
       "Simple Sticks (No Balls)",
       lambda func: map(lambda imol: set_draw_stick_mode_atoms(imol, 0), molecule_number_list()))


     add_simple_coot_menu_menuitem(
       submenu_representation,
       "Clear Ball & Stick...",
       lambda func: molecule_chooser_gui(
         "Choose a molecule from which to clear Ball&Stick objects",
         lambda imol: clear_ball_and_stick(imol)))


     add_simple_coot_menu_menuitem(
       submenu_representation,
       "Electrostatic Surface...",
       lambda func: molecule_chooser_gui(
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
       
     add_simple_coot_menu_menuitem(
       submenu_representation,
       "Clipped Surface Here (This Residue)",
       lambda func:
         surface_func1(1))

     add_simple_coot_menu_menuitem(
       submenu_representation,
       "Full Surface Around Here (This Residue)",
       lambda func:
         surface_func1())
     

     add_simple_coot_menu_menuitem(
       submenu_representation,
       "Un-Surface...",
       lambda func: molecule_chooser_gui(
          "Choose a molecule to represent conventionally...",
          lambda imol: do_surface(imol, 0)))

     def hilight_site_func():
       active_atom = active_residue()
       if (active_atom):
         imol = active_atom[0]
         centre_residue_spec = [active_atom[1],
                                active_atom[2],
                                active_atom[3]]
         hilight_binding_site(imol, centre_residue_spec, 230,4)

     add_simple_coot_menu_menuitem(
       submenu_representation,
       "Highlight Interesting Site (here)...",
       lambda func: hilight_site_func())
     

     def make_dot_surf_func(imol,text):
        # I think a single colour is better than colour by atom
        set_dots_colour(imol, 0.5, 0.5, 0.5)
	dots_handle = dots(imol, text, text, 2, 1)
	print "dots handle: ", dots_handle

     add_simple_coot_menu_menuitem(
       submenu_representation,
       "Dotted Surface...",
       lambda func: generic_chooser_and_entry("Surface for molecule", 
		"Atom Selection:", "//A/1-2", 
		lambda imol, text: make_dot_surf_func(imol, text)))


     def clear_dot_surf_func(imol,text):
	try:
          n = int(text)
          clear_dots(imol,n)
	except:
          print "BL WARNING:: dots handle number shall be an integer!!"

     add_simple_coot_menu_menuitem(
       submenu_representation,
       "Clear Surface Dots...",
       lambda func: generic_chooser_and_entry("Molecule with Dotted Surface", 
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
         
     add_simple_coot_menu_menuitem(
       submenu_representation,
       "Limit Model Display Radius...",
       lambda func: generic_single_entry("Display Radius Limit (0 for \'no limit\') ",
                                         #  "15.0" ;; maybe this should be the map radius
                                         # BL says:: I think it should be the current one
                                         str(get_map_radius()),
                                         "Set: ",
                                         lambda text: limit_model_disp_func(text)))

     
     add_simple_coot_menu_menuitem(
         submenu_representation,
         "HOLE...",
         lambda func: hole_ify())

     add_simple_coot_menu_menuitem(
       submenu_representation,
       "Label All CAs...",
       lambda func: molecule_chooser_gui("Choose a molecule to label",
                                         lambda imol: label_all_CAs(imol)))
 
     # Views submenu
     submenu = gtk.Menu()
     menuitem2 = gtk.MenuItem("Views")
 
     menuitem2.set_submenu(submenu)
     menu.append(menuitem2)
     menuitem2.show()

     add_simple_coot_menu_menuitem(
       submenu,
       "Add View...",
       lambda func: view_saver_gui())

     # BL says:: maybe check if number at some point
     add_simple_coot_menu_menuitem(
       submenu,
       "Add a Spin View...",
       lambda func: generic_double_entry("Number of Steps", 
                         "Number of Degrees (total)", "3600", "360", 
			 False, False, 		#check button text and callback
			 "  Add Spin  ",
                         lambda text_1, text_2: add_spin_view("Spin", int(text_1), float(text_2))))

     add_simple_coot_menu_menuitem(
       submenu,
       "Views Panel...",
       lambda func: views_panel_gui())
 
     add_simple_coot_menu_menuitem(
       submenu,
       "Play Views",
       lambda func: map(eval,["go_to_first_view(1)",
                              "time.sleep(1)", "play_views()"]))
 
     # BL says:: maybe check if number at some point
     add_simple_coot_menu_menuitem(
       submenu, "Set Views Play Speed...",
       lambda func: generic_single_entry("Set Views Play Speed",
			str(views_play_speed()), "  Set it  ",
			lambda text: set_views_play_speed(float(text))))


     add_simple_coot_menu_menuitem(
       submenu, "Save Views...",
       lambda func: generic_single_entry("Save Views",
                                         "coot-views.py", " Save ",
                                         lambda txt: save_views(txt)))


     #---------------------------------------------------------------------
     #     3D annotations
     #---------------------------------------------------------------------

     submenu = gtk.Menu()
     menuitem2 = gtk.MenuItem("3D Annotations...")
 
     menuitem2.set_submenu(submenu)
     submenu_representation.append(menuitem2)
     menuitem2.show()

     add_simple_coot_menu_menuitem(
       submenu,
       "Annotate position...",
       lambda func: generic_single_entry("Annotation: ", "",
                                         "Make Annotation",
                                         lambda txt: add_annotation_here(txt)))


     # BL says:: maybe this (and next) should have a file chooser/selector!?
     add_simple_coot_menu_menuitem(
       submenu,
       "Save Annotations...",
       lambda func: generic_single_entry("Save Annotations",
                                         "coot_annotations.py",
                                         " Save ",
                                         lambda file_name: save_annotations(file_name)))
       

     add_simple_coot_menu_menuitem(
       submenu,
       "Load Annotations...",
       lambda func: generic_single_entry("Load Annotations",
                                         "coot_annotations.py",
                                         " Load ",
                                         lambda file_name: load_annotations(file_name)))

     
     add_simple_coot_menu_menuitem(
       submenu,
       "Remove annotation here",
       lambda func: remove_annotation_here())

     
     add_simple_coot_menu_menuitem(
       submenu,
       "Remove annotation near click",
       lambda func: remove_annotation_at_click())
     

     #---------------------------------------------------------------------
     #     Other Representation Programs
     #
     #---------------------------------------------------------------------

     # shall use subprocess at some point
     def ccp4mg_func1():
	import os
	pd_file_name = "1.mgpic.py"
	write_ccp4mg_picture_description(pd_file_name)
	if os.name == 'nt':
          ccp4mg_exe = "winccp4mg.exe"
	else:
          ccp4mg_exe = "ccp4mg"
	if command_in_path_qm(ccp4mg_exe):
          ccp4mg_file_exe = find_exe(ccp4mg_exe, "PATH")
          pd_file_name = os.path.abspath(pd_file_name)
          args = [ccp4mg_file_exe, "-pict", pd_file_name]
          try:
            import subprocess
            subprocess.Popen(args).pid
            print "BL DEBUG:: new subprocess"
          except:
            # no subprocess, use old style
            os.spawnv(os.P_NOWAIT, ccp4mg_file_exe, args)
	else:
          print "BL WARNING:: sorry cannot find %s in $PATH" %ccp4mg_exe

     add_simple_coot_menu_menuitem(
       submenu_representation, "CCP4MG...",
       lambda func: ccp4mg_func1())


     # ---------------------------------------------------------------------
     #     PISA Interface and Assemblies
     # ---------------------------------------------------------------------

     add_simple_coot_menu_menuitem(
       submenu_pisa, "PISA Assemblies...",
       lambda func:
       molecule_chooser_gui("Choose molecule for PISA assembly analysis",
                            lambda imol:
                            pisa_assemblies(imol)))

     add_simple_coot_menu_menuitem(
       submenu_pisa, "PISA Interfaces...",
       lambda func:
       molecule_chooser_gui("Choose molecule for PISA interface analysis",
                            lambda imol:
                            pisa_interfaces(imol)))


     # ---------------------------------------------------------------------
     #     Modules
     # ---------------------------------------------------------------------

     add_simple_coot_menu_menuitem(
       submenu_modules, "CCP4...",
       lambda func: add_module_ccp4())

     add_simple_coot_menu_menuitem(
       submenu_modules, "SHELX...",
       lambda func: add_module_shelx())

     add_simple_coot_menu_menuitem(
       submenu_modules, "User-defined Restraints...",
       lambda func: add_module_user_defined_restraints())

     add_simple_coot_menu_menuitem(
         submenu_modules, "ProSMART",
         lambda func: add_module_prosmart())

     add_simple_coot_menu_menuitem(
         submenu_modules, "Carbohydrate",
         lambda func: add_module_carbohydrate_gui())
     
     add_simple_coot_menu_menuitem(
         submenu_modules, "Cryo-EM",
         lambda func: add_module_cryo_em())

     
     # ---------------------------------------------------------------------
     #     Settings
     # ---------------------------------------------------------------------

     submenu = gtk.Menu()
     menuitem2 = gtk.MenuItem("Rotate Translate Zone Mode...")
 
     menuitem2.set_submenu(submenu)
     submenu_settings.append(menuitem2)
     menuitem2.show()

     add_simple_coot_menu_menuitem(
       submenu, "Rotate About Fragment Centre",
       lambda func: set_rotate_translate_zone_rotates_about_zone_centre(1))


     add_simple_coot_menu_menuitem(
       submenu, "Rotate About Second Clicked Atom",
       lambda func: set_rotate_translate_zone_rotates_about_zone_centre(0))


     add_simple_coot_menu_menuitem(
       submenu_settings,
       "Set Density Fit Graph Weight...",
       lambda func: generic_single_entry("set weight (smaller means apparently better fit)",
		str("%.2f" %residue_density_fit_scale_factor()), "Set it",
		lambda text: set_den_gra_func(text)))

     # BL says:: maybe check if number at some point
     add_simple_coot_menu_menuitem(
       submenu_settings,
       "Set Spin Speed",
       lambda func: generic_single_entry("Set Spin Speed (smaller is slower)",
			str(idle_function_rotate_angle()), "Set it",
			lambda text: set_idle_function_rotate_angle(float(text))))


     add_simple_coot_menu_menuitem(
       submenu_settings, "Nudge Centre...",
       lambda func: nudge_screen_centre_gui())


     def all_mol_symm_func():
	for imol in molecule_number_list():
		if valid_model_molecule_qm(imol):
			set_symmetry_whole_chain(imol, 1)

     add_simple_coot_menu_menuitem(
       submenu_settings, "All Molecules use \"Near Chains\" Symmetry", 
       lambda func: map(lambda imol: valid_model_molecule_qm(imol) and
                         set_symmetry_whole_chain(imol, 1),
                         molecule_number_list()))


     add_simple_coot_menu_menuitem(
       submenu_refine, "Question Accept Refinement", 
       lambda func: set_refinement_immediate_replacement(0))


     add_simple_coot_menu_menuitem(
       submenu_settings, "Save Graphics Size and Positions",
       lambda func: graphics_window_size_and_position_to_preferences())


     def save_dialog_func():
       post_model_fit_refine_dialog()
       post_go_to_atom_window()

       def delete_event(*args):
         window.destroy()
         return False
       
       window = gtk.Window(gtk.WINDOW_TOPLEVEL)
       label_text = "   When happy, press \"Save\" to save   \n" + "   dialog positions"
       label = gtk.Label(label_text)
       h_sep = gtk.HSeparator()
       cancel_button = gtk.Button("  Cancel  ")
       go_button = gtk.Button("  Save  ")
       vbox = gtk.VBox(False, 4)
       hbox = gtk.HBox(False, 4)

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

     add_simple_coot_menu_menuitem(
       submenu_settings, "Save Dialog Positions...",
       lambda func: save_dialog_func())


     add_simple_coot_menu_menuitem(
       submenu_settings, "Key Bindings...",
       lambda func: key_bindings_gui())

     def install_and_show_key_bindings():
          file_to_preferences("template_key_bindings.py") # copy and evaluate
          key_bindings_gui()
       
     add_simple_coot_menu_menuitem(
       submenu_settings, "Python: Install Template Keybindings",
       lambda func: install_and_show_key_bindings())

     def quick_save_func(txt):
       try:
         n = int(txt)
         gobject.timeout_add(1000*n, quick_save)
       except:
         print "BL INFO:: could not add timer for auto save!"

     add_simple_coot_menu_menuitem(
       submenu_settings, "Enable Quick-Save checkpointing...",
       lambda func:
          generic_single_entry("Checkpoint interval (seconds)",
                               "30",
                               " Start Auto-saving ",
                               lambda txt:
                                  quick_save_func(txt)))

     # Doesnt seem to be working right currently, so comment out?! Not any more?!
     # add to validate menu
     menu = coot_menubar_menu("Validate")

     add_simple_coot_menu_menuitem(menu, "Pukka Puckers...?",
                                   lambda func: molecule_chooser_gui(
       "Choose a molecule for ribose pucker analysis",
       lambda imol: pukka_puckers_qm(imol)))

     add_simple_coot_menu_menuitem(
       menu,
       "Alignment vs PIR...",
       lambda func: molecule_chooser_gui("Alignment vs PIR info for molecule:",
                                         lambda imol: wrapper_alignment_mismatches_gui(imol)))
          

  else:
	print "BL WARNING:: could not find the main_menubar! Sorry, no extensions menu!"

