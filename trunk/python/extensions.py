
# Copyright 2007, 2008 by Bernhard Lohkamp
# Copyright 2006, 2007, 2008 by The University of York

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

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
   
have_coot_python = False
try: 
  import coot_python
  have_coot_python = True
except:
  print """BL WARNING:: could not import coot_python module!!
Some things, esp. extensions, may be crippled!"""

if (have_coot_python):
  if coot_python.main_menubar():

     # ---------------------------------------------
     #           coot news dialog
     # ---------------------------------------------
     menu = coot_menubar_menu("About")
     if (menu):
       add_simple_coot_menu_menuitem(menu, "Coot News",
                                     lambda func: whats_new_dialog())

     
     # ---------------------------------------------
     #           extensions
     # ---------------------------------------------

     menu = coot_menubar_menu("E_xtensions")

     #---------------------------------------------------------------------
     #     Post MR
     #
     #---------------------------------------------------------------------

     add_simple_coot_menu_menuitem(menu,"[Post MR] Fill Partial Residues...",
	lambda func: molecule_chooser_gui("Find and Fill residues with missing atoms",
		lambda imol: fill_partial_residues(imol)))


     add_simple_coot_menu_menuitem(menu,"[Post MR] Fit Protein...",
	lambda func: molecule_chooser_gui("Fit Protein using Rotamer Search",
		lambda imol: (imol_refinement_map() == -1 and
                              add_status_bar_text("oops. Must set a map to fit") or
                              fit_protein(imol))[0]))                                          


     add_simple_coot_menu_menuitem(menu,"[Post MR] Stepped Refine...",
	lambda func: molecule_chooser_gui("Stepped Refine: ",
		lambda imol: (imol_refinement_map() == -1 and
                              add_status_bar_text("oops. Must set a map to fit") or
                              stepped_refine_protein(imol))[0]))


     add_coot_menu_separator(menu)

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
       
     add_simple_coot_menu_menuitem(menu,"Mask Map by Atom Selection...",
         lambda func: molecule_chooser_gui("Define the molecule that has atoms to mask the map",
             lambda imol: generic_multiple_entries_with_check_button(
                         [[" Map molecule number: ", mask_map_func()],
                          [" Atom selection: ", "//A/1"],
                          ["Radius around atoms: ", mask_map_radius_func()]],
                         [" Invert Masking? ", lambda active_state: mask_map_func1(active_state)],
                         "  Mask Map  ", lambda text_list, invert: mask_map_func2(imol, text_list, invert))))


     add_simple_coot_menu_menuitem(menu,"Copy Map Molecule...", 
	lambda func: map_molecule_chooser_gui("Map to Copy...", 
		lambda imol: copy_molecule(imol)))


     add_simple_coot_menu_menuitem(menu, "Make a Difference Map...",
                                   lambda func: make_difference_map_gui())

     add_simple_coot_menu_menuitem(menu, "Transform map by LSQ model fit...",
                                   lambda func: transform_map_using_lsq_matrix_gui())


     add_simple_coot_menu_menuitem(menu, "Export map...",
        lambda func: generic_chooser_and_file_selector("Export Map: ",
                valid_map_molecule_qm, "File_name: ", "",
                lambda imol, text: (export_map(imol, text) == 1 and
                                    add_status_bar_text("Map " + str(imol) + " exported to " + text))))
                                    

     add_simple_coot_menu_menuitem(menu,"Brighten Maps",
               lambda func: brighten_maps())


     def set_diff_map_func(imol):
	print "setting map number %s to be a difference map" %imol
        set_map_is_difference_map(imol)
        
     add_simple_coot_menu_menuitem(menu,"Set map is a difference map...",
	lambda func: map_molecule_chooser_gui("Which map should be considered a difference map?",
		lambda imol: set_diff_map_func(imol)))


     add_simple_coot_menu_menuitem(menu, "Another (contour) level...",
	lambda func: another_level())


     add_simple_coot_menu_menuitem(menu, "Multi-chicken...",
        lambda func: map_molecule_chooser_gui("Choose a molecule for multiple contouring",
		lambda imol: (set_map_displayed(imol, 0), multi_chicken(imol))))

     add_coot_menu_separator(menu)
     
     #---------------------------------------------------------------------
     #     Molecule functions
     #
     #---------------------------------------------------------------------

     add_simple_coot_menu_menuitem(menu,"Copy Fragment...", 
	lambda func: generic_chooser_and_entry("Create a new Molecule\n \
                                  From which molecule shall we seed?", 
                                 "Atom selection for fragment", "//A/1-10", 
		lambda imol, text: new_molecule_by_atom_selection(imol,text)))


     add_simple_coot_menu_menuitem(menu,"Copy Coordinates Molecule...", 
	lambda func: molecule_chooser_gui("Molecule to Copy...", 
		lambda imol: copy_molecule(imol)))


     add_simple_coot_menu_menuitem(menu,"Residue Type Selection",
	lambda func: generic_chooser_and_entry("Choose a molecule to select residues from: ","Residue Type:","",
                                               lambda imol, text: (new_molecule_by_residue_type_selection(imol, text),
                                                                   update_go_to_atom_window_on_new_mol())))


     def new_mol_sphere_func1(imol, text):
       try:
         radius = float(text)
       except:
         print "WARNING:: no valid radius", text
       args = [imol] + rotation_centre() + [radius]
       new_molecule_by_sphere_selection(*args)

     add_simple_coot_menu_menuitem(menu, "New Molecule by Sphere...",
                                   lambda func: generic_chooser_and_entry(
       "Choose a molecule from which to select a sphere of atoms:",
       "Radius:", "10.0",
       lambda imol, text: new_mol_sphere_func1(imol, text)))


# BL says:: may work, not sure about function entirely
     add_simple_coot_menu_menuitem(menu, "Replace Fragment",
	lambda func: molecule_chooser_gui("Define the molecule that needs updating",
		lambda imol_base: generic_chooser_and_entry(
				"Molecule that contains the new fragment:",
				"Atom Selection","//",
				lambda imol_fragment, atom_selection_str:
				replace_fragment(imol_base, imol_fragment, atom_selection_str))))


     add_simple_coot_menu_menuitem(menu, "Renumber Waters...",
	lambda func: molecule_chooser_gui("Renumber waters of which molecule?",
		lambda imol: renumber_waters(imol)))


     def apply_ncs_loop_func(imol, chain_id, text1, text2, dum):
       try:
         r1 = int(text1)
         r2 = int(text2)
         copy_residue_range_from_ncs_master_to_others(imol, chain_id, r1, r2)
       except:
         print "BL WARNING:: no valid number input"
       
     add_simple_coot_menu_menuitem(menu, "Apply NCS loop...",
                                   lambda func: generic_chooser_and_entry("Apply NCS Range from Master",
                                               "Master Chain ID",
                                               "",
                lambda imol, chain_id: generic_double_entry("Start of Residue Number Range",
                                       "End of Residue Number Range",
                                       "", "", False, False,
                                       lambda text1, text2, dum: apply_ncs_loop_func(imol, chain_id, text1, text2, dum))))


     add_simple_coot_menu_menuitem(menu, "Residues with Alt Confs...",
	lambda func: molecule_chooser_gui(
		"Which molecule to check for Alt Confs?",
		lambda imol: alt_confs_gui(imol)))


     add_simple_coot_menu_menuitem(menu, "Cis Peptides...",
        lambda func: molecule_chooser_gui("Choose a molecule for checking for Cis Peptides",
                lambda imol: cis_peptides_gui(imol)))


     add_simple_coot_menu_menuitem(menu, "Phosphorylate this residue",
	lambda func: phosphorylate_active_residue())


     add_simple_coot_menu_menuitem(menu,"Use SEGIDs...",
	lambda func: molecule_chooser_gui("Exchange the Chain IDs, replace with SEG IDs",
		lambda imol: exchange_chain_ids_for_seg_ids(imol)))


     add_coot_menu_separator(menu)

     #---------------------------------------------------------------------
     #     Other Programs (?)
     #
     #---------------------------------------------------------------------

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
		os.spawnv(os.P_NOWAIT, ccp4mg_file_exe, args)
	else:
		print "BL WARNING:: sorry cannot find %s in $PATH" %ccp4mg_exe

     add_simple_coot_menu_menuitem(menu, "CCP4MG...",
	lambda func: ccp4mg_func1())


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
       options_menu_mol_list_pair = generic_molecule_chooser(vbox, chooser_hint_text)
       entry = file_selector_entry(vbox, entry_hint_text)

       def shelx_delete_event(*args):
         window.destroy()
         return False

       def shelx_go_funcn_event(*args):
         import operator
         txt = entry.get_text()
         imol = get_option_menu_active_molecule(option_menu_list_pair)
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
       
     add_simple_coot_menu_menuitem(menu, "SHELXL Refine...", 
                                   lambda func: shelx_ref_func())


     add_coot_menu_separator(menu)

     # ---------------------------------------------------------------------
     #     Building
     # ---------------------------------------------------------------------

     add_simple_coot_menu_menuitem(menu, "Dock Sequence...", 
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
     
     add_simple_coot_menu_menuitem(menu, "Associate Sequence...",
                                   lambda func: generic_chooser_entry_and_file_selector(
       "Associate Sequence to Model: ",
       valid_model_molecule_qm,
       "Chain ID",
       "",
       "Select PIR file",
       lambda imol, chain_id, seq_file_name: associate_seq_func(imol, chain_id, seq_file_name)))

     add_simple_coot_menu_menuitem(menu,"Add Strand Here...",
                                   lambda func: place_strand_here_gui())


     add_coot_menu_separator(menu)

     # ---------------------------------------------------------------------
     #     Settings
     # ---------------------------------------------------------------------
     #


     add_simple_coot_menu_menuitem(menu, "Set Refinement options...",
                                   lambda func: refinement_options_gui())
 
     submenu = gtk.Menu()
     menuitem2 = gtk.MenuItem("Peptide Restraints...")

     menuitem2.set_submenu(submenu)
     menu.append(menuitem2)
     menuitem2.show()

     def add_restr_func1():
	print 'Planar Peptide Restraints added'
	add_planar_peptide_restraints()

     add_simple_coot_menu_menuitem(submenu, "Add Planar Peptide Restraints",
			lambda func: add_restr_func1())


     def add_restr_func2():
	print 'Planar Peptide Restraints removed'
	remove_planar_peptide_restraints()

     add_simple_coot_menu_menuitem(submenu, "Remove Planar Peptide Restraints",
	lambda func: add_restr_func2())


     # An example with a submenu:
     #
     submenu = gtk.Menu()
     menuitem2 = gtk.MenuItem("Refinement Speed...")
     menuitem2.set_submenu(submenu)
     menu.append(menuitem2)
     menuitem2.show()

     add_simple_coot_menu_menuitem(submenu, "Molasses Refinement mode", 
	lambda func: (printf("Molasses..."),
                      set_dragged_refinement_steps_per_frame(4)))


     add_simple_coot_menu_menuitem(submenu, "Crocodile Refinement mode", 
	lambda func: (printf("Crock..."),
                      set_dragged_refinement_steps_per_frame(120)))


     add_simple_coot_menu_menuitem(submenu, "Normal Refinement mode (1 Emsley)", 
	lambda func: (printf("Default Refinement mode (1 Emsley)..."),
                      set_dragged_refinement_steps_per_frame(50)))


# BL says:: maybe check if number at some point
     add_simple_coot_menu_menuitem(menu, "Set Spin Speed",
		lambda func: generic_single_entry("Set Spin Speed (smaller is slower)",
			str(idle_function_rotate_angle()), "Set it",
			lambda text: set_idle_function_rotate_angle(float(text))))


     add_simple_coot_menu_menuitem(menu, "Nudge Centre...",
	lambda func: nudge_screen_centre_gui())


     add_simple_coot_menu_menuitem(menu,"Set Undo Molecule...",
	lambda func: molecule_chooser_gui("Set the Molecule for 'Undo' Operations",
		lambda imol: set_undo_molecule(imol)))


# BL says: has no checking for text = number yet
     add_simple_coot_menu_menuitem(menu, "B factor bonds scale factor...",
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

     add_simple_coot_menu_menuitem(menu,"Set Matrix (Refinement Weight)...",
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

     add_simple_coot_menu_menuitem(menu,"Set Density Fit Graph Weight...",
	lambda func: generic_single_entry("set weight (smaller means apparently better fit)", 
		"1.0", "Set it", 
		lambda text: set_den_gra_func(text)))


#
     submenu = gtk.Menu()
     menuitem2 = gtk.MenuItem("Rotate Translate Zone Mode...")
 
     menuitem2.set_submenu(submenu)
     menu.append(menuitem2)
     menuitem2.show()

     add_simple_coot_menu_menuitem(submenu, "Rotate About Fragment Centre",
	lambda func: set_rotate_translate_zone_rotates_about_zone_centre(1))


     add_simple_coot_menu_menuitem(submenu, "Rotate About Second Clicked Atom",
	lambda func: set_rotate_translate_zone_rotates_about_zone_centre(0))


     def all_mol_symm_func():
	for imol in molecule_number_list():
		if valid_model_molecule_qm(imol):
			set_symmetry_whole_chain(imol, 1)

     add_simple_coot_menu_menuitem(menu, "All Molecules use \"Near Chains\" Symmetry", 
	lambda func: map(lambda imol: valid_model_molecule_qm(imol) and
                         set_symmetry_whole_chain(imol, 1),
                         molecule_number_list()))
#	lambda func: all_mol_symm_func())


     add_simple_coot_menu_menuitem(menu, "Question Accept Refinement", 
	lambda func: set_refinement_immediate_replacement(0))


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

     add_simple_coot_menu_menuitem(menu, "Save Dialog Positions...",
                                   lambda func: save_dialog_func())


     # ---------------------------------------------------------------------
     #     Views/Representations
     # ---------------------------------------------------------------------
     #
     add_coot_menu_separator(menu)

     def make_ball_n_stick_func(imol, text):
       bns_handle = make_ball_and_stick(imol, text, 0.18, 0.3, 1)
       print "handle: ", bns_handle

     add_simple_coot_menu_menuitem(menu, "Ball & Stick...",
        lambda func: generic_chooser_and_entry("Ball & Stick",
                                               "Atom Selection:",
                                               "//A/1-2",
                                               lambda imol, text: make_ball_n_stick_func(imol, text)))


     add_simple_coot_menu_menuitem(menu, "Clear Ball & Stick...",
                                   lambda func: molecule_chooser_gui(
       "Choose a molecule from which to clear Ball&Stick objects",
       lambda imol: clear_ball_and_stick(imol)))


     def hilight_site_func():
       active_atom = active_residue()
       if (active_atom):
         imol = active_atom[0]
         centre_residue_spec = [active_atom[1],
                                active_atom[2],
                                active_atom[3]]
         hilight_binding_site(imol, centre_residue_spec, 230,4)

     add_simple_coot_menu_menuitem(menu, "Highlight Interesting Site (here)...",
                                   lambda func: hilight_site_func())
     

     def make_dot_surf_func(imol,text):
	dots_handle = dots(imol, text, 1, 1)
	print "dots handle: ", dots_handle

     add_simple_coot_menu_menuitem(menu,"Dotted Surface...",
	lambda func: generic_chooser_and_entry("Surface for molecule", 
		"Atom Selection:", "//A/1-2", 
		lambda imol, text: make_dot_surf_func(imol, text)))


     def clear_dot_surf_func(imol,text):
	try:
		n = int(text)
		clear_dots(imol,n)
	except:
		print "BL WARNING:: dots handle number shall be an integer!!"

     add_simple_coot_menu_menuitem(menu,"Clear Surface Dots...",
	lambda func: generic_chooser_and_entry("Molecule with Dotted Surface", 
		"Dots Handle Number:", "0", 
		lambda imol, text: clear_dot_surf_func(imol, text)))


     add_simple_coot_menu_menuitem(menu, "Label All CAs...",
        lambda func: molecule_chooser_gui("Choose a molecule to label",
                                          lambda imol: label_all_CAs(imol)))
 
# Views submenu
     submenu = gtk.Menu()
     menuitem2 = gtk.MenuItem("Views")
 
     menuitem2.set_submenu(submenu)
     menu.append(menuitem2)
     menuitem2.show()

     add_simple_coot_menu_menuitem(submenu, "Add View...",
                        lambda func: view_saver_gui())

# BL says:: maybe check if number at some point
     add_simple_coot_menu_menuitem(submenu, "Add a Spin View...",
                  lambda func: generic_double_entry("Number of Steps", 
                         "Number of Degrees (total)", "3600", "360", 
			 False, False, 		#check button text and callback
			 "  Add Spin  ",
                         lambda text_1, text_2: add_spin_view("Spin", int(text_1), float(text_2))))

     add_simple_coot_menu_menuitem(submenu, "Views Panel...",
                        lambda func: views_panel_gui())
 
     add_simple_coot_menu_menuitem(submenu, "Play Views",
                        lambda func: map(eval,["go_to_first_view(1)",
                        "time.sleep(1)", "play_views()"]))
 
# BL says:: maybe check if number at some point
     add_simple_coot_menu_menuitem(submenu, "Set Views Play Speed...",
		lambda func: generic_single_entry("Set Views Play Speed",
			str(views_play_speed()), "  Set it  ",
			lambda text: set_views_play_speed(float(text))))

     add_simple_coot_menu_menuitem(submenu, "Save Views...",
		lambda func: generic_single_entry("Save Views",
			"coot-views.py", " Save ",
			lambda txt: save_views(txt)))

     # add to validate menu
     menu = coot_menubar_menu("Validate")

     add_simple_coot_menu_menuitem(menu, "Pukka Puckers...?",
                                   lambda func: molecule_chooser_gui(
       "Choose a molecule for ribose pucker analysis",
       lambda imol: pukka_puckers_qm(imol)))

  else:
	print "BL WARNING:: could not find the main_menubar! Sorry, no extensions menu!"

