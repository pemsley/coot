# Copyright 2006, 2007 by Bernhard Lohkamp
# Copyright 2006 by Paul Emsley, The University of York

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc.,  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import pygtk, gtk, pango
import time

# BL says: Havent been able to get the proper extensions to run,
# so for now we just make a new window with a menubar
# may be fixed at some point....
# try to program it, so that we can easily switch to proper extension when
# eventualy functional!?

try: 
  import coot_python
  if coot_python.main_menubar():

     menu = coot_menubar_menu("Extensions")

     add_simple_coot_menu_menuitem(menu,"[Post MR] Fill Partial Residues...",
	lambda func: molecule_chooser_gui("Find and Fill residues with missing atoms",
		lambda imol: fill_partial_residues(imol)))

     def fit_prot_func(imol):
	if imol_refinement_map()==-1:
	   add_status_bar_text("oops. Must set a map to fit")
	else:
	   fit_protein(imol)
     add_simple_coot_menu_menuitem(menu,"[Post MR] Fit Protein...",
	lambda func: molecule_chooser_gui("Fit Protein using Rotamer Search",
		lambda imol: fit_prot_func(imol)))

     def step_ref_func(imol):
	if imol_refinement_map()==-1:
		add_status_bar_text("oops. Must set a map to fit")
	else:
		stepped_refine_protein(imol)
     add_simple_coot_menu_menuitem(menu,"[Post MR] Stepped Refine...",
	lambda func: molecule_chooser_gui("Stepped Refine: ",
		lambda imol: step_ref_func(imol)))

     add_simple_coot_menu_menuitem(menu,"Residue Type Selection",
	lambda func: generic_chooser_and_entry("Choose a molecule to select residues from: ","Residue Type:","",
		lambda imol,text: map(eval,("new_molecule_by_residue_type_selection(imol,text)", "update_go_to_atom_window_on_new_mol()"))))


# BL says: currently doesnt check for valid molecules etc.!!
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
	return f
     def mask_map_func1(active_state):
	print "changed active_state to ", active_state
     def mask_map_func2(text_1, text_2, imol, invert_mask_qm):
	import operator
	n = int(text_1)
	if invert_mask_qm:
		invert = 1
	else:
		invert = 0
	if (operator.isNumberType(n)):
		mask_map_by_atom_selection(n, imol, text_2, invert)

     add_simple_coot_menu_menuitem(menu,"Mask Map by Atom Selection...", 
             lambda func: molecule_chooser_gui("Define the molecule that has atoms to mask the map", 
                  lambda imol: generic_double_entry("Map molecule number: ", 
			"Atom selection: ", mask_map_func(), "//A/1", 
			" Invert Masking? ", lambda active_state:
			mask_map_func1(active_state),
			"Mask Map", 
			lambda text_1, text_2, invert: mask_map_func2(text_1, text_2, imol, invert))))

     add_simple_coot_menu_menuitem(menu,"Copy Coordinates Molecule...", 
	lambda func: molecule_chooser_gui("Molecule to Copy...", 
		lambda imol: copy_molecule(imol)))

     add_simple_coot_menu_menuitem(menu,"Copy Map Molecule...", 
	lambda func: map_molecule_chooser_gui("Molecule to Copy...", 
		lambda imol: copy_molecule(imol)))

     add_simple_coot_menu_menuitem(menu,"Copy Fragment...", 
	lambda func: generic_chooser_and_entry("Create a new Molecule\nFrom which molecule shall we seed?", 
	"Atom selection for fragment", "//A/1-10", 
		lambda imol, text: new_molecule_by_atom_selection(imol,text)))

# BL says:: may work, not sure about function entirely
     add_simple_coot_menu_menuitem(menu, "Replace Fragment",
	lambda func: molecule_chooser_gui("Define the molecule that needs updating",
		lambda imol_base: generic_chooser_and_entry(
				"Molecule that contains the new fragment:",
				"Atom Selection","//",
				lambda imol_fragment, atom_selection_str:
				replace_fragment(imol_base, imol_fragment, atom_selection_str))))

#
#
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
#

     add_simple_coot_menu_menuitem(menu,"Set Undo Molecule...",
	lambda func: molecule_chooser_gui("Set the Molecule for 'Undo' Operations",
		lambda imol: set_undo_molecule(imol)))

     def set_diff_map_func(imol):
	print "setting map number %s to be a difference map" %imol
	set_map_is_difference_map(imol)
     add_simple_coot_menu_menuitem(menu,"Set map is a difference map...",
	lambda func: map_molecule_chooser_gui("Which map should be considered a difference map?",
		lambda imol: set_diff_map_func(imol)))

     add_simple_coot_menu_menuitem(menu,"Use SEGIDs...",
	lambda func: molecule_chooser_gui("Exchange the Chain IDs, replace with SEG IDs",
		lambda imol: exchange_chain_ids_for_seg_ids(imol)))

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

# BL says: has no checking for text = number yet
     add_simple_coot_menu_menuitem(menu, "B factor bonds scale factor...",
	lambda func: generic_chooser_and_entry("Choose a molecule to which the B-factor colour scale is applied:",
		"B-factor scale:", "1.0", 
		lambda imol, text: set_b_factor_bonds_scale_factor(imol,float(text))))

     add_simple_coot_menu_menuitem(menu,"Renumber Waters...",
	lambda func: molecule_chooser_gui("Renumber waters of which molecule?",
		lambda imol: renumber_waters(imol)))

     add_simple_coot_menu_menuitem(menu, "Residues with Alt Confs...",
	lambda func: molecule_chooser_gui(
		"Which molecule to check for Alt Confs?",
		lambda imol: alt_confs_gui(imol)))

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

# An example with a submenu:
#
     submenu = gtk.Menu()
     menuitem2 = gtk.MenuItem("Refinement Speed...")
     menuitem2.set_submenu(submenu)
     menu.append(menuitem2)
     menuitem2.show()

     def ref_mode_func1():
	print "Molasses..."
	set_dragged_refinement_steps_per_frame(4)
     add_simple_coot_menu_menuitem(submenu, "Molasses Refinement mode", 
	lambda func: ref_mode_func1())

     def ref_mode_func2():
	print "Crock..."
	set_dragged_refinement_steps_per_frame(120)
     add_simple_coot_menu_menuitem(submenu, "Crocodile Refinement mode", 
	lambda func: ref_mode_func2())

     def ref_mode_func2():
	print "Back to normal (1 Emsley)..."
	set_dragged_refinement_steps_per_frame(50)
     add_simple_coot_menu_menuitem(submenu, "Normal speed Refinement mode (1 Emsley)", 
	lambda func: ref_mode_func2())

#
# another submenu
     submenu = gtk.Menu()
     menuitem2 = gtk.MenuItem("Rotate Translate Zone Mode...")
 
     menuitem2.set_submenu(submenu)
     menu.append(menuitem2)
     menuitem2.show()

     add_simple_coot_menu_menuitem(submenu, "Rotate About Fragment Centre",
	lambda func: set_rotate_translate_zone_rotates_about_zone_centre(1))

     add_simple_coot_menu_menuitem(submenu, "Rotate About Second Clicked Atom",
	lambda func: set_rotate_translate_zone_rotates_about_zone_centre(0))

#
# back to 'normal' extensions
     def all_mol_symm_func():
	for imol in molecule_number_list():
		if valid_model_molecule_qm(imol):
			set_symmetry_whole_chain(imol, 1)
     add_simple_coot_menu_menuitem(menu, "All Molecules use \"Near Chains\" Symmetry", 
	lambda func: all_mol_symm_func())

     add_simple_coot_menu_menuitem(menu, "Question Accept Refinement", 
	lambda func: set_refinement_immediate_replacement(0))

     add_simple_coot_menu_menuitem(menu, "Another (contour) level...",
	lambda func: another_level())

     def ccp4mg_func1():
	import os
	pd_file_name = "1.mgpic.py"
	write_ccp4mg_picture_description(pd_file_name)
	if os.name == 'nt':
		ccp4mg_exe = "winccp4mg.exe"
	else:
		ccp4mg_exe = "ccp4mg"
	if command_in_path_qm(ccp4mg_exe):
		ccp4mg_file_exe = find_exe(ccp4mg_exe, "PATH", "")
		pd_file_name = os.path.abspath(pd_file_name)
		args = [ccp4mg_file_exe, "-pict", pd_file_name]
		os.spawnv(os.P_NOWAIT, ccp4mg_file_exe, args)
	else:
		print "BL WARNING:: sorry cannot find %s in $PATH" %ccp4mg_exe
     add_simple_coot_menu_menuitem(menu, "CCP4MG...",
	lambda func: ccp4mg_func1())

     def shelx_ref_func(imol,text):
	if len(text)>0:
		shelxl_refine(imol,text)
	else:
		shelxl_refine(imol)
     add_simple_coot_menu_menuitem(menu, "SHELXL Refine...", 
	lambda func: generic_chooser_and_entry("Molecule for refinement:", 
		"HKL data filename (leave blank for default)", "", 
		lambda imol, text: shelx_ref_func(imol,text)))

     add_simple_coot_menu_menuitem(menu, "Phosphorylate this residue",
	lambda func: phosphorylate_active_residue())

# BL says:: maybe check if number at some point
     add_simple_coot_menu_menuitem(menu, "Set Spin Speed",
		lambda func: generic_single_entry("Set Spin Speed (smaller is slower)",
			str(idle_function_rotate_angle()), "Set it",
			lambda text: set_idle_function_rotate_angle(float(text))))

     add_simple_coot_menu_menuitem(menu,"Brighten Maps",
               lambda func: brighten_maps())

     add_simple_coot_menu_menuitem(menu, "Nudge Centre...",
	lambda func: nudge_screen_centre_gui())

# Not this.  Instead fire up a new top level, where we have a molecule chooser, 
# an entry for the chain spec and a text window where can paste in a sequence.
# 
# No, that doesn't work either.  We need a set of pairs of entries
# and text boxes.
# 
# Hmm..

#      (add-simple-coot-menu-menuitem
#       menu "Cootaneer this fragment [try sequence assignment]"
#       (lambda ()
#	 (let ((imol-map (imol-refinement-map)))
#	   (if (= imol-map -1)
#	       (info-dialog "Need to assign a map to fit against")
#	       (let ((active-atom (active-residue)))
#		 (if (list? active-atom)
#		     (let ((imol     (list-ref active-atom 0))
#			   (chain-id (list-ref active-atom 1))
#			   (resno    (list-ref active-atom 2))
#			   (inscode  (list-ref active-atom 3))
#			   (at-name  (list-ref active-atom 4))
#			   (alt-conf (list-ref active-atom 5)))
#		       (cootaneer imol-map imol (list chain-id resno inscode 
#						      at-name alt-conf)))))))))

#     submenu = gtk.Menu()
#     menuitem2 = gtk.MenuItem("Cootaneering")

#     menuitem2.set_submenu(submenu)
#     menu.append(menuitem2)
#     menuitem2.show()

#     def cootaneer_func1(imol, chain_id, pir_file_name):
#	print "assoc seq: ", imol, chain_id, pir_file_name
#	import os
#	if os.path.isfile(pir_file_name):
#		fin = open(pir_file_name, 'r')
#		lines = fin.readlines()
#		for line in lines:
#			seq_text += line
#		fin.close()
#		assign_pir_sequence(imol, chain_id, seq_text)
#	else:
#		print "BL WARINING: no such file ", pir_file_name
#
#     add_simple_coot_menu_menuitem(submenu, "Associate Sequence....", 
#	lambda func: generic_chooser_entry_and_file_selector(
#		"Associate Sequence to Model: ",
#		"Chain ID", "", "Select PIR file", 
#		lambda imol, chain_id, pir_file_name: 
#			cootaneer_func1(imol, chain_id, pir_file_name)))

#     add_simple_coot_menu_menuitem(submenu, "Cootaneer this fragment...",
#	lambda func: molecule_chooser_gui("Choose a molecule to apply sequence assignment", 
#		lambda imol: coontaneer_gui(imol)))


     add_simple_coot_menu_menuitem(menu,"Add Strand Here...",
                        lambda func: place_strand_here_gui())

 
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

  else:
	print "BL WARNING:: could not find the main_menubar! Sorry, no extensions!"

except: print """BL WARNING:: could not import coot_python module!!
Some things, esp. extensions, may be crippled!"""

# finally, show it

