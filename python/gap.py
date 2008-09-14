# gap.py 
# Copyright 2004, 2005, 2006 by Paul Emsley, The University of York
# Copyright 2005, 2006 Bernhard Lohkamp
# Copyright 2008 by Bernhard Lohkamp, The University of York
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


# Fit missing loop in protein
# fit both directions if end residues exist and start gui to select between
# both models if they are different
# if only one end, do 'normal' loop fit
# direction can be forwards or backwards
#

def fit_gap(imol, chain_id, start_resno, stop_resno, sequence=""):
	
	imol_map = imol_refinement_map()
	if (imol_map == -1):
		info_dialog("Need to set a map to fit a loop")
	else:
		# normal execution
		if (stop_resno < start_resno):
			res_limits = [stop_resno - 1, start_resno + 1]
		else:
			res_limits = [start_resno - 1, stop_resno + 1]

		if all(map(lambda resno: residue_exists_qm(imol, chain_id, resno,""), res_limits)):
			# build both ways
			imol_copy = copy_molecule(imol)
			# build A:
			fit_gap_generic(imol, chain_id, start_resno, stop_resno, sequence)
			# build B:
			fit_gap_generic(imol_copy, chain_id, stop_resno, start_resno, sequence)

			# get the fit result:
			result_a = low_density_average(imol_map, imol, chain_id, start_resno, stop_resno)
			result_b = low_density_average(imol_map, imol_copy, chain_id, start_resno, stop_resno)

			# print "BL DEBUG:: fit a, b:", result_a, result_b
			cut_off = 0.90
			if (min(result_a,  result_b) / max(result_a, result_b) >= cut_off):
				# solutions are identical?! (cut-off 99%)
				# we keep only the first one
				close_molecule(imol_copy)
			else:
				# different by cut-off -> display both options if pygtk
				# otherwise use the 'better' one
				if ('pygtk' in sys.modules.keys()):
					# have pygtk
					atom_selection = "//" + chain_id + "/" + str(min(start_resno, stop_resno) - 1) + \
									 "-"  + str(max(start_resno, stop_resno) + 1)
					# we make a fragment for each loop
					imol_fragment_a = new_molecule_by_atom_selection(imol     , atom_selection)
					imol_fragment_b = new_molecule_by_atom_selection(imol_copy, atom_selection)
					close_molecule(imol_copy)															
					
					buttons = [[" Loop A ", "(replace_fragment(" + str(imol) + ", " + str(imol_fragment_a) + ", \""\
								+ atom_selection + "\"), set_mol_displayed(" + str(imol_fragment_a) + \
								", 0), set_mol_displayed(" + str(imol_fragment_b) + ", 0))"],
							   [" Loop B ", "(replace_fragment(" + str(imol) + ", " + str(imol_fragment_b) + ", \""\
								+ atom_selection + "\"), set_mol_displayed(" + str(imol_fragment_a) + \
								", 0), set_mol_displayed(" + str(imol_fragment_b) + ", 0))"]]
					selected_button = 0
					if (result_b > result_a):
						selected_button = 1

					go_function = "(close_molecule(" + str(imol_fragment_a) + "), close_molecule(" + str(imol_fragment_b) + "))"
					dialog_box_of_radiobuttons("Select Loop", [200, 100],
											   buttons, "  Accept  ",
											   go_function, selected_button)
										
				else:
					# no pygtk take better one
					if (result_a > result_b):
						close_molecule(imol_copy)
					else:
						atom_selection = "//" + str(min(start_resno, stop_resno) - 1) + \
										 "-"  + str(max(start_resno, stop_resno) + 1)
						replace_fragment(imol, imol_copy, atom_selection)
						close_molecule(imol_copy)
						

		else:
			# either end residue is missing -> single build
			fit_gap_generic(imol, chain_id, start_resno, stop_resno, sequence)
			

# Fit missing loop in protein.
#
# direction is either 'forwards' or 'backwards'
#
# start-resno is higher than stop-resno if we are building backwards
# 
# fit_gap_generic(0,"A",23,26)   ; we'll build forwards
#
# fit_gap_generic(0,"A",26,23)   ; we'll build backwards
#
def fit_gap_generic(imol, chain_id, start_resno, stop_resno, sequence=""):

 import string
 sequence = string.upper(sequence)

 if (valid_model_molecule_qm(imol) == 0):
  print "Molecule number %(a)i is not a valid model molecule" %{"a":imol}
 else: 
  backup_mode = backup_state(imol)
  turn_off_backup(imol)

  # -----------------------------------------------
  # Make poly ala
  # -----------------------------------------------
  set_residue_selection_flash_frames_number(0)

  immediate_refinement_mode = refinement_immediate_replacement_state()
#  print " BL DEBUG:: start_resno, stop_resno",start_resno,stop_resno
  if (stop_resno < start_resno):
     direction = "backwards"
  else:
     direction = "forwards"

  print "direction is ", direction
	     
  set_refinement_immediate_replacement(1)
	     
  # recur over residues:
  if direction == "forwards":
        resno = start_resno - 1
  else:
        resno = start_resno + 1
  for i in range(abs(start_resno - stop_resno) + 1):
	       
     print "add-terminal-residue: residue number: ",resno
     status = add_terminal_residue(imol, chain_id, resno, "ALA", 1)
     if status:
        # first do a refinement of what we have 
        refine_auto_range(imol, chain_id, resno, "")
        accept_regularizement()
        if direction == "forwards":
           resno = resno + 1
        else:
           resno = resno - 1
     else:
        print "Failure in fit-gap at residue ",resno

  # -----------------------------------------------
  # From poly ala to sequence (if given):
  # -----------------------------------------------

  if (not sequence == ""):
	  print "mutate-and-autofit-residue-range ",imol, chain_id, start_resno,stop_resno, sequence
          if direction == "forwards":
              mutate_and_autofit_residue_range(imol, chain_id,
                                               start_resno, stop_resno,
                                               sequence)
          else:
              mutate_and_autofit_residue_range(imol, chain_id,
                                               stop_resno, start_resno,
                                               sequence)

  # -----------------------------------------------
  # Refine new zone
  # -----------------------------------------------

  if residue_exists_qm(imol,chain_id,start_resno - 1,""):
         print "Test finds"
  else:
         print "Test: not there"

  if residue_exists_qm(imol,chain_id,start_resno - 1,""):
         low_end = start_resno - 1
  else:
         low_end = start_resno
  if residue_exists_qm(imol,chain_id,stop_resno + 1,""):
         high_end = stop_resno + 1
  else:
         high_end = stop_resno
  if direction == "forwards":
         final_zone = [low_end,high_end]
  else:
         final_zone = [high_end,low_end]

 # we also need to check that start-resno-1 exists and
 # stop-resno+1 exists.

  refine_zone(imol,chain_id,final_zone[0],final_zone[1],"")
	# set the refinement dialog flag back to what it was:
  if immediate_refinement_mode == 0:
        set_refinement_immediate_replacement(0)
  accept_regularizement()

  if (backup_mode==1):
     turn_on_backup(imol)


# For Kay Diederichs, autofit without a map (find rotamer with best
# clash score). This ignores alt conformations and residues with
# insertion codes.
#
def de_clash (imol,chain_id,resno_start,resno_end):

    resno = resno_start
    while resno <= resno_end:
	   auto_fit_best_rotamer(resno,"","",chain_id,imol,-1,1,0.1)
           resno = resno + 1

# calculate the average of 20% lowest density at all atom_positions in fragment
def low_density_average(imol_map, imol, chain_id, start_resno, stop_resno):

	map_coords = []
	map_density = []
	for resno in range(start_resno, stop_resno + 1):
		atom_ls = residue_info(imol, chain_id, resno, "")  # ignoring ins_code
		for atom in atom_ls:
			# we only take main chain + CB into account
			if atom[0][0] in [' N  ', ' CA ', ' CB ', ' C  ', ' O  ']:  
				map_coords.append(atom[2])

	for [x, y, z] in map_coords:
		map_density.append(density_at_point(imol_map, x, y,  z))

	map_density.sort()  # sort ascending

	cut_off = len(map_density) / 5
	if (cut_off <= 0):
		cut_off = 1   # take at least one point

	map_average = sum(map_density[0:cut_off]) / float(cut_off)  # make it a float
	return map_average
