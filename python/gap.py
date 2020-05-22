# gap.py
# Copyright 2004, 2005, 2006 by Paul Emsley, The University of York
# Copyright 2005, 2006 Bernhard Lohkamp
# Copyright 2008 by Bernhard Lohkamp, The University of York
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# Fit missing loop in protein
# fit both directions if end residues exist and start gui to select between
# both models if they are different
# if only one end, do 'normal' loop fit
# direction can be forwards or backwards
#

def fit_gap(imol, chain_id, start_resno, stop_resno,
            sequence="", use_rama_restraints=1):

   imol_map = imol_refinement_map()
   if (imol_map == -1):
      info_dialog("Need to set a map to fit a loop")
   else:
      # normal execution
      backup_mode = backup_state(imol)
      make_backup(imol)
      turn_off_backup(imol)

      # backup the rama_status
      rama_status = refine_ramachandran_angles_state()
      set_refine_ramachandran_angles(use_rama_restraints)

      if (stop_resno < start_resno):
         res_limits = [stop_resno - 1, start_resno + 1]
      else:
         res_limits = [start_resno - 1, stop_resno + 1]

      if all([residue_exists_qm(imol, chain_id, resno,"") for resno in res_limits]):
         # build both ways
         imol_backwards = copy_molecule(imol)
         loop_len = abs(start_resno - stop_resno) + 1
         if (loop_len >= 6):
            imol_both = copy_molecule(imol)

         # make a backup copy of the original terminal residues
         atom_selection = "//" + chain_id + "/" + str(min(start_resno, stop_resno) - 1) + \
                          "-"  + str(max(start_resno, stop_resno) + 1)
         imol_fragment_backup = new_molecule_by_atom_selection(imol, atom_selection)
         set_mol_displayed(imol_fragment_backup, 0)

         # build A:
         fit_gap_generic(imol, chain_id, start_resno, stop_resno, sequence)
         # build B:
         fit_gap_generic(imol_backwards, chain_id, stop_resno, start_resno, sequence)

         # get the fit result:
         result_a = low_density_average(imol_map, imol, chain_id, start_resno, stop_resno)
         result_b = low_density_average(imol_map, imol_backwards, chain_id, start_resno, stop_resno)
         loop_list = [[imol, result_a], [imol_backwards, result_b]]

         # if longer loop (>=6) build half from both sides
         if (loop_len >= 6):
            start_resno1 = min([start_resno, stop_resno])
            stop_resno2  = max([start_resno, stop_resno])
            stop_resno1  = start_resno1 + loop_len/2 - 1
            start_resno2 = stop_resno1 + 1
            sequence1    = sequence[0:loop_len//2]
            sequence2    = sequence[loop_len//2:len(sequence)]
            fit_gap_generic(imol_both, chain_id, start_resno1, stop_resno1, sequence1)
            fit_gap_generic(imol_both, chain_id, stop_resno2, start_resno2, sequence2)
            # finally refine the 'gap'; check refinement immediate status
            immediate_refinement_mode = refinement_immediate_replacement_state()
            set_refinement_immediate_replacement(1)
            refine_zone(imol_both, chain_id, stop_resno1 - loop_len//3,
                                            start_resno2 + loop_len//3, "")
            accept_regularizement()
            set_refinement_immediate_replacement(immediate_refinement_mode)
            result_c = low_density_average(imol_map, imol_both, chain_id, start_resno, stop_resno)
            loop_list.append([imol_both, result_c])

         #print "BL DEBUG:: fit a, b:", result_a, result_b, result_c
         cut_off = 0.90
         # filter out redundant results based on cut-off
         i = 0
         while i < (len(loop_list) -1):
            j = i + 1
            while j < len(loop_list):
               if (min(loop_list[i][1], loop_list[j][1]) / max(loop_list[i][1], loop_list[j][1]) > cut_off):
                  # solutions are identical?! (cut-off 90%)
                  close_molecule(loop_list[j][0])
                  loop_list.pop(j)
               j += 1
            i += 1

         # different by cut-off -> display both options if pygtk
         # otherwise use the 'better' one
         if ('pygtk' in list(sys.modules.keys())):
            # have pygtk
            # we make a fragment for each loop
            fragment_list = []
            for i, (imol_loop, result) in enumerate(loop_list):
               imol_fragment = new_molecule_by_atom_selection(imol_loop, atom_selection)
               fragment_list.append([imol_fragment, result])
               set_mol_displayed(imol_fragment, 0)
               # we close all mols and work with the fragments (except the original one)
               if (i > 0):
                  close_molecule(imol_loop)

            buttons = []
            selected_button = 0
            max_result = fragment_list[0][1]
            for i, (imol_fragment, result) in enumerate(fragment_list):
               label_str   = " Loop " + chr(65 + i) + " "
               replace_str = "replace_fragment(" + str(imol) + ", " + str(imol_fragment) + \
                                                 ", \"" + atom_selection + "\"), "
               #display_ls  = map(lambda (x, y): "set_mol_displayed(" + str(x) + ", 0)", fragment_list)

               #buttons.append([label_str, "(" + replace_str + ', '.join(display_ls) + ")"])
               buttons.append([label_str, "(" + replace_str + ")"])

               # activate if better result than previous
               if (i < len(fragment_list) - 1):
                  if (max_result < fragment_list[i+1][1]):
                     selected_button = i + 1
                     max_result = fragment_list[i+1][1]

            fragment_list.append([imol_fragment_backup, -99999.])
            close_ls    = ["close_molecule(" + str(x_y[0]) + ")" for x_y in fragment_list]
            go_function = "(" + ', '.join(close_ls) + ")"
            cancel_function = "(delete_residue_range(" + str(imol) + ", \"" + str(chain_id) + \
                              "\", " + str(min(start_resno, stop_resno)) + \
                              ", " + str(max(start_resno, stop_resno)) + \
                              "), replace_fragment(" + str(imol) + ", " + \
                              str(imol_fragment_backup) + ", \"" + atom_selection + "\"), " + \
                              ', '.join(close_ls) + ")"


            # only show if more than 1 loop left
            if (len(buttons) >1):
               dialog_box_of_radiobuttons("Select Loop", [200, 100],
                                          buttons, "  Accept  ",
                                          go_function, selected_button,
                                          "  Reject  ", cancel_function)

         else:
            # no pygtk take best one
            if (result_a > result_b):
               close_molecule(imol_copy)
            else:
               replace_fragment(imol, imol_copy, atom_selection)
               close_molecule(imol_copy)


      else:
         # either end residue is missing -> single build
         fit_gap_generic(imol, chain_id, start_resno, stop_resno, sequence)

      set_refine_ramachandran_angles(rama_status)
      if (backup_mode==1):
         turn_on_backup(imol)



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
      print("Molecule number %(a)i is not a valid model molecule" %{"a":imol})
   else:

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

      print("direction is ", direction)

      set_refinement_immediate_replacement(1)

      # recur over residues:
      if direction == "forwards":
         resno = start_resno - 1
      else:
         resno = start_resno + 1

      for i in range(abs(start_resno - stop_resno) + 1):

         print("add-terminal-residue: residue number: ",resno)
         status = add_terminal_residue(imol, chain_id, resno, "auto", 1)
         if status:
            # first do a refinement of what we have
            refine_auto_range(imol, chain_id, resno, "")
            accept_regularizement()
            if direction == "forwards":
               resno = resno + 1
            else:
               resno = resno - 1
         else:
            print("Failure in fit-gap at residue ",resno)



      # -----------------------------------------------
      # From poly ala to sequence (if given):
      # -----------------------------------------------
      # only if sequence is hasnt been assigned

      if (not sequence == "" and not has_sequence_qm(imol, chain_id)):
         print("mutate-and-autofit-residue-range ",imol, chain_id, start_resno,stop_resno, sequence)
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
         print("Test finds")
      else:
         print("Test: not there")

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


# helper function to see if a sequence has been assigned to a chain in imol
# return True if sequence is there, False otherwise
#
def has_sequence_qm(imol, chain_id_ref):
   if sequence_info(imol):
      for item in sequence_info(imol):
         chain_id = item[0]
         sequence = item[1]
         if (chain_id_ref == chain_id and len(sequence) > 0):
            return True
   return False


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

	cut_off = len(map_density) // 5
	if (cut_off <= 0):
		cut_off = 1   # take at least one point

	map_average = sum(map_density[0:cut_off]) / float(cut_off)  # make it a float
	return map_average
