# translated by BL

# Copyright 2004, 2005 by The University of York
# Copyright 2005 by Bernhard Lohkamp
 
#;;; This program is free software; you can redistribute it and/or modify
#;;; it under the terms of the GNU General Public License as published by
#;;; the Free Software Foundation; either version 2 of the License, or (at
#;;; your option) any later version.
 
#;;; This program is distributed in the hope that it will be useful, but
#;;; WITHOUT ANY WARRANTY; without even the implied warranty of
#;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#;;; General Public License for more details.
 
#;;; You should have received a copy of the GNU General Public License
#;;; along with this program; if not, write to the Free Software
#;;; Foundation, Inc.,  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#;; direction is either 'forwards or 'backwards
#;; 
#;; start-resno is higher than stop-resno if we are building backwards
#;; 
#;; (fit-gap 0 "A"  23 26)   ; we'll build forwards
#;; (fit-gap 0 "A"  26 23)   ; we'll build backwards
#;;
def fit_gap(imol, chain_id, start_resno, stop_resno, sequence):

 if (valid_model_molecule(imol) == 0):
	
  print "Molecule number %(a)i is not a valid model molecule" %{"a":imol}

 else: 

  with_no_backups(imol)

#   ;; -----------------------------------------------
#   ;; Make poly ala
#   ;; -----------------------------------------------
	   
#				;(set-terminal-residue-do-rigid-body-refine 0)
#				;(set-add-terminal-residue-n-phi-psi-trials 2000)
  set_residue_selection_flash_frames_number = 0

  immediate_refinement_mode = refinement_immediate_replacement_state
  print " BL DEBUG:: start_resno, stop_resno",start_resno,stop_resno
  if stop_resno < start_resno:
     direction = "backwards"
  else:
     direction = "forwards"

  print "direction is ", direction
	     
  set_refinement_immediate_replacement = 1
	     
#     ;; recur over residues:
#  print "BL DEBUG:: we do adding so many times =", abs(start_resno-stop_resno)+1
  if direction == "forwards":
        resno = start_resno - 1
  else:
        resno = start_resno + 1
  for i in range(abs(start_resno - stop_resno) + 1):
	       
     print "add-terminal-residue: residue number: ",resno
     status = add_terminal_residue(imol,chain_id,resno,"ALA",1)
     if status:
#					; first do a refinement of what we have 
        refine_auto_range(imol,chain_id,resno,"")
        accept_regularizement()
        if direction == "forwards":
           resno = resno + 1
        else:
           resno = resno - 1
     else:
        print "Failure in fit-gap at residue ",resno

#     ;; -----------------------------------------------
#     ;; From poly ala to sequence (if given):
#     ;; -----------------------------------------------

  if (not sequence == ""):
	  print "mutate-and-autofit-residue-range ",imol, chain_id, start_resno,stop_resno, sequence
          mutate_and_autofit_residue_range(imol,chain_id,
                                          start_resno,stop_resno,sequence)

#     ;; -----------------------------------------------
#     ;; Refine new zone
#     ;; -----------------------------------------------

  if residue_exists(imol,chain_id,start_resno - 1,""):
         print "Test finds"
  else:
         print "Test: not there"

  if residue_exists(imol,chain_id,start_resno - 1,""):
         low_end = start_resno - 1
  else:
         low_end = start_resno
  if residue_exists(imol,chain_id,stop_resno + 1,""):
         high_end = stop_resno + 1
  else:
         high_end = stop_resno
  if direction == "forwards":
         final_zone = [low_end,high_end]
  else:
         final_zone = [high_end,low_end]

#       ;; we also need to check that start-resno-1 exists and
#       ;; stop-resno+1 exists.

  refine_zone(imol,chain_id,final_zone[0],final_zone[1],"")
#				; set the refinement dialog flag back to what it was:
  if immediate_refinement_mode == 0:
        set_refinement_immediate_replacement = 0
  accept_regularizement()

       

#;;; For Kay Diederichs, autofit without a map (best clash score)
#;;; This ignores alt conformations and residues with insertion codes
#;;;
def de_clash (imol,chain_id,resno_start,resno_end):

    resno = resno_start
    while resno <= resno_end:
	   auto_fit_best_rotamer(resno,"","",chain_id,imol,-1,1,0.1)
           resno = resno + 1


