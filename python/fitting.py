# translated to by Bernhard Lohkamp

#;;;; Copyright 2004 The University of York
#;;;; Copyright 2004, 2005, 2006, 2007 Bernhard Lohkamp
 
#;;;; This program is free software; you can redistribute it and/or modify
#;;;; it under the terms of the GNU General Public License as published by
#;;;; the Free Software Foundation; either version 2 of the License, or (at
#;;;; your option) any later version.
 
#;;;; This program is distributed in the hope that it will be useful, but
#;;;; WITHOUT ANY WARRANTY; without even the implied warranty of
#;;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#;;;; General Public License for more details.
 
#;;;; You should have received a copy of the GNU General Public License
#;;;; along with this program; if not, write to the Free Software
#;;;; Foundation, Inc.,  59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


#;;; 
def fit_protein(imol):

# for testing we set chain_ids
#    chain_ids=["A","B"]

    set_go_to_atom_molecule(imol)
    make_backup(imol) # do a backup first
    backup_mode = backup_state(imol)
    imol_map  = imol_refinement_map()
    alt_conf = ""
    replacement_state = refinement_immediate_replacement_state()

    turn_off_backup(imol)
    set_refinement_immediate_replacement(1)
	  
    for chain_id in chain_ids(imol):
	     n_residues = chain_n_residues(chain_id,imol)
             print "There are %(a)i residues in chain %(b)s" % {"a":n_residues,"b":chain_id}
	       
	     for serial_number in range(n_residues):
		  
                res_name = resname_from_serial_number(imol,chain_id,serial_number)
                res_no = seqnum_from_serial_number(imol,chain_id,serial_number)
                ins_code = insertion_code_from_serial_number(imol,chain_id,serial_number)
                print "centering on ",chain_id,res_no,"CA"
                set_go_to_atom_chain_residue_atom_name(chain_id,res_no,"CA")
                rotate_y_scene(30,0.3) # n-frames frame-interval(degrees)
#                print" BL DEBUG:: params :",res_no,alt_conf,ins_code,chain_id,imol,imol_map,1,0.1
                auto_fit_best_rotamer(res_no,alt_conf,ins_code,chain_id,imol,imol_map,1,0.1)
                if (imol_map >= 0):
#        	  refine_auto_range(imol,chain_id,res_no,alt_conf)
                  if (res_no != n_residues):
 		      refine_zone(imol,chain_id,res_no-1,res_no+1,alt_conf)
                  else:
 		      refine_zone(imol,chain_id,res_no-1,res_no,alt_conf)
		  accept_regularizement()
                rotate_y_scene(30,0.3)
      
    if (replacement_state == 0):
	  set_refinement_immediate_replacement(0)
    if (backup_mode == 1):
	  turn_on_backup(imol)

def fit_waters(imol):
    imol_map = imol_refinement_map()
#    print "BL DEBUG:: imol_map is:", imol_map
    if (imol_map != -1):
       replacement_state = refinement_immediate_replacement_state()
       backup_mode = backup_state(imol)
       alt_conf = ""

       print "BL DEBUG:: replacement_state is :",replacement_state

       turn_off_backup(imol)
       set_refinement_immediate_replacement(1)
       
       for chain_id in chain_ids(imol):
#            print "BL DEBUG:: is_solven_chain:",is_solvent_chain_qm(imol,chain_id)
            if (is_solvent_chain_qm(imol,chain_id)):
                print"BL DEBUG:: we found a solvent chain: ",chain_id
                n_residues = chain_n_residues(chain_id,imol)
                print "There are %(a)i residues in chain %(b)s" % {"a":n_residues,"b":chain_id}
                for serial_number in range(n_residues):
                     res_no = seqnum_from_serial_number(imol,chain_id,serial_number)
                     print "BL DEBUG:: imol, chain_id, res_no, alt_conf",imol,chain_id,res_no,alt_conf
                     refine_zone(imol,chain_id,res_no,res_no,alt_conf)
                     accept_regularizement()
 
       if (replacement_state == 0):
          set_refinement_immediate_replacement(0)
       if (backup_mode == 1):
          turn_on_backup(imol)
