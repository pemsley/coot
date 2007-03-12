#print " we load mutate, maybe..."

def mutate_chain(imol,chain_id,sequence):

   if (len(sequence) != chain_n_iresidues(chain_id, imol)) :
       print "sequence mismatch: molecule", chain_n_residues(chain_id, imol), "new sequences:", len(sequence)
   else:
       make_backup(imol)
       backup_mode = backup_state(imol)
       turn_off_backup(imol)
       baddies = 0
       for ires in sequence : 
           res = mutate_single_residue_by_serial_number(ires, chain_id, imol, sequence_list[ires])
           if (res != 1) :
               baddies += 1

       set_have_unsaved_changes(imol)
       if (backup_mode == 1) :
           turn_on_backup(imol)
       update_go_to_atom_window_on_changed_mol(imol)
       graphics_draw()

def mutate_residue_range(imol,chain_id,start_res_no,stop_res_no,sequence):

   make_backup(imol)
   n_residues = stop_res_no - start_res_no + 1
   if (len(sequence) != n_residues) :
       print "sequence length mismatch:", len(sequence), n_residues
   else:
       backup_mode = backup_state(imol)
       turn_off_backup(imol)
       multi_mutate(mutate_single_residue_by_seqno, imol, start_res_no, chain_id, sequence)
       set_have_unsaved_changes(imol)
       if (backup_mode == 1) :
           turn_on_backup(imol)
       update_go_to_atom_window_on_changed_mol(imol)
       graphics_draw()

def mutate_and_autofit_residue_range(imol,chain_id,start_res_no,stop_res_no,sequence):
   
   mutate_residue_range(imol,chain_id,start_res_no,stop_res_no,sequence)
   mol_for_map = apply(imol_refinement_map)
#   print "mol_for_map is : ",mol_for_map
   if (mol_for_map >= 0) :
       backup_mode = backup_state(imol)
       turn_off_backup(imol)
       for ires in range(len(sequence)) :
          clash = 1
          altloc = ""
          inscode = ""
          resno = ires + start_res_no
          print "auto-fit-best-rotamer ",resno,altloc,inscode,chain_id,imol,mol_for_map,clash
          score = auto_fit_best_rotamer(resno,altloc,inscode,chain_id,imol,mol_for_map,clash,0.5)
          print "   Best score: ",score
#          number_list(start_res_no,stop_res_no)
       if (backup_mode == 1) :
           turn_on_backup(imol)
   else:
       print "WARNING:: no map set for refinement.  Can't fit"

def mutate_and_auto_fit(residue_number,chain_id,mol,mol_for_map,residue_type):

   mutate(residue_number,chain_id,mol,residue_type)
   auto_fit_best_rotamer(residue_number,"","",chain_id,mol,mol_for_map,0,0.5)

def multi_mutate(mutate_function,imol,start_res_no,chain_id,sequence):

   baddies = 0
   for ires in range(len(sequence)) :
#       print "start_res_no,ires, chain_id,imol,sequence[ires] are",start_res_no, ires,chain_id,imol,sequence[ires]
       result = mutate_function(start_res_no+ires,"",chain_id,imol,sequence[ires])
       if (result != 1) :
           baddies += 1
       print "multi_mutate of", len(sequence), "residues had", baddies, "errors"
