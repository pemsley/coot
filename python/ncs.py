# ncs.py
# Copyright 2007, 2008 by Bernhard Lohkamp
# copyright 2008 The University of York
# Copyright 2007 by Paul Emsley, The University of York
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


def find_first_model_molecule():

    imols = molecule_number_list()
    if not imols:
        return False
    else:
        for molecule in imols:
            if valid_model_molecule_qm(molecule):
                return molecule

# Skip the residue in the next chain (typically of a molecule with
# NCS) with the same residue number.  If on the last chain, then wrap
# to beginning.  If it can't find anything then don't move (and put a
# message in the status bar)
#
def skip_to_next_ncs_chain(direction):

  import types

  # Given a chain-id and a list of chain-ids, return the chain-id of
  # the next chain to be jumped to (use wrapping).  If the list of
  # chain-ids is less then length 2, return False.
  #
  def skip_to_chain_internal(this_chain_id, chain_id_list):
    # print "this_chain_id: ", this_chain_id
    if len(chain_id_list) < 2:
       return False
    else:
        # we do it differnt to Paul's function, as I dont understand it and
        # feel that this is easier and equally good!?
        # Agreed, it does the same thing and more obviously correct 20100126-PE
        current_chain_index = chain_id_list.index(this_chain_id)
	if current_chain_index == len(chain_id_list)-1:	# last chain
            # return the first chain
            return chain_id_list[0]
	else:
            # return the next chain
            return chain_id_list[current_chain_index + 1]


  def skip_to_chain(imol, this_chain_id, chain_id_list):

      # Given a chain-id and a list of chain-ids, return the chain-id of
      # the next chain to be jumped to (use wrapping).  If the list of
      # chain-ids is less then length 2, return False.
      #
      chain_guess = skip_to_chain_internal(this_chain_id, chain_id_list)

      if ((not type(chain_guess) is types.StringType)):
          return chain_guess
      elif (is_solvent_chain_qm(imol, chain_guess)):
          skip_to_chain(imol, chain_guess, chain_id_list)
      else:
          return chain_guess

  def get_chain_id_list(imol, this_chain_id):
      att = ncs_chain_ids(imol)
      if (not att):
          return chain_ids(imol)
      else:
          for attempt in att:
              if (this_chain_id in attempt):
                  return attempt
          # we havent found this_chain_id in ncs chains, so return all chains
          return chain_ids(imol)

  # First, what is imol? imol is the go to atom molecule
  imol = go_to_atom_molecule_number()
  this_chain_id = go_to_atom_chain_id()
  chains = get_chain_id_list(imol, this_chain_id)
  # try to get ghosts to be able to apply orientations
  make_ncs_ghosts_maybe(imol)
  found_atom_state = 0
  if (not chains):
      print "BL WARNING:: empty set of chains!!! This should never happen."
      msg_txt = "BL WARNING:: somehow there are no chains in mol " + str(imol)
      msg_txt+= "\n Try \"p\" to update the Go To Atom molecule and then\n"
      msg_txt+= "skip again. Good luck skipping!\n"
      info_dialog(msg_txt)
      # what shall we do now? Bail out I guess. But this should never happen!
      return
  next_chain = skip_to_chain(imol, this_chain_id, chains)

  try_next_chain = next_chain
  while (not try_next_chain == this_chain_id):

      # OK, stop trying for next chain if we have looped round
      # (as it were) so that we are back at the starting chain:
      # e.g. consider the case: ["A" is protein, "B" is water,
      # "C" is ligand]
      #
      if (not try_next_chain):
          add_status_bar_text("No 'NCS Next Chain' found")
          break
      else:
          if (not (try_next_chain == this_chain_id)):
              found_atom_state = set_go_to_atom_chain_residue_atom_name_no_redraw(
                  try_next_chain,
                  go_to_atom_residue_number(),
                  go_to_atom_atom_name(),
                  0)
              
              # now, did that set-go-to-atom function work (was there a
              # real atom)?  If not, then that could have been the ligand
              # chain or the water chain that we tried to go to.  We want
              # to try again, and we shbould keep trying again until we get
              # back to this-chain-id - in which case we have a "No NCS
              # Next Chain atom" status-bar message.

          if (found_atom_state == 0):
              # then we did *not* find the atom, e.g. next-chain was
              # the water chain
              try_next_chain = skip_to_chain(imol, try_next_chain, chains)
          else:
              # otherwise all was hunkey-dorey
              # set the orientation
              forward_flag = 0
              if (direction == "forward"):
                  forward_flag = 1
              apply_ncs_to_view_orientation_and_screen_centre(imol, this_chain_id, next_chain, forward_flag)
              break
  
            

# A function inspired by a question from Bill Scott.  He wanted to
# RNA ghosts.  Because RNA does not work with SSM, we need to define
# the matrix manually.  Let's make a copy of given imol and get
# the rtop from that.  Typical usage manual_ncs_ghosts(0, 1, 10, "A", "C")
# 
def single_manual_ncs_ghosts(imol, resno_start, resno_end, ref_chain, peer_chain):

	imol_copy = copy_molecule(imol)
	clear_lsq_matches()
	add_lsq_match(resno_start, resno_end, ref_chain, 
                      resno_start, resno_end, peer_chain, 0) # ref mov - all atoms
	rtop = apply_lsq_matches(imol_copy, imol_copy)
	close_molecule(imol_copy)
	if (not rtop):
            print "Failed to get matching matrix"
	else:
            clear_ncs_ghost_matrices(imol)
            set_draw_ncs_ghosts(imol, 1)
            args = [imol, peer_chain, ref_chain] + rtop[0] + rtop[1]
            add_ncs_matrix(*args)

# chain-id-list is ["A", "B", "C", "D"], i.e. the
# reference/target/master chain-id first and then the peers.  This
# allows us to add many peers at the same time (unlike above
# function).
#
def manual_ncs_ghosts(imol, resno_start, resno_end, chain_id_list):

    from types import ListType

    if (type(chain_id_list) is ListType):
        if (len(chain_id_list) > 1):
            # OK, OK, to the standard SSM-based NCS matrices are bad for this
            # molecule, lets use LSQ.
            clear_ncs_ghost_matrices(imol)
            imol_copy = copy_molecule(imol)
            for chain_id in chain_id_list[1:]:       # I dont think we need to superpose A onto A?!
                clear_lsq_matches()
                add_lsq_match(resno_start, resno_end, chain_id_list[0],
                              resno_start, resno_end, chain_id, 1)
                rtop = apply_lsq_matches(imol, imol_copy)
                if (not rtop):
                    print "Failed to get LSQ matching matrix", chain_id
                else:
                    master = chain_id_list[0]
                    print "chain_id %s master %s rtop: %s" %(chain_id, master, rtop)
                    if master:  # should be string or False, so ok I belive
                        args = [imol, chain_id, master] + rtop[0] + rtop[1]
                        add_ncs_matrix(*args)
            close_molecule(imol_copy)
            set_draw_ncs_ghosts(imol, 1)

# Update NCS ghosts based on local environment (residues within 6A of
# (and including) the active residue).
#
# Typically one would bind this function to a key.
#
def update_ncs_ghosts_by_local_sphere():
    """Update NCS ghosts based on local environment (residues within 6A of
    (and including) the active residue).
    
    Typically one would bind this function to a key."""

    with UsingActiveAtom() as [aa_imol, aa_chain_id, aa_res_no, aa_ins_code,
                               aa_atom_name, aa_alt_conf]:
        # clear_ncs_ghost_matrices(aa_imol)
        ghost_ncs_chain_ids = ncs_chain_ids(aa_imol)
        if (isinstance(ghost_ncs_chain_ids, list)):
            ghost_chain_id_list = ghost_ncs_chain_ids[0]
            if (isinstance(ghost_chain_id_list, list)):
                clear_ncs_ghost_matrices(aa_imol)
                for chain_id in ghost_chain_id_list[1]:
                    imol_copy = copy_molecule(aa_imol)
                    set_mol_displayed(imol_copy, 0)
                    set_mol_active(imol_copy, 0)
                    clear_lsq_matches()
                    active_residue_spec = [aa_chain_id, aa_res_no, aa_ins_code]
                    near_residues = residues_near_residue(aa_imol,
                                                          active_residue_spec,
                                                          6)
                    sphere_residues = [active_residue_spec, near_residues]
                    for residue_spec in sphere_residues:
                        res_no = residue_spec_to_res_no(residue_spec)
                        add_lsq_match(res_no, res_no, ghost_chain_id_list[0],
                                      res_no, res_no, chain_id, 1)

                    rtop = apply_lsq_matches(aa_imol, imol_copy)
                    close_molecule(imol_copy)
                    if not rtop:
                        print "Failed to get LSQ matching matrix", chain_id
                    else:
                        master = ghost_chain_id_list[0]
                        if (isinstance(master, str)):
                            args = [aa_imol, chain_id, master] + rtop[0] + rtop[1]
                            add_ncs_matrix(*args)

            
# Return the first master chain id (usually there is only one of course) or False
#
def ncs_master_chain_id(imol):

    from types import ListType
    cids_ls = ncs_chain_ids(imol)
    if type(cids_ls) is not ListType:
        return False
    else:
        if (len(cids_ls) == 0):
            return False
        else:
            cids = cids_ls[0]
            if (len(cids) == 0):
                return False
            else:
                return cids[0]

# This was designed to create an NCS copy of a ligand (or range of
# residues) in the active site of one chain to the as yet unoccupied
# active site of another, i.e. it makes a NCS ligand "D"1 that is a NCS
# copy of ligand "C"1 using an NCS operator that maps protein chain "A"
# onto chain "B".
# 
def ncs_ligand(imol_protein, ncs_master_chain_id,
               imol_ligand, chain_id_ligand,
               resno_ligand_start, resno_ligand_stop):
    
    # find ghost in ghosts that has a chain-id matching
    # chain-id-protein and get its rtop.  Return #f on not finding the
    # ghost
    def rtop_from_ghost_with_chain_id(ghosts, chain_id):
        for ghost in ghosts:
            if (ghost[1] == chain_id):
                return ghost[3]
        return False

    chain_ids_from_ncs = ncs_chain_ids(imol_protein)
    if chain_ids_from_ncs:
        ligand_selection_string = "//" + chain_id_ligand + "/" + str(resno_ligand_start) + \
                                  "-" + str(resno_ligand_stop)
        imol_ligand_fragment = new_molecule_by_atom_selection(imol_ligand, ligand_selection_string)
        ghosts = ncs_ghosts(imol_protein)
        for chain_ids in chain_ids_from_ncs:
            if (chain_ids[0] == ncs_master_chain_id):
                peer_chains = chain_ids[1:len(chain_ids)]
                candidate_name = "Candidate NCS-related ligand"
                for chain_id_protein in peer_chains:
                    rtop = rtop_from_ghost_with_chain_id(ghosts, chain_id_protein)
                    if (not rtop):
                        print "Opps - ncs-ligand: Missing ghost rt-op!"
                        info_dialog("Opps - ncs-ligand: Missing ghost rt-op!")
                    else:
                        new_lig_mol = copy_molecule(imol_ligand_fragment)
                        transform_coords_molecule(new_lig_mol, inverse_rtop(rtop))
                        set_molecule_name(new_lig_mol,
                                          str(new_lig_mol) +
                                          ": " +
                                          candidate_name +
                                          " to protein chain " +
                                          chain_id_protein)

        def test_func(imol):
            if (not valid_model_molecule_qm(imol)):
                return False
            else:
                name = molecule_name(imol)
                if (candidate_name in name):
                    ls = [name]
                    for i in molecule_centre(imol):
                        ls.append(i)
                    return ls
                else:
                    return False
        molecules_matching_criteria(lambda imol: test_func(imol))     
                
                
                
