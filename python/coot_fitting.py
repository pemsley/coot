# fitting.py
# Copyright 2005, 2006, 2007 by Bernhard Lohkamp 
# Copyright 2004, 2005, 2006, 2008, 2009 by The University of York
# Copyright 2009 by Bernhard Lohkamp
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

import coot
import coot_utils

# import coot_gui - coot_gui - this file is for scripting
# import coot_toolbuttons # 20220910-PE FIXME later

# For each residue in the protein (molecule number @var{imol}), do a
# rotamer fit and real-space refinement.  Update the graphics and
# rotate the scene at each residue for eye candy goodness.
#
# Note that residue with alt confs do not undergo auto-fit-rotamer.
# This is because that autofit-rotamer then refine will tend to put
# both rotamers into the same place.  Not good.  It seems a
# reasonable expectation that residues with an alternate conformation
# are already resonably well-fitted.  So residues with alternate
# conformations undergo only real space refinement.
#
# This is simple-minded and outdated now we have the interruptible
# version (below).
#
def fit_protein(imol, rotamer_only=False, animate=True):

    coot.set_go_to_atom_molecule(imol)
    coot.make_backup(imol) # do a backup first
    backup_mode = coot.backup_state(imol)
    imol_map  = coot.imol_refinement_map()
    replacement_state = coot.refinement_immediate_replacement_state()

    if imol_map == -1:
        coot.info_dialog("Oops.  Must set a map to fit")
    else:

        coot.turn_off_backup(imol)
        coot.set_refinement_immediate_replacement(1)
          
        for chain_id in coot_utils.chain_ids(imol):
         if (not coot_utils.is_solvent_chain_qm(imol,chain_id)):
             n_residues = coot.chain_n_residues(chain_id,imol)
             print("There are %(a)i residues in chain %(b)s" % {"a":n_residues,"b":chain_id})

             for serial_number in range(n_residues):
                  
                res_name = coot.resname_from_serial_number(imol, chain_id, serial_number)
                res_no = coot.seqnum_from_serial_number(imol, chain_id, serial_number)
                ins_code = coot.insertion_code_from_serial_number(imol, chain_id, serial_number)
                res_atoms = coot.residue_info_py(imol, chain_id, res_no, ins_code)
                
                if (ins_code is not None):
                    if (len(res_atoms) > 3):
                    #if (not res_name=="HOH"): not needed if only refining 3 or  more atoms
                        for alt_conf in coot_utils.residue_alt_confs(imol, chain_id, res_no, ins_code):
                            print("centering on ",chain_id,res_no," CA")
                            coot.set_go_to_atom_chain_residue_atom_name(chain_id,res_no,"CA")
                            if animate:
                                coot.rotate_y_scene(30, 0.3) # n-frames frame-interval(degrees)
                            if (alt_conf == ""):
                                coot.auto_fit_best_rotamer(res_no, alt_conf, ins_code, chain_id, imol,
                                                      imol_map, 1, 0.1)
                            if (imol_map >= 0 and
                                not rotamer_only):
                                coot.refine_zone(imol, chain_id, res_no, res_no, alt_conf)
                                coot.accept_regularizement()
                            if animate:
                                coot.rotate_y_scene(30, 0.3)
      
    if (replacement_state == 0):
        coot.set_refinement_immediate_replacement(0)
    if (backup_mode == 1):
        coot.turn_on_backup(imol)

# Paul: 20090517: thinking about making the fit-protein function
# interruptible with a toolbar button press.  How do we do that? 
# fit_protein needs to be split into 2 parts, one, that generates a
# list of residues specs the other that does a refinement given a
# residue spec, then we run and idle function that calls
# fit_residue_by_spec or each spec in turn and at the end return
#

# These 2 variables are used by multi-refine function(s), called by
# idle functions to refine just one residue.
# 
global continue_multi_refine
global multi_refine_spec_list
global multi_refine_idle_proc
global multi_refine_stop_button
global multi_refine_cancel_button
global multi_refine_continue_button
global multi_refine_separator
continue_multi_refine = False
multi_refine_spec_list = []
multi_refine_idle_proc = False
multi_refine_stop_button = False
multi_refine_cancel_button = False
multi_refine_continue_button = False
multi_refine_separator = False

# Return a list of residue specs
#
# chain-specifier can be a string, where it is the chain of interest.
# or 'all-chains, where all chains are chosen.
#
def fit_protein_make_specs(imol, chain_specifier):
    specs_list = []
    if coot.is_valid_model_molecule(imol):
        pass
    chain_list = []
    if chain_specifier == 'all-chains':
        chain_list = coot_utils.chain_ids(imol)
    else:
        chain_list = [chain_specifier]
    for ch_id in chain_list:
        serial_number_list = list(range(coot.chain_n_residues(ch_id, imol)))
        for serial_number in serial_number_list:
            res_no   = coot.seqnum_from_serial_number(imol, ch_id, serial_number)
            ins_code = coot.insertion_code_from_serial_number(imol, ch_id, serial_number)
            if ins_code is not None:
                specs_list.append([imol, ch_id, res_no, ins_code])
        
    return specs_list

    
def fit_protein_make_specs_from_residue_range(imol, chain_id,
                                              res_no_start, res_no_end):
    if coot.valid_model_molecule_qm(imol):
        if res_no_end > res_no_start:
            ret = [[chain_id, res_no, ""] for res_no in range(res_no_start, res_no_end + 1 )]
            return ret
    return []
    
    
def fit_protein_fit_function(res_spec, imol_map):

    imol     = res_spec[0]
    chain_id = res_spec[1]
    res_no   = res_spec[2]
    ins_code = res_spec[3]

    res_name = coot.residue_name(imol, chain_id, res_no, ins_code)
    if isinstance(res_name, str):
        if (res_name != "HOH"):
            for alt_conf in coot_utils.residue_alt_confs(imol, chain_id, res_no, ins_code):
                print("centering on", chain_id, res_no, "CA")
                coot.set_go_to_atom_chain_residue_atom_name(chain_id, res_no, "CA")
                coot.rotate_y_scene(10, 0.3) # n_frames frame_interval(degrees)        
                res_atoms = coot.residue_info_py(imol, chain_id, res_no, ins_code)
                if (len(res_atoms) > 3):
                    # if (not res_name == "HOH"):
                    # not needed as we only refine more than 3 atom res
                    if (alt_conf == ""):
                        with coot_utils.NoBackups(imol):
                            coot.auto_fit_best_rotamer(res_no, alt_conf, ins_code, chain_id, imol,
                                                  imol_map, 1, 0.1)
                    if (coot_utils.valid_map_molecule_qm(imol_map)):
                        with coot_utils.NoBackups(imol):
                            with coot_utils.AutoAccept():
                                coot.refine_zone(imol, chain_id, res_no, res_no, alt_conf)
                    coot.rotate_y_scene(10, 0.3)


def fit_protein_stepped_refine_function(res_spec, imol_map, use_rama = False):

    imol     = res_spec[0]
    chain_id = res_spec[1]
    res_no   = res_spec[2]
    ins_code = res_spec[3]
    current_rama_state = coot.refine_ramachandran_angles_state()
    
    for alt_conf in coot_utils.residue_alt_confs(imol, chain_id, res_no, ins_code):
        if use_rama:
            coot.set_refine_ramachandran_angles(1)
        res_name = coot.residue_name(imol, chain_id, res_no, ins_code)
        if (not res_name == "HOH"):
            print("centering on", chain_id, res_no, "CA")
            coot.set_go_to_atom_chain_residue_atom_name_full(chain_id, res_no,
                                                        ins_code, "CA",
                                                        alt_conf)
            coot.rotate_y_scene(10, 0.3) # n_frames frame_interval(degrees)
            with coot_utils.NoBackups(imol):
                with coot_utils.AutoAccept():
                    coot.refine_auto_range(imol, chain_id, res_no, alt_conf)
            coot.rotate_y_scene(10, 0.3)    

    coot.set_refine_ramachandran_angles(current_rama_state)

    
def fit_protein_rama_fit_function(res_spec, imol_map):
    # BL says: make it more generic
    fit_protein_stepped_refine_function(res_spec, imol_map, True)

# func is a refinement function that takes 2 args, one a residue
# spec, the other the coot.imol_refinement_map.  e.g. fit_protein_fit_function
#
def interruptible_fit_protein(imol, func):

    # import gobject
    global multi_refine_spec_list
    global continue_multi_refine
    global multi_refine_idle_proc
    global multi_refine_stop_button
    global multi_refine_separator
    specs = fit_protein_make_specs(imol, 'all-chains')
    if specs:
        # lets make a backup before we start
        coot.make_backup(imol)

        if False:
           # 20230501-PE old stuff commented out on gtk3 branch merge
           # multi_refine_separator = coot_toolbuttons.add_coot_toolbar_separator()
           # multi_refine_stop_button = coot_gui.coot_toolbar_button("Stop", "stop_interruptible_fit_protein()", "gtk-stop")
           pass

        multi_refine_spec_list = specs
        def idle_func():
            global multi_refine_spec_list
            global continue_multi_refine
            global multi_refine_idle_proc
            global multi_refine_stop_button
            global multi_refine_separator
            if not multi_refine_spec_list:
                #set_visible_toolbar_multi_refine_stop_button(0)
                #set_visible_toolbar_multi_refine_continue_button(0)
                multi_refine_stop_button.destroy()
                multi_refine_stop_button = False
                multi_refine_separator.destroy()
                multi_refine_separator = False
                return False
            if continue_multi_refine:
                imol_map = coot.imol_refinement_map()
                func(multi_refine_spec_list[0], imol_map)
                del multi_refine_spec_list[0]
                return True
            else:
                # finish what we have been doing first before we stop
                coot.accept_regularizement()
                return False
        # gobject.idle_add(idle_func)
        multi_refine_idle_proc = idle_func
        coot.set_go_to_atom_molecule(imol)


# For each residue in chain chain-id of molecule number imol, do a
# rotamer fit and real space refinement of each residue.  Don't
# update the graphics while this is happening (which makes it faster
# than fit-protein, but much less interesting to look at).
#
def fit_chain(imol, chain_id):

    coot.make_backup(imol)
    backup_mode = coot.backup_state(imol)
    imol_map = coot.imol_refinement_map()
    alt_conf = ""
    replacement_state = coot.refinement_immediate_replacement_state()

    coot.turn_off_backup(imol)
    coot.set_refinement_immediate_replacement(1)

    if imol_map == -1:
       print("WARNING:: fit-chain undefined imol-map. Skipping!!")
    else:
       n_residues = coot.chain_n_residues(chain_id,imol)
       for serial_number in range(n_residues):
           res_name = coot.resname_from_serial_number(imol,chain_id,serial_number)
           res_no = coot.seqnum_from_serial_number(imol,chain_id,serial_number)
           ins_code = coot.insertion_code_from_serial_number(imol,chain_id,serial_number)
           if ins_code is not None:
               res_atoms = coot.residue_info_py(imol, chain_id, res_no, ins_code)
               if (len(res_atoms) > 3):  # actually then we dont need the water check any more?!
                   #if (not res_name == "HOH"):
                       print("centering on ", chain_id, res_no, " CA")
                       coot.set_go_to_atom_chain_residue_atom_name(chain_id,res_no,"CA")
                       coot.auto_fit_best_rotamer(res_no, alt_conf, ins_code,
                                                  chain_id, imol, imol_map, 1, 0.1)
                       if (imol_map >= 1):
                           coot.refine_zone(imol,chain_id,res_no,res_no,alt_conf)
                           coot.accept_moving_atoms()

    if replacement_state == 0:
        coot.set_refinement_immediate_replacement(0)
    if backup_mode == 1:
        coot.turn_on_backup(imol)

# As fit_chain, but only for a residue range from resno_start to resno_end
def fit_residue_range(imol, chain_id, resno_start, resno_end):

    coot.make_backup(imol)
    backup_mode = coot.backup_state(imol)
    imol_map = coot.imol_refinement_map()
    alt_conf = ""
    replacement_state = coot.refinement_immediate_replacement_state()

    coot.turn_off_backup(imol)
    coot.set_refinement_immediate_replacement(1)
    
    if (imol_map == -1):
       print("WARNING:: fit-chain undefined imol-map. Skipping!!")
    else:
       n_residues = coot.chain_n_residues(chain_id,imol)
       ins_code = ""
       for res_no in coot_utils.number_list(resno_start, resno_end):
           print("centering on ", chain_id, res_no, " CA")
           coot.set_go_to_atom_chain_residue_atom_name(chain_id, res_no, "CA")
           res_atoms = coot.residue_info_py(imol, chain_id, res_no, ins_code)
           if (len(res_atoms) > 3):
               coot.auto_fit_best_rotamer(res_no, alt_conf, ins_code, chain_id,
                                          imol, imol_map, 1, 0.1)
               if (imol_map >= 1):
                   coot.refine_zone(imol,chain_id,res_no,res_no,alt_conf)
                   coot.accept_regularizement()

    if (replacement_state == 0):
        coot.set_refinement_immediate_replacement(0)
    if (backup_mode == 1):
          coot.turn_on_backup(imol)


# For each residue in the solvent chains of molecule number
# @var{imol}, do a rigid body fit of the water to the density.
#
# BL says: we pass *args where args[0]=imol and args[1]=animate_qm (if there)
def fit_waters(imol, animate_qm = False):

    print("animate?:", animate_qm)
    imol_map = coot.imol_refinement_map()
    do_animate_qm = False
    if (animate_qm):
        do_animate_qm = True

    print("do_animate?: ", do_animate_qm)
 
    if (imol_map != -1):
        replacement_state = coot.refinement_immediate_replacement_state()
        backup_mode = coot.backup_state(imol)
        alt_conf = ""

        coot.turn_off_backup(imol)
        coot.set_refinement_immediate_replacement(1)
        coot.set_go_to_atom_molecule(imol)

        # refine waters
        for chain_id in coot_utils.chain_ids(imol):
            if (coot_utils.is_solvent_chain_qm(imol, chain_id)):
                n_residues = coot.chain_n_residues(chain_id, imol)
                print("There are %(a)i residues in chain %(b)s" % {"a":n_residues,"b":chain_id})
                for serial_number in range(n_residues):
                    res_no = coot.seqnum_from_serial_number(imol,chain_id,serial_number)
                    if do_animate_qm:
                        res_info = coot.residue_info_py(imol, chain_id, res_no, "")
                        if not res_info == []:
                            atom = res_info[0]
                            atom_name = atom[0]
                            coot.set_go_to_atom_chain_residue_atom_name(chain_id, res_no, atom_name)
                            coot.refine_zone(imol, chain_id, res_no, res_no, alt_conf)
                            coot.rotate_y_scene(30, 0.6)	# n-frames frame-interval(degrees)
                    else:
                        coot.refine_zone(imol,chain_id,res_no,res_no,alt_conf)
                    coot.accept_regularizement()

        if (replacement_state == 0):
            coot.set_refinement_immediate_replacement(0)
        if (backup_mode == 1):
            coot.turn_on_backup(imol)

    else:
        coot.add_status_bar_text("You need to define a map to fit the waters")


# BL thingy:
# as fit_waters, but will only fit a range of waters (start to end).
# speeds things up when refining a lot of waters
#
def fit_waters_range(imol, chain_id, start, end):

    imol_map = coot.imol_refinement_map()
    if (imol_map != -1):
       replacement_state = coot.refinement_immediate_replacement_state()
       backup_mode = coot.backup_state(imol)
       alt_conf = ""

       coot.turn_off_backup(imol)
       coot.set_refinement_immediate_replacement(1)

       if (coot_utils.is_solvent_chain_qm(imol,chain_id)):
          serial_number = start
          while serial_number <= end:
                res_no = serial_number
                coot.refine_zone(imol,chain_id,res_no,res_no,alt_conf)
                coot.accept_regularizement()
                serial_number += 1

       if (replacement_state == 0):
          coot.set_refinement_immediate_replacement(0)
       if (backup_mode == 1):
          coot.turn_on_backup(imol)

# Step through the residues of molecule number imol and at each step
# do a residue range refinement (unlike fit-protein for example,
# which does real-space refinement for every residue).
#
# The step is set internally to 2.
#
def stepped_refine_protein(imol, res_step = 2):

    import types
    from types import IntType

    imol_map = coot.imol_refinement_map()
    if (not coot_utils.valid_map_molecule_qm(imol_map)):            
        coot.info_dialog("Oops, must set map to refine to")
    else:
        def refine_func(chain_id, res_no):
            #print "centering on ",chain_id,res_no," CA"
            coot.set_go_to_atom_chain_residue_atom_name(chain_id, res_no, "CA")
            coot.rotate_y_scene(30, 0.3) # n-frames frame-interval(degrees)
            coot.refine_auto_range(imol, chain_id, res_no, "")
            coot.accept_regularizement()
            coot.rotate_y_scene(30,0.3)
        stepped_refine_protein_with_refine_func(imol, refine_func, res_step)


# refine each residue with ramachandran restraints
#
def stepped_refine_protein_for_rama(imol):

    imol_map = coot.imol_refinement_map()
    if (not coot_utils.valid_map_molecule_qm(imol_map)):
        coot.info_dialog("Oops, must set map to refine to")
    else:
        current_rama_state  = coot.refine_ramachandran_angles_state()
        def refine_func(chain_id, res_no):
            coot.set_go_to_atom_chain_residue_atom_name(chain_id, res_no, "CA")
            coot.refine_auto_range(imol, chain_id, res_no, "")
            coot.accept_regularizement()

        coot.set_refine_ramachandran_angles(1)
        stepped_refine_protein_with_refine_func(imol, refine_func, 1)
        coot.set_refine_ramachandran_angles(current_rama_state)

#
def stepped_refine_protein_with_refine_func(imol, refine_func, res_step):

    coot.set_go_to_atom_molecule(imol)
    coot.make_backup(imol)
    backup_mode = coot.backup_state(imol)
    alt_conf = ""
    imol_map = coot.imol_refinement_map()
    replacement_state = coot.refinement_immediate_replacement_state()

    if (imol_map == -1):
        # actually shouldnt happen as we set check the map earlier...
        coot.add_status_bar_text("Oops.  Must set a map to fit")
    else:
        coot.turn_off_backup(imol)
        coot.set_refinement_immediate_replacement(1)
        res_step = int(res_step)
        if (res_step <= 1):
            coot.set_refine_auto_range_step(1)
            res_step = 1
        else:
            coot.set_refine_auto_range_step(int(res_step / 2))

        for chain_id in coot_utils.chain_ids(imol):
            n_residues = coot.chain_n_residues(chain_id,imol)
            print("There are %(a)i residues in chain %(b)s" % {"a":n_residues,"b":chain_id})
            
            for serial_number in range(0, n_residues, res_step):
                res_name = coot.resname_from_serial_number(imol, chain_id, serial_number)
                res_no = coot.seqnum_from_serial_number(imol, chain_id, serial_number)
                ins_code = coot.insertion_code_from_serial_number(imol, chain_id, serial_number)
                if ins_code is not None:
                    print("centering on ", chain_id, res_no, " CA")
                    refine_func(chain_id, res_no)
                
        if (replacement_state == 0):
            coot.set_refinement_immediate_replacement(0)
        if (backup_mode == 1):
            coot.turn_on_backup(imol)


# The gui that you see after ligand finding. 
# 
def post_ligand_fit_gui():

     def test(imol):
       if not coot.is_valid_model_molecule(imol):
          return False
       else:
          name = coot.molecule_name(imol)
          if "Fitted ligand" in name:
             list = [name]
             # BL says: I think there must be a cleverer way!
             for i in coot_utils.molecule_centre(imol):
                list.append(i)
             return list
          else:
             return False 
     molecules_matching_criteria(lambda imol: test(imol))

# test_func is a function given one argument (a molecule number) that
# returns either False if the condition is not satisfied or something
# else if it is.  And that "something else" can be a list like
# [label, x, y, z]
# or 
# ["Bad Chiral", 0, "A", 23, "", "CA", "A"]
# 
# It is used in the create a button label and "what to do when the
# button is pressed".
# 
def molecules_matching_criteria(test_func):

    # 2026-02-18-PE function deleted
    pass

# This totally ignores insertion codes.  A clever algorithm would
# need a re-write, I think.  Well, we'd have at this end a function
# that took a chain_id res_no_1 ins_code_1 res_no_2 ins_code_2 
# 
# And refine-zone would need to be re-written too, of course.  So
# let's save that for a rainy day (days... (weeks)).
# 
# BL says:: this seems to ignore backup!
def refine_active_residue_generic(side_residue_offset):

    active_atom = active_residue()

    if not active_atom:
       print("No active atom")
    else:
       imol       = active_atom[0]
       chain_id   = active_atom[1]
       res_no     = active_atom[2]
       ins_code   = active_atom[3]
       atom_name  = active_atom[4]
       alt_conf   = active_atom[5]
    
       print("active-atom:", active_atom)
       imol_map = coot.imol_refinement_map()
       replacement_state = coot.refinement_immediate_replacement_state()
       if imol_map == -1:
          coot.info_dialog("Oops.  Must Select Map to fit to!")
       else:
          coot.set_refinement_immediate_replacement(1)
          coot.refine_zone(imol, chain_id, res_no - side_residue_offset,
                                      res_no + side_residue_offset,
                                      alt_conf)
          coot.accept_regularizement()
       if replacement_state == 0:
          coot.set_refinement_immediate_replacement(0)

# Function for keybinding. Speaks for itself
def refine_active_residue():
    refine_active_residue_generic(0)

# And another one
def refine_active_residue_triple():
    refine_active_residue_generic(1)


# For just one (this) residue, side-residue-offset is 0.
# 
def manual_refine_residues(side_residue_offset):

    active_atom = coot.active_residue()

    if not active_atom:
       print("No active atom")
    else:
       imol       = active_atom[0]
       chain_id   = active_atom[1]
       res_no     = active_atom[2]
       ins_code   = active_atom[3]
       atom_name  = active_atom[4]
       alt_conf   = active_atom[5]

    imol_map = coot.coot.imol_refinement_map()

    if (imol_map == -1):
        coot.info_dialog("Oops.  Must Select Map to fit to!")
    else:
        coot.refine_zone(imol, chain_id,
                         res_no - side_residue_offset,
                         res_no + side_residue_offset,
                         alt_conf)

# generic spherical refinement (use_map) (or regularization, dont use map):
def sphere_refine_regularize_generic(use_map=True, radius=3, expand=False):
    # from types import ListType
    active_atom = coot.active_residue_py()
    if (not active_atom):
        coot.add_status_bar_text("No active residue")
    else:
        if (use_map and not coot_utils.valid_map_molecule_qm(coot.imol_refinement_map())):
            coot.show_select_map_dialog()
        else:
            imol      = active_atom[0]
            chain_id  = active_atom[1]
            res_no    = active_atom[2]
            ins_code  = active_atom[3]
            atom_name = active_atom[4]
            alt_conf  = active_atom[5]
            centred_residue = active_atom[1:4]
            other_residues = coot.residues_near_residue_py(imol, centred_residue, radius)
            coot_utils.all_residues = [centred_residue]
            coot_utils.all_residues += other_residues

            # extend?
            if expand:
                print("in sphere_refine_regularize_generic, all_residues is", all_residues, "using radius", radius)
                coot_utils.all_residues.sort()
                tmp_ls = coot_utils.all_residues[:]
                for res in tmp_ls:
                    before_res = res[:]
                    after_res = res[:]
                    before_res[1] = before_res[1]-1
                    after_res[1] = after_res[1]+1
                    if not before_res in coot_utils.all_residues:
                        coot_utils.all_residues.append(before_res)
                    if not after_res in coot_utils.all_residues:
                        coot_utils.all_residues.append(after_res)
                coot_utils.all_residues.sort()  # not needed

            print("imol: %s residues: %s" %(imol, coot_utils.all_residues))
            if use_map:
                # don't use 'soft-mode/hard-mode' at the moment
                # (not sure how to integrate weight change into dragged refinement)
                coot.refine_residues_with_modes_with_alt_conf_py(imol, coot_utils.all_residues, "",
                                                                 '#soft-mode/hard-mode', False, False)
            else:
                coot.regularize_residues_py(imol, coot_utils.all_residues)

                # Sphere refinement (around radius)
#
def sphere_refine(radius=4.5, expand=False):
    sphere_refine_regularize_generic(True, radius, expand)

def sphere_refine_plus(radius=4.5):
    sphere_refine(radius, True)

# Sphere regularization (as above)
def sphere_regularize(radius=4.5, expand=False):
    sphere_refine_regularize_generic(False, radius, expand)

def sphere_regularize_plus(radius=4.5):
    sphere_regularize(radius, True)


def refine_tandem_residues():
    active_atom = coot.closest_atom_simple_py() # active_atom returns the CA if it can
    if not active_atom:
       print("No active atom")
    else:
       imol       = active_atom[0]
       chain_id   = active_atom[1]
       res_no     = active_atom[2]
       ins_code   = active_atom[3]
       atom_name  = active_atom[4]
       alt_conf   = active_atom[5]
       specs = []
       for ires in range(res_no-3, res_no+4):
           try:
              # test if the residue exists by looking for a residue name
              rn = coot.residue_name(imol, chain_id, ires, ins_code)
              if len(rn) > 0:
                  specs.append([chain_id, ires, ins_code])
           except TypeError as e:
               pass # no need to tell us
       coot.refine_residues(imol, specs)


# Pepflip the active residue - needs a key binding
#
def pepflip_active_residue():
    active_atom = coot.closest_atom_simple_py() # active_atom returns the CA if it can
    if not active_atom:
       print("No active atom")
    else:
       imol       = active_atom[0]

       ca = coot.closest_atom_raw() # don't map to CA

       chain_id   = ca[1]
       res_no     = ca[2]
       ins_code   = ca[3]
       atom_name  = ca[4]
       alt_conf   = ca[5]

       if (atom_name == " N  "): # PDBv3 fixme
           res_no -= 1;
       coot.pepflip(imol, chain_id, res_no, ins_code, alt_conf)

    

# Another cool function that needs a key binding
#
def auto_fit_rotamer_active_residue():

    active_atom = coot.active_residue()

    if not active_atom:
       print("No active atom")
    else:
       imol       = active_atom[0]
       chain_id   = active_atom[1]
       res_no     = active_atom[2]
       ins_code   = active_atom[3]
       atom_name  = active_atom[4]
       alt_conf   = active_atom[5]
   
       print("active-atom:", active_atom)
       imol_map = coot.imol_refinement_map()
       replacement_state = coot.refinement_immediate_replacement_state()
       if imol_map == -1:
          coot.info_dialog("Oops.  Must Select Map to fit to!")
       else:
          coot.auto_fit_best_rotamer(res_no, alt_conf, ins_code, chain_id, imol, imol_map, 1, 0.1)


# Backrub rotamers for chain. After alignment mutation we should run this.
#
def backrub_rotamers_for_chain(imol, ch_id):

    """Backrub rotamers for chain. After alignment mutation we should run this."""

    coot.set_rotamer_search_mode(2) # ROTAMERSEARCHLOWRES
    coot.make_backup(imol)

    with coot_utils.NoBackups(imol):
        n_times = 2
        imol_map = coot.imol_refinement_map()
        if coot_utils.valid_map_molecule_qm(imol_map):
            n_res = coot.chain_n_residues(ch_id, imol)
            for i_round in range(n_times):
                for serial_number in range(n_res):
                    res_name = coot.resname_from_serial_number(imol, ch_id, serial_number)
                    res_no = coot.seqnum_from_serial_number(imol, ch_id, serial_number)
                    ins_code = coot.insertion_code_from_serial_number(imol, ch_id, serial_number)
                    print("debug res_name", res_name)
                    print("debug res_no", res_no)
                    print("debug ins_code", ins_code)
                    print("debug ch_id", ch_id)
                    if isinstance(ins_code, str):   # valid residue check :-)
                        if not res_name == "HOH":
                            # coot.auto_fit_best_rotamer(res_no, "", ins_code, ch_id, imol, imol_map, 1, 0.1)
                            coot.auto_fit_best_rotamer(imol, ch_id, res_no, ins_code, "", imol_map, 1, 0.1)


# Restrain the atoms in imol (in give range selection) to
# corresponding atoms in imol_ref.
# 
# atom_sel_type is either 'all' 'main-chain' or 'ca'
#
def add_extra_restraints_to_other_molecule(imol, chain_id,
                                           resno_range_start, resno_range_end,
                                           atom_sel_type, imol_ref):

    for res_no in range(resno_range_start, resno_range_end):
        pass # guess should do seomething?!?!

def add_extra_start_pos_restraints(imol, residue_spec, esd):

    ri = coot.residue_info_py(imol,
                              coot_utils.residue_spec_to_chain_id(residue_spec),
                              coot_utils.residue_spec_to_res_no(residue_spec),
                              coot_utils.residue_spec_to_ins_code(residue_spec))
    for atom_info in ri:
        atom_name = coot_utils.residue_atom2atom_name(atom_info)
        alt_conf  = coot_utils.residue_atom2alt_conf(atom_info)
        coot.add_extra_start_pos_restraint(imol,
                                      coot_utils.residue_spec_to_chain_id(residue_spec),
                                      coot_utils.residue_spec_to_res_no(residue_spec),
                                      coot_utils.residue_spec_to_ins_code(residue_spec),
                                      atom_name, alt_conf, esd)

        
        pass # fill me
