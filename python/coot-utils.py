# Makefile.am
# 
# Copyright 2005, 2006, 2007 The University of York
# Copyright 2005, 2006, 2007 Bernhard Lohkamp
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA

#;; A silly (but somewhat amusing) demo that changes the background colour
#;; over a range of colours.
#;;
#;; The scene is spun too, for extra eye-candy (so have a molecule
#;; loaded before you run this).
#;;
#;; What's the point?  It's to show that we can display the result of
#;; an arbitarily complex computation, i.e, we have a real extension
#;; language, not just a list of commands (like some other molecular
#;; graphics programs).

# Copyright Bernhard Lohkamp
# adapted from coot-utils.scm
# so far only the things we really need

def set_virtual_trackball_type(type):
    if (type == "flat"):
        vt_surface(1)
    elif (type == "spherical-surface"):
        vt_surface(0)
    else:
        print "virtual trackball type",type,"not understood"

def number_list(a,b):
    result = []
    if a == b: 
       result.append(a)
       return result
    elif a > b : return result
    else :  
        while a <=b :
           result.append(a)
           a = a + 1
        return result

def rotation_centre():
   return [rotation_centre_position(0),
           rotation_centre_position(1),
           rotation_centre_position(2)]

def molecule_centre(imol):
   return [molecule_centre_internal(imol,0),
           molecule_centre_internal(imol,1),
           molecule_centre_internal(imol,2)]


def move_molecule_to_screen_centre(imol):
    #; We need to know what the current molecule centre for imol is.
    #; That is not available for export at the moment, but it should
    #; be.
    current_mol_centre = molecule_centre(imol)
    rotate_centre = rotation_centre()
    translate_molecule_by(imol,(rotate_centre[0]-current_mol_centre[0]),
                               (rotate_centre[1]-current_mol_centre[1]),
                               (rotate_centre[2]-current_mol_centre[2]))

def transform_map(imol,mat,trans,about_pt):

    transform_map_raw(imol,mat[0],mat[1],mat[2],
                           mat[3],mat[4],mat[5],
                           mat[6],mat[7],mat[8],
                           trans[0],trans[1],trans[2],
                           about_pt[0],about_pt[1],about_pt[2])


def chain_ids(imol):

    chain_id_is = []
    number_of_chains = n_chains(imol)
    for chain_no in range(number_of_chains):
        chain_id_is.append(chain_id(imol,chain_no))
    return chain_id_is

def is_solvent_chain_qm(imol,chain_id):
    if (is_solvent_chain_p(imol,chain_id)==1): return True

def valid_model_molecule_qm(imol): 
    if (is_valid_model_molecule(imol)==1): return True

def residue_exists_qm(imol,chain_id,resno,ins_code): 
    if (does_residue_exist_p(imol,chain_id,resno,ins_code)==1): return True

#not in 0.0.31
def centre_of_mass(imol):
    centre = centre_of_mass_string(imol)
    if (centre == 0):
       print "molecule number",imol,"is not valid"
#BL and now Paul has some call-with-input-string which I dunno whatto do with
#or need to do???

#Paul says:
#; backups wrapper: doesn't work currently, I think.  More cleverness required.
def with_no_backups(imol):
    backup_mode = backup_state(imol)
    turn_off_backup(imol)
    if backup_mode == 1:
       turn_on_backup(imol)

#Paul says:
#;;; This is not really a util, perhaps it should be somewhere else?
#not in 0.0.31
def print_sequence(imol):
    for chain in chain_ids(imol):
       print_sequence(imol,chain)





