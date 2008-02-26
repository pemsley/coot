# coot-utils.py 
# adapted from coot-utils.scm

# Copyright 2004, 2005, 2006, 2007 by Bernhard Lohkamp
# Copyright 2000 by Paul Emsley
# Copyright 2004, 2005, 2006, 2007 by Paul Emsley, The University of York

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

# 'Macro' to tidy up a set of functions to be run with automatic
# accepting of the refinement
# 
# funcs is a normal set of functions (not a thunk)
# 
def with_auto_accept(*funcs):

    replace_state = refinement_immediate_replacement_state()
    set_refinement_immediate_replacement(1)
    for f in funcs:
        func = f[0]
        args = f[1:len(f)]
        #print "BL DEBUG:: func and args", func, args
        func(*args)
        accept_regularizement()
    
    if (replace_state == 0):
        set_refinement_immediate_replacement(0)


# return a list of molecule numbers (closed and open)
# The elements of the returned list need to be tested against is_valid_model_molecule_qm
def molecule_number_list():
    n_molecules = graphics_n_molecules()
    if (n_molecules == 0):
       return []
    else:
       return number_list(0,(n_molecules-1))

# Return a list of molecules that are maps
# 
def map_molecule_list():

    map_list = []
    for i in range(graphics_n_molecules()):
       if is_valid_map_molecule(i):
          map_list.append(i)
    return map_list

# Set the virtual trackball behaviour.
#
# trackball @var{type} is a string: either 'flat' or 'spherical-surface'
#
def set_virtual_trackball_type(type):
    if (type == "flat"):
        vt_surface(1)
    elif (type == "spherical-surface"):
        vt_surface(0)
    else:
        print "virtual trackball type",type,"not understood"

#Is ls a list of strings? Return True or False
#
def list_of_strings_qm(ls):
    import types
    not_str =0
    if type(ls) is not ListType:
       return False
    else:
       for item in ls:
           if isinstance(item,types.StringTypes): pass
           else: not_str += 1
       if not_str == 0:
          return True
       else:
          return False

# string concat with spaces, @var{ls} must be a list of strings.
#
def string_append_with_spaces(ls):
    
    import string
    if ls:
       return string.join(ls)
    else:
       return [""] 

# The screen centre.
#
# return the rotation centre as a 3 membered list of numbers
# is python list [...] !!!
#
def rotation_centre():
   return [rotation_centre_position(0),
           rotation_centre_position(1),
           rotation_centre_position(2)]

# this is actually not essentail since python has these funtion(s)
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

# Return True or False
# adapted from find_exe
def command_in_path_qm(cmd):

  # test for command (see goosh-command-with-file-input description)
  # 
    import os,string

    file_ext = file_name_extension(cmd)
    if ((file_ext <> 'exe') and (os.name == 'nt')):
    	cmd += ".exe"
    # we shall check for full path names first
    if (os.path.isfile(cmd)):
     return True
    else:
     try:
       primary_path = os.environ["PATH"]
       for path in string.split(primary_path, os.pathsep):
           program_exe = os.path.join(path,cmd)
#           print "BL DEBUG:: program_exe is", program_exe
           if (os.path.isfile(program_exe)):
               return True
               break
           else:
               pass
       return False
     except:
       print "BL WARNING:: couldnt open $PATH"  # this shouldnt happen
  
# Where cmd is e.g. "refmac" 
#       args is ["HKLIN","thing.mtz"]
#       log_file_name is "refmac.log"      
#       data_list is ["HEAD","END"]
# 
# Return the exist status e.g. 0 or 1.
# 
def popen_command(cmd, args, data_list, log_file, screen_flag):
    import string, os

    if not(command_in_path_qm(cmd)):
       print "command ", cmd, " not found in $PATH!"
       print "BL INFO:: Maybe we'll find it somewhere else later..."
    else:
       cmd_execfile = find_exe(cmd,"CCP4_BIN","PATH")
       data_list_file = "data_list_file_tmp.txt"

       # make args string
       args_string = string.join(args)
    
       # write tmp input file
       input = file (data_list_file,'w')
       for data in data_list:
           input.write(data + '\n')
       input.close()

#       print "BL DEBUG:: popen command is", cmd_execfile + " " + args_string + " < " + data_list_file + " > " + log_file
       status = os.popen(cmd_execfile + " " + args_string + " < " + data_list_file + " > " + log_file, 'r')
       cmd_finished = status.close()
       if (cmd_finished != None):
          print "running command ", cmd, " failed!"
          return 0
       else:
#          print "cmd ", cmd ," seems to have run ok!"
          return 1
       os.remove(data_list_file)

# example usage:
# popen_command("mtzdump",["HKLIN","a.mtz"],["HEAD","END"],"test.log",0)
  
# "a.b.res" -> "a.b"
# file_name_sans_extension
def strip_extension(s):
    import os
    head, tail = os.path.splitext(s)
    return head

# What is the extension of file_name?
# 
# "a.pdb" -> "pdb"
# "" -> ""
#
def file_name_extension(file_name):
    
    import os, string
    root, ext = os.path.splitext(file_name)
    if ext:
       ext = string.lstrip(ext,'.')
       return ext
    else:
       return ""

# e.g. "a.pdb" -> "a-tmp.pdb"
#
def add_tmp_extension_to(file_name):
  
    import types
    if isinstance(file_name,types.StringTypes):
       root, ext = os.path.splitext(file_name)
       f = root + "-tmp" + ext
       return f
    else:
       return "tmp"

# /a/b.t -> b.t   d/e.ext -> e.ext
# file-name-sans-path
def strip_path(s):
    import os
    head, tail = os.path.split(s)
    return tail

# does s start with a "/" ?
# return True or False
def slash_start_qm(s):
    import types
    if isinstance(s, types.StringTypes):
       if len(s) > 0:
          if s.startswith("/"): return True
          else: return False
       else: return False
    return False

# return a string that contains the date/time
# e.g. "2006-01-02_2216.03"
#
def unique_date_time_str():
    import time
    lt = time.strftime("%Y-%m-%d_%H%M", time.localtime())
    return lt

# return a list that has only every-nth members; 
# e.g. @code{every_nth ([0,1,2,3,4,5,6,7,8],2)} -> [0,2,4,6,8]
#      @code{every_nth ([0,1,2,3,4,5,6,7,8],3)} -> [0,3,6]
# 
# @var{n} must be positive
def every_nth(ls, n):

    elements = range(0,len(ls),n)
    a =[]
    for i in elements:
        a.append(ls[i])
    return a


# multi_read pdb reads all the files matching
# @code{@emph{glob_pattern}} in
# directory @code{@emph{dir}}.  Typical usage of this might be:
# @code{multi_read_pdb("a*.pdb",".")}
# BL says: in windows dir needs the 'C:/' pattern, '/c/'won't work
def multi_read_pdb(glob_pattern, dir):
    import glob, os
    patt = os.path.normpath(dir+'/*.'+glob_pattern)
    all_files = glob.glob(patt)
    for file in all_files:
        print "BL INFO:: reading ", file
        read_pdb(file)

# read_pdb_all reads all the "*.pdb" files in the current directory.
#
def read_pdb_all():
    import glob, os
    recentre_status = recentre_on_read_pdb()
    set_recentre_on_read_pdb(0)
    patt = os.path.normpath(os.path.abspath(".")+'/*.pdb')
    all_files = glob.glob(patt)
    for file in all_files:
        print "BL INFO:: reading ", file
        read_pdb(file)
    set_recentre_on_read_pdb(recentre_status)

# return False if dir_name is a file or we can't do the mkdir
def coot_mkdir(dir_name):
  import os
  if (os.path.isfile(dir_name)):
     return False
  else:
     if (os.path.isdir(dir_name)):
       return True
     else:
       os.mkdir(dir_name)
       return True

# return the view matrix (useful for molscript, perhaps).
# BL says: like all matrices is a python list [...]
def view_matrix():
    return [get_view_matrix_element(row_number,column_number) for row_number in range(3) for column_number in range(3)]

# return the transposed view matrix (useful for molscript, perhaps).
# BL says: like all matrices is a python list [...]
def view_matrix_transp():
    return [get_view_matrix_element(column_number,row_number) for row_number in range(3) for column_number in range(3)]

# return the view quaternion
def view_quaternion():

	ret = map(get_view_quaternion_internal,[0,1,2,3])
	return ret

# Return the view number 
#
def add_view(position, quaternion, zoom, view_name):

	args = position + quaternion
	args.append(zoom)
	args.append(view_name)
	ret = add_view_raw(*args)
	return ret

# Convert a view matrix to a view quaternion to set Coot view internals.
#
def matrix2quaternion(m00, m10, m20, m01, m11, m21, m02, m12, m22):

    import math

    # From an idea by "Christian" at euclidianspace.com.  The
    # rotation matrix is special orthogonal, so (1 + trace) > 0. So
    # we can simply do a sqrt on the 4 trace variants.  Most of the
    # code here is simply recovering the sign.
    
    # return x with the sign of y
    def convert_sign(x, y):
        if (x > 0 and y > 0):
            return x
        elif (x < 0 and y > 0):
            return -x
        elif (x > 0 and y <0):
            return -x
        else:
            return x

    pw = 1 + m00 + m11 + m22
    px = 1 + m00 - m11 - m22
    py = 1 - m00 + m11 - m22
    pz = 1 - m00 - m11 + m22

    pr = []
    for v in [pw, px, py, pz]:
        if v < 0:
            v1 = 0
        else:
            v1 = v
        pr.append(math.sqrt(v1) / 2)

    ls = map(convert_sign, pr[1:], [m21 - m12, m02 - m20, m10 - m01])

    ls.append(pr[0])

    return ls

# e.g
# matrix2quaternion(0.0347695872187614, 0.773433089256287, 0.632923781871796,
#                   0.774806916713715, 0.379149734973907, -0.505885183811188,
#                  -0.631241261959076, 0.507983148097992, -0.586078405380249)
# -> 
# [-0.55715757608, -0.694704711, -7.549694273e-4, 0.45492890477] or similar

# Set the view matrix using matrix->quaternion. 
#
# Useful for using a view matrix from another program, perhaps.
#
def set_view_matrix(m00, m10, m20, m01, m11, m21, m02, m12, m22):

    set_view_quaternion(matrix2quaternion(m00, m10, m20,
                                          m01, m11, m21,
                                          m02, m12, m22))

# Miguel's molecular orientation axes
def miguels_axes():
    set_axis_orientation_matrix(*view_matrix())
    set_axis_orientation_matrix_usage(1)

# Return the molecule centre as a list of 3 numbers.
#
#  Note: mol_cen could contain values less than -9999.
#
def molecule_centre(imol):
   return [molecule_centre_internal(imol,0),
           molecule_centre_internal(imol,1),
           molecule_centre_internal(imol,2)]

# Move the centre of molecule number imol to the current screen centre
#
def move_molecule_to_screen_centre(imol):
  if valid_model_molecule_qm(imol):
    rotate_centre = rotation_centre()
    translate_molecule_by(imol,(rotate_centre[0]-molecule_centre(imol)[0]),
                               (rotate_centre[1]-molecule_centre(imol)[1]),
                               (rotate_centre[2]-molecule_centre(imol)[2]))
# This is a short name for the above.
# deftexi move_molecule_here
move_molecule_here = move_molecule_to_screen_centre

# Return a nine-membered list of numbers.
#
def identity_matrix():
    return [1,0,0,0,1,0,0,0,1]

# e.g. translation('x',2)
#  -> [2, 0, 0]
# return False on error
def translation(axis,length):
    import operator
# BL says: we dont check if axis is string, yet at least not directly
    if (operator.isNumberType(length)):
       if (axis=="x"):
          return [length,0,0]
       elif (axis=="y"):
          return [0,length,0]
       elif (axis=="z"):
          return [0,0,length]
       else:
          print "symbol axis: ", axis, " incomprehensible"
          return False
    else:
       print "incomprehensible length argument: ",length
       return False

# Rotate degrees about screen axis, where axis is either 'x', 'y' or 'z'.
# 
def rotate_about_screen_axis(axis,degrees):
    import math, operator

    def deg_to_rad(degs):
        return (degs * 3.1415926 /180.0)

    def simple_rotation_x(alpha):
        cos_alpha = math.cos(alpha)
        sin_alpha = math.sin(alpha)
        return [1,0,0,0,cos_alpha,-sin_alpha,0,sin_alpha,cos_alpha]

    def simple_rotation_y(alpha):
        cos_alpha = math.cos(alpha)
        sin_alpha = math.sin(alpha)
        return [cos_alpha,0,sin_alpha,0,1,0,-sin_alpha,0,cos_alpha]

    def simple_rotation_z(alpha):
        cos_alpha = math.cos(alpha)
        sin_alpha = math.sin(alpha)
        return [cos_alpha,-sin_alpha,0,sin_alpha,cos_alpha,0,0,0,]

# BL says: I dont know what the next 2 defines are for... 
# looks not used and/or useless to me
# seems that only 2nd matrix is used and not view_matrix!
    def vm():view_matrix()
    def mult(mat1,mat2):mat2
# end of uselessness...

    if (operator.isNumberType(degrees)): 
       if (axis=="x"):
          mult(view_matrix(),simple_rotation_x(deg_to_rad(degrees)))
       elif (axis=="y"):
          mult(view_matrix(),simple_rotation_y(deg_to_rad(degrees)))
       elif (axis=="z"):
          mult(view_matrix(),simple_rotation_z(deg_to_rad(degrees)))
       else:
          print "symbol axis: ", axis, " incomprehensible"
    else:
       print "incomprehensible length argument: ", degrees

# transform a coordinates molecule by a coot-rtop (which is a Python
# expression of a clipper::RTop)
#
def transform_coord_molecule(imol, rtop):

    transform_molecule_by()

# @code{transform_map(imol,mat,trans,about_pt,radius)}
#
# or @code{transform_map(imol,trans,about_pt,radius)} for a simple translation
#
# or @code{transform_map(imol,trans,radius)} when using the default 
# rotation-centre as the about-pt
#
def transform_map(*args):

    def tf(imol,mat,trans,about_pt,radius):
        transform_map_raw(imol,mat[0],mat[1],mat[2],
                               mat[3],mat[4],mat[5],
                               mat[6],mat[7],mat[8],
                               trans[0],trans[1],trans[2],
                               about_pt[0],about_pt[1],about_pt[2],
                               radius)
    if (len(args)==5):
       tf(args[0],args[1],args[2],args[3],args[4])
#no matrix specified:
    elif (len(args)==4):
       tf(args[0],identity_matrix(),args[1],args[2],args[3])
#no matrix or about point specified:
    elif (len(args)==3):
       tf(args[0],identity_matrix(),args[1],rotation_centre(),args[2])
    else:
       print "arguments to transform-map incomprehensible: args: ",args

# Define a map transformation function that obeys Lapthorns Law of
# NCS Handling Programs
# 
# typical usage: transform_map_using_lsq_matrix(1, "A", 10, 30, 0, "A", 10, 30, 2, rotation_centre(), 6)
#
# Remember, that now the about-pt is the "to" point, i.e. the maps are brought from 
# somewhere else and generated about the about-pt.
#
def transform_map_using_lsq_matrix(imol_ref, ref_chain, ref_resno_start, ref_resno_end, imol_mov, mov_chain, mov_resno_start, mov_resno_end, imol_map, about_pt, radius):

	clear_lsq_matches()
	add_lsq_match(ref_resno_start, ref_resno_end, ref_chain,
		      mov_resno_start, mov_resno_end, mov_chain, 1)
	rtop = apply_lsq_matches(imol_ref, imol_mov)
	transform_map(imol_map, rtop[0], rtop[1], about_pt, radius)

# Make the imol-th map brighter.
#
# Scale_factor > 1 makes brighter...
#
def brighten_map(imol, scale_factor):

    from types import ListType

    if valid_map_molecule_qm(imol):
       current_colour = map_colour_components(imol)
       if type(current_colour) is ListType:
           new_v = []
           for i in range(len(current_colour)):
              new_v.append(current_colour[i] * float(scale_factor))
              if new_v[i] < 0.05:
                  new_v[i] = 0.05
              elif new_v[i] > 1.0:
                  new_v[i] = 1.0
              else:
                  pass
           set_map_colour(imol, *new_v)
       else:
          print "bad non-list current-colour ", current_colour
       graphics_draw()

# Make all maps brighter
#
def brighten_maps():
    map(lambda imap: brighten_map(imap, 1.25), map_molecule_list())

# Make all maps darker
#
def darken_maps():
    map(lambda imap: brighten_map(imap, 0.8), map_molecule_list())

# return a list of chain ids for given molecule number @var{imol}.
# return empty list on error
#
def chain_ids(imol):

    chain_id_is = []
    number_of_chains = n_chains(imol)
    for chain_no in range(number_of_chains):
        chain_id_is.append(chain_id(imol,chain_no))
    return chain_id_is

# convert from interface name to schemisch name to be equivalent to Paul's naming
#
def is_solvent_chain_qm(imol,chain_id):
    if (is_solvent_chain_p(imol,chain_id)==1): return True

# schemey interface to eponymous scripting interface function!?
def valid_model_molecule_qm(imol): 
    if (is_valid_model_molecule(imol)==1): return True
    else: return False

# schemey interface to eponymous scripting interface function.
def valid_map_molecule_qm(imol):
    if (is_valid_map_molecule(imol)==1): return True
    else: return False

# Does residue resno with insertion code ins_code of chain chain_id
# and in molecule number imol exist?
#
# Return True or False
def residue_exists_qm(imol,chain_id,resno,ins_code): 
    if (does_residue_exist_p(imol,chain_id,resno,ins_code)==1): return True
    else: return False

# Return a list of 3 float for the centre of mas of molecule number imol.
#
# on failure return False.
#
def centre_of_mass(imol):
    centre = eval(centre_of_mass_string(imol))
    if (centre == 0):
       print "molecule number",imol,"is not valid"
       return False
    else:
       print "Centre of mass for molecule %s is %s" % (imol, centre)
       # for use somewhere let's return the centre
       return centre

# Return as a list the occupancy, temperature_factor, x y z coordinates
# of the given atom.
# (the x,y,z are in Cartesian Angstroms).
# 
# on error (e.g. atom not found) return False
# 
def atom_specs(imol, chain_id, resno, ins_code, atom_name, alt_conf):

  # not sure what paul is on about port, so just print it!
  # at leats not entirely sure...
  try:
      v = atom_info_string(imol, chain_id, resno, ins_code, atom_name, alt_conf)
  except:
      v = False

  return v

# Paul says:
# backups wrapper: doesn't work currently, I think.  More cleverness required.
def with_no_backups(imol):

    backup_mode = backup_state(imol)
    turn_off_backup(imol)
    if backup_mode == 1:
        turn_on_backup(imol)

# return a guess at the map to be refined (usually called after
# imol_refinement_map returns -1)
#
# Basically uses the first not difference map we find!
#
def guess_refinement_map():

    map_list = map_molecule_list()
    if map_list == []:
        return -1       # failed to find a map
    else:
        for map_mol in map_list:
            if map_is_difference_map(map_mol) == 0:
                return map_mol
            else:
                pass
        print "BL WARNING:: we couldnt find a non difference map for fitting!"
        return -1
        

# Print the sequence of molecule number @var{imol}
#
# This is not really a util, perhaps it should be somewhere else?
#
def print_sequence(imol):
    
    for chain in chain_ids(imol):
       print_sequence_chain(imol,chain)

# comma key hook
def graphics_comma_key_pressed_hook():
	pass

# dot key hook
def graphics_dot_key_pressed_hook():
	pass

# general key press hook
def graphics_general_key_press_hook(key):
	print "Key %s was pressed" %(key)
#	return False
# example

# The P key is already bound
# skip_to_next_ncs_chain is already bound to O key.

# 	if (key==80): 		# bind the P key
# 		set_pointer_atom_is_dummy(1)
# 		place_atom_at_pointer()
# 	if (key==75):		# bind the K key
# 		skip_to_next_ncs_chain()
# 	else:
# 		print "No python binding for key ", key


# def read_vu_file(filename, obj_name):
# BL says:: should do that at some point.

# residue-test-function is a function that takes 4 arguments, the
# chain-id, resno, inscode and residue-serial-number (should it be
# needed) and returns either False or return something interesting
# (e.g. text for a button label).
#
# Return a list of residues, each of which has a return value at the
# start, ie. [return_value, chain_id, res_no, ins_code]
def residue_matching_criteria(imol, residue_test_func):

    chain_list = chain_ids(imol)
    seq_residue_list = []
    alt_conf_residue_list = []

    for chain_id in chain_list:
        for serial_number in range(chain_n_residues(chain_id,imol)):
            res_no = seqnum_from_serial_number(imol, chain_id, serial_number)
            ins_code = insertion_code_from_serial_number(imol, chain_id, serial_number)
            r = residue_test_func(chain_id, res_no, ins_code, serial_number)
            if r:
                seq_residue_list.append([chain_id, res_no, ins_code])
    if seq_residue_list:
        return seq_residue_list
    else:
        return False

# Return a list of all residues that have alt confs: where a residue
# is specified thusly: [[chain_id, resno, ins_code], [...] ]
# 
def residues_with_alt_confs(imol):

    # return False if there are no atoms with alt-confs, else return
    # a list of the residue's spec [chain_id, res_no, ins_code]
    def alt_func1(chain_id, res_no, ins_code, res_serial_no):
        r = False
        atom_ls = residue_info(imol, chain_id, res_no, ins_code)
        for i in range(len(atom_ls)-1):
            alt_conf_str = atom_ls[i][0][1]
            if alt_conf_str:
                r = True
                break
        return r

    return residue_matching_criteria(imol, lambda chain_id, res_no, ins_code, res_serial_no: alt_func1(chain_id, res_no, ins_code, res_serial_no))

# Return a list of all the altconfs in the residue. 
# Typically this will return [] or ["A","B"]
# 
def residue_alt_confs(imol, chain_id, res_no, ins_code):

	atom_ls = residue_info(imol, chain_id, res_no, ins_code)
        alt_confs = []

	if atom_ls:
		for i in range(len(atom_ls)-1):
			alt_conf_str = atom_ls[i][0][1]
			if alt_conf_str:
				if not alt_conf_str in alt_confs:
					alt_confs.append(alt_conf_str)
	return alt_confs

# Return False if no atom can be found given the spec else return a list
# consisting of the atom name and alt-conf specifier.  
# 
# Choose an atom that is called " CA ".  Failing that choose the
# first atom.
# 
def residue_spec2atom_for_centre(imol, chain_id, res_no, ins_code):

    from types import ListType
    # residue-info can return False
    atom_ls = residue_info(imol, chain_id, res_no, ins_code)

    centre_atom_name_alt_conf = False

    if type(atom_ls) is ListType:
        for i in range(len(atom_ls)-1):
            alt_conf_str = atom_ls[i][0][1]
            atom = atom_ls[i][0][0]
            if (atom == " CA "):
                centre_atom_name_alt_conf = [atom, alt_conf_str]
                break
        if (not centre_atom_name_alt_conf) and (len(atom_ls)>0):
            # take the first atom
            centre_atom_name_alt_conf = atom_ls[0][0][0:2]
                 
    return centre_atom_name_alt_conf
 

def update_go_to_atom_from_current_atom():

    active_atom = active_residue()
    if active_atom:
	imol      = active_atom[0]
	chain_id  = active_atom[1]
	resno     = active_atom[2]
	ins_code  = active_atom[3]
	atom_name = active_atom[4]
	alt_conf  = active_atom[5]
        go_to_atom_imol_current = go_to_atom_molecule_number()
        set_go_to_atom_molecule(imol)
        # if imol != goto_atom_imol_current
        update_go_to_atom_window_on_other_molecule_chosen(imol)
        set_go_to_atom_chain_residue_atom_name(chain_id, resno, atom_name)

# Typically one might want to use this on a water, but it deletes the
# nearest CA currently...  Needs a re-think.  Should active-atom just
# return the nearest atom and not be clever about returning a CA.
# 
def delete_atom_by_active_residue():

	active_atom = active_residue()
	if active_atom:
		delete_atom(active_atom)

# general mutate
# 
# typically:
# 
# overlay PTY onto given TYR
#
# delete speced TYR
# 
# merge molecules PTY molecule int molecule number imol
# 
# change resno of PTY to that of the speced TYR
# 
# change chain id of PTY to that of speced TYR
#
# change chain ids with residue range for the PTY
# 
def mutate_by_overlap(imol, chain_id, resno, tlc):

	def mutate_it():
		# BL says:: currently we dont get an imol back from get_monomer (python thing)
		# so we run get_monomer and use the last imol as a ligand imol (may work..)
		imol_ligand = get_monomer(tlc)
		mol_ls = molecule_number_list()
		imol_ligand = mol_ls[len(mol_ls)-1]
		while imol_ligand == -1:
 			print "BL DEBUG:: still waitin", imol_ligand
		if not valid_model_molecule_qm(imol_ligand):
			s = " Oops.  Failed to get monomer " + str(tlc)
			add_status_bar_text(s)
		else:
			delete_residue_hydrogens(imol_ligand, "A", 1, "", "")
			delete_atom(imol_ligand, "A", 1, "", " OXT", "")
			overlap_ligands(imol_ligand, imol, chain_id, resno)
			delete_residue(imol, chain_id, resno, "")
			new_chain_id_info = merge_molecules([imol_ligand], imol)
			print "INFO:: new_chain_id_info: ", new_chain_id_info
			merge_status = new_chain_id_info[0]
			if merge_status == 1:
				new_chain_id = new_chain_id_info[1]
				change_residue_number(imol, new_chain_id, 1, "", resno, "")
				change_chain_id(imol, new_chain_id, chain_id, 1, resno, resno)

				replacement_state = refinement_immediate_replacement_state()
				imol_map = imol_refinement_map()
				set_refinement_immediate_replacement(1)
				if imol_map == -1:
					regularize_zone(imol, chain_id, resno, resno, "")	
					# not sure where to continue
				else:
					spin_atoms = [" P  ", " O1P", " O2P", " O3P"]
					phos_dir = {
						'PTR': [" CZ ", " OH "],
						'SEP': [" CB ", " OG "],
						'TPO': [" CB ", " OG1"] }
					if tlc in phos_dir.keys():
						dir_atoms = phos_dir[tlc]
					else:
						dir_atoms = False
					print ".... spining atoms ", spin_atoms
					if dir_atoms:
						spin_search(imol_map, imol, chain_id, resno, "", dir_atoms, spin_atoms)
					refine_zone(imol, chain_id, resno, resno, "")
				accept_regularizement()
				set_refinement_immediate_replacement(replacement_state)

				set_mol_displayed(imol_ligand, 0)
				set_mol_active(imol_ligand, 0)

			else:
				# guess merge failed!?!
				print "BL WARNING:: merge failed!?"

	# First, if there are multiple maps, force the user to choose one,
	# rather than continuing.
	imol_map = imol_refinement_map()
	import threading
	if (imol_map == -1):
		map_mols = map_molecule_list()
		if len(map_mols) > 1:
# BL says:: I dunno how to wait for the input from map selection since there is no return
# waiting for mouse/keyboard input may be an option but no curses on windows, maybe
# fix later, for now we use the first map
#			show_select_map_dialog()
			print "BL INFO:: we use the first map! If you wanna use another one, please select from Model/Fit/Refine"
			mutate_it()
		else:
			mutate_it()
	elif map_molecule_list() == []:
		print "BL WARNING:: no valid maps around!"
	else:
		mutate_it()

# A bit of fun
#
def phosphorylate_active_residue():

	active_atom = active_residue()
	imol       = active_atom[0]
	chain_id   = active_atom[1]
	resno      = active_atom[2]
	inscode    = active_atom[3]
	res_name = residue_name(imol, chain_id, resno, inscode)

	if res_name == 'TYR':
		mutate_by_overlap(imol, chain_id, resno, "PTR")
	elif res_name == 'SER':
		mutate_by_overlap(imol, chain_id, resno, "SEP")
	elif res_name == 'THR':
		mutate_by_overlap(imol, chain_id, resno, "TPO")
	else:
		s = "Can't Phosphorylate residue of type " + res_name
		info_dialog(s)

def label_all_CAs(imol):
    
    for chain_id in chain_ids(imol):
        if not (is_solvent_chain_qm(imol, chain_id)):
            n_residues = chain_n_residues(chain_id, imol)
            for serial_number in number_list(0, n_residues):
                res_name = resname_from_serial_number(imol, chain_id, serial_number)
                res_no = seqnum_from_serial_number(imol, chain_id, serial_number)
                ins_code = insertion_code_from_serial_number(imol, chain_id, serial_number)
                add_atom_label(imol, chain_id, res_no, " CA ")
    graphics_draw()


def label_all_atoms_in_residue(imol, chain_id, resno, inscode):
	import types

	atom_list = residue_info(imol, chain_id, resno, inscode)
	if type(atom_list) is ListType:
		for atom_info in atom_list:
                        add_atom_label(imol,chain_id, resno, atom_info[0][0])
                graphics_draw()

def label_all_active_residue_atoms():

	active_atom = active_residue()
        imol       = active_atom[0]
        chain_id   = active_atom[1]
        resno      = active_atom[2]
        inscode    = active_atom[3]

	atom_list = residue_info(imol, chain_id, resno, inscode)
	if type(atom_list) is ListType:
		for atom_info in atom_list:
			add_atom_label(imol,chain_id, resno, atom_info[0][0])
		graphics_draw()

#############
# some re-definitions from coot python functions
############

sequence_info    = sequence_info_py
atom_info_string = atom_info_string_py
active_residue   = active_residue_py
closest_atom     = closest_atom_py
residue_info     = residue_info_py
residue_name     = residue_name_py

generic_string_vector_to_list_internal = generic_string_vector_to_list_internal_py
generic_list_to_string_vector_internal = generic_list_to_string_vector_internal_py
generic_int_vector_to_list_internal = generic_int_vector_to_list_internal_py

get_symmetry           = get_symmetry_py
dictionaries_read      = dictionaries_read_py
monomer_restraints     = monomer_restraints_py
set_monomer_restraints = set_monomer_restraints_py
map_colour_components  = map_colour_components_py
cis_peptides           = cis_peptides_py
view_name              = view_name_py
view_description       = view_description_py
go_to_view             = go_to_view_py
merge_molecules        = merge_molecules_py
clear_and_update_molecule = clear_and_update_molecule_py
add_molecule           = add_molecule_py
spin_search            = spin_search_py
merge_molecules        = merge_molecules_py
chain_id_for_shelxl_residue_number = chain_id_for_shelxl_residue_number_py
overlap_ligands        = overlap_ligands_py
analyse_ligand_differences = analyse_ligand_differences_py
execute_ligand_search  = execute_ligand_search_py
movie_file_name_prefix = movie_file_name_prefix_py
apply_lsq_matches      = apply_lsq_matches_py
rotamer_graphs         = rotamer_graphs_py
cootaneer              = cootaneer_py

## and some extra ones
show_set_undo_molecule_chooser = show_set_undo_molecule_chooser_py
centre_of_mass_string          = centre_of_mass_string_py
make_image_raster3d            = make_image_raster3d_py
make_image_povray              = make_image_povray_py
run_state_file                 = run_state_file_py
wrapped_create_run_state_file_dialog =  wrapped_create_run_state_file_dialog_py
raster_screen_shot             = raster_screen_shot_py
handle_read_draw_probe_dots_unformatted = handle_read_draw_probe_dots_unformatted_py
probe_available_p              = probe_available_p_py
gln_asn_b_factor_outliers      = gln_asn_b_factor_outliers_py

############################################################################################
# end of Paul's scripting
############################################################################################
#
# some BL functions
#
############################################################################################
#

# for easier switching on of GL lighting on surfaces:
def GL_light_on():
    set_do_GL_lighting(1)
    do_GL_lighting_state()

# and to turn it off
def GL_light_off():
    set_do_GL_lighting(0)
    do_GL_lighting_state()


# BL module to find exe files
# we need this for popen as it requires the full path of the exe file
# we use arguments:
#
# program_name	: name of exe to find
#
# first_path	: path name to serach first (usually "PATH")
#
# second_path	: next path to serach in (e.g. "CCP4_BIN")
#
# then we search everywhere
# 
# returns full path of exe when successful, False otherwise
#
def find_exe(program_name, first_path, second_path):

    import os, string, fnmatch

    global search_disk
    search_disk = None

    def search_disk_dialog():

        import pygtk, gtk, pango

        dialog = gtk.Dialog("Search whole disk dialog", None,
                 gtk.DIALOG_MODAL | gtk.DIALOG_NO_SEPARATOR,
                 (gtk.STOCK_YES, gtk.RESPONSE_ACCEPT,
                 gtk.STOCK_NO, gtk.RESPONSE_REJECT))
        ifont = gtk.gdk.Font("fixed")
	label_text = "Couldn't find %s in default path" %(program_name)
        if (first_path):
            label_text += " and "
            label_text += first_path
        if (second_path):
            label_text += " and "
            label_text += second_path
        label_text += "\n\nShall we search the whole disk?\n"
        label = gtk.Label(label_text)

        dialog.vbox.pack_end(label, True, True, 0)
	dialog.show_all()
	result = dialog.run()
        if result == gtk.RESPONSE_ACCEPT:
           ret = True
        else:
           ret = False
        dialog.destroy()
        return ret

    if (os.name == 'nt'):
       drive = "C:\\"
       file_ext = file_name_extension(program_name)
       if (file_ext <> 'exe'):
       		extension = ".exe"
       else:
       		extension = ""
    else:
       extension = ""
       drive = "/"

    program_name_noext = program_name
    program_name += extension

    if (first_path):
       try:
          primary_path = os.environ[first_path]
#          print "BL DEBUG:: primary path is ", primary_path
          for path in string.split(primary_path, os.pathsep):
              program_exe = os.path.join(path,program_name)
#              print "BL DEBUG:: prog exe is ", program_exe
              if (os.path.isfile(program_exe)):
                 print "BL INFO:: We found ", program_exe
                 return program_exe
                 break
       except:
          print "BL WARNING:: %s not defined!" %first_path
        
    if (second_path):
       try:
          secondary_path = os.environ[second_path]
#          print "BL DEBUG:: secondary path is ", secondary_path
          for path in string.split(secondary_path, os.pathsep):
              program_exe = os.path.join(path,program_name)
#              print "BL DEBUG:: prog exe is ", program_exe
              if (os.path.isfile(program_exe)):
                 print "BL INFO:: We found ", program_exe
                 return program_exe
                 break
       except:
          print "BL WARNING:: %s not defined!" %first_path

    # BL says: before we search everywhere we might want to ask
    # the user if he actually wishes to do so!
    # lets insert a small pygtk dialog and ask!
    search_disk = search_disk_dialog()
#    print "BL DEBUG:: search_disk is", search_disk
    if search_disk:
       # search everywhere 
       for root, dir, file in os.walk(drive):
           program_exe = os.path.join(root,program_name)
           if (os.path.isfile(program_exe)):
#              print "BL DEBUG:: have found  ", program_exe
              return program_exe
              break
    else:
       print "BL INFO:: we don't search the whole disk for", program_name_noext

    print "BL WARNING:: We cannot find %s anywhere! Program %s won't run!" %(program_name_noext, program_name_noext)
    return False

# for running online docs
def open_url(url):
    import webbrowser

    try:
      webbrowser.open(url,1,1)
    except:
      print "BL WARNING:: Cannot open the URL %s in webbrowser %s!" %(url,webbrowser.get())


# to reload modules
def reload_module(name):
	import os
	path = os.getenv('COOT_PYTHON_DIR')
	file = os.path.join(path, name)
	execfile(file)
